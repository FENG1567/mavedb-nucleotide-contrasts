import argparse,csv,gzip,json,math
from collections import defaultdict,Counter
from datetime import datetime,timezone
from pathlib import Path

REPO_ROOT=Path(__file__).resolve().parents[1]
parser=argparse.ArgumentParser(description='Run the v2 fixed-pair splicing technical audit from local source tables.')
parser.add_argument('--project-root',type=Path,required=True,
    help='source tree containing outputs/opensplice_* and incoming/splice_opensplice')
parser.add_argument('--output-dir',type=Path,default=REPO_ROOT/'results/splicing_audit_v2')
args=parser.parse_args()
ROOT=args.project_root.expanduser().resolve()
PAIR=args.output_dir.expanduser().resolve()
PAIR.mkdir(parents=True,exist_ok=True)
AA=ROOT/'outputs'/'opensplice_compass_aa.csv'
OV=ROOT/'outputs'/'opensplice_compass_overlap.csv'
OP=ROOT/'outputs'/'opensplice_overlap_psi.csv'
CP=ROOT/'incoming'/'splice_opensplice'/'compass'/'compass_psi.csv.gz'
MASTER=ROOT/'incoming'/'splice_opensplice'/'data'/'opensplice_master.tsv.gz'
MAN=PAIR/'splicing_fixed_pair_manifest_20261123_v2.csv'
OUTG=PAIR/'splicing_fixed_pairs_20261123_v2.csv'
OUTV=PAIR/'splicing_fixed_pair_variants_20261123_v2.csv'
OUTJ=PAIR/'splicing_fixed_pair_summary_20261123_v2.json'
OUTM=PAIR/'splicing_fixed_pair_methods_20261123_v2.md'
CELLS=['HeLa','HEK','K562','MCF7','HMC3']

def num(x):
    try:
        y=float(str(x).strip())
        return y if math.isfinite(y) else None
    except Exception:
        return None

def nv(x):
    if not x: return ''
    p=[z.strip() for z in str(x).strip().split(':')]
    if len(p)==3 and '>' in p[2]:
        ref,alt=p[2].split('>',1)
    elif len(p)==4:
        ref,alt=p[2],p[3]
    else:
        return ''
    c,pos=p[0],p[1]
    if not c.lower().startswith('chr'): c='chr'+c
    try: pos=str(int(float(pos)))
    except Exception: return ''
    return f'{c}:{pos}:{ref.upper()}:{alt.upper()}'

def rk(r):
    return nv(':'.join([r.get('CHROM',''),r.get('POS',''),r.get('REF',''),r.get('ALT','')]))

def vk(v):
    c,p,ref,alt=v.split(':')
    c=c.lower().removeprefix('chr')
    rank=int(c) if c.isdigit() else {'x':23,'y':24,'m':25,'mt':25}.get(c,26)
    return rank,int(p),ref,alt

def tn(x):
    return ' '.join(str(x or '').strip().lower().split())

def ranks(x):
    ix=sorted(range(len(x)),key=lambda i:x[i]); out=[0.0]*len(x); i=0
    while i<len(ix):
        j=i+1
        while j<len(ix) and x[ix[j]]==x[ix[i]]: j+=1
        q=(i+1+j)/2.0
        for z in ix[i:j]: out[z]=q
        i=j
    return out

def corr(x,y):
    if len(x)<2: return None
    mx=sum(x)/len(x); my=sum(y)/len(y)
    a=[z-mx for z in x]; b=[z-my for z in y]
    d=math.sqrt(sum(z*z for z in a)*sum(z*z for z in b))
    return sum(i*j for i,j in zip(a,b))/d if d else None

def fisher(r,n):
    if r is None or n<=3 or abs(r)>=1: return None
    z=math.atanh(r); se=1/math.sqrt(n-3)
    return [math.tanh(z-1.95996398454*se),math.tanh(z+1.95996398454*se)]

def wilson(k,n):
    if not n: return None
    z=1.95996398454; p=k/n; d=1+z*z/n
    c=(p+z*z/(2*n))/d; h=z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/d
    return [max(0,c-h),min(1,c+h)]

def bp(k,n):
    if not n: return {'one_sided':None,'two_sided':None}
    d=2**n
    up=sum(math.comb(n,i) for i in range(k,n+1))/d
    lo=sum(math.comb(n,i) for i in range(k+1))/d
    return {'one_sided':up,'two_sided':min(1.0,2*min(up,lo))}

print('Fixing all candidate groups and pair members before reading score fields.',flush=True)
gs=defaultdict(lambda:{'vars':set(),'tx':set()})
with open(AA,newline='',encoding='utf-8-sig') as f:
    for r in csv.DictReader(f):
        v=rk(r)
        try: pos=int(float(r.get('protein_pos','')))
        except Exception: continue
        if not v or not r.get('gene') or not r.get('ref_aa') or not r.get('alt_aa'): continue
        k=(r['gene'].strip(),pos,r['ref_aa'].strip(),r['alt_aa'].strip())
        gs[k]['vars'].add(v)
        if r.get('transcript'): gs[k]['tx'].add(r['transcript'].strip())
candidates=[]
for k,d in sorted(gs.items()):
    vs=sorted(d['vars'],key=vk)
    if len(vs)>=2:
        candidates.append({'gene':k[0],'pos':k[1],'ref':k[2],'alt':k[3],
          'n':len(vs),'tx':sorted(d['tx']),'all':vs,'pair':vs[:2]})
for i,g in enumerate(candidates,1): g['id']=f'SG{i:03d}'
fixed={v for g in candidates for v in g['pair']}
allvars={v for g in candidates for v in g['all']}
assert len(candidates)==22, f'expected 22 candidate AA groups, observed {len(candidates)}'
assert len(allvars)==46, f'expected 46 pre-score candidate variants, observed {len(allvars)}'
assert len(fixed)==44, f'expected 44 fixed pair members, observed {len(fixed)}'

manifest=[]
for g in candidates:
    for rank,v in enumerate(g['all'],1):
        manifest.append({'group_id':g['id'],'gene':g['gene'],'protein_pos':g['pos'],'ref_aa':g['ref'],'alt_aa':g['alt'],
          'candidate_variant_n':g['n'],'candidate_rank_by_coordinate':rank,'variant_hg38':v,
          'selected_for_fixed_pair':v in g['pair'],'selected_pair_rank':(g['pair'].index(v)+1 if v in g['pair'] else ''),
          'transcripts_seen':'|'.join(g['tx']),
          'pair_rule':'first two unique DNA keys after chromosome/position/ref/alt sorting; fixed before score lookup'})
with open(MAN,'w',newline='',encoding='utf-8') as f:
    w=csv.DictWriter(f,fieldnames=list(manifest[0])); w.writeheader(); w.writerows(manifest)
print(f'Frozen manifest: {len(candidates)} groups, {len(allvars)} original candidates, {len(fixed)} fixed pair members.',flush=True)

source_sample=''
with gzip.open(CP,'rt',newline='',encoding='utf-8') as f:
    for r in csv.DictReader(f):
        if r.get('variant_hg38'):
            source_sample=str(r['variant_hg38']).strip()
            break
assert source_sample and '>' in source_sample and nv(source_sample), f'unexpected actual COMPASS key format: {source_sample!r}'
parts=source_sample.split(':')
assert len(parts)==3 and '>' in parts[2]
sr,sa=parts[2].split('>',1)
assert nv(source_sample)==nv(':'.join([parts[0],parts[1],sr,sa])), 'source-format parser round-trip failed'
fixed_example=sorted(fixed,key=vk)[0]
fc,fp,fr,fa=fixed_example.split(':')
fixed_source_form=f'{fc}:{fp}:{fr}>{fa}'
assert nv(fixed_source_form)==fixed_example, 'fixed DNA key source-format round-trip failed'
print(json.dumps({'actual_COMPASS_key_sample':source_sample,'normalized_sample':nv(source_sample),
                  'fixed_key_example':fixed_example,'fixed_source_format_example':fixed_source_form},sort_keys=True),flush=True)

exons=defaultdict(set)
overlap_keys=set()
with open(OV,newline='',encoding='utf-8-sig') as f:
    for r in csv.DictReader(f):
        v=rk(r)
        if v: overlap_keys.add(v)
        if v in allvars and r.get('compass_exon'): exons[v].add(r['compass_exon'].strip())
assert fixed.issubset(overlap_keys), f'fixed keys absent from exact overlap map: {sorted(fixed-overlap_keys)}'
oprows=defaultdict(list)
with open(OP,newline='',encoding='utf-8-sig') as f:
    for r in csv.DictReader(f):
        v=rk(r)
        if v in allvars: oprows[v].append(r)

cprows=defaultdict(list)
print('Scanning COMPASS once for all 46 pre-score candidates and all cell-line columns.',flush=True)
with gzip.open(CP,'rt',newline='',encoding='utf-8') as f:
    for r in csv.DictReader(f):
        raw=str(r.get('variant_hg38','')).strip()
        v=nv(raw)
        if v in allvars:
            cprows[v].append(r)
print('COMPASS scan finished.',flush=True)

variant_rows=[]; vlookup={}
for v in sorted(allvars,key=vk):
    ors=oprows.get(v,[])
    if len(ors)==1:
        op=ors[0]; oscore=num(op.get('delta_logit'))
        ost='unique_row_with_score' if oscore is not None else 'unique_row_missing_delta_logit'
    elif not ors: op={}; oscore=None; ost='no_exact_variant_row'
    else: op={}; oscore=None; ost='multiple_derived_rows_for_variant'
    es=sorted(exons.get(v,set())); raw=cprows.get(v,[])
    if len(es)!=1:
        matches=[]; cst='no_compass_exon_mapping' if not es else 'multiple_compass_exon_mappings'
    else:
        matches=[r for r in raw if tn(r.get('gene_exon'))==tn(es[0])]
        if not raw: cst='exact_variant_absent_in_COMPASS'
        elif not matches: cst='gene_exon_mismatch'
        elif len(matches)>1: cst='multiple_rows_for_variant_and_gene_exon'
        else: cst='unique_event_row'
    cp=matches[0] if len(matches)==1 else {}
    vr={'variant_hg38':v,'candidate_for_fixed_pair':v in fixed,
      'OpenSplice_rows':len(ors),'OpenSplice_gene':op.get('gene',''),'OpenSplice_exon_id':op.get('exon_id',''),
      'OpenSplice_region':op.get('region',''),'OpenSplice_psi':num(op.get('psi')) if len(ors)==1 else None,
      'OpenSplice_delta_psi':num(op.get('delta_psi')) if len(ors)==1 else None,'OpenSplice_delta_logit':oscore,
      'OpenSplice_status':ost,'COMPASS_exon_labels':'|'.join(es),'COMPASS_raw_rows_for_variant':len(raw),
      'COMPASS_raw_event_ids':'|'.join(sorted({str(r.get('event_id','')) for r in raw if r.get('event_id')})),
      'COMPASS_raw_gene_exons':'|'.join(sorted({str(r.get('gene_exon','')) for r in raw if r.get('gene_exon')})),
      'COMPASS_matching_event_rows':len(matches),'COMPASS_matching_event_id':cp.get('event_id',''),
      'COMPASS_matching_gene_exon':cp.get('gene_exon',''),'COMPASS_matching_transcript_id':cp.get('transcript_id',''),
      'COMPASS_status':cst}
    for cell in CELLS:
        vr[f'COMPASS_{cell}_delta_logit_pooled']=num(cp.get(f'{cell}_delta_logit_pooled')) if cp else None
        vr[f'COMPASS_{cell}_dpsi_pooled']=num(cp.get(f'{cell}_dpsi_pooled')) if cp else None
        vr[f'COMPASS_{cell}_pooled_psi_raw']=num(cp.get(f'{cell}_pooled_psi_raw')) if cp else None
        vr[f'COMPASS_{cell}_total_pooled']=num(cp.get(f'{cell}_total_pooled')) if cp else None
        vr[f'COMPASS_{cell}_status']='unique_event_row_with_score' if cp and num(cp.get(f'{cell}_delta_logit_pooled')) is not None else ('unique_event_row_missing_score' if cp else cst)
    variant_rows.append(vr); vlookup[v]=vr

group_rows=[]
for g in candidates:
    v1,v2=g['pair']; a,b=vlookup[v1],vlookup[v2]
    row={'group_id':g['id'],'gene':g['gene'],'protein_pos':g['pos'],'ref_aa':g['ref'],'alt_aa':g['alt'],
      'candidate_variant_n':g['n'],'transcripts_seen':'|'.join(g['tx']),
      'pair_rule':'first two unique DNA keys by chromosome/position/ref/alt, fixed before score lookup',
      'variant_1':v1,'variant_2':v2,'variant_1_OpenSplice_status':a['OpenSplice_status'],
      'variant_2_OpenSplice_status':b['OpenSplice_status'],'variant_1_COMPASS_event_status':a['COMPASS_status'],
      'variant_2_COMPASS_event_status':b['COMPASS_status']}
    dop=b['OpenSplice_delta_logit']-a['OpenSplice_delta_logit'] if a['OpenSplice_delta_logit'] is not None and b['OpenSplice_delta_logit'] is not None else None
    row['OpenSplice_delta_logit_diff_v2_minus_v1']=dop
    for cell in ['HeLa','HEK']:
        x=a[f'COMPASS_{cell}_delta_logit_pooled']; y=b[f'COMPASS_{cell}_delta_logit_pooled']
        dcp=y-x if x is not None and y is not None else None
        row[f'COMPASS_{cell}_delta_logit_diff_v2_minus_v1']=dcp
        complete=dop is not None and dcp is not None
        if not complete: direction='unavailable'
        elif abs(dop)<=1e-12 or abs(dcp)<=1e-12: direction='tie'
        else: direction='concordant' if dop*dcp>0 else 'discordant'
        row[f'{cell}_score_status']='complete' if complete else ('missing_OpenSplice' if dop is None else f'missing_COMPASS_{cell}')
        row[f'{cell}_direction_status']=direction
    group_rows.append(row)

def metrics_for(cell):
    scored=[r for r in group_rows if r[f'{cell}_score_status']=='complete']
    non=[r for r in scored if r[f'{cell}_direction_status'] in ('concordant','discordant')]
    k=sum(r[f'{cell}_direction_status']=='concordant' for r in non)
    d=sum(r[f'{cell}_direction_status']=='discordant' for r in non)
    ties=sum(r[f'{cell}_direction_status']=='tie' for r in scored)
    missing=len(group_rows)-len(scored)
    x=[r['OpenSplice_delta_logit_diff_v2_minus_v1'] for r in scored]
    y=[r[f'COMPASS_{cell}_delta_logit_diff_v2_minus_v1'] for r in scored]
    pr=corr(x,y); sr=corr(ranks(x),ranks(y)) if len(x)>1 else None
    return {'complete_fixed_pairs':len(scored),'non_tie_pairs':len(non),'concordant':k,'discordant':d,'ties':ties,'unavailable':missing,
      'R_among_scored_non_tie_pairs':k/len(non) if non else None,'R_Wilson95CI_conditional':wilson(k,len(non)),
      'exact_binomial_vs_0.5':bp(k,len(non)),'Pearson_r':pr,'Pearson_Fisher95CI':fisher(pr,len(scored)),
      'Spearman_rho':sr,'n_for_correlations':len(scored),
      'all_22_candidate_R_lower_upper_bounds':[k/len(group_rows),(k+ties+missing)/len(group_rows)] if group_rows else [None,None],
      'status':'provisional descriptive only; OpenSplice event provenance is incomplete'}

with open(OUTG,'w',newline='',encoding='utf-8') as f:
    w=csv.DictWriter(f,fieldnames=list(group_rows[0])); w.writeheader(); w.writerows(group_rows)
with open(OUTV,'w',newline='',encoding='utf-8') as f:
    w=csv.DictWriter(f,fieldnames=list(variant_rows[0])); w.writeheader(); w.writerows(variant_rows)

size=MASTER.stat().st_size if MASTER.exists() else None
magic=b''
if MASTER.exists():
    with open(MASTER,'rb') as f: magic=f.read(8)
summary={'created_utc':datetime.now(timezone.utc).isoformat(),
 'pre_score_manifest':str(MAN),'candidate_groups_expected':22,'candidate_groups_found':len(candidates),
 'original_candidates_expected':46,'original_candidates_found':len(allvars),'fixed_pair_members':len(fixed),
 'candidate_group_size_counts':dict(Counter(str(g['n']) for g in candidates)),
 'selection_rule':'Within each gene + protein position + ref AA + alt AA group, deduplicate DNA keys; select the first two by chromosome/position/ref/alt order before score lookup.',
 'source_key_preflight':{'actual_COMPASS_source_value':source_sample,'normalized_value':nv(source_sample),
   'source_value_format':'chr:position:REF>ALT','fixed_key_example':fixed_example,'fixed_key_source_format_example':fixed_source_form,
   'fixed_key_roundtrip_passed':nv(fixed_source_form)==fixed_example,'fixed_keys_present_in_exact_overlap_map':len(fixed.intersection(overlap_keys))},
 'source_provenance':{'candidate_AA_map':str(AA),'exact_overlap_event_map':str(OV),'OpenSplice_derived_score_extract':str(OP),
   'OpenSplice_master_path':str(MASTER),'OpenSplice_master_bytes':size,'OpenSplice_master_gzip_magic':magic[:2]==bytes([31,139]),
   'OpenSplice_master_HTML_prefix':magic.lstrip().lower().startswith(b'<html'),
   'OpenSplice_access':'Saved master is 146-byte HTTP 403 HTML, not gzip; current official Figshare article API call returned HTTP 403 text/html. Event-level OpenSplice identity remains unverified.',
   'COMPASS_path':str(CP),'COMPASS_bytes':CP.stat().st_size,'COMPASS_join_fields':['variant_hg38','gene_exon','event_id','transcript_id'],
   'COMPASS_cells_retained':CELLS},
 'score_definitions':{'OpenSplice':'derived delta_logit matched by exact DNA key; source extract lacks transcript, strand, event, and condition.',
   'COMPASS':'GSE307247 pooled delta_logit for all five cell lines, matched by exact DNA key plus unique gene_exon label.',
   'primary_reproduction':'HeLa endpoint, matching the original internal concordance script.',
   'secondary_condition_sensitivity':'HEK endpoint retained from the official OpenSplice external-comparison R script; condition matching is not claimed because OpenSplice context is absent.'},
 'coverage':{'all_original_candidate_variants':len(allvars),'fixed_pair_members':len(fixed),
   'OpenSplice_unique_score_variants':sum(v['OpenSplice_delta_logit'] is not None for v in variant_rows),
   'COMPASS_unique_event_rows':sum(v['COMPASS_matching_event_rows']==1 for v in variant_rows),
   'COMPASS_HeLa_score_variants':sum(v['COMPASS_HeLa_delta_logit_pooled'] is not None for v in variant_rows),
   'COMPASS_HEK_score_variants':sum(v['COMPASS_HEK_delta_logit_pooled'] is not None for v in variant_rows)},
 'primary_HeLa_metrics':metrics_for('HeLa'),'secondary_HEK_metrics':metrics_for('HEK'),
 'interpretation':'The fixed-pair table corrects the earlier 0-match parser error and retains all 22 groups, all 46 pre-score candidate variants, and all 44 fixed pair members. Any R, p-value, or correlation is provisional descriptive bookkeeping only; do not present the former complete-case result as inferential cross-study validation until OpenSplice event/transcript/strand/condition identity is recovered.',
 'run1_parser_error_note':'Run 1 outputs are preserved separately. Its COMPASS join parser incorrectly required four colon-delimited fields; the observed source format is chr:pos:REF>ALT. Run 1 zero matches are a parser failure, not assay/clinical missingness.'}
OUTJ.write_text(json.dumps(summary,indent=2,ensure_ascii=False,sort_keys=True)+chr(10),encoding='utf-8')
notes=[
 '# OpenSplice × COMPASS fixed-pair selection audit (corrected parser)',
 '',
 'The score-independent manifest records all 46 candidate DNA variants in 22 groups and marks the 44 fixed pair members. The two variants per group were selected by genomic coordinate order before any score field was read.',
 '',
 'The first pass used an incorrect COMPASS key parser and produced zero joins. Those run 1 outputs remain preserved as a parser-failure record and are not evidence of missing assay or clinical data. The corrected parser handles an observed raw source value in chr:position:REF>ALT format, round-trips a fixed DNA key into that format, and asserts all 44 fixed pair members occur in the independent 2,073-key exact-overlap table.',
 '',
 'COMPASS event_id, gene_exon, transcript_id, and all five cell-line score columns are retained for the 46 original candidates. HeLa reproduces the prior internal endpoint; HEK is retained as a secondary condition sensitivity from the official OpenSplice comparison script. OpenSplice scores still come from the derived extract because the master file is a 146-byte 403 HTML response, not a valid gzip. The extract lacks transcript, strand, event, and condition fields.',
 '',
 'R, binomial p, and correlations are provisional descriptive values only. Do not claim inferential cross-study validation unless OpenSplice source identity is recovered.'
]
OUTM.write_text(chr(10).join(notes)+chr(10),encoding='utf-8')
print(json.dumps({'summary':str(OUTJ),'groups':len(group_rows),'candidates':len(allvars),'fixed_members':len(fixed),
 'coverage':summary['coverage'],'HeLa':summary['primary_HeLa_metrics'],'HEK':summary['secondary_HEK_metrics'],
 'source_key_preflight':summary['source_key_preflight']},sort_keys=True),flush=True)


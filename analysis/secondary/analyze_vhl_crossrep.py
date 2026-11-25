#!/usr/bin/env python3
"""Prespecified within-study VHL RNA/function cross-replicate analysis."""
import argparse
from pathlib import Path
from hashlib import sha256
import io, json, re, zipfile
import numpy as np
import pandas as pd
import openpyxl

REPO_ROOT = Path(__file__).resolve().parents[2]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--input-zip', type=Path, required=True,
                    help='article supplementary ZIP containing the ST1 workbook')
parser.add_argument('--output-dir', type=Path, default=REPO_ROOT / 'results' / 'secondary' / 'vhl',
                    help='directory for generated analysis outputs')
args = parser.parse_args()
ZIP = args.input_zip.resolve()
OUT = args.output_dir.resolve()
OUT.mkdir(parents=True, exist_ok=True)
MEMBER = '41588_2024_1800_MOESM4_ESM.xlsx'
BOOT = 2000
SEED = 20261122

with zipfile.ZipFile(ZIP) as zf:
    payload=zf.read(MEMBER)
    wb=openpyxl.load_workbook(io.BytesIO(payload),read_only=True,data_only=True)
    ws=wb['ST1']
    header=[x for x in next(ws.iter_rows(min_row=3,max_row=3,values_only=True))]
    rows=list(ws.iter_rows(min_row=4,values_only=True))
    wb.close()
df=pd.DataFrame(rows,columns=header).dropna(how='all')
input_sha=sha256(payload).hexdigest()

required=['hg38_pos','ref','alt','cHGVS','pHGVS','conseq','protPos',
          'tHDR_pre','tHDR_post','tHDR_rna','rLD2_tHDR_pre','rLD2_tHDR_post','rLD2_tHDR_rna']
missing_cols=[c for c in required if c not in df.columns]
if missing_cols: raise RuntimeError(f'Missing required columns: {missing_cols}')
for c in ['hg38_pos','protPos','tHDR_pre','tHDR_post','tHDR_rna','rLD2_tHDR_pre','rLD2_tHDR_post','rLD2_tHDR_rna']:
    df[c]=pd.to_numeric(df[c],errors='coerce')
for c in ['ref','alt','cHGVS','pHGVS','conseq']:
    df[c]=df[c].astype('string').str.strip()

before=len(df)
missense=df[df['conseq'].str.upper().eq('NON_SYNONYMOUS') & df['cHGVS'].notna() & df['pHGVS'].notna() & df['hg38_pos'].notna() & df['ref'].notna() & df['alt'].notna()].copy()
missense['dna_key']=missense.apply(lambda r:f"hg38:{int(r.hg38_pos)}:{r.ref}>{r.alt}",axis=1)

# Collapse only repeated exact-DNA records with consistent annotations and identical counts.
count_cols=['tHDR_pre','tHDR_post','tHDR_rna','rLD2_tHDR_pre','rLD2_tHDR_post','rLD2_tHDR_rna']
dedup_rows=[]; inconsistent=[]; duplicated_keys=0; duplicate_conflicts=0
for key,g in missense.groupby('dna_key',sort=True,dropna=False):
    if g['cHGVS'].nunique(dropna=False)!=1 or g['pHGVS'].nunique(dropna=False)!=1:
        inconsistent.append(key); continue
    if len(g)>1:
        duplicated_keys+=1
        if not all(g[c].nunique(dropna=False)==1 for c in count_cols):
            duplicate_conflicts+=1; continue
    dedup_rows.append(g.sort_values(['sge_region'] if 'sge_region' in g.columns else ['dna_key'],kind='stable').iloc[0])
u=pd.DataFrame(dedup_rows).reset_index(drop=True)
if u.empty: raise RuntimeError('No unique consistent missense DNA rows remain')

# Eligibility uses day-6 DNA depth only, in both source transfection replicates.
u['DNA6_1']=u['tHDR_pre']; u['DNA20_1']=u['tHDR_post']; u['RNA6_1']=u['tHDR_rna']
u['DNA6_2']=u['rLD2_tHDR_pre']; u['DNA20_2']=u['rLD2_tHDR_post']; u['RNA6_2']=u['rLD2_tHDR_rna']
depth_ok=u['DNA6_1'].ge(10) & u['DNA6_2'].ge(10)
u['depth_eligible']=depth_ok
u['pHGVS_group']=u['pHGVS']
group_counts=u[depth_ok].groupby('pHGVS_group')['dna_key'].nunique()
eligible_groups=set(group_counts[group_counts.ge(2)].index)
u['in_eligible_group']=depth_ok & u['pHGVS_group'].isin(eligible_groups)
elig=u[u['in_eligible_group']].copy()
if elig.empty: raise RuntimeError('No eligible same-pHGVS groups after prespecified DNA depth filter')

def logdiff(num,den):
    return np.log2(np.asarray(num,dtype=float)+1)-np.log2(np.asarray(den,dtype=float)+1)
elig['F1']=logdiff(elig['DNA20_1'],elig['DNA6_1'])
elig['R1']=logdiff(elig['RNA6_1'],elig['DNA6_1'])
elig['F2']=logdiff(elig['DNA20_2'],elig['DNA6_2'])
elig['R2']=logdiff(elig['RNA6_2'],elig['DNA6_2'])
elig['RNA6_1_zero']=elig['RNA6_1'].eq(0); elig['RNA6_2_zero']=elig['RNA6_2'].eq(0)

PAIRINGS=[('within_rep1','F1','R1'),('within_rep2','F2','R2'),('cross_F1_R2','F1','R2'),('cross_F2_R1','F2','R1')]
groups={k:g.sort_values('cHGVS',kind='stable').reset_index(drop=True) for k,g in elig.groupby('pHGVS_group',sort=True)}
group_pos={k:int(g['protPos'].iloc[0]) for k,g in groups.items()}
pos_groups={}
for k,pos in group_pos.items(): pos_groups.setdefault(pos,[]).append(k)
positions=np.array(sorted(pos_groups.keys()),dtype=int)

def describe(group_keys,xcol,ycol):
    covs=[]; cx=[]; cy=[]; concordant=0; sign_n=0; ties=0
    for key in group_keys:
        g=groups[key]
        x=g[xcol].to_numpy(float); y=g[ycol].to_numpy(float)
        if len(x)<2: continue
        covs.append(float(np.cov(x,y,ddof=1)[0,1]))
        cx.extend((x-x.mean()).tolist()); cy.extend((y-y.mean()).tolist())
        dx=x[0]-x[1]; dy=y[0]-y[1]
        if dx==0 or dy==0: ties+=1
        else:
            sign_n+=1
            concordant+=int(np.sign(dx)==np.sign(dy))
    if not covs: return {'mean_group_covariance':None,'groupcentered_r':None,'first2_sign_concordance':None,'sign_groups_n':0,'sign_ties':ties}
    if np.std(cx)==0 or np.std(cy)==0: corr=None
    else: corr=float(np.corrcoef(cx,cy)[0,1])
    return {'mean_group_covariance':float(np.mean(covs)),'groupcentered_r':corr,
            'first2_sign_concordance':(concordant/sign_n if sign_n else None),
            'sign_groups_n':sign_n,'sign_ties':ties}

point=[]
for name,x,y in PAIRINGS:
    point.append({'pairing':name,'function_measure':x,'rna_measure':y,**describe(list(groups),x,y)})

# Position-cluster bootstrap; every sampled position carries all its pHGVS groups.
rng=np.random.default_rng(SEED)
boot={p:{'mean_group_covariance':[],'groupcentered_r':[],'first2_sign_concordance':[]} for p,_,_ in PAIRINGS}
for _ in range(BOOT):
    sampled=rng.choice(positions,size=len(positions),replace=True)
    keys=[]
    for pos in sampled: keys.extend(pos_groups[int(pos)])
    for name,x,y in PAIRINGS:
        stat=describe(keys,x,y)
        for metric in boot[name]:
            val=stat.get(metric)
            if val is not None and np.isfinite(val): boot[name][metric].append(val)

summary=[]
for item in point:
    name=item['pairing']; b=boot[name]
    for metric in ('mean_group_covariance','groupcentered_r','first2_sign_concordance'):
        vals=np.asarray(b[metric],dtype=float)
        item[metric+'_ci95_low']=float(np.quantile(vals,.025)) if len(vals) else None
        item[metric+'_ci95_high']=float(np.quantile(vals,.975)) if len(vals) else None
        item[metric+'_bootstrap_n']=int(len(vals))
    summary.append(item)

group_audit=[]
for key,g in u.groupby('pHGVS_group',sort=True):
    group_audit.append({'pHGVS':key,'protPos':g['protPos'].iloc[0], 'unique_DNA_pre_filter':int(g['dna_key'].nunique()),
                        'eligible_DNA_after_both_pre10':int(g.loc[g['depth_eligible'],'dna_key'].nunique()),
                        'retained_for_analysis':key in eligible_groups})

variant_cols=['dna_key','hg38_pos','ref','alt','cHGVS','pHGVS','protPos','function_class','DNA6_1','DNA20_1','RNA6_1','DNA6_2','DNA20_2','RNA6_2','F1','R1','F2','R2','RNA6_1_zero','RNA6_2_zero']
elig[[c for c in variant_cols if c in elig.columns]].sort_values(['pHGVS','cHGVS']).to_csv(OUT/'variant_level.csv',index=False)
pd.DataFrame(group_audit).to_csv(OUT/'group_audit.csv',index=False)
pd.DataFrame(summary).to_csv(OUT/'crossrep_summary.csv',index=False)

counts={'source_workbook_member':MEMBER,'source_member_sha256':input_sha,'source_bundle_sha256':sha256(ZIP.read_bytes()).hexdigest(),
        'source_rows':int(before),'missense_rows_with_exact_ids':int(len(missense)),'unique_consistent_exact_DNA_rows':int(len(u)),
        'duplicate_DNA_keys_collapsed_identical_counts':int(duplicated_keys-duplicate_conflicts),
        'duplicate_DNA_keys_with_conflicting_counts_excluded':int(duplicate_conflicts),
        'inconsistent_cHGVS_or_pHGVS_keys_excluded':int(len(inconsistent)),
        'unique_DNA_rows_below_day6_DNA10_either_replicate':int((~depth_ok).sum()),
        'unique_DNA_rows_meeting_both_day6_DNA10':int(depth_ok.sum()),
        'same_pHGVS_groups_before_depth':int(u['pHGVS_group'].nunique()),
        'same_pHGVS_groups_retained':int(len(eligible_groups)),
        'retained_DNA_variants':int(len(elig)),'retained_protein_positions':int(elig['protPos'].nunique()),
        'RNA6_rep1_zero_count':int(elig['RNA6_1_zero'].sum()),'RNA6_rep2_zero_count':int(elig['RNA6_2_zero'].sum()),
        'bootstrap':{'cluster':'protPos','iterations':BOOT,'seed':SEED,'positions_sampled_each_iteration':int(len(positions))},
        'summary':summary}
(OUT/'summary.json').write_text(json.dumps(counts,indent=2,ensure_ascii=False),encoding='utf-8')
(OUT/'methods.md').write_text(
    '# VHL cross-replicate analysis result\n\n'
    'The prespecified within-study analysis uses exact missense cHGVS/genomic variants grouped by identical pHGVS. Variant inclusion is determined only by day-6 gDNA counts >=10 in both transfection replicates; each retained pHGVS group has >=2 distinct DNA variants. RNA and day-20 zeros remain, using log2(count+1).\n\n'
    'F1/F2 are day-20/day-6 gDNA log differences; R1/R2 are day-6 RNA/day-6 gDNA log differences for each replicate. Four pairings (within replicate 1/2 and both cross-replicate directions) are reported without selecting a preferred direction. Outcomes are group-balanced mean covariance, variant-level group-centered Pearson r, and sign concordance for the first two cHGVS-sorted DNA variants in each pHGVS group. 95% percentile intervals use 2,000 bootstrap resamples of protein-position clusters.\n\n'
    'No pooled `function_score_final`/`rna_score` correlation was used. This is an internal replicate measurement check in one study; it cannot establish causality, an independent-study replication, or a new mechanism. See CSVs and `summary.json` for attrition and estimates.\n',encoding='utf-8')
print(json.dumps(counts,ensure_ascii=False))

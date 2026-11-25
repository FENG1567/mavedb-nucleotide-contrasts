#!/usr/bin/env python3
"""Exact group reconciliation from frozen hash-collapsed outputs and bounded 27-set source check."""
import argparse,csv,json,math,hashlib
from collections import defaultdict,Counter
from pathlib import Path
from datetime import datetime,timezone

REPO_ROOT=Path(__file__).resolve().parents[1]
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--project-root',type=Path,required=True,
    help='original analysis tree with revision_20261122/dna_audit, work/multigene_mapped and incoming/multigene_scores')
parser.add_argument('--output-dir',type=Path,default=REPO_ROOT/'results/collection_raw_audit')
args=parser.parse_args()
ROOT=args.project_root.expanduser().resolve()
AUD=ROOT/'revision_20261122/dna_audit'; MAP=ROOT/'work/multigene_mapped'
RAW=ROOT/'incoming/multigene_scores'; OUT=args.output_dir.expanduser().resolve()
OUT.mkdir(parents=True,exist_ok=True)
def rows(path):
    with open(path,newline='',encoding='utf-8-sig') as f: yield from csv.DictReader(f)
def read(path): return list(rows(path))
def write(path,data):
    with open(path,'w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(data[0]) if data else [],extrasaction='ignore');w.writeheader();w.writerows(data)
def yes(v):return str(v).strip().lower() in {'true','1','yes'}
def number(v):
    try:
        x=float(v);return x if math.isfinite(x) else None
    except:return None
def score_col(fields):
    for c in ('score','Score','SCORE'):
        if c in fields:return c
    return next((c for c in fields if 'score' in c.lower() and 'error' not in c.lower()),None)
def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()
def ids(s):return sorted(set(x for x in str(s or '').split(';') if x))

print('STAGE 1 load crosswalk and precomputed qualification groups',flush=True)
legacy=read(AUD/'legacy_key_crosswalk.csv')
old_by_group=defaultdict(list);key_to_groups=defaultdict(list);set_gene={}
for r in legacy:
    g=(r['score_set'],r['gene'],str(r['protein_pos']),r['ref_aa'],r['alt_aa'])
    old_by_group[g].append(r);key_to_groups[(r['score_set'],r['legacy_key'])].append(g);set_gene[r['score_set']]=r['gene']
ref_groups=defaultdict(dict);strict_groups=defaultdict(dict)
def load_qual(path,nfield,out):
    for r in rows(path):
        g=(r['score_set'],r['gene'],str(r['protein_pos']),r['ref_aa'],r['alt_aa'])
        acc=r['transcript_accession'];vals=set(ids(r.get('dna_identities','')))
        n=int(r[nfield])
        if len(vals)!=n or n<2:raise ValueError(f'identity-list mismatch or <2 threshold: {path} {g} {acc}')
        if acc in out[g]:raise ValueError(f'duplicate transcript-qualified group row: {path} {g} {acc}')
        out[g][acc]=vals
load_qual(AUD/'reference_supported_groups.csv','unique_reference_supported_DNA_count',ref_groups)
load_qual(AUD/'strict_groups.csv','unique_verified_DNA_count',strict_groups)
per={r['score_set']:r for r in read(AUD/'per_score_set.csv')}
old_per=read(AUD/'legacy_per_set.csv')
eligible_sets={r['score_set'] for r in old_per if int(r.get('legacy_groups_ge2_unique_keys') or 0)>0}
target_sets={r['score_set'] for r in read(AUD/'mavedb_bare_c_construct_by_set.csv')}
print(f'STAGE 1 done: {len(old_by_group)} groups, {len(eligible_sets)} eligible sets, {sum(len(x) for x in ref_groups.values())} ref-qualified accession groups, {sum(len(x) for x in strict_groups.values())} versioned-request accession groups',flush=True)

print('STAGE 2 stream target-sequence compatibility results',flush=True)
target_ids=defaultdict(set);target_raw=defaultdict(set);target_conflicts=Counter()
for r in rows(AUD/'mavedb_bare_c_construct_validation.csv'):
    gs=key_to_groups.get((r['score_set'],r['raw_hgvs_nt']),[])
    if not gs or r.get('coordinate_check_status')!='target_sequence_ref_and_codon_pass':continue
    aa=(str(r.get('dna_protein_pos','')),r.get('dna_ref_aa',''),r.get('dna_alt_aa',''))
    if r.get('p_annotation_status')=='p_annotation_conflict':
        target_conflicts.update(gs);continue
    for g in gs:
        if aa!=(g[2],g[3],g[4]):continue
        ident='SEQ:'+r.get('target_sequence_sha256','')+'|c.'+str(r.get('c_position',''))+str(r.get('c_ref',''))+'>'+str(r.get('c_alt',''))
        target_ids[g].add(ident);target_raw[g].add(r['raw_hgvs_nt'])
print(f'STAGE 2 done: {sum(map(len,target_ids.values()))} conditional target identities across {len(target_ids)} historical groups',flush=True)

print('STAGE 3 compare raw score exports with mapped rows for 27 eligible sets',flush=True)
keys_by_set=defaultdict(set);legacy_pairs_by_set=defaultdict(set);groups_by_set=defaultdict(set)
for g,rs in old_by_group.items():
    if g[0] not in eligible_sets:continue
    groups_by_set[g[0]].add(g)
    for r in rs:
        if r.get('key_category')!='protein_only':
            keys_by_set[g[0]].add(r['legacy_key'])
            legacy_pairs_by_set[g[0]].add((g[2],g[3],g[4],r['legacy_key']))
input_hashes={r['path']:r.get('sha256','') for r in rows(AUD/'input_sha256.csv')}
raw_hashes_by_name={Path(path).name:digest for path,digest in input_hashes.items()
                    if 'multigene_scores/' in path.replace('\\','/')}
source_metrics={};group_source_rows=Counter();group_source_scores=defaultdict(set)
for urn in sorted(eligible_sets):
    raw_path=RAW/(urn.replace(':','_')+'.csv');mapped_path=MAP/(urn.replace(':','_')+'.csv')
    if not raw_path.exists() or not mapped_path.exists():raise FileNotFoundError(urn)
    expected_raw_sha=raw_hashes_by_name.get(raw_path.name)
    if not expected_raw_sha:raise ValueError(f'Missing frozen source SHA256 for {urn}')
    if sha(raw_path)!=expected_raw_sha:raise ValueError(f'Frozen source SHA256 mismatch for {urn}')
    src_counter=Counter();src_keys=set();src_rows=0
    with open(raw_path,newline='',encoding='utf-8-sig') as f:
        rd=csv.DictReader(f);fields=rd.fieldnames or [];sc=score_col(fields)
        if sc is None or 'hgvs_nt' not in fields:raise ValueError(f'missing source score/hgvs_nt field: {urn}')
        for r in rd:
            x=number(r.get(sc))
            if x is None:continue
            nt=str(r.get('hgvs_nt','') or '')
            if nt not in keys_by_set[urn]:continue
            src_rows+=1;src_keys.add(nt)
            src_counter[(str(r.get('accession','') or ''),nt,str(r.get('hgvs_pro','') or ''),x)]+=1
    map_counter=Counter();mapped_rows=0;mapped_pairs=set()
    with open(mapped_path,newline='',encoding='utf-8-sig') as f:
        for r in csv.DictReader(f):
            if r.get('class')!='missense':continue
            nt=str(r.get('hgvs_nt','') or '')
            if nt not in keys_by_set[urn]:continue
            x=number(r.get('score'))
            if x is None:continue
            map_counter[(str(r.get('variant_id','') or ''),nt,str(r.get('hgvs_pro','') or ''),x)]+=1
            mapped_rows+=1
            g=(urn,set_gene.get(urn,''),str(r.get('protein_pos','')),r.get('ref_aa',''),r.get('alt_aa',''))
            if g in groups_by_set[urn]:
                mapped_pairs.add((g[2],g[3],g[4],nt));group_source_rows[g]+=1;group_source_scores[g].add(x)
    matched=sum((src_counter & map_counter).values())
    pair_found=len(legacy_pairs_by_set[urn] & mapped_pairs)
    source_metrics[urn]={
      'score_set':urn,'gene':set_gene.get(urn,''),'legacy_multikey_groups':len(groups_by_set[urn]),
      'legacy_author_nt_group_keys':len(legacy_pairs_by_set[urn]),'legacy_distinct_nt_strings':len(keys_by_set[urn]),
      'raw_numeric_rows_for_legacy_nt_keys':src_rows,'raw_distinct_nt_keys_found':len(src_keys),
      'raw_legacy_nt_keys_missing':len(keys_by_set[urn]-src_keys),
      'mapped_missense_rows_for_legacy_nt_keys':mapped_rows,'mapped_legacy_group_key_pairs_found':pair_found,
      'mapped_legacy_group_key_pairs_missing':len(legacy_pairs_by_set[urn])-pair_found,
      'raw_rows_reconciled_to_mapped_numeric_records':matched,'raw_numeric_rows_unmatched':sum(src_counter.values())-matched,
      'mapped_numeric_rows_without_raw_match':mapped_rows-matched,'exact_numeric_record_multiset_match':src_counter==map_counter,
      'raw_input_path':str(raw_path),'raw_input_sha256_from_frozen_manifest':expected_raw_sha,
      'mapped_input_path':str(mapped_path)
    }
    state='SOURCE OK' if src_counter==map_counter and pair_found==len(legacy_pairs_by_set[urn]) else 'SOURCE CHECK'
    print(f'{state} {urn}: groups={len(groups_by_set[urn])} group_keys={pair_found}/{len(legacy_pairs_by_set[urn])} numeric_rows={mapped_rows}',flush=True)
print('STAGE 3 done',flush=True)

# Precomputed strict/reference group files are complete for their stated transcript-qualified >=2 tier.
# If a group is absent, crosswalk key flags distinguish zero support; a one-key result is inferred only
# when the entire score set has one transcript accession. Multi-accession residuals stay unknown.
group_rows=[];per_set_groups=defaultdict(list)
for g,rs in sorted(old_by_group.items()):
    urn,gene,pos,ref,alt=g
    ref_raw=sum(yes(x.get('reference_supported_DNA')) for x in rs)
    strict_raw=sum(yes(x.get('versioned_identifier_request_returned_sequence')) for x in rs)
    rg=ref_groups.get(g,{});sg=strict_groups.get(g,{})
    rn=sum(len(v) for v in rg.values());sn=sum(len(v) for v in sg.values())
    rmax=max((len(v) for v in rg.values()),default=0);smax=max((len(v) for v in sg.values()),default=0)
    nacc=int(per.get(urn,{}).get('transcript_accession_count') or 0)
    def classify(nraw,ids_n,max_one):
        if max_one>=2:return ('ge2',ids_n if nacc==1 else None)
        if nraw==0:return ('zero',0)
        if nacc==1:return ('one',1)
        return ('unknown_0_or_1',None)
    rstate,rknown=classify(ref_raw,rn,rmax);sstate,sknown=classify(strict_raw,sn,smax)
    # A precomputed >=2 group row must be represented in the original key-level qualification flags.
    if rmax>=2 and ref_raw<rmax:raise ValueError(f'reference group/key crosswalk inconsistency: {g}')
    if smax>=2 and strict_raw<smax:raise ValueError(f'versioned group/key crosswalk inconsistency: {g}')
    t_ids=target_ids.get(g,set());t_raw=target_raw.get(g,set())
    if gene=='TPK1':tstate='tpk1_full_construct_unresolved'
    elif urn not in target_sets:tstate='not_assessed_target_sequence'
    elif len(t_ids)>=2:tstate='conditional_ge2'
    elif len(t_ids)==1:tstate='conditional_one'
    else:tstate='zero_target_compatible_under_assumption'
    categories=Counter(x.get('key_category','') for x in rs)
    # Group key provenance is checked by the complete raw/mapped per-set signature and separately counted above.
    row={
      'score_set':urn,'gene':gene,'protein_group':f'{pos}{ref}>{alt}','legacy_unique_keys':len(rs),
      'legacy_author_nt_score_keys':sum(x.get('key_category')!='protein_only' for x in rs),
      'legacy_protein_fallback_keys':sum(x.get('key_category')=='protein_only' for x in rs),
      'legacy_key_categories_json':json.dumps(dict(sorted(categories.items())),sort_keys=True),
      'numeric_mapped_score_rows_for_legacy_keys':group_source_rows.get(g,0),
      'distinct_numeric_score_values_diagnostic_only':len(group_source_scores.get(g,set())),
      'accession_reference_supported_raw_keys':ref_raw,
      'accession_reference_supported_unique_DNA_identities_known':rknown if rknown is not None else '',
      'accession_reference_group_state':rstate,
      'accession_reference_alias_keys_collapsed_known':ref_raw-rknown if rknown is not None else '',
      'versioned_request_supported_raw_keys':strict_raw,
      'versioned_request_supported_unique_DNA_identities_known':sknown if sknown is not None else '',
      'versioned_request_group_state':sstate,
      'versioned_request_alias_keys_collapsed_known':strict_raw-sknown if sknown is not None else '',
      'target_sequence_conditional_raw_keys':len(t_raw),
      'target_sequence_conditional_unique_DNA_identities':len(t_ids),
      'target_sequence_conditional_alias_keys_collapsed':len(t_raw)-len(t_ids),
      'target_sequence_group_state':tstate,'target_annotation_conflict_rows':target_conflicts.get(g,0),
      'empirical_independent_DNA_measurement_provenance':'not established by score/reference outputs',
      'qualified_reference_identity_ids_by_transcript_json':json.dumps({a:sorted(v) for a,v in rg.items()},sort_keys=True),
      'qualified_versioned_identity_ids_by_transcript_json':json.dumps({a:sorted(v) for a,v in sg.items()},sort_keys=True),
      'target_conditional_identity_ids_json':json.dumps(sorted(t_ids))
    }
    group_rows.append(row);per_set_groups[urn].append(row)
write(OUT/'group_level_empirical_testability.csv',group_rows)
write(OUT/'source_provenance_per_score_set.csv',list(source_metrics.values()))

per_set=[]
for urn in sorted(set(r['score_set'] for r in old_per)|set(per)):
    old=next((x for x in old_per if x['score_set']==urn),{});p=per.get(urn,{});gg=per_set_groups.get(urn,[])
    non_tpk=[x for x in gg if x['gene']!='TPK1']
    def ct(gg,field,state):return sum(x[field]==state for x in gg)
    per_set.append({
      'score_set':urn,'gene':p.get('gene',old.get('gene','')),'manifest_label_prefix':p.get('label_prefix',old.get('label_prefix','')),
      'legacy_missense_rows':int(old.get('legacy_missense_rows') or 0),
      'legacy_groups_ge2_unique_keys':int(old.get('legacy_groups_ge2_unique_keys') or 0),
      'legacy_unique_keys_in_ge2_groups':int(old.get('legacy_unique_keys_in_ge2_groups') or 0),
      'accession_reference_groups_ge2_all':ct(gg,'accession_reference_group_state','ge2'),
      'accession_reference_groups_ge2_non_TPK1':ct(non_tpk,'accession_reference_group_state','ge2'),
      'accession_reference_groups_dropping_to_one':ct(gg,'accession_reference_group_state','one'),
      'accession_reference_groups_zero_supported':ct(gg,'accession_reference_group_state','zero'),
      'accession_reference_groups_unknown_0_or_1':ct(gg,'accession_reference_group_state','unknown_0_or_1'),
      'versioned_request_groups_ge2_all':ct(gg,'versioned_request_group_state','ge2'),
      'versioned_request_groups_ge2_non_TPK1':ct(non_tpk,'versioned_request_group_state','ge2'),
      'versioned_request_groups_dropping_to_one':ct(gg,'versioned_request_group_state','one'),
      'versioned_request_groups_zero_supported':ct(gg,'versioned_request_group_state','zero'),
      'versioned_request_groups_unknown_0_or_1':ct(gg,'versioned_request_group_state','unknown_0_or_1'),
      'target_conditional_groups_ge2_non_TPK1':ct(non_tpk,'target_sequence_group_state','conditional_ge2'),
      'target_conditional_groups_one_non_TPK1':ct(non_tpk,'target_sequence_group_state','conditional_one'),
      'target_conditional_groups_zero_supported_non_TPK1':ct(non_tpk,'target_sequence_group_state','zero_target_compatible_under_assumption'),
      'target_conditional_groups_not_assessed':ct(gg,'target_sequence_group_state','not_assessed_target_sequence'),
      'TPK1_full_construct_unresolved_groups':ct(gg,'target_sequence_group_state','tpk1_full_construct_unresolved'),
      'accession_reference_raw_keys':sum(int(x['accession_reference_supported_raw_keys']) for x in gg),
      'accession_reference_unique_DNA_ids_known':sum(int(x['accession_reference_supported_unique_DNA_identities_known'] or 0) for x in gg),
      'accession_reference_alias_keys_collapsed_known':sum(int(x['accession_reference_alias_keys_collapsed_known'] or 0) for x in gg),
      'versioned_request_raw_keys':sum(int(x['versioned_request_supported_raw_keys']) for x in gg),
      'versioned_request_unique_DNA_ids_known':sum(int(x['versioned_request_supported_unique_DNA_identities_known'] or 0) for x in gg),
      'versioned_request_alias_keys_collapsed_known':sum(int(x['versioned_request_alias_keys_collapsed_known'] or 0) for x in gg),
      'non_reference_supported_legacy_keys':sum(int(x['legacy_unique_keys'])-int(x['accession_reference_supported_raw_keys']) for x in gg)
    })
write(OUT/'per_score_set_reconciliation.csv',per_set)

expected=json.load(open(AUD/'legacy_reproduction.json'))['observed']
cats=Counter(r.get('key_category','') for r in legacy)
dist_ref=Counter(r['accession_reference_group_state'] for r in group_rows)
dist_strict=Counter(r['versioned_request_group_state'] for r in group_rows)
src_ok=sum(bool(x['exact_numeric_record_multiset_match']) for x in source_metrics.values())
source_keypairs=sum(len(legacy_pairs_by_set[u]) for u in eligible_sets)
mapped_keypairs=sum(int(x['mapped_legacy_group_key_pairs_found']) for x in source_metrics.values())
all_mapped_keypairs=sum(len(legacy_pairs_by_set[u]) for u in eligible_sets)
obs={
 'historical_manifest_score_sets':92,
 'historical_analyzable_score_sets':sum(int(r.get('legacy_missense_rows') or 0)>0 for r in old_per),
 'historical_multi_key_groups':len(group_rows),'historical_keys':sum(int(r['legacy_unique_keys']) for r in group_rows),
 'eligible_score_sets_with_multi_key_groups':len(eligible_sets),
 'historical_key_categories':dict(sorted(cats.items())),
 'accession_reference_group_states_all':dict(sorted(dist_ref.items())),
 'versioned_request_group_states_all':dict(sorted(dist_strict.items())),
 'accession_reference_groups_ge2_non_TPK1':sum(r['gene']!='TPK1' and r['accession_reference_group_state']=='ge2' for r in group_rows),
 'versioned_request_groups_ge2_non_TPK1':sum(r['gene']!='TPK1' and r['versioned_request_group_state']=='ge2' for r in group_rows),
 'accession_reference_groups_one':dist_ref.get('one',0),'accession_reference_groups_zero':dist_ref.get('zero',0),
 'accession_reference_groups_unknown_0_or_1':dist_ref.get('unknown_0_or_1',0),
 'versioned_request_groups_one':dist_strict.get('one',0),'versioned_request_groups_zero':dist_strict.get('zero',0),
 'versioned_request_groups_unknown_0_or_1':dist_strict.get('unknown_0_or_1',0),
 'accession_reference_supported_raw_keys':sum(int(r['accession_reference_supported_raw_keys']) for r in group_rows),
 'accession_reference_unique_DNA_ids_known':sum(int(r['accession_reference_supported_unique_DNA_identities_known'] or 0) for r in group_rows),
 'accession_reference_alias_keys_collapsed_known':sum(int(r['accession_reference_alias_keys_collapsed_known'] or 0) for r in group_rows),
 'versioned_request_supported_raw_keys':sum(int(r['versioned_request_supported_raw_keys']) for r in group_rows),
 'versioned_request_unique_DNA_ids_known':sum(int(r['versioned_request_supported_unique_DNA_identities_known'] or 0) for r in group_rows),
 'versioned_request_alias_keys_collapsed_known':sum(int(r['versioned_request_alias_keys_collapsed_known'] or 0) for r in group_rows),
 'target_conditional_groups_ge2_non_TPK1':sum(r['gene']!='TPK1' and r['target_sequence_group_state']=='conditional_ge2' for r in group_rows),
 'target_conditional_unique_DNA_ids_all':sum(int(r['target_sequence_conditional_unique_DNA_identities']) for r in group_rows),
 'target_conditional_TPK1_ids_excluded_from_group_eligibility':sum(int(r['target_sequence_conditional_unique_DNA_identities']) for r in group_rows if r['gene']=='TPK1'),
 'TPK1_legacy_keys':sum(int(r['legacy_unique_keys']) for r in group_rows if r['gene']=='TPK1'),
 'TPK1_historical_unresolved_keys_held_unqualified':sum(int(r['legacy_unique_keys']) for r in group_rows if r['gene']=='TPK1'),
 'source_sets_checked':len(source_metrics),'source_sets_exact_numeric_multiset_match':src_ok,
 'legacy_group_key_pairs_in_27_sets':source_keypairs,'mapped_group_key_pairs_found':mapped_keypairs,
 'mapped_group_key_pairs_missing':all_mapped_keypairs-mapped_keypairs,
 'source_to_mapped_numeric_rows_reconciled':sum(int(x['raw_rows_reconciled_to_mapped_numeric_records']) for x in source_metrics.values()),
 'source_numeric_rows_unmatched':sum(int(x['raw_numeric_rows_unmatched']) for x in source_metrics.values()),
 'mapped_rows_without_raw_source_match':sum(int(x['mapped_numeric_rows_without_raw_match']) for x in source_metrics.values())
}
checks={
 'historical_groups_equal_expected':obs['historical_multi_key_groups']==expected['groups_ge2'],
 'historical_keys_equal_expected':obs['historical_keys']==expected['unique_keys_in_ge2_groups'],
 'score_sets_with_legacy_group_equal_27':len(eligible_sets)==27,
 'legacy_key_crosswalk_rows_equal_expected':len(legacy)==expected['unique_keys_in_ge2_groups'],
 'source_nt_rows_cover_all_historical_group_key_pairs':mapped_keypairs==source_keypairs and mapped_keypairs==all_mapped_keypairs,
 'raw_input_and_mapped_numeric_records_match_for_all_27':src_ok==27 and obs['source_numeric_rows_unmatched']==0 and obs['mapped_rows_without_raw_source_match']==0
}
hashes=[]
for p in [AUD/'legacy_reproduction.json',AUD/'legacy_per_set.csv',AUD/'legacy_key_crosswalk.csv',AUD/'reference_supported_groups.csv',AUD/'strict_groups.csv',AUD/'per_score_set.csv',AUD/'mavedb_bare_c_construct_validation.csv',AUD/'mavedb_bare_c_construct_methods.md',AUD/'input_sha256.csv',Path(__file__).resolve().parent/'legacy_mapping/map_all_scores.py']:
    hashes.append({'path':str(p),'sha256':sha(p),'size_bytes':p.stat().st_size})
report={
 'run_at_utc':datetime.now(timezone.utc).isoformat(),
 'source_tier_definitions':{
  'accession_reference_supported':'Submitted hgvs_nt SNV with pinned CDS hash, reference match, DNA-derived translation, and supplied p annotation absent or consistent; versioned-request status is separate.',
  'versioned_request_supported':'Narrower archived sequence-request tier; returned accession/version identity was not independently verified because FASTA deflines were not preserved.',
  'target_sequence_conditional':'Bare c.HGVS SNV matched supplied MaveDB targetSequence under the explicit, unverified c.1-to-target-base-1 assumption. Kept separate from accession-qualified and empirical-measurement counts.',
  'reported_nucleotide_resolved_score_record':'Original incoming score export carries hgvs_nt and finite numeric score. map_all_scores.py preserves that raw field and has no protein-to-DNA reverse-translation path. This does not prove physically independent DNA measurements.'
 },
 'historical_expected':expected,'observed':obs,'reconciliation_checks':checks,'input_output_hashes':hashes,
 'precomputed_group_list_validation':'The frozen audit script builds strict_groups.csv and reference_supported_groups.csv from every qualifying canonical identity, stratified by score set and transcript accession, and writes only groups with >=2 unique identities. The absence of a group is therefore a complete >=2 absence for that tier; 0 versus 1 is called exact only when key-level flags show zero or the score set has exactly one transcript accession. Other residuals are unknown.',
 'independent_measurement_boundary':'Independent empirical DNA measurement provenance is not identified in these files. Reported nucleotide score records are not treated as proof of separate clones/readouts; distinct numeric scores are diagnostics only.',
 'tpk1_boundary':'The historical 6,913 unresolved TPK1 keys remain unqualified for full-construct group eligibility. Any local target-sequence compatibility is shown conditionally and excluded.',
 'historic_support_count_boundary':'The previous 24,810/31,723 sequence-support tally validates only its historical numerator under mixed reference tiers; it is not revised collection coverage. The present output keeps accession and conditional target tiers separate.'
}
checks['all_reconciliation_checks_pass']=all(checks.values())
report['reconciliation_checks']=checks
(OUT/'collection_group_reconciliation.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
print('STAGE 4 outputs written',flush=True)
print(json.dumps(checks,sort_keys=True),flush=True)
print(json.dumps(obs,sort_keys=True),flush=True)


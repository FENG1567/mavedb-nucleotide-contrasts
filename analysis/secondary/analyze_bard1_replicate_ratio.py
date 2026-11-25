#!/usr/bin/env python3
"""BARD1 same-target/same-AA replicate-ratio analysis under frozen protocol."""
from pathlib import Path
from datetime import datetime, timezone
import argparse, hashlib, json, re
import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT_DIR = REPO_ROOT / 'data' / 'external' / 'bard1'
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--legacy-counts', type=Path, default=DEFAULT_INPUT_DIR / '20250922_BARD1RNAscores.tsv',
                    help='pinned legacy RNA/D05/D13 count table')
parser.add_argument('--annotation', type=Path, default=DEFAULT_INPUT_DIR / '20260807_BARD1.allscores.tsv',
                    help='same-commit annotation table')
parser.add_argument('--current-counts', type=Path, default=DEFAULT_INPUT_DIR / '20260807_BARD1.snvcounts.tsv',
                    help='newer count table used only for D05 version comparison')
parser.add_argument('--output-dir', type=Path, default=REPO_ROOT / 'results' / 'secondary' / 'bard1',
                    help='directory for generated analysis outputs')
args = parser.parse_args()

LEGACY = args.legacy_counts.resolve()
ANNOT = args.annotation.resolve()
CURRENT = args.current_counts.resolve()
OUT = args.output_dir.resolve()
OUT.mkdir(parents=True, exist_ok=True)
COMMIT = '72ba50c8480a8dd88528c811ccc906d41fb3d971'
BOOTSTRAPS = 2000
BOOTSTRAP_SEED = 20261122



def sha256(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):
            h.update(block)
    return h.hexdigest()


def normalize(frame):
    frame=frame.copy()
    frame['chrom']=frame['chrom'].astype(str).str.replace('chr','',regex=False)
    frame['pos']=pd.to_numeric(frame['pos'],errors='coerce').astype('Int64')
    frame['ref']=frame['ref'].astype(str).str.upper()
    frame['alt']=frame['alt'].astype(str).str.upper()
    return frame


def float_or_none(x):
    return None if x is None or not np.isfinite(float(x)) else float(x)


def position_from_aa(value):
    m=re.search(r'(\d+)',str(value))
    return int(m.group(1)) if m else None

# Read only the legacy raw count table for RNA/D05/D13 analysis. The newer count
# file is read only below to document version discordance; no new-version counts
# enter any ratio calculation.
legacy=pd.read_csv(LEGACY,sep='\t',low_memory=False)
annot=pd.read_csv(ANNOT,sep='\t',low_memory=False)
current=pd.read_csv(CURRENT,sep='\t',low_memory=False)
legacy=normalize(legacy); annot=normalize(annot); current=normalize(current)
variant_keys=['chrom','pos','ref','alt']
full_keys=variant_keys+['target']
if legacy.duplicated(full_keys).any():
    raise SystemExit('Legacy target-specific variant key is not unique.')
if annot.duplicated(variant_keys).any():
    raise SystemExit('Annotation genomic variant key is not unique; mapping ambiguous.')

annotation_columns=variant_keys+['target','consequence','amino_acid_change','hgvs_p','hgvs_c','CDS_position','cDNA_position']
ann=annot[annotation_columns].rename(columns={
    'target':'target_annotation','consequence':'consequence_annotation'
})
data=legacy.merge(ann,on=variant_keys,how='left',indicator=True,validate='many_to_one')
data['annotation_mapped']=data['_merge'].eq('both')
data['dna_key']=data.apply(lambda r:f"{r['chrom']}:{int(r['pos'])}:{r['ref']}>{r['alt']}" if pd.notna(r['pos']) else '',axis=1)

count_cols=['D05_R1','D05_R2','D05_R3','D13_R1','D13_R2','D13_R3','rna_count_R1','rna_count_R2','rna_count_R3']
missing_cols=[c for c in count_cols if c not in data.columns]
if missing_cols: raise SystemExit(f'Missing legacy count columns: {missing_cols}')
for c in count_cols: data[c]=pd.to_numeric(data[c],errors='coerce')

# Exact version comparison: same genomic SNV and measured target, for D05 only.
version_cmp=legacy[full_keys+['D05_R1','D05_R2','D05_R3']].merge(
    current[full_keys+['D05_R1','D05_R2','D05_R3']],on=full_keys,how='inner',suffixes=('_legacy','_current'),validate='one_to_one')
version_rows=[]
for c in ['D05_R1','D05_R2','D05_R3']:
    a=pd.to_numeric(version_cmp[c+'_legacy'],errors='coerce'); b=pd.to_numeric(version_cmp[c+'_current'],errors='coerce'); ok=a.notna()&b.notna()
    version_rows.append({'field':c,'n_overlap_nonmissing':int(ok.sum()),'n_equal':int((a[ok]==b[ok]).sum()),'n_different':int((a[ok]!=b[ok]).sum())})
pd.DataFrame(version_rows).to_csv(OUT/'legacy_vs_current_D05_comparison.csv',index=False)

# Eligibility: annotation-backed missense, all three D05 baseline counts >= 10,
# then at least two distinct genomic SNVs within each measured target x AA group.
missense=data['consequence_annotation'].astype(str).str.contains('missense',case=False,na=False)
missense &= data['amino_acid_change'].notna() & data['amino_acid_change'].astype(str).str.strip().ne('')
d05_values=data[['D05_R1','D05_R2','D05_R3']].to_numpy(dtype=float)
d05_complete=np.isfinite(d05_values).all(axis=1)
d05_ge10=d05_complete & (d05_values>=10).all(axis=1)
candidates=data.loc[missense.to_numpy() & d05_ge10].copy()
candidates['AA_group']=candidates['amino_acid_change'].astype(str).str.strip()
candidates['AA_position']=candidates['AA_group'].map(position_from_aa)
candidates['target_group']=candidates['target'].astype(str)+'||'+candidates['AA_group']
if candidates['AA_position'].isna().any():
    raise SystemExit('Some eligible missense annotations have no protein position for position-cluster bootstrap.')

# Unsplit same-AA inventory is for reporting only; statistical groups include target.
unsplit=(candidates.groupby('AA_group',as_index=False)
         .agg(n_distinct_DNA=('dna_key','nunique'),n_legacy_targets=('target','nunique'),
              legacy_targets=('target',lambda x:';'.join(sorted(set(map(str,x)))))))
unsplit_ge2=unsplit[unsplit['n_distinct_DNA']>=2].copy()
unsplit_ge2['spans_multiple_targets']=unsplit_ge2['n_legacy_targets']>1

group_eligibility=(candidates.groupby(['target','AA_group','AA_position'],as_index=False)
                   .agg(n_distinct_DNA=('dna_key','nunique'),n_rows=('dna_key','size')))
group_eligibility=group_eligibility[group_eligibility['n_distinct_DNA']>=2].copy()
group_eligibility['target_group']=group_eligibility['target'].astype(str)+'||'+group_eligibility['AA_group'].astype(str)
# Set the group key on the full annotation/count table consistently before selection.
data['AA_group']=data['amino_acid_change'].astype(str).str.strip()
data['AA_position']=data['AA_group'].map(position_from_aa)
data['target_group']=data['target'].astype(str)+'||'+data['AA_group'].astype(str)
eligible=data.loc[data['target_group'].isin(set(group_eligibility['target_group'])) & missense.to_numpy() & d05_ge10].copy()
eligible=eligible.sort_values(['target','AA_group','dna_key'],kind='stable').reset_index(drop=True)
if eligible.duplicated(full_keys).any(): raise SystemExit('Duplicate target-specific rows survived into eligible data.')

# The exact DNA allele can be measured in more than one amplicon. These are
# repeated assay-context rows, not new DNA alleles; they remain in target-specific
# groups and are clustered by protein position for uncertainty.
key_target_counts=eligible.groupby('dna_key')['target'].nunique()
multitarget_dna_keys=int((key_target_counts>1).sum())
extra_target_rows=int((key_target_counts-1).clip(lower=0).sum())

# Count missing/zero values without filtering on RNA or D13 observations.
miss_rows=[]
for c in count_cols:
    x=data[c]; y=eligible[c]
    miss_rows.append({'count_field':c,'legacy_rows':int(len(data)),'legacy_missing':int(x.isna().sum()),'legacy_zero':int((x==0).sum()),
                      'eligible_rows':int(len(eligible)),'eligible_missing':int(y.isna().sum()),'eligible_zero':int((y==0).sum())})
pd.DataFrame(miss_rows).to_csv(OUT/'count_missingness_and_zeros.csv',index=False)

# Source annotation/coverage summary, including target membership in the newer
# semicolon-delimited annotation target labels.
def target_in_annotation(row):
    if pd.isna(row['target_annotation']): return False
    return str(row['target']) in set(str(row['target_annotation']).split(';'))
annotation_mapped=data['annotation_mapped']
target_match=int(data.loc[annotation_mapped].apply(target_in_annotation,axis=1).sum())

# Construct fixed +1 raw-count log ratios (zero counts remain valid).
for i in (1,2,3):
    data[f'R{i}']=np.log2(data[f'rna_count_R{i}']+1)-np.log2(data[f'D05_R{i}']+1)
    data[f'F{i}']=np.log2(data[f'D13_R{i}']+1)-np.log2(data[f'D05_R{i}']+1)
for i in (1,2,3):
    eligible[f'R{i}']=np.log2(eligible[f'rna_count_R{i}']+1)-np.log2(eligible[f'D05_R{i}']+1)
    eligible[f'F{i}']=np.log2(eligible[f'D13_R{i}']+1)-np.log2(eligible[f'D05_R{i}']+1)

ratio_cols=[f'{p}{i}' for p in ['R','F'] for i in (1,2,3)]
export_cols=['target','target_annotation','AA_group','AA_position','dna_key','hgvs_p','hgvs_c','consequence_annotation',
             'D05_R1','D05_R2','D05_R3','rna_count_R1','rna_count_R2','rna_count_R3','D13_R1','D13_R2','D13_R3']+ratio_cols
eligible[[c for c in export_cols if c in eligible.columns]].to_csv(OUT/'eligible_variants_and_ratios.csv',index=False)
group_eligibility.sort_values(['AA_position','target','AA_group']).to_csv(OUT/'target_AA_group_eligibility.csv',index=False)
unsplit_ge2.sort_values('AA_group').to_csv(OUT/'unsplit_AA_target_spans.csv',index=False)

# Compute within-group sample moments once for each direction and group.
directions=[(1,1,'within'),(2,2,'within'),(3,3,'within'),
            (1,2,'cross'),(1,3,'cross'),(2,1,'cross'),(2,3,'cross'),(3,1,'cross'),(3,2,'cross')]
group_rows=[]; pair_rows=[]
for (target,aa),g in eligible.groupby(['target','AA_group'],sort=True):
    g=g.sort_values('dna_key',kind='stable').reset_index(drop=True)
    pos=int(g['AA_position'].iloc[0])
    first=g.iloc[0]; second=g.iloc[1]
    for ri,fj,kind in directions:
        x=pd.to_numeric(g[f'R{ri}'],errors='coerce').to_numpy(dtype=float)
        y=pd.to_numeric(g[f'F{fj}'],errors='coerce').to_numpy(dtype=float)
        ok=np.isfinite(x)&np.isfinite(y); x=x[ok]; y=y[ok]
        n=int(len(x))
        if n>=2:
            xc=x-x.mean(); yc=y-y.mean()
            cov=float(np.dot(xc,yc)/(n-1)); varx=float(np.dot(xc,xc)/(n-1)); vary=float(np.dot(yc,yc)/(n-1))
            rg=cov/float(np.sqrt(varx*vary)) if varx>0 and vary>0 else np.nan
        else:
            cov=varx=vary=rg=np.nan
        group_rows.append({'target':target,'AA_group':aa,'AA_position':pos,'direction':f'R{ri}:F{fj}',
                           'replicate_relation':kind,'n_complete_variants':n,
                           'sample_covariance':cov,'sample_variance_R':varx,'sample_variance_F':vary,'group_centered_r':rg})
        x1=float(first[f'R{ri}']) if pd.notna(first[f'R{ri}']) else np.nan
        x2=float(second[f'R{ri}']) if pd.notna(second[f'R{ri}']) else np.nan
        y1=float(first[f'F{fj}']) if pd.notna(first[f'F{fj}']) else np.nan
        y2=float(second[f'F{fj}']) if pd.notna(second[f'F{fj}']) else np.nan
        if not np.isfinite([x1,x2,y1,y2]).all(): status='missing'
        else:
            dx=x2-x1; dy=y2-y1
            status='concordant' if dx*dy>0 else 'discordant' if dx*dy<0 else 'tie'
        pair_rows.append({'target':target,'AA_group':aa,'AA_position':pos,'direction':f'R{ri}:F{fj}',
                          'variant_1':first['dna_key'],'variant_2':second['dna_key'],
                          'R_1':x1,'R_2':x2,'delta_R_2_minus_1':x2-x1 if np.isfinite([x1,x2]).all() else np.nan,
                          'F_1':y1,'F_2':y2,'delta_F_2_minus_1':y2-y1 if np.isfinite([y1,y2]).all() else np.nan,
                          'first_pair_status':status})

group_stats=pd.DataFrame(group_rows)
pairs=pd.DataFrame(pair_rows)
group_stats.to_csv(OUT/'group_direction_moments.csv',index=False)
pairs.to_csv(OUT/'first_two_DNA_lexicographic_signs.csv',index=False)

# Position-cluster bootstrap using cached per-group sample moments. Each cluster
# contributes all target-AA groups at a protein position. A multinomial weight
# matrix resamples positions; no per-draw pandas or covariance calls are made.
rng=np.random.default_rng(BOOTSTRAP_SEED)
summary=[]; boot_rows=[]
for ri,fj,kind in directions:
    direction=f'R{ri}:F{fj}'
    gs=group_stats[(group_stats['direction']==direction)&group_stats['sample_covariance'].notna()].copy()
    if gs.empty:
        summary.append({'direction':direction,'replicate_relation':kind,'n_groups':0,'n_target_variant_observations':0,'n_positions':0,
                        'mean_group_sample_covariance':np.nan,'group_equal_centered_r':np.nan,
                        'cov_ci_lower':np.nan,'cov_ci_upper':np.nan,'r_ci_lower':np.nan,'r_ci_upper':np.nan,
                        'bootstrap_valid_draws':0,'bootstrap_requested_draws':BOOTSTRAPS})
        continue
    point_cov=float(gs['sample_covariance'].mean())
    sumxy=float(gs['sample_covariance'].sum()); sumxx=float(gs['sample_variance_R'].sum()); sumyy=float(gs['sample_variance_F'].sum())
    point_r=sumxy/float(np.sqrt(sumxx*sumyy)) if sumxx>0 and sumyy>0 else np.nan
    pos_agg=gs.groupby('AA_position').agg(n_groups=('sample_covariance','size'),sum_cov=('sample_covariance','sum'),
                                          sum_varR=('sample_variance_R','sum'),sum_varF=('sample_variance_F','sum')).sort_index()
    pos_ids=pos_agg.index.to_numpy(dtype=int); k=len(pos_ids)
    W=rng.multinomial(k,np.full(k,1.0/k),size=BOOTSTRAPS).astype(float)
    n_boot=W@pos_agg['n_groups'].to_numpy(float)
    cov_num=W@pos_agg['sum_cov'].to_numpy(float)
    varR_num=W@pos_agg['sum_varR'].to_numpy(float)
    varF_num=W@pos_agg['sum_varF'].to_numpy(float)
    boot_cov=np.divide(cov_num,n_boot,out=np.full(BOOTSTRAPS,np.nan),where=n_boot>0)
    boot_r=np.divide(cov_num,np.sqrt(varR_num*varF_num),out=np.full(BOOTSTRAPS,np.nan),where=(varR_num>0)&(varF_num>0))
    valid_cov=np.isfinite(boot_cov); valid_r=np.isfinite(boot_r)
    cov_lo,cov_hi=np.nanpercentile(boot_cov[valid_cov],[2.5,97.5]) if valid_cov.any() else (np.nan,np.nan)
    r_lo,r_hi=np.nanpercentile(boot_r[valid_r],[2.5,97.5]) if valid_r.any() else (np.nan,np.nan)
    # Descriptive first-pair sign concordance; non-ties form the rate denominator.
    ps=pairs[pairs.direction==direction]
    ps_valid=ps[ps.first_pair_status.isin(['concordant','discordant','tie'])]
    pos_pair=ps_valid.groupby('AA_position').agg(n_pairs=('first_pair_status','size'),
        n_concordant=('first_pair_status',lambda x:int((x=='concordant').sum())),
        n_discordant=('first_pair_status',lambda x:int((x=='discordant').sum())),
        n_tie=('first_pair_status',lambda x:int((x=='tie').sum()))).sort_index()
    pair_ids=pos_pair.index.to_numpy(dtype=int)
    if len(pair_ids):
        Wp=rng.multinomial(len(pair_ids),np.full(len(pair_ids),1.0/len(pair_ids)),size=BOOTSTRAPS).astype(float)
        conc=Wp@pos_pair['n_concordant'].to_numpy(float)
        discord=Wp@pos_pair['n_discordant'].to_numpy(float)
        nontie=conc+discord
        boot_rate=np.divide(conc,nontie,out=np.full(BOOTSTRAPS,np.nan),where=nontie>0)
        valid_rate=np.isfinite(boot_rate)
        rate_lo,rate_hi=np.nanpercentile(boot_rate[valid_rate],[2.5,97.5]) if valid_rate.any() else (np.nan,np.nan)
    else:
        rate_lo=rate_hi=np.nan;valid_rate=np.zeros(BOOTSTRAPS,dtype=bool);boot_rate=np.full(BOOTSTRAPS,np.nan)
    concordant=int((ps.first_pair_status=='concordant').sum())
    discordant=int((ps.first_pair_status=='discordant').sum())
    ties=int((ps.first_pair_status=='tie').sum())
    pair_missing=int((ps.first_pair_status=='missing').sum())
    non_tie=concordant+discordant
    summary.append({'direction':direction,'replicate_relation':kind,'n_groups':int(len(gs)),
                    'n_target_variant_observations':int(gs['n_complete_variants'].sum()),'n_positions':int(k),
                    'mean_group_sample_covariance':point_cov,'group_equal_centered_r':point_r,
                    'cov_ci_lower':float(cov_lo),'cov_ci_upper':float(cov_hi),
                    'r_ci_lower':float(r_lo),'r_ci_upper':float(r_hi),
                    'bootstrap_valid_cov_draws':int(valid_cov.sum()),'bootstrap_valid_r_draws':int(valid_r.sum()),
                    'bootstrap_requested_draws':BOOTSTRAPS,
                    'first2_groups_concordant':concordant,'first2_groups_discordant':discordant,
                    'first2_groups_tie':ties,'first2_groups_missing':pair_missing,
                    'first2_concordant_fraction_non_ties':concordant/non_tie if non_tie else np.nan,
                    'first2_concordance_ci_lower':float(rate_lo),'first2_concordance_ci_upper':float(rate_hi),
                    'first2_bootstrap_valid_draws':int(valid_rate.sum())})
    for b in range(BOOTSTRAPS):
        boot_rows.append({'direction':direction,'bootstrap_draw':b+1,'mean_group_sample_covariance':boot_cov[b],
                          'group_equal_centered_r':boot_r[b],
                          'first2_concordant_fraction_non_ties':boot_rate[b]})
summary_df=pd.DataFrame(summary)
summary_df.to_csv(OUT/'direction_summary_and_position_cluster_CI.csv',index=False)
# Keep complete bootstrap draws for reproducibility/audit; these are summaries only.
pd.DataFrame(boot_rows).to_csv(OUT/'position_cluster_bootstrap_draws.csv',index=False)

# Matched count source/version and field availability details.
legacy_sha=sha256(LEGACY); annot_sha=sha256(ANNOT); current_sha=sha256(CURRENT)
legacy_count_field_summary=pd.DataFrame(miss_rows)
route={
 'gene':'BARD1','route_status':'viable_target_AA_group_comparison' if len(group_eligibility)>0 else 'stop_no_eligible_group',
 'pinned_repository_commit':COMMIT,'legacy_table_rows':int(len(legacy)),
 'legacy_unique_target_variant_rows':int(legacy[full_keys].drop_duplicates().shape[0]),
 'legacy_unique_genomic_variants':int(legacy[variant_keys].drop_duplicates().shape[0]),
 'newer_core_snvcount_rows':int(len(current)),'annotation_rows':int(len(annot)),
 'legacy_rows_with_unique_genomic_annotation':int(annotation_mapped.sum()),
 'legacy_rows_without_annotation':int((~annotation_mapped).sum()),
 'legacy_rows_with_HGVS_p':int(data['hgvs_p'].notna().sum()),'legacy_rows_with_HGVS_c':int(data['hgvs_c'].notna().sum()),
 'legacy_rows_with_both_HGVS_p_c':int((data['hgvs_p'].notna()&data['hgvs_c'].notna()).sum()),
 'legacy_rows_with_source_target':int(legacy['target'].notna().sum()),'legacy_target_annotation_membership_matches':target_match,
 'legacy_rows_baseline_D05_all_reps_ge10':int(d05_ge10.sum()),
 'annotated_missense_rows_baseline_eligible':int(len(candidates)),
 'unsplit_AA_groups_ge2_DNA':int(len(unsplit_ge2)),
 'unsplit_AA_groups_spanning_targets':int(unsplit_ge2['spans_multiple_targets'].sum()),
 'target_AA_groups_ge2_DNA':int(len(group_eligibility)),
 'eligible_target_AA_variant_rows':int(len(eligible)),
 'eligible_unique_DNA_variants':int(eligible['dna_key'].nunique()),
 'eligible_DNA_variants_measured_in_multiple_targets':multitarget_dna_keys,
 'eligible_extra_target_observation_rows':extra_target_rows,
 'legacy_sha256':legacy_sha,'annotation_sha256':annot_sha,'current_snvcounts_sha256_for_version_comparison_only':current_sha,
 'bootstrap_draws':BOOTSTRAPS,'bootstrap_seed':BOOTSTRAP_SEED,
 'analysis_completed_utc':datetime.now(timezone.utc).isoformat()}
(OUT/'route_status.json').write_text(json.dumps(route,indent=2)+'\n',encoding='utf-8')
pd.DataFrame([route]).to_csv(OUT/'route_status.csv',index=False)

# Compact report. Exact means/intervals are descriptive summaries of ratio
# associations across target-AA groups, with position-cluster bootstrap CIs.
version_diff_text=', '.join(f"{r['field']} {r['n_different']:,}/{r['n_overlap_nonmissing']:,}" for r in version_rows)
lines=['# BARD1 same-amino-acid cross-replicate ratio analysis','',
       '## Source and eligibility','',
       f'- Count source: `20250922_BARD1RNAscores.tsv`, pinned repository commit `{COMMIT}`; all RNA, D05, and D13 counts come from this one legacy version.',
       f'- Legacy table rows: {len(legacy):,}; newer core SNV count rows: {len(current):,}. The legacy table has {int((~annotation_mapped).sum()):,} rows without a unique annotation match in the same-commit final annotation table.',
       f'- Legacy table has `{int(legacy["target"].notna().sum()):,}` measured-target labels, `pos_id` and `cdna_pos`, but no HGVS columns. Same-commit annotation maps {int(annotation_mapped.sum()):,} rows; {int(data["hgvs_p"].notna().sum()):,} have HGVS.p, {int(data["hgvs_c"].notna().sum()):,} have HGVS.c, and {int((data["hgvs_p"].notna()&data["hgvs_c"].notna()).sum()):,} have both.',
       f'- The legacy cohort has {int(d05_ge10.sum()):,} rows with all three D05 counts present and at least 10. It contributes {int(len(candidates)):,} annotated missense rows; {int(len(unsplit_ge2)):,} unsplit amino-acid substitution groups contain at least two distinct DNA variants, {int(unsplit_ge2["spans_multiple_targets"].sum()):,} of those span multiple legacy targets.',
       f'- For count-ratio centering, groups are split by measured `(target, amino-acid substitution)`: {int(len(group_eligibility)):,} eligible groups and {int(len(eligible)):,} target-specific variant rows. {multitarget_dna_keys:,} exact DNA variants occur in more than one eligible target, adding {extra_target_rows:,} repeated target observations; these are assay-context repeats, not new DNA alleles.',
       '- Newer annotation target strings may list overlapping amplicons separated by semicolons; the legacy row target is the measured amplicon used for grouping. The paper methods describe variant counts per target, replicate, and time point; its resource statement links the analysis repo, IGVF raw reads, and MaveDB final scores.',
       '- No RNA score, functional score, RNA depth, or D13 depth filter was applied. Zero RNA and D13 counts were retained via +1; missing counts remain missing and are tabulated separately.',
       '', '## Definitions','',
       '- For RNA replicate i, `R_i = log2(rna_count_Ri + 1) - log2(D05_Ri + 1)`.',
       '- For function replicate j, `F_j = log2(D13_Rj + 1) - log2(D05_Rj + 1)`.',
       '- Covariance is the arithmetic mean of within-group sample covariances (`ddof=1`), equal-weighting each target-AA group. Centered r is `mean(group covariance) / sqrt(mean(group R variance) * mean(group F variance))`.',
       '- 95% percentile intervals use 2,000 protein-position cluster-bootstrap draws, carrying all target-AA groups at a sampled position; seed 20261122. No p-values are reported.',
       '', '## All nine replicate directions','',
       '| Direction | Relation | Groups | Target-specific variant observations | Positions | Mean sample covariance (95% CI) | Group-equal centered r (95% CI) | First-pair concordant / discordant / tie / missing | First-pair concordance fraction among non-ties (95% CI) |',
       '|---|---|---:|---:|---:|---:|---:|---:|---:|']
for _,r in summary_df.iterrows():
    ci_cov=f"{r['mean_group_sample_covariance']:.5g} ({r['cov_ci_lower']:.5g}, {r['cov_ci_upper']:.5g})" if pd.notna(r['mean_group_sample_covariance']) else 'NA'
    ci_r=f"{r['group_equal_centered_r']:.5g} ({r['r_ci_lower']:.5g}, {r['r_ci_upper']:.5g})" if pd.notna(r['group_equal_centered_r']) else 'NA'
    p=f"{int(r['first2_groups_concordant'])} / {int(r['first2_groups_discordant'])} / {int(r['first2_groups_tie'])} / {int(r['first2_groups_missing'])}"
    ci_pair=f"{r['first2_concordant_fraction_non_ties']:.4f} ({r['first2_concordance_ci_lower']:.4f}, {r['first2_concordance_ci_upper']:.4f})" if pd.notna(r['first2_concordant_fraction_non_ties']) else 'NA'
    lines.append(f"| {r['direction']} | {r['replicate_relation']} | {int(r['n_groups'])} | {int(r['n_target_variant_observations'])} | {int(r['n_positions'])} | {ci_cov} | {ci_r} | {p} | {ci_pair} |")
lines += ['', 'First-pair signs use the two lexicographically first exact DNA keys in each target-AA group. These group-level concordance counts are descriptive; the bootstrap unit for association uncertainty is protein position.',
          '', '## Version and measurement-context checks','',
          f'- Legacy versus newer core D05 count mismatches on exact target-specific overlapping rows: {version_diff_text} (see `legacy_vs_current_D05_comparison.csv`). All RNA and function ratios therefore use the legacy table only.',
          f'- Cross-target amino-acid substitutions were split into target-specific groups before centering; position bootstrap resamples all target-AA groups at a residue together. This handles amplicon-specific library totals and avoids treating repeated variants across amplicons as new DNA alleles.',
          '- The legacy table contains fewer rows than the newer core SNV count table, and the author repository does not specify why the legacy RNA table has this subset. Results therefore describe eligible missense SNVs represented in that legacy table, not all assayed BARD1 variants.',
          '- The nine directions are related comparisons, not nine independent validations; the results are descriptive and do not establish causality or mechanism.',
          '', '## Files and reproducibility','',
          '- `eligible_variants_and_ratios.csv`: selected variant rows, count inputs, and all nine log-ratio coordinates.',
          '- `target_AA_group_eligibility.csv`: common group eligibility and DNA-variant counts.',
          '- `group_direction_moments.csv`: cached within-group moments for every direction.',
          '- `direction_summary_and_position_cluster_CI.csv`: all nine estimates and position-cluster intervals.',
          '- `first_two_DNA_lexicographic_signs.csv`: group-level first-pair records for all directions.',
          '- `position_cluster_bootstrap_draws.csv`: bootstrap estimates for reproducibility.',
          '- `count_missingness_and_zeros.csv`, `legacy_vs_current_D05_comparison.csv`, `route_status.csv`, and `input_manifest.csv`: input checks and source accounting.',
          f'- SHA256 legacy count table: `{legacy_sha}`; annotation table: `{annot_sha}`; current count table used only for comparison: `{current_sha}`.',
          f'- Pinned code repository: `https://github.com/ivanw314/BARD1_SGE_analysis/tree/{COMMIT}`.']
(OUT/'report.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')

manifest=[]
for p,usage in [(LEGACY,'only RNA/D05/D13 count source'),(ANNOT,'annotation only; score values unused'),(CURRENT,'version-comparison count table only; not used in ratios')]:
    manifest.append({'file':str(p),'bytes':p.stat().st_size,'sha256':sha256(p),'repository_commit':COMMIT,'usage':usage})
pd.DataFrame(manifest).to_csv(OUT/'input_manifest.csv',index=False)
print(json.dumps(route,indent=2))
print(summary_df.to_string(index=False))
print('outputs='+str(OUT))

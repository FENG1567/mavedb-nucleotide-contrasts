import argparse, json, hashlib
from collections import defaultdict, Counter
from pathlib import Path
import numpy as np
import pandas as pd

REPO_ROOT=Path(__file__).resolve().parents[1]
DATA_ROOT=REPO_ROOT/'data/derived'
EXTERNAL_ROOT=REPO_ROOT/'data/external/brca2_cc0'
parser=argparse.ArgumentParser(description='Recompute BRCA2 direction/tie estimators from the two mapped score tables.')
parser.add_argument('--hap1',type=Path,default=EXTERNAL_ROOT/'mapped_BRCA2_00001225-a-1.csv')
parser.add_argument('--mes',type=Path,default=EXTERNAL_ROOT/'mapped_BRCA2_00001242-a-1.csv')
parser.add_argument('--phasea-summary',type=Path,default=DATA_ROOT/'work/phaseA_summary.json')
parser.add_argument('--holdout-groups',type=Path,default=DATA_ROOT/'work/holdout_paired_groups.csv')
parser.add_argument('--covariance-summary',type=Path,default=DATA_ROOT/'revision_20261122/statistics/brca2_covariance_summary.json')
parser.add_argument('--output-dir',type=Path,default=REPO_ROOT/'results/brca2_direction')
args=parser.parse_args()
HAP1,MES=args.hap1.expanduser().resolve(),args.mes.expanduser().resolve()
PHASEA,HOLD,OLD_STATS=(args.phasea_summary.expanduser().resolve(),
                       args.holdout_groups.expanduser().resolve(),
                       args.covariance_summary.expanduser().resolve())
OUT=args.output_dir.expanduser().resolve()
OUT.mkdir(parents=True,exist_ok=True)
for required in (HAP1,MES,PHASEA,HOLD,OLD_STATS):
    if not required.is_file(): raise FileNotFoundError(f'Required BRCA2 input is missing: {required}')

def load_mapped(path):
    d=pd.read_csv(path,dtype=str,keep_default_na=False)
    d['score']=pd.to_numeric(d['score'],errors='coerce')
    d['protein_pos']=pd.to_numeric(d['protein_pos'],errors='coerce')
    d=d.dropna(subset=['score','protein_pos'])
    d=d[d['class'].eq('missense')].copy()
    d['z']=(d.score-d.score.mean())/d.score.std(ddof=0)
    return d

def map_groups(df):
    g=defaultdict(dict)
    for _,r in df.iterrows():
        key=(int(r.protein_pos),str(r.ref_aa),str(r.alt_aa))
        g[key][str(r.hgvs_nt)]={'score':float(r.score),'z':float(r.z),'cds_pos':str(r.get('cds_pos',''))}
    return g

def p_diff_any(v):
    n=len(v)
    if n<2: return 0.0
    c=Counter(float(x) for x in v)
    equal=sum(k*(k-1)//2 for k in c.values())
    return (n*(n-1)//2-equal)/(n*(n-1)/2)

def selected(x,y,keys,mode):
    n=len(x)
    if mode=='hap1_extreme':
        if n==2: i,j=0,1; da=x[i]-x[j]; db=y[i]-y[j]
        else:
            order=np.argsort(x); i,j=int(order[0]),int(order[-1]); da=x[j]-x[i]; db=y[j]-y[i]
    elif mode=='mes_extreme':
        if n==2: i,j=0,1; da=x[i]-x[j]; db=y[i]-y[j]
        else:
            order=np.argsort(y); i,j=int(order[0]),int(order[-1]); da=x[j]-x[i]; db=y[j]-y[i]
    elif mode=='fixed_first2':
        i,j=0,1; da=x[i]-x[j]; db=y[i]-y[j]
    else: raise ValueError(mode)
    return i,j,float(da),float(db),int(float(da)*float(db)>0)

def endpoint_tie_summary(axis,x,y):
    if len(axis)==2:
        k=int(float(x[0]-x[1])*float(y[0]-y[1])>0)
        return (False,False,k,k)
    if len(axis)<2: return (False,False,0,0)
    lo=float(np.min(axis)); hi=float(np.max(axis))
    low=np.flatnonzero(np.asarray(axis)==lo); high=np.flatnonzero(np.asarray(axis)==hi)
    signs=[]; pairn=0
    for i in low:
        for j in high:
            if int(i)==int(j): continue
            pairn+=1
            if axis is x: da=float(x[j]-x[i]); db=float(y[j]-y[i])
            else: da=float(x[j]-x[i]); db=float(y[j]-y[i])
            signs.append(int(da*db>0))
    return (len(low)>1 or len(high)>1, len(set(signs))>1 if signs else False,
            min(signs) if signs else 0, max(signs) if signs else 0)

def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''): h.update(b)
    return h.hexdigest()

A=load_mapped(HAP1)
B=load_mapped(MES)
ga,gb=map_groups(A),map_groups(B)
rows=[]
for group in sorted(set(ga)&set(gb)):
    keys=sorted(set(ga[group])&set(gb[group]))
    if len(keys)<2: continue
    x=np.array([ga[group][k]['z'] for k in keys],float)
    y=np.array([gb[group][k]['z'] for k in keys],float)
    rawx=np.array([ga[group][k]['score'] for k in keys],float)
    rawy=np.array([gb[group][k]['score'] for k in keys],float)
    one={}
    for mode,label in [('hap1_extreme','primary'),('mes_extreme','reverse_source'),('fixed_first2','fixed_pair')]:
        i,j,da,db,k=selected(x,y,keys,mode)
        axis=x if mode!='mes_extreme' else y
        endpoint_tie,sign_ambig,tie_min,tie_max=endpoint_tie_summary(axis,x,y)
        if mode=='mes_extreme': p_nonzero=p_diff_any(x); selected_nonzero=(db!=0)
        else: p_nonzero=p_diff_any(y); selected_nonzero=(da!=0)
        exp=0.5*p_nonzero if selected_nonzero else 0.0
        one[label]=(i,j,da,db,k,exp,endpoint_tie,sign_ambig,tie_min,tie_max)
    i,j,da,db,k,exp,et,amb,tmin,tmax=one['primary']
    ri,rj,rda,rdb,rk,rexp,ret,ramb,rtmin,rtmax=one['reverse_source']
    fi,fj,fda,fdb,fk,fexp,fet,famb,ftmin,ftmax=one['fixed_pair']
    cx=Counter(rawx.tolist()); cy=Counter(rawy.tolist())
    hx_pairs=sum(v*(v-1)//2 for v in cx.values())
    my_pairs=sum(v*(v-1)//2 for v in cy.values())
    rows.append({
        'protein_pos':group[0],'ref_aa':group[1],'alt_aa':group[2],'group_id':f'{group[1]}{group[0]}{group[2]}','n_variants':len(keys),
        'shared_hgvs_keys':'|'.join(keys),
        'Hap1_any_exact_score_tie':len(set(rawx.tolist()))<len(keys),'mES_any_exact_score_tie':len(set(rawy.tolist()))<len(keys),
        'Hap1_tied_variant_pairs':hx_pairs,'mES_tied_variant_pairs':my_pairs,
        'primary_Hap1_min_key':keys[i],'primary_Hap1_max_key':keys[j],'primary_d_Hap1_z':da,'primary_d_mES_z':db,
        'primary_Hap1_axis_zero':da==0,'primary_mES_axis_zero':db==0,'primary_concordant':k,'primary_tie_adjusted_chance':exp,
        'primary_score_extreme_endpoint_tie':et,'primary_tie_break_sign_ambiguous':amb,'primary_tie_break_min_concordant':tmin,'primary_tie_break_max_concordant':tmax,
        'reverse_mES_min_key':keys[ri],'reverse_mES_max_key':keys[rj],'reverse_d_Hap1_z':rda,'reverse_d_mES_z':rdb,
        'reverse_Hap1_axis_zero':rda==0,'reverse_mES_axis_zero':rdb==0,'reverse_concordant':rk,'reverse_tie_adjusted_chance':rexp,
        'reverse_score_extreme_endpoint_tie':ret,'reverse_tie_break_sign_ambiguous':ramb,'reverse_tie_break_min_concordant':rtmin,'reverse_tie_break_max_concordant':rtmax,
        'fixed_pair_key1':keys[fi],'fixed_pair_key2':keys[fj],'fixed_d_Hap1_z':fda,'fixed_d_mES_z':fdb,
        'fixed_Hap1_axis_zero':fda==0,'fixed_mES_axis_zero':fdb==0,'fixed_concordant':fk,'fixed_tie_adjusted_chance':fexp,
        'Hap1_all_pairwise_nonzero_probability':p_diff_any(rawx),'mES_all_pairwise_nonzero_probability':p_diff_any(rawy)
    })

df=pd.DataFrame(rows)
assert len(df)==388, f'expected 388 groups, got {len(df)}'
assert int(df.n_variants.sum())==823, f'expected 823 shared HGVS keys, got {df.n_variants.sum()}'
assert int(df.protein_pos.nunique())==367
phase=json.loads(PHASEA.read_text())
hold=pd.read_csv(HOLD)
old_phi=float(phase['holdout']['phi_A'])
old_ci=list(phase['holdout']['phi_A_ci'])
old_r=float(phase['holdout']['replication']['all']['R'])
assert abs(float(df.primary_concordant.mean())-old_r)<1e-12, (df.primary_concordant.mean(),old_r)
old_stats=json.loads(OLD_STATS.read_text())
old_cluster_ci=old_stats['position_cluster_bootstrap_5000_seed_20261122']['direction_R_ci']

def cluster_ci(sub, col, seed):
    if sub.empty: return {'n_positions':0,'ci95_percentile':[None,None]}
    by=sub.groupby('protein_pos',sort=False)[col].agg(['count','sum'])
    P=len(by); rng=np.random.default_rng(seed)
    W=rng.multinomial(P,np.full(P,1/P),size=5000)
    den=W@by['count'].to_numpy(float); num=W@by['sum'].to_numpy(float)
    vals=num/den
    return {'n_positions':P,'ci95_percentile':[float(np.quantile(vals,.025)),float(np.quantile(vals,.975))]}

summary={'source_files':{
  'mapped_1225':{'path':str(HAP1),'sha256':sha(HAP1)},
  'mapped_1242':{'path':str(MES),'sha256':sha(MES)},
  'holdout_paired_groups':{'path':str(HOLD),'sha256':sha(HOLD)}},
 'group_universe':{'groups':len(df),'shared_HGVS_variants':int(df.n_variants.sum()),'protein_positions':int(df.protein_pos.nunique()),
   'group_size_counts':{str(k):int(v) for k,v in df.n_variants.value_counts().sort_index().items()}},
 'preserved_primary':{'Phi_A':old_phi,'Phi_A_original_group_bootstrap_CI95':old_ci,'existing_position_cluster_R_CI95':old_cluster_ci,
   'primary_R_reconstructed':old_r,'primary_R_reconstruction_matches_phaseA':True},
 'tie_definitions':'Exact equality of the exported score values (equality is invariant under per-assay z-scoring); tied variant-pair count is the number of within-group unordered pairs having equal exported score.',
 'chance_benchmark_definition':'For each method, conditional on one assay and the selected pair under that method, uniformly reassign the other assay’s observed scores over the group’s variant keys. The exact expected concordance is 0.5 times the probability that the selected endpoints receive unequal scores, with all observed tie multiplicities preserved. This is a conditional exchangeability benchmark, not a p-value or proof that variants are exchangeable.',
 'position_cluster_bootstrap':'For new sensitivity estimators, resample observed protein-position clusters with replacement (same number of nonempty positions for the estimator), retain every group within each selected position, 5,000 draws, seed 20261123 + method index. Percentile 95% intervals; descriptive dependence-robust intervals.',
 'estimators':{}}
for idx,(name,col,chance,axis_a,axis_b) in enumerate([
 ('primary_Hap1_extremes','primary_concordant','primary_tie_adjusted_chance','primary_Hap1_axis_zero','primary_mES_axis_zero'),
 ('reverse_mES_extremes','reverse_concordant','reverse_tie_adjusted_chance','reverse_Hap1_axis_zero','reverse_mES_axis_zero'),
 ('fixed_first_two_HGVS_keys','fixed_concordant','fixed_tie_adjusted_chance','fixed_Hap1_axis_zero','fixed_mES_axis_zero')]):
    summary['estimators'][name]={}
    for label,sub in [('all',df),('n2',df[df.n_variants==2]),('n3plus',df[df.n_variants>=3])]:
        n=len(sub); k=int(sub[col].sum()) if n else 0
        summary['estimators'][name][label]={
          'groups':n,'protein_positions':int(sub.protein_pos.nunique()),'concordant':k,'R':float(k/n) if n else None,
          'Hap1_selected_difference_zero':int(sub[axis_a].sum()),'mES_selected_difference_zero':int(sub[axis_b].sum()),
          'either_selected_axis_zero':int((sub[axis_a]|sub[axis_b]).sum()),'both_selected_axes_zero':int((sub[axis_a]&sub[axis_b]).sum()),
          'any_Hap1_score_ties_within_group':int(sub.Hap1_any_exact_score_tie.sum()),'any_mES_score_ties_within_group':int(sub.mES_any_exact_score_tie.sum()),
          'groups_with_score_extreme_endpoint_ties':int(sub[f'{name.split("_")[0]}_score_extreme_endpoint_tie'].sum()) if name!='fixed_first_two_HGVS_keys' else 0,
          'groups_with_tie_break_sign_ambiguity':int(sub[f'{name.split("_")[0]}_tie_break_sign_ambiguous'].sum()) if name!='fixed_first_two_HGVS_keys' else 0,
          'tie_adjusted_expected_R':float(sub[chance].mean()) if n else None,
          'R_minus_tie_adjusted_expected_R':float(sub[col].mean()-sub[chance].mean()) if n else None,
          'position_cluster_bootstrap_95':cluster_ci(sub,col,20261123+idx*10+(0 if label=='all' else 1 if label=='n2' else 2))}
# Ensure known legacy interpretation remains exact for primary.
summary['selected_pair_notes']={
 'primary_Hap1_extremes':'Matches run_phaseA.py: n=2 uses sorted-HGVS variants 0 minus 1; n>=3 uses np.argsort(Hap1 z) extrema and reads mES values for the same two variants.',
 'reverse_mES_extremes':'Sensitivity reverses the data-dependent selection source for n>=3; n=2 pair ordering remains the same.',
 'fixed_first_two_HGVS_keys':'Uses the first two shared exact HGVS keys sorted lexicographically in every group, independent of scores; this matches the fixed key order already used by the two-variant estimator.'}
summary['diagnostics']={
 'primary_extreme_tie_break_ambiguity_groups_all':int(df.primary_tie_break_sign_ambiguous.sum()),
 'primary_extreme_tie_break_ambiguity_groups_n3plus':int(df.loc[df.n_variants>=3,'primary_tie_break_sign_ambiguous'].sum()),
 'primary_tie_break_R_min_max_all':[float(df.primary_tie_break_min_concordant.mean()),float(df.primary_tie_break_max_concordant.mean())],
 'reverse_tie_break_ambiguity_groups_all':int(df.reverse_tie_break_sign_ambiguous.sum()),
 'reverse_tie_break_ambiguity_groups_n3plus':int(df.loc[df.n_variants>=3,'reverse_tie_break_sign_ambiguous'].sum()),
 'reverse_tie_break_R_min_max_all':[float(df.reverse_tie_break_min_concordant.mean()),float(df.reverse_tie_break_max_concordant.mean())]}
assert summary['estimators']['primary_Hap1_extremes']['all']['R']==old_r
outcsv=OUT/'brca2_direction_group_sensitivities_20261123.csv'
df.to_csv(outcsv,index=False)
(OUT/'brca2_direction_summary_20261123.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False,sort_keys=True)+'\n',encoding='utf-8')
(OUT/'brca2_direction_methods_20261123.md').write_text(
'''# BRCA2 paired direction and tie audit\n\nInputs are the mapped BRCA2 HAP1 (`urn:mavedb:00001225-a-1`) and mES (`urn:mavedb:00001242-a-1`) tables. For each assay, scores were z-standardized across all missense rows with population SD, matching `run_phaseA.py`; same-position, same-reference-amino-acid, same-alternate-amino-acid groups were intersected by exact `hgvs_nt`.\n\nThe primary direction estimator is reproduced exactly: two-member groups compare the lexicographically sorted HGVS keys; groups with three or more members select HAP1-score extremes via `np.argsort`, then compare mES scores for those same variants. Sensitivities use mES-score extremes (reverse-source selection) or the first two shared HGVS keys independent of scores. Ties are exact equalities of the exported MaveDB `score` values; z-scoring preserves equalities.\n\nThe tie-adjusted chance rate is a conditional expectation under uniformly reassignment of the non-selection assay's observed scores across variant keys within each group, preserving score ties. It is reported as an expectation under an exchangeability benchmark, not a p-value. Position-cluster percentile intervals resample protein positions and retain all groups per sampled position.\n\nThe original Phase-A Phi and its interval are copied into the summary as preserved values. This audit does not revise Phi, infer clinical independence, or establish that variants are exchangeable within groups. See the JSON and group-level CSV for counts, selected endpoints, tie structures, and interval values.\n''',encoding='utf-8')
print(json.dumps({'groups':len(df),'group_sizes':summary['group_universe']['group_size_counts'],'primary_R':old_r,'summary':str(OUT/'brca2_direction_summary_20261123.json'),'groups_csv':str(outcsv)},sort_keys=True),flush=True)

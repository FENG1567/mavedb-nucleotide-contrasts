#!/usr/bin/env python3
"""Independent QC for VHL exact-DNA counts and SGE target-region mixing."""
from pathlib import Path
from hashlib import sha256
import argparse
import io, json, zipfile
import numpy as np
import pandas as pd
import openpyxl

REPO_ROOT = Path(__file__).resolve().parents[2]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--input-zip', type=Path, required=True,
                    help='same article supplementary ZIP used by the main analysis')
parser.add_argument('--output-dir', type=Path, default=REPO_ROOT / 'results' / 'secondary' / 'vhl',
                    help='directory containing main-analysis outputs; QC writes to its qc/ subdirectory')
parser.add_argument('--main-script', type=Path, default=Path(__file__).with_name('analyze_vhl_crossrep.py'),
                    help='main-analysis script whose actual SHA-256 is recorded')
args = parser.parse_args()
ZIP = args.input_zip.resolve()
OUT = args.output_dir.resolve()
MAIN_SCRIPT = args.main_script.resolve()
QC = OUT / 'qc'
QC.mkdir(parents=True, exist_ok=True)
MEMBER = '41588_2024_1800_MOESM4_ESM.xlsx'
BOOT = 2000
SEED = 20261123

# Load the exact retained variant set and the already-downloaded ST1 region metadata.
eligible=pd.read_csv(OUT/'variant_level.csv')
with zipfile.ZipFile(ZIP) as zf:
    member_bytes=zf.read(MEMBER)
    wb=openpyxl.load_workbook(io.BytesIO(member_bytes),read_only=True,data_only=True)
    ws=wb['ST1']
    headers=[x for x in next(ws.iter_rows(min_row=3,max_row=3,values_only=True))]
    raw=pd.DataFrame(list(ws.iter_rows(min_row=4,values_only=True)),columns=headers).dropna(how='all')
    wb.close()

required=['hg38_pos','ref','alt','cHGVS','pHGVS','protPos','sge_region',
          'tHDR_pre','tHDR_post','tHDR_rna','rLD2_tHDR_pre','rLD2_tHDR_post','rLD2_tHDR_rna']
missing=[x for x in required if x not in raw.columns]
if missing: raise RuntimeError(f'Missing ST1 columns: {missing}')
raw['dna_key']=raw.apply(lambda r:f"hg38:{int(r.hg38_pos)}:{str(r.ref).strip()}>{str(r.alt).strip()}",axis=1)
if raw['dna_key'].duplicated().any(): raise AssertionError('Unexpected duplicate genomic DNA keys in ST1; review before grouping')
st1=raw[raw['dna_key'].isin(eligible['dna_key'])].copy()
if len(st1)!=len(eligible): raise AssertionError(f"Join mismatch: ST1 eligible rows={len(st1)}, saved main eligible variants={len(eligible)}")
st1=st1.merge(eligible[['dna_key','F1','R1','F2','R2']],on='dna_key',how='left',validate='one_to_one')
st1['protPos']=pd.to_numeric(st1['protPos'],errors='coerce')
for c in ['tHDR_pre','tHDR_post','tHDR_rna','rLD2_tHDR_pre','rLD2_tHDR_post','rLD2_tHDR_rna']:
    st1[c]=pd.to_numeric(st1[c],errors='coerce')

count_map={'DNA6_1':'tHDR_pre','DNA20_1':'tHDR_post','RNA6_1':'tHDR_rna',
           'DNA6_2':'rLD2_tHDR_pre','DNA20_2':'rLD2_tHDR_post','RNA6_2':'rLD2_tHDR_rna'}
issues=[]
for name,col in count_map.items():
    x=st1[col].to_numpy(dtype=float)
    finite=np.isfinite(x)
    nonnegative=finite & (x>=0)
    integer=finite & (x==np.floor(x))
    if not finite.all() or not nonnegative.all() or not integer.all():
        issues.append({'analysis_name':name,'source_column':col,'n':len(x),
                       'nonfinite':int((~finite).sum()),'negative':int((finite & (x<0)).sum()),
                       'noninteger':int((finite & (x!=np.floor(x))).sum())})
if issues: raise AssertionError(f'Count validation failed: {issues}')
prot_counts=eligible.groupby('pHGVS')['protPos'].nunique(dropna=False)
if not prot_counts.eq(1).all(): raise AssertionError('At least one exact pHGVS group maps to multiple protPos values')

# Record the actual main-analysis implementation used for this run.
script_sha = sha256(MAIN_SCRIPT.read_bytes()).hexdigest()

eligible_groups=eligible[['pHGVS','protPos','cHGVS','dna_key','F1','R1','F2','R2']].merge(
    st1[['dna_key','sge_region']],on='dna_key',how='left',validate='one_to_one')
if eligible_groups['sge_region'].isna().any(): raise AssertionError('Missing target region for retained DNA variants')
reg_audit=[]
for p,g in eligible_groups.groupby('pHGVS',sort=True):
    regs=sorted(g['sge_region'].astype(str).unique())
    reg_audit.append({'pHGVS':p,'protPos':int(g['protPos'].iloc[0]),'dna_variants':int(g['dna_key'].nunique()),
                      'n_sge_regions':len(regs),'sge_regions':'|'.join(regs),'region_constant':len(regs)==1})
pd.DataFrame(reg_audit).to_csv(QC/'region_by_group.csv',index=False)
mixed=set(pd.DataFrame(reg_audit).loc[lambda x:~x['region_constant'],'pHGVS'])

def group_metrics(groups, keys, xcol, ycol):
    covs=[]; xc=[]; yc=[]; conc=0; sign_n=0; ties=0
    for k in keys:
        g=groups[k].sort_values('cHGVS',kind='stable')
        x=g[xcol].to_numpy(float); y=g[ycol].to_numpy(float)
        if len(x)<2: continue
        covs.append(float(np.cov(x,y,ddof=1)[0,1]))
        xc.extend((x-x.mean()).tolist()); yc.extend((y-y.mean()).tolist())
        dx=x[0]-x[1]; dy=y[0]-y[1]
        if dx==0 or dy==0: ties+=1
        else: sign_n+=1; conc+=int(np.sign(dx)==np.sign(dy))
    r=float(np.corrcoef(xc,yc)[0,1]) if len(xc)>1 and np.std(xc)>0 and np.std(yc)>0 else None
    return {'mean_group_covariance':float(np.mean(covs)) if covs else None,'variant_weighted_groupcentered_r':r,
            'first2_sign_concordance':conc/sign_n if sign_n else None,'sign_groups_n':sign_n,'sign_ties':ties}

sens=[]; sensitivity_groups=0; sensitivity_variants=0
if mixed:
    # Regroup by exact pHGVS and source-reported target region; retain only groups with >=2 DNA variants.
    eligible_groups['sens_group']=eligible_groups['pHGVS'].astype(str)+'||'+eligible_groups['sge_region'].astype(str)
    counts=eligible_groups.groupby('sens_group')['dna_key'].nunique()
    retained=set(counts[counts>=2].index)
    data=eligible_groups[eligible_groups['sens_group'].isin(retained)].copy()
    groups={k:g.reset_index(drop=True) for k,g in data.groupby('sens_group',sort=True)}
    sensitivity_groups=len(groups); sensitivity_variants=len(data)
    pos_groups={}
    for k,g in groups.items(): pos_groups.setdefault(int(g['protPos'].iloc[0]),[]).append(k)
    positions=np.array(sorted(pos_groups),dtype=int)
    pairs=[('within_rep1','F1','R1'),('within_rep2','F2','R2'),('cross_F1_R2','F1','R2'),('cross_F2_R1','F2','R1')]
    rng=np.random.default_rng(SEED)
    bstore={p:{m:[] for m in ('mean_group_covariance','variant_weighted_groupcentered_r','first2_sign_concordance')} for p,_,_ in pairs}
    for name,x,y in pairs:
        sens.append({'pairing':name,'function_measure':x,'rna_measure':y,
                     'eligible_same_region_groups':sensitivity_groups,'eligible_DNA_variants':sensitivity_variants,
                     **group_metrics(groups,list(groups),x,y)})
    if len(positions) and groups:
        for _ in range(BOOT):
            sampled=rng.choice(positions,size=len(positions),replace=True)
            keys=[]
            for pos in sampled: keys.extend(pos_groups[int(pos)])
            for name,x,y in pairs:
                vals=group_metrics(groups,keys,x,y)
                for metric in bstore[name]:
                    if vals[metric] is not None and np.isfinite(vals[metric]): bstore[name][metric].append(vals[metric])
        for row in sens:
            b=bstore[row['pairing']]
            for metric in b:
                vals=np.asarray(b[metric],dtype=float)
                row[metric+'_ci95_low']=float(np.quantile(vals,.025)) if len(vals) else None
                row[metric+'_ci95_high']=float(np.quantile(vals,.975)) if len(vals) else None
                row[metric+'_bootstrap_n']=int(len(vals))
    pd.DataFrame(sens).to_csv(QC/'same_region_sensitivity.csv',index=False)

qc={'main_analysis_script_sha256':script_sha,
    'main_analysis_script_path':str(MAIN_SCRIPT),'main_summary_sha256':sha256((OUT/'summary.json').read_bytes()).hexdigest(),
    'main_crossrep_csv_sha256':sha256((OUT/'crossrep_summary.csv').read_bytes()).hexdigest(),
    'source_workbook_member_sha256':sha256(member_bytes).hexdigest(),'retained_variants':int(len(st1)),
    'eligible_groups':int(eligible['pHGVS'].nunique()),'count_columns':{},'count_validation_passed':not issues,
    'pHGVS_group_protPos_unique_assertion_passed':True,'n_sge_region_column_values':int(st1['sge_region'].nunique()),
    'region_mixed_groups_n':len(mixed),'region_constant_groups_n':int(len(reg_audit)-len(mixed)),
    'region_groups':reg_audit,'same_region_sensitivity_run':bool(mixed),'same_region_sensitivity_groups':sensitivity_groups,
    'same_region_sensitivity_variants':sensitivity_variants,'same_region_sensitivity':sens,
    'sensitivity_groupcentered_r_weighting':'variant-weighted: each retained DNA variant contributes once after centering within pHGVS+region group; larger groups contribute more',
    'covariance_weighting':'equal weight per pHGVS+region group','bootstrap':{'cluster':'protPos','iterations':BOOT if mixed else 0,'seed':SEED}}
for name,col in count_map.items():
    x=st1[col].to_numpy(dtype=float)
    qc['count_columns'][name]={'source_column':col,'n':len(x),'finite':int(np.isfinite(x).sum()),
                               'nonnegative':int((x>=0).sum()),'integer':int((x==np.floor(x)).sum()),
                               'min':float(x.min()),'max':float(x.max())}
(QC/'qc_summary.json').write_text(json.dumps(qc,indent=2,ensure_ascii=False),encoding='utf-8')
print(json.dumps({k:v for k,v in qc.items() if k!='region_groups'},ensure_ascii=False))

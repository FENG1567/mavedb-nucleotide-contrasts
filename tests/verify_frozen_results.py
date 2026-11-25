#!/usr/bin/env python3
"""Recalculate manuscript-facing aggregates; no raw-data or bootstrap claim."""
from pathlib import Path
import csv
import gzip
import json
import math
import sys
from collections import Counter, defaultdict

ROOT = Path(__file__).resolve().parents[1]
D = ROOT / 'data' / 'derived'
V4 = D / 'revision_20261123_review'
checks = []

def check(name, passed, observed=None):
    checks.append({'name': name, 'passed': bool(passed), 'observed': observed})
    if not passed:
        raise AssertionError(f'{name}: {observed}')

def rows(p):
    opener = gzip.open if str(p).endswith('.gz') else open
    with opener(p, 'rt', encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))

def truth(x):
    return str(x).lower() in ('true', '1')

def close(a,b):
    return math.isclose(float(a), float(b), rel_tol=1e-10, abs_tol=1e-11)

def main():
    groups = rows(V4 / 'collection/group_level_empirical_testability.csv')
    per_set = rows(V4 / 'collection/per_score_set_reconciliation.csv')
    check('collection historical groups', len(groups)==13041, len(groups))
    check('collection frozen sets', len(per_set)==92, len(per_set))
    for field, state, expected in [
        ('accession_reference_group_state','ge2',3836),
        ('versioned_request_group_state','ge2',2679),
        ('target_sequence_group_state','conditional_ge2',1013)]:
        n=sum(r[field]==state for r in groups)
        check(field, n==expected, n)
    check('historical keys from groups',sum(int(r['legacy_unique_keys']) for r in groups)==31723)

    bg = rows(V4 / 'pairs/brca2_direction_group_sensitivities_20261123.csv')
    check('BRCA2 group count',len(bg)==388,len(bg))
    check('BRCA2 shared keys',sum(int(r['n_variants']) for r in bg)==823)
    check('BRCA2 positions',len(set(r['protein_pos'] for r in bg))==367)
    check('BRCA2 size distribution',Counter(int(r['n_variants']) for r in bg)==Counter({2:341,3:47}))
    check('BRCA2 concordance',sum(int(r['primary_concordant']) for r in bg)==193)
    for size,n in [(2,171),(3,22)]:
        check(f'BRCA2 size {size} concordance',sum(int(r['primary_concordant']) for r in bg if int(r['n_variants'])==size)==n)
    h=rows(D / 'work/holdout_paired_groups.csv')
    phi=1-sum(float(r['cov']) for r in h)/sum(float(r['mean_var']) for r in h)
    check('BRCA2 Phi recomputed from 388 group moments',close(phi,.8750878921486126),phi)
    check('BRCA2 Phi differs from class-control estimate',not close(phi,.8758893944894707))

    m=rows(V4 / 'msh2/output/primary_metrics.csv')
    e=rows(V4 / 'msh2/output/primary_group_errors.csv.gz')
    expected={1:(8547,24619,887),2:(8416,23997,900),3:(8474,24156,892)}
    for held,cohort in expected.items():
        metric=next(r for r in m if int(r['held_rep'])==held and r['endpoint']=='6TG_mock')
        sub=[r for r in e if int(r['held_rep'])==held and r['endpoint']=='6TG_mock']
        observed=(len(sub),sum(int(r['n_variants']) for r in sub),len(set(r['pos'] for r in sub)))
        check(f'MSH2 held {held} cohort',observed==cohort,observed)
        a=sum(float(r['mse_aa']) for r in sub)/len(sub)
        f=sum(float(r['mse_full']) for r in sub)/len(sub)
        gain=1-math.sqrt(f/a)
        check(f'MSH2 held {held} equal-group MSE',close(a,metric['mse_aa']) and close(f,metric['mse_full']))
        check(f'MSH2 held {held} RMSE gain',close(gain,metric['relative_rmse_gain']),gain)
    sp=rows(D/'context/splicing_fixed_pair_manifest_20261123_v2.csv')
    pair_members=Counter(r['group_id'] for r in sp if truth(r['selected_for_fixed_pair']))
    check('splicing fixed 22 pairs',len(pair_members)==22 and set(pair_members.values())=={2},dict(pair_members))
    report={'scope':'Frozen-input aggregate recalculation; does not replace raw-source reconstruction or bootstrap reruns.',
            'passed':True,'checks':checks}
    out=ROOT/'results/verification'
    out.mkdir(parents=True,exist_ok=True)
    (out/'frozen_aggregate_checks.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps(report,indent=2))

if __name__=='__main__':
    main()

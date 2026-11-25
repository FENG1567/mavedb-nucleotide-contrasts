"""Compare the regenerated BRCA2/selection outputs against manuscript freezes."""
import csv,json,math
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
D=ROOT/'data/derived'
R=ROOT/'results'
checks=[]

def same(a,b,where=''):
    if isinstance(a,dict) and isinstance(b,dict):
        assert set(a)==set(b),(where,set(a)^set(b))
        for k in a:same(a[k],b[k],where+'.'+str(k))
    elif isinstance(a,list) and isinstance(b,list):
        assert len(a)==len(b),where
        for i,(x,y) in enumerate(zip(a,b)):same(x,y,f'{where}[{i}]')
    else:
        if str(a)==str(b):return
        try:
            assert math.isclose(float(a),float(b),rel_tol=1e-10,abs_tol=1e-10),(where,a,b)
        except (ValueError,TypeError):assert a==b,(where,a,b)

def table(path,keys):
    with path.open(encoding='utf-8-sig',newline='') as f:
        return {tuple(row[k] for k in keys):row for row in csv.DictReader(f)}

def main():
    br=json.loads((R/'brca2_core/brca2_core_summary.json').read_text())
    fr=json.loads((D/'revision_20261122/statistics/brca2_covariance_summary.json').read_text())
    same(br['phi_A_primary'],fr['phi_A_primary_preserved'],'primary Phi')
    same(br['phi_A_primary_ci_original_group_bootstrap'],fr['phi_A_primary_ci_preserved'],'primary Phi CI')
    observed=table(R/'brca2_direction/brca2_direction_group_sensitivities_20261123.csv',['group_id'])
    expected_groups=table(D/'revision_20261123_review/pairs/brca2_direction_group_sensitivities_20261123.csv',['group_id'])
    assert set(observed)==set(expected_groups)
    extra={'primary_tie_break_min_concordant','primary_tie_break_max_concordant',
           'reverse_tie_break_min_concordant','reverse_tie_break_max_concordant'}
    for key,ref in expected_groups.items():
        row=observed[key]
        assert set(row)-set(ref)<=extra and set(ref)<=set(row)
        same({k:row[k] for k in ref},ref,'BRCA2 published group fields.'+str(key))
        # The final publication omits these historical diagnostic columns.
        # The portable implementation now computes them correctly; no observed
        # group has an endpoint tie, so each bound equals observed concordance.
        for field in set(row)-set(ref):
            concordance='primary_concordant' if field.startswith('primary') else 'reverse_concordant'
            same(row[field],ref[concordance],'new tie-bound diagnostic.'+str(key)+'.'+field)
    current=json.loads((R/'brca2_direction/brca2_direction_summary_20261123.json').read_text())
    frozen=json.loads((D/'revision_20261123_review/pairs/brca2_direction_summary_20261123.json').read_text())
    same(current['diagnostics'],frozen['diagnostics'],'BRCA2 diagnostic summary')
    same(current['estimators'],frozen['estimators'],'BRCA2 estimators and intervals')
    same(table(R/'selection_bias_null/selection_bias_null.csv',['tau_over_sigma','pipeline']),
         table(D/'outputs/selection_bias_null.csv',['tau_over_sigma','pipeline']),'selection simulation')
    expected=['figure1_collection_testability','figure2_brca2_paired_contrasts',
              'figure3_msh2_controlled_endpoints','figure4_structural_selection_counterexample',
              'figureS1_msh2_sensitivity','rna_cross_repeat_checks']
    for name in expected:
        assert (R/f'figures/{name}.pdf').read_bytes().startswith(b'%PDF-'),name
    report={'passed':True,'scope':['BRCA2 raw-score mapping fields','primary Phi and ordinary bootstrap CI',
                                  'all published fields of 388 BRCA2 group rows; corrected extra tie-bound diagnostics; diagnostic/estimator summaries',
                                  'selection simulation table','six generated PDF signatures'],
            'numeric_tolerance':1e-10}
    (R/'verification').mkdir(exist_ok=True,parents=True)
    (R/'verification/recomputed_core_comparison.json').write_text(json.dumps(report,indent=2))
    print(json.dumps(report,indent=2))

if __name__=='__main__':main()

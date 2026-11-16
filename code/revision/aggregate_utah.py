"""Summarize the new cohort with one-bone-per-donor confidence intervals."""
from __future__ import annotations
from collections import Counter
import csv
import hashlib
import json
import os
from pathlib import Path
import numpy as np
from scipy.stats import beta

R = Path(__file__).resolve().parents[2]
E = R / 'external_extension_2026_11_12'


def load(p): return json.loads(p.read_text())
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def dump(p, d): p.write_text(json.dumps(d, indent=2, ensure_ascii=False) + '\n')


def ci(k, n):
    if not n: return None
    return [0.0 if k == 0 else float(beta.ppf(.025, k, n-k+1)),
            1.0 if k == n else float(beta.ppf(.975, k+1, n-k))]


def describe(rows, independent):
    accepted = [r for r in rows if r['accepted']]
    evaluable = [r for r in rows if r['source_evaluable']]
    a = {'attempted': len(rows), 'source_evaluable': len(evaluable), 'accepted': len(accepted),
         'statuses': dict(Counter(r['outcome'] for r in rows)),
         'donors': len({r['donor_id'] for r in rows}),
         'minimum_continuous_clearance_mm': min((r['continuous_minimum_mm'] for r in accepted), default=None),
         'contracted_outputs': sum(r['radial_scale'] < 1 for r in accepted),
         'confidence_interval_unit': 'one selected humerus per source-defined donor' if independent else 'not computed for correlated bilateral bones'}
    if independent:
        assert len(rows) == a['donors']
        a['acceptance_per_attempted_exact_95_ci'] = ci(len(accepted), len(rows))
        a['acceptance_per_evaluable_exact_95_ci'] = ci(len(accepted), len(evaluable))
    if accepted:
        a['continuous_clearance_quartiles_mm'] = np.quantile([r['continuous_minimum_mm'] for r in accepted], [.25,.5,.75]).tolist()
        a['wall_seconds_quartiles'] = np.quantile([r['wall_seconds'] for r in accepted], [.25,.5,.75]).tolist()
    return a


def main():
    inv = load(E / 'inventory.json')
    assert sha(E / 'inventory.json') == load(E / 'selection_freeze.json')['inventory_sha256']
    prep = {x['case_id']: x for x in load(E / 'preparation_status.json')}
    runs = {x['case_id']: x for x in load(E / 'execution_status.json')}
    frozen = load(E / 'freeze.json')
    assert len(prep) == len(runs) == len(inv['cases']) == 125
    rows = []
    for c in inv['cases']:
        case = c['case_id']; p = prep[case]; run = runs[case]
        row = {k: c[k] for k in ['case_id','donor_id','analysis','sex','side','age_years','cadaver','pathology','notes','source_sha256','archive_member','source_faces']}
        row.update({'preparation_status': p['status'], 'source_evaluable': p['status'] == 'PREPARED',
                    'accepted': False, 'outcome': p['status'], 'continuous_minimum_mm': None,
                    'sampled_minimum_mm': None, 'radial_scale': None, 'cavity_volume_mm3': None,
                    'wall_seconds': run.get('wall_seconds'), 'reason': p.get('reason', '')})
        if p['status'] == 'PREPARED':
            d = load(E / 'runs' / case / 'decision.json')
            assert d['source_sha256'] == p['prepared_sha256']
            assert d['config_sha256'] == p['config_sha256']
            assert d['pipeline_sha256'] == frozen['files']['pipeline']['sha256']
            assert d['auditor_sha256'] == frozen['files']['auditor']['sha256']
            row.update({'accepted': bool(d['accepted']) and run['returncode'] == 0,
                        'outcome': d['status'], 'reason': d.get('reason', '')})
            if d['status'] == 'WITHHELD_OR_ERROR':
                reason = d['reason']
                if reason.startswith('Fail-closed source validation:'):
                    row.update(outcome='SOURCE_REJECTED_IN_BLENDER', source_evaluable=False)
                elif any(reason.startswith(x) for x in ['Containment guard', 'Fail-closed output validation:', 'Could not find an interior center', 'No outer vertices near validated spine section']):
                    row['outcome'] = 'METHOD_REJECTED'
                else:
                    row['outcome'] = 'TECHNICAL_FAILURE'
            meta = E / 'runs' / case / 'candidate/regional_run_metadata.json'
            if meta.exists():
                m = load(meta)
                assert m['source_topology_valid'] and m['output_all_valid']
                row.update(sampled_minimum_mm=m['relation_diagnostics']['nearest_unsigned_min_mm'],
                           radial_scale=m['containment_guard']['selected_radial_scale'],
                           cavity_volume_mm3=m['inner']['volume_mm3'])
            audit = E / 'runs' / case / 'continuous_audit.json'
            if audit.exists():
                a = load(audit)
                row['continuous_minimum_mm'] = a['continuous_triangle_surface_minimum_mm']
                if d['status'] == 'WITHHELD_CONTINUOUS_CLEARANCE':
                    row['reason'] = f"Continuous minimum {row['continuous_minimum_mm']:.12f} mm below the 0.150-mm final threshold."
                assert a['required_clearance_mm'] == .15
                assert a['continuous_clearance_pass'] == (row['continuous_minimum_mm'] + 1e-9 >= .15)
                if row['accepted']:
                    assert a['continuous_clearance_pass']
                    for name, key in [('outer.stl','outer_sha256'), ('cavity.stl','inner_sha256')]:
                        assert sha(E / 'runs' / case / 'accepted' / name) == a[key]
            assert (E / 'runs' / case / 'accepted').exists() == row['accepted']
        rows.append(row)
    primary = [r for r in rows if r['analysis'] == 'primary']
    left = [r for r in rows if r['analysis'] == 'secondary_left']
    donor_pairs = sorted({r['donor_id'] for r in left})
    summary = {'source_doi': inv['source_doi'], 'source_version': inv['source_version'],
               'protocol_freeze_utc': frozen['frozen_at_utc'],
               'selection_freeze_utc': load(E / 'selection_freeze.json')['frozen_at_utc'],
               'primary': describe(primary, True), 'secondary_left': describe(left, True),
               'all_bones_descriptive': describe(rows, False),
               'paired_donors': {'attempted': len(donor_pairs),
                                  'both_sides_source_evaluable': sum(all(r['source_evaluable'] for r in rows if r['donor_id']==donor) for donor in donor_pairs),
                                  'both_sides_accepted': sum(all(r['accepted'] for r in rows if r['donor_id']==donor) for donor in donor_pairs)},
               'primary_demographics': {'source_defined_donors':len(primary),
                                       'sex_counts':dict(Counter(r['sex'] for r in primary)),
                                       'age_min_years':min(r['age_years'] for r in primary),
                                       'age_max_years':max(r['age_years'] for r in primary),
                                       'age_median_years':float(np.median([r['age_years'] for r in primary])),
                                       'under_18_years':sum(r['age_years']<18 for r in primary),
                                       'all_cadaveric':all(r['cadaver']=='1' for r in primary)},
               'failures':[{k:r[k] for k in ['case_id','analysis','outcome','reason']} for r in rows if not r['accepted']],
               'environment_setup_failures':125,
               'environment_setup_failures_note':'Missing SciPy before source diagnostics; retained logs and identical retry. Not method failures.',
               'public_release':False, 'anatomical_or_physical_validation':False}
    dump(E / 'results.json', rows); dump(E / 'summary.json', summary)
    with (E / 'results.csv').open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == '__main__': main()

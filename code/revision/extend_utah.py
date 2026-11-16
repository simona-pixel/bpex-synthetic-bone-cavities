"""Acquire metadata, prepare and evaluate the fixed Utah full-humerus cohort.

Run phases sequentially: inventory, prepare, run. Source archive is obtained from
Zenodo separately. The numerical generator, auditor and profile are not modified.
"""
from __future__ import annotations

import argparse
import concurrent.futures
import csv
from datetime import datetime, timezone
import hashlib
import io
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import time
from zipfile import ZipFile

R = Path(__file__).resolve().parents[2]
E = R / 'external_extension_2026_11_12'
H = R.parent / 'benchmark/external_validation_v1/scripts'
PROFILE = R.parent / 'benchmark/external_validation_v1/protocol/EXTERNAL_PARAMETER_PROFILE_v1.json'


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def dump(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False) + '\n')


def verify_freeze():
    f = json.loads((E / 'freeze.json').read_text())
    for item in f['files'].values():
        assert sha(R.parent / item['path']) == item['sha256'], item['path']


def inventory():
    import numpy as np
    import trimesh
    verify_freeze()
    target = E / 'inventory.json'
    if target.exists():
        raise RuntimeError('Inventory already exists; refusing to overwrite it.')
    archive = E / 'acquisition/database_v1.11.zip'
    md5 = hashlib.md5()
    with archive.open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            md5.update(block)
    assert md5.hexdigest() == 'fd87208ea7484b27b6ce0ebfba6711fb'
    cases = []
    with ZipFile(archive) as z:
        for group, expected in [('single', 5), ('paired', 120)]:
            csvs = [n for n in z.namelist() if n.startswith(group + '_humeri_model_data') and n.endswith('.csv')]
            assert len(csvs) == 1
            csvname = csvs[0]
            csvbytes = z.read(csvname)
            (E / 'acquisition' / Path(csvname).name).write_bytes(csvbytes)
            records = list(csv.DictReader(io.StringIO(csvbytes.decode('utf-8-sig'))))
            records = [r for r in records if r['Model ID'].strip()]
            by_model = {r['Model ID']: r for r in records}
            assert len(by_model) == len(records) == expected
            names = sorted(n for n in z.namelist() if n.startswith(group + '_humeri_STLs/') and n.lower().endswith('.stl'))
            assert len(names) == expected
            for name in names:
                model = Path(name).stem
                m = re.fullmatch(r'(single|paired)_humerus_(\d+)_([MF])_(\d+)_([LR])', model)
                assert m and m.group(1) == group, model
                row = by_model[model]
                _, number, sex, age, side = m.groups()
                assert int(row['Age (years)']) == int(age)
                assert int(row['M/F 1/0']) == (sex == 'M')
                assert int(row['R/L 1/0']) == (side == 'R')
                notes = row['Notes'].strip()
                exclusion = None
                if row['Whole Humerus 1/0'] != '1':
                    exclusion = 'metadata_marks_incomplete_humerus'
                if re.search(r'cropp|clipp|incomplete|partial', notes, re.I):
                    exclusion = 'metadata_explicitly_describes_incomplete_geometry'
                raw = z.read(name)
                h = hashlib.sha256(raw).hexdigest()
                mesh = trimesh.load(io.BytesIO(raw), file_type='stl', process=False, maintain_order=True)
                n = len(mesh.faces)
                points = np.unique(np.asarray(mesh.vertices), axis=0)
                assert np.isfinite(points).all()
                geometry_hash = hashlib.sha256(points.astype('<f8').tobytes()).hexdigest()
                case_id = 'utah_' + group + '_' + number + '_humerus_' + side.lower()
                source = E / 'inputs/source' / (case_id + '.stl')
                source.parent.mkdir(parents=True, exist_ok=True)
                source.write_bytes(raw)
                cases.append({'case_id': case_id, 'donor_id': 'utah_' + group + '_' + number,
                              'analysis': 'primary' if group == 'single' or side == 'R' else 'secondary_left',
                              'side': side, 'sex': sex, 'age_years': int(age),
                              'cadaver': row['Cadaver 1/0'], 'pathology': row['Pathology'], 'notes': notes,
                              'archive_member': name, 'metadata_file': csvname, 'model_id': model,
                              'source_relative': str(source.relative_to(E)), 'source_sha256': h,
                              'raw_unique_vertex_set_sha256': geometry_hash,
                              'source_faces': n, 'metadata_exclusion': exclusion, 'metadata': row})
    assert len(cases) == 125
    hashes = {}
    for c in cases:
        key = c['raw_unique_vertex_set_sha256']
        assert key not in hashes, ('Duplicate raw vertex geometry requires review before evaluation', hashes.get(key), c['case_id'])
        hashes[key] = c['case_id']
    donors = sorted({c['donor_id'] for c in cases})
    assert len(donors) == 65
    primary = [c for c in cases if c['analysis'] == 'primary']
    assert len(primary) == len({c['donor_id'] for c in primary}) == 65
    for donor in donors:
        rows = [c for c in cases if c['donor_id'] == donor]
        if 'paired' in donor:
            assert len(rows) == 2 and {c['side'] for c in rows} == {'L', 'R'}
            assert len({(c['age_years'], c['sex']) for c in rows}) == 1
    result = {'created_at_utc': datetime.now(timezone.utc).isoformat(),
              'source_doi': '10.5281/zenodo.19077748', 'source_version': 'v1.11',
              'licence': 'CC-BY-4.0', 'archive_bytes': archive.stat().st_size,
              'archive_md5': md5.hexdigest(), 'archive_sha256': sha(archive),
              'expected_cases': 125, 'source_defined_donors': 65,
              'duplicate_vertex_sets': 0, 'cases': cases}
    dump(target, result)
    with (E / 'cases.csv').open('w', newline='') as f:
        fields = [k for k in cases[0] if k != 'metadata']
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows({k: c[k] for k in fields} for c in cases)
    dump(E / 'selection_freeze.json', {'frozen_at_utc': datetime.now(timezone.utc).isoformat(),
                                      'inventory_sha256': sha(target), 'cases_csv_sha256': sha(E / 'cases.csv'),
                                      'method_outcomes_inspected': False,
                                      'metadata_excluded': [c['case_id'] for c in cases if c['metadata_exclusion']]})
    print(json.dumps({'cases': len(cases), 'donors': len(donors), 'primary': len(primary),
                      'metadata_excluded': [c['case_id'] for c in cases if c['metadata_exclusion']],
                      'ages': [min(c['age_years'] for c in primary), max(c['age_years'] for c in primary)],
                      'notes': [{k: c[k] for k in ['case_id', 'pathology', 'notes']} for c in cases
                                if c['notes'] != 'NA' or c['pathology'] != 'NA']}, indent=2))


def check_selection():
    verify_freeze()
    f = json.loads((E / 'selection_freeze.json').read_text())
    assert sha(E / 'inventory.json') == f['inventory_sha256']
    return json.loads((E / 'inventory.json').read_text())['cases']


def prepare():
    import scipy  # fail before scheduling if the preparation runtime is incomplete
    cases = check_selection()
    assert not (E / 'preparation_status.json').exists()

    def job(c):
        case = c['case_id']
        out = E / 'inputs/prepared' / case
        out.mkdir(parents=True, exist_ok=True)
        if c['metadata_exclusion']:
            return {'case_id': case, 'status': 'METADATA_EXCLUDED', 'reason': c['metadata_exclusion']}
        source = E / c['source_relative']
        assert sha(source) == c['source_sha256']
        command = [sys.executable, str(H / 'prepare_mesh_v1.py'), '--case-id', case,
                   '--source', str(source), '--output', str(out / 'prepared.stl'),
                   '--report', str(out / 'preparation.json')]
        start = time.perf_counter()
        with (out / 'preparation.log').open('w') as log:
            proc = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT)
        result = {'case_id': case, 'returncode': proc.returncode, 'preparation_wall_seconds': time.perf_counter() - start}
        if proc.returncode:
            log = (out / 'preparation.log').read_text()
            result.update(status='SOURCE_REJECTED' if 'Fail-closed external source:' in log or 'Prepared STL failed post-export validation' in log else 'TECHNICAL_FAILURE', reason=log[-3000:])
        else:
            report = json.loads((out / 'preparation.json').read_text())
            cfg = E / 'configs' / (case + '.json')
            p = subprocess.run([sys.executable, str(H / 'make_case_config_v1.py'), '--case-id', case,
                                '--prepared-report', str(out / 'preparation.json'), '--profile', str(PROFILE),
                                '--output', str(cfg)], capture_output=True, text=True)
            if p.returncode:
                result.update(status='TECHNICAL_FAILURE', reason=p.stdout + p.stderr)
            else:
                result.update(status='PREPARED', longitudinal_length_mm=report['longitudinal_length_mm'],
                              prepared_sha256=sha(out / 'prepared.stl'), config_sha256=sha(cfg),
                              selected_faces=report['prepared_after_export']['faces'],
                              preparation_report_sha256=sha(out / 'preparation.json'))
        dump(out / 'status.json', result)
        print(json.dumps({k: result[k] for k in ['case_id', 'status']}), flush=True)
        return result

    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        results = list(pool.map(job, cases))
    dump(E / 'preparation_status.json', results)
    dump(E / 'execution_freeze.json', {'frozen_at_utc': datetime.now(timezone.utc).isoformat(),
                                     'preparation_status_sha256': sha(E / 'preparation_status.json'),
                                     'cases': results, 'method_outcomes_inspected': False})


def run(blender):
    cases = check_selection()
    f = json.loads((E / 'execution_freeze.json').read_text())
    assert sha(E / 'preparation_status.json') == f['preparation_status_sha256']
    statuses = {c['case_id']: c for c in f['cases']}
    assert not (E / 'execution_status.json').exists()

    def job(c):
        case = c['case_id']
        s = statuses[case]
        if s['status'] != 'PREPARED':
            return {'case_id': case, 'status': 'NOT_RUN', 'preparation_status': s['status']}
        source = E / 'inputs/prepared' / case / 'prepared.stl'
        config = E / 'configs' / (case + '.json')
        assert sha(source) == s['prepared_sha256'] and sha(config) == s['config_sha256']
        out = E / 'runs' / case
        start = time.perf_counter()
        proc = subprocess.run([sys.executable, str(R / 'code/run_checked.py'), '--blender', blender,
                               '--source', str(source), '--config', str(config), '--case-id', case,
                               '--out', str(out)], capture_output=True, text=True)
        result = {'case_id': case, 'returncode': proc.returncode, 'wall_seconds': time.perf_counter() - start}
        record = out / 'decision.json'
        if record.is_file():
            result['decision'] = json.loads(record.read_text())
            result['status'] = result['decision']['status']
        else:
            result.update(status='TECHNICAL_FAILURE', reason=(proc.stdout + proc.stderr)[-3000:])
        dump(out / 'launcher_status.json', result)
        print(json.dumps({'case_id': case, 'status': result['status'], 'wall_seconds': result['wall_seconds'],
                          'minimum_mm': result.get('decision', {}).get('continuous_minimum_mm')}), flush=True)
        return result

    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        results = list(pool.map(job, cases))
    dump(E / 'execution_status.json', results)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('phase', choices=['inventory', 'prepare', 'run'])
    p.add_argument('--blender', default='blender')
    args = p.parse_args()
    if args.phase == 'inventory': inventory()
    elif args.phase == 'prepare': prepare()
    else: run(args.blender)

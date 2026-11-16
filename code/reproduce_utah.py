"""Portable rerun of the fixed Utah cohort from the versioned source archive.

Requires external Python with numpy, trimesh, scipy and networkx, plus Blender.
The original acquisition/run logs remain separate from this reproduction route.
"""
from __future__ import annotations
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
from pathlib import Path
import subprocess
import sys
from zipfile import ZipFile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent


def digest_file(path, kind='sha256'):
    h = hashlib.new(kind)
    with path.open('rb') as f:
        for data in iter(lambda: f.read(1024 * 1024), b''): h.update(data)
    return h.hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--archive', type=Path, required=True)
    p.add_argument('--out', type=Path)
    p.add_argument('--blender', default='blender')
    p.add_argument('--workers', type=int, choices=[1,2], default=2)
    p.add_argument('--subset', choices=['all','primary','secondary_left'], default='all')
    p.add_argument('--case-id', action='append')
    p.add_argument('--dry-run', action='store_true')
    a = p.parse_args()
    local = ROOT / 'external_extension_2026_11_12/inventory.json'
    manifest = local if local.is_file() else ROOT / 'results/utah/inventory.json'
    inv = json.loads(manifest.read_text())
    assert digest_file(a.archive, 'md5') == inv['archive_md5']
    cases = [c for c in inv['cases'] if a.subset == 'all' or c['analysis'] == a.subset]
    if a.case_id:
        chosen = set(a.case_id); cases = [c for c in cases if c['case_id'] in chosen]
        assert {c['case_id'] for c in cases} == chosen
    print(json.dumps({'selected_bones':len(cases), 'source_defined_donors':len({c['donor_id'] for c in cases}),
                      'subset':a.subset, 'archive_md5_verified':True, 'dry_run':a.dry_run}), flush=True)
    if a.dry_run: return 0
    if a.out is None: p.error('--out is required except with --dry-run')
    out = a.out.resolve()
    if out.exists() and any(out.iterdir()): p.error('Output directory must be empty.')
    import numpy, scipy, trimesh, networkx
    out.mkdir(parents=True, exist_ok=True)
    history = HERE / 'historical'
    if not history.is_dir(): history = ROOT.parent / 'benchmark/external_validation_v1/scripts'
    profile = ROOT / 'configs/EXTERNAL_PARAMETER_PROFILE_v1.json'
    if not profile.is_file(): profile = ROOT.parent / 'benchmark/external_validation_v1/protocol/EXTERNAL_PARAMETER_PROFILE_v1.json'
    assert digest_file(profile) == '6f19a8831457d25a56ef63b99cde3b3443b9d0e8023241d71ab7f03eedafa5fd'
    for filename, expected in [('prepare_mesh_v1.py','69bb7573930123ee01f0a8211f702f4286e9d1461d837ed547fb3849bf88592a'),
                               ('make_case_config_v1.py','bd929eb041f75d54da662d48627fa10ba862a47d34f12d18541df816e50036c8')]:
        assert digest_file(history / filename) == expected

    def job(c):
        case = c['case_id']; d = out / case; d.mkdir()
        result = {'case_id':case, 'analysis':c['analysis'], 'donor_id':c['donor_id'], 'accepted':False}
        if c['metadata_exclusion']:
            result.update(status='METADATA_EXCLUDED', reason=c['metadata_exclusion'])
        else:
            with ZipFile(a.archive) as z: raw = z.read(c['archive_member'])
            assert hashlib.sha256(raw).hexdigest() == c['source_sha256']
            (d / 'source.stl').write_bytes(raw)
            commands = [
                [sys.executable,str(history / 'prepare_mesh_v1.py'),'--case-id',case,'--source',str(d / 'source.stl'),'--output',str(d / 'prepared.stl'),'--report',str(d / 'preparation.json')],
                [sys.executable,str(history / 'make_case_config_v1.py'),'--case-id',case,'--prepared-report',str(d / 'preparation.json'),'--profile',str(profile),'--output',str(d / 'config.json')],
                [sys.executable,str(HERE / 'run_checked.py'),'--blender',a.blender,'--source',str(d / 'prepared.stl'),'--config',str(d / 'config.json'),'--case-id',case,'--out',str(d / 'evaluation')]]
            for step, cmd in zip(['prepare','configure','evaluate'], commands):
                with (d / (step + '.log')).open('w') as log:
                    proc = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT)
                if proc.returncode:
                    result.update(status='STOPPED', step=step, returncode=proc.returncode)
                    break
            decision = d / 'evaluation/decision.json'
            if decision.is_file(): result.update(json.loads(decision.read_text()))
        (d / 'result.json').write_text(json.dumps(result, indent=2)+'\n')
        print(json.dumps({k:result[k] for k in ['case_id','accepted','status']}), flush=True)
        return result

    with ThreadPoolExecutor(max_workers=a.workers) as pool: results = list(pool.map(job, cases))
    (out / 'results.json').write_text(json.dumps(results, indent=2)+'\n')
    return 0 if all(r['accepted'] for r in results) else 2


if __name__ == '__main__': raise SystemExit(main())

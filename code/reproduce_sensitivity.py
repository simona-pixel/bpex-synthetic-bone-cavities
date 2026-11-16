"""Rerun the 90 archived configurations through both checking stages.

Prepared meshes must be acquired and aligned using the documented source
protocol. They are not redistributed. This launcher creates a new run and
does not overwrite the recorded study results.
"""
from __future__ import annotations
import argparse
import concurrent.futures
import json
import subprocess
import sys
from pathlib import Path


def main():
    root = Path(__file__).resolve().parents[1]
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--prepared-inputs', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--blender', default='blender')
    p.add_argument('--workers', type=int, default=2)
    p.add_argument('--case-id', help='Optional single case')
    p.add_argument('--variant', help='Optional single configuration stem')
    p.add_argument('--dry-run', action='store_true')
    a = p.parse_args()
    configs = sorted((root / 'configs/sensitivity').glob('*/*.json'))
    configs = [c for c in configs if (not a.case_id or c.parent.name == a.case_id)
               and (not a.variant or c.stem == a.variant)]
    if not configs:
        p.error('No matching configuration')
    commands = []
    for c in configs:
        case = c.parent.name
        source = a.prepared_inputs / case / 'prepared.stl'
        if not source.is_file():
            source = a.prepared_inputs / (case + '.stl')
        if not source.is_file():
            p.error('Missing prepared mesh for ' + case)
        commands.append((case, c.stem, [sys.executable, str(root / 'code/run_checked.py'),
            '--blender', a.blender, '--source', str(source.resolve()),
            '--config', str(c.resolve()), '--case-id', case,
            '--out', str((a.out / case / c.stem).resolve())]))
    if a.dry_run:
        print(json.dumps([{'case_id': c, 'variant': v, 'command': cmd}
                          for c, v, cmd in commands], indent=2))
        return 0
    if a.out.exists() and any(a.out.iterdir()):
        p.error('Output directory must be empty')
    a.out.mkdir(parents=True, exist_ok=True)

    def run(item):
        case, variant, cmd = item
        r = subprocess.run(cmd, text=True, capture_output=True)
        out = a.out / case / variant
        out.mkdir(parents=True, exist_ok=True)
        (out / 'launcher_stdout.txt').write_text(r.stdout)
        (out / 'launcher_stderr.txt').write_text(r.stderr)
        return {'case_id': case, 'variant': variant, 'returncode': r.returncode}

    with concurrent.futures.ThreadPoolExecutor(max_workers=max(1, a.workers)) as pool:
        rows = list(pool.map(run, commands))
    (a.out / 'execution_summary.json').write_text(json.dumps(rows, indent=2) + '\n')
    print(json.dumps(rows, indent=2))
    return 1 if any(r['returncode'] not in (0, 2) for r in rows) else 0


if __name__ == '__main__':
    raise SystemExit(main())

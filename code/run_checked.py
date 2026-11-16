"""Portable command-line entry point requiring sampled and continuous checks."""
from __future__ import annotations
import argparse,json,shutil,subprocess,sys
from pathlib import Path

def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--blender',default='blender');p.add_argument('--source',required=True);p.add_argument('--config',required=True);p.add_argument('--case-id',required=True);p.add_argument('--out',required=True);a=p.parse_args()
 out=Path(a.out).resolve()
 if out.exists() and any(out.iterdir()):p.error('Output directory must be empty; previous results are never overwritten.')
 out.mkdir(parents=True,exist_ok=True)
 command=[a.blender,'--background','--factory-startup','--threads','2','--python-exit-code','1','--python',str(Path(__file__).with_name('checked_worker_blender.py')),'--','--source',str(Path(a.source).resolve()),'--config',str(Path(a.config).resolve()),'--case-id',a.case_id,'--out',str(out)]
 (out/'command.json').write_text(json.dumps(command,indent=2)+'\n')
 with (out/'execution.log').open('w') as log:r=subprocess.run(command,stdout=log,stderr=subprocess.STDOUT)
 record=out/'decision.json'
 if not record.is_file():
  print(json.dumps({'status':'TECHNICAL_FAILURE','returncode':r.returncode,'log':str(out/'execution.log')}));return 1
 decision=json.loads(record.read_text());print(json.dumps(decision,indent=2))
 return 0 if decision['accepted'] and r.returncode==0 else 2
if __name__=='__main__':raise SystemExit(main())

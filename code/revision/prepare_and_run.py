"""Freeze and run the predeclared exploratory BPEX revision analysis.

Run using an external Python interpreter, not Blender. Historical inputs are
read only; all products are written to this revision directory.
"""
from __future__ import annotations
import argparse, concurrent.futures, csv, hashlib, json, shutil, subprocess, sys, time
from datetime import datetime, timezone
import os
from pathlib import Path

REV = Path(__file__).resolve().parents[2]
ARTICLE = Path(os.environ.get('BPEX_STUDY_ROOT', REV.parent))
BENCH = ARTICLE / 'benchmark'
BLENDER = os.environ.get('BLENDER_BIN', 'blender')

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(p, obj):
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(obj, indent=2, ensure_ascii=False) + '\n')

def prepare():
    core = REV / 'code/core'; core.mkdir(parents=True, exist_ok=True)
    source_paths = [BENCH/'pipeline_v3.py', BENCH/'revision_validation_v1/exact_surface_clearance_blender.py']
    for p in source_paths: shutil.copy2(p, core/p.name)
    selected = [
      ('20260825_vsd_v3_holdout_v1','vsd_002_femur_r','002','femur'),
      ('20260825_vsd_v3_holdout_v1','vsd_z064_femur_r','z064','femur'),
      ('20260825_vsd_v3_holdout_v1','vsd_006_tibia_r','006','tibia'),
      ('20260825_vsd_v3_holdout_v1','vsd_z062_tibia_r','z062','tibia'),
      ('20260825_zenodo_multi_os_v3_development','zenodo_6614647_humerus_r','ABH074','humerus'),
      ('20260825_zenodo_multi_os_v3_development','zenodo_6614647_ulna_r','ABH074','ulna'),
    ]
    variants = [('baseline',None,None)]
    for field, vals in [
      ('canal_ratio_isthme',[0.40,0.60]),('canal_ratio_diaphyse',[0.50,0.70]),
      ('canal_ratio_metaphyse',[0.74,0.94]),
      ('validated_loft_radial_safety_mm',[0.20,0.60]),
      ('containment_guard_minimum_clearance_mm',[0.075,0.300]),
      ('validated_loft_z_margin_fraction',[0.20,0.10]),
      ('resolution',[(201,90),(801,360)])]:
        for val in vals: variants.append((field+'_'+str(val).replace(' ','').replace('(','').replace(')','').replace(',','x'),field,val))
    cases=[]
    for run, case, subj, bone in selected:
        src = BENCH/'external_validation_v1/runs'/run/'cases'/case
        inp = REV/'analysis/inputs'/case; inp.mkdir(parents=True,exist_ok=True)
        for name in ['prepared.stl','prepared_report.json','case_config.json']:
            shutil.copy2(src/name,inp/name)
        base=json.loads((inp/'case_config.json').read_text())
        case_jobs=[]
        for label,field,value in variants:
            cfg=json.loads(json.dumps(base))
            cfg['revision_analysis_role']='exploratory parameter sensitivity; reused subjects'
            if field=='resolution':
                cfg['benchmark']['validated_loft_sections'],cfg['benchmark']['validated_loft_angles']=value
            elif field and field.startswith('canal_ratio'):
                cfg['bones'][case]['parameters'][field]=value
            elif field: cfg['benchmark'][field]=value
            out=REV/'analysis/runs'/case/label
            dump(out/'config.json',cfg)
            case_jobs.append({'label':label,'field':field,'value':value,'config':str(out/'config.json'),'output':str(out)})
        cases.append({'case_id':case,'subject_id':subj,'bone':bone,'historical_run':run,'prepared':str(inp/'prepared.stl'),'input_sha256':sha(inp/'prepared.stl'),'historical_output':str(src/'rep_1/output'),'jobs':case_jobs})
    freeze={'frozen_at_utc':datetime.now(timezone.utc).isoformat(),'role':'exploratory post-review analysis; not an independent validation','cases':cases,'case_count':6,'donor_count':5,'configurations_per_case':len(variants),'historical_code':{p.name:sha(p) for p in source_paths},'historical_profile_sha256':sha(BENCH/'external_validation_v1/protocol/EXTERNAL_PARAMETER_PROFILE_v1.json'),'protocol_sha256':sha(REV/'PLAN_REVISION.md'),'runner_sha256':sha(__file__)}
    dump(REV/'analysis/FROZEN_PLAN.json',freeze)
    return freeze

def execute(case):
    log=REV/'analysis/logs'/(case['case_id']+'.log'); log.parent.mkdir(parents=True,exist_ok=True)
    started=time.perf_counter()
    with log.open('w') as f:
        cp=subprocess.run([BLENDER,'--background','--factory-startup','--threads','2','--python',str(REV/'code/run_case_blender.py'),'--','--plan',str(REV/'analysis/FROZEN_PLAN.json'),'--case',case['case_id']],stdout=f,stderr=subprocess.STDOUT,timeout=14400)
    return {'case_id':case['case_id'],'returncode':cp.returncode,'wall_seconds':time.perf_counter()-started,'log':str(log)}

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--prepare-only',action='store_true');ap.add_argument('--resume',action='store_true');ap.add_argument('--workers',type=int,default=2);a=ap.parse_args()
    freeze=json.loads((REV/'analysis/FROZEN_PLAN.json').read_text()) if a.resume else prepare()
    if a.prepare_only: print('Frozen',len(freeze['cases']),'cases x',freeze['configurations_per_case'],'configurations',flush=True);return
    results=[]
    with concurrent.futures.ThreadPoolExecutor(max_workers=a.workers) as pool:
        for future in concurrent.futures.as_completed([pool.submit(execute,c) for c in freeze['cases']]):
            result=future.result(); results.append(result); print(json.dumps(result),flush=True)
    dump(REV/'analysis/execution_status.json',results)
    if any(x['returncode'] for x in results): raise SystemExit(1)
if __name__=='__main__':main()

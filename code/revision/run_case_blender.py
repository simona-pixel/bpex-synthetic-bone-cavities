"""Execute one public case and all frozen configurations, retaining failures."""
from __future__ import annotations
import argparse, importlib.util, json, sys, time, traceback
import os
from pathlib import Path
from types import SimpleNamespace
import bpy

HERE=Path(__file__).resolve().parent
def load(name,p):
    spec=importlib.util.spec_from_file_location(name,p);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
pipe=load('historical_pipeline_v3',HERE/'core/pipeline_v3.py')
audit=load('continuous_audit',HERE/'core/exact_surface_clearance_blender.py')
def dump(p,o):p.write_text(json.dumps(o,indent=2,allow_nan=False)+'\n')

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--plan');ap.add_argument('--case');a=ap.parse_args(sys.argv[sys.argv.index('--')+1:])
    plan=json.loads(Path(a.plan).read_text());case=next(c for c in plan['cases'] if c['case_id']==a.case)
    results=[]
    for job in case['jobs']:
        out=Path(job['output']); result_file=out/'result.json'
        if result_file.exists():
            print('REUSE',a.case,job['label'],flush=True);results.append(json.loads(result_file.read_text()));continue
        result={'case_id':a.case,'subject_id':case['subject_id'],'bone':case['bone'],'variant':job['label'],'field':job['field'],'value':job['value'],'status':'STARTED','blender':bpy.app.version_string,'blender_python':sys.version}
        start=time.perf_counter()
        try:
            args=SimpleNamespace(config=job['config'],bone=a.case,method='regional',source=case['prepared'],output=str(out),thickness=None,regional_postprocess='none',regional_kernel='validated_loft',allow_invalid_source=False,allow_invalid_output=True,containment_guard=True)
            meta=pipe.run(args)
            cfg=json.loads(Path(job['config']).read_text());threshold=cfg['benchmark']['containment_guard_minimum_clearance_mm']
            result.update({'historical_criteria_pass':meta['output_all_valid'],'volume_mm3':meta['inner']['volume_mm3'],'sampled_minimum_mm':meta['relation_diagnostics']['nearest_unsigned_min_mm'],'radial_scale':meta['containment_guard']['selected_radial_scale'],'generation_internal_seconds':meta['timing_seconds']['total_internal'],'outer_path':meta['outer_path'],'inner_path':meta['inner_path'],'inner_signature':meta['inner_mesh_signature'],'threshold_mm':threshold})
            # A changed acceptance threshold never substitutes for the fixed 0.150 mm reference.
            if meta['output_geometry_valid']:
                check=audit.audit(Path(meta['outer_path']),Path(meta['inner_path']),a.case+'__'+job['label'],threshold)
                dump(out/'continuous_audit.json',check)
                result.update({'continuous_minimum_mm':check['continuous_triangle_surface_minimum_mm'],'continuous_audit_seconds':check['elapsed_seconds'],'final_criteria_pass':check['continuous_clearance_pass'] and meta['output_all_valid'],'fixed_0150_pass':check['continuous_triangle_surface_minimum_mm']+1e-9>=.15,'outer_sha256':check['outer_sha256'],'inner_sha256':check['inner_sha256']})
            else: result.update({'final_criteria_pass':False,'fixed_0150_pass':False})
            result['status']='COMPLETED'
        except Exception as exc:
            result.update({'status':'FAILED','error':str(exc),'traceback':traceback.format_exc(),'final_criteria_pass':False})
            print(traceback.format_exc(),flush=True)
        result['wall_seconds']=time.perf_counter()-start
        dump(result_file,result);results.append(result)
        print('REVISION_RESULT='+json.dumps({k:result.get(k) for k in ['case_id','variant','status','final_criteria_pass','continuous_minimum_mm','wall_seconds']}),flush=True)
    dump(Path(case['jobs'][0]['output']).parent/'case_summary.json',results)
    print('CASE_COMPLETE',a.case,len(results),flush=True)
if __name__=='__main__':main()

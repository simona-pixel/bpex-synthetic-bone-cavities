"""Blender worker. Candidate files never imply final acceptance."""
import argparse,hashlib,importlib.util,json,shutil,sys,traceback
from pathlib import Path
from types import SimpleNamespace
HERE=Path(__file__).resolve().parent
def load(name,p):
 s=importlib.util.spec_from_file_location(name,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 p=argparse.ArgumentParser();p.add_argument('--source');p.add_argument('--config');p.add_argument('--case-id');p.add_argument('--out');a=p.parse_args(sys.argv[sys.argv.index('--')+1:]);out=Path(a.out)
 pipe=load('pipeline',HERE/'core/pipeline_v3.py');aud=load('auditor',HERE/'core/exact_surface_clearance_blender.py');decision={'case_id':a.case_id,'accepted':False,'status':'STARTED','source_sha256':sha(a.source),'config_sha256':sha(a.config),'pipeline_sha256':sha(HERE/'core/pipeline_v3.py'),'auditor_sha256':sha(HERE/'core/exact_surface_clearance_blender.py')}
 try:
  cfg=json.loads(Path(a.config).read_text());candidate=out/'candidate'
  args=SimpleNamespace(config=a.config,bone=a.case_id,method='regional',source=a.source,output=str(candidate),thickness=None,regional_postprocess='none',regional_kernel='validated_loft',allow_invalid_source=False,allow_invalid_output=False,containment_guard=True)
  meta=pipe.run(args);threshold=float(cfg['benchmark']['containment_guard_minimum_clearance_mm']);r=aud.audit(Path(meta['outer_path']),Path(meta['inner_path']),a.case_id,threshold)
  (out/'continuous_audit.json').write_text(json.dumps(r,indent=2)+'\n');decision['continuous_minimum_mm']=r['continuous_triangle_surface_minimum_mm'];decision['threshold_mm']=threshold
  decision['accepted']=bool(meta['output_all_valid'] and r['continuous_clearance_pass']);decision['status']='ACCEPTED' if decision['accepted'] else 'WITHHELD_CONTINUOUS_CLEARANCE'
  if decision['accepted']:
   target=out/'accepted';target.mkdir()
   for key,name in [('outer_path','outer.stl'),('inner_path','cavity.stl')]:shutil.copy2(meta[key],target/name)
 except Exception as exc:
  decision.update({'status':'WITHHELD_OR_ERROR','reason':str(exc),'traceback':traceback.format_exc()})
 (out/'decision.json').write_text(json.dumps(decision,indent=2)+'\n')
 print('CHECKED_DECISION='+json.dumps({k:decision[k] for k in ['case_id','accepted','status']}),flush=True)
if __name__=='__main__':main()

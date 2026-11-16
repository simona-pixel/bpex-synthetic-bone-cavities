"""Build analysis tables from executed records without pooling reused subjects."""
import collections,csv,hashlib,json,platform,statistics,sys
import numpy as np
import os
from pathlib import Path
R=Path(__file__).resolve().parents[2]
def dump(p,o):p.write_text(json.dumps(o,indent=2,ensure_ascii=False)+'\n')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def csvout(p,rows):
 fields=list(dict.fromkeys(k for r in rows for k in r))
 with p.open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields,extrasaction='raise');w.writeheader();w.writerows(rows)
plan=json.loads((R/'analysis/FROZEN_PLAN.json').read_text());rows=[];replay=[]
for c in plan['cases']:
 base=json.loads((Path(c['jobs'][0]['output'])/'result.json').read_text())
 original=Path(c['historical_output'])/(c['case_id']+'_inner_negative.stl')
 a=original.read_bytes();b=Path(base['inner_path']).read_bytes();stl_dtype=np.dtype([('n','<f4',(3,)),('v','<f4',(3,3)),('attribute','<u2')]);oldtri=np.frombuffer(a,stl_dtype,offset=84);newtri=np.frombuffer(b,stl_dtype,offset=84)
 replay.append({'case_id':c['case_id'],'historical_inner_sha256':sha(original),'revision_baseline_inner_sha256':base['inner_sha256'],'byte_identical':a==b,'triangle_coordinates_and_order_identical':bool(np.array_equal(oldtri['v'],newtri['v'])),'maximum_coordinate_difference_mm':float(np.max(np.abs(oldtri['v']-newtri['v']))),'maximum_stored_normal_component_difference':float(np.max(np.abs(oldtri['n']-newtri['n']))),'headers_identical':a[:80]==b[:80]})
 for job in c['jobs']:
  x=json.loads((Path(job['output'])/'result.json').read_text())
  if x['status']=='FAILED' and 'containment guard could not find a valid radial scale' in x.get('error',''):
   x['status']='METHOD_REJECTED'
  if x.get('continuous_minimum_mm') is not None:x['fixed_0150_pass']=x['continuous_minimum_mm']+1e-9>=.15 and x['historical_criteria_pass']
  x['volume_percent_change_from_baseline']=100*(x['volume_mm3']/base['volume_mm3']-1) if 'volume_mm3' in x else None
  rows.append({k:v for k,v in x.items() if k not in ['traceback','outer_path','inner_path','blender_python','inner_signature']})
if len(rows)!=90:raise ValueError(f'Expected 90 results, found {len(rows)}')
groups=[]
for variant in dict.fromkeys(x['variant'] for x in rows):
 a=[x for x in rows if x['variant']==variant];d=[x['continuous_minimum_mm'] for x in a if x.get('continuous_minimum_mm') is not None];vol=[x['volume_percent_change_from_baseline'] for x in a if x.get('volume_percent_change_from_baseline') is not None]
 groups.append({'variant':variant,'attempted':len(a),'completed':sum(x['status']=='COMPLETED' for x in a),'historical_pass':sum(bool(x.get('historical_criteria_pass')) for x in a),'final_pass_requested_threshold':sum(bool(x.get('final_criteria_pass')) for x in a),'final_pass_fixed_0150':sum(bool(x.get('fixed_0150_pass')) for x in a),'continuous_minimum_mm':min(d) if d else None,'continuous_maximum_mm':max(d) if d else None,'volume_change_min_percent':min(vol) if vol else None,'volume_change_max_percent':max(vol) if vol else None,'contractions':sum(x.get('radial_scale',1)<1 for x in a)})
res=json.loads((R/'analysis/resolution_comparison.json').read_text());wind=json.loads((R/'analysis/qualification/independent_winding_summary.json').read_text());fixture=json.loads((R/'analysis/qualification/threshold_fixtures.json').read_text())
summary={'role':plan['role'],'case_count':6,'donor_count':5,'executions':len(rows),'completed':sum(x['status']=='COMPLETED' for x in rows),'method_rejections':sum(x['status']=='METHOD_REJECTED' for x in rows),'technical_failures':sum(x['status']=='FAILED' for x in rows),'historical_passes':sum(bool(x.get('historical_criteria_pass')) for x in rows),'final_passes':sum(bool(x.get('final_criteria_pass')) for x in rows),'baseline_replay_all_byte_identical':all(x['byte_identical'] for x in replay),'all_continuous_minimum_mm':min(x['continuous_minimum_mm'] for x in rows if 'continuous_minimum_mm' in x),'smallest_passing_clearance_mm':min(x['continuous_minimum_mm'] for x in rows if x['final_criteria_pass']),'groups':groups,'failures':[x for x in rows if not x['final_criteria_pass']],'resolution':res,'winding':{'total_queries':sum(x['queries'] for x in wind),'known_analytic_queries':sum(x['known_expected_queries'] for x in wind),'parity_winding_disagreements':sum(x['parity_winding_disagreements'] for x in wind),'analytic_errors':sum(x['winding_expected_errors']+x['parity_expected_errors'] for x in wind),'nearest_face_winding_disagreements':sum(x['nearest_face_winding_disagreements'] for x in wind)},'threshold_fixture_count':len(fixture),'threshold_expected_decisions_match':all(x['expected_continuous_pass']==x['observed_pass'] for x in fixture),'threshold_maximum_abs_error_mm':max(x['absolute_error_mm'] for x in fixture),'new_environment':{'external_python':sys.version,'platform':platform.platform(),'blender':rows[0]['blender'],'blender_python':json.loads((Path(plan['cases'][0]['jobs'][0]['output'])/'result.json').read_text())['blender_python'],'blender_threads':2,'concurrent_sensitivity_processes':2,'repetitions_per_configuration':1}}
dump(R/'analysis/sensitivity_summary.json',summary);dump(R/'analysis/baseline_replay.json',replay);dump(R/'analysis/results_derived.json',rows);csvout(R/'analysis/sensitivity_results.csv',rows);csvout(R/'analysis/sensitivity_by_configuration.csv',groups)
print(json.dumps({k:v for k,v in summary.items() if k not in ['groups','failures','resolution','new_environment']},indent=2))
print(json.dumps(groups,indent=2))

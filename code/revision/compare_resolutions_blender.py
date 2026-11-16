"""Bidirectional sampled vertex-to-surface comparisons of resolution variants."""
import importlib.util,json
import os
from pathlib import Path
import numpy as np
from mathutils.bvhtree import BVHTree
R=Path(__file__).resolve().parents[2]
s=importlib.util.spec_from_file_location('audit',R/'code/core/exact_surface_clearance_blender.py');m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
plan=json.loads((R/'analysis/FROZEN_PLAN.json').read_text());rows=[]
for case in plan['cases']:
 results={j['label']:json.loads((Path(j['output'])/'result.json').read_text()) for j in case['jobs']}
 fine=results['resolution_801x360'];fv,ff=m.read_stl(Path(fine['inner_path']));ft=BVHTree.FromPolygons(fv,ff,all_triangles=True)
 for label in ['resolution_201x90','baseline']:
  r=results[label];v,f=m.read_stl(Path(r['inner_path']));tree=BVHTree.FromPolygons(v,f,all_triangles=True)
  d1=np.array([ft.find_nearest(p)[3] for p in v]);d2=np.array([tree.find_nearest(p)[3] for p in fv])
  rows.append({'case_id':case['case_id'],'comparison':label+'_versus_801x360','variant_vertices':len(v),'reference_vertices':len(fv),'variant_to_fine_mean_mm':float(d1.mean()),'fine_to_variant_mean_mm':float(d2.mean()),'bidirectional_vertex_max_mm':float(max(d1.max(),d2.max())),'mean_of_directional_means_mm':float((d1.mean()+d2.mean())/2),'volume_percent_difference_from_fine':100*(r['volume_mm3']/fine['volume_mm3']-1),'clearance_difference_from_fine_mm':r['continuous_minimum_mm']-fine['continuous_minimum_mm'],'scope':'vertex-sampled distance; not a certified continuous Hausdorff distance'})
 print('RESOLUTION_COMPARED',case['case_id'],flush=True)
(R/'analysis/resolution_comparison.json').write_text(json.dumps(rows,indent=2)+'\n')

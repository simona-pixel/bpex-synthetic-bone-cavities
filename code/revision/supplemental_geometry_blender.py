"""Generate independent-predicate queries, threshold fixtures and a hollow coupon."""
import importlib.util,json,sys
import os
from pathlib import Path
import bpy,bmesh,numpy as np
from mathutils import Vector
R=Path(__file__).resolve().parents[2]
def load(name,p):
 s=importlib.util.spec_from_file_location(name,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
pipe=load('core',R/'code/core/pipeline_v3.py');aud=load('audit',R/'code/core/exact_surface_clearance_blender.py')
O=R/'analysis/qualification';O.mkdir(parents=True,exist_ok=True)
def dump(p,o):p.write_text(json.dumps(o,indent=2)+'\n')
def obj(name,v,f):
 m=bpy.data.meshes.new(name);m.from_pydata(v,[],f);m.update();o=bpy.data.objects.new(name,m);bpy.context.collection.objects.link(o);pipe.recalc_normals_outside(o);return o
def cube(name,half):
 bpy.ops.mesh.primitive_cube_add(size=2*half);o=bpy.context.object;o.name=name;bm=bmesh.new();bm.from_mesh(o.data);bmesh.ops.triangulate(bm,faces=list(bm.faces));bm.to_mesh(o.data);bm.free();pipe.recalc_normals_outside(o);return o
def lprism():
 pts=[(0,0),(6,0),(6,2),(2,2),(2,6),(0,6)];v=[(x,y,z) for z in [-3,3] for x,y in pts];t=[(0,1,3),(1,2,3),(0,3,5),(3,4,5)];f=[tuple(reversed(q)) for q in t]+[tuple(i+6 for i in q) for q in t]
 for i in range(6):j=(i+1)%6;f.extend([(i,j,j+6),(i,j+6,i+6)])
 return obj('L_prism',v,f)
def tunnel():
 outer=[(-5,-5),(5,-5),(5,5),(-5,5)];hole=[(-2,-2),(2,-2),(2,2),(-2,2)];v=[(x,y,z) for z in [-3,3] for points in [outer,hole] for x,y in points];f=[]
 for i in range(4):
  j=(i+1)%4;f.extend([(i,j,8+j),(i,8+j,8+i),(4+i,12+i,12+j),(4+i,12+j,4+j),(8+i,8+j,12+j),(8+i,12+j,12+i),(i,4+i,4+j),(i,4+j,j)])
 return obj('square_tunnel',v,f)
def arrays(o):
 return np.array([tuple(o.matrix_world@v.co) for v in o.data.vertices]),np.array([tuple(p.vertices) for p in o.data.polygons],dtype=int)
def save_queries(name,o,q,expected,kinds):
 v,f=arrays(o);bm,tree=pipe.bvh_for(o);maxd=float(np.linalg.norm(v.max(0)-v.min(0))*4)
 votes=[pipe.point_inside_votes(tree,Vector(p),maxd) for p in q]
 nearest=[]
 for p in q:
  hit=tree.find_nearest(Vector(p));nearest.append((Vector(p)-hit[0]).dot(hit[1])<0)
 bm.free();np.savez_compressed(O/(name+'.npz'),vertices=v,faces=f,queries=np.array(q),ray_votes=votes,nearest_face_inside=nearest,expected=np.array(expected),kinds=np.array(kinds))
def threshold_tests():
 rows=[]
 for d in [.149,.1498,.14995,.15,.15005,.1502,.151]:
  pipe.clear_scene();outer=cube('outer',5);inner=cube('inner',5-d)
  op=O/'fixture_outer.stl';ip=O/('fixture_inner_'+str(d)+'.stl');pipe.export_stl(outer,str(op));pipe.export_stl(inner,str(ip))
  a=aud.audit(op,ip,'clearance_'+str(d),.15);rows.append({'nominal_mm':d,'continuous_mm':a['continuous_triangle_surface_minimum_mm'],'expected_continuous_pass':d>=.15,'expected_sampled_guard_pass':d>=.15-.0001,'observed_pass':a['continuous_clearance_pass'],'absolute_error_mm':abs(d-a['continuous_triangle_surface_minimum_mm'])})
 dump(O/'threshold_fixtures.json',rows)
def analytic_queries():
 pts=[(x,y,z) for x in [-6.25,-4.75,-3.25,-1.75,-.25,1.25,2.75,4.25,5.75] for y in [-6.25,-4.75,-3.25,-1.75,-.25,1.25,2.75,4.25,5.75] for z in [-4.5,-1.5,1.5,4.5]]
 for name,make,pred in [
 ('cube',lambda:cube('cube',5),lambda x,y,z:abs(x)<5 and abs(y)<5 and abs(z)<5),
 ('L_prism',lprism,lambda x,y,z:abs(z)<3 and ((0<x<6 and 0<y<2) or (0<x<2 and 0<y<6))),
 ('square_tunnel',tunnel,lambda x,y,z:abs(z)<3 and abs(x)<5 and abs(y)<5 and not(abs(x)<2 and abs(y)<2))]:
  pipe.clear_scene();o=make();save_queries('analytic_'+name,o,pts,[pred(*p) for p in pts],['analytic_known']*len(pts))
def external_queries(plan):
 for case in plan['cases']:
  pipe.clear_scene();base=Path(case['jobs'][0]['output']);res=json.loads((base/'result.json').read_text())
  outer=pipe.import_stl(res['outer_path'],'outer');inner=pipe.import_stl(res['inner_path'],'inner');v,f=arrays(outer);iv,_=arrays(inner)
  q=list(iv[np.linspace(0,len(iv)-1,128,dtype=int)]);kinds=['generated_vertex']*128;expected=[-1]*128
  for idx in np.linspace(0,len(f)-1,64,dtype=int):
   tri=v[f[idx]];n=np.cross(tri[1]-tri[0],tri[2]-tri[0]);n/=np.linalg.norm(n);c=tri.mean(0)
   for sign,kind in [(-1,'near_surface_inward'),(1,'near_surface_outward')]:
    q.append(c+sign*.01*n);kinds.append(kind);expected.append(-1)
  save_queries(case['case_id'],outer,q,expected,kinds)
def apply_boolean(target,cutter,operation):
 pipe.select_only(target);m=target.modifiers.new('demonstration_'+operation,'BOOLEAN');m.operation=operation;m.solver='EXACT';m.object=cutter;bpy.ops.object.modifier_apply(modifier=m.name);bpy.data.objects.remove(cutter,do_unlink=True)
def digital_coupon(plan,explicit_case=None,explicit_result=None):
 case=explicit_case or plan['cases'][0];res=explicit_result or json.loads((Path(case['jobs'][0]['output'])/'result.json').read_text());pipe.clear_scene()
 outer=pipe.import_stl(res['outer_path'],'outer');inner=pipe.import_stl(res['inner_path'],'cavity');v,_=arrays(outer);z0=float((v[:,2].min()+v[:,2].max())/2);zlo,zhi=z0-20,z0+20
 bpy.ops.mesh.primitive_cube_add(size=2,location=((v[:,0].min()+v[:,0].max())/2,(v[:,1].min()+v[:,1].max())/2,z0));box=bpy.context.object;box.dimensions=(float(np.ptp(v[:,0])+20),float(np.ptp(v[:,1])+20),40);bpy.ops.object.transform_apply(location=False,rotation=False,scale=True)
 apply_boolean(outer,box,'INTERSECT');solid_volume=pipe.mesh_diagnostics(outer)['volume_mm3'];apply_boolean(outer,inner,'DIFFERENCE')
 bm=bmesh.new();bm.from_mesh(outer.data);bmesh.ops.triangulate(bm,faces=list(bm.faces));bm.to_mesh(outer.data);bm.free();pipe.recalc_normals_outside(outer)
 diagnostics=pipe.mesh_diagnostics(outer);vv,ff=arrays(outer);euler=int(len(vv)-diagnostics['edges']+len(ff));crossings=pipe.self_intersection_count(outer)
 path=O/'digital_hollow_diaphyseal_coupon.stl';pipe.export_stl(outer,str(path));np.savez_compressed(O/'digital_coupon_mesh.npz',vertices=vv,faces=ff)
 # A connected orientable closed surface with Euler characteristic zero has
 # genus one. In this construction that handle is the through cavity.
 report={'case_id':case['case_id'],'length_mm':40,'z_limits_mm':[zlo,zhi],'solid_segment_volume_mm3':solid_volume,'hollow_material_volume_mm3':diagnostics['volume_mm3'],'removed_volume_mm3':solid_volume-diagnostics['volume_mm3'],'diagnostics':diagnostics,'non_adjacent_self_intersection_pairs':crossings,'euler_characteristic':euler,'genus_one_through_cavity':euler==0 and diagnostics['components']==1 and diagnostics['closed'],'parent_pair_clearance_mm':res['continuous_minimum_mm'],'scope':'digital Boolean and topology demonstration only; no printing, instrumentation or physical validation','output_stl':str(path)}
 if pipe.diagnostic_failures(diagnostics) or crossings or euler!=0:raise RuntimeError('Coupon failed the intended digital checks: '+json.dumps(report))
 dump(O/'digital_coupon_report.json',report)
if __name__=='__main__':
 plan=json.loads((R/'analysis/FROZEN_PLAN.json').read_text());threshold_tests();analytic_queries();external_queries(plan);digital_coupon(plan);print('SUPPLEMENTAL_GEOMETRY_COMPLETE',flush=True)

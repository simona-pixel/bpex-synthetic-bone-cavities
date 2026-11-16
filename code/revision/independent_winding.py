"""Independent solid-angle sum in external NumPy Python, no Blender predicates.

Implements the oriented solid-angle definition of generalized winding number.
This targets closed consistently oriented triangulations, not adversarial soups.
"""
import json,time
import os
from pathlib import Path
import numpy as np
R=Path(__file__).resolve().parents[2];O=R/'analysis/qualification'
def winding(v,f,q):
 tri=v[f];out=[]
 for start in range(0,len(q),8):
  a=tri[None,:,0,:]-q[start:start+8,None,:];b=tri[None,:,1,:]-q[start:start+8,None,:];c=tri[None,:,2,:]-q[start:start+8,None,:]
  la=np.linalg.norm(a,axis=2);lb=np.linalg.norm(b,axis=2);lc=np.linalg.norm(c,axis=2)
  num=np.sum(a*np.cross(b,c),axis=2)
  den=la*lb*lc+np.sum(a*b,axis=2)*lc+np.sum(b*c,axis=2)*la+np.sum(c*a,axis=2)*lb
  out.extend(np.sum(2*np.arctan2(num,den),axis=1)/(4*np.pi))
 return np.array(out)
rows=[];details=[]
for file in sorted(O.glob('*.npz')):
 if file.name=='digital_coupon_mesh.npz':continue
 z=np.load(file);t=time.perf_counter();w=winding(z['vertices'],z['faces'],z['queries']);inside=np.abs(w)>.5;parity=z['ray_votes']>=2;expected=z['expected'];known=expected>=0
 row={'case_id':file.stem,'queries':len(w),'parity_winding_disagreements':int(np.count_nonzero(inside!=parity)),'nearest_face_winding_disagreements':int(np.count_nonzero(inside!=z['nearest_face_inside'])),'known_expected_queries':int(known.sum()),'winding_expected_errors':int(np.count_nonzero(inside[known]!=expected[known])),'parity_expected_errors':int(np.count_nonzero(parity[known]!=expected[known])),'numpy_seconds':time.perf_counter()-t,'maximum_distance_to_integer_winding':float(np.max(np.abs(np.abs(w)-np.round(np.abs(w)))))}
 rows.append(row)
 for idx in range(len(w)):details.append({'case_id':file.stem,'index':idx,'kind':str(z['kinds'][idx]),'winding':float(w[idx]),'inside_winding':bool(inside[idx]),'inside_parity':bool(parity[idx]),'inside_nearest_face':bool(z['nearest_face_inside'][idx]),'expected':int(expected[idx])})
 print(json.dumps(row),flush=True)
(O/'independent_winding_summary.json').write_text(json.dumps(rows,indent=2)+'\n');(O/'independent_winding_queries.json').write_text(json.dumps(details,indent=2)+'\n')

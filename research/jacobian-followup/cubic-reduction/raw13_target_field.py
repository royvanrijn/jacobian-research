#!/usr/bin/env python3
"""Select an exact specialization rank certificate for the full cubic target field."""
from pathlib import Path
import json,hashlib,time
B=Path(__file__).resolve().parent;O=B/'raw13-target-field'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
P=1000000007
BASIS=[(i,j,d-i-j) for d in range(8) for i in range(d,-1,-1) for j in range(d-i,-1,-1)]
def add(a,b):return tuple((x+y)%P for x,y in zip(a,b))
def scale(a,c):return tuple(x*c%P for x in a)
def mul(a,b,st,q):
 z=[0]*5
 for i,x in enumerate(a):
  for j,y in enumerate(b):z[i+j]=(z[i+j]+x*y)%P
 for k in [4,3]:
  x=z[k];z[k]=0;z[k-1]=(z[k-1]+x)%P;z[k-2]=(z[k-2]-st*x)%P;z[k-3]=(z[k-3]+q*x)%P
 return tuple(z[:3])
def solve(M,b):
 A=[list(row)+[x] for row,x in zip(M,b)]
 for i in range(3):
  pivot=next((j for j in range(i,3) if A[j][i]),None)
  if pivot is None:return None
  A[i],A[pivot]=A[pivot],A[i];v=pow(A[i][i],P-2,P);A[i]=[x*v%P for x in A[i]]
  for j in range(3):
   if j!=i:
    x=A[j][i];A[j]=[(a-x*b)%P for a,b in zip(A[j],A[i])]
 return tuple(row[-1] for row in A)
def point_rows(point):
 s,t,u=point;st=s*t%P;q=2*s*s*u%P;one=(1,0,0);w=(0,1,0);e=(st,-2,3)
 columns=[mul(e,z,st,q) for z in [one,w,(0,0,1)]];einverse=solve(list(map(list,zip(*columns))),one)
 if einverse is None or s%P==0:return None
 invs=pow(s,P-2,P);inv2=pow(2,P-2,P);a=scale(einverse,s);b=scale(add(w,scale(e,-1)),inv2*invs%P);e2=mul(e,e,st,q);e3=mul(e2,e,st,q);c=scale(add(add(scale(e2,5),scale(mul(w,e,st,q),-3)),scale(e3,-2)),inv2*invs*invs%P)
 powers=[]
 for x in [a,b,c]:
  z=[one]
  for k in range(7):z.append(mul(z[-1],x,st,q))
  powers.append(z)
 vals=[mul(mul(powers[0][i],powers[1][j],st,q),powers[2][k],st,q) for i,j,k in BASIS]
 return [[z[j] for z in vals] for j in [1,2]]
def main():
 start=time.monotonic();O.mkdir(exist_ok=True);raw=B/'raw13-specialization/input.json';D=json.loads(raw.read_text());positions={e:i for i,e in enumerate(BASIS)};known=[]
 for row in range(3):
  v=[0]*len(BASIS);ee=[0]*3;ee[row]=1;v[positions[tuple(ee)]]=1
  for atom in D['atoms']:
   if atom['row']==row:v[positions[tuple(atom['exponent'])]]=int(atom['coefficient'])%P
  known.append(v)
 pivot_basis={};selected=[];points=[];all_rows=[];seed=20261115;skipped=[];cp=O/'checkpoint.json'
 def checkpoint(status):
  z={'status':status,'source_sha256':sha(Path(__file__)),'source_input_sha256':sha(raw),'modulus':P,'source_degree_cap':7,'columns':len(BASIS),'rank':len(pivot_basis),'points':len(points),'selected_rows':len(selected),'primary_method':'Incremental row elimination;certificate checker uses companion matrices and column elimination.','target_degree_or_denominator_cap':None,'candidate_coefficient_search':False,'seconds':time.monotonic()-start}
  cp.write_text(json.dumps(z,indent=2)+'\n')
 checkpoint('running')
 for attempt in range(240):
  pt=[]
  for j in range(3):seed=(1664525*seed+1013904223)%2**32;pt.append(seed%97+1)
  if pt in points:continue
  rows=point_rows(pt)
  if rows is None:skipped.append(pt);continue
  pi=len(points);points.append(pt)
  for component,row in enumerate(rows,1):
   assert row[0]==0 and all(sum(x*y for x,y in zip(row,v))%P==0 for v in known)
   index=len(all_rows);all_rows.append(row);red=list(row)
   for pivot,old in sorted(pivot_basis.items()):
    fac=red[pivot]
    if fac:red=[(a-fac*b)%P for a,b in zip(red,old)]
   pivot=next((i for i,x in enumerate(red) if x),None)
   if pivot is not None:
    inv=pow(red[pivot],P-2,P);pivot_basis[pivot]=[x*inv%P for x in red];selected.append({'point':pi,'component':component,'row_index':index,'pivot_column':pivot})
  assert len(pivot_basis)<=116
  if len(pivot_basis)==116:break
  if len(points)%16==0:checkpoint('running');print('RANK',len(points),len(pivot_basis),flush=True)
 if len(pivot_basis)!=116:
  checkpoint('terminal-incomplete');print('INCOMPLETE rank',len(pivot_basis),flush=True);return
 payload={'source_sha256':sha(Path(__file__)),'source_input_sha256':sha(raw),'modulus':P,'source_basis':[list(e) for e in BASIS],'target_points':points,'excluded_singular_points':skipped,'selected_rows':selected,'minor_columns':[z['pivot_column'] for z in selected],'expected_rank':116,'kernel_dimension':4,'cubic':'W^3-W^2+s*t*W-2*s^2*u','derivative':'E=3*W^2-2*W+s*t','source_recovery':['a=s/E','b=(W-E)/(2*s)','c=(5*E^2-3*W*E-2*E^3)/(2*s^2)'],'scope':'Necessary full target-field membership functionals on source degree<=7;rank uses exact rational specialization,not candidate-map sampling.'}
 ip=O/'input.json';ip.write_text(json.dumps(payload,indent=2)+'\n');checkpoint('complete');z=json.loads(cp.read_text());z['input_sha256']=sha(ip);cp.write_text(json.dumps(z,indent=2)+'\n');print('COMPLETE rank116/120,',len(points),'target evaluation points;full field membership,no target-degree cap',flush=True)
if __name__=='__main__':main()

#!/usr/bin/env python3
"""Select an exact specialization rank certificate for the full cubic target field."""
from pathlib import Path
import json,hashlib,time
B=Path(__file__).resolve().parent;O=B/'negative-source-quadratic/base-field14'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
P=1000000007
BASIS=[(i,j,d-i-j) for d in range(15) for i in range(d,-1,-1) for j in range(d-i,-1,-1)]
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
  for k in range(14):z.append(mul(z[-1],x,st,q))
  powers.append(z)
 vals=[mul(mul(powers[0][i],powers[1][j],st,q),powers[2][k],st,q) for i,j,k in BASIS]
 return [[z[j] for z in vals] for j in [1,2]]
def main():
 start=time.monotonic();O.mkdir(parents=True,exist_ok=True);raw=B/'raw13-specialization/input.json';D=json.loads(raw.read_text());positions={e:i for i,e in enumerate(BASIS)};known=[]
 from verify import variable,add as pa,mul as pm,scale as ps,compose
 from fractions import Fraction
 from itertools import product
 F=[variable(i,3) for i in range(3)]
 for atom in D['atoms']:F[atom['row']]=pa(F[atom['row']],{tuple(atom['exponent']):Fraction(atom['coefficient'])})
 def power(f,k):
  v={(0,0,0):Fraction(1)}
  for _ in range(k):v=pm(v,f)
  return v
 S=pa(pa(power(F[1],3),ps(pm(F[0],power(F[2],2)),27)),ps(pm(F[1],F[2]),-9));T=pa(power(F[1],2),ps(pm(F[0],S),-4));gens=F+[S,T];degrees=[4,6,7,8,8];exact=[];targets=[];pivots={}
 for ex in sorted(product(range(4),range(3),range(3),range(2),range(2)),key=lambda e:(sum(d*k for d,k in zip(degrees,e)),e)):
  if sum(d*k for d,k in zip(degrees,ex))>14:continue
  v={(0,0,0):Fraction(1)}
  for f,k in zip(gens,ex):v=pm(v,power(f,k))
  red=dict(v)
  for pivot,old in sorted(pivots.items()):
   fac=red.get(pivot,0)
   if fac:red=pa(red,ps(old,-fac))
  if not red:continue
  pivot=min(red);pivots[pivot]=ps(red,1/red[pivot]);exact.append(v);targets.append(list(ex))
 for v in exact:known.append([int(v.get(e,0))%P for e in BASIS])
 from verify import encode
 (O/'known.json').write_text(json.dumps({'generators':['H1','H2','H3','S','T'],'target_generator_exponents':targets,'source_polynomials':[encode(v) for v in exact],'dimension':len(exact)},indent=2)+'\n')
 expected=len(BASIS)-len(exact)
 print('KNOWN',len(exact),'expectedrank',expected,flush=True)
 pivot_basis={};selected=[];points=[];all_rows=[];seed=20261115;skipped=[];cp=O/'checkpoint.json'
 def checkpoint(status):
  z={'status':status,'source_sha256':sha(Path(__file__)),'source_input_sha256':sha(raw),'modulus':P,'source_degree_cap':14,'columns':len(BASIS),'rank':len(pivot_basis),'points':len(points),'selected_rows':len(selected),'primary_method':'Incremental row elimination;certificate checker uses companion matrices and column elimination.','target_degree_or_denominator_cap':None,'candidate_coefficient_search':False,'seconds':time.monotonic()-start}
  cp.write_text(json.dumps(z,indent=2)+'\n')
 checkpoint('running')
 for attempt in range(400):
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
  assert len(pivot_basis)<=expected
  if len(pivot_basis)==expected:break
  if len(points)%16==0:checkpoint('running');print('RANK',len(points),len(pivot_basis),flush=True)
 if len(pivot_basis)!=expected:
  checkpoint('terminal-incomplete');print('INCOMPLETE rank',len(pivot_basis),flush=True);return
 payload={'source_sha256':sha(Path(__file__)),'source_input_sha256':sha(raw),'modulus':P,'source_basis':[list(e) for e in BASIS],'target_points':points,'excluded_singular_points':skipped,'selected_rows':selected,'minor_columns':[z['pivot_column'] for z in selected],'expected_rank':expected,'kernel_dimension':len(exact),'cubic':'W^3-W^2+s*t*W-2*s^2*u','derivative':'E=3*W^2-2*W+s*t','source_recovery':['a=s/E','b=(W-E)/(2*s)','c=(5*E^2-3*W*E-2*E^3)/(2*s^2)'],'scope':'Necessary full target-field membership functionals on source degree<=14;rank uses exact rational specialization,not candidate-map sampling.'}
 ip=O/'input.json';ip.write_text(json.dumps(payload,indent=2)+'\n');checkpoint('complete');z=json.loads(cp.read_text());z['input_sha256']=sha(ip);cp.write_text(json.dumps(z,indent=2)+'\n');print('COMPLETE rank',expected,'/',len(BASIS),,len(points),'target evaluation points;full field membership,no target-degree cap',flush=True)
if __name__=='__main__':main()

#!/usr/bin/env python3
"""Exact finite E6/E7 Cartan root data ONLY; not an all-order geometric checker."""
from fractions import Fraction as F
import argparse
from pathlib import Path
import json,hashlib
def inverse(A):
 n=len(A);M=[[F(x) for x in A[i]]+[F(i==j) for j in range(n)] for i in range(n)]
 for i in range(n):
  r=next(j for j in range(i,n) if M[j][i]);M[i],M[r]=M[r],M[i]
  c=M[i][i];M[i]=[x/c for x in M[i]]
  for j in range(n):
   if j!=i:
    c=M[j][i];M[j]=[a-c*b for a,b in zip(M[j],M[i])]
 return [r[n:] for r in M]
def model(n,edges,omega):
 A=[[2*int(i==j) for j in range(n)] for i in range(n)]
 for i,j in edges:A[i][j]=A[j][i]=-1
 Ai=inverse(A);start=tuple(int(i==omega) for i in range(n));orb={start};todo=[start]
 while todo:
  v=todo.pop()
  for i in range(n):
   w=tuple(v[j]-v[i]*A[j][i] for j in range(n))
   if w not in orb:orb.add(w);todo.append(w)
 roots=set();todo=[]
 for i in range(n):
  v=tuple(int(j==i) for j in range(n));roots.add(v);todo.append(v)
 while todo:
  v=todo.pop()
  for i in range(n):
   c=sum(A[i][j]*v[j] for j in range(n));w=tuple(v[j]-(c if j==i else 0) for j in range(n))
   if w not in roots:roots.add(w);todo.append(w)
 pos=[v for v in roots if min(v)>=0]
 norm=Ai[omega][omega];dots=sorted({sum(Ai[omega][j]*v[j] for j in range(n)) for v in orb})
 return {'rank':n,'edges':edges,'minuscule_node':omega,'cartan':A,'inverse_cartan':[[str(x) for x in r] for r in Ai],'weight_orbit_size':len(orb),'roots':len(roots),'positive_roots':len(pos),'positive_roots_outside_centralizer':sum(v[omega]!=0 for v in pos),'norm':str(norm),'dot_values':[str(x) for x in dots],'weight_orbit':sorted(orb),'root_set':sorted(roots),'positive_root_set':sorted(pos),'centralizer_root_set':sorted(v for v in roots if v[omega]==0),'root_reflection_closed':all(tuple(v[j]-(sum(A[i][q]*v[q] for q in range(n)) if j==i else 0) for j in range(n)) in roots for v in roots for i in range(n)),'weight_reflection_closed':all(tuple(v[j]-v[i]*A[j][i] for j in range(n)) in orb for v in orb for i in range(n)),'all_roots_squared_length_two':all(sum(v[i]*A[i][j]*v[j] for i in range(n) for j in range(n))==2 for v in roots)}
d={'EIII':model(6,[(0,1),(1,2),(2,3),(3,4),(2,5)],0),'EVII':model(7,[(0,1),(1,2),(2,3),(3,4),(4,5),(2,6)],5)}
assert (d['EIII']['weight_orbit_size'],d['EIII']['roots'],d['EIII']['positive_roots_outside_centralizer'])==(27,72,16)
assert d['EIII']['dot_values']==['-2/3','1/3','4/3']
assert (d['EVII']['weight_orbit_size'],d['EVII']['roots'],d['EVII']['positive_roots_outside_centralizer'])==(56,126,27)
assert d['EVII']['dot_values']==['-3/2','-1/2','1/2','3/2']
d['status']='PASS_EXACT_FINITE_ROOT_DATA_ONLY'
d['assurance_flags']={k:False for k in ['external_review','formal_verification','independent_replay','worldwide_priority']}
ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',default=str(Path(__file__).with_name('piecewise_coadjoint_exceptional_root_data.json')));args=ap.parse_args()
assert all(d[k][f] for k in ['EIII','EVII']for f in ['root_reflection_closed','weight_reflection_closed','all_roots_squared_length_two'])
out=Path(args.output);out.write_text(json.dumps(d,indent=2)+'\n')
print(str(out),hashlib.sha256(out.read_bytes()).hexdigest())
print({k:{q:v[q] for q in ['weight_orbit_size','roots','positive_roots_outside_centralizer','norm','dot_values']} for k,v in d.items() if k in ['EIII','EVII']})



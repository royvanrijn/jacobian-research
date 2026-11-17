#!/usr/bin/env python3
"""Independent full-map and branch-certificate audit of the19-atom weighted tail family."""
from pathlib import Path
from fractions import Fraction as Q
from math import comb
from collections import Counter
import hashlib,json,time
from verify import add,mul,variable,derivative,determinant3,evaluate
from verify_power_single_compress import decode,scale
from weighted_verify_utils import power,substitution,expression,equal_generators,rank,coefficient_buckets
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def canonical(p):return tuple(sorted(p.items()))
def main():
 started=time.monotonic();source=HERE/'weighted-tail7/input.json';d=json.loads(source.read_text());names=d['variables'];n=len(names);assert n==19;index={name:i for i,name in enumerate(names)};atoms=d['tail_atoms'];assert len(atoms)==17
 assert d['source_sha256']==sha(HERE/'weighted_tail7.py')
 expected=set()
 for row in [1,2]:
  for a in range(8):
   for b in range(8):
    for c in range(8):
     if 2<=a+b+c<=7 and a-b-2*c==-row:expected.add((row,(a,b,c)))
 assert {(a['row'],tuple(a['exponent'])) for a in atoms}==expected and len(expected)==17
 # Direct universal22-variable determinant:3 source coordinates and19 parameters.
 N=3+n;F=[variable(i,N) for i in range(3)];fixed=[((2,1,0),Q(-3)),((3,0,1),Q(-1))]
 for e,c in fixed:F[0]=add(F[0],{e+(0,)*n:c})
 for atom in atoms:
  e=list(atom['exponent'])+[0]*n;e[3+index[atom['coefficient']]]=1;F[atom['row']]=add(F[atom['row']],{tuple(e):Q(1)})
 determinant=add(determinant3([[derivative(f,j) for j in range(3)] for f in F],N),{(0,)*N:Q(1)},Q(-1));buckets=coefficient_buckets(determinant,3)
 assert all(e[0]==e[1]+2*e[2] for e in buckets)
 jac=list(buckets.values());one={(0,)*n:Q(1)};beta=variable(index['beta'],n);sat=variable(index['sat'],n)
 boundaries=[{1:dict(one)},{0:dict(one),1:scale(one,Q(-3))}]
 for atom in atoms:
  row=atom['row']-1;b,c=atom['quotient_exponent'];assert [b,c]==atom['exponent'][1:]
  v=variable(index[atom['coefficient']],n)
  for j in range(c+1):boundaries[row][b+j]=add(boundaries[row].get(b+j,{}),v,Q(comb(c,j)*((-3)**j)))
 boundary=[f for k,f in boundaries[0].items() if k>=2 and f]+[f for k,f in boundaries[1].items() if k>=3 and f]
 assert Counter(map(canonical,boundary))==Counter(canonical(decode(f,n)) for f in d['boundary_equations'])
 collision={}
 for atom in atoms:
  if atom['row']!=1:continue
  b,c=atom['quotient_exponent'];i=index[atom['coefficient']]
  for j in range(c+1):
   e=[0]*n;e[i]=1;e[index['beta']]+=b+j;collision=add(collision,{tuple(e):Q(((-1)**b)*comb(c,j)*(3**j))})
 saturation=add(one,mul(beta,sat),Q(-1));core=jac+[collision,saturation];actual=core+boundary;original=[decode(f,n) for f in d['equations']]
 assert Counter(map(canonical,actual))==Counter(map(canonical,original))
 lookup={(a['row'],tuple(a['quotient_exponent'])):index[a['coefficient']] for a in atoms}
 branch_zeros={
  'small-h':[index[a['coefficient']] for a in atoms if a['row']==2 and a['quotient_exponent'] not in [[2,0],[1,1]]],
  'high-w':[lookup[key] for key in [(1,(1,2)),(1,(2,1)),(2,(4,0)),(2,(3,1))]],
  'mid-w':[lookup[key] for key in [(1,(1,2)),(2,(0,3)),(2,(1,2)),(2,(4,0)),(2,(3,1))]]}
 out=HERE/'weighted-unit-certificates';manifest=json.loads((out/'manifest.json').read_text());assert manifest['status']=='three-complete-branch-unit-certificates';caseinfo={c['tag']:c for c in manifest['cases']};results=[];maximum=0;terms=0
 paths={'small-h':HERE/'weighted-tail7/small-h/input.json','high-w':HERE/'weighted-tail7/high-w/generic/input.json','mid-w':HERE/'weighted-tail7/mid-w/input.json'}
 for tag in ['small-h','high-w','mid-w']:
  data=json.loads(paths[tag].read_text());free=data['variables'];m=len(free);values=[expression(data['substitution'][name],free) for name in names]
  linear=boundary+[variable(i,n) for i in branch_zeros[tag]]
  assert all(sum(e)==1 for f in linear for e in f)
  oldfree=[name for name in free if name in index]
  for i,name in enumerate(names):
   if name in free:assert values[i]==variable(free.index(name),m)
  assert all(not substitution(f,values,m) for f in linear)
  assert rank([[f.get(tuple(int(k==j) for k in range(n)),Q(0)) for j in range(n)] for f in linear],n)==n-len(oldfree)
  derived=[substitution(f,values,m) for f in original];derived=[f for f in derived if f]
  if tag!='small-h':
   inv=variable(free.index('extra_inv'),m);K=values[lookup[1,(0,2)]];product=mul(inv,K)
   if tag=='high-w':product=mul(product,values[lookup[2,(0,3)]])
   derived.append(add({(0,)*m:Q(1)},product,Q(-1)))
  reduced=[decode(f,m) for f in data['equations']];equal_generators(derived,reduced)
  ip=out/(tag+'-input.json');inputs=json.loads(ip.read_text());assert inputs['input_sha256']==sha(paths[tag]) and inputs['variables']==free
  generators=[decode(f,m) for f in inputs['generators']];assert Counter(map(canonical,generators))==Counter(map(canonical,reduced))
  cp=out/(tag+'.json');cert=json.loads(cp.read_text());assert caseinfo[tag]['certificate_sha256']==sha(cp) and cert['input_sha256']==sha(ip) and cert['source_sha256']==sha(HERE/'certify_weighted_units.py')
  assert cert['singular_source_sha256']==sha(out/(tag+'.sing'));multipliers=[decode(f,m) for f in cert['multipliers']];assert len(multipliers)==len(generators)
  result={}
  for a,b in zip(generators,multipliers):result=add(result,mul(a,b))
  assert result=={(0,)*m:Q(1)},tag;maximum=max(maximum,max((sum(e) for f in multipliers for e in f),default=0));terms+=sum(len(f) for f in multipliers);results.append(dict(tag=tag,certificate_sha256=sha(cp),input_sha256=sha(paths[tag])))
 # Independently check the top-v coefficient identities used in the structural proof.
 k=19;u=variable(0,k);v=variable(1,k);z=[variable(i+2,k) for i in range(17)];oneK={(0,)*k:Q(1)}
 def poly(coefs):
  p={}
  for i,c in enumerate(coefs):p=add(p,mul(c,power(u,i,k)))
  return p
 Gs=[poly(z[0:2]),poly(z[2:5]),poly(z[5:7])];Hs=[poly(z[7:10]),poly(z[10:14]),poly(z[14:16]),z[16]]
 G={};H={}
 for i,f in enumerate(Gs):G=add(G,mul(f,power(v,i,k)))
 for i,f in enumerate(Hs):H=add(H,mul(f,power(v,i,k)))
 PDE=add(mul(add(G,mul(v,derivative(G,1))),derivative(H,0)),mul(derivative(G,0),add(scale(H,Q(2)),mul(v,derivative(H,1)))),Q(-1))
 def cv(p,a):
  out={}
  for e,c in p.items():
   if e[1]==a:
    ee=list(e);ee[1]=0;out[tuple(ee)]=c
  return out
 assert cv(PDE,5)==scale(mul(derivative(Gs[2],0),Hs[3]),Q(-5))
 c4=add(add(scale(mul(Gs[2],derivative(Hs[2],0)),Q(3)),mul(derivative(Gs[2],0),Hs[2]),Q(-4)),mul(derivative(Gs[1],0),Hs[3]),Q(-5));assert cv(PDE,4)==c4
 c3=add(add(scale(mul(Gs[2],derivative(Hs[1],0)),Q(3)),mul(derivative(Gs[2],0),Hs[1]),Q(-3)),scale(mul(Gs[1],derivative(Hs[2],0)),Q(2)));c3=add(add(c3,mul(derivative(Gs[1],0),Hs[2]),Q(-4)),mul(derivative(Gs[0],0),Hs[3]),Q(-5));assert cv(PDE,3)==c3
 # Linear-g lemma coordinate change: V(q)=kappa*q, V(v)=-kappa*v.
 t=variable(2,3);uu=variable(0,3);vv=variable(1,3);one3={(0,)*3:Q(1)};kap=add(one3,t,Q(-3));q=add(mul(kap,uu),mul(t,add(one3,vv,Q(-1))));A=add(mul(kap,uu),mul(t,add(one3,vv,Q(-2))))
 vector=lambda f:add(mul(A,derivative(f,0)),mul(mul(kap,vv),derivative(f,1)),Q(-1))
 assert vector(q)==mul(kap,q) and vector(vv)==scale(mul(kap,vv),Q(-1))
 # Reduction to the already independently certified normalized raw13 equations.
 oldraw=json.loads((HERE/'raw13-specialization/input.json').read_text());oldnormal=json.loads((HERE/'raw13-normalized/input.json').read_text());oldv=json.loads((HERE/'raw13-rigidity/validation.json').read_text())
 assert oldv['status']=='pass' and oldv['source_sha256']==sha(HERE/'verify_raw13_rigidity.py') and oldv['normalized_input_sha256']==sha(HERE/'raw13-normalized/input.json')
 oldmap={(a['row'],tuple(a['exponent'])):j for j,a in enumerate(oldraw['atoms'][2:])};oldvalues=[{} for i in range(n)];known=[Q(0)]*n
 for atom in atoms:
  i=index[atom['coefficient']];key=(atom['row'],tuple(atom['exponent']))
  if key in oldmap:oldvalues[i]=variable(oldmap[key],13);known[i]=Q(oldraw['atoms'][2+oldmap[key]]['coefficient'])
 oldvalues[index['beta']]=variable(11,13);oldvalues[index['sat']]=variable(12,13);known[index['beta']]=known[index['sat']]=1
 reducedcore=[substitution(f,oldvalues,13) for f in core];equal_generators(reducedcore,[decode(f,13) for f in oldnormal['equations']])
 assert all(evaluate(f,known)==0 for f in original) and 2+sum(v!=0 for v in known[:17])==13
 receipt=dict(status='pass',source_sha256=sha(Path(__file__)),arithmetic_sha256=sha(HERE/'verify.py'),utility_sha256=sha(HERE/'weighted_verify_utils.py'),input_sha256=sha(source),raw13_validation_sha256=sha(HERE/'raw13-rigidity/validation.json'),
  nonlinear_atoms=19,tail_atoms=[8,9],branch_unit_identities=3,max_multiplier_degree=maximum,multiplier_monomials=terms,cases=results,
  full_three_variable_determinant_reconstructed=True,necessary_boundary_equations_reconstructed=True,linear_parameterizations_complete=True,all_branch_identities_multiplied=True,top_v_identities_checked=True,linear_g_coordinate_change_checked=True,old13_reduction_checked=True,
  unique_nonlinear_total=13,seconds=round(time.monotonic()-started,3),scope='With first output a-3a^2b-a^3c, every degree<=7 weight-compatible tail pair of source weights(1,-1,-2), identity linear part and specified beta!=0 collision is the known representative. All six added tail atoms vanish. This does not cover changing the first output, arbitrary weights, collision families or global raw support minima.')
 (HERE/'weighted-tail7/validation.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n');print(json.dumps(receipt),flush=True)
if __name__=='__main__':main()

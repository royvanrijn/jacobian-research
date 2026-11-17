#!/usr/bin/env python3
"""Independent live-substitution replay of all512 linear carrier-subspace charts."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import combinations
from collections import defaultdict,Counter
import hashlib,json
from verify import decode,variable,compose,add,scale,mul,evaluate
HERE=Path(__file__).resolve().parent
DEGREES={3:6,4:6,5:5,6:4,7:4,8:3,9:3,10:3,11:3,12:3}
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def param_decode(row,np):
 p={tuple(t['exponent']):Q(t['coefficient']) for t in row};assert len(p)==len(row) and all(len(e)==np and all(k>=0 for k in e) and c for e,c in p.items());return p
def source_evaluate(p,point,np):
 out={}
 for e,c in p.items():
  for a,k in zip(point,e[:12]):c*=a**k
  v=e[12:];out[v]=out.get(v,Q(0))+c
 return {e:c for e,c in out.items() if c}
def main():
 rawpath=HERE/'cubic12-integer-witness.json';raw=json.loads(rawpath.read_text());F=[decode(p,12) for p in raw['coefficients']];points=[[Q(v) for v in p] for p in raw['points']];image=list(map(Q,raw['common_image']))
 assert points[0]!=points[1] and all([evaluate(p,pt) for p in F]==image for pt in points)
 cp=HERE/'linear-carrier-subspaces12/certificates.json';cert=json.loads(cp.read_text());assert cert['status']=='complete-primary-certificates' and cert['charts']==512 and cert['unit_obstructions']==508
 assert cert['source_sha256']==sha(HERE/'linear_carrier_certificates12.py') and cert['model_producer_sha256']==sha(HERE/'linear_carrier_subspaces12.py')
 records={(r['dimension'],r['case']):r for r in cert['cases']};assert len(records)==512;units=0;linear_units=0;fixed_bounds=[];families={};old_witnesses={};model_hashes={};total_coefficients=0
 for n,D in DEGREES.items():
  choices=list(combinations(range(3,12),n-3));folder=HERE/('linear-subspaces-n%02d-d%d'%(n,D));checkpoint=json.loads((folder/'checkpoint.json').read_text())
  assert checkpoint['status']=='complete' and checkpoint['charts_total']==len(choices)==checkpoint['charts_completed'] and checkpoint['source_sha256']==sha(HERE/'linear_carrier_subspaces12.py') and checkpoint['input_sha256']==sha(rawpath)
  bestpath=HERE/('carrier-subset-frontier12/best-n%02d.json'%n);best=json.loads(bestpath.read_text());old_witnesses[str(n)]=sha(bestpath)
  for case,kept in enumerate(choices):
   deleted=[j for j in range(3,12) if j not in kept];pairs=[(j,k) for j in deleted for k in kept if k<j];np=len(pairs);nv=12+np;X=[variable(j,nv) for j in range(nv)];embed=lambda p:{e+(0,)*np:c for e,c in p.items()};current=[embed(p) for p in F];E=X[:12]
   # Reconstruct the LIVE graph equations, instead of composing a precomputed embedding.
   for j in deleted:
    h=current[j]
    for z,(a,k) in enumerate(pairs):
     if a==j:h=add(h,mul(X[12+z],current[k]))
    graph=add(X[j],h,Q(-1));assert all(e[j]==0 and not any(e[k] for k in range(j+1,12)) for e in graph)
    S=list(X);S[j]=graph;current=[compose(p,S,nv) for p in current];E=[compose(p,S,nv) for p in E]
   retained=[0,1,2]+list(kept);G={i:current[i] for i in retained}
   for j in deleted:
    h=current[j]
    for z,(a,k) in enumerate(pairs):
     if a==j:h=add(h,mul(X[12+z],current[k]))
    assert not h
   assert all(not any(e[j] for j in deleted) for p in G.values() for e in p)
   assert all({e:c for e,c in p.items() if sum(e[:12])<=1}==X[i] for i,p in G.items())
   for pt in points:
    for j,p in enumerate(E):assert source_evaluate(p,pt,np)==({(0,)*np:pt[j]} if pt[j] else {})
    for i,p in G.items():assert source_evaluate(p,pt,np)==({(0,)*np:image[i]} if image[i] else {})
   assert [points[0][j] for j in retained]!=[points[1][j] for j in retained]
   mp=folder/('case-%03d.json'%case);model=json.loads(mp.read_text());rec=records[(n,case)]
   assert rec['model']==str(mp.relative_to(HERE)) and rec['model_sha256']==sha(mp);model_hashes[str(mp.relative_to(HERE))]=sha(mp)
   assert model['input_sha256']==sha(rawpath) and model['source_sha256']==sha(HERE/'linear_carrier_subspaces12.py')
   assert model['parameter_pairs']==[[j+1,k+1] for j,k in pairs] and model['kept']==[j+1 for j in kept] and model['deleted']==[j+1 for j in deleted]
   assert [decode(p,nv) for p in model['source_embedding']]==E and {int(i)-1:decode(p,nv) for i,p in model['outputs'].items()}==G
   slots={}
   for i,p in G.items():
    groups=defaultdict(dict)
    for e,c in p.items():
     if sum(e[:12])>=2:groups[e[:12]][e[12:]]=c
    slots.update({(i,e):q for e,q in groups.items()})
   high={k:q for k,q in slots.items() if sum(k[1])>D};listed=model['high_coefficients'];assert [(int(t['output'])-1,tuple(t['source_exponent'])) for t in listed]==sorted(high)
   high_polys=[param_decode(t['parameter_polynomial'],np) for t in listed];assert high_polys==[high[k] for k in sorted(high)];total_coefficients+=len(high)
   def combination(items):
    out={};assert len({t['high_index'] for t in items})==len(items)
    for t in items:assert Q(t['coefficient']);out=add(out,high_polys[t['high_index']],Q(t['coefficient']))
    return out
   if rec['status']=='unit-obstruction':
    assert combination(rec['unit'])=={(0,)*np:Q(1)};units+=1
    if model['status']!='constant-obstruction':linear_units+=1
   else:
    assert rec['status']=='surviving-family';families[n]=(case,G,slots,np,X,rec)
    for t in rec['generators']:assert combination(t['linear_combination'])==param_decode(t['polynomial'],np)
   fixed=[k for k,q in slots.items() if set(q)=={(0,)*np}];fixed_bounds.append((n,case,len(fixed)))
   if len(fixed)<13:
    assert n==6 and case==22 and len(fixed)==12 and rec['status']=='unit-obstruction'
    assert all(set(high_polys[t['high_index']])!={(0,)*np} for t in rec['unit'])
   if model['deleted']==best['deleted']:
    compressed=[{tuple(e[j] for j in retained):c for e,c in p.items() if not any(e[12:])} for p in G.values()]
    assert compressed==[decode(p,n) for p in best['coefficients']]
    assert max(sum(e) for p in compressed for e in p)==best['degree']
    assert sum(sum(e)>=2 for p in compressed for e in p)==best['nonlinear_total']
  print('DIMENSION',n,'PASS',len(choices),'charts',flush=True)
 assert units==508 and linear_units==9 and set(families)=={4,5,7,12}
 # Read and pin the previously independently verified positive frontier.
 oldpath=HERE/'carrier-subset-frontier12/validation.json';old=json.loads(oldpath.read_text());assert old['status']=='pass' and old['all_subsets_replayed']==512 and old['source_sha256']==sha(HERE/'verify_carrier_subsets12.py') and old['input_sha256']==sha(rawpath)
 for n,digest in old['witness_hashes'].items():assert old_witnesses[n]==digest
 family_records=[]
 for n,(case,G,slots,np,X,rec) in families.items():
  T=X[12:];one={(0,)*np:Q(1)};pt=[variable(j,np) for j in range(np)];pair=lambda a,b:mul(pt[a],pt[b]);expected=[]
  base={i:{e:c for e,c in p.items() if not any(e[12:])} for i,p in G.items()}
  if n==4:
   expected=[pt[2]];S=X[:12]+T[:];S[12+2]={};actual={i:compose(p,S,12+np) for i,p in G.items()};U=G[4];UU=mul(U,U)
   wanted=dict(base)
   for i,c,a,b in [(0,1,0,1),(1,-12,1,6),(2,-4,1,5)]:wanted[i]=add(wanted[i],scale(mul(mul(T[a],T[b]),UU),c))
   assert actual==wanted and len(UU)==3
   assert all(not(set(base[i])&set(UU)) for i in [0,1,2])
   e=[0]*12;e[0]=3;e[1]=2;e[2]=1;assert slots[(1,tuple(e))]=={(0,)*np:Q(12)}
   support_formula='15+3[t0*t1!=0]+3[t1*t6!=0]+3[t1*t5!=0]'
  elif n==5:
   expected=[pair(0,1),pair(1,8),pair(1,6)];U=G[4];V=G[7];UU=mul(U,U);UV=mul(U,V);wanted=dict(base)
   wanted[0]=add(wanted[0],mul(mul(T[0],T[1]),UU))
   wanted[1]=add(add(wanted[1],scale(mul(mul(T[1],T[8]),UU),-12)),scale(mul(mul(T[1],T[9]),UV),-12))
   wanted[2]=add(add(wanted[2],scale(mul(mul(T[1],T[6]),UU),-4)),scale(mul(mul(T[1],T[7]),UV),-4))
   assert G==wanted and len(UV)==4 and all(not(set(base[i])&set(UV)) for i in [1,2])
   e=[0]*12;e[0]=3;e[2]=1;e[7]=1;assert slots[(1,tuple(e))]=={(0,)*np:Q(-12)}
   support_formula='17+4[t1*t9!=0]+4[t1*t7!=0], subject to t0*t1=t1*t8=t1*t6=0'
  elif n==7:
   expected=[pt[4],pt[5]]
   for a,b in [(0,1),(6,7),(10,11)]:expected.extend([pair(2,a),add(pair(2,b),pair(3,a)),pair(3,b)])
   for j in [8,9,12,13]:expected.extend([pair(2,j),pair(3,j)])
   assert len(expected)==19
   for zero_indices in [[2,3,4,5],[0,1,4,5,6,7,8,9,10,11,12,13]]:
    S=list(X)
    for j in zero_indices:S[12+j]={}
    assert {i:compose(p,S,len(X)) for i,p in G.items()}==base
   e=[0]*12;e[0]=3;e[2]=1;assert slots[(0,tuple(e))]=={(0,)*np:Q(-1)}
   support_formula='23 on both complete branches'
  else:assert n==12 and np==0 and G=={i:p for i,p in enumerate(F)};support_formula='36'
  assert expected==[param_decode(t['polynomial'],np) for t in rec['generators']]
  counts=sum(sum(e[:12])>=2 for p in base.values() for e in p);assert counts=={4:15,5:17,7:23,12:36}[n]
  family_records.append({'dimension':n,'case':case,'necessary_generators':len(expected),'normal_form_checked':True,'minimum_nonlinear_occurrences_at_degree_cap':counts,'support_formula':support_formula})
 receipt={'status':'pass','claim':'KCUB29','checker_sha256':sha(Path(__file__)),'arithmetic_sha256':sha(HERE/'verify.py'),'input_sha256':sha(rawpath),'producer_sha256':sha(HERE/'linear_carrier_subspaces12.py'),'certificate_producer_sha256':sha(HERE/'linear_carrier_certificates12.py'),'certificates_sha256':sha(cp),'old_subset_validation_sha256':sha(oldpath),'charts':512,'unit_obstructions':508,'multi_coefficient_constant_units':9,'all_generic_collisions_and_identity_jets':True,'model_hashes':model_hashes,'upper_witness_hashes':old_witnesses,'fixed_nonlinear_coefficient_minima':{str(n):min(k for a,b,k in fixed_bounds if a==n) for n in DEGREES},'one_extra_coefficient_unit':{'dimension':6,'case':22,'fixed':12,'extra_at_least':1},'minimum_occurrences_in_entire_class':13,'minimum_degree_by_dimension':{'3':7,'4':6,'5':5,'6':5,'7':4,'8':4,'9':4,'10':4,'11':4,'12':3},'families':family_records,'scope':'All homogeneous linear carrier-output relation subspaces, normalized highest-pivot graphs with raw source complement and raw retained outputs; no arbitrary affine term bound, nonlinear graphs or target repairs. No global optimum claim.','goal_complete':False}
 proof=HERE/'LINEAR_CARRIER_SUBSPACE_LIMITS.md'
 if proof.exists():receipt['proof_sha256']=sha(proof)
 (HERE/'linear-carrier-subspaces12/validation.json').write_text(json.dumps(receipt,indent=2)+'\n');print('PASS512 charts;508 unit obstructions; raw occurrence minimum13; degree minima verified; goal open',flush=True)
if __name__=='__main__':main()

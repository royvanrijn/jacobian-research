#!/usr/bin/env python3
"""Independent arrangement coverage and exact coefficient-power replay."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import permutations,combinations,product
from math import gcd,lcm
import json,hashlib,time
from verify import decode,variable,add,mul,scale,compose,derivative,determinant3
HERE=Path(__file__).resolve().parent;OUT=HERE/'weight-resonance7';ROOT=HERE.parents[1]
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def primitive(v):
 g=gcd(*v);v=tuple(z//g for z in v)
 if next(z for z in v if z)<0:v=tuple(-z for z in v)
 return v

def canon(v):return min(primitive(p) for p in permutations(v))
def native():
 return [(i,(a,b,d-a-b)) for i in range(3) for d in range(2,8) for a in range(d+1) for b in range(d-a+1)]
def relation(i,e):return tuple(e[j]-int(j==i) for j in range(3))
def proportional(v,w):return all(v[i]*w[j]==v[j]*w[i] for i in range(3) for j in range(i+1,3))
def kernel(v,w):
 A=[list(map(Q,v)),list(map(Q,w))];pivots=[];row=0
 for col in range(3):
  k=next((j for j in range(row,2) if A[j][col]),None)
  if k is None:continue
  A[row],A[k]=A[k],A[row];c=A[row][col];A[row]=[z/c for z in A[row]]
  for j in range(2):
   if j!=row:
    c=A[j][col];A[j]=[x-c*y for x,y in zip(A[j],A[row])]
  pivots.append(col);row+=1
  if row==2:break
 if row<2:return None
 free=next(j for j in range(3) if j not in pivots);z=[Q(0)]*3;z[free]=Q(1)
 for i,p in enumerate(pivots):z[p]=-A[i][free]
 den=lcm(*(x.denominator for x in z));return canon(tuple(int(x*den) for x in z))
def graph_acyclic(edges):return any(all(p.index(j)<p.index(i) for j,i in edges) for p in permutations(range(3)))
def high_rank(support):
 rows=[{e for i,e in support if i==j and sum(e)>=4} for j in range(3)];owner={}
 def augment(i,seen):
  for e in rows[i]:
   if e in seen:continue
   seen.add(e)
   if e not in owner or augment(owner[e],seen):owner[e]=i;return True
  return False
 return sum(augment(i,set()) for i in range(3))
def mixed(w):
 if 0 in w or not min(w)<0<max(w):return None
 if sum(z>0 for z in w)==2:w=tuple(-z for z in w)
 return (next(z for z in w if z>0),)+tuple(sorted(-z for z in w if z<0))
def prior(m):
 if m is None:return None
 p,q,r=m
 if q==r and (p,q)!=(1,1):return 'KCUB32'
 if p>=2 and (p==q or p==r):return 'KCUB33'
 if p==q==1:return {2:'KCUB15:raw13-minimum-only',3:'KCUB18'}.get(r,'KCUB34' if r>=4 else None)
 if (p,q,r)==(1,2,3):return 'KCUB16'
 if (p,q,r)==(2,1,3):return 'KCUB17'
 if p==2 and q==1 and r>=4:return 'KCUB35'
 return None

def coefficient_equations(basis):
 # Direct Leibniz expansion of derivative atoms, with integer source/parameter exponents.
 N=len(basis);zero=(0,)*N;entries=[[[] for _ in range(3)] for _ in range(3)]
 for i in range(3):entries[i][i].append(((0,0,0),zero,1))
 for k,(i,e) in enumerate(basis):
  par=tuple(int(j==k) for j in range(N))
  for j in range(3):
   if e[j]:entries[i][j].append((tuple(e[a]-int(a==j) for a in range(3)),par,e[j]))
 result={}
 for p in permutations(range(3)):
  sign=(-1)**sum(p[i]>p[j] for i in range(3) for j in range(i+1,3))
  for terms in product(*(entries[i][p[i]] for i in range(3))):
   se=tuple(sum(t[0][j] for t in terms) for j in range(3));pe=tuple(sum(t[1][j] for t in terms) for j in range(N));c=sign
   for t in terms:c*=t[2]
   row=result.setdefault(se,{});row[pe]=row.get(pe,0)+c
 result[(0,0,0)][zero]-=1
 return [{e:Q(c) for e,c in row.items() if c} for se,row in sorted(result.items(),reverse=True) if any(row.values())]
def total(*ps):
 out={}
 for p in ps:out=add(out,p)
 return out

def power(p,k,n):
 out={(0,)*n:Q(1)}
 for _ in range(k):out=mul(out,p)
 return out

def con(c,n):return {(0,)*n:Q(c)} if c else {}
def main():
 start=time.monotonic();cp=OUT/'checkpoint.json';meta=json.loads(cp.read_text());src=OUT/'census.json';d=json.loads(src.read_text());assert meta['status']=='complete' and meta['census_sha256']==sha(src) and d['source_sha256']==meta['source_sha256']==sha(HERE/'weight_resonance7.py');atoms=native();relations=sorted({relation(i,e) for i,e in atoms});planes={primitive(v) for v in relations};rays={z for v,w in combinations(relations,2) if (z:=kernel(v,w)) is not None};assert (len(atoms),len(relations),len(planes),len(rays))==(348,182,154,896);assert rays=={tuple(r['weight']) for r in d['records']};assert planes=={tuple(r['relation_line']) for r in d['rank1_planes']};records={};gates={}
 for rec in d['records']:
  w=tuple(rec['weight']);support=sorted((i,e) for i,e in atoms if sum(a*b for a,b in zip(relation(i,e),w))==0);assert support==[(z['output'],tuple(z['exponent'])) for z in rec['support']];assert len(support)==rec['nonlinear_slots'];edges={(j,i) for i,e in support for j in range(3) if i!=j and e[j]};assert sorted(edges)==list(map(tuple,rec['cross_edges']));rank=high_rank(support);assert rank==rec['high_degree_structural_rank'];m=mixed(w);assert rec['mixed_positive_negative_magnitudes']==(list(m) if m else None);gate='acyclic' if graph_acyclic(edges) else 'GLDC10:high-rank<3' if rank<3 else prior(m) or 'unresolved';assert gate==rec['gate'];gates[gate]=gates.get(gate,0)+1;records[w]=rec
 assert gates==d['gate_counts'];assert gates['unresolved']==463
 # Universal rank-one negative relation is one elementary monomial shear.
 rank1_shears=0
 for rec in d['rank1_planes']:
  v=tuple(rec['relation_line']);support=sorted((i,e) for i,e in atoms if proportional(v,relation(i,e)));assert support==[(z['output'],tuple(z['exponent'])) for z in rec['support']];v=tuple(-z for z in v) if sum(v)<0 else v;assert sum(v)>0
  if min(v)<0:
   j=next(i for i,z in enumerate(v) if z<0);assert v[j]==-1 and sum(z<0 for z in v)==1 and len(support)==1;assert support[0]==(j,tuple(v[i]+int(i==j) for i in range(3)));n=4;x=[variable(i,n) for i in range(3)];e=support[0][1];term=mul({e+(0,):Q(1)},variable(3,n));F=list(x);F[j]=add(F[j],term);inv=list(x);inv[j]=add(inv[j],term,Q(-1));assert [compose(z,inv+[variable(3,n)],n) for z in F]==x and [compose(z,F+[variable(3,n)],n) for z in inv]==x and determinant3([[derivative(z,k) for k in range(3)] for z in F],n)==con(1,n);rank1_shears+=1
  else:assert all(z>=0 for z in v)
 assert rank1_shears==99
 # Formal logarithmic Jacobian rank-one update and every possible degree <=6 triple.
 n=10;t=variable(0,n);vs=[variable(i,n) for i in range(1,4)];fs=[variable(i,n) for i in range(4,7)];ds=[variable(i,n) for i in range(7,10)];M=[[add(fs[i] if i==j else {},mul(t,mul(ds[i],vs[j]))) for j in range(3)] for i in range(3)];rhs=mul(fs[0],mul(fs[1],fs[2]))
 for i in range(3):
  z=mul(t,mul(vs[i],ds[i]))
  for j in range(3):
   if i!=j:z=mul(z,fs[j])
  rhs=add(rhs,z)
 assert determinant3(M,n)==rhs;leading_cases=0
 for degrees in product(range(7),repeat=3):
  n=7;T=variable(0,n);VV=[variable(i,n) for i in range(1,4)];AA=[variable(i,n) for i in range(4,7)];FF=[con(1,n) if k==0 else add(con(1,n),mul(AA[i],power(T,k,n))) for i,k in enumerate(degrees)];expr=mul(FF[0],mul(FF[1],FF[2]))
  for i in range(3):
   z=mul(T,mul(VV[i],derivative(FF[i],0)))
   for j in range(3):
    if i!=j:z=mul(z,FF[j])
   expr=add(expr,z)
  lead={e[1:]:c for e,c in expr.items() if e[0]==sum(degrees)};target=con(1,6)
  for i,k in enumerate(degrees):target=add(target,scale(variable(i,6),k))
  for i,k in enumerate(degrees):
   if k:target=mul(target,variable(3+i,6))
  assert lead==target;leading_cases+=1
 # Generic determinants behind both additional elementary zero-weight reductions.
 n=7;za,zb,zc,zd,ze,zf,zg=[variable(i,n) for i in range(n)];assert determinant3([[za,{},{}],[zb,zc,zd],[ze,zf,zg]],n)==mul(za,add(mul(zc,zg),mul(zd,zf),Q(-1)))
 # All mixed native source models, checked by independent derivative-atom expansion.
 lc=OUT/'mixed-lift-checkpoint.json';lifts=json.loads(lc.read_text());assert lifts['status']=='complete' and lifts['source_sha256']==sha(HERE/'certify_weight_resonance7.py') and lifts['census_sha256']==sha(src);expected={tuple(r['mixed_positive_negative_magnitudes']) for r in d['records'] if r['gate']=='unresolved' and r['mixed_positive_negative_magnitudes'] and r['mixed_positive_negative_magnitudes']!=[1,1,1]};assert expected=={tuple(z['magnitudes']) for z in lifts['cases']} and len(expected)==451;counts={};terms=0;targets=0;maxdeg=0;maxpower=0;verified=[];timeouts=[];model_hashes={};cert_hashes={}
 for j,rec in enumerate(lifts['cases']):
  p,q,r=rec['magnitudes'];weight=(p,-q,-r);name=rec['name'];ip=OUT/'mixed-cases'/(name+'.json');inp=json.loads(ip.read_text());assert sha(ip)==rec['input_sha256'] and inp['source_sha256']==sha(HERE/'certify_weight_resonance7.py');basis=sorted((i,e) for i,e in atoms if sum(e[k]*weight[k] for k in range(3))==weight[i]);assert basis==[(z['output'],tuple(z['exponent'])) for z in inp['basis']];N=len(basis);assert inp['parameter_names']==['c'+str(k) for k in range(N)];equations=coefficient_equations(basis);assert equations==[decode(z,N) for z in inp['generators']];first=[k for k,(i,e) in enumerate(basis) if i==0];assert first==inp['positive_output_parameters'];script=OUT/'mixed-cases'/(name+'.sing');assert rec['script_sha256']==sha(script);model_hashes[name]=sha(ip);counts[rec['status']]=counts.get(rec['status'],0)+1
  if rec['status']=='timeout':timeouts.append({'magnitudes':[p,q,r],'name':name,'status':'unchanged-new5-second-timeout'});continue
  assert rec['status']=='complete-coefficient-powers';certp=OUT/'mixed-cases'/(name+'-certificate.json');cert=json.loads(certp.read_text());assert sha(certp)==rec['certificate_sha256'] and cert['source_sha256']==sha(HERE/'certify_weight_resonance7.py') and cert['input_sha256']==sha(ip) and cert['script_sha256']==sha(script);assert cert['unresolved_parameters']==[] and [z['parameter'] for z in cert['targets']]==first
  localterms=0
  for z in cert['targets']:
   k=z['parameter'];ex=z['power'];assert ex in [1,2,4,8,16,32];mm=[decode(z,N) for z in z['multipliers']];assert len(mm)==len(equations);assert total(*(mul(a,b) for a,b in zip(mm,equations)))=={tuple(ex if i==k else 0 for i in range(N)):Q(1)};targets+=1;localterms+=sum(map(len,mm));maxdeg=max(maxdeg,max((sum(e) for p in mm for e in p),default=0));maxpower=max(maxpower,ex)
  assert localterms==rec['multiplier_terms'];terms+=localterms;verified.append([p,q,r]);cert_hashes[name]=sha(certp)
  if len(verified)%50==0:print('CERTIFIED',len(verified),'models',targets,'powers',terms,'multiplier terms',flush=True)
 assert counts=={'complete-coefficient-powers':447,'timeout':4} and (targets,terms)==(594,153878)
 zero_closed={(0,0,1),(0,1,1)};remaining=[]
 for rec in d['records']:
  if rec['gate']!='unresolved':continue
  w=tuple(rec['weight']);m=rec['mixed_positive_negative_magnitudes']
  if w in zero_closed:continue
  if m and m in verified:continue
  remaining.append({'weight':list(w),'mixed_positive_negative_magnitudes':m,'nonlinear_slots':rec['nonlinear_slots'],'source':'new5-second-timeout' if m and m!=[1,1,1] else 'retained-broad-family-unresolved'})
 assert len(remaining)==14
 imports=['GLDC10','KCUB15','KCUB16','KCUB17','KCUB18','KCUB32','KCUB33','KCUB34','KCUB35'];ledger=json.loads((ROOT/'MATH_STATUS.json').read_text());E={e['id']:e for e in ledger['entries']};pins={}
 for k in imports:
  e=E[k];source=ROOT/e['canonical_source'];assert source.is_file();pins[k]={'state':e['state'],'scope_fingerprint':hashlib.sha256((e['state']+'\0'+e['title']+'\0'+e['scope']).encode()).hexdigest()[:16],'canonical_sha256':sha(source)};assert e['state']=='proved'
 proof=HERE/'WEIGHT_RESONANCE_CENSUS7.md';receipt={'status':'pass','claim':'KCUB36','checker_sha256':sha(Path(__file__)),'arithmetic_sha256':sha(HERE/'verify.py'),'census_sha256':sha(src),'census_producer_sha256':sha(HERE/'weight_resonance7.py'),'census_checkpoint_sha256':sha(cp),'coefficient_producer_sha256':sha(HERE/'certify_weight_resonance7.py'),'coefficient_checkpoint_sha256':sha(lc),'labelled_atoms':348,'distinct_relation_vectors':182,'rank1_relation_planes':154,'rank1_generic_shear_controls':rank1_shears,'rank1_leading_cases':leading_cases,'rank2_weight_rays':896,'primary_gate_counts':gates,'complete_native_mixed_models':451,'certified_positive_first_linear_models':len(verified),'coefficient_power_targets':targets,'multiplier_monomials':terms,'max_multiplier_degree':maxdeg,'max_target_power':maxpower,'preserved_timeouts':timeouts,'additional_zero_weight_unit_reductions':2,'remaining_weight_types':remaining,'imports':pins,'model_hashes':model_hashes,'certificate_hashes':cert_hashes,'raw_term_bound_scope':'Every degree<=7 diagonally Gm-equivariant3D in counted source coordinates Keller map with an algebraic-closure collision and<=12 nonlinear occurrences must have one of the14 retained weights, up to coordinate permutation and overall weight sign. General source linear diagonalization is outside the raw-count theorem. Repeated-weight output normalization may change counts, but all newly discarded repeated-weight cases have full no-collision closures;112 is distinct-weight and its sharp13 import is explicit.','no_new_positive_bound_or_full_profile_count':True,'goal_complete':False,'trust_boundary':'Independent exact arrangement/derivative-atom reconstruction and rational certificate multiplication plus written universal rank-one, normalization, module/unit and finite-coverage arguments; imports retain their own written/external boundaries. Not formalized or externally reviewed.','seconds':time.monotonic()-start}
 if proof.exists():receipt['proof_sha256']=sha(proof)
 (OUT/'validation.json').write_text(json.dumps(receipt,indent=2)+'\n');print('PASS universal rank-one lemma and full degree7 weighted frontend:896 rays ->14 retained weights;447 exact full-model closures;594 powers,153878 terms; goal open',flush=True)
if __name__=='__main__':main()

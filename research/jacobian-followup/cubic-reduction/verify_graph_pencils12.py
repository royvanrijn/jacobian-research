#!/usr/bin/env python3
"""Independent uniform dual checks for all weighted carrier graph charts."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import combinations_with_replacement
from collections import defaultdict
import hashlib,json
from verify import decode,variable,mul,add,scale,compose,derivative,evaluate
HERE=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
W=[1,-1,-2,0,-3,2,-1,-2,-2,-1,-1,0]
def trunc(p,D):return {e:c for e,c in p.items() if sum(e[:12])<=D}
def make_graph(F,pivot,others):
 n=12+len(others);X=[variable(i,n) for i in range(n)];embed=lambda p:{e+(0,)*len(others):c for e,c in p.items()}
 h=embed(F[pivot]);graph=add(h,X[pivot],Q(-1))
 for j,k in enumerate(others):graph=add(graph,mul(X[12+j],embed(F[k])));h=add(h,mul(X[12+j],embed(F[k])))
 graph=scale(graph,-1);assert all(e[pivot]==0 for e in graph)
 sub=X[:12];sub[pivot]=graph
 G={i:compose(p,sub,n) for i,p in enumerate(F) if i!=pivot}
 # The full source straightening and target row change have exact inverses.
 forward=X[:];forward[pivot]=add(X[pivot],graph);back=X[:];back[pivot]=add(X[pivot],graph,Q(-1))
 assert compose(forward[pivot],back,n)==X[pivot] and compose(back[pivot],forward,n)==X[pivot]
 assert compose(h,forward,n)==X[pivot]
 target=X[pivot]
 for j,k in enumerate(others):target=add(target,mul(X[12+j],X[k]))
 ti=X[:];ti[pivot]=add(scale(target,-1),X[pivot],2)
 assert compose(target,ti,n)==X[pivot]
 return G,graph,X

def source_evaluate(p,point):
 out={}
 for e,c in p.items():
  v=c
  for a,k in zip(point,e[:12]):v*=a**k
  key=e[12:];out[key]=out.get(key,Q(0))+v
 return {e:c for e,c in out.items() if c}
def collision_check(G,points,image,pivot,np):
 assert tuple(v for i,v in enumerate(points[0]) if i!=pivot)!=tuple(v for i,v in enumerate(points[1]) if i!=pivot)
 for point in points:
  for j,p in G.items():assert source_evaluate(p,point)==({(0,)*np:image[j]} if image[j] else {})

def main():
 path=HERE/'cubic12-integer-witness.json';raw=json.loads(path.read_text());F=[decode(p,12) for p in raw['coefficients']]
 points=[[Q(v) for v in p] for p in raw['points']];image=list(map(Q,raw['common_image']))
 assert points[0]!=points[1] and all([evaluate(p,pt) for p in F]==image for pt in points)
 assert all(all(sum(a*b for a,b in zip(e,W))==W[i] for e in p) for i,p in enumerate(F))
 groups=defaultdict(list)
 for j in range(3,12):groups[W[j]].append(j+1)
 charts={j:[a for a in group if a<j] for group in groups.values() for j in group}
 assert charts=={4:[],12:[4],5:[],6:[],7:[],10:[7],11:[7,10],8:[],9:[8]}
 # Exact transfer identities to the old literal graphs, no numerical parameter samples.
 transfers=[]
 for pivot,others,bad in [(10,[9,6],2),(11,[3],1)]:
  G,graph,X=make_graph(F,pivot,others);n=len(X);np=n-12
  base,_,_=make_graph(F,pivot,[]);base={i:{e+(0,)*np:c for e,c in p.items()} for i,p in base.items()}
  assert all(G[j]==base[j] for j in G if j!=bad)
  if pivot==10:
   expected=add(scale(mul(X[12],mul(G[6],G[9])),4),scale(mul(X[13],mul(G[6],G[6])),4))
  else:expected=scale(mul(X[12],mul(G[6],G[3])),12)
  assert add(G[bad],base[bad],Q(-1))==expected
  collision_check(G,points,image,pivot,np);transfers.append({'pivot':pivot+1,'bad':bad+1,'identity_checked':True})
 records=[];model_hashes={};dual_hashes={}
 for pivot,other,bad,dual_name in [(9,6,1,'dual-adaptive.json'),(8,7,2,'dual.json')]:
  folder=HERE/'graph-pencils12'/('pivot-%02d-D7'%(pivot+1));model_path=folder/'model.json';dual_path=folder/dual_name
  model=json.loads(model_path.read_text());dual=json.loads(dual_path.read_text());D=model['source_order_cap'];assert D==7 and model['pivot']==pivot+1 and model['other']==other+1 and model['bad']==bad+1 and dual['model_sha256']==sha(model_path)
  G,graph,X=make_graph(F,pivot,[other]);assert graph==decode(model['graph'],13)
  assert G=={int(j)-1:decode(p,13) for j,p in model['raw_outputs'].items()}
  for j,p in G.items():
   assert not any(c for e,c in p.items() if sum(e[:12])==0)
   assert {e:c for e,c in p.items() if sum(e[:12])==1}==X[j]
   if j!=bad:assert max(sum(e[:12]) for e in p)<=3
  collision_check(G,points,image,pivot,1)
  C={int(a)-1:decode(p,13) for a,p in model['derivation_coefficients'].items()}
  assert set(C)==set(G) and all(e[pivot]==0 and sum(e[:12])<=D for p in C.values() for e in p)
  def apply(p):
   out={}
   for a,q in C.items():out=add(out,trunc(mul(q,derivative(p,a)),D))
   return out
  one={(0,)*13:Q(1)}
  # These identities certify the UNIQUE actual inverse-Jacobian jets because JG(0)=I.
  for j,p in G.items():assert apply(p)==(one if j==bad else {})
  for a,p in C.items():assert {e:c for e,c in p.items() if sum(e[:12])==0}==(one if a==bad else {})
  defect=add(one,C[bad],Q(-1));assert defect==decode(model['defect'],13)
  weights={tuple(w['source_exponent']):{int(k):Q(c) for k,c in w['parameter_coefficients'].items()} for w in dual['weights']}
  assert len(weights)==len(dual['weights']) and all(len(e)==12 and e[pivot]==0 and sum(e)<=7 and w and all(c and 0<=k<=2 for k,c in w.items()) for e,w in weights.items())
  def lam(p):
   out={}
   for e,c in p.items():
    if e[:12] not in weights:continue
    for k,v in weights[e[:12]].items():out[e[12]+k]=out.get(e[12]+k,Q(0))+c*v
   return {k:c for k,c in out.items() if c}
  source=[]
  for degree in (2,3):
   for word in combinations_with_replacement([a for a in range(12) if a!=pivot],degree):
    e=[0]*13
    for a in word:e[a]+=1
    source.append(tuple(e[:12]));assert lam(apply({tuple(e):Q(1)}))=={}
  assert len(source)==352 and source==list(map(tuple,model['source_monomials']))
  assert lam(defect)=={0:Q(1)}
  # The uniform identity covers the formerly omitted rational value as well.
  assert sum(c*Q(3,8)**k for k,c in lam(defect).items())==1
  record={'pivot':pivot+1,'bad':bad+1,'source_order':7,'basis_annihilations':352,'unit_pairing':True,'generic_parameter_collision_checked':True,'inverse_jet_identities':11,'functional_source_support':len(weights),'parameter_terms':sum(len(w) for w in weights.values()),'maximum_parameter_degree':max(k for w in weights.values() for k in w)};records.append(record)
  model_hashes[str(pivot+1)]=sha(model_path);dual_hashes[str(pivot+1)]=sha(dual_path)
  print('PENCIL',pivot+1,'PASS',record,flush=True)
 assert sum(r['functional_source_support'] for r in records)==101 and sum(r['parameter_terms'] for r in records)==144
 # Retained all-degree literal theorem is a pinned premise, not a newly replayed search.
 old_path=HERE/'graph36-validation.json';old=json.loads(old_path.read_text())
 assert old['status']=='pass' and old['all_nine_arbitrary_degree_certificates'] and old['input_sha256']==sha(path) and old['source_sha256']==sha(HERE/'verify_graph36.py') and old['arithmetic_sha256']==sha(HERE/'verify.py')
 assert {r['pivot'] for r in old['arbitrary_target_degree_cases']}==set(range(4,13))
 for j,pins in old['case_hashes'].items():
  assert pins['raw_sha256']==sha(HERE/'graph36'/('pivot-%02d.json'%int(j)))
  assert pins['arbitrary_degree_sha256']==sha(HERE/'adaptive36-wide'/('pivot-%02d.json'%int(j)))
 receipt={'status':'pass','claim':'KCUB27','complete_projective_weight_charts':9,'weight_blocks':{str(k):v for k,v in groups.items()},'literal_KCUB4_charts':[4,5,6,7,8],'transfer_charts':transfers,'essential_pencils':records,'new_basis_annihilations':704,'new_unit_pairings':2,'parameter_denominators_absent':True,'polynomial_target_degree_unbounded_for_all_charts':True,'rational_target_degree_unbounded_for_two_essential_pencils':True,
 'checker_sha256':sha(Path(__file__)),'arithmetic_sha256':sha(HERE/'verify.py'),'proof_sha256':sha(HERE/'WEIGHTED_GRAPH_COMPLETION.md'),'input_sha256':sha(path),'literal_validation_sha256':sha(old_path),'literal_checker_sha256':sha(HERE/'verify_graph36.py'),'model_sha256':model_hashes,'dual_sha256':dual_hashes,
 'boundary':'New graphs, collisions, inverse jets, polynomial duals and transfer identities are independently reconstructed with Fraction arithmetic. Projective coverage, graph Keller inheritance and the uncapped target consequence are written. KCUB4 literal certificates are imported as a pinned premise. No general coordinated/ordered target or global dimension exclusion.','positive_bounds_unchanged':True,'goal_complete':False}
 (HERE/'graph-pencils12/validation.json').write_text(json.dumps(receipt,indent=2)+'\n')
 print('PASS nine weighted projective charts;704 annihilations and2 unit pairings;all polynomial target degrees;global bounds unchanged.')
if __name__=='__main__':main()

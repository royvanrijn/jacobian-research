#!/usr/bin/env python3
"""Signed generic rank-one extension controls; standard library by default.

The frozen373 package supplies pinned exact Laurent/PD helpers and its first
core contradiction. This checker verifies three additional cores, three raw
local folds, both signed basis determinants and12 actual signed transfers.
Optional --census re-reads the independent named census and actual signed
Spherogram knot_group model on all16 literal diagrams. No parameter sampling.
The formal coefficient families are independent roots and independent
covectors; they need not be bases of the ambient space/dual. No inverse
of their evaluation matrix occurs. The paired-span and nondegenerate-form
consequences also use their separate written linear-algebra arguments.
"""
import argparse
from copy import deepcopy
import hashlib
import importlib.util
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
INPUT=HERE/'generic_kch_extensions_inputs.json'
DATA=json.loads(INPUT.read_text())
for path,sha in DATA['source_sha256'].items():
 assert hashlib.sha256((HERE/path).read_bytes()).hexdigest()==sha,path
spec=importlib.util.spec_from_file_location('generic_kch_frozen_base',HERE/'verify_generic_kch_obstruction.py')
BASE=importlib.util.module_from_spec(spec);spec.loader.exec_module(BASE)
Ring=BASE.Ring
SIGNED={c['name']:{'seed_arcs':c['seed_arcs'],'signed_crossings':c['canonical_signed_relations']} for c in DATA['cases']}
EXPECTED_CHANGED={'13a_373':[1,2,6,8],'13a_375':[1,2,6,9],
 '13a_391':[1,2,6,9],'13a_392':[1,2,6,12],
 '13a_3264':[1,3,6,10],'13a_3266':[1,3,6,9],
 '13a_3275':[1,3,6,9],'13a_3276':[1,3,6,9],
 '13a_304':[0,1,6,9],'13a_362':[0,1,6,11],
 '13a_623':[0,1,5,8],'13a_629':[0,1,5,8],
 '13a_3232':[0,1,6,9],'13a_3259':[0,1,6,11],
 '13a_3374':[0,1,5,8],'13a_3375':[0,1,5,8]}

def setup(names,size):
 R=Ring(names);c=R.c;mu=R.v('mu');mi=R.v('mu',-1);d=R.sub(c(1),mu)
 U=[[c(int(i==j)) for j in range(size)] for i in range(size)]
 return R,c,mu,mi,d,U

def engine(R,B):
 add,sub,mul=R.add,R.sub,R.mul
 def dot(f,u):return add(*(mul(f[i],B[i][j],u[j]) for i in range(len(B)) for j in range(len(B))))
 def conj(o,p,e):
  u,f=o;v,g=p;a=dot(f,v);b=dot(g,u);mi=R.v('mu',-1)
  if e==1:out=([sub(x,mul(a,y)) for x,y in zip(v,u)],[add(x,mul(mi,b,y)) for x,y in zip(g,f)])
  else:out=([add(x,mul(mi,a,y)) for x,y in zip(v,u)],[sub(x,mul(b,y)) for x,y in zip(g,f)])
  assert dot(out[1],out[0])==B[0][0]
  return out
 return dot,conj

def signed_step(name,o,a,b):
 rel=next(r for r in SIGNED[name]['signed_crossings'] if r[0]==o and {r[1],r[2]}=={a,b})
 return o,a,b,rel[3] if rel[1]==a else -rel[3]

def paired_siblings():
 R,c,mu,mi,d,E=setup('mu x X y Y z Z',3);x,X,y,Y,z,Z=[R.v(s) for s in ['x','X','y','Y','z','Z']]
 add,sub,mul,neg=R.add,R.sub,R.mul,R.neg
 B=[[d,x,y],[X,d,z],[Y,Z,d]];dot,conj=engine(R,B);U,V,W=E
 pair=lambda u:(u,[dict(a) for a in u])
 plans={375:[(6,9,10),(10,7,8),(4,10,11),(11,4,5),(8,11,12)],391:[(9,6,5),(6,9,10),(4,10,11),(5,11,12),(11,7,8)],392:[(12,6,5),(5,12,11),(4,11,10),(6,10,9),(10,7,8)]}
 reports=[]
 for index,plan in plans.items():
  name='13a_'+str(index);roots={0:pair(U),4:pair([neg(a) for a in U]),6:pair(V),7:pair(V),12 if index==392 else 9:pair(W)}
  actual=[]
  for step in plan:
   o,a,b,e=signed_step(name,*step);assert o in roots and a in roots and b not in roots
   roots[b]=conj(roots[o],roots[a],e);actual.append((o,a,b,e))
  if index==375:
   T=conj(pair(V),pair(W),-1);A=conj(pair(U),T,1)
   t=dot(U,T[0]);tau=dot(T[1],U);prod=add(c(-1),neg(mul(mi,t,tau)))
   assert roots[10]==T and roots[11]==A and roots[8]==conj(T,pair(V),-1)
   assert roots[5]==([add(mul(prod,a),mul(mi,tau,b)) for a,b in zip(U,T[0])],[sub(mul(prod,a),mul(t,b)) for a,b in zip(U,T[1])])
   assert roots[12][0][0]==neg(t) and roots[12][1][0]==mul(mi,tau)
  elif index==391:
   T=conj(pair(V),pair(W),1);A=conj(pair(U),T,-1);D=conj(pair(W),pair(V),-1)
   t=dot(U,T[0]);tau=dot(T[1],U)
   assert roots[10]==T and roots[11]==A and roots[5]==D
   assert roots[12]==conj(D,A,1) and roots[8]==conj(A,pair(V),1)
   assert roots[12][0][0]==mul(mi,t) and roots[12][1][0]==neg(tau)
   assert conj(A,roots[8],-1)==pair(V)
  else:
   h=add(c(1),mul(mi,z,Z));t=add(neg(mul(mi,z,x)),mul(h,y));tau=add(mul(Z,X),mul(h,Y))
   k=add(mul(mi,t,X),mul(z,sub(h,mul(mi,d))));K=add(neg(mul(tau,x)),mul(Z,add(d,h)))
   er={11:[{},neg(mul(mi,z)),h],10:[mul(mi,t),neg(mul(mi,z)),h],9:[mul(mi,t),sub(neg(mul(mi,z)),k),h],8:[mul(mi,mi,K,t),sub(c(1),mul(mi,mi,K,z)),mul(mi,K,h)]}
   ef={11:[{},Z,h],10:[neg(tau),Z,h],9:[neg(tau),add(Z,mul(mi,K)),h],8:[mul(k,tau),sub(c(1),mul(k,Z)),neg(mul(k,h))]}
   for i in er:assert roots[i]==(er[i],ef[i])
  reports.append({'name':name,'paired_trace':actual,'all_raw_directed_pairings_free':True,'generic_root_and_covector_formulas':True})
 # Auxiliary constant-one certificate for392 after the legal nonzero cancellations.
 R=Ring('mu beta h z Z k K p');c=R.c;mu,bi,h,z,Z,k,K,p=[R.v(n) for n in R.names];mi=R.v('mu',-1);beta_inv=R.v('beta',-1);add,sub,mul=R.add,R.sub,R.mul
 fs=[sub(sub(mul(mu,h),mu),mul(z,Z)),add(k,mul(mi,z)),add(K,mul(mu,Z)),sub(K,mul(mu,bi)),add(k,beta_inv),sub(mul(p,h),c(1))]
 hs=[mul(p,mi),mul(p,Z),R.neg(mul(p,mi,k)),mul(p,mi,k),mul(p,bi),c(-1)]
 assert add(*(mul(a,b) for a,b in zip(hs,fs)))==c(1)
 corrupted=[dict(h) for h in hs];corrupted[0]=add(corrupted[0],c(1))
 assert add(*(mul(a,b) for a,b in zip(corrupted,fs)))!=c(1)
 return {'systems':reports,'corrupted_auxiliary_multiplier_rejected':True,'392_auxiliary_unit_identity':1,'only_units':['mu','beta']}

def local_folds():
 R,c,mu,mi,d,E=setup('mu a b x X z Z kap',3);a,b,x,X,z,Z,kap=[R.v(s) for s in ['a','b','x','X','z','Z','kap']];ki=R.v('kap',-1)
 add,sub,mul,neg=R.add,R.sub,R.mul,R.neg;A,B,C=E;pair=lambda u:(u,[dict(a) for a in u])
 raw=[[d,a,x],[b,d,z],[X,Z,d]];dot,conj=engine(R,raw)
 p0=conj(pair(B),pair(A),1);p4=conj(pair(A),pair(B),-1);p2=conj(p0,pair(B),1);k=add(c(1),mul(mi,a,b))
 assert p2==([neg(mul(mi,a)),k,{}],[b,k,{}])
 p7=conj(p2,pair(C),1);r=dot(p2[1],C);s=dot(C,p2[0]);p=dot(p7[1],A);q=dot(A,p7[0]);closing=conj(p7,pair(A),1)
 assert p7==([mul(mi,a,r),neg(mul(k,r)),c(1)],[mul(mi,b,s),mul(mi,k,s),c(1)])
 assert closing==([sub(c(1),mul(p,mi,a,r)),mul(p,k,r),neg(p)],[add(c(1),mul(mi,mi,q,b,s)),mul(mi,mi,q,k,s),mul(mi,q)])
 special=[[d,neg(mul(mu,ki)),{}],[kap,d,z],[{},Z,d]];dd,cj=engine(R,special)
 q0=cj(pair(B),pair(A),1);q4=cj(pair(A),pair(B),-1);q2=cj(q0,pair(B),1);q7=cj(q2,pair(C),1)
 assert q2==([mul(ki,u) for u in A],[mul(kap,u) for u in A]) and q7==pair(C)
 assert q4==([neg(mul(ki,u)) for u in q0[0]],[neg(mul(kap,u)) for u in q0[1]])
 for typ in ['I','II']:
  p2=conj(pair(A),pair(B),-1) if typ=='I' else conj(pair(C),pair(B),-1)
  p3=conj(pair(C),p2,-1) if typ=='I' else conj(pair(A),p2,-1)
  closing=conj(p3,pair(A),-1)
  if typ=='I':
   t=add(Z,mul(mi,a,X));T=sub(z,mul(b,x));H=sub(mul(mu,b),mul(T,X));K=mul(mi,add(a,mul(t,x)))
   assert p3==([mul(mi,a),c(1),mul(mi,t)],[neg(b),c(1),neg(T)])
   assert closing==([add(c(1),mul(mi,mi,a,H)),mul(mi,H),mul(mi,mi,t,H)],[add(c(1),mul(b,K)),neg(K),mul(K,T)])
   special=[[d,neg(mul(mu,ki)),x],[kap,d,mul(kap,x)],[X,mul(ki,X),d]]
  else:
   t=add(a,mul(mi,Z,x));T=sub(b,mul(z,X))
   assert p3==([mul(mi,t),c(1),mul(mi,Z)],[neg(T),c(1),neg(z)])
   assert closing==([add(c(1),mul(mi,t,T)),T,mul(mi,T,Z)],[add(c(1),mul(mi,t,T)),neg(mul(mi,t)),mul(mi,t,z)])
   special=[[d,neg(mul(mu,ki)),x],[kap,d,{}],[X,{},d]]
  dd,cj=engine(R,special)
  q2=cj(pair(A),pair(B),-1) if typ=='I' else cj(pair(C),pair(B),-1)
  q3=cj(pair(C),q2,-1) if typ=='I' else cj(pair(A),q2,-1)
  q4=cj(pair(B),q3,-1);qj=cj(q2,pair(C),-1)
  assert q4==([neg(mul(ki,u)) for u in A],[neg(mul(kap,u)) for u in A])
  assert qj==pair(C)
  assert cj(q3,pair(A),-1)==([mul(kap,u) for u in B],[mul(ki,u) for u in B])
 return {'alternate_old_four_fold':True,'A6_negative_local_types':['I','II'],'only_inverted_parameters':['mu','kap'],'all_directed_pairings_independent':True}

def generic_basis_changes():
 specs=[{'name':c['name'],'original_seeds':c['seed_arcs'],'changed_seeds':c['changed_seed_arcs'],'literal_pair_derivations':[t[:3] for t in c['signed_basis_derivations']]} for c in DATA['cases']]
 names='mu '+' '.join('a%d%d'%(i,j) for i in range(4) for j in range(4) if i!=j)
 R,c,mu,mi,d,E=setup(names,4);B=[[d if i==j else R.v('a%d%d'%(i,j)) for j in range(4)] for i in range(4)];dot,conj=engine(R,B);reports=[]
 for spec in specs:
  name=spec['name'];old=spec['original_seeds'];new=spec['changed_seeds'];assert old==SIGNED[name]['seed_arcs']
  roots={a:(E[i],[dict(x) for x in E[i]]) for i,a in enumerate(old)};actual=[]
  for o,a,b in spec['literal_pair_derivations']:
   o,a,b,e=signed_step(name,o,a,b);assert o in roots and a in roots and b not in roots
   roots[b]=conj(roots[o],roots[a],e);actual.append((o,a,b,e))
  us=[roots[a][0] for a in new];fs=[roots[a][1] for a in new]
  assert us[:2]==E[:2] and fs[:2]==E[:2]
  assert R.sub(R.mul(us[2][2],us[3][3]),R.mul(us[2][3],us[3][2]))==c(1)
  assert R.sub(R.mul(fs[2][2],fs[3][3]),R.mul(fs[2][3],fs[3][2]))==c(1)
  case=next(c for c in DATA['cases'] if c['name']==name)
  assert [list(t) for t in actual]==case['signed_basis_derivations']
  assert new==EXPECTED_CHANGED[name]
  reports.append({'name':name,'original_seeds':old,'changed_seeds':new,'actual_signed_derivations':actual,'root_determinant':1,'covector_determinant':1})
 return reports

def signed_core_transfers(records=None):
 if records is None:records=DATA['signed_reduced_transfers']
 assert len(records)==12 and {r['new'] for r in records}==set(EXPECTED_CHANGED)-{'13a_373','13a_375','13a_391','13a_392'}
 def key(r):
  o,a,b,e=r;return o,min(a,b),max(a,b),e if a<b else -e
 def core(name,j=6,jj=7):
  fold=lambda i:0 if i==4 else j if i==jj else i
  return {key((fold(o),fold(a),fold(b),e)) for o,a,b,e in SIGNED[name]['signed_crossings'] if not {o,a,b}&{1,2,3}}
 for r in records:
  name=r['new'];j,jj=(8,9) if name in ['13a_623','13a_629','13a_3374','13a_3375'] else(6,7)
  source=core(name,j,jj);m=r['relabel'];image={key((int(m[str(o)]),int(m[str(a)]),int(m[str(b)]),e)) for o,a,b,e in source} if m else source
  target={(o,a,b,e*r['global_sign']) for o,a,b,e in core(r['target'])}
  assert image==target and r['global_sign'] in [1,-1]
 # Bind each local generic fold to actual signed crossings, not to old unsigned data.
 for name in ['13a_373','13a_375','13a_391','13a_392']:
  steps=[(2,6,7),(7,2,3),(3,1,0),(1,3,4),(0,1,2)]
  assert [signed_step(name,*t)[3] for t in steps]==[1,1,1,-1,-1]
 for name in ['13a_3264','13a_3266','13a_3275','13a_3276']:
  steps=[(3,1,0),(1,3,4),(0,3,2),(2,6,7),(7,1,2)]
  assert [signed_step(name,*t)[3] for t in steps]==[1,-1,1,1,1]
 for name in ['13a_304','13a_362','13a_623','13a_629','13a_3232','13a_3259','13a_3374','13a_3375']:
  j,jj=(8,9) if name in ['13a_623','13a_629','13a_3374','13a_3375'] else (6,7)
  typ1=name in ['13a_304','13a_362','13a_623','13a_629']
  steps=[(0,1,2),(j,2,3),(3,0,1),(1,3,4),(2,j,jj)] if typ1 else [(j,1,2),(0,2,3),(3,0,1),(1,3,4),(2,j,jj)]
  assert [signed_step(name,*t)[3] for t in steps]==[-1]*5
 return records


def case_binding(case):
 name=case['name'];assert name in EXPECTED_CHANGED
 source=json.loads((HERE/case['source_input']).read_text())
 original=next(c for c in source['cases'] if c['name']==name)
 for key in ['name','pd','seed_arcs','seed_trace','census_pd','census_DT','census_map_witness']:
  assert case[key]==original[key],(name,key)
 signed,_=BASE.diagram(case['pd']);assert len(signed)==13
 assert signed==case['signed_crossings']
 assert [[r['over'],r['input'],r['output'],r['exponent']] for r in signed]==case['canonical_signed_relations']
 BASE.census_bijection(case['pd'],case['census_pd'],case['census_map_witness'])
 known=set(case['seed_arcs']);used=set();assert len(known)==4
 for o,a,b,e in case['signed_trace']:
  assert o in known and a in known and b not in known
  hits=[r for r in signed if r['over']==o and {r['input'],r['output']}=={a,b}]
  assert len(hits)==1;r=hits[0]
  assert e==(r['exponent'] if a==r['input'] else -r['exponent'])
  known.add(b);used.add(r['crossing'])
 assert known==set(range(13)) and len(used)==9
 left=[[r['over'],r['input'],r['output'],r['exponent']] for r in signed if r['crossing'] not in used]
 assert left==case['remaining_signed_crossings'] and len(left)==4
 return {'name':name,'signed_crossings':13,'signed_seed_steps':9,'named_planar_map_checked':True}


def actual_signed_model(case,Link):
 # Independent port-map recovery by incidence propagation from the constructor.
 pd=case['pd'];link=Link(pd);ends={}
 for i,row in enumerate(pd):
  for j,e in enumerate(row):ends.setdefault(e,[]).append((i,j))
 alpha={a:b for pair in ends.values() for a,b in (pair,pair[::-1])}
 solutions=[]
 for start in (0,2):
  rotations={0:start};queue=[0];valid=True
  for i in queue:
   for j,(other,k) in enumerate(link.crossings[i].adjacent):
    h,l=alpha[i,(j+rotations[i])%4]
    if h!=other.label:valid=False;break
    r=(l-k)%4
    if r not in (0,2):valid=False;break
    if h in rotations:
     if rotations[h]!=r:valid=False;break
    else:rotations[h]=r;queue.append(h)
   if not valid:break
  if valid and len(rotations)==13:solutions.append(rotations)
 assert len(solutions)==1;rotation=solutions[0]
 parent={e:e for e in ends}
 def find(e):
  while parent[e]!=e:e=parent[e]
  return e
 for a,b,c,d in pd:parent[find(b)]=find(d)
 components={}
 for e in parent:components.setdefault(find(e),[]).append(e)
 components=sorted(sorted(c) for c in components.values())
 arc={e:i for i,c in enumerate(components) for e in c}
 pieces=link._pieces();piece_arc={}
 for i,piece in enumerate(pieces):
  labels={arc[pd[c.label][(j+rotation[c.label])%4]] for c,j in piece}
  assert len(labels)==1;piece_arc[i]=labels.pop()
 actual=[]
 for crossing in link.crossings:
  found={}
  for i,piece in enumerate(pieces):
   for position,(c,j) in enumerate(piece):
    if c!=crossing:continue
    typ='output' if position==0 else 'input' if position==len(piece)-1 else 'over'
    if typ in found:assert found[typ]==i
    found[typ]=i
  assert set(found)=={'input','output','over'}
  actual.append([piece_arc[found['over']],piece_arc[found['input']],piece_arc[found['output']],-crossing.sign])
 assert actual==case['canonical_signed_relations']
 return True


def census_replay(data):
 import importlib.metadata as metadata
 import snappy_15_knots
 from spherogram import Link
 source=data['census_provenance']
 assert metadata.version('snappy_15_knots')==source['version']
 assert metadata.version('spherogram')==source['spherogram_version']
 root=Path(snappy_15_knots.__file__).parent
 for pin in source['source_files']:
  raw=(root/pin['path']).read_bytes()
  assert len(raw)==pin['bytes'] and hashlib.sha256(raw).hexdigest()==pin['sha256']
 for case in data['cases']:
  link=Link('K'+case['name'].replace('_',''))
  assert [list(row) for row in link.PD_code()]==case['census_pd']
  assert json.loads(json.dumps(link.DT_code()))==case['census_DT']
  actual_signed_model(case,Link)
 return {'named_census_rereads':16,'actual_signed_knot_group_models':16}


def main():
 parser=argparse.ArgumentParser();parser.add_argument('--census',action='store_true');args=parser.parse_args()
 assert len(DATA['names'])==len(set(DATA['names']))==16
 assert set(DATA['names'])==set(EXPECTED_CHANGED)=={c['name'] for c in DATA['cases']}
 bindings=[case_binding(c) for c in DATA['cases']]
 base_data=json.loads((HERE/'generic_kch_obstruction_inputs.json').read_text())
 first_core={}
 for fn in (BASE.raw_first_loop,BASE.paired_core,BASE.unit_certificate):first_core.update(fn(base_data))
 cores=paired_siblings();folds=local_folds();changes=generic_basis_changes();transfers=signed_core_transfers()
 rejected=[]
 for kind in ('crossing_sign','named_map'):
  bad=deepcopy(DATA['cases'][0])
  if kind=='crossing_sign':bad['signed_trace'][0][3]*=-1
  else:bad['census_map_witness']['dart_images'][0]=bad['census_map_witness']['dart_images'][1]
  try:case_binding(bad)
  except AssertionError:rejected.append(kind)
  else:raise AssertionError('corrupted case accepted: '+kind)
 bad=deepcopy(DATA['signed_reduced_transfers']);bad[0]['global_sign']*=-1
 try:signed_core_transfers(bad)
 except AssertionError:rejected.append('global_transfer_sign')
 else:raise AssertionError('corrupted transfer sign accepted')
 census=census_replay(DATA) if args.census else None
 files=[Path(__file__),INPUT,HERE/'GENERIC_KCH_EXTENSIONS.tex']
 out={'status':'PASS','names':DATA['names'],'literal_signed_bindings':bindings,
      'frozen373_first_core_controls':first_core,'additional_core_controls':cores,
      'raw_local_folds':folds,'generic_signed_basis_changes':changes,
      'signed_reduced_transfers':transfers,'corruptions_rejected':rejected,
      'optional_census_replay':census,
      'representation_scope':DATA['representation_scope'],'augmentation_scope':DATA['augmentation_scope'],
      'paired_span_scope':DATA['paired_span_scope'],'form_preserving_scope':DATA['form_preserving_scope'],
      'separate_corollary_dependencies':DATA['separate_corollary_dependencies'],
      'input_sha256':{p.relative_to(HERE).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in files}|DATA['source_sha256'],
      'assurance':'Exact internal integer/Laurent and literal oriented-PD controls, including both basis determinants for separately independent root and covector families and no inverse evaluation matrix. Written span, independence, nonzero-cancellation and nondegenerate-form adjoint arguments remain essential. Positive attainment and meridional-rank statements are separate mathematical imports for their own corollaries; no formal verification, independent specialist review or worldwide-priority claim.'}
 receipt_path=HERE/'generic_kch_extensions_verification.json'
 if receipt_path.is_file():
  previous=json.loads(receipt_path.read_text())
  if 'internal_pdf_reviews' in previous:
   # Preserve each review's own source/PDF hashes; never transfer its
   # assurance automatically to a later source or PDF edition.
   out['internal_pdf_reviews']=previous['internal_pdf_reviews']
 receipt_path.write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps({'status':'PASS','names':16,'additional_generic_cores':3,'raw_local_fold_types':3,
                   'both_signed_basis_determinants':16,'signed_reduced_transfers':12,
                   'representation_fields':'all','meridian_determinant':'every nonzero scalar',
                   'written_paired_span_consequence':'min(r,q)<=3 and r+q<=7',
                   'optional_census_replay':census},indent=2))

if __name__=='__main__':main()

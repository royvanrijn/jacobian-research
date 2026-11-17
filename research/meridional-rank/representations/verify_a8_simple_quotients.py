#!/usr/bin/env python3
"""Twenty-four positive A8 double-transposition epimorphisms.

Default replay uses only the standard library. Optional --census rereads
all24 literal named PD/DT records from the pinned external census.
No cover, homology, older affine bound or full293 negative scan is used.
"""
import argparse
from copy import deepcopy
import hashlib
import itertools
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
INPUT=HERE/'a8_simple_quotient_inputs.json'
NAMES=['13a_373','13a_375','13a_391','13a_392','13a_3264','13a_3266','13a_3275','13a_3276',
       '13n_1005','13n_1006','13n_1007','13n_1011','13n_1012','13n_1013','13n_1043','13n_1044',
       '13n_1045','13n_1046','13n_1047','13n_1048','13n_3842','13n_3844','13n_3855','13n_3856']
ID=tuple(range(8))

def compose(a,b):return tuple(a[b[i]] for i in range(8))
def inverse(a):return tuple(a.index(i) for i in range(8))
def even(a):return sum(a[i]>a[j] for i in range(8) for j in range(i+1,8))%2==0

def valid_perm(a):
 assert len(a)==8 and sorted(a)==list(range(8)) and even(a)

def class_member(a):
 valid_perm(a)
 assert compose(a,a)==ID and sum(a[i]!=i for i in range(8))==4

def all_double_transpositions():
 out=set()
 for a,b,c,d in itertools.combinations(range(8),4):
  for pairs in [[(a,b),(c,d)],[(a,c),(b,d)],[(a,d),(b,c)]]:
   p=list(ID)
   for i,j in pairs:p[i],p[j]=p[j],p[i]
   out.add(tuple(p))
 assert len(out)==210
 return out

_GROUP_CACHE={}
_ELEMENTS_CACHE={}
def generated_group(generators):
 gens=tuple(sorted(set(tuple(p) for p in generators)))
 for p in gens:valid_perm(p)
 if gens in _GROUP_CACHE:return _GROUP_CACHE[gens]
 steps=tuple(sorted(set(gens+tuple(inverse(g) for g in gens))))
 seen={ID};queue=[ID]
 for a in queue:
  for b in steps:
   c=compose(a,b)
   if c not in seen:seen.add(c);queue.append(c)
 assert len(seen)<=20160
 result={'order':len(seen),'generators':len(gens),'steps_including_inverses':len(steps)}
 _GROUP_CACHE[gens]=result
 _ELEMENTS_CACHE[gens]=tuple(seen)
 return result

def simplicity_control(pair):
 gens=tuple(sorted(set(tuple(p) for p in pair)))
 generated_group(pair)
 unseen=set(_ELEMENTS_CACHE[gens]);sizes=[]
 steps=tuple(sorted(set(gens+tuple(inverse(g) for g in gens))))
 while unseen:
  first=min(unseen);orbit={first};queue=[first]
  for x in queue:
   for a in steps:
    c=compose(compose(a,x),inverse(a))
    if c not in orbit:orbit.add(c);queue.append(c)
  assert orbit<=unseen
  unseen-=orbit;sizes.append(len(orbit))
 sizes.sort()
 assert sizes==[1,105,112,210,1120,1260,1344,1344,1344,1680,2520,2880,2880,3360]
 sums={0}
 for n in sizes[1:]:sums|={x+n for x in sums}
 proper_divisor_unions=[1+x for x in sums if 1<1+x<20160 and 20160%(1+x)==0]
 assert not proper_divisor_unions
 return {'conjugacy_class_sizes':sizes,'proper_nontrivial_normal_union_divisors':proper_divisor_unions,
         'argument':'A normal subgroup is a union of conjugacy classes containing identity, and its order divides20160.'}

def diagram(pd):
 ends={};adj={x:set() for row in pd for x in row}
 for i,row in enumerate(pd):
  assert len(row)==4 and len(set(row))==4
  for j,e in enumerate(row):ends.setdefault(e,[]).append(4*i+j)
  b,d=row[1],row[3];adj[b].add(d);adj[d].add(b)
 assert len(ends)==2*len(pd) and all(len(v)==2 for v in ends.values())
 alpha={a:b for pair in ends.values() for a,b in (pair,pair[::-1])}
 unseen=set(alpha);faces=0
 while unseen:
  start=x=min(unseen)
  while True:
   unseen.remove(x);y=alpha[x];x=4*(y//4)+(y+1)%4
   if x==start:break
  faces+=1
 assert faces==len(pd)+2
 x=0;strand=set()
 for _ in range(2*len(pd)):
  assert x not in strand;strand.add(x);y=4*(x//4)+(x+2)%4;x=alpha[y]
 assert x==0 and len(strand)==2*len(pd) and len({x//4 for x in strand})==len(pd)
 comps=[];unused=set(adj)
 while unused:
  cc=set();queue=[min(unused)]
  for x in queue:
   if x in cc:continue
   cc.add(x);queue.extend(adj[x]-cc)
  comps.append(sorted(cc));unused-=cc
 comps.sort();assert len(comps)==len(pd)
 arc={v:i for i,c in enumerate(comps) for v in c}
 return [(arc[b],arc[a],arc[c]) for a,b,c,d in pd],alpha

def census_bijection(pd,census_pd,cert):
 _,a=diagram(pd);_,b=diagram(census_pd)
 f=cert['dart_images'];orientation=cert['orientation'];swap=cert['over_under_swap']
 assert orientation in (-1,1) and swap in (0,1)
 assert len(pd)==len(census_pd) and sorted(f)==list(range(4*len(pd)))
 for x,y in enumerate(f):
  assert f[a[x]]==b[y]
  assert (x%2)^(y%2)==swap
  assert f[4*(x//4)+(x+1)%4]==4*(y//4)+(y+orientation)%4

def group_control(data):
 assert data['target_group']=='A8' and data['target_order']==20160
 assert data['class_size']==210 and data['ordinary_rank']==2 and data['meridian_class_rank']==4
 pair=data['ordinary_generators'];assert len(pair)==2
 ordinary=generated_group(pair);assert ordinary['order']==20160
 cls=all_double_transpositions();first=min(cls);orbit={first};queue=[first]
 for x in queue:
  for a in pair:
   a=tuple(a);c=compose(compose(a,x),inverse(a))
   if c not in orbit:orbit.add(c);queue.append(c)
 assert orbit==cls
 quartet=data['attaining_class_quartet'];assert len(quartet)==4
 for p in quartet:class_member(p)
 attained=generated_group(quartet);assert attained['order']==20160
 return {'target_order':20160,'class_size':210,'single_A8_conjugacy_orbit_checked':210,
         'ordinary_pair_group_order':ordinary['order'],'attaining_class_quartet_order':attained['order'],
         'simplicity':simplicity_control(pair),
         'class_rank_lower':'Any3 double transpositions supply at most6 edges; a connected graph on8 vertices needs at least7.'}

def check_case(case):
 assert case['name'] in NAMES and len(case['pd'])==13
 crossings,_=diagram(case['pd'])
 census_bijection(case['pd'],case['census_pd'],case['census_map_witness'])
 seeds=case['seed_arcs'];known=set(seeds);assert len(seeds)==len(known)==4
 assert len(case['seed_trace'])==9
 for o,a,b in case['seed_trace']:
  assert o in known and a in known and b not in known
  assert (o,a,b) in crossings or (o,b,a) in crossings
  known.add(b)
 assert known==set(range(13))
 perms=[tuple(p) for p in case['permutations']];assert len(perms)==13
 for p in perms:class_member(p)
 for o,a,b in crossings:assert compose(compose(perms[o],perms[a]),perms[o])==perms[b]
 gens=[perms[i] for i in seeds];group=generated_group(gens)
 assert case['expected_image_order']==group['order']==20160
 orbit={0};queue=[0]
 for i in queue:
  for g in gens:
   if g[i] not in orbit:orbit.add(g[i]);queue.append(g[i])
 assert len(orbit)==8
 return {'name':case['name'],'literal_crossings_verified':13,'class_meridians_verified':13,
         'four_seed_coloring_steps':9,'sphere_embedding_and_one_component':True,
         'named_census_dart_map':True,'image_order':20160,'transitive_degree':8,
         'lower_bound':4,'bridge_upper_via_BKVV':4}

def optional_census(data):
 import importlib.metadata as metadata
 import snappy_15_knots
 from spherogram import Link
 src=data['census_provenance']
 assert metadata.version('snappy_15_knots')==src['version']
 assert metadata.version('spherogram')==src['spherogram_version']
 root=Path(snappy_15_knots.__file__).parent
 for pin in src['source_files']:
  raw=(root/pin['path']).read_bytes()
  assert len(raw)==pin['bytes'] and hashlib.sha256(raw).hexdigest()==pin['sha256']
 for c in data['cases']:
  link=Link('K'+c['name'].replace('_',''))
  assert [list(row) for row in link.PD_code()]==c['census_pd']
  assert json.loads(json.dumps(link.DT_code()))==c['census_DT']
 return {'fresh_named_PD_DT_rereads':24,'package':src['package'],'version':src['version']}

def rejected_mutations(data):
 rejected=[]
 for kind in ['arc_permutation','seed_trace','named_dart_map']:
  bad=deepcopy(data['cases'][0])
  if kind=='arc_permutation':bad['permutations'][0]=list(ID)
  elif kind=='seed_trace':bad['seed_trace'][0][2]=bad['seed_arcs'][0]
  else:bad['census_map_witness']['dart_images'][0]=bad['census_map_witness']['dart_images'][1]
  try:check_case(bad)
  except AssertionError:rejected.append(kind)
  else:raise AssertionError('Corruption accepted: '+kind)
 for kind in ['ordinary_pair','class_quartet']:
  bad=deepcopy(data)
  if kind=='ordinary_pair':bad['ordinary_generators'][1]=bad['ordinary_generators'][0]
  else:bad['attaining_class_quartet']=[bad['attaining_class_quartet'][0]]*4
  try:group_control(bad)
  except AssertionError:rejected.append(kind)
  else:raise AssertionError('Corruption accepted: '+kind)
 return rejected

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--census',action='store_true');args=parser.parse_args()
 data=json.loads(INPUT.read_text())
 assert len(data['cases'])==len(NAMES)==24 and {c['name'] for c in data['cases']}==set(NAMES)
 for filename,pin in data['source_sha256'].items():assert hashlib.sha256((HERE/filename).read_bytes()).hexdigest()==pin
 group=group_control(data);cases=[check_case(c) for c in data['cases']];rejected=rejected_mutations(data)
 census=optional_census(data) if args.census else None
 files=[Path(__file__),INPUT,HERE/'A8_SIMPLE_QUOTIENTS.tex']
 receipt={'status':'PASS','group':group,'named_cases':cases,'corruptions_rejected':rejected,
          'optional_census':census,'scope':data['scope'],
          'input_sha256':{p.relative_to(HERE).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in files}|data['source_sha256'],
          'assurance':'Exact internal positive finite-permutation and literal-PD replay. The named census and published BKVV bridge-Wirtinger theorem remain external imports. No older affine lower bound or cover calculation is a premise. No formal verification, external specialist review, new named count or worldwide-priority claim.'}
 path=HERE/'a8_simple_quotient_verification.json'
 if path.is_file():
  old=json.loads(path.read_text())
  if 'internal_pdf_reviews' in old:
   # Each record keeps its own source/PDF hashes, without assurance
   # automatically propagating to a later edition.
   receipt['internal_pdf_reviews']=old['internal_pdf_reviews']
 path.write_text(json.dumps(receipt,indent=2)+'\n')
 print(json.dumps({'status':'PASS','named_epimorphisms':24,'target_order':20160,
                   'class_size':210,'class_rank':4,'ordinary_group_rank':2,
                   'literal_crossing_relations':312,'corruptions_rejected':rejected,
                   'optional_census':census},indent=2))

if __name__=='__main__':main()

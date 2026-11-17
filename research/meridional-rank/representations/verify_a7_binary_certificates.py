#!/usr/bin/env python3
"""Three exact binary affine sharp-four certificates; standard library by default.

Optional --census re-reads the pinned independent named census. Optional
--cover-homology reconstructs complete directed-edge branched-cover complexes
at characteristic2. Neither control is a premise of the direct matrix bound.
"""
import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
NAMES=['13n_4641','13n_4642','13n_4645']
ID=tuple(range(7));UNIT=(ID,(0,)*7)

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

def compose(a,b):return tuple(a[b[i]] for i in range(7))

def rank2(rows):
 basis={}
 for row in rows:
  while row:
   k=row.bit_length()-1
   if k in basis:row^=basis[k]
   else:basis[k]=row;break
 return len(basis)

def full_edge_chain(perms,crossings):
 count=len(perms);n=7;edges=[(g,s) for g in range(count) for s in range(n)]
 parent=list(range(n))
 def root(i):
  while parent[i]!=i:i=parent[i]
  return i
 tree=[]
 # First-edge tree, independent of reverse-tree existing checker.
 for j,(g,s) in enumerate(edges):
  t=perms[g][s]
  if root(s)!=root(t):parent[root(s)]=root(t);tree.append(j)
 assert len(tree)==6
 cells=[]
 for o,a,b in crossings:
  for start in range(n):
   s=start;bits=0
   for g in (o,a,o,b):bits^=1<<(n*g+s);s=perms[g][s]
   assert s==start;cells.append(bits)
 for g,p in enumerate(perms):
  unused=set(range(n))
  while unused:
   start=s=min(unused);bits=0
   while True:
    unused.remove(s);bits^=1<<(n*g+s);s=p[s]
    if s==start:break
   cells.append(bits)
 chords=[j for j in range(len(edges)) if j not in tree]
 collapsed=[sum(((row>>j)&1)<<k for k,j in enumerate(chords)) for row in cells]
 rank=rank2(collapsed)
 return {'edges':len(edges),'vertices':7,'tree_edges':tree,'columns_after_tree_collapse':len(chords),'filling_and_crossing_cells':len(cells),'boundary_rank_F2':rank,'H1_F2_dimension':len(chords)-rank}


def census_bijection(pd,census_pd,cert):
 _,a=diagram(pd);_,b=diagram(census_pd)
 f=cert['dart_images'];orientation=cert['orientation'];swap=cert['over_under_swap']
 assert orientation in (-1,1) and swap in (0,1)
 assert len(pd)==len(census_pd) and sorted(f)==list(range(4*len(pd)))
 for x,y in enumerate(f):
  assert f[a[x]]==b[y]
  assert (x%2)^(y%2)==swap
  assert f[4*(x//4)+(x+1)%4]==4*(y//4)+(y+orientation)%4

def rank_rows(rows):return rank2([sum((v%2)<<j for j,v in enumerate(row)) for row in rows])

def mm(a,b):return [[sum(x*y for x,y in zip(row,col))%2 for col in zip(*b)] for row in a]

def affine_multiply(a,b):
 p=compose(a[0],b[0]);v=list(a[1])
 for i in range(7):v[a[0][i]]^=b[1][i]
 return p,tuple(v)

def affine_inverse(a):
 p=tuple(a[0].index(i) for i in range(7));v=[0]*7
 for i in range(7):v[p[i]]=a[1][i]
 return p,tuple(v)

def affine_word(word,generators):
 result=UNIT
 for i in word:result=affine_multiply(result,generators[i])
 return result

def schreier_structure(generators):
 # Choose lifts t_p of every projected element p by a fixed-order BFS.
 sections={ID:UNIT};queue=[ID]
 for p in queue:
  for s in generators:
   value=affine_multiply(sections[p],s)
   if value[0] not in sections:sections[value[0]]=value;queue.append(value[0])
 # These loops generate the entire kernel, by telescoping any kernel word.
 # Negative edges are inverses of positive edges at their other endpoint.
 rows=[]
 for p in queue:
  for s in generators:
   value=affine_multiply(sections[p],s)
   loop=affine_multiply(value,affine_inverse(sections[value[0]]))
   assert loop[0]==ID and sum(loop[1])%2==0
   rows.append(loop[1])
 dimension=rank_rows(rows)
 return {'projection_order':len(queue),'kernel_dimension':dimension,
         'image_order':len(queue)*2**dimension,'schreier_loops_checked':len(rows)},queue

def check_case(case,cover_homology=False):
 assert case['name'] in NAMES
 pd=case['pd'];crossings,_=diagram(pd)
 assert len(pd)==13
 census_bijection(pd,case['census_pd'],case['census_map_witness'])
 seeds=case['seed_arcs'];known=set(seeds)
 assert len(known)==len(seeds)==4
 for o,a,b in case['seed_trace']:
  assert o in known and a in known and b not in known
  assert (o,a,b) in crossings or (o,b,a) in crossings
  known.add(b)
 assert known==set(range(13))
 ps=case['permutations'];vs=case['arc_translations_F2'];ms=case['homogeneous_W_matrices_F2']
 assert len(ps)==len(vs)==len(ms)==13
 identity=[[int(i==j) for j in range(7)] for i in range(7)]
 displacements=[]
 for p,v,m in zip(ps,vs,ms):
  assert len(p)==len(v)==7 and sorted(p)==list(range(7))
  assert all(x in (0,1) for x in v) and sum(v)%2==0
  assert sum(p[i]!=i for i in range(7))==4 and compose(p,p)==ID
  assert sum(p[i]>p[j] for i in range(7) for j in range(i+1,7))%2==0
  expected=[[0]*7 for _ in range(7)]
  for i in range(6):
   for j in range(6):expected[i][j]=(int(p[j]==i)+int(p[6]==i))%2
   expected[i][6]=v[i]
  expected[6][6]=1
  assert expected==m and mm(m,m)==identity
  delta=[[(a-b)%2 for a,b in zip(r,s)] for r,s in zip(m,identity)]
  assert rank_rows(delta)==2
  displacements.append(delta)
 for o,a,b in crossings:assert mm(mm(ms[o],ms[a]),ms[o])==ms[b]
 assert rank_rows([r for i in seeds for r in displacements[i]])==7
 assert rank_rows([r for delta in displacements for r in delta])==7
 generators=[(tuple(ps[i]),tuple(vs[i])) for i in seeds]
 full,projected=schreier_structure(generators)
 assert full['projection_order']==2520 and full['kernel_dimension']==6
 assert full['image_order']==161280
 pair_words=case['ordinary_pair_words_on_seed_indices']
 assert len(pair_words)==2
 pair,_=schreier_structure([affine_word(word,generators) for word in pair_words])
 for key in ('projection_order','kernel_dimension','image_order'):assert pair[key]==full[key]
 # Explicit A7 conjugation class, without a class-table import.
 sigma=tuple(ps[0]);conjugates=set()
 for p in projected:
  inverse=tuple(p.index(i) for i in range(7))
  conjugates.add(compose(compose(p,sigma),inverse))
 assert len(conjugates)==105 and all(tuple(p) in conjugates for p in ps)
 cover=None
 if cover_homology:
  cover=full_edge_chain(ps,crossings)
  assert cover['H1_F2_dimension']==case['expected_H1_F2_dimension']==1
 return {'name':case['name'],**full,'ordinary_pair_projection_order':pair['projection_order'],
         'ordinary_pair_kernel_dimension':pair['kernel_dimension'],
         'ordinary_pair_words_on_seed_indices':pair_words,
         'ordinary_generating_rank':2,'meridian_class_generating_rank':4,
         'meridian_class_size':105*2**2,'representation_dimension':7,
         'meridian_displacement_rank':2,'common_fixed_dimension':0,
         'bridge_upper':4,'meridional_lower':4,'optional_full_cover_F2':cover}

def check_external_census(data):
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
 return len(data['cases'])

def main():
 parser=argparse.ArgumentParser()
 parser.add_argument('--census',action='store_true')
 parser.add_argument('--cover-homology',action='store_true')
 args=parser.parse_args()
 path=HERE/'a7_binary_inputs.json';data=json.loads(path.read_text())
 assert data['field']==2 and data['dimension']==7
 assert [c['name'] for c in data['cases']]==NAMES
 for relative,digest in data['source_sha256'].items():
  assert hashlib.sha256((HERE/relative).read_bytes()).hexdigest()==digest
 results=[check_case(case,args.cover_homology) for case in data['cases']]
 corruptions=[]
 for kind in ('identity_matrix','seed_trace','census_map','ordinary_pair'):
  bad=deepcopy(data['cases'][0])
  if kind=='identity_matrix':bad['homogeneous_W_matrices_F2'][0]=[[int(i==j) for j in range(7)] for i in range(7)]
  elif kind=='seed_trace':bad['seed_trace'][0][2]=bad['seed_arcs'][0]
  elif kind=='census_map':bad['census_map_witness']['dart_images'][0]=bad['census_map_witness']['dart_images'][1]
  else:bad['ordinary_pair_words_on_seed_indices']=[[0],[0]]
  try:check_case(bad)
  except (AssertionError,KeyError,ValueError):corruptions.append(kind)
  else:raise AssertionError('corrupted certificate accepted: '+kind)
 census=check_external_census(data) if args.census else None
 files=[Path(__file__),path,HERE/'A7_BINARY_SHARP_FOUR.tex']+[HERE/p for p in data['source_sha256']]
 receipt={'status':'PASS','scope':'Three retained named knots have mr=b=4 by direct binary affine matrices. The common full finite image is W_F2 semidirect A7 of order161280, ordinary rank2 and split meridian class rank4, class size420. No branched-cover topology classification, search completeness, formal verification, external review or global priority is inferred.',
          'field':2,'results':results,'optional_named_census_reread_count':census,
          'optional_full_cover_homology_checked':args.cover_homology,
          'corruptions_rejected':corruptions,
          'input_sha256':{p.relative_to(HERE).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in files}}
 (HERE/'a7_binary_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
 print(json.dumps({'status':'PASS','sharp_four':3,'group_order':161280,'class_size':420,
                   'ordinary_rank':2,'class_rank':4,'census_rereads':census,
                   'full_cover_homology_checked':args.cover_homology,'corruptions_rejected':len(corruptions)},indent=2))

if __name__=='__main__':main()

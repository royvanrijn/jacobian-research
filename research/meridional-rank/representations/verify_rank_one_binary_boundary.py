#!/usr/bin/env python3
"""Complete binary rank-one boundary on293 named four-meridian diagrams.

Default: independently normalize the four seed covectors to the coordinate
functionals and exhaust exactly7^4=2401 models per case. Standard library only.
--full-class-control also replays the independent105-class/192-centralizer scan.
--census re-reads the pinned independently packaged named census.
"""
import argparse
from copy import deepcopy
import hashlib
from itertools import product
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
INPUT=HERE/'rank_one_binary_boundary_inputs.json'
IDENTITY=(1,2,4,8)

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


def rank_rows(rows):
 basis={}
 for row in rows:
  while row:
   k=row.bit_length()-1
   if k in basis:row^=basis[k]
   else:basis[k]=row;break
 return len(basis)

def matrix_product(a,b):
 answer=[]
 for row in a:
  value=0
  for j in range(4):
   if row>>j&1:value^=b[j]
  answer.append(value)
 return tuple(answer)

def matrix_inverse(a):
 rows=[row|(1<<(i+4)) for i,row in enumerate(a)]
 for j in range(4):
  k=next(k for k in range(j,4) if rows[k]>>j&1)
  rows[j],rows[k]=rows[k],rows[j]
  for k in range(4):
   if k!=j and rows[k]>>j&1:rows[k]^=rows[j]
 assert tuple(row&15 for row in rows)==IDENTITY
 return tuple(row>>4 for row in rows)

def class_model(data):
 pairs=[(u,f) for u in range(1,16) for f in range(1,16) if (u&f).bit_count()%2==0]
 assert len(pairs)==105 and pairs==list(map(tuple,data['class_pairs_u_phi']))
 matrices=[tuple((1<<i)^(f if u>>i&1 else 0) for i in range(4)) for u,f in pairs]
 assert matrices==list(map(tuple,data['class_matrices_row_masks']))
 lookup={pair:i for i,pair in enumerate(pairs)}
 table=[]
 for i,(u,f) in enumerate(pairs):
  assert matrix_product(matrices[i],matrices[i])==IDENTITY
  assert rank_rows([a^b for a,b in zip(matrices[i],IDENTITY)])==1
  row=[]
  for j,(v,g) in enumerate(pairs):
   root=v^(u if (f&v).bit_count()%2 else 0)
   form=g^(f if (g&u).bit_count()%2 else 0)
   k=lookup[root,form]
   assert matrix_product(matrix_product(matrices[i],matrices[j]),matrices[i])==matrices[k]
   row.append(k)
  table.append(row)
 assert table==data['conjugation_table']
 choices=[[lookup[u,1<<i] for u in range(1,16) if not(u>>i&1)] for i in range(4)]
 assert list(map(len,choices))==[7]*4
 models=list(product(*choices));assert len(models)==data['normalized_covector_seed_count']==2401
 roots=data['positive_control_seed_roots']
 control=tuple(lookup[u,1<<i] for i,u in enumerate(roots))
 assert control in models and rank_rows([pairs[x][1] for x in control])==4
 return pairs,matrices,table,models

def validate_case(case,expected_crossings=13,require_census=True):
 crossings,_=diagram(case['pd']);n=len(case['pd'])
 assert n==expected_crossings
 if require_census:census_bijection(case['pd'],case['census_pd'],case['census_map_witness'])
 key=lambda step:(step[0],)+tuple(sorted(step[1:]))
 known=set(case['seed_arcs']);assert len(known)==len(case['seed_arcs'])==4
 used=set();allcross=set(map(key,crossings))
 for o,a,b in case['seed_trace']:
  assert o in known and a in known and b not in known
  assert key((o,a,b)) in allcross
  known.add(b);used.add(key((o,a,b)))
 assert known==set(range(n)) and len(case['seed_trace'])==len(used)==n-4
 assert len(case['remaining_crossings'])==4
 assert set(map(key,case['remaining_crossings']))==allcross-used


def evaluate_models(case,table,models):
 seeds=case['seed_arcs'];trace=case['seed_trace'];left=case['remaining_crossings']
 valid=0;first=None
 for model in models:
  arcs=[0]*len(case['pd'])
  for a,label in zip(seeds,model):arcs[a]=label
  for o,a,b in trace:arcs[b]=table[arcs[o]][arcs[a]]
  if all(arcs[b]==table[arcs[o]][arcs[a]] for o,a,b in left):
   valid+=1
   if first is None:first={'seed_class_ids':list(model),'arc_class_ids':arcs[:]}
 return valid,first


def direct_search(case,table,models):
 valid,_=evaluate_models(case,table,models)
 assert valid==0,(case['name'],'rank-four-covector representation found')
 return {'name':case['name'],'ordered_covector_basis_models_checked':len(models),
         'valid_rank_four_covector_models':valid,'four_meridian_upper':4,
         'all_dimension_binary_rank_one_fixed_codimension_upper':3}


def positive_knot_control(data,pairs,matrices,table,models):
 control=data['positive_knot_control']
 assert control['name'] not in {c['name'] for c in data['cases']}
 assert control['spherogram_constructor']=='T(4,5)'
 validate_case(control,expected_crossings=15,require_census=False)
 count,witness=evaluate_models(control,table,models)
 assert count==control['expected_valid_rank_four_covector_models']==1
 assert witness is not None
 assert witness['seed_class_ids']==control['representative']['seed_class_ids']
 assert witness['arc_class_ids']==control['representative']['arc_class_ids']
 images=[matrices[i] for i in witness['arc_class_ids']]
 crossings,_=diagram(control['pd'])
 for t in images:
  assert matrix_product(t,t)==IDENTITY
  assert rank_rows([a^b for a,b in zip(t,IDENTITY)])==1
 for o,a,b in crossings:assert matrix_product(matrix_product(images[o],images[a]),images[o])==images[b]
 seeds=control['seed_arcs'];generators=[images[i] for i in seeds]
 assert rank_rows([r^i for t in generators for r,i in zip(t,IDENTITY)])==4
 group={IDENTITY};queue=[IDENTITY]
 for t in queue:
  for g in generators:
   v=matrix_product(t,g)
   if v not in group:group.add(v);queue.append(v)
 assert len(group)==control['expected_image_order']==120
 return {'name':'T(4,5)','outside_negative_panel':True,'crossings':15,
         'ordered_models_checked':len(models),'accepted_models':count,
         'common_fixed_codimension':4,'image_order':len(group),
         'all_literal_crossing_matrices_verified':True,**witness}



def full_class_control(data,pairs,matrices,table):
 control=data['optional_full_class_control']
 group=[a for a in product(range(16),repeat=4) if rank_rows(a)==4]
 assert len(group)==control['group_order']==20160
 centralizer=[g for g in group if matrix_product(g,matrices[0])==matrix_product(matrices[0],g)]
 assert len(centralizer)==192 and centralizer==list(map(tuple,control['centralizer_matrices_row_masks']))
 lookup={a:i for i,a in enumerate(matrices)};actions=[]
 for g in centralizer:
  inverse=matrix_inverse(g)
  action=[lookup[matrix_product(matrix_product(g,t),inverse)] for t in matrices]
  assert action[0]==0 and sorted(action)==list(range(105));actions.append(action)
 assert actions==control['centralizer_class_actions']
 # Inner conjugation sends the first transvection onto every class member.
 assert len({matrix_product(matrix_product(g,matrices[0]),matrix_inverse(g)) for g in group})==105
 tails=list(map(tuple,control['normalized_tails']));assert len(tails)==8176
 covered=bytearray(105**3)
 code=lambda t:(t[0]*105+t[1])*105+t[2]
 for t in tails:
  orbit={(a[t[0]],a[t[1]],a[t[2]]) for a in actions}
  assert min(orbit)==t
  indices=[code(v) for v in orbit]
  assert not any(covered[i] for i in indices)
  for i in indices:covered[i]=1
 assert all(covered)
 valid_total=0;rank_four=0;results=[]
 for case in data['cases']:
  valid=0;seeds=case['seed_arcs'];trace=case['seed_trace'];left=case['remaining_crossings']
  for tail in tails:
   arcs=[0]*13
   for a,label in zip(seeds,(0,)+tail):arcs[a]=label
   for o,a,b in trace:arcs[b]=table[arcs[o]][arcs[a]]
   if not all(arcs[b]==table[arcs[o]][arcs[a]] for o,a,b in left):continue
   valid+=1
   if rank_rows([pairs[arcs[a]][1] for a in seeds])==4:rank_four+=1
  assert valid==case['expected_valid_full_class_orbits']
  valid_total+=valid;results.append({'name':case['name'],'valid_inner_colourings':valid})
 assert valid_total==control['expected_valid_orbits']==440 and rank_four==0
 return {'group_order':20160,'class_size':105,'centralizer_order':192,
         'raw_first_normalized_tails':105**3,'tail_orbits':len(tails),
         'valid_inner_colouring_orbits':valid_total,'rank_four_covector_orbits':rank_four,'cases':results}


def reduction_controls():
 # q=1 must be excepted: a single2D transvection becomes identity modulo its fixed line.
 transvection=(3,2);assert rank_rows([transvection[0]^1,transvection[1]^2])==1
 assert transvection[0]&1 and transvection[1]==2
 # Trefoil/S3 control, with an added trivial line: q=2 survives on the quotient.
 a=(3,2,4);b=(1,3,4)
 assert matrix_product(matrix_product(a,b),a)==matrix_product(matrix_product(b,a),b)
 assert rank_rows([a[i]^(1<<i) for i in range(3)]+[b[i]^(1<<i) for i in range(3)])==2
 assert rank_rows([a[i]^(1<<i) for i in range(2)])==rank_rows([b[i]^(1<<i) for i in range(2)])==1
 return {'codimension_one_exception_checked':True,'codimension_two_reduction_control_checked':True}


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
 control=data['positive_knot_control']
 assert metadata.version('spherogram')==control['spherogram_version']
 assert [list(row) for row in Link(control['spherogram_constructor']).PD_code()]==control['pd']
 return len(data['cases'])


def main():
 parser=argparse.ArgumentParser();parser.add_argument('--census',action='store_true');parser.add_argument('--full-class-control',action='store_true');args=parser.parse_args()
 data=json.loads(INPUT.read_text());cases=data['cases']
 assert data['field']==2 and data['quotient_dimension']==4
 assert len(cases)==len({c['name'] for c in cases})==293
 source=data['panel_source'];snapshot=HERE/source['snapshot'];raw=snapshot.read_bytes()
 assert hashlib.sha256(raw).hexdigest()==source['sha256']
 assert hashlib.sha1(b'blob '+str(len(raw)).encode()+b'\0'+raw).hexdigest()==source['git_blob']
 assert len(json.loads(raw))==293 and set(json.loads(raw))=={c['name'] for c in cases}
 pairs,matrices,table,models=class_model(data)
 for case in cases:validate_case(case)
 results=[direct_search(case,table,models) for case in cases]
 controls=reduction_controls()
 knot_control=positive_knot_control(data,pairs,matrices,table,models)
 bad=deepcopy(cases[0]);bad['seed_trace'][0][2]=bad['seed_arcs'][0]
 try:validate_case(bad)
 except AssertionError:pass
 else:raise AssertionError('corrupted seed trace accepted')
 bad=deepcopy(cases[0]);bad['census_map_witness']['dart_images'][0]=bad['census_map_witness']['dart_images'][1]
 try:validate_case(bad)
 except AssertionError:pass
 else:raise AssertionError('corrupted census map accepted')
 full=full_class_control(data,pairs,matrices,table) if args.full_class_control else None
 census=check_external_census(data) if args.census else None
 files=[Path(__file__),INPUT,HERE/'RANK_ONE_BINARY_BOUNDARY.tex',snapshot]
 receipt={'status':'PASS','named_cases':293,'direct_models_per_case':2401,'direct_models_total':293*2401,
          'valid_four_dimensional_zero_fixed_transvection_models':0,'all_dimension_binary_rank_one_fixed_codimension_upper':3,
          'dual_displacement_span_upper':3,'irreducible_binary_rank_one_dimension_upper':3,
          'results':results,'reduction_controls':controls,'positive_knot_control':knot_control,'positive_rank_four_seed_control_verified':True,
          'corruptions_rejected':['seed_trace','census_map'],'optional_full_class_control':full,'optional_named_census_reread_count':census,
          'input_sha256':{p.relative_to(HERE).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
          'scope':'Exactly the frozen293 named four-seed diagrams. EveryF2 representation in any finite ambient dimension with meridian displacement rank1 has fixed-space codimension<=3. Main finite proof directly exhausts coordinate-covector models, including reducible images; independent full-class orbit control is optional. The written contragredient corollary additionally bounds displacement span and irreducible dimension by3. No conclusion for larger fields, displacement rank>=2, generic eigenvalues or unrestricted meridional rank.'}
 (HERE/'rank_one_binary_boundary_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
 print(json.dumps({'status':'PASS','named_cases':293,'direct_models_per_case':2401,'total_direct_models':293*2401,
                   'valid_rank_four_models':0,'all_dimension_binary_rank_one_codimension_upper':3,
                   'full_class_control':full is not None,'census_rereads':census,'positive_T4_5_models':knot_control['accepted_models']},indent=2))

if __name__=='__main__':main()

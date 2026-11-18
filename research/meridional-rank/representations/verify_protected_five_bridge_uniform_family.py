#!/usr/bin/env python3
"""Exact uniform five-bridge braid/surgery controls and optional parent geometry.

The all-parameter surgery and filling arguments are written in the companion
note. Default mode checks exact group words and literal planar wiring. --sage
reruns rigorous interval hyperbolicity of the two-cusped parent.
"""
from __future__ import annotations
import argparse
from collections import Counter,defaultdict
from copy import deepcopy
from pathlib import Path
import hashlib,importlib,json,platform,random,sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
INPUT=HERE/'protected_five_bridge_uniform_family_inputs.json'
from verify_unbounded_prime_augmentation_gap_audit import (
 compose,inverse,cycle,perm_hurwitz,free_hurwitz,word_mul,word_inv)
from verify_protected_five_bridge_prime import exact_certificate as base_certificate


def require(value,message):
 if not value:raise AssertionError(message)


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def exact_certificate(data):
 U=[1,2,-3,-3,2,-3,-3,2,1,1,2,2,2,-3]
 Q=[4]+U+[-4]+[2,3]*4
 P=[5,4,3,2,2,3,4,-5]
 require(data['prefix_U']==U and data['base_braid_Q']==Q
         and data['circle_suffix_P']==P and data['literal_braid_word']==Q+P,
         'literal augmented braid and slab binding')
 require(data['braid_strands']==6 and len(data['pd'])==len(Q+P)==32,
         'six strands and thirty-two presented crossings')
 pd=data['pd'];ends={}
 for i,row in enumerate(pd):
  require(len(row)==len(set(row))==4,'four distinct half-edge labels')
  for j,e in enumerate(row):ends.setdefault(e,[]).append(4*i+j)
 require(len(ends)==64 and all(len(x)==2 for x in ends.values()),'paired diagram edges')
 alpha={a:b for pair in ends.values()for a,b in(pair,pair[::-1])}
 unused=set(alpha);faces=0
 while unused:
  start=x=min(unused)
  while True:
   unused.remove(x);y=alpha[x];x=4*(y//4)+(y+1)%4
   if x==start:break
  faces+=1
 require(faces==34,'Euler characteristic two for the connected PD projection')
 graph={i:set()for i in range(32)}
 for a,b in alpha.items():graph[a//4].add(b//4)
 seen={0};todo=[0]
 for i in todo:
  for j in graph[i]-seen:seen.add(j);todo.append(j)
 require(len(seen)==32,'connected projection graph')
 # Wire every literal crossing to the PD, tracking the added sixth component.
 wire=[None]*6;top=[None]*6;colors=[0]*5+[1];edgecolors=defaultdict(set)
 axis_self_crossings=0
 for letter,row in zip(Q+P,pd):
  j=abs(letter)-1;a,b=colors[j:j+2]
  if a==b==1:axis_self_crossings+=1
  if letter>0:ports=[a,b,a,b];incoming=[0,3];outgoing=[1,2]
  else:ports=[b,a,b,a];incoming=[1,0];outgoing=[2,3]
  for col,port in[(j,incoming[0]),(j+1,incoming[1])]:
   if wire[col]is None:wire[col]=row[port];top[col]=row[port]
   require(wire[col]==row[port],'literal incoming wiring')
  wire[j:j+2]=[row[k]for k in outgoing];colors[j:j+2]=[b,a]
  for edge,color in zip(row,ports):edgecolors[edge].add(color)
 require(wire==top and colors==[0]*5+[1],'six closure identifications')
 require(all(len(x)==1 for x in edgecolors.values()),'consistent component labels')
 circle=sorted(e for e,x in edgecolors.items()if x=={1})
 knot=sorted(e for e,x in edgecolors.items()if x=={0})
 components=data['oriented_component_edges']
 require(components[data['circle_component']]==circle and len(circle)==8,
         'added sixth strand is the retained eight-edge circle component')
 require(components[data['knot_component']]==knot and len(knot)==56,
         'retained fifty-six-edge knot component')
 require(axis_self_crossings==0,'circle has no self crossings and zero blackboard framing')
 # Artin pure-braid factorization of the local point-pushing circle.
 def A(i):return list(range(5,i,-1))+[i,i]+[-j for j in range(i+1,6)]
 def reduce_indices(word):
  out=[]
  for j in word:
   if out and out[-1]==-j:out.pop()
   else:out.append(j)
  return out
 require(reduce_indices(A(2)+A(3)+A(4))==P,'P=A26*A36*A46')
 labels=list(range(1,7));roles=[];over=[];linking_twice=0
 for letter in P:
  j=abs(letter)-1;a,b=labels[j:j+2];roles.append([a,b]);over.append(b if letter>0 else a)
  if 6 in(a,b):linking_twice += 1 if letter>0 else -1
  labels[j:j+2]=[b,a]
 require(labels==list(range(1,7)),'circle suffix is pure')
 require(roles==[[5,6],[4,6],[3,6],[2,6],[6,2],[6,3],[6,4],[6,5]],'literal selected-strand roles')
 require(over==[6,6,6,6,2,3,4,6] and linking_twice==6,'circle over/under oval and linking three')
 require(data['linking_matrix']==[[0,3],[3,0]],'retained oriented linking matrix')
 for module,name in [('snappy','uniform5-snappy332-init.py'),
                     ('spherogram.links.links_base','uniform5-spherogram241-links-base.py'),
                     ('spherogram.links.invariants','uniform5-spherogram241-invariants.py')]:
  require(sha(HERE/'sources'/name)==data['geometry']['runtime_topology_source_hashes'][module],
          'retained topology source edition: '+module)
 # Universal free-group action; c is the circle meridian, w=a*b*d.
 free=[(j,)for j in range(1,7)];c=free[5];w=word_mul(*free[1:4])
 cout=word_mul(word_inv(word_mul(w,c)),c,w,c)
 expected=[free[0]]+[word_mul(word_inv(c),x,c)for x in free[1:4]]
 expected +=[word_mul(cout,word_inv(c),free[4],c,word_inv(cout)),cout]
 require(free_hurwitz(free,P)==expected,'universal circle slab action and preferred longitude word')
 triple=[(1,),(2,),(3,)];product=word_mul(*triple)
 require(free_hurwitz(triple,[1,2]*3)==
         [word_mul(word_inv(product),x,product)for x in triple],
         'positive full twist equals conjugation by inverse ordered product')
 # The conjugated subbraid really is delta in the protected subgroup.
 five=[(j,)for j in range(1,6)]
 require(free_hurwitz(five,[2,3,4,-3])==free_hurwitz(five,[-4,2,3,4]),
         'delta=sigma4^-1*(sigma2*sigma3)*sigma4 Artin identity')
 old=json.loads((HERE/'protected_five_bridge_prime_inputs.json').read_text())
 base_certificate(old)
 seed=[tuple(p)for p in old['seed_permutations']]
 pfx=perm_hurwitz(seed,U);C=pfx[3]
 framed=[pfx[1],pfx[2],compose(compose(C,pfx[4]),inverse(C))]
 p=compose(compose(framed[0],framed[1]),framed[2]);power=tuple(range(10));order=0
 while True:
  order+=1;power=compose(power,p)
  if power==tuple(range(10)):break
 require(order==12 and order<100,'framed alternating product has exact order twelve')
 delta=[2,3,4,-3]
 require(perm_hurwitz(pfx,delta*36)==pfx,'thirty-six delta powers fix the prefix tuple')
 require(perm_hurwitz(pfx,delta*4)==seed,'base returning alternating tuple')
 strand=list(range(5))
 for j in delta:at=abs(j)-1;strand[at],strand[at+1]=strand[at+1],strand[at]
 require(strand==[0,2,4,3,1],'delta strand permutation has order three')
 return {'literal_braid_PD_binding':True,'connected_spherical_PD':True,
         'knot_and_circle_components':[56,8],'circle_component_index':data['circle_component'],
         'circle_no_self_crossings':True,'circle_slabrules_and_pure_factorization':True,
         'circle_disk_piercings':3,'circle_linking_number':3,
         'universal_surgery_group_identity':True,'positive_full_twist_group_identity':True,
         'framed_A10_product_order':12,'delta_coloring_period':36,
         'all_k_nonnegative_family':'closure(U*delta^(4+36*k))',
         'all_k_protected_augmentation_upper':4,'all_k_mr_and_bridge':5,
         'surgery_slope_in_disk_marking':[1,'-12*k'],
         'infinite_hyperbolicity_and_distinctness_require_written_filling_and_volume_proofs':True}


def sage_certificate(data):
 for p in(ROOT/'.scratch/mrc-sage-snappy',Path('/tmp/mrc-knot-audit-venv/lib/python3.12/site-packages')):
  if p.exists():sys.path.append(str(p))
 import snappy,spherogram
 from sage.version import version as sage_version
 from sage.all import ComplexIntervalField,QQ
 from snappy.verify.hyperbolicity import check_logarithmic_gluing_equations_and_positively_oriented_tets
 g=data['geometry'];versions={'sage':sage_version,'snappy':snappy.__version__,
  'spherogram':spherogram.__version__,'python':sys.version.split()[0]}
 require(versions==g['runtime_versions'],'pinned analytic versions')
 sources={name:sha(importlib.import_module(name).__file__)for name in g['runtime_analytic_source_hashes']}
 require(sources==g['runtime_analytic_source_hashes'],'analytic Python source pins')
 topology_sources={name:sha(importlib.import_module(name).__file__)for name in g['runtime_topology_source_hashes']}
 require(topology_sources==g['runtime_topology_source_hashes'],'topology Python source pins')
 binary=Path(snappy.__file__).parent/g['runtime_binary']['file']
 require(sha(binary)==g['runtime_binary']['sha256'],'analytic kernel binary pin')
 L=spherogram.Link(braid_closure=data['literal_braid_word'])
 require([list(row)for row in L.PD_code()]==data['pd'],'fresh augmented braid gives retained PD')
 require(json.loads(json.dumps(L.DT_code()))==data['DT'],'fresh DT binding')
 freshL=spherogram.Link(data['pd']);index={C:i for i,C in enumerate(freshL.crossings)}
 components=[sorted(data['pd'][index[ce.crossing]][ce.strand_index]for ce in comp)for comp in freshL.link_components]
 require(components==data['oriented_component_edges'],'fresh PD component ordering and cusp label binding')
 random.seed(g['random_seed']);snappy.set_rand_seed(g['random_seed']);fresh=freshL.exterior()
 for _ in range(g['randomizations']):fresh.randomize()
 require(fresh.triangulation_isosig()==g['decorated_triangulation_isosig'],'fresh PD exterior triangulation binding')
 M=snappy.Manifold(g['triangulation_snappy'])
 require(M.triangulation_isosig()==g['decorated_triangulation_isosig'],'retained parent triangulation binding')
 require(M.num_tetrahedra()==g['tetrahedra']==31 and M.num_cusps()==g['cusps']==2,'thirty-one tetrahedra and two cusps')
 require(json.loads(json.dumps(M._get_tetrahedra_gluing_data()))==g['face_pairings'],'retained face pairings')
 require(Counter(tuple(map(int,row))for row in M.gluing_equations())==Counter(map(tuple,g['log_gluing_equations'])),'retained log-gluing rows up to edge ordering')
 CIF=ComplexIntervalField(g['bits_prec']);RIF=CIF.real_field()
 retained=[CIF(RIF(QQ(a),QQ(b)),RIF(QQ(c),QQ(d)))for a,b,c,d in g['shape_intervals_rational']]
 check_logarithmic_gluing_equations_and_positively_oriented_tets(M,retained)
 okay,shapes=M.verify_hyperbolicity(bits_prec=g['bits_prec']);require(okay,'fresh parent Krawczyk existence proof')
 check_logarithmic_gluing_equations_and_positively_oriented_tets(M,shapes)
 fresh_okay,_=fresh.verify_hyperbolicity(bits_prec=g['bits_prec']);require(fresh_okay,'fresh literal PD exterior interval proof')
 def ri(x):return[str(x.lower().exact_rational()),str(x.upper().exact_rational())]
 return {'runtime_versions':versions,'analytic_source_hashes':sources,'topology_source_hashes':topology_sources,'kernel_sha256':sha(binary),
 'literal_augmented_braid_PD_exterior_binding':True,'fresh_component_cusp_ordering':True,
 'retained_interval_log_equations_and_orientations':True,'fresh_Krawczyk_parent_hyperbolicity':True,
 'fresh_PD_exterior_hyperbolicity':True,'tetrahedra':31,'cusps':2,'bits_prec':g['bits_prec'],
 'decorated_triangulation_isosig':M.triangulation_isosig(),
 'fresh_shape_intervals_rational':[ri(z.real())+ri(z.imag())for z in shapes]}


def main():
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument('--sage',action='store_true');parser.add_argument('--write-receipt',type=Path)
 args=parser.parse_args();data=json.loads(INPUT.read_text())
 report={'schema_version':1,'python':platform.python_version(),'checker_sha256':sha(__file__),
 'input_sha256':sha(INPUT),'dependency_sha256':{name:sha(HERE/name)for name in(
 'verify_unbounded_prime_augmentation_gap_audit.py','verify_protected_five_bridge_prime.py',
 'verify_new_seventeen_crossing_knot.py',
 'protected_five_bridge_prime_inputs.json')},'finite_certificate':exact_certificate(data)}
 for field in('circle_suffix_P','circle_component'):
  bad=deepcopy(data)
  if field=='circle_suffix_P':bad[field][-1]*=-1
  else:bad[field]=data['knot_component']
  try:exact_certificate(bad)
  except AssertionError:pass
  else:raise AssertionError('corrupted surgery marking accepted: '+field)
 report['negative_controls_rejected']=2
 if args.sage:report['rigorous_parent_geometry']=sage_certificate(data)
 report['scope']={'formal_verification':False,'independent_replay':False,'external_review':False}
 report['result']='PASS'
 if args.write_receipt:args.write_receipt.write_text(json.dumps(report,indent=2)+'\n')
 print(json.dumps({k:v for k,v in report.items()if k!='rigorous_parent_geometry'},indent=2))
 if args.sage:print('Rigorous two-cusped parent geometry: PASS.')

if __name__=='__main__':main()
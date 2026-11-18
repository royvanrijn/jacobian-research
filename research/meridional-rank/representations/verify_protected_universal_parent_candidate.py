#!/usr/bin/env python3
"""Exact three-component universal-parent candidate controls and optional parent geometry.

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
INPUT=HERE/'protected_universal_parent_candidate_inputs.json'
from verify_unbounded_prime_augmentation_gap_audit import (
 compose,inverse,cycle,perm_hurwitz,free_hurwitz,word_mul,word_inv)



def require(value,message):
 if not value:raise AssertionError(message)


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def exact_certificate(data):
 U=[1,2,-3,-3,2,-3,-3,2,1,1,2,2,2,-3]
 R=[4,3,2,2,3,4,-5,4]+U+[-4,2,3,5];D=[6,5,4,4,5,-6];word=R+D;pd=data['pd'];N=len(word)
 require(data['prefix_U']==U and data['base_braid_R']==R and data['new_circle_suffix_D']==D
  and data['literal_braid_word']==word,'literal universal candidate braid')
 require(data['braid_strands']==7 and len(pd)==N==32,'seven strands and thirty-two crossings')
 ends={}
 for i,row in enumerate(pd):
  require(len(row)==len(set(row))==4,'distinct PD crossing ports')
  for j,e in enumerate(row):ends.setdefault(e,[]).append(4*i+j)
 require(len(ends)==2*N and all(len(x)==2 for x in ends.values()),'paired PD edges')
 alpha={a:b for pair in ends.values()for a,b in(pair,pair[::-1])};unused=set(alpha);faces=0
 while unused:
  start=x=min(unused)
  while True:
   unused.remove(x);y=alpha[x];x=4*(y//4)+(y+1)%4
   if x==start:break
  faces+=1
 require(faces==N+2,'spherical projection Euler characteristic')
 graph={i:set()for i in range(N)}
 for a,b in alpha.items():graph[a//4].add(b//4)
 reached={0};todo=[0]
 for i in todo:
  for j in graph[i]-reached:reached.add(j);todo.append(j)
 require(len(reached)==N,'connected PD projection')
 colors=[0,0,0,0,1,0,2];wire=[None]*7;top=wire[:];edgecolors=defaultdict(set)
 lk2=defaultdict(int);selfcross=defaultdict(int)
 for letter,row in zip(word,pd):
  j=abs(letter)-1;a,b=colors[j:j+2]
  if a==b:selfcross[a]+=1
  else:lk2[tuple(sorted((a,b)))]+=1 if letter>0 else-1
  if letter>0:ports=[a,b,a,b];incoming=[0,3];outgoing=[1,2]
  else:ports=[b,a,b,a];incoming=[1,0];outgoing=[2,3]
  for col,port in[(j,incoming[0]),(j+1,incoming[1])]:
   if wire[col]is None:wire[col]=row[port];top[col]=row[port]
   require(wire[col]==row[port],'literal candidate incoming wiring')
  wire[j:j+2]=[row[k]for k in outgoing];colors[j:j+2]=[b,a]
  for edge,color in zip(row,ports):edgecolors[edge].add(color)
 require(wire==top and colors==[0,0,0,0,1,0,2],'all seven closure identifications')
 require(all(len(v)==1 for v in edgecolors.values()),'consistent edge colors')
 components=data['oriented_component_edges'];ids=data['component_indices']
 for color,name,count in[(0,'knot_K',48),(1,'old_circle_C',10),(2,'new_circle_D',6)]:
  edges=sorted(e for e,v in edgecolors.items()if v=={color})
  require(edges==components[ids[name]] and len(edges)==count,'marked component '+name)
 require(selfcross[1]==selfcross[2]==0,'both circle components have no self crossings')
 require(dict(lk2)=={(0,1):6,(1,2):2,(0,2):2},'old linking three and new circle linking one with each')
 require(data['linking_matrix']==[[0,1,1],[1,0,3],[1,3,0]],'retained linking matrix')
 # D=A47*A57, so it is a local disk circle around positions4,5.
 free=[(i,)for i in range(1,8)];a,b,c,e=free[3],free[4],free[6],free[5]
 w=word_mul(a,b);cout=word_mul(word_inv(word_mul(w,c)),c,w,c)
 expected=free[:3]+[word_mul(word_inv(c),a,c),word_mul(word_inv(c),b,c),
  word_mul(cout,word_inv(c),e,c,word_inv(cout)),cout]
 require(free_hurwitz(free,D)==expected,'two-punctured D-disk free-group action')
 require(data['disk_selected_positions']==[4,5] and data['D_filling_slope']==[1,'-q'],
  'geometric disk slope convention')
 for module,name in [('snappy','uniform5-snappy332-init.py'),('spherogram.links.links_base','uniform5-spherogram241-links-base.py'),('spherogram.links.invariants','uniform5-spherogram241-invariants.py')]:
  require(sha(HERE/'sources'/name)==data['geometry']['runtime_topology_source_hashes'][module],'topology source pin')
 require(data['uniform_identification_certified_by_this_geometry_checker'] is False,'geometry checker coverage stays explicit')
 return {'literal_braid_PD_binding':True,'connected_spherical_PD':True,
  'component_edge_counts':{'D':6,'C':10,'K':48},'new_circle_D_links_C_and_K_once':True,
  'D_disk_selected_strands':[4,5],'D_disk_surgery_identity':True,
  'proved_D_filling_braid':'closure(R*sigma4^(2*q))',
  'uniform_M_b_identification_certified_by_this_geometry_checker':False}

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
 require(M.num_tetrahedra()==g['tetrahedra']==35 and M.num_cusps()==g['cusps']==3,'thirty-five tetrahedra and three cusps')
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
 'fresh_PD_exterior_hyperbolicity':True,'tetrahedra':35,'cusps':3,'bits_prec':g['bits_prec'],
 'decorated_triangulation_isosig':M.triangulation_isosig(),
 'fresh_shape_intervals_rational':[ri(z.real())+ri(z.imag())for z in shapes]}


def main():
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument('--sage',action='store_true');parser.add_argument('--write-receipt',type=Path)
 args=parser.parse_args();data=json.loads(INPUT.read_text())
 report={'schema_version':1,'python':platform.python_version(),'checker_sha256':sha(__file__),
 'input_sha256':sha(INPUT),'dependency_sha256':{name:sha(HERE/name)for name in(
 'verify_unbounded_prime_augmentation_gap_audit.py',)},'finite_certificate':exact_certificate(data)}
 for field in('new_circle_suffix_D','component_indices'):
  bad=deepcopy(data)
  if field=='new_circle_suffix_D':bad[field][-1]*=-1
  else:bad[field]['new_circle_D']=data['component_indices']['old_circle_C']
  try:exact_certificate(bad)
  except AssertionError:pass
  else:raise AssertionError('corrupted surgery marking accepted: '+field)
 report['negative_controls_rejected']=2
 if args.sage:report['rigorous_parent_geometry']=sage_certificate(data)
 report['scope']={'formal_verification':False,'independent_replay':False,'external_review':False}
 report['result']='PASS'
 if args.write_receipt:args.write_receipt.write_text(json.dumps(report,indent=2)+'\n')
 print(json.dumps({k:v for k,v in report.items()if k!='rigorous_parent_geometry'},indent=2))
 if args.sage:print('Rigorous three-cusped universal candidate geometry: PASS.')

if __name__=='__main__':main()
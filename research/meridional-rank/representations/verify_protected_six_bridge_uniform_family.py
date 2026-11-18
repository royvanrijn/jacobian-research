#!/usr/bin/env python3
"""Exact all-bridge coloring controls and six-bridge surgery controls and optional parent geometry.

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
INPUT=HERE/'protected_six_bridge_uniform_family_inputs.json'
from verify_unbounded_prime_augmentation_gap_audit import (
 compose,inverse,cycle,perm_hurwitz,free_hurwitz,word_mul,word_inv)



def require(value,message):
 if not value:raise AssertionError(message)


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def uniform_algebra_certificate():
 U=[1,2,-3,-3,2,-3,-3,2,1,1,2,2,2,-3]
 checked=[]
 for b in range(4,25):
  d=b-2;n=2*b;rho=list(range(b-1,3,-1))
  delta=[-j for j in rho[::-1]]+list(range(2,b-1))+rho
  seed=[cycle(n,(0,1),(2,3)),cycle(n,(1,4),(2,3)),cycle(n,(0,5),(1,6)),cycle(n,(1,2),(6,7))]
  if b>=5:seed.append(cycle(n,(0,8),(7,9)))
  seed +=[cycle(n,(2*j-2,2*j),(2*j-1,2*j+1))for j in range(5,b)]
  pfx=perm_hurwitz(seed,U);C=pfx[3]
  framed=[pfx[1],pfx[2]]+[compose(compose(C,z),inverse(C))for z in pfx[4:]]
  p=tuple(range(n))
  for z in framed:p=compose(p,z)
  require(p==cycle(n,tuple([0]+list(range(8,n,2))+[5]),
   tuple([1,6]+list(range(9,n,2))+[4]),(2,3)),'uniform framed product cycle formula')
  require(perm_hurwitz(seed,U+delta*(d+1))==seed,'uniform base return exact control')
  X,Y=seed[1],seed[2];a,h=pfx[1],pfx[2]
  require(compose(compose(p,X),inverse(p))==h,'first uniform return identity')
  pd=tuple(range(n))
  for _ in range(d):pd=compose(pd,p)
  require(compose(compose(pd,Y),inverse(pd))==a,'last uniform return identity')
  z=Y
  for actual in framed[2:]:
   z=compose(compose(p,z),inverse(p));require(z==actual,'uniform new generator recurrence')
  strand=list(range(b))
  for letter in U+delta*(d+1):
   i=abs(letter)-1;strand[i],strand[i+1]=strand[i+1],strand[i]
  visited=set();at=0
  while at not in visited:visited.add(at);at=strand[at]
  require(len(visited)==b and at==0,'uniform base closure is a knot')
  checked.append(b)
 return {'bounded_formula_controls_b':checked,'uniform_written_seed_and_period_theorem':True,
  'framed_product_cycle_lengths':['d','d+1',2,1],
  'framed_product_order':'d*(d+1)','coloring_period':'d*d*(d+1)',
  'family_exponent':'(d+1)*(1+d*d*k)',
  'all_b_all_k_gap':'ar_C<=b-1<mr=br=b',
  'unrestricted_primality_is_not_inferred_from_bounded_controls':True}


def exact_certificate(data):
 b=6;d=b-2;U=[1,2,-3,-3,2,-3,-3,2,1,1,2,2,2,-3]
 rho=list(range(b-1,3,-1));Q=rho+U+[-j for j in rho[::-1]]+list(range(2,b-1))*(d+1)
 P=list(range(b,1,-1))+[2]+list(range(3,b))+[-b];word=Q+P;pd=data['pd'];N=len(word)
 require(data['bridge_parameter_b']==b and data['prefix_U']==U and data['base_braid_Q']==Q
  and data['circle_suffix_P']==P and data['literal_braid_word']==word,'literal six-bridge parent binding')
 require(data['braid_strands']==b+1 and len(pd)==N==43,'seven strands and forty-three crossings')
 ends={}
 for i,row in enumerate(pd):
  require(len(row)==len(set(row))==4,'PD crossing labels')
  for j,e in enumerate(row):ends.setdefault(e,[]).append(4*i+j)
 require(len(ends)==2*N and all(len(x)==2 for x in ends.values()),'paired PD half edges')
 alpha={a:c for pair in ends.values()for a,c in(pair,pair[::-1])};unused=set(alpha);faces=0
 while unused:
  start=x=min(unused)
  while True:
   unused.remove(x);y=alpha[x];x=4*(y//4)+(y+1)%4
   if x==start:break
  faces+=1
 require(faces==N+2,'spherical PD Euler characteristic')
 graph={i:set()for i in range(N)}
 for a,c in alpha.items():graph[a//4].add(c//4)
 seen={0};todo=[0]
 for i in todo:
  for j in graph[i]-seen:seen.add(j);todo.append(j)
 require(len(seen)==N,'connected projection')
 wire=[None]*(b+1);top=wire[:];colors=[0]*b+[1];edgecolors=defaultdict(set);self_crossings=0
 for letter,row in zip(word,pd):
  j=abs(letter)-1;a,c=colors[j:j+2]
  if a==c==1:self_crossings+=1
  if letter>0:ports=[a,c,a,c];incoming=[0,3];outgoing=[1,2]
  else:ports=[c,a,c,a];incoming=[1,0];outgoing=[2,3]
  for col,port in[(j,incoming[0]),(j+1,incoming[1])]:
   if wire[col]is None:wire[col]=row[port];top[col]=row[port]
   require(wire[col]==row[port],'literal incoming braid wiring')
  wire[j:j+2]=[row[k]for k in outgoing];colors[j:j+2]=[c,a]
  for edge,color in zip(row,ports):edgecolors[edge].add(color)
 require(wire==top and colors==[0]*b+[1],'all seven closure identifications')
 require(all(len(x)==1 for x in edgecolors.values()),'consistent edge components')
 circle=sorted(e for e,x in edgecolors.items()if x=={1});knot=sorted(e for e,x in edgecolors.items()if x=={0})
 require(data['oriented_component_edges'][data['circle_component']]==circle and len(circle)==10,'circle cusp label')
 require(data['oriented_component_edges'][data['knot_component']]==knot and len(knot)==76,'knot cusp label')
 require(self_crossings==0,'unknotted zero-self-crossing circle')
 # Literal pure point-pushing suffix and its geometric oval roles.
 def A(i):return list(range(b,i,-1))+[i,i]+[-j for j in range(i+1,b+1)]
 reduced=[]
 for letter in sum((A(i)for i in range(2,b)),[]):
  if reduced and reduced[-1]==-letter:reduced.pop()
  else:reduced.append(letter)
 require(reduced==P,'circle point-pushing product')
 labels=list(range(1,b+2));roles=[];over=[];lk2=0
 for letter in P:
  j=abs(letter)-1;a,c=labels[j:j+2];roles.append([a,c]);over.append(c if letter>0 else a)
  if b+1 in(a,c):lk2+=1 if letter>0 else -1
  labels[j:j+2]=[c,a]
 expected_roles=[[j,b+1]for j in range(b,1,-1)]+[[b+1,j]for j in range(2,b+1)]
 require(roles==expected_roles and over==[b+1]*(b-1)+list(range(2,b))+[b+1],
  'circle detour and over/under oval')
 require(labels==list(range(1,b+2))and lk2==2*d and data['linking_matrix']==[[0,d],[d,0]],'circle linking four')
 free=[(j,)for j in range(1,b+2)];c=free[-1];w=word_mul(*free[1:b-1]);e=free[b-1]
 cout=word_mul(word_inv(word_mul(w,c)),c,w,c)
 expected=[free[0]]+[word_mul(word_inv(c),x,c)for x in free[1:b-1]]
 expected +=[word_mul(cout,word_inv(c),e,c,word_inv(cout)),cout]
 require(free_hurwitz(free,P)==expected,'four-punctured disk slab free-group identity')
 triple=[(j,)for j in range(1,d+1)];product=word_mul(*triple)
 require(free_hurwitz(triple,list(range(1,d))*d)==
  [word_mul(word_inv(product),x,product)for x in triple],'positive four-strand full twist identity')
 for module,name in [('snappy','uniform5-snappy332-init.py'),('spherogram.links.links_base','uniform5-spherogram241-links-base.py'),('spherogram.links.invariants','uniform5-spherogram241-invariants.py')]:
  require(sha(HERE/'sources'/name)==data['geometry']['runtime_topology_source_hashes'][module],'retained topology source pin')
 return {'uniform_algebra':uniform_algebra_certificate(),'literal_braid_PD_binding':True,
  'connected_spherical_PD':True,'circle_component_index':data['circle_component'],
  'circle_disk_piercings':4,'circle_linking_number':4,'circle_self_crossings':0,
  'exact_point_pushing_and_surgery_identities':True,'six_bridge_exponent':'5+80*k',
  'disk_framed_circle_slope':[1,'-20*k'],'all_k_augmentation_upper':5,
  'all_k_meridional_rank_and_bridge':6,
  'prime_tail_and_infinite_distinctness_require_written_filling_argument':True}

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
 require(M.num_tetrahedra()==g['tetrahedra']==34 and M.num_cusps()==g['cusps']==2,'thirty-four tetrahedra and two cusps')
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
 'fresh_PD_exterior_hyperbolicity':True,'tetrahedra':34,'cusps':2,'bits_prec':g['bits_prec'],
 'decorated_triangulation_isosig':M.triangulation_isosig(),
 'fresh_shape_intervals_rational':[ri(z.real())+ri(z.imag())for z in shapes]}


def main():
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument('--sage',action='store_true');parser.add_argument('--write-receipt',type=Path)
 args=parser.parse_args();data=json.loads(INPUT.read_text())
 report={'schema_version':1,'python':platform.python_version(),'checker_sha256':sha(__file__),
 'input_sha256':sha(INPUT),'dependency_sha256':{name:sha(HERE/name)for name in(
 'verify_unbounded_prime_augmentation_gap_audit.py',)},'finite_certificate':exact_certificate(data)}
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
 if args.sage:print('Rigorous six-bridge two-cusped parent geometry: PASS.')

if __name__=='__main__':main()
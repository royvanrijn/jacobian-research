#!/usr/bin/env python3
"""Exact K_* certificate; stdlib by default, optional --topology / --sage.

--topology needs SnapPy/Spherogram/Regina. --sage needs a Sage Python
with SnapPy. Each optional replay is stored against its own source edition.
Default replay does not independently certify analytic interval assertions.
"""
import argparse
from collections import defaultdict
from copy import deepcopy
from fractions import Fraction
import hashlib
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
INPUT=HERE/'new_seventeen_crossing_knot_inputs.json'
EXPECTED_BETA=[1,2,-3,-3,2,-3,-3,2,1,1,2,2,2,-3,2,2,2]
ID=tuple(range(8))

def compose(a,b):return tuple(a[b[i]] for i in range(8))
def group_order(generators):
 seen={ID};queue=[ID]
 for a in queue:
  for b in generators:
   c=compose(a,b)
   if c not in seen:seen.add(c);queue.append(c)
 return len(seen)

def diagram(pd):
 ends={};adj={e:set() for row in pd for e in row}
 for i,row in enumerate(pd):
  assert len(row)==len(set(row))==4
  for j,e in enumerate(row):ends.setdefault(e,[]).append(4*i+j)
  adj[row[1]].add(row[3]);adj[row[3]].add(row[1])
 assert len(ends)==2*len(pd) and all(len(x)==2 for x in ends.values())
 alpha={a:b for pair in ends.values() for a,b in (pair,pair[::-1])}
 unused=set(alpha);faces=0
 while unused:
  start=x=min(unused)
  while True:
   unused.remove(x);y=alpha[x];x=4*(y//4)+(y+1)%4
   if x==start:break
  faces+=1
 assert faces==len(pd)+2
 incoming={};x=0;seen=set()
 while x not in seen:
  seen.add(x);i,j=divmod(x,4);incoming.setdefault(i,[]).append(j);x=alpha[4*i+(j+2)%4]
 assert x==0 and len(seen)==2*len(pd) and all(len(v)==2 for v in incoming.values())
 components=[];unused=set(adj)
 while unused:
  cc={min(unused)};todo=list(cc)
  for e in todo:
   for f in adj[e]-cc:cc.add(f);todo.append(f)
  components.append(sorted(cc));unused-=cc
 components.sort();assert len(components)==len(pd)
 arc={e:i for i,cc in enumerate(components) for e in cc};signed=[]
 for i,row in enumerate(pd):
  ui=next(j for j in incoming[i] if j%2==0);oi=next(j for j in incoming[i] if j%2==1)
  sign=1 if (oi-ui)%4==3 else -1
  signed.append((arc[row[1]],arc[row[ui]],arc[row[(ui+2)%4]],-sign))
 return signed,components,arc

def knot_certificate(data):
 case=data['case'];beta=case['literal_braid_word'];pd=case['pd']
 assert beta==EXPECTED_BETA and case['hurwitz_word']==[-k for k in beta]
 assert case['braid_strands']==4 and len(pd)==len(beta)==17
 signed,components,arc=diagram(pd)
 images=[tuple(p) for p in case['arc_permutations']];assert len(images)==17
 for p in images:
  assert sorted(p)==list(range(8)) and compose(p,p)==ID and sum(p[i]!=i for i in range(8))==4
  assert sum(p[i]>p[j] for i in range(8) for j in range(i+1,8))%2==0
 for o,a,b,e in signed:assert compose(compose(images[o],images[a]),images[o])==images[b]
 seeds=case['seed_arcs'];known=set(seeds);assert len(seeds)==len(known)==4 and len(case['seed_trace'])==13
 for o,a,b in case['seed_trace']:
  assert o in known and a in known and b not in known
  assert any(oo==o and {aa,bb}=={a,b} for oo,aa,bb,e in signed);known.add(b)
 assert known==set(range(17))
 initial=[images[i] for i in seeds];assert group_order(initial)==case['expected_image_order']==20160
 # Independent literal braid wiring and coloring: both input and output ports
 # are checked, including allfour closure identifications.
 current=initial[:];wire=[None]*4;top=[None]*4;tags=list(range(4));seededges=[None]*4;edgecolors=defaultdict(set);perm=list(range(4))
 for letter,row in zip(beta,pd):
  i=abs(letter)-1;a,b=current[i:i+2]
  if letter>0:
   c=compose(compose(b,a),b);ports=[a,b,c,b];incoming=[0,3];outgoing=[1,2];new=[b,c]
  else:
   c=compose(compose(a,b),a);ports=[b,a,c,a];incoming=[1,0];outgoing=[2,3];new=[c,a]
  for col,port in [(i,incoming[0]),(i+1,incoming[1])]:
   if wire[col] is None:wire[col]=row[port];top[col]=row[port]
   assert wire[col]==row[port]
   if tags[col] is not None:seededges[tags[col]]=row[port]
  wire[i:i+2]=[row[p] for p in outgoing];current[i:i+2]=new;tags[i]=tags[i+1]=None
  perm[i],perm[i+1]=perm[i+1],perm[i]
  for e,p in zip(row,ports):edgecolors[e].add(p)
 assert wire==top and current==initial and all(x is not None for x in seededges)
 assert [arc[e] for e in seededges]==seeds
 assert len(edgecolors)==34 and all(len(v)==1 for v in edgecolors.values())
 for i,cc in enumerate(components):assert {next(iter(edgecolors[e])) for e in cc}=={images[i]}
 x=0;orbit=[]
 while x not in orbit:orbit.append(x);x=perm[x]
 assert len(orbit)==4 and x==0
 return {'literal_braid_and_PD_wiring':True,'one_component':True,'presented_crossings':17,
         'literal_crossing_relations':17,'image_order':20160,'meridian_class':'double transpositions',
         'seed_coloring_steps':13,'four_cycle_braid_permutation':perm,'meridional_lower':4,'bridge_upper_via_BKVV':4}

def determinant(a):
 a=[row[:] for row in a];n=len(a);sgn=1;previous=1
 for k in range(n-1):
  pivot=next((i for i in range(k,n) if a[i][k]),None)
  if pivot is None:return 0
  if pivot!=k:a[k],a[pivot]=a[pivot],a[k];sgn=-sgn
  value=a[k][k]
  for i in range(k+1,n):
   for j in range(k+1,n):
    numerator=value*a[i][j]-a[i][k]*a[k][j];assert numerator%previous==0;a[i][j]=numerator//previous
   a[i][k]=0
  previous=value
 return sgn*a[-1][-1]

def interpolate_integer_values(values):
 # Newtonforward-difference interpolation in thebasis binomial(t,k).
 differences=[];row=list(values)
 while row:
  differences.append(row[0]);row=[b-a for a,b in zip(row,row[1:])]
 answer=[Fraction(0)]*len(values);basis=[Fraction(1)]
 for k,value in enumerate(differences):
  for i,x in enumerate(basis):answer[i]+=value*x
  nextbasis=[Fraction(0)]*(len(basis)+1)
  for i,x in enumerate(basis):nextbasis[i]-=k*x;nextbasis[i+1]+=x
  basis=[x/Fraction(k+1) for x in nextbasis]
 assert all(x.denominator==1 for x in answer)
 return [int(x) for x in answer]

def alexander_certificate(data):
 case=data['case'];signed,_,_=diagram(case['pd']);n=len(signed);a0=[[0]*n for _ in range(n)];a1=[[0]*n for _ in range(n)]
 for i,(o,a,b,e) in enumerate(signed):
  a0[i][o]+=1;a1[i][o]-=1
  if e==1:a1[i][a]+=1;a0[i][b]-=1
  else:a0[i][a]-=1;a1[i][b]+=1
 values=[determinant([[a0[i][j]+t*a1[i][j] for j in range(n-1)] for i in range(n-1)]) for t in range(n)]
 poly=interpolate_integer_values(values)
 while poly and poly[-1]==0:poly.pop()
 while poly and poly[0]==0:poly.pop(0)
 assert poly
 if poly[-1]<0:poly=[-x for x in poly]
 assert poly==poly[::-1]==case['alexander_coefficients_ascending']
 assert abs(sum(poly))==1
 det=abs(sum(c*(-1)**i for i,c in enumerate(poly)))
 assert len(poly)-1==case['alexander_breadth']==14 and det==case['alexander_determinant']==117
 assert case['expected_genus']==7 and case['expected_bridge_number']==case['expected_braid_index']==4 and case['expected_crossing_number']==17
 return {'method':'Literal orientedPD Foxminor; exact17point integerdeterminants/interpolation, degreebound16',
         'coefficients_ascending':poly,'breadth':14,'determinant':det,'genus_lower_via_Seifert_matrix':7,
         'canonical_four_braid_genus_upper':7,'crossing_number_via_Yamada':17}

def interval(value):
 assert len(value)==2;a,b=map(Fraction,value);assert a<=b;return a,b

def analytic_certificate_structure(data):
 c=data['closed_representation_certificate'];h=data['cusped_hyperbolicity_certificate'];n=len(c['shape_intervals_rational']);assert n==52 and len(h['shape_intervals_rational'])==26
 negative=0
 for box in c['shape_intervals_rational']:
  assert len(box)==4;interval(box[:2]);lo,hi=interval(box[2:]);assert hi<0 or lo>0;negative+=hi<0
 assert negative==c['negative_oriented_shape_intervals']==14
 for box in h['shape_intervals_rational']:interval(box[:2]);assert interval(box[2:])[0]>0
 hnf=c['exact_edge_equation_HNF'];full=hnf['full_rectangular_equations'];selected=hnf['selected_rectangular_equations'];columns=hnf['unimodular_transform_columns']
 assert len(full)==n+1 and len(selected)==n and len(columns)==n and all(len(col)==n for col in columns)
 udet=determinant([list(row) for row in zip(*columns)]);assert abs(udet)==1 and udet==hnf['determinant']
 augmented=[]
 for a,b,sgn in full[:n]:assert len(a)==len(b)==n and sgn in (-1,1);augmented.append(a+b+[(1-sgn)//2])
 transformed=[[sum(augmented[k][j]*col[k] for k in range(n)) for j in range(2*n+1)] for col in columns]
 assert transformed[0]==hnf['discarded_identity_equation'] and all(x==0 for x in transformed[0][:-1]) and transformed[0][-1]%2==0
 derived=[(r[:n],r[n:2*n],(-1)**r[-1]) for r in transformed[1:]]
 assert [[a,b,c] for a,b,c in derived]+[full[-1]]==selected
 logs=c['log_gluing_rows'];assert len(logs)==n+1 and all(len(row)==3*n for row in logs)
 for row,(a,b,sgn) in zip(logs,full):
  assert a==[row[3*i]-row[3*i+2] for i in range(n)]
  assert b==[row[3*i+2]-row[3*i+1] for i in range(n)]
  assert sgn==(-1)**sum(row[3*i+2] for i in range(n))
 # The raw rectangular sign and the z'' substitution sign cancel.
 # The ACTUAL filled peripheral multiplier has RHS +1.
 assert full[-1][2]*((-1)**sum(logs[-1][3*i+2] for i in range(n)))==1
 assert len(c['log_equation_integer_lifts'])==n+1 and c['log_equation_integer_lifts'][-1]==-2
 lo,hi=interval(c['longitude_log_interval'][:2]);assert hi<0 or lo>0
 lo,hi=interval(c['signed_Bloch_Wigner_sum_rational_interval']);assert Fraction(1037,100)<lo<=hi<Fraction(1039,100)
 top=data['actual_double_cover'];assert top['cover_info']['degree']==2 and top['cover_info']['type']=='cyclic' and top['lifted_filling']==[[1.0,0.0]] and top['closed_tetrahedra']==25 and top['closed_H1']=='Z_117'
 return {'exact_unimodular_edge_equation_HNF':True,'full_edge_equations':52,'independent_edge_equations':51,
         'rational_cusped_shape_boxes':26,'rational_closed_representation_shape_boxes':52,
         'negative_closed_shapes':14,'nonzero_longitude_real_log_recorded':True,
         'raw_filling_rectangular_constant':full[-1][2],'actual_filling_multiplier_RHS':1,
         'signed_volume_rational_interval_positive':True,'analytic_existence_and_transcendentals_rechecked_by_default':False}


def rank_mod2(rows):
 basis={}
 for row in rows:
  while row:
   i=row.bit_length()-1
   if i in basis:row^=basis[i]
   else:basis[i]=row;break
 return len(basis)

def mod2_kernel_certificate(data):
 k=data['actual_double_cover']['mod2_orbifold_kernel'];gens=k['generators'];n=len(gens);phi=k['generator_images_mod2']
 assert n==len(set(gens))==len(phi)==4 and all(x in (0,1) for x in phi)
 def mask(word):
  answer=0
  for c in word:answer^=1<<gens.index(c.lower())
  return answer
 rows=[mask(w) for w in k['relators']];assert rank_mod2(rows)==k['relator_matrix_mod2_rank']==n-1
 value=sum(b<<i for i,b in enumerate(phi))
 assert all((row&value).bit_count()%2==0 for row in rows)
 assert (mask(k['meridian_word'])&value).bit_count()%2==1
 assert (mask(k['longitude_word'])&value).bit_count()%2==0
 assert rank_mod2(rows+[mask(k['meridian_word'])])==n
 assert k['boundary_cover_lattice_columns_in_base_meridian_longitude']==[[2,0],[0,1]]
 return {'unique_meridional_mod2_character':phi,'orbifold_generators':n,'mod2_relator_rank':n-1,
         'meridian_character':1,'longitude_character':0,'primitive_lifted_meridian_squared':True}


def trefoil_sum_control(data):
 # Completing meridian squares gives S3; its meridional-parity kernel is C3.
 a=(1,0,2);b=(0,2,1);identity=(0,1,2)
 def cp(p,q):return tuple(p[q[i]] for i in range(3))
 group={identity};queue=[identity]
 for p in queue:
  for q in (a,b):
   r=cp(p,q)
   if r not in group:group.add(r);queue.append(r)
 assert len(group)==6 and cp(cp(cp(a,b),cp(a,b)),cp(a,b))==identity
 kernel=[p for p in group if sum(p[i]>p[j] for i in range(3) for j in range(i+1,3))%2==0]
 assert len(kernel)==3 and all(cp(cp(p,p),p)==identity for p in kernel)
 sigma=data['case']['arc_permutations'][data['case']['seed_arcs'][0]];unseen=set(range(8));cycles=[]
 while unseen:
  start=x=min(unseen);cycle=[]
  while True:
   unseen.remove(x);cycle.append(x);x=sigma[x]
   if x==start:break
  cycles.append(cycle)
 assert sorted(map(len,cycles))==[1,1,1,1,2,2]
 return {'seed_degree':8,'meridian_cycle_count':6,'index_delta':2,
         'trefoil_double_cover_group_order':3,'F3_H1_per_double_component':1,
         'per_trefoil_added_F3_H1':2,'sphere_gluing_graph_vertices_per_step':7,
         'sphere_gluing_graph_edges_per_step':6,'all_s_lower_formula':'ceil((7+h0+2s)/2)>=4+s, h0>=0'}

def check_backend_sources(data,snappy):
 root=Path(snappy.__file__).parent
 for name,pin in data['backend_source_sha256'].items():assert hashlib.sha256((root/name).read_bytes()).hexdigest()==pin
 return dict(data['backend_source_sha256'])

def topology_replay(data):
 import random,warnings
 warnings.filterwarnings('ignore',message='Plink failed.*')
 import snappy,spherogram,regina
 from importlib.metadata import version
 assert version('snappy')==data['required_optional_versions']['snappy']
 assert version('spherogram')==data['required_optional_versions']['spherogram']
 assert regina.versionString()==data['required_optional_versions']['regina']
 binary=data['trusted_software_imports']['topological_cover_and_peripheral_lifting']['runtime_binary']
 assert hashlib.sha256((Path(snappy.__file__).parent/binary['file']).read_bytes()).hexdigest()==binary['sha256']
 c=data['case'];t=data['actual_double_cover'];seed=t['random_seed'];random.seed(seed);snappy.set_rand_seed(seed);regina.RandomEngine.reseedWithDefault()
 braid=spherogram.Link(braid_closure=c['literal_braid_word'])
 assert [list(row) for row in braid.PD_code()]==c['pd']
 assert json.loads(json.dumps(braid.DT_code()))==c['DT'] and len(braid.link_components)==1
 # Reset before the actual exterior construction, exactly as in the retained state.
 random.seed(seed);snappy.set_rand_seed(seed);E=spherogram.Link(c['pd']).exterior()
 assert E.num_cusps()==1 and E.triangulation_isosig()==t['exterior_isosig_with_peripheral']
 # Use the retained tetrahedron ordering for the separately pinned group words.
 # The decorated isomorphism above already identifies meridian and longitude.
 E=snappy.Manifold(t['exterior_snappy'])
 assert E.triangulation_isosig()==t['exterior_isosig_with_peripheral']
 E.dehn_fill((2,0));G=E.fundamental_group();k=t['mod2_orbifold_kernel']
 assert G.generators()==k['generators'] and G.relators()==k['relators']
 assert list(G.peripheral_curves()[0])==[k['meridian_word'],k['longitude_word']]
 covers=E.covers(2);assert len(covers)==1;C=covers[0]
 assert C.cover_info()['degree']==2 and C.cover_info()['type']=='cyclic'
 assert C.num_cusps()==1 and C.cusp_info('filling')==[(1.0,0.0)]
 assert C.triangulation_isosig()==t['cover_isosig_with_peripheral']
 explicit=E.cover([[1,0] if b else [0,1] for b in k['generator_images_mod2']])
 assert explicit.triangulation_isosig()==C.triangulation_isosig()
 assert explicit.cusp_info('filling')==[(1.0,0.0)]
 T=regina.Triangulation3(C.filled_triangulation()._to_string());T.simplify()
 assert T.isValid() and T.isClosed() and T.isConnected() and T.isOrientable()
 assert T.isoSig()==t['closed_isosig'] and T.size()==t['closed_tetrahedra'] and str(T.homology())==t['closed_H1']
 return {'literal_braid_PD_DT_freshly_rebuilt':True,'explicit_mod2_kernel_cover_matches_unique_cover':True,
         'base_filling':[2,0],'tracked_primitive_lifted_filling':[1,0],
         'actual_closed_double_cover_isosig':T.isoSig(),'closed_tetrahedra':T.size(),'closed_H1':str(T.homology()),
         'versions':{'snappy':snappy.__version__,'spherogram':version('spherogram'),'regina':regina.versionString()}}

def sage_replay(data):
 import sage.all as sage
 from sage.version import version as sageversion
 import snappy
 from snappy.verify import hyperbolicity
 from snappy.verify.volume import _volume_from_shape
 from snappy.snap.shapes import enough_gluing_equations
 assert sageversion==data['required_optional_versions']['sage'] and snappy.__version__==data['required_optional_versions']['snappy']
 binary=data['trusted_software_imports']['verified_analytic_backend']['runtime_binary']
 assert hashlib.sha256((Path(snappy.__file__).parent/binary['file']).read_bytes()).hexdigest()==binary['sha256']
 sources=check_backend_sources(data,snappy);top=data['actual_double_cover'];h=data['cusped_hyperbolicity_certificate'];c=data['closed_representation_certificate']
 M=snappy.Manifold(top['exterior_snappy']);assert M.num_cusps()==1 and M.triangulation_isosig()==top['exterior_isosig_with_peripheral']
 shapes=M.tetrahedra_shapes('rect',bits_prec=h['bits_prec'],intervals=True)
 hyperbolicity.check_logarithmic_gluing_equations_and_positively_oriented_tets(M,shapes)
 def ri(a):return [str(a.lower().exact_rational()),str(a.upper().exact_rational())]
 def ci(z):return ri(z.real())+ri(z.imag())
 assert [ci(z) for z in shapes]==h['shape_intervals_rational']
 cusped={'verified_complete_cusped_hyperbolicity':True,'bits_prec':h['bits_prec'],
         'shape_intervals_rational':[ci(z) for z in shapes],'triangulation_isosig':M.triangulation_isosig()}
 C=snappy.Manifold(top['filled_cover_snappy']);assert C.triangulation_isosig()==top['cover_isosig_with_peripheral']
 assert C.cusp_info('filling')==[(1.0,0.0)]
 assert [list(map(int,row)) for row in C.gluing_equations()]==c['log_gluing_rows']
 rect=[(list(map(int,a)),list(map(int,b)),int(sign)) for a,b,sign in C.gluing_equations('rect')]
 assert json.loads(json.dumps(rect))==c['exact_edge_equation_HNF']['full_rectangular_equations']
 assert json.loads(json.dumps(enough_gluing_equations(C)))==c['exact_edge_equation_HNF']['selected_rectangular_equations']
 shapes=C.tetrahedra_shapes('rect',bits_prec=c['bits_prec'],intervals=True)
 CIF=shapes[0].parent();RIF=CIF.real_field();two_pi=RIF(2)*RIF.pi()
 assert all(z.imag()>0 or z.imag()<0 for z in shapes)  # In particular, no z is0 or1.
 logs=[x for z in shapes for x in (z.log(),(1/(1-z)).log(),((z-1)/z).log())]
 lifts=[]
 for row in c['log_gluing_rows']:
  lhs=sum((n*z for n,z in zip(row,logs)),CIF(0));ratio=lhs.imag()/two_pi;k=sage.ZZ(ratio.center().round())
  assert abs(lhs.real())<RIF(1)/100 and abs(ratio-k)<RIF(1)/100;lifts.append(int(k))
 assert lifts==c['log_equation_integer_lifts'] and lifts[-1]==-2
 # The full shape multiplier equals +1: the two rectangular minus signs cancel.
 assert rect[-1][2]*((-1)**sum(c['log_gluing_rows'][-1][3*i+2] for i in range(len(shapes))))==1
 unfilled=C.copy();unfilled.dehn_fill((0,0));peripheral=list(unfilled.gluing_equations())[-2:]
 longitude=sum((int(n)*z for n,z in zip(peripheral[1],logs)),CIF(0))
 assert longitude.real()>0 or longitude.real()<0
 volume=sum((_volume_from_shape(z) for z in shapes),RIF(0))
 assert volume>RIF(1037)/100 and volume<RIF(1039)/100
 assert [ci(z) for z in shapes]==c['shape_intervals_rational']
 closed={'full_rectangular_shape_solution_certified':True,'bits_prec':c['bits_prec'],
         'shape_intervals_rational':[ci(z) for z in shapes],'full_log_integer_lifts':lifts,
         'actual_filling_multiplier_RHS':1,'longitude_log_interval':ci(longitude),
         'signed_representation_volume_rational_interval':ri(volume),
         'positive_signed_representation_volume':True,'closed_cover_hyperbolicity_asserted':False}
 return {'cusped':cusped,'closed_representation':closed,'backend_source_sha256':sources,
         'versions':{'sage':sageversion,'snappy':snappy.__version__}}

def rejected_mutations(data):
 rejected=[]
 for kind in ['literal_PD','arc_permutation','seed_trace','Alexander_coefficient','HNF_transform','signed_volume_sign','mod2_meridian']:
  bad=deepcopy(data)
  if kind=='literal_PD':bad['case']['pd'][0][0]=bad['case']['pd'][0][1];call=knot_certificate
  elif kind=='arc_permutation':bad['case']['arc_permutations'][0]=list(ID);call=knot_certificate
  elif kind=='seed_trace':bad['case']['seed_trace'][0][2]=bad['case']['seed_arcs'][0];call=knot_certificate
  elif kind=='Alexander_coefficient':bad['case']['alexander_coefficients_ascending'][7]+=1;call=alexander_certificate
  elif kind=='HNF_transform':bad['closed_representation_certificate']['exact_edge_equation_HNF']['unimodular_transform_columns'][0][0]+=1;call=analytic_certificate_structure
  elif kind=='signed_volume_sign':
   box=bad['closed_representation_certificate']['signed_Bloch_Wigner_sum_rational_interval'];bad['closed_representation_certificate']['signed_Bloch_Wigner_sum_rational_interval']=[str(-Fraction(box[1])),str(-Fraction(box[0]))];call=analytic_certificate_structure
  else:bad['actual_double_cover']['mod2_orbifold_kernel']['generator_images_mod2']=[0]*4;call=mod2_kernel_certificate
  try:call(bad)
  except (AssertionError,ValueError,ZeroDivisionError,KeyError):rejected.append(kind)
  else:raise AssertionError('Corruption accepted: '+kind)
 return rejected

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--topology',action='store_true');parser.add_argument('--sage',action='store_true');args=parser.parse_args();assert not(args.topology and args.sage),'Run topology and Sage modes in separate dependency environments.'
 data=json.loads(INPUT.read_text())
 for name,pin in data['source_sha256'].items():assert hashlib.sha256((HERE/name).read_bytes()).hexdigest()==pin
 finite=knot_certificate(data);alex=alexander_certificate(data);structure=analytic_certificate_structure(data);kernel=mod2_kernel_certificate(data);trefoils=trefoil_sum_control(data);corrupted=rejected_mutations(data)
 topology=topology_replay(data) if args.topology else None
 analytic=sage_replay(data) if args.sage else None
 files=[Path(__file__),INPUT,HERE/'NEW_SEVENTEEN_CROSSING_KNOT.tex'];hashes={p.relative_to(HERE).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in files}|data['source_sha256']
 receipt={'status':'PASS','finite_knot_certificate':finite,'Alexander_and_crossing_certificate':alex,
          'analytic_input_structure':structure,'mod2_cover_kernel':kernel,'trefoil_sum_control':trefoils,'corruptions_rejected':corrupted,
          'input_sha256':hashes,'optional_topology_replay':topology,'optional_Sage_replay':analytic,
          'assurance':'Default is exact stdlib finite and rational-certificate structure replay. Shape-solution existence, positivity and dilogarithm bounds are separately checked in Sage mode. Actual PD-to-double-cover binding is separately checked in topology mode. Published BKVV, Yamada, Francaviglia and LTY statements remain distinct written imports; no formal verification or external human review.'}
 path=HERE/'new_seventeen_crossing_knot_verification.json';previous=json.loads(path.read_text()) if path.is_file() else {}
 if 'internal_pdf_reviews' in previous:receipt['internal_pdf_reviews']=previous['internal_pdf_reviews']
 # Preserve explicit edition-bound runs, never converting an old backend result
 # into assurance about a later source, input or PDF edition.
 runs=previous.get('edition_replays',[])
 if args.topology or args.sage:
  record={'mode':'topology' if args.topology else 'sage','input_sha256':hashes,
          'result':topology if args.topology else analytic}
  runs=[r for r in runs if not(r['mode']==record['mode'] and r['input_sha256']==hashes)]+[record]
 if runs:receipt['edition_replays']=runs
 path.write_text(json.dumps(receipt,indent=2)+'\n')
 print(json.dumps({'status':'PASS','mr':4,'bridge':4,'genus':7,'braid_index':4,'crossing':17,
                   'default_analytic_existence_check':False,'topology_replayed':bool(topology),'Sage_replayed':bool(analytic),
                   'corruptions_rejected':corrupted},indent=2))

if __name__=='__main__':main()

#!/usr/bin/env python3
"""Exact finite controls for the all-parameter actual-cover theorem.
Default: stdlib algebra, literal PD, oriented basis chains and source pins.
--topology: fresh Spherogram/SnapPy actual cover and direct marked parent.
--sage: separate interval geometry verifier; no numerical isometry is used.
"""
from pathlib import Path
from copy import deepcopy
import argparse,json,hashlib,random,itertools,sys,time,importlib.util
HERE=Path(__file__).resolve().parent
INPUT=HERE/'all_parameter_hyperbolic_covers_inputs.json'
RECEIPT=HERE/'all_parameter_hyperbolic_covers_verification.json'
B=[1,2,-3,-3,2,-3,-3,2,1,1,2,2,2,-3,2,2,2]
P=[4,3,2,2,3,-4]
def red(w):
 out=[]
 for x in w:
  if out and out[-1]==-x:out.pop()
  else:out.append(x)
 return out
def inv(w):return [-x for x in reversed(w)]
def mul(*ws):return red(sum(ws,[]))
def hurwitz(top,word):
 cur=[a[:] for a in top]
 for e in word:
  i=abs(e)-1;a,b=cur[i:i+2]
  cur[i:i+2]=[b,mul(inv(b),a,b)] if e>0 else [mul(a,b,inv(a)),a]
 return cur
def algebra():
 assert red([4,3,2,2,-3,-4]+[4,3,3,-4])==P
 top=[[i] for i in range(1,6)];a,b,d,c=top[1:];w=mul(a,b)
 cout=mul(inv(mul(w,c)),c,w,c)
 expected=[top[0],mul(inv(c),a,c),mul(inv(c),b,c),mul(cout,inv(c),d,c,inv(cout)),cout]
 assert hurwitz(top,P)==expected
 assert hurwitz([a,b],[1,1])==[mul(inv(w),a,w),mul(inv(w),b,w)]
 assert mul(*hurwitz([a,b],[1]))==w
 # Geometric crossing roles: all crossings involve the added circle.
 roles=list(range(1,6));cross=[]
 for e in P:
  i=abs(e)-1;x,y=roles[i:i+2];assert 5 in (x,y)
  over=y if e>0 else x
  cross.append({'letter':e,'roles':[x,y],'over':over})
  roles[i],roles[i+1]=roles[i+1],roles[i]
 assert roles==list(range(1,6))
 assert cross[0]['over']==cross[-1]['over']==5
 assert [r['over']==5 for r in cross[1:-1]]==[True,True,False,False]
 return {'pure_factorization':True,'literal_free_word_Hurwitz_identity':True,'T2_conjugates_by_ab_inverse':True,'all_suffix_crossings_involve_circle':cross,'written_induction':'T^(2q) gives conjugation by(ab)^(-q); circle longitude is ab by the actual crossing disk.'}
def character(M):
 G=M.fundamental_group();gens=G.generators();rels=G.relators();per=G.peripheral_curves()
 def parity(w,phi):return sum(phi[gens.index(c.lower())] for c in w)%2
 chars=[phi for phi in itertools.product([0,1],repeat=len(gens)) if not any(parity(r,phi) for r in rels) and [parity(m,phi) for m,l in per]==[0,1]]
 assert len(chars)==1
 phi=chars[0];assert [parity(l,phi) for m,l in per]==[0,0]
 originals=[sum(phi[abs(i)-1] for i in w)%2 for w in G.original_generators(as_int_list=True)]
 return [[1,0] if v else [0,1] for v in phi],originals


def diagram(data):
 f=data['family'];assert f['base_word']==B and f['augmentation_suffix']==P
 pd=f['augmented_PD'];assert len(pd)==23
 ends={}
 for j,r in enumerate(pd):
  assert len(set(r))==4
  for h,e in enumerate(r):ends.setdefault(e,[]).append(4*j+h)
 assert set(ends)==set(range(46)) and all(len(v)==2 for v in ends.values())
 alpha={a:b for v in ends.values() for a,b in (v,v[::-1])}
 todo=set(alpha);faces=0
 while todo:
  start=x=min(todo)
  while True:
   todo.remove(x);y=alpha[x];x=4*(y//4)+(y+1)%4
   if x==start:break
  faces+=1
 assert faces==25
 # Components from actual opposite-port strands, without orientation choices.
 parent={e:e for e in ends}
 def find(e):
  while parent[e]!=e:parent[e]=parent[parent[e]];e=parent[e]
  return e
 for r in pd:
  for i in (0,1):parent[find(r[i])]=find(r[i+2])
 comps={}
 for e in ends:comps.setdefault(find(e),set()).add(e)
 groups=sorted(comps.values(),key=len);assert list(map(len,groups))==[6,40]
 circle=groups[0]
 assert all(not(circle&set(r)) for r in pd[:17])
 assert all(len(circle&set(r))==2 for r in pd[17:])
 over=[r[1] in circle for r in pd[17:]];assert over==[True,True,True,False,False,True]
 return {'literal_PD_crossings':23,'Euler_faces':faces,'circle_segment_count':6,'knot_segment_count':40,'circle_has_no_self_crossings':True,'circle_overpass_pattern':over}

def basis_static(data):
 b=data['basis_certificate'];ci=b['fixed_cusp_indices'];cur=b['fixed_curves']
 assert b['original_face_character']==[1,0,0,0,1,1,1,1,0,0,1,0,1,1]
 assert len(ci)==56 and all(len(r)==4 for r in ci) and set(x for r in ci for x in r)=={0,1,2}
 assert len(cur)==4*56 and all(len(r)==16 and all(type(x)is int for x in r) for r in cur)
 for q in (0,1):
  di,dcur=b['basis_curve_covers'][str(q)]
  assert di==ci and len(dcur)==len(cur)
  for j,(a,bb) in enumerate([[-1,2-q],[-1,2-q],[1,0]]):
   for t in range(56):
    for sheet in range(2):
     for v in range(4):
      if ci[t][v]!=j:continue
      for f in range(4):
       k=4*v+f
       assert dcur[4*t+sheet][k]==a*cur[4*t+sheet][k]+bb*cur[4*t+2+sheet][k]
 mat=[[-1,0],[2,-1]];assert b['circle_lift_matrix']==mat
 assert mat[0][0]*mat[1][1]-mat[0][1]*mat[1][0]==1
 par=data['parent'];s=[mat[0][0],mat[1][0]];fiber=[mat[0][1],mat[1][1]]
 assert par['base_circle_fill_coefficients']==s==[-1,2]
 assert par['slope_constant']==s and par['slope_n_coefficient']==[-6*x for x in fiber]==[0,6]
 assert par['base_circle_fill_cusp']==1 and par['permanently_filled_cusps']==[1,2] and par['surviving_source_cusp']==0
 assert par['fixed_fillings']==[[0,0],[0,0],[1,0]]
 return {'two_oriented_basis_chains_verified_entrywise':True,'fixed_circle_lift_matrix':mat,'determinant':1,'all_n_slope_constant':s,'all_n_slope_coefficient':[0,6],'basis_orientation_premise':b['oriented_meridian_rule'],'uniform_proof':'Degree1 torus lift linearity and relative pants parameter sum; not fitted family values.'}

def topology(data):
 def exact_fillings(M):
  out=[]
  for pair in M.cusp_info('filling'):
   vals=[int(x) for x in pair];assert all(x==y for x,y in zip(pair,vals));out.append(vals)
  return out
 import importlib.metadata
 from spherogram import Link
 import snappy
 assert snappy.__version__==data['software']['snappy']=='3.3.2'
 assert hashlib.sha256(Path(snappy.SnapPy.__file__).read_bytes()).hexdigest() in data['software']['allowed_topology_kernel_sha256']
 assert importlib.metadata.version('spherogram')==data['software']['spherogram']=='2.4.1'
 random.seed(1729);snappy.set_rand_seed(1729)
 L=Link(braid_closure=B+P);assert [list(r) for r in L.PD_code()]==data['family']['augmented_PD']
 assert L.linking_matrix()==[[0,2],[2,0]] and len(L.sublink([0]).crossings)==0
 E=L.exterior();E.dehn_fill((2,0),1);phi,orig=character(E);C=E.cover(phi)
 fixed=snappy.Manifold(str(HERE/data['parent']['fixed_cover']))
 assert exact_fillings(C)==exact_fillings(fixed)==data['parent']['fixed_fillings']==[[0,0],[0,0],[1,0]]
 assert orig==data['basis_certificate']['original_face_character']
 assert C._get_tetrahedra_gluing_data()==fixed._get_tetrahedra_gluing_data()
 assert C._get_cusp_indices_and_peripheral_curve_data()==fixed._get_cusp_indices_and_peripheral_curve_data()
 ci,cur=C._get_cusp_indices_and_peripheral_curve_data()
 assert [ci,cur]==[data['basis_certificate']['fixed_cusp_indices'],data['basis_certificate']['fixed_curves']]
 for q in (0,1):
  Eq=E.copy();Eq.dehn_fill((1,q),0);phiq,oq=character(Eq);D=Eq.cover(phiq)
  assert oq==orig and D._get_tetrahedra_gluing_data()==C._get_tetrahedra_gluing_data()
  di,dcur=D._get_cusp_indices_and_peripheral_curve_data()
  assert exact_fillings(D)==[[1,0],[1,0],[1,0]]
  assert [di,dcur]==data['basis_certificate']['basis_curve_covers'][str(q)]
 # The direct parent is generated afresh with the pinned seed and source marking.
 random.seed(1729);snappy.set_rand_seed(1729)
 F=snappy.Manifold(str(HERE/data['parent']['fixed_cover']));F.dehn_fill((-1,2),1)
 N=F.filled_triangulation([1,2]);saved=snappy.Manifold(str(HERE/data['parent']['parent']))
 assert N._get_tetrahedra_gluing_data()==saved._get_tetrahedra_gluing_data()
 assert N._get_cusp_indices_and_peripheral_curve_data()==saved._get_cusp_indices_and_peripheral_curve_data()
 assert N.triangulation_isosig()==saved.triangulation_isosig()==data['parent']['decorated_signature']
 assert N.num_tetrahedra()==14 and N.num_cusps()==1 and str(N.homology())=='Z'
 assert exact_fillings(N)==exact_fillings(saved)==[[0,0]]
 assert N.cusp_info(0)['is_complete'] and saved.cusp_info(0)['is_complete']
 return {'fresh_augmented_literal_PD':True,'unique_character_and_oriented_original_face_images':orig,'raw_fixed_marked_cover_identical':True,'fresh_and_saved_cover_exact_fillings':exact_fillings(C),'basis_curve_cover_primitive_fillings_verified':[[1,0],[1,0],[1,0]],'fresh_and_saved_parent_complete_unfilled':True,'basis_lifts_freshly_reconstructed':[0,1],'fresh_parent_raw_gluings_and_peripheral_chains_identical':True,'parent_tetrahedra':14,'parent_signature':N.triangulation_isosig(),'numerical_isometry_used':False,'separate_randomize_call_used':False,'snappy_version':snappy.__version__,'spherogram_version':importlib.metadata.version('spherogram')}

def corruptions(data):
 bads=[]
 tests=[('PD_segment',lambda d:d['family']['augmented_PD'][0].__setitem__(0,1000),diagram),('suffix_sign',lambda d:d['family']['augmentation_suffix'].__setitem__(0,-4),diagram),('oriented_lifted_chain',lambda d:d['basis_certificate']['basis_curve_covers']['1'][1][0].__setitem__(0,d['basis_certificate']['basis_curve_covers']['1'][1][0][0]+1),basis_static),('covering_matrix_sign',lambda d:d['basis_certificate']['circle_lift_matrix'][1].__setitem__(1,1),basis_static),('base_filling',lambda d:d['parent'].__setitem__('base_circle_fill_coefficients',[1,0]),basis_static),('doubling_step',lambda d:d['parent'].__setitem__('slope_n_coefficient',[0,3]),basis_static)]
 for name,mut,check in tests:
  d=deepcopy(data);mut(d)
  try:check(d)
  except (AssertionError,KeyError,IndexError):bads.append(name)
  else:raise AssertionError('corrupt control accepted: '+name)
 return bads

def save(out):
 for _ in range(5):
  oldbytes=RECEIPT.read_bytes() if RECEIPT.exists() else None
  old=json.loads(oldbytes) if oldbytes else {}
  if 'internal_pdf_reviews' in old:out['internal_pdf_reviews']=old['internal_pdf_reviews']
  if old.get('edition_sha256')==out['edition_sha256']:
   out['mode_replays']={**old.get('mode_replays',{}),**out['mode_replays']}
  tmp=RECEIPT.with_suffix('.json.tmp');tmp.write_text(json.dumps(out,indent=2)+'\n')
  if (RECEIPT.read_bytes() if RECEIPT.exists() else None)==oldbytes:tmp.replace(RECEIPT);return
 raise RuntimeError('Concurrent review update; no metadata overwritten')

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--topology',action='store_true');parser.add_argument('--sage',action='store_true');args=parser.parse_args()
 start=time.monotonic();data=json.loads(INPUT.read_text())
 for name,h in data['source_sha256'].items():assert hashlib.sha256((HERE/name).read_bytes()).hexdigest()==h,name
 finite={'algebra':algebra(),'literal_diagram':diagram(data),'oriented_basis':basis_static(data),'corruptions_rejected':corruptions(data)}
 modes={'default':finite}
 if args.topology:modes['topology']=topology(data)
 if args.sage:
  spec=importlib.util.spec_from_file_location('all_parameter_cover_geometry',HERE/data['geometry']['verifier']);mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
  modes['sage_geometry']=mod.run_verification()
 names=[INPUT.name,Path(__file__).name,'ALL_PARAMETER_HYPERBOLIC_DOUBLE_COVERS.tex']+list(data['source_sha256'])
 edition={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in dict.fromkeys(names)}
 out={'status':'PASS exact finite all-parameter cover controls; universal topology and filling proof written separately','mode_replays':modes,'this_invocation':{'topology':args.topology,'sage':args.sage},'edition_sha256':edition,'python_version':sys.version.split()[0],'elapsed_seconds':round(time.monotonic()-start,3),'scope':{'all_n_topological_identity':'Sigma2(K_n)=N(-1,2+6n)','all_n_geometry':'Separate verified area11 cusp and length>10+19n>6; Lackenby plus geometrization gives actual closed hyperbolicity','old_n0_representation_volume_premise':False,'numerical_isometry_premise':False,'formal_verification':False,'external_specialist_review':False},'written_premises':['Explicit point-push crossing disk and signed disk-twist homeomorphism','Annulus zero framing and pants relative solid-torus parameter sum','Correctness of pinned exact diagram, cover and marked permanent-filling software','Separate interval geometry certificate when Sage replayed','Primary six-theorem and geometrization statement','Independent prime-family rank/genus/crossing and trefoil consequences']}
 save(out);print(json.dumps(out,indent=2))
if __name__=='__main__':main()

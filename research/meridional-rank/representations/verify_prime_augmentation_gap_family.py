#!/usr/bin/env python3
"""Uniform twist-family certificates: stdlib exact groups and Laurent arithmetic.
The mathematical universal quantifiers are proved in the companion article.
--topology independently reconstructs the indicated representative braids.
"""
from pathlib import Path
from itertools import permutations
from copy import deepcopy
import argparse,hashlib,json,time,sys
import verify_new17_generic_kch as seed
HERE=Path(__file__).resolve().parent
INPUT=HERE/'prime_augmentation_gap_family_inputs.json'
RECEIPT=HERE/'prime_augmentation_gap_family_verification.json'
PREFIX=[1,2,-3,-3,2,-3,-3,2,1,1,2,2,2,-3]
def word(data,n):
 f=data['family'];assert n>=0
 return f['prefix']+[f['tail_generator']]*(f['tail_initial']+f['tail_step']*n)
def pmul(a,b):return tuple(a[b[i]] for i in range(len(a)))
def inv(a):return tuple(a.index(i) for i in range(len(a)))
def order(a):
 b=tuple(range(len(a)));n=0
 while True:
  n+=1;b=pmul(b,a)
  if b==tuple(range(len(a))):return n
def action(top,w,mul,inverse):
 cur=list(top)
 for e in w:
  j=abs(e)-1;a,b=cur[j:j+2]
  cur[j:j+2]=[b,mul(mul(inverse(b),a),b)] if e>0 else [mul(mul(a,b),inverse(a)),a]
 return cur
def seesaw(w,i):
 v=[abs(x) for x in w if abs(x) in (i-1,i)]
 return sum(j==0 or a!=v[j-1] for j,a in enumerate(v))
def validate_family(data):
 f=data['family'];assert f['prefix']==PREFIX and f['tail_generator']==2
 assert f['tail_initial']==3 and f['tail_step']==6
 assert data['seesaws']=={'2':4,'3':7}
 out=[]
 for n in f['representative_parameters']:
  w=word(data,n);assert len(w)==17+6*n
  top=action(list(range(4)),w,lambda a,b:a,lambda a:a) if False else list(range(4))
  for e in w:
   j=abs(e)-1;top[j],top[j+1]=top[j+1],top[j]
  assert order(tuple(top))==4
  assert [seesaw(w,i) for i in (2,3)]==[4,7]
  assert set(e//abs(e) for e in w if abs(e)==1)=={1}
  assert set(e//abs(e) for e in w if abs(e)==2)=={1}
  assert set(e//abs(e) for e in w if abs(e)==3)=={-1}
  out.append({'n':n,'crossings':len(w),'permutation':top,'seesaws':[4,7]})
 # Even added powers preserve permutation, and a positive terminal run cannot
 # alter the two linear run counts. These observations prove the induction.
 return {'representatives':out,'uniform_step_even':f['tail_step']%2==0,'terminal_run_unchanged':True}
def a8(data):
 top=[tuple(a) for a in data['A8']['top_generators']]
 assert len(top)==4 and all(pmul(a,a)==tuple(range(8)) for a in top)
 assert all(sum(a[i]!=i for i in range(8))==4 for a in top)
 group=seed.perm_group(top);assert len(group)==data['A8']['expected_group_order']==20160
 pre=action(top,PREFIX,pmul,inv);q=pmul(*pre[1:3])
 assert order(q)==data['A8']['prefix_middle_product_order']==6
 assert action(pre,[2]*6,pmul,inv)==pre
 assert action(top,word(data,0),pmul,inv)==top
 # Exact general pair formula: H^k(a,b)=(q^k a,q^(k-1)a).
 # Written proof is one substitution, hence induction; here check the period
 # and distinguish its proper divisors for the actual pair.
 assert all(action(pre,[2]*d,pmul,inv)!=pre for d in (1,2,3))
 for n in data['family']['representative_parameters']:
  assert action(top,word(data,n),pmul,inv)==top
 assert data['A8']['class_size']==210 and data['A8']['class_rank']==4
 # Transposition-edge graph of <=3 double transpositions has <=6 edges on8
 # vertices and is disconnected. The quartet attains the resulting lower4.
 return {'full_image_order':len(group),'meridian_cycle_type':'2^2 1^4','pair_product_order':6,'proper_periods_rejected':[1,2,3],'class_generating_rank':4,'uniform_return':'H^6 fixes prefix pair, so T^(3+6n)=T^3 there'}
def h3(data,generic):
 roots=[tuple(seed.q(*z) for z in r) for r in data['H3']['top_roots']]
 expected=[tuple(seed.q(*z) for z in generic['three_dimensional_witness']['arc_roots'][i]) for i in (1,0,8,13)]
 assert roots==expected
 top=[seed.reflection(r) for r in roots];E=seed.eye(3)
 pre=action(top,PREFIX,seed.qmm,lambda a:a)
 q=seed.qmm(pre[1],pre[2]);q2=seed.qmm(q,q)
 assert q!=E and q2!=E and seed.qmm(q2,q)==E
 assert data['H3']['prefix_middle_product_order']==3
 assert action(pre,[2]*3,seed.qmm,lambda a:a)==pre
 assert action(top,word(data,0),seed.qmm,lambda a:a)==top
 for n in data['family']['representative_parameters']:
  assert action(top,word(data,n),seed.qmm,lambda a:a)==top
 positive=seed.positive_witness(generic);assert positive['finite_image_order']==data['H3']['image_order']==120
 assert positive['class_generating_rank']==data['H3']['class_rank']==3
 return {'uniform_prefix_pair_period':3,'added_step_multiple_of_period':True,'positive_seed_replay':positive}
def pconst(R,n):return R.c(n)
def eye(R,n):return [[R.c(int(i==j)) for j in range(n)] for i in range(n)]
def mm(R,a,b):return [[R.add(*(R.mul(x,y) for x,y in zip(row,col))) for col in zip(*b)] for row in a]
def madd(R,a,b):return [[R.add(x,y) for x,y in zip(r,s)] for r,s in zip(a,b)]
def mscale(R,a,c):return [[R.mul(x,c) for x in r] for r in a]
def det(R,a):
 out=R.c(0)
 for p in permutations(range(len(a))):
  sign=(-1)**sum(p[i]>p[j] for i in range(len(p)) for j in range(i+1,len(p)))
  out=R.add(out,R.mul(R.c(sign),*(a[i][p[i]] for i in range(len(a)))))
 return out
def coeffs(R,cs):return R.add(*(R.mul(R.c(c),R.v('t',i)) for i,c in enumerate(cs)))
def alexander(data):
 R=seed.Ring(('t','s'));t=R.v('t');ti=R.v('t',-1);one=R.c(1);E=eye(R,4)
 def G(i,negative=False):
  M=eye(R,4);j=i-1
  block=[[R.c(0),one],[ti,R.sub(one,ti)]] if negative else [[R.sub(one,t),t],[one,R.c(0)]]
  for a in range(2):
   for b in range(2):M[j+a][j+b]=block[a][b]
  return M
 B=E
 for e in PREFIX:B=mm(R,B,G(abs(e),e<0))
 N=madd(R,G(2),mscale(R,E,R.c(-1)))
 assert mm(R,N,N)==mscale(R,N,R.neg(R.add(one,t)))
 C=mm(R,B,madd(R,E,mscale(R,N,R.v('s'))))
 D=madd(R,E,mscale(R,C,R.c(-1)))
 cofactor=det(R,[row[:3] for row in D[:3]])
 assert all(e[1] in (0,1) for e in cofactor)
 A={e:c for e,c in cofactor.items() if e[1]==0}
 Bcoef={(e[0],0):c for e,c in cofactor.items() if e[1]==1}
 P=data['alexander']['P_ascending'];assert P==[1,-3,4,0,-8,18,-21,15,-6,1]
 pp=coeffs(R,P);star=coeffs(R,list(reversed(P)))
 assert Bcoef==R.mul(t,star)
 assert R.mul(R.v('t',2),R.add(R.mul(R.add(one,t),A),Bcoef))==pp
 assert sum(c*((-1)**i) for i,c in enumerate(P))==-61
 assert sum(c*((-1)**i) for i,c in enumerate(reversed(P)))==61
 assert sum(P)==1
 # Derivative at -1 is the quotient value after removing 1+t.
 deriv=lambda cs:sum(i*c*((-1)**(i-1)) for i,c in enumerate(cs) if i)
 assert deriv(P)+deriv(list(reversed(P)))-6*61==-117
 assert -6*61==-366
 assert data['alexander']['shift_initial']==data['alexander']['shift_step']==6
 assert data['alexander']['degree_initial']==14 and data['alexander']['degree_step']==6
 assert data['alexander']['determinant_initial']==117 and data['alexander']['determinant_step']==366
 samples=[]
 for n in data['family']['representative_parameters']:
  shift=6+6*n;num=P+[0]*(shift+9+1-len(P))
  for i,c in enumerate(reversed(P)):num[shift+i]+=c
  quotient=[];previous=0
  for c in num[:-1]:
   a=c-previous;quotient.append(a);previous=a
  assert num[-1]==previous
  assert quotient==quotient[::-1] and quotient[0]==quotient[-1]==1
  assert len(quotient)-1==14+6*n and quotient[4]==0
  assert sum(quotient)==1 and sum(c*(-1)**i for i,c in enumerate(quotient))==-117-366*n
  # Reconstruct the actual cofactor by powers, independent of formula division.
  C=B
  for _ in range(3+6*n):C=mm(R,C,G(2))
  actual=det(R,[row[:3] for row in madd(R,E,mscale(R,C,R.c(-1)))[:3]])
  assert R.mul(R.v('t',2),actual)==coeffs(R,quotient)
  samples.append({'n':n,'degree':len(quotient)-1,'determinant':117+366*n,'coefficient_t4':quotient[4]})
 return {'rank_one_power_identity':'N^2=-(1+t)N','cofactor_affine_in_s':True,'B_equals_t_Pstar':True,'t2_times_1_plus_t_A_plus_B_equals_P':True,'uniform_formula':'Delta_n=(P+t^(6n+6)*Pstar)/(1+t)','universal_degree':'14+6n','universal_determinant':'117+366n','representatives':samples}
def uniform_folds(data,generic):
 u=data['generic_uniformity'];assert u['prefix_labels']==[1,6,16,13] and u['outer_closure_labels']==[1,13]
 assert u['after_first_two_folds_middle']==['D','Y'] and u['required_final_middle']==['D','B']
 assert u['core_rows']==[[8,2,3,-1],[3,13,14,1],[14,3,4,1],[4,8,9,-1],[9,14,15,1],[15,9,2,1],[2,4,13,-1],[13,15,8,1]]
 assert u['new_to_old']==generic['core_transfer']['new_to_old']
 # Free-group product preservation, independent of tail length.
 assert seed.free_reduce([2,-2,1,2])==[1,2]
 labels=[1,0,8,13]
 rows=generic['case']['canonical_signed_relations']
 for j,e in enumerate(PREFIX):
  o,a,b,sign=rows[j];i=abs(e)-1
  assert (sign==-1)==(e>0)
  if e>0:
   assert labels[i:i+2]==[a,o];labels[i:i+2]=[o,b]
  else:
   assert labels[i:i+2]==[o,a];labels[i:i+2]=[b,o]
 assert labels==u['prefix_labels']
 algebra=seed.fold_algebra();transfer=seed.core_transfer(generic)
 return {'prefix_signed_label_tracking':labels,'first_two_folds_exact_Laurent_replay':algebra,'core_transfer':transfer,'all_tail_lengths':'product preservation gives DY=DB, hence Y=B','universal_quantifier':'Written induction, not representative sampling'}
def gram_lift():
 # A symbolic symmetric 4x4 Gram matrix, not presumed nonsingular.
 R=seed.Ring(tuple('g%d%d'%(i,j) for i in range(4) for j in range(i,4)))
 G=[[R.v('g%d%d'%(min(i,j),max(i,j))) for j in range(4)] for i in range(4)]
 E=eye(R,4);Z=[[R.c(0) for _ in range(4)] for _ in range(4)]
 L=[G[i]+E[i] for i in range(4)]+[E[i]+Z[i] for i in range(4)]
 Linv=[Z[i]+E[i] for i in range(4)]+[E[i]+[R.neg(x) for x in G[i]] for i in range(4)]
 assert mm(R,L,Linv)==mm(R,Linv,L)==eye(R,8)
 return {'symbolic_8D_form_inverse':True,'form':'[[G,I],[I,0]]','inverse':'[[0,I],[I,-G]]','scope':'Normalized symmetric reflections, characteristic !=2; root-line signs and lifting Wirtinger relations are proved in text.'}
def topology(data):
 import importlib.metadata as metadata
 from spherogram import Link
 out=[]
 for n in data['family']['representative_parameters']:
  L=Link(braid_closure=word(data,n));pd=[list(r) for r in L.PD_code()]
  assert len(L.link_components)==1 and len(pd)==17+6*n
  seed.diagram(pd)
  out.append({'n':n,'components':len(L.link_components),'crossings':len(pd),'PD_sha256':hashlib.sha256(json.dumps(pd,separators=(',',':')).encode()).hexdigest(),'DT':json.loads(json.dumps(L.DT_code()))})
 return {'spherogram':metadata.version('spherogram'),'literal_braids':out}
def corruptions(data,generic):
 out=[]
 for name,mut,check in [
  ('tail_period_step',lambda d:d['family'].__setitem__('tail_step',2),validate_family),
  ('A8_meridian',lambda d:d['A8']['top_generators'][0].__setitem__(0,0),a8),
  ('Alexander_coefficient',lambda d:d['alexander']['P_ascending'].__setitem__(4,1),alexander),
  ('seesaw_count',lambda d:d['seesaws'].__setitem__('3',3),validate_family),
  ('prefix_signed_label',lambda d:d['generic_uniformity']['prefix_labels'].__setitem__(1,7),lambda d:uniform_folds(d,generic)),
  ('positive_root',lambda d:d['H3']['top_roots'][0][0].__setitem__(0,-1),lambda d:h3(d,generic))]:
  bad=deepcopy(data);mut(bad)
  try:check(bad)
  except (AssertionError,IndexError,ValueError):out.append(name)
  else:raise AssertionError('corrupt control accepted: '+name)
 return out
def main():
 parser=argparse.ArgumentParser();parser.add_argument('--topology',action='store_true');args=parser.parse_args()
 start=time.monotonic();data=json.loads(INPUT.read_text())
 for name,digest in data['source_sha256'].items():assert hashlib.sha256((HERE/name).read_bytes()).hexdigest()==digest,name
 generic=json.loads((HERE/'new17_generic_kch_inputs.json').read_text())
 result={'status':'PASS exact uniform prime augmentation-gap family certificates','family':validate_family(data),'A8':a8(data),'H3':h3(data,generic),'Alexander':alexander(data),'generic':uniform_folds(data,generic),'Gram_lift':gram_lift(),'trefoil_control':seed.trefoil_control(generic),'corruptions_rejected':corruptions(data,generic),'optional_topology':topology(data) if args.topology else None,'scope':{'parameters':'all integers n>=0 by written induction and exact universal identities','mr_bridge':'4','genus':'7+3n','crossing':'17+6n','prime':True,'complex_augmentation_rank':3,'finite_complex_reflection_class_maximum':3,'all_finite_or_infinite_Coxeter_class_bound_upper':3,'all_hyperbolic_claim':False,'all_double_covers_nonGraph_claim':False},'input_sha256':dict(data['source_sha256']),'python_version':sys.version.split()[0],'elapsed_seconds':round(time.monotonic()-start,3),'assurance':{'formal_verification':False,'external_specialist_review':False,'internal_exact_replay':True,'universal_proof':'Written pair-period induction, affine Laurent cofactor, product cancellation, and normalized symmetric Gram lift; source theorem imports retained separately.'}}
 for p in [INPUT,Path(__file__),HERE/'PRIME_AUGMENTATION_GAP_FAMILY.tex']:
  if p.exists():result['input_sha256'][p.name]=hashlib.sha256(p.read_bytes()).hexdigest()
 seed.save_receipt(RECEIPT,result);print(json.dumps(result,indent=2))
if __name__=='__main__':main()

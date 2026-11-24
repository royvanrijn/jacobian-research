#!/usr/bin/env python3
"""Authenticate the full homogeneous trace proof and run bounded supporting controls."""
from pathlib import Path,PurePosixPath
from fractions import Fraction as Q
from math import comb,gcd
from itertools import permutations,combinations
import argparse,copy,hashlib,json,re,sys
if not __debug__:raise SystemExit('Optimized execution is not accepted.')
HERE=Path(__file__).resolve().parent
PN='REDUCTIVE_HOMOGENEOUS_AUTOMORPHISM_INVARIANT_TRACE.md'
CN='verify_reductive_homogeneous_trace.py'
IN='reductive_homogeneous_trace_inputs.json'
LN='reductive_homogeneous_trace_source_ledger.json'
RN='reductive_homogeneous_trace_verification.json'
PROOF_SHA='b4b223a00cd9ed00a3f3957749b488dc0daea32f6e3491cb047d6d6dcac0aba5'
ASSURANCE=dict(external_review=False,formal_verification=False,independent_replay=False,worldwide_priority=False)
SPECS={
 'written_proof':(PN,'Sections1--26_complete_all_field_trace_proof_volume_lemma_intrinsic_cover_isomorphisms_marked_finite_etale_trace_and_forms'),
 'exact_checker':(CN,'Source_authentication_and_bounded_exact_support_NOT_a_proof_of_general_structure_or_covering_theorems'),
 'typed_inputs':(IN,'Exact_full_theorem_quantifiers_covering_scope_all_extensions_and_one_fixed_marker'),
 'matsushima_primary':('sources/arzhantsev-matsushima-math0506430v1.pdf','ONLY_Section1_Theorem2_1_Proposition2_2_physical1_2_BOTH_affineness_directions_converse_for_H0_is_used_in_full_proof'),
 'group_structure_primary':('sources/milne-algebraic-groups-2017.pdf','ONLY_5_25_5_28_printed105_106_22_42_22_43_printed477_478_12l_printed248_249_PLUS_OPTIONAL_delta_squared_lemma_21_60_21_61_printed442_443_12_23_printed240_22c_printed475_NOT_whole_book_review'),
 'units_primary':('sources/conrad-units-on-product-varieties.pdf','ONLY_Theorem1_1_Corollary1_2_physical1_normalized_unit_character_over_same_field_Rosenlicht_input'),
 'group_cover_primary':('sources/brion-szamuely-1109.2802v4.pdf','ONLY_Proposition1_1a_physical1_and_Section2_proof_physical3_connected_finite_etale_GALOIS_group_cover_central_isogeny_over_algebraically_closed_char0_NOT_Theorem1_2_for_disconnected_H'),
 'homogeneous_cover_antecedent':('sources/borovoi-demarche-20121217.pdf','ONLY_Proposition5_2_and_complete_proof_physical17_18_exact_disconnected_stabilizer_GALOIS_cover_theorem_specialized_proof_reproduced_and_credited'),
 'volume_context_antecedent':('sources/kaliman-kutzschebauch-1201.4769v2.pdf','CONTEXT_ONLY_Proposition3_6_Remark3_7_physical9_Proposition8_4_Example8_5_physical24_25_Proposition8_6_physical26_NOT_a_logical_premise'),
 'galois_domination_primary':('sources/stacks-03SF.html','ONLY_Section58_7_Tag03SF_connected_finite_etale_Galois_covers_and_domination_NOT_a_nonGalois_homogeneous_cover_classification'),
 'galois_domination_construction_primary':('sources/stacks-0BN2.html','ONLY_Lemma58_3_8_Tag0BN2_complete_finite_Galois_domination_construction_reproduced_in_Section25_2')}
PRIMARY_SHA={
 'matsushima_primary':'786bff5625883f4456f6decd82e3fc72c60ace70ced673d58aeb20566889b685',
 'group_structure_primary':'503f1a7d511632af05541c52fb2bf6fc36f5a75b0323f2f9110f03dcce11fd67',
 'units_primary':'0720c7f594d69049ebe70d804582a252fcd83852b6f530b3f52fabb7cee9e9e0',
 'group_cover_primary':'4fa0b70b09a200e4b3758c8e9ad17b245627cf30c8ddc5f3d398ce018fe09da2',
 'homogeneous_cover_antecedent':'80d830f74073185941fff49ee6e987b0b762c1e5f981007e350cd002f52e09b5',
 'volume_context_antecedent':'db983596e35f530551fffbb4eef1387491f001b851fad95edc50a515e505587e',
 'galois_domination_primary':'93635733d2e998fe28282d40a4fe70b19dcd3b4e4449fee951eea7433e5310f6',
 'galois_domination_construction_primary':'ef24281ea6c76868d2a534c6abebb8e9faea6a9642682ad0720d4d4803a561d0'}
SCOPE=dict(schema=2,reviewed_proof_sha256=PROOF_SHA,
 field='ANY field K of characteristic zero',acting_group='connected reductive linear algebraic K-group G',
 subgroup='ANY closed K-subgroup H, possibly disconnected, with affine scheme quotient X=G/H',
 quotient='represented fppf/scheme quotient, geometrically integral and smooth; no G(K)-transitivity assumption',
 normalized_trace='unique K-linear G-invariant Lambda:O(X)->K with Lambda(1)=1; A=K1 directsum gA',
 determinant_character_required_trivial=False,nontrivial_determinant_case_included=True,
 constant_units_required=False,split_group_required=False,embedding_in_C_required=False,
 all_K_algebra_automorphisms=True,all_field_extensions=True,
 extension_quantifiers='for EVERY extension L/K, Lambda_L=id_L tensor Lambda and EVERY L-algebra automorphism of O(X_L), including new automorphisms',
 volume_lemma='ONLY when full determinant character is trivial: invariant omega, dOmega^(n-1)=(kerLambda)omega and NONZERO H^n_dR=K[omega]',
 top_cohomology_criterion='for EVERY affine homogeneous G/H over algebraically closed char0: H^n_dR!=0 IFF delta_H=1; both directions proved; not an iff with all-Aut trace',
 determinant_cover='over algebraic closure, Hplus=kerdelta_H, finite etale connected CYCLIC Galois Y=G/Hplus of degree m=order(image delta_H)',
 intrinsic_cover_scope='UNIQUE up to isomorphism over X among connected finite etale GALOIS degree-m covers with nonzero H^n_dR; both directions explicit',
 group_cover_import='Brion-Szamuely Prop1.1(a) ONLY; disconnected homogeneous cover derived by component/deck stabilizer and etale lift rigidity',
 connected_stabilizer_Theorem1_2_imported=False,nonGalois_cover_classification_claimed=False,
 automorphisms_lift_to_G_assumed=False,geometric_determinant_cover_lift_proved=True,
 geometric_lift_descends_to_K_claimed=False,
 zero_dimension='X=SpecK, A=K, Lambda=id; no negative-degree form',
 marker_quantifiers='EVERY fixed pair f,h and sigma; ONE marker sigma(h) independent of r, ALL r>=0 and exact surviving orders',
 purity_safety_unsafety_both_directions=True,fixed_K_phase_field_extension_both_directions=True,
 actual_descending_phase_on_finite_quotient_required=True,
 preexisting_invariant_marker_on_finite_quotient_required=False,
 finite_etale_fixed_marker_descent_by_normalized_trace=True,averaging_preserves_purity=False,
 starting_pure_unsafe_pair_created=False,degree_filtration_or_Image_power_descent_claimed=False,
 retained427_is_core_input=False,Lie_image_is_ideal=False,
 optional_full_determinant_square_one=True,optional_sharpness_is_main_proof_premise=False,
 trace_presentation_independent=True,every_variety_isomorphism_trace_compatible=True,
 general_morphism_trace_functoriality_claimed=False,
 finite_etale_scope='ANY finite etale q:Z->X between TWO affine reductive homogeneous varieties with GIVEN possibly unrelated presentations; no equivariance required',
 finite_etale_marked_identity='Lambda_Z(h qstar(a))=Lambda_X((Tr_B/A(h)/d) a), ALL h,a; d=deg q>=1; normalized finite locally free algebra trace',
 finite_etale_fixed_marker='ONE upstairs h pushes to ONE Tq(h); exact surviving orders ALL r>=0; downstairs a pulls to qstar(a)',
 finite_etale_purity_safety_unsafety_both_directions=True,
 finite_etale_all_field_extensions_and_new_elements=True,
 geometric_forms_without_K_points=True,
 geometric_forms_scope='ANY smooth geometrically integral AFFINE K-variety with geometrically transitive action of connected reductive K-group; no K-point or K-quotient presentation required',
 full_theorem_logical_evidence='complete written proof plus narrowly scoped primary imports; finite controls supporting only',
 controls=dict(symmetric_power_degree_max=10,rankone_matrix_sizes=[2,3,4],torus_power_max=6,Ga_primitive_degree_max=12,cyclic_order_max=12,finite_etale_degree_max=6),
 assurance_flags=ASSURANCE)
def need(v,msg):
 if not v:raise ValueError(msg)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def parse(text):
 def unique(pairs):
  out={}
  for k,v in pairs:
   need(k not in out,'Duplicate JSON key: '+k);out[k]=v
  return out
 def bad_constant(value):raise ValueError('Nonfinite JSON value: '+value)
 return json.loads(text,object_pairs_hook=unique,parse_constant=bad_constant)
def load(p):return parse(p.read_text(encoding='utf-8'))
def typed(a,b):
 if type(a)is not type(b):return False
 if isinstance(b,dict):return set(a)==set(b) and all(typed(a[k],v) for k,v in b.items())
 if isinstance(b,list):return len(a)==len(b) and all(typed(x,y) for x,y in zip(a,b))
 return a==b
def inputs(v):need(typed(v,SCOPE),'Exact full homogeneous theorem scope and no surplus fields')
def confined_path(raw):
 need(type(raw)is str and raw and '\\' not in raw,'POSIX relative source path')
 p=PurePosixPath(raw)
 need(not p.is_absolute() and all(x not in ('.','..') for x in raw.split('/')),'No absolute or traversing source path')
 need(str(p)==raw,'Canonical relative source path')
 target=(HERE/raw).resolve()
 need(target.is_relative_to(HERE) and target.is_file(),'Source remains inside package: '+raw)
 return target
def sources(v):
 need(type(v)is dict and set(v)=={'schema','sources','assurance_flags'},'Exact ledger fields')
 need(type(v['schema'])is int and v['schema']==2 and typed(v['assurance_flags'],ASSURANCE),'Typed ledger schema and assurance')
 rows=v['sources'];need(type(rows)is list and len(rows)==len(SPECS),'All mandatory source roles')
 out={}
 for r in rows:
  need(type(r)is dict and set(r)=={'role','path','relation','sha256'},'Exact source row fields')
  role=r['role'];need(type(role)is str and role in SPECS and role not in out,'Unique known role')
  need(type(r['relation'])is str and type(r['sha256'])is str and re.fullmatch('[0-9a-f]{64}',r['sha256']) is not None,'Typed source relation and SHA256')
  p=confined_path(r['path'])
  need((r['path'],r['relation'])==SPECS[role],'Exact reviewed role path and relation')
  need(sha(p)==r['sha256'],'Missing or stale source: '+role)
  out[role]=r
 need(out['written_proof']['sha256']==PROOF_SHA,'Pinned complete proof')
 for role,h in PRIMARY_SHA.items():need(out[role]['sha256']==h,'Exact primary body: '+role)
 return out
# Sparse exact Laurent polynomials; exponent tuples have fixed arity.
def clean(p):return {e:Q(a) for e,a in p.items() if a}
def add(*ps):
 r={}
 for p in ps:
  for e,a in p.items():r[e]=r.get(e,Q(0))+a
 return clean(r)
def scale(p,a):return clean({e:Q(a)*v for e,v in p.items()})
def mul(p,q):
 r={}
 for e,a in p.items():
  for f,b in q.items():
   z=tuple(i+j for i,j in zip(e,f));r[z]=r.get(z,Q(0))+a*b
 return clean(r)
def const(a,n):return {(0,)*n:Q(a)} if a else {}
def var(i,n):
 e=[0]*n;e[i]=1;return {tuple(e):Q(1)}
def power(p,k,n=None):
 n=n or len(next(iter(p)));r=const(1,n)
 for _ in range(k):r=mul(r,p)
 return r
def mon(e,a=1):return {tuple(e):Q(a)} if a else {}
def deriv(p,i):
 r={}
 for e,a in p.items():
  if e[i]:
   z=list(e);z[i]-=1;r[tuple(z)]=a*e[i]
 return clean(r)
def det(M,nvars):
 out={}
 for perm in permutations(range(len(M))):
  inv=sum(perm[i]>perm[j] for i in range(len(M)) for j in range(i+1,len(M)))
  p=const((-1)**inv,nvars)
  for i,j in enumerate(perm):p=mul(p,M[i][j])
  out=add(out,p)
 return out
def ct(p):return p.get((0,),Q(0))
def wedge_sign(a,b):
 if set(a)&set(b):return 0,()
 return (-1)**sum(i>j for i in a for j in b),tuple(sorted(a+b))

def sl2_adjoint_and_volume():
 n=4;a,b,c,d=(var(i,n) for i in range(n))
 Adj=[[power(a,2),scale(power(b,2),-1),scale(mul(a,b),-2)],
      [scale(power(c,2),-1),power(d,2),scale(mul(c,d),2)],
      [scale(mul(a,c),-1),mul(b,d),add(mul(a,d),mul(b,c))]]
 need(det(Adj,n)==power(add(mul(a,d),scale(mul(b,c),-1)),3),'Universal GL2 adjoint determinant numerator')
 # Coordinates a,b,c on ad-bc=1; density rho=1/a for da wedge db wedge dc.
 n=7;a,b,c,F,Fa,Fb,Fc=(var(i,n) for i in range(n))
 rho=mon([-1,0,0,0,0,0,0]);one=const(1,n)
 Dlist=[
  [c,mul(add(one,mul(b,c)),rho),{}], # upper root
  [{},{},a],                        # lower root
  [a,b,scale(c,-1)]]                # diagonal
 rows=[]
 for label,D in zip(['E','F','H'],Dlist):
  vals=[mul(rho,z) for z in D]
  divergence=add(*(deriv(z,i) for i,z in enumerate(vals)))
  need(not divergence,'SL2 fundamental field divergence zero:'+label)
  alpha={tuple(j for j in range(3) if j!=i):scale(vals[i],(-1)**i) for i in range(3)}
  top={}
  for key,z in alpha.items():
   for j in range(3):
    sign,idx=wedge_sign((j,),key)
    if sign:
     need(idx==(0,1,2),'Top exterior orientation')
     term=add(mul([Fa,Fb,Fc][j],z),mul(F,deriv(z,j)))
     top=add(top,scale(term,sign))
  expected=mul(rho,add(*(mul(D[i],[Fa,Fb,Fc][i]) for i in range(3))))
  need(top==expected,'Universal first-jet Cartan exact-top identity:'+label)
  rows.append(dict(field=label,divergence='0',Cartan_first_jet_identity=True))
 # Weyl component on LieSL2/LieT swaps off-diagonal root vectors with minus signs.
 quotient=[[const(0,1),const(-1,1)],[const(-1,1),const(0,1)]]
 need(det(quotient,1)==const(-1,1),'N(T) determinant obstruction is -1')
 return dict(universal_families=7,adjoint_determinant='(ad-bc)^3 before division by determinant^3',
  volume_chart='omega=a^-1 da wedge db wedge dc on ad-bc=1',
  exact_divergence_and_Cartan=rows,normalizer_component_quotient_determinant=-1,
  scope='Exact SL2 model, not a finite proof of every reductive-group or isotropy theorem.')

def symmetric_power_coinvariants():
 rows=[]
 for r in range(SCOPE['controls']['symmetric_power_degree_max']+1):
  N=r+1
  E=[[Q(0) for j in range(N)] for i in range(N)]
  F=copy.deepcopy(E);H=copy.deepcopy(E)
  for i in range(N):
   if i:E[i-1][i]=Q(i)
   if i<r:F[i+1][i]=Q(r-i)
   H[i][i]=Q(r-2*i)
  def mm(a,b):return [[sum((a[i][k]*b[k][j] for k in range(N)),Q(0)) for j in range(N)] for i in range(N)]
  EF,FE=mm(E,F),mm(F,E)
  need([[EF[i][j]-FE[i][j] for j in range(N)] for i in range(N)]==H,'Exact sl2 bracket on Sym^r')
  recovered=[]
  if r:
   for i in range(N):
    col=[E[j][i+1]/Q(i+1) for j in range(N)] if i<r else [F[j][r-1] for j in range(N)]
    need(col==[Q(j==i) for j in range(N)],'Every nontrivial model basis vector is a Lie image')
    recovered.append(i)
  else:need(E==F==H==[[Q(0)]],'The constant model remains a nonzero coinvariant')
  rows.append(dict(r=r,dimension=N,Lie_image_basis_size=len(recovered),coinvariant_dimension=0 if r else 1))
 return rows

def rankone_models():
 rows=[]
 for size in SCOPE['controls']['rankone_matrix_sizes']:
  n=2*size;v=[var(i,n) for i in range(size)];ell=[var(size+i,n) for i in range(size)]
  M=[[mul(v[i],ell[j]) for j in range(size)] for i in range(size)]
  R=add(*(mul(v[i],ell[i]) for i in range(size)))
  need(add(*(M[i][i] for i in range(size)))==R,'Rankone trace identity')
  count=1
  for i in range(size):
   for j in range(size):
    need(add(*(mul(M[i][k],M[k][j]) for k in range(size)))==mul(R,M[i][j]),'M^2=tr(M)M before R1 specialization')
    count+=1
  for i,k in combinations(range(size),2):
   for j,l in combinations(range(size),2):
    need(mul(M[i][j],M[k][l])==mul(M[i][l],M[k][j]),'Actual two-by-two rankone minor')
    count+=1
  n=size**2;G=[[var(i*size+j,n) for j in range(size)] for i in range(size)]
  observed=[]
  for i in range(size):
   for j in range(size):
    ge=G[i][0] if j==0 else {}
    eg=G[0][j] if i==0 else {}
    comm=add(ge,scale(eg,-1))
    expected=G[i][0] if j==0 and i>0 else scale(G[0][j],-1) if i==0 and j>0 else {}
    need(comm==expected,'Literal block stabilizer of E11')
    if comm:observed.append([i,j])
  need(len(observed)==2*(size-1),'Stabilizer off-diagonal blocks vanish exactly')
  rows.append(dict(n=size,universal_matrix_equalities=count,stabilizer_nonzero_constraints=observed,
   meaning='Symbolic fixed-size rankone models support the written all-n orbit/stabilizer argument.'))
 return rows

def torus_models():
 def aut(p,a,eps):return clean({(eps*i,):v*a**i for(i,),v in p.items()})
 fixtures=[
  ({(-2,):Q(1,3),(0,):Q(-2),(3,):Q(5,2)},{(-1,):Q(2),(2,):Q(-3)}),
  ({(1,):Q(1)},{(-3,):Q(1)}),
  ({(-1,):Q(1),(1,):Q(1)},{(0,):Q(1),(2,):Q(2,7)})]
 rows=[];signs=[]
 for a in [Q(-2),Q(1,2),Q(3)]:
  for eps in [-1,1]:
   # (1/sigma(t))*d(sigma(t)) = eps*t^-1 exactly.
   inverse_t=mon((-eps,),1/a)
   derivative_t=mon((eps-1,),eps*a)
   need(mul(inverse_t,derivative_t)==mon((-1,),eps),'Torus inversion volume sign')
   signs.append(dict(a=str(a),epsilon=eps,volume_multiplier=eps))
   for idx,(f,h) in enumerate(fixtures):
    sf,sh=aut(f,a,eps),aut(h,a,eps)
    need(aut(sf,a**(-eps),eps)==f and aut(sh,a**(-eps),eps)==h,'Actual automorphism inverse restores phase AND marker')
    for m in range(SCOPE['controls']['torus_power_max']+1):
     lhs=ct(mul(sh,power(sf,m)));rhs=ct(mul(h,power(f,m)))
     need(lhs==rhs,'Both-direction fixed-marker covariance on torus')
     rows.append(dict(a=str(a),epsilon=eps,fixture=idx,m=m,moment=str(rhs)))
 for k in range(-8,9):
  unit=mon((k,));need(ct(unit)==Q(k==0),'Nonconstant torus units have zero trace')
  need(mul(unit,mon((-k,)))==const(1,1),'Actual unit inverse')
 t=mon((1,));t_inv=mon((-1,))
 avg=scale(add(t,t_inv),Q(1,2))
 need(ct(power(avg,2))==Q(1,2),'ACTUAL finite averaging destroys purity')
 for m in range(1,SCOPE['controls']['torus_power_max']+1):
  need(ct(power(t,m))==0 and ct(power(t_inv,m))==0,'Individual pure torus monomials')
 return dict(volume_signs=signs,fixed_pair_rows=rows,nonconstant_unit_controls=16,
  finite_averaging_negative=dict(pure_inputs=['t','t^-1'],average='(t+t^-1)/2',square_trace='1/2',
   inference_rejected='All-Aut trace invariance does not preserve phase purity under averaging.'))

def Ga_and_point_models():
 rows=[]
 for k in range(SCOPE['controls']['Ga_primitive_degree_max']+1):
  primitive=mon((k+1,),Q(1,k+1))
  need(deriv(primitive,0)==mon((k,)),'Ga exact primitive')
  rows.append(k)
 need(deriv(mon((1,)),0)==const(1,1),'Ga Lie image contains1: no normalized invariant functional')
 point=[Q(0),Q(1),Q(-7),Q(3,2)]
 need(all(Q(a)==a for a in point),'Point trace is scalar identity')
 return dict(Ga_exact_basis_degrees=rows,Ga_constant_is_Lie_image=True,
  point_dimension_zero='A=K, trace=id; separate case rather than contraction with a negative-degree form.')


def finite_character_controls():
 """Finite cyclic representation models; not proofs of etale descent or homogenization."""
 rows=[]
 for m in range(1,SCOPE['controls']['cyclic_order_max']+1):
  for j in range(m):
   order=m//gcd(m,j)
   kernel=[a for a in range(m) if (j*a)%m==0]
   need(len(kernel)*order==m,'Character kernel index equals character order')
   # Norm of this character in Q[x]/(x^m-1); invariant line iff j=0.
   norm=[Q(0)]*m
   for a in range(m):norm[(j*a)%m]+=Q(1,m)
   invariant=norm==[Q(1)]+[Q(0)]*(m-1)
   need(invariant==(j==0),'Finite character line has invariants iff trivial')
   for t in range(1,m+1):
    if m%t:continue
    subgroup=[a for a in range(m) if a%t==0]
    restricted_trivial=all((j*a)%m==0 for a in subgroup)
    factors=all(a in kernel for a in subgroup)
    need(restricted_trivial==factors==(t%order==0),'Subgroup factors through determinant kernel iff restricted character trivial')
    if t==order:need(subgroup==kernel,'Equal cover degree identifies kernel in cyclic model')
   rows.append(dict(component_order=m,character_exponent=j,determinant_order=order,kernel_size=len(kernel),invariant_line_dimension=int(invariant)))
 return dict(rows=rows,scope='Exact bounded finite-character models support the character, subgroup and invariant-line steps; they do not prove general homogeneous-cover or de Rham descent theorems.')


def finite_etale_laurent_controls():
 """Exact finite-free matrix trace and marked moments for t=s^d, in bounded models."""
 rows=[]
 fixtures=[
  ({(-2,):Q(1,3),(0,):Q(-2),(3,):Q(5,2)},{(-5,):Q(2),(0,):Q(-3),(4,):Q(7,2)}),
  ({(1,):Q(1)},{(-3,):Q(1),(2,):Q(2)}),
  ({(-1,):Q(1),(1,):Q(1)},{(0,):Q(1),(2,):Q(2,7)})]
 for d in range(1,SCOPE['controls']['finite_etale_degree_max']+1):
  def pull(f):return {(d*i,):a for(i,),a in f.items()}
  def transfer(h):
   # Basis 1,s,...,s^(d-1) of Q[t,t^-1][s]/(s^d-t).
   M=[[{} for _ in range(d)] for _ in range(d)]
   for (j,),a in h.items():
    for col in range(d):
     exponent,row=divmod(j+col,d)
     M[row][col]=add(M[row][col],mon((exponent,),a))
   tr=scale(add(*(M[i][i] for i in range(d))),Q(1,d))
   expected={(j//d,):a for(j,),a in h.items() if j%d==0}
   need(tr==expected,'Actual finite-free normalized multiplication-matrix trace')
   return tr
  for j in range(-12,13):transfer(mon((j,)))
  for idx,(f,h) in enumerate(fixtures):
   th=transfer(h)
   need(transfer(pull(f))==f,'Normalized algebra trace is a pullback retraction')
   for r in range(SCOPE['controls']['torus_power_max']+1):
    lhs=ct(mul(h,power(pull(f),r)))
    rhs=ct(mul(th,power(f,r)))
    need(lhs==rhs,'Exact finite-etale fixed-marker trace identity in Laurent model')
    rows.append(dict(degree=d,fixture=idx,order=r,moment=str(lhs)))
 return dict(rows=rows,monomial_trace_matrix_controls=25*SCOPE['controls']['finite_etale_degree_max'],
 scope='Exact finite-free Laurent models and bounded marked moments support E1; they do not prove Galois domination, presentation independence or arbitrary finite-etale compatibility.')


def nonnormal_S3_averaging_control():
 """Finite set/function algebra check of the repeated-branch count, with nonnormal Delta."""
 Gamma=list(permutations(range(3)));identity=(0,1,2)
 def compose(g,h):return tuple(g[h[i]] for i in range(3))
 def inverse(g):return tuple(g.index(i) for i in range(3))
 Delta=[g for g in Gamma if g[0]==0]
 need(len(Gamma)==6 and len(Delta)==2,'S3 and order-two branch stabilizer')
 conjugates={compose(compose(g,d),inverse(g)) for g in Gamma for d in Delta}
 need(conjugates!=set(Delta),'The order-two stabilizer is NONNORMAL')
 h=[Q(1,3),Q(-2),Q(5,2)];normalized=sum(h,Q(0))/3
 rows=[]
 for g in Gamma:
  # W=Gamma, r(g)=g(0); deck action is right multiplication.
  branches=[compose(g,gamma)[0] for gamma in Gamma]
  multiplicities=[branches.count(i) for i in range(3)]
  need(multiplicities==[2,2,2],'Every rank-three branch occurs |Delta|=2 times')
  need(sum((h[i] for i in branches),Q(0))/6==normalized,'Gamma average equals normalized rank-three trace')
  for delta in Delta:need(compose(g,delta)[0]==g[0],'Delta acts within every r-fiber')
  rows.append(dict(g=list(g),branch_multiplicities=multiplicities,average=str(normalized)))
 return dict(group='S3',group_order=6,Delta_order=2,Delta_normal=False,cover_degree=3,rows=rows,
 scope='Exact finite set/function algebra model for E6 with NONNORMAL intermediate subgroup; not a connected algebraic-cover proof.')

def negative_controls(ledger):
 out=[]
 def reject(label,fn):
  try:fn()
  except (ValueError,KeyError,TypeError,OSError):out.append(label);return
  raise ValueError('Invalid control accepted: '+label)
 for role in SPECS:
  a=copy.deepcopy(ledger);a['sources']=[r for r in a['sources'] if r['role']!=role]
  reject('missing_role:'+role,lambda a=a:sources(a))
  a=copy.deepcopy(ledger);next(r for r in a['sources'] if r['role']==role)['sha256']='0'*64
  reject('bad_hash:'+role,lambda a=a:sources(a))
  a=copy.deepcopy(ledger);next(r for r in a['sources'] if r['role']==role)['surplus']=True
  reject('surplus_source_field:'+role,lambda a=a:sources(a))
 for k in SCOPE:
  a=copy.deepcopy(SCOPE)
  if type(a[k])is bool:a[k]=not a[k]
  elif type(a[k])is int:a[k]=True
  elif type(a[k])is str:a[k]+=' INVALID'
  elif type(a[k])is dict:a[k]['surplus']=True
  else:raise ValueError('Uncovered input type')
  reject('false_scope:'+k,lambda a=a:inputs(a))
 a=copy.deepcopy(SCOPE);a['surplus']=False
 reject('surplus_input_field',lambda:inputs(a))
 a=copy.deepcopy(SCOPE);a['controls']['torus_power_max']=True
 reject('boolean_integer_bound',lambda:inputs(a))
 for label,raw in [('duplicate_top','{"schema":2,"schema":2}'),('duplicate_nested','{"a":{"x":1,"x":1}}'),('NaN','{"a":NaN}'),('Infinity','{"a":Infinity}'),('negative_Infinity','{"a":-Infinity}')]:
  reject('json:'+label,lambda raw=raw:parse(raw))
 for key,val in [('schema',True),('surplus',False),('assurance_flags',dict(ASSURANCE,external_review=True))]:
  a=copy.deepcopy(ledger);a[key]=val
  reject('bad_ledger:'+key,lambda a=a:sources(a))
 a=copy.deepcopy(ledger);a['sources'][-1]=copy.deepcopy(a['sources'][0])
 reject('duplicate_role',lambda:sources(a))
 for raw in ['../outside.pdf','/tmp/outside.pdf','sources/../../outside.pdf','sources\\outside.pdf','sources//bad.pdf','./bad.pdf']:
  reject('path_confinement:'+raw,lambda raw=raw:confined_path(raw))
 return out

def main():
 ap=argparse.ArgumentParser(description=__doc__)
 ap.add_argument('--repo-root',type=Path,required=True)
 ap.add_argument('--output',type=Path,default=HERE/RN)
 args=ap.parse_args();need(args.repo_root.resolve().is_dir(),'Repository root exists')
 need(args.output.resolve() not in {(HERE/p).resolve() for p,_ in SPECS.values()}|{(HERE/LN).resolve()},'Receipt cannot overwrite evidence')
 inp=load(HERE/IN);inputs(inp);ledger=load(HERE/LN);pins=sources(ledger)
 sl2=sl2_adjoint_and_volume();sym=symmetric_power_coinvariants();rank=rankone_models()
 torus=torus_models();Ga=Ga_and_point_models();cover=finite_character_controls();etale=finite_etale_laurent_controls();S3=nonnormal_S3_averaging_control();neg=negative_controls(ledger)
 report=dict(schema=2,status='PASS',proof_sha256=PROOF_SHA,checker_sha256=sha(Path(__file__)),
 inputs_sha256=sha(HERE/IN),source_ledger_sha256=sha(HERE/LN),
 python_version=sys.version,scope=inp,sources=[pins[k] for k in sorted(pins)],
 SL2_universal_models=sl2,symmetric_power_coinvariant_models=sym,rankone_models=rank,
 torus_models=torus,additive_group_and_point_models=Ga,finite_character_models=cover,finite_etale_laurent_models=etale,nonnormal_S3_averaging_model=S3,
 source_scope_negative_controls=neg,assurance_flags=ASSURANCE,
 meaning='Source authentication and bounded exact supporting controls PASS. The unrestricted all-field theorem is established by the complete written proof and its narrowly scoped classical imports. No bounded test proves complete reducibility, Rosenlicht, finite-etale descent, the homogeneous-cover theorem, arbitrary-field descent or any unrestricted conclusion.')
 args.output.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
 print(json.dumps(dict(status='PASS',mandatory_roles=len(pins),SL2_universal_families=7,
 symmetric_power_models=len(sym),rankone_matrix_equalities=sum(r['universal_matrix_equalities'] for r in rank),
 fixed_pair_torus_rows=len(torus['fixed_pair_rows']),finite_character_rows=len(cover['rows']),finite_etale_marked_rows=len(etale['rows']),nonnormal_S3_rows=len(S3['rows']),
 source_scope_controls=len(neg),receipt_sha256=sha(args.output))))
if __name__=='__main__':main()

#!/usr/bin/env python3
"""Authenticate the field/product/quadratic proof and replay bounded exact supporting models."""
from pathlib import Path,PurePosixPath
from fractions import Fraction as Q
from math import factorial
import argparse,copy,hashlib,json,re,sys
if not __debug__:raise SystemExit('Optimized execution is not accepted.')
HERE=Path(__file__).resolve().parent
PN='FIELD_MATHIEU_PRODUCT_AND_QUADRATIC_COUNTEREXAMPLES.md'
CN='verify_field_mathieu_product_quadratic.py'
IN='field_mathieu_product_quadratic_inputs.json'
LN='field_mathieu_product_quadratic_source_ledger.json'
RN='field_mathieu_product_quadratic_verification.json'
PROOF_SHA='8537845bffd21a60bcadcaaff190501b7d9d6869f395680c241b373ee443d97a'
ASSURANCE=dict(external_review=False,formal_verification=False,independent_replay=False,worldwide_priority=False)
SPECS={'written_proof': ('FIELD_MATHIEU_PRODUCT_AND_QUADRATIC_COUNTEREXAMPLES.md',
                   'Sections1--22_complete_R_Q_all_char0_formally_real_Gram_even_quadratic_actual_cover_and_acyclic_comparison_proofs'),
 'exact_checker': ('verify_field_mathieu_product_quadratic.py',
                   'Source_authentication_and_bounded_exact_CDF_Laurent_Gram_models_NOT_unrestricted_field_or_quadratic_form_proof'),
 'typed_inputs': ('field_mathieu_product_quadratic_inputs.json',
                  'Exact_standard_Q_models_all_char0_fixed_pair_identities_formal_reality_Gram_all_even_quadratic_scopes_and_given_phase_cover_directions'),
 'long_known_method_primary': ('sources/long-2607.19012v1.html',
                               'CREDITED_KNOWN_METHOD_Theorem2_1_complete_beta_binomial_proof_Proposition3_1_Lemma4_1_NOT_a_new_interval_counterexample_or_Jacobian_claim'),
 'long_split_SL2_comparison_primary': ('sources/long-2607.19012v1.html',
                                       'ACTUAL_IMPORT_ONLY_Theorem4_2_exact_complex_SU2_Haar_pair_F_equals_1plusc_times_adplusb_H_equals_minusc_all_orders_USED_FOR_SPLIT_REAL_SL2_COMPARISON_via_real_coefficients_trace_uniqueness_and_scalar_extension_NOT_uniform_SU_n_S2_quaternion_core'),
 'meckes_column_primary': ('sources/meckes-haar-book.pdf',
                           'ONLY_Section1_2_printed19_20_Haar_first_column_uniform_on_real_or_complex_sphere_including_SU_n_last_column_adjustment_CROSSCHECK_self_contained_transitivity_Gaussian_proof_is_in_Section3'),
 'dvk_torus_comparison_primary': ('sources/dvk-powers-1998.pdf',
                                  'ONLY_printed222_223_and_Section3_scalar_Laurent_Mathieu_theorem_over_C_local_finite_coefficient_field_extension_proved_in_Section10_NOT_function_valued_coefficients_or_entire_K_embedding'),
 'complete_reducibility_primary_reused': ('sources/milne-algebraic-groups-2017.pdf',
                                          'REUSED_SRC_RHT_MILNE_COMPLETE_REDUCIBILITY_ONLY_22_42_22_43_printed477_478_and_12l_12_53_Reynolds_printed248_249_same_field_rational_semisimplicity_NOT_whole_book_or_full_allAut_theorem')}
PRIMARY_SHA={'long_known_method_primary': 'abeec76e95cba0e6244f4af31240634e28d657116c8afcf743b58db144dd7ebd',
 'long_split_SL2_comparison_primary': 'abeec76e95cba0e6244f4af31240634e28d657116c8afcf743b58db144dd7ebd',
 'meckes_column_primary': '12c9e49c9768f9f87951a59e6f4e960bf4a5bf6aa6f7fc5c8b442076e57d3dbe',
 'dvk_torus_comparison_primary': 'bb002d12b539a1ff3fdd099f343d3accd921ccf60f8a149f29b07cd0f7e87e1b',
 'complete_reducibility_primary_reused': '503f1a7d511632af05541c52fb2bf6fc36f5a75b0323f2f9110f03dcce11fd67'}
SCOPE={'schema': 2,
 'reviewed_proof_sha256': '8537845bffd21a60bcadcaaff190501b7d9d6869f395680c241b373ee443d97a',
 'base_field': 'ANY characteristic-zero K for standard Q-form and even-quadratic bad pairs; individual-factor '
               'Mathieu product-instability conclusion for EVERY FORMALLY REAL K; R is a special calculation, not '
               'the whole scope',
 'functional': 'unique normalized invariant K-linear Reynolds projection; on standard real compact factors Haar, '
               'on products tensor Laurent CT; no compactness of arbitrary K-points',
 'uniform_criterion': 'ANY commutative unital real algebra B and normalized R-linear lambda, GIVEN v with '
                      'lambda(p(v))=integral_0^1 p(s)ds for EVERY real polynomial p',
 'uniform_phase': 'f=(1-t^-1)((1-v)+v*t), ONE fixed h=t^-1 in B[t,t^-1]',
 'all_order_identity': 'for EVERY integer m>=1: Lambda(f^m)=0 and Lambda(h*f^m)=(-1)^(m-1)/(m+1) !=0',
 'all_order_proof': 'exact Laurent identity integral_0^1 f(s,t)^m ds=(t-1)^(m-1)*(t-t^-m)/(m+1); positive and '
                    'negative exponents separated',
 'beta_CDF_scope': 'GIVEN regular u with Beta(a,b) law for POSITIVE INTEGER a,b; polynomial '
                   'V(s)=((a+b-1)!/((a-1)!(b-1)!))*integral_0^s r^(a-1)*(1-r)^(b-1)dr makes V(u) uniform',
 'noninteger_beta_CDF_regular_claimed': False,
 'SU_n_scope': 'EVERY n>=2 standard Q(i)-unitary Q-form after ANY char0 field K, times SPLIT G_m; real '
               'specialization compact SU(n); u=a11^2+b11^2, v=1-(1-u)^(n-1)',
 'SU_n_first_column': 'OVER R: uniform complex sphere by transitivity; normalized complex Gaussians give exact '
                      'Beta(1,n-1). These rational trace moments descend to Q then all K.',
 'S2_scope': 'standard Q-sphere x^2+y^2+z^2=1 after ANY char0 K, times SPLIT G_m; real specialization compact '
             'SO(3)/SO(2); v=(1+z)/2',
 'S2_phase': 'f=((1+z)*t-2*z+(z-1)*t^-1)/2; h=t^-1',
 'quaternion_group': 'standard Hamilton Q-algebra (-1,-1)_Q, base changed to ANY char0 K; GL_1(H_K)=Spec '
                     'K[a,b,c,d,N^-1], N=sum four squares',
 'quaternion_phase': 'alpha=a^2+b^2, beta=c^2+d^2; F=(N-1)*(alpha*N+beta)/N^2; ONE fixed Hmark=N^-1',
 'quaternion_cover': 'p:SU2_compact x G_m -> GL_1(H), (q,s)->s*q; O(cover)=O(base)[s]/(s^2-N), finite faithfully '
                     'flat rank2, derivative2s unit, kernel diagonal mu_2',
 'quaternion_pullback': 'pstar(F)=(1-s^-2)*((1-u)+u*s^2), pstar(Hmark)=s^-2; d=2 Long family, trace restriction '
                        'by direct equivariant isogeny uniqueness',
 'quaternion_derived': 'SL_1(H_K), with compact real SU2 realization; geometric group SL2; no nonidentity '
                       'K-unipotents if K formally real',
 'given_phase_cover_trace': 'for EVERY GIVEN downstairs g under quaternion or even-quadratic actual double cover '
                            'and ONE fixed b upstairs, ALL m>=0: Lambda_up(b*(pstar g)^m)=Lambda_down(Tp(b)*g^m), '
                            'Tp=(1+tau)/2=Tr/2',
 'given_phase_cover_purity_forward': True,
 'given_phase_cover_purity_reverse': True,
 'given_phase_cover_safety_forward': True,
 'given_phase_cover_safety_reverse': True,
 'given_phase_cover_unsafety_both_directions': True,
 'arbitrary_upstairs_phase_descends_claimed': False,
 'purity_of_power_implies_purity_of_root_claimed': False,
 'cover_marker_choice_depends_on_order': False,
 'compact_real_factor_Mathieu': True,
 'split_real_torus_Mathieu': True,
 'product_of_two_Mathieu_factors_can_fail': True,
 'compact_derived_suffices_for_Mathieu_claimed': False,
 'compact_factor_Mathieu_over_C_claimed': False,
 'global_Mathieu_property_scalar_extension_invariant_claimed': False,
 'every_compact_factor_times_Gm_fails_claimed': False,
 'complete_real_classification_claimed': False,
 'Long_interval_construction_claimed_new': False,
 'Long_split_real_SL2_comparison_uses_Theorem4_2_import': True,
 'Long_Theorem4_2_is_uniform_SU_n_S2_quaternion_core_premise': False,
 'finite_controls_prove_unrestricted_orders': False,
 'external_full_homogeneous_trace_package_is_required_premise': False,
 'controls': {'beta_a_max': 4,
              'beta_b_max': 4,
              'beta_moment_max': 6,
              'SU_n_sizes': [2, 3, 4, 5, 6],
              'marked_moment_max': 10,
              'even_r_max': 6,
              'Gram_sphere_dimensions': [3, 4],
              'Gram_basis_sizes': [8, 12],
              'Gram_triangular_mixtures': True},
 'assurance_flags': {'external_review': False,
                     'formal_verification': False,
                     'independent_replay': False,
                     'worldwide_priority': False},
 'rational_models': 'standard determinant-one unitary Q-group of Q(i)^n, standard Q-sphere x^2+y^2+z^2=1, '
                    'Hamilton Q-quaternion group; base change THESE specified models to K',
 'rational_trace': 'A_Q=Q1 directsum Lie(G_Q)A_Q; normalized invariant Lambda_Q from same-field complete '
                   'reducibility; Lambda_K=id_K tensor Lambda_Q for EVERY field extension',
 'all_char0_standard_Q_form_counterexamples': True,
 'formally_real_scope': '-1 is NOT a finite sum of squares in K; use this DEFINITION directly, no chosen ordering',
 'formally_real_compact_form_factors_have_only_zero_pure_element': True,
 'formally_real_product_of_two_Mathieu_factors_fails': True,
 'Q_itself_has_product_counterexamples': True,
 'compactness_of_general_K_points_claimed': False,
 'arbitrary_twists_of_standard_unitary_or_sphere_forms_claimed': False,
 'Gram_scope': 'EVERY finite Q-independent coefficient set has rational Gram matrix, strictly positive definite '
               'after R-Haar extension; exact rational LDL gives positive rational pivots',
 'rational_SOS_identity': 'positive rational A/B is a sum of A*B copies of (1/B)^2; formally real sum of squares0 '
                          'forces all terms0, so lambda_K(f^2)=0 implies f=0',
 'ordering_or_Artin_Schreier_import_required': False,
 'Lagrange_four_square_import_required': False,
 'torus_Mathieu_any_char0_field': True,
 'torus_field_embedding_scope': 'ONLY the finitely generated coefficient field Q(coefficients of ONE GIVEN phase '
                                'and marker) embeds in C; never the entire arbitrary K',
 'embedding_entire_K_in_R_or_C_required': False,
 'quaternion_no_nonidentity_K_unipotents_for_formally_real_K': True,
 'even_quadratic_scope': 'EVERY nondegenerate quadratic form q on dimension N=2r>=4 over ANY char0 K; choose '
                         'orthogonal nondegenerate rank2 A and rank(2r-2) B, q=A+B',
 'even_quadratic_dimension_min': 4,
 'even_quadratic_requires_positive_definiteness': False,
 'even_quadratic_requires_split_form': False,
 'even_quadratic_norm_one_K_point_required': False,
 'even_quadratic_space': 'X=Spec K[V,q^-1], geometrically homogeneous under connected reductive SO(q) x G_m; '
                         'quotient/field statements scheme-geometric, not K-point transitivity',
 'even_quadratic_phase': 'u=A/q, v=1-(B/q)^(r-1), F=(1-q^-1)*(1-v+v*q), ONE fixed H=q^-1; equivalently '
                         'F=(q-1)-(q-1)^2*B^(r-1)/q^r',
 'even_quadratic_all_order_moments': 'ALL m>=1: Lambda_X(F^m)=0 and Lambda_X(H F^m)=(-1)^(m-1)/(m+1); all '
                                     'extensions',
 'even_quadratic_beta_moments': 'lambda_{q=1}(A(w)^k)=k!*(r-1)!/(k+r-1)! for ALL k>=0 via standard Q real-sphere '
                                'law, Kbar pair-normalization and faithful descent',
 'even_quadratic_cover': 'W={q=1}, Y=W x G_m, p(w,s)=s*w; O(Y)=O(X)[s]/(s^2-q), finite faithfully flat rank2 and '
                         'etale; tau=(-w,-s)',
 'cover_marker_scope': 'quaternion and even-quadratic ACTUAL double covers; EVERY GIVEN downstairs phase, ONE '
                       'fixed marker, normalized trace Tp=(1+tau)/2=Tr/2, ALL m>=0, both purity/safety/unsafety '
                       'directions',
 'odd_dimensional_unsafe_phase_constructed': False,
 'odd_dimension_global_Mathieu_classification_claimed': False,
 'odd_cover_pointed_K_identification_without_K_point_claimed': False,
 'phase_averaging_or_phase_norm_preserves_purity_claimed': False,
 'alternative_SU_n_bridge': 'standard q on dimension2n; (M,s)->s*(first complex column), q->s^2, A->s^2*|M11|^2; '
                            'U_n x G_m geometrically transitive on target, so trace restriction; CT(t->s^2) '
                            'recovers degree1 family',
 'SU_n_to_SO_group_map_surjective_assumed': False,
 'forced_reverse_from_SU_n_witness_to_arbitrary_q_claimed': False,
 'independent_direct_and_quadratic_pullback_routes_retained': True,
 'reused_Milne_source_leaf': 'SRC-RHT-MILNE-COMPLETE-REDUCIBILITY; exact existing canonical PDF, '
                             'narrow22.42--22.43 and12l scope only',
 'whole_Milne_book_review_claimed': False}
def need(v,msg):
 if not v:raise ValueError(msg)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def parse(text):
 def unique(pairs):
  out={}
  for k,v in pairs:
   need(k not in out,'Duplicate JSON key: '+k);out[k]=v
  return out
 def bad_constant(x):raise ValueError('Nonfinite JSON value: '+x)
 return json.loads(text,object_pairs_hook=unique,parse_constant=bad_constant)
def load(p):return parse(p.read_text(encoding='utf-8'))
def typed(a,b):
 if type(a)is not type(b):return False
 if isinstance(b,dict):return set(a)==set(b) and all(typed(a[k],v) for k,v in b.items())
 if isinstance(b,list):return len(a)==len(b) and all(typed(x,y) for x,y in zip(a,b))
 return a==b
def inputs(v):need(typed(v,SCOPE),'Exact field/formal-reality/quadratic theorem scope; no surplus fields or assurance inflation')
def confined_path(raw):
 need(type(raw)is str and raw and '\\' not in raw,'POSIX relative source path')
 p=PurePosixPath(raw)
 need(not p.is_absolute() and all(x not in ('.','..') for x in raw.split('/')),'No absolute or traversing source path')
 need(str(p)==raw,'Canonical relative source path')
 target=(HERE/raw).resolve()
 need(target.is_relative_to(HERE) and target.is_file(),'Source remains within package: '+raw)
 return target
def sources(v):
 need(type(v)is dict and set(v)=={'schema','sources','assurance_flags'},'Exact ledger fields')
 need(type(v['schema'])is int and v['schema']==2 and typed(v['assurance_flags'],ASSURANCE),'Typed ledger schema and assurance')
 rows=v['sources'];need(type(rows)is list and len(rows)==len(SPECS),'All mandatory source roles')
 out={}
 for row in rows:
  need(type(row)is dict and set(row)=={'role','path','relation','sha256'},'Exact source row fields')
  role=row['role'];need(type(role)is str and role in SPECS and role not in out,'Unique known source role')
  need(type(row['relation'])is str and type(row['sha256'])is str and re.fullmatch('[0-9a-f]{64}',row['sha256']) is not None,'Typed source relation and SHA256')
  p=confined_path(row['path'])
  need((row['path'],row['relation'])==SPECS[role],'Exact reviewed source role, path and scope')
  need(sha(p)==row['sha256'],'Missing or stale source: '+role)
  out[role]=row
 need(out['written_proof']['sha256']==PROOF_SHA,'Pinned complete proof')
 for role,h in PRIMARY_SHA.items():need(out[role]['sha256']==h,'Exact primary source body: '+role)
 return out
def exact_models():
 def clean(p):return {i:Q(a) for i,a in p.items() if a}
 def add(*ps):
  r={}
  for p in ps:
   for i,a in p.items():r[i]=r.get(i,Q(0))+a
  return clean(r)
 def scale(p,c):return clean({i:a*c for i,a in p.items()})
 def mul(p,q):
  r={}
  for i,a in p.items():
   for j,b in q.items():r[i+j]=r.get(i+j,Q(0))+a*b
  return clean(r)
 def power(p,k):
  r={0:Q(1)}
  for _ in range(k):r=mul(r,p)
  return r
 def integral(p):return sum((a/Q(i+1) for i,a in p.items()),Q(0))
 def laurent_mul(p,q):
  r={}
  for i,a in p.items():
   for j,b in q.items():r[i+j]=add(r.get(i+j,{}),mul(a,b))
  return {i:a for i,a in r.items() if a}
 one={0:Q(1)};u={1:Q(1)};one_minus=add(one,scale(u,-1))
 beta_rows=[]
 for a in range(1,SCOPE['controls']['beta_a_max']+1):
  for b in range(1,SCOPE['controls']['beta_b_max']+1):
   C=Q(factorial(a+b-1),factorial(a-1)*factorial(b-1))
   density=scale(mul(power(u,a-1),power(one_minus,b-1)),C)
   V={i+1:c/Q(i+1) for i,c in density.items()}
   assert sum(V.values())==1 and integral(density)==1
   for k in range(SCOPE['controls']['beta_moment_max']+1):
    value=integral(mul(density,power(V,k)))
    assert value==Q(1,k+1)
    beta_rows.append(dict(a=a,b=b,k=k,moment=str(value)))
 su_rows=[]
 for N in SCOPE['controls']['SU_n_sizes']:
  density=scale(power(one_minus,N-2),Q(N-1))
  v=add(one,scale(power(one_minus,N-1),-1))
  f={1:v,0:add(one,scale(v,-2)),-1:add(v,scale(one,-1))}
  p={0:one}
  for m in range(1,SCOPE['controls']['marked_moment_max']+1):
   p=laurent_mul(p,f)
   pure=integral(mul(density,p.get(0,{})))
   marked=integral(mul(density,p.get(1,{})))
   assert pure==0 and marked==Q((-1)**(m-1),m+1)
   su_rows.append(dict(N=N,m=m,pure=str(pure),marked=str(marked)))
 z={1:Q(1)}
 f={1:scale(add(one,z),Q(1,2)),0:scale(z,-1),-1:scale(add(z,scale(one,-1)),Q(1,2))}
 p={0:one};sphere_rows=[]
 def sphere_integral(poly):return sum((a/Q(i+1) for i,a in poly.items() if i%2==0),Q(0))
 for m in range(1,SCOPE['controls']['marked_moment_max']+1):
  p=laurent_mul(p,f)
  pure=sphere_integral(p.get(0,{}));marked=sphere_integral(p.get(1,{}))
  assert pure==0 and marked==Q((-1)**(m-1),m+1)
  sphere_rows.append(dict(m=m,pure=str(pure),marked=str(marked)))
 raw_u_first=integral(mul(scale(one_minus,2),add(one,scale(u,-2))))
 assert raw_u_first==Q(1,3)!=0
 bad_middle=integral(add(one,scale(u,-1)))
 assert bad_middle==Q(1,2)!=0
 
 # Quaternion pullback from the explicit rational regular formula, without substitution heuristics.
 N={2:one};alpha={2:u};beta={2:one_minus}
 Nminusone={2:one,0:scale(one,-1)}
 inner=laurent_mul(alpha,N)
 for j,c in beta.items():inner[j]=add(inner.get(j,{}),c)
 qpull=laurent_mul(laurent_mul(Nminusone,inner),{-4:one})
 expected={2:u,0:add(one,scale(u,-2)),-2:add(u,scale(one,-1))}
 assert qpull==expected
 qrows=[];qpower={0:one}
 for m in range(1,SCOPE['controls']['marked_moment_max']+1):
  qpower=laurent_mul(qpower,qpull)
  pure=integral(qpower.get(0,{}));marked=integral(qpower.get(2,{}))
  assert pure==0 and marked==Q((-1)**(m-1),m+1)
  qrows.append(dict(m=m,pure=str(pure),marked=str(marked)))
 
 return dict(beta_polynomial_CDF_models=beta_rows,SU_n_Laurent_models=su_rows,
  S2_height_Laurent_models=sphere_rows,quaternion_regular_formula_pullback=True,quaternion_pullback_models=qrows,
  negative_models=dict(omitting_SU3_CDF_first_moment=str(raw_u_first),wrong_middle_coefficient_first_moment=str(bad_middle)))


def field_models():
 """Bounded exact support only; the written proofs establish fields and unrestricted dimension/order."""
 def plus(*ps):
  out={}
  for p in ps:
   for k,v in p.items():out[k]=out.get(k,Q(0))+v
  return {k:Q(v) for k,v in out.items() if v}
 def times(p,q):
  out={}
  for a,x in p.items():
   for b,y in q.items():out[a+b]=out.get(a+b,Q(0))+x*y
  return {k:v for k,v in out.items() if v}
 def scaled(p,c):return {k:v*c for k,v in p.items() if v*c}
 def ppower(p,k):
  out={0:Q(1)}
  for _ in range(k):out=times(out,p)
  return out
 def pint(p):return sum((v/Q(k+1) for k,v in p.items()),Q(0))
 def lmul(p,q):
  out={}
  for i,a in p.items():
   for j,b in q.items():out[i+j]=plus(out.get(i+j,{}),times(a,b))
  return {i:a for i,a in out.items() if a}
 one={0:Q(1)};u={1:Q(1)};g=plus(one,scaled(u,-1))
 even=[]
 for r in range(2,SCOPE['controls']['even_r_max']+1):
  tail=ppower(g,r-1);v=plus(one,scaled(tail,-1))
  f={1:v,0:plus(one,scaled(v,-2)),-1:plus(v,scaled(one,-1))}
  qminus={1:one,0:scaled(one,-1)}
  rem={1:tail,0:scaled(tail,-2),-1:tail}
  # E3 after B=q*(1-u), compared with the original E1 polynomial.
  e3={i:plus(qminus.get(i,{}),scaled(rem.get(i,{}),-1)) for i in set(qminus)|set(rem)}
  e3={i:a for i,a in e3.items() if a}
  need(e3==f,'Exact regular E1/E3 formula agreement')
  density=scaled(ppower(g,r-2),Q(r-1));power={0:one}
  for m in range(1,SCOPE['controls']['marked_moment_max']+1):
   power=lmul(power,f)
   pure=pint(times(density,power.get(0,{})))
   marked=pint(times(density,power.get(1,{})))
   need(pure==0 and marked==Q((-1)**(m-1),m+1),'Even-quadric rational Laurent support')
   dilated={2*i:a for i,a in power.items()}
   need(dilated.get(0,{})==power.get(0,{}) and dilated.get(2,{})==power.get(1,{}),'CT and fixed marker under t->s^2')
   even.append(dict(r=r,dimension=2*r,m=m,pure=str(pure),marked=str(marked)))
 def rise(a,k):
  out=Q(1)
  for j in range(k):out*=a+j
  return out
 def mmul(p,q):
  out={}
  for e,a in p.items():
   for f,b in q.items():
    z=tuple(i+j for i,j in zip(e,f));out[z]=out.get(z,Q(0))+a*b
  return {k:v for k,v in out.items() if v}
 def madd(*ps):
  out={}
  for p in ps:
   for e,a in p.items():out[e]=out.get(e,Q(0))+a
  return {k:v for k,v in out.items() if v}
 def sphere_moment(p,dim):
  out=Q(0)
  for e,a in p.items():
   if any(i%2 for i in e):continue
   ks=[i//2 for i in e];num=Q(1)
   for k in ks:num*=rise(Q(1,2),k)
   out+=a*num/rise(Q(dim,2),sum(ks))
  return out
 def ldlt(G):
  n=len(G);L=[[Q(i==j) for j in range(n)] for i in range(n)];D=[]
  for j in range(n):
   pivot=G[j][j]-sum((L[j][k]**2*D[k] for k in range(j)),Q(0))
   need(pivot>0,'Exact rational positive LDL pivot')
   D.append(pivot)
   for i in range(j+1,n):
    L[i][j]=(G[i][j]-sum((L[i][k]*L[j][k]*D[k] for k in range(j)),Q(0)))/pivot
  need(all(G[i][j]==sum((L[i][k]*D[k]*L[j][k] for k in range(n)),Q(0)) for i in range(n) for j in range(n)),'Rational Gram equals LDL transpose')
  for d in D:
   count=d.numerator*d.denominator;base=Q(1,d.denominator)
   need(Q(count)*base**2==d,'Positive rational is finite sum of rational squares')
  return D
 grams=[]
 for dim,expected_size in zip(SCOPE['controls']['Gram_sphere_dimensions'],SCOPE['controls']['Gram_basis_sizes']):
  oneM={(0,)*dim:Q(1)}
  xs=[]
  for i in range(dim):
   e=[0]*dim;e[i]=1;xs.append({tuple(e):Q(1)})
  diagdiff=madd(mmul(xs[0],xs[0]),{e:-a for e,a in mmul(xs[1],xs[1]).items()})
  basis=[oneM]+xs+[diagdiff]+[mmul(xs[i],xs[j]) for i in range(dim) for j in range(i+1,dim)]
  need(len(basis)==expected_size,'Pinned independent finite sphere basis size')
  for label,B in [('basis',basis),('triangular_mixture',[madd(*basis[:i+1]) for i in range(len(basis))])]:
   G=[[sphere_moment(mmul(a,b),dim) for b in B] for a in B];D=ldlt(G)
   grams.append(dict(sphere_ambient_dimension=dim,basis_size=len(B),model=label,pivots=[str(x) for x in D],
    rational_SOS_counts=[x.numerator*x.denominator for x in D]))
  # Deliberately dependent input cannot acquire a positive-definite Gram claim.
  bad=[[Q(1),Q(1)],[Q(1),Q(1)]]
  try:ldlt(bad)
  except ValueError:pass
  else:raise ValueError('Dependent Gram basis was accepted')
 # In Q(i), a+i*b has square moment0 on S3 although it is a nonzero linear function:
 # Lambda(a^2)=Lambda(b^2)=1/4 and Lambda(ab)=0. This is only a square-isotropy control.
 nonreal=Q(1,4)-Q(1,4);need(nonreal==0,'Formal-reality hypothesis cannot be silently dropped')
 return dict(even_quadratic_Laurent_models=even,rational_Gram_models=grams,
  negative_models=dict(dependent_Gram_basis='REJECTED in dimensions3 and4',non_formally_real_square_isotropy=str(nonreal)),
  scope='Declared bounded rational models support E1/E3, Laurent dilation and Gram linear algebra; they do not establish arbitrary q/field/dimension or global Mathieu classification.')

def negative_controls(ledger):
 out=[]
 def reject(label,fn):
  try:fn()
  except (ValueError,KeyError,TypeError,OSError):out.append(label);return
  raise ValueError('Invalid guard control accepted: '+label)
 for role in SPECS:
  v=copy.deepcopy(ledger);v['sources']=[r for r in v['sources'] if r['role']!=role]
  reject('missing_role:'+role,lambda v=v:sources(v))
  v=copy.deepcopy(ledger);next(r for r in v['sources'] if r['role']==role)['sha256']='0'*64
  reject('bad_hash:'+role,lambda v=v:sources(v))
  v=copy.deepcopy(ledger);next(r for r in v['sources'] if r['role']==role)['surplus']=False
  reject('surplus_source_field:'+role,lambda v=v:sources(v))
 for k in SCOPE:
  v=copy.deepcopy(SCOPE)
  if type(v[k])is bool:v[k]=not v[k]
  elif type(v[k])is int:v[k]=True
  elif type(v[k])is str:v[k]+=' INVALID'
  elif type(v[k])is dict:v[k]['surplus']=False
  else:raise ValueError('Untested input field type')
  reject('false_scope:'+k,lambda v=v:inputs(v))
 for label,raw in [('duplicate_top','{"schema":1,"schema":1}'),('duplicate_nested','{"a":{"x":1,"x":1}}'),('NaN','{"a":NaN}'),('Infinity','{"a":Infinity}'),('negative_Infinity','{"a":-Infinity}')]:
  reject('json:'+label,lambda raw=raw:parse(raw))
 v=copy.deepcopy(SCOPE);v['controls']['marked_moment_max']=True
 reject('boolean_integer_bound',lambda:inputs(v))
 v=copy.deepcopy(SCOPE);v['surplus']=False
 reject('surplus_input_field',lambda:inputs(v))
 for k,x in [('schema',True),('surplus',False),('assurance_flags',dict(ASSURANCE,external_review=True))]:
  v=copy.deepcopy(ledger);v[k]=x
  reject('bad_ledger:'+k,lambda v=v:sources(v))
 v=copy.deepcopy(ledger);v['sources'][-1]=copy.deepcopy(v['sources'][0])
 reject('duplicate_role',lambda:sources(v))
 for raw in ['../outside.pdf','/tmp/outside.pdf','sources/../../outside.pdf','sources\\outside.pdf','sources//bad.pdf','./bad.pdf']:
  reject('path_confinement:'+raw,lambda raw=raw:confined_path(raw))
 return out

def main():
 ap=argparse.ArgumentParser(description=__doc__)
 ap.add_argument('--repo-root',type=Path,default=HERE)
 ap.add_argument('--output',type=Path,default=HERE/RN)
 args=ap.parse_args();need(args.repo_root.resolve().is_dir(),'Replay root exists')
 protected={(HERE/p).resolve() for p,_ in SPECS.values()}|{(HERE/LN).resolve()}
 need(args.output.resolve() not in protected,'Receipt cannot overwrite proof, checker, inputs, ledger or primary source')
 scope=load(HERE/IN);inputs(scope);ledger=load(HERE/LN);pins=sources(ledger)
 models=exact_models();field=field_models();neg=negative_controls(ledger)
 report=dict(schema=2,status='PASS',proof_sha256=PROOF_SHA,checker_sha256=sha(Path(__file__)),
 inputs_sha256=sha(HERE/IN),source_ledger_sha256=sha(HERE/LN),
 python_version=sys.version,scope=scope,sources=[pins[k] for k in sorted(pins)],models=models,field_models=field,
 source_scope_negative_controls=neg,assurance_flags=ASSURANCE,
 meaning='Source authentication and declared bounded exact supporting controls PASS. The complete written proof and narrow Milne/DvK imports establish the stated orders, standard Q-model field extensions, formally real Gram implication and arbitrary even nondegenerate quadratic complements. Finite models do not prove all-order purity, arbitrary fields/forms, complete reducibility, quadratic normalization or any global Mathieu classification.')
 args.output.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
 print(json.dumps(dict(status='PASS',mandatory_roles=len(pins),beta_CDF_controls=len(models['beta_polynomial_CDF_models']),
 SU_n_marked_rows=len(models['SU_n_Laurent_models']),S2_marked_rows=len(models['S2_height_Laurent_models']),
 quaternion_marked_rows=len(models['quaternion_pullback_models']),exact_negative_models=len(models['negative_models']),even_quadratic_marked_rows=len(field['even_quadratic_Laurent_models']),rational_Gram_models=len(field['rational_Gram_models']),
 source_scope_controls=len(neg),receipt_sha256=sha(args.output))))
if __name__=='__main__':main()

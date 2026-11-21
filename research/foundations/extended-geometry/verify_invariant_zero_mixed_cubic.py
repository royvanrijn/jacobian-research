#!/usr/bin/env python3
"""Bounded exact identities, retained certificate replay and narrow source authentication."""
from pathlib import Path
from fractions import Fraction as Q
from math import comb
import argparse,copy,hashlib,json
if not __debug__:raise SystemExit('Optimized execution is not accepted.')
HERE=Path(__file__).resolve().parent
PN='INVARIANT_ZERO_MIXED_CUBIC_LOW_DEGREE_CLASSIFICATION.md'
CN='verify_invariant_zero_mixed_cubic.py';IN='invariant_zero_mixed_cubic_inputs.json'
LN='invariant_zero_mixed_cubic_source_ledger.json';RN='invariant_zero_mixed_cubic_verification.json'
PROOF_SHA='0b65f1da0010eaa646355670099cdd0bbd111ace80abb85f76eeccc1e5547fc9'
PREFIX='research/foundations/extended-geometry/'
ASSURANCE=dict(external_review=False,formal_verification=False,independent_replay=False,worldwide_priority=False)
SPECS={
 'written_proof':('local',PN,'Sections1--10_complete_depression_central_classification_entire_I0_and_negative_center_height3_nonpower'),
 'exact_checker':('local',CN,'bounded_algebra_certificate_and_source_controls_NOT_prefix_inference_of_rational_u'),
 'typed_inputs':('local',IN,'exact_central_unit_entire_I0_markers_and_nonpower_scope'),
 'odd_low_degree':('repo',PREFIX+'ODD_CUBIC_HEIGHT_LOW_DEGREE_CLASSIFICATION.md','ONLY_unit_center_degree_at_most5_and_all_order_rational_u_endpoint_reconstruction_NOT_old_sharpness_witness'),
 'odd_exact_certificate':('repo',PREFIX+'odd_cubic_height_low_degree_certificates.json','ACTUAL_four_Bezout_certificates_with_moment_generators_and_scalings_recomputed'),
 'odd_seven':('repo',PREFIX+'ODD_CUBIC_DEGREE_SEVEN_SAFETY.md','ONLY_exact_degree7_unit_center_exclusion_and_pole_root_cases_with_new_complete_curve_separation'),
 'one_pole':('repo',PREFIX+'ONE_POLE_ODD_CUBIC_CLASSIFICATION.md','ONLY_constant_coefficient_boundary_one_geometric_root_degrees_and_N1_fixed_zero_or_inversion_branches'),
 'mixed425':('repo',PREFIX+'MIXED_CUBIC_ZERO_CENTER_AND_CONDITIONAL_OFFCENTER.md','ONLY_Section2_algebraic_first4_and_Sections3--4_E0unit_C0_plus_or_minus3B0_I0_branch_NOT_FPYZ_or_all_height'),
 'milne_galois':('repo',PREFIX+'sources/milne_fields_galois_v510.pdf','ACTUAL_Theorem3_17_printed39_finite_Galois_correspondence_used_in_imported_one_pole_cyclic_step_NOT_Remark5_47_norm_branch')}
REVIEWED_RAW={
 'odd_low_degree':'51e8ccec1e80166c1c70f87e5cc83b0725a8d3bb04eb697931da18b11a06a826',
 'odd_exact_certificate':'78341369e8dc6135c9af33a3e32bf746abd49c2174a3cdab7bcee873f9bfa0b5',
 'odd_seven':'af8611c27bd832dedcbb77587fedd604d811c93f2ff79a210937bbf6ab03a958',
 'one_pole':'9cc0bd3906988d71d7b47211949823df1c4696183c28e69cf8d2299141db7bf6',
 'mixed425':'c4fee78a3b99b286e22818f63c90f8da450293bde974b354fb22bd2866bf1dc6',
 'milne_galois':'5c43ea0bf4ec190b819727fe720414a80e09f9e762533d08bad5b8a0b6de7273'}
SCOPE=dict(schema=1,reviewed_proof_sha256=PROOF_SHA,characteristic_zero=True,
 regularity=['C0+D0=0','B0+E0=0'],
 invariant='I=27*C*E^2-9*B*D*E+2*D^3',
 central_unit_hypotheses=['I=0','E0!=0','C0=D0=0','B0=-E0'],
 central_classification_max_least_degree=7,safe_central_exact_degrees=[2,5],
 unsafe_central_exact_degrees=[6],odd_degree6_unique='H=b*(1+a*x)^2,E=-b*(1+a*x)^4,a*b!=0',
 rational_depression_polynomial_H=True,depression_does_not_increase_least_degree=True,
 exact_moment_transfer='EVERY m>=1 and EVERY polynomial-x monomial marker x^ell,ell>=0',
 arbitrary_marker_rational_shift_safety_transfer=False,
 entire_invariant_zero_degree5_and7_safe=True,E0_zero_case_requires_I_zero=True,
 E0_zero_case_all_center_coefficients_zero=True,full_zero_center_uses='425 algebraic first FOUR moments only',
 nonpower_exact_J_height=3,nonpower_negative_center='(C0 != 0) OR (B0 != 0)',
 nonpower_over_every_extension=True,pure_unsafe_height3_nonpower=True,
 original_lift_power_argument_uses_both_gradings=True,
 pure_phase_forced_I_zero=False,generic_I_nonzero_classified=False,
 unrestricted_global_degrees5_or7_resolved=False,prime11_minimality_claimed=False,
 rational_u_inferred_from_finite_controls=False,full_shape430_needed=False,
 FPYZ_needed=False,all_Aut_or_424_needed=False,
 controls=dict(raw_m_max=6,polynomial_x_markers=[0,1,2,3],degree6_boxes=6,
 retained_Bezout_cases=['h1','h3_j1','h0_d_nonzero','zero_center_core']),
 assurance_flags=ASSURANCE)
def need(v,msg):
 if not v:raise ValueError(msg)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):
 def unique(pairs):
  out={}
  for k,v in pairs:
   need(k not in out,'Duplicate JSON key');out[k]=v
  return out
 return json.loads(p.read_text(),object_pairs_hook=unique)
def typed(a,b):
 if type(a)is not type(b):return False
 if isinstance(b,dict):return set(a)==set(b) and all(typed(a[k],v) for k,v in b.items())
 if isinstance(b,list):return len(a)==len(b) and all(typed(x,y) for x,y in zip(a,b))
 return a==b
def inputs(v):need(typed(v,SCOPE),'Exact typed invariant-zero/height3 scope')
def sources(v,repo):
 need(isinstance(v,dict) and set(v)=={'schema','sources','assurance_flags'},'Exact source fields')
 need(type(v['schema'])is int and v['schema']==1 and typed(v['assurance_flags'],ASSURANCE),'Typed schema/assurance')
 rows=v['sources'];need(type(rows)is list and len(rows)==len(SPECS),'Every exact mandatory role')
 out={}
 for row in rows:
  need(type(row)is dict and set(row)=={'role','kind','path','relation','sha256'},'Exact source row')
  role=row['role'];need(type(role)is str and role in SPECS and role not in out,'Unique known role')
  need((row['kind'],row['path'],row['relation'])==SPECS[role],'Exact role/relation/path')
  root=HERE if row['kind']=='local' else repo;p=(root/row['path']).resolve()
  need(p.is_relative_to(root) and p.is_file() and sha(p)==row['sha256'],'Missing/stale source:'+role)
  out[role]=row
 need(out['written_proof']['sha256']==PROOF_SHA,'Exact reviewed proof')
 for role,h in REVIEWED_RAW.items():need(out[role]['sha256']==h,'Reviewed raw mathematical body:'+role)
 return out
def add(*ps):
 r={}
 for p in ps:
  for e,a in p.items():r[e]=r.get(e,0)+a
 return {e:a for e,a in r.items() if a}
def mul(p,q):
 r={}
 for e,a in p.items():
  for f,b in q.items():
   k=tuple(x+y for x,y in zip(e,f));r[k]=r.get(k,0)+a*b
 return {e:a for e,a in r.items() if a}
def scale(p,c):return {e:a*c for e,a in p.items() if a*c}
def power(p,k,n=2):
 r={(0,)*n:Q(1)}
 for _ in range(k):r=mul(r,p)
 return r
def avg(p,l=0):
 return sum((Q(a,j+1) for(i,j),a in p.items() if i+l==0 and j%2==0),Q(0))
def divide_x(p):return {(i-1,j):a for(i,j),a in p.items()}
x={(1,0):Q(1)};t={(0,1):Q(1)};one={(0,0):Q(1)}
def phase_type2(b,r,s):
 c={(1,0):Q(r),(2,0):Q(s)}
 z=add(t,scale(c,-1));return divide_x(scale(add(z,scale(power(z,3),-1)),Q(b)))
def phase_type3(a,b,p0,p1):
 U=add(one,scale(x,Q(a)))
 L=add(scale(x,Q(a)*Q(p0)+Q(p1)),scale(power(x,2),Q(a)*Q(p1)))
 C=scale(mul(U,add(power(L,3),scale(L,-1))),Q(b))
 B=scale(mul(power(U,2),add(one,scale(power(L,2),-3))),Q(b))
 D=scale(mul(power(U,3),L),3*Q(b));E=scale(power(U,4),-Q(b))
 P=add(C,mul(B,t),mul(D,power(t,2)),mul(E,power(t,3)))
 H=scale(power(U,2),Q(b))
 G=add(mul(H,t),mul(E,power(t,3)))
 return divide_x(P),divide_x(G)

def central_controls():
 # Generic polynomial identities in variables E,H,c,t.
 n=4
 def v(i):
  e=[0]*n;e[i]=1;return {tuple(e):Q(1)}
 E,H,c,T=(v(i) for i in range(4))
 C=add(scale(mul(E,power(c,3,n)),-1),scale(mul(H,c),-1))
 B=add(H,scale(mul(E,power(c,2,n)),3));D=scale(mul(E,c),-3)
 I=add(scale(mul(C,power(E,2,n)),27),scale(mul(mul(B,D),E),-9),scale(power(D,3,n),2))
 assert not I
 P=add(C,mul(B,T),mul(D,power(T,2,n)),mul(E,power(T,3,n)))
 z=add(T,scale(c,-1))
 assert P==add(mul(E,power(z,3,n)),mul(H,z))
 assert add(scale(mul(E,B),3),scale(power(D,2,n),-1))==scale(mul(E,H),3)
 raw=[];moment_points=0
 # Polynomial-shear safe models and genuine rational-center unsafe models.
 cases=[('safe',1,1,1),('safe',Q(3,2),0,1),
  ('unsafe',1,1,1,0),('unsafe',1,1,0,1),('unsafe',-2,Q(3,2),2,-3)]
 for case in cases:
  if case[0]=='safe':
   _,b,r,s=case;f=phase_type2(b,r,s)
   g=divide_x(scale(add(t,scale(power(t,3),-1)),Q(b)))
  else:
   _,a,b,p0,p1=case;f,g=phase_type3(a,b,p0,p1)
  fp=one;gp=one;rows=[]
  for m in range(1,SCOPE['controls']['raw_m_max']+1):
   fp=mul(fp,f);gp=mul(gp,g)
   values=[]
   for ell in SCOPE['controls']['polynomial_x_markers']:
    a1,a2=avg(fp,ell),avg(gp,ell);assert a1==a2
    values.append(str(a1));moment_points+=1
   assert avg(fp)==0
   if case[0]=='unsafe' and m%2==0:
    k=sum((Q((-1)**r*comb(m,r),m+2*r+1) for r in range(m+1)),Q(0))
    assert avg(fp,1)==Q(b)**m*Q(a)**(m-1)*k and avg(fp,1)!=0
   rows.append(dict(m=m,markers_x0_to_x3=values))
  degree=max(i+j for i,j in f)
  assert degree==(5 if case[0]=='safe' else 6)
  raw.append(dict(parameters=[str(v) for v in case],least_degree=degree,rows=rows))
 boxes=[(c,e) for c in range(7) for e in range(5) if c%2==e%2 and max(c,e+2)==6]
 assert boxes==[(0,4),(2,4),(4,4),(6,0),(6,2),(6,4)]
 # A double center pole at U=0 cannot be silently accepted.
 U=add(one,x);L=add(one,scale(power(U,2),-1))
 # c=1/U^2-1; U^2*C = U^-? cleared exact numerator below.
 badC_numerator=add(power(L,3),scale(mul(L,power(U,2)),-1))
 value_at_minus1=sum(a*(-1)**i for(i,j),a in badC_numerator.items())
 assert value_at_minus1==1

 return dict(universal_identities=['depressed_invariant_zero','literal_shift_identity','polynomial_H_cleared_identity'],degree6_boxes=boxes,raw_parameter_cases=raw,exact_xmarker_moment_comparisons=moment_points,double_pole_negative=True)

def cert_clean(p):
    return {k: Q(v) for k, v in p.items() if v}

def cert_add(p, q, scale=1):
    p = dict(p)
    for k, v in q.items():
        p[k] = p.get(k, Q(0)) + scale * v
    return cert_clean(p)

def cert_scale(p, v):
    return cert_clean({k: v * c for k, c in p.items()})

def cert_mul(p, q):
    out = {}
    for a, v in p.items():
        for b, w in q.items():
            k = tuple((i + j for i, j in zip(a, b)))
            out[k] = out.get(k, Q(0)) + v * w
    return cert_clean(out)

def cert_power(p, n):
    out = {(0,) * len(next(iter(p))): Q(1)}
    for _ in range(n):
        out = cert_mul(out, p)
    return out

def cert_evaluate(p, values):
    return sum((v * cert_prod((a ** k for a, k in zip(values, powers))) for powers, v in p.items()), Q(0))

def cert_prod(values):
    out = 1
    for v in values:
        out *= v
    return out

def cert_decode(rows):
    return cert_clean({tuple(r['powers']): Q(r['coefficient']) for r in rows})

def cert_decode_uni(rows):
    return cert_clean({(r['power'],): Q(r['coefficient']) for r in rows})

def cert_h01_moment(m, linear_H):
    out = {}
    N = 3 * m + 1
    for h in range(m + 1) if linear_H else [0]:
        for j in range(m - h + 1):
            for k in range((m - h - j) // 2 + 1):
                ell = m - h - j - 2 * k
                mon = (j, k, ell)
                v = comb(m, h) * comb(N, j) * comb(N - j, k) * (-1) ** ell * comb(m + ell, ell)
                out[mon] = out.get(mon, 0) + v
    return cert_clean(out)

def cert_h3_moment(m, linear_J):
    out = {}
    for j in range(m + 1):
        for k in range((m - j) // 2 + 1):
            for ell in range((m - j - 2 * k) // 3 + 1):
                if j + k + ell > m:
                    continue
                r = m - j - 2 * k - 3 * ell
                if not linear_J and r:
                    continue
                v = comb(m, j) * comb(m - j, k) * comb(m - j - k, ell) * (comb(3 * m + 1, r) if linear_J else 1)
                mon = (j, k, ell)
                out[mon] = out.get(mon, 0) + v
    return cert_clean(out)

def cert_solve_c(f2):
    alpha = f2[0, 1, 0]
    assert all((j == 0 or (i, j, k) == (0, 1, 0) for i, j, k in f2))
    return cert_clean({(i, k): -v / alpha for (i, j, k), v in f2.items() if not j})

def cert_substitute_c(p, cstar):
    out = {}
    for (i, j, k), v in p.items():
        out = cert_add(out, cert_mul({(i, k): v}, cert_power(cstar, j)))
    return out

def cert_set_d_one(p):
    out = {}
    for (i, j), v in p.items():
        out[i,] = out.get((i,), Q(0)) + v
    return cert_clean(out)

def cert_check_unit_certificate(case, fn):
    cstar = cert_solve_c(fn(2))
    generators = []
    for m, factor, encoded in zip(case['generator_moments'], case['generator_scalars_relative_to_substituted_moment'], case['primitive_generators']):
        got = cert_scale(cert_substitute_c(fn(m), cstar), Q(factor))
        assert got == cert_decode(encoded)
        generators.append(got)
    multipliers = list(map(cert_decode, case['multipliers']))
    total = {}
    for a, g in zip(multipliers, generators):
        total = cert_add(total, cert_mul(a, g))
    assert total == {(0, 0): Q(1)}
    altered = multipliers[0].copy()
    altered[0, 0] = altered.get((0, 0), Q(0)) + 1
    bad = cert_mul(altered, generators[0])
    for a, g in zip(multipliers[1:], generators[1:]):
        bad = cert_add(bad, cert_mul(a, g))
    assert bad != {(0, 0): Q(1)}
    return {'degree_bound': case['degree_bound'], 'multiplier_terms': list(map(len, multipliers)), 'identity': 'exactly1', 'corrupted_multiplier_rejected': True}

def cert_check_univariate_certificate(case, expected):
    gens = list(map(cert_decode_uni, case['generators']))
    multipliers = list(map(cert_decode_uni, case['multipliers']))
    assert gens == expected
    total = {}
    for a, g in zip(multipliers, gens):
        total = cert_add(total, cert_mul(a, g))
    assert total == {(0,): Q(1)}
    return {'generator_degrees': [max((m[0] for m in g)) for g in gens], 'multiplier_degrees': [max((m[0] for m in a)) for a in multipliers], 'identity': 'exactly1'}

def cert_h0_controls(case):
    fn = lambda m: cert_h01_moment(m, False)
    cstar = cert_solve_c(fn(2))
    g4 = cert_substitute_c(fn(4), cstar)
    g6 = cert_substitute_c(fn(6), cstar)
    result = cert_check_univariate_certificate(case, [cert_set_d_one(g4), cert_set_d_one(g6)])
    assert {(i,): v for (i, j), v in g4.items() if not j} == {(4,): -Q(1157)}
    result['d_zero_restriction'] = '-1157*b^4'
    return result

def cert_constant_J_controls():
    fn = lambda m: cert_h3_moment(m, False)
    cstar = cert_solve_c(fn(2))
    g4 = cert_substitute_c(fn(4), cstar)
    g6 = cert_substitute_c(fn(6), cstar)
    assert g4 == {(4, 0): -Q(7, 2), (1, 1): Q(12)}
    assert {(j,): v for (i, j), v in g6.items() if not i} == {(2,): Q(15)}
    assert cert_evaluate(g6, [Q(1), Q(7, 24)]) == Q(1397, 192)
    return {'first_two_equations': 'b=-a^2/2; a*(24c-7a^3)/2=0', 'a_zero_sixth_equation': '15c^2=0', 'a_nonzero_normalized_sixth': str(Q(1397, 192))}

def cert_zero_center_controls(case):
    expected = []
    for m in [2, 4]:
        expected.append({(m - j,): Q(comb(m, j), m + 2 * j + 1) for j in range(m + 1)})
    return cert_check_univariate_certificate(case, expected)

def cert_unary_mul(a, b, order):
    out = [Q(0)] * min(len(a) + len(b) - 1, order + 1)
    for i, v in enumerate(a):
        for j, w in enumerate(b):
            if i + j < len(out):
                out[i + j] += v * w
    return out

def cert_unary_power(a, n, order):
    out = [Q(1)]
    for _ in range(n):
        out = cert_unary_mul(out, a, order)
    return out

def cert_raw_even_moment(C, D, m):
    return sum((Q(comb(m, j), m + 2 * j + 1) * (cert_unary_mul(cert_unary_power(C, m - j, m), cert_unary_power(D, j, m), m) + [Q(0)] * (m + 1))[m] for j in range(m + 1)), Q(0))

def cert_kappa(m):
    if m % 2:
        return Q(0)
    return sum((Q((-1) ** j * comb(m, j), m + 2 * j + 1) for j in range(m + 1)), Q(0))

def cert_raw_controls():
    count = 0
    for linearH in [False, True]:
        for b, c, d in [(0, 0, 0), (1, 2, 3), (-2, 1, 2), (Q(1, 2), -1, Q(2, 3))]:
            H = [1, 1] if linearH else [1]
            J = [1, b, c]
            K = [1, d]
            C = cert_unary_mul(H, cert_unary_power(J, 2, 5), 5)
            D = [-v for v in cert_unary_mul(H, cert_unary_power(K, 2, 3), 3)]
            for m in [2, 4, 6, 8]:
                assert cert_raw_even_moment(C, D, m) == cert_kappa(m) * cert_evaluate(cert_h01_moment(m, linearH), [b, c, d])
                count += 1
    for b, c, d in [(0, 0, 0), (1, 2, 3), (-2, 1, 2), (Q(1, 2), -1, Q(2, 3))]:
        H = [1, b, c, d]
        C = cert_unary_mul(H, [1, 2, 1], 5)
        D = [-v for v in H]
        for m in [2, 4, 6, 8]:
            assert cert_raw_even_moment(C, D, m) == cert_kappa(m) * cert_evaluate(cert_h3_moment(m, True), [b, c, d])
            count += 1
    return count

def invariant_and_grading_controls():
 n=2;C={(1,0):Q(1)};B={(0,1):Q(1)}
 D=scale(C,-1);E=scale(B,-1)
 value=add(scale(mul(C,power(E,2,n)),27),scale(mul(mul(B,D),E),-9),scale(power(D,3,n),2))
 need(value==scale(mul(C,add(scale(power(B,2,n),9),scale(power(C,2,n),-1))),2),
      'Universal regular-center invariant factorization')
 n=3;C={(1,0,0):Q(1)};B={(0,1,0):Q(1)};T={(0,0,1):Q(1)}
 negative=mul(add({(0,0,0):Q(1)},scale(power(T,2,n),-1)),add(C,mul(B,T)))
 need(negative=={(1,0,0):1,(0,1,1):1,(1,0,2):-1,(0,1,3):-1},
      'Negative center vanishes iff BOTH center coefficients vanish')
 # E0=B0=0 and C0=-D0 gives I0=2*D0^3, not automatic zero.
 n=1;D={(1,):Q(1)}
 need(scale(power(D,3,n),2)=={(3,):Q(2)},'E0zero still requires I0')
 # Exact-height3 does not alone imply nonpower: t^3 is a literal cube.
 T={(0,1):Q(1)};cube=power(T,3)
 need(cube=={(0,3):1},'Negative-center hypothesis cannot be omitted')
 # In FULL original variables H=xi1*(1+z1); H^3 has one dual degree
 # but FOUR coordinate degrees. One grading alone cannot descend a power.
 H={(1,0,0,0):Q(1),(1,0,1,0):Q(1)};H3=power(H,3,4)
 dual=sorted({e[0]+e[1] for e in H3});coord=sorted({e[2]+e[3] for e in H3})
 need(dual==[3] and coord==[0,1,2,3],'Both original gradings are indispensable')
 return dict(universal=['regular_center_I_factorization','negative_center_joint_nonzero_coefficients','E0zero_I0_constraint'],
   actual_negative_controls=['t^3 is exact height3 but a cube with zero center',
      'E0=B0=0,D0=1,C0=-1 gives I0=2',
      'xi1*(1+z1) cubed is dual-homogeneous but not coordinate-homogeneous'],
   full_original_dual_degrees=dual,full_original_coordinate_degrees=coord)

def retained_certificates(data):
 need(type(data)is dict and set(SCOPE['controls']['retained_Bezout_cases'])<=set(data),'Actual retained certificate cases')
 out=dict(
  h1=cert_check_unit_certificate(data['h1'],lambda m:cert_h01_moment(m,True)),
  h3_j1=cert_check_unit_certificate(data['h3_j1'],lambda m:cert_h3_moment(m,True)),
  h0_d_nonzero=cert_h0_controls(data['h0_d_nonzero']),
  zero_center_core=cert_zero_center_controls(data['zero_center_core']),
  constant_J=cert_constant_J_controls(),
  raw_interval_generator_comparisons=cert_raw_controls())
 need(len(SCOPE['controls']['retained_Bezout_cases'])==4,'Four actual Bezout identities')
 return out

def negatives(ledger,repo):
 out=[]
 def reject(label,fn):
  try:fn()
  except (ValueError,KeyError,TypeError,OSError):out.append(label);return
  raise ValueError('Negative accepted:'+label)
 for role in SPECS:
  a=copy.deepcopy(ledger);a['sources']=[r for r in a['sources'] if r['role']!=role]
  reject('missing_role:'+role,lambda a=a:sources(a,repo))
  a=copy.deepcopy(ledger)
  next(r for r in a['sources'] if r['role']==role)['sha256']='0'*64
  reject('bad_hash:'+role,lambda a=a:sources(a,repo))
 a=copy.deepcopy(ledger);a['schema']=True
 reject('Boolean_ledger_schema',lambda:sources(a,repo))
 a=copy.deepcopy(ledger);a['sources'].append(copy.deepcopy(a['sources'][0]))
 reject('duplicate_source_role',lambda:sources(a,repo))
 changes={
  'schema':True,'central_classification_max_least_degree':True,
  'nonpower_exact_J_height':2,'nonpower_negative_center':'(C0 != 0) AND (B0 != 0)',
  'E0_zero_case_requires_I_zero':False,
  'arbitrary_marker_rational_shift_safety_transfer':True,
  'pure_phase_forced_I_zero':True,'generic_I_nonzero_classified':True,
  'unrestricted_global_degrees5_or7_resolved':True,
  'prime11_minimality_claimed':True,
  'original_lift_power_argument_uses_both_gradings':False,
  'rational_u_inferred_from_finite_controls':True,
  'unsafe_central_exact_degrees':[5,6,7],
  'full_shape430_needed':True,'FPYZ_needed':True,'all_Aut_or_424_needed':True}
 for key,value in changes.items():
  a=copy.deepcopy(SCOPE);a[key]=value
  reject('false_scope:'+key,lambda a=a:inputs(a))
 a=copy.deepcopy(SCOPE);a['assurance_flags']['external_review']=True
 reject('unsupported_assurance',lambda:inputs(a))
 return out

def main():
 ap=argparse.ArgumentParser(description=__doc__)
 ap.add_argument('--repo-root',type=Path,required=True)
 ap.add_argument('--output',type=Path,default=HERE/RN)
 args=ap.parse_args();repo=args.repo_root.resolve()
 need(repo.is_dir(),'Repository root exists')
 inp=load(HERE/IN);inputs(inp);ledger=load(HERE/LN);pins=sources(ledger,repo)
 controls=central_controls()
 additional=invariant_and_grading_controls()
 cert=load(repo/SPECS['odd_exact_certificate'][1])
 certificates=retained_certificates(cert)
 controls_negative=negatives(ledger,repo)
 report=dict(schema=1,status='PASS',proof_sha256=PROOF_SHA,checker_sha256=sha(Path(__file__)),
  inputs_sha256=sha(HERE/IN),source_ledger_sha256=sha(HERE/LN),scope=inp,
  sources=[pins[k] for k in sorted(pins)],central_controls=controls,
  additional_universal_and_boundary_controls=additional,retained_certificates=certificates,
  source_scope_negative_controls=controls_negative,
  assurance_flags=ASSURANCE,
  meaning='Exact bounded controls and source authentication support the written all-order proof. Neither rational-u reconstruction, complete-curve separation nor all-order classification is inferred from a finite moment prefix.')
 args.output.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
 print(json.dumps(dict(status='PASS',mandatory_roles=len(pins),universal=6,
  Bezout_identities=4,raw_moment_comparisons=controls['exact_xmarker_moment_comparisons'],
  certificate_generator_comparisons=certificates['raw_interval_generator_comparisons'],
  degree6_boxes=len(controls['degree6_boxes']),source_scope_controls=len(controls_negative),
  receipt_sha256=sha(args.output))))
if __name__=='__main__':main()

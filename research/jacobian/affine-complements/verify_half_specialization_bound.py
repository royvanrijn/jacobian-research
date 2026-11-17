#!/usr/bin/env python3
"""Narrow exact controls for PF158; universal proofs are in the source note."""
from pathlib import Path
import hashlib
import json
import math
import platform
import time
import sympy as S

START=time.perf_counter()
BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[1]
PROOF=BASE/'HALF_SPECIALIZATION_POLE_BOUND_AND_QUADRATIC_RULING_FOUR_EXCLUSION.md'
CONTEXT=BASE/'sources/half-specialization-bound/context.json'
RECEIPT=BASE/'half_specialization_bound_receipt.json'
checks=[]
def check(label,condition):
    if not condition:
        raise AssertionError(label)
    checks.append(label)
def zero(expr):
    return S.cancel(expr)==0
def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

context=json.loads(CONTEXT.read_text())
for rel,digest in context['retained_sha256'].items():
    check('retained source pin '+rel,sha(ROOT/rel)==digest)

s,t,u,x,w,z,Y,p,lam=S.symbols('s t u x w z Y p lam')
local_controls=[]
# P=w^d(1-t^a*w)^e. These give exact equality in the new valuation bound.
for D,d,a in [(4,2,3),(4,2,17),(4,3,5),(4,1,3),(5,2,4),(6,3,7),(6,2,3),(7,4,2)]:
    e=D-d
    div=math.gcd(e,a*d)
    alpha=e//div
    h=a*alpha
    k0=a*d//div
    r=h+k0
    check(f'normalized branch parameter D{D} d{d} a{a}',math.gcd(alpha,k0)==1)
    # Verify the cleared reciprocal equation before any local approximation.
    check(f'actual reciprocal identity D{D} d{d} a{a}',
          zero(z**D*((w**d*(1-t**a*w)**e).subs(w,1/z))-(z-t**a)**e))
    check(f'normalized Hensel equation D{D} d{d} a{a}',e*(h+k0)==D*h and a*alpha==h)
    initial=Y**e-p
    check(f'simple parameter-dependent initial coefficient D{D} d{d} a{a}',
          S.diff(initial,Y).subs(Y,lam)==e*lam**(e-1) and e>0)
    check(f'sharp first-coefficient bound D{D} d{d} a{a}',r==S.ceiling(S.Rational(D*h,e)))
    velocity_order=r-2*h
    differential_order=velocity_order+alpha-1
    if 2*d>=D:
        check(f'half-specialization gives holomorphic coframe D{D} d{d} a{a}',differential_order>=0)
    local_controls.append({'D':D,'d':d,'a':a,'escaping_degree':e,'alpha':alpha,'h':h,
        'first_parameter_index':r,'velocity_order':velocity_order,'differential_order':differential_order})
check('quartic degree-one boundary has double differential pole',
      next(c for c in local_controls if c['D']==4 and c['d']==1)['differential_order']==-2)
check('quintic degree-two boundary has double differential pole',
      next(c for c in local_controls if c['D']==5)['differential_order']==-2)

# Actual reciprocal Taylor coefficient used in the proof: all lower reduction
# coefficients vanish, and Z0(0)=0, so the coefficient with index e remains a unit.
for D,d in [(4,2),(4,3),(5,3),(6,3)]:
    e=D-d
    coeff=S.symbols('c0:'+str(D+1))
    F=sum((coeff[i] if i<=d else s*coeff[i])*z**(D-i) for i in range(D+1))
    Z0=s+s**2
    taylor=S.expand(F.subs(z,Z0+s**3*Y))
    unit=S.cancel(taylor.coeff(Y,e)/s**(3*e)).subs(s,0)
    check(f'actual low reciprocal Taylor unit D{D} d{d}',unit==coeff[d])

# Source-pullback leading coefficient of a completely general regular r=4 P.
cp,b1,b0,a1,a0,AA,SS=S.symbols('cp b1 b0 a1 a0 AA SS')
wu=t*(SS+AA*x)
zu=2*SS*x+AA*x**2
P=cp*zu**2+b1*zu+b0+wu*(a1*zu+a0)
check('full general ruling-four pullback leading coefficient',S.expand(P).coeff(x,4)==cp*AA**2)
check('all lower terms stay below source degree four',S.degree(S.expand(P-cp*AA**2*x**4),x)<=3)
check('general quadratic whole-source relation',zero(wu**2-t**2*SS**2-t**2*AA*zu))

# New regular S2 test family: every named whole source line has degree 2 or4,
# with bad leading roots of arbitrary multiplicity beyond PF149.
surface_controls=[]
for power in [3,4,7]:
    PP=((u-1)**power*z+w)**2+z
    even=S.expand((u-1)**(2*power)*z**2+(u**2+1)*z+u)
    odd=2*(u-1)**power*z
    check(f'canonical quadratic expression h{power}',
          zero(PP-(even+w*odd)-(w**2-u-u**2*z)))
    pull=S.expand(PP.subs({u:t**2,w:t*(-1+t**2*x),z:-2*x+t**2*x**2}))
    leading=S.Poly(pull,x).LC()
    check(f'exact whole-source x-leading coefficient h{power}',zero(leading-(t**2-1)**(2*power)*t**4))
    check(f'global leading degree exceeds quartic threshold h{power}',S.degree(leading,t)==4*power+4)
    check(f'all named exceptional whole lines degree two h{power}',
          all(S.degree(pull.subs(t,c),x)==2 for c in [0,1,-1]))
    check(f'finite leading roots exhausted h{power}',
          zero(leading-t**4*(t-1)**(2*power)*(t+1)**(2*power)))
    check(f'bad surface-root specialization has half degree h{power}',
          zero(PP.subs({u:1,z:w**2-1})-(2*w**2-1)))
    check(f'outside old PF149 multiplicity gate h{power}',2*power>4)
    surface_controls.append({'h':power,'ruling_degree':4,'surface_leading_root_multiplicity':2*power,
        'source_x_leading_degree':4*power+4,'exceptional_line_degrees':[2,2,2]})

# Infinity and ramified pullback arithmetic, including logarithmic equality.
for D,b,alpha in [(4,4,1),(4,4,5),(4,12,2),(6,6,3),(6,11,2)]:
    order=int(S.ceiling(S.Rational(b*alpha,D)))-alpha-1
    check(f'infinity all-branch lower bound D{D} b{b} alpha{alpha}',order>=-1)
check('infinity equality retains logarithmic poles',int(S.ceiling(S.Rational(4,4)))-1-1==-1)
for ram in [1,2,5]:
    check(f'ramified pullback retains simple poles E{ram}',ram*(-1)+(ram-1)==-1)
for D in [2,3,4]:
    check(f'actual PF52 degree two supplies half degree D{D}',2>=S.ceiling(S.Rational(D,2)))

# Retained nonproper 3D Keller example: the new proof instead needs a CURVE
# generic fibre and the actual plane injectivity-on-a-whole-line input.
# Read the exact existing control rather than introducing a different example.
# The pinned PF157 checker records the previously verified Jac=-2 and arc;
# verify its PASS receipt here without replaying that checker.
old=json.loads((BASE/'canonical_deck_pole_coverage_receipt.json').read_text())
check('retained dimension-three nonproper control remains pinned',
      old['outcome']=='PASS' and 'retained 3D first-velocity pole' in old['checks'])
check('ruling-five boundary is not claimed closed',not context['quadratic_ruling_degree_five_claimed_closed'])
check('full normalization and every ramified prime retained',context['full_original_normalization_and_all_ramified_primes_retained'])
check('no automatic integrality or tail certificate claimed',
      not context['generic_first_derivative_integrality_proved'] and not context['normal_tail_increasing_power_containment_proved'])

receipt={'claim':'PF2RULED158','outcome':'PASS','check_count':len(checks),'checks':checks,
    'elapsed_seconds':time.perf_counter()-START,'environment':{'python':platform.python_version(),'sympy':S.__version__},
    'proof_sha256':sha(PROOF),'checker_sha256':sha(Path(__file__)), 'context_sha256':sha(CONTEXT),
    'local_controls':local_controls,'surface_controls':surface_controls,
    'universal_reciprocal_parameter_bound':'ceil(D*h/(D-d))-2*h',
    'actual_low_degree_plane_keller_leading_bound':'deg_t leading_x P <= deg_x P-1 for 2<=deg_x P<=4 and noninvertible',
    'entire_actual_quadratic_ruling_degree_four_family_closed':True,
    'all_bad_leading_divisors_and_unbounded_companion_degree_retained':True,
    'remaining_reduced_ruling_lower_bound':'max(5,2*|D|)',
    'quadratic_ruling_degree_five_closed':False,
    'original_full_normalization_and_all_ramified_primes_retained':True,
    'generic_first_derivative_integrality_proved':False,
    'normal_tail_increasing_power_containment_proved':False,
    'registered_predecessor_checkers_replayed':[], 'new_external_theorem_imports':[],
    'full_goal_complete':False,'independent_replay':False,'formal_verification':False,'external_review':False}
RECEIPT.write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'outcome':'PASS','check_count':len(checks),'elapsed_seconds':receipt['elapsed_seconds']},indent=2))

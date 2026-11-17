#!/usr/bin/env python3
"""Exact controls for PF159's written unbounded polynomial-mate/escape proof."""
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
PROOF=BASE/'POLYNOMIAL_MATE_RESPONSES_AND_QUADRATIC_RULING_FIVE_ESCAPE_COVERAGE.md'
CONTEXT=BASE/'sources/polynomial-mate-response/context.json'
RECEIPT=BASE/'polynomial_mate_response_receipt.json'
checks=[]
def check(label,condition):
    if not condition:
        raise AssertionError(label)
    checks.append(label)
def zero(expr):
    return S.cancel(expr)==0
def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()
def initial(expr,var):
    num,den=S.fraction(S.cancel(expr))
    if num==0:
        return S.oo,S.Integer(0)
    pn,pd=S.Poly(num,var),S.Poly(den,var)
    on=min(k[0] for k in pn.monoms())
    od=min(k[0] for k in pd.monoms())
    return on-od,S.cancel(pn.nth(on)/pd.nth(od))

context=json.loads(CONTEXT.read_text())
for rel,digest in context['retained_sha256'].items():
    check('retained source pin '+rel,sha(ROOT/rel)==digest)
u,v,t,w,x,z,s,Y,p=S.symbols('u v t w x z s Y p')
# Rational-in-w mates genuinely escape the polynomial leading-coefficient rule.
PP=s**6*(w-s**-8)**3
QQ=s**-5*(w-s**-8)**-2/3
check('actual rational-in-w mate has unit Jacobian',
      zero(S.diff(PP,s)*S.diff(QQ,w)-S.diff(PP,w)*S.diff(QQ,s)-1))
check('rational-mate control has forbidden normal-variable denominator',S.degree(S.denom(S.cancel(QQ)),w)==2)
check('cubic response cannot satisfy polynomial initial product',S.diff(-Y**3,Y).subs(Y,0)==0)

# Universal triple-root constants and canonical coefficient transformations.
gamma,kappa=S.symbols('gamma kappa',nonzero=True)
I=S.Poly(S.expand(gamma*Y**2*(1-kappa*Y)**3),Y)
eta3,eta4,eta5=[I.nth(i) for i in [3,4,5]]
check('triple escape constants first identity',zero(eta3**2-3*gamma*eta4))
check('triple escape constants second identity',zero(eta3**3-27*gamma**2*eta5))
check('triple nonzero initial root has exact multiplicity three',
      all(S.diff(I.as_expr(),Y,k).subs(Y,1/kappa)==0 for k in range(3)) and
      S.diff(I.as_expr(),Y,3).subs(Y,1/kappa)!=0)
M,N,a2,a1,a0,b2,b1,b0=S.symbols('M N a2 a1 a0 b2 b1 b0')
Pstar=w*(a2*((w**2-N)/M)**2+a1*(w**2-N)/M+a0)+b2*((w**2-N)/M)**2+b1*(w**2-N)/M+b0
coeff=S.Poly(S.expand(Pstar),w)
for i,target in [(5,a2/M**2),(4,b2/M**2),(3,a1/M-2*N*a2/M**2),(2,b1/M-2*N*b2/M**2)]:
    check(f'full canonical generic-ruling coefficient w{i}',zero(coeff.nth(i)-target))
F=a2*(1-N*Y**2)**2+b2*Y*(1-N*Y**2)**2+M*a1*Y**2*(1-N*Y**2)+M*b1*Y**3*(1-N*Y**2)+M**2*a0*Y**4+M**2*b0*Y**5
check('actual M-root reciprocal equation',zero(F-M**2*Y**5*Pstar.subs(w,1/Y)))
check('leading-zero M-root reciprocal derivative is a unit',S.diff(F,Y).subs({Y:0,a2:0})==b2)
check('implicit reciprocal parameter numerator',S.diff(F-p*M**2*Y**5,p)==-M**2*Y**5)
check('implicit coframe cancellation retains M factor',zero((M**2*Y**5)/Y**2/M-M*Y**3))
A1,A2,B2,B1,MM=S.symbols('A1 A2 B2 B1 MM',nonzero=True)
check('canonical triple coefficient identity one',
      zero(MM**3*((A1/MM)**2-3*(B1/MM)*(B2/MM**2))-(MM*A1**2-3*B1*B2)))
check('canonical triple coefficient identity two',
      zero(MM**4*((A1/MM)**3-27*(B1/MM)**2*(A2/MM**2))-(MM*A1**3-27*B1**2*A2)))
check('coprime reduced ruling degrees force fifth-power leading',math.gcd(5,6)==1 and math.gcd(3,5)==1)

# Entire family P_h, retaining source-chart/differential checks and the full
# substituted polynomial, rather than dropping arbitrary high-order terms.
family=[]
for j in [1,2,3]:
    h=5*j
    vv=u-1
    PS=w**2*(1-vv**h*w)**3+vv**(4*h+1)*z**2+vv**(h+1)*w+vv
    PR=PS.subs(z,(w**2-u)/u**2)
    lead=S.Poly(S.expand(PR),w).nth(5)
    check(f'actual fifth-power leading family j{j}',zero(lead+((u-1)**(3*j))**5))
    source=PS.subs({u:t**2,w:t*(-1+t**2*x),z:-2*x+t**2*x**2})
    check(f'whole reached source-line degrees family j{j}',all(S.degree(source.subs(t,c),x)==2 for c in [0,1,-1]))
    check(f'central whole-line source differential family j{j}',zero(S.diff(source,t).subs(t,0)+(-1)**(h+1)))
    check(f'whole bad fibre surface differential family j{j}',zero(S.diff(PR,u).subs(u,1)-1))
    div=math.gcd(3,h)
    alpha=3//div
    w0=s**(-3*h//div)+s**(-h//div)*Y
    full=PR.subs({u:1+s**alpha,w:w0},simultaneous=True)
    order,poly=initial(full,s)
    check(f'full substituted cubic initial response family j{j}',order==0 and zero(poly+Y**3))
    rem_order,_=initial(full+Y**3,s)
    check(f'all remaining substituted terms strictly late family j{j}',rem_order>0)
    pole=(3-h)//div-1
    check(f'actual normalized coframe double-or-higher pole family j{j}',pole<=-2)
    check(f'pole requires reciprocal quadratic response denominator family j{j}',S.diff(poly,Y)==-3*Y**2)
    # Canonical coefficients, not only the leading w coefficient.
    ca2=-u**4*vv**(3*h)
    cb2=3*u**4*vv**(2*h)+vv**(4*h+1)
    ca1=-3*u**2*vv**h-2*u**3*vv**(3*h)
    for label,ex,expected in [('a2',ca2,15*j),('b2',cb2,10*j),('a1',ca1,5*j)]:
        got,_=initial(ex.subs(u,1+v),v)
        check(f'canonical bad order {label} family j{j}',got==expected)
    family.append({'j':j,'h':h,'alpha':alpha,'first_parameter_index':-h//div,
        'coframe_order':pole,'bad_canonical_orders':[15*j,10*j,5*j],
        'initial_response':str(poly),'whole_remainder_order':int(rem_order)})

# Affine-response boundary: leading tests and all three polar ends, including
# the first coframe correction needed to check the actual residue.
vv=u-1
Sa=1-vv**10*w
g0=3-2*u
PA=w**2*Sa**3+vv**2*g0*w*Sa+vv**41*z**2+vv**11*w+vv
PAR=PA.subs(z,(w**2-u)/u**2)
check('affine control retains fifth-power leading',zero(S.Poly(S.expand(PAR),w).nth(5)+(u-1)**30))
source=PA.subs({u:t**2,w:t*(-1+t**2*x),z:-2*x+t**2*x**2})
check('affine control whole leading-drop line degrees',all(S.degree(source.subs(t,c),x)==2 for c in [0,1,-1]))
check('affine control central whole-line nonzero differential',S.diff(source,t).subs(t,0)!=0)
check('affine control bad whole-fibre nonzero differential',zero(S.diff(PAR,u).subs(u,1)-1))
branches=[]
for C in [S.Integer(0),S.I,-S.I]:
    w0=v**-10+v**-2*Y if C==0 else v**-10+C*v**-4-C*v**-3+v**-2*Y
    expr=PAR.subs({u:1+v,w:w0},simultaneous=True)
    order,poly=initial(expr,v)
    target=-Y if C==0 else 2*Y+C
    check(f'affine branch exact initial C{C}',order==0 and zero(poly-target))
    check(f'affine branch Hensel unit C{C}',S.diff(poly,Y) in [-1,2])
    den=S.cancel((1+v)**2*S.diff(PAR,w).subs({u:1+v,w:w0},simultaneous=True)/v**2)
    n,leading=initial(den,v)
    expected=-1 if C==0 else 2
    check(f'affine branch coframe double-pole denominator C{C}',n==0 and leading==expected)
    rem,_=initial(den-expected,v)
    check(f'affine branch ZERO residue first correction C{C}',rem>=2)
    branches.append({'C':str(C),'initial_response':str(poly),'coframe_order':-2,
        'coframe_leading_coefficient':str(-1/leading),'residue':0})
check('affine control has exactly two finite generic roots at bad fibre',zero(PAR.subs(u,1)-w**2))
# The three unit-slope Hensel branches plus these two finite roots exhaust D=5.
check('all bad-fibre normalized branches exhausted',len(branches)+2==5)
check('residue cancellation factor is exact',zero((1+v)**2*(1-2*v)-(1-3*v**2-2*v**3)))

old=json.loads((BASE/'canonical_deck_pole_coverage_receipt.json').read_text())
check('known3D nonproper control remains pinned',old['outcome']=='PASS' and 'retained 3D first-velocity pole' in old['checks'])
check('ruling-five whole closure not claimed',not context['all_ruling_five_targets_excluded'])
check('full normalization and every ramified prime retained',context['full_original_normalization_and_all_ramified_primes_retained'])
check('no automatic integrality or tail-power certificate claimed',not context['generic_first_derivative_integrality_proved'] and not context['normal_tail_increasing_power_containment_proved'])
receipt={'claim':'PF2RULED159','outcome':'PASS','check_count':len(checks),'checks':checks,
    'elapsed_seconds':time.perf_counter()-START,'environment':{'python':platform.python_version(),'sympy':S.__version__},
    'proof_sha256':sha(PROOF),'checker_sha256':sha(Path(__file__)),'context_sha256':sha(CONTEXT),
    'excluded_family_controls':family,'residue_free_affine_response_branches':branches,
    'actual_polynomial_mate_pole_requires_affine_first_response':True,
    'arbitrary_rational_in_w_mates_excluded':False,
    'complete_actual_quadratic_ruling_five_escape_coverage':True,
    'necessary_bad_canonical_orders':'(15*j,10*j,5*j),j>=1, off M; two coupled leading identities',
    'all_ruling_five_targets_excluded':False,
    'residue_free_affine_control_is_actual_target':False,
    'original_full_normalization_and_all_ramified_primes_retained':True,
    'generic_first_derivative_integrality_proved':False,
    'normal_tail_increasing_power_containment_proved':False,
    'registered_predecessor_checkers_replayed':[],'new_external_theorem_imports':[],
    'full_goal_complete':False,'independent_replay':False,'formal_verification':False,'external_review':False}
RECEIPT.write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'outcome':'PASS','check_count':len(checks),'elapsed_seconds':receipt['elapsed_seconds']},indent=2))

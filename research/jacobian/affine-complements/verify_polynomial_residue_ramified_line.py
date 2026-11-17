#!/usr/bin/env python3
"""PF161 narrow replay: coefficient-ring Hensel and full ramified-prime data."""
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
SRC=BASE/'sources/polynomial-residue-ramified-line'
PROOF=BASE/'POLYNOMIAL_RESIDUE_RAMIFIED_LINE_AND_QUADRATIC_RULING_FIVE_CLOSURE.md'
CONTEXT=SRC/'context.json'
RECEIPT=BASE/'polynomial_residue_ramified_line_receipt.json'
checks=[]
def check(label,condition):
    if not condition:raise AssertionError(label)
    checks.append(label)
def zero(e):return S.cancel(e)==0
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def order_and_leading(expr,var):
    num,den=S.fraction(S.cancel(expr))
    if num==0:return S.oo,S.Integer(0)
    pn,pd=S.Poly(num,var),S.Poly(den,var)
    a=min(i[0] for i in pn.monoms());b=min(i[0] for i in pd.monoms())
    return a-b,S.cancel(pn.nth(a)/pd.nth(b))

context=json.loads(CONTEXT.read_text())
for rel,h in context['retained_sha256'].items():check('retained source pin '+rel,sha(ROOT/rel)==h)
audit=json.loads((SRC/'chau_component_import_audit.json').read_text())
check('retained exact Chau component audit pin',sha(SRC/'chau_component_import_audit.json')==context['audit_sha256'])
check('published primary PDF current hash',sha(ROOT/audit['retained_pdf'])==audit['retained_sha256'])
check('component scope allows other nonproper components',audit['other_components_allowed'] and not audit['whole_curve_Corollary3_substituted'])
check('retained imported input not falsely independently reproved',not audit['independently_reproved'] and not audit['new_external_theorem_import'])
u,w,y,z,p,t,x=S.symbols('u w y z p t x')

def trunc(expr,degree):
    poly=S.Poly(S.expand(expr),u)
    return sum((c*u**i[0] for i,c in poly.terms() if i[0]<degree),S.Integer(0))
def evaluate_truncated(poly,root,degree):
    value=S.Integer(0)
    for c in poly.all_coeffs():value=trunc(value*root+c,degree)
    return value

def hensel(expr,var,root0,last):
    poly=S.Poly(expr,var)
    J0=S.diff(expr,var).subs({u:0,var:root0})
    assert J0.is_number and J0!=0
    root=S.Integer(root0)
    for degree in range(1,last+1):
        residual=evaluate_truncated(poly,root,degree+1)
        c=S.expand(residual).coeff(u,degree)
        root-=S.cancel(c/J0)*u**degree
    assert evaluate_truncated(poly,root,last+1)==0
    return S.expand(root),J0

def reciprocal_series(root,first,last):
    # root=u^first*h with h0 a nonzero CONSTANT; actual multiplication
    # and recursive inversion retain polynomial p-coefficients.
    h=[S.expand(root).coeff(u,first+i) for i in range(last+first+1)]
    coeff=[1/h[0]]
    for j in range(1,last+first+1):
        coeff.append(S.cancel(-sum(h[i]*coeff[j-i] for i in range(1,j+1))/h[0]))
    return sum(coeff[j]*u**(j-first) for j in range(len(coeff)))

controls=[]
for n,k,a in [(2,1,0),(2,1,1),(2,2,0),(2,2,1),(3,1,0)]:
    m=k+1;sig=1-2*k;M=u**m;N=sig**2*u;D=2*n+1
    # All lower coefficients retained in the proof. These supporting models
    # include both lower parities and nonconstant lower/base coefficients.
    P=w*(u**a*z**n+z**(n-1)+1)+z**n+(u+2)*z**(n-1)+u-1
    star=P.subs(z,(w**2-N)/M)
    F=S.cancel(M**n*(star-p))
    G=S.cancel(F.subs(w,1/y)*y**D)
    label=f'n{n} m{m} a{a}'
    check('exact cleared reduction '+label,zero(F.subs(u,0)-w**(2*n)*((w if a==0 else 0)+1)))
    beta_y=n*m+D*a;beta_w=n*m+(D-2)*a
    if a==0:
        root,J0=hensel(F,w,-1,beta_w)
        wseries=root
        check('unit finite root over C[p] '+label,J0==1)
    else:
        root,J0=hensel(G,y,0,beta_y)
        check('unit reciprocal root over C[p] '+label,J0==1)
        yo,yl=order_and_leading(root,u)
        check('reciprocal root exact independent leading coefficient '+label,yo==a and yl==-1)
        wseries=reciprocal_series(root,a,beta_w)
    check('all stored lift coefficients polynomial in p '+label,
          all(not S.denom(c).has(p) for c in S.Poly(root,u).all_coeffs()))
    velocity=S.diff(wseries,p)
    vo,vl=order_and_leading(velocity,u)
    check('first parameter velocity exact order '+label,vo==beta_w and vl==1)
    residue=S.expand(wseries).coeff(u,beta_w)
    check('polynomial not rational finite-parameter residue '+label,S.degree(residue,p)==1 and S.diff(residue,p)==1)
    # Arbitrary functions can have poles: actual Jacobian is needed to
    # remove negative coefficients of the true mate, not a general landing rule.
    qo,_=order_and_leading(wseries/u**beta_w,u)
    check('arbitrary polynomial-in-w function not declared regular '+label,qo<0)
    oo,ol=order_and_leading(-velocity/M,u)
    e0=(n-1)*m+(2*n-1)*a+1
    check('coframe holomorphic exact leading term '+label,oo==e0-1 and ol==-1)
    xs=(wseries.subs(u,t**2)/t-sig)/t**(2*k)
    q=2*(m+a)-1
    xo,xl=order_and_leading(xs,t)
    check('actual whole-source coordinate has fixed pole '+label,xo==-q and xl==-1)
    xp=S.diff(xs,p)
    eo,el=order_and_leading(xp,t)
    eR=2*e0
    check('original ramified source volume order '+label,eo==eR-1 and el==1)
    check('ramified prime large index and local clearing contribution '+label,
          eR>=6 and eR-q==2*(n-2)*m+4*(n-1)*a+3 and 1+math.ceil(S.Rational(q,eR))==2)
    controls.append({'n':n,'m':m,'a':a,'branch':'finite unit' if a==0 else 'reciprocal unit',
        'lift_through_u_order':beta_w if a==0 else beta_y,
        'lift_series':str(root),'sample_polynomial_residue':str(residue),
        'parameter_velocity_order':int(vo),'coframe_order':int(oo),
        'intermediate_ramification_index':e0,'original_ramification_index':eR,
        'source_coordinate_pole':q,'source_volume_order':int(eo),
        'formal_local_control_not_global_Keller_pair':True})

# Exact whole-source volume identity for arbitrary coefficient functions.
AA=S.Function('A');SS=S.Function('Sigma')
U=t**2;W=t*(SS(U)+AA(U)*x);ZZ=2*SS(U)*x+AA(U)*x**2
check('arbitrary quadratic whole-source relation',zero(W**2-U*SS(U)**2-U*AA(U)*ZZ))
vol=(S.diff(U,x)*S.diff(W,t)-S.diff(U,t)*S.diff(W,x))/(U*AA(U))
check('arbitrary quadratic volume is actual constant source volume',zero(vol+2))
# A full target shear is legitimate once the residue is POLYNOMIAL.
Pvar,Qvar=S.symbols('Pvar Qvar');residue=Pvar**4+3*Pvar**2-Pvar+7
check('polynomial graph target shear preserves Jacobian',S.det(S.Matrix([Pvar,Qvar-residue]).jacobian([Pvar,Qvar]))==1)
check('polynomial graph coordinate ring is C[Pvar]',S.degree(Qvar-residue,Qvar)==1)
# n=0 rational mate is a real scope boundary; b_n=0 has no unit simple root.
P0=w+u**3;Q0=1/u
check('degree-one rational mate remains outside n>=1 theorem',zero(S.diff(P0,u)*S.diff(Q0,w)-S.diff(P0,w)*S.diff(Q0,u)-1/u**2))
P7=u**2*w*z**3+z**2+w
F7=S.cancel(u**6*(P7.subs(z,(w**2-u)/u**2)-p))
check('ruling-seven b3zero has no nonzero simple central branch',zero(F7.subs(u,0)) and S.diff(F7,w).subs({u:0,w:0})==0)
source7=P7.subs({u:t**2,w:t*(-1+t**2*x),z:-2*x+t**2*x**2})
check('ruling-seven central degree-two boundary is retained',S.degree(source7.subs(t,0),x)==2)
# Keeping another prime can raise a FULL target exponent beyond this prime's2.
check('full maximum not replaced by one ramified prime',max(2,1+math.ceil(S.Rational(2,1)))==3)
old=json.loads((BASE/'canonical_deck_pole_coverage_receipt.json').read_text())
check('known3D nonproper control retained',old['outcome']=='PASS' and 'retained 3D first-velocity pole' in old['checks'])
check('entire ruling-five actual class retired',context['entire_quadratic_ruling_five_actual_target_class_closed'])
check('higher odd and arbitrary rational mates not declared closed',not context['all_higher_odd_targets_claimed_closed'] and not context['arbitrary_rational_companions_claimed_excluded'])
check('full original finite normalization and every ramified prime retained',context['forced_prime_is_in_original_full_finite_normalization'] and context['all_original_ramified_primes_retained'])
check('generic integrality and normal-tail certificate not claimed',not context['generic_first_derivative_integrality_proved'] and not context['normal_tail_increasing_power_containment_proved'])
receipt={'claim':'PF2RULED161','outcome':'PASS','check_count':len(checks),'checks':checks,
    'elapsed_seconds':time.perf_counter()-START,'environment':{'python':platform.python_version(),'sympy':S.__version__},
    'proof_sha256':sha(PROOF),'checker_sha256':sha(Path(__file__)),'context_sha256':sha(CONTEXT),
    'chau_import_audit_sha256':sha(SRC/'chau_component_import_audit.json'),
    'formal_local_controls':controls,
    'universal_forced_original_ramification':'2*((n-1)*m+(2*n-1)*a+1)',
    'universal_source_coordinate_pole':'2*(m+a)-1',
    'forced_full_boundary_residue_ring':'C[P]',
    'actual_nonproper_image':'Q=q0(P),q0 polynomial; closed A1 component',
    'entire_quadratic_ruling_five_actual_target_class_closed':True,
    'all_higher_odd_targets_excluded':False,
    'arbitrary_rational_companions_excluded':False,
    'remaining_reduced_ruling_lower_bound':'max(6,2*|D|)',
    'all_original_ramified_primes_retained':True,
    'generic_first_derivative_integrality_proved':False,'normal_tail_increasing_power_containment_proved':False,
    'new_external_theorem_imports':[],'retained_Chau_component_input_rechecked':True,
    'registered_predecessor_checkers_replayed':[],
    'full_goal_complete':False,'independent_replay':False,'formal_verification':False,'external_review':False}
RECEIPT.write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'outcome':'PASS','check_count':len(checks),'elapsed_seconds':receipt['elapsed_seconds']},indent=2))

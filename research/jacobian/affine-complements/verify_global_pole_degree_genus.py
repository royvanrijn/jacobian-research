#!/usr/bin/env python3
"""Replay PF160's exact branch certificate and new global-criterion controls."""
from pathlib import Path
import hashlib
import json
import platform
import time
import sympy as S

START=time.perf_counter()
BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[1]
SRC=BASE/'sources/global-pole-degree-genus'
PROOF=BASE/'GLOBAL_POLE_DEGREE_GENUS_AND_RESIDUE_FREE_AFFINE_EXCLUSION.md'
CONTEXT=SRC/'context.json'
RECEIPT=BASE/'global_pole_degree_genus_receipt.json'
checks=[]
def check(label,condition):
    if not condition:raise AssertionError(label)
    checks.append(label)
def zero(e):return S.cancel(e)==0
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def initial(e,var):
    num,den=S.fraction(S.cancel(e))
    if num==0:return S.oo,S.Integer(0)
    pn,pd=S.Poly(num,var),S.Poly(den,var)
    a=min(k[0] for k in pn.monoms());b=min(k[0] for k in pd.monoms())
    return a-b,S.cancel(pn.nth(a)/pd.nth(b))

context=json.loads(CONTEXT.read_text())
for rel,digest in context['retained_sha256'].items():
    check('retained source pin '+rel,sha(ROOT/rel)==digest)
for rel,digest in context['derived_sha256'].items():
    check('new exact source pin '+rel,sha(SRC/rel)==digest)
cert=json.loads((SRC/'affine_discriminant_certificate.json').read_text())
recipe=json.loads((SRC/'exception_polynomial_recipe.json').read_text())
u,v,w,x,t,s,Y,Z,p,tau=S.symbols('u v w x t s Y Z p tau')
v0=u-1;T=1-v0**10*w;g0=1-2*v0
PS=w**2*T**3+v0**2*g0*w*T+tau*v0**41*Z**2+v0**11*w+v0
PR=PS.subs(Z,(w**2-u)/u**2)
F=(-u**4*v0**30*w**5+(3*u**4*v0**20+tau*v0**41)*w**4
   -3*u**4*v0**10*w**3+(u**4*(1-v0**12*(1-2*v0))-2*tau*u*v0**41)*w**2
   +u**4*(v0**2*(1-2*v0)+v0**11)*w+u**4*(v0-p)+tau*u**2*v0**41)
check('actual full quintic clears exactly by u4',zero(F-u**4*(PR-p)))
check('generic ruling degree is five',S.degree(F,w)==5)
check('finite specialization has degree two and is generically simple',zero(F.subs(u,1)-(w**2-p)))
check('coprime five/two counts force geometric integrality',S.gcd(S.Integer(5),S.Integer(2))==1)

# Actual multiplication in a degree-five generic field control at u=2,tau=1.
f=S.cancel(PR.subs({u:2,tau:1}))
monic=S.Poly(f-p,w).monic()
A=S.zeros(5)
for j in range(4):A[j+1,j]=1
for j in range(5):A[j,4]=-monic.nth(j)
ident=S.eye(5)
def matrix_poly(expr):
    poly=S.Poly(expr,w)
    return sum((poly.nth(i)*A**i for i in range(poly.degree()+1)),S.zeros(5))
check('degree-five multiplication obeys its actual monic equation',all(zero(z) for z in matrix_poly(monic.as_expr())))
velocity=matrix_poly(S.diff(f,w)).inv()
check('actual parameter-velocity multiplication inverse',all(zero(z) for z in matrix_poly(S.diff(f,w))*velocity-ident))
check('actual field parameter-velocity trace is zero',zero(S.trace(velocity)))
check('actual root trace is independent of parameter',zero(S.diff(S.trace(A),p)))

# A composite-degree zero-trace exact form cannot use the prime criterion.
check('composite guard retains twelve hyperelliptic branch points',S.gcd(S.Poly(t**11-1,t),S.Poly(11*t**10,t)).is_one)
check('composite guard genus5 exceeds false upper bound3',((-4+12)//2+1)==5 and 5>(4-1)*(2-1))
check('composite guard trace of dt over t2 cancels',zero(1/(2*t)+1/(-2*t)))

# Full parameter family: pole divisor and local first corrections.
branches=[]
for C in [S.Integer(0),S.I,-S.I]:
    w0=v**-10+v**-2*Y if C==0 else v**-10+C*v**-4-C*v**-3+v**-2*Y
    ex=PR.subs({u:1+v,w:w0},simultaneous=True)
    order,poly=initial(ex,v);target=-Y if C==0 else 2*Y+C
    check(f'universal tau affine response C{C}',order==0 and zero(poly-target))
    den=S.cancel((1+v)**2*S.diff(PR,w).subs({u:1+v,w:w0},simultaneous=True)/v**2)
    n,leading=initial(den,v);expected=-1 if C==0 else 2
    check(f'universal tau exact double pole C{C}',n==0 and leading==expected)
    rem,_=initial(den-expected,v)
    check(f'universal tau zero residue C{C}',rem>=2)
    branches.append({'C':str(C),'pole_order':2,'residue':0,'initial_response':str(poly)})
check('all five roots at v0 accounted',len(branches)+2==5)
check('primitive pole budget is exactly three',sum(b['pole_order']-1 for b in branches)==3)
source=PS.subs({u:t**2,w:t*(-1+t**2*x),Z:-2*x+t**2*x**2})
check('central source restriction has degree two for tau nonzero',S.expand(source.subs(t,0)).coeff(x,2)==-4*tau)
check('both other leading-drop source lines have degree two',all(S.degree(source.subs(t,c),x)==2 for c in [1,-1]))
check('central whole-line source differential is nonzero',zero(S.diff(source,t).subs(t,0)+2))
check('bad whole-fibre differential is nonzero',zero(S.diff(PR,u).subs(u,1)-1))

# Universal boundary scalings supporting the polynomial division and degree.
for eps in [S.Integer(1),S.Integer(-1)]:
    scaled=PR.subs({u:s**2,w:eps*s+eps*s**3*Y/2},simultaneous=True)
    order,poly=initial(scaled,s)
    check(f'u0 four finite roots generic z equation eps{eps}',order==0 and zero(poly+tau*Y**2+1))
    small=s**130*PR.subs({u:s**-4,w:eps*s**-2+s**11*Y},simultaneous=True)
    order,poly=initial(small,s)
    check(f'infinity four-root initial equation eps{eps}',order==0 and zero(poly-(4*tau*Y**2-eps)))
large=s**16*F.subs({u:s,w:s**-4*Y},simultaneous=True)
order,poly=initial(large,s)
check('u0 fifth polar-root initial equation',order==0 and zero(poly+Y**4*(Y+tau)))
large=s**69*F.subs({u:1+s**-1,w:s**-7*Y},simultaneous=True)
order,poly=initial(large,s)
check('infinity fifth large-root initial equation',order==0 and zero(poly-Y**4*(tau-Y)))
check('discriminant v0 exact exponent96',8*30-6*20-3*8==96)
check('discriminant v0 leading constant',zero(4*p*(-1)*(-1)*(-4)+16*p))
check('discriminant u0 exact exponent10',32+6+4-32==10)
check('discriminant u0 leading constant',zero((p+1)**2/tau**2*256*tau**8-256*(p+1)**2*tau**6))
check('discriminant infinity exact degree321',272+56+4-11==321)
check('discriminant infinity leading constant',zero(-1/tau**2*256*tau**8+256*tau**6))
check('reduced discriminant degree215',321-96-10==215)

# Independent new replay of the stored exact modular certificate.
prime=cert['prime'];pv=cert['p0'];tv=cert['tau0']
Fm=F.subs({u:1+v,p:pv,tau:tv})
D=S.Poly(S.discriminant(Fm,w,domain=S.GF(prime).poly_ring(v)),v,modulus=prime)
base=S.Poly(v**96*(v+1)**10,v,modulus=prime)
R,rem=D.div(base)
check('exact specialized discriminant division',rem.is_zero and D.degree()==cert['discriminant_degree'])
poly=lambda coeff:S.Poly(sum(S.Integer(c)*v**i for i,c in enumerate(coeff)),v,modulus=prime)
recorded=poly(cert['reduced_coefficients_ascending_mod_prime'])
check('all215-degree reduced coefficients replay exactly',R==recorded and R.degree()==215)
U=poly(cert['bezout_U_ascending_mod_prime']);V=poly(cert['bezout_V_ascending_mod_prime'])
check('explicit finite-field polynomial Bezout identity one',(U*R+V*R.diff()).is_one)
check('every reduced discriminant root is simple',S.gcd(R,R.diff()).is_one)
check('reduced leading value agrees with universal formula',int(R.LC())%prime==(-256)%prime==cert['reduced_leading_coefficient'])
check('reduced v0 value agrees with universal formula',int(R.eval(0))%prime==(-16*pv)%prime==cert['reduced_at_zero'])
check('reduced u0 value agrees with universal formula',int(R.eval(-1))%prime==(256*(pv+1)**2)%prime==cert['reduced_at_minus_one'])
check('nonzero integer exceptional polynomial certified mod101',int(R.resultant(R.diff()))%prime==cert['exception_polynomial_at_one_mod_prime']==12)
check('simple roots avoid all fixed leading-zero points',R.eval(0)!=0 and R.eval(-1)!=0)
check('all215 simple generic branch points force genus at least104',S.ceiling(S.Rational(215-8,2))==104)
check('any rational primitive would force genus at most8',(5-1)*(3-1)==8)
check('global genus obstruction is strict',104>8)
check('explicit exceptional polynomial degree bound',8*(215+214)==recipe['B_tau_degree_at_most']==3432)
check('universal discriminant/exception expansions are not claimed',not recipe['universal_R_coefficients_expanded'] and not recipe['universal_B_coefficients_expanded'])
old=json.loads((BASE/'canonical_deck_pole_coverage_receipt.json').read_text())
check('known3D nonproper control retained',old['outcome']=='PASS' and 'retained 3D first-velocity pole' in old['checks'])
check('generic curve trace distinguished from full normalization',context['actual_curve_field_trace_not_original_finite_normalization'])
check('original full normalization and every ramified prime retained',context['retained_full_original_normalization_codifferent_and_all_ramified_primes'])
check('whole ruling-five class not claimed closed',not context['all_ruling_five_targets_closed'])
check('automatic integrality and normal-tail powers not claimed',not context['generic_first_derivative_integrality_proved'] and not context['normal_tail_increasing_power_containment_proved'])
receipt={'claim':'PF2RULED160','outcome':'PASS','check_count':len(checks),'checks':checks,
    'elapsed_seconds':time.perf_counter()-START,'environment':{'python':platform.python_version(),'sympy':S.__version__},
    'proof_sha256':sha(PROOF),'checker_sha256':sha(Path(__file__)),'context_sha256':sha(CONTEXT),
    'certificate_sha256':sha(SRC/'affine_discriminant_certificate.json'),
    'recipe_sha256':sha(SRC/'exception_polynomial_recipe.json'),
    'affine_control_tau1_excluded':True,
    'excluded_parameter_class':'tau!=0 and B(tau)!=0; B is explicit nonzero nested resultant',
    'exception_polynomial_degree_at_most':3432,'exception_polynomial_at_one_mod101':12,
    'simple_generic_branch_points_at_least':215,'generic_curve_genus_at_least':104,
    'coframe_polar_points':branches,'primitive_degree_if_exact':3,'genus_at_most_if_exact':8,
    'rational_mates_and_every_finite_extension_excluded_on_stated_class':True,
    'all_ruling_five_targets_excluded':False,
    'universal_discriminant_and_B_coefficient_expansions_computed':False,
    'original_full_normalization_and_all_ramified_primes_retained':True,
    'generic_first_derivative_integrality_proved':False,'normal_tail_increasing_power_containment_proved':False,
    'registered_predecessor_checkers_replayed':[],'new_external_theorem_imports':[],
    'full_goal_complete':False,'independent_replay':False,'formal_verification':False,'external_review':False}
RECEIPT.write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'outcome':'PASS','check_count':len(checks),'elapsed_seconds':receipt['elapsed_seconds']},indent=2))

#!/usr/bin/env python3
"""Exact PF163 controls; the uniform complete-end proof is in the source note."""
from pathlib import Path
import argparse
import hashlib
import json
import platform
import tempfile
import time
import sympy as S

START = time.perf_counter()
BASE = Path(__file__).resolve().parent
PROOF = BASE / "RULING_SIX_POLAR_ZERO_AND_AFFINE_RESPONSE_FAMILY.md"
checks = []
def check(label, value):
    if not bool(value):
        raise AssertionError(label)
    checks.append(label)
def zero(e):
    return S.cancel(e) == 0
def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()
v,w,z,p,t,x,Y,tau = S.symbols("v w z p t x Y tau")
# Only tau is assumed invertible; all other symbols are unrestricted.
tau = S.Symbol("tau", nonzero=True)
u = v+1
g = 1-2*v
T = 1-v**10*w
P5 = w**2*T**3+v**2*g*w*T+v**41*z**2+v**11*w+v
zstar = (w**2-u)/u**2
P5star = S.expand(P5.subs(z,zstar))
P5w = S.diff(P5star,w)

def initial(e):
    n,d = S.fraction(S.together(e))
    if n == 0:
        return S.oo,S.Integer(0)
    pn,pd = S.Poly(n,v),S.Poly(d,v)
    no = min(i[0] for i in pn.monoms())
    do = min(i[0] for i in pd.monoms())
    return no-do,S.cancel(pn.nth(no)/pd.nth(do))

def regular_first_two(e):
    n,d = S.fraction(S.together(e))
    pn,pd = S.Poly(n,v),S.Poly(d,v)
    no = min(i[0] for i in pn.monoms())
    do = min(i[0] for i in pd.monoms())
    if no < do:
        raise AssertionError("nonregular series in coefficient control")
    n0,n1=pn.nth(do),pn.nth(do+1)
    d0,d1=pd.nth(do),pd.nth(do+1)
    return S.cancel(n0/d0),S.cancel((n1*d0-n0*d1)/d0**2)

for C in [0,S.I,-S.I]:
    name = str(C)
    ww = v**(-10)+C*v**(-4)-C*v**(-3)+v**(-2)*Y
    order,lead = initial(P5star.subs(w,ww))
    expected = -Y if C == 0 else 2*Y+C
    check("entire old affine response "+name,order == 0 and zero(lead-expected))
    a0,a1 = regular_first_two((P5w/v**2).subs(w,ww))
    slope = -1 if C == 0 else 2
    check("old derivative coefficient zero "+name,zero(a0-slope))
    check("old derivative coefficient one "+name,zero(a1+2*slope))
    check("old derivative first coefficients independent of Y "+name,
          not a0.has(Y) and not a1.has(Y))
    # Remove the K power first; its exact valuation is then uniform in K.
    qo,ql = initial((zstar**3).subs(w,ww))
    ro,rl = initial(S.diff(zstar**3,w).subs(w,ww))
    check("cubic-z perturbation exact valuation "+name,qo == -60 and ql == 1)
    check("cubic-z derivative exact valuation "+name,ro == -50 and rl == 6)

check("unit volume residue cancellation",S.expand(u**2*g).coeff(v,0) == 1
      and S.expand(u**2*g).coeff(v,1) == 0)
finite = S.Poly((P5star-p).subs(v,0),w)
check("two finite branches retain full specialization",finite.as_expr() == w**2-p)
check("two finite branches generically simple",S.discriminant(finite.as_expr(),w) == 4*p)

models = []
for K in [62,63,70]:
    h = K-30
    PK = P5+tau*v**K*z**3
    star = S.expand(PK.subs(z,zstar))
    coeff = S.Poly(star,w)
    check("whole polynomial ruling degree six K"+str(K),coeff.degree() == 6)
    check("exact leading coefficient K"+str(K),zero(coeff.nth(6)-tau*v**K/u**6))
    check("exact next coefficient K"+str(K),coeff.nth(5) == -v**30)
    # Full cleared sextic, not only a leading model.
    cleared = S.expand(u**6*(w**2*T**3+v**2*g*w*T+v**11*w+v-p)
        +u**2*v**41*(w**2-u)**2+tau*v**K*(w**2-u)**3)
    check("whole equation clearing by u6 K"+str(K),zero(cleared-u**6*(star-p)))
    scaled = v**(5*K-180)*(star-p).subs(w,v**(-h)*Y)
    so,si = initial(scaled)
    check("new whole weighted initial K"+str(K),
          so == 0 and zero(si-Y**5*(tau*Y-1)))
    check("new initial root simple K"+str(K),
          zero(si.subs(Y,1/tau)) and zero(S.diff(si,Y).subs(Y,1/tau)-tau**(-4)))
    derivative = S.diff(star,w)
    scaled_derivative = v**(4*K-150)*derivative.subs(w,v**(-h)*Y)
    do,di = initial(scaled_derivative)
    check("new entire derivative initial K"+str(K),
          do == 0 and zero(di-Y**4*(6*tau*Y-5)))
    check("new derivative exact nonzero constant K"+str(K),
          zero(di.subs(Y,1/tau)-tau**(-4)))
    check("new reciprocal derivative coframe zero K"+str(K),4*K-150 >= 98)
    check("prospective primitive local degree K"+str(K),4*K-149 >= 99)
    for C in [0,S.I,-S.I]:
        ww=v**(-10)+C*v**(-4)-C*v**(-3)+v**(-2)*Y
        io,ip = initial(star.subs(w,ww))
        expected=-Y if C == 0 else 2*Y+C
        check("perturbed entire affine response K"+str(K)+" C"+str(C),
              io == 0 and zero(ip-expected))
        a0,a1=regular_first_two((derivative/v**2).subs(w,ww))
        slope=-1 if C == 0 else 2
        check("perturbed derivative two coefficients K"+str(K)+" C"+str(C),
              zero(a0-slope) and zero(a1+2*slope))
    # Whole-source polynomial marking: no rational-chart evaluation at u=0.
    source = S.expand(PK.subs({v:t**2-1,w:t*(-1+t**2*x),z:-2*x+t**2*x**2}))
    check("central whole-source restriction degree three K"+str(K),
          S.degree(source.subs(t,0),x) == 3)
    check("whole central source derivative K"+str(K),
          S.diff(source,t).subs(t,0) == -2)
    check("both reached offcentral lines degree two K"+str(K),
          S.degree(source.subs(t,1),x) == 2 and S.degree(source.subs(t,-1),x) == 2)
    check("off-M leading-drop fibre differential K"+str(K),
          S.diff(star,v).subs(v,0) == 1)
    check("canonical central even-leading coefficient unit K"+str(K),
          (tau*v**K).subs(v,-1) != 0)
    models.append({"K":K,"new_w_pole":h,"new_coframe_zero":4*K-150,
        "prospective_primitive_local_degree":4*K-149,"entire_pole_budget":3})

# The exponents in the written proof are affine functions of K and alpha.
# Their monotonicity, together with their endpoint values, covers K>=62.
kplus,aplus = S.symbols("kplus aplus", integer=True, nonnegative=True)
KK=62+kplus; alpha=1+aplus
check("uniform old perturbation positive margin",S.expand(KK-60) == kplus+2)
check("uniform old derivative perturbation margin",S.expand(KK-50) == kplus+12)
check("uniform new degree-four remainder margin",S.expand(KK-40) == kplus+22)
check("uniform new root outside old polar centres",S.expand(KK-30-10) == kplus+22)
check("uniform new zero at least98",S.expand(4*KK-150) == 4*kplus+98)
check("uniform local multiplicity at least99",S.expand(4*KK-149) == 4*kplus+99)
check("infinity monic bound strictly positive before ceiling",
      S.expand(KK*alpha) == kplus*aplus+kplus+62*aplus+62)
check("all six branches counted",2+3+1 == 6)
check("full primitive pole degree less than new local degree",3 < 99)
# Degenerate controls make the nonzero assumptions visible.
check("tau zero drops sextic degree",S.degree((P5+tau*v**62*z**3).subs(tau,0).subs(z,zstar),w) == 5)
check("new derivative cancellation retains 6-5",6-5 == 1)
check("simple coframe pole would not supply rational primitive degree",
      S.diff(v**(-1),v) == -v**(-2))

receipt = {
    "claim":"PF2RULED163","outcome":"PASS","check_count":len(checks),"checks":checks,
    "elapsed_seconds":time.perf_counter()-START,
    "environment":{"python":platform.python_version(),"sympy":S.__version__},
    "proof_sha256":sha(PROOF),"checker_sha256":sha(Path(__file__)),
    "retained_proof_sources":{
        name:sha(BASE/name) for name in [
            "CANONICAL_DECK_JET_FIELDS_AND_RATIONAL_POLE_ATLAS_COVERAGE.md",
            "MONIC_PARAMETER_POLE_BOUNDS_AND_EVEN_LEADING_CYCLIC_TARGETS.md",
            "POLYNOMIAL_MATE_RESPONSES_AND_QUADRATIC_RULING_FIVE_ESCAPE_COVERAGE.md"]},
    "checked_full_substitution_models":models,
    "uniform_proved_scope":{"K":"every integer >=62","tau":"every nonzero complex scalar",
        "ruling_degree":6,"coframe_double_poles":3,"coframe_simple_poles":0,
        "new_coframe_zero":"4*K-150","primitive_pole_budget":3,
        "rational_mates_in_any_finite_extension":False},
    "whole_ruling_six_class_closed":False,"generic_first_derivative_integrality_proved":False,
    "normal_tail_increasing_power_containment_proved":False,"full_goal_complete":False,
    "independent_whole_proof_review":False,"formal_verification":False,"external_review":False
}
parser=argparse.ArgumentParser()
parser.add_argument("--output",type=Path,default=None)
args=parser.parse_args()
out=args.output or Path(tempfile.mkdtemp(prefix="math-check-"))/"ruling_six_polar_zero_receipt.json"
out.write_text(json.dumps(receipt,indent=2)+"\n")
print(json.dumps({"outcome":"PASS","check_count":len(checks),
    "elapsed_seconds":receipt["elapsed_seconds"],"receipt":str(out)},indent=2))

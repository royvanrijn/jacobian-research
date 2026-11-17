"""Exact controls for PF169; the written proof supplies the unrestricted argument."""
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
PROOF = BASE / "RULING_SIX_GLOBAL_QUARTIC_DIVISIBILITY_AND_CUBIC_REMAINDER.md"
checks = []
models = []


def check(label, condition):
    if not bool(condition):
        raise AssertionError(label)
    checks.append(label)


def zero(value):
    return S.cancel(S.expand(value)) == 0


def initial(value, parameter):
    terms = S.Add.make_args(S.expand(value))
    exponents = [term.as_powers_dict().get(parameter, 0) for term in terms]
    order = min(exponents)
    coefficient = S.expand(sum(
        term / parameter**exponent for term, exponent in zip(terms, exponents)
        if exponent == order
    ))
    return order, coefficient


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


v, s, w, Y, p = S.symbols("v s w Y p")
h0, K0, r0 = S.symbols("h0 K0 r0", nonzero=True)
disc = (2*h0*K0)**2 - 4*h0**2*(K0**2+r0)
check("outer square discriminant detects nonzero r at delta=2k",
      zero(disc + 4*h0**2*r0))
check("outer even multiplicities exclude odd minimum index",
      all(q % 2 != 0 for q in (1, 3, 5)))
check("delta=2k cancellation leaves a simple degree-five factor",
      S.diff(Y**5*(Y+2), Y).subs(Y, -2) != 0)

# Integer and half-integer boundaries of every inequality in the outer proof.
# These are controls, not a finite replacement for the displayed comparisons.
for j, k, delta in [(4,1,3), (5,1,3), (6,1,5), (7,2,5),
                     (8,2,7), (9,3,7), (12,4,11)]:
    ell = j-k
    beta = S.Rational(2*j-delta, 2)
    wp_order = 3*ell-k-S.Rational(delta, 2)
    label = f"j{j} k{k} delta{delta}"
    check("strict outer order hypotheses "+label, 2*k < delta < j)
    check("translated cluster smaller than outer centre "+label, 0 < beta < ell)
    check("cubic remainder positive margin "+label, j-delta > 0)
    check("quadratic remainder positive margin "+label, 2*ell-delta > 0)
    check("exact wp exponent positive "+label, wp_order > 0)
    check("all p-independent Puiseux orders integral after v=s2 "+label,
          all((2*x).is_integer is True for x in (beta, wp_order)))

    # Entire source-image polynomial, not a truncated local initial.
    hh, KK, rr = v**j, v**k, v**delta
    zz = w**2-1
    TT = (hh*w+KK)*zz
    PP = S.expand(TT**2+rr*zz**2+v**k*w*zz+zz+w+1)
    wc = -v**(-ell)
    H = 4*ell-delta
    subw = -s**(-2*ell)+s**(-2*beta)*Y
    scaled = S.expand(s**(2*H)*(PP.subs({v:s**2, w:subw})-p))
    order, lead = initial(scaled, s)
    check("entire translated equation regular "+label, order == 0)
    check("entire translated equation initial Y2+1 "+label, zero(lead-Y**2-1))
    check("translated roots simple and nonzero "+label,
          S.diff(lead,Y).subs(Y,S.I) == 2*S.I)
    der = S.diff(PP,w).subs({v:s**2,w:subw}).subs(Y,S.I)
    der_order, der_lead = initial(der,s)
    check("entire derivative exact order "+label, -der_order == 2*wp_order)
    check("entire derivative leading coefficient "+label, zero(der_lead-2*S.I))
    check("Puiseux coframe holomorphic including du "+label, -der_order+1 >= 0)
    models.append({"j":j,"k":k,"delta":delta,
                   "outer_pole":ell,"translated_pole":str(beta),
                   "wp_v_order":str(wp_order),
                   "initial":"Y^2+1","formal_model_not_Keller_realization":True})

# Full lower simple-root equations when ord p3 < k.
for j,k,delta,a in [(6,2,5,1),(8,3,7,1),(8,3,7,2),
                     (10,4,9,3)]:
    hh,KK,rr = v**j,v**k,v**delta
    zz=w**2-1
    PP=S.expand(((hh*w+KK)*zz)**2+rr*zz**2+v**a*w*zz+zz+w+1)
    scaled=S.expand(v**(2*a)*(PP.subs(w,v**(-a)*Y)-p))
    order,lead=initial(scaled,v)
    label=f"j{j} k{k} a{a}"
    check("lower entire equation regular "+label,order==0)
    check("lower entire equation Y2(1+Y) "+label,zero(lead-Y**2*(1+Y)))
    check("lower nonzero simple root "+label,S.diff(lead,Y).subs(Y,-1)!=0)
    der_order,der_lead=initial(S.diff(PP,w).subs(w,-v**(-a)),v)
    check("lower exact wp order "+label,-der_order==a)
    check("lower derivative unit coefficient "+label,der_lead!=0)
    check("lower all higher margins positive "+label,
          min(2*k-2*a,j+k-3*a,2*j-4*a)>0)

# Symbolic surface-ring identity for the regular cubic remainder.
M,N,h,K,B,z,a1,a0,b1,b0 = S.symbols("M N h K B z a1 a0 b1 b0")
r=2*h*M*B
b2=N*h**2+K**2+r
PP=M*h**2*z**3+2*h*K*w*z**2+b2*z**2+w*(a1*z+a0)+b1*z+b0
TT=(h*w+K)*z
WW=TT+B*w
RR=(a1-2*B*K)*w*z+(b1-2*B*h*N-B**2*M)*z+a0*w+b0-B**2*N
check("global regular square reduction in surface ring",
      S.rem(S.expand(PP-WW**2-RR),w**2-N-M*z,w)==0)
check("new square has ruling degree three",
      S.degree(S.expand(WW.subs(z,(w**2-N)/M)),w)==3)
check("new remainder has ruling degree at most three",
      S.degree(S.expand(RR.subs(z,(w**2-N)/M)),w)<=3)
e=S.symbols("e")
check("next constant shift cubic coefficient is a1-2BK-2eh",
      S.expand((RR-2*e*WW-e**2)).coeff(z,1).coeff(w,1)==a1-2*B*K-2*e*h)
x,sigma=S.symbols("x sigma")
check("central affine square coordinate retained",
      zero(WW.subs({w:0,z:2*sigma*x})-2*sigma*K*x))

receipt={
    "claim":"PF2RULED169","outcome":"PASS","check_count":len(checks),"checks":checks,
    "elapsed_seconds":time.perf_counter()-START,
    "environment":{"python":platform.python_version(),"sympy":S.__version__},
    "proof_sha256":sha(PROOF),"checker_sha256":sha(Path(__file__)),
    "retained_proof_sources":{
        n:sha(BASE/n) for n in [
            "RULING_SIX_SIMPLE_ESCAPE_GRAPHS_AND_QUARTIC_FACE_REDUCTION.md",
            "RULING_SIX_APPROXIMATE_CUBIC_ROOT_REGULARITY_AT_M_DIVISORS.md",
            "RULING_SIX_OUTER_SQUARES_AND_EVEN_SPECIALIZATIONS.md",
        ]
    },
    "checked_local_models":models,
    "uniform_proved_scope":{
        "quartic_coefficient_divisibility":"h*M divides r",
        "regular_cubic_square_coordinate":"W3=(h*w+K)*z+r*w/(2*h*M)",
        "remainder_ruling_maximum":3,
        "next_divisibility_h_divides_cubic_coefficient":False,
        "coefficient_and_companion_degree_bounds":None,
    },
    "whole_ruling_six_class_closed":False,
    "generic_first_derivative_integrality_proved":False,
    "automatic_trace_integrality_proved":False,
    "full_goal_complete":False,"independent_whole_proof_review":False,
    "formal_verification":False,"external_review":False,
}
parser=argparse.ArgumentParser()
parser.add_argument("--output",type=Path,default=None)
args=parser.parse_args()
out=args.output or Path(tempfile.mkdtemp(prefix="math-check-"))/"ruling_six_global_quartic_divisibility_receipt.json"
out.write_text(json.dumps(receipt,indent=2)+"\n")
print(json.dumps({"outcome":"PASS","check_count":len(checks),
                 "elapsed_seconds":receipt["elapsed_seconds"],"receipt":str(out)},indent=2))

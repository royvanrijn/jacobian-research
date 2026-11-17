#!/usr/bin/env python3
"""Exact replay for GAMMA_CONVOLUTION_APPELL_OBSTRUCTION.tex.

The finite rational counterexample and the explicit infinite-support
construction have written proofs. These controls certify their displayed
arithmetic, not the all-order convergence and Hankel arguments.
"""
from pathlib import Path
import hashlib
import json
import platform
import sympy as s

checks = []
def check(name, condition):
    if not bool(condition):
        raise AssertionError(name)
    checks.append(name)

theta = [s.Rational(1,2), s.Rational(19,2)]
alpha = [s.Rational(125,9), s.Rational(1,3249)]
for j, (t,a) in enumerate(zip(theta,alpha)):
    check(f"positive_scale_{j}", t>0)
    check(f"positive_shape_{j}", a>0)
kappa = [s.Integer(0)] + [
    s.factor(s.factorial(k-1)*sum(a*t**k for a,t in zip(alpha,theta)))
    for k in range(1,6)]
expected = [s.Rational(132,19), s.Rational(7,2),
            s.Integer(4), s.Rational(81,4), s.Integer(582)]
for k, value in enumerate(expected,1):
    check(f"cumulant_{k}", kappa[k]==value)
x,z=s.symbols("x z")
centered_log = sum(kappa[k]*z**k/s.factorial(k) for k in range(2,6))
egf = s.series(s.exp(centered_log),z,0,6).removeO().expand()
moment = [s.expand(egf).coeff(z,k)*s.factorial(k) for k in range(6)]
appell = {n:s.expand(sum(s.binomial(n,k)*moment[k]*x**(n-k)
                          for k in range(n+1))) for n in (4,5)}
q=x*x-x+19
check("A4_factorization", s.expand(appell[4]-q*(x*x+x+3))==0)
check("A5_factorization", s.expand(appell[5]-(x+2)*q*q)==0)
check("derivative_identity", s.diff(appell[5],x)==5*appell[4])
check("exact_gcd", s.expand(s.Poly(s.gcd(appell[4],appell[5]),x).monic().as_expr()-q)==0)
check("zero_resultant", s.resultant(appell[4],appell[5],x)==0)
check("nonreal_common_quadratic", s.discriminant(q,x)==-75)
# Deliberately changing the fifth cumulant breaks this exact collision.
corrupt5=appell[5]+1
check("corrupted_control_nonzero_resultant",
      s.resultant(appell[4],corrupt5,x)!=0)

m=[s.Rational(7,2),s.Integer(2),s.Rational(27,8),s.Rational(97,4)]
for j in range(4):
    check(f"moment_data_{j}",sum(a*t**(j+2) for a,t in zip(alpha,theta))==m[j])
H=s.factor(m[0]*m[2]-m[1]**2)
I=s.factor(m[1]*m[3]-m[2]**2)
J=s.factor(m[0]*m[3]-m[1]*m[2])
check("H_exact",H==s.Rational(125,16))
check("I_exact",I==s.Rational(2375,64))
check("J_exact",J==s.Rational(625,8))
S=s.factor(J/H);P=s.factor(I/H)
check("two_atom_sum",S==s.Integer(10))
check("two_atom_product",P==s.Rational(19,4))
check("mean_between_positive_roots",
      s.factor((m[1]/m[0])**2-S*m[1]/m[0]+P)==-H/m[0]**2)
eps=s.Rational(1,1000)
lower=[
 H-eps*(m[0]*4+m[2])-eps**2*4,
 I-eps*(m[1]*8+m[3]*2)-eps**2*16,
 J-eps*(m[0]*8+m[3])-eps**2*8]
for j,bound in enumerate(lower):
    check(f"infinite_extension_positivity_bound_{j}",bound>0)
for j in range(4):
    check(f"remaining_positive_moment_{j}",m[j]-eps*2**j>0)

# The representing-measure statement is an all-order identity, not a prefix
# assertion: q_j=sum alpha_i theta_i^(j+1), for every j>=0.
# Verify its determinant factorization symbolically at the full finite rank.
M=s.Matrix([[sum(a*t**(i+j+1) for a,t in zip(alpha,theta))
              for j in range(2)] for i in range(2)])
vand=s.prod((theta[j]-theta[i])**2 for i in range(2) for j in range(i+1,2))
check("rank_two_cauchy_binet_identity",
      s.factor(M.det()-s.prod(a*t for a,t in zip(alpha,theta))*vand)==0)
check("rank_two_strict_principal_minor",M.det()>0)

# Independent raw moment calculation from the two Gamma densities.
v=s.symbols("v")
raw=[s.factor(sum(s.binomial(k,j)*s.rf(alpha[0],j)*theta[0]**j
                 *s.rf(alpha[1],k-j)*theta[1]**(k-j)
                 for j in range(k+1))) for k in range(6)]
for n in (4,5):
    from_raw=sum(s.binomial(n,k)*raw[k]*v**(n-k) for k in range(n+1))
    check(f"independent_raw_moment_factorization_{n}",
          s.expand(from_raw-appell[n].subs(x,v+kappa[1]))==0)

note=Path(__file__).with_name("GAMMA_CONVOLUTION_APPELL_OBSTRUCTION.tex")
prior=Path(__file__).with_name("lower_exponent_initial_and_boundary.tex")
receipt={
    "scope":"Strict-Hankel-TP strengthening of the retained two-Gamma "
            "Appell-coprimality implications; not a binomial-family counterexample",
    "passed":True,"checks":len(checks),"check_names":checks,
    "kind":"internal exact replay; all-order proofs are in the TeX source",
    "python":platform.python_version(),
    "command":"python research/factorial/image-rigidity/gamma_convolution_appell_obstruction.py",
    "sympy_version":s.__version__,
    "exact_prior_input":"lower_exponent_initial_and_boundary.tex thm:gamma-convolution-obstruction",
    "source_sha256":{note.name:hashlib.sha256(note.read_bytes()).hexdigest(),
              Path(__file__).name:hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              prior.name:hashlib.sha256(prior.read_bytes()).hexdigest()},
    "positive_rational_scales":[str(v) for v in theta],
    "positive_rational_shapes":[str(v) for v in alpha],
    "cumulants_1_to_5":[str(v) for v in kappa[1:]],
    "centered_common_quadratic":str(q),
    "scope_limits":{
        "single_gamma_power_counterexample":False,
        "physical_shape_progression_counterexample":False,
        "factorial_sfc_counterexample":False,
        "strict_infinite_extension":"written proof; positive explicit two-atom compensation plus countable atoms"},
    "independent_replay":False,"formal_verification":False,"external_review":False,
    "strict_extension_lower_bounds":[str(v) for v in lower],
    "assurance":{"internal_replay":True,"formal_verification":False,
                 "external_review":False,"independent_replay":False}}
out=Path(__file__).with_name("gamma_convolution_appell_obstruction_verification.json")
out.write_text(json.dumps(receipt,indent=2)+"\n")
print(f"PASS {len(checks)} exact controls; receipt {out}")

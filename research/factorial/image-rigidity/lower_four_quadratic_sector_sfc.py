#!/usr/bin/env python3
"""Exact symbolic controls for the written all-real-shape sector proof."""
from __future__ import annotations
import hashlib
import json
from pathlib import Path
import platform
import sys
import sympy as sp

HERE = Path(__file__).resolve().parent
PROOF = HERE/"LOWER_FOUR_QUADRATIC_SECTOR_SFC.tex"
IMPORTS = {
    "LOWER_FOUR_QUADRATIC_GAP_SFC.tex":
        "ea2c5885f6f03a06ca2aa12840eddb191c128c20848ca8cf83d472268b89aca4",
    "LOWER_FIVE_QUADRATIC_SECTOR_SFC.tex":
        "150d433c66c8bd1ac0c82bbe0fff6a53ad701c2aede630e961c783f27c8facce",
}
for name,digest in IMPORTS.items():
    actual=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
    if actual != digest:
        raise AssertionError(f"changed exact proof input {name}: {actual}")
checks = []
def check(name, expression):
    value = sp.cancel(sp.expand(expression))
    if value != 0:
        raise AssertionError(f"{name}: nonzero exact remainder {sp.factor(value)}")
    checks.append(name)

def positive_polynomial(name, expression, variable):
    poly = sp.Poly(sp.expand(expression), variable)
    if not poly or any(c < 0 for c in poly.all_coeffs()) or not any(c > 0 for c in poly.all_coeffs()):
        raise AssertionError(f"{name}: not a nonzero polynomial with nonnegative coefficients")
    checks.append(name)

n, M, u = sp.symbols("n M u")
gamma = M - 2*n
s = 2*(n+1)
A = sp.rf(M+2, 4)
k = (M+4)*(M+5)
b = s*(2*M+7)
U = A + b*u
W = s*u*(k+u)

# Independent rising-jet Ore reduction for the general shape-shift-four input.
def E(vector):
    r, e, j = vector
    return [sp.expand(-2*u*sp.diff(r,u)-gamma*u*j),
            sp.expand(r-e-2*u*sp.diff(e,u)+u*j),
            sp.expand(e+M*j-2*u*sp.diff(j,u))]
jet = [sp.Integer(1), sp.Integer(0), sp.Integer(0)]
jets = {0:jet}
for j in range(6):
    acted = E(jet)
    jet = [sp.expand(a+j*v) for a,v in zip(acted,jet)]
    jets[j+1] = jet
upper = [sp.expand(u*a-c) for a,c in zip(jets[4],jets[6])]
check("full symbolic offset-four K1 coefficient", upper[1]+W)
check("full symbolic offset-four K2 coefficient", upper[2]+U)

# The imaginary part is checked with independent real scalar placeholders.
x,y,AA,bb,kk,ss = sp.symbols("x y AA bb kk ss", real=True)
uu = x+sp.I*y
uu_bar = x-sp.I*y
UU = AA+bb*uu
WW = ss*uu*(kk+uu)
UU_bar = AA+bb*uu_bar
WW_bar = ss*uu_bar*(kk+uu_bar)
imag_over_y = sp.cancel((WW*UU_bar-WW_bar*UU)/(2*sp.I*y))
check("full imaginary numerator identity",
      imag_over_y-ss*(AA*kk+2*AA*x+bb*(x*x+y*y)))

Q = 2*(n+1)*(2*M+7)/((M+2)*(M+3))
check("necessary angle ratio normalization", b*k/A-Q)
check("negative derivative factor",
      sp.diff(Q,M)+2*(n+1)*(2*M*M+14*M+23)/((M+2)**2*(M+3)**2))
positive_polynomial("Q derivative numerator positive", 2*M*M+14*M+23, M)

# Exact all-n endpoint certificates.  Positivity follows coefficientwise.
numerator = 2*(n+1)*(2*M+7)
denominator = (M+2)*(M+3)
endpoints = [
    ("beta at n/2", sp.Rational(13,2)*n, sp.Rational(8,13), 2*(130*n+67)/13),
    ("beta at 2n", 8*n, sp.Rational(1,2), 26*n+11),
    ("beta at 8n", 14*n, sp.Rational(2,7), 2*(175*n+43)/7),
]
for name, endpoint, target, margin in endpoints:
    actual = (numerator-target*denominator).subs(M,endpoint)
    check(name+" exact cleared margin", actual-margin)
    positive_polynomial(name+" margin positive for all n>0", margin,n)

a = sp.symbols("a")
R = 4*(n-1)*(a+1)/(a-n+1)**2
check("Gamma-J ratio decreasing factor",
      sp.diff(R,a)+4*(n-1)*(a+n+1)/(a-n+1)**3)
ratio_cases = [
    ("first region",2*n,sp.Integer(8),4*(5*n+3),(n+1)**2),
    ("second region",sp.Rational(9,4)*n,sp.Rational(144,25),
     16*(485*n+244),25*(5*n+4)**2),
    ("third region",3*n,sp.Integer(3),20*n+7,(2*n+1)**2),
    ("fourth region",6*n,sp.Integer(1),n*n+30*n+5,(5*n+1)**2),
]
for name, alpha, bound, num, den in ratio_cases:
    check(name+" exact Gamma-J margin", bound-R.subs(a,alpha)-num/den)
    positive_polynomial(name+" Gamma-J numerator positive",num,n)

angular_cases = [
    ("first angle",sp.Integer(8),sp.Rational(8,13),sp.Rational(7,9),sp.Rational(11,1053)),
    ("second angle",sp.Rational(144,25),sp.Rational(1,2),sp.Rational(119,169),sp.Rational(239,57122)),
    ("third angle",sp.Integer(3),sp.Rational(2,7),sp.Rational(1,2),sp.Rational(1,28)),
]
for name, q, common_bound, cosine, gap in angular_cases:
    check(name+" double-angle cosine", (q-1)/(q+1)-cosine)
    check(name+" exact strict angular gap",common_bound-cosine*cosine-gap)
    if gap <= 0:
        raise AssertionError(name+" gap nonpositive")
    checks.append(name+" positive angular gap")
check("fourth angle right half-plane", (sp.Integer(1)-1)/(sp.Integer(1)+1))

# Reciprocal normalization in the Gamma-product sector lemma, exactly at
# symbolic shapes and small degrees as regression controls, not a prefix theorem.
t,ga = sp.symbols("t ga")
for degree in range(1,5):
    aa=ga/2
    C=sum(sp.binomial(degree,j)*sp.rf(aa,j)*sp.rf(aa+sp.Rational(1,2),j)
          *u**(degree-j) for j in range(degree+1))
    T=sum(sp.binomial(degree,j)*sp.rf(ga,2*j)*t**j for j in range(degree+1))
    check(f"Gamma duplication/reversal degree {degree}",
          sp.expand((4*t)**degree*C.subs(u,1/(4*t))-T))

# Generic nonreal partial fraction: each positive Pick coefficient remains so
# in both half-planes, not merely in the right half-plane.
aa = sp.symbols("aa",real=True)
pick = uu/(aa*(uu+aa))
pick_bar = uu_bar/(aa*(uu_bar+aa))
check("Pick elementary imaginary identity",
      (pick-pick_bar)/(2*sp.I*y)-1/((x+aa)**2+y*y))

# Corruption controls ensure that the detector does not accept the plausible
# missing k factor or the wrong shift coefficient.
bad_ratio = b/A
if sp.cancel(bad_ratio-Q) == 0:
    raise AssertionError("corrupted angle normalization accepted")
checks.append("corrupted missing-k normalization rejected")
if sp.cancel(upper[2]+U+s*u) == 0:
    raise AssertionError("corrupted offset-four coefficient accepted")
checks.append("corrupted upper coefficient rejected")

inputs = {
    PROOF.name: PROOF,
    **{name:HERE/name for name in IMPORTS},
}
for key,path in inputs.items():
    if not path.is_file():
        raise AssertionError(f"missing proof input {key}")
# Exact named imported statements are present; their written proofs are not
# certified by this existence check.
four = inputs["LOWER_FOUR_QUADRATIC_GAP_SFC.tex"].read_text()
five = inputs["LOWER_FIVE_QUADRATIC_SECTOR_SFC.tex"].read_text()
for label,text in [
    ("lem:offset-four-transform",four),
    ("lem:pick-general",five),
    ("lem:moment-sector",five),
]:
    if "\\label{"+label+"}" not in text:
        raise AssertionError(f"missing imported statement {label}")

sha = lambda path:hashlib.sha256(path.read_bytes()).hexdigest()
receipt = {
    "passed":True,
    "status":"PASS",
    "symbolic_controls":len(checks),
    "total_checks":len(checks),
    "control_names":checks,
    "proof_scope":"all integers n>=1 and all real beta>0 for support {4,6}; written sector argument supplies the quantifiers",
    "finite_normalization_controls_are_not_an_index_theorem":True,
    "scope":{
        "all_positive_real_Gamma_beta":True,
        "lower_four_quadratic_gap_all_indices":True,
        "ordinary_support_four_six_all_indices":True,
        "lower_moment_uniform_Hurwitz_assumption":False,
        "algebraic_or_p_integral_shape_restriction":False,
        "all_degree_seven_binomial_SFC":False,
        "unrestricted_FC_or_SFC":False,
    },
    "source_sha256":{**{key:sha(path) for key,path in inputs.items()},
                     Path(__file__).name:sha(Path(__file__))},
    "exact_imports":{
        "general_shape_shift_four":"LOWER_FOUR_QUADRATIC_GAP_SFC.tex lem:offset-four-transform",
        "strict_Pick":"LOWER_FIVE_QUADRATIC_SECTOR_SFC.tex lem:pick-general",
        "Gamma_product_sector":"LOWER_FIVE_QUADRATIC_SECTOR_SFC.tex lem:moment-sector",
    },
    "checker_sha256":sha(Path(__file__)),
    "environment":{"python":platform.python_version(),"sympy":sp.__version__},
    "assurance":{"independent_replay":False,"formal_verification":False,"external_review":False},
    "independent_replay":False,"formal_verification":False,"external_review":False,
    "novelty_claim":False,"priority_claim":False,
    "command":".venv/bin/python research/factorial/image-rigidity/lower_four_quadratic_sector_sfc.py",
}
(HERE/"lower_four_quadratic_sector_sfc_verification.json").write_text(json.dumps(receipt,indent=2)+"\n")
print(f"PASS {len(checks)} exact symbolic controls")
for key,path in inputs.items():
    print(sha(path),key)
print(sha(Path(__file__)),Path(__file__).name)

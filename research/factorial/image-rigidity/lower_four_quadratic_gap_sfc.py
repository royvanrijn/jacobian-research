"""Exact internal replay for the ordinary support {4,6} scalar proof.

The all-index conclusion is supplied by the written proof. These controls
verify its generic differential algebra and complete finite specializations.
Run with the repository virtual environment, which supplies SymPy.
"""
from pathlib import Path
import hashlib
import json
import platform
import sympy as sp

HERE = Path(__file__).resolve().parent
PROOF = HERE / "LOWER_FOUR_QUADRATIC_GAP_SFC.tex"
IMPORTS = {
    "CUBIC_ADJOINT_STABILITY_SFC.tex":
        "ea2510e511545b9dd7a49a0796d75c45575687e60770c3359d38c5bc99cb4125",
    "lower_exponent_initial_and_boundary.tex":
        "e66e02d0354de5cd92d33eef54efea15c387b9f6a0c1b62ed379295e13f07d45",
}
for name, wanted in IMPORTS.items():
    assert hashlib.sha256((HERE/name).read_bytes()).hexdigest() == wanted, (name, "exact import changed")
counts = {}
def count(key, amount=1):
    counts[key] = counts.get(key, 0) + amount

M, n, z = sp.symbols("M n z")
def reduced_derivative(v):
    a, b, c = v
    return [
        sp.diff(a,z)-2*n*c/z,
        a+sp.diff(b,z)+(1-M*(M-1)/z**2)*c,
        b+sp.diff(c,z)+(2*M-2)*c/z,
    ]

v = [sp.Integer(1),sp.Integer(0),sp.Integer(0)]
vectors = {0:v}
for j in range(6):
    dv = reduced_derivative(v)
    v = [sp.factor((M+j)*v[k]-z*dv[k]) for k in range(3)]
    vectors[j+1] = v
remainder = [sp.factor(z**2*vectors[4][k]-vectors[6][k]) for k in range(3)]
A = (M+2)*(M+3)*(M+4)*(M+5)
s = 2*(n+1)
U = A+s*(2*M+7)*z**2
W = s*z**2*((M+4)*(M+5)+z**2)
q = -M*(M+1)*A+(M+2)*(M+3)*(M*(M+1)-4*n*(2*M+5))*z**2-4*n*(n+1)*z**4
b, c = z*(2*M*U+W), -z**2*U
for got, wanted in zip(remainder,[q,b,c]):
    assert sp.expand(got-wanted) == 0
    count("generic_full_Ore_remainder_coefficients")
# Both gauge-stripped expressions, with B and B' independent.
for sign in [1,-1]:
    derivative_coefficient = sign*c
    constant_coefficient = b-sign*c+(M-1)*c/z
    wanted_derivative = -sign*z**2*U
    wanted_constant = z*((sign*z+M+1)*U+W)
    assert sp.expand(derivative_coefficient-wanted_derivative) == 0
    assert sp.expand(constant_coefficient-wanted_constant) == 0
    count("generic_stripped_gauge_coefficients",2)
# Corrupting a genuinely essential shift-four coefficient is detected.
assert sp.expand(remainder[2]-(-z**2*(A+s*(2*M+5)*z**2))) != 0
count("corrupted_generic_coefficient_rejections")

def moment_R(index, gamma):
    return sp.Poly(sum(
        (-1)**j*sp.binomial(index,j)*sp.rf(gamma,2*j)*z**(2*(index-j))
        for j in range(index+1)),z).as_expr()
for index in range(1,7):
    gamma = 4*index+1
    MM = gamma+2*index
    R = moment_R(index,gamma)
    E = lambda f,j: (MM+j)*f-z*sp.diff(f,z)
    vv, seq = R,{0:R}
    for j in range(6):
        vv = sp.expand(E(vv,j))
        seq[j+1] = vv
    lhs = sp.rf(gamma,4)*moment_R(index+1,gamma+4)
    rhs = z**2*seq[4]-seq[6]
    delta = sp.Poly(sp.expand(lhs-rhs),z)
    for exponent in range(2*(index+1)+1):
        assert delta.nth(exponent) == 0
        count("literal_shape_shift_coefficients")
    count("literal_shape_shift_identities")
    ode = z**2*sp.diff(R,z,3)+(2-2*MM)*z*sp.diff(R,z,2)+(MM*(MM-1)-z**2)*sp.diff(R,z)+2*index*z*R
    assert sp.expand(ode) == 0
    count("literal_cubic_ODE_identities")

# Complete reduced operator on k[z]/z^2 for actual gamma odd.
A0 = [[0,0],[0,0]]
for degree in range(2):
    A0[degree][degree] = (-degree)%2
    if degree+1<2:
        A0[degree+1][degree] = 1
A02 = [[sum(A0[i][k]*A0[k][j] for k in range(2))%2 for j in range(2)] for i in range(2)]
for i in range(2):
    for j in range(2):
        assert A02[i][j] == A0[i][j]
        count("complete_reduced_idempotent_entries")
basis = [[1,A0[0][0]],[0,A0[1][0]]]
assert (basis[0][0]*basis[1][1]-basis[0][1]*basis[1][0])%2 == 1
count("complete_cyclic_basis_determinants")
for vector in [[A0[0][0],A0[1][0]],
               [(A0[0][0]*0+A0[0][1]*1)%2,(A0[1][0]*0+A0[1][1]*1)%2]]:
    assert vector == [0,1]
    count("complete_positive_power_seed_vectors")

# Exact valuation differences; the written inequalities use L>=1 and
# respectively L<a+1 or L>a+1. Boundary coefficient gains are positive.
aa, bb, L = sp.symbols("a b0 L")
ee = aa+bb+2
w1, w2 = 2*bb+2+L, bb+1+2*L
for got,wanted in [
    (w1-(bb+1+L),bb+1),(w2-(bb+1+L),L),
    (w1-ee,bb+L-aa),(w2-ee,2*L-aa-1),
    (w1.subs(L,aa+1)-ee,bb+1),(w2.subs(L,aa+1)-ee,aa+1),
]:
    assert sp.expand(got-wanted) == 0
    count("symbolic_valuation_difference_identities")

def trim(p):
    while len(p)>1 and p[-1]==0:
        p.pop()
    return p
def add(p,q):
    r = [0]*max(len(p),len(q))
    for i,v in enumerate(p): r[i] += v
    for i,v in enumerate(q): r[i] += v
    return trim(r)
def mul(p,q):
    r = [0]*(len(p)+len(q)-1)
    for i,a in enumerate(p):
        for j,b in enumerate(q):
            r[i+j] += a*b
    return trim(r)
def scale_variable(p,c):
    return [a*c**i for i,a in enumerate(p)]
def auxiliary(index,gamma):
    p = [1]
    for j in range(index):
        r = [0]*(len(p)+1)
        for k,a in enumerate(p):
            r[k] += (2*j+gamma+1-k)*a
            r[k+1] += a
        p = trim(r)
    assert p[-1] == 1
    return p
def valuation2(x):
    assert x != 0
    return (abs(x)&-abs(x)).bit_length()-1
def affine_reduction_ok(p):
    return len(p)>=2 and [c%2 for c in p] == [1,1]+[0]*(len(p)-2)

profiles = []
for index in [1,2,3,4,6,10,18,42,85,127]:
    gamma = 4*index+1
    MM, ss = 6*index+1,2*(index+1)
    AA = (MM+2)*(MM+3)*(MM+4)*(MM+5)
    aa, bb = valuation2(3*index+2),valuation2(index+1)
    ee, q0 = aa+bb+2,2**(aa+1)
    assert valuation2(AA) == ee
    assert valuation2(ss) == bb+1
    assert valuation2(MM+4) == 0
    assert valuation2(MM+5) == bb+1
    assert valuation2(2*MM+7) == 0
    count("literal_parameter_valuation_identities",5)
    B,C = auxiliary(index,gamma),auxiliary(index+1,gamma)
    B1,C1 = B[1::2] or [0],C[1::2] or [0]
    Up = [AA,ss*(2*MM+7)]
    Wp = [0,ss*(MM+4)*(MM+5),ss]
    numerator = add(
        mul(scale_variable(Up,q0),scale_variable(C1,q0)),
        mul(scale_variable(Wp,q0),scale_variable(B1,q0)))
    assert all(c%(2**ee)==0 for c in numerator)
    F = [c//2**ee for c in numerator]
    assert affine_reduction_ok(F)
    # Independently form H(z) before extracting the odd component.
    Uz = [AA,0,ss*(2*MM+7)]
    Wz = [0,0,ss*(MM+4)*(MM+5),0,ss]
    Hz = add(mul(Uz,C),mul(Wz,B))
    odd = scale_variable(Hz[1::2],q0)
    assert odd == numerator
    count("full_literal_odd_scalar_polynomials")
    count("literal_scalar_coefficients_checked",len(F))
    profiles.append({
        "n":index,"a":aa,"b0":bb,"e":ee,"q0":q0,"scalar_degree":len(F)-1,
        "whole_affine_reduction":[1,1],"all_coefficients_integral":True,
        "odd_section_before_and_after_substitution_agree":True,
        "full_integer_scalar_coefficients_sha256":
            hashlib.sha256(json.dumps(F,separators=(",",":")).encode()).hexdigest(),
    })
for corrupted in [[1,0,2],[1,1,1]]:
    assert not affine_reduction_ok(corrupted)
    count("corrupted_scalar_coefficient_rejections")

receipt = {
    "kind":"internal exact replay; all-index quantifiers are proved in the retained TeX",
    "passed":True,
    "scope":{
        "ordinary_support_four_six_all_indices":True,
        "all_positive_Gamma_shapes":False,
        "all_degree_six_supports":False,
        "all_ordinary_univariate_binomials":False,
        "unrestricted_FC_or_SFC":False,
    },
    "controls":counts,
    "full_scalar_profiles":profiles,
    "exact_imports":{
        "simplicity":"CUBIC_ADJOINT_STABILITY_SFC.tex thm:simple",
        "zero_concomitant":"CUBIC_ADJOINT_STABILITY_SFC.tex equation(10)",
        "operator_normalization":"lower_exponent_initial_and_boundary.tex equations(G1)-(G2)",
    },
    "source_sha256":{
        PROOF.name:hashlib.sha256(PROOF.read_bytes()).hexdigest(),
        Path(__file__).name:hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        **IMPORTS,
    },
    "python":platform.python_version(),"sympy":sp.__version__,
    "independent_replay":False,"formal_verification":False,"external_review":False,
    "novelty_claim":False,"priority_claim":False,
    "command":".venv/bin/python research/factorial/image-rigidity/lower_four_quadratic_gap_sfc.py",
}
out = HERE/"lower_four_quadratic_gap_sfc_verification.json"
out.write_text(json.dumps(receipt,indent=2)+"\n")
print(json.dumps({"passed":True,"controls":counts,"receipt":str(out)},indent=2))

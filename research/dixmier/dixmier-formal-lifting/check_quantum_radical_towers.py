#!/usr/bin/env python3
"""Exact controls for QRG1; the internal-tower theorem is a written proof."""
from pathlib import Path
from math import comb
import hashlib
import json
import platform
import tempfile
import sympy as s

HERE = Path(__file__).resolve().parent
h, r, a, ss, t, u, x = s.symbols("h r a s t u x")

def clean(op):
    return {key: s.factor(val) for key, val in op.items()
            if s.factor(val) != 0}

def add(*ops):
    out = {}
    for op in ops:
        for key, value in op.items():
            out[key] = out.get(key, 0)+value
    return clean(out)

def scale(op, c):
    return clean({key: c*value for key, value in op.items()})

def mul(left, right):
    """Laurent Ore normal product [b,v]=h; keys (v exponent,b order)."""
    out = {}
    for (e, i), c in left.items():
        for (f, j), d in right.items():
            derivative = s.Integer(1)
            for k in range(i+1):
                key = (e+f-k, i+j-k)
                value = c*d*comb(i,k)*derivative*h**k
                out[key] = out.get(key, 0)+value
                derivative *= f-k
    return clean(out)

def power(op, k):
    out = {(0,0): s.Integer(1)}
    for _ in range(k):
        out = mul(out, op)
    return out

def comm(left, right):
    return add(mul(left,right), scale(mul(right,left),-1))

def h_coefficient(op, n):
    return clean({key: s.expand(value).coeff(h,n)
                  for key, value in op.items()})

def serial(op):
    return [[e,j,str(value)] for (e,j),value in sorted(op.items())]

# Full classical Long extension: discriminants and squarefree certificate.
f = r*a**3-2*a*a+ss*a+4*t
delta = s.discriminant(f,a)
expected = 4*ss**2-4*r*ss**3+128*t-144*r*ss*t-432*r*r*t*t
assert s.expand(delta-expected) == 0
disc_t = s.factor(s.discriminant(delta,t))
assert s.factor(disc_t-256*(4-3*r*ss)**3) == 0
field = s.QQ.frac_field(r,ss)
p_delta = s.Poly(delta,t,domain=field)
assert s.gcd(p_delta,p_delta.diff()).degree() == 0
assert p_delta.degree() == 2

# Actual finite positive controls: monomial coverings, with all mates finite Laurent.
coverings = []
for k in range(2,8):
    P = {(k,0):s.Integer(1)}
    Q = {(1-k,1):s.Rational(1,k)}
    assert comm(Q,P) == {(0,0):h}
    assert power({(1,0):s.Integer(1)},k) == P
    assert scale(mul({(k-1,0):s.Integer(1)},Q),k) == {(0,1):s.Integer(1)}
    # A parameter-pole shear of b; its radical presentation is still actual.
    Qpole = add(Q,{(4-k,0):s.Rational(1,k)/h**3})
    assert comm(Qpole,P) == {(0,0):h}
    coverings.append({
        "k":k, "P":serial(P), "Q":serial(Q),
        "parameter_pole_shear_Q":serial(Qpole),
        "canonical_commutator":"h",
        "internal_root":"v^k=P",
        "ambient_b_recovery":"b=k*v^(k-1)*Q",
    })

# Noncommutative Hensel control: A=(v+h*b)^2 and first correction from v.
v = {(1,0):s.Integer(1)}
b = {(0,1):s.Integer(1)}
theta = add(v,scale(b,h))
A = power(theta,2)
error = add(A,scale(power(v,2),-1))
first_error = h_coefficient(error,1)
assert first_error == {(1,1):s.Integer(2)}
assert power(add(v,scale(b,h)),2) == A
assert comm(v,b) != {}
assert comm(A,v) != {}  # The radicand need not be central.

# The first differing coefficient of two kth powers is k*v^(k-1)*U.
U = {(2,1):s.Integer(1),(-1,0):s.Integer(1)}
uniqueness = []
for k in range(2,7):
    for n in range(1,4):
        Y = add(theta,scale(U,h**n))
        difference = add(power(Y,k),scale(power(theta,k),-1))
        expected_leading = scale(mul({(k-1,0):s.Integer(1)},U),k)
        assert h_coefficient(difference,n) == expected_leading
        assert all(h_coefficient(difference,j) == {} for j in range(n))
        uniqueness.append({"k":k,"first_difference_order":n,
                           "leading_coefficient":"k*v^(k-1)*U","passed":True})

# Actual quadratic extension hidden by completion; no false descent to E.
jets = []
for N in range(2,11):
    root = sum(s.binomial(s.Rational(1,2),j)*(4*h*u)**j
               for j in range(N+1))
    square_error = s.series(root*root-(1+4*h*u),h,0,N+1).removeO()
    identity_error = s.series(root.subs(u,x+h*x*x)-(1+2*h*x),
                             h,0,N+1).removeO()
    assert s.expand(square_error) == 0
    assert s.expand(identity_error) == 0
    jets.append(N)
assert s.expand((1+2*h*x)**2-(1+4*h*(x+h*x*x))) == 0
assert s.diff(1+4*h*u,u) == 4*h

# Essential residue-commutativity boundary.
X = s.diag(1,-1)
N12 = s.Matrix([[0,1],[0,0]])
Y = X+h*N12
assert X*X == s.eye(2) and Y*Y == s.eye(2)
assert (Y-X).subs(h,0) == s.zeros(2)
assert Y != X and X*N12+N12*X == s.zeros(2)

names = [
    "QUANTUM_RADICAL_TOWERS_AND_LONG_OBSTRUCTION.md",
    "check_quantum_radical_towers.py",
    "QUANTUM_DIVISION_AND_COMPLETION.md",
]
receipt = {
    "claim":"QRG1", "status":"PASS",
    "environment":{"python":platform.python_version(),"sympy":s.__version__},
    "Long_cubic":{
        "polynomial":str(f),
        "discriminant":str(s.factor(delta)),
        "discriminant_as_quadratic_in_t":str(disc_t),
        "delta_squarefree_in_t":True,
        "geometric_normal_closure":"S3, by written cubic argument",
    },
    "actual_finite_monomial_covers":coverings,
    "noncommutative_Hensel":{
        "root":"theta=v+h*b", "radicand":serial(A),
        "initial_lift":"v", "first_correction":"h*b",
        "radicand_noncentral":True,
        "finite_uniqueness_controls":uniqueness,
    },
    "completion_boundary":{
        "u":"x+h*x^2", "theta":"1+2*h*x", "theta_squared":"1+4*h*u",
        "binomial_jet_orders":jets,
        "actual_extension":"quadratic",
        "residue_extension_degree":1,
        "conclusion":"Root lies in completed image, not necessarily actual image.",
    },
    "noncommutative_residue_boundary":{
        "X":"diag(1,-1)", "Y":"X+h*e12",
        "squares_equal":"I", "same_residue":True, "distinct":True,
        "conclusion":"Commutative residue is essential to root uniqueness.",
    },
    "hashes":{name:hashlib.sha256((HERE/name).read_bytes()).hexdigest()
              for name in names},
    "coverage_boundary":(
        "Exact cubic discriminants and squarefree certificate; finite monomial "
        "coverings and parameter-pole shears; noncentral Weyl root and finite "
        "first-difference controls; binomial completion jets and matrix residue "
        "countercontrol. Universal noncommutative Hensel existence/uniqueness, "
        "internal tower induction, scalar descent and all-coordinate-change "
        "transport use the written proof. No mixed Long mate or bounded lift."
    ),
    "independent_replay":False,"formal_verification":False,"external_review":False,
}
out=Path(tempfile.mkdtemp(prefix="qrg1-check-"))/"quantum_radical_tower_controls.json"
out.write_text(json.dumps(receipt,indent=2)+"\n")
print(json.dumps({"claim":"QRG1","status":"PASS","coverings":len(coverings),
                  "uniqueness_cases":len(uniqueness),"receipt":str(out)},indent=2))


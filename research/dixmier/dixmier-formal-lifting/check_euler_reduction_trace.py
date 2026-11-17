#!/usr/bin/env python3
"""Exact finite controls for QER1; the universal trace theorem is written."""
from itertools import product
from math import comb, factorial
from pathlib import Path
import hashlib
import json
import platform
import tempfile
import sympy as s

HERE = Path(__file__).resolve().parent
h, eta, lam, z, q = s.symbols("h eta lam z q")
W = (1, 1, 2)
ZERO = (0, 0, 0)

def falling(n, k):
    return factorial(n) // factorial(n-k) if n >= k else 0

def clean(op):
    out = {}
    for key, val in op.items():
        val = s.factor(val)
        if val != 0:
            out[key] = val
    return out

def add(*ops):
    out = {}
    for op in ops:
        for key, val in op.items():
            out[key] = out.get(key, 0) + val
    return clean(out)

def scale(op, c):
    return clean({key: c*val for key, val in op.items()})

def monomial(a=ZERO, b=ZERO):
    return {(tuple(a), tuple(b)): s.Integer(1)}

def mul(left, right):
    out = {}
    for (a, b), c in left.items():
        for (d, e), f in right.items():
            for k in product(*(range(min(b[i], d[i])+1) for i in range(3))):
                aa = tuple(a[i]+d[i]-k[i] for i in range(3))
                bb = tuple(b[i]+e[i]-k[i] for i in range(3))
                v = c*f*h**sum(k)
                for i in range(3):
                    v *= comb(b[i], k[i])*falling(d[i], k[i])
                out[(aa, bb)] = out.get((aa, bb), 0)+v
    return clean(out)

def comm(left, right):
    return add(mul(left, right), scale(mul(right, left), -1))

def tau_diagonal(a):
    """8 times the (-1)^N coefficient polynomial from the pole at t=-1."""
    aa = sum(a)
    M = a[0]+a[1]+2*a[2]
    B = a[2]+1
    F = h**aa
    for i in a:
        F *= factorial(i)
    F *= (q-1)**M/(2-q)**(aa+3)
    out = 0
    for k in range(1, B+1):
        j = B-k
        ck = s.diff(F, q, j).subs(q, 0)/factorial(j)
        # The binomial is polynomial in z; expand without assumptions.
        bk = s.prod(z+i for i in range(1, k))/s.Integer(factorial(k-1))
        out += ck*bk
    return s.factor(8*out)

def tau(op):
    out = 0
    for (a, b), coeff in op.items():
        assert sum(W[i]*(a[i]-b[i]) for i in range(3)) == 0
        if a == b:
            out += coeff*tau_diagonal(a).subs(z, eta/h)
    return s.factor(out)

def finite_diagonal_trace(a, N):
    total = 0
    for n3 in range(N//2+1):
        for n1 in range(N-2*n3+1):
            n2 = N-2*n3-n1
            total += falling(n1, a[0])*falling(n2, a[1])*falling(n3, a[2])
    return total

# Independent interpolation from finite-dimensional weighted polynomial modules.
continuation = []
for a in product(range(3), repeat=3):
    degree = sum(a)+2  # pole at +1 has order sum(a)+3
    p0 = s.interpolate([(2*j, finite_diagonal_trace(a, 2*j))
                        for j in range(degree+1)], z)
    p1 = s.interpolate([(2*j+1, finite_diagonal_trace(a, 2*j+1))
                        for j in range(degree+1)], z)
    actual = s.factor(4*(p0-p1)*h**sum(a))
    predicted = tau_diagonal(a)
    assert s.factor(actual-predicted) == 0, a
    for N in [2*degree+2, 2*degree+3]:
        expected = (p0 if N % 2 == 0 else p1).subs(z, N)
        assert expected == finite_diagonal_trace(a, N)
    continuation.append({'alpha': list(a), 'normalized_trace': str(predicted)})

one = monomial()
Hs = [monomial(tuple(int(i == j) for i in range(3)),
               tuple(int(i == j) for i in range(3))) for j in range(3)]
E = add(Hs[0], Hs[1], scale(Hs[2], 2))
assert tau(one) == 1
assert tau(Hs[0]) == -h/2 and tau(Hs[1]) == -h/2
assert s.factor(tau(Hs[2])-(eta+h)/2) == 0
assert tau(mul(Hs[0], Hs[0])) == 0
assert tau(mul(Hs[0], Hs[1])) == h*h/4

generators = [
    monomial((1,0,0), (0,1,0)), monomial((0,1,0), (1,0,0)),
    monomial((2,0,0), (0,0,1)), monomial((0,0,1), (2,0,0)),
    monomial((0,2,0), (0,0,1)), monomial((0,0,1), (0,2,0)),
    monomial((1,1,0), (0,0,1)), monomial((0,0,1), (1,1,0)),
] + Hs
nonzero_commutators = 0
for i, U in enumerate(generators):
    for V in generators[i+1:]:
        C = comm(U, V)
        if C:
            nonzero_commutators += 1
        assert tau(C) == 0
for D in generators + [one, mul(Hs[0], Hs[1]), mul(generators[2], generators[3])]:
    assert tau(add(mul(E, D), scale(D, -eta))) == 0

# Fourier moment shift is checked by the actual normal relation p*x=x*p+h.
# Fourier(x*p)=p*(-x)=-(x*p+h).
fourier_negative_terms = add(Hs[0], add(Hs[1], scale(one, h)),
                             scale(add(Hs[2], scale(one, h)), 2))
assert fourier_negative_terms == add(E, scale(one, 3*h))

# Retained normalized Keller map: weights, determinant and exact Euler row.
x, y, zz = s.symbols("x y zz")
u0 = 1+x*y
f = s.Matrix([
    (2*x-3*x*x*y-x**3*zz)/2,
    y+3*x*u0**2*zz+3*x*y*y*(4+3*x*y),
    u0**3*zz+y*y*u0*(4+3*x*y),
])
weights = (1,-1,-2)
variables = (x,y,zz)
for i, fi in enumerate(f):
    for exponents, coeff in s.Poly(fi, variables).terms():
        assert sum(weights[j]*exponents[j] for j in range(3)) == weights[i]
J = f.jacobian(variables)
assert s.factor(J.det()) == 1
weighted_source = s.Matrix([x,-y,-2*zz])
weighted_target = s.Matrix([f[0],-f[1],-2*f[2]])
assert all(s.expand(t) == 0 for t in J*weighted_source-weighted_target)
assert all(s.factor(t) == 0 for t in
           (weighted_target.T*J.adjugate().T-weighted_source.T))

# Position quotient and its nonunit Jacobian.
u, v = s.symbols("u v")
A = (2-3*u-v)/2
B = u+3*v*(1+u)**2+3*u*u*(4+3*u)
C = v*(1+u)**3+u*u*(1+u)*(4+3*u)
quotient = s.Matrix([A*B,A*A*C])
assert s.factor(quotient.jacobian((u,v)).det()-A*A) == 0
assert s.factor(f[0]-x*A.subs({u:x*y, v:x*x*zz})) == 0
assert s.factor(f[0]*f[1]-(A*B).subs({u:x*y,v:x*x*zz})) == 0
assert s.factor(f[0]**2*f[2]-(A*A*C).subs({u:x*y,v:x*x*zz})) == 0

names = [
    "EULER_REDUCTION_TRACE_OBSTRUCTION.md",
    "check_euler_reduction_trace.py",
    "IDENTITY_LIFTS_AND_UNBOUNDED_INVERSES.md",
]
receipt = {
    "claim": "QER1",
    "status": "PASS",
    "environment": {"python": platform.python_version(), "sympy": s.__version__},
    "trace_continuations": continuation,
    "finite_trace_interpolation": {
        "alpha_box": "0<=alpha_i<=2",
        "cases": len(continuation),
        "independent_finite_sums": True,
        "interpolation_degree_bound": "|alpha|+2 from pole at +1",
        "additional_points_per_case": 2,
    },
    "sparse_Weyl_trace_controls": {
        "invariant_operators": len(generators),
        "pair_commutators": len(generators)*(len(generators)-1)//2,
        "nonzero_commutators_with_trace_zero": nonzero_commutators,
        "Euler_level_relation": "PASS",
        "normalized_identity": "1",
        "H1": "-h/2", "H2": "-h/2", "H3": "(eta+h)/2",
        "H1_squared": "0", "H1_H2": "h^2/4",
    },
    "cotangent_descent_controls": {
        "weights": [1,-1,-2],
        "normalized_Keller_determinant": "1",
        "weighted_Euler_identity": "PASS",
        "inverse_Jacobian_moment_preservation": "PASS",
        "Fourier_moment": "E_plus+3*h",
        "original_to_positive_level": "eta=lambda-3*h",
        "quotient_positions": ["A*B", "A^2*C"],
        "quotient_Jacobian": "(2-3*u-v)^2/4=A^2",
        "localized_first_position": "x*A",
    },
    "hashes": {name: hashlib.sha256((HERE/name).read_bytes()).hexdigest()
               for name in names},
    "coverage_boundary": (
        "Exact finite trace interpolation for the alpha box, sparse invariant "
        "Weyl commutators and Euler level relation, normalized Keller weights "
        "and moment identity, Fourier shift, and quotient Jacobian. The all-"
        "operator continuation, all scalar levels, matrix and full-corner "
        "obstructions, global chart faithfulness and reduced nonsurjectivity "
        "are written proofs. No bounded Weyl lift or DC1/DC2 resolution."
    ),
    "independent_replay": False,
    "formal_verification": False,
    "external_review": False,
}
out = Path(tempfile.mkdtemp(prefix="qer1-check-"))/"euler_reduction_controls.json"
out.write_text(json.dumps(receipt, indent=2)+"\n")
print(json.dumps({"claim": "QER1", "status": "PASS",
                  "trace_cases": len(continuation),
                  "commutators": len(generators)*(len(generators)-1)//2,
                  "receipt": str(out)}, indent=2))


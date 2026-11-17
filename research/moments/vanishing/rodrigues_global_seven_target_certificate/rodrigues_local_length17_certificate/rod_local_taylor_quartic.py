"""Exact local sphere-moment jets; scratch analysis, not a deliverable.

q=2*z*w+t**2=1, F0=w+z*t-z**3/2.
Solve moments 1..10 formally on an explicit transverse slice through degree 3.
Then compute residual moments 11..14 through ordinary degree 4.
No quotient multiplication table, numerical approximation, or weighted initial
ideal is used.
"""
import sympy as sp
from functools import lru_cache

z, t, w, a, b = sp.symbols("z t w a b")
F0 = w + z*t - z**3/2
basis = [1, z, t, w, z*z, z*t, t*t-z*w, w*t, w*w,
         z**3, z*z*t, z*(t*t-z*w/2), t**3-3*z*w*t,
         w*(t*t-z*w/2), w*w*t, w**3]

def poly_dict(f):
    return sp.Poly(sp.expand(f), z, t, w).as_dict()

powers = [poly_dict(F0**n) for n in range(15)]

@lru_cache(None)
def marked_monomial(n, iz, it, iw):
    out = sp.S.Zero
    for (jz, jt, jw), coeff in powers[n].items():
        Z, T, W = iz+jz, it+jt, iw+jw
        if Z == W and T % 2 == 0:
            out += (coeff * sp.factorial(Z) * sp.factorial2(T-1)
                    / sp.factorial2(2*Z+T+1))
    return out

@lru_cache(None)
def H(n, marker):
    return sum(coeff * marked_monomial(n, *ijk)
               for ijk, coeff in poly_dict(marker).items())

def parameter_H(n, marker):
    return sum(H(n, coeff) * a**i * b**j
               for (i, j), coeff in
               sp.Poly(sp.expand(marker), a, b).as_dict().items())

Qa = w*w - sp.Rational(1, 10)*(t**3-3*z*w*t) + sp.Rational(17, 70)*t
Qb = w*t + sp.Rational(2, 5)*z*(t*t-z*w/2) - z/5
T1 = a*Qa+b*Qb

J = sp.Matrix([[n*H(n-1, q) for q in basis] for n in range(1, 11)])
_, columns = J.rref()
assert columns == (0, 1, 2, 3, 6, 7, 8, 13, 14, 15)
directions = [basis[i] for i in columns]
J0 = J[:, columns]
Jinv = J0.inv()
assert all(n*H(n-1, Q) == 0 for n in range(1, 11) for Q in (Qa, Qb))
print("Jacobian determinant:", sp.factor(J0.det()))

v2 = sp.Matrix([sp.binomial(n, 2)*parameter_H(n-2, T1*T1)
                if n >= 2 else 0 for n in range(1, 11)])
T2 = sp.expand(sum(c*q for c, q in zip(-Jinv*v2, directions)))
v3 = sp.Matrix([
    (sp.binomial(n, 3)*parameter_H(n-3, T1**3) if n >= 3 else 0)
    + (n*(n-1)*parameter_H(n-2, T1*T2) if n >= 2 else 0)
    for n in range(1, 11)])
T3 = sp.expand(sum(c*q for c, q in zip(-Jinv*v3, directions)))

# The defining equalities below directly certify the implicit corrections.
assert all(sp.expand(n*parameter_H(n-1, T2) + v2[n-1]) == 0
           for n in range(1, 11))
assert all(sp.expand(n*parameter_H(n-1, T3) + v3[n-1]) == 0
           for n in range(1, 11))

for n in (11, 12, 13, 14):
    # The unknown T4 term is multiplied by D mu_n(F0), which vanishes.
    assert all(n*H(n-1, q) == 0 for q in basis)
    r2 = sp.factor(sp.binomial(n, 2)*parameter_H(n-2, T1*T1))
    r3 = sp.factor(sp.binomial(n, 3)*parameter_H(n-3, T1**3)
                   + n*(n-1)*parameter_H(n-2, T1*T2))
    r4 = sp.factor(sp.binomial(n, 4)*parameter_H(n-4, T1**4)
                   + 3*sp.binomial(n, 3)*parameter_H(n-3, T1*T1*T2)
                   + sp.binomial(n, 2)*parameter_H(n-2, 2*T1*T3+T2*T2))
    assert r2 == r3 == 0
    if n in (13, 14):
        assert r4 == 0
    print(f"R{n}: quadratic={r2}; cubic={r3}; quartic={r4}")

print("T2:", sp.collect(T2, [a*a, a*b, b*b]))
print("T3:", [(ij, str(sp.factor(c)))
              for ij, c in sp.Poly(T3, a, b).as_dict().items()])

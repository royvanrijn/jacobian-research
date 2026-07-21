"""Exact relative Groebner classification of the three-transfer block Z_3."""

import sympy as sp

Z = sp.symbols("Z")
p, q, r, X, Y, T = sp.symbols("p q r X Y T")
u = sp.symbols("u0:9")

S = Z**3+p*Z**2+q*Z+r
V = sp.expand(S**2+X*Z**2+Y*Z+T)
U = Z**9+sum(u[index]*Z**index for index in range(9))
difference = sp.Poly(sp.expand(U**2-V**3), Z)

# Degrees 17,...,9 eliminate the nine coefficients of U triangularly.
solution = {}
for degree, variable in zip(range(17, 8, -1), reversed(u)):
    equation = sp.expand(difference.coeff_monomial(Z**degree).subs(solution))
    roots = sp.solve(equation, variable, dict=False)
    assert len(roots) == 1
    solution[variable] = sp.factor(roots[0])

remaining = [
    sp.factor(difference.coeff_monomial(Z**degree).subs(solution))
    for degree in range(8, -1, -1)
]

# Affine difference omits degrees one and zero.  They are nevertheless
# consequences of the degree 8,...,2 equations.
coefficient_domain = sp.QQ.poly_ring(p, q, r)
affine_groebner = sp.groebner(
    remaining[:7], T, Y, X, order="grevlex", domain=coefficient_domain
)
strong_groebner = sp.groebner(
    remaining, T, Y, X, order="grevlex", domain=coefficient_domain
)
assert all(affine_groebner.reduce(equation)[1] == 0 for equation in remaining)
assert len(affine_groebner.polys) == len(strong_groebner.polys) == 7

expected = [
    T**3-6*T**2*r**2-6*X**2*p*r**3+12*X*Y*r**3,
    T**2*X-2*T**2*p*r-2*X**2*p**2*r**2+4*X*Y*p*r**2,
    T*X**2-2*T**2*q-2*X**2*p*q*r+4*X*Y*q*r,
    X**2*Y-4*T**2*p-4*X**2*p**2*r+8*X*Y*p*r,
    X**3-6*T**2-6*X**2*p*r+12*X*Y*r,
    2*T*Y+X**2*p*q-X**2*r-2*X*Y*q,
    Y**2+2*T*X+X**2*p**2-X**2*q-2*X*Y*p,
]
expected_groebner = sp.groebner(
    expected, T, Y, X, order="grevlex", domain=coefficient_domain
)
assert all(expected_groebner.reduce(equation)[1] == 0
           for equation in remaining[:7])
assert all(affine_groebner.reduce(equation)[1] == 0 for equation in expected)

# Monic leading terms give a free rank-eight module over k[p,q,r].
basis = (
    sp.Integer(1), T, T**2, X, T*X, X**2, Y, X*Y,
)
assert all(expected_groebner.reduce(element)[1] == element for element in basis)

# Coincident-root fiber and its two independent socle elements T^2 and X*Y.
special_relations = [equation.subs({p: 0, q: 0, r: 0}) for equation in expected]
special = sp.groebner(special_relations, T, Y, X, order="grevlex", domain=sp.QQ)
assert all(special.reduce(element)[1] == element for element in basis)
for socle_element in (T**2, X*Y):
    for generator in (T, X, Y):
        assert special.reduce(socle_element*generator)[1] == 0

print("PASS: the nine high coefficients eliminate U triangularly")
print("PASS: affine difference and strong equality give the same Z_3 ideal")
print("PASS: Z_3 is finite flat of rank 8 over k[p,q,r]")
print("PASS: its coincident fiber has Hilbert function (1,3,3,1) and socle dimension 2")

"""Exact relative Groebner classification of the four-transfer block Z_4."""

import sympy as sp

Z = sp.symbols("Z")
p, q, r, t = sp.symbols("p q r t")
A, B, C, D = sp.symbols("A B C D")
u = sp.symbols("u0:12")

S = Z**4+p*Z**3+q*Z**2+r*Z+t
V = sp.expand(S**2+A*Z**3+B*Z**2+C*Z+D)
U = Z**12+sum(u[index]*Z**index for index in range(12))
difference = sp.Poly(sp.expand(U**2-V**3), Z)

# Degrees 23,...,12 eliminate all twelve coefficients of U triangularly.
solution = {}
for degree, variable in zip(range(23, 11, -1), reversed(u)):
    equation = sp.expand(difference.coeff_monomial(Z**degree).subs(solution))
    roots = sp.solve(equation, variable, dict=False)
    assert len(roots) == 1
    solution[variable] = sp.factor(roots[0])

remaining = [
    sp.factor(difference.coeff_monomial(Z**degree).subs(solution))
    for degree in range(11, -1, -1)
]

coefficient_domain = sp.QQ.poly_ring(p, q, r, t)
affine_groebner = sp.groebner(
    remaining[:10], D, C, B, A, order="grevlex", domain=coefficient_domain
)
strong_groebner = sp.groebner(
    remaining, D, C, B, A, order="grevlex", domain=coefficient_domain
)
assert len(affine_groebner.polys) == len(strong_groebner.polys) == 13
assert all(affine_groebner.reduce(equation)[1] == 0 for equation in remaining)

expected_leading_exponents = (
    (1, 1, 0, 2), (3, 0, 0, 0), (2, 1, 0, 0),
    (2, 0, 1, 0), (1, 0, 2, 0), (0, 0, 3, 0),
    (2, 0, 0, 1), (1, 0, 1, 1), (0, 0, 2, 1),
    (0, 0, 1, 2), (0, 0, 0, 3), (0, 2, 0, 0),
    (0, 1, 1, 0),
)
assert tuple(
    polynomial.LM(order=affine_groebner.order).exponents
    for polynomial in affine_groebner.polys
) == expected_leading_exponents
assert all(polynomial.LC(order=affine_groebner.order) == 1
           for polynomial in affine_groebner.polys)

# The monic leading terms are constant over the S-parameter space.
standard_basis = tuple(map(sp.sympify, (
    1, A, B, C, D,
    A**2, A*B, B**2, A*C, A*D, B*D, C*D, D**2,
    A**2*C, A**2*D, A*C*D,
)))
assert len(standard_basis) == 16
assert all(affine_groebner.reduce(element)[1] == element
           for element in standard_basis)

# Analyze the coincident-root fiber p=q=r=t=0.
special_relations = [
    polynomial.as_expr().subs({p: 0, q: 0, r: 0, t: 0})
    for polynomial in affine_groebner.polys
]
special = sp.groebner(special_relations, D, C, B, A,
                      order="grevlex", domain=sp.QQ)
assert all(special.reduce(element)[1] == element for element in standard_basis)

monomials = [sp.Poly(element, D, C, B, A).monoms()[0]
             for element in standard_basis]
index = {monomial: position for position, monomial in enumerate(monomials)}


def coordinate_vector(polynomial):
    vector = sp.zeros(16, 1)
    remainder = sp.Poly(special.reduce(polynomial)[1], D, C, B, A)
    for monomial, coefficient in remainder.terms():
        vector[index[monomial]] = coefficient
    return vector


multiplication = [
    sp.Matrix.hstack(*[
        coordinate_vector(generator*element) for element in standard_basis
    ])
    for generator in (A, B, C, D)
]

# Socle is the common kernel of multiplication by the maximal ideal.
socle_dimension = 16-sp.Matrix.vstack(*multiplication).rank()
assert socle_dimension == 4

# Compute dimensions of successive maximal-ideal powers.
power_space = sp.Matrix.hstack(*[
    coordinate_vector(element) for element in standard_basis[1:]
])
power_dimensions = []
while True:
    power_space = sp.Matrix.hstack(*power_space.columnspace())
    power_dimensions.append(power_space.rank())
    if power_space.rank() == 0:
        break
    power_space = sp.Matrix.hstack(*[
        matrix*column
        for matrix in multiplication
        for column in power_space.columnspace()
    ])
assert power_dimensions == [15, 11, 5, 1, 0]
assert (1, 4, 6, 4, 1) == (
    16-power_dimensions[0],
    *(power_dimensions[index]-power_dimensions[index+1]
      for index in range(len(power_dimensions)-1)),
)

print("PASS: the twelve high coefficients eliminate U triangularly")
print("PASS: affine difference and strong equality give the same Z_4 ideal")
print("PASS: Z_4 is finite flat of rank 16 over the monic quartic S-space")
print("PASS: its coincident fiber has Hilbert function (1,4,6,4,1) and socle dimension 4")

#!/usr/bin/env python3
"""Exact checks for the foundational map's first arithmetic-dynamics note.

The all-iterate degree theorem is proved in the accompanying note by tracking
the unique top-total-degree summand.  This checker verifies the polynomial
inputs to that induction, the resulting matrix recurrence, the first two
literal iterates, and the exceptional two-line orbit.
"""

import sympy as sp


x, y, z = sp.symbols("x y z")
u = 1 + x*y
F = sp.Matrix(
    [
        u**3*z + y**2*u*(4 + 3*x*y),
        y + 3*x*u**2*z + 3*x*y**2*(4 + 3*x*y),
        2*x - 3*x**2*y - x**3*z,
    ]
)
variables = (x, y, z)


def total_degree(polynomial):
    return sp.Poly(sp.expand(polynomial), *variables).total_degree()


def top_homogeneous_part(polynomial):
    expanded = sp.Poly(sp.expand(polynomial), *variables)
    degree = expanded.total_degree()
    return sp.expand(
        sum(
            coefficient*x**exponents[0]*y**exponents[1]*z**exponents[2]
            for exponents, coefficient in expanded.terms()
            if sum(exponents) == degree
        )
    )


# The base multidegree and its unique top homogeneous summands.
assert [total_degree(component) for component in F] == [7, 6, 4]
assert [top_homogeneous_part(component) for component in F] == [
    x**3*y**3*z,
    3*x**3*y**2*z,
    -x**3*z,
]

# If (A,B,C) has component degrees (alpha,beta,gamma) and
# Delta=alpha-beta+gamma>0, the selected summands exceed their only
# competitors by Delta.  The new multidegree is M(alpha,beta,gamma).
alpha, beta, gamma = sp.symbols(
    "alpha beta gamma", integer=True, positive=True
)
delta = alpha - beta + gamma
selected_degrees = sp.Matrix(
    [
        3*alpha + 3*beta + gamma,
        3*alpha + 2*beta + gamma,
        3*alpha + gamma,
    ]
)
competitor_degrees = sp.Matrix(
    [
        2*alpha + 4*beta,
        2*alpha + 3*beta,
        2*alpha + beta,
    ]
)
assert sp.simplify(selected_degrees - competitor_degrees) == sp.Matrix(
    [delta, delta, delta]
)
assert sp.simplify(
    selected_degrees[0] - selected_degrees[1] + selected_degrees[2]
) == 3*alpha + beta + gamma

M = sp.Matrix([[3, 3, 1], [3, 2, 1], [3, 0, 1]])
assert selected_degrees == M*sp.Matrix([alpha, beta, gamma])
T = sp.symbols("T")
assert sp.factor(M.charpoly(T).as_expr()) == T*(T**2 - 6*T - 1)

# Literal composition checks the first nontrivial transition without relying
# on the abstract recurrence.
F2 = [
    sp.expand(component.subs(dict(zip(variables, F)), simultaneous=True))
    for component in F
]
assert [total_degree(component) for component in F2] == [43, 37, 25]

multidegrees = []
vector = sp.Matrix([1, 1, 1])
for _ in range(12):
    vector = M*vector
    multidegrees.append(tuple(int(entry) for entry in vector))

expected_first_six = [
    (7, 6, 4),
    (43, 37, 25),
    (265, 228, 154),
    (1633, 1405, 949),
    (10063, 8658, 5848),
    (62011, 53353, 36037),
]
assert multidegrees[:6] == expected_first_six
degrees = [1] + [entry[0] for entry in multidegrees]
assert all(
    degrees[index + 2] == 6*degrees[index + 1] + degrees[index]
    for index in range(len(degrees) - 2)
)

lambda_plus = 3 + sp.sqrt(10)
lambda_minus = 3 - sp.sqrt(10)
closed_form = (
    (4 + sp.sqrt(10))*lambda_plus**T
    + (sp.sqrt(10) - 4)*lambda_minus**T
) / (2*sp.sqrt(10))
assert sp.simplify(closed_form.subs(T, 0) - 1) == 0
assert sp.simplify(closed_form.subs(T, 1) - 7) == 0

# The Jacobian of every iterate is constant by the chain rule.
assert sp.factor(F.jacobian(variables).det()) == -2

# The coordinate lines L_x and L_z form an invariant two-cycle.  This gives
# exact escaping rational orbits with linear logarithmic-height growth.
q = sp.symbols("q")
Lx_image = tuple(
    sp.expand(component.subs({x: q, y: 0, z: 0})) for component in F
)
Lz_image = tuple(
    sp.expand(component.subs({x: 0, y: 0, z: q})) for component in F
)
assert Lx_image == (0, 0, 2*q)
assert Lz_image == (q, 0, 0)

# The first pullback of the nonproperness quartic is an exact small target for
# the next phase of the programme.  It has degree eight, far below the naive
# degree bound 4*deg(F)=28, but no invariant-factor claim is made.
a, b, c = F
discriminant = 27*a**2*c**2 - 18*a*b*c + 16*a + b**3*c - b**2
assert total_degree(discriminant) == 8
expected_pullback = (
    -9*x**4*y**2*z**2
    - 54*x**3*y**3*z
    - 18*x**3*y*z**2
    - 81*x**2*y**4
    - 72*x**2*y**2*z
    - 9*x**2*z**2
    - 54*x*y**3
    + 6*x*y*z
    + 63*y**2
    + 16*z
)
assert sp.factor(discriminant) == expected_pullback

print("PASS foundational dynamics: exact all-iterate degree recurrence")
print("PASS foundational dynamics: dynamical degree is 3 + sqrt(10)")
print("PASS foundational dynamics: generic degree tower is addressed in the note")
print("PASS foundational dynamics: exceptional two-line escaping orbit")
print("PASS foundational dynamics: first discriminant pullback has degree 8")

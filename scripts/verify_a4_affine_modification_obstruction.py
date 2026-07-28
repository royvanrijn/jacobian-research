#!/usr/bin/env python3
"""Exact audit of the forced A4 exceptional-quotient affine modification."""

from itertools import combinations

import sympy as sp


a, b, c, T = sp.symbols("a b c T")
q, w = sp.symbols("q w")
u, k, r, s = sp.symbols("u k r s")

rho = b**2 + 3 * b + 9
B = (
    a**3
    - 3 * a * b**2
    + 2 * b**3
    - 9 * a * b
    + 9 * b**2
    - 27 * a
    + 27 * b
    + 27
)
A = a**3 - b**3 - 9 * b**2 - 27 * b - 54
C = a**3 - b**3 + 27
P = (
    T**4
    - 6 * A * B * T**2
    - 8 * B**3 * T
    + B**2 * (9 * A**2 - 12 * C * B)
)

# Use z=rho and c=2b+3.  The etale local coefficient relation is
# c^2+27=4z, and c is a unit at either conjugate cluster point.
b_in_c = (c - 3) / 2
P_c = sp.expand(16 * P.subs(b, b_in_c))


def exact_divide_monomial(expression, variable, exponent):
    """Divide by a verified common monomial factor."""

    polynomial = sp.Poly(expression, variable)
    assert min(power[0] for power, coefficient in polynomial.terms()) >= exponent
    quotient = sp.cancel(expression / variable**exponent)
    assert sp.denom(quotient) == 1
    return sp.expand(quotient)


def jacobian_minors(equations, variables):
    """Return the maximal minors of the equation Jacobian."""

    jacobian = sp.Matrix(equations).jacobian(variables)
    size = len(equations)
    return [
        sp.expand(jacobian[:, columns].det())
        for columns in combinations(range(len(variables)), size)
    ]


def assert_zero_mod(expressions, generators, variables):
    """Verify ideal membership with one exact Groebner basis."""

    basis = sp.groebner(
        generators,
        *variables,
        order="lex",
        domain=sp.QQ,
    )
    for expression in expressions:
        assert basis.reduce(sp.expand(expression))[1] == 0


# ---------------------------------------------------------------------------
# 1. The direct affine contraction q=(T+a^3)/a^3 is nonnormal
# ---------------------------------------------------------------------------

direct_equation = sp.expand(P_c.subs(T, a**3 * (q - 1)))
cluster_relation = c**2 + 27

# At a=0 the equation has fourth order along the cluster and is independent
# of q.  The cluster cylinder is contained in the full hypersurface
# singular locus.
assert sp.factor(direct_equation.subs(a, 0)) == (
    3
    * c**2
    * cluster_relation**4
    * (11 * c**2 - 18 * c + 243)
    / 64
)
assert_zero_mod(
    [direct_equation]
    + [
        sp.diff(direct_equation, variable)
        for variable in (a, c, q)
    ],
    [a, cluster_relation],
    (q, a, c),
)

# This is a codimension-one singular cylinder in a hypersurface surface.
# Hence the direct affine-modification ring fails R1 and is not normal.
print("PASS: q=(T+a^3)/a^3 gives a codimension-one singular cylinder")
print("PASS: the direct exceptional-quotient affine modification is nonnormal")


# ---------------------------------------------------------------------------
# 2. The E3 integral chart is still singular
# ---------------------------------------------------------------------------

# Adjoin w=rho/a^3.  Reduction by c^2+27=4a^3w and division by the exact
# total order a^12 gives the E3 strict-transform equation.
e3_relation = c**2 + 27 - 4 * a**3 * w
e3_total = sp.rem(direct_equation, e3_relation, c)
e3_equation = exact_divide_monomial(e3_total, a, 12)

# The two conjugate points
#
#     a=0, q=1, 27w=c, c^2+27=0
#
# lie in the complete-intersection singular locus.
e3_point = [a, q - 1, 27 * w - c, cluster_relation]
assert_zero_mod(
    [e3_relation, e3_equation]
    + jacobian_minors(
        (e3_relation, e3_equation),
        (a, c, w, q),
    ),
    e3_point,
    (q, w, a, c),
)
print("PASS: adjoining w=rho/a^3 gives the exact E3 integral chart")
print("PASS: the E3 chart retains a conjugate pair of surface singularities")


# ---------------------------------------------------------------------------
# 3. The F chart is nonnormal along the forced triple line
# ---------------------------------------------------------------------------

# On the F chart,
#
#     a=u^2 k,  rho=u^3 k,  T=u^3 r-u^6 k^3.
#
# The last formula is G=T+a^3=u^3r.
f_relation = c**2 + 27 - 4 * u**3 * k
f_total = sp.rem(
    sp.expand(
        P_c.subs(
            {
                a: u**2 * k,
                T: u**3 * r - u**6 * k**3,
            }
        )
    ),
    f_relation,
    c,
)
f_equation = exact_divide_monomial(f_total, u, 12) / 8
assert sp.denom(f_equation) == 1
f_equation = sp.expand(f_equation)

f_exceptional = sp.expand(f_equation.subs(u, 0))
simple_line = 2 * r - (27 - 3 * c) * k
triple_line = 2 * r - (c - 9) * k
assert sp.rem(
    sp.expand(8 * f_exceptional - simple_line * triple_line**3),
    cluster_relation,
    c,
) == 0

# The total F surface is singular at the generic point of the triple line:
# all maximal Jacobian minors vanish modulo
# (u, triple_line, c^2+27).  This is again codimension one.
f_triple_prime = [u, triple_line, cluster_relation]
assert_zero_mod(
    [f_relation, f_equation]
    + jacobian_minors(
        (f_relation, f_equation),
        (u, k, c, r),
    ),
    f_triple_prime,
    (r, k, u, c),
)
print("PASS: the F exceptional equation is one simple line plus one triple line")
print("PASS: the F affine chart is nonnormal along the triple line")


# ---------------------------------------------------------------------------
# 4. Normalizing the triple line forces another singular center
# ---------------------------------------------------------------------------

# Its Newton edge is l^3+u^2*l+u^3, so the next integral quotient is
#
#     s=(2r-(c-9)k)/u.
#
# Substitute r=((c-9)k+us)/2 and divide the exact u^3 total order.
normalized_total = sp.rem(
    sp.together(
        f_equation.subs(
            r,
            ((c - 9) * k + u * s) / 2,
        )
    )
    * 16,
    f_relation,
    c,
)
normalized_equation = exact_divide_monomial(
    sp.expand(normalized_total),
    u,
    3,
)

# The new chart is still singular at u=k=s=0 over c^2+27=0.  This is the
# adjacent fan center which a single affine chart cannot contain smoothly.
normalized_point = [u, k, s, cluster_relation]
assert_zero_mod(
    [f_relation, normalized_equation]
    + jacobian_minors(
        (f_relation, normalized_equation),
        (u, k, c, s),
    ),
    normalized_point,
    (s, k, u, c),
)
print("PASS: normalization forces s=(2r-(c-9)k)/u")
print("PASS: the normalized triple chart reaches another singular fan center")


# ---------------------------------------------------------------------------
# 5. The only smooth completion is the nonaffine full resolution
# ---------------------------------------------------------------------------

# In chain order E1-F-E2-E3, the four point blowups have intersection
# matrix below.  It is negative definite and contracts back to the cluster,
# but the smooth resolution contains four complete exceptional curves and
# therefore is not affine.
intersection_matrix = sp.Matrix(
    [
        [-3, 1, 0, 0],
        [1, -1, 1, 0],
        [0, 1, -3, 1],
        [0, 0, 1, -1],
    ]
)
leading_minors = tuple(
    intersection_matrix[:size, :size].det()
    for size in range(1, 5)
)
assert leading_minors == (-3, 2, -3, 1)
assert intersection_matrix.det() == 1

print("PASS: the full exceptional chain is negative definite")
print("OBSTRUCTION: affine quotient charts are singular or nonnormal")
print("OBSTRUCTION: the smooth full-chain modification is not affine")
print("NOTE: any further affine model must delete a codimension-one divisor")

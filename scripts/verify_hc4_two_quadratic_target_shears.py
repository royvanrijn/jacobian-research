#!/usr/bin/env python3
"""Exclude length-two quadratic triangular target words over the reals.

For the normalized foundational map F=(P,B,C), take distinct i,j and let k
be the remaining target index.  Consider

    A = F_i + Q(F_j,F_k),
    L = F_j + R(A,F_k),

where Q and R have positive total degree at most two.  This checker verifies
an exact coefficient proof that the bordered-Hessian invariant

    K(L) = -grad(L)^T adj(Hess(L)) grad(L)

never vanishes identically for any of the six ordered words with real
coefficients.

The proof splits according to the genuine coupling in the second shear:

* r20 != 0, where R contains A^2;
* r20 = 0 and r11 != 0, where R contains A*F_k;
* r20 = r11 = 0, which is a single quadratic triangular shear after a
  target rescaling and is already covered by the cubic single-shear
  obstruction.

Only the first two strata require new calculations below.  Positivity of two
univariate factors is checked by exact real-root counts, not numerically.
"""

from __future__ import annotations

import sympy as sp


x, y, z = sp.symbols("x y z")
source_variables = (x, y, z)
t = 1 + x * y
q = t**2 * z + y**2 * (1 + 3 * t)
P = sp.expand(t * q)
B = sp.expand(y + 3 * x * q)
C = sp.expand(x * (5 - 3 * t) - x**3 * z)
mapping = (P, B, C)

(
    q10,
    q01,
    q20,
    q11,
    q02,
    r10,
    r01,
    r20,
    r11,
    r02,
) = sp.symbols("q10 q01 q20 q11 q02 r10 r01 r20 r11 r02")
parameters = (q10, q01, q20, q11, q02, r10, r01, r20, r11, r02)


def add_jets(
    left: tuple[sp.Expr, sp.Matrix, sp.Matrix],
    right: tuple[sp.Expr, sp.Matrix, sp.Matrix],
) -> tuple[sp.Expr, sp.Matrix, sp.Matrix]:
    return left[0] + right[0], left[1] + right[1], left[2] + right[2]


def source_jet(
    polynomial: sp.Expr,
    substitutions: dict[sp.Symbol, sp.Expr],
) -> tuple[sp.Expr, sp.Matrix, sp.Matrix]:
    return (
        sp.expand(polynomial.subs(substitutions)),
        sp.Matrix(
            [sp.diff(polynomial, variable).subs(substitutions) for variable in source_variables]
        ),
        sp.Matrix(
            3,
            3,
            lambda row, column: sp.diff(
                polynomial,
                source_variables[row],
                source_variables[column],
            ).subs(substitutions),
        ),
    )


def quadratic_jet(
    coefficients: tuple[sp.Expr, sp.Expr, sp.Expr, sp.Expr, sp.Expr],
    first: tuple[sp.Expr, sp.Matrix, sp.Matrix],
    second: tuple[sp.Expr, sp.Matrix, sp.Matrix],
) -> tuple[sp.Expr, sp.Matrix, sp.Matrix]:
    a10, a01, a20, a11, a02 = coefficients
    first_value, second_value = first[0], second[0]
    derivative_first = a10 + 2 * a20 * first_value + a11 * second_value
    derivative_second = a01 + a11 * first_value + 2 * a02 * second_value
    value = (
        a10 * first_value
        + a01 * second_value
        + a20 * first_value**2
        + a11 * first_value * second_value
        + a02 * second_value**2
    )
    gradient = derivative_first * first[1] + derivative_second * second[1]
    hessian = (
        derivative_first * first[2]
        + derivative_second * second[2]
        + 2 * a20 * first[1] * first[1].T
        + a11 * (first[1] * second[1].T + second[1] * first[1].T)
        + 2 * a02 * second[1] * second[1].T
    )
    return value, gradient, hessian


def invariant_from_jet(
    jet: tuple[sp.Expr, sp.Matrix, sp.Matrix],
    *,
    expand_result: bool = True,
) -> sp.Expr:
    gradient = jet[1]
    hessian = jet[2]
    h11, h12, h13 = hessian[0, 0], hessian[0, 1], hessian[0, 2]
    h22, h23, h33 = hessian[1, 1], hessian[1, 2], hessian[2, 2]
    gx, gy, gz = gradient
    result = -(
        gx**2 * (h22 * h33 - h23**2)
        + 2 * gx * gy * (h13 * h23 - h12 * h33)
        + 2 * gx * gz * (h12 * h23 - h13 * h22)
        + gy**2 * (h11 * h33 - h13**2)
        + 2 * gy * gz * (h12 * h13 - h11 * h23)
        + gz**2 * (h11 * h22 - h12**2)
    )
    return sp.expand(result) if expand_result else result


def word_invariant(
    indices: tuple[int, int, int],
    axis: sp.Symbol,
    parameter_substitutions: dict[sp.Symbol, sp.Expr] | None = None,
) -> sp.Poly:
    """Return K on one source axis, substituting parameters before expansion."""

    source_substitutions = {
        variable: sp.Integer(0)
        for variable in source_variables
        if variable != axis
    }
    jets = tuple(source_jet(polynomial, source_substitutions) for polynomial in mapping)
    substitutions = parameter_substitutions or {}
    coefficients = tuple(sp.sympify(entry).subs(substitutions) for entry in parameters)
    first_index, second_index, remaining_index = indices
    intermediate = add_jets(
        jets[first_index],
        quadratic_jet(
            coefficients[:5],  # type: ignore[arg-type]
            jets[second_index],
            jets[remaining_index],
        ),
    )
    retained = add_jets(
        jets[second_index],
        quadratic_jet(
            coefficients[5:],  # type: ignore[arg-type]
            intermediate,
            jets[remaining_index],
        ),
    )
    return sp.Poly(invariant_from_jet(retained), axis)


def direct_invariant(polynomial: sp.Expr) -> sp.Expr:
    gradient = sp.Matrix([sp.diff(polynomial, variable) for variable in source_variables])
    hessian = sp.hessian(polynomial, source_variables)
    return invariant_from_jet(
        (polynomial, gradient, hessian),
        expand_result=False,
    )


def restricted_polynomial(
    polynomial: sp.Expr,
    variable: sp.Symbol,
    substitutions: dict[sp.Symbol, sp.Expr],
) -> sp.Poly:
    return sp.Poly(sp.expand(direct_invariant(polynomial).subs(substitutions)), variable)


def coefficient(polynomial: sp.Poly, degree: int) -> sp.Expr:
    return sp.factor(polynomial.coeff_monomial(polynomial.gens[0] ** degree))


words = (
    (0, 1, 2),
    (0, 2, 1),
    (1, 0, 2),
    (1, 2, 0),
    (2, 0, 1),
    (2, 1, 0),
)


# ---------------------------------------------------------------------------
# Quadratic coupling r20 != 0.
# ---------------------------------------------------------------------------

x_axis = {word: word_invariant(word, x) for word in words}
assert coefficient(x_axis[(0, 1, 2)], 20) == 9_437_184 * q02**8 * r20**4
assert coefficient(x_axis[(0, 2, 1)], 20) == 9_437_184 * q20**8 * r20**4
assert coefficient(x_axis[(1, 0, 2)], 20) == 9_437_184 * q02**8 * r20**4
assert coefficient(x_axis[(1, 2, 0)], 20) == 9_437_184 * q20**8 * r20**4
assert coefficient(x_axis[(2, 0, 1)], 12) == 2_304 * r20**4
assert coefficient(x_axis[(2, 1, 0)], 12) == 2_304 * r20**4

# The two words retaining C close after q20=0 forces q10=0.
for word in ((0, 2, 1), (1, 2, 0)):
    branch = word_invariant(word, x, {q20: 0})
    assert coefficient(branch, 12) == 2_304 * q10**8 * r20**4
    assert coefficient(word_invariant(word, x, {q20: 0, q10: 0}), 8) == 9

# Word A=P+Q(B,C), L=B+R(A,C).
S = q01**2 * r20 + q01 * r11 + r02
branch_01 = word_invariant((0, 1, 2), x, {q02: 0})
assert coefficient(branch_01, 12) == 2_304 * S**4
positive_quartic = 48 * q20**4 + 48 * q20**3 + 836 * q20**2 - 268 * q20 + 99
assert sp.Poly(positive_quartic, q20).count_roots(-sp.oo, sp.oo) == 0
y_01 = word_invariant((0, 1, 2), y, {q02: 0, r02: -q01**2 * r20 - q01 * r11})
assert coefficient(y_01, 12) == (
    48 * r20**4 * (q20 + 4) ** 4 * positive_quartic
)
y_01_final = word_invariant(
    (0, 1, 2),
    y,
    {
        q02: 0,
        r02: -q01**2 * r20 - q01 * r11,
        q20: -4,
        q10: 0,
    },
)
assert coefficient(
    word_invariant(
        (0, 1, 2),
        y,
        {
            q02: 0,
            r02: -q01**2 * r20 - q01 * r11,
            q20: -4,
        },
    ),
    8,
) == 1_140_624 * q10**4 * r20**4
assert sp.expand(
    coefficient(y_01_final, 4)
    - 3 * r10 * (23_763 * r10**3 - 920 * r20)
) == 0
assert sp.expand(
    coefficient(y_01_final, 3) + 6 * (4_005 * r10**3 - 8 * r20)
) == 0
# For r20 != 0 the degree-three equation makes r10 nonzero; eliminating
# r20 between the last two displayed equations leaves a nonzero multiple.
assert sp.factor(
    (23_763 * r10**3 - 920 * r20).subs(
        r20,
        sp.Rational(4_005, 8) * r10**3,
    )
) != 0

# Word A=B+Q(P,C), L=P+R(A,C).
y_10 = word_invariant((1, 0, 2), y)
assert coefficient(y_10, 28) == 1_275_605_286_912 * q20**8 * r20**4
x_10_branch = word_invariant(
    (1, 0, 2),
    x,
    {q02: 0, r02: -q01**2 * r20 - q01 * r11, q20: 0},
)
assert coefficient(x_10_branch, 8) == 9 * (q01 * r10 + r01) ** 4
y_10_branch = word_invariant(
    (1, 0, 2),
    y,
    {
        q02: 0,
        r02: -q01**2 * r20 - q01 * r11,
        q20: 0,
        r01: -q01 * r10,
    },
)
assert coefficient(y_10_branch, 12) == 1_216_512 * q10**8 * r20**4
y_10_reduced = word_invariant(
    (1, 0, 2),
    y,
    {
        q02: 0,
        r02: -q01**2 * r20 - q01 * r11,
        q20: 0,
        r01: -q01 * r10,
        q10: 0,
    },
)
groebner_variables = (sp.Symbol("w"), q01, r11, r10, q11, r20)
inverse = groebner_variables[0]
groebner_basis = sp.groebner(
    [entry for (_, entry) in y_10_reduced.terms()]
    + [inverse * r20 - 1],
    *groebner_variables,
    order="grevlex",
)
assert groebner_basis.reduce(r10**3)[1] == 0
assert groebner_basis.reduce((r11 + 2 * q01 * r20) ** 2)[1] == 0

# Those forced equalities cancel q01 and reduce the word to
# L=P+r(B+u*P*C)^2.
u, r = sp.symbols("u r")
normal_quadratic = sp.expand(P + r * (B + u * P * C) ** 2)
normal_z = restricted_polynomial(normal_quadratic, z, {x: 0, y: 0})
normal_y = restricted_polynomial(normal_quadratic, y, {x: 0, z: 0})
assert coefficient(normal_z, 2) == -64 * r * u**2 - 168 * r * u - 108 * r + 9
normal_y_coefficient = coefficient(normal_y, 4)
forced_r = sp.Rational(9, 4) / ((8 * u + 9) * (2 * u + 3))
normal_line = restricted_polynomial(normal_quadratic, y, {x: 1, z: 0})
assert sp.factor(sp.cancel(coefficient(normal_line, 36).subs(r, forced_r))) == (
    -910_050_728_661
    * u**8
    / (4 * (2 * u + 3) ** 4 * (8 * u + 9) ** 4)
)
assert normal_y_coefficient.subs({u: 0, r: sp.Rational(1, 12)}) == sp.Rational(
    35_641,
    144,
)


# ---------------------------------------------------------------------------
# Bilinear coupling r20=0, r11 != 0.
# ---------------------------------------------------------------------------

r20_zero = {r20: 0}
x_bilinear = {word: word_invariant(word, x, r20_zero) for word in words}
assert coefficient(x_bilinear[(0, 1, 2)], 16) == 186_624 * q02**4 * r11**4
assert coefficient(x_bilinear[(1, 0, 2)], 16) == 186_624 * q02**4 * r11**4
assert coefficient(x_bilinear[(0, 2, 1)], 12) == 2_304 * q20**4 * r10**4
assert coefficient(x_bilinear[(1, 2, 0)], 12) == 2_304 * q20**4 * r10**4
assert coefficient(x_bilinear[(2, 0, 1)], 8) == 9 * r10**4
assert coefficient(x_bilinear[(2, 1, 0)], 8) == 9 * r10**4

# A=P+Q(B,C), L=B+R(A,C).
y_01_bilinear = word_invariant((0, 1, 2), y, {r20: 0, q02: 0})
assert coefficient(y_01_bilinear, 4) == 3 * r10**4 * positive_quartic
assert coefficient(
    word_invariant((0, 1, 2), y, {r20: 0, q02: 0, r10: 0}),
    0,
) == (2 * r11 + 3) ** 2
x_01_bilinear = word_invariant(
    (0, 1, 2),
    x,
    {r20: 0, q02: 0, r10: 0, r11: sp.Rational(-3, 2)},
)
assert coefficient(x_01_bilinear, 12) == 144 * (3 * q01 - 2 * r02) ** 4
a, b, c, d = sp.symbols("a b c d")
normal_01 = sp.expand(
    B
    + d * C
    - sp.Rational(3, 2) * P * C
    - sp.Rational(3, 2) * a * B * C
    - sp.Rational(3, 2) * b * B**2 * C
    - sp.Rational(3, 2) * c * B * C**2
)
line_01_y = restricted_polynomial(normal_01, y, {x: 1, z: 0})
line_01_x = restricted_polynomial(normal_01, x, {y: 1, z: 0})
assert coefficient(line_01_y, 24) == -sp.Rational(15_533_624_506_455, 16) * b**4
assert coefficient(line_01_x, 24).subs(b, 0) == -sp.Rational(
    3_486_784_401 * 19,
    16,
) * c**4
assert coefficient(
    restricted_polynomial(normal_01.subs({b: 0, c: 0}), y, {x: 1, z: 0}),
    16,
) == -sp.Rational(1_712_421, 4)

# A=B+Q(P,C), L=P+R(A,C).
assert coefficient(
    word_invariant((1, 0, 2), y, {r20: 0, q02: 0, r10: 0}),
    6,
) == 8_192 * q20**2 * r11**2
assert coefficient(
    word_invariant(
        (1, 0, 2),
        y,
        {r20: 0, q02: 0, r10: 0, q20: 0},
    ),
    4,
) == 297
assert coefficient(
    word_invariant(
        (1, 0, 2),
        y,
        {r20: 0, q02: 0, q20: 0},
    ),
    4,
) == 297 * (q10 * r10 + 1) ** 4
assert coefficient(
    word_invariant(
        (1, 0, 2),
        z,
        {r20: 0, q02: 0, q20: 0, q10: -1 / r10},
    ),
    0,
) == (2 * q11 * r10**2 + 3 * r10**2 - 2 * r11) ** 2
normal_10 = sp.expand(
    P
    + r
    * (B - P / r + q01 * C + u * P * C)
    + r01 * C
    + r**2
    * (2 * u + 3)
    / 2
    * (B - P / r + q01 * C + u * P * C)
    * C
    + r02 * C**2
)
line_10 = restricted_polynomial(normal_10, y, {x: 1, z: 0})
assert coefficient(line_10, 20) == (
    -sp.Rational(4_782_969, 4) * r**8 * u**4 * (2 * u + 3) ** 4
)
assert coefficient(
    restricted_polynomial(normal_10.subs(u, 0), y, {x: 1, z: 0}),
    16,
) == -sp.Rational(1_712_421, 4) * r**4

# The two words retaining C.  The branch r10=0 is immediate; the branch
# q20=0 first forces q10*r10=-1.
positive_quadratic = 9 * q02**2 + 27 * q02 + 53
assert sp.Poly(positive_quadratic, q02).count_roots(-sp.oo, sp.oo) == 0
assert coefficient(
    word_invariant((0, 2, 1), y, {r20: 0, r10: 0}),
    8,
) == 9 * r11**4 * positive_quadratic**2
assert coefficient(
    word_invariant((1, 2, 0), x, {r20: 0, r10: 0}),
    8,
) == 9
for word in ((0, 2, 1), (1, 2, 0)):
    assert coefficient(
        word_invariant(word, x, {r20: 0, q20: 0}),
        8,
    ) == 9 * (q10 * r10 + 1) ** 4

# A=P+Q(C,B), L=C+R(A,B), after q10*r10=-1.
assert coefficient(
    word_invariant(
        (0, 2, 1),
        y,
        {r20: 0, q20: 0, q10: -1 / r10},
    ),
    8,
) == 9 * r11**4 * positive_quadratic**2

# A=B+Q(C,P), L=C+R(A,P), after q10*r10=-1.
normal_12_intermediate = sp.expand(B - C / r + q01 * P + u * C * P + q02 * P**2)
normal_12 = sp.expand(
    C
    + r * normal_12_intermediate
    + r01 * P
    + r11 * normal_12_intermediate * P
    + r02 * P**2
)
assert coefficient(
    restricted_polynomial(normal_12, y, {x: 0, z: 0}),
    20,
) == 1_576_599_552 * r11**4 * q02**4
assert coefficient(
    restricted_polynomial(normal_12.subs(q02, 0), x, {y: 0, z: 0}),
    0,
) == (-2 * r11 + 2 * r**2 * u + 3 * r**2) ** 2
normal_12_forced = normal_12.subs(
    {
        q02: 0,
        r11: r**2 * (2 * u + 3) / 2,
    }
)
line_12 = restricted_polynomial(normal_12_forced, y, {x: 1, z: 0})
assert coefficient(line_12, 32) == (
    -sp.Rational(52_612_659, 16) * r**8 * u**4 * (2 * u + 3) ** 4
)
line_12_u0 = restricted_polynomial(normal_12_forced.subs(u, 0), y, {x: 1, z: 0})
assert coefficient(line_12_u0, 28) == -14_580 * (
    2 * r02 + 3 * q01 * r**2
) ** 4
assert coefficient(
    restricted_polynomial(
        normal_12_forced.subs(
            {
                u: 0,
                r02: -sp.Rational(3, 2) * q01 * r**2,
            }
        ),
        y,
        {x: 1, z: 0},
    ),
    24,
) == -sp.Rational(291_761_109, 4) * r**8

# A=C+Q(P,B), L=P+R(A,B).
assert coefficient(
    word_invariant((2, 0, 1), x, {r20: 0, r10: 0}),
    4,
) == 144 * r11**4

# A=C+Q(B,P), L=B+R(A,P).
assert coefficient(
    word_invariant((2, 1, 0), x, {r20: 0, r10: 0}),
    0,
) == (2 * r11 + 3) ** 2
assert coefficient(
    word_invariant(
        (2, 1, 0),
        y,
        {r20: 0, r10: 0, r11: sp.Rational(-3, 2)},
    ),
    20,
) == 7_981_535_232 * q02**4
assert coefficient(
    word_invariant(
        (2, 1, 0),
        y,
        {
            r20: 0,
            r10: 0,
            r11: sp.Rational(-3, 2),
            q02: 0,
        },
    ),
    16,
) == 7_776 * 21_304 * q11**4
normal_21 = sp.expand(
    B
    + d * P
    - sp.Rational(3, 2) * C * P
    - sp.Rational(3, 2) * a * B * P
    - sp.Rational(3, 2) * b * B**2 * P
    + c * P**2
)
assert coefficient(
    restricted_polynomial(normal_21, y, {x: 1, z: 0}),
    36,
) == -sp.Rational(8_755_315_630_911, 4) * b**4
assert coefficient(
    restricted_polynomial(normal_21.subs(b, 0), z, {x: 0, y: 0}),
    6,
) == 144 * c**4
normal_21_reduced = normal_21.subs({b: 0, c: 0})
assert coefficient(
    restricted_polynomial(normal_21_reduced, y, {x: 1, z: 0}),
    24,
) == -sp.Rational(291_761_109, 4) * a**4
assert coefficient(
    restricted_polynomial(normal_21_reduced.subs(a, 0), z, {x: 0, y: 0}),
    2,
) == 9 * d**4
assert coefficient(
    restricted_polynomial(
        normal_21_reduced.subs({a: 0, d: 0}),
        y,
        {x: 1, z: 0},
    ),
    16,
) == -sp.Rational(1_712_421, 4)


print("PASS: all six r20 != 0 quadratic-coupling words are excluded")
print("PASS: all six r20 = 0, r11 != 0 bilinear-coupling words are excluded")
print("PASS: the remaining r20 = r11 = 0 words reduce to one triangular shear")
print("PASS: no real length-two quadratic triangular target word has K(L)=0")

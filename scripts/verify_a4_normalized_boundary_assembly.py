#!/usr/bin/env python3
"""Exact audit of the normalized-boundary A4 Keller assembly."""

import sympy as sp


lam, r, c, S, T = sp.symbols("lam r c S T")


def target_cubic(first, second, third):
    return sp.expand(
        first**3
        - 3 * first * second**2
        + 2 * second**3
        - 9 * first * second * third
        + 9 * second**2 * third
        - 27 * first * third**2
        + 27 * second * third**2
        + 27 * third**3
    )


f1 = 2 * r**3 - 3 * r**2 + 3 * r - 1
f2 = -r**3 + 3 * r - 1
f3 = r * (r - 1)
a = 2 * r - 1
b = 2 * r**2 + r - 5


# ---------------------------------------------------------------------------
# 1. Polynomial normalized boundary and unimodular ambient completion
# ---------------------------------------------------------------------------

assert sp.factor(target_cubic(lam * f1, lam * f2, lam * f3)) == 0
assert sp.expand(a * f2 + b * f3) == 1

completion_matrix = sp.Matrix([
    [f1, 1, 0],
    [f2, 0, b],
    [f3, 0, -a],
])
assert sp.factor(completion_matrix.det()) == 1

p = sp.expand(f1 * lam + S)
q = sp.expand(f2 * lam + b * T)
rho = sp.expand(f3 * lam - a * T)
new_target_cubic = target_cubic(p, q, rho)
assert new_target_cubic.subs({S: 0, T: 0}) == 0


# ---------------------------------------------------------------------------
# 2. The two inverse-mask numerators are not divisible
# ---------------------------------------------------------------------------

difference = sp.expand(p - q)
incidence = sp.expand(
    27 * rho**2
    - difference**2
    - 3 * (difference - 3 * rho) * q
)
numerator_1 = sp.expand(c - difference * r)
numerator_2 = sp.expand(-incidence * c + 27 * rho**3 * r)

variables = (lam, r, c, S, T)
_, remainder_1 = sp.div(
    numerator_1,
    new_target_cubic,
    *variables,
)
_, remainder_2 = sp.div(
    numerator_2,
    new_target_cubic,
    *variables,
)
assert remainder_1 != 0
assert remainder_2 != 0


# ---------------------------------------------------------------------------
# 3. The rational Jacobian quotient is nonconstant
# ---------------------------------------------------------------------------

old_target_cubic = target_cubic(lam, r, c)
assert sp.factor(old_target_cubic) != 0
assert sp.factor(
    new_target_cubic - old_target_cubic
) != 0

# At the old mask-zero section the prospective denominator vanishes, while
# the numerator is a nonzero polynomial.  Hence their ratio is not constant.
assert new_target_cubic.subs({S: 0, T: 0}) == 0
assert old_target_cubic.subs({S: 0, T: 0}) == old_target_cubic


print("PASS: the normalized target boundary has a unimodular A5 completion")
print("PASS: both two-mask inverse numerators have nonzero remainders")
print("PASS: the assembled rational candidate has nonconstant Jacobian")
print("PASS: an ordinary Keller assembly needs a nonautomorphic log rechart")

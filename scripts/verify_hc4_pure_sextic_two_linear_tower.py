#!/usr/bin/env python3
"""Verify the complete two-linear-form closure on the pure-sixth boundary.

HC4RSD31 leaves

    c_6=x^6,
    c_5=x^4 L,
    c_4=x^3 M+x^2 L^2/4

up to unary terms and translation in x.  This checker treats the four
constant-rank possibilities for the passive linear forms L,M.  It verifies
the two exact cubic packets, the immutable descendants which kill their
curved members, and the all-degree cylinder identities used at the final
linear tails.
"""

from __future__ import annotations

import hashlib
import json
from functools import cache
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "hc4_pure_sextic_two_linear_tower.json"
)

x, y, z = sp.symbols("x y z")
variables = (x, y, z)


def bordered_invariant(polynomial: sp.Expr) -> sp.Expr:
    gradient = sp.Matrix([sp.diff(polynomial, variable) for variable in variables])
    hessian = sp.hessian(polynomial, variables)
    return sp.expand((gradient.T * hessian.adjugate() * gradient)[0])


@cache
def as_polynomial(polynomial: sp.Expr) -> sp.Poly:
    return sp.Poly(polynomial, *variables)


def coefficient(polynomial: sp.Expr, monomial: tuple[int, int, int]) -> sp.Expr:
    return sp.factor(as_polynomial(polynomial).coeff_monomial(monomial))


def homogeneous_face(polynomial: sp.Expr, degree: int) -> sp.Expr:
    result = 0
    for monomial, value in as_polynomial(polynomial).terms():
        if sum(monomial) == degree:
            result += value * x**monomial[0] * y**monomial[1] * z**monomial[2]
    return sp.expand(result)


def same(left: sp.Expr, right: sp.Expr) -> bool:
    return sp.expand(left - right) == 0


b = sp.symbols("b0:10")
q = sp.symbols("q0:6")
ell = sp.symbols("ell0:3")
cubic_monomials = (
    x**3,
    x**2 * y,
    x**2 * z,
    x * y**2,
    x * y * z,
    x * z**2,
    y**3,
    y**2 * z,
    y * z**2,
    z**3,
)
quadratic_monomials = (x**2, x * y, x * z, y**2, y * z, z**2)
tail = (
    sum(value * monomial for value, monomial in zip(b, cubic_monomials))
    + sum(value * monomial for value, monomial in zip(q, quadratic_monomials))
    + ell[0] * x
    + ell[1] * y
    + ell[2] * z
)


# There are two cubic coefficient packets.  Their reduced radical bases
# give elementary two-chart covers: b_3 != 0 is the finite chart and
# b_3=0 is completed by the b_5 chart at infinity.
shifted_degree_twelve = (
    36 * b[3] * b[5] - 9 * b[4] ** 2 + 6 * b[4] - 1,
    36 * b[3] * b[8]
    - 36 * b[4] * b[7]
    + 108 * b[5] * b[6]
    - b[5]
    + 12 * b[7],
    9 * b[3] * b[9] - 3 * b[4] * b[8] + 3 * b[5] * b[7] + b[8],
    108 * b[6] * b[8] - 36 * b[7] ** 2 - b[8],
    108 * b[6] * b[9] - 12 * b[7] * b[8] - b[9],
    3 * b[7] * b[9] - b[8] ** 2,
)
shifted_basis = (
    -b[8] ** 2 + 3 * b[7] * b[9],
    -12 * b[7] * b[8] + 108 * b[6] * b[9] - b[9],
    -36 * b[7] ** 2 + 108 * b[6] * b[8] - b[8],
    -2 * b[5] * b[8] + 3 * b[4] * b[9] - b[9],
    -6 * b[5] * b[7] + 3 * b[4] * b[8] - b[8],
    -108 * b[5] * b[6] + 18 * b[4] * b[7] + b[5] - 6 * b[7],
    -b[5] * b[7] + 3 * b[3] * b[9],
    -108 * b[5] * b[6] + 36 * b[3] * b[8] + b[5],
    -324 * b[4] * b[6]
    + 216 * b[3] * b[7]
    + 3 * b[4]
    + 108 * b[6]
    - 1,
    -9 * b[4] ** 2 + 36 * b[3] * b[5] + 6 * b[4] - 1,
)
shifted_original_gb = sp.groebner(shifted_degree_twelve, *b[3:10], order="grevlex")
shifted_expected_gb = sp.groebner(shifted_basis, *b[3:10], order="grevlex")
assert all(shifted_original_gb.reduce(value**2)[1] == 0 for value in shifted_basis)
assert all(
    shifted_expected_gb.reduce(value)[1] == 0 for value in shifted_degree_twelve
)

unshifted_degree_twelve = (
    4 * b[3] * b[5] - b[4] ** 2,
    b[3] * b[8] - b[4] * b[7] + 3 * b[5] * b[6],
    3 * b[3] * b[9] - b[4] * b[8] + b[5] * b[7],
    3 * b[6] * b[8] - b[7] ** 2,
    9 * b[6] * b[9] - b[7] * b[8],
    3 * b[7] * b[9] - b[8] ** 2,
)
unshifted_basis = (
    -b[8] ** 2 + 3 * b[7] * b[9],
    -b[7] * b[8] + 9 * b[6] * b[9],
    -b[7] ** 2 + 3 * b[6] * b[8],
    -2 * b[5] * b[8] + 3 * b[4] * b[9],
    -2 * b[5] * b[7] + b[4] * b[8],
    -6 * b[5] * b[6] + b[4] * b[7],
    -b[5] * b[7] + 3 * b[3] * b[9],
    -3 * b[5] * b[6] + b[3] * b[8],
    -3 * b[4] * b[6] + 2 * b[3] * b[7],
    -b[4] ** 2 + 4 * b[3] * b[5],
)
unshifted_original_gb = sp.groebner(
    unshifted_degree_twelve, *b[3:10], order="grevlex"
)
unshifted_expected_gb = sp.groebner(unshifted_basis, *b[3:10], order="grevlex")
assert all(unshifted_original_gb.reduce(value**2)[1] == 0 for value in unshifted_basis)
assert all(
    unshifted_expected_gb.reduce(value)[1] == 0
    for value in unshifted_degree_twelve
)

a, p, u = sp.symbols("a p u")
finite_shifted = {
    b[3]: p,
    b[4]: sp.Rational(1, 3) + 2 * a * p,
    b[5]: a**2 * p,
    b[6]: (u + 1) / 108,
    b[7]: a * u / 36,
    b[8]: a**2 * u / 36,
    b[9]: a**3 * u / 108,
}
h, r, v = sp.symbols("h r v")
infinite_shifted = {
    b[3]: h**2 * r,
    b[4]: sp.Rational(1, 3) + 2 * h * r,
    b[5]: r,
    b[6]: (1 + h**3 * v) / 108,
    b[7]: h**2 * v / 36,
    b[8]: h * v / 36,
    b[9]: v / 108,
}
finite_unshifted = {
    b[3]: p,
    b[4]: 2 * a * p,
    b[5]: a**2 * p,
    b[6]: u,
    b[7]: 3 * a * u,
    b[8]: 3 * a**2 * u,
    b[9]: a**3 * u,
}
infinite_unshifted = {
    b[3]: h**2 * r,
    b[4]: 2 * h * r,
    b[5]: r,
    b[6]: h**3 * v,
    b[7]: 3 * h**2 * v,
    b[8]: 3 * h * v,
    b[9]: v,
}
for packet, substitutions in (
    (shifted_degree_twelve, finite_shifted),
    (shifted_degree_twelve, infinite_shifted),
    (unshifted_degree_twelve, finite_unshifted),
    (unshifted_degree_twelve, infinite_unshifted),
):
    assert all(sp.expand(value.subs(substitutions)) == 0 for value in packet)


# Independent L,M.  Normalize L=y and M=z.  The finite cubic packet first
# aligns its moving direction, then degree nine kills both curvature
# parameters.  The infinity chart has direct v^2 and r^2 obstructions.
independent = bordered_invariant(
    x**6 + x**4 * y + x**2 * y**2 / 4 + x**3 * z + tail
)
assert sp.expand(homogeneous_face(independent, 12).subs(finite_shifted)) == 0
E = -2 * a**2 * b[1] + 12 * a**2 * q[3] + 2 * a * b[2] - 12 * a * q[4] + 12 * q[5] - 1
assert same(coefficient(independent.subs(finite_shifted), (11, 0, 0)), 12 * p * E)
assert same(coefficient(independent.subs(finite_shifted), (10, 1, 0)), u * E / 3)
q5_from_E = (
    1 + 2 * a**2 * b[1] - 12 * a**2 * q[3] - 2 * a * b[2] + 12 * a * q[4]
) / 12
finite_collision = sp.expand(independent.subs(finite_shifted).subs(q[5], q5_from_E))
R = -2 * a * b[1] + 12 * a * q[3] + b[2] - 6 * q[4]
assert same(coefficient(finite_collision, (10, 0, 0)), -R**2)
assert same(coefficient(finite_collision, (6, 4, 0)), -a**2 * u**2 / 1296)
assert sp.expand(
    coefficient(finite_collision, (8, 2, 0))
    - (a * u * R / 18 - 4 * a**2 * p**2)
) == 0
aligned_collision = finite_collision.subs({a: 0, q[4]: b[2] / 6})
assert same(coefficient(aligned_collision, (8, 0, 1)), -4 * p)
assert same(coefficient(aligned_collision, (7, 1, 1)), -u / 9)

infinite_independent = sp.expand(independent.subs(infinite_shifted))
assert sp.expand(homogeneous_face(infinite_independent, 12)) == 0
assert same(coefficient(infinite_independent, (6, 0, 4)), -v**2 / 1296)
assert same(coefficient(infinite_independent.subs(v, 0), (8, 0, 2)), -4 * r**2)

base_independent = independent.subs(
    {
        b[3]: 0,
        b[4]: sp.Rational(1, 3),
        b[5]: 0,
        b[6]: sp.Rational(1, 108),
        b[7]: 0,
        b[8]: 0,
        b[9]: 0,
    }
)
assert sp.expand(
    coefficient(base_independent, (4, 4, 0)) - 5 * (12 * q[5] - 1) / 1296
) == 0
assert sp.expand(
    coefficient(base_independent, (0, 6, 0)) - (9 * q[5] - 1) / 11664
) == 0
assert same(
    coefficient(base_independent.subs(q[5], sp.Rational(1, 12)), (0, 6, 0)),
    -sp.Rational(1, 46656),
)


# Dependent L != 0.  Translation in x removes the component of M along L,
# so normalize L=y and M=0.  Five square coefficients make the cubic tail
# binary.  The nonzero z^2 branch dies by -4*q5^2; in the zero branch two
# incompatible resonances force the remaining x^2*z coefficient to vanish.
dependent = bordered_invariant(x**6 + x**4 * y + x**2 * y**2 / 4 + tail)
assert same(coefficient(dependent, (6, 0, 4)), -9 * b[9] ** 2)
assert same(coefficient(dependent.subs(b[9], 0), (6, 2, 2)), -112 * b[8] ** 2)
assert same(
    coefficient(dependent.subs({b[8]: 0, b[9]: 0}), (6, 4, 0)),
    -109 * b[7] ** 2,
)
assert same(
    coefficient(dependent.subs({b[7]: 0, b[8]: 0, b[9]: 0}), (8, 0, 2)),
    -4 * b[5] ** 2,
)
assert same(
    coefficient(
        dependent.subs({b[5]: 0, b[7]: 0, b[8]: 0, b[9]: 0}),
        (12, 0, 0),
    ),
    -36 * b[4] ** 2,
)
dependent_reduced = dependent.subs(
    {b[4]: 0, b[5]: 0, b[7]: 0, b[8]: 0, b[9]: 0}
)
assert same(coefficient(dependent_reduced, (11, 0, 0)), 144 * b[3] * q[5])
assert same(
    coefficient(dependent_reduced, (10, 1, 0)),
    4 * q[5] * (108 * b[6] - 1),
)
dependent_nonzero_scalar = dependent_reduced.subs(
    {b[3]: 0, b[6]: sp.Rational(1, 108)}
)
assert same(coefficient(dependent_nonzero_scalar, (6, 0, 2)), -4 * q[5] ** 2)
dependent_zero_scalar = dependent_reduced.subs({q[5]: 0, q[4]: b[2] / 6})
assert same(
    coefficient(dependent_zero_scalar, (4, 1, 1)),
    -b[2] ** 3 * (72 * b[6] - 1) / 2,
)
assert same(
    coefficient(dependent_zero_scalar, (2, 2, 1)),
    -b[2] ** 3 * (108 * b[6] - 1) / 18,
)


# The last transverse form is h(x,y)+z*(alpha*x+beta).  Its z coefficient
# is a product in the polynomial ring and forces alpha=0 because h_yy has
# the nonzero leading term x^2/2.
alpha, beta = sp.symbols("alpha beta")
hx, hy, hxx, hxy, hyy = sp.symbols("hx hy hxx hxy hyy")
generic_gradient = sp.Matrix([hx + alpha * z, hy, alpha * x + beta])
generic_hessian = sp.Matrix(
    [[hxx, hxy, alpha], [hxy, hyy, 0], [alpha, 0, 0]]
)
generic_linear_tail = sp.expand(
    (generic_gradient.T * generic_hessian.adjugate() * generic_gradient)[0]
)
assert same(
    sp.Poly(generic_linear_tail, z).coeff_monomial(z),
    -2 * alpha**2 * (alpha * x + beta) * hyy,
)


# L=0, M!=0.  Normalize M=z.  The unshifted cubic packet has the same two
# charts.  Its finite curved chart has incompatible 1/12 and 3/20 faces.
transverse_quartic = bordered_invariant(x**6 + x**3 * z + tail)
assert sp.expand(
    homogeneous_face(transverse_quartic, 12).subs(finite_unshifted)
) == 0
S = a**2 * q[3] - a * q[4] + q[5]
assert same(
    coefficient(transverse_quartic.subs(finite_unshifted), (11, 0, 0)),
    12 * p * (12 * S - 1),
)
finite_transverse_collision = transverse_quartic.subs(
    finite_unshifted
).subs(q[5], sp.Rational(1, 12) - a**2 * q[3] + a * q[4])
assert same(coefficient(finite_transverse_collision, (4, 4, 0)), -36 * u**2)
assert same(
    coefficient(finite_transverse_collision.subs(u, 0), (6, 2, 0)),
    -12 * p**2,
)

S_infinity = h**2 * q[5] - h * q[4] + q[3]
infinite_transverse = transverse_quartic.subs(infinite_unshifted)
assert same(
    coefficient(infinite_transverse, (10, 0, 1)),
    36 * v * (12 * S_infinity - h**2),
)
infinite_transverse_collision = infinite_transverse.subs(
    q[3], h**2 / 12 - h**2 * q[5] + h * q[4]
)
assert same(
    coefficient(infinite_transverse_collision, (4, 0, 4)), -36 * h**2 * v**2
)
assert same(
    coefficient(infinite_transverse_collision.subs(v, 0), (6, 0, 2)),
    -12 * h**2 * r**2,
)

infinite_endpoint = transverse_quartic.subs(
    {
        b[3]: 0,
        b[4]: 0,
        b[5]: r,
        b[6]: 0,
        b[7]: 0,
        b[8]: 0,
        b[9]: v,
        q[3]: 0,
        q[4]: 0,
    }
)
assert same(coefficient(infinite_endpoint, (9, 0, 0)), 12 * b[1] ** 2 * r)
assert same(coefficient(infinite_endpoint, (8, 0, 1)), 36 * b[1] ** 2 * v)
endpoint_b1_zero = infinite_endpoint.subs(b[1], 0)
assert same(coefficient(endpoint_b1_zero, (7, 0, 0)), 36 * q[1] ** 2 * r)
assert same(coefficient(endpoint_b1_zero, (0, 0, 4)), -9 * q[1] ** 2 * v**2)
endpoint_q1_zero = endpoint_b1_zero.subs(q[1], 0)
assert same(coefficient(endpoint_q1_zero, (5, 0, 0)), 60 * ell[1] ** 2 * r)
assert same(coefficient(endpoint_q1_zero, (4, 0, 1)), 180 * ell[1] ** 2 * v)

base_transverse = transverse_quartic.subs({value: 0 for value in b[3:10]})
E0 = 12 * q[3] * q[5] - q[3] - 3 * q[4] ** 2
E1 = 20 * q[3] * q[5] - 2 * q[3] - 5 * q[4] ** 2
assert same(coefficient(base_transverse, (10, 0, 0)), 12 * E0)
assert same(coefficient(base_transverse, (7, 0, 1)), 12 * E1)
assert same(E1 - sp.Rational(5, 3) * E0, -q[3] / 3)
base_passive_linear = base_transverse.subs({q[3]: 0, q[4]: 0})
assert same(
    coefficient(base_passive_linear, (8, 0, 0)), b[1] ** 2 * (12 * q[5] - 1)
)
assert same(coefficient(base_passive_linear, (5, 0, 1)), -4 * b[1] ** 2 * q[5])
base_b1_zero = base_passive_linear.subs(b[1], 0)
assert same(
    coefficient(base_b1_zero, (6, 0, 0)), 4 * q[1] ** 2 * (9 * q[5] - 1)
)
assert same(coefficient(base_b1_zero, (0, 0, 2)), -4 * q[1] ** 2 * q[5] ** 2)


# L=M=0.  A nonzero cubic packet and the passive quadratic align to one
# constant passive form.  The remaining possible direction break is
# h(x,w)+v*D(x), whose v coefficient is exact.  If the passive quadratic
# also vanishes, the whole invariant is a Wronskian square.
assert sp.expand(
    4 * q[3] * S - (2 * a * q[3] - q[4]) ** 2
    - (4 * q[3] * q[5] - q[4] ** 2)
) == 0

D, Dp, Dpp = sp.symbols("D Dp Dpp")
hw, hxw = sp.symbols("hw hxw")
generic_D_gradient = sp.Matrix([hx + z * Dp, hw, D])
generic_D_hessian = sp.Matrix(
    [[hxx + z * Dpp, hxw, Dp], [hxw, hyy, 0], [Dp, 0, 0]]
)
generic_D_invariant = sp.expand(
    (generic_D_gradient.T * generic_D_hessian.adjugate() * generic_D_gradient)[0]
)
assert same(
    sp.Poly(generic_D_invariant, z).coeff_monomial(z),
    D * hyy * (D * Dpp - 2 * Dp**2),
)

gamma, delta, epsilon = sp.symbols("gamma delta epsilon")
quadratic_D = gamma * x**2 + delta * x + epsilon
quadratic_resonance = sp.Poly(
    sp.expand(quadratic_D * sp.diff(quadratic_D, x, 2) - 2 * sp.diff(quadratic_D, x) ** 2),
    x,
)
assert quadratic_resonance.coeff_monomial(x**2) == -6 * gamma**2
assert quadratic_resonance.coeff_monomial(1).subs(gamma, 0) == -2 * delta**2

f, fp, g, gp = sp.symbols("f fp g gp")
wronskian_gradient = sp.Matrix([hx + y * fp + z * gp, f, g])
wronskian_hessian = sp.Matrix(
    [[hxx, fp, gp], [fp, 0, 0], [gp, 0, 0]]
)
wronskian_invariant = sp.factor(
    (wronskian_gradient.T * wronskian_hessian.adjugate() * wronskian_gradient)[0]
)
assert wronskian_invariant == -(f * gp - fp * g) ** 2


payload = {
    "format": "hc4-pure-sextic-two-linear-tower-v1",
    "status": {
        "id": "HC4RSD32",
        "kind": "hybrid theorem",
        "scope": "the complete two-linear-form tower left by HC4RSD31",
        "result": (
            "the independent packet is empty and every dependent packet is "
            "a fixed cylinder; hence every scalar pure-sixth boundary closes"
        ),
    },
    "cubic_packets": {
        "independent_shifted": "prime dimension-three two-chart packet",
        "dependent_unshifted": "prime dimension-three two-chart packet",
    },
    "independent_terminal_coefficients": [
        "-a^2*u^2/1296",
        "-4*a^2*p^2",
        "-4*p",
        "-u/9",
        "-v^2/1296",
        "-4*r^2",
        "-1/46656",
    ],
    "global_obstructions": [
        "D*h_ww*(D*D''-2*(D')^2)",
        "-(f*g'-f'*g)^2",
    ],
    "conclusion": "all scalar degree-six leading directions are fixed cylinders",
}

serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
OUTPUT.write_text(serialized)
digest = hashlib.sha256(serialized.encode()).hexdigest()

print("PASS: classified both cubic two-chart packets")
print("PASS: closed the independent two-linear-form tower")
print("PASS: closed all three dependent-rank endpoints")
print("PASS: verified the global D-resonance and Wronskian-square obstructions")
print("THEOREM: every scalar pure-sixth boundary is a fixed cylinder")
print(f"PASS: wrote {OUTPUT.relative_to(ROOT)}")
print(f"SHA256: {digest}")

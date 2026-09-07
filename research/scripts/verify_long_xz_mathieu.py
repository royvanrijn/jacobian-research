#!/usr/bin/env python3
"""Exact, dependency-free checks for Long's (xz), SU(2), and SO(3) witnesses.

Source: Christopher D. Long, "Counterexamples to the (xz)-Conjecture and
the Mathieu Conjecture for SU(2)", arXiv:2607.19012v1 (21 July 2026).

The (xz) calculation is self-contained.  This script verifies the SU(2)
algebraic substitution and monomial identities on the right side of the
Mueger--Tuset integration formula quoted by Long.  The companion
verify_long_su2_haar.py supplies the independent Haar-measure proof.

The SO(3) calculation replays Long's 28 July 2026 announcement.  It uses
only that the third column of a Haar rotation is uniform on S^2, where
U*V+T^2=1, the phase of U is uniform, and T is uniform on [-1,1].
"""

from fractions import Fraction
from math import comb, factorial


Laurent = dict[tuple[int, int, int], Fraction]  # powers of x, z1, z2
SphereLaurent = dict[tuple[int, int], Fraction]  # powers of U, T


def add(*polys: Laurent) -> Laurent:
    out: Laurent = {}
    for poly in polys:
        for monomial, coefficient in poly.items():
            out[monomial] = out.get(monomial, Fraction(0)) + coefficient
            if not out[monomial]:
                del out[monomial]
    return out


def multiply(left: Laurent, right: Laurent) -> Laurent:
    out: Laurent = {}
    for (ax, az1, az2), ac in left.items():
        for (bx, bz1, bz2), bc in right.items():
            monomial = (ax + bx, az1 + bz1, az2 + bz2)
            out[monomial] = out.get(monomial, Fraction(0)) + ac * bc
    return {monomial: coefficient for monomial, coefficient in out.items() if coefficient}


def power(poly: Laurent, exponent: int) -> Laurent:
    out: Laurent = {(0, 0, 0): Fraction(1)}
    base = poly
    n = exponent
    while n:
        if n & 1:
            out = multiply(out, base)
        base = multiply(base, base)
        n //= 2
    return out


def integral_constant_term(poly: Laurent) -> Fraction:
    """Apply integral_0^1 CT_(z1,z2), term by term."""
    total = Fraction(0)
    for (x_degree, z1_degree, z2_degree), coefficient in poly.items():
        if z1_degree == z2_degree == 0:
            total += coefficient * Fraction(1, x_degree + 1)
    return total


ONE: Laurent = {(0, 0, 0): Fraction(1)}
X: Laurent = {(1, 0, 0): Fraction(1)}
Z1: Laurent = {(0, 1, 0): Fraction(1)}
Z1_INV: Laurent = {(0, -1, 0): Fraction(1)}
Z2: Laurent = {(0, 0, 1): Fraction(1)}
Z2_INV: Laurent = {(0, 0, -1): Fraction(1)}

# f=(1-z^-1)((1-x)+xz)
f = multiply(
    add(ONE, {(0, -1, 0): Fraction(-1)}),
    add(ONE, {(1, 0, 0): Fraction(-1)}, multiply(X, Z1)),
)


def beta_integral(n: int, k: int) -> Fraction:
    """Integral of x^k(1-x)^(n-k), expanded exactly."""
    return sum(
        Fraction((-1) ** j * comb(n - k, j), k + j + 1)
        for j in range(n - k + 1)
    )


def mueger_tuset_monomial_rhs(r: int, s: int, t: int, u: int) -> Fraction:
    """Algebraic beta-map torus/Beta integral for a^r b^s c^t d^u."""
    if r != u or s != t:
        return Fraction(0)
    return Fraction((-1) ** t * factorial(r) * factorial(s), factorial(r + s + 1))


def sphere_add(*polys: SphereLaurent) -> SphereLaurent:
    out: SphereLaurent = {}
    for poly in polys:
        for monomial, coefficient in poly.items():
            out[monomial] = out.get(monomial, Fraction(0)) + coefficient
            if not out[monomial]:
                del out[monomial]
    return out


def sphere_scale(coefficient: int, poly: SphereLaurent) -> SphereLaurent:
    return {
        monomial: Fraction(coefficient) * value
        for monomial, value in poly.items()
        if coefficient * value
    }


def sphere_multiply(
    left: SphereLaurent, right: SphereLaurent
) -> SphereLaurent:
    out: SphereLaurent = {}
    for (left_u, left_t), left_coefficient in left.items():
        for (right_u, right_t), right_coefficient in right.items():
            monomial = (left_u + right_u, left_t + right_t)
            out[monomial] = (
                out.get(monomial, Fraction(0))
                + left_coefficient * right_coefficient
            )
    return {
        monomial: coefficient
        for monomial, coefficient in out.items()
        if coefficient
    }


def sphere_power(poly: SphereLaurent, exponent: int) -> SphereLaurent:
    out: SphereLaurent = {(0, 0): Fraction(1)}
    base = poly
    n = exponent
    while n:
        if n & 1:
            out = sphere_multiply(out, base)
        base = sphere_multiply(base, base)
        n //= 2
    return out


def sphere_integral(poly: SphereLaurent) -> Fraction:
    """Apply uniform phase extraction and normalized height integration."""
    total = Fraction(0)
    for (u_degree, t_degree), coefficient in poly.items():
        if u_degree == 0 and t_degree % 2 == 0:
            total += coefficient * Fraction(1, t_degree + 1)
    return total


def long_so3_witness() -> tuple[SphereLaurent, SphereLaurent]:
    """Return P,Q after eliminating V through U*V+T^2=1."""
    one: SphereLaurent = {(0, 0): Fraction(1)}
    u: SphereLaurent = {(1, 0): Fraction(1)}
    u_inverse: SphereLaurent = {(-1, 0): Fraction(1)}
    t_squared: SphereLaurent = {(0, 2): Fraction(1)}

    # V=(1-T^2)/U on the sphere.
    v = sphere_add(
        u_inverse,
        sphere_scale(-1, sphere_multiply(u_inverse, t_squared)),
    )
    displayed = sphere_multiply(
        sphere_add(one, u),
        sphere_add(
            v,
            sphere_scale(
                -1,
                sphere_multiply(sphere_add(sphere_scale(2, one), u), t_squared),
            ),
        ),
    )

    # Equivalent endpoint-contact form:
    # P=(1+U)/U * (1-T^2(1+U)^2).
    endpoint_form = sphere_multiply(
        sphere_multiply(sphere_add(one, u), u_inverse),
        sphere_add(
            one,
            sphere_scale(
                -1,
                sphere_multiply(t_squared, sphere_power(sphere_add(one, u), 2)),
            ),
        ),
    )
    assert displayed == endpoint_form
    return displayed, u


def long_so3_height_integral(order: int) -> Fraction:
    """Integral_0^1 (1-v^2)^order dv, evaluated termwise."""
    return sum(
        Fraction((-1) ** index * comb(order, index), 2 * index + 1)
        for index in range(order + 1)
    )


def long_so3_jet(order: int, degree: int) -> Fraction:
    """Coefficient of s^degree in H_m(1+s)."""
    return sum(
        Fraction(
            (-1) ** index
            * comb(order, index)
            * comb(order + 2 * index, degree),
            2 * index + 1,
        )
        for index in range(order + 1)
    )


def main() -> None:
    # Proof-oriented beta/binomial identity.
    for n in range(1, 21):
        for k in range(n + 1):
            assert comb(n, k) * beta_integral(n, k) == Fraction(1, n + 1)

    # Direct Laurent-polynomial regression of the two moments.
    f_power = ONE
    for n in range(1, 16):
        f_power = multiply(f_power, f)
        assert integral_constant_term(f_power) == 0
        assert integral_constant_term(multiply(Z1_INV, f_power)) == Fraction(
            (-1) ** (n - 1), n + 1
        )

    # Mueger--Tuset beta substitution, in Long's coordinate order (a,b,c,d):
    # ((1-x)z2, xz1, -z1^-1, z2^-1).
    a = add(Z2, multiply({(0, 0, 1): Fraction(-1)}, X))
    b = multiply(X, Z1)
    c = {(0, -1, 0): Fraction(-1)}
    d = Z2_INV
    F_beta = multiply(add(ONE, c), add(multiply(a, d), b))
    G_beta = {(0, -1, 0): Fraction(1)}  # -c
    assert F_beta == f
    assert G_beta == Z1_INV

    # Check the quoted formula's monomial right side in a useful exact box.
    for r in range(5):
        for s in range(5):
            for t in range(5):
                for u in range(5):
                    beta_image = multiply(
                        multiply(power(a, r), power(b, s)),
                        multiply(power(c, t), power(d, u)),
                    )
                    assert integral_constant_term(beta_image) == mueger_tuset_monomial_rhs(
                        r, s, t, u
                    )

    # Long's SO(3) witness depends only on the third rotation column.
    # Haar pushforward makes that column uniform on S^2.  Direct phase/height
    # integration checks the displayed moments without using SU(2).
    so3_p, so3_q = long_so3_witness()
    so3_power: SphereLaurent = {(0, 0): Fraction(1)}
    for order in range(1, 16):
        so3_power = sphere_multiply(so3_power, so3_p)
        expected = Fraction(
            4**order * factorial(order) ** 2,
            factorial(2 * order + 1),
        )
        assert sphere_integral(so3_power) == 0
        assert sphere_integral(sphere_multiply(so3_q, so3_power)) == expected

    # Replay the all-order endpoint-jet mechanism through a much larger
    # exact cutoff.  The note proves it uniformly: after w=vX,
    # H_m(X)=X^(m-1) J_m(X), and J_m'(X)=(1-X^2)^m has an m-fold zero at 1.
    for order in range(1, 101):
        height_integral = long_so3_height_integral(order)
        expected = Fraction(
            4**order * factorial(order) ** 2,
            factorial(2 * order + 1),
        )
        assert height_integral == expected
        assert long_so3_jet(order, order) == 0
        assert long_so3_jet(order, order - 1) == height_integral

    print("PASS Long xz: beta/binomial identity n=1..20")
    print("PASS Long xz: exact Laurent moments n=1..15")
    print("PASS Long SU(2): beta substitution and monomial RHS in degrees 0..4")
    print("PASS Long SU(2): combine with verify_long_su2_haar.py for the full Haar proof")
    print("PASS Long SO(3): exact spherical moments n=1..15")
    print("PASS Long SO(3): endpoint-jet and beta identities n=1..100")


if __name__ == "__main__":
    main()

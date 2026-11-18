"""Replay integral moment congruences and the full moving cubic source."""

from __future__ import annotations

import argparse
from collections import defaultdict
import hashlib
import json
from math import factorial
from pathlib import Path

import sympy as sp


if not __debug__:
    raise RuntimeError("Exact assertions must be enabled")


def truncated_product(left, right, prime):
    result = defaultdict(int)
    for (s, y), a in left.items():
        for (u, v), b in right.items():
            if s + u < prime and y + v < prime:
                result[s + u, y + v] = (result[s + u, y + v] + a * b) % prime
    return {key: value for key, value in result.items() if value}


def wick_mod(polynomial, prime):
    return sum(
        a * factorial(s) * (factorial(y) // (2 ** (y // 2) * factorial(y // 2)))
        for (s, y), a in polynomial.items()
        if y % 2 == 0
    ) % prime


def frobenius_controls():
    cases = [
        (1, {0: 2, 2: 1, 3: 1}),
        (1, {1: 2, 3: 1}),
        (2, {2: 1, 3: 2}),
        (3, {2: 1, 4: 2}),
        (4, {0: 2, 1: 1, 3: 2}),
    ]
    residue_checks = 0
    recurrence_checks = 0
    failed_order_checks = 0
    for d, h in [*cases, (2, {3: 1, 4: 2})]:
        for prime in (7, 11, 13, 17, 19):
            n, slack = divmod(prime - 1, d)
            phase = {(d, j): a for j, a in h.items()} | {(d + 1, 0): 3, (d + 2, 1): 1}
            h_poly = {(0, j): a for j, a in h.items()}
            power = h_power = {(0, 0): 1}
            for _ in range(n):
                power = truncated_product(power, phase, prime)
                h_power = truncated_product(h_power, h_poly, prime)

            def functional(j):
                return wick_mod(truncated_product(h_power, {(0, j): 1}, prime), prime)

            observed = []
            for j in range(max(h) + 1):
                original = wick_mod(truncated_product(power, {(slack, j): 1}, prime), prime)
                assert original == factorial(prime - 1) * functional(j) % prime
                observed.append(original)
                residue_checks += 1
            if min(h) <= d:
                assert any(observed)
            else:
                assert not h_power and not any(observed)
                failed_order_checks += 1
            for j in range(2 * prime):
                degree_lowering = sum(
                    a * ((j + (n + 1) * k) * (functional(j + k - 1) if j + k else 0)
                         - functional(j + k + 1))
                    for k, a in h.items()
                )
                assert degree_lowering % prime == 0
                recurrence_checks += 1
            if d == 1 and min(h) == 0:
                assert sum(a * functional(j) for j, a in h.items()) % prime == h[0] ** prime % prime
            if d == 1 and min(h) == 1:
                assert observed[0] == (-h[1] ** (prime - 1) * int(sp.factorial2(prime - 2))) % prime
    return dict(primes=[7, 11, 13, 17, 19], radial_orders=[1, 2, 3, 4],
                exact_residue_comparisons=residue_checks,
                triangular_recurrence_comparisons=recurrence_checks,
                vanishing_truncations_outside_order_hypothesis=failed_order_checks)


def cubic_radial_tail_controls():
    s, y = sp.symbols("s y")
    primes = (13, 17, 19, 23, 29)

    def powers(poly, order, prime):
        polynomial = {(i, j): int(c) % prime for (i, j), c in sp.Poly(poly, s, y).terms()}
        result = {(0, 0): 1}
        for _ in range(order):
            result = truncated_product(result, polynomial, prime)
        return result

    checked_faces = 0
    fixtures = [
        (s**3 + (2 + y**2) * s**2 + (y**3 + y**4) * s, 6, 4),
        (s**3 + (y + y**2) * s**2 + (y**3 + y**4) * s, 5, 3),
        (s**3 + y**2 * s**2 + (y**3 + y**4) * s, 9, 5),
        (s**3 + y**4 * s, 12, 6),
        (s**3 + y**3 * s**2 + y**7 * s, 9, 4),
        (s**3 + y**2 * s**2 + y**8 * s, 14, 7),
        (s**3 + y**4 * s**2, 12, 5),
    ]
    for phase, l, times in fixtures:
        for p in primes:
            q, slack = divmod(p - 1, l)
            power = powers(phase, times * q, p)
            assert wick_mod(truncated_product(power, {(slack, slack): 1}, p), p)
            checked_faces += 1
    for j in (2, 3, 4):
        for p in primes:
            n = (j + 1) * (p - 1) // (3 * j) - 1
            k0 = 2 * n - (p - 1) // j
            marker_powers = [(p - 1 - n - k, p - 1 - 2 * j * n + j * k) for k in (k0, k0 + 1)]
            assert all(2 <= i <= 5 and 0 <= l < 2 * j for i, l in marker_powers)
            power = powers(s**3 + 2 * y**j * s**2 + y**(2*j) * s, n, p)
            assert any(wick_mod(truncated_product(power, {marker: 1}, p), p) for marker in marker_powers)
            checked_faces += 1

    shifted_checks = 0
    for extra, q in ((y**3, 3), (y**4, 4)):
        j = (s + y + y**2) ** 2 + extra
        # Construct an actual critical-section approximation in the original ring.
        section = sp.Integer(0)
        for k in range(1, q + 1):
            value = sp.expand(sp.diff(j, y).subs(y, section)).coeff(s, k)
            section -= value * s**k / 2
        assert sp.expand(j.subs(y, section)).coeff(s, q) != 0
        assert all(sp.expand(j.subs(y, section)).coeff(s, k) == 0 for k in range(q))
        for p in primes:
            count, slack = divmod(p - 1, 2 * (q + 1))
            # In these examples the critical-section coefficients are localized at 2.
            fixed_marker = sp.Poly(s**slack * (y - section)**slack, s, y)
            residues = {(i, k): int(c.p) * pow(int(c.q), -1, p) % p
                        for (i, k), c in fixed_marker.terms()}
            power = powers(s*j, (q + 2) * count, p)
            assert wick_mod(truncated_product(power, residues, p), p)
            shifted_checks += 1

    square_checks = 0
    for p in primes:
        for quadratic, cubic in ((2, 0), (pow(2, -1, p), 0), (2, 1)):
            b = y + quadratic * y**2 + cubic * y**3
            bp = {(0, 1): 1, (0, 2): quadratic, (0, 3): cubic}
            n = (p - 1) // 2
            actual = powers(s * (s+b)**2, n, p)
            b_power = powers(b, n, p)
            exp_b = {(0, 0): 1}
            b_to_i = {(0, 0): 1}
            for i in range(1, p):
                b_to_i = truncated_product(b_to_i, bp, p)
                for exponents, c in b_to_i.items():
                    exp_b[exponents] = (exp_b.get(exponents, 0) + c * pow(factorial(i), -1, p)) % p
            shifted = truncated_product(b_power, exp_b, p)
            first_marker_values = []
            top_degree = 5 if cubic else 2 if 2 * quadratic % p == 1 else 3
            for j in range(top_degree):
                marker = {(0, j): 1}
                expected = factorial(p - 1) * (-1)**n * wick_mod(truncated_product(shifted, marker, p), p)
                observed = wick_mod(truncated_product(actual, marker, p), p)
                assert observed == expected % p
                first_marker_values.append(observed)
                square_checks += 1
            assert any(first_marker_values)
    return dict(weighted_face_checks=checked_faces,
                polynomial_critical_section_checks=shifted_checks,
                shifted_gamma_square_checks=square_checks)


def cubic_controls():
    s, u, y, a, r, k = sp.symbols("s u y a r k", nonzero=True)
    A, H, G = sp.symbols("A H G")
    C = A / (3 * a)
    x = k * C
    xi = (H - r - A**2 / (3 * a)) / (r * k)
    offset = G + (2 * A**3 / (27 * a**2) - A * (H - r) / (3 * a)) / r
    cubic = (a * s**3 + A * s**2 + (H - r) * s) / r + G
    assert sp.expand((offset + u**3 / 3 + xi * u).subs(u, k * (s + C))).subs(k**3, 3*a/r).expand() == cubic.expand()
    assert sp.expand(offset + x**3 / 3 + xi * x).subs(k**3, 3*a/r).expand() == G

    x, xi, scale, dx, dl, dxi = sp.symbols("x xi scale dx dl dxi")
    e_B = -(u**2 + xi)

    def derivative(f):
        return sp.diff(f, u) + (u**2 + xi) * f

    def reduce(f):
        f = sp.Poly(sp.expand(f), u)
        while f.degree() >= 3:
            power = u ** (f.degree() - 3)
            f = sp.Poly(f.as_expr() - f.LC() * derivative((u - x) * power), u)
        return f.as_expr()

    moving = dx + scale * (u - x)
    dG = dl + x * dxi + (x**2 + xi) * dx
    for f, df, wanted in (
        (sp.Integer(1), 0, dl + dxi * u - dx * e_B),
        (u, 0, dl * u - xi * dxi - (x * dx + dxi) * e_B),
        (e_B, -dxi, dG * e_B),
    ):
        actual = df + (dl + u * dxi) * f + derivative(moving * f)
        assert sp.expand(reduce(actual - wanted)) == 0

    x, xi, offset = [sp.Function(name)(y, r) for name in ("x", "xi", "Lambda")]
    endpoint = offset + x**3 / 3 + xi * x

    def connection(var):
        Dx = sp.diff(x, var)
        Da = sp.diff(offset, var)
        Dz = sp.diff(xi, var)
        return sp.Matrix([
            [sp.diff(endpoint, var), -Dx, -x * Dx - Dz],
            [0, Da, -xi * Dz],
            [0, Dz, Da],
        ])

    Dy, Dr = connection(y), connection(r)
    curvature = Dr.diff(y) - Dy.diff(r) + Dy * Dr - Dr * Dy
    assert all(sp.expand(entry) == 0 for entry in curvature)
    return dict(cubic_completion=True, moving_endpoint_connections=3,
                flat_curvature_entries=9, original_seed_factor="1/kappa")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    here = Path(__file__).resolve().parent
    root = here.parents[2]
    paths = [Path(__file__).resolve(), *[here / name for name in (
        "GAUSSIAN_RADIAL_FROBENIUS_AND_CUBIC_FIBRES.md",
        "GAUSSIAN_LOGARITHM_EXCLUSION_AND_RATIONAL_PRIMITIVES.md",
        "GAUSSIAN_INVARIANT_CORE_BOUNDARY_REDUCTION.md",
        "GAUSSIAN_ANGULAR_CORE_FULL_IMAGE_INVARIANCE.md",
        "GAUSSIAN_QUADRATIC_RADIAL_UNRESTRICTED_FINITE_IMAGES.md",
        "GAUSSIAN_VARIABLE_RADIAL_QUADRATIC_FIBRE_FINITE_IMAGES.md",
        "sources/stacks-054K-radial-reductions.html",
    )], root / "unpublished/gaussian-boundary/analytic_core.tex"]
    report = dict(status="passed", sympy=sp.__version__, frobenius=frobenius_controls(),
                  cubic_radial_tails=cubic_radial_tail_controls(),
                  cubic=cubic_controls(), sha256={str(path.relative_to(root)):
                      hashlib.sha256(path.read_bytes()).hexdigest() for path in paths},
                  scope="Bounded exact controls. The fixed-marker, source and uniform fibre estimates are written proofs.",
                  independent_replay=False, formal_verification=False, external_review=False)
    result = json.dumps(report, indent=2) + "\n"
    if args.output:
        args.output.write_text(result, encoding="utf-8", newline="\n")
    print(result, end="")


if __name__ == "__main__":
    main()

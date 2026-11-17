"""Bounded exact controls for the W-affine c=0 generic original image.

The unrestricted selected-field, original-frame and kernel statements are
proved in the companion note. These are finite controls, not degree-prefix
evidence for the all-degree theorem.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from math import factorial
from pathlib import Path
import platform

import sympy as sp

z, w, y, t = sp.symbols("z w y t")


def zero(expr):
    assert sp.cancel(sp.expand(expr)) == 0, expr


def trunc(expr, order):
    return sp.expand(sp.series(expr, t, 0, order + 1).removeO())


def wick(expr):
    total = 0
    for (iz, iw, iy), coefficient in sp.Poly(sp.expand(expr), z, w, y).terms():
        if iz == iw and iy % 2 == 0:
            total += coefficient * factorial(iz) * sp.Rational(
                factorial(iy), 2 ** (iy // 2) * factorial(iy // 2))
    return sp.expand(total)


def raw_periods(phase, markers, order):
    rows = [sp.Integer(0) for _ in markers]
    power = sp.Integer(1)
    for n in range(order + 1):
        if n:
            power = sp.expand(power * phase)
        for j, marker in enumerate(markers):
            rows[j] += wick(marker * power) * t**n / factorial(n)
    return [sp.expand(row) for row in rows]


def selected_root(h, order):
    value = sp.Integer(0)
    for _ in range(order + 1):
        value = trunc(t * h.subs(z, value), order)
    zero(trunc(value - t * h.subs(z, value), order))
    return value


def normalized_monomial(i, j, k, h, u, kk, aa):
    d = 1 - 2 * aa * t
    jac = 1 - t * sp.diff(h, z)
    mu, variance = t * u / d, 1 / d
    moment = sum(
        sp.Rational(factorial(k), 2**v * factorial(v) * factorial(k - 2*v))
        * variance**v * mu**(k - 2*v)
        for v in range(k // 2 + 1)
    )
    exponent_z = t * sp.diff(kk, z) + t**2 * u * sp.diff(u, z) / d
    correction = exponent_z / jac + t * sp.diff(h, z, 2) / jac**2
    value = sp.cancel(z**i * moment)
    for _ in range(j):
        value = sp.cancel(sp.diff(value, z) / jac + correction * value)
    return value


def normalized_marker(marker, h, u, kk, aa):
    return sp.cancel(sum(
        coefficient * normalized_monomial(i, j, k, h, u, kk, aa)
        for (i, j, k), coefficient in sp.Poly(marker, z, w, y).terms()
    ))


def field_reduce(expr, h):
    """Reduce in the SELECTED field, never the full reducible quotient."""
    if h.subs(z, 0) == 0:
        return sp.cancel(expr.subs(z, 0))
    q = sp.Poly(z - t*h, z, domain=sp.QQ.frac_field(t))
    numerator, denominator = sp.fraction(sp.cancel(expr))
    n = sp.Poly(numerator, z, domain=q.domain).rem(q)
    den = sp.Poly(denominator, z, domain=q.domain).rem(q)
    assert not den.is_zero
    inverse = sp.invert(den, q)
    return sp.cancel((n*inverse).rem(q).as_expr())


def check_original_sources():
    markers = [sp.Integer(1), z, z**2, y, y**2, w, z*w, w*y, w**2]
    cases = [
        ("zero-h", sp.Integer(0), sp.Rational(1, 3), 1+z, z**2),
        ("selected-zero-factor", z+z**3, sp.Integer(0), 1+z, z+z**2),
        ("nonzero-constant-h", sp.Integer(2), sp.Integer(0), z+1, z**2-z),
        ("rational-degree-one", 2+3*z, sp.Rational(1, 4), z+1, z**2-z),
        ("quadratic-field", 1+z**2, sp.Rational(1, 3), 1+z, z),
        ("cubic-field", 1+z+z**3, sp.Integer(0), 1+z**2, z**2),
    ]
    order = 4
    results = []
    for name, h, aa, u, kk in cases:
        phase = w*h + aa*y**2 + u*y + kk
        root = selected_root(h, order)
        d = 1 - 2*aa*t
        jac = trunc(1 - t*sp.diff(h, z).subs(z, root), order)
        exponent = trunc(t*kk.subs(z, root) + t**2*u.subs(z, root)**2/(2*d), order)
        f = trunc(sp.exp(exponent)/(jac*sp.sqrt(d)), order)
        actual = raw_periods(phase, markers, order)
        for marker, row in zip(markers, actual):
            ratio = normalized_marker(marker, h, u, kk, aa)
            ratio_series = trunc(ratio.subs(z, root), order)
            zero(row - trunc(f*ratio_series, order))
        # Derivative compatibility is exact in the selected algebraic field.
        jac_z = 1 - t*sp.diff(h, z)
        exponent_z = t*kk + t**2*u**2/(2*d)
        delta = lambda value: sp.diff(value, t) + h/jac_z * sp.diff(value, z)
        tau = delta(exponent_z) - delta(jac_z)/jac_z - sp.diff(d, t)/(2*d)
        zero(field_reduce(normalized_marker(phase, h, u, kk, aa)-tau, h))
        rank = 1 if h.subs(z, 0) == 0 else max(1, int(sp.degree(h, z)))
        if h.subs(z, 0) == 0:
            assert root == 0
        else:
            factors = sp.factor_list(z-t*h, z)[1]
            assert len(factors) == 1 and factors[0][1] == 1
            assert int(sp.degree(factors[0][0], z)) == rank
        results.append(dict(name=name, h=str(h), a=str(aa), u=str(u), k=str(kk),
                            selected_root_through_t4=str(root), selected_generic_rank=rank,
                            original_markers=[str(q) for q in markers],
                            Wick_orders=[0, order], exact_connection_identity=True))
    return results


def check_kernel_and_frame():
    h = 1+z**2
    zero_args = (h, sp.Integer(0), sp.Integer(0), sp.Integer(0))
    kernel_marker = w*(1-z**2)
    kernel_value = normalized_marker(kernel_marker, *zero_args)
    zero(field_reduce(kernel_value, h))
    nonzero_value = field_reduce(normalized_marker(z, *zero_args), h)
    zero(nonzero_value-z)
    assert nonzero_value != 0
    # The entire defining polynomial is reducible when h(0)=0; its other
    # factors must not be included in the selected actual generic image.
    h_zero = z+z**3
    q_full = sp.Poly(z-t*h_zero, z, domain=sp.QQ.frac_field(t))
    assert q_full.degree() == 3
    assert field_reduce(z, h_zero) == 0
    assert sp.Poly(z, z, domain=q_full.domain).rem(q_full).as_expr() == z
    return dict(exact_kernel_marker=str(kernel_marker), exact_kernel_reduction="0",
                exact_nonzero_marker="Z -> z in the quadratic selected field",
                selected_zero_negative_control="whole degree-three quotient is not the selected rank-one field")


def check_unit_rank_defect():
    h = 1+z**2
    phase = w*h
    order = 8
    f, phi_z, kernel = raw_periods(phase, [sp.Integer(1), z, w*(1-z**2)], order)
    delta = 1-4*t**2
    expected_f = trunc(delta**sp.Rational(-1, 2), order)
    zero(f-expected_f)
    zero(f-2*t*phi_z-1)
    zero(trunc(sp.diff(f, t)-4*t/delta*f, order-1))
    zero(kernel)
    # Each finite zero of Delta has odd order, so Delta is not a rational square.
    factors = sp.factor_list(delta, t)[1]
    assert len(factors) == 2 and all(multiplicity == 1 for _, multiplicity in factors)
    # Verify the defining field relation and two-element basis degree.
    q = sp.Poly(z-t*h, z, domain=sp.QQ.frac_field(t))
    assert q.degree() == 2 and q.is_irreducible
    return dict(phase=str(phase), original_Wick_orders=[0, order],
                unit_ODE="(1-4t^2) f' - 4t f = 0",
                unit_cyclic_rank=1, full_actual_generic_rank=2,
                frame=["Phi[1]", "Phi[Z]"],
                exact_time_amplitude="Phi[1-2tZ]=1",
                kernel_marker="Phi[W(1-Z^2)]=0",
                generic_not_global_scope_boundary="P=Y^2 has rank one and unbounded normalized finite-pole orders")


def checks():
    return dict(
        status="passed",
        original_source_controls=check_original_sources(),
        exact_selected_field_kernel_controls=check_kernel_and_frame(),
        unit_vs_full_rank_control=check_unit_rank_defect(),
        scope=("Bounded exact Wick, field, kernel and connection controls only. "
               "The companion proof establishes all h,u,k, all constant a, all "
               "characteristic-zero fields, and all original polynomial markers. "
               "No global K[t]-lattice or literal filtered-image classification is inferred.")
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = checks()
    directory = Path(__file__).resolve().parent
    names = ["GAUSSIAN_W_AFFINE_C_ZERO_GENERIC_IMAGE.md", Path(__file__).name,
             "GAUSSIAN_MULTIBASE_QUADRATIC_FIBRE_PRIMITIVITY.md",
             "GAUSSIAN_QUADRATIC_W_FINITE_IMAGE_CLASSIFICATION.md"]
    result.update(python=platform.python_version(), sympy=sp.__version__,
                  sha256={name: hashlib.sha256((directory/name).read_bytes()).hexdigest()
                          for name in names},
                  independent_replay=False, formal_verification=False, external_review=False)
    report = json.dumps(result, indent=2)+"\n"
    if args.output is not None:
        args.output.write_text(report, encoding="utf-8", newline="\n")
    print(report, end="")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Exact and numerical controls for the second-line counting theorem.

Finite content controls and area quadrature are explicitly not proofs of the
unrestricted asymptotic. Its analytic input is the retained pair estimate.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
from pathlib import Path

import mpmath as mp
import sympy as sp


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    checks: list[dict[str, object]] = []

    def require(name: str, condition: bool, **details) -> None:
        if not condition:
            raise AssertionError(name)
        checks.append({"name": name, "passed": True, **details})

    def equal(name: str, left, right=0) -> None:
        difference = sp.factor(left - right)
        if difference != 0:
            raise AssertionError(f"{name}: {difference}")
        checks.append({"name": name, "passed": True})

    D, K = 32768, 27
    m, n, r = sp.symbols("m n r")
    A = D * m**2 + K * n**2
    rr = A / (D * m * n)
    forms = [4608 * m * n, 9216 * m * n,
             16 * (5376 * m * n - A), -3 * (1792 * m * n + A)]
    equal("coordinate_u", forms[1], 2 * forms[0])
    equal("coordinate_b", forms[2] / forms[0], sp.Rational(56, 3) - sp.Rational(1024, 9) * rr)
    equal("coordinate_c", forms[3] / forms[0], -sp.Rational(7, 6) - sp.Rational(64, 3) * rr)
    equal("coercive_combination", forms[2] + 16 * forms[3], -64 * A)

    def valuation(value: int, prime: int) -> int:
        value = abs(value)
        if not value:
            raise ValueError("valuation of zero was not requested")
        order = 0
        while value % prime == 0:
            order += 1
            value //= prime
        return order

    def content_control(mm: int, nn: int) -> None:
        aa = D * mm**2 + K * nn**2
        entries = [4608 * mm * nn, 9216 * mm * nn,
                   16 * (5376 * mm * nn - aa), -3 * (1792 * mm * nn + aa)]
        actual = math.gcd(*entries)
        expected = 2**min(2 * valuation(nn, 2), 15) * 3**min(2 * valuation(mm, 3), 3)
        if actual != expected:
            raise AssertionError((mm, nn, actual, expected))

    finite_pairs = 0
    for mm in range(-100, 101):
        for nn in range(1, 101):
            if mm and math.gcd(mm, nn) == 1:
                content_control(mm, nn)
                finite_pairs += 1
    require("primitive_content_bounded_control", finite_pairs == 12174, coprime_pairs=finite_pairs)
    valuation_pairs = 0
    for ee in range(16):
        for ff in range(10):
            for sign in (-1, 1):
                for unit_m, unit_n in ((1, 1), (5, 7), (7, 5)):
                    for mm, nn in ((sign * 3**ff * unit_m, 2**ee * unit_n),
                                   (sign * 2**ee * unit_m, 3**ff * unit_n)):
                        if math.gcd(mm, nn) == 1:
                            content_control(mm, nn)
                            valuation_pairs += 1
    require("primitive_content_large_valuation_controls", valuation_pairs > 1000, coprime_pairs=valuation_pairs)

    W2 = 1 + sp.Rational(1, 7) + 8 + 64 + sp.Rational(512, 7)
    W3 = 1 + sp.Rational(27, 26) + sp.Rational(1, 26)
    equal("two_adic_weight", W2, sp.Rational(1024, 7))
    equal("three_adic_weight", W3, sp.Rational(27, 13))
    equal("total_valuation_weight", W2 * W3, sp.Rational(27648, 91))
    equal("full_euler_prefactor", sp.Rational(7, 36) / 3 * W2 * W3, sp.Rational(256, 13))
    equal("log_and_involution_factors_cancel", sp.Rational(2, 1) * sp.Rational(1, 2), 1)

    R0 = 768 * sp.sqrt(6)
    lower, upper = sp.Rational(80640, 19), sp.Rational(91392, 13)
    equal("AM_GM_threshold", R0**2, 4 * D * K)
    equal("lower_breakpoint", 16 * (5376 - lower), 3 * (lower + 1792))
    equal("upper_breakpoint", 16 * (upper - 5376), 3 * (upper + 1792))
    require("breakpoint_order", bool(1792 < R0 < lower < 5376 < upper))
    require("second_coordinate_dominated_positive", bool(R0 > 1280))
    equal("logarithmic_integral_discriminant", 5376**2 - R0**2, (768 * sp.sqrt(43))**2)
    equal("arctangent_integral_discriminant", R0**2 - 1792**2, (256 * sp.sqrt(5))**2)

    # Independent numerical area controls: the R integral derived in the
    # source and a polar-coordinate integral directly from the height Q.
    mp.mp.dps = 60
    root = 768 * mp.sqrt(6)
    lo, hi = mp.mpf(80640) / 19, mp.mpf(91392) / 13
    q = mp.sqrt((5376 - root) / (5376 + root))
    u_lo = mp.sqrt((lo - root) / (lo + root))
    u_hi = mp.sqrt((hi - root) / (hi + root))
    eta = mp.sqrt((root - 1792) / (root + 1792))
    closed = (mp.log(((u_hi + q) / (u_hi - q)) * ((q + u_lo) / (q - u_lo))) / (12288 * mp.sqrt(43))
              + (mp.atan(u_hi * eta) - mp.atan(u_lo * eta)) / (384 * mp.sqrt(5)))
    # R=R0*cosh(s) removes the integrable endpoint singularity.
    slo, shi = mp.acosh(lo / root), mp.acosh(hi / root)
    neg_area = mp.quad(lambda s: 1 / (16 * (root * mp.cosh(s) + 5376)), [0, slo, shi, mp.inf])
    pos_area = (mp.quad(lambda s: 1 / (16 * (5376 - root * mp.cosh(s))), [0, slo])
                + mp.quad(lambda s: 1 / (3 * (root * mp.cosh(s) + 1792)), [slo, shi])
                + mp.quad(lambda s: 1 / (16 * (root * mp.cosh(s) - 5376)), [shi, mp.inf]))
    r_integral = neg_area + pos_area

    def q_polar(theta):
        xx, yy = mp.cos(theta), mp.sin(theta)
        aa = D * xx**2 + K * yy**2
        return max(9216 * abs(xx * yy), 16 * abs(5376 * xx * yy - aa), 3 * abs(1792 * xx * yy + aa))

    t_breaks = [(v + sign * mp.sqrt(v**2 - root**2)) / (2 * D)
                for v in (lo, hi) for sign in (-1, 1)]
    theta_breaks = sorted([mp.mpf(0), mp.pi / 2, mp.pi]
                          + [mp.atan(1 / tt) for tt in t_breaks])
    polar_integral = mp.quad(lambda theta: 1 / (2 * q_polar(theta)), theta_breaks)
    require("closed_area_matches_R_integral", abs(closed - r_integral) < mp.mpf("1e-45"))
    require("independent_polar_area", abs(polar_integral - r_integral) < mp.mpf("1e-45"))
    require("both_parameter_signs_positive", neg_area > 0 and pos_area > 0)

    # Rational-point parameter height controls, for the proved upper bound.
    cubic_controls = 0
    for mm in range(-40, 41):
        for nn in range(1, 41):
            if math.gcd(mm, nn) != 1:
                continue
            numerator = 32 * mm**3 - 9 * mm * nn**2
            denominator = 32 * nn**3
            common = math.gcd(numerator, denominator)
            maximum = max(abs(mm), nn)
            require_bound = (32 % common == 0
                             and 23 * maximum**3 <= max(abs(numerator), denominator)
                             and max(abs(numerator), denominator) <= 41 * maximum**3)
            if not require_bound:
                raise AssertionError((mm, nn, common))
            cubic_controls += 1
    require("rational_point_cubic_height_controls", cubic_controls > 1000, coprime_pairs=cubic_controls)

    # The injectivity obstruction and exact primitive content for all
    # globally soluble targets. These bounded controls supplement the proof.
    require("injectivity_normalized_quadratic_form_odd",
            all((u*u + u*v + v*v) % 2 == 1
                for u, v in ((1, 0), (0, 1), (1, 1))))
    require("injectivity_RHS_valuation_odd", -5 % 2 == 1)
    global_forms = [18*n**3, 36*n**3,
                    336*n**3 - 2048*m**3 + 576*m*n**2,
                    -21*n**3 - 384*m**3 + 108*m*n**2]
    root_r = m**3/n**3 - sp.Rational(9, 32)*m/n
    equal("global_target_coordinate_b", global_forms[2]/global_forms[0],
          sp.Rational(56, 3) - sp.Rational(1024, 9)*root_r)
    equal("global_target_coordinate_c", global_forms[3]/global_forms[0],
          -sp.Rational(7, 6) - sp.Rational(64, 3)*root_r)
    global_finite = 0
    for mm in range(-100, 101):
        for nn in range(1, 101):
            if math.gcd(mm, nn) != 1:
                continue
            values = [18*nn**3, 36*nn**3,
                      336*nn**3 - 2048*mm**3 + 576*mm*nn**2,
                      -21*nn**3 - 384*mm**3 + 108*mm*nn**2]
            order_two = valuation(nn, 2)
            ee = 0 if order_two == 0 else 3 if order_two == 1 else 7
            expected = 2**ee * (3 if mm % 3 == 0 else 1)
            if math.gcd(*values) != expected:
                raise AssertionError((mm, nn, math.gcd(*values), expected))
            global_finite += 1
    require("global_primitive_content_controls", global_finite == 12175,
            coprime_pairs=global_finite)
    local_two = sp.Rational(1, 2) + sp.Rational(1, 8)*4 + sp.Rational(1, 8)*2**sp.Rational(14, 3)
    local_three = sp.Rational(2, 3) + sp.Rational(2, 9)*3**sp.Rational(2, 3)
    equal("global_two_adic_weight", local_two, 1 + 2**sp.Rational(5, 3))
    equal("global_three_adic_weight", local_three,
          sp.Rational(2, 9)*(3 + 3**sp.Rational(2, 3)))
    equal("global_constant_weight", sp.Rational(3, 2)*local_two*local_three,
          (1 + 2**sp.Rational(5, 3))*(3 + 3**sp.Rational(2, 3))/3)
    require("rational_zero_target_height", max(18, 36, 336, 21)//3 == 112)

    # Independent cubic archimedean area evaluations, in t and polar angle.
    cubic_f = lambda tt: tt**3 - mp.mpf(9)/32*tt
    global_lower_r = mp.mpf(315)/2432
    global_upper_r = mp.mpf(357)/1664
    global_lo = mp.findroot(lambda tt: cubic_f(tt) - global_lower_r,
                            (mp.mpf(".5"), mp.mpf(".8")))
    global_hi = mp.findroot(lambda tt: cubic_f(tt) - global_upper_r,
                            (mp.mpf(".5"), mp.mpf(".9")))
    global_area_t = (
        mp.quad(lambda tt: (336 - 2048*cubic_f(tt))**(-mp.mpf(2)/3),
                [-mp.inf, 0, global_lo])
        + mp.quad(lambda tt: (21 + 384*cubic_f(tt))**(-mp.mpf(2)/3),
                  [global_lo, global_hi])
        + mp.quad(lambda tt: (2048*cubic_f(tt) - 336)**(-mp.mpf(2)/3),
                  [global_hi, mp.inf]))/2

    def cubic_q_polar(theta):
        xx, yy = mp.cos(theta), mp.sin(theta)
        return max(36*yy**3, abs(336*yy**3 - 2048*xx**3 + 576*xx*yy**2),
                   abs(21*yy**3 + 384*xx**3 - 108*xx*yy**2))

    global_theta_breaks = [mp.mpf(0), mp.atan(1/global_hi),
                           mp.atan(1/global_lo), mp.pi/2, mp.pi]
    global_area_polar = mp.quad(
        lambda theta: cubic_q_polar(theta)**(-mp.mpf(2)/3)/2,
        global_theta_breaks)
    require("independent_global_cubic_area",
            abs(global_area_t-global_area_polar) < mp.mpf("1e-45"))
    global_constant = (global_area_t*(1+2**(mp.mpf(5)/3))
                       *(3+3**(mp.mpf(2)/3))/(3*mp.zeta(2)))
    # A finite Euler prefix is only a numerical orientation for the exact
    # positive product in the theorem. The following tail bound is rigorous.
    cutoff = 100000
    product = mp.mpf(1)
    for prime in sp.primerange(2, cutoff + 1):
        if prime % 3 == 2:
            term = mp.mpf(prime)**-3
            product *= (1 + term) / (1 - term)
    const_prefix = 256 * closed * product / (13 * mp.sqrt(3) * mp.zeta(2))
    relative_tail_bound = mp.expm1(1 / (mp.mpf(cutoff)**2 * (1 - mp.mpf(cutoff)**-3)))
    source_path = Path(__file__).with_name("FIXED_MAP_SECOND_LINE_COUNT.md")
    report = {
        "status": "passed",
        "scope": "Exact coefficient and constant controls, finite content controls, and numerical independent area checks; not a full analytic proof replay or external review.",
        "python": platform.python_version(),
        "sympy": sp.__version__,
        "mpmath": mp.__version__,
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "source_sha256": hashlib.sha256(source_path.read_bytes()).hexdigest(),
        "check_count": len(checks),
        "checks": checks,
        "numerical_controls": {
            "precision_decimal_digits": mp.mp.dps,
            "area_closed": mp.nstr(closed, 35),
            "area_R_integral": mp.nstr(r_integral, 35),
            "area_polar": mp.nstr(polar_integral, 35),
            "negative_parameter_area": mp.nstr(neg_area, 35),
            "positive_parameter_area": mp.nstr(pos_area, 35),
            "global_cubic_area_t": mp.nstr(global_area_t, 35),
            "global_cubic_area_polar": mp.nstr(global_area_polar, 35),
            "global_leading_constant": mp.nstr(global_constant, 30),
            "Euler_prime_cutoff": cutoff,
            "constant_Euler_prefix": mp.nstr(const_prefix, 30),
            "Euler_relative_tail_bound": mp.nstr(relative_tail_bound, 15),
        },
    }
    if args.output:
        args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": report["status"], "check_count": len(checks),
                      "constant_approximation": mp.nstr(const_prefix, 13),
                      "output": str(args.output) if args.output else None}))


if __name__ == "__main__":
    main()

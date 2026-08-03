#!/usr/bin/env python3
"""Exact regressions for the ternary cusp-profile suspension theorem.

The companion note gives the all-order proof.  This checker uses only
integer and rational arithmetic to verify:

* the cusp-defect and homogeneous-divisibility identities;
* the full phase-coefficient ladder for several windings and profiles;
* top Laplacian contractions, exact trace depth, and shifted-power
  detectors in bounded non-power examples;
* the original, endpoint-power, and radial-padding specializations.
"""

from __future__ import annotations

import importlib.util
import json
from fractions import Fraction
from math import comb, factorial
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BASE_PATH = ROOT / "scripts" / "verify_gvc3_homogeneous_counterexample.py"
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "gvc3_cusp_profile_suspension.json"
)

SPEC = importlib.util.spec_from_file_location("gvc3_base", BASE_PATH)
assert SPEC and SPEC.loader
base = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(base)


def polynomial_sum(*polynomials):
    """Call the base sparse sum, including for an empty input."""
    return base.add(*polynomials) if polynomials else {}


def polynomial_product(*polynomials):
    """Multiply a nonempty list of base sparse polynomials."""
    assert polynomials
    result = {base.ZERO: 1}
    for polynomial in polynomials:
        result = base.multiply(result, polynomial)
    return result


def double_factorial(value: int) -> int:
    result = 1
    for entry in range(value, 0, -2):
        result *= entry
    return result


def univariate_multiply(left: list[int], right: list[int]) -> list[int]:
    result = [0] * (len(left) + len(right) - 1)
    for left_degree, left_coefficient in enumerate(left):
        for right_degree, right_coefficient in enumerate(right):
            result[left_degree + right_degree] += (
                left_coefficient * right_coefficient
            )
    return result


def univariate_power(polynomial: list[int], exponent: int) -> list[int]:
    result = [1]
    factor = polynomial
    remaining = exponent
    while remaining:
        if remaining & 1:
            result = univariate_multiply(result, factor)
        remaining //= 2
        if remaining:
            factor = univariate_multiply(factor, factor)
    return result


def endpoint_profile_power(r: int, coefficients: list[int], m: int) -> list[int]:
    """Coefficients of ((1-z)^(2r) S(z))^m."""
    endpoint = [(-1) ** j * comb(2 * r, j) for j in range(2 * r + 1)]
    return univariate_power(
        univariate_multiply(endpoint, coefficients),
        m,
    )


def endpoint_moment(r: int, coefficients: list[int], m: int) -> Fraction:
    profile = endpoint_profile_power(r, coefficients, m)
    return sum(
        (Fraction(coefficient, 2 * degree + 1)
         for degree, coefficient in enumerate(profile)),
        Fraction(0),
    )


def phase_coefficient(
    r: int,
    coefficients: list[int],
    m: int,
    u_degree: int,
) -> Fraction:
    """Coefficient of u^d in the phase-integrated H_m(1+u)."""
    profile = endpoint_profile_power(r, coefficients, m)
    return sum(
        (
            Fraction(coefficient, 2 * degree + 1)
            * comb(r * m + 2 * degree, u_degree)
            for degree, coefficient in enumerate(profile)
            if r * m + 2 * degree >= u_degree
        ),
        Fraction(0),
    )


def build_geometry():
    rho = base.add(base.monomial((0, 0, 2)), base.monomial((1, 1, 0)))
    a = base.add(rho, base.monomial((2, 0, 0)))
    c = base.add(
        base.multiply(base.Y, base.power(rho, 2)),
        base.scale(
            base.multiply(base.multiply(base.X, base.power(base.T, 2)), rho),
            -2,
        ),
        base.scale(base.monomial((3, 0, 2)), -1),
    )
    delta = base.add(
        base.monomial((1, 1, 0), 4),
        base.monomial((0, 0, 2)),
    )
    return rho, a, c, delta


def build_profile(
    r: int,
    coefficients: list[int],
    h: int,
):
    assert r >= 1 and h >= 0 and coefficients and any(coefficients)
    rho, a, c, _ = build_geometry()
    e = len(coefficients) - 1
    t2a2 = base.multiply(base.power(base.T, 2), base.power(a, 2))
    homogenized_terms = []
    for j, coefficient in enumerate(coefficients):
        if not coefficient:
            continue
        term = base.multiply(
            base.power(t2a2, j),
            base.power(rho, 3 * (e - j)),
        )
        homogenized_terms.append(base.scale(term, coefficient))
    s_hom = polynomial_sum(*homogenized_terms)
    p = polynomial_product(
        base.power(rho, h),
        base.power(a, r),
        base.power(c, 2 * r),
        s_hom,
    )
    n = 6 * r + 3 * e + h
    if e:
        assert base.homogeneous_degree(s_hom) == 6 * e
    else:
        assert s_hom == {base.ZERO: coefficients[0]}
    assert base.homogeneous_degree(p) == 2 * n
    return p, s_hom, n


def expected_top_contraction(
    n: int,
    r: int,
    m: int,
    ell: int,
    c_m: Fraction,
) -> Fraction:
    return (
        2 ** (n * m + ell)
        * factorial(n * m + ell)
        * double_factorial(2 * n * m + 2 * ell + 1)
        * comb(r * m - 1, ell - 1)
        * c_m
    )


def expected_terminal_trace(n: int, m: int, c_m: Fraction) -> Fraction:
    return (
        2 ** (n * m - 3)
        * factorial(n * m - 1)
        * double_factorial(2 * n * m + 3)
        * c_m
    )


def expected_shifted_detector(
    n: int,
    r: int,
    m: int,
    d: int,
    c_shifted: Fraction,
) -> Fraction:
    return (
        2 ** (n * m + (n - 3) * d)
        * factorial(n * m + (n - 1) * d)
        * double_factorial(2 * n * m + 2 * (n + 1) * d + 1)
        * comb(r * (m + d) - 1, d - 1)
        * c_shifted
    )


def verify_phase_ladders():
    records = []
    profiles = [
        (1, [1], "constant"),
        (1, [1, 1], "positive non-power 1+z"),
        (1, [1, -1, 1], "positive quadratic 1-z+z^2"),
        (2, [1], "winding two"),
        (2, [2, 1], "winding two with profile 2+z"),
    ]
    for r, coefficients, label in profiles:
        for m in range(1, 6):
            c_m = endpoint_moment(r, coefficients, m)
            assert c_m != 0
            pure = phase_coefficient(r, coefficients, m, r * m)
            assert pure == 0
            ladder = []
            for ell in range(1, r * m + 1):
                value = phase_coefficient(
                    r,
                    coefficients,
                    m,
                    r * m - ell,
                )
                expected = comb(r * m - 1, ell - 1) * c_m
                assert value == expected
                ladder.append(str(value))
            records.append(
                {
                    "profile": label,
                    "r": r,
                    "S_coefficients": coefficients,
                    "m": m,
                    "c_m": str(c_m),
                    "pure_coefficient": str(pure),
                    "lower_ladder": ladder,
                }
            )
    return records


def verify_polynomial_specializations():
    rho, a, c, _ = build_geometry()
    base_p = base.multiply(a, base.power(c, 2))

    original, _, original_n = build_profile(1, [1], 0)
    assert original_n == 6 and original == base_p

    padded = []
    for h in range(4):
        p, _, n = build_profile(1, [1], h)
        assert p == base.multiply(base.power(rho, h), base_p)
        assert n == 6 + h
        padded.append({"h": h, "N": n, "term_count": len(p)})

    endpoint_powers = []
    for s in range(2, 6):
        e = s - 2
        coefficients = [(-1) ** j * comb(e, j) for j in range(e + 1)]
        p, s_hom, n = build_profile(1, coefficients, 0)
        expected_hom = base.power(
            base.add(
                base.power(rho, 3),
                base.scale(
                    base.multiply(base.power(base.T, 2), base.power(a, 2)),
                    -1,
                ),
            ),
            e,
        )
        assert s_hom == expected_hom
        expected = polynomial_product(
            a,
            base.monomial((e, 0, 0)),
            base.power(c, s),
        )
        assert p == expected
        assert n == 3 * s
        endpoint_powers.append(
            {"s": s, "N": n, "term_count": len(p)}
        )

    return {
        "original": {"r": 1, "S": "1", "h": 0, "N": original_n},
        "radial_padding": padded,
        "endpoint_power_family": endpoint_powers,
    }


def verify_cusp_identity():
    rho, a, c, _ = build_geometry()
    defect = base.add(
        base.multiply(base.power(base.T, 2), base.power(a, 2)),
        base.scale(base.power(rho, 3), -1),
    )
    assert defect == base.scale(base.multiply(base.X, c), -1)
    return {
        "U": "rho=t^2+x*y",
        "V": "t*(rho+x^2)",
        "identity": "V^2-U^3=-x*C",
        "x_zero_normalization": "(U,V)=(t^2,t^3)",
    }


def verify_suspension_strata():
    records = []
    for k in range(6, 19):
        strata = []
        for r in range(1, k // 6 + 1):
            for e in range((k - 6 * r) // 3 + 1):
                h = k - 6 * r - 3 * e
                assert h >= 0
                strata.append({"r": r, "e": e, "h": h})
        assert {"r": 1, "e": 0, "h": k - 6} in strata
        if k == 6:
            assert strata == [{"r": 1, "e": 0, "h": 0}]
        records.append({"k": k, "strata": strata})

    # A trailing zero at profile level e is exactly three radial suspensions.
    with_trailing_zero, _, n_left = build_profile(1, [1, 1, 0], 0)
    with_radial_shift, _, n_right = build_profile(1, [1, 1], 3)
    assert n_left == n_right
    assert with_trailing_zero == with_radial_shift
    return {
        "powers": records,
        "trailing_profile_zero_equals_radial_shift": True,
    }


def verify_differential_replays():
    _, _, _, delta = build_geometry()
    cases = [
        (1, [1], 0, 2, "original"),
        (1, [1, 1], 0, 2, "positive non-power profile"),
        (1, [1], 1, 2, "radial suspension"),
        (2, [1], 0, 1, "winding two"),
    ]
    records = []
    for r, coefficients, h, max_m, label in cases:
        p, _, n = build_profile(r, coefficients, h)
        p_power = {base.ZERO: 1}
        for m in range(1, max_m + 1):
            p_power = base.multiply(p_power, p)
            c_m = endpoint_moment(r, coefficients, m)
            assert c_m != 0
            pure = base.apolar_scalar(base.power(delta, n * m), p_power)
            assert pure == 0
            ladder = []
            for ell in range(1, min(r * m, 2) + 1):
                mixed_input = base.multiply(
                    base.monomial((2 * ell, 0, 0)),
                    p_power,
                )
                actual = base.apolar_scalar(
                    base.power(delta, n * m + ell),
                    mixed_input,
                )
                expected = expected_top_contraction(n, r, m, ell, c_m)
                assert expected.denominator == 1
                assert actual == expected.numerator
                assert actual != 0
                ladder.append({"ell": ell, "scalar": str(actual)})

            terminal = base.apply_operator(
                base.power(delta, n * m - 1),
                p_power,
            )
            assert terminal
            terminal_detector = base.apolar_scalar(
                base.monomial((0, 2, 0)),
                terminal,
            )
            terminal_expected = expected_terminal_trace(n, m, c_m)
            assert terminal_expected.denominator == 1
            assert terminal_detector == terminal_expected.numerator
            assert terminal_detector != 0
            records.append(
                {
                    "case": label,
                    "r": r,
                    "S_coefficients": coefficients,
                    "h": h,
                    "N": n,
                    "m": m,
                    "term_count_P": len(p),
                    "pure": pure,
                    "top_ladder": ladder,
                    "terminal_trace_terms": len(terminal),
                    "terminal_dy2_detector": str(terminal_detector),
                }
            )

    # Replay the shifted-power statement on a genuinely non-power profile.
    r, coefficients, h = 1, [1, 1], 0
    p, _, n = build_profile(r, coefficients, h)
    shifted_records = []
    for m, d in ((1, 1), (1, 2)):
        shifted_power = base.power(p, m + d)
        output = base.apply_operator(base.power(delta, n * m), shifted_power)
        assert output
        detector = base.multiply(
            base.power(delta, (n - 1) * d),
            base.monomial((0, 2 * d, 0)),
        )
        actual = base.apolar_scalar(detector, output)
        expected = expected_shifted_detector(
            n,
            r,
            m,
            d,
            endpoint_moment(r, coefficients, m + d),
        )
        assert expected.denominator == 1
        assert actual == expected.numerator
        assert actual != 0
        shifted_records.append(
            {
                "m": m,
                "d": d,
                "output_terms": len(output),
                "detector": str(actual),
            }
        )

    return records, shifted_records


def main() -> None:
    phase_ladders = verify_phase_ladders()
    specializations = verify_polynomial_specializations()
    cusp = verify_cusp_identity()
    strata = verify_suspension_strata()
    differential, shifted = verify_differential_replays()

    artifact = {
        "format": "gvc3-cusp-profile-suspension-v1",
        "field": "characteristic zero",
        "family": {
            "rho": "t^2+x*y",
            "A": "rho+x^2",
            "C": "(rho^3-t^2*A^2)/x",
            "S_hom": "sum_j s_j*(t^2*A^2)^j*rho^(3*(e-j))",
            "P": "rho^h*A^r*C^(2r)*S_hom",
            "N": "6r+3e+h",
            "operator": "Lambda=(4*d_x*d_y+d_t^2)^N",
            "c_m": "integral_0^1 (1-v^2)^(2rm)*S(v^2)^m dv",
        },
        "all_order_formulas": {
            "pure": "Delta^(Nm)(P^m)=0",
            "ladder": (
                "Delta^(Nm+ell)(x^(2ell)*P^m)="
                "2^(Nm+ell)*(Nm+ell)!*(2Nm+2ell+1)!!*"
                "binom(rm-1,ell-1)*c_m"
            ),
            "terminal_trace": (
                "d_y^2*Delta^(Nm-1)(P^m)="
                "2^(Nm-3)*(Nm-1)!*(2Nm+3)!!*c_m"
            ),
            "shifted_power_detector": (
                "Delta^((N-1)d)*d_y^(2d)*Lambda^m(P^(m+d))="
                "2^(Nm+(N-3)d)*(Nm+(N-1)d)!*"
                "(2Nm+2(N+1)d+1)!!*binom(r(m+d)-1,d-1)*c_(m+d)"
            ),
        },
        "cusp_geometry": cusp,
        "rank_power_profile_strata": strata,
        "phase_ladder_replays": phase_ladders,
        "polynomial_specializations": specializations,
        "differential_replays": differential,
        "shifted_power_replays": shifted,
        "beta_rank": {
            str(n): max(n - 2, 0) for n in range(1, 8)
        },
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(artifact, indent=2) + "\n")

    print("PASS cusp-defect and homogeneous specialization identities")
    print("PASS admissible rank/power profile strata through k=18")
    print("PASS full phase ladder for five profiles through m=5")
    print("PASS top contractions and exact trace depths")
    print("PASS shifted-power detectors for a non-power profile")
    print(f"PASS wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

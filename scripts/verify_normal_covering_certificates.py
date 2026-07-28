#!/usr/bin/env python3
"""Verify faithful normal-cover certificates and their local witnesses."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
from typing import Any

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.normal_covering import analyze_component_action


DEFAULT_CERTIFICATES = (
    ROOT / "arithmetic/certificates/normal_cover_s3_quintic.json",
    ROOT / "arithmetic/certificates/normal_cover_v4_sextic.json",
)
SCHEMA = "normal-covering-certificate/v1"
T = sp.Symbol("T")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "certificates",
        nargs="*",
        type=Path,
        default=list(DEFAULT_CERTIFICATES),
    )
    parser.add_argument("--json", action="store_true")
    return parser.parse_args()


def polynomial_from_coefficients(coefficients: list[int]) -> sp.Poly:
    expression = sum(
        sp.Integer(coefficient) * T**degree
        for degree, coefficient in enumerate(coefficients)
    )
    return sp.Poly(expression, T, domain=sp.QQ)


def valuation(value: int, prime: int) -> int:
    value = abs(value)
    if value == 0:
        raise ValueError("the exact Hensel witnesses must have nonzero values")
    exponent = 0
    while value % prime == 0:
        exponent += 1
        value //= prime
    return exponent


def factor_record_map(arithmetic: dict[str, Any]) -> dict[str, sp.Poly]:
    factors: dict[str, sp.Poly] = {}
    for record in arithmetic["factors"]:
        name = record["name"]
        if name in factors:
            raise ValueError(f"duplicate factor name {name}")
        factor = polynomial_from_coefficients(
            record["coefficients_low_to_high"]
        )
        assert factor.degree() >= 2
        assert factor.is_irreducible
        assert factor.discriminant() == record["discriminant"]
        factors[name] = factor
    return factors


def verify_arithmetic(certificate: dict[str, Any]) -> dict[str, Any]:
    arithmetic = certificate["arithmetic"]
    polynomial = polynomial_from_coefficients(
        arithmetic["polynomial_coefficients_low_to_high"]
    )
    factors = factor_record_map(arithmetic)
    product = sp.Poly(1, T, domain=sp.QQ)
    for factor in factors.values():
        product *= factor
    assert product == polynomial
    assert polynomial.gcd(polynomial.diff()).degree() == 0
    assert all(factor.degree() >= 2 for factor in factors.values())

    ramified_primes = arithmetic["splitting_field_ramified_primes"]
    assert ramified_primes == sorted(set(ramified_primes))
    witnessed_primes: list[int] = []
    for witness in arithmetic["local_witnesses"]:
        prime = witness["prime"]
        assert prime in ramified_primes
        factor = factors[witness["factor"]]
        approximation = witness["approximation"]
        value = int(factor.eval(approximation))
        derivative_value = int(factor.diff().eval(approximation))
        value_valuation = valuation(value, prime)
        derivative_valuation = valuation(derivative_value, prime)
        assert value_valuation == witness["value_valuation"]
        assert derivative_valuation == witness["derivative_valuation"]
        if witness["criterion"] == "simple_hensel":
            assert value_valuation >= 1
            assert derivative_valuation == 0
        elif witness["criterion"] == "strong_hensel":
            assert value_valuation > 2 * derivative_valuation
        else:
            raise ValueError(f"unknown Hensel criterion {witness['criterion']}")
        witnessed_primes.append(prime)
    assert sorted(witnessed_primes) == ramified_primes

    real_witness = arithmetic["real_witness"]
    real_factor = factors[real_witness["factor"]]
    left, right = real_witness["interval"]
    left_value = real_factor.eval(left)
    right_value = real_factor.eval(right)
    assert left_value == 0 or right_value == 0 or left_value * right_value < 0
    return {
        "polynomial_degree": polynomial.degree(),
        "factor_degrees": sorted(factor.degree() for factor in factors.values()),
        "ramified_primes": ramified_primes,
        "local_witness_count": len(witnessed_primes),
    }


def verify_certificate(path: Path) -> dict[str, Any]:
    certificate = json.loads(path.read_text())
    assert certificate["schema"] == SCHEMA
    action = certificate["action"]
    analysis = analyze_component_action(
        degree=action["degree"],
        generators=action["generators"],
        components=action["components"],
    )
    expected = certificate["expected"]
    assert len(analysis.group) == expected["group_order"]
    assert list(analysis.factorization_shape) == expected["factorization_shape"]
    assert len(analysis.components) == expected["component_count"]
    assert analysis.is_normal_cover is expected["normal_cover"]
    assert len(analysis.common_core) == expected["common_core_order"]
    assert analysis.normal_covering_number == expected["normal_covering_number"]
    assert (
        len(analysis.components) == analysis.normal_covering_number
    ) is expected["covering_is_minimal"]
    assert analysis.index_sum == expected["index_sum"]
    assert analysis.index_sum == action["degree"]

    component_factors = {
        component["factor"].replace("^", "**")
        for component in action["components"]
    }
    arithmetic_factors = {
        factor["name"]: factor
        for factor in certificate["arithmetic"]["factors"]
    }
    assert component_factors == {
        "".join(str(polynomial_from_coefficients(record[
            "coefficients_low_to_high"
        ]).as_expr()).split())
        for record in arithmetic_factors.values()
    }
    arithmetic_summary = verify_arithmetic(certificate)
    assert arithmetic_summary["polynomial_degree"] == action["degree"]
    assert (
        arithmetic_summary["factor_degrees"]
        == expected["factorization_shape"]
    )
    return {
        "certificate": str(path.relative_to(ROOT)),
        "group_label": action["group_label"],
        "group_order": len(analysis.group),
        "factorization_shape": list(analysis.factorization_shape),
        "component_count": len(analysis.components),
        "normal_covering_number": analysis.normal_covering_number,
        "normal_cover": analysis.is_normal_cover,
        "common_core_order": len(analysis.common_core),
        "index_sum": analysis.index_sum,
        **arithmetic_summary,
    }


def main() -> None:
    args = parse_args()
    summaries = [
        verify_certificate(path.resolve())
        for path in args.certificates
    ]
    if args.json:
        print(json.dumps(summaries, indent=2, sort_keys=True))
        return
    for summary in summaries:
        print(
            "PASS "
            f"{summary['group_label']}: shape={summary['factorization_shape']} "
            f"r=gamma={summary['normal_covering_number']} "
            f"core={summary['common_core_order']} "
            f"ramified={summary['ramified_primes']}"
        )


if __name__ == "__main__":
    main()

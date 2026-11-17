#!/usr/bin/env python3
"""Minimal independent replay of the degree-42 ``c6`` certificate.

This checker uses only Python's standard library.  It does not import SymPy,
Singular, python-flint, the Macaulay builder, or the block-Wiedemann solver.
"""

from __future__ import annotations

import hashlib
import json
import sys
from fractions import Fraction
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CERTIFICATE = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "degree42_c6_macaulay_certificate.json"
)
# Filled after the canonical problem has been derived once from the committed
# Ritt residual construction.  The checker refuses a self-consistent but
# substituted polynomial system.
EXPECTED_PROBLEM_SHA256 = (
    "b5c115b3f2efb61c54017f251274cf7f537d650e5333a039949b4d5833fc4218"
)


def digest(value: object) -> str:
    return hashlib.sha256(
        json.dumps(
            value,
            sort_keys=True,
            separators=(",", ":"),
        ).encode()
    ).hexdigest()


def coefficient(value: object) -> Fraction:
    if isinstance(value, list):
        return Fraction(int(value[0]), int(value[1]))
    return Fraction(int(value))


def polynomial(data: list[list[object]]) -> dict[tuple[int, ...], Fraction]:
    return {
        tuple(int(exponent) for exponent in monomial): coefficient(value)
        for monomial, value in data
        if coefficient(value)
    }


def add(
    destination: dict[tuple[int, ...], Fraction],
    source: dict[tuple[int, ...], Fraction],
    scale: Fraction,
) -> None:
    for monomial, value in source.items():
        updated = destination.get(monomial, Fraction(0)) + scale * value
        if updated:
            destination[monomial] = updated
        else:
            destination.pop(monomial, None)


def monomial_multiple(
    source: dict[tuple[int, ...], Fraction],
    multiplier: tuple[int, ...],
    normal_cutoff: int | None,
) -> dict[tuple[int, ...], Fraction]:
    result = {}
    for monomial, value in source.items():
        product = tuple(
            left + right
            for left, right in zip(monomial, multiplier, strict=True)
        )
        if normal_cutoff is not None and sum(product[:2]) >= normal_cutoff:
            continue
        result[product] = value
    return result


def multiply_variable(
    source: dict[tuple[int, ...], Fraction],
    variable: int,
) -> dict[tuple[int, ...], Fraction]:
    result = {}
    for monomial, value in source.items():
        product = list(monomial)
        product[variable] += 1
        result[tuple(product)] = value
    return result


def replay_annihilator(
    generators: list[dict[tuple[int, ...], Fraction]],
    target: dict[tuple[int, ...], Fraction],
    entries: list[dict[str, Any]],
    normal_cutoff: int,
) -> None:
    residual = {monomial: -value for monomial, value in target.items()}
    for entry in entries:
        generator = generators[int(entry["generator"])]
        multiplier = tuple(int(value) for value in entry["monomial"])
        add(
            residual,
            monomial_multiple(generator, multiplier, normal_cutoff),
            coefficient(entry["coefficient"]),
        )
    assert not residual


def pairing(
    functional: dict[tuple[int, ...], Fraction],
    source: dict[tuple[int, ...], Fraction],
) -> Fraction:
    return sum(
        value * source.get(monomial, Fraction(0))
        for monomial, value in functional.items()
    )


def relevant_multipliers(
    support: tuple[tuple[int, ...], ...],
    generator: dict[tuple[int, ...], Fraction],
) -> set[tuple[int, ...]]:
    result = set()
    for support_monomial in support:
        for generator_monomial in generator:
            if all(
                left >= right
                for left, right in zip(
                    support_monomial,
                    generator_monomial,
                    strict=True,
                )
            ):
                result.add(
                    tuple(
                        left - right
                        for left, right in zip(
                            support_monomial,
                            generator_monomial,
                            strict=True,
                        )
                    )
                )
    return result


def main() -> None:
    path = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_CERTIFICATE
    data = json.loads(path.read_text())
    assert data["schema"] == "degree42-c6-macaulay-certificate.v1"
    claimed_certificate_digest = data.pop("certificate_sha256")
    assert digest(data) == claimed_certificate_digest
    problem = data["problem"]
    claimed_problem_digest = problem.pop("sha256")
    assert digest(problem) == claimed_problem_digest
    assert claimed_problem_digest == EXPECTED_PROBLEM_SHA256
    assert problem["variables"] == [
        "u",
        "v",
        "sync42p_w0",
        "sync42p_w1",
        "sync42p_w2",
    ]
    assert problem["normal_cutoff"] == 6
    assert data["macaulay"] == {
        "membership": {
            "rows": 646,
            "columns": 2400,
            "rank": 547,
            "normal_multiplier_degree": 4,
            "boundary_multiplier_degree": 2,
        },
        "dual": {
            "equations": 1249,
            "functional_coordinates": 210,
            "rank": 134,
            "boundary_support_degree": 2,
        },
    }
    for run_kind in ("membership_runs", "dual_runs"):
        runs = data["modular"][run_kind]
        assert len(runs) == 8
        assert len({run["prime"] for run in runs}) == 8
        assert all(int(run["prime"]).bit_length() == 31 for run in runs)
        assert all(
            run["algorithm"] == "two-sided-block-krylov-wiedemann"
            for run in runs
        )
    generators = [polynomial(item) for item in problem["generators"]]
    c6 = polynomial(problem["c6"])
    assert len(generators) == 16
    w0_c6 = multiply_variable(c6, 2)
    w2_c6 = multiply_variable(c6, 4)
    replay_annihilator(
        generators,
        w0_c6,
        data["annihilation_certificates"]["w0"],
        6,
    )
    replay_annihilator(
        generators,
        w2_c6,
        data["annihilation_certificates"]["w2"],
        6,
    )
    functional = polynomial(data["nonmembership_functional"])
    assert pairing(functional, c6) == 1
    support = tuple(functional)
    assert all(sum(monomial[:2]) < 6 for monomial in support)
    for generator in generators:
        for multiplier in relevant_multipliers(support, generator):
            product = monomial_multiple(generator, multiplier, None)
            assert pairing(functional, product) == 0
    # Extending the functional by zero outside its finite support therefore
    # kills every generator multiple.  It also kills (u,v)^6 by the support
    # condition above, while taking c6 to one.
    print("PASS: w0*c6 and w2*c6 vanish modulo the exact rational J6")
    print("PASS: finite Macaulay dual functional proves c6 is nonzero")
    print(
        "PASS: characteristic-zero degree-42 embedded-support "
        "obstruction certified"
    )


if __name__ == "__main__":
    main()

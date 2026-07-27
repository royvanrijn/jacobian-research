#!/usr/bin/env python3
"""Verify finite permutation actions of the Roe--Turturean G_Q2 presentation.

The checker is deliberately dependency-free.  A permutation is stored by its
zero-based image list, and multiplication is functional composition:

    (a * b)(i) = a(b(i)).

Thus the paper's conventions are evaluated literally:

    x^g = g^-1 * x * g,    [x,y] = x^-1 * y^-1 * x * y.

This verifies presentation admissibility and finite action invariants.  An
optional comparison block can add a narrowly supported polynomial-to-action
certificate; it is not inferred from orbit sizes.
"""

from __future__ import annotations

import argparse
import json
from collections import deque
from math import gcd, lcm
from pathlib import Path
from typing import Any


Permutation = tuple[int, ...]
SCHEMA = "gq2-permutation-action/v1"
CONVENTION = {
    "points": "zero-based",
    "permutation_encoding": "image-list",
    "multiplication": "(a*b)(i)=a(b(i))",
    "conjugation": "x^g=g^-1*x*g",
    "commutator": "[x,y]=x^-1*y^-1*x*y",
    "frobenius": "sigma maps to geometric Frobenius",
}


def identity(degree: int) -> Permutation:
    return tuple(range(degree))


def validate_permutation(values: Any, degree: int, name: str) -> Permutation:
    assert isinstance(values, list), f"{name}: expected an image list"
    result = tuple(values)
    assert len(result) == degree, f"{name}: wrong degree"
    assert all(isinstance(value, int) for value in result), (
        f"{name}: images must be integers"
    )
    assert set(result) == set(range(degree)), f"{name}: not a permutation"
    return result


def multiply(left: Permutation, right: Permutation) -> Permutation:
    """Return ``left * right`` (apply ``right`` first)."""

    assert len(left) == len(right)
    return tuple(left[right[index]] for index in range(len(left)))


def inverse(value: Permutation) -> Permutation:
    result = [0] * len(value)
    for source, target in enumerate(value):
        result[target] = source
    return tuple(result)


def power(value: Permutation, exponent: int) -> Permutation:
    assert exponent >= 0
    result = identity(len(value))
    base = value
    while exponent:
        if exponent & 1:
            result = multiply(result, base)
        base = multiply(base, base)
        exponent //= 2
    return result


def conjugate(value: Permutation, by: Permutation) -> Permutation:
    return multiply(multiply(inverse(by), value), by)


def commutator(left: Permutation, right: Permutation) -> Permutation:
    return multiply(
        multiply(multiply(inverse(left), inverse(right)), left),
        right,
    )


def cycle_lengths(value: Permutation, *, include_fixed: bool = True) -> list[int]:
    seen: set[int] = set()
    lengths: list[int] = []
    for start in range(len(value)):
        if start in seen:
            continue
        current = start
        length = 0
        while current not in seen:
            seen.add(current)
            current = value[current]
            length += 1
        if include_fixed or length != 1:
            lengths.append(length)
    return sorted(lengths, reverse=True)


def permutation_order(value: Permutation) -> int:
    result = 1
    for length in cycle_lengths(value):
        result = lcm(result, length)
    return result


def omega2_exponent(order: int) -> int:
    """CRT representative: 1 on the 2-part and 0 on the odd part."""

    assert order >= 1
    two_part = 1
    odd_part = order
    while odd_part % 2 == 0:
        two_part *= 2
        odd_part //= 2
    if two_part == 1:
        return 0
    return (odd_part * pow(odd_part, -1, two_part)) % order


def power_omega2(value: Permutation) -> tuple[Permutation, int, int]:
    order = permutation_order(value)
    exponent = omega2_exponent(order)
    return power(value, exponent), order, exponent


def generated_group(
    generators: list[Permutation],
    *,
    max_group_order: int,
) -> tuple[Permutation, ...]:
    assert generators
    degree = len(generators[0])
    one = identity(degree)
    steps = tuple(dict.fromkeys(generators + [inverse(g) for g in generators]))
    seen = {one}
    queue: deque[Permutation] = deque([one])
    while queue:
        current = queue.popleft()
        for step in steps:
            candidate = multiply(current, step)
            if candidate in seen:
                continue
            seen.add(candidate)
            assert len(seen) <= max_group_order, (
                "generated image exceeds --max-group-order"
            )
            queue.append(candidate)
    return tuple(sorted(seen))


def normal_closure(
    group: tuple[Permutation, ...],
    generators: list[Permutation],
    *,
    max_group_order: int,
) -> tuple[Permutation, ...]:
    conjugates = [
        conjugate(generator, element)
        for element in group
        for generator in generators
    ]
    return generated_group(conjugates, max_group_order=max_group_order)


def is_power_of_two(value: int) -> bool:
    return value > 0 and value & (value - 1) == 0


def action_orbits(
    group: tuple[Permutation, ...],
    degree: int,
) -> list[list[int]]:
    remaining = set(range(degree))
    result: list[list[int]] = []
    while remaining:
        representative = min(remaining)
        orbit = sorted({element[representative] for element in group})
        result.append(orbit)
        remaining.difference_update(orbit)
    return sorted(result, key=lambda orbit: (len(orbit), orbit), reverse=True)


def v_p(value: int, prime: int) -> int:
    assert value != 0
    value = abs(value)
    result = 0
    while value % prime == 0:
        result += 1
        value //= prime
    return result


def cubic_discriminant(coefficients: list[int]) -> int:
    """Discriminant for coefficients in low-to-high order."""

    assert len(coefficients) == 4
    d, c, b, a = coefficients
    return (
        b * b * c * c
        - 4 * a * c * c * c
        - 4 * b * b * b * d
        - 27 * a * a * d * d
        + 18 * a * b * c * d
    )


def verify_tame_eisenstein_cubic_comparison(
    comparison: dict[str, Any],
    generators: dict[str, Permutation],
    report: dict[str, Any],
) -> dict[str, Any]:
    """Check the supported Q_2 comparison for the tame S_3 cubic."""

    coefficients = comparison.get("polynomial_coefficients_low_to_high")
    assert isinstance(coefficients, list) and len(coefficients) == 4
    assert all(isinstance(value, int) for value in coefficients)
    prime = comparison.get("prime")
    assert prime == 2, "the v1 named comparison is specifically dyadic"
    assert coefficients[-1] % prime != 0
    assert all(value % prime == 0 for value in coefficients[:-1])
    assert coefficients[0] % (prime * prime) != 0

    discriminant = cubic_discriminant(coefficients)
    disc_v2 = v_p(discriminant, 2)
    disc_unit_mod_8 = (discriminant // (2**disc_v2)) % 8
    assert disc_v2 % 2 == 0 and disc_unit_mod_8 != 1, (
        "cubic discriminant must be nonsquare in Q_2"
    )

    one = identity(3)
    assert generators["x0"] == one and generators["x1"] == one
    assert cycle_lengths(generators["sigma"]) == [2, 1]
    assert cycle_lengths(generators["tau"]) == [3]
    assert report["image_group_order"] == 6
    assert report["orbit_sizes"] == [3]

    return {
        "type": comparison["type"],
        "prime": 2,
        "eisenstein": True,
        "discriminant": discriminant,
        "discriminant_v2": disc_v2,
        "discriminant_unit_mod_8": disc_unit_mod_8,
        "conclusion": (
            "irreducible tame cubic with S_3 splitting action; "
            "sigma is geometric Frobenius and tau is tame inertia"
        ),
    }


def verify_expected(actual: dict[str, Any], expected: dict[str, Any]) -> None:
    for key, value in expected.items():
        assert key in actual, f"expected field {key!r} is not reported"
        assert actual[key] == value, (
            f"expected {key}={value!r}, got {actual[key]!r}"
        )


def verify_certificate(
    certificate: dict[str, Any],
    *,
    max_group_order: int = 2_000_000,
) -> dict[str, Any]:
    assert certificate.get("schema") == SCHEMA
    assert certificate.get("conventions") == CONVENTION
    degree = certificate.get("degree")
    assert isinstance(degree, int) and degree >= 1

    raw_generators = certificate.get("generators")
    assert isinstance(raw_generators, dict)
    assert set(raw_generators) == {"sigma", "tau", "x0", "x1"}
    generators = {
        name: validate_permutation(raw_generators[name], degree, name)
        for name in ("sigma", "tau", "x0", "x1")
    }
    sigma = generators["sigma"]
    tau = generators["tau"]
    x0 = generators["x0"]
    x1 = generators["x1"]
    one = identity(degree)

    sigma2, sigma_order, sigma_omega_exponent = power_omega2(sigma)
    x0tau = multiply(x0, tau)
    x1tau = multiply(x1, tau)
    u0, x0tau_order, x0tau_omega_exponent = power_omega2(x0tau)
    u1, x1tau_order, x1tau_omega_exponent = power_omega2(x1tau)
    d0 = multiply(u0, inverse(x0))
    z0 = conjugate(x0, sigma2)
    c0 = commutator(d0, z0)
    g0 = power(sigma2, 2)
    dg = conjugate(d0, g0)
    hc = commutator(dg, d0)
    h0 = one
    for factor in (
        conjugate(x0, g0),
        x0,
        dg,
        d0,
        power(d0, 2),
        hc,
    ):
        h0 = multiply(h0, factor)

    tame_left = conjugate(tau, sigma)
    tame_right = power(tau, 2)
    tame_relator = multiply(tame_left, inverse(tame_right))
    wild_relator = one
    for factor in (h0, inverse(u1), conjugate(x1, sigma), c0):
        wild_relator = multiply(wild_relator, factor)

    image = generated_group(
        [sigma, tau, x0, x1],
        max_group_order=max_group_order,
    )
    wild = normal_closure(
        image,
        [x0, x1],
        max_group_order=max_group_order,
    )
    orbits = action_orbits(image, degree)
    orbit_sizes = sorted((len(orbit) for orbit in orbits), reverse=True)
    stabilizer_orders = sorted(
        (len(image) // len(orbit) for orbit in orbits),
        reverse=True,
    )

    report: dict[str, Any] = {
        "schema": SCHEMA,
        "degree": degree,
        "admissible": (
            tame_relator == one
            and wild_relator == one
            and is_power_of_two(len(wild))
        ),
        "tame_relation": tame_relator == one,
        "wild_relation": wild_relator == one,
        "wild_normal_closure_is_2_group": is_power_of_two(len(wild)),
        "image_group_order": len(image),
        "wild_normal_closure_order": len(wild),
        "orbits": orbits,
        "orbit_sizes": orbit_sizes,
        "stabilizer_orders": stabilizer_orders,
        "generator_cycle_types": {
            name: cycle_lengths(value)
            for name, value in generators.items()
        },
        "omega2_evaluations": {
            "sigma": {
                "order": sigma_order,
                "exponent": sigma_omega_exponent,
            },
            "x0*tau": {
                "order": x0tau_order,
                "exponent": x0tau_omega_exponent,
            },
            "x1*tau": {
                "order": x1tau_order,
                "exponent": x1tau_omega_exponent,
            },
        },
    }
    assert report["admissible"], "the marked permutation action is not admissible"

    comparison = certificate.get("comparison")
    if comparison is not None:
        assert isinstance(comparison, dict)
        comparison_type = comparison.get("type")
        if comparison_type == "tame-eisenstein-cubic-q2/v1":
            report["comparison"] = verify_tame_eisenstein_cubic_comparison(
                comparison,
                generators,
                report,
            )
        else:
            raise AssertionError(f"unsupported comparison type: {comparison_type!r}")

    expected = certificate.get("expected")
    if expected is not None:
        assert isinstance(expected, dict)
        verify_expected(report, expected)
    return report


def verify_certificate_path(
    path: Path,
    *,
    max_group_order: int = 2_000_000,
) -> dict[str, Any]:
    return verify_certificate(
        json.loads(path.read_text()),
        max_group_order=max_group_order,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("certificate", type=Path)
    parser.add_argument("--max-group-order", type=int, default=2_000_000)
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()
    report = verify_certificate_path(
        args.certificate,
        max_group_order=args.max_group_order,
    )
    if args.json:
        print(json.dumps(report, indent=2, sort_keys=True))
    else:
        print(
            "PASS G_Q2 action:",
            f"degree={report['degree']}",
            f"image_order={report['image_group_order']}",
            f"orbits={report['orbit_sizes']}",
            f"wild_core_order={report['wild_normal_closure_order']}",
        )
        if "comparison" in report:
            print("PASS polynomial comparison:", report["comparison"]["conclusion"])


if __name__ == "__main__":
    main()

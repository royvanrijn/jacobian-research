#!/usr/bin/env python3
"""Exact bounded construction fingerprints for ICARM curves 273 and 281--286.

The family-recognition screen is complete for the 2,329 normalized,
generically nonsingular six-root Mestre tuples of diameter at most 300, plus
the larger normalized Fermigier tuple.  It solves the exact j-equation in
z=T^2.  Modular no-root witnesses reject almost all pairs; the survivors are
factored over QQ and every reported match is checked by exact substitution.

This is a bounded family census, not a proof that an unmatched curve has no
Mestre, K3, Nagao, Kihara, isogenous, or private-family construction.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction
from math import gcd, isqrt
from pathlib import Path
from typing import Any, Iterable, Sequence


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
EC_GENERATED = GENERATED / "elliptic-curves"
DEFAULT_OUTPUT = EC_GENERATED / "icarm_construction_fingerprints_v1.json"
Q = Fraction

import sys

sys.path.insert(0, str(ROOT / "elliptic-curves/cas"))

from icarm_curve273 import GENERAL_WEIERSTRASS_COEFFICIENTS, POINTS  # noqa: E402
import icarm_curve302  # noqa: E402
from mestre_root_tuples import SixRootMestreConstruction  # noqa: E402


TARGET_SOURCE = EC_GENERATED / "icarm_7fff_zip_public_source_281_282_285_286.json"
CENSUS_SOURCES = (
    ROOT / "archive/elliptic-curves/artifacts/generated-results/elliptic_mestre_root_tuple_scale_max200_census.json",
    ROOT / "archive/elliptic-curves/artifacts/generated-results/elliptic_mestre_root_tuple_scale_max300_census.json",
)
FERMIGIER_NORMALIZED_ROOTS = (0, 29, 658, 722, 981, 1036)
MODULAR_PRIMES = (101, 103, 107, 109, 127, 131, 137, 139)


def text(value: Fraction | int) -> str:
    value = Q(value)
    return str(value.numerator) if value.denominator == 1 else str(value)


def weierstrass_invariants(
    coefficients: Sequence[Fraction | int],
) -> tuple[Fraction, Fraction, Fraction, Fraction]:
    a1, a2, a3, a4, a6 = (Q(value) for value in coefficients)
    b2 = a1 * a1 + 4 * a2
    b4 = 2 * a4 + a1 * a3
    b6 = a3 * a3 + 4 * a6
    b8 = a1 * a1 * a6 + 4 * a2 * a6 - a1 * a3 * a4 + a2 * a3 * a3 - a4 * a4
    c4 = b2 * b2 - 24 * b4
    c6 = -b2**3 + 36 * b2 * b4 - 216 * b6
    discriminant = -b2 * b2 * b8 - 8 * b4**3 - 27 * b6**2 + 9 * b2 * b4 * b6
    if discriminant == 0:
        raise ValueError("singular Weierstrass vector")
    return c4, c6, discriminant, c4**3 / discriminant


def load_targets(*, include_curve302: bool = False) -> dict[int, dict[str, Any]]:
    source = json.loads(TARGET_SOURCE.read_text())
    targets: dict[int, dict[str, Any]] = {}
    for record in source["curves"]:
        curve_id = int(record["id"])
        coefficients = tuple(Q(value) for value in record["ainvs"])
        targets[curve_id] = {
            "ainvs": coefficients,
            "points": tuple(
                (Q(point[0]), Q(point[1])) for point in record["points"]
            ),
        }
    targets[273] = {
        "ainvs": tuple(GENERAL_WEIERSTRASS_COEFFICIENTS),
        "points": tuple(POINTS),
    }
    if include_curve302:
        targets[302] = {
            "ainvs": tuple(icarm_curve302.GENERAL_WEIERSTRASS_COEFFICIENTS),
            "points": tuple(icarm_curve302.POINTS),
        }
    for target in targets.values():
        c4, c6, discriminant, j_value = weierstrass_invariants(target["ainvs"])
        target.update(c4=c4, c6=c6, discriminant=discriminant, j=j_value)
    return dict(sorted(targets.items()))


def load_census_roots() -> tuple[tuple[int, ...], ...]:
    max200 = json.loads(CENSUS_SOURCES[0].read_text())
    max300 = json.loads(CENSUS_SOURCES[1].read_text())
    roots = [
        tuple(values)
        for values in max200["tuple_populations"][
            "generically_nonsingular_nonreflection_roots"
        ]
    ]
    roots.extend(
        tuple(values)
        for values in max300["tuple_populations"]["genuinely_new_nonsingular_roots"]
    )
    if len(roots) != 2329 or len(set(roots)) != 2329:
        raise AssertionError("the complete diameter-at-most-300 census changed")
    if FERMIGIER_NORMALIZED_ROOTS in roots:
        raise AssertionError("the external Fermigier control unexpectedly entered the census")
    return tuple(roots) + (FERMIGIER_NORMALIZED_ROOTS,)


def multiply_mod(left: Sequence[int], right: Sequence[int], prime: int) -> list[int]:
    answer = [0] * (len(left) + len(right) - 1)
    for i, left_value in enumerate(left):
        for j, right_value in enumerate(right):
            answer[i + j] = (answer[i + j] + left_value * right_value) % prime
    return answer


def family_equation_mod(
    roots: Sequence[int], parameter: int, target_j: Fraction, prime: int
) -> int:
    parameter %= prime
    if parameter == 0 or target_j.denominator % prime == 0:
        raise ValueError("bad modular evaluation point")
    product = [1]
    for root in roots:
        product = multiply_mod(product, (-(root + parameter), 1), prime)
    for root in roots:
        product = multiply_mod(product, (-(root - parameter), 1), prime)
    approximant = [0] * 7
    approximant[6] = 1
    inverse_two = pow(2, -1, prime)
    for index in range(5, -1, -1):
        square = multiply_mod(approximant, approximant, prime)
        degree = 6 + index
        approximant[index] = (product[degree] - square[degree]) * inverse_two % prime
    square = multiply_mod(approximant, approximant, prime)
    inverse_t2 = pow(parameter * parameter % prime, -1, prime)
    quartic = [
        (square[index] - product[index]) * inverse_t2 % prime
        for index in range(5)
    ]
    e, d, c, b, a = quartic
    invariant_i = (12 * a * e - 3 * b * d + c * c) % prime
    invariant_j = (
        72 * a * c * e
        + 9 * b * c * d
        - 27 * a * d * d
        - 27 * b * b * e
        - 2 * c**3
    ) % prime
    a4 = -27 * invariant_i % prime
    a6 = -27 * invariant_j % prime
    c4 = -48 * a4 % prime
    discriminant = -16 * (4 * a4**3 + 27 * a6**2) % prime
    j_mod = target_j.numerator * pow(target_j.denominator, -1, prime) % prime
    return (c4**3 - j_mod * discriminant) % prime


def add_mod(left: Sequence[int], right: Sequence[int], prime: int) -> list[int]:
    answer = [0] * max(len(left), len(right))
    for index in range(len(answer)):
        answer[index] = (
            (left[index] if index < len(left) else 0)
            + (right[index] if index < len(right) else 0)
        ) % prime
    while len(answer) > 1 and answer[-1] == 0:
        answer.pop()
    return answer


def interpolate_mod(points: Sequence[tuple[int, int]], prime: int) -> list[int]:
    answer = [0]
    for index, (x_value, y_value) in enumerate(points):
        basis = [1]
        denominator = 1
        for other_index, (other_x, _) in enumerate(points):
            if index == other_index:
                continue
            basis = multiply_mod(basis, (-other_x, 1), prime)
            denominator = denominator * (x_value - other_x) % prime
        scale = y_value * pow(denominator, -1, prime) % prime
        answer = add_mod(answer, [scale * value % prime for value in basis], prime)
    return answer


def evaluate_mod(coefficients: Sequence[int], value: int, prime: int) -> int:
    answer = 0
    for coefficient in reversed(coefficients):
        answer = (answer * value + coefficient) % prime
    return answer


def modular_root_witness(
    roots: Sequence[int], target_j: Fraction, prime: int
) -> tuple[list[int], tuple[int, ...]] | None:
    if target_j.denominator % prime == 0:
        return None
    samples = []
    occupied: set[int] = set()
    for parameter in range(1, 14):
        z_value = parameter * parameter % prime
        if z_value in occupied:
            return None
        occupied.add(z_value)
        samples.append(
            (z_value, family_equation_mod(roots, parameter, target_j, prime))
        )
    polynomial = interpolate_mod(samples, prime)
    if len(polynomial) != 13 or polynomial[-1] == 0:
        return None
    for parameter in (14, 15, prime - 14):
        expected = family_equation_mod(roots, parameter, target_j, prime)
        if evaluate_mod(polynomial, parameter * parameter % prime, prime) != expected:
            raise AssertionError("the j-equation exceeded degree twelve in z=T^2")
    roots_mod_prime = tuple(
        value for value in range(prime) if evaluate_mod(polynomial, value, prime) == 0
    )
    return polynomial, roots_mod_prime


def multiply_fraction(
    left: Sequence[Fraction], right: Sequence[Fraction]
) -> list[Fraction]:
    answer = [Q(0)] * (len(left) + len(right) - 1)
    for i, left_value in enumerate(left):
        for j, right_value in enumerate(right):
            answer[i + j] += left_value * right_value
    return answer


def add_fraction(
    left: Sequence[Fraction], right: Sequence[Fraction]
) -> list[Fraction]:
    answer = [Q(0)] * max(len(left), len(right))
    for index in range(len(answer)):
        answer[index] = (
            (left[index] if index < len(left) else Q(0))
            + (right[index] if index < len(right) else Q(0))
        )
    return answer


def family_equation_exact(
    construction: SixRootMestreConstruction,
    parameter: Fraction,
    target_j: Fraction,
) -> Fraction:
    invariant_i, invariant_j = construction.binary_invariants(parameter)
    a4 = -27 * invariant_i
    a6 = -27 * invariant_j
    c4 = -48 * a4
    discriminant = -16 * (4 * a4**3 + 27 * a6**2)
    return c4**3 - target_j * discriminant


def interpolate_exact_j_equation(
    roots: Sequence[int], target_j: Fraction
) -> tuple[Fraction, ...]:
    construction = SixRootMestreConstruction(tuple(Q(value) for value in roots))
    points = [
        (Q(parameter * parameter), family_equation_exact(construction, Q(parameter), target_j))
        for parameter in range(1, 14)
    ]
    answer = [Q(0)]
    for index, (x_value, y_value) in enumerate(points):
        basis = [Q(1)]
        denominator = Q(1)
        for other_index, (other_x, _) in enumerate(points):
            if index == other_index:
                continue
            basis = multiply_fraction(basis, (-other_x, Q(1)))
            denominator *= x_value - other_x
        answer = add_fraction(answer, [y_value * value / denominator for value in basis])
    while len(answer) > 1 and answer[-1] == 0:
        answer.pop()
    if len(answer) != 13:
        raise AssertionError("the exact j-equation did not have degree twelve")
    for parameter in (Q(14), Q(15), Q(-14)):
        expected = family_equation_exact(construction, parameter, target_j)
        observed = sum(
            coefficient * (parameter * parameter) ** index
            for index, coefficient in enumerate(answer)
        )
        if observed != expected:
            raise AssertionError("exact j-equation interpolation failed")
    return tuple(answer)


def rational_square_root(value: Fraction) -> Fraction | None:
    if value <= 0:
        return None
    numerator = isqrt(value.numerator)
    denominator = isqrt(value.denominator)
    if numerator * numerator != value.numerator or denominator * denominator != value.denominator:
        return None
    return Q(numerator, denominator)


def exact_rational_square_roots(coefficients: Sequence[Fraction]) -> tuple[Fraction, ...]:
    try:
        import sympy as sp
    except ImportError as error:  # pragma: no cover - dependency is pinned
        raise RuntimeError("sympy is required for the surviving exact factorization") from error
    variable = sp.symbols("z")
    polynomial = sp.Poly(
        sum(
            sp.Rational(value.numerator, value.denominator) * variable**index
            for index, value in enumerate(coefficients)
        ),
        variable,
        domain=sp.QQ,
    )
    answer = []
    for root in sp.polys.polytools.ground_roots(polynomial):
        value = Q(int(root.p), int(root.q))
        square_root = rational_square_root(value)
        if square_root is not None:
            answer.append(square_root)
    return tuple(sorted(set(answer)))


def recognize_families(
    targets: dict[int, dict[str, Any]], roots_census: Sequence[tuple[int, ...]]
) -> dict[str, Any]:
    records = []
    for curve_id, target in targets.items():
        matches = []
        exact_factorization_survivor_count = 0
        survivors_with_rational_square_parameter = []
        witness_histogram: dict[str, int] = {}
        for roots in roots_census:
            excluded = False
            local_profile = []
            for prime in MODULAR_PRIMES:
                witness = modular_root_witness(roots, target["j"], prime)
                if witness is None:
                    continue
                _, roots_mod_prime = witness
                local_profile.append([prime, len(roots_mod_prime)])
                if not roots_mod_prime:
                    witness_histogram[str(prime)] = witness_histogram.get(str(prime), 0) + 1
                    excluded = True
                    break
            if excluded:
                continue
            coefficients = interpolate_exact_j_equation(roots, target["j"])
            parameters = exact_rational_square_roots(coefficients)
            exact_factorization_survivor_count += 1
            if parameters:
                survivors_with_rational_square_parameter.append(
                    {
                        "roots": list(roots),
                        "modular_root_counts": local_profile,
                        "rational_square_parameters": [text(value) for value in parameters],
                    }
                )
            construction = SixRootMestreConstruction(tuple(Q(value) for value in roots))
            for parameter in parameters:
                _, _, family_discriminant, family_j = weierstrass_invariants(
                    construction.primitive_jacobian_coefficients(parameter)
                )
                if family_discriminant and family_j == target["j"]:
                    matches.append({"roots": list(roots), "parameter_T": text(parameter)})
        records.append(
            {
                "curve_id": curve_id,
                "families_tested": len(roots_census),
                "modular_no_root_witness_histogram": witness_histogram,
                "exact_factorization_survivor_count": exact_factorization_survivor_count,
                "survivors_with_rational_square_parameter": survivors_with_rational_square_parameter,
                "exact_j_matches": matches,
            }
        )
    return {
        "scope": "complete normalized nonsingular six-root census through diameter 300 plus the larger Fermigier control tuple",
        "census_family_count": len(roots_census),
        "diameter_at_most_300_family_count": len(roots_census) - 1,
        "modular_primes": list(MODULAR_PRIMES),
        "targets": records,
        "boundary": "No-match results exclude only this fixed-root Mestre census. They do not exclude larger root tuples, generalized Mestre constructions, Kihara/Nagao families, K3 descendants, isogenous images, or private families.",
    }


def parse_vector(values: Sequence[Any]) -> tuple[Fraction, ...] | None:
    if len(values) != 5:
        return None
    answer = []
    for value in values:
        if isinstance(value, bool) or isinstance(value, (list, dict)):
            return None
        try:
            answer.append(Q(str(value)))
        except (ValueError, ZeroDivisionError):
            return None
    return tuple(answer)


def vectors_in_json(value: Any, pointer: str = "") -> Iterable[tuple[str, tuple[Fraction, ...]]]:
    if isinstance(value, list):
        vector = parse_vector(value)
        if vector is not None:
            try:
                weierstrass_invariants(vector)
            except ValueError:
                pass
            else:
                yield pointer or "/", vector
        for index, child in enumerate(value):
            yield from vectors_in_json(child, f"{pointer}/{index}")
    elif isinstance(value, dict):
        for key, child in value.items():
            yield from vectors_in_json(child, f"{pointer}/{key}")


def repository_j_scan(
    targets: dict[int, dict[str, Any]], *, excluded_paths: Sequence[Path] = ()
) -> dict[str, Any]:
    candidate_count = 0
    files_with_candidates = 0
    matches: dict[int, list[dict[str, str]]] = {curve_id: [] for curve_id in targets}
    for path in sorted(GENERATED.rglob("*.json")):
        if path in excluded_paths:
            continue
        try:
            payload = json.loads(path.read_text())
        except (OSError, json.JSONDecodeError):
            continue
        file_count = 0
        for pointer, vector in vectors_in_json(payload):
            file_count += 1
            candidate_count += 1
            _, _, _, j_value = weierstrass_invariants(vector)
            for curve_id, target in targets.items():
                if j_value == target["j"]:
                    matches[curve_id].append(
                        {
                            "path": str(path.relative_to(ROOT)),
                            "json_pointer": pointer,
                        }
                    )
        files_with_candidates += bool(file_count)
    return {
        "json_files_with_nonsingular_five_vectors": files_with_candidates,
        "nonsingular_five_vectors_tested": candidate_count,
        "matches_by_exact_j": {str(key): value for key, value in matches.items()},
        "interpretation": "Equal j is necessary but not sufficient for Q-isomorphism. Absence of an equal-j vector excludes a repository-stored Q-isomorphic copy; listed equal-j records still require twist/isomorphism checking.",
    }


def prime_factors(value: int) -> list[int]:
    value = abs(value)
    answer = []
    divisor = 2
    while divisor * divisor <= value:
        if value % divisor == 0:
            answer.append(divisor)
            while value % divisor == 0:
                value //= divisor
        divisor = 3 if divisor == 2 else divisor + 2
    if value > 1:
        answer.append(value)
    return answer


def denominator_fingerprint(points: Sequence[tuple[Fraction, Fraction]]) -> dict[str, Any]:
    square_roots = []
    for x_value, _ in points:
        denominator = x_value.denominator
        root = isqrt(denominator)
        if root * root != denominator:
            raise AssertionError("an x-coordinate denominator is not a square")
        square_roots.append(root)
    common = 0
    for value in square_roots:
        common = gcd(common, value)
    distinct = sorted(set(square_roots))
    return {
        "point_count": len(points),
        "integral_x_count": square_roots.count(1),
        "distinct_sqrt_x_denominator_count": len(distinct),
        "gcd_sqrt_x_denominators": common,
        "sqrt_x_denominators": distinct,
        "prime_support": sorted({prime for value in distinct for prime in prime_factors(value)}),
        "interpretation_boundary": "A reduced Mordell--Weil basis can obscure the raw construction sections; this fingerprint is diagnostic, not a family certificate.",
    }


def build_result(*, include_curve302: bool = False) -> dict[str, Any]:
    targets = load_targets(include_curve302=include_curve302)
    roots_census = load_census_roots()
    target_records = {
        str(curve_id): {
            "ainvs": [text(value) for value in target["ainvs"]],
            "j": text(target["j"]),
            "point_denominators": denominator_fingerprint(target["points"]),
            "rational_torsion": "trivial, by the separately pinned exact curve certificate",
        }
        for curve_id, target in targets.items()
    }
    result = {
        "schema": (
            "elliptic-curves.icarm-construction-fingerprints.v2"
            if include_curve302
            else "elliptic-curves.icarm-construction-fingerprints.v1"
        ),
        "status": "complete bounded exact construction-recognition screen",
        "inputs": {
            "target_source": str(TARGET_SOURCE.relative_to(ROOT)),
            "target_source_sha256": hashlib.sha256(TARGET_SOURCE.read_bytes()).hexdigest(),
            "census_sources": [str(path.relative_to(ROOT)) for path in CENSUS_SOURCES],
            "fermigier_control_roots": list(FERMIGIER_NORMALIZED_ROOTS),
        },
        "targets": target_records,
        "repository_model_scan": repository_j_scan(
            targets,
            excluded_paths=(
                (EC_GENERATED / "icarm_construction_fingerprints_v2.json",)
                if include_curve302
                else ()
            ),
        ),
        "six_root_mestre_recognition": recognize_families(targets, roots_census),
        "forced_torsion_exclusion": {
            "curves": sorted(targets),
            "result": "none is a direct specialization of the implemented Elkies--Klagsbrun K3 model y^2=x*(x^2+2A*x+B)",
            "reason": "that family has the rational 2-torsion point (0,0), while every target has certified trivial rational torsion",
            "boundary": "This does not exclude an isogenous quotient, a different K3 fibration, or another family on the same surface.",
        },
        "reproducing_command": (
            "python3 elliptic-curves/cas/analyze_icarm_construction_fingerprints.py --include-curve302"
            if include_curve302
            else "python3 elliptic-curves/cas/analyze_icarm_construction_fingerprints.py"
        ),
    }
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--include-curve302", action="store_true")
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()
    output = arguments.output
    if output is None:
        output = (
            EC_GENERATED / "icarm_construction_fingerprints_v2.json"
            if arguments.include_curve302
            else DEFAULT_OUTPUT
        )
    result = build_result(include_curve302=arguments.include_curve302)
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if arguments.check:
        if not output.exists() or output.read_text() != rendered:
            raise SystemExit(f"stale or missing artifact: {output}")
        print(f"PASS {output}")
        return
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(rendered)
    print(output)


if __name__ == "__main__":
    main()

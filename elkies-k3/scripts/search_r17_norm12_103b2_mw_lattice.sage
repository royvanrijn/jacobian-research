#!/usr/bin/env sage-python
"""Sieve sparse signed combinations in the ``0x103b2`` rank-17 subgroup."""

from __future__ import annotations

import argparse
import hashlib
from itertools import combinations, product
import json
from pathlib import Path
import runpy

from sage.all import PolynomialRing, QQ


ROOT = Path(__file__).resolve().parents[2]
HELPER_PATH = ROOT / "elkies-k3/scripts/certify_r17_norm12_103b2_jacobian.sage"
OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-norm12-103b2-mw-lattice-sieve-v1.json"
)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path):
    return str(path.resolve().relative_to(ROOT))


def reduce_rational(value, prime):
    value = QQ(value)
    return int(value.numerator()) % prime * pow(int(value.denominator()) % prime, -1, prime) % prime


def negate(point, a1, a3, prime):
    if point is None:
        return None
    x_value, y_value = point
    return x_value, (-y_value - a1 * x_value - a3) % prime


def add_points(left, right, coefficients, inverses, prime):
    if left is None:
        return right
    if right is None:
        return left
    a1, a2, a3, a4, a6 = coefficients
    x1, y1 = left
    x2, y2 = right
    if x1 != x2:
        denominator = (x2 - x1) % prime
        slope = (y2 - y1) * inverses[denominator] % prime
        intercept = (y1 * x2 - y2 * x1) * inverses[denominator] % prime
    else:
        if (y1 + y2 + a1 * x1 + a3) % prime == 0:
            return None
        denominator = (2 * y1 + a1 * x1 + a3) % prime
        if not denominator:
            return None
        slope = (
            3 * x1 * x1 + 2 * a2 * x1 + a4 - a1 * y1
        ) * inverses[denominator] % prime
        intercept = (
            -x1 * x1 * x1 + a4 * x1 + 2 * a6 - a3 * y1
        ) * inverses[denominator] % prime
    x3 = (slope * slope + a1 * slope - a2 - x1 - x2) % prime
    y3 = (-(slope + a1) * x3 - intercept - a3) % prime
    return x3, y3


def evaluate_quartic(coefficients, parameter, prime):
    value = 0
    for coefficient in reversed(coefficients):
        value = (value * parameter + coefficient) % prime
    return value


def build_prime_contexts(
    helper, records, pointed_curve, pointed_basis, constants, coefficient_radius,
    prime_count,
):
    integer_curves = [helper["normalized_integer_quartic"](record) for record in records]
    unused_a, unused_b, c_value, d_value, unused_e, v0 = constants
    rationals = [*pointed_curve.a_invariants(), c_value, d_value, v0]
    for point in pointed_basis:
        rationals.extend((point[0], point[1]))
    primes = []
    for candidate in helper["prime_range"](19, 20_000):
        prime = int(candidate)
        if any(int(value.denominator()) % prime == 0 for value in rationals):
            continue
        if int(pointed_curve.discriminant().numerator()) % prime == 0:
            continue
        if any(scalar % prime == 0 for scalar, unused in integer_curves):
            continue
        primes.append(prime)
        if len(primes) == prime_count:
            break
    if len(primes) != prime_count:
        raise ArithmeticError("not enough primes good for all seventeen basis points")
    target_index = next(
        index for index, record in enumerate(records)
        if record["label"] == helper["TARGET_LABEL"]
    )
    other_indices = [index for index in range(len(records)) if index != target_index]
    contexts = []
    for prime in primes:
        inverses = [0] + [pow(value, prime - 2, prime) for value in range(1, prime)]
        curve_coefficients = tuple(
            reduce_rational(value, prime) for value in pointed_curve.a_invariants()
        )
        a1, unused_a2, a3, unused_a4, unused_a6 = curve_coefficients
        positive_points = [
            (reduce_rational(point[0], prime), reduce_rational(point[1], prime))
            for point in pointed_basis
        ]
        multiple_points = []
        for positive_point in positive_points:
            values = [None] * (2 * coefficient_radius + 1)
            current = None
            for coefficient in range(1, coefficient_radius + 1):
                current = add_points(
                    current, positive_point, curve_coefficients, inverses, prime
                )
                values[coefficient_radius + coefficient] = current
                values[coefficient_radius - coefficient] = negate(
                    current, a1, a3, prime
                )
            multiple_points.append(values)
        field = helper["GF"](prime)
        reference_curve = helper["EllipticCurve"](field, curve_coefficients)
        for basis_index, positive_point in enumerate(positive_points):
            reference_point = reference_curve(*positive_point)
            for coefficient in range(-coefficient_radius, coefficient_radius + 1):
                if not coefficient:
                    continue
                expected = coefficient * reference_point
                observed = multiple_points[basis_index][coefficient_radius + coefficient]
                if expected.is_zero():
                    if observed is not None:
                        raise ArithmeticError("finite-field scalar multiplication mismatch")
                elif observed != (int(expected[0]), int(expected[1])):
                    raise ArithmeticError("finite-field scalar multiplication mismatch")
        for basis_index in range(len(positive_points)):
            next_index = (basis_index + 1) % len(positive_points)
            observed = add_points(
                positive_points[basis_index], positive_points[next_index],
                curve_coefficients, inverses, prime,
            )
            expected = reference_curve(*positive_points[basis_index]) + reference_curve(
                *positive_points[next_index]
            )
            if expected.is_zero():
                if observed is not None:
                    raise ArithmeticError("finite-field point addition mismatch")
            elif observed != (int(expected[0]), int(expected[1])):
                raise ArithmeticError("finite-field point addition mismatch")
        reduced_c = reduce_rational(c_value, prime)
        reduced_d = reduce_rational(d_value, prime)
        reduced_v0 = reduce_rational(v0, prime)
        base_parameter = reduce_rational(helper["KNOWN_PARAMETER"], prime)
        numerator_scale = 4 * reduced_v0 * reduced_v0 % prime
        denominator_scale = 2 * reduced_v0 % prime

        square_masks = []
        for parameter in range(prime):
            mask = 0
            for offset, record_index in enumerate(other_indices):
                scalar, coefficients = integer_curves[record_index]
                reduced_coefficients = tuple(value % prime for value in coefficients)
                value = scalar % prime * evaluate_quartic(
                    reduced_coefficients, parameter, prime
                ) % prime
                if value == 0 or pow(value, (prime - 1) // 2, prime) == 1:
                    mask |= 1 << offset
            square_masks.append(mask)
        contexts.append({
            "prime": prime,
            "inverses": inverses,
            "curve_coefficients": curve_coefficients,
            "coefficient_radius": coefficient_radius,
            "multiple_points": multiple_points,
            "c": reduced_c,
            "d": reduced_d,
            "v0": reduced_v0,
            "base_parameter": base_parameter,
            "numerator_scale": numerator_scale,
            "denominator_scale": denominator_scale,
            "square_masks": square_masks,
        })
    return contexts, other_indices


def parameter_mod_prime(point, context):
    if point is None:
        return None
    x_value, y_value = point
    prime = context["prime"]
    if not y_value:
        return None
    numerator = (
        context["numerator_scale"] * (x_value + context["c"])
        - context["d"] * context["d"]
    ) % prime
    denominator = context["denominator_scale"] * y_value % prime
    if not denominator:
        return None
    return (
        numerator * context["inverses"][denominator] + context["base_parameter"]
    ) % prime


def coefficient_vectors(rank, max_support, coefficient_radius):
    nonzero_coefficients = tuple(range(-coefficient_radius, 0)) + tuple(
        range(1, coefficient_radius + 1)
    )
    for support in range(1, max_support + 1):
        for indices in combinations(range(rank), support):
            for coefficients in product(nonzero_coefficients, repeat=support):
                yield indices, coefficients


def local_sieve(contexts, other_count, rank, max_support, coefficient_radius):
    full_mask = (1 << other_count) - 1
    candidate_count = 0
    survivors = []
    first_empty_histogram = [0] * len(contexts)
    for indices, coefficients in coefficient_vectors(
        rank, max_support, coefficient_radius
    ):
        candidate_count += 1
        mask = full_mask
        for context_index, context in enumerate(contexts):
            point = None
            for basis_index, coefficient in zip(indices, coefficients):
                point = add_points(
                    point,
                    context["multiple_points"][basis_index][
                        context["coefficient_radius"] + coefficient
                    ],
                    context["curve_coefficients"],
                    context["inverses"],
                    context["prime"],
                )
            parameter = parameter_mod_prime(point, context)
            if parameter is not None:
                mask &= context["square_masks"][parameter]
            if not mask:
                first_empty_histogram[context_index] += 1
                break
        if mask:
            survivors.append((indices, coefficients, mask))
    return candidate_count, first_empty_histogram, survivors


def exact_replay(helper, records, other_indices, primitive, pointed_curve, pointed_basis,
                 constants, survivors):
    ring = primitive.parent()
    results = []
    for indices, coefficients, mask in survivors:
        point = pointed_curve(0)
        vector = [0] * len(pointed_basis)
        for basis_index, coefficient in zip(indices, coefficients):
            point += coefficient * pointed_basis[basis_index]
            vector[basis_index] = coefficient
        if point.is_zero() or not point[1]:
            results.append({
                "coefficient_vector": vector,
                "status": "EXCEPTIONAL_POINT_FOR_INVERSE_MAP",
            })
            continue
        parameter, cover_coordinate = helper["inverse_parameter"](
            point, constants, helper["KNOWN_PARAMETER"]
        )
        if cover_coordinate**2 != primitive(parameter):
            raise ArithmeticError("exact inverse map failed")
        locally_surviving_labels = []
        exact_split_labels = []
        for offset, record_index in enumerate(other_indices):
            if not (mask >> offset) & 1:
                continue
            label = records[record_index]["label"]
            locally_surviving_labels.append(label)
            branch = helper["polynomial"](
                records[record_index]["branch_polynomial_q_coefficients_low_to_high"],
                ring,
            )
            if helper["rational_square_root"](branch(parameter)) is not None:
                exact_split_labels.append(label)
        results.append({
            "coefficient_vector": vector,
            "t": helper["rational_text"](parameter),
            "s_on_primitive_quartic": helper["rational_text"](cover_coordinate),
            "locally_surviving_cover_labels": locally_surviving_labels,
            "exact_split_cover_labels": exact_split_labels,
            "status": "EXACT_REPLAYED",
        })
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-support", type=int, default=5)
    parser.add_argument("--coefficient-radius", type=int, default=2)
    parser.add_argument("--prime-count", type=int, default=32)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if (
        not 1 <= args.max_support <= 17
        or not 1 <= args.coefficient_radius <= 8
        or args.prime_count < 8
    ):
        parser.error(
            "max support must be 1..17, coefficient radius 1..8, and at least eight primes are required"
        )

    helper = runpy.run_path(str(HELPER_PATH), run_name="r17_jacobian_helpers")
    splitting = json.loads(helper["SPLITTING"].read_text())
    records = splitting["construction"]["records"]
    target_index = next(
        index for index, record in enumerate(records)
        if record["label"] == helper["TARGET_LABEL"]
    )
    ring = PolynomialRing(QQ, "t")
    published = helper["polynomial"](
        records[target_index]["branch_polynomial_q_coefficients_low_to_high"], ring
    )
    primitive, unused_scale = helper["primitive_integral_quartic"](published)
    pointed_curve, unused_generator, constants = helper["pointed_jacobian"](
        primitive, helper["KNOWN_PARAMETER"], QQ(80653002864) / 625
    )
    minimal = pointed_curve.global_minimal_model()
    basis, point_search = helper["discover_rank_basis"](
        primitive, pointed_curve, constants
    )
    inverse_isomorphism = minimal.isomorphism_to(pointed_curve)
    pointed_basis = [inverse_isomorphism(point) for point in basis]
    contexts, other_indices = build_prime_contexts(
        helper, records, pointed_curve, pointed_basis, constants,
        args.coefficient_radius, args.prime_count,
    )
    candidate_count, first_empty_histogram, survivors = local_sieve(
        contexts, len(other_indices), len(pointed_basis), args.max_support,
        args.coefficient_radius,
    )
    exact_results = exact_replay(
        helper, records, other_indices, primitive, pointed_curve, pointed_basis,
        constants, survivors,
    )
    exact_hits = [
        result for result in exact_results if result.get("exact_split_cover_labels")
    ]
    result = {
        "schema": "elkies-k3.r17-norm12-103b2-mw-lattice-sieve.v1",
        "status": "PASS_BOUNDED_SPARSE_SIGNED_LATTICE_SIEVE",
        "inputs": {
            relative(helper["SPLITTING"]): digest(helper["SPLITTING"]),
            relative(HELPER_PATH): digest(HELPER_PATH),
        },
        "source_cover": helper["TARGET_LABEL"],
        "primitive_integral_quartic_coefficients_low_to_high": [
            str(coefficient) for coefficient in primitive
        ],
        "jacobian": {
            "global_minimal_model_a_invariants": [
                str(value) for value in minimal.a_invariants()
            ],
            "rational_torsion_order": int(minimal.torsion_subgroup().order()),
            "known_subgroup_rank": len(basis),
            "rank_claim": "at least 17; no rank upper bound is used",
            "independent_generators_on_global_minimal_model": [
                [str(point[0]), str(point[1])] for point in basis
            ],
        },
        "point_search": point_search,
        "coefficient_population": {
            "coefficient_radius": args.coefficient_radius,
            "allowed_coefficients": list(range(-args.coefficient_radius, args.coefficient_radius + 1)),
            "maximum_support": args.max_support,
            "candidate_count": candidate_count,
        },
        "other_cover_count": len(other_indices),
        "sieve_primes": [context["prime"] for context in contexts],
        "first_empty_prime_histogram": {
            str(contexts[index]["prime"]): count
            for index, count in enumerate(first_empty_histogram) if count
        },
        "local_survivor_count": len(survivors),
        "exact_results": exact_results,
        "exact_simultaneous_splits": exact_hits,
        "proof_boundary": (
            f"The sieve is exhaustive only for coefficient vectors in [-{args.coefficient_radius},"
            f"{args.coefficient_radius}]^17 of support at most {args.max_support}. "
            "A nonsquare reduction proves "
            "a rational nonsquare; every all-prime survivor is replayed exactly."
        ),
        "reproducing_command": (
            "sage -python elkies-k3/scripts/search_r17_norm12_103b2_mw_lattice.sage "
            f"--max-support {args.max_support} --coefficient-radius {args.coefficient_radius} "
            f"--prime-count {args.prime_count}"
        ),
    }
    serialized = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.check:
        if args.output.read_text() != serialized:
            raise ArithmeticError("stored lattice-sieve certificate differs from replay")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(serialized)
    print(
        f"R17103B2LATTICE|candidates={candidate_count}|local_survivors={len(survivors)}|"
        f"exact_hits={len(exact_hits)}|output={args.output}"
    )


if __name__ == "__main__":
    main()

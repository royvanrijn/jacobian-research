#!/usr/bin/env sage-python
"""Compute the Jacobian of the pointed norm-12 ``0x103b2`` quartic.

The exact rank computation is separated from the targeted parameter sweep:

* PARI ``ellrank`` supplies rigorous lower and upper Mordell--Weil bounds;
* eclib saturation proves that the visible rank-one point is primitive;
* the inverse pointed-quartic map defines exact rational parameters ``t_n``
  from multiples ``n*G``;
* reductions of those exact parameters give a rational nonsquare witness for
  every other compiled cover.  This avoids expanding enormous rational
  numerators while retaining an exact, replayable exclusion certificate.
"""

from __future__ import annotations

import argparse
from array import array
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import sys

from sage.all import EllipticCurve, GF, PolynomialRing, QQ, ZZ, gcd, lcm, pari, prime_range


ROOT = Path(__file__).resolve().parents[2]
SPLITTING = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-genus-one-bisection-splitting-search-v1.json"
)
OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-norm12-103b2-jacobian-v1.json"
)
TARGET_LABEL = "norm12-orbit-103b2"
KNOWN_PARAMETER = QQ(1) / 25
KNOWN_COVER_COORDINATE = QQ("3521934804796232704/643125")


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def rational_text(value) -> str:
    return str(QQ(value))


def polynomial(coefficients, ring):
    return ring([QQ(value) for value in coefficients])


def rational_square_root(value):
    value = QQ(value)
    if value < 0:
        return None
    numerator = math.isqrt(int(value.numerator()))
    denominator = math.isqrt(int(value.denominator()))
    if numerator * numerator != value.numerator():
        return None
    if denominator * denominator != value.denominator():
        return None
    return QQ(numerator) / denominator


def primitive_integral_quartic(q):
    denominator = lcm([coefficient.denominator() for coefficient in q])
    integral = [ZZ(coefficient * denominator) for coefficient in q]
    content = gcd(integral)
    primitive = q.parent()([coefficient // content for coefficient in integral])
    if primitive.leading_coefficient() < 0:
        primitive = -primitive
        content = -content
    scale_square = QQ(content) / denominator
    scale = rational_square_root(scale_square)
    if scale is None:
        raise ArithmeticError("the target quartic clearing scale is not a rational square")
    if q != scale**2 * primitive:
        raise ArithmeticError("primitive quartic normalization failed")
    return primitive, scale


def pointed_jacobian(quartic, parameter, cover_coordinate):
    local_ring = PolynomialRing(QQ, "z")
    z = local_ring.gen()
    shifted = local_ring(quartic(z + parameter))
    e, d, c, b, a = (QQ(shifted[index]) for index in range(5))
    v0 = QQ(cover_coordinate)
    if not v0 or v0**2 != e:
        raise ArithmeticError("invalid pointed-quartic base point")
    curve = EllipticCurve(QQ, [
        d / v0,
        c - d**2 / (4 * v0**2),
        2 * v0 * b,
        -4 * v0**2 * a,
        a * (d**2 - 4 * v0**2 * c),
    ])
    opposite_x = d**2 / (4 * v0**2) - c
    candidates = [point for point in curve.lift_x(opposite_x, all=True) if point[1]]
    if len(candidates) != 1 or candidates[0].has_finite_order():
        raise ArithmeticError("the pointed quartic did not supply one visible infinite-order point")
    return curve, candidates[0], (a, b, c, d, e, v0)


def inverse_parameter(point, constants, parameter):
    unused_a, unused_b, c, d, unused_e, v0 = constants
    x_value, y_value = point[:2]
    local_parameter = (4 * v0**2 * (x_value + c) - d**2) / (2 * v0 * y_value)
    cover_coordinate = (x_value * local_parameter**2 - d * local_parameter) / (2 * v0) - v0
    return QQ(local_parameter + parameter), QQ(cover_coordinate)


def certify_inverse_map_identity(quartic, curve, constants):
    unused_a, unused_b, c, d, unused_e, v0 = constants
    affine = PolynomialRing(QQ, names=("X", "Y"))
    X, Y = affine.gens()
    field = affine.fraction_field()
    a1, a2, a3, a4, a6 = curve.a_invariants()
    equation = Y**2 + a1 * X * Y + a3 * Y - X**3 - a2 * X**2 - a4 * X - a6
    local_parameter = field((4 * v0**2 * (X + c) - d**2) / (2 * v0 * Y))
    cover_coordinate = field(
        (X * local_parameter**2 - d * local_parameter) / (2 * v0) - v0
    )
    local_ring = PolynomialRing(QQ, "z")
    z = local_ring.gen()
    shifted = local_ring(quartic(z + KNOWN_PARAMETER))
    difference = field(cover_coordinate**2 - shifted(local_parameter))
    numerator = affine(difference.numerator())
    quotient, remainder = numerator.quo_rem(equation)
    if remainder:
        raise ArithmeticError("the inverse pointed-quartic map identity failed")
    return {
        "status": "EXACT_POLYNOMIAL_IDENTITY",
        "numerator_total_degree": int(numerator.total_degree()),
        "elliptic_equation_total_degree": int(equation.total_degree()),
        "quotient_sha256": hashlib.sha256(str(quotient).encode()).hexdigest(),
    }


def normalized_integer_quartic(record):
    values = [Fraction(value) for value in record["branch_polynomial_q_coefficients_low_to_high"]]
    denominator = 1
    for value in values:
        denominator = math.lcm(denominator, value.denominator)
    coefficients = tuple(
        value.numerator * (denominator // value.denominator) for value in values
    )
    return denominator, coefficients


def mod_rational(value, field, prime):
    value = QQ(value)
    denominator = int(value.denominator()) % prime
    if not denominator:
        raise ZeroDivisionError
    return field(int(value.numerator()) % prime) / field(denominator)


def select_good_primes(curves, curve, generator, constants, count):
    unused_a, unused_b, c, d, unused_e, v0 = constants
    rationals = [*curve.a_invariants(), generator[0], generator[1], c, d, v0]
    selected = []
    for prime in prime_range(19, 20_000):
        prime = int(prime)
        if any(int(value.denominator()) % prime == 0 for value in rationals):
            continue
        if int(curve.discriminant().numerator()) % prime == 0:
            continue
        if any(scalar % prime == 0 for scalar, unused in curves):
            continue
        selected.append(prime)
        if len(selected) == count:
            return selected
    raise ArithmeticError("not enough uniformly good witness primes")


def targeted_modular_sweep(records, source_index, curve, generator, constants, sample_count, prime_count):
    curves = [normalized_integer_quartic(record) for record in records]
    other_indices = [index for index in range(len(records)) if index != source_index]
    other_labels = [records[index]["label"] for index in other_indices]
    primes = select_good_primes(curves, curve, generator, constants, prime_count)
    witnesses = [[0] * len(other_indices) for unused in range(sample_count)]
    remaining = sample_count * len(other_indices)
    witness_histogram = {prime: 0 for prime in primes}
    unused_a, unused_b, c, d, unused_e, v0 = constants

    for prime in primes:
        field = GF(prime)
        reduced_curve = EllipticCurve(field, [
            mod_rational(value, field, prime) for value in curve.a_invariants()
        ])
        reduced_generator = reduced_curve(
            mod_rational(generator[0], field, prime),
            mod_rational(generator[1], field, prime),
        )
        reduced_c = mod_rational(c, field, prime)
        reduced_d = mod_rational(d, field, prime)
        reduced_v0 = mod_rational(v0, field, prime)
        reduced_parameter = mod_rational(KNOWN_PARAMETER, field, prime)
        reduced_curves = [
            (
                field(scalar % prime),
                tuple(field(coefficient % prime) for coefficient in coefficients),
            )
            for scalar, coefficients in curves
        ]

        current = reduced_generator
        for sample_offset in range(sample_count):
            current += reduced_generator
            if current.is_zero() or not current[1]:
                continue
            local_parameter = (
                4 * reduced_v0**2 * (current[0] + reduced_c) - reduced_d**2
            ) / (2 * reduced_v0 * current[1])
            target_parameter = local_parameter + reduced_parameter
            row = witnesses[sample_offset]
            for other_offset, record_index in enumerate(other_indices):
                if row[other_offset]:
                    continue
                scalar, coefficients = reduced_curves[record_index]
                value = scalar * sum(
                    coefficients[index] * target_parameter**index for index in range(5)
                )
                if value and not value.is_square():
                    row[other_offset] = prime
                    witness_histogram[prime] += 1
                    remaining -= 1
        if not remaining:
            break

    unresolved = [
        {
            "multiple": sample_offset + 2,
            "cover_labels": [
                other_labels[index] for index, witness in enumerate(row) if not witness
            ],
        }
        for sample_offset, row in enumerate(witnesses)
        if 0 in row
    ]
    if unresolved:
        raise ArithmeticError(
            f"{sum(len(row['cover_labels']) for row in unresolved)} pairs lack a local nonsquare witness: "
            f"{unresolved[:5]}"
        )

    matrix_hash = hashlib.sha256()
    for row in witnesses:
        packed = array("H", row)
        if sys.byteorder != "little":
            packed.byteswap()
        matrix_hash.update(packed.tobytes())
    by_cover = {
        label: {
            str(prime): sum(row[offset] == prime for row in witnesses)
            for prime in primes
            if any(row[offset] == prime for row in witnesses)
        }
        for offset, label in enumerate(other_labels)
    }
    return {
        "multiple_range_inclusive": [2, sample_count + 1],
        "exact_rational_parameter_count": sample_count,
        "other_cover_count": len(other_indices),
        "parameter_cover_pairs": sample_count * len(other_indices),
        "uniformly_good_witness_primes": primes,
        "witness_prime_histogram": {
            str(prime): count for prime, count in witness_histogram.items() if count
        },
        "witness_histogram_by_cover": by_cover,
        "first_witness_matrix": {
            "row_order": "multiples 2 through sample_count+1",
            "column_labels": other_labels,
            "entry_format": "little-endian uint16 first nonsquare-witness prime",
            "sha256": matrix_hash.hexdigest(),
        },
        "unresolved_pair_count": 0,
        "simultaneous_split_count": 0,
        "conclusion": (
            "For every listed n and every other compiled cover, the stored first witness prime "
            "has nonsquare reduction. Hence no other cover splits at t_n over QQ."
        ),
    }


def explicit_prefix(quartic, curve, generator, constants, count):
    rows = []
    manifest = hashlib.sha256()
    current = generator
    for multiple in range(2, count + 2):
        current += generator
        parameter, cover_coordinate = inverse_parameter(current, constants, KNOWN_PARAMETER)
        if cover_coordinate**2 != quartic(parameter):
            raise ArithmeticError(f"exact inverse-map check failed at multiple {multiple}")
        line = f"{multiple}\t{parameter}\t{cover_coordinate}\n"
        manifest.update(line.encode())
        rows.append({
            "multiple": multiple,
            "t": rational_text(parameter),
            "s_on_primitive_quartic": rational_text(cover_coordinate),
            "t_numerator_bits": int(abs(parameter.numerator()).nbits()),
            "t_denominator_bits": int(parameter.denominator().nbits()),
        })
    return rows, manifest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-count", type=int, default=10_000)
    parser.add_argument("--explicit-prefix", type=int, default=8)
    parser.add_argument("--prime-count", type=int, default=48)
    parser.add_argument("--pari-stack-gb", type=int, default=2)
    parser.add_argument("--skip-rank-descent", action="store_true")
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.sample_count < 100 or args.explicit_prefix < 1 or args.prime_count < 16:
        parser.error("use at least 100 samples, one explicit point, and 16 witness primes")

    splitting = json.loads(SPLITTING.read_text())
    records = splitting["construction"]["records"]
    by_label = {record["label"]: index for index, record in enumerate(records)}
    source_index = by_label[TARGET_LABEL]
    target = records[source_index]
    ring = PolynomialRing(QQ, "t")
    q = polynomial(target["branch_polynomial_q_coefficients_low_to_high"], ring)
    primitive, scale = primitive_integral_quartic(q)
    if scale * (QQ(80653002864) / 625) != KNOWN_COVER_COORDINATE:
        raise ArithmeticError("the normalized known point changed")

    curve, generator, constants = pointed_jacobian(
        primitive, KNOWN_PARAMETER, QQ(80653002864) / 625
    )
    minimal = curve.global_minimal_model()
    isomorphism = curve.isomorphism_to(minimal)
    minimal_generator = isomorphism(generator)
    if minimal_generator.height() <= 0:
        raise ArithmeticError("the visible Jacobian point is not nontorsion")
    torsion = minimal.torsion_subgroup()
    if torsion.order() != 1:
        raise ArithmeticError("unexpected rational torsion")
    saturated_points, saturation_index, regulator = minimal.saturation(
        [minimal_generator], max_prime=370
    )
    if saturation_index != 1 or saturated_points != [minimal_generator]:
        raise ArithmeticError("the visible point is not primitive in its rank-one span")

    if args.skip_rank_descent:
        if not args.output.exists():
            raise FileNotFoundError("--skip-rank-descent requires an existing certificate")
        stored = json.loads(args.output.read_text())
        rank_bounds = stored["jacobian"]["mordell_weil_rank_bounds"]
        descent_record = stored["jacobian"]["pari_2_descent"]
    else:
        pari.allocatemem(args.pari_stack_gb * 1024**3)
        descent = pari(minimal).ellrank(
            0,
            pari([[
                minimal_generator[0],
                minimal_generator[1],
            ]]),
        )
        rank_bounds = [int(descent[0]), int(descent[1])]
        descent_record = {
            "algorithm": "PARI ellrank with the visible point supplied",
            "raw_result": str(descent),
            "lower_bound": rank_bounds[0],
            "upper_bound": rank_bounds[1],
            "sha_2_information": int(descent[2]),
        }
    if rank_bounds != [1, 1]:
        raise ArithmeticError(f"the exact Jacobian rank is not one: {rank_bounds}")

    map_identity = certify_inverse_map_identity(primitive, curve, constants)
    prefix, prefix_hash = explicit_prefix(
        primitive, curve, generator, constants, args.explicit_prefix
    )
    sweep = targeted_modular_sweep(
        records, source_index, curve, generator, constants,
        args.sample_count, args.prime_count,
    )

    a, b, c, d, e, v0 = constants
    result = {
        "schema": "elkies-k3.r17-norm12-103b2-jacobian.v1",
        "status": "PASS_EXACT_RANK_ONE_AND_TARGETED_NO_OTHER_SPLITS",
        "inputs": {relative(SPLITTING): digest(SPLITTING)},
        "source_cover": {
            "label": TARGET_LABEL,
            "published_branch_coefficients_low_to_high": [str(value) for value in q],
            "primitive_integral_quartic_coefficients_low_to_high": [
                str(value) for value in primitive
            ],
            "published_to_primitive_cover_coordinate_scale": rational_text(scale),
            "known_point_published": {
                "t": rational_text(KNOWN_PARAMETER),
                "s": rational_text(KNOWN_COVER_COORDINATE),
            },
            "known_point_primitive": {
                "t": rational_text(KNOWN_PARAMETER),
                "s": rational_text(v0),
            },
        },
        "jacobian": {
            "pointed_model_a_invariants": [rational_text(value) for value in curve.a_invariants()],
            "global_minimal_model_a_invariants": [
                rational_text(value) for value in minimal.a_invariants()
            ],
            "pointed_to_minimal_isomorphism_u_r_s_t": [
                rational_text(value) for value in isomorphism.tuple()
            ],
            "generator_on_pointed_model": [
                rational_text(generator[0]), rational_text(generator[1])
            ],
            "generator_on_global_minimal_model": [
                rational_text(minimal_generator[0]), rational_text(minimal_generator[1])
            ],
            "generator_canonical_height": str(minimal_generator.height(precision=256)),
            "torsion_order": int(torsion.order()),
            "mordell_weil_rank_bounds": rank_bounds,
            "mordell_weil_group": "Z generated by the displayed point",
            "pari_2_descent": descent_record,
            "rank_one_saturation": {
                "index": int(saturation_index),
                "index_bound": 370,
                "regulator_approx": str(regulator),
            },
            "root_number": int(minimal.root_number()),
            "conductor": str(minimal.conductor()),
        },
        "inverse_pointed_quartic_map": {
            "shifted_coefficients_a_b_c_d_e": [rational_text(value) for value in (a, b, c, d, e)],
            "formula": {
                "z": "(4*v0^2*(X+c)-d^2)/(2*v0*Y)",
                "t": "z+1/25",
                "s": "(X*z^2-d*z)/(2*v0)-v0",
            },
            "identity_certificate": map_identity,
            "exception_boundary": (
                "The denominator Y vanishes at -G, not at nG for n>=2, because G is nontorsion."
            ),
        },
        "targeted_parameter_sweep": {
            **sweep,
            "parameter_definition": "t_n=t(n*G) under the displayed exact inverse map",
            "distinctness": (
                "The degree-two t-map pairs nG with (1-n)G. For distinct n,m>=2, "
                "neither n=m nor n=1-m, so the listed parameters are distinct."
            ),
            "explicit_prefix": prefix,
            "explicit_prefix_manifest_sha256": prefix_hash,
            "storage_boundary": (
                "Only the prefix is expanded: numerator and denominator heights grow quadratically in n. "
                "The full exact sequence is stored by generator, map, and multiple interval; every cross-cover "
                "exclusion has an independently replayed finite-field nonsquare witness."
            ),
        },
        "proof_boundary": (
            "PARI's equal 2-descent bounds prove rank one; torsion and saturation make the displayed point a generator. "
            "The targeted sweep is complete only for the declared multiple interval and the 142 other compiled covers."
        ),
        "reproducing_command": (
            "sage -python elkies-k3/scripts/certify_r17_norm12_103b2_jacobian.sage"
        ),
    }
    serialized = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.check:
        if args.output.read_text() != serialized:
            raise ArithmeticError("stored Jacobian certificate differs from replay")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(serialized)
    print(
        f"R17103B2JAC|rank={rank_bounds[0]}|torsion={torsion.order()}|"
        f"samples={args.sample_count}|other_covers={sweep['other_cover_count']}|"
        f"unresolved={sweep['unresolved_pair_count']}|output={args.output}"
    )


if __name__ == "__main__":
    main()

#!/usr/bin/env sage-python
"""Audit alternate-Q80 visibility in seven high-rank public fibres.

For the two 11952 fibres, three 08f72 fibres, and two 0e80b/103b2 fibres,
this replay specializes an exact saturated generic MW17 basis, proves the
displayed public points independent by finite reductions, recovers candidate
integral coordinates from a high-precision height solve, and then verifies
every relation by exact elliptic-curve group law.  It also exhausts every
constructed alternate-Q80 bisection cover at each fibre and computes the
exact Smith-coordinate span of every rational split.

Numerical heights are used only to propose integer relations.  Every retained
relation and every quotient calculation is exact and fail-closed.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys

from sage.all import EllipticCurve, Matrix, PolynomialRing, QQ, ZZ, matrix
from sage.env import SAGE_VERSION


sys.set_int_max_str_digits(0)
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "elliptic-curves/cas"))

from mod2_reduction_independence import (  # noqa: E402
    combined_mod2_rank,
    find_mod2_reduction_certificate,
    find_two_torsion_certificate_prime,
)


SWEEP = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-icarm-database-sweep-v1.json"
PUBLIC = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-icarm-public-fibres-v1.json"
LOCAL = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-icarm-local-fingerprints-v1.json"
CURVE12 = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-curve12-alternate-q80-quotient-v1.json"
OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-native-icarm-quotient-audit-v1.json"
PRECISION = 180
CURVE12_PREFERRED_BASIS = [
    "P2", "P11", "P4", "P3", "P6", "P8", "P17", "P10", "P28", "P24", "P19", "P15"
]

CHARTS = (
    {
        "source_chart": "norm12-orbit-11952",
        "representative": "norm12-orbit-11952",
        "curve_ids": (12, 395),
        "direct": ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json",
        "covers": ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-11952-alternate-bisections-full-v1.json",
    },
    {
        "source_chart": "norm12-orbit-08f72",
        "representative": "norm12-orbit-08f72",
        "curve_ids": (363, 364, 378),
        "direct": ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-orbit08f72-direct-fibration-v1.json",
        "covers": ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-08f72-alternate-bisections-full-v1.json",
    },
    {
        "source_chart": "norm12-orbit-103b2",
        "representative": "norm12-orbit-0e80b",
        "curve_ids": (393, 404),
        "direct": ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-orbit103b2-direct-fibration-v1.json",
        "covers": ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-103b2-bisections-full-v1.json",
    },
)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def rational_text(value) -> str:
    value = QQ(value)
    if value.denominator() == 1:
        return str(value.numerator())
    return f"{value.numerator()}/{value.denominator()}"


def polynomial(ring, values):
    return ring([QQ(value) for value in values])


def rational_function(ring, record):
    return ring.fraction_field()(
        polynomial(ring, record["numerator_coefficients_low_to_high"])
        / polynomial(ring, record["denominator_coefficients_low_to_high"])
    )


def short_invariants(ainvs):
    a1, a2, a3, a4, a6 = map(QQ, ainvs)
    b2 = a1**2 + 4 * a2
    b4 = a1 * a3 + 2 * a4
    b6 = a3**2 + 4 * a6
    c4 = b2**2 - 24 * b4
    c6 = -b2**3 + 36 * b2 * b4 - 216 * b6
    return a1, a3, b2, -c4 / 48, -c6 / 864


def python_fraction(value):
    value = QQ(value)
    return Fraction(int(value.numerator()), int(value.denominator()))


def finite_reduction_certificate(curve, points):
    coefficients = [
        Fraction(0),
        Fraction(0),
        Fraction(0),
        python_fraction(curve.a4()),
        python_fraction(curve.a6()),
    ]
    affine = [
        (python_fraction(point[0]), python_fraction(point[1])) for point in points
    ]
    signatures = find_mod2_reduction_certificate(
        coefficients, affine, prime_bound=500
    )
    rank = combined_mod2_rank(signatures, len(points))
    if rank != len(points):
        raise ArithmeticError(
            f"finite reductions reached rank {rank}, not {len(points)}"
        )
    no_two_torsion_prime = find_two_torsion_certificate_prime(
        coefficients, prime_bound=500
    )
    return {
        "proof": (
            "Full F2 rank in a product of E(F_p)/2E(F_p), plus the separate "
            "no-rational-2-torsion prime, proves independence by infinite descent."
        ),
        "combined_exact_rank_over_F2": rank,
        "two_torsion_certificate_prime": no_two_torsion_prime,
        "certificate_primes": [signature.prime for signature in signatures],
        "signatures": [
            {
                "prime": signature.prime,
                "group_order": signature.group_order,
                "doubled_subgroup_order": signature.doubled_subgroup_order,
                "quotient_dimension": signature.quotient_dimension,
                "rows": [list(row) for row in signature.rows],
            }
            for signature in signatures
        ],
    }


def exact_linear_combination(curve, coefficients, points):
    answer = curve(0)
    for coefficient, point in zip(coefficients, points):
        answer += int(coefficient) * point
    return answer


def recover_coordinates(curve, public_points, target_points):
    """Propose integral coordinates numerically and prove them exactly."""

    all_points = list(public_points) + list(target_points)
    heights = matrix(
        curve.height_pairing_matrix(all_points, precision=PRECISION)
    )
    count = len(public_points)
    public_gram = heights[:count, :count]
    cross = heights[:count, count:]
    real_coordinates = public_gram.solve_right(cross)
    coordinates = matrix(
        ZZ,
        count,
        len(target_points),
        lambda row, column: ZZ(real_coordinates[row, column].round()),
    )
    maximum_error = max(
        abs(real_coordinates[row, column] - coordinates[row, column])
        for row in range(count)
        for column in range(len(target_points))
    )
    if maximum_error >= 2 ** (-100):
        raise ArithmeticError(
            "height-coordinate recovery did not clear the 2^-100 separation gate"
        )
    for column, target in enumerate(target_points):
        if exact_linear_combination(
            curve, coordinates.column(column), public_points
        ) != target:
            raise ArithmeticError(f"exact relation {column + 1} failed")
    return coordinates


def evaluate_cover_splits(covers, parameter, ring):
    split_records = []
    indicator = hashlib.sha256()
    for cover in covers["bisections"]:
        branch = cover["branch"]
        q_value = (
            polynomial(ring, branch["numerator_coefficients"])(parameter)
            / polynomial(ring, branch["denominator_coefficients"])(parameter)
        )
        split = bool(q_value.is_square())
        indicator.update(
            f"{cover['label']}\t{int(split)}\n".encode()
        )
        if not split:
            continue
        if q_value == 0:
            raise ArithmeticError(f"ramified specialization at {cover['label']}")
        square_root = q_value.sqrt()
        lift = cover["lifted_section"]
        values = {}
        for coordinate in ("x", "y"):
            values[f"{coordinate}_plus"] = (
                polynomial(ring, lift[f"{coordinate}0_coefficients"])(parameter)
                + square_root
                * polynomial(ring, lift[f"{coordinate}1_coefficients"])(parameter)
            )
            values[f"{coordinate}_minus"] = (
                polynomial(ring, lift[f"{coordinate}0_coefficients"])(parameter)
                - square_root
                * polynomial(ring, lift[f"{coordinate}1_coefficients"])(parameter)
            )
        split_records.append(
            {
                "cover": cover,
                "q_value": q_value,
                "square_root": square_root,
                **values,
            }
        )
    return split_records, indicator.hexdigest()


def special_fibre(config, hit, public_record, direct, covers, ring):
    curve_id = int(public_record["id"])
    source_chart = config["source_chart"]
    native = next(
        record
        for record in hit["native_chart_twists"]
        if record["chart"] == source_chart
    )
    if native["twist"]["status"] != "QQ_ISOMORPHIC_UNTWISTED":
        raise ArithmeticError(f"curve {curve_id} is twisted in {source_chart}")
    parameter = QQ(native["native_parameter"]["numerator"]) / QQ(
        native["native_parameter"]["denominator"]
    )
    projective_scale_q = QQ(native["twist"]["quadratic_twist_parameter_q"])
    projective_scale_s = QQ(native["twist"]["qq_isomorphism_scale_s_with_s_squared_q"])
    if projective_scale_s**2 != projective_scale_q:
        raise ArithmeticError(f"curve {curve_id} has inconsistent projective scales")

    direct_A = polynomial(
        ring, direct["weierstrass_model"]["A_coefficients_low_to_high"]
    )
    direct_B = polynomial(
        ring, direct["weierstrass_model"]["B_coefficients_low_to_high"]
    )
    a1, a3, b2, target_A, target_B = short_invariants(public_record["ainvs"])
    fibre_A = direct_A(parameter)
    fibre_B = direct_B(parameter)
    scale_q = target_B * fibre_A / (fibre_B * target_A)
    if not scale_q.is_square():
        raise ArithmeticError(f"curve {curve_id} has a nontrivial affine-chart twist")
    scale_s = scale_q.sqrt()
    parameter_denominator = QQ(native["native_parameter"]["denominator"])
    if projective_scale_q * parameter_denominator**4 != scale_q:
        raise ArithmeticError(f"curve {curve_id} projective/affine q scales disagree")
    if projective_scale_s * parameter_denominator**2 not in (scale_s, -scale_s):
        raise ArithmeticError(f"curve {curve_id} projective/affine s scales disagree")
    if target_A != scale_q**2 * fibre_A:
        raise ArithmeticError(f"curve {curve_id} failed the exact A-isomorphism")
    if target_B != scale_q**3 * fibre_B:
        raise ArithmeticError(f"curve {curve_id} failed the exact B-isomorphism")
    curve = EllipticCurve(QQ, [target_A, target_B])

    public_points = [
        curve(
            QQ(x) + b2 / 12,
            QQ(y) + (a1 * QQ(x) + a3) / 2,
        )
        for x, y in public_record["points"]
    ]
    if len(public_points) != int(hit["snapshot_rank_lower_bound"]):
        raise ArithmeticError(f"curve {curve_id} public point count changed")
    independence = finite_reduction_certificate(curve, public_points)

    generic_points = []
    for record in direct["sections"]["records"]:
        x_value = rational_function(ring, record["X"])(parameter)
        y_value = rational_function(ring, record["Y"])(parameter)
        generic_points.append(
            curve(scale_q * x_value, scale_s**3 * y_value)
        )
    if len(generic_points) != 17:
        raise ArithmeticError("the saturated generic basis no longer has rank 17")

    split_records, split_digest = evaluate_cover_splits(covers, parameter, ring)
    split_plus = []
    split_minus = []
    for split in split_records:
        plus = curve(
            scale_q * split["x_plus"], scale_s**3 * split["y_plus"]
        )
        minus = curve(
            scale_q * split["x_minus"], scale_s**3 * split["y_minus"]
        )
        trace_word = split["cover"]["section_basis_w"]
        trace_point = exact_linear_combination(curve, trace_word, generic_points)
        if plus + minus != trace_point:
            raise ArithmeticError(
                f"split trace relation failed for {split['cover']['label']}"
            )
        split_plus.append(plus)
        split_minus.append(minus)

    all_target_points = generic_points + split_plus
    coordinates = recover_coordinates(curve, public_points, all_target_points)
    generic_coordinates = coordinates[:, :17]
    split_coordinates = coordinates[:, 17:]
    if generic_coordinates.rank() != 17:
        raise ArithmeticError(f"curve {curve_id} generic specialization lost rank")
    smith_diagonal, smith_left, _ = generic_coordinates.smith_form()
    elementary = generic_coordinates.elementary_divisors()
    quotient_rank = len(public_points) - 17
    if elementary != [ZZ(1)] * 17 + [ZZ(0)] * quotient_rank:
        raise ArithmeticError(
            f"curve {curve_id} displayed quotient is not free rank {quotient_rank}"
        )
    if smith_diagonal[:17, :17] != Matrix.identity(ZZ, 17):
        raise ArithmeticError("unexpected Smith diagonal normalization")

    public_quotient_vectors = []
    for index in range(len(public_points)):
        unit = matrix(ZZ, len(public_points), 1)
        unit[index, 0] = 1
        transformed = smith_left * unit
        public_quotient_vectors.append(
            [int(value) for value in transformed.column(0)[17:]]
        )
    simple_tail_basis = (
        generic_coordinates[17:, :] == matrix(ZZ, quotient_rank, 17)
        and abs(generic_coordinates[:17, :].det()) == 1
    )
    preferred_public_basis = (
        [f"P{index}" for index in range(18, len(public_points) + 1)]
        if simple_tail_basis
        else (
            CURVE12_PREFERRED_BASIS
            if curve_id == 12
            else None
        )
    )

    split_output = []
    quotient_columns = []
    for index, split in enumerate(split_records):
        public_vector = split_coordinates.column(index)
        quotient_vector = smith_left * matrix(ZZ, len(public_points), 1, list(public_vector))
        quotient_vector = quotient_vector.column(0)[17:]
        quotient_columns.append(list(quotient_vector))
        split_output.append(
            {
                "label": split["cover"]["label"],
                "priority_rank": int(split["cover"]["priority_rank"]),
                "branch_square_value": rational_text(split["q_value"]),
                "branch_square_root": rational_text(split["square_root"]),
                "trace_word_in_saturated_generic_basis": [
                    int(value) for value in split["cover"]["section_basis_w"]
                ],
                "plus_point_coordinates_in_ordered_public_points": [
                    int(value) for value in public_vector
                ],
                "plus_quotient_vector_in_deterministic_smith_basis": [
                    int(value) for value in quotient_vector
                ],
                "minus_quotient_vector_in_deterministic_smith_basis": [
                    -int(value) for value in quotient_vector
                ],
                "plus_vector_in_preferred_public_quotient_basis": (
                    [int(value) for value in public_vector[17:]]
                    if simple_tail_basis
                    else None
                ),
                "both_branches_and_trace_verified_by_exact_group_law": True,
            }
        )

    if quotient_columns:
        split_span = matrix(ZZ, quotient_columns).transpose()
        split_rank = split_span.rank()
        split_elementary = split_span.elementary_divisors()
        nonzero_split_invariants = [
            abs(int(value)) for value in split_elementary if value != 0
        ]
        annihilator = split_span.transpose().right_kernel_matrix()
        annihilator_rows = [
            [int(value) for value in row] for row in annihilator.rows()
        ]
        if simple_tail_basis:
            preferred_span = matrix(
                ZZ,
                [
                    [int(value) for value in split_coordinates.column(index)[17:]]
                    for index in range(len(split_records))
                ],
            ).transpose()
            preferred_annihilator_rows = [
                [int(value) for value in row]
                for row in preferred_span.transpose().right_kernel_matrix().rows()
            ]
        else:
            preferred_annihilator_rows = None
    else:
        split_rank = 0
        nonzero_split_invariants = []
        annihilator_rows = [
            [int(row == column) for column in range(quotient_rank)]
            for row in range(quotient_rank)
        ]
        preferred_annihilator_rows = None

    return {
        "curve_id": curve_id,
        "native_chart": source_chart,
        "representative_class": config["representative"],
        "native_parameter": rational_text(parameter),
        "snapshot_rank_lower_bound": int(hit["snapshot_rank_lower_bound"]),
        "displayed_point_count": len(public_points),
        "public_point_independence": independence,
        "specialized_generic_subgroup": {
            "rank": 17,
            "coordinate_matrix_rows_in_ordered_public_points": [
                [int(value) for value in row] for row in generic_coordinates.rows()
            ],
            "all_relations_verified_by_exact_group_law": True,
            "height_recovery_separation_gate": "maximum rounding error < 2^-100",
            "numerical_heights_used_in_proof": False,
        },
        "displayed_exceptional_quotient": {
            "quotient": f"Z^{quotient_rank}",
            "free_rank": quotient_rank,
            "smith_nonzero_invariant_factors": [1] * 17,
            "preferred_public_quotient_basis": preferred_public_basis,
            "preferred_basis_reason": (
                "The first seventeen public points contain the generic subgroup "
                "unimodularly and every later public point has zero generic coordinate."
                if simple_tail_basis
                else (
                    "The previously certified curve-12 unimodular complement."
                    if curve_id == 12
                    else "No unit-vector public complement was selected; use the Smith basis."
                )
            ),
            "public_point_images_in_deterministic_smith_basis": {
                f"P{index + 1}": vector
                for index, vector in enumerate(public_quotient_vectors)
            },
            "scope": (
                "quotient of the subgroup generated by the displayed public points; "
                "not a proved quotient of the full Mordell--Weil group"
            ),
        },
        "alternate_q80_cover_audit": {
            "covers_evaluated": len(covers["bisections"]),
            "complete_inventory_status": covers["status"],
            "split_indicator_sha256_in_inventory_order": split_digest,
            "rational_split_count": len(split_records),
            "splits": split_output,
            "exact_split_span_rank_in_exceptional_quotient": int(split_rank),
            "split_span_nonzero_smith_invariant_factors": nonzero_split_invariants,
            "split_span_is_primitive": all(value == 1 for value in nonzero_split_invariants),
            "split_span_annihilator_rows_in_deterministic_smith_dual": annihilator_rows,
            "split_span_annihilator_rows_in_preferred_public_quotient_dual": preferred_annihilator_rows,
        },
    }


def build():
    sweep = json.loads(SWEEP.read_text())
    public = json.loads(PUBLIC.read_text())
    local = json.loads(LOCAL.read_text())
    curve12 = json.loads(CURVE12.read_text())
    if local["status"] != "PASS_EXACT_LOCAL_FINGERPRINTS_FOR_ALL_69_RECOGNIZED_FIBRES":
        raise ArithmeticError("the all-69 local fingerprint input is not exact")
    if curve12["status"] != "PROVED_CURVE12_NATIVE_ALTERNATE_Q80_AND_DISPLAYED_QUOTIENT":
        raise ArithmeticError("the original curve-12 quotient certificate changed")
    if curve12["displayed_exceptional_quotient"]["free_basis_modulo_specialized_generic"] != CURVE12_PREFERRED_BASIS:
        raise ArithmeticError("the original curve-12 preferred quotient basis changed")
    hit_records = {
        int(record["curve_id"]): record
        for record in sweep["rational_j_hits_and_twists"]
    }
    public_records = {int(record["id"]): record for record in public["records"]}
    ring = PolynomialRing(QQ, "u")
    fibres = []
    inputs = {
        relative(SWEEP): digest(SWEEP),
        relative(PUBLIC): digest(PUBLIC),
        relative(LOCAL): digest(LOCAL),
        relative(CURVE12): digest(CURVE12),
    }
    for config in CHARTS:
        direct_path = config["direct"]
        covers_path = config["covers"]
        direct = json.loads(direct_path.read_text())
        covers = json.loads(covers_path.read_text())
        if direct["sections"]["status"] != "PASS_EXACT_SATURATED_RANK17_BASIS":
            raise ArithmeticError(f"{config['source_chart']} generic basis is not saturated")
        if len(covers["bisections"]) != int(covers["construction"]["record_count"]):
            raise ArithmeticError(f"{config['source_chart']} cover inventory is incomplete")
        inputs[relative(direct_path)] = digest(direct_path)
        inputs[relative(covers_path)] = digest(covers_path)
        for curve_id in config["curve_ids"]:
            print(
                f"R17Q80AUDIT|curve={curve_id}|stage=exact_specialization",
                flush=True,
            )
            fibres.append(
                special_fibre(
                    config,
                    hit_records[curve_id],
                    public_records[curve_id],
                    direct,
                    covers,
                    ring,
                )
            )

    if [record["curve_id"] for record in fibres] != [12, 395, 363, 364, 378, 393, 404]:
        raise ArithmeticError("the priority fibre order changed")
    if fibres[0]["displayed_exceptional_quotient"]["free_rank"] != 12:
        raise ArithmeticError("curve 12 lost its exact displayed Z^12 quotient")
    return {
        "schema": "elkies-k3.r17-norm12-native-icarm-quotient-audit.v1",
        "status": "PASS_EXACT_NATIVE_ALTERNATE_Q80_QUOTIENT_AUDIT_FOR_SEVEN_FIBRES",
        "summary": {
            "priority_fibres_audited": 7,
            "curve_ids": [12, 395, 363, 364, 378, 393, 404],
            "all_public_subgroups_independent_by_exact_finite_reductions": True,
            "all_generic_specializations_primitive_in_displayed_subgroups": True,
            "all_constructed_covers_exhaustively_evaluated": True,
            "all_69_recognized_fibres_have_local_fingerprints": True,
        },
        "fibres": fibres,
        "claim_boundary": {
            "proved": [
                "the exact displayed-subgroup quotient by generic MW17 for each of the seven fibres",
                "the complete rational-split inventory and its exact quotient span for each available native cover table",
                "the stated public point independence and every recovered group-law relation",
            ],
            "not_proved": [
                "that any displayed subgroup is the full Mordell--Weil group",
                "an unconditional upper bound or exact rank for any fibre",
                "that absence of a split in a fixed cover inventory means geometric invisibility to every alternate-Q80 or norm-eight construction",
                "that local/Nagao features predict rank jumps",
            ],
        },
        "inputs": inputs,
        "software_assumptions": {
            "sage_version": SAGE_VERSION,
            "required_features": [
                "exact QQ elliptic-curve group law",
                "exact Smith normal form",
                "canonical height pairing used only for candidate recovery",
            ],
        },
        "reproducing_command": (
            "PYTHONPATH=elliptic-curves/cas sage -python elkies-k3/scripts/"
            "certify_r17_norm12_native_icarm_quotient_audit.sage"
        ),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    output = args.output.resolve()
    serialized = json.dumps(build(), indent=2, sort_keys=True) + "\n"
    if args.check:
        if not output.exists() or output.read_text() != serialized:
            raise ArithmeticError("stored quotient audit differs from exact replay")
    else:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(serialized)
    print(
        "R17Q80AUDIT|fibres=7|curve12_quotient=Z^12|status=PASS|"
        f"output={relative(output)}",
        flush=True,
    )


if __name__ == "__main__":
    main()

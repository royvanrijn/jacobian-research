#!/usr/bin/env sage-python
"""Certify specialized MW17 quotients on every refreshed atlas hit.

The 573-curve atlas refresh supplies seventeen new untwisted norm-twelve
fibres, including eleven rank-at-least-25 priority fibres.  For each one this
replay specializes the exact saturated generic MW17 basis from its native
chart, independently proves the displayed public points independent, recovers
the generic subgroup coordinates by a high-precision height solve, verifies
every relation by exact group law, and computes the Smith quotient.  Numerical
heights discover relations only.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import runpy
import sys

from sage.all import EllipticCurve, GF, PolynomialRing, QQ, ZZ, matrix
from sage.env import SAGE_VERSION


sys.set_int_max_str_digits(0)
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "elliptic-curves/cas"))
SWEEP = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-icarm-database-sweep-v2.json"
OVERVIEW = ROOT / "artifacts/generated-results/elliptic-curves/icarm_curve_refresh_475_573_overview_v1.json"
COMMON = ROOT / "elkies-k3/scripts/certify_r17_norm12_native_icarm_quotient_audit.sage"
MOD3_SOURCE = ROOT / "elliptic-curves/cas/search_mestre_root_tuple_scale_max200.py"
OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-refresh-priority-quotients-v1.json"

from search_mestre_root_tuple_scale_max200 import (  # noqa: E402
    mod3_independence_certificate,
)

PINNED_HASHES = {
    SWEEP: "77a3c051111e7ead5ee2a6f88df4a975c2f5bdb87be1bfe4d88b195f293da50c",
    OVERVIEW: "1db137c4c006f774ad653b41b8c04ecc7b332d1104905dcb3f5eb23732904e3c",
    COMMON: "878981317bcc71f72aabfd5e88ce3051a629a55bfd793156305d764adf44516c",
    MOD3_SOURCE: "405a2b9f7653c89af0e3e6caf2e77765cb4bfc88fccf88edffa67d3435aebf24",
}

CONFIGS = (
    {
        "source_chart": "norm12-orbit-07ca9",
        "representative": "norm12-orbit-07ca9",
        "curve_ids": (543, 544, 545, 499),
        "direct": ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-orbit07ca9-direct-fibration-v1.json",
        "sha256": "a9c21568aa5f909013a951924f78f6a222f59b38c12f3048b8a3fb5febc1871b",
    },
    {
        "source_chart": "norm12-orbit-08234",
        "representative": "norm12-orbit-08234",
        "curve_ids": (531, 534, 535, 536, 537, 546, 478),
        "direct": ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-orbit08234-direct-fibration-v1.json",
        "sha256": "650c200300f884b266c316eae5bd6c7567c1707dfc21196827b8c18d84a16ddd",
    },
    {
        "source_chart": "norm12-orbit-11952",
        "representative": "norm12-orbit-11952",
        "curve_ids": (540, 532),
        "direct": ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json",
        "sha256": "76c54483c93c7090def42a8dad256838eb9510cd8479d07c5e3123eefa5cfe66",
    },
    {
        "source_chart": "norm12-orbit-103b2",
        "representative": "norm12-orbit-0e80b",
        "curve_ids": (541, 539),
        "direct": ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-orbit103b2-direct-fibration-v1.json",
        "sha256": "3f676dd0ce76da7f3092b073519b11af964c6982bfbc7262057d0f5c66234b9f",
    },
    {
        "source_chart": "norm12-orbit-08f72",
        "representative": "norm12-orbit-08f72",
        "curve_ids": (498,),
        "direct": ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-orbit08f72-direct-fibration-v1.json",
        "sha256": "5091c7fd2199692b538ecfd64b2bb3a2a208ac01ae28db082ca0dff6e21d1b90",
    },
    {
        "source_chart": "norm12-orbit-074d9",
        "representative": "norm12-orbit-074d9",
        "curve_ids": (538,),
        "direct": ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-wgxli-lineage-fibres-v1.json",
        "sha256": "15aeb70029a078dd51024a2e4c8d75a336a7197abe9471d79151a91178a59ec5",
        "format": "wgxli-lineage",
    },
)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def direct_payload(config, raw):
    """Put the older 074d9 lineage certificate in the direct-chart schema."""

    if config.get("format") != "wgxli-lineage":
        return raw
    representative = raw["representative"]
    if representative["chart"] != config["source_chart"]:
        raise ArithmeticError("the 074d9 lineage representative changed")
    generic = raw["generic_basis"]
    if generic["rank"] != 17 or not generic["saturated"]:
        raise ArithmeticError("the 074d9 lineage basis is not saturated rank 17")
    records = []
    for record in representative["sections"]:
        records.append(
            {
                "X": {
                    "numerator_coefficients_low_to_high": record[
                        "representative_x_coefficients_low_to_high"
                    ],
                    "denominator_coefficients_low_to_high": ["1"],
                },
                "Y": {
                    "numerator_coefficients_low_to_high": record[
                        "representative_y_coefficients_low_to_high"
                    ],
                    "denominator_coefficients_low_to_high": ["1"],
                },
            }
        )
    return {
        "weierstrass_model": {
            "A_coefficients_low_to_high": representative[
                "A_coefficients_low_to_high"
            ],
            "B_coefficients_low_to_high": representative[
                "B_coefficients_low_to_high"
            ],
        },
        "sections": {
            "status": "PASS_EXACT_SATURATED_RANK17_BASIS",
            "records": records,
        },
    }


def finite_reduction_certificate_with_mod3(curve, points, common):
    coefficients = [
        Fraction(0),
        Fraction(0),
        Fraction(0),
        common["python_fraction"](curve.a4()),
        common["python_fraction"](curve.a6()),
    ]
    affine = [
        (common["python_fraction"](point[0]), common["python_fraction"](point[1]))
        for point in points
    ]
    signatures = common["find_mod2_reduction_certificate"](
        coefficients, affine, prime_bound=500
    )
    mod2_rank = common["combined_mod2_rank"](signatures, len(points))
    no_two_torsion_prime = common["find_two_torsion_certificate_prime"](
        coefficients, prime_bound=500
    )
    mod3 = None
    certified_rank = mod2_rank
    if mod2_rank != len(points):
        mod3 = mod3_independence_certificate(
            coefficients, tuple(affine), prime_bound=499
        )
        certified_rank = int(mod3["combined_exact_rank_over_F3"])
    if certified_rank != len(points):
        raise ArithmeticError(
            f"finite reductions reached rank {certified_rank}, not {len(points)}"
        )
    return {
        "proof": (
            "Full rank in an exact product of good-reduction mod-2 or mod-3 "
            "quotients proves Z-independence by infinite descent."
        ),
        "certified_rank": certified_rank,
        "combined_exact_rank_over_F2": mod2_rank,
        "mod3_fallback_certificate": mod3,
        "two_torsion_certificate_prime": no_two_torsion_prime,
        "certificate_primes": [signature.prime for signature in signatures],
        "mod2_signatures": [
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


def curve499_noninclusion(config, hit, public_record, direct, ring, common):
    """Prove why the displayed curve-499 quotient is not defined."""

    curve_id = 499
    source_chart = config["source_chart"]
    native = next(
        record
        for record in hit["native_chart_twists"]
        if record["chart"] == source_chart
    )
    if native["twist"]["status"] != "QQ_ISOMORPHIC_UNTWISTED":
        raise ArithmeticError("curve 499 lost its untwisted native chart")
    parameter = QQ(native["native_parameter"]["numerator"]) / QQ(
        native["native_parameter"]["denominator"]
    )
    projective_q = QQ(native["twist"]["quadratic_twist_parameter_q"])
    projective_s = QQ(
        native["twist"]["qq_isomorphism_scale_s_with_s_squared_q"]
    )
    denominator = QQ(native["native_parameter"]["denominator"])
    direct_a = common["polynomial"](
        ring, direct["weierstrass_model"]["A_coefficients_low_to_high"]
    )
    direct_b = common["polynomial"](
        ring, direct["weierstrass_model"]["B_coefficients_low_to_high"]
    )
    a1, a3, b2, target_a, target_b = common["short_invariants"](
        public_record["ainvs"]
    )
    fibre_a = direct_a(parameter)
    fibre_b = direct_b(parameter)
    scale_q = target_b * fibre_a / (fibre_b * target_a)
    if not scale_q.is_square():
        raise ArithmeticError("curve 499 acquired a nontrivial affine twist")
    scale_s = scale_q.sqrt()
    if projective_q * denominator**4 != scale_q:
        raise ArithmeticError("curve 499 projective and affine q scales disagree")
    if projective_s * denominator**2 not in (scale_s, -scale_s):
        raise ArithmeticError("curve 499 projective and affine s scales disagree")
    if target_a != scale_q**2 * fibre_a or target_b != scale_q**3 * fibre_b:
        raise ArithmeticError("curve 499 failed the exact short-model isomorphism")
    curve = EllipticCurve(QQ, [target_a, target_b])
    public_points = [
        curve(
            QQ(x_value) + b2 / 12,
            QQ(y_value) + (a1 * QQ(x_value) + a3) / 2,
        )
        for x_value, y_value in public_record["points"]
    ]
    independence = finite_reduction_certificate_with_mod3(
        curve, public_points, common
    )
    generic_points = []
    for record in direct["sections"]["records"]:
        x_value = common["rational_function"](ring, record["X"])(parameter)
        y_value = common["rational_function"](ring, record["Y"])(parameter)
        generic_points.append(
            curve(scale_q * x_value, scale_s**3 * y_value)
        )
    heights = matrix(
        curve.height_pairing_matrix(public_points + generic_points, precision=300)
    )
    public_count = len(public_points)
    real_coordinates = heights[:public_count, :public_count].solve_right(
        heights[:public_count, public_count:]
    )
    triple_coordinates = matrix(
        ZZ,
        public_count,
        17,
        lambda row, column: ZZ((3 * real_coordinates[row, column]).round()),
    )
    maximum_error = max(
        abs(3 * real_coordinates[row, column] - triple_coordinates[row, column])
        for row in range(public_count)
        for column in range(17)
    )
    if maximum_error >= 2 ** (-100):
        raise ArithmeticError("curve 499 third-coordinate recovery did not separate")
    for column, generic_point in enumerate(generic_points):
        if 3 * generic_point != common["exact_linear_combination"](
            curve, triple_coordinates.column(column), public_points
        ):
            raise ArithmeticError("a curve 499 triple relation failed exact group law")
    if triple_coordinates.rank() != 17:
        raise ArithmeticError("curve 499 generic specialization lost rank")
    mod3_rank = matrix(GF(3), triple_coordinates).rank()
    if mod3_rank != 1:
        raise ArithmeticError("the curve 499 commensurability index changed")
    return {
        "curve_id": curve_id,
        "status": "PROVED_GENERIC_SUBGROUP_NOT_CONTAINED_IN_DISPLAYED_SUBGROUP",
        "native_chart": source_chart,
        "representative_class": config["representative"],
        "native_parameter": common["rational_text"](parameter),
        "displayed_subgroup_rank": public_count,
        "specialized_generic_subgroup_rank": 17,
        "public_point_independence": independence,
        "three_times_generic_coordinates_in_ordered_public_points_rows": [
            [int(value) for value in row] for row in triple_coordinates.rows()
        ],
        "all_triple_relations_verified_by_exact_group_law": True,
        "height_recovery_separation_gate": "maximum error in 3*C is below 2^-100",
        "numerical_heights_used_in_proof": False,
        "rank_of_triple_coordinate_matrix_modulo_3": int(mod3_rank),
        "overgroup_generated_by_displayed_and_generic_modulo_displayed": "Z/3Z",
        "displayed_subgroup_contains_specialized_generic_subgroup": False,
        "displayed_quotient_by_generic_subgroup": "NOT_DEFINED_NONINCLUSION",
        "proof": (
            "Exact relations 3*G_j=sum_i N_ij*P_i put the generic subgroup in "
            "the rational span of the independent displayed points. Since N mod 3 "
            "has rank one, adjoining the generic subgroup enlarges the displayed "
            "subgroup by Z/3Z, so the generic subgroup is not contained in it."
        ),
    }


def build():
    for path, expected in PINNED_HASHES.items():
        observed = digest(path)
        if observed != expected:
            raise ArithmeticError(f"pinned input changed: {relative(path)}")
    common = runpy.run_path(str(COMMON))
    special_fibre = common["special_fibre"]
    special_fibre.__globals__["finite_reduction_certificate"] = (
        lambda curve, points: finite_reduction_certificate_with_mod3(
            curve, points, common
        )
    )

    sweep = json.loads(SWEEP.read_text())
    overview = json.loads(OVERVIEW.read_text())
    if sweep["status"] != "PASS_EXACT_COMPLETE_PINNED_ICARM_J_PREIMAGE_AND_TWIST_SWEEP":
        raise ArithmeticError("the refresh sweep is not exact")
    if overview["status"] != "PASS_EXACT_OVERVIEW_OF_ICARM_CURVES_475_THROUGH_573":
        raise ArithmeticError("the refresh overview is not exact")
    priority_ids = [543, 531, 534, 535, 536, 544, 545, 537, 540, 541, 546]
    if overview["summary"]["priority_atlas_fibre_ids"] != priority_ids:
        raise ArithmeticError("the rank-at-least-24 atlas-priority list changed")
    expected_ids = priority_ids + [498, 539, 538, 478, 499, 532]
    overview_hits = sorted(
        record["curve_id"]
        for record in overview["curve_overview"]
        if record["atlas"]["status"] == "EXACT_UNTWISTED_NORM12_ATLAS_FIBRE"
    )
    if overview_hits != sorted(expected_ids):
        raise ArithmeticError("the refreshed atlas-hit inventory changed")

    hits = {
        int(record["curve_id"]): record
        for record in sweep["rational_j_hits_and_twists"]
    }
    public = {
        int(record["id"]): record for record in overview["snapshot"]["records"]
    }
    ring = PolynomialRing(QQ, "u")
    fibres = []
    noninclusive = []
    inputs = {relative(path): expected for path, expected in PINNED_HASHES.items()}
    for config in CONFIGS:
        direct_path = config["direct"]
        if digest(direct_path) != config["sha256"]:
            raise ArithmeticError(f"native direct fibration changed: {relative(direct_path)}")
        direct = direct_payload(config, json.loads(direct_path.read_text()))
        if direct["sections"]["status"] != "PASS_EXACT_SATURATED_RANK17_BASIS":
            raise ArithmeticError(f"{config['source_chart']} generic basis is not saturated")
        inputs[relative(direct_path)] = config["sha256"]
        for curve_id in config["curve_ids"]:
            hit = hits.get(curve_id)
            if hit is None or hit["representative"] != config["representative"]:
                raise ArithmeticError(f"curve {curve_id} lost its exact atlas class")
            print(f"R17REFRESHQ|curve={curve_id}|stage=specialize-and-recover", flush=True)
            if curve_id == 499:
                noninclusive.append(
                    curve499_noninclusion(
                        config,
                        hit,
                        public[curve_id],
                        direct,
                        ring,
                        common,
                    )
                )
                continue
            fibres.append(
                special_fibre(
                    config,
                    hit,
                    public[curve_id],
                    direct,
                    None,
                    ring,
                )
            )

    fibres.sort(key=lambda record: (-record["snapshot_rank_lower_bound"], record["curve_id"]))
    quotient_ids = [curve_id for curve_id in expected_ids if curve_id != 499]
    if [record["curve_id"] for record in fibres] != quotient_ids:
        raise ArithmeticError("the refreshed quotient output order changed")
    if [record["curve_id"] for record in noninclusive] != [499]:
        raise ArithmeticError("the refreshed noninclusion inventory changed")
    quotient_by_curve = {
        str(record["curve_id"]): record["displayed_exceptional_quotient"]["quotient"]
        for record in fibres
    }
    expected_quotients = {
        "543": "Z^12",
        "531": "Z^11",
        "534": "Z^11",
        "535": "Z^11",
        "536": "Z^11",
        "544": "Z^11",
        "545": "Z^11",
        "537": "Z^10",
        "540": "Z^8",
        "541": "Z^8",
        "546": "Z^8",
        "498": "Z^6",
        "539": "Z^6",
        "538": "Z^5",
        "478": "Z^4",
        "532": "Z^3",
    }
    if quotient_by_curve != expected_quotients:
        raise ArithmeticError("a refreshed displayed quotient changed")
    return {
        "schema": "elkies-k3.r17-norm12-refresh-priority-quotients.v1",
        "status": "PASS_EXACT_REFRESH_ATLAS_HIT_SPECIALIZATION_AUDIT",
        "summary": {
            "new_atlas_hit_count": len(expected_ids),
            "new_atlas_hit_ids": expected_ids,
            "quotient_curve_count": len(fibres),
            "quotient_curve_ids": quotient_ids,
            "priority_curve_count": len(priority_ids),
            "priority_curve_ids": priority_ids,
            "nonpriority_curve_ids": [498, 539, 538, 478, 499, 532],
            "all_public_subgroups_independent_by_exact_finite_reductions": True,
            "all_generic_specializations_have_rank_17": True,
            "all_quotient_generic_subgroups_primitive_in_displayed_subgroups": True,
            "noninclusive_displayed_subgroup_curve_ids": [499],
            "displayed_exceptional_quotients": quotient_by_curve,
            "curve543_rank_lower_bound": 29,
            "curve543_displayed_exceptional_quotient": "Z^12",
        },
        "fibres": fibres,
        "noninclusive_fibres": noninclusive,
        "claim_boundary": {
            "proved": [
                "all seventeen refreshed atlas hits are untwisted fibres of the stated native norm-twelve charts",
                "their displayed public points are independent by exact finite reductions",
                "on sixteen hits the saturated generic MW17 bases specialize primitively into the displayed subgroups and the displayed-subgroup quotients are the stated free abelian groups",
                "on curve 499 the specialized generic subgroup is not contained in the displayed subgroup and adjoining it enlarges that subgroup by Z/3Z",
            ],
            "not_proved": [
                "that any displayed subgroup is the full Mordell-Weil group",
                "an unconditional rank upper bound or exact rank for any fibre",
                "the geometric origin of every exceptional quotient direction",
                "a cover-visibility result, because no cover inventory is evaluated here",
            ],
        },
        "inputs": inputs,
        "software_assumptions": {
            "sage_version": SAGE_VERSION,
            "required_features": [
                "exact QQ elliptic-curve group law",
                "exact Smith normal form",
                "canonical heights used only for candidate relation recovery",
            ],
        },
        "reproducing_command": (
            "PYTHONPATH=elliptic-curves/cas sage -python elkies-k3/scripts/"
            "certify_r17_norm12_refresh_priority_quotients.sage"
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
        if not output.is_file() or output.read_text() != serialized:
            raise ArithmeticError("the stored refresh-priority quotient certificate changed")
    else:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(serialized)
    print(
        "R17REFRESHQ|curves=17|priority=11|curve543=Z^12|rank28_fibres=Z^11|"
        f"status=PROVED|output={relative(output)}",
        flush=True,
    )


if __name__ == "__main__":
    main()

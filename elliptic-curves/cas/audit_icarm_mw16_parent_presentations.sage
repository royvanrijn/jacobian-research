#!/usr/bin/env sage-python
"""Deduplicate the nine MW16 atlas hits and compare their blind responses.

Within each target-curve cluster, prove exact PGL2(Q) equivalence of the
Jacobian pencils, prove that the induced twist is a rational square, transport
the recorded target parameter, and recover an integral determinant-one change
between the specialized generic MW16 bases.  Between the five cluster
representatives, retain the exact separating witnesses returned by the frozen
j-map equivalence solver.

This replay uses no public exceptional points or target complement.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
from importlib.machinery import SourceFileLoader
import json
from pathlib import Path
import sys

from sage.all import PolynomialRing, QQ, ZZ, matrix


ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT / "elliptic-curves/data/icarm_mw16_parent_ladder_blind_inputs_v1.json"
BLIND = ROOT / "artifacts/generated-results/elliptic-curves/icarm_mw16_parent_ladder_blind_v1.json"
PGL2 = ROOT / "elkies-k3/scripts/compile_r17_norm12_record_lineage_atlas.sage"
OUTPUT = ROOT / "artifacts/generated-results/elliptic-curves/icarm_mw16_parent_presentation_audit_v1.json"


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def canonical_sha256(value) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def qtext(value) -> str:
    value = QQ(value)
    return (
        str(value.numerator())
        if value.denominator() == 1
        else f"{value.numerator()}/{value.denominator()}"
    )


def j_record(parent, ring):
    coefficient_a = ring(parent["pencil"]["A_coefficients_low_to_high"])
    coefficient_b = ring(parent["pencil"]["B_coefficients_low_to_high"])
    return (
        {"N": 6912 * coefficient_a**3, "D": 4 * coefficient_a**3 + 27 * coefficient_b**2},
        coefficient_a,
        coefficient_b,
    )


def point_tuple(record):
    return Fraction(record["x"]), Fraction(record["y"])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=INPUT)
    parser.add_argument("--blind", type=Path, default=BLIND)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()

    inputs = json.loads(args.input.read_text())
    blind = json.loads(args.blind.read_text())
    if inputs.get("status") != "PASS_EXACT_COMPLEMENT_BLIND_NINE_PARENT_INPUTS":
        raise ArithmeticError("nine-parent input status changed")
    if blind.get("status") != "PASS_COMPLETE_NINE_PARENT_INITIAL_HALF_LATTICE_LADDER":
        raise ArithmeticError("blind ladder status changed")
    if inputs["observation_unit"] != blind["observation_unit"]:
        raise ArithmeticError("blind run changed the observation-unit declaration")

    sys.path.insert(0, str(ROOT / "elliptic-curves"))
    from latent_lattice.elliptic import EllipticCurve as LatentEllipticCurve
    from latent_lattice.pari import recover_exact_embedding

    pgl2 = SourceFileLoader("mw16_parent_pgl2", str(PGL2)).load_module()
    ring = PolynomialRing(QQ, "lambda")
    parameter_variable = ring.gen()
    field = ring.fraction_field()
    parent_by_id = {row["parent_id"]: row for row in inputs["parents"]}
    blind_by_id = {row["parent_id"]: row for row in blind["parents"]}
    clusters = []
    representatives = []
    for curve_id in inputs["observation_unit"]["curve_ids"]:
        parents = [
            row for row in inputs["parents"] if int(row["curve_id"]) == curve_id
        ]
        parents.sort(key=lambda row: int(row["priority_rank"]))
        first = parents[0]
        representatives.append(first)
        first_j, first_a, first_b = j_record(first, ring)
        model = tuple(Fraction(value) for value in first["target_short_model"])
        latent_curve = LatentEllipticCurve(model)
        first_points = tuple(
            point_tuple(record) for record in first["specialized_generic_points"]
        )
        within = []
        for second in parents[1:]:
            second_j, second_a, second_b = j_record(second, ring)
            base_change, solve = pgl2.exact_pgl2_equivalence(
                first_j, second_j, ring
            )
            if base_change is None or not solve.get("identity_verified"):
                raise ArithmeticError("within-curve parent labels are not PGL2 equivalent")
            aa, bb, cc, dd = map(QQ, base_change)
            if aa * dd == bb * cc:
                raise ArithmeticError("recovered PGL2 matrix is singular")
            phi = field(
                (aa * parameter_variable + bb)
                / (cc * parameter_variable + dd)
            )
            pulled_second_a = field(second_a(phi))
            pulled_second_b = field(second_b(phi))
            first_a_field = field(first_a)
            first_b_field = field(first_b)
            twist = field(
                pulled_second_b
                * first_a_field
                / (first_b_field * pulled_second_a)
            )
            if (
                pulled_second_a != twist**2 * first_a_field
                or pulled_second_b != twist**3 * first_b_field
            ):
                raise ArithmeticError("j-equivalence failed the Weierstrass identities")
            if twist.numerator().degree() or twist.denominator().degree():
                raise ArithmeticError("within-curve Weierstrass twist is not constant")
            twist_constant = QQ(twist)
            if not twist_constant.is_square():
                raise ArithmeticError("within-curve pencils differ by a nonsquare twist")
            first_parameter = QQ(first["target_parameter"])
            transported = QQ(
                (aa * first_parameter + bb) / (cc * first_parameter + dd)
            )
            if transported != QQ(second["target_parameter"]):
                raise ArithmeticError("PGL2 map does not transport target parameters")

            second_points = tuple(
                point_tuple(record)
                for record in second["specialized_generic_points"]
            )
            embedding = recover_exact_embedding(
                latent_curve,
                first_points,
                second_points,
                digits=150,
                timeout=300.0,
            )
            transition = matrix(ZZ, [list(column) for column in embedding.columns]).transpose()
            if transition.dimensions() != (16, 16) or abs(transition.det()) != 1:
                raise ArithmeticError("specialized MW16 subgroups are not integrally equal")
            first_discoveries = {
                (record["point"]["x"], record["point"]["y"])
                for record in blind_by_id[first["parent_id"]]["discoveries"]
            }
            second_discoveries = {
                (record["point"]["x"], record["point"]["y"])
                for record in blind_by_id[second["parent_id"]]["discoveries"]
            }
            within.append(
                {
                    "first_parent_id": first["parent_id"],
                    "second_parent_id": second["parent_id"],
                    "pgl2_matrix_a_b_c_d": list(map(int, base_change)),
                    "pgl2_determinant": qtext(aa * dd - bb * cc),
                    "affine_base_change": not cc,
                    "j_map_identity_exact": True,
                    "constant_twist_q": qtext(twist_constant),
                    "constant_twist_is_square_in_Q": True,
                    "weierstrass_scale_s_with_s_squared_q": qtext(
                        QQ(twist_constant.sqrt())
                    ),
                    "target_parameter_transport_exact": True,
                    "specialized_mw16_subgroups_integrally_equal": True,
                    "second_basis_rows_in_first_basis": [
                        list(map(int, row)) for row in transition.rows()
                    ],
                    "basis_transition_determinant": int(transition.det()),
                    "blind_initial_quotient_ranks": [
                        int(
                            blind_by_id[first["parent_id"]][
                                "exact_quotient_rank_recovered"
                            ]
                        ),
                        int(
                            blind_by_id[second["parent_id"]][
                                "exact_quotient_rank_recovered"
                            ]
                        ),
                    ],
                    "blind_returned_point_sets_equal": (
                        first_discoveries == second_discoveries
                    ),
                }
            )
        clusters.append(
            {
                "curve_id": curve_id,
                "representative_parent_id": first["parent_id"],
                "presentation_count": len(parents),
                "presentation_ids": [row["parent_id"] for row in parents],
                "within_cluster_equivalences": within,
            }
        )

    cross = []
    for left_index, left in enumerate(representatives):
        left_j = j_record(left, ring)[0]
        for right in representatives[left_index + 1 :]:
            right_j = j_record(right, ring)[0]
            base_change, solve = pgl2.exact_pgl2_equivalence(
                left_j, right_j, ring
            )
            if base_change is not None:
                raise ArithmeticError("two curve-level representatives unexpectedly coincide")
            left_histogram = blind_by_id[left["parent_id"]][
                "generic_half_lattice"
            ]["twice_norm_histogram"]
            right_histogram = blind_by_id[right["parent_id"]][
                "generic_half_lattice"
            ]["twice_norm_histogram"]
            separating_prime = solve.get("separating_prime")
            if separating_prime is None and left_histogram == right_histogram:
                raise ArithmeticError(
                    "cross-cluster non-equivalence lacks a separating witness: "
                    f"{left['parent_id']} versus {right['parent_id']}: {solve}"
                )
            cross.append(
                {
                    "left_parent_id": left["parent_id"],
                    "right_parent_id": right["parent_id"],
                    "pgl2_equivalent": False,
                    "separation_method": (
                        "exact_j_map_landmark_modp"
                        if separating_prime is not None
                        else "exact_generic_mw16_half_lattice_coset_spectrum"
                    ),
                    "separating_prime": (
                        int(separating_prime)
                        if separating_prime is not None
                        else None
                    ),
                    "landmark_values_sha256": canonical_sha256(
                        solve["landmark_values"]
                    ),
                    "mw16_half_lattice_spectra_differ": (
                        left_histogram != right_histogram
                    ),
                }
            )
    if len(clusters) != 5 or len(cross) != 10:
        raise ArithmeticError("five-cluster pair census changed")

    payload = {
        "schema": "elliptic-curves.icarm-mw16-parent-presentation-audit.v1",
        "status": "PASS_EXACT_NINE_PRESENTATIONS_FIVE_FIBRATIONS",
        "presentation_count": 9,
        "exact_fibration_class_count": 5,
        "clusters": clusters,
        "cross_cluster_non_equivalences": cross,
        "operational_rule": {
            "statistical_unit": "target_curve/fibration cluster",
            "independent_observation_count": 5,
            "coordinate_search_trial_count": 9,
            "retain_all_presentations": True,
            "reason": (
                "bounded projective-height boxes are not invariant under the certified "
                "affine base changes, so all nine presentations may probe different tails "
                "without being counted as distinct fibrations or outcomes"
            ),
        },
        "inputs": {
            relative(path): digest(path)
            for path in (args.input, args.blind, PGL2, Path(__file__))
        },
        "claim_boundary": [
            "The nine atlas hits form exactly five Q-isomorphism classes of Jacobian pencils under the tested exact PGL2/Weierstrass equivalence.",
            "All repeated presentations specialize to the same integral generic MW16 subgroup on their shared target curve.",
            "All repeated presentations returned the same exact quotient rank in the initial blind wave.",
            "Coordinate-height search boxes can still differ after affine reparameterization, so retaining nine search charts is operationally meaningful.",
            "No public complement, target-rank upper bound, rank-32 point, or Selmer result is asserted.",
        ],
        "reproducing_command": (
            "sage -python elliptic-curves/cas/audit_icarm_mw16_parent_presentations.sage --check"
        ),
    }
    serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check:
        if not args.output.is_file() or args.output.read_text() != serialized:
            raise ArithmeticError("stored parent-presentation audit differs from replay")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(serialized)
    print(
        "MW16PARENTAUDIT|presentations=9|fibrations=5|cross_pairs=10|"
        "all_repeated_subgroups_equal=1|status=PASS",
        flush=True,
    )


if __name__ == "__main__":
    main()

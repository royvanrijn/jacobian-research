#!/usr/bin/env sage-python
"""Retrospective 512-parity M26 fibre scan for the five remaining strict 302 labels.

This repeats the declared M25 fibre diagnostic after the independently
certified rank-26 recovery.  It holds fixed each direction's original M17
vetted shell parity and enumerates all nine parity bits newly available in
M26.  It is diagnostic-only: targets are known here and no point search or
prospective selection follows from the result automatically.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.machinery
import json
from fractions import Fraction
from pathlib import Path

from sage.all import EllipticCurve, QQ, ZZ, matrix, pari, vector


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
ART = ROOT / "artifacts/generated-results/elliptic-curves"
VISIBILITY = ART / "curve302_residual_visibility_geometry_v1.json"
M26 = ART / "curve302_m25_fibre_orbit64677_extension103_mod2_v1.json"
HELPER_SOURCE = CAS / "audit_curve302_m25_remaining_strict_cvp.sage"
MAPPER = CAS / "factor_free_pari_mapping.sage"
CVP_SOURCE = CAS / "audit_curve302_residual_visibility_geometry.sage"
OUTPUT = ART / "curve302_m26_remaining_strict_fibre_scan_v1.json"
M17 = 17
RANK = 26
EXTENSION_DIMENSION = 9
PRECISION = 192
EXCLUDED = {"residual-strict-04", "residual-strict-06"}


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rel(path: Path) -> str:
    return str(path.relative_to(ROOT))


def load(path: Path):
    return json.loads(path.read_text())


def ints(row):
    return [int(value) for value in row]


def parity_mask(row) -> int:
    return sum((int(value) & 1) << index for index, value in enumerate(row))


def build() -> dict:
    visibility, m26 = load(VISIBILITY), load(M26)
    if visibility["status"] != "PASS_RETROSPECTIVE_VETTED_VISIBILITY_DIAGNOSTIC":
        raise ArithmeticError("M17 visibility input changed")
    if m26["status"] != "COMPLETE_DECLARED_FINITE_AUDIT" or m26["rank_lower_bound"] != RANK:
        raise ArithmeticError("M26 certified input changed")
    helpers = importlib.machinery.SourceFileLoader("curve302_m26_helpers", str(HELPER_SOURCE)).load_module()
    cvp = importlib.machinery.SourceFileLoader("curve302_m26_cvp", str(CVP_SOURCE)).load_module()
    mapper = importlib.machinery.SourceFileLoader("curve302_m26_mapper", str(MAPPER)).load_module()
    from alternate_quartic_covers import short_add
    from half_lattice_pointed_sieve import linear_combination

    model = tuple(Fraction(value) for value in m26["curve"])
    curve = EllipticCurve(QQ, [QQ(value) for value in m26["curve"]])
    basis = tuple(tuple(QQ(value) for value in point) for point in m26["independent_points"])
    map_basis = tuple(tuple(Fraction(value) for value in point) for point in m26["independent_points"])
    if len(basis) != RANK or any(point not in curve for point in basis):
        raise ArithmeticError("M26 basis is malformed")
    fixed = {
        row["id"]: row for row in visibility["directions"]
        if row["cohort"] == "residual_strict" and row["id"] not in EXCLUDED
    }
    if len(fixed) != 5:
        raise ArithmeticError("remaining residual strict panel differs")
    prepared = []
    for direction_id, entry in sorted(fixed.items()):
        old_parity = entry["minimum_finite_coordinate_parity_in_vetted_set"]
        chart = next(row for row in entry["exact_pointed_quartic_rows"] if row["parity"] == old_parity)
        target = tuple(QQ(value) for value in chart["exact_target_point_short_model"])
        if target not in curve:
            raise ArithmeticError("retrospective target is off the M26 curve")
        prepared.append((direction_id, entry, chart, target, tuple(Fraction(value) for value in chart["exact_target_point_short_model"])))
    heights, asymmetry = helpers.rounded_gram(curve, (*basis, *(row[3] for row in prepared)))
    gram = heights[:RANK, :RANK]
    change = matrix(ZZ, pari(gram).qflllgram()).transpose()
    if abs(change.det()) != 1:
        raise ArithmeticError("M26 LLL change is not unimodular")
    inverse = change.inverse()
    reduced = change * gram * change.transpose()
    if not reduced.is_positive_definite():
        raise ArithmeticError("M26 reduced metric is not positive definite")
    dd = cvp.RoundedMetricCVP(reduced, "dd", PRECISION)
    mpfr = cvp.RoundedMetricCVP(reduced, "mpfr", PRECISION)
    directions = []
    for target_index, (direction_id, entry, chart, target, target_for_map) in enumerate(prepared):
        target_cross = vector(ZZ, [2 * heights[row, RANK + target_index] for row in range(RANK)])
        projection_original = gram.solve_right(target_cross)
        projection = vector(QQ, (matrix(QQ, 1, RANK, projection_original).row(0) * inverse))
        base_parity = [int(value) & 1 for value in chart["centre_m17_word"]]
        trials, exact = [], {}
        for extension in range(1 << EXTENSION_DIMENSION):
            original_parity = vector(ZZ, base_parity + [(extension >> bit) & 1 for bit in range(EXTENSION_DIMENSION)])
            reduced_parity = vector(ZZ, (matrix(ZZ, 1, RANK, original_parity) * inverse).row(0))
            reduced_parity = vector(ZZ, [int(value) % 2 for value in reduced_parity])
            u, centre_distance = cvp.cross_precision_closest(dd, mpfr, -vector(QQ, reduced_parity) / 2)
            centre = vector(ZZ, (matrix(ZZ, 1, RANK, reduced_parity + 2 * u) * change).row(0))
            z, target_distance = cvp.cross_precision_closest(dd, mpfr, vector(QQ, [(projection[index] - reduced_parity[index]) / 2 for index in range(RANK)]))
            translation = vector(ZZ, (matrix(ZZ, 1, RANK, u - z) * change).row(0))
            added = linear_combination(model, map_basis, ints(translation))
            moved = target_for_map if added is None else short_add(model, target_for_map, added)
            if moved is None:
                raise ArithmeticError("target translate became infinity")
            mapped = helpers.coordinate(mapper, model, map_basis, centre, moved)
            row = {
                "extension_mask_in_new_M26_coordinates": extension,
                "M26_parity_mask": parity_mask(centre),
                "centre_M26_word": ints(centre),
                "centre_metric_norm": int(centre * gram * centre),
                "centre_CVP_distance_in_reduced_metric": str(centre_distance),
                "target_translation_M26_word": ints(translation),
                "target_CVP_distance_in_reduced_metric": str(target_distance),
                "reduced_coordinate": mapped["reduced_coordinate"],
                "reduced_coordinate_height": mapped["reduced_coordinate_height"],
                "reduced_coordinate_decimal_digits": mapped["reduced_coordinate_decimal_digits"],
            }
            trials.append(row)
            exact[extension] = mapped
        trials.sort(key=lambda row: (int(row["reduced_coordinate_height"]), row["extension_mask_in_new_M26_coordinates"]))
        winner = trials[0]
        directions.append({
            "id": direction_id,
            "fixed_M17_chart": {
                "parity": chart["parity"],
                "centre_M17_word": chart["centre_m17_word"],
                "degree_two_translation_orbit": chart["degree_two_translation_orbit"],
                "old_vetted_coordinate_height": entry["minimum_finite_reduced_coordinate_height_in_vetted_set"],
                "old_vetted_coordinate_decimal_digits": entry["minimum_finite_reduced_coordinate_decimal_digits_in_vetted_set"],
            },
            "M26_extension_count": len(trials),
            "M26_extension_winner": {**winner, **exact[winner["extension_mask_in_new_M26_coordinates"]]},
            "M26_extension_coordinate_trials": trials,
            "inside_historical_height_125000": int(winner["reduced_coordinate_height"]) <= 125000,
        })
        print("CURVE302M26FIBRE|{}|winner_digits={}|extension={}".format(direction_id, winner["reduced_coordinate_decimal_digits"], winner["extension_mask_in_new_M26_coordinates"]), flush=True)
    return {
        "schema": "elliptic-curves.curve302-m26-remaining-strict-fibre-scan.v1",
        "status": "PASS_RETROSPECTIVE_M26_PARITY_FIBRE_SCAN",
        "inputs": {rel(path): sha(path) for path in (VISIBILITY, M26, HELPER_SOURCE, MAPPER, CVP_SOURCE, Path(__file__))},
        "protocol": {
            "directions": 5,
            "parities_per_direction": 1 << EXTENSION_DIMENSION,
            "fixed_M17_shell_projection": "hold each original M17 chart parity fixed and enumerate every newly available M26 parity bit",
            "centre_and_target_reduction": "DD/MPFR-agreeing CVP in a rank-26 rounded 384-bit canonical-height metric",
            "boundary": "Known targets are diagnostic inputs. The result is not a full M26 parity scan, a global coordinate minimum, a point search, or a prospective selection claim.",
        },
        "maximum_height_matrix_asymmetry": asymmetry,
        "directions": directions,
        "reproducing_command": "sage -python elliptic-curves/cas/audit_curve302_m26_remaining_strict_fibre_scan.sage --check",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build", action="store_true")
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()
    if arguments.build == arguments.check:
        parser.error("choose exactly one of --build or --check")
    payload = build()
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if arguments.build:
        if OUTPUT.exists():
            raise FileExistsError("preserve immutable M26 parity-fibre diagnostic")
        OUTPUT.write_text(rendered)
    elif OUTPUT.read_text() != rendered:
        raise ArithmeticError("stored M26 parity-fibre diagnostic did not replay")
    print("CURVE302M26FIBRE|directions=5|parities_per_direction=512|status=PASS", flush=True)


if __name__ == "__main__":
    main()

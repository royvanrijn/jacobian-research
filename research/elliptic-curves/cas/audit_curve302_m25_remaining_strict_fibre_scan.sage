#!/usr/bin/env sage-python
"""Retrospectively scan M25 parity extensions of six residual strict charts.

For each remaining residual strict target, retain its fixed M17 chart parity
from the previously declared 32-chart diagnostic.  Enumerate exactly the 256
extensions of that parity in the eight newly available M25 coordinates.  In
each extension choose the metric-shortest centre and the target-relative M25
translate by cross-precision CVP, then verify the pointed-quartic coordinate
exactly.  This is a finite retrospective fibre scan, not a point search nor a
global M25/2M25 visibility theorem.
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
M25 = ART / "curve302_m24_bisection_orbit117420_mod2_v1.json"
M25_CVP = ART / "curve302_m25_remaining_strict_cvp_v1.json"
M25_CVP_SOURCE = CAS / "audit_curve302_m25_remaining_strict_cvp.sage"
MAPPER = CAS / "factor_free_pari_mapping.sage"
CVP_SOURCE = CAS / "audit_curve302_residual_visibility_geometry.sage"
OUTPUT = ART / "curve302_m25_remaining_strict_fibre_scan_v1.json"
M17 = 17
RANK = 25
EXTENSION_DIMENSION = 8
PRECISION = 192


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
    visibility, m25, prior = load(VISIBILITY), load(M25), load(M25_CVP)
    if visibility["status"] != "PASS_RETROSPECTIVE_VETTED_VISIBILITY_DIAGNOSTIC":
        raise ArithmeticError("residual M17 visibility input is unavailable")
    if m25["rank_lower_bound"] != RANK or m25["status"] != "COMPLETE_DECLARED_FINITE_AUDIT":
        raise ArithmeticError("certified M25 input changed")
    if prior["status"] != "PASS_RETROSPECTIVE_FIXED_CHART_M25_CVP":
        raise ArithmeticError("fixed-chart M25 CVP diagnostic is unavailable")
    cvp = importlib.machinery.SourceFileLoader("curve302_m25_fibre_cvp", str(CVP_SOURCE)).load_module()
    helpers = importlib.machinery.SourceFileLoader("curve302_m25_fibre_helpers", str(M25_CVP_SOURCE)).load_module()
    mapper = importlib.machinery.SourceFileLoader("curve302_m25_fibre_mapper", str(MAPPER)).load_module()
    from alternate_quartic_covers import short_add
    from half_lattice_pointed_sieve import linear_combination

    model = tuple(Fraction(value) for value in m25["curve"])
    curve = EllipticCurve(QQ, [QQ(value) for value in m25["curve"]])
    basis = tuple(tuple(QQ(value) for value in point) for point in m25["independent_points"])
    map_basis = tuple(tuple(Fraction(value) for value in point) for point in m25["independent_points"])
    if len(basis) != RANK or any(point not in curve for point in basis):
        raise ArithmeticError("M25 basis is malformed")
    gram = matrix(ZZ, prior["M25_rounded_height_gram"])
    change = matrix(ZZ, prior["M25_LLL_change"])
    if gram.nrows() != RANK or change.nrows() != RANK or abs(change.det()) != 1:
        raise ArithmeticError("retained M25 metric transport changed")
    reduced = change * gram * change.transpose()
    inverse = change.inverse()
    if not reduced.is_positive_definite():
        raise ArithmeticError("M25 reduced metric is not positive definite")
    dd = cvp.RoundedMetricCVP(reduced, "dd", PRECISION)
    mpfr = cvp.RoundedMetricCVP(reduced, "mpfr", PRECISION)

    fixed_rows = {
        row["id"]: row for row in visibility["directions"]
        if row["cohort"] == "residual_strict" and row["id"] != "residual-strict-06"
    }
    if len(fixed_rows) != 6:
        raise ArithmeticError("remaining strict panel changed")
    prepared = []
    for direction_id, entry in sorted(fixed_rows.items()):
        old_parity = entry["minimum_finite_coordinate_parity_in_vetted_set"]
        chart = next(row for row in entry["exact_pointed_quartic_rows"] if row["parity"] == old_parity)
        target = tuple(QQ(value) for value in chart["exact_target_point_short_model"])
        if target not in curve:
            raise ArithmeticError("retrospective target is off curve")
        prepared.append((direction_id, entry, chart, target, tuple(Fraction(value) for value in chart["exact_target_point_short_model"])))
    heights, asymmetry = helpers.rounded_gram(curve, (*basis, *(row[3] for row in prepared)))
    if heights[:RANK, :RANK] != gram:
        raise ArithmeticError("rank-25 height Gram differs from the pinned fixed-chart diagnostic")

    directions = []
    for target_index, (direction_id, entry, chart, target, target_for_map) in enumerate(prepared):
        base_parity = [int(value) & 1 for value in chart["centre_m17_word"]]
        target_cross = vector(ZZ, [2 * heights[row, RANK + target_index] for row in range(RANK)])
        projection_original = gram.solve_right(target_cross)
        projection = vector(QQ, (matrix(QQ, 1, RANK, projection_original).row(0) * inverse))
        trials = []
        exact_by_extension = {}
        for extension in range(1 << EXTENSION_DIMENSION):
            original_parity = vector(ZZ, base_parity + [(extension >> bit) & 1 for bit in range(EXTENSION_DIMENSION)])
            reduced_parity = vector(ZZ, (matrix(ZZ, 1, RANK, original_parity) * inverse).row(0))
            reduced_parity = vector(ZZ, [int(value) % 2 for value in reduced_parity])
            u, centre_distance = cvp.cross_precision_closest(dd, mpfr, -vector(QQ, reduced_parity) / 2)
            reduced_centre = reduced_parity + 2 * u
            centre = vector(ZZ, (matrix(ZZ, 1, RANK, reduced_centre) * change).row(0))
            z, target_distance = cvp.cross_precision_closest(
                dd, mpfr, vector(QQ, [(projection[index] - reduced_parity[index]) / 2 for index in range(RANK)])
            )
            translation = vector(ZZ, (matrix(ZZ, 1, RANK, u - z) * change).row(0))
            added = linear_combination(model, map_basis, ints(translation))
            moved = target_for_map if added is None else short_add(model, target_for_map, added)
            if moved is None:
                raise ArithmeticError("target translate became infinity")
            exact = helpers.coordinate(mapper, model, map_basis, centre, moved)
            trial = {
                "extension_mask_in_new_M25_coordinates": extension,
                "M25_parity_mask": parity_mask(centre),
                "centre_M25_word": ints(centre),
                "centre_metric_norm": int(centre * gram * centre),
                "centre_CVP_distance_in_reduced_metric": str(centre_distance),
                "target_translation_M25_word": ints(translation),
                "target_CVP_distance_in_reduced_metric": str(target_distance),
                "reduced_coordinate": exact["reduced_coordinate"],
                "reduced_coordinate_height": exact["reduced_coordinate_height"],
                "reduced_coordinate_decimal_digits": exact["reduced_coordinate_decimal_digits"],
            }
            trials.append(trial)
            exact_by_extension[extension] = exact
        trials.sort(key=lambda row: (int(row["reduced_coordinate_height"]), row["extension_mask_in_new_M25_coordinates"]))
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
            "M25_extension_count": len(trials),
            "M25_extension_winner": {**winner, **exact_by_extension[winner["extension_mask_in_new_M25_coordinates"]]},
            "M25_extension_coordinate_trials": trials,
            "inside_historical_height_125000": int(winner["reduced_coordinate_height"]) <= 125000,
        })
        print("CURVE302M25FIBRE|{}|winner_digits={}|extension={}".format(direction_id, winner["reduced_coordinate_decimal_digits"], winner["extension_mask_in_new_M25_coordinates"]), flush=True)
    return {
        "schema": "elliptic-curves.curve302-m25-remaining-strict-fibre-scan.v1",
        "status": "PASS_RETROSPECTIVE_M25_PARITY_FIBRE_SCAN",
        "inputs": {rel(path): sha(path) for path in (VISIBILITY, M25, M25_CVP, M25_CVP_SOURCE, MAPPER, CVP_SOURCE, Path(__file__))},
        "protocol": {
            "directions": 6,
            "parities_per_direction": 1 << EXTENSION_DIMENSION,
            "fixed_M17_shell_projection": "each direction retains one already-vetted M17 centre parity; all eight newly available M25 parity bits are enumerated",
            "centre_and_target_reduction": "DD/MPFR-agreeing CVP in the pinned rank-25 rounded 384-bit height metric",
            "boundary": "Known target points are retrospective inputs. This does not scan every M25 parity, prove a global minimum coordinate, select a future execution chart, or search for points.",
        },
        "maximum_height_matrix_asymmetry": asymmetry,
        "directions": directions,
        "reproducing_command": "sage -python elliptic-curves/cas/audit_curve302_m25_remaining_strict_fibre_scan.sage --check",
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
            raise FileExistsError("preserve immutable M25 parity-fibre diagnostic")
        OUTPUT.write_text(rendered)
    elif OUTPUT.read_text() != rendered:
        raise ArithmeticError("stored M25 parity-fibre diagnostic did not replay")
    print("CURVE302M25FIBRE|directions=6|parities_per_direction=256|status=PASS", flush=True)


if __name__ == "__main__":
    main()

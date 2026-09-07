#!/usr/bin/env sage-python
"""Exact M28 parity-fibre diagnostic for one remaining curve302 strict label."""

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
M28 = ART / "curve302_m27_fibre_orbit106210_extension673_mod2_v1.json"
HELPERS = CAS / "audit_curve302_m25_remaining_strict_cvp.sage"
MAPPER = CAS / "factor_free_pari_mapping.sage"
CVP = CAS / "audit_curve302_residual_visibility_geometry.sage"
RANK = 28
M17 = 17
EXTENSION_DIMENSION = 11
PRECISION = 192
ALLOWED = {"residual-strict-03", "residual-strict-05", "residual-strict-07"}


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rel(path: Path) -> str:
    return str(path.relative_to(ROOT))


def output(direction: str) -> Path:
    return ART / ("curve302_m28_" + direction.replace("-", "_") + "_fibre_v1.json")


def build(direction: str):
    visibility, m28 = json.loads(VISIBILITY.read_text()), json.loads(M28.read_text())
    if direction not in ALLOWED:
        raise ValueError("not a remaining M28 residual strict direction")
    if m28["status"] != "COMPLETE_DECLARED_FINITE_AUDIT" or m28["rank_lower_bound"] != RANK:
        raise ArithmeticError("rank-28 source cloud changed")
    entry = next(row for row in visibility["directions"] if row["id"] == direction)
    old_parity = entry["minimum_finite_coordinate_parity_in_vetted_set"]
    chart = next(row for row in entry["exact_pointed_quartic_rows"] if row["parity"] == old_parity)
    helpers = importlib.machinery.SourceFileLoader("curve302_m28_helpers", str(HELPERS)).load_module()
    cvp = importlib.machinery.SourceFileLoader("curve302_m28_cvp", str(CVP)).load_module()
    mapper = importlib.machinery.SourceFileLoader("curve302_m28_mapper", str(MAPPER)).load_module()
    from alternate_quartic_covers import short_add
    from half_lattice_pointed_sieve import linear_combination
    model = tuple(Fraction(value) for value in m28["curve"])
    curve = EllipticCurve(QQ, [QQ(value) for value in m28["curve"]])
    basis = tuple(tuple(QQ(value) for value in point) for point in m28["independent_points"])
    map_basis = tuple(tuple(Fraction(value) for value in point) for point in m28["independent_points"])
    target = tuple(QQ(value) for value in chart["exact_target_point_short_model"])
    target_map = tuple(Fraction(value) for value in chart["exact_target_point_short_model"])
    if len(basis) != RANK or target not in curve or any(point not in curve for point in basis):
        raise ArithmeticError("rank-28 points are malformed")
    heights, asymmetry = helpers.rounded_gram(curve, (*basis, target))
    gram = heights[:RANK, :RANK]
    change = matrix(ZZ, pari(gram).qflllgram()).transpose()
    inverse = change.inverse()
    reduced = change * gram * change.transpose()
    dd, mpfr = cvp.RoundedMetricCVP(reduced, "dd", PRECISION), cvp.RoundedMetricCVP(reduced, "mpfr", PRECISION)
    target_cross = vector(ZZ, [2 * heights[row, RANK] for row in range(RANK)])
    projection = vector(QQ, (matrix(QQ, 1, RANK, gram.solve_right(target_cross)).row(0) * inverse))
    base_parity = [int(value) & 1 for value in chart["centre_m17_word"]]
    trials, exact = [], {}
    for extension in range(1 << EXTENSION_DIMENSION):
        original = vector(ZZ, base_parity + [(extension >> bit) & 1 for bit in range(EXTENSION_DIMENSION)])
        parity = vector(ZZ, (matrix(ZZ, 1, RANK, original) * inverse).row(0))
        parity = vector(ZZ, [int(value) % 2 for value in parity])
        u, centre_distance = cvp.cross_precision_closest(dd, mpfr, -vector(QQ, parity) / 2)
        centre = vector(ZZ, (matrix(ZZ, 1, RANK, parity + 2 * u) * change).row(0))
        z, target_distance = cvp.cross_precision_closest(dd, mpfr, vector(QQ, [(projection[index] - parity[index]) / 2 for index in range(RANK)]))
        translation = vector(ZZ, (matrix(ZZ, 1, RANK, u - z) * change).row(0))
        added = linear_combination(model, map_basis, [int(value) for value in translation])
        moved = target_map if added is None else short_add(model, target_map, added)
        if moved is None:
            raise ArithmeticError("still-missing target translate became infinity")
        mapped = helpers.coordinate(mapper, model, map_basis, centre, moved)
        record = {"extension_mask_in_new_M28_coordinates": extension, "M28_parity_mask": sum((int(value) & 1) << index for index, value in enumerate(centre)), "centre_M28_word": [int(value) for value in centre], "centre_metric_norm": int(centre * gram * centre), "centre_CVP_distance_in_reduced_metric": str(centre_distance), "target_translation_M28_word": [int(value) for value in translation], "target_CVP_distance_in_reduced_metric": str(target_distance), "reduced_coordinate": mapped["reduced_coordinate"], "reduced_coordinate_height": mapped["reduced_coordinate_height"], "reduced_coordinate_decimal_digits": mapped["reduced_coordinate_decimal_digits"]}
        trials.append(record)
        exact[extension] = mapped
    trials.sort(key=lambda row: (int(row["reduced_coordinate_height"]), row["extension_mask_in_new_M28_coordinates"]))
    winner = trials[0]
    return {"schema": "elliptic-curves.curve302-m28-one-residual-fibre.v1", "status": "PASS_RETROSPECTIVE_M28_ONE_RESIDUAL_FIBRE_SCAN", "inputs": {rel(path): sha(path) for path in (VISIBILITY, M28, HELPERS, MAPPER, CVP, Path(__file__))}, "protocol": {"direction": direction, "parities": 2048, "fixed_M17_chart": chart["parity"], "boundary": "Known target only; no point search, global M28 scan, or prospective selection claim."}, "maximum_height_matrix_asymmetry": asymmetry, "fixed_M17_chart": {"parity": chart["parity"], "centre_M17_word": chart["centre_m17_word"], "degree_two_translation_orbit": chart["degree_two_translation_orbit"], "old_vetted_coordinate_height": entry["minimum_finite_reduced_coordinate_height_in_vetted_set"]}, "winner": {**winner, **exact[winner["extension_mask_in_new_M28_coordinates"]]}, "coordinate_trials": trials, "inside_historical_height_125000": int(winner["reduced_coordinate_height"]) <= 125000, "reproducing_command": "sage -python elliptic-curves/cas/audit_curve302_m28_one_residual_fibre.sage --check --direction " + direction}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--direction", choices=sorted(ALLOWED), required=True)
    parser.add_argument("--build", action="store_true")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.build == args.check:
        parser.error("choose exactly one of --build or --check")
    payload = build(args.direction)
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    path = output(args.direction)
    if args.build:
        if path.exists():
            raise FileExistsError("preserve immutable M28 one-direction diagnostic")
        path.write_text(rendered)
    elif path.read_text() != rendered:
        raise ArithmeticError("stored M28 diagnostic did not replay")
    print("CURVE302M28ONE|{}|parities=2048|winner_digits={}".format(args.direction, payload["winner"]["reduced_coordinate_decimal_digits"]), flush=True)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Frozen one-chart M28 arm: norm-ten orbit4761, extension1772."""

from __future__ import annotations

import argparse
import copy
from pathlib import Path

import certify_compact_r17_candidates as cert
import run_curve302_m27_fibre_orbit106210_extension673 as imported


ROOT = Path(__file__).resolve().parents[2]
ART = ROOT / "artifacts/generated-results/elliptic-curves"
SEED = ART / "curve302_m27_fibre_orbit106210_extension673_mod2_v1.json"
LOCAL = ROOT / "artifacts/local/elliptic-curves/curve302-m28-fibre-orbit4761-extension1772-v1"
OUTPUT = ART / "curve302_m28_fibre_orbit4761_extension1772_v1.json"
MOD2 = ART / "curve302_m28_fibre_orbit4761_extension1772_mod2_v1.json"
MODL = ART / "curve302_m28_fibre_orbit4761_extension1772_modl_v1.json"
RANK, PARITY, EXTENSION, ORBIT, PARITY_MASK = 28, 8515, 1772, 4761, 232268099
CENTER = [-5, 3, 10, 0, 0, 4, -5, 6, 1, 4, -8, 0, 0, -3, -4, 0, 0, 2, 2, -1, -1, 0, -1, 1, -1, 2, 1, -1]


base = imported.base
middle = imported.template
inner = middle.template
for module in (middle,):
    module.SEED_CLOUD, module.LOCAL, module.MOD2, module.MODL = SEED, LOCAL, MOD2, MODL
    module.M17_PARITY, module.EXTENSION_MASK, module.INITIAL_RANK = PARITY, EXTENSION, RANK
    module.EXPECTED_CENTRE, module.EXPECTED_PARITY = CENTER, PARITY_MASK
inner.SEED_CLOUD, inner.INITIAL_RANK = SEED, RANK
base.M24, base.LOCAL, base.OUTPUT, base.MOD2, base.MODL, base.M24_DIMENSION, base.ORBIT_MASK = SEED, LOCAL, OUTPUT, MOD2, MODL, RANK, ORBIT


def hashes():
    result = inner.BASE_SOURCES()
    for path in (Path(inner.__file__), Path(middle.__file__), Path(imported.__file__), Path(__file__)):
        result[base.rel(path)] = cert.hashed(path)
    return result


def centre():
    row = base.orbit_row()
    if row["orbit_mask"] != ORBIT or row["minimum_generic_MW17_norm"] != 10 or row["category"] != "rational":
        raise ArithmeticError("frozen norm-ten orbit changed")
    if len(CENTER) != RANK or sum((value & 1) << index for index, value in enumerate(CENTER)) != PARITY_MASK:
        raise ArithmeticError("frozen M28 parity changed")
    return {"base_M17_parity_mask": PARITY, "new_M25_extension_mask": EXTENSION, "M25_parity_mask": PARITY_MASK, "representative": CENTER, "selection_rule": "metric-shortest representative in the fixed M28 parity, rederived in Sage before mapping", "generic_degree_two_orbit": row}


def protocol():
    result = copy.deepcopy(inner.BASE_PROTOCOL())
    result.update(schema="elliptic-curves.curve302-m28-fibre-orbit4761-extension1772.v1", initial_rank=RANK, centre_selection={"rule": "fix M17 parity8515 (generic rational norm-ten degree-two orbit4761), extend it by M28 bits1772, then select the shortest representative in the rounded specialized M28 metric", "base_M17_parity_mask": PARITY, "frozen_M28_extension_mask": EXTENSION, "generic_degree_two_orbit_mask": ORBIT, "calibration_status": "retrospectively calibrated; execution-blind", "minimum_change": "one additional M28 centre/chart", "selection_inputs": [base.rel(path) for path in (base.PARENT, base.LATTICE, base.ORBITS, SEED)]}, execution_blindness={"allowed_artifact_inputs": [base.rel(path) for path in (SEED, base.PARENT, base.LATTICE, base.ORBITS)], "forbidden_inputs": ["the residual-visibility diagnostics", "all public missing-point coordinates and public rank-31 point lists", "post-search finite-rank outcomes"], "enforcement": "geometry, worker, and replay install an artifact-read guard before arithmetic imports"}, boundary="One retrospectively calibrated, execution-blind M28 chart; only a completed-cloud lower bound may result.")
    return result


def freeze(): middle.freeze()
def geometry(): middle.geometry()
def worker(): middle.worker()
def replay(): middle.replay()
def certify():
    middle.certify()
    if base.read(MOD2)["rank_lower_bound"] < 29 or any(row["finite_column_rank"] < 29 for row in base.read(MODL)["audits"]):
        raise ArithmeticError("frozen M28 arm did not certify rank at least 29")


inner.source_hashes = hashes
middle.source_hashes = hashes
base.source_hashes = hashes
inner.frozen_centre = centre
middle.frozen_centre = centre
base.centre_from_orbit = centre
base.protocol = protocol


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("freeze", "geometry", "worker", "replay", "certify"))
    args = parser.parse_args()
    globals()[args.stage]()

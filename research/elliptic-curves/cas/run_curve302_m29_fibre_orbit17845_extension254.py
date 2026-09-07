#!/usr/bin/env python3
"""Frozen one-chart M29 arm: nearby norm-eight orbit17845, extension254."""

from __future__ import annotations

import argparse
import copy
from pathlib import Path

import certify_compact_r17_candidates as cert
import run_curve302_m28_fibre_orbit4761_extension1772 as template


ROOT = Path(__file__).resolve().parents[2]
ART = ROOT / "artifacts/generated-results/elliptic-curves"
SEED = ART / "curve302_m28_fibre_orbit4761_extension1772_mod2_v1.json"
LOCAL = ROOT / "artifacts/local/elliptic-curves/curve302-m29-fibre-orbit17845-extension254-v1"
OUTPUT = ART / "curve302_m29_fibre_orbit17845_extension254_v1.json"
MOD2 = ART / "curve302_m29_fibre_orbit17845_extension254_mod2_v1.json"
MODL = ART / "curve302_m29_fibre_orbit17845_extension254_modl_v1.json"
RANK, PARITY, EXTENSION, ORBIT, PARITY_MASK = 29, 42625, 254, 17845, 33334913
CENTER = [-7, 4, 12, 2, 0, 6, -4, 7, 2, 5, -11, 0, 2, -7, -6, 1, 2, 2, 3, -1, -1, -1, -1, 1, -1, 2, 0, 0, 0]


base, middle, inner = template.base, template.middle, template.inner
middle.SEED_CLOUD, middle.LOCAL, middle.MOD2, middle.MODL = SEED, LOCAL, MOD2, MODL
middle.M17_PARITY, middle.EXTENSION_MASK, middle.INITIAL_RANK = PARITY, EXTENSION, RANK
middle.EXPECTED_CENTRE, middle.EXPECTED_PARITY = CENTER, PARITY_MASK
inner.SEED_CLOUD, inner.INITIAL_RANK = SEED, RANK
base.M24, base.LOCAL, base.OUTPUT, base.MOD2, base.MODL, base.M24_DIMENSION, base.ORBIT_MASK = SEED, LOCAL, OUTPUT, MOD2, MODL, RANK, ORBIT


def hashes():
    result = inner.BASE_SOURCES()
    for path in (Path(inner.__file__), Path(middle.__file__), Path(template.__file__), Path(__file__)):
        result[base.rel(path)] = cert.hashed(path)
    return result


def nearby_orbit():
    lattice = base.read(base.LATTICE)
    if lattice["status"] != "PASS_COMPLETE_DEGREE2_TRANSLATION_QUOTIENT" or base.sha(base.ORBITS) != lattice["orbits_tsv_sha256"]:
        raise ArithmeticError("complete degree-two table changed")
    row = next(line.split("\t") for line in base.ORBITS.read_text().splitlines()[1:] if int(line.split("\t", 1)[0]) == ORBIT)
    if row[1] != "genus_one" or int(row[2]) != 8:
        raise ArithmeticError("frozen nearby-shell orbit changed")
    return {"orbit_mask": ORBIT, "category": "genus_one_bisection_lattice_candidate", "minimum_generic_MW17_norm": 8, "degree_two_LLL_minimum_word": [int(value) for value in row[3].split()], "parent_MW17_representative": [int(value) for value in row[4].split()], "degree_two_divisor_class": row[5]}


def centre():
    orbit = nearby_orbit()
    if len(CENTER) != RANK or sum((value & 1) << index for index, value in enumerate(CENTER)) != PARITY_MASK:
        raise ArithmeticError("frozen M29 parity changed")
    return {"base_M17_parity_mask": PARITY, "new_M25_extension_mask": EXTENSION, "M25_parity_mask": PARITY_MASK, "representative": CENTER, "selection_rule": "metric-shortest representative in the fixed M29 parity, rederived in Sage before mapping", "generic_degree_two_orbit": orbit}


def protocol():
    result = copy.deepcopy(inner.BASE_PROTOCOL())
    result.update(schema="elliptic-curves.curve302-m29-fibre-orbit17845-extension254.v1", initial_rank=RANK, centre_selection={"rule": "fix M17 parity42625 (nearby generic norm-eight degree-two orbit17845), extend it by M29 bits254, then select the shortest representative in the rounded specialized M29 metric", "base_M17_parity_mask": PARITY, "frozen_M29_extension_mask": EXTENSION, "generic_degree_two_orbit_mask": ORBIT, "calibration_status": "retrospectively calibrated; execution-blind", "minimum_change": "one additional M29 centre/chart", "selection_inputs": [base.rel(path) for path in (base.PARENT, base.LATTICE, base.ORBITS, SEED)]}, execution_blindness={"allowed_artifact_inputs": [base.rel(path) for path in (SEED, base.PARENT, base.LATTICE, base.ORBITS)], "forbidden_inputs": ["the residual-visibility diagnostics", "all public missing-point coordinates and public rank-31 point lists", "post-search finite-rank outcomes"], "enforcement": "geometry, worker, and replay install an artifact-read guard before arithmetic imports"}, boundary="One retrospectively calibrated, execution-blind M29 chart; only a completed-cloud lower bound may result.")
    return result


def freeze(): middle.freeze()
def geometry(): middle.geometry()
def worker(): middle.worker()
def replay(): middle.replay()
def certify():
    middle.certify()
    if base.read(MOD2)["rank_lower_bound"] < 30 or any(row["finite_column_rank"] < 30 for row in base.read(MODL)["audits"]):
        raise ArithmeticError("frozen M29 arm did not certify rank at least 30")


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

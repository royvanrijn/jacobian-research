#!/usr/bin/env python3
"""Final frozen M30 arm: pad the execution-blind M29 norm-eight fibre centre."""

from __future__ import annotations

import argparse
import copy
from pathlib import Path

import certify_compact_r17_candidates as cert
import run_curve302_m24_bisection_orbit117420 as base


ROOT = Path(__file__).resolve().parents[2]
ART = ROOT / "artifacts/generated-results/elliptic-curves"
SEED = ART / "curve302_m29_fibre_orbit17845_extension254_mod2_v1.json"
LOCAL = ROOT / "artifacts/local/elliptic-curves/curve302-m30-padded-fibre-orbit114326-extension1502-v1"
OUTPUT = ART / "curve302_m30_padded_fibre_orbit114326_extension1502_v1.json"
MOD2 = ART / "curve302_m30_padded_fibre_orbit114326_extension1502_mod2_v1.json"
MODL = ART / "curve302_m30_padded_fibre_orbit114326_extension1502_modl_v1.json"
RANK, M17_PARITY, ORBIT, M29_EXTENSION = 30, 123536, 114326, 1502
CENTER = [6, -4, -12, -2, 1, -6, 4, -7, -2, -5, 10, 0, -2, 7, 7, -1, -1, -2, -3, 1, 1, 1, 0, -1, 1, -1, 0, 1, 0, 0]
PARITY_MASK = 196993680


base.M24, base.LOCAL, base.OUTPUT, base.MOD2, base.MODL, base.M24_DIMENSION, base.ORBIT_MASK = SEED, LOCAL, OUTPUT, MOD2, MODL, RANK, ORBIT
BASE_SOURCES, BASE_PROTOCOL = base.source_hashes, base.protocol


def hashes():
    result = BASE_SOURCES()
    result[base.rel(Path(__file__))] = cert.hashed(Path(__file__))
    return result


def seed_from_m30():
    cloud = base.read(SEED)
    if cloud["status"] != "COMPLETE_DECLARED_FINITE_AUDIT" or cloud["rank_lower_bound"] != RANK:
        raise ArithmeticError("rank-30 source cloud changed")
    points, proof = cloud["independent_points"], cloud["rank_certificate"]
    if len(points) != RANK:
        raise ArithmeticError("rank-30 source cloud lacks its basis")
    actual = base.checked_rank(tuple(map(cert.F, cloud["curve"])), [tuple(map(cert.F, point)) for point in points], [signature["prime"] for signature in proof["signatures"]], proof["no_rational_2_torsion_prime"])
    if base.digest(actual) != base.digest(proof):
        raise ArithmeticError("rank-30 seed certificate did not replay")
    return {"family": cloud["family"], "parameter": cloud["parameter"], "curve": cloud["curve"], "points": points, "generic_points": points[:17], "rank_certificate": proof, "input_cloud": base.rel(SEED), "input_cloud_sha256": base.sha(SEED)}


def nearby_orbit():
    lattice = base.read(base.LATTICE)
    if lattice["status"] != "PASS_COMPLETE_DEGREE2_TRANSLATION_QUOTIENT" or base.sha(base.ORBITS) != lattice["orbits_tsv_sha256"]:
        raise ArithmeticError("complete degree-two table changed")
    row = next(line.split("\t") for line in base.ORBITS.read_text().splitlines()[1:] if int(line.split("\t", 1)[0]) == ORBIT)
    if row[1] != "genus_one" or int(row[2]) != 8:
        raise ArithmeticError("frozen nearby-shell orbit changed")
    return {"orbit_mask": ORBIT, "category": "genus_one_bisection_lattice_candidate", "minimum_generic_MW17_norm": 8, "degree_two_LLL_minimum_word": [int(value) for value in row[3].split()], "parent_MW17_representative": [int(value) for value in row[4].split()], "degree_two_divisor_class": row[5]}


def centre():
    if len(CENTER) != RANK or sum((value & 1) << index for index, value in enumerate(CENTER)) != PARITY_MASK:
        raise ArithmeticError("padded M30 centre parity changed")
    return {"base_M17_parity_mask": M17_PARITY, "M29_extension_mask": M29_EXTENSION, "M30_padding_coordinate": 0, "M30_parity_mask": PARITY_MASK, "representative": CENTER, "selection_rule": "the fixed M29 parity-fibre representative, padded by zero in the newly certified M30 coordinate", "generic_degree_two_orbit": nearby_orbit()}


def protocol():
    result = copy.deepcopy(BASE_PROTOCOL())
    result.update(schema="elliptic-curves.curve302-m30-padded-fibre-orbit114326-extension1502.v1", initial_rank=RANK, centre_selection={"rule": "fix the prior M29 parity-fibre centre over M17 parity123536 (nearby norm-eight degree-two orbit114326) with M29 extension1502, then append zero in the M30 coordinate", "base_M17_parity_mask": M17_PARITY, "frozen_M29_extension_mask": M29_EXTENSION, "M30_padding_coordinate": 0, "generic_degree_two_orbit_mask": ORBIT, "calibration_status": "retrospectively calibrated; execution-blind", "minimum_change": "one additional M30 chart", "selection_inputs": [base.rel(path) for path in (base.PARENT, base.LATTICE, base.ORBITS, SEED)]}, execution_blindness={"allowed_artifact_inputs": [base.rel(path) for path in (SEED, base.PARENT, base.LATTICE, base.ORBITS)], "forbidden_inputs": ["the residual-visibility diagnostics", "all public missing-point coordinates and public rank-31 point lists", "post-search finite-rank outcomes"], "enforcement": "geometry, worker, and replay install an artifact-read guard before arithmetic imports"}, boundary="One retrospectively calibrated, execution-blind M30 chart; only a completed-cloud lower bound may result.")
    return result


def freeze(): base.freeze()
def geometry(): base.geometry()
def worker(): base.worker()
def replay(): base.replay()
def certify():
    base.certify()
    if base.read(MOD2)["rank_lower_bound"] < 31 or any(row["finite_column_rank"] < 31 for row in base.read(MODL)["audits"]):
        raise ArithmeticError("final frozen M30 arm did not certify rank at least 31")


base.source_hashes = hashes
base.protocol = protocol
base.seed_from_m24 = seed_from_m30
base.centre_from_orbit = centre


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("freeze", "geometry", "worker", "replay", "certify"))
    args = parser.parse_args()
    globals()[args.stage]()

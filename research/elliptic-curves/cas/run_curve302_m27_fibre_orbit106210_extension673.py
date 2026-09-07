#!/usr/bin/env python3
"""One frozen M27 parity-fibre arm for curve302 residual strict-01."""

from __future__ import annotations

import argparse
import copy
import json
import sys
from pathlib import Path

import certify_compact_r17_candidates as cert
import run_curve302_m26_fibre_orbit58145_extension472 as template


ROOT = Path(__file__).resolve().parents[2]
ART = ROOT / "artifacts/generated-results/elliptic-curves"
SEED_CLOUD = ART / "curve302_m26_fibre_orbit58145_extension472_mod2_v1.json"
LOCAL = ROOT / "artifacts/local/elliptic-curves/curve302-m27-fibre-orbit106210-extension673-v1"
OUTPUT = ART / "curve302_m27_fibre_orbit106210_extension673_v1.json"
MOD2 = ART / "curve302_m27_fibre_orbit106210_extension673_mod2_v1.json"
MODL = ART / "curve302_m27_fibre_orbit106210_extension673_modl_v1.json"
M17_PARITY = 103240
EXTENSION_MASK = 673
INITIAL_RANK = 27
EXPECTED_CENTRE = [-6, 4, 10, 1, 0, 4, -5, 6, 1, 3, -8, 0, 1, -4, -4, 1, 1, 3, 2, -2, 0, 0, -1, 0, -1, 2, 1]
EXPECTED_PARITY = 88314696


base = template.base
inner = template.template
template.SEED_CLOUD = SEED_CLOUD
template.LOCAL = LOCAL
template.OUTPUT = OUTPUT
template.MOD2 = MOD2
template.MODL = MODL
template.M17_PARITY = M17_PARITY
template.EXTENSION_MASK = EXTENSION_MASK
template.INITIAL_RANK = INITIAL_RANK
template.EXPECTED_CENTRE = EXPECTED_CENTRE
template.EXPECTED_PARITY = EXPECTED_PARITY
inner.SEED_CLOUD = SEED_CLOUD
inner.INITIAL_RANK = INITIAL_RANK
base.M24 = SEED_CLOUD
base.LOCAL = LOCAL
base.OUTPUT = OUTPUT
base.MOD2 = MOD2
base.MODL = MODL
base.M24_DIMENSION = INITIAL_RANK
base.ORBIT_MASK = 106210


def source_hashes() -> dict[str, str]:
    result = inner.BASE_SOURCES()
    for path in (Path(inner.__file__), Path(template.__file__), Path(__file__)):
        result[base.rel(path)] = cert.hashed(path)
    return result


def frozen_centre() -> dict:
    orbit = base.orbit_row()
    if orbit["orbit_mask"] != 106210 or orbit["minimum_generic_MW17_norm"] != 10 or orbit["category"] != "rational":
        raise ArithmeticError("frozen parent bisection orbit changed")
    if len(EXPECTED_CENTRE) != INITIAL_RANK or sum((value & 1) << index for index, value in enumerate(EXPECTED_CENTRE)) != EXPECTED_PARITY:
        raise ArithmeticError("frozen M27 centre parity changed")
    return {"base_M17_parity_mask": M17_PARITY, "new_M25_extension_mask": EXTENSION_MASK, "M25_parity_mask": EXPECTED_PARITY, "representative": EXPECTED_CENTRE, "selection_rule": "metric-shortest representative in the fixed M27 parity; the stored word is rederived in the Sage geometry stage", "generic_degree_two_orbit": orbit}


def protocol() -> dict:
    result = copy.deepcopy(inner.BASE_PROTOCOL())
    result.update(schema="elliptic-curves.curve302-m27-fibre-orbit106210-extension673.v1", initial_rank=INITIAL_RANK, centre_selection={"rule": "fix M17 parity103240 (generic rational norm-ten degree-two orbit106210), extend it by M27 bits673, then select the shortest representative in the rounded specialized M27 metric", "base_M17_parity_mask": M17_PARITY, "frozen_M27_extension_mask": EXTENSION_MASK, "generic_degree_two_orbit_mask": 106210, "calibration_status": "retrospectively calibrated; execution-blind", "minimum_change": "one additional M27 centre/chart; earlier sealed policies remain unmodified", "selection_inputs": [base.rel(path) for path in (base.PARENT, base.LATTICE, base.ORBITS, SEED_CLOUD)]}, execution_blindness={"allowed_artifact_inputs": [base.rel(path) for path in (SEED_CLOUD, base.PARENT, base.LATTICE, base.ORBITS)], "forbidden_inputs": ["the residual-visibility diagnostics", "all public missing-point coordinates and public rank-31 point lists", "post-search finite-rank outcomes"], "enforcement": "geometry, worker, and replay install an artifact-read guard before their arithmetic imports and serialize accepted read paths"}, boundary="A single retrospectively calibrated target-blind M27 parity-fibre arm. It can certify only the completed point cloud's lower bound; it is neither a prospective ranking validation nor an exact-rank proof.")
    return result


def freeze() -> None:
    template.freeze()


def geometry() -> None:
    template.geometry()


def worker() -> None:
    template.worker()


def replay() -> None:
    template.replay()


def certify() -> None:
    template.certify()
    mod2, modl = base.read(MOD2), base.read(MODL)
    if mod2["rank_lower_bound"] < 28 or any(audit["finite_column_rank"] < 28 for audit in modl["audits"]):
        raise ArithmeticError("the single frozen M27 parity-fibre arm did not certify rank at least 28")


def report() -> None:
    p = base.assert_protocol()
    if OUTPUT.exists():
        raise FileExistsError("preserve immutable M27 fibre report")
    names = ("seed.json", "centre.json", "maps.json", "result.json", "replay.json", "certification-ledger.json")
    local = [LOCAL / name for name in names]
    if any(not path.exists() for path in local) or not MOD2.exists() or not MODL.exists():
        raise ArithmeticError("complete M27 arm evidence is required")
    centre, maps, result, replay, certificate = (base.read(path) for path in local[1:])
    mod2, modl, old = base.read(MOD2), base.read(MODL), base.read(base.OLD_M17_MAPS)
    if any(row["parity"] == M17_PARITY for row in old["sample"]) or any(row["parity"] == M17_PARITY for row in old["centres"]):
        raise ArithmeticError("historical deep arm unexpectedly contains the M17 base parity")
    if certificate["status"] != "PASS" or mod2["rank_lower_bound"] < 28 or any(audit["finite_column_rank"] < 28 for audit in modl["audits"]):
        raise ArithmeticError("rank-28 certificate missing")
    payload = {"schema": "elliptic-curves.curve302-m27-fibre-orbit106210-extension673-report.v1", "status": "PASS_CERTIFIED_RANK_AT_LEAST_28", "inputs": {**base.frozen_inputs(), base.rel(base.OLD_M17_MAPS): base.sha(base.OLD_M17_MAPS), **{base.rel(path): base.sha(path) for path in (*local, MOD2, MODL)}}, "sources": source_hashes(), "frozen_rule": p["centre_selection"], "exact_M27_centre": {**centre, "centre_point_short_model": result["search"]["base_point"], "reduced_quartic_coefficients_ascending": result["search"]["coefficients"], "coordinate_matrix": maps["mapping"]["matrix"]}, "execution_oracle_audit": {"rule_is_retrospectively_calibrated": True, "execution_is_target_blind": True, "allowed_inputs": p["execution_blindness"]["allowed_artifact_inputs"], "geometry_reads": base.read(LOCAL / "geometry-data-access.json"), "worker_reads": base.read(LOCAL / "worker-data-access.json"), "replay_reads": base.read(LOCAL / "replay-data-access.json"), "excluded_from_execution": p["execution_blindness"]["forbidden_inputs"]}, "completed_search": {"rank_before": INITIAL_RANK, "rank_after_worker_mod2_admission": result["rank_lower_bound"], "rank_after_exact_map_replay": replay["rank_lower_bound"], "mod_2_certified_rank_lower_bound": mod2["rank_lower_bound"], "mod_3_5_ranks": {str(audit["modulus"]): audit["finite_column_rank"] for audit in modl["audits"]}, "height": p["height"], "seconds_per_chart": p["seconds_per_chart"], "search_status": result["search"]["status"], "returned_finite_curve_points": result["search"]["finite_curve_points"]}, "boundary": "This execution was target-blind but its single rule was retrospectively calibrated. It proves a lower bound for the completed cloud only.", "reproducing_command": "python3 elliptic-curves/cas/run_curve302_m27_fibre_orbit106210_extension673.py launch"}
    OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print("REPORTED CURVE302 M27 FIBRE|rank>={}".format(mod2["rank_lower_bound"]), flush=True)


def launch() -> None:
    freeze()
    ledger = {"schema": "elliptic-curves.curve302-m27-fibre-orbit106210-extension673-launch.v1", "status": "RUNNING", "stages": []}
    base.checkpoint(LOCAL / "ledger.json", ledger)
    for name, command, seconds in (("geometry", [base.SAGE, str(Path(__file__).resolve()), "geometry"], 180), ("worker", [sys.executable, str(Path(__file__).resolve()), "worker"], 120), ("replay", [sys.executable, str(Path(__file__).resolve()), "replay"], 120), ("certify", [sys.executable, str(Path(__file__).resolve()), "certify"], 900), ("report", [sys.executable, str(Path(__file__).resolve()), "report"], 120)):
        base.launch_stage(name, command, seconds, ledger)
    ledger.update(status="PASS", report=base.rel(OUTPUT), report_sha256=base.sha(OUTPUT))
    base.checkpoint(LOCAL / "ledger.json", ledger)


template.source_hashes = source_hashes
inner.source_hashes = source_hashes
base.source_hashes = source_hashes
template.frozen_centre = frozen_centre
inner.frozen_centre = frozen_centre
base.centre_from_orbit = frozen_centre
base.protocol = protocol


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("freeze", "geometry", "worker", "replay", "certify", "report", "launch"))
    args = parser.parse_args()
    globals()[args.stage]()

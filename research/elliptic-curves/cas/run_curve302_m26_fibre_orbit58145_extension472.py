#!/usr/bin/env python3
"""One frozen M26 parity-fibre arm for curve302 residual strict-02.

The predeclared target-blind rule fixes M17 parity28399 (parent rational
norm-ten degree-two orbit58145), fixes its nine-bit M26 extension472, and
derives the metric-shortest centre from the rank-26 seed alone.
"""

from __future__ import annotations

import argparse
import copy
import json
import sys
from decimal import Decimal, getcontext
from importlib.machinery import SourceFileLoader
from pathlib import Path

import certify_compact_r17_candidates as cert
import run_curve302_m25_fibre_orbit64677_extension103 as template


ROOT = Path(__file__).resolve().parents[2]
ART = ROOT / "artifacts/generated-results/elliptic-curves"
SEED_CLOUD = ART / "curve302_m25_fibre_orbit64677_extension103_mod2_v1.json"
LOCAL = ROOT / "artifacts/local/elliptic-curves/curve302-m26-fibre-orbit58145-extension472-v1"
OUTPUT = ART / "curve302_m26_fibre_orbit58145_extension472_v1.json"
MOD2 = ART / "curve302_m26_fibre_orbit58145_extension472_mod2_v1.json"
MODL = ART / "curve302_m26_fibre_orbit58145_extension472_modl_v1.json"
M17_PARITY = 28399
EXTENSION_MASK = 472
INITIAL_RANK = 26
EXPECTED_CENTRE = [5, -3, -7, -3, 0, -5, 1, -3, -2, -3, 7, -1, -2, 7, 5, -2, 0, 0, -2, 0, 1, 1, 0, -1, 1, -1]
EXPECTED_PARITY = 61894383


# The preceding arm supplies the generic fixed-map/certificate machinery.  It
# is re-bound here to a new seed and a new frozen parity, never to a target.
base = template.base
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
base.M24 = SEED_CLOUD
base.LOCAL = LOCAL
base.OUTPUT = OUTPUT
base.MOD2 = MOD2
base.MODL = MODL
base.M24_DIMENSION = INITIAL_RANK
base.ORBIT_MASK = 58145


def source_hashes() -> dict[str, str]:
    result = template.BASE_SOURCES()
    result[base.rel(Path(template.__file__))] = cert.hashed(Path(template.__file__))
    result[base.rel(Path(__file__))] = cert.hashed(Path(__file__))
    return result


def frozen_centre() -> dict:
    orbit = base.orbit_row()
    if orbit["orbit_mask"] != 58145 or orbit["minimum_generic_MW17_norm"] != 10 or orbit["category"] != "rational":
        raise ArithmeticError("frozen parent bisection orbit changed")
    if len(EXPECTED_CENTRE) != INITIAL_RANK or sum((value & 1) << index for index, value in enumerate(EXPECTED_CENTRE)) != EXPECTED_PARITY:
        raise ArithmeticError("frozen M26 centre parity changed")
    return {"base_M17_parity_mask": M17_PARITY, "new_M25_extension_mask": EXTENSION_MASK, "M25_parity_mask": EXPECTED_PARITY, "representative": EXPECTED_CENTRE, "selection_rule": "metric-shortest representative in the fixed M26 parity; the stored word is rederived in the Sage geometry stage", "generic_degree_two_orbit": orbit}


def protocol() -> dict:
    result = copy.deepcopy(template.BASE_PROTOCOL())
    result.update(
        schema="elliptic-curves.curve302-m26-fibre-orbit58145-extension472.v1",
        initial_rank=INITIAL_RANK,
        centre_selection={"rule": "fix M17 parity28399 (generic rational norm-ten degree-two orbit58145), extend it by M26 bits472, then select the shortest representative in the rounded specialized M26 metric", "base_M17_parity_mask": M17_PARITY, "frozen_M26_extension_mask": EXTENSION_MASK, "generic_degree_two_orbit_mask": 58145, "calibration_status": "retrospectively calibrated; execution-blind", "minimum_change": "one additional M26 centre/chart; earlier sealed policies remain unmodified", "selection_inputs": [base.rel(path) for path in (base.PARENT, base.LATTICE, base.ORBITS, SEED_CLOUD)]},
        execution_blindness={"allowed_artifact_inputs": [base.rel(path) for path in (SEED_CLOUD, base.PARENT, base.LATTICE, base.ORBITS)], "forbidden_inputs": ["the residual-visibility diagnostics", "all public missing-point coordinates and public rank-31 point lists", "post-search finite-rank outcomes"], "enforcement": "geometry, worker, and replay install an artifact-read guard before their arithmetic imports and serialize accepted read paths"},
        boundary="A single retrospectively calibrated target-blind M26 parity-fibre arm. It can certify only the completed point cloud's lower bound; it is neither a prospective ranking validation nor an exact-rank proof.",
    )
    return result


def freeze() -> None:
    template.freeze()


def geometry() -> None:
    base.install_execution_guard()
    p = base.assert_protocol()
    out = LOCAL / "maps.json"
    if out.exists():
        raise FileExistsError("preserve frozen M26 fibre map")
    seed, centre = base.read(LOCAL / "seed.json"), base.read(LOCAL / "centre.json")
    if base.parent_prefix_in_short_model() != seed["generic_points"]:
        raise ArithmeticError("rank-26 prefix is not the certified M17 parent specialization")
    from sage.all import EllipticCurve, QQ, ZZ, matrix, pari, vector

    cvp = SourceFileLoader("curve302_m26_execution_cvp", str(ROOT / "elliptic-curves/cas/audit_curve302_residual_visibility_geometry.sage")).load_module()
    mapper = SourceFileLoader("curve302_m26_execution_mapper", str(ROOT / "elliptic-curves/cas/factor_free_pari_mapping.sage")).load_module()
    model = tuple(map(cert.F, seed["curve"]))
    points = tuple(tuple(map(cert.F, point)) for point in seed["points"])
    getcontext().prec = 110
    pari.default("realprecision", 110)
    pari.allocatemem(256_000_000, 1_073_741_824, silent=True)
    height = pari(EllipticCurve(QQ, [QQ(value) for value in seed["curve"]])).ellheightmatrix([list(point) for point in points], precision=384)
    decimal = [[Decimal(str(height[i, j])) for j in range(INITIAL_RANK)] for i in range(INITIAL_RANK)]
    if max(abs(decimal[i][j] - decimal[j][i]) for i in range(INITIAL_RANK) for j in range(INITIAL_RANK)) > Decimal("1e-90"):
        raise ArithmeticError("M26 height matrix asymmetry exceeds guard")
    gram = matrix(ZZ, INITIAL_RANK, INITIAL_RANK, [int((entry * Decimal(1_000_000)).to_integral_value()) for row in decimal for entry in row])
    change = matrix(ZZ, pari(gram).qflllgram()).transpose()
    if abs(change.det()) != 1:
        raise ArithmeticError("M26 LLL change is not unimodular")
    inverse = change.inverse()
    original_parity = vector(ZZ, [(M17_PARITY >> index) & 1 for index in range(17)] + [(EXTENSION_MASK >> bit) & 1 for bit in range(INITIAL_RANK - 17)])
    reduced_parity = vector(ZZ, (matrix(ZZ, 1, INITIAL_RANK, original_parity) * inverse).row(0))
    reduced_parity = vector(ZZ, [int(value) % 2 for value in reduced_parity])
    reduced_gram = change * gram * change.transpose()
    dd = cvp.RoundedMetricCVP(reduced_gram, "dd", 192)
    mpfr = cvp.RoundedMetricCVP(reduced_gram, "mpfr", 192)
    u, distance = cvp.cross_precision_closest(dd, mpfr, -vector(reduced_parity) / 2)
    derived = vector(ZZ, (matrix(ZZ, 1, INITIAL_RANK, reduced_parity + 2 * u) * change).row(0))
    if list(map(int, derived)) != centre["representative"] or sum((int(value) & 1) << index for index, value in enumerate(derived)) != EXPECTED_PARITY:
        raise ArithmeticError("frozen parity no longer derives the pinned M26 centre")
    mapping = mapper.mapping(model, points, centre)
    base.checkpoint(out, {"schema": "elliptic-curves.curve302-m26-fibre-map.v1", "status": "COMPLETE_DECLARED_SINGLE_MAP", "protocol_sha256": base.sha(LOCAL / "protocol.json"), "seed_sha256": base.sha(LOCAL / "seed.json"), "centre_sha256": base.sha(LOCAL / "centre.json"), "metric_gram": [list(map(int, row)) for row in gram.rows()], "change_of_basis": [list(map(int, row)) for row in change.rows()], "centre_CVP_distance": str(distance), "centre_derived_from_frozen_parity": True, "mapping": mapping})
    base.checkpoint(LOCAL / "geometry-data-access.json", sorted(base.READS))
    print("FROZEN CURVE302 M26 FIBRE MAP|orbit=58145|extension=472", flush=True)


def worker() -> None:
    template.worker()


def replay() -> None:
    template.replay()


def certify() -> None:
    template.certify()
    mod2, modl = base.read(MOD2), base.read(MODL)
    if mod2["rank_lower_bound"] < 27 or any(audit["finite_column_rank"] < 27 for audit in modl["audits"]):
        raise ArithmeticError("the single frozen M26 parity-fibre arm did not certify rank at least 27")


def report() -> None:
    p = base.assert_protocol()
    if OUTPUT.exists():
        raise FileExistsError("preserve immutable M26 fibre report")
    names = ("seed.json", "centre.json", "maps.json", "result.json", "replay.json", "certification-ledger.json")
    local = [LOCAL / name for name in names]
    if any(not path.exists() for path in local) or not MOD2.exists() or not MODL.exists():
        raise ArithmeticError("complete M26 arm evidence is required")
    centre, maps, result, replay, certificate = (base.read(path) for path in local[1:])
    mod2, modl, old = base.read(MOD2), base.read(MODL), base.read(base.OLD_M17_MAPS)
    if any(row["parity"] == M17_PARITY for row in old["sample"]) or any(row["parity"] == M17_PARITY for row in old["centres"]):
        raise ArithmeticError("historical deep arm unexpectedly contains the M17 base parity")
    if certificate["status"] != "PASS" or mod2["rank_lower_bound"] < 27 or any(audit["finite_column_rank"] < 27 for audit in modl["audits"]):
        raise ArithmeticError("rank-27 certificate missing")
    payload = {
        "schema": "elliptic-curves.curve302-m26-fibre-orbit58145-extension472-report.v1",
        "status": "PASS_CERTIFIED_RANK_AT_LEAST_27",
        "inputs": {**base.frozen_inputs(), base.rel(base.OLD_M17_MAPS): base.sha(base.OLD_M17_MAPS), **{base.rel(path): base.sha(path) for path in (*local, MOD2, MODL)}},
        "sources": source_hashes(),
        "frozen_rule": p["centre_selection"],
        "exact_M26_centre": {**centre, "centre_point_short_model": result["search"]["base_point"], "reduced_quartic_coefficients_ascending": result["search"]["coefficients"], "coordinate_matrix": maps["mapping"]["matrix"]},
        "historical_policy_exclusion": {"base_M17_parity": M17_PARITY, "old_initial_SHA_sample_contains_base_parity": False, "old_initial_selected_centres_contain_base_parity": False, "later_M22_floor_excludes_M17_supported_base_parity": True},
        "execution_oracle_audit": {"rule_is_retrospectively_calibrated": True, "execution_is_target_blind": True, "allowed_inputs": p["execution_blindness"]["allowed_artifact_inputs"], "geometry_reads": base.read(LOCAL / "geometry-data-access.json"), "worker_reads": base.read(LOCAL / "worker-data-access.json"), "replay_reads": base.read(LOCAL / "replay-data-access.json"), "excluded_from_execution": p["execution_blindness"]["forbidden_inputs"]},
        "completed_search": {"rank_before": INITIAL_RANK, "rank_after_worker_mod2_admission": result["rank_lower_bound"], "rank_after_exact_map_replay": replay["rank_lower_bound"], "mod_2_certified_rank_lower_bound": mod2["rank_lower_bound"], "mod_3_5_ranks": {str(audit["modulus"]): audit["finite_column_rank"] for audit in modl["audits"]}, "height": p["height"], "seconds_per_chart": p["seconds_per_chart"], "search_status": result["search"]["status"], "returned_finite_curve_points": result["search"]["finite_curve_points"]},
        "boundary": "This successful execution was target-blind but the rule was retrospectively calibrated. It proves a lower bound for the completed cloud only; it is not a prospective selector validation or exact-rank proof.",
        "reproducing_command": "python3 elliptic-curves/cas/run_curve302_m26_fibre_orbit58145_extension472.py launch",
    }
    OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print("REPORTED CURVE302 M26 FIBRE|rank>={}".format(mod2["rank_lower_bound"]), flush=True)


def launch() -> None:
    freeze()
    ledger = {"schema": "elliptic-curves.curve302-m26-fibre-orbit58145-extension472-launch.v1", "status": "RUNNING", "stages": []}
    base.checkpoint(LOCAL / "ledger.json", ledger)
    for name, command, seconds in (("geometry", [base.SAGE, str(Path(__file__).resolve()), "geometry"], 180), ("worker", [sys.executable, str(Path(__file__).resolve()), "worker"], 120), ("replay", [sys.executable, str(Path(__file__).resolve()), "replay"], 120), ("certify", [sys.executable, str(Path(__file__).resolve()), "certify"], 900), ("report", [sys.executable, str(Path(__file__).resolve()), "report"], 120)):
        base.launch_stage(name, command, seconds, ledger)
    ledger.update(status="PASS", report=base.rel(OUTPUT), report_sha256=base.sha(OUTPUT))
    base.checkpoint(LOCAL / "ledger.json", ledger)


template.source_hashes = source_hashes
base.source_hashes = source_hashes
template.frozen_centre = frozen_centre
base.centre_from_orbit = frozen_centre
base.protocol = protocol


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("freeze", "geometry", "worker", "replay", "certify", "report", "launch"))
    args = parser.parse_args()
    globals()[args.stage]()

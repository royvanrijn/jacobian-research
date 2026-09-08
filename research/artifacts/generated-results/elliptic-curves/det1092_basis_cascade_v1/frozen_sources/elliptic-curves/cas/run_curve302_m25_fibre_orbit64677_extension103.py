#!/usr/bin/env python3
"""One frozen M25 parity-fibre arm for the next curve-302 strict direction.

The rule is retrospectively calibrated but execution-blind: retain M17 parity
11144 (whose parent degree-two label is rational norm-ten orbit64677), append
the fixed eight-bit M25 extension103, and choose the metric-shortest centre in
that parity.  Geometry, search, and replay never read a residual target or the
retrospective selector diagnostics.
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
import run_curve302_m24_bisection_orbit117420 as base


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
ART = ROOT / "artifacts/generated-results/elliptic-curves"
SEED_CLOUD = ART / "curve302_m24_bisection_orbit117420_mod2_v1.json"
LOCAL = ROOT / "artifacts/local/elliptic-curves/curve302-m25-fibre-orbit64677-extension103-v1"
OUTPUT = ART / "curve302_m25_fibre_orbit64677_extension103_v1.json"
MOD2 = ART / "curve302_m25_fibre_orbit64677_extension103_mod2_v1.json"
MODL = ART / "curve302_m25_fibre_orbit64677_extension103_modl_v1.json"
M17_PARITY = 11144
EXTENSION_MASK = 103
INITIAL_RANK = 25
EXPECTED_CENTRE = [4, -2, -6, 1, 0, -2, 4, -3, -1, -3, 4, -1, 0, 1, 2, 0, 0, -1, -1, 1, 0, 0, 1, -1, 0]
EXPECTED_PARITY = 13511560


# Reuse the separately replayed backend, map validation, finite certificates,
# and source manifest.  Only the seed, frozen selector, and Sage geometry are
# changed for this one centre arm.
BASE_SOURCES = base.source_hashes
BASE_PROTOCOL = base.protocol
base.M24 = SEED_CLOUD
base.LOCAL = LOCAL
base.OUTPUT = OUTPUT
base.MOD2 = MOD2
base.MODL = MODL
base.M24_DIMENSION = INITIAL_RANK
base.ORBIT_MASK = 64677


def source_hashes() -> dict[str, str]:
    result = BASE_SOURCES()
    result[base.rel(Path(__file__))] = cert.hashed(Path(__file__))
    return result


def seed_from_m25() -> dict:
    cloud = base.read(SEED_CLOUD)
    if cloud["status"] != "COMPLETE_DECLARED_FINITE_AUDIT" or cloud["rank_lower_bound"] != INITIAL_RANK:
        raise ArithmeticError("certified rank-25 source cloud changed")
    points = cloud["independent_points"]
    if len(points) != INITIAL_RANK:
        raise ArithmeticError("rank-25 cloud lacks its independent basis")
    model = tuple(map(cert.F, cloud["curve"]))
    proof = cloud["rank_certificate"]
    actual = base.checked_rank(
        model,
        [tuple(map(cert.F, point)) for point in points],
        [signature["prime"] for signature in proof["signatures"]],
        proof["no_rational_2_torsion_prime"],
    )
    if base.digest(actual) != base.digest(proof):
        raise ArithmeticError("rank-25 source independent certificate did not replay")
    return {
        "family": cloud["family"],
        "parameter": cloud["parameter"],
        "curve": cloud["curve"],
        "points": points,
        "generic_points": points[:17],
        "rank_certificate": proof,
        "input_cloud": base.rel(SEED_CLOUD),
        "input_cloud_sha256": base.sha(SEED_CLOUD),
    }


def frozen_centre() -> dict:
    orbit = base.orbit_row()
    if orbit["orbit_mask"] != 64677 or orbit["minimum_generic_MW17_norm"] != 10 or orbit["category"] != "rational":
        raise ArithmeticError("frozen parent bisection orbit changed")
    if len(EXPECTED_CENTRE) != INITIAL_RANK or sum((value & 1) << index for index, value in enumerate(EXPECTED_CENTRE)) != EXPECTED_PARITY:
        raise ArithmeticError("frozen M25 centre parity changed")
    return {
        "base_M17_parity_mask": M17_PARITY,
        "new_M25_extension_mask": EXTENSION_MASK,
        "M25_parity_mask": EXPECTED_PARITY,
        "representative": EXPECTED_CENTRE,
        "selection_rule": "metric-shortest representative in the fixed M25 parity; the stored word is rederived in the Sage geometry stage",
        "generic_degree_two_orbit": orbit,
    }


def protocol() -> dict:
    result = copy.deepcopy(BASE_PROTOCOL())
    result.update(
        schema="elliptic-curves.curve302-m25-fibre-orbit64677-extension103.v1",
        initial_rank=INITIAL_RANK,
        centre_selection={
            "rule": "fix M17 parity11144 (generic rational norm-ten degree-two orbit64677), extend it by M25 bits103, then select the shortest representative in the rounded specialized M25 metric",
            "base_M17_parity_mask": M17_PARITY,
            "frozen_M25_extension_mask": EXTENSION_MASK,
            "generic_degree_two_orbit_mask": 64677,
            "calibration_status": "retrospectively calibrated; execution-blind",
            "minimum_change": "one additional M25 centre/chart; earlier sealed policies remain unmodified",
            "selection_inputs": [base.rel(path) for path in (base.PARENT, base.LATTICE, base.ORBITS, SEED_CLOUD)],
        },
        execution_blindness={
            "allowed_artifact_inputs": [base.rel(path) for path in (SEED_CLOUD, base.PARENT, base.LATTICE, base.ORBITS)],
            "forbidden_inputs": [
                "the residual-visibility diagnostics",
                "all public missing-point coordinates and public rank-31 point lists",
                "post-search finite-rank outcomes",
            ],
            "enforcement": "geometry, worker, and replay install an artifact-read guard before their arithmetic imports and serialize accepted read paths",
        },
        boundary="A single retrospectively calibrated target-blind M25 parity-fibre arm. It can certify only the completed point cloud's lower bound; it is neither a prospective ranking validation nor an exact-rank proof.",
    )
    return result


def geometry() -> None:
    """Derive the frozen parity's shortest M25 centre without target data."""

    base.install_execution_guard()
    p = base.assert_protocol()
    out = LOCAL / "maps.json"
    if out.exists():
        raise FileExistsError("preserve frozen M25 fibre map")
    seed, centre = base.read(LOCAL / "seed.json"), base.read(LOCAL / "centre.json")
    if base.parent_prefix_in_short_model() != seed["generic_points"]:
        raise ArithmeticError("rank-25 prefix is not the certified M17 parent specialization")
    from sage.all import EllipticCurve, QQ, ZZ, matrix, pari, vector

    cvp = SourceFileLoader("curve302_fibre_execution_cvp", str(CAS / "audit_curve302_residual_visibility_geometry.sage")).load_module()
    mapper = SourceFileLoader("curve302_fibre_execution_mapper", str(CAS / "factor_free_pari_mapping.sage")).load_module()
    model = tuple(map(cert.F, seed["curve"]))
    points = tuple(tuple(map(cert.F, point)) for point in seed["points"])
    getcontext().prec = 110
    pari.default("realprecision", 110)
    pari.allocatemem(256_000_000, 1_073_741_824, silent=True)
    height = pari(EllipticCurve(QQ, [QQ(value) for value in seed["curve"]])).ellheightmatrix([list(point) for point in points], precision=384)
    decimal = [[Decimal(str(height[i, j])) for j in range(INITIAL_RANK)] for i in range(INITIAL_RANK)]
    if max(abs(decimal[i][j] - decimal[j][i]) for i in range(INITIAL_RANK) for j in range(INITIAL_RANK)) > Decimal("1e-90"):
        raise ArithmeticError("M25 height matrix asymmetry exceeds guard")
    gram = matrix(ZZ, INITIAL_RANK, INITIAL_RANK, [int((entry * Decimal(1_000_000)).to_integral_value()) for row in decimal for entry in row])
    change = matrix(ZZ, pari(gram).qflllgram()).transpose()
    if abs(change.det()) != 1:
        raise ArithmeticError("M25 LLL change is not unimodular")
    inverse = change.inverse()
    original_parity = vector(ZZ, [(M17_PARITY >> index) & 1 for index in range(17)] + [(EXTENSION_MASK >> bit) & 1 for bit in range(8)])
    reduced_parity = vector(ZZ, (matrix(ZZ, 1, INITIAL_RANK, original_parity) * inverse).row(0))
    reduced_parity = vector(ZZ, [int(value) % 2 for value in reduced_parity])
    reduced_gram = change * gram * change.transpose()
    dd = cvp.RoundedMetricCVP(reduced_gram, "dd", 192)
    mpfr = cvp.RoundedMetricCVP(reduced_gram, "mpfr", 192)
    u, distance = cvp.cross_precision_closest(dd, mpfr, -vector(reduced_parity) / 2)
    derived = vector(ZZ, (matrix(ZZ, 1, INITIAL_RANK, reduced_parity + 2 * u) * change).row(0))
    if list(map(int, derived)) != centre["representative"] or sum((int(value) & 1) << index for index, value in enumerate(derived)) != EXPECTED_PARITY:
        raise ArithmeticError("frozen parity no longer derives the pinned M25 centre")
    mapping = mapper.mapping(model, points, centre)
    checkpoint = base.checkpoint
    checkpoint(out, {
        "schema": "elliptic-curves.curve302-m25-fibre-map.v1",
        "status": "COMPLETE_DECLARED_SINGLE_MAP",
        "protocol_sha256": base.sha(LOCAL / "protocol.json"),
        "seed_sha256": base.sha(LOCAL / "seed.json"),
        "centre_sha256": base.sha(LOCAL / "centre.json"),
        "metric_gram": [list(map(int, row)) for row in gram.rows()],
        "change_of_basis": [list(map(int, row)) for row in change.rows()],
        "centre_CVP_distance": str(distance),
        "centre_derived_from_frozen_parity": True,
        "mapping": mapping,
    })
    checkpoint(LOCAL / "geometry-data-access.json", sorted(base.READS))
    print("FROZEN CURVE302 M25 FIBRE MAP|orbit=64677|extension=103", flush=True)


def freeze() -> None:
    base.freeze()


def worker() -> None:
    base.worker()


def replay() -> None:
    base.replay()


def certify() -> None:
    base.certify()
    mod2, modl = base.read(MOD2), base.read(MODL)
    if mod2["rank_lower_bound"] < 26 or any(audit["finite_column_rank"] < 26 for audit in modl["audits"]):
        raise ArithmeticError("the single frozen M25 parity-fibre arm did not certify rank at least 26")


def report() -> None:
    p = base.assert_protocol()
    if OUTPUT.exists():
        raise FileExistsError("preserve immutable M25 fibre report")
    names = ("seed.json", "centre.json", "maps.json", "result.json", "replay.json", "certification-ledger.json")
    local = [LOCAL / name for name in names]
    if any(not path.exists() for path in local) or not MOD2.exists() or not MODL.exists():
        raise ArithmeticError("complete execution and certificates are required")
    seed, centre, maps, result, replay, certificate = (base.read(path) for path in local)
    mod2, modl, old = base.read(MOD2), base.read(MODL), base.read(base.OLD_M17_MAPS)
    if any(row["parity"] == M17_PARITY for row in old["sample"]) or any(row["parity"] == M17_PARITY for row in old["centres"]):
        raise ArithmeticError("historical deep arm unexpectedly included the M17 base parity")
    if mod2["rank_lower_bound"] < 26 or any(audit["finite_column_rank"] < 26 for audit in modl["audits"]):
        raise ArithmeticError("rank-26 certificate missing")
    report_data = {
        "schema": "elliptic-curves.curve302-m25-fibre-orbit64677-extension103-report.v1",
        "status": "PASS_CERTIFIED_RANK_AT_LEAST_26",
        "inputs": {
            **base.frozen_inputs(),
            base.rel(base.OLD_M17_MAPS): base.sha(base.OLD_M17_MAPS),
            **{base.rel(path): base.sha(path) for path in (*local, MOD2, MODL)},
        },
        "sources": source_hashes(),
        "frozen_rule": p["centre_selection"],
        "exact_M25_centre": {
            **centre,
            "centre_point_short_model": result["search"]["base_point"],
            "reduced_quartic_coefficients_ascending": result["search"]["coefficients"],
            "coordinate_matrix": maps["mapping"]["matrix"],
        },
        "historical_policy_exclusion": {
            "base_M17_parity": M17_PARITY,
            "old_initial_SHA_sample_contains_base_parity": False,
            "old_initial_selected_centres_contain_base_parity": False,
            "later_M22_floor_excludes_M17_supported_base_parity": True,
        },
        "execution_oracle_audit": {
            "rule_is_retrospectively_calibrated": True,
            "execution_is_target_blind": True,
            "allowed_inputs": p["execution_blindness"]["allowed_artifact_inputs"],
            "geometry_reads": base.read(LOCAL / "geometry-data-access.json"),
            "worker_reads": base.read(LOCAL / "worker-data-access.json"),
            "replay_reads": base.read(LOCAL / "replay-data-access.json"),
            "excluded_from_execution": p["execution_blindness"]["forbidden_inputs"],
        },
        "completed_search": {
            "rank_before": INITIAL_RANK,
            "rank_after_worker_mod2_admission": result["rank_lower_bound"],
            "rank_after_exact_map_replay": replay["rank_lower_bound"],
            "mod_2_certified_rank_lower_bound": mod2["rank_lower_bound"],
            "mod_3_5_ranks": {str(audit["modulus"]): audit["finite_column_rank"] for audit in modl["audits"]},
            "height": p["height"],
            "seconds_per_chart": p["seconds_per_chart"],
            "search_status": result["search"]["status"],
            "returned_finite_curve_points": result["search"]["finite_curve_points"],
        },
        "boundary": "This single rule was chosen after retrospective analysis, so the successful target-blind execution is not prospective selection validation. It proves only a rank lower bound for its completed cloud and does not alter the earlier sealed campaigns or prove exact rank.",
        "reproducing_command": "python3 elliptic-curves/cas/run_curve302_m25_fibre_orbit64677_extension103.py launch",
    }
    OUTPUT.write_text(json.dumps(report_data, indent=2, sort_keys=True) + "\n")
    print("REPORTED CURVE302 M25 FIBRE|rank>={}".format(mod2["rank_lower_bound"]), flush=True)


def launch() -> None:
    freeze()
    commands = (("geometry", [base.SAGE, str(Path(__file__).resolve()), "geometry"], 180), ("worker", [sys.executable, str(Path(__file__).resolve()), "worker"], 120), ("replay", [sys.executable, str(Path(__file__).resolve()), "replay"], 120), ("certify", [sys.executable, str(Path(__file__).resolve()), "certify"], 900), ("report", [sys.executable, str(Path(__file__).resolve()), "report"], 120))
    ledger = {"schema": "elliptic-curves.curve302-m25-fibre-orbit64677-extension103-launch.v1", "status": "RUNNING", "stages": []}
    base.checkpoint(LOCAL / "ledger.json", ledger)
    for name, command, seconds in commands:
        base.launch_stage(name, command, seconds, ledger)
    ledger.update(status="PASS", report=base.rel(OUTPUT), report_sha256=base.sha(OUTPUT))
    base.checkpoint(LOCAL / "ledger.json", ledger)


base.source_hashes = source_hashes
base.protocol = protocol
base.seed_from_m24 = seed_from_m25
base.centre_from_orbit = frozen_centre


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("freeze", "geometry", "worker", "replay", "certify", "report", "launch"))
    args = parser.parse_args()
    globals()[args.stage]()

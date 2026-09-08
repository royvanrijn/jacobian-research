#!/usr/bin/env python3
"""Freeze the exact, adaptive curve-302 strict-residual recovery chain.

This is deliberately an audit manifest, not a replacement protocol.  The
historical rank-17 -> 19 -> 22 -> 24 calibration, and each subsequently
sealed one-chart arm, remain separately addressable.  The manifest records
their hashes and the exact finite certificates that establish every recorded
rank lower bound.

It therefore makes no prospective-selection claim: the post-M24 arm choices
were retrospectively calibrated, although their execution stages were guarded
against the diagnostic and target point files.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
ART = ROOT / "artifacts/generated-results/elliptic-curves"
LOCAL = ROOT / "artifacts/local/elliptic-curves"
OUTPUT = ART / "curve302_residual_strict_bootstrap_chain_v1.json"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rel(path: Path) -> str:
    return str(path.relative_to(ROOT))


def read(path: Path) -> dict:
    return json.loads(path.read_text())


def hashed_paths(paths: list[Path]) -> dict[str, str]:
    if len(paths) != len(set(paths)) or any(not path.is_file() for path in paths):
        raise ArithmeticError("bootstrap-chain input path is missing or duplicated")
    return {rel(path): sha(path) for path in paths}


def check_certificate(mod2: dict, modl: dict, rank: int) -> dict:
    if mod2.get("status") != "COMPLETE_DECLARED_FINITE_AUDIT":
        raise ArithmeticError("mod-2 finite audit is not complete")
    if mod2.get("rank_lower_bound") != rank:
        raise ArithmeticError("mod-2 certificate rank differs from chain rank")
    proof = mod2.get("rank_certificate", {})
    if proof.get("rank_lower_bound") != rank or len(mod2.get("independent_points", [])) != rank:
        raise ArithmeticError("mod-2 independent basis is incomplete")
    if modl.get("status") != "COMPLETE_BOUNDED_QUOTIENT_AUDIT":
        raise ArithmeticError("odd-prime finite audit is not complete")
    ranks = {str(row["modulus"]): row["finite_column_rank"] for row in modl.get("audits", [])}
    if ranks.get("3") != rank or ranks.get("5") != rank:
        raise ArithmeticError("independent mod-3/mod-5 replay does not match chain rank")
    return {
        "mod2_rank_lower_bound": rank,
        "mod3_rank_lower_bound": ranks["3"],
        "mod5_rank_lower_bound": ranks["5"],
        "no_rational_2_torsion_prime": proof["no_rational_2_torsion_prime"],
        "certificate_argument": proof["argument"],
    }


def calibration_states() -> list[dict]:
    parent = ART / "curve302_recovered_mw17_parent_proof_v1.json"
    stages = [
        ("M17", 17, parent),
        ("M19", 19, ART / "curve302_focused_curve302_generic17_mod2_v1.json"),
        ("M22", 22, ART / "curve302_recovered_followup_wave_01_mod2_v1.json"),
        ("M24", 24, ART / "curve302_recovered_followup_wave_02_mod2_v1.json"),
        ("M24-policy-stall", 24, ART / "curve302_recovered_followup_wave_03_mod2_v1.json"),
    ]
    rows = []
    for label, rank, path in stages:
        data = read(path)
        if label == "M17":
            if data.get("status") != "PASS_FULL_ARITHMETIC_MW17_PARENT" or data.get("generic_arithmetic_MW_rank") != rank:
                raise ArithmeticError("certified M17 parent changed")
            basis_size = data["full_basis_size"]
        else:
            if data.get("status") != "COMPLETE_DECLARED_FINITE_AUDIT" or data.get("rank_lower_bound") != rank:
                raise ArithmeticError(f"{label} calibration state changed")
            basis_size = len(data["independent_points"])
            if basis_size != rank:
                raise ArithmeticError(f"{label} calibration basis changed")
        rows.append({
            "stage": label,
            "rank_lower_bound": rank,
            "basis_size": basis_size,
            "certificate_or_cloud": rel(path),
            "sha256": sha(path),
        })
    return rows


ARM_SPECS = (
    # rank_before, rank_after, direction, orbit, shell, extension/padding, local directory, generated artifacts, optional report
    (24, 25, "residual-strict-06", 117420, "rational norm-10", "append seven zero M24 coordinates", "curve302-m24-bisection-orbit117420-v1", "curve302_m24_bisection_orbit117420_mod2_v1.json", "curve302_m24_bisection_orbit117420_modl_v1.json", "curve302_m24_bisection_orbit117420_v1.json"),
    (25, 26, "residual-strict-04", 64677, "rational norm-10", "M25 extension mask 103", "curve302-m25-fibre-orbit64677-extension103-v1", "curve302_m25_fibre_orbit64677_extension103_mod2_v1.json", "curve302_m25_fibre_orbit64677_extension103_modl_v1.json", "curve302_m25_fibre_orbit64677_extension103_v1.json"),
    (26, 27, "residual-strict-02", 58145, "rational norm-10", "M26 extension mask 472", "curve302-m26-fibre-orbit58145-extension472-v1", "curve302_m26_fibre_orbit58145_extension472_mod2_v1.json", "curve302_m26_fibre_orbit58145_extension472_modl_v1.json", "curve302_m26_fibre_orbit58145_extension472_v1.json"),
    (27, 28, "residual-strict-01", 106210, "rational norm-10", "M27 extension mask 673", "curve302-m27-fibre-orbit106210-extension673-v1", "curve302_m27_fibre_orbit106210_extension673_mod2_v1.json", "curve302_m27_fibre_orbit106210_extension673_modl_v1.json", "curve302_m27_fibre_orbit106210_extension673_v1.json"),
    (28, 29, "residual-strict-05", 4761, "rational norm-10", "M28 extension mask 1772", "curve302-m28-fibre-orbit4761-extension1772-v1", "curve302_m28_fibre_orbit4761_extension1772_mod2_v1.json", "curve302_m28_fibre_orbit4761_extension1772_modl_v1.json", None),
    (29, 30, "residual-strict-03", 17845, "genus-one nearby norm-8", "M29 extension mask 254", "curve302-m29-fibre-orbit17845-extension254-v1", "curve302_m29_fibre_orbit17845_extension254_mod2_v1.json", "curve302_m29_fibre_orbit17845_extension254_modl_v1.json", None),
    (30, 31, "residual-strict-07", 114326, "genus-one nearby norm-8", "M29 extension mask 1502; M30 padding coordinate 0", "curve302-m30-padded-fibre-orbit114326-extension1502-v1", "curve302_m30_padded_fibre_orbit114326_extension1502_mod2_v1.json", "curve302_m30_padded_fibre_orbit114326_extension1502_modl_v1.json", None),
)


CHECKPOINT_FILES = (
    "protocol.json", "seed.json", "centre.json", "maps.json", "result.json", "replay.json",
    "geometry-data-access.json", "worker-data-access.json", "replay-data-access.json",
    "certification-ledger.json", "independent-rank.log", "independent-rank.supervisor.json",
    "mod2-build.log", "mod2-build.supervisor.json", "mod2-check.log", "mod2-check.supervisor.json",
    "modl-build.log", "modl-build.supervisor.json", "modl-check.log", "modl-check.supervisor.json",
)


def arm_record(spec: tuple) -> dict:
    before, after, direction, orbit, shell, extension, local_name, mod2_name, modl_name, report_name = spec
    local = LOCAL / local_name
    paths = [local / name for name in CHECKPOINT_FILES]
    protocol, centre, maps, result, replay = (read(local / name) for name in ("protocol.json", "centre.json", "maps.json", "result.json", "replay.json"))
    if protocol.get("initial_rank") != before or protocol.get("chart_count") != 1:
        raise ArithmeticError("sealed arm rank or chart count changed")
    if not protocol.get("execution_blindness", {}).get("forbidden_inputs"):
        raise ArithmeticError("sealed arm is missing its execution-blindness declaration")
    if result.get("rank_lower_bound", before) < before or replay.get("rank_lower_bound", before) < before:
        raise ArithmeticError("arm worker/replay state regressed")
    if maps.get("protocol_sha256") != sha(local / "protocol.json") or result.get("protocol_sha256") != sha(local / "protocol.json"):
        raise ArithmeticError("map or worker no longer binds to its frozen protocol")
    mod2_path, modl_path = ART / mod2_name, ART / modl_name
    mod2, modl = read(mod2_path), read(modl_path)
    certificate = check_certificate(mod2, modl, after)
    if len(mod2["independent_points"]) != after:
        raise ArithmeticError("chain witness cloud has wrong size")
    report = None if report_name is None else ART / report_name
    all_paths = paths + [mod2_path, modl_path] + ([] if report is None else [report])
    for path, expected in protocol.get("inputs", {}).items():
        candidate = ROOT / path
        if not candidate.is_file() or sha(candidate) != expected:
            raise ArithmeticError("frozen arm input hash changed: " + path)
    return {
        "direction": direction,
        "rank_before": before,
        "rank_after": after,
        "degree_two_orbit_mask": orbit,
        "generic_shell": shell,
        "frozen_extension_rule": extension,
        "exact_new_independent_witness": mod2["independent_points"][-1],
        "certificate": certificate,
        "frozen_protocol": protocol,
        "exact_centre": centre,
        "map_status": maps["status"],
        "worker_status": result["status"],
        "replay_status": replay["status"],
        "execution_reads": {
            "geometry": read(local / "geometry-data-access.json"),
            "worker": read(local / "worker-data-access.json"),
            "replay": read(local / "replay-data-access.json"),
        },
        "immutable_checkpoints": hashed_paths(all_paths),
    }


def build() -> dict:
    calibration = calibration_states()
    arms = [arm_record(spec) for spec in ARM_SPECS]
    if [row["rank_before"] for row in arms] != [24, 25, 26, 27, 28, 29, 30] or [row["rank_after"] for row in arms] != [25, 26, 27, 28, 29, 30, 31]:
        raise ArithmeticError("rank transition sequence changed")
    if any(row["rank_after"] != row["rank_before"] + 1 for row in arms):
        raise ArithmeticError("this is no longer a one-direction-at-a-time chain")
    return {
        "schema": "elliptic-curves.curve302-residual-strict-bootstrap-chain.v1",
        "status": "PASS_EXACT_CHAIN_PRESERVED_RANK_AT_LEAST_31",
        "audit_source_sha256": sha(Path(__file__)),
        "calibration_states": calibration,
        "recovery_arms": arms,
        "chain_interpretation": {
            "selection_status": "Every post-M24 arm is retrospectively calibrated. Its geometry, worker, and replay were instead target-blind under the separately frozen execution guard recorded above.",
            "rank_claim": "Each step is an independently certified rank lower bound, by its retained mod-2 certificate and independent mod-3/mod-5 finite checks.",
            "not_claimed": [
                "a prospective target-free selector has been validated",
                "an unconditional exact-rank upper bound for curve 302",
                "a theorem about every parity class or every degree-two multisection orbit",
            ],
        },
        "reproducing_command": "python3 elliptic-curves/cas/audit_curve302_residual_strict_bootstrap_chain.py --check",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build", action="store_true")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.build == args.check:
        parser.error("choose exactly one of --build or --check")
    payload = build()
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.build:
        if OUTPUT.exists():
            raise FileExistsError("preserve existing bootstrap-chain manifest")
        OUTPUT.write_text(rendered)
    elif OUTPUT.read_text() != rendered:
        raise ArithmeticError("bootstrap-chain manifest did not replay")
    print("CURVE302BOOTSTRAP|states=5|arms=7|rank_at_least=31|status=PASS", flush=True)


if __name__ == "__main__":
    main()

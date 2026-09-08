#!/usr/bin/env python3
"""Matched, generic-only norm-ten bisection-orbit exposure on det-1092 fibres.

The discovery path receives a redacted MW17 parent, the complete generic
degree-two orbit table, and a frozen score-only candidate set.  It does not
read the earlier point-exposure ledger, a catalogue, curve-302 exceptional
points, or any public rank outcome.  The old policy is opened only by the
separate post-terminal comparison stage.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

import certify_compact_r17_candidates as cert
import pari_pointed_backend as backend
from memory_rank_certificate import checked_rank
from pointed_quartic_search import PointedQuarticSearch
from research_runtime.memory_store import MemoryFactStore
from research_runtime.quotient_only_reduction import QuotientOnlyReductionCache as ReductionCache
from research_runtime.search_state import raw_state
from research_runtime.store import checkpoint, digest
from research_runtime.supervisor import Limits, run


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
ART = ROOT / "artifacts/generated-results/elliptic-curves"
CAMPAIGN_TAG = os.environ.get("DET1092_BISECTION_CAMPAIGN", "v1")
if not re.fullmatch(r"v[1-9][0-9]*", CAMPAIGN_TAG):
    raise ValueError("DET1092_BISECTION_CAMPAIGN must be a version tag such as v2")
LOCAL = ROOT / "artifacts/local/elliptic-curves" / ("det1092-bisection-orbit-holdout-" + CAMPAIGN_TAG)
# The candidate roster is parameterized in this already-certified generic
# determinant-reduced chart.  It contains the same ordered seventeen generic
# sections as the normalized parent, transported by an exact generic
# isomorphism; it contains no exceptional-fibre point data.
PARENT = ART / "curve302_recovered_mw17_parent_v1.json"
REDUCED_PARENT = ART / "det1092_reduced_parameter_chart_v1/reduced-parent.json"
LATTICE = ART / "curve302_parent_degree2_multisection_lattice_v1.json"
ORBITS = ART / "curve302_parent_degree2_multisection_orbits_v1.tsv"
SELECTION = ROOT / "artifacts/local/elliptic-curves/det1092-record-scale-selection-v1/selection-result.json"
OLD_LEDGER = ROOT / "artifacts/local/elliptic-curves/det1092-record-scale-points-v1/ledger.json"
SAGE = "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python"
OUTPUT_COMPARISON = ART / ("det1092_bisection_orbit_holdout_comparison_" + CAMPAIGN_TAG + ".json")
CHARTS = 49
INITIAL_RANK = 17
GENERIC_NORM = 10


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rel(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def sources() -> dict[str, str]:
    names = [
        "det1092_bisection_orbit_holdout.py",
        "prepare_det1092_bisection_orbit_seed.sage",
        "prepare_det1092_bisection_orbit_maps.sage",
        "audit_recorded_point_mod2_rank_v3.py",
        "audit_retained_cloud_modl.py",
        "verify_factor_free_rank.sage",
    ]
    return {
        **backend.sources(),
        **{rel(CAS / name): cert.hashed(CAS / name) for name in names},
    }


def read(path: Path):
    return cert.read(path)


def parent_redaction() -> dict:
    parent, reduced = read(PARENT), read(REDUCED_PARENT)
    required = ("generic_height_gram",)
    if int(parent["height_determinant"]) != 1092 or any(key not in parent for key in required):
        raise ArithmeticError("full determinant-1092 MW17 parent input changed")
    reduced_required = ("a_invariants", "basis_weierstrass_coordinates", "source_sha256")
    if any(key not in reduced for key in reduced_required):
        raise ArithmeticError("generic determinant-reduced chart input changed")
    if any(item["numerator"] for item in reduced["a_invariants"][:3]):
        raise ArithmeticError("generic determinant-reduced chart is not short")
    sections = reduced["basis_weierstrass_coordinates"]
    if len(sections) != INITIAL_RANK or any(len(row) != 2 for row in sections):
        raise ArithmeticError("reduced chart no longer exports exactly seventeen sections")
    return {
        "a_invariants": reduced["a_invariants"],
        "basis_weierstrass_coordinates": sections,
        "generic_height_gram": parent["generic_height_gram"],
        "normalized_parent_sha256": sha(PARENT),
        "reduced_chart_parent_sha256": sha(REDUCED_PARENT),
    }


def candidate_set() -> list[dict]:
    selection = read(SELECTION)
    if selection["status"] != "PASS" or len(selection["selected"]) != 48:
        raise ArithmeticError("the score-only held-out selection is incomplete")
    rows = []
    for row in selection["selected"]:
        item = {key: row[key] for key in ("id", "family", "parameter", "model", "height_bin", "stratum", "score_units")}
        if item["family"] != "det1092-reduced" or len(item["model"]) != 5:
            raise ArithmeticError("unexpected held-out determinant-1092 candidate")
        rows.append(item)
    if len({row["id"] for row in rows}) != 48:
        raise ArithmeticError("held-out identifiers are not unique")
    return rows


def freeze() -> None:
    if LOCAL.exists() or OUTPUT_COMPARISON.exists():
        raise FileExistsError("preserve frozen bisection-orbit holdout evidence")
    LOCAL.mkdir(parents=True)
    redacted = parent_redaction()
    rows = candidate_set()
    parent_file = LOCAL / "parent-sections.json"
    candidates_file = LOCAL / "candidates.json"
    checkpoint(parent_file, redacted)
    checkpoint(candidates_file, {"schema": "elliptic-curves.det1092-bisection-heldout-candidates.v1", "rows": rows,
                                 "selection_sha256": sha(SELECTION)})
    lattice = read(LATTICE)
    if lattice["status"] != "PASS_COMPLETE_DEGREE2_TRANSLATION_QUOTIENT" or sha(ORBITS) != lattice["orbits_tsv_sha256"]:
        raise ArithmeticError("complete generic degree-two orbit data changed")
    protocol = {
        "schema": "elliptic-curves.det1092-bisection-orbit-holdout.v1",
        "campaign_tag": CAMPAIGN_TAG,
        "sources": sources(),
        "inputs": {rel(path): sha(path) for path in (parent_file, candidates_file, LATTICE, ORBITS)},
        "rows": rows,
        "initial_rank": INITIAL_RANK,
        "generic_norm": GENERIC_NORM,
        "complete_generic_rational_norm10_orbits": lattice["rational_bisections"]["translation_orbits"],
        "charts_per_fibre": CHARTS,
        "height": 125000,
        "seconds_per_chart": 10,
        "rss_bytes": 2147483648,
        "geometry_seconds": 600,
        "worker_seconds": 1200,
        "replay_seconds": 1200,
        "certificate_seconds": 600,
        "maximum_workers": 1,
        "gp_sha256": cert.hashed(Path("/usr/bin/gp")),
        "centre_policy": {
            "enumeration": "all rational degree-two M/2M orbits with exact generic MW17 minimum norm 10",
            "ranking": "descending exact specialized rounded-height norm; then ascending generic l1, generic linfinity, orbit mask",
            "metric": "384-bit canonical-height Gram rounded at 10^6; every selected generic orbit word is rechecked against the exact determinant-1092 Gram",
            "selection": "top 49 from the full norm-10 orbit list on each fibre",
            "no_sampling": "No SHA parity sample or 2048-candidate bottleneck.",
            "no_novelty_floor": "Every chart uses the original seventeen generic sections only; there are no adaptive waves or subgroup-coordinate floor.",
        },
        "execution_blindness": {
            "discovery_inputs": [rel(path) for path in (parent_file, candidates_file, LATTICE, ORBITS)],
            "forbidden": [
                "curve302 exceptional points or recovered M24/M25/M31 clouds",
                "the completed old point-exposure ledger or rank outcomes",
                "catalogue or novelty queries",
                "outcome-dependent candidate replacement or chart selection",
            ],
            "enforcement": "geometry, worker and replay install an artifact-read guard; the old-policy comparison is a separate terminal-only stage",
        },
        "equal_budget": "49 factor-free charts per fixed held-out fibre, height 125000 and ten seconds per chart; all maps freeze before any point worker starts",
        "certification": "Every returned cloud receives independent mod-2, mod-3, mod-5 and standalone finite-rank certificates. A gain proves only a rank lower bound.",
        "boundary": "This is a finite prospective bisection-orbit exposure experiment, not an exact-rank proof, a general bisection-solubility theorem, or a claim about unsearched norm-ten orbits beyond the declared top-49 rule.",
    }
    checkpoint(LOCAL / "protocol.json", protocol)
    print("FROZEN DET1092 BISECTION HOLDOUT|fibres=48|charts_per_fibre=49|norm10_orbits={}".format(
        protocol["complete_generic_rational_norm10_orbits"]), flush=True)


def campaign() -> dict:
    protocol = read(LOCAL / "protocol.json")
    if protocol["sources"] != sources() or any(sha(ROOT / path) != value for path, value in protocol["inputs"].items()):
        raise ArithmeticError("frozen bisection holdout source or input changed")
    if len(protocol["rows"]) != 48 or protocol["charts_per_fibre"] != CHARTS:
        raise ArithmeticError("held-out roster or equal budget changed")
    return protocol


def install_execution_guard() -> None:
    """Reject all artifact reads outside the generic frozen discovery inputs."""
    approved = {LOCAL.resolve(), LATTICE.resolve(), ORBITS.resolve()}

    def guard(event, args):
        if event != "open" or not args or not isinstance(args[0], (str, bytes)):
            return
        path = Path(args[0]).resolve()
        artifact_root = ROOT / "artifacts"
        if not path.is_relative_to(artifact_root):
            return
        if path.is_relative_to(LOCAL) or path in approved:
            return
        raise PermissionError("bisection discovery rejected artifact input: " + rel(path))

    sys.addaudithook(guard)


def configure(index: int) -> tuple[dict, Path]:
    protocol = campaign()
    if not 0 <= index < len(protocol["rows"]):
        raise ValueError("frozen held-out index required")
    row = protocol["rows"][index]
    return row, LOCAL / row["id"]


def initial_state(seed: dict):
    model = tuple(map(cert.F, seed["curve"]))
    points = tuple(tuple(map(cert.F, point)) for point in seed["points"])
    cache = ReductionCache(MemoryFactStore())
    state = raw_state(model, points, cache=cache, prime_bound=1000)
    if (state.rank != INITIAL_RANK
            or tuple(tuple(cert.F(value) for value in point) for point in state.basis) != points):
        raise ArithmeticError("specialized generic MW17 seed changed")
    proof = checked_rank(model, points, state.reductions.primes, state.no_two_torsion_prime)
    if digest(proof) != digest(seed["rank_certificate"]):
        raise ArithmeticError("seed finite-independence certificate changed")
    return model, points, state


def worker(index: int) -> None:
    install_execution_guard()
    row, directory = configure(index)
    protocol, seed, maps = campaign(), read(directory / "seed.json"), read(directory / "maps.json")
    output = directory / "result.json"
    if output.exists():
        raise FileExistsError("preserve completed fixed-budget bisection worker")
    model, points, state = initial_state(seed)
    if maps["status"] != "COMPLETE_DECLARED_BISECTION_MAPS" or maps["protocol_sha256"] != sha(LOCAL / "protocol.json"):
        raise ArithmeticError("all frozen bisection maps are required")
    if len(maps["rows"]) != CHARTS or maps["seed_sha256"] != sha(directory / "seed.json"):
        raise ArithmeticError("fixed equal chart roster changed")
    data = {
        "schema": "elliptic-curves.det1092-bisection-orbit-worker.v1", "status": "RUNNING_FIXED_EQUAL_BUDGET",
        "protocol_hash": digest(protocol), "seed_sha256": sha(directory / "seed.json"), "maps_sha256": sha(directory / "maps.json"),
        "family": row["family"], "parameter": row["parameter"], "curve": seed["curve"], "generic_points": seed["points"],
        "initial_state": state.record(), "final_state": state.record(), "charts": [], "rank_lower_bound": INITIAL_RANK,
        "policy": "all 49 maps use the fixed original MW17 seed; no admission, adaptive wave, or novelty floor operates during search",
    }
    checkpoint(output, data)
    incomplete = False
    for chart in maps["rows"]:
        centre = chart["centre"]
        search = PointedQuarticSearch(state=state, centre={"coefficients": centre["representative"]},
                                      coordinate_policy=chart["coordinate_policy"])
        transcript, found = backend.execute(search, chart, protocol["height"], protocol["seconds_per_chart"], protocol["gp_sha256"])
        incomplete |= transcript["status"] != "bounded_search_complete"
        data["charts"].append({"index": chart["index"], "centre": centre, "mapping": chart,
                               "search": transcript, "returned_point_count": len(found)})
        checkpoint(output, data)
        print("BISECTION ORBIT", row["id"], chart["index"] + 1, "/", CHARTS, transcript["status"], flush=True)
    data["status"] = "CENSORED_FIXED_EQUAL_BUDGET" if incomplete else "COMPLETE_FIXED_EQUAL_BUDGET"
    checkpoint(output, data)
    if incomplete:
        raise RuntimeError("preserved censored bisection chart; certificates are not launched")


def replay(index: int) -> None:
    install_execution_guard()
    row, directory = configure(index)
    protocol, seed, maps, data = campaign(), read(directory / "seed.json"), read(directory / "maps.json"), read(directory / "result.json")
    model, points, state = initial_state(seed)
    if data["protocol_hash"] != digest(protocol) or data["maps_sha256"] != sha(directory / "maps.json"):
        raise ArithmeticError("fixed worker binding changed")
    if data["initial_state"] != state.record() or data["final_state"] != state.record() or data["curve"] != seed["curve"]:
        raise ArithmeticError("fixed seed state changed")
    if len(data["charts"]) != CHARTS or data["status"] != "COMPLETE_FIXED_EQUAL_BUDGET":
        raise ArithmeticError("only a complete fixed equal-budget run can replay")
    for expected, actual in zip(maps["rows"], data["charts"]):
        if actual["index"] != expected["index"] or actual["centre"] != expected["centre"] or actual["mapping"] != expected:
            raise ArithmeticError("chart/map roster differs")
        search = PointedQuarticSearch(state=state, centre={"coefficients": expected["centre"]["representative"]},
                                      coordinate_policy=expected["coordinate_policy"])
        found = backend.replay(search, expected, actual["search"])
        if len(found) != actual["returned_point_count"]:
            raise ArithmeticError("replayed point cloud cardinality differs")
    print("REPLAYED FIXED BISECTION ORBIT EXPOSURE", row["id"], CHARTS, flush=True)


def certify(index: int) -> dict:
    row, directory = configure(index)
    protocol = campaign()
    result = directory / "result.json"
    ledger = directory / "certification-ledger.json"
    if ledger.exists():
        raise FileExistsError("preserve independent bisection certification")
    prefix = "det1092_bisection_orbit_" + row["id"].replace("-", "_")
    mod2, modl = ART / (prefix + "_mod2_v1.json"), ART / (prefix + "_modl_v1.json")
    if mod2.exists() or modl.exists():
        raise FileExistsError("preserve generated point-cloud certificates")
    certificate_protocol = {
        "schema": "elliptic-curves.det1092-bisection-orbit-certification.v1",
        "inputs": {rel(path): sha(path) for path in (result, directory / "maps.json", directory / "seed.json", LOCAL / "protocol.json")},
        "sources": {rel(CAS / name): cert.hashed(CAS / name) for name in (
            "audit_recorded_point_mod2_rank_v3.py", "audit_retained_cloud_modl.py", "verify_factor_free_rank.sage")},
        "independence": "independent whole-cloud finite quotient certificates modulo 2, 3 and 5",
    }
    checkpoint(directory / "certification-protocol.json", certificate_protocol)
    jobs = [
        ("mod2-build", [sys.executable, str(CAS / "audit_recorded_point_mod2_rank_v3.py"), "--input", str(result),
                        "--input-sha256", sha(result), "--output", str(mod2), "--prime-bound", "997"], ROOT),
        ("mod2-check", [sys.executable, str(CAS / "audit_recorded_point_mod2_rank_v3.py"), "--check", str(mod2)], ROOT),
        ("modl-build", [sys.executable, str(CAS / "audit_retained_cloud_modl.py"), "--input", str(mod2), "--output", str(modl)], ROOT),
        ("modl-check", [sys.executable, str(CAS / "audit_retained_cloud_modl.py"), "--check", str(modl)], ROOT),
    ]
    record = {"status": "RUNNING", "stages": []}
    checkpoint(ledger, record)
    for name, command, cwd in jobs:
        supervision = run(command, limits=Limits(protocol["certificate_seconds"], protocol["rss_bytes"]),
                          log_path=directory / (name + ".log"), checkpoint_path=directory / (name + ".supervisor.json"), cwd=cwd)
        ok = supervision["outcome"] == "completed" and supervision["returncode"] == 0
        record["stages"].append({"name": name, "status": "PASS" if ok else "FAILED_OR_CENSORED", "supervision": supervision})
        checkpoint(ledger, record)
        if not ok:
            record["status"] = "FAILED_OR_CENSORED"
            checkpoint(ledger, record)
            raise ArithmeticError("preserve failed independent finite certificate")
    fresh = directory / "independent-rank"
    fresh.mkdir(exist_ok=False)
    for path in (CAS / "verify_factor_free_rank.sage", mod2):
        shutil.copy2(path, fresh / path.name)
    checkpoint(fresh / "protocol.json", {"files": {path.name: sha(path) for path in (CAS / "verify_factor_free_rank.sage", mod2)},
                                            "seconds": protocol["certificate_seconds"], "rss_bytes": protocol["rss_bytes"]})
    supervision = run([SAGE, str(fresh / "verify_factor_free_rank.sage"), "--input", str(fresh / mod2.name)],
                      limits=Limits(protocol["certificate_seconds"], protocol["rss_bytes"]), log_path=directory / "independent-rank.log",
                      checkpoint_path=directory / "independent-rank.supervisor.json", cwd=fresh)
    ok = supervision["outcome"] == "completed" and supervision["returncode"] == 0
    record["stages"].append({"name": "independent-rank", "status": "PASS" if ok else "FAILED_OR_CENSORED", "supervision": supervision})
    if not ok:
        record["status"] = "FAILED_OR_CENSORED"
        checkpoint(ledger, record)
        raise ArithmeticError("independent rank replay failed")
    cloud, odd = read(mod2), read(modl)
    if cloud["points"][:INITIAL_RANK] != read(directory / "seed.json")["points"]:
        raise ArithmeticError("whole-cloud certificate lost generic prefix")
    odd_ranks = {str(item["modulus"]): item["finite_column_rank"] for item in odd["audits"]}
    if any(rank < cloud["rank_lower_bound"] for rank in odd_ranks.values()):
        raise ArithmeticError("odd-modulus certificate fails to support mod-2 cloud rank")
    record.update(status="PASS", rank_lower_bound=cloud["rank_lower_bound"], discovered_rank_gain=cloud["rank_lower_bound"] - INITIAL_RANK,
                  completed_boxes=CHARTS, point_count=len(cloud["points"]), odd_modulus_ranks=odd_ranks,
                  certificates={rel(path): sha(path) for path in (mod2, modl)})
    checkpoint(ledger, record)
    return record


def launch() -> None:
    protocol = campaign()
    ledger_path = LOCAL / "ledger.json"
    if ledger_path.exists():
        raise FileExistsError("preserve frozen bisection campaign ledger")
    ledger = {"status": "RUNNING_GEOMETRY", "rows": [], "completed_boxes": 0,
              "discovery_scope": "generic parent + frozen score-only candidates + generic norm-ten orbit table only"}
    checkpoint(ledger_path, ledger)
    # All 48 map rosters are frozen before the first point-search worker.
    for index, candidate in enumerate(protocol["rows"]):
        directory = LOCAL / candidate["id"]
        directory.mkdir(exist_ok=False)
        entry = {"id": candidate["id"], "status": "RUNNING_GEOMETRY", "height_bin": candidate["height_bin"],
                 "stratum": candidate["stratum"], "stages": []}
        ledger["rows"].append(entry)
        checkpoint(ledger_path, ledger)
        supervision = run([SAGE, str(CAS / "prepare_det1092_bisection_orbit_seed.sage"), "--index", str(index)],
                          limits=Limits(protocol["geometry_seconds"], protocol["rss_bytes"]), log_path=directory / "seed.log",
                          checkpoint_path=directory / "seed.supervisor.json", cwd=ROOT)
        ok = supervision["outcome"] == "completed" and supervision["returncode"] == 0
        entry["stages"].append({"name": "seed", "status": "PASS" if ok else "FAILED_OR_CENSORED", "supervision": supervision})
        checkpoint(ledger_path, ledger)
        if not ok:
            raise ArithmeticError("preserve failed generic-only seed")
        supervision = run([SAGE, str(CAS / "prepare_det1092_bisection_orbit_maps.sage"), "--index", str(index)],
                          limits=Limits(protocol["geometry_seconds"], protocol["rss_bytes"]), log_path=directory / "maps.log",
                          checkpoint_path=directory / "maps.supervisor.json", cwd=ROOT)
        ok = supervision["outcome"] == "completed" and supervision["returncode"] == 0
        entry["stages"].append({"name": "maps", "status": "PASS" if ok else "FAILED_OR_CENSORED", "supervision": supervision})
        checkpoint(ledger_path, ledger)
        if not ok:
            raise ArithmeticError("preserve failed full-orbit geometry")
        entry["status"] = "MAPS_FROZEN"
        checkpoint(ledger_path, ledger)
        print("BISECTION MAPS FROZEN", candidate["id"], index + 1, "/", len(protocol["rows"]), flush=True)
    ledger["status"] = "RUNNING_FIXED_EQUAL_BUDGET_SEARCH"
    checkpoint(ledger_path, ledger)
    for index, candidate in enumerate(protocol["rows"]):
        entry, directory = ledger["rows"][index], LOCAL / candidate["id"]
        entry["status"] = "RUNNING"
        checkpoint(ledger_path, ledger)
        for name, command, seconds in (
            ("worker", [sys.executable, str(Path(__file__).resolve()), "worker", "--index", str(index)], protocol["worker_seconds"]),
            ("replay", [sys.executable, str(Path(__file__).resolve()), "replay", "--index", str(index)], protocol["replay_seconds"]),
        ):
            supervision = run(command, limits=Limits(seconds, protocol["rss_bytes"]), log_path=directory / (name + ".log"),
                              checkpoint_path=directory / (name + ".supervisor.json"), cwd=ROOT)
            ok = supervision["outcome"] == "completed" and supervision["returncode"] == 0
            entry["stages"].append({"name": name, "status": "PASS" if ok else "FAILED_OR_CENSORED", "supervision": supervision})
            checkpoint(ledger_path, ledger)
            if not ok:
                raise ArithmeticError("preserve failed/censored equal-budget bisection worker")
        certificate = certify(index)
        entry.update(status="PASS", rank_lower_bound=certificate["rank_lower_bound"], discovered_rank_gain=certificate["discovered_rank_gain"],
                     certification_sha256=sha(directory / "certification-ledger.json"), completed_boxes=CHARTS)
        ledger["completed_boxes"] += CHARTS
        checkpoint(ledger_path, ledger)
        print("BISECTION FIBRE COMPLETE", candidate["id"], "rank>=", entry["rank_lower_bound"], flush=True)
    ledger["status"] = "PASS"
    checkpoint(ledger_path, ledger)


def detach() -> None:
    campaign()
    record, log = LOCAL / "detached-launch.json", LOCAL / "launch.log"
    if record.exists() or log.exists():
        raise FileExistsError("preserve existing detached bisection campaign launch")
    with log.open("x") as stream:
        process = subprocess.Popen([sys.executable, str(Path(__file__).resolve()), "launch"], cwd=ROOT,
                                   stdout=stream, stderr=stream, start_new_session=True)
    checkpoint(record, {"schema": "elliptic-curves.det1092-bisection-orbit-detached-launch.v1", "status": "STARTED", "pid": process.pid,
                        "started_at_unix": int(time.time()), "scope": "frozen 48 fibres, 49 norm-ten bisection charts each"})
    print("STARTED", process.pid, flush=True)


def compare() -> None:
    """Open the old policy only after the new discovery campaign is terminal."""
    protocol, ledger = campaign(), read(LOCAL / "ledger.json")
    if ledger["status"] != "PASS" or len(ledger["rows"]) != len(protocol["rows"]):
        raise ArithmeticError("new bisection campaign is not terminal")
    old = read(OLD_LEDGER)
    if old["status"] != "PASS" or [row["id"] for row in old["rows"]] != [row["id"] for row in protocol["rows"]]:
        raise ArithmeticError("matched old equal-budget control is unavailable")
    new_gains = [row for row in ledger["rows"] if row["rank_lower_bound"] > INITIAL_RANK]
    old_gains = [row for row in old["rows"] if row["rank_lower_bound"] > INITIAL_RANK]
    if OUTPUT_COMPARISON.exists():
        raise FileExistsError("preserve bisection coverage comparison")
    checkpoint(OUTPUT_COMPARISON, {
        "schema": "elliptic-curves.det1092-bisection-orbit-holdout-comparison.v1", "status": "PASS_TERMINAL_MATCHED_COMPARISON",
        "candidate_count": len(protocol["rows"]), "equal_boxes_per_candidate": CHARTS,
        "old_policy": {"name": "2048-SHA-parity then top-49 specialized norms", "completed_boxes": old["completed_boxes"],
                       "gain_count": len(old_gains), "gain_rate": str(len(old_gains)) + "/" + str(len(old["rows"]))},
        "bisection_policy": {"name": "all generic norm-ten rational-bisection orbits then top-49 specialized norms", "completed_boxes": ledger["completed_boxes"],
                             "gain_count": len(new_gains), "gain_rate": str(len(new_gains)) + "/" + str(len(ledger["rows"])),
                             "gain_rows": [{"id": row["id"], "rank_lower_bound": row["rank_lower_bound"]} for row in new_gains]},
        "inputs": {rel(path): sha(path) for path in (LOCAL / "ledger.json", OLD_LEDGER, LOCAL / "protocol.json")},
        "boundary": "This matched finite portfolio comparison measures certified gain incidence under two fixed 49-chart policies only. It is not a population-level rate estimate or a claim that unchanged rank is an upper bound.",
    })


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("stage", choices=("freeze", "launch", "detach", "worker", "replay", "compare"))
    parser.add_argument("--index", type=int)
    args = parser.parse_args()
    if args.stage in ("worker", "replay"):
        if args.index is None:
            parser.error("--index is required for worker/replay")
        globals()[args.stage](args.index)
    else:
        globals()[args.stage]()

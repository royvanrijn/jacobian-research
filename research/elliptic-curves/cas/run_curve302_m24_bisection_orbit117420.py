#!/usr/bin/env python3
"""One frozen, execution-blind M24 pointed-quartic arm on curve 302.

The sole centre is the canonical norm-ten representative of the complete
degree-two translation orbit ``117420``.  The label was chosen after a
retrospective diagnostic, so this is a calibrated one-chart recovery test,
not evidence that it was a prospective ranking rule.  Once ``freeze`` has
written the protocol, geometry, point search, and replay receive neither a
public residual-point source nor the retrospective diagnostic.

All writes are new immutable local checkpoints or new generated certificates.
``launch`` is a convenience wrapper around the separately replayable stages.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import sys
from importlib.machinery import SourceFileLoader
from pathlib import Path

import certify_compact_r17_candidates as cert
import pari_pointed_backend as backend
from memory_rank_certificate import checked_rank
from pointed_quartic_search import PointedQuarticSearch
from research_runtime.memory_store import MemoryFactStore
from research_runtime.pointed_orbit_compression import compress
from research_runtime.quotient_only_reduction import QuotientOnlyReductionCache as ReductionCache
from research_runtime.search_state import raw_state
from research_runtime.store import checkpoint, digest
from research_runtime.supervisor import Limits, run


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
ART = ROOT / "artifacts/generated-results/elliptic-curves"
LOCAL = ROOT / "artifacts/local/elliptic-curves/curve302-m24-bisection-orbit117420-v1"

M24 = ART / "curve302_recovered_followup_wave_03_mod2_v1.json"
PARENT = ART / "curve302_recovered_mw17_parent_v1.json"
LATTICE = ART / "curve302_parent_degree2_multisection_lattice_v1.json"
ORBITS = ART / "curve302_parent_degree2_multisection_orbits_v1.tsv"
OLD_M17_MAPS = (
    ROOT
    / "artifacts/local/elliptic-curves/curve302-focused-point-exposure-v2"
    / "curve302-generic17/maps.json"
)
OUTPUT = ART / "curve302_m24_bisection_orbit117420_v1.json"
MOD2 = ART / "curve302_m24_bisection_orbit117420_mod2_v1.json"
MODL = ART / "curve302_m24_bisection_orbit117420_modl_v1.json"
SAGE = "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python"
ORBIT_MASK = 117420
M17_DIMENSION = 17
M24_DIMENSION = 24

READS: set[str] = set()


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rel(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def source_hashes() -> dict[str, str]:
    direct = [
        Path(__file__),
        CAS / "factor_free_pari_mapping.sage",
        CAS / "load_curve302_recovered_parent.sage",
        CAS / "audit_recorded_point_mod2_rank_v3.py",
        CAS / "audit_retained_cloud_modl.py",
        CAS / "verify_factor_free_rank.sage",
    ]
    return {
        **backend.sources(),
        **{rel(path): cert.hashed(path) for path in direct},
    }


def frozen_inputs() -> dict[str, str]:
    return {rel(path): sha(path) for path in (M24, PARENT, LATTICE, ORBITS)}


def read(path: Path):
    return cert.read(path)


def install_execution_guard() -> None:
    """Allow only frozen inputs and this arm's own checkpoints under artifacts.

    The compiled sieve cache is generated solely from the source hash checked by
    the backend, rather than an arithmetic input.  It is allowed so a worker
    can be reused without widening the data inputs.
    """

    approved = {path.resolve() for path in (M24, PARENT, LATTICE, ORBITS)}
    sieve_cache = ROOT / "artifacts/local/elliptic-curves/pointed-sieve-build"

    def guard(event, args):
        if event != "open" or not args or not isinstance(args[0], (str, bytes)):
            return
        path = Path(args[0]).resolve()
        if not path.is_relative_to(ROOT / "artifacts"):
            return
        if path in approved or path.is_relative_to(LOCAL) or path.is_relative_to(sieve_cache):
            READS.add(rel(path))
            return
        raise PermissionError("execution input rejected: " + rel(path))

    sys.addaudithook(guard)


def orbit_row() -> dict:
    lattice = read(LATTICE)
    if lattice["status"] != "PASS_COMPLETE_DEGREE2_TRANSLATION_QUOTIENT":
        raise ArithmeticError("complete degree-two quotient certificate is unavailable")
    if sha(ORBITS) != lattice["orbits_tsv_sha256"]:
        raise ArithmeticError("degree-two orbit table hash differs")
    rows = [line.split("\t") for line in ORBITS.read_text().splitlines()[1:]]
    matches = [row for row in rows if int(row[0]) == ORBIT_MASK]
    if len(matches) != 1:
        raise ArithmeticError("frozen degree-two orbit label is not unique")
    row = matches[0]
    if len(row) != 6 or row[1] != "rational" or int(row[2]) != 10:
        raise ArithmeticError("frozen orbit stopped being a norm-ten rational bisection")
    word = [int(value) for value in row[4].split()]
    if len(word) != M17_DIMENSION:
        raise ArithmeticError("degree-two orbit has a non-M17 representative")
    parent = read(PARENT)
    gram = parent["generic_height_gram"]
    norm = sum(word[i] * int(gram[i][j]) * word[j] for i in range(M17_DIMENSION) for j in range(M17_DIMENSION))
    if norm != 10:
        raise ArithmeticError("degree-two representative lost generic norm ten")
    return {
        "orbit_mask": ORBIT_MASK,
        "category": row[1],
        "minimum_generic_MW17_norm": int(row[2]),
        "degree_two_LLL_minimum_word": [int(value) for value in row[3].split()],
        "parent_MW17_representative": word,
        "degree_two_divisor_class": row[5],
        "parent_MW17_norm_rechecked": norm,
    }


def centre_from_orbit() -> dict:
    row = orbit_row()
    representative = row["parent_MW17_representative"] + [0] * (M24_DIMENSION - M17_DIMENSION)
    parity = sum((value & 1) << index for index, value in enumerate(representative))
    if len(representative) != M24_DIMENSION or parity != 95909:
        raise ArithmeticError("frozen orbit no longer gives the declared M24 parity")
    return {
        **row,
        "representative": representative,
        "M24_parity_mask": parity,
        "M24_extension_rule": "append seven zero coordinates to the certified M17 representative",
    }


def seed_from_m24() -> dict:
    cloud = read(M24)
    if cloud["status"] != "COMPLETE_DECLARED_FINITE_AUDIT" or cloud["rank_lower_bound"] != M24_DIMENSION:
        raise ArithmeticError("the frozen M24 lower-bound cloud changed")
    points = cloud["independent_points"]
    if len(points) != M24_DIMENSION or cloud["points"] != points:
        raise ArithmeticError("M24 cloud no longer contains exactly its certified basis")
    model = tuple(map(cert.F, cloud["curve"]))
    proof = cloud["rank_certificate"]
    actual = checked_rank(
        model,
        [tuple(map(cert.F, point)) for point in points],
        [signature["prime"] for signature in proof["signatures"]],
        proof["no_rational_2_torsion_prime"],
    )
    if digest(actual) != digest(proof):
        raise ArithmeticError("M24 independent rank certificate does not replay")
    return {
        "family": cloud["family"],
        "parameter": cloud["parameter"],
        "curve": cloud["curve"],
        "points": points,
        "generic_points": points[:M17_DIMENSION],
        "rank_certificate": proof,
        "input_cloud": rel(M24),
        "input_cloud_sha256": sha(M24),
    }


def protocol() -> dict:
    return {
        "schema": "elliptic-curves.curve302-m24-bisection-orbit117420.v1",
        "sources": source_hashes(),
        "inputs": frozen_inputs(),
        "initial_rank": M24_DIMENSION,
        "chart_count": 1,
        "height": 125000,
        "seconds_per_chart": 10,
        "rss_bytes": 2147483648,
        "gp_sha256": cert.hashed(Path("/usr/bin/gp")),
        "centre_selection": {
            "rule": "one predeclared complete-degree-two rational-bisection translation orbit; take its canonical exported M17 minimum and append seven M24 zeros",
            "frozen_orbit_mask": ORBIT_MASK,
            "calibration_status": "retrospectively calibrated; execution-blind",
            "minimum_change": "one additional centre/chart; the sealed 17-to-24 deep-centre policy is not modified",
            "selection_inputs": [rel(PARENT), rel(LATTICE), rel(ORBITS)],
        },
        "execution_blindness": {
            "allowed_artifact_inputs": [rel(path) for path in (M24, PARENT, LATTICE, ORBITS)],
            "forbidden_inputs": [
                "the residual-visibility diagnostic",
                "any public missing-point coordinate or public rank-31 point list",
                "post-search finite-rank outcomes",
            ],
            "enforcement": "geometry, worker, and replay install an artifact-read guard before their arithmetic imports; each serializes the accepted read paths",
        },
        "boundary": "A single retrospectively calibrated target-blind execution arm. A rank gain is an exact lower bound after independent finite certificates, not a prospective selector validation, an exact-rank proof, or an assertion about any other bisection orbit.",
    }


def assert_protocol() -> dict:
    value = read(LOCAL / "protocol.json")
    if value != protocol():
        raise ArithmeticError("frozen protocol or one of its source/input hashes changed")
    return value


def freeze() -> None:
    if LOCAL.exists() or OUTPUT.exists() or MOD2.exists() or MODL.exists():
        raise FileExistsError("preserve existing orbit117420 evidence")
    LOCAL.mkdir(parents=True)
    seed = seed_from_m24()
    centre = centre_from_orbit()
    checkpoint(LOCAL / "protocol.json", protocol())
    checkpoint(LOCAL / "seed.json", seed)
    checkpoint(LOCAL / "centre.json", centre)
    print("FROZEN CURVE302 M24 BISECTION ARM|orbit=117420|charts=1|rank=24", flush=True)


def parent_prefix_in_short_model():
    """Specialize the parent basis and apply the certified public short map."""

    from sage.all import QQ

    loader = SourceFileLoader("curve302_parent_loader", str(CAS / "load_curve302_recovered_parent.sage")).load_module()
    curve, basis, t0 = loader.load_curve302_recovered_parent(PARENT)
    if t0 != 0 or tuple(curve.a_invariants())[:3] != (1, 1, 1):
        raise ArithmeticError("parent specialization presentation changed")
    rows = []
    for point in basis:
        x, y = QQ(point[0](t0)), QQ(point[1](t0))
        # X=36x+15 and Y=108(2y+x+1), as certified for the literal 302 fibre.
        rows.append([str(cert.F(str(36 * x + 15))), str(cert.F(str(108 * (2 * y + x + 1))) )])
    return rows


def geometry() -> None:
    install_execution_guard()
    p = assert_protocol()
    out = LOCAL / "maps.json"
    if out.exists():
        raise FileExistsError("preserve frozen bisection map")
    seed, centre = read(LOCAL / "seed.json"), read(LOCAL / "centre.json")
    if parent_prefix_in_short_model() != seed["generic_points"]:
        raise ArithmeticError("M24 prefix is not the specialized certified M17 parent basis")
    mapper = SourceFileLoader("curve302_bisection_mapper", str(CAS / "factor_free_pari_mapping.sage")).load_module()
    mapping = mapper.mapping(tuple(map(cert.F, seed["curve"])), tuple(tuple(map(cert.F, point)) for point in seed["points"]), centre)
    if mapping["centre"] != centre or mapping["coordinate_policy"] != {"kind": "raw", "matrix": mapping["matrix"]}:
        raise ArithmeticError("frozen centre/map binding differs")
    checkpoint(out, {
        "schema": "elliptic-curves.curve302-m24-bisection-map.v1",
        "status": "COMPLETE_DECLARED_SINGLE_MAP",
        "protocol_sha256": sha(LOCAL / "protocol.json"),
        "seed_sha256": sha(LOCAL / "seed.json"),
        "centre_sha256": sha(LOCAL / "centre.json"),
        "parent_prefix_verified": True,
        "mapping": mapping,
    })
    checkpoint(LOCAL / "geometry-data-access.json", sorted(READS))
    print("FROZEN CURVE302 M24 BISECTION MAP|orbit=117420", flush=True)


def state_from_seed(seed, cache):
    state = raw_state(
        tuple(map(cert.F, seed["curve"])),
        tuple(tuple(map(cert.F, point)) for point in seed["points"]),
        cache=cache,
        prime_bound=1000,
    )
    if state.rank != M24_DIMENSION:
        raise ArithmeticError("M24 seed did not reproduce its certified independent rank")
    return state


def worker() -> None:
    install_execution_guard()
    p = assert_protocol()
    out = LOCAL / "result.json"
    if out.exists():
        raise FileExistsError("preserve completed point attempt")
    seed, centre, maps = read(LOCAL / "seed.json"), read(LOCAL / "centre.json"), read(LOCAL / "maps.json")
    if maps["status"] != "COMPLETE_DECLARED_SINGLE_MAP" or maps["protocol_sha256"] != sha(LOCAL / "protocol.json"):
        raise ArithmeticError("single map is not frozen against this protocol")
    cache = ReductionCache(MemoryFactStore())
    state = state_from_seed(seed, cache)
    mapping = maps["mapping"]
    search = PointedQuarticSearch(
        state=state,
        centre={"coefficients": centre["representative"]},
        coordinate_policy=mapping["coordinate_policy"],
    )
    transcript, points = backend.execute(search, mapping, p["height"], p["seconds_per_chart"], p["gp_sha256"])
    compression = compress(tuple(map(cert.F, seed["curve"])), state.basis, centre["representative"], points)
    for index in compression["kept_indices"]:
        state = state.adjoin(points[index], cache=cache)
    final = state.record()
    checkpoint(out, {
        "schema": "elliptic-curves.curve302-m24-bisection-worker.v1",
        "status": "COMPLETE_DECLARED_SINGLE_CHART" if transcript["status"] == "bounded_search_complete" else "CENSORED_SINGLE_CHART",
        "protocol_sha256": sha(LOCAL / "protocol.json"),
        "seed_sha256": sha(LOCAL / "seed.json"),
        "maps_sha256": sha(LOCAL / "maps.json"),
        "family": seed["family"],
        "parameter": seed["parameter"],
        "curve": seed["curve"],
        "generic_points": seed["generic_points"],
        "initial_state": state_from_seed(seed, ReductionCache(MemoryFactStore())).record(),
        "centre": centre,
        "mapping": mapping,
        "search": transcript,
        "charts": [{
            "index": 0,
            "centre": centre,
            "search": transcript,
            "admission_compression": compression,
        }],
        "admission_compression": compression,
        "final_state": final,
        "rank_lower_bound": state.rank,
    })
    checkpoint(LOCAL / "worker-data-access.json", sorted(READS))
    print("CURVE302 M24 BISECTION WORKER|status={}|rank={}".format(transcript["status"], state.rank), flush=True)


def replay() -> None:
    install_execution_guard()
    p = assert_protocol()
    seed, centre, maps, data = (read(LOCAL / name) for name in ("seed.json", "centre.json", "maps.json", "result.json"))
    if data["protocol_sha256"] != sha(LOCAL / "protocol.json") or data["maps_sha256"] != sha(LOCAL / "maps.json"):
        raise ArithmeticError("worker transcript is not bound to the frozen execution")
    cache = ReductionCache(MemoryFactStore())
    state = state_from_seed(seed, cache)
    initial = state.record()
    if data["initial_state"] != initial or data["centre"] != centre or data["mapping"] != maps["mapping"]:
        raise ArithmeticError("worker initial state or exact map differs")
    search = PointedQuarticSearch(state=state, centre={"coefficients": centre["representative"]}, coordinate_policy=maps["mapping"]["coordinate_policy"])
    points = backend.replay(search, maps["mapping"], data["search"])
    compression = compress(tuple(map(cert.F, seed["curve"])), state.basis, centre["representative"], points)
    if compression != data["admission_compression"]:
        raise ArithmeticError("replayed involution compression differs")
    for index in compression["kept_indices"]:
        state = state.adjoin(points[index], cache=cache)
    if state.record() != data["final_state"] or state.rank != data["rank_lower_bound"]:
        raise ArithmeticError("replayed point admissions differ")
    if data["search"]["status"] != "bounded_search_complete":
        raise ArithmeticError("censored single chart cannot be called a completed experiment")
    checkpoint(LOCAL / "replay.json", {
        "schema": "elliptic-curves.curve302-m24-bisection-replay.v1",
        "status": "PASS_EXACT_MAP_AND_POINT_REPLAY",
        "result_sha256": sha(LOCAL / "result.json"),
        "rank_lower_bound": state.rank,
    })
    checkpoint(LOCAL / "replay-data-access.json", sorted(READS))
    print("REPLAYED CURVE302 M24 BISECTION|rank={}".format(state.rank), flush=True)


def certificate_stage(name: str, command: list[str], seconds: int, ledger: dict) -> None:
    record = run(
        command,
        limits=Limits(seconds, 2147483648),
        log_path=LOCAL / (name + ".log"),
        checkpoint_path=LOCAL / (name + ".supervisor.json"),
        cwd=ROOT,
    )
    passed = record["outcome"] == "completed" and record["returncode"] == 0
    ledger["stages"].append({"name": name, "status": "PASS" if passed else "FAILED_OR_CENSORED", "supervision": record})
    checkpoint(LOCAL / "certification-ledger.json", ledger)
    if not passed:
        ledger["status"] = "FAILED_OR_CENSORED"
        checkpoint(LOCAL / "certification-ledger.json", ledger)
        raise ArithmeticError("certificate stage failed: " + name)


def certify() -> None:
    assert_protocol()
    if not (LOCAL / "replay.json").exists() or MOD2.exists() or MODL.exists() or (LOCAL / "certification-ledger.json").exists():
        raise FileExistsError("preserve existing certification evidence")
    result = LOCAL / "result.json"
    ledger = {"schema": "elliptic-curves.curve302-m24-bisection-certification.v1", "status": "RUNNING", "stages": []}
    checkpoint(LOCAL / "certification-ledger.json", ledger)
    certificate_stage("mod2-build", [sys.executable, str(CAS / "audit_recorded_point_mod2_rank_v3.py"), "--input", str(result), "--input-sha256", sha(result), "--output", str(MOD2), "--prime-bound", "997"], 180, ledger)
    certificate_stage("mod2-check", [sys.executable, str(CAS / "audit_recorded_point_mod2_rank_v3.py"), "--check", str(MOD2)], 180, ledger)
    certificate_stage("modl-build", [sys.executable, str(CAS / "audit_retained_cloud_modl.py"), "--input", str(MOD2), "--output", str(MODL)], 180, ledger)
    certificate_stage("modl-check", [sys.executable, str(CAS / "audit_retained_cloud_modl.py"), "--check", str(MODL)], 180, ledger)
    copied = LOCAL / "independent-rank"
    copied.mkdir()
    verifier = copied / "verify_factor_free_rank.sage"
    cloud = copied / MOD2.name
    shutil.copy2(CAS / "verify_factor_free_rank.sage", verifier)
    shutil.copy2(MOD2, cloud)
    checkpoint(copied / "protocol.json", {"files": {path.name: sha(path) for path in (verifier, cloud)}, "seconds": 180, "rss_bytes": 2147483648})
    certificate_stage("independent-rank", [SAGE, str(verifier), "--input", str(cloud)], 180, ledger)
    mod2 = read(MOD2)
    modl = read(MODL)
    if mod2["rank_lower_bound"] < 24 or any(audit["finite_column_rank"] < 25 for audit in modl["audits"]):
        raise ArithmeticError("the frozen arm did not certify rank at least 25 in each odd-modulus descent audit")
    ledger.update(status="PASS", rank_lower_bound=mod2["rank_lower_bound"], odd_modulus_ranks={str(audit["modulus"]): audit["finite_column_rank"] for audit in modl["audits"]})
    checkpoint(LOCAL / "certification-ledger.json", ledger)
    print("CERTIFIED CURVE302 M24 BISECTION|rank>={}".format(mod2["rank_lower_bound"]), flush=True)


def report() -> None:
    p = assert_protocol()
    if OUTPUT.exists():
        raise FileExistsError("preserve immutable generated report")
    needed = [LOCAL / name for name in ("seed.json", "centre.json", "maps.json", "result.json", "replay.json", "certification-ledger.json")]
    if any(not path.exists() for path in needed) or not MOD2.exists() or not MODL.exists():
        raise ArithmeticError("complete execution/replay/certificates are required before reporting")
    old = read(OLD_M17_MAPS)
    centre, maps, result, replay, certificate, mod2, modl = (read(path) for path in (LOCAL / "centre.json", LOCAL / "maps.json", LOCAL / "result.json", LOCAL / "replay.json", LOCAL / "certification-ledger.json", MOD2, MODL))
    parity = centre["M24_parity_mask"]
    if any(row["parity"] == parity for row in old["sample"]) or any(row["parity"] == parity for row in old["centres"]):
        raise ArithmeticError("old M17 deep-centre arm unexpectedly contained this parity")
    if parity >> 22:
        raise ArithmeticError("padded M17 parity unexpectedly survives the M22-floor exclusion")
    odd_rank = min(audit["finite_column_rank"] for audit in modl["audits"])
    if odd_rank < 25 or certificate["status"] != "PASS":
        raise ArithmeticError("report requires a certified rank-25 outcome")
    report_data = {
        "schema": "elliptic-curves.curve302-m24-bisection-orbit117420-report.v1",
        "status": "PASS_CERTIFIED_RANK_AT_LEAST_25",
        "inputs": {
            **frozen_inputs(),
            rel(OLD_M17_MAPS): sha(OLD_M17_MAPS),
            **{rel(path): sha(path) for path in (*needed, MOD2, MODL)},
        },
        "sources": source_hashes(),
        "frozen_rule": p["centre_selection"],
        "exact_half_lattice_degree_two_reconstruction": {
            **centre,
            "pointed_quartic_centre_short_model": result["search"]["base_point"],
            "reduced_quartic_coefficients_ascending": result["search"]["coefficients"],
            "raw_coordinate_matrix": maps["mapping"]["matrix"],
            "quartic_coefficient_bits": result["search"]["maximum_coefficient_bits"],
        },
        "why_sealed_17_to_24_policy_excluded_it": {
            "initial_M17_deep_sample_contains_parity": False,
            "initial_M17_selected_centres_contain_parity": False,
            "initial_policy": "2048 SHA masks; select the 49 largest specialized rounded-height norms",
            "later_policy": "each follow-up wave requires a nonzero coordinate above the preceding subgroup rank",
            "padded_M24_parity_mask": parity,
            "padded_parity_shift_right_22": parity >> 22,
            "consequence": "This M17-supported parity was not sampled initially and is arithmetically ineligible in the M24 wave because the M22 floor rejects every M17-supported word.",
        },
        "execution_oracle_audit": {
            "rule_is_retrospectively_calibrated": True,
            "execution_is_target_blind": True,
            "allowed_inputs": p["execution_blindness"]["allowed_artifact_inputs"],
            "geometry_reads": read(LOCAL / "geometry-data-access.json"),
            "worker_reads": read(LOCAL / "worker-data-access.json"),
            "replay_reads": read(LOCAL / "replay-data-access.json"),
            "excluded_from_execution": p["execution_blindness"]["forbidden_inputs"],
        },
        "completed_search": {
            "charts": 1,
            "height": p["height"],
            "seconds_per_chart": p["seconds_per_chart"],
            "search_status": result["search"]["status"],
            "returned_finite_curve_points": result["search"]["finite_curve_points"],
            "rank_before": M24_DIMENSION,
            "rank_after_worker": result["rank_lower_bound"],
            "rank_after_exact_replay": replay["rank_lower_bound"],
            "mod_2_certified_rank_lower_bound": mod2["rank_lower_bound"],
            "certified_rank_lower_bound_by_odd_descent": odd_rank,
            "mod_3_5_ranks": {str(audit["modulus"]): audit["finite_column_rank"] for audit in modl["audits"]},
        },
        "boundary": "This recovers an additional direction from an execution that did not read its point. The one-chart rule was nonetheless chosen retrospectively from orbit117420, so it is not a prospective visibility ranking validation. It gives a certified lower bound only and changes neither the sealed historical 17-to-24 policy nor the curve's already larger public rank lower bound.",
        "reproducing_command": "python3 elliptic-curves/cas/run_curve302_m24_bisection_orbit117420.py launch",
    }
    OUTPUT.write_text(json.dumps(report_data, indent=2, sort_keys=True) + "\n")
    print("REPORTED CURVE302 M24 BISECTION|rank>={}".format(mod2["rank_lower_bound"]), flush=True)


def launch_stage(name: str, command: list[str], seconds: int, ledger: dict) -> None:
    record = run(command, limits=Limits(seconds, 2147483648), log_path=LOCAL / (name + ".log"), checkpoint_path=LOCAL / (name + ".supervisor.json"), cwd=ROOT)
    passed = record["outcome"] == "completed" and record["returncode"] == 0
    ledger["stages"].append({"name": name, "status": "PASS" if passed else "FAILED_OR_CENSORED", "supervision": record})
    checkpoint(LOCAL / "ledger.json", ledger)
    if not passed:
        ledger["status"] = "FAILED_OR_CENSORED"
        checkpoint(LOCAL / "ledger.json", ledger)
        raise ArithmeticError("launch stage failed: " + name)


def launch() -> None:
    freeze()
    ledger = {"schema": "elliptic-curves.curve302-m24-bisection-orbit117420-launch.v1", "status": "RUNNING", "stages": []}
    checkpoint(LOCAL / "ledger.json", ledger)
    for name, command, seconds in (
        ("geometry", [SAGE, str(Path(__file__).resolve()), "geometry"], 180),
        ("worker", [sys.executable, str(Path(__file__).resolve()), "worker"], 120),
        ("replay", [sys.executable, str(Path(__file__).resolve()), "replay"], 120),
        ("certify", [sys.executable, str(Path(__file__).resolve()), "certify"], 900),
        ("report", [sys.executable, str(Path(__file__).resolve()), "report"], 120),
    ):
        launch_stage(name, command, seconds, ledger)
    ledger.update(status="PASS", report=rel(OUTPUT), report_sha256=sha(OUTPUT))
    checkpoint(LOCAL / "ledger.json", ledger)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("freeze", "geometry", "worker", "replay", "certify", "report", "launch"))
    args = parser.parse_args()
    globals()[args.stage]()

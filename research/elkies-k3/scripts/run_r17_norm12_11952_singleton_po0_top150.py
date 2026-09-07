#!/usr/bin/env python3
"""Run a checkpointed exact isolated-branch P.O=0 screen on top singleton twists."""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
LOCAL = ROOT / "artifacts/local/elkies-k3/twist-polynomial-sections"
CENSUS = GENERATED / (
    "elkies-k3-r17-norm12-11952-singleton-twist-census-top200-v1.json"
)
BISECTIONS = GENERATED / (
    "elkies-k3-r17-norm12-11952-alternate-bisections-full-v1.json"
)
MODEL = GENERATED / (
    "elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json"
)
EXPORTER = ROOT / "elkies-k3/scripts/export_elkies_2026_twist_polynomial_sections_modp.sage"
DERIVER = ROOT / "elkies-k3/scripts/derive_elkies_2026_singleton_twist_section.sage"
ENUMERATOR = ROOT / "elkies-k3/scripts/bruteforce_twist_polynomial_sections_modp.cpp"
BRUTE_WRAPPER = ROOT / "elkies-k3/scripts/run_twist_polynomial_sections_bruteforce.py"
LIFTER = ROOT / "elkies-k3/scripts/lift_r17_norm12_direct_singleton_po0_bruteforce.sage"
OUTPUT = GENERATED / (
    "elkies-k3-r17-norm12-11952-singleton-twist-po0-top150-campaign-v1.json"
)
SCHEMA = "elkies-k3.r17-norm12-11952-singleton-po0-topn-campaign.v1"
COMPLETE_EXPORT = "PASS_EXACT_MODP_REDUCED_POLYNOMIAL_SECTION_EXPORT"
COMPLETE_BRUTE = "PASS_EXHAUSTIVE_FINITE_FIELD_ENUMERATION_OF_EXPORTED_BLOCKS"
COMPLETE_LIFT = "PASS_EXACT_BRUTEFORCE_SEED_HENSEL_AUDIT"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def load(path: Path):
    return json.loads(path.read_text())


def pins_current(payload) -> bool:
    return all(
        (ROOT / name).is_file() and digest(ROOT / name) == expected
        for name, expected in payload.get("inputs", {}).items()
    )


def paths(orbit: str, prime: int) -> tuple[Path, Path, Path, Path]:
    export = LOCAL / f"direct-singleton-alternate-orbit-{orbit}/p{prime}/export.json"
    section = GENERATED / (
        "elkies-k3-r17-norm12-11952-singleton-twist-section-"
        f"orbit{orbit}-p{prime}-v1.json"
    )
    brute = GENERATED / (
        "elkies-k3-r17-norm12-11952-singleton-twist-po0-"
        f"orbit{orbit}-p{prime}-v1.json"
    )
    lift = GENERATED / (
        "elkies-k3-r17-norm12-11952-singleton-twist-po0-"
        f"orbit{orbit}-p{prime}-hensel-v1.json"
    )
    return export, section, brute, lift


def valid_brute(orbit: str, prime: int, path: Path) -> bool:
    if not path.is_file():
        return False
    payload = load(path)
    candidate = payload.get("candidate", {})
    return (
        payload.get("status") == COMPLETE_BRUTE
        and int(payload.get("prime", 0)) == prime
        and candidate.get("key") == f"alternate-orbit-{orbit}"
        and candidate.get("chi") == 3
        and candidate.get("x_degree_bound") == 6
        and candidate.get("y_degree_bound") == 9
        and len(payload.get("known_section_match_indices", [])) == 1
        and payload.get("inputs", {}).get(relative(ENUMERATOR)) == digest(ENUMERATOR)
        and pins_current(payload)
    )


def run(command: list[str]) -> tuple[bool, str]:
    completed = subprocess.run(
        command,
        cwd=ROOT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    return completed.returncode == 0, completed.stdout


def prepare_shell(
    sage: str,
    label: str,
    orbit: str,
    prime: int,
) -> tuple[Path, Path, Path] | None:
    export, section, brute, unused_lift = paths(orbit, prime)
    ok, exporter_output = run(
        [
            sage,
            "-python",
            str(EXPORTER),
            "--direct-label",
            label,
            "--prime",
            str(prime),
            "--bisections",
            str(BISECTIONS),
            "--model",
            str(MODEL),
        ]
    )
    if not ok:
        print(
            f"TOP150SKIP|label={label}|p={prime}|stage=export|"
            f"tail={exporter_output.splitlines()[-1] if exporter_output.splitlines() else 'error'}",
            flush=True,
        )
        return None
    if load(export).get("status") != COMPLETE_EXPORT:
        return None
    ok, deriver_output = run(
        [
            sage,
            "-python",
            str(DERIVER),
            "--direct-label",
            label,
            "--prime",
            str(prime),
            "--bisections",
            str(BISECTIONS),
            "--model",
            str(MODEL),
            "--export",
            str(export),
            "--output",
            str(section),
        ]
    )
    if not ok:
        print(
            f"TOP150SKIP|label={label}|p={prime}|stage=known|"
            f"tail={deriver_output.splitlines()[-1] if deriver_output.splitlines() else 'error'}",
            flush=True,
        )
        return None
    ok, brute_output = run(
        [
            sys.executable,
            str(BRUTE_WRAPPER),
            "--export",
            str(export),
            "--known-section",
            str(section),
            "--output",
            str(brute),
        ]
    )
    if not ok or not valid_brute(orbit, prime, brute):
        print(brute_output, flush=True)
        raise ArithmeticError(f"failed exhaustive brute shell for {label} at p={prime}")
    return export, section, brute


def obstruction_exponent(message: str) -> int:
    match = re.fullmatch(r"no lift modulo (?:p|\d+)\^(\d+)", message)
    if match is None:
        raise ArithmeticError(f"unexpected obstruction message {message}")
    return int(match.group(1))


def write_manifest(path: Path, result) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")


def verify_manifest(path: Path) -> None:
    payload = load(path)
    if payload.get("schema") != SCHEMA:
        raise ValueError("unexpected singleton campaign schema")
    for name, expected in payload["inputs"].items():
        resolved = ROOT / name
        if not resolved.is_file() or digest(resolved) != expected:
            raise ArithmeticError(f"campaign input digest changed: {name}")
    census = load(CENSUS)
    candidates = census["top"][: payload["candidate_limit"]]
    excluded_name = payload.get("excluded_prime_manifest")
    excluded_primes = {}
    if excluded_name is not None:
        excluded_payload = load(ROOT / excluded_name)
        if excluded_payload.get("schema") != SCHEMA:
            raise ArithmeticError("excluded-prime campaign schema changed")
        excluded_primes = {
            str(item["label"]): int(item["prime"])
            for item in excluded_payload["records"]
        }
    obstruction_histogram = Counter()
    total_tested = 0
    total_solutions = 0
    total_extras = 0
    unresolved = []
    exact_new_sections = []
    for rank, (record, candidate_record) in enumerate(
        zip(payload["records"], candidates, strict=True), start=1
    ):
        if (
            record["discovery_rank"] != rank
            or record["label"] != candidate_record["label"]
            or record["orbit_hex"] != candidate_record["orbit_hex"]
            or record.get("excluded_prime") != excluded_primes.get(record["label"])
            or (
                record.get("excluded_prime") is not None
                and record["prime"] == record["excluded_prime"]
            )
        ):
            raise ArithmeticError(f"campaign rank {rank} no longer matches the census")
        for field in ("export", "known_section", "bruteforce", "hensel"):
            name = record.get(field)
            if name is None:
                continue
            resolved = ROOT / name
            expected = record[f"{field}_sha256"]
            if not resolved.is_file() or digest(resolved) != expected:
                raise ArithmeticError(f"campaign output digest changed: {name}")
        export = load(ROOT / record["export"])
        section = load(ROOT / record["known_section"])
        brute_path = ROOT / record["bruteforce"]
        brute = load(brute_path)
        prime = int(record["prime"])
        label = str(record["label"])
        orbit = label.removeprefix("alternate-orbit-")
        if (
            export.get("status") != COMPLETE_EXPORT
            or export.get("candidate", {}).get("key") != label
            or int(export.get("prime", 0)) != prime
            or section.get("status") != "PASS_EXACT_DESCENDED_SINGLETON_TWIST_SECTION"
            or section.get("modular_identification", {}).get("export_sha256")
            != digest(ROOT / record["export"])
            or section.get("modular_identification", {}).get(
                "exported_equation_residuals_zero"
            )
            is not True
            or not valid_brute(orbit, prime, brute_path)
            or brute.get("export_status") != COMPLETE_EXPORT
            or brute.get("inputs", {}).get(record["export"])
            != digest(ROOT / record["export"])
            or brute.get("inputs", {}).get(record["known_section"])
            != digest(ROOT / record["known_section"])
        ):
            raise ArithmeticError(f"campaign rank {rank} exact shell metadata changed")
        solutions = brute["solutions"]
        known = {int(index) for index in brute["known_section_match_indices"]}
        extras = [
            index
            for index, solution in enumerate(solutions)
            if int(solution["full_shell_tangent_rank"]) == 8 and index not in known
        ]
        if (
            len(known) != 1
            or record["known_solution_index"] != min(known)
            or record["isolated_extra_solution_indices"] != extras
            or record["x_polynomials_tested"]
            != int(brute["enumeration"]["x_polynomials_tested"])
            or record["representative_sign_solution_count"] != len(solutions)
            or brute["enumeration"]["representative_sign_solutions"] != len(solutions)
        ):
            raise ArithmeticError(f"campaign rank {rank} solution counts changed")

        recorded_obstructions = {
            int(item["solution_index"]): int(item["first_impossible_exponent"])
            for item in record["exact_local_obstructions"]
        }
        recomputed_unresolved = []
        recomputed_exact = []
        if extras:
            hensel = load(ROOT / record["hensel"])
            if (
                hensel.get("status") != COMPLETE_LIFT
                or hensel.get("inputs", {}).get(record["bruteforce"])
                != digest(brute_path)
                or not pins_current(hensel)
            ):
                raise ArithmeticError(f"campaign rank {rank} Hensel pins changed")
            lifts = {int(item["solution_index"]): item for item in hensel["lifts"]}
            if not set(extras).issubset(lifts):
                raise ArithmeticError(f"campaign rank {rank} omitted an isolated branch")
            for index in extras:
                item = lifts[index]
                obstruction = item.get("exact_local_obstruction")
                exact = item.get("exact_rational_reconstruction")
                if obstruction:
                    exponent = obstruction_exponent(obstruction)
                    obstruction_histogram[exponent] += 1
                    if recorded_obstructions.get(index) != exponent:
                        raise ArithmeticError(
                            f"campaign rank {rank} obstruction exponent changed"
                        )
                elif exact and exact.get("literal_curve_substitution") is True:
                    recomputed_exact.append(index)
                    exact_new_sections.append(
                        {
                            "discovery_rank": rank,
                            "label": label,
                            "prime": prime,
                            "solution_index": index,
                        }
                    )
                else:
                    recomputed_unresolved.append(index)
                    unresolved.append(
                        {
                            "discovery_rank": rank,
                            "label": label,
                            "prime": prime,
                            "solution_index": index,
                        }
                    )
        elif record.get("hensel") is not None or recorded_obstructions:
            raise ArithmeticError(f"campaign rank {rank} has a spurious Hensel record")
        if (
            set(recorded_obstructions) != set(extras) - set(recomputed_unresolved) - set(recomputed_exact)
            or record["unresolved_isolated_solution_indices"] != recomputed_unresolved
            or record["exact_new_section_indices"] != recomputed_exact
        ):
            raise ArithmeticError(f"campaign rank {rank} lift classification changed")
        total_tested += record["x_polynomials_tested"]
        total_solutions += len(solutions)
        total_extras += len(extras)
    if len(payload["records"]) != payload["candidate_limit"]:
        raise ArithmeticError("campaign record count changed")
    expected_summary = {
        "completed_character_count": len(payload["records"]),
        "x_polynomials_tested_across_distinct_shells": total_tested,
        "representative_sign_solution_count": total_solutions,
        "isolated_extra_branch_count": total_extras,
        "exact_local_obstruction_count": sum(obstruction_histogram.values()),
        "first_impossible_exponent_histogram": {
            str(exponent): obstruction_histogram[exponent]
            for exponent in sorted(obstruction_histogram)
        },
        "unresolved_isolated_branch_count": len(unresolved),
        "exact_new_section_count": len(exact_new_sections),
        "unresolved_isolated_branches": unresolved,
        "exact_new_sections": exact_new_sections,
    }
    if payload["summary"] != expected_summary:
        raise ArithmeticError("campaign summary differs from exact shell replay")
    expected_status = (
        "PASS_EXACT_NEW_RATIONAL_SECTION_CANDIDATE"
        if exact_new_sections
        else "INCOMPLETE_PADIC_SURVIVORS_REQUIRE_FURTHER_LIFTING"
        if unresolved
        else "PASS_BOUNDED_TOPN_ALL_ISOLATED_EXTRA_BRANCHES_LOCALLY_OBSTRUCTED"
    )
    if payload["status"] != expected_status:
        raise ArithmeticError("campaign status differs from exact shell replay")
    print(
        "R17SINGLETONPO0TOPNCHECK"
        f"|characters={len(payload['records'])}"
        f"|isolated_extras={payload['summary']['isolated_extra_branch_count']}"
        f"|unresolved={payload['summary']['unresolved_isolated_branch_count']}"
        f"|status={payload['status']}",
        flush=True,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate-limit", type=int, default=150)
    parser.add_argument("--primes", type=int, nargs="+", default=[17, 23, 29, 31, 37, 41])
    parser.add_argument(
        "--exclude-primes-from",
        type=Path,
        help="campaign manifest whose chosen prime is excluded character by character",
    )
    parser.add_argument("--hensel-depth", type=int, default=128)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    output = args.output if args.output.is_absolute() else ROOT / args.output
    if args.check:
        verify_manifest(output)
        return
    if args.candidate_limit < 1 or args.candidate_limit > 200:
        parser.error("--candidate-limit must be in 1..200")
    if args.hensel_depth < 2:
        parser.error("--hensel-depth must be at least two")
    sage = shutil.which("sage")
    if sage is None:
        raise FileNotFoundError("the Sage launcher is required")

    census = load(CENSUS)
    candidates = census["top"][: args.candidate_limit]
    excluded_manifest_path = (
        args.exclude_primes_from.resolve() if args.exclude_primes_from is not None else None
    )
    excluded_primes = {}
    if excluded_manifest_path is not None:
        excluded_manifest = load(excluded_manifest_path)
        if excluded_manifest.get("schema") != SCHEMA:
            raise ValueError("--exclude-primes-from has the wrong campaign schema")
        excluded_primes = {
            str(record["label"]): int(record["prime"])
            for record in excluded_manifest["records"]
        }
        if not all(str(record["label"]) in excluded_primes for record in candidates):
            raise ValueError("excluded-prime manifest does not cover the requested prefix")
    inputs = (
        Path(__file__).resolve(),
        CENSUS,
        BISECTIONS,
        MODEL,
        EXPORTER,
        DERIVER,
        ENUMERATOR,
        BRUTE_WRAPPER,
        LIFTER,
        *((excluded_manifest_path,) if excluded_manifest_path is not None else ()),
    )
    base = {
        "schema": SCHEMA,
        "status": "RUNNING_CHECKPOINTED_EXACT_PO0_CAMPAIGN",
        "candidate_limit": args.candidate_limit,
        "prime_preference": args.primes,
        "hensel_depth": args.hensel_depth,
        "excluded_prime_manifest": (
            relative(excluded_manifest_path) if excluded_manifest_path is not None else None
        ),
        "inputs": {relative(path): digest(path) for path in inputs},
        "records": [],
        "summary": {},
        "reproducing_command": shlex.join(sys.argv),
        "proof_boundary": (
            "Each chosen finite-field shell is exhaustive for polynomial sections "
            "with P.O=0, deg(X)<=6, deg(Y)<=9. Exact prime-power obstruction excludes "
            "only the indicated isolated modular branch. Singular modular branches, "
            "sections with P.O>0, candidates beyond the declared census prefix, and "
            "non-singleton characters remain open."
        ),
    }
    obstruction_histogram = Counter()
    total_tested = 0
    total_solutions = 0
    total_extras = 0
    unresolved = []
    exact_new_sections = []

    for rank, candidate_record in enumerate(candidates, start=1):
        label = str(candidate_record["label"])
        orbit = label.removeprefix("alternate-orbit-")
        excluded_prime = excluded_primes.get(label)
        candidate_primes = [prime for prime in args.primes if prime != excluded_prime]
        if not candidate_primes:
            raise ValueError(f"all requested primes are excluded for {label}")
        chosen = None
        # Reuse a hash-current complete artifact before creating another shell.
        for prime in candidate_primes:
            unused_export, unused_section, brute, unused_lift = paths(orbit, prime)
            if valid_brute(orbit, prime, brute):
                existing = load(brute)
                export_names = [
                    name for name in existing["inputs"] if name.endswith("/export.json")
                ]
                section_names = [
                    name
                    for name in existing["inputs"]
                    if "singleton-twist-section" in name
                ]
                if len(export_names) != 1 or len(section_names) != 1:
                    raise ArithmeticError(
                        f"cannot infer exact inputs for existing shell {label} p={prime}"
                    )
                export = ROOT / export_names[0]
                section = ROOT / section_names[0]
                chosen = (prime, export, section, brute)
                break
        if chosen is None:
            for prime in candidate_primes:
                prepared = prepare_shell(sage, label, orbit, prime)
                if prepared is not None:
                    export, section, brute = prepared
                    chosen = (prime, export, section, brute)
                    break
        if chosen is None:
            raise ArithmeticError(f"no usable complete shell found for {label}")
        prime, export, section, brute = chosen
        brute_payload = load(brute)
        solutions = brute_payload["solutions"]
        known = {int(index) for index in brute_payload["known_section_match_indices"]}
        extras = [
            index
            for index, solution in enumerate(solutions)
            if int(solution["full_shell_tangent_rank"]) == 8 and index not in known
        ]

        lift_path = paths(orbit, prime)[3]
        lift_payload = None
        if extras:
            reuse_lift = False
            if lift_path.is_file():
                possible = load(lift_path)
                by_index = {int(item["solution_index"]): item for item in possible.get("lifts", [])}
                reuse_lift = (
                    possible.get("status") == COMPLETE_LIFT
                    and possible.get("inputs", {}).get(relative(brute)) == digest(brute)
                    and set(extras).issubset(by_index)
                    and pins_current(possible)
                )
                if reuse_lift:
                    lift_payload = possible
            if not reuse_lift:
                command = [
                    sage,
                    "-python",
                    str(LIFTER),
                    "--bruteforce",
                    str(brute),
                    "--bisections",
                    str(BISECTIONS),
                    "--model",
                    str(MODEL),
                    "--hensel-depth",
                    str(args.hensel_depth),
                    "--output",
                    str(lift_path),
                ]
                for index in extras:
                    command.extend(("--solution-index", str(index)))
                ok, lift_output = run(command)
                if not ok:
                    print(lift_output, flush=True)
                    raise ArithmeticError(f"Hensel audit failed for {label} at p={prime}")
                lift_payload = load(lift_path)

        obstruction_records = []
        unresolved_indices = []
        exact_indices = []
        if lift_payload is not None:
            lifts = {int(item["solution_index"]): item for item in lift_payload["lifts"]}
            for index in extras:
                item = lifts[index]
                obstruction = item.get("exact_local_obstruction")
                exact = item.get("exact_rational_reconstruction")
                if obstruction:
                    exponent = obstruction_exponent(obstruction)
                    obstruction_histogram[exponent] += 1
                    obstruction_records.append(
                        {"solution_index": index, "first_impossible_exponent": exponent}
                    )
                elif exact and exact.get("literal_curve_substitution") is True:
                    exact_indices.append(index)
                    exact_new_sections.append(
                        {"discovery_rank": rank, "label": label, "prime": prime, "solution_index": index}
                    )
                else:
                    unresolved_indices.append(index)
                    unresolved.append(
                        {"discovery_rank": rank, "label": label, "prime": prime, "solution_index": index}
                    )

        record = {
            "discovery_rank": rank,
            "label": label,
            "orbit_hex": candidate_record["orbit_hex"],
            "prime": prime,
            "excluded_prime": excluded_prime,
            "x_polynomials_tested": int(brute_payload["enumeration"]["x_polynomials_tested"]),
            "representative_sign_solution_count": len(solutions),
            "known_solution_index": min(known),
            "isolated_extra_solution_indices": extras,
            "exact_local_obstructions": obstruction_records,
            "unresolved_isolated_solution_indices": unresolved_indices,
            "exact_new_section_indices": exact_indices,
            "export": relative(export),
            "export_sha256": digest(export),
            "known_section": relative(section),
            "known_section_sha256": digest(section),
            "bruteforce": relative(brute),
            "bruteforce_sha256": digest(brute),
            "hensel": relative(lift_path) if extras else None,
            "hensel_sha256": digest(lift_path) if extras else None,
        }
        base["records"].append(record)
        total_tested += record["x_polynomials_tested"]
        total_solutions += len(solutions)
        total_extras += len(extras)
        base["summary"] = {
            "completed_character_count": len(base["records"]),
            "x_polynomials_tested_across_distinct_shells": total_tested,
            "representative_sign_solution_count": total_solutions,
            "isolated_extra_branch_count": total_extras,
            "exact_local_obstruction_count": sum(obstruction_histogram.values()),
            "first_impossible_exponent_histogram": {
                str(exponent): obstruction_histogram[exponent]
                for exponent in sorted(obstruction_histogram)
            },
            "unresolved_isolated_branch_count": len(unresolved),
            "exact_new_section_count": len(exact_new_sections),
            "unresolved_isolated_branches": unresolved,
            "exact_new_sections": exact_new_sections,
        }
        write_manifest(output, base)
        print(
            "TOP150PROGRESS"
            f"|rank={rank}/{args.candidate_limit}|label={label}|p={prime}"
            f"|tested={record['x_polynomials_tested']}|solutions={len(solutions)}"
            f"|extras={len(extras)}|unresolved={len(unresolved_indices)}",
            flush=True,
        )

    base["status"] = (
        "PASS_EXACT_NEW_RATIONAL_SECTION_CANDIDATE"
        if exact_new_sections
        else "INCOMPLETE_PADIC_SURVIVORS_REQUIRE_FURTHER_LIFTING"
        if unresolved
        else "PASS_BOUNDED_TOPN_ALL_ISOLATED_EXTRA_BRANCHES_LOCALLY_OBSTRUCTED"
    )
    write_manifest(output, base)
    print(
        "R17SINGLETONPO0TOPN"
        f"|characters={len(base['records'])}|tested={total_tested}"
        f"|isolated_extras={total_extras}|obstructed={sum(obstruction_histogram.values())}"
        f"|unresolved={len(unresolved)}|exact_new={len(exact_new_sections)}"
        f"|status={base['status']}",
        flush=True,
    )


if __name__ == "__main__":
    main()

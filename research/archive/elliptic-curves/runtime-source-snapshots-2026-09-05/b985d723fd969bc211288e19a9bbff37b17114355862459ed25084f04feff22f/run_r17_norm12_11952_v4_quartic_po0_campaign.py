#!/usr/bin/env python3
"""Run a checkpointed bounded P.O=0 sieve on shortlisted V4 product twists.

status: ACTIVE_SEARCH
claim: bounded modular polynomial-section search for declared quartic twists
inputs: exact alternate-Q80 V4 pair shortlist and direct model certificates
outputs: artifacts/generated-results/elkies-k3-r17-norm12-11952-v4-quartic-po0-campaign-v1.json
supersedes: none; consumes the intersection-one rational-base shortlist

The default is intentionally a one-pair, one-prime, one-system-group pilot.
Larger runs require explicit limits.  Modular solutions are discovery inputs;
empty finite-field systems and timeouts are not characteristic-zero rank
statements.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import shlex
import shutil
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[2]
SHORTLIST = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-norm12-11952-v4-pair-shortlist-64-v1.json"
)
RANK_SCREEN = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-norm12-11952-v4-base-rank-screen-v1.json"
)
BISECTIONS = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-norm12-11952-alternate-bisections-cheapest-1024-v1.json"
)
DIRECT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json"
)
EXPORTER = ROOT / "elkies-k3/scripts/export_elkies_2026_twist_polynomial_sections_modp.sage"
SOLVER = ROOT / "elkies-k3/scripts/run_elkies_2026_twist_polynomial_sections_msolve.py"
OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-norm12-11952-v4-quartic-po0-campaign-v1.json"
)
SCHEMA = "elkies-k3.r17-norm12-11952-v4-quartic-po0-campaign.v1"


def digest(path: Path) -> str:
    result = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            result.update(block)
    return result.hexdigest()


def display_path(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(ROOT))
    except ValueError:
        return str(path.resolve())


def candidate_tag(pair_key: str) -> str:
    return f"direct-product-{pair_key.replace(':', '--')}"


def verify_manifest(path: Path) -> None:
    payload = json.loads(path.read_text())
    if payload.get("schema") != SCHEMA:
        raise ValueError("unexpected campaign schema")
    for input_path, expected in payload["inputs"].items():
        resolved = ROOT / input_path
        if digest(resolved) != expected:
            raise ArithmeticError(f"campaign input digest changed: {resolved}")
    for job in payload["jobs"]:
        export = ROOT / job["export"]
        if digest(export) != job["export_sha256"]:
            raise ArithmeticError(f"export digest changed: {export}")
        summary_name = job.get("summary")
        if summary_name:
            summary = ROOT / summary_name
            if digest(summary) != job["summary_sha256"]:
                raise ArithmeticError(f"solver summary digest changed: {summary}")
    print(
        "ALTV4PO0CHECK|"
        f"jobs={len(payload['jobs'])}|status={payload['status']}|output={display_path(path)}",
        flush=True,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--shortlist", type=Path, default=SHORTLIST)
    parser.add_argument("--rank-screen", type=Path, default=RANK_SCREEN)
    parser.add_argument("--bisections", type=Path, default=BISECTIONS)
    parser.add_argument("--direct", type=Path, default=DIRECT)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--pair-limit", type=int, default=1)
    parser.add_argument(
        "--selection",
        choices=("exact-rank-one", "shortlist"),
        default="exact-rank-one",
        help="prefer certified rank-one base Jacobians or retain shortlist order",
    )
    parser.add_argument("--primes-per-pair", type=int, default=1)
    parser.add_argument("--export-only", action="store_true")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--max-groups", type=int, default=1)
    group.add_argument("--all-groups", action="store_true")
    parser.add_argument("--timeout", type=float, default=60.0)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.check:
        verify_manifest(args.output)
        return
    if (
        args.pair_limit < 1
        or args.primes_per_pair < 1
        or args.timeout <= 0
        or args.threads < 1
        or args.jobs < 1
        or (args.max_groups is not None and args.max_groups < 1)
    ):
        parser.error("all campaign limits must be positive")

    shortlist = json.loads(args.shortlist.read_text())
    if shortlist.get("schema") != "elkies-k3.r17-norm12-11952-v4-pair-shortlist.v1":
        raise ValueError("unexpected V4 shortlist schema")
    if shortlist.get("status") != "PASS_EXACT_BOUNDED_RATIONAL_V4_PAIR_SHORTLIST":
        raise ValueError("the V4 shortlist is not an exact rational-base certificate")
    if args.selection == "exact-rank-one":
        rank_screen = json.loads(args.rank_screen.read_text())
        if rank_screen.get("status") != "PASS_BOUNDED_EXACT_BASE_JACOBIAN_RANK_INTERVAL_SCREEN":
            raise ValueError("the base-Jacobian rank screen is not certified")
        exact_rank_one_keys = [
            row["pair_key"]
            for row in rank_screen["results"]
            if row.get("status") == "completed"
            and int(row["rank_lower_bound"]) == 1
            and int(row["rank_upper_bound"]) == 1
        ]
        by_key = {row["pair_key"]: row for row in shortlist["pairs"]}
        pairs = [by_key[key] for key in exact_rank_one_keys[: args.pair_limit]]
    else:
        rank_screen = None
        pairs = shortlist["pairs"][: args.pair_limit]
    if len(pairs) != args.pair_limit:
        raise ValueError("--pair-limit exceeds the exact shortlist")

    sage = shutil.which("sage")
    if sage is None:
        raise FileNotFoundError("the Sage launcher is required")
    jobs = []
    for pair in pairs:
        pair_key = str(pair["pair_key"])
        primes = pair["recommended_complete_sieve_primes"][: args.primes_per_pair]
        if len(primes) != args.primes_per_pair:
            raise ValueError(f"pair {pair_key} has too few recommended primes")
        for prime in primes:
            export_command = [
                sage,
                "-python",
                str(EXPORTER),
                "--direct-product-key",
                pair_key,
                "--prime",
                str(prime),
                "--bisections",
                str(args.bisections),
                "--model",
                str(args.direct),
            ]
            subprocess.run(export_command, cwd=ROOT, check=True)
            tag = candidate_tag(pair_key)
            export_path = (
                ROOT
                / "artifacts/local/elkies-k3/twist-polynomial-sections"
                / tag
                / f"p{prime}"
                / "export.json"
            )
            export = json.loads(export_path.read_text())
            if export["candidate"]["kind"] != "direct_product" or export["candidate"]["key"] != pair_key:
                raise ArithmeticError("export candidate does not match the shortlisted pair")
            if int(export["prime"]) != int(prime) or int(export["candidate"]["chi"]) != 4:
                raise ArithmeticError("export lost the declared prime or chi=4 degree box")

            job = {
                "pair_key": pair_key,
                "shortlist_rank": int(pair["shortlist_rank"]),
                "prime": int(prime),
                "export_command": shlex.join(export_command),
                "export": display_path(export_path),
                "export_sha256": digest(export_path),
                "export_status": export["status"],
                "polynomial_section_degree_bounds": {
                    "X": int(export["candidate"]["x_degree_bound"]),
                    "Y": int(export["candidate"]["y_degree_bound"]),
                },
                "total_blocks": len(export["systems"]),
            }
            if not args.export_only:
                summary_path = (
                    ROOT
                    / "artifacts/generated-results"
                    / f"elkies-k3-r17-norm12-11952-v4-po0-{tag}-p{prime}-msolve.json"
                )
                solver_command = [
                    sys.executable,
                    str(SOLVER),
                    "--export",
                    str(export_path),
                    "--threads",
                    str(args.threads),
                    "--jobs",
                    str(args.jobs),
                    "--timeout",
                    str(args.timeout),
                    "--summary",
                    str(summary_path),
                ]
                if not args.all_groups:
                    solver_command.extend(("--max-groups", str(args.max_groups)))
                subprocess.run(solver_command, cwd=ROOT, check=True)
                summary = json.loads(summary_path.read_text())
                job.update(
                    {
                        "solver_command": shlex.join(solver_command),
                        "summary": display_path(summary_path),
                        "summary_sha256": digest(summary_path),
                        "solver_status": summary["status"],
                        "attempted_distinct_systems": int(
                            summary["attempted_distinct_systems"]
                        ),
                        "classification_counts_by_distinct_system": summary[
                            "classification_counts_by_distinct_system"
                        ],
                    }
                )
            jobs.append(job)

    complete_modp = bool(jobs) and not args.export_only and args.all_groups and all(
        job.get("solver_status") == "PASS_COMPLETE_MODP_POLYNOMIAL_SECTION_SCHEME_SOLVED"
        for job in jobs
    )
    result = {
        "schema": SCHEMA,
        "status": (
            "PASS_COMPLETE_DECLARED_MODP_PO0_SCHEMES"
            if complete_modp
            else "INCOMPLETE_BOUNDED_V4_PO0_CAMPAIGN"
        ),
        "inputs": {
            display_path(path): digest(path)
            for path in (
                Path(__file__).resolve(),
                args.shortlist,
                *((args.rank_screen,) if rank_screen is not None else ()),
                args.bisections,
                args.direct,
                EXPORTER,
                SOLVER,
            )
        },
        "limits": {
            "pair_limit": args.pair_limit,
            "base_selection": args.selection,
            "primes_per_pair": args.primes_per_pair,
            "export_only": args.export_only,
            "all_distinct_system_groups": args.all_groups,
            "maximum_distinct_system_groups_per_job": (
                None if args.all_groups else args.max_groups
            ),
            "timeout_seconds_per_distinct_system": args.timeout,
            "threads_per_msolve_process": args.threads,
            "concurrent_msolve_processes": args.jobs,
        },
        "jobs": jobs,
        "reproducing_command": shlex.join(sys.argv),
        "proof_boundary": (
            "Each export is the exact chi=4 polynomial P.O=0 system in the displayed "
            "finite-field chart. Only jobs marked complete solve every distinct exported "
            "system. Timeouts, empty modular systems, and finite-field solutions are not "
            "characteristic-zero Mordell--Weil conclusions. A QQ section must be lifted "
            "and verified exactly before any rank-20 claim."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(
        "ALTV4PO0CAMPAIGN|"
        f"pairs={len(pairs)}|jobs={len(jobs)}|status={result['status']}|"
        f"output={display_path(args.output)}",
        flush=True,
    )


if __name__ == "__main__":
    main()

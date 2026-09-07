#!/usr/bin/env python3
"""Run one deterministic MW<=3 source scout in every eligible foundry class.

The child calculation is delegated to ``hunt_lattice_foundry_rootful_source.sage``.
This driver supplies coverage and provenance: it chooses the lowest-MW
catalogued target in the requested MW15--17 band, assigns a deterministic seed,
runs independent classes concurrently, and collects their exact child
certificates into one bounded-search ledger.

The current child Niemeier gluing certificate supports cyclic discriminant
groups.  Noncyclic classes are reported as unsupported, never as negative
search results.  A missed MW<=3 target is complete only for the declared beam,
samples, primes, seed, and generation bound.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DATABASE = ROOT / "artifacts/generated-results/elkies-k3-lattice-foundry-v1.json"
HUNTER = HERE / "hunt_lattice_foundry_rootful_source.sage"
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-lattice-foundry-mw3-broad-scout-v1.json"
)
DEFAULT_PRIMES = "3,7,11,13,17,23"
PINNED_SAGE_PYTHON = Path(
    "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python"
)
DEFAULT_SAGE_PYTHON = (
    PINNED_SAGE_PYTHON
    if PINNED_SAGE_PYTHON.is_file()
    else Path(sys.executable)
)


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def numeric_ns_id(ns_id: str) -> int:
    assert ns_id.startswith("NS") and ns_id[2:].isdigit()
    return int(ns_id[2:])


def artifact_paths(ns_id: str) -> tuple[Path, Path]:
    stem = f"elkies-k3-lattice-foundry-{ns_id.lower()}-mw3-broad-scout-v1"
    output = ROOT / "artifacts/generated-results" / f"{stem}.json"
    frame = ROOT / "artifacts/generated-results" / f"{stem}-root-adapted.txt"
    return output, frame


def selected_target(ns_row: dict, mw_min: int, mw_max: int) -> dict | None:
    eligible = [
        frame
        for frame in ns_row["frames"]
        if mw_min <= int(frame["mw_rank_for_rho_19"]) <= mw_max
    ]
    if not eligible:
        return None
    return min(
        eligible,
        key=lambda frame: (
            int(frame["mw_rank_for_rho_19"]),
            -int(frame["signed_root_count"]),
            frame["frame_id"],
        ),
    )


def child_command(
    ns_id: str,
    frame_id: str,
    output: Path,
    frame_output: Path,
    arguments: argparse.Namespace,
) -> list[str]:
    return [
        str(arguments.sage_python),
        str(HUNTER),
        "--database",
        str(arguments.database.resolve()),
        "--ns-id",
        ns_id,
        "--target-frame-id",
        frame_id,
        "--generations",
        str(arguments.generations),
        "--beam",
        str(arguments.beam),
        "--samples-per-parent",
        str(arguments.samples_per_parent),
        "--primes",
        arguments.primes,
        "--seed",
        str(arguments.seed_base + numeric_ns_id(ns_id)),
        "--target-root-rank",
        str(arguments.target_root_rank),
        "--allow-below-target",
        "--output",
        str(output),
        "--root-adapted-frame-output",
        str(frame_output),
    ]


def payload_matches_run(
    payload: dict,
    ns_id: str,
    frame_id: str,
    arguments: argparse.Namespace,
) -> bool:
    search = payload.get("search", {})
    return (
        payload.get("schema") == "elkies-k3.lattice-foundry-rootful-source.v1"
        and payload.get("target", {}).get("ns_id") == ns_id
        and payload.get("target", {}).get("frame_id") == frame_id
        and int(search.get("generations_bound", -1)) == arguments.generations
        and int(search.get("beam", -1)) == arguments.beam
        and int(search.get("samples_per_parent", -1)) == arguments.samples_per_parent
        and int(search.get("seed", -1))
        == arguments.seed_base + numeric_ns_id(ns_id)
        and int(search.get("target_root_rank", -1)) == arguments.target_root_rank
        and search.get("allow_below_target") is True
    )


def run_one(specification: dict, arguments: argparse.Namespace) -> dict:
    ns_id = specification["ns_id"]
    frame_id = specification["frame_id"]
    output = specification["output"]
    frame_output = specification["frame_output"]
    if output.is_file() and frame_output.is_file() and not arguments.rebuild:
        payload = json.loads(output.read_text())
        if payload_matches_run(payload, ns_id, frame_id, arguments):
            return {"status": "REUSED_MATCHING_CHILD", **specification}
    command = child_command(ns_id, frame_id, output, frame_output, arguments)
    completed = subprocess.run(
        command,
        cwd=ROOT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        timeout=arguments.timeout_seconds,
        check=False,
    )
    if completed.returncode != 0:
        return {
            "status": "FAILED_CHILD",
            "returncode": completed.returncode,
            "output_tail": completed.stdout[-4000:],
            **specification,
        }
    return {
        "status": "PASS_NEW_CHILD",
        "output_tail": completed.stdout[-1000:],
        **specification,
    }


def source_row(specification: dict) -> dict:
    output = specification["output"]
    frame_output = specification["frame_output"]
    payload = json.loads(output.read_text())
    source = payload["source"]
    search = payload["search"]
    mw_rank = int(source["mw_rank_for_rho_19"])
    components = source["root_components"]
    return {
        "ns_id": specification["ns_id"],
        "determinant": specification["determinant"],
        "target_frame_id": specification["frame_id"],
        "target_mw_rank": specification["target_mw_rank"],
        "child_certificate": relative(output),
        "child_certificate_sha256": digest(output),
        "root_adapted_frame": relative(frame_output),
        "root_adapted_frame_sha256": digest(frame_output),
        "source_root_type": source["root_type"],
        "source_root_rank": int(source["root_rank"]),
        "source_mw_rank": mw_rank,
        "mw_at_most_three_hit": mw_rank <= 3,
        "reducible_fibre_support_count": len(components),
        "semistable_configuration_compatible": all(
            component["type"].startswith("A") for component in components
        ),
        "target_root_rank_reached": bool(search["target_root_rank_reached"]),
        "generations_used": int(search["generations_used"]),
        "visited_reduced_keys": int(search["visited_reduced_keys"]),
        "kneser_discovery_edges": int(
            payload["kneser_p_neighbor_provenance"]["edge_count"]
        ),
        "niemeier_ambient": payload["niemeier_certificate"]["ambient_label"],
    }


def build_ledger(
    database_path: Path,
    arguments: argparse.Namespace,
    specifications: list[dict],
    unsupported: list[dict],
    failures: list[dict],
) -> dict:
    rows = [source_row(specification) for specification in specifications]
    rows.sort(
        key=lambda row: (
            0 if row["source_mw_rank"] <= 3 else 1,
            row["source_mw_rank"],
            row["reducible_fibre_support_count"],
            0 if row["semistable_configuration_compatible"] else 1,
            row["ns_id"],
        )
    )
    return {
        "schema": "elkies-k3.lattice-foundry-mw3-broad-scout.v1",
        "status": (
            "PASS_EXACT_CHILD_CERTIFICATES_BOUNDED_DISCOVERY_LEDGER"
            if not failures
            else "PARTIAL_CHILD_FAILURES_RECORDED"
        ),
        "objective": (
            "For each eligible foundry NS class, start at a catalogued MW15--17 "
            "frame and seek an exact source fibration of MW rank at most three."
        ),
        "proof_boundary": {
            "proved": (
                "Every retained child row has an exact primitive root lattice, MW "
                "rank, Niemeier gluing certificate, and deterministic run provenance."
            ),
            "not_proved": (
                "A miss is only a negative result inside its declared Kneser beam. "
                "Semistable compatibility does not construct an equation or prove "
                "rational descent. Unsupported noncyclic classes were not searched."
            ),
        },
        "search": {
            "target_mw_range": [arguments.target_mw_min, arguments.target_mw_max],
            "target_root_rank": arguments.target_root_rank,
            "generations": arguments.generations,
            "beam": arguments.beam,
            "samples_per_parent": arguments.samples_per_parent,
            "primes": [int(value) for value in arguments.primes.split(",")],
            "seed_base": arguments.seed_base,
            "selection_rule": (
                "minimum target MW rank, then maximum signed-root count, then frame id"
            ),
        },
        "coverage": {
            "database_ns_classes": len(
                json.loads(database_path.read_text())["ns_classes"]
            ),
            "searched_cyclic_classes": len(rows),
            "unsupported_noncyclic_classes": len(unsupported),
            "failed_children": len(failures),
            "mw_at_most_three_hits": sum(
                row["mw_at_most_three_hit"] for row in rows
            ),
        },
        "unsupported": unsupported,
        "failures": [
            {
                key: relative(value) if isinstance(value, Path) else value
                for key, value in failure.items()
            }
            for failure in failures
        ],
        "candidates": rows,
        "inputs": {relative(database_path): digest(database_path)},
        "reproduce": (
            f"{sys.executable} {relative(Path(__file__))} "
            f"--sage-python {arguments.sage_python} "
            f"--workers {arguments.workers} --generations {arguments.generations} "
            f"--beam {arguments.beam} --samples-per-parent "
            f"{arguments.samples_per_parent} --target-root-rank "
            f"{arguments.target_root_rank} --seed-base {arguments.seed_base}"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--database", type=Path, default=DATABASE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--ns-id", action="append", default=[])
    parser.add_argument("--target-mw-min", type=int, default=15)
    parser.add_argument("--target-mw-max", type=int, default=17)
    parser.add_argument("--target-root-rank", type=int, default=14)
    parser.add_argument("--generations", type=int, default=10)
    parser.add_argument("--beam", type=int, default=10)
    parser.add_argument("--samples-per-parent", type=int, default=50)
    parser.add_argument("--primes", default=DEFAULT_PRIMES)
    parser.add_argument("--sage-python", type=Path, default=DEFAULT_SAGE_PYTHON)
    parser.add_argument("--seed-base", type=int, default=20263100)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--timeout-seconds", type=int, default=900)
    parser.add_argument("--rebuild", action="store_true")
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()
    if arguments.workers <= 0 or arguments.generations <= 0:
        parser.error("workers and generations must be positive")
    if arguments.beam <= 0 or arguments.samples_per_parent <= 0:
        parser.error("beam and samples per parent must be positive")
    if not arguments.sage_python.is_file():
        parser.error(f"Sage Python executable not found: {arguments.sage_python}")

    database_path = arguments.database.resolve()
    database = json.loads(database_path.read_text())
    wanted = set(arguments.ns_id)
    known = {row["ns_id"] for row in database["ns_classes"]}
    if wanted - known:
        parser.error(f"unknown NS ids: {sorted(wanted - known)}")

    specifications = []
    unsupported = []
    for ns_row in database["ns_classes"]:
        ns_id = ns_row["ns_id"]
        if wanted and ns_id not in wanted:
            continue
        target = selected_target(
            ns_row, arguments.target_mw_min, arguments.target_mw_max
        )
        if target is None:
            unsupported.append({"ns_id": ns_id, "reason": "NO_TARGET_IN_MW_BAND"})
            continue
        if int(ns_row["discriminant_length"]) != 1:
            unsupported.append(
                {
                    "ns_id": ns_id,
                    "determinant": int(ns_row["determinant"]),
                    "discriminant_length": int(ns_row["discriminant_length"]),
                    "reason": "NONCYCLIC_DISCRIMINANT_NOT_SUPPORTED_BY_CHILD_GLUE_V1",
                }
            )
            continue
        output, frame_output = artifact_paths(ns_id)
        specifications.append(
            {
                "ns_id": ns_id,
                "determinant": int(ns_row["determinant"]),
                "frame_id": target["frame_id"],
                "target_mw_rank": int(target["mw_rank_for_rho_19"]),
                "output": output,
                "frame_output": frame_output,
            }
        )

    if arguments.check:
        missing = [
            specification["ns_id"]
            for specification in specifications
            if not specification["output"].is_file()
            or not specification["frame_output"].is_file()
        ]
        if missing:
            raise SystemExit(f"missing child artifacts for: {missing}")
        mismatched = []
        for specification in specifications:
            payload = json.loads(specification["output"].read_text())
            if not payload_matches_run(
                payload,
                specification["ns_id"],
                specification["frame_id"],
                arguments,
            ):
                mismatched.append(specification["ns_id"])
        if mismatched:
            raise SystemExit(f"child artifacts do not match run: {mismatched}")
        failures = []
        successful = specifications
    else:
        results = []
        with ThreadPoolExecutor(max_workers=arguments.workers) as executor:
            futures = {
                executor.submit(run_one, specification, arguments): specification
                for specification in specifications
            }
            for future in as_completed(futures):
                result = future.result()
                results.append(result)
                print(
                    "FOUNDRYMW3BROAD|"
                    f"ns={result['ns_id']}|status={result['status']}",
                    flush=True,
                )
        failures = [row for row in results if row["status"] == "FAILED_CHILD"]
        successful_ids = {
            row["ns_id"] for row in results if row["status"] != "FAILED_CHILD"
        }
        successful = [
            row for row in specifications if row["ns_id"] in successful_ids
        ]

    ledger = build_ledger(
        database_path, arguments, successful, unsupported, failures
    )
    serialized = json.dumps(ledger, indent=2, sort_keys=True) + "\n"
    output_path = arguments.output.resolve()
    if arguments.check:
        if output_path.read_text() != serialized:
            raise SystemExit("MW3 broad-scout ledger is stale")
    else:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(serialized)
    print(
        "FOUNDRYMW3BROAD|"
        f"searched={ledger['coverage']['searched_cyclic_classes']}|"
        f"unsupported={ledger['coverage']['unsupported_noncyclic_classes']}|"
        f"hits={ledger['coverage']['mw_at_most_three_hits']}|"
        f"failures={ledger['coverage']['failed_children']}|status=PASS",
        flush=True,
    )


if __name__ == "__main__":
    main()

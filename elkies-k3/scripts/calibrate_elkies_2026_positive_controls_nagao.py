#!/usr/bin/env python3
"""Calibrate a complete compact-``t`` Nagao scan on the four published fibres.

The prime set is split round-robin into three pairwise-disjoint ensembles.
Each local contribution is centered and population-standardized over the good
fibres of P^1(F_p); a parameter's primary score is its weakest normalized
ensemble sum.  The C++ hot loop enumerates every primitive ``t=a/b`` in the
requested box and reports exact population ranks for the four controls.

The method is accepted only if every control lies in the top one percent of
the complete H <= 10,000 population.  This is a heuristic acceptance test,
not a rank bound and not a substitute for residual descent.
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
from time import perf_counter


SCRIPTS = Path(__file__).resolve().parent
ROOT = SCRIPTS.parents[1]
sys.path.insert(0, str(SCRIPTS))

from search_h92_q12o5867_rootless_nagao import (  # noqa: E402
    build_residue_tables,
    export_cpp_tables,
    is_prime,
    load_family_model,
)


CPP = SCRIPTS / "scan_h92_q12o5867_rootless_nagao.cpp"
CONTROLS = ((-2, 377), (-308, 251), (2456, 135), (-9529, 5471))
CONTROL_CERTIFICATE = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "elkies_2026_high_rank_positive_controls_v2.json"
)
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "elkies_2026_compact_t_nagao_positive_control_h10000_v1.json"
)


def file_sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def round_robin_prime_ensembles(
    lower: int, upper: int, count: int = 3
) -> tuple[tuple[int, ...], ...]:
    primes = [value for value in range(lower, upper + 1) if is_prime(value)]
    if count < 3 or len(primes) < count:
        raise ValueError("at least three nonempty prime ensembles are required")
    return tuple(tuple(primes[index::count]) for index in range(count))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--height", type=int, default=10_000)
    parser.add_argument("--prime-lower", type=int, default=19)
    parser.add_argument("--prime-upper", type=int, default=599)
    parser.add_argument("--ensemble-count", type=int, default=3)
    parser.add_argument("--finalists", type=int, default=1000)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--local-directory",
        type=Path,
        default=ROOT / "artifacts/local/elkies-k3/elkies-2026-control-nagao",
    )
    args = parser.parse_args()
    if args.height < 1 or args.finalists < 1:
        raise SystemExit("--height and --finalists must be positive")

    control_document = json.loads(CONTROL_CERTIFICATE.read_text())
    if control_document.get("status") != (
        "PASS_EXACT_ELKIES_2026_HIGH_RANK_POSITIVE_CONTROLS"
    ):
        raise SystemExit("the exact four-fibre positive-control certificate is missing")
    observed_gains = [
        row["public_positive_control"]["quotient_gain_beyond_generic_rank_17"]
        for row in control_document["fibres"]
    ]
    if observed_gains != [8, 9, 10, 11]:
        raise SystemExit("the exact positive-control quotient gains changed")

    compiler = shutil.which("g++")
    if compiler is None:
        raise SystemExit("g++ is required for the complete H=10,000 scan")
    ensembles = round_robin_prime_ensembles(
        args.prime_lower, args.prime_upper, args.ensemble_count
    )
    model = load_family_model()
    started = perf_counter()
    table_blocks, rejected = build_residue_tables(model, ensembles)
    if rejected:
        raise SystemExit(f"unexpected rejected calibration primes: {rejected}")

    args.local_directory.mkdir(parents=True, exist_ok=True)
    table_path = args.local_directory / "compact-t-tables.txt"
    binary_path = args.local_directory / "scan-worst-block"
    export_cpp_tables(table_path, model, table_blocks)
    compile_command = [
        compiler,
        "-O3",
        "-std=c++17",
        "-Wall",
        "-Wextra",
        "-pedantic",
        str(CPP),
        "-o",
        str(binary_path),
    ]
    subprocess.run(compile_command, check=True)
    controls_text = ",".join(f"{a}/{b}" for a, b in CONTROLS)
    scan_command = [
        str(binary_path),
        str(table_path),
        str(args.height),
        str(args.height),
        "1",
        ",".join("1" for _ in ensembles),
        str(args.finalists),
        str(args.output),
        "1",
        controls_text,
    ]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(scan_command, check=True)

    result = json.loads(args.output.read_text())
    if result.get("status") != "PASS_POSITIVE_CONTROL_SCORING_GATE":
        raise SystemExit("the scoring method failed its positive-control gate")
    result["orchestration"] = {
        "script": str(Path(__file__).resolve()),
        "script_sha256": file_sha256(Path(__file__).resolve()),
        "scanner_source": str(CPP.resolve()),
        "scanner_source_sha256": file_sha256(CPP),
        "positive_control_certificate": str(CONTROL_CERTIFICATE.resolve()),
        "positive_control_certificate_sha256": file_sha256(CONTROL_CERTIFICATE),
        "table_sha256": file_sha256(table_path),
        "compile_command": shlex.join(compile_command),
        "scan_command": shlex.join(scan_command),
        "wall_seconds_including_tables_and_compile": perf_counter() - started,
    }
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    ranks = [row["population_rank"] for row in result["positive_controls"]]
    print(
        f"ELKIES2026NAGAO|height={args.height}|population={result['population_count']}|"
        f"control_ranks={','.join(map(str, ranks))}|status={result['status']}|"
        f"output={args.output}"
    )


if __name__ == "__main__":
    main()

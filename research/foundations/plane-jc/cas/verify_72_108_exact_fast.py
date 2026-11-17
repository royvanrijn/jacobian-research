#!/usr/bin/env python3
"""Portable parallel replay of the pinned (72,108) exact certificates."""

from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import os
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_REPLAY = (
    REPO_ROOT
    / "plane-jc"
    / "external"
    / "zenodo-21479814"
    / "bilLkarkariy-jc2-72-108-exact-certificates-d9ea4fd"
    / "release_bundle"
    / "exact_replay"
)
FAST_HARD = Path(__file__).with_name("verify_72_108_hard_fast.py")


@dataclass(frozen=True)
class Stage:
    label: str
    scripts: tuple[Path, ...]


@dataclass(frozen=True)
class StageResult:
    label: str
    seconds: float
    stdout: str
    stderr: str
    returncode: int


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(8 * 1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def verify_manifest(replay: Path, name: str) -> None:
    manifest = replay / name
    for line_number, raw_line in enumerate(
        manifest.read_text(encoding="utf-8").splitlines(), 1
    ):
        if not raw_line:
            continue
        try:
            expected, relative = raw_line.split(maxsplit=1)
        except ValueError as error:
            raise RuntimeError(
                f"malformed {name} line {line_number}"
            ) from error
        relative = relative.removeprefix("*")
        target = replay / relative
        actual = sha256(target)
        if actual != expected:
            raise RuntimeError(
                f"{name}: checksum mismatch for {relative}: "
                f"expected {expected}, got {actual}"
            )
        print(f"{relative}: OK")


def run_stage(stage: Stage, replay: Path, environment: dict[str, str]) -> StageResult:
    started = time.perf_counter()
    stdout_parts = []
    stderr_parts = []
    returncode = 0
    for script in stage.scripts:
        completed = subprocess.run(
            [sys.executable, str(script)],
            cwd=replay,
            env=environment,
            text=True,
            capture_output=True,
            check=False,
        )
        stdout_parts.append(completed.stdout)
        stderr_parts.append(completed.stderr)
        if completed.returncode:
            returncode = completed.returncode
            break
    return StageResult(
        label=stage.label,
        seconds=time.perf_counter() - started,
        stdout="".join(stdout_parts),
        stderr="".join(stderr_parts),
        returncode=returncode,
    )


def run_wave(
    stages: tuple[Stage, ...],
    replay: Path,
    environment: dict[str, str],
    jobs: int,
) -> None:
    with concurrent.futures.ThreadPoolExecutor(
        max_workers=min(jobs, len(stages))
    ) as executor:
        futures = {
            stage.label: executor.submit(run_stage, stage, replay, environment)
            for stage in stages
        }
        results = {label: future.result() for label, future in futures.items()}

    failures = []
    for stage in stages:
        result = results[stage.label]
        if result.stdout:
            print(result.stdout, end="")
        if result.stderr:
            print(result.stderr, end="", file=sys.stderr)
        print(f"TIMING {result.label} {result.seconds:.3f}s")
        if result.returncode:
            failures.append(f"{result.label} (exit {result.returncode})")
    if failures:
        raise RuntimeError("failed replay stages: " + ", ".join(failures))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--replay-root",
        type=Path,
        default=DEFAULT_REPLAY,
        help="path to the extracted exact_replay directory",
    )
    parser.add_argument(
        "--jobs",
        type=int,
        default=min(4, os.cpu_count() or 1),
        help="maximum concurrent CPU stages (default: min(4, CPU count))",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.jobs < 1:
        raise SystemExit("--jobs must be positive")
    replay = args.replay_root.resolve()
    if not (replay / "hard" / "h_certificate_exact.txt").is_file():
        raise SystemExit(f"missing hard certificate under {replay}")

    started = time.perf_counter()
    verify_manifest(replay, "EXACT_SHA256SUMS.txt")
    verify_manifest(replay, "RECONSTRUCTED_CERTIFICATES.sha256")

    environment = os.environ.copy()
    environment["PYTHONHASHSEED"] = "0"
    environment["JC2_72_108_REPLAY_ROOT"] = str(replay)

    run_wave(
        (
            Stage("laurent", (replay / "verify_laurent_reduction.py",)),
            Stage("case1-reduction", (replay / "verify_case1_reduction.py",)),
        ),
        replay,
        environment,
        args.jobs,
    )
    run_wave(
        (
            Stage("firstblock", (replay / "verify_firstblock_exact.py",)),
            Stage("case1-cascade", (replay / "case1_cascade_machine.py",)),
            Stage(
                "branch-one",
                (replay / "derive_hne0.py", replay / "build_degree5.py"),
            ),
            Stage(
                "branch-two",
                (
                    replay / "derive_hne0_branch2.py",
                    replay / "build_degree5_branch2.py",
                ),
            ),
        ),
        replay,
        environment,
        args.jobs,
    )
    run_wave(
        (
            Stage(
                "serialized",
                (replay / "verify_serialized_certificates.py",),
            ),
            Stage("hard-and-symmetry", (FAST_HARD,)),
        ),
        replay,
        environment,
        args.jobs,
    )

    print(f"TOTAL_SECONDS={time.perf_counter() - started:.3f}")
    print("JC2_72_108_EXACT_REPLAY_PASS")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Run a bounded PARI ellrank diagnostic on one conductor-first target.

PARI's cubic-field BNF is provisional under GRH unless separately certified.
This runner therefore never promotes the returned upper endpoint to an
unconditional rank theorem.  Rational points returned by PARI are checked on
the exact descent model and greedily admitted only when an independent exact
mod-2 finite-reduction certificate proves a genuine rank gain.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import shutil
import signal
import subprocess
import sys
import tempfile
import time
from typing import Any, Sequence


ROOT = Path(__file__).resolve().parents[2]
CAS = Path(__file__).resolve().parent
sys.path.insert(0, str(CAS))

from build_conductor_first_near_miss_magma import (  # noqa: E402
    DEFAULT_MANIFEST,
    load_target,
)
from build_conductor_first_near_miss_targets import mod2_certificate  # noqa: E402
from elliptic_candidate_record import is_on_weierstrass_curve  # noqa: E402


Q = Fraction
PROTOCOL = "CFNMPARI"
SCHEMA = "elliptic-curves.conductor-first-pari-diagnostic.v1"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def qtext(value: Fraction | int | str) -> str:
    value = Q(value)
    return str(value.numerator) if value.denominator == 1 else str(value)


def load_exact_input(
    manifest_path: Path, target_id: str
) -> tuple[dict[str, Any], dict[str, Any], tuple[Q, ...], tuple[tuple[Q, Q], ...]]:
    manifest, target = load_target(manifest_path, target_id)
    model = tuple(Q(value) for value in target["descent_model"])
    points = tuple(tuple(Q(value) for value in point) for point in target["known_basis"])
    rank = int(target["certified_known_rank"])
    if len(model) != 5:
        raise ValueError("descent_model must contain five Weierstrass coefficients")
    if len(points) != rank or len(set(points)) != rank:
        raise ValueError("known_basis does not contain the certified number of distinct points")
    if any(not is_on_weierstrass_curve(model, point) for point in points):
        raise ArithmeticError("known_basis contains a point off descent_model")
    return manifest, target, model, points


def gp_program(
    target_id: str,
    model: Sequence[Q],
    points: Sequence[tuple[Q, Q]],
    *,
    stack_bytes: int,
    effort: int,
) -> str:
    if stack_bytes < 64_000_000 or effort < 0:
        raise ValueError("invalid PARI stack or effort")
    model_text = ",".join(qtext(value) for value in model)
    point_text = ",".join(
        "[" + ",".join(qtext(value) for value in point) + "]" for point in points
    )
    return f'''default(parisizemax,{stack_bytes});
E=ellinit([{model_text}]);
P=[{point_text}];
print("{PROTOCOL}|version=1|target={target_id}|stage=input|status=complete|pari=",version(),"|known=",#P);
print("{PROTOCOL}|target={target_id}|stage=ellrank|status=start|effort={effort}");
gettime();
iferr(R=ellrank(E,{effort},P),ERR,print("{PROTOCOL}|target={target_id}|stage=ellrank|status=error|pari_error=",ERR);quit(2));
print("{PROTOCOL}|target={target_id}|stage=ellrank|status=complete|milliseconds=",gettime(),"|lower=",R[1],"|upper=",R[2],"|sha=",R[3],"|point_count=",#R[4]);
for(i=1,#R[4],print("{PROTOCOL}|target={target_id}|stage=point|index=",i,"|point=",R[4][i]));
print("{PROTOCOL}|target={target_id}|stage=all|status=complete");
quit;
'''


def marker_fields(line: str) -> dict[str, str] | None:
    if not line.startswith(f"{PROTOCOL}|"):
        return None
    fields: dict[str, str] = {}
    for item in line.rstrip("\n").split("|")[1:]:
        if "=" in item:
            key, value = item.split("=", 1)
            fields[key] = value
    return fields


POINT_RE = re.compile(r"^\[([^,]+),\s*([^,]+)\]$")


def parse_point(text: str) -> tuple[Q, Q]:
    match = POINT_RE.fullmatch(text.strip())
    if match is None:
        raise ValueError(f"cannot parse PARI point {text!r}")
    return Q(match.group(1)), Q(match.group(2))


def parse_output(stdout: str) -> dict[str, Any]:
    result: dict[str, Any] = {
        "pari_version": None,
        "completed": False,
        "rank_lower": None,
        "rank_upper": None,
        "sha_indicator": None,
        "elapsed_milliseconds": None,
        "reported_point_count": None,
        "points": [],
        "last_stage": None,
    }
    for line in stdout.splitlines():
        fields = marker_fields(line)
        if fields is None:
            continue
        result["last_stage"] = fields.get("stage", result["last_stage"])
        if "pari" in fields:
            result["pari_version"] = fields["pari"]
        if fields.get("stage") == "ellrank" and fields.get("status") == "complete":
            result["rank_lower"] = int(fields["lower"])
            result["rank_upper"] = int(fields["upper"])
            result["sha_indicator"] = int(fields["sha"])
            result["elapsed_milliseconds"] = int(fields["milliseconds"])
            result["reported_point_count"] = int(fields["point_count"])
        elif fields.get("stage") == "point":
            result["points"].append(parse_point(fields["point"]))
        elif fields.get("stage") == "all" and fields.get("status") == "complete":
            result["completed"] = True
    count = result["reported_point_count"]
    if count is not None and count != len(result["points"]):
        raise RuntimeError("PARI point marker count is inconsistent")
    return result


def read_rss_bytes(pid: int) -> int:
    for line in Path(f"/proc/{pid}/status").read_text().splitlines():
        if line.startswith("VmRSS:"):
            return int(line.split()[1]) * 1024
    raise RuntimeError(f"VmRSS is absent for pid {pid}")


def terminate_owned(process: subprocess.Popen[str], sig: signal.Signals) -> None:
    if process.poll() is None:
        group = os.getpgid(process.pid)
        if group != process.pid:
            raise RuntimeError("refusing to signal an unexpected process group")
        os.killpg(group, sig)


def run_process(
    executable: str,
    program: str,
    *,
    stack_bytes: int,
    timeout_seconds: float,
    rss_limit_bytes: int,
) -> dict[str, Any]:
    with tempfile.TemporaryFile(mode="w+") as stdout_file, tempfile.TemporaryFile(
        mode="w+"
    ) as stderr_file:
        process = subprocess.Popen(
            [executable, "-f", "-q", "-s", str(stack_bytes)],
            stdin=subprocess.PIPE,
            stdout=stdout_file,
            stderr=stderr_file,
            text=True,
            start_new_session=True,
        )
        assert process.stdin is not None
        process.stdin.write(program)
        process.stdin.close()
        started = time.monotonic()
        peak_rss = 0
        outcome = "running"
        while process.poll() is None:
            elapsed = time.monotonic() - started
            if elapsed >= timeout_seconds:
                outcome = "strict_wall_timeout"
                terminate_owned(process, signal.SIGTERM)
                break
            try:
                rss = read_rss_bytes(process.pid)
            except (FileNotFoundError, ProcessLookupError):
                break
            peak_rss = max(peak_rss, rss)
            if rss > rss_limit_bytes:
                outcome = "strict_rss_limit"
                terminate_owned(process, signal.SIGTERM)
                break
            time.sleep(0.25)
        if process.poll() is None:
            try:
                process.wait(timeout=10)
            except subprocess.TimeoutExpired:
                terminate_owned(process, signal.SIGKILL)
                process.wait()
        wall_seconds = time.monotonic() - started
        stdout_file.seek(0)
        stderr_file.seek(0)
        stdout = stdout_file.read()
        stderr = stderr_file.read()
    parsed = parse_output(stdout)
    if outcome == "running":
        outcome = (
            "completed"
            if process.returncode == 0 and parsed["completed"]
            else "pari_failure"
        )
    return {
        "outcome": outcome,
        "returncode": process.returncode,
        "wall_seconds": wall_seconds,
        "peak_observed_rss_bytes": peak_rss,
        "stdout": stdout,
        "stderr": stderr,
        "parsed": parsed,
    }


def negative(model: Sequence[Q], point: tuple[Q, Q]) -> tuple[Q, Q]:
    a1, _, a3, _, _ = model
    x, y = point
    return x, -a1 * x - a3 - y


def exact_point_gains(
    model: Sequence[Q],
    known: Sequence[tuple[Q, Q]],
    returned: Sequence[tuple[Q, Q]],
) -> dict[str, Any]:
    if any(not is_on_weierstrass_curve(model, point) for point in returned):
        raise ArithmeticError("PARI returned a point off the descent model")
    accepted = list(known)
    new_points: list[tuple[Q, Q]] = []
    final_certificate = None
    for point in returned:
        if point in accepted or negative(model, point) in accepted:
            continue
        try:
            certificate = mod2_certificate(model, accepted + [point])
        except ArithmeticError:
            continue
        accepted.append(point)
        new_points.append(point)
        final_certificate = certificate
    return {
        "exact_on_curve": True,
        "certified_new_points": [[qtext(x), qtext(y)] for x, y in new_points],
        "certified_rank_lower_bound": len(accepted),
        "combined_mod2_certificate": final_certificate,
    }


def run(args: argparse.Namespace) -> int:
    _, target, model, known = load_exact_input(args.manifest, args.target)
    executable = str(args.gp) if args.gp is not None else shutil.which("gp")
    if executable is None or not Path(executable).is_file():
        raise FileNotFoundError("PARI/GP executable was not found")
    program = gp_program(
        target["id"], model, known, stack_bytes=args.stack_bytes, effort=args.effort
    )
    print(
        f"{PROTOCOL}|target={target['id']}|stage=supervisor|status=start"
        f"|timeout={args.timeout}|rss_limit={args.rss_limit_bytes}",
        flush=True,
    )
    process = run_process(
        executable,
        program,
        stack_bytes=args.stack_bytes,
        timeout_seconds=args.timeout,
        rss_limit_bytes=args.rss_limit_bytes,
    )
    parsed = process.pop("parsed")
    gains = exact_point_gains(model, known, parsed["points"])
    status = (
        "COMPLETED_GRH_CONDITIONAL_UPPER_BOUND"
        if process["outcome"] == "completed"
        else "NO_COMPLETE_SELMER_DIAGNOSTIC"
    )
    payload = {
        "schema": SCHEMA,
        "status": status,
        "mathematical_status": "exact_point_lower_bound_only; upper endpoint GRH-conditional",
        "input": {
            "manifest": str(args.manifest.resolve()),
            "manifest_sha256": sha256_file(args.manifest),
            "target": target["id"],
            "descent_model": [qtext(value) for value in model],
            "known_rank_lower_bound": len(known),
            "known_point_membership_checked": True,
        },
        "backend": {
            "engine": "PARI/GP ellrank",
            "executable": str(Path(executable).resolve()),
            "pari_version": parsed["pari_version"],
            "python_version": platform.python_version(),
            "program_sha256": hashlib.sha256(program.encode()).hexdigest(),
            "script_sha256": sha256_file(Path(__file__).resolve()),
        },
        "bounds": {
            "effort": args.effort,
            "wall_timeout_seconds": args.timeout,
            "stack_bytes": args.stack_bytes,
            "rss_limit_bytes": args.rss_limit_bytes,
        },
        "process": process,
        "provisional_ellrank": {
            key: parsed[key]
            for key in (
                "completed",
                "rank_lower",
                "rank_upper",
                "sha_indicator",
                "elapsed_milliseconds",
                "reported_point_count",
                "last_stage",
            )
        },
        "exact_point_recovery": gains,
        "claim_boundary": [
            "PARI ellrank uses provisional cubic-field BNF data; no upper endpoint is unconditional without a separate bnfcertify result.",
            "A timeout, resource stop, or PARI failure is not a Selmer or rank bound.",
            "Only returned rational points with an independent full mod-2 certificate increase the unconditional rank lower bound.",
        ],
        "reproducing_command": (
            f".venv/bin/python {Path(__file__).relative_to(ROOT)} "
            f"--manifest {args.manifest} --target {target['id']} --output {args.output} "
            f"--gp {executable} --effort {args.effort} --timeout {args.timeout} "
            f"--stack-bytes {args.stack_bytes} --rss-limit-bytes {args.rss_limit_bytes}"
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    mode = "w" if args.overwrite else "x"
    with args.output.open(mode) as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(f"status={status}")
    print(f"outcome={process['outcome']}")
    print(f"certified_rank_lower_bound={gains['certified_rank_lower_bound']}")
    print(f"output={args.output.resolve()}")
    return 0 if process["outcome"] in {
        "completed",
        "strict_wall_timeout",
        "strict_rss_limit",
    } else 1


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--target", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--gp", type=Path)
    parser.add_argument("--effort", type=int, default=0)
    parser.add_argument("--timeout", type=float, default=3600.0)
    parser.add_argument("--stack-bytes", type=int, default=8_000_000_000)
    parser.add_argument("--rss-limit-bytes", type=int, default=10_000_000_000)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    if args.effort < 0:
        parser.error("--effort must be nonnegative")
    if args.timeout <= 0 or args.timeout > 86400:
        parser.error("--timeout must lie in (0,86400]")
    if min(args.stack_bytes, args.rss_limit_bytes) < 64_000_000:
        parser.error("stack and RSS limits must be at least 64 MB")
    raise SystemExit(run(args))


if __name__ == "__main__":
    main()

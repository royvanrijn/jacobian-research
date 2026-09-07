#!/usr/bin/env python3
"""Run PARI's unconditional u=42 2-descent with strict resource controls.

PARI documents the rank interval returned by ``ellrank`` as unconditional.
The expensive step here is initialization of the irreducible cubic field.
This runner keeps the process in its own recorded process group, enforces both
a wall-clock cap and an observed-RSS cap, and supplies the exact 17-point basis.
It writes no files; stdout is suitable for capture in a generated artifact.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time


REPOSITORY = Path(__file__).resolve().parents[2]
SOURCE = (
    REPOSITORY
    / "archive/elliptic-curves/artifacts/generated-results/elliptic_nagao_u42_height_10000000.json"
)
SOURCE_SHA256 = "4fea0207fd637988bcc1147143657cbec5c2404cb81b4c4a487e2dde20cc43b8"
A = Fraction(-74879150695093957092648257365083, 92236816)
B = Fraction(
    121839825430716337244033564674334552301153773691,
    442921190432,
)


def gp_rational(value: Fraction | str) -> str:
    rational = Fraction(value)
    if rational.denominator == 1:
        return str(rational.numerator)
    return f"({rational.numerator}/{rational.denominator})"


def load_points() -> list[tuple[str, str]]:
    raw = SOURCE.read_bytes()
    digest = hashlib.sha256(raw).hexdigest()
    if digest != SOURCE_SHA256:
        raise SystemExit(f"refusing changed input: expected {SOURCE_SHA256}, got {digest}")
    data = json.loads(raw)
    points = [
        (record["jacobian_x"], record["jacobian_y"])
        for record in data["small_prime_saturation"]["saturated_basis"]
    ]
    if len(points) != 17:
        raise SystemExit(f"expected 17 points, got {len(points)}")
    return points


def gp_program(stack_bytes: int) -> str:
    points = ",".join(
        f"[{gp_rational(x_value)},{gp_rational(y_value)}]"
        for x_value, y_value in load_points()
    )
    return f"""default(parisizemax,{stack_bytes});
default(realprecision,100);
E=ellinit([0,0,0,{gp_rational(A)},{gp_rational(B)}]);
P=[{points}];
print("PARI_VERSION=",version());
print("SOURCE_SHA256={SOURCE_SHA256}");
print("EXACT_POINTS_ON_CURVE=",vecsum(vector(#P,i,ellisoncurve(E,P[i]))));
gettime();
R=ellrank(E,0,P);
print("ELLRANK_DONE_MS=",gettime());
print("RANK_LOWER=",R[1]);
print("RANK_UPPER=",R[2]);
print("SHA_PAIRING_RANK=",R[3]);
print("RETURNED_POINTS=",#R[4]);
"""


def observed_rss_bytes(pid: int) -> int:
    result = subprocess.run(
        ["/bin/ps", "-p", str(pid), "-o", "rss="],
        check=True,
        capture_output=True,
        text=True,
    )
    value = result.stdout.strip()
    if not value:
        raise RuntimeError(f"ps returned no RSS for owned pid {pid}")
    return int(value) * 1024


def terminate_owned(process: subprocess.Popen[str], sig: signal.Signals) -> None:
    if process.poll() is not None:
        return
    process_group = os.getpgid(process.pid)
    if process_group != process.pid:
        raise RuntimeError(
            f"refusing unexpected process group {process_group} for pid {process.pid}"
        )
    os.killpg(process_group, sig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--timeout", type=float, default=600.0)
    parser.add_argument("--stack-bytes", type=int, default=8_000_000_000)
    parser.add_argument("--rss-limit-bytes", type=int, default=8_000_000_000)
    parser.add_argument("--heartbeat", type=float, default=45.0)
    args = parser.parse_args()
    if args.timeout <= 0 or args.timeout > 3600:
        raise SystemExit("--timeout must be in (0,3600]")
    if args.heartbeat <= 0 or args.heartbeat > 60:
        raise SystemExit("--heartbeat must be in (0,60]")
    if min(args.stack_bytes, args.rss_limit_bytes) < 64_000_000:
        raise SystemExit("stack and RSS limits must each be at least 64MB")

    process = subprocess.Popen(
        ["gp", "-fq"],
        stdin=subprocess.PIPE,
        text=True,
        start_new_session=True,
    )
    print(f"OWNED_GP_PID={process.pid}", file=sys.stderr, flush=True)
    assert process.stdin is not None
    process.stdin.write(gp_program(args.stack_bytes))
    process.stdin.close()
    started = time.monotonic()
    next_heartbeat = args.heartbeat

    try:
        while process.poll() is None:
            elapsed = time.monotonic() - started
            if elapsed >= args.timeout:
                terminate_owned(process, signal.SIGTERM)
                try:
                    process.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    terminate_owned(process, signal.SIGKILL)
                    process.wait()
                print(
                    f"STRICT_TIMEOUT_SECONDS={args.timeout:g}",
                    file=sys.stderr,
                    flush=True,
                )
                raise SystemExit(124)
            if elapsed >= next_heartbeat:
                rss = observed_rss_bytes(process.pid)
                print(
                    f"HEARTBEAT_SECONDS={int(elapsed)} RSS_BYTES={rss}",
                    file=sys.stderr,
                    flush=True,
                )
                if rss > args.rss_limit_bytes:
                    terminate_owned(process, signal.SIGTERM)
                    process.wait(timeout=5)
                    print(
                        f"STRICT_RSS_LIMIT_BYTES={args.rss_limit_bytes}",
                        file=sys.stderr,
                        flush=True,
                    )
                    raise SystemExit(125)
                next_heartbeat += args.heartbeat
            time.sleep(1)
    except BaseException:
        if process.poll() is None:
            terminate_owned(process, signal.SIGTERM)
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                terminate_owned(process, signal.SIGKILL)
                process.wait()
        raise
    raise SystemExit(process.returncode)


if __name__ == "__main__":
    main()

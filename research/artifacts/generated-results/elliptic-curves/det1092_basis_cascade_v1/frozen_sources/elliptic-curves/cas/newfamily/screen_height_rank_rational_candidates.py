#!/usr/bin/env sage -python
"""Fast height-matrix screen for rational newfamily specializations.

This is a triage stage, not an exact Mordell--Weil rank proof.  The eleven
hidden generic sections are known to be generically independent.  For each
rational T=a/b this script specializes them on the fixed short finite-minimal
family and computes their canonical-height Gram matrix.  A positive-definite
11x11 matrix with stable numerical rank 11 is enough to keep a specialization
for later extra-point searches; exact eclib work is deferred until a genuinely
new point is found.

Every candidate runs in an isolated Sage subprocess so pathological height
calls cannot stall the whole shortlist.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction
import json
from pathlib import Path
import subprocess
import sys
import time

from sage.all import EllipticCurve, QQ, RR, RealField, ZZ, load, matrix

ROOTS = (-47, -43, -31, 30, 45, 46)


@dataclass(frozen=True)
class Candidate:
    numerator: int
    denominator: int
    discovery: float
    held: float

    @property
    def total(self) -> float:
        return self.discovery + self.held

    @property
    def minimum(self) -> float:
        return min(self.discovery, self.held)

    @property
    def balanced(self) -> float:
        return 0.35 * self.discovery + 0.65 * self.held

    @property
    def name(self) -> str:
        return f"{self.numerator}/{self.denominator}"


def parse_scan(paths: list[Path]) -> list[Candidate]:
    unique: dict[tuple[int, int], Candidate] = {}
    for path in paths:
        for line in path.read_text().splitlines():
            if not line.startswith("C "):
                continue
            fields = line.split()
            c = Candidate(
                int(fields[1]), int(fields[2]), float(fields[3]), float(fields[4])
            )
            unique[(c.numerator, c.denominator)] = c
    return list(unique.values())


def diversify(candidates: list[Candidate], quota: int, small_denominator: int) -> list[Candidate]:
    selected: dict[tuple[int, int], Candidate] = {}

    def add(rows):
        for c in rows[:quota]:
            selected[(c.numerator, c.denominator)] = c

    add(sorted(candidates, key=lambda c: c.total, reverse=True))
    add(sorted(candidates, key=lambda c: c.held, reverse=True))
    add(sorted(candidates, key=lambda c: c.minimum, reverse=True))
    add(sorted(candidates, key=lambda c: c.balanced, reverse=True))
    add(
        sorted(
            (c for c in candidates if c.denominator <= small_denominator),
            key=lambda c: (c.total, c.held),
            reverse=True,
        )
    )
    return sorted(
        selected.values(),
        key=lambda c: (c.minimum, c.total, -c.denominator),
        reverse=True,
    )


def load_builder():
    paths = [
        Path("elliptic-curves/cas/newfamily/newfamily_rank11_minimal_common.py"),
        Path("elliptic-curves/cas/newfamily/archive/newfamily_rank11_minimal_common.py"),
        Path("/tmp/newfamily_rank11_minimal_common.py"),
    ]
    for path in paths:
        if path.exists():
            sys.path.insert(0, str(path.parent.resolve()))
            from newfamily_rank11_minimal_common import build_finite_minimal_family
            return build_finite_minimal_family
    raise RuntimeError("missing newfamily_rank11_minimal_common.py")


def load_sections(path: str):
    recovered = load(path)
    sections = recovered["sections"]
    if set(sections) != set(range(11)):
        raise RuntimeError("section sobj does not contain U0..U10")
    return [sections[i] for i in range(11)]


def qbits(value) -> int:
    value = QQ(value)
    return max(ZZ(abs(value.numerator())).nbits(), ZZ(value.denominator()).nbits())


def run_single(args) -> int:
    started = time.monotonic()
    t = QQ(args.numerator) / QQ(args.denominator)
    print("PHASE load", flush=True)
    sections = load_sections(args.sections_sobj)
    build = load_builder()
    family = build(ROOTS)

    print("PHASE evaluate_family", flush=True)
    A = QQ(family["Amin"](t))
    B = QQ(family["Bmin"](t))
    E = EllipticCurve(QQ, [0, 0, 0, A, B])
    if E.discriminant() == 0:
        print(json.dumps({"status": "singular", "parameter": str(t)}))
        return 0

    print("PHASE evaluate_sections", flush=True)
    points = []
    bits = []
    for xfun, yfun in sections:
        x = QQ(xfun(t))
        y = QQ(yfun(t))
        P = E([x, y])
        points.append(P)
        bits.append(max(qbits(x), qbits(y)))
    if len(set(points)) != 11:
        raise RuntimeError("hidden sections collide at this specialization")

    print(
        f"PHASE height_start precision={args.precision} max_bits={max(bits)} median_bits={sorted(bits)[5]}",
        flush=True,
    )
    hs = time.monotonic()
    G = E.height_pairing_matrix(points, precision=args.precision)
    height_seconds = time.monotonic() - hs
    print(f"PHASE height_done seconds={height_seconds:.6f}", flush=True)

    RF = RealField(args.precision)
    GR = matrix(RF, G)
    eigenvalues = sorted(RF(x) for x in GR.eigenvalues())
    largest = max(abs(x) for x in eigenvalues)
    tolerance = max(RF(2) ** (-args.rank_bits), largest * RF(2) ** (-args.rank_bits))
    numerical_rank = sum(abs(x) > tolerance for x in eigenvalues)
    positive_definite = all(x > tolerance for x in eigenvalues)
    regulator = GR.det()

    result = {
        "status": "completed",
        "parameter": f"{args.numerator}/{args.denominator}",
        "numerator": args.numerator,
        "denominator": args.denominator,
        "discovery": args.discovery,
        "held": args.held,
        "total": args.discovery + args.held,
        "minimum": min(args.discovery, args.held),
        "height_rank": int(numerical_rank),
        "positive_definite": bool(positive_definite),
        "height_precision": args.precision,
        "rank_tolerance": str(tolerance),
        "smallest_eigenvalue": str(eigenvalues[0]),
        "largest_eigenvalue": str(eigenvalues[-1]),
        "regulator": str(regulator),
        "regulator_log2": float(RR(abs(regulator)).log2()),
        "generator_bits_min": min(bits),
        "generator_bits_median": sorted(bits)[5],
        "generator_bits_max": max(bits),
        "height_seconds": height_seconds,
        "wall_seconds": time.monotonic() - started,
    }
    print(json.dumps(result, sort_keys=True))
    return 0


def run_parent(args) -> int:
    candidates = parse_scan([Path(x) for x in args.scan])
    shortlist = diversify(candidates, args.quota, args.small_denominator)
    if args.limit is not None:
        shortlist = shortlist[: args.limit]

    print(
        f"retained={len(candidates)} diversified={len(shortlist)} timeout={args.timeout} precision={args.precision}",
        flush=True,
    )
    records = []
    for pos, c in enumerate(shortlist, 1):
        print(
            f"[{pos}/{len(shortlist)}] T={c.name} D={c.discovery:.4f} H={c.held:.4f}",
            flush=True,
        )
        cmd = [
            sys.executable,
            str(Path(__file__).resolve()),
            "--single",
            "--sections-sobj", args.sections_sobj,
            "--numerator", str(c.numerator),
            "--denominator", str(c.denominator),
            "--discovery", repr(c.discovery),
            "--held", repr(c.held),
            "--precision", str(args.precision),
            "--rank-bits", str(args.rank_bits),
        ]
        try:
            completed = subprocess.run(
                cmd,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                timeout=args.timeout,
                check=False,
            )
        except subprocess.TimeoutExpired as exc:
            tail = (exc.stdout or "")
            if isinstance(tail, bytes):
                tail = tail.decode(errors="replace")
            phases = [line for line in tail.splitlines() if line.startswith("PHASE ")]
            record = {
                "status": "timeout",
                "parameter": c.name,
                "numerator": c.numerator,
                "denominator": c.denominator,
                "discovery": c.discovery,
                "held": c.held,
                "timeout_seconds": args.timeout,
                "phase_tail": phases[-8:],
            }
            records.append(record)
            print("  TIMEOUT last phases:", flush=True)
            for line in phases[-8:]:
                print("   ", line, flush=True)
            continue

        lines = [line for line in completed.stdout.splitlines() if line.strip()]
        try:
            record = json.loads(lines[-1])
        except Exception:
            record = {
                "status": "error",
                "parameter": c.name,
                "returncode": completed.returncode,
                "output_tail": "\n".join(lines[-20:]),
            }
        records.append(record)
        print(
            "  status=%s hrank=%s pd=%s bits=%s height_s=%s wall=%s" % (
                record.get("status"),
                record.get("height_rank"),
                record.get("positive_definite"),
                record.get("generator_bits_max"),
                record.get("height_seconds"),
                record.get("wall_seconds"),
            ),
            flush=True,
        )

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(records, indent=2, sort_keys=True) + "\n")
    completed = [r for r in records if r.get("status") == "completed"]
    good = [
        r for r in completed
        if r.get("height_rank") == 11 and r.get("positive_definite") is True
    ]
    print(json.dumps({
        "attempted": len(records),
        "completed": len(completed),
        "height_rank_11_positive_definite": len(good),
        "timeouts": sum(r.get("status") == "timeout" for r in records),
        "output": str(out),
    }, sort_keys=True))
    return 0


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--scan", action="append")
    p.add_argument("--sections-sobj", default="/tmp/newfamily_hidden_sections_complete.sobj")
    p.add_argument("--quota", type=int, default=8)
    p.add_argument("--small-denominator", type=int, default=100)
    p.add_argument("--limit", type=int)
    p.add_argument("--timeout", type=int, default=120)
    p.add_argument("--precision", type=int, default=180)
    p.add_argument("--rank-bits", type=int, default=80)
    p.add_argument(
        "--output",
        default="artifacts/local/elliptic-curves/newfamily/height_rank_rational_screen.json",
    )
    p.add_argument("--single", action="store_true", help=argparse.SUPPRESS)
    p.add_argument("--numerator", type=int, help=argparse.SUPPRESS)
    p.add_argument("--denominator", type=int, help=argparse.SUPPRESS)
    p.add_argument("--discovery", type=float, default=0.0, help=argparse.SUPPRESS)
    p.add_argument("--held", type=float, default=0.0, help=argparse.SUPPRESS)
    args = p.parse_args()
    if args.single:
        if args.numerator is None or args.denominator is None:
            p.error("single mode requires numerator and denominator")
        return run_single(args)
    if not args.scan:
        p.error("provide at least one --scan file")
    return run_parent(args)


if __name__ == "__main__":
    raise SystemExit(main())

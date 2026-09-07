#!/usr/bin/env sage -python
"""Seeded arithmetic screen for rational newfamily specializations.

The scanner only ranks local Nagao signal.  This replay takes a diversified
shortlist, specializes the eleven exact hidden generic sections, transports
them to a global minimal model, and asks eclib for the rank of that *known*
subgroup with saturation disabled.  A candidate is useful for the next phase
only if the baseline subgroup is recovered as rank 11 cheaply.

Every specialization is isolated in a child Sage process.  This prevents one
large global-minimal-model or eclib call from stalling the entire shortlist.
No bounded point search is performed here; this stage establishes a trustworthy
baseline for later rank-gain searches.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction
import json
from math import gcd
from pathlib import Path
import subprocess
import sys
import time

from sage.all import EllipticCurve, QQ, ZZ, load
from sage.libs.eclib.interface import mwrank_EllipticCurve, mwrank_MordellWeil

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
    def parameter(self) -> Fraction:
        return Fraction(self.numerator, self.denominator)

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
            candidate = Candidate(
                int(fields[1]), int(fields[2]), float(fields[3]), float(fields[4])
            )
            unique[(candidate.numerator, candidate.denominator)] = candidate
    return list(unique.values())


def diversify(candidates: list[Candidate], quota: int, small_denominator: int) -> list[Candidate]:
    selected: dict[tuple[int, int], Candidate] = {}

    def add(rows):
        for candidate in rows[:quota]:
            selected[(candidate.numerator, candidate.denominator)] = candidate

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
    candidates = [
        Path("elliptic-curves/cas/newfamily/newfamily_rank11_minimal_common.py"),
        Path("elliptic-curves/cas/newfamily/archive/newfamily_rank11_minimal_common.py"),
        Path("/tmp/newfamily_rank11_minimal_common.py"),
    ]
    for path in candidates:
        if path.exists():
            sys.path.insert(0, str(path.parent.resolve()))
            from newfamily_rank11_minimal_common import build_finite_minimal_family
            return build_finite_minimal_family
    raise RuntimeError("missing newfamily_rank11_minimal_common.py")


def qbits(value) -> int:
    value = QQ(value)
    return max(
        ZZ(abs(value.numerator())).nbits(),
        ZZ(value.denominator()).nbits(),
    )


def to_mwrank_triple(point) -> list:
    if point.is_zero():
        return [ZZ(0), ZZ(1), ZZ(0)]
    x = QQ(point[0])
    y = QQ(point[1])
    z = ZZ(x.denominator()).lcm(ZZ(y.denominator()))
    X = ZZ(x * z)
    Y = ZZ(y * z)
    common = gcd(gcd(abs(int(X)), abs(int(Y))), abs(int(z)))
    if common:
        X //= common
        Y //= common
        z //= common
    if z < 0:
        X, Y, z = -X, -Y, -z
    if QQ(X) / z != x or QQ(Y) / z != y:
        raise AssertionError("mwrank projective conversion failed")
    return [X, Y, z]


def load_sections_from_sobj(path: str):
    recovered = load(path)
    sections = recovered["sections"]
    if set(sections) != set(range(11)):
        raise RuntimeError("section sobj does not contain U0..U10")
    return [sections[index] for index in range(11)]


def run_single(args) -> int:
    started = time.monotonic()
    parameter = QQ(args.numerator) / QQ(args.denominator)
    sections = load_sections_from_sobj(args.sections_sobj)
    build = load_builder()
    family = build(ROOTS)

    A = QQ(family["Amin"](parameter))
    B = QQ(family["Bmin"](parameter))
    E = EllipticCurve(QQ, [0, 0, 0, A, B])
    if E.discriminant() == 0:
        print(json.dumps({"status": "singular", "parameter": str(parameter)}))
        return 0

    specialized = []
    for index, (x_function, y_function) in enumerate(sections):
        x = QQ(x_function(parameter))
        y = QQ(y_function(parameter))
        point = E([x, y])
        specialized.append(point)

    if len(set(specialized)) != 11:
        raise RuntimeError("hidden sections collide at this specialization")

    raw_bits = [max(qbits(P[0]), qbits(P[1])) for P in specialized]

    minimal_started = time.monotonic()
    Emin = E.global_minimal_model()
    minimal_seconds = time.monotonic() - minimal_started
    iso = E.isomorphism_to(Emin)
    points = [iso(P) for P in specialized]
    minimal_bits = [max(qbits(P[0]), qbits(P[1])) for P in points]

    order = sorted(range(11), key=lambda index: minimal_bits[index])
    mwcurve = mwrank_EllipticCurve([ZZ(value) for value in Emin.ainvs()])
    mw = mwrank_MordellWeil(mwcurve, verbose=False, pp=1, maxr=32)

    steps = []
    for index in order:
        step_started = time.monotonic()
        before = len(mw.points())
        mw.process([to_mwrank_triple(points[index])], saturation_bound=0)
        after = len(mw.points())
        steps.append(
            {
                "section": index,
                "bits": minimal_bits[index],
                "rank_before": before,
                "rank_after": after,
                "seconds": time.monotonic() - step_started,
            }
        )

    baseline_rank = len(mw.points())
    root_number = int(Emin.root_number())
    result = {
        "status": "completed",
        "parameter": f"{args.numerator}/{args.denominator}",
        "numerator": args.numerator,
        "denominator": args.denominator,
        "discovery": args.discovery,
        "held": args.held,
        "total": args.discovery + args.held,
        "minimum": min(args.discovery, args.held),
        "baseline_rank": baseline_rank,
        "root_number": root_number,
        "raw_generator_bits_min": min(raw_bits),
        "raw_generator_bits_median": sorted(raw_bits)[5],
        "raw_generator_bits_max": max(raw_bits),
        "minimal_generator_bits_min": min(minimal_bits),
        "minimal_generator_bits_median": sorted(minimal_bits)[5],
        "minimal_generator_bits_max": max(minimal_bits),
        "minimal_discriminant_bits": ZZ(abs(Emin.discriminant())).nbits(),
        "minimal_seconds": minimal_seconds,
        "generator_steps": steps,
        "wall_seconds": time.monotonic() - started,
    }
    print(json.dumps(result, sort_keys=True))
    return 0


def run_parent(args) -> int:
    scan_paths = [Path(value) for value in args.scan]
    candidates = parse_scan(scan_paths)
    shortlist = diversify(candidates, args.quota, args.small_denominator)
    if args.limit is not None:
        shortlist = shortlist[: args.limit]

    print(
        f"retained={len(candidates)} diversified={len(shortlist)} "
        f"timeout={args.timeout}",
        flush=True,
    )

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    records = []

    for position, candidate in enumerate(shortlist, 1):
        command = [
            sys.executable,
            str(Path(__file__).resolve()),
            "--single",
            "--sections-sobj", args.sections_sobj,
            "--numerator", str(candidate.numerator),
            "--denominator", str(candidate.denominator),
            "--discovery", repr(candidate.discovery),
            "--held", repr(candidate.held),
        ]
        print(
            f"[{position}/{len(shortlist)}] T={candidate.name} "
            f"D={candidate.discovery:.4f} H={candidate.held:.4f}",
            flush=True,
        )
        try:
            completed = subprocess.run(
                command,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                timeout=args.timeout,
                check=False,
            )
        except subprocess.TimeoutExpired:
            record = {
                "status": "timeout",
                "parameter": candidate.name,
                "numerator": candidate.numerator,
                "denominator": candidate.denominator,
                "discovery": candidate.discovery,
                "held": candidate.held,
                "timeout_seconds": args.timeout,
            }
            records.append(record)
            print("  TIMEOUT", flush=True)
            continue

        lines = [line for line in completed.stdout.splitlines() if line.strip()]
        try:
            record = json.loads(lines[-1])
        except Exception:
            record = {
                "status": "error",
                "parameter": candidate.name,
                "returncode": completed.returncode,
                "output_tail": "\n".join(lines[-20:]),
            }
        records.append(record)
        print(
            "  status=%s baseline=%s root=%s bits=%s wall=%s" % (
                record.get("status"),
                record.get("baseline_rank"),
                record.get("root_number"),
                record.get("minimal_generator_bits_max"),
                record.get("wall_seconds"),
            ),
            flush=True,
        )

    output.write_text(json.dumps(records, indent=2, sort_keys=True) + "\n")

    completed = [record for record in records if record.get("status") == "completed"]
    good = [record for record in completed if record.get("baseline_rank") == 11]
    print(
        json.dumps(
            {
                "output": str(output),
                "attempted": len(records),
                "completed": len(completed),
                "baseline_rank_11": len(good),
                "timeouts": sum(record.get("status") == "timeout" for record in records),
            },
            sort_keys=True,
        )
    )
    return 0


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--scan", action="append")
    parser.add_argument(
        "--sections-sobj",
        default="/tmp/newfamily_hidden_sections_complete.sobj",
    )
    parser.add_argument("--quota", type=int, default=8)
    parser.add_argument("--small-denominator", type=int, default=100)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--timeout", type=int, default=120)
    parser.add_argument(
        "--output",
        default="artifacts/local/elliptic-curves/newfamily/seeded_rational_screen.json",
    )

    parser.add_argument("--single", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--numerator", type=int, help=argparse.SUPPRESS)
    parser.add_argument("--denominator", type=int, help=argparse.SUPPRESS)
    parser.add_argument("--discovery", type=float, default=0.0, help=argparse.SUPPRESS)
    parser.add_argument("--held", type=float, default=0.0, help=argparse.SUPPRESS)
    args = parser.parse_args()

    if args.single:
        if args.numerator is None or args.denominator is None:
            parser.error("single mode requires numerator and denominator")
        return run_single(args)

    if not args.scan:
        parser.error("provide at least one --scan file")
    return run_parent(args)


if __name__ == "__main__":
    raise SystemExit(main())

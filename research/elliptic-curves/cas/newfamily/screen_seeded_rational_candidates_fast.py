#!/usr/bin/env sage -python
"""Fast seeded baseline screen for rational newfamily specializations.

Unlike ``screen_seeded_rational_candidates.py``, this stage deliberately avoids
``global_minimal_model()`` and root-number computation.  For T=a/b it uses the
short-model isomorphism

    X = b^4 x,  Y = b^6 y,

so that

    Y^2 = X^3 + b^8 A(a/b) X + b^12 B(a/b).

The finite-minimal family has degree(A)<=8 and degree(B)<=12, so this is the
natural homogeneous integral specialization when the polynomial coefficients
are integral.  If a fixed rational coefficient denominator remains, a small
extra integral scale c is applied:

    X <- c^2 X, Y <- c^3 Y,
    A <- c^4 A, B <- c^6 B.

The goal is only to establish the exact known subgroup rank 11 cheaply.  Global
minimalization, root number, conductor and deep searches belong to later phases.
Each candidate runs in a child process; phase markers are flushed so a timeout
reports where it happened.
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
    def total(self):
        return self.discovery + self.held

    @property
    def minimum(self):
        return min(self.discovery, self.held)

    @property
    def balanced(self):
        return 0.35 * self.discovery + 0.65 * self.held

    @property
    def name(self):
        return f"{self.numerator}/{self.denominator}"


def parse_scan(paths):
    unique = {}
    for path in paths:
        for line in Path(path).read_text().splitlines():
            if not line.startswith("C "):
                continue
            f = line.split()
            c = Candidate(int(f[1]), int(f[2]), float(f[3]), float(f[4]))
            unique[(c.numerator, c.denominator)] = c
    return list(unique.values())


def diversify(candidates, quota, small_denominator):
    selected = {}

    def add(rows):
        for c in rows[:quota]:
            selected[(c.numerator, c.denominator)] = c

    add(sorted(candidates, key=lambda c: c.total, reverse=True))
    add(sorted(candidates, key=lambda c: c.held, reverse=True))
    add(sorted(candidates, key=lambda c: c.minimum, reverse=True))
    add(sorted(candidates, key=lambda c: c.balanced, reverse=True))
    add(sorted(
        (c for c in candidates if c.denominator <= small_denominator),
        key=lambda c: (c.total, c.held), reverse=True,
    ))
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


def qbits(value):
    value = QQ(value)
    return max(ZZ(abs(value.numerator())).nbits(), ZZ(value.denominator()).nbits())


def to_mwrank_triple(point):
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
    assert QQ(X) / z == x and QQ(Y) / z == y
    return [X, Y, z]


def load_sections(path):
    recovered = load(path)
    sections = recovered["sections"]
    if set(sections) != set(range(11)):
        raise RuntimeError("section sobj does not contain U0..U10")
    return [sections[i] for i in range(11)]


def phase(name, **extra):
    fields = " ".join(f"{k}={v}" for k, v in extra.items())
    print(f"PHASE {name}" + (" " + fields if fields else ""), flush=True)


def extra_integral_scale(Ah, Bh):
    """Return c with Ah*c^4 and Bh*c^6 integral.

    Usually c=1 because b^8 and b^12 already homogenize the polynomial family.
    We intentionally avoid factoring denominators: lcm is a safe fallback and
    these are fixed coefficient denominators, not the large b-powers.
    """
    da = ZZ(QQ(Ah).denominator())
    db = ZZ(QQ(Bh).denominator())
    if da == 1 and db == 1:
        return ZZ(1)
    return da.lcm(db)


def run_single(args):
    started = time.monotonic()
    a = ZZ(args.numerator)
    b = ZZ(args.denominator)
    t = QQ(a) / QQ(b)

    phase("load")
    sections = load_sections(args.sections_sobj)
    build = load_builder()
    family = build(ROOTS)

    phase("evaluate_family")
    A = QQ(family["Amin"](t))
    B = QQ(family["Bmin"](t))

    # Homogeneous specialization of the short model.
    Ah = A * b**8
    Bh = B * b**12
    c = extra_integral_scale(Ah, Bh)
    Aint_q = Ah * c**4
    Bint_q = Bh * c**6
    if Aint_q.denominator() != 1 or Bint_q.denominator() != 1:
        raise RuntimeError("failed to integralize homogeneous short model")
    Aint = ZZ(Aint_q)
    Bint = ZZ(Bint_q)

    phase(
        "integral_model",
        b_bits=b.nbits(),
        extra_scale_bits=c.nbits(),
        A_bits=abs(Aint).nbits(),
        B_bits=abs(Bint).nbits(),
    )
    E = EllipticCurve(QQ, [0, 0, 0, Aint, Bint])
    if E.discriminant() == 0:
        print(json.dumps({"status": "singular", "parameter": str(t)}), flush=True)
        return 0

    # Total point scaling from original rational short model.
    xscale = b**4 * c**2
    yscale = b**6 * c**3

    phase("evaluate_sections")
    points = []
    bits = []
    for index, (xf, yf) in enumerate(sections):
        x = QQ(xf(t)) * xscale
        y = QQ(yf(t)) * yscale
        P = E([x, y])
        points.append(P)
        bits.append(max(qbits(x), qbits(y)))
    if len(set(points)) != 11:
        raise RuntimeError("hidden sections collide at this specialization")

    order = sorted(range(11), key=lambda i: bits[i])
    phase("eclib_init", max_bits=max(bits), median_bits=sorted(bits)[5])
    mwcurve = mwrank_EllipticCurve([ZZ(v) for v in E.ainvs()])
    mw = mwrank_MordellWeil(mwcurve, verbose=False, pp=1, maxr=32)

    steps = []
    for position, index in enumerate(order, 1):
        step_started = time.monotonic()
        before = len(mw.points())
        phase("eclib_point_start", position=position, section=index, bits=bits[index], rank=before)
        mw.process([to_mwrank_triple(points[index])], saturation_bound=0)
        after = len(mw.points())
        seconds = time.monotonic() - step_started
        phase("eclib_point_done", position=position, section=index, rank=after, seconds=f"{seconds:.6f}")
        steps.append({
            "section": index,
            "bits": bits[index],
            "rank_before": before,
            "rank_after": after,
            "seconds": seconds,
        })

    baseline_rank = len(mw.points())
    result = {
        "status": "completed",
        "model": "homogeneous_integral_short",
        "parameter": f"{args.numerator}/{args.denominator}",
        "numerator": args.numerator,
        "denominator": args.denominator,
        "discovery": args.discovery,
        "held": args.held,
        "total": args.discovery + args.held,
        "minimum": min(args.discovery, args.held),
        "baseline_rank": baseline_rank,
        "generator_bits_min": min(bits),
        "generator_bits_median": sorted(bits)[5],
        "generator_bits_max": max(bits),
        "integral_a4_bits": abs(Aint).nbits(),
        "integral_a6_bits": abs(Bint).nbits(),
        "extra_scale": str(c),
        "generator_steps": steps,
        "wall_seconds": time.monotonic() - started,
    }
    print(json.dumps(result, sort_keys=True), flush=True)
    return 0


def partial_text(value):
    if value is None:
        return ""
    if isinstance(value, bytes):
        return value.decode(errors="replace")
    return str(value)


def run_parent(args):
    candidates = parse_scan(args.scan)
    shortlist = diversify(candidates, args.quota, args.small_denominator)
    if args.limit is not None:
        shortlist = shortlist[:args.limit]

    print(f"retained={len(candidates)} diversified={len(shortlist)} timeout={args.timeout}", flush=True)
    records = []
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)

    for position, candidate in enumerate(shortlist, 1):
        print(
            f"[{position}/{len(shortlist)}] T={candidate.name} "
            f"D={candidate.discovery:.4f} H={candidate.held:.4f}",
            flush=True,
        )
        command = [
            sys.executable, str(Path(__file__).resolve()), "--single",
            "--sections-sobj", args.sections_sobj,
            "--numerator", str(candidate.numerator),
            "--denominator", str(candidate.denominator),
            "--discovery", repr(candidate.discovery),
            "--held", repr(candidate.held),
        ]
        try:
            completed = subprocess.run(
                command,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                timeout=args.timeout,
                check=False,
            )
            lines = [line for line in completed.stdout.splitlines() if line.strip()]
            try:
                record = json.loads(lines[-1])
            except Exception:
                record = {
                    "status": "error",
                    "parameter": candidate.name,
                    "returncode": completed.returncode,
                    "output_tail": "\n".join(lines[-30:]),
                }
        except subprocess.TimeoutExpired as exc:
            tail = partial_text(exc.stdout).splitlines()[-30:]
            record = {
                "status": "timeout",
                "parameter": candidate.name,
                "numerator": candidate.numerator,
                "denominator": candidate.denominator,
                "discovery": candidate.discovery,
                "held": candidate.held,
                "timeout_seconds": args.timeout,
                "phase_tail": tail,
            }
            print("  TIMEOUT last phases:", flush=True)
            for line in tail[-8:]:
                print("   ", line, flush=True)
            records.append(record)
            output.write_text(json.dumps(records, indent=2, sort_keys=True) + "\n")
            continue

        records.append(record)
        output.write_text(json.dumps(records, indent=2, sort_keys=True) + "\n")
        print(
            "  status=%s baseline=%s bits=%s wall=%s" % (
                record.get("status"), record.get("baseline_rank"),
                record.get("generator_bits_max"), record.get("wall_seconds"),
            ),
            flush=True,
        )

    completed = [r for r in records if r.get("status") == "completed"]
    good = [r for r in completed if r.get("baseline_rank") == 11]
    print(json.dumps({
        "output": str(output),
        "attempted": len(records),
        "completed": len(completed),
        "baseline_rank_11": len(good),
        "timeouts": sum(r.get("status") == "timeout" for r in records),
    }, sort_keys=True), flush=True)
    return 0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--scan", action="append")
    parser.add_argument("--sections-sobj", default="/tmp/newfamily_hidden_sections_complete.sobj")
    parser.add_argument("--quota", type=int, default=8)
    parser.add_argument("--small-denominator", type=int, default=100)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--timeout", type=int, default=120)
    parser.add_argument(
        "--output",
        default="artifacts/local/elliptic-curves/newfamily/seeded_rational_screen_fast.json",
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
        parser.error("provide at least one --scan")
    return run_parent(args)


if __name__ == "__main__":
    raise SystemExit(main())

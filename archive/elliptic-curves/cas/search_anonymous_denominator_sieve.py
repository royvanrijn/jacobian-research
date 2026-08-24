#!/usr/bin/env python3
from __future__ import annotations

import argparse
from fractions import Fraction
from math import gcd, isqrt
from pathlib import Path
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent))

from anonymous_rank_candidate import A4 as A, A6 as B


PROTOCOL = "ANONSIEVE"

# The square condition is valid at bad primes too:
# an integer square remains a square mod p.
SIEVE_PRIMES = (
    3, 5, 7, 11, 13, 17, 19, 23, 29, 31,
    37, 41, 43, 47, 53, 59, 61, 67, 71, 73,
    79, 83, 89, 97, 101, 103, 107, 109,
    113, 127, 131, 137, 139, 149,
)


def icbrt(n: int) -> int:
    """Exact floor cube root."""
    if n < 0:
        return -icbrt(-n)
    if n < 2:
        return n

    x = 1 << ((n.bit_length() + 2) // 3)

    while True:
        y = (2 * x + n // (x * x)) // 3
        if y >= x:
            break
        x = y

    while (x + 1) ** 3 <= n:
        x += 1
    while x ** 3 > n:
        x -= 1

    return x


def f_integer(x: int) -> int:
    """
    For z = 2y+x:

        z^2 = 4x^3 + x^2 + 4Ax + 4B
    """
    return 4 * x**3 + x*x + 4*A*x + 4*B


def real_root_center(scale: int) -> int:
    """
    Locate the real zero of F(x) by integer bisection.
    """
    lo = -8 * scale
    hi = 8 * scale

    while f_integer(lo) > 0:
        lo *= 2

    while f_integer(hi) < 0:
        hi *= 2

    while hi - lo > 1:
        mid = (lo + hi) // 2

        if f_integer(mid) <= 0:
            lo = mid
        else:
            hi = mid

    if abs(f_integer(lo)) <= abs(f_integer(hi)):
        return lo

    return hi


def intrinsic_centers():
    """
    Geometry-derived x centers.

    F'(x) = 12x^2 + 2x + 4A

    so its stationary points are approximately

        (-1 +/- sqrt(1-48A))/12.
    """

    bscale = icbrt(abs(B))
    rr = real_root_center(max(1, bscale))

    c4 = 1 - 48*A
    s = isqrt(c4)

    # nearest integer sqrt
    if (s + 1) * (s + 1) - c4 < c4 - s*s:
        s += 1

    stat_pos = (-1 + s) // 12
    stat_neg = (-1 - s) // 12

    centers = [
        ("zero", 0),
        ("realroot", rr),
        ("stationary_neg", stat_neg),
        ("stationary_pos", stat_pos),
        ("bcuberoot_neg", -bscale),
        ("bcuberoot_pos", bscale),
    ]

    seen = set()
    out = []

    for label, center in centers:
        if center in seen:
            continue

        seen.add(center)
        out.append((label, center))

    return out


def repeat_pattern(pattern: int, p: int, length: int) -> int:
    """
    Repeat a p-bit residue pattern over `length` positions.

    This is much faster than iterating over every k.
    """

    reps = (length + p - 1) // p

    repunit = (
        ((1 << (p * reps)) - 1)
        // ((1 << p) - 1)
    )

    return (
        pattern * repunit
    ) & ((1 << length) - 1)


def allowed_mask(
    center: int,
    scale: int,
    d: int,
    p: int,
    radius: int,
    length: int,
) -> int:
    """
    Parameterization

        x = center + scale*k/d^2
          = a/d^2

        a = center*d^2 + scale*k.

    We require

        H = d^6 z^2
          = 4a^3
            + a^2 d^2
            + 4A a d^4
            + 4B d^6

    to be a square modulo p.
    """

    dmod = d % p

    d2 = dmod * dmod % p
    d4 = d2 * d2 % p
    d6 = d4 * d2 % p

    center_mod = center % p
    scale_mod = scale % p

    a0 = center_mod * d2 % p

    A_mod = A % p
    B_mod = B % p

    squares = {
        x*x % p
        for x in range(p)
    }

    pattern = 0
    radius_mod = radius % p

    # index = k + radius
    for index_residue in range(p):
        k_residue = (
            index_residue - radius_mod
        ) % p

        a = (
            a0
            + scale_mod * k_residue
        ) % p

        h = (
            4 * a*a % p * a
            + a*a * d2
            + 4 * A_mod * a * d4
            + 4 * B_mod * d6
        ) % p

        if h in squares:
            pattern |= 1 << index_residue

    return repeat_pattern(
        pattern,
        p,
        length,
    )


def exact_square_for(
    center: int,
    scale: int,
    d: int,
    k: int,
):
    d2 = d*d

    a = (
        center * d2
        + scale * k
    )

    d4 = d2*d2
    d6 = d4*d2

    h = (
        4*a**3
        + a*a*d2
        + 4*A*a*d4
        + 4*B*d6
    )

    if h < 0:
        return None

    r = isqrt(h)

    if r*r != h:
        return None

    return a, d, r


def map_point(
    a: int,
    d: int,
    r: int,
):
    x = Fraction(a, d*d)
    z = Fraction(r, d*d*d)

    p1 = (
        x,
        (-x + z) / 2,
    )

    p2 = (
        x,
        (-x - z) / 2,
    )

    for xx, yy in (p1, p2):
        assert (
            yy*yy + xx*yy
            ==
            xx**3 + A*xx + B
        )

    return p1, p2


def parse_ints(text: str):
    return [
        int(v.strip())
        for v in text.split(",")
        if v.strip()
    ]


def main():
    ap = argparse.ArgumentParser(
        description=(
            "Denominator-normalized quadratic-residue sieve "
            "for the anonymous high-rank candidate."
        )
    )

    ap.add_argument(
        "--d-min",
        type=int,
        default=1,
    )

    ap.add_argument(
        "--d-max",
        type=int,
        default=500,
    )

    ap.add_argument(
        "--radius",
        type=int,
        default=50000,
    )

    ap.add_argument(
        "--scales",
        default=(
            "1,"
            "1000000,"
            "1000000000000,"
            "1000000000000000000,"
            "1000000000000000000000000,"
            "1000000000000000000000000000"
        ),
    )

    ap.add_argument(
        "--centers",
        default="all",
        help=(
            "all or comma-separated intrinsic labels: "
            "zero,realroot,stationary_neg,stationary_pos,"
            "bcuberoot_neg,bcuberoot_pos"
        ),
    )

    ap.add_argument(
        "--prime-count",
        type=int,
        default=len(SIEVE_PRIMES),
        help=(
            f"prefix of the {len(SIEVE_PRIMES)} "
            "pinned residue-sieve primes"
        ),
    )

    ap.add_argument(
        "--out",
        default=(
            "artifacts/local/elliptic-curves/"
            "anonymous-denominator-sieve-points.tsv"
        ),
    )

    args = ap.parse_args()

    if args.d_min < 1:
        raise SystemExit("d-min must be >= 1")

    if args.d_max < args.d_min:
        raise SystemExit("bad denominator interval")

    if args.radius < 1:
        raise SystemExit("radius must be positive")

    if not (
        1
        <= args.prime_count
        <= len(SIEVE_PRIMES)
    ):
        raise SystemExit("bad --prime-count")

    scales = parse_ints(args.scales)

    if (
        not scales
        or any(s <= 0 for s in scales)
    ):
        raise SystemExit(
            "scales must all be positive"
        )

    centers = intrinsic_centers()

    if args.centers != "all":
        wanted = {
            value.strip()
            for value
            in args.centers.split(",")
            if value.strip()
        }

        centers = [
            (label, center)
            for label, center in centers
            if label in wanted
        ]

        missing = (
            wanted
            - {label for label, _ in centers}
        )

        if missing:
            raise SystemExit(
                f"unknown centers: {sorted(missing)}"
            )

    if not centers:
        raise SystemExit("no centers selected")

    primes = SIEVE_PRIMES[
        :args.prime_count
    ]

    radius = args.radius
    length = 2*radius + 1

    full_mask = (
        1 << length
    ) - 1

    # k=0 was already covered by the center
    # itself; remove it.
    zero_bit = 1 << radius

    print(
        f"{PROTOCOL}|stage=start"
        f"|d={args.d_min}..{args.d_max}"
        f"|radius={radius}"
        f"|scales={scales}"
        f"|prime_count={len(primes)}"
        f"|centers="
        + ",".join(
            label
            for label, _ in centers
        ),
        flush=True,
    )

    for label, center in centers:
        print(
            f"{PROTOCOL}|stage=center"
            f"|label={label}"
            f"|x={center}"
            f"|F={f_integer(center)}",
            flush=True,
        )

    raw_per_box = length - 1

    denominator_count = (
        args.d_max
        - args.d_min
        + 1
    )

    boxes = (
        len(centers)
        * len(scales)
        * denominator_count
    )

    declared = (
        boxes
        * raw_per_box
    )

    print(
        f"{PROTOCOL}|stage=space"
        f"|boxes={boxes}"
        f"|declared_candidates={declared}",
        flush=True,
    )

    seen_x = set()
    found = []

    exact_tests = 0
    boxes_nonempty = 0

    t0 = time.monotonic()
    last_report = t0

    done_boxes = 0

    for label, center in centers:
        for scale in scales:
            for d in range(
                args.d_min,
                args.d_max + 1,
            ):
                mask = (
                    full_mask
                    ^ zero_bit
                )

                for p in primes:
                    mask &= allowed_mask(
                        center,
                        scale,
                        d,
                        p,
                        radius,
                        length,
                    )

                    if not mask:
                        break

                done_boxes += 1

                if mask:
                    boxes_nonempty += 1

                # Enumerate only the surviving bits.
                while mask:
                    low = mask & -mask

                    index = (
                        low.bit_length() - 1
                    )

                    mask ^= low

                    k = index - radius

                    # If gcd(k,d)>1 then this is an
                    # unnecessarily unreduced copy.
                    if gcd(abs(k), d) != 1:
                        continue

                    exact_tests += 1

                    sq = exact_square_for(
                        center,
                        scale,
                        d,
                        k,
                    )

                    if sq is None:
                        continue

                    a, dd, r = sq

                    points = map_point(
                        a,
                        dd,
                        r,
                    )

                    x = points[0][0]

                    if x in seen_x:
                        continue

                    seen_x.add(x)

                    record = (
                        label,
                        center,
                        scale,
                        d,
                        k,
                        points[0],
                        points[1],
                    )

                    found.append(record)

                    print(
                        f"{PROTOCOL}|stage=FOUND"
                        f"|center={label}"
                        f"|scale={scale}"
                        f"|d={d}"
                        f"|k={k}"
                        f"|x={points[0][0]}"
                        f"|y={points[0][1]}",
                        flush=True,
                    )

                now = time.monotonic()

                if now - last_report >= 5:
                    print(
                        f"{PROTOCOL}|stage=progress"
                        f"|boxes={done_boxes}/{boxes}"
                        f"|nonempty={boxes_nonempty}"
                        f"|exact_tests={exact_tests}"
                        f"|found={len(found)}"
                        f"|seconds={now-t0:.2f}",
                        flush=True,
                    )

                    last_report = now

    out = Path(args.out)

    out.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    with out.open("w") as h:
        h.write(
            "center\tscale\td\tk\tx\ty\n"
        )

        for (
            label,
            center,
            scale,
            d,
            k,
            p1,
            p2,
        ) in found:

            h.write(
                f"{label}\t"
                f"{scale}\t"
                f"{d}\t"
                f"{k}\t"
                f"{p1[0]}\t"
                f"{p1[1]}\n"
            )

            h.write(
                f"{label}\t"
                f"{scale}\t"
                f"{d}\t"
                f"{k}\t"
                f"{p2[0]}\t"
                f"{p2[1]}\n"
            )

    elapsed = (
        time.monotonic() - t0
    )

    print(
        f"{PROTOCOL}|stage=done"
        f"|status=PASS"
        f"|boxes={boxes}"
        f"|declared_candidates={declared}"
        f"|boxes_nonempty={boxes_nonempty}"
        f"|exact_tests={exact_tests}"
        f"|distinct_x={len(found)}"
        f"|seconds={elapsed:.3f}"
        f"|out={out}",
        flush=True,
    )


if __name__ == "__main__":
    main()

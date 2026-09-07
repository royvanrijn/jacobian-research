#!/usr/bin/env python3

from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
from pathlib import Path
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent))

from search_icarm_curve273_rank31 import (
    Q,
    PUBLISHED_POINTS,
    QuarticChart,
    affine_substitute,
    discover_relation,
    exact_linear_combination,
    point_negate,
    search_chart,
    short_weierstrass_coefficients,
    slope_discriminant,
    slope_to_points,
    to_short_point,
    x_chart_to_points,
    x_polynomial,
)

from mod2_reduction_independence import (
    combined_mod2_rank,
    find_mod2_reduction_certificate,
)


PROTOCOL = "R31LATTICE"

ROOT = Path(__file__).resolve().parents[2]

DEFAULT_VECTORS = (
    ROOT /
    "artifacts/local/elliptic-curves/curve273-rank30/"
    "short-coefficient-vectors.tsv"
)


@dataclass(frozen=True)
class Plan:
    chart: QuarticChart
    map_kind: str
    base_point: tuple[Fraction, Fraction] | None


def rat_bits(q):
    q = Q(q)
    return max(
        abs(q.numerator).bit_length(),
        q.denominator.bit_length(),
    )


def poly_bits(poly):
    return max(rat_bits(x) for x in poly)


def canonical_point(P):
    N = point_negate(P)
    return min(P, N)


def load_vectors(path, count):
    rows = []

    lines = path.read_text().splitlines()

    for line in lines[1:]:
        if not line.strip():
            continue

        parts = line.split("\t")

        height = float(parts[0])
        coeffs = tuple(int(x) for x in parts[1:])

        if len(coeffs) != 30:
            raise ValueError("short-vector row width changed")

        rows.append((height, coeffs))

        if len(rows) >= count:
            break

    return rows


def exact_subgroup_bases(rows, max_coordinate_bits):
    result = []

    for height, coeffs in rows:
        P = exact_linear_combination(coeffs)

        if P is None:
            raise AssertionError(
                "nonzero coefficient vector vanished despite rank-30 certificate"
            )

        bits = max(
            rat_bits(P[0]),
            rat_bits(P[1]),
        )

        if bits > max_coordinate_bits:
            continue

        result.append(
            (height, coeffs, P, bits)
        )

    return result


def secant_slope(P, Qp):
    if P[0] == Qp[0]:
        return None

    return (
        Q(Qp[1]) - Q(P[1])
    ) / (
        Q(Qp[0]) - Q(P[0])
    )


def signed_basis_points():
    out = []

    for P in PUBLISHED_POINTS:
        out.append(P)
        out.append(point_negate(P))

    return tuple(out)


SIGNED_BASIS = signed_basis_points()


def compression_score(values):
    hs = sorted(rat_bits(x) for x in values)

    return (
        hs[len(hs) // 4],
        hs[len(hs) // 2],
        hs[-1],
        sum(hs),
    )


def plans_for_base(
    base_id,
    base,
    *,
    x_pair_count,
    slope_anchor_count,
    slope_pair_count,
    x_offset_height,
    x_pair_height,
    slope_offset_height,
    slope_pair_height,
):
    plans = []

    xp = x_polynomial()

    # --------------------------------------------------------
    # Direct x-offset chart at the short subgroup point.
    # --------------------------------------------------------

    chart = QuarticChart(
        identifier=f"{base_id}_xoffset",
        kind="short_x_offset",
        polynomial=affine_substitute(
            xp,
            base[0],
            Q(1),
        ),
        center=base[0],
        scale=Q(1),
        base_index=None,
        seed_parameters=(Q(0),),
        height_bound=x_offset_height,
    )

    plans.append(
        Plan(chart, "x", None)
    )

    # --------------------------------------------------------
    # x-pair charts between the subgroup base and original
    # generators. Rank them by transformed coefficient size.
    # --------------------------------------------------------

    xranks = []

    for index, P in enumerate(PUBLISHED_POINTS):
        scale = Q(P[0]) - Q(base[0])

        if scale == 0:
            continue

        poly = affine_substitute(
            xp,
            base[0],
            scale,
        )

        xranks.append(
            (
                poly_bits(poly),
                index,
                scale,
                poly,
            )
        )

    xranks.sort(
        key=lambda x: (x[0], x[1])
    )

    for _, index, scale, poly in xranks[:x_pair_count]:
        chart = QuarticChart(
            identifier=(
                f"{base_id}_xpair_p{index+1:02d}"
            ),
            kind="short_x_pair",
            polynomial=poly,
            center=base[0],
            scale=scale,
            base_index=None,
            seed_parameters=(Q(0), Q(1)),
            height_bound=x_pair_height,
        )

        plans.append(
            Plan(chart, "x", None)
        )

    # --------------------------------------------------------
    # Slope quartic at this exact subgroup point.
    # --------------------------------------------------------

    D = slope_discriminant(base)

    anchors = []

    for index, P in enumerate(SIGNED_BASIS):
        m = secant_slope(base, P)

        if m is None:
            continue

        anchors.append(
            (
                rat_bits(m),
                index,
                m,
            )
        )

    anchors.sort(
        key=lambda x: (x[0], x[1])
    )

    # Remove duplicate slope values.
    distinct = []
    seen = set()

    for _, index, m in anchors:
        if m in seen:
            continue
        seen.add(m)
        distinct.append((index, m))

    pool = distinct[:slope_anchor_count]

    # Individual slope-offset charts.
    for index, m in pool:
        chart = QuarticChart(
            identifier=(
                f"{base_id}_soffset_s{index+1:02d}"
            ),
            kind="short_slope_offset",
            polynomial=affine_substitute(
                D,
                m,
                Q(1),
            ),
            center=m,
            scale=Q(1),
            base_index=None,
            seed_parameters=(Q(0),),
            height_bound=slope_offset_height,
        )

        plans.append(
            Plan(chart, "slope", base)
        )

    # Pair-normalized slope charts.
    ranked_pairs = []

    all_slopes = [
        m for _, m in distinct
    ]

    for (
        (i1, m1),
        (i2, m2),
    ) in combinations(pool, 2):

        scale = m2 - m1

        if scale == 0:
            continue

        normalized = tuple(
            (m - m1) / scale
            for m in all_slopes
        )

        poly = affine_substitute(
            D,
            m1,
            scale,
        )

        score = (
            compression_score(normalized),
            poly_bits(poly),
            i1,
            i2,
        )

        ranked_pairs.append(
            (
                score,
                i1,
                i2,
                m1,
                scale,
                poly,
            )
        )

    ranked_pairs.sort(
        key=lambda x: x[0]
    )

    for (
        _,
        i1,
        i2,
        center,
        scale,
        poly,
    ) in ranked_pairs[:slope_pair_count]:

        chart = QuarticChart(
            identifier=(
                f"{base_id}_spair_s{i1+1:02d}_s{i2+1:02d}"
            ),
            kind="short_slope_pair",
            polynomial=poly,
            center=center,
            scale=scale,
            base_index=None,
            seed_parameters=(Q(0), Q(1)),
            height_bound=slope_pair_height,
        )

        plans.append(
            Plan(chart, "slope", base)
        )

    return plans


def map_images(plan, quartic_points):
    images = []

    chart = plan.chart

    for parameter, square_root in quartic_points:
        if parameter in chart.seed_parameters:
            continue

        coordinate = (
            chart.center
            + chart.scale * parameter
        )

        if plan.map_kind == "x":
            mapped = x_chart_to_points(
                coordinate,
                square_root,
            )

        elif plan.map_kind == "slope":
            mapped = slope_to_points(
                plan.base_point,
                coordinate,
                square_root,
            )

        else:
            raise AssertionError(plan.map_kind)

        images.extend(mapped)

    return tuple(images)


def classify_candidate(
    point,
    *,
    relation_timeout,
    stack_bytes,
    certificate_prime_bound,
):
    relation = discover_relation(
        point,
        timeout=relation_timeout,
        stack_bytes=stack_bytes,
    )

    if relation is not None:
        return (
            "rank30_subgroup",
            30,
            relation,
            (),
        )

    augmented = tuple(
        to_short_point(P)
        for P in PUBLISHED_POINTS
    ) + (
        to_short_point(point),
    )

    signatures = find_mod2_reduction_certificate(
        short_weierstrass_coefficients(),
        augmented,
        prime_bound=certificate_prime_bound,
    )

    rank = combined_mod2_rank(
        signatures,
        31,
    )

    return (
        (
            "EXACT_INDEPENDENT_31ST"
            if rank == 31
            else "unresolved"
        ),
        rank,
        None,
        tuple(s.prime for s in signatures),
    )


def retry_old_timeouts(
    *,
    timeout,
    stack_bytes,
    x_offset_height,
):
    xp = x_polynomial()
    found = []

    # Previous tier timed out only on p28 and p29.
    for index in (27, 28):
        P = PUBLISHED_POINTS[index]

        chart = QuarticChart(
            identifier=f"retry_xoffset_p{index+1:02d}",
            kind="retry_basis_xoffset",
            polynomial=affine_substitute(
                xp,
                P[0],
                Q(1),
            ),
            center=P[0],
            scale=Q(1),
            base_index=None,
            seed_parameters=(Q(0),),
            height_bound=x_offset_height,
        )

        try:
            qpts, ms, wall = search_chart(
                chart,
                timeout=timeout,
                stack_bytes=stack_bytes,
            )
        except subprocess.TimeoutExpired:
            print(
                f"{PROTOCOL}|retry"
                f"|p={index+1}"
                f"|status=TIMEOUT",
                flush=True,
            )
            continue

        images = []

        for u, z in qpts:
            if u == 0:
                continue

            x = P[0] + u

            images.extend(
                x_chart_to_points(x, z)
            )

        print(
            f"{PROTOCOL}|retry"
            f"|p={index+1}"
            f"|images={len(images)}"
            f"|pari_ms={ms}"
            f"|wall={wall:.3f}",
            flush=True,
        )

        found.extend(images)

    return found


def main():
    ap = argparse.ArgumentParser()

    ap.add_argument(
        "--vectors",
        type=Path,
        default=DEFAULT_VECTORS,
    )

    ap.add_argument(
        "--base-count",
        type=int,
        default=120,
    )

    ap.add_argument(
        "--candidate-pool",
        type=int,
        default=600,
    )

    ap.add_argument(
        "--max-coordinate-bits",
        type=int,
        default=1400,
    )

    ap.add_argument(
        "--x-pair-count",
        type=int,
        default=4,
    )

    ap.add_argument(
        "--slope-anchor-count",
        type=int,
        default=8,
    )

    ap.add_argument(
        "--slope-pair-count",
        type=int,
        default=12,
    )

    ap.add_argument(
        "--x-offset-height",
        type=int,
        default=200000,
    )

    ap.add_argument(
        "--x-pair-height",
        type=int,
        default=20000,
    )

    ap.add_argument(
        "--slope-offset-height",
        type=int,
        default=20000,
    )

    ap.add_argument(
        "--slope-pair-height",
        type=int,
        default=20000,
    )

    ap.add_argument(
        "--chart-timeout",
        type=float,
        default=4.0,
    )

    ap.add_argument(
        "--retry-timeout",
        type=float,
        default=30.0,
    )

    ap.add_argument(
        "--relation-timeout",
        type=float,
        default=90.0,
    )

    ap.add_argument(
        "--stack-bytes",
        type=int,
        default=1_000_000_000,
    )

    ap.add_argument(
        "--certificate-prime-bound",
        type=int,
        default=3000,
    )

    args = ap.parse_args()

    rows = load_vectors(
        args.vectors,
        args.candidate_pool,
    )

    bases = exact_subgroup_bases(
        rows,
        args.max_coordinate_bits,
    )

    bases = bases[:args.base_count]

    print(
        f"{PROTOCOL}|stage=bases"
        f"|loaded_vectors={len(rows)}"
        f"|usable_bases={len(bases)}"
        f"|requested={args.base_count}",
        flush=True,
    )

    for i, (h, c, P, bits) in enumerate(bases[:10], 1):
        print(
            f"{PROTOCOL}|base={i}"
            f"|height={h:.12g}"
            f"|bits={bits}"
            f"|support={sum(x != 0 for x in c)}"
            f"|l1={sum(abs(x) for x in c)}",
            flush=True,
        )

    plans = []

    for index, (height, coeffs, base, bits) in enumerate(bases, 1):
        plans.extend(
            plans_for_base(
                f"q{index:03d}",
                base,
                x_pair_count=args.x_pair_count,
                slope_anchor_count=args.slope_anchor_count,
                slope_pair_count=args.slope_pair_count,
                x_offset_height=args.x_offset_height,
                x_pair_height=args.x_pair_height,
                slope_offset_height=args.slope_offset_height,
                slope_pair_height=args.slope_pair_height,
            )
        )

    print(
        f"{PROTOCOL}|stage=start"
        f"|bases={len(bases)}"
        f"|charts={len(plans)}",
        flush=True,
    )

    discovered = {}

    # First close the two holes in the previous run.
    for point in retry_old_timeouts(
        timeout=args.retry_timeout,
        stack_bytes=args.stack_bytes,
        x_offset_height=1_000_000,
    ):
        discovered.setdefault(
            canonical_point(point),
            set(),
        ).add("retry")

    started = time.monotonic()
    completed = 0
    timeouts = 0

    for index, plan in enumerate(plans, 1):
        try:
            qpts, ms, wall = search_chart(
                plan.chart,
                timeout=args.chart_timeout,
                stack_bytes=args.stack_bytes,
            )

        except subprocess.TimeoutExpired:
            timeouts += 1
            continue

        completed += 1

        for point in map_images(plan, qpts):
            key = canonical_point(point)

            discovered.setdefault(
                key,
                set(),
            ).add(plan.chart.identifier)

        if index % 100 == 0:
            print(
                f"{PROTOCOL}|progress"
                f"|charts={index}/{len(plans)}"
                f"|completed={completed}"
                f"|timeouts={timeouts}"
                f"|images={len(discovered)}"
                f"|seconds={time.monotonic()-started:.1f}",
                flush=True,
            )

    print(
        f"{PROTOCOL}|stage=classify"
        f"|unique_images={len(discovered)}",
        flush=True,
    )

    hit = False

    # Raw submitted points and their negatives are known subgroup points.
    known = {
        canonical_point(P)
        for P in PUBLISHED_POINTS
    }

    for point, sources in sorted(
        discovered.items(),
        key=lambda item: (item[0][0], item[0][1]),
    ):
        if point in known:
            continue

        (
            classification,
            rank,
            relation,
            primes,
        ) = classify_candidate(
            point,
            relation_timeout=args.relation_timeout,
            stack_bytes=args.stack_bytes,
            certificate_prime_bound=args.certificate_prime_bound,
        )

        print(
            f"{PROTOCOL}|candidate"
            f"|classification={classification}"
            f"|rank={rank}"
            f"|sources={len(sources)}"
            f"|x={point[0]}"
            f"|y={point[1]}"
            f"|primes={','.join(map(str,primes))}",
            flush=True,
        )

        if classification == "EXACT_INDEPENDENT_31ST":
            hit = True

            print(
                f"{PROTOCOL}|RANK31"
                f"|status=EXACT_UNCONDITIONAL"
                f"|x={point[0]}"
                f"|y={point[1]}",
                flush=True,
            )

            break

    print(
        f"{PROTOCOL}|stage=done"
        f"|completed={completed}/{len(plans)}"
        f"|timeouts={timeouts}"
        f"|unique_images={len(discovered)}"
        f"|rank31_hit={str(hit).lower()}"
        f"|seconds={time.monotonic()-started:.3f}",
        flush=True,
    )


if __name__ == "__main__":
    main()

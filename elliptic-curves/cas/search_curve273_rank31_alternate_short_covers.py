#!/usr/bin/env python3

from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
import json
from pathlib import Path
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent))

import search_icarm_curve273_rank31 as engine

Q = engine.Q
QuarticChart = engine.QuarticChart
affine_substitute = engine.affine_substitute
discover_relation = engine.discover_relation
exact_linear_combination = engine.exact_linear_combination
gp_rational = engine.gp_rational
gp_vector = engine.gp_vector
point_negate = engine.point_negate
poly_add = engine.poly_add
poly_evaluate = engine.poly_evaluate
poly_multiply = engine.poly_multiply
poly_scale = engine.poly_scale
search_chart = engine.search_chart
run_gp = engine.run_gp

PUBLISHED_POINTS = ()
short_weierstrass_coefficients = None
to_short_point = None
A1 = A2 = A3 = Q(0)

from mod2_reduction_independence import (
    combined_mod2_rank,
    find_mod2_reduction_certificate,
)


PROTOCOL = "R31ALT"

ROOT = Path(__file__).resolve().parents[2]

DEFAULT_VECTORS = (
    ROOT
    / "artifacts/local/elliptic-curves/curve273-rank30"
    / "short-coefficient-vectors.tsv"
)


RationalPoint = tuple[Fraction, Fraction]
Polynomial = tuple[Fraction, ...]


def configure_curve(curve_id, curve_json=None):
    global PUBLISHED_POINTS, short_weierstrass_coefficients, to_short_point
    global A1, A2, A3, SHORT_BASIS, SIGNED_SHORT_BASIS

    engine.load_curve_data(curve_id, curve_json)
    PUBLISHED_POINTS = engine.PUBLISHED_POINTS
    short_weierstrass_coefficients = engine.short_weierstrass_coefficients
    to_short_point = engine.to_short_point
    A1, A2, A3 = engine.A1, engine.A2, engine.A3
    SHORT_BASIS = tuple(to_short_point(point) for point in PUBLISHED_POINTS)
    SIGNED_SHORT_BASIS = tuple(
        oriented
        for point in SHORT_BASIS
        for oriented in (point, short_negate(point))
    )


@dataclass(frozen=True)
class CoverBase:
    index: int
    predicted_height: float
    coefficients: tuple[int, ...]
    general_point: RationalPoint
    short_point: RationalPoint
    coordinate_bits: int


@dataclass(frozen=True)
class AltPlan:
    chart: QuarticChart
    cover: CoverBase


def rat_bits(value):
    value = Q(value)
    return max(
        abs(value.numerator).bit_length(),
        value.denominator.bit_length(),
    )


def parameter_height(value):
    value = Q(value)
    return max(
        abs(value.numerator),
        value.denominator,
    )


def polynomial_bits(poly):
    return max(rat_bits(v) for v in poly)


def canonical_general_point(point):
    inverse = point_negate(point)
    return min(point, inverse)


# ------------------------------------------------------------
# Short model helpers
#
# X = 36*x + 3*b2
# Y = 108*(2*y+a1*x+a3)
# ------------------------------------------------------------

def from_short_point(point):
    X, Y = map(Q, point)

    b2 = A1*A1 + 4*A2
    x = (X - 3*b2) / 36
    y = Y / 216 - (A1*x + A3) / 2

    return x, y


def short_negate(point):
    return point[0], -point[1]


def point_on_short(point):
    x, y = map(Q, point)

    coeffs = short_weierstrass_coefficients()

    return (
        y*y
        ==
        x**3 + coeffs[3]*x + coeffs[4]
    )


SHORT_BASIS = ()
SIGNED_SHORT_BASIS = ()


# ------------------------------------------------------------
# Alternate degree-two cover
#
# For Q=(Xq,Yq):
#
#     t_Q(P) = (Y(P)+Yq)/(X(P)-Xq)
#
# Lines through -Q give:
#
#     w^2 = D_Q(t)
#
# ------------------------------------------------------------

def alternate_discriminant(short_q):
    xq, yq = map(Q, short_q)

    Ashort = Q(
        short_weierstrass_coefficients()[3]
    )

    t: Polynomial = (Q(0), Q(1))

    # Line:
    #
    #   y = t*x - yq - t*xq
    #
    intercept: Polynomial = (
        -yq,
        -xq,
    )

    t2 = poly_multiply(t, t)

    alpha = poly_add(
        (xq,),
        poly_scale(t2, Q(-1)),
    )

    c1 = poly_add(
        poly_scale(
            poly_multiply(t, intercept),
            Q(2),
        ),
        (-Ashort,),
    )

    D = poly_add(
        poly_add(
            poly_multiply(alpha, alpha),
            poly_scale(
                poly_multiply(alpha, (xq,)),
                Q(-4),
            ),
        ),
        poly_scale(c1, Q(4)),
    )

    if len(D) != 5 or D[-1] != 1:
        raise AssertionError(
            "alternate discriminant is not monic quartic"
        )

    return D


def alternate_parameter(short_q, short_p):
    xq, yq = map(Q, short_q)
    xp, yp = map(Q, short_p)

    if xp == xq:
        raise ValueError("vertical alternate slope")

    return (yp + yq) / (xp - xq)


def alternate_to_short_points(
    short_q,
    parameter,
    square_root,
):
    xq, yq = map(Q, short_q)

    parameter = Q(parameter)
    square_root = Q(square_root)

    alpha = xq - parameter*parameter

    x1 = (-alpha + square_root) / 2
    x2 = (-alpha - square_root) / 2

    intercept = (
        -yq
        - parameter*xq
    )

    points = (
        (
            x1,
            parameter*x1 + intercept,
        ),
        (
            x2,
            parameter*x2 + intercept,
        ),
    )

    if not all(
        point_on_short(P)
        for P in points
    ):
        raise AssertionError(
            "alternate cover mapped off short curve"
        )

    return points


# ------------------------------------------------------------
# Load short MW vectors produced in previous stage
# ------------------------------------------------------------

def load_vectors(path, count):
    rows = []

    lines = path.read_text().splitlines()

    for line in lines[1:]:
        if not line.strip():
            continue

        parts = line.split("\t")

        height = float(parts[0])
        coeffs = tuple(
            int(x)
            for x in parts[1:]
        )

        if len(coeffs) != len(PUBLISHED_POINTS):
            raise ValueError(
                f"short-vector width != {len(PUBLISHED_POINTS)}"
            )

        rows.append(
            (height, coeffs)
        )

        if len(rows) >= count:
            break

    return rows


def build_cover_bases(
    rows,
    *,
    count,
    max_coordinate_bits,
):
    result = []

    seen = set()

    for height, coeffs in rows:
        P = exact_linear_combination(coeffs)

        if P is None:
            raise AssertionError(
                "nonzero coefficient vector vanished"
            )

        key = canonical_general_point(P)

        if key in seen:
            continue

        seen.add(key)

        short = to_short_point(P)

        if not point_on_short(short):
            raise AssertionError(
                "short transport failed"
            )

        bits = max(
            rat_bits(P[0]),
            rat_bits(P[1]),
        )

        if bits > max_coordinate_bits:
            continue

        result.append(
            CoverBase(
                index=len(result) + 1,
                predicted_height=height,
                coefficients=coeffs,
                general_point=P,
                short_point=short,
                coordinate_bits=bits,
            )
        )

        if len(result) >= count:
            break

    return result


# ------------------------------------------------------------
# Parameter compression
# ------------------------------------------------------------

def signed_public_parameters(base):
    values = []

    for index, P in enumerate(
        SIGNED_SHORT_BASIS
    ):
        try:
            t = alternate_parameter(
                base.short_point,
                P,
            )
        except ValueError:
            continue

        values.append(
            (index, t)
        )

    return values


def distinct_parameter_values(parameters):
    seen = set()
    answer = []

    for index, value in sorted(
        parameters,
        key=lambda item: (
            parameter_height(item[1]),
            item[0],
        ),
    ):
        if value in seen:
            continue

        seen.add(value)
        answer.append(
            (index, value)
        )

    return answer


def compression_signature(values):
    heights = sorted(
        rat_bits(v)
        for v in values
    )

    if not heights:
        return (10**9,)

    return (
        heights[
            min(4, len(heights)-1)
        ],
        heights[len(heights)//4],
        heights[len(heights)//2],
        heights[-1],
        sum(heights),
    )


def build_plans_for_cover(
    base,
    *,
    offset_count,
    normalization_pool,
    affine_count,
    offset_height,
    affine_height,
):
    D = alternate_discriminant(
        base.short_point
    )

    parameters = signed_public_parameters(
        base
    )

    distinct = distinct_parameter_values(
        parameters
    )

    if len(distinct) < 4:
        return []

    plans = []

    # --------------------------------------------------------
    # Offset charts around smallest known alternate parameters
    # --------------------------------------------------------

    for offset_index, (
        source_index,
        center,
    ) in enumerate(
        distinct[:offset_count],
        1,
    ):
        poly = affine_substitute(
            D,
            center,
            Q(1),
        )

        plans.append(
            AltPlan(
                QuarticChart(
                    identifier=(
                        f"q{base.index:03d}"
                        f"_alt_offset"
                        f"{offset_index:02d}"
                    ),
                    kind="alternate_short_offset",
                    polynomial=poly,
                    center=center,
                    scale=Q(1),
                    base_index=None,
                    seed_parameters=(Q(0),),
                    height_bound=offset_height,
                ),
                base,
            )
        )

    # --------------------------------------------------------
    # Affine normalizations t = center + scale*u
    #
    # Normalize two known alternate parameters to 0 and 1.
    # --------------------------------------------------------

    pool = distinct[
        :min(
            normalization_pool,
            len(distinct),
        )
    ]

    ranked = []

    all_values = tuple(
        value
        for _, value in distinct
    )

    for (
        (index1, t1),
        (index2, t2),
    ) in combinations(pool, 2):

        scale = t2 - t1

        if scale == 0:
            continue

        normalized = tuple(
            (value - t1) / scale
            for value in all_values
        )

        poly = affine_substitute(
            D,
            t1,
            scale,
        )

        score = (
            compression_signature(
                normalized
            ),
            polynomial_bits(poly),
            index1,
            index2,
        )

        ranked.append(
            (
                score,
                t1,
                scale,
                poly,
                index1,
                index2,
            )
        )

    ranked.sort(
        key=lambda item: item[0]
    )

    for rank, (
        score,
        center,
        scale,
        poly,
        index1,
        index2,
    ) in enumerate(
        ranked[:affine_count],
        1,
    ):

        plans.append(
            AltPlan(
                QuarticChart(
                    identifier=(
                        f"q{base.index:03d}"
                        f"_alt_affine"
                        f"{rank:02d}"
                    ),
                    kind="alternate_short_affine",
                    polynomial=poly,
                    center=center,
                    scale=scale,
                    base_index=None,
                    seed_parameters=(
                        Q(0),
                        Q(1),
                    ),
                    height_bound=affine_height,
                ),
                base,
            )
        )

    return plans


# ------------------------------------------------------------
# Map quartic images back to curve 273
# ------------------------------------------------------------

def map_plan_images(
    plan,
    quartic_points,
):
    chart = plan.chart

    images = []

    for local_t, square_root in quartic_points:

        if (
            poly_evaluate(
                chart.polynomial,
                local_t,
            )
            !=
            square_root*square_root
        ):
            raise AssertionError(
                "PARI point missed alternate quartic"
            )

        if local_t in chart.seed_parameters:
            continue

        parameter = (
            chart.center
            + chart.scale*local_t
        )

        short_points = (
            alternate_to_short_points(
                plan.cover.short_point,
                parameter,
                square_root,
            )
        )

        for short_point in short_points:
            general = from_short_point(
                short_point
            )

            # Verify by transporting back again.
            replay = to_short_point(
                general
            )

            if replay != short_point:
                raise AssertionError(
                    "short/general transport mismatch"
                )

            images.append(general)

    return tuple(images)


# ------------------------------------------------------------
# Rank-31 classification
# ------------------------------------------------------------

def certificate_rank(
    point,
    prime_bound,
):
    augmented = (
        tuple(
            to_short_point(P)
            for P in PUBLISHED_POINTS
        )
        +
        (
            to_short_point(point),
        )
    )

    signatures = (
        find_mod2_reduction_certificate(
            short_weierstrass_coefficients(),
            augmented,
            prime_bound=prime_bound,
        )
    )

    rank = combined_mod2_rank(
        signatures,
        len(augmented),
    )

    return (
        rank,
        tuple(
            s.prime
            for s in signatures
        ),
    )


def classify(
    point,
    *,
    prime_bound,
    relation_timeout,
    stack_bytes,
):
    # First try the cheap exact rank-31 proof.
    rank, primes = certificate_rank(
        point,
        prime_bound,
    )

    if rank == len(PUBLISHED_POINTS) + 1:
        return (
            "EXACT_INDEPENDENT_NEXT_POINT",
            rank,
            primes,
        )

    # If the mod-2 image lies in the known span, try to
    # prove that the point really belongs to the subgroup.
    relation = discover_relation(
        point,
        timeout=relation_timeout,
        stack_bytes=stack_bytes,
    )

    if relation is not None:
        return (
            "certified_subgroup",
            rank,
            primes,
        )

    return (
        "UNRESOLVED_MOD2_DEPENDENT",
        rank,
        primes,
    )


def discover_relations_batch(
    points,
    *,
    timeout,
    stack_bytes,
):
    """Recover many subgroup relations with one height-matrix initialization."""

    if not points:
        return {}
    curve = ",".join(gp_rational(value) for value in engine.GENERAL_WEIERSTRASS_COEFFICIENTS)
    basis = ",".join(gp_vector(value) for value in PUBLISHED_POINTS)
    commands = [
        "default(realprecision,140);",
        f"E=ellinit([{curve}]);",
        f"B=[{basis}];",
        "H=ellheightmatrix(E,B);",
    ]
    for index, point in enumerate(points):
        commands.extend(
            [
                f"Q={gp_vector(point)};",
                "V=vector(#B,j,ellheight(E,B[j],Q))~;",
                "C=round(matsolve(H,V));",
                "S=[0];for(j=1,#B,S=elladd(E,S,ellmul(E,B[j],C[j])));",
                f'print("BATCHREL|{index}|",Vec(C),"|",S==Q);',
            ]
        )
    commands.append("quit")
    try:
        output, _wall = run_gp(
            "\n".join(commands) + "\n",
            timeout=timeout,
            stack_bytes=stack_bytes,
        )
    except (RuntimeError, subprocess.TimeoutExpired):
        return {}

    relations = {}
    for line in output.splitlines():
        if not line.startswith("BATCHREL|"):
            continue
        _tag, raw_index, raw_vector, exact = line.split("|", 3)
        if exact != "1":
            continue
        coefficients = tuple(
            int(value.strip())
            for value in raw_vector.strip()[1:-1].split(",")
        )
        index = int(raw_index)
        if exact_linear_combination(coefficients) != points[index]:
            raise AssertionError("batched PARI relation failed exact replay")
        relations[index] = coefficients
    return relations


def main():
    ap = argparse.ArgumentParser()

    ap.add_argument("--curve-id", type=int, default=273)
    ap.add_argument("--curve-json", type=Path)
    ap.add_argument("--output", type=Path)

    ap.add_argument(
        "--vectors",
        type=Path,
        default=None,
    )

    ap.add_argument(
        "--vector-pool",
        type=int,
        default=600,
    )

    ap.add_argument(
        "--cover-count",
        type=int,
        default=120,
    )

    ap.add_argument(
        "--max-coordinate-bits",
        type=int,
        default=1400,
    )

    ap.add_argument(
        "--offset-count",
        type=int,
        default=4,
    )

    ap.add_argument(
        "--normalization-pool",
        type=int,
        default=10,
    )

    ap.add_argument(
        "--affine-count",
        type=int,
        default=8,
    )

    ap.add_argument(
        "--offset-height",
        type=int,
        default=200000,
    )

    ap.add_argument(
        "--affine-height",
        type=int,
        default=100000,
    )

    ap.add_argument(
        "--chart-timeout",
        type=float,
        default=6.0,
    )

    ap.add_argument(
        "--relation-timeout",
        type=float,
        default=30.0,
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

    configure_curve(args.curve_id, args.curve_json)
    if args.vectors is None:
        if args.curve_id == 273:
            args.vectors = DEFAULT_VECTORS
        elif args.curve_id == 245:
            args.vectors = (
                ROOT
                / "artifacts/local/elliptic-curves/curve245-rank20/"
                "short-coefficient-vectors.tsv"
            )
        else:
            raise SystemExit("--vectors is required for this ICARM curve id")

    rows = load_vectors(
        args.vectors,
        args.vector_pool,
    )

    covers = build_cover_bases(
        rows,
        count=args.cover_count,
        max_coordinate_bits=args.max_coordinate_bits,
    )

    print(
        f"{PROTOCOL}|stage=bases"
        f"|vectors={len(rows)}"
        f"|covers={len(covers)}",
        flush=True,
    )

    for cover in covers[:10]:
        print(
            f"{PROTOCOL}|base={cover.index}"
            f"|height={cover.predicted_height:.12g}"
            f"|bits={cover.coordinate_bits}"
            f"|support="
            f"{sum(c != 0 for c in cover.coefficients)}"
            f"|l1="
            f"{sum(abs(c) for c in cover.coefficients)}",
            flush=True,
        )

    plans = []

    for cover in covers:
        plans.extend(
            build_plans_for_cover(
                cover,
                offset_count=args.offset_count,
                normalization_pool=args.normalization_pool,
                affine_count=args.affine_count,
                offset_height=args.offset_height,
                affine_height=args.affine_height,
            )
        )

    print(
        f"{PROTOCOL}|stage=start"
        f"|covers={len(covers)}"
        f"|charts={len(plans)}"
        f"|offset_height={args.offset_height}"
        f"|affine_height={args.affine_height}",
        flush=True,
    )

    known_basis = {
        canonical_general_point(P)
        for P in PUBLISHED_POINTS
    }

    discovered = {}

    completed = 0
    timeouts = 0

    total_quartic_points = 0
    total_images = 0
    known_basis_images = 0

    started = time.monotonic()

    for index, plan in enumerate(
        plans,
        1,
    ):

        try:
            (
                quartic_points,
                milliseconds,
                wall,
            ) = search_chart(
                plan.chart,
                timeout=args.chart_timeout,
                stack_bytes=args.stack_bytes,
            )

        except subprocess.TimeoutExpired:
            timeouts += 1

            print(
                f"{PROTOCOL}|timeout"
                f"|chart={plan.chart.identifier}",
                flush=True,
            )

            continue

        completed += 1

        total_quartic_points += len(
            quartic_points
        )

        images = map_plan_images(
            plan,
            quartic_points,
        )

        total_images += len(images)

        for point in images:
            key = canonical_general_point(
                point
            )

            if key in known_basis:
                known_basis_images += 1
                continue

            discovered.setdefault(
                key,
                set(),
            ).add(
                plan.chart.identifier
            )

        if (
            index % 100 == 0
            or images
        ):
            print(
                f"{PROTOCOL}|progress"
                f"|charts={index}/{len(plans)}"
                f"|completed={completed}"
                f"|timeouts={timeouts}"
                f"|quartic_points={total_quartic_points}"
                f"|images={total_images}"
                f"|known_basis_images={known_basis_images}"
                f"|unique_nonbasis={len(discovered)}"
                f"|seconds={time.monotonic()-started:.1f}",
                flush=True,
            )

    print(
        f"{PROTOCOL}|stage=classify"
        f"|unique_nonbasis={len(discovered)}"
        f"|total_images={total_images}",
        flush=True,
    )

    rank31_hit = False
    subgroup_count = 0
    unresolved_count = 0
    ordered = sorted(
        discovered.items(),
        key=lambda item: (
            item[0][0],
            item[0][1],
        ),
    )
    all_points = [point for point, _sources in ordered]
    batched = discover_relations_batch(
        all_points,
        timeout=max(60.0, args.relation_timeout),
        stack_bytes=args.stack_bytes,
    )
    candidate_rows = []
    for index, (point, sources) in enumerate(ordered):
        relation = batched.get(index)
        if relation is not None:
            classification = "certified_subgroup"
            rank = len(PUBLISHED_POINTS)
            primes = ()
        else:
            # Only points not replayed in the known integer span pay for the
            # exhaustive finite-reduction calculation.
            relation = discover_relation(
                point,
                timeout=args.relation_timeout,
                stack_bytes=args.stack_bytes,
            )
            if relation is not None:
                classification = "certified_subgroup"
                rank = len(PUBLISHED_POINTS)
                primes = ()
            else:
                rank, primes = certificate_rank(point, args.certificate_prime_bound)
                if rank == len(PUBLISHED_POINTS) + 1:
                    classification = "EXACT_INDEPENDENT_NEXT_POINT"
                    rank31_hit = True
                else:
                    classification = "UNRESOLVED_MOD2_DEPENDENT"
        candidate_rows.append(
            {
                "point": point,
                "sources": sources,
                "classification": classification,
                "rank": rank,
                "primes": primes,
                "relation": relation,
            }
        )

    for row in candidate_rows:
        point = row["point"]
        sources = row["sources"]
        classification = row["classification"]
        rank = row["rank"]
        primes = row["primes"]

        if classification == "certified_subgroup":
            subgroup_count += 1

        if classification == "UNRESOLVED_MOD2_DEPENDENT":
            unresolved_count += 1

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

        if classification == "EXACT_INDEPENDENT_NEXT_POINT":
            print(
                f"{PROTOCOL}|NEXT_POINT"
                f"|status=EXACT_UNCONDITIONAL"
                f"|x={point[0]}"
                f"|y={point[1]}",
                flush=True,
            )

    elapsed = time.monotonic() - started

    print(
        f"{PROTOCOL}|stage=done"
        f"|completed={completed}/{len(plans)}"
        f"|timeouts={timeouts}"
        f"|quartic_points={total_quartic_points}"
        f"|total_images={total_images}"
        f"|known_basis_images={known_basis_images}"
        f"|unique_nonbasis={len(discovered)}"
        f"|subgroup={subgroup_count}"
        f"|unresolved={unresolved_count}"
        f"|rank31_hit={str(rank31_hit).lower()}"
        f"|seconds={elapsed:.3f}",
        flush=True,
    )

    if args.output is not None:
        artifact = {
            "schema_version": 1,
            "artifact_kind": "bounded_alternate_degree_two_cover_search",
            "curve_id": args.curve_id,
            "basis_rank": len(PUBLISHED_POINTS),
            "claim_scope": (
                "Exact curve maps, memberships, subgroup relations, and any "
                "finite-reduction independence certificate; bounded chart boxes only."
            ),
            "parameters": {
                "vector_pool": args.vector_pool,
                "cover_count": args.cover_count,
                "offset_count": args.offset_count,
                "normalization_pool": args.normalization_pool,
                "affine_count": args.affine_count,
                "offset_height": args.offset_height,
                "affine_height": args.affine_height,
                "chart_timeout": args.chart_timeout,
                "relation_timeout": args.relation_timeout,
                "certificate_prime_bound": args.certificate_prime_bound,
            },
            "summary": {
                "completed": completed,
                "chart_count": len(plans),
                "timeouts": timeouts,
                "quartic_points": total_quartic_points,
                "total_images": total_images,
                "known_basis_images": known_basis_images,
                "unique_nonbasis": len(discovered),
                "certified_subgroup": subgroup_count,
                "unresolved": unresolved_count,
                "independent_next_point_hit": rank31_hit,
                "wall_seconds": elapsed,
            },
            "candidates": [
                {
                    "x": str(row["point"][0]),
                    "y": str(row["point"][1]),
                    "sources": sorted(row["sources"]),
                    "classification": row["classification"],
                    "augmented_mod2_rank": row["rank"],
                    "certificate_primes": list(row["primes"]),
                    "basis_relation": (
                        None
                        if row["relation"] is None
                        else list(row["relation"])
                    ),
                }
                for row in candidate_rows
            ],
        }
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
        print(f"{PROTOCOL}|saved={args.output}", flush=True)


if __name__ == "__main__":
    main()

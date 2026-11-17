#!/usr/bin/env python3
"""Exact bounded search for failures of five-to-four-channel descent.

Fix the two-endpoint operator

    A(X,Y) = X^r - r!/s! Y^s

and normalize the endpoint coefficients of P at (r,0) and (0,s) to one.
For every choice of at most three additional support points in the
requested box, this script constructs the scalar moments

    M_m = [A(partial)^m(P^m)]_(0,0)

over QQ.  It saturates the ideal (M_2,...,M_N) by the product of all
additional coefficients, so a reported survivor has exactly the stated
channels nonzero.  When the box multiple is at least N, points outside the
box cannot contribute through order N; the search then covers arbitrary
nonnegative support with at most five channels for this finite moment
truncation.

This is an exact bounded experiment, not an all-support theorem and not a
search for the full polynomial identities required by GVC.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from itertools import combinations, permutations
from math import factorial
from pathlib import Path

import sympy as sp


def compositions(
    total: int, parts: int, prefix: tuple[int, ...] = ()
):
    if parts == 1:
        yield prefix + (total,)
        return
    for value in range(total + 1):
        yield from compositions(total - value, parts - 1, prefix + (value,))


def scalar_moments(
    r: int,
    s: int,
    support: tuple[tuple[int, int], ...],
    coefficients: tuple[sp.Symbol, ...],
    order: int,
) -> list[sp.Expr]:
    q = sp.Rational(-factorial(r), factorial(s))
    answer: list[sp.Expr] = []
    for moment_order in range(2, order + 1):
        moment = sp.Integer(0)
        for counts in compositions(moment_order, len(support)):
            x_exp = sum(
                counts[index] * support[index][0]
                for index in range(len(support))
            )
            y_exp = sum(
                counts[index] * support[index][1]
                for index in range(len(support))
            )
            if x_exp % r or y_exp % s:
                continue
            operator_y_count = y_exp // s
            if not (
                0 <= operator_y_count <= moment_order
                and x_exp == r * (moment_order - operator_y_count)
            ):
                continue

            multinomial = factorial(moment_order)
            for count in counts:
                multinomial //= factorial(count)
            coefficient = (
                multinomial
                * sp.binomial(moment_order, operator_y_count)
                * q**operator_y_count
                * factorial(x_exp)
                * factorial(y_exp)
            )
            for variable, power in zip(coefficients, counts[2:]):
                coefficient *= variable**power
            moment += coefficient
        answer.append(sp.factor(moment))
    return answer


def balanced_selections(
    r: int,
    s: int,
    support: tuple[tuple[int, int], ...],
    moment_order: int,
) -> tuple[tuple[int, ...], ...]:
    """Return ``(operator-y count, polynomial counts...)`` balance rows."""

    rows: list[tuple[int, ...]] = []
    for counts in compositions(moment_order, len(support)):
        x_exp = sum(
            counts[index] * support[index][0]
            for index in range(len(support))
        )
        y_exp = sum(
            counts[index] * support[index][1]
            for index in range(len(support))
        )
        if x_exp % r or y_exp % s:
            continue
        operator_y_count = y_exp // s
        if (
            0 <= operator_y_count <= moment_order
            and x_exp == r * (moment_order - operator_y_count)
        ):
            rows.append((operator_y_count, *counts))
    return tuple(sorted(rows))


def canonical_return_signature(
    r: int,
    s: int,
    support: tuple[tuple[int, int], ...],
    order: int,
) -> list[list[list[int]]]:
    """Canonicalize balance rows under endpoint swap and extra permutations."""

    extra_count = len(support) - 2
    raw = {
        moment_order: balanced_selections(r, s, support, moment_order)
        for moment_order in range(2, order + 1)
    }
    candidates: list[tuple[tuple[tuple[int, ...], ...], ...]] = []
    for permutation in permutations(range(extra_count)):
        for swap_endpoints in (False, True):
            levels: list[tuple[tuple[int, ...], ...]] = []
            for moment_order in range(2, order + 1):
                transformed: list[tuple[int, ...]] = []
                for row in raw[moment_order]:
                    operator_y_count, endpoint_a, endpoint_b, *extra = row
                    if swap_endpoints:
                        prefix = (
                            moment_order - operator_y_count,
                            endpoint_b,
                            endpoint_a,
                        )
                    else:
                        prefix = (
                            operator_y_count,
                            endpoint_a,
                            endpoint_b,
                        )
                    transformed.append(
                        (*prefix, *(extra[index] for index in permutation))
                    )
                levels.append(tuple(sorted(transformed)))
            candidates.append(tuple(levels))
    canonical = min(candidates)
    return [
        [list(row) for row in level]
        for level in canonical
    ]


def is_unit_groebner(basis: sp.GroebnerBasis) -> bool:
    return (
        len(basis.polys) == 1
        and basis.polys[0].as_expr() == 1
    )


def search(r: int, s: int, order: int, box_multiple: int) -> dict[str, object]:
    if r == s:
        raise ValueError("choose unequal endpoint orders")
    if order < 2:
        raise ValueError("order must be at least two")
    if box_multiple < 1:
        raise ValueError("box multiple must be positive")

    coefficient_pool = sp.symbols("c0 c1 c2")
    inverse = sp.symbols("inverse")
    endpoint_a = (r, 0)
    endpoint_b = (0, s)
    points = [
        (x_exp, y_exp)
        for x_exp in range(box_multiple * r + 1)
        for y_exp in range(box_multiple * s + 1)
        if (x_exp, y_exp) not in (endpoint_a, endpoint_b)
    ]

    tested_by_extra_count: dict[int, int] = {}
    first_contradiction = {
        moment_order: 0 for moment_order in range(2, order + 1)
    }
    pivot_examples: dict[
        int, list[tuple[tuple[tuple[int, int], ...], list[sp.Expr]]]
    ] = {
        moment_order: [] for moment_order in range(2, order + 1)
    }
    survivors: list[
        tuple[tuple[tuple[int, int], ...], list[sp.Expr]]
    ] = []
    final_pivots: list[dict[str, object]] = []
    for extra_count in range(4):
        tested_by_extra_count[extra_count] = 0
        coefficients = coefficient_pool[:extra_count]
        for extra_support in combinations(points, extra_count):
            tested_by_extra_count[extra_count] += 1
            support = (endpoint_a, endpoint_b) + extra_support
            moments = scalar_moments(r, s, support, coefficients, order)

            contradicted = False
            for index, moment_order in enumerate(range(2, order + 1)):
                prefix = moments[: index + 1]
                if any(
                    moment.is_number and moment != 0 for moment in prefix
                ):
                    is_unit = True
                elif coefficients:
                    coefficient_product = sp.prod(coefficients)
                    saturated_basis = sp.groebner(
                        prefix + [inverse * coefficient_product - 1],
                        inverse,
                        *coefficients,
                        order="grevlex",
                    )
                    is_unit = is_unit_groebner(saturated_basis)
                else:
                    is_unit = False
                if is_unit:
                    first_contradiction[moment_order] += 1
                    if len(pivot_examples[moment_order]) < 3:
                        pivot_examples[moment_order].append(
                            (extra_support, moments)
                        )
                    if moment_order == order:
                        final_pivots.append(
                            {
                                "endpoint_pair": [r, s],
                                "extra_support": [
                                    list(point) for point in extra_support
                                ],
                                "extra_count": extra_count,
                                "moments": [str(moment) for moment in moments],
                                "return_signature": canonical_return_signature(
                                    r,
                                    s,
                                    support,
                                    order,
                                ),
                            }
                        )
                    contradicted = True
                    break
            if not contradicted:
                survivors.append((extra_support, moments))

    print(
        f"r={r}, s={s}, moments=1..{order}, "
        f"box=[0,{box_multiple * r}]x[0,{box_multiple * s}]"
    )
    print(f"supports tested by added-channel count: {tested_by_extra_count}")
    print(f"supports tested in total: {sum(tested_by_extra_count.values())}")
    print(f"first contradiction by moment: {first_contradiction}")
    print(f"nonzero-channel torus survivors: {len(survivors)}")
    if pivot_examples.get(order):
        print(f"sample supports first killed at M_{order}:")
        for support, moments in pivot_examples[order]:
            print(f"  extra support {support}")
            for moment_order, moment in enumerate(moments, 2):
                print(f"    M_{moment_order} = {moment}")
    for support, moments in survivors[:10]:
        print(f"  extra support {support}")
        for moment_order, moment in enumerate(moments, 2):
            print(f"    M_{moment_order} = {moment}")
    if box_multiple >= order:
        print(
            "coverage: all nonnegative supports with at most five channels "
            f"for moments 1..{order}"
        )
    else:
        print("coverage: bounded support box only")
    print("STATUS: exact finite-moment computation; no full GVC claim")
    return {
        "endpoint_pair": [r, s],
        "order": order,
        "box_multiple": box_multiple,
        "tested_by_extra_count": tested_by_extra_count,
        "tested_total": sum(tested_by_extra_count.values()),
        "first_contradiction": first_contradiction,
        "survivor_count": len(survivors),
        "final_pivots": final_pivots,
    }


def cluster_final_pivots(
    results: list[dict[str, object]],
) -> list[dict[str, object]]:
    clusters: dict[str, dict[str, object]] = {}
    for result in results:
        for pivot in result["final_pivots"]:
            signature = pivot["return_signature"]
            encoded = json.dumps(signature, separators=(",", ":"))
            cluster = clusters.setdefault(
                encoded,
                {
                    "return_signature": signature,
                    "count": 0,
                    "endpoint_pair_counts": {},
                    "representatives": [],
                },
            )
            cluster["count"] += 1
            pair_key = ",".join(str(value) for value in pivot["endpoint_pair"])
            pair_counts = cluster["endpoint_pair_counts"]
            pair_counts[pair_key] = pair_counts.get(pair_key, 0) + 1
            if len(cluster["representatives"]) < 4:
                cluster["representatives"].append(
                    {
                        "endpoint_pair": pivot["endpoint_pair"],
                        "extra_support": pivot["extra_support"],
                        "moments": pivot["moments"],
                    }
                )

    answer: list[dict[str, object]] = []
    for encoded, cluster in clusters.items():
        cluster["id"] = "R-" + hashlib.sha256(
            encoded.encode("utf-8")
        ).hexdigest()[:12]
        cluster["endpoint_pair_counts"] = dict(
            sorted(cluster["endpoint_pair_counts"].items())
        )
        answer.append(cluster)
    answer.sort(key=lambda item: (-item["count"], item["id"]))
    return answer


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--r", type=int)
    parser.add_argument("--s", type=int)
    parser.add_argument(
        "--frontier-suite",
        action="store_true",
        help="run (1,2), (1,3), (1,4), and (2,3)",
    )
    parser.add_argument("--order", type=int, default=4)
    parser.add_argument("--box-multiple", type=int, default=4)
    parser.add_argument(
        "--json-output",
        type=Path,
        help="write exact pivot records and canonical return clusters",
    )
    return parser.parse_args()


def main() -> None:
    arguments = parse_args()
    if arguments.r is None and arguments.s is None:
        endpoint_pairs = (
            ((1, 2), (1, 3), (1, 4), (2, 3))
            if arguments.frontier_suite
            else ((1, 2), (1, 3))
        )
    elif arguments.r is None or arguments.s is None:
        raise ValueError("--r and --s must be supplied together")
    elif arguments.frontier_suite:
        raise ValueError("--frontier-suite cannot be combined with --r and --s")
    else:
        endpoint_pairs = ((arguments.r, arguments.s),)
    results = [
        search(
            r=r,
            s=s,
            order=arguments.order,
            box_multiple=arguments.box_multiple,
        )
        for r, s in endpoint_pairs
    ]
    clusters = cluster_final_pivots(results)
    print(f"canonical final-pivot return clusters: {len(clusters)}")
    for cluster in clusters[:12]:
        print(
            f"  {cluster['id']}: count={cluster['count']}, "
            f"pairs={cluster['endpoint_pair_counts']}"
        )

    if arguments.json_output is not None:
        payload = {
            "schema_version": 1,
            "status": (
                "exact finite-moment computation; canonical return "
                "clustering is not an all-support GVC theorem"
            ),
            "command": " ".join(
                ["python", "scripts/search_binary_gvc_five_channel_descent.py"]
                + sys.argv[1:]
            ),
            "software": {
                "python": sys.version.split()[0],
                "sympy": sp.__version__,
                "coefficient_domain": "QQ",
            },
            "results": results,
            "cluster_count": len(clusters),
            "clusters": clusters,
        }
        arguments.json_output.parent.mkdir(parents=True, exist_ok=True)
        arguments.json_output.write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        print(f"wrote {arguments.json_output}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Exact bounded LNED search with positive-dimensional plinth residuals.

This continues ``search_lnd_nonprincipal_plinth.py`` for

    D = u*d/dx + v*d/dy,       w = u*y-v*x,

but replaces zero-dimensional quotients by ideals whose quotient retains a
free y-direction.  The images of u, v, and w are nilpotent in every tested
quotient, so the image of Q[u,v,w] still has an exact finite vector-space
basis; no y-degree truncation is used in a membership decision.

Pure powers through six and mixed powers four through six remain a bounded
counterexample search.
"""

from __future__ import annotations

from itertools import product

import sympy as sp

import search_lnd_nonprincipal_plinth as npl


u, v, x, y = npl.u, npl.v, npl.x, npl.y
invariant = npl.invariant

CHARTS = (
    npl.PlinthIdeal("free-line", (u, v, x)),
    npl.PlinthIdeal(
        "first-plinth-jet",
        (u**2, u * v, v**2, x),
    ),
    npl.PlinthIdeal(
        "split-plinth-jet",
        (u**2, v**2, x),
    ),
    npl.PlinthIdeal(
        "tilted-first-jet",
        (u**2, u * v, v**2, x - u * y),
    ),
    npl.PlinthIdeal(
        "tilted-split-jet",
        (u**2, v**2, x - u * y),
    ),
)


def exponent_tuples(total_degree: int) -> tuple[tuple[int, ...], ...]:
    return tuple(
        exponents
        for exponents in product(
            range(total_degree + 1), repeat=len(npl.variables)
        )
        if sum(exponents) == total_degree
    )


def coefficient_vector(
    f: sp.Expr, exponents: tuple[tuple[int, ...], ...]
) -> tuple[sp.Rational, ...]:
    poly = sp.Poly(sp.expand(f), *npl.variables, domain=sp.QQ)
    return tuple(poly.coeff_monomial(term) for term in exponents)


def verify_free_line_characterization(maximum_degree: int = 8) -> None:
    """Check D(I_n)=im(D_n) intersect ker[v*y^(n-1)] exactly."""
    for total_degree in range(maximum_degree + 1):
        exponents = exponent_tuples(total_degree)
        monomials = tuple(
            sp.prod(
                variable**power
                for variable, power in zip(
                    npl.variables, powers, strict=True
                )
            )
            for powers in exponents
        )
        image_matrix = sp.Matrix.hstack(
            *(
                sp.Matrix(coefficient_vector(npl.derivation(m), exponents))
                for m in monomials
            )
        )
        ideal_columns = [
            sp.Matrix(coefficient_vector(npl.derivation(m), exponents))
            for m, powers in zip(monomials, exponents, strict=True)
            if sum(powers[:3]) > 0
        ]
        ideal_image_matrix = (
            sp.Matrix.hstack(*ideal_columns)
            if ideal_columns
            else sp.zeros(len(exponents), 0)
        )
        if total_degree == 0:
            functional = sp.zeros(1, len(exponents))
        else:
            detected = (0, 1, 0, total_degree - 1)
            functional = sp.zeros(1, len(exponents))
            functional[0, exponents.index(detected)] = 1

        assert functional * ideal_image_matrix == sp.zeros(
            1, ideal_image_matrix.cols
        )
        restricted_rank = (functional * image_matrix).rank()
        assert ideal_image_matrix.rank() == (
            image_matrix.rank() - restricted_rank
        )


def verify_plinth_filtration(maximum_degree: int = 7) -> None:
    """Check im(D) intersect (u,v)^N = D((u,v)^(N-1))."""
    for total_degree in range(maximum_degree + 1):
        exponents = exponent_tuples(total_degree)
        monomials = tuple(
            sp.prod(
                variable**power
                for variable, power in zip(
                    npl.variables, powers, strict=True
                )
            )
            for powers in exponents
        )
        image_columns = tuple(
            sp.Matrix(coefficient_vector(npl.derivation(m), exponents))
            for m in monomials
        )
        image_matrix = sp.Matrix.hstack(*image_columns)
        image_rank = image_matrix.rank()
        for plinth_order in range(1, total_degree + 2):
            filtered_columns = [
                column
                for column, powers in zip(
                    image_columns, exponents, strict=True
                )
                if sum(powers[:2]) >= plinth_order - 1
            ]
            filtered_image = (
                sp.Matrix.hstack(*filtered_columns)
                if filtered_columns
                else sp.zeros(len(exponents), 0)
            )
            low_rows = [
                row
                for row, powers in enumerate(exponents)
                if sum(powers[:2]) < plinth_order
            ]
            low_projection = image_matrix[low_rows, :]
            intersection_dimension = (
                image_rank - low_projection.rank()
            )
            assert filtered_image[low_rows, :] == sp.zeros(
                len(low_rows), filtered_image.cols
            )
            assert filtered_image.rank() == intersection_dimension


def main() -> None:
    pure_bound = 6
    mixed_start = 4
    seeds = npl.seed_family()
    multipliers = (sp.Integer(1), u, v, x, y, invariant)

    pure_queries = tuple(
        sp.expand(seed**exponent)
        for seed in seeds
        for exponent in range(1, pure_bound + 1)
    )
    pure_primitives = npl.particular_primitives(pure_queries)
    verify_free_line_characterization()
    verify_plinth_filtration()

    print(
        "SEARCH:",
        "D=u*d_x+v*d_y;",
        f"{len(seeds)} seeds on positive-dimensional residuals",
    )
    print(
        "THEOREM-REGRESSION:",
        "D((u,v,x))=im(D) intersect ker([v^1]h(0,v,0,y))",
        "through total degree 8",
    )
    print(
        "THEOREM-REGRESSION:",
        "im(D) intersect (u,v)^N=D((u,v)^(N-1))",
        "through total degree 7",
    )
    total_survivors = 0
    total_tail_obstructions = 0
    distinct_survivors: dict[str, sp.Expr] = {}
    chart_survivor_counts: list[int] = []
    chart_escape_prefixes: list[int] = []

    for chart in CHARTS:
        quotient_data = npl.quotient_kernel_span(chart)
        groebner, kernel_basis = quotient_data
        assert not groebner.is_zero_dimensional

        membership_rows = [
            tuple(
                npl.primitive_enters_ideal(
                    pure_primitives[
                        seed_index * pure_bound + exponent - 1
                    ],
                    quotient_data,
                )
                for exponent in range(1, pure_bound + 1)
            )
            for seed_index in range(len(seeds))
        ]
        survivors = [
            seed
            for seed_index, seed in enumerate(seeds)
            if all(membership_rows[seed_index])
        ]
        for seed in survivors:
            distinct_survivors[sp.srepr(seed)] = seed

        escape_prefixes = [
            (
                next(
                    (
                        exponent
                        for exponent, enters in enumerate(row)
                        if not enters
                    ),
                    pure_bound,
                ),
                seed,
            )
            for seed, row in zip(seeds, membership_rows, strict=True)
            if not npl.in_base_support_cone(seed)
        ]
        best_escape_prefix = max(
            prefix for prefix, _ in escape_prefixes
        )
        best_escapes = [
            seed
            for prefix, seed in escape_prefixes
            if prefix == best_escape_prefix
        ]

        mixed_queries = tuple(
            sp.expand(multiplier * seed**exponent)
            for seed in survivors
            for multiplier in multipliers
            for exponent in range(mixed_start, pure_bound + 1)
        )
        mixed_primitives = (
            npl.particular_primitives(mixed_queries)
            if mixed_queries
            else ()
        )
        tail_obstructions: list[tuple[sp.Expr, sp.Expr]] = []
        query_index = 0
        for seed in survivors:
            for multiplier in multipliers:
                tail = mixed_primitives[query_index : query_index + 3]
                query_index += 3
                if tail and all(
                    not npl.primitive_enters_ideal(
                        primitive, quotient_data
                    )
                    for primitive in tail
                ):
                    tail_obstructions.append((seed, multiplier))

        total_survivors += len(survivors)
        total_tail_obstructions += len(tail_obstructions)
        chart_survivor_counts.append(len(survivors))
        chart_escape_prefixes.append(best_escape_prefix)
        print(
            "CHART:",
            chart.name,
            f"kernel-image dimension={len(kernel_basis)}",
            f"pure-prefix survivors={len(survivors)}",
            f"bounded-tail obstructions={len(tail_obstructions)}",
            f"best outside-support prefix={best_escape_prefix}",
        )
        for seed in survivors[:10]:
            print("  SURVIVOR:", f"f={seed}")
        for seed, multiplier in tail_obstructions[:10]:
            print(
                "  CANDIDATE-ONLY:",
                f"f={seed}",
                f"g={multiplier}",
            )
        for seed in best_escapes[:3]:
            print("  NEAREST-ESCAPE:", f"f={seed}")

    print(f"SUMMARY: pure-prefix survivors={total_survivors}")
    print(
        "SUMMARY:",
        f"distinct pure-prefix survivors={len(distinct_survivors)}",
    )
    print(
        f"SUMMARY: bounded-tail obstructions={total_tail_obstructions}"
    )
    assert total_tail_obstructions == 0
    assert chart_survivor_counts == [75, 54, 54, 54, 52]
    assert chart_escape_prefixes == [0, 0, 0, 0, 0]
    assert len(distinct_survivors) == 75
    assert all(
        npl.in_base_support_cone(survivor)
        for survivor in distinct_survivors.values()
    ), "a survivor escaped the (u,v)-support cone"
    print("NOTE: membership is exact; the exponent range is bounded")


if __name__ == "__main__":
    main()

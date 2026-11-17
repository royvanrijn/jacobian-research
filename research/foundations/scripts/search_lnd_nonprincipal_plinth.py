#!/usr/bin/env python3
"""Exact bounded LNED search for a nonprincipal plinth ideal.

Work over Q[u,v,x,y] with

    D = u*d/dx + v*d/dy,       w = u*y-v*x in ker(D).

The two local slices x/u and y/v meet over the codimension-two plinth
locus (u,v).  With wt(u)=wt(v)=0 and wt(x)=wt(y)=1, D lowers weight by
one.  Singular computes one exact bounded primitive over Q[u,v].
Modulo a zero-dimensional ideal, a finite closure under u, v, and w
then decides whether a kernel correction puts that primitive in the
ideal.

Individual membership decisions are exact.  Pure powers through six and
mixed powers four through six form only a bounded counterexample search.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
import re
import subprocess

import sympy as sp


u, v, x, y = sp.symbols("u v x y")
variables = (u, v, x, y)
quotient_order = (y, x, v, u)
invariant = sp.expand(u * y - v * x)


def derivation(f: sp.Expr) -> sp.Expr:
    return sp.expand(u * sp.diff(f, x) + v * sp.diff(f, y))


def weight_monomials(target_weight: int) -> tuple[sp.Expr, ...]:
    if target_weight < 0:
        return ()
    return tuple(
        x ** (target_weight - y_power) * y**y_power
        for y_power in range(target_weight + 1)
    )


def mixed_weight_basis(maximum_weight: int) -> tuple[sp.Expr, ...]:
    return tuple(
        monomial
        for target_weight in range(maximum_weight + 1)
        for monomial in weight_monomials(target_weight)
    )


def mixed_weight(f: sp.Expr) -> int:
    return max(
        x_power + y_power
        for (_, _, x_power, y_power), coefficient in sp.Poly(
            sp.expand(f), *variables
        ).terms()
        if coefficient
    )


def vector_over_base(
    f: sp.Expr, basis: tuple[sp.Expr, ...]
) -> tuple[sp.Expr, ...]:
    poly = sp.Poly(sp.expand(f), x, y)
    result = tuple(
        sp.expand(poly.coeff_monomial(monomial))
        for monomial in basis
    )
    assert all(not coefficient.has(x, y) for coefficient in result)
    return result


def singular_expr(f: sp.Expr) -> str:
    return str(sp.expand(f)).replace("**", "^")


def singular_vector(vector: tuple[sp.Expr, ...]) -> str:
    return "[" + ",".join(singular_expr(entry) for entry in vector) + "]"


def bounded_primitive_module(
    maximum_target_weight: int,
) -> tuple[
    tuple[sp.Expr, ...],
    tuple[sp.Expr, ...],
    tuple[tuple[sp.Expr, ...], ...],
]:
    target_basis = mixed_weight_basis(maximum_target_weight)
    source_candidates = mixed_weight_basis(maximum_target_weight + 1)
    source_basis: list[sp.Expr] = []
    image_vectors: list[tuple[sp.Expr, ...]] = []
    for source in source_candidates:
        image = derivation(source)
        if image != 0:
            source_basis.append(source)
            image_vectors.append(vector_over_base(image, target_basis))
    return target_basis, tuple(source_basis), tuple(image_vectors)


def particular_primitives(
    expressions: tuple[sp.Expr, ...],
) -> tuple[sp.Expr | None, ...]:
    grouped: dict[int, list[tuple[int, sp.Expr]]] = {}
    for query_id, expression in enumerate(expressions):
        grouped.setdefault(mixed_weight(expression), []).append(
            (query_id, expression)
        )

    code = ["ring r=0,(u,v),dp;", "short=0;"]
    source_bases: dict[int, tuple[sp.Expr, ...]] = {}
    for maximum_weight, queries in sorted(grouped.items()):
        target_basis, source_basis, image_vectors = (
            bounded_primitive_module(maximum_weight)
        )
        source_bases[maximum_weight] = source_basis
        module_name = f"PM{maximum_weight}"
        standard_name = f"PG{maximum_weight}"
        code.append(
            f"module {module_name}="
            + ",".join(
                singular_vector(vector) for vector in image_vectors
            )
            + ";"
        )
        code.append(f"module {standard_name}=std({module_name});")
        for query_id, expression in queries:
            target_vector = vector_over_base(expression, target_basis)
            vector_name = f"pv{query_id}"
            lift_name = f"pl{query_id}"
            code.append(
                f"vector {vector_name}="
                f"{singular_vector(target_vector)};"
            )
            code.append(
                f"if (reduce({vector_name},{standard_name})==0)"
                "{"
                f"matrix {lift_name}=lift({module_name},{vector_name});"
                f'print("Q{query_id}=1");'
            )
            for source_index in range(len(source_basis)):
                code.append(
                    f'print("C{query_id}_{source_index}="'
                    f"+string({lift_name}[{source_index + 1},1]));"
                )
            code.append(f'}}else{{print("Q{query_id}=0");}};')
    code.append("quit;")

    completed = subprocess.run(
        ["Singular", "-q"],
        input="\n".join(code),
        text=True,
        capture_output=True,
        check=True,
    )
    membership = {
        int(query_id): value == "1"
        for query_id, value in re.findall(
            r"Q(\d+)=(0|1)", completed.stdout
        )
    }
    assert len(membership) == len(expressions), (
        completed.stdout + completed.stderr
    )
    coefficients: dict[tuple[int, int], sp.Expr] = {}
    for query_id, source_index, value in re.findall(
        r"C(\d+)_(\d+)=(.*)", completed.stdout
    ):
        cleaned = value.strip().strip("`").replace("^", "**")
        coefficients[(int(query_id), int(source_index))] = sp.sympify(
            cleaned, locals={"u": u, "v": v}
        )

    primitives: list[sp.Expr | None] = []
    for query_id, expression in enumerate(expressions):
        if not membership[query_id]:
            primitives.append(None)
            continue
        source_basis = source_bases[mixed_weight(expression)]
        primitive = sp.expand(
            sum(
                coefficients[(query_id, source_index)] * source
                for source_index, source in enumerate(source_basis)
            )
        )
        assert derivation(primitive) == expression
        primitives.append(primitive)
    return tuple(primitives)


@dataclass(frozen=True)
class PlinthIdeal:
    name: str
    generators: tuple[sp.Expr, ...]


left_primary = (u, v**2, x, y)
right_primary = (u**2, v, x - u, y - v)
primary_product = tuple(
    sp.expand(left * right)
    for left in left_primary
    for right in right_primary
)

CHARTS = (
    PlinthIdeal("plinth-radical", (u, v, x, y)),
    PlinthIdeal("double-plinth", (u**2, v**2, x, y)),
    PlinthIdeal(
        "balanced-plinth-jet",
        (
            u**2,
            v**2,
            u * x,
            v * x,
            u * y,
            v * y,
            x**2,
            x * y,
            y**2,
        ),
    ),
    PlinthIdeal(
        "tilted-plinth-jet",
        (
            u**2,
            v**2,
            u * x,
            v * x,
            u * y,
            v * y,
            x**2 + u * v,
            x * y + u * v,
            y**2,
        ),
    ),
    PlinthIdeal("plinth-primary-product", primary_product),
)


def polynomial_span_contains(
    basis: tuple[sp.Expr, ...], candidate: sp.Expr
) -> bool:
    expressions = basis + (sp.expand(candidate),)
    monomials = sorted(
        {
            exponents
            for expression in expressions
            for exponents, coefficient in sp.Poly(
                expression, *quotient_order, domain=sp.QQ
            ).terms()
            if coefficient
        }
    )
    matrix = sp.Matrix(
        [
            [
                sp.Poly(
                    expression, *quotient_order, domain=sp.QQ
                ).coeff_monomial(exponents)
                for expression in expressions
            ]
            for exponents in monomials
        ]
    )
    return matrix[:, :-1].rank() == matrix.rank()


def quotient_kernel_span(chart: PlinthIdeal) -> tuple[
    sp.GroebnerBasis, tuple[sp.Expr, ...]
]:
    groebner = sp.groebner(
        chart.generators,
        *quotient_order,
        order="lex",
        domain=sp.QQ,
    )

    def remainder(f: sp.Expr) -> sp.Expr:
        return sp.expand(groebner.reduce(sp.expand(f))[1])

    basis: list[sp.Expr] = [remainder(sp.Integer(1))]
    frontier = list(basis)
    for _ in range(128):
        next_frontier: list[sp.Expr] = []
        for element in frontier:
            for multiplier in (u, v, invariant):
                candidate = remainder(element * multiplier)
                if candidate != 0 and not polynomial_span_contains(
                    tuple(basis), candidate
                ):
                    basis.append(candidate)
                    next_frontier.append(candidate)
        if not next_frontier:
            return groebner, tuple(basis)
        frontier = next_frontier
    raise AssertionError("kernel-image closure did not stabilize")


def primitive_enters_ideal(
    primitive: sp.Expr | None,
    quotient_data: tuple[sp.GroebnerBasis, tuple[sp.Expr, ...]],
) -> bool:
    if primitive is None:
        return False
    groebner, kernel_basis = quotient_data
    remainder = sp.expand(groebner.reduce(primitive)[1])
    return polynomial_span_contains(kernel_basis, remainder)


def in_base_support_cone(f: sp.Expr) -> bool:
    return all(
        u_power + v_power >= 1
        for (u_power, v_power, _, _), coefficient in sp.Poly(
            f, *variables
        ).terms()
        if coefficient
    )


def seed_family() -> tuple[sp.Expr, ...]:
    coefficients = (sp.Integer(1), u, v, u - v, u * v)
    seeds: list[sp.Expr] = []
    for target_weight in range(1, 4):
        monomials = weight_monomials(target_weight)
        atoms = tuple(
            sp.expand(coefficient * monomial)
            for coefficient in coefficients
            for monomial in monomials
        )
        seeds.extend(atoms)
        for left, right in combinations(atoms, 2):
            seeds.append(sp.expand(left + right))
            seeds.append(sp.expand(left - right))

    mixed_atoms = tuple(
        sp.expand(coefficient * monomial)
        for coefficient in coefficients
        for monomial in (x, y, x**2, x * y, y**2)
    )
    for left, right in combinations(mixed_atoms, 2):
        if mixed_weight(left) != mixed_weight(right):
            seeds.append(sp.expand(left + right))
            seeds.append(sp.expand(left - right))
    seeds.extend(
        (
            invariant,
            u * invariant,
            v * invariant,
            u * x + v * y,
            u * y + v * x,
            x**2 + u,
            y**2 + v,
        )
    )
    for chart in CHARTS:
        for generator in chart.generators:
            for multiplier in (sp.Integer(1), u, v, x, y):
                image = derivation(sp.expand(generator * multiplier))
                if image != 0 and mixed_weight(image) <= 3:
                    seeds.append(image)

    unique: dict[str, sp.Expr] = {}
    for seed in seeds:
        expanded = sp.expand(seed)
        if expanded != 0:
            unique[sp.srepr(expanded)] = expanded
    return tuple(unique.values())


def main() -> None:
    pure_bound = 6
    mixed_start = 4
    seeds = seed_family()
    multipliers = (sp.Integer(1), u, v, x, y, invariant)

    assert derivation(invariant) == 0
    pure_queries = tuple(
        sp.expand(seed**exponent)
        for seed in seeds
        for exponent in range(1, pure_bound + 1)
    )
    pure_primitives = particular_primitives(pure_queries)

    print(
        "SEARCH:",
        "D=u*d_x+v*d_y;",
        f"{len(seeds)} nonprincipal-plinth sparse seeds",
    )
    total_survivors = 0
    total_tail_obstructions = 0
    distinct_survivors: dict[str, sp.Expr] = {}
    chart_survivor_counts: list[int] = []
    chart_escape_prefixes: list[int] = []

    for chart in CHARTS:
        quotient_data = quotient_kernel_span(chart)
        membership_rows = [
            tuple(
                primitive_enters_ideal(
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
            if not in_base_support_cone(seed)
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
            particular_primitives(mixed_queries)
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
                    not primitive_enters_ideal(
                        primitive, quotient_data
                    )
                    for primitive in tail
                ):
                    tail_obstructions.append((seed, multiplier))

        total_survivors += len(survivors)
        chart_survivor_counts.append(len(survivors))
        chart_escape_prefixes.append(best_escape_prefix)
        total_tail_obstructions += len(tail_obstructions)
        print(
            "CHART:",
            chart.name,
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
    assert chart_survivor_counts == [100, 100, 98, 98, 98]
    assert chart_escape_prefixes == [0, 0, 0, 0, 0]
    assert len(distinct_survivors) == 100
    assert total_tail_obstructions == 0
    assert all(
        in_base_support_cone(survivor)
        for survivor in distinct_survivors.values()
    ), "a survivor escaped the (u,v)-support cone"
    print("NOTE: membership is exact; the exponent range is bounded")


if __name__ == "__main__":
    main()

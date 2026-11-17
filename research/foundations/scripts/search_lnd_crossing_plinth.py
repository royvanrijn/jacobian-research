#!/usr/bin/env python3
"""Exact bounded LNED search at an intersecting plinth divisor.

Work over Q[u,v,y,z] with

    D = u*v*d/dy + (y^2+u)*d/dz.

The plinth divisor u*v=0 has two components meeting at u=v=0.  Relative
to wt(u)=wt(v)=0, wt(y)=1, wt(z)=3, the leading derivation lowers weight
by one and the perturbation u*d/dz lowers weight by three.  Every
primitive can be normalized modulo

    Q[u,v, 3*u*v*z-y^3-3*u*y]

to weight at most one above its target.  Singular computes one exact
bounded primitive over Q[u,v]; a finite quotient/kernel calculation then
decides whether a kernel correction places it in the chosen ideal.

Individual membership decisions are exact.  Pure powers through six and
mixed powers four through six are only a bounded counterexample search.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
import re
import subprocess

import sympy as sp


u, v, y, z = sp.symbols("u v y z")
base_product = u * v
variables = (u, v, y, z)
quotient_order = (z, y, v, u)
z_weight = 3
invariant = sp.expand(3 * u * v * z - y**3 - 3 * u * y)


def derivation(f: sp.Expr) -> sp.Expr:
    return sp.expand(
        base_product * sp.diff(f, y)
        + (y**2 + u) * sp.diff(f, z)
    )


def weight_monomials(target_weight: int) -> tuple[sp.Expr, ...]:
    if target_weight < 0:
        return ()
    return tuple(
        y ** (target_weight - z_weight * z_power) * z**z_power
        for z_power in range(target_weight // z_weight + 1)
    )


def mixed_weight_basis(maximum_weight: int) -> tuple[sp.Expr, ...]:
    return tuple(
        monomial
        for target_weight in range(maximum_weight + 1)
        for monomial in weight_monomials(target_weight)
    )


def mixed_weight(f: sp.Expr) -> int:
    return max(
        y_power + z_weight * z_power
        for (_, _, y_power, z_power), coefficient in sp.Poly(
            sp.expand(f), *variables
        ).terms()
        if coefficient
    )


def vector_over_base(
    f: sp.Expr, basis: tuple[sp.Expr, ...]
) -> tuple[sp.Expr, ...]:
    poly = sp.Poly(sp.expand(f), y, z)
    result = tuple(
        sp.expand(poly.coeff_monomial(monomial))
        for monomial in basis
    )
    assert all(not coefficient.has(y, z) for coefficient in result)
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
        maximum_weight = mixed_weight(expression)
        source_basis = source_bases[maximum_weight]
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
class CrossingIdeal:
    name: str
    generators: tuple[sp.Expr, ...]


left_primary = (u, v**2, y, z)
right_primary = (u**2, v, y - u, z - y)
primary_product = tuple(
    sp.expand(left * right)
    for left in left_primary
    for right in right_primary
)

CHARTS = (
    CrossingIdeal("crossing-radical", (u, v, y, z)),
    CrossingIdeal("double-crossing", (u**2, v**2, y, z)),
    CrossingIdeal(
        "balanced-crossing-jet",
        (
            u**2,
            v**2,
            u * y,
            v * y,
            y**2 + u * v * z,
            y * z,
            z**2,
        ),
    ),
    CrossingIdeal(
        "tilted-crossing-jet",
        (
            u**2,
            v**2,
            u * y,
            v * y,
            y**2 + u * v * z + u * v,
            y * z + u * v,
            z**2,
        ),
    ),
    CrossingIdeal("crossing-primary-product", primary_product),
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


def quotient_kernel_span(chart: CrossingIdeal) -> tuple[
    sp.GroebnerBasis, tuple[sp.Expr, ...]
]:
    groebner = sp.groebner(
        chart.generators,
        *quotient_order,
        order="lex",
        domain=sp.QQ,
    )
    assert groebner.is_zero_dimensional

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


def divisible_by_full_plinth(f: sp.Expr) -> bool:
    return all(
        u_power >= 1 and v_power >= 1
        for (u_power, v_power, _, _), coefficient in sp.Poly(
            f, *variables
        ).terms()
        if coefficient
    )


def seed_family() -> tuple[sp.Expr, ...]:
    coefficients = (
        sp.Integer(1),
        u,
        v,
        u - v,
        u + v,
        u**2,
        u * v,
        v**2,
    )
    seeds: list[sp.Expr] = []
    for target_weight in range(1, 5):
        monomials = weight_monomials(target_weight)
        atoms = tuple(
            sp.expand(coefficient * monomial)
            for coefficient in coefficients
            for monomial in monomials
        )
        seeds.extend(atoms)
        if len(monomials) > 1:
            for left, right in combinations(atoms, 2):
                seeds.append(sp.expand(left + right))
                seeds.append(sp.expand(left - right))

    mixed_atoms = tuple(
        sp.expand(coefficient * monomial)
        for coefficient in coefficients
        for monomial in (y, y**2, z)
    )
    for left, right in combinations(mixed_atoms, 2):
        if mixed_weight(left) != mixed_weight(right):
            seeds.append(sp.expand(left + right))
            seeds.append(sp.expand(left - right))
    seeds.extend(
        (
            y**2 + u,
            u * v * (y**2 + u),
            y**2 + u + u * v * y,
            y**2 + u - u * v * y,
            z + y**2 + u,
            z - y**2 - u,
            invariant,
            u * v * invariant,
        )
    )
    for chart in CHARTS:
        for generator in chart.generators:
            for multiplier in (sp.Integer(1), u, v, y):
                image = derivation(sp.expand(generator * multiplier))
                if image != 0 and mixed_weight(image) <= 4:
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
    multipliers = (sp.Integer(1), u, v, y, z)

    assert derivation(invariant) == 0
    pure_queries = tuple(
        sp.expand(seed**exponent)
        for seed in seeds
        for exponent in range(1, pure_bound + 1)
    )
    pure_primitives = particular_primitives(pure_queries)

    print(
        "SEARCH:",
        "D=u*v*d_y+(y^2+u)*d_z;",
        f"{len(seeds)} crossing-aware sparse seeds",
    )
    total_survivors = 0
    total_tail_obstructions = 0
    distinct_survivors: dict[str, sp.Expr] = {}
    chart_survivor_counts: list[int] = []

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
            if not divisible_by_full_plinth(seed)
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
        total_tail_obstructions += len(tail_obstructions)
        print(
            "CHART:",
            chart.name,
            f"pure-prefix survivors={len(survivors)}",
            f"bounded-tail obstructions={len(tail_obstructions)}",
            f"best nonplinth prefix={best_escape_prefix}",
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
    assert chart_survivor_counts == [36, 36, 35, 35, 35]
    assert len(distinct_survivors) == 36
    assert total_tail_obstructions == 0
    assert all(
        divisible_by_full_plinth(survivor)
        for survivor in distinct_survivors.values()
    ), "a survivor escaped full-plinth divisibility"
    print("NOTE: membership is exact; the exponent range is bounded")


if __name__ == "__main__":
    main()

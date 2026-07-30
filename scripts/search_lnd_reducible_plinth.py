#!/usr/bin/env python3
"""Exact reducible-plinth LNED searches over a one-dimensional base.

The locally nilpotent derivations include

    D_r = x*(x-1)*d/dy + y^r*d/dz,  r=1,2.

Give x weight zero, y weight one, and z weight r+1.  Then D_r lowers
weight by exactly one, so membership of a weight-W polynomial in D_r(I)
is an exact finite module problem over Q[x], not a degree-truncated
approximation.  Singular computes the required Q[x]-module standard
bases.

The final profile adds the branch-asymmetric term x*d/dz.  Its grading is
broken, but every primitive can be normalized modulo the full kernel to
weight at most one above the target.  Exact Q[x]-module lifts followed by
finite quotient/kernel tests therefore retain exact membership decisions.

The ideals below couple the two plinth fibers x=0 and x=1.  Sparse seeds
carry coefficients on one or both branches.  Six pure powers and the mixed
tail for five multipliers are candidate generation only; promotion still
requires an all-order certificate.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
import re
import subprocess

import sympy as sp


x, y, z = sp.symbols("x y z")
a = sp.expand(x * (x - 1))
variables = (x, y, z)


def derivation(f: sp.Expr, chain_power: int) -> sp.Expr:
    return sp.expand(
        a * sp.diff(f, y) + y**chain_power * sp.diff(f, z)
    )


def weight_monomials(
    target_weight: int, z_weight: int
) -> tuple[sp.Expr, ...]:
    if target_weight < 0:
        return ()
    return tuple(
        y ** (target_weight - z_weight * z_power) * z**z_power
        for z_power in range(target_weight // z_weight + 1)
    )


def weight(f: sp.Expr, z_weight: int) -> int:
    weights = {
        y_power + z_weight * z_power
        for (_, y_power, z_power), coefficient in sp.Poly(
            sp.expand(f), *variables
        ).terms()
        if coefficient
    }
    assert len(weights) == 1
    return weights.pop()


def vector_over_qx(
    f: sp.Expr, target_weight: int, z_weight: int
) -> tuple[sp.Expr, ...]:
    poly = sp.Poly(sp.expand(f), y, z)
    result = tuple(
        sp.expand(poly.coeff_monomial(monomial))
        for monomial in weight_monomials(target_weight, z_weight)
    )
    assert all(not coefficient.has(y, z) for coefficient in result)
    return result


def singular_expr(f: sp.Expr) -> str:
    return str(sp.expand(f)).replace("**", "^")


def singular_vector(vector: tuple[sp.Expr, ...]) -> str:
    return "[" + ",".join(singular_expr(entry) for entry in vector) + "]"


@dataclass(frozen=True)
class WeightedIdeal:
    name: str
    generators: tuple[sp.Expr, ...]


def charts(chain_power: int) -> tuple[WeightedIdeal, ...]:
    q0 = (x**2, y, z)
    q1 = (x - 1, y**2, z)
    asymmetric_product = tuple(
        sp.expand(left * right)
        for left in q0
        for right in q1
    )
    invariant = (chain_power + 1) * a * z - y ** (chain_power + 1)
    return (
        WeightedIdeal("radical-fiber-pair", (a, y, z)),
        WeightedIdeal("double-base-fiber", (a**2, y, z)),
        WeightedIdeal(
            "invariant-coupling",
            (a**2, a * y, invariant, y * z, z**2),
        ),
        WeightedIdeal(
            "tilted-cross-jet",
            (
                a**2,
                a * y,
                y ** (chain_power + 1) + a * z,
                y * z,
                z**2,
            ),
        ),
        WeightedIdeal("asymmetric-primary-pair", asymmetric_product),
    )


def image_module_generators(
    chart: WeightedIdeal,
    target_weight: int,
    chain_power: int,
) -> tuple[tuple[sp.Expr, ...], ...]:
    z_weight = chain_power + 1
    source_weight = target_weight + 1
    vectors: list[tuple[sp.Expr, ...]] = []
    for generator in chart.generators:
        generator_weight = weight(generator, z_weight)
        for multiplier in weight_monomials(
            source_weight - generator_weight, z_weight
        ):
            source = sp.expand(generator * multiplier)
            image = derivation(source, chain_power)
            assert image == 0 or weight(image, z_weight) == target_weight
            if image != 0:
                vectors.append(
                    vector_over_qx(image, target_weight, z_weight)
                )
    if not vectors:
        vectors.append(
            tuple(
                sp.Integer(0)
                for _ in weight_monomials(target_weight, z_weight)
            )
        )
    return tuple(vectors)


def exact_memberships(
    chart: WeightedIdeal,
    expressions: tuple[sp.Expr, ...],
    chain_power: int,
) -> tuple[bool, ...]:
    z_weight = chain_power + 1
    grouped: dict[int, list[tuple[int, sp.Expr]]] = {}
    for query_id, expression in enumerate(expressions):
        grouped.setdefault(weight(expression, z_weight), []).append(
            (query_id, expression)
        )

    code = ["ring r=0,x,dp;", "short=0;"]
    for target_weight, queries in sorted(grouped.items()):
        generators = image_module_generators(
            chart, target_weight, chain_power
        )
        module_name = f"M{target_weight}"
        standard_name = f"G{target_weight}"
        code.append(
            f"module {module_name}="
            + ",".join(singular_vector(vector) for vector in generators)
            + ";"
        )
        code.append(f"module {standard_name}=std({module_name});")
        for query_id, expression in queries:
            vector = vector_over_qx(
                expression, target_weight, z_weight
            )
            vector_name = f"v{query_id}"
            code.append(f"vector {vector_name}={singular_vector(vector)};")
            code.append(
                f'if (reduce({vector_name},{standard_name})==0) '
                f'{{print("Q{query_id}=1");}} '
                f'else {{print("Q{query_id}=0");}};'
            )
    code.append("quit;")

    completed = subprocess.run(
        ["Singular", "-q"],
        input="\n".join(code),
        text=True,
        capture_output=True,
        check=True,
    )
    answers = {
        int(query_id): value == "1"
        for query_id, value in re.findall(r"Q(\d+)=(0|1)", completed.stdout)
    }
    assert len(answers) == len(expressions), completed.stdout + completed.stderr
    return tuple(answers[index] for index in range(len(expressions)))


def seed_family(chain_power: int) -> tuple[sp.Expr, ...]:
    z_weight = chain_power + 1
    branch_coefficients = (sp.Integer(1), x, x - 1, 2 * x - 1, a)
    weight_one = tuple(coefficient * y for coefficient in branch_coefficients)

    higher_weight: list[sp.Expr] = []
    for target_weight in range(2, z_weight + 2):
        monomials = weight_monomials(target_weight, z_weight)
        atoms = tuple(
            coefficient * monomial
            for coefficient in branch_coefficients
            for monomial in monomials
        )
        higher_weight.extend(atoms)
        if len(monomials) > 1:
            for left, right in combinations(atoms, 2):
                higher_weight.append(sp.expand(left + right))
                higher_weight.append(sp.expand(left - right))
    return weight_one + tuple(higher_weight)


def search_profile(chain_power: int) -> tuple[int, int]:
    pure_bound = 6
    mixed_start = 4
    z_weight = chain_power + 1
    seeds = seed_family(chain_power)
    multipliers = (sp.Integer(1), x, x - 1, y, z)
    invariant = (chain_power + 1) * a * z - y ** (chain_power + 1)

    assert derivation(x, chain_power) == 0
    assert derivation(y, chain_power) == a
    assert derivation(z, chain_power) == y**chain_power
    assert derivation(invariant, chain_power) == 0
    assert weight(derivation(z, chain_power), z_weight) == chain_power

    print(
        "SEARCH:",
        f"D=x*(x-1)*d_y+y^{chain_power}*d_z;",
        f"{len(seeds)} branch-aware sparse seeds",
    )
    total_survivors = 0
    total_tail_obstructions = 0

    for chart in charts(chain_power):
        pure_queries = tuple(
            sp.expand(seed**exponent)
            for seed in seeds
            for exponent in range(1, pure_bound + 1)
        )
        pure_answers = exact_memberships(
            chart, pure_queries, chain_power
        )
        survivors = [
            seed
            for seed_index, seed in enumerate(seeds)
            if all(
                pure_answers[seed_index * pure_bound + exponent - 1]
                for exponent in range(1, pure_bound + 1)
            )
        ]
        assert all(
            sp.expand(sp.div(seed, a, x)[1]) == 0
            for seed in survivors
        )

        mixed_queries = tuple(
            sp.expand(multiplier * seed**exponent)
            for seed in survivors
            for multiplier in multipliers
            for exponent in range(mixed_start, pure_bound + 1)
        )
        mixed_answers = (
            exact_memberships(chart, mixed_queries, chain_power)
            if mixed_queries
            else ()
        )
        tail_obstructions: list[tuple[sp.Expr, sp.Expr]] = []
        query_index = 0
        for seed in survivors:
            for multiplier in multipliers:
                tail = mixed_answers[query_index : query_index + 3]
                query_index += 3
                if tail and not any(tail):
                    tail_obstructions.append((seed, multiplier))

        total_survivors += len(survivors)
        total_tail_obstructions += len(tail_obstructions)
        print(
            "CHART:",
            chart.name,
            f"pure-prefix survivors={len(survivors)}",
            f"bounded-tail obstructions={len(tail_obstructions)}",
        )
        for seed in survivors[:8]:
            print("  SURVIVOR:", f"f={seed}")
        for seed, multiplier in tail_obstructions[:8]:
            print(
                "  CANDIDATE-ONLY:",
                f"f={seed}",
                f"g={multiplier}",
            )

    print(
        f"SUMMARY r={chain_power}: "
        f"pure-prefix survivors={total_survivors}"
    )
    print(
        f"SUMMARY r={chain_power}: "
        f"bounded-tail obstructions={total_tail_obstructions}"
    )
    return total_survivors, total_tail_obstructions


def mixed_weight(f: sp.Expr, z_weight: int) -> int:
    return max(
        y_power + z_weight * z_power
        for (_, y_power, z_power), coefficient in sp.Poly(
            sp.expand(f), *variables
        ).terms()
        if coefficient
    )


def mixed_weight_basis(
    maximum_weight: int, z_weight: int
) -> tuple[sp.Expr, ...]:
    return tuple(
        monomial
        for target_weight in range(maximum_weight + 1)
        for monomial in weight_monomials(target_weight, z_weight)
    )


def mixed_vector_over_qx(
    f: sp.Expr, basis: tuple[sp.Expr, ...]
) -> tuple[sp.Expr, ...]:
    poly = sp.Poly(sp.expand(f), y, z)
    result = tuple(
        sp.expand(poly.coeff_monomial(monomial))
        for monomial in basis
    )
    assert all(not coefficient.has(y, z) for coefficient in result)
    return result


def asymmetric_derivation(f: sp.Expr) -> sp.Expr:
    return sp.expand(
        a * sp.diff(f, y) + (y**2 + x) * sp.diff(f, z)
    )


asymmetric_invariant = sp.expand(3 * a * z - y**3 - 3 * x * y)


def bounded_primitive_modules(
    maximum_target_weight: int,
) -> tuple[
    tuple[sp.Expr, ...],
    tuple[sp.Expr, ...],
    tuple[tuple[sp.Expr, ...], ...],
]:
    z_weight = 3
    target_basis = mixed_weight_basis(maximum_target_weight, z_weight)
    source_candidates = mixed_weight_basis(
        maximum_target_weight + 1, z_weight
    )
    source_basis: list[sp.Expr] = []
    image_vectors: list[tuple[sp.Expr, ...]] = []
    for source in source_candidates:
        image = asymmetric_derivation(source)
        if image != 0:
            source_basis.append(source)
            image_vectors.append(
                mixed_vector_over_qx(image, target_basis)
            )
    return target_basis, tuple(source_basis), tuple(image_vectors)


def particular_primitives(
    expressions: tuple[sp.Expr, ...],
) -> tuple[sp.Expr | None, ...]:
    z_weight = 3
    grouped: dict[int, list[tuple[int, sp.Expr]]] = {}
    for query_id, expression in enumerate(expressions):
        grouped.setdefault(
            mixed_weight(expression, z_weight), []
        ).append((query_id, expression))

    code = ["ring r=0,x,dp;", "short=0;"]
    source_bases: dict[int, tuple[sp.Expr, ...]] = {}
    for maximum_weight, queries in sorted(grouped.items()):
        target_basis, source_basis, image_vectors = (
            bounded_primitive_modules(maximum_weight)
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
            target_vector = mixed_vector_over_qx(
                expression, target_basis
            )
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
            cleaned
        )

    primitives: list[sp.Expr | None] = []
    for query_id, expression in enumerate(expressions):
        if not membership[query_id]:
            primitives.append(None)
            continue
        maximum_weight = mixed_weight(expression, z_weight)
        source_basis = source_bases[maximum_weight]
        primitive = sp.expand(
            sum(
                coefficients[(query_id, source_index)] * source
                for source_index, source in enumerate(source_basis)
            )
        )
        assert asymmetric_derivation(primitive) == expression
        primitives.append(primitive)
    return tuple(primitives)


def polynomial_span_contains(
    basis: tuple[sp.Expr, ...], candidate: sp.Expr
) -> bool:
    expressions = basis + (sp.expand(candidate),)
    monomials = sorted(
        {
            exponents
            for expression in expressions
            for exponents, coefficient in sp.Poly(
                expression, z, y, x, domain=sp.QQ
            ).terms()
            if coefficient
        }
    )
    matrix = sp.Matrix(
        [
            [
                sp.Poly(expression, z, y, x, domain=sp.QQ)
                .coeff_monomial(exponents)
                for expression in expressions
            ]
            for exponents in monomials
        ]
    )
    return matrix[:, :-1].rank() == matrix.rank()


def quotient_kernel_span(chart: WeightedIdeal) -> tuple[
    sp.GroebnerBasis, tuple[sp.Expr, ...]
]:
    groebner = sp.groebner(
        chart.generators, z, y, x, order="lex", domain=sp.QQ
    )

    def remainder(f: sp.Expr) -> sp.Expr:
        return sp.expand(groebner.reduce(sp.expand(f))[1])

    basis: list[sp.Expr] = [remainder(sp.Integer(1))]
    frontier = list(basis)
    for _ in range(64):
        next_frontier: list[sp.Expr] = []
        for element in frontier:
            for multiplier in (x, asymmetric_invariant):
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


def asymmetric_charts() -> tuple[WeightedIdeal, ...]:
    q0 = (x**2, y, z)
    q1 = ((x - 1) ** 2, y - (x - 1), z - y)
    nonhomogeneous_product = tuple(
        sp.expand(left * right)
        for left in q0
        for right in q1
    )
    return charts(2) + (
        WeightedIdeal(
            "mixed-weight-jet",
            (
                a**2,
                a * y,
                y**2 + a * z + a * (2 * x - 1),
                y * z + a,
                z**2,
            ),
        ),
        WeightedIdeal(
            "asymmetric-offset-jet",
            (
                a**2,
                a * y,
                y**2 + a * z + x * a,
                y * z + (x - 1) * a,
                z**2,
            ),
        ),
        WeightedIdeal(
            "nonhomogeneous-primary-pair",
            nonhomogeneous_product,
        ),
    )


def asymmetric_seed_family() -> tuple[sp.Expr, ...]:
    seeds = list(seed_family(2))
    branch_coefficients = (sp.Integer(1), x, x - 1, 2 * x - 1, a)
    atoms = tuple(
        sp.expand(coefficient * monomial)
        for coefficient in branch_coefficients
        for monomial in (y, y**2, z)
    )
    seeds.extend(atoms)
    for left, right in combinations(atoms, 2):
        if mixed_weight(left, 3) != mixed_weight(right, 3):
            seeds.append(sp.expand(left + right))
            seeds.append(sp.expand(left - right))
    seeds.extend(
        (
            y**2 + x,
            a * (y**2 + x),
            y**2 + x + a * y,
            y**2 + x - a * y,
            z + y**2 + x,
            z - y**2 - x,
        )
    )
    unique: dict[str, sp.Expr] = {}
    for seed in seeds:
        expanded = sp.expand(seed)
        if expanded != 0:
            unique[sp.srepr(expanded)] = expanded
    return tuple(unique.values())


def search_asymmetric_profile() -> tuple[int, int]:
    pure_bound = 6
    mixed_start = 4
    seeds = asymmetric_seed_family()
    multipliers = (sp.Integer(1), x, x - 1, y, z)
    profile_charts = asymmetric_charts()

    assert asymmetric_derivation(asymmetric_invariant) == 0
    pure_queries = tuple(
        sp.expand(seed**exponent)
        for seed in seeds
        for exponent in range(1, pure_bound + 1)
    )
    pure_primitives = particular_primitives(pure_queries)

    print(
        "SEARCH:",
        "D=x*(x-1)*d_y+(y^2+x)*d_z;",
        f"{len(seeds)} branch-aware sparse seeds",
    )
    total_survivors = 0
    total_tail_obstructions = 0
    for chart in profile_charts:
        quotient_data = quotient_kernel_span(chart)
        survivors = [
            seed
            for seed_index, seed in enumerate(seeds)
            if all(
                primitive_enters_ideal(
                    pure_primitives[
                        seed_index * pure_bound + exponent - 1
                    ],
                    quotient_data,
                )
                for exponent in range(1, pure_bound + 1)
            )
        ]
        assert all(
            sp.expand(sp.div(seed, a, x)[1]) == 0
            for seed in survivors
        )
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
        total_tail_obstructions += len(tail_obstructions)
        print(
            "CHART:",
            chart.name,
            f"pure-prefix survivors={len(survivors)}",
            f"bounded-tail obstructions={len(tail_obstructions)}",
        )
        for seed in survivors[:8]:
            print("  SURVIVOR:", f"f={seed}")
        for seed, multiplier in tail_obstructions[:8]:
            print(
                "  CANDIDATE-ONLY:",
                f"f={seed}",
                f"g={multiplier}",
            )

    print(
        "SUMMARY asymmetric:",
        f"pure-prefix survivors={total_survivors}",
    )
    print(
        "SUMMARY asymmetric:",
        f"bounded-tail obstructions={total_tail_obstructions}",
    )
    return total_survivors, total_tail_obstructions


def main() -> None:
    results = tuple(search_profile(chain_power) for chain_power in (1, 2))
    assert results[0] == (14, 0)
    asymmetric_result = search_asymmetric_profile()
    assert asymmetric_result == (48, 0)
    print("NOTE: module membership is exact; the exponent range is bounded")


if __name__ == "__main__":
    main()

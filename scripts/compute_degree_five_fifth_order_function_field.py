#!/usr/bin/env python3
"""Compute the generic bounded hbar^5 period by a fraction-free kernel.

This is the characteristic-zero continuation of
``explore_degree_five_fifth_order_period_constraints.py``.  It constructs
the universal sparse hbar^3 lift and the five lower-kernel directions which
give the stable 15-condition pivot set.  SymPy is used for the sparse Weyl
calculus over ``QQ(a,tau)``.  The final 15-by-16 kernel is delegated to
Singular's syzygy engine, avoiding rational Gaussian-elimination swell.

The five selectable support charts produce *bounded* period candidates.  A
complete Laurent certificate still requires checking all additional Laurent
lower-lift directions.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import tempfile
from itertools import combinations_with_replacement
from pathlib import Path
from time import perf_counter

from sympy.polys.domains import QQ

from explore_degree_five_quantum_residue import (
    add,
    degree_five_family,
    laurent_monomials,
    pi_power,
    poisson,
    scale,
    solve_affine,
    third_order_family,
)
from verify_degree_five_laurent_quantum_obstruction import SUPPORT
from verify_degree_five_third_order_function_field import (
    S_SUPPORT,
    T_SUPPORT,
    clear_coefficient,
)


SELECTED_DIRECTIONS = (14, 18, 22, 33, 40)
FULL_DIRECTION_GROUPS = (
    (14, 15),
    (16, 17),
    (18, 19),
    (20, 21),
    (22, 23),
    (29, 30),
    (31, 32),
    (33,),
    (40,),
)
SECOND_SUPPORT = tuple(
    monomial
    for monomial in SUPPORT
    if monomial != (19, 7, 0)
) + ((20, 0, 4),)
THIRD_SUPPORT = tuple(
    monomial
    for monomial in SECOND_SUPPORT
    if monomial != (19, 5, 1)
) + ((19, 7, 0),)
FOURTH_SUPPORT = tuple(
    monomial
    for monomial in SUPPORT
    if monomial != (19, 7, 0)
) + ((20, 2, 3),)
FIFTH_SUPPORT = tuple(
    monomial
    for monomial in SUPPORT
    if monomial != (19, 7, 0)
) + ((20, 4, 2),)
D5_S_MONOMIALS = (
    (10, 0, 0),
    (11, 1, 0),
    (12, 2, 0),
    (13, 3, 0),
    (14, 4, 0),
    (15, 5, 0),
    (12, 0, 1),
    (13, 1, 1),
    (14, 2, 1),
    (15, 3, 1),
)


def split_solution(vector, s_support, t_support):
    split = len(s_support)
    return (
        {
            s_support[index]: coefficient
            for index, coefficient in vector.items()
            if index < split
        },
        {
            t_support[index - split]: coefficient
            for index, coefficient in vector.items()
            if index >= split
        },
    )


def to_singular(value) -> str:
    return str(value).replace("**", "^")


def construct_system(
    full_bounded: bool = True,
    functional_support=SUPPORT,
):
    # Only the supports of these numerical kernel vectors are used.  Each
    # support has function-field nullity one below.
    sample_s, sample_t = degree_five_family(
        QQ,
        -QQ(2) / QQ(3),
        QQ.one,
    )
    sample_family = third_order_family(sample_s, sample_t, QQ)
    direction_supports = {
        index: (
            sorted(sample_family.kernel[index][0]),
            sorted(sample_family.kernel[index][1]),
        )
        for index in set().union(*map(set, FULL_DIRECTION_GROUPS))
    }

    field = QQ.frac_field("a", "tau")
    a, tau = field.gens
    ring = a.numer.ring
    S, T = degree_five_family(
        field,
        a,
        tau,
        verify_canonical=False,
    )
    a_poly = a.numer
    denominator_s = a_poly**4 * (a_poly + ring.one) ** 10
    denominator_t = a_poly**4 * (a_poly + ring.one) ** 8
    A = {
        monomial: clear_coefficient(
            coefficient,
            denominator_s,
            field,
            ring,
        )
        for monomial, coefficient in S.items()
    }
    B = {
        monomial: clear_coefficient(
            coefficient,
            denominator_t,
            field,
            ring,
        )
        for monomial, coefficient in T.items()
    }

    def as_field(poly):
        return {
            monomial: field(coefficient)
            for monomial, coefficient in poly.items()
        }

    base_columns = [
        as_field(
            scale(
                poisson({monomial: ring.one}, B),
                denominator_s,
            )
        )
        for monomial in S_SUPPORT
    ]
    base_columns += [
        as_field(
            scale(
                poisson(A, {monomial: ring.one}),
                denominator_t,
            )
        )
        for monomial in T_SUPPORT
    ]
    rhs = as_field(
        scale(
            pi_power(A, B, 3),
            -ring.domain.one / ring.domain(24),
        )
    )
    base, base_kernel, base_rank = solve_affine(
        base_columns,
        rhs,
        field,
    )
    if base_rank != 42 or base_kernel:
        raise AssertionError((base_rank, len(base_kernel)))
    base_s, base_t = split_solution(base, S_SUPPORT, T_SUPPORT)
    print("lifted universal 42-term hbar^3 solution", flush=True)

    cached_s = {}
    cached_t = {}
    for s_support, t_support in direction_supports.values():
        for monomial in s_support:
            if monomial not in cached_s:
                cached_s[monomial] = as_field(
                    scale(
                        poisson({monomial: ring.one}, B),
                        denominator_s,
                    )
                )
        for monomial in t_support:
            if monomial not in cached_t:
                cached_t[monomial] = as_field(
                    scale(
                        poisson(A, {monomial: ring.one}),
                        denominator_t,
                    )
                )

    directions = {}
    if not full_bounded:
        for index in SELECTED_DIRECTIONS:
            s_support, t_support = direction_supports[index]
            columns = [cached_s[monomial] for monomial in s_support]
            columns += [cached_t[monomial] for monomial in t_support]
            _, kernel, rank = solve_affine(columns, {}, field)
            if rank != len(columns) - 1 or len(kernel) != 1:
                raise AssertionError((index, len(columns), rank, len(kernel)))
            directions[index] = split_solution(
                kernel[0],
                s_support,
                t_support,
            )

    def restrict(poly):
        return [
            poly.get(monomial, field.zero)
            for monomial in functional_support
        ]

    constraints = []
    if full_bounded:
        for s_degree, z_order, left in (
            (21, 1, True),
            (17, 0, False),
        ):
            for monomial in laurent_monomials(
                s_degree,
                z_order,
                0,
                3,
            ):
                image = (
                    as_field(poisson({monomial: ring.one}, B))
                    if left
                    else as_field(poisson(A, {monomial: ring.one}))
                )
                row = restrict(image)
                if any(row):
                    constraints.append(row)
    else:
        constraints = [
            restrict(as_field(poisson({monomial: ring.one}, B)))
            for monomial in D5_S_MONOMIALS
        ]

    if full_bounded:
        full_directions = []
        for group in FULL_DIRECTION_GROUPS:
            s_support = sorted(
                set().union(
                    *(set(direction_supports[index][0]) for index in group)
                )
            )
            t_support = sorted(
                set().union(
                    *(set(direction_supports[index][1]) for index in group)
                )
            )
            columns = [cached_s[monomial] for monomial in s_support]
            columns += [cached_t[monomial] for monomial in t_support]
            _, kernel, rank = solve_affine(columns, {}, field)
            if (
                rank != len(columns) - len(group)
                or len(kernel) != len(group)
            ):
                raise AssertionError(
                    (group, len(columns), rank, len(kernel))
                )
            full_directions.extend(
                split_solution(vector, s_support, t_support)
                for vector in kernel
            )
            print(
                "lifted bounded kernel group",
                group,
                f"({len(columns)} columns, nullity {len(kernel)})",
                flush=True,
            )

        for direction_s, direction_t in full_directions:
            variation = add(
                poisson(direction_s, base_t),
                poisson(base_s, direction_t),
            )
            variation = add(
                variation,
                pi_power(direction_s, T, 3),
                field.one / field(24),
            )
            variation = add(
                variation,
                pi_power(S, direction_t, 3),
                field.one / field(24),
            )
            row = restrict(variation)
            if any(row):
                constraints.append(row)

        for left, right in combinations_with_replacement(
            range(len(full_directions)),
            2,
        ):
            left_s, left_t = full_directions[left]
            right_s, right_t = full_directions[right]
            variation = (
                poisson(left_s, left_t)
                if left == right
                else add(
                    poisson(left_s, right_t),
                    poisson(right_s, left_t),
                )
            )
            row = restrict(variation)
            if any(row):
                constraints.append(row)
    else:
        full_directions = None

    if full_bounded:
        if len(full_directions) != 16:
            raise AssertionError(len(full_directions))
    else:
        # These are the five conditions in the stable pivot subsystem.
        pass

    if not full_bounded:
        for index in (18, 33):
            direction_s, direction_t = directions[index]
            variation = add(
                poisson(direction_s, base_t),
                poisson(base_s, direction_t),
            )
            variation = add(
                variation,
                pi_power(direction_s, T, 3),
                field.one / field(24),
            )
            variation = add(
                variation,
                pi_power(S, direction_t, 3),
                field.one / field(24),
            )
            constraints.append(restrict(variation))
        for left, right in ((14, 22), (18, 33), (22, 40)):
            left_s, left_t = directions[left]
            right_s, right_t = directions[right]
            constraints.append(
                restrict(
                    add(
                        poisson(left_s, right_t),
                        poisson(right_s, left_t),
                    )
                )
            )

    defect = poisson(base_s, base_t)
    defect = add(
        defect,
        pi_power(base_s, T, 3),
        field.one / field(24),
    )
    defect = add(
        defect,
        pi_power(S, base_t, 3),
        field.one / field(24),
    )
    defect = add(
        defect,
        pi_power(S, T, 5),
        field.one / field(1920),
    )
    return constraints, restrict(defect)


def singular_program(
    constraints,
    defect,
    functional_support=SUPPORT,
    full_factors: bool = False,
) -> str:
    # Singular modules store generators as columns.  The 16 generators
    # below are therefore the columns of the 15-by-16 condition matrix.
    generators = []
    for coordinate in range(len(functional_support)):
        entries = [
            to_singular(row[coordinate])
            for row in constraints
        ]
        generators.append("[" + ",".join(entries) + "]")

    lines = [
        "ring r=(0,a,tau),(x),dp;",
        "module M=" + ",".join(generators) + ";",
        "module K=syz(M);",
        'print("KERNEL_GENERATORS="+string(ncols(K)));',
        'print("KERNEL_ROWS="+string(nrows(K)));',
        "vector v=K[1];",
        "number period=0;",
    ]
    for coordinate, value in enumerate(defect, start=1):
        lines.append(
            "period=period+leadcoef(v["
            + str(coordinate)
            + "])*("
            + to_singular(value)
            + ");"
        )
    lines += [
        "number period_numerator=numerator(period);",
        "number period_denominator=denominator(period);",
        'print("PERIOD_DENOMINATOR="+string(period_denominator));',
        "ring polynomial_ring=0,(a,tau),dp;",
        "short=0;",
        "poly P=imap(r,period_numerator);",
        "list factors=factorize(P);",
        "poly large_factor=1;",
        'print("PERIOD_TERMS="+string(size(P)));',
        'print("PERIOD_TOTAL_DEGREE="+string(deg(P)));',
        'print("FACTOR_COUNT="+string(size(factors[1])));',
        "int factor_index;",
        "for(factor_index=1;"
        "factor_index<=size(factors[1]);"
        "factor_index++)",
        "{",
    ]
    if full_factors:
        lines += [
            "  if(size(factors[1][factor_index])>5)",
            "  {",
            "    large_factor=factors[1][factor_index];",
            "  }",
            '  print("FACTOR="+string(factors[1][factor_index])'
            '    +", EXP="+string(factors[2][factor_index])'
            '    +", TERMS="+string(size(factors[1][factor_index]))'
            '    +", DEG="+string(deg(factors[1][factor_index])));',
        ]
    else:
        lines += [
            "  if(size(factors[1][factor_index])<=5)",
            "  {",
            '    print("FACTOR="+string(factors[1][factor_index])'
            '      +", EXP="+string(factors[2][factor_index])'
            '      +", TERMS="+string(size(factors[1][factor_index]))'
            '      +", DEG="+string(deg(factors[1][factor_index])));',
            "  }",
            "  else",
            "  {",
            "    large_factor=factors[1][factor_index];",
            '    print("FACTOR=<omitted>, EXP="'
            '      +string(factors[2][factor_index])'
            '      +", TERMS="+string(size(factors[1][factor_index]))'
            '      +", DEG="+string(deg(factors[1][factor_index])));',
            "  }",
        ]
    lines += [
        "}",
        "poly kappa_one_slice=subst(large_factor,a,-2/3);",
        'print("LARGE_FACTOR_KAPPA_ONE="+string(kappa_one_slice));',
        "quit;",
    ]
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--full-factors",
        action="store_true",
        help="print the 381-term irreducible period factor",
    )
    parser.add_argument(
        "--pivot-only",
        action="store_true",
        help="use only the stable 15-condition subsystem",
    )
    parser.add_argument(
        "--second-support",
        action="store_true",
        help="use the support found on the P_41 finite-field slice",
    )
    parser.add_argument(
        "--third-support",
        action="store_true",
        help="use the support found on V(P_41,Q_32) modulo 32003",
    )
    parser.add_argument(
        "--fourth-support",
        action="store_true",
        help="use the complementary support from V(P_41,Q_32) modulo 32003",
    )
    parser.add_argument(
        "--fifth-support",
        action="store_true",
        help="use the vertical support found at (a,tau)=(-1/2,-3)",
    )
    parser.add_argument(
        "--factor-output",
        help="write the large irreducible factor in Singular syntax",
    )
    args = parser.parse_args()
    singular = shutil.which("Singular")
    if singular is None:
        raise SystemExit("Singular is required for the fraction-free kernel")

    started = perf_counter()
    if sum(
        map(
            int,
            (
                args.second_support,
                args.third_support,
                args.fourth_support,
                args.fifth_support,
            ),
        )
    ) > 1:
        parser.error("choose at most one alternate support")
    functional_support = (
        FIFTH_SUPPORT
        if args.fifth_support
        else FOURTH_SUPPORT
        if args.fourth_support
        else THIRD_SUPPORT
        if args.third_support
        else SECOND_SUPPORT
        if args.second_support
        else SUPPORT
    )
    constraints, defect = construct_system(
        full_bounded=not args.pivot_only,
        functional_support=functional_support,
    )
    print(
        "constructed function-field system:",
        f"{len(constraints)}x{len(SUPPORT)}",
        f"in {perf_counter() - started:.2f}s",
        flush=True,
    )
    program = singular_program(
        constraints,
        defect,
        functional_support=functional_support,
        full_factors=args.full_factors or bool(args.factor_output),
    )
    with tempfile.TemporaryDirectory(
        prefix="degree-five-qper-",
    ) as directory:
        path = Path(directory) / "period.sing"
        path.write_text(program)
        result = subprocess.run(
            [singular, "-q", str(path)],
            check=True,
            capture_output=True,
            text=True,
        )
    print(result.stdout, end="")
    if result.stderr:
        print(result.stderr, end="")
    if args.factor_output:
        large_lines = [
            line
            for line in result.stdout.splitlines()
            if line.startswith("FACTOR=")
            and "TERMS=" in line
            and "TERMS=1," not in line
            and "TERMS=2," not in line
        ]
        if len(large_lines) != 1:
            raise AssertionError(
                f"expected one large factor, got {len(large_lines)}"
            )
        factor_text = large_lines[0].split(", EXP=", 1)[0][len("FACTOR=") :]
        Path(args.factor_output).write_text(factor_text + "\n")


if __name__ == "__main__":
    main()

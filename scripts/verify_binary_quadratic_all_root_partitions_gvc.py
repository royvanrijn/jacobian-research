#!/usr/bin/env python3
"""Verify the complete binary (r,deg(P))=(2,6) GVC row.

The proof is in
``extended-geometry/BINARY_QUADRATIC_ALL_ROOT_PARTITIONS_GVC.md``.
Singular is required for the exact characteristic-zero radicals.
"""

from __future__ import annotations

import shutil
import subprocess

import sympy as sp

from verify_binary_quartic_triple_simple_root_gvc import (
    apply_operator,
    moment,
    singular_expression,
)


ROOT_PARTITIONS = ((2,), (1, 1))


def compositions(total: int, length: int):
    if length == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for tail in compositions(total - first, length - 1):
            yield (first,) + tail


def has_matching(derivatives, annihilators) -> bool:
    def extend(position: int, used: frozenset[int]) -> bool:
        if position == len(derivatives):
            return True
        direction = derivatives[position]
        return any(
            factor not in used
            and annihilator != direction
            and extend(position + 1, used | {factor})
            for factor, annihilator in enumerate(annihilators)
        )

    return extend(0, frozenset())


def verify_hall_locus() -> None:
    for partition in ROOT_PARTITIONS:
        derivatives = tuple(
            direction
            for direction, multiplicity in enumerate(partition)
            for _ in range(multiplicity)
        )
        for counts in compositions(6, len(partition) + 1):
            annihilators = tuple(
                direction
                for direction, count in enumerate(counts[:-1])
                for _ in range(count)
            ) + (-1,) * counts[-1]
            expected_failure = any(
                counts[direction] >= 7 - multiplicity
                for direction, multiplicity in enumerate(partition)
            )
            assert has_matching(derivatives, annihilators) != expected_failure


def radical_equal(
    variables: tuple[sp.Symbol, ...],
    equations: tuple[sp.Expr, ...],
    expected: tuple[sp.Expr, ...],
) -> None:
    """Check equality with the exact radical over Q in Singular."""

    singular = shutil.which("Singular")
    if singular is None:
        raise RuntimeError("Singular is required for the radical replay")
    program = f"""
ring rr=0,({",".join(map(str, variables))}),dp;
ideal I={",".join(map(singular_expression, equations))};
ideal E={",".join(map(singular_expression, expected))};
LIB "primdec.lib";
ideal R=std(radical(I));
ideal L=minbase(reduce(R,std(E)));
ideal T=minbase(reduce(std(E),R));
if ((size(L)==0) && (size(T)==0)) {{ print("PASS"); }}
else {{ print("FAIL"); R; L; T; }}
quit;
"""
    result = subprocess.run(
        [singular, "-q"],
        input=program,
        text=True,
        capture_output=True,
        check=True,
    )
    assert "PASS" in result.stdout and "FAIL" not in result.stdout, (
        result.stdout,
        result.stderr,
    )


def all_moment_coefficients(polynomial, operator, last: int):
    return tuple(
        coefficient
        for order in range(1, last + 1)
        for coefficient in moment(polynomial, operator, order).values()
    )


def verify_distinct_root() -> None:
    """Replay the full first-equation reduction and second-moment ladder."""

    degree = 6
    p = {
        (i, j): sp.symbols(f"distinct_p{i}{j}")
        for i in range(degree + 1)
        for j in range(degree + 1 - i)
    }
    p[(6, 0)] = sp.Integer(1)
    p[(0, 6)] = sp.Integer(0)
    a = {order: sp.symbols(f"distinct_a{order}") for order in range(3, 13)}
    b = {order: sp.symbols(f"distinct_b{order}") for order in range(3, 13)}
    first_operator = {
        (1, 1): 1,
        **{(order, 0): a[order] for order in range(3, 7)},
        **{(0, order): b[order] for order in range(3, 7)},
    }

    # The XY coefficient is a unit, so the first equation solves every
    # mixed coefficient of P triangularly, from high total degree down.
    substitution: dict[sp.Symbol, sp.Expr] = {}
    for total in range(4, -1, -1):
        for output_x in range(total + 1):
            output_y = total - output_x
            target = (output_x + 1, output_y + 1)
            current = {
                exponent: sp.expand(
                    coefficient.subs(substitution)
                    if hasattr(coefficient, "subs")
                    else coefficient
                )
                for exponent, coefficient in p.items()
            }
            equation = apply_operator(current, first_operator).get(
                (output_x, output_y), 0
            )
            target_coefficient = p[target]
            if (
                isinstance(target_coefficient, sp.Symbol)
                and target_coefficient not in substitution
            ):
                substitution[target_coefficient] = sp.factor(
                    sp.solve(equation, target_coefficient)[0]
                )

    reduced_polynomial = {
        exponent: sp.expand(
            coefficient.subs(substitution)
            if hasattr(coefficient, "subs")
            else coefficient
        )
        for exponent, coefficient in p.items()
    }
    assert all(
        sp.factor(coefficient) == 0
        for coefficient in apply_operator(
            reduced_polynomial, first_operator
        ).values()
    )

    full_operator = {
        (1, 1): 1,
        **{(order, 0): a[order] for order in a},
        **{(0, order): b[order] for order in b},
    }
    second = moment(reduced_polynomial, full_operator, 2)

    p05 = p[(0, 5)]
    p04 = p[(0, 4)]
    p03 = p[(0, 3)]
    p02 = p[(0, 2)]
    assert sp.factor(second[(4, 3)]) == 1200 * p05
    assert sp.factor(second[(4, 2)].subs(p05, 0)) == 720 * p04
    assert sp.factor(second[(6, 0)].subs({p05: 0, p04: 0})) == (
        201600 * a[3] ** 2
    )
    first_zeros = {p05: 0, p04: 0, a[3]: 0}
    assert sp.factor(second[(4, 1)].subs(first_zeros)) == 360 * p03
    second_zeros = first_zeros | {p03: 0}
    assert sp.factor(
        second[(4, 0)].subs(second_zeros)
        - 120 * (113040 * a[4] ** 2 + p02)
    ) == 0
    assert sp.factor(
        second[(1, 1)].subs(second_zeros)
        + 2880 * a[4] * (15840 * a[4] ** 2 + p02)
    ) == 0
    third_zeros = second_zeros | {a[4]: 0, p02: 0}
    assert sp.factor(second[(2, 0)].subs(third_zeros)) == (
        213580800 * a[5] ** 2
    )
    fourth_zeros = third_zeros | {a[5]: 0}
    assert sp.factor(second[(0, 0)].subs(fourth_zeros)) == (
        466560000 * a[6] ** 2
    )


def verify_nonpure_double_line() -> None:
    """Check the four complete equality faces over P_6=x*y^5."""

    # The first two half-integral faces have one polynomial companion.
    for suffix, y_order, correction_exponent in (
        ("32", 3, (3, 2)),
        ("52", 5, (3, 0)),
    ):
        a, c = sp.symbols(f"nonpure_half_a{suffix} nonpure_half_c{suffix}")
        operator = {(2, 0): 1, (0, y_order): a}
        polynomial = {(1, 5): 1, correction_exponent: c}
        equations = all_moment_coefficients(polynomial, operator, 6)
        radical_equal((a, c), equations, (a, c))

    # The remaining half-integral crossings are killed by an extremal
    # coefficient before another polynomial companion can enter.
    a = sp.symbols("nonpure_half_a72")
    assert sp.factor(
        moment({(1, 5): 1}, {(2, 0): 1, (0, 7): a}, 2)[(0, 3)]
    ) == 2419200 * a
    a = sp.symbols("nonpure_half_a92")
    assert sp.factor(
        moment({(1, 5): 1}, {(2, 0): 1, (0, 9): a}, 2)[(0, 1)]
    ) == 14515200 * a

    for slope in range(2, 6):
        h, z = sp.symbols(f"nonpure_h{slope} nonpure_z{slope}")
        polynomial_support = [
            (x_degree, y_degree)
            for x_degree in range(7)
            for y_degree in range(7 - x_degree)
            if slope * x_degree + y_degree == slope + 5
            and (x_degree, y_degree) != (1, 5)
        ]
        corrections = sp.symbols(
            f"nonpure_c{slope}_0:{len(polynomial_support)}"
        )
        polynomial = {
            (1, 5): 1,
            **{
                exponent: coefficient
                for exponent, coefficient in zip(
                    polynomial_support, corrections, strict=True
                )
            },
        }
        operator = {(2, 0): 1, (1, slope): h, (0, 2 * slope): z}
        equations = all_moment_coefficients(polynomial, operator, 6)
        radical_equal(
            (h, z, *corrections),
            equations,
            (h, z, *corrections),
        )

    # Record the three terminal coefficients after moment two solves the
    # pure-Y companion on the last three faces.
    h = sp.symbols("nonpure_terminal_h")
    terminal_data = (
        (3, -30, -sp.Rational(2, 7), 2, (1, 1), -3888000),
        (4, -60, -sp.Rational(37, 84), 3, (0, 3), -2193454080000),
        (5, -60, -sp.Rational(83, 168), 3, (0, 0), -15435360000000),
    )
    for slope, correction, companion, order, exponent, coefficient in (
        terminal_data
    ):
        polynomial = {
            (1, 5): 1,
            (2, 5 - slope): correction * h,
        }
        operator = {
            (2, 0): 1,
            (1, slope): h,
            (0, 2 * slope): companion * h**2,
        }
        assert sp.factor(moment(polynomial, operator, order)[exponent]) == (
            coefficient * h**3
        )


def verify_pure_double_line() -> None:
    """Verify every primary and secondary Newton face over P_6=y^6."""

    # Slope 3/2: the pure Y^3 migration and both polynomial companions die.
    v, p, q = sp.symbols("pure_32_v pure_32_p pure_32_q")
    operator = {(2, 0): 1, (0, 3): v}
    polynomial = {(0, 6): 1, (2, 3): p, (4, 0): q}
    equations = all_moment_coefficients(polynomial, operator, 5)
    radical_equal((v, p, q), equations, (v, p, q))

    # Slope 2: exactly two coordinate axes survive.
    A, B, z, q, r = sp.symbols(
        "pure_2_A pure_2_B pure_2_z pure_2_q pure_2_r"
    )
    operator = {(2, 0): 1, (1, 2): B, (0, 4): A}
    polynomial = {(0, 6): 1, (1, 4): z, (2, 2): q, (3, 0): r}
    equations = all_moment_coefficients(polynomial, operator, 8)
    radical_equal((A, B, z, q, r), equations, (r, q, A, B * z))

    # On the B-axis, slope 3 removes the first migration.
    A, z, q = sp.symbols("pure_B3_A pure_B3_z pure_B3_q")
    operator = {(1, 2): 1, (0, 5): A}
    polynomial = {(0, 6): 1, (1, 3): z, (2, 0): q}
    equations = all_moment_coefficients(polynomial, operator, 7)
    radical_equal((A, z, q), equations, (A, z, q))

    # Its slope-4 common-threshold face is also the origin.
    A, z = sp.symbols("pure_B4_A pure_B4_z")
    operator = {(1, 2): 1, (0, 6): A}
    polynomial = {(0, 6): 1, (1, 2): z}
    equations = all_moment_coefficients(polynomial, operator, 7)
    radical_equal((A, z), equations, (A, z))

    # On the z-axis, the first secondary face dies.
    u, v, p = sp.symbols("pure_z_secondary_u pure_z_secondary_v pure_z_p")
    operator = {(2, 0): 1, (1, 3): u, (0, 5): v}
    polynomial = {(0, 6): 1, (1, 4): 1, (2, 1): p}
    equations = all_moment_coefficients(polynomial, operator, 7)
    radical_equal((u, v, p), equations, (u, v, p))

    # The intervening pure Y^6 and Y^7 crossings are triangular.
    a = sp.symbols("pure_z_a6")
    assert sp.factor(
        moment(
            {(0, 6): 1, (1, 4): 1},
            {(2, 0): 1, (0, 6): a},
            1,
        )[(0, 0)]
    ) == 720 * a
    a = sp.symbols("pure_z_a7")
    assert sp.factor(
        moment(
            {(0, 6): 1, (1, 4): 1},
            {(2, 0): 1, (0, 7): a},
            2,
        )[(0, 1)]
    ) == 161280 * a

    # Its final weight-four common-threshold face dies as well.
    u, v, q = sp.symbols("pure_z4_u pure_z4_v pure_z4_q")
    operator = {(2, 0): 1, (1, 4): u, (0, 8): v}
    polynomial = {(1, 4): 1, (2, 0): q}
    equations = all_moment_coefficients(polynomial, operator, 7)
    radical_equal((u, v, q), equations, (u, v, q))

    # On the intersection of the slope-two axes, slope 5/2 dies.
    A, q = sp.symbols("pure_52_A pure_52_q")
    operator = {(2, 0): 1, (0, 5): A}
    polynomial = {(0, 6): 1, (2, 1): q}
    equations = all_moment_coefficients(polynomial, operator, 5)
    radical_equal((A, q), equations, (A, q))

    # The final slope-three face has exactly two one-sided axes.
    A, B, z, q = sp.symbols(
        "pure_3_A pure_3_B pure_3_z pure_3_q"
    )
    operator = {(2, 0): 1, (1, 3): B, (0, 6): A}
    polynomial = {(0, 6): 1, (1, 3): z, (2, 0): q}
    equations = all_moment_coefficients(polynomial, operator, 7)
    radical_equal((A, B, z, q), equations, (q, A, B * z))


def main() -> None:
    verify_hall_locus()
    verify_distinct_root()
    verify_nonpure_double_line()
    verify_pure_double_line()
    print("verified complete binary quadratic-leading sextic GVC row")


if __name__ == "__main__":
    main()

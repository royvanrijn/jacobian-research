#!/usr/bin/env python3
"""Universal canonical-different complex for quartic kernel tensors.

Let mu_ij and s_ij be the trace-free and scalar parts of multiplication in
the generalized cubic algebra.  This checker forms the explicit 4-by-7
matrix

    d1 = [(0,z,-y,x), (s_ij,2*mu_ij)_{0<=i<=j<=2}].

For each squarefree cubic symbol and the full 24-parameter quartic kernel,
it constructs a canonical 7-by-3 matrix d2.  Its six lower rows are the
fixed incidence syzygies obtained by multiplying (z,-y,x) by e_j.  Its
top row is central-quadratic plus parameter-linear cubic.  The checker
verifies d1*d2=0 and the Buchsbaum--Eisenbud grade bounds, so

    0 -> R^3 -> R^7 -> R^4

is exact.  Transposing d2 computes Ext^2 of coker(d1).  The six linear rows
already present a length-six module killed by (x,y,z)^2, so the top row is
redundant and the associated parameter module is free of rank six.

The separate dense-plane tail checker verifies that coker(d1) equals the
actual ramification-support module on seven full-support planes.  Extending
that annihilator/different equality universally is the remaining step; this
checker does not silently identify the two modules.
"""

from __future__ import annotations

import itertools
import multiprocessing
import shutil
import subprocess
from typing import Any

import sympy as sp

import research_universal_cubic_quartic_kernel_saturation as frontier
import verify_cubic_symbol_double_saturation as cubic_audit


PAIRS = tuple(itertools.combinations_with_replacement(range(3), 2))
RELATION_COEFFICIENTS = (
    cubic_audit.z,
    -cubic_audit.y,
    cubic_audit.x,
)


def relation_eigenvalue(
    vector: sp.Matrix,
) -> sp.Expr:
    """Return q from vector=q*(z,-y,x), asserting polynomial divisibility."""

    if vector == sp.zeros(3, 1):
        return sp.Integer(0)
    for entry, coefficient in zip(
        vector, RELATION_COEFFICIENTS
    ):
        if entry == 0:
            continue
        quotient = sp.cancel(entry / coefficient)
        if sp.denom(quotient) == 1:
            result = sp.expand(quotient)
            assert (
                vector - result * cubic_audit.RELATION
            ).applyfunc(sp.expand) == sp.zeros(3, 1)
            return result
    raise AssertionError("relation eigenvalue is not polynomial")


def fixed_linear_tail() -> sp.Matrix:
    """Return the six incidence rows below the varying top row."""

    tail = sp.zeros(6, 3)
    for column in range(3):
        for first, coefficient in enumerate(
            RELATION_COEFFICIENTS
        ):
            pair = tuple(sorted((first, column)))
            tail[PAIRS.index(pair), column] += coefficient
    return tail


def audit_stratum(
    task: tuple[
        str,
        tuple[sp.Symbol, ...],
        dict[tuple[int, int, int], sp.Expr],
    ],
) -> dict[str, Any]:
    """Verify the universal canonical-different complex for one symbol."""

    stratum, parameters, tensor = task
    cubic_audit.FACTOR_SINGULAR_EXPRESSIONS = False
    trace_free_products, scalar_products = (
        cubic_audit.multiplication_table(
            cubic_audit.CUBIC_STRATA[stratum], tensor
        )
    )

    different_columns = [
        sp.Matrix((0, *RELATION_COEFFICIENTS))
    ]
    different_columns.extend(
        sp.Matrix(
            (
                scalar_products[pair],
                *(
                    sp.expand(2 * entry)
                    for entry in trace_free_products[pair]
                ),
            )
        )
        for pair in PAIRS
    )
    first_differential = sp.Matrix.hstack(*different_columns)

    second_differential = sp.zeros(7, 3)
    second_differential[1:, :] = fixed_linear_tail()
    eigenvalues: list[sp.Expr] = []
    variables = (*parameters, *cubic_audit.BASE_VARIABLES)
    bidegrees: set[tuple[int, int]] = set()
    for column in range(3):
        scalar_relation = sp.expand(
            sum(
                RELATION_COEFFICIENTS[first]
                * scalar_products[tuple(sorted((first, column)))]
                for first in range(3)
            )
        )
        assert scalar_relation == 0
        vector_relation = sp.Matrix(
            [
                sp.expand(
                    sum(
                        RELATION_COEFFICIENTS[first]
                        * 2
                        * trace_free_products[
                            tuple(sorted((first, column)))
                        ][row]
                        for first in range(3)
                    )
                )
                for row in range(3)
            ]
        )
        eigenvalue = relation_eigenvalue(vector_relation)
        eigenvalues.append(eigenvalue)
        second_differential[0, column] = -eigenvalue
        for monomial, _coefficient in sp.Poly(
            eigenvalue, *variables
        ).terms():
            parameter_degree = sum(
                monomial[: len(parameters)]
            )
            collision_degree = sum(
                monomial[len(parameters) :]
            )
            assert (parameter_degree, collision_degree) in {
                (0, 2),
                (1, 3),
            }
            bidegrees.add((parameter_degree, collision_degree))

    assert (
        first_differential * second_differential
    ).applyfunc(sp.expand) == sp.zeros(4, 3)

    # I_3(d2) contains x^3, y^3, and z^3 from the fixed lower rows,
    # so its grade is three.  A nonzero central 4-by-4 minor of d1
    # proves generic rank four and grade I_4(d1)>=1.
    linear_tail = fixed_linear_tail()
    maximal_tail_minors = {
        sp.expand(linear_tail[list(rows), :].det())
        for rows in itertools.combinations(range(6), 3)
    }
    assert any(
        sign * cubic_audit.x**3 in maximal_tail_minors
        for sign in (1, -1)
    )
    assert any(
        sign * cubic_audit.y**3 in maximal_tail_minors
        for sign in (1, -1)
    )
    assert any(
        sign * cubic_audit.z**3 in maximal_tail_minors
        for sign in (1, -1)
    )

    central_differential = first_differential.subs(
        {parameter: 0 for parameter in parameters}
    )
    nonzero_minor: tuple[int, ...] | None = None
    for columns in itertools.combinations(range(7), 4):
        if sp.expand(
            central_differential[:, list(columns)].det()
        ) != 0:
            nonzero_minor = columns
            break
    assert nonzero_minor is not None

    return {
        "stratum": stratum,
        "relation_eigenvalue_bidegrees": sorted(bidegrees),
        "nonzero_different_minor_columns": list(nonzero_minor),
        "eigenvalue_term_counts": [
            len(sp.Poly(value, *variables).terms())
            for value in eigenvalues
        ],
    }


def audit_linear_tail() -> None:
    """Check length six, square annihilation, and constant Fittings."""

    linear_tail = fixed_linear_tail()
    coefficient_matrix = sp.zeros(9, 6)
    for relation in range(6):
        for generator in range(3):
            polynomial = sp.Poly(
                linear_tail[relation, generator],
                *cubic_audit.BASE_VARIABLES,
            )
            for variable_index, variable in enumerate(
                cubic_audit.BASE_VARIABLES
            ):
                coefficient_matrix[
                    3 * generator + variable_index, relation
                ] = polynomial.coeff_monomial(variable)
    assert coefficient_matrix.rank() == 6

    singular = shutil.which("Singular")
    assert singular is not None, "Singular is required"
    relation_columns = [
        list(linear_tail[row, :])
        for row in range(linear_tail.rows)
    ]
    program = f"""
ring R=0,(x,y,z),dp;
module linear_tail=
{",".join(map(cubic_audit.singular_vector, relation_columns))};
linear_tail=std(linear_tail);
print("LENGTH="+string(vdim(linear_tail)));
module maximal_square_free=
  (x2+xy+xz+y2+yz+z2)*freemodule(3);
module square_action=simplify(
  reduce(maximal_square_free,linear_tail),2
);
print("SQUARE_ACTION="+string(size(square_action)));
quit;
"""
    result = subprocess.run(
        [singular, "-q"],
        input=program,
        text=True,
        capture_output=True,
        check=True,
    )
    assert "LENGTH=6" in result.stdout
    assert "SQUARE_ACTION=0" in result.stdout


def main() -> None:
    cubic_audit.FACTOR_SINGULAR_EXPRESSIONS = False
    parameters, tensor = frontier.universal_tensor()
    audit_linear_tail()
    tasks = [
        (stratum, parameters, tensor)
        for stratum in sorted(cubic_audit.SQUAREFREE_STRATA)
    ]
    context = multiprocessing.get_context("fork")
    with context.Pool(processes=4) as pool:
        results = list(pool.imap_unordered(audit_stratum, tasks))
    results.sort(key=lambda item: item["stratum"])
    for result in results:
        print(
            "PASS:",
            result["stratum"],
            "universal different complex;",
            "eigenvalue bidegrees",
            result["relation_eigenvalue_bidegrees"],
        )

    print(
        "PASS: every squarefree universal canonical-different complex "
        "is Buchsbaum--Eisenbud exact"
    )
    print(
        "PASS: its Ext^2 is the constant six-dimensional two-layer "
        "module, so the parameter Fittings are Fitt_6=(1), Fitt_5=(0)"
    )
    print(
        "OPEN: identify the universal annihilator Ann(Omega) with the "
        "canonical different module; this equality holds on the "
        "seven checked full-support planes"
    )


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Exact checks for the index-three shallow Segre rigidity theorem.

The proof itself is written in
``extended-geometry/INDEX_THREE_SHALLOW_SEGRE_RIGIDITY.md``.  This script
checks:

* the full-output and rectangular block-cube formulas;
* the integrability identity whose curl is a coefficient commutator;
* a seven-variable, generic-rank-four, exact-index-three calibration;
* an exact polynomial inverse for that calibration.
"""

from __future__ import annotations

import json
from itertools import combinations
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "index_three_shallow_segre_rigidity.json"
)


def full_output_block_check() -> None:
    """Verify formulas (5)--(6) with generic commutative entries."""

    size = 2
    t = sp.symbols("t")
    A = sp.Matrix(size, size, sp.symbols(f"a0:{size * size}"))
    B = sp.Matrix(size, size, sp.symbols(f"b0:{size * size}"))
    q = sp.Matrix(sp.symbols(f"q0:{size}"))
    y = sp.Matrix(sp.symbols(f"y0:{size}"))
    zero = sp.zeros(size)
    M = (t * A).row_join(-t**2 * sp.eye(size)).col_join(
        B.row_join(zero)
    )
    expected = (
        t**3 * (A**3 - A * B - B * A)
    ).row_join(-t**4 * (A**2 - B)).col_join(
        (t**2 * (B * A**2 - B**2)).row_join(-t**3 * B * A)
    )
    assert sp.simplify(M**3 - expected) == sp.zeros(2 * size)

    u = (q - 2 * t * y).col_join(sp.zeros(size, 1))
    expected_column = (
        t**2 * (A**2 - B) * (q - 2 * t * y)
    ).col_join(t * B * A * (q - 2 * t * y))
    assert sp.simplify(M**2 * u - expected_column) == sp.zeros(
        2 * size, 1
    )


def rectangular_block_check() -> None:
    """Verify the rank-compressed residual system (12)."""

    n, rank = 2, 1
    t = sp.symbols("rect_t")
    A = sp.Matrix(n, n, sp.symbols(f"ra0:{n * n}"))
    U = sp.Matrix(n, rank, sp.symbols(f"ru0:{n * rank}"))
    D = sp.Matrix(rank, n, sp.symbols(f"rd0:{rank * n}"))
    q = sp.Matrix(sp.symbols(f"rq0:{n}"))
    y = sp.Matrix(sp.symbols(f"ry0:{rank}"))
    R = U * D
    E = A**2 - R
    M = (t * A).row_join(t**2 * U).col_join(
        (-D).row_join(sp.zeros(rank))
    )
    expected = (
        t**3 * (A**3 - R * A - A * R)
    ).row_join(t**4 * E * U).col_join(
        (-t**2 * D * E).row_join(-t**3 * D * A * U)
    )
    assert sp.simplify(M**3 - expected) == sp.zeros(n + rank)

    u = (q + 2 * t * U * y).col_join(sp.zeros(rank, 1))
    expected_column = (t**2 * E * (q + 2 * t * U * y)).col_join(
        -t * D * A * (q + 2 * t * U * y)
    )
    assert sp.simplify(M**2 * u - expected_column) == sp.zeros(
        n + rank, 1
    )


def integrability_commutator_check() -> None:
    """Check curl(A^2)=[A_i,A_j]x for a generic quadratic map."""

    dimension = 3
    x = sp.Matrix(sp.symbols(f"cx0:{dimension}"))
    q = []
    for output in range(dimension):
        expression = 0
        for left in range(dimension):
            for right in range(left, dimension):
                coefficient = sp.symbols(f"s{output}_{left}_{right}")
                monomial = x[left] * x[right]
                expression += coefficient * monomial
        q.append(expression)
    A = sp.Matrix(q).jacobian(tuple(x))
    coefficient_matrices = [A.diff(variable) for variable in x]
    square = A**2
    for left, right in combinations(range(dimension), 2):
        curl = square.diff(x[left])[:, right] - square.diff(x[right])[:, left]
        commutator = (
            coefficient_matrices[left] * coefficient_matrices[right]
            - coefficient_matrices[right] * coefficient_matrices[left]
        ) * x
        assert sp.simplify(curl - commutator) == sp.zeros(dimension, 1)


def calibration_check() -> dict[str, object]:
    """Verify the seven-variable exact-index-three positive control."""

    x1, x2, x3, y1, y2, y3, t = variables = sp.symbols(
        "x1 x2 x3 y1 y2 y3 t"
    )
    x = sp.Matrix([x1, x2, x3])
    y = sp.Matrix([y1, y2, y3])
    q = sp.Matrix([0, x1**2 / 2, x1 * x2])
    c = sp.Matrix([0, 0, x1**3 / 3])
    h = (t * q - t**2 * y).col_join(c).col_join(sp.Matrix([0]))
    jacobian = h.jacobian(variables)
    assert jacobian**2 != sp.zeros(7)
    assert jacobian**3 == sp.zeros(7)
    assert sp.factor((sp.eye(7) + jacobian).det()) == 1
    assert jacobian.rank() == 4

    # The base scaled map E_t=x+tQ+t^2C has this exact inverse.
    u1, u2, u3, v1, v2, v3, tau = outputs = sp.symbols(
        "u1 u2 u3 v1 v2 v3 tau"
    )
    w1 = u1 + tau**2 * v1
    w2 = u2 + tau**2 * v2
    w3 = u3 + tau**2 * v3
    inverse_x = sp.Matrix(
        [
            w1,
            w2 - tau * w1**2 / 2,
            w3 - tau * w1 * w2 + tau**2 * w1**3 / 6,
        ]
    )
    inverse_y = sp.Matrix(
        [
            v1,
            v2,
            v3 - inverse_x[0] ** 3 / 3,
        ]
    )
    inverse = inverse_x.col_join(inverse_y).col_join(sp.Matrix([tau]))
    forward = sp.Matrix(variables) + h
    substitutions = dict(zip(variables, inverse))
    assert [
        sp.expand(component.subs(substitutions, simultaneous=True))
        for component in forward
    ] == list(outputs)

    return {
        "dimension": 7,
        "generic_rank_JH": 4,
        "nilpotency_index_JH": 3,
        "jacobian_determinant": 1,
        "Q": ["0", "x1^2/2", "x1*x2"],
        "C": ["0", "0", "x1^3/3"],
        "inverse_verified": True,
    }


def main() -> None:
    full_output_block_check()
    rectangular_block_check()
    integrability_commutator_check()
    calibration = calibration_check()
    payload = {
        "format": "index-three-shallow-segre-rigidity-v1",
        "field": "characteristic zero",
        "classified_class": (
            "H(x,y,t)=(t*Q(x)-t^2*y,C(x),0), "
            "Q homogeneous quadratic and C homogeneous cubic"
        ),
        "exact_index_three_equivalence": (
            "(JH)^3=0 iff JC=(JQ)^2 and (JQ)^3=0"
        ),
        "integrability_consequence": (
            "the constant coefficient matrices of JQ commute and are "
            "simultaneously strictly triangular"
        ),
        "classification_conclusion": (
            "every map in the classified index-three class is a polynomial "
            "automorphism; hence the class contains no collision"
        ),
        "collision_gauge": (
            "distinct collision points are linearly independent and may be "
            "sent to e1,e2 by linear conjugacy; translations are not used"
        ),
        "rank_compressed_residual": [
            "E*U=0",
            "D*E=0",
            "D*A*U=0",
            "A^3-R*A-A*R=0",
            "E*Q=0",
            "D*A*Q=0",
            "where R=U*D and E=A^2-R",
        ],
        "calibration": calibration,
        "status": (
            "exhaustive theorem for the displayed shallow circuit class; "
            "not a lower bound for arbitrary cubic-homogeneous maps"
        ),
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    print("PASS shallow Segre: exact full and rectangular block cubes")
    print("PASS shallow Segre: integrability curl equals coefficient commutator")
    print("PASS shallow Segre: 7D calibration has rank 4 and exact index 3")
    print("PASS shallow Segre: exact polynomial inverse verified")
    print(f"PASS shallow Segre: wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

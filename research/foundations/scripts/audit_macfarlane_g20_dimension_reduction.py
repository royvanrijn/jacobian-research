#!/usr/bin/env python3
"""Exact reduction audit for MacFarlane's 20-variable Keller map.

The external map is

    G20(x,w,t) = (x + t R(x) + t^2 B w, w - gamma(x), t),

where its dehomogenization is F13=x+R+B gamma.  This checker does not claim
to construct a 19-variable map.  It certifies the linear obstructions and
the two numerical gates that a successful rank-compressed construction must
beat.
"""

from __future__ import annotations

import hashlib
import json
from collections import defaultdict
from itertools import combinations_with_replacement
from pathlib import Path

import sympy as sp

from rank_compressed_bcw_homogenization import constant_jacobian_kernel


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "macfarlane_g20_dimension_reduction_audit.json"
)
SOURCE_COMMIT = "dad6090bf4f01b3cdad04048fbe16f3be52b485c"
SOURCE_URL = (
    "https://github.com/Amacfa/keller-counterexamples-13-20/tree/"
    + SOURCE_COMMIT
)


def coefficient_matrix(
    expressions: list[sp.Expr], variables: tuple[sp.Symbol, ...]
) -> sp.Matrix:
    polynomials = [sp.Poly(value, *variables, domain=sp.QQ) for value in expressions]
    monomials = sorted(
        {
            exponents
            for polynomial in polynomials
            for exponents, coefficient in polynomial.terms()
            if coefficient
        }
    )
    return sp.Matrix(
        [
            [polynomial.coeff_monomial(exponents) for exponents in monomials]
            for polynomial in polynomials
        ]
    )


def sparse_vector(vector: sp.Matrix) -> list[list[object]]:
    return [
        [index + 1, str(sp.cancel(value))]
        for index, value in enumerate(vector)
        if value
    ]


def build_maps() -> dict[str, object]:
    x = tuple(sp.symbols("x1:14"))
    w = tuple(sp.symbols("w1:7"))
    tau = sp.Symbol("tau")
    (
        x1,
        x2,
        x3,
        x4,
        x5,
        x6,
        x7,
        x8,
        x9,
        x10,
        x11,
        x12,
        x13,
    ) = x

    quadratic = [
        -x11 * x12,
        3 * x1 * x3 - x8 * x9 - 3 * x5 * x6,
        -x8 * x10 + 4 * x2**2 - x4 * x5 - x6 * x7,
        2 * x12 * x13,
        3 * x2**2,
        sp.Integer(0),
        3 * x2 * x3 - x2 * x5,
        x1 * x2,
        6 * x1 * x3 - 3 * x1 * x5 - 3 * x3 * x6,
        -x1 * x7 + 7 * x2**2 - x3 * x4,
        x1 * x3,
        -sp.Rational(1, 2) * x1**2,
        x2**2,
    ]
    gamma = [
        -2 * x1 * x3 * x8
        + x1 * x5 * x8
        - sp.Rational(1, 3) * x1 * x2 * x9
        + 4 * x1 * x2**2
        + x3 * x6 * x8
        - 3 * x2**2 * x6,
        x1 * x7 * x8
        - x1 * x2 * x10
        - 7 * x2**2 * x8
        + x3 * x4 * x8
        - 3 * x2**2 * x4
        - 3 * x2 * x3 * x6
        + x2 * x5 * x6,
        x1**2 * x13 - 2 * x12 * x2**2,
        -sp.Rational(1, 2) * x1**2 * x11 + x1 * x12 * x3,
        x1 * x2 * x3,
        x1**2 * x2,
    ]
    B = sp.zeros(13, 6)
    B[0, 3] = -1
    B[0, 5] = -sp.Rational(3, 2)
    B[1, 0] = 3
    B[2, 1] = 1
    B[2, 4] = 3
    B[3, 2] = -1
    B[4, 4] = 1
    B[5, 5] = 1

    cubic = list(B * sp.Matrix(gamma))
    nonlinear13 = [sp.expand(q + c) for q, c in zip(quadratic, cubic)]
    f13 = [sp.expand(variable + correction) for variable, correction in zip(x, nonlinear13)]
    variables20 = x + w + (tau,)
    h20 = [
        sp.expand(tau * quadratic[index] + tau**2 * (B * sp.Matrix(w))[index])
        for index in range(13)
    ] + [-value for value in gamma] + [sp.Integer(0)]
    g20 = [
        sp.expand(variable + correction)
        for variable, correction in zip(variables20, h20)
    ]

    p13 = sp.Matrix([0, 0, -sp.Rational(1, 4)] + [0] * 10)
    q13 = sp.Matrix(
        [
            1,
            -sp.Rational(3, 2),
            sp.Rational(13, 2),
            -sp.Rational(9, 4),
            3,
            sp.Rational(3, 2),
            sp.Rational(99, 4),
            sp.Rational(3, 2),
            -sp.Rational(3, 4),
            -sp.Rational(45, 8),
            -sp.Rational(13, 2),
            sp.Rational(1, 2),
            -sp.Rational(9, 4),
        ]
    )
    gamma_q = sp.Matrix(gamma).subs(dict(zip(x, q13)))
    p20 = p13.col_join(sp.zeros(6, 1)).col_join(sp.Matrix([1]))
    q20 = q13.col_join(gamma_q).col_join(sp.Matrix([1]))

    return {
        "x": x,
        "w": w,
        "tau": tau,
        "variables20": variables20,
        "R": quadratic,
        "gamma": gamma,
        "B": B,
        "C": cubic,
        "N13": nonlinear13,
        "F13": f13,
        "H20": h20,
        "G20": g20,
        "p13": p13,
        "q13": q13,
        "p20": p20,
        "q20": q20,
    }


def fixed_space_through_degree_three(
    f13: list[sp.Expr], variables: tuple[sp.Symbol, ...]
) -> dict[str, object]:
    # The exact torus grading is recovered from the component monomials.
    nonlinear = [sp.expand(value - variables[index]) for index, value in enumerate(f13)]
    grading_rows: list[list[int]] = []
    for index, expression in enumerate(nonlinear):
        for exponents, coefficient in sp.Poly(
            expression, *variables, domain=sp.QQ
        ).terms():
            if coefficient:
                row = list(exponents)
                row[index] -= 1
                grading_rows.append(row)
    grading_matrix = sp.Matrix(grading_rows)
    grading_kernel = grading_matrix.nullspace()
    assert len(grading_kernel) == 1
    grading = grading_kernel[0]
    denominator = sp.ilcm(*[value.q for value in grading])
    weights = [int(value * denominator) for value in grading]
    divisor = sp.igcd(*weights)
    weights = [value // divisor for value in weights]
    assert weights == [-1, 1, 2, 0, 2, -1, 3, 0, 1, 2, 1, -2, 2]

    sectors: dict[int, list[sp.Expr]] = defaultdict(list)
    for degree in range(4):
        for indices in combinations_with_replacement(range(len(variables)), degree):
            monomial = sp.Mul(*(variables[index] for index in indices))
            sectors[sum(weights[index] for index in indices)].append(monomial)

    substitution = dict(zip(variables, f13))
    prime = 1_000_003

    def residue(value: sp.Expr) -> int:
        rational = sp.Rational(value)
        return (
            int(rational.p) % prime
        ) * pow(int(rational.q) % prime, prime - 2, prime) % prime

    def sparse_modular_rank(
        columns: list[dict[tuple[int, ...], sp.Expr]]
    ) -> tuple[int, list[tuple[int, ...]]]:
        pivots: dict[tuple[int, ...], dict[tuple[int, ...], int]] = {}
        for rational_column in columns:
            column = {
                monomial: coefficient
                for monomial, value in rational_column.items()
                if (coefficient := residue(value))
            }
            while column:
                pivot = min(column)
                if pivot not in pivots:
                    inverse = pow(column[pivot], prime - 2, prime)
                    column = {
                        monomial: coefficient * inverse % prime
                        for monomial, coefficient in column.items()
                    }
                    pivots[pivot] = column
                    break
                factor = column[pivot]
                for monomial, coefficient in pivots[pivot].items():
                    updated = (column.get(monomial, 0) - factor * coefficient) % prime
                    if updated:
                        column[monomial] = updated
                    else:
                        column.pop(monomial, None)
        return len(pivots), sorted(pivots)

    sector_records: list[dict[str, int]] = []
    total_nullity = 0
    for weight, monomials in sorted(sectors.items()):
        columns: list[dict[tuple[int, ...], sp.Expr]] = []
        row_monomials: set[tuple[int, ...]] = set()
        for monomial in monomials:
            defect = sp.Poly(
                sp.expand(monomial.subs(substitution, simultaneous=True) - monomial),
                *variables,
                domain=sp.QQ,
            )
            column = dict(defect.terms())
            columns.append(column)
            row_monomials.update(column)
        rank, pivots = sparse_modular_rank(columns)
        nullity = len(monomials) - rank
        total_nullity += nullity
        sector_records.append(
            {
                "weight": weight,
                "columns": len(monomials),
                "rows": len(row_monomials),
                "rank_mod_1000003": rank,
                "nullity": nullity,
                "pivot_count": len(pivots),
            }
        )
        assert nullity == (1 if weight == 0 else 0)
    # Reduction modulo a good prime has full column rank after removing the
    # literal constant column.  The same nonzero minors are nonzero over QQ.
    assert total_nullity == 1
    return {
        "torus_weights": weights,
        "rank_certificate_prime": prime,
        "sector_records": sector_records,
        "fixed_basis": ["1"],
    }


def main() -> None:
    data = build_maps()
    x = data["x"]
    variables20 = data["variables20"]
    quadratic = data["R"]
    B = data["B"]
    cubic = data["C"]
    nonlinear13 = data["N13"]
    f13 = data["F13"]
    h20 = data["H20"]
    g20 = data["G20"]
    p13 = data["p13"]
    q13 = data["q13"]
    p20 = data["p20"]
    q20 = data["q20"]

    # Reproduce the two collisions from the external certificate.
    f13_matrix = sp.Matrix(f13)
    g20_matrix = sp.Matrix(g20)
    assert f13_matrix.subs(dict(zip(x, p13))) == p13
    assert f13_matrix.subs(dict(zip(x, q13))) == p13
    assert g20_matrix.subs(dict(zip(variables20, p20))) == p20
    assert g20_matrix.subs(dict(zip(variables20, q20))) == p20
    assert p13 != q13 and p20 != q20
    assert all(
        not polynomial.is_zero and polynomial.total_degree() == 3
        for polynomial in [
            sp.Poly(value, *variables20, domain=sp.QQ)
            for value in h20
            if value
        ]
    )

    # The two rank-compressed gates: 13+6+1=20 at present.
    cubic_output_rank = coefficient_matrix(cubic, x).rank()
    assert cubic_output_rank == 6
    assert B.rank() == 6

    h20_polys = tuple(sp.Poly(value, *variables20, domain=sp.QQ) for value in h20)
    n13_polys = tuple(sp.Poly(value, *x, domain=sp.QQ) for value in nonlinear13)
    kernel20 = constant_jacobian_kernel(h20_polys, variables20)
    kernel13 = constant_jacobian_kernel(n13_polys, x)
    assert kernel20.cols == 0
    assert kernel13.cols == 0

    h20_coefficients = coefficient_matrix(h20, variables20)
    n13_coefficients = coefficient_matrix(nonlinear13, x)
    assert h20_coefficients.rank() == 19
    assert n13_coefficients.rank() == 13
    h20_left_kernel = h20_coefficients.T.nullspace()
    assert len(h20_left_kernel) == 1
    assert h20_left_kernel[0] == sp.eye(20)[:, 19]

    # Exact linear algebra used in the hyperplane no-go proofs.
    r_coefficients = coefficient_matrix(quadratic, x)
    r_left_kernel = r_coefficients.T.nullspace()
    expected_r_kernel = [
        sp.eye(13)[:, 5],
        -sp.Rational(1, 3) * sp.eye(13)[:, 4] + sp.eye(13)[:, 12],
    ]
    assert sp.Matrix.hstack(*r_left_kernel).columnspace() == sp.Matrix.hstack(
        *expected_r_kernel
    ).columnspace()
    assert sp.Matrix.vstack(r_coefficients.T, B.T).nullspace() == []

    # If a.R is a nonzero square (a.x)^2, square-monomial support forces
    # a into span(e1,e2); B then kills neither direction.
    square_variables = sorted(
        {
            index
            for expression in quadratic
            for exponents, coefficient in sp.Poly(
                expression, *x, domain=sp.QQ
            ).terms()
            for index, exponent in enumerate(exponents)
            if coefficient and exponent == 2
        }
    )
    assert square_variables == [0, 1]
    assert B[:2, :].rank() == 2

    # For F13, the collision conditions cut ker(a.R) to the displayed line.
    collision_rows = sp.Matrix.vstack(p13.T, q13.T)
    restricted_r_kernel = sp.Matrix.hstack(*r_left_kernel)
    collision_kernel = (collision_rows * restricted_r_kernel).nullspace()
    assert len(collision_kernel) == 1
    collision_covector = restricted_r_kernel * collision_kernel[0]
    collision_covector *= sp.ilcm(*[value.q for value in collision_covector])
    divisor = sp.igcd(*[int(value) for value in collision_covector])
    collision_covector /= divisor
    assert collision_covector == sp.Matrix(
        [0, 0, 0, 0, -2, 13, 0, 0, 0, 0, 0, 0, 6]
    )
    residual_cubic = sp.expand(collision_covector.dot(sp.Matrix(cubic)))
    assert residual_cubic == sp.expand(
        sp.Rational(13, 1) * x[0] ** 2 * x[1]
        - 2 * x[0] * x[1] * x[2]
    )

    artifact = {
        "format": "macfarlane-g20-dimension-reduction-audit-v1",
        "external_source": SOURCE_URL,
        "status": (
            "exact obstruction audit; no 19-variable counterexample is constructed"
        ),
        "verified_external_structure": {
            "F13_dimension": 13,
            "G20_dimension": 20,
            "cubic_output_rank_of_F13": cubic_output_rank,
            "rank_compressed_count": "13+6+1=20",
            "F13_collision": True,
            "G20_collision": True,
            "G20_cubic_homogeneous": True,
        },
        "linear_audit": {
            "G20_constant_JH_kernel_dimension": kernel20.cols,
            "G20_H_output_span_dimension": h20_coefficients.rank(),
            "G20_fixed_linear_covectors": [
                sparse_vector(vector) for vector in h20_left_kernel
            ],
            "F13_constant_nonlinear_Jacobian_kernel_dimension": kernel13.cols,
            "F13_nonlinear_output_span_dimension": n13_coefficients.rank(),
            "quadratic_R_left_kernel": [
                sparse_vector(vector) for vector in expected_r_kernel
            ],
            "ker_left_R_intersection_ker_left_B_dimension": 0,
            "G20_Keller_linear_hyperplane_restrictions_through_collision": 0,
            "F13_Keller_affine_hyperplane_restrictions_through_collision": 0,
        },
        "dimension_19_success_gates": [
            "a collision-preserving degree-three source in dimension 12 with cubic-output rank at most 6",
            "a 13-variable stable-equivalent source with cubic-output rank at most 5",
            "a genuinely nonlinear quotient or restriction outside the audited linear and degree-at-most-three invariant mechanisms",
        ],
    }
    OUTPUT.write_text(json.dumps(artifact, indent=2) + "\n")
    digest = hashlib.sha256(OUTPUT.read_bytes()).hexdigest()

    print("PASS MacFarlane G20: reproduced exact 13D and 20D collisions")
    print("PASS MacFarlane G20: cubic-output rank is 6, giving 13+6+1=20")
    print("PASS linear audit: G20 and F13 have zero constant input kernel")
    print("PASS linear audit: no Keller hyperplane restriction through the collision")
    print("OPEN dimension 19: the two rank-compressed construction gates remain")
    print(f"PASS wrote {OUTPUT.relative_to(ROOT)}")
    print(f"SHA256 {digest}")


if __name__ == "__main__":
    main()

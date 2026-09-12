#!/usr/bin/env python3
"""Verify the corner-derived F2 d=3 modified-chart obstruction.

This checker separates two assertions which were previously conflated.

* The F2 complete chain itself fixes the selected root multiplicity
  ``gamma=2``.  The monomial chart

      X = xi*y^2,   Y = xi^(-2)*y^(-7)

  then sends the certified bracket ``X^4`` to ``-3*xi^2*y^2`` and sends
  the terminal vertices to x-degrees ``3*r`` and ``3*(2*r-1)``.  The
  corner/source envelope gives finite coefficient supports in this chart.

* Passing from the full Laurent pullback to a polynomial modified pair is
  still a cutting/normal-form statement.  This checker does not assert it.

The support box is only an envelope: coefficients on one translated source
band satisfy exact binomial-jet relations.  The checker constructs their
linear image, rather than treating the displayed support positions as
independent variables.  In the literal nonnegative-xi projection, the r=3
top diagonal then has an exact unit-ideal obstruction.  This kills every
branch of that *conditional polynomial projection*, but does not supply the
missing Laurent-to-polynomial cutting theorem.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import hashlib
import json
from pathlib import Path
import shutil
import subprocess

import sympy as sp


ROOT = Path(__file__).resolve().parents[2]
ARTIFACT = (
    ROOT
    / "artifacts/generated-results/jc2_f2_modified_chart_bridge.json"
)


def ceil_div(numerator: int, denominator: int) -> int:
    return -((-numerator) // denominator)


def source_band_parameters(
    degree: int,
    terminal_height: int,
    r: int,
    layer: int,
) -> tuple[tuple[int, ...], int, int]:
    """Return source exponents, first surviving exponent, and jet order."""

    source_i_min = max(0, ceil_div(layer, 5))
    source_i_max = (degree + layer) // 6
    if source_i_max < source_i_min:
        return (), 0, 0
    source_j = tuple(
        5 * source_i - layer
        for source_i in range(source_i_min, source_i_max + 1)
    )
    denominator = 5 * r - 3
    t_min = max(
        layer,
        ceil_div((7 * r - 4) * layer - terminal_height, denominator),
    )
    jet_order = t_min - layer
    if jet_order >= len(source_j):
        return (), t_min, jet_order
    return source_j, t_min, jet_order


def source_band_support(
    degree: int,
    terminal_height: int,
    r: int,
    layer: int,
) -> tuple[int, ...]:
    """Return the *possible* active t exponents on one translated band.

    An original monomial x^i*y^j has layer ``ell=5*i-j`` and becomes
    ``t^ell*(1+t)^j*z^ell``.  The total-degree and terminal-halfspace
    inequalities give the source interval and the initial jet order below.

    Every returned exponent can occur, but exponents on the same band are
    generally linked.  ``exact_source_projection`` retains those linear
    binomial-jet relations.
    """

    source_j, t_min, _jet_order = source_band_parameters(
        degree, terminal_height, r, layer
    )
    if not source_j:
        return ()

    assert len(set(source_j)) == len(source_j)
    return tuple(range(t_min, layer + max(source_j) + 1))


def pre_shift_supports(r: int, side: str) -> dict[int, set[int]]:
    m = 2 * r - 1
    if side == "P":
        degree = 25 * r
        height = r
        expected_x_degree = 3 * r
        maximum_layer = 5 * r
    elif side == "Q":
        degree = 25 * m
        height = m
        expected_x_degree = 3 * m
        maximum_layer = 5 * m
    else:
        raise ValueError(side)

    supports: dict[int, set[int]] = defaultdict(set)
    # Every source band lies in [-degree,5*degree].  Empty or terminal-killed
    # bands are discarded by source_band_support.
    # The translated type-II edge is the second supporting halfspace:
    # layer=X_power-Y_power<=5*r (respectively 5*m).
    for layer in range(-degree, maximum_layer + 1):
        for t_exponent in source_band_support(
            degree, height, r, layer
        ):
            # t=xi^(-1)y^(-5), z=xi^2*y^7.
            x_exponent = 2 * layer - t_exponent
            y_exponent = 7 * layer - 5 * t_exponent
            if x_exponent >= 0:
                supports[x_exponent].add(y_exponent)

    assert max(supports) == expected_x_degree
    assert supports[expected_x_degree] == {0}
    assert supports[expected_x_degree - 1] == {-2, -5}
    return dict(supports)


def power_support(exponents: set[int], power: int) -> set[int]:
    if power == 0:
        return {0}
    result = {0}
    for _ in range(power):
        result = {left + right for left in result for right in exponents}
    return result


def post_shift_supports(
    pre_shift: dict[int, set[int]], degree: int
) -> dict[int, list[int]]:
    """Return the coefficient envelope after xi=x-G and p_1=0.

    The leading-minus-one coefficient proves
    ``Supp(G) subset {-2,-5}``.  Coefficient index i means x^(degree-i).
    """

    g_support = {-2, -5}
    result: dict[int, list[int]] = {}
    for index in range(2, degree + 1):
        target_x = degree - index
        support: set[int] = set()
        for old_x in range(target_x, degree + 1):
            for old_y in pre_shift.get(old_x, set()):
                support.update(
                    old_y + shift_y
                    for shift_y in power_support(
                        g_support, old_x - target_x
                    )
                )
        result[index] = sorted(support, reverse=True)
    return result


def exact_source_projection(
    r: int, side: str
) -> tuple[list[tuple[int, int]], sp.Matrix]:
    """Return the exact nonnegative-xi binomial-jet image over QQ.

    Columns are a basis of source coefficients satisfying the killed initial
    jets.  Rows are the possible modified-chart monomials.  Different layer
    blocks have disjoint rows because the monomial bridge is invertible.
    """

    m = 2 * r - 1
    if side == "P":
        degree, height, maximum_layer = 25 * r, r, 5 * r
    elif side == "Q":
        degree, height, maximum_layer = 25 * m, m, 5 * m
    else:
        raise ValueError(side)

    positions: list[tuple[int, int]] = []
    blocks: list[sp.Matrix] = []
    for layer in range(-degree, maximum_layer + 1):
        source_j, t_min, jet_order = source_band_parameters(
            degree, height, r, layer
        )
        if not source_j:
            continue
        if jet_order:
            killed = sp.Matrix(
                [
                    [sp.binomial(j, order) for j in source_j]
                    for order in range(jet_order)
                ]
            )
            kernel_vectors = killed.nullspace()
            if not kernel_vectors:
                continue
            kernel = sp.Matrix.hstack(*kernel_vectors)
        else:
            kernel = sp.eye(len(source_j))

        t_exponents = [
            exponent
            for exponent in range(t_min, layer + max(source_j) + 1)
            if 2 * layer - exponent >= 0
        ]
        block = sp.Matrix(
            [
                [
                    sum(
                        kernel[source_index, column]
                        * sp.binomial(
                            source_j[source_index], exponent - layer
                        )
                        for source_index in range(len(source_j))
                    )
                    for column in range(kernel.cols)
                ]
                for exponent in t_exponents
            ]
        )
        nonzero_rows = [
            row
            for row in range(block.rows)
            if any(block[row, column] for column in range(block.cols))
        ]
        if not nonzero_rows:
            continue
        block = block[nonzero_rows, :]
        t_exponents = [t_exponents[row] for row in nonzero_rows]
        positions.extend(
            (2 * layer - exponent, 7 * layer - 5 * exponent)
            for exponent in t_exponents
        )
        blocks.append(block)

    assert len(positions) == len(set(positions))
    return positions, sp.diag(*blocks)


def top_band_projection_matrix(n: int) -> sp.Matrix:
    """Return the (3*n+1)-by-(2*n+1) top-band source image."""

    source_j = tuple(range(0, 20 * n + 1, 5))
    killed = sp.Matrix(
        [
            [sp.binomial(j, order) for j in source_j]
            for order in range(2 * n)
        ]
    )
    kernel = sp.Matrix.hstack(*killed.nullspace())
    matrix = sp.Matrix(
        [
            [
                sum(
                    kernel[source_index, column]
                    * sp.binomial(source_j[source_index], 2 * n + index)
                    for source_index in range(len(source_j))
                )
                for column in range(kernel.cols)
            ]
            for index in range(3 * n + 1)
        ]
    )
    assert matrix.shape == (3 * n + 1, 2 * n + 1)
    assert matrix.rank() == 2 * n + 1
    return matrix


def full_laurent_seed_audit() -> dict[str, object]:
    """Verify the exact two-parameter full common-root top band."""

    t, u, rho1, rho2 = sp.symbols("t u rho1 rho2")
    v = (1 + t) ** 5 - 1
    phi = sp.cancel(v / t)
    r_polynomial = sp.Rational(1, 25) + rho1 * v + rho2 * v**2
    h_polynomial = sp.expand(phi**2 * r_polynomial)
    a0 = sp.expand(t**2 * h_polynomial)
    assert sp.expand(a0 - v**2 * r_polynomial) == 0
    assert sp.Poly(h_polynomial, t).degree() == 18
    assert h_polynomial.subs(t, 0) == 1
    assert sp.Poly(h_polynomial, t).LC() == rho2

    a0_in_u = sp.Poly(sp.expand(a0.subs(t, u - 1)), u)
    assert all(exponent[0] % 5 == 0 for exponent, _ in a0_in_u.terms())
    assert sp.expand(a0_in_u.as_expr().subs(u, 1)) == 0

    # First term lost by the positive projection.  If A=H^r and
    # B=trunc_(<=3r)(A), comparison with the identically zero Wronskian of
    # H^r,H^m gives
    #   [s^(3r)](m*B'*D-r*B*D')=-m*(3r+1)*[s^(3r+1)]H^r.
    # Verify the r=3 instance directly in the two-parameter source seed.
    h_coefficients = [h_polynomial.coeff(t, index) for index in range(11)]
    p_power = truncated_power(h_coefficients, 3, 11)
    q_power = truncated_power(h_coefficients, 5, 11)
    positive_p = sum(p_power[index] * t**index for index in range(10))
    positive_q = sum(q_power[index] * t**index for index in range(11))
    projected_wronskian = sp.expand(
        5 * sp.diff(positive_p, t) * positive_q
        - 3 * positive_p * sp.diff(positive_q, t)
    )
    first_tail_coefficient = sp.expand(p_power[10])
    assert first_tail_coefficient != 0
    assert sp.expand(
        projected_wronskian.coeff(t, 9) + 50 * first_tail_coefficient
    ) == 0

    # The complete tail is the image of an honest source polynomial band.
    # Here y_source is the coordinate before the Puiseux translation.
    x_source, y_source = sp.symbols("x_source y_source")
    source_v = x_source * y_source**5 - 1
    source_root = sp.expand(
        x_source
        * source_v**2
        * (sp.Rational(1, 25) + rho1 * source_v + rho2 * source_v**2)
    )
    assert sp.Poly(source_root, x_source, y_source).total_degree() == 25
    assert sp.Poly(source_root, x_source, y_source).coeff_monomial(
        x_source**5 * y_source**20
    ) == rho2

    return {
        "common_root": "C_top=xi^3*H(s), s=xi^-1*y^-5",
        "Phi": "((1+s)^5-1)/s",
        "v": "(1+s)^5-1",
        "H": "Phi(s)^2*(1/25+rho1*v+rho2*v^2)",
        "free_top_parameters": ["rho1", "rho2"],
        "endpoint_condition": "rho2!=0",
        "endpoint_monomial": "rho2*xi^-15*y^-90",
        "source_invariance": (
            "s^2*H(s)=v^2*(1/25+rho1*v+rho2*v^2) belongs to QQ[(1+s)^5]"
        ),
        "exact_polynomial_source_lift": (
            "x*(x*y^5-1)^2*(1/25+rho1*(x*y^5-1)+"
            "rho2*(x*y^5-1)^2)"
        ),
        "corrected_top_tangent": (
            "the exact source bands do not force the former divisible "
            "common-root slice; for r=3 the first-five kernel dimensions "
            "are 6,6,7,7,10, and the extra descent-5 mode is the commuting "
            "C0^4 term"
        ),
        "lambda_location": (
            "C0^-1 is the descent-30/layer-10 formal resonance; it is not "
            "an independent homogeneous source-band kernel"
        ),
        "first_positive_projection_tail_correction": (
            "[s^(3*r)]((2*r-1)*B'*D-r*B*D')="
            "-(2*r-1)*(3*r+1)*[s^(3*r+1)]H^r"
        ),
        "r3_first_correction_is_nonzero_polynomial": True,
        "claim_boundary": (
            "this is the exact full-Laurent top seed, not a finite plane-map "
            "candidate; all lower positive/negative bands remain to be glued"
        ),
    }


def source_projection_audit() -> tuple[
    dict[str, object], list[tuple[int, int]], sp.Matrix
]:
    """Audit exact r=3 projection ranks and terminal non-membership."""

    p_positions, p_matrix = exact_source_projection(3, "P")
    q_positions, q_matrix = exact_source_projection(3, "Q")
    assert len(p_positions) == 83
    assert p_matrix.rank() == 74
    assert len(q_positions) == 215
    assert q_matrix.rank() == 196

    top = top_band_projection_matrix(3)
    left_kernel = top.T.nullspace()
    assert len(left_kernel) == 3
    relation_vector = left_kernel[0]
    denominator = sp.ilcm(*(entry.q for entry in relation_vector))
    relation = [int(entry * denominator) for entry in relation_vector]
    common_divisor = sp.igcd(*relation)
    relation = [entry // common_divisor for entry in relation]
    if relation[0] < 0:
        relation = [-entry for entry in relation]
    expected_relation = [
        3470256,
        -1118621,
        299936,
        -65511,
        11250,
        -1430,
        120,
        -5,
        0,
        0,
    ]
    assert relation == expected_relation
    assert (sp.Matrix(relation).T * top).is_zero_matrix

    # P_T=x^9+a*x^2*y has only its leading monomial on this diagonal.
    terminal_top = sp.Matrix([1] + [0] * 9)
    terminal_evaluation = int((sp.Matrix(relation).T * terminal_top)[0])
    assert terminal_evaluation == 3470256

    payload = {
        "r3_nonnegative_projection": {
            "P_support_positions": len(p_positions),
            "P_source_image_rank": p_matrix.rank(),
            "P_linear_relations": len(p_positions) - p_matrix.rank(),
            "Q_support_positions": len(q_positions),
            "Q_source_image_rank": q_matrix.rank(),
            "Q_linear_relations": len(q_positions) - q_matrix.rank(),
        },
        "P_top_diagonal": {
            "positions": [
                [3 * 3 - index, -5 * index] for index in range(10)
            ],
            "source_image_shape": list(top.shape),
            "source_image_rank": top.rank(),
            "primitive_relation": relation,
            "terminal_resonance_evaluation": terminal_evaluation,
            "terminal_resonance_in_source_image": False,
        },
        "claim_boundary": (
            "support positions are an envelope, not independent scalar "
            "coordinates; these ranks retain all exact binomial-jet links"
        ),
    }
    return payload, p_positions, p_matrix


def truncated_product(
    left: list[sp.Expr], right: list[sp.Expr], order: int
) -> list[sp.Expr]:
    """Multiply coefficient lists modulo v**order."""

    result = [sp.S.Zero] * order
    for left_index, left_coefficient in enumerate(left[:order]):
        if left_coefficient == 0:
            continue
        for right_index, right_coefficient in enumerate(
            right[: order - left_index]
        ):
            if right_coefficient != 0:
                result[left_index + right_index] += (
                    left_coefficient * right_coefficient
                )
    return [sp.expand(coefficient) for coefficient in result]


def truncated_power(
    coefficients: list[sp.Expr], exponent: int, order: int
) -> list[sp.Expr]:
    """Raise a coefficient list to a nonnegative power modulo v**order."""

    result = [sp.S.One] + [sp.S.Zero] * (order - 1)
    for _ in range(exponent):
        result = truncated_product(result, coefficients, order)
    return result


def top_gap_equations(
    n: int,
    coefficients: tuple[sp.Symbol, sp.Symbol, sp.Symbol],
) -> list[sp.Expr]:
    """Return the n exact source gaps for a common cubic raised to n.

    With v=(1+t)^5-1, tau=(1+v)^(1/5)-1, and

        L(v)=H(tau)*(tau/v)^2,

    top-band membership of H(t)^n is equivalent to the vanishing of
    [v^q]L(v)^n for 2*n < q <= 3*n.
    """

    c1, c2, c3 = coefficients
    v = sp.symbols("v")
    order = 3 * n + 1
    tau = sp.series(
        (1 + v) ** sp.Rational(1, 5) - 1,
        v,
        0,
        order + 2,
    ).removeO()
    h_tau = 1 + c1 * tau + c2 * tau**2 + c3 * tau**3
    l_series = sp.series(
        h_tau * (tau / v) ** 2, v, 0, order
    ).removeO().expand()
    l_coefficients = [l_series.coeff(v, index) for index in range(order)]
    power = truncated_power(l_coefficients, n, order)
    return [
        sp.expand(power[index]) for index in range(2 * n + 1, 3 * n + 1)
    ]


def top_diagonal_gap_audit() -> dict[str, object]:
    """Prove the r=3 projected top-diagonal ideal is the unit ideal."""

    c1, c2, c3 = sp.symbols("c1 c2 c3")
    variables = (c1, c2, c3)

    r2_equations = top_gap_equations(2, variables)
    m3_equations = top_gap_equations(3, variables)
    r2_combined = sp.groebner(
        [*r2_equations, *m3_equations],
        *variables,
        order="lex",
        method="f5b",
    )
    assert [polynomial.as_expr() for polynomial in r2_combined.polys] == [1]

    r3_equations = top_gap_equations(3, variables)
    # Verify directly that the v-gap equations are just an invertible linear
    # re-basing of the three binomial-source relations evaluated on H(t)^3.
    t = sp.symbols("t")
    h_polynomial = 1 + c1 * t + c2 * t**2 + c3 * t**3
    h_cube_coefficients = sp.Matrix(
        [
            sp.expand(h_polynomial**3).coeff(t, index)
            for index in range(10)
        ]
    )
    source_relations = [
        sp.expand((relation.T * h_cube_coefficients)[0])
        for relation in top_band_projection_matrix(3).T.nullspace()
    ]
    monomials = sorted(
        set().union(
            *(
                set(sp.Poly(equation, *variables).monoms())
                for equation in [*source_relations, *r3_equations]
            )
        )
    )

    def coefficient_matrix(equations: list[sp.Expr]) -> sp.Matrix:
        return sp.Matrix(
            [
                [
                    sp.Poly(equation, *variables).coeff_monomial(monomial)
                    for equation in equations
                ]
                for monomial in monomials
            ]
        )

    source_relation_matrix = coefficient_matrix(source_relations)
    gap_relation_matrix = coefficient_matrix(r3_equations)
    assert source_relation_matrix.rank() == 3
    assert gap_relation_matrix.rank() == 3
    assert source_relation_matrix.row_join(gap_relation_matrix).rank() == 3
    change_matrix = source_relation_matrix.gauss_jordan_solve(
        gap_relation_matrix
    )[0]
    assert sp.factor(change_matrix.det()) == sp.Rational(1, 5**42)

    r3_basis = sp.groebner(
        r3_equations, *variables, order="lex", method="f5b"
    )
    leading_exponents = [
        polynomial.LM(order=r3_basis.order).exponents
        for polynomial in r3_basis.polys
    ]
    assert leading_exponents == [(1, 0, 0), (0, 1, 0), (0, 0, 27)]
    artinian_length = 27

    # For r=3, m=5.  The very first Q gap is already a unit in the
    # 27-dimensional P-gap algebra.  In triangular coordinates its residue
    # is a univariate polynomial of degree 26, coprime to the degree-27
    # eliminant.  Their resultant is the multiplication determinant/norm.
    q_first_gap = top_gap_equations(5, variables)[0]
    q_remainder = sp.cancel(r3_basis.reduce(q_first_gap)[1])
    eliminant = next(
        polynomial.as_expr()
        for polynomial in r3_basis.polys
        if polynomial.as_expr().free_symbols <= {c3}
    )
    assert q_remainder.free_symbols <= {c3}
    assert sp.degree(eliminant, c3) == 27
    assert sp.degree(q_remainder, c3) == 26
    gcd = sp.gcd(sp.Poly(eliminant, c3), sp.Poly(q_remainder, c3))
    assert gcd.degree() == 0
    resultant = sp.cancel(sp.resultant(eliminant, q_remainder, c3))
    assert resultant != 0
    numerator, denominator = sp.fraction(resultant)
    resultant_text = f"{numerator}/{denominator}"
    resultant_digest = hashlib.sha256(resultant_text.encode()).hexdigest()
    assert resultant_digest == (
        "cd610d23ca92bc10c7788b72a603d3f71353d36f107bb87695669f4decb27314"
    )

    # At r=4, P membership alone is already inconsistent.  This is an exact
    # family-dependence sample; it is not promoted to an all-r induction.
    r4_basis = sp.groebner(
        top_gap_equations(4, variables),
        *variables,
        order="grevlex",
        method="f5b",
    )
    assert [polynomial.as_expr() for polynomial in r4_basis.polys] == [1]

    return {
        "derivation": {
            "top_band_membership": (
                "B_n(t)=Phi(t)^(2*n)*R_n((1+t)^5) mod t^(3*n+1), "
                "deg(R_n)<=2*n"
            ),
            "Phi": "((1+t)^5-1)/t",
            "common_power": (
                "top bracket zero implies m*B'*D-r*B*D'=0; "
                "gcd(r,2*r-1)=1 gives B=H^r, D=H^(2*r-1), deg(H)<=3"
            ),
            "gap_function": (
                "L(v)=H((1+v)^(1/5)-1)*(((1+v)^(1/5)-1)/v)^2"
            ),
            "gap_ideal": (
                "I_r=([v^q]L^r, q=2*r+1..3*r; "
                "[v^q]L^(2*r-1), q=4*r-1..6*r-3)"
            ),
        },
        "r2": {
            "equation_count": 5,
            "groebner_basis": ["1"],
            "verdict": "projected top diagonal empty",
        },
        "r3": {
            "P_gap_equation_count": 3,
            "source_relation_to_gap_change_determinant": "5^-42",
            "P_gap_leading_ideal": ["c1", "c2", "c3^27"],
            "P_gap_artinian_length": artinian_length,
            "first_Q_gap": "[v^11]L(v)^5",
            "first_Q_gap_residue_degree": 26,
            "P_eliminant_degree": 27,
            "residue_gcd_degree": gcd.degree(),
            "multiplication_determinant_nonzero": True,
            "resultant_sign": 1 if resultant > 0 else -1,
            "resultant_numerator_digits": len(str(abs(int(numerator)))),
            "resultant_denominator_digits": len(str(abs(int(denominator)))),
            "resultant_sha256": resultant_digest,
            "combined_gap_ideal": "unit ideal",
            "verdict": "all literal polynomial-projection branches empty",
        },
        "r4": {
            "P_gap_equation_count": 4,
            "P_gap_groebner_basis": ["1"],
            "verdict": "P projection alone empty",
        },
        "claim_boundary": (
            "the obstruction is unconditional for the literal nonnegative-"
            "xi source projection; applying it to an F2 Keller map still "
            "requires a cutting theorem showing that this projection is a "
            "polynomial bracket-preserving modified pair"
        ),
    }


def singular_top_gap_unit_audit(n: int) -> None:
    """Run an optional exact QQ top-gap sample with Singular."""

    singular = shutil.which("Singular")
    if singular is None:
        raise RuntimeError("Singular is required for the extended family sample")
    c1, c2, c3 = sp.symbols("c1 c2 c3")
    equations = top_gap_equations(n, (c1, c2, c3))
    primitive_equations: list[sp.Expr] = []
    for equation in equations:
        _, integral = sp.Poly(equation, c1, c2, c3).clear_denoms(
            convert=True
        )
        _, primitive = integral.primitive()
        primitive_equations.append(primitive.as_expr())

    def singular_expression(expression: sp.Expr) -> str:
        return sp.sstr(expression).replace("**", "^")

    lines = ["ring R=0,(c1,c2,c3),dp;", "option(redSB);"]
    lines.extend(
        f"poly f{index}={singular_expression(equation)};"
        for index, equation in enumerate(primitive_equations, start=1)
    )
    generators = ",".join(
        f"f{index}" for index in range(1, len(primitive_equations) + 1)
    )
    lines.extend(
        [
            f"ideal I={generators};",
            "ideal G=std(I);",
            'if ((size(G)==1) && (G[1]==1)) { print("UNIT"); }',
        ]
    )
    completed = subprocess.run(
        [singular, "-q"],
        input="\n".join(lines) + "\n",
        text=True,
        capture_output=True,
        check=True,
        timeout=120,
    )
    if "UNIT" not in completed.stdout:
        raise AssertionError(
            f"the n={n} Singular top-gap basis was not the unit ideal"
        )


def monomial_bracket(
    left: list[tuple[sp.Expr, sp.Expr, sp.Expr]],
    right: list[tuple[sp.Expr, sp.Expr, sp.Expr]],
) -> dict[tuple[sp.Expr, sp.Expr], sp.Expr]:
    result: dict[tuple[sp.Expr, sp.Expr], sp.Expr] = defaultdict(
        lambda: sp.Integer(0)
    )
    for left_x, left_y, left_coefficient in left:
        for right_x, right_y, right_coefficient in right:
            exponent = (
                sp.expand(left_x + right_x - 1),
                sp.expand(left_y + right_y - 1),
            )
            result[exponent] += sp.expand(
                left_coefficient
                * right_coefficient
                * (left_x * right_y - left_y * right_x)
            )
    return {
        exponent: sp.factor(coefficient)
        for exponent, coefficient in result.items()
        if sp.factor(coefficient) != 0
    }


def terminal_resonance_audit() -> dict[str, object]:
    r = sp.symbols("r", integer=True, positive=True)
    m = 2 * r - 1
    a = sp.symbols("a", nonzero=True)
    binomial_two = sp.factor(m * (r - 1) / (2 * r**2))
    terminal_p = [(3 * r, 0, 1), (2, 1, a)]
    terminal_q = [
        (3 * m, 0, 1),
        (3 * r - 1, 1, m * a / r),
        (1, 2, binomial_two * a**2),
    ]
    bracket = monomial_bracket(terminal_p, terminal_q)
    expected_coefficient = sp.factor(
        3 * a**3 * (r - 1) * m / (2 * r**2)
    )
    assert bracket == {(sp.Integer(2), sp.Integer(2)): expected_coefficient}

    resonance_index = 3 * r - 2
    leading_q_degree = 3 * m
    assert sp.expand(leading_q_degree - 2 * resonance_index) == 1
    assert sp.expand(leading_q_degree - 3 * resonance_index) == 3 - 3 * r
    assert sp.expand(3 * r + (3 - 3 * r) - 1) == 2

    # The first omitted binomial coefficient is cancelled by F.  This also
    # checks the bracket coefficient independently through 3*r*F'.
    alpha = sp.factor(m / r)
    omitted = sp.factor(alpha * (alpha - 1) * (alpha - 2) / 6)
    f_coefficient = sp.factor(-omitted)
    assert sp.factor(
        3 * r * 3 * f_coefficient * a**3 - expected_coefficient
    ) == 0

    # At r=3 this is B=1+a*y*t^7.  Its formal cube root has t-support in
    # multiples of seven, so its t^18 coefficient (x^-15 in C) vanishes.
    assert resonance_index.subs(r, 3) == 7
    assert 18 % 7 != 0

    return {
        "family_pair": {
            "P": "x^(3*r)+a*x^2*y",
            "Q": (
                "x^(3*(2*r-1))+((2*r-1)/r)*a*x^(3*r-1)*y+"
                "((2*r-1)*(r-1)/(2*r^2))*a^2*x*y^2"
            ),
        },
        "bracket": (
            "3*(r-1)*(2*r-1)*a^3/(2*r^2) * x^2*y^2"
        ),
        "polynomial_coordinate": (
            "B(t)=1+a*y*t^(3*r-2); Q is the polynomial part of "
            "P^((2*r-1)/r)"
        ),
        "first_omitted_x_exponent": "3-3*r",
        "first_visible_F_coefficient": (
            "-binom((2*r-1)/r,3)*a^3*y^3"
        ),
        "r3_specialization": {
            "B": "1+a*y*t^7",
            "P": "x^9+a*x^2*y",
            "Q": "x^15+(5*a/3)*x^8*y+(5*a^2/9)*x*y^2",
            "bracket": "5*a^3*x^2*y^2/3",
            "common_top_test": (
                "the formal cube root is supported in t-degrees 7*k, "
                "so [t^18]B^(1/3)=0"
            ),
        },
        "claim_boundary": (
            "this is the exact formal terminal block in the ambient support "
            "box; it is outside the exact nonnegative-xi source image, omits "
            "the upper common-root endpoint, and is not a Keller counterexample"
        ),
    }


def tangent_space_audit(
    p_supports: dict[int, list[int]],
    source_positions: list[tuple[int, int]],
    source_matrix: sp.Matrix,
) -> dict[str, object]:
    """Compute the ambient support-box tangent space at the terminal point."""

    X, z = sp.symbols("X z")
    p0 = X**9 + z**2 * X**2
    q0 = (
        X**15
        + sp.Rational(5, 3) * z**2 * X**8
        + sp.Rational(5, 9) * z**4 * X
    )

    def bracket_operator(left: sp.Expr, right: sp.Expr) -> sp.Expr:
        return 3 * (
            5 * sp.diff(left, X) * right
            - 3 * left * sp.diff(right, X)
            - z
            * (
                sp.diff(left, X) * sp.diff(right, z)
                - sp.diff(left, z) * sp.diff(right, X)
            )
        )

    def linearized(delta_p: sp.Expr, delta_q: sp.Expr) -> sp.Expr:
        return sp.expand(
            bracket_operator(delta_p, q0)
            + bracket_operator(p0, delta_q)
        )

    names: list[str] = []
    columns: list[dict[tuple[int, int], sp.Expr]] = []
    row_keys: set[tuple[int, int]] = set()

    for index, y_exponents in sorted(p_supports.items()):
        x_exponent = 9 - index
        for y_exponent in y_exponents:
            assert (index - y_exponent) % 3 == 0
            z_exponent = (index - y_exponent) // 3
            delta_p = X**x_exponent * z**z_exponent
            # d(pol_X(P^(5/3)))=(5/3)*pol_X(P^(2/3)*dP).
            # At P0, only the first two binomial terms can have nonnegative
            # X-degree for these supports.
            delta_q = (
                sp.Rational(5, 3)
                * X ** (x_exponent + 6)
                * z**z_exponent
            )
            if x_exponent >= 1:
                delta_q += (
                    sp.Rational(10, 9)
                    * X ** (x_exponent - 1)
                    * z ** (z_exponent + 2)
                )
            polynomial = sp.Poly(linearized(delta_p, delta_q), X, z)
            column = {
                (int(monomial[0]), int(monomial[1])): coefficient
                for monomial, coefficient in polynomial.terms()
                if coefficient != 0
            }
            names.append(f"a{index}_{z_exponent}")
            columns.append(column)
            row_keys.update(column)

    # Linearize E-mu*z^6*(X-g2*z-g5*z^2)^2 at
    extra_columns = [
        ("mu", {(2, 6): sp.Integer(-1)}),
        ("g2", {(1, 7): sp.Rational(10, 3)}),
        ("g5", {(1, 8): sp.Rational(10, 3)}),
    ]
    for name, column in extra_columns:
        names.append(name)
        columns.append(column)
        row_keys.update(column)

    ordered_rows = sorted(row_keys)
    row_index = {key: index for index, key in enumerate(ordered_rows)}
    matrix = sp.zeros(len(ordered_rows), len(columns))
    for column_index, column in enumerate(columns):
        for key, coefficient in column.items():
            matrix[row_index[key], column_index] = coefficient

    rank = matrix.rank()
    nullspace = matrix.nullspace()
    assert len(names) == 83  # 80 P coefficients plus mu,g2,g5.
    assert rank == 77
    assert len(nullspace) == 6
    modes = [
        [
            [names[index], str(sp.factor(value))]
            for index, value in enumerate(vector)
            if value != 0
        ]
        for vector in nullspace
    ]
    assert modes == [
        [["a4_2", "1"]],
        [["a7_3", "1"]],
        [["a2_1", "6"], ["a9_3", "1"]],
        [["a7_2", "1/5"], ["mu", "1"]],
        [["a8_3", "-2"], ["g2", "1"]],
        [["a8_4", "1"], ["g5", "1"]],
    ]

    # The 80 support-box coordinates are not independent source directions.
    # Pull the six ambient modes back across xi=x-G.  At first order,
    # delta(p_pre)=delta(p_post)+(g2*z+g5*z^2)*d_X(p0).
    pre_mode_polynomials = [
        X**5 * z**2,
        X**2 * z**3,
        6 * X**7 * z + z**3,
        sp.Rational(1, 5) * X**2 * z**2,
        -2 * X * z**3 + z * sp.diff(p0, X),
        X * z**4 + z**2 * sp.diff(p0, X),
    ]
    mode_terms: list[dict[tuple[int, int], sp.Expr]] = []
    outside_positions: set[tuple[int, int]] = set()
    source_position_set = set(source_positions)
    for polynomial in pre_mode_polynomials:
        terms: dict[tuple[int, int], sp.Expr] = {}
        for (x_exponent, z_exponent), coefficient in sp.Poly(
            sp.expand(polynomial), X, z
        ).terms():
            position = (
                int(x_exponent),
                int(9 - x_exponent - 3 * z_exponent),
            )
            terms[position] = coefficient
            if position not in source_position_set:
                outside_positions.add(position)
        mode_terms.append(terms)

    all_positions = [*source_positions, *sorted(outside_positions)]
    extended_source = sp.zeros(len(all_positions), source_matrix.cols)
    extended_source[: len(source_positions), :] = source_matrix
    position_index = {
        position: index for index, position in enumerate(all_positions)
    }
    mode_matrix = sp.zeros(len(all_positions), len(pre_mode_polynomials))
    for mode_index, terms in enumerate(mode_terms):
        for position, coefficient in terms.items():
            mode_matrix[position_index[position], mode_index] = coefficient
    source_relations = extended_source.T.nullspace()
    compatibility = sp.Matrix(
        [
            [
                (relation.T * mode_matrix[:, mode_index])[0]
                for mode_index in range(mode_matrix.cols)
            ]
            for relation in source_relations
        ]
    )
    compatible_modes = compatibility.nullspace()
    expected_compatible_modes = [
        sp.Matrix([1, 0, 0, 0, 0, 0]),
        sp.Matrix([0, 1, 0, 0, 0, 0]),
        sp.Matrix([0, 0, 0, 1, 0, 0]),
    ]
    assert compatible_modes == expected_compatible_modes

    deepest_names = [f"a{index}_{2 * index}" for index in range(2, 10)]
    for name in deepest_names:
        coordinate = names.index(name)
        assert all(vector[coordinate] == 0 for vector in nullspace)

    canonical_matrix = [
        [
            [int(entry.p), int(entry.q)]
            if isinstance(entry, sp.Rational)
            else [int(entry), 1]
            for entry in matrix.row(row)
        ]
        for row in range(matrix.rows)
    ]
    matrix_digest = hashlib.sha256(
        json.dumps(canonical_matrix, separators=(",", ":")).encode()
    ).hexdigest()
    return {
        "equation": (
            "3*(5*p_X*q-3*p*q_X-z*(p_X*q_z-p_z*q_X))="
            "mu*z^6*(X-g2*z-g5*z^2)^2"
        ),
        "coordinates": (
            "x=y*X, z=y^-3, P=y^9*p(X,z), Q=y^15*q(X,z)"
        ),
        "support_box_coefficient_positions": 80,
        "support_box_coordinates_are_independent_source_variables": False,
        "target_parameters": ["mu", "g2", "g5"],
        "matrix_shape": [matrix.rows, matrix.cols],
        "matrix_rank": rank,
        "tangent_dimension": len(nullspace),
        "nullspace_modes": modes,
        "ambient_modes_lying_in_exact_source_linear_image": [
            "a4_2",
            "a7_3",
            "(1/5)*a7_2+mu",
        ],
        "source_image_compatible_mode_dimension": len(compatible_modes),
        "deepest_common_top_directions": deepest_names,
        "deepest_directions_in_tangent_space": False,
        "matrix_sha256": matrix_digest,
        "claim_boundary": (
            "the six-dimensional calculation is the tangent space in the "
            "80-position support-box over-approximation.  The terminal base "
            "point itself is outside the exact source image, so the three "
            "listed source-image directions are not tangent branches based "
            "at a source-compatible point"
        ),
    }


def build_payload() -> dict[str, object]:
    r_symbol = sp.symbols("r", integer=True, positive=True)
    m_symbol = 2 * r_symbol - 1

    # The complete-chain generated-corner formula is
    # A1=A0'+gamma*(1/5,1).  Its second coordinate fixes gamma=2.
    gamma = sp.solve(
        [sp.Eq(1 + sp.Symbol("gamma") / 5, sp.Rational(7, 5)),
         sp.Eq(sp.Symbol("gamma"), 2)],
        [sp.Symbol("gamma")],
        dict=True,
    )
    assert gamma == [{sp.Symbol("gamma"): 2}]

    X, Y, xi, y = sp.symbols("X Y xi y", nonzero=True)
    X_sub = xi * y**2
    Y_sub = xi**-2 * y**-7
    coordinate_jacobian = sp.factor(
        sp.diff(X_sub, xi) * sp.diff(Y_sub, y)
        - sp.diff(X_sub, y) * sp.diff(Y_sub, xi)
    )
    assert coordinate_jacobian == -3 * xi**-2 * y**-6
    assert sp.factor(X_sub**4 * coordinate_jacobian) == -3 * xi**2 * y**2

    # Check the general terminal vertices under (a,b)->(a-2b,2a-7b).
    exponent_map = lambda a, b: (
        sp.expand(a - 2 * b),
        sp.expand(2 * a - 7 * b),
    )
    assert exponent_map(7 * r_symbol, 2 * r_symbol) == (3 * r_symbol, 0)
    assert exponent_map(4, 1) == (2, 1)
    assert exponent_map(7 * m_symbol, 2 * m_symbol) == (3 * m_symbol, 0)
    assert exponent_map(7 * r_symbol - 3, m_symbol) == (3 * r_symbol - 1, 1)
    assert exponent_map(1, 0) == (1, 2)

    # The common root t^7*H(t)*z^5 becomes xi^3*H(xi^-1*y^-5).
    t_sub = sp.factor(X_sub * Y_sub)
    z_sub = sp.factor(Y_sub**-1)
    assert t_sub == 1 / (xi * y**5)
    assert z_sub == xi**2 * y**7
    h = sp.Function("H")
    common_root = sp.factor(t_sub**7 * h(t_sub) * z_sub**5)
    assert common_root == xi**3 * h(1 / (xi * y**5))

    pre_p = pre_shift_supports(3, "P")
    pre_q = pre_shift_supports(3, "Q")
    p_masks = post_shift_supports(pre_p, 9)
    q_masks = post_shift_supports(pre_q, 15)
    expected_p_masks = {
        2: [-1, -4, -7, -10],
        3: [0, -3, -6, -9, -12, -15],
        4: [-2, -5, -8, -11, -14, -17, -20],
        5: [-1, -4, -7, -10, -13, -16, -19, -22, -25],
        6: [0, -3, -6, -9, -12, -15, -18, -21, -24, -27, -30],
        7: [1, -2, -5, -8, -11, -14, -17, -20, -23, -26, -29, -32, -35],
        8: [-1, -4, -7, -10, -13, -16, -19, -22, -25, -28, -31, -34, -37, -40],
        9: [0, -3, -6, -9, -12, -15, -18, -21, -24, -27, -30, -33, -36, -39, -42, -45],
    }
    assert p_masks == expected_p_masks
    assert sum(map(len, p_masks.values())) == 80
    assert sum(map(len, q_masks.values())) == 212

    # Sample enough r to guard the symbolic family pattern of the leading
    # chart degree and G support.
    family_samples = []
    for r in range(2, 9):
        for side in ("P", "Q"):
            support = pre_shift_supports(r, side)
            expected_degree = 3 * (r if side == "P" else 2 * r - 1)
            assert max(support) == expected_degree
            assert support[expected_degree - 1] == {-2, -5}
        family_samples.append(
            {
                "r": r,
                "P_x_degree": 3 * r,
                "Q_x_degree": 3 * (2 * r - 1),
                "leading_minus_one_support": [-2, -5],
            }
        )

    projection, p_source_positions, p_source_matrix = source_projection_audit()
    full_laurent_seed = full_laurent_seed_audit()
    tangent = tangent_space_audit(
        p_masks, p_source_positions, p_source_matrix
    )
    top_diagonal = top_diagonal_gap_audit()
    terminal = terminal_resonance_audit()
    return {
        "schema": "plane-jc.f2-modified-chart-bridge.v3",
        "status": "exact-source-projection-obstruction-conditional-on-cut",
        "forced_chain_data": {
            "A0": [5, 20],
            "A0_prime": [1, 0],
            "A1": ["7/5", 2],
            "generated_corner_formula": "A1=A0'+gamma*(1/5,1)",
            "gamma": 2,
            "historical_d": 3,
            "consequence": (
                "the d=2 historical branch is not an F2 branch; it belongs "
                "to a different generated corner/orientation"
            ),
        },
        "monomial_chart": {
            "substitution": {"X": "xi*y^2", "Y": "xi^-2*y^-7"},
            "exponent_map": "(a,b)->(a-2*b,2*a-7*b)",
            "coordinate_jacobian": "-3*xi^-2*y^-6",
            "bracket": "X^4 -> -3*xi^2*y^2",
            "terminal_degrees": {
                "P": "3*r",
                "Q": "3*(2*r-1)",
            },
            "common_root": "xi^3*H(xi^-1*y^-5)",
            "tschirnhaus_shift_support": ["y^-2", "y^-5"],
            "shift": "xi=x-G, G=g_-2*y^-2+g_-5*y^-5",
        },
        "support_derivation": {
            "source_constraints": (
                "ell=5*i-j, i>=0, j>=0, i+j<=D"
            ),
            "terminal_t_bound": (
                "k>=ceil(((7*r-4)*ell-height)/(5*r-3))"
            ),
            "bridge_exponents": "x=2*ell-k, y=7*ell-5*k",
            "binomial_jet_statement": (
                "the killed initial jets are imposed before projecting the "
                "remaining binomial-span image; positions on one band are "
                "not treated as independent"
            ),
            "r3_pre_shift_leading": {
                "P_x9": [0],
                "P_x8": [-2, -5],
                "Q_x15": [0],
                "Q_x14": [-2, -5],
            },
            "r3_P_polynomial_masks": {
                f"p{index}": exponents
                for index, exponents in p_masks.items()
            },
            "r3_Q_mask_summary": {
                f"q{index}": {
                    "maximum": max(exponents),
                    "minimum": min(exponents),
                    "count": len(exponents),
                }
                for index, exponents in q_masks.items()
            },
            "r3_P_support_box_positions_after_shift": 80,
            "r3_Q_support_box_positions_after_shift": 212,
            "family_samples": family_samples,
            "claim_boundary": (
                "these are exact coefficient envelopes for the nonnegative "
                "xi projection; the corner chain alone does not prove that "
                "discarding the negative xi tail preserves the bracket"
            ),
        },
        "exact_source_projection": projection,
        "full_laurent_residual_seed": full_laurent_seed,
        "obsolete_sharp_candidate": {
            "old_candidate_bracket": "(35/9)*y^6*x^2",
            "forced_F2_chart_bracket_weight": "y^2*(x-G)^2",
            "verdict": (
                "excluded from the forced gamma=2 chart by its leading "
                "y-weight"
            ),
        },
        "terminal_resonance": terminal,
        "ambient_terminal_tangent_space": tangent,
        "projected_top_diagonal_obstruction": top_diagonal,
        "remaining_bridge": (
            "prove that the literal nonnegative-xi source projection is the "
            "polynomial bracket-preserving modified pair (or derive the "
            "precise negative-tail correction).  The r=3 projected residue "
            "itself is the unit ideal, so no lower projected branch remains"
        ),
        "software": {"sympy": sp.__version__, "arithmetic": "exact QQ"},
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--refresh", action="store_true")
    parser.add_argument(
        "--extended-r5",
        action="store_true",
        help="also prove the r=5 P-only top-gap ideal is (1) with Singular",
    )
    parser.add_argument("--artifact", type=Path, default=ARTIFACT)
    args = parser.parse_args()

    payload = build_payload()
    if args.extended_r5:
        singular_top_gap_unit_audit(5)
    artifact = args.artifact.resolve()
    if args.refresh:
        artifact.parent.mkdir(parents=True, exist_ok=True)
        artifact.write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        try:
            display = artifact.relative_to(ROOT)
        except ValueError:
            display = artifact
        print(f"WROTE {display}")
    else:
        expected = json.loads(artifact.read_text(encoding="utf-8"))
        current_claim = {
            key: value for key, value in payload.items() if key != "software"
        }
        pinned_claim = {
            key: value for key, value in expected.items() if key != "software"
        }
        if current_claim != pinned_claim:
            raise AssertionError(
                "pinned F2 modified-chart bridge artifact is stale; "
                "inspect before --refresh"
            )

    print("F2_MODIFIED_GAMMA2_BRIDGE_PASS")
    print("F2_MODIFIED_CORNER_SUPPORT_MASKS_PASS")
    print("F2_MODIFIED_EXACT_SOURCE_PROJECTION_RANK_PASS")
    print("F2_MODIFIED_TERMINAL_OUTSIDE_SOURCE_IMAGE_PASS")
    print("F2_MODIFIED_TERMINAL_RESONANCE_PASS")
    print("F2_MODIFIED_TERMINAL_TANGENT_RANK_PASS")
    print("F2_R3_PROJECTED_TOP_GAP_ARTINIAN_LENGTH=27")
    print("F2_R3_PROJECTED_TOP_GAP_RESULTANT_NONZERO_PASS")
    print("F2_R3_PROJECTED_TOP_GAP_UNIT_IDEAL_PASS")
    if args.extended_r5:
        print("F2_R5_PROJECTED_P_TOP_GAP_UNIT_IDEAL_PASS")


if __name__ == "__main__":
    main()

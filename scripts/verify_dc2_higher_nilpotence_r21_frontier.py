#!/usr/bin/env python3
"""Exact higher-nilpotence family and reciprocal-R21 admission frontier.

The theorem part constructs an all-degree Hamiltonian Hessian family with
``N^4=0`` and ``N^2 != 0`` whose Cayley transforms are integrable. It also
proves that this family is a triangular symplectic automorphism control and
is Moyal-flat. A differential calculation explains why its surrounding
two-variable triangular ansatz cannot make ``N^2`` nonconstant.

The frontier part reconstructs the type-(2,1) reciprocal cancellation map
over ``QQ[q]/(q^2-4q+6)``, verifies its Keller and cotangent-graph identities,
checks the natural unimodular stable charts supplied by ``1+x*y^2``, and
performs two exact low-complexity admission tests.  No nonzero constant
target two-form has constant pullback, excluding the affine-contact lift,
while one elementary tangent-normalized shear removes the complete linear
graph defect.  Exact tame corrections continue through degree four, but the
reducible factorization of R excludes every fiber-preserving completion.
The necessary stable-mixed U_2 chart is then normalized through degree six
by twenty-six exact factors, leaving the first defect in degree seven.  An
Euler-homotopy recurrence proves that every subsequent homogeneous defect
can be removed in the completed local ring.  This is formal Darboux
triviality, not polynomial R21 rank-two admission.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import sys
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from master_cancellation import fiber_antiderivative, reduce_q  # noqa: E402


q1, q2, p1, p2 = DARBOUX_VARIABLES = sp.symbols("q1 q2 p1 p2")
parameter_c = sp.symbols("c")

# Bracket convention {p_i,q_j}=delta_ij.
POISSON = sp.Matrix(
    (
        (0, 0, -1, 0),
        (0, 0, 0, -1),
        (1, 0, 0, 0),
        (0, 1, 0, 0),
    )
)


def is_zero_matrix(matrix: sp.MatrixBase) -> bool:
    return all(sp.expand(entry) == 0 for entry in matrix)


def rows_are_closed(jacobian: sp.MatrixBase) -> bool:
    return all(
        sp.expand(
            sp.diff(jacobian[row, left], DARBOUX_VARIABLES[right])
            - sp.diff(jacobian[row, right], DARBOUX_VARIABLES[left])
        )
        == 0
        for row in range(4)
        for left in range(4)
        for right in range(left + 1, 4)
    )


def pi_power(left: sp.Expr, right: sp.Expr, power: int) -> sp.Expr:
    """Apply the constant Poisson bidifferential ``power`` times."""

    terms = tuple(
        (i, j, POISSON[i, j])
        for i in range(4)
        for j in range(4)
        if POISSON[i, j]
    )
    zero = (0, 0, 0, 0)
    states = {(zero, zero): sp.Integer(1)}
    for _ in range(power):
        next_states: dict[
            tuple[tuple[int, ...], tuple[int, ...]], sp.Expr
        ] = {}
        for (left_orders, right_orders), coefficient in states.items():
            for left_index, right_index, sign in terms:
                next_left = list(left_orders)
                next_right = list(right_orders)
                next_left[left_index] += 1
                next_right[right_index] += 1
                key = (tuple(next_left), tuple(next_right))
                next_states[key] = next_states.get(key, 0) + sign * coefficient
        states = next_states

    answer = 0
    for (left_orders, right_orders), coefficient in states.items():
        left_derivative = left
        right_derivative = right
        for variable, order in zip(DARBOUX_VARIABLES, left_orders, strict=True):
            if order:
                left_derivative = sp.diff(left_derivative, variable, order)
        for variable, order in zip(DARBOUX_VARIABLES, right_orders, strict=True):
            if order:
                right_derivative = sp.diff(right_derivative, variable, order)
        answer += coefficient * left_derivative * right_derivative
    return sp.expand(answer)


def verify_higher_nilpotence_family() -> dict[str, object]:
    """Verify the all-degree triangular regular-[4] Cayley family."""

    g_coefficients = sp.symbols("g0:5")
    h_coefficients = sp.symbols("h0:6")
    g = sum(value * p1**degree for degree, value in enumerate(g_coefficients))
    h = sum(value * p1**degree for degree, value in enumerate(h_coefficients))
    potential = sp.expand(
        p1 * q2 + parameter_c * p2**2 / 2 + g * p2 + h
    )
    nilpotent = POISSON * sp.hessian(potential, DARBOUX_VARIABLES)
    square = (nilpotent * nilpotent).applyfunc(sp.expand)
    cube = (square * nilpotent).applyfunc(sp.expand)
    fourth = (square * square).applyfunc(sp.expand)
    expected_square = sp.Matrix(
        (
            (0, 0, 0, parameter_c),
            (0, 0, -parameter_c, 0),
            (0, 0, 0, 0),
            (0, 0, 0, 0),
        )
    )
    expected_cube = sp.Matrix(
        (
            (0, 0, parameter_c, 0),
            (0, 0, 0, 0),
            (0, 0, 0, 0),
            (0, 0, 0, 0),
        )
    )
    assert square == expected_square
    assert cube == expected_cube
    assert is_zero_matrix(fourth)

    cayley = sp.eye(4) + nilpotent + square / 2 + cube / 4
    assert rows_are_closed(cayley)
    outputs = (
        sp.expand(
            q1
            - q2
            - sp.diff(g, p1) * p2
            - sp.diff(h, p1)
            + parameter_c * p1 / 4
            + parameter_c * p2 / 2
        ),
        sp.expand(q2 - g - parameter_c * p1 / 2 - parameter_c * p2),
        p1,
        p1 + p2,
    )
    jacobian = sp.Matrix(outputs).jacobian(DARBOUX_VARIABLES)
    assert jacobian == cayley
    assert sp.expand(jacobian.det()) == 1
    assert (jacobian * POISSON * jacobian.T).applyfunc(sp.expand) == POISSON

    # The inverse is triangular after reading p1 and p2 from the last two
    # outputs. Check both compositions with independent target symbols.
    Q1, Q2, P1, P2 = sp.symbols("Q1 Q2 P1 P2")
    g_target = g.subs(p1, P1)
    h_prime_target = sp.diff(h, p1).subs(p1, P1)
    g_prime_target = sp.diff(g, p1).subs(p1, P1)
    p1_inverse = P1
    p2_inverse = P2 - P1
    q2_inverse = sp.expand(
        Q2
        + g_target
        + parameter_c * P1 / 2
        + parameter_c * p2_inverse
    )
    q1_inverse = sp.expand(
        Q1
        + q2_inverse
        + g_prime_target * p2_inverse
        + h_prime_target
        - parameter_c * P1 / 4
        - parameter_c * p2_inverse / 2
    )
    inverse = (q1_inverse, q2_inverse, p1_inverse, p2_inverse)
    forward_substitution = dict(zip((Q1, Q2, P1, P2), outputs, strict=True))
    for recovered, expected in zip(inverse, DARBOUX_VARIABLES, strict=True):
        assert sp.expand(
            recovered.subs(forward_substitution, simultaneous=True) - expected
        ) == 0

    # Every higher odd Moyal term vanishes. The displayed coefficient-generic
    # degrees cover the only possible powers for this family.
    for left, right in itertools.combinations(range(4), 2):
        assert pi_power(outputs[left], outputs[right], 3) == 0
        assert pi_power(outputs[left], outputs[right], 5) == 0

    # Rigidity of the surrounding triangular ansatz. For
    # A=p1*q2+p2^2+f(p1,p2), these are the only two nonzero row curls.
    arbitrary_f = sp.Function("f")(p1, p2)
    general_potential = p1 * q2 + p2**2 + arbitrary_f
    general_nilpotent = POISSON * sp.hessian(
        general_potential, DARBOUX_VARIABLES
    )
    general_cayley = (
        sp.eye(4)
        + general_nilpotent
        + general_nilpotent**2 / 2
        + general_nilpotent**3 / 4
    )
    nonzero_curls = []
    for row in range(4):
        for left in range(4):
            for right in range(left + 1, 4):
                value = sp.simplify(
                    sp.diff(general_cayley[row, left], DARBOUX_VARIABLES[right])
                    - sp.diff(
                        general_cayley[row, right], DARBOUX_VARIABLES[left]
                    )
                )
                if value != 0:
                    nonzero_curls.append((row, left, right, value))
    f_222 = sp.diff(arbitrary_f, p2, 3)
    f_122 = sp.diff(arbitrary_f, p1, p2, 2)
    assert nonzero_curls == [
        (0, 2, 3, f_222 / 4 - f_122 / 2),
        (1, 2, 3, -f_222 / 2),
    ]

    return {
        "potential": "p1*q2+(c/2)*p2^2+g(p1)*p2+h(p1)",
        "coefficient_base": "Q[c,g_0,...,g_4,h_0,...,h_5]",
        "nilpotence": {
            "N_squared": str(expected_square.tolist()),
            "N_cubed": str(expected_cube.tolist()),
            "N_fourth_power_zero": True,
            "regular_index_four_open": "c != 0",
        },
        "cayley_map": [str(component) for component in outputs],
        "polynomial_inverse": [str(component) for component in inverse],
        "status": (
            "exact symplectic polynomial automorphism; all higher odd "
            "Moyal brackets vanish"
        ),
        "triangular_rigidity": {
            "ansatz": "p1*q2+p2^2+f(p1,p2)",
            "nonzero_row_curls": ["f_222/4-f_122/2", "-f_222/2"],
            "conclusion": (
                "row closure forces f_222=f_122=0, hence f_22 is "
                "constant and N^2 has no nonconstant part"
            ),
        },
    }


def verify_r21_frontier() -> dict[str, object]:
    """Verify exact R21 data and the natural stable denominator charts."""

    x, y, z, e, q = sp.symbols("x y z e q")
    variables = (x, y, z, e)
    modulus = q**2 - 4 * q + 6

    def red(expression: sp.Expr) -> sp.Expr:
        return reduce_q(sp.expand(expression), q, modulus)

    A = 1 + x * y**2
    h = q + (4 * q - 6) * A
    B = sp.expand(A**2 * z + y**3 * h)
    P = sp.expand(A * B)
    Q = sp.expand(y + x * B)
    s = x / A
    R = red(sp.cancel(fiber_antiderivative(2, 1, s, P, Q)))
    assert sp.denom(R) == 1
    r_cofactor = red(sp.cancel(R / x))
    assert sp.denom(r_cofactor) == 1
    assert red(R - x * r_cofactor) == 0
    assert red(r_cofactor.subs(x, 0) - 1) == 0
    assert red(r_cofactor - 1) != 0
    base_jacobian = sp.Matrix((P, Q, R)).jacobian((x, y, z))
    assert red(base_jacobian.det()) == -1

    # An affine contact symplectization would require a nonzero constant
    # target two-form whose pullback is again constant.  Write its three
    # coefficients over K as a_i+q*b_i and kill every positive-degree source
    # coefficient.  The resulting rational linear system has full rank six.
    rational_parts = sp.symbols("a0:3")
    q_parts = sp.symbols("b0:3")
    two_form_coefficients = tuple(
        rational_parts[index] + q * q_parts[index] for index in range(3)
    )
    target_two_form = sp.Matrix(
        (
            (0, two_form_coefficients[0], two_form_coefficients[1]),
            (-two_form_coefficients[0], 0, two_form_coefficients[2]),
            (-two_form_coefficients[1], -two_form_coefficients[2], 0),
        )
    )
    pulled_two_form = (
        base_jacobian.T * target_two_form * base_jacobian
    ).applyfunc(red)
    contact_equations = []
    for row in range(3):
        for column in range(row + 1, 3):
            for monomial, coefficient in sp.Poly(
                pulled_two_form[row, column], x, y, z, q
            ).terms():
                if sum(monomial[:3]) > 0:
                    contact_equations.append(coefficient)
    contact_matrix, _ = sp.linear_eq_to_matrix(
        contact_equations, rational_parts + q_parts
    )
    assert contact_matrix.shape == (262, 6)
    assert contact_matrix.rank() == 6

    T = sp.symbols("T")
    incidence = (
        T
        - Q**2 * T**2 / 2
        + 2 * P * Q * T**3 / 3
        - P**2 * T**4 / 4
        - R
    )
    assert red(sp.cancel(incidence.subs(T, s))) == 0
    assert sp.factor(
        sp.diff(incidence, T) - (1 - T * (Q - P * T) ** 2)
    ) == 0

    # The incidence root generates the complete source function field.  This
    # also makes the generic degree exactly four: over K(P,Q), the polynomial
    # R-g(P,Q,T) is primitive and linear in the independent variable R, hence
    # irreducible, while its T-degree is four.
    # Check the reconstruction in stages to avoid a useless expansion of the
    # final nested rational expression.
    assert red(sp.cancel(Q - P * s - y)) == 0
    assert sp.cancel(1 - s * y**2 - 1 / A) == 0
    assert sp.cancel(s / (1 - s * y**2) - x) == 0
    assert red(sp.cancel(P / A - B)) == 0
    assert red(sp.cancel((B - y**3 * h) / A**2 - z)) == 0

    # The graph form is dP^dQ+de^dR. Its coefficient matrix is a polynomial
    # symplectic form with determinant one.
    target_poisson = sp.Matrix(
        (
            (0, -1, 0, 0),
            (1, 0, 0, 0),
            (0, 0, 0, -1),
            (0, 0, 1, 0),
        )
    )
    graph_map = sp.Matrix((P, Q, e, R))
    graph_jacobian = graph_map.jacobian(variables)
    graph_form = (
        graph_jacobian.T * (-target_poisson) * graph_jacobian
    ).applyfunc(red)
    assert red(graph_form.det()) == 1

    # Put the graph map in tangent-identity coordinates
    # (u0,u1,u2,u3)=(z,y,e,x).  The only linear bracket defect is u0 in the
    # (1,2) entry.  The elementary automorphism u1 -> u1-u0*u3 removes it,
    # and the exact corrected defect has no terms of degrees one or two.
    u = sp.symbols("u0:4")
    tangent_substitution = {x: u[3], y: u[1], z: u[0], e: u[2]}
    normalized_graph = sp.Matrix(
        [
            red(component.subs(tangent_substitution, simultaneous=True))
            for component in graph_map
        ]
    )
    origin = {variable: 0 for variable in u}
    assert normalized_graph.jacobian(u).subs(origin) == sp.eye(4)

    def homogeneous_part(expression: sp.Expr, degree: int) -> sp.Expr:
        terms = sp.Poly(red(expression), *u, q).terms()
        return red(
            sum(
                coefficient
                * sp.prod(
                    variable**power
                    for variable, power in zip((*u, q), monomial, strict=True)
                )
                for monomial, coefficient in terms
                if sum(monomial[:4]) == degree
            )
        )

    raw_jacobian = normalized_graph.jacobian(u)
    raw_defect = (
        raw_jacobian * target_poisson * raw_jacobian.T - target_poisson
    ).applyfunc(red)
    raw_linear_defects = {
        (row, column): homogeneous_part(raw_defect[row, column], 1)
        for row in range(4)
        for column in range(row + 1, 4)
        if homogeneous_part(raw_defect[row, column], 1) != 0
    }
    assert raw_linear_defects == {(1, 2): u[0]}

    tangent_shear = sp.Matrix(
        (u[0], u[1] - u[0] * u[3], u[2], u[3])
    )
    tangent_shear_inverse = sp.Matrix(
        (u[0], u[1] + u[0] * u[3], u[2], u[3])
    )
    assert sp.expand(tangent_shear.jacobian(u).det()) == 1
    assert sp.Matrix(
        [
            component.subs(
                dict(zip(u, tangent_shear_inverse, strict=True)),
                simultaneous=True,
            )
            for component in tangent_shear
        ]
    ).applyfunc(sp.expand) == sp.Matrix(u)
    corrected_graph = sp.Matrix(
        [
            red(
                component.subs(
                    dict(zip(u, tangent_shear, strict=True)),
                    simultaneous=True,
                )
            )
            for component in normalized_graph
        ]
    )
    corrected_jacobian = corrected_graph.jacobian(u)
    corrected_defect = (
        corrected_jacobian
        * target_poisson
        * corrected_jacobian.T
        - target_poisson
    ).applyfunc(red)
    for row in range(4):
        for column in range(row + 1, 4):
            assert homogeneous_part(corrected_defect[row, column], 1) == 0
            assert homogeneous_part(corrected_defect[row, column], 2) == 0
    cubic_defects = {
        (0, 1): -3 * u[1] ** 2 * u[3],
        (0, 2): (21 - 15 * q) * u[0] * u[1] ** 2,
        (0, 3): 3 * u[1] * u[3] ** 2,
        (1, 2): (5 * q - 6) * u[1] ** 3,
        (1, 3): sp.Integer(0),
        (2, 3): 3 * u[1] ** 2 * u[3],
    }
    for pair, expected in cubic_defects.items():
        actual = homogeneous_part(corrected_defect[pair], 3)
        assert red(actual - expected) == 0, (pair, actual, expected)

    # A layer of degree-four coordinate shears cannot remove the cubic
    # paired defect.  For V0 independent of u0 and V1 independent of u1,
    # the linearized change in the (0,1) bracket is
    # -partial_0(V0)-partial_1(V1)=0.
    quartic_monomials_away_from_u0 = [
        u[1] ** left * u[2] ** middle * u[3] ** (4 - left - middle)
        for left in range(5)
        for middle in range(5 - left)
    ]
    quartic_monomials_away_from_u1 = [
        u[0] ** left * u[2] ** middle * u[3] ** (4 - left - middle)
        for left in range(5)
        for middle in range(5 - left)
    ]
    shear_coefficients = sp.symbols("v0_0:15") + sp.symbols("v1_0:15")
    vector_0 = sum(
        coefficient * monomial
        for coefficient, monomial in zip(
            shear_coefficients[:15],
            quartic_monomials_away_from_u0,
            strict=True,
        )
    )
    vector_1 = sum(
        coefficient * monomial
        for coefficient, monomial in zip(
            shear_coefficients[15:],
            quartic_monomials_away_from_u1,
            strict=True,
        )
    )
    assert -sp.diff(vector_0, u[0]) - sp.diff(vector_1, u[1]) == 0

    # Coupled shears do succeed.  The sparse degree-four correction vector is
    #
    #   (-3*u0*u1^2*u3, 0, (6-5q)*u0*u1^3, 3*u1^2*u3^2/2).
    #
    # Its (u0,u3)-part is Hamiltonian with
    # -3*u1^2*u0*u3^2/2.  Decompose u0*u3^2 into three cubes of linear forms;
    # each summand then integrates to an exact polynomial shear.  The next
    # degree-five correction is handled in the same way using a five-term
    # Waring decomposition of u0^2*u3^3.  We check exact inverses factor by
    # factor and use degree-six jets only to keep the certificate compact.
    def truncate_in_u(expression: sp.Expr, bound: int) -> sp.Expr:
        return red(
            sum(
                coefficient
                * sp.prod(
                    variable**power
                    for variable, power in zip((*u, q), monomial, strict=True)
                )
                for monomial, coefficient in sp.Poly(
                    red(expression), *u, q
                ).terms()
                if sum(monomial[:4]) <= bound
            )
        )

    def compose_jet(
        polynomial_map: sp.MatrixBase,
        source_change: sp.MatrixBase,
        bound: int,
    ) -> sp.Matrix:
        substitution = dict(zip(u, source_change, strict=True))
        return sp.Matrix(
            [
                truncate_in_u(
                    component.subs(substitution, simultaneous=True), bound
                )
                for component in polynomial_map
            ]
        )

    def linear_hamiltonian_shear(
        linear_form: sp.Expr,
        coefficient: sp.Expr,
        invariant_prefactor: sp.Expr,
        power: int,
    ) -> sp.Matrix:
        hamiltonian = coefficient * invariant_prefactor * linear_form**power
        shear = sp.Matrix(
            (
                u[0] + sp.diff(hamiltonian, u[3]),
                u[1],
                u[2],
                u[3] - sp.diff(hamiltonian, u[0]),
            )
        )
        inverse_hamiltonian = -hamiltonian
        inverse = sp.Matrix(
            (
                u[0] + sp.diff(inverse_hamiltonian, u[3]),
                u[1],
                u[2],
                u[3] - sp.diff(inverse_hamiltonian, u[0]),
            )
        )
        assert sp.expand(shear.jacobian(u).det()) == 1
        for left, right in ((shear, inverse), (inverse, shear)):
            substitution = dict(zip(u, right, strict=True))
            assert sp.Matrix(
                [
                    sp.expand(
                        component.subs(substitution, simultaneous=True)
                    )
                    for component in left
                ]
            ) == sp.Matrix(u)
        return shear

    def elementary_u2_shear(increment: sp.Expr) -> sp.Matrix:
        shear = sp.Matrix((u[0], u[1], u[2] + increment, u[3]))
        inverse = sp.Matrix((u[0], u[1], u[2] - increment, u[3]))
        assert sp.expand(shear.jacobian(u).det()) == 1
        substitution = dict(zip(u, inverse, strict=True))
        assert sp.Matrix(
            [
                sp.expand(component.subs(substitution, simultaneous=True))
                for component in shear
            ]
        ) == sp.Matrix(u)
        return shear

    cubic_correction_factors = [
        linear_hamiltonian_shear(
            u[0] + u[3], -sp.Rational(1, 4), u[1] ** 2, 3
        ),
        linear_hamiltonian_shear(
            u[0] - u[3], -sp.Rational(1, 4), u[1] ** 2, 3
        ),
        linear_hamiltonian_shear(
            u[0], sp.Rational(1, 2), u[1] ** 2, 3
        ),
        elementary_u2_shear(
            (6 - 5 * q) * u[0] * u[1] ** 3
        ),
    ]

    jet_map = sp.Matrix(
        [truncate_in_u(component, 6) for component in normalized_graph]
    )
    jet_map = compose_jet(jet_map, tangent_shear, 6)
    for factor in cubic_correction_factors:
        jet_map = compose_jet(jet_map, factor, 6)
    cubic_corrected_jacobian = jet_map.jacobian(u)
    cubic_corrected_defect = (
        cubic_corrected_jacobian
        * target_poisson
        * cubic_corrected_jacobian.T
        - target_poisson
    ).applyfunc(lambda entry: truncate_in_u(entry, 4))
    for row in range(4):
        for column in range(row + 1, 4):
            for degree in range(1, 4):
                assert (
                    homogeneous_part(
                        cubic_corrected_defect[row, column], degree
                    )
                    == 0
                )
    quartic_defects = {
        (0, 1): 8 * u[0] * u[1] * u[3] ** 2,
        (0, 2): (30 * q - 48) * u[0] ** 2 * u[1] * u[3],
        (0, 3): -sp.Rational(8, 3) * u[0] * u[3] ** 3,
        (1, 2): (40 - 30 * q) * u[0] * u[1] ** 2 * u[3],
        (1, 3): sp.Rational(8, 3) * u[1] * u[3] ** 3,
        (2, 3): -8 * u[0] * u[1] * u[3] ** 2,
    }
    for pair, expected in quartic_defects.items():
        assert red(
            homogeneous_part(cubic_corrected_defect[pair], 4) - expected
        ) == 0

    quartic_correction_factors = [
        linear_hamiltonian_shear(
            u[0] + u[3], -sp.Rational(1, 45), u[1], 5
        ),
        linear_hamiltonian_shear(
            u[0] - u[3], sp.Rational(1, 45), u[1], 5
        ),
        linear_hamiltonian_shear(
            u[0] + 2 * u[3], sp.Rational(1, 90), u[1], 5
        ),
        linear_hamiltonian_shear(
            u[0] - 2 * u[3], -sp.Rational(1, 90), u[1], 5
        ),
        linear_hamiltonian_shear(
            u[3], -sp.Rational(2, 3), u[1], 5
        ),
        elementary_u2_shear(
            (15 * q - 20) * u[0] ** 2 * u[1] ** 2 * u[3]
        ),
    ]
    for factor in quartic_correction_factors:
        jet_map = compose_jet(jet_map, factor, 6)
    quartic_corrected_jacobian = jet_map.jacobian(u)
    quartic_corrected_defect = (
        quartic_corrected_jacobian
        * target_poisson
        * quartic_corrected_jacobian.T
        - target_poisson
    ).applyfunc(lambda entry: truncate_in_u(entry, 5))
    for row in range(4):
        for column in range(row + 1, 4):
            for degree in range(1, 5):
                assert (
                    homogeneous_part(
                        quartic_corrected_defect[row, column], degree
                    )
                    == 0
                )
    quintic_defects = {
        (0, 1): -5 * u[0] ** 2 * u[3] ** 3,
        (0, 2): (27 - 15 * q) * u[0] ** 3 * u[3] ** 2
        + (330 - 111 * q) * u[1] ** 5,
        (0, 3): sp.Integer(0),
        (1, 2): (45 * q - 66) * u[0] ** 2 * u[1] * u[3] ** 2,
        (1, 3): -sp.Rational(5, 2) * u[0] * u[3] ** 4,
        (2, 3): 5 * u[0] ** 2 * u[3] ** 3,
    }
    for pair, expected in quintic_defects.items():
        assert red(
            homogeneous_part(quartic_corrected_defect[pair], 5) - expected
        ) == 0

    # Every isolated power A^n has an elementary one-variable stabilization:
    # A^n-x*y^2*(1+A+...+A^(n-1))=1.
    stable_charts = []
    for exponent in range(1, 5):
        quotient = y**2 * sum(A**degree for degree in range(exponent))
        matrix = sp.Matrix(((A**exponent, x), (quotient, 1)))
        inverse = sp.Matrix(((1, -x), (-quotient, A**exponent)))
        assert sp.expand(matrix.det()) == 1
        assert (matrix * inverse).applyfunc(sp.expand) == sp.eye(2)
        stable_charts.append(
            {
                "power": exponent,
                "matrix": str(matrix.tolist()),
                "inverse": str(inverse.tolist()),
            }
        )

    # The two natural affine stable coordinates absorb B and P. Verify that
    # they are coordinates, but also that neither inverse chart alone is a
    # Darboux trivialization of the graph form.
    natural_chart_defects = []
    for exponent, translation, label in (
        (2, y**3 * h, "B+x*e"),
        (3, A * y**3 * h, "P+x*e"),
    ):
        quotient = y**2 * sum(A**degree for degree in range(exponent))
        new_z = sp.expand(A**exponent * z + x * e + translation)
        new_e = sp.expand(quotient * z + e)
        old_z = sp.expand((new_z - translation) - x * new_e)
        old_e = sp.expand(
            -quotient * (new_z - translation) + A**exponent * new_e
        )
        assert sp.expand(old_z - z) == 0
        assert sp.expand(old_e - e) == 0

        X, Y, Z, E = sp.symbols("X Y Z E")
        A_new = 1 + X * Y**2
        h_new = q + (4 * q - 6) * A_new
        translation_new = (
            Y**3 * h_new if exponent == 2 else A_new * Y**3 * h_new
        )
        quotient_new = Y**2 * sum(
            A_new**degree for degree in range(exponent)
        )
        inverse_substitution = {
            x: X,
            y: Y,
            z: sp.expand((Z - translation_new) - X * E),
            e: sp.expand(
                -quotient_new * (Z - translation_new)
                + A_new**exponent * E
            ),
        }
        pulled_outputs = sp.Matrix(
            [
                red(component.subs(inverse_substitution, simultaneous=True))
                for component in (P, Q, e, R)
            ]
        )
        pulled_jacobian = pulled_outputs.jacobian((X, Y, Z, E))
        pulled_poisson = (
            pulled_jacobian * target_poisson * pulled_jacobian.T
        ).applyfunc(red)
        defect_count = sum(
            red(pulled_poisson[row, column] - target_poisson[row, column])
            != 0
            for row in range(4)
            for column in range(row + 1, 4)
        )
        assert defect_count == 6
        natural_chart_defects.append(
            {
                "stable_coordinate": label,
                "nonzero_symplectic_pair_defects": defect_count,
            }
        )

    # The fiber-preserving no-go makes stable mixing mandatory.  On the U_2
    # chart, tangent-normalize by (u0,u1,u2,u3)=(Z,Y,E,X).  Two elementary
    # stable-mixing shears kill all linear defects, one u2-shear kills the
    # quadratic defect, and the three cubic Hamiltonian shears already used
    # above kill the cubic defect.  The first remaining term is quartic.
    X, Y, Z, E = sp.symbols("X Y Z E")
    A_mixed = 1 + X * Y**2
    h_mixed = q + (4 * q - 6) * A_mixed
    quotient_mixed = Y**2 * (1 + A_mixed)
    translation_mixed = Y**3 * h_mixed
    mixed_inverse_substitution = {
        x: X,
        y: Y,
        z: sp.expand(Z - translation_mixed - X * E),
        e: sp.expand(
            -quotient_mixed * (Z - translation_mixed)
            + A_mixed**2 * E
        ),
    }
    mixed_outputs = sp.Matrix(
        [
            red(
                component.subs(
                    mixed_inverse_substitution, simultaneous=True
                ).subs(
                    {X: u[3], Y: u[1], Z: u[0], E: u[2]},
                    simultaneous=True,
                )
            )
            for component in (P, Q, e, R)
        ]
    )
    assert mixed_outputs.jacobian(u).subs(origin) == sp.eye(4)
    mixed_raw_jacobian = mixed_outputs.jacobian(u)
    mixed_raw_defect = (
        mixed_raw_jacobian
        * target_poisson
        * mixed_raw_jacobian.T
        - target_poisson
    ).applyfunc(red)

    # The straight-line Moser path from the standard form W=-Pi to the exact
    # U_2 pullback form Omega has
    #
    #   pf((1-s)W+s*Omega)=1+s(1-s)*Delta.
    #
    # Compute Delta without introducing s: it is the polarization of the
    # four-dimensional Pfaffian at (W, Omega-W).  A nonzero Delta proves that
    # the canonical scalar geometric series does not terminate, even though
    # it may still admit a special closed-form resummation.
    mixed_standard_form = -target_poisson

    def pfaffian_four(matrix: sp.MatrixBase) -> sp.Expr:
        return red(
            matrix[0, 1] * matrix[2, 3]
            - matrix[0, 2] * matrix[1, 3]
            + matrix[0, 3] * matrix[1, 2]
        )

    mixed_pullback_form = (
        mixed_raw_jacobian.T
        * mixed_standard_form
        * mixed_raw_jacobian
    ).applyfunc(red)
    assert pfaffian_four(mixed_standard_form) == 1
    assert pfaffian_four(mixed_pullback_form) == 1
    mixed_form_difference = (
        mixed_pullback_form - mixed_standard_form
    ).applyfunc(red)
    moser_delta = red(
        mixed_standard_form[0, 1] * mixed_form_difference[2, 3]
        + mixed_form_difference[0, 1] * mixed_standard_form[2, 3]
        - mixed_standard_form[0, 2] * mixed_form_difference[1, 3]
        - mixed_form_difference[0, 2] * mixed_standard_form[1, 3]
        + mixed_standard_form[0, 3] * mixed_form_difference[1, 2]
        + mixed_form_difference[0, 3] * mixed_standard_form[1, 2]
    )
    moser_quadratic_coefficient = pfaffian_four(mixed_form_difference)
    assert red(moser_quadratic_coefficient + moser_delta) == 0
    assert moser_delta != 0
    moser_delta_terms = sp.Poly(moser_delta, *u, q).terms()
    moser_delta_min_degree = min(
        sum(monomial[:4]) for monomial, _ in moser_delta_terms
    )
    moser_delta_max_degree = max(
        sum(monomial[:4]) for monomial, _ in moser_delta_terms
    )
    moser_delta_leading_part = homogeneous_part(
        moser_delta, moser_delta_min_degree
    )

    # A reciprocal-compatible constant-Pfaffian path is supplied by source
    # dilation.  Since the U_2 graph map K is tangent to the identity, put
    # K_s(u)=K(su)/s.  Its Jacobian is J_K(su), hence still has determinant
    # one.  The canonical trivializing field
    #
    #   V_s=-(J K_s)^(-1) partial_s K_s
    #
    # is therefore polynomial: the inverse Jacobian is its adjugate.  Its
    # formal flow satisfies K_s o H_s=id.  At s=1 this is the formal inverse
    # of K, so it cannot be a polynomial time-one map because the incidence
    # equation makes K generically degree four.
    path_parameter = sp.symbols("s")
    dilation_substitution = {
        variable: path_parameter * variable for variable in u
    }
    dilation_map = sp.Matrix(
        [
            red(
                component.subs(
                    dilation_substitution, simultaneous=True
                )
                / path_parameter
            )
            for component in mixed_outputs
        ]
    )
    assert dilation_map.subs(path_parameter, 0) == sp.Matrix(u)
    dilation_jacobian = dilation_map.jacobian(u)
    assert red(dilation_jacobian.det()) == 1
    dilation_parameter_derivative = dilation_map.diff(path_parameter)
    dilation_field = (
        -dilation_jacobian.adjugate() * dilation_parameter_derivative
    ).applyfunc(red)
    assert (
        dilation_jacobian * dilation_field
        + dilation_parameter_derivative
    ).applyfunc(red) == sp.zeros(4, 1)
    dilation_field_terms = [
        sp.Poly(component, *u, path_parameter, q).terms()
        for component in dilation_field
    ]
    dilation_field_term_counts = [
        len(terms) for terms in dilation_field_terms
    ]
    dilation_field_source_degrees = [
        max(sum(monomial[:4]) for monomial, _ in terms)
        for terms in dilation_field_terms
    ]
    dilation_field_parameter_degrees = [
        max(monomial[4] for monomial, _ in terms)
        for terms in dilation_field_terms
    ]
    assert dilation_field_term_counts == [301, 142, 303, 149]
    assert dilation_field_source_degrees == [54, 40, 53, 43]
    assert dilation_field_parameter_degrees == [52, 38, 51, 41]

    # All twenty-six explicit corrections preserve the b=u1 fibration.
    # Exclude the entire polynomial automorphism subgroup preserving K[b],
    # rather than just those factors.  Such an automorphism sends b to
    # lambda*b+mu.  On a b-fiber, Omega has rank two in
    # (a,c,d)=(u0,u2,u3).
    # Relative to da^dc^dd its primitive kernel derivation is
    #
    #   delta=(Omega_cd, -Omega_ad, Omega_ac).
    #
    # If a b-fibration-preserving polynomial Darboux normalizer existed, its
    # restriction from the fiber b=-mu/lambda to b=0 would conjugate this
    # derivation, up to a nonzero constant, to partial_a and hence make it
    # locally nilpotent.  Specialization at b=0 contradicts that: d divides
    # delta(d)=d^4(a-cd)/6 but delta(d) is nonzero.  For an LND on a domain,
    # f|delta(f) implies delta(f)=0.
    fiber_kernel = sp.Matrix(
        (
            mixed_pullback_form[2, 3],
            -mixed_pullback_form[0, 3],
            mixed_pullback_form[0, 2],
        )
    ).applyfunc(red)
    fiber_kernel_at_zero = fiber_kernel.subs(u[1], 0).applyfunc(red)
    expected_fiber_kernel_at_zero = sp.Matrix(
        (
            1
            - u[0] * u[3]
            + u[2] * u[3] ** 2
            - u[0] ** 2 * u[3] ** 3 / 3
            + sp.Rational(5, 6) * u[0] * u[2] * u[3] ** 4
            - u[2] ** 2 * u[3] ** 5 / 2,
            -u[0] + u[2] * u[3],
            u[3] ** 4 * (u[0] - u[2] * u[3]) / 6,
        )
    )
    assert (
        fiber_kernel_at_zero - expected_fiber_kernel_at_zero
    ).applyfunc(red) == sp.zeros(3, 1)
    assert fiber_kernel_at_zero[0].subs(
        {u[0]: 0, u[2]: 0, u[3]: 0}
    ) == 1
    assert fiber_kernel_at_zero[2] != 0
    assert sp.rem(
        fiber_kernel_at_zero[2], u[3], u[0], u[2], u[3]
    ) == 0

    # Strengthen non-local-nilpotence to absence of any polynomial slice.
    # Put t=a-cd.  The kernel derivation has two invariant/semi-invariant
    # equations I=d^4*t^2-12d and J=d^3*c-2.  On I=1,J=0, the coordinate
    # y=d^2*t identifies the curve with
    # K[y,1/(y^2-1)] and delta(y)=(y^2-1)^2/144.  A slice would give a
    # rational primitive of 144/(y^2-1)^2, impossible because this rational
    # differential has nonzero residues at y=1 and y=-1.
    fiber_variables = (u[0], u[2], u[3])

    def fiber_delta(expression: sp.Expr) -> sp.Expr:
        return red(
            sum(
                fiber_kernel_at_zero[index]
                * sp.diff(expression, fiber_variables[index])
                for index in range(3)
            )
        )

    fiber_t = u[0] - u[2] * u[3]
    fiber_invariant = u[3] ** 4 * fiber_t**2 - 12 * u[3]
    fiber_semi_invariant = u[3] ** 3 * u[2] - 2
    fiber_curve_coordinate = u[3] ** 2 * fiber_t
    assert fiber_delta(fiber_invariant) == 0
    assert red(
        fiber_delta(fiber_semi_invariant)
        - u[3] ** 3 * fiber_t * fiber_semi_invariant / 2
    ) == 0
    assert red(fiber_delta(fiber_curve_coordinate) - u[3] ** 2) == 0
    residue_variable = sp.symbols("y_curve")
    curve_time_form = 144 / (residue_variable**2 - 1) ** 2
    curve_time_residues = [
        sp.residue(curve_time_form, residue_variable, point)
        for point in (1, -1)
    ]
    assert curve_time_residues == [-36, 36]

    # The no-slice condition is attainable by tame coordinates, so it is a
    # genuine selection rule rather than a disguised impossibility.  For
    # k>=1, compose the two elementary maps
    #
    #   G=b+a*d,        F_k=a+c^k*G.
    #
    # On F_k=0, with (G,c,d) as coordinates, the Hamiltonian derivation is
    #
    #   D_k(G)=k*G^2*c^(2k-1)-1,
    #   D_k(c)=-G*c^(2k),
    #   D_k(d)=-k*G*c^(k-1).
    #
    # Put t=G*c^k, I=t^2-2c and J=1+c^k*d.  Then D_k(t)=-c^k,
    # D_k(I)=0 and D_k(J)=-k*t*c^(k-1)*J.  On I=1,J=0 the time form is
    # -2^k*dt/(t^2-1)^k and has nonzero residues.  The note proves these
    # formulas for arbitrary k; verify the first three exact rows here.
    tame_no_slice_rows = []
    calibration_coordinates = sp.symbols("F_cal G_cal c_cal d_cal")
    calibration_f, calibration_g, calibration_c, calibration_d = (
        calibration_coordinates
    )
    calibration_curve_variable = sp.symbols("t_cal")
    for exponent in (1, 2, 3):
        tame_g = u[1] + u[0] * u[3]
        tame_f = u[0] + u[2] ** exponent * tame_g
        tame_map = sp.Matrix((tame_f, tame_g, u[2], u[3]))
        tame_jacobian = tame_map.jacobian(u)
        assert sp.expand(tame_jacobian.det()) == 1
        tame_brackets = (tame_jacobian * target_poisson * tame_jacobian.T)
        tame_inverse = {
            u[0]: calibration_f
            - calibration_c**exponent * calibration_g,
            u[1]: calibration_g
            - (
                calibration_f
                - calibration_c**exponent * calibration_g
            )
            * calibration_d,
            u[2]: calibration_c,
            u[3]: calibration_d,
        }
        tame_fiber_kernel = [
            sp.expand(
                tame_brackets[0, column]
                .subs(tame_inverse, simultaneous=True)
                .subs(calibration_f, 0)
            )
            for column in range(1, 4)
        ]
        expected_tame_fiber_kernel = [
            exponent
            * calibration_g**2
            * calibration_c ** (2 * exponent - 1)
            - 1,
            -calibration_g * calibration_c ** (2 * exponent),
            -exponent
            * calibration_g
            * calibration_c ** (exponent - 1),
        ]
        assert tame_fiber_kernel == expected_tame_fiber_kernel

        def tame_fiber_delta(expression: sp.Expr) -> sp.Expr:
            return sp.expand(
                sum(
                    tame_fiber_kernel[index]
                    * sp.diff(
                        expression,
                        (
                            calibration_g,
                            calibration_c,
                            calibration_d,
                        )[index],
                    )
                    for index in range(3)
                )
            )

        tame_t = calibration_g * calibration_c**exponent
        tame_invariant = tame_t**2 - 2 * calibration_c
        tame_semi_invariant = (
            1 + calibration_c**exponent * calibration_d
        )
        assert tame_fiber_delta(tame_t) == -calibration_c**exponent
        assert tame_fiber_delta(tame_invariant) == 0
        assert sp.expand(
            tame_fiber_delta(tame_semi_invariant)
            + exponent
            * tame_t
            * calibration_c ** (exponent - 1)
            * tame_semi_invariant
        ) == 0
        tame_time_form = (
            -(2**exponent)
            / (calibration_curve_variable**2 - 1) ** exponent
        )
        tame_residues = [
            sp.residue(tame_time_form, calibration_curve_variable, point)
            for point in (1, -1)
        ]
        assert all(value != 0 for value in tame_residues)
        tame_m = (
            calibration_g**2
            * calibration_c ** (2 * exponent - 1)
            - 2
        )
        tame_constant_i = calibration_c * tame_m
        tame_constant_j = (
            1 + calibration_c**exponent * calibration_d
        )
        tame_constant_s = sp.expand(
            tame_m**exponent * tame_constant_j
        )
        tame_constant_sigma = (-2) ** exponent
        tame_constant_u = sp.cancel(
            tame_constant_s
            * (tame_constant_s - tame_constant_sigma)
            / tame_constant_i**exponent
        )
        assert sp.denom(tame_constant_u) == 1
        tame_constant_u = sp.expand(tame_constant_u)
        for invariant in (
            tame_constant_i,
            tame_constant_s,
            tame_constant_u,
        ):
            assert tame_fiber_delta(invariant) == 0
        assert sp.expand(
            tame_constant_i**exponent * tame_constant_u
            - tame_constant_s
            * (tame_constant_s - tame_constant_sigma)
        ) == 0
        tame_no_slice_rows.append(
            {
                "exponent": exponent,
                "coordinate_degree": exponent + 2,
                "time_form_pole_order": exponent,
                "transverse_weight": exponent,
                "constant_ring_danielewski_exponent": exponent,
                "constant_ring_special_value": tame_constant_sigma,
                "time_form_residues_at_1_and_minus_1": [
                    str(value) for value in tame_residues
                ],
            }
        )

    # Bounded Darboux-polynomial calibration for the best tangential match
    # k=2.  For a fixed logarithmic cofactor q, solve D(p)=weight*q*p in the
    # space of polynomials of total degree at most seven.  This is only a
    # finite exact census, not a classification theorem.  It shows that the
    # first extra semi-invariant beside powers of the boundary coordinate
    # occurs at weight two for the tame model, but at weight three for R21.
    def bounded_weight_space_dimension(
        variables: tuple[sp.Symbol, ...],
        components: tuple[sp.Expr, ...],
        logarithmic_cofactor: sp.Expr,
        weight: int,
        degree_bound: int,
    ) -> int:
        monomials = sorted(
            sp.itermonomials(variables, degree_bound),
            key=lambda monomial: (
                sp.total_degree(monomial),
                sp.default_sort_key(monomial),
            ),
        )
        coefficients = sp.symbols(f"weight_{weight}_0:{len(monomials)}")
        polynomial = sum(
            coefficient * monomial
            for coefficient, monomial in zip(
                coefficients, monomials, strict=True
            )
        )
        equation = sp.Poly(
            sp.expand(
                sum(
                    components[index]
                    * sp.diff(polynomial, variables[index])
                    for index in range(3)
                )
                - weight * logarithmic_cofactor * polynomial
            ),
            *variables,
        )
        matrix, _ = sp.linear_eq_to_matrix(
            equation.coeffs(), coefficients
        )
        return len(coefficients) - matrix.rank()

    eigenspace_degree_bound = 7
    eigenspace_weights = (1, 2, 3, 4)
    tame_two_kernel = (
        2 * calibration_g**2 * calibration_c**3 - 1,
        -calibration_g * calibration_c**4,
        -2 * calibration_g * calibration_c,
    )

    def tame_two_delta(expression: sp.Expr) -> sp.Expr:
        return sp.expand(
            sum(
                tame_two_kernel[index]
                * sp.diff(
                    expression,
                    (calibration_g, calibration_c, calibration_d)[index],
                )
                for index in range(3)
            )
        )

    tame_two_logarithmic_cofactor = -calibration_g * calibration_c**3
    tame_two_weight_dimensions = [
        bounded_weight_space_dimension(
            (calibration_g, calibration_c, calibration_d),
            tame_two_kernel,
            tame_two_logarithmic_cofactor,
            weight,
            eigenspace_degree_bound,
        )
        for weight in eigenspace_weights
    ]
    target_logarithmic_cofactor = (
        u[3] ** 3 * (u[0] - u[2] * u[3]) / 6
    )
    target_weight_dimensions = [
        bounded_weight_space_dimension(
            fiber_variables,
            tuple(fiber_kernel_at_zero),
            target_logarithmic_cofactor,
            weight,
            eigenspace_degree_bound,
        )
        for weight in eigenspace_weights
    ]
    assert tame_two_weight_dimensions == [2, 2, 2, 3]
    assert target_weight_dimensions == [1, 1, 2, 2]
    tame_two_extra_semi_invariant = 1 + calibration_c**2 * calibration_d
    assert sp.expand(
        sum(
            tame_two_kernel[index]
            * sp.diff(
                tame_two_extra_semi_invariant,
                (calibration_g, calibration_c, calibration_d)[index],
            )
            for index in range(3)
        )
        - 2
        * tame_two_logarithmic_cofactor
        * tame_two_extra_semi_invariant
    ) == 0
    assert red(
        fiber_delta(fiber_semi_invariant)
        - 3 * target_logarithmic_cofactor * fiber_semi_invariant
    ) == 0

    # The degree-seven weight census is not an all-degree prime
    # classification.  At degree eight the k=2 tame derivation already has
    # another irreducible weight-one Darboux polynomial.  The correct global
    # invariant is its complete polynomial constant ring.
    tame_two_m = calibration_g**2 * calibration_c**3 - 2
    tame_two_i = calibration_c * tame_two_m
    tame_two_j = 1 + calibration_c**2 * calibration_d
    tame_two_s = sp.expand(tame_two_m**2 * tame_two_j)
    tame_two_sigma = 4
    tame_two_u = sp.cancel(
        tame_two_s * (tame_two_s - tame_two_sigma) / tame_two_i**2
    )
    assert sp.denom(tame_two_u) == 1
    tame_two_u = sp.expand(tame_two_u)
    for invariant in (tame_two_i, tame_two_s, tame_two_u):
        assert tame_two_delta(invariant) == 0
    assert sp.expand(
        tame_two_i**2 * tame_two_u
        - tame_two_s * (tame_two_s - tame_two_sigma)
    ) == 0

    tame_two_degree_eight_prime = sp.expand(
        tame_two_m * tame_two_j + calibration_c
    )
    assert sp.Poly(
        tame_two_degree_eight_prime,
        calibration_g,
        calibration_c,
        calibration_d,
    ).total_degree() == 8
    assert sp.gcd(
        sp.Poly(
            tame_two_degree_eight_prime,
            calibration_d,
        ).coeff_monomial(calibration_d),
        sp.Poly(
            tame_two_degree_eight_prime,
            calibration_d,
        ).coeff_monomial(1),
    ) == 1
    assert sp.expand(
        tame_two_delta(tame_two_degree_eight_prime)
        - tame_two_logarithmic_cofactor
        * tame_two_degree_eight_prime
    ) == 0

    # R21 has the analogous constant ring with exponent three.  Put
    # N=d^3*(a-c*d)^2-12, I=d*N, J=d^3*c-2 and S=N^3*J.  The quotient
    # S*(S-3456)/I^3 is again polynomial and invariant.  The note proves by
    # localization and a two-branch valuation/intersection argument that
    # these displayed generators are the complete kernels:
    #
    #   ker(D_2)=K[I,S,U]/(I^2*U-S*(S-4)),
    #   ker(delta_0)=K[I,S,U]/(I^3*U-S*(S-3456)).
    #
    # They are Danielewski surfaces with different intrinsic exponents.
    target_n = u[3] ** 3 * fiber_t**2 - 12
    target_i = sp.expand(u[3] * target_n)
    target_j = fiber_semi_invariant
    target_s = sp.expand(target_n**3 * target_j)
    target_sigma = 3456
    target_u = sp.cancel(
        target_s * (target_s - target_sigma) / target_i**3
    )
    assert sp.denom(target_u) == 1
    target_u = sp.expand(target_u)
    assert red(target_i - fiber_invariant) == 0
    for invariant in (target_i, target_s, target_u):
        assert fiber_delta(invariant) == 0
    assert sp.expand(
        target_i**3 * target_u
        - target_s * (target_s - target_sigma)
    ) == 0

    # The k=2 tame model is nevertheless exactly conjugate to R21 after
    # deleting I=0.  This is the closest positive construction found so far.
    # In target variables set G=t/36, C=6d and
    # W=-1/(36*d^2)+I*(J/d^3)/864.  W is polynomial, and the map intertwines
    # delta_0 with -D_2/36.  Its Jacobian is -I/5184, so it becomes an
    # isomorphism on I!=0 but cannot extend to a polynomial automorphism.
    dense_tame_g = fiber_t / 36
    dense_tame_c = 6 * u[3]
    dense_tame_w = sp.cancel(
        -sp.Rational(1, 36) / u[3] ** 2
        + target_i * (target_j / u[3] ** 3) / 864
    )
    assert sp.denom(dense_tame_w) == 1
    dense_tame_w = sp.expand(dense_tame_w)
    dense_tame_coordinates = sp.Matrix(
        (dense_tame_g, dense_tame_c, dense_tame_w)
    )
    assert red(
        fiber_delta(dense_tame_g)
        + (2 * dense_tame_g**2 * dense_tame_c**3 - 1) / 36
    ) == 0
    assert red(
        fiber_delta(dense_tame_c)
        - dense_tame_g * dense_tame_c**4 / 36
    ) == 0
    assert red(
        fiber_delta(dense_tame_w)
        - dense_tame_g * dense_tame_c / 18
    ) == 0
    dense_conjugacy_jacobian = red(
        dense_tame_coordinates.jacobian(fiber_variables).det()
    )
    assert red(dense_conjugacy_jacobian + target_i / 5184) == 0
    assert red(
        dense_tame_g**2 * dense_tame_c**4
        - 2 * dense_tame_c
        - target_i
    ) == 0
    assert red(
        (1 + dense_tame_c**2 * dense_tame_w) / dense_tame_c**2
        - target_i * (target_j / u[3] ** 3) / 864
    ) == 0

    # The split signature is in fact attained by a new degree-seven
    # coordinate.  Work first in a standard symplectic four-space
    # (G,P,C,Q), put Q'=Q+G^2/2, and use the Bezout pair
    #
    #   A=G*C^4/6, B=G^2*C^3/3-1,
    #   R=2*G^3*C^2/3, S=-(B+2),
    #   R*A+S*B=1.
    #
    # Then F=A*Q'-B*P and W=R*P+S*Q' complete to a determinant-one
    # polynomial coordinate system (F,G,C,W).  On F=0 the affine map
    # (G,C,W)->(a,c,d)=(G+C*W,W,C) pulls the complete R21 fiber form back
    # to the restricted standard form, not merely its characteristic line.
    symbol_g, symbol_p, symbol_c, symbol_q = sp.symbols("G P C Q")
    symbol_variables = (symbol_g, symbol_p, symbol_c, symbol_q)
    symbol_a = symbol_g * symbol_c**4 / 6
    symbol_b = symbol_g**2 * symbol_c**3 / 3 - 1
    symbol_r = 2 * symbol_g**3 * symbol_c**2 / 3
    symbol_s = -(symbol_b + 2)
    shifted_q = symbol_q + symbol_g**2 / 2
    assert sp.expand(symbol_r * symbol_a + symbol_s * symbol_b) == 1
    admission_f = sp.expand(symbol_a * shifted_q - symbol_b * symbol_p)
    admission_w = sp.expand(symbol_r * symbol_p + symbol_s * shifted_q)
    admission_coordinates = sp.Matrix(
        (admission_f, symbol_g, symbol_c, admission_w)
    )
    assert sp.expand(
        admission_coordinates.jacobian(symbol_variables).det()
    ) == 1
    assert sp.Poly(admission_f, *symbol_variables).total_degree() == 7

    capital_f, capital_w = sp.symbols("F W")
    inverse_p = sp.expand(
        symbol_a * capital_w + (symbol_b + 2) * capital_f
    )
    inverse_q = sp.expand(
        symbol_b * capital_w + symbol_r * capital_f - symbol_g**2 / 2
    )
    inverse_substitution = {
        symbol_p: inverse_p,
        symbol_q: inverse_q,
    }
    assert sp.expand(
        admission_f.subs(inverse_substitution, simultaneous=True)
        - capital_f
    ) == 0
    assert sp.expand(
        admission_w.subs(inverse_substitution, simultaneous=True)
        - capital_w
    ) == 0

    # Classify the complete reciprocal-compatible freedom in this Bezout
    # ansatz.  Replacing Q'=Q+G^2/2 by Q+H(G,C) changes the fiber bracket to
    # {F,W}=-H_G.  Matching R21 forces H=G^2/2+k(C).  Such k(C) is merely the
    # ambient symplectic shear Q->Q+k(C), and the full source form in
    # (F,G,C,W) is unchanged.  Every other Bezout pair is
    # (R+B*T,S-A*T), which only replaces W by W-T*F.
    arbitrary_hamiltonian_shift = sp.Function("H")(symbol_g, symbol_c)
    arbitrary_shifted_q = symbol_q + arbitrary_hamiltonian_shift
    arbitrary_admission_f = sp.expand(
        symbol_a * arbitrary_shifted_q - symbol_b * symbol_p
    )
    arbitrary_admission_w = sp.expand(
        symbol_r * symbol_p + symbol_s * arbitrary_shifted_q
    )
    arbitrary_inverse_q = sp.expand(
        symbol_b * capital_w
        + symbol_r * capital_f
        - arbitrary_hamiltonian_shift
    )
    arbitrary_c_shift = sp.Function("k")(symbol_c)
    general_shifted_q = shifted_q + arbitrary_c_shift
    general_admission_f = sp.expand(
        symbol_a * general_shifted_q - symbol_b * symbol_p
    )
    general_admission_w = sp.expand(
        symbol_r * symbol_p + symbol_s * general_shifted_q
    )
    general_inverse_q = sp.expand(
        symbol_b * capital_w
        + symbol_r * capital_f
        - symbol_g**2 / 2
        - arbitrary_c_shift
    )
    general_inverse_source = sp.Matrix(
        (symbol_g, inverse_p, symbol_c, general_inverse_q)
    )
    general_source_form = (
        general_inverse_source.jacobian(
            (capital_f, symbol_g, symbol_c, capital_w)
        ).T
        * mixed_standard_form
        * general_inverse_source.jacobian(
            (capital_f, symbol_g, symbol_c, capital_w)
        )
    ).applyfunc(sp.expand)
    bezout_parameter = sp.Function("T")(symbol_g, symbol_c)
    general_bezout_r = symbol_r + symbol_b * bezout_parameter
    general_bezout_s = symbol_s - symbol_a * bezout_parameter
    general_bezout_w = sp.expand(
        general_bezout_r * symbol_p
        + general_bezout_s * shifted_q
    )
    assert sp.expand(
        general_bezout_r * symbol_a
        + general_bezout_s * symbol_b
    ) == 1
    assert sp.expand(
        general_bezout_w
        - admission_w
        + bezout_parameter * admission_f
    ) == 0

    def symbol_bracket(left: sp.Expr, right: sp.Expr) -> sp.Expr:
        return sp.expand(
            (
                sp.Matrix([left]).jacobian(symbol_variables)
                * target_poisson
                * sp.Matrix([right]).jacobian(symbol_variables).T
            )[0]
        )

    arbitrary_fiber_bracket = sp.expand(
        symbol_bracket(arbitrary_admission_f, arbitrary_admission_w)
        .subs(
            {
                symbol_p: inverse_p,
                symbol_q: arbitrary_inverse_q,
            },
            simultaneous=True,
        )
        .subs(capital_f, 0)
    )
    assert arbitrary_fiber_bracket == -sp.diff(
        arbitrary_hamiltonian_shift, symbol_g
    )

    admission_brackets = (
        symbol_bracket(admission_f, symbol_g),
        symbol_bracket(admission_f, symbol_c),
        sp.expand(
            symbol_bracket(admission_f, admission_w)
            .subs(inverse_substitution, simultaneous=True)
            .subs(capital_f, 0)
        ),
    )
    assert admission_brackets == (
        -symbol_b,
        symbol_a,
        -symbol_g,
    )

    new_variables = (capital_f, symbol_g, symbol_c, capital_w)
    old_source_coordinates = sp.Matrix(
        (symbol_g, inverse_p, symbol_c, inverse_q)
    )
    source_form_in_new_coordinates = (
        old_source_coordinates.jacobian(new_variables).T
        * mixed_standard_form
        * old_source_coordinates.jacobian(new_variables)
    ).applyfunc(sp.expand)
    assert general_source_form == source_form_in_new_coordinates
    general_fiber_bracket = sp.expand(
        symbol_bracket(general_admission_f, general_admission_w)
        .subs(
            {
                symbol_p: inverse_p,
                symbol_q: general_inverse_q,
            },
            simultaneous=True,
        )
        .subs(capital_f, 0)
    )
    assert general_fiber_bracket == -symbol_g
    fiber_target_map = sp.Matrix(
        (symbol_g + symbol_c * capital_w, capital_w, symbol_c)
    )
    assert sp.expand(
        fiber_target_map.jacobian(
            (symbol_g, symbol_c, capital_w)
        ).det()
    ) == -1
    fiber_target_substitution = {
        u[0]: fiber_target_map[0],
        u[1]: 0,
        u[2]: fiber_target_map[1],
        u[3]: fiber_target_map[2],
    }
    target_form_on_admission_fiber = mixed_pullback_form.subs(
        fiber_target_substitution, simultaneous=True
    ).applyfunc(red)
    fiber_embedding_jacobian = sp.Matrix(
        (
            (1, capital_w, symbol_c),
            (0, 0, 0),
            (0, 0, 1),
            (0, 1, 0),
        )
    )
    pulled_fiber_form = (
        fiber_embedding_jacobian.T
        * target_form_on_admission_fiber
        * fiber_embedding_jacobian
    ).applyfunc(red)
    source_fiber_form = source_form_in_new_coordinates.extract(
        (1, 2, 3), (1, 2, 3)
    ).subs(capital_f, 0)
    fiber_form_difference = (
        pulled_fiber_form - source_fiber_form
    ).applyfunc(red)
    assert fiber_form_difference == sp.zeros(3), fiber_form_difference

    # The same statement can be read at the derivation level.  Under
    # (G,C,W)=(a-c*d,d,c), delta_0 is the bracket row of F.
    pulled_target_kernel = (
        fiber_delta(u[0] - u[2] * u[3]),
        fiber_delta(u[3]),
        fiber_delta(u[2]),
    )
    assert tuple(
        red(
            component.subs(
                fiber_target_substitution, simultaneous=True
            )
        )
        for component in pulled_target_kernel
    ) == admission_brackets

    # Solve the first transverse polynomial-automorphism compatibility
    # equation.  For
    #
    #   D(G)=B, D(C)=-A, D(W)=G,
    #
    # the free normal coefficient must satisfy D(s)=-C^2*G/3.  The obvious
    # primitive -2/C is not polynomial.  However I=C^4*G^2-12*C and
    # L=W-2/C^3 are D-invariant, and its pole cancels exactly against
    # I^2*L/144.
    def admission_delta(expression: sp.Expr) -> sp.Expr:
        return sp.expand(
            symbol_b * sp.diff(expression, symbol_g)
            - symbol_a * sp.diff(expression, symbol_c)
            + symbol_g * sp.diff(expression, capital_w)
        )

    admission_i = symbol_c**4 * symbol_g**2 - 12 * symbol_c
    admission_l = capital_w - 2 / symbol_c**3
    assert admission_delta(admission_i) == 0
    assert sp.cancel(admission_delta(admission_l)) == 0
    normal_s = sp.cancel(
        -2 / symbol_c - admission_i**2 * admission_l / 144
    )
    assert sp.denom(normal_s) == 1
    normal_s = sp.expand(normal_s)
    assert sp.expand(
        admission_delta(normal_s) + symbol_c**2 * symbol_g / 3
    ) == 0

    normal_d = sp.expand(
        symbol_g * symbol_c**3 / 3
        + symbol_c**4 * symbol_g * normal_s / 6
    )
    normal_c = sp.expand(-symbol_g * normal_s)
    normal_a = sp.expand(
        symbol_c**3 * symbol_g * capital_w / 3
        - 2 * symbol_c**2 * symbol_g**2 / 3
        + normal_s
        * (
            symbol_c**4 * symbol_g * capital_w
            - 2 * symbol_c**3 * symbol_g**2
            - 6 * symbol_c * symbol_g
            + 6
        )
        / 6
    )
    normal_column = sp.Matrix((normal_a, 1, normal_c, normal_d))
    full_first_jet_jacobian = sp.Matrix.hstack(
        normal_column, fiber_embedding_jacobian
    )
    source_form_at_zero = source_form_in_new_coordinates.subs(
        capital_f, 0
    )
    assert (
        full_first_jet_jacobian.T
        * target_form_on_admission_fiber
        * full_first_jet_jacobian
        - source_form_at_zero
    ).applyfunc(red) == sp.zeros(4)

    # In source fiber coordinates the induced first-order vector field is
    # divergence-free, so it is tangent to the volume-preserving polynomial
    # automorphism condition.  It is not itself locally nilpotent: C divides
    # its nonzero C-component.  Thus its single exponential does not provide
    # a finite resummation.
    normal_vector = sp.Matrix(
        (
            sp.expand(
                normal_a - capital_w * normal_d - symbol_c * normal_c
            ),
            normal_d,
            normal_c,
        )
    )
    normal_divergence = sp.expand(
        sp.diff(normal_vector[0], symbol_g)
        + sp.diff(normal_vector[1], symbol_c)
        + sp.diff(normal_vector[2], capital_w)
    )
    assert normal_divergence == 0
    assert normal_vector[1] != 0
    assert sp.rem(
        normal_vector[1], symbol_c, symbol_g, symbol_c, capital_w
    ) == 0

    # Advance one transverse order using the formal-flow second jet of the
    # divergence-free normal vector.  All tangent--tangent coefficients of
    # the order-F residual vanish.  The three remaining normal coefficients
    # are polynomial and q-independent.  Bezout gives an explicit polynomial
    # correction Z, unique modulo the characteristic field D.  Requiring Z
    # to preserve fiber volume asks for D(lambda)=-div(Z).  In the canonical
    # grading this equation has a unique infinite N-series and no polynomial
    # solution; hence this particular first jet cannot be completed by a
    # volume-compatible polynomial second jet.
    def normal_vector_delta(expression: sp.Expr) -> sp.Expr:
        return sp.expand(
            sum(
                normal_vector[index]
                * sp.diff(
                    expression,
                    (symbol_g, symbol_c, capital_w)[index],
                )
                for index in range(3)
            )
        )

    normal_target = sp.Matrix((normal_a, normal_c, normal_d))
    flow_second_target = sp.Matrix(
        [normal_vector_delta(component) / 2 for component in normal_target]
    ).applyfunc(sp.expand)
    target_form_first_variation = sp.zeros(4)
    for row in range(4):
        for column in range(4):
            variation = sp.diff(mixed_pullback_form[row, column], u[1])
            for target_index, normal_component in zip(
                (0, 2, 3), normal_target, strict=True
            ):
                variation += normal_component * sp.diff(
                    mixed_pullback_form[row, column], u[target_index]
                )
            target_form_first_variation[row, column] = red(
                variation.subs(
                    fiber_target_substitution, simultaneous=True
                )
            )

    first_jacobian_variation = sp.zeros(4)
    first_jacobian_variation[:, 0] = sp.Matrix(
        (
            2 * flow_second_target[0],
            0,
            2 * flow_second_target[1],
            2 * flow_second_target[2],
        )
    )
    fiber_symbols = (symbol_g, symbol_c, capital_w)
    first_jacobian_variation[0, 1:4] = sp.Matrix(
        [normal_a]
    ).jacobian(fiber_symbols)
    first_jacobian_variation[2, 1:4] = sp.Matrix(
        [normal_c]
    ).jacobian(fiber_symbols)
    first_jacobian_variation[3, 1:4] = sp.Matrix(
        [normal_d]
    ).jacobian(fiber_symbols)
    second_transverse_residual = (
        first_jacobian_variation.T
        * target_form_on_admission_fiber
        * full_first_jet_jacobian
        + full_first_jet_jacobian.T
        * target_form_first_variation
        * full_first_jet_jacobian
        + full_first_jet_jacobian.T
        * target_form_on_admission_fiber
        * first_jacobian_variation
        - sp.diff(source_form_in_new_coordinates, capital_f)
    ).applyfunc(red)
    assert second_transverse_residual + second_transverse_residual.T == sp.zeros(4)
    assert all(
        second_transverse_residual[left, right] == 0
        for left in range(1, 4)
        for right in range(left + 1, 4)
    )
    second_normal_residual = sp.Matrix(
        [second_transverse_residual[0, index] for index in range(1, 4)]
    )
    assert all(
        sp.diff(component, q) == 0 for component in second_normal_residual
    )
    assert [
        len(sp.Poly(component, *fiber_symbols).terms())
        for component in second_normal_residual
    ] == [9, 11, 7]

    second_bezout_correction = sp.Matrix(
        (
            -second_normal_residual[2] * symbol_r / 2,
            -second_normal_residual[2] * symbol_s / 2,
            (
                second_normal_residual[0] * symbol_r
                + second_normal_residual[1] * symbol_s
            )
            / 2,
        )
    ).applyfunc(sp.expand)
    assert (
        source_fiber_form * second_bezout_correction
        - second_normal_residual / 2
    ).applyfunc(sp.expand) == sp.zeros(3, 1)
    second_correction_divergence = sp.factor(
        sum(
            sp.diff(second_bezout_correction[index], fiber_symbols[index])
            for index in range(3)
        )
    )
    admission_n = symbol_c**3 * symbol_g**2
    admission_j = symbol_c**3 * capital_w
    second_divergence_polynomial = (
        admission_j * admission_n**4
        - 2 * admission_n**4
        - 42 * admission_j * admission_n**3
        - 24 * admission_n**3
        + 594 * admission_j * admission_n**2
        - 432 * admission_n**2
        - 1728 * admission_j * admission_n
        - 2592 * admission_n
        - 12960 * admission_j
        + 46656
    )
    assert sp.expand(
        second_correction_divergence
        - symbol_c * second_divergence_polynomial / 15552
    ) == 0

    # Give an all-degree non-solvability proof for the volume equation.
    # With wt(C)=2, wt(G)=-3, wt(W)=-6, D has weight 3 and the divergence
    # has weight 2.  The only relevant lambda-monomials have weight -1 and
    # are G*C*N^i*J^j.  Thus lambda=G*C*p(N,J).  The coefficient recurrence
    # for
    #
    #   L(p)=(N-6)p+N(N-12)p_N+3N(2-J)p_J=-H/2592
    #
    # reaches the nonzero constant p_5=-17/5388768.  Thereafter
    # p_n=n*p_(n-1)/(6*(2n+1)) is nonzero for every n, so p cannot be a
    # polynomial.
    recurrence_n, recurrence_j = sp.symbols("N J")
    recurrence_h = (
        recurrence_j * recurrence_n**4
        - 2 * recurrence_n**4
        - 42 * recurrence_j * recurrence_n**3
        - 24 * recurrence_n**3
        + 594 * recurrence_j * recurrence_n**2
        - 432 * recurrence_n**2
        - 1728 * recurrence_j * recurrence_n
        - 2592 * recurrence_n
        - 12960 * recurrence_j
        + 46656
    )
    recurrence_rhs = sp.Poly(
        -recurrence_h / 2592, recurrence_n
    )
    recurrence_coefficients = []
    for degree in range(6):
        forcing = recurrence_rhs.coeff_monomial(recurrence_n**degree)
        if degree == 0:
            coefficient = -forcing / 6
        else:
            previous = recurrence_coefficients[-1]
            coefficient = sp.expand(
                (
                    degree * previous
                    + 3
                    * (2 - recurrence_j)
                    * sp.diff(previous, recurrence_j)
                    - forcing
                )
                / (6 * (2 * degree + 1))
            )
        recurrence_coefficients.append(sp.factor(coefficient))
    assert recurrence_coefficients[5] == -sp.Rational(17, 5388768)

    # The invariant freedom in the first jet does not rescue the coordinate.
    # Write s=s_min+h with D(h)=0.  The second-order divergence is affine
    # linear in an arbitrary h (the quadratic coefficient vanishes before
    # imposing invariance).  The grading
    #
    #   wt(C)=2, wt(G)=-3, wt(W)=-6
    #
    # makes D homogeneous of weight three.  A shift of weight w changes the
    # divergence in weight w+4, so only ker(D)_{-2} can affect the baseline
    # weight-two obstruction.  From
    #
    #   ker(D)=K[I,S,U]/(I^3 U-S(S-3456)),
    #
    # its weight-minus-two piece is I^2*U*K[S].  For the top monomial
    # I^2*U*S^m, put r=m+2.  The exact leading J^r coefficient of p_5 is
    #
    # (-1728)^r (r-1)(3r-4)(3r-2)(3r-1)(9r^2-24r-5)/716636160,
    #
    # nonzero for every integer r>=2.  Descending on deg K excludes every
    # nonzero K; K=0 leaves the already nonzero baseline p_5.
    reduced_target_b_derivative = sp.diff(
        mixed_pullback_form, u[1]
    ).subs(fiber_target_substitution, simultaneous=True).applyfunc(red)
    reduced_source_form_derivative = sp.diff(
        source_form_in_new_coordinates, capital_f
    )
    reduced_n = symbol_c**3 * symbol_g**2
    reduced_j = symbol_c**3 * capital_w
    reduced_sigma_min = (
        -2 - (reduced_n - 12) ** 2 * (reduced_j - 2) / 144
    )

    def second_divergence_for_sigma(sigma: sp.Expr) -> sp.Expr:
        reduced_normal_a = sp.expand(
            symbol_g
            * (reduced_j * (sigma + 2) / 6 - sigma)
            + (sigma - (sigma + 2) * reduced_n / 3) / symbol_c
        )
        reduced_normal_c = sp.expand(-symbol_g * sigma / symbol_c)
        reduced_normal_d = sp.expand(
            symbol_g * symbol_c**3 * (sigma + 2) / 6
        )
        reduced_normal = sp.Matrix(
            (reduced_normal_a, reduced_normal_c, reduced_normal_d)
        )
        reduced_vector = sp.Matrix(
            (
                (sigma - (sigma + 2) * reduced_n / 3) / symbol_c,
                reduced_normal_d,
                reduced_normal_c,
            )
        )

        def reduced_vector_delta(expression: sp.Expr) -> sp.Expr:
            return sp.expand(
                sum(
                    reduced_vector[index]
                    * sp.diff(expression, fiber_symbols[index])
                    for index in range(3)
                )
            )

        reduced_second_target = sp.Matrix(
            [
                reduced_vector_delta(component) / 2
                for component in reduced_normal
            ]
        )
        reduced_form_variation = (
            reduced_target_b_derivative
            + target_form_on_admission_fiber.applyfunc(reduced_vector_delta)
        )
        reduced_jacobian_zero = sp.zeros(4)
        reduced_jacobian_zero[:, 0] = sp.Matrix(
            (
                reduced_normal_a,
                1,
                reduced_normal_c,
                reduced_normal_d,
            )
        )
        reduced_jacobian_zero[:, 1:4] = fiber_embedding_jacobian
        reduced_jacobian_one = sp.zeros(4)
        reduced_jacobian_one[:, 0] = sp.Matrix(
            (
                2 * reduced_second_target[0],
                0,
                2 * reduced_second_target[1],
                2 * reduced_second_target[2],
            )
        )
        reduced_jacobian_one[0, 1:4] = sp.Matrix(
            [reduced_normal_a]
        ).jacobian(fiber_symbols)
        reduced_jacobian_one[2, 1:4] = sp.Matrix(
            [reduced_normal_c]
        ).jacobian(fiber_symbols)
        reduced_jacobian_one[3, 1:4] = sp.Matrix(
            [reduced_normal_d]
        ).jacobian(fiber_symbols)
        reduced_residual = (
            reduced_jacobian_one.T
            * target_form_on_admission_fiber
            * reduced_jacobian_zero
            + reduced_jacobian_zero.T
            * reduced_form_variation
            * reduced_jacobian_zero
            + reduced_jacobian_zero.T
            * target_form_on_admission_fiber
            * reduced_jacobian_one
            - reduced_source_form_derivative
        )
        reduced_residual_normal = sp.Matrix(
            [reduced_residual[0, index] for index in range(1, 4)]
        )
        reduced_correction = sp.Matrix(
            (
                -reduced_residual_normal[2] * symbol_r / 2,
                -reduced_residual_normal[2] * symbol_s / 2,
                (
                    reduced_residual_normal[0] * symbol_r
                    + reduced_residual_normal[1] * symbol_s
                )
                / 2,
            )
        )
        return sp.expand(
            sum(
                sp.diff(reduced_correction[index], fiber_symbols[index])
                for index in range(3)
            )
        )

    shift_epsilon = sp.symbols("epsilon")
    arbitrary_shift = sp.Function("h")(
        symbol_g, symbol_c, capital_w
    )
    arbitrary_shift_divergence = second_divergence_for_sigma(
        reduced_sigma_min + shift_epsilon * symbol_c * arbitrary_shift
    )
    assert sp.diff(
        arbitrary_shift_divergence, shift_epsilon, 2
    ) == 0

    shift_power = sp.symbols("r", integer=True, positive=True)
    top_weight_minus_two_shift = (
        (reduced_n - 12) ** (3 * shift_power - 1)
        * reduced_j**shift_power
    )
    top_shift_response = sp.factor(
        sp.diff(
            second_divergence_for_sigma(
                reduced_sigma_min
                + shift_epsilon * top_weight_minus_two_shift
            ),
            shift_epsilon,
        ).subs(shift_epsilon, 0)
    )
    top_response_polynomial = (
        reduced_n**3
        - 30 * reduced_n**2
        + (108 * shift_power + 126) * reduced_n
        + 648 * shift_power
        + 432
    )
    expected_top_shift_response = sp.factor(
        -symbol_c
        * reduced_j ** (shift_power - 1)
        * (reduced_n - 12) ** (3 * shift_power - 2)
        * (
            reduced_j * top_response_polynomial
            + 6 * shift_power * reduced_n**3
            - 126 * shift_power * reduced_n**2
            + 648 * shift_power * reduced_n
        )
        / 108
    )
    assert sp.simplify(
        sp.powsimp(
            top_shift_response / expected_top_shift_response,
            force=True,
        )
    ) == 1

    top_forcing_coefficients = [
        648 * shift_power + 432,
        108 * shift_power + 126,
        -30,
        1,
    ]

    def shifted_power_coefficient(degree: int) -> sp.Expr:
        return (
            sp.binomial(3 * shift_power - 2, degree)
            * (-12) ** (3 * shift_power - 2 - degree)
        )

    def top_forcing_coefficient(degree: int) -> sp.Expr:
        return sp.simplify(
            sum(
                top_forcing_coefficients[index]
                * shifted_power_coefficient(degree - index)
                for index in range(4)
                if degree >= index
            )
            / 18
        )

    top_recurrence_coefficients = []
    for degree in range(6):
        forcing = top_forcing_coefficient(degree)
        if degree == 0:
            coefficient = -forcing / 6
        else:
            coefficient = (
                (degree - 3 * shift_power)
                * top_recurrence_coefficients[-1]
                - forcing
            ) / (6 * (2 * degree + 1))
        top_recurrence_coefficients.append(sp.factor(coefficient))
    top_p5_formula = (
        (-1728) ** shift_power
        * (shift_power - 1)
        * (3 * shift_power - 4)
        * (3 * shift_power - 2)
        * (3 * shift_power - 1)
        * (9 * shift_power**2 - 24 * shift_power - 5)
        / 716636160
    )
    assert sp.simplify(
        sp.powsimp(
            top_recurrence_coefficients[5] - top_p5_formula,
            force=True,
        )
    ) == 0
    mixed_linear_defects = {
        (row, column): homogeneous_part(mixed_raw_defect[row, column], 1)
        for row in range(4)
        for column in range(row + 1, 4)
        if homogeneous_part(mixed_raw_defect[row, column], 1) != 0
    }
    assert mixed_linear_defects == {
        (0, 2): -u[2],
        (0, 3): u[3],
        (1, 2): u[0],
    }

    stable_base_shear = sp.Matrix(
        (u[0] + u[2] * u[3], u[1], u[2], u[3])
    )
    stable_base_shear_inverse = sp.Matrix(
        (u[0] - u[2] * u[3], u[1], u[2], u[3])
    )
    stable_completion_shear = sp.Matrix(
        (u[0], u[1], u[2] - u[0] ** 2 / 2, u[3])
    )
    stable_completion_shear_inverse = sp.Matrix(
        (u[0], u[1], u[2] + u[0] ** 2 / 2, u[3])
    )
    for shear, inverse in (
        (stable_base_shear, stable_base_shear_inverse),
        (stable_completion_shear, stable_completion_shear_inverse),
    ):
        assert sp.expand(shear.jacobian(u).det()) == 1
        substitution = dict(zip(u, inverse, strict=True))
        assert sp.Matrix(
            [
                sp.expand(component.subs(substitution, simultaneous=True))
                for component in shear
            ]
        ) == sp.Matrix(u)

    mixed_correction_factors = [
        stable_base_shear,
        stable_completion_shear,
        elementary_u2_shear(2 * u[0] * u[1] ** 2),
        *cubic_correction_factors[:3],
    ]
    mixed_jet_map = sp.Matrix(
        [truncate_in_u(component, 8) for component in mixed_outputs]
    )
    for factor in mixed_correction_factors:
        mixed_jet_map = compose_jet(mixed_jet_map, factor, 8)
    mixed_corrected_jacobian = mixed_jet_map.jacobian(u)
    mixed_corrected_defect = (
        mixed_corrected_jacobian
        * target_poisson
        * mixed_corrected_jacobian.T
        - target_poisson
    ).applyfunc(lambda entry: truncate_in_u(entry, 4))
    for row in range(4):
        for column in range(row + 1, 4):
            for degree in range(1, 4):
                assert (
                    homogeneous_part(
                        mixed_corrected_defect[row, column], degree
                    )
                    == 0
                )
    mixed_quartic_defects = {
        (0, 1): -u[0] * u[1] * u[3] ** 2,
        (0, 2): (60 - 50 * q) * u[1] ** 4
        - 6 * u[0] ** 2 * u[1] * u[3],
        (0, 3): sp.Rational(1, 3) * u[0] * u[3] ** 3,
        (1, 2): 7 * u[0] * u[1] ** 2 * u[3],
        (1, 3): -sp.Rational(1, 3) * u[1] * u[3] ** 3,
        (2, 3): u[0] * u[1] * u[3] ** 2,
    }
    for pair, expected in mixed_quartic_defects.items():
        assert red(
            homogeneous_part(mixed_corrected_defect[pair], 4) - expected
        ) == 0

    mixed_quartic_correction_factors = [
        linear_hamiltonian_shear(
            u[0] + u[3], sp.Rational(1, 360), u[1], 5
        ),
        linear_hamiltonian_shear(
            u[0] - u[3], -sp.Rational(1, 360), u[1], 5
        ),
        linear_hamiltonian_shear(
            u[0] + 2 * u[3], -sp.Rational(1, 720), u[1], 5
        ),
        linear_hamiltonian_shear(
            u[0] - 2 * u[3], sp.Rational(1, 720), u[1], 5
        ),
        linear_hamiltonian_shear(
            u[3], sp.Rational(1, 12), u[1], 5
        ),
        elementary_u2_shear(
            -sp.Rational(7, 2) * u[0] ** 2 * u[1] ** 2 * u[3]
            + (12 - 10 * q) * u[1] ** 5
        ),
    ]
    for factor in mixed_quartic_correction_factors:
        mixed_jet_map = compose_jet(mixed_jet_map, factor, 8)
    mixed_quartic_corrected_jacobian = mixed_jet_map.jacobian(u)
    mixed_quartic_corrected_defect = (
        mixed_quartic_corrected_jacobian
        * target_poisson
        * mixed_quartic_corrected_jacobian.T
        - target_poisson
    ).applyfunc(lambda entry: truncate_in_u(entry, 5))
    for row in range(4):
        for column in range(row + 1, 4):
            for degree in range(1, 5):
                assert (
                    homogeneous_part(
                        mixed_quartic_corrected_defect[row, column], degree
                    )
                    == 0
                )
    mixed_quintic_defects = {
        (0, 1): -sp.Rational(1, 3) * u[0] ** 2 * u[3] ** 3,
        (0, 2): (12 - 10 * q) * u[1] ** 5
        + 4 * u[0] * u[1] ** 3 * u[3],
        (0, 3): sp.Integer(0),
        (1, 2): u[0] ** 2 * u[1] * u[3] ** 2
        - u[1] ** 4 * u[3],
        (1, 3): -sp.Rational(1, 6) * u[0] * u[3] ** 4,
        (2, 3): sp.Rational(1, 3) * u[0] ** 2 * u[3] ** 3,
    }
    for pair, expected in mixed_quintic_defects.items():
        assert red(
            homogeneous_part(mixed_quartic_corrected_defect[pair], 5)
            - expected
        ) == 0

    # Solve the quintic correction equation.  The sparse degree-six vector
    # is
    #
    # (-a^3*d^3/9, 0,
    #  -a^3*b*d^2/3+a*b^4*d+(6-5q)*b^6/3, a^2*d^4/12).
    #
    # Its (a,d)-part is generated by H=-a^3*d^4/36.  A seven-term Waring
    # decomposition of a^3*d^4 gives exact linear-Hamiltonian shears.
    mixed_quintic_correction_factors = [
        linear_hamiltonian_shear(
            u[0] + u[3], sp.Rational(13, 60480), 1, 7
        ),
        linear_hamiltonian_shear(
            u[0] - u[3], sp.Rational(13, 60480), 1, 7
        ),
        linear_hamiltonian_shear(
            u[0] + 2 * u[3], -sp.Rational(1, 15120), 1, 7
        ),
        linear_hamiltonian_shear(
            u[0] - 2 * u[3], -sp.Rational(1, 15120), 1, 7
        ),
        linear_hamiltonian_shear(
            u[0] + 3 * u[3], sp.Rational(1, 181440), 1, 7
        ),
        linear_hamiltonian_shear(
            u[0] - 3 * u[3], sp.Rational(1, 181440), 1, 7
        ),
        linear_hamiltonian_shear(
            u[0], -sp.Rational(1, 3240), 1, 7
        ),
        elementary_u2_shear(
            -sp.Rational(1, 3) * u[0] ** 3 * u[1] * u[3] ** 2
            + u[0] * u[1] ** 4 * u[3]
            + (6 - 5 * q) * u[1] ** 6 / 3
        ),
    ]
    for factor in mixed_quintic_correction_factors:
        mixed_jet_map = compose_jet(mixed_jet_map, factor, 8)
    mixed_quintic_corrected_jacobian = mixed_jet_map.jacobian(u)
    mixed_quintic_corrected_defect = (
        mixed_quintic_corrected_jacobian
        * target_poisson
        * mixed_quintic_corrected_jacobian.T
        - target_poisson
    ).applyfunc(lambda entry: truncate_in_u(entry, 6))
    for row in range(4):
        for column in range(row + 1, 4):
            for degree in range(1, 6):
                assert (
                    homogeneous_part(
                        mixed_quintic_corrected_defect[row, column], degree
                    )
                    == 0
                )
    mixed_sextic_defects = {
        (0, 1): -sp.Rational(27, 4) * u[0] ** 2 * u[1] ** 4
        - sp.Rational(9, 2) * u[0] * u[1] ** 4 * u[3]
        - sp.Rational(3, 4) * u[1] ** 4 * u[3] ** 2,
        (0, 2): 50 * q * u[0] * u[1] ** 4 * u[3]
        + sp.Rational(9, 4) * u[0] ** 2 * u[1] ** 4
        - sp.Rational(117, 2) * u[0] * u[1] ** 4 * u[3]
        - sp.Rational(27, 4) * u[1] ** 4 * u[3] ** 2,
        (0, 3): -9 * u[0] ** 3 * u[1] ** 3
        + 27 * u[0] ** 2 * u[1] ** 3 * u[3]
        + 9 * u[0] * u[1] ** 3 * u[3] ** 2
        + u[1] ** 3 * u[3] ** 3,
        (1, 2): (12 - 10 * q) * u[1] ** 5 * u[3]
        + sp.Rational(1, 3) * u[0] ** 3 * u[3] ** 3,
        (1, 3): sp.Rational(27, 4) * u[0] ** 2 * u[1] ** 4
        - sp.Rational(27, 2) * u[0] * u[1] ** 4 * u[3]
        - sp.Rational(9, 4) * u[1] ** 4 * u[3] ** 2,
        (2, 3): sp.Rational(27, 4) * u[0] ** 2 * u[1] ** 4
        + sp.Rational(9, 2) * u[0] * u[1] ** 4 * u[3]
        + sp.Rational(3, 4) * u[1] ** 4 * u[3] ** 2,
    }
    for pair, expected in mixed_sextic_defects.items():
        assert red(
            homogeneous_part(mixed_quintic_corrected_defect[pair], 6)
            - expected
        ) == 0

    # The sextic equation is again solvable.  A sparse degree-seven vector is
    # generated in the (a,d)-plane by
    #
    # b^4*(9*a^4/16-9*a^3*d/4-9*a^2*d^2/8-a*d^3/4),
    #
    # plus one elementary u2-shear.  Five fourth-power Hamiltonians integrate
    # the binary quartic exactly.
    mixed_sextic_correction_factors = [
        linear_hamiltonian_shear(
            u[0], sp.Rational(571, 384), u[1] ** 4, 4
        ),
        linear_hamiltonian_shear(
            u[0] + u[3], -sp.Rational(49, 32), u[1] ** 4, 4
        ),
        linear_hamiltonian_shear(
            u[0] + 2 * u[3], sp.Rational(59, 64), u[1] ** 4, 4
        ),
        linear_hamiltonian_shear(
            u[0] + 3 * u[3], -sp.Rational(37, 96), u[1] ** 4, 4
        ),
        linear_hamiltonian_shear(
            u[0] + 4 * u[3], sp.Rational(9, 128), u[1] ** 4, 4
        ),
        elementary_u2_shear(
            -sp.Rational(1, 12) * u[0] ** 4 * u[3] ** 3
            - (12 - 10 * q) * u[0] * u[1] ** 5 * u[3]
            - sp.Rational(27, 20) * u[1] ** 5 * u[3] ** 2
        ),
    ]
    for factor in mixed_sextic_correction_factors:
        mixed_jet_map = compose_jet(mixed_jet_map, factor, 8)
    mixed_sextic_corrected_jacobian = mixed_jet_map.jacobian(u)
    mixed_sextic_corrected_defect = (
        mixed_sextic_corrected_jacobian
        * target_poisson
        * mixed_sextic_corrected_jacobian.T
        - target_poisson
    ).applyfunc(lambda entry: truncate_in_u(entry, 7))
    for row in range(4):
        for column in range(row + 1, 4):
            for degree in range(1, 7):
                assert (
                    homogeneous_part(
                        mixed_sextic_corrected_defect[row, column], degree
                    )
                    == 0
                )
    mixed_septic_defects = {
        (0, 1): -2 * u[0] * u[1] ** 3 * u[3] ** 3,
        (0, 2): (126 - 91 * q) * u[1] ** 6 * u[3]
        + 24 * u[0] ** 2 * u[1] ** 3 * u[3] ** 2,
        (0, 3): sp.Rational(3, 2) * u[0] * u[1] ** 2 * u[3] ** 4,
        (1, 2): -sp.Rational(21, 2) * u[0] * u[1] ** 4 * u[3] ** 2,
        (1, 3): -sp.Rational(1, 2) * u[1] ** 3 * u[3] ** 4,
        (2, 3): 2 * u[0] * u[1] ** 3 * u[3] ** 3,
    }
    for pair, expected in mixed_septic_defects.items():
        assert red(
            homogeneous_part(mixed_sextic_corrected_defect[pair], 7)
            - expected
        ) == 0

    # General formal Darboux recurrence.  If D_m is the degree-m Poisson
    # defect, the corresponding form defect is E_m=W*D_m*W, W=-Pi.  Since
    # E_m is closed and homogeneous, alpha=i_Euler(E_m)/(m+2) satisfies
    # d(alpha)=E_m.  The vector V=Pi*alpha then has linearized Poisson change
    # -D_m.  This inductively kills every finite homogeneous order formally.
    standard_form = mixed_standard_form

    def radial_form_correction(
        defects: dict[tuple[int, int], sp.Expr], degree: int
    ) -> dict[str, object]:
        poisson_defect = sp.zeros(4)
        for (row, column), value in defects.items():
            poisson_defect[row, column] = value
            poisson_defect[column, row] = -value
        form_defect = (
            standard_form * poisson_defect * standard_form
        ).applyfunc(red)
        for left in range(4):
            for middle in range(left + 1, 4):
                for right in range(middle + 1, 4):
                    assert red(
                        sp.diff(form_defect[middle, right], u[left])
                        - sp.diff(form_defect[left, right], u[middle])
                        + sp.diff(form_defect[left, middle], u[right])
                    ) == 0
        radial_primitive = sp.Matrix(
            [
                red(
                    sum(
                        u[index] * form_defect[index, column]
                        for index in range(4)
                    )
                    / (degree + 2)
                )
                for column in range(4)
            ]
        )
        radial_vector = (target_poisson * radial_primitive).applyfunc(red)
        derivative = radial_vector.jacobian(u)
        linearized_change = (
            derivative * target_poisson
            + target_poisson * derivative.T
        ).applyfunc(red)
        assert all(
            red(linearized_change[row, column] + poisson_defect[row, column])
            == 0
            for row in range(4)
            for column in range(4)
        )
        return {
            "defect_degree": degree,
            "radial_denominator": degree + 2,
            "vector_component_term_counts": [
                len(sp.Poly(component, *u, q).terms())
                if component != 0
                else 0
                for component in radial_vector
            ],
            "closed_form_defect": True,
            "linearized_cancellation": True,
            "radial_vector": [str(component) for component in radial_vector],
        }

    radial_recurrence_checks = [
        radial_form_correction(mixed_sextic_defects, 6),
        radial_form_correction(mixed_septic_defects, 7),
    ]

    return {
        "coefficient_base": "Q[q]/(q^2-4q+6)",
        "map_term_counts": {
            name: len(sp.Poly(value, x, y, z, q).terms())
            for name, value in (("P", P), ("Q", Q), ("R", R))
        },
        "jacobian_determinant": "-1",
        "inverse_incidence_degree": sp.Poly(incidence, T).degree(),
        "inverse_incidence_derivative": "1-T*(Q-P*T)^2",
        "fiber_preserving_stabilization_no_go": {
            "factorization": "R=x*S with S nonconstant and S|_(x=0)=1",
            "cofactor_source_term_count": len(
                sp.Poly(r_cofactor, x, y, z, q).terms()
            ),
            "subgroup": (
                "source changes (a,b,c,d)->(h(a,b,d),lambda*c+G(a,b,d)) "
                "with h a polynomial base automorphism and lambda a unit"
            ),
            "conclusion": (
                "a Darboux trivialization in this subgroup would force "
                "R after a base automorphism to be affine-linear in d, "
                "hence R would be a coordinate; the displayed nontrivial "
                "factorization excludes this. Stable-variable mixing or a "
                "different target polarization is necessary"
            ),
        },
        "graph_form_determinant": "1",
        "affine_contact_screen": {
            "positive_degree_equations": len(contact_equations),
            "coefficient_dimension_over_Q": 6,
            "coefficient_rank_over_Q": contact_matrix.rank(),
            "kernel_dimension_over_K": 0,
            "conclusion": (
                "no nonzero constant target two-form has constant pullback; "
                "affine-contact symplectization with constant unit multiplier "
                "is excluded"
            ),
        },
        "tangent_normalization": {
            "linear_coordinates": "(u0,u1,u2,u3)=(z,y,e,x)",
            "unique_linear_defect": "D_12=u0",
            "polynomial_shear": "u1 -> u1-u0*u3",
            "shear_inverse": "u1 -> u1+u0*u3",
            "corrected_defect_order": 3,
            "cubic_defects": {
                f"D_{left}{right}": str(value)
                for (left, right), value in cubic_defects.items()
            },
            "degree_four_coordinate_shear_screen": (
                "the paired cubic D_01 cannot be removed by one layer of "
                "coordinate shears; a coupled or commutator correction is "
                "required"
            ),
            "tame_jet_correction": {
                "cubic_layer_factor_count": len(cubic_correction_factors),
                "quartic_layer_factor_count": len(
                    quartic_correction_factors
                ),
                "factor_type": (
                    "linear-Hamiltonian polynomial shears with displayed "
                    "polynomial inverses, plus elementary u2-shears"
                ),
                "corrected_through_degree": 4,
                "first_remaining_degree": 5,
                "quintic_defects": {
                    f"D_{left}{right}": str(value)
                    for (left, right), value in quintic_defects.items()
                },
                "scope": (
                    "finite tame jet normalization, not a finite global "
                    "Darboux trivialization"
                ),
            },
        },
        "stable_unimodular_charts": stable_charts,
        "natural_chart_test": natural_chart_defects,
        "stable_mixed_u2_jet": {
            "stable_coordinate": "B+x*e",
            "unique_linear_defects": {
                f"D_{left}{right}": str(value)
                for (left, right), value in mixed_linear_defects.items()
            },
            "correction_factor_count": len(mixed_correction_factors)
            + len(mixed_quartic_correction_factors)
            + len(mixed_quintic_correction_factors)
            + len(mixed_sextic_correction_factors),
            "stable_mixing_factors": [
                "u0 -> u0+u2*u3",
                "u2 -> u2-u0^2/2",
            ],
            "additional_factors": (
                "u2 -> u2+2*u0*u1^2, three cubic Hamiltonian shears, "
                "five quintic Hamiltonian shears, seven septic Hamiltonian "
                "shears, five quartic Hamiltonian shears, and three further "
                "u2-shears"
            ),
            "corrected_through_degree": 6,
            "first_remaining_degree": 7,
            "septic_defects": {
                f"D_{left}{right}": str(value)
                for (left, right), value in mixed_septic_defects.items()
            },
            "conclusion": (
                "the U_2 chart escapes the fiber-preserving no-go and has "
                "an explicit exact tame correction through sextic order; "
                "the nonzero septic row is the next coefficient of an "
                "all-order formal normalization rather than finite closure"
            ),
        },
        "formal_darboux_recurrence": {
            "form_defect": "E_m=(-Pi)*D_m*(-Pi)",
            "radial_primitive": (
                "alpha_m=i_Euler(E_m)/(m+2), with d(alpha_m)=E_m"
            ),
            "correction_vector": "V_(m+1)=Pi*alpha_m",
            "checks": radial_recurrence_checks,
            "conclusion": (
                "induction kills every fixed homogeneous order in the "
                "formal power-series topology"
            ),
            "polynomial_scope": (
                "the infinite formal composition need not terminate or "
                "resum to a polynomial automorphism; formal elimination is "
                "not R21 admission"
            ),
            "straight_line_moser_path": {
                "pfaffian": "1+s*(1-s)*Delta",
                "delta_nonzero": True,
                "delta_term_count": len(moser_delta_terms),
                "delta_minimum_source_degree": moser_delta_min_degree,
                "delta_maximum_source_degree": moser_delta_max_degree,
                "delta_leading_part": str(moser_delta_leading_part),
                "conclusion": (
                    "the canonical geometric series in "
                    "s*(1-s)*Delta is genuinely infinite; a special "
                    "resummation is still possible"
                ),
            },
            "constant_pfaffian_dilation_path": {
                "path": "K_s(u)=K(s*u)/s",
                "jacobian_determinant": "1",
                "trivializing_field": (
                    "V_s=-(Jac K_s)^(-1)*partial_s(K_s), polynomial "
                    "because det(Jac K_s)=1"
                ),
                "field_component_term_counts": dilation_field_term_counts,
                "field_component_source_degrees": (
                    dilation_field_source_degrees
                ),
                "field_component_parameter_degrees": (
                    dilation_field_parameter_degrees
                ),
                "formal_flow_identity": "K_s composed with H_s = id",
                "graph_generic_degree": 4,
                "degree_proof": (
                    "the quartic incidence is irreducible over K(P,Q,R), "
                    "and x,y,z are rationally reconstructed from P,Q,T"
                ),
                "time_one_scope": (
                    "H_1 is the formal inverse of the generically "
                    "degree-four U_2 graph map, so this path has no "
                    "polynomial time-one automorphism"
                ),
            },
            "target_b_component_obstruction": {
                "necessary_condition": (
                    "for H^*Omega=omega, the Hamiltonian kernel of the "
                    "source hypersurface H_1=0 must have no polynomial "
                    "slice"
                ),
                "fiber": "b=0",
                "primitive_kernel_at_b_zero": [
                    str(component) for component in fiber_kernel_at_zero
                ],
                "non_lnd_witness": (
                    "d divides delta(d)=d^4*(a-c*d)/6, which is nonzero"
                ),
                "no_slice_curve": {
                    "equations": [
                        "d^3*c-2=0",
                        "d^4*(a-c*d)^2-12*d=1",
                    ],
                    "coordinate": "y=d^2*(a-c*d)",
                    "induced_derivation": (
                        "delta(y)=(y^2-1)^2/144 on "
                        "K[y,1/(y^2-1)]"
                    ),
                    "time_form_residues_at_1_and_minus_1": [
                        str(value) for value in curve_time_residues
                    ],
                    "conclusion": (
                        "delta has no polynomial slice because a rational "
                        "derivative has zero residue at every pole"
                    ),
                },
                "conclusion": (
                    "the target-b component H_1 cannot be any elementary "
                    "coordinate lambda*u_i+F(other variables), nor any "
                    "coordinate having a polynomial Poisson mate modulo "
                    "H_1; this includes all b-fibration-preserving "
                    "normalizers, elementary b-shears, and the twenty-six "
                    "explicit tame factors"
                ),
                "degree_seven_fiber_admission": {
                    "source_coordinates": "(G,P,C,Q)",
                    "bezout_data": (
                        "A=G*C^4/6, B=G^2*C^3/3-1, "
                        "R=2*G^3*C^2/3, S=-(B+2), R*A+S*B=1"
                    ),
                    "coordinate": "F=A*(Q+G^2/2)-B*P",
                    "companion": "W=R*P+S*(Q+G^2/2)",
                    "coordinate_degree": sp.Poly(
                        admission_f, *symbol_variables
                    ).total_degree(),
                    "coordinate_term_count": len(
                        sp.Poly(admission_f, *symbol_variables).terms()
                    ),
                    "coordinate_jacobian": "1 for (G,P,C,Q)->(F,G,C,W)",
                    "inverse": [
                        "P=A*W+(B+2)*F",
                        "Q=B*W+R*F-G^2/2",
                    ],
                    "bezout_ansatz_classification": {
                        "general_shift": (
                            "Q'=Q+H(G,C) gives {F,W}|_(F=0)=-H_G"
                        ),
                        "fiber_condition": (
                            "R21 requires H_G=G, hence "
                            "H=G^2/2+k(C)"
                        ),
                        "c_shift": (
                            "k(C) is the ambient symplectic shear "
                            "Q->Q+k(C); the full source form in "
                            "(F,G,C,W) is unchanged"
                        ),
                        "all_bezout_pairs": (
                            "(R_T,S_T)=(R+B*T,S-A*T), so "
                            "W_T=W-T*F"
                        ),
                        "conclusion": (
                            "the second-order exclusion applies to the "
                            "complete reciprocal-compatible Bezout ansatz, "
                            "not only the displayed companion choice"
                        ),
                    },
                    "fiber_map": "(G,C,W)->(a,c,d)=(G+C*W,W,C)",
                    "fiber_map_jacobian": "-1",
                    "fiber_bracket_row": [
                        str(value) for value in admission_brackets
                    ],
                    "result": (
                        "the pullback of the complete R21 b=0 two-form "
                        "equals the restricted standard source form exactly"
                    ),
                    "normal_jet": {
                        "invariants": (
                            "I=C^4*G^2-12*C, L=W-2/C^3"
                        ),
                        "pole_cancellation": (
                            "s=-2/C-I^2*L/144 is polynomial"
                        ),
                        "s_polynomial": str(normal_s),
                        "transport_equation": "D(s)=-C^2*G/3",
                        "component_degrees": [
                            sp.Poly(
                                component, symbol_g, symbol_c, capital_w
                            ).total_degree()
                            for component in (normal_a, normal_c, normal_d)
                        ],
                        "full_form_match_on_fiber": True,
                        "source_fiber_vector_divergence": str(
                            normal_divergence
                        ),
                        "locally_nilpotent": False,
                        "non_lnd_witness": (
                            "C divides the nonzero C-component"
                        ),
                        "canonical_second_order_test": {
                            "flow_choice": (
                                "the F^2/2 formal-flow jet of the displayed "
                                "divergence-free normal vector"
                            ),
                            "tangent_tangent_residual_zero": True,
                            "normal_residual_term_counts": [
                                len(
                                    sp.Poly(
                                        component, *fiber_symbols
                                    ).terms()
                                )
                                for component in second_normal_residual
                            ],
                            "normal_residual_degrees": [
                                sp.Poly(
                                    component, *fiber_symbols
                                ).total_degree()
                                for component in second_normal_residual
                            ],
                            "coefficient_base": "Q (independent of q)",
                            "bezout_correction_degrees": [
                                sp.Poly(
                                    component, *fiber_symbols
                                ).total_degree()
                                for component in second_bezout_correction
                            ],
                            "bezout_correction_divergence": str(
                                second_correction_divergence
                            ),
                            "volume_equation": (
                                "D(lambda)=-div(Z), with every solution of "
                                "the form equation equal to Z+lambda*D"
                            ),
                            "grading": (
                                "wt(C)=2, wt(G)=-3, wt(W)=-6; the required "
                                "weight forces lambda=G*C*p(N,J), "
                                "N=C^3*G^2, J=C^3*W"
                            ),
                            "coefficient_recurrence_p_0_through_p_5": [
                                str(value)
                                for value in recurrence_coefficients
                            ],
                            "nontermination_witness": "p_5=-17/5388768",
                            "conclusion": (
                                "the canonical minimal first jet has no "
                                "polynomial volume-compatible second jet"
                            ),
                            "all_invariant_shifts": {
                                "general_first_jet": (
                                    "s=s_min+h with h in ker(D)"
                                ),
                                "affine_linearity": (
                                    "the second-order divergence has zero "
                                    "quadratic variation for arbitrary h"
                                ),
                                "grading": (
                                    "wt(C)=2, wt(G)=-3, wt(W)=-6; a "
                                    "weight-w shift changes the divergence "
                                    "in weight w+4"
                                ),
                                "relevant_kernel_piece": (
                                    "ker(D)_(-2)=I^2*U*K[S]"
                                ),
                                "top_monomial": (
                                    "for I^2*U*S^m put r=m+2"
                                ),
                                "top_p_5_coefficient": str(top_p5_formula),
                                "nonvanishing_range": "every integer r>=2",
                                "conclusion": (
                                    "descending on deg(K) excludes every "
                                    "nonzero I^2*U*K[S] shift; K=0 leaves "
                                    "p_5=-17/5388768. No polynomial "
                                    "invariant first-jet shift admits a "
                                    "volume-compatible second completion"
                                ),
                            },
                        },
                    },
                    "scope": (
                        "exact polynomial hypersurface admission and "
                        "volume-compatible first transverse jet, followed "
                        "by an all-invariant second-order obstruction. No "
                        "polynomial R21 normalizer can have this displayed "
                        "F as its target-b component"
                    ),
                },
                "sharp_tame_no_slice_family": {
                    "coordinates": (
                        "G=b+a*d, F_k=a+c^k*G; "
                        "(F_k,G,c,d) is a tame determinant-one automorphism"
                    ),
                    "fiber_invariant_data": (
                        "t=G*c^k, I=t^2-2*c, J=1+c^k*d; "
                        "D_k(t)=-c^k, D_k(I)=0, "
                        "D_k(J)=-k*t*c^(k-1)*J"
                    ),
                    "curve": "I=1, J=0",
                    "curve_coordinate_ring": "K[t,1/(t^2-1)]",
                    "time_form": "-2^k*dt/(t^2-1)^k",
                    "verified_rows": tame_no_slice_rows,
                    "all_degree_conclusion": (
                        "the residues are nonzero for every k>=1, so this "
                        "is an infinite tame family of no-slice coordinates; "
                        "the degree-three row shows the target-coordinate "
                        "obstruction is sharp in degree"
                    ),
                    "r21_comparison": (
                        "the R21 invariant curve has pole order 2 but "
                        "transverse logarithmic weight 3; this family ties "
                        "both values to k. The complete constant ring has "
                        "intrinsic Danielewski exponent k, whereas R21 has "
                        "exponent 3; hence k=3 is forced globally, but the "
                        "generic time-form divisor forces k=2. No member "
                        "of the family can be the target coordinate"
                    ),
                    "bounded_semi_invariant_census": {
                        "degree_bound": eigenspace_degree_bound,
                        "weights": list(eigenspace_weights),
                        "tame_k_2_dimensions": (
                            tame_two_weight_dimensions
                        ),
                        "r21_dimensions": target_weight_dimensions,
                        "first_extra_irreducible_weights": {
                            "tame_k_2": 2,
                            "r21": 3,
                        },
                        "status": (
                            "exact bounded experiment only; a new "
                            "irreducible weight-one Darboux polynomial "
                            "appears in degree eight for the tame k=2 row"
                        ),
                    },
                    "degree_eight_correction": {
                        "polynomial": str(tame_two_degree_eight_prime),
                        "weight": 1,
                        "irreducibility_check": (
                            "primitive and linear in d over K[G,c]"
                        ),
                        "conclusion": (
                            "the degree-seven weight census does not "
                            "classify all Darboux primes"
                        ),
                    },
                    "polynomial_constant_ring_obstruction": {
                        "tame_k_2": {
                            "generators": (
                                "I=c*(G^2*c^3-2), "
                                "S=(G^2*c^3-2)^2*(1+c^2*d), "
                                "U=S*(S-4)/I^2"
                            ),
                            "relation": "I^2*U=S*(S-4)",
                            "u_term_count": len(
                                sp.Poly(
                                    tame_two_u,
                                    calibration_g,
                                    calibration_c,
                                    calibration_d,
                                ).terms()
                            ),
                            "u_total_degree": sp.Poly(
                                tame_two_u,
                                calibration_g,
                                calibration_c,
                                calibration_d,
                            ).total_degree(),
                        },
                        "r21": {
                            "generators": (
                                "I=d*(d^3*(a-c*d)^2-12), "
                                "S=(d^3*(a-c*d)^2-12)^3*(d^3*c-2), "
                                "U=S*(S-3456)/I^3"
                            ),
                            "relation": "I^3*U=S*(S-3456)",
                            "u_term_count": len(
                                sp.Poly(target_u, *fiber_variables).terms()
                            ),
                            "u_total_degree": sp.Poly(
                                target_u, *fiber_variables
                            ).total_degree(),
                        },
                        "all_k_formula": (
                            "ker(D_k)=K[I,S,U]/"
                            "(I^k*U-S*(S-(-2)^k))"
                        ),
                        "intrinsic_exponent": (
                            "for k>=2 the Makar-Limanov invariant is K[I] "
                            "and the canonical plinth ideal is (I^k), so "
                            "the Danielewski exponent k is preserved by "
                            "ring isomorphism"
                        ),
                        "conclusion": (
                            "the constant rings force k=3, while the "
                            "generic time-form divisor forces k=2; the "
                            "entire tame no-slice family is excluded"
                        ),
                    },
                    "dense_open_k_2_conjugacy": {
                        "target_coordinates": "t=a-c*d, I=d^4*t^2-12*d",
                        "map": (
                            "G=t/36, C=6*d, "
                            "W=-1/(36*d^2)+I*(c-2/d^3)/864"
                        ),
                        "w_polynomial": str(dense_tame_w),
                        "intertwining": "delta_0 composed with Psi=(-1/36)*D_2",
                        "jacobian": "-I/5184",
                        "expanded_jacobian": str(
                            dense_conjugacy_jacobian
                        ),
                        "constant_field_map": (
                            "I_2 maps to I, L_2 maps to I*L_R21/864"
                        ),
                        "scope": (
                            "an exact isomorphism after localizing I; the "
                            "Jacobian factor I prevents a global polynomial "
                            "fiber isomorphism"
                        ),
                    },
                },
            },
        },
        "rank_two_status": (
            "not admitted: the exact graph form is polynomial and the "
            "isolated reciprocal powers are stably unimodular; affine "
            "contact and every fiber-preserving stabilization are excluded, "
            "the stable-mixed chart is formally Darboux-trivial to all "
            "orders, but no finite polynomial automorphism trivializing the "
            "coupled graph form is certified. A new degree-seven coordinate "
            "does pass the complete b=0 two-form gate and has a polynomial, "
            "volume-compatible first transverse jet, but its canonical "
            "normal vector is not locally nilpotent"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    certificate = {
        "format": "dc2-higher-nilpotence-r21-frontier-v13",
        "higher_nilpotence_family": verify_higher_nilpotence_family(),
        "r21_frontier": verify_r21_frontier(),
        "conclusion": (
            "N^2 != 0 is achieved by an exact regular-index-four family, "
            "but triangular rigidity makes N^2 constant and the family is "
            "an automorphism/Moyal-flat control. For R21, affine contact is "
            "excluded and exact tame shears normalize the graph through "
            "degree four, but the whole fiber-preserving stabilization "
            "subgroup is excluded by the reducible factorization of R. "
            "The stable-mixed U_2 chart survives and is normalized through "
            "degree six with an explicit septic remainder. Euler homotopy "
            "then eliminates every homogeneous defect to all formal orders. "
            "The exact straight-line Moser Pfaffian has nonzero Delta, so "
            "its canonical geometric series does not terminate. "
            "Source dilation gives a constant-Pfaffian path with a "
            "polynomial trivializing field, but its time-one map is the "
            "nonpolynomial inverse of the generic-degree-four graph map. "
            "The exact fiber-kernel derivation also excludes every "
            "polynomial normalizer whose target-b component has a "
            "polynomial Poisson mate on its zero fiber. In particular, "
            "that component cannot be an elementary coordinate; this "
            "includes the complete twenty-six-factor correction subgroup. "
            "An explicit infinite tame family F_k=a+c^k*(b+a*d) has no "
            "polynomial Poisson mate, showing that the necessary condition "
            "is attainable and already sharp in degree three. Its curve "
            "pole order and transverse weight are both k, whereas R21 "
            "requires the split signature (2,3). The complete polynomial "
            "constant rings are Danielewski surfaces with exponents k and "
            "3, while the generic time-form divisors have pole orders k "
            "and 2. These two intrinsic requirements exclude every F_k. "
            "For k=2 there is an exact polynomial conjugacy on I!=0, but "
            "its Jacobian is -I/5184, locating the global obstruction on "
            "the affine-modification divisor I=0. A new degree-seven "
            "Bezout coordinate F=A*(Q+G^2/2)-B*P now crosses that divisor: "
            "its zero hypersurface carries exactly the R21 b=0 two-form. "
            "The first transverse automorphism equation has a polynomial "
            "solution because the apparent pole -2/C cancels against the "
            "constant-ring invariant I^2*L/144. The resulting normal jet "
            "is divergence-free but not locally nilpotent, so a single "
            "finite exponential does not finish the construction. At the "
            "next order its formal-flow jet leaves only three q-independent "
            "normal defects. Bezout solves the form equation polynomially, "
            "but the required volume correction D(lambda)=-div(Z) has a "
            "nonterminating coefficient recurrence, beginning with the "
            "nonzero tail p_5=-17/5388768. Thus the canonical minimal first "
            "jet has no polynomial second completion. More generally, the "
            "obstruction is affine-linear in every invariant shift. Grading "
            "reduces the only relevant freedom to I^2*U*K[S], and the top "
            "J^(m+2) coefficient for I^2*U*S^m is nonzero for every m>=0. "
            "Hence no polynomial invariant shift works: this degree-seven "
            "target-b coordinate is excluded at second order. Replacing "
            "Q+G^2/2 by Q+H(G,C) forces H=G^2/2+k(C), an ambient "
            "symplectic shear, while every other Bezout companion is "
            "W->W-T*F. Thus the whole reciprocal-compatible Bezout ansatz "
            "is excluded, not only one companion choice. "
            "Whether that formal normalizer resums to a finite polynomial "
            "automorphism remains the coupled graph-Darboux problem."
        ),
    }
    rendered = json.dumps(certificate, indent=2, sort_keys=True) + "\n"
    digest = hashlib.sha256(rendered.encode()).hexdigest()
    if args.output:
        output = args.output.resolve()
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(rendered)
    else:
        print(rendered, end="")

    print("PASS: constructed an all-degree Cayley-integrable N^2 != 0 family")
    print("PASS: triangular row closure forces the variable part of N^2 to vanish")
    print("PASS: verified the exact R21 Keller map and polynomial graph form")
    print("PASS: excluded the constant-two-form affine-contact lift for R21")
    print("PASS: a polynomial shear raises the first R21 graph defect from degree 1 to 3")
    print("PASS: ten exact tame shears remove the R21 graph defect through degree 4")
    print("PASS: reducibility of R excludes the full fiber-preserving stabilization subgroup")
    print("PASS: the stable-mixed U_2 chart has an exact tame correction through degree 6")
    print("PASS: verified the Euler-homotopy formal Darboux recurrence in degrees 6 and 7")
    print("PASS: proved the straight-line Moser geometric series is nonterminating")
    print("PASS: the constant-Pfaffian dilation path has an exact polynomial field")
    print("PASS: excluded every elementary target-b coordinate for a Darboux normalizer")
    print("PASS: the b=0 fiber kernel has no polynomial slice")
    print("PASS: constructed an infinite tame family of no-slice coordinates")
    print("PASS: compared tame and R21 semi-invariants through degree seven")
    print("PASS: computed the tame and R21 polynomial constant rings")
    print("PASS: excluded the full tame no-slice family by split invariants")
    print("PASS: conjugated the k=2 and R21 fibers exactly on I nonzero")
    print("PASS: a degree-seven coordinate admits the complete R21 b=0 form")
    print("PASS: its first transverse volume-compatible jet is polynomial")
    print("PASS: its canonical second volume equation has no polynomial solution")
    print("PASS: grading excludes every invariant shift of that first jet")
    print("PASS: classified and excluded the full reciprocal Bezout ansatz")
    print("PASS: verified four reciprocal-power stable unimodular charts")
    print("OPEN: no polynomial Darboux trivialization of the R21 graph is certified")
    print(f"SHA256: {digest}")


if __name__ == "__main__":
    main()

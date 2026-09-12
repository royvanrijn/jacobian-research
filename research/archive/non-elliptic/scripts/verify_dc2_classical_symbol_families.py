#!/usr/bin/env python3
"""Exact construction checks for four two-parameter DC_2 symbol families.

The families deliberately separate a genuine marked-root search branch from
three controls:

* the degree-six weighted marked-root family;
* a two-parameter target-symplectic orbit of the normalized-factorization
  (``c=-9``) completion;
* a low-degree Hamiltonian suspension; and
* a rank-one fibre shear over a central parameter, including its connection
  coordinate.

The last two families are polynomial symplectic automorphisms with exact
Weyl lifts.  They are useful zero-obstruction controls, not DC_2
counterexample candidates.  The normalized-factorization orbit transports
the known order-five obstruction and is likewise not a new moduli branch.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import sympy as sp
from sympy.polys.domains import QQ


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


from jcsearch.deformation_complex import relative_fiber_connection_complex
from verify_degree_six_relative_quantization_obstruction import degree_six_pair


def canonical_bracket(
    left: sp.Expr,
    right: sp.Expr,
    positions: tuple[sp.Symbol, sp.Symbol],
    momenta: tuple[sp.Symbol, sp.Symbol],
) -> sp.Expr:
    """Poisson bracket with ``{p_i,q_j}=delta_ij``."""

    return sp.expand(
        sum(
            sp.diff(left, momentum) * sp.diff(right, position)
            - sp.diff(left, position) * sp.diff(right, momentum)
            for position, momentum in zip(positions, momenta)
        )
    )


def pi_power_constant(
    left: sp.Expr,
    right: sp.Expr,
    power: int,
    positions: tuple[sp.Symbol, sp.Symbol],
    momenta: tuple[sp.Symbol, sp.Symbol],
) -> sp.Expr:
    """Apply the constant Poisson bidifferential ``power`` times."""

    variables = positions + momenta
    variable_index = {variable: index for index, variable in enumerate(variables)}
    terms = []
    for position, momentum in zip(positions, momenta):
        terms.append((momentum, position, 1))
        terms.append((position, momentum, -1))

    states = {((0,) * 4, (0,) * 4): sp.Integer(1)}
    for _ in range(power):
        next_states: dict[tuple[tuple[int, ...], tuple[int, ...]], sp.Expr] = {}
        for (left_degrees, right_degrees), coefficient in states.items():
            for left_variable, right_variable, sign in terms:
                new_left = list(left_degrees)
                new_right = list(right_degrees)
                new_left[variable_index[left_variable]] += 1
                new_right[variable_index[right_variable]] += 1
                key = (tuple(new_left), tuple(new_right))
                next_states[key] = next_states.get(key, 0) + sign * coefficient
        states = next_states

    value = 0
    for (left_degrees, right_degrees), coefficient in states.items():
        left_derivative = left
        right_derivative = right
        for variable, degree in zip(variables, left_degrees):
            if degree:
                left_derivative = sp.diff(left_derivative, variable, degree)
        for variable, degree in zip(variables, right_degrees):
            if degree:
                right_derivative = sp.diff(right_derivative, variable, degree)
        value += coefficient * left_derivative * right_derivative
    return sp.expand(value)


def verify_canonical_tuple(
    outputs: tuple[sp.Expr, sp.Expr, sp.Expr, sp.Expr],
    positions: tuple[sp.Symbol, sp.Symbol],
    momenta: tuple[sp.Symbol, sp.Symbol],
) -> None:
    """Check the six canonical brackets in ``(q1,q2,p1,p2)`` order."""

    q1, q2, p1, p2 = outputs
    expected = {
        (0, 1): 0,
        (0, 2): -1,
        (0, 3): 0,
        (1, 2): 0,
        (1, 3): -1,
        (2, 3): 0,
    }
    for (left, right), value in expected.items():
        assert canonical_bracket(
            outputs[left], outputs[right], positions, momenta
        ) == value


def verify_exact_moyal_relations(
    outputs: tuple[sp.Expr, sp.Expr, sp.Expr, sp.Expr],
    positions: tuple[sp.Symbol, sp.Symbol],
    momenta: tuple[sp.Symbol, sp.Symbol],
) -> None:
    """Check that every higher odd Moyal term through the degree cutoff is zero."""

    for left in range(4):
        for right in range(left + 1, 4):
            for power in (3, 5, 7):
                assert (
                    pi_power_constant(
                        outputs[left],
                        outputs[right],
                        power,
                        positions,
                        momenta,
                    )
                    == 0
                )


def marked_root_family() -> dict[str, object]:
    field = QQ.frac_field("sigma", "tau")
    sigma, tau = field.gens
    S, T = degree_six_pair(field, sigma, tau)
    assert len(S) == 75
    assert len(T) == 59
    return {
        "key": "MR6",
        "construction": "degree-six weighted marked-root incidence",
        "parameter_ring": "Q[sigma,tau,1/sigma] on the exact-degree chart",
        "symbol_terms": {"S": len(S), "T": len(T)},
        "order_five_locus": (
            "generic obstruction; reconstructed quartic closed component "
            "is checked by verify_degree_six_order_five_survivor.py"
        ),
        "candidate_status": "genuine noninjective classical branch",
    }


def normalized_factorization_family() -> dict[str, object]:
    # Abstract target coordinates for the exact c=-9 completion.  The two
    # linear shears are exact symplectic target automorphisms, so composing
    # them with the normalized-factorization completion neither changes its
    # collision degree nor removes its transported restricted obstruction.
    R, T, D, S, rho, eta = sp.symbols("R T D S rho eta")
    outputs = (R, T, D + rho * R, S + eta * T)
    verify_canonical_tuple(outputs, (R, T), (D, S))
    verify_exact_moyal_relations(outputs, (R, T), (D, S))
    return {
        "key": "NF3",
        "construction": (
            "two-parameter target-symplectic orbit of the normalized "
            "linear-quadratic factorization completion"
        ),
        "parameter_ring": "Q[rho,eta]",
        "target_shear": ["D -> D+rho*R", "S -> S+eta*T"],
        "order_five_locus": "empty in the transported c=-9 filtration",
        "order_five_value": "-49",
        "candidate_status": "gauge-equivalent obstruction control",
    }


def hamiltonian_suspension_family() -> dict[str, object]:
    q1, q2, p1, p2, lam, mu = sp.symbols("q1 q2 p1 p2 lambda mu")
    hamiltonian = lam * q1**3 / 3 + mu * q1 * q2**2
    outputs = (
        q1,
        q2,
        p1 + sp.diff(hamiltonian, q1),
        p2 + sp.diff(hamiltonian, q2),
    )
    verify_canonical_tuple(outputs, (q1, q2), (p1, p2))
    verify_exact_moyal_relations(outputs, (q1, q2), (p1, p2))
    return {
        "key": "HS3",
        "construction": "low-degree Hamiltonian suspension",
        "parameter_ring": "Q[lambda,mu]",
        "hamiltonian": "lambda*q1^3/3+mu*q1*q2^2",
        "order_five_locus": "all of A^2",
        "order_seven_locus": "all of A^2",
        "candidate_status": "exact Weyl automorphism control",
    }


def central_fibre_family() -> dict[str, object]:
    v, R, P, U, lam, mu = sp.symbols("v R P U lambda mu")
    potential = lam * R * v**3 / 3 + mu * R**2 * v**2 / 2
    outputs = (
        v,
        R,
        P + sp.diff(potential, v),
        U + sp.diff(potential, R),
    )
    verify_canonical_tuple(outputs, (v, R), (P, U))
    verify_exact_moyal_relations(outputs, (v, R), (P, U))
    return {
        "key": "RF3",
        "construction": "rank-one fibre map over central R with connection",
        "parameter_ring": "Q[lambda,mu]",
        "potential": "lambda*R*v^3/3+mu*R^2*v^2/2",
        "connection_coordinate": "U+lambda*v^3/3+mu*R*v^2",
        "order_five_locus": "all of A^2",
        "order_seven_locus": "all of A^2",
        "candidate_status": "exact Weyl automorphism control",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    relative = relative_fiber_connection_complex().complex
    assert relative.dimensions == (35, 60, 26)
    assert relative.ranks == (34, 26)
    assert relative.cohomology_dimensions == (1, 0, 0)
    assert relative.dual_obstruction_cocycles() == ()

    families = [
        marked_root_family(),
        normalized_factorization_family(),
        hamiltonian_suspension_family(),
        central_fibre_family(),
    ]
    assert len({family["key"] for family in families}) == 4

    certificate = {
        "scope": (
            "four explicit two-parameter classical families; exact Moyal "
            "control through order seven for the tame controls; no DC_2 claim"
        ),
        "relative_fibre_connection_complex": {
            "dimensions": list(relative.dimensions),
            "ranks": list(relative.ranks),
            "cohomology_dimensions": list(relative.cohomology_dimensions),
            "fitting_0_obstruction_module": "unit ideal",
        },
        "families": families,
    }
    if args.output:
        output = args.output.resolve()
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(json.dumps(certificate, indent=2, sort_keys=True) + "\n")

    print("PASS: constructed four explicit two-parameter symplectic families")
    print("PASS: the canonical relative complex is 35 -> 60 -> 26 and exact in H^2")
    print("PASS: HS3 and RF3 have zero Moyal cocycles through order seven")
    print("PASS: NF3 is a transported c=-9 obstruction orbit, not new moduli")
    print("SCOPE: only MR6 advances to a nonautomorphism PBW search")


if __name__ == "__main__":
    main()

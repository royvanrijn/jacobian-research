#!/usr/bin/env sage-python
"""Upgrade the marked G720-S0128 point to a formal Z_7 family.

The tangent artifact has 55 equations in 46 variables and a unit rank-45
minor, so a formal implicit-function argument applies only after explaining
the ten nonpivot equations.  For a short Weierstrass cubic

    f(X) = X^3 + A X + B,  D = 4 A^3 + 27 B^2,
    N = 2 A X + 3 B,

the exact identity

    8 A^3 f(X) = D (N-B) - 9 B N^2 + N^3

shows, wherever A is a unit, that ord(D)>=6 and ord(N),ord(Y)>=d force
ord(Y^2-f(X))>=min(6,2d).  Applied at 0, 1, and infinity with depths 1, 1,
and 3, the nonidentity-section residual is divisible by t^2(t-1)^2 and has
degree at most six.  It therefore has only three free quotient coefficients;
the pivot equations at residual coefficients 2,3,4 kill all three.

Consequently the full 55-equation germ equals the 45-equation pivot germ.
The stored unit minor then gives a one-parameter formal Z_7 solution by the
multivariate formal implicit-function theorem.  This is local 7-adic/formal
evidence, not a rational family over Q.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_LIFT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-golay-det720-3a5-marked-gf7-lift-v1.json"
)
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-golay-det720-3a5-formal-smoothness-v1.json"
)


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def universal_identity(base_ring) -> bool:
    ring = PolynomialRing(base_ring, names=("A", "B", "X", "Y"))
    A, B, X, unused_Y = ring.gens()
    discriminant = 4 * A**3 + 27 * B**2
    node = 2 * A * X + 3 * B
    cubic = X**3 + A * X + B
    return bool(
        8 * A**3 * cubic
        == discriminant * (node - B) - 9 * B * node**2 + node**3
    )


def localized_jet_membership(base_ring, depth: int) -> dict:
    names = [f"{letter}{index}" for letter in "abxy" for index in range(6)]
    names.append("z")
    coefficient_ring = PolynomialRing(base_ring, names=names, order="degrevlex")
    variables = list(coefficient_ring.gens())
    a = variables[0:6]
    b = variables[6:12]
    x = variables[12:18]
    y = variables[18:24]
    z = variables[24]
    function_ring = PolynomialRing(coefficient_ring, "u")
    u = function_ring.gen()

    def polynomial(coefficients):
        return sum(value * u**index for index, value in enumerate(coefficients))

    A = polynomial(a)
    B = polynomial(b)
    X = polynomial(x)
    Y = polynomial(y)
    discriminant = 4 * A**3 + 27 * B**2
    node = 2 * A * X + 3 * B
    residual = Y**2 - X**3 - A * X - B
    generators = (
        [discriminant[index] for index in range(6)]
        + [node[index] for index in range(depth)]
        + y[:depth]
        + [z * a[0] - 1]
    )
    ideal = coefficient_ring.ideal(generators)
    basis = ideal.groebner_basis()
    target_orders = range(2 * depth)
    zero_remainders = [
        bool(ideal.reduce(residual[index]) == 0) for index in target_orders
    ]
    if not all(zero_remainders):
        raise ArithmeticError("localized component-jet implication failed")
    return {
        "depth": depth,
        "discriminant_order": 6,
        "forced_residual_order": 2 * depth,
        "localized_at": "A(0)",
        "groebner_basis_size": len(basis),
        "checked_residual_coefficients": list(target_orders),
        "all_remainders_zero": True,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lift", type=Path, default=DEFAULT_LIFT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()
    lift_path = arguments.lift.resolve()
    output_path = arguments.output.resolve()
    lift = json.loads(lift_path.read_text())
    if lift.get("schema") != "elkies-k3.golay-det720-3a5-marked-gf7-lift.v1":
        raise ValueError("unexpected marked-lift schema")
    certificate = lift["jacobian_certificate"]
    if (
        int(certificate["rank_mod_7"]) != 45
        or int(certificate["tangent_dimension"]) != 1
        or int(certificate["pivot_minor_determinant_mod_7"]) % 7 == 0
    ):
        raise ValueError("marked point lacks the required unit rank-45 minor")

    # Rows 29..41 are residual coefficients 0..12 for the nonidentity section.
    # The pivot retains coefficients 2,3,4 and omits exactly the ten locally
    # forced coefficients.
    pivot_rows = set(map(int, certificate["pivot_row_indices"]))
    selected_residual_coefficients = sorted(
        row - 29 for row in pivot_rows if 29 <= row <= 41
    )
    omitted_residual_coefficients = sorted(
        set(range(13)) - set(selected_residual_coefficients)
    )
    if selected_residual_coefficients != [2, 3, 4]:
        raise ArithmeticError("the stored pivot section equations changed")
    if omitted_residual_coefficients != [0, 1, 5, 6, 7, 8, 9, 10, 11, 12]:
        raise ArithmeticError("unexpected omitted section residuals")

    seed = lift["seed"]["coordinates_mod_7"]
    a = seed[:9]
    support_a_values = {
        "zero": a[0] % 7,
        "one": sum(a) % 7,
        "infinity_scaled": a[8] % 7,
    }
    if not all(support_a_values.values()):
        raise ArithmeticError("A is not a unit at every marked support")

    identity_checks = {
        "over_Q": universal_identity(QQ),
        "over_GF7": universal_identity(GF(7)),
    }
    if not all(identity_checks.values()):
        raise ArithmeticError("universal node/discriminant identity failed")
    jet_checks = {
        "over_Q_depth_1": localized_jet_membership(QQ, 1),
        "over_Q_depth_3": localized_jet_membership(QQ, 3),
        "over_GF7_depth_1": localized_jet_membership(GF(7), 1),
        "over_GF7_depth_3": localized_jet_membership(GF(7), 3),
    }

    output = {
        "schema": "elkies-k3.golay-det720-3a5-formal-smoothness.v1",
        "status": "PASS_ONE_PARAMETER_FORMAL_Z7_MARKED_FAMILY",
        "prime": 7,
        "inputs": {relative(lift_path): digest(lift_path)},
        "support_units_mod_7": support_a_values,
        "universal_identity": {
            "formula": (
                "8*A^3*(X^3+A*X+B) = D*(N-B)-9*B*N^2+N^3, "
                "D=4*A^3+27*B^2, N=2*A*X+3*B"
            ),
            "symbolically_checked": identity_checks,
        },
        "localized_jet_certificates": jet_checks,
        "global_residual_reduction": {
            "finite_support_orders": {"zero": 2, "one": 2},
            "infinity_scaled_order": 6,
            "consequence": (
                "R(t) is divisible by t^2*(t-1)^2 and has degree at most 6, "
                "so R=t^2*(t-1)^2*(c0+c1*t+c2*t^2)."
            ),
            "selected_residual_coefficients": selected_residual_coefficients,
            "quotient_triangular_relations": [
                "R[2]=c0",
                "R[3]=c1-2*c0",
                "R[4]=c2-2*c1+c0",
            ],
            "omitted_residual_coefficients_forced": omitted_residual_coefficients,
        },
        "formal_implicit_function_certificate": {
            "ambient_variables": 46,
            "independent_equations": 45,
            "free_parameter": certificate["omitted_free_variable_names"][0],
            "jacobian_minor_mod_7": int(
                certificate["pivot_minor_determinant_mod_7"]
            ),
            "formal_relative_dimension": 1,
            "conclusion": (
                "The 45 pivot equations define a unique solution in Z_7[[s6-s6_0]] "
                "for the other 45 coordinates. The residual reduction makes all ten "
                "omitted equations automatic, hence the full marked system has the "
                "same one-parameter formal germ."
            ),
        },
        "proof_boundary": {
            "proved": (
                "The complete normalized marked system has a one-parameter formal "
                "Z_7 solution through the pinned GF(7) point. The ten nonpivot section "
                "equations follow from an exact localized Weierstrass identity and the "
                "three retained residual coefficients."
            ),
            "not_proved": (
                "Formal Z_7 existence is not algebraization over Q. No Q-rational "
                "parameterization, physical elliptic-neighbour corridor, effective "
                "target multisection, or specialization rank jump is proved."
            ),
        },
        "reproduce": (
            "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
            "elkies-k3/scripts/certify_golay_det720_3a5_formal_smoothness.sage"
        ),
    }
    serialized = json.dumps(output, indent=2, sort_keys=True) + "\n"
    if arguments.check:
        if output_path.read_text() != serialized:
            raise SystemExit("Golay-720 formal-smoothness artifact is stale")
    else:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(serialized)
    print(
        "GOLAY7203A5FORMAL|independent_equations=45|variables=46|"
        "formal_dimension=1|omitted_equations_forced=10|status=PASS",
        flush=True,
    )


if __name__ == "__main__":
    main()

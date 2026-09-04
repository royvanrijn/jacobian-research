#!/usr/bin/env sage-python
"""Prove formal smoothness of the marked NS0031 model-157 branch.

The normalized presentation has 52 variables and 59 displayed equations.
Its unit Jacobian minor retains 51 rows: every fibre/component row, residual
coefficients 2..8 for the pole-zero section, and residual coefficients 0..16
for the pole-one section.  The exact discriminant/node identity forces the
eight omitted section rows.  Hence the full marked germ is a one-parameter
formally smooth Z_7 branch.

This is a formal local theorem.  It does not algebraize the branch, produce a
Q-rational point, or prove a rational NS0031 marking.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import GF, PolynomialRing, ZZ


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_HENSEL = (
    ROOT
    / "artifacts/generated-results/"
    "elkies-k3-lattice-foundry-ns0031-marked-gf7-hensel-v1.json"
)
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/"
    "elkies-k3-lattice-foundry-ns0031-marked-formal-smoothness-v1.json"
)


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hensel", type=Path, default=DEFAULT_HENSEL)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()
    hensel_path = arguments.hensel.resolve()
    output_path = arguments.output.resolve()
    hensel = json.loads(hensel_path.read_text())

    if hensel.get("schema") != (
        "elkies-k3.lattice-foundry-ns0031-marked-gf7-hensel.v1"
    ):
        raise ValueError("unexpected NS0031 marked-lift schema")
    if hensel.get("status") != (
        "PASS_ONE_DIMENSIONAL_MARKED_TANGENT_AND_EXPLICIT_Z7_"
        "LIFT_TO_REQUESTED_PRECISION"
    ):
        raise ValueError("formal certificate requires the smooth GF(7) seed")
    if hensel["system"]["equation_block_sizes"] != {
        "normalization": 1,
        "fibre_at_zero": 2,
        "fibre_at_one": 8,
        "fibre_at_infinity": 8,
        "component_marking": 8,
        "pole_zero_section": 13,
        "pole_one_section": 19,
    }:
        raise ValueError("marked system block structure changed")

    # Rows 27..39 are pole-zero residual coefficients 0..12; rows 40..58
    # are pole-one residual coefficients 0..18.
    expected_pivot_rows = (
        list(range(27)) + list(range(29, 36)) + list(range(40, 57))
    )
    certificate = hensel["jacobian_certificate"]
    pivot_rows = list(map(int, certificate["pivot_row_indices"]))
    minor = int(certificate["pivot_minor_determinant_mod_7"])
    if pivot_rows != expected_pivot_rows or minor % 7 != 1:
        raise ArithmeticError("unit minor no longer retains the expected rows")
    if certificate["omitted_free_variable_names"] != ["m9"]:
        raise ArithmeticError("the formal free coordinate changed")

    coefficient_ring = PolynomialRing(ZZ, names=("A", "B", "X", "C", "H"))
    A, B, X, C, H = coefficient_ring.gens()
    D = 4 * A**3 + 27 * B**2
    F = X**3 + A * X * C**4 + B * C**6
    identity_right = D * C**4 * (H - B * C**2) - 9 * B * H**2 * C**2 + H**3
    identity_left = 8 * A**3 * F
    if identity_left != identity_right.subs({H: 2 * A * X + 3 * B * C**2}):
        raise ArithmeticError("discriminant/node identity failed")

    field = GF(7)
    coordinates = hensel["seed"]["coordinates_mod_7"]
    names = hensel["seed"]["coordinate_names"]
    seed = dict(zip(names, coordinates))
    a = [field(seed[f"a{index}"]) for index in range(9)]
    support_a_values = {
        "zero": int(a[0]),
        "one": int(sum(a, field.zero())),
        "infinity_scaled": int(a[8]),
    }
    if not all(support_a_values.values()):
        raise ArithmeticError("A is not a unit at every marked support")
    c_at_one = int(field(seed["c0"]) + 1)
    if not c_at_one:
        raise ArithmeticError("the pole denominator meets the I8 support")

    # Exact coefficient bookkeeping for the two global residual arguments.
    p_coefficients = PolynomialRing(ZZ, names=[f"u{i}" for i in range(7)])
    p_t_ring = PolynomialRing(p_coefficients, "t")
    p_t = p_t_ring.gen()
    p_quotient = sum(p_coefficients.gen(i) * p_t**i for i in range(7))
    p_residual = p_t**2 * p_quotient
    if p_residual.degree() != 8:
        raise ArithmeticError("pole-zero residual support calculation changed")
    if [p_residual[i] for i in range(2, 9)] != list(p_coefficients.gens()):
        raise ArithmeticError("pole-zero quotient coefficients are not triangular")

    r_coefficients = PolynomialRing(ZZ, names=("v17", "v18"))
    v17, v18 = r_coefficients.gens()
    r_t_ring = PolynomialRing(r_coefficients, "t")
    r_t = r_t_ring.gen()
    r_tail = v17 * r_t**17 + v18 * r_t**18
    r_value_relations = [r_tail(1), r_tail.derivative()(1)]
    r_relation_matrix = [[1, 1], [17, 18]]
    if r_value_relations != [v17 + v18, 17 * v17 + 18 * v18]:
        raise ArithmeticError("pole-one tail relations changed")
    if ZZ(r_relation_matrix[0][0] * r_relation_matrix[1][1]
          - r_relation_matrix[0][1] * r_relation_matrix[1][0]) != 1:
        raise ArithmeticError("pole-one tail relation matrix is not unimodular")

    payload = {
        "schema": "elkies-k3.ns0031-marked-formal-smoothness.v1",
        "status": "PASS_ONE_DIMENSIONAL_FORMALLY_SMOOTH_Z7_MARKED_FAMILY",
        "prime": 7,
        "inputs": {relative(hensel_path): digest(hensel_path)},
        "identity": {
            "definitions": {
                "D": "4*A^3+27*B^2",
                "H": "2*A*X+3*B*C^2",
                "F": "X^3+A*X*C^4+B*C^6",
            },
            "formula": (
                "8*A^3*F = D*C^4*(H-B*C^2) - "
                "9*B*H^2*C^2 + H^3"
            ),
            "verified_symbolically_over_Z": True,
        },
        "support_localization": {
            "A_at_supports_mod_7": support_a_values,
            "all_A_values_are_units": True,
            "pole_denominator_at_one_mod_7": c_at_one,
            "pole_denominator_is_a_unit_at_one": True,
        },
        "order_arguments": {
            "pole_zero_at_zero": {
                "inputs": [
                    "ord(D)>=2",
                    "ord(H)>=1",
                    "ord(Y)>=1",
                    "A(0) unit",
                ],
                "conclusion": "ord(Y^2-F)>=2",
            },
            "pole_zero_at_infinity": {
                "inputs": [
                    "ord_infinity(D)>=8 in weight 24",
                    "ord_infinity(H)>=2 in weight 12",
                    "ord_infinity(Y)>=2 in weight 6",
                    "A(infinity) unit",
                ],
                "conclusion": "ord_infinity(Y^2-F)>=4 in weight 12",
            },
            "pole_zero_global": (
                "The degree-at-most-12 residual is divisible by t^2 and has "
                "degree at most 8. Its coefficients 2..8 kill the seven "
                "quotient coefficients triangularly."
            ),
            "pole_one_at_one": {
                "inputs": [
                    "ord(D)>=8",
                    "ord(H)>=1",
                    "ord(M)>=1",
                    "A(1) and C(1) units",
                ],
                "conclusion": "ord(M^2-F)>=2",
            },
            "pole_one_global": (
                "After coefficients 0..16 vanish, the degree-at-most-18 "
                "residual is v17*t^17+v18*t^18. Its value and derivative at "
                "t=1 give a unimodular 2-by-2 system, so v17=v18=0."
            ),
        },
        "independent_system": {
            "variable_count": 52,
            "displayed_equation_count": 59,
            "independent_equation_count": 51,
            "independent_displayed_row_indices": expected_pivot_rows,
            "retained_pole_zero_residual_coefficients": list(range(2, 9)),
            "forced_pole_zero_residual_coefficients": [0, 1, 9, 10, 11, 12],
            "retained_pole_one_residual_coefficients": list(range(17)),
            "forced_pole_one_residual_coefficients": [17, 18],
            "unit_minor_mod_7": minor,
            "free_parameter": "m9",
            "formal_relative_dimension": 1,
        },
        "proof_boundary": {
            "proved": (
                "In the localization where A is a unit at all three supports "
                "and C(1) is a unit, the eight omitted section equations "
                "follow from the exact identity, fibre orders, component "
                "jets, and retained residual coefficients. The retained 51 "
                "equations have a unit Jacobian minor, so the full marked "
                "branch is formally smooth of relative dimension one over Z_7."
            ),
            "not_proved": (
                "Formal smoothness does not algebraize or rationally "
                "parameterize the branch over Q, produce a Q-rational point, "
                "prove geometric Picard rank 19, or descend a rational NS0031 "
                "marking."
            ),
        },
        "reproduce": (
            "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
            "elkies-k3/scripts/"
            "certify_lattice_foundry_ns0031_marked_formal_smoothness.sage"
        ),
    }
    serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if arguments.check:
        if not output_path.exists() or output_path.read_text() != serialized:
            raise SystemExit(f"stale artifact: {output_path}")
    else:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(serialized)
    print(
        "FOUNDRYNS0031FORMAL|independent=51|variables=52|minor=1|"
        "forced=8|dimension=1|status=PASS",
        flush=True,
    )


if __name__ == "__main__":
    main()

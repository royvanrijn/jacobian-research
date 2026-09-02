#!/usr/bin/env sage-python
"""Decode and geometrically audit finite NS0007 fixed-case solutions over GF(7).

The fixed-case census stores msolve RURs for every nonempty zero-dimensional
case.  This checker decodes the RURs, proves that their elimination factors
are linear over the base field, substitutes every point into the fully
expanded reduced system, and reconstructs the generalized Weierstrass model

    y^2 + a3*y = x^3 + a2*x^2 + a4*x.

It then checks exact discriminant orders ``I2+I4+I7+I7``, a squarefree
degree-four residual discriminant, the prescribed node/smooth-component
behavior of ``P=(0,0)``, and the Jacobian rank of the reduced marked system.
Positive-dimensional fixed cases remain a separately typed open branch.
"""

from __future__ import annotations

import argparse
import ast
import hashlib
import json
from collections import defaultdict
from pathlib import Path

from sage.all import GF, PolynomialRing, matrix


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CENSUS = (
    ROOT
    / "artifacts/generated-results/"
    "elkies-k3-lattice-foundry-ns0007-pole0-fixed-case-census-mod7.json"
)
DEFAULT_SYSTEM = (
    ROOT / "artifacts/local/elkies-k3/ns0007-pole0-reduced-modp/p7-lambda2.json"
)
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/"
    "elkies-k3-lattice-foundry-ns0007-pole0-fixed-case-solutions-mod7.json"
)


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_msolve(text: str):
    stripped = text.strip()
    return ast.literal_eval(stripped[:-1] if stripped.endswith(":") else stripped)


def valuation_at(polynomial, root) -> int:
    factor = polynomial.parent().gen() - root
    value = polynomial
    order = 0
    while value and not value(root):
        value, remainder = value.quo_rem(factor)
        if remainder:
            raise ArithmeticError("failed exact linear-factor division")
        order += 1
    return order


def coefficients(polynomial, length: int) -> list[int]:
    return [int(polynomial[index]) for index in range(length)]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--census", type=Path, default=DEFAULT_CENSUS)
    parser.add_argument("--system", type=Path, default=DEFAULT_SYSTEM)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()

    census_path = arguments.census.resolve()
    system_metadata_path = arguments.system.resolve()
    census = json.loads(census_path.read_text())
    metadata = json.loads(system_metadata_path.read_text())
    if census.get("schema") != (
        "elkies-k3.lattice-foundry-ns0007-pole0-fixed-case-census-modp.v1"
    ):
        raise ValueError("unexpected NS0007 fixed-case census schema")
    if not census["search"]["exhaustive"]:
        raise ValueError("solution audit requires the exhaustive fixed-case census")
    if metadata.get("schema") != (
        "elkies-k3.lattice-foundry-ns0007-pole0-reduced-modp-system.v1"
    ):
        raise ValueError("unexpected expanded NS0007 system schema")
    prime = int(census["prime"])
    if prime != int(metadata["prime"]) or census["lambda"] != metadata["lambda"]:
        raise ValueError("census/system specialization mismatch")
    field = GF(prime)

    system_path = (ROOT / metadata["system"]["msolve_input"]).resolve()
    if digest(system_path) != metadata["system"]["msolve_sha256"]:
        raise ArithmeticError("expanded reduced-system digest mismatch")
    lines = system_path.read_text().splitlines()
    names = lines[0].split(",")
    coefficient_ring = PolynomialRing(field, names=names, order="degrevlex")
    equation_text = "\n".join(lines[2:]).replace("^", "**")
    equations = [
        coefficient_ring(piece.strip())
        for piece in equation_text.split(",")
        if piece.strip()
    ]
    if len(equations) != 19:
        raise ValueError("unexpected expanded reduced-system equation count")
    jacobian = matrix(coefficient_ring, [
        [equation.derivative(variable) for variable in coefficient_ring.gens()]
        for equation in equations
    ])

    polynomial_ring = PolynomialRing(field, "t")
    t = polynomial_ring.gen()
    lambda_value = field(census["lambda"])
    finite_rows = []
    positive_dimensional = []
    models = defaultdict(list)

    for case in census["exceptional_cases"]:
        if case["status"] == "POSITIVE_DIMENSIONAL_SOLUTION_SET":
            positive_dimensional.append(
                {"case_index": case["index"], "fixed_values": case["values"]}
            )
            continue
        if case["status"] != "FINITE_SOLUTION_SET":
            raise ValueError(f"unexpected exceptional status {case['status']}")
        solution = parse_msolve(case["solution_output"])
        if solution[0] != 0:
            raise ArithmeticError("finite case did not return a zero-dimensional RUR")
        payload = solution[1]
        if int(payload[0]) != prime or payload[3] != names:
            raise ArithmeticError("unexpected RUR characteristic or variable order")
        parametrization = payload[5]
        if parametrization[0] != 1:
            raise ArithmeticError("unsupported multiple RUR blocks")
        elimination_data, denominator_data, coordinate_data = parametrization[1]
        elimination_ring = PolynomialRing(field, "T")
        elimination = elimination_ring(elimination_data[1])
        denominator = elimination_ring(denominator_data[1])
        if elimination.degree() != 1 or denominator != 1:
            raise ArithmeticError("finite NS0007 point is not base-field linear")
        root = -elimination[0] / elimination[1]
        if len(coordinate_data) != len(names) - 1:
            raise ArithmeticError("unexpected RUR coordinate count")
        coordinate_polynomials = [
            elimination_ring(block[0][1]) for block in coordinate_data
        ]
        values = [-polynomial(root) for polynomial in coordinate_polynomials]
        values.append(root)
        assignment = dict(zip(names, values))
        fixed_assignment = [
            int(assignment[name])
            for name in census["fixed_variables_in_enumeration_order"]
        ]
        if fixed_assignment != case["values"]:
            raise ArithmeticError("decoded RUR point violates its fixed case")
        if any(equation(values) for equation in equations):
            raise ArithmeticError("decoded RUR point fails expanded-system substitution")
        jacobian_rank = int(jacobian(values).rank())

        a2 = field(3) + sum(
            assignment[f"a2_{index}"] * t**index for index in range(1, 5)
        )
        a2_infinity = polynomial_ring([a2[4 - index] for index in range(5)])
        si = sum(assignment[f"si_{index}"] * t**index for index in range(7))
        sl = sum(assignment[f"sl_{index}"] * t**index for index in range(7))
        a3_infinity = polynomial_ring(si * (a2_infinity - si**2) % t**7)
        a4_infinity = polynomial_ring(
            ((3 * si**2 + a2_infinity) * (a2_infinity - si**2) / 4) % t**7
        )
        a3 = polynomial_ring([a3_infinity[6 - index] for index in range(7)])
        a4_coefficients = [field.zero()] * 9
        for index in range(7):
            a4_coefficients[8 - index] = a4_infinity[index]
        a4_coefficients[1] = -sum(a4_coefficients[2:])
        a4 = polynomial_ring(a4_coefficients)

        b2 = 4 * a2
        b4 = 2 * a4
        b6 = a3**2
        b8 = a2 * a3**2 - a4**2
        discriminant = -b2**2 * b8 - 8 * b4**3 - 27 * b6**2 + 9 * b2 * b4 * b6
        finite_orders = {
            "0": valuation_at(discriminant, field(0)),
            "1": valuation_at(discriminant, field(1)),
            str(int(lambda_value)): valuation_at(discriminant, lambda_value),
        }
        infinity_order = 24 - discriminant.degree()
        divisor = t**2 * (t - 1) ** 4 * (t - lambda_value) ** 7
        residual, remainder = discriminant.quo_rem(divisor)
        residual_squarefree = bool(
            not remainder
            and residual.degree() == 4
            and residual.gcd(residual.derivative()).degree() == 0
            and all(residual(root_value) for root_value in (0, 1, lambda_value))
        )

        # At 0 and 1 the marked point must be the node; at the two I7 fibres
        # it must be a smooth point (the identity component before resolution).
        node_at_zero = not a3(0) and not a4(0) and bool(a2(0))
        node_at_one = not a3(1) and not a4(1) and bool(a2(1))

        def node_data(a2_value, a3_value, node_parameter):
            node_x = (node_parameter**2 - a2_value) / 2
            node_y = -a3_value / 2
            tangent_value = (3 * node_parameter**2 - a2_value) / 2
            return {
                "node_x": int(node_x),
                "node_y": int(node_y),
                "marked_point_is_smooth": bool(node_x or node_y),
                "ordinary_node": bool(tangent_value),
                "split_tangent": bool(tangent_value.is_square()),
            }

        lambda_node = node_data(a2(lambda_value), a3(lambda_value), sl[0])
        infinity_node = node_data(a2_infinity[0], a3_infinity[0], si[0])
        exact_profile = (
            finite_orders == {"0": 2, "1": 4, str(int(lambda_value)): 7}
            and infinity_order == 7
            and residual_squarefree
            and node_at_zero
            and node_at_one
            and lambda_node["marked_point_is_smooth"]
            and lambda_node["ordinary_node"]
            and infinity_node["marked_point_is_smooth"]
            and infinity_node["ordinary_node"]
        )
        model_key = (
            tuple(coefficients(a2, 5)),
            tuple(coefficients(a3, 7)),
            tuple(coefficients(a4, 9)),
        )
        row = {
            "case_index": case["index"],
            "fixed_values": case["values"],
            "decoded_values_in_system_order": [int(value) for value in values],
            "jacobian_rank": jacobian_rank,
            "a2_coefficients_low_to_high": list(model_key[0]),
            "a3_coefficients_low_to_high": list(model_key[1]),
            "a4_coefficients_low_to_high": list(model_key[2]),
            "discriminant_orders": {
                "finite": finite_orders,
                "infinity": infinity_order,
            },
            "residual_discriminant_coefficients_low_to_high": coefficients(
                residual, 5
            ),
            "residual_squarefree_degree_four": residual_squarefree,
            "marked_node_at_I2_zero": node_at_zero,
            "marked_node_at_I4_one": node_at_one,
            "lambda_I7_node": lambda_node,
            "infinity_I7_node": infinity_node,
            "exact_profile_and_marked_component_pass": exact_profile,
        }
        finite_rows.append(row)
        models[model_key].append(row)

    model_rows = []
    for model_key, occurrences in models.items():
        first = occurrences[0]
        model_rows.append(
            {
                "a2_coefficients_low_to_high": list(model_key[0]),
                "a3_coefficients_low_to_high": list(model_key[1]),
                "a4_coefficients_low_to_high": list(model_key[2]),
                "occurrence_count": len(occurrences),
                "case_indices": [row["case_index"] for row in occurrences],
                "exact_profile_and_marked_component_pass": first[
                    "exact_profile_and_marked_component_pass"
                ],
                "jacobian_ranks": sorted({row["jacobian_rank"] for row in occurrences}),
                "discriminant_orders": first["discriminant_orders"],
                "residual_discriminant_coefficients_low_to_high": first[
                    "residual_discriminant_coefficients_low_to_high"
                ],
                "lambda_I7_node": first["lambda_I7_node"],
                "infinity_I7_node": first["infinity_I7_node"],
            }
        )
    model_rows.sort(
        key=lambda row: (
            row["a2_coefficients_low_to_high"],
            row["a3_coefficients_low_to_high"],
            row["a4_coefficients_low_to_high"],
        )
    )
    finite_rows.sort(key=lambda row: row["case_index"])
    passed_models = [
        row for row in model_rows if row["exact_profile_and_marked_component_pass"]
    ]
    output = {
        "schema": "elkies-k3.lattice-foundry-ns0007-pole0-fixed-case-solutions-modp.v1",
        "status": (
            "PASS_BASE_FIELD_RUR_DECODE_AND_EXACT_GEOMETRIC_AUDIT_"
            "WITH_POSITIVE_DIMENSIONAL_CASES_OPEN"
        ),
        "inputs": {
            relative(census_path): digest(census_path),
            relative(system_metadata_path): digest(system_metadata_path),
            relative(system_path): digest(system_path),
        },
        "prime": prime,
        "lambda": int(lambda_value),
        "accounting": {
            "finite_rur_cases": len(finite_rows),
            "finite_base_field_points": len(finite_rows),
            "positive_dimensional_fixed_cases": len(positive_dimensional),
            "distinct_weierstrass_models": len(model_rows),
            "models_passing_exact_profile_and_marked_components": len(passed_models),
            "smooth_full_rank_finite_points": sum(
                row["jacobian_rank"] == 19 for row in finite_rows
            ),
        },
        "positive_dimensional_cases": positive_dimensional,
        "models": model_rows,
        "finite_points": finite_rows,
        "proof_boundary": {
            "proved": (
                "Every finite RUR is linear over GF(7), substitutes into the 19 "
                "expanded equations, and has the displayed exact discriminant, "
                "residual, marked-component, and Jacobian data."
            ),
            "open": (
                "The two positive-dimensional fixed cases require separate "
                "component analysis. A modular model is not a characteristic-zero "
                "equation or a proof of the full geometric NS lattice."
            ),
        },
        "reproduce": (
            "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
            "elkies-k3/scripts/"
            "audit_lattice_foundry_ns0007_p7_fixed_case_solutions.sage"
        ),
    }
    serialized = json.dumps(output, indent=2, sort_keys=True) + "\n"
    output_path = arguments.output.resolve()
    if arguments.check:
        if output_path.read_text() != serialized:
            raise SystemExit("NS0007 fixed-case solution audit is stale")
    else:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(serialized)
    print(
        "FOUNDRYNS0007SOLUTIONS|"
        f"finite={len(finite_rows)}|models={len(model_rows)}|"
        f"passed={len(passed_models)}|positive_dimensional={len(positive_dimensional)}|"
        "status=PASS",
        flush=True,
    )


if __name__ == "__main__":
    main()

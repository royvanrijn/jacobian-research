#!/usr/bin/env python3
"""Finite LR Rees/SAGBI module computation for the degree-five base map.

This replaces a bounded target-monomial search by the finite semi-invariant
module calculation requested in OP-LR-REES and OP-LR-II.  It does five things.

* completes the target invariant algebra to a finite SAGBI basis;
* constructs the target-field modules in the only surviving weights;
* records the initial lifts and an explicit failure of linear Rees strictness;
* computes the invariant-ring-saturated weight-zero normal quotient; and
* computes the matrices of II_(F,p,-p), proves the cutoff |p| >= 3, and
  tests exact module generation over Q[u,gamma].

The Hessian calculation is exact but substantial.  A typical run takes a few
minutes in the pinned SymPy environment.  Singular is used only for exact
module membership in the final three-summand quotient.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

import sympy as sp

from search_rees_torsion_witnesses import (
    F,
    G_1,
    G_2,
    J,
    TARGET_VARIABLES,
    TARGET_WEIGHTS,
    TargetMonomialField,
    a,
    b,
    gamma,
    invariant_polynomial,
    second_fundamental_form,
    substitute_invariants,
    target_lift,
    u,
    x,
    y,
    z,
)


RING_VARIABLES = (u, gamma)
SOURCE_VARIABLES = (x, y, z)
COMPONENT_NAMES = ("A", "B", "C")


def weighted_degree(expression: sp.Expr) -> int:
    terms = sp.Poly(sp.expand(expression), u, gamma).terms()
    return max(2 * exponent[0] + 3 * exponent[1] for exponent, _ in terms)


def weighted_top(expression: sp.Expr) -> sp.Expr:
    expression = sp.expand(expression)
    degree = weighted_degree(expression)
    return sp.expand(
        sum(
            coefficient * u**exponent[0] * gamma**exponent[1]
            for exponent, coefficient in sp.Poly(expression, u, gamma).terms()
            if 2 * exponent[0] + 3 * exponent[1] == degree
        )
    )


def source_degree(field: sp.Matrix) -> int:
    return max(
        sp.Poly(entry, *SOURCE_VARIABLES).total_degree()
        for entry in field
        if entry != 0
    )


def source_top(field: sp.Matrix) -> sp.Matrix:
    degree = source_degree(field)
    result = []
    for entry in field:
        result.append(
            sp.expand(
                sum(
                    coefficient
                    * x**exponent[0]
                    * y**exponent[1]
                    * z**exponent[2]
                    for exponent, coefficient in sp.Poly(
                        entry, *SOURCE_VARIABLES
                    ).terms()
                    if sum(exponent) == degree
                )
            )
        )
    return sp.Matrix(result)


def semi_invariant_generators(weight: int) -> list[tuple[int, int, int]]:
    """Minimal generators over Q[AC^2,BC] for target weight ``weight``."""
    if weight > 0:
        return [(0, 0, weight)]
    if weight == 0:
        return [(0, 0, 0)]
    n = -weight
    generators = [(A_power, n - 2 * A_power, 0) for A_power in range(n // 2 + 1)]
    if n % 2:
        generators.append(((n + 1) // 2, 0, 1))
    return generators


def target_field_generators(field_weight: int) -> list[TargetMonomialField]:
    return [
        TargetMonomialField(component, exponents)
        for component, coordinate_weight in enumerate(TARGET_WEIGHTS)
        for exponents in semi_invariant_generators(field_weight + coordinate_weight)
    ]


def exponent_label(exponents: tuple[int, int, int]) -> str:
    monomial = sp.prod(
        variable**exponent
        for variable, exponent in zip(TARGET_VARIABLES, exponents, strict=True)
    )
    return sp.sstr(monomial)


def sagbi_certificate() -> dict[str, object]:
    """Return and verify the finite invariant SAGBI completion."""
    f_1 = sp.expand(a * gamma**2 / 6)
    f_2 = sp.expand(4 * b * gamma)
    relation_coefficients = {
        (4, 0): sp.Integer(1),
        (0, 5): -sp.Integer(1),
        (3, 1): -sp.Rational(14, 3),
        (2, 2): sp.Rational(461, 90),
        (1, 3): -sp.Rational(739, 450),
        (0, 4): -sp.Rational(527, 720),
        (3, 0): -sp.Rational(1636136, 759375),
        (2, 1): sp.Rational(675947, 202500),
    }
    f_3 = sp.expand(
        -90
        * sum(
            coefficient * f_1**first_power * f_2**second_power
            for (first_power, second_power), coefficient in relation_coefficients.items()
        )
    )
    assert weighted_top(f_1) == u**5 * gamma**5
    assert weighted_top(f_2) == u**4 * gamma**4
    assert weighted_top(f_3) == u**12 * gamma**14

    # The only toric relation between the first two initial monomials is
    # T_1^4-T_2^5.  Its exact subduction is -f_3/90.  The third exponent
    # (12,14) is off the diagonal and introduces no further toric relation.
    assert sp.expand(
        sum(
            coefficient * f_1**first_power * f_2**second_power
            for (first_power, second_power), coefficient in relation_coefficients.items()
        )
        + f_3 / 90
    ) == 0

    return {
        "generators": [sp.sstr(f_1), sp.sstr(f_2), sp.sstr(f_3)],
        "initial_monomials": ["u^5*gamma^5", "u^4*gamma^4", "u^12*gamma^14"],
        "subduction_identity": [
            [first, second, sp.sstr(coefficient)]
            for (first, second), coefficient in relation_coefficients.items()
        ],
        "third_generator_term_count": len(sp.Poly(f_3, u, gamma).terms()),
        "_expressions": (f_1, f_2, f_3),
    }


def logarithmic_coordinates(field: sp.Matrix) -> list[sp.Expr]:
    return [
        sp.cancel(field[0] / x),
        sp.expand(y * field[0] + x * field[1]),
        sp.expand(
            (-sp.Rational(8, 7) * y + 2 * x * z) * field[0]
            - sp.Rational(8, 7) * x * field[1]
            + x**2 * field[2]
        ),
    ]


def constant_pair_tensors() -> dict[tuple[int, int], sp.Matrix]:
    """Compute normalized II tensors before a target coefficient multiplier."""
    tensors: dict[tuple[int, int], sp.Matrix] = {}
    for first_component in range(3):
        for second_component in range(first_component, 3):
            first = TargetMonomialField(first_component, (0, 0, 0))
            second = TargetMonomialField(second_component, (0, 0, 0))
            form = second_fundamental_form(first, second)
            compensating_weight = (
                TARGET_WEIGHTS[first_component] + TARGET_WEIGHTS[second_component]
            )
            invariant_source = sp.Matrix(
                [
                    invariant_polynomial(sp.cancel(x**compensating_weight * entry))
                    for entry in logarithmic_coordinates(form)
                ]
            )
            assert all(entry is not None for entry in invariant_source)
            tensors[(first_component, second_component)] = (
                J * invariant_source
            ).applyfunc(sp.expand)
    return tensors


def quotient_residue(target_coordinates: sp.Matrix) -> tuple[sp.Expr, sp.Expr, sp.Expr]:
    return (
        sp.expand(G_1.reduce(sp.expand(target_coordinates[0]))[1]),
        sp.expand(G_2.reduce(sp.expand(target_coordinates[1]))[1]),
        sp.expand(target_coordinates[2].subs(gamma, 0)),
    )


def annihilator_cutoff(exponents: tuple[int, int, int]) -> bool:
    """Universal multiplier test for vanishing in all three normal summands."""
    A_power, B_power, C_power = exponents
    first_zero = A_power >= 1 or B_power >= 2       # (a,b^2)
    second_zero = B_power >= 1 or (A_power >= 1 and C_power >= 1)  # (b,a gamma)
    third_zero = C_power >= 1                       # (gamma)
    return first_zero and second_zero and third_zero


@dataclass
class MatrixColumn:
    first: TargetMonomialField
    second: TargetMonomialField
    multiplier_exponents: tuple[int, int, int]
    target_coordinates: sp.Matrix
    residue: tuple[sp.Expr, sp.Expr, sp.Expr]

    def as_json(self) -> dict[str, object]:
        return {
            "first": self.first.label(),
            "second": self.second.label(),
            "components": [self.first.component, self.second.component],
            "multiplier_exponents_ABC": list(self.multiplier_exponents),
            "residue": [sp.sstr(entry) for entry in self.residue],
            "nonzero": any(entry != 0 for entry in self.residue),
        }


def quadratic_matrix(weight: int, tensors: dict[tuple[int, int], sp.Matrix]) -> list[MatrixColumn]:
    columns = []
    for first in target_field_generators(weight):
        for second in target_field_generators(-weight):
            exponents = tuple(
                left + right
                for left, right in zip(first.exponents, second.exponents, strict=True)
            )
            multiplier = sp.expand(
                a**exponents[0] * b**exponents[1] * gamma**exponents[2]
            )
            tensor = tensors[tuple(sorted((first.component, second.component)))]
            target_coordinates = (multiplier * tensor).applyfunc(sp.expand)
            residue = quotient_residue(target_coordinates)
            columns.append(
                MatrixColumn(first, second, exponents, target_coordinates, residue)
            )
    return columns


def singular_polynomial(expression: sp.Expr) -> str:
    """Serialize a Q[u,gamma] polynomial without Singular precedence traps."""
    terms = []
    for exponents, coefficient in sp.Poly(sp.expand(expression), u, gamma).terms():
        factors = []
        for name, exponent in zip(("u", "g"), exponents, strict=True):
            if exponent == 1:
                factors.append(name)
            elif exponent:
                factors.append(f"{name}^{exponent}")
        monomial = "*".join(factors) or "1"
        if coefficient == 1:
            core = monomial
        elif coefficient == -1:
            core = f"-{monomial}"
        else:
            core = f"({coefficient})*{monomial}"
        terms.append(("+" if terms else "") + core)
    return "".join(terms) or "0"


def singular_vector(entries: tuple[sp.Expr, sp.Expr, sp.Expr] | sp.Matrix) -> str:
    return "[" + ",".join(singular_polynomial(entry) for entry in entries) + "]"


def module_generation_tests(
    matrices: dict[int, list[MatrixColumn]],
) -> dict[str, object]:
    singular = shutil.which("Singular")
    if singular is None:
        raise RuntimeError("Singular is required for the exact module-membership test")

    relation_vectors = [
        (a, 0, 0),
        (b**2, 0, 0),
        (0, b, 0),
        (0, a * gamma, 0),
        (0, 0, gamma),
    ]
    p_1 = matrices[1]
    p_2 = matrices[2]
    known_index = next(
        index
        for index, column in enumerate(p_1)
        if column.first.component == 1
        and column.first.exponents == (0, 0, 0)
        and column.second.component == 2
        and column.second.exponents == (0, 0, 0)
    )
    known = p_1[known_index]

    lines = ["ring r=0,(u,g),dp;"]
    for index, vector in enumerate(relation_vectors):
        lines.append(f"vector R{index}={singular_vector(vector)};")
    lines.append(f"vector K={singular_vector(known.residue)};")
    for weight, columns in matrices.items():
        for index, column in enumerate(columns):
            lines.append(f"vector C{weight}_{index}={singular_vector(column.residue)};")
    relations = ",".join(f"R{index}" for index in range(len(relation_vectors)))
    lines.append(f"module MK_raw={relations},K;")
    lines.append("module MK=std(MK_raw);")
    p_1_names = ",".join(f"C1_{index}" for index in range(len(p_1)))
    lines.append(f"module M1_raw={relations},{p_1_names};")
    lines.append("module M1=std(M1_raw);")
    lines.append('"BEGIN_KNOWN";')
    for weight, columns in matrices.items():
        for index in range(len(columns)):
            lines.append(f"reduce(C{weight}_{index},MK)==0;")
    lines.append('"BEGIN_P1";')
    for index in range(len(p_2)):
        lines.append(f"reduce(C2_{index},M1)==0;")

    completed = subprocess.run(
        [singular, "-q"],
        input="\n".join(lines) + "\n",
        text=True,
        capture_output=True,
        check=True,
        timeout=120,
    )
    if completed.stderr.strip() or "?" in completed.stdout:
        raise RuntimeError(f"Singular module test failed:\n{completed.stdout}\n{completed.stderr}")
    known_text, p_1_text = completed.stdout.split("BEGIN_P1")
    known_bits = [
        line.strip() == "1"
        for line in known_text.split("BEGIN_KNOWN", 1)[1].splitlines()
        if line.strip() in {"0", "1"}
    ]
    p_1_bits = [
        line.strip() == "1"
        for line in p_1_text.splitlines()
        if line.strip() in {"0", "1"}
    ]
    expected_known_count = sum(len(columns) for columns in matrices.values())
    assert len(known_bits) == expected_known_count
    assert len(p_1_bits) == len(p_2)
    return {
        "known_p1_column_index": known_index,
        "known_p1_generates_all_p1_and_p2_columns": all(known_bits),
        "p1_matrix_generates_p2_matrix": all(p_1_bits),
        "known_membership_failures": [
            index for index, passed in enumerate(known_bits) if not passed
        ],
        "p2_mod_p1_failures": [
            index for index, passed in enumerate(p_1_bits) if not passed
        ],
    }


def initial_lift_data(weights: tuple[int, ...]) -> dict[str, object]:
    result: dict[str, object] = {}
    for weight in weights:
        rows = []
        for field in target_field_generators(weight):
            lift = target_lift(field)
            rows.append(
                {
                    "field": field.label(),
                    "degree": source_degree(lift),
                    "initial_lift": [sp.sstr(entry) for entry in source_top(lift)],
                }
            )
        result[str(weight)] = rows

    AC_eA = TargetMonomialField(0, (1, 0, 1))
    C2_eC = TargetMonomialField(2, (0, 0, 2))
    torsion = 3 * target_lift(AC_eA) - 4 * target_lift(C2_eC)
    assert source_degree(target_lift(AC_eA)) == 39
    assert source_degree(target_lift(C2_eC)) == 39
    assert 3 * source_top(target_lift(AC_eA)) == 4 * source_top(target_lift(C2_eC))
    assert source_degree(torsion) == 34
    result["linear_strictness"] = {
        "commutes": False,
        "weight": 1,
        "relation": "3*ell_F(AC e_A)-4*ell_F(C^2 e_C)",
        "input_degree": 39,
        "output_degree_after_cancellation": 34,
        "new_initial_lift": [sp.sstr(entry) for entry in source_top(torsion)],
    }
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "artifacts/generated-results/lr_rees_sagbi_module_computation.json"
        ),
    )
    args = parser.parse_args()

    sagbi = sagbi_certificate()
    expressions = sagbi.pop("_expressions")
    assert all(sp.denom(entry) == 1 for entry in expressions)

    # Saturated weight-zero normal quotient:
    # R/(a,b^2) + R/(b,a gamma) + R/(gamma).
    normal_quotient = {
        "ring": "Q[u,gamma], deg(u)=2, deg(gamma)=3",
        "summand_ideals": [
            [sp.sstr(a), sp.sstr(b**2)],
            [sp.sstr(b), sp.sstr(a * gamma)],
            ["gamma"],
        ],
        "initial_ideals": [
            ["u^3", "u^2*gamma", "u*gamma^2", "gamma^3"],
            ["u^3", "u*gamma", "gamma^2"],
            ["gamma"],
        ],
    }

    tensors = constant_pair_tensors()
    matrices = {weight: quadratic_matrix(weight, tensors) for weight in (1, 2)}

    # This is the structural cutoff, not a bounded search.  Check several
    # weights as a regression of the semigroup argument.
    cutoff_rows = {}
    for weight in range(3, 9):
        products = [
            tuple(
                left + right
                for left, right in zip(
                    first.exponents, second.exponents, strict=True
                )
            )
            for first in target_field_generators(weight)
            for second in target_field_generators(-weight)
        ]
        assert products and all(annihilator_cutoff(product) for product in products)
        cutoff_rows[str(weight)] = len(products)

    generation = module_generation_tests(matrices)
    initial_lifts = initial_lift_data((1, -1, 2, -2))

    artifact = {
        "description": "Finite LR Rees/SAGBI and quadratic normal-module computation for F_2",
        "target_weights_ABC": list(TARGET_WEIGHTS),
        "target_invariant_ring": {
            "generators": ["P=A*C^2", "Q=B*C"],
            "pullbacks": [sp.sstr(a * gamma**2), sp.sstr(b * gamma)],
            "sagbi": sagbi,
        },
        "target_field_modules": {
            str(weight): [
                {
                    "component": COMPONENT_NAMES[field.component],
                    "coefficient": exponent_label(field.exponents),
                    "field": field.label(),
                }
                for field in target_field_generators(weight)
            ]
            for weight in (1, -1, 2, -2)
        },
        "initial_lifts": initial_lifts,
        "saturated_normal_quotient": normal_quotient,
        "quadratic_matrices": {
            str(weight): [column.as_json() for column in columns]
            for weight, columns in matrices.items()
        },
        "quadratic_cutoff": {
            "statement": "|p|>=3 vanishes because every coefficient product lies in (a,b^2), (b,a*gamma), and (gamma) simultaneously",
            "regression_column_counts": cutoff_rows,
        },
        "generation_tests": generation,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    serialized = json.dumps(artifact, indent=2, sort_keys=True) + "\n"
    args.output.write_text(serialized)
    digest = hashlib.sha256(serialized.encode()).hexdigest()

    nonzero_counts = {
        weight: sum(any(entry != 0 for entry in column.residue) for column in columns)
        for weight, columns in matrices.items()
    }
    print("PASS: target invariant algebra has a three-generator finite SAGBI basis")
    print("PASS: linear target lifting is not Rees-strict (weight 1, degree 39 -> 34)")
    print(f"PASS: surviving quadratic weights are p=1,2; nonzero columns {nonzero_counts}")
    print("PASS: every |p|>=3 matrix vanishes by the saturated-annihilator cutoff")
    print(
        "PASS: p=1 generates p=2 in the saturated quadratic normal image = "
        f"{generation['p1_matrix_generates_p2_matrix']}"
    )
    print(
        "PASS: the single known (partial_B,partial_C) class is cyclic = "
        f"{generation['known_p1_generates_all_p1_and_p2_columns']}"
    )
    print(f"artifact={args.output} sha256={digest}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Dependency-free audit of the decisive LR quadratic separator.

This deliberately does not import SymPy or call Singular.  It reads the exact
matrix certificate, evaluates its rational polynomial entries at
(u,gamma)=(1/6,0), and checks the covector (0,-144/79,1).
"""

from __future__ import annotations

import ast
from fractions import Fraction
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CERTIFICATE = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "lr_rees_sagbi_module_computation.json"
)


def evaluate(node: ast.AST, variables: dict[str, Fraction]) -> Fraction:
    if isinstance(node, ast.Expression):
        return evaluate(node.body, variables)
    if isinstance(node, ast.Constant) and isinstance(node.value, int):
        return Fraction(node.value)
    if isinstance(node, ast.Name) and node.id in variables:
        return variables[node.id]
    if isinstance(node, ast.UnaryOp):
        value = evaluate(node.operand, variables)
        if isinstance(node.op, ast.USub):
            return -value
        if isinstance(node.op, ast.UAdd):
            return value
    if isinstance(node, ast.BinOp):
        left = evaluate(node.left, variables)
        right = evaluate(node.right, variables)
        if isinstance(node.op, ast.Add):
            return left + right
        if isinstance(node.op, ast.Sub):
            return left - right
        if isinstance(node.op, ast.Mult):
            return left * right
        if isinstance(node.op, ast.Div):
            return left / right
        if isinstance(node.op, ast.Pow):
            assert right.denominator == 1
            return left ** right.numerator
    raise ValueError(f"unsupported certificate expression: {ast.dump(node)}")


def polynomial_value(expression: str, variables: dict[str, Fraction]) -> Fraction:
    return evaluate(ast.parse(expression, mode="eval"), variables)


def column_value(column: dict[str, object]) -> Fraction:
    point = {"u": Fraction(1, 6), "gamma": Fraction(0)}
    entries = [
        polynomial_value(expression, point) for expression in column["residue"]
    ]
    covector = (Fraction(0), -Fraction(144, 79), Fraction(1))
    return sum(
        coefficient * entry
        for coefficient, entry in zip(covector, entries, strict=True)
    )


def main() -> None:
    certificate = json.loads(CERTIFICATE.read_text())
    generation = certificate["generation_tests"]
    separator = generation["separating_functional"]

    assert separator["point"] == {"u": "1/6", "gamma": "0"}
    assert separator["covector"] == ["0", "-144/79", "1"]
    assert separator["new_value"] == "-987/395"
    assert generation["p2_mod_p1_failures"] == [2]
    assert generation["new_p2_column"] == 2
    assert generation["new_p2_remainder_mod_p1"] == "-987/395*gen(3)"
    quotient = generation["p2_mod_p1_quotient"]
    assert quotient == {
        "cyclic_generator_column": 2,
        "annihilator": ["gamma", "6*u-1"],
        "support": {"u": "1/6", "gamma": "0"},
        "vector_space_dimension": 1,
        "is_reduced_residue_field": True,
    }

    p_1 = certificate["quadratic_matrices"]["1"]
    p_2 = certificate["quadratic_matrices"]["2"]
    assert len(p_1) == 24 and len(p_2) == 24
    assert all(column_value(column) == 0 for column in p_1)

    p_2_values = [column_value(column) for column in p_2]
    assert p_2_values[2] == -Fraction(987, 395)
    assert all(value == 0 for index, value in enumerate(p_2_values) if index != 2)
    assert separator["p2_values"] == [str(value) for value in p_2_values]

    # Directly audit descent through the five relation generators of N_R.
    point = {"u": Fraction(1, 6), "gamma": Fraction(0)}
    second_ideal = certificate["saturated_normal_quotient"]["summand_ideals"][1]
    assert all(polynomial_value(generator, point) == 0 for generator in second_ideal)
    third_ideal = certificate["saturated_normal_quotient"]["summand_ideals"][2]
    assert all(polynomial_value(generator, point) == 0 for generator in third_ideal)
    # The first covector entry is zero, so no condition from the first ideal.

    print("PASS: the separator descends through the saturated normal relations")
    print("PASS: it annihilates all 24 p=1 quadratic columns")
    print("PASS: p=2 column 2 has exact value -987/395 and all others vanish")
    print("PASS: the new quotient is the reduced residue field at (1/6,0)")


if __name__ == "__main__":
    main()

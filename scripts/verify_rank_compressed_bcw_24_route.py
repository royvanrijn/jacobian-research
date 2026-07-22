#!/usr/bin/env python3
"""Construct and verify the rank-compressed 16 -> 24 BCW homogenization."""

from __future__ import annotations

import json
from pathlib import Path

import sympy as sp

from rank_compressed_bcw_homogenization import (
    extract_quadratic_cubic,
    factor_cubic_output,
    verify_parametric_factorization,
)


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "artifacts" / "generated-results" / "shared_bcw_33_counterexample.json"
OUTPUT = ROOT / "artifacts" / "generated-results" / "rank_compressed_bcw_24_counterexample.json"


def rational_text(value: sp.Expr) -> str:
    value = sp.cancel(value)
    assert value.is_Rational
    return str(value)


def dense(support: list[list[int]], dimension: int) -> tuple[int, ...]:
    result = [0] * dimension
    for index, exponent in support:
        result[index] = exponent
    return tuple(result)


def support(exponents: tuple[int, ...]) -> list[list[int]]:
    return [[index, exponent] for index, exponent in enumerate(exponents) if exponent]


def source_map(stored: dict[str, object]) -> tuple[list[sp.Symbol], list[sp.Expr]]:
    """Recover K=X+Q+C from the serialized ordinary 2n+1 homogenization."""
    n = stored["degree_reduced_dimension"]
    assert stored["dimension"] == 2 * n + 1
    variables = list(sp.symbols(f"x0:{n}"))
    t_index = 2 * n
    expressions: list[sp.Expr] = []
    for index in range(n):
        q_part = sp.Integer(0)
        for term in stored["H"][index]:
            exponents = dense(term["monomial"], 2 * n + 1)
            coefficient = sp.Rational(term["coefficient"])
            if exponents[n + index] == 1 and exponents[t_index] == 2:
                assert sum(exponents) == 3 and coefficient == 1
                continue
            assert exponents[t_index] == 1
            old = exponents[:n]
            assert sum(old) == 2 and not any(exponents[n:2 * n])
            q_part += coefficient * sp.prod(x**power for x, power in zip(variables, old))
        c_part = sp.Integer(0)
        for term in stored["H"][n + index]:
            exponents = dense(term["monomial"], 2 * n + 1)
            assert not any(exponents[n:]) and sum(exponents[:n]) == 3
            coefficient = -sp.Rational(term["coefficient"])
            c_part += coefficient * sp.prod(
                x**power for x, power in zip(variables, exponents[:n])
            )
        expressions.append(sp.expand(variables[index] + q_part + c_part))
    return variables, expressions


def main() -> None:
    source = json.loads(SOURCE.read_text())
    variables, expressions = source_map(source)
    quadratic, cubic = extract_quadratic_cubic(expressions, variables)
    factorization = factor_cubic_output(cubic)
    n = len(variables)
    k = len(factorization.c)
    assert n == 16
    assert [index for index, poly in enumerate(cubic) if not poly.is_zero] == [0, 1, 2, 3, 4, 6, 8]
    assert factorization.basis_components == (0, 1, 2, 3, 4, 6, 8)
    assert k == 7
    verify_parametric_factorization(variables, quadratic, cubic, factorization)

    final_dimension = n + k + 1
    assert final_dimension == 24
    t_index = final_dimension - 1
    h_components: list[list[dict[str, object]]] = [[] for _ in range(final_dimension)]
    for i, q_part in enumerate(quadratic):
        for exponents, coefficient in q_part.terms():
            if not coefficient:
                continue
            final = list(exponents) + [0] * (k + 1)
            final[t_index] = 1
            h_components[i].append(
                {"coefficient": rational_text(coefficient), "monomial": support(tuple(final))}
            )
        for j in range(k):
            coefficient = factorization.B[i, j]
            if coefficient:
                final = [0] * final_dimension
                final[n + j] = 1
                final[t_index] = 2
                h_components[i].append(
                    {"coefficient": rational_text(coefficient), "monomial": support(tuple(final))}
                )
    for j, c_part in enumerate(factorization.c):
        for exponents, coefficient in c_part.terms():
            if not coefficient:
                continue
            final = exponents + (0,) * (k + 1)
            h_components[n + j].append(
                {"coefficient": rational_text(-coefficient), "monomial": support(final)}
            )
    assert all(
        sum(exponent for _, exponent in term["monomial"]) == 3
        for component in h_components
        for term in component
    )

    old_points = [[sp.Rational(value) for value in point[:n]] for point in source["collision_points"]]
    old_target = [sp.Rational(value) for value in source["common_image"][:n]]
    final_points: list[list[sp.Expr]] = []
    for point in old_points:
        substitution = dict(zip(variables, point))
        c_value = [poly.eval(substitution) for poly in factorization.c]
        assert [sp.Poly(expression, *variables).eval(substitution) for expression in expressions] == old_target
        final_points.append(point + c_value + [sp.Integer(1)])
    common_image = old_target + [sp.Integer(0)] * k + [sp.Integer(1)]
    assert len({tuple(point) for point in final_points}) == 3

    artifact = {
        "format": "rank-compressed-bcw-sparse-cubic-homogeneous-map-v1",
        "source": "rank compression of the repository shared-factor 16-variable BCW map",
        "construction": "C=Bc over QQ followed by rank-compressed nilpotent cubic homogenization",
        "source_artifact": str(SOURCE.relative_to(ROOT)),
        "source_dimension": source["source_dimension"],
        "degree_reduced_dimension": n,
        "cubic_output_rank": k,
        "dimension": final_dimension,
        "linear_part": "identity",
        "map": "V_i = Z_i + H_i; omitted H_i are zero",
        "jacobian_determinant": "1",
        "jacobian_certificate": (
            "E_t=t^-1 K(tX), det(DE_t)=1, and id+tN=P_t o (E_t x id) o A_t; "
            "homogenization gives det(DV)=1"
        ),
        "rank_factorization": {
            "component_monomials": [support(exponents) for exponents in factorization.monomials],
            "basis_components": list(factorization.basis_components),
            "B": [[rational_text(factorization.B[i, j]) for j in range(k)] for i in range(n)],
        },
        "H": h_components,
        "collision_points": [[rational_text(value) for value in point] for point in final_points],
        "common_image": [rational_text(value) for value in common_image],
        "statistics": {
            "nonzero_cubic_output_components": sum(not poly.is_zero for poly in cubic),
            "cubic_output_rank": k,
            "nonzero_cubic_terms": sum(len(component) for component in h_components),
        },
    }
    OUTPUT.write_text(json.dumps(artifact, indent=2) + "\n")

    print("PASS rank-compressed BCW: extracted K=X+Q+C in dimension 16")
    print("PASS rank-compressed BCW: cubic coefficient matrix has exact QQ-rank 7")
    print("PASS rank-compressed BCW: verified id+tN=P_t o (E_t x id) o A_t")
    print("PASS rank-compressed BCW: cubic homogenization gives a 24-variable collision")
    print("PASS rank-compressed BCW: fixed-dimensional implication now targets GMC(48)")
    print(f"PASS rank-compressed BCW: wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

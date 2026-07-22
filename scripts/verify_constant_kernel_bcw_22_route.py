#!/usr/bin/env python3
"""Remove the constant-Jacobian-kernel directions from the 24D BCW map.

If F=X+H and K is the constant kernel of JH, then H is constant on cosets
of K.  For a quotient B with ker(B)=K and a section C with BC=I, the map

    f(q) = q + B H(Cq)

is a cubic-homogeneous Keller map.  In coordinates X=Cq+Kr the original
map is the triangular extension (f(q), r+g(q)), so determinants agree and
every collision descends under B.
"""

from __future__ import annotations

import json
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "artifacts" / "generated-results" / "rank_compressed_bcw_24_counterexample.json"
OUTPUT = ROOT / "artifacts" / "generated-results" / "constant_kernel_bcw_22_counterexample.json"


def qtext(value: sp.Expr) -> str:
    value = sp.cancel(value)
    assert value.is_Rational
    return str(value)


def sparse_vector(vector: sp.Matrix) -> list[list[object]]:
    return [[i, qtext(value)] for i, value in enumerate(vector) if value]


def sparse_matrix_rows(matrix: sp.Matrix) -> list[list[list[object]]]:
    return [sparse_vector(matrix.row(i).T) for i in range(matrix.rows)]


def monomial(variables: list[sp.Symbol], support: list[list[int]]) -> sp.Expr:
    return sp.prod(variables[index] ** exponent for index, exponent in support)


def decode_h(stored: dict[str, object]) -> tuple[list[sp.Symbol], list[sp.Expr]]:
    variables = list(sp.symbols(f"x0:{stored['dimension']}"))
    components = [
        sp.expand(
            sum(
                sp.Rational(term["coefficient"])
                * monomial(variables, term["monomial"])
                for term in component
            )
        )
        for component in stored["H"]
    ]
    return variables, components


def encode_h(components: list[sp.Expr], variables: list[sp.Symbol]) -> list[list[dict[str, object]]]:
    encoded: list[list[dict[str, object]]] = []
    for expression in components:
        terms: list[dict[str, object]] = []
        for exponents, coefficient in sp.Poly(expression, *variables, domain=sp.QQ).terms():
            if coefficient:
                terms.append(
                    {
                        "coefficient": qtext(coefficient),
                        "monomial": [
                            [index, exponent]
                            for index, exponent in enumerate(exponents)
                            if exponent
                        ],
                    }
                )
        encoded.append(terms)
    return encoded


def constant_kernel(components: list[sp.Expr], variables: list[sp.Symbol]) -> sp.Matrix:
    jacobian = sp.Matrix(components).jacobian(variables)
    monomials = sorted(
        {
            exponents
            for entry in jacobian
            for exponents, coefficient in sp.Poly(entry, *variables, domain=sp.QQ).terms()
            if coefficient
        }
    )
    coefficient_matrix = sp.Matrix(
        [
            [
                sp.Poly(jacobian[i, j], *variables, domain=sp.QQ).coeff_monomial(exponents)
                for j in range(len(variables))
            ]
            for i in range(len(components))
            for exponents in monomials
        ]
    )
    basis = coefficient_matrix.nullspace()
    return sp.Matrix.hstack(*basis) if basis else sp.zeros(len(variables), 0)


def dependency_sccs(components: list[sp.Expr], variables: list[sp.Symbol]) -> list[list[int]]:
    jacobian = sp.Matrix(components).jacobian(variables)
    edges = {
        i: {j for j in range(len(variables)) if jacobian[i, j] != 0}
        for i in range(len(components))
    }
    indices: dict[int, int] = {}
    low: dict[int, int] = {}
    stack: list[int] = []
    on_stack: set[int] = set()
    result: list[list[int]] = []

    def visit(vertex: int) -> None:
        indices[vertex] = low[vertex] = len(indices)
        stack.append(vertex)
        on_stack.add(vertex)
        for target in edges[vertex]:
            if target not in indices:
                visit(target)
                low[vertex] = min(low[vertex], low[target])
            elif target in on_stack:
                low[vertex] = min(low[vertex], indices[target])
        if low[vertex] == indices[vertex]:
            component: list[int] = []
            while True:
                target = stack.pop()
                on_stack.remove(target)
                component.append(target)
                if target == vertex:
                    break
            result.append(sorted(component))

    for vertex in range(len(variables)):
        if vertex not in indices:
            visit(vertex)
    return result


def main() -> None:
    source = json.loads(SOURCE.read_text())
    variables, h = decode_h(source)
    n = len(variables)
    assert n == 24

    kernel = constant_kernel(h, variables)
    assert kernel.shape == (24, 2)
    expected_kernel = sp.zeros(24, 2)
    expected_kernel[6, 0] = -3
    expected_kernel[10, 0] = 1
    expected_kernel[12, 1] = -sp.Rational(2, 11)
    expected_kernel[13, 1] = -sp.Rational(3, 11)
    expected_kernel[15, 1] = 1
    assert kernel.columnspace() == expected_kernel.columnspace()
    kernel = expected_kernel

    # Rows of B are a basis of the annihilator of K.  A pivot-coordinate
    # section gives a sparse rational C with BC=I.
    B = sp.Matrix.hstack(*kernel.T.nullspace()).T
    assert B.shape == (22, 24) and B * kernel == sp.zeros(22, 2)
    pivots = list(B.rref()[1])
    C = sp.zeros(24, 22)
    pivot_inverse = B[:, pivots].inv()
    for local_index, ambient_index in enumerate(pivots):
        C[ambient_index, :] = pivot_inverse[local_index, :]
    assert B * C == sp.eye(22)
    change = C.row_join(kernel)
    assert change.det() != 0
    inverse_change = change.inv()
    assert inverse_change[:22, :] == B

    # Exact GZ-type quotient identity: H(X)=H(CBX), hence f=BFC.
    cb_substitution = dict(zip(variables, list(C * B * sp.Matrix(variables))))
    assert all(sp.expand(left - right.subs(cb_substitution)) == 0 for left, right in zip(h, h))

    quotient_variables = list(sp.symbols("q0:22"))
    section_substitution = dict(zip(variables, list(C * sp.Matrix(quotient_variables))))
    quotient_h = [sp.expand(value) for value in B * sp.Matrix(h).subs(section_substitution)]
    assert all(
        sp.Poly(value, *quotient_variables, domain=sp.QQ).total_degree() in (-sp.oo, 3)
        for value in quotient_h
    )
    assert constant_kernel(quotient_h, quotient_variables).cols == 0

    # In (q,r) coordinates the source map is a skew product with a literal
    # two-dimensional identity linear factor.  This is the determinant bridge.
    fiber_h = [
        sp.expand(value)
        for value in inverse_change[22:, :] * sp.Matrix(h).subs(section_substitution)
    ]
    assert len(fiber_h) == 2

    source_points = [
        sp.Matrix([sp.Rational(value) for value in point])
        for point in source["collision_points"]
    ]
    source_image = sp.Matrix([sp.Rational(value) for value in source["common_image"]])
    quotient_points = [B * point for point in source_points]
    quotient_image = B * source_image
    assert len({tuple(point) for point in quotient_points}) == 3
    for point in quotient_points:
        substitution = dict(zip(quotient_variables, point))
        assert point + sp.Matrix(quotient_h).subs(substitution) == quotient_image

    encoded_h = encode_h(quotient_h, quotient_variables)
    used_variables = sorted(
        {
            index
            for component in encoded_h
            for term in component
            for index, _ in term["monomial"]
        }
    )
    nonzero_outputs = [index for index, component in enumerate(encoded_h) if component]
    sccs = dependency_sccs(quotient_h, quotient_variables)
    assert used_variables == list(range(22))
    assert nonzero_outputs == list(range(21))
    assert sorted(sccs, key=len) == [[21], list(range(21))]

    artifact = {
        "format": "constant-kernel-quotient-bcw-sparse-cubic-homogeneous-map-v1",
        "source": "constant-JH-kernel quotient of the rank-compressed 24-variable BCW map",
        "source_artifact": str(SOURCE.relative_to(ROOT)),
        "dimension": 22,
        "linear_part": "identity",
        "map": "f(q)=q+B H(Cq)",
        "jacobian_determinant": "1",
        "jacobian_certificate": (
            "ker(B)=constant kernel(JH), BC=I, and in X=Cq+Kr coordinates "
            "the 24D map is (f(q),r+g(q))"
        ),
        "quotient_factorization": {
            "kernel_basis_columns": sparse_matrix_rows(kernel.T),
            "B_shape": [22, 24],
            "B_rows": sparse_matrix_rows(B),
            "C_shape": [24, 22],
            "C_rows": sparse_matrix_rows(C),
            "identities": ["B K=0", "B C=I_22", "H=H C B", "f=B F C"],
        },
        "H": encoded_h,
        "collision_points": [
            [qtext(value) for value in point] for point in quotient_points
        ],
        "common_image": [qtext(value) for value in quotient_image],
        "diagnostics": {
            "source_constant_kernel_dimension": 2,
            "quotient_constant_kernel_dimension": 0,
            "variables_in_coordinate_support": used_variables,
            "nonzero_H_outputs": nonzero_outputs,
            "coordinate_dependency_sccs": sccs,
            "triangular_identity_fiber_dimension_removed": 2,
            "remaining_identity_output_coordinates": [21],
        },
        "statistics": {
            "nonzero_cubic_terms": sum(len(component) for component in encoded_h),
            "dimension": 22,
            "fixed_dimensional_GMC_target": 44,
        },
    }
    OUTPUT.write_text(json.dumps(artifact, indent=2) + "\n")

    print("PASS constant-kernel BCW: dim ker_Q(JH)=2 in the 24D map")
    print("PASS constant-kernel BCW: exact B,C quotient and triangular-fiber identities")
    print("PASS constant-kernel BCW: quotient has zero constant JH-kernel")
    print("PASS constant-kernel BCW: cubic collision descends to dimension 22")
    print("PASS constant-kernel BCW: fixed-dimensional implication targets GMC(44)")
    print(f"PASS constant-kernel BCW: wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

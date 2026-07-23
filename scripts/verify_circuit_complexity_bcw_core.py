#!/usr/bin/env python3
"""Verify the output-feedback core behind the circuit-level BCW endpoints.

For a polynomial correction H with constant output-space factorization

    H = B p,

put R = Jp B.  Then JH = B Jp and, for k >= 1,

    (JH)^k = B R^(k-1) Jp.

The script verifies this factorization for the two current foundational
rank/index endpoints, computes the exact polynomial nilpotency index of R,
and audits the feedback graph.  It also records a compact ten-multiplication,
depth-four circuit for the foundational map and demonstrates why the naive
quadratic graph-map compilation is not determinant preserving off the circuit
graph.
"""

from __future__ import annotations

import json
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
ARTIFACTS = ROOT / "artifacts" / "generated-results"
OUTPUT = ARTIFACTS / "circuit_complexity_bcw_core.json"

ENDPOINTS = (
    (
        "rank_endpoint",
        "rank_reduced_bcw_24_counterexample.json",
        17,
        18,
    ),
    (
        "index_endpoint",
        "index_reduced_bcw_22_counterexample.json",
        18,
        18,
    ),
)


def foundational_circuit() -> dict[str, object]:
    """Verify a ten-gate, multiplication-depth-four factorization."""

    x, y, z = sp.symbols("x y z")
    gate_depths = {
        "a": 1,
        "r": 1,
        "x2": 1,
        "p0": 2,
        "v2": 2,
        "v2z": 3,
        "y2": 1,
        "y2b": 2,
        "xw": 4,
        "vw": 4,
    }
    assert len(gate_depths) == 10
    assert max(gate_depths.values()) == 4
    a = x * y
    r = x * z
    x2 = x * x
    p0 = x2 * (3 * y + r)
    v = 1 + 2 * a
    b = 2 + 3 * a
    v2 = v * v
    v2z = v2 * z
    y2 = y * y
    y2b = y2 * b
    w = v2z + 4 * y2b
    xw = x * w
    vw = v * w
    circuit_map = sp.Matrix((x - p0, y + 3 * xw, vw))

    long_map = sp.Matrix(
        (
            v**3 * z + 4 * y**2 * v * b,
            y + 3 * x * v**2 * z + 12 * x * y**2 * b,
            -x + 3 * x**2 * y + x**3 * z,
        )
    )
    normalized = sp.Matrix((-long_map[2], long_map[1], long_map[0]))
    assert all(sp.expand(left - right) == 0 for left, right in zip(circuit_map, normalized))
    assert sp.factor(circuit_map.jacobian((x, y, z)).det()) == 1

    points = (
        (sp.Integer(0), sp.Integer(0), -sp.Rational(1, 8)),
        (sp.Integer(1), -sp.Rational(3, 4), sp.Rational(13, 4)),
        (-sp.Integer(1), sp.Rational(3, 4), sp.Rational(13, 4)),
    )
    images = tuple(
        tuple(circuit_map.subs(dict(zip((x, y, z), point))))
        for point in points
    )
    assert images[0] == images[1] == images[2]

    # The tempting graph construction uses one variable per multiplication
    # gate.  It agrees with F on the circuit graph but is not Keller on the
    # full stabilized space.  Singularity of its affine jet is an exact,
    # inexpensive obstruction.
    u = list(sp.symbols("graph_u0:10"))
    gate_relations = (
        x * y,
        x * z,
        x * x,
        u[2] * (3 * y + u[1]),
        (1 + 2 * u[0]) ** 2,
        u[4] * z,
        y * y,
        u[6] * (2 + 3 * u[0]),
        x * (u[5] + 4 * u[7]),
        (1 + 2 * u[0]) * (u[5] + 4 * u[7]),
    )
    graph_variables = (x, y, z, *u)
    graph_map = sp.Matrix(
        (
            x - u[3],
            y + 3 * u[8],
            u[9],
            *(u[index] - gate_relations[index] for index in range(10)),
        )
    )
    origin = dict.fromkeys(graph_variables, 0)
    affine_jet = graph_map.jacobian(graph_variables).subs(origin)
    assert affine_jet.rank() == 12
    assert affine_jet.det() == 0

    return {
        "multiplication_gates": 10,
        "multiplication_depth": 4,
        "output_roots": 3,
        "factorization": {
            "a": "x*y",
            "r": "x*z",
            "x2": "x*x",
            "p0": "x2*(3*y+r)",
            "v2": "(1+2*a)^2",
            "v2z": "v2*z",
            "y2": "y*y",
            "y2b": "y2*(2+3*a)",
            "w": "v2z+4*y2b",
            "xw": "x*w",
            "vw": "(1+2*a)*w",
            "outputs": ["x-p0", "y+3*xw", "vw"],
        },
        "determinant": "1",
        "collision_size": 3,
        "naive_graph_compilation": {
            "dimension": 13,
            "affine_jet_rank": 12,
            "affine_jet_determinant": "0",
            "conclusion": (
                "agreement on the circuit graph does not preserve the "
                "Jacobian determinant on the ambient stabilization"
            ),
        },
    }


def parse_h(artifact: dict[str, object]) -> tuple[tuple[sp.Symbol, ...], tuple[sp.Poly, ...]]:
    dimension = int(artifact["dimension"])
    variables = tuple(sp.symbols(f"core_x0:{dimension}"))
    components: list[sp.Poly] = []
    for raw_component in artifact["H"]:
        expression = sum(
            sp.Rational(term["coefficient"])
            * sp.prod(variables[index] ** exponent for index, exponent in term["monomial"])
            for term in raw_component
        )
        components.append(sp.Poly(expression, *variables, domain=sp.QQ))
    return variables, tuple(components)


def output_factorization(
    components: tuple[sp.Poly, ...],
) -> tuple[sp.Matrix, tuple[sp.Poly, ...], tuple[int, ...]]:
    variables = components[0].gens
    monomials = tuple(
        sorted(
            {
                monomial
                for component in components
                for monomial, coefficient in component.terms()
                if coefficient
            }
        )
    )
    coefficients = sp.Matrix(
        [
            [component.coeff_monomial(monomial) for monomial in monomials]
            for component in components
        ]
    )
    pivots = tuple(coefficients.T.rref()[1])
    basis_coefficients = coefficients[list(pivots), :]
    rows = []
    for row in range(coefficients.rows):
        solution, parameters = basis_coefficients.T.gauss_jordan_solve(
            coefficients[row, :].T
        )
        assert parameters.rows == 0
        rows.append(list(solution))
    B = sp.Matrix(rows)
    basis = tuple(components[index] for index in pivots)
    assert B.rank() == len(pivots)
    reconstructed = B * sp.Matrix([component.as_expr() for component in basis])
    assert all(
        sp.Poly(
            reconstructed[index] - components[index].as_expr(),
            *variables,
            domain=sp.QQ,
        ).is_zero
        for index in range(len(components))
    )
    return B, basis, pivots


def strongly_connected_component_sizes(adjacency: tuple[tuple[int, ...], ...]) -> list[int]:
    indices = [-1] * len(adjacency)
    lowlinks = [0] * len(adjacency)
    stack: list[int] = []
    on_stack: set[int] = set()
    components: list[list[int]] = []
    counter = 0

    def visit(vertex: int) -> None:
        nonlocal counter
        indices[vertex] = counter
        lowlinks[vertex] = counter
        counter += 1
        stack.append(vertex)
        on_stack.add(vertex)
        for target in adjacency[vertex]:
            if indices[target] < 0:
                visit(target)
                lowlinks[vertex] = min(lowlinks[vertex], lowlinks[target])
            elif target in on_stack:
                lowlinks[vertex] = min(lowlinks[vertex], indices[target])
        if lowlinks[vertex] == indices[vertex]:
            component = []
            while True:
                target = stack.pop()
                on_stack.remove(target)
                component.append(target)
                if target == vertex:
                    break
            components.append(component)

    for vertex in range(len(adjacency)):
        if indices[vertex] < 0:
            visit(vertex)
    return sorted((len(component) for component in components), reverse=True)


def audit_endpoint(
    label: str,
    filename: str,
    expected_rank: int,
    expected_index: int,
) -> dict[str, object]:
    artifact = json.loads((ARTIFACTS / filename).read_text())
    variables, components = parse_h(artifact)
    B, basis, pivots = output_factorization(components)
    Jp = sp.Matrix([component.as_expr() for component in basis]).jacobian(variables)
    JH = sp.Matrix([component.as_expr() for component in components]).jacobian(variables)
    assert JH == B * Jp
    core = Jp * B

    # The graph uses the convention input core coordinate j -> output core
    # coordinate i when R_ij is a nonzero polynomial.
    adjacency = tuple(
        tuple(
            row
            for row in range(core.rows)
            if not sp.Poly(core[row, column], *variables, domain=sp.QQ).is_zero
        )
        for column in range(core.cols)
    )
    scc_sizes = strongly_connected_component_sizes(adjacency)

    power = core
    power_profile = []
    core_index = None
    for exponent in range(1, core.rows + 2):
        nonzero_entries = 0
        scalar_terms = 0
        for entry in power:
            polynomial = sp.Poly(entry, *variables, domain=sp.QQ)
            if not polynomial.is_zero:
                nonzero_entries += 1
                scalar_terms += len(polynomial.terms())
        power_profile.append(
            {
                "power": exponent,
                "nonzero_entries": nonzero_entries,
                "scalar_terms": scalar_terms,
            }
        )
        if nonzero_entries == 0:
            core_index = exponent
            break
        power = (power * core).applyfunc(sp.expand)

    assert core_index == 17
    assert power_profile[-2]["nonzero_entries"] > 0
    statistics = artifact["statistics"]
    assert int(statistics["generic_rank_JH_over_QQ_x"]) == expected_rank
    assert int(statistics["nilpotency_index_JH"]) == expected_index
    assert int(statistics["exact_certificate"]["nonzero_entries_JH_power_17"]) > 0
    assert int(statistics["exact_certificate"]["nonzero_entries_JH_power_18"]) == 0

    # Since H=Bp, JH=BJp and (JH)^k=B R^(k-1) Jp.  R^17=0 therefore
    # supplies the exact upper certificate (JH)^18=0; the stored nonzero
    # seventeenth power supplies the lower certificate.
    return {
        "label": label,
        "source_artifact": filename,
        "ambient_dimension": int(artifact["dimension"]),
        "generic_rank_JH": expected_rank,
        "nilpotency_index_JH": expected_index,
        "output_core_dimension": len(pivots),
        "basis_components": list(pivots),
        "feedback_edges": sum(len(targets) for targets in adjacency),
        "feedback_scc_sizes": scc_sizes,
        "feedback_core_nilpotency_index": core_index,
        "identity": "(JH)^k = B*(Jp*B)^(k-1)*Jp",
        "power_profile": power_profile,
    }


def main() -> None:
    payload = {
        "format": "circuit-complexity-bcw-output-core-v1",
        "theorem": {
            "factorization": "H=Bp with B constant of full column rank s",
            "feedback_core": "R=Jp*B",
            "power_identity": "(JH)^k=B*R^(k-1)*Jp for k>=1",
            "rank_bound": "rank(JH)<=s",
            "index_bound": "nilpotency_index(JH)<=nilpotency_index(R)+1<=s+1",
        },
        "foundational_circuit": foundational_circuit(),
        "endpoint_cores": [
            audit_endpoint(label, filename, rank, index)
            for label, filename, rank, index in ENDPOINTS
        ],
        "conclusion": (
            "nilpotency is controlled by the post-placement output-feedback "
            "core; the raw multiplication DAG alone omits this feedback data"
        ),
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    print("PASS circuit complexity: foundational map has a 10-gate depth-4 circuit")
    print("PASS circuit complexity: naive one-variable-per-gate graph map has singular affine jet")
    for endpoint in payload["endpoint_cores"]:
        print(
            "PASS output core:"
            f" {endpoint['label']} has core dimension"
            f" {endpoint['output_core_dimension']} and exact core index"
            f" {endpoint['feedback_core_nilpotency_index']}"
        )
    print("PASS output core: both ambient endpoint Jacobians have exact index 18")
    print(f"PASS wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

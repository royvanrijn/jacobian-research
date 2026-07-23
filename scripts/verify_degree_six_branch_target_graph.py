#!/usr/bin/env python3
"""Verify the complete stable-target graph in the degree-six cluster chart.

For three moving branch values with 0 and infinity fixed, the open target
moduli is the complement of six lines in P^2.  Blowing up its four triple
points gives Mbar_0,5.  The first three blowups are the toric/radial
permutohedral layer; the fourth is the triple-Maxwell refinement.
"""

from __future__ import annotations

import json
import sys
from collections import deque
from itertools import combinations, product
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = (
    ROOT / "artifacts" / "generated-results" / "branch_target_graph_degree6.json"
)


def primitive_projective(vector: sp.Matrix) -> tuple[int, int, int]:
    """Normalize a nonzero rational P^2 vector to a primitive integer tuple."""

    entries = [sp.Rational(entry) for entry in vector]
    common_denominator = sp.ilcm(*(entry.q for entry in entries))
    integers = [int(entry * common_denominator) for entry in entries]
    common_divisor = abs(sp.igcd(*integers))
    integers = [entry // common_divisor for entry in integers]
    first_nonzero = next(entry for entry in integers if entry)
    if first_nonzero < 0:
        integers = [-entry for entry in integers]
    return tuple(integers)


# Homogeneous branch-value coordinates [A:B:C] represent the five marked
# target points (0, infinity, A, B, C) modulo target scaling.
A, B, C = sp.symbols("A B C")
target_coordinates = sp.Matrix([A, B, C])
line_equations = {
    "A=0": sp.Matrix([1, 0, 0]),
    "B=0": sp.Matrix([0, 1, 0]),
    "C=0": sp.Matrix([0, 0, 1]),
    "A=B": sp.Matrix([1, -1, 0]),
    "A=C": sp.Matrix([1, 0, -1]),
    "B=C": sp.Matrix([0, 1, -1]),
}

# Recover every intersection point and the full incidence directly from the
# six line equations.
intersection_points: set[tuple[int, int, int]] = set()
for left, right in combinations(line_equations, 2):
    point = primitive_projective(
        line_equations[left].cross(line_equations[right])
    )
    intersection_points.add(point)

incidence = {}
for point in sorted(intersection_points):
    point_vector = sp.Matrix(point)
    incidence[point] = sorted(
        name
        for name, coefficients in line_equations.items()
        if coefficients.dot(point_vector) == 0
    )

triple_points = {
    point: lines for point, lines in incidence.items() if len(lines) == 3
}
double_points = {
    point: lines for point, lines in incidence.items() if len(lines) == 2
}
assert triple_points == {
    (0, 0, 1): ["A=0", "A=B", "B=0"],
    (0, 1, 0): ["A=0", "A=C", "C=0"],
    (1, 0, 0): ["B=0", "B=C", "C=0"],
    (1, 1, 1): ["A=B", "A=C", "B=C"],
}
assert double_points == {
    (0, 1, 1): ["A=0", "B=C"],
    (1, 0, 1): ["A=C", "B=0"],
    (1, 1, 0): ["A=B", "C=0"],
}

# Every arrangement line passes through two triple points and one residual
# double point.  After blowing up the four triple points, its strict transform
# has self-intersection 1-2=-1.  The four exceptional curves also have
# self-intersection -1, giving the ten boundary curves of Mbar_0,5.
line_triple_incidence = {
    line: sorted(point for point, lines in triple_points.items() if line in lines)
    for line in line_equations
}
line_double_incidence = {
    line: sorted(point for point, lines in double_points.items() if line in lines)
    for line in line_equations
}
assert all(len(points) == 2 for points in line_triple_incidence.values())
assert all(len(points) == 1 for points in line_double_incidence.values())

line_vertices = tuple(f"L:{name}" for name in line_equations)
exceptional_vertices = tuple(f"E:{point}" for point in triple_points)
boundary_vertices = line_vertices + exceptional_vertices
boundary_adjacency = {vertex: set() for vertex in boundary_vertices}

# Strict transforms meet over the three double points.
for lines in double_points.values():
    left, right = (f"L:{line}" for line in lines)
    boundary_adjacency[left].add(right)
    boundary_adjacency[right].add(left)

# Each exceptional curve meets the strict transforms of the three lines
# through its blown-up point.
for point, lines in triple_points.items():
    exceptional = f"E:{point}"
    for line in lines:
        line_vertex = f"L:{line}"
        boundary_adjacency[exceptional].add(line_vertex)
        boundary_adjacency[line_vertex].add(exceptional)

assert len(boundary_vertices) == 10
assert all(len(neighbors) == 3 for neighbors in boundary_adjacency.values())
assert sum(map(len, boundary_adjacency.values())) // 2 == 15


def graph_distance(start: str, finish: str) -> int:
    queue = deque([(start, 0)])
    visited = {start}
    while queue:
        vertex, distance = queue.popleft()
        if vertex == finish:
            return distance
        for neighbor in boundary_adjacency[vertex] - visited:
            visited.add(neighbor)
            queue.append((neighbor, distance + 1))
    raise AssertionError("boundary graph is disconnected")


# The dual boundary graph is the Petersen graph: it is connected,
# three-regular, triangle-free, has girth five, and diameter two.
for first, second in combinations(boundary_vertices, 2):
    assert graph_distance(first, second) <= 2
for first in boundary_vertices:
    for second in boundary_adjacency[first]:
        assert not (
            boundary_adjacency[first] & boundary_adjacency[second]
        ), "a triangle survived the four blowups"

has_five_cycle = False
for start in boundary_vertices:
    for second in boundary_adjacency[start]:
        for third in boundary_adjacency[second] - {start}:
            for fourth in boundary_adjacency[third] - {second}:
                for fifth in boundary_adjacency[fourth] - {third}:
                    if (
                        start in boundary_adjacency[fifth]
                        and len({start, second, third, fourth, fifth}) == 5
                    ):
                        has_five_cycle = True
assert has_five_cycle

# Blowing up four points in P^2 gives K^2=9-4=5.  Together with the ten
# (-1)-curves and the arrangement presentation, this is the degree-five del
# Pezzo/Kapranov model of Mbar_0,5.
assert 9 - len(triple_points) == 5
boundary_self_intersections = {
    **{vertex: -1 for vertex in line_vertices},
    **{vertex: -1 for vertex in exceptional_vertices},
}

# The first three target points are coordinate points.  Their blowups give
# Bl_{p_A,p_B,p_C} P^2, the toric permutohedral surface.  The diagonal point
# is the unique additional blowup required for the stable-target model.
coordinate_points = {(1, 0, 0), (0, 1, 0), (0, 0, 1)}
diagonal_point = (1, 1, 1)
assert set(triple_points) == coordinate_points | {diagonal_point}

# Pull the four target centers back through the leading branch-scale map
# [x:y:z] -> [-9/4*x^2 : -y^2 : -9*z^2].
x, y, z = sp.symbols("x y z")
leading_branch_values = (
    -sp.Rational(9, 4) * x**2,
    -y**2,
    -9 * z**2,
)
diagonal_preimages = []
for y_sign, z_sign in product((-1, 1), repeat=2):
    point = (sp.Integer(1), y_sign * sp.Rational(3, 2), z_sign * sp.Rational(1, 2))
    values = tuple(
        sp.factor(value.subs({x: point[0], y: point[1], z: point[2]}))
        for value in leading_branch_values
    )
    assert values[0] == values[1] == values[2] == -sp.Rational(9, 4)
    diagonal_preimages.append(point)
assert len(set(diagonal_preimages)) == 4

# The tangent cone of the exact triple-Maxwell locus is reduced at those four
# directions.  Hence after the radial blowup its four formal branches are
# separated and smooth.
tangent_equations = (4 * y**2 - 9 * x**2, x**2 - 4 * z**2)
for point in diagonal_preimages:
    assert all(equation.subs({x: point[0], y: point[1], z: point[2]}) == 0
               for equation in tangent_equations)
    affine_jacobian = sp.Matrix(
        [
            [sp.diff(equation.subs(x, 1), variable) for variable in (y, z)]
            for equation in tangent_equations
        ]
    )
    assert affine_jacobian.det().subs({y: point[1], z: point[2]}) != 0

# On the etale critical-point chart, each exact moving branch value is its
# cluster scale squared times a unit.  This validates the Kummer pullback,
# rather than only its leading valuation.
W, R = sp.symbols("W R")
family = sp.expand(
    W * (W - x) * (W - 1) * (W - 1 - y) * (W - 3) * (W - 3 - z)
)
cluster_data = (
    (sp.Integer(0), x, -sp.Rational(9, 4)),
    (sp.Integer(1), y, -sp.Integer(1)),
    (sp.Integer(3), z, -sp.Integer(9)),
)
for center, scale, expected_value_unit in cluster_data:
    derivative_quotient = sp.cancel(
        sp.diff(family, W).subs(W, center + scale * R) / scale
    )
    value_quotient = sp.cancel(
        family.subs(W, center + scale * R) / scale**2
    )
    assert scale not in sp.denom(derivative_quotient).free_symbols
    assert scale not in sp.denom(value_quotient).free_symbols
    limiting_derivative = sp.factor(
        derivative_quotient.subs({x: 0, y: 0, z: 0})
    )
    limiting_value = sp.factor(value_quotient.subs({x: 0, y: 0, z: 0}))
    assert limiting_derivative.subs(R, sp.Rational(1, 2)) == 0
    assert sp.diff(limiting_derivative, R).subs(R, sp.Rational(1, 2)) != 0
    assert limiting_value.subs(R, sp.Rational(1, 2)) == expected_value_unit

artifact = {
    "experiment": "complete degree-six stable-target graph",
    "open_target_moduli": {
        "coordinates": "[lambda_x:lambda_y:lambda_z] in P^2",
        "removed_lines": list(line_equations),
    },
    "arrangement": {
        "triple_points": {
            str(point): lines for point, lines in sorted(triple_points.items())
        },
        "double_points": {
            str(point): lines for point, lines in sorted(double_points.items())
        },
    },
    "compactifications": {
        "radial": "Blow up the three coordinate points of P^2",
        "stable_target": "Blow up the three coordinate points and [1:1:1]",
        "stable_target_identification": "Mbar_0,5",
        "canonical_square": 5,
        "boundary_curve_count": 10,
        "boundary_dual_graph": "Petersen graph",
    },
    "source_pullback": {
        "leading_map": "[x:y:z] -> [-9/4*x^2:-y^2:-9*z^2]",
        "diagonal_target_point": list(diagonal_point),
        "diagonal_preimages": [
            [str(coordinate) for coordinate in point]
            for point in diagonal_preimages
        ],
        "exact_resonance_ideal_on_critical_chart": [
            "lambda_x-lambda_y",
            "lambda_z-lambda_y",
        ],
        "conclusion": (
            "after the radial modification, blow up the four separated "
            "triple-Maxwell branches to obtain the complete stable-target graph"
        ),
    },
    "boundary_adjacency": {
        vertex: sorted(neighbors)
        for vertex, neighbors in sorted(boundary_adjacency.items())
    },
    "boundary_self_intersections": boundary_self_intersections,
}

expected_artifact = json.dumps(artifact, indent=2, sort_keys=True) + "\n"
if "--emit-json" in sys.argv:
    print(expected_artifact, end="")
    raise SystemExit(0)
assert ARTIFACT.read_text() == expected_artifact, (
    f"{ARTIFACT.relative_to(ROOT)} is stale; regenerate it from this script"
)

print("PASS target arrangement: six lines, four triple points, three double points")
print("PASS Kapranov surface: Bl_4(P^2), K^2=5, ten (-1)-curves")
print("PASS boundary incidence: Petersen graph")
print("PASS source pullback: four reduced triple-Maxwell directions")
print("PASS exact critical chart: lambda_i=x_i^2*unit")
print("DEGREE_SIX_BRANCH_TARGET_GRAPH_PASS")

#!/usr/bin/env python3
"""Verify the second-order obstruction in the final HC4 affine-plane frame.

This research script starts with the exact first-order system used by
``verify_hc4_affine_plane_bridge.py``.  It parameterizes every admissible
connection jet, differentiates that linear system in all four frame
directions, imposes zero curvature of the ambient affine connection, and
eliminates the derivative parameters.  The output is the complete quadratic
compatibility ideal on the first-order connection coefficients.

The first flatness prolongation forces the non-Schubert scalar to one of two
signs.  The HC4 maximal-motion identity makes the determinant ``p*q`` of the
selected Gauss-kernel-line differential a nonzero constant.  Prolonging with
``d(p*q)=0`` then makes the compatibility ideal the unit ideal on ``a != 0``.
This excludes the Gauss-rank-two maximal-motion branch.  The script also
certifies the exact lower flag-motion tensors and the rank-one split used by
the global hyperplane-pencil argument in ``HC4RSD80``.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "artifacts" / "generated-results" / "hc4_affine_plane_prolongation.json"
EXPECTED_OUT_SHA256 = "d258b2b9be7a0906bae70a317044ba2011ed5dea6d2b8d765f8d0571a0217cf0"


def audit_existing() -> None:
    actual = hashlib.sha256(OUT.read_bytes()).hexdigest()
    assert actual == EXPECTED_OUT_SHA256, (actual, EXPECTED_OUT_SHA256)
    payload = json.loads(OUT.read_text(encoding="utf-8"))
    assert payload["scope"] == (
        "second-order flatness and lower flag geometry after HC4RSD77"
    )
    assert payload["first_order_unknowns"] == 64
    assert payload["first_order_rank"] == 47
    assert payload["first_order_parameters"] == 17
    assert payload["curvature_equations"] == 96
    assert payload["derivative_unknowns"] == 68
    assert payload["derivative_rank"] == 48
    assert payload["compatibility_equations"] == 4
    assert payload["first_order_witness_survives_flatness"] is True
    assert payload["maximal_motion"]["saturation_by_a"] == "unit ideal"
    assert payload["maximal_motion"]["flat_witness_survives"] is False
    assert payload["canonical_frame_gauge"]["residual_group"] == "+/- I"
    assert "local rank-one split" in payload["proof_boundary"]
    assert "proved in the companion note" in payload["proof_boundary"]
    print(
        "PASS: committed HC4 affine-plane prolongation artifact is intact and "
        "retains its local-proof boundary; no symbolic replay or rewrite"
    )


parser = argparse.ArgumentParser()
parser.add_argument(
    "--audit-existing-only",
    action="store_true",
    help="validate the committed artifact without symbolic replay or rewriting it",
)
arguments = parser.parse_args()
if arguments.audit_existing_only:
    audit_existing()
    raise SystemExit(0)

import sympy as sp

n = 4
Gamma: dict[tuple[int, int, int], sp.Symbol] = {}
variables: list[sp.Symbol] = []
for i in range(n):
    for j in range(n):
        for k in range(n):
            symbol = sp.symbols(f"g{i+1}{j+1}{k+1}")
            Gamma[i, j, k] = symbol
            variables.append(symbol)

S = sp.zeros(n)
for i in range(n):
    S[i, n - 1 - i] = 1
N = sp.zeros(n)
for j in range(1, n):
    N[j - 1, j] = 1
T = S * N


def cov(matrix: sp.Matrix, i: int, j: int, k: int) -> sp.Expr:
    return -sum(
        Gamma[i, j, a] * matrix[a, k]
        + Gamma[i, k, a] * matrix[j, a]
        for a in range(n)
    )


first_order_equations: list[sp.Expr] = []
for matrix in (S, T):
    for i in range(n):
        for j in range(n):
            for k in range(n):
                component = cov(matrix, i, j, k)
                first_order_equations.append(component - cov(matrix, j, i, k))
                first_order_equations.append(component - cov(matrix, i, k, j))

# Frobenius of ker N^2 and ker N^3.
for k in (2, 3):
    first_order_equations.append(Gamma[0, 1, k] - Gamma[1, 0, k])
for i in range(3):
    for j in range(i + 1, 3):
        first_order_equations.append(Gamma[i, j, 3] - Gamma[j, i, 3])

# Primitive quasi-translation kernel and constant affine frame volume.
for k in range(n):
    first_order_equations.append(Gamma[0, 0, k])
for i in range(n):
    first_order_equations.append(sum(Gamma[i, j, j] for j in range(n)))

first_order_equations = [
    sp.expand(e) for e in first_order_equations if sp.expand(e) != 0
]
linear_matrix, _ = sp.linear_eq_to_matrix(first_order_equations, variables)
nullspace = linear_matrix.nullspace()
parameters = sp.symbols(f"p0:{len(nullspace)}")

# Gamma is a row convention: Gamma_i[j,k] is the e_k coefficient of
# nabla_{e_i} e_j.
gamma_vector = sp.zeros(len(variables), 1)
for parameter, basis_vector in zip(parameters, nullspace):
    gamma_vector += parameter * basis_vector
gamma_substitution = {
    variable: sp.expand(gamma_vector[index])
    for index, variable in enumerate(variables)
}

connection = []
for i in range(n):
    connection.append(
        sp.Matrix(
            n,
            n,
            lambda j, k: gamma_substitution[Gamma[i, j, k]],
        )
    )

# Directional derivatives of the 17 free first-order parameters.
dparameters: dict[tuple[int, int], sp.Symbol] = {}
derivative_variables: list[sp.Symbol] = []
for i in range(n):
    for alpha in range(len(parameters)):
        symbol = sp.symbols(f"dp{i+1}_{alpha}")
        dparameters[i, alpha] = symbol
        derivative_variables.append(symbol)

derivative_connection = []
for i in range(n):
    directional = []
    for j in range(n):
        matrix = sp.zeros(n)
        for alpha, basis_vector in enumerate(nullspace):
            for row in range(n):
                for column in range(n):
                    variable_index = variables.index(Gamma[j, row, column])
                    matrix[row, column] += (
                        dparameters[i, alpha] * basis_vector[variable_index]
                    )
        directional.append(matrix)
    derivative_connection.append(directional)

curvature_equations: list[sp.Expr] = []
for i in range(n):
    for j in range(i + 1, n):
        bracket = [
            connection[i][j, a] - connection[j][i, a]
            for a in range(n)
        ]
        curvature = (
            derivative_connection[i][j]
            - derivative_connection[j][i]
            + connection[j] * connection[i]
            - connection[i] * connection[j]
            - sum((bracket[a] * connection[a] for a in range(n)), sp.zeros(n))
        )
        curvature_equations.extend(sp.expand(entry) for entry in curvature)

def parameter_expression(symbol: sp.Symbol) -> sp.Expr:
    return gamma_substitution[symbol]


a = sp.factor(parameter_expression(Gamma[2, 0, 1]))
b = sp.factor(parameter_expression(Gamma[3, 0, 1]))
q = sp.factor(parameter_expression(Gamma[3, 1, 3]))
r = sp.factor(parameter_expression(Gamma[3, 1, 2]))
p0_geometric = sp.factor(parameter_expression(Gamma[2, 2, 3]))
s = sp.factor(parameter_expression(Gamma[3, 2, 3]))

# Exact flag-motion tensors.  Rows in source_kernel_projective_derivative are
# the e2,e3,e4 components and columns are the e1,...,e4 directions.  The A_i
# matrices are the derivatives of E2=<e1,e2> in Hom(E2,V/E2).  The E3 matrix
# is the affine second fundamental form of E3=<e1,e2,e3> along E3.
source_kernel_projective_derivative = sp.Matrix(
    3,
    4,
    lambda row, direction: connection[direction][0, row + 1],
)
middle_plane_A3 = sp.Matrix(
    2,
    2,
    lambda row, column: connection[2][column, row + 2],
)
middle_plane_A4 = sp.Matrix(
    2,
    2,
    lambda row, column: connection[3][column, row + 2],
)
top_hyperplane_second_fundamental = sp.Matrix(
    3,
    3,
    lambda left, right: connection[left][right, 3],
)
top_hyperplane_transverse_derivative = sp.Matrix(
    1,
    3,
    lambda _, column: connection[3][column, 3],
)

assert source_kernel_projective_derivative == sp.Matrix([
    [0, 0, a, b],
    [0, 0, 0, a],
    [0, 0, 0, 0],
])
assert middle_plane_A3 == sp.Matrix([
    [0, -(a + q) / 2],
    [0, 0],
])
assert middle_plane_A4 == sp.Matrix([
    [a, r],
    [0, q],
])
assert top_hyperplane_second_fundamental == sp.diag(0, 0, p0_geometric)
assert top_hyperplane_transverse_derivative == sp.Matrix([[0, q, s]])


def derivative_of(expression: sp.Expr, direction: int) -> sp.Expr:
    return sp.expand(
        sum(
            sp.diff(expression, parameter) * dparameters[direction, alpha]
            for alpha, parameter in enumerate(parameters)
        )
    )


def eliminate_derivatives(equations: list[sp.Expr]):
    matrix, rhs = sp.linear_eq_to_matrix(equations, derivative_variables)
    left_kernel = matrix.T.nullspace()
    conditions = []
    for vector in left_kernel:
        equation = sp.factor((vector.T * rhs)[0])
        if equation != 0 and equation not in conditions and -equation not in conditions:
            conditions.append(equation)
    return matrix, conditions


derivative_matrix, compatibility = eliminate_derivatives(curvature_equations)

# For ell=S*e1, the fixed components of S give
#
#   nabla_i ell = -sum_j Gamma^4_{i,j} theta^j.
#
# On the two transverse directions e3,e4 and modulo ell, the two relevant
# target directions theta^3,theta^2 therefore give the triangular matrix
# displayed below.  Its determinant is p*q.  HC4RSD72 identifies this
# determinant, up to the fixed nonzero canonical-frame factor, with the
# nonzero constant Hessian determinant.
selected_line_derivative = sp.Matrix([
    [-connection[2][2, 3], -connection[3][2, 3]],
    [-connection[2][1, 3], -connection[3][1, 3]],
])
assert selected_line_derivative[0, 0] == -p0_geometric
assert selected_line_derivative[1, 0] == 0
assert selected_line_derivative[1, 1] == -q
assert sp.factor(selected_line_derivative.det() - p0_geometric * q) == 0

maximal_motion_equations = curvature_equations + [
    derivative_of(p0_geometric * q, direction) for direction in range(n)
]
maximal_derivative_matrix, maximal_compatibility = eliminate_derivatives(
    maximal_motion_equations
)

# On a generic rank-one stratum of the projective source-kernel differential,
# a=0 and b!=0.  Flatness already forces q=0 after a=0; impose that stratum
# identity and its first derivatives before eliminating the remaining jets.
rank_one_substitution = {a: 0, q: 0}
rank_one_equations = [
    sp.expand(equation.subs(rank_one_substitution))
    for equation in curvature_equations
] + [
    derivative_of(a, direction) for direction in range(n)
] + [
    derivative_of(q, direction) for direction in range(n)
]
rank_one_derivative_matrix, rank_one_compatibility = eliminate_derivatives(
    rank_one_equations
)

# Three displayed generators suffice for the contradiction on a != 0.
flatness_gate = [
    sp.expand(a * (p0_geometric - q)),
    sp.expand(
        4 * p0_geometric * a - 3 * a**2 - 4 * a * q + 3 * q**2
    ),
    sp.expand(
        p0_geometric
        * (2 * p0_geometric * a - a * q + 3 * q**2)
    ),
]


def is_nonzero_scalar_multiple(left: sp.Expr, right: sp.Expr) -> bool:
    if right == 0:
        return False
    ratio = sp.cancel(left / right)
    return ratio != 0 and not ratio.free_symbols


for expected in flatness_gate:
    assert any(
        is_nonzero_scalar_multiple(expected, equation)
        for equation in maximal_compatibility
    )

inverse_a = sp.symbols("inverse_a")
unit_basis = sp.groebner(
    flatness_gate + [inverse_a * a - 1],
    inverse_a,
    p0_geometric,
    q,
    a,
    order="lex",
)
assert list(unit_basis.polys) == [sp.Poly(1, *unit_basis.gens)]

rank_one_gate = sp.expand(p0_geometric * r)
assert any(
    is_nonzero_scalar_multiple(rank_one_gate, equation)
    for equation in rank_one_compatibility
)
assert all(
    not (equation.free_symbols & b.free_symbols)
    for equation in rank_one_compatibility
)

# Lower-motion geometry used by HC4RSD80.  When a=q=0, b=0 fixes the source
# kernel line projectively.  If b!=0, flatness gives p0*r=0.  On r=0 the
# middle plane E2 is parallel in every ambient direction.  On p0=0 the
# distribution E3 is autoparallel, so its leaves are affine hyperplanes.
rank_zero_substitution = {a: 0, b: 0, q: 0}
assert source_kernel_projective_derivative.subs(rank_zero_substitution) == sp.zeros(3, 4)

parallel_middle_substitution = {a: 0, q: 0, r: 0}
# For every ambient direction and both E2 basis vectors, both normal
# components vanish.
assert all(
    sp.expand(connection[direction][vector, normal].subs(parallel_middle_substitution))
    == 0
    for direction in range(4)
    for vector in range(2)
    for normal in (2, 3)
)

affine_hyperplane_substitution = {p0_geometric: 0}
assert (
    top_hyperplane_second_fundamental.subs(affine_hyperplane_substitution)
    == sp.zeros(3)
)

rank_one_flag_substitution = {a: 0, q: 0, p0_geometric: 0}
assert source_kernel_projective_derivative.subs(
    rank_one_flag_substitution
)[:, 3] == sp.Matrix([b, 0, 0])
assert middle_plane_A4.subs(rank_one_flag_substitution) == sp.Matrix([
    [0, r],
    [0, 0],
])

# The residual gauge of an S-adapted regular-nilpotent frame is only a sign.
# A centralizer element is a polynomial in N; S-orthogonality is G^2=I.
c0, c1, c2, c3 = sp.symbols("c0 c1 c2 c3")
centralizer = c0 * sp.eye(n) + c1 * N + c2 * N**2 + c3 * N**3
centralizer_square = sp.expand(centralizer**2 - sp.eye(n))
centralizer_equations = []
for entry in centralizer_square:
    entry = sp.factor(entry)
    if entry != 0 and entry not in centralizer_equations:
        centralizer_equations.append(entry)
expected_centralizer_equations = [
    (c0 - 1) * (c0 + 1),
    2 * c0 * c1,
    2 * c0 * c2 + c1**2,
    2 * (c0 * c3 + c1 * c2),
]
assert len(centralizer_equations) == len(expected_centralizer_equations)
assert all(
    sp.expand(actual - expected) == 0
    for actual, expected in zip(
        centralizer_equations, expected_centralizer_equations
    )
)

# Test the formal first-order witness a=q=1, r=0 with every other nullspace
# parameter zero.  In the present canonical nullspace a and q are individual
# parameters; solve the linear parameter conditions rather than relying on
# their positions.
witness_conditions = [a - 1, q - 1, r, p0_geometric - 1]
witness_solution = sp.linsolve(witness_conditions, parameters)
witness_tuple = next(iter(witness_solution))
witness_substitution = {}
free_witness_symbols = set().union(*(entry.free_symbols for entry in witness_tuple))
for symbol in free_witness_symbols:
    if symbol in parameters:
        witness_substitution[symbol] = 0
witness_tuple = tuple(sp.expand(entry.subs(witness_substitution)) for entry in witness_tuple)
witness_parameter_substitution = dict(zip(parameters, witness_tuple))
witness_obstructions = [
    sp.factor(equation.subs(witness_parameter_substitution))
    for equation in compatibility
]
witness_obstructions = [equation for equation in witness_obstructions if equation != 0]
maximal_witness_obstructions = [
    sp.factor(equation.subs(witness_parameter_substitution))
    for equation in maximal_compatibility
]
maximal_witness_obstructions = [
    equation for equation in maximal_witness_obstructions if equation != 0
]

equations_involving_aq = []
a_symbols = a.free_symbols
q_symbols = q.free_symbols
for equation in compatibility:
    if equation.free_symbols & a_symbols and equation.free_symbols & q_symbols:
        equations_involving_aq.append(str(equation))

result = {
    "scope": "second-order flatness and lower flag geometry after HC4RSD77",
    "status": (
        "Gauss-rank-two maximal motion excluded; local lower-motion split "
        "and flag tensors certified"
    ),
    "first_order_unknowns": len(variables),
    "first_order_rank": int(linear_matrix.rank()),
    "first_order_parameters": len(parameters),
    "curvature_equations": len(curvature_equations),
    "derivative_unknowns": len(derivative_variables),
    "derivative_rank": int(derivative_matrix.rank()),
    "compatibility_equations": len(compatibility),
    "distinguished_parameters": {
        "a": str(a),
        "b": str(b),
        "q": str(q),
        "r": str(r),
        "p0": str(p0_geometric),
        "s": str(s),
    },
    "compatibility_involving_a_and_q": equations_involving_aq,
    "first_order_witness_survives_flatness": not witness_obstructions,
    "first_order_witness_obstructions": [str(e) for e in witness_obstructions],
    "compatibility_generators": [str(e) for e in compatibility],
    "maximal_motion": {
        "selected_line_derivative": str(selected_line_derivative),
        "selected_line_determinant": str(p0_geometric * q),
        "additional_equations": (
            "d(p0*q)=0 from p0*q=delta up to a fixed nonzero frame factor"
        ),
        "derivative_rank": int(maximal_derivative_matrix.rank()),
        "compatibility_equations": len(maximal_compatibility),
        "compatibility_generators": [str(e) for e in maximal_compatibility],
        "flat_witness_survives": not maximal_witness_obstructions,
        "flat_witness_obstructions": [str(e) for e in maximal_witness_obstructions],
        "contradiction_generators": [str(e) for e in flatness_gate],
        "saturation_by_a": "unit ideal",
        "conclusion": (
            "a!=0 gives p0=q and q^2=a^2 from flatness, while "
            "d(p0*q)=0 gives a=-3q; characteristic zero forces "
            "q=0, a contradiction"
        ),
    },
    "canonical_frame_gauge": {
        "centralizer_form": "c0 I+c1 N+c2 N^2+c3 N^3",
        "orthogonality_equations": [str(e) for e in centralizer_equations],
        "residual_group": "+/- I",
    },
    "source_kernel_motion_rank_one": {
        "stratum": "a=0, q=0, b!=0",
        "derivative_rank": int(rank_one_derivative_matrix.rank()),
        "compatibility_equations": len(rank_one_compatibility),
        "compatibility_generators": [str(e) for e in rank_one_compatibility],
        "compatibility_involving_b": [
            str(e) for e in rank_one_compatibility if e.free_symbols & b.free_symbols
        ],
        "conclusion": "flatness leaves p0*r=0 and no equation involving b",
    },
    "lower_motion_flag_geometry": {
        "source_kernel_projective_derivative": str(
            source_kernel_projective_derivative
        ),
        "middle_plane_A3": str(middle_plane_A3),
        "middle_plane_A4": str(middle_plane_A4),
        "top_hyperplane_second_fundamental": str(
            top_hyperplane_second_fundamental
        ),
        "top_hyperplane_transverse_derivative": str(
            top_hyperplane_transverse_derivative
        ),
        "rank_zero": "a=b=q=0 fixes the projective source-kernel line",
        "r_zero": "a=q=r=0 makes E2 parallel in every ambient direction",
        "p0_zero": "p0=0 makes E3 autoparallel with affine-hyperplane leaves",
        "moving_p0_zero_flag": (
            "for a=q=p0=0, the e4 derivatives have e1->b*e2 and "
            "E2/E1->r*(E3/E2)"
        ),
    },
    "proof_boundary": (
        "The checker supplies the local rank-one split and flag tensors. "
        "HC4RSD80 closes p0=0 globally by the degree-one hyperplane-incidence "
        "lemma; that projective argument is proved in the companion note."
    ),
}
OUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
print(json.dumps(result, indent=2))

#!/usr/bin/env sage -python
"""Compile an exact QQ 3A3 -> A3+2A2 replacement for q4/o323.

The lifted branch 16 has P.O=1.  Smooth saturation at its simple pole gives
a 5 -> 3 linear Riemann--Roch space.  At the compact t=0 I4, the constant
term of the chord quartic splits as two squared linear traces.  Exactly one
trace leaves two square linear factors and a quartic squareclass.  Binary
quartic invariants then give I4+2I3+14I1.  No elimination or Groebner basis
is used.
"""

import hashlib
import json
import time
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, ZZ, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
GENERATED = ROOT / "artifacts/generated-results"
COMPACT = LOCAL / "q4o208-compact-weierstrass-qq.json"
LIFTS = LOCAL / "q4o208-q4o323-horizontal-resolved-qq.json"
SHELL = LOCAL / "q4o208-q4o323-horizontal-marking-qq.json"
PHYSICAL = GENERATED / "elkies-k3-h3-q4o208-physical-3a3-marking.json"
EDGE = GENERATED / "elkies-k3-h3-q4o208-physical-q4o1599-a3-2a2-certificate.json"
OUTPUT = LOCAL / "q4o208-q4o1599-a3-2a2-rr-qq.json"
INPUTS = (COMPACT, LIFTS, SHELL, PHYSICAL, EDGE)
started = time.monotonic()


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def bits(values):
    answer = 0
    for value in values:
        value = QQ(value)
        answer = max(answer, abs(ZZ(value.numerator())).nbits(), ZZ(value.denominator()).nbits())
    return int(answer)


def coefficients(poly):
    return [str(value) for value in poly.list()]


compact = json.loads(COMPACT.read_text())
lifts = json.loads(LIFTS.read_text())
shell = json.loads(SHELL.read_text())
physical = json.loads(PHYSICAL.read_text())
edge = json.loads(EDGE.read_text())
assert compact["status"] == "PASS_EXACT_QQ_Q4O208_COMPACT_WEIERSTRASS_NORMALIZATION"
assert lifts["status"] == "PASS_EXACT_QQ_Q4O323_RESOLVED_SIMPLE_POLE_HORIZONTAL"
assert shell["status"] == "PASS_EXACT_QQ_Q4O323_LIFTED_SHELL_EXCLUDES_TARGET"
assert physical["status"] == "PASS_EXACT_Q4O208_PHYSICAL_3A3_MARKING"
assert edge["status"] == "PASS_EXACT_MARKED_DEGREE_TWO_CANDIDATE_CERTIFICATE"
assert edge["candidate_id"]["label"] == "q4o1599-exact-qq-a3-2a2"
assert edge["child"]["root_data"] == [7, 24, 36]

R = PolynomialRing(QQ, "t")
t = R.gen()
A = R([QQ(value) for value in compact["compact_model"]["A_coefficients_low_to_high"]])
B = R([QQ(value) for value in compact["compact_model"]["B_coefficients_low_to_high"]])
record = next(
    item for item in lifts["exact_QQ_horizontal_sections"]
    if int(item["branch_index"]) == 16
)
Z = R([QQ(value) for value in record["Z_coefficients_low_to_high"]])
X = R([QQ(value) for value in record["X_coefficients_low_to_high"]])
Y = R([QQ(value) for value in record["Y_coefficients_low_to_high"]])
assert Z.degree() == 1 and X.degree() == 6 and Y.degree() == 8
assert Y**2 == X**3 + A*X*Z**4 + B*Z**6

# Branch 16 has the same sign and marked NS class in every complete graph
# solution.  The residual twofold ambiguity concerns only unused shell
# sections, so it does not propagate to this replacement edge.
branch_matches = [
    next(
        item for item in graph_solution["branch_matches"]
        if int(item["stored_branch_index"]) == 16
    )
    for graph_solution in shell["lattice_match"]["all_graph_solution_branch_matches"]
]
assert len(branch_matches) == shell["lattice_match"]["complete_graph_solutions_before_target_pin"] == 2
assert all(item["equation_Y_sign_relative_to_stored"] == 1 for item in branch_matches)
branch_classes = {tuple(item["matched_NS_coordinates"]) for item in branch_matches}
assert len(branch_classes) == 1
branch_class = list(branch_classes)[0]
raw_to_physical = matrix(ZZ, physical["C5_child_basis_in_physical_3A3"])
branch_class_physical = vector(ZZ, branch_class) * raw_to_physical
physical_fibre = vector(ZZ, edge["source_to_child_basis"][0])
old_zero = vector(ZZ, [-1, 1] + [0] * 17)
vertical_residual = physical_fibre - old_zero - branch_class_physical
# For an old-fibre degree-two neighbour, D=O+P+V.  This checks the exact
# horizontal section used by the equation against the certified lattice
# fibre, including equality of their Mordell--Weil projections.
assert vertical_residual[1] == 0
assert not any(vertical_residual[11:])

# Ambient f=(a0+a1*t+a2*t^2+a3*t^3+b*Z*m)/Z^2 with
# m=(y+yP)/(x-xP).  Smoothness at Z=0 is a*X-b*Y == 0 mod Z^2.
smooth_rows = matrix(QQ, [
    [(t**degree*X).mod(Z**2)[order] for degree in range(4)]
    + [-Y.mod(Z**2)[order]]
    for order in range(2)
])
assert smooth_rows.rank() == 2
smooth_kernel = smooth_rows.right_kernel_matrix()
assert smooth_kernel.nrows() == 3

# Derive, rather than guess, the resolved t=0 trace.  The constant term of
# the chord radicand has exactly two squared linear branches.
P = PolynomialRing(QQ, ("p", "q", "r"))
p, q, r = P.gens()
parameters = (p, q, r)
L = PolynomialRing(P, "l")
l = L.gen()
ambient = [
    sum(parameters[row]*smooth_kernel[row, column] for row in range(3))
    for column in range(5)
]
a_general = sum(ambient[degree]*l**degree for degree in range(4))
b_general = ambient[4]


def lift_poly(poly):
    return L([P(value) for value in poly.list()])


A_l, X_l, Y_l, Z_l = map(lift_poly, (A, X, Y, Z))
raw_general = (
    a_general**4 - 6*X_l*a_general**2*b_general**2
    + 8*Y_l*a_general*b_general**3 - 3*X_l**2*b_general**4
    - 4*A_l*b_general**4*Z_l**4
)
general_after_pole, remainder = raw_general.quo_rem(Z_l**4)
assert not remainder
constant_factors = list(P(general_after_pole[0]).factor())
linear_traces = [factor for factor, exponent in constant_factors if factor.degree() == 1 and exponent == 2]
assert len(linear_traces) == 2


def compile_trace(trace):
    trace_row = vector(QQ, [trace.monomial_coefficient(value) for value in parameters])
    trace_kernel = matrix(QQ, [list(trace_row)]).right_kernel_matrix()
    h0 = trace_kernel*smooth_kernel
    assert h0.nrows() == 2 and h0.rank() == 2
    S = PolynomialRing(QQ, "u")
    u = S.gen()
    T = PolynomialRing(S, "t")
    tt = T.gen()
    a_u = sum(S(h0[0, degree] + u*h0[1, degree])*tt**degree for degree in range(4))
    b_u = S(h0[0, 4] + u*h0[1, 4])
    A_u, X_u, Y_u, Z_u = map(T, (A, X, Y, Z))
    raw = (
        a_u**4 - 6*X_u*a_u**2*b_u**2 + 8*Y_u*a_u*b_u**3
        - 3*X_u**2*b_u**4 - 4*A_u*b_u**4*Z_u**4
    )
    after_pole, pole_remainder = raw.quo_rem(Z_u**4)
    assert not pole_remainder
    factors = list(after_pole.factor())
    odd = [factor for factor, exponent in factors if exponent % 2]
    quartics = [factor for factor in odd if factor.degree() == 4]
    return h0, factors, quartics


compiled = [compile_trace(trace) for trace in linear_traces]
quartic_hits = [item for item in compiled if len(item[2]) == 1 and len(item[2][0].list()) == 5]
assert len(quartic_hits) == 1
h0, radicand_factors, quartic_list = quartic_hits[0]
quartic = quartic_list[0]
assert sorted((factor.degree(), int(exponent)) for factor, exponent in radicand_factors) == [
    (1, 2), (1, 2), (4, 1),
]

S = quartic.base_ring()
u = S.gen()
quartic_coefficients = list(quartic.list())
e, d, c, b, a = quartic_coefficients
I = 12*a*e - 3*b*d + c**2
J = 72*a*c*e + 9*b*c*d - 27*a*d**2 - 27*b**2*e - 2*c**3
A_child = S(-27*I)
B_child = S(-27*J)
Delta_child = S(-16*(4*A_child**3 + 27*B_child**2))
assert (A_child.degree(), B_child.degree(), Delta_child.degree()) == (8, 12, 20)
assert A_child.leading_coefficient() and B_child.leading_coefficient()

delta_factors = list(Delta_child.factor())
cubic_factors = [(factor, exponent) for factor, exponent in delta_factors if exponent == 3]
nodal_factors = [(factor, exponent) for factor, exponent in delta_factors if exponent == 1]
assert len(cubic_factors) == 2 and all(factor.degree() == 1 for factor, unused in cubic_factors)
assert len(nodal_factors) == 1 and nodal_factors[0][0].degree() == 14
assert nodal_factors[0][0].is_squarefree()
assert all(A_child.gcd(factor) == B_child.gcd(factor) == 1 for factor, unused in delta_factors)
# Degree 20 with nonzero leading A,B gives an I4 fibre at infinity.
assert (24-Delta_child.degree(), 8-A_child.degree(), 12-B_child.degree()) == (4, 0, 0)

payload = {
    "schema": "elkies-k3.h3-q4o208-q4o1599-a3-2a2-rr-qq.v1",
    "status": "PASS_EXACT_QQ_Q4O1599_A3_2A2_RR_AND_JACOBIAN",
    "replacement": {
        "replaces_failed_edge": "q4/o323",
        "orbit": 1599,
        "horizontal_stored_branch": 16,
        "equation_Y_sign_relative_to_stored": 1,
        "P_dot_O": 1,
        "physical_fibre": edge["source_to_child_basis"][0],
        "lattice_child_root_data": edge["child"]["root_data"],
        "complete_graph_solution_count": len(branch_matches),
        "branch_class_constant_across_complete_graph_solutions": True,
        "horizontal_NS_class_in_C5_child": list(map(int, branch_class)),
        "horizontal_NS_class_in_physical_3A3": list(map(int, branch_class_physical)),
        "physical_divisor_decomposition": {
            "identity": "D = old_zero + horizontal_section + vertical_residual",
            "vertical_residual_in_physical_3A3": list(map(int, vertical_residual)),
            "MW_tail_zero": True,
        },
        "marking_ambiguity": None,
    },
    "horizontal": {
        "Z_coefficients_low_to_high": coefficients(Z),
        "X_coefficients_low_to_high": coefficients(X),
        "Y_coefficients_low_to_high": coefficients(Y),
        "degrees_X_Y_Z": [int(X.degree()), int(Y.degree()), int(Z.degree())],
        "maximum_rational_bits": bits(Z.list()+X.list()+Y.list()),
        "exact_compact_weierstrass_identity": True,
    },
    "resolved_RR": {
        "ambient_basis": ["1/Z^2", "t/Z^2", "t^2/Z^2", "t^3/Z^2", "m/Z"],
        "chord": "m=(y+y_P)/(x-x_P)",
        "smooth_saturation_rows": [[str(value) for value in row] for row in smooth_rows.rows()],
        "selected_resolved_trace_method": "unique t=0 squared linear trace with quartic squareclass",
        "kernel_basis": [[str(value) for value in row] for row in h0.rows()],
        "dimensions": {"ambient": 5, "smooth": 3, "h0": 2},
        "raw_square_factor_degrees_and_exponents": [
            [int(factor.degree()), int(exponent)] for factor, exponent in radicand_factors
        ],
    },
    "quartic": {
        "coefficients_in_t_low_to_high": [coefficients(value) for value in quartic_coefficients],
        "coefficient_degrees_in_u": [int(value.degree()) for value in quartic_coefficients],
        "degree_in_t": 4,
        "maximum_rational_bits": bits(
            coefficient for value in quartic_coefficients for coefficient in value.list()
        ),
    },
    "child": {
        "minimal_A_coefficients_low_to_high": coefficients(A_child),
        "minimal_B_coefficients_low_to_high": coefficients(B_child),
        "degrees_A_B_Delta": [8, 12, 20],
        "finite_reducible_fibres": [
            {"factor": str(factor.monic()), "kodaira": "I3", "delta_order": 3}
            for factor, unused in cubic_factors
        ],
        "finite_nodal_factor": str(nodal_factors[0][0].monic()),
        "finite_nodal_factor_degree": 14,
        "infinity": {"kodaira": "I4", "orders_A_B_Delta": [0, 0, 4]},
        "fibre_profile": "I4+2I3+14I1",
        "ADE": "A3+2A2",
        "root_data": [7, 24, 36],
        "MW_rank_if_rho19": 10,
        "euler_number": 24,
        "maximum_A_B_rational_bits": bits(A_child.list()+B_child.list()),
    },
    "method": {
        "large_Groebner_required": False,
        "full_discriminant_factorization_required": False,
        "runtime_seconds": time.monotonic()-started,
    },
    "proof_boundary": (
        "Exact QQ h0=2 pencil, quartic, Jacobian, I4+2I3+14I1 fibre profile, "
        "Euler 24 and A3+2A2/MW10 conditional on rho=19. The q4/o1599 lattice "
        "edge and bidirectional NS transports are exact. Branch 16 has the same sign "
        "and NS class in both complete lifted-section graph solutions, so the residual "
        "ambiguity among unused shell sections does not affect pinned-R17 alignment."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in INPUTS],
        "sha256": {str(path.relative_to(ROOT)): sha256(path) for path in INPUTS},
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
print(
    "Q4O1599RRQQ|ambient=5|smooth=3|h0=2|quartic=4|"
    "fibres=I4+2I3+14I1|ADE=A3+2A2|status={}|output={}".format(
        payload["status"], OUTPUT
    )
)

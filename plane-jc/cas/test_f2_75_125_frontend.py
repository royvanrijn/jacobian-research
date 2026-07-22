#!/usr/bin/env python3
"""Regression tests for the exact F2 j=1 source skeleton."""

import sympy as sp

from f2_75_125_frontend import (
    chain_data,
    forced_edges,
    machine_certificate,
    normalized_terminal_edge,
)


chain = chain_data()
assert chain["degree_pair"] == [75, 125]
assert chain["initial_direction"] == [5, -1]
assert chain["translation_root_multiplicities"] == {"P": 6, "Q": 10}
assert chain["final_corner"] == ["7/5", 2]

initial, terminal = forced_edges()
assert initial.P_start == (75, 60)
assert initial.P_end == terminal.P_start == (21, 6)
assert initial.Q_end == terminal.Q_start == (35, 10)
assert terminal.P_end == (4, 1)
assert terminal.Q_end == (1, 0)

X, y = sp.symbols("X y")
P, Q, bracket = normalized_terminal_edge()
assert bracket == X**4
assert sp.Poly(P, X, y).terms() == [((21, 6), 1), ((4, 1), 1)]
assert sp.Poly(Q, X, y).terms() == [
    ((35, 10), sp.Rational(-9, 5)),
    ((18, 5), -3),
    ((1, 0), -1),
]

certificate = machine_certificate()
assert certificate["schema"] == "plane-jc.f2-75-125-residual.v1"
assert certificate["frontend_complete"] is False
assert len(certificate["residual_obligations"]) == 3
assert certificate["terminal_edge_normalization"]["de_rham_obstruction_rank"] == 0
assert certificate["terminal_edge_normalization"]["complete_supports"]["Q"] == [
    [1, 0], [18, 5], [35, 10]
]
assert len(certificate["terminal_edge_normalization"]["coefficient_system"]["equations"]) == 3
assert certificate["common_power_top_band"]["unresolved_layer_gap"] == 35
assert certificate["laurent_polygon_branches"]["known_branch_count"] is None

print("PASS: F2 chain arithmetic and Puiseux translation are exact")
print("PASS: all forced edge vertices are nonzero and lattice-consistent")
print("PASS: the normalized terminal type-I bracket is exactly X^4")
print("PASS: the machine certificate remains explicitly non-exhaustive")

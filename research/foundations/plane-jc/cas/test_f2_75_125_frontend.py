#!/usr/bin/env python3
"""Regression tests for the exact F2 j=1 source skeleton."""

import sympy as sp

from f2_75_125_frontend import (
    chain_data,
    forced_edges,
    machine_certificate,
    normalized_terminal_edge,
    terminal_kummer_characters,
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

characters = terminal_kummer_characters()
assert characters["P_X_characters"] == [1, 4]
assert characters["Q_X_characters"] == [0, 1, 3]
assert characters["target_character"] == 4
assert characters["trivial_character_descent"] is False

certificate = machine_certificate()
assert certificate["schema"] == "plane-jc.f2-75-125-residual.v3"
assert certificate["frontend_complete"] is False
assert len(certificate["residual_obligations"]) == 3
assert certificate["terminal_edge_normalization"]["de_rham_obstruction_rank"] == 0
assert certificate["terminal_edge_normalization"]["complete_supports"]["Q"] == [
    [1, 0], [18, 5], [35, 10]
]
assert len(certificate["terminal_edge_normalization"]["coefficient_system"]["equations"]) == 3
assert certificate["common_power_top_band"]["unresolved_layer_gap"] == 35
assert certificate["common_power_top_band"]["formal_top_layer"] == 40
assert certificate["common_power_top_band"]["missing_zero_layers"] == [39, 5]
assert certificate["common_power_top_band"]["band_chart"]["jacobian"] == (
    "[t,z]_(X,y)=-z"
)
assert certificate["common_power_top_band"]["Q_top_band"].startswith("-9/5")
assert certificate["laurent_polygon_branches"]["known_branch_count"] is None
assert certificate["terminal_kummer_characters"] == characters

print("PASS: F2 chain arithmetic and Puiseux translation are exact")
print("PASS: all forced edge vertices are nonzero and lattice-consistent")
print("PASS: the normalized terminal type-I bracket is exactly X^4")
print("PASS: the terminal mu_5 character profile blocks trivial Kummer descent")
print("PASS: the corrected Laurent chart has [t,z]=-z and top layer 40")
print("PASS: the machine certificate remains explicitly non-exhaustive")

#!/usr/bin/env python3
"""Regression checks for the F2 carrier-Wronskian exclusion."""

from __future__ import annotations

import importlib.util
from pathlib import Path


MODULE_PATH = Path(__file__).with_name(
    "verify_f2_75_125_carrier_wronskian.py"
)
SPEC = importlib.util.spec_from_file_location("f2_carrier_wronskian", MODULE_PATH)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError("cannot load the F2 carrier-Wronskian checker")
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


payload = MODULE.build_payload()
assert payload["status"] == "exact-finite-carrier-classification-row-still-live"
assert payload["target_monomial_change"]["forced_first_nonshear_descent"] == 36
assert payload["pre_target_shear_normalization"][
    "removable_target_shear_descents"
] == [5, 10, 15, 20, 25, 30, 35]
assert payload["carrier_target_completion"]["carrier_target_ray"] == [5, 36]
assert payload["carrier_target_completion"]["blowup_count"] == 12
assert payload["carrier_target_completion"][
    "residue_orders_on_left_carrier_right"
] == [-1, 0, 1]

cases = {row["case"]: row for row in payload["carrier_cases"]}
assert cases["squarefree_R"]["divisor_shape"]["finite_H_orders"] == [
    -7,
    -15,
    -7,
    -7,
]
assert cases["squarefree_R"]["forced_R"] == "(v^2-3*v+3)/25"
assert cases["squarefree_R"]["residue_degree"] == 3
assert cases["squarefree_R"]["normalized_carrier_residue_map"] == "1+1/(v-1)^3"
assert cases["double_root_R"]["divisor_shape"]["finite_H_orders"] == [
    -7,
    -15,
    -15,
]
assert cases["double_root_R"]["forced_double_root_equation"] == "rho^2-3*rho+1=0"
assert cases["double_root_R"]["residue_degree"] == 6
assert cases["double_root_R"]["map_identity"] == "g(v)=(729/125)*h(s)"

print("PASS: target monomials force the first nonshear coefficient at descent 36")
print("PASS: seven earlier rational coefficients are removable target shears")
print("PASS: twelve target blowups extract the generic carrier ray (5,36)")
print("PASS: the squarefree Wronskian has one cyclic-cubic carrier row")
print("PASS: the double-root Wronskian has two conjugate terminal-Belyi rows")

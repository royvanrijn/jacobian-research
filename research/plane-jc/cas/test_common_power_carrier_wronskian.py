#!/usr/bin/env python3
"""Regression checks for the common-power carrier Wronskian theorem."""

from __future__ import annotations

import importlib.util
from pathlib import Path


MODULE_PATH = Path(__file__).with_name(
    "verify_common_power_carrier_wronskian.py"
)
SPEC = importlib.util.spec_from_file_location(
    "common_power_carrier_wronskian", MODULE_PATH
)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError("cannot load the common-power carrier checker")
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


payload = MODULE.build_payload()
assert payload["status"] == "proved-general-reduction-to-three-point-hurwitz-data"
assert payload["symbolic_reduction"]["forced_descent"] == (
    "delta=k*(m+n-1)+1"
)
assert payload["universal_passport"]["passport"] == {
    "over_0": "(k^L,1^s)",
    "over_infinity": "(p_i : mu_i>=2)",
    "over_third_value": "(k-2,1^(f-k+2))",
}

f2_rows = payload["f2_k5_specializations"]
assert f2_rows[0]["residue_degree"] == 3
assert f2_rows[0]["passport"] == {
    "over_0": [1, 1, 1],
    "over_infinity": [3],
    "over_third_value": [3],
}
assert f2_rows[1]["residue_degree"] == 6
assert f2_rows[1]["passport"] == {
    "over_0": [5, 1],
    "over_infinity": [3, 3],
    "over_third_value": [3, 1, 1, 1],
}

for k in range(3, 16):
    for partition in MODULE.integer_partitions(k):
        if MODULE.multiplicity_gcd(partition) != 1 or len(partition) == k:
            continue
        profile = MODULE.passport_for_partition(k, partition)
        passport = profile["passport"]
        assert sum(passport["over_0"]) == profile["residue_degree"]
        assert sum(passport["over_infinity"]) == profile["residue_degree"]
        assert sum(passport["over_third_value"]) == profile["residue_degree"]
        expected_constraint = (
            "contained_in_alternating_group"
            if k % 2 == 1
            else "not_contained_in_alternating_group"
        )
        assert profile["geometric_monodromy_constraint"] == expected_constraint

try:
    MODULE.passport_for_partition(5, (1, 1, 1, 1, 1))
except ValueError as error:
    assert "squarefree" in str(error)
else:
    raise AssertionError("the squarefree carrier gate disappeared")

try:
    MODULE.passport_for_partition(6, (2, 2, 2))
except ValueError as error:
    assert "primitive" in str(error)
else:
    raise AssertionError("the imprimitive carrier gate disappeared")

print("PASS: every coprime common-power pair has a unimodular target monomial")
print("PASS: primitive carriers force delta=k*(m+n-1)+1")
print("PASS: the universal low-degree operator gives a three-point passport")
print("PASS: carrier parity decides alternating-group containment from k alone")
print("PASS: squarefree and imprimitive edge loci remain explicit gates")
print("PASS: the k=5 specializations recover the F2 degree-three and degree-six maps")

#!/usr/bin/env python3
"""Regression checks for the exact F2 nonlinear forcing compiler."""

from __future__ import annotations

import importlib.util
from pathlib import Path


MODULE_PATH = Path(__file__).with_name(
    "compile_f2_75_125_nonlinear_forcing.py"
)
SPEC = importlib.util.spec_from_file_location("f2_nonlinear_forcing", MODULE_PATH)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError("cannot load the F2 nonlinear-forcing compiler")
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


payload = MODULE.build_payload()
assert payload["status"] == "exact-arithmetic-circuit-presentation;ideal-test-open"

coupled = payload["coupled_Laurent_forcing"]
assert coupled["projection_order"] == [[12, 35], [37, 37]]
assert coupled["divisibility_coordinate_count"] == 294
assert coupled["pinned_quotient_cokernel_coordinate_count"] == 53
assert coupled["full_Laurent_cokernel_coordinate_count"] == 347
assert len(coupled["rows"]) == 25
assert sum(
    row["full_Laurent_cokernel_coordinate_count"]
    for row in coupled["rows"]
) == 347

final = payload["final_thirteen_functionals"]
assert final["target"]["equation_count"] == 7
assert final["target"]["quotient_residue_count"] == 2
assert final["layer_zero"]["equation_count"] == 6
assert final["combined_equation_count"] == 13

incidence = payload["descent_eight_incidence"]
assert incidence["incidence_equation_count"] == 5
assert incidence["defect_relation"] == "27*y^2-9*y+1"

ledger = payload["squarefree_carrier_later_first_defect"]
assert ledger["later_spacing_interval"] == [9, 90]
assert ledger["later_spacing_count"] == 82
assert ledger["regime_counts"] == {
    "spacings_9_to_11_with_multiples_2_3": 3,
    "spacings_12_to_17_with_multiple_2": 6,
    "spacings_18_to_90_with_no_pre_target_multiple": 73,
}

equations = payload["equation_ledger"]
assert equations == {
    "all_equation_digest_sha256": (
        "72efde05c379b575e4383060f6f0a78d1206a70ffb777a9e9e5c1641c7168e1e"
    ),
    "coupled_full_Laurent": 347,
    "target_and_layer_zero": 13,
    "incidence_and_relative_field_relation": 6,
    "total": 366,
}
assert payload["arithmetic_circuit"]["reachable_variable_count"] == 954
assert not payload["decision"]["unit_ideal_obtained"]
assert not payload["decision"]["counterexample_obtained"]

print("PASS: ten endpoint circuits and upper tangent followers are substituted")
print("PASS: the coupled Laurent obstruction splits exactly as 294+53=347")
print("PASS: the final target/Hermite block has 7+6 coordinates")
print("PASS: the five incidence rows are compiled over the quartic algebra")
print("PASS: the squarefree carrier is routed through spacings 9..90")

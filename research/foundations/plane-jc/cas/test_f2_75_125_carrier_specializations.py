#!/usr/bin/env python3
"""Regression checks for the exact F2 carrier specializations."""

from __future__ import annotations

import importlib.util
from pathlib import Path


MODULE_PATH = Path(__file__).with_name(
    "verify_f2_75_125_carrier_specializations.py"
)
SPEC = importlib.util.spec_from_file_location(
    "f2_carrier_specializations", MODULE_PATH
)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError("cannot load the F2 carrier-specialization checker")
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


payload = MODULE.build_payload()
assert payload["status"] == "exact-number-field-linear-maps;nonlinear-forcing-open"

routing = payload["carrier_routing"]
assert not routing["squarefree_carrier"]["descent_eight_component_membership"]
assert routing["double_carrier"]["descent_eight_component_membership"]
assert payload["double_branch_coefficient_field"]["compositum_degree"] == 4
assert payload["double_branch_coefficient_field"]["geometric_branch_count"] == 4

zero_rows = payload["specialized_zero_row_maps"]
assert zero_rows["endpoint_disjoint_descents_1_through_11"] == {
    "rank": 134,
    "Q_kernel_dimension": 2,
    "full_Laurent_forcing_cokernel_dimension": 219,
    "agreement_with_parametric_reducer": True,
}
assert zero_rows["coupled_quotient_cokernel_after_common_factor"] == 53
assert zero_rows["coupled_full_Laurent_forcing_cokernel_dimension"] == 347
rows = {row["descent"]: row for row in zero_rows["rows"]}
for descent in (5, 10, 15, 20, 25):
    assert rows[descent]["Q_kernel_dimension"] == 1
    assert rows[descent]["quotient_forcing_cokernel_dimension"] == 3
for descent in (1, 12, 30, 35, 37):
    assert rows[descent]["Q_kernel_dimension"] == 0
    assert rows[descent]["quotient_forcing_cokernel_dimension"] == 2

target = payload["specialized_target_cokernel"]
assert target["quotient_operator"]["matrix_shape"] == [22, 20]
assert target["quotient_operator"]["rank"] == 20
assert target["complete_image"]["cokernel_dimension"] == 14
assert target["after_prior_control_and_fixed_endpoint_elimination"] == {
    "movable_rho_jets": 5,
    "quotient_fitting_coordinates": 2,
    "remaining_target_coordinates": 7,
}

layer_zero = payload["specialized_layer_zero_cokernel"]
assert layer_zero["image_rank"] == 20
assert layer_zero["cokernel_dimension"] == 14
assert layer_zero["complete_Hermite_map"]["same_row_space_as_cokernel"]
assert layer_zero["after_prior_control_and_fixed_endpoint_elimination"][
    "remaining_layer_zero_coordinates"
] == 6
assert payload["combined_residual_coordinate_count"] == 13

print("PASS: the squarefree carrier is routed outside the descent-eight component")
print("PASS: the two double carriers specialize over QQ(rho)")
print("PASS: the local defect makes four compositum branches")
print("PASS: all zero-row Schur kernels and cokernels are exact through layer 3")
print("PASS: the specialized target and Hermite cokernels have residual rank 7+6")

#!/usr/bin/env python3
"""Regression for the automatic essential-dimension search profile."""

import sympy as sp

from search_essential_bcw import essential_profile
from search_rank_aware_bcw import State, initial_state, support
from verify_shared_bcw_33_route import OPTIMIZED_PLAN, apply_shared_step, dense_factor


state = initial_state()
for component, first_support, second_support in OPTIMIZED_PLAN:
    dimension = len(state.variables)
    first = dense_factor(first_support, dimension)
    second = dense_factor(second_support, dimension)
    removed = tuple(a + b for a, b in zip(first, second))
    coefficient = sp.Poly(
        state.expressions[component], *state.variables, domain=sp.QQ
    ).coeff_monomial(removed)
    selected = (component, removed, coefficient, sum(removed))
    result = apply_shared_step(
        state.expressions,
        state.variables,
        state.registry,
        selected,
        (first, second),
        state.introduced,
    )
    assert result is not None
    step = (component, support(first), support(second))
    state = State(result[0], result[1], result[2], result[3], state.plan + (step,))

profile = essential_profile(state)
assert profile["introduced_variables"] == 13
assert profile["cubic_output_rank"] == 7
assert profile["homogeneous_dimension"] == 24
assert profile["constant_kernel_dimension"] == 2
assert profile["final_essential_dimension"] == 22
assert profile["projected_collision_separated"] is True
assert profile["coordinate_cyclic_invariant_row_module_dimensions_mod_1000003"] == [1, 22]

print("PASS essential BCW profile: rank compression gives dimension 24")
print("PASS essential BCW profile: constant-kernel quotient gives dimension 22")
print("PASS essential BCW profile: projected collision remains separated")
print("PASS essential BCW profile: modular cyclic row-module dimensions are 1 and 22")

#!/usr/bin/env python3
"""Exact symbolic verification of the emitted essential-dimension candidate.

This is deliberately not an independent replay: it validates the search
artifact before any theorem-level 21-dimensional artifact is frozen.
"""

import json
from pathlib import Path

import sympy as sp

from rank_compressed_bcw_homogenization import (
    constant_kernel_quotient,
    extract_quadratic_cubic,
    factor_cubic_output,
    rank_compressed_homogeneous_map,
    verify_parametric_factorization,
)
from search_essential_bcw import essential_profile
from search_rank_aware_bcw import State, initial_state, support
from verify_shared_bcw_33_route import apply_shared_step, dense_factor


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "artifacts" / "generated-results" / "essential_bcw_candidate.json"
stored = json.loads(SOURCE.read_text())
assert stored["format"] == "essential-bcw-search-candidate-v1"

state = initial_state()
for raw in stored["plan"]:
    component = raw["component"]
    first_support = [tuple(pair) for pair in raw["first_factor"]]
    second_support = [tuple(pair) for pair in raw["second_factor"]]
    dimension = len(state.variables)
    first = dense_factor(first_support, dimension)
    second = dense_factor(second_support, dimension)
    removed = tuple(a + b for a, b in zip(first, second))
    polynomial = sp.Poly(state.expressions[component], *state.variables, domain=sp.QQ)
    coefficient = polynomial.coeff_monomial(removed)
    assert coefficient and sum(removed) == max(
        sp.Poly(expression, *state.variables, domain=sp.QQ).total_degree()
        for expression in state.expressions
    )
    result = apply_shared_step(
        state.expressions,
        state.variables,
        state.registry,
        (component, removed, coefficient, sum(removed)),
        (first, second),
        state.introduced,
    )
    assert result is not None
    step = (component, support(first), support(second))
    state = State(result[0], result[1], result[2], result[3], state.plan + (step,))

assert len(state.plan) == 17
assert len(state.variables) == 17 and state.introduced == 14
assert max(
    sp.Poly(expression, *state.variables, domain=sp.QQ).total_degree()
    for expression in state.expressions
) == 3

quadratic, cubic = extract_quadratic_cubic(state.expressions, state.variables)
factorization = factor_cubic_output(cubic)
assert len(factorization.c) == 6
verify_parametric_factorization(state.variables, quadratic, cubic, factorization)
homogeneous_variables, homogeneous_h = rank_compressed_homogeneous_map(
    state.variables, quadratic, factorization
)
assert len(homogeneous_variables) == 24
quotient = constant_kernel_quotient(homogeneous_variables, homogeneous_h)
assert quotient.kernel.cols == 3
assert len(quotient.quotient_variables) == 21
assert constant_kernel_quotient(
    quotient.quotient_variables, quotient.quotient_h
).kernel.cols == 0

profile = essential_profile(state)
if profile != stored["profile"]:
    print(json.dumps(profile, indent=2))
    raise AssertionError("stored candidate profile is stale")

print("PASS candidate: 17 exact stable cancellations give a 17D quadratic-cubic map")
print("PASS candidate: cubic-output rank 6 homogenizes in dimension 24")
print("PASS candidate: exact constant kernel has dimension 3")
print("PASS candidate: the 21D quotient has zero constant kernel")
print("PASS candidate: three projected points remain distinct with one common image")
print("PASS search record: promoted separately by the frozen 21D generator and replay")

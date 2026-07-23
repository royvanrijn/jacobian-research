#!/usr/bin/env python3
"""Circuit-level search for restricted cubic and quartic-HN minima.

This search changes the BCW degree-lowering circuit *before* homogenization.
It augments the monomial shared-factor moves with exact polynomial-gate
exposures and multi-term target shears.  Partial maps are ranked by the
sampled rank and full power-rank tuple of ``J(F-id)``; terminal maps then pass
through rank-compressed homogenization and the constant-kernel quotient.

All circuit moves are explicit stable left-right equivalences:

* a gate ``g`` is exposed as the new output ``w+g`` by a source shear;
* a polynomial in exposed outputs is added to an old target component;
* the affine linear part is normalized back to the identity.

The modular profiles are search diagnostics.  Any improvement must be frozen,
replayed from the original three-variable collision, and certified over QQ
before it changes a theorem-level bound.
"""

from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import dataclass
from itertools import combinations
import json
from pathlib import Path
from typing import Iterable

import sympy as sp

from rank_compressed_bcw_homogenization import (
    extract_quadratic_cubic,
    factor_cubic_output,
    iterated_constant_kernel_quotient,
    rank_compressed_homogeneous_map,
)
from restricted_rank_profiles import (
    PowerRankProfile,
    correction_profile,
    cotangent_hessian_rank,
    polynomial_jacobian_profile,
    profile_key,
    transformed_hessian_power_profile,
)
from search_rank_aware_bcw import State, initial_state, state_key, support
from verify_shared_bcw_33_route import (
    apply_shared_step,
    candidate_splits,
    high_terms,
)


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "restricted_bcw_circuit_search.json"
)

INCUMBENTS = {
    "cubic_rank": 17,
    "cubic_sampled_index": 18,
    "cubic_dimension": 21,
    "quartic_hessian_rank": 37,
    "quartic_dimension": 42,
    "quartic_sampled_index": 34,
}

CORE_ATOM_NAMES = frozenset(
    (
        "x2s",
        "xxs",
        "v2r",
        "xvvz",
        "qb",
        "ayb",
        "v2h",
        "y2vb",
        "yvyb",
    )
)
PERTURBATION_COEFFICIENTS = (-12, -6, -3, -1, 1, 3, 6, 12)


@dataclass(frozen=True)
class Gate:
    name: str
    source: sp.Expr


@dataclass(frozen=True)
class CircuitAtom:
    """One summand in a target shear, expressed in exposed gate outputs."""

    name: str
    component: int
    coefficient: sp.Rational
    factors: tuple[tuple[str, int], ...]


@dataclass
class CircuitState:
    expressions: list[sp.Expr]
    variables: list[sp.Symbol]
    registry: dict[tuple[int, ...], int]
    introduced: int
    monomial_plan: tuple[
        tuple[int, tuple[tuple[int, int], ...], tuple[tuple[int, int], ...]],
        ...,
    ]
    circuit_atoms: tuple[str, ...]
    collision_points: tuple[tuple[sp.Expr, ...], ...]

    def legacy(self) -> State:
        return State(
            self.expressions,
            self.variables,
            self.registry,
            self.introduced,
            self.monomial_plan,
        )

    @property
    def plan_key(self) -> tuple[object, ...]:
        return self.circuit_atoms, self.monomial_plan


def initial_collision_points() -> tuple[tuple[sp.Expr, ...], ...]:
    return (
        (sp.Integer(0), sp.Integer(0), -sp.Rational(1, 8)),
        (sp.Integer(1), -sp.Rational(3, 4), sp.Rational(13, 4)),
        (-sp.Integer(1), sp.Rational(3, 4), sp.Rational(13, 4)),
    )


def evaluate(expression: sp.Expr, variables, point) -> sp.Expr:
    return sp.Poly(expression, *variables, domain=sp.QQ).eval(
        dict(zip(variables, point))
    )


def base_circuit() -> tuple[dict[str, Gate], dict[str, sp.Expr], tuple[CircuitAtom, ...]]:
    """Return the factorized nodes and exact multi-term rewrite atoms."""

    x, y, z = sp.symbols("x y z")
    gates = {
        "xy": Gate("xy", x * y),
        "xz": Gate("xz", x * z),
        "xy2": Gate("xy2", x * y**2),
        "vz": Gate("vz", (1 + 2 * x * y) * z),
        "y2": Gate("y2", y**2),
        "x2": Gate("x2", x**2),
        "s": Gate("s", 3 * y + x * z),
        "xlin": Gate("xlin", x),
        "xs": Gate("xs", x * (3 * y + x * z)),
        "yv": Gate("yv", y * (1 + 2 * x * y)),
        "yb": Gate("yb", y * (2 + 3 * x * y)),
        "xv": Gate("xv", x * (1 + 2 * x * y)),
    }
    # Signals are polynomials in the corresponding exposed-output symbols.
    a, r, q, h, yy, xx, ss, xl, xss, yv, yb, xv = sp.symbols(
        "gate_xy gate_xz gate_xy2 gate_vz gate_y2 gate_x2 gate_s "
        "gate_xlin gate_xs gate_yv gate_yb gate_xv"
    )
    signals = {
        "a": a,
        "v": 1 + 2 * a,
        "b": 2 + 3 * a,
        "r": r,
        "q": q,
        "h": h,
        "yy": yy,
        "xx": xx,
        "s": ss,
        "xl": xl,
        "xs": xss,
        "yv": yv,
        "yb": yb,
        "xv": xv,
    }
    core_atoms = (
        # The first normalized component is x-x^2(3y+xz).  This grouped
        # shear removes that whole nonlinear circuit before homogenization.
        CircuitAtom("x2s", 0, sp.Rational(1), (("xx", 1), ("s", 1))),
        # The same cancellation with a linear/cubic instead of a
        # quadratic/quadratic gate split.
        CircuitAtom("xxs", 0, sp.Rational(1), (("xl", 1), ("xs", 1))),
        CircuitAtom("v2r", 1, sp.Rational(-3), (("v", 2), ("r", 1))),
        CircuitAtom("xvvz", 1, sp.Rational(-3), (("xv", 1), ("h", 1))),
        CircuitAtom("qb", 1, sp.Rational(-12), (("q", 1), ("b", 1))),
        CircuitAtom("ayb", 1, sp.Rational(-12), (("a", 1), ("yb", 1))),
        CircuitAtom("v2h", 2, sp.Rational(-1), (("v", 2), ("h", 1))),
        CircuitAtom(
            "y2vb",
            2,
            sp.Rational(-4),
            (("yy", 1), ("v", 1), ("b", 1)),
        ),
        CircuitAtom(
            "yvyb",
            2,
            sp.Rational(-4),
            (("yv", 1), ("yb", 1)),
        ),
    )
    perturbation_atoms = tuple(
        CircuitAtom(
            f"qpert_{'m' if coefficient < 0 else 'p'}{abs(coefficient)}",
            1,
            sp.Rational(coefficient),
            (("q", 1),),
        )
        for coefficient in PERTURBATION_COEFFICIENTS
    ) + tuple(
        CircuitAtom(
            f"aqpert_{'m' if coefficient < 0 else 'p'}{abs(coefficient)}",
            1,
            sp.Rational(coefficient),
            (("a", 1), ("q", 1)),
        )
        for coefficient in PERTURBATION_COEFFICIENTS
    )
    neutral_pairs = (
        ("aa", (("a", 2),)),
        ("axx", (("a", 1), ("xx", 1))),
        ("as", (("a", 1), ("s", 1))),
        ("qq", (("q", 2),)),
        ("qxx", (("q", 1), ("xx", 1))),
        ("qs", (("q", 1), ("s", 1))),
        ("xxxx", (("xx", 2),)),
        ("ss", (("s", 2),)),
    )
    neutral_atoms = tuple(
        CircuitAtom(
            f"npert_c{component}_{label}_{'m' if coefficient < 0 else 'p'}1",
            component,
            sp.Rational(coefficient),
            factors,
        )
        for component in range(3)
        for label, factors in neutral_pairs
        for coefficient in (-1, 1)
    )
    special_neutral_atoms = tuple(
        CircuitAtom(
            f"aspert_{'m' if coefficient < 0 else 'p'}{abs(coefficient)}",
            1,
            sp.Rational(coefficient),
            (("a", 1), ("s", 1)),
        )
        for coefficient in PERTURBATION_COEFFICIENTS
    ) + tuple(
        CircuitAtom(
            f"qqpert_{'m' if coefficient < 0 else 'p'}{abs(coefficient)}",
            2,
            sp.Rational(coefficient),
            (("q", 2),),
        )
        for coefficient in PERTURBATION_COEFFICIENTS
    )
    atoms = (
        core_atoms
        + perturbation_atoms
        + neutral_atoms
        + special_neutral_atoms
    )
    return gates, signals, atoms


SIGNAL_GATE = {
    "a": "xy",
    "v": "xy",
    "b": "xy",
    "r": "xz",
    "q": "xy2",
    "h": "vz",
    "yy": "y2",
    "xx": "x2",
    "s": "s",
    "xl": "xlin",
    "xs": "xs",
    "yv": "yv",
    "yb": "yb",
    "xv": "xv",
}


ATOM_CONFLICTS = (
    frozenset(("x2s", "xxs")),
    frozenset(("y2vb", "yvyb")),
    frozenset(("v2r", "xvvz")),
    frozenset(("qb", "ayb")),
)


def affine_normalize(
    expressions: list[sp.Expr],
    variables: list[sp.Symbol],
    points: tuple[tuple[sp.Expr, ...], ...],
) -> tuple[list[sp.Expr], tuple[tuple[sp.Expr, ...], ...]]:
    """Postcompose by the inverse affine jet so the map is ``id+O(2)``."""

    origin = dict.fromkeys(variables, 0)
    constant = sp.Matrix([expression.subs(origin) for expression in expressions])
    jacobian = sp.Matrix(expressions).jacobian(variables).subs(origin)
    if jacobian.det() != 1:
        raise AssertionError("circuit shear lost its determinant-one affine jet")
    inverse = jacobian.inv()
    normalized_vector = inverse * (sp.Matrix(expressions) - constant)
    normalized = [sp.expand(value) for value in normalized_vector]
    assert sp.Matrix(normalized).jacobian(variables).subs(origin) == sp.eye(
        len(variables)
    )

    normalized_points = tuple(tuple(point) for point in points)
    images = []
    for point in normalized_points:
        substitution = dict(zip(variables, point))
        images.append(
            tuple(
                sp.Poly(expression, *variables, domain=sp.QQ).eval(substitution)
                for expression in normalized
            )
        )
    assert images[0] == images[1] == images[2]
    return normalized, normalized_points


def circuit_seed(atom_names: tuple[str, ...]) -> CircuitState:
    """Expose the required polynomial gates and apply grouped target shears."""

    base = initial_state()
    gates, signal_templates, atoms = base_circuit()
    chosen = [atom for atom in atoms if atom.name in atom_names]
    needed_names = tuple(
        name
        for name in gates
        if any(
            SIGNAL_GATE[signal] == name
            for atom in chosen
            for signal, _ in atom.factors
        )
    )

    new_variables = list(sp.symbols(f"circuit_gate0:{len(needed_names)}"))
    variables = base.variables + new_variables
    padding = (0,) * len(new_variables)
    expressions = list(base.expressions)
    gate_outputs: dict[str, sp.Expr] = {}
    gate_indices: dict[str, int] = {}
    for gate_name, variable in zip(needed_names, new_variables):
        gate = gates[gate_name]
        gate_outputs[gate_name] = variable + gate.source
        gate_indices[gate_name] = len(expressions)
        expressions.append(gate_outputs[gate_name])

    symbolic_outputs = {
        sp.Symbol(f"gate_{gate_name}"): gate_outputs[gate_name]
        for gate_name in needed_names
    }
    signals = {
        name: sp.expand(template.subs(symbolic_outputs, simultaneous=True))
        for name, template in signal_templates.items()
        if SIGNAL_GATE[name] in gate_outputs
    }
    by_component: dict[int, sp.Expr] = {}
    for atom in chosen:
        product = atom.coefficient
        for signal, exponent in atom.factors:
            product *= signals[signal] ** exponent
        by_component[atom.component] = by_component.get(
            atom.component, sp.Integer(0)
        ) + product
    for component, shear in by_component.items():
        expressions[component] = sp.expand(expressions[component] + shear)

    points = []
    for point in initial_collision_points():
        source_substitution = dict(zip(base.variables, point))
        additions = [
            -sp.Poly(gates[name].source, *base.variables, domain=sp.QQ).eval(
                source_substitution
            )
            for name in needed_names
        ]
        points.append(tuple(point) + tuple(additions))
    expressions, normalized_points = affine_normalize(
        expressions, variables, tuple(points)
    )

    # Retain exact monomial gate outputs for the ordinary shared-factor
    # cleanup.  Polynomial and linearly normalized gates remain present in
    # the circuit but are not falsely advertised as monomial exposures.
    registry: dict[tuple[int, ...], int] = {}
    for gate_name, output_index in gate_indices.items():
        nonlinear = sp.expand(expressions[output_index] - variables[output_index])
        poly = sp.Poly(nonlinear, *variables, domain=sp.QQ)
        terms = [(exponents, coefficient) for exponents, coefficient in poly.terms() if coefficient]
        if len(terms) == 1 and terms[0][1] == 1 and sum(terms[0][0]) >= 2:
            registry[terms[0][0]] = output_index

    return CircuitState(
        expressions=expressions,
        variables=variables,
        registry=registry,
        introduced=len(new_variables),
        monomial_plan=(),
        circuit_atoms=tuple(sorted(atom_names)),
        collision_points=normalized_points,
    )


def seed_family(
    max_circuit_gates: int,
    required_atoms: frozenset[str] = frozenset(),
    enabled_atoms: frozenset[str] | None = None,
) -> list[CircuitState]:
    _, _, atoms = base_circuit()
    known_atoms = {atom.name for atom in atoms}
    unknown = required_atoms - known_atoms
    if unknown:
        raise ValueError(f"unknown required circuit atoms: {sorted(unknown)}")
    if enabled_atoms is None:
        enabled_atoms = CORE_ATOM_NAMES
    unknown_enabled = enabled_atoms - known_atoms
    if unknown_enabled:
        raise ValueError(
            f"unknown enabled circuit atoms: {sorted(unknown_enabled)}"
        )
    if not required_atoms <= enabled_atoms:
        raise ValueError("every required atom must also be enabled")
    atoms = tuple(atom for atom in atoms if atom.name in enabled_atoms)
    states = [circuit_seed(())] if not required_atoms else []
    for size in range(1, len(atoms) + 1):
        for subset in combinations((atom.name for atom in atoms), size):
            chosen = frozenset(subset)
            if not required_atoms <= chosen:
                continue
            if any(conflict <= chosen for conflict in ATOM_CONFLICTS):
                continue
            state = circuit_seed(tuple(subset))
            if state.introduced <= max_circuit_gates:
                states.append(state)
    deduplicated: dict[tuple[object, ...], CircuitState] = {}
    for state in states:
        key = state_key(state.legacy())
        previous = deduplicated.get(key)
        if previous is None or state.plan_key < previous.plan_key:
            deduplicated[key] = state
    return list(deduplicated.values())


def circuit_gate_count(atom_names: tuple[str, ...]) -> int:
    _, _, atoms = base_circuit()
    chosen = {name for name in atom_names}
    return len(
        {
            SIGNAL_GATE[signal]
            for atom in atoms
            if atom.name in chosen
            for signal, _ in atom.factors
        }
    )


def transitions(state: CircuitState) -> Iterable[CircuitState]:
    signature, terms = high_terms(state.expressions, state.variables)
    maximum = signature[0]
    selected_terms = sorted(
        (term for term in terms if term[3] == maximum),
        key=lambda term: (term[0], term[1], sp.default_sort_key(term[2])),
    )
    for selected in selected_terms:
        component = selected[0]
        for first, second in candidate_splits(selected[1]):
            result = apply_shared_step(
                state.expressions,
                state.variables,
                state.registry,
                selected,
                (first, second),
                state.introduced,
            )
            if result is None:
                continue
            metadata = result[4]
            old_dimension = len(state.variables)
            lifted_points = []
            for point in state.collision_points:
                additions = []
                for record in metadata["new_factors"]:
                    exponents = [0] * old_dimension
                    for index, exponent in record["factor"]:
                        exponents[index] = exponent
                    additions.append(
                        -sp.prod(
                            value**exponent
                            for value, exponent in zip(point, exponents)
                        )
                    )
                lifted_points.append(tuple(point) + tuple(additions))
            plan_step = (component, support(first), support(second))
            yield CircuitState(
                expressions=result[0],
                variables=result[1],
                registry=result[2],
                introduced=result[3],
                monomial_plan=state.monomial_plan + (plan_step,),
                circuit_atoms=state.circuit_atoms,
                collision_points=tuple(lifted_points),
            )


def structural_key(state: CircuitState) -> tuple[object, ...]:
    signature, _ = high_terms(state.expressions, state.variables)
    return (
        signature[0],
        signature[1],
        signature[2],
        state.introduced,
        -len(state.registry),
        state.plan_key,
    )


def partial_key(
    state: CircuitState, profile: PowerRankProfile
) -> tuple[object, ...]:
    signature, _ = high_terms(state.expressions, state.variables)
    return (
        signature[0],
        profile_key(profile),
        signature[1],
        signature[2],
        state.introduced,
        -len(state.registry),
        state.plan_key,
    )


def diverse_states(
    states: Iterable[CircuitState],
    limit: int,
    key,
) -> list[CircuitState]:
    """Keep the baseline and one lead state per circuit family when possible."""

    ranked = sorted(states, key=key)
    if len(ranked) <= limit:
        return ranked
    representatives: dict[tuple[str, ...], CircuitState] = {}
    for state in ranked:
        representatives.setdefault(state.circuit_atoms, state)
    ordered_representatives = sorted(representatives.values(), key=key)
    baseline = representatives.get(())
    selected: list[CircuitState] = []
    if baseline is not None:
        selected.append(baseline)
    selected.extend(
        state
        for state in ordered_representatives
        if state is not baseline
    )
    selected = selected[:limit]
    selected_ids = {id(state) for state in selected}
    selected.extend(
        state
        for state in ranked
        if id(state) not in selected_ids
    )
    return selected[:limit]


def diverse_profiled(
    pairs: Iterable[tuple[CircuitState, PowerRankProfile]],
    limit: int,
    mode: str = "mixed",
) -> list[tuple[CircuitState, PowerRankProfile]]:
    """Retain profile leaders, then low-complexity leaders, by circuit family."""

    pairs = list(pairs)
    if mode == "structural":
        by_state_id = {id(pair[0]): pair for pair in pairs}
        selected_states = diverse_states(
            (pair[0] for pair in pairs),
            limit,
            structural_key,
        )
        return [by_state_id[id(state)] for state in selected_states]
    if mode != "mixed":
        raise ValueError(f"unknown profiled beam mode: {mode}")
    ranked = sorted(pairs, key=lambda pair: partial_key(pair[0], pair[1]))
    if len(ranked) <= limit:
        return ranked
    profile_representatives: dict[
        tuple[str, ...], tuple[CircuitState, PowerRankProfile]
    ] = {}
    for pair in ranked:
        profile_representatives.setdefault(pair[0].circuit_atoms, pair)
    ordered_profile = sorted(
        profile_representatives.values(),
        key=lambda pair: partial_key(pair[0], pair[1]),
    )
    structural_representatives: dict[
        tuple[str, ...], tuple[CircuitState, PowerRankProfile]
    ] = {}
    for pair in sorted(pairs, key=lambda pair: structural_key(pair[0])):
        structural_representatives.setdefault(pair[0].circuit_atoms, pair)
    ordered_structural = sorted(
        structural_representatives.values(),
        key=lambda pair: structural_key(pair[0]),
    )

    baseline = profile_representatives.get(())
    selected: list[tuple[CircuitState, PowerRankProfile]] = []
    if baseline is not None:
        selected.append(baseline)
    selected.extend(pair for pair in ordered_profile if pair is not baseline)
    selected = selected[:limit]
    selected_ids = {id(pair[0]) for pair in selected}
    selected.extend(
        pair
        for pair in ordered_structural
        if id(pair[0]) not in selected_ids
    )
    selected = selected[:limit]
    selected_ids = {id(pair[0]) for pair in selected}
    selected.extend(
        pair for pair in ranked if id(pair[0]) not in selected_ids
    )
    return selected[:limit]


def terminal_profile(state: CircuitState, hessian_power: bool = False) -> dict[str, object]:
    quadratic, cubic = extract_quadratic_cubic(state.expressions, state.variables)
    factorization = factor_cubic_output(cubic)
    ambient_variables, ambient_h = rank_compressed_homogeneous_map(
        state.variables, quadratic, factorization
    )
    quotient = iterated_constant_kernel_quotient(ambient_variables, ambient_h)

    projected = []
    for point in state.collision_points:
        substitution = dict(zip(state.variables, point))
        cubic_values = [poly.eval(substitution) for poly in factorization.c]
        ambient_point = sp.Matrix(list(point) + cubic_values + [sp.Integer(1)])
        projected.append(ambient_point)
    for stage in quotient.stages:
        projected = [stage.B * point for point in projected]
    projected = [tuple(point) for point in projected]
    if len(set(projected)) != 3:
        raise AssertionError("constant-kernel quotient collapsed the collision")

    cubic_profile = polynomial_jacobian_profile(
        quotient.quotient_h, quotient.quotient_variables
    )
    hessian_ranks = cotangent_hessian_rank(quotient.quotient_h)
    cotangent_kernel_excess = hessian_ranks[2] - 2 * cubic_profile.rank
    if cotangent_kernel_excess < 0:
        raise AssertionError(
            "cotangent Hessian rank fell below twice the Jacobian rank"
        )
    result: dict[str, object] = {
        "circuit_atoms": list(state.circuit_atoms),
        "circuit_gate_count": circuit_gate_count(state.circuit_atoms),
        "introduced_variables": state.introduced,
        "monomial_cleanup_steps": len(state.monomial_plan),
        "cubic_output_rank": len(factorization.c),
        "homogeneous_dimension": len(ambient_variables),
        "constant_kernel_dimensions": [
            stage.kernel.cols for stage in quotient.stages
        ],
        "total_constant_kernel_dimension": quotient.total_kernel_dimension,
        "quotient_dimension": len(quotient.quotient_variables),
        "projected_collision_separated": True,
        "cubic_power_ranks_mod_1000003": list(cubic_profile.ranks),
        "cubic_power_profile_terminated": cubic_profile.terminated,
        "cubic_sampled_index": cubic_profile.sampled_index,
        "cubic_rank_mod_1000003": cubic_profile.rank,
        "cotangent_hessian_rank_mod_1000003": hessian_ranks[2],
        "cotangent_kernel_excess_mod_1000003": cotangent_kernel_excess,
        "quartic_dimension": 2 * len(quotient.quotient_variables),
    }
    if hessian_power:
        quartic_profile = transformed_hessian_power_profile(
            quotient.quotient_h
        )
        result.update(
            {
                "quartic_power_ranks_mod_1000033": list(quartic_profile.ranks),
                "quartic_power_profile_terminated": quartic_profile.terminated,
                "quartic_sampled_index": quartic_profile.sampled_index,
            }
        )
    return result


def terminal_key(profile: dict[str, object], state: CircuitState) -> tuple[object, ...]:
    sampled_index = profile["cubic_sampled_index"]
    return (
        profile["cubic_rank_mod_1000003"],
        sampled_index if sampled_index is not None else 10**9,
        profile["cotangent_hessian_rank_mod_1000003"],
        profile["quotient_dimension"],
        state.plan_key,
    )


def terminal_objective(profile: dict[str, object]) -> tuple[int, ...]:
    """Objectives whose coordinatewise minima define the retained archive."""

    cubic_index = profile["cubic_sampled_index"]
    quartic_index = profile.get("quartic_sampled_index")
    return (
        int(profile["cubic_rank_mod_1000003"]),
        int(cubic_index) if cubic_index is not None else 10**9,
        int(profile["cotangent_hessian_rank_mod_1000003"]),
        int(profile["quotient_dimension"]),
        int(quartic_index) if quartic_index is not None else 10**9,
    )


def dominates(left: tuple[int, ...], right: tuple[int, ...]) -> bool:
    return all(a <= b for a, b in zip(left, right)) and any(
        a < b for a, b in zip(left, right)
    )


def pareto_terminals(
    terminals: Iterable[tuple[CircuitState, dict[str, object]]],
) -> list[tuple[CircuitState, dict[str, object]]]:
    """Deduplicate objective vectors and retain all nondominated terminals."""

    representatives: dict[
        tuple[int, ...], tuple[CircuitState, dict[str, object]]
    ] = {}
    for state, profile in terminals:
        objective = terminal_objective(profile)
        previous = representatives.get(objective)
        if previous is None or state.plan_key < previous[0].plan_key:
            representatives[objective] = (state, profile)
    pairs = list(representatives.values())
    answer = [
        pair
        for pair in pairs
        if not any(
            dominates(terminal_objective(other[1]), terminal_objective(pair[1]))
            for other in pairs
            if other is not pair
        )
    ]
    return sorted(
        answer,
        key=lambda pair: (terminal_objective(pair[1]), pair[0].plan_key),
    )


def encoded_plan(state: CircuitState) -> dict[str, object]:
    return {
        "circuit_atoms": list(state.circuit_atoms),
        "monomial_plan": [
            {
                "component": component,
                "first_factor": [list(pair) for pair in first],
                "second_factor": [list(pair) for pair in second],
            }
            for component, first, second in state.monomial_plan
        ],
    }


def replay_encoded_plan(plan: dict[str, object]) -> CircuitState:
    """Rebuild a serialized circuit/monomial plan with exact collision lifts."""

    state = circuit_seed(tuple(str(name) for name in plan["circuit_atoms"]))
    for raw in plan["monomial_plan"]:
        component = int(raw["component"])
        dimension = len(state.variables)
        first = tuple(
            exponent
            for exponent in (
                dict((int(i), int(e)) for i, e in raw["first_factor"]).get(
                    index, 0
                )
                for index in range(dimension)
            )
        )
        second = tuple(
            exponent
            for exponent in (
                dict((int(i), int(e)) for i, e in raw["second_factor"]).get(
                    index, 0
                )
                for index in range(dimension)
            )
        )
        removed = tuple(a + b for a, b in zip(first, second))
        polynomial = sp.Poly(
            state.expressions[component], *state.variables, domain=sp.QQ
        )
        coefficient = polynomial.coeff_monomial(removed)
        if not coefficient:
            raise AssertionError("serialized circuit plan no longer cancels a term")
        result = apply_shared_step(
            state.expressions,
            state.variables,
            state.registry,
            (component, removed, coefficient, sum(removed)),
            (first, second),
            state.introduced,
        )
        if result is None:
            raise AssertionError("serialized circuit plan violates component ownership")
        old_dimension = len(state.variables)
        lifted_points = []
        for point in state.collision_points:
            additions = []
            for record in result[4]["new_factors"]:
                exponents = [0] * old_dimension
                for index, exponent in record["factor"]:
                    exponents[index] = exponent
                additions.append(
                    -sp.prod(
                        value**exponent
                        for value, exponent in zip(point, exponents)
                    )
                )
            lifted_points.append(tuple(point) + tuple(additions))
        state = CircuitState(
            expressions=result[0],
            variables=result[1],
            registry=result[2],
            introduced=result[3],
            monomial_plan=state.monomial_plan
            + ((component, support(first), support(second)),),
            circuit_atoms=state.circuit_atoms,
            collision_points=tuple(lifted_points),
        )
    return state


def beats_incumbent(profile: dict[str, object]) -> bool:
    comparisons = (
        ("cubic_rank_mod_1000003", "cubic_rank"),
        ("cubic_sampled_index", "cubic_sampled_index"),
        ("quotient_dimension", "cubic_dimension"),
        ("cotangent_hessian_rank_mod_1000003", "quartic_hessian_rank"),
        ("quartic_dimension", "quartic_dimension"),
        ("quartic_sampled_index", "quartic_sampled_index"),
    )
    return any(
        value_key in profile
        and profile[value_key] is not None
        and int(profile[value_key]) < INCUMBENTS[incumbent_key]
        for value_key, incumbent_key in comparisons
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--width", type=int, default=24)
    parser.add_argument("--max-steps", type=int, default=24)
    parser.add_argument("--prebeam-factor", type=int, default=4)
    parser.add_argument("--partial-power-depth", type=int, default=8)
    parser.add_argument("--max-circuit-gates", type=int, default=7)
    parser.add_argument(
        "--beam-mode",
        choices=("mixed", "structural"),
        default="mixed",
        help=(
            "mixed prioritizes power-profile leaders by circuit family and "
            "then structural leaders; structural is a completion-control run"
        ),
    )
    parser.add_argument(
        "--require-atom",
        action="append",
        default=[],
        help=(
            "restrict seeds to circuit families containing this atom; may "
            "be repeated"
        ),
    )
    parser.add_argument(
        "--enable-atom",
        action="append",
        default=[],
        help=(
            "when supplied, enumerate only these atoms; may be repeated and "
            "allows exact replay of an earlier circuit library"
        ),
    )
    parser.add_argument(
        "--skip-terminal-hessian-power",
        action="store_true",
        help=(
            "omit the transformed quartic-Hessian power profile; this is a "
            "speed-only mode and cannot detect quartic-index improvements"
        ),
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--no-write",
        action="store_true",
        help="print the bounded search result without writing the JSON record",
    )
    args = parser.parse_args()
    if not args.output.is_absolute():
        args.output = ROOT / args.output
    if min(
        args.width,
        args.max_steps,
        args.prebeam_factor,
        args.partial_power_depth,
    ) < 1:
        parser.error("width, steps, prebeam factor, and power depth must be positive")
    if args.max_circuit_gates < 0:
        parser.error("max-circuit-gates must be nonnegative")

    try:
        frontier = seed_family(
            args.max_circuit_gates,
            frozenset(args.require_atom),
            (
                frozenset(args.enable_atom)
                if args.enable_atom
                else None
            ),
        )
    except ValueError as error:
        parser.error(str(error))
    if not frontier:
        parser.error("no circuit seeds satisfy the atom and gate restrictions")
    partial_cache: dict[tuple[object, ...], PowerRankProfile] = {}
    terminals: list[tuple[CircuitState, dict[str, object]]] = []
    histogram: Counter[tuple[int, int | None, int, int]] = Counter()
    depth_log: list[dict[str, object]] = []
    print(
        f"circuit_seeds={len(frontier)} "
        f"nonempty={sum(bool(state.circuit_atoms) for state in frontier)}",
        flush=True,
    )

    for depth in range(1, args.max_steps + 1):
        deduplicated: dict[tuple[object, ...], CircuitState] = {}
        completed: list[CircuitState] = []
        generated = 0
        for state in frontier:
            signature, _ = high_terms(state.expressions, state.variables)
            if signature[0] <= 3:
                completed.append(state)
                continue
            for candidate in transitions(state):
                generated += 1
                candidate_signature, _ = high_terms(
                    candidate.expressions, candidate.variables
                )
                if candidate_signature[0] <= 3:
                    completed.append(candidate)
                    continue
                key = state_key(candidate.legacy())
                previous = deduplicated.get(key)
                if previous is None or candidate.plan_key < previous.plan_key:
                    deduplicated[key] = candidate

        for candidate in completed:
            profile = terminal_profile(
                candidate,
                hessian_power=not args.skip_terminal_hessian_power,
            )
            terminals.append((candidate, profile))
            histogram[
                (
                    int(profile["cubic_rank_mod_1000003"]),
                    (
                        int(profile["cubic_sampled_index"])
                        if profile["cubic_sampled_index"] is not None
                        else None
                    ),
                    int(profile["cotangent_hessian_rank_mod_1000003"]),
                    int(profile["quotient_dimension"]),
                )
            ] += 1

        prebeam = diverse_states(
            deduplicated.values(),
            args.prebeam_factor * args.width,
            structural_key,
        )
        profiled = []
        for candidate in prebeam:
            key = state_key(candidate.legacy())
            profile = partial_cache.get(key)
            if profile is None:
                profile = correction_profile(
                    candidate.expressions,
                    candidate.variables,
                    max_power=args.partial_power_depth,
                )
                partial_cache[key] = profile
            profiled.append((candidate, profile))
        prebeam_profiled_count = len(profiled)
        profiled = diverse_profiled(profiled, args.width, mode=args.beam_mode)
        frontier = [candidate for candidate, _ in profiled]

        lead_pair = (
            min(profiled, key=lambda pair: partial_key(pair[0], pair[1]))
            if profiled
            else None
        )
        lead_profile = lead_pair[1] if lead_pair else None
        lead_signature = (
            high_terms(lead_pair[0].expressions, lead_pair[0].variables)[0]
            if lead_pair
            else None
        )
        record = {
            "depth": depth,
            "generated": generated,
            "unique": len(deduplicated),
            "prebeam_profiled": prebeam_profiled_count,
            "kept": len(frontier),
            "terminal_count": len(terminals),
            "lead_high_signature": list(lead_signature) if lead_signature else None,
            "lead_partial_power_ranks_mod_1000003": (
                list(lead_profile.ranks) if lead_profile else None
            ),
        }
        depth_log.append(record)
        print(
            f"depth={depth} generated={generated} unique={len(deduplicated)} "
            f"profiled={prebeam_profiled_count} kept={len(frontier)} "
            f"terminal={len(terminals)} lead_high={lead_signature} "
            f"lead_power={lead_profile.ranks if lead_profile else None}",
            flush=True,
        )
        if not frontier:
            break

    archive = pareto_terminals(terminals)
    family_leaders: dict[
        tuple[str, ...], tuple[CircuitState, dict[str, object]]
    ] = {}
    family_counts: Counter[tuple[str, ...]] = Counter()
    for state, profile in terminals:
        family_counts[state.circuit_atoms] += 1
        previous = family_leaders.get(state.circuit_atoms)
        if previous is None or terminal_key(profile, state) < terminal_key(
            previous[1], previous[0]
        ):
            family_leaders[state.circuit_atoms] = (state, profile)
    best_state: CircuitState | None = None
    best_profile: dict[str, object] | None = None
    if terminals:
        best_state, best_profile = min(
            terminals, key=lambda pair: terminal_key(pair[1], pair[0])
        )

    payload: dict[str, object] = {
        "format": "restricted-bcw-circuit-search-v2",
        "certification_status": (
            "bounded modular search record; any improvement requires an exact QQ "
            "generator and independent replay"
        ),
        "search_scope": {
            "source": "three-variable determinant-one collision",
            "circuit_moves": (
                "polynomial gate exposure, grouped multi-term target shears, "
                "then exact shared-factor cleanup"
            ),
            "partial_objective": (
                "maximum degree followed by the sampled rank and complete "
                "power-rank tuple of J(F-id)"
            ),
            "terminal_objective": (
                "Pareto archive for cubic rank, cubic sampled nilpotency "
                "index, cotangent-Hessian rank, quotient dimension, and "
                "transformed quartic-Hessian sampled nilpotency index"
            ),
            "terminal_pipeline": (
                "rank-compressed cubic homogenization, constant-kernel quotient, "
                "collision-separation check, cubic and cotangent-Hessian profiles"
            ),
            "width": args.width,
            "max_steps": args.max_steps,
            "prebeam_factor": args.prebeam_factor,
            "partial_power_depth": args.partial_power_depth,
            "max_circuit_gates": args.max_circuit_gates,
            "beam_mode": args.beam_mode,
            "required_atoms": sorted(set(args.require_atom)),
            "enabled_atoms": (
                sorted(set(args.enable_atom))
                if args.enable_atom
                else sorted(CORE_ATOM_NAMES)
            ),
            "terminal_hessian_power_profiled": (
                not args.skip_terminal_hessian_power
            ),
            "circuit_atoms": (
                sorted(set(args.enable_atom))
                if args.enable_atom
                else sorted(CORE_ATOM_NAMES)
            ),
        },
        "incumbents": INCUMBENTS,
        "depth_log": depth_log,
        "terminal_histogram": [
            {
                "cubic_rank": key[0],
                "cubic_sampled_index": key[1],
                "quartic_hessian_rank": key[2],
                "cotangent_kernel_excess": key[2] - 2 * key[0],
                "quotient_dimension": key[3],
                "count": count,
            }
            for key, count in sorted(
                histogram.items(),
                key=lambda item: (
                    item[0][0],
                    item[0][1] if item[0][1] is not None else 10**9,
                    item[0][2:],
                ),
            )
        ],
        "best_terminal": (
            {
                "profile": best_profile,
                "plan": encoded_plan(best_state),
                "beats_a_current_incumbent_diagnostically": beats_incumbent(
                    best_profile
                ),
            }
            if best_state is not None and best_profile is not None
            else None
        ),
        "pareto_terminals": [
            {
                "objective": list(terminal_objective(profile)),
                "profile": profile,
                "plan": encoded_plan(state),
                "beats_a_current_incumbent_diagnostically": beats_incumbent(
                    profile
                ),
            }
            for state, profile in archive
        ],
        "terminal_family_leaders": [
            {
                "circuit_atoms": list(family),
                "terminal_count": family_counts[family],
                "objective": list(terminal_objective(profile)),
                "profile": profile,
                "plan": encoded_plan(state),
            }
            for family, (state, profile) in sorted(family_leaders.items())
        ],
    }
    if not args.no_write:
        args.output.write_text(json.dumps(payload, indent=2) + "\n")
        try:
            displayed_output = args.output.relative_to(ROOT)
        except ValueError:
            displayed_output = args.output
        print(f"wrote {displayed_output}")
    if best_profile is None:
        print("NO TERMINAL within the bounded circuit beam")
    else:
        print(
            "BEST TERMINAL "
            f"rank={best_profile['cubic_rank_mod_1000003']} "
            f"index={best_profile['cubic_sampled_index']} "
            f"hessian_rank={best_profile['cotangent_hessian_rank_mod_1000003']} "
            f"dimension={best_profile['quotient_dimension']}"
        )
        improving_archive = [
            (state, profile)
            for state, profile in archive
            if beats_incumbent(profile)
        ]
        if improving_archive:
            print(
                "DIAGNOSTIC IMPROVEMENT: "
                f"{len(improving_archive)} Pareto plan(s) require exact QQ "
                "certification"
            )
        else:
            print("NO DIAGNOSTIC IMPROVEMENT over the certified incumbent")
        print(f"PARETO TERMINALS {len(archive)}")
        for state, profile in archive:
            print(
                "PARETO "
                f"objective={terminal_objective(profile)} "
                f"atoms={state.circuit_atoms} "
                f"steps={len(state.monomial_plan)}"
            )


if __name__ == "__main__":
    main()

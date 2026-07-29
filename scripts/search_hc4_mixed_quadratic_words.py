#!/usr/bin/env python3
"""Length-two mixed canonical-word searches for HC(4).

This extends ``search_hc4_mixed_canonical_pivots.py`` from one Hamiltonian
shear to ordered noncommuting words.  It contains three finite families:

* two signed quadratic letters;
* one unit-time quadratic and one unit-time cubic in both orders;
* the fixed-order signed word ``T_H2 o T_H1`` with quadratic H1 and cubic
  H2, classified by its factored Poisson bracket before expansion.

Every letter has the form ``tau*(q_i + epsilon*p_j)^d`` and mixes source and
dual variables.  Commuting words and duplicate polynomial maps are removed.

The quadratic--quadratic maps are linear symplectic and preserve the parent
constant Hessian determinant by congruence.  A word containing a cubic
letter is nonlinear, so its parent determinant is audited separately.
The general searches use the same staged gates as the one-letter search:

1. complete-support affine pivot detection;
2. D(mu+lambda*A,w) in the scalar repair box;
3. the simultaneous rank-at-most-two budget;
4. every complete descended determinant in the small two-pivot repair box.

The fixed-order signed family is pair-focused: only its words with a
two-dimensional coordinate-affine block proceed to the reduced-rank and
complete-determinant gates.

Unequal values or excess rank modulo the good prime 1000003 are exact
characteristic-zero rejection witnesses.  Agreement on the deterministic
point set is reported only as a modular survivor.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from dataclasses import dataclass
from itertools import combinations, product
from pathlib import Path

import sympy as sp

import search_hc4_mixed_canonical_pivots as base


@dataclass(frozen=True)
class Letter:
    letter_id: str
    hamiltonian: sp.Expr
    flow: tuple[sp.Expr, ...]
    metadata: dict[str, int]


@dataclass(frozen=True)
class Word:
    word_id: str
    first: Letter
    second: Letter
    flow: tuple[sp.Expr, ...]


def flow_of(hamiltonian: sp.Expr) -> tuple[sp.Expr, ...]:
    return tuple(
        [
            sp.expand(q_i + sp.diff(hamiltonian, p_i))
            for q_i, p_i in zip(base.q, base.p, strict=True)
        ]
        + [
            sp.expand(p_i - sp.diff(hamiltonian, q_i))
            for q_i, p_i in zip(base.q, base.p, strict=True)
        ]
    )


def letters(
    degree: int = 2, tau_box: tuple[int, ...] = (-1, 1)
) -> tuple[Letter, ...]:
    rows: list[Letter] = []
    for source, dual, epsilon, tau in product(
        range(3), range(3), (-1, 1), tau_box
    ):
        linear = base.q[source] + epsilon * base.p[dual]
        hamiltonian = tau * linear**degree
        rows.append(
            Letter(
                letter_id=(
                    f"d{degree}-q{source}-p{dual}"
                    f"-e{epsilon:+d}-t{tau:+d}"
                ),
                hamiltonian=hamiltonian,
                flow=flow_of(hamiltonian),
                metadata={
                    "source": source,
                    "dual": dual,
                    "epsilon": epsilon,
                    "tau": tau,
                    "degree": degree,
                },
            )
        )
    assert len(rows) == 18 * len(tau_box)
    return tuple(rows)


def linear_map_key(flow: tuple[sp.Expr, ...]) -> tuple[tuple[int, ...], ...]:
    return tuple(
        tuple(
            int(sp.Poly(expression, *base.variables).coeff_monomial(variable))
            for variable in base.variables
        )
        for expression in flow
    )


def polynomial_map_key(
    flow: tuple[sp.Expr, ...],
) -> tuple[tuple[tuple[tuple[int, ...], int], ...], ...]:
    return tuple(
        tuple(
            (monomial, int(coefficient))
            for monomial, coefficient in sp.Poly(
                expression, *base.variables
            ).terms()
        )
        for expression in flow
    )


def is_pure_source_cotangent(flow: tuple[sp.Expr, ...]) -> bool:
    """A linear cotangent lift has no dual variable in its source outputs."""

    return not any(expression.has(*base.p) for expression in flow[:3])


def compose_word_sets(
    alphabet_pairs: tuple[
        tuple[tuple[Letter, ...], tuple[Letter, ...]], ...
    ],
    *,
    expected_noncommuting: int,
    expected_pure_source: int,
    expected_unique: int,
) -> tuple[Word, ...]:
    unique: dict[object, Word] = {}
    noncommuting = 0
    pure_source = 0
    for first_alphabet, second_alphabet in alphabet_pairs:
        for first, second in product(first_alphabet, second_alphabet):
            if base.poisson_bracket(
                first.hamiltonian, second.hamiltonian
            ) == 0:
                continue
            noncommuting += 1

            # Pullback order is Phi o T_first o T_second.
            second_substitution = dict(
                zip(base.variables, second.flow, strict=True)
            )
            flow = tuple(
                sp.expand(
                    expression.subs(
                        second_substitution, simultaneous=True
                    )
                )
                for expression in first.flow
            )
            if is_pure_source_cotangent(flow):
                pure_source += 1
                continue
            key = polynomial_map_key(flow)
            unique.setdefault(
                key,
                Word(
                    word_id=f"{first.letter_id}__{second.letter_id}",
                    first=first,
                    second=second,
                    flow=flow,
                ),
            )

    result = tuple(unique.values())
    assert noncommuting == expected_noncommuting
    assert pure_source == expected_pure_source
    assert len(result) == expected_unique
    return result


def words() -> tuple[Word, ...]:
    quadratic = letters()
    return compose_word_sets(
        ((quadratic, quadratic),),
        expected_noncommuting=648,
        expected_pure_source=48,
        expected_unique=600,
    )


def quadratic_cubic_words() -> tuple[Word, ...]:
    quadratic = letters(degree=2, tau_box=(1,))
    cubic = letters(degree=3, tau_box=(1,))
    return compose_word_sets(
        ((quadratic, cubic), (cubic, quadratic)),
        expected_noncommuting=324,
        expected_pure_source=0,
        expected_unique=324,
    )


def signed_quadratic_cubic_words() -> tuple[Word, ...]:
    """Return the signed words ``T_H2 o T_H1`` with degrees ``(2, 3)``.

    ``compose_word_sets`` stores the outer letter first because its flow is
    the expression into which the inner flow is substituted.  Thus the
    cubic alphabet is first here, while ``word.second`` is the mathematical
    first letter H1 and ``word.first`` is H2.
    """

    quadratic = letters(degree=2, tau_box=(-1, 1))
    cubic = letters(degree=3, tau_box=(-1, 1))
    return compose_word_sets(
        ((cubic, quadratic),),
        expected_noncommuting=648,
        expected_pure_source=0,
        expected_unique=648,
    )


def bracket_signature(h1: Letter, h2: Letter) -> dict[str, object]:
    """Classify ``{H1,H2}`` without expanding either Hamiltonian."""

    first_hits_second = h1.metadata["source"] == h2.metadata["dual"]
    second_hits_first = h1.metadata["dual"] == h2.metadata["source"]
    kappa = (
        h2.metadata["epsilon"] * int(first_hits_second)
        - h1.metadata["epsilon"] * int(second_hits_first)
    )
    assert kappa
    if first_hits_second and second_hits_first:
        incidence = "reciprocal"
    elif first_hits_second:
        incidence = "h1_source_hits_h2_dual"
    else:
        assert second_hits_first
        incidence = "h1_dual_hits_h2_source"

    coefficient = (
        h1.metadata["degree"]
        * h2.metadata["degree"]
        * h1.metadata["tau"]
        * h2.metadata["tau"]
        * kappa
    )
    return {
        "incidence": incidence,
        "linear_pairing": kappa,
        "bracket_coefficient": coefficient,
        "factorization": (
            f"{coefficient}*L1^{h1.metadata['degree'] - 1}"
            f"*L2^{h2.metadata['degree'] - 1}"
        ),
    }


def coordinate_affine_blocks(
    polynomial: sp.Poly,
) -> tuple[int, tuple[tuple[int, ...], ...]]:
    """Return the largest jointly affine coordinate blocks.

    This is deliberately a coordinate-subspace census.  It does not claim
    to classify oblique constant affine directions.
    """

    affine = base.affine_coordinates(polynomial)
    for size in range(len(affine), 0, -1):
        blocks = tuple(
            block
            for block in combinations(affine, size)
            if all(
                sum(bool(monomial[index]) for index in block) <= 1
                for monomial, _ in polynomial.terms()
            )
        )
        if blocks:
            return size, blocks
    return 0, ()


def reduced_rank_audit(
    evaluations: tuple[
        tuple[
            tuple[int, list[int], list[list[int]]],
            tuple[int, list[int], list[list[int]]],
            tuple[int, list[int], list[list[int]]],
        ],
        ...,
    ],
) -> dict[str, object]:
    """Give exact modular lower-rank witnesses for a reduced pencil."""

    maximum = -1
    maximum_witness: dict[str, object] | None = None
    budget_witness: dict[str, object] | None = None
    for point_index, point_evaluations in enumerate(evaluations):
        _, _, h1 = point_evaluations[0]
        _, _, h2 = point_evaluations[1]
        _, _, hb = point_evaluations[2]
        parameters = (
            (0, 0),
            (1, 0),
            (0, 1),
            (1, 1),
            (point_index + 2, 3 * point_index + 1),
        )
        for s1, s2 in parameters:
            reduced = base.matrix_add(
                hb,
                base.matrix_scale(s1, h1),
                base.matrix_scale(s2, h2),
            )
            rank = base.rank_mod(reduced)
            witness = {
                "point": point_index,
                "s": [s1, s2],
                "rank": rank,
            }
            if rank > maximum:
                maximum = rank
                maximum_witness = witness
            if rank > 2 and budget_witness is None:
                budget_witness = witness

    assert maximum_witness is not None
    result: dict[str, object] = {
        "generic_rank_lower_bound": maximum,
        "corank_upper_bound": 4 - maximum,
        "maximum_rank_witness_mod_p": maximum_witness,
    }
    if maximum == 4:
        # A nonzero specialized determinant proves that the generic rank is
        # exactly four over Q, not merely at least four modulo the prime.
        result["generic_rank"] = 4
        result["generic_corank"] = 0
    if budget_witness is not None:
        result["rank_at_least_three_witness_mod_p"] = budget_witness
    return result


def transformed_potential(word: Word) -> sp.Expr:
    substitution = dict(zip(base.variables, word.flow, strict=True))
    return sp.expand(
        base.base_potential.subs(substitution, simultaneous=True)
    )


def split_scalar_poly(
    polynomial: sp.Poly, pivot_index: int
) -> tuple[sp.Poly, sp.Poly]:
    retained = tuple(
        variable
        for index, variable in enumerate(base.variables)
        if index != pivot_index
    )
    a_terms: dict[tuple[int, ...], int] = {}
    b_terms: dict[tuple[int, ...], int] = {}
    for monomial, coefficient in polynomial.terms():
        reduced = monomial[:pivot_index] + monomial[pivot_index + 1 :]
        if monomial[pivot_index] == 1:
            a_terms[reduced] = int(coefficient)
        else:
            assert monomial[pivot_index] == 0
            b_terms[reduced] = int(coefficient)
    return (
        sp.Poly.from_dict(a_terms, retained, domain=sp.ZZ),
        sp.Poly.from_dict(b_terms, retained, domain=sp.ZZ),
    )


def split_pair_poly(
    polynomial: sp.Poly, pair: tuple[int, int]
) -> tuple[sp.Poly, sp.Poly, sp.Poly]:
    retained_indices = tuple(
        index for index in range(6) if index not in pair
    )
    retained = tuple(base.variables[index] for index in retained_indices)
    left_terms: dict[tuple[int, ...], int] = {}
    right_terms: dict[tuple[int, ...], int] = {}
    base_terms: dict[tuple[int, ...], int] = {}
    for monomial, coefficient in polynomial.terms():
        reduced = tuple(monomial[index] for index in retained_indices)
        left_power = monomial[pair[0]]
        right_power = monomial[pair[1]]
        if left_power:
            assert left_power == 1 and right_power == 0
            left_terms[reduced] = int(coefficient)
        elif right_power:
            assert right_power == 1
            right_terms[reduced] = int(coefficient)
        else:
            base_terms[reduced] = int(coefficient)
    return (
        sp.Poly.from_dict(left_terms, retained, domain=sp.ZZ),
        sp.Poly.from_dict(right_terms, retained, domain=sp.ZZ),
        sp.Poly.from_dict(base_terms, retained, domain=sp.ZZ),
    )


def run_search(family: str = "quadratic-quadratic") -> dict[str, object]:
    if family == "quadratic-quadratic":
        all_words = words()
        summary = {
            "raw_letters": 36,
            "raw_ordered_words": 36**2,
            "noncommuting_ordered_words": 648,
            "pure_source_cotangent_words_removed": 48,
            "unique_mixed_words": len(all_words),
        }
        parent_constant_by_congruence = True
    elif family == "quadratic-cubic":
        all_words = quadratic_cubic_words()
        summary = {
            "quadratic_letters": 18,
            "cubic_letters": 18,
            "raw_ordered_words": 2 * 18**2,
            "noncommuting_ordered_words": 324,
            "pure_source_cotangent_words_removed": 0,
            "unique_mixed_words": len(all_words),
        }
        parent_constant_by_congruence = False
    else:
        raise ValueError(f"unknown word family: {family}")

    # Seven independent points are enough for fast exact rejections: a
    # single unequal pair proves nonconstancy.  Agreement remains labelled
    # only as a modular survivor.
    scalar_points = base.deterministic_points(5, count=7)
    pair_points = base.deterministic_points(4, count=7)
    parent_points = base.deterministic_points(6, count=7)

    summary |= {
        "charts_with_affine_pivot": 0,
        "parent_constant_hessian": 0,
        "parent_nonconstant_hessian": 0,
        "parent_hessian_unresolved": 0,
        "scalar_affine_pivots": 0,
        "scalar_remainder_trials": 0,
        "scalar_remainder_exact_survivors": 0,
        "scalar_remainder_modular_survivors": 0,
        "simultaneous_affine_pairs": 0,
        "corank_budget_modular_survivors": 0,
        "complete_determinant_trials": 0,
        "complete_determinant_modular_survivors": 0,
    }
    chart_rows: list[dict[str, object]] = []

    for word_index, word in enumerate(all_words, start=1):
        if word_index == 1 or word_index % 25 == 0:
            print(
                f"progress={word_index}/{len(all_words)} "
                f"word={word.word_id}",
                flush=True,
            )

        transformed = transformed_potential(word)
        polynomial = sp.Poly(
            transformed, *base.variables, domain=sp.ZZ
        )
        affine = base.affine_coordinates(polynomial)
        pairs = base.affine_pairs(polynomial, affine)
        if affine:
            summary["charts_with_affine_pivot"] += 1
        summary["scalar_affine_pivots"] += len(affine)
        summary["simultaneous_affine_pairs"] += len(pairs)

        row: dict[str, object] = {
            "word_id": word.word_id,
            "first": word.first.metadata,
            "second": word.second.metadata,
            "map_coordinate_degrees": [
                sp.Poly(
                    expression, *base.variables
                ).total_degree()
                for expression in word.flow
            ],
            "terms": len(polynomial.terms()),
            "affine_coordinates": [
                str(base.variables[index]) for index in affine
            ],
            "affine_pairs": [
                [str(base.variables[left]), str(base.variables[right])]
                for left, right in pairs
            ],
        }
        if parent_constant_by_congruence:
            row["linear_map"] = linear_map_key(word.flow)

        remainder_survivors: list[dict[str, object]] = []
        remainder_witnesses: list[dict[str, object]] = []
        for pivot_index in affine:
            A, B0 = split_scalar_poly(polynomial, pivot_index)
            jets = (base.Jet.from_poly(A), base.Jet.from_poly(B0))
            evaluations = tuple(
                (jets[0].evaluate(point), jets[1].evaluate(point))
                for point in scalar_points
            )
            for lam, mu in product(
                base.SCALAR_LAMBDAS, base.SCALAR_MUS
            ):
                summary["scalar_remainder_trials"] += 1
                values = base.scalar_remainder_values(
                    evaluations, lam, mu
                )
                witness = base.first_difference(values)
                trial = {
                    "pivot": str(base.variables[pivot_index]),
                    "lambda": lam,
                    "mu": mu,
                }
                if witness is not None:
                    if len(remainder_witnesses) < 3:
                        remainder_witnesses.append(trial | witness)
                    continue
                structural_zero = (
                    family == "quadratic-quadratic"
                    and
                    base.variables[pivot_index] == base.z
                    and not word.first.hamiltonian.has(base.z, base.w)
                    and not word.second.hamiltonian.has(base.z, base.w)
                )
                if structural_zero:
                    summary["scalar_remainder_exact_survivors"] += 1
                    remainder_survivors.append(
                        trial
                        | {
                            "constant": "0",
                            "proof": (
                                "doubled row-count plus constant "
                                "symplectic congruence on (x,y,u,v)"
                            ),
                        }
                    )
                else:
                    summary["scalar_remainder_modular_survivors"] += 1
                    remainder_survivors.append(
                        trial | {"values_mod_p": sorted(set(values))}
                    )
        if remainder_survivors:
            row["scalar_remainder_survivors"] = remainder_survivors
        if remainder_witnesses:
            row["sample_scalar_remainder_witnesses_mod_p"] = (
                remainder_witnesses
            )

        corank_survivors: list[list[str]] = []
        determinant_survivors: list[dict[str, object]] = []
        determinant_witnesses: list[dict[str, object]] = []
        for pair in pairs:
            A1, A2, B0 = split_pair_poly(polynomial, pair)
            jets = (
                base.Jet.from_poly(A1),
                base.Jet.from_poly(A2),
                base.Jet.from_poly(B0),
            )
            evaluations = tuple(
                (
                    jets[0].evaluate(point),
                    jets[1].evaluate(point),
                    jets[2].evaluate(point),
                )
                for point in pair_points
            )

            corank_modular = True
            for point_index, point_evaluations in enumerate(evaluations):
                _, _, h1 = point_evaluations[0]
                _, _, h2 = point_evaluations[1]
                _, _, hb = point_evaluations[2]
                s1 = (point_index + 2) % base.PRIME
                s2 = (3 * point_index + 1) % base.PRIME
                reduced = base.matrix_add(
                    hb,
                    base.matrix_scale(s1, h1),
                    base.matrix_scale(s2, h2),
                )
                if base.rank_mod(reduced) > 2:
                    corank_modular = False
                    break
            if corank_modular:
                summary["corank_budget_modular_survivors"] += 1
                corank_survivors.append(
                    [
                        str(base.variables[pair[0]]),
                        str(base.variables[pair[1]]),
                    ]
                )

            for repair, mu in product(
                base.PAIR_LAMBDAS, base.PAIR_MUS
            ):
                summary["complete_determinant_trials"] += 1
                values = base.pair_repair_values(
                    evaluations, repair, mu
                )
                witness = base.first_difference(values)
                trial = {
                    "pair": [
                        str(base.variables[pair[0]]),
                        str(base.variables[pair[1]]),
                    ],
                    "lambda": repair,
                    "mu": mu,
                }
                if witness is not None:
                    if len(determinant_witnesses) < 3:
                        determinant_witnesses.append(trial | witness)
                    continue
                summary[
                    "complete_determinant_modular_survivors"
                ] += 1
                determinant_survivors.append(
                    trial | {"values_mod_p": sorted(set(values))}
                )

        if corank_survivors:
            row["corank_budget_modular_survivors"] = corank_survivors
        if determinant_survivors:
            row["complete_determinant_modular_survivors"] = (
                determinant_survivors
            )
        if determinant_witnesses:
            row["sample_complete_determinant_witnesses_mod_p"] = (
                determinant_witnesses
            )

        if parent_constant_by_congruence:
            parent_constant: bool | None = True
            parent_witness = None
        else:
            parent_constant, parent_witness = base.audit_parent_hessian(
                transformed,
                parent_points,
                quadratic_generator=False,
            )
        if parent_constant is True:
            summary["parent_constant_hessian"] += 1
            row["parent_constant_hessian"] = True
            row["parent_constant_proof"] = (
                "constant congruence by a determinant-one linear "
                "symplectic map"
            )
        elif parent_constant is False:
            summary["parent_nonconstant_hessian"] += 1
            row["parent_constant_hessian"] = False
            row["parent_hessian_witness_mod_p"] = parent_witness
        else:
            summary["parent_hessian_unresolved"] += 1
            row["parent_constant_hessian"] = (
                "unresolved_after_modular_points"
            )
        chart_rows.append(row)

    return {
        "status": "bounded_search",
        "scope": {
            "base": (
                "collision-centred foundational cubic Keller doubling"
            ),
            "prime": base.PRIME,
            "family": family,
            "letter": (
                "tau*(q_i+epsilon*p_j)^d with the degree/time "
                "box specified by family"
            ),
            "epsilon_box": (-1, 1),
            "scalar_lambda_box": base.SCALAR_LAMBDAS,
            "scalar_mu_box": base.SCALAR_MUS,
            "pair_lambda_box": base.PAIR_LAMBDAS,
            "pair_mu_box": base.PAIR_MUS,
            "commuting_words": "excluded",
            "duplicate_maps": "deduplicated by exact linear matrix",
            "pure_source_transformations": "excluded",
        },
        "summary": summary,
        "charts": chart_rows,
    }


def signed_qc_preclassification() -> tuple[dict[str, object], ...]:
    """Classify the signed degree-(2,3) box before composing its flows."""

    census: Counter[tuple[str, int, int]] = Counter()
    commuting = 0
    for h1, h2 in product(
        letters(degree=2, tau_box=(-1, 1)),
        letters(degree=3, tau_box=(-1, 1)),
    ):
        first_hits_second = (
            h1.metadata["source"] == h2.metadata["dual"]
        )
        second_hits_first = (
            h1.metadata["dual"] == h2.metadata["source"]
        )
        kappa = (
            h2.metadata["epsilon"] * int(first_hits_second)
            - h1.metadata["epsilon"] * int(second_hits_first)
        )
        if not kappa:
            commuting += 1
            continue
        signature = bracket_signature(h1, h2)
        census[
            (
                str(signature["incidence"]),
                int(signature["linear_pairing"]),
                int(signature["bracket_coefficient"]),
            )
        ] += 1

    assert commuting == 648
    assert sum(census.values()) == 648
    return tuple(
        {
            "degree_pair": [2, 3],
            "incidence": incidence,
            "linear_pairing": kappa,
            "bracket_coefficient": coefficient,
            "words": count,
        }
        for (incidence, kappa, coefficient), count in sorted(census.items())
    )


def run_signed_quadratic_cubic_search() -> dict[str, object]:
    """Run the pair-focused signed ``T_H2 o T_H1`` canonical-word box."""

    preclassification = signed_qc_preclassification()
    all_words = signed_quadratic_cubic_words()
    pair_points = base.deterministic_points(4, count=7)
    parent_points = base.deterministic_points(6, count=7)

    summary: dict[str, int] = {
        "quadratic_letters": 36,
        "cubic_letters": 36,
        "raw_ordered_words": 36**2,
        "commuting_words_excluded": 648,
        "noncommuting_ordered_words": 648,
        "duplicate_maps_removed": 0,
        "unique_noncommuting_words": len(all_words),
        "words_with_two_coordinate_affine_pivots": 0,
        "coordinate_affine_pivot_blocks": 0,
        "reduced_pencils_with_rank_at_least_three": 0,
        "reduced_pencils_with_generic_rank_four": 0,
        "parent_constant_hessian": 0,
        "parent_nonconstant_hessian": 0,
        "parent_hessian_unresolved": 0,
        "complete_determinant_trials": 0,
        "complete_determinant_modular_survivors": 0,
    }
    pivot_dimension_census: Counter[int] = Counter()
    classification_census: Counter[
        tuple[str, int, int, int, str]
    ] = Counter()
    chart_rows: list[dict[str, object]] = []

    for word_index, word in enumerate(all_words, start=1):
        if word_index == 1 or word_index % 25 == 0:
            print(
                f"progress={word_index}/{len(all_words)} "
                f"word={word.word_id}",
                flush=True,
            )

        # compose_word_sets stores the outer cubic letter first.
        h1 = word.second
        h2 = word.first
        assert h1.metadata["degree"] == 2
        assert h2.metadata["degree"] == 3
        signature = bracket_signature(h1, h2)

        transformed = transformed_potential(word)
        polynomial = sp.Poly(
            transformed, *base.variables, domain=sp.ZZ
        )
        pivot_dimension, blocks = coordinate_affine_blocks(polynomial)
        pivot_dimension_census[pivot_dimension] += 1
        if pivot_dimension < 2:
            continue

        summary["words_with_two_coordinate_affine_pivots"] += 1
        summary["coordinate_affine_pivot_blocks"] += len(blocks)
        row: dict[str, object] = {
            "word_id": word.word_id,
            "composition": "T_H2 o T_H1",
            "H1": h1.metadata,
            "H2": h2.metadata,
            "poisson_bracket": signature,
            "map_coordinate_degrees": [
                sp.Poly(expression, *base.variables).total_degree()
                for expression in word.flow
            ],
            "transformed_potential_terms": len(polynomial.terms()),
            "coordinate_affine_pivot_dimension": pivot_dimension,
            "maximal_coordinate_affine_blocks": [
                [str(base.variables[index]) for index in block]
                for block in blocks
            ],
        }

        parent_constant, parent_witness = base.audit_parent_hessian(
            transformed,
            parent_points,
            quadratic_generator=False,
        )
        if parent_constant is True:
            parent_class = "constant"
            summary["parent_constant_hessian"] += 1
        elif parent_constant is False:
            parent_class = "nonconstant"
            summary["parent_nonconstant_hessian"] += 1
            row["parent_hessian_witness_mod_p"] = parent_witness
        else:
            parent_class = "unresolved"
            summary["parent_hessian_unresolved"] += 1
        row["parent_hessian"] = parent_class

        block_rows: list[dict[str, object]] = []
        for block in blocks:
            assert len(block) == 2
            A1, A2, B0 = split_pair_poly(polynomial, block)
            jets = (
                base.Jet.from_poly(A1),
                base.Jet.from_poly(A2),
                base.Jet.from_poly(B0),
            )
            evaluations = tuple(
                (
                    jets[0].evaluate(point),
                    jets[1].evaluate(point),
                    jets[2].evaluate(point),
                )
                for point in pair_points
            )
            rank_audit = reduced_rank_audit(evaluations)
            if int(rank_audit["generic_rank_lower_bound"]) >= 3:
                summary[
                    "reduced_pencils_with_rank_at_least_three"
                ] += 1
            if rank_audit.get("generic_rank") == 4:
                summary["reduced_pencils_with_generic_rank_four"] += 1

            determinant_witnesses: list[dict[str, object]] = []
            determinant_survivors: list[dict[str, object]] = []
            for repair, mu in product(
                base.PAIR_LAMBDAS, base.PAIR_MUS
            ):
                summary["complete_determinant_trials"] += 1
                values = base.pair_repair_values(
                    evaluations, repair, mu
                )
                witness = base.first_difference(values)
                trial = {
                    "lambda": repair,
                    "mu": mu,
                }
                if witness is not None:
                    if len(determinant_witnesses) < 3:
                        determinant_witnesses.append(trial | witness)
                    continue
                summary[
                    "complete_determinant_modular_survivors"
                ] += 1
                determinant_survivors.append(
                    trial | {"values_mod_p": sorted(set(values))}
                )

            block_row: dict[str, object] = {
                "pivots": [
                    str(base.variables[index]) for index in block
                ],
                "transformed_reduced_hessian": rank_audit,
                "sample_complete_determinant_witnesses_mod_p": (
                    determinant_witnesses
                ),
            }
            if determinant_survivors:
                block_row[
                    "complete_determinant_modular_survivors"
                ] = determinant_survivors
            block_rows.append(block_row)

            classification_census[
                (
                    str(signature["incidence"]),
                    abs(int(signature["linear_pairing"])),
                    pivot_dimension,
                    int(rank_audit["generic_rank_lower_bound"]),
                    parent_class,
                )
            ] += 1
        row["blocks"] = block_rows
        chart_rows.append(row)

    classification = tuple(
        {
            "degree_pair": [2, 3],
            "poisson_bracket_type": incidence,
            "absolute_linear_pairing": absolute_pairing,
            "coordinate_affine_pivot_dimension": pivot_dimension,
            "transformed_reduced_hessian_generic_rank": rank,
            "transformed_reduced_hessian_corank": 4 - rank,
            "parent_hessian": parent_class,
            "words": count,
        }
        for (
            incidence,
            absolute_pairing,
            pivot_dimension,
            rank,
            parent_class,
        ), count in sorted(classification_census.items())
    )

    assert summary["words_with_two_coordinate_affine_pivots"] == 216
    assert summary["coordinate_affine_pivot_blocks"] == 216
    assert summary["reduced_pencils_with_generic_rank_four"] == 216
    assert summary["parent_nonconstant_hessian"] == 216
    assert summary["complete_determinant_trials"] == 34_992
    assert summary["complete_determinant_modular_survivors"] == 0

    return {
        "status": "bounded_search",
        "scope": {
            "base": (
                "collision-centred foundational cubic Keller doubling"
            ),
            "prime": base.PRIME,
            "composition": "T_H2 o T_H1",
            "degree_pair": [2, 3],
            "H1": "tau1*(q_i+epsilon1*p_j)^2",
            "H2": "tau2*(q_k+epsilon2*p_l)^3",
            "tau_box": [-1, 1],
            "epsilon_box": [-1, 1],
            "commuting_words": "excluded before composition",
            "duplicate_maps": "deduplicated by exact polynomial map",
            "search_gate": (
                "only words with a two-dimensional jointly affine "
                "coordinate block"
            ),
            "affine_subspaces": (
                "coordinate subspaces only; oblique constant directions "
                "are not classified"
            ),
            "pair_lambda_box": base.PAIR_LAMBDAS,
            "pair_mu_box": base.PAIR_MUS,
        },
        "pre_expansion_poisson_bracket_census": preclassification,
        "coordinate_affine_pivot_dimension_census": {
            str(dimension): count
            for dimension, count in sorted(pivot_dimension_census.items())
        },
        "classification": classification,
        "summary": summary,
        "charts": chart_rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        help="write the complete JSON census to this path",
    )
    parser.add_argument(
        "--family",
        choices=(
            "quadratic-quadratic",
            "quadratic-cubic",
            "signed-quadratic-cubic",
        ),
        default="quadratic-quadratic",
    )
    args = parser.parse_args()

    if args.family == "signed-quadratic-cubic":
        result = run_signed_quadratic_cubic_search()
    else:
        result = run_search(args.family)
    encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(encoded)
        digest = hashlib.sha256(encoded.encode()).hexdigest()
        print(f"artifact={args.output}")
        print(f"artifact_sha256={digest}")

    print(f"HC4_MIXED_WORD_SUMMARY family={args.family}")
    print(json.dumps(result["summary"], indent=2, sort_keys=True))
    print(
        f"SCOPE: finite exact length-two mixed word box only "
        f"({args.family})"
    )


if __name__ == "__main__":
    main()

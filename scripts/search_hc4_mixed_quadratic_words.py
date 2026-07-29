#!/usr/bin/env python3
"""Length-two mixed quadratic canonical-word search for HC(4).

This extends ``search_hc4_mixed_canonical_pivots.py`` from one Hamiltonian
shear to ordered words of two noncommuting quadratic shears

    H = tau*(q_i + epsilon*p_j)^2,

with tau,epsilon in {-1,1}.  Every letter mixes source and dual variables.
Commuting words, duplicate symplectic maps, and the exceptional compositions
whose source block is independent of all dual variables are removed.

The surviving maps are linear symplectic, so the collision-centred parent
potential retains its constant Hessian determinant by constant congruence.
Every chart is tested in the same order as the one-letter search:

1. complete-support affine pivot detection;
2. D(mu+lambda*A,w) in the scalar repair box;
3. the simultaneous rank-at-most-two budget;
4. every complete descended determinant in the small two-pivot repair box.

Unequal values or excess rank modulo the good prime 1000003 are exact
characteristic-zero rejection witnesses.  Agreement on the deterministic
point set is reported only as a modular survivor.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from dataclasses import dataclass
from itertools import product
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


def letters() -> tuple[Letter, ...]:
    rows: list[Letter] = []
    for source, dual, epsilon, tau in product(
        range(3), range(3), (-1, 1), (-1, 1)
    ):
        linear = base.q[source] + epsilon * base.p[dual]
        hamiltonian = tau * linear**2
        rows.append(
            Letter(
                letter_id=(
                    f"q{source}-p{dual}-e{epsilon:+d}-t{tau:+d}"
                ),
                hamiltonian=hamiltonian,
                flow=flow_of(hamiltonian),
                metadata={
                    "source": source,
                    "dual": dual,
                    "epsilon": epsilon,
                    "tau": tau,
                },
            )
        )
    assert len(rows) == 36
    return tuple(rows)


def linear_map_key(flow: tuple[sp.Expr, ...]) -> tuple[tuple[int, ...], ...]:
    return tuple(
        tuple(
            int(sp.Poly(expression, *base.variables).coeff_monomial(variable))
            for variable in base.variables
        )
        for expression in flow
    )


def is_pure_source_cotangent(flow: tuple[sp.Expr, ...]) -> bool:
    """A linear cotangent lift has no dual variable in its source outputs."""

    return not any(expression.has(*base.p) for expression in flow[:3])


def words() -> tuple[Word, ...]:
    unique: dict[tuple[tuple[int, ...], ...], Word] = {}
    alphabet = letters()
    for first, second in product(alphabet, repeat=2):
        if base.poisson_bracket(
            first.hamiltonian, second.hamiltonian
        ) == 0:
            continue

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
            continue
        key = linear_map_key(flow)
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
    assert len(result) == 600
    return result


def transformed_potential(word: Word) -> sp.Expr:
    substitution = dict(zip(base.variables, word.flow, strict=True))
    return sp.expand(
        base.base_potential.subs(substitution, simultaneous=True)
    )


def run_search() -> dict[str, object]:
    all_words = words()
    # Seven independent points are enough for fast exact rejections: a
    # single unequal pair proves nonconstancy.  Agreement remains labelled
    # only as a modular survivor.
    scalar_points = base.deterministic_points(5, count=7)
    pair_points = base.deterministic_points(4, count=7)

    summary = {
        "raw_letters": 36,
        "raw_ordered_words": 36**2,
        "noncommuting_ordered_words": 648,
        "pure_source_cotangent_words_removed": 48,
        "unique_mixed_words": len(all_words),
        "charts_with_affine_pivot": 0,
        "scalar_affine_pivots": 0,
        "scalar_remainder_trials": 0,
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
            "linear_map": linear_map_key(word.flow),
            "terms": len(polynomial.terms()),
            "affine_coordinates": [
                str(base.variables[index]) for index in affine
            ],
            "affine_pairs": [
                [str(base.variables[left]), str(base.variables[right])]
                for left, right in pairs
            ],
            "parent_constant_hessian": True,
            "parent_constant_proof": (
                "constant congruence by a determinant-one linear "
                "symplectic map"
            ),
        }

        remainder_survivors: list[dict[str, object]] = []
        remainder_witnesses: list[dict[str, object]] = []
        for pivot_index in affine:
            A, B0, _ = base.split_scalar(transformed, pivot_index)
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
                summary["scalar_remainder_modular_survivors"] += 1
                remainder_survivors.append(
                    trial | {"values_mod_p": sorted(set(values))}
                )
        if remainder_survivors:
            row["scalar_remainder_modular_survivors"] = (
                remainder_survivors
            )
        if remainder_witnesses:
            row["sample_scalar_remainder_witnesses_mod_p"] = (
                remainder_witnesses
            )

        corank_survivors: list[list[str]] = []
        determinant_survivors: list[dict[str, object]] = []
        determinant_witnesses: list[dict[str, object]] = []
        for pair in pairs:
            A1, A2, B0, _ = base.split_pair(transformed, pair)
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
        chart_rows.append(row)

    return {
        "status": "bounded_search",
        "scope": {
            "base": (
                "collision-centred foundational cubic Keller doubling"
            ),
            "prime": base.PRIME,
            "letter": "tau*(q_i+epsilon*p_j)^2",
            "tau_box": (-1, 1),
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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        help="write the complete JSON census to this path",
    )
    args = parser.parse_args()

    result = run_search()
    encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(encoded)
        digest = hashlib.sha256(encoded.encode()).hexdigest()
        print(f"artifact={args.output}")
        print(f"artifact_sha256={digest}")

    print("HC4_MIXED_QUADRATIC_WORD_SUMMARY")
    print(json.dumps(result["summary"], indent=2, sort_keys=True))
    print(
        "SCOPE: finite exact length-two mixed quadratic word box only"
    )


if __name__ == "__main__":
    main()

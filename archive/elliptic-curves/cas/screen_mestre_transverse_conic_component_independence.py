#!/usr/bin/env python3
"""Bounded exact finite-reduction quotient screen for the second component.

This uses the conic-rational parameter ``s`` directly.  A strict increase in
the mod-3 quotient rank at one fibre would disprove a generic relation with
the visible subgroup.  Conversely, a no-gain result is only the stated
bounded search result.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from math import gcd
import json
from pathlib import Path
from typing import Iterable

from mestre_root_tuples import SixRootMestreConstruction
from search_mestre_dsquare_four import rational_square_root
from search_mestre_root_tuple_scale import (
    primitive_visible_points,
    quartic_point_to_jacobian,
    quartic_value,
)
from search_mestre_root_tuple_scale_max200 import mod3_independence_certificate
from verify_mestre_transverse_two_section_conic_component import (
    component_coordinates,
    split_roots,
)


Q = Fraction


def rational_text(value: Fraction) -> str:
    value = Q(value)
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def bounded_rationals(height: int, *, nonzero: bool = False) -> tuple[Fraction, ...]:
    values = set()
    for denominator in range(1, height + 1):
        for numerator in range(-height, height + 1):
            if (nonzero and numerator == 0) or gcd(numerator, denominator) != 1:
                continue
            values.add(Q(numerator, denominator))
    return tuple(sorted(values, key=lambda value: (abs(value.numerator), value.denominator, value.numerator)))


def specialization_record(s_value: Fraction, t_value: Fraction, *, prime_bound: int) -> dict[str, object] | None:
    roots = split_roots(s_value)
    if len(set(roots)) != 6:
        return None
    construction = SixRootMestreConstruction(roots)
    coordinates = component_coordinates(s_value)
    quartic = construction.primitive_quartic_coefficients(t_value)
    visible = tuple(
        quartic_point_to_jacobian(construction, t_value, point)
        for point in primitive_visible_points(construction, t_value)
    )
    affine = []
    for intercept, slope in ((coordinates[4], coordinates[5]), (coordinates[6], coordinates[7])):
        x_value = intercept + slope * t_value
        y_value = rational_square_root(quartic_value(quartic, x_value))
        if y_value in (None, 0):
            return None
        affine.append(quartic_point_to_jacobian(construction, t_value, (x_value, y_value)))
    coefficients = construction.primitive_jacobian_coefficients(t_value)
    visible_certificate = mod3_independence_certificate(coefficients, visible, prime_bound=prime_bound)
    augmented_certificate = mod3_independence_certificate(
        coefficients, (*visible, *affine), prime_bound=prime_bound
    )
    return {
        "s": rational_text(s_value),
        "T": rational_text(t_value),
        "visible_rank": visible_certificate["combined_exact_rank_over_F3"],
        "augmented_rank": augmented_certificate["combined_exact_rank_over_F3"],
        "visible_pivots": visible_certificate["independent_subset_indices_one_based"],
        "augmented_pivots": augmented_certificate["independent_subset_indices_one_based"],
        "certificate_primes": augmented_certificate["certificate_primes"],
    }


def screen(
    *, root_height: int, parameter_height: int, prime_bound: int,
    start: int = 0, count: int | None = None,
) -> dict[str, object]:
    roots = tuple(value for value in bounded_rationals(root_height) if value not in (-7, 7))
    parameters = bounded_rationals(parameter_height, nonzero=True)
    pairs = tuple((s_value, t_value) for s_value in roots for t_value in parameters)
    selected = pairs[start:] if count is None else pairs[start : start + count]
    records = []
    for s_value, t_value in selected:
        try:
            record = specialization_record(s_value, t_value, prime_bound=prime_bound)
        except (AssertionError, ValueError, ZeroDivisionError):
            record = None
        if record is not None:
            records.append(record)
    gains = [record for record in records if record["augmented_rank"] > record["visible_rank"]]
    best = max(
        records,
        default=None,
        key=lambda record: (record["augmented_rank"], record["visible_rank"], record["s"], record["T"]),
    )
    result = {
        "status": "bounded exact second-component independence screen completed",
        "root_height": root_height,
        "mestre_parameter_height": parameter_height,
        "prime_bound": prime_bound,
        "candidate_pair_count": len(pairs),
        "admissible_specialization_count": len(records),
        "strict_quotient_rank_gain_count": len(gains),
        "strict_quotient_rank_gains": gains,
        "best_observed": best,
        "scope": "a bounded finite-reduction screen; no failure is a rank upper bound",
    }
    if start or count is not None:
        result.update({"pair_offset": start, "pair_batch_count": len(selected)})
    return result


def merge_batches(paths: Iterable[Path]) -> dict[str, object]:
    batches = [json.loads(path.read_text()) for path in paths]
    if not batches:
        raise ValueError("at least one bounded result is required")
    reference = batches[0]
    keys = ("root_height", "mestre_parameter_height", "prime_bound", "candidate_pair_count", "scope")
    if any(any(batch[key] != reference[key] for key in keys) for batch in batches[1:]):
        raise ValueError("the bounded results do not have common screen parameters")
    covered = set()
    for batch in batches:
        offset, batch_count = int(batch["pair_offset"]), int(batch["pair_batch_count"])
        interval = set(range(offset, offset + batch_count))
        if covered & interval:
            raise ValueError("the bounded results overlap")
        covered |= interval
    total = int(reference["candidate_pair_count"])
    if covered != set(range(total)):
        raise ValueError("the bounded results do not cover every parameter pair")
    gains = [record for batch in batches for record in batch["strict_quotient_rank_gains"]]
    observed = [batch["best_observed"] for batch in batches if batch["best_observed"] is not None]
    return {
        "status": "bounded exact second-component independence screen completed",
        "root_height": reference["root_height"],
        "mestre_parameter_height": reference["mestre_parameter_height"],
        "prime_bound": reference["prime_bound"],
        "candidate_pair_count": total,
        "admissible_specialization_count": sum(int(batch["admissible_specialization_count"]) for batch in batches),
        "strict_quotient_rank_gain_count": len(gains),
        "strict_quotient_rank_gains": gains,
        "best_observed": max(
            observed, default=None,
            key=lambda record: (record["augmented_rank"], record["visible_rank"], record["s"], record["T"]),
        ),
        "scope": reference["scope"],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root-height", type=int, default=5)
    parser.add_argument("--parameter-height", type=int, default=5)
    parser.add_argument("--prime-bound", type=int, default=101)
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("--count", type=int)
    parser.add_argument("--merge-batches", nargs="+", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = (
        merge_batches(args.merge_batches)
        if args.merge_batches is not None
        else screen(
            root_height=args.root_height, parameter_height=args.parameter_height,
            prime_bound=args.prime_bound, start=args.start, count=args.count,
        )
    )
    rendered = json.dumps(
        result,
        indent=2,
        sort_keys=True,
    ) + "\n"
    if args.output is None:
        print(rendered, end="")
    else:
        args.output.write_text(rendered)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Exact branch-incidence codes for quadratic covers of ``P^1``.

The input is a finite list of square conditions ``z_i^2=g_i(u)`` over
``QQ(u)``.  The module reduces each divisor modulo two, records the resulting
binary incidence row, and enumerates every elementary quadratic quotient of
the associated multiquadratic cover.  A closed point of degree ``d`` is kept
as one Galois-orbit column of weight ``d``; all of its geometric branch
places have the same incidence column.  This compression is lossless for the
branch counts and genera used here.

Constants are deliberately discarded: the analysis is geometric over
``Qbar``.  Consequently a zero incidence word is an unramified, geometrically
split square condition rather than a genus-zero double cover.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
import hashlib
from itertools import combinations
from math import factorial
from typing import Any, Iterable, Mapping, Sequence

import sympy as sp


Q = Fraction


def rational_text(value: Any) -> str:
    """Write a rational value in the stable notation used by JSON artifacts."""

    value = sp.Rational(value)
    return str(value.p) if value.q == 1 else f"{value.p}/{value.q}"


def polynomial_coefficients(polynomial: sp.Poly) -> list[str]:
    """Return ascending, primitive monic-normalized coefficients."""

    return [rational_text(value) for value in reversed(polynomial.all_coeffs())]


def polynomial_key(polynomial: sp.Poly) -> tuple[str, ...]:
    """Identify a QQ-place by its monic irreducible polynomial."""

    monic = polynomial.monic()
    return tuple(polynomial_coefficients(monic))


def sha256_lines(lines: Iterable[str]) -> str:
    digest = hashlib.sha256()
    for line in lines:
        digest.update((line + "\n").encode())
    return digest.hexdigest()


def polynomial_sha256(polynomial: sp.Poly) -> str:
    """Stable digest of a polynomial's ascending rational coefficients."""

    return sha256_lines(polynomial_coefficients(polynomial))


def gf2_rank(rows: Iterable[int]) -> int:
    """Rank of bit vectors represented by non-negative Python integers."""

    pivots: dict[int, int] = {}
    for row in rows:
        while row:
            pivot = row.bit_length() - 1
            previous = pivots.get(pivot)
            if previous is None:
                pivots[pivot] = row
                break
            row ^= previous
    return len(pivots)


def quotient_genus(geometric_branch_count: int) -> int | None:
    """Genus of a connected quadratic quotient of ``P^1``.

    A nontrivial quadratic extension of ``Qbar(u)`` has an even positive
    number of branch places.  Returning ``None`` for zero makes the constant
    squareclass case explicit rather than accidentally calling it genus zero.
    """

    if geometric_branch_count == 0:
        return None
    if geometric_branch_count < 2 or geometric_branch_count % 2:
        raise AssertionError("a quadratic branch divisor must have positive even degree")
    return (geometric_branch_count - 2) // 2


@dataclass(frozen=True)
class SquareCondition:
    """One exact condition ``z^2=numerator/denominator`` in ``QQ(u)``."""

    label: str
    numerator: sp.Poly
    denominator: sp.Poly

    @classmethod
    def polynomial(cls, label: str, polynomial: sp.Poly) -> "SquareCondition":
        return cls(label, polynomial, sp.Poly(1, polynomial.gens[0], domain=sp.QQ))


def conditions_from_json_records(
    records: Sequence[Mapping[str, Any]], *, parameter_name: str = "u"
) -> list[SquareCondition]:
    """Decode exact square conditions from coefficient-list records.

    Each record has a distinct ``label`` and a nonempty
    ``numerator_coefficients_ascending`` list.  The optional denominator list
    has the same convention and defaults to ``[1]``.  Coefficients are parsed
    by :class:`sympy.Rational`, so JSON strings such as ``"-17/35"`` remain
    exact.  This is the interchange format for later ICARM and K3 inputs.
    """

    variable = sp.Symbol(parameter_name)

    def polynomial(values: Any, field: str) -> sp.Poly:
        if not isinstance(values, Sequence) or isinstance(values, (str, bytes)) or not values:
            raise ValueError(f"{field} must be a nonempty coefficient list")
        expression = sum(sp.Rational(value) * variable**index for index, value in enumerate(values))
        result = sp.Poly(expression, variable, domain=sp.QQ)
        if result.is_zero:
            raise ValueError(f"{field} must not be zero")
        return result

    conditions = []
    for record in records:
        try:
            label = record["label"]
            numerator = polynomial(record["numerator_coefficients_ascending"], "numerator_coefficients_ascending")
            denominator = polynomial(record.get("denominator_coefficients_ascending", [1]), "denominator_coefficients_ascending")
        except KeyError as error:
            raise ValueError(f"square-condition record is missing {error.args[0]}") from error
        if not isinstance(label, str) or not label:
            raise ValueError("square-condition label must be a nonempty string")
        conditions.append(SquareCondition(label, numerator, denominator))
    return conditions


def _odd_factor_places(condition: SquareCondition) -> dict[tuple[str, ...], sp.Poly]:
    """Return the finite support of ``div(g) mod 2``.

    Numerator and denominator valuations have the same parity modulo two, so
    their odd irreducible factors are toggled into a single support set.
    """

    answer: dict[tuple[str, ...], sp.Poly] = {}
    for polynomial in (condition.numerator, condition.denominator):
        for factor, exponent in sp.factor_list(polynomial)[1]:
            if exponent % 2:
                key = polynomial_key(factor)
                if key in answer:
                    del answer[key]
                else:
                    answer[key] = factor.monic()
    return answer


def _mask_labels(mask: int, labels: Sequence[str]) -> list[str]:
    return [label for index, label in enumerate(labels) if mask & (1 << index)]


def analyze_square_conditions(
    conditions: Sequence[SquareCondition],
    *,
    maximum_records: int | None = None,
) -> dict[str, Any]:
    """Build the exact weighted binary branch-incidence code.

    ``maximum_records`` may cap only the verbose codeword list; all aggregate
    minima, histograms, dimensions, and the low-genus test are still exact.
    The applications in this repository are small enough to retain every
    codeword.
    """

    if not conditions:
        raise ValueError("at least one square condition is required")
    labels = [condition.label for condition in conditions]
    if len(set(labels)) != len(labels):
        raise ValueError("square-condition labels must be distinct")
    variable = conditions[0].numerator.gens[0]
    if any(
        condition.numerator.gens != (variable,)
        or condition.denominator.gens != (variable,)
        for condition in conditions
    ):
        raise ValueError("all conditions must use one common parameter")

    condition_places = [_odd_factor_places(condition) for condition in conditions]
    places: dict[tuple[str, ...], sp.Poly] = {}
    for factors in condition_places:
        places.update(factors)
    ordered_keys = sorted(places, key=lambda key: (len(key) - 1, key))
    finite_place_index = {key: index for index, key in enumerate(ordered_keys)}

    finite_rows: list[int] = []
    infinity_bits: list[int] = []
    for factors in condition_places:
        row = 0
        finite_degree = 0
        for key, factor in factors.items():
            row |= 1 << finite_place_index[key]
            finite_degree += factor.degree()
        finite_rows.append(row)
        infinity_bits.append(finite_degree % 2)

    use_infinity = any(infinity_bits)
    infinity_index = len(ordered_keys) if use_infinity else None
    rows = [
        row | ((infinity_bit << infinity_index) if use_infinity else 0)
        for row, infinity_bit in zip(finite_rows, infinity_bits, strict=True)
    ]
    place_degrees = [places[key].degree() for key in ordered_keys]
    place_labels = [f"B{index + 1}" for index in range(len(ordered_keys))]
    if use_infinity:
        place_degrees.append(1)
        place_labels.append("infinity")

    condition_records = []
    for condition, factors, row, infinity_bit in zip(
        conditions, condition_places, rows, infinity_bits, strict=True
    ):
        support = [
            place_labels[index]
            for index in range(len(place_labels))
            if row & (1 << index)
        ]
        branch_degree = sum(
            place_degrees[index]
            for index in range(len(place_labels))
            if row & (1 << index)
        )
        condition_records.append(
            {
                "label": condition.label,
                "numerator_sha256": polynomial_sha256(condition.numerator),
                "denominator_sha256": polynomial_sha256(condition.denominator),
                "finite_odd_factor_count": len(factors),
                "infinity_branched": bool(infinity_bit),
                "branch_support": support,
                "geometric_branch_count": branch_degree,
                "quadratic_quotient_genus": quotient_genus(branch_degree),
            }
        )

    # Every elementary quadratic quotient is one nonzero coefficient vector.
    codeword_histogram: Counter[tuple[int, int | None]] = Counter()
    low_genus: list[dict[str, Any]] = []
    cancellations: list[dict[str, Any]] = []
    records: list[dict[str, Any]] = []
    codeword_manifest_lines: list[str] = []
    minimum_branch_count: int | None = None
    minimum_masks: list[int] = []
    for mask in range(1, 1 << len(rows)):
        word = 0
        for index, row in enumerate(rows):
            if mask & (1 << index):
                word ^= row
        branch_count = sum(
            place_degrees[index]
            for index in range(len(place_labels))
            if word & (1 << index)
        )
        genus = quotient_genus(branch_count)
        codeword_histogram[branch_count, genus] += 1
        if minimum_branch_count is None or branch_count < minimum_branch_count:
            minimum_branch_count = branch_count
            minimum_masks = [mask]
        elif branch_count == minimum_branch_count:
            minimum_masks.append(mask)
        record = {
            "combination": _mask_labels(mask, labels),
            "geometric_branch_count": branch_count,
            "quadratic_quotient_genus": genus,
        }
        if branch_count == 0:
            cancellations.append(record)
        elif genus is not None and genus <= 1:
            low_genus.append(record)
        codeword_manifest_lines.append(
            "|".join(
                [
                    ",".join(record["combination"]),
                    str(record["geometric_branch_count"]),
                    str(record["quadratic_quotient_genus"]),
                ]
            )
        )
        if maximum_records is None or len(records) < maximum_records:
            records.append(record)

    # Columns with the same incidence vector cannot be distinguished by any
    # quotient.  First record this on the compressed closed-place code.  Then
    # expand their weights formally over Qbar: an irreducible factor of degree
    # d comprises d geometric branch places, all with the same incidence
    # pattern.  The latter is the coding-theoretic branch-place identification
    # relevant to geometric code automorphisms.
    column_classes: dict[tuple[int, int], list[str]] = defaultdict(list)
    geometric_classes: dict[int, list[tuple[str, int]]] = defaultdict(list)
    for index, label in enumerate(place_labels):
        pattern = sum(
            ((row >> index) & 1) << direction for direction, row in enumerate(rows)
        )
        column_classes[pattern, place_degrees[index]].append(label)
        geometric_classes[pattern].append((label, place_degrees[index]))
    indistinguishable_columns = [
        {
            "incidence_pattern_binary": format(pattern, f"0{len(rows)}b"),
            "geometric_degree_per_place": degree,
            "branch_places": members,
        }
        for (pattern, degree), members in sorted(column_classes.items())
        if len(members) > 1
    ]
    geometric_incidence_classes = [
        {
            "incidence_pattern_binary": format(pattern, f"0{len(rows)}b"),
            "closed_branch_places": [label for label, _ in members],
            "geometric_branch_place_count": sum(degree for _, degree in members),
        }
        for pattern, members in sorted(geometric_classes.items())
        if sum(degree for _, degree in members) > 1
    ]
    geometric_code_automorphism_order = 1
    for members in geometric_classes.values():
        geometric_code_automorphism_order *= factorial(
            sum(degree for _, degree in members)
        )

    pair_records = []
    for left, right in combinations(range(len(rows)), 2):
        shared = rows[left] & rows[right]
        shared_degree = sum(
            place_degrees[index]
            for index in range(len(place_labels))
            if shared & (1 << index)
        )
        if shared_degree:
            pair_records.append(
                {
                    "directions": [labels[left], labels[right]],
                    "shared_geometric_branch_count": shared_degree,
                    "shared_branch_places": [
                        place_labels[index]
                        for index in range(len(place_labels))
                        if shared & (1 << index)
                    ],
                }
            )

    rank = gf2_rank(rows)
    full_union = 0
    for row in rows:
        full_union |= row
    full_branch_count = sum(
        place_degrees[index]
        for index in range(len(place_labels))
        if full_union & (1 << index)
    )
    full_cover_genus = (
        1 + (2 ** (rank - 2)) * (full_branch_count - 4)
        if rank >= 2
        else quotient_genus(full_branch_count)
    )

    return {
        "definition": (
            "Rows encode div(g_i) modulo 2.  Irreducible QQ factors are closed "
            "branch places of geometric weight equal to their degree; infinity is "
            "included when it has odd valuation."
        ),
        "parameter": str(variable),
        "direction_count": len(rows),
        "geometric_branch_place_count": sum(place_degrees),
        "closed_branch_place_count": len(place_labels),
        "code_dimension": rank,
        "conditions": condition_records,
        "branch_places": [
            {
                "label": place_labels[index],
                "geometric_degree": place_degrees[index],
                "polynomial_sha256": (
                    polynomial_sha256(places[ordered_keys[index]])
                    if index < len(ordered_keys)
                    else None
                ),
            }
            for index in range(len(place_labels))
        ],
        "minimum_quadratic_quotient": {
            "geometric_branch_count": minimum_branch_count,
            "genus": quotient_genus(minimum_branch_count or 0),
            "combinations": [_mask_labels(mask, labels) for mask in minimum_masks],
        },
        "quadratic_quotient_histogram": [
            {
                "geometric_branch_count": branch_count,
                "genus": genus,
                "count": count,
            }
            for (branch_count, genus), count in sorted(
                codeword_histogram.items(), key=lambda item: (item[0][0], str(item[0][1]))
            )
        ],
        "all_quadratic_quotients": records,
        "all_quadratic_quotients_recorded": len(records) == (1 << len(rows)) - 1,
        "low_genus_quadratic_quotients": low_genus,
        "complete_branch_cancellations": cancellations,
        "shared_branch_pairs": pair_records,
        "incidence_indistinguishable_closed_branch_places": indistinguishable_columns,
        "geometric_incidence_indistinguishable_branch_places": geometric_incidence_classes,
        "geometric_incidence_code_automorphisms": {
            "coordinate_permutation_group_order": geometric_code_automorphism_order,
            "meaning": (
                "Product of symmetric groups on geometric branch places with "
                "the same binary incidence column.  This is a code automorphism "
                "count only; it does not assert that these permutations are "
                "induced by PGL_2 or by automorphisms of the cover."
            ),
        },
        "multiquadratic_cover": {
            "connected_degree": 2**rank,
            "geometric_branch_count": full_branch_count,
            "genus": full_cover_genus,
            "genus_formula": (
                "1 + 2^(r-2)*(B-4) for a connected (Z/2)^r cover of P1, "
                "where r is code dimension and B is the branch-union size"
            ),
        },
        "codeword_manifest_sha256": sha256_lines(codeword_manifest_lines),
    }

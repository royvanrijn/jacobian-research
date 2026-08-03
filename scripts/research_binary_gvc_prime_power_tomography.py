#!/usr/bin/env python3
"""Prime-power tomography for a finite binary Hall-packet surrogate.

This script studies the two-colour partition configuration

    R_i -> (1, 0, i),    B_j -> (0, 1, j),    0 <= i,j <= r.

A fibre fixes the two colour counts and the total level.  Its primitive
semigroup relations are the Graver basis of this configuration.  The
already-closed circuit range has support at most four, so the pilot packet
envelope consists of mixed-colour Graver relations with support at least
five.  This is the projected-scroll quotient model appearing in Section 7
of BINARY_GVC_UNIFORM_FACE_TERMINATION.md; it is not an enumeration of all
Hall--jet shells.

For a fibre state x at scale N=q*p**e, the script records

* Legendre valuations of both multinomials and the common radial factorial;
* every Kummer carry position on each marked side;
* p-free factorial units modulo p**k;
* low base-p digits of the side partitions; and
* exact group-ring data recovering all C_t character traces.

The v2 output keeps the scalar factorial-weight tomography separate from
the marked-side decorations.  It reports the first valuation or unit
separator, an exact Stirling separator, the first marked/character
separator, and any finite-window collision not explained by exact
scaled-factorial equality.  ``--primitive-only`` omits global repeated-state
grouping but retains every primitive relation, every projected support basis,
first-separator certificates, and full probe rows for surviving collisions.

Normaliz computes the Graver basis exactly through a Lawrence lifting.  A
single full-support computation contains the basis of every coordinate
subconfiguration: filter for relations whose support lies in that
subconfiguration.  The finite prime/exponent window and the projected
two-colour model make the collision census a bounded experiment, not the
proof of unrestricted GVC(2), and a surviving signature is not by itself a
GVC counterexample.  The unrestricted proof is the separate Hall-envelope
theorem.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
import shutil
import subprocess
from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from functools import cache
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any


Vector = tuple[int, ...]


@dataclass(frozen=True, order=True)
class State:
    operator: Vector
    polynomial: Vector


@dataclass(frozen=True, order=True)
class Relation:
    left: State
    right: State


def is_prime(number: int) -> bool:
    if number < 2:
        return False
    divisor = 2
    while divisor * divisor <= number:
        if number % divisor == 0:
            return False
        divisor += 1
    return True


def parse_csv_integers(text: str) -> tuple[int, ...]:
    values = tuple(int(value) for value in text.split(",") if value)
    if not values:
        raise ValueError("expected a nonempty comma-separated integer list")
    return values


def normaliz_version(executable: str) -> str:
    completed = subprocess.run(
        [executable, "--version"],
        check=True,
        capture_output=True,
        text=True,
    )
    return completed.stdout.splitlines()[0].strip()


def normaliz_hilbert_basis(
    equations: tuple[Vector, ...],
    executable: str,
) -> tuple[Vector, ...]:
    """Hilbert basis of the nonnegative integer solutions of equations."""

    if not equations or not equations[0]:
        raise ValueError("Normaliz equations must be a nonempty matrix")
    width = len(equations[0])
    if any(len(row) != width for row in equations):
        raise ValueError("ragged Normaliz equation matrix")

    lines = [f"amb_space {width}", f"inequalities {width}"]
    for row in range(width):
        lines.append(
            " ".join("1" if row == column else "0" for column in range(width))
        )
    lines.append(f"equations {len(equations)}")
    lines.extend(" ".join(str(value) for value in row) for row in equations)

    with TemporaryDirectory(prefix="binary-gvc-tomography-") as directory:
        project = Path(directory) / "basis"
        project.with_suffix(".in").write_text(
            "\n".join(lines) + "\n",
            encoding="utf-8",
        )
        completed = subprocess.run(
            [executable, "-N", "-f", "-x=1", str(project)],
            capture_output=True,
            text=True,
        )
        if completed.returncode:
            raise RuntimeError(
                "Normaliz failed:\n"
                + completed.stdout
                + completed.stderr
            )
        tokens = tuple(
            int(token)
            for token in project.with_suffix(".gen").read_text(
                encoding="utf-8"
            ).split()
        )

    row_count, column_count = tokens[:2]
    if column_count != width:
        raise AssertionError((column_count, width))
    entries = tokens[2 : 2 + row_count * column_count]
    if len(entries) != row_count * column_count:
        raise AssertionError("truncated Normaliz generator matrix")
    return tuple(
        tuple(entries[index * width : (index + 1) * width])
        for index in range(row_count)
    )


def state_from_dense(vector: Vector, radial_degree: int) -> State:
    width = radial_degree + 1
    if len(vector) != 2 * width:
        raise ValueError("wrong dense state width")
    return State(vector[:width], vector[width:])


def orient_relation(left: State, right: State) -> Relation:
    return Relation(left, right) if left <= right else Relation(right, left)


def colored_graver_basis(
    radial_degree: int,
    executable: str,
) -> tuple[Relation, ...]:
    """Exact Graver basis via the Hilbert basis of the Lawrence lifting."""

    width = radial_degree + 1
    operator_count = (1,) * width + (0,) * width
    polynomial_count = (0,) * width + (1,) * width
    total_level = tuple(range(width)) + tuple(range(width))
    configuration = (operator_count, polynomial_count, total_level)
    equations = tuple(row + tuple(-value for value in row) for row in configuration)
    lifted = normaliz_hilbert_basis(equations, executable)

    relations: set[Relation] = set()
    state_width = 2 * width
    for element in lifted:
        positive = element[:state_width]
        negative = element[state_width:]
        if any(left and right for left, right in zip(positive, negative)):
            # The diagonal e_i+e_i generators remove every overlap.  A
            # nontrivial Lawrence generator is therefore disjoint.
            continue
        if not any(positive) or not any(negative):
            continue
        left = state_from_dense(positive, radial_degree)
        right = state_from_dense(negative, radial_degree)
        relation = orient_relation(left, right)
        verify_relation(relation)
        relations.add(relation)
    return tuple(sorted(relations))


def state_counts(state: State) -> tuple[int, int]:
    return sum(state.operator), sum(state.polynomial)


def state_level(state: State) -> int:
    return sum(
        level * multiplicity
        for side in (state.operator, state.polynomial)
        for level, multiplicity in enumerate(side)
    )


def verify_relation(relation: Relation) -> None:
    if state_counts(relation.left) != state_counts(relation.right):
        raise AssertionError("colour counts differ across relation")
    if state_level(relation.left) != state_level(relation.right):
        raise AssertionError("total levels differ across relation")
    dense_left = relation.left.operator + relation.left.polynomial
    dense_right = relation.right.operator + relation.right.polynomial
    if any(left and right for left, right in zip(dense_left, dense_right)):
        raise AssertionError("Graver relation is not sign-disjoint")


def relation_support_size(relation: Relation) -> int:
    dense_left = relation.left.operator + relation.left.polynomial
    dense_right = relation.right.operator + relation.right.polynomial
    return sum(bool(left or right) for left, right in zip(dense_left, dense_right))


def active_levels(relation: Relation) -> tuple[int, ...]:
    width = len(relation.left.operator)
    return tuple(
        level
        for level in range(width)
        if any(
            side[level]
            for state in (relation.left, relation.right)
            for side in (state.operator, state.polynomial)
        )
    )


def transform_relation(
    relation: Relation,
    *,
    reverse: bool,
    swap_colours: bool,
    swap_sides: bool,
) -> Relation:
    def transform_state(state: State) -> State:
        operator = state.operator[::-1] if reverse else state.operator
        polynomial = state.polynomial[::-1] if reverse else state.polynomial
        if swap_colours:
            operator, polynomial = polynomial, operator
        return State(operator, polynomial)

    left = transform_state(relation.left)
    right = transform_state(relation.right)
    if swap_sides:
        left, right = right, left
    return Relation(left, right)


def primitive_normalization(relation: Relation) -> Relation:
    """Remove common level translation/dilation, then quotient symmetries."""

    levels = active_levels(relation)
    first = levels[0]
    differences = tuple(level - first for level in levels[1:])
    divisor = math.gcd(*differences) if differences else 1
    divisor = max(divisor, 1)
    span = (levels[-1] - first) // divisor

    def normalize_state(state: State) -> State:
        operator = [0] * (span + 1)
        polynomial = [0] * (span + 1)
        for level, value in enumerate(state.operator):
            if value:
                operator[(level - first) // divisor] = value
        for level, value in enumerate(state.polynomial):
            if value:
                polynomial[(level - first) // divisor] = value
        return State(tuple(operator), tuple(polynomial))

    normalized = Relation(
        normalize_state(relation.left),
        normalize_state(relation.right),
    )
    variants = tuple(
        transform_relation(
            normalized,
            reverse=reverse,
            swap_colours=swap_colours,
            swap_sides=swap_sides,
        )
        for reverse in (False, True)
        for swap_colours in (False, True)
        for swap_sides in (False, True)
    )
    answer = min(variants)
    verify_relation(answer)
    return answer


def sparse_side(side: Vector) -> list[list[int]]:
    return [
        [level, multiplicity]
        for level, multiplicity in enumerate(side)
        if multiplicity
    ]


def relation_record(relation: Relation) -> dict[str, Any]:
    return {
        "span": len(relation.left.operator) - 1,
        "support_size": relation_support_size(relation),
        "colour_counts": list(state_counts(relation.left)),
        "total_level": state_level(relation.left),
        "left": {
            "operator": sparse_side(relation.left.operator),
            "polynomial": sparse_side(relation.left.polynomial),
        },
        "right": {
            "operator": sparse_side(relation.right.operator),
            "polynomial": sparse_side(relation.right.polynomial),
        },
    }


def relation_digest(relation: Relation) -> str:
    payload = json.dumps(
        relation_record(relation),
        sort_keys=True,
        separators=(",", ":"),
    ).encode()
    return hashlib.sha256(payload).hexdigest()


def positive_partition(side: Vector) -> tuple[int, ...]:
    return tuple(sorted((value for value in side if value), reverse=True))


def side_metadata(state: State, span: int) -> dict[str, Any]:
    operator_count, polynomial_count = state_counts(state)
    total_level = state_level(state)
    total_capacity = span * (operator_count + polynomial_count)
    operator_level = sum(
        level * value for level, value in enumerate(state.operator)
    )
    polynomial_level = sum(
        level * value for level, value in enumerate(state.polynomial)
    )
    return {
        "radial_vector": [total_level, total_capacity - total_level],
        "operator_partition": list(positive_partition(state.operator)),
        "polynomial_partition": list(positive_partition(state.polynomial)),
        "marked_operator_partition": sparse_side(state.operator),
        "marked_polynomial_partition": sparse_side(state.polynomial),
        "marked_radial_levels": [operator_level, polynomial_level],
    }


def base_digits(number: int, prime: int) -> tuple[int, ...]:
    if number == 0:
        return (0,)
    answer = []
    remaining = number
    while remaining:
        answer.append(remaining % prime)
        remaining //= prime
    return tuple(answer)


def kummer_carries(addends: tuple[int, ...], prime: int) -> tuple[int, ...]:
    target = sum(addends)
    if not addends:
        return ()
    digit_rows = tuple(base_digits(value, prime) for value in addends)
    target_digits = base_digits(target, prime)
    width = max(max(map(len, digit_rows)), len(target_digits)) + 1
    incoming = 0
    carries = []
    for position in range(width):
        column = incoming + sum(
            digits[position] if position < len(digits) else 0
            for digits in digit_rows
        )
        target_digit = (
            target_digits[position] if position < len(target_digits) else 0
        )
        if (column - target_digit) % prime:
            raise AssertionError("invalid base-p carry column")
        outgoing = (column - target_digit) // prime
        carries.append(outgoing)
        incoming = outgoing
    while carries and carries[-1] == 0:
        carries.pop()
    return tuple(carries)


@cache
def factorial_valuation(number: int, prime: int) -> int:
    answer = 0
    quotient = number
    while quotient:
        quotient //= prime
        answer += quotient
    return answer


def legendre_layers(number: int, prime: int) -> tuple[int, ...]:
    answer = []
    quotient = number
    while quotient:
        quotient //= prime
        if quotient:
            answer.append(quotient)
    return tuple(answer)


_FACTORIAL_UNIT_TABLES: dict[tuple[int, int], list[int]] = {}


def factorial_unit(number: int, prime: int, unit_power: int) -> int:
    """Return the p-free part of n! modulo p**unit_power.

    The strengthened p^3 probes reach substantially larger factorial
    arguments than the pilot run.  Extending one prefix table per modulus
    avoids recomputing every smaller factorial for each packet.
    """

    modulus = prime**unit_power
    table = _FACTORIAL_UNIT_TABLES.setdefault((prime, unit_power), [1])
    for value in range(len(table), number + 1):
        unit = value
        while unit % prime == 0:
            unit //= prime
        table.append(table[-1] * (unit % modulus) % modulus)
    return table[number]


def multinomial_data(
    partition: tuple[int, ...],
    scale: int,
    prime: int,
    unit_power: int,
) -> tuple[int, int, tuple[int, ...]]:
    scaled = tuple(scale * value for value in partition)
    total = sum(scaled)
    valuation = factorial_valuation(total, prime) - sum(
        factorial_valuation(value, prime) for value in scaled
    )
    carries = kummer_carries(scaled, prime)
    if sum(carries) != valuation:
        raise AssertionError((partition, scale, prime, carries, valuation))
    modulus = prime**unit_power
    denominator = math.prod(
        factorial_unit(value, prime, unit_power) for value in scaled
    ) % modulus
    unit = (
        factorial_unit(total, prime, unit_power)
        * pow(denominator, -1, modulus)
    ) % modulus
    return valuation, unit, carries


def group_ring_trace(side: Vector, scale: int, order: int) -> list[int]:
    """C_order group-ring coefficients; Fourier evaluation gives traces."""

    histogram = [0] * order
    for level, multiplicity in enumerate(side):
        histogram[level % order] += scale * multiplicity
    return histogram


def torsion_snapshot(state: State, scale: int, order: int) -> dict[str, Any]:
    operator_level = sum(
        level * value for level, value in enumerate(state.operator)
    )
    polynomial_level = sum(
        level * value for level, value in enumerate(state.polynomial)
    )
    return {
        "operator_group_ring": group_ring_trace(state.operator, scale, order),
        "polynomial_group_ring": group_ring_trace(
            state.polynomial, scale, order
        ),
        "marked_monomial_class": [
            scale * operator_level % order,
            scale * polynomial_level % order,
        ],
    }


def probe_side(
    state: State,
    span: int,
    *,
    prime: int,
    exponent: int,
    quotient: int,
    unit_power: int,
    torsion_orders: tuple[int, ...],
) -> dict[str, Any]:
    scale = quotient * prime**exponent
    operator_partition = positive_partition(state.operator)
    polynomial_partition = positive_partition(state.polynomial)
    operator_valuation, operator_unit, operator_carries = multinomial_data(
        operator_partition,
        scale,
        prime,
        unit_power,
    )
    polynomial_valuation, polynomial_unit, polynomial_carries = multinomial_data(
        polynomial_partition,
        scale,
        prime,
        unit_power,
    )

    operator_count, polynomial_count = state_counts(state)
    total_level = state_level(state)
    total_capacity = span * (operator_count + polynomial_count)
    radial = (scale * total_level, scale * (total_capacity - total_level))
    radial_valuations = tuple(
        factorial_valuation(value, prime) for value in radial
    )
    radial_valuation = sum(radial_valuations)
    modulus = prime**unit_power
    radial_unit = math.prod(
        factorial_unit(value, prime, unit_power) for value in radial
    ) % modulus
    total_unit = operator_unit * polynomial_unit * radial_unit % modulus

    return {
        "scale": scale,
        "low_digits": {
            "operator_partition": [
                list(base_digits(scale * value, prime))
                for value in operator_partition
            ],
            "polynomial_partition": [
                list(base_digits(scale * value, prime))
                for value in polynomial_partition
            ],
        },
        "valuations": {
            "operator_multinomial": operator_valuation,
            "polynomial_multinomial": polynomial_valuation,
            "radial_factorials": list(radial_valuations),
            "total": operator_valuation
            + polynomial_valuation
            + radial_valuation,
        },
        "carry_positions": {
            "operator": [
                [position, carry]
                for position, carry in enumerate(operator_carries)
                if carry
            ],
            "polynomial": [
                [position, carry]
                for position, carry in enumerate(polynomial_carries)
                if carry
            ],
        },
        "legendre_layers": {
            "radial_x": list(legendre_layers(radial[0], prime)),
            "radial_y": list(legendre_layers(radial[1], prime)),
        },
        "factorial_units": {
            "modulus": modulus,
            "operator_multinomial": operator_unit,
            "polynomial_multinomial": polynomial_unit,
            "radial_factorials": radial_unit,
            "total": total_unit,
        },
        "torsion_character_traces": {
            str(order): torsion_snapshot(state, scale, order)
            for order in torsion_orders
        },
    }


def comparison_signature(
    metadata: dict[str, Any],
    probes: list[dict[str, Any]],
) -> dict[str, Any]:
    """Requested invariant, excluding exact marked levels and their sums."""

    return {
        "radial_vector": metadata["radial_vector"],
        "operator_partition": metadata["operator_partition"],
        "polynomial_partition": metadata["polynomial_partition"],
        "probes": probes,
    }


def signature_digest(signature: dict[str, Any]) -> str:
    payload = json.dumps(
        signature,
        sort_keys=True,
        separators=(",", ":"),
    ).encode()
    return hashlib.sha256(payload).hexdigest()


def state_identifier(state: State, span: int) -> str:
    payload = json.dumps(
        {
            "span": span,
            "operator": list(state.operator),
            "polynomial": list(state.polynomial),
        },
        sort_keys=True,
        separators=(",", ":"),
    ).encode()
    return hashlib.sha256(payload).hexdigest()


def difference_relation(left: State, right: State) -> Relation:
    positive_operator = tuple(
        max(left_value - right_value, 0)
        for left_value, right_value in zip(
            left.operator, right.operator, strict=True
        )
    )
    positive_polynomial = tuple(
        max(left_value - right_value, 0)
        for left_value, right_value in zip(
            left.polynomial, right.polynomial, strict=True
        )
    )
    negative_operator = tuple(
        max(right_value - left_value, 0)
        for left_value, right_value in zip(
            left.operator, right.operator, strict=True
        )
    )
    negative_polynomial = tuple(
        max(right_value - left_value, 0)
        for left_value, right_value in zip(
            left.polynomial, right.polynomial, strict=True
        )
    )
    relation = orient_relation(
        State(positive_operator, positive_polynomial),
        State(negative_operator, negative_polynomial),
    )
    verify_relation(relation)
    return relation


def first_difference(left: Any, right: Any, path: str = "") -> str | None:
    if type(left) is not type(right):
        return path or "$"
    if isinstance(left, dict):
        keys = sorted(set(left) | set(right))
        for key in keys:
            child = f"{path}.{key}" if path else key
            if key not in left or key not in right:
                return child
            difference = first_difference(left[key], right[key], child)
            if difference is not None:
                return difference
        return None
    if isinstance(left, list):
        if len(left) != len(right):
            return f"{path}.length"
        for index, (left_value, right_value) in enumerate(zip(left, right)):
            difference = first_difference(
                left_value,
                right_value,
                f"{path}[{index}]",
            )
            if difference is not None:
                return difference
        return None
    return None if left == right else (path or "$")


def support_pair(relation: Relation) -> tuple[tuple[int, ...], tuple[int, ...]]:
    operator = tuple(
        level
        for level in range(len(relation.left.operator))
        if relation.left.operator[level] or relation.right.operator[level]
    )
    polynomial = tuple(
        level
        for level in range(len(relation.left.polynomial))
        if relation.left.polynomial[level] or relation.right.polynomial[level]
    )
    return operator, polynomial


def support_is_subset(
    relation: Relation,
    operator_support: tuple[int, ...],
    polynomial_support: tuple[int, ...],
) -> bool:
    operator_allowed = set(operator_support)
    polynomial_allowed = set(polynomial_support)
    return all(
        not (relation.left.operator[level] or relation.right.operator[level])
        or level in operator_allowed
        for level in range(len(relation.left.operator))
    ) and all(
        not (
            relation.left.polynomial[level]
            or relation.right.polynomial[level]
        )
        or level in polynomial_allowed
        for level in range(len(relation.left.polynomial))
    )


def exact_scalar_factorial_collision(relation: Relation) -> bool:
    """All-scale scalar factorial equality, by scaled-factorial rigidity."""

    left = sorted(
        positive_partition(relation.left.operator)
        + positive_partition(relation.left.polynomial),
        reverse=True,
    )
    right = sorted(
        positive_partition(relation.right.operator)
        + positive_partition(relation.right.polynomial),
        reverse=True,
    )
    return left == right


def exact_marked_partition_collision(relation: Relation) -> bool:
    return (
        positive_partition(relation.left.operator)
        == positive_partition(relation.right.operator)
        and positive_partition(relation.left.polynomial)
        == positive_partition(relation.right.polynomial)
    )


def combined_partition(state: State) -> tuple[int, ...]:
    return tuple(
        sorted(
            positive_partition(state.operator)
            + positive_partition(state.polynomial),
            reverse=True,
        )
    )


def fraction_record(value: Fraction) -> list[int]:
    return [value.numerator, value.denominator]


def asymptotic_partition_separator(
    relation: Relation,
) -> dict[str, Any] | None:
    """First exact Stirling component separating the scalar factorial rays."""

    left = combined_partition(relation.left)
    right = combined_partition(relation.right)
    if left == right:
        return None

    components = (
        ("total_mass", sum(left), sum(right)),
        (
            "entropy_exponential",
            math.prod(value**value for value in left),
            math.prod(value**value for value in right),
        ),
        ("part_count", len(left), len(right)),
        ("part_product", math.prod(left), math.prod(right)),
    )
    for component, left_value, right_value in components:
        if left_value != right_value:
            return {
                "component": component,
                "left": left_value,
                "right": right_value,
            }

    distinct_parts = len(set(left) | set(right))
    for index in range(1, distinct_parts + 1):
        power = 2 * index - 1
        left_value = sum(
            (Fraction(1, value**power) for value in left),
            start=Fraction(),
        )
        right_value = sum(
            (Fraction(1, value**power) for value in right),
            start=Fraction(),
        )
        if left_value != right_value:
            return {
                "component": f"reciprocal_power_sum_{power}",
                "left": fraction_record(left_value),
                "right": fraction_record(right_value),
            }

    raise AssertionError("unequal partitions have identical Stirling data")


def scalar_valuation_signature(probes: list[dict[str, Any]]) -> list[int]:
    return [probe["valuations"]["total"] for probe in probes]


def scalar_unit_signature(probes: list[dict[str, Any]]) -> list[list[int]]:
    return [
        [
            probe["factorial_units"]["modulus"],
            probe["factorial_units"]["total"],
        ]
        for probe in probes
    ]


def marked_digit_carry_signature(
    probes: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    return [
        {
            "low_digits": probe["low_digits"],
            "side_valuations": {
                "operator_multinomial": probe["valuations"][
                    "operator_multinomial"
                ],
                "polynomial_multinomial": probe["valuations"][
                    "polynomial_multinomial"
                ],
            },
            "carry_positions": probe["carry_positions"],
        }
        for probe in probes
    ]


def marked_unit_signature(probes: list[dict[str, Any]]) -> list[list[int]]:
    return [
        [
            probe["factorial_units"]["operator_multinomial"],
            probe["factorial_units"]["polynomial_multinomial"],
        ]
        for probe in probes
    ]


def first_adelic_probe_separator(
    left_probes: list[dict[str, Any]],
    right_probes: list[dict[str, Any]],
) -> dict[str, Any] | None:
    for index, (left, right) in enumerate(
        zip(left_probes, right_probes, strict=True)
    ):
        left_value = left["valuations"]["total"]
        right_value = right["valuations"]["total"]
        if left_value != right_value:
            return {
                "layer": "valuation",
                "probe_offset": index,
                "left": left_value,
                "right": right_value,
            }
    for index, (left, right) in enumerate(
        zip(left_probes, right_probes, strict=True)
    ):
        left_value = left["factorial_units"]["total"]
        right_value = right["factorial_units"]["total"]
        if left_value != right_value:
            if (
                left["factorial_units"]["modulus"]
                != right["factorial_units"]["modulus"]
            ):
                raise AssertionError("probe moduli differ across a relation")
            return {
                "layer": "factorial_unit",
                "probe_offset": index,
                "modulus": left["factorial_units"]["modulus"],
                "left": left_value,
                "right": right_value,
            }
    return None


def first_packet_decoration_separator(
    relation: Relation,
    left_probes: list[dict[str, Any]],
    right_probes: list[dict[str, Any]],
    torsion_orders: tuple[int, ...],
) -> dict[str, Any] | None:
    left_marked = marked_digit_carry_signature(left_probes)
    right_marked = marked_digit_carry_signature(right_probes)
    for index, (left, right) in enumerate(
        zip(left_marked, right_marked, strict=True)
    ):
        difference = first_difference(left, right)
        if difference is not None:
            return {
                "layer": "marked_digit_or_carry",
                "probe_offset": index,
                "first_field": difference,
            }
    left_units = marked_unit_signature(left_probes)
    right_units = marked_unit_signature(right_probes)
    for index, (left, right) in enumerate(
        zip(left_units, right_units, strict=True)
    ):
        if left != right:
            return {
                "layer": "marked_factorial_unit",
                "probe_offset": index,
                "left": left,
                "right": right,
            }
    if not exact_marked_partition_collision(relation):
        return {
            "layer": "marked_partition",
            "left": {
                "operator": list(positive_partition(relation.left.operator)),
                "polynomial": list(
                    positive_partition(relation.left.polynomial)
                ),
            },
            "right": {
                "operator": list(
                    positive_partition(relation.right.operator)
                ),
                "polynomial": list(
                    positive_partition(relation.right.polynomial)
                ),
            },
        }
    for order in torsion_orders:
        left = torsion_snapshot(relation.left, 1, order)
        right = torsion_snapshot(relation.right, 1, order)
        if left != right:
            return {
                "layer": "torsion_character",
                "order": order,
                "first_field": first_difference(left, right),
                "left": left,
                "right": right,
            }
    return None


def six_step_collision_relation(
    radial_degree: int,
    operator_low: int,
    polynomial_low: int,
) -> Relation:
    if not 0 <= operator_low <= radial_degree - 6:
        raise ValueError("six-step operator pair lies outside the span")
    if not 0 <= polynomial_low <= radial_degree - 4:
        raise ValueError("shifted adjacent polynomial pairs lie outside the span")
    operator_left = [0] * (radial_degree + 1)
    operator_right = [0] * (radial_degree + 1)
    polynomial_left = [0] * (radial_degree + 1)
    polynomial_right = [0] * (radial_degree + 1)
    operator_left[operator_low + 6] = 1
    operator_right[operator_low] = 1
    polynomial_left[polynomial_low] = 1
    polynomial_left[polynomial_low + 1] = 1
    polynomial_right[polynomial_low + 3] = 1
    polynomial_right[polynomial_low + 4] = 1
    return orient_relation(
        State(tuple(operator_left), tuple(polynomial_left)),
        State(tuple(operator_right), tuple(polynomial_right)),
    )


def expected_six_step_collision_orbits(
    radial_degree: int,
) -> set[Relation]:
    return {
        primitive_normalization(
            six_step_collision_relation(
                radial_degree,
                operator_low,
                polynomial_low,
            )
        )
        for operator_low in range(radial_degree - 5)
        for polynomial_low in range(radial_degree - 3)
    }


def configured_torsion_collision(
    relation: Relation,
    torsion_orders: tuple[int, ...],
) -> bool:
    return all(
        torsion_snapshot(relation.left, 1, order)
        == torsion_snapshot(relation.right, 1, order)
        for order in torsion_orders
    )


def first_separating_torsion_order(
    relation: Relation,
    limit: int,
) -> int | None:
    for order in range(2, limit + 1):
        if (
            torsion_snapshot(relation.left, 1, order)
            != torsion_snapshot(relation.right, 1, order)
        ):
            return order
    return None


def weak_compositions(total: int, length: int):
    if length == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for tail in weak_compositions(total - first, length - 1):
            yield (first,) + tail


def exact_fiber_states(relation: Relation) -> tuple[State, ...]:
    operator_support, polynomial_support = support_pair(relation)
    operator_count, polynomial_count = state_counts(relation.left)
    target_level = state_level(relation.left)
    width = len(relation.left.operator)
    states = []
    for operator_values in weak_compositions(
        operator_count, len(operator_support)
    ):
        operator_level = sum(
            level * value
            for level, value in zip(
                operator_support, operator_values, strict=True
            )
        )
        for polynomial_values in weak_compositions(
            polynomial_count, len(polynomial_support)
        ):
            polynomial_level = sum(
                level * value
                for level, value in zip(
                    polynomial_support, polynomial_values, strict=True
                )
            )
            if operator_level + polynomial_level != target_level:
                continue
            operator = [0] * width
            polynomial = [0] * width
            for level, value in zip(
                operator_support, operator_values, strict=True
            ):
                operator[level] = value
            for level, value in zip(
                polynomial_support, polynomial_values, strict=True
            ):
                polynomial[level] = value
            states.append(State(tuple(operator), tuple(polynomial)))
    return tuple(sorted(states))


def circuit_connected_in_exact_fiber(
    relation: Relation,
) -> tuple[int, bool]:
    states = exact_fiber_states(relation)
    frontier = [relation.left]
    visited = {relation.left}
    while frontier:
        state = frontier.pop()
        if state == relation.right:
            return len(states), True
        for neighbor in states:
            if neighbor in visited:
                continue
            difference = difference_relation(state, neighbor)
            if relation_support_size(difference) <= 4:
                visited.add(neighbor)
                frontier.append(neighbor)
    return len(states), False


def build_experiment(
    *,
    radial_degree: int,
    max_support: int | None,
    primitive_only: bool,
    primes: tuple[int, ...],
    max_exponent: int,
    max_quotient: int,
    unit_power: int,
    torsion_orders: tuple[int, ...],
    normaliz: str,
) -> dict[str, Any]:
    raw_basis = colored_graver_basis(radial_degree, normaliz)
    normalized_basis = tuple(
        sorted({primitive_normalization(relation) for relation in raw_basis})
    )
    raw_ids = {
        relation: f"G{index:04d}"
        for index, relation in enumerate(raw_basis, start=1)
    }

    candidates = tuple(
        relation
        for relation in normalized_basis
        if relation_support_size(relation) >= 5
        and (
            max_support is None
            or relation_support_size(relation) <= max_support
        )
        and all(state_counts(relation.left))
    )
    packet_ids = {
        relation: f"P{index:04d}"
        for index, relation in enumerate(candidates, start=1)
    }

    probe_index = [
        {"prime": prime, "e": exponent, "q": quotient,
         "scale": quotient * prime**exponent}
        for prime in primes
        for exponent in range(1, max_exponent + 1)
        for quotient in range(1, max_quotient + 1)
    ]

    packets = []
    state_signatures: dict[str, dict[str, Any]] = {}
    state_values: dict[str, State] = {}
    state_records: dict[str, dict[str, Any]] = {}
    state_occurrences: dict[str, list[dict[str, str]]] = {}
    scalar_collisions = []
    marked_partition_collisions = []
    full_collisions = []
    exact_configured_collisions = []
    valuation_collision_ids = []
    finite_adelic_collision_ids = []
    finite_adelic_false_collision_ids = []
    first_adelic_counts: Counter[str] = Counter()
    first_adelic_probe_counts: Counter[
        tuple[str, int, int, int]
    ] = Counter()
    asymptotic_component_counts: Counter[str] = Counter()
    exact_scalar_decoration_counts: Counter[str] = Counter()

    for relation in candidates:
        span = len(relation.left.operator) - 1
        left_metadata = side_metadata(relation.left, span)
        right_metadata = side_metadata(relation.right, span)
        left_probes = []
        right_probes = []
        for probe in probe_index:
            arguments = {
                "prime": probe["prime"],
                "exponent": probe["e"],
                "quotient": probe["q"],
                "unit_power": unit_power,
                "torsion_orders": torsion_orders,
            }
            left_probes.append(
                probe_side(relation.left, span, **arguments)
            )
            right_probes.append(
                probe_side(relation.right, span, **arguments)
            )

        left_signature = comparison_signature(left_metadata, left_probes)
        right_signature = comparison_signature(right_metadata, right_probes)
        left_digest = signature_digest(left_signature)
        right_digest = signature_digest(right_signature)
        collision = left_digest == right_digest
        if collision and left_signature != right_signature:
            raise AssertionError("SHA-256 collision in tomography signatures")
        packet_id = packet_ids[relation]
        scalar_collision = exact_scalar_factorial_collision(relation)
        marked_collision = exact_marked_partition_collision(relation)
        asymptotic_separator = asymptotic_partition_separator(relation)
        if (asymptotic_separator is None) != scalar_collision:
            raise AssertionError(
                "Stirling separator and scaled-factorial rigidity disagree"
            )
        first_adelic_separator = first_adelic_probe_separator(
            left_probes, right_probes
        )
        if first_adelic_separator is not None:
            first_adelic_separator = {
                **first_adelic_separator,
                "probe": probe_index[
                    first_adelic_separator["probe_offset"]
                ],
            }
        first_adelic_layer = (
            first_adelic_separator["layer"]
            if first_adelic_separator is not None
            else None
        )
        valuation_collision = first_adelic_layer != "valuation"
        finite_adelic_collision = first_adelic_separator is None
        decoration_separator = first_packet_decoration_separator(
            relation,
            left_probes,
            right_probes,
            torsion_orders,
        )
        if (
            decoration_separator is not None
            and "probe_offset" in decoration_separator
        ):
            decoration_separator = {
                **decoration_separator,
                "probe": probe_index[
                    decoration_separator["probe_offset"]
                ],
            }
        decoration_layer = (
            decoration_separator["layer"]
            if decoration_separator is not None
            else None
        )
        torsion_collision = configured_torsion_collision(
            relation, torsion_orders
        )
        exact_configured_collision = marked_collision and torsion_collision
        if collision != exact_configured_collision:
            raise AssertionError(
                "finite probes disagree with the exact partition/character test"
            )
        first_torsion_separator = first_separating_torsion_order(
            relation, span + 1
        )
        if scalar_collision:
            scalar_collisions.append(packet_id)
            exact_scalar_decoration_counts[
                decoration_layer or "configured_indistinguishable"
            ] += 1
        asymptotic_component_counts[
            (
                asymptotic_separator["component"]
                if asymptotic_separator is not None
                else "exact_scalar_factorial_collision"
            )
        ] += 1
        first_adelic_counts[
            first_adelic_layer or "finite_adelic_collision"
        ] += 1
        if first_adelic_separator is not None:
            first_probe = first_adelic_separator["probe"]
            first_adelic_probe_counts[
                (
                    first_adelic_layer,
                    first_probe["prime"],
                    first_probe["e"],
                    first_probe["q"],
                )
            ] += 1
        if valuation_collision:
            valuation_collision_ids.append(packet_id)
        if finite_adelic_collision:
            finite_adelic_collision_ids.append(packet_id)
            if not scalar_collision:
                finite_adelic_false_collision_ids.append(packet_id)
        if marked_collision:
            marked_partition_collisions.append(packet_id)
        if collision:
            full_collisions.append(packet_id)
        if exact_configured_collision:
            exact_configured_collisions.append(packet_id)
        collision_fiber = None
        if exact_configured_collision:
            fiber_size, circuit_connected = circuit_connected_in_exact_fiber(
                relation
            )
            collision_fiber = {
                "state_count": fiber_size,
                "endpoints_connected_by_support_at_most_four_moves": (
                    circuit_connected
                ),
            }

        left_adelic_weight_signature = {
            "valuations": scalar_valuation_signature(left_probes),
            "units": scalar_unit_signature(left_probes),
        }
        right_adelic_weight_signature = {
            "valuations": scalar_valuation_signature(right_probes),
            "units": scalar_unit_signature(right_probes),
        }
        collision_tomography = None
        if exact_configured_collision or (
            finite_adelic_collision and not scalar_collision
        ):
            collision_tomography = [
                {
                    "probe": probe,
                    "left": left_probe,
                    "right": right_probe,
                }
                for probe, left_probe, right_probe in zip(
                    probe_index,
                    left_probes,
                    right_probes,
                    strict=True,
                )
            ]

        left_state_id = state_identifier(relation.left, span)
        right_state_id = state_identifier(relation.right, span)
        if not primitive_only:
            for label, state, metadata, probes, digest, signature, state_id in (
                (
                    "left",
                    relation.left,
                    left_metadata,
                    left_probes,
                    left_digest,
                    left_signature,
                    left_state_id,
                ),
                (
                    "right",
                    relation.right,
                    right_metadata,
                    right_probes,
                    right_digest,
                    right_signature,
                    right_state_id,
                ),
            ):
                previous_signature = state_signatures.setdefault(
                    state_id, signature
                )
                if previous_signature != signature:
                    raise AssertionError(
                        "one state acquired two tomography signatures"
                    )
                state_values.setdefault(state_id, state)
                state_records.setdefault(
                    state_id,
                    {
                        "id": state_id,
                        "span": span,
                        "metadata": metadata,
                        "signature_sha256": digest,
                        "tomography": [
                            {**probe, "data": probe_data}
                            for probe, probe_data in zip(
                                probe_index,
                                probes,
                                strict=True,
                            )
                        ],
                    },
                )
                state_occurrences.setdefault(state_id, []).append(
                    {"packet": packet_id, "side": label, "digest": digest}
                )

        packets.append(
            {
                "id": packet_id,
                "relation": relation_record(relation),
                "relation_sha256": relation_digest(relation),
                "left_metadata": left_metadata,
                "right_metadata": right_metadata,
                "exact_scalar_factorial_collision": scalar_collision,
                "asymptotic_partition_separator": asymptotic_separator,
                "all_configured_total_valuations_equal": (
                    valuation_collision
                ),
                "finite_adelic_weight_collision": finite_adelic_collision,
                "finite_adelic_false_collision": (
                    finite_adelic_collision and not scalar_collision
                ),
                "first_adelic_probe_separator": first_adelic_separator,
                "first_packet_decoration_separator": decoration_separator,
                "left_adelic_weight_signature_sha256": signature_digest(
                    left_adelic_weight_signature
                ),
                "right_adelic_weight_signature_sha256": signature_digest(
                    right_adelic_weight_signature
                ),
                "surviving_collision_tomography": collision_tomography,
                "exact_marked_partition_collision": marked_collision,
                "configured_torsion_collision": torsion_collision,
                "exact_all_scale_configured_signature_collision": (
                    exact_configured_collision
                ),
                "first_separating_torsion_order": first_torsion_separator,
                "collision_exact_fiber": collision_fiber,
                "full_probe_collision": collision,
                "first_separating_field": first_difference(
                    left_signature, right_signature
                ),
                "left_state_id": left_state_id,
                "right_state_id": right_state_id,
                "left_signature_sha256": left_digest,
                "right_signature_sha256": right_digest,
            }
        )

    collision_groups: dict[str, list[str]] = {}
    for state_id, signature in state_signatures.items():
        collision_groups.setdefault(signature_digest(signature), []).append(
            state_id
        )
    nontrivial_global_groups = []
    not_pairwise_circuit_groups = []
    for digest, state_ids in sorted(collision_groups.items()):
        if len(state_ids) < 2:
            continue
        pair_relations = []
        all_closed_circuits = True
        for left_index, left_id in enumerate(sorted(state_ids)):
            for right_id in sorted(state_ids)[left_index + 1 :]:
                relation = difference_relation(
                    state_values[left_id], state_values[right_id]
                )
                support_size = relation_support_size(relation)
                normalized = primitive_normalization(relation)
                primitive = normalized in normalized_basis
                closed_circuit = support_size <= 4
                all_closed_circuits &= closed_circuit
                pair_relations.append(
                    {
                        "left_state_id": left_id,
                        "right_state_id": right_id,
                        "support_size": support_size,
                        "graver_primitive": primitive,
                        "closed_circuit_range": closed_circuit,
                        "relation": relation_record(relation),
                    }
                )
        group = {
            "signature_sha256": digest,
            "state_ids": sorted(state_ids),
            "all_pairs_explained_by_closed_circuits": all_closed_circuits,
            "pair_relations": pair_relations,
            "occurrences": [
                occurrence
                for state_id in sorted(state_ids)
                for occurrence in state_occurrences[state_id]
            ],
        }
        nontrivial_global_groups.append(group)
        if not all_closed_circuits:
            not_pairwise_circuit_groups.append(group)

    candidate_supports = sorted({support_pair(relation) for relation in candidates})

    if (
        radial_degree in (6, 7)
        and (max_support is None or max_support >= 6)
        and set(torsion_orders) == {2, 3}
    ):
        expected_collisions = expected_six_step_collision_orbits(
            radial_degree
        )
        actual_collisions = {
            relation
            for relation in candidates
            if packet_ids[relation] in exact_configured_collisions
        }
        if actual_collisions != expected_collisions:
            raise AssertionError(
                "six-step configured collision family changed: "
                f"{actual_collisions} != {expected_collisions}"
            )

    raw_ids_by_support: dict[
        tuple[tuple[int, ...], tuple[int, ...]], list[str]
    ] = {}
    for relation in raw_basis:
        raw_ids_by_support.setdefault(support_pair(relation), []).append(
            raw_ids[relation]
        )
    packet_ids_by_support: dict[
        tuple[tuple[int, ...], tuple[int, ...]], list[str]
    ] = {}
    for relation in candidates:
        packet_ids_by_support.setdefault(support_pair(relation), []).append(
            packet_ids[relation]
        )

    semigroups = []
    for index, (operator_support, polynomial_support) in enumerate(
        candidate_supports, start=1
    ):
        basis_ids = []
        for operator_mask in range(1 << len(operator_support)):
            operator_subset = tuple(
                level
                for bit, level in enumerate(operator_support)
                if operator_mask & (1 << bit)
            )
            for polynomial_mask in range(1 << len(polynomial_support)):
                polynomial_subset = tuple(
                    level
                    for bit, level in enumerate(polynomial_support)
                    if polynomial_mask & (1 << bit)
                )
                basis_ids.extend(
                    raw_ids_by_support.get(
                        (operator_subset, polynomial_subset), []
                    )
                )
        basis_ids.sort(key=lambda identifier: int(identifier[1:]))
        exact_packet_ids = packet_ids_by_support[
            (operator_support, polynomial_support)
        ]
        semigroups.append(
            {
                "id": f"S{index:04d}",
                "operator_support": list(operator_support),
                "polynomial_support": list(polynomial_support),
                "graver_basis_ids": basis_ids,
                "graver_basis_size": len(basis_ids),
                "exact_unresolved_packet_ids": exact_packet_ids,
            }
        )

    result = {
        "schema": "binary-gvc-prime-power-tomography-v2",
        "status": {
            "kind": "bounded_experiment",
            "claim": (
                "projected two-colour Graver census only; not the "
                "Hall-envelope proof and not a GVC counterexample"
            ),
        },
        "parameters": {
            "radial_degree": radial_degree,
            "max_support": max_support,
            "primitive_only": primitive_only,
            "primes": list(primes),
            "max_exponent": max_exponent,
            "max_quotient": max_quotient,
            "unit_power": unit_power,
            "torsion_orders": list(torsion_orders),
            "normaliz": normaliz_version(normaliz),
        },
        "model": {
            "columns": "R_i=(1,0,i), B_j=(0,1,j)",
            "closed_support_bound": 4,
            "candidate_filter": (
                "mixed-colour primitive Graver relations of support at least "
                "5, optionally truncated by max_support"
            ),
            "symmetry_quotient": [
                "common level translation",
                "common level dilation",
                "level reversal",
                "operator/polynomial colour exchange",
                "relation-side exchange",
            ],
            "torus_character_note": (
                "the return matrix fixes colour counts and total level, so "
                "the global torus character agrees across every relation"
            ),
            "signature_note": (
                "finite adelic weight collisions compare only total "
                "valuations and total p-free units; marked digits, carries, "
                "side units, partitions, and torsion characters form a "
                "separate decoration funnel"
            ),
            "record_scope": (
                "primitive relations only; global repeated-state groups "
                "omitted"
                if primitive_only
                else "primitive relations and global repeated-state groups"
            ),
        },
        "probe_index": probe_index,
        "summary": {
            "raw_full_support_graver_basis": len(raw_basis),
            "symmetry_normalized_graver_basis": len(normalized_basis),
            "candidate_packets": len(candidates),
            "candidate_support_size_distribution": {
                str(size): count
                for size, count in sorted(
                    Counter(
                        relation_support_size(relation)
                        for relation in candidates
                    ).items()
                )
            },
            "candidate_nonfree_support_semigroups": len(semigroups),
            "asymptotic_partition_classification": dict(
                sorted(asymptotic_component_counts.items())
            ),
            "first_adelic_probe_classification": dict(
                sorted(first_adelic_counts.items())
            ),
            "first_adelic_separator_probe_distribution": [
                {
                    "layer": layer,
                    "prime": prime,
                    "e": exponent,
                    "q": quotient,
                    "count": count,
                }
                for (
                    layer,
                    prime,
                    exponent,
                    quotient,
                ), count in sorted(first_adelic_probe_counts.items())
            ],
            "configured_total_valuation_collisions": len(
                valuation_collision_ids
            ),
            "finite_adelic_weight_collisions": len(
                finite_adelic_collision_ids
            ),
            "finite_adelic_false_collisions": len(
                finite_adelic_false_collision_ids
            ),
            "exact_scalar_factorial_collisions": len(scalar_collisions),
            "exact_scalar_collision_decoration_classification": dict(
                sorted(exact_scalar_decoration_counts.items())
            ),
            "exact_marked_partition_collisions": len(
                marked_partition_collisions
            ),
            "within_relation_full_probe_collisions": len(full_collisions),
            "within_relation_exact_all_scale_configured_collisions": len(
                exact_configured_collisions
            ),
            "global_distinct_state_full_probe_collision_groups": (
                None
                if primitive_only
                else len(nontrivial_global_groups)
            ),
            "global_collision_groups_not_pairwise_support_four": (
                None
                if primitive_only
                else len(not_pairwise_circuit_groups)
            ),
        },
        "collision_ids": {
            "configured_total_valuation": valuation_collision_ids,
            "finite_adelic_weight": finite_adelic_collision_ids,
            "finite_adelic_false": finite_adelic_false_collision_ids,
            "exact_scalar_factorial": scalar_collisions,
            "exact_marked_partition": marked_partition_collisions,
            "within_relation_full_probe": full_collisions,
            "within_relation_exact_all_scale_configured": (
                exact_configured_collisions
            ),
        },
        "global_full_probe_collision_groups": nontrivial_global_groups,
        "global_collision_groups_not_pairwise_support_four": (
            not_pairwise_circuit_groups
        ),
        "raw_graver_basis": [
            {"id": raw_ids[relation], **relation_record(relation)}
            for relation in raw_basis
        ],
        "nonfree_support_semigroups": semigroups,
        "states": [state_records[state_id] for state_id in sorted(state_records)],
        "packets": packets,
    }
    result["result_sha256"] = hashlib.sha256(
        json.dumps(result, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    return result


def compact_summary(result: dict[str, Any]) -> dict[str, Any]:
    collision_ids = set(
        result["collision_ids"][
            "within_relation_exact_all_scale_configured"
        ]
    )
    collision_packets = [
        packet for packet in result["packets"] if packet["id"] in collision_ids
    ]
    false_collision_ids = set(
        result["collision_ids"]["finite_adelic_false"]
    )
    collision_state_ids = {
        state_id
        for packet in collision_packets
        for state_id in (packet["left_state_id"], packet["right_state_id"])
    }
    return {
        "schema": result["schema"] + "-summary",
        "status": result["status"],
        "parameters": result["parameters"],
        "model": result["model"],
        "summary": result["summary"],
        "collision_ids": result["collision_ids"],
        "configured_collision_packets": collision_packets,
        "finite_adelic_false_collision_packets": [
            packet
            for packet in result["packets"]
            if packet["id"] in false_collision_ids
        ],
        "configured_collision_states": [
            state
            for state in result["states"]
            if state["id"] in collision_state_ids
        ],
        "global_full_probe_collision_groups": result[
            "global_full_probe_collision_groups"
        ],
        "global_collision_groups_not_pairwise_support_four": result[
            "global_collision_groups_not_pairwise_support_four"
        ],
        "full_result_sha256": result["result_sha256"],
    }


def write_json(path: Path, value: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    encoded = (json.dumps(value, indent=2, sort_keys=True) + "\n").encode()
    if path.suffix == ".gz":
        path.write_bytes(gzip.compress(encoded, compresslevel=9, mtime=0))
    else:
        path.write_bytes(encoded)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--radial-degree", type=int, default=6)
    parser.add_argument("--max-support", type=int)
    parser.add_argument("--primitive-only", action="store_true")
    parser.add_argument("--primes", default="2,3,5,7,11")
    parser.add_argument("--max-exponent", type=int, default=2)
    parser.add_argument("--max-quotient", type=int, default=3)
    parser.add_argument("--unit-power", type=int, default=2)
    parser.add_argument("--torsion-orders", default="2,3")
    parser.add_argument("--normaliz", default="normaliz")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--summary-output", type=Path)
    arguments = parser.parse_args()

    if arguments.radial_degree < 1:
        parser.error("--radial-degree must be positive")
    if arguments.max_support is not None and arguments.max_support < 5:
        parser.error("--max-support must be at least five")
    if arguments.max_exponent < 1:
        parser.error("--max-exponent must be positive")
    if arguments.max_quotient < 1:
        parser.error("--max-quotient must be positive")
    if arguments.unit_power < 1:
        parser.error("--unit-power must be positive")

    primes = parse_csv_integers(arguments.primes)
    if any(not is_prime(prime) for prime in primes):
        parser.error("--primes must contain primes only")
    torsion_orders = parse_csv_integers(arguments.torsion_orders)
    if any(order < 2 for order in torsion_orders):
        parser.error("--torsion-orders must be at least two")
    normaliz = shutil.which(arguments.normaliz)
    if normaliz is None:
        parser.error(f"Normaliz executable not found: {arguments.normaliz}")

    result = build_experiment(
        radial_degree=arguments.radial_degree,
        max_support=arguments.max_support,
        primitive_only=arguments.primitive_only,
        primes=primes,
        max_exponent=arguments.max_exponent,
        max_quotient=arguments.max_quotient,
        unit_power=arguments.unit_power,
        torsion_orders=torsion_orders,
        normaliz=normaliz,
    )
    if arguments.output is not None:
        write_json(arguments.output, result)
    if arguments.summary_output is not None:
        write_json(arguments.summary_output, compact_summary(result))

    summary = result["summary"]
    print(f"Normaliz: {result['parameters']['normaliz']}")
    print(
        "Graver census: "
        f"raw={summary['raw_full_support_graver_basis']}, "
        f"normalized={summary['symmetry_normalized_graver_basis']}, "
        f"candidate packets={summary['candidate_packets']}, "
        f"support semigroups={summary['candidate_nonfree_support_semigroups']}"
    )
    print(
        "collision funnel: "
        "valuation="
        f"{summary['configured_total_valuation_collisions']}, "
        "finite adelic="
        f"{summary['finite_adelic_weight_collisions']}, "
        "finite false="
        f"{summary['finite_adelic_false_collisions']}, "
        f"scalar-factorial={summary['exact_scalar_factorial_collisions']}, "
        f"marked-partition={summary['exact_marked_partition_collisions']}, "
        f"within-relation full={summary['within_relation_full_probe_collisions']}, "
        "exact all-scale configured="
        f"{summary['within_relation_exact_all_scale_configured_collisions']}, "
        "global full groups="
        f"{summary['global_distinct_state_full_probe_collision_groups']}, "
        "not pairwise support-four="
        f"{summary['global_collision_groups_not_pairwise_support_four']}"
    )
    print(f"result sha256: {result['result_sha256']}")
    if arguments.output is not None:
        print(f"wrote {arguments.output}")
    if arguments.summary_output is not None:
        print(f"wrote {arguments.summary_output}")
    print(
        "STATUS: bounded projected-scroll experiment; not the Hall-envelope "
        "GVC(2) theorem and not a counterexample"
    )


if __name__ == "__main__":
    main()

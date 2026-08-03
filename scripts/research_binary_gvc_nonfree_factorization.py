#!/usr/bin/env python3
"""Second-order factorial tomography for binary projected return fibers.

The first prime-power tomography census studies primitive relations between
two states in a colored partition fiber.  This script studies the next
obstruction: the affine semigroup of *all scaled states* in each such fiber
can have more than one factorization into Hilbert atoms.

For a normalized primitive packet with operator support ``R``, polynomial
support ``B``, color counts ``(c_R,c_B)``, level ``w``, and span ``s``, the
scaled-state semigroup is

    sum x_R = c_R N,
    sum x_B = c_B N,
    sum i*x_i = w N,
    N,x_i >= 0.

Normaliz computes its Hilbert basis.  Each Hilbert atom ``h`` has a
side-refined signed factorial vector encoding

    (c_R N)! / prod(x_R!),
    (c_B N)! / prod(x_B!),
    (w N)! ((s(c_R+c_B)-w)N)!.

The matrix consisting of the Hilbert atoms and these signed vectors detects
factorizations which give both the same semigroup state and the same exact
factorial product at every dilation.  Its Lawrence lifting gives the Graver
basis of the factorial-compatible factorization lattice.  Equality of the
signed vectors is stronger than equality in any finite prime-power scalar
window: the product has the same Legendre valuation and factorial unit at
every prime and every scale.  Individual marked carry positions are retained
as separate packet data and can still distinguish two atoms.

The default span-four census is complete at the Hilbert-basis and exact-rank
levels.  Full Graver bases are computed for collision profiles with at most
``--max-graver-atoms`` atoms.  Larger collision lattices are retained with an
exact integer kernel basis and are explicitly marked as not Graver-complete.

This is an exact bounded theorem about the projected packet model, not the
proof of unrestricted GVC(2).  A surviving factorization relation is a
candidate for the parked Hall-shell inheritance problem, not a binary GVC
counterexample.  The unrestricted proof is the separate Hall-envelope
theorem.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from collections import Counter
from pathlib import Path
from typing import Any

import sympy

from research_binary_gvc_prime_power_tomography import (
    Relation,
    colored_graver_basis,
    normaliz_hilbert_basis,
    normaliz_version,
    primitive_normalization,
    relation_record,
    relation_support_size,
    state_counts,
    state_level,
    support_pair,
)


Vector = tuple[int, ...]
FactorKey = tuple[str, int]


def residue_histogram_matrix(
    span: int,
    orders: tuple[int, ...],
) -> sympy.Matrix:
    """Matrix of all residue-class sums on levels zero through the span."""

    return sympy.Matrix(
        [
            [int(level % order == residue) for level in range(span + 1)]
            for order in orders
            for residue in range(order)
        ]
    )


def verify_consecutive_residue_reconstruction(max_order: int = 16) -> None:
    """Replay the consecutive-modulus theorem through a configurable order."""

    for order in range(2, max_order + 1):
        orders = (order, order + 1)

        odd_span = 2 * order - 1
        odd_matrix = residue_histogram_matrix(odd_span, orders)
        if odd_matrix.rank() != odd_span + 1 or odd_matrix.nullspace():
            raise AssertionError(
                f"C{order},C{order + 1} failed on odd span {odd_span}"
            )

        even_span = 2 * order
        even_matrix = residue_histogram_matrix(even_span, orders)
        kernel = integer_kernel_basis(even_matrix)
        expected = (1,) * order + (0,) + (-1,) * order
        if even_matrix.rank() != even_span or kernel != (expected,):
            raise AssertionError(
                f"unexpected span-{even_span} residue kernel: "
                f"rank={even_matrix.rank()}, kernel={kernel}"
            )
        if sum(expected) != 0 or sum(
            level * value for level, value in enumerate(expected)
        ) != -order * (order + 1):
            raise AssertionError(
                f"incorrect count or weighted level at order {order}"
            )

        safe_sum = [0] * (even_span + 1)
        for level in range(order):
            safe_sum[level] += 1
            safe_sum[level + order + 1] -= 1
        if tuple(safe_sum) != expected:
            raise AssertionError(
                f"the order-{order} kernel did not split into safe swaps"
            )


def scaled_fiber_hilbert_basis(
    relation: Relation,
    normaliz: str,
) -> tuple[tuple[int, ...], tuple[int, ...], tuple[Vector, ...]]:
    """Return the supports and exact Hilbert basis of one scaled fiber."""

    operator_support, polynomial_support = support_pair(relation)
    operator_count, polynomial_count = state_counts(relation.left)
    level = state_level(relation.left)
    equations = (
        (-operator_count,)
        + (1,) * len(operator_support)
        + (0,) * len(polynomial_support),
        (-polynomial_count,)
        + (0,) * len(operator_support)
        + (1,) * len(polynomial_support),
        (-level,) + operator_support + polynomial_support,
    )
    atoms = tuple(sorted(normaliz_hilbert_basis(equations, normaliz)))
    if not atoms:
        raise AssertionError("a primitive packet must have Hilbert atoms")
    return operator_support, polynomial_support, atoms


def add_factor(
    vector: dict[FactorKey, int],
    side: str,
    slope: int,
    multiplicity: int,
) -> None:
    if not slope or not multiplicity:
        return
    key = side, slope
    vector[key] = vector.get(key, 0) + multiplicity
    if not vector[key]:
        del vector[key]


def atom_factorial_vector(
    atom: Vector,
    *,
    operator_width: int,
    operator_count: int,
    polynomial_count: int,
    level: int,
    complementary_level: int,
) -> dict[FactorKey, int]:
    """Canonical side-refined signed factorial slopes of one atom."""

    scale = atom[0]
    operator = atom[1 : 1 + operator_width]
    polynomial = atom[1 + operator_width :]
    vector: dict[FactorKey, int] = {}
    add_factor(vector, "operator", operator_count * scale, 1)
    for part in operator:
        add_factor(vector, "operator", part, -1)
    add_factor(vector, "polynomial", polynomial_count * scale, 1)
    for part in polynomial:
        add_factor(vector, "polynomial", part, -1)
    add_factor(vector, "radial_x", level * scale, 1)
    add_factor(vector, "radial_y", complementary_level * scale, 1)
    return vector


def augmented_signature_matrix(
    atoms: tuple[Vector, ...],
    *,
    operator_width: int,
    operator_count: int,
    polynomial_count: int,
    level: int,
    complementary_level: int,
) -> tuple[sympy.Matrix, tuple[FactorKey, ...], tuple[dict[FactorKey, int], ...]]:
    """Matrix for equal state plus equal side-refined factorial vector."""

    factorial_vectors = tuple(
        atom_factorial_vector(
            atom,
            operator_width=operator_width,
            operator_count=operator_count,
            polynomial_count=polynomial_count,
            level=level,
            complementary_level=complementary_level,
        )
        for atom in atoms
    )
    factor_keys = tuple(
        sorted({key for vector in factorial_vectors for key in vector})
    )
    rows = [
        [atom[row] for atom in atoms]
        for row in range(len(atoms[0]))
    ]
    rows.extend(
        [vector.get(key, 0) for vector in factorial_vectors]
        for key in factor_keys
    )
    return sympy.Matrix(rows), factor_keys, factorial_vectors


def primitive_integer_vector(vector: sympy.Matrix) -> Vector:
    denominator = math.lcm(*(entry.q for entry in vector))
    entries = [int(entry * denominator) for entry in vector]
    divisor = math.gcd(*entries)
    entries = [entry // abs(divisor) for entry in entries]
    first = next(entry for entry in entries if entry)
    if first < 0:
        entries = [-entry for entry in entries]
    return tuple(entries)


def integer_kernel_basis(matrix: sympy.Matrix) -> tuple[Vector, ...]:
    return tuple(
        primitive_integer_vector(vector)
        for vector in matrix.nullspace()
    )


def reversal_permutation(
    atoms: tuple[Vector, ...],
    *,
    operator_support: tuple[int, ...],
    polynomial_support: tuple[int, ...],
    span: int,
    level: int,
    operator_count: int,
    polynomial_count: int,
) -> tuple[int, ...] | None:
    """Internal level reversal, when it preserves this exact profile."""

    capacity = span * (operator_count + polynomial_count)
    if 2 * level != capacity:
        return None
    if tuple(span - value for value in reversed(operator_support)) != operator_support:
        return None
    if (
        tuple(span - value for value in reversed(polynomial_support))
        != polynomial_support
    ):
        return None

    operator_width = len(operator_support)
    atom_index = {atom: index for index, atom in enumerate(atoms)}
    permutation = []
    for atom in atoms:
        reversed_atom = (
            (atom[0],)
            + tuple(reversed(atom[1 : 1 + operator_width]))
            + tuple(reversed(atom[1 + operator_width :]))
        )
        if reversed_atom not in atom_index:
            raise AssertionError("reversal did not preserve the Hilbert basis")
        permutation.append(atom_index[reversed_atom])
    return tuple(permutation)


def reversal_constraint_matrix(permutation: tuple[int, ...]) -> sympy.Matrix:
    """Equations z_i+z_rev(i)=0, with z_i=0 at fixed atoms."""

    rows = []
    seen: set[int] = set()
    for index, image in enumerate(permutation):
        if index in seen:
            continue
        row = [0] * len(permutation)
        row[index] = 1
        if image != index:
            row[image] = 1
        rows.append(row)
        seen.update((index, image))
    return sympy.Matrix(rows)


def kernel_is_reversal_only(
    matrix: sympy.Matrix,
    permutation: tuple[int, ...] | None,
) -> bool:
    if permutation is None or matrix.cols == matrix.rank():
        return False
    constraints = reversal_constraint_matrix(permutation)
    return matrix.col_join(constraints).rank() == matrix.rank()


def orient_factorization(left: Vector, right: Vector) -> tuple[Vector, Vector]:
    return (left, right) if left < right else (right, left)


def factorial_graver_basis(
    matrix: sympy.Matrix,
    normaliz: str,
) -> tuple[tuple[Vector, Vector], ...]:
    """Graver basis of the augmented columns by a Lawrence lifting."""

    rows = tuple(
        tuple(int(matrix[row, column]) for column in range(matrix.cols))
        for row in range(matrix.rows)
    )
    equations = tuple(
        row + tuple(-entry for entry in row)
        for row in rows
    )
    lifted = normaliz_hilbert_basis(equations, normaliz)
    relations: set[tuple[Vector, Vector]] = set()
    for element in lifted:
        left = element[: matrix.cols]
        right = element[matrix.cols :]
        if not any(left) or not any(right):
            continue
        if any(a and b for a, b in zip(left, right, strict=True)):
            continue
        relations.add(orient_factorization(left, right))
    return tuple(
        sorted(
            relations,
            key=lambda pair: (
                max(sum(pair[0]), sum(pair[1])),
                sum(pair[0]) + sum(pair[1]),
                pair,
            ),
        )
    )


def sparse_vector(vector: Vector) -> list[list[int]]:
    return [
        [index, value]
        for index, value in enumerate(vector)
        if value
    ]


def atom_packet_signature(
    atom: Vector,
    *,
    operator_support: tuple[int, ...],
    polynomial_support: tuple[int, ...],
    level: int,
    complementary_level: int,
    torsion_orders: tuple[int, ...],
) -> tuple[Any, ...]:
    """Packet data short of inserting the exact marked atom by definition."""

    operator_width = len(operator_support)
    scale = atom[0]
    operator = atom[1 : 1 + operator_width]
    polynomial = atom[1 + operator_width :]
    signature: list[Any] = [
        (level * scale, complementary_level * scale),
        tuple(sorted((part for part in operator if part), reverse=True)),
        tuple(sorted((part for part in polynomial if part), reverse=True)),
    ]
    for order in torsion_orders:
        operator_histogram = [0] * order
        polynomial_histogram = [0] * order
        for marked_level, multiplicity in zip(
            operator_support,
            operator,
            strict=True,
        ):
            operator_histogram[marked_level % order] += multiplicity
        for marked_level, multiplicity in zip(
            polynomial_support,
            polynomial,
            strict=True,
        ):
            polynomial_histogram[marked_level % order] += multiplicity
        signature.extend(
            (tuple(operator_histogram), tuple(polynomial_histogram))
        )
    return tuple(signature)


def packet_signature_is_injective(
    atoms: tuple[Vector, ...],
    *,
    operator_support: tuple[int, ...],
    polynomial_support: tuple[int, ...],
    level: int,
    complementary_level: int,
    torsion_orders: tuple[int, ...],
) -> bool:
    signatures = {
        atom_packet_signature(
            atom,
            operator_support=operator_support,
            polynomial_support=polynomial_support,
            level=level,
            complementary_level=complementary_level,
            torsion_orders=torsion_orders,
        )
        for atom in atoms
    }
    return len(signatures) == len(atoms)


def factorization_packet_signature(
    side: Vector,
    atoms: tuple[Vector, ...],
    *,
    operator_support: tuple[int, ...],
    polynomial_support: tuple[int, ...],
    level: int,
    complementary_level: int,
    torsion_orders: tuple[int, ...],
) -> tuple[tuple[Any, ...], ...]:
    expanded = []
    for index, multiplicity in enumerate(side):
        signature = atom_packet_signature(
            atoms[index],
            operator_support=operator_support,
            polynomial_support=polynomial_support,
            level=level,
            complementary_level=complementary_level,
            torsion_orders=torsion_orders,
        )
        expanded.extend((signature,) * multiplicity)
    return tuple(sorted(expanded))


def first_packet_separator(
    relation: tuple[Vector, Vector],
    atoms: tuple[Vector, ...],
    *,
    operator_support: tuple[int, ...],
    polynomial_support: tuple[int, ...],
    level: int,
    complementary_level: int,
) -> str | None:
    left, right = relation
    stages = (
        ("packet_partition", ()),
        ("C2", (2,)),
        ("C3", (2, 3)),
    )
    for name, orders in stages:
        left_signature = factorization_packet_signature(
            left,
            atoms,
            operator_support=operator_support,
            polynomial_support=polynomial_support,
            level=level,
            complementary_level=complementary_level,
            torsion_orders=orders,
        )
        right_signature = factorization_packet_signature(
            right,
            atoms,
            operator_support=operator_support,
            polynomial_support=polynomial_support,
            level=level,
            complementary_level=complementary_level,
            torsion_orders=orders,
        )
        if left_signature != right_signature:
            return name
    return None


def graver_record(
    relation: tuple[Vector, Vector],
    reversal: tuple[int, ...] | None,
    atoms: tuple[Vector, ...],
    *,
    operator_support: tuple[int, ...],
    polynomial_support: tuple[int, ...],
    level: int,
    complementary_level: int,
) -> dict[str, Any]:
    left, right = relation
    reversed_left = None
    if reversal is not None:
        reversed_left = tuple(left[reversal[index]] for index in range(len(left)))
    return {
        "left": sparse_vector(left),
        "right": sparse_vector(right),
        "factorization_degrees": [sum(left), sum(right)],
        "reversal_relation": reversed_left == right,
        "first_packet_separator": first_packet_separator(
            relation,
            atoms,
            operator_support=operator_support,
            polynomial_support=polynomial_support,
            level=level,
            complementary_level=complementary_level,
        ),
    }


def factorial_vector_record(vector: dict[FactorKey, int]) -> list[list[Any]]:
    return [
        [side, slope, multiplicity]
        for (side, slope), multiplicity in sorted(vector.items())
    ]


def profile_record(
    profile_id: str,
    relation: Relation,
    *,
    normaliz: str,
    max_graver_atoms: int,
) -> dict[str, Any]:
    operator_support, polynomial_support, atoms = scaled_fiber_hilbert_basis(
        relation,
        normaliz,
    )
    operator_count, polynomial_count = state_counts(relation.left)
    level = state_level(relation.left)
    span = len(relation.left.operator) - 1
    complementary_level = span * (operator_count + polynomial_count) - level
    matrix, factor_keys, factorial_vectors = augmented_signature_matrix(
        atoms,
        operator_width=len(operator_support),
        operator_count=operator_count,
        polynomial_count=polynomial_count,
        level=level,
        complementary_level=complementary_level,
    )
    rank = matrix.rank()
    nullity = matrix.cols - rank
    kernel_basis = integer_kernel_basis(matrix)
    if len(kernel_basis) != nullity:
        raise AssertionError("incorrect exact kernel dimension")

    reversal = reversal_permutation(
        atoms,
        operator_support=operator_support,
        polynomial_support=polynomial_support,
        span=span,
        level=level,
        operator_count=operator_count,
        polynomial_count=polynomial_count,
    )
    reversal_only = kernel_is_reversal_only(matrix, reversal)

    graver_complete = bool(nullity and len(atoms) <= max_graver_atoms)
    graver = (
        factorial_graver_basis(matrix, normaliz)
        if graver_complete
        else ()
    )
    if graver_complete and not graver:
        raise AssertionError("nonzero lattice has an empty Graver basis")
    graver_records = tuple(
        graver_record(
            relation,
            reversal,
            atoms,
            operator_support=operator_support,
            polynomial_support=polynomial_support,
            level=level,
            complementary_level=complementary_level,
        )
        for relation in graver
    )
    if reversal_only and graver_complete:
        if not all(record["reversal_relation"] for record in graver_records):
            raise AssertionError("a reversal-only lattice has a non-reversal move")

    if nullity == 0:
        classification = "factorial_signature_injective"
    elif reversal_only:
        classification = "reversal_only"
    else:
        classification = "same_vector_candidate"

    partition_injective = packet_signature_is_injective(
        atoms,
        operator_support=operator_support,
        polynomial_support=polynomial_support,
        level=level,
        complementary_level=complementary_level,
        torsion_orders=(),
    )
    c2_injective = packet_signature_is_injective(
        atoms,
        operator_support=operator_support,
        polynomial_support=polynomial_support,
        level=level,
        complementary_level=complementary_level,
        torsion_orders=(2,),
    )
    c2_c3_injective = packet_signature_is_injective(
        atoms,
        operator_support=operator_support,
        polynomial_support=polynomial_support,
        level=level,
        complementary_level=complementary_level,
        torsion_orders=(2, 3),
    )
    if partition_injective:
        first_atom_separator = "packet_partition"
    elif c2_injective:
        first_atom_separator = "C2"
    elif c2_c3_injective:
        first_atom_separator = "C3"
    else:
        first_atom_separator = None

    return {
        "id": profile_id,
        "classification": classification,
        "seed_packet": relation_record(relation),
        "operator_support": list(operator_support),
        "polynomial_support": list(polynomial_support),
        "color_counts": [operator_count, polynomial_count],
        "radial_vector": [level, complementary_level],
        "hilbert_basis_size": len(atoms),
        "hilbert_basis": [list(atom) for atom in atoms],
        "factorial_vectors": [
            factorial_vector_record(vector)
            for vector in factorial_vectors
        ],
        "factorial_row_keys": [list(key) for key in factor_keys],
        "augmented_matrix_rank": rank,
        "factorial_collision_lattice_nullity": nullity,
        "integer_kernel_basis": [sparse_vector(vector) for vector in kernel_basis],
        "atom_packet_signature": {
            "partition_injective": partition_injective,
            "C2_injective": c2_injective,
            "C2_C3_injective": c2_c3_injective,
            "first_injective_stage": first_atom_separator,
        },
        "internal_reversal": list(reversal) if reversal is not None else None,
        "collision_lattice_is_reversal_only": reversal_only,
        "factorial_graver_complete": graver_complete,
        "factorial_graver_basis_size": len(graver),
        "factorial_graver_basis": list(graver_records),
    }


def first_candidate_relation(profiles: list[dict[str, Any]]) -> dict[str, Any]:
    candidates = []
    for profile in profiles:
        if profile["classification"] != "same_vector_candidate":
            continue
        if not profile["factorial_graver_complete"]:
            continue
        for relation in profile["factorial_graver_basis"]:
            if relation["reversal_relation"]:
                continue
            candidates.append(
                (
                    max(relation["factorization_degrees"]),
                    sum(relation["factorization_degrees"]),
                    profile["id"],
                    relation,
                    profile,
                )
            )
    if not candidates:
        raise AssertionError("the configured census has no candidate relation")
    _, _, _, relation, profile = min(candidates, key=lambda item: item[:3])
    return {
        "profile_id": profile["id"],
        "operator_support": profile["operator_support"],
        "polynomial_support": profile["polynomial_support"],
        "color_counts": profile["color_counts"],
        "radial_vector": profile["radial_vector"],
        "hilbert_basis": profile["hilbert_basis"],
        "relation": relation,
    }


def logical_digest(result: dict[str, Any]) -> str:
    payload = dict(result)
    payload.pop("result_sha256", None)
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def signature_only_profile_record(
    profile_id: str,
    relation: Relation,
    *,
    normaliz: str,
    torsion_orders: tuple[int, ...],
) -> dict[str, Any]:
    operator_support, polynomial_support, atoms = scaled_fiber_hilbert_basis(
        relation,
        normaliz,
    )
    operator_count, polynomial_count = state_counts(relation.left)
    level = state_level(relation.left)
    span = len(relation.left.operator) - 1
    complementary_level = span * (operator_count + polynomial_count) - level

    stages: list[tuple[str, tuple[int, ...]]] = [("packet_partition", ())]
    prefix: list[int] = []
    for order in torsion_orders:
        prefix.append(order)
        stages.append((f"C{order}", tuple(prefix)))

    stage_records = {}
    first_injective_stage = None
    final_groups: dict[tuple[Any, ...], list[int]] = {}
    for name, orders in stages:
        groups: dict[tuple[Any, ...], list[int]] = {}
        for index, atom in enumerate(atoms):
            signature = atom_packet_signature(
                atom,
                operator_support=operator_support,
                polynomial_support=polynomial_support,
                level=level,
                complementary_level=complementary_level,
                torsion_orders=orders,
            )
            groups.setdefault(signature, []).append(index)
        collisions = tuple(
            tuple(indices)
            for indices in groups.values()
            if len(indices) > 1
        )
        injective = not collisions
        stage_records[name] = injective
        if first_injective_stage is None and injective:
            first_injective_stage = name
        final_groups = groups

    final_collisions = [
        indices
        for indices in final_groups.values()
        if len(indices) > 1
    ]
    return {
        "id": profile_id,
        "seed_packet": relation_record(relation),
        "operator_support": list(operator_support),
        "polynomial_support": list(polynomial_support),
        "color_counts": [operator_count, polynomial_count],
        "radial_vector": [level, complementary_level],
        "hilbert_basis_size": len(atoms),
        "hilbert_basis": [list(atom) for atom in atoms],
        "atom_packet_signature": {
            "injective_by_stage": stage_records,
            "first_injective_stage": first_injective_stage,
            "final_collision_atom_indices": final_collisions,
        },
    }


def build_signature_only_census(
    *,
    radial_degree: int,
    torsion_orders: tuple[int, ...],
    normaliz: str,
) -> dict[str, Any]:
    raw_basis = colored_graver_basis(radial_degree, normaliz)
    normalized_basis = tuple(
        sorted({primitive_normalization(relation) for relation in raw_basis})
    )
    candidates = tuple(
        relation
        for relation in normalized_basis
        if relation_support_size(relation) >= 5
        and all(state_counts(relation.left))
    )
    profile_relations: dict[
        tuple[Any, ...], Relation
    ] = {}
    for relation in candidates:
        key = (
            support_pair(relation),
            state_counts(relation.left),
            state_level(relation.left),
            len(relation.left.operator) - 1,
        )
        profile_relations.setdefault(key, relation)

    profiles = [
        signature_only_profile_record(
            f"F{index:04d}",
            relation,
            normaliz=normaliz,
            torsion_orders=torsion_orders,
        )
        for index, relation in enumerate(profile_relations.values(), start=1)
    ]
    injectivity = Counter(
        profile["atom_packet_signature"]["first_injective_stage"]
        or "unresolved"
        for profile in profiles
    )
    unresolved = [
        profile
        for profile in profiles
        if profile["atom_packet_signature"]["first_injective_stage"] is None
    ]
    result: dict[str, Any] = {
        "schema": "binary-gvc-nonfree-factorization-signature-v1",
        "status": {
            "kind": "exact bounded computation",
            "claim": (
                "Complete Hilbert-basis and cumulative packet-signature "
                "injectivity census in the configured projected span. "
                "It omits the factorial collision ranks and Graver bases."
            ),
        },
        "parameters": {
            "radial_degree": radial_degree,
            "torsion_orders": list(torsion_orders),
            "normaliz": normaliz_version(normaliz),
        },
        "summary": {
            "raw_projected_graver_basis": len(raw_basis),
            "symmetry_normalized_graver_basis": len(normalized_basis),
            "candidate_packets": len(candidates),
            "candidate_profiles": len(profiles),
            "atom_signature_injectivity": dict(sorted(injectivity.items())),
            "maximum_hilbert_basis_size": max(
                profile["hilbert_basis_size"] for profile in profiles
            ),
            "all_signature_collision_profiles": len(unresolved),
        },
        "unresolved_profiles": [profile["id"] for profile in unresolved],
        "profiles": profiles,
    }

    if radial_degree == 5 and torsion_orders == (2, 3, 4):
        expected = {
            "raw_projected_graver_basis": 2225,
            "symmetry_normalized_graver_basis": 460,
            "candidate_packets": 404,
            "candidate_profiles": 400,
            "atom_signature_injectivity": {
                "C2": 80,
                "C3": 143,
                "C4": 8,
                "packet_partition": 169,
            },
            "maximum_hilbert_basis_size": 155,
            "all_signature_collision_profiles": 0,
        }
        for key, value in expected.items():
            if result["summary"][key] != value:
                raise AssertionError(
                    f"span-five signature census changed at {key}: "
                    f"{result['summary'][key]} != {value}"
                )

    if radial_degree == 6 and torsion_orders == (2, 3, 4):
        expected = {
            "raw_projected_graver_basis": 8559,
            "symmetry_normalized_graver_basis": 1584,
            "candidate_packets": 1490,
            "candidate_profiles": 1469,
            "atom_signature_injectivity": {
                "C2": 358,
                "C3": 599,
                "C4": 130,
                "packet_partition": 382,
            },
            "maximum_hilbert_basis_size": 445,
            "all_signature_collision_profiles": 0,
        }
        for key, value in expected.items():
            if result["summary"][key] != value:
                raise AssertionError(
                    f"span-six signature census changed at {key}: "
                    f"{result['summary'][key]} != {value}"
                )

    result["result_sha256"] = logical_digest(result)
    return result


def build_census(
    *,
    radial_degree: int,
    max_graver_atoms: int,
    normaliz: str,
) -> dict[str, Any]:
    raw_basis = colored_graver_basis(radial_degree, normaliz)
    normalized_basis = tuple(
        sorted({primitive_normalization(relation) for relation in raw_basis})
    )
    candidates = tuple(
        relation
        for relation in normalized_basis
        if relation_support_size(relation) >= 5
        and all(state_counts(relation.left))
    )
    profile_keys = {
        (
            support_pair(relation),
            state_counts(relation.left),
            state_level(relation.left),
            len(relation.left.operator) - 1,
        )
        for relation in candidates
    }
    if len(profile_keys) != len(candidates):
        raise AssertionError(
            "the configured range has duplicate normalized fiber profiles"
        )

    profiles = [
        profile_record(
            f"F{index:04d}",
            relation,
            normaliz=normaliz,
            max_graver_atoms=max_graver_atoms,
        )
        for index, relation in enumerate(candidates, start=1)
    ]
    classifications = Counter(
        profile["classification"]
        for profile in profiles
    )
    graver_profiles = [
        profile for profile in profiles if profile["factorial_graver_complete"]
    ]
    graver_relations = sum(
        profile["factorial_graver_basis_size"]
        for profile in graver_profiles
    )
    graver_reversal_relations = sum(
        relation["reversal_relation"]
        for profile in graver_profiles
        for relation in profile["factorial_graver_basis"]
    )
    result: dict[str, Any] = {
        "schema": "binary-gvc-nonfree-factorization-v1",
        "status": {
            "kind": "exact bounded computation",
            "claim": (
                "Complete Hilbert-basis and side-refined factorial-rank "
                "census in the configured projected span; Graver completeness "
                "is recorded profile by profile. Survivors are Hall-promotion "
                "candidates, not GVC counterexamples."
            ),
        },
        "parameters": {
            "radial_degree": radial_degree,
            "max_graver_atoms": max_graver_atoms,
            "normaliz": normaliz_version(normaliz),
            "sympy": sympy.__version__,
        },
        "model": {
            "candidate_filter": (
                "symmetry-normalized two-color Graver packets with mixed "
                "colors and support at least five"
            ),
            "factorial_signature": (
                "exact side-refined signed slopes for operator multinomial, "
                "polynomial multinomial, and the two radial factorials"
            ),
            "safe_symmetry": "internal level reversal",
        },
        "summary": {
            "raw_projected_graver_basis": len(raw_basis),
            "symmetry_normalized_graver_basis": len(normalized_basis),
            "candidate_profiles": len(candidates),
            "classification": dict(sorted(classifications.items())),
            "factorial_graver_complete_profiles": len(graver_profiles),
            "factorial_graver_deferred_profiles": sum(
                profile["factorial_collision_lattice_nullity"] > 0
                and not profile["factorial_graver_complete"]
                for profile in profiles
            ),
            "computed_factorial_graver_relations": graver_relations,
            "computed_reversal_relations": graver_reversal_relations,
            "computed_nonreversal_relations": (
                graver_relations - graver_reversal_relations
            ),
        },
        "first_factorial_only_candidate": first_candidate_relation(profiles),
        "profiles": profiles,
    }

    atom_injectivity = Counter(
        profile["atom_packet_signature"]["first_injective_stage"]
        or "unresolved"
        for profile in profiles
    )
    collision_atom_injectivity = Counter(
        profile["atom_packet_signature"]["first_injective_stage"]
        or "unresolved"
        for profile in profiles
        if profile["factorial_collision_lattice_nullity"]
    )
    primitive_separators = Counter(
        relation["first_packet_separator"] or "unresolved"
        for profile in graver_profiles
        for relation in profile["factorial_graver_basis"]
    )
    result["summary"].update(
        {
            "atom_signature_injectivity": dict(sorted(atom_injectivity.items())),
            "factorial_collision_atom_injectivity": dict(
                sorted(collision_atom_injectivity.items())
            ),
            "computed_primitive_packet_separators": dict(
                sorted(primitive_separators.items())
            ),
            "all_profiles_C2_C3_atom_injective": all(
                profile["atom_packet_signature"]["C2_C3_injective"]
                for profile in profiles
            ),
            "all_signature_collision_profiles": sum(
                not profile["atom_packet_signature"]["C2_C3_injective"]
                for profile in profiles
            ),
        }
    )

    if radial_degree == 4 and max_graver_atoms >= 20:
        expected = {
            "raw_projected_graver_basis": 426,
            "symmetry_normalized_graver_basis": 90,
            "candidate_profiles": 65,
            "classification": {
                "factorial_signature_injective": 52,
                "reversal_only": 1,
                "same_vector_candidate": 12,
            },
            "factorial_graver_complete_profiles": 11,
            "factorial_graver_deferred_profiles": 2,
            "atom_signature_injectivity": {
                "C2": 14,
                "C3": 16,
                "packet_partition": 35,
            },
            "factorial_collision_atom_injectivity": {"C3": 13},
            "computed_primitive_packet_separators": {
                "C2": 207,
                "C3": 86,
                "packet_partition": 15,
            },
            "all_profiles_C2_C3_atom_injective": True,
            "all_signature_collision_profiles": 0,
        }
        for key, value in expected.items():
            if result["summary"][key] != value:
                raise AssertionError(
                    f"span-four census changed at {key}: "
                    f"{result['summary'][key]} != {value}"
                )
        first = result["first_factorial_only_candidate"]
        if (
            first["operator_support"] != [0, 4]
            or first["polynomial_support"] != [0, 1, 2, 4]
            or first["color_counts"] != [1, 3]
            or first["radial_vector"] != [8, 8]
            or first["relation"]["factorization_degrees"] != [2, 2]
        ):
            raise AssertionError("the first span-four collision candidate changed")

    result["result_sha256"] = logical_digest(result)
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--radial-degree", type=int, default=4)
    parser.add_argument("--max-graver-atoms", type=int, default=20)
    parser.add_argument("--normaliz", default="normaliz")
    parser.add_argument(
        "--signature-only",
        action="store_true",
        help="compute all Hilbert atoms and packet-signature injectivity only",
    )
    parser.add_argument(
        "--verify-consecutive-residues",
        action="store_true",
        help=(
            "replay the all-span consecutive-residue theorem through "
            "cyclic order 16, without invoking Normaliz"
        ),
    )
    parser.add_argument(
        "--torsion-orders",
        default="2,3",
        help="cumulative torsion orders used by --signature-only",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "artifacts/generated-results/"
            "binary_gvc_nonfree_factorization_span4.json"
        ),
    )
    arguments = parser.parse_args()
    if arguments.radial_degree < 1:
        parser.error("--radial-degree must be positive")
    if arguments.max_graver_atoms < 1:
        parser.error("--max-graver-atoms must be positive")

    torsion_orders = tuple(
        int(value)
        for value in arguments.torsion_orders.split(",")
        if value
    )
    if not torsion_orders or any(order < 2 for order in torsion_orders):
        parser.error("--torsion-orders must list integers at least two")
    if len(set(torsion_orders)) != len(torsion_orders):
        parser.error("--torsion-orders must not contain duplicates")

    verify_consecutive_residue_reconstruction()
    if arguments.verify_consecutive_residues:
        print(
            "PASS consecutive-residue reconstruction through C16,C17: "
            "odd spans injective; even kernels split into safe beta swaps"
        )
        print(
            "STATUS: bounded exact regression for the all-order incidence-"
            "forest proof in the canonical note"
        )
        return

    if arguments.signature_only:
        result = build_signature_only_census(
            radial_degree=arguments.radial_degree,
            torsion_orders=torsion_orders,
            normaliz=arguments.normaliz,
        )
    else:
        result = build_census(
            radial_degree=arguments.radial_degree,
            max_graver_atoms=arguments.max_graver_atoms,
            normaliz=arguments.normaliz,
        )
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    summary = result["summary"]
    if arguments.signature_only:
        print(
            "PASS span-"
            f"{arguments.radial_degree} Hilbert-atom signature census: "
            f"{summary['candidate_profiles']} profiles; "
            f"{summary['atom_signature_injectivity']}"
        )
        print(
            "all-signature collision profiles: "
            f"{summary['all_signature_collision_profiles']}"
        )
        print(f"logical result sha256: {result['result_sha256']}")
        print(
            "STATUS: exact bounded projected-semigroup computation; "
            "factorial collision ranks were intentionally omitted"
        )
        return
    print(
        "PASS span-"
        f"{arguments.radial_degree} nonfree factorization census: "
        f"{summary['candidate_profiles']} profiles; "
        f"{summary['classification']}"
    )
    print(
        "factorial-compatible Graver bases: "
        f"{summary['factorial_graver_complete_profiles']} complete, "
        f"{summary['factorial_graver_deferred_profiles']} deferred; "
        f"{summary['computed_factorial_graver_relations']} primitive moves"
    )
    first = result["first_factorial_only_candidate"]
    print(
        "first factorial-only candidate: "
        f"R={first['operator_support']}, B={first['polynomial_support']}, "
        f"counts={first['color_counts']}, radial={first['radial_vector']}, "
        f"degrees={first['relation']['factorization_degrees']}"
    )
    print(
        "packet signature: all Hilbert atoms are injective after C2,C3; "
        "no all-signature factorization collision in the configured span"
    )
    print(f"logical result sha256: {result['result_sha256']}")
    print(
        "STATUS: exact bounded projected-semigroup computation; "
        "affine Hall-shell promotion remains unproved in the parked route, "
        "while Hall-envelope separation proves unrestricted GVC(2)"
    )


if __name__ == "__main__":
    main()

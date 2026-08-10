#!/usr/bin/env python3
"""Compile the exact global-attachment obligations for F2 ``(75,125)``.

The terminal residue theorem supplies one explicit boundary row

    source ray (12,-17) -> target ray (5,2),   (e,f)=(1,6),

but it does not supply a completion of either affine plane.  This checker
separates consequences of that certified row from data that only a global
completion can provide.  In particular it

* orients the two target nodes in the regular fan around ``(5,2)``;
* compiles all endpoint and interior attachment slots;
* derives the source-boundary tree lower bounds;
* writes the two live finite-normalization degree equations and retains the
  distinct-target double row as an explicitly excluded audit case; and
* emits a machine-readable contract for the missing class-group, unit,
  canonical, purity, spectator, and global-meridian ledgers.

The pinned default artifact is an exact *incomplete attachment audit*.  It is
not a constructed global completion and not an exclusion of ``(75,125)``.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
import hashlib
import json
from math import gcd
from pathlib import Path
from typing import Sequence


ROOT = Path(__file__).resolve().parents[2]
ARTIFACT = (
    ROOT
    / "artifacts/generated-results/jc2_f2_75_125_global_attachment.json"
)


def determinant(left: tuple[int, int], right: tuple[int, int]) -> int:
    return left[0] * right[1] - left[1] * right[0]


def bareiss_determinant(matrix: Sequence[Sequence[int]]) -> int:
    """Return an exact determinant using fraction-free elimination."""

    size = len(matrix)
    if size == 0:
        return 1
    if any(len(row) != size for row in matrix):
        raise ValueError("determinant requires a square matrix")
    work = [list(map(int, row)) for row in matrix]
    sign = 1
    previous = 1
    for column in range(size - 1):
        pivot = next(
            (row for row in range(column, size) if work[row][column]),
            None,
        )
        if pivot is None:
            return 0
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            sign *= -1
        value = work[column][column]
        for row in range(column + 1, size):
            for inner in range(column + 1, size):
                numerator = (
                    work[row][inner] * value
                    - work[row][column] * work[column][inner]
                )
                if numerator % previous:
                    raise AssertionError("Bareiss division ceased to be exact")
                work[row][inner] = numerator // previous
            work[row][column] = 0
        previous = value
    return sign * work[-1][-1]


def solve_linear_system(
    matrix: Sequence[Sequence[int]],
    vector: Sequence[int],
) -> tuple[Fraction, ...]:
    """Solve a nonsingular square system over the rationals."""

    size = len(matrix)
    if len(vector) != size or any(len(row) != size for row in matrix):
        raise ValueError("linear system dimensions do not match")
    augmented = [
        [Fraction(value) for value in row] + [Fraction(vector[index])]
        for index, row in enumerate(matrix)
    ]
    for column in range(size):
        pivot = next(
            (row for row in range(column, size) if augmented[row][column]),
            None,
        )
        if pivot is None:
            raise ValueError("linear system is singular")
        augmented[column], augmented[pivot] = (
            augmented[pivot],
            augmented[column],
        )
        scale = augmented[column][column]
        augmented[column] = [value / scale for value in augmented[column]]
        for row in range(size):
            if row == column:
                continue
            scale = augmented[row][column]
            if scale:
                augmented[row] = [
                    left - scale * right
                    for left, right in zip(
                        augmented[row], augmented[column]
                    )
                ]
    return tuple(augmented[index][-1] for index in range(size))


def symmetric_inertia(
    matrix: Sequence[Sequence[int | Fraction]],
) -> tuple[int, int, int]:
    """Return positive, negative, and zero inertia by exact congruence."""

    work = [[Fraction(value) for value in row] for row in matrix]
    if any(len(row) != len(work) for row in work):
        raise ValueError("inertia requires a square matrix")
    positive = negative = zero = 0
    while work:
        size = len(work)
        diagonal = next(
            (index for index in range(size) if work[index][index]),
            None,
        )
        if diagonal is not None:
            if diagonal:
                work[0], work[diagonal] = work[diagonal], work[0]
                for row in work:
                    row[0], row[diagonal] = row[diagonal], row[0]
            pivot = work[0][0]
            if pivot > 0:
                positive += 1
            else:
                negative += 1
            work = [
                [
                    work[row][column]
                    - work[row][0] * work[0][column] / pivot
                    for column in range(1, size)
                ]
                for row in range(1, size)
            ]
            continue

        off_diagonal: tuple[int, int] | None = None
        for row in range(size):
            for column in range(row + 1, size):
                if work[row][column]:
                    off_diagonal = (row, column)
                    break
            if off_diagonal is not None:
                break
        if off_diagonal is None:
            zero += size
            break

        first, second = off_diagonal
        order = [first, second, *(
            index for index in range(size) if index not in (first, second)
        )]
        work = [[work[row][column] for column in order] for row in order]
        block = work[0][1]
        # The leading block [[0,block],[block,0]] has inertia (1,1).
        positive += 1
        negative += 1
        inverse = (
            (Fraction(0), Fraction(1, 1) / block),
            (Fraction(1, 1) / block, Fraction(0)),
        )
        reduced: list[list[Fraction]] = []
        for row in range(2, size):
            reduced_row: list[Fraction] = []
            for column in range(2, size):
                correction = sum(
                    work[row][left]
                    * inverse[left][right]
                    * work[right][column]
                    for left in range(2)
                    for right in range(2)
                )
                reduced_row.append(work[row][column] - correction)
            reduced.append(reduced_row)
        work = reduced
    return positive, negative, zero


def regular_subdivision(
    target_coordinates: tuple[int, int],
    ambient_basis: tuple[tuple[int, int], tuple[int, int]],
    *,
    left_boundary: bool,
    right_boundary: bool,
) -> dict[str, object]:
    """Build the minimal Stern--Brocot fan containing one primitive ray.

    The coordinate cone has endpoints ``(1,0),(0,1)``.  Boundary endpoints
    start with self-intersection one, as for a projective-plane boundary
    line; inserted exceptional rays start with self-intersection minus one.
    The nonboundary endpoint models the closure of an affine coordinate line.
    """

    if target_coordinates[0] <= 0 or target_coordinates[1] <= 0:
        raise ValueError("the target ray must lie in the open coordinate cone")
    if gcd(*target_coordinates) != 1:
        raise ValueError("the target ray must be primitive")
    if determinant(*ambient_basis) != 1:
        raise ValueError("the ambient cone basis must be positively regular")
    left = (1, 0)
    right = (0, 1)
    weights: dict[tuple[int, int], int | None] = {
        left: 1 if left_boundary else None,
        right: 1 if right_boundary else None,
    }
    insertion_order: list[tuple[int, int]] = []
    while True:
        middle = (left[0] + right[0], left[1] + right[1])
        insertion_order.append(middle)
        for endpoint in (left, right):
            if weights[endpoint] is not None:
                weights[endpoint] -= 1
        weights[middle] = -1
        if middle == target_coordinates:
            break
        if (
            target_coordinates[1] * middle[0]
            < middle[1] * target_coordinates[0]
        ):
            right = middle
        else:
            left = middle

    def ambient_ray(coordinates: tuple[int, int]) -> tuple[int, int]:
        return (
            coordinates[0] * ambient_basis[0][0]
            + coordinates[1] * ambient_basis[1][0],
            coordinates[0] * ambient_basis[0][1]
            + coordinates[1] * ambient_basis[1][1],
        )

    ordered = sorted(
        weights,
        key=lambda ray: (
            Fraction(ray[1], ray[0]) if ray[0] else Fraction(10**9)
        ),
    )
    if any(determinant(ordered[index], ordered[index + 1]) != 1 for index in range(len(ordered) - 1)):
        raise AssertionError("the Stern--Brocot subdivision is not regular")
    boundary_coordinates = [ray for ray in ordered if weights[ray] is not None]
    boundary_weights = [int(weights[ray]) for ray in boundary_coordinates]
    matrix = [
        [
            boundary_weights[row]
            if row == column
            else int(abs(row - column) == 1)
            for column in range(len(boundary_coordinates))
        ]
        for row in range(len(boundary_coordinates))
    ]
    matrix_determinant = bareiss_determinant(matrix)
    inertia = symmetric_inertia(matrix)
    canonical = solve_linear_system(
        matrix,
        [-2 - weight for weight in boundary_weights],
    )
    if (
        abs(matrix_determinant) != 1
        or inertia != (1, len(matrix) - 1, 0)
        or any(
        value.denominator != 1 for value in canonical
        )
    ):
        raise AssertionError("the minimal fan lost its A2 boundary ledger")
    target_ambient = ambient_ray(target_coordinates)
    target_index = boundary_coordinates.index(target_coordinates)
    return {
        "ambient_basis": [list(ray) for ray in ambient_basis],
        "target_coordinates_in_basis": list(target_coordinates),
        "target_ambient_ray": list(target_ambient),
        "insertion_order_coordinates": [list(ray) for ray in insertion_order],
        "ordered_fan_coordinates": [list(ray) for ray in ordered],
        "ordered_fan_ambient_rays": [list(ambient_ray(ray)) for ray in ordered],
        "boundary_component_coordinates": [
            list(ray) for ray in boundary_coordinates
        ],
        "boundary_component_ambient_rays": [
            list(ambient_ray(ray)) for ray in boundary_coordinates
        ],
        "boundary_self_intersections": boundary_weights,
        "target_boundary_component_index": target_index,
        "target_self_intersection": boundary_weights[target_index],
        "boundary_intersection_matrix": matrix,
        "boundary_intersection_determinant": matrix_determinant,
        "boundary_intersection_inertia_positive_negative_zero": list(inertia),
        "canonical_coefficients": [int(value) for value in canonical],
    }


def source_proximity_fan_audit() -> dict[str, object]:
    """Recover the original nonmonomial valuation and its two-stage fan."""

    normal_t_z = (12, -17)
    # The terminal ``y`` is y_tr=y_old-X^-1, not the original affine y_old.
    # Since t=X*y_tr=X*y_old-1 and z=y_tr^-1=X/t, one has X=t*z.
    order_X = normal_t_z[0] + normal_t_z[1]
    order_y_translated = -normal_t_z[1]
    order_X_inverse = -order_X
    order_y_old = min(order_X_inverse, order_y_translated)
    order_x = 5 * order_X
    order_v_minus_one = normal_t_z[0]
    if (
        order_X,
        order_y_translated,
        order_y_old,
        order_x,
        order_v_minus_one,
    ) != (-5, 17, 5, -25, 12):
        raise AssertionError("the original-coordinate valuation changed")

    # The monomial carrier is the normalized ray (-5,1).  In the P2 cone
    # from affine y_old=0 to the line at infinity it is
    # (-5,1)=6*(0,1)+5*(-1,-1).
    carrier_ray = (-5, 1)
    left = (0, 1)
    infinity = (-1, -1)
    carrier_coordinates = (6, 5)
    if (
        carrier_coordinates[0] * left[0]
        + carrier_coordinates[1] * infinity[0],
        carrier_coordinates[0] * left[1]
        + carrier_coordinates[1] * infinity[1],
    ) != carrier_ray:
        raise AssertionError("the carrier ray is not in the declared P2 cone")
    carrier_fan = regular_subdivision(
        carrier_coordinates,
        (left, infinity),
        left_boundary=False,
        right_boundary=True,
    )
    if len(carrier_fan["insertion_order_coordinates"]) != 6:
        raise AssertionError("the source carrier blowup length changed")

    # At the carrier point v=x*y_old^5=1, q=y_old is transverse and
    # r=v-1 is tangential.  Their terminal orders are (5,12).
    principal_arm = regular_subdivision(
        (5, 12),
        ((1, 0), (0, 1)),
        left_boundary=True,
        right_boundary=False,
    )
    if len(principal_arm["insertion_order_coordinates"]) != 6:
        raise AssertionError("the principal proximity-arm length changed")
    if principal_arm["target_self_intersection"] != -1:
        raise AssertionError("the terminal principal divisor is no longer a -1 curve")

    return {
        "terminal_chart_coordinates": {
            "t": "X*y_translated=X*y_old-1",
            "z": "y_translated^-1=X/t",
            "y_translated": "y_old-X^-1",
        },
        "orders_on_kummer_chart": {
            "X": order_X,
            "y_translated": order_y_translated,
            "t": normal_t_z[0],
            "z": normal_t_z[1],
        },
        "orders_on_original_source": {
            "x": order_x,
            "y_old": order_y_old,
            "v_minus_1_where_v=x*y_old^5": order_v_minus_one,
        },
        "value_group_gcd": gcd(
            abs(order_x), order_y_old, order_v_minus_one
        ),
        "monomial_carrier_ray": list(carrier_ray),
        "carrier_center_coordinate": "v=x*y_old^5",
        "fixed_principal_center": "v=1",
        "standard_P2_carrier_fan": carrier_fan,
        "principal_center_local_orders_q_r": [5, 12],
        "principal_center_local_coordinates": [
            "q=y_old (carrier transverse)",
            "r=v-center (carrier tangential)",
        ],
        "terminal_residue_initial_monomial": (
            "s_terminal=unit*r^5/q^12; its orders on adjacent rays "
            "(3,7),(2,5) are -1,+1"
        ),
        "terminal_endpoint_orientation": {
            "s=infinity": "adjacent local ray (3,7)",
            "s=0": "adjacent local ray (2,5)",
        },
        "principal_arm_fan": principal_arm,
        "exact_consequence": (
            "six blowups extract the carrier ray (-5,1), and six further "
            "blowups at v=1 extract the terminal valuation with "
            "nu(x),nu(y_old),nu(v-1)=(-25,5,12)"
        ),
        "cofactor_centers": {
            "squarefree": (
                "the two simple roots rho_1,rho_2 of R mark two distinct "
                "spectator points v=rho_i on the carrier"
            ),
            "double_root": (
                "the nonzero double root rho != 1 marks a second principal "
                "center v=rho with the same local (5,12) arm"
            ),
        },
    }


def minimal_source_principal_boundary(case: str) -> dict[str, object]:
    """Compile the minimal source graph forced by principal packets."""

    if case not in CASES:
        raise ValueError(f"unknown F2 attachment case: {case}")
    packet_count = PACKET_COUNTS[case]
    proximity = source_proximity_fan_audit()
    carrier_fan = proximity["standard_P2_carrier_fan"]
    carrier_weights = list(carrier_fan["boundary_self_intersections"])
    carrier_index = int(carrier_fan["target_boundary_component_index"])
    names = [f"carrier_chain_{index}" for index in range(len(carrier_weights))]
    weights = dict(zip(names, carrier_weights))
    edges = [
        (names[index], names[index + 1])
        for index in range(len(names) - 1)
    ]
    carrier_name = names[carrier_index]

    arm_fan = proximity["principal_arm_fan"]
    arm_weights = list(arm_fan["boundary_self_intersections"])
    if int(arm_fan["target_boundary_component_index"]) != 4:
        raise AssertionError("the principal-arm terminal index changed")
    # The first local boundary component is the pre-existing carrier.  Its
    # local weight changes from 1 to 0, hence every arm lowers carrier^2 by 1.
    if arm_weights[0] != 0:
        raise AssertionError("the principal arm no longer starts with one carrier blowup")
    new_arm_weights = arm_weights[1:]
    terminal_components: list[str] = []
    packet_centers: dict[str, str] = {}
    attachment_neighbors: dict[str, dict[str, str]] = {}
    for packet_index in range(1, packet_count + 1):
        packet = f"packet_{packet_index}"
        center = "1" if packet_index == 1 else "rho_double"
        packet_centers[packet] = center
        weights[carrier_name] -= 1
        previous = carrier_name
        arm_names: list[str] = []
        for index, weight in enumerate(new_arm_weights, start=1):
            name = f"{packet}_arm_{index}"
            names.append(name)
            arm_names.append(name)
            weights[name] = weight
            edges.append((previous, name))
            previous = name
        terminal_name = arm_names[3]
        terminal_components.append(terminal_name)
        weights[terminal_name] -= 3
        marked: dict[str, str] = {
            "endpoint_s=0": arm_names[4],
            "endpoint_s=infinity": arm_names[2],
        }
        for attachment_index, slot in enumerate(REQUIRED_SLOT_NAMES[2:], start=1):
            leaf = f"{packet}_attachment_{attachment_index}"
            names.append(leaf)
            weights[leaf] = -1
            edges.append((terminal_name, leaf))
            marked[slot] = leaf
        attachment_neighbors[packet] = marked

    edge_set = {tuple(sorted(edge)) for edge in edges}
    matrix = [
        [
            weights[left]
            if left == right
            else int(tuple(sorted((left, right))) in edge_set)
            for right in names
        ]
        for left in names
    ]
    matrix_determinant = bareiss_determinant(matrix)
    inertia = symmetric_inertia(matrix)
    canonical = solve_linear_system(
        matrix,
        [-2 - weights[name] for name in names],
    )
    adjacency = {name: set() for name in names}
    for left_name, right_name in edges:
        adjacency[left_name].add(right_name)
        adjacency[right_name].add(left_name)
    leaf_count = sum(len(adjacency[name]) == 1 for name in names)
    expected = (16, 6) if packet_count == 1 else (25, 10)
    if (len(names), leaf_count) != expected:
        raise AssertionError("the minimal source principal topology changed")
    if abs(matrix_determinant) != 1 or inertia != (1, len(names) - 1, 0):
        raise AssertionError("the source principal boundary lost its A2 ledger")
    if any(value.denominator != 1 for value in canonical):
        raise AssertionError("the source principal canonical class is not integral")
    if any(len(adjacency[name]) != 5 for name in terminal_components):
        raise AssertionError("a terminal component lost its five marked neighbors")
    return {
        "case": case,
        "component_order": names,
        "self_intersections_in_component_order": [weights[name] for name in names],
        "edges": [list(edge) for edge in edges],
        "carrier_component": carrier_name,
        "carrier_valency": len(adjacency[carrier_name]),
        "terminal_components": terminal_components,
        "terminal_valencies": [len(adjacency[name]) for name in terminal_components],
        "packet_centers_on_carrier_v": packet_centers,
        "attachment_neighbors": attachment_neighbors,
        "component_count": len(names),
        "leaf_count": leaf_count,
        "intersection_matrix": matrix,
        "intersection_determinant": matrix_determinant,
        "intersection_inertia_positive_negative_zero": list(inertia),
        "canonical_coefficients_in_component_order": [
            int(value) for value in canonical
        ],
        "unresolved_additions": (
            ["two spectator branches at v=rho_1,rho_2", "purity/other map-resolution branches"]
            if case == "squarefree_one_packet"
            else ["purity/other map-resolution branches"]
        ),
    }


def target_minimal_completion_audit() -> dict[str, object]:
    """Compile the canonical P2 blowup chain extracting target ray (5,2)."""

    fan = regular_subdivision(
        (5, 2),
        ((1, 0), (0, 1)),
        left_boundary=True,
        right_boundary=False,
    )
    if len(fan["insertion_order_coordinates"]) != 4:
        raise AssertionError("the target extraction length changed")
    if fan["boundary_self_intersections"] != [-2, -2, -1, -3, -2]:
        raise AssertionError("the target boundary weights changed")
    return {
        "local_coordinates": ["a=(-Q)^-1", "b=P/(-Q)"],
        "fan": fan,
        "exact_consequence": (
            "four point blowups extract (5,2); together with the transformed "
            "line at infinity the target boundary chain has weights "
            "(-2,-2,-1,-3,-2), determinant one, and integral canonical "
            "coefficients (-3,-6,-9,-4,-2)"
        ),
    }


def target_valuation_uniqueness_audit() -> dict[str, object]:
    """Prove that all principal packets restrict to one target divisor."""

    order_a = 5
    order_b = 2
    order_pi = 3 * order_b - order_a
    order_eta = 2 * order_a - 5 * order_b
    if (gcd(order_a, order_b), order_pi, order_eta) != (1, 1, 0):
        raise AssertionError("the target valuation coordinates changed")
    return {
        "global_target_coordinates": {
            "a": "(-Q)^-1",
            "b": "P/(-Q)",
        },
        "orders_a_b": [order_a, order_b],
        "uniformizer": "pi=b^3/a",
        "uniformizer_order": order_pi,
        "residue_coordinate": "eta=a^2/b^5",
        "residue_coordinate_order": order_eta,
        "residue_pullback": "eta^-1=h(s) is nonconstant of degree 6",
        "conclusion": (
            "every principal source valuation is centered at the generic "
            "point of the unique target divisor for ray (5,2)"
        ),
        "double_root_consequence": (
            "the two principal packets necessarily lie over the same target "
            "divisor; the distinct-target case is excluded"
        ),
    }


def multiply_permutations(
    left: Sequence[int], right: Sequence[int]
) -> tuple[int, ...]:
    """Compose image-list permutations as ``left after right``."""

    if len(left) != len(right):
        raise ValueError("permutation degrees differ")
    return tuple(left[right[index] - 1] for index in range(len(left)))


def permutation_is_valid(permutation: Sequence[int]) -> bool:
    return sorted(permutation) == list(range(1, len(permutation) + 1))


def permutation_orbit(
    generators: Sequence[Sequence[int]], start: int
) -> frozenset[int]:
    seen = {start}
    frontier = [start]
    while frontier:
        point = frontier.pop()
        for generator in generators:
            image = generator[point - 1]
            if image not in seen:
                seen.add(image)
                frontier.append(image)
    return frozenset(seen)


def nontrivial_cycle_lengths(permutation: Sequence[int]) -> tuple[int, ...]:
    """Return the decreasing lengths of the nontrivial cycles."""

    seen: set[int] = set()
    lengths: list[int] = []
    for start in range(1, len(permutation) + 1):
        if start in seen:
            continue
        point = start
        length = 0
        while point not in seen:
            seen.add(point)
            length += 1
            point = permutation[point - 1]
        if length > 1:
            lengths.append(length)
    return tuple(sorted(lengths, reverse=True))


@dataclass(frozen=True)
class AttachmentSlot:
    name: str
    source_location: str
    target_location: str
    tangential_index: int
    residue_different: int
    forces_new_interior_branch: bool


CASES = (
    "squarefree_one_packet",
    "double_same_target",
    "double_distinct_targets",
)

PACKET_COUNTS = {
    "squarefree_one_packet": 1,
    "double_same_target": 2,
    "double_distinct_targets": 2,
}

REQUIRED_SLOT_NAMES = (
    "endpoint_s=0",
    "endpoint_s=infinity",
    "interior_s=-1",
    "interior_denominator_root_1",
    "interior_denominator_root_2",
)


def target_fan_audit() -> dict[str, object]:
    left = (3, 1)
    terminal = (5, 2)
    right = (2, 1)
    eta = (2, -5)
    if determinant(left, terminal) != 1:
        raise AssertionError("the left target cone is not regular")
    if determinant(terminal, right) != 1:
        raise AssertionError("the right target cone is not regular")
    if tuple(x + y for x, y in zip(left, right)) != terminal:
        raise AssertionError("the target-ray self-intersection relation changed")
    eta_orders = tuple(
        ray[0] * eta[0] + ray[1] * eta[1]
        for ray in (left, terminal, right)
    )
    if eta_orders != (1, 0, -1):
        raise AssertionError("the target residue-coordinate orientation changed")
    return {
        "regular_fan_rays_left_terminal_right": [
            list(left),
            list(terminal),
            list(right),
        ],
        "adjacent_determinants": [1, 1],
        "ray_relation": "(3,1)+(2,1)=(5,2)",
        "terminal_self_intersection_in_this_local_fan": -1,
        "residue_coordinate": "eta=a^2/b^5",
        "eta_orders_on_left_terminal_right": list(eta_orders),
        "residue_map": "h=eta^-1",
        "node_orientation": {
            "h=infinity": "intersection of rays (3,1) and (5,2)",
            "h=0": "intersection of rays (5,2) and (2,1)",
        },
        "warning": (
            "the self-intersection is local-fan data before any further "
            "global blowups at the marked attachment points"
        ),
    }


def attachment_slots() -> tuple[AttachmentSlot, ...]:
    slots = (
        AttachmentSlot(
            name="endpoint_s=0",
            source_location="toric endpoint s=0",
            target_location="target node h=0",
            tangential_index=1,
            residue_different=0,
            forces_new_interior_branch=False,
        ),
        AttachmentSlot(
            name="endpoint_s=infinity",
            source_location="toric endpoint s=infinity",
            target_location="smooth target point h=125/729",
            tangential_index=3,
            residue_different=2,
            forces_new_interior_branch=False,
        ),
        AttachmentSlot(
            name="interior_s=-1",
            source_location="interior point s=-1",
            target_location="target node h=0",
            tangential_index=5,
            residue_different=4,
            forces_new_interior_branch=True,
        ),
        AttachmentSlot(
            name="interior_denominator_root_1",
            source_location="first root of 9*s^2+15*s+5",
            target_location="target node h=infinity",
            tangential_index=3,
            residue_different=2,
            forces_new_interior_branch=True,
        ),
        AttachmentSlot(
            name="interior_denominator_root_2",
            source_location="second root of 9*s^2+15*s+5",
            target_location="target node h=infinity",
            tangential_index=3,
            residue_different=2,
            forces_new_interior_branch=True,
        ),
    )
    if tuple(slot.name for slot in slots) != REQUIRED_SLOT_NAMES:
        raise AssertionError("the terminal attachment-slot order changed")
    if sum(slot.residue_different for slot in slots) != 10:
        raise AssertionError("the residue different no longer totals ten")
    if sum(slot.forces_new_interior_branch for slot in slots) != 3:
        raise AssertionError("the interior attachment count changed")
    return slots


def topology_bounds(packet_count: int) -> dict[str, int | str]:
    if packet_count == 1:
        minimum_components = 16
        minimum_leaves = 6
    elif packet_count == 2:
        minimum_components = 25
        minimum_leaves = 10
    else:
        raise ValueError("only one or two terminal packets occur in F2")
    return {
        "terminal_vertex_valency_lower_bound": 5,
        "minimum_source_boundary_components": minimum_components,
        "minimum_source_boundary_leaves": minimum_leaves,
        "proof": (
            "six blowups extract the shared (-5,1) carrier; each principal "
            "center contributes a six-blowup (5,12) arm and three interior "
            "attachment blowups; the leaf count follows from the compiled tree"
        ),
    }


def finite_normalization_case(case: str) -> dict[str, object]:
    if case not in CASES:
        raise ValueError(f"unknown F2 attachment case: {case}")
    packets = PACKET_COUNTS[case]
    if case == "squarefree_one_packet":
        target_equations = [
            {
                "target_component": "T_terminal",
                "known_boundary_rows": [[1, 6]],
                "known_contribution": 6,
                "equation": "d=6+rho_T",
                "remainder_constraint": "rho_T>=0",
            }
        ]
        degree_floor = 6
        spectator_count = 2
        packet_placement = "one terminal target component"
    elif case == "double_same_target":
        target_equations = [
            {
                "target_component": "T_terminal",
                "known_boundary_rows": [[1, 6], [1, 6]],
                "known_contribution": 12,
                "equation": "d=12+rho_T",
                "remainder_constraint": "rho_T>=0",
            }
        ]
        degree_floor = 12
        spectator_count = 0
        packet_placement = "two distinct source valuations over one target divisor"
    else:
        target_equations = [
            {
                "target_component": "T_terminal_1",
                "known_boundary_rows": [[1, 6]],
                "known_contribution": 6,
                "equation": "d=6+rho_1",
                "remainder_constraint": "rho_1>=0",
            },
            {
                "target_component": "T_terminal_2",
                "known_boundary_rows": [[1, 6]],
                "known_contribution": 6,
                "equation": "d=6+rho_2",
                "remainder_constraint": "rho_2>=0 and rho_1=rho_2",
            },
        ]
        degree_floor = 6
        spectator_count = 0
        packet_placement = "one source valuation over each of two target divisors"
    result = {
        "case": case,
        "principal_packet_count": packets,
        "packet_placement": packet_placement,
        "simple_R_spectator_orbit_count": spectator_count,
        "geometric_degree_floor": degree_floor,
        "target_component_equations": target_equations,
        "source_topology": topology_bounds(packets),
        "minimal_source_principal_boundary": minimal_source_principal_boundary(case),
        "purity_obligation": {
            "terminal_rows_have_transverse_different": 0,
            "required_elsewhere": (
                "a separate row over an affine nonproperness curve with "
                "e>1, together with its positive affine companion"
            ),
            "does_not_raise_current_degree_floor": True,
        },
        "status": "global_completion_input_required",
    }
    if case == "double_distinct_targets":
        result.update(
            geometric_degree_floor=12,
            status="excluded_by_target_valuation_uniqueness",
            exclusion_reason=(
                "both packets have global target orders (a,b)=(5,2), "
                "uniformizer b^3/a, and nonconstant residue a^2/b^5; "
                "both therefore dominate the same extracted target divisor"
            ),
        )
    return result


def ledger_contract() -> dict[str, object]:
    return {
        "source_boundary": {
            "required": [
                "a complete geometric SNC boundary component list",
                "the dual-tree edges and every self-intersection",
                "one terminal-component label for each principal packet",
                "the original orders nu(x),nu(y_old),nu(x*y_old^5-center)=(-25,5,12) and the carrier center for every terminal component",
                "a certified translated/Kummer identification and complete proximity chain for every terminal component",
                "five distinct neighboring components for each terminal packet, indexed by the compiled slots",
            ],
            "hard_gates": [
                "the dual graph is a connected tree",
                "each terminal component has valency at least five",
                "the complete boundary intersection matrix has determinant plus or minus one",
                "the intersection form has Hodge inertia (1,n-1,0)",
            ],
        },
        "class_group_and_units": {
            "forced_affine_output": "Cl(A2)=0 and O(A2)^*=k^*",
            "required_input": (
                "the complete source boundary intersection/valuation matrix; "
                "a partial local fan cannot certify either exact sequence"
            ),
            "candidate_gate": (
                "for a declared complete smooth A2 boundary, its component "
                "classes form the Picard basis and the intersection matrix "
                "must be unimodular"
            ),
        },
        "canonical": {
            "required_input": "all source self-intersections and genera",
            "equation": "M*kappa=(-2-E_i^2)_i for rational boundary components",
            "candidate_gate": "the solved canonical coefficient vector is integral",
        },
        "finite_normalization": {
            "required_input": [
                "one common geometric field degree d",
                "an exhaustive boundary/affine pullback ledger over every target curve used",
                "finite-flatness and source-to-target transfer certificates for every ledger",
                "certified placement of each principal packet and spectator orbit",
            ],
            "already_compiled": "the three symbolic target-component equations in cases[]",
        },
        "spectators": {
            "squarefree_required_count": 2,
            "missing_for_each": [
                "target component and target branch value",
                "residue degree and transverse index",
                "actual inertia permutation on global sheets",
                "node/endpoint incidence and Kummer-torsor identification",
            ],
        },
        "global_meridians": {
            "already_certified": (
                "the local A6 triple has passport (5,1)|(3,3)|(3,1,1,1) "
                "and product one"
            ),
            "required_input": (
                "global sheet degree and every additional branch cycle on "
                "each compact target component"
            ),
            "hard_gates": [
                "every listed permutation has the declared degree",
                "the ordered product is the identity",
                "the generated action is transitive for a connected cover",
                "one product-one meridian system is supplied for every target ledger",
            ],
        },
        "refusal": (
            "without these declarations the compiler reports incomplete; it "
            "does not promote cover contacts to ramification rows or infer a "
            "global completion from the local toric fan"
        ),
    }


def source_boundary_audit(
    source: object,
    case: str,
) -> dict[str, object]:
    expected_packets = PACKET_COUNTS[case]
    reasons: list[str] = []
    exclusions: list[str] = []
    if not isinstance(source, dict):
        return {
            "status": "incomplete",
            "reasons": ["source_boundary is missing"],
        }
    components_data = source.get("components")
    edges_data = source.get("edges")
    terminal_components = source.get("terminal_components")
    neighbors = source.get("attachment_neighbors")
    complete = source.get("complete") is True
    if not complete:
        reasons.append("the source boundary is not declared complete")
    if source.get("open_surface") != "A2":
        reasons.append("the source open surface is not certified as A2")
    for field in (
        "smooth_completion_certified",
        "snc_certified",
        "geometric_components_certified",
    ):
        if source.get(field) is not True:
            reasons.append(f"{field} is missing")
    if not isinstance(components_data, list) or not components_data:
        reasons.append("source components are missing")
    if not isinstance(edges_data, list):
        reasons.append("source dual-tree edges are missing")
    if not isinstance(terminal_components, list):
        reasons.append("terminal component labels are missing")
    if not isinstance(neighbors, dict):
        reasons.append("attachment-neighbor labels are missing")
    if reasons:
        return {"status": "incomplete", "reasons": reasons}

    assert isinstance(components_data, list)
    assert isinstance(edges_data, list)
    assert isinstance(terminal_components, list)
    assert isinstance(neighbors, dict)
    names: list[str] = []
    self_intersections: dict[str, int] = {}
    genera: dict[str, int] = {}
    component_records: dict[str, dict[str, object]] = {}
    for component in components_data:
        if not isinstance(component, dict):
            exclusions.append("a source component is not an object")
            continue
        name = component.get("name")
        self_intersection = component.get("self_intersection")
        genus = component.get("genus")
        if not isinstance(name, str) or not name:
            exclusions.append("a source component has no valid name")
            continue
        if not isinstance(self_intersection, int):
            exclusions.append(f"{name} has no integral self-intersection")
            continue
        if not isinstance(genus, int) or genus < 0:
            exclusions.append(f"{name} has no valid genus")
            continue
        names.append(name)
        self_intersections[name] = self_intersection
        genera[name] = genus
        component_records[name] = component
        if genus != 0:
            exclusions.append(f"{name} is not a rational A2 boundary component")
    if len(set(names)) != len(names):
        exclusions.append("source component names are not distinct")
    name_set = set(names)
    if len(terminal_components) != expected_packets:
        exclusions.append("the number of terminal components does not match the F2 case")
    if len(set(terminal_components)) != len(terminal_components):
        exclusions.append("terminal component labels are not distinct")
    if any(name not in name_set for name in terminal_components):
        exclusions.append("a terminal label is not a source component")

    edge_set: set[tuple[str, str]] = set()
    for edge in edges_data:
        if (
            not isinstance(edge, list)
            or len(edge) != 2
            or not all(isinstance(name, str) for name in edge)
        ):
            exclusions.append("a source edge is malformed")
            continue
        left, right = edge
        if left == right or left not in name_set or right not in name_set:
            exclusions.append("a source edge has an invalid endpoint")
            continue
        normalized = tuple(sorted((left, right)))
        if normalized in edge_set:
            exclusions.append("the source dual graph has a repeated edge")
        edge_set.add(normalized)

    adjacency = {name: set() for name in names}
    for left, right in edge_set:
        adjacency[left].add(right)
        adjacency[right].add(left)
    if names:
        seen = {names[0]}
        frontier = [names[0]]
        while frontier:
            point = frontier.pop()
            for other in adjacency[point]:
                if other not in seen:
                    seen.add(other)
                    frontier.append(other)
        if seen != name_set:
            exclusions.append("the source dual graph is disconnected")
        if len(edge_set) != len(names) - 1:
            exclusions.append("the source dual graph is not a tree")

    for index, terminal in enumerate(terminal_components, start=1):
        if terminal not in name_set:
            continue
        packet_key = f"packet_{index}"
        marking = neighbors.get(packet_key)
        if not isinstance(marking, dict):
            exclusions.append(f"{packet_key} has no attachment-neighbor map")
            continue
        marked_neighbors = [marking.get(slot) for slot in REQUIRED_SLOT_NAMES]
        if not all(isinstance(name, str) for name in marked_neighbors):
            exclusions.append(f"{packet_key} does not label all five slots")
            continue
        if len(set(marked_neighbors)) != 5:
            exclusions.append(f"{packet_key} does not have five distinct neighbors")
        if any(name not in adjacency[terminal] for name in marked_neighbors):
            exclusions.append(f"{packet_key} has a marked non-neighbor")
        if len(adjacency[terminal]) < 5:
            exclusions.append(f"{terminal} has valency below five")
        if genera.get(terminal) != 0:
            exclusions.append(f"{terminal} is not rational")
        terminal_record = component_records.get(terminal, {})
        original_orders = terminal_record.get("original_valuation_orders")
        expected_orders = {
            "x": -25,
            "y_old": 5,
            "x*y_old^5-center": 12,
        }
        if original_orders is None:
            reasons.append(f"{terminal} has no original valuation-order declaration")
        elif original_orders != expected_orders:
            exclusions.append(f"{terminal} does not have original orders (-25,5,12)")
        expected_center = "1" if index == 1 else "rho_double"
        if terminal_record.get("carrier_center") != expected_center:
            reasons.append(f"{terminal} has no certified carrier center {expected_center}")
        if terminal_record.get("translated_kummer_identification_certified") is not True:
            reasons.append(f"{terminal} is not identified with the translated Kummer valuation")
        if terminal_record.get("proximity_chain_certified") is not True:
            reasons.append(f"{terminal} has no certified proximity chain")

    matrix = [
        [
            self_intersections[left]
            if left == right
            else int(tuple(sorted((left, right))) in edge_set)
            for right in names
        ]
        for left in names
    ] if len(self_intersections) == len(names) else []
    matrix_determinant: int | None = None
    canonical: list[int] | None = None
    if matrix:
        matrix_determinant = bareiss_determinant(matrix)
        if abs(matrix_determinant) != 1:
            exclusions.append("the complete A2 boundary matrix is not unimodular")
        else:
            adjunction = [
                2 * genera[name] - 2 - self_intersections[name]
                for name in names
            ]
            solution = solve_linear_system(matrix, adjunction)
            if any(value.denominator != 1 for value in solution):
                exclusions.append("the canonical boundary coefficients are not integral")
            else:
                canonical = [int(value) for value in solution]

    inertia = symmetric_inertia(matrix) if matrix else None
    if matrix and inertia != (1, len(matrix) - 1, 0):
        exclusions.append("the complete A2 boundary matrix has the wrong Hodge signature")

    if exclusions:
        status = "excluded_candidate"
    elif reasons:
        status = "incomplete"
    else:
        status = "passes_source_boundary_gates"
    leaf_count = sum(len(adjacency[name]) == 1 for name in names)
    bound = topology_bounds(expected_packets)
    if names and len(names) < int(bound["minimum_source_boundary_components"]):
        exclusions.append("the source boundary has too few components")
        status = "excluded_candidate"
    if names and leaf_count < int(bound["minimum_source_boundary_leaves"]):
        exclusions.append("the source boundary has too few leaves")
        status = "excluded_candidate"
    return {
        "status": status,
        "reasons": [*exclusions, *reasons],
        "terminal_packet_components": {
            f"packet_{index}": component
            for index, component in enumerate(terminal_components, start=1)
        },
        "component_order": names,
        "component_count": len(names),
        "leaf_count": leaf_count,
        "intersection_matrix": matrix,
        "intersection_determinant": matrix_determinant,
        "intersection_inertia_positive_negative_zero": (
            list(inertia) if inertia is not None else None
        ),
        "canonical_coefficients_in_component_order": canonical,
    }


def target_ledger_audit(
    ledgers: object,
    degree: object,
    case: str,
    packet_sources: object,
) -> dict[str, object]:
    reasons: list[str] = []
    exclusions: list[str] = []
    if not isinstance(degree, int) or degree <= 0:
        reasons.append("a positive geometric_degree is missing")
    if not isinstance(ledgers, list) or not ledgers:
        reasons.append("target_ledgers are missing")
    if reasons:
        return {"status": "incomplete", "reasons": reasons}
    assert isinstance(degree, int)
    assert isinstance(ledgers, list)
    packet_targets: dict[str, str] = {}
    ledger_names: list[str] = []
    purity_witnesses: list[str] = []
    compiled_ledgers: list[dict[str, object]] = []
    for ledger in ledgers:
        if not isinstance(ledger, dict):
            exclusions.append("a target ledger is not an object")
            continue
        name = ledger.get("name")
        center = ledger.get("center")
        boundary_rows = ledger.get("boundary_rows")
        affine_rows = ledger.get("affine_rows")
        exhaustive = ledger.get("exhaustive") is True
        if not isinstance(name, str) or not name:
            exclusions.append("a target ledger has no name")
            continue
        ledger_names.append(name)
        if not isinstance(boundary_rows, list) or not isinstance(affine_rows, list):
            exclusions.append(f"{name} has malformed pullback rows")
            continue
        contribution = 0
        ramified_affine_center = False
        for row in boundary_rows:
            if not isinstance(row, dict):
                exclusions.append(f"{name} has a malformed boundary row")
                continue
            e = row.get("e")
            f = row.get("f")
            packet = row.get("principal_packet")
            if not isinstance(e, int) or e <= 0 or not isinstance(f, int) or f <= 0:
                exclusions.append(f"{name} has a nonpositive boundary row")
                continue
            contribution += e * f
            if isinstance(packet, str):
                if packet in packet_targets:
                    exclusions.append(f"principal packet {packet} is listed twice")
                packet_targets[packet] = name
                if (e, f) != (1, 6):
                    exclusions.append(f"principal packet {packet} lost row (1,6)")
                if center != "infinity":
                    exclusions.append(f"principal packet {packet} is not centered at infinity")
                if (
                    isinstance(packet_sources, dict)
                    and row.get("source_component") != packet_sources.get(packet)
                ):
                    exclusions.append(
                        f"principal packet {packet} is not attached to its declared source component"
                    )
            if center == "affine_nonproperness" and e > 1:
                ramified_affine_center = True
        affine_contribution = 0
        for row in affine_rows:
            if not isinstance(row, dict) or not isinstance(row.get("f"), int) or row["f"] <= 0:
                exclusions.append(f"{name} has a malformed affine row")
                continue
            affine_contribution += row["f"]
        contribution += affine_contribution
        if ledger.get("finite_flat_certified") is not True:
            reasons.append(f"{name} is not certified finite flat")
        if ledger.get("target_transfer_certified") is not True:
            reasons.append(f"{name} has no certified source-to-target transfer")
        if contribution > degree:
            exclusions.append(f"{name} already exceeds the geometric degree")
        if not exhaustive:
            reasons.append(f"{name} is not declared exhaustive")
        elif contribution != degree:
            exclusions.append(f"{name} does not satisfy the degree identity")
        if ramified_affine_center and affine_contribution > 0:
            purity_witnesses.append(name)
        compiled_ledgers.append(
            {
                "name": name,
                "center": center,
                "boundary_contribution": contribution - affine_contribution,
                "affine_contribution": affine_contribution,
                "total_contribution": contribution,
                "degree_identity_holds": exhaustive and contribution == degree,
            }
        )

    if len(set(ledger_names)) != len(ledger_names):
        exclusions.append("target ledger names are not distinct")

    required_packets = {f"packet_{index}" for index in range(1, PACKET_COUNTS[case] + 1)}
    if set(packet_targets) != required_packets:
        exclusions.append("the principal-packet target placement is not exhaustive")
    elif case == "double_same_target" and len(set(packet_targets.values())) != 1:
        exclusions.append("the same-target case places its packets on distinct ledgers")
    elif case == "double_distinct_targets" and len(set(packet_targets.values())) != 2:
        exclusions.append("the distinct-target case does not use two ledgers")
    if not purity_witnesses:
        reasons.append("no certified affine purity row with e>1 and an affine companion is present")

    if exclusions:
        status = "excluded_candidate"
    elif reasons:
        status = "incomplete"
    else:
        status = "passes_finite_normalization_and_purity_gates"
    return {
        "status": status,
        "reasons": [*exclusions, *reasons],
        "geometric_degree": degree,
        "packet_targets": packet_targets,
        "purity_witness_ledgers": purity_witnesses,
        "ledgers": compiled_ledgers,
    }


def spectator_audit(
    spectators: object,
    case: str,
    degree: object,
    target_ledger_names: set[str],
) -> dict[str, object]:
    required = 2 if case == "squarefree_one_packet" else 0
    if required == 0:
        if spectators not in (None, []):
            return {
                "status": "excluded_candidate",
                "reasons": ["the double-root case has no simple-R spectator orbit"],
            }
        return {"status": "not_applicable", "required_count": 0, "reasons": []}
    if not isinstance(spectators, list) or len(spectators) != required:
        return {
            "status": "incomplete",
            "required_count": required,
            "reasons": ["two squarefree simple-R spectator orbits are not classified"],
        }
    missing: list[str] = []
    names: list[str] = []
    required_fields = (
        "name",
        "target_ledger",
        "target_branch_value",
        "transverse_index",
        "residue_degree",
        "inertia_permutation",
        "geometry_certified",
    )
    for index, spectator in enumerate(spectators, start=1):
        if not isinstance(spectator, dict):
            missing.append(f"spectator {index} is malformed")
            continue
        absent = [field for field in required_fields if field not in spectator]
        if absent:
            missing.append(f"spectator {index} lacks {', '.join(absent)}")
        if spectator.get("geometry_certified") is not True:
            missing.append(f"spectator {index} geometry is not certified")
        name = spectator.get("name")
        if isinstance(name, str):
            names.append(name)
        if spectator.get("target_ledger") not in target_ledger_names:
            missing.append(f"spectator {index} refers to an unknown target ledger")
        e = spectator.get("transverse_index")
        f = spectator.get("residue_degree")
        if not isinstance(e, int) or e <= 0 or not isinstance(f, int) or f <= 0:
            missing.append(f"spectator {index} has invalid (e,f) data")
        inertia = spectator.get("inertia_permutation")
        if (
            not isinstance(degree, int)
            or not isinstance(inertia, list)
            or len(inertia) != degree
            or not permutation_is_valid(inertia)
        ):
            missing.append(f"spectator {index} has no valid global inertia permutation")
    if len(set(names)) != required:
        missing.append("spectator orbit names are not distinct")
    return {
        "status": "complete" if not missing else "incomplete",
        "required_count": required,
        "reasons": missing,
    }


def meridian_audit(
    systems: object,
    degree: object,
    expected_packets: int,
    expected_target_ledgers: set[str],
) -> dict[str, object]:
    if not isinstance(degree, int) or degree <= 0:
        return {"status": "incomplete", "reasons": ["geometric degree is missing"]}
    if not isinstance(systems, list) or not systems:
        return {
            "status": "incomplete",
            "reasons": ["global meridian systems are missing"],
        }
    exclusions: list[str] = []
    records: list[dict[str, object]] = []
    terminal_packet_occurrences: dict[str, str] = {}
    system_names: list[str] = []
    identity = tuple(range(1, degree + 1))
    for system in systems:
        if not isinstance(system, dict) or not isinstance(system.get("cycles"), list):
            exclusions.append("a meridian system is malformed")
            continue
        name = system.get("name", "unnamed")
        if not isinstance(name, str):
            exclusions.append("a meridian system has no valid target-ledger name")
            continue
        system_names.append(name)
        cycles = system["cycles"]
        parsed: list[tuple[int, ...]] = []
        for cycle in cycles:
            if not isinstance(cycle, list) or len(cycle) != degree or not permutation_is_valid(cycle):
                exclusions.append(f"{name} contains an invalid degree-{degree} permutation")
                continue
            parsed.append(tuple(cycle))
        if len(parsed) != len(cycles):
            continue
        packet_indices = system.get("terminal_packets")
        if not isinstance(packet_indices, dict):
            exclusions.append(f"{name} has no terminal-packet cycle index map")
            packet_indices = {}
        for packet, indices in packet_indices.items():
            if packet in terminal_packet_occurrences:
                exclusions.append(f"{packet} occurs in more than one meridian system")
                continue
            terminal_packet_occurrences[packet] = str(name)
            if (
                not isinstance(indices, list)
                or len(indices) != 3
                or not all(isinstance(index, int) for index in indices)
                or any(index < 0 or index >= len(parsed) for index in indices)
            ):
                exclusions.append(f"{packet} has malformed terminal-cycle indices")
                continue
            triple = tuple(parsed[index] for index in indices)
            types = tuple(nontrivial_cycle_lengths(cycle) for cycle in triple)
            if types != ((5,), (3, 3), (3,)):
                exclusions.append(f"{packet} does not have the terminal A6 passport")
            local_product = identity
            for cycle in triple:
                local_product = multiply_permutations(local_product, cycle)
            support = {
                point
                for cycle in triple
                for point in range(1, degree + 1)
                if cycle[point - 1] != point
            }
            local_transitive = bool(support) and len(
                permutation_orbit(triple, min(support))
            ) == len(support)
            if local_product != identity:
                exclusions.append(f"{packet} terminal triple does not have product one")
            if len(support) != 6 or not local_transitive:
                exclusions.append(f"{packet} is not a transitive six-sheet terminal packet")
        product = identity
        for cycle in parsed:
            product = multiply_permutations(product, cycle)
        product_one = product == identity
        transitive = bool(parsed) and len(permutation_orbit(parsed, 1)) == degree
        if not product_one:
            exclusions.append(f"{name} does not have product one")
        if system.get("connected") is True and not transitive:
            exclusions.append(f"{name} is declared connected but is not transitive")
        records.append(
            {
                "name": name,
                "cycle_count": len(parsed),
                "product_one": product_one,
                "transitive": transitive,
                "terminal_packets": sorted(packet_indices),
            }
        )
    required_packets = {
        f"packet_{index}" for index in range(1, expected_packets + 1)
    }
    if set(terminal_packet_occurrences) != required_packets:
        exclusions.append("the global meridian systems do not place every principal packet")
    if len(set(system_names)) != len(system_names):
        exclusions.append("global meridian system names are not distinct")
    if set(system_names) != expected_target_ledgers:
        exclusions.append("global meridian systems do not cover every target ledger")
    return {
        "status": "excluded_candidate" if exclusions else "passes_declared_meridian_gates",
        "reasons": exclusions,
        "systems": records,
    }


def audit_candidate(data: object) -> dict[str, object]:
    if not isinstance(data, dict):
        raise ValueError("candidate JSON must contain an object")
    case = data.get("case")
    if case not in CASES:
        raise ValueError(f"candidate case must be one of {', '.join(CASES)}")
    assert isinstance(case, str)
    source = source_boundary_audit(data.get("source_boundary"), case)
    target = target_ledger_audit(
        data.get("target_ledgers"),
        data.get("geometric_degree"),
        case,
        source.get("terminal_packet_components"),
    )
    target_ledger_names = {
        ledger.get("name")
        for ledger in data.get("target_ledgers", [])
        if isinstance(ledger, dict) and isinstance(ledger.get("name"), str)
    } if isinstance(data.get("target_ledgers"), list) else set()
    spectators = spectator_audit(
        data.get("spectator_orbits"),
        case,
        data.get("geometric_degree"),
        target_ledger_names,
    )
    meridians = meridian_audit(
        data.get("global_meridian_systems"),
        data.get("geometric_degree"),
        PACKET_COUNTS[case],
        target_ledger_names,
    )
    case_gate = {
        "status": (
            "excluded_candidate"
            if case == "double_distinct_targets"
            else "applicable"
        ),
        "reason": (
            "the two principal packets necessarily restrict to the same "
            "target divisor (5,2)"
            if case == "double_distinct_targets"
            else "the attachment case survives target-valuation uniqueness"
        ),
    }
    audits = (source, target, spectators, meridians, case_gate)
    if any(audit["status"] == "excluded_candidate" for audit in audits):
        status = "excluded_candidate"
    elif any(audit["status"] == "incomplete" for audit in audits):
        status = "incomplete"
    else:
        status = "passes_compiled_necessary_gates_not_an_existence_proof"
    return {
        "schema": "plane-jc.f2-75-125-global-attachment-candidate-audit.v1",
        "case": case,
        "status": status,
        "case_gate": case_gate,
        "source_boundary": source,
        "target_ledgers": target,
        "spectators": spectators,
        "global_meridians": meridians,
        "claim_boundary": (
            "passing proves only compatibility with these necessary compiled "
            "gates; it does not construct a polynomial Keller map"
        ),
    }


def build_payload() -> dict[str, object]:
    slots = attachment_slots()
    cases = [finite_normalization_case(case) for case in CASES]
    payload = {
        "schema": "plane-jc.f2-75-125-global-attachment.v1",
        "status": "exact-obligation-compiler-global-completion-still-missing",
        "certified_input": {
            "source_ray": [12, -17],
            "source_ray_coordinates": ["t=X*y_translated", "z=y_translated^-1"],
            "original_source_valuation_orders_x_y_v_minus_center": [-25, 5, 12],
            "source_pole_orders_P_Q": [3, 5],
            "target_ray": [5, 2],
            "transverse_index": 1,
            "residue_degree": 6,
            "residue_map": "125*s*(s+1)^5/(9*s^2+15*s+5)^3",
            "passport": [[5, 1], [3, 3], [3, 1, 1, 1]],
            "geometric_monodromy": "A6",
            "target_fixed_deck_group": "trivial",
        },
        "source_original_proximity_resolution": source_proximity_fan_audit(),
        "target_local_fan": target_fan_audit(),
        "target_minimal_completion": target_minimal_completion_audit(),
        "target_valuation_uniqueness": target_valuation_uniqueness_audit(),
        "attachment_slots": [asdict(slot) for slot in slots],
        "attachment_summary": {
            "toric_endpoint_neighbor_count_per_packet": 2,
            "new_interior_branch_count_per_packet": 3,
            "terminal_valency_lower_bound": 5,
            "node_h=0": ["endpoint_s=0", "interior_s=-1"],
            "node_h=infinity": [
                "interior_denominator_root_1",
                "interior_denominator_root_2",
            ],
            "smooth_h=125/729": ["endpoint_s=infinity"],
        },
        "cases": cases,
        "ledger_contract": ledger_contract(),
        "global_result": {
            "new_exact_gates": [
                "the translated terminal chart restricts to the original nonmonomial valuation nu(x),nu(y_old),nu(x*y_old^5-center)=(-25,5,12)",
                "six blowups extract the shared carrier ray (-5,1), and six further blowups at each marked carrier center extract a principal terminal arm",
                "the target ray (5,2) is extracted by four blowups and gives the unimodular boundary chain (-2,-2,-1,-3,-2)",
                "one terminal packet forces a 16-component, 6-leaf minimal principal boundary; two packets sharing the carrier force 25 components and 10 leaves",
                "squarefree spectator centers are the two simple R roots on the carrier; the double R root is the second principal center",
                "both double-root packets restrict to the unique target divisor (5,2), excluding the distinct-target case and forcing d>=12",
                "a declared complete A2 source boundary must have a unimodular intersection matrix and integral adjunction solution",
            ],
            "not_resolved": [
                "additional source blowups required by spectators, purity, or other map-resolution rows beyond the compiled minimal principal graph",
                "the two squarefree simple-R spectator orbit profiles",
                "the purity-forced affine ramification row",
                "the global sheet degree and global meridian systems",
            ],
            "verdict": "(75,125) remains unexcluded",
        },
        "candidate_mode": {
            "command": (
                ".venv/bin/python plane-jc/cas/"
                "verify_f2_75_125_global_attachment.py --candidate FILE.json"
            ),
            "behavior": (
                "audit a supplied global completion against the source-tree, "
                "class/unit/canonical, finite-normalization, purity, spectator, "
                "and meridian gates"
            ),
        },
        "reproduction_command": (
            ".venv/bin/python plane-jc/cas/"
            "verify_f2_75_125_global_attachment.py"
        ),
        "software": {"python": "standard library"},
    }
    return payload


def artifact_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=ARTIFACT)
    parser.add_argument("--refresh", action="store_true")
    parser.add_argument("--candidate", type=Path)
    args = parser.parse_args()

    if args.candidate is not None:
        candidate = json.loads(args.candidate.read_text())
        print(json.dumps(audit_candidate(candidate), indent=2, sort_keys=True))
        return

    payload = build_payload()
    artifact = args.artifact.resolve()
    if args.refresh:
        artifact.parent.mkdir(parents=True, exist_ok=True)
        artifact.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        try:
            display = artifact.relative_to(ROOT)
        except ValueError:
            display = artifact
        print(f"WROTE {display}")
    else:
        expected = json.loads(artifact.read_text())
        if expected != payload:
            raise AssertionError(
                "the pinned F2 global-attachment artifact is stale; inspect "
                "the change before using --refresh"
            )

    print("F2_75_125_TARGET_FAN=(3,1),(5,2),(2,1)")
    print("F2_75_125_ATTACHMENTS_PER_TERMINAL=5")
    print("F2_75_125_INTERIOR_ATTACHMENTS_PER_TERMINAL=3")
    print("F2_75_125_ONE_PACKET_TREE_MIN=components:16,leaves:6")
    print("F2_75_125_TWO_PACKET_TREE_MIN=components:25,leaves:10")
    print("F2_75_125_GLOBAL_ATTACHMENT_STATUS=INCOMPLETE")
    print(f"F2_75_125_GLOBAL_ATTACHMENT_ARTIFACT_SHA256={artifact_sha256(artifact)}")


if __name__ == "__main__":
    main()

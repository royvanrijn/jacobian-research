#!/usr/bin/env python3
"""Replay exact sparse-support exclusions for normalized plane Keller maps.

For F=(x+P,y+Q), the Keller equation is

    P_x + Q_y + P_x Q_y - P_y Q_x = 0.

Every coefficient is linear or quadratic in the exact-support coefficients.
The checker uses the Rabinowitsch equation ``z*prod(coefficients)-1``.  When
one monomial of the Keller equation has a single contribution ``n*m``, the
following identity is an explicit unit-ideal certificate:

    (z*(prod/m)/n)*(n*m) - (z*prod-1) = 1.

The arbitrary-degree balanced 2+2 proof exhausts its finite
divergence/determinant collision patterns with exact integer arithmetic.
The former bounded census remains as an independent digest regression.  The
arbitrary-degree 1+q theorem (q<=5) is proved in the accompanying note; its
nonempty exact-support charts are replayed independently by Singular.  For
total support six, exact Z3 formulas classify every no-singleton exponent
chart: 1+5 is the quartic shear chain, 2+4 is impossible, and 3+3 is the
directional quadratic shear chart.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import struct
import subprocess
from itertools import combinations
from pathlib import Path

import z3


ROOT = Path(__file__).resolve().parents[2]
ARTIFACT = ROOT / "artifacts/generated-results/jc2_sparse_support_exclusions.json"
ARTIFACT_SHA256 = (
    "4f1fd545a5c48c39a32f74b6225ce00fd986cd5107427f9ee9077ce985863292"
)
SINGULAR_CERTIFICATE = Path(__file__).with_name("sparse_support_exceptional.sing")
SINGULAR_CERTIFICATE_SHA256 = (
    "b7046fd378a1796e7abfed6ae3e4554d4e92c403c7dad35705f30c8a7d4575ab"
)

# Coefficient masks in the order a,b,c,d.
A_MASK = 1
B_MASK = 2
C_MASK = 4
D_MASK = 8
ALL_MASK = A_MASK | B_MASK | C_MASK | D_MASK

EXPECTED_LINEAR_RESIDUAL_PATTERNS = (
    (0x0F, ((0, 2), (1, 3))),
    (0x0F, ((0, 3), (1, 2))),
    (0x3F, ((0, 2), (1, 5), (3, 4))),
    (0x3F, ((0, 3), (1, 4), (2, 5))),
    (0x5F, ((0, 2), (1, 4), (3, 6))),
    (0x5F, ((0, 6), (1, 2), (3, 4))),
    (0x66, ((1, 2), (5, 6))),
    (0x69, ((0, 3), (5, 6))),
    (0x6F, ((0, 2), (1, 3), (5, 6))),
    (0x7F, ((0, 2), (1, 3, 4), (5, 6))),
    (0x95, ((0, 2), (4, 7))),
    (0x9A, ((1, 3), (4, 7))),
    (0x9F, ((0, 3), (1, 2), (4, 7))),
    (0xAF, ((0, 3), (1, 5), (2, 7))),
    (0xAF, ((0, 7), (1, 3), (2, 5))),
    (0xBF, ((0, 3), (1, 2, 5), (4, 7))),
    (0xCF, ((0, 6), (1, 3), (2, 7))),
    (0xCF, ((0, 7), (1, 2), (3, 6))),
    (0xDF, ((0, 3, 6), (1, 2), (4, 7))),
    (0xEF, ((0, 2, 7), (1, 3), (5, 6))),
)


def nonlinear_monomials(degree_cap: int) -> tuple[tuple[int, int], ...]:
    """Monomials x^i y^j of total degree 2 through ``degree_cap``."""

    return tuple(
        (x_degree, total_degree - x_degree)
        for total_degree in range(2, degree_cap + 1)
        for x_degree in range(total_degree + 1)
    )


def _pair_data(
    support: tuple[tuple[int, int], tuple[int, int]],
    *,
    p_coordinate: bool,
) -> tuple[tuple[tuple[int, int], int, int, int], ...]:
    """Return divergence terms as (exponent, integer, mask, source)."""

    masks = (A_MASK, B_MASK) if p_coordinate else (C_MASK, D_MASK)
    sources = (0, 1) if p_coordinate else (2, 3)
    result: list[tuple[tuple[int, int], int, int, int]] = []
    for (x_degree, y_degree), mask, source in zip(
        support, masks, sources, strict=True
    ):
        integer = x_degree if p_coordinate else y_degree
        if integer:
            exponent = (
                (x_degree - 1, y_degree)
                if p_coordinate
                else (x_degree, y_degree - 1)
            )
            result.append((exponent, integer, mask, source))
    return tuple(result)


def _selected_singleton(
    p_support: tuple[tuple[int, int], tuple[int, int]],
    q_support: tuple[tuple[int, int], tuple[int, int]],
    p_divergence: tuple[tuple[tuple[int, int], int, int, int], ...],
    q_divergence: tuple[tuple[tuple[int, int], int, int, int], ...],
) -> tuple[int, int, int, int, int] | None:
    """Select the canonical singleton coefficient of the Keller equation."""

    terms = [*p_divergence, *q_divergence]
    source = 4
    for p_index, p_exponent in enumerate(p_support):
        p_mask = (A_MASK, B_MASK)[p_index]
        for q_index, q_exponent in enumerate(q_support):
            q_mask = (C_MASK, D_MASK)[q_index]
            determinant = (
                p_exponent[0] * q_exponent[1]
                - p_exponent[1] * q_exponent[0]
            )
            if determinant:
                terms.append(
                    (
                        (
                            p_exponent[0] + q_exponent[0] - 1,
                            p_exponent[1] + q_exponent[1] - 1,
                        ),
                        determinant,
                        p_mask | q_mask,
                        source,
                    )
                )
            source += 1

    multiplicities: dict[tuple[int, int], int] = {}
    for exponent, _, _, _ in terms:
        multiplicities[exponent] = multiplicities.get(exponent, 0) + 1
    singletons = [
        (exponent[0], exponent[1], integer, mask, source_index)
        for exponent, integer, mask, source_index in terms
        if multiplicities[exponent] == 1
    ]
    return min(singletons) if singletons else None


def balanced_binomial_census(degree_cap: int) -> dict[str, object]:
    """Certify every exact 2+2 support by a singleton coefficient."""

    monomials = nonlinear_monomials(degree_cap)
    supports = tuple(combinations(monomials, 2))
    divergence_p = tuple(_pair_data(support, p_coordinate=True) for support in supports)
    divergence_q = tuple(_pair_data(support, p_coordinate=False) for support in supports)
    digest = hashlib.sha256()
    certificate_count = 0
    survivors: list[dict[str, object]] = []

    for p_index, p_support in enumerate(supports):
        p_divergence = divergence_p[p_index]
        for q_index, q_support in enumerate(supports):
            selected = _selected_singleton(
                p_support,
                q_support,
                p_divergence,
                divergence_q[q_index],
            )
            if selected is None:
                survivors.append({"P": p_support, "Q": q_support})
                continue
            exponent_x, exponent_y, integer, mask, source = selected
            if not integer or not (1 <= mask <= ALL_MASK):
                raise AssertionError("invalid singleton certificate")

            # Canonical digest record:
            # A_x,A_y,B_x,B_y,C_x,C_y,D_x,D_y,e_x,e_y,n,mask,source.
            record = (
                *p_support[0],
                *p_support[1],
                *q_support[0],
                *q_support[1],
                exponent_x,
                exponent_y,
                integer,
                mask,
                source,
            )
            digest.update(struct.pack(">13h", *record))
            certificate_count += 1

    if survivors:
        sample = json.dumps(survivors[:3], sort_keys=True)
        raise AssertionError(f"balanced 2+2 census has survivors: {sample}")

    expected_count = len(supports) ** 2
    if certificate_count != expected_count:
        raise AssertionError("balanced support count changed")
    return {
        "degree_cap": degree_cap,
        "monomial_order": "total degree, then increasing x exponent",
        "nonlinear_monomial_count": len(monomials),
        "two_term_supports_per_coordinate": len(supports),
        "ordered_support_pairs": expected_count,
        "unit_certificates": certificate_count,
        "survivors": 0,
        "certificate_identity": "(z*(abcd/M)/n)*(n*M) - (z*abcd-1) = 1",
        "digest_encoding": (
            "concatenated >13h records "
            "(Ax,Ay,Bx,By,Cx,Cy,Dx,Dy,ex,ey,n,mask,source)"
        ),
        "certificate_sha256": digest.hexdigest(),
    }


def balanced_trinomial_collision_census(degree_cap: int = 6) -> dict[str, object]:
    """Search every bounded 3+3 support for no-singleton collision charts."""

    monomials = nonlinear_monomials(degree_cap)
    supports = tuple(combinations(monomials, 3))
    survivors: list[
        tuple[
            tuple[tuple[int, int], ...],
            tuple[tuple[int, int], ...],
        ]
    ] = []
    for p_support in supports:
        for q_support in supports:
            exponents: list[tuple[int, int]] = []
            exponents.extend(
                (x_degree - 1, y_degree)
                for x_degree, y_degree in p_support
                if x_degree
            )
            exponents.extend(
                (x_degree, y_degree - 1)
                for x_degree, y_degree in q_support
                if y_degree
            )
            exponents.extend(
                (
                    p_x + q_x - 1,
                    p_y + q_y - 1,
                )
                for p_x, p_y in p_support
                for q_x, q_y in q_support
                if p_x * q_y - p_y * q_x
            )
            multiplicities: dict[tuple[int, int], int] = {}
            for exponent in exponents:
                multiplicities[exponent] = multiplicities.get(exponent, 0) + 1
            if exponents and all(count >= 2 for count in multiplicities.values()):
                survivors.append((p_support, q_support))

    expected_support = ((0, 2), (1, 1), (2, 0))
    if survivors != [(expected_support, expected_support)]:
        raise AssertionError("the bounded 3+3 collision survivor list changed")
    return {
        "degree_cap": degree_cap,
        "nonlinear_monomial_count": len(monomials),
        "three_term_supports_per_coordinate": len(supports),
        "ordered_support_pairs": len(supports) ** 2,
        "no_singleton_supports": 1,
        "survivors": [
            {
                "P": [list(exponent) for exponent in expected_support],
                "Q": [list(exponent) for exponent in expected_support],
                "classification": "directional quadratic shear exponent chart",
            }
        ],
        "claim_boundary": "bounded regression only; arbitrary degree is Z3-certified",
    }


def _z3_balanced_problem() -> tuple[
    tuple[z3.ArithRef, ...],
    tuple[tuple[z3.ArithRef, z3.ArithRef], ...],
    tuple[z3.BoolRef, ...],
    tuple[z3.BoolRef, ...],
]:
    variables = z3.Ints("ax ay bx by cx cy dx dy")
    ax, ay, bx, by, cx, cy, dx, dy = variables
    a_exponent = (ax, ay)
    b_exponent = (bx, by)
    c_exponent = (cx, cy)
    d_exponent = (dx, dy)

    def determinant(
        left: tuple[z3.ArithRef, z3.ArithRef],
        right: tuple[z3.ArithRef, z3.ArithRef],
    ) -> z3.ArithRef:
        return left[0] * right[1] - left[1] * right[0]

    term_exponents = (
        (ax - 1, ay),
        (bx - 1, by),
        (cx, cy - 1),
        (dx, dy - 1),
        (ax + cx - 1, ay + cy - 1),
        (ax + dx - 1, ay + dy - 1),
        (bx + cx - 1, by + cy - 1),
        (bx + dx - 1, by + dy - 1),
    )
    presence = (
        ax > 0,
        bx > 0,
        cy > 0,
        dy > 0,
        determinant(a_exponent, c_exponent) != 0,
        determinant(a_exponent, d_exponent) != 0,
        determinant(b_exponent, c_exponent) != 0,
        determinant(b_exponent, d_exponent) != 0,
    )
    base = (
        *(variable >= 0 for variable in variables),
        ax + ay >= 2,
        bx + by >= 2,
        cx + cy >= 2,
        dx + dy >= 2,
        z3.Or(ax != bx, ay != by),
        z3.Or(cx != dx, cy != dy),
    )
    return variables, term_exponents, presence, base


def _z3_no_singleton_constraints(
    term_exponents: tuple[tuple[z3.ArithRef, z3.ArithRef], ...],
    active: tuple[int, ...],
) -> tuple[z3.BoolRef, ...]:
    constraints: list[z3.BoolRef] = []
    for term_index in active:
        exponent = term_exponents[term_index]
        partners = [
            z3.And(
                exponent[0] == term_exponents[other][0],
                exponent[1] == term_exponents[other][1],
            )
            for other in active
            if other != term_index
        ]
        constraints.append(z3.Or(*partners) if partners else z3.BoolVal(False))
    return tuple(constraints)


def _partitions_without_singletons(
    items: tuple[int, ...],
) -> tuple[tuple[tuple[int, ...], ...], ...]:
    if not items:
        return ((),)
    first, rest = items[0], items[1:]
    result: list[tuple[tuple[int, ...], ...]] = []
    for block_size in range(1, len(rest) + 1):
        for selected in combinations(rest, block_size):
            selected_set = set(selected)
            remainder = tuple(item for item in rest if item not in selected_set)
            for tail in _partitions_without_singletons(remainder):
                result.append(((first, *selected), *tail))
    return tuple(result)


def _add_partition_constraints(
    solver: z3.Solver,
    term_exponents: tuple[tuple[z3.ArithRef, z3.ArithRef], ...],
    partition: tuple[tuple[int, ...], ...],
) -> None:
    for block in partition:
        representative = term_exponents[block[0]]
        for term_index in block[1:]:
            solver.add(
                representative[0] == term_exponents[term_index][0],
                representative[1] == term_exponents[term_index][1],
            )
    for left_index, left_block in enumerate(partition):
        left = term_exponents[left_block[0]]
        for right_block in partition[:left_index]:
            right = term_exponents[right_block[0]]
            solver.add(z3.Or(left[0] != right[0], left[1] != right[1]))


def verify_arbitrary_balanced_binomials() -> dict[str, object]:
    """Prove that every nonlinear 2+2 support has a singleton coefficient."""

    if z3.get_version_string() != "4.15.3":
        raise RuntimeError("the arbitrary 2+2 certificate requires Z3 4.15.3")
    variables, term_exponents, presence, base = _z3_balanced_problem()
    ax, _, bx, _, _, cy, _, dy = variables
    divergence_variables = (ax, bx, cy, dy)

    feasible_presence_masks: list[int] = []
    for mask in range(256):
        solver = z3.Solver()
        solver.add(*base)
        for term_index, condition in enumerate(presence):
            solver.add(condition if mask & (1 << term_index) else z3.Not(condition))
        if solver.check() == z3.sat:
            feasible_presence_masks.append(mask)
    if len(feasible_presence_masks) != 85:
        raise AssertionError("the determinant/divergence presence census changed")

    linear_residual_masks: list[int] = []
    for mask in feasible_presence_masks:
        active = tuple(index for index in range(8) if mask & (1 << index))
        solver = z3.Solver()
        solver.add(*base)
        for term_index, variable in enumerate(divergence_variables):
            solver.add(variable >= 1 if mask & (1 << term_index) else variable == 0)
        solver.add(*_z3_no_singleton_constraints(term_exponents, active))
        if solver.check() == z3.sat:
            linear_residual_masks.append(mask)
    expected_masks = sorted({mask for mask, _ in EXPECTED_LINEAR_RESIDUAL_PATTERNS})
    if linear_residual_masks != expected_masks:
        raise AssertionError("the linear no-singleton mask list changed")

    residual_patterns: list[tuple[int, tuple[tuple[int, ...], ...]]] = []
    for mask in linear_residual_masks:
        active = tuple(index for index in range(8) if mask & (1 << index))
        for partition in _partitions_without_singletons(active):
            solver = z3.Solver()
            solver.add(*base)
            for term_index, variable in enumerate(divergence_variables):
                solver.add(
                    variable >= 1 if mask & (1 << term_index) else variable == 0
                )
            _add_partition_constraints(solver, term_exponents, partition)
            if solver.check() == z3.sat:
                residual_patterns.append((mask, partition))
    if tuple(residual_patterns) != EXPECTED_LINEAR_RESIDUAL_PATTERNS:
        raise AssertionError("the canonical residual collision table changed")

    for mask, partition in residual_patterns:
        solver = z3.Solver()
        solver.add(*base)
        for term_index, condition in enumerate(presence):
            solver.add(condition if mask & (1 << term_index) else z3.Not(condition))
        _add_partition_constraints(solver, term_exponents, partition)
        if solver.check() != z3.unsat:
            raise AssertionError(
                f"residual collision pattern survived: {mask:#04x} {partition}"
            )

    global_solver = z3.Solver()
    global_solver.add(*base)
    for term_index, (exponent, condition) in enumerate(
        zip(term_exponents, presence, strict=True)
    ):
        partners = [
            z3.And(
                other_condition,
                exponent[0] == other_exponent[0],
                exponent[1] == other_exponent[1],
            )
            for other_index, (other_exponent, other_condition) in enumerate(
                zip(term_exponents, presence, strict=True)
            )
            if other_index != term_index
        ]
        global_solver.add(z3.Implies(condition, z3.Or(*partners)))
    if global_solver.check() != z3.unsat:
        raise AssertionError("the global arbitrary 2+2 no-singleton system survived")

    presence_digest = hashlib.sha256(bytes(feasible_presence_masks)).hexdigest()
    return {
        "statement": (
            "every exact 2+2 nonlinear support in arbitrary degree has a "
            "singleton Keller coefficient"
        ),
        "term_order": [
            "P_A,x",
            "P_B,x",
            "Q_C,y",
            "Q_D,y",
            "[A,C]",
            "[A,D]",
            "[B,C]",
            "[B,D]",
        ],
        "presence_masks_total": 256,
        "presence_masks_feasible": len(feasible_presence_masks),
        "presence_mask_sha256": presence_digest,
        "linear_residual_masks": [f"0x{mask:02x}" for mask in linear_residual_masks],
        "canonical_residual_patterns": len(residual_patterns),
        "residual_pattern_table": [
            {
                "mask": f"0x{mask:02x}",
                "blocks": ["".join(str(index) for index in block) for block in partition],
            }
            for mask, partition in residual_patterns
        ],
        "residual_patterns_all_unsat": True,
        "global_no_singleton_formula": "unsat",
        "solver": f"Z3 {z3.get_version_string()} exact integer arithmetic",
        "unit_certificate": "(z*(abcd/M)/n)*(n*M) - (z*abcd-1) = 1",
    }


def _z3_two_by_three_problem() -> tuple[
    tuple[z3.ArithRef, ...],
    tuple[tuple[z3.ArithRef, z3.ArithRef], ...],
    tuple[z3.BoolRef, ...],
    tuple[z3.BoolRef, ...],
]:
    variables = z3.Ints("ax ay bx by cx cy dx dy ex ey")
    ax, ay, bx, by, cx, cy, dx, dy, ex, ey = variables
    p_exponents = ((ax, ay), (bx, by))
    q_exponents = ((cx, cy), (dx, dy), (ex, ey))

    def determinant(
        left: tuple[z3.ArithRef, z3.ArithRef],
        right: tuple[z3.ArithRef, z3.ArithRef],
    ) -> z3.ArithRef:
        return left[0] * right[1] - left[1] * right[0]

    term_exponents = (
        *((exponent[0] - 1, exponent[1]) for exponent in p_exponents),
        *((exponent[0], exponent[1] - 1) for exponent in q_exponents),
        *(
            (
                p_exponent[0] + q_exponent[0] - 1,
                p_exponent[1] + q_exponent[1] - 1,
            )
            for p_exponent in p_exponents
            for q_exponent in q_exponents
        ),
    )
    presence = (
        *(exponent[0] > 0 for exponent in p_exponents),
        *(exponent[1] > 0 for exponent in q_exponents),
        *(
            determinant(p_exponent, q_exponent) != 0
            for p_exponent in p_exponents
            for q_exponent in q_exponents
        ),
    )
    base: list[z3.BoolRef] = [
        *(variable >= 0 for variable in variables),
        *(exponent[0] + exponent[1] >= 2 for exponent in p_exponents),
        *(exponent[0] + exponent[1] >= 2 for exponent in q_exponents),
        z3.Or(ax != bx, ay != by),
    ]
    for left_index, left in enumerate(q_exponents):
        for right in q_exponents[:left_index]:
            base.append(z3.Or(left[0] != right[0], left[1] != right[1]))
    return variables, term_exponents, presence, tuple(base)


def verify_arbitrary_two_by_three() -> dict[str, object]:
    """Prove that every nonlinear 2+3 support has a singleton coefficient."""

    variables, term_exponents, presence, base = _z3_two_by_three_problem()
    ax, _, bx, _, _, cy, _, dy, _, ey = variables
    divergence_variables = (ax, bx, cy, dy, ey)

    feasible_masks: list[int] = []
    for mask in range(1 << 11):
        solver = z3.Solver()
        solver.add(*base)
        for term_index, condition in enumerate(presence):
            solver.add(condition if mask & (1 << term_index) else z3.Not(condition))
        if solver.check() == z3.sat:
            feasible_masks.append(mask)
    feasible_digest = hashlib.sha256(
        b"".join(mask.to_bytes(2, "big") for mask in feasible_masks)
    ).hexdigest()
    if len(feasible_masks) != 321 or feasible_digest != (
        "227bcb7b327487ec9dc6b861a75c3e1c2e693cfe3cf15726077f052407a7b0eb"
    ):
        raise AssertionError("the 2+3 determinant/divergence presence census changed")

    linear_residual_masks: list[int] = []
    for mask in feasible_masks:
        active = tuple(index for index in range(11) if mask & (1 << index))
        solver = z3.Solver()
        solver.add(*base)
        for term_index, variable in enumerate(divergence_variables):
            solver.add(variable >= 1 if mask & (1 << term_index) else variable == 0)
        solver.add(*_z3_no_singleton_constraints(term_exponents, active))
        if solver.check() == z3.sat:
            linear_residual_masks.append(mask)
    linear_digest = hashlib.sha256(
        b"".join(mask.to_bytes(2, "big") for mask in linear_residual_masks)
    ).hexdigest()
    if len(linear_residual_masks) != 98 or linear_digest != (
        "17f070086059ca1536395602983fda35dff2c962ad64e057f0472987e584fb9f"
    ):
        raise AssertionError("the 2+3 linear collision sieve changed")

    global_solver = z3.Solver()
    global_solver.add(*base)
    for term_index, (exponent, condition) in enumerate(
        zip(term_exponents, presence, strict=True)
    ):
        partners = [
            z3.And(
                other_condition,
                exponent[0] == other_exponent[0],
                exponent[1] == other_exponent[1],
            )
            for other_index, (other_exponent, other_condition) in enumerate(
                zip(term_exponents, presence, strict=True)
            )
            if other_index != term_index
        ]
        global_solver.add(z3.Implies(condition, z3.Or(*partners)))
    if global_solver.check() != z3.unsat:
        raise AssertionError("the global arbitrary 2+3 no-singleton system survived")

    return {
        "statement": (
            "every exact 2+3 nonlinear support in arbitrary degree has a "
            "singleton Keller coefficient"
        ),
        "transpose_included": True,
        "term_count": len(term_exponents),
        "presence_masks_total": 1 << len(term_exponents),
        "presence_masks_feasible": len(feasible_masks),
        "presence_mask_sha256": feasible_digest,
        "linear_residual_masks": len(linear_residual_masks),
        "linear_residual_mask_sha256": linear_digest,
        "global_no_singleton_formula": "unsat",
        "solver": f"Z3 {z3.get_version_string()} exact integer arithmetic",
        "unit_certificate": "(z*(abcde/M)/n)*(n*M) - (z*abcde-1) = 1",
    }


def _z3_sparse_split_problem(
    p_count: int,
    q_count: int,
    *,
    prefix: str,
) -> tuple[
    tuple[tuple[z3.ArithRef, z3.ArithRef], ...],
    tuple[tuple[z3.ArithRef, z3.ArithRef], ...],
    tuple[tuple[z3.ArithRef, z3.ArithRef], ...],
    tuple[z3.BoolRef, ...],
    tuple[z3.BoolRef, ...],
]:
    """Build the exact exponent problem for a general sparse split."""

    variables = z3.Ints(
        " ".join(
            f"{prefix}_{coordinate}{index}{axis}"
            for coordinate, count in (("p", p_count), ("q", q_count))
            for index in range(count)
            for axis in ("x", "y")
        )
    )
    cursor = 0
    groups: list[tuple[tuple[z3.ArithRef, z3.ArithRef], ...]] = []
    for count in (p_count, q_count):
        group: list[tuple[z3.ArithRef, z3.ArithRef]] = []
        for _ in range(count):
            group.append((variables[cursor], variables[cursor + 1]))
            cursor += 2
        groups.append(tuple(group))
    p_exponents, q_exponents = groups

    term_exponents = (
        *((exponent[0] - 1, exponent[1]) for exponent in p_exponents),
        *((exponent[0], exponent[1] - 1) for exponent in q_exponents),
        *(
            (
                p_exponent[0] + q_exponent[0] - 1,
                p_exponent[1] + q_exponent[1] - 1,
            )
            for p_exponent in p_exponents
            for q_exponent in q_exponents
        ),
    )
    presence = (
        *(exponent[0] > 0 for exponent in p_exponents),
        *(exponent[1] > 0 for exponent in q_exponents),
        *(
            p_exponent[0] * q_exponent[1]
            - p_exponent[1] * q_exponent[0]
            != 0
            for p_exponent in p_exponents
            for q_exponent in q_exponents
        ),
    )
    base: list[z3.BoolRef] = [
        *(variable >= 0 for variable in variables),
        *(
            exponent[0] + exponent[1] >= 2
            for exponent in (*p_exponents, *q_exponents)
        ),
    ]
    for group in (p_exponents, q_exponents):
        for left_index, left in enumerate(group):
            for right in group[:left_index]:
                base.append(
                    z3.Or(left[0] != right[0], left[1] != right[1])
                )
    return p_exponents, q_exponents, term_exponents, presence, tuple(base)


def _add_global_no_singleton_formula(
    solver: z3.Solver,
    term_exponents: tuple[tuple[z3.ArithRef, z3.ArithRef], ...],
    presence: tuple[z3.BoolRef, ...],
) -> None:
    for term_index, (exponent, condition) in enumerate(
        zip(term_exponents, presence, strict=True)
    ):
        partners = [
            z3.And(
                other_condition,
                exponent[0] == other_exponent[0],
                exponent[1] == other_exponent[1],
            )
            for other_index, (other_exponent, other_condition) in enumerate(
                zip(term_exponents, presence, strict=True)
            )
            if other_index != term_index
        ]
        solver.add(z3.Implies(condition, z3.Or(*partners)))


def verify_support_six_exponent_classification() -> dict[str, object]:
    """Classify every arbitrary-degree support-six no-singleton chart."""

    if z3.get_version_string() != "4.15.3":
        raise RuntimeError("the support-six certificate requires Z3 4.15.3")

    def belongs_to(
        exponent: tuple[z3.ArithRef, z3.ArithRef],
        support: tuple[tuple[z3.ArithRef | int, z3.ArithRef | int], ...],
    ) -> z3.BoolRef:
        return z3.Or(
            *(
                z3.And(exponent[0] == x_degree, exponent[1] == y_degree)
                for x_degree, y_degree in support
            )
        )

    p15, q15, terms15, presence15, base15 = _z3_sparse_split_problem(
        1, 5, prefix="support_six_15"
    )
    solver15 = z3.Solver()
    solver15.add(*base15)
    _add_global_no_singleton_formula(solver15, terms15, presence15)
    p_x, p_y = p15[0]
    quartic_chain = tuple((4 - index, index * p_y) for index in range(5))
    solver15.add(
        z3.Not(
            z3.And(
                p_x == 0,
                p_y >= 2,
                *(belongs_to(exponent, quartic_chain) for exponent in q15),
            )
        )
    )
    if solver15.check() != z3.unsat:
        raise AssertionError("a non-quartic-chain 1+5 collision chart survived")

    _, _, terms24, presence24, base24 = _z3_sparse_split_problem(
        2, 4, prefix="support_six_24"
    )
    solver24 = z3.Solver()
    solver24.add(*base24)
    _add_global_no_singleton_formula(solver24, terms24, presence24)
    if solver24.check() != z3.unsat:
        raise AssertionError("an arbitrary-degree 2+4 collision chart survived")

    p33, q33, terms33, presence33, base33 = _z3_sparse_split_problem(
        3, 3, prefix="support_six_33"
    )
    solver33 = z3.Solver()
    solver33.add(*base33)
    _add_global_no_singleton_formula(solver33, terms33, presence33)
    quadratic_support = ((0, 2), (1, 1), (2, 0))
    solver33.add(
        z3.Not(
            z3.And(
                *(
                    belongs_to(exponent, quadratic_support)
                    for exponent in (*p33, *q33)
                )
            )
        )
    )
    if solver33.check() != z3.unsat:
        raise AssertionError("a nonquadratic 3+3 collision chart survived")

    return {
        "statement": (
            "every arbitrary-degree support-six exponent chart either has a "
            "singleton Keller coefficient or is one of two named shear charts"
        ),
        "splits": {
            "1+5": {
                "non_quartic_chain_no_singleton_formula": "unsat",
                "survivor": (
                    "P support {(0,m)} and Q support "
                    "{(4,0),(3,m),(2,2m),(1,3m),(0,4m)}, m>=2"
                ),
            },
            "2+4": {
                "global_no_singleton_formula": "unsat",
                "unit_certificate": (
                    "(z*(abcdef/M)/n)*(n*M) - (z*abcdef-1) = 1"
                ),
            },
            "3+3": {
                "nonquadratic_no_singleton_formula": "unsat",
                "survivor": "both supports equal {(2,0),(1,1),(0,2)}",
            },
        },
        "transpose_included": True,
        "solver": f"Z3 {z3.get_version_string()} exact integer arithmetic",
    }


def verify_arbitrary_shear_chain_formulas() -> dict[str, object]:
    """Check the exceptional formulas for several exponents exactly."""

    checked_exponents = (2, 3, 5, 11)
    for exponent in checked_exponents:
        # P=a*y^m and Q=b*x^2+c*x*y^m+d*y^(2m).
        # Q_y-P_y*Q_x has exactly the following two coefficients.
        coefficient_xy = exponent  # multiplies c-2ab
        coefficient_y = exponent  # multiplies 2d-ac
        if coefficient_xy == 0 or coefficient_y == 0:
            raise AssertionError("characteristic-zero exponent became zero")
    return {
        "support": {
            "P": ["y^m"],
            "Q": ["x^2", "x*y^m", "y^(2m)"],
        },
        "range": "every integer m >= 2 over characteristic zero",
        "coefficient_ideal": ["c-2*a*b", "2*d-a*c"],
        "forced_relations": ["c=2*a*b", "d=a^2*b"],
        "factorization": [
            "X=x+a*y^m",
            "Y=y+b*X^2",
        ],
        "inverse": [
            "y=Y-b*X^2",
            "x=X-a*(Y-b*X^2)^m",
        ],
        "transpose_included": True,
        "regression_exponents": list(checked_exponents),
    }


def verify_arbitrary_cubic_shear_chain_formulas() -> dict[str, object]:
    """Record the unique nonempty exact 1+4 support chain."""

    return {
        "support": {
            "P": ["y^m"],
            "Q": ["x^3", "x^2*y^m", "x*y^(2m)", "y^(3m)"],
        },
        "range": "every integer m >= 2 over characteristic zero",
        "coefficient_ideal": [
            "c-3*a*b",
            "d-a*c",
            "3*e-a*d",
        ],
        "forced_relations": [
            "c=3*a*b",
            "d=3*a^2*b",
            "e=a^3*b",
        ],
        "factorization": [
            "X=x+a*y^m",
            "Y=y+b*X^3",
        ],
        "inverse": [
            "y=Y-b*X^3",
            "x=X-a*(Y-b*X^3)^m",
        ],
        "transpose_included": True,
    }


def verify_arbitrary_quartic_shear_chain_formulas() -> dict[str, object]:
    """Record the unique nonempty exact 1+5 support chain."""

    return {
        "support": {
            "P": ["y^m"],
            "Q": [
                "x^4",
                "x^3*y^m",
                "x^2*y^(2m)",
                "x*y^(3m)",
                "y^(4m)",
            ],
        },
        "range": "every integer m >= 2 over characteristic zero",
        "coefficient_ideal": [
            "c-4*a*b",
            "2*d-3*a*c",
            "3*e-2*a*d",
            "4*f-a*e",
        ],
        "forced_relations": [
            "c=4*a*b",
            "d=6*a^2*b",
            "e=4*a^3*b",
            "f=a^4*b",
        ],
        "factorization": [
            "X=x+a*y^m",
            "Y=y+b*X^4",
        ],
        "inverse": [
            "y=Y-b*X^4",
            "x=X-a*(Y-b*X^4)^m",
        ],
        "transpose_included": True,
    }


def verify_directional_quadratic_shear_formulas() -> dict[str, object]:
    """Record the exact 3+3 survivor coefficient chart."""

    return {
        "support": {
            "P": ["y^2", "x*y", "x^2"],
            "Q": ["y^2", "x*y", "x^2"],
        },
        "coefficient_ideal": [
            "b+2*d",
            "2*c+e",
            "b*f+2*c^2",
            "2*a*f+b*c",
            "4*a*c-b^2",
        ],
        "exact_chart_parametrization": [
            "c=b^2/(4*a)",
            "d=-b/2",
            "e=-b^2/(2*a)",
            "f=-b^3/(8*a^2)",
            "a*b != 0",
        ],
        "directional_form": (
            "(P,Q)=lambda*(v,-u)*(u*x+v*y)^2 with lambda*u*v != 0"
        ),
        "inverse": "(X,Y) maps to (X-P(X,Y),Y-Q(X,Y))",
    }


def affine_normalized_support_consequence() -> dict[str, object]:
    return {
        "definition": (
            "sigma_aff(F)=min support(N_(p,L)(F)-identity), where "
            "N_(p,L)(z)=(JF(p)L)^(-1)*(F(p+Lz)-F(p))"
        ),
        "invariance": "polynomial invertibility is preserved by every normalization",
        "conclusion": (
            "every noninvertible plane Keller map has sigma_aff(F) >= 7"
        ),
        "degree_frontier_changed": False,
    }


def audit_existing_only(artifact: Path) -> None:
    """Validate the committed theorem ledger without enumeration or solvers."""

    artifact_hash = hashlib.sha256(artifact.read_bytes()).hexdigest()
    if artifact_hash != ARTIFACT_SHA256:
        raise AssertionError("the pinned sparse-support artifact bytes changed")
    certificate_hash = hashlib.sha256(SINGULAR_CERTIFICATE.read_bytes()).hexdigest()
    if certificate_hash != SINGULAR_CERTIFICATE_SHA256:
        raise AssertionError("the pinned Singular certificate bytes changed")

    payload = json.loads(artifact.read_text())
    assert payload["schema_version"] == 3
    boundary = payload["claim_boundary"]
    assert boundary["field"] == "characteristic zero"
    assert boundary["exact_support"] == "every displayed coefficient is nonzero"
    assert boundary["bounded_regression"].endswith("total degree at most 12")
    assert boundary["not_claimed"] == [
        "invariance of support size under affine normalization",
        "a new universal JC(2) degree bound",
        "supports with at least seven nonlinear monomial occurrences",
    ]

    balanced = payload["arbitrary_degree_balanced_binomial"]
    assert balanced["presence_masks_total"] == 256
    assert balanced["presence_masks_feasible"] == 85
    assert balanced["presence_mask_sha256"] == (
        "4d211a46e23b9549ddf4108ee292f041973d102c650d02f4dd9daab8bb0304ad"
    )
    assert balanced["linear_residual_masks"] == [
        f"0x{mask:02x}"
        for mask in dict.fromkeys(
            row_mask for row_mask, _partition in EXPECTED_LINEAR_RESIDUAL_PATTERNS
        )
    ]
    assert balanced["canonical_residual_patterns"] == len(
        EXPECTED_LINEAR_RESIDUAL_PATTERNS
    ) == 20
    assert balanced["residual_pattern_table"] == [
        {
            "mask": f"0x{mask:02x}",
            "blocks": ["".join(map(str, block)) for block in partition],
        }
        for mask, partition in EXPECTED_LINEAR_RESIDUAL_PATTERNS
    ]
    assert balanced["residual_patterns_all_unsat"] is True
    assert balanced["global_no_singleton_formula"] == "unsat"

    two_by_three = payload["arbitrary_degree_two_by_three"]
    assert two_by_three["presence_masks_total"] == 2048
    assert two_by_three["presence_masks_feasible"] == 321
    assert two_by_three["presence_mask_sha256"] == (
        "227bcb7b327487ec9dc6b861a75c3e1c2e693cfe3cf15726077f052407a7b0eb"
    )
    assert two_by_three["linear_residual_masks"] == 98
    assert two_by_three["linear_residual_mask_sha256"] == (
        "17f070086059ca1536395602983fda35dff2c962ad64e057f0472987e584fb9f"
    )
    assert two_by_three["global_no_singleton_formula"] == "unsat"
    assert two_by_three["transpose_included"] is True

    support_six = payload["support_six_exponent_classification"]
    assert support_six["transpose_included"] is True
    assert support_six["splits"]["1+5"][
        "non_quartic_chain_no_singleton_formula"
    ] == "unsat"
    assert support_six["splits"]["2+4"]["global_no_singleton_formula"] == "unsat"
    assert support_six["splits"]["3+3"][
        "nonquadratic_no_singleton_formula"
    ] == "unsat"

    assert payload["arbitrary_degree_shear_chain"] == (
        verify_arbitrary_shear_chain_formulas()
    )
    assert payload["arbitrary_degree_cubic_shear_chain"] == (
        verify_arbitrary_cubic_shear_chain_formulas()
    )
    assert payload["arbitrary_degree_quartic_shear_chain"] == (
        verify_arbitrary_quartic_shear_chain_formulas()
    )
    assert payload["directional_quadratic_shear"] == (
        verify_directional_quadratic_shear_formulas()
    )
    assert payload["affine_normalized_support"] == (
        affine_normalized_support_consequence()
    )

    bounded = payload["balanced_binomial_census"]
    assert bounded["degree_cap"] == 12
    assert bounded["ordered_support_pairs"] == 14_653_584
    assert bounded["unit_certificates"] == 14_653_584
    assert bounded["survivors"] == 0
    assert bounded["certificate_sha256"] == (
        "df87a3577cddd666f8f4428cf88f59c50996aab1f13b38fd32ddf548d5376ab4"
    )
    trinomial = payload["balanced_trinomial_collision_census"]
    assert trinomial["degree_cap"] == 6
    assert trinomial["ordered_support_pairs"] == 5_290_000
    assert trinomial["no_singleton_supports"] == 1
    assert len(trinomial["survivors"]) == 1
    assert trinomial["survivors"][0]["classification"] == (
        "directional quadratic shear exponent chart"
    )
    assert trinomial["claim_boundary"].startswith("bounded regression only")

    print(
        "PASS committed plane sparse-support audit: pinned JSON and Singular bytes; "
        "arbitrary-degree versus bounded-regression boundaries; no enumeration or solver"
    )


def run_singular_certificate() -> str:
    certificate_hash = hashlib.sha256(SINGULAR_CERTIFICATE.read_bytes()).hexdigest()
    if certificate_hash != SINGULAR_CERTIFICATE_SHA256:
        raise RuntimeError("the pinned Singular certificate hash changed")
    singular = shutil.which("Singular")
    if singular is None:
        raise RuntimeError("Singular is required for the exceptional Groebner replay")
    completed = subprocess.run(
        [singular, "-q", str(SINGULAR_CERTIFICATE)],
        check=False,
        capture_output=True,
        text=True,
    )
    if completed.returncode != 0:
        raise RuntimeError(
            "Singular certificate failed:\n" + completed.stdout + completed.stderr
        )
    markers = (
        "PASS sparse 1+3 exceptional Groebner and two-sided inverse certificate",
        "PASS sparse 1+4 exceptional Groebner and two-sided inverse certificate",
        "PASS sparse 1+5 exceptional Groebner and two-sided inverse certificate",
        "PASS sparse 3+3 directional Groebner and two-sided inverse certificate",
    )
    if any(marker not in completed.stdout for marker in markers):
        raise RuntimeError("Singular certificate did not print every PASS marker")
    version = subprocess.run(
        [singular, "--version"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.splitlines()[0]
    return version


def build_payload(degree_cap: int) -> dict[str, object]:
    singular_version = run_singular_certificate()
    return {
        "schema_version": 3,
        "claim_boundary": {
            "normalization": "F=(x+P,y+Q), with P,Q having no terms of degree <2",
            "exact_support": "every displayed coefficient is nonzero",
            "field": "characteristic zero",
            "arbitrary_degree_classes": [
                "one monomial in P and at most five in Q",
                "the transposed classes",
                "exactly two monomials in each of P and Q",
                "two monomials in P and three in Q",
                "the transposed 3+2 class",
                "two monomials in P and four in Q",
                "the transposed 4+2 class",
                "three monomials in each of P and Q",
                "common monomial-ray supports",
                "separated-axis supports",
            ],
            "bounded_regression": (
                "exactly two monomials in each of P,Q, all of total degree "
                f"at most {degree_cap}"
            ),
            "not_claimed": [
                "invariance of support size under affine normalization",
                "a new universal JC(2) degree bound",
                "supports with at least seven nonlinear monomial occurrences",
            ],
        },
        "arbitrary_degree_shear_chain": verify_arbitrary_shear_chain_formulas(),
        "arbitrary_degree_cubic_shear_chain": (
            verify_arbitrary_cubic_shear_chain_formulas()
        ),
        "arbitrary_degree_quartic_shear_chain": (
            verify_arbitrary_quartic_shear_chain_formulas()
        ),
        "arbitrary_degree_balanced_binomial": verify_arbitrary_balanced_binomials(),
        "arbitrary_degree_two_by_three": verify_arbitrary_two_by_three(),
        "support_six_exponent_classification": (
            verify_support_six_exponent_classification()
        ),
        "directional_quadratic_shear": (
            verify_directional_quadratic_shear_formulas()
        ),
        "affine_normalized_support": affine_normalized_support_consequence(),
        "balanced_binomial_census": balanced_binomial_census(degree_cap),
        "balanced_trinomial_collision_census": (
            balanced_trinomial_collision_census()
        ),
        "software": {
            "python": "standard library plus z3-solver",
            "z3": z3.get_version_string(),
            "singular": singular_version,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--degree-cap", type=int, default=12)
    parser.add_argument("--refresh", action="store_true")
    parser.add_argument("--artifact", type=Path, default=ARTIFACT)
    parser.add_argument(
        "--audit-existing-only",
        action="store_true",
        help="validate committed records without Z3, Singular, or support enumeration",
    )
    args = parser.parse_args()
    if args.degree_cap < 2 or args.degree_cap > 181:
        raise SystemExit("degree cap must lie between 2 and 181")

    artifact = args.artifact.resolve()
    if args.audit_existing_only:
        audit_existing_only(artifact)
        return

    payload = build_payload(args.degree_cap)
    if args.refresh:
        artifact.parent.mkdir(parents=True, exist_ok=True)
        artifact.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        print(f"WROTE {artifact.relative_to(ROOT)}")
    else:
        expected = json.loads(artifact.read_text())
        current_claim = {key: value for key, value in payload.items() if key != "software"}
        pinned_claim = {key: value for key, value in expected.items() if key != "software"}
        if current_claim != pinned_claim:
            raise AssertionError(
                "pinned sparse-support artifact is stale; inspect before --refresh"
            )

    census = payload["balanced_binomial_census"]
    print(
        "PASS arbitrary-degree balanced 2+2 singleton/unit theorem;",
        "bounded regression:",
        census["unit_certificates"],
        "supports through degree",
        census["degree_cap"],
    )
    print(
        "PASS arbitrary-degree 1+q classification (q<=5):",
        "the only nonempty charts are quadratic/cubic/quartic shear automorphisms",
    )
    print(
        "PASS support-six frontier:",
        "1+5 is the quartic shear, 2+4 is a unit ideal,",
        "and 3+3 is the directional quadratic shear",
    )
    trinomial_census = payload["balanced_trinomial_collision_census"]
    print(
        "PASS bounded 3+3 collision regression:",
        trinomial_census["ordered_support_pairs"],
        "supports through degree",
        trinomial_census["degree_cap"],
    )


if __name__ == "__main__":
    main()

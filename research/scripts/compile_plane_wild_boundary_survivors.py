#!/usr/bin/env python3
"""Compile the proved wild-boundary rows and the remaining search frontier.

The output deliberately separates four statuses:

* ``known_keller`` is reserved for the proved characteristic-two plane row;
* ``obstructed`` means that a proved different, unit, Euler, or class-group
  gate rejects affine-plane reconstruction;
* ``needs_reconstruction`` is a finite plane cover whose boundary ledger is
  known but whose candidate affine source has not been constructed; and
* ``local_comparison_only`` is an Artin--Schreier, Witt, or Kummer local row,
  not a global plane-cover reconstruction.

Passing the multiple-fibre packet gate is only a necessary condition.  It is
never reported as existence of a polynomial Keller map.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from collections import Counter
from itertools import combinations_with_replacement
from pathlib import Path
from typing import Any

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
CAS = ROOT / "plane-jc" / "cas"
sys.path.insert(0, str(CAS))

from boundary_lattice_prefilter import multiple_fiber_invariants


OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "plane_wild_boundary_survivor_atlas.json"
)

HIDDEN_ORDER_ROWS = (
    (2, 2),
    (2, 4),
    (2, 6),
    (3, 3),
    (3, 6),
    (5, 5),
    (7, 7),
)

THICKENED_ROWS = (
    (2, 2, 2),
    (2, 2, 4),
    (2, 2, 8),
    (3, 3, 2),
    (3, 3, 3),
    (5, 5, 2),
)

ODD_PRIMES = (3, 5, 7)
MAX_COVER_DEGREE = 15
MAX_PACKET_MULTIPLICITY = 12
MAX_PACKET_LENGTH = 3


def residue_parts(prime: int, degree: int) -> tuple[int, int]:
    """Split ``degree`` into prime-to-p and p-power factors."""

    inseparable = 1
    separable = degree
    while separable % prime == 0:
        separable //= prime
        inseparable *= prime
    return separable, inseparable


def conductor(kind: str, value: int | list[int] | None) -> dict[str, Any]:
    return {"kind": kind, "value": value}


def different(
    value: int | None,
    *,
    certification: str = "exact",
) -> dict[str, Any]:
    return {"certification": certification, "value": value}


def ledger(
    ramification_index: int,
    separable_residue_degree: int,
    inseparable_residue_degree: int,
    different_value: int | None,
    conductor_kind: str,
    conductor_value: int | list[int] | None,
    sheet_loss: int,
    *,
    different_certification: str = "exact",
) -> dict[str, Any]:
    """Return the normalized six-entry boundary ledger."""

    return {
        "conductor": conductor(conductor_kind, conductor_value),
        "different": different(
            different_value,
            certification=different_certification,
        ),
        "e": ramification_index,
        "f_insep": inseparable_residue_degree,
        "f_sep": separable_residue_degree,
        "sheet_loss": sheet_loss,
    }


def p_free_part(value: int, prime: int) -> tuple[int, int]:
    """Return ``(p_power, prime_to_p_part)``."""

    p_power = 1
    remainder = value
    while remainder % prime == 0:
        p_power *= prime
        remainder //= prime
    return p_power, remainder


def verify_balanced_support_identity() -> dict[str, Any]:
    """Replay the multi-retained support identity and coefficient kernel."""

    source_p, source_q, source_t = sp.symbols("P Q T")
    checks = 0
    admissible_degrees: dict[str, list[int]] = {}
    for prime in ODD_PRIMES:
        degree = prime
        allowed: list[int] = []
        for retained in range(1, MAX_COVER_DEGREE - degree + 1):
            boundary_factor = source_t**degree - source_p ** (degree + 1)
            retained_polynomial = 1 + source_t**retained
            order = (
                retained_polynomial * boundary_factor
                + source_p ** (degree - 1) * source_q * source_t
            )
            identity = (
                retained_polynomial * sp.diff(order, source_t)
                - sp.diff(retained_polynomial, source_t) * order
                - source_p ** (degree - 1)
                * source_q
                * (
                    retained_polynomial
                    - source_t * sp.diff(retained_polynomial, source_t)
                )
            )
            assert sp.Poly(
                sp.expand(identity),
                source_p,
                source_q,
                source_t,
                modulus=prime,
            ).is_zero
            kernel_exponents = [
                exponent
                for exponent in range(1, retained + 1)
                if (1 - exponent) % prime == 0
            ]
            assert kernel_exponents == list(range(1, retained + 1, prime))
            coefficient_kernel = retained in kernel_exponents
            assert coefficient_kernel == (retained % prime == 1)
            if coefficient_kernel:
                allowed.append(retained)
                support_factor = sp.Poly(
                    retained_polynomial
                    - source_t * sp.diff(retained_polynomial, source_t),
                    source_t,
                    modulus=prime,
                )
                assert support_factor.degree() == 0
                assert int(support_factor.LC()) % prime == 1
            checks += 1
        admissible_degrees[str(prime)] = allowed
    return {
        "admissible_retained_degrees_through_bound": admissible_degrees,
        "coefficient_characterization": "A=a0+T*B(T^p)",
        "identity": "A*H_T-A'*H=P^(N-1)*Q*(A-T*A')",
        "rows_checked": checks,
    }


def verify_balanced_root_count_identity() -> dict[str, Any]:
    """Replay the birational chart and the geometric point-count formula."""

    source_x, source_u = sp.symbols("x u")
    checks = 0
    for prime in ODD_PRIMES:
        degree = prime
        for retained in range(1, MAX_COVER_DEGREE - degree + 1, prime):
            target_p = source_x * source_u
            target_t = source_x**2 * source_u
            retained_polynomial = 1 + target_t**retained
            target_q = retained_polynomial * (
                source_u - source_x ** (degree - 1)
            )
            hidden_order = (
                retained_polynomial
                * (target_t**degree - target_p ** (degree + 1))
                + target_p ** (degree - 1) * target_q * target_t
            )
            assert sp.Poly(
                sp.expand(hidden_order),
                source_x,
                source_u,
                modulus=prime,
            ).is_zero
            checks += 1

    field_order, root_count = sp.symbols(
        "q n", integer=True, positive=True
    )
    chart_points = field_order**2 - root_count * (field_order - 1)
    root_fibre_points = root_count * (2 * field_order - 1)
    cover_points = sp.expand(chart_points + root_fibre_points)
    open_points = sp.expand(cover_points - field_order)
    assert sp.expand(
        cover_points - (field_order**2 + root_count * field_order)
    ) == 0
    assert sp.expand(
        open_points - (field_order**2 + (root_count - 1) * field_order)
    ) == 0
    return {
        "boundary_count": "q",
        "chart": (
            "P=x*u,T=x^2*u,Q=A(x^2*u)*(u-x^(N-1)); "
            "isomorphism over D(A)"
        ),
        "chart_count": "q^2-n_q(A)*(q-1)",
        "cover_count": "q^2+n_q(A)*q",
        "geometric_consequence": (
            "after a finite splitting field n_q(A)=deg(A), so C-E is "
            "geometrically A2 only if deg(A)=1"
        ),
        "open_count": "q^2+(n_q(A)-1)*q",
        "root_fibre_count": "2*q-1 per rational root",
        "rows_checked": checks,
    }


def hidden_order_row(
    prime: int,
    degree: int,
    exponent_a: int,
    exponent_b: int = 0,
) -> dict[str, Any]:
    """Compile one proved monomial hidden-order row."""

    if degree % prime:
        raise ValueError("hidden order requires p dividing N")
    if exponent_a < 0 or exponent_b < 0:
        raise ValueError("gluing exponents must be nonnegative")

    separable, inseparable = residue_parts(prime, degree)
    q_exponent = exponent_b + 1
    fierce = ledger(
        1,
        separable,
        inseparable,
        q_exponent,
        "boundary_normalization_exponent",
        degree * (degree - 2),
        degree,
    )
    row: dict[str, Any] = {
        "architecture": "target_supported_hidden_order",
        "characteristic": prime,
        "cover_degree": degree + 1,
        "fierce_ledger": fierce,
        "gluing": {"a": exponent_a, "b": exponent_b},
        "id": f"hidden-p{prime}-N{degree}-a{exponent_a}-b{exponent_b}",
        "omitted_boundary_count": 1,
        "primitive_order": (
            "(T-1)*(T^N-P^(N+1))+P^a*Q^(b+1)*T"
        ),
        "wild_degree": degree,
    }

    if exponent_b > 0:
        if exponent_a != degree - 1:
            raise ValueError("stored thickened controls are balanced")
        p_power, tame_degree = p_free_part(q_exponent, prime)
        row["base_change"] = {
            "degree": q_exponent,
            "purely_inseparable_degree": p_power,
            "separable_degree": tame_degree,
        }
        if tame_degree > 1:
            euler = (degree + 1) * tame_degree - degree
            row.update(
                {
                    "decision": "obstructed",
                    "decisive_gate": "compactly_supported_euler_characteristic",
                    "obstruction": {"chi_c": euler, "required_for_A2": 1},
                }
            )
        elif degree > 2:
            row.update(
                {
                    "decision": "obstructed",
                    "decisive_gate": "purely_inseparable_class_transfer",
                    "obstruction": {
                        "class_group_contains": f"Z/{degree - 1}",
                        "named_class": "L1",
                        "named_class_order": degree - 1,
                    },
                }
            )
        else:
            profile = multiple_fiber_invariants(
                (2, 2, 2),
                reduced_sum_principal=True,
                generic_unit_rank=1,
                generic_class_group_trivial=True,
                other_vertical_classes_trivial=True,
            )
            row.update(
                {
                    "decision": "obstructed",
                    "decisive_gate": "full_source_class_group",
                    "obstruction": {
                        "class_group": [2, 2],
                        "core": "D(x*u)",
                        "multiple_fiber_multiplicities": [2, 2, 2],
                        "smith_diagonal": list(profile.smith_diagonal),
                        "source_fill_matrix": [[1, 0], [-1, 1]],
                    },
                }
            )
        return row

    if exponent_a == 0:
        row.update(
            {
                "decision": "obstructed",
                "decisive_gate": "unit_and_class_group",
                "obstruction": {
                    "class_group": f"Z/{degree + 1}",
                    "unit_rank": 1,
                },
            }
        )
    elif exponent_a < degree - 1:
        common = math.gcd(exponent_a, degree - 1)
        tame_index = (degree - 1) // common
        companion = ledger(
            tame_index,
            common,
            1,
            tame_index - 1,
            "none",
            None,
            degree - 1,
        )
        row.update(
            {
                "companion_ledger": companion,
                "decision": "obstructed",
                "decisive_gate": "companion_tame_different",
                "obstruction": {
                    "different": tame_index - 1,
                    "ramification_index": tame_index,
                },
            }
        )
    elif exponent_a == degree - 1:
        row["companion_ledgers"] = [
            ledger(1, 1, 1, 0, "none", None, 1),
            ledger(1, degree - 1, 1, 0, "none", None, degree - 1),
        ]
        if degree == 2:
            row.update(
                {
                    "decision": "known_keller",
                    "decisive_gate": "all_known_gates_pass",
                    "reconstruction": "characteristic_two_affine_plane",
                }
            )
        else:
            row.update(
                {
                    "decision": "obstructed",
                    "decisive_gate": "balanced_source_class_group",
                    "obstruction": {
                        "class_group": f"Z/{degree - 1}",
                        "named_class": "L1",
                        "named_class_order": degree - 1,
                    },
                }
            )
    else:
        row.update(
            {
                "companion_ledger": ledger(
                    degree,
                    1,
                    1,
                    None,
                    "none",
                    None,
                    degree,
                    different_certification="proved_positive_not_computed",
                ),
                "decision": "obstructed",
                "decisive_gate": "wild_ramification_over_P",
                "obstruction": {
                    "different": "positive",
                    "ramification_index": degree,
                },
            }
        )
    return row


def hidden_order_rows() -> list[dict[str, Any]]:
    """Compile one representative of every proved monomial regime."""

    rows: list[dict[str, Any]] = []
    for prime, degree in HIDDEN_ORDER_ROWS:
        for exponent_a in range(degree + 1):
            rows.append(hidden_order_row(prime, degree, exponent_a))
    rows.extend(
        hidden_order_row(prime, degree, degree - 1, q_exponent - 1)
        for prime, degree, q_exponent in THICKENED_ROWS
    )
    return rows


def prescribed_degree_rows() -> list[dict[str, Any]]:
    """Compile the theorem-backed odd-characteristic prime-degree queue."""

    rows: list[dict[str, Any]] = []
    for prime in ODD_PRIMES:
        degree = prime
        separable, inseparable = residue_parts(prime, degree)
        for cover_degree in range(degree + 1, MAX_COVER_DEGREE + 1):
            if cover_degree % prime == 0:
                continue
            retained = cover_degree - degree
            row = {
                "architecture": "prescribed_degree_hidden_order",
                "characteristic": prime,
                "cover_degree": cover_degree,
                "fierce_ledger": ledger(
                    1,
                    separable,
                    inseparable,
                    1,
                    "boundary_normalization_exponent",
                    degree * (degree - 2),
                    degree,
                ),
                "id": f"prescribed-p{prime}-N{degree}-d{cover_degree}",
                "node_count": retained * (degree + 1),
                "retained_sheet_count": retained,
                "wild_degree": degree,
            }
            if retained == 1:
                row.update(
                    {
                        "decision": "obstructed",
                        "decisive_gate": "linear_retained_case_companion_different",
                        "obstruction": {"different": degree - 2},
                    }
                )
            else:
                row.update(
                    {
                        "decision": "needs_reconstruction",
                        "missing_certificates": [
                            "normalization of the proposed source open",
                            "complete different away from the omitted boundary",
                            "generic-fibre unit and class groups",
                            "source-fill valuation matrix",
                            "constant-Jacobian coordinates",
                        ],
                    }
                )
            rows.append(row)
    return rows


def balanced_prescribed_degree_rows() -> list[dict[str, Any]]:
    """Compile the exact support sieve for balanced multi-retained orders.

    For ``H=A(T)F+P^(N-1)QT`` one has
    ``A*H_T-A'*H=P^(N-1)Q(A-TA')``.  With ``A`` monic, squarefree, and
    ``A(0) != 0``, target different support inside ``V(PQ)`` therefore
    forces ``A-TA'`` to be the nonzero constant ``A(0)``.  In characteristic
    ``p`` this is equivalent to ``A=A(0)+T*B(T^p)``.
    """

    rows: list[dict[str, Any]] = []
    for prime in ODD_PRIMES:
        degree = prime
        separable, inseparable = residue_parts(prime, degree)
        for cover_degree in range(degree + 1, MAX_COVER_DEGREE + 1):
            if cover_degree % prime == 0:
                continue
            retained = cover_degree - degree
            row: dict[str, Any] = {
                "architecture": "balanced_prescribed_degree_hidden_order",
                "characteristic": prime,
                "cover_degree": cover_degree,
                "equation": "A(T)*(T^N-P^(N+1))+P^(N-1)*Q*T",
                "fierce_ledger": ledger(
                    1,
                    separable,
                    inseparable,
                    1,
                    "boundary_normalization_exponent",
                    degree * (degree - 2),
                    degree,
                ),
                "id": f"balanced-prescribed-p{prime}-N{degree}-d{cover_degree}",
                "retained_polynomial_degree": retained,
                "support_identity": "A*H_T-A'*H=P^(N-1)*Q*(A-T*A')",
                "wild_degree": degree,
            }
            if retained % prime != 1:
                row.update(
                    {
                        "decision": "obstructed",
                        "decisive_gate": "extra_target_different",
                        "obstruction": {
                            "leading_coefficient_of_A_minus_TA_prime_mod_p": (
                                1 - retained
                            )
                            % prime,
                            "reason": "A-T*A' cannot be constant",
                        },
                    }
                )
            elif retained == 1:
                row.update(
                    {
                        "admissible_polynomial_shape": "A=a0+a1*T",
                        "decision": "obstructed",
                        "decisive_gate": "balanced_source_class_group",
                        "obstruction": {
                            "class_group": f"Z/{degree - 1}",
                            "named_class_order": degree - 1,
                        },
                    }
                )
            else:
                row.update(
                    {
                        "admissible_polynomial_shape": "A=a0+T*B(T^p)",
                        "canonical_squarefree_control": f"A=1+T^{retained}",
                        "decision": "obstructed",
                        "decisive_gate": (
                            "geometric_point_count_after_retained_roots_split"
                        ),
                        "obstruction": {
                            "open_count_over_Fq": "q^2+(n_q(A)-1)*q",
                            "open_count_over_splitting_field": (
                                f"q^2+{retained - 1}*q"
                            ),
                            "retained_root_count_over_splitting_field": retained,
                        },
                    }
                )
            rows.append(row)
    return rows


def artin_schreier_witt_different(
    prime: int,
    jumps: tuple[int, ...],
) -> int:
    return sum(
        (prime**index - prime ** (index - 1)) * (jump + 1)
        for index, jump in enumerate(jumps, start=1)
    )


def comparison_rows() -> list[dict[str, Any]]:
    """Compile representative additive, AS, AS--Witt, and Kummer rows."""

    rows: list[dict[str, Any]] = []
    for prime in ODD_PRIMES:
        for derivative_order in (1, 2):
            rows.append(
                {
                    "architecture": "fierce_additive_block",
                    "characteristic": prime,
                    "cover_degree": prime,
                    "decision": "needs_reconstruction",
                    "equation": "Z^p+Q^m*Z-P",
                    "id": f"additive-p{prime}-m{derivative_order}",
                    "ledger": ledger(
                        1,
                        1,
                        prime,
                        derivative_order,
                        "none",
                        None,
                        prime,
                    ),
                    "missing_certificates": [
                        "global omitted-boundary presentation",
                        "affine-plane source reconstruction",
                        "determinant cancellation ledger",
                    ],
                }
            )
        for pole_order in (1, 2):
            if pole_order % prime == 0:
                continue
            rows.append(
                {
                    "architecture": "artin_schreier_local",
                    "characteristic": prime,
                    "decision": "local_comparison_only",
                    "equation": "Y^p-Y=c*Q^(-m)",
                    "id": f"as-p{prime}-m{pole_order}",
                    "ledger": ledger(
                        prime,
                        1,
                        1,
                        (prime - 1) * (pole_order + 1),
                        "artin",
                        pole_order + 1,
                        prime,
                    ),
                }
            )
        jumps = (1, prime + 1)
        rows.append(
            {
                "architecture": "artin_schreier_witt_local",
                "characteristic": prime,
                "decision": "local_comparison_only",
                "id": f"asw-p{prime}-j1-{prime + 1}",
                "jumps": list(jumps),
                "ledger": ledger(
                    prime**2,
                    1,
                    1,
                    artin_schreier_witt_different(prime, jumps),
                    "artin_witt_vector",
                    [jump + 1 for jump in jumps],
                    prime**2,
                ),
            }
        )
        for index in range(2, 6):
            if math.gcd(index, prime) != 1:
                continue
            rows.append(
                {
                    "architecture": "tame_kummer_local",
                    "characteristic": prime,
                    "decision": "local_comparison_only",
                    "id": f"kummer-p{prime}-e{index}",
                    "ledger": ledger(
                        index,
                        1,
                        1,
                        index - 1,
                        "none",
                        None,
                        index,
                    ),
                }
            )
    return rows


def packet_scan() -> dict[str, Any]:
    """Enumerate bounded odd-characteristic packets surviving vertical torsion."""

    by_prime: dict[str, Any] = {}
    for prime in ODD_PRIMES:
        survivors: list[dict[str, Any]] = []
        rejected = 0
        tested = 0
        torsion_orders: Counter[int] = Counter()
        for length in range(2, MAX_PACKET_LENGTH + 1):
            for multiplicities in combinations_with_replacement(
                range(2, MAX_PACKET_MULTIPLICITY + 1),
                length,
            ):
                if not any(value % prime == 0 for value in multiplicities):
                    continue
                tested += 1
                profile = multiple_fiber_invariants(
                    multiplicities,
                    reduced_sum_principal=True,
                    generic_unit_rank=1,
                )
                torsion_order = math.prod(profile.class_torsion)
                if torsion_order > 1:
                    rejected += 1
                    torsion_orders[torsion_order] += 1
                    assert any(
                        math.gcd(left, right) > 1
                        for index, left in enumerate(multiplicities)
                        for right in multiplicities[index + 1 :]
                    )
                    continue
                assert all(
                    math.gcd(left, right) == 1
                    for index, left in enumerate(multiplicities)
                    for right in multiplicities[index + 1 :]
                )
                survivors.append(
                    {
                        "multiplicities": list(multiplicities),
                        "smith_diagonal": list(profile.smith_diagonal),
                        "status": "packet_gate_only",
                    }
                )
        survivors.sort(
            key=lambda row: (
                sum(row["multiplicities"]),
                len(row["multiplicities"]),
                row["multiplicities"],
            )
        )
        by_prime[str(prime)] = {
            "minimal_survivors": survivors[:20],
            "rejected_by_vertical_torsion": rejected,
            "survivor_count": len(survivors),
            "tested_count": tested,
            "torsion_order_distribution": {
                str(order): count
                for order, count in sorted(torsion_orders.items())
            },
        }
        assert rejected + len(survivors) == tested
    return {
        "bounds": {
            "max_length": MAX_PACKET_LENGTH,
            "max_multiplicity": MAX_PACKET_MULTIPLICITY,
        },
        "by_characteristic": by_prime,
        "interpretation": (
            "These are only packets not rejected by the vertical class-group "
            "gate. Pairwise coprimality is necessary and sufficient inside "
            "the certified principal-sum, generic-unit-rank-one model."
        ),
    }


def characteristic_zero_module_template() -> dict[str, Any]:
    """Compile the class-presentation input preceding the Case-1 residue."""

    packet_presentations = []
    for multiplicities in ((2, 3), (2, 5), (2, 7)):
        profile = multiple_fiber_invariants(
            multiplicities,
            reduced_sum_principal=True,
            generic_unit_rank=1,
        )
        packet_presentations.append(
            {
                "multiplicities": list(multiplicities),
                "relation_matrix_R": [list(row) for row in profile.relation_matrix],
                "smith_diagonal": list(profile.smith_diagonal),
                "vertical_class_group": list(profile.class_torsion),
            }
        )
        assert not profile.class_torsion
    return {
        "class_group_block": "[[V,A],[0,R]]",
        "finite_support_residue": (
            "H^0_Z(coker(Phi))=(im(Phi):I_Z^infinity)/im(Phi)"
        ),
        "minimal_packet_presentations": packet_presentations,
        "required_order_of_gates": [
            "different vanishes away from the omitted boundary",
            "compiled determinant class vanishes in Cl(U) or Pic(U)",
            "finite-support local-cohomology residue vanishes",
            "affineness and polynomial-coordinate reconstruction",
        ],
        "status": "coherent_template_not_instantiated",
    }


def compile_atlas() -> dict[str, Any]:
    support_certificate = verify_balanced_support_identity()
    root_count_certificate = verify_balanced_root_count_identity()
    hidden = hidden_order_rows()
    prescribed = prescribed_degree_rows()
    balanced_prescribed = balanced_prescribed_degree_rows()
    comparisons = comparison_rows()
    packets = packet_scan()

    hidden_decisions = Counter(row["decision"] for row in hidden)
    prescribed_decisions = Counter(row["decision"] for row in prescribed)
    balanced_prescribed_decisions = Counter(
        row["decision"] for row in balanced_prescribed
    )
    comparison_decisions = Counter(row["decision"] for row in comparisons)

    known = [row for row in hidden if row["decision"] == "known_keller"]
    assert len(known) == 1
    assert known[0]["characteristic"] == 2
    assert known[0]["wild_degree"] == 2
    assert not any(
        row["decision"] == "known_keller" and row["characteristic"] % 2
        for row in hidden
    )
    assert all(
        row["cover_degree"] % row["characteristic"] != 0
        for row in prescribed
    )
    support_only_queue = [
        row
        for row in balanced_prescribed
        if row["retained_polynomial_degree"] > 1
        if row["retained_polynomial_degree"] % row["characteristic"] == 1
    ]
    assert [
        (row["characteristic"], row["cover_degree"])
        for row in support_only_queue
    ] == [(3, 7), (3, 10), (3, 13), (5, 11), (7, 15)]
    assert all(row["decision"] == "obstructed" for row in support_only_queue)
    support_only_rows = []
    for row in support_only_queue:
        support_only_rows.append(
            {
                "candidate_id": row["id"],
                "characteristic": row["characteristic"],
                "cover_degree": row["cover_degree"],
                "decision": "obstructed",
                "former_status": "support_only_survivor",
                "obstruction": (
                    "geometric point count after the retained polynomial "
                    "splits"
                ),
                "retained_polynomial_degree": row[
                    "retained_polynomial_degree"
                ],
            }
        )

    return {
        "balanced_prescribed_degree_rows": balanced_prescribed,
        "balanced_root_count_certificate": root_count_certificate,
        "balanced_support_certificate": support_certificate,
        "characteristic_zero_boundary_module": (
            characteristic_zero_module_template()
        ),
        "comparison_rows": comparisons,
        "format": "plane-wild-boundary-survivor-atlas-v2",
        "hidden_order_rows": hidden,
        "odd_characteristic_packet_scan": packets,
        "odd_characteristic_reconstruction_queue": [],
        "odd_characteristic_support_only_rows": support_only_rows,
        "prescribed_degree_rows": prescribed,
        "scope": (
            "Exact compilation of the proved target-supported hidden-order "
            "rows, theorem-backed prescribed-degree covers through the stated "
            "bound, representative additive/AS/AS-Witt/Kummer comparison "
            "rows, the geometric retained-root point-count obstruction, and "
            "a bounded abstract multiple-fibre packet scan. No nonlinear "
            "balanced retained-polynomial row survives affine-plane "
            "reconstruction. A bounded packet survivor is not a cover and "
            "not a Keller map."
        ),
        "summary": {
            "balanced_prescribed_degree_decisions": dict(
                sorted(balanced_prescribed_decisions.items())
            ),
            "balanced_prescribed_degree_row_count": len(balanced_prescribed),
            "comparison_decisions": dict(sorted(comparison_decisions.items())),
            "comparison_row_count": len(comparisons),
            "hidden_decisions": dict(sorted(hidden_decisions.items())),
            "hidden_row_count": len(hidden),
            "odd_hidden_keller_rows": 0,
            "prescribed_degree_decisions": dict(
                sorted(prescribed_decisions.items())
            ),
            "prescribed_degree_row_count": len(prescribed),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--write",
        action="store_true",
        help="rewrite the pinned generated artifact",
    )
    parser.add_argument(
        "--emit-json",
        action="store_true",
        help="print the compiled atlas instead of checking the artifact",
    )
    args = parser.parse_args()

    artifact = compile_atlas()
    serialized = json.dumps(artifact, indent=2, sort_keys=True) + "\n"
    if args.emit_json:
        print(serialized, end="")
        return
    if args.write:
        OUTPUT.write_text(serialized)
    else:
        assert OUTPUT.exists(), (
            f"missing {OUTPUT.relative_to(ROOT)}; regenerate with --write"
        )
        assert OUTPUT.read_text() == serialized, (
            f"{OUTPUT.relative_to(ROOT)} is stale; regenerate with --write"
        )

    digest = hashlib.sha256(serialized.encode()).hexdigest()
    print(
        "PASS plane wild-boundary survivor atlas: "
        f"hidden={artifact['summary']['hidden_row_count']}, "
        f"prescribed={artifact['summary']['prescribed_degree_row_count']}, "
        "balanced_prescribed="
        f"{artifact['summary']['balanced_prescribed_degree_row_count']}, "
        f"comparisons={artifact['summary']['comparison_row_count']}"
    )
    for prime, profile in artifact["odd_characteristic_packet_scan"][
        "by_characteristic"
    ].items():
        print(
            f"PASS p={prime} packet sieve: tested={profile['tested_count']}, "
            f"rejected={profile['rejected_by_vertical_torsion']}, "
            f"survivors={profile['survivor_count']}"
        )
    print(f"SHA256: {digest}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Shared retained-root Euler gate for boundary-search compilers.

The gate applies only to an explicitly certified balanced presentation

    A(T)*(T^N-P^(N+1)) + P^(N-1)*Q*T

with squarefree ``A``, ``A(0) != 0``, and exactly one omitted fierce
boundary.  It never infers those hypotheses from a coarse branch ledger.
When applicable, the normalized complement has geometric compactly supported
Euler characteristic ``deg(A)`` and therefore cannot be an affine plane when
``deg(A) > 1``.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from enum import Enum
from typing import Any


class RetainedRootEulerStatus(str, Enum):
    """Four-way result preserving inapplicability and missing proof data."""

    NOT_APPLICABLE = "not_applicable"
    UNCERTIFIED = "uncertified"
    PASSES = "passes"
    OBSTRUCTED = "obstructed"


@dataclass(frozen=True)
class RetainedRootEulerDatum:
    """Proof-bearing input for the geometric retained-root theorem."""

    retained_degree: int
    squarefree: bool
    nonzero_constant_term: bool
    omitted_fierce_boundary_count: int
    balanced_chart_certificate: str
    different_support_certificate: str
    root_fibre_certificate: str
    omitted_boundary_certificate: str
    theorem_source: str

    def __post_init__(self) -> None:
        if (
            isinstance(self.retained_degree, bool)
            or not isinstance(self.retained_degree, int)
            or self.retained_degree <= 0
        ):
            raise ValueError("retained degree must be positive")
        if not isinstance(self.squarefree, bool) or not isinstance(
            self.nonzero_constant_term, bool
        ):
            raise ValueError("squarefree and nonzero_constant_term must be booleans")
        if (
            isinstance(self.omitted_fierce_boundary_count, bool)
            or not isinstance(self.omitted_fierce_boundary_count, int)
            or self.omitted_fierce_boundary_count < 0
        ):
            raise ValueError("omitted-boundary count must be nonnegative")
        certificate_values = (
            self.balanced_chart_certificate,
            self.different_support_certificate,
            self.root_fibre_certificate,
            self.omitted_boundary_certificate,
            self.theorem_source,
        )
        if not all(isinstance(value, str) for value in certificate_values):
            raise ValueError("retained-root certificate references must be strings")

    @property
    def certificates_complete(self) -> bool:
        return all(
            value.strip()
            for value in (
                self.balanced_chart_certificate,
                self.different_support_certificate,
                self.root_fibre_certificate,
                self.omitted_boundary_certificate,
                self.theorem_source,
            )
        )


@dataclass(frozen=True)
class RetainedRootEulerAudit:
    """Machine-readable result of the retained-root Euler gate."""

    status: RetainedRootEulerStatus
    applicable: bool
    obstruction: bool
    retained_degree: int | None
    geometric_euler_characteristic: int | None
    expected_affine_plane_euler_characteristic: int
    geometric_open_class: str | None
    finite_field_open_count: str | None
    reason: str
    theorem_source: str | None

    @property
    def allows_continuation(self) -> bool:
        """Only a proved obstruction stops the downstream search."""

        return not self.obstruction

    def to_dict(self) -> dict[str, Any]:
        result = asdict(self)
        result["status"] = self.status.value
        return result


def audit_retained_root_euler(
    datum: RetainedRootEulerDatum | None,
) -> RetainedRootEulerAudit:
    """Evaluate the theorem without promoting uncertified hypotheses."""

    if datum is None:
        return RetainedRootEulerAudit(
            status=RetainedRootEulerStatus.NOT_APPLICABLE,
            applicable=False,
            obstruction=False,
            retained_degree=None,
            geometric_euler_characteristic=None,
            expected_affine_plane_euler_characteristic=1,
            geometric_open_class=None,
            finite_field_open_count=None,
            reason="no retained-polynomial presentation was declared",
            theorem_source=None,
        )
    if (
        not datum.squarefree
        or not datum.nonzero_constant_term
        or datum.omitted_fierce_boundary_count != 1
    ):
        failed = []
        if not datum.squarefree:
            failed.append("the retained polynomial is not certified squarefree")
        if not datum.nonzero_constant_term:
            failed.append("the retained constant term may vanish")
        if datum.omitted_fierce_boundary_count != 1:
            failed.append("the presentation does not omit exactly one fierce boundary")
        return RetainedRootEulerAudit(
            status=RetainedRootEulerStatus.NOT_APPLICABLE,
            applicable=False,
            obstruction=False,
            retained_degree=datum.retained_degree,
            geometric_euler_characteristic=None,
            expected_affine_plane_euler_characteristic=1,
            geometric_open_class=None,
            finite_field_open_count=None,
            reason="; ".join(failed),
            theorem_source=datum.theorem_source or None,
        )
    if not datum.certificates_complete:
        return RetainedRootEulerAudit(
            status=RetainedRootEulerStatus.UNCERTIFIED,
            applicable=False,
            obstruction=False,
            retained_degree=datum.retained_degree,
            geometric_euler_characteristic=None,
            expected_affine_plane_euler_characteristic=1,
            geometric_open_class=None,
            finite_field_open_count=None,
            reason=(
                "balanced chart, different support, root fibres, omitted "
                "boundary, and theorem source must all be certified"
            ),
            theorem_source=datum.theorem_source or None,
        )

    degree = datum.retained_degree
    obstruction = degree > 1
    return RetainedRootEulerAudit(
        status=(
            RetainedRootEulerStatus.OBSTRUCTED
            if obstruction
            else RetainedRootEulerStatus.PASSES
        ),
        applicable=True,
        obstruction=obstruction,
        retained_degree=degree,
        geometric_euler_characteristic=degree,
        expected_affine_plane_euler_characteristic=1,
        geometric_open_class=f"L^2+{degree - 1}*L",
        finite_field_open_count="q^2+(n_q(A)-1)*q",
        reason=(
            f"geometric Euler characteristic {degree} differs from 1"
            if obstruction
            else "the Euler gate is neutral in retained degree one"
        ),
        theorem_source=datum.theorem_source,
    )


def retained_root_euler_datum_from_dict(
    data: dict[str, object],
) -> RetainedRootEulerDatum:
    """Parse the shared JSON input block used by the log-boundary compiler."""

    retained_degree = data.get("retained_degree")
    squarefree = data.get("squarefree", False)
    nonzero_constant_term = data.get("nonzero_constant_term", False)
    omitted_count = data.get("omitted_fierce_boundary_count", 0)
    if isinstance(retained_degree, bool) or not isinstance(retained_degree, int):
        raise ValueError("retained_degree must be an integer")
    if not isinstance(squarefree, bool) or not isinstance(
        nonzero_constant_term, bool
    ):
        raise ValueError("retained-root hypothesis flags must be booleans")
    if isinstance(omitted_count, bool) or not isinstance(omitted_count, int):
        raise ValueError("omitted_fierce_boundary_count must be an integer")

    certificate_fields = (
        "balanced_chart_certificate",
        "different_support_certificate",
        "root_fibre_certificate",
        "omitted_boundary_certificate",
        "theorem_source",
    )
    certificate_values = {
        field: data.get(field, "") for field in certificate_fields
    }
    if not all(
        isinstance(value, str) for value in certificate_values.values()
    ):
        raise ValueError("retained-root certificate references must be strings")
    return RetainedRootEulerDatum(
        retained_degree=retained_degree,
        squarefree=squarefree,
        nonzero_constant_term=nonzero_constant_term,
        omitted_fierce_boundary_count=omitted_count,
        balanced_chart_certificate=certificate_values[
            "balanced_chart_certificate"
        ],
        different_support_certificate=certificate_values[
            "different_support_certificate"
        ],
        root_fibre_certificate=certificate_values["root_fibre_certificate"],
        omitted_boundary_certificate=certificate_values[
            "omitted_boundary_certificate"
        ],
        theorem_source=certificate_values["theorem_source"],
    )


__all__ = [
    "RetainedRootEulerAudit",
    "RetainedRootEulerDatum",
    "RetainedRootEulerStatus",
    "audit_retained_root_euler",
    "retained_root_euler_datum_from_dict",
]

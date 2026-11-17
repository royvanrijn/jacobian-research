#!/usr/bin/env python3
"""Proof-bearing conductor/contact-loss jet truncation audit.

For a finite normalization ``A -> B`` of a reduced curve, the conductor
annihilates ``Q=B/A``.  On a completed normalization branch with parameter
``t`` and conductor exponent ``c``, a section of ``Q`` therefore depends
only on its image modulo ``t^c``.  If the expression producing that section
can lose at most ``lambda`` orders of contact, input jets modulo
``t^(c+lambda)`` suffice.

This module records the safe additive bound

    lambda = differential_order + pole_order + additional_contact_loss

and never infers any of its terms from a coarse Newton polygon.  It is a
truncation audit, not a vanishing or existence test for the resulting
local-cohomology class.

The dependency-sensitive interface propagates a separate loss along every
named input-to-output path.  Its optional valuation adapter converts a
certified valuation frontier of omitted Newton support into the corresponding
available normal-jet order.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from enum import Enum
from typing import Any, Iterable


class ConductorJetStatus(str, Enum):
    """Four-way result which never promotes missing data to a proof."""

    NOT_DECLARED = "not_declared"
    UNCERTIFIED = "uncertified"
    PASSES = "passes"
    INSUFFICIENT = "insufficient"


@dataclass(frozen=True)
class ConductorBranchJetDatum:
    """One completed normalization-branch contact ledger.

    ``available_jet_order=n`` means that the inputs are known modulo
    ``t^n``.  Thus the known coefficients have exponents ``0,...,n-1``.
    """

    name: str
    conductor_exponent: int
    differential_order: int
    pole_order: int
    additional_contact_loss: int
    available_jet_order: int | None
    conductor_certificate: str
    expression_certificate: str
    valuation_certificate: str

    def __post_init__(self) -> None:
        if not isinstance(self.name, str) or not self.name.strip():
            raise ValueError("a conductor branch needs a nonempty name")
        for field_name in (
            "conductor_exponent",
            "differential_order",
            "pole_order",
            "additional_contact_loss",
        ):
            value = getattr(self, field_name)
            if isinstance(value, bool) or not isinstance(value, int) or value < 0:
                raise ValueError(f"{field_name} must be a nonnegative integer")
        if self.conductor_exponent == 0:
            raise ValueError(
                "only branches in the finite conductor support should be listed"
            )
        if self.available_jet_order is not None and (
            isinstance(self.available_jet_order, bool)
            or not isinstance(self.available_jet_order, int)
            or self.available_jet_order < 0
        ):
            raise ValueError(
                "available_jet_order must be null or a nonnegative integer"
            )
        certificates = (
            self.conductor_certificate,
            self.expression_certificate,
            self.valuation_certificate,
        )
        if not all(isinstance(value, str) for value in certificates):
            raise ValueError("conductor-jet certificate references must be strings")

    @property
    def contact_loss(self) -> int:
        return (
            self.differential_order
            + self.pole_order
            + self.additional_contact_loss
        )

    @property
    def required_jet_order(self) -> int:
        return self.conductor_exponent + self.contact_loss

    @property
    def certificates_complete(self) -> bool:
        return all(
            value.strip()
            for value in (
                self.conductor_certificate,
                self.expression_certificate,
                self.valuation_certificate,
            )
        )

    @property
    def margin(self) -> int | None:
        if self.available_jet_order is None:
            return None
        return self.available_jet_order - self.required_jet_order

    def to_dict(self) -> dict[str, Any]:
        result = asdict(self)
        result.update(
            contact_loss=self.contact_loss,
            required_jet_order=self.required_jet_order,
            margin=self.margin,
        )
        return result


@dataclass(frozen=True)
class ConductorJetAudit:
    """Combined truncation result over every conductor branch."""

    status: ConductorJetStatus
    truncation_certified: bool
    reason: str
    branches: tuple[ConductorBranchJetDatum, ...]

    @property
    def requires_band_recovery(self) -> bool:
        return self.status in {
            ConductorJetStatus.UNCERTIFIED,
            ConductorJetStatus.INSUFFICIENT,
        }

    def to_dict(self) -> dict[str, Any]:
        return {
            "status": self.status.value,
            "truncation_certified": self.truncation_certified,
            "requires_band_recovery": self.requires_band_recovery,
            "reason": self.reason,
            "branches": [branch.to_dict() for branch in self.branches],
        }


def audit_conductor_jet_truncation(
    branches: Iterable[ConductorBranchJetDatum] | None,
) -> ConductorJetAudit:
    """Apply the conductor/contact-loss theorem to a finite branch ledger."""

    if branches is None:
        return ConductorJetAudit(
            status=ConductorJetStatus.NOT_DECLARED,
            truncation_certified=False,
            reason="no conductor-jet ledger was declared",
            branches=(),
        )
    branch_tuple = tuple(branches)
    if not branch_tuple:
        return ConductorJetAudit(
            status=ConductorJetStatus.NOT_DECLARED,
            truncation_certified=False,
            reason="the declared conductor-jet ledger has no supported branches",
            branches=(),
        )
    names = [branch.name for branch in branch_tuple]
    if len(set(names)) != len(names):
        raise ValueError("conductor branch names must be distinct")

    uncertified = [
        branch.name
        for branch in branch_tuple
        if not branch.certificates_complete
        or branch.available_jet_order is None
    ]
    if uncertified:
        return ConductorJetAudit(
            status=ConductorJetStatus.UNCERTIFIED,
            truncation_certified=False,
            reason="missing jet or proof data on: " + ", ".join(uncertified),
            branches=branch_tuple,
        )

    insufficient = [
        branch.name
        for branch in branch_tuple
        if branch.margin is not None and branch.margin < 0
    ]
    if insufficient:
        return ConductorJetAudit(
            status=ConductorJetStatus.INSUFFICIENT,
            truncation_certified=False,
            reason=(
                "available jets do not reach the conductor/contact-loss bound on: "
                + ", ".join(insufficient)
            ),
            branches=branch_tuple,
        )

    return ConductorJetAudit(
        status=ConductorJetStatus.PASSES,
        truncation_certified=True,
        reason="every branch reaches its conductor/contact-loss bound",
        branches=branch_tuple,
    )


def conductor_branch_jet_datum_from_dict(
    data: dict[str, object],
) -> ConductorBranchJetDatum:
    """Parse one JSON branch without truthiness-based type coercions."""

    integer_fields = (
        "conductor_exponent",
        "differential_order",
        "pole_order",
        "additional_contact_loss",
    )
    parsed_integers: dict[str, int] = {}
    for field_name in integer_fields:
        value = data.get(field_name, 0)
        if isinstance(value, bool) or not isinstance(value, int):
            raise ValueError(f"{field_name} must be an integer")
        parsed_integers[field_name] = value
    available = data.get("available_jet_order")
    if available is not None and (
        isinstance(available, bool) or not isinstance(available, int)
    ):
        raise ValueError("available_jet_order must be null or an integer")
    name = data.get("name")
    if not isinstance(name, str):
        raise ValueError("conductor branch name must be a string")
    certificate_fields = (
        "conductor_certificate",
        "expression_certificate",
        "valuation_certificate",
    )
    certificates = {
        field_name: data.get(field_name, "") for field_name in certificate_fields
    }
    if not all(isinstance(value, str) for value in certificates.values()):
        raise ValueError("conductor-jet certificate references must be strings")
    return ConductorBranchJetDatum(
        name=name,
        available_jet_order=available,
        conductor_certificate=certificates["conductor_certificate"],
        expression_certificate=certificates["expression_certificate"],
        valuation_certificate=certificates["valuation_certificate"],
        **parsed_integers,
    )


@dataclass(frozen=True)
class NormalJetInputDatum:
    """Available normal-parameter jet order for one named input series."""

    name: str
    available_jet_order: int | None
    valuation_certificate: str

    def __post_init__(self) -> None:
        if not isinstance(self.name, str) or not self.name.strip():
            raise ValueError("a normal-jet input needs a nonempty name")
        if self.available_jet_order is not None and (
            isinstance(self.available_jet_order, bool)
            or not isinstance(self.available_jet_order, int)
            or self.available_jet_order < 0
        ):
            raise ValueError(
                "available_jet_order must be null or a nonnegative integer"
            )
        if not isinstance(self.valuation_certificate, str):
            raise ValueError("valuation_certificate must be a string")


@dataclass(frozen=True)
class ContactExpression:
    """A small expression tree with mechanically propagated contact loss."""

    operation: str
    operands: tuple["ContactExpression", ...] = ()
    input_name: str | None = None
    loss_order: int = 0

    _NULLARY = frozenset({"constant", "input"})
    _NARY = frozenset({"add", "multiply"})
    _UNARY_ZERO_LOSS = frozenset({"unit_change", "invert_unit"})
    _UNARY_LOSS = frozenset({"derivative", "pole", "certified_loss"})

    def __post_init__(self) -> None:
        allowed = (
            self._NULLARY
            | self._NARY
            | self._UNARY_ZERO_LOSS
            | self._UNARY_LOSS
        )
        if self.operation not in allowed:
            raise ValueError(f"unknown contact-expression operation {self.operation!r}")
        if (
            isinstance(self.loss_order, bool)
            or not isinstance(self.loss_order, int)
            or self.loss_order < 0
        ):
            raise ValueError("contact-expression loss_order must be nonnegative")
        if self.operation == "input":
            if (
                not isinstance(self.input_name, str)
                or not self.input_name.strip()
                or self.operands
                or self.loss_order != 0
            ):
                raise ValueError("an input node needs only a nonempty input_name")
        elif self.operation == "constant":
            if self.input_name is not None or self.operands or self.loss_order != 0:
                raise ValueError("a constant node has no input, operands, or loss")
        elif self.operation in self._NARY:
            if self.input_name is not None or not self.operands or self.loss_order != 0:
                raise ValueError(
                    f"{self.operation} needs operands and has no direct loss"
                )
        elif self.operation in self._UNARY_ZERO_LOSS:
            if (
                self.input_name is not None
                or len(self.operands) != 1
                or self.loss_order != 0
            ):
                raise ValueError(f"{self.operation} is a zero-loss unary operation")
        elif self.input_name is not None or len(self.operands) != 1:
            raise ValueError(f"{self.operation} is a unary loss operation")

    @classmethod
    def input(cls, name: str) -> "ContactExpression":
        return cls("input", input_name=name)

    @classmethod
    def constant(cls) -> "ContactExpression":
        return cls("constant")

    @classmethod
    def combine(
        cls, operation: str, *operands: "ContactExpression"
    ) -> "ContactExpression":
        return cls(operation, operands=tuple(operands))

    @classmethod
    def shift(
        cls, operation: str, operand: "ContactExpression", order: int
    ) -> "ContactExpression":
        return cls(operation, operands=(operand,), loss_order=order)

    @property
    def sensitivities(self) -> dict[str, int]:
        """Maximum certified contact loss along every input dependency path."""

        if self.operation == "constant":
            return {}
        if self.operation == "input":
            assert self.input_name is not None
            return {self.input_name: 0}
        merged: dict[str, int] = {}
        for operand in self.operands:
            for input_name, loss in operand.sensitivities.items():
                merged[input_name] = max(merged.get(input_name, -1), loss)
        if self.operation in self._UNARY_LOSS:
            merged = {
                input_name: loss + self.loss_order
                for input_name, loss in merged.items()
            }
        return merged

    def to_dict(self) -> dict[str, Any]:
        result: dict[str, Any] = {"operation": self.operation}
        if self.input_name is not None:
            result["input_name"] = self.input_name
        if self.operands:
            result["operands"] = [operand.to_dict() for operand in self.operands]
        if self.loss_order:
            result["loss_order"] = self.loss_order
        return result


@dataclass(frozen=True)
class BoundaryOutputExpressionDatum:
    """One matching-matrix entry or distinguished residue component."""

    name: str
    expression: ContactExpression
    expression_certificate: str

    def __post_init__(self) -> None:
        if not isinstance(self.name, str) or not self.name.strip():
            raise ValueError("a boundary output needs a nonempty name")
        if not isinstance(self.expression, ContactExpression):
            raise ValueError("boundary output expression has the wrong type")
        if not isinstance(self.expression_certificate, str):
            raise ValueError("expression_certificate must be a string")


@dataclass(frozen=True)
class ConductorBranchSensitivityDatum:
    """Dependency-sensitive conductor ledger on one normalization branch."""

    name: str
    conductor_exponent: int
    inputs: tuple[NormalJetInputDatum, ...]
    outputs: tuple[BoundaryOutputExpressionDatum, ...]
    conductor_certificate: str
    dependency_completeness_certificate: str

    def __post_init__(self) -> None:
        if not isinstance(self.name, str) or not self.name.strip():
            raise ValueError("a conductor branch needs a nonempty name")
        if (
            isinstance(self.conductor_exponent, bool)
            or not isinstance(self.conductor_exponent, int)
            or self.conductor_exponent <= 0
        ):
            raise ValueError("conductor_exponent must be positive")
        if not self.outputs:
            raise ValueError("a sensitivity branch needs at least one output")
        input_names = [item.name for item in self.inputs]
        output_names = [item.name for item in self.outputs]
        if len(set(input_names)) != len(input_names):
            raise ValueError("normal-jet input names must be distinct on a branch")
        if len(set(output_names)) != len(output_names):
            raise ValueError("boundary output names must be distinct on a branch")
        declared = set(input_names)
        referenced = {
            input_name
            for output in self.outputs
            for input_name in output.expression.sensitivities
        }
        unknown = referenced - declared
        if unknown:
            raise ValueError(
                "boundary expressions reference undeclared inputs: "
                + ", ".join(sorted(unknown))
            )
        if not isinstance(self.conductor_certificate, str) or not isinstance(
            self.dependency_completeness_certificate, str
        ):
            raise ValueError("branch certificate references must be strings")

    @property
    def used_input_names(self) -> frozenset[str]:
        return frozenset(
            input_name
            for output in self.outputs
            for input_name in output.expression.sensitivities
        )

    @property
    def unused_input_names(self) -> tuple[str, ...]:
        return tuple(
            item.name for item in self.inputs if item.name not in self.used_input_names
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "conductor_exponent": self.conductor_exponent,
            "inputs": [asdict(item) for item in self.inputs],
            "outputs": [
                {
                    "name": output.name,
                    "expression": output.expression.to_dict(),
                    "sensitivities": output.expression.sensitivities,
                    "expression_certificate": output.expression_certificate,
                }
                for output in self.outputs
            ],
            "conductor_certificate": self.conductor_certificate,
            "dependency_completeness_certificate": (
                self.dependency_completeness_certificate
            ),
        }


@dataclass(frozen=True)
class DependencyRequirement:
    branch: str
    output: str
    input_name: str
    contact_loss: int
    conductor_exponent: int
    required_jet_order: int
    available_jet_order: int | None
    margin: int | None
    deficit: int | None
    certified: bool

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class ConductorSensitivityAudit:
    """Exact branch/output/input truncation decision and deficit report."""

    status: ConductorJetStatus
    truncation_certified: bool
    reason: str
    branches: tuple[ConductorBranchSensitivityDatum, ...]
    requirements: tuple[DependencyRequirement, ...]
    unused_inputs: tuple[tuple[str, tuple[str, ...]], ...]

    @property
    def requires_band_recovery(self) -> bool:
        return self.status in {
            ConductorJetStatus.UNCERTIFIED,
            ConductorJetStatus.INSUFFICIENT,
        }

    @property
    def maximum_deficit(self) -> int | None:
        deficits = [
            item.deficit for item in self.requirements if item.deficit is not None
        ]
        return max(deficits, default=None)

    def to_dict(self) -> dict[str, Any]:
        return {
            "mode": "dependency_sensitive",
            "status": self.status.value,
            "truncation_certified": self.truncation_certified,
            "requires_band_recovery": self.requires_band_recovery,
            "maximum_deficit": self.maximum_deficit,
            "reason": self.reason,
            "branches": [branch.to_dict() for branch in self.branches],
            "requirements": [item.to_dict() for item in self.requirements],
            "unused_inputs": {
                branch: list(names) for branch, names in self.unused_inputs
            },
        }


def audit_conductor_sensitivity_ledger(
    branches: Iterable[ConductorBranchSensitivityDatum] | None,
) -> ConductorSensitivityAudit:
    """Compile the full branch/output/input contact-loss inequalities."""

    if branches is None:
        return ConductorSensitivityAudit(
            status=ConductorJetStatus.NOT_DECLARED,
            truncation_certified=False,
            reason="no dependency-sensitive conductor ledger was declared",
            branches=(),
            requirements=(),
            unused_inputs=(),
        )
    branch_tuple = tuple(branches)
    if not branch_tuple:
        return ConductorSensitivityAudit(
            status=ConductorJetStatus.NOT_DECLARED,
            truncation_certified=False,
            reason="the dependency-sensitive conductor ledger has no branches",
            branches=(),
            requirements=(),
            unused_inputs=(),
        )
    branch_names = [branch.name for branch in branch_tuple]
    if len(set(branch_names)) != len(branch_names):
        raise ValueError("sensitivity-ledger branch names must be distinct")

    requirements: list[DependencyRequirement] = []
    incomplete_labels: list[str] = []
    for branch in branch_tuple:
        branch_complete = bool(
            branch.conductor_certificate.strip()
            and branch.dependency_completeness_certificate.strip()
        )
        if not branch_complete:
            incomplete_labels.append(branch.name)
        inputs = {item.name: item for item in branch.inputs}
        for output in branch.outputs:
            output_complete = bool(output.expression_certificate.strip())
            if not output_complete:
                incomplete_labels.append(f"{branch.name}/{output.name}")
            for input_name, contact_loss in sorted(
                output.expression.sensitivities.items()
            ):
                input_datum = inputs[input_name]
                required = branch.conductor_exponent + contact_loss
                available = input_datum.available_jet_order
                margin = None if available is None else available - required
                certified = bool(
                    branch_complete
                    and output_complete
                    and input_datum.valuation_certificate.strip()
                    and available is not None
                )
                if not certified:
                    incomplete_labels.append(
                        f"{branch.name}/{output.name}<-{input_name}"
                    )
                requirements.append(
                    DependencyRequirement(
                        branch=branch.name,
                        output=output.name,
                        input_name=input_name,
                        contact_loss=contact_loss,
                        conductor_exponent=branch.conductor_exponent,
                        required_jet_order=required,
                        available_jet_order=available,
                        margin=margin,
                        deficit=(
                            None if margin is None else max(0, -margin)
                        ),
                        certified=certified,
                    )
                )

    unused = tuple(
        (branch.name, branch.unused_input_names) for branch in branch_tuple
    )
    if incomplete_labels:
        return ConductorSensitivityAudit(
            status=ConductorJetStatus.UNCERTIFIED,
            truncation_certified=False,
            reason=(
                "missing jet or proof data on: "
                + ", ".join(dict.fromkeys(incomplete_labels))
            ),
            branches=branch_tuple,
            requirements=tuple(requirements),
            unused_inputs=unused,
        )
    short = [
        f"{item.branch}/{item.output}<-{item.input_name}"
        for item in requirements
        if item.margin is not None and item.margin < 0
    ]
    if short:
        return ConductorSensitivityAudit(
            status=ConductorJetStatus.INSUFFICIENT,
            truncation_certified=False,
            reason="available jets miss the detailed bound on: " + ", ".join(short),
            branches=branch_tuple,
            requirements=tuple(requirements),
            unused_inputs=unused,
        )
    return ConductorSensitivityAudit(
        status=ConductorJetStatus.PASSES,
        truncation_certified=True,
        reason="every used input reaches every output-specific conductor bound",
        branches=branch_tuple,
        requirements=tuple(requirements),
        unused_inputs=unused,
    )


@dataclass(frozen=True)
class NormalValuationRule:
    """Monomial valuation and trivialization shift on one normal branch."""

    coordinate_names: tuple[str, ...]
    coordinate_orders: tuple[int, ...]
    normalization_shift: int
    certificate: str

    def __post_init__(self) -> None:
        if not self.coordinate_names or len(self.coordinate_names) != len(
            self.coordinate_orders
        ):
            raise ValueError("valuation coordinates and orders must have equal length")
        if len(set(self.coordinate_names)) != len(self.coordinate_names) or not all(
            isinstance(name, str) and name.strip() for name in self.coordinate_names
        ):
            raise ValueError("valuation coordinate names must be nonempty and distinct")
        for value in (*self.coordinate_orders, self.normalization_shift):
            if isinstance(value, bool) or not isinstance(value, int):
                raise ValueError("normal valuation orders and shift must be integers")
        if not isinstance(self.certificate, str):
            raise ValueError("normal valuation certificate must be a string")

    def order(self, exponent_vector: tuple[int, ...]) -> int:
        if len(exponent_vector) != len(self.coordinate_orders):
            raise ValueError("Newton exponent vector has the wrong dimension")
        if any(
            isinstance(value, bool) or not isinstance(value, int)
            for value in exponent_vector
        ):
            raise ValueError("Newton exponents must be integers")
        return self.normalization_shift + sum(
            exponent * weight
            for exponent, weight in zip(exponent_vector, self.coordinate_orders)
        )


@dataclass(frozen=True)
class OmittedNewtonSupportDatum:
    """Certified valuation frontier of the unknown Newton support.

    The completeness certificate must prove that every omitted monomial has
    normal order at least the minimum attained on the displayed vectors.  A
    coordinatewise antichain alone is enough only when the valuation is
    monotone on the omitted-support cone.
    """

    input_name: str
    first_omitted_exponents: tuple[tuple[int, ...], ...]
    completeness_certificate: str

    def __post_init__(self) -> None:
        if not isinstance(self.input_name, str) or not self.input_name.strip():
            raise ValueError("omitted Newton support needs an input name")
        if not self.first_omitted_exponents:
            raise ValueError("list at least one first omitted Newton monomial")
        dimensions = {len(vector) for vector in self.first_omitted_exponents}
        if len(dimensions) != 1:
            raise ValueError("omitted Newton exponent vectors need one dimension")
        if not isinstance(self.completeness_certificate, str):
            raise ValueError("support completeness certificate must be a string")


@dataclass(frozen=True)
class NewtonSupportJetCompilation:
    input: NormalJetInputDatum
    coordinate_names: tuple[str, ...]
    coordinate_orders: tuple[int, ...]
    normalization_shift: int
    first_omitted_exponents: tuple[tuple[int, ...], ...]
    first_omitted_normal_orders: tuple[int, ...]

    def to_dict(self) -> dict[str, Any]:
        return {
            "input": asdict(self.input),
            "coordinate_names": list(self.coordinate_names),
            "coordinate_orders": list(self.coordinate_orders),
            "normalization_shift": self.normalization_shift,
            "first_omitted_exponents": [
                list(vector) for vector in self.first_omitted_exponents
            ],
            "first_omitted_normal_orders": list(
                self.first_omitted_normal_orders
            ),
        }


def compile_newton_support_jet(
    support: OmittedNewtonSupportDatum,
    valuation: NormalValuationRule,
) -> NewtonSupportJetCompilation:
    """Convert a certified omitted-support valuation frontier to a jet order."""

    orders = tuple(
        valuation.order(vector) for vector in support.first_omitted_exponents
    )
    if min(orders) < 0:
        raise ValueError(
            "an omitted monomial remains polar after the declared trivialization"
        )
    certificate = ""
    if valuation.certificate.strip() and support.completeness_certificate.strip():
        certificate = (
            f"{valuation.certificate}; {support.completeness_certificate}"
        )
    return NewtonSupportJetCompilation(
        input=NormalJetInputDatum(
            name=support.input_name,
            available_jet_order=min(orders),
            valuation_certificate=certificate,
        ),
        coordinate_names=valuation.coordinate_names,
        coordinate_orders=valuation.coordinate_orders,
        normalization_shift=valuation.normalization_shift,
        first_omitted_exponents=support.first_omitted_exponents,
        first_omitted_normal_orders=orders,
    )


def contact_expression_from_dict(data: dict[str, object]) -> ContactExpression:
    operation = data.get("operation")
    if not isinstance(operation, str):
        raise ValueError("contact-expression operation must be a string")
    raw_operands = data.get("operands", [])
    if not isinstance(raw_operands, list) or not all(
        isinstance(item, dict) for item in raw_operands
    ):
        raise ValueError("contact-expression operands must be JSON objects")
    input_name = data.get("input_name")
    if input_name is not None and not isinstance(input_name, str):
        raise ValueError("contact-expression input_name must be null or a string")
    loss_order = data.get("loss_order", 0)
    if isinstance(loss_order, bool) or not isinstance(loss_order, int):
        raise ValueError("contact-expression loss_order must be an integer")
    return ContactExpression(
        operation=operation,
        operands=tuple(contact_expression_from_dict(item) for item in raw_operands),
        input_name=input_name,
        loss_order=loss_order,
    )


def conductor_branch_sensitivity_from_dict(
    data: dict[str, object],
) -> ConductorBranchSensitivityDatum:
    """Parse one dependency-sensitive branch ledger from JSON."""

    name = data.get("name")
    conductor = data.get("conductor_exponent")
    if not isinstance(name, str):
        raise ValueError("sensitivity branch name must be a string")
    if isinstance(conductor, bool) or not isinstance(conductor, int):
        raise ValueError("sensitivity conductor_exponent must be an integer")
    raw_inputs = data.get("inputs", [])
    raw_outputs = data.get("outputs", [])
    if not isinstance(raw_inputs, list) or not all(
        isinstance(item, dict) for item in raw_inputs
    ):
        raise ValueError("sensitivity inputs must be JSON objects")
    if not isinstance(raw_outputs, list) or not all(
        isinstance(item, dict) for item in raw_outputs
    ):
        raise ValueError("sensitivity outputs must be JSON objects")
    raw_valuation = data.get("normal_valuation")
    valuation_rule: NormalValuationRule | None = None
    if raw_valuation is not None:
        if not isinstance(raw_valuation, dict):
            raise ValueError("normal_valuation must be one JSON object")
        coordinate_names = raw_valuation.get("coordinate_names")
        coordinate_orders = raw_valuation.get("coordinate_orders")
        normalization_shift = raw_valuation.get("normalization_shift", 0)
        valuation_certificate = raw_valuation.get("certificate", "")
        if (
            not isinstance(coordinate_names, list)
            or not all(isinstance(item, str) for item in coordinate_names)
            or not isinstance(coordinate_orders, list)
            or not all(
                not isinstance(item, bool) and isinstance(item, int)
                for item in coordinate_orders
            )
            or isinstance(normalization_shift, bool)
            or not isinstance(normalization_shift, int)
            or not isinstance(valuation_certificate, str)
        ):
            raise ValueError("normal_valuation has invalid coordinates or orders")
        valuation_rule = NormalValuationRule(
            coordinate_names=tuple(coordinate_names),
            coordinate_orders=tuple(coordinate_orders),
            normalization_shift=normalization_shift,
            certificate=valuation_certificate,
        )
    inputs = []
    for item in raw_inputs:
        input_name = item.get("name")
        if not isinstance(input_name, str):
            raise ValueError("sensitivity input names must be strings")
        raw_omitted = item.get("first_omitted_exponents")
        if raw_omitted is not None:
            if "available_jet_order" in item:
                raise ValueError(
                    "derive a sensitivity input from omitted support or declare "
                    "available_jet_order, not both"
                )
            if valuation_rule is None:
                raise ValueError(
                    "first_omitted_exponents require a branch normal_valuation"
                )
            completeness = item.get("support_completeness_certificate", "")
            if (
                not isinstance(raw_omitted, list)
                or not raw_omitted
                or not all(
                    isinstance(vector, list)
                    and all(
                        not isinstance(value, bool) and isinstance(value, int)
                        for value in vector
                    )
                    for vector in raw_omitted
                )
                or not isinstance(completeness, str)
            ):
                raise ValueError("first omitted Newton support is malformed")
            compiled_input = compile_newton_support_jet(
                OmittedNewtonSupportDatum(
                    input_name=input_name,
                    first_omitted_exponents=tuple(
                        tuple(vector) for vector in raw_omitted
                    ),
                    completeness_certificate=completeness,
                ),
                valuation_rule,
            )
            inputs.append(compiled_input.input)
            continue
        available = item.get("available_jet_order")
        certificate = item.get("valuation_certificate", "")
        if not isinstance(certificate, str):
            raise ValueError("sensitivity input valuation certificate is a string")
        inputs.append(NormalJetInputDatum(input_name, available, certificate))
    outputs = []
    for item in raw_outputs:
        output_name = item.get("name")
        expression = item.get("expression")
        certificate = item.get("expression_certificate", "")
        if (
            not isinstance(output_name, str)
            or not isinstance(expression, dict)
            or not isinstance(certificate, str)
        ):
            raise ValueError("sensitivity output needs name, expression, certificate")
        outputs.append(
            BoundaryOutputExpressionDatum(
                name=output_name,
                expression=contact_expression_from_dict(expression),
                expression_certificate=certificate,
            )
        )
    conductor_certificate = data.get("conductor_certificate", "")
    dependency_certificate = data.get(
        "dependency_completeness_certificate", ""
    )
    if not isinstance(conductor_certificate, str) or not isinstance(
        dependency_certificate, str
    ):
        raise ValueError("sensitivity branch certificates must be strings")
    return ConductorBranchSensitivityDatum(
        name=name,
        conductor_exponent=conductor,
        inputs=tuple(inputs),
        outputs=tuple(outputs),
        conductor_certificate=conductor_certificate,
        dependency_completeness_certificate=dependency_certificate,
    )


__all__ = [
    "BoundaryOutputExpressionDatum",
    "ConductorBranchJetDatum",
    "ConductorBranchSensitivityDatum",
    "ConductorJetAudit",
    "ConductorJetStatus",
    "ConductorSensitivityAudit",
    "ContactExpression",
    "DependencyRequirement",
    "NewtonSupportJetCompilation",
    "NormalJetInputDatum",
    "NormalValuationRule",
    "OmittedNewtonSupportDatum",
    "audit_conductor_sensitivity_ledger",
    "audit_conductor_jet_truncation",
    "compile_newton_support_jet",
    "conductor_branch_sensitivity_from_dict",
    "conductor_branch_jet_datum_from_dict",
    "contact_expression_from_dict",
]

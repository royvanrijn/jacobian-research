#!/usr/bin/env python3
"""Exact regressions for the conductor/contact-loss truncation theorem."""

from __future__ import annotations

import sys
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from conductor_jet_truncation import (  # noqa: E402
    BoundaryOutputExpressionDatum,
    ConductorBranchJetDatum,
    ConductorBranchSensitivityDatum,
    ContactExpression,
    NormalJetInputDatum,
    NormalValuationRule,
    OmittedNewtonSupportDatum,
    audit_conductor_sensitivity_ledger,
    audit_conductor_jet_truncation,
    compile_newton_support_jet,
    conductor_branch_sensitivity_from_dict,
    conductor_branch_jet_datum_from_dict,
)


t = sp.symbols("t")


def order(expression: sp.Expr) -> int:
    """Return the exact t-adic order of a nonzero polynomial."""

    polynomial = sp.Poly(sp.expand(expression), t)
    if polynomial.is_zero:
        return sp.oo
    return min(exponent[0] for exponent, _ in polynomial.terms())


def truncated_coefficients(expression: sp.Expr, length: int) -> tuple[sp.Expr, ...]:
    polynomial = sp.Poly(sp.expand(expression), t)
    return tuple(polynomial.nth(index) for index in range(length))


def contact_expression(
    expression: sp.Expr, differential_order: int, pole_order: int
) -> sp.Expr:
    return sp.cancel(sp.diff(expression, t, differential_order) / t**pole_order)


def numerical_semigroup(generators: tuple[int, ...]) -> tuple[int, tuple[int, ...]]:
    """Compute the conductor and gaps of a gcd-one numerical semigroup."""

    if not generators or any(value <= 0 for value in generators):
        raise ValueError("semigroup generators must be positive")
    if sp.gcd_list(generators) != 1:
        raise ValueError("semigroup generators must have gcd one")
    largest = max(generators)
    reachable = {0}
    bound = largest
    while True:
        for value in range(bound + 1):
            if any(value - generator in reachable for generator in generators):
                reachable.add(value)
        for candidate in range(bound + 1):
            if all(candidate + offset in reachable for offset in range(largest)):
                gaps = tuple(value for value in range(candidate) if value not in reachable)
                return candidate, gaps
        bound *= 2


def assert_threshold_is_sufficient_and_sharp(
    conductor: int, differential_order: int, pole_order: int
) -> None:
    """Check sufficiency and the one-order-short counterexample over Q."""

    contact_loss = differential_order + pole_order
    required = conductor + contact_loss
    high_perturbation = t**required * (1 + 2 * t + 3 * t**2)
    high_output = contact_expression(
        high_perturbation, differential_order, pole_order
    )
    assert order(high_output) >= conductor

    low_perturbation = t ** (required - 1)
    low_output = contact_expression(
        low_perturbation, differential_order, pole_order
    )
    assert order(low_output) == conductor - 1


# The two basic conductor quotients.
node_class = lambda expression: sp.expand(expression.subs(t, 1) - expression.subs(t, 0))
cusp_class = lambda expression: sp.Poly(sp.expand(expression), t).nth(1)

node_conductor = t * (t - 1)
assert node_class(node_conductor * (3 + 5 * t)) == 0
assert node_class(1 + 7 * t) == 7
assert cusp_class(t**2 * (3 + 5 * t)) == 0
assert cusp_class(1 + 7 * t) == 7
print("PASS: node and cusp conductor ideals annihilate their normalization quotients")


# Contact loss: c+d+ell is sufficient and is sharp in general.
for conductor in (1, 2, 6):
    for differential_order in range(3):
        for pole_order in range(3):
            assert_threshold_is_sufficient_and_sharp(
                conductor, differential_order, pole_order
            )
print("PASS: the c+d+pole jet bound is sufficient and one-order sharp over Q")


# Arbitrary unibranch monomial curves: gaps give a basis for B/A, while all
# exponents at least the numerical-semigroup conductor disappear.
semigroup_expectations = {
    (2, 3): (2, (1,)),
    (3, 4): (6, (1, 2, 5)),
    (4, 6, 9): (12, (1, 2, 3, 5, 7, 11)),
}
for generators, expected in semigroup_expectations.items():
    conductor, gaps = numerical_semigroup(generators)
    assert (conductor, gaps) == expected
    coefficients = sp.symbols(f"a0:{conductor + 3}")
    expression = sum(
        coefficient * t**exponent
        for exponent, coefficient in enumerate(coefficients)
    )
    quotient_vector = tuple(
        truncated_coefficients(expression, conductor)[gap] for gap in gaps
    )
    perturbed = expression + t**conductor * (1 + t + t**2)
    perturbed_vector = tuple(
        truncated_coefficients(perturbed, conductor)[gap] for gap in gaps
    )
    assert quotient_vector == perturbed_vector
print("PASS: numerical-semigroup conductor jets recover every monomial-curve quotient")


def branch(
    name: str,
    conductor: int,
    derivative: int,
    pole: int,
    available: int | None,
    *,
    additional: int = 0,
    certificates: bool = True,
) -> ConductorBranchJetDatum:
    certificate = "exact fixture" if certificates else ""
    return ConductorBranchJetDatum(
        name=name,
        conductor_exponent=conductor,
        differential_order=derivative,
        pole_order=pole,
        additional_contact_loss=additional,
        available_jet_order=available,
        conductor_certificate=certificate,
        expression_certificate=certificate,
        valuation_certificate=certificate,
    )


passing = audit_conductor_jet_truncation(
    (
        branch("node-left", 1, 1, 2, 4),
        branch("node-right", 1, 0, 0, 1),
        branch("cusp", 2, 2, 1, 6, additional=1),
    )
)
assert passing.status.value == "passes"
assert passing.truncation_certified
assert [item.margin for item in passing.branches] == [0, 0, 0]

insufficient = audit_conductor_jet_truncation(
    (branch("cusp", 2, 2, 1, 4),)
)
assert insufficient.status.value == "insufficient"
assert insufficient.requires_band_recovery
assert insufficient.branches[0].margin == -1

uncertified = audit_conductor_jet_truncation(
    (branch("unknown", 2, 1, 0, None, certificates=False),)
)
assert uncertified.status.value == "uncertified"
assert uncertified.requires_band_recovery

absent = audit_conductor_jet_truncation(None)
assert absent.status.value == "not_declared"
assert not absent.requires_band_recovery
print("PASS: the proof-bearing audit separates pass, insufficiency, and missing data")


# The JSON parser rejects booleans masquerading as integers and preserves the
# exact threshold convention.
parsed = conductor_branch_jet_datum_from_dict(
    {
        "name": "semigroup-branch",
        "conductor_exponent": 6,
        "differential_order": 2,
        "pole_order": 3,
        "additional_contact_loss": 1,
        "available_jet_order": 12,
        "conductor_certificate": "semigroup gaps",
        "expression_certificate": "expression tree",
        "valuation_certificate": "normal valuation ledger",
    }
)
assert parsed.contact_loss == 6
assert parsed.required_jet_order == 12
assert parsed.margin == 0
try:
    conductor_branch_jet_datum_from_dict(
        {
            "name": "bad",
            "conductor_exponent": True,
        }
    )
except ValueError:
    pass
else:
    raise AssertionError("boolean conductor exponent was accepted")
print("PASS: JSON conductor-jet records preserve strict integer semantics")


# A cokernel class only sees the quotient coordinates.  This finite matrix
# fixture represents two matching columns and a distinguished determinant
# residue in the cusp quotient with gap basis [t].
u0, u1, v0, v1 = sp.symbols("u0 u1 v0 v1")
source = u0 + u1 * t
target = v0 + v1 * t
matching_columns = sp.Matrix([[cusp_class(source), cusp_class(target)]])
residue = sp.Matrix([cusp_class(source - target)])
source_perturbed = source + t**2 * (1 + t)
target_perturbed = target + t**3
perturbed_columns = sp.Matrix(
    [[cusp_class(source_perturbed), cusp_class(target_perturbed)]]
)
perturbed_residue = sp.Matrix(
    [cusp_class(source_perturbed - target_perturbed)]
)
assert matching_columns == perturbed_columns
assert residue == perturbed_residue
print("PASS: matching-map cokernels and distinguished classes use only conductor jets")


# Dependency-sensitive compilation avoids coupling the shortest input jet to
# the largest loss of an unrelated output path.  P is differentiated once;
# Q is differentiated twice and then divided by t^2.  R is certified unused
# by the complete expression graph and therefore needs no jet certificate.
p_input = ContactExpression.input("P")
q_input = ContactExpression.input("Q")
rho_expression = ContactExpression.combine(
    "add",
    ContactExpression.shift("derivative", p_input, 1),
    ContactExpression.shift(
        "pole",
        ContactExpression.shift("derivative", q_input, 2),
        2,
    ),
)
detailed_branch = ConductorBranchSensitivityDatum(
    name="case1-shaped-cusp",
    conductor_exponent=2,
    inputs=(
        NormalJetInputDatum("P", 3, "P known modulo t^3"),
        NormalJetInputDatum("Q", 6, "Q known modulo t^6"),
        NormalJetInputDatum("R", None, ""),
    ),
    outputs=(
        BoundaryOutputExpressionDatum(
            name="rho",
            expression=rho_expression,
            expression_certificate="complete synthetic P/Q residue tree",
        ),
    ),
    conductor_certificate="cusp conductor (t^2)",
    dependency_completeness_certificate="rho has exactly the displayed P/Q paths",
)
detailed = audit_conductor_sensitivity_ledger((detailed_branch,))
assert detailed.status.value == "passes"
assert detailed.truncation_certified
assert detailed.unused_inputs == (("case1-shaped-cusp", ("R",)),)
assert [
    (item.input_name, item.contact_loss, item.required_jet_order, item.margin)
    for item in detailed.requirements
] == [("P", 1, 3, 0), ("Q", 4, 6, 0)]

# The old scalar compression couples min(3,6) to max(1,4), so it fails even
# though every actual dependency passes.  This is the optimization supplied
# by the matrix ledger.
scalar_compression = audit_conductor_jet_truncation(
    (branch("scalar-compression", 2, 2, 2, 3),)
)
assert scalar_compression.status.value == "insufficient"

detailed_short = audit_conductor_sensitivity_ledger(
    (
        ConductorBranchSensitivityDatum(
            name=detailed_branch.name,
            conductor_exponent=detailed_branch.conductor_exponent,
            inputs=(
                detailed_branch.inputs[0],
                NormalJetInputDatum("Q", 5, "Q known modulo t^5"),
                detailed_branch.inputs[2],
            ),
            outputs=detailed_branch.outputs,
            conductor_certificate=detailed_branch.conductor_certificate,
            dependency_completeness_certificate=(
                detailed_branch.dependency_completeness_certificate
            ),
        ),
    )
)
assert detailed_short.status.value == "insufficient"
assert detailed_short.maximum_deficit == 1
q_requirement = next(
    item for item in detailed_short.requirements if item.input_name == "Q"
)
assert q_requirement.deficit == 1
print("PASS: dependency sensitivities avoid false coupling and report exact deficits")


# A toroidal valuation turns a certified valuation frontier of omitted Newton
# support into the available normal-jet order.  No coefficient nonvanishing is
# assumed: the minimum possible omitted order is the safe bound.
valuation_rule = NormalValuationRule(
    coordinate_names=("X", "z"),
    coordinate_orders=(1, 4),
    normalization_shift=0,
    certificate="ord_t(X)=1 and ord_t(z)=4 on the selected branch",
)
p_support = compile_newton_support_jet(
    OmittedNewtonSupportDatum(
        input_name="P",
        first_omitted_exponents=((2, 1), (6, 0)),
        completeness_certificate="all omitted P monomials have weight at least six",
    ),
    valuation_rule,
)
q_support = compile_newton_support_jet(
    OmittedNewtonSupportDatum(
        input_name="Q",
        first_omitted_exponents=((1, 2), (9, 0)),
        completeness_certificate="all omitted Q monomials have weight at least nine",
    ),
    valuation_rule,
)
assert p_support.input.available_jet_order == 6
assert p_support.first_omitted_normal_orders == (6, 6)
assert q_support.input.available_jet_order == 9
assert q_support.first_omitted_normal_orders == (9, 9)

uncertified_frontier = compile_newton_support_jet(
    OmittedNewtonSupportDatum(
        input_name="P",
        first_omitted_exponents=((2, 1),),
        completeness_certificate="",
    ),
    valuation_rule,
)
frontier_audit = audit_conductor_sensitivity_ledger(
    (
        ConductorBranchSensitivityDatum(
            name="uncertified-frontier",
            conductor_exponent=1,
            inputs=(uncertified_frontier.input,),
            outputs=(
                BoundaryOutputExpressionDatum(
                    "rho",
                    ContactExpression.input("P"),
                    "rho is the P conductor class",
                ),
            ),
            conductor_certificate="node conductor",
            dependency_completeness_certificate="rho is the only output",
        ),
    )
)
assert frontier_audit.status.value == "uncertified"
print("PASS: certified omitted-support frontiers compile to exact normal jets")


parsed_sensitivity = conductor_branch_sensitivity_from_dict(
    {
        "name": "parsed-node",
        "conductor_exponent": 1,
        "conductor_certificate": "node conductor",
        "dependency_completeness_certificate": "complete expression tree",
        "inputs": [
            {
                "name": "P",
                "available_jet_order": 2,
                "valuation_certificate": "P known modulo t^2",
            }
        ],
        "outputs": [
            {
                "name": "phi_1",
                "expression_certificate": "first derivative",
                "expression": {
                    "operation": "derivative",
                    "loss_order": 1,
                    "operands": [
                        {"operation": "input", "input_name": "P"}
                    ],
                },
            }
        ],
    }
)
parsed_sensitivity_audit = audit_conductor_sensitivity_ledger(
    (parsed_sensitivity,)
)
assert parsed_sensitivity_audit.status.value == "passes"
assert parsed_sensitivity_audit.requirements[0].required_jet_order == 2
print("PASS: dependency-sensitive JSON expressions compile without scalar coercion")

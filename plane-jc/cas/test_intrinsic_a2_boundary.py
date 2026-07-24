#!/usr/bin/env python3
"""Regression tests for the intrinsic affine-plane boundary gates."""

import sympy as sp

from boundary_lattice_prefilter import (
    boundary_intersection_matrix,
    standard_completion,
)
from intrinsic_a2_boundary import (
    IntrinsicA2Boundary,
    KellerDicriticalDatum,
    audit_a2_boundary,
    audit_keller_pole_vector,
    audit_keller_residual_different,
    contract_keller_vertical_boundary,
    infer_finite_model_dicritical_projection_budget,
    infer_keller_dicritical_budget,
    symmetric_inertia,
)


def package(configuration, initial_form):
    return IntrinsicA2Boundary(
        names=configuration.names,
        intersection_matrix=boundary_intersection_matrix(
            configuration, initial_form
        ),
        genera=(0,) * len(configuration.names),
    )


assert symmetric_inertia(sp.Matrix([[0, 1], [1, 0]])) == (1, 1, 0)
assert symmetric_inertia(sp.diag(1, -1, -2, 0)) == (1, 2, 1)


# Every honest boundary blowup of the standard completion satisfies not only
# SNF/unimodularity but also the adjunction and rational-surface identities.
base, initial_form, _ = standard_completion("A2")
configurations = [base]
current = base
for center in (("L",), ("E1",), ("E1", "E2"), ("E3",)):
    current = current.blow_up(center)
    configurations.append(current)

for configuration in configurations:
    audit = audit_a2_boundary(package(configuration, initial_form))
    assert audit.passes, audit.failures
    assert audit.canonical_square == 10 - len(configuration.names)


# Unimodularity and Hodge signature alone do not characterize an A2
# boundary.  This rational star has det=-1 and inertia (1,3), but adjunction
# gives K^2=-2 instead of 6.
fake_star = IntrinsicA2Boundary(
    names=("C", "A", "B", "D"),
    intersection_matrix=sp.Matrix(
        [
            [-1, 1, 1, 1],
            [1, -5, 0, 0],
            [1, 0, -3, 0],
            [1, 0, 0, -2],
        ]
    ),
    genera=(0, 0, 0, 0),
)
fake_audit = audit_a2_boundary(fake_star)
assert fake_audit.determinant == -1
assert fake_audit.inertia == (1, 3, 0)
assert fake_audit.canonical_square == -2
assert not fake_audit.passes
assert any("K_X^2=-2" in failure for failure in fake_audit.failures)


# The identity P2 -> P2 has pole vector (1), degree one, and zero ordinary
# and logarithmic ramification.  It is proper, so it has no dicritical.
base_audit = audit_a2_boundary(package(base, initial_form))
identity = audit_keller_pole_vector(
    base_audit, (1,), require_nonproper=False
)
assert identity.passes
assert identity.geometric_degree == 1
assert identity.ramification_coefficients == (0,)
assert identity.log_ramification_coefficients == (0,)
assert identity.dicritical_candidates == ()
assert not audit_keller_pole_vector(
    base_audit, (1,), require_nonproper=True
).passes


# Canonical coefficients under free boundary blowups begin
# -3,-2,-1,0.  Hence fewer than three successive one-parent blowups cannot
# contain a dicritical prime of a nonproper Keller resolution.
free = base
free_audits = []
for parent in ("L", "E1", "E2"):
    free = free.blow_up((parent,))
    free_audits.append(audit_a2_boundary(package(free, initial_form)))

assert free_audits[0].canonical_coefficients == (-3, -2)
assert free_audits[1].canonical_coefficients == (-3, -2, -1)
assert not free_audits[1].can_support_nonproper_keller_boundary
assert free_audits[2].canonical_coefficients == (-3, -2, -1, 0)
assert free_audits[2].nonnegative_canonical_components == ("E3",)


# The first numerically possible free-depth package has a pole vector meeting
# every exact canonical, nefness, degree, ramification, and dicritical gate.
# This is a consistency witness, not an existence claim for a Keller map.
first_possible = audit_keller_pole_vector(
    free_audits[2], (3, 2, 1, 0), require_nonproper=True
)
assert first_possible.passes, first_possible.failures
assert first_possible.hyperplane_intersections == (2, 0, 0, 1)
assert first_possible.geometric_degree == 6
assert first_possible.ramification_coefficients == (6, 4, 2, 0)
assert first_possible.log_ramification_coefficients == (4, 3, 2, 1)
assert first_possible.dicritical_candidates == ("E3",)

# The complete graph now feeds the finite-normalization residual-different
# audit.  Here (f^*L).E3=1 forces residue degree one over a line, while the
# adjacent coefficient two forces companion intersection two.  Declaring
# no companion intersection therefore fails the exact identity.
first_budget = infer_keller_dicritical_budget(first_possible, "E3", 1)
assert first_budget.ramification_index == 1
assert first_budget.residue_degree == 1
assert first_budget.available_residual_intersection == 2
assert first_budget.forced_companion_intersection == 2
assert first_budget.feasible

first_contraction = contract_keller_vertical_boundary(first_possible)
assert first_contraction.passes
assert first_contraction.contracted_names == ("E1", "E2")
assert first_contraction.surviving_names == ("L", "E3")
assert first_contraction.contracted_inertia == (0, 2, 0)

# On the finite Stein model the surface different records the singularity
# created where the contracted H-null chain met E3.  For a target line the
# normalization correction is zero and the projection and corrected
# residual budgets agree exactly.
first_projection = infer_finite_model_dicritical_projection_budget(
    first_possible, "E3", 1
)
assert first_projection.budgets_match
assert first_projection.target_normalization_correction == 0
assert (
    first_projection.residual_forced_companion_intersection
    == first_projection.projection_forced_companion_intersection
)

# Scaling the same pole ray makes (f^*L).E3=3.  Interpreting the image as a
# rational cubic leaves a naive budget gap two, exactly the degree of the
# normalization conductor, (c-1)(c-2).  Including that target correction
# restores equality.
cubic_image_package = audit_keller_pole_vector(
    free_audits[2], (9, 6, 3, 0), require_nonproper=True
)
assert cubic_image_package.passes
cubic_projection = infer_finite_model_dicritical_projection_budget(
    cubic_image_package, "E3", 3
)
assert cubic_projection.residue_degree == 1
assert cubic_projection.target_normalization_correction == 2
assert (
    cubic_projection.projection_forced_companion_intersection
    - cubic_projection.residual_forced_companion_intersection
    == 2
)
assert cubic_projection.budgets_match

unpaid_first_possible = audit_keller_residual_different(
    first_possible,
    (KellerDicriticalDatum("E3", 1, 0),),
)
assert not unpaid_first_possible.all_identities_hold
assert unpaid_first_possible.total_available_intersection == 2
assert unpaid_first_possible.total_required_intersection == 0

paid_first_possible = audit_keller_residual_different(
    first_possible,
    (KellerDicriticalDatum("E3", 1, 2),),
)
assert paid_first_possible.all_identities_hold
assert paid_first_possible.component_audits[0].transverse_coefficient_matches

try:
    audit_keller_residual_different(
        first_possible,
        (KellerDicriticalDatum("E2", 1, 0),),
    )
except ValueError as error:
    assert "not dicritical" in str(error)
else:
    raise AssertionError("a nondicritical component entered the residual audit")

try:
    audit_keller_residual_different(
        first_possible,
        (KellerDicriticalDatum("E3", 2, 0),),
    )
except ValueError as error:
    assert "not divisible" in str(error)
else:
    raise AssertionError("a nonintegral residue degree entered the residual audit")


print("PASS: exact inertia handles isotropic affine-plane boundary lattices")
print("PASS: every compiled boundary blowup satisfies adjunction and K^2+rho=10")
print("PASS: the Noether gate rejects a unimodular Hodge-signature fake tree")
print("PASS: pole vectors reconstruct ordinary and logarithmic ramification")
print("PASS: nonproper Keller resolutions require canonical free depth at least three")
print("PASS: intrinsic dicriticals feed the finite residual-different budget")

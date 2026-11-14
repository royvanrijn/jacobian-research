# Compact-group Mathieu classification: source audit

Reviewed 14 November 2026. Christopher D. Long's
[The Mathieu Property for Compact Connected Lie Groups, v1](https://arxiv.org/abs/2609.16178v1)
proves that the kernel of normalized Haar integration on representative
functions is Mathieu--Zhao exactly when the compact connected Lie group is
a torus. This replaces the existence question for simple higher-rank groups
in the archived [SO(3) discussion](../../archive/non-elliptic/extended-geometry/LONG_SU2_AND_BCW_REPRODUCTIONS.md).

The audited dependency chain is Sections 3--7: a homogeneous phase-balanced
Hopf pair; averaging onto invariant radial measures; a fundamental-weight
root doublet and equivariant orthogonal projection; descent using the
center-trivial representation rho tensor conjugate(rho); and pullback from
an adjoint simple quotient. Positivity follows from the projection being
nonzero near the identity. This constructs representative functions on the
whole group, resolving the subgroup-extension problem. The torus direction
uses Duistermaat--van der Kallen. Local review is neither peer review nor a
formal verification of the Lie-theoretic argument.

The independent [checker](verify_moment_transfers.py) verifies 132 exact
sphere/marker identities through m=6 in complex dimensions 2,3,4,6, including
zero tails. Its polynomial expansion agrees with the source's Pascal-row
formula and defining-representation beta factor. This supports transcription;
the all-order and all-group conclusions use the written proof above.

The [componentwise finite-trace deduction](DISCONNECTED_COMPACT_GROUPS.md)
now settles disconnected compact Lie groups as well. Group-specific minimal
representation degree remains outside these classifications. It changes no Gaussian or differential-operator
dimension bound.

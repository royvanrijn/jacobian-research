# Cancellation continuation pointers

This file no longer maintains an independent roadmap.  The sole continuation
queue is generated in [`STATUS.md`](../STATUS.md) from
[`MATH_STATUS.json`](../MATH_STATUS.json).

The primary cancellation continuations are `OP-CR`, `OP-SUSP`, and
`OP-UG3`.
The latter is now the
[minimal-boundary gateway and classification conjecture](MINIMAL_BOUNDARY_CLASSIFICATION.md),
with geometric degree three as its first target.  `MBP1` now formalizes
selected critical boundary, saturated link, boundary monotonicity, ledger
completeness, puncture rank, primitive conormal, noncontraction, and chart
straightenability as predicates on a finite-normalization diagram.  The open
step is to construct that diagram canonically from a numerically minimal
unmarked normalization and prove its predicates.  The subsequent cubic work
is organized as two alternatives: either extract the marked suspension
package, or prove intrinsic
curvilinearity (hence point-flatness) and straighten the
Deligne--Faddeev coefficient orbit.  Completing either cubic route gives the
foundational map and reconstructs the other package, so they should not be
pursued as cumulative independent requirements.
The second unramified nonproperness divisor is no longer an open cubic
certificate: the foundational competitor has one irreducible target
boundary component, so lexicographic boundary minimality forces the same
for every boundary-minimal cubic.  This does not settle an arbitrary
ungraded cubic: no reduction preserving torus-freeness from a larger
boundary ledger to the minimal stratum is known.  The global degree-minimum
target is therefore tracked separately as `OP-UG3`, and still includes the
phantom-divisor unit test.
For a reduced minimal point defect, fiber-minimality is now the concrete
problem of excluding a square-zero length-four collision of the ramified
and affine sheets; the foundational triple-root collision is instead
curvilinear of length three and must remain allowed.
More generally, proving that every intrinsic collision fiber is
curvilinear closes point-flatness outright by the local monogenicity
criterion, including all nonreduced and higher determinantal defects.
Equivalently, it suffices to prove that the relative cotangent module of
every collision fiber is cyclic; this is the unit-`Fitt_1` condition already
visible in the intrinsic scheme package.
In the existing conormal language, the same condition says that the
primitive conormal class generates the full nilradical after every
closed-point specialization.  This closed-point saturation statement is
the common missing lemma in the normalization and branchwise frontends.
It follows from the exact Hartogs package: pure two-dimensional `S_2`
scheme-theoretic ramification support, a rank-one full-support `S_1`
cotangent module, and codimension-one primitive generation.  Equivalently,
the only closed-point obstructions are
`Ext_A^2(T,A)` for the ramification support and `Ext_A^3(Omega_{B/A},A)`
for the primitive-generation cokernel.  Proposition 1.15 replaces them by
the two finite double-saturation quotients
`T^[2]/T` and `Omega_{B/A}/T tau`; their canonical duals are exactly those
two `Ext` modules.  Proposition 1.16 shows that after the first quotient
vanishes, the second is just `H_Z^0(Omega_{B/A})`.  Proposition 1.17 turns
this into the presentation saturation `N:I^infinity=N`, with
`I=Fitt_3(B)`.
The tame local structure theorem removes the simple-normal-crossing branch
locus from these tests: the normalization there is the free sum of a
quadratic Kummer sheet and a trivial sheet.  Thus the remaining saturation
calculations first reduce to closed non-SNC points of the critical
discriminant.  Ordinary-cusp braid monodromy removes the clean cusp locus
as well: its only three-sheet actions are the finite-free `2+1` Kummer and
transitive cubic root covers.  The surviving E target is therefore a
closed worse-than-ordinary-cusp point.  `KBD6` quantifies the reduced row:
the ternary-cubic symbol has degree-six line-section discriminant, so a
squarefree symbol forces branch multiplicity six and a non-squarefree
symbol forces multiplicity at least seven.  The next geometric step must
exclude this multiplicity jump using the marked Keller boundary; the
three-coordinate invariant `mu` does not see it.
The ranked execution plan and the closed phantom-boundary plus remaining
graded gauge certificates are maintained in
[`CUBIC_CLOSURE_ATTACKS.md`](CUBIC_CLOSURE_ATTACKS.md).
Cancellation arithmetic is parked as `OP-ARITH`.  The former roadmap is preserved in
[`archive/legacy-notes`](../archive/legacy-notes/CANCELLATION_RESEARCH_ROADMAP_2026-07-22.md).

The geometric deck group of the cancellation inverse polynomial is treated
separately in [`INVERSE_MONODROMY.md`](INVERSE_MONODROMY.md).  Its generic
`A_N/S_N` theorem and branch-cycle search are independent of the parameter-
polynomial Galois problem in `OP-ARITH`.

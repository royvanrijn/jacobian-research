# Cancellation continuation pointers

This file no longer maintains an independent roadmap.  The sole continuation
queue is generated in [`STATUS.md`](../STATUS.md) from
[`MATH_STATUS.json`](../MATH_STATUS.json).

The primary cancellation continuations are `OP-CR` and `OP-SUSP`.
The latter is now the
[minimal-boundary classification conjecture](MINIMAL_BOUNDARY_CLASSIFICATION.md),
with geometric degree three as its first target.  Its cubic work is organized
as two alternative gateways: either extract the marked suspension package,
or prove intrinsic curvilinearity (hence point-flatness), exclude a second
unramified nonproperness divisor, and straighten the
Deligne--Faddeev coefficient orbit.  Completing either gateway gives the
foundational map and reconstructs the other package, so they should not be
pursued as cumulative independent requirements.
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
the common missing lemma in the normalization and branchwise gateways.
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
The ranked execution plan and the matching phantom-boundary and graded
gauge certificates are maintained in
[`CUBIC_CLOSURE_ATTACKS.md`](CUBIC_CLOSURE_ATTACKS.md).
Cancellation arithmetic is parked as `OP-ARITH`.  The former roadmap is preserved in
[`archive/legacy-notes`](../archive/legacy-notes/CANCELLATION_RESEARCH_ROADMAP_2026-07-22.md).

The geometric deck group of the cancellation inverse polynomial is treated
separately in [`INVERSE_MONODROMY.md`](INVERSE_MONODROMY.md).  Its generic
`A_N/S_N` theorem and branch-cycle search are independent of the parameter-
polynomial Galois problem in `OP-ARITH`.

# Mathieu kernels of quadratic norm complements

This Lean 4.32.1 project follows the [uniform counterexample](../../extended-geometry/QUADRATIC_COMPLEMENT_UNIFORM_NORMALIZER_COUNTEREXAMPLES.md)
and the [sharp dimension classification](../../extended-geometry/QUADRATIC_COMPLEMENT_SHARP_DIMENSION_CLASSIFICATION.md).

Run both commands from this directory.

    lake build
    lake env lean Audit.lean

The dependency lock pins Mathlib v4.32.1. The audit prints the transitive
Lean axioms of the principal declarations.

## Statements in the current draft

The [paper draft](../../../../unpublished/mathieu-quadratic-complements/main.tex)
opens over the rational sum of squares with
\(p=(1+x)(1+2x-y^2-z^2)/q\). The main theorem uses the distinct product
family \(F_{α,β,γ}=q^{-α}(1+q)^β W^γ\).

- [SimplePole.lean](MathieuQuadratic/SimplePole.lean) proves the opening
  pure and marked moments, with marked value \(2^m/(m+1)\). On a
  normalized diagonal chart, it gives the dimension-dependent formula
  for the same construction. It derives the moment formula from
  scalar weights, rank-two rotations, and finite polynomial
  coefficient extraction. The geometric restriction at \(x=-1\)
  excludes proper powers up to units. The denominator order is one,
  every pure polynomial phase is safe, and a constant multiplier
  cannot detect any pure phase.
- [Moments.lean](MathieuQuadratic/Moments.lean) and
  [Draft.lean](MathieuQuadratic/Draft.lean) give the mixed-moment
  proposition for arbitrary product parameters, the zero for a negative
  extraction index, and the first omitted coefficient. They identify the
  rational beta constants with the integrals in the draft. The uniform
  prescription \(F_{M,M,M+1}\) has the fixed marker \(q^M\).
  [Nonpower.lean](MathieuQuadratic/Nonpower.lean) checks both
  restriction orders over the algebraic closure.
- [BlockPolynomials.lean](MathieuQuadratic/BlockPolynomials.lean) proves
  the block-invariant polynomial description and degree bound. It obtains
  a block average from invariant complements and the trace comparison.
  [DraftBounds.lean](MathieuQuadratic/DraftBounds.lean) clears the norm
  denominator for regular multipliers and proves the sharp cutoff and least
  detecting polynomial degrees from the block averaging input.
- [SharpMultiplier.lean](MathieuQuadratic/SharpMultiplier.lean) proves
  the last-order detection by \(B^J\), and the purity and detection
  alternatives for the product family given the block expansions.
  [MultiplierTail.lean](MathieuQuadratic/MultiplierTail.lean) gives the
  rational test for a finite expansion and proves both eventual outcomes.
  It checks the cancellation example for the uniform ternary phase.
- [Perturbation.lean](MathieuQuadratic/Perturbation.lean) gives the
  first-order cancellation and the all-order witness for \(F_*\).
  [PoleOrder.lean](MathieuQuadratic/PoleOrder.lean) proves the exact
  denominator orders, the sharp minimum within the product family, and
  the absolute nonpower assertion for the perturbation.
- [MarkedRecurrence.lean](MathieuQuadratic/MarkedRecurrence.lean)
  proves that a supplied nonzero polynomial recurrence gives a finite
  window for eventual vanishing. The proof chooses a threshold beyond
  the integer roots of both endpoint coefficients.
- [Instability.lean](MathieuQuadratic/Instability.lean) gives the
  eventual-vanishing bound from a positive weight filtration. It proves
  purity and safety for \(F_{M+1,M,2M+1}\) from the block average and
  shows that its two extreme moments give a nonzero invariant under
  the stated covariance. A closed invariant fiber then excludes zero
  from its orbit closure.
- [NormOneChart.lean](MathieuQuadratic/NormOneChart.lean) constructs the
  hyperbolic norm-one quotient. The torus action gives unequal-exponent
  vanishing, and the same rotation identity as on the complement gives the
  height moment formula. Together they give the exact all-order witness.
  A geometric simple zero verifies nonpower.
- [SplitNormTraces.lean](MathieuQuadratic/SplitNormTraces.lean) proves
  transitivity on a nonzero vector–covector pairing, the split
  quadratic étale model. It also proves that square detection implies
  the Mathieu property, as used for the real sphere.
- [AlgebraicTrace.lean](MathieuQuadratic/AlgebraicTrace.lean) constructs
  the normalized functional from polynomial moments. The radial
  recurrence gives a trace independent of the denominator
  presentation. The polynomial chain rule proves its scalar and
  orthogonal invariance. The same moment formula works over the dual
  numbers, giving infinitesimal invariance in
  [DualTrace.lean](MathieuQuadratic/DualTrace.lean).
- [OrthogonalTrace.lean](MathieuQuadratic/OrthogonalTrace.lean) verifies
  all the geometric trace hypotheses for the constructed functional
  in the plane-height chart. [Draft.lean](MathieuQuadratic/Draft.lean)
  obtains the uniform witness in every dimension at least three by
  diagonalization and transport. The marker remains \(q^{N-2}\).
- [GeometricAction.lean](MathieuQuadratic/GeometricAction.lean) and
  [Transitivity.lean](MathieuQuadratic/Transitivity.lean) prove geometric
  transitivity and the constant invariant subspace.
  [GeometricProjection.lean](MathieuQuadratic/GeometricProjection.lean)
  proves that the algebraic trace commutes with extension of the
  ground field.
  [TraceUniqueness.lean](MathieuQuadratic/TraceUniqueness.lean) proves
  uniqueness on geometric points by testing powers of linear forms.
  [CanonicalTrace.lean](MathieuQuadratic/CanonicalTrace.lean) transports
  this trace through arbitrary diagonal charts and identifies it with
  the normalized invariant functional over the algebraic closure.
- [AlgebraicLaurentTrace.lean](MathieuQuadratic/AlgebraicLaurentTrace.lean)
  identifies the algebraic trace with Laurent constant coefficient in
  dimensions one and two. The binary proof conjugates the split torus
  into the diagonal isometry group and uses compatibility with scalar
  extension. The Mathieu conclusion on each of these coordinate models
  uses Laurent necessity, verified here in rank one.
- [LaurentTrace.lean](MathieuQuadratic/LaurentTrace.lean) derives the
  weights and constant-term comparison from the actual diagonal torus
  pullbacks in the split low-dimensional charts.
  [QuadraticCoordinates.lean](MathieuQuadratic/QuadraticCoordinates.lean)
  gives the norm-preserving localization equivalence for the rank-two chart.
- [Draft.lean](MathieuQuadratic/Draft.lean) gives the dimension cutoff
  and field-extension result from the geometric trace comparisons
  and the complex Laurent theorem. Its high-dimensional conclusion
  uses the draft's product phase and prescribed norm-power marker.

The canonical trace theorem in [Draft.lean](MathieuQuadratic/Draft.lean)
proves the sharp dimension cutoff from the rank-two interior case of
Laurent necessity. Dimension one is verified in full. The counterexample in
every dimension at least three has the prescribed norm-power marker and uses
the canonical trace.

[RankTwoLaurent.lean](MathieuQuadratic/RankTwoLaurent.lean) verifies the
supporting-line reduction in Section 3 of Duistermaat–van der Kallen,
including lower-dimensional Newton polygons. It uses
[LaurentFace.lean](MathieuQuadratic/LaurentFace.lean) to preserve the constant
terms of all powers when restricting to the zero-weight face and projecting
that face to one Laurent coordinate.
[TorusLaurent.lean](MathieuQuadratic/TorusLaurent.lean) checks the initial
torus integral for the constant-term generating function. The interior case
still requires the compactification and continuation argument.

[OneVariableLaurent.lean](MathieuQuadratic/OneVariableLaurent.lean) derives
rank-one Laurent necessity from the one-variable coefficient extraction
proof in [GMC-2-lean at commit 1782de7](https://github.com/MurrellGroup/GMC-2-lean/tree/1782de7ff6c97eb1d98e63e7ff34df18b9cd322e).
[The adapted MIT-licensed modules](MathieuQuadratic/External/DvK)
retain the Laurent root-sum proof, omit the Gaussian supporting-face
construction and use the pinned Mathlib's coefficient and valuation
APIs. The audit checks the imported coefficient extraction and the
resulting rank-one theorem on the current pin.

The recurrence proposition still needs the residue-elimination proof
of existence and effectivity. The nullcone classification still needs
the geometric Hilbert–Mumford comparison and the two-term
ordered-field positivity argument. The Hermitian corollary still
needs the algebra splitting, the quaternionic symplectic case, and
the canonical trace comparisons. The norm-one corollary uses torus
and infinitesimal invariance as hypotheses. The real sphere observation
uses square detection as a hypothesis.

[AngularEndpoints.lean](MathieuQuadratic/AngularEndpoints.lean) verifies
the draft's angular endpoint proposition using termwise integration over
the ground field. The endpoint at height ratio one is zero or has strictly
one-sided Laurent support. In even dimension the same conclusion holds at
height ratio zero. [Draft.lean](MathieuQuadratic/Draft.lean) transports both
tests to the canonical trace through its orthogonal block chart.

[MarkedRecurrence.lean](MathieuQuadratic/MarkedRecurrence.lean) verifies the
single-order test once the finite spanning family of markers and invertible
moment recurrence are available. Establishing the twisted de Rham dimension
bound and effective basis construction remains open in the formalization.
The quadratic-numerator and affine-angular classifications also remain open.

## Supporting squared-normalizer construction

The coordinate ring in [Complement.lean](MathieuQuadratic/Complement.lean)
is the localization of the multivariate polynomial ring away from the norm.
[Diagonal.lean](MathieuQuadratic/Diagonal.lean) obtains the rank-two and
rank-\(M\) orthogonal splitting from a nondegenerate quadratic form of
dimension \(M+2\) and identifies its norm polynomial.

- [Scalar.lean](MathieuQuadratic/Scalar.lean) constructs scalar pullback on
  the localization and derives the vanishing at nonzero scalar weights.
  [Rotations.lean](MathieuQuadratic/Rotations.lean) defines the rank-two
  derivations, shows that they kill the norm, and proves the summed
  differential identity and height recurrence.
  [Infinitesimal.lean](MathieuQuadratic/Infinitesimal.lean) realizes each
  rotation on the norm localization over the dual numbers.
- [Endpoint.lean](MathieuQuadratic/Endpoint.lean) proves the rational beta
  integral, positivity and factorial formula, and the endpoint divisibility.
  [Moments.lean](MathieuQuadratic/Moments.lean) extracts the mixed moments
  with the generalized binomial for every \(\ell\ge-k\), including \(k<M\).
  Scalar weights give zero for \(\ell>k\).
- [Uniform.lean](MathieuQuadratic/Uniform.lean) proves the pure-power
  threshold, the eventual and all-order radial marker thresholds, and the
  exact marked sum for \(F=P^M(1+P)\) with \(H=q^M\).
- [Nonpower.lean](MathieuQuadratic/Nonpower.lean) verifies the simple zero in
  the algebraically closed scalar extension using dual numbers. The
  uniform phase has a simple zero at norm \(-1/2\), and the direct phase
  \(P\) has one at norm \(-1\). A unit times a power of exponent at least two
  cannot have either image.
- [Counterexample.lean](MathieuQuadratic/Counterexample.lean) assembles the
  explicit uniform and ternary witnesses on the localized ring.
  Its trace hypothesis requires normalization, invariance under scalar two,
  and annihilation of the rank-two rotation images with norm-power
  denominators. [Reynolds.lean](MathieuQuadratic/Reynolds.lean) constructs
  the unique normalized linear projection from a direct sum of the constants
  and an augmentation subspace, and derives the trace hypotheses when the
  augmentation contains the scalar differences and rotation images.
- [Extension.lean](MathieuQuadratic/Extension.lean) identifies the
  scalar-extended localization and verifies that the same phase and marker
  formulas hold there.

## Dimension classification

[BaseChange.lean](MathieuQuadratic/BaseChange.lean) proves that purity,
safety, and the all-order witness property survive scalar extension in both
directions for a fixed phase and marker from the original field. It proves
descent of the Mathieu property by tensor product induction on the markers
after extension.

[Laurent.lean](MathieuQuadratic/Laurent.lean) proves Laurent purity and
eventual safety from a positive integral separator, with a cutoff for each
marker. [CoefficientField.lean](MathieuQuadratic/CoefficientField.lean)
transfers the complex necessity theorem to every characteristic-zero
field using the field generated by the coefficients of one phase.
[Newton.lean](MathieuQuadratic/Newton.lean) relates integral separation
to exclusion of zero from the real Newton polytope, so the complex
necessity hypothesis has the form stated by Duistermaat and van der Kallen.

[SplitLocalization.lean](MathieuQuadratic/SplitLocalization.lean) identifies
the localization away from a monomial with the Laurent polynomial algebra.
[LowDimension.lean](MathieuQuadratic/LowDimension.lean) constructs the
unary and split binary coordinates. [LaurentTrace.lean](MathieuQuadratic/LaurentTrace.lean)
uses the scalar character in dimension one and the two characters
\(i+j\), \(i-j\) in dimension two to derive the constant-term trace and its
Mathieu property from character invariance and complex Laurent necessity.
[BinarySupport.lean](MathieuQuadratic/BinarySupport.lean) proves the
scalar-weight sign criterion for a support invariant under swapping the
two Laurent coordinates.
[Classification.lean](MathieuQuadratic/Classification.lean) transports
the resulting low-dimensional criterion and the higher-dimensional
witness through the indicated algebra and trace identifications.

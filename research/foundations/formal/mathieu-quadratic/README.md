# Mathieu kernels of quadratic norm complements

This Lean 4.32.1 project follows the [uniform counterexample](../../extended-geometry/QUADRATIC_COMPLEMENT_UNIFORM_NORMALIZER_COUNTEREXAMPLES.md)
and the [sharp dimension classification](../../extended-geometry/QUADRATIC_COMPLEMENT_SHARP_DIMENSION_CLASSIFICATION.md).

From this directory, run:

    lake build
    lake env lean Audit.lean

The dependency lock pins Mathlib v4.32.1. The audit prints the transitive
Lean axioms of the principal declarations.

## Statements in the current draft

The [paper draft](../../../../unpublished/mathieu-quadratic-complements/main.tex)
uses the product phase \(F_{α,β,γ}=q^{-α}(1+q)^β W^γ\). The formal
statements follow this choice of phase.

- [Moments.lean](MathieuQuadratic/Moments.lean) and
  [Draft.lean](MathieuQuadratic/Draft.lean) give the mixed-moment
  proposition for arbitrary product parameters, the zero for a negative
  extraction index, and the first omitted coefficient. They identify the
  rational beta constants with the integrals in the draft.
- [Draft.lean](MathieuQuadratic/Draft.lean) gives the opening ternary
  witness \(F_{1,1,1}\), its first marked moment \(2/3\), and the
  prescription \(F_{M,M,M+1}\) with marker \(q^M\). These statements use
  a normalized trace with scalar and geometric infinitesimal invariance.
  [Nonpower.lean](MathieuQuadratic/Nonpower.lean) checks the two
  restriction orders over the algebraic closure.
- [SharpMultiplier.lean](MathieuQuadratic/SharpMultiplier.lean) proves
  the sharp order cutoff on the span of the averaged monomials, the
  last-order detection by \(B^J\), and the purity and detection
  alternatives for the product family given the block expansions.
  [MultiplierTail.lean](MathieuQuadratic/MultiplierTail.lean) gives the
  cleared rational function for a finite expansion and proves the
  equivalence with eventual vanishing.
- [Perturbation.lean](MathieuQuadratic/Perturbation.lean) gives the
  first-order cancellation, all positive pure and marked moments,
  and below-threshold vanishing on the block span for \(F_*\).
  [PoleOrder.lean](MathieuQuadratic/PoleOrder.lean) proves the exact
  clearing powers for both phases and assembles the perturbation witness
  with its absolute nonpower property.
- [OrthogonalTrace.lean](MathieuQuadratic/OrthogonalTrace.lean) realizes
  the rank-two rotations as special orthogonal dual-number points.
  [GeometricAction.lean](MathieuQuadratic/GeometricAction.lean) and
  [Transitivity.lean](MathieuQuadratic/Transitivity.lean) prove geometric
  transitivity and the constant invariant subspace.
  [QuadraticCoordinates.lean](MathieuQuadratic/QuadraticCoordinates.lean)
  gives the norm-preserving localization equivalence for the rank-two chart.
- [Draft.lean](MathieuQuadratic/Draft.lean) proves the stated dimension
  cutoff and field-extension result given the trace comparisons and
  the complex Laurent theorem. The high-dimensional conclusion uses
  the draft's product phase and prescribed norm-power marker.

The geometric Reynolds decomposition, its trace comparisons, the
unconditional block averaging proposition, the isotropic ternary linear
multiplier example, and the complex Laurent necessity theorem remain
formalization obligations.

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

The geometric Reynolds decomposition for the \(SO(q)\times\mathbb G_m\)
action remains an input. The geometric trace comparisons must identify
the canonical group trace with the normalized projections and torus
character traces in these modules. The complex Duistermaat–van der Kallen
necessity direction remains an explicit hypothesis of the classification.

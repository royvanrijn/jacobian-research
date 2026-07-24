# Lean formalization of GMC(2)

This package formalizes the complete two-variable theorem proved in
[`TWO_REAL_GMC_LOWER_FACE_THEOREM.md`](../../extended-geometry/TWO_REAL_GMC_LOWER_FACE_THEOREM.md).

Build it with:

```sh
lake build
```

The checked modules are:

- `SupportingFace`: the rational lower-face certificate;
- `WeightedInitial`, `CircularFlatten`, and `LowerFaceExtraction`: the
  associated-graded lower-face calculation and exposed lowest coefficient;
- `DuistermaatVanDerKallen`: the imported constant-term theorem;
- `Specialization`: the good-reduction interface for finite-type
  `ℤ`-domains;
- `CoefficientRing` and `SpecializationArgument`: descent to the finite
  coefficient ring and the full reduction-modulo-\(p\) contradiction;
- `FactorialQuotient` and `PrimeIsolation`: the factorial divisibility and
  characteristic-`p` isolation argument;
- `ConstantTerm`: Frobenius commuting with angular constant-term extraction;
- `OneSided` and `LowerFace`: the one-sided-support conclusion and exported
  circular theorem;
- `BivariateGaussian`: the ordinary bivariate-polynomial statement and
  Wick monomial formula.

There are no `sorry` declarations.  The rational supporting-face extraction
is proved directly from the finite weighted support.  Two deep inputs are
explicitly marked as `axiom`: DvdK and finite-type specialization.  They are
isolated interfaces rather than hidden proof gaps and can be replaced
locally as the relevant mathlib APIs mature.

The public theorem is
`GMC2.gaussianMomentsConjecture_two_variables_bivariate`.  Its axiom audit
contains exactly those two mathematical inputs, in addition to Lean's
standard logical choice, quotient, and propositional-extensionality
principles.

# Lean formalization of the GMC(2) lower-face argument

This package formalizes the modular core of
[`TWO_REAL_GMC_LOWER_FACE_THEOREM.md`](../../extended-geometry/TWO_REAL_GMC_LOWER_FACE_THEOREM.md).

Build it with:

```sh
lake build
```

The checked modules are:

- `SupportingFace`: the rational lower-face certificate;
- `DuistermaatVanDerKallen`: the imported constant-term theorem;
- `Specialization`: the good-reduction interface for finite-type
  `ℤ`-domains;
- `FactorialQuotient` and `PrimeIsolation`: the factorial divisibility and
  characteristic-`p` isolation argument;
- `ConstantTerm`: Frobenius commuting with angular constant-term extraction;
- `LowerFace`: the final contradiction and eventual one-sided-weight lemma.

There are no `sorry` declarations.  Three deep inputs are explicitly marked
as `axiom`: DvdK, finite-type specialization, and rational supporting-face
extraction.  The latter two are isolated interfaces rather than hidden proof
gaps, and can be replaced locally as mathlib's polyhedral and
scheme-specialization APIs mature.


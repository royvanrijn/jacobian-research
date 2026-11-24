# The Mathieu property for quadratic complements

**Unpublished manuscript, ready for mathematical review**

The normalized rotation-and-scaling trace on a nondegenerate quadratic
complement in characteristic zero has a Mathieu–Zhao kernel exactly in
dimensions one and two. In every higher dimension, a linear multiplier
detects every order of an explicit pure phase. The paper also gives the
sharp degree and cutoff for its original uniform witness.

Polynomial recurrences decide eventual vanishing for each fixed phase in
the block-invariant ring and each fixed regular multiplier. The geometric
nullcone is strictly smaller than the pure phases with eventual marked
vanishing in dimensions at least three.

The same trace classification gives explicit witnesses for quadratic-étale
and quaternionic Hermitian norm complements. Isotropic norm-one quadrics
also have witnesses over the ground field. For positive-definite real
norm-one quadrics of dimension at least two, the real trace is Mathieu–Zhao
and the complex trace fails.

[PDF](mathieu-quadratic-complements.pdf) · [LaTeX source](main.tex) ·
[Metadata](metadata.json) · [Verification record](verification.json)

## Build

Run from this directory with TeX Live.

```sh
latexmk -pdf -interaction=nonstopmode -halt-on-error -jobname=mathieu-quadratic-complements main.tex
```

The single source file contains the bibliography and builds on its own.

## Proof and review

The written arguments are in the paper and the supporting research notes
on the [dimension classification](../../research/foundations/extended-geometry/QUADRATIC_COMPLEMENT_SHARP_DIMENSION_CLASSIFICATION.md),
[rank-two moments](../../research/foundations/extended-geometry/QUADRATIC_COMPLEMENT_UNIFORM_NORMALIZER_COUNTEREXAMPLES.md),
[sharp multiplier bounds](../../research/foundations/extended-geometry/QUADRATIC_COMPLEMENT_SHARP_MULTIPLIER_BOUNDS.md),
and [isotropic counterexamples](../../research/foundations/extended-geometry/ISOTROPIC_QUADRATIC_COMPLEMENT_COUNTEREXAMPLES.md).
The verification record pins the revised edition and records the bounded
arithmetic checks and PDF-only internal review. External specialist review
and independent external replay remain open.

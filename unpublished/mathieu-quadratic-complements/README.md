# The Mathieu property for quadratic complements

**Unpublished manuscript, ready for mathematical review**

The normalized rotation-and-scaling trace on a nondegenerate quadratic
complement in characteristic zero has a Mathieu–Zhao kernel exactly in
dimensions one and two. The paper opens with a double-pole ternary
example and gives a uniform witness in higher dimensions. For this
witness the least polynomial degree detecting infinitely many orders is
\(2N-4\). The paper gives the sharp lower-degree cutoff, a test for fixed
regular multipliers, and a pole-reducing perturbation in dimensions at
least six. The introduction credits Long's spherical classification
for the complex case. An isotropic ternary example illustrates how the
least detecting polynomial degree can change with the choice of witness.

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

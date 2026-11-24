# The Mathieu property for quadratic complements

**Unpublished manuscript, ready for mathematical review**

For every nondegenerate quadratic form of positive dimension over a
characteristic-zero field, the kernel of the normalized rotation-and-scaling
trace on its nonzero locus is a Mathieu–Zhao subspace exactly in dimensions
one and two. The manuscript opens with an explicit ternary counterexample
and gives a uniform regular witness detected at every order by one
fixed multiplier in higher dimensions.

[PDF](mathieu-quadratic-complements.pdf) · [LaTeX source](main.tex) ·
[Metadata](metadata.json) · [Verification record](verification.json)

## Build

Run from this directory with TeX Live.

```sh
latexmk -pdf -interaction=nonstopmode -halt-on-error -jobname=mathieu-quadratic-complements main.tex
```

The single source file contains the bibliography and builds on its own.

## Proof and review

The paper contains the full written argument. Its source proofs are the
[dimension classification](../../research/foundations/extended-geometry/QUADRATIC_COMPLEMENT_SHARP_DIMENSION_CLASSIFICATION.md) and the
[uniform moment formula](../../research/foundations/extended-geometry/QUADRATIC_COMPLEMENT_UNIFORM_NORMALIZER_COUNTEREXAMPLES.md). The exact replay compares
selected moment coefficients with direct Laurent expansions and includes
a negative control. The universal theorem follows from the written proof.
The verification record identifies the reviewed edition, the replay scope,
and the revisions from the isolated PDF-only review. External specialist
review and independent replay remain open.

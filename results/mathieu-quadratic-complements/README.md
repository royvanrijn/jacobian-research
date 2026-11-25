# The Mathieu property for quadratic complements

**Under review at Hexagon.** Reserved identifier `hexagon:2611.01050v2`. Updated November 25, 2026.

The normalized rotation-and-scaling trace on a nondegenerate quadratic
complement in characteristic zero has a Mathieu–Zhao kernel exactly in
dimensions one and two. In every higher dimension, one phase attains pole
order one, multiplier degree one, and numerator degree three. These are
the three minima. For a numerator of degree at most two over \(q\), the
first \(N\) pure moments test purity, and every pure phase contracts over
the ground field.

Eventual vanishing against every fixed multiplier for a block-invariant
phase affine in \(S=B/q\) is equivalent to one-sided scalar support,
apart from the zero phase. The angular endpoint test also gives necessary
conditions for general block-invariant phases. For any nonzero regular
phase \(f=p/q^a\), a twisted de Rham computation gives a list of at most
\(D^N\) multipliers, where \(D\) is the degree of the squarefree part of
\(qp\). Testing that list at one computable order decides eventual
vanishing against every fixed regular multiplier.

The same trace classification gives explicit witnesses for quadratic-étale
and quaternionic Hermitian norm complements. Isotropic norm-one quadrics
also have witnesses over the ground field. Over the reals, the kernel of
the special orthogonal Reynolds projection on the complement is
Mathieu–Zhao exactly when the form is definite or \(N\le2\).

[PDF](mathieu-quadratic-complements.pdf) · [LaTeX source](main.tex) ·
[Metadata](metadata.json) · [Verification record](verification.json) ·
[HexagonMath source ZIP](mathieu-quadratic-complements-hexagon-source.zip)

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

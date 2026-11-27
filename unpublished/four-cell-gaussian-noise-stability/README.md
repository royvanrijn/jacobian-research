# The tetrahedral partition maximizes Gaussian noise stability for small positive correlations

Roy van Rijn · prepared 27 November 2026 · ready for specialist review

The paper proves a qualitative, dimension-independent interval of positive correlations on which the central regular tetrahedral partition maximizes Gaussian noise stability among all balanced ordinary measurable four-cell partitions. Equality is classified up to orthogonal maps, labels, cylindrical extension and null sets. The main statement is for ambient dimension at least three. The threshold is not explicit, and the paper does not address every positive correlation.

[Read the paper](four-cell-gaussian-noise-stability.pdf) · [LaTeX source](main.tex) · [Source archive](source.zip) · [Metadata](metadata.json) · [Verification and blind-reader record](verification.json) · [Checksums](SHA256SUMS.txt)

## Proof and dependencies

The article contains every new argument needed for its theorem. Its chain is:

1. Zhang's sharp equal-cell first-moment theorem and balanced covariance deficit identity, stated precisely in Proposition 2.1; Milman–Neeman supplies the Gaussian multi-bubble perimeter input used there.
2. Quantitative covariance, moment and partition alignment (Proposition 2.2), including exact orthogonality to the closest regular moment orbit.
3. The first-moment interface form and its exact rotational kernel (Lemma 3.1).
4. Saturated tetrahedral Ornstein–Uhlenbeck score margins, a Gaussian weighted exterior estimate and the normalized Mehler kernel bound (Section 4).
5. Regret-normalized normal-flux compactness, its sharp one-sided cost and interaction limit (Proposition 5.1).
6. Three-dimensional rigidity (Theorem 6.1), followed by Heilman–Tarter maximizer existence and dimension reduction to obtain Theorem 1.1 in every ambient dimension at least three.

The stronger maximizer localization and analytic interface chart remain in the [research companion](../../research/plurality-stability/MAXIMIZER_PERTURBATIVE_RIGIDITY.md). They are not premises of the shorter proof in this article. Section 7 retains the tail-swap obstruction to a noise-stability deficit proportional to the first-moment deficit.

The current research proofs and exact dependency scopes remain in [quantitative first-moment stability](../../research/plurality-stability/QUANTITATIVE_SIMPLEX_FIRST_MOMENT_STABILITY.md), [interface rigidity](../../research/plurality-stability/MAXIMIZER_PERTURBATIVE_RIGIDITY.md), [the small-correlation proof](../../research/plurality-stability/SMALL_POSITIVE_CORRELATION_INTERVAL.md), and [MATH_STATUS.json](../../research/MATH_STATUS.json). The article is their standalone presentation, not an independent verification of those arguments.

Zhang's input is a frozen external preprint, release `gaussian-balanced-simplex-all-k-v1`, commit `ebfd90558439afca21babeb7b77c8d8837ea745f`. The [retained text](../../research/plurality-stability/sources/zhang-equal-gaussian-cells-ebfd905.md) and [source/application ledger](../../research/plurality-stability/sources/zhang-equal-gaussian-cells-sources-ebfd905.md) distinguish the frozen project author field from the signed companion's personal attribution. No peer-review status is inferred for that preprint.

## Build

From this directory, with TeX Live, AMS packages, `microtype`, `mathtools`, `xcolor`, `hyperref` and `latexmk` installed:

```sh
SOURCE_DATE_EPOCH=1795737600 FORCE_SOURCE_DATE=1 latexmk -pdf -interaction=nonstopmode -halt-on-error -jobname=four-cell-gaussian-noise-stability main.tex
```

The bibliography is included in `main.tex`; no repository files or network access are needed to compile the source archive. The designated PDF has ten pages. The epoch fixes PDF creation metadata for artifact reproduction; it has no mathematical role.

## Review status

The current edition has a warning-free build, an all-page visual review and an isolated internal AI PDF-only reading. The reader found no concrete gap in the score-regret normalization, weighted tightness, flux coefficient, interaction limit or removal of rotations. The verification record preserves the edition hash, five reader answers, limitations and optional exposition suggestions.

External specialist review, formal verification and independent mathematical replay have not occurred. The imported theorem applications and the new analytic proof remain subject to that review. This paper has not been submitted, released or published. No larger-correlation research was undertaken in preparing this edition.

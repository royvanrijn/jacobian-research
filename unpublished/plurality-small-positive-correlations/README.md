# Plurality Is Stablest for Small Positive Correlations: A Gaussian Simplex Theorem

Roy van Rijn · revised 27 November 2026 · prepared for specialist review

For every fixed number of cells $k\ge2$, the paper proves that the central
regular-simplex partition maximizes Gaussian noise stability among all
balanced measurable partitions for a qualitative interval
$0<\rho\le\rho_k$. The threshold may depend on $k$ but is independent of
the ambient dimension $d\ge k-1$. Equality is classified up to orthogonal
maps, labels, cylindrical extension and null sets. The paper does not
assert a numerical threshold, a threshold uniform in $k$, or optimality
throughout $0<\rho<1$.

[Read the paper](plurality-small-positive-correlations.pdf) ·
[LaTeX source](main.tex) · [Source archive](source.zip) ·
[Metadata](metadata.json) · [Verification and reader record](verification.json) ·
[Checksums](SHA256SUMS.txt)

## Argument and mathematical dependencies

The article contains the complete essential argument for every fixed $k$.
Section 2 reproduces Zhang's covariance-path proof with attribution,
including balancing prices, Gaussian flux, the ordinary-perimeter
normalization, singular covariance endpoints and equality. The geometric
input is Milman–Neeman's Gaussian multi-bubble theorem in dimension $k-1$.
The quantitative alignment proof uses compactness of a uniformly positive
covariance set instead of unnecessary explicit constants.

Sections 3–5 identify the rotational interface kernel and prove weighted
normal-flux compactness for arbitrary measurable transfers. The sharp
one-sided cost and compact-cutoff interaction limit are displayed
explicitly. Section 6 gives the contradiction, proves existence with the
mass constraints retained, and uses Heilman–Tarter dimension reduction
only for actual maximizers. Section 7 replaces a four-cell tail box by a
single reflected remote ball, showing for every $k$ why first-moment
deficit cannot be the uniform scale of a full noise-value gap.

Zhang's argument is credited to the frozen external preprint
`gaussian-balanced-simplex-all-k-v1`, commit
`ebfd90558439afca21babeb7b77c8d8837ea745f`. The
[retained source](../../research/plurality-stability/sources/zhang-equal-gaussian-cells-ebfd905.md)
and [attribution ledger](../../research/plurality-stability/sources/zhang-equal-gaussian-cells-sources-ebfd905.md)
give its exact edition and authorship evidence. The paper supplies that
argument for reader verification; it does not claim authorship of the
first-moment comparison.

The [first-moment stability proof](../../research/plurality-stability/QUANTITATIVE_SIMPLEX_FIRST_MOMENT_STABILITY.md)
and [fixed-cell normal-flux proof](../../research/plurality-stability/SMALL_POSITIVE_CORRELATION_INTERVAL.md)
remain the research sources. The separate
[continuation results](../../research/plurality-stability/POSITIVE_CORRELATION_CONTINUATION.md)
are not premises of this article. The full positive-correlation conjecture
remains open in the current mathematical account.

## Build and review

With TeX Live, AMS packages, `geometry`, `microtype`, `mathtools`, `xcolor`,
`hyperref` and `latexmk` installed, run from this directory:

```sh
SOURCE_DATE_EPOCH=1795737600 FORCE_SOURCE_DATE=1 latexmk -pdf -interaction=nonstopmode -halt-on-error -jobname=plurality-small-positive-correlations main.tex
```

The bibliography is included in `main.tex`. The source archive builds
without repository files or network access. The designated PDF has nine
pages. The fixed epoch controls PDF metadata for artifact reproduction;
it has no mathematical role.

The current source has a warning-free build and an all-page visual
inspection. The exact PDF-only reader test, its edition, findings and
any resulting revisions are recorded in [verification.json](verification.json),
together with the clean archive build and repository checks. Earlier
four-cell review outcomes do not certify this enlarged edition.

External specialist review, formal verification and independent
mathematical replay have not occurred. No submission, release or
publication is asserted. The present edition is the single current paper
home; the previous four-cell edition remains in Git history.

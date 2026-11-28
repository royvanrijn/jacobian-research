# Quantitative Gaussian Simplex Stability at Small Positive Correlations

Roy van Rijn · revised 28 November 2026 · prepared for specialist review

For every fixed $k\ge2$, the partition $T$ of a centered regular simplex maximizes
Gaussian noise stability among balanced measurable $k$-cell partitions for
all sufficiently small positive correlations. More precisely,
\[
 S_\rho(T)-S_\rho(A)\ge c_k\rho\,d_k(A)^2,
\]
where $d_k(A)$ is the infimum of the sum of Gaussian symmetric differences
over the labeled simplex orbit. Both $c_k>0$ and the correlation threshold
are independent of the ambient dimension $d\ge k-1$. Both constants may
depend on $k$ and are left unevaluated. Equality of noise stabilities
identifies the simplex orbit.

[Read the paper](plurality-small-positive-correlations.pdf) ·
[LaTeX source](main.tex) · [Source archive](source.zip) ·
[Metadata](metadata.json) · [Verification and reader record](verification.json) ·
[Checksums](SHA256SUMS.txt) · [Hexagon submission fields](#hexagon-submission)

## Proof

The paper proves the comparison directly for arbitrary measurable
competitors. Section 2 presents Yongxian Zhang's first-moment comparison along an
additive covariance path using the Milman–Neeman Gaussian multi-bubble
theorem. The subsequent compactness and score argument gives uniform
localization in the nearest simplex's moment coordinates.
[Signed first-moment preprint](https://doi.org/10.5281/zenodo.23250730)

Section 3 derives the explicit aligned facet inequality from Prékopa's
marginal theorem. Sections 4 and 5 use the exact Ornstein–Uhlenbeck score
loss to bound arbitrary label changes and identify the limiting facet
densities. Section 6 controls the higher Hermite degrees uniformly across
the active and transverse coordinates and combines the local and global
estimates. The argument avoids dimension reduction for maximizers.
Section 7 explains the score-loss normalization using a reflected
exchange far inside two cells.

The [underlying analytic result](../../research/plurality-stability/SMALL_POSITIVE_CORRELATION_INTERVAL.md#dimension-independent-local-comparison)
contains the direct local and global comparisons. The
[frozen Zhang source](../../research/plurality-stability/sources/zhang-equal-gaussian-cells-ebfd905.md)
and [attribution ledger](../../research/plurality-stability/sources/zhang-equal-gaussian-cells-sources-ebfd905.md)
identify the first-moment preprint. The introduction records Heilman's
balanced three-cell comparison through $\rho=2/5$. The present theorem gives a full label-distance deficit and
applies to every fixed $k$, including $k\ge4$.
[Heilman, Theorems 1.2 and A.1](https://arxiv.org/abs/2608.00333v1)

## Build and review

With TeX Live, AMS packages, `geometry`, `microtype`, `mathtools`, `xcolor`,
`hyperref` and `latexmk` installed, run from this directory or from a clean
extraction of the source archive:

```sh
SOURCE_DATE_EPOCH=1795824000 FORCE_SOURCE_DATE=1 latexmk -pdf -interaction=nonstopmode -halt-on-error -jobname=plurality-small-positive-correlations main.tex
```

The source archive contains the complete LaTeX and bibliography. The
designated PDF is nine pages and reproduces byte for byte with the command
above on the recorded toolchain. The fixed epoch makes the PDF reproducible.

The [verification record](verification.json) documents the current build,
layout inspection and the contributor-reported review of commit
`dd646006`, with its findings and the changes made in response.
Specialist mathematical review, independent mathematical replay, formal
verification and submission remain outstanding.

## Hexagon submission

Open the [Hexagon submission form](https://hexagonmath.org/submit) with ORCID
and upload [source.zip](source.zip). Its only entry is `main.tex`, with
the bibliography included. The archive reproduces the nine-page review PDF.
The same form values are available as [Hexagon metadata](hexagon-metadata.json).
[Hexagon deposit guidance](https://hexagonmath.org/about)

| Field | Value |
|---|---|
| Title | Quantitative Gaussian Simplex Stability at Small Positive Correlations |
| Contributor | Roy van Rijn |
| ORCID | `0009-0000-2808-0645` |
| Affiliation | Independent researcher |
| Primary category | `math.PR` — Probability |
| Cross-list | `math.MG` — Metric Geometry |
| Primary MSC 2020 | `60G15` — Gaussian processes |
| Secondary MSC 2020 | `49K21` — Optimality conditions for problems involving relations other than differential equations<br>`28A75` — Length, area, volume, other geometric measure theory |
| Keywords | Gaussian noise stability, regular simplex, standard simplex conjecture, first moments |
| AI tools | ChatGPT (OpenAI), Codex (OpenAI) |
| Author comments | 9 pages. |
| Version note, if offered | Initial deposit |

### Abstract

```text
For each fixed number of cells $k\ge2$, we prove that the partition defined by a centered regular simplex maximizes Gaussian noise stability among equal-measure measurable $k$-cell partitions at sufficiently small positive correlations $\rho$ in dimensions $d\ge k-1$. The noise stability deficit is bounded below by a positive constant times $\rho$ times the square of the summed Gaussian symmetric differences from the closest simplex partition. The constant and the correlation threshold are independent of the ambient dimension. Zhang's first-moment theorem localizes near maximizers. To prove the positive-correlation bound, we normalize arbitrary label changes by their exact loss against the simplex's Ornstein--Uhlenbeck scores. A Gaussian covariance estimate on the facets bounds the limiting first-moment interaction, and a Hermite expansion controls all remaining directions uniformly in the dimension.
```

### AI provenance

```text
ChatGPT (OpenAI) assisted with mathematical arguments, literature checks, exposition, proof review and verification code. OpenAI Codex performed a source-aware proof audit and an isolated reading of the compiled manuscript.
```

Confirm the contributor and authorship entries, text origin, extent of
human understanding, and rights at deposit. Choose the distribution
license there. CC BY 4.0 matches the author's earlier Hexagon deposits.
[Prior deposit](https://hexagonmath.org/2611.00596v1)

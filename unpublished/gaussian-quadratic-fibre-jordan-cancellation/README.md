# Gaussian Moment Vanishing with One Base Variable and Constant Quadratic Fibres

**Unpublished · partial draft.**

Classifies Gaussian phases with one unrestricted ordinary base variable and a constant quadratic fibre matrix, including arbitrary Jordan blocks. Purity, safety and a polynomial unit period are equivalent to finite rational identities. The determinant-one precision matrix gives a sufficient mixed-moment bound.

[PDF](gaussian_moments_one_base_constant_quadratic_fibres.pdf) · [LaTeX](main.tex) · [Metadata](metadata.json) · [Verification](verification.json) · [Source bundle](gaussian_moments_one_base_source.zip)

The [dimension-degree boundary paper](../gaussian-boundary/) now focuses on ternary cubics. Its former optional semisimple-fibre discussion is covered here by the stronger all-Jordan theorem and a short scalar pole-cancellation corollary. The [marker-reduction paper](../gaussian-constant-fibre-base-marker-certificates/) combines the complementary base-marker theorem and exact coefficient-ring obstruction.

Build from this directory:

```sh
latexmk -pdf -jobname=gaussian_moments_one_base_constant_quadratic_fibres -interaction=nonstopmode -halt-on-error main.tex
```

Only `main.tex` is needed to compile. The verification record identifies the current source and PDF, the tool versions, and the exact scope of the internal checks. Historical checker receipts retained in the source bundle keep their original scope.

Independent mathematical review remains pending. The purity-to-safety question for one constant quadratic fibre over two ordinary base variables remains open.

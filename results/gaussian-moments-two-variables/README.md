# The Gaussian Moments Conjecture in Two Variables

Proves the Gaussian Moments Conjecture for two real variables by constant-term extraction and prime isolation.

**Published preprint · 2026-07-23 · [10.5281/zenodo.21534699](https://doi.org/10.5281/zenodo.21534699)**

[PDF](gaussian-moments-two-variables.pdf) · [LaTeX](main.tex) · [Zenodo](https://zenodo.org/records/21534699) · [Submission metadata](zenodo-metadata.json) · [Publication metadata](metadata.json)

The PDF and deposited bundle files are the exact Zenodo edition. [Checksums](deposited-files.json) record their provenance.

A [working revision](../../unpublished/gaussian-moments-two-variables-revision/) is kept separately; it is not the deposited edition.

## Verification

Complete bivariate Lean development without sorry, conditional on two named mathematical inputs: Duistermaat-van der Kallen and finite-type specialization. These two inputs remain explicit axioms; independent external review is not claimed.

- [Evidence 1](../../research/foundations/formal/gmc2/README.md)

Build the deposited source from this directory by running `pdflatex -interaction=nonstopmode -halt-on-error main.tex` twice, using the supplied main.bbl. The deposited source archive does not include references.bib. Formal packages and broader research notes remain under `research/`.

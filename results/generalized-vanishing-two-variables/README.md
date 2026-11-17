# The Generalized Vanishing Conjecture: The Two-Variable Theorem and the First Failing Dimension

Proves GVC in two variables and gives a homogeneous counterexample in dimension three, determining the exact dimensional boundary.

**Published preprint · 2026-08-05 · [10.5281/zenodo.21803325](https://doi.org/10.5281/zenodo.21803325)**

[PDF](generalized-vanishing-two-variables.pdf) · [LaTeX](main.tex) · [Zenodo](https://zenodo.org/records/21803325) · [Hexagon](https://hexagonmath.org/2610.00580v2) · [Submission metadata](zenodo-metadata.json) · [Publication metadata](metadata.json)

The PDF and deposited bundle files are the exact Zenodo edition. [Checksums](deposited-files.json) record their provenance.

A [working revision](../../unpublished/generalized-vanishing-two-variables-revision/) is kept separately; it is not the deposited edition.

## Verification

Written proof with a claim-by-claim audit. The ternary counterexample and full profile family have Lean proofs. Binary GVC is only partly formalized. The complete binary theorem is not fully Lean-verified; independent external review is not claimed.

- [Evidence 1](../../unpublished/generalized-vanishing-two-variables-revision/VERIFICATION.md)
- [Evidence 2](../../research/foundations/formal/gvc/README.md)

Build the deposited source from this directory with `latexmk -pdf main.tex`. Formal packages and broader research notes remain under `research/`.

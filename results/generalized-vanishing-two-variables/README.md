# The Generalized Vanishing Conjecture: The Two-Variable Theorem and the First Failing Dimension

The current repository edition proves binary GVC for arbitrary constant-coefficient operators, extracts finite support certificates and a sharp operator-independent cutoff, and gives homogeneous counterexamples in every dimension at least three.

[Current PDF](generalized-vanishing-two-variables.pdf) | [LaTeX](main.tex) | [Source bundle](source.zip) | [Metadata](metadata.json) | [Verification](verification.json) | [Proof and formalization scope](VERIFICATION.md)

This is a local revision of the [published edition](published/). The exact Zenodo files and their [original checksums](published/deposited-files.json) are preserved there; the current source and PDF differ from the Zenodo and Hexagon edition. The revision incorporates the former unpublished manuscript and its prior corrections.

Build from this directory:

~~~sh
latexmk -pdf -interaction=nonstopmode -halt-on-error -jobname=generalized-vanishing-two-variables main.tex
~~~

The written argument proves the binary theorem. The retained Lean development covers the ternary counterexample and admissible profile family; the complete binary theorem remains partly formalized. The current verification record distinguishes source and layout checks from mathematical certification, a fresh Lean replay, independent replay and external review.

The exposition now states the main result before the literature discussion, explains support separation through a simple derivative example and identifies the prime-dilated endpoint argument as the difficult step. The abstract and metadata describe the same theorem scope. Required AI-assistance wording is included verbatim in the manuscript.

An isolated internal AI PDF-only reading of the current edition is recorded with its hash and five reader answers in verification.json. Independent mathematical review, including the support-certificate and optimal-cutoff arguments, remains pending. No new public deposit is claimed.

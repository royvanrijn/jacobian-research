# A Filtered Frobenius Proof of the Two-Variable Gaussian Moments Conjecture

An alternative proof of the bivariate support classification and sharp mixed-moment cutoff m > deg Q. The introduction attributes the same classification and cutoff to Wilson and compares Long's stronger theorem with one radial variable and an arbitrary-rank torus.

[Current PDF](gaussian-moments-two-variables.pdf) · [LaTeX](main.tex) · [Source bundle](source.zip) · [Metadata](metadata.json) · [Verification record](verification.json) · [Conditional Lean supplement](SUPPLEMENT.md)

This directory contains the current local revision. The [published edition](published/) preserves the exact files deposited on 23 July 2026 at [Zenodo](https://zenodo.org/records/21534699); its [checksums](published/deposited-files.json) remain intact. The current PDF differs from that deposited edition and has not been deposited as a new version.

## Argument and verification

Circular coordinates separate angular constant terms from the radial factorial functional. A lower Newton face supplies a nonzero lowest coefficient; filtered Frobenius reduction then prevents every pure moment from vanishing unless all rotational weights have one strict sign. This gives the mixed-moment conclusion and its optimal cutoff.

The companion Lean development has no sorry declarations but uses two explicit mathematical axioms: the Duistermaat--van der Kallen theorem and finite-type specialization. The explicit cutoff and its sharpness are written results, outside the exported formal endpoint. This revision does not claim a fresh Lean replay or external review.

The current source combines the existing revised proof with a worked circular-coordinate example, a precise source for finite-field specialization, and an explanation of why factorial normalization precedes reduction. Build and rendered-PDF checks are recorded in [verification.json](verification.json). An isolated internal AI PDF-only reading and targeted wording follow-up are recorded with their hashes in verification.json. Independent mathematical review remains open.

## Build

From this directory:

~~~sh
latexmk -pdf -interaction=nonstopmode -halt-on-error -jobname=gaussian-moments-two-variables main.tex
~~~

The source bundle includes main.tex, references.bib, SUPPLEMENT.md and this build guide. Formal toolchain, endpoint and axiom-report instructions are in the supplement.

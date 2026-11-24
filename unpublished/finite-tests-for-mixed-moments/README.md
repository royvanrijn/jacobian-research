# Finite Tests for Mixed Moments of Regular Functions

**Unpublished companion draft**

[PDF](finite-tests-for-mixed-moments.pdf) · [LaTeX source](main.tex) ·
[Metadata](metadata.json) · [Verification record](verification.json)

The paper gives a finite test for cofinal mixed moment vanishing
against every fixed multiplier under an algebraic integration-by-parts
hypothesis. Independent formal variables retain all offsets from one
diagonal order. An indicial polynomial bounds the number of ordinary
moments needed for the test. For Laurent constant terms a common
linear weight gives an exact criterion and a positive integer ray
detects every failure.

## Build

Run from this directory with TeX Live.

    latexmk -pdf -interaction=nonstopmode -halt-on-error -jobname=finite-tests-for-mixed-moments main.tex

The single LaTeX source includes the bibliography.

## Review

The verification record pins the PDF and records the exact arithmetic,
source inspection and PDF-only internal review. Specialist review of
the effective logarithmic computation remains open. The paper gives a
general decision procedure and a complete worked Laurent example.

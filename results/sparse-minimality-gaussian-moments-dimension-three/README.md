# Sparse Minimality of Gaussian-Moments Counterexamples in Dimension Three

**Current repository revision · 17 November 2026 · not deposited.**

[Current PDF](sparse-minimality-gaussian-moments-dimension-three.pdf) ·
[LaTeX](main.tex) · [Source bundle](source.zip) ·
[Verification](verification.json) · [Metadata](metadata.json)

An exhaustive exact calculation excludes Gaussian-moment failures with **both**
degree at most four and at most four monomials in Gaussian Witt coordinates.
All-order arguments treat the three surviving coefficient families. Together
with Long's example, this establishes five-term minimality among failures of
degree at most four. It does not establish global degree or Pareto minimality,
or exclude higher-degree four-term failures.

The revision explains the central finite-to-all-order argument, proves the
one-dimensional charge enumeration bound directly, expands Gaussian invariance
and shear calculations, and removes the deposited edition's unsupported
product-order minimality inference. It supplies a standalone exact checker
supplement and the repository's standard AI attribution.

The exact [published edition](published/) and its [checksums](published/deposited-files.json)
are preserved separately: [Zenodo record](https://zenodo.org/records/21534917),
deposited 24 July 2026. The current PDF and source differ from that edition.

## Build and replay

Build from this directory, including after extracting the source bundle:

~~~sh
latexmk -pdf -interaction=nonstopmode -halt-on-error -jobname=sparse-minimality-gaussian-moments-dimension-three main.tex
~~~

The source includes its bibliography; no research-worktree file is needed to
compile. Run the exact computation from [the supplement](supplement/README.md)
with Python 3.12, SymPy 1.14.0 and Singular 4.3.2.

The retained research checker was replayed successfully for this revision.
Its two source files are included byte-for-byte in the supplement. The
[verification record](verification.json) distinguishes exact replay,
source-aware proof review, bibliography checks, source-bundle build and PDF
inspection from independent mathematical review.

An isolated internal AI reader read the PDF and rechecked the corrected supplement identifiers and filter cutoff. External mathematical review of the computational coverage and all-order deductions remains open.

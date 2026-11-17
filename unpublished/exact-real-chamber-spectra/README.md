# Exact Real Chamber Spectra of Explicit Keller Maps

**Unpublished · parked.**

For every N >= 3, one explicit rational Keller map realizes every parity-compatible real fiber count on the open target locus C != 0 with squarefree degree-N inverse polynomial. The folds belong to the inverse-incidence projection, while the Keller map has determinant one everywhere.

[PDF](paper.pdf) · [LaTeX](main.tex) · [Metadata](metadata.json) · [Verification](verification.json) · [Source bundle](source.zip) · [Finite checker](verify.py)

The 17 November 2026 revision corrects the transverse fold parameter, gives the polynomial-divisibility and chain-rule Jacobian calculations, and states the open-locus scope in the abstract. It credits Gallagher's all-degree construction and Migus's generic-degree classification; historical priority for the spectrum statement is not claimed as settled.

Build from this directory:

```sh
latexmk -pdf -jobname=paper -interaction=nonstopmode -halt-on-error main.tex
```

With Python 3 and SymPy (tested on Python 3.12.3 and SymPy 1.14.0), run from the repository root:

```sh
.venv/bin/python unpublished/exact-real-chamber-spectra/verify.py
```

The checker verifies rational seeds, endpoint and divisibility identities and rational reconstruction for N=3 through 8. For N=3 through 6 it isolates the discriminant levels and counts real roots exactly over rational sample levels, recovering the full parity spectrum. Its output is [symbolic_checks.json](symbolic_checks.json). These finite controls supplement the written all-degree proof; they do not establish its universal scope or historical novelty.

The repository construction source is pinned at [c32b9b4](https://github.com/royvanrijn/jacobian-research/tree/c32b9b4d273cf7988e85184a0d5eae9c488727fd/unpublished/exact-real-chamber-spectra). The supplied source bundle contains this revised manuscript and the checker; it has not been publicly deposited. Fresh compilation and source-aware internal visual review are recorded in verification.json. An isolated internal AI PDF-only reading and a targeted PDF follow-up are recorded in verification.json. The follow-up confirmed the added finite-normalization argument; external mathematical review remains unperformed.

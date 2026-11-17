# Generic Discriminants of Polynomial Tangent Pencils

**Unpublished · parked.**

A self-contained contact-stratum proof of the classical cusp-node specialization for the reduced discriminant curve of H(W)-sW+t, with explicit genericity for polynomial graphs and smoothness at the unique point at infinity. Priority for the numerical counts is not established.

[PDF](paper.pdf) · [LaTeX](main.tex) · [Metadata](metadata.json) · [Audit](AUDIT.md) · [Verification](verification.json) · [Source bundle](source.zip)

The 17 November 2026 revision explains the thin-family contribution, compares the precise sufficient hypotheses in Yuran's Theorems 1.3 and 1.16, gives the exact cubic example H=W^3, keeps Theorem 1.1 together, and replaces stale audit descriptions of removed relative-infinity and Jacobian-map material.

Build from this directory:

```sh
latexmk -pdf -jobname=paper -interaction=nonstopmode -halt-on-error main.tex
```

The [partial Lean source](https://github.com/royvanrijn/jacobian-research/tree/c32b9b4d273cf7988e85184a0d5eae9c488727fd/research/foundations/formal/discriminant-pencils) is pinned at repository commit c32b9b4d273cf7988e85184a0d5eae9c488727fd, with Lean 4.32.1 and the project's mathlib lockfile. From that project directory, the recorded command is:

```sh
lake build
```

The project covers tangent parametrization and jets, contact algebra, formal infinity calculations and stated affine-exhaustion interfaces. The global dimension argument, full normalization and degree theorem, scheme-theoretic cusp/node criteria and genus-delta step are written arguments outside its scope. No fresh Lean build or transitive axiom audit was run during this revision.

The supplied source bundle contains the revised note and edition-specific audit; it has not been publicly deposited. Fresh LaTeX compilation and source-aware internal visual review are recorded in verification.json. The revised edition has received an isolated internal AI PDF-only reading, recorded in verification.json; external human review remains unperformed. Keep parked until the intended standalone or companion contribution is settled.

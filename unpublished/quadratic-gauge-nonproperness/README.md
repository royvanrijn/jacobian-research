# The Exact Nonproperness Locus of the Quadratic-Gauge Keller Maps

**Unpublished · parked.**

Computes the exact reduced nonproperness locus, complete boundary fiber table and discriminant order for the quadratic-gauge maps.

[PDF](paper.pdf) · [LaTeX](main.tex) · [Metadata](metadata.json) · [Verification](verification.json) · [Source bundle](source.zip) · [Finite checker](verify.py)

The 17 November 2026 revision includes all three map coordinates, the source variables t=1+xy and q, derivative-unit reconstruction, the Jacobian and generic-degree argument, and the identities on t=0 and q=0. The boundary table can therefore be checked without access to the companion. Jelonek's Proposition 6, p.261, was checked in the primary 1993 text; the paper states its characteristic-zero transfer explicitly. The precise [companion edition](https://github.com/royvanrijn/jacobian-research/tree/c32b9b4d273cf7988e85184a0d5eae9c488727fd/unpublished/common-arithmetic-fibers) is cited for context.

Build from this directory:

```sh
latexmk -pdf -jobname=paper -interaction=nonstopmode -halt-on-error main.tex
```

With Python 3 and SymPy (tested on Python 3.12.3 and SymPy 1.14.0), run from the repository root:

```sh
.venv/bin/python unpublished/quadratic-gauge-nonproperness/verify.py
```

The checker verifies the generic cubic map determinant and the boundary identities with symbolic g1,g2,g3. It checks the degree-three discriminant specialization, the symbolic degree-four boundary coefficient, and discriminant orders for G=sum S^k in degrees four through seven at (B,C)=(2,3),(3,2),(1,-1). It records [symbolic_checks.json](symbolic_checks.json). These checks do not verify the all-degree Newton argument, graph-boundary descent or Jelonek criterion; none of those is part of the companion's Lean formalization.

The supplied source bundle contains the revised standalone note and checker. It has not been publicly deposited. Fresh compilation and source-aware internal visual inspection are recorded in verification.json. An isolated internal AI PDF-only reading and a targeted PDF follow-up are recorded in verification.json. The follow-up confirmed the explicit determinant-one output normalization; external specialist review remains unperformed.

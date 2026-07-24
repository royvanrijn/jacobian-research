# Finite Étale Algebras as Keller Fibers

This paper introduces **Keller fibers** and proves a complete rank
classification: every nonzero finite étale algebra of rank other than two
occurs as a full fiber of a polynomial Keller map. For rank at least three,
the realization is explicit in affine three-space, has Jacobian determinant
`1`, and has coordinate degree at most `6N+2` in rank `N`.

The arithmetic applications include:

- an explicit degree-five Keller fiber that is everywhere locally soluble
  over `Q` but has no rational point, with degree five proved optimal;
- one fixed Keller map with infinitely many such Hasse-failing fibers;
- exact transfer of connectedness, signatures, splitting fields, local
  factorization data, and intersectivity;
- compatibility with extension of the ground field.

The directory name is retained as a stable repository path from the earlier
draft.

The four focused audits accompanying the active draft are:

1. the exact arbitrary-characteristic-zero-field scope of the degree-two
   Galois exclusion, with Campbell--Razar--Wright provenance and an explicit
   descent proof;
2. two-sided coordinate-ring reconstruction using a Bézout inverse of `E'`;
3. the scaling identity `F_displayed = diag(1,19,19) F_normalized` for the
   optimal quintic example;
4. the dated and qualified [literature audit](LITERATURE_AUDIT.md).

Run the exact checker from the repository root:

```bash
.venv/bin/python scripts/verify_finite_etale_keller_fibers.py
```

Build the staged Lean certificate with:

```bash
cd formal/finite-etale-keller
lake build
```

Build the paper with:

```bash
latexmk -pdf -interaction=nonstopmode -halt-on-error \
  -outdir=tmp/pdfs/common-arithmetic-fibers \
  papers/common-arithmetic-fibers/main.tex
```

The active-paper build copies the PDF to:

```text
output/pdf/common-arithmetic-fibers.pdf
```

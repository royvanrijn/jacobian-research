# Finite Étale Algebras as Keller Fibers

This paper introduces **Keller fibers** and proves a complete rank
classification: every nonzero finite étale algebra of rank other than two
occurs as a full fiber of a polynomial Keller map. For rank at least three,
the realization is explicit in affine three-space and has determinant `-2`.

The arithmetic applications include:

- an explicit degree-five Keller fiber that is everywhere locally soluble
  over `Q` but has no rational point, with degree five proved optimal;
- one fixed Keller map with infinitely many such Hasse-failing fibers;
- exact transfer of connectedness, signatures, splitting fields, local
  factorization data, and intersectivity.

The directory name is retained as a stable repository path from the earlier
draft.

Run the exact checker from the repository root:

```bash
.venv/bin/python scripts/verify_common_arithmetic_fibers.py
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

# Common Arithmetic Fibers of Stably Inequivalent Keller Maps

This paper proves that two fixed Keller maps over `Q` share infinitely many
complete connected arithmetic fibers while remaining stably inequivalent.
Over `Q(sqrt(-2))`, the same statement holds for three fixed pairwise
inequivalent quartic maps.

Run the exact synthesis checker from the repository root:

```bash
.venv/bin/python scripts/verify_common_arithmetic_fibers.py
```

Build the paper with:

```bash
latexmk -pdf -interaction=nonstopmode -halt-on-error \
  -outdir=tmp/pdfs/common-arithmetic-fibers \
  papers/common-arithmetic-fibers/main.tex
```

The verified stable PDF is copied to:

```text
output/pdf/common-arithmetic-fibers.pdf
```

# Marked-root Keller maps and degreewise stable multiplicity

This is the canonical repository statement and proof of the `tau(N-1)` lower
bound, organized through five separately auditable lemmas.  The foundational map, cubic marked-root
model, exact image theorem, and weighted theorem appear only as motivation; the
construction formulas and boundary arguments used by the theorem are
included in the paper.

Build from this directory with:

```bash
latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
```

The noncanonical proof-audit companion with repository reproduction hooks is
[`../../DEGREEWISE_MULTIPLICITY_AUDIT.md`](../../DEGREEWISE_MULTIPLICITY_AUDIT.md).

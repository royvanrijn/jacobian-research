# Marked-root Keller maps and degreewise stable multiplicity

This is the canonical repository statement and proof of two complementary
degreewise results: the `(N-3)`-dimensional weighted stable-moduli theorem and
the `tau(N-1)` cross-construction lower bound, the latter organized through
five separately auditable lemmas.  The foundational map, cubic marked-root
model, exact image theorem, and weighted theorem appear only as motivation; the
construction formulas and boundary arguments used by the theorem are
included in the paper.

The proof-critical normalization, stabilization, and thick-contact
calculations are isolated as
[`boundary-exhaustion.tex`](boundary-exhaustion.tex),
[`stable-functoriality.tex`](stable-functoriality.tex), and
[`thick-intersection.tex`](thick-intersection.tex).  They are included by the
main paper and are not competing theorem sources.

Build from this directory with:

```bash
latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
```

From the repository root, the focused regressions for this section are:

```bash
python3 scripts/audit_boundary_exhaustion_independent.py
python3 scripts/audit_thick_intersection_local.py
.venv/bin/python scripts/verify_weighted_seed_schema.py
.venv/bin/python scripts/verify_scheme_boundary_all_parameters.py
.venv/bin/python scripts/verify_quartic_nonproperness_paths.py
```

The noncanonical proof-audit companion with repository reproduction hooks is
[`../../DEGREEWISE_MULTIPLICITY_AUDIT.md`](../../DEGREEWISE_MULTIPLICITY_AUDIT.md).

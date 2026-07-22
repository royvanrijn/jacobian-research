# Marked-root Keller maps and degreewise stable multiplicity

This is the canonical repository statement and proof of the
`1+(N-1)tau(N-1)-sigma(N-1)` cross-construction and parameter-branch
lower bound, organized through the weighted and cancellation constructions,
complete boundary exhaustion, stable normalization functoriality, thick
intersection, and cancellation-parameter faithfulness.  The foundational map, cubic marked-root
model, exact image theorem, and weighted theorem appear only as motivation; the
construction formulas and boundary arguments used by the theorem are
included in the paper.

The former decorated-normalization and Hurwitz--LL material has been split
into [`../decorated-discriminant-normalization/main.tex`](../decorated-discriminant-normalization/main.tex)
and [`../hurwitz-ll-rerooting/main.tex`](../hurwitz-ll-rerooting/main.tex).

The proof-critical normalization, stabilization, and thick-contact
calculations are isolated as
[`boundary-exhaustion.tex`](boundary-exhaustion.tex),
[`stable-functoriality.tex`](stable-functoriality.tex),
[`cancellation-parameter-faithfulness.tex`](cancellation-parameter-faithfulness.tex), and
[`thick-intersection.tex`](thick-intersection.tex).  They are included by the
main paper and are not competing theorem sources.

The construction-independent canonical source for stabilization, including
normalization, valuation, intersection, nilradical, relative-differential,
Fitting, and conductor base change, is
[`../../verified/STABLE_NORMALIZATION_FUNCTORIALITY.md`](../../verified/STABLE_NORMALIZATION_FUNCTORIALITY.md).
The TeX section reproduces the part needed by this paper.

Build from this directory with:

```bash
latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
```

From the repository root, the focused regressions for this section are:

```bash
python3 scripts/audit_boundary_exhaustion_independent.py
python3 scripts/audit_thick_intersection_local.py
.venv/bin/python scripts/verify_weighted_seed_schema.py
.venv/bin/python scripts/verify_cancellation_parameter_faithfulness.py
.venv/bin/python scripts/verify_scheme_boundary_all_parameters.py
.venv/bin/python scripts/verify_quartic_nonproperness_paths.py
```

The noncanonical proof-audit companion with repository reproduction hooks is
[`../../DEGREEWISE_MULTIPLICITY_AUDIT.md`](../../DEGREEWISE_MULTIPLICITY_AUDIT.md).

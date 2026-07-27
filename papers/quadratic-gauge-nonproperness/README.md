# The exact nonproperness locus of the quadratic-gauge Keller maps

This geometric companion contains the material removed from the focused
finite-étale-fiber paper:

- graph-boundary base change and descent to Jelonek's complex criterion;
- the exact reduced nonproperness locus;
- the complete fiber table on `Pi = 0`; and
- the global `Pi`-adic discriminant order from the Newton polygon.

These results are logically independent of prescribed-fiber fullness, which
is already proved by equality of the literal fiber rank and geometric degree.
They are not part of the Lean formalization.

Build with:

```bash
latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
```

The exact symbolic audit is:

```bash
.venv/bin/python scripts/verify_quadratic_gauge_nonproperness.py
```

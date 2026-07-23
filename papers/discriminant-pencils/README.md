# Generic discriminants of polynomial tangent pencils

This directory contains the short replacement paper about the reduced
discriminant curve of

\[
H(W)-sW+t.
\]

The paper proves only the generic cusp--node theorem. The former real-chamber
theorem and repository-specific endpoint slice have deliberately been left
out so that the manuscript has one question and one main result.

Build with:

```bash
latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
```

The internal adversarial proof audit remains in `AUDIT.md`.

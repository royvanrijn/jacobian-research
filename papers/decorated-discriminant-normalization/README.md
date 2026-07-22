# Decorated discriminant normalization and stable moduli

This paper contains the scheme-theoretic discriminant-normalization package:
the full Fitting divisor, node pairing and conductor, boundary marks, the
generic degree-`N-2` decorated map, and the `(N-3)`-dimensional stable-moduli
theorem.  Adding the unique unramified affine root sheet makes the normalized
seed open the normalization of the marked decorated image.  The selected
root already extends on the marked admissible-cover stack, including the
normalized-Stein and conductor comparisons at arbitrary collisions.  Its
descent after contraction is also complete: the finite birational closure is
the normal marked coarse compactification and has one point over every DVR
limit.

All stabilization and Fitting/base-change steps cite the construction-independent
[`../../verified/STABLE_NORMALIZATION_FUNCTORIALITY.md`](../../verified/STABLE_NORMALIZATION_FUNCTORIALITY.md).

Build from this directory with:

```bash
latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
```

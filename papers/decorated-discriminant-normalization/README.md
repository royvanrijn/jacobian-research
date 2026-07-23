# Decorated discriminant normalization and stable moduli

The main moduli argument uses the minimal marked Hessian-divisor invariant
`(P^1; div(H''), 0, infinity)`.  Its quotient map has generic degree `N-2`
and image dimension `N-3`; adding the unique unramified affine root sheet
recovers the normalized seed.  Node pairing and conductor are not inputs to
those D1/F2 proofs.  They are retained as the stronger scheme-theoretic D2
decoration, together with the full Fitting divisor and boundary marks.  The selected
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

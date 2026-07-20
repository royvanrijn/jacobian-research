# Weighted-seed theorem schema and bounded scan

This is an exploratory search ledger. The universal inverse, discriminant, and
monodromy statements are now proved in `WEIGHTED_SEED_THEOREM.md`. The scan is
still not theorem-level evidence for exact images, boundary strata, or
finite-field formulas; those require the seed-specific analyses separated
there.

## Universal one-variable model

For an admissible seed `p(w)`, put

\[
H(w)=\int_0^w p(s)\,ds.
\]

Writing a target as `(A,B,C)` and retaining `c_0` for the seed normalization,
the inverse pencil is

\[
E_{A,B,C}(w)=H(w)-BCw+c_0AC^2.
\]

On `C!=0`, a simple root reconstructs one finite source point and

\[
\gamma=-E'(w)/c_0,\qquad x=-c_0C/E'(w).
\]

The repeated-root discriminant therefore has the rational parameterization

\[
s=p(w),\qquad t=wp(w)-H(w),
\]

where `s=BC` and `t=c_0AC^2`. This is the common algebraic spine for image,
nonproperness, monodromy, and finite-field fiber statistics.

The affine `x=0` map is always triangular. If
`kappa=p'(1)/c_0`, its second-coordinate linear coefficient is
`-c_0/(kappa+2)` and its first-coordinate `z` coefficient is
`b(kappa+2)`. Both are nonzero by the construction hypotheses, so every target
on `C=0` has a unique finite boundary preimage.

## Scan executed

`scripts/scan_weighted_seeds.py` scanned 627 distinct admissible seeds through
degree six. It used the canonical primitives

\[
H_d(w)=w^d(1-w),\qquad 2\le d\le6,
\]

and all deformations

\[
p(w)=2w-3w^2+\sum_{j=1}^{d-2}\theta_j
w(1-w)\left(w^j-\frac{6}{(j+2)(j+3)}\right)
\]

with `theta_j` in `{-2,-1,0,1,2}` and nonzero top coefficient. Simple-root
histograms for `H(w)-sw+t` were recorded over `F_5`, `F_7`, and `F_11` when
the seed retained its degree in that characteristic.

Results:

- exactly 5 of 627 candidates had no primitive zeros beyond the distinguished
  roots `0` and `1`; these were precisely the five canonical `H_d` models;
- 624 of 627 had squarefree `p'(w)` in characteristic zero;
- the three failures of squarefree critical behavior were canonical degrees
  `d=4,5,6`, where the high-multiplicity zero of `H_d` forces higher
  ramification at `w=0`; and
- noncanonical degree-three seeds all introduced one additional root of `H`.
  The already-tested seed `p=w-2w^3` gives
  `H=w^2(1-w)(1+w)/2`, exposing the extra `w=-1` boundary branch.

Thus the scan reveals a real tradeoff. The canonical family gives the cleanest
`C=0` geometry, but from degree four onward it has nongeneric higher
ramification. Generic deformations restore simple critical behavior at the
cost of additional `C=0` branches.

The finite-field histograms use the number of simple roots of the inverse
pencil, so their possible fiber sizes line up with fixed-point counts in the
natural permutation action. This is evidence for an `S_n` theorem schema, not
a proof of geometric monodromy or an exact all-`q` distribution.

Reproduce with:

```bash
make scan-weighted-seeds
```

The detailed generated output is `results/weighted_seed_scan.json`.

# The dimension-38 quartic HN frontier

## Two independent constructions

There are now two exact routes to a homogeneous quartic
Hessian-nilpotent counterexample in dimension \(38\).

The public construction is Mohammad Traboulsi's
[`Nested-Tail-components-for-Keller-Maps`](https://github.com/mtraboulsi689/Nested-Tail-components-for-Keller-Maps)
at commit `73635c96034bb8364c036fbec2e366224e601b40`.  Its chain is

\[
3\longrightarrow 11\longrightarrow 18\longrightarrow 19
\longrightarrow 38.
\]

The eleven-variable map is degree three.  After normalization its cubic
output space has rank seven, with active coordinates
\(1,2,3,4,5,6,9\).  Rank-compressed companion cancellation therefore uses
\(11+7=18\) variables, cubic homogenization adds one variable, and the
homogeneous cotangent lift doubles \(19\) to \(38\).

The independent repository construction starts from MacFarlane's
\(F_{13}\).  A source coordinate
\[
s=F_{13,13}=x_{13}+x_2^2
\]
and the target completion \(y_4\mapsto y_4-y_8^2\) produce the explicit
twelve-variable degree-three map in
[the canonical note](../verified/TWELVE_VARIABLE_DEGREE_THREE_KELLER_COUNTEREXAMPLE.md).
Its cubic output rank is six.  Thus its chain has the same compressed cost:

\[
(n,r)=(12,6),\qquad n+r=18,\qquad
12+6+1=19,\qquad 2\cdot19=38.
\]

The constructions are not the same reduction in different notation.  They
give two points on the same \((n,r)\) tradeoff line:

\[
(11,7)\quad\text{and}\quad(12,6).
\]

This changes the immediate search objective.  With the ordinary cotangent
lift, dimension \(37\) cannot occur because the lift doubles dimension.
The next record available to this pipeline is \(36\), and it requires

\[
n+r\le17
\]

before the final homogenizing coordinate is added.

## What the public construction contributes

The public paper's nested-tail companion is more general than the
degree-three rank compression previously used here.  If

\[
\Phi=X+H_2+\cdots+H_D
\]

and \(r_d\) is the output coefficient rank of the tail
\(H_d+\cdots+H_D\), its companion cost is

\[
n+\sum_{d=3}^{D}r_d.
\]

This is the correct score for a graph-deletion route that raises degree.
Variable count alone, terminal degree alone, and the rank of only the
highest homogeneous piece can all select the wrong route.

The public exact frontier also prevents several duplicate searches:

- among its 256 fixed-order graph-coordinate deletions, the unique minimum
  companion cost is \(14\), but it is a three-variable degree-seven map
  with tail ranks \(3,3,2,2,1\); it yields degree-eight generalized
  vanishing examples, not a quartic HN improvement;
- the only degree-three route in that deletion family is the original
  eleven-variable map;
- a six-variable affine feature-graph repair has inconsistent first-order
  determinant systems of ranks \((17,18)\), \((17,17)\), and \((17,18)\);
- after deleting \(h,k\), the standard one-auxiliary cubicizations based
  on \(x^2\) or \(xz\) both have obstruction ranks \((87,88)\);
- deleting the zero homogenizing component from the nineteen-variable
  parent degenerates the ordinary Laplacian and does not give dimension
  \(37\).

These are bounded obstructions to the stated ansatzes, not lower bounds.
In particular, they do not cover the source-coordinate/target-completion
move that produced the repository's \(13\to12\) reduction.

## Cross-construction search

The accompanying checker applies techniques learned from each construction
to the other.

### Completing a pivot before deletion

For each source-affine graph coordinate \(v\) of the public \(G_{11}\), it
tests target coordinates

\[
g(Y)=Y_v+P(Y),
\]

where \(P\) is an arbitrary linear or quadratic polynomial in target
components whose pullbacks are independent of \(v\).  The completion is
shifted to vanish at the common collision image.  Solving
\(g(G_{11}(X))=0\) for \(v\) is then polynomial and all coefficients above
degree three are linear in the coefficients of \(P\).

For the seven source-affine pivots

\[
a,b,c,q,s,h,k,
\]

every exact system is inconsistent: its augmented rank is one larger than
its full column rank.  The coordinate \(d\) occurs quadratically in other
components and is deliberately outside this linear audit.

### The nonlinear \(z_8\) pivot

The local pivot \(K_8=z_8+z_1z_2\) is more subtle because \(K_4\) is
quadratic in \(z_8\).  Let \(P\) be a target polynomial in the eight output
components whose pullbacks are independent of \(z_8\), and complete the
pivot to \(g=Y_8+P\).  On the slice \(g=0\), writing
\(R=P(K_{\mathrm{allowed}})\), the fourth component becomes

\[
K_4=z_4+(z_1z_2)^2-R^2.                                        \tag{10}
\]

If (10) has degree at most three, the highest homogeneous part of \(R^2\)
first forces \(\deg R\le2\).  Unique factorization then forces

\[
R_2=\pm z_1z_2.
\]

For all target completions \(P\) of degree at most two, the exact condition
that \(P(K_{\mathrm{allowed}})\) have degree at most two is a
\(206\times44\) linear system of rank \(39\), with a five-dimensional
kernel.  The quadratic parts of that kernel have rank five, while adjoining
\(z_1z_2\) raises the rank to six.  Thus no degree-at-most-two target
completion of this nonlinear pivot can preserve degree three.

### Coordinated source shears

On the local \(F_{12}\), every single-coordinate quadratic source shear
that preserves degree three is trivial.  This does **not** remain true for
coordinated shears.  Simultaneously shear

\[
(z_4,z_7,z_9,z_{10},z_{11})
\]

by arbitrary quadratic polynomials in

\[
(z_1,z_2,z_3,z_5,z_6,z_8,z_{12}).
\]

These five coordinates form a jointly affine block in the map.  The exact
degree-at-most-three system has size \(146\times140\), rank \(126\), and
therefore a new fourteen-dimensional kernel.

This is a real extension of the current search language: separate
coordinate tests miss all fourteen directions.  It does not yet lower the
record.  On the cubic coefficient matrix, a fixed \(5\times5\) minor has
determinant \(6\).  The Schur-complement conditions for total cubic-output
rank at most five reduce to six affine equations in the fourteen shear
parameters.  Their coefficient and augmented ranks are

\[
(5,6).
\]

Hence no member of this entire coordinated quadratic family lowers the
cubic rank from six to five.

## Current conclusion and next search

No dimension below \(38\) is proved here.  The exact new information is:

1. both independent records saturate the same cost \(n+r=18\);
2. seven natural completed-pivot reductions of the public lift are
   obstructed;
3. the closest nonlinear local pivot has no quadratic target completion;
4. the local map has a previously unused fourteen-parameter family of
   coordinated degree-preserving source automorphisms;
5. that whole family is disjoint from the cubic-rank-at-most-five locus.

The most credible remaining \(36\)-dimensional searches are now:

- nonlinear pivot completion at the public coordinate \(d\), or a
  target-degree-at-least-three completion of the local \(z_8\) pivot;
- coordinated source **and** target automorphisms, rather than a shear on
  only one side;
- an earlier coordinate-pair restriction of the public eleven-variable
  stable lift or the local twelve-variable map outside the graph-coordinate
  basis;
- a non-nested state realization with total tail-rank score at most
  seventeen;
- a symmetric Hessian-nilpotent lift that does not double every cubic
  state.

The first two are the nearest extensions of the successful \(F_{13}\)
square-completion mechanism.  The last one is structurally different and
is the only listed route that could naturally address an odd target such
as \(37\).

## Subsequent closure

The follow-up
[nonlinear-pivot and coordinated left-right audit](HVC38_GAP_CLOSURE.md)
excludes the public \(d\)-pivot and local \(z_8\)-pivot through target degree
eight.  It also derives a seventeen-parameter genuine triangular
source-target family and proves that its exact degree-three locus never has
cubic-output rank below six.  Thus the corresponding open items above
should now be read with those stronger bounds.

## Verification

Run

```bash
.venv/bin/python scripts/audit_hvc38_cross_construction_frontier.py
```

The generated exact record is
[`hvc38_cross_construction_frontier.json`](../artifacts/generated-results/hvc38_cross_construction_frontier.json).
The calculation proves only the finite rank statements described above;
it is not a global minimality claim and has no Lean formalization or
external review.

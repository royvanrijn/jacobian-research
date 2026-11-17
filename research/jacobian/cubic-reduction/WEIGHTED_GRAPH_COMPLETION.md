# Nonliteral weighted carrier graphs: no one-stage cubic completion

KCUB27, 15 November 2026. Let K be the fixed integral twelve-variable cubic
of [KCUB3](PROOF12.md), and Cj=Kj for j=4,...,12. Take any nonzero linear
combination H of carriers of one common weight, normalize its highest
nonzero-index coefficient to1, solve H=0 for that source coordinate, and
delete that coordinate and output. In every resulting eleven-variable raw
map, a specified bad output cannot be made cubic by subtracting any
polynomial of the other ten raw outputs, in any target degree.

Thus none of these complete projective weighted linear graph families
produces an eleven-variable cubic by the specified one-stage elementary
completion. This is not a minimum-dimension theorem, nor an exclusion of
cross-weight/base-output combinations, nonlinear graph coordinates,
coordinated target changes using the bad output, or ordered stages.
No smaller positive witness is found. Global cubic dimension[4,12],
ordinary cubic12/36, cubic-linear89 and all other simplicity bounds stay
unchanged. The full objective remains open.

## All nine projective charts

K is weighted for

```
(1,-1,-2,0,-3,2,-1,-2,-2,-1,-1,0).
```

The carrier weight blocks are

```
weight0: C4,C12;  weight-3: C5;  weight2: C6;
weight-1: C7,C10,C11;  weight-2: C8,C9.
```

Normalize the last nonzero coefficient in each block. This gives exactly
nine charts: literal C4,C5,C6,C7,C8; the two pencils C10+t*C7 and C9+t*C8;
and C11+s*C10+t*C7 and C12+s*C4. All zero-coefficient boundaries occur
among these charts. No generic nonzero parameter stratum is dropped.

Each H is monic affine in its highest source coordinate, since the
carriers are triangular. Source graph straightening and the linear target
replacement of the deleted row are polynomial automorphisms of determinant1.
Their composition has the deleted output equal to the new source coordinate.
Restriction to its zero level is Keller. The two original collision points
lie on every such level because all nine Cj vanish there; their retained
base coordinates are distinct. These facts are independently reconstructed
as polynomial identities. The obstacle is cubic degree after completion,
not determinant or collision inheritance.

The five literal charts are excluded by the retained arbitrary-target-degree
[KCUB4 certificates](GRAPH_DELETION_OBSTRUCTION.md).

## Two charts inherit the old obstruction exactly

Only K3 and C11 involve source x11. On H=C11+s*C10+t*C7=0, all ten outputs
other than the deleted row and K3 are unchanged from the literal C11 graph.
Writing G for the restriction gives the exact identity

```
G3(s,t)=G3(0,0)+4*s*C7*C10+4*t*C7^2.
```

The added terms are polynomials in unchanged retained outputs. Any polynomial
completion of G3(s,t) would therefore give a completion of the literal
G3(0,0), contradicting KCUB4.

Similarly, only K2 and C12 involve source x12. On H=C12+s*C4=0,

```
G2(s)=G2(0)+12*s*C7*C4,
```

and the other ten outputs are unchanged. KCUB4 again transfers in every
polynomial target degree. These statements concern the fixed raw output
algebra; they do not exclude first mixing the bad output into other rows.

## The two genuine pencils

The remaining families are

| Deleted output | Source graph | Bad retained row |
| --- | --- | ---: |
|C10+t*C7|x10=-x1*x8-t*(x7+x1*x3)|2|
|C9+t*C8|x9=-2*x1*x5+3*x2^2-t*(x8+x2^2)|3|

In both, the other ten raw outputs are cubic for every t, and the retained
map G has zero constant part and identity linear part. We prove a stronger
statement for these two pencils: even a rational function of the other
ten outputs cannot repair the bad row to a polynomial cubic.

Let i be the bad row and let Delta be the inverse-Jacobian column derivation
of the eleven-variable G, so Delta(G_j)=delta_ij. If
f=G_i-q(G_j:j!=i) is a polynomial cubic, then Delta f=1, including rational q.
This is a necessary condition, not a sufficiency claim about target membership.
Since the other ten G_j are cubic with distinct identity linear jets,
subtracting a constant and their scalar multiples normalizes f to

```
f=X_i+Q2+Q3,
```

with Q2,Q3 homogeneous of degrees2,3 in the eleven retained source variables.
There are66+286=352 possible monomials. The bad linear coefficient is1
because Delta at0 is partial_i. Consequently a repair requires

```
Delta(Q2+Q3)=1-Delta(X_i).
```

Only source coefficients through order7 are needed to obstruct this
finite-dimensional source equation. Supplied order7 polynomial jets D of
Delta satisfy D(G_j)=delta_ij modulo source order>=8. Because JG(0)=I,
the inverse-Jacobian jet is unique, so these are the actual order7 jets,
for every t. The independent checker verifies those identities directly;
it does not trust the primary inverse-substitution algorithm.

## Polynomial dual identities, valid at every parameter

For each pencil there is a coefficient functional

```
Lambda_t(P)=sum_e w_e(t)*[X^e]P,
```

where all w_e are rational polynomials in t of degree<=2, and all sampled
source monomials have degree<=7. Exact reconstruction gives

```
Lambda_t(D(Q))=0  for every one of the352 quadratic/cubic basis monomials Q,
Lambda_t(1-D(X_i))=1.
```

| Pencil | Source coefficient weights | Nonzero parameter coefficient terms | Parameter degree |
| --- | ---: | ---: | ---: |
|C10+t*C7|42|60|2|
|C9+t*C8|59|84|2|

Thus applying Lambda_t to the necessary equation gives0=1. The identities
are polynomial in t, with no parameter denominator and no omitted roots.
They hold over every characteristic-zero field and its extensions. This
proves the two pencil exclusions in arbitrary polynomial/rational target
degree. The order7 truncation is a finite test of the desired degree<=3 source
equation; it places no degree bound on the unknown target q.

## Independent replay and retained discovery limits

Run

```
.venv/bin/python research/jacobian/cubic-reduction/verify_graph_pencils12.py
```

The checker uses standard-library Fraction arithmetic. It independently
reconstructs both source graphs and all retained outputs, verifies graph
straightening, generic parameter collisions and identity jets, checks the
supplied inverse-Jacobian jets on every raw output, and multiplies all704
basis annihilations plus the two unit pairings exactly. It also reconstructs
the two target-product transfer identities and all nine projective chart
labels. The old literal exclusions are an explicit pinned KCUB4 premise;
they are not silently re-proved by the new pencil check.

- [New independent receipt](graph-pencils12/validation.json)
- [Pencil10 order7 model](graph-pencils12/pivot-10-D7/model.json)
- [Pencil10 uniform dual](graph-pencils12/pivot-10-D7/dual-adaptive.json)
- [Pencil9 order7 model](graph-pencils12/pivot-09-D7/model.json)
- [Pencil9 uniform dual](graph-pencils12/pivot-09-D7/dual.json)

The order5 system was insufficient for pencil9; its apparent cubic-slice
candidates were necessary-system false positives, not witnesses. Increasing
to order7 removes them. For pencil10, the first restricted-support dual had
pairing t-3/8, so it did not cover t=3/8. The exact system at that point
remained inconsistent. Adding its coefficient support produced the uniform
unit pairing above. The preliminary bounded dual misses and gcd receipt
are retained; no generic denominator was discarded and no finite sampling
was promoted to a theorem. The universal conclusion depends on the final
polynomial identities, not on the discovery probes or their caps.

The older BCR11/12 full linear-graph theorem concerns the different six-cubic
BCR2 map. Its cubic-output rank is6, whereas the present map has rank5;
that older result was not reused as a nonliteral-graph exclusion here.

## Remaining construction routes

This closes same-weight carrier graph combinations with one-stage repair
of a bad raw output. It leaves different weights mixed together, combinations
involving base outputs, nonlinear graph coordinates, source changes, higher
or ordered target stages that first alter the raw output algebra, and
other constructions. It proves no global eleven-dimensional exclusion,
no full cubic target-completion theorem, and no sparsity or special-form
optimality beyond its stated family. The strict3D census stays62 excluded
and12 remaining, and the global3D degree interval stays[5,7].

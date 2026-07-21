# Boundary geometry and stable distinction

This document contains the third cancellation result: the canonical boundary
object and its application separating every noncubic cancellation map from a
generic weighted seed of the same inverse degree.  It also records the exact
hypotheses needed to make the object canonical.

## 1. Compactification setup

Let `F:U->Y` be a dominant quasi-finite morphism of normal affine varieties,
with `Y` normal.  Write `K=k(Y)` and `L=k(U)`, and let

\[
 \bar X_F=\operatorname{Norm}_Y(L),\qquad
 j_F:U\longrightarrow\bar X_F,
 \qquad \partial_F=(\bar X_F\setminus U)_{\mathrm{red}}.
\]

The map `j_F` is an open immersion: it is the factorization supplied by the
normalization, is birational onto its image, and is quasi-finite; the normal
form of Zariski's Main Theorem identifies it with an open subvariety.  This
point is essential—the boundary is taken in the canonical finite
normalization, not in a chosen projective compactification.

For a boundary prime `E`, retain its image `Z=\overline{F(E)}`, the valuation
extension, ramification index `e(E/Z)`, residue degree `f(E/Z)`, and sheet
loss.  When every `E` dominates a divisor, these data form a labelled
bipartite diagram

\[
 \{E\subset\partial_F\}\longrightarrow\{Z\subset Y\}.
\]

The divisorial-image hypothesis must not be hidden.  In the general
quasi-finite setup a boundary divisor may map into codimension at least two.
The divisorial invariant therefore either assumes boundary noncontraction or
is understood as the codimension-one truncation of a larger stratified
object retaining all images.  Keller maps between affine spaces satisfy the
needed purity/noncontraction statement in the applications here.

## 2. Complete boundary list

Completeness is certified locally.  For a target prime divisor `Z` and a
finite cover of generic degree `nu`, all primes over the generic point obey

\[
 \sum_{E\mid Z}e(E/Z)f(E/Z)=\nu.                           \tag{1}
\]

Thus a proposed list is exhaustive only after it includes both affine and
boundary primes and reaches the full sum (1).  Chart visibility alone is not
a proof.

For a generic weighted seed of degree `N`, all divisorial boundary images are
the discriminant `Delta_H` and the second divisor `C=0`.  For the cancellation
map with

\[
 N=r(m+1)+1,
\]

all divisorial boundary images are the discriminant `Delta_{m,r}` and `P=0`.
Over the discriminant there is one boundary prime with `(e,f)=(r+1,1)`.
Over `P=0` there are `mr-1` geometric unramified boundary primes with
`(e,f)=(1,1)`; the remaining cluster is the affine divisor `B=0`.  Consequently
`P=0` is a boundary-image divisor exactly when `(m,r)!=(1,1)`.

The chart confinement and both degree sums prove that these are all
divisorial components of `partial_F` and all their target images.  In
particular, no unseen component can duplicate a marking, alter the incidence
diagram, or create an extra automorphism.

## 3. Normalization and divisorial images

The generic inverse is controlled by

\[
 \Psi(T)=C\int_0^T\{1-t(Q-Pt)^m\}^r\,dt-R.
\]

Its critical divisor is

\[
 1-T(Q-PT)^m=0.
\]

Putting `Y=Q-PT` normalizes that divisor by

\[
 T=Y^{-m},\qquad P=(Q-Y)Y^m.
\]

It maps birationally to a reduced prime discriminant divisor.  The exponent
`r` gives intrinsic ramification index `r+1`; it does not thicken the
discriminant generator.  The discriminant vertex is therefore intrinsically
marked as the unique target divisor receiving ramification, while `P=0` is
the unique vertex receiving only unramified boundary primes in every
noncubic case.

## 4. Scheme-theoretic intersections

Eliminating the normalized critical coordinate and restricting the reduced
discriminant to `P=0` gives, up to a nonzero scalar,

\[
 Q^{mr(m+1)}\bigl((r+1)RQ^m-C\bigr).                       \tag{2}
\]

Hence its scheme-theoretic intersection ring is

\[
 k[Q,R]/\left(Q^{mr(m+1)}((r+1)RQ^m-C)\right)
 \simeq k[Q,R]/(Q^{mr(m+1)})\times k[Q,Q^{-1}].            \tag{3}
\]

The nilradical has exact index `mr(m+1)`.  Its reduction is

\[
 k[R]\times k[Q,Q^{-1}],                                  \tag{4}
\]

so the reduced intersection is `A^1 disjoint-union G_m`.

Reduced and scheme-theoretic intersections must be distinguished.  The
reduced object records incidence topology; (3) additionally records contact
multiplicity and is strictly stronger.  For `(m,r)=(1,1)`, formulas (2)--(4)
remain true as a coordinate trace, but `P=0` is not a boundary vertex, so no
intrinsic boundary claim is made from that trace.

## 5. The invariant ladder and functoriality

There is a natural hierarchy

\[
 \mathcal I^{\mathrm{formal}}(F)
 \longrightarrow \mathcal I^{\mathrm{sch}}(F)
 \longrightarrow \mathcal I^{\mathrm{red}}(F).            \tag{5}
\]

- `I^red` retains the upstairs and downstairs primes, their labels, and all
  reduced multiple intersections.
- `I^sch` retains the corresponding scheme-theoretic multiple
  intersections and multiplicities.
- `I^formal` retains completed local maps along strata, conductor ideals,
  differents, and valuation filtrations.

The arrows in (5) forget structure.  Conductor data belongs to the formal
layer and must not be claimed to follow from a reduced incidence diagram.
The upstairs vertices are necessary: target diagrams alone forget multiple
primes above one divisor, residue-field covers, upstairs intersections, and
specialization relations.

A polynomial left--right equivalence induces an isomorphism of function-field
extensions and, by functoriality of normalization, an isomorphism of finite
covers carrying the distinguished affine open to the distinguished affine
open.  It therefore identifies all three layers whenever the indicated
formal data are defined.  Stabilization by identity variables is flat base
change: it tensors intersection rings with a polynomial ring, preserves
reduction, nilpotency indices, valuation labels, and the bipartite incidence
diagram.  These statements require the canonical normalization boundary;
they are not assertions about arbitrary compactifications.

## 6. Stable comparison with weighted seeds

For a generic weighted seed, the saturated intersection of its two boundary
images has coordinate ring `k[B]`, hence is reduced.  For every noncubic
cancellation pair, (3) is nonreduced.  Therefore:

### Boundary distinction theorem

For `(m,r)!=(1,1)`, the cancellation map is not polynomially left--right
equivalent, even after stabilization, to a generic weighted seed of inverse
degree `r(m+1)+1`.

The theorem uses the complete boundary list, not merely the components visible
in a resolvent chart.  At the reduced layer one also sees
`A^1 disjoint-union G_m` versus `A^1`; the thick layer supplies the uniform
nilpotent obstruction.

## 7. Independent audit notes

The boundary exhaustion was independently rederived from the generic inverse,
the normalized critical divisor, both degenerate charts, squarefreeness, and
the local degree sums.  The clean-room checker uses only Python's standard
library and is run as

```bash
python3 scripts/audit_boundary_exhaustion_independent.py
```

The scheme-intersection formula is checked in bounded symbolic ranges by
`scripts/verify_c24_scheme_boundary_all_parameters.py`; its all-parameter
proof is the monic-resultant and beta-integral calculation summarized in
(2)--(3).  A genuinely independent audit of that all-parameter thickening is
still desirable.

Detailed former component notes are retained in
[the cancellation archive](../../archive/cancellation-components/).

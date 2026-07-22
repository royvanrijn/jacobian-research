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

The divisorial-image point is a theorem here, not a hidden assumption.  The
finite normalization is affine.  If a component of the complement of the
dense affine open `U` had codimension at least two, normality and extension
across codimension two would identify the coordinate rings on a principal
neighborhood and its intersection with `U`, contradicting the missing
component.  Thus `partial_F` is pure of codimension one.  If `E` corresponds
to a height-one prime `q` in the normalization and `p=q cap k[Y]`, then
`p!=0`; going down for the integral extension over normal `k[Y]` gives
`ht(p)<=1`, hence `ht(p)=1`.  Therefore every boundary prime maps onto a
target prime divisor.  Boundary divisors can contract in arbitrary projective
compactifications or blowups, but not in this canonical affine finite-
normalization setup.  Without affineness, normality, or finiteness, one must
instead impose noncontraction or retain all higher-codimension images in a
stratified object.

## 2. Complete boundary list

This section is a summary.  The canonical paper contains the complete local
proof, including the four generic DVRs, explicit uniformizers, every affine
and boundary prime, geometric and arithmetic residue degrees, reduced prime
discriminant equations, and the no-residual-component saturation argument;
see [the extracted boundary-exhaustion section](../papers/marked-root-multiplicity/boundary-exhaustion.tex).
Its exact weighted-discriminant lemma proves for the displayed degreewise
seed that `Disc_W(E_N)=lambda_N C^2 delta_N`, with `delta_N` prime, reduced,
and nonzero modulo `C`.

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
diagram.  These assertions are instances of the construction-independent
[stable normalization functoriality theorem](../verified/STABLE_NORMALIZATION_FUNCTORIALITY.md),
including its Fitting and conductor base-change clauses.  They require the
canonical normalization boundary;
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

## 7. Degreewise consequence of the boundary calculation

The boundary-exhaustion and contact calculations above supply the boundary inputs to the degreewise result: cancellation type `(m,r)` has `e_Delta=r+1` and `mu=mr(m+1)`, whereas the split weighted contact is reduced.  The marked reconstruction open additionally separates all `mr` parameter roots within a type.  Combined with the two construction lemmas, this yields at least `1+(N-1)tau(N-1)-sigma(N-1)` stable classes in every generic degree `N>=4`.

This section is only a component summary.  The canonical statement and proof are in [Marked-Root Keller Maps and Degreewise Stable Multiplicity](../papers/marked-root-multiplicity/main.tex); the [four-input audit](../DEGREEWISE_MULTIPLICITY_AUDIT.md) records the independent proof obligations.

## 8. Independent audit notes

The boundary exhaustion was independently rederived from the generic inverse,
the normalized critical divisor, both degenerate charts, squarefreeness, and
the local degree sums.  The clean-room checker uses only Python's standard
library and is run as

```bash
python3 scripts/audit_boundary_exhaustion_independent.py
```

The scheme-intersection formula now has two independent all-parameter proofs
in the canonical paper's
[thick-intersection section](../papers/marked-root-multiplicity/thick-intersection.tex):
a primitive monic-resultant calculation and a completed-local-ring length
calculation.  The former is checked in bounded symbolic ranges by
`scripts/verify_scheme_boundary_all_parameters.py`; the latter has a
standard-library exact-arithmetic audit in
`scripts/audit_thick_intersection_local.py`.

The canonical proof and the independent audit serve different purposes.  The
former is the source of the theorem and its local-ring hypotheses; the latter
is a clean-room regression against omitted branches and transcription errors.

Detailed former component notes are retained in
[the cancellation archive](../archive/cancellation-components/).

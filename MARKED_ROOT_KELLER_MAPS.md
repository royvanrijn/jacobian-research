# Marked-root Keller maps

The weighted and cancellation constructions share one inverse-geometric
framework.  Their source formulas differ, but both affine sources are
regular-reconstruction opens inside finite normalized covers obtained by
marking one inverse root.

## A unifying thesis

The repository's organizing mechanism is developed in the standalone
[unifying thesis](UNIFYING_THESIS.md).  In brief, a ramified marked-root
incidence is suspended so that its core Jacobian zero is cancelled by a
vertical or rational source-chart factor.  The resulting polynomial map is
etale, but affine source space is only the regular-reconstruction open in
the finite normalized root cover.  A boundary divisor is omitted exactly
when some source reconstruction coordinate has negative valuation there.

The two principal realizations are:

- a polynomial weighted suspension of the ramified tangent incidence;
- a birational cancellation suspension whose source chart contributes the
  reciprocal boundary power.

The foundational cubic is the smallest point at which the two realizations
are polynomially left--right equivalent.  The standalone note proves the
boundary--reconstruction criterion, identifies the nonproperness locus with
the finite image of the normalization boundary, works out both cubic
ledgers, and separates direct consequences from family-specific theorems and
later analogies.

The deepest recurring dichotomy is

\[
\boxed{\text{formal/local solvability}
\quad\text{versus}\quad
\text{global polynomial algebraization}.}
\]

The thesis is proved for the displayed weighted and cancellation families.
Its use as an exhaustive description is the scoped
[minimal-boundary classification conjecture](cancellation/MINIMAL_BOUNDARY_CLASSIFICATION.md).

## 1. Marked-root presentations

Let `Y` be a normal affine variety over a characteristic-zero field `k`.  A
degree-`N` marked-root presentation consists of:

1. a separable irreducible equation
   \[
   \Psi(T;y)=0\quad\text{over }k(Y)
   \]
   of degree `N` in `T`;
2. the normalization `Xbar_Psi -> Y` of `Y` in
   `k(Y)[T]/(Psi)`;
3. rational reconstruction functions
   \[
   X_1(T,y),\ldots,X_d(T,y)
   \]
   on `Xbar_Psi`;
4. the regular-reconstruction open `U_Psi`, where these functions and their
   inverse identities are regular;
5. an isomorphism `A^d -> U_Psi` under which the finite projection restricts
   to a polynomial map `F:A^d->Y`; and
6. an argument that `det DF` is a nonzero constant.

This definition is intrinsic to the finite function-field extension.  A
displayed affine polynomial in `T` may lose degree on a target divisor, and a
different primitive element may introduce spurious discriminant factors.
The finite normalization and its distinguished regular-reconstruction open,
not the raw resolvent presentation, are the canonical objects.

On the generic etale locus, the `N` roots of `Psi` correspond exactly to the
`N` source points.  Critical root strata describe colliding sheets.  When a
critical root makes a reconstruction coordinate singular, the corresponding
valuation lies on

\[
\partial_F=(\bar X_\Psi\setminus U_\Psi)_{\mathrm{red}}.
\]

The resulting boundary primes, their ramification and residue degrees, and
their upstairs and downstairs intersections form the invariant developed in
[BOUNDARY_GEOMETRY.md](cancellation/BOUNDARY_GEOMETRY.md).

## 2. The two principal realizations

### Weighted marked roots

For an admissible degree-`N` seed `H`, the inverse equation is

\[
\Psi_H(W;A,B,C)=H(W)-BCW+cAC^2.
\]

The inverse-pencil, Jacobian, discriminant-normalization, reconstruction-pole,
and Hessian-Fitting assertions for this equation are all consequences of the
[tangent-map core theorem](verified/TANGENT_MAP_CORE.md).

The generic cover has monodromy `S_N`.  Its discriminant has simple fold
inertia, so its canonical ramified boundary prime has index two.  The global
normalized incidence and reconstruction charts are summarized in the
[weighted theorem](verified/WEIGHTED_SEED_THEOREM.md).

### Cancellation marked roots

For parameters `m,r>=1`, the inverse equation is

\[
\Psi_{m,r}(T;P,Q,R)
=C\int_0^T\{1-t(Q-Pt)^m\}^r\,dt-R.
\]

Its degree is `r(m+1)+1`, and its derivative is

\[
\partial_T\Psi_{m,r}=C\{1-T(Q-PT)^m\}^r.
\]

Thus the distinguished critical divisor has intrinsic ramification index
`r+1`.  Exact reconstruction and polynomial cancellation are proved in
[the cancellation construction](cancellation/CONSTRUCTION.md).

The foundational Keller map is the minimal cubic marked-root map.  It is the degree-three point at
which the weighted and cancellation descriptions coincide up to linear
reparametrization.

## 3. Degreewise stable-multiplicity consequence

The marked-root framework explains the sharp repository lower bound.  For a
fixed proper divisor `r|(N-1)`, put `m=(N-1)/r-1`.  The reconstruction open
recovers the selected cancellation root, so that type contributes all
`mr=N-1-r` parameter-root classes, not merely one representative.  Summing
over the proper divisors and adding one split weighted class gives

\[
 \boxed{1+(N-1)\tau(N-1)-\sigma(N-1)}
\]

pairwise stably inequivalent degree-`N` maps for every `N>=4`.  Distinct
cancellation types are separated by `(e_Delta,mu)`, parameter roots within a
type by reconstruction-open faithfulness, and the weighted class by its
reduced boundary contact.

This is a **finite degreewise count** of explicitly separated classes.  It is
not the dimension of a moduli space.  The separate weighted-seed theorem
gives an `(N-3)`-dimensional family of stable classes; neither statement
subsumes the other because the finite count also distinguishes the
cancellation branches from the weighted construction.

The canonical statement and proof are the standalone paper [Marked-Root Keller Maps and Degreewise Stable Multiplicity](papers/marked-root-multiplicity/main.tex).  The [five-lemma audit](DEGREEWISE_MULTIPLICITY_AUDIT.md) is its verification companion.  This document is canonical for the marked-root framework, not for the degreewise theorem.

## 4. The degree-three `A_2` picture

For a cubic, ordering the three roots gives the type-`A_2` reflection
arrangement: its three Weyl walls are the pairwise diagonals.  Marking one
root distinguishes the two walls on which the marked root collides from the
wall on which only the unmarked pair collides.  The first two form the
repeated-marked-root divisor removed from the affine source, and their
valuation descends to the dicritical normalization boundary.  The unmarked
collision wall remains in the regular-reconstruction open and supplies the
simple marked point over a double-root cubic.

All three walls map to the cubic discriminant after forgetting the ordering,
but the marking determines which sheets stay affine and which escape.  This
is the degree-three model of the general mechanism: collision walls give
discriminant inertia and become canonical normalization boundary precisely
when reconstruction ceases to be regular.

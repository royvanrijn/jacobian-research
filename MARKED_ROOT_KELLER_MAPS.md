# Marked-root Keller maps

The weighted and cancellation constructions share one inverse-geometric
framework.  Their source formulas differ, but both affine sources are
regular-reconstruction opens inside finite normalized covers obtained by
marking one inverse root.

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

The marked-root framework explains the degreewise lower bound: for every `N>=4`, one split weighted map and one cancellation type for each proper divisor of `N-1` give at least `tau(N-1)` stable left--right classes.  Distinct cancellation types are separated by `(e_Delta,mu)`, while reducedness separates the weighted class.

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

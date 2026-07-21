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
[BOUNDARY_GEOMETRY.md](experimental/cancellation/BOUNDARY_GEOMETRY.md).

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
[the cancellation construction](experimental/cancellation/CONSTRUCTION.md).

C01 is the minimal cubic marked-root map.  It is the degree-three point at
which the weighted and cancellation descriptions coincide up to linear
reparametrization.

## 3. Degreewise divisor-count theorem

### Theorem

For every integer `N>=4`, there are at least

\[
\boxed{\tau(N-1)}
\]

pairwise stably inequivalent polynomial maps `C^3->C^3` with nonzero constant
Jacobian, nontrivial fibers, and generic fiber degree `N`.  Here `tau` is the
positive-divisor counting function.

For the weighted map one may take the split admissible seed

\[
H_N(W)=W^2(1-W)(1+W^{N-3}).                                \tag{1}
\]

Its extra roots are the distinct roots of `1+W^{N-3}`; they are nonzero and
different from `1`.  Moreover,

\[
H_N'(1)=-2,
\qquad
\frac{H_N''(1)}{2}=-(N+1)\ne-2,                            \tag{2}
\]

so (1) is admissible.  The weighted theorem gives constant nonzero Jacobian
after taking, for example, `b_0=1` and the forced
`a_0=-(1+kappa)/(2+kappa)`.  It has generic degree `N`; its generic etale
fiber therefore contains `N` distinct points.

Put `n=N-1`.  Every proper positive divisor `r` of `n` determines

\[
m=\frac nr-1\ge1,
\qquad
r(m+1)+1=N.                                                \tag{3}
\]

Conversely, every positive pair satisfying the degree equation arises this
way.  There are therefore `tau(n)-1` such cancellation types.

For each pair, the cancellation parameter polynomial is monic of degree
`mr>=1`.  It has a root over `C`, and uniform separability makes every root
simple, so the finite recurrence constructs a valid polynomial cancellation
jet.  The resulting map has nonzero constant Jacobian, generic degree `N`,
and an explicit `N`-point fiber.

Because `N>=4`, none of these pairs is the exceptional cubic `(m,r)=(1,1)`.
The uniform boundary-exhaustion theorem therefore applies to every pair.
The complete normalization boundary has two intrinsically different target
vertices.  The discriminant vertex is the unique one receiving a ramified
boundary prime, of index

\[
e_\Delta=r+1,                                               \tag{4}
\]

while the `P=0` vertex receives only unramified boundary primes.  They cannot
be exchanged by a left--right equivalence.  Their scheme-theoretic
intersection has nilradical of exact index

\[
\mu_{m,r}=mr(m+1)=m(N-1)
=(N-1)\left(\frac{N-1}{r}-1\right).                        \tag{5}
\]

Distinct proper divisors `r` give distinct ramification indices (4), and
equivalently distinct values of (5).  Functoriality of the finite
normalization preserves these labelled data under arbitrary polynomial
left--right equivalence; stabilization tensors the intersection rings with a
polynomial ring and preserves both ramification and nilpotency index.  Hence
the `tau(N-1)-1` cancellation types are pairwise stably inequivalent.

The weighted map contributes one further class: its corresponding canonical
boundary intersection is reduced, whereas every cancellation intersection
above is nonreduced.  This gives

\[
(\tau(N-1)-1)+1=\tau(N-1)
\]

classes.  The former degreewise doubling theorem is the weakest case of this
divisor-count result.  The boundary proof is also recorded in
[BOUNDARY_GEOMETRY.md](experimental/cancellation/BOUNDARY_GEOMETRY.md).

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

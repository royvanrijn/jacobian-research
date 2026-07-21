# Boundary-exhaustion certificates

This note supplies the reusable completeness criterion needed by intrinsic
boundary signatures.  A resolvent chart may exhibit several affine and
boundary branches, but displaying them does not by itself prove that no
other branch exists.  Exhaustion follows from two independent checks:

1. a support argument confines every divisorial boundary image to a finite
   list of target divisors; and
2. a DVR degree count proves that every branch over each listed divisor has
   been found.

Use the admissible boundary setup of
[the canonical incidence theorem](CANONICAL_BOUNDARY_INTERSECTION_INVARIANT.md):
`X` and `Y` are normal integral affine varieties, `F:X -> Y` is dominant and
quasi-finite, and

\[
 X\hookrightarrow\overline X_F\xrightarrow{\pi_F}Y
\]

is the canonical open immersion followed by the finite normalization.  Put
`K=k(Y)`, `L=k(X)`, and `nu=[L:K]`.

## 1. The local degree certificate

Let `Z` be a prime divisor of `Y`.  Since `Y` is normal, the local ring
`R=O_(Y,eta_Z)` is a DVR.  Let `S` be its integral closure in `L`.  The height-
one primes `E_1,...,E_t` of `bar(X)_F` over `Z` are the maximal ideals of the
semilocal Dedekind ring `S`, and

\[
 \boxed{\qquad
 \nu=\sum_{i=1}^t e(E_i/Z)f(E_i/Z).
 \qquad}                                                       \tag{1}
\]

### Lemma 1.1 (local exhaustion)

Suppose distinct primes `E_1,...,E_s` over `Z` have been constructed and
their indices and residue degrees satisfy

\[
 \sum_{i=1}^s e(E_i/Z)f(E_i/Z)=\nu.                            \tag{2}
\]

Then `s=t`: there is no further prime over `Z`.

**Proof.**  The finite torsion-free `R`-module `S` is free of rank `nu`.
If `pi` is a uniformizer of `R`, then `dim_(kappa(Z))(S/pi S)=nu`.  Localizing
the Artinian ring `S/pi S` at its maximal ideals decomposes this dimension
into the positive contributions `e(E_i/Z)f(E_i/Z)`, proving (1).  Equality
(2) leaves no positive contribution for another prime.  QED

The sum must include both boundary primes and primes whose generic points lie
in the affine-source open.  Counting only the missing branches generally
does not reach `nu` and cannot prove exhaustion.

## 2. The global certificate

### Theorem 2.1 (support plus degree accounting)

Let `Z_1,...,Z_s` be distinct target prime divisors.  Assume:

1. every divisorial component of `partial_F` maps into one of the `Z_i`
   (the stronger statement that the canonical open immersion is surjective
   over `Y minus union_i Z_i` is a convenient sufficient check); and
2. for each `Z_i`, a list of distinct affine and boundary primes over `Z_i`
   has been constructed whose `sum ef` is `nu`.

Then the displayed boundary primes are all divisorial components of
`partial_F`, their images are exactly the members of the candidate list which
receive a boundary prime, and their ramification profiles are exhaustive.

**Proof.**  The support hypothesis rules out a divisorial boundary image away
from the candidate list.  Lemma 1.1 proves completeness of the prime list
over every candidate divisor.  Classifying the generic point of each listed
prime as lying inside or outside the distinguished open `X` then gives the
complete boundary list and profile.  QED

A convenient way to prove the first hypothesis is to give regular
reconstruction formulas on the complement of the candidate divisors and to
show that all `nu` inverse branches reconstruct there.  A denominator list
alone is insufficient: apparent poles may be removable chart poles or finite
source divisors.

## 3. Generic weighted marked-root theorem seeds

Let the inverse degree be `n`, and let

\[
 H(W)=hW^2(W-1)\prod_{j=1}^{n-3}(W-\rho_j)
\]

with distinct nonzero `rho_j`.  The inverse reconstruction and the two root
charts confine divisorial boundary images to the discriminant `Delta_H` and
`C=0`: when `C` and the discriminant are nonzero, all `n` inverse roots are
simple and every reconstruction denominator is a unit, so every normalized
branch lies in the affine-source open.

At the generic point of `Delta_H`, one boundary prime contributes `e f=2`
and the other `n-2` inverse branches are affine and simple.  The ramified
prime cannot lie in `X`, because the Keller map is étale there.  Their total
is

\[
 2+(n-2)=n.                                                   \tag{3}
\]

At the generic point of `C=0`, the double-zero cluster contributes two affine
sheets through one étale affine prime of `(e,f)=(1,2)`, the distinguished root
`W=1` contributes one affine `(1,1)` prime, and the `n-3` extra simple roots
give boundary primes.  Thus

\[
 2+1+(n-3)=n.                                                 \tag{4}
\]

Lemma 1.1 proves that both lists are exhaustive.  In particular, the complete
boundary profiles are one `(e,f)=(2,1)` prime over `Delta_H` and `n-3`
geometric `(1,1)` primes over `C=0`.

## 4. Cancellation construction

For `cancellation type (m,r)`, put `N=r(m+1)+1`.  The affine and projective resolvent charts
and the reconstruction formulas confine divisorial boundary images to the
discriminant `Delta_(m,r)` and `P=0`: off these divisors, every inverse root
is noncritical, hence has `D!=0`, and the reconstruction formulas recover all
`N` roots inside the affine-source open.

At the generic point of the discriminant, the integral critical divisor gives
one boundary prime with `e f=r+1`; the other `N-(r+1)` sheets are affine.
Again the ramified prime cannot lie in the étale affine-source open.  Hence

\[
 (r+1)+N-(r+1)=N.                                            \tag{5}
\]

At `P=0`, the `U=0` cluster is the affine divisor `B=0`.  Étaleness gives
`e=1`, while its generic residue degree is `r+1`, so it contributes `r+1`
sheets.  The distinguished nonzero root is the affine divisor `A=0` with
`(e,f)=(1,1)`.  The remaining `mr-1` simple roots are boundary primes.  The
exact count is

\[
 (r+1)+1+(mr-1)=r(m+1)+1=N.                                  \tag{6}
\]

Lemma 1.1 therefore proves the projective two-chart list exhaustive.  The
complete profiles are one `(r+1,1)` boundary prime over the discriminant and
`mr-1` geometric `(1,1)` primes over `P=0`.

For `r=1,m>1`, equations (3)--(6) also prove the uniqueness needed by the
weighted--cancellation boundary-intersection theorem: the discriminant divisor is the
unique candidate with a ramified boundary prime of index two, and the second
divisor is the unique candidate with positive unramified boundary
contribution `m-1`.

The arithmetic version groups geometric primes and retains their actual
residue degrees.  Equality (1), not a root-orbit heuristic, remains the
exhaustion certificate.

## 5. Canonical boundary-exhaustion theorem for weighted marked-root theorem and cancellation construction

The preceding calculations can now be stated at the level required by the
canonical boundary invariant, rather than at the level of a chosen
resolvent.

### Theorem 5.1 (complete normalization boundary)

Work geometrically over an algebraic closure of the characteristic-zero
coefficient field.

1. For a generic split weighted marked-root theorem seed of inverse degree `n>=3`, with
   
   \[
    H(W)=hW^2(W-1)\prod_{j=1}^{n-3}(W-\rho_j),
   \]
   
   the irreducible components of
   `bar(X)_F minus A^3` are exactly
   
   \[
    E_\Delta,\quad E_{\rho_1},\ldots,E_{\rho_{n-3}}.
   \]
   
   Their target images are
   
   \[
    \pi_F(E_\Delta)=\Delta_H,\qquad
    \pi_F(E_{\rho_j})=V(C).
   \]
   
   The first component has `(e,f)=(2,1)`; every `E_(rho_j)` has
   `(e,f)=(1,1)`.  When `n=3`, the second list is empty and `V(C)` is not a
   boundary-image divisor.

2. For `cancellation type (m,r)`, with `N=r(m+1)+1`, the irreducible components of
   `bar(X)_F minus A^3` are exactly
   
   \[
    E_\Delta,\quad E_\omega
    \quad
    (\omega\text{ a nonzero root of }K_{m,r}
       \text{ other than the distinguished affine root}).
   \]
   
   Their target images are
   
   \[
    \pi_F(E_\Delta)=\Delta_{m,r},\qquad
    \pi_F(E_\omega)=V(P).
   \]
   
   The first component has `(e,f)=(r+1,1)`, and the `mr-1` components
   `E_omega` have `(e,f)=(1,1)`.  When `(m,r)=(1,1)`, the latter list is
   empty and `V(P)` is not a boundary-image divisor.

Over the original coefficient field, the components in either explicit root
list descend in Galois orbits.  In these charts their residue extensions are
the constant extensions generated by the corresponding `rho_j` or `omega`,
so the arithmetic residue degrees are the degrees of their irreducible
factors.  No additional arithmetic boundary prime is created.

**Proof.**  Proposition 1.2 of the canonical boundary theorem says that the
reduced complement of `A^3` in `bar(X)_F` is pure of codimension one and that
every component maps onto a target divisor.  It is therefore enough to prove
support containment and then exhaust every prime over the candidate target
divisors.

For weighted marked-root theorem, normalize the finite marked-root incidence.  Its function field is
`k(A^3)`, so uniqueness of normalization identifies it with `bar(X)_F`; its
regular-reconstruction open is exactly the distinguished source `A^3`.  On
`D(C Delta_H)`, the incidence is finite etale, every marked root is simple,
and `E_W=-c gamma` together with the reconstruction formulas makes `C` and
every reconstruction denominator a unit.  Hence

\[
 \pi_F^{-1}(D(C\Delta_H))\subset A^3,                       \tag{7}
\]

which confines every boundary image to `Delta_H` or `C=0`.  The normalized
double-root chart over the generic point of `Delta_H` and the normalized
zero/root charts over the generic point of `C=0` construct precisely the
primes listed in part 1.  Their contributions are (3)--(4), each totaling
`n`; Lemma 1.1 rules out every further prime over either target divisor.

For cancellation construction, on `D(P Delta_(m,r))` the degree-`N` resolvent is monic after
inverting `P`, is etale, and satisfies `D(T)!=0` at every root.  Formula (20)
of the master construction reconstructs every root regularly and uniquely,
so

\[
 \pi_F^{-1}(D(P\Delta_{m,r}))\subset A^3.                   \tag{8}
\]

Thus every boundary image is `Delta_(m,r)` or `P=0`.  The normalized critical
chart constructs `E_Delta`; the regular projective chart `T=U/P` constructs
the `U=0` affine prime, the distinguished nonzero affine prime, and the
`mr-1` displayed boundary primes.  The contributions (5)--(6) total `N`, so
Lemma 1.1 rules out every additional prime over both candidate divisors.

Finally, the explicit `C=0` and `P=0` chart equations descend over the
coefficient field.  Factoring their squarefree residual root polynomials
groups the displayed geometric primes into arithmetic primes, with residue
degree equal to the corresponding factor degree.  Applying Lemma 1.1 over
the original field again exhausts the total degree, so descent cannot create
an additional arithmetic prime.  QED

### Corollary 5.2 (marking safety)

For every noncubic case in Theorem 5.1, the two displayed target divisors are
the only target vertices of the canonical boundary diagram.  The discriminant
vertex is intrinsically the unique one receiving a ramified boundary prime;
the other is the unique vertex receiving only unramified boundary primes.
Consequently no unlisted boundary divisor can duplicate a profile, change the
incidence diagram, or enlarge the automorphism group used in the weighted--cancellation
comparison.

This conclusion concerns the complete canonical normalization.  The root
equations are used only to prove (7)--(8) and to construct enough normalized
DVR branches to saturate the local degree formula.

## 6. Independent audit

Theorem 2.1 and Lemma 1.1 are abstract algebraic statements independent of
the two families.  The family-specific inputs to Theorem 5.1 have now been
rederived in the
[independent boundary-exhaustion audit](BOUNDARY_EXHAUSTION_INDEPENDENT_AUDIT.md).
That proof starts from the defining inverse incidences and independently
checks these four points:

1. the resolvent incidence has the same function field as the polynomial map;
2. its normalization is the canonical finite normalization over the target;
3. the displayed reconstruction open is exactly `A^3`, especially in the
   degenerate `C=0` and `P=0` charts; and
4. the local branches used in (3)--(6) are distinct normalization primes with
   the asserted `e` and `f`.

Its standard-library checker verifies the inverse derivative and degree,
the projective chart scaling and `P=0` factorization, squarefreeness of the
boundary polynomial, the parameter-root identity, and both local degree
saturations without importing a project module.  The all-parameter proof is
the valuation and normalization argument in the audit, not extrapolation from
the bounded regression grid.

This closes the previously identified independent boundary-exhaustion risk.
It does not independently reprove every other part of the cancellation construction, such as polynomial
cancellation, monodromy, parameter arithmetic, or unrestricted equivalence.

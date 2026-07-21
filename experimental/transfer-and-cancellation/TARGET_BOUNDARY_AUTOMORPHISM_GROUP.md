# Target automorphisms of the C24 marked boundary pair

The target-fixed theorem reduces equivalence of distinct C24 parameter roots
to automorphisms of the common finite normalization covering nonidentity
target automorphisms.  This note begins the resulting, smaller problem by
computing the restriction of every labelled-boundary automorphism to the
distinguished plane `P=0`.

Fix a noncubic `C24_(m,r)` map over a characteristic-zero field `k`, and put

\[
 Y=\mathbb A^3_{P,Q,R},\qquad D_P=V(P),\qquad
 D_\Delta=V(\Delta_{m,r}).                                 \tag{1}
\]

The labels distinguish these divisors intrinsically: `D_Delta` is the unique
target divisor receiving a ramified boundary prime, while `D_P` receives only
unramified boundary primes.  Let

\[
 \operatorname{Aut}(Y;D_P,D_\Delta)
\]

mean polynomial target automorphisms preserving both labelled divisors.  Such
an automorphism automatically preserves their full scheme-theoretic
intersection.

## 1. Restriction to the boundary plane

The canonical intersection on `D_P` has coordinate ring

\[
 k[Q,R]/\left(
 Q^M\{(r+1)RQ^m-C\}
 \right),
 \qquad M=mr(m+1)>0.                                      \tag{2}
\]

Its reduction is the disjoint union

\[
 V(Q)\simeq\mathbb A^1_R,
 \qquad
 V((r+1)RQ^m-C)\simeq\mathbb G_m.                         \tag{3}
\]

### Theorem 1.1 (boundary-plane restriction)

For every
`beta in Aut(Y;D_P,D_Delta)`, there are unique `a,u in k^*` such that

\[
 \beta^*P=aP,                                             \tag{4}
\]

and modulo `P`,

\[
 \boxed{\quad
 \beta^*Q=uQ,\qquad \beta^*R=u^{-m}R.
 \quad}                                                   \tag{5}
\]

In particular, the restriction of the labelled-pair automorphism group to
`D_P` is contained in a one-dimensional torus; no translation, triangular
shear, or interchange of the two intersection components is possible.

**Proof.**  Preservation of the prime divisor `V(P)` gives (4), because the
only units of `k[P,Q,R]` are constants.  Let `bar(beta)` be the induced
automorphism of `k[Q,R]`.

The two components in (3) cannot be exchanged: their coordinate rings have
different unit groups, so `A^1` and `G_m` are not isomorphic.  Hence
`bar(beta)` preserves the prime `(Q)`, giving

\[
 \bar\beta^*Q=uQ,qquad u\in k^*.
\]

After composing with this scaling in `Q`, an automorphism of the polynomial
ring `k[Q][R]` over `k[Q]` is affine linear in `R`.  Therefore

\[
 \bar\beta^*R=vR+f(Q),qquad v\in k^*,\quad f\in k[Q].     \tag{6}
\]

Apply (5)--(6) provisionally to the second factor of (2).  Preservation of
its prime ideal gives

\[
 (r+1)(vR+f(Q))(uQ)^m-C
 =b\{(r+1)RQ^m-C\}                                       \tag{7}
\]

for some `b in k^*`.  Comparing constant terms gives `b=1`; comparing the
`RQ^m` term gives `vu^m=1`; and the remaining term gives `f=0`.  This proves
(5).  Uniqueness is immediate.  QED

The nonreduced exponent `M` is not needed to derive (5), but its preservation
confirms that the entire scheme intersection, rather than only its reduced
support, is respected.

## 2. The visible scaling subgroup

For every `lambda in k^*`, define

\[
 \sigma_\lambda(P,Q,R)
 =\left(\lambda^{m+1}P,\lambda Q,\lambda^{-m}R\right).      \tag{8}
\]

### Proposition 2.1 (cover-preserving torus)

The maps `sigma_lambda` form a subgroup `G_m` of the labelled-boundary
automorphism group.  They lift to the common C24 finite normalization by

\[
 T\longmapsto\lambda^{-m}T.                               \tag{9}
\]

They preserve every C24 affine reconstruction marking and act trivially on
the parameter roots.

**Proof.**  Under (8)--(9),

\[
 Q-PT\longmapsto\lambda(Q-PT),
 \qquad T(Q-PT)^m\longmapsto T(Q-PT)^m.                   \tag{10}
\]

Changing the integration variable in the inverse resolvent gives

\[
 C\int_0^{\lambda^{-m}T}
 \{1-t(\lambda Q-\lambda^{m+1}Pt)^m\}^r\,dt
 =\lambda^{-m}C\int_0^T\{1-t(Q-Pt)^m\}^r\,dt.             \tag{11}
\]

Thus (8) preserves the cover and its discriminant.  In the projective chart
`U=PT`, one has `U -> lambda U`; hence

\[
 w={U\over Q}\longmapsto w.                              \tag{12}
\]

The filled branch is `w=-q/(1-q)`, so every `q` is fixed.  The source lift is
explicitly

\[
 (x,y,z)\longmapsto
 (\lambda^{-m}x,\lambda y,\lambda^{m+1}z).                 \tag{13}
\]

It fixes `A=1+xy^m` and sends `B` to `lambda^(m+1)B` for the same cancellation
polynomial `h`.  It therefore preserves the distinguished affine open and
does not change the parameter branch.  QED

## 3. Exact reduction of the remaining group problem

Theorem 1.1 defines a restriction character

\[
 \rho:\operatorname{Aut}(Y;D_P,D_\Delta)\longrightarrow\mathbb G_m,
 \qquad \rho(\beta)=u.                                    \tag{14}
\]

Proposition 2.1 splits it.  After composing any `beta` with
`sigma_(rho(beta))^(-1)`, the resulting automorphism lies in the congruence
kernel

\[
 \mathcal K=left\{\beta:
 \beta^*Q\equiv Q,
 \ \beta^*R\equiv R\pmod P,
 \ \beta^*P=aP
 \right\}.                                                \tag{15}
\]

Consequently

\[
 \operatorname{Aut}(Y;D_P,D_\Delta)
 =\mathcal K\rtimes\mathbb G_m                            \tag{16}
\]

whenever `rho` is read as the split restriction map.  The scaling factor acts
trivially on the parameter roots.  Therefore **all possible motion of the
parameter roots is concentrated in `mathcal K`**.

For the equivalence problem one needs the still smaller subgroup
`K_cov subset K` whose elements lift to the common finite normalization and
preserve its labelled upstairs--downstairs boundary object.  Such a lift acts
on the finite set of nonzero roots of `K_(m,r)(w)`, and hence gives a
homomorphism

\[
 \mathcal K_{\rm cov}longrightarrow
 \operatorname{Perm}\{w:K_{m,r}(w)=0\}.                   \tag{17}
\]

The parameter root `q` is recovered injectively from its filled root by

\[
 q=-{w\over1-w}.                                          \tag{18}
\]

Thus any equivalence between distinct normalized parameter roots must come
from a nontrivial element of `K_cov` acting nontrivially in (16).

This replaces unrestricted left--right equivalence by two concrete tasks:

1. prove `K_cov=1`, which would finish parameter-root inequivalence; or
2. compute its finite permutation image on the roots of `K_(m,r)`.

The full target-only kernel `K` could conceivably contain automorphisms that
do not lift to the inverse cover.  Such automorphisms are irrelevant to C24
equivalence and should not be included merely because they preserve the two
target equations.

## 4. What remains

Theorem 1.1 rules out all boundary-plane freedom except weighted scaling, and
that scaling fixes `q`.  The unresolved issue is now an infinitesimal one:
can a polynomial automorphism congruent to the identity on `P=0`, possibly
with `P` rescaled, preserve the full discriminant and lift so as to permute
the projective normalization branches?

The next useful calculation is the first formal neighborhood of `P=0` in the
normalized resolvent chart.  It should determine the action of `K_cov` on
`w` and is strictly smaller than another cancellation classification.

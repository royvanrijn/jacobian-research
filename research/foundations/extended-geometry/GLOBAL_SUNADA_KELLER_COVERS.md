# Global Sunada pairs of relative Keller maps

This note promotes the degree-seven Davenport calculation from isolated
vertical specializations to a single pair of nonisomorphic covers over a
common two-dimensional target.  There are two complementary Keller
realizations.  One primitive Cox-ledger coordinate makes the two covers
determinant-one morphisms of smooth affine threefolds to the exact same
target.  After a finite tangent-mark base change, the covers are also literal
inverse-cover pullbacks of two stably inequivalent relative weighted Keller
maps.

Neither realization is an absolute polynomial self-map of affine
three-space: the direct Cox version lives on boundary complements, and the
weighted version is relative over a tangent-mark curve.  Constructing two
new polynomial self-maps of one affine space whose full inverse covers form
a Gassmann pair remains a controlled-boundary suspension problem.

## 1. The common global Sunada base

Let

\[
 K=\mathbb Q(a),\qquad a^2+a+2=0,
\]

and let \(g_T(Y)\) and \(h_T(Z)=\overline{g_T}(Z)\) be the degree-seven
Davenport polynomials from
[the arithmetic-indistinguishability note](ARITHMETIC_INDISTINGUISHABILITY.md).
Use \(u\) for the polynomial value.  Exact elimination gives

\[
 \operatorname {Res}_Y(g_T'(Y),g_T(Y)-u)
 =\operatorname {Res}_Z(h_T'(Z),h_T(Z)-u)
 =7^{-6}\Delta(T,u)^2,                                  \tag{1.1}
\]

where

\[
\begin{aligned}
\Delta(T,u)={}&17920T^{10}+10472T^9-2464T^8
 -1792T^7u-1728T^7\\
&+11956T^6u+7056T^5u-4802T^3u^2+343u^3.
\end{aligned}                                           \tag{1.2}
\]

The cubic \(\Delta\) is irreducible in \(K(T)[u]\).  Define

\[
 B=\operatorname {Spec}K[T,u,\Delta^{-1}]
\]

and the two finite etale degree-seven covers

\[
 X_g=V(g_T(Y)-u)\longrightarrow B,\qquad
 X_h=V(h_T(Z)-u)\longrightarrow B.                     \tag{1.3}
\]

These are global over the surface \(B\): \(T\) is not fixed.

The cubic

\[
\begin{aligned}
Q_T(Y,Z)={}&Y^3+aY^2Z+Y^2Z+aYZ^2-Z^3\\
 &+(5+3a)TY+(-2+3a)TZ+(2a+1)T
\end{aligned}                                           \tag{1.4}
\]

divides \(g_T(Y)-h_T(Z)\).  The common Galois closure has group

\[
 \Gamma=\operatorname {GL}_3(\mathbb F_2)
\]

and the two degree-seven actions are its actions on Fano points and Fano
lines.  If \(P\) and \(L\) are a point and a line, respectively, then

\[
 H_P=\operatorname {Stab}_\Gamma(P),\qquad
 H_L=\operatorname {Stab}_\Gamma(L)                    \tag{1.5}
\]

have order \(24\), are not conjugate, and have equal induced permutation
characters.  Therefore \(X_g\) and \(X_h\) are nonisomorphic \(B\)-covers
with one common Galois closure.

This is the global Sunada statement missing from the earlier
specialization-only formulation.

## 2. Equality of every good fiber zeta function

For a finite etale zero-dimensional scheme \(X/\mathbb F_q\), if geometric
Frobenius has cycle lengths \(d_1,\ldots,d_j\), then

\[
 Z(X,v)=\prod_{i=1}^j(1-v^{d_i})^{-1}.                 \tag{2.1}
\]

For every \(\gamma\in\Gamma\), the point and line permutations have the same
cycle partition.  One way to see the strengthening from character equality
to cycle equality is to apply the common character to every power
\(\gamma^e\); the fixed-point counts of all powers recover the number of
cycles of each length.  The exact checker also compares the cycle partitions
directly for all \(168\) elements.

Spread the common Galois closure out after inverting a finite set of primes
containing \(2\) and \(7\).  A **good fiber** means a closed point of this
model lying off \(\Delta=0\), at a prime where that common
\(\Gamma\)-torsor model is finite etale.  Its Frobenius is a conjugacy class
in \(\Gamma\), so

\[
 \boxed{Z((X_g)_b,v)=Z((X_h)_b,v)}                     \tag{2.2}
\]

for every good closed point \(b\).  Equivalently, the two fibers have the
same residue-degree multiset and the same number of points over every finite
extension of \(\kappa(b)\).  This is pointwise equality, not merely equality
of Chebotarev densities.

As a finite regression, the checker factors both polynomials on every
unramified rational \((T,u)\)-fiber over \(\mathbb F_{11}\) and
\(\mathbb F_{23}\), the first two split good residue fields.  All \(620\)
tested fibers have identical factor-degree multisets.

## 3. Direct Cox-ledger Kellerization over one target

There is a general one-unit observation behind the first
[three-factor Cox ledger](COX_LEDGER_THREE_FACTOR.md).  If

\[
 f:X\longrightarrow Y
\]

is a finite etale morphism of smooth affine \(d\)-folds and chosen volume
forms give

\[
 f^*\omega_Y=j\omega_X,\qquad j\in\mathcal O(X)^*,
\]

then

\[
 \widehat f:X\times\mathbb A^1_z\longrightarrow
 Y\times\mathbb A^1_Z,\qquad
 (x,z)\longmapsto\left(f(x),\frac z{j(x)}\right)       \tag{3.1}
\]

has determinant one.  Indeed,

\[
 f^*\omega_Y\wedge d(z/j)
 =j\omega_X\wedge\left(\frac{dz}{j}
                 -\frac{z\,dj}{j^2}\right)
 =\omega_X\wedge dz,                                  \tag{3.2}
\]

because \(\omega_X\wedge dj=0\).  The extra target coordinate determines
\(z=jZ\), so the suspension does not alter any finite fiber.

For (1.3), use source coordinates \((T,Y)\) and \((T,Z)\).  Their Jacobian
units are

\[
 j_g=g_T'(Y),\qquad j_h=h_T'(Z).                       \tag{3.3}
\]

These derivatives are units after restricting to \(B\).  This follows
scheme-theoretically from the Bezout identity for the resultants (1.1), not
only pointwise from the absence of repeated roots.  Hence

\[
\begin{aligned}
\widehat\pi_g:X_g\times\mathbb A^1_z&\longrightarrow
 B\times\mathbb A^1_Z,
 &(T,Y,z)&\longmapsto(T,g_T(Y),z/g_T'(Y)),\\
\widehat\pi_h:X_h\times\mathbb A^1_z&\longrightarrow
 B\times\mathbb A^1_Z,
 &(T,Z,z)&\longmapsto(T,h_T(Z),z/h_T'(Z))
\end{aligned}                                         \tag{3.4}
\]

are finite etale degree-seven morphisms with determinant one and the exact
same smooth affine threefold target.  Their common Galois closure is the
pullback of the common \(\Gamma\)-closure in Section 1, their fibers retain
the zeta equality (2.2), and nonconjugacy in (1.5) makes them nonisomorphic
over the target.

These are **Cox-Keller morphisms**, in the terminology of the three-factor
ledger: constant-Jacobian maps between smooth affine boundary complements.
They solve the common-target determinant problem without pretending that
the source and target complements are affine spaces.

## 4. Simultaneous relative weighted-Keller realization

The remaining issue is that a weighted seed needs a distinguished affine
tangent mark.  Introduce nonzero \(r,\rho\) satisfying

\[
\begin{aligned}
J_g(T,r)&=
\frac{g_T(r)-g_T(0)-r g_T'(0)}{r^2}=0,\\
J_h(T,\rho)&=
\frac{h_T(\rho)-h_T(0)-\rho h_T'(0)}{\rho^2}=0.
\end{aligned}                                          \tag{4.1}
\]

Both equations have degree five in their marking variable.  Let
\(\mathcal U\) be the open part of their fiber product over the \(T\)-line
on which

\[
r\rho c_gc_h(\kappa_g+2)(\kappa_h+2)\ne0,              \tag{4.2}
\]

where the quantities below define \(c_\bullet,\kappa_\bullet\).  Put

\[
\begin{aligned}
H_g(W)&=g_T(rW)-g_T(0)-r g_T'(0)W,\\
H_h(W)&=h_T(\rho W)-h_T(0)-\rho h_T'(0)W.
\end{aligned}                                          \tag{4.3}
\]

Equation (4.1) gives

\[
H_\bullet(0)=H_\bullet'(0)=H_\bullet(1)=0.
\]

Set

\[
c_\bullet=-H_\bullet'(1),\qquad
\kappa_\bullet=\frac{H_\bullet''(1)}{c_\bullet},\qquad
b_\bullet=c_\bullet^{-1}.                              \tag{4.4}
\]

The universal weighted-seed theorem now constructs

\[
F_g,F_h:\mathbb A^3_{\mathcal U}\longrightarrow
              \mathbb A^3_{\mathcal U}
\]

with relative Jacobians

\[
\det D_{\mathcal U}F_g=b_gc_g=1,\qquad
\det D_{\mathcal U}F_h=b_hc_h=1.                       \tag{4.5}
\]

Thus these are relative Keller maps, and the maps obtained by adjoining the
identity on \(\mathcal U\) are etale wherever \(\mathcal U\) is smooth.

Let

\[
 \mathcal L=(\mathcal U\times\mathbb A^1_u)
             \cap(\Delta\ne0).
\]

For weighted target coordinates \((A,B,C)\), embed this one common abstract
target into the two Keller targets by

\[
\begin{array}{lll}
\iota_g:& C=1,&
B=-r g_T'(0),\quad c_gA=g_T(0)-u,\\[2mm]
\iota_h:& C=1,&
B=-\rho h_T'(0),\quad c_hA=h_T(0)-u.
\end{array}                                             \tag{4.6}
\]

The two inverse pencils restrict exactly to

\[
\begin{aligned}
H_g(W)-BCW+c_gAC^2&=g_T(rW)-u,\\
H_h(W)-BCW+c_hAC^2&=h_T(\rho W)-u.                     \tag{4.7}
\end{aligned}
\]

Since \(r,\rho\) are units and \(\Delta\ne0\), every root is simple and the
weighted reconstruction formula is regular.  Consequently

\[
F_g^{-1}(\iota_g(\mathcal L))\simeq X_g\times_B\mathcal L,
\qquad
F_h^{-1}(\iota_h(\mathcal L))\simeq X_h\times_B\mathcal L. \tag{4.8}
\]

Equations (1.5), (2.2), and (4.8) are the requested global Sunada pair of
Keller inverse covers over one positive-dimensional target.

## 5. Stable inequivalence

For \(T\ne0\), both \(g_T''\) and \(h_T''\) have leading coefficient \(6\)
and zero degree-four coefficient.  If their effective divisors were related
by \(Y=\alpha Z+\beta\), coefficient comparison would first give
\(\beta=0\).  The degree-three and degree-two coefficients would then give

\[
\alpha^2=-\frac{1+a}{a},\qquad
\alpha^3=-\frac{1+a}{a}.                               \tag{5.1}
\]

Thus \(\alpha=1\) and \(2a+1=0\), contradicting \(a^2+a+2=0\).  Affine
normalization by \(r\) and \(\rho\) cannot change this conclusion.

At the generic point of \(\mathcal U\), the two normalized Hessian divisors
are therefore inequivalent.  Stable normalization functoriality and the
weighted Hessian obstruction imply that \(F_g\) and \(F_h\) are not stably
polynomially left-right equivalent over the generic base.  Any relative
stable equivalence would restrict to such an equivalence at the generic
point, so the relative families themselves are stably inequivalent.

## 6. Why this does not solve the absolute affine-space problem

The weighted construction uses \(T\) as a base coordinate and adjoins the
finite tangent marks \(r,\rho\).  It gives an etale morphism of affine
three-space bundles over \(\mathcal U\).  The direct construction (3.4)
has one common target and no tangent extension, but its source and target
are affine boundary complements.  Neither is a new polynomial self-map of
\(\mathbb A^3_K\) with a full-target \(\Gamma\)-inverse cover.

The marking extension itself is now explicit: its normalization is a
rational conic, but the nonzero-mark locus removes two rational points.
Therefore it has no nonconstant affine-line parametrization.  See the
[tangent-mark audit](DAVENPORT_TANGENT_MARK_CURVE.md).

That obstruction applies to the fixed mark and to marks affine-linear in
`T`, but not to marks linear in a square-root parameter.  The
[proportional-section audit](DAVENPORT_PROPORTIONAL_TANGENT_SECTIONS.md)
finds four exact sections `T=tau*s^2`; the conjugation-invariant choice
`T=-s^2` marks both Davenport covers simultaneously.  After removing a
common `s^6` factor, the weighted seeds extend across `s=0`, leaving the
three explicit determinant/endpoint divisors recorded there.

There are two independent obstructions to claiming more:

1. allowing the second tangent-pencil parameter restores universal
   \(S_7\)-monodromy;
2. turning a prescribed Davenport core into an absolute Keller map requires
   a polynomial suspension chart whose divisor ledger cancels all six
   critical branches.

The direct two-boundary multiplicative ansatz acquires the extra divisor
\(N_{d_1,d_2}=0\).  The new three-factor Cox ledger and the suspension
(3.1) show how a surviving boundary unit can instead be absorbed by a
primitive coordinate.  What they do not provide is an affine-space chart:
the source and target remain boundary complements.  For the Davenport
derivative, which has six generic simple critical points grouped over three
values, an absolute construction must simultaneously straighten those
complements and retain the unit cancellation.  This is now a sharply stated
factoriality/coordinate problem rather than a missing group-theory problem.
The complete height-one calculation and the stable unit-rank obstruction are
proved in the
[Davenport Cox-boundary audit](DAVENPORT_COX_BOUNDARY_OBSTRUCTION.md).

The cancellation parameter polynomial does not by itself solve this
remaining problem.  Its truncated-binomial, geometric-derivative, Shabat,
and Schubert incarnations govern the field of the cancellation coefficient
branch.  The generic inverse covers of the cancellation maps still have
\(A_N\) or \(S_N\) monodromy.  A useful dessin input would have to prescribe
the **inverse-core** passport and arrive with a compatible suspension ledger,
not only prove irreducibility of the coefficient parameter.

## 7. Reproduction

Run

```bash
.venv/bin/python scripts/verify_global_sunada_keller.py
```

The checker verifies:

1. the global cubic correspondence;
2. the common squared branch resultant (1.1);
3. the determinant-one direct Cox suspensions (3.4);
4. the finite tangent-mark equations and exact pullback identities (4.7);
5. the Hessian-divisor obstruction;
6. equal point/line cycle partitions for all \(168\) elements of
   \(\operatorname {GL}_3(\mathbb F_2)\);
7. nonconjugacy of the two order-\(24\) stabilizers; and
8. all unramified rational fibers over \(\mathbb F_{11}\) and
   \(\mathbb F_{23}\).

The symbolic and permutation data are reusable from
[`jcsearch/sunada.py`](../jcsearch/sunada.py).

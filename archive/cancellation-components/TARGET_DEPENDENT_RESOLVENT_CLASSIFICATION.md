# Classification of target-dependent polynomial resolvents

This note gives the maximal derivative classification inside the generalized
cancellation coordinate skeleton.  Instead of prescribing a power or a product of
normalized factors, allow the inverse-resolvent derivative to be an arbitrary
polynomial in the inverse variable and the target coordinates.  The Keller
condition forces the original cancellation power identically.

Let `k` have characteristic zero, let `f in k[y]` be nonconstant, let
`a,b>=1`, and put

\[
 A=1+xf(y),\qquad B=A^bz+g(y,A),                              \tag{1}
\]

\[
 P=A^aB,\qquad Q=y+xA^{a-1}B,
 \qquad s={x\over A},\qquad n=a+b-2.                          \tag{2}
\]

For an arbitrary nonzero polynomial `H in k[T,P,Q]`, define in the localized
ring

\[
 R=C\int_0^s H(t,P,Q)dt,
 \qquad C\in k^*.                                             \tag{3}
\]

The lower endpoint only fixes the additive target normalization of `R`.

## 1. Algebraic independence

The standard localized calculation gives

\[
 \det {\partial(s,P,Q)\over\partial(x,y,z)}=-A^n.             \tag{4}
\]

In particular, `s,P,Q` are algebraically independent.  Equivalently, the
localized reconstruction is

\[
 y=Q-sP,qquad
 D=1-sf(Q-sP),qquad
 A=D^{-1},qquad x=sD^{-1},                                  \tag{5}
\]

followed by `B=PD^a` and reconstruction of `z` from (1).  Thus substitution

\[
 k[T,P,Q]\longrightarrow k(s,P,Q),\qquad T\longmapsto s       \tag{6}
\]

is injective.

## 2. Universal derivative rigidity

### Theorem 2.1

The localized map `(P,Q,R)` has nonzero constant Jacobian if and only if

\[
 \boxed{\quad
 H(T,P,Q)=\lambda\{1-Tf(Q-PT)\}^n,
 \qquad \lambda\in k^*.
 \quad}                                                       \tag{7}
\]

Hence target-dependent polynomial coefficients, sums, and arbitrary finite
factorizations produce no derivative beyond the original cancellation power.

**Proof.**  Holding `(P,Q)` fixed in (3) gives

\[
 R_s=C H(s,P,Q).
\]

Combining this with (4), the Jacobian is

\[
 -C A^nH(s,P,Q).                                             \tag{8}
\]

It is a nonzero constant exactly when

\[
 H(s,P,Q)=\lambda A^{-n}.
\]

By (5), `A^(-1)=1-sf(Q-Ps)`.  Therefore

\[
 H(s,P,Q)=\lambda\{1-sf(Q-Ps)\}^n.                           \tag{9}
\]

The injectivity in (6) upgrades (9) to the polynomial identity (7).
Conversely, substituting (7) in (8) gives the constant `-C lambda`.  QED

### Corollary 2.2 (complete polynomial classification in the skeleton)

If the three coordinates `(P,Q,R)` are polynomial and have nonzero constant
Jacobian, then the map is, up to translations, nonzero scalings, and
polynomial left--right equivalence, a cancellation map.

**Proof.**  Theorem 2.1 forces the derivative used in the
[generalized cancellation theorem](GENERALIZED_CANCELLATION_MECHANISM.md).
If `n=0`, then `H=lambda` and `R=C lambda x/A`, which is not polynomial.
Thus `n>=1`.  The generalized theorem forces `a=1`, `b=n+1`, classifies the
polynomial cancellation jet spectrally, and removes its tail by a polynomial
source automorphism.  Its resulting map is cancellation-equivalent.  QED

The
[finite-factor theorem](TWO_FACTOR_RESOLVENT_CLASSIFICATION.md) is now also a
direct corollary: any factorized polynomial derivative must multiply out to
(7).  Its Laurent-unit proof remains useful because it classifies the
individual factors without first expanding them.

This theorem completes Option B within the coordinate skeleton (1)--(2).
A genuinely new branch must change the reconstruction skeleton itself, for
example by introducing another source function or another inverse variable.

## 3. Exact regression

Run

```bash
.venv/bin/python scripts/verify_target_dependent_resolvent.py
```

The script differentiates nonlinear representatives and independently solves
bounded generic coefficient systems for `H`.  The all-degree result is the
algebraic-independence proof above.

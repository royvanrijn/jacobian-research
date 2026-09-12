# F2 `k=1` affine target collision theorem

> **Status.**  On the first F2 nonproperness normalization chart
> `(deg p,deg q)=(3,5)`, affine reparametrization and affine target changes
> give
> `p=t^3+a*t` and `q=t^5+b*t^4+c*t^2+d*t`.  All normalization collisions
> are indexed by the roots of one explicit quartic
> `R(u)=u^4+b*u^3+a*u^2+(2*a*b-c)*u-(a^2+d)`.  A root is diagonal exactly
> when `3*u^2+4*a=0`; otherwise it gives the unordered pair of roots of
> `z^2-u*z+(u^2+a)`.  There is a nonempty open subchart with four distinct
> ordinary affine nodes.  Its unique infinity branch has type `(2,5)` and
> delta invariant `2`, while the four nodes contribute `4`; together they
> exhaust the arithmetic genus `6` of a rational plane quintic.  This is an
> exact conductor/collision atlas for `k=1`, not an F2 exclusion or a
> construction of the source pullback.

The formulas, generic witness, tangent test, and genus ledger are replayed by
[`verify_f2_affine_target_k1_collision.py`](../scripts/verify_f2_affine_target_k1_collision.py).

## 1. Four-parameter normal form

Let `nu(t)=(p(t),q(t))` be the normalization parametrization on the `k=1`
chart of the
[`target-curve atlas`](F2_AFFINE_TARGET_CURVE_ATLAS.md).  Thus

\[
 \deg p=3,\qquad \deg q=5.                       \tag{1.1}
\]

An affine change of parameter makes `p` monic and removes its quadratic
term.  Target translations remove the two constants, independent target
scalings normalize both leading coefficients, and the target shear
`q -> q-lambda*p` removes the cubic term of `q`.  Hence

\[
 \boxed{
 p(t)=t^3+at,\qquad
 q(t)=t^5+bt^4+ct^2+dt.}                         \tag{1.2}
\]

All these operations preserve the F2 coordinate-degree ratio and the Keller
condition up to a nonzero constant rescaling.

## 2. Exact collision quartic

For two normalization parameters `s,t`, put

\[
 u=s+t,\qquad v=st.                              \tag{2.1}
\]

The divided difference of `p` is

\[
 \Delta_p=s^2+st+t^2+a=u^2-v+a.                 \tag{2.2}
\]

Thus a collision must satisfy

\[
 v=u^2+a.                                        \tag{2.3}
\]

Substitution in the divided difference of `q` gives

\[
 \Delta_q=-R(u),                                 \tag{2.4}
\]

where

\[
 \boxed{
 R(u)=u^4+bu^3+au^2+(2ab-c)u-(a^2+d).}           \tag{2.5}
\]

Conversely, for every root of `R`, the roots `s,t` of

\[
 z^2-uz+(u^2+a)=0                                \tag{2.6}
\]

have the same `p` and `q` values.  The discriminant of (2.6) is

\[
 (s-t)^2=-3u^2-4a.                               \tag{2.7}
\]

Therefore the root is a diagonal critical point exactly when
`3u^2+4a=0`; every other root is an off-diagonal normalization collision.

Reducing (1.2) modulo (2.6) gives the common target point explicitly:

\[
 \boxed{
 x(u)=-u(u^2+a),\qquad
 y(u)=(u^2+a)(u^3+2au+ab-c).}                    \tag{2.8}
\]

Equations (2.5)--(2.8) replace a two-variable collision elimination by one
quartic and two evaluation formulas.

## 3. The generic four-node packet

The following conditions are Zariski open in `(a,b,c,d)`:

1. `R` is squarefree;
2. `Res(R,3u^2+4a)` is nonzero, so every pair is off diagonal;
3. the two branch tangent vectors at every pair are distinct; and
4. the four target points (2.8) are distinct.

The open set is nonempty.  At

\[
 (a,b,c,d)=(1,0,0,0),                            \tag{3.1}
\]

one has

\[
 R=u^4+u^2-1,qquad
 \operatorname{Disc}(R)=-400,qquad
 \operatorname{Res}(R,3u^2+4)=25.               \tag{3.2}
\]

The branch-tangent resultant is `-10000`, and saturation by `u-v` proves
that no two distinct roots of `R` have the same pair `(x(u),y(u))`.
Consequently this member, and hence a dense open subchart, has four distinct
ordinary affine nodes.  It is generically injective and therefore is indeed
a normalization parametrization of a rational quintic.

## 4. Complete generic conductor budget

At infinity, in the target chart `q=1`, the normalization has local orders

\[
 \operatorname{ord}_\infty(p/q)=2,qquad
 \operatorname{ord}_\infty(1/q)=5.              \tag{4.1}
\]

The unique infinity branch is the `(2,5)` plane cusp, with

\[
 \delta_\infty=\frac{(2-1)(5-1)}2=2.            \tag{4.2}
\]

A plane quintic has arithmetic genus

\[
 p_a=\frac{(5-1)(5-2)}2=6.                      \tag{4.3}
\]

The normalization has genus zero.  Four ordinary nodes contribute four more
delta units, so

\[
 \boxed{p_a=\delta_\infty+4\delta_{\rm node}=2+4=6.} \tag{4.4}
\]

There are no further singularities on the generic `k=1` target.  The
degenerate complement is explicitly reached when the collision quartic has
multiple roots, a root becomes diagonal, branch tangents coincide, or target
collision values merge.  Those strata contain the cusp, tacnode, and
higher-multiple-fiber packets that a complete pullback compiler must retain.

## 5. Consequence and next target

For `k=1`, the target conductor is no longer unspecified: generically it is
the direct sum of four affine nodal conductor quotients, with the separate
fixed cusp at infinity.  The remaining source calculation is to substitute
the exposed F2 Laurent data in the implicit equation and factor the pullback
while retaining the quartic roots as marked target points.  The subsequent
[`implicit-conductor theorem`](F2_AFFINE_TARGET_K1_IMPLICIT_CONDUCTOR.md)
completes the first of those operations: it gives the exact twelve-support
quintic and packages all eight nodal preimages in one conductor polynomial.

<!-- status-consumer: PF2K1I1 a7582c1e36140840 -->

The later
[`fixed-coordinate Keller-pullback theorem`](F2_AFFINE_K1_KELLER_PULLBACK.md)
transports these four node values back to the fixed F2 target chart and
proves that their affine source pullbacks are reduced ordinary nodes with
normalization-defect length one and conductor-divisor degree two.  Their four
fiber counts, rather than their local
types, are the remaining affine calculation.

<!-- status-consumer: PF2K1PB1 6f837229017243c4 -->

The
[`all-stratum conductor theorem`](F2_AFFINE_TARGET_K1_CONDUCTOR_CONSERVATION.md)
shows that the discriminant complement does not create a new total-conductor
problem: every degeneration has affine delta four and the same exact
degree-eight conductor divisor, with multiplicities recording the collided
branches.

<!-- status-consumer: PF2K1CC1 f152c82ef2d54c32 -->

The subsequent
[`generic affine-row Chern theorem`](F2_AFFINE_K1_LOG_CH2.md) already
computes the divisorial conormal contribution from `(e,f,E^2)` and the
target carrier contact.  In the special-residue chart that contact is the
truncated jet order `b=min(ord_u(w|_C),8)`; the remaining boundary
factorization must now supply the source integers and boundary-supported
`Fitt_1` corrections.

This theorem does not prove that the purity-forced component uses `k=1`, and
it does not determine the boundary indices above the four nodes.  The other
23 degree charts and all degenerate strata remain possible.  No exclusion of
`(75,125)` or proof of `JC(2)` is claimed.

The subsequent puncture theorem locates this quintic transversely on the
target `(5,2)` divisor.  The value `125/729` selects the special carrier
point, but the already resolved terminal neighborhood cannot extract a
divisor dominating the quintic and therefore does not determine its purity
index.

<!-- status-consumer: PF2PPA1 b24c4d80c2f8230e -->

On the generic cusp degeneration, conditional minimal transverse SNC
boundary attachments contribute total point length `2f`; this is a boundary
correction, not the affine cusp conductor length or raw jet corank.

<!-- status-consumer: PF2K1L1 5221f5659fc19729 -->

The complement-monodromy theorem now covers every immersed, distinct-image
root partition of this quartic, not only the four-node open set.  It also
covers the generic `A_2+3A_1` and `2A_2+2A_1` cusp strata.  All seven
affine complements are `Z`; their fixed affine sheets force a second
ramified target component and a second new source-boundary divisor.  The
first noncyclic escape is the `E_6+A_1` stratum, whose exact degree-six
fixed-sheet action shows that topology alone cannot finish the cusp attack.

<!-- status-consumer: PF2K1M1 fafcbb3c2e6ceb2b -->

## Reproduction

```bash
.venv/bin/python scripts/verify_f2_affine_target_k1_collision.py
```

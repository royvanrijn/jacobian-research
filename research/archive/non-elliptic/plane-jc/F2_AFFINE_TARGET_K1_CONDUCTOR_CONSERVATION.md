# F2 `k=1` all-stratum conductor-conservation theorem

> **Status.**  Every degree-`(3,5)` target normal form is a birational
> parametrization of a rational plane quintic.  Its infinity branch always
> has type `(2,5)` and delta invariant two, so its complete affine
> normalization defect has length four on every collision degeneration.
> The degree-eight resultant `C(t)` is the conductor divisor on the affine
> normalization on every stratum, not only on the generic four-node locus.
> Under a Keller pullback, affine normalization-defect length is at most
> `4(d-1)` and affine conductor-divisor degree is at most `8(d-1)`.

The universal resultant, two degeneration witnesses, and degree bounds are
replayed by
[`verify_f2_affine_target_k1_conductor_conservation.py`](../scripts/verify_f2_affine_target_k1_conductor_conservation.py).

## 1. Birationality and the fixed genus budget

Consider the complete `k=1` normal form

\[
 p=t^3+at,\qquad q=t^5+bt^4+ct^2+dt.             \tag{1.1}
\]

Every dominant polynomial parametrization factors through the normalization
by a polynomial self-map of `A1`.  Its degree must divide both coordinate
degrees.  Since `gcd(3,5)=1`, (1.1) is automatically birational onto its
image for every choice of `a,b,c,d`.

The projective image is consequently a rational plane quintic.  At its
unique infinity point, the two local coordinate orders are `(2,5)`.  After
taking a square-root uniformizer for the order-two coordinate, the
order-five term remains nonzero, so this branch has characteristic pair
`(2,5)` on every parameter stratum.  Therefore

\[
 p_a=\frac{(5-1)(5-2)}2=6,
 \qquad \delta_\infty=\frac{(2-1)(5-1)}2=2.      \tag{1.2}
\]

The normalization has genus zero.  Hence the complete affine defect is

\[
 \boxed{
 \sum_{z\in\operatorname{Sing}(C)\cap\mathbb A^2}\delta_z=4.} \tag{1.3}
\]

Nodes may collide into cusps, tacnodes, or higher multiple fibers, but they
cannot change (1.3).

## 2. The resultant is the conductor on every stratum

Let

\[
 R(u)=u^4+bu^3+au^2+(2ab-c)u-(a^2+d)             \tag{2.1}
\]

and

\[
 C(t)=\operatorname{Res}_u
 \left(R(u),t^2-ut+(u^2+a)\right).               \tag{2.2}
\]

Direct elimination gives a monic polynomial of degree eight for all values
of `a,b,c,d`; its degree does not drop on the discriminant.  The exact
gradient identities remain

\[
 F_P(p,q)=q'C,
 \qquad
 F_Q(p,q)=-p'C.                                  \tag{2.3}
\]

For a plane curve, the dualizing differential is represented rationally by

\[
 \frac{dP}{F_Q}=-\frac{dQ}{F_P}.                 \tag{2.4}
\]

Pulling (2.4) to the normalization and using (2.3) gives

\[
 \nu^*\!\left(\frac{dP}{F_Q}\right)
 =-\frac{dt}{C(t)}.                              \tag{2.5}
\]

The adjoint/conductor description of the dualizing module for a reduced
Gorenstein curve therefore identifies `div(C)` with the conductor divisor
on the affine normalization.  Equivalently, at each affine singular point
`z`,

\[
 \deg C|_{\nu^{-1}(z)}=2\delta_z.                \tag{2.6}
\]

Summing (2.6) and using (1.3) gives

\[
 \boxed{\deg C=8=2\delta_{\rm aff}}              \tag{2.7}
\]

on every `k=1` stratum.  Thus the collision resultant is not merely a
set-theoretic marker: its root multiplicities are the exact conductor
exponents.

## 3. Exact degeneration witnesses

At `(a,b,c,d)=(0,0,0,0)`,

\[
 p=t^3,\qquad q=t^5,\qquad F=P^5-Q^3,qquad C=t^8. \tag{3.1}
\]

All four affine delta units concentrate in the `(3,5)` cusp, whose conductor
exponent is eight.

At `(a,b,c,d)=(0,0,1,0)`,

\[
 p=t^3,\qquad q=t^5+t^2,qquad
 C=t^2(t^3+1)^2.                                 \tag{3.2}
\]

The point `t=0` is a `(2,3)` cusp with delta one and conductor exponent two.
The three roots of `t^3=-1` map to one ordinary triple point.  That point has
delta three and conductor exponent two on each branch.  Thus both ledgers
remain

\[
 1+3=4,qquad 2+3\cdot2=8.                       \tag{3.3}
\]

This explicitly demonstrates how the generic eight simple conductor points
specialize without losing length.

## 4. Keller pullback on every stratum

Let the target affine singular points be `z_i`, with delta invariants
`delta_i`, and let `N_i` be the number of affine source points over `z_i`.
The fixed-coordinate Keller-pullback theorem proves that every source germ is
an étale copy of the target germ.  Hence

\[
 \ell(\widetilde{D}/D)_{\rm aff}
 =\sum_iN_i\delta_i,                              \tag{4.1}
\]

while the conductor divisor on the source normalization has affine degree

\[
 \deg\operatorname{Cond}_{\rm aff}
 =2\sum_iN_i\delta_i.                            \tag{4.2}
\]

Every `z_i` lies on the nonproperness set, so its finite normalization-cover
fiber contains a nonempty boundary part.  At geometric degree `d`, therefore,
`0<=N_i<=d-1`.  Combining this with `sum delta_i=4` gives the
all-stratum bounds

\[
 \boxed{
 \ell(\widetilde D/D)_{\rm aff}\le4(d-1),
 \qquad
 \deg\operatorname{Cond}_{\rm aff}\le8(d-1).}   \tag{4.3}
\]

At the squarefree/double degree floors these are respectively `20/44` for
the normalization quotient and `40/88` for the conductor divisor.

## 5. What this closes

The degenerate `k=1` target no longer requires a separate total-conductor
calculation.  The single polynomial `C(t)` carries its exact conductor
multiplicities through every node collision, diagonal critical point, cusp,
tacnode, and higher multiple fiber.  Étale pullback also fixes every affine
source local type once its target value is known.

What can still depend on the degeneration is the distribution of the four
delta units among target values and branches, the corresponding finite fiber
counts, and their attachments to the unresolved source boundary.  Those
data, together with `(e,f,E^2)` and the boundary logarithmic matrices, still
require the complete F2 source pair.  This theorem does not exclude
`(75,125)` or prove `JC(2)`.

The logarithmic interpretation is now separated from this ordinary
conductor ledger.  The
[`tame-node packet theorem`](F2_AFFINE_K1_TAME_NODE_PACKET.md) proves that an
fs tame Kummer toroidal source packet over any resolved target node has zero
logarithmic cotangent cokernel and zero localized `ch_2`, even in the
collided cyclic model with its full exceptional chain.  Consequently the
degree-eight divisor computed here cannot be charged as a logarithmic point
correction without additional non-toroidal source data.  Even there,
full-rank SNC exponent data kill the cokernel, and a rank-one packet must
pass two explicit logarithmic-unit first-jet equations before its
determinant support can be singular.

<!-- status-consumer: PF2K1TN1 521fb57f7e6abc1f -->

The stronger
[`affine strict-log-étale resolution theorem`](AFFINE_KELLER_STRICT_LOG_ETALE_RESOLUTION.md)
applies before any singularity classification: every embedded resolution of
this affine conductor curve pulls back strictly étale through the Keller
map.  Hence the delta and conductor lengths computed here are ordinary
curve/fiber invariants, never direct finite-length summands of the relative
logarithmic cotangent module.  Only their missing boundary attachments can
carry such a class.

<!-- status-consumer: PAER1 60eb24b2232d159e -->

## Source

- Ernst Kunz, *Introduction to Plane Algebraic Curves*, Chapter 17, for the
  conductor/value-semigroup and adjoint-differential description of reduced
  plane curve singularities.

## Reproduction

```bash
.venv/bin/python scripts/verify_f2_affine_target_k1_conductor_conservation.py
```

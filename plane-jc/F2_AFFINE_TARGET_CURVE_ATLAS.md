# F2 affine target-curve and collision atlas

> **Status.**  Let `C` be any irreducible component of the affine
> nonproperness set of a hypothetical F2 `(75,125)` Keller map.  Its
> normalization has a polynomial parametrization
> `(p(t),q(t))` with exact degrees `(3k,5k)` for one of the 24 integers
> `1<=k<=24`.  The projective closure of `C` has degree exactly `5k`, its
> unique point at infinity is `[0:1:0]`, and a normalized implicit equation
> has top homogeneous part `P^(5k)`.  The affine curve `C` is singular.
> Equivalently, the two divided differences of `p` and `q` have a common
> zero: every candidate lies on the collision/critical ideal, not on the
> embedding open.  This reduces the target search from arbitrary degree-124
> parametrizations to 24 finite coefficient charts of at most 194 raw
> coefficients.  It does not identify which chart or component occurs and
> does not construct a Keller map.

The arithmetic reduction and collision/critical examples are replayed by
[`verify_f2_affine_target_curve_atlas.py`](../scripts/verify_f2_affine_target_curve_atlas.py).

## 1. Normalization and the degree ratio

Work first over `C`.  Let

\[
 F=(P,Q):\mathbb A^2\longrightarrow\mathbb A^2,
 \qquad (\deg P,\deg Q)=(75,125).
\]

Jelonek--Lasoń give a polynomial parametrization of every irreducible
component `C` of the nonproperness set with degree at most
`max(75,125)-1=124`.  Nguyen Van Chau's dicritical-series theorem gives a
polynomial parametrization `h=(h_1,h_2)` of the same component satisfying

\[
 \frac{\deg h_1}{\deg h_2}
 =\frac{\deg P}{\deg Q}=\frac35.                \tag{1.1}
\]

The one-place-at-infinity theorem and rational parametrizability identify
the affine normalization with `A^1`.  Write its normalization map as

\[
 \nu(t)=(p(t),q(t)): \mathbb A^1\longrightarrow C. \tag{1.2}
\]

Every dominant polynomial parametrization of `C` factors uniquely through
`nu` by a nonconstant polynomial self-map of `A^1`.  Thus if
`h=nu o r`, then

\[
 \deg h_i=(\deg \nu_i)(\deg r).                 \tag{1.3}
\]

Equation (1.1) therefore holds for the normalization degrees themselves.
The bounded parametrization also factors through `nu`, so its degree bounds
the degrees of `p,q`.  Consequently

\[
 \boxed{(\deg p,\deg q)=(3k,5k),\qquad1\le k\le24.} \tag{1.4}
\]

The factorization argument is invariant under extension of the
characteristic-zero ground field to `C`, so the conclusion descends by the
Lefschetz principle.

## 2. Exact projective degree and the infinity form

The normalization extends to a birational morphism

\[
 \bar\nu:\mathbb P^1\longrightarrow\bar C\subset\mathbb P^2.
\]

Homogenizing both coordinates to degree `5k` gives

\[
 [S:T]\longmapsto
 [T^{2k}p_h(S,T):q_h(S,T):T^{5k}].              \tag{2.1}
\]

There is no base point, and the pullback of a generic target line has degree
`5k`.  Since `bar nu` is birational,

\[
 \boxed{\deg\bar C=5k\in\{5,10,\ldots,120\}.}   \tag{2.2}
\]

At `T=0`, (2.1) lands at `[0:1:0]`.  Chau's theorem says that the complete
nonproperness curve has only this point at infinity.  If `G(P,Q,W)` is an
irreducible homogeneous equation of degree `5k`, its restriction to `W=0`
is therefore a binary form supported at the single point `P=0`.  After
scaling `G`,

\[
 G(P,Q,0)=P^{5k}.                               \tag{2.3}
\]

Thus an affine equation may be normalized as

\[
 g(P,Q)=P^{5k}
   +\sum_{i+j\le5k-1}c_{ij}P^iQ^j.              \tag{2.4}
\]

This is stronger than the previous undirected bound `deg C<=124`.

## 3. Singularity and the universal collision ideal

Chau's affine-line obstruction says that a nonsingular polynomial map of
the plane cannot have a nonproperness component isomorphic to `A^1`.  Since
the normalization of `C` is `A^1`, the finite birational map `nu` cannot be
an isomorphism.  Hence `C` has at least one affine singular point.

For independent variables `s,t`, define the divided differences

\[
 \Delta_p(s,t)=\frac{p(s)-p(t)}{s-t},\qquad
 \Delta_q(s,t)=\frac{q(s)-q(t)}{s-t}.            \tag{3.1}
\]

They are polynomials, with diagonal restrictions

\[
 \Delta_p(t,t)=p'(t),\qquad
 \Delta_q(t,t)=q'(t).                            \tag{3.2}
\]

At a singular point of `C`, either the normalization has at least two
preimages, giving an off-diagonal common zero of (3.1), or it has one
preimage and is nonimmersive there, giving a diagonal common zero by (3.2).
Therefore

\[
 \boxed{(\Delta_p,\Delta_q)\subset k[s,t]
        \text{ is a proper ideal}.}              \tag{3.3}
\]

This is the first exact coefficient condition on the unknown purity target
curve.  It is precisely a collision/critical condition.  A candidate with
unit divided-difference ideal would be a finite injective immersion and hence
a closed embedding of `A^1`; it cannot occur for a Keller nonproperness
component.

After translating one singular point to the target origin, (2.4) also obeys

\[
 c_{00}=c_{10}=c_{01}=0.                         \tag{3.4}
\]

For fixed `k`, the raw normalization-parametrization chart has only

\[
 (3k+1)+(5k+1)=8k+2                             \tag{3.5}
\]

coefficients before reparametrization and target normalizations; its maximum
is `194`, not the thousands of coefficients in an unrestricted implicit
degree-124 equation.

## 4. What this does and does not close

Combining this theorem with the affine-purity frontier gives a finite first
target for the missing row:

1. enumerate `k=1,...,24`;
2. normalize `deg(p,q)=(3k,5k)`;
3. impose the proper ideal `(Delta_p,Delta_q)` and irreducibility/birationality;
4. form the implicit equation with top term `P^(5k)`;
5. factor `g(P,Q)` on the affine source and compute all boundary valuations;
6. attach the resulting component to the common resolution and evaluate its
   logarithmic Chern module.

The theorem supplies Steps 1--3 as necessary geometry.  It does not select
`k`, solve the collision ideal, factor the pullback, determine `(e_i,f_i)`,
or prove that the new component cannot cancel the root contribution `27`.
Thus `(75,125)` and `JC(2)` remain open.

The subsequent
[`k=1 collision theorem`](F2_AFFINE_TARGET_K1_COLLISION.md) solves the first
chart's divided-difference equations by one quartic and identifies its
generic conductor as four affine nodes plus the fixed `(2,5)` infinity cusp.

<!-- status-consumer: PF2K1C1 358a6ba820e8b2f1 -->

The
[`puncture-attachment theorem`](F2_AFFINE_PURITY_PUNCTURE_ATTACHMENT.md)
also shows that every chart meets the target divisor `(5,2)` with contact
`k` and leading residue `p_lead^5/(-q_lead)^3`.

<!-- status-consumer: PF2PPA1 b24c4d80c2f8230e -->

## Sources

- Nguyen Van Chau,
  [*Non-proper value set and the Jacobian condition*](https://arxiv.org/abs/math/0305088),
  for the one-point-at-infinity theorem.
- Nguyen Van Chau,
  [*Plane Jacobian conjecture for simple polynomials*](https://arxiv.org/abs/0711.3894),
  Theorems 3 and 4, for the coordinate-degree ratio and affine-line
  obstruction.
- Z. Jelonek and M. Lasoń,
  [*Quantitative properties of the non-properness set of a polynomial map*](https://arxiv.org/abs/1411.5011),
  for the degree-`124` parametrization bound.

## Reproduction

```bash
.venv/bin/python scripts/verify_f2_affine_target_curve_atlas.py
```

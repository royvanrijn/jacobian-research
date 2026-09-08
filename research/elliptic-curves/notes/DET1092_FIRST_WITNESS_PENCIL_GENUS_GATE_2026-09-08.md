# First302 witness: the whole next RR pencil has genus two over Q

## Result and scope

**Retrospective verified application and new deduction.** Every rational
affine member of the first-centre RR pencil through the historical first
exceptional witness has geometric genus2. There is no rationally defined
genus-zero or genus-one member in that pencil. Its single projective member
at infinity is `C_min+F_0`, the old bisection plus the entire302 fibre.

This closes the low-genus explanation in the **whole fixed next RR system
through this point**, rather than testing one member at a time. It does not
exclude another centre, another parity class, another generic translate of
the exceptional point, or a larger/different divisor system. The full
rank-gain mechanism remains unknown.

This is not an oracle-free construction. The previously saved successful
chart coordinate `z0=-1714/2373` selects the diagnostic pencil. No new point
search or member selection is performed, and the exceptional point's saved
coordinates are not used.

Exact data and replay:

- [sextic and modular exposure certificate](../../artifacts/generated-results/elliptic-curves/det1092_first_witness_pencil_genus_gate_v1.json);
- [independent proof certificate](../../artifacts/generated-results/elliptic-curves/det1092_first_witness_pencil_genus_gate_replay_v1.json);
- [audit](../cas/audit_det1092_first_witness_pencil.sage) and
  [checker](../cas/verify_det1092_first_witness_pencil.sage).

## Uniform equation, without member-by-member searches

**Verified application.** The prior full-chart bridge gives the exact rational
`u0=u(z0)`. Write the pencil as

\[
f_i(t,v)=B_i(t)+(u_0+vt)A_i(t).
\]

In the standard short model let `C=P_w=(c_x,c_y)` and `a=-c4/48`, and put
`h^2=den(x_C)`. The checker verifies the identity

\[
f_2^4H\left(t,-f_1/f_2+a_1/2\right)
=h^6\,c\,q(t,v),
\]

where `H` is the already checked slope quartic, `c` is a nonzero rational
constant retained explicitly, and `q` is a primitive integral polynomial.
Its degrees are at most6 in `t` and4 in `v`. The factors `f2^4,h^6` are
squares; retaining `c` prevents a normalization from silently changing the
rational squareclass. This is a uniform branch equation, not a separately
factored example. Also `f2(0,v)` is a fixed nonzero constant, so no finite
rational `v` makes the displayed map identically undefined.

Let `d(v)` be the universal degree-six discriminant of `q(t,v)` in `t`.
It is an integer polynomial. Its nonvanishing is the squarefreeness condition
for the homogeneous binary sextic, including the infinity chart. A simple
root at infinity (affine degree5) does not change genus2.

## Why the discriminant has degree exactly28

**New deduction with exact coefficient hypotheses checked.** Write

\[
q(t,v)=\sum_{j=0}^4 t^jv^j q_j(t).
\]

The support condition is checked coefficientwise. The top term is
`t^4 v^4 q4(t)`, with `q4` a squarefree quadratic and `q4(0)!=0`.
The quartic `sum_j s^j q_j(0)` is also squarefree of degree4.

Put `epsilon=1/v` and `q_epsilon=epsilon^4 q(t,1/epsilon)`.
Two roots tend to the distinct nonzero roots of `q4`; four roots have the
form `epsilon*s_i+O(epsilon^2)` for distinct `s_i`. Thus the product of
squared root differences has valuation `2*binom(4,2)=12`. The leading
sextic coefficient is an epsilon-unit. Since the sextic discriminant is
homogeneous of degree10 in its coefficients,

\[
\operatorname{disc}(q_\epsilon)=\epsilon^{40}d(1/\epsilon),
\qquad \deg d=40-12=28.
\]

This proves the degree over Q; it is not inferred from a modular degree.
The discriminant/resultant convention agrees with the
[Sage polynomial documentation](https://doc.sagemath.org/html/en/reference/polynomial_rings/sage/rings/polynomial/polynomial_element.html).

## The exact modulo191 obstruction

**Verified application.** The discriminant modulo191, in ascending order,
has coefficients

```text
71,92,169,118,63,95,176,88,76,48,18,141,155,139,2,
27,141,173,139,155,62,9,21,2,90,38,154,119,180.
```

The leading coefficient is nonzero and

\[
\gcd\bigl(d(v)\bmod191,\ v^{191}-v\bigr)=1.
\]

An independent checker reconstructs this degree28 polynomial from29
evaluations of explicit11-by11 Sylvester determinants over `F191`.
It does not trust the constructor's modular resultant call.

Because the global discriminant is integral and its leading coefficient is
a191-adic unit, every rational root would be191-integral. Reduction would
then give a root modulo191, which the gcd excludes. Therefore

\[
d(v)\ne0\qquad\text{for every }v\in\mathbf Q.
\]

**Established literature applied to verified data.** Every such binary
sextic has six simple geometric branch points, so its normalized double
cover has genus2 by
[Riemann--Hurwitz](https://stacks.math.columbia.edu/tag/0C1B).
A rational or genus-one curve cannot map nonconstantly to this genus2
normalization. This is a genuine obstruction to a low-genus base change
through these particular members, not merely a failure to parametrize them.

## What remains, and computation limits

The projective pencil has only one member outside the affine chart:
`tA`, namely `C_min+F_0`. Its vertical302 component is not a point
construction. The preceding
[reducible-locus proof](DET1092_RR_NET_REDUCIBLE_LOCUS_2026-09-08.md)
already excluded known-section pairs through the independent witness.

The useful conclusion is to stop seeking a rational or elliptic member
inside this fixed pencil. A constructive explanation must change the centre,
the representative being diagnosed, or the divisor system, or explain the
arithmetic on a genuinely higher-genus base. This note selects none of those
alternatives by using the exceptional point as an execution oracle.

Limits: one exact sextic and at most18 modular discriminants at fixed primes
through197. The certificate retains all17 checks through the successful191
obstruction, including degree drops and roots at earlier primes. No global
discriminant factorization, rational point search, parameter sweep, class-group
computation, or V3 change occurred. Both jobs finished in seconds.

```sh
sage -python research/elliptic-curves/cas/verify_det1092_first_witness_pencil.sage
```

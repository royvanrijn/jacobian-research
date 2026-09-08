# Singular members of the first-centre RR net: a genus-nine gate

## Result and boundary

**Verified application and new deduction.** On the explicitly specified open
part of the [first-centre RR net](DET1092_FIRST_UNLOCK_RR_NET_2026-09-08.md),
singular members are parametrized by the curve of halves `2Q=P_w(t)`.
This curve is geometrically connected of degree4 over the parameter line and
has genus9. Its map to the singular-member parameter locus is birational onto
its image. Consequently that component of the singular-member locus has
genus9, not genus0 or1.

At `t=0` there is an exact local obstruction to halving: the slope quartic
reduces modulo17 to `m^4+2m^2+8m+14`, with no root. This excludes a rational
singularity above zero in the stated open part. It does not exclude a member
that is singular above a different parameter and passes smoothly through the
first302 point.

**Unknown.** Exceptional members, reducible members, multiple tangencies,
and singularities over excluded fibres remain unclassified. The genus-nine
gate does not rule out isolated rational singular members or special rational
normalizations. It does not explain the independent302 point. In particular,
the halving curve here describes **nodes of the multisections**, not the
exceptional Mordell--Weil directions as literal halves of MW17.

The [generic-input certificate](../../artifacts/generated-results/elliptic-curves/det1092_rr_net_halving_gate_v1.json)
and [independent replay](../../artifacts/generated-results/elliptic-curves/det1092_rr_net_halving_gate_replay_v1.json)
contain all rational-function coefficients and the local obstruction.

## The explicit halving curve

**Verified application.** Keep the historical generic norm10 word `w` from
the RR-net certificate, and write `C=P_w`. Put the generic parent in the
standard short model

\[
Y^2=X^3+aX+b,\qquad
X=x+b_2/12,\quad Y=y+(a_1x+a_3)/2,
\quad a=-c_4/48,\quad b=-c_6/864.
\]

Let `C=(c_x,c_y)` in this model. The line of slope `m` through `-C` has
two residual intersections whose discriminant is

\[
H(t,m)=m^4-6c_xm^2-8c_ym-3c_x^2-4a.
\]

The curve `H=0` has the exact map

\[
X_Q=(m^2-c_x)/2,\qquad Y_Q=m(X_Q-c_x)-c_y.
\]

The checker verifies, modulo `H`, both the curve equation and the tangent
identity `2mY_Q=3X_Q^2+a`, and then the duplication identities `2Q=C`.
Conversely, on the finite nondegenerate chart, the tangent at any half of `C`
passes through `-C` and gives a root of `H`.

The important coefficientwise identity is

\[
\boxed{\operatorname{disc}_m H=256\Delta_E(t).}
\]

This is checked exactly for the supplied parent, not inferred from numerical
root clustering or finite-field patterns.

## Why its genus is nine

**Verified application of established surface theory.** The parent is the
proved24-I1 K3 with smooth infinity. Every nonzero geometric section has
height `4+2(P.O)`, so there is no nonzero torsion and no generic half of a
height10 section: such a half would have height `10/4 < 4`.

These facts also prove geometric connectedness of the degree-four halving
cover. A degree-one component would give a generic half. A degree-two
component with conjugate halves `Q,Q'` would give a rational nonzero
2-torsion point

\[
Q+Q'-C=Q'-Q.
\]

Every partition of4 into multiple components has a part of size1 or2.
Thus no disconnected decomposition is possible. This argument uses the
positive height bound and does not require a new saturation computation.

**Established literature.** Multiplication by2 is etale on a smooth elliptic
curve in characteristic zero, with the same tangent-space argument applying
to its smooth family; see [Stacks, Lemma39.9.9](https://stacks.math.columbia.edu/tag/0BFH).
The height formula is already applied and certified in the
[parent note](CURVE302_RECOVERED_MW17_PARENT_2026-09-07.md), using
[Schuett--Shioda](https://arxiv.org/abs/0907.0298).

**Verified application.** The checker proves that `c_x,c_y` have no poles
at the24 simple discriminant roots. There the monic quartic discriminant has
valuation1, giving one simple ramification point and profile `(2,1,1)`.
There is no ramification over any smooth fibre, including apparent poles of
the slope chart and infinity. Consequently the total ramification is24.

**Established literature applied to these data.**
[Riemann--Hurwitz](https://stacks.math.columbia.edu/tag/0C1B) gives

\[
2g-2=-2\cdot4+24=16,\qquad g=9.
\]

## Exact map to singular members

**New deduction, with rational identities independently checked.** Write the
net as `B+(u+vt)A`. In literal coordinates the slope is `m-a1/2`, so define

\[
r(t,m)=\frac{N(t,m)}{L(t,m)}
=-\frac{B_1+(m-a_1/2)B_2}{A_1+(m-a_1/2)A_2}.
\]

`N,L` are explicit linear polynomials in `m` with nonzero Mobius determinant.
On `H=0`, the branch derivative is `dm/dt=-H_t/H_m`. Therefore set

\[
v=r_t-r_mH_t/H_m,\qquad u=r-tv.
\]

These rational maps are defined where the parent fibre is smooth, coefficients
are finite, and `L H_m r_m != 0`. Substituting `r=u+vt` into the local
residual-discriminant equation gives both a zero and a zero first derivative
at the selected branch point. This is the singular-member incidence condition.
In local square coordinates `s^2=D(t)`, it is `D=D'=0`.

The map is not secretly a high-degree parametrization of a rational curve.
Along `H`, write `d/dt` for the function-field derivation. Then `dr/dt=v`, so

\[
du=-t\,dv.
\]

If `v` were constant, `r-vt` would be constant, making `m` rational in `t`
by the inverse Mobius map, contradicting the degree-four connected cover.
Thus `dv != 0`. In characteristic zero the differential ratio belongs to the
image function field `Q(u,v)`, yielding

\[
t=-du/dv,\qquad r=u+vt,
\]

and then `m` by the inverse Mobius map. The map is consequently birational
onto its image. This proves the genus-nine assertion for the dense-open
singular-member component without computing a large resultant.

A nonconstant rational or elliptic parametrization of this component is
impossible by Riemann--Hurwitz. This says nothing about the existence of
individual rational points on it or exceptional fibres of the normalization
map. In particular, multiple tangencies and reducible members need separate
analysis; they are not discarded by this argument.

## Exact local obstruction at302

**Verified application.** The independent replay reduces the monic halving
quartic at zero modulo17 and checks

\[
H(0,m)\bmod17=m^4+2m^2+8m+14
\]

has no root in `F17`. Its coefficients are17-integral and its leading
coefficient is1. A rational root would therefore be17-integral and reduce to
a root. Hence there is no rational half above zero. The missing slope at
infinity cannot supply one: the centre is finite and nonzero, and a vertical
tangent doubles to the identity, not to `C(0)`.

This is fully compatible with the successful pointed search, which found a
**nonzero square value** of a quartic, not a zero of its branch polynomial.
The distinction is arithmetic, not merely a choice of chart.

## Limits, next unresolved pieces and replay

The constructor used one quartic discriminant and a fixed list of at most18
primes through197; the first prime17 supplied the obstruction. The independent
checker uses direct coefficient and finite-field identities. No general
function-field genus algorithm, class-group calculation, parameter sweep,
point search, or V3 pilot modification was performed. Jobs finished in seconds.

The subsequent [reducible-locus calculation](DET1092_RR_NET_REDUCIBLE_LOCUS_2026-09-08.md)
classifies all reducible members:22 pairs of known generic sections and the
old-bisection-plus-fibre pencil. Irreducible curves with lower-genus
normalizations remain unresolved. No unbounded discriminant-locus search
is authorized by these calculations.

```sh
sage -python research/elliptic-curves/cas/verify_det1092_rr_net_halving_gate.sage
```

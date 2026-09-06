# Fixing the descent cubic requires genus28 or31 in the actual families

The fixed-incidence twist control cannot be transferred to the five relevant
K3 families by a rational or elliptic parameter change which keeps the
entire cubic descent field constant. The minimal constant-cubic isomorphism
cover has **degree6 and genus31** for the three R17 families, and **degree6
and genus28** for the two MW16 families. These are equation-only geometric
computations, independent of exceptional points and rank labels.

There is also a direct paired obstruction: **all eight frozen high/low pairs
have different discriminant squareclasses**. Their cubic fields are therefore
nonisomorphic. They cannot be treated as different solubility states of one
unchanged cubic class group.

This closes a direct transfer route. It does not exclude simultaneous
solubility mechanisms in which the cubic field varies, nor low-degree
constructions which transport less than the full two-torsion representation.

## Exact geometric panel

The [input](../../artifacts/generated-results/elliptic-curves/rank_jump_fixed_field_transfer_geometry_inputs_v1.json)
contains only the five compact A(t),B(t) arrays and one previously retained
smooth parameter per family. That parameter is used only to prove generic
cubic irreducibility by reduction modulo a good prime. No point, Kummer
class or rank label is projected.

For f_t(x)=x³+A(t)x+B(t), the
[certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_fixed_field_transfer_geometry_v1.json)
verifies degree A≤8, degree B≤12, discriminant degree24, gcd(c4,D)=1,
complete squarefree multiplicities and good reduction at infinity.
Let r be the number of geometric discriminant roots of odd multiplicity.

| Family | r | Genus of degree3 root cover | Genus of degree2 resolvent cover | Genus of degree6 constant-cubic carrier |
|---|---:|---:|---:|---:|
| 074d9 | 24 | 10 | 11 | 31 |
| 103b2 | 24 | 10 | 11 | 31 |
| 11952 | 24 | 10 | 11 | 31 |
| a1-fibration-01 | 22 | 9 | 10 | 28 |
| a1-fibration-05 | 22 | 9 | 10 | 28 |

The R17 discriminant divisors are squarefree. Each MW16 divisor has one
double root and22 simple roots. These are family-wide geometric facts, not
features separating high and low rational specializations.

## Why the monodromy and genus calculations are exact

The specialized cubic is irreducible modulo a good prime19,31,31,7,23
respectively. Since it is monic in x, this proves irreducibility over Q(t).
The discriminant has an odd-order zero, so its squareclass is nonconstant
even over Qbar(t). The arithmetic Galois group is S3. The geometric group
is a normal subgroup of S3 and contains an odd permutation; it must therefore
also be S3. In particular the splitting cover is geometrically connected.

At a finite discriminant root, c4 is a unit, so exactly two cubic roots
coalesce and the third remains simple. The discriminant of the local
quadratic factor has the same valuation m as D. Over Qbar((u)), its unit
part is a square; adjoining its roots requires sqrt(u^m). Thus odd m gives
a transposition and even m gives no ramification. Infinity is unramified:
after x=t⁴X and y=t⁶Y, the discriminant at t=infinity is the verified
nonzero leading degree24 coefficient.

Riemann–Hurwitz now gives

\[
 2g_3-2=-6+r,\qquad
 2g_2-2=-4+r,\qquad
 2g_6-2=-12+3r.
\]

The contributions per branch point are1 in the three-letter action,
1 in the sign action, and3 in the regular six-letter action. No branch
point has been omitted, including the MW16 double root and infinity.

## The explicit isomorphism carrier: a condition on t

Fix any separable cubic field K0=Q(eta) with eta³+a eta+b=0. Every trace-zero
element of K0 has the form

\[
 q=v\eta+w(\eta^2+2a/3).
\]

Its characteristic polynomial is x³+A0(v,w)x+B0(v,w), where

\[
 A_0=av^2+3bvw-\frac{a^2}{3}w^2,
\]
\[
 B_0=bv^3-\frac{2a^2}{3}v^2w-abvw^2
       -\left(\frac{2a^3}{27}+b^2\right)w^3.
\]

Consequently, at a smooth rational t,

\[
 \boxed{K_t\simeq K_0
 \quad\Longleftrightarrow\quad
 \exists(v,w)\in\mathbb Q^2:
 A(t)=A_0(v,w),\ B(t)=B_0(v,w).}
\tag{1}
\]

Indeed, an isomorphism sends the trace-zero cubic generator to q. Conversely,
the two coefficient identities give its characteristic polynomial f_t.
Nonzero discriminant ensures q is a primitive generator and gives an
isomorphism. The
[carrier artifact](../../artifacts/generated-results/elliptic-curves/rank_jump_cubic_isomorphism_carrier_v1.json)
verifies the universal characteristic-polynomial identity and supplies exact
coefficients for the five fixed fields obtained at the retained parameters.
Each has the identity lift v=1,w=0 at its defining parameter; this uses no
elliptic point.

The determinant of the basis(1,q,q²) relative to(1,eta,eta²) is

\[
 J(v,w)=v^3+avw^2+bw^3,
\]

and the exact discriminant identity is

\[
 -4A_0^3-27B_0^2=(-4a^3-27b^2)J(v,w)^2.
\tag{2}
\]

Thus the familiar discriminant-square condition is a necessary projection
of (1), not a sufficient field-isomorphism condition. Neither condition
asserts rational solubility of an exceptional Mordell–Weil class.

Over the smooth t-line, (1) is the scheme of isomorphisms between two cubic
étale algebras. Its geometric fibre has6 points. Geometrically it is the
S3 splitting cover, twisted over Q by K0. Its smooth projective normalization
therefore has the genus g6 in the table, independently of the chosen K0.

## Minimality and finiteness

Let C be a smooth projective geometrically integral curve over Q, with a
nonconstant t-map to P1, such that the pulled-back cubic algebra becomes
the constant algebra K0 tensor Q(C). After extending constants to Qbar,
the cubic splits. Therefore Qbar(C) contains its geometric splitting field,
and C geometrically dominates the degree6 splitting cover.

Writing n=degree(C→P1), we obtain

\[
 6\mid n,\qquad
 g(C)-1\ge\frac n6(g_6-1).
\]

In particular, genus(C)≥31 for these R17 families and≥28 for these MW16
families. The isomorphism cover itself attains degree6 and the stated genus.
No nonconstant rational substitution t=R(s), and no genus-one parameter
curve, can keep the full cubic algebra constant. This remains true after
any finite extension of the constant number field.

There is also an arithmetic consequence: for any one fixed cubic field K0,
only **finitely many rational smooth parameters t** in each of these families
can satisfy K_t≅K0. The isomorphism carrier has genus>1, so this follows from
Faltings's theorem; see
[Milne, IV§1](https://www.jmilne.org/math/CourseNotes/AV.pdf).
This is a finiteness theorem, not an effective enumeration or a numerical
bound. Its rational points are not computed here.

## The original matched pairs already have different fields

The paired arithmetic reads only frozen equations and the original pair
indices. Every cubic is independently certified irreducible with S3 Galois
group. For each pair it computes the exact rational ratio of cubic
discriminants. A negative ratio or an integer floor-square certificate
proves that ratio nonsquare. No factorization or exceptional point is needed.

| High token | Matched low token | Discriminant ratio | Cubic fields |
|---|---|---|---|
| case-00 | case-01 | negative | nonisomorphic |
| case-02 | case-03 | positive nonsquare | nonisomorphic |
| case-04 | case-05 | negative | nonisomorphic |
| case-06 | case-07 | negative | nonisomorphic |
| case-08 | case-05 | positive nonsquare | nonisomorphic |
| case-09 | case-05 | negative | nonisomorphic |
| case-10 | case-05 | positive nonsquare | nonisomorphic |
| case-11 | case-12 | negative | nonisomorphic |

Tokens retain the [original panel identities](FRESH_RANK27_GOVERNING_AND_CT_COMPARISON.md).
Repeated controls are not independent replications. These tests say why one
cannot transport the same cubic ideal classes directly across the pairs;
they do not show that a particular field invariant discriminates jump size.

## Consequence for the next experiment

The [six-direction twist switch](FIXED_INCIDENCE_SIX_DIRECTION_SOLUBILITY_SWITCH.md)
remains a valid controlled solubility mechanism. What fails is a direct
transfer through a low-genus base change that preserves its entire descent
field inside these actual families. The obstruction is geometric and exists
before exceptional points are supplied.

A useful replacement must allow the cubic field to vary while defining
the proposed additional classes coherently, for example as2-covers over a
family or an auxiliary parameter curve. Before trying to construct such a
block, the next test should bound how many independent classes can extend
over the generic family with the necessary geometric local conditions.
This distinguishes classes already present in a global cover construction
from arithmetic classes which arise only after specialization.

No claim is made that all mechanisms require a constant field. A common
auxiliary curve for several rational points may keep neither E[2] nor the
cubic field constant; the genus bounds above do not apply to that weaker
requirement. No visibility feature or new candidate score is proposed.

## Reproduction

All five30-second workers completed. The
[portable geometric verifier](../../artifacts/generated-results/elliptic-curves/rank_jump_fixed_field_transfer_geometry_verification_v1.json)
checks rational polynomial identities, squarefreeness/coprimality, finite
irreducibility and all Hurwitz calculations. The
[portable carrier verifier](../../artifacts/generated-results/elliptic-curves/rank_jump_cubic_isomorphism_carrier_verification_v1.json)
checks245 rational trace/norm/discriminant identities on grids sufficient
for exact interpolation of each specified carrier, individual cubic S3
witnesses and all eight discriminant-ratio certificates.

```sh
timeout 60 python3 elliptic-curves/rank-jump/verify_fixed_field_transfer_geometry.py check
timeout 30 python3 elliptic-curves/rank-jump/verify_cubic_isomorphism_carrier.py check
```

The protocols are [geometry](FIXED_FIELD_TRANSFER_GEOMETRY_PROTOCOL.json)
and [carrier/pairs](CUBIC_ISOMORPHISM_CARRIER_PROTOCOL.json). No new
parameter search, class-group computation, elliptic point search or active
search modification was performed. Mathematical-status entries are untouched.

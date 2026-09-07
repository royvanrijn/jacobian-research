# Additive-component compatibility closes the arithmetic mixed ranks

Follow-up: the [equal-class Picard test](EQUAL_CLASS_COMPATIBILITY_IS_NOT_SUFFICIENT.md)
proves geometric mixed rank zero for the compatible low-control triple.
Compatibility is therefore not sufficient even for geometric incidence;
the compatible high-control triple remains at rank `0–1`.

**The paired constructions now have exact arithmetic generic ranks: two
for the two-point base and three for the three-point genus-one base.**
Every non-singleton arithmetic character space is zero on both anchors.

This strengthens the earlier arithmetic intervals; their geometric Picard
bounds remain valid. The new proof does not need another surface point count.
It isolates a rational-descent obstruction: a mixed section would require
different anchor Kummer classes to coincide at rational additive fibres.

| Retained MW16-04 anchor | Original generic rank / witness quotient | Two-point base arithmetic rank | Three-point base arithmetic rank |
|---|---|---:|---:|
| `-1647/91` | `16 / 9` | 2 | 3 |
| `-2177/2397` | `16 / observed 0`, censored | 2 | 3 |

These are ranks over the **new function fields**, not bounds on either
production curve. At `u=0`, the guaranteed sections are the already selected
generic points. No new exceptional quotient direction is produced.

## 1. The necessary compatibility theorem

Let \(f(T)=T^3+AT+B\) be irreducible over \(\mathbb Q\), with Galois group
`S3` and \(B\operatorname{disc}(f)\ne0\). Put \(K=\mathbb Q(\theta)\).
Use the fixed-cubic pencil
\[
E_u:\ y^2=x^3+2Au x^2+(A+3Bu+A^2u^2)x+B+ABu^2-B^2u^3.
\]
Choose distinct nonzero rational \(a_1,\ldots,a_m\), with \(f(a_i)\ne0\),
and set \(d(u)=\prod_i(1-a_i u)\), \(m\ge1\).
Write \(E_u^d\) in the conventional twist coordinates, whose roots are
\(d(\theta_i+u\theta_i^2)\).

Define
\[
\begin{aligned}
\gamma&=1-u\theta,\\
\kappa&=1+u\theta+u^2(A+\theta^2),\\
D&=1+Au^2+Bu^3=\gamma\kappa,\\
\lambda(a)&=[f(a)(a-\theta)]\in K^\times/K^{\times2},\\
\lambda_\infty&=[-B\theta].
\end{aligned}
\]
All the \(\lambda\)'s have square norm. Then:

1. \(\operatorname{rank}E_u^d(\mathbb Q(u))\le1\).
2. Positive rank requires **all** \(\lambda(a_i)\) to be equal.
3. If `m` is even, their common value must also equal
   \(\lambda_\infty\).
4. These conditions are necessary, not sufficient.

When \((a_i,y_i)\in E_0(\mathbb Q)\), the finite condition simplifies to
\[
\boxed{[a_1-\theta]=\cdots=[a_m-\theta].}
\]
Thus a product character involving two different anchor Kummer classes has
arithmetic rank zero. This statement concerns **solubility/rational descent**
of the chosen construction; it does not assert that its geometric character
space vanishes.

The conditions are available from the equation and selected branch positions.
In the present experiment those positions are generic-point abscissas, with
no exceptional coordinates. Their availability does not make this a
validated predictor of large specialization rank.

## 2. Rational additive fibres kill the constant ambiguity

At each rational root `u=1/a_i` of `d`, the twist has fibre type `I0*`.
Its geometric component group is \((\mathbb Z/2)^2\), with its three nonzero
elements permuted by the constant `S3` two-torsion action. It therefore
has no nonzero rational element.

Equivalently, the three nonidentity outer components are labelled by the
roots of the reduced cubic. Its residue algebra is still `K`:
the determinant of the basis change from
\((1,\theta,\theta^2)\) to
\((1,\theta+u\theta^2,(\theta+u\theta^2)^2)\) is `D(u)`,
and \(D(1/a_i)=f(a_i)/a_i^3\ne0\).
If `m` is even, infinity is also `I0*`, with residue roots proportional
to \(\theta_i^2\). The corresponding basis determinant is `B != 0`.

A rational section therefore meets the identity component at each of
these rational additive places. The identity-component group over
\(\mathbb Q((v))\) is 2-divisible: its reduction is the additive group
\(\mathbb G_a(\mathbb Q)\), and its formal kernel is 2-divisible by Hensel
lifting. Consequently
\[
E_u^d(\mathbb Q((v)))/2E_u^d(\mathbb Q((v)))=0
\]
at every such place.

The component-group facts and the injection of prime-to-characteristic
torsion into an additive component group are standard; see
[Schütt–Shioda, section 7, especially Lemma 7.8](https://arxiv.org/pdf/0907.0298).
Here the residue field has characteristic zero, and its `S3` action has
zero invariants. The places in this argument are **places of the parameter
line**, with residue field \(\mathbb Q\), not the bad rational primes of a
specialized production curve.

## 3. Only one nonconstant Kummer pattern is possible

For \(P\in E_u^d(\mathbb Q(u))\), its Kummer class is represented by
\[
\beta_P=x(P)-d(\theta+u\theta^2)\in K(u)^\times/K(u)^{\times2}.
\]
The usual norm-kernel two-descent identifies this with an injective map
from \(E_u^d(\mathbb Q(u))/2\); a square representative means divisibility
by two over \(\mathbb Q(u)\), not merely over the splitting field.
At every good finite place in `u` its valuation is even. Indeed, a pole of
the elliptic x-coordinate has even order; near a smooth two-torsion point,
the equation gives even valuation of the corresponding root difference.

At the roots of `d` the entire local Kummer image is zero by the preceding
argument, so these places also contribute no odd valuation.

Over `K`, the remaining bad polynomial factors as `D=gamma*kappa`.
The quadratic `kappa` is irreducible: its discriminant
\(-3\theta^2-4A\) cuts out the quadratic splitting-field extension of the
non-Galois cubic `K`.

At `gamma=0`, the root labelled by `theta` is the **single**, noncolliding
root of the nodal cubic. The other two roots coalesce. Thus
`x(P)-d(theta+u*theta^2)` again has even valuation: it is a unit at the
node, has even valuation near the single smooth two-torsion point, and
has even pole order near zero. The only possible odd finite divisor is
therefore `kappa`.

Unique factorization in `K[u]` proves
\[
\boxed{[\beta_P]=[\beta_0\kappa^e],\qquad
       \beta_0\in K^\times/K^{\times2},\quad e\in\{0,1\}.}
\]
Also \(N(\kappa)=D^2\), so \(\beta_0\) has square norm.
In a splitting field, the three valuation rows of `kappa` are
\[
\begin{pmatrix}0&1&1\\1&0&1\\1&1&0\end{pmatrix};
\]
the diagonal entries correspond to the single roots just described.
This is a divisor statement, not an inference from finitely many local tests.

Now evaluate at a rational additive place `u=1/a`. The local class is zero,
and the universal identity
\[
a^2\kappa(1/a)(a-\theta)=f(a)
\]
gives
\[
[\beta_0]=\lambda(a)^e.
\]
If `e=0`, this forces \(\beta_0\) to be a square, hence \(P\in2E_u^d(\mathbb Q(u))\).
If `e=1`, every additive place forces the **same** constant class
\(\beta_0=\lambda(a_i)\). At an additive infinity fibre, the leading
coefficient of `kappa` has squareclass `[-B*theta]`, since
\[
\theta^2(A+\theta^2)=-B\theta.
\]
This proves the extra even-degree condition.

The map from \(E_u^d(\mathbb Q(u))/2\) to the bit `e` is injective.
The group is finitely generated and has no rational two-torsion, so its
rank is at most one. If the required constant classes disagree, this
quotient is zero and its rank is zero. This proves all parts of the theorem.

## 4. Exact application to the paired controls

The [frozen protocol](COMPONENT_KUMMER_GATE_PROTOCOL.json) tests only the
three previously selected generic points, using odd primes at most 97.
For a prime ideal \((p,\theta-r)\), a nonsquare value of
\((a_i-r)(a_j-r)\) proves that \([a_i-\theta]\ne[a_j-\theta]\).
All displayed residues are nonzero and `p` avoids the cubic discriminant.

| Anchor | Classes distinguished | One exact prime/root witness |
|---|---|---|
| high | `0 vs 1` | `p=53, r=10` |
| high | `0 vs 2` and `1 vs 2` | `p=37, r=9` |
| low | `0 vs 1` and `0 vs 2` | `p=19, r=5` |
| low | `1 vs 2` | `p=41, r=10` |

The [certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_component_kummer_gate_v1.json)
retains all 27 tested prime/root characters and the exact witnesses.
Their combined signatures have rank three on each anchor. Hence every
mask `3,5,6,7` contains incompatible classes, and all four arithmetic
mixed ranks are exactly zero.

The singleton ranks are one, and the original arithmetic generic rank is
zero because `B` is nonsquare. Character decomposition consequently gives
the exact ranks two and three in the opening table.

More generally, for **any number** of distinct point abscissas with pairwise
different anchor Kummer classes, adjoining their linear square roots adds
exactly their singleton arithmetic directions:
\[
\operatorname{rank}E_u\bigl(\mathbb Q(u)(\sqrt{1-a_1u},\ldots,\sqrt{1-a_nu})\bigr)
=n+\mathbf1_{\{B\text{ is a rational square}\}}.
\]
Every higher mixed character contains an incompatible pair. This corollary
does not require its higher-degree surfaces to be K3.
It explains why enlarging this particular construction using distinct
classes cannot create bonus arithmetic directions.

## 5. A control that passes compatibility but does not prove a section

On each anchor form, using only the same three generic points,
\[
R_0=P_0,\qquad R_1=P_0+2P_1,\qquad R_2=P_0+2P_2.
\]
These remain three independent rational directions, with index four in
the original three-point lattice. Their Kummer classes, however, all equal
\([P_0]\). The distinction between free rank and mod-two class rank matters.

This equality is verified by explicit cubic-field square roots. If the
chords for `P+Q` and `(P+Q)+Q` are \(\ell_1,\ell_2\), then
\[
\frac{x(P+2Q)-\theta}{x(P)-\theta}
 =\left(\frac{\ell_2(\theta)}{\ell_1(\theta)}\right)^2.
\]
All four resulting square roots and chord coefficients are retained.
An independent Sage group-law calculation and polynomial reduction check
the identities.

The triple made from `R0,R1,R2` passes the **finite additive-component**
gate; its odd degree has no additive condition at infinity.
Its global rational-section existence remains **UNKNOWN**. A good fibre
at infinity and the other good places can still impose solubility
conditions, and compatibility does not solve the two-cover.

Each paired subproduct still fails the additive condition at infinity:
the certificate distinguishes `[P0]` from `[-B*theta]` on both anchors.
Thus this equal-class triple's new function field has arithmetic rank
**three or four**, with only its triple character unresolved. The rank of
its new parameter curve remains UNKNOWN; it is a different curve from
the previously certified positive-rank base.

This control also prevents misusing the criterion as an intrinsic rank
score: such matching classes can be manufactured from generic points on
either control. The next falsifiable question is whether this compatible
triple actually supplies its permitted mixed section, and what that section
specializes to modulo the original generic subgroup.

## 6. Reproduction, limits, and remaining mechanisms

The [independent replay](../../artifacts/generated-results/elliptic-curves/rank_jump_component_kummer_gate_verification_v1.json)
checks the universal identities, residue-algebra determinants, all finite
characters, and all four chord square roots. No new surface counts,
specialized curves, exceptional point searches, or class-group computations
were needed.

```sh
python3 elliptic-curves/rank-jump/component_kummer_gate.py check
sage -python elliptic-curves/rank-jump/verify_component_kummer_gate.py
python3 -m unittest discover -s elliptic-curves/rank-jump -p test_component_kummer_gate.py
```

The resulting mechanism ranking is:

1. **Solubility/rational descent:** additive fibres impose equality of a
   single constant Kummer class. This is an exact obstruction and now closes
   the arithmetic ranks of the previous constructions.
2. **Incidence:** a compatible triple leaves at most one possible mixed
   arithmetic direction. A rational section is still needed to realize it;
   geometric Picard rank alone would not suffice.
3. **Weak explanation:** many distinct promoted classes do not create
   additional mixed arithmetic directions in this pencil, even on an
   infinite rational genus-one base. The construction does not account for
   the original +9-versus-zero difference.
4. **Visibility:** no visibility statistic was tested. Even if a future mixed
   section specializes to a class already in the generic mod-two image,
   that alone would not prove membership in the generic rational span:
   an exceptional quotient vector can be twice another vector.

The remaining implication is now more specific:
\[
\text{compatible additive Kummer data}
\;\stackrel{?}{\Longrightarrow}\;
\text{a rational mixed section}
\;\stackrel{?}{\Longrightarrow}\;
\text{an exceptional specialization direction}.
\]
Agent 1 receives a construction constraint, not a candidate selector.
Existing search protocols, candidates, certificates and mathematical-status
entries are unchanged.

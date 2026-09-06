# Mixed-character bounds on a matched high/zero-gain pair

**Result:** the same two generic-point constructions leave **at most one**
additional mixed-character direction on both retained MW16-04 fibres.
The frozen two-prime experiment does **not** decide whether that direction
exists. The relevant reductions have matching discriminant squareclasses,
so the previous branch-splitting argument cannot be reused to assert rank zero.

This is a retrospective incidence experiment. It uses the first two
generic points on each curve and no exceptional coordinates. It gives
function-field bounds, not rank upper bounds for the production curves.

| Retained MW16-04 fibre | Marked generic rank | Certified witness quotient | Mixed-character geometric rank |
|---|---:|---:|---|
| `-1647/91` | 16 | 9 | 0 or 1 |
| `-2177/2397` | 16 | observed 0, censored | 0 or 1 |

The first three columns inherit the exact subgroup accounting in
[the original paired panel](ANALYSIS.md). The low endpoint is not an
exact rank-16 assertion. No currently changing search output was edited.

## 1. Why this test

The [full branch-splitting test](FULL_BRANCH_SPLITTING_AND_THREE_CHARACTER_GATE.md)
closed one proposed source of extra generic directions for the rank-20 anchor.
The mixed-character possibility remains different: two covers with individual
sections can, in principle, contribute a further section in their product
character. This is a simultaneous construction question, not a measure of
how easily a chart exposes points.

The [protocol](MIXED_CHARACTER_PROTOCOL.json) selects the second original
matched pair and generic-point indices 0 and 1 on each row. Selection used
only a good-reduction preflight: the first pair lacked a second eligible
high-fibre prime at most 19. No Frobenius counts were inspected before freezing.
The selected primes are `13,17` on the high fibre and `17,19` on the low fibre.
Degrees are at most three, with 40 seconds per prime and no point searches.

The [input projection](../../artifacts/generated-results/elliptic-curves/rank_jump_mixed_character_pair_inputs_v1.json)
contains just each short equation and those two generic points, plus provenance.
The high/low choice is retrospective; it is not an independently validated
prospective rank selector.

## 2. The base is rational, and there are four character spaces

For each anchor write \(E_0:y^2=f(x)=x^3+Ax+B\), and call the two
generic points \((a,y_a),(b,y_b)\). They have distinct nonzero abscissas.
Use the fixed-cubic deformation
\[
E_u:\ y^2=x^3+2Au x^2+(A+3Bu+A^2u^2)x+B+ABu^2-B^2u^3,
\quad D(u)=1+Au^2+Bu^3.
\]
Let \(g_a=1-au,\ g_b=1-bu\). The two-cover compositum is
\[
C_{a,b}:\quad z_a^2=g_a,\qquad z_b^2=g_b.
\]

**Correction to the initial progress message:** this base has genus zero,
not genus one. The covers share their branch point at infinity. Their
three inertia vectors are \(10,01,11\), and Riemann–Hurwitz gives
\(2g-2=-8+6=-2\). The rational point above `u=0` makes it a rational curve.
An exact parametrization is
\[
z_a=\frac{at^2-2at+b}{at^2-b},\qquad
z_b=1+t(z_a-1),\qquad u=\frac{1-z_a^2}{a}.
\]
The symbolic replay checks both square equations.

The four character spaces over \(\overline{\mathbb Q}(u)\) are
\(1,g_a,g_b,g_ag_b\). The first has geometric rank one and arithmetic
rank zero since `B` is nonsquare on both anchors. Each singleton has
geometric and arithmetic rank one, by the exact
[linear-twist theorem](LINEAR_TWIST_SOLUBLE_BLOCKS.md); its rational section
uses \(y_a^2=f(a)\), respectively \(y_b^2=f(b)\).

Thus the pair contributes two independent arithmetic directions before
the mixed character is considered. Its generic rank is not the original
MW16 family rank. The reference family has changed.

## 3. The mixed twist is a K3 surface

Put \(d=g_ag_b\). The conventional twist of \(E_u\) has coefficients
\((d a_2,d^2a_4,d^3a_6)\), discriminant
\[
16\delta D^2d^6,\qquad \delta=-4A^3-27B^2,
\]
and \(c_4=16d^2(A^2u^2-9Bu-3A)\).
Its singular fibres are
\[
3I_2+3I_0^*.
\]
The `I2` fibres occur at `D=0`; the `I0*` fibres occur at `1/a,1/b,∞`.
The declared nonvanishing conditions keep them disjoint. At infinity
the discriminant and `c4` valuations are six and two. Hence the minimal
surface is K3, with trivial lattice
\[
U\oplus A_1^3\oplus D_4^3,\qquad \operatorname{rank}=17.
\]
Shioda–Tate gives mixed-character rank \(\rho-17\).
The surface/fibre formulas and rank theorem are standard; see
[Schütt–Shioda, sections 5–6](https://arxiv.org/pdf/0907.0298).
The new application is checked symbolically in
[the geometry verifier](verify_mixed_character_geometry.py).

Before counting, this leaves the geometric rank between zero and three.
After counting, the upper bound is one on each curve.

## 4. Point counts and the five-dimensional complement

The good-prime conditions are
\[
p>3,\qquad p\nmid B\delta\,ab(a-b)f(a)f(b),
\]
including all coefficient denominators. They preserve the displayed
minimal configurations and their resolutions.

At \(q=p^n\), let \(r_q\) be the number of roots of `f` in \(\mathbb F_q\).
The trivial lattice has Frobenius representation
\[
\mathbb Q_\ell(-1)^5\oplus
  \bigl(\mathbb Q_\ell[\{\theta_1,\theta_2,\theta_3\}](-1)\bigr)^4.
\]
Indeed, the three `A1` components form one root permutation module.
Each `D4` contributes a fixed central component and its three outer
components, labelled by the nonzero two-torsion. Fibre and zero section
give the other two fixed classes. Unit twists do not change this root
permutation action on an `I0*` fibre.

Writing \(F_u(x)\) for the untwisted cubic, the resolved count is therefore
\[
\#X(\mathbb F_q)
=q^2+1+q(5+4r_q)+T_q,\qquad
T_q=\sum_{u,x\in\mathbb F_q}\chi(d(u)F_u(x)).
\]
The five-dimensional complementary Frobenius polynomial has the form
\[
(X-\epsilon p)
\left(X^4-sX^3+cX^2-sp^2X+p^4\right).
\]
Three power traces determine the sign and coefficients in every retained
case. A separate companion-matrix replay verifies those three traces.
Factoring the normalized polynomial removes **all** cyclotomic factors
of degree at most five; no root-of-unity test is based on floating point.

| Anchor | p | \(T_p,T_{p^2},T_{p^3}\) | Geometric Picard rank of reduction | NS squareclass used by endpoint |
|---|---:|---|---:|---:|
| high | 13 | `7,543,-4103` | 18 | −1 |
| high | 17 | `-1,239,-6427` | 18 | −1 |
| low | 17 | `-9,-497,21813` | 20 | unused |
| low | 19 | `11,1109,-10525` | 18 | −1 |

For the upper bound in characteristic zero, the cycle-class injection and
good-reduction specialization suffice; see
[van Luijk, Proposition 2.2 and Corollary 2.3](https://arxiv.org/pdf/math/0506416).
Equality for the tabulated reduction ranks uses Tate for elliptic K3 surfaces.

For the three rank-18 reductions, all algebraic classes are fixed over
\(\mathbb F_{p^6}\). Artin–Tate then gives the signed discriminant
squareclass from
\[
-\frac{\det(1-(F_{\mathrm{trans}}/p)^6)}{p^6}.
\]
The determinant is computed once by a resultant and independently by a
four-by-four matrix determinant. All three squareclasses are −1.
The reduction at `p=17` on the low anchor has an additional cyclotomic
factor of order three and is not used for the rank-17 endpoint.

There is **no discriminant mismatch**. Thus neither characteristic-zero
Picard rank 17 nor rank 18 is proved. In particular, matching reduction
squareclasses are not evidence sufficient to assert an extra divisor.

The resulting unconditional intervals are
\[
\begin{array}{c|c}
\text{object}&\text{rank interval on either anchor}\\ \hline
E_u^{g_ag_b}/\overline{\mathbb Q}(u)&[0,1]\\
E_u/\overline{\mathbb Q}(C_{a,b})&[3,4]\\
E_u/\mathbb Q(C_{a,b})&[2,3].
\end{array}
\]
These do not bound the ranks at particular rational points of `C_{a,b}`.
A possible extra mixed section could also specialize into the original
generic subgroup. No exceptional quotient direction is explained yet.

## 5. Independent verification and retained failures

[Complete counts](../../artifacts/generated-results/elliptic-curves/rank_jump_mixed_character_counts_v2.json)
use Sage 10.9 / PARI elliptic cardinalities and explicit Frobenius orbits:
6,984 smooth cardinality calls and seven singular checks.
The [independent NumPy replay](../../artifacts/generated-results/elliptic-curves/rank_jump_mixed_character_verification_v1.json)
directly enumerates character sums for **all 20,056 base parameters**,
including singular fibres, without an elliptic point counter or orbit reuse.
All twelve field arrays agree exactly.

The [analysis](../../artifacts/generated-results/elliptic-curves/rank_jump_mixed_character_v1.json)
and [geometry/bounds](../../artifacts/generated-results/elliptic-curves/rank_jump_mixed_character_geometry_v1.json)
bind their inputs and producers. Initial partial counts remain in
[version 1](../../artifacts/generated-results/elliptic-curves/rank_jump_mixed_character_counts_v1.json).
The [first source and failure](../../artifacts/generated-results/elliptic-curves/rank_jump_mixed_character_attempt_v1.json)
record a singular-fibre handler that incorrectly expected the two coalescing
roots to be separately rational. The corrected handler uses the rational
double root directly, and the same four primes were recounted.
The [count-producer snapshot](../../artifacts/generated-results/elliptic-curves/rank_jump_mixed_character_count_source_v1.json)
also preserves a subsequent Sage-integer JSON serialization failure.
That failure did not change counts. The final analyzer verifies the sealed
producer binding and explicitly converts its integer dimension.

Replay, from the repository root:

```sh
sage -python elliptic-curves/rank-jump/mixed_character.py check
sage -python elliptic-curves/rank-jump/verify_mixed_character_geometry.py
python3 -m unittest discover -s elliptic-curves/rank-jump -p test_mixed_character.py
python3 elliptic-curves/rank-jump/verify_mixed_character.py --case 0 --prime 13
python3 elliptic-curves/rank-jump/verify_mixed_character.py --case 0 --prime 17
python3 elliptic-curves/rank-jump/verify_mixed_character.py --case 1 --prime 17
python3 elliptic-curves/rank-jump/verify_mixed_character.py --case 1 --prime 19
```

No prospective parameter population, scoring policy, search script, or
Agent 1 output was modified.

## 6. Mechanisms and the next discriminating step

1. **Incidence — mixed characters on a small base remain viable, bounded.**
   On this pair they can supply at most one extra geometric direction.
   The exact rank and rational descent of that direction remain unknown.
   A positive mixed-section identity would be meaningful new evidence;
   equal discriminant squareclasses are not a substitute.
2. **Solubility — simultaneous rational lifts remain a separate gate.**
   The pair base is rational and its two singleton directions are explicit.
   Neither fact proves a large specialized quotient, particularly at a
   fibre where their values already lie in the marked generic subgroup.
3. **Weak explanation — this two-point construction does not distinguish
   the retained +9 fibre from its observed-zero control.**
   Both have the same certified intervals. This excludes neither all
   mixed-character constructions nor further invisible points.
4. **Visibility — no new feature was measured.** Finite-field point counts
   determine cohomological constraints; they are not rational point-search
   yield or chart scores.

The next useful small-base experiment is to add the **third generic point**
on this same pair. Three linear characters share infinity and have four
branch points in total, so their compositum really has genus one. Its four
non-singleton characters are still K3 twists; an explicit section or a
certified bound in each would determine whether this construction can carry
a larger block before exceptional points are supplied. Freeze those classes
and bounds before counting. A successful construction must then pass two
further tests: rational descent and independence modulo the original generic
subgroup after specialization. No claim here makes either implication.

# One square condition creates two independent rational directions

There is now a complete small mechanism of simultaneous rational solubility:

\[
 \boxed{2m+2+c^2=n^2
 \ \Longrightarrow\ \text{two sections over one rational double cover}
 \ \Longrightarrow\ \text{arithmetic generic rank }1\to3.}
\]

This is proved for **c=1 and c=3** in the explicit family

\[
 E_{m,c}:\quad y^2=x^3+mx^2-(m+3)x+c^2.
\]

The c=1 case explains the rank-three member of the earlier exact
rational-versus-Sha control. The c=3 case has generic two-division Galois
group S₃, so this construction does not require a cyclic cubic field.
These are small mechanism controls, not new high-rank candidates or an
explanation yet of the production +8,+10,+14 jumps.

The result is about function-field ranks. It does **not** say that every
specialization keeps the three directions. The experiment found and proves
an explicit specialization locus where they become dependent.

## The common arithmetic event

Write F(x)=x³+mx²−(m+3)x+c² and d=2m+2+c². Direct identities give

\[
 F(0)=c^2,\qquad F(-1)=F(2)=d.
\]

Thus R=(0,c) is already rational over Q(m). The two degree-two
multisections at x=−1 and x=2 have the **same** equation y²=d over
the parameter line. A single base change

\[
 B_c:\ n^2=2m+2+c^2,\qquad m=(n^2-c^2-2)/2
\]

makes both sections rational:

\[
 P_- =(-1,n),\qquad P_+=(2,n).
\]

The auxiliary base B_c is a rational curve, with a degree-two map to the
m-line. This is simultaneous **global solubility**, supplied by explicit
sections; no Selmer radical is being promoted to a point. The incidence
calculation below proves that the two sections add independent generic
directions. No point-visibility claim is involved.

The field shared here is the field of definition of the multisections,
Q(m)(√d). It is not a common halving field for independent rational
Mordell–Weil directions; the earlier independent-halving-field results
are unaffected.

## Geometric incidence: exactly three new geometric directions

Both parent elliptic surfaces are rational. Their minimal singular fibres
and Shioda–Tate ranks are:

| Fixed c | Finite parent fibres | Parent infinity | Parent geometric rank | Finite fibres after base change | New infinity | Base-changed geometric rank |
|---|---|---|---:|---|---|---:|
| 1 | 2 II | I₂* | 2 | 4 II | I₄ | 5 |
| 3 | 4 I₁ | I₂* | 2 | 8 I₁ | I₄ | 5 |

For c=1, put q=m²+3m+9. The invariants are

\[
 c_4=16q,\quad\Delta=16q^2.
\]

The two simple roots of q have orders (v(c₄),v(Δ))=(1,2), hence
type II. At infinity, the integral chart x=v⁻²X, y=v⁻³Y gives orders
(2,8), hence I₂*. After m=(n²−3)/2,

\[
 c_4=4(n^4+27),\quad\Delta=(n^4+27)^2.
\]

There are four type-II fibres and an I₄ fibre at infinity.

For c=3 the parent discriminant is

\[
 \Delta=16(m^4-26m^3-117m^2-378m-2079),
\]

while c₄ remains 16q. Exact gcd checks show that this quartic is
squarefree and coprime to c₄. Its pullback after m=(n²−11)/2 is
likewise squarefree of degree eight and coprime to the pulled-back c₄.
The branch point n=0 lies over a smooth parent fibre. At infinity the
orders are again (2,8) before base change and (0,4) afterwards.

Each surface has Euler number 12, hence χ=1 and geometric Picard rank
10. The parent reducible-fibre root lattice is D₆, of rank 6; after
base change it is A₃, of rank 3. Consequently the geometric ranks are
10−2−6=2 and 10−2−3=5. This uses the standard rational-surface and
Shioda–Tate framework in
[Schütt–Shioda, §§6 and 8](https://arxiv.org/pdf/0907.0298).
The generic-rank calculation precedes any claim about a specialization.

## Galois structure: exactly two of the new directions are rational

Let D=c²−2. Two additional identities are

\[
 F(1)=D,\qquad
 F(c^2/2)=\frac{c^2D}{8}\,d.
\]

They give five geometric sections on the base-changed surface:

\[
 R=(0,c),\quad P_-=(-1,n),\quad P_+=(2,n),\quad
 T=(1,\sqrt D),\quad
 U=\left(c^2/2,\frac{cn\sqrt{2D}}4\right).
\]

The certificate proves these five sections independent. The first three
are rational. T has constant-field character χ_D and U has χ_(2D).
For c=1 these are χ_−1 and χ_−2; for c=3 they are χ₇ and χ₁₄.
They are distinct nontrivial quadratic characters. Since the five
sections reach the geometric upper bound, their rational span is the
whole geometric Mordell–Weil vector space. Therefore

\[
 \boxed{\operatorname{rank}E_{m,c}(\mathbf Q(m))=1,
 \qquad\operatorname{rank}E_{(n^2-c^2-2)/2,c}(\mathbf Q(n))=3.}
\]

The deck involution n↦−n fixes R,T and negates P₋,P₊,U.
Thus the new geometric character has dimension 3, with rational dimension
2. The extra geometric direction U does not become rational merely
because d becomes square.

For the independence proof, exact good-prime Kummer fingerprints and
real signs certify rank three for R,P₋,P₊ at the fixed anchors
(c,n)=(1,5) and (3,7). A relation among the generic sections would
specialize to one there. T and U are independently proved nontorsion at
those anchors using two split good-prime reduction bounds over their
constant quadratic fields. Their distinct Galois characters then make
them independent of the three rational sections and of each other.

Both two-division cubics are generically irreducible, witnessed by monic
irreducible specializations modulo 2. For c=1 the discriminant is a
square, giving C₃. For c=3 it is a squarefree nonconstant quartic in m,
giving S₃. The n-family remains S₃ in the c=3 case, since its degree-eight
discriminant is squarefree and the anchor again certifies irreducibility.

## Connection to the strict rational-versus-Sha block

At c=1,n=5, m=11 and the three rational points are
(0,1), (−1,5), (2,5) on the exact rank-three control E₊.
The earlier [strict-class certificate](NORM_LIFTS_CAN_BE_ENTIRELY_SHA.md)
identifies

\[
 \beta_0=\delta(P_-+R),\qquad \beta_1=\delta(P_+).
\]

Thus the two strict classes previously constructed from cubic arithmetic
are represented by the two new section directions, after one generic
correction. On the −1 twist E₀ they are the entire two-dimensional
Sha[2] block. The common-value condition explains their rational
representatives on E₊ without declaring norm lifts or Selmer lifts
automatically rational.

This condition was recognized retrospectively. It is now an
equation-defined construction, with an exact proof, but it has not been
used as a prospective selector on production curves.

## Falsification controls: rational points can collapse together

The original S₃ independence anchor was c=3,n=5, or m=7. Its three
point fingerprints had rank two. The
[bounded halving experiment](SHARED_VALUE_HALVING_PROTOCOL.json) tested
only the seven nonzero subset sums of those prescribed points. It found

\[
 R+P_-=-2P_+.
\]

This is genuine dependence, not a missing saturation direction. Two
independent fingerprints prove that the retained subgroup has rank
exactly two. Its whole curve rank is not determined here.

More strongly, symbolic group-law calculation proves the identity on the
entire locus

\[
 \boxed{n=c+2,\quad m=2c+1:qquad R+P_-+2P_+=0.}
\]

This explains why the first anchor failed. One predeclared replacement
anchor, c=3,n=7, passed the rank-three gate; no further anchors or
parameter sweep were used.

There is an even stronger collapse at c=n=1, m=−1:

\[
 P_+=2R,\qquad P_-=-3R.
\]

The retained image has rank exactly one. R is nontorsion by good-prime
counts 15 and 16 at 11 and 13. Therefore the existence of the three
explicit points alone does not prove a rank-three specialization.
At any fibre, independence must still be checked relative to the
specified generic subgroup. Once the double cover is adopted as the
base, its rank-three subgroup is generic and those two directions must
no longer be counted as an additional jump.

## Reproducible bounded experiment

The [first protocol](SHARED_VALUE_SOLUBLE_BLOCK_PROTOCOL.json) fixed c=1,3,
one initial independence anchor per family, and a collapse control.
It allowed 30 seconds per worker, fingerprint primes at most 199 and
torsion-bound primes at most 43. It performed symbolic identities,
exact finite counts and prescribed point arithmetic—no rank search,
descent, class-group computation or new point enumeration.

The failed S₃ gate led to the single-layer exact halving check, followed
by [one replacement-anchor protocol](SHARED_VALUE_COMPLETION_PROTOCOL.json).
Original sources, the successful first checkpoint and the failed log
remain preserved. A read-only adapter records the replacement's effective
protocol explicitly. The portable input retains the original checkpoint
and failure log, so replay does not need ignored local files.

Artifacts:
[initial evidence](../../artifacts/generated-results/elliptic-curves/rank_jump_shared_value_initial_inputs_v1.json),
[exact halving and dependence](../../artifacts/generated-results/elliptic-curves/rank_jump_shared_value_halving_v1.json),
[complete mechanism certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_shared_value_soluble_block_completion_v1.json).

Replay from the repository root:

```sh
sage -python elliptic-curves/rank-jump/replay_shared_value_block.py check
sage -python elliptic-curves/rank-jump/shared_value_halving.py check
python3 -m unittest discover -s elliptic-curves/rank-jump -p test_shared_value_block.py
```

## Ranked implications for the main rank-jump question

1. **Constructively verified mechanism:** several multisections sharing
   one nontrivial value squareclass can become rational together over a
   low-degree base. Here geometry and Galois action prove a two-direction
   rational block, including in an S₃ family.
2. **Necessary companion test:** common solubility does not ensure
   specialization independence. The explicit collapse locus gives a
   counterexample and a model for the relations that a candidate mechanism
   must exclude.
3. **Weak explanations:** simply sharing a cubic field, having many
   local Selmer classes, or finding many point representations remains
   insufficient. This mechanism succeeds because it supplies both actual
   sections and an independent rank calculation.
4. **Missing production computation:** identify an equation-defined common
   squareclass or low-degree auxiliary base carrying several genuinely new
   multisections of an existing production family. Its generic rank and
   independence must be proved before any claim about high-jump incidence.
   No such production construction is certified by this control.
5. **Potential information for Agent 1:** a proven common-cover condition
   could become a simultaneous-solubility selector. A generic-coordinate
   visibility score could not substitute for it. No current search policy,
   candidate set or worker configuration was changed.

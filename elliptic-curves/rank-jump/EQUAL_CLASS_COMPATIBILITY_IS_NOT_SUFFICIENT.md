# Equal Kummer classes need not create a geometric direction

The compatible three-point construction has **geometric mixed rank zero on
the retained low-gain anchor**. On the matched high-gain anchor the mixed
rank is still `0–1`. Thus the necessary additive compatibility condition
does not imply even geometric incidence, let alone a new rational direction.

| MW16-04 anchor | Original generic rank / witness quotient | Equal-class triple: geometric mixed rank | Arithmetic rank over the full new base |
|---|---|---:|---:|
| `-1647/91` | `16 / 9` | `0–1`, unresolved | `3–4` |
| `-2177/2397` | `16 / observed 0`, censored | **0** | **3** |

These are function-field ranks of a new pencil through each anchor.
Neither row gives a rank upper bound for its production curve. The low
row's original observed quotient remains censored.

## Construction and necessary gate

Start from the three independent generic points `P0,P1,P2` in the pinned
[paired inputs](../../artifacts/generated-results/elliptic-curves/rank_jump_triple_character_inputs_v1.json).
Use
\[
R_0=P_0,\qquad R_1=P_0+2P_1,\qquad R_2=P_0+2P_2.
\]
Their change-of-basis determinant is four, so they remain independent in
the free group but have the same Kummer class. Exact point coordinates,
chord identities, and square roots of their class ratios are retained in
the [component certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_component_kummer_gate_v1.json).
No exceptional point enters this construction.

Write `Ri=(ai,yi)`, `f(T)=T^3+AT+B`, and use
\[
E_u:\ y^2=x^3+2Au x^2+(A+3Bu+A^2u^2)x+B+ABu^2-B^2u^3.
\]
The full base is `C: zi^2=1-ai*u`; the triple character is the twist
\(d=(1-a_0u)(1-a_1u)(1-a_2u)\). Its three finite additive places pass
the [necessary Kummer compatibility theorem](ADDITIVE_COMPONENT_KUMMER_COMPATIBILITY.md):
all required classes are `[a0-theta]`. Infinity is smooth, so there is no
additional additive condition there. This is not a proof of all-place
local solubility for a particular two-cover.

The three pair characters still have arithmetic rank zero: their common
finite class differs from `[-B*theta]`, the required infinity class.
The three singleton characters each have arithmetic rank one; the trivial
character has arithmetic rank zero because `B` is nonsquare. Hence
\[
\operatorname{rank}E_u(\mathbb Q(C))
=3+\operatorname{rank}E_u^d(\mathbb Q(u)),
\qquad 0\le\operatorname{rank}E_u^d(\mathbb Q(u))\le1.
\]
The earlier non-torsion proofs for the *distinct-class* auxiliary bases
do not prove positive rank for these changed parameter curves. Their
rational rank remains unknown.

## Frozen two-place Picard test

The [protocol](EQUAL_CLASS_PICARD_PROTOCOL.json) selects the first two
eligible retained primes, adding the first good prime above 19 and at
most 31 only when necessary. Eligibility depends on the equation and
branch positions; it was checked before inspecting the weighted traces.
This gives `13,17` on the high anchor and `17,31` on the low anchor.
Only the last prime requires new untwisted counts.

The triple twist has three `I2` fibres and three `I0*` fibres, with smooth
infinity. Its trivial lattice has rank 17. The geometry, trace formula,
and two-place discriminant argument are established in the
[three-character proof](TRIPLE_CHARACTER_GENUS_ONE_BOUNDS.md).
For each good reduction let `T1,T2,T3` be the complementary five-dimensional
Frobenius traces over degrees one, two and three. The computation gives:

| Anchor | Prime | `(T1,T2,T3)` | Geometric reduction Picard rank | NS discriminant squareclass |
|---|---:|---|---:|---:|
| high | 13 | `(11,719,-463)` | 18 | `-1` |
| high | 17 | `(27,559,17133)` | 18 | `-1` |
| low | 17 | `(-1,239,-6427)` | 18 | `-1` |
| low | 31 | `(3,2055,-85905)` | 18 | `-113` |

In every row the normalized complementary polynomial has exactly one
cyclotomic factor, `X-1`. The noncyclotomic quartic supplies the Artin–Tate
determinant. Equal-rank specialization would preserve the Néron–Severi
discriminant squareclass. The low anchor's mismatch therefore forces
characteristic-zero Picard rank 17, and Shioda–Tate gives geometric mixed
rank zero. Arithmetic mixed rank is then zero too.

The high anchor's matching squareclasses do **not** prove Picard rank 18.
They merely leave the existing interval unresolved. There is no certified
high/low difference in positive mixed rank.

## Independent verification and limits

The [input projection](../../artifacts/generated-results/elliptic-curves/rank_jump_equal_class_inputs_v1.json),
[raw arrays](../../artifacts/generated-results/elliptic-curves/rank_jump_equal_class_counts_v1.json),
[Picard report](../../artifacts/generated-results/elliptic-curves/rank_jump_equal_class_picard_v1.json),
and [independent certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_equal_class_verification_v1.json)
bind the source versions and earlier proofs.

Sage/PARI supplied the new untwisted fibre counts. Independent NumPy field
arithmetic checked all 30,783 new base parameters at prime 31, by direct
character summation on 10,478 independently computed Frobenius orbits and
comparison with every stored orbit entry. This includes singular fibres.
The 12,817 reused parameters are bound to the prior independent replay.
Every new twist weight and every smooth infinity contribution was recomputed.
Companion matrices independently check the three traces, cyclotomic
dimensions, and degree-six Artin–Tate determinants.

The new count had a 40-second cap; each independent prime replay had a
60-second cap. The complete first independent verification took 36.415
seconds. A fresh replay uses temporary checkpoints and does not overwrite
the retained evidence:

```sh
sage -python elliptic-curves/rank-jump/equal_class_picard.py check
sage -python elliptic-curves/rank-jump/replay_equal_class_picard.py
python3 -m unittest discover -s elliptic-curves/rank-jump -p test_equal_class_picard.py
```

No parameter sweep, point search, class-group computation, production
selector, or active search output was changed.

## Mechanisms, missing implications, and handoff

1. **Incidence:** a geometric Picard increase is a precise necessary event
   for a mixed section in this construction. The low control proves that
   it is not forced by coincident Kummer classes. The high control still
   needs an exact geometric rank determination.
2. **Solubility:** additive-component compatibility is an exact exclusion
   theorem for rational sections. Passing it neither creates the geometric
   class nor proves rational descent. A section, its field of definition,
   and its specialization modulo the original generic subgroup are all
   still missing on the high control.
3. **Weakened explanations:** common Kummer classes, a genus-one parameter
   curve with a rational origin, and several branch characters do not
   suffice. Coincident classes can be manufactured from generic points on
   either anchor. Matching reduction discriminants are also insufficient.
4. **Visibility:** nothing here measures search visibility. The geometric
   zero is not a failure to find a point and cannot veto the low production
   fibre's further point search.

For the larger rank-jump question, the strongest unresolved route remains
an equation-defined Selmer incidence block together with a simultaneous
global-solubility certificate. The retained
[strict Selmer class-field blocks](STRICT_SELMER_AND_ARTIN_BLOCKS.md) and
[cup-product obstruction](CUP_IDEAL_AND_STRICT_LIFTING_OBSTRUCTION.md)
describe those two layers separately. Their missing implication is a
criterion forcing several classes to be rational rather than Sha, followed
by exact quotient independence. This equal-class test does not close it.

Agent 1 can eventually use a proved incidence-and-solubility criterion
available from equations and generic data. For now these are construction
constraints, not a validated score. The next small falsification test is
one additional coefficient-selected good reduction of the high equal-class
surface: a discriminant mismatch would close its geometric rank at zero;
another match would retain `UNKNOWN`. Any later section certificate must
also address its image modulo the original rank-16 subgroup; a generic
Kummer image alone does not settle free quotient membership.

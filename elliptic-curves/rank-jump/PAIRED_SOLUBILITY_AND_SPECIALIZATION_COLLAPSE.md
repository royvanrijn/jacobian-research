# Four soluble covers give exactly three directions on both high fibres

Both successful four-cover systems give **exactly three** directions modulo
the marked generic rank-17 subgroup. All eight constructed points belong
integrally to the retained independent witness groups. On the +8 fibre,
the two covers `0911e` and `1795d` give opposite representatives of the
**same quotient line**. Their simultaneous rational solubility therefore
supplies one displayed direction, despite their two independent generic
characters. On the +7 fibre, every pair is independent but a triple is
dependent modulo the generic subgroup.

This closes the earlier interval 3..4 without changing any full fibre rank
claim. The known directions still unexplained by these quartets number
four and five. A separate, fixed calculation bounds all 46 character twists
of the previously frozen cover systems. In particular, the `0911e` twist
has generic rank at most six: at least two directions in the +8 witness
group lie outside the specialized full generic group of that single cover,
even allowing undisplayed sections of arbitrarily large height.

## Exact quotient accounting

Let G1,...,G17 be the generic prefix of the retained independent witness
basis, and let Q1,...,Qj be its remaining generators modulo G. Labels below
omit `orbit-`. The points use the positive square-root branches fixed in
the original solubility certificate.

| +7 cover | Its class modulo G | +8 cover | Its class modulo G |
|---|---|---|---|
| 01333 | -Q5 | 0911e | Q4 |
| 0b2d0 | Q3 | 0a037 | -Q2 |
| 13109 | -Q1+Q3 | 1795d | -Q4 |
| 19e45 | Q3+Q5 | 18f5d | -Q5 |

The sole rational relation in each quartet, up to scaling, is

\[
P_{01333}-P_{0b2d0}+P_{19e45}\in G,\qquad
P_{0911e}+P_{1795d}\in G.
\]

The certificate retains the exact integral words in G, as well as all
eight complete integral witness coordinates. Conjugating a cover replaces
its point by a generic translate of its negative, so signs depend on the
branch choice but the quotient lines and their dependencies do not.

| Fibre | Published R17 parameter | Witness rank | Generic rank | Witness quotient rank | Quartet quotient rank | Witness quotient beyond quartet |
|---|---:|---:|---:|---:|---:|---:|
| 08234-003 | -288/65 | 24 | 17 | 7 | **3** | **4** |
| 08234-009 | -4112/1937 | 25 | 17 | 8 | **3** | **5** |

These are ranks of specified subgroups and their rational spans. Full
ranks of the two original fibres remain **UNKNOWN**.

The [relation protocol](PAIRED_QUARTET_RELATIONS_PROTOCOL.json) permits only
these two frozen fibres. Canonical-height matrices at 100 decimal digits
propose witness coordinates with denominator at most 64; both workers are
bounded at 60 seconds. Every proposal is then checked by exact rational
elliptic group addition. All denominators are one. The independent verifier
uses a separate pure-Python rational group law and finite Kummer signatures
to reprove witness independence and every relation; it uses no numerical
heights. The retained exceptional basis is an **oracle retrospective
diagnostic**, never a prospective selector.

## Generic capacity versus global solubility

The [character protocol](PAIRED_CHARACTER_MOMENTS_PROTOCOL.json) freezes
the 15 nonempty products in each of the two positive quartets and the
previously obstructed quartet, plus the cross-group pair F,D: 46 twists.
There are no new rational parameters or new fibre counts. All traces are
reweighted from the already verified complete ledger over F131 and F131².

Every product of k quadratics remains squarefree of degree 2k, coprime to
the original degree-24 discriminant at 131. Thus all 46 pass the same
good-reduction argument as the
[single-twist bound](NATIVE_SINGLE_COVER_CANNOT_EXPLAIN_THE_WHOLE_PLUS8.md).
The nontrivial cohomology dimension is N=20+4k. Exact first and second
Frobenius moments give an upper bound on the central eigenvalue
multiplicity; the functional-equation sign refines that upper bound using
the parity of the **analytic multiplicity**, without invoking a
Mordell–Weil parity conjecture. The independent verifier recomputes the
characters in actual finite fields and checks local root numbers in the
residue fields of the discriminant factors.

| Frozen system | Singleton twist upper bounds, in stored label order | Pair bounds | Triple bounds | Quartet bound |
|---|---|---|---|---:|
| +7 positive quartet | 7, 8, 8, 8 | 6..9 | 8..9 | 12 |
| +8 positive quartet | 6, 8, 7, 8 | 8..10 | 9..10 | 10 |
| Obstructed quartet A,B,C,D | 6, 7, 8, 7 | 7..9 | 8..10 | 9 |

Each singleton has a known non-torsion generic section, hence lower bound
one. Higher products have only lower bound zero here. The displayed upper
bounds are neither rank estimates nor evidence that their capacity is used.

Three fixed pair carriers now have the following comparison. The middle
rank column concerns the **original elliptic curve over the carrier's
function field**; the Jacobian column concerns a different elliptic curve.

| Carrier | Original E generic rank over carrier | Carrier Jacobian rank | Carrier rational solubility | Retained specialization evidence |
|---|---:|---:|---|---|
| F,G = 1795d,0911e | 19..40 | 2 | Yes | Both split on +8, but displayed quotient rank **1** |
| F,D = 1795d,11278 | 19..38 | 3 | Yes | No simultaneous split in the 32-fibre sample |
| A,D = 030cb,11278 | 19..38 | 2 | **No: nonzero Sha class**, despite every local point | No rational simultaneous specialization anywhere |

The global claims are inherited from the independently verified
[positive carriers](GLOBAL_CARRIER_SOLUBILITY_AND_SPECIALIZATION.md) and
[labelled A,D Cassels–Tate obstruction](NATIVE_PAIR_CARRIER_HAS_A_SHA_OBSTRUCTION.md).
All three generic lower bounds are 17+1+1 by distinct quadratic characters.
The equal intervals for F,D and A,D do **not** prove equal actual generic
ranks. They show only that these certified capacity bounds do not separate
global solubility. Even a globally empty carrier has a function field and
these generic sections.

For F,G, the existing generic independence theorem and the new exact
specialized relation demonstrate that generic independence need not survive
at the selected rational lift. This failure occurs after global solubility;
it is not a point-visibility failure.

For the single `0911e` cover, character decomposition gives full generic
rank at most 17+6=23. Its rational lift on the rank-25 witness fibre implies
that at least 25-23=2 witness directions lie outside the rational span of
the specialized full generic cover group. This is stronger than the
previous bound of one from `1795d`. It does not identify the outside points
or compute the exact twist rank.

## What remains of the proposed mechanism

The current partial chain is

\[
\text{simultaneous native square conditions}
\Longrightarrow\text{four rational points}
\Longrightarrow\text{exactly three quotient directions on each high fibre}.
\]

The first implication has explicit generic maps, and the second now has
exact specialization certificates. The missing positive implication is a
simpler pre-point arithmetic condition forcing a sufficiently large
**independent** soluble block. Both full quartet carriers have
[genus 17 and gonality 8](SOLUBLE_QUARTETS_REQUIRE_HIGHER_GENUS_LIFTS.md),
so a rational or elliptic parametrization of the entire quartet is excluded.
An auxiliary elliptic point on F,G does not by itself lift through the
remaining conditions, and here even its two displayed cover points collapse
to one quotient line.

For example, an identity Pj=S±Pi with S generic forces another rational
point but adds no new quotient line. Such group-law constructions must be
accounted for before counting simultaneous square conditions as dimensions.
The relations found here were discovered retrospectively; no equation for
their specialization locus, independent of the supplied points, is claimed.

The ranked implications for further work are:

1. **Solubility:** a labelled global obstruction or explicit simultaneous
   lift remains the strongest demonstrated discriminator. The native A,D
   example proves that local solubility and substantial generic capacity
   can coexist with a global obstruction. Higher-genus lifting is still
   unresolved for general positive specializations.
2. **Incidence:** independent generic characters must be followed by an
   independence check after specialization. The +8 pair gives a concrete
   counterexample to counting two soluble native covers as two new
   directions. It would be useful to derive its collapse locus directly
   from generic lift maps, before supplying exceptional points.
3. **Weak explanations:** large generic upper bounds, distinct cover
   labels, and pair-carrier Jacobian rank do not yet predict a large jump.
   A single `0911e` pullback, even with all hidden generic sections, cannot
   account for the whole observed +8 witness group.
4. **Missing computations:** identify the remaining four and five witness
   directions by pre-point constructions; separate rationally soluble
   higher-genus lifts from Sha obstructions; and prove independence for
   the resulting block. These are distinct obligations.
5. **Visibility:** no new visibility feature is proposed. Agent 1 can
   eventually use a certified simultaneous-lift condition with a proved
   quotient-dimension guarantee. Raw counts of split charts have no such
   guarantee, and these retrospective quartets must not be used as an
   oracle-trained prospective selector.

## Replay and immutable evidence

From the repository root (Sage 10.9 for the finite-field checks):

```sh
sage -python elliptic-curves/rank-jump/paired_character_moments.py check
sage -python elliptic-curves/rank-jump/verify_paired_character_moments.py check
python3 elliptic-curves/rank-jump/verify_paired_quartet_relations.py check
```

The frozen [relation input](../../artifacts/generated-results/elliptic-curves/rank_jump_paired_quartet_relations_inputs_v1.json),
[exact coordinates](../../artifacts/generated-results/elliptic-curves/rank_jump_paired_quartet_relations_v1.json),
[independent group-law replay](../../artifacts/generated-results/elliptic-curves/rank_jump_paired_quartet_relations_verification_v1.json),
[46 character bounds](../../artifacts/generated-results/elliptic-curves/rank_jump_paired_character_moments_v1.json),
and [independent arithmetic replay](../../artifacts/generated-results/elliptic-curves/rank_jump_paired_character_moments_verification_v1.json)
bind their scripts and mathematical inputs by hash. No live search output,
protocol, selector, worker setting, or mathematical status entry was changed.

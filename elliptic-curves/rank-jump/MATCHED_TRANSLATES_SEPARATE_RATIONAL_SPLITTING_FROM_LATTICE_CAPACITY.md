# Matched translates separate rational splitting from lattice capacity

Three translates chosen by the generic lattice alone give the same
intersection degree twelve on the same successful native triple. Two have
irreducible degree-twelve intersection polynomials and no rational points.
The third has the known rational component and an irreducible degree-eleven
remainder. Their Galois actions are respectively **S12, 1+S11, S12**.

The negative schemes already have local obstructions at 61 and 53. Thus the
generic trace norms and intersection capacity do not force simultaneous
rational solubility. In this comparison, the distinguishing event is a
rational fixed point of the Galois action on a finite intersection scheme.
The original genus-five carrier is fixed and globally soluble throughout;
the extra group-relation condition is what fails on the controls.

## Fixed carrier and a selector without the successful word

The [preceding minimal-carrier result](MINIMAL_CARRIER_AND_RATIONAL_SPLITTING_OF_A_TWO_DIRECTION_BLOCK.md)
fixes A=`01333`, B=`0b2d0`, D=`19e45` and their native traces
\(\tau_A,\tau_B,\tau_D\). It studies

\[
P_A-P_B+P_D=S,
\]

as the intersection of A's rational bisection with the curve traced by
\(S+P_B-P_D\). Set

\[
w=\tau_A-\tau_B+\tau_D,\qquad R=2S-w.
\]

The previous proof gives intersection degree \(h(R)+2\). The
[new protocol](TRIPLE_TRANSLATE_CONTROL_PROTOCOL.json) selects translates
using only the generic height Gram matrix and these three trace vectors.
The selection function receives no successful word, exceptional point,
parameter, or rank outcome.

It enumerates the vectors of norm at most ten in \(2M+\mathbf Z w\).
Since M has minimum four, every nonzero vector this short lies in the
coset \(w+2M\). There are exactly **34 signed vectors**, all of norm ten.
Retaining \(S=(R+w)/2\) of height ten leaves 14 words. Simultaneous cover
conjugation identifies \(S\) with \(w-S\), leaving seven classes.

Each class is represented by the least admissible word under coefficient
l1 norm, then maximum absolute coefficient, then lexicographic order.
The first three classes are frozen before the known successful word is
read for comparison. The second selected class happens to be the conjugate
of that word. It is reported as a calibration, not as independent positive
evidence or a new specialization.

The trace triple itself was selected retrospectively in the preceding
experiment. This selector isolates the effect of the translate within that
fixed triple; it does not establish an oracle-free choice of the triple.

## Exact matched results

Use P1,...,P17 for the published generic basis.

| Control | Generic translate S | h(S) | h(2S-w) | Image D.O | Intersection factors over Q | Rational intersection |
|---|---|---:|---:|---:|---|---|
| 0 | P4-P8-P12+P14 | 10 | 10 | 8 | 12 | None |
| 1 | P7-P8+P9-P12+P14-P16+P17 | 10 | 10 | 10 | 1+11 | t=-288/65 |
| 2 | -P2+P4-P5+P7-P10+P11+P13+P14-P16+P17 | 10 | 10 | 10 | 12 | None |

All three image curves have degree four over the original parameter line,
self-intersection 16, arithmetic genus nine, and normalization genus one.
All intersect A's bisection in twelve points. In particular controls 1
and 2 also match in their intersection with the zero section. This is a
comparison of relation divisors on one fixed carrier, not a measurement
of low ranks on newly chosen curves.

The existing exact triple compiler was imported and used without changing
it. Each worker had a 60-second limit. All completed. Their affine-open
intersection polynomials are squarefree of degree twelve, and exact
rational functions provide the three cover roots and points over each
finite parameter algebra. A separate explicit group law verifies the
relation on every geometric point. The twelve distinct constructed points
exhaust the proper intersection by its independently calculated degree.
There are no unexamined extra points at excluded denominators or infinity.

Control 1 has exactly the same intersection polynomial as the previous
successful word. Its three root functions are negatives of the previous
ones in that finite algebra. Hence it is precisely the simultaneous-sign
conjugate construction. Among the 32 frozen parameters, its only hit is
the retained `08234-003`; controls 0 and 2 have none. No new parameter was
evaluated and no original-curve point search was run.

## Local obstructions before rational factor extraction

At the following primes, the primitive integral intersection polynomial
has unit leading coefficient and squarefree reduction:

| Control | Prime | Factor degrees modulo p | Consequence |
|---|---:|---|---|
| 0 | 61 | 12 | No point over Q61 |
| 2 | 53 | 2+3+7 | No point over Q53 |

Neither reduction has an Fp root. Any Qp root of this integral polynomial
would be integral because its leading coefficient is a p-adic unit, and
would reduce to an Fp root. Thus these are exact local obstructions, not
bounded failures to find rational points.

The finite intersection algebra already includes every characteristic-zero
intersection point. These obstructions therefore apply to the entire
specified relation scheme. They are **not** obstructions to the fixed
genus-five carrier, which has the known rational point for control 1, and
they are not Sha assertions about the original elliptic fibre.

This supplies an interpretable pre-point **solubility** test for a proposed
relation construction: a good-prime reduction without a degree-one factor
excludes rational solubility of that construction. Passing several such
tests does not prove global solubility or rank gain. The primes above are
the first certified obstructions in the fixed list of primes at most 131,
not a claim about the smallest possible obstructing primes under every
integral presentation.

## Exact Galois certificates

These are Galois groups of the finite intersection schemes, not of the
original elliptic curve's torsion field. All modular factorizations used
below preserve degree and are squarefree, so they give Frobenius cycle
types in the characteristic-zero permutation action.

For **control 0**, irreducibility modulo 61 proves irreducibility over Q
and hence transitivity. Modulo 79 the type 1+11 supplies an 11-cycle,
which excludes every nontrivial block system in degree twelve. Modulo 73
the type 1+2+9 supplies a transposition by taking the ninth power. A
primitive group containing a transposition is S12.

For **control 2**, no single irreducible reduction occurs in the bounded
prime list, but three reductions prove rational irreducibility. Any proper
rational factor degree must be a subset sum of every modular factor-degree
list:

| Prime | Factor degrees | Remaining possible proper Q-factor degrees |
|---:|---|---|
| 53 | 2+3+7 | 2,3,5,7,9,10 |
| 73 | 1+5+6 | 5,7 |
| 89 | 1+2+9 | None |

Thus the action is transitive. The type 2+3+7 at 53 supplies a 7-cycle by
taking the sixth power, and a transposition by taking the twenty-first
power. A prime cycle longer than half the degree excludes imprimitivity:
in a nontrivial block system both the number and size of the blocks are
at most six, and the corresponding wreath product has no element of
order seven. The group is therefore S12.

For the **degree-eleven factor of control 1**, irreducibility modulo 73
gives a transitive action of prime degree, hence a primitive action.
Modulo 79 that factor has type 1+2+3+5; its fifteenth power is a
transposition. Its group is S11. The rational root is fixed, giving the
full scheme action 1+S11.

For completeness, the primitive-plus-transposition step is elementary:
the edges given by all conjugates of the transposition form a
group-invariant graph. Its connected components are blocks. Primitivity
makes it connected, and transpositions along a connected graph generate
the whole symmetric group. No statistical Galois-group estimate is used.

## Implications for the rank-jump mechanism

1. **Solubility:** the successful construction has a rational Galois-fixed
   component in a degree-twelve intersection. The matched failures have
   full S12 action and explicit local obstructions. This isolates an
   arithmetic splitting event beyond the generic trace lattice.
2. **Weak explanation, now falsified:** height-ten translates, a minimal
   norm-ten residual trace, intersection degree twelve, and the same
   globally soluble carrier do not force a rational relation point.
   Controls 1 and 2 retain even the same zero-section intersection.
3. **Incidence:** on the successful component the retained independent
   certificate gives two quotient directions. Solubility of the relation
   scheme and independence remain separate obligations. The experiment
   does not explain the full +7 quotient or the four remaining witness
   directions beyond its rank-three quartet.
4. **Missing implication:** find a condition on a triple selected without
   exceptional data that forces a rational intersection component and
   preserves multiple independent directions. The lattice-only translate
   selector is now separated from the oracle word, but the input triple
   still comes from a successful retrospective block.
5. **Information for Agent 1:** local degree-one-factor tests can exclude
   a proposed simultaneous relation construction. Exact rational-factor
   extraction plus lift maps can certify its global solubility. These are
   solubility features of a construction, not visibility scores or rank
   predictors for arbitrary parameters. No search settings were changed.

The next useful question is whether this local obstruction can be obtained
directly from the native cover and relation equations, before expensive
characteristic-zero elimination, while accounting for every chart and bad
fibre. A successful finite-field implementation would give a cheap necessary
condition for these arithmetic blocks; passing it would still leave the
global rational-component condition open.

## Reproduction

Sage 10.9 and PARI 2.17.3 were used. The independent coset replay performs
exact rational LDL enumeration, visits 1,878 nodes, and recovers the same
34 signed vectors and seven eligible conjugacy classes. It uses no
numerical-height assumption or PARI short-vector enumeration.

```sh
sage -python elliptic-curves/rank-jump/triple_translate_controls.py check --case 0
sage -python elliptic-curves/rank-jump/triple_translate_controls.py check --case 1
sage -python elliptic-curves/rank-jump/triple_translate_controls.py check --case 2
sage -python elliptic-curves/rank-jump/verify_triple_translate_controls.py check
```

The immutable [selection](../../artifacts/generated-results/elliptic-curves/rank_jump_triple_translate_selection_v1.json),
[completed control manifest](../../artifacts/generated-results/elliptic-curves/rank_jump_triple_translate_controls_v1.json),
and [independent verification](../../artifacts/generated-results/elliptic-curves/rank_jump_triple_translate_controls_verification_v1.json)
bind all three case inputs, complete intersection polynomials, lift maps,
source scripts, and protocols. Active search files and mathematical status
entries were not modified.

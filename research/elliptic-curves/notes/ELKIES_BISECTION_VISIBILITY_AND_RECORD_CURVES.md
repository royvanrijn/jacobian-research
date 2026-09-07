# Bisection visibility and exact published-R17 fibre exclusions

<!-- status-consumer: EC-K3-ELKIES-2026-BISECTION-VISIBILITY-RECORD-CURVES 1c39220ee5fedc77 -->

## Outcome

The complete specialization census changes the next search decision in two
ways.

1. There is no new bisection trace shell to construct.  The existing 39,120
   records already represent every section-nonnegative degree-two `(-2)`-curve
   modulo generic-section translation and sign.  The underlying complete
   norm-ten shell has 806,238 unoriented representatives.  Replacing a retained
   trace by a higher-height translate gives the same bisection orbit, the same
   quadratic extension and the same class modulo the generic seventeen.
2. The 2024 Elkies--Klagsbrun rank-29 curve and ICARM curves 273, 302, and
   398--400 cannot be inserted as rational fibres of the published `R17`
   fibration.  Their exact degree-24 recognition equations against its `j`-map
   are irreducible over `QQ`.  This remains true after quadratic twisting,
   because twisting preserves `j`.

The right continuation is therefore a new geometric or arithmetic mechanism
targeted at the missing rank-28 quotient classes, while the existing
bisections remain a strong conductor-search filter.

## Exact visibility spaces

Row reduction of every split-bisection class in the stored ordered public
complement gives:

| parameter | known complement | visible span | invisible quotient |
| --- | ---: | ---: | ---: |
| `-2/377` | 8 | 5 | 3 |
| `-308/251` | 9 | 3 | 6 |
| `2456/135` | 10 | 2 | 8 |
| `-9529/5471` | 11 | 1 | 10 |
| `3/8` | 4 | 4 | 0 |

For the rank-28 fibre, the sole visible class is

```text
01011001010 = Q2 + Q4 + Q5 + Q8 + Q10.
```

Its first pivot is `Q2`.  The deterministic row-reduction complement is

```text
Q1, Q3, Q4, Q5, Q6, Q7, Q8, Q9, Q10, Q11.
```

In the original public point list these are points

```text
1, 3, 4, 7, 8, 9, 11, 15, 19, 22.
```

This complement is a reproducible coordinate choice, not an intrinsic direct
sum decomposition of the Mordell--Weil group.  It is nevertheless the exact
ten-class packet against which a different point-construction mechanism can
be scored.

At `t=3/8` the invisible quotient is zero.  The current atlas should therefore
be kept as an inexpensive pre-descent gate for low-conductor searches: a
candidate on which the split classes reproduce the expected exceptional span
deserves the more expensive exact descent and saturation work.

## Why another trace shell cannot help

The complete bisection theorem is a quotient statement, not a bounded
height-ten point search.  The 806,238 norm-ten representatives map onto all
39,120 surviving translation classes, and the equation atlas contains one
certified record for every class.  Section translation acts within one such
class.  At a rational specialization it translates a split point by a member
of the generic `R17` subgroup, so its finite-quotient class modulo `R17` is
unchanged.

Thus enlarging the trace height while retaining degree-two rational
bisections can only repeat the existing visibility space.  To reach one of
the ten missing rank-28 classes, at least one ingredient must change.  The
useful choices are:

- higher-degree multisections whose splitting fields are not in the complete
  quadratic bisection atlas;
- explicit residual 2-Selmer covers, after the repository's Selmer gate;
- direct point searches whose candidates are rejected immediately unless
  their exact finite-quotient class meets the ten-class target packet.

The third route is the cheapest immediate experiment.  The first two are the
structural routes if direct searches do not hit the target packet.

## Exact fibre-recognition test

For the published short family

```text
y^2 = x^3 + A(t)x + B(t),
c4_R17(t) = -48 A(t),
Delta_R17(t) = -16(4 A(t)^3 + 27 B(t)^2),
```

the `j`-map is reduced and has numerator/denominator degrees `24/24`.  For a
target curve with invariants `(c4,Delta)`, a necessary specialization equation
is

```text
c4_R17(t)^3 Delta - c4^3 Delta_R17(t) = 0.
```

After removing integer content, all six target equations are primitive
degree-24 polynomials:

| target | certified rank lower bound | finite-field witness | conclusion |
| --- | ---: | --- | --- |
| 2024 Elkies--Klagsbrun / ICARM 12 | 29 | irreducible of degree 24 modulo `461` | no rational `t` |
| ICARM curve 273 | 30 | irreducible of degree 24 modulo `367` | no rational `t` |
| ICARM curve 302 | 31 | irreducible of degree 24 modulo `397` | no rational `t` |
| ICARM curve 398 | 30 | irreducible of degree 24 modulo `1009` | no rational `t` |
| ICARM curve 399 | 29 | irreducible of degree 24 modulo `83` | no rational `t` |
| ICARM curve 400 | 28 | irreducible of degree 24 modulo `157` | no rational `t` |

Each reduction retains degree 24.  Gauss's lemma therefore proves that each
integer polynomial is irreducible over `QQ`; in particular none has a rational
affine root.  Their leading coefficients are nonzero, so the point at infinity
is not a solution either.

The historical rank-28 fibre is a calibration control.  At its published
parameter the same construction factors exactly as

```text
(5471*t + 9529) * (irreducible degree-23 cofactor).
```

The cofactor is irreducible modulo `197`, and the linear factor recovers
`t=-9529/5471`.  Thus the exclusion test detects a known rational fibre rather
than rejecting it through a normalization error.

This rules out rational fibres in the published fibration, not merely failure
of one affine chart, and the conclusion is twist-stable.  It does not rule out
an isogenous family, another elliptic fibration on the same K3, or an unrelated
construction.  The rank-29 provenance makes the second possibility the main
structural target; the exact marking of that other fibration has not been
published.

The public generalized Weierstrass equations and current metadata for ICARM
398--400 are pinned in
[`elkies_2026_r17_j_recognition_targets.json`](../data/elkies_2026_r17_j_recognition_targets.json).
This snapshot is input to `j`-recognition only; it does not duplicate the
ICARM point-independence certificates.

## Replay and claim boundary

```bash
.venv/bin/python \
  elliptic-curves/scripts/analyze_elkies_bisection_visibility_and_record_curves.py \
  --check
```

The pinned certificate is
[`elkies_2026_bisection_visibility_record_curves_v1.json`](../../artifacts/generated-results/elliptic-curves/elkies_2026_bisection_visibility_record_curves_v1.json).
Its SHA-256 is

```text
4e7b9fcbae2eae2d950fba867b9a6a27fe7352e9af15b09eb5361d9486f8c6ca
```

The visibility calculation concerns finite-quotient classes inside already
known public complements.  It supplies no upper bound for any Mordell--Weil
rank.  The provenance calculation is an exact non-membership result for the
published `R17` fibration only.

# Bisection visibility, the rank-28 target, and the rank-30/31 records

<!-- status-consumer: EC-K3-ELKIES-2026-BISECTION-VISIBILITY-RECORD-CURVES c3c7fbac71f4624f -->

## Outcome

The complete specialization census changes the next search decision in two
ways.

1. There is no new bisection trace shell to construct.  The existing 39,120
   records already represent every section-nonnegative degree-two `(-2)`-curve
   modulo generic-section translation and sign.  The underlying complete
   norm-ten shell has 806,238 unoriented representatives.  Replacing a retained
   trace by a higher-height translate gives the same bisection orbit, the same
   quadratic extension and the same class modulo the generic seventeen.
2. ICARM curves 273 and 302 cannot be inserted as additional fibres of this
   atlas.  Their exact degree-24 recognition equations against the published
   `R17` `j`-map are irreducible over `QQ`, so neither curve is a direct rational
   specialization of the published family.

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

## Exact record-curve provenance test

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

After removing integer content, both target equations are primitive degree-24
polynomials:

| target | certified rank lower bound | finite-field witness | conclusion |
| --- | ---: | --- | --- |
| ICARM curve 273 | 30 | irreducible of degree 24 modulo `367` | no rational `t` |
| ICARM curve 302 | 31 | irreducible of degree 24 modulo `397` | no rational `t` |

The reductions retain degree 24.  Gauss's lemma therefore proves that the
integer polynomials are irreducible over `QQ`; in particular they have no
rational affine roots.  Their leading coefficients are nonzero, so the point
at infinity is not a solution either.

This rules out a direct rational specialization in the published chart.  It
does not rule out an isogenous family, a different K3 fibration, or an
unrelated construction.  Curves 273 and 302 should consequently remain
external high-rank controls.  Their exact finite-quotient certificates are
still useful implementation regressions, but their 30 and 31 points cannot be
labelled `R17 + 13` and `R17 + 14` without a different exact family
identification.

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
112cbe9030fc8657a0622ae34e2681328f91c8cfcb390d4f1d5acdbf88099bc9
```

The visibility calculation concerns finite-quotient classes inside already
known public complements.  It supplies no upper bound for any Mordell--Weil
rank.  The provenance calculation is an exact non-membership result for the
published `R17` family only.

# ICARM curve 356: rank at least 29, size record, and construction fingerprint

Status: **exact unconditional rank lower bound and exact local arithmetic;
construction provenance unresolved**.

The curve-351 comparison below is a bounded numerical fingerprint.  It is not
a family-recognition theorem and does not identify either curve with the
Elkies H3/rootless-`MW17` construction.

## Bottom line

[ICARM curve 356](https://elliptic-rank.icarm.cloud/curve/356), submitted by
`wgxli` on 2026-08-25 at 14:20:06 UTC, is

```text
y^2 = x^3 + x^2 + A*x + B

A = -24391876744717707263532695900840552395172973498186560300
B = 46943906433780620456844832699051340439698711588743845207309557656274241785479710000.
```

The 29 displayed points are independently replayed in this repository, so

```text
rank E(Q) >= 29
```

unconditionally.  No unconditional upper bound is known here.  A 55-second
PARI `ellrank(E,0,points)` attempt returned no bound and was stopped; that
timeout has no mathematical meaning.

This is not a new raw-rank record: curve 302 has rank at least 31.  It is a
new public *size* record for rank at least 29.  At the retrieval time above,
the complete ICARM database gave:

| curve | rank lower bound | `log N` | naive height | Faltings height | `log |Delta|` |
|---|---:|---:|---:|---:|---:|
| curve 356 | 29 | 317.8863 | 394.2398 | 30.7273 | 383.0718 |
| 2024 Elkies--Klagsbrun, ICARM 12 | 29 | 343.7201 | 436.0125 | 34.2354 | 425.4423 |
| curve 273 | 30 | 339.3479 | 442.0854 | 34.7705 | 432.1249 |
| curve 302 | 31 | 375.2224 | 468.2771 | 36.7425 | 453.0469 |

Thus curve 356 improves all four public size coordinates for the
rank-at-least-29 threshold, not only conductor.

## Exact rank-lower-bound replay

The exact data are in
[`icarm_curve356.py`](../cas/icarm_curve356.py), and the checker is
[`verify_icarm_curve356_rank29.py`](../cas/verify_icarm_curve356_rank29.py).
It verifies the following.

1. All 29 rational pairs satisfy the public equation.
2. The rational change

   ```text
   X = 36*x + 12,
   Y = 216*y
   ```

   transports the curve and all points to the integral short model

   ```text
   Y^2 = X^3 - 27*c4*X - 54*c6.
   ```

3. The short 2-division cubic has no root modulo 17, so `E(Q)[2]=0`.
4. Exact point counts give `#E(F_7)=12` and `#E(F_17)=25`.  Their gcd is one,
   so the complete rational torsion subgroup is trivial.
5. Exhaustive finite-group calculations at

   ```text
   11,53,73,79,83,89,101,109,127,137,
   157,163,173,179,191,193,229,251,263,281
   ```

   produce 30 binary rows in the product of the quotients
   `E(F_p)/2E(F_p)`.  The 29-column matrix has rank 29.

Any integral relation among the displayed points therefore has all
coefficients even.  Dividing such a relation by two produces rational
2-torsion; because `E(Q)[2]=0`, infinite descent forces the relation to be
zero.  This exact finite-quotient implementation is separate from the
Brumer--Cremona quadratic-character verifier used by the ICARM website.

## Minimal model, bad fibres, and conductor economy

The public integral model has

```text
c4 = 1170810083746449948649569403240346514968302727912954894416

c6 = -40559535158786456074713935459005218642378386504572098678557536904830761270132198806464.
```

The reduced `j=c4^3/Delta` has SHA-256

```text
74999308db4bee588f5d8b7fe0c6a8cf5f581d5313f35c91305e1d86a790e25e.
```

Its exact discriminant factorization is

```text
Delta = -2^8 * 3^10 * 5^7 * 13^6 * 23^2 * 29^3 * 37^4
          * 41^2 * 139^2 * 751 * 28960331
          * 1204882855601765528877267647500895974865482613
          * 197980272243427555346397293722916980361535459279712115031762027678304939.
```

Every displayed discriminant valuation is below 12, which already proves
global minimality.  PARI/GP 2.15.4 independently returns the identity minimal
change at every bad prime and the fibre profile

```text
I1*@2,
I10@3, I7@5, I6@13, I2@23, I3@29, I4@37,
I2@41, I2@139,
I1@751, I1@28960331, and I1 at each of the two large primes.
```

At 2 the minimal discriminant valuation is 8 rather than the tame `I1*`
value because of wild reduction; the conductor exponent is 3.  Every odd bad
prime is multiplicative and has conductor exponent one.  Consequently

```text
N = 2^3 * 3 * 5 * 13 * 23 * 29 * 37 * 41 * 139 * 751 * 28960331
    * 1204882855601765528877267647500895974865482613
    * 197980272243427555346397293722916980361535459279712115031762027678304939.
```

This is the mechanism behind the small conductor: many powers in `Delta`
collapse to one copy of each odd prime in `N`.  The global root number is
`-1`, consistent with the proved odd rank but not an upper bound.

## What can be recovered about the construction

The public curve page contains no commentary, family equation, parameter,
search bounds, or generic-section formulas.  At the 2026-08-25 snapshot used
for the original replay, the submitter's only other public ICARM entry was
[curve 351](https://elliptic-rank.icarm.cloud/curve/351), a rank-at-least-25
curve submitted 27 minutes 44 seconds earlier.  Comparing those two records
gave a strong common-lineage signal.

The later 2026-09-01 sweep finds curves 376, 377, and 385 in the same bounded
ordered-height component and exports all five canonical short fibres as an
inverse-interpolation input.  See
[`ICARM_WGXLI_RANK17_LINEAGE.md`](ICARM_WGXLI_RANK17_LINEAGE.md).  This update
strengthens the lineage evidence but remains numerical family-recognition
evidence, not a construction certificate.

<!-- status-consumer: EC-ICARM-WGXLI-R17-LINEAGE 90790392f558f0a0 -->

### Ordered 17-point fingerprint

For a rational point, write the reduced `x` denominator as `d^2`.  At the
following positions among the first seventeen displayed points, curves 351
and 356 have exactly the same `d`:

```text
position:  2  4   5  11   13  15  16  17
d:         1  1  71   5  679   1   7  41.
```

The nontrivial ordered matches `71,5,679,7,41` are difficult to explain as
an unordered rank coincidence.  More strongly, PARI canonical-height Grams
of the two ordered 17-tuples satisfy the best scalar fit

```text
H_356 approximately 1.4208782482875446 * H_351,
relative Frobenius residual = 0.11220111822209557,
Pearson correlation of all 289 entries = 0.9748839794656168.
```

The source hashes, denominator check, and 80-digit PARI computation are
replayed by
[`analyze_icarm_curve356_lineage.py`](../cas/analyze_icarm_curve356_lineage.py).

The same command also runs the repository's exact fixed-root Mestre
recognizer.  It tests all 2,329 normalized, generically nonsingular six-root
tuples of diameter at most 300, plus the larger Fermigier control tuple.  Of
the 2,330 target `j` equations, 111 survive the modular filters to exact
factorization, but none has a rational-square parameter and there is no exact
`j` match.  This excludes precisely that bounded census.  It does not exclude
larger root tuples, generalized Mestre constructions, Nagao or Kihara
families, K3 descendants, isogenous images, or a private family.

The most economical interpretation is that the displayed ordering preserves
seventeen corresponding sections of a common non-isotrivial family.  Under
that still-unproved identification, curve 351 has at least `25-17=8` and
curve 356 at least `29-17=12` exceptional specialization directions.

### Why this is not yet a construction certificate

The inference stops at a **common ordered 17-section template**.  It does not
establish any of the following:

- that the generic Mordell--Weil rank is exactly 17;
- the family equation or the two specialization parameters;
- that the first seventeen submitted points are literally evaluations of
  rational functions over the family base;
- that the family is the Elkies--Klagsbrun K3 family;
- that it is the repository's H3/rootless-`MW17` model;
- that it lies outside every Mestre-type construction (only the declared
  2,330-family fixed-root census was excluded);
- or that the repeated discriminant valuations came from a particular CRT
  search rather than another conductor-shaping method.

The repository's calibrated height-lattice audit already shows that
unlabelled rank-17 Gram fitting can mistake a known Fermigier--Mestre control
for the R17 lattice.  The direct ordering and denominator agreement here are
stronger evidence of common lineage between 351 and 356, but they do not fix
the lineage's name.  The rank-mutation theorem likewise says what must happen
*after* a marked K3 identification; it cannot supply that identification from
an isolated rational curve.

## Exact next gates

The construction becomes reproducible only after one of the following is
available.

1. The submitter supplies the family equation, parameters for curves 351 and
   356, and the seventeen section functions.
2. An independently proposed family has an exact `j` match at rational
   parameters, followed by a `Q`-isomorphism and transport of all seventeen
   labeled points.
3. The repository's rootless H3 equation is completed and its exact `j`-map
   can be solved against the curve-356 `j` above.  Equality of `j` alone is
   insufficient because twists must still be separated.

Only then is it meaningful to study the twelve exceptional directions in the
family's Mordell--Weil quotient or to reproduce the conductor-directed search
near the curve-356 parameter.

## Reproduction

From the repository root:

```bash
PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/verify_icarm_curve356_rank29.py \
  --verify-primality

PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/analyze_icarm_curve356_lineage.py
```

The first command is exact.  The second requires network access to the two
hash-pinned ICARM JSON records and PARI/GP; its height comparison is numerical.

## Public sources

- [ICARM curve 356](https://elliptic-rank.icarm.cloud/curve/356): public
  minimal model, 29-point witness, invariants, submitter, and timestamps.
- [ICARM curve 351](https://elliptic-rank.icarm.cloud/curve/351): the earlier
  same-submitter comparison curve.
- [ICARM API documentation](https://elliptic-rank.icarm.cloud/api): exact
  independence-verification semantics and global-minimal storage convention.
- [ICARM database JSON](https://elliptic-rank.icarm.cloud/database.json): the
  size-record comparison snapshot.
- [Elkies--Klagsbrun construction account](https://mathoverflow.net/questions/477849/background-for-the-elkies-klagsbrun-curve-of-rank-29/478509): the public
  generic-`Z^17` precedent; it is context, not an identification of curve 356.

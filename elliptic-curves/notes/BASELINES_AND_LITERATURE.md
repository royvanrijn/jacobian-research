# Baselines and literature

Checked on 2026-08-20. Current-record statements are external and should be
rechecked before publication.

## Numerical baselines

Jonathan Bober's [conditional analytic-rank paper](https://arxiv.org/abs/1112.1503)
lists natural-log conductors `196.68` for the historical E21 curve and `182.72`
for Fermigier's E22 curve.  Exact PARI/GP recomputation for E22 gives

```text
N = 22720638514787473197194583889675055980109503436060704437972911338086049759883790
log(N) = 182.72491095063742879610833035152407869...
```

Thus the programme uses the user's literal strict cutoff `log(N) < 182.72`,
not the rounded E22 value.  The known E22 model and its 22 independent points
come from Fermigier's
[1997 paper](https://matwbn.icm.edu.pl/ksiazki/aa/aa82/aa8243.pdf).
Klagsbrun--Sherman--Weigandt prove exact arithmetic ranks for the historical
E20--E24 and E28 curves [subject to GRH for associated number-field zeta
functions](https://arxiv.org/abs/1606.07178).  Independence alone supplies the
unconditional lower bounds.

The current public general-rank record found in 2026 is rank at least 30. The
[maintained equation and 30 points](https://web.math.pmf.unizg.hr/~duje/tors/rk30.html)
give the public data. The repository checks those data exactly and proves the
points independent without an analytic-rank assumption. The previous 2024
[rank-at-least-29 record](https://web.math.pmf.unizg.hr/~duje/tors/rk29.html)
remains a historical calibration; its public exact-rank statement is
conditional on GRH.

## Local exact baseline replays

For the 2026 record, all 30 public points satisfy the curve equation exactly.
After transport to an integral short model, their images in products of
`E(F_p)/2E(F_p)` form a binary matrix of rank 30, and the modulo-23
2-division polynomial has no root. Infinite descent proves the points
independent. The
[pinned certificate](../../artifacts/generated-results/elliptic-curves/icarm_curve273_rank30_v1.json)
is independently replayed using Sage invariant factors and discrete
logarithms. It proves `rank >= 30` unconditionally, not `rank = 30`.

The repository now checks the public rank-at-least-29 curve without relying on
a floating-point height determinant.  At 29 cyclic good reductions, the
discrete logarithms of the 29 displayed points form a full-rank matrix modulo
2.  A separate good prime, (p=67), has group order 83 and therefore rules out
rational 2-torsion.  Infinite descent proves the points independent.  The
[pinned certificate](../../artifacts/generated-results/elliptic-curves/elkies_klagsbrun_e29_independence_v1.json)
is an unconditional rank lower bound only: it neither supplies a thirtieth
point nor replays the conditional rank upper bound.

Kihara's fully published arithmetic rank-at-least-14 family is also reproduced
as a pipeline-development fallback.  At (t=2), the exact quartic construction
gives the paper's fifteen points; taking the fifteenth as origin, a full-rank
finite-reduction matrix modulo 5 proves the other fourteen independent.  The
torsion witness is (p=11), where the group order is 18.  Independence at this
defined specialization also rules out a relation between the fourteen
rational-function sections.  The
[family record](../families/kihara_rank14.json) and
[pinned certificate](../../artifacts/generated-results/elliptic-curves/kihara_rank14_t2_v1.json)
make this replay exact.  It is not a rank-30 candidate and no fifteenth
independent section is claimed.

## Closest search precedents

- Elkies--Watkins,
  [*Elliptic curves of large rank and small conductor*](https://arxiv.org/abs/math/0403374),
  searches minimal models and integral points while tracking conductor and
  Szpiro ratio.
- Fermigier's E22 search starts from a Mestre family of generic rank at least
  12 and uses staged Nagao-score cutoffs.  This is an important precedent for
  deterministic rare-event levels.
- Elkies--Klagsbrun,
  [*New Rank Records for Elliptic Curves Having Rational Torsion*](https://arxiv.org/abs/2003.00077),
  precomputes local tables, sieves rational parameters, uses staged cutoffs and
  skew regions, and explicitly asks for Bayesian rank probabilities and better
  treatment of bad primes.  Those components are therefore precedent, not new
  contributions here.
- Elkies's [K3/high-rank lectures](https://arxiv.org/abs/0709.2908) describe
  the rank-17 fibration behind the rank-28 search.  The current programme has
  an exact Fermigier rank-at-least-12 adapter and a public Kihara
  rank-at-least-14 fallback, but not the unpublished rank-17 Weierstrass
  fibration and sections behind the current rank-29 record.
- Bai--Brent--Thome's
  [number-field-sieve root optimization](https://arxiv.org/abs/1212.1958)
  is the closest cross-field analogue: it uses roots modulo small prime powers,
  Hensel lifting, and sublattice sieving to improve binary-form values.

In the primary sources checked, we did not find the full combination of
deliberate discriminant prime powers, CRT selection, and shortest rational
reconstruction used as the core of a rank/conductor record search.  This is a
literature-search observation, not a novelty theorem.

## Heuristic context

Park--Poonen--Voight--Wood's
[random-matrix model](https://arxiv.org/abs/1602.01431) heuristically predicts
only finitely many curves over \(\mathbb Q\) of rank greater than 21.  That is
motivation for rare-event methods, not evidence that any particular score or
candidate has high rank.

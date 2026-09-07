# A new compact rank-at-least-26 elliptic curve

The compact R17 family `07ca9` at parameter `-2507/3068` gives the global
minimal integral equation

```text
y^2 + x*y + y = x^3 - x^2
  - 27129978154374720446538711466583970429714689912*x
  + 1828432266054662976336361254643265057329703273521399600629500669704699.
```

Its **26 independent rational points** are contained in the
[standalone proof](../../artifacts/generated-results/elliptic-curves/new_compact_rank26_proof_v1.json)
and [Sage curve file](../../artifacts/generated-results/elliptic-curves/new_compact_rank26_curve.sage).
The stable inventory ID is `new-20260905-37`. The
[Python checker](../cas/certify_new_compact_rank26_v2.py) verifies the rank
lower bound, coordinate transport, global minimality and finite novelty
comparison without Sage, numerical heights, factorization or point search.

There is no rational-isomorphism match in the pinned 586-equation
[ICARM catalogue](https://elliptic-rank.icarm.cloud/curves), or in 275
previously retained equation records. The catalogue was fetched on
2026-09-05 at 18:57:46 UTC; its SHA256 is
`1ec915b1d108f906791f5361f8150d328ce96e5f41d95d9e78d9a354e175e53a`.
ICARM is supported by NSF grant DMS2425401. This is a finite comparison,
not proof of universal novelty. Exact rank and conductor remain unknown.
The new rank-at-least-28 and rank-at-least-32 targets remain open.

## Exact proof and usable model

The discovered short model has

```text
A = -434079650469995527144619383465343526875435038587/16
B = 58509832513749215242763343108759246836786932442993054548380583713031083/32.
```

The isomorphism from the displayed integral model is
`x_short = X - 1/4`, `y_short = Y + (X+1)/2`. The checker verifies every
point on the integral model and its image on the short model. The finite
mod-2 quotient columns of the 26 short-model points are independent, and
a separate good-prime witness excludes rational 2-torsion. Infinite descent
on any integral relation proves independence; the rational group isomorphism
transports this proof to the integral model.

For that integral model, `gcd(c4,c6)=27`. Nonminimality at a prime would
require its fourth power to divide `c4` and its sixth power to divide `c6`,
hence its fourth power to divide their gcd. No fourth prime power divides
27. This proves global minimality without a discriminant factorization.
The first attempted exporter called a short-model-only rank helper on the
integral equation and failed explicitly. Its source and failed log are
retained; version 2 checks independence on the short model and proves the
transport exactly. No certificate was emitted by the failed version.

```sh
python3 elliptic-curves/cas/certify_new_compact_rank26_v2.py --check \
  artifacts/generated-results/elliptic-curves/new_compact_rank26_proof_v1.json
python3 elliptic-curves/cas/export_new_compact_rank26_sage.py --check
```

## The balanced six-family experiment

This addresses two previously unrun machinery gaps: fresh generic parity
censuses for all six compact R17 presentations, and a balanced wider
parameter search. It uses the [exact six-family atlas](COMPACT_SIX_R17_ATLAS_2026-09-05.md)
with its 102 generic sections. The earlier prospective rank-25 yield supplies
the scheduling gate; known-record points and parameters are excluded from
candidate selection and point search. The initial protocol's reference to
cross-family incidence was clarified before launch: that incidence exclusion
applies to the **first 32 measured curves**, not every future fibre.

Each fresh census visits all `2^17` generic parity classes. Across six
families, **786,432** integral representatives, parity identities and exact
quadratic norms passed independent replay, including positive definiteness
of every Gram matrix by exact leading principal minors. Every parity class
has a recorded representative of generic norm at most 12. The 43 selected
mask sets agree with the historical sets; that comparison was made only
after the fresh censuses finished. The computed representative norms are
not certified coset minima. Floating CVP output does not establish exact
covering radii or rank jumps.

The selector scores all **122,400,468** signed primitive nonzero addresses
with `abs(n), d <= 4096`, equally across the six families. Every score uses
all 562 primes from 5 through 4093. Four addresses per family are fixed
before any point search or catalogue comparison. All 24 receive fresh
attempts, including earlier addresses; there is no replacement of known
curves. Each worker has 300 seconds and 1.5 GiB, at most four concurrently.
The 43 generic charts per curve use height 100,000 and four seconds each.
Numerical canonical heights explicitly use PARI precision 384 bits and only
schedule searches; exact certificates establish the rank lower bounds.

All **1,032** retained charts and admission histories replayed successfully.
All **24 complete point clouds**, totalling **6,845** point occurrences up
to sign, passed independent finite-reduction audits through prime 997.
These audits did not raise any worker's lower bound. All 24-result rank
certificates and their post-batch equation comparisons passed standalone
replay. Chart completion means the declared bounded attempt completed;
it does not imply every chart's entire height box was exhausted.

The independent selection audit recomputes all 24 selected scores using
every retained prime and checks 144 finite-field traces by direct point
counts, at primes 5, 7, 11, 13, 997 and 4093. A separate totient sieve
checks the population size. It does not independently re-enumerate the
entire score ordering; the retained scanner transcripts provide that part
of the computational evidence.

The three additions with lower bound at least 22 are:

| Stable ID | Family | Parameter | Certified lower bound |
|---|---|---|---:|
| `new-20260905-37` | `07ca9` | `-2507/3068` | 26 |
| `new-20260905-38` | `074d9` | `-1525/1388` | 25 |
| `new-20260905-39` | `07ca9` | `-951/1696` | 25 |

The [batch certificate](../../artifacts/generated-results/elliptic-curves/compact_r17_wide_results_v1.json)
retains all 24 measurements, known matches and earlier-equation matches.
The [stable-ID inventory](../../artifacts/generated-results/elliptic-curves/new_high_rank_curve_index_v3.json)
now contains 39 curves: one with certified lower bound 26, five with 25,
eight with 24, eleven with 23 and fourteen with 22. These buckets count each
curve once. All 39 were rechecked against the same 586-equation snapshot.
The [CSV](../../artifacts/generated-results/elliptic-curves/new_high_rank_curve_index_v3.csv)
provides their equations.

The [evidence manifest](../../artifacts/generated-results/elliptic-curves/compact_r17_wide_evidence_v1.json)
and adjacent ZIP retain census, selection, chart, proof and source inputs.
The separately bounded adaptive follow-up on the new rank-26 curve has its
own protocol and checkpoints under
`artifacts/local/elliptic-curves/compact-r17-new26-followup-v1/`.
It targets two further independent directions using the nine discovered
directions beyond the generic rank-17 group. Its finite outcome is separate
from the completed initial experiment and the rank-26 proof above.

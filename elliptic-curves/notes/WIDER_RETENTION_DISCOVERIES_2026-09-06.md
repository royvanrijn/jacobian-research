# Nineteen more curves from wider candidate retention

The fixed wider-retention experiment adds **19 catalogue-unmatched curves**,
including two with 27 exactly independent rational points and three with 26.
The [89-curve inventory](../../artifacts/generated-results/elliptic-curves/new_high_rank_curve_index_v7.json)
and [equation CSV](../../artifacts/generated-results/elliptic-curves/new_high_rank_curve_index_v7.csv)
preserve all seventy earlier IDs. Their lower-bound distribution is
5×27, 8×26, 21×25, 19×24, 19×23 and 17×22. These are lower bounds, not exact ranks.

Both new 27-point curves have proved global minimal integral equations in
the convention `y² + a1*x*y + a3*y = x³ + a2*x² + a4*x + a6`:

| ID | Family | Parameter | Global minimal coefficients `[a1,a2,a3,a4,a6]` |
|---|---|---|---|
| new-20260906-71 |103b2|3726/881|`[1,0,0,-206146472796129715074018059573733454563478276563335,1105311419338647732650534549181813400819706534279793673221000968292289116297]`|
| new-20260906-72 |11952|2012/211|`[1,0,0,-63732593342850967458504535648490033346020923285,7968611834872631407039058304455964638234186241293355884337746122912097]`|

For both, `gcd(c4,c6)=1` excludes nonminimality at every prime. The
[minimal-model certificate](../../artifacts/generated-results/elliptic-curves/retention_high_rank_minimal_proofs_v1.json)
replays exact point transports and finite-reduction independence. Ready-to-load
Sage files retain the [103b2 points](../../artifacts/generated-results/elliptic-curves/new_retention_rank27_curve_103b2.sage)
and [11952 points](../../artifacts/generated-results/elliptic-curves/new_retention_rank27_curve_11952.sage).

## What changed in selection

The same H4096 parameter box contains 122,400,468 signed addresses across six
families. The earlier pipeline saved128 finalists per sign and retained128
overall per family. The [previous discarded-shard experiment](TRANSLATED_HEIGHT_AND_DISCARDED_SHARDS_2026-09-06.md)
showed that this cutoff lost useful candidates. A predetermined single-shard
benchmark retained512 per sign, preserved the old128 prefix exactly, and
recomputed every returned score from all562 canonical trace tables.

The subsequent twelve-shard run retained6144 addresses, again checking every
score and all twelve old128 prefixes. It reused one benchmark shard and took
119.467 seconds for the other eleven. The parameter height, scanner and table
bytes were unchanged. Complete population ranking still trusts the pinned
scanner; replay independently checks the returned addresses and their scores.

Excluding all1536 previously extended addresses left a fixed4608. Each received
5978 further prime traces, a total27,546,624. The score through32749 selected
exactly four per family; primes32771–65521 stayed outside selection, including
tie-breaking. The run took416.336 seconds and its exact roster/score replay
21.039 seconds. A premature launch before protocol preparation failed without
creating any score output; its original log and sequencing-correction record
remain alongside the successful separately named launch.

All24 finalists started from their17 generic sections and received every43
or49 maximum generic parity class at height125000. The budget came from the
separately completed known-control calibration. All1080 boxes completed in
661.229 seconds with two workers; all admission/archive histories and full
point-cloud proofs passed in187.257 seconds. The
[batch certificate](../../artifacts/generated-results/elliptic-curves/retention24_r17_results_v1.json)
retains every attempt, including the three catalogue-unmatched bounds19 or21
below the inventory threshold and the known-curve rediscoveries.

One illustrative cutoff loss is `08f72` at `1103/375`: it was721st by short
score among the1024 retained family addresses but first by the extended score.
It blindly supplied27 independent points and was then identified as known
ICARM curve363. It is a control rediscovery, not a new curve. The other known
match is394. No catalogue or published points entered selection or searches.
Two new27 curves and nineteen additions demonstrate productive finite output;
they do not establish a general score-to-rank theorem or close new rank28/32.

## Fresh comparison and remaining boundaries

The public [ICARM database](https://elliptic-rank.icarm.cloud/) was fetched again
at `2026-09-06T02:16:21.873568+00:00`:593 equations, SHA256
`4ab1b29bbcb90ba329c280b26924313685a98b30891baefafd98dfa4f25199ed`.
The [exact refreshed comparison](../../artifacts/generated-results/elliptic-curves/retention_refreshed_catalogue_comparison_v1.json)
finds no rational-isomorphism match for any of the89 inventory equations.
The frozen original586-row comparisons remain intact. ICARM's highest reported
lower bound is31; its public metadata is not independently recertified here.
Catalogue absence never proves that no one has seen a curve before.
Acknowledgement: ICARM and NSF Grant DMS2425401.

The new nineteen also complete228 exact incidence checks against the same
twelve presentations. There are206 rational-preimage exclusions and22
complete preimage factorizations. Three extra publishedR17 presentations are
the proved08234 equivalence `s=-26t-50`, with the same integral generic subgroup.
The [aggregate binding](../../artifacts/generated-results/elliptic-curves/inventory89_cross_family_incidence_v1.json)
now covers1068 pairs:958 exclusions,89 original presentations and21 duplicates.
No additional recorded generic directions are supplied. Other families and
nongeneric rational points remain open.

A separate98-word retrospective diagnostic translates the known control's
last published direction by its blindly recovered28-point subgroup. Exact
group and chart replay still gives918522 as the smallest height among those
proposals, outside the original100000 boxes. This is no global minimum or
prospective point oracle, and it triggers no automatic height expansion.

The [portable evidence supplement](../../artifacts/generated-results/elliptic-curves/retention24_discovery_evidence_v1.json)
retains all raw traces, charts, state archives, protocols, failed launches,
source dependencies and exact proofs. Its isolated verifier declares39 stages
covering the scores,89 point certificates andCSV, all1068 incidence checks,
generic section transport and refreshed catalogue comparison. Local admission
histories are retained and already passed, without claiming an isolated repeat.
Separate301-centre follow-ups on the two new27 subgroups have their own fixed
protocols; they do not alter this initial experiment or its certificates.

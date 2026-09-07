# Six additional compact-atlas curves, including rank at least 25

The balanced pilot on the [six compact R17 models](COMPACT_SIX_R17_ATLAS_2026-09-05.md)
produced six additional curves with exact rank lower bounds 22–25. All have
distinct `j`-invariants, also distinct from the fifteen earlier compact-R17
curves, and no rational-isomorphism match in the pinned 584-curve ICARM
snapshot. The combined discovery result is twenty-one curves with lower
bounds at least 22, including one at least 25. The subsequent
[MW16 searches](PROSPECTIVE_COMPACT_MW16_PILOT_2026-09-05.md) add eleven more,
including two further rank-at-least-25 examples. No new curve yet meets the
rank-at-least-28 near-record or rank-at-least-32 record target.

| Compact family | New compact parameter | Certified rank at least | Gain beyond the seventeen generic sections |
|---|---|---:|---:|
| `07ca9` | `505/794` | 25 | 8 |
| `08f72` | `-704/93` | 24 | 7 |
| `074d9` | `-839/505` | 23 | 6 |
| `08234` | `-506/9` | 23 | 6 |
| `074d9` | `148/225` | 22 | 5 |
| `11952` | `744/5` | 22 | 5 |

Parameters refer to the **new compact coordinate in the atlas**, not the old
native chart or the parameter of the published compact equation used in the
earlier fifteen-curve experiment. The atlas records every exact base map.

The [point certificates](../../artifacts/generated-results/elliptic-curves/compact_atlas_new_curves_v1.json)
contain all models, independent points, finite quotient matrices, family
transports and comparison equations. Replay without Sage or discovery caches:

```sh
python3 elliptic-curves/cas/certify_compact_atlas_candidates.py --check \
  artifacts/generated-results/elliptic-curves/compact_atlas_new_curves_v1.json
```

## The rank-at-least-25 curve

A convenient integral equation for the `07ca9`, `505/794` specialization is

```text
y^2 + x*y + y = x^3 - x^2
 - 1059556401049745356603514361839762881637078027*x
 + 14325193788023719095355453992974755532155118126141941235002031329435.
```

From its short certificate model `(X,Y)`, set `x=X+1/4` and
`y=Y-(x+1)/2`. All 25 certified points transport to this displayed equation.
No global minimality, exact rank or conductor assertion is needed.

The separate Sage-free checker recomputes membership and full column rank
in a retained product of `E(F_p)/2E(F_p)`, and verifies a modular witness for
`E(Q)[2]=0`. Every integral relation must then have even coefficients;
division by two and infinite descent prove independence. The first seventeen
points also match the exactly specialized atlas sections after their common
rational scale and sign. This uses no numerical height or analytic-rank
assumption. The discovery and replay implementations share group-law
primitives; external or formal verification is not claimed.

## Completed initial experiment and limits

Each family scored 1,275,854 signed primitive height-1024 parameters with all
562 primes from 5 through 4093 before retention. Four finalists per family
were frozen. Five catalogue matches were skipped without refilling. All
nineteen fresh workers completed their declared 43-chart plans, with measured
lower-bound counts `17:5, 19:2, 20:4, 21:2, 22:2, 23:2, 24:1, 25:1`.
The chart height was 100000 and the allowance four seconds per chart;
completed plans are not assertions that all their boxes were exhaustive.

The source, protocols, population files, exact post-freeze deduplication and
worker results are retained in
`artifacts/local/elliptic-curves/compact-six-r17-h1024-v2/`.
The independent point export uses only terminal initial measurements.
The trace-table generation was checked against Sage cardinalities at every
projective residue for primes 5, 7, 11 and 13 on all six families.
All nineteen initial transcripts, comprising 817 chart and admission records,
passed replay without resieving. The
[replay manifest](../../artifacts/generated-results/elliptic-curves/compact_six_r17_initial_replay_v1.json)
and [portable bundle](../../artifacts/generated-results/elliptic-curves/compact_six_r17_initial_replay_v1.zip)
retain the complete result transcripts and selection, launch and verification
ledgers. Each result replays with `replay_compact_atlas_search.py RESULT.json`.

A separate 301-chart follow-up starts with the new rank-25 subgroup, using
only points recovered by the prospective search. It retains the earlier
four-second chart allowance, one 1500-second/1.5-GiB worker and the rank-28
stop target. Its first startup accidentally used a fixed 17-dimensional
oracle and was rejected before any chart ran. That failed attempt is
preserved as `compact-atlas-first25-followup-v1`. The corrected v2 uses the
existing general-dimension geometry and checks every centre's length and
parity; exact norm regressions at dimensions 25 and 28 passed. It hit the
1500-second limit after 262 charts. A separately declared five-minute
continuation completed the remaining 39 frozen charts, with no further gain.

The three other high-gain follow-ups retained 257 charts on `08f72,-704/93`,
215 on `074d9,-839/505` and 153 on `08234,-506/9`, keeping lower bounds
24, 23 and 23. The first two hit their 1500-second wall limits; the third
hit its 1.5-GiB RSS limit. All four retained transcripts, 926 charts in total,
passed exact replay without resieving. The last replay required a separately
declared 600-second attempt after its original 300-second allowance expired.
Only the rank-25 follow-up completed its 301-chart plan. The
[follow-up evidence manifest](../../artifacts/generated-results/elliptic-curves/compact_atlas_followup_replay_v1.json)
and [portable bundle](../../artifacts/generated-results/elliptic-curves/compact_atlas_followup_replay_v1.zip)
preserve source, protocols, startup failure, stopped snapshots, continuations
and verification outcomes.

Separate exact mod-3 and mod-5 audits of the four initial transcripts' unique
ambiguous points retained ranks 23/23, 25/25, 23/23 and 24/24 respectively
for `074d9,-839/505`, `07ca9,505/794`, `08234,-506/9` and `08f72,-704/93`.
All four certificate replays passed. The audit tested 567, 442, 1063 and 433
points up to sign; it found no additional independent direction and proves
neither dependence nor a rank upper bound.

A bounded conductor preflight on the displayed rank-25 equation removed all
discriminant prime factors through 10000 and retained a 318-bit composite
cofactor. This exceeded its declared 192-bit composite-factorization gate;
the conductor remains unknown. This failed gate is included in the bundle.
The [separate MW16 pilot](PROSPECTIVE_COMPACT_MW16_PILOT_2026-09-05.md) uses
the five newly compacted inputs to broaden the prospective population.

## Novelty boundary

The ICARM comparison uses the 584-equation snapshot with SHA-256
`7e80549befa11a07422a3960967f4cd80264d8675cb3e0a99f0c9c5afb340f72`,
dated 2026-09-05. The certificate also retains the equations of all fifteen
earlier prospective curves and excludes rational isomorphisms to them.
All twenty-one exact `j`-invariants were compared across the five exports.
Absence from these finite sets does not prove absence from every publication
or unpublished search. This work used ICARM, supported by NSF Grant DMS 2425401.

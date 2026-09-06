# All seventy curves checked against the twelve recorded presentations

The [23-curve extension](../../artifacts/generated-results/elliptic-curves/latest23_cross_family_j_incidence_v1.json)
closes the missing incidence coverage in the
[70-curve inventory](../../artifacts/generated-results/elliptic-curves/new_high_rank_curve_index_v6.json).
Its 276 pairs comprise 251 exact rational-preimage exclusions and 25 complete
rational-root factorizations. No case is unresolved and no target occurs at
parameter infinity. Building took 1.963 seconds; exact rational replay took
0.379 seconds. These are finite, exact incidence tests, not point searches.

The two extra presentations are:

| Curve | Compact08234 parameter s | PublishedR17 parameter t |
|---|---:|---:|
| new-20260906-56 | 1504/3 | -827/39 |
| new-20260906-58 | -2387/360 | -1201/720 |

Both obey `s = -26t - 50`. The existing
[generic transport proof](COMPACT_CROSS_FAMILY_INCIDENCE_2026-09-05.md)
identifies the models by scale26 and all seventeen sections by an integral
matrix of determinant -1. Its exact function-field replay passes without
numerical heights. Thus these presentations give the same generic subgroup;
they supply no additional independent directions. Equal j alone would not
have established this conclusion.

The [aggregate certificate](../../artifacts/generated-results/elliptic-curves/inventory70_cross_family_incidence_v1.json)
binds all four original cohorts to the same current equations and family
models: **384 + 84 + 96 + 276 = 840 = 70 × 12** checks, with 752 exclusions,
70 original presentations and 18 duplicates of the proved08234/publishedR17
equivalence. The twelve presentations consist of six compactMW17, five
compactMW16 and publishedR17. The four incidence replays and the generic
transport replay remain separate proof obligations.

The self-contained [evidence archive](../../artifacts/generated-results/elliptic-curves/inventory70_incidence_evidence_v1.json)
retains those obligations and their inputs. Its
[isolated replay](../../artifacts/generated-results/elliptic-curves/inventory70_incidence_portable_replay_v1.json)
checks all four cohorts, the function-field transport and the aggregate binding.
It does not replay the separate independent-point certificates.

This closes this route to an additional **recorded generic subgroup** on these
seventy curves. Other families, future fibres, nongeneric points, saturation
and rank upper bounds remain open. The 27-point lower bounds on the strongest
new curves are unchanged.

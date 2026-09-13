# The201-curve model and conductor audit —7 September2026

This is the frozen201-row audit. Its counts and table layout are historical.
The [current inventory](../INVENTORY.md), [JSON](../data/research_curves/database.json)
and [CSV](../data/research_curves/database.csv) select later certified additions;
this note does not prescribe the current database selection or rendering.

## Exact models and numerical columns

The [metric certificate](../../artifacts/generated-results/elliptic-curves/inventory201_table_metrics_v1.json)
checks201 global minimal models and4796 exact point transports. Local minimality
is checked independently at2,3 and the prime divisors of the invariant gcd,
which is at most10^7. Outside that set an invariant is a unit, excluding a
smaller integral model. Exact isomorphisms transfer the original independence
certificates to the displayed models; complete discriminant factorization is
not required. The rank28 row reproduces public points, and seven rows match
the pinned ICARM catalogue.

The historical table uses natural logarithms. Naive height is
`log max(abs(c4)^3,c6^2)` on the minimal model. Its period-area Faltings value
is `-1/2 log(area)`, matching the ICARM convention and Sage `stable=False`;
this differs from the stable-height option on additive curves. Numerical
columns were checked at96 and160 bits, displayed to two decimals, and matched
all seven public rows within10^-9. Full exact integers and transported points
remain in the certificates. Conductor divisors and upper bounds occupy separate
fields; an unresolved exact conductor stays null and a partial prime list is
not a complete submission list.

## Completed conductor continuation

The [201-row supplement](../../artifacts/generated-results/elliptic-curves/inventory201_conductors_v2/summary.json)
contains129 exact conductors,100 more than the earlier29. The other72 remain
UNKNOWN in this batch. Both public rows600 and619 gain local exact proofs.
All129 exact values exceed their recorded rank-threshold benchmark; this pass
proves no conductor record. Large upper bounds do not exclude the unresolved curves.

The frozen protocol covers172 rows without a local exact proof, including two
with previously reported public values. Each receives at most15 seconds of
PARI partial factoring and30 seconds of deterministic ECM, with four workers
and per-attempt checkpoints. A repository relocation interrupted the pass;
retained checkpoints resumed without restarting completed curves. Prior factors
were candidate inputs, checked by product. Exact prime certificates and
independent Sage generic Tate/PARI local exponents were required for completion.
No point search or rank change occurred. The [earlier24-curve screen](NEW_CURVE_CONDUCTOR_RECORD_SCREEN_2026-09-07.md)
and the [four submission certificates](SUBMITTED627_630_BAD_PRIMES_2026-09-07.md)
remain distinct evidence.

## Replay and historical commands

`certify_inventory201_conductors_v2.sage --check` reuses saved factors and prime
certificates but repeats primality-certificate validation and both local
reduction implementations on all201 rows. Its summary also reads the retained
`artifacts/local/elliptic-curves/conductor-inventory-continuation-v2/factor_protocol.json`.
It performs no factor search; missing inputs do not authorize reconstruction.

The [original note](../../archive/elliptic-curves/notes/INVENTORY201_TABLE_AND_CONDUCTORS_2026-09-07.before-2026-09-13.md.txt)
preserves the old database-selection and README-rendering commands. Those
commands can select the historical201-row batch; use the [current replay guide](../../REPRODUCE.md)
for maintained entry points. Further factorization requires a new bounded
protocol and versioned certificates; the completed pass retains partial factors
and exhausted attempts.

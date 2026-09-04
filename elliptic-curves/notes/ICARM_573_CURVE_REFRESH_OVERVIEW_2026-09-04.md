# ICARM 573-curve refresh and missing-curve overview

<!-- status-consumer: EC-K3-R17-NORM12-ICARM-573-REFRESH a93ce35de34fde21 -->

Status: **exact hash-pinned ids-1--573 norm-twelve atlas sweep; exact
ids-475--573 equation/point/discriminant overview; exact displayed-point
independence and MW17 specialization audit for all seventeen new atlas
hits**.

## Outcome

The 2026-09-04 ICARM response has SHA-256
`e57d991894722f0e5ab2f548b77f09064a46ec926c93ef3730f47685e016aab0`
and contains ids 1 through 573.  The original 474-curve artifact remains an
immutable historical theorem.  A new v2 sweep decides all `573*6=3,438`
curve/class pairs against the same six exact rational-`PGL2` norm-twelve
`j`-classes.

There are 86 hits and 3,352 misses.  Relative to the 474-curve snapshot, 17 of
the 99 appended curves are new hits and 82 miss all six classes.  All 103 new
native chart/fibre comparisons are untwisted over `QQ`; together with the old
results this gives 479 untwisted comparisons.

| representative | old hits | refreshed hits | new hit ids |
| --- | ---: | ---: | --- |
| `074d9` | 5 | 6 | 538 |
| `08234` | 54 | 61 | 478, 531, 534, 535, 536, 537, 546 |
| `0e80b` | 2 | 4 | 539, 541 |
| `11952` | 2 | 4 | 532, 540 |
| `07ca9` | 3 | 7 | 499, 543, 544, 545 |
| `08f72` | 3 | 4 | 498 |

The companion overview retains the complete public projection for ids 475
through 573.  All displayed points were checked exactly on their equations,
and all 99 stored discriminants were recomputed from the displayed
Weierstrass coefficients.  Commentary-only provenance tags give 36 declared
Elkies--Klagsbrun rank-9-with-2-torsion fibres, six declared published-R17
fibres, three declared R17 two-neighbour fibres, two declared
Mestre--Fermigier fibres, five other commented records, and 47 records with no
public construction commentary.  These tags organize the intake; they do not
independently prove provenance.

## Incorporated atlas hits and priority misses

Every appended curve of rank at least 24 received an independent exact rank
lower-bound replay and a trivial-torsion certificate.  Twelve close with
products of `E(F_p)/2E(F_p)`.  Curve 542's displayed subgroup is not
2-saturated enough for that lightweight proof: the mod-2 image has rank 25,
but the declared exact mod-3 fallback has rank 26 and closes all 26 points.

Eleven of the thirteen priority curves are exact untwisted atlas fibres.  The
specialized saturated generic MW17 basis was recovered inside every displayed
public subgroup, all relations were rechecked by exact group law, and Smith
normal form proves the following displayed-subgroup quotients:

| curve | rank at least | native class/chart | displayed quotient by MW17 |
| ---: | ---: | --- | --- |
| [543](https://elliptic-rank.icarm.cloud/curve/543) | 29 | `07ca9` | `Z^12` |
| [531](https://elliptic-rank.icarm.cloud/curve/531) | 28 | `08234` | `Z^11` |
| [534](https://elliptic-rank.icarm.cloud/curve/534) | 28 | `08234` | `Z^11` |
| [535](https://elliptic-rank.icarm.cloud/curve/535) | 28 | `08234` | `Z^11` |
| [536](https://elliptic-rank.icarm.cloud/curve/536) | 28 | `08234` | `Z^11` |
| [544](https://elliptic-rank.icarm.cloud/curve/544) | 28 | `07ca9` | `Z^11` |
| [545](https://elliptic-rank.icarm.cloud/curve/545) | 28 | `07ca9` | `Z^11` |
| [537](https://elliptic-rank.icarm.cloud/curve/537) | 27 | `08234` | `Z^10` |
| [540](https://elliptic-rank.icarm.cloud/curve/540) | 25 | `11952` | `Z^8` |
| [541](https://elliptic-rank.icarm.cloud/curve/541) | 25 | `103b2` in class `0e80b` | `Z^8` |
| [546](https://elliptic-rank.icarm.cloud/curve/546) | 25 | `08234` | `Z^8` |
| [498](https://elliptic-rank.icarm.cloud/curve/498) | 23 | `08f72` | `Z^6` |
| [539](https://elliptic-rank.icarm.cloud/curve/539) | 23 | `103b2` in class `0e80b` | `Z^6` |
| [538](https://elliptic-rank.icarm.cloud/curve/538) | 22 | `074d9` | `Z^5` |
| [478](https://elliptic-rank.icarm.cloud/curve/478) | 21 | `08234` | `Z^4` |
| [532](https://elliptic-rank.icarm.cloud/curve/532) | 20 | `11952` | `Z^3` |

For every quotient row the simple public complement is `P18,...,P_r`.  Thus
the public statement that curve 543 consists of seventeen generic directions and twelve
covering directions is now an exact internal specialization theorem for the
displayed subgroup.  It is not an assertion that the displayed subgroup is
the full Mordell--Weil group.

Curve 499, the remaining lower-rank hit, has a different exact outcome.  Its
twenty displayed points are independent and the generic MW17 specialization
has rank 17, but the latter is not contained in the displayed subgroup.  If
`G_j` denotes the specialized generic basis and `P_i` the displayed basis,
the certificate verifies exact relations `3 G_j = sum_i N_ij P_i`; the matrix
`N mod 3` has rank one.  Consequently adjoining the generic subgroup enlarges
the displayed subgroup by `Z/3Z`.  The expression "displayed subgroup modulo
MW17" is therefore not defined for curve 499, and the artifact records the
commensurability obstruction rather than inventing a quotient.

The other two priority curves are exact six-class misses:

- curve 542 has rank at least 26 by the exact mod-3 fallback and no public
  construction commentary;
- curve 548 has rank at least 24 by the exact mod-2 replay and no public
  construction commentary.

Their absence is only from this 43-chart norm-twelve atlas, not from every K3
or high-rank family.

## Remaining appended curves

The other 82 appended curves have one exact no-root decision for each of the
six classes.  Their public point lists were checked for membership, but the
overview does not independently reproduce point independence below the
priority threshold unless they are one of the six atlas hits audited above.
Forty-seven of all 99 appended rows have no public
construction commentary, so those rows remain construction `UNKNOWN` rather
than being assigned to a family by appearance.

The public conductors are absent for curves 537, 543, 545, and 568.  A bounded
Sage factorization attempt on curve 537 was stopped after 90 seconds without a
factorization, so no conductor or local-reduction claim is made for those four
rows.

## Replay

```bash
sage -python elkies-k3/scripts/certify_r17_norm12_icarm_database_sweep.sage \
  --refresh-573 --check

.venv/bin/python \
  elliptic-curves/cas/audit_icarm_curve_refresh_overview.py --check

PYTHONPATH=elliptic-curves/cas sage -python \
  elkies-k3/scripts/certify_r17_norm12_refresh_priority_quotients.sage --check
```

The first two artifacts embed the claim-sufficient source projections, so the
default replays are offline.  Reconstructing them from the captured raw ICARM
response requires its exact byte hash.

## Claim boundary

This refresh proves exact snapshot-bounded recognition, exact priority rank
lower bounds, sixteen displayed MW17 quotients, and the curve-499
commensurability obstruction.  It proves no exact
rank upper bound, no full Mordell--Weil group, no construction for the 82
atlas misses, no completeness beyond the six norm-twelve classes, and no
geometric origin for every exceptional quotient direction.

# Elliptic curves over `Q`

This is the active arithmetic programme for exceptional Mordell--Weil rank and
small conductor. The original target was to find either

- `rank E(Q) >= 21` with `log N(E) < 182.72`; or
- `rank E(Q) >= 30`.

The second branch is proved and has since been improved to rank at least 31.
The first branch is now completely replayed, and the current repository-local
rank-at-least-21 conductor anchor is ICARM curve 394 at
`log(N)=166.252098...`.

[`MATH_STATUS.json`](../MATH_STATUS.json) is the sole status authority. This
page is a map of that record, not a second status database.

## Current status

| Curve or family | What is proved here | Boundary |
| --- | --- | --- |
| ICARM curve 302 | `rank E(Q) >= 31`, trivial torsion, global minimality, exact discriminant/conductor and local fibre data; two exact independence implementations | No unconditional rank upper bound; no K3-family identification |
| ICARM curve 356 | `rank E(Q) >= 29`, trivial torsion, global minimality, exact conductor and local fibre data; current size record at the rank-at-least-29 threshold | No unconditional rank upper bound; strong common-17-section fingerprint with curve 351, but no family or K3 identification |
| ICARM curve 273 | `rank E(Q) >= 30`, independently replayed | No unconditional exact-rank statement |
| ICARM curves 285 and 286 | Twenty-one displayed points on each curve are independent; trivial torsion; global minimality and exact conductor reconstructed from local Tate data | Rank lower bounds only; no unconditional upper bounds |
| ICARM curve 394 | Compact Elkies `t=3/8` fibre; generic seventeen plus four public directions certify `rank E(Q) >= 21`; exact local conductor replay gives `log N=166.252098...` | Rank lower bound only; still above curve 245's `150.668907...` conductor line |
| Fermigier `E22` | `rank E(Q) >= 22` | `log N=182.724910...`, so it misses the strict cutoff |
| ICARM curve 245 | `rank E(Q) >= 20` and exact `log N=150.668907...` | One point short; no rank upper bound |
| Fermigier--Mestre `u=28917/20` | `rank E(Q) >= 20` and exact `log N=159.934825...` | One point short; no rank upper bound |
| Split-infinity Mestre frontiers | Two exact rank-at-least-19 curves below the cutoff | Exact rank 19 only conditionally under the hypotheses stated in the note |
| New six-root family at `T=83/6` | Exact rank 14 from a rank-14 subgroup and PARI interval `[14,14]` | A calibration specialization, not a target curve |

The canonical statements are:

- `ECR31`: [curve 302 rank-at-least-31 certificate](notes/ICARM_CURVE302_RANK31.md);
- `ECR29-IC356`: [curve 356 rank-at-least-29 size record and construction fingerprint](notes/ICARM_CURVE356_RANK29_AND_CONSTRUCTION.md);
- `ECR30`: [curve 273 rank-at-least-30 certificate](notes/ICARM_CURVE273_RANK30.md);
- `EC-R21-ICARM`: [curves 285/286 point-independence replay](notes/ICARM_7FFF_ZIP_SEQUENCE.md);
- `EC-R21-ICARM394`: [curve 394 compact specialization and conductor replay](notes/ICARM_CURVE394_RANK21.md);
- `EC-R20-IC245`: [curve 245 low-conductor rank-20 replay](notes/ICARM_CURVE245_RANK20.md);
- [conductor-first descent inputs for the four rank-19/20 near misses](notes/CONDUCTOR_FIRST_NEAR_MISS_DESCENT.md);
<!-- status-consumer: EC-CF-NEARMISS-DESCENT-INPUTS 25c9f212e5162216 -->
- `EC-NF-R14`: [new-family exact rank-14 specialization](notes/NEWFAMILY_RANK14_T83_6.md);
- `OP-EC-NEXT`: rank 32, exact-rank, low-conductor certification, and K3 continuation.

## Evidence labels

These labels are deliberately non-interchangeable.

| Label | Meaning |
| --- | --- |
| **Exact lower bound** | Displayed rational points are checked and their independence is certified. It does not assert equality. |
| **Exact rank** | A certified lower bound and an unconditional upper bound agree. |
| **Conditional closure** | An upper bound uses named hypotheses such as GRH or BSD. |
| **Exact structural computation** | A symbolic identity, finite quotient, local classification, or construction is checked exactly, but it need not imply a rank claim. |
| **Bounded experiment** | Only the declared finite search was completed. A miss is not an upper bound. |
| **Public-source snapshot** | External data are pinned for provenance; they are not automatically an independent replay. |

In particular, a point list is not a rank lower bound until independence is
proved, a PARI score is not a rank estimate, and a discriminant radical is not
the conductor. See [THEORY.md](THEORY.md) for the exact independence and local
arithmetic used by the active certificates.

## Active work

There are five live arithmetic gates.

1. In the compact published Elkies `t` coordinate, retain the exact rank-25--28
   positive controls, rank by weakest performance across at least three
   disjoint prime ensembles, and compute the actual residual 2-Selmer quotient
   before any expensive point search. The complete rank-28 discriminant
   factorization and all bad-place known Kummer images are now pinned; the
   original PARI, factored-PARI, and eclib attempts still have no completed
   Selmer bound. The factor-supplied class-group worker now isolates its stall
   in relation generation. An exact reduced cubic lowers the defining-order
   index by 27 but reaches the same class-relation plateau, while the exact
   factor-base-1000 BNF-free pilot is far below a valid generation bound and
   has no noncanonical quotient gain.
   The exact rank-28 local positive control further shows that all eleven
   exceptional directions add zero bad-place signature rank. Full local-image
   coverage is certified at only four odd places and infinity; a bounded
   resumable norm-one-cover pilot certifies 60 selected local witnesses but no
   everywhere-local class or obstruction. The eleven public complement points
   are now also materialized as exact Kummer classes and two-covers with
   rational witnesses, certifying the expected residual Selmer lower bound 11
   while leaving the ambient upper bound open. None authorizes search. The
   unconditional Magma job is generated but not run locally.
   The same controls now have a reproducible
   [escape-from-R17 fingerprint](notes/ELKIES_RANK_JUMP_FINGERPRINTS.md):
   the `t=3/8` mechanism control has quotient rank 4 and the high-rank controls
   have exact displayed quotient ranks `8,9,10,11`, unit Smith saturation indices,
   complete projected-height successive minima, and degree-two visible spans
   `4` and `5,3,2,1`. Degree-three and degree-four visibility remain explicitly
   missing rather than zero.
2. Prove an unconditional upper bound for curve 302 or find a rank-at-least-32
   curve.
3. Compute residual descent quotients for the four low-conductor rank-19/20
   near misses, construct any surviving covers, and search their surrounding
   families in conductor-first order.
4. Complete the residual 2-Selmer calculation for curve 273; Selmer classes
   must not be promoted to Mordell--Weil directions without the remaining
   global argument.
5. Preserve the now-complete H3/rootless-MW17 equation transport as the exact
   construction source for the compact arithmetic programme. Alternate suffix,
   reverse-lift, and family-inference routes are parked unless they directly
   support specialization.

The low-conductor Fermigier/Mestre search remains useful, but its accumulated
negative scans are historical calibration. They are indexed in the archive
and are not part of the active search command surface.  Their recoverable
score and outcome data are now joined, with censored rather than fabricated
negative labels, by the
[`Fermigier labelled-corpus protocol`](notes/FERMIGIER_LABELLED_CORPUS.md).

The [rank-jump laboratory](notes/RANK_JUMP_LABORATORY.md) inventories
certified within-family positives and censored controls separately. Its first
executable corpus replays the complete compact-R17 height-10,000 Nagao
ranking. Its complete
[Fermigier replay](notes/FERMIGIER_RANK_JUMP_REPLAY.md) now shows that all four
frozen local orderings miss both certified quotient jumps through a
100,000-candidate budget; it is a retrospective negative retrieval result, not
a prospective holdout. The normalized
[Nagao section-7 replay](notes/NAGAO_SECTION7_RANK_JUMP_REPLAY.md) adds a
complete 18,244,819-parameter population and an exact quotient fingerprint:
free rank 8, mod-`2/3/5` dimensions `19/8/8`, and bounded degree-two searches
spanning all eight free directions. Its frozen bands rank the positive at
9,041,935 and 755,065, so this remains development evidence rather than a
search launch. Mestre and E29 remain behind explicit provenance gates.

## Repository layout

- [`scripts/`](scripts/) contains stable generation and replay entry points.
- [`cas/`](cas/) contains the exact arithmetic modules, specialized CAS
  checkers, and the few active search drivers used by the five gates above.
- [`ecsearch/`](ecsearch/) is the dependency-free exact arithmetic layer.
- [`tests/`](tests/) contains current regression tests.
- [`families/`](families/) contains normalized family equations and provenance.
- [`notes/`](notes/) contains canonical mathematical explanations.
- [`../artifacts/generated-results/elliptic-curves/`](../artifacts/generated-results/elliptic-curves/)
  contains the compact current artifact set and its evidence catalogue.
- Ignored raw runs and checkpoints belong in
  `artifacts/local/elliptic-curves/`.
- [`../archive/elliptic-curves/`](../archive/elliptic-curves/) contains the
  preserved bounded-search history: scripts, their matching tests, generated
  outputs, original commands, hashes, and old paths.

The archive is for provenance and idea recovery. It is not a current source of
mathematical status.

## Reproduction

The short command catalogue is [REPRODUCE.md](REPRODUCE.md). The principal
checks are:

```sh
make verify-elliptic-curves

PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/check_icarm_curve302_rank31_pinned.py

PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/verify_icarm_curve356_rank29.py

.venv/bin/python elliptic-curves/cas/verify_icarm_curve273_rank30.py --check

.venv/bin/python elliptic-curves/cas/analyze_icarm_7fff_zip_sequence.py --check
```

Commands requiring Sage, Singular, Magma, PARI/GP, eclib, or long searches are
identified individually in the reproduction catalogue. Do not refresh a
pinned artifact merely because a newer software version renders it
differently; record the command, version, old and new hashes, and the reason.

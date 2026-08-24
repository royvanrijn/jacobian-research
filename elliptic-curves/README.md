# Elliptic curves over `Q`

This is the active arithmetic programme for exceptional Mordell--Weil rank and
small conductor. The original target was to find either

- `rank E(Q) >= 21` with `log N(E) < 182.72`; or
- `rank E(Q) >= 30`.

The second branch is proved and has since been improved to rank at least 31.
The first branch has two public candidates whose point independence is replayed
exactly here, but whose conductor/minimality calculations have not yet been
independently reconstructed in this repository.

[`MATH_STATUS.json`](../MATH_STATUS.json) is the sole status authority. This
page is a map of that record, not a second status database.

## Current status

| Curve or family | What is proved here | Boundary |
| --- | --- | --- |
| ICARM curve 302 | `rank E(Q) >= 31`, trivial torsion, global minimality, exact discriminant/conductor and local fibre data; two exact independence implementations | No unconditional rank upper bound; no K3-family identification |
| ICARM curve 273 | `rank E(Q) >= 30`, independently replayed | No unconditional exact-rank statement |
| ICARM curves 285 and 286 | Twenty-one displayed points on each curve are independent; trivial torsion | Their sub-threshold conductors are imported public data, not yet a repository-local Tate-algorithm replay |
| Fermigier `E22` | `rank E(Q) >= 22` | `log N=182.724910...`, so it misses the strict cutoff |
| ICARM curve 245 | `rank E(Q) >= 20` and exact `log N=150.668907...` | One point short; no rank upper bound |
| Fermigier--Mestre `u=28917/20` | `rank E(Q) >= 20` and exact `log N=159.934825...` | One point short; no rank upper bound |
| Split-infinity Mestre frontiers | Two exact rank-at-least-19 curves below the cutoff | Exact rank 19 only conditionally under the hypotheses stated in the note |
| New six-root family at `T=83/6` | Exact rank 14 from a rank-14 subgroup and PARI interval `[14,14]` | A calibration specialization, not a target curve |

The canonical statements are:

- `ECR31`: [curve 302 rank-at-least-31 certificate](notes/ICARM_CURVE302_RANK31.md);
- `ECR30`: [curve 273 rank-at-least-30 certificate](notes/ICARM_CURVE273_RANK30.md);
- `EC-R21-ICARM`: [curves 285/286 point-independence replay](notes/ICARM_7FFF_ZIP_SEQUENCE.md);
- `EC-R20-IC245`: [curve 245 low-conductor rank-20 replay](notes/ICARM_CURVE245_RANK20.md);
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

There are four live arithmetic gates.

1. Prove an unconditional upper bound for curve 302 or find a rank-at-least-32
   curve.
2. Reconstruct global minimality and every local conductor exponent for ICARM
   curves 285 and 286. Their 21-point independence certificates are already
   exact.
3. Complete the residual 2-Selmer calculation for curve 273; Selmer classes
   must not be promoted to Mordell--Weil directions without the remaining
   global argument.
4. Continue the H3/rootless-MW17 equation transport documented in
   [the curve-273 construction investigation](notes/ICARM_CURVE273_CONSTRUCTION_INVESTIGATION.md).

The low-conductor Fermigier/Mestre search remains useful, but its accumulated
negative scans are historical calibration. They are indexed in the archive
and are not part of the active command surface.

## Repository layout

- [`scripts/`](scripts/) contains stable generation and replay entry points.
- [`cas/`](cas/) contains the exact arithmetic modules, specialized CAS
  checkers, and the few active search drivers used by the four gates above.
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

.venv/bin/python elliptic-curves/cas/verify_icarm_curve273_rank30.py --check

.venv/bin/python elliptic-curves/cas/analyze_icarm_7fff_zip_sequence.py --check
```

Commands requiring Sage, Singular, Magma, PARI/GP, eclib, or long searches are
identified individually in the reproduction catalogue. Do not refresh a
pinned artifact merely because a newer software version renders it
differently; record the command, version, old and new hashes, and the reason.

# Rootless J2 completeness track — 2026-08-31

## Outcome

The rootless-frame classification track is now initialized with an exact
genus and corpus audit.  It is **not complete**.

The two mandatory positive controls pass:

| frame | determinant | minimum | norm-four pairs | automorphism order |
| --- | ---: | ---: | ---: | ---: |
| published R17 | 948 | 4 | 1,311 | 2 |
| alternate Q80-derived | 948 | 4 | 1,313 | 4 |

They have the same local genus and PARI `qfisom` proves that they are not
integrally isometric.  Thus the exact lower bound remains two rootless `J2`
frame classes.

The proposed unfiltered genus traversal is not the cheap formulation for this
lattice.  The target genus is

```text
signature: (17,0)
2-adic:    1^-16:[4^1]_1
3-adic:    1^-16 3^-1
79-adic:   1^16 79^-1
```

and its exact mass is

```text
77731517730627488307787 / 925557271717478400
  = 83983.476880245...
```

Every positive-definite lattice has `+/-identity` in its automorphism group,
so each mass summand is at most `1/2`.  The full genus therefore contains at
least

```text
ceil(2*mass) = 167967
```

integral isometry classes.  A mass certificate would require accounting for
all of their reciprocal automorphism orders.  It cannot certify completeness
of the rootless subset from the two rootless controls alone.

## Required correction to the old neighbour corpus

The 65 stored files under
[`seeds/target-genus-rootless-pneighbor`](seeds/target-genus-rootless-pneighbor)
are rootless, even, rank 17, and determinant 948, but they are **not** in the
R17/Q80 target genus.  Their local symbols are

```text
2-adic:    1^16:[4^-1]_5
3-adic:    1^-16 3^-1
79-adic:   1^-16 79^1
```

Exact isometry deduplication gives nineteen classes inside that rejected
corpus.  None counts toward the present `J2` lower bound.  The directory name
and the comment in the old beam script were based on determinant/root-shell
filters rather than a full local-genus check.

This is a reusable guard: rank, determinant, parity, rootlessness, and even a
close theta-shell fingerprint do not determine the discriminant form or local
genus.  Every future neighbour or Niemeier complement must pass the exact
local-symbol gate before isometry deduplication or equation-cost ranking.

## Consequence for the classification route

The operational order should be reversed from the initial cheap-track
proposal:

1. Pin the rank-seven auxiliary lattice with discriminant form opposite to the
   H3 Neron--Severi lattice.
2. Enumerate primitive embeddings into all 24 Niemeier lattices modulo Weyl
   and glue automorphisms.
3. Compute primitive-closure complements and immediately reject roots,
   determinant failures, and local-genus failures.
4. Deduplicate the survivors by integral frame isometry, with the published
   R17 and alternate Q80 frames as mandatory controls.
5. Retain the full pinned `U` marking and Niemeier embedding provenance for
   every survivor.
6. Only then rank survivors by actual equation cost and the rank-29
   fingerprint.

This Niemeier-first route uses the geometric filter before entering a genus
with at least 167,967 classes.  It also supplies the requested provenance as
part of enumeration rather than as a later reconstruction problem.

No third target-genus rootless frame is certified by the present corpus, and
no theorem that there are exactly two is claimed.  The alternate Q80 frame
therefore remains the first construction priority after the published R17
exclusion: it is the only non-published target-genus control with a pinned
degree-two final neighbour and an existing `h0=5` / three-transform compiler
plan.  This completeness track does not block its characteristic-zero lift.

## Replay

The exact audit, including LLL normalization before `qfisom`, local genus
symbols, mass, automorphism orders, and all source hashes, is
[`elkies-k3-rootless-j2-completeness-track.json`](../artifacts/generated-results/elkies-k3-rootless-j2-completeness-track.json).

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/audit_rootless_j2_completeness_track.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/audit_rootless_j2_completeness_track.sage --check
```

The first command writes the deterministic generated artifact; the second
checks it byte-for-byte.  Neither command enumerates the full genus or proves
a complete Kneser--Nishiyama classification.

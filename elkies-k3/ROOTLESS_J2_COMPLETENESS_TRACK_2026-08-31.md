# Rootless J2 completeness track — 2026-08-31

<!-- status-consumer: EC-K3-H3-ROOTLESS-J2-COMPLETE c6f054948b04b507 -->

## Outcome

The Niemeier-first rootless-frame classification is now **complete at the
`J2` (frame-isometry) level**.  There are exactly two rootless rank-17 frame
classes of determinant 948 for the pinned Picard-rank-19 lattice:

The two mandatory positive controls pass:

| frame | determinant | minimum | norm-four pairs | automorphism order |
| --- | ---: | ---: | ---: | ---: |
| published R17 | 948 | 4 | 1,311 | 2 |
| alternate Q80-derived | 948 | 4 | 1,313 | 4 |

They have the same local genus and PARI `qfisom` proves that they are not
integrally isometric.  Exhaustive primitive-embedding enumeration in the 24
Niemeier lattices proves that these are the only two rootless `J2` frame
classes.  This does not identify `J1` orbits under the surface automorphism
group and does not construct the still-missing characteristic-zero equation
for the alternate frame.

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

## Niemeier-first classification route

The successful operational order reverses the initial cheap-track proposal:

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

## Complete enumeration

The exact rooted-Niemeier catalogue contains 23 hash-pinned unimodular Gram
matrices.  Exact root decompositions exclude ten of them because they contain
no `D5` subsystem; the Leech lattice is excluded because the auxiliary itself
has roots.  The remaining thirteen rooted Niemeier classes contain sixteen
`D5` anchor orbits under their full automorphism groups.  Direct exceptional
root-system enumeration gives one internal `D5` Weyl orbit in each of `E6`,
`E7`, and `E8` (15, 120, and 1,260 fixed-root subsystems before stabilizer
reduction).

After a `D5` anchor is fixed, the auxiliary has a standardized sixth vector
of norm 12 and a seventh vector of norm 24.  An explicit determinant-one
isometry checks that this standardized Gram is the pinned auxiliary Gram.
Residual-Weyl dominance and exact fixed-line enumeration give 3,220 primitive
sixth-vector representatives across all sixteen anchors.

For a fixed sixth vector, move the seventh vector into the residual Weyl
chamber.  Its complement is rootless exactly when every residual Dynkin label
is strictly positive: a zero label makes the corresponding simple root
orthogonal, while positive labels pair positively with every positive root.
The forced `D5` and sixth-vector pairings consume norm `17/4`, leaving the
exact projected norm budget `79/4`.  This gate reduces 3,220 sixth vectors to
21 and leaves only 167 positive-label tuples.  Exact rational LDL enumeration
of the remaining fixed-space ellipsoids gives 25 solutions, of which 12 are
integral ambient vectors.  All twelve full auxiliary embeddings are primitive
and rootless.

Exact integral-isometry deduplication of their saturated rank-17 complements
gives precisely:

| class | cover representatives | Niemeier provenance | norm-four pairs | automorphism order |
| --- | ---: | --- | ---: | ---: |
| published R17 | 8 | `N(2A7+2D5)` | 1,311 | 2 |
| alternate Q80-derived | 4 | `N(2A7+2D5)` | 1,313 | 4 |

The cover counts are not claimed as full-automorphism embedding-orbit counts:
the enumeration deliberately retains possible duplicates under anchor
stabilizers.  They are enough for completeness of the complement isometry
classes.  No other Niemeier ambient produces a rootless complement.  Since
the two survivors are exactly the already pinned controls, their existing
primitive `U` markings and equation-cost data attach without a new marking
ambiguity.  The alternate Q80 frame remains the only open equation lift.

## Rank-seven auxiliary and positive-control lift

The required auxiliary lattice is now pinned at
[`data/lattice/rootless_j2_auxiliary_rank7_gram.txt`](data/lattice/rootless_j2_auxiliary_rank7_gram.txt).
It is even, positive definite, rank seven, and has determinant 948.  Its
cyclic discriminant module has quadratic Gram

```text
[1267/948] on Z/948.
```

This is anti-isometric to the discriminant form of each mandatory rank-17
control.  For the published generator `[773/948]`, multiplication by `43`
is one explicit anti-isometry.  Among the four even positive-definite
rank-seven genera of determinant 948, exactly one has this anti-discriminant
form.  Its exact mass is `80119/20736`, its automorphism group has order
7,680, and its 40 roots span `D5`.

This correction is not recoverable from the old seed labels: all 228 retained
rank-seven determinant-948 Gram matrices in the seed files have the wrong
discriminant form.  None can serve as the auxiliary for the target R17 genus.

Gluing the pinned auxiliary to either control is completely explicit because
both discriminant groups are cyclic.  There are eight anti-isometry units for
each control.  Every one produces an even unimodular rank-24 lattice whose
192 roots split as

```text
(rank, signed roots) = (7,56), (7,56), (5,40), (5,40),
```

so the ambient Niemeier lattice is `N(2A7+2D5)`.  In every case the auxiliary
embedding is primitive and its saturated orthogonal complement is rootless,
has determinant 948, passes the exact target-genus gate, and is integrally
isometric to the requested control.  Double-coset reduction by
`Aut(K) x Aut(R)` gives two primitive embedding orbits for the published R17
control and one for the alternate Q80 control.

This certifies the two required Niemeier positive controls independently of
the complete orbit-reduced enumeration above.  A direct vector enumeration
is excluded by an exact cost probe: after fixing the
control `D5`, its rank-19 determinant-4 orthogonal lattice contains
329,206,692 sign-pairs of norm 12.  The complete traversal therefore has to
work with Weyl/stabilizer orbits or equivalent glue-code combinatorics.

The complete enumeration proves that there is no third target-genus rootless
frame.  The alternate Q80 frame therefore remains the first construction
priority after the published R17 exclusion: it is the unique non-published
rootless `J2` class, with a pinned degree-two final neighbour and an existing
`h0=5` / three-transform compiler plan.  This completeness theorem does not
block its characteristic-zero lift.

## Replay

The exact audit, including LLL normalization before `qfisom`, local genus
symbols, mass, automorphism orders, and all source hashes, is
[`elkies-k3-rootless-j2-completeness-track.json`](../artifacts/generated-results/elkies-k3-rootless-j2-completeness-track.json).

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/audit_rootless_j2_completeness_track.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/audit_rootless_j2_completeness_track.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_rootless_j2_niemeier_controls.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/build_rooted_niemeier_catalog.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/enumerate_niemeier_d5_anchor_orbits.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/enumerate_niemeier_auxiliary_sixth_dominant.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/classify_rootless_j2_niemeier_first.sage --check
```

The generation form of each command omits `--check`; check mode recomputes and
compares its deterministic artifact byte-for-byte.  The final artifact is
[`elkies-k3-rootless-j2-niemeier-first.json`](../artifacts/generated-results/elkies-k3-rootless-j2-niemeier-first.json).
Together the catalogue, anchor, sixth-vector, and final seventh-vector passes
prove the complete rootless Kneser--Nishiyama `J2` classification.  They do
not quotient primitive embeddings to `J1` surface-automorphism orbits.

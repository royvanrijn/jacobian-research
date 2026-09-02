# The determinant-78 frame genus: complete Niemeier classification

<!-- status-consumer: EC-K3-E6-RANK4-DET78-GLOBAL-ROOTFUL 648ec884ce7152bb -->

## Result

Let `F_src` be the positive-definite rank-17 frame in
[`data/lattice/e6_rank4_det78_frame.txt`](data/lattice/e6_rank4_det78_frame.txt).
It is the frame of the explicit `A1+2E6/MW4` fibration constructed in
[`E6_RANK4_LINEAR_CHORD_INCIDENCE_2026-09-02.md`](E6_RANK4_LINEAR_CHORD_INCIDENCE_2026-09-02.md).

The exact all-Niemeier computation in
[`scripts/classify_e6_rank4_det78_niemeier_frames.sage`](scripts/classify_e6_rank4_det78_niemeier_frames.sage)
proves the following J2/frame-level statement.

> **Determinant-78 classification.** The genus of `F_src` contains exactly
> 1,549 integral-isometry classes. Every class contains roots. Their root-rank
> distribution is
>
> | frame root rank | MW rank for rho = 19 | classes |
> |---:|---:|---:|
> | 10 | 7 | 1 |
> | 11 | 6 | 45 |
> | 12 | 5 | 249 |
> | 13 | 4 | 543 |
> | 14 | 3 | 477 |
> | 15 | 2 | 200 |
> | 16 | 1 | 33 |
> | 17 | 0 | 1 |

In particular, this exact determinant-78 Neron--Severi lattice admits **no
rootless MW17 frame**. The proposed corridor

```text
explicit MW4 source -> same NS -> rootless MW17
```

does not exist, so the request to optimize a corridor back from such an
endpoint is vacuous. No neighbour walk was run. The largest Mordell--Weil rank
available in this frame genus is instead 7: the unique class 1058 has root
type `A4+A6`, representative provenance `2A7+2D5`, and frame Gram SHA-256
`e9a8cd27517a7c87e09ad81302d2932ecba2e4df4ce9c10eb8ea9d9135cd0414`.

This class count is J2, namely frame isometry. It is not a J1 quotient by
surface automorphisms, and it neither supplies equations for all classes nor
asserts that the 1,549 rows are distinct surface-automorphism orbits.

## Auxiliary lattice and completeness reduction

The positive rank-seven Nishiyama auxiliary used by the computation is

```text
[ 2 -1  0  0  0  0 -1]
[-1  2 -1  0  0  0  0]
[ 0 -1  2  0  0  0  0]
[ 0  0  0  2 -1  0  0]
[ 0  0  0 -1  2  0  0]
[ 0  0  0  0  0  2  0]
[-1  0  0  0  0  0  4]
```

It has determinant 78 and root system `A3+A2+A1`. The last generator has norm
4, pairs with an **endpoint** of the displayed `A3`, and has residual Schur
norm `13/4`. Its Gram SHA-256 is
`470802fbec420262d1d020160ef4a935eb42837f2d603450523dd755543558de`.
The exact finite quadratic-form check identifies its discriminant form with
the negative of the source-frame form. Thus every frame class in the target
genus can be glued primitively to this auxiliary inside a positive even
unimodular rank-24 lattice.

The completeness argument implemented in the classifier is:

1. The auxiliary has roots, so the ambient unimodular lattice cannot be the
   Leech lattice. The computation processes all 23 rooted Niemeier lattices in
   the pinned full-glue Gram catalogue.
2. It enumerates `A3` root-sublattice orbits, then `A2` and `A1` in successive
   orthogonal root systems. Stabilizers are exact residual Weyl groups; keeping
   extra representatives can only duplicate the cover.
3. For each primitive root anchor it moves the last vector to the residual
   dominant chamber, enumerates every Dynkin-label tuple within the exact
   `13/4` budget, and solves the remaining shifted ellipsoid exactly.
4. Every integral primitive auxiliary embedding yields its saturated
   rank-17 orthogonal complement. Complements are LLL-rebased only by recorded
   unimodular transformations.
5. Frames are deduplicated by an exact root-span/glue isometry test. After
   quotienting by root reflections, it tests Dynkin graph maps and the
   low-rank root-orthogonal isometry group; integrality is an exact finite
   congruence. Rootless inputs would fall back to full PARI integral isometry.

The cover contains 1,591 primitive root anchors, 737,139 dominant label
tuples, 231,160 exact ellipsoid solutions, and 37,397 integral primitive
auxiliary embeddings. These deduplicate to the 1,549 classes above.

## Independent mass closure

For each deduplicated frame `F`, the computation obtains the exact order of
`Aut(F)` as

```text
|W(Root(F))| * number of glue-extending (Dynkin outer map, quotient automorphism) pairs.
```

Every accepted extension is checked as an integral unimodular Gram isometry.
The resulting weighted class sum is

```text
sum_F 1/|Aut(F)| = 1463420154787 / 4131952105881600.
```

Sage's exact genus mass for `Genus(F_src)` is independently

```text
1463420154787 / 4131952105881600.
```

The equality is asserted before the full artifact is written. Together with
the exhaustive primitive-embedding cover and exact isometry deduplication,
this closes the genus and rules out an omitted rootless class.

The explicit source is recovered exactly once, as class 222. It has root type
`A1+2E6`, MW rank 4, representative provenance `D10+2E7`, and automorphism
group order 42,998,169,600. This is the mandatory positive control connecting
the classification to the explicit source family.

## Reproduction and pinned outputs

Run the complete classification with Sage 10.9:

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/classify_e6_rank4_det78_niemeier_frames.sage
```

Then replay it byte-for-byte:

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/classify_e6_rank4_det78_niemeier_frames.sage --check
```

The complete run took about 20 minutes and 444 MB peak resident memory on the
2026-09-03 workstation. Its machine-readable output is
[`../artifacts/generated-results/elkies-k3-e6-rank4-det78-niemeier-frames-v1.json`](../artifacts/generated-results/elkies-k3-e6-rank4-det78-niemeier-frames-v1.json),
SHA-256
`be549ce3ca1a74fd3c4df5133521b54f67612a75b10ea635c5671ba5881f3fc6`.
It stores every class Gram, root data, automorphism order, representative
primitive embedding, ambient provenance counts, and all global accounting.

Pinned inputs are:

- source frame SHA-256
  `b441ba2d8440ddd226bcb839a065a2a25ffbd99f58e7a1223088479402080f04`;
- rooted Niemeier catalogue SHA-256
  `abbfb31ca37af8fb739fbbd0bd2575a3b63d8d1031ef4f962d5a79e674d2ce53`.

## Literature context

This is the standard Nishiyama/Kneser strategy carried out with full Niemeier
glue and an explicit class-mass closure. The relevant general references are
K. Nishiyama, *The Jacobian fibrations on some K3 surfaces and their
Mordell--Weil groups*, Japan. J. Math. 22 (1996), 293--347,
<https://doi.org/10.4099/math1924.22.293>; M. Schuett and T. Shioda,
*Elliptic surfaces*, <https://arxiv.org/abs/0907.0298>; and A. Braun,
Y. Kimura, and T. Watari, *The Noether--Lefschetz problem and gauge-group
resolved landscapes*, <https://arxiv.org/abs/1312.4421>, especially the
primitive-embedding double quotient underlying J2. The need to retain the full
Niemeier overlattice rather than only its root system is illustrated in
M.-J. Bertin and O. Lecacheux, *Apéry--Fermi pencil of K3-surfaces and their
20 elliptic fibrations*, <https://doi.org/10.5802/pmb.44>. A targeted search
found no published determinant-78 classification matching this exact source
frame.

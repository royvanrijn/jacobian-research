# Ignition-first rank-growth search

> **Historical pre-endpoint workflow.** The rank-17 model and its seventeen
> generic sections have since been recovered exactly.  Current rank-32 work
> must first pass the fail-closed residual quotient
> `Sel_2(E_t)/<P1,...,P17>` on the same parameter and minimal curve.  This
> height-Gram ignition/cascade interface remains useful only after that gate or
> for replaying historical controls; it is not authorization for a raw point
> search.  See [`README.md`](README.md) and
> [`../elliptic-curves/README.md`](../elliptic-curves/README.md).

This workflow implements the search strategy learned from the controlled
rank-21 experiments:

1. find a genuine rank-17 -> rank-18 **ignition** event by minimizing the
   positive Schur-complement height;
2. after ignition, search sequentially for rank 19, 20 and 21;
3. at each cascade step, prefer candidates whose projection onto the added
   rank-growth block is strongly aligned with the **newest orthogonal
   increment**.

The controlled rank-21 experiment motivating this had median newest-increment
correlations approximately

- rank 18 -> 19: `1.000`
- rank 19 -> 20: `0.958`
- rank 20 -> 21: `0.946`

The invariant ignition quantity for a candidate point `P` above a current
height Gram `G` is

```text
delta(P) = h(P) - b^T G^-1 b,
```

where `b_i = <P, P_i>`.  A positive `delta` is the squared height of the
component transverse to the current Mordell-Weil span.

## Superseded reconstruction limitation

At the time of this workflow, `data/k3-model/README.md` still listed the
canonical Elkies rank-17 Weierstrass model and generic Mordell--Weil sections
as artifacts to recover.  They are now pinned in
[`data/fibrations/elkies_2026_published_r17_model.json`](data/fibrations/elkies_2026_published_r17_model.json)
and the endpoint certificates described in [`README.md`](README.md).
The code still consumes height-Gram records rather than inventing a
specialization or point-finding formula, but present use is subordinate to the
residual 2-Selmer gate above.

### Historical reconstruction checkpoint

The model-recovery work is now active rather than purely archival.  See
[`RECONSTRUCTION_PROGRESS.md`](RECONSTRUCTION_PROGRESS.md) for the current
checkpoint. As of the consolidated 2026-08-20 checkpoint, the repository has
an exact Coxeter-9 reduction of the recovered rank-17 lattice and an optimal
determinant-2 extension chain back to the full saturation. The old numerical
rank-10 hit on `root-000029` has been downgraded to a near-collision. The
primary route now passes through an E6/MW3 neighbor: one smooth `GF(31)`
one-section seed is exact, and seven complete split-root cores remain for the
reduced canonical `P2` search. The characteristic-zero rank-17 model is still
pending.

## Candidate TSV contract

Both search drivers consume a tab-separated file with at least:

```text
candidate_id    gram
```

Optional columns are:

```text
parameter    points    new_point    source
```

Relative file paths are resolved relative to the TSV file.

For ignition, each Gram matrix must contain the known 17 specialized basis
vectors first, followed by one or more point candidates.

For a cascade step, each Gram matrix must contain the *current canonical
basis* first (`17 + all previously accepted new directions`), followed by one
or more point candidates.  The driver checks this leading Gram block against
the parent hit and rejects mismatches.

A manifest can be generated from result directories with:

```bash
python elkies-k3/scripts/build_rank_growth_candidate_manifest.py \
  --gram-glob 'elkies-k3/results/specializations/*/height-gram.txt' \
  --out elkies-k3/results/rank-growth-candidates.tsv
```

## Self-test

```bash
python elkies-k3/scripts/test_rank_growth_geometry.py
```

Expected output begins with `PASS`.

## Rank-18 ignition

```bash
python elkies-k3/scripts/search_rank18_ignition.py \
  --candidates elkies-k3/results/rank-growth-candidates.tsv \
  --out elkies-k3/results/rank18-ignition
```

Every retained hit gets a directory containing:

```text
hit-000001/
    height-gram.txt
    pairing-vector.txt
    projection-coefficients.txt
    parameter.txt
    points.txt          # when supplied
    new-point.txt       # when supplied
    ignition.json
```

`hits.tsv` is ordered by increasing positive transverse height.

## Rank-19 cascade

Prepare candidate Grams whose first 18 rows are the canonical basis from the
chosen ignition hit, then run:

```bash
python elkies-k3/scripts/search_rank_cascade.py \
  --current-hit elkies-k3/results/rank18-ignition/hit-000001 \
  --candidates elkies-k3/results/rank19-candidates.tsv \
  --out elkies-k3/results/rank19-cascade
```

The default `--min-alignment 0.85` is intentionally below the ~0.95 median
observed in the controlled experiment.  `candidates.tsv` retains all scored
independent candidates; `hits.tsv` contains candidates passing the alignment
prior.

Each cascade hit writes `current-height-gram.txt`, the canonical basis for the
next step.

## Rank 20 and rank 21

The same driver is recursive:

```bash
python elkies-k3/scripts/search_rank_cascade.py \
  --current-hit elkies-k3/results/rank19-cascade/hit-000001 \
  --candidates elkies-k3/results/rank20-candidates.tsv \
  --out elkies-k3/results/rank20-cascade

python elkies-k3/scripts/search_rank_cascade.py \
  --current-hit elkies-k3/results/rank20-cascade/hit-000001 \
  --candidates elkies-k3/results/rank21-candidates.tsv \
  --out elkies-k3/results/rank21-cascade
```

The drivers are numerical discovery tools.  Any claimed rank increase still
needs the repository's exact independence/saturation/certification path before
being treated as a proof.

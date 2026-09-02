# Fermigier rank-20 direct direction search

This experiment avoids general 2-descent on the Fermigier--Mestre anchor
`u=28917/20` (`T=28917/10`). It searches the original quartic and alternate
exact degree-two charts of the same elliptic curve, then filters returned
points against the pinned rank-20 subgroup.

The runner never turns bounded point-search absence into a rank upper bound.
Rank is promoted only when an augmented finite-quotient certificate is replayed
successfully.

## First run: smoke stage

```bash
mkdir -p artifacts/local/elliptic-curves

PYTHONUNBUFFERED=1 caffeinate -i \
  sage -python elliptic-curves/cas/search_fermigier_rank20_direction.py \
  --stage smoke \
  2>&1 | tee artifacts/local/elliptic-curves/fermigier_rank20_direction_smoke.log
```

The smoke stage uses a few low-weight alternate covers and complementary skew
boxes. The default timeout is 60 seconds per PARI call; timed-out jobs are
recorded rather than retried.

## Full resumable cover frontier

```bash
PYTHONUNBUFFERED=1 caffeinate -i \
  sage -python elliptic-curves/cas/search_fermigier_rank20_direction.py \
  --stage frontier \
  2>&1 | tee artifacts/local/elliptic-curves/fermigier_rank20_frontier.log
```

This scans all `2^20-1` represented nonzero classes in Gray-code order. The
atomic plan checkpoint is written after every 65,536 classes. Interrupting and
rerunning the same command resumes from the last complete checkpoint.

Then search the selected covers:

```bash
PYTHONUNBUFFERED=1 caffeinate -i \
  sage -python elliptic-curves/cas/search_fermigier_rank20_direction.py \
  --stage search \
  2>&1 | tee artifacts/local/elliptic-curves/fermigier_rank20_direction_full.log
```

## Sharding

Build the frontier once, copy its plan file to each machine, then run disjoint
search shards. For example, with four shards:

```bash
sage -python elliptic-curves/cas/search_fermigier_rank20_direction.py \
  --stage search --shard-count 4 --shard-index 0
```

Use indices `0`, `1`, `2`, and `3`. With the default output path, each shard
writes a distinct `.shard_NNN_of_NNN.json` artifact.

## Result interpretation

- A height-pairing relation counts only after exact `Fraction` group-law replay.
- Overcomplete mod-2/mod-3/mod-5 quotient-bank escape is exact prioritization
  evidence, but not yet a promoted Mordell--Weil rank claim.
- `certified_rank_lower_bound_after_search: 21` or `22` is emitted only after
  the augmented point set receives and replays a complete finite-quotient
  infinite-descent certificate.
- A result remaining at 20 says only that this bounded search did not certify a
  new direction.

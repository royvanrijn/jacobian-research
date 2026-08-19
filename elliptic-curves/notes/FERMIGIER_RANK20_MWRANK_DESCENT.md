# Fermigier rank-20 eclib/mwrank descent

This is the second Magma-free backend for the Fermigier rank-20 relative-rank
experiment. It avoids PARI's `ellrankinit` cubic class-group initialization and
runs Cremona's independent eclib/mwrank 2-descent implementation.

The default mode is **Selmer-only**. mwrank does not accept the repository's
known rank-20 basis as descent generators, so the final interval combines:

1. the repository's exact finite-quotient certificate proving `rank >= 20`;
2. mwrank's rigorous 2-descent upper bound.

Run from the repository root:

```bash
PYTHONUNBUFFERED=1 caffeinate -i \
  sage -python elliptic-curves/cas/run_fermigier_rank20_mwrank_descent.py \
  2>&1 | tee artifacts/local/elliptic-curves/fermigier_rank20_mwrank_descent.log
```

The program emits `R20MWRANK|...` records immediately before and after the
expensive descent and writes a deterministic JSON summary under
`artifacts/local/elliptic-curves/`.

Classification:

- **M0**: combined interval `[20,20]`; exact rank 20 follows from the independent
  lower certificate plus mwrank's upper bound.
- **M2**: combined interval `[20,R]`, `R>20`; residual 2-Selmer room remains.
- **M3**: mwrank itself finds a lower bound at least 21. In Selmer-only mode this
  is unlikely; rerun with `--search-points` only if explicit point discovery is
  desired.
- an upper bound below 20 is a contradiction/failure, not a mathematical result.

Useful options:

```bash
# Quiet machine-readable mode
sage -python elliptic-curves/cas/run_fermigier_rank20_mwrank_descent.py --quiet

# Permit quartic point searches after the Selmer computation
sage -python elliptic-curves/cas/run_fermigier_rank20_mwrank_descent.py \
  --search-points --n-aux 22
```

For this curve the rational 2-torsion rank is expected to be zero. In that case
mwrank's rank upper bound equals its computed 2-Selmer rank. The script records
both values rather than assuming equality silently.

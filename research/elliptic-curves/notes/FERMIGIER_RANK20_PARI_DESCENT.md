# Fermigier rank-20 Sage/PARI descent

This is the Magma-free first-line continuation of the relative-descent experiment.
It uses the same corrected bounded-saturation rank-20 basis and exact mod-2
finite-quotient certificate as the Magma job.

Run from the repository root:

```bash
sage -python elliptic-curves/cas/run_fermigier_rank20_pari_descent.py \
  2>&1 | tee artifacts/local/elliptic-curves/fermigier_rank20_pari_descent.log
```

The default effort sequence is `0,1,2`. A more expensive continuation is:

```bash
sage -python elliptic-curves/cas/run_fermigier_rank20_pari_descent.py \
  --efforts 0,1,2,4,8 \
  2>&1 | tee artifacts/local/elliptic-curves/fermigier_rank20_pari_descent.log
```

PARI `ellrank(E, effort, points)` accepts the twenty known rational points
explicitly. It performs 2-descent and Cassels-pairing restrictions and returns
`[r1,r2,s,L]`, where `r1 <= rank(E) <= r2` and `s` is the obstruction rank
that PARI detects in `Sha[2]/2Sha[4]`. The repository's exact independent
20-point certificate is combined with this to replace a weak `r1` by 20.
For this cubic 2-division field, PARI constructs provisional BNF data under
GRH.  Its source does not call `bnfcertify`, so an upper endpoint is
GRH-conditional unless a separate BNF certification is supplied.  Rational
points returned in `L` can still yield unconditional lower bounds after exact
independent verification.

Classification:

- **P0**: combined interval `[20,20]`; this conditionally closes the curve
  under GRH, but is not an unconditional exact-rank certificate.
- **P2**: interval `[20,R]` with `R>20`; the conditional 2-descent still leaves
  residual room.
- **P3**: PARI itself returns lower bound at least 21; this is a genuine new-rank
  signal and the returned points should immediately be extracted and independently
  checked against the rank-20 subgroup.
- errors, timeouts, or an upper bound below 20 are failures/contradictions, not
  mathematical results.

Unlike the Magma experiment this does not expose individual 2- and 4-covers.
Its advantage is that it is free, accepts the known points directly, and its
provisional upper bound incorporates the Cassels-pairing restriction. If P2
remains after serious effort, explicit cover-level work is still useful.

## eclib fallback

If `ellrankinit` remains inside PARI's cubic class-group computation for several
hours, use the independent eclib/mwrank runner documented in
`FERMIGIER_RANK20_MWRANK_DESCENT.md`. It performs a different 2-descent and does
not invoke PARI's `ellrankinit`.

# Newfamily specialization `T = 83/6`: exact rank 14

## Status

For the six-root quartic family with

```text
(-47,-43,-31,30,45,46)
```

and specialization

\[
T=\frac{83}{6},
\]

the eleven recovered generic hidden sections specialize to rank 11. Three additional rational points enlarge the subgroup to rank 14, and PARI `ellrank` at effort 0 returns the unconditional rank interval `[14,14]`.

Therefore

\[
\boxed{\operatorname{rank} E_{83/6}(\mathbf Q)=14.}
\]

The exact-rank certificate is

```text
artifacts/generated-results/elliptic-curves/newfamily_rank14_t83_6_pari_exact_rank_v1.json
```

The earlier lower-bound certificate remains

```text
artifacts/generated-results/elliptic-curves/newfamily_rank14_t83_6_v1.json
```

and the git-only replay driver is

```text
elliptic-curves/cas/newfamily/certify_rank_t83_6.py
```

## Lower-bound replay

The baseline-first eclib verification gives

```text
known hidden sections : rank 11
Q2                    : 11 -> 12
Q3                    : 12 -> 13
Q4                    : 13 -> 14
```

On the global minimal model the three rank-increasing points are

\[
Q_2=(-16951807,\ 1852437557223),
\]

\[
Q_3=(-13752559,\ 1821538211319),
\]

\[
Q_4=(-4217857,\ 1724943248823).
\]

Their high-precision Schur residuals during discovery were approximately

```text
Q2  0.43971405055972258...
Q3  0.23146760463960279...
Q4  0.23425871898606799...
```

The Schur test was only triage; the rank increases above were checked by exact eclib group processing.

## Exact upper bound

The committed replay then feeds the resulting 14-point subgroup to PARI `ellrank`:

```text
T83RANK|stage=eclib|status=complete|baseline=11|lower=14
T83RANK|stage=ellrankinit|status=complete|seconds=3.797025
T83RANK|stage=ellrank|status=complete|effort=0|pari_lower=14|pari_upper=14|effective_lower=14|sha_pairing_rank=0|returned_points=14|seconds=1.274649
T83RANK|stage=done|interval=14,14|classification=exact_rank_14
```

Thus the exact eclib lower bound and the PARI 2-descent/Cassels-pairing upper bound coincide at 14.

## Discovery context

This specialization was found after changing the search objective from large rational Nagao winners to low-projective-height rational parameters where the specialized known sections remain accessible to ordinary point search.

The initial unseeded search produced thirteen distinct Schur-triage hits. Exact processing showed that only three contributed new independent Mordell-Weil directions; the remaining ten lie in the resulting rank-14 subgroup.

The discovery-run Nagao scores were

```text
discovery = 5.132090535805
held      = 3.650137875399
```

The global minimal model has root number `+1` and a 169-bit absolute discriminant.

## Current follow-up

A deeper H16/H18/H20 search is being used as a control and neighborhood probe. Since the exact rank is now proved to be 14, no further rational point on this same curve can raise its rank: any new points found at larger search height must lie in the established rank-14 Mordell-Weil group.

The more promising next rank-search targets are the six specializations in `newfamily_exact_subgroup_rank_gain_batch_v1.json` with certified subgroup rank at least 13, especially those with small point sizes and many unused Schur hits.

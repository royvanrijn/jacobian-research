# Newfamily specialization `T = 83/6`: exact subgroup rank at least 14

## Status

This note records an exact lower-bound result for the six-root quartic family with

```text
(-47,-43,-31,30,45,46).
```

At the rational specialization

\[
T=\frac{83}{6},
\]

the eleven recovered generic hidden sections specialize to an exact rank-11 subgroup. Three additional rational points found by the low-projective-height unseeded search independently enlarge that subgroup to rank 14.

Therefore

\[
\boxed{\operatorname{rank} E_{83/6}(\mathbf Q)\ge 14.}
\]

This is a lower bound only. No upper bound or full saturation statement is claimed here.

The pinned compact certificate is

```text
artifacts/generated-results/elliptic-curves/newfamily_rank14_t83_6_v1.json
```

and the exact batch verifier is

```text
elliptic-curves/cas/newfamily/batch_verify_v2_rank_gain_hits.py
```

## Discovery path

The result came from a rational-parameter search that was deliberately changed from a pure Nagao-score objective to a low-projective-height search. Large rational parameters had strong local scores but specialized known-section naive heights typically around 35--50 or higher, making ordinary eclib point search ineffective. Restricting projective parameter height exposed a much more searchable population.

For `T=83/6`, the unseeded search produced thirteen distinct Schur-triage hits. Exact eclib processing showed that only three of those add independent Mordell--Weil directions; the remaining ten lie in the resulting rank-14 subgroup.

The Nagao scores attached to this specialization in the discovery run were

```text
discovery = 5.132090535805
held      = 3.650137875399
```

The root number computed on the global minimal model is `+1`.

## Exact baseline

The verifier processes the known specialized sections first, sorted only within the known-section phase by coordinate bit size. Their exact eclib rank growth is

```text
U1   0 -> 1
U2   1 -> 2
U5   2 -> 3
U7   3 -> 4
U3   4 -> 5
U4   5 -> 6
U6   6 -> 7
U9   7 -> 8
U10  8 -> 9
U8   9 -> 10
U0  10 -> 11
```

Thus the specialized known subgroup is exactly rank 11 as processed by eclib.

## Three new independent points

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

Their exact rank growth is

```text
known subgroup : rank 11
Q2             : 11 -> 12
Q3             : 12 -> 13
Q4             : 13 -> 14
```

The corresponding points on the fixed homogeneous short model used during discovery are

```text
Q2 = (-2722114156253184, 3769442154946420392591360)
Q3 = (-2208380243300352, 3706569428268172097617920)
Q4 = ( -677301762576384, 3510021470235349179432960)
```

The high-precision relative Schur residuals that triggered exact verification were approximately

```text
Q2  0.43971405055972258...
Q3  0.23146760463960279...
Q4  0.23425871898606799...
```

For comparison, the `T=11` rank-11 calibration produced dependent-point residual noise around `4e-76`, so these were intentionally treated only as triage signals until exact eclib processing confirmed the three independent directions.

## Computational profile

The global minimal model has a 169-bit absolute discriminant. In the exact verification run, point-coordinate sizes ranged from 39 to 216 bits, with median 41 bits. The full baseline-first verification completed in about 3.1 seconds on the development machine, including about 0.014 seconds for global minimalization.

This makes `T=83/6` a particularly convenient specialization for deeper work: the three new generators are small, the model is inexpensive to rebuild, and exact subgroup processing is fast.

## Next steps

The immediate follow-up should be reproducible and committed:

1. attempt PARI/Sage rank bounds using the 14-point subgroup;
2. run controlled saturation diagnostics on the rank-14 subgroup;
3. search for a 15th point at increasing eclib heights on this small specialization;
4. preserve any successful upper-bound or further rank-gain result as a new versioned certificate under `artifacts/generated-results/elliptic-curves/`;
5. finish the corrected baseline-first batch over the other low-height specializations and promote that compact batch result separately.

Do not infer `rank = 14` from this note. The theorem-strength statement currently recorded is only

\[
\operatorname{rank}E_{83/6}(\mathbf Q)\ge14.
\]

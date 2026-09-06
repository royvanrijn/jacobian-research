# The six rank-27 search subgroups are saturated at 2, 3 and 5

The [five proper endpoint section indices](ENDPOINT_SECTION_SPANS_2026-09-06.md)
do not explain the null specialized-parity experiment on the six new rank-27
curves. The [exact certificate](../../artifacts/generated-results/elliptic-curves/new27_small_prime_saturation_v1.json)
proves that each particular 27-point subgroup used in that experiment is
saturated at 2, 3 and 5. Full saturation and whole-curve rank remain unknown.

The six V12 inventory IDs are 40, 41, 48, 71, 72 and 90. Their exact equations,
27-point generators, finite quotient matrices and no-torsion witnesses are
included in the certificate. The checker recomputes the finite reductions on
the 27 seed points themselves; it does not infer their primitivity merely from
the rank of a larger returned point cloud. No new point search is involved.

## Why the finite matrices prove this particular saturation claim

Let `H = <P_1,...,P_r>` in `E(Q)`. Suppose the images of these generators are
independent in a product of finite quotients `E(F_p)/ell E(F_p)`, and suppose
`E(Q)[ell] = 0`. If a rational point `Q` satisfies

```
ell Q = a_1 P_1 + ... + a_r P_r,
```

reduction in those quotients forces every `a_i` to be divisible by `ell`.
Thus `ell(Q - sum (a_i/ell) P_i) = 0`; the no-torsion witness gives `Q in H`.
This proves `ell`-saturation of `H` in `E(Q)`, even if the whole curve has more
independent directions than `H`.

All six seed subgroups have the required independent columns and no-torsion
witnesses for `ell = 2, 3, 5`. The odd-prime checks use the same retained prime
sets as the preceding cloud certificates, restricted to their first 27 columns
and independently recomputed. Consequently no index-three or index-five
correction of the type found at the endpoints can enlarge these six subgroups.

This is a limited consequence of exact finite data, not a whole-curve rank
upper bound. It leaves divisions at larger primes, full saturation, independent
directions beyond 27 and the effectiveness of different chart policies open.
The existing inventory remains 101 curves, including six with certified lower
bound 27.

Both build and exact replay pass:

```sh
python3 elliptic-curves/cas/certify_new27_small_prime_saturation.py --check
```

The source and certificate are retained for the next evidence supplement;
this note does not claim an isolated portable replay for this additional result.

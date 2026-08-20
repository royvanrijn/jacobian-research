# ICARM curve 273: construction and family investigation

Status: **ongoing source audit and bounded computation**.  This note does not
identify a specialization parameter for curve 273 and does not prove that it
belongs to the Elkies--Klagsbrun rank-17 K3 family.

## Bottom line

The curve itself and 30 independent points are reproducible: see
[`ICARM_CURVE273_RANK30.md`](ICARM_CURVE273_RANK30.md).  The construction is not
yet reproducible.  The current public ICARM entry attributes the find to
Claude, with Levent Alpöge and Ava Howell, but neither that entry nor the public
leaderboard discussion gives a family equation, search parameter, or
specialization certificate.

The strongest working hypothesis is that curve 273 came from the same broad
rank-17 K3 specialization programme as the public rank-28 and rank-29 records.
This is plausible and useful, but it remains a hypothesis.  A proof requires
one of the following:

1. a discoverer-supplied construction record;
2. an explicit rank-17 family and a rational parameter whose specialization is
   isomorphic over `Q` to curve 273; or
3. an equivalent exact specialization certificate, including the transport of
   the generic sections.

## Public provenance recovered on 2026-08-20

- The [ICARM curve 273 page](https://elliptic-rank.icarm.cloud/curve/273)
  records the equation, 30 points, the submitter name `ranksunbounded`, and the
  attribution "Claude, with Levent Alpöge and Ava Howell".
- The [public ICARM Zulip topic](https://icarm.zulipchat.com/#narrow/channel/519875-general/topic/Elliptic.20Curve.20Rank.20Leaderboard/near/603443505)
  records that the initial 29-point display was a parser omission, after which
  the thirtieth submitted point was restored.  It also records the conditional
  analytic-rank discussion, but no construction recipe.
- Dujella's maintained [rank-30 page](https://web.math.pmf.unizg.hr/~duje/tors/rk30.html)
  reproduces the curve and its 30 points.

The earlier pinned source audit
[`elliptic_rank30_public_source_audit.json`](../../artifacts/generated-results/elliptic_rank30_public_source_audit.json)
predates curve 273.  It remains useful for the rank-29 construction history,
but its conclusion that no public rank-30 curve was available has been
superseded by the 2026-08-20 record.

## What is known about the preceding record family

The public rank-29 announcement says that its search used a rank-17 elliptic
fibration on the same K3 surface used for the rank-28 record, sieved rational
specializations, and searched outside the generic `Z^17`.  The announcement
did not publish the family equation or the 17 generic sections.

The independent reconstruction programme in [`../../elkies-k3/`](../../elkies-k3/)
has identified the recovered rank-17 Mordell--Weil lattice with the Shimura
datum

```text
quaternion discriminant D = 6
level M = 79
rank-17 lattice determinant = 948 = 2^2 * 3 * 79.
```

The most economical current explicit route is the `E6/MW3` neighbour with

```text
ADE = E6 + A3^2 + A1^2
MW rank = 3
fibres = IV* + I4 + I4 + I2 + I2 + 4 I1
reduced MW determinant = 79/16.
```

See [`../../elkies-k3/RECONSTRUCTION_PROGRESS.md`](../../elkies-k3/RECONSTRUCTION_PROGRESS.md)
and [`../../elkies-k3/E6_P2_REDUCTION_2026-08-20.md`](../../elkies-k3/E6_P2_REDUCTION_2026-08-20.md).
Recovering this family remains more valuable than repeating large blind
specialization searches: it would expose the actual parameter geometry and
make searches beyond rank 30 reproducible.

## Structural comparison with the rank-28 and rank-29 curves

Exact invariant calculations give the following small-prime valuations in the
displayed integral discriminants.  These are fingerprints, not family
certificates; minimalization and the distinction between discriminant and
conductor remain essential.

| curve | selected discriminant valuations |
|---|---|
| rank 28 | `2^15 3^6 5^6 7^4 11^2 13^4 17^5` |
| rank 29 | `2^19 3^7 5^7 7^4 11^5 13^3 17^4 31^3 41^2` |
| curve 273 | `2^16 3^12 5^8 7^5 13^5 31^2 41^2 47^4 53^3 67^3` |

Thus all three share the bad-prime pattern `2,3,5,7,13`, and the last two also
share `31,41`.  For the rank-29 and rank-30 displayed models,

```text
gcd(|Delta_29|, |Delta_30|)
  = 95418385098324986880000000
  = 2^16 * 3^7 * 5^7 * 7^4 * 13^3 * 31^2 * 41^2.
```

On the other hand, the three `j`-invariants are pairwise different.  For the
rank-29 and rank-30 models, `gcd(c4_29,c4_30)=1` and
`gcd(c6_29,c6_30)=1`.  In particular curve 273 is not merely the same curve in
a scaled Weierstrass presentation, and the invariant comparison gives no
simple twist explanation.

Interpretation: the shared discriminant support is consistent with a common
CRT-shaped search lineage, but it does not establish membership in `X(6,79)`.

A calibrated bounded run of
[`../../elkies-k3/scripts/search_rank17_embedding_graph_v2.py`](../../elkies-k3/scripts/search_rank17_embedding_graph_v2.py)
used 12,000 ambient short-vector lines and the first eight eligible norm shells.
Both the rank-29 control and curve 273 reached only partial graph depth `2`.
The result is therefore non-discriminating and does not justify a larger
heuristic search.  The local checkpoints are
`artifacts/local/elkies-k3/rank17-E29-control-bounded-20260820-best-partial.txt`
and `artifacts/local/elkies-k3/rank17-E30-bounded-20260820-best-partial.txt`.
The common options were

```text
--limit 12000 --max-shells 8 --max-shell-lines 1000
--node-limit 100000 --seconds-per-shell 4
```

## Exact rank and conditional upper bound

The repository certificate proves unconditionally

```text
rank E(Q) >= 30.
```

It does not prove exact rank 30.  The ICARM discussion reports an analytic
upper bound of 31 under GRH and uses root number `+1` with BSD parity to obtain
conditional exact rank 30.  The hypotheses are indispensable; the repository
therefore retains only the unconditional lower-bound claim.

## Residual 2-descent: exact bounded progress

The relation search began with a six-large-prime target and has now produced
the exact support-size chain

```text
6 -> 4 -> 2.
```

The successful stage-two target was

```text
505724623:356162826
84664160213:21346805921
541738517197:261717997519
28691731813798755604363789:17957201189903465826327159
```

Forcing the first three ideals with `--top 2000` gave 399 exactly factored
candidates (the remaining shortlisted candidates exceeded the 160-bit
factorization guard), 268 improvements, and the two-ideal residual

```text
28691731813798755604363789:17957201189903465826327159
7159638381133483906634203654283170391:12780381281373253031851035100853459
```

The new relation is represented by

```text
m = 14332057143548341066300258343667194241.
```

Raw log:
`artifacts/local/elliptic-curves/crt-cycle-stage2.log`.

Two bounded continuations explain why the same one-dimensional search should
not simply be enlarged:

| forced ideal | shortlist | exactly tested | improvement |
|---|---:|---:|---:|
| old 85-bit ideal | 2000 | 2000 | none |
| new 123-bit ideal | 2000 | 6 | none |

The 123-bit modulus is already larger than the real-root scale.  Almost all
points in the integral arithmetic progression then have residual cofactors
above the factorization guard.  Raising the guard would mostly buy harder
factorizations, not better geometry.

The clean continuation is to leave the restricted elements `m-theta` and work
in the full cubic field.  Construct the product of the two target prime ideals,
reduce its three-dimensional Minkowski lattice (with the bad-prime `S` support
handled explicitly), enumerate short elements, and factor their exact norms.
This searches all three coefficient directions and aligns the enumeration with
algebraic norm.  A naive unweighted two-coordinate `a-b*theta` lattice was
tested and rejected: the large curve coefficients make its norm residuals
worse, not better.

## Efficient next gates

1. **Family certificate first.**  Recover the explicit `E6/MW3` model and its
   rank-jump locus, or obtain the discoverers' construction data.  Then test
   curve 273 by exact specialization and `Q`-isomorphism.
2. **Lattice fingerprint second.**  Once an explicit generic rank-17 height
   matrix is available, find an exact isometric embedding into the curve-273
   Mordell--Weil height lattice and transport the corresponding points.  A
   heuristic short-vector graph match alone is not a certificate.
3. **Full-ideal descent, not wider integral CRT.**  Use reduced ideal lattices
   and `S`-unit normalization for the current two-ideal target.  Keep the
   factorization bound and enumeration radius explicit so every negative run
   remains a bounded experiment.
4. **Only then search for improvements.**  With the family parameter exposed,
   optimize specializations by exact local conditions and finite-quotient
   escape of new points.  Re-minimize every specialization before interpreting
   discriminant valuations or conductor.

These gates distinguish two separate objectives: reconstructing how the curve
was found, and proving whether its algebraic rank is exactly 30.  Progress on
one does not automatically settle the other.

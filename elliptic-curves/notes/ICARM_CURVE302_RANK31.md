# ICARM curve 302: unconditional rank at least 31

Status: **exact unconditional rank lower bound, independently replayed**.
No unconditional exact-rank statement is claimed.

## Public record

On 2026-08-23 David Renshaw pointed to ICARM curve 302:

- <https://x.com/dwrensha/status/2091619906701689228>
- <https://elliptic-rank.icarm.cloud/curve/302>

The ICARM entry attributes the find to Claude, Levent Alpöge, and Ava Howell
and gives the integral model

```text
y^2 + x*y + y = x^3 + x^2 + A*x + B

A = -1284727764113567728281797636015784768866707681415849262157224232063
B = 560368321454261339256859338901915312332769858684945406858043869199456710681989058863306170127006181.
```

The page supplies 31 rational points and reports a conditional BSD+GRH closure
at rank 31. The conditional upper-bound calculation is not reproduced here;
this note and its artifact use only the unconditional lower-bound data.

## Exact independent certificate

The checker is
[`verify_icarm_curve302_rank31.py`](../cas/verify_icarm_curve302_rank31.py).
It performs the following exact replay.

1. All 31 public rational points satisfy the displayed Weierstrass equation.
2. The rational change of variables

   ```text
   X = 36*x + 15,
   Y = 108*(2*y + x + 1)
   ```

   transports the curve and every point to the integral short model

   ```text
   Y^2 = X^3
       - 1665007182291183775853209736276457060451253155114940643755762604753675*X
       + 26144544405770017044368029315807785787305444894561450700514938908225759062376208254436069009884671694150.
   ```

3. The short 2-division cubic has no root modulo 31, so `E(Q)[2]=0`.
4. Independently, exact point counts give

   ```text
   #E(F_17) = 26,
   #E(F_31) = 43.
   ```

   These coprime orders force the complete rational torsion subgroup to be
   trivial.
5. Exact exhaustive finite-group calculations at

   ```text
   17,47,53,61,67,71,79,83,89,101,107,113,127,137,
   149,179,191,197,211,233,241,263,269,281,283,293,311
   ```

   give 32 rows in the product of `E(F_p)/2E(F_p)`. Their 31-column binary
   matrix has rank 31.
6. A separate Brumer--Cremona quadratic-character implementation reaches
   binary rank 31 with 32 rows from 26 primes, ending at 283. This independently
   cross-checks both the point set and the finite-quotient certificate.

If an integral relation among the points existed, its images in every finite
quotient would force all coefficients to be even. Dividing by two produces a
rational 2-torsion point; because `E(Q)[2]=0`, the relation itself divides by
two. Infinite descent forces every coefficient to vanish. Hence

```text
rank E(Q) >= 31
```

unconditionally.

The finite-quotient proof is structurally different from the
Brumer--Cremona method used by the ICARM submission verifier. The local checker
also contains an independently written Brumer--Cremona replay as a second exact
certificate.

## Exact arithmetic fingerprint

The public equation has

```text
c4 = 61666932677451250957526286528757668905601968707960764583546763139025

c6 = -484158229736481797117926468811255292357508238788175012972498868670847390044003856563630907590456883225
```

and discriminant

```text
56986667162894850943626331069759320443067423220527305082244725946223212417010574363800489149279120357391198820328667683955054948591887684940751932057558936585975074155693048247805850350226862080000.
```

The exact `j`-invariant has SHA-256

```text
5939208330113d89ae063d62053f0c8383e18b3a564919b86f86a02a4d13a550.
```

The listed discriminant and conductor factorizations multiply back exactly.
All 26 distinct factors occurring in `c4` or the discriminant also pass
SymPy's primality checker. The only discriminant valuation at least 12 is
`v_2(Delta)=15`, while `c4` is odd, so the displayed integral equation is
already globally minimal.

The local fibre profile is particularly simple:

```text
p=5:  IV
all other bad p: I_{v_p(Delta)}
```

Thus the model is multiplicative at all bad primes except 5. The resulting
local conductor exponents recover the public conductor exactly: exponent 2 at
5 and exponent 1 at every other bad prime.

## Reproduction

From the repository root:

```bash
PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/verify_icarm_curve302_rank31.py \
  --verify-primality

PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/verify_icarm_curve302_rank31.py \
  --check --verify-primality
```

The pinned artifact is
[`icarm_curve302_rank31_v1.json`](../../artifacts/generated-results/elliptic-curves/icarm_curve302_rank31_v1.json).
Its SHA-256 is

```text
3be0d6fe82c58e0f9284df5d9340332944a1d906508ea986d4abe00357036991.
```

## Claim boundary

Proved here:

```text
rank E(Q) >= 31,
E(Q)_tors = 0,
the public model is globally minimal,
the published discriminant/conductor products and local fibre profile are exact.
```

Not proved here:

```text
rank E(Q) = 31 unconditionally,
curve 302 belongs to the reconstructed Elkies H3/rootless-MW17 family,
a rank-32 curve or a smaller rank-31 curve.
```

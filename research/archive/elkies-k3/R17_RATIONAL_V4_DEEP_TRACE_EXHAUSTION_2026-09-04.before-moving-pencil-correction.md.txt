# Complete rational-V4 deep-trace exhaustion (2026-09-04)

<!-- status-consumer: EC-K3-R17-NORM12-11952-COMPLETE-RATIONAL-V4-DEEP-TRACE-EXHAUSTION 8a0ccf906acc72a5 -->

## Status

This is an exact negative result for the remaining norm-twelve **integral
coboundary** layer of the alternate-Q80 rational `V4` construction.  It does
not construct a product-character section and therefore does not prove generic
rank at least 20.

The complete native atlas contains 39,147 smooth rational bisection classes.
For norm-ten vectors `w_i,w_j`, the corresponding bisections meet once exactly
when `<w_i,w_j>=7`.  Exact blockwise lattice pairing finds 4,358,409 such
unordered pairs.  Their distinct irreducible quadratic branch divisors are
coprime, so every pair gives a connected genus-one `V4` base.  The unique
intersection point is rational because both curves and their degree-one
intersection cycle are defined over `QQ`.

The product-Tate parity certificate leaves 49 minimum-norm-twelve trace
parities.  For each trace, every possible height-ten half-point carrier has

```text
M = M0 + lambda*h^2,       lambda in P1,
```

and degree-at-most-eight branch polynomial `q_lambda`.  It has the target
quartic character `d` precisely when `q_lambda=d*r^2` for a polynomial `r` of
degree at most two, including the rational constant squareclass.

The checker streams the literal branch triples from the 364-MiB equation
atlas; it does not primitive-normalize away their contents.  At `p=131` all 49
deep traces have good, nonzero reduction.  Factoring every member of their
projective parameter lines gives 49 distinct quartic squareclass keys.  An
exact vectorized scan compares those keys with all 4,358,409 rational-V4
products, or

```text
4,358,409 * 49 = 213,562,041
```

target--trace comparisons.  No pair survives.  A synthetic quartic constructed
from the first deep trace is recovered through the identical scalar-sensitive
hash path at its expected parameter, so the empty result is not caused by a
chart or constant-factor mismatch.  A second complete replay reproduces the
stored artifact byte for byte.

The generated certificate is
[`elkies-k3-r17-norm12-11952-v4-deep-trace-inversion-full-v1.json`](../artifacts/generated-results/elkies-k3-r17-norm12-11952-v4-deep-trace-inversion-full-v1.json),
SHA-256
`b3c81f6e5bb9a60b4169c857d8747bb606ee6424d5b3031096ff949a49ec4ccf`.

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 sage -python \
  elkies-k3/scripts/search_r17_norm12_11952_all_rational_v4_deep_trace_inversion.sage \
  --primes 131,137 \
  --output artifacts/generated-results/elkies-k3-r17-norm12-11952-v4-deep-trace-inversion-full-v1.json \
  --check
```

## Exact boundary

The theorem excludes a third character arising from any of the 49 deep
norm-twelve integral trace carriers for every rational intersection-one pair in
the complete smooth atlas.  Together with the separate norm-eight inversion,
this exhausts the zero-Tate-class height-eight carriers in those rational
`V4` targets.

It does not compute the product-twist Mordell--Weil groups, their Tate
quotients, or any nonzero Tate class.  A non-coboundary height-eight section, a
higher-height section, or a `V4` base built from rational bisections outside the
smooth native atlas remains open.  In particular,

```text
rank E/QQ(C) >= 20
```

is still `UNKNOWN`.

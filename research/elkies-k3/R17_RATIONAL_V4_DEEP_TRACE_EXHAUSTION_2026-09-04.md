# Regular-chart deep-trace exclusion for all rational-V4 targets

<!-- status-consumer: EC-K3-R17-NORM12-11952-COMPLETE-RATIONAL-V4-DEEP-TRACE-EXHAUSTION ff8d6c796122fe34 -->

## Status

This is an exact negative result for the **fixed regular-slope norm-twelve
families** on every target in the alternate-Q80 rational `V4` atlas. It does
not exhaust all norm-twelve integral coboundaries or construct a third
character. The [moving-pencil correction](Q80_ALL_SMOOTH_GENUS_ONE_BISECTIONS_2026-09-13.md)
identifies the missing finite zero-contact parameter. The
[earlier source](../archive/elkies-k3/R17_RATIONAL_V4_DEEP_TRACE_EXHAUSTION_2026-09-04.before-moving-pencil-correction.md.txt)
is preserved; the original computation and certificate are unchanged.

The complete native atlas contains 39,147 smooth rational bisection classes.
For norm-ten vectors `w_i,w_j`, the corresponding bisections meet once exactly
when `<w_i,w_j>=7`.  Exact blockwise lattice pairing finds 4,358,409 such
unordered pairs.  Their distinct irreducible quadratic branch divisors are
coprime, so every pair gives a connected genus-one `V4` base.  The unique
intersection point is rational because both curves and their degree-one
intersection cycle are defined over `QQ`.

The product-Tate parity certificate leaves 49 minimum-norm-twelve trace
parities. For each trace, the calculation tests the regular-slope family

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

The theorem excludes a third character arising from these 49 fixed
regular-slope families for every rational intersection-one pair in the
complete smooth atlas. The target coverage is complete; the norm-twelve
carrier coverage is not. In particular, the full pencil of `O+tau-F` has
slope `M0/h-c*h/(t-u)`, with `c=coeff(M0,t^7)`, and allows its unique zero
contact to move to finite t. Those members are not covered by a constant
lambda in the displayed regular-slope formula. Thus this calculation and
the separate norm-eight inversion do not by themselves exhaust every
zero-Tate-class height-eight carrier in the rational `V4` atlas.

It does not compute the product-twist Mordell--Weil groups, their Tate
quotients, or any nonzero Tate class. Moving-contact coboundaries,
non-coboundary height-eight sections, higher-height sections and other
bisection atlases remain outside this comparison. The separate
[arithmetic rank-zero theorem for seventeen selected products](R17_PRODUCT_REGULATOR_OBSTRUCTION_SWEEP_2026-09-05.md)
already closes those particular targets at every height; this correction
does not reopen them. For the general integral `V4` target,

```text
rank E/QQ(C) >= 20
```

is still `UNKNOWN`.

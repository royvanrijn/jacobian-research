# Rational-normalization boundary beyond the smooth bisection atlas (2026-09-04)

<!-- status-consumer: EC-K3-R17-NORM12-SINGULAR-GENUS1-RATIONAL-NORMALIZATION-EXHAUSTION bf05d9b06ccc1502 -->
<!-- status-consumer: EC-K3-R17-NORM12-GENUS2-RATIONAL-NORMALIZATION-SCREEN 08357122de3c43ff -->

## Result

The branch-character search has no failure in the first singular row on either
direct norm-twelve rank-17 chart.

For `norm12-orbit-11952` (alternate Q80), all 63,917 minimum norm-eight
translation classes were checked.  For the hidden `norm12-orbit-103b2` chart,
all 63,925 classes were checked.  Neither search contains a nonsplit rational
singular member, hence neither produces a rational quadratic normalization,
a character already present in the complete smooth atlas, a collision between
new characters, or a three-character relation.

Combined with the complete smooth-atlas injectivity certificates, this proves
injectivity through the arithmetic-genus-one row having rational
normalization.  The least possible arithmetic genus for a failure in this
rational-normalization filtration is therefore **at least two**.  It is not
currently known to equal two: the genus-two computation below is an exact
multi-prime and bounded-CRT miss, not a global characteristic-zero
nonexistence theorem.

The smooth elliptic-normalization part of arithmetic genus one has quartic
branch support and is not a rational quadratic base.  It is outside the
rational-normalization assertion above.

## Exact singular genus-one certificate

Let `m` be the exact number of minimum norm-eight representatives up to sign
in a parity class.  The regular chord gauge puts one known split member at
infinity, with `q_infinity=h^2`; the other `m-1` split members are finite.  For
every class, the finite pencil discriminant has degree 22 and its even part
has degree exactly `2*(m-1)`.  Thus the known split members exhaust the entire
even-multiplicity part, including possible cuspidal or tangential cases that
would otherwise disappear on taking the polynomial squareclass.

The remaining odd part has no rational projective root.  This is certified
class by class by exact factorization, exact irreducibility, or a good-prime
projective-root obstruction that also retains the required full and odd
discriminant degrees.  Any accepted rational quadratic normalization would
additionally have had to pass exact characteristic-zero factorization,
quadratic squarefreeness, the two section identities, and comparison with the
complete smooth squareclass manifest.

The merged counts are:

| chart | classes | exact factorization | exact odd irreducibility | ordinary modular obstruction | full-discriminant modular obstruction | candidates |
|---|---:|---:|---:|---:|---:|---:|
| alternate Q80 | 63,917 | 6,517 | 53,730 | 3,520 | 150 | 0 |
| hidden 103b2 | 63,925 | 1,288 | 1 | 2 | 62,634 | 0 |

The hidden run had three long-lived modular survivors.  Exact resolution
excluded indices 3,091 and 43,478 at primes 1013 and 1009, respectively, and
proved that the degree-20 odd discriminant at index 13,206 is irreducible over
`QQ`.

Primary merged certificates:

- [`elkies-k3-r17-norm12-11952-singular-bisection-search-complete-v1.json`](../artifacts/generated-results/elkies-k3-r17-norm12-11952-singular-bisection-search-complete-v1.json), SHA-256 `2e20206614e17cad992e98742f917c6a50abccd39799a618ad6126fc62733f13`;
- [`elkies-k3-r17-norm12-103b2-singular-bisection-search-complete-v1.json`](../artifacts/generated-results/elkies-k3-r17-norm12-103b2-singular-bisection-search-complete-v1.json), SHA-256 `a35bf3c64fa889036d41097f193026409dde3ee44882e978badd69f21e5b7a13`.

The hidden multiplicity table is generated independently by
[`rank_r17_norm12_103b2_norm8_pencils.sage`](scripts/rank_r17_norm12_103b2_norm8_pencils.sage).
Its exact Fincke--Pohst traversal counts 516,046 signed vectors through norm
eight, of which 460,080 lie on the norm-eight shell, and leaves precisely
63,925 parity classes after removing lower-norm cosets.

The search and merge programs are
[`search_r17_norm12_direct_singular_bisections.sage`](scripts/search_r17_norm12_direct_singular_bisections.sage),
[`search_r17_norm12_direct_norm8_singular_modp.sage`](scripts/search_r17_norm12_direct_norm8_singular_modp.sage),
and
[`merge_r17_norm12_direct_norm8_singular_search.py`](scripts/merge_r17_norm12_direct_norm8_singular_search.py).
The exact command lines and input hashes are stored in every shard and merged
certificate.

## Two-node genus-two normalization screen

The norm-six trace shells contain 26,645 classes on alternate Q80 and 26,672
on hidden 103b2.  The four-essential-variable systems were screened at
`p=17,23,29,31,37`, with full finite-field factorization and scalar
squareclass retained.

Intersecting cover identities in the order `23,29,31,17,37` gives:

| chart | smooth-versus-norm-six intersections | distinct norm-six intersections |
|---|---|---|
| alternate Q80 | 3,372,368 -> 10,748 -> 36 -> 0 -> 0 | 3,136,622 -> 30,837 -> 264 -> 4 -> 0 |
| hidden 103b2 | 3,011,174 -> 7,732 -> 19 -> 1 -> 0 | 2,993,900 -> 18,789 -> 113 -> 1 -> 0 |

CRT reconstruction modulo 13,005,833 then tested every reconstructed rational
parameter exactly over `QQ`.  Alternate Q80 had 502 common traces, 314,458 CRT
tuples, and 78,860 distinct reconstructed parameter pairs.  Hidden 103b2 had
610 common traces, 427,063 tuples, and 107,633 pairs.  Exact factorization and
section identities accepted zero candidates on both charts.

The compact certificates are:

- alternate cover intersection,
  [`elkies-k3-r17-norm12-11952-genus2-cover-collision-intersection-v1.json`](../artifacts/generated-results/elkies-k3-r17-norm12-11952-genus2-cover-collision-intersection-v1.json), SHA-256 `21cb2d2ee98ab1fc160e176e2374b73a3eda09c4afa5e2aaa0bcaa344253a994`;
- hidden cover intersection,
  [`elkies-k3-r17-norm12-103b2-genus2-cover-collision-intersection-v1.json`](../artifacts/generated-results/elkies-k3-r17-norm12-103b2-genus2-cover-collision-intersection-v1.json), SHA-256 `c111cd8cdf0f013f6f98be70cbd1d5b0c1ef3ef37cc2fe2a9d02f094a1b722d3`;
- alternate CRT reconstruction,
  [`elkies-k3-r17-norm12-11952-genus2-normalization-reconstruction-full-v1.json`](../artifacts/generated-results/elkies-k3-r17-norm12-11952-genus2-normalization-reconstruction-full-v1.json), SHA-256 `4618a2b9086de0c2802b083d5afe296819bd26918d4b3cf845cd157e5bdf7229`;
- hidden CRT reconstruction,
  [`elkies-k3-r17-norm12-103b2-genus2-normalization-reconstruction-full-v1.json`](../artifacts/generated-results/elkies-k3-r17-norm12-103b2-genus2-normalization-reconstruction-full-v1.json), SHA-256 `c9fcb1ac3f15e0835761d2eb0f0411e5f444f42673da1adb9d24fb7ca6be75d8`.

This genus-two result is deliberately fail-closed.  The empty multi-prime
intersection excludes characteristic-zero collisions only under simultaneous
good reduction in the displayed affine charts.  The reconstruction excludes
only rational parameters inside the standard reconstruction bound and
integral at all five primes.  Parameter-at-infinity, bad-reduction charts, and
unbounded rational heights remain open.  Accordingly no unconditional
arithmetic-genus-two injectivity theorem, rank-19 rational quadratic base, or
rank-20 `V4` construction is claimed.

## Current boundary

The answer presently certified is

\[
  p_a^{\rm first\ possible\ failure}\ge 2,
\]

for the rational-normalization branch-character filtration on both direct
rank-17 charts.  Whether equality holds, or whether injectivity continues
globally through arithmetic genus two, remains `UNKNOWN`.

# R17 multisection diversity beyond raw counts (2026-09-02)

## Result

Raw bisection and trisection counts discard substantial structure.  On the
published rootless R17 frame, the first orbit/metric/overlap profile gives
three concrete corrections to the count-only picture:

1. `Aut(R17)={+I,-I}` is too small to compress the bisection atlas: all
   degree-two cosets are fixed.  The 39,120 rational bisection vertices are
   nevertheless organized by one connected zero-intersection graph with
   8,895,801 edges and 157,553,175 triangles.
2. Threshold counts miss deep cosets.  In degree two there are 43 cosets of
   minimum norm 12 beyond the 63,925 genus-one frontier vertices of norm 8.
   The complete degree-three census similarly has 320 rational vertices of
   norm 26 and 3,388 genus-one vertices of norm 24 beyond its first qualifying
   shells.
3. Degree overlap changes geometric labels.  Under the natural inclusion
   `M/2M -> M/4M`, all 39,120 rational bisection vertices become genus-one
   quadrisection vertices of minimum norm 40.  They are therefore absent from
   the old degree-four low-height cap 38 even though they are the most direct
   inherited degree-four mechanisms.
4. The complete comparison with `NS0032-F011` and `NS0028-F005` reverses
   several count-only rankings.  In particular, NS0032 has the most rational
   vertices but only 20,933 automorphism orbits, while NS0028 has 41,376;
   NS0028 has more vertices than R17 but fewer zero-intersection triangles.

This supports replacing the scalar pair `(N_2,N_3)` by a **multisection
diversity profile**.  A single new scalar is premature: on R17, orbit entropy
is almost vacuous because the automorphism group has order two, whereas
quotient separation, deep-hole mass, equation cost, and local squareclass
entropy remain nontrivial.

The degree-two results below are complete computational lattice results.  The
degree-three one-vertex spectrum is the existing complete census, but its
graph is sampled.  Degree four is sampled except for the exact embedded
two-torsion calculation.  Only the already certified rational bisection atlas
has equation and squareclass meaning.

## Intrinsic graph

For `c in M/dM`, put

```text
mu_d(c) = min { (w,w) : w == c mod dM }.
```

Proposition F6 in
[`RANK_MUTATION_AND_LIFT_THEOREMS.md`](RANK_MUTATION_AND_LIFT_THEOREMS.md)
proves that degree-`d` vertices of arithmetic genera `g,h` have minimum
intersection, after independent section translations,

```text
mu_d(c-c')/2 + g + h - 2.
```

This is the correct quotient substitute for a raw angle.  It is invariant
under section translation and `Aut(M)`.  By contrast, an angle between two
chosen shortest vectors changes when a different shortest lift is selected.
The artifact retains an angle histogram in the deterministic
equation-priority gauge, but labels it as representative-dependent.

The first graph uses an edge when the minimum intersection is at most one.
For the rational degree-two subgraph, the sharper zero-intersection graph is
the previously studied norm-four XOR graph.

## Complete degree-two structure

The minimum-norm distribution on all `2^17` cosets is

| minimum norm | cosets |
|---:|---:|
| 0 | 1 |
| 4 | 1,311 |
| 6 | 26,672 |
| 8 | 63,925 |
| 10 | 39,120 |
| 12 | 43 |

Thus the discrete covering radius is exactly

```text
rho_2 = sqrt(12)/2 = sqrt(3).
```

The qualifying vertex sets are:

| label | threshold | vertices | deep vertices | `Aut(M)` orbits | span over `F_2` |
|---|---:|---:|---:|---:|---:|
| rational | 10 | 39,120 | 0 | 39,120 | 17 |
| genus one | 8 | 63,968 | 43 at norm 12 | 63,968 | 17 |

The older value 63,925 counted only the minimum-height norm-eight frontier.
For diversity, all 63,968 section-nonnegative integral genus-one cosets should
be retained, with depth as a weight rather than as an exclusion.

Among the `binomial(39120,2)=765,167,640` rational pairs, the exact quotient
separation/intersection distribution is

| `mu_2(c-c')` | minimum intersection | pairs |
|---:|---:|---:|
| 4 | 0 | 8,895,801 |
| 6 | 1 | 161,556,656 |
| 8 | 2 | 372,003,497 |
| 10 | 3 | 222,477,240 |
| 12 | 4 | 234,446 |

The zero-intersection graph is connected.  Its degrees range from 342 to 854,
with median 454, so even one connected mechanism space is highly
nonuniform.  It contains exactly 157,553,175 graph triangles.  A triangle
means that every pair has quotient minimum intersection zero; it does **not**
prove that one simultaneous choice of three translated representatives is
pairwise disjoint.

This is already a stronger discriminator than `N_2`: two lattices with the
same vertex count can have different separation entropy, degree distribution,
connectivity, clique counts, and deep-hole mass.

## Degree three and four

The complete degree-three covering radius is `sqrt(26)/3`.  Its low-genus
sets split as follows:

| label | vertices | first shell | deeper shells | `Aut(M)` orbits |
|---|---:|---:|---:|---:|
| rational | 18,024,616 | 18,024,296 at norm 20 | 320 at norm 26 | 9,012,308 |
| genus one | 33,484,468 | 33,481,080 at norm 18 | 3,388 at norm 24 | 16,742,234 |

Since inversion has no fixed nonzero point modulo three, every orbit has size
two.  A deterministic inversion-closed sample of 1,025 cosets spans all 17
dimensions over `F_3` in both qualifying labels.  Its graph statistics are
sampling diagnostics only; the full graph would require a separate streamed
or Fourier calculation rather than materializing tens of millions of
vertices.

The corresponding degree-four sample is also inversion-closed and spans all
17 dimensions modulo two in each of the rational, genus-one, and genus-two
labels.  Its observed maximum coset minimum is 40, but the exact degree-two
embedding below improves the covering-radius lower bound to `sqrt(3)` by
exhibiting norm-48 degree-four cosets.  No complete degree-four covering radius
is claimed.

## Exact overlap across degrees

View `M/dM` as the `d`-torsion subgroup `(1/d)M/M` of the lattice torus.  Then

```text
T_d intersect T_e = T_gcd(d,e).
```

Consequently the displayed positive-threshold degree-two and degree-three
sets, and the degree-three and degree-four sets, have no literal common
vertices: the corresponding torsion subgroups meet only at zero.  There is no
canonical reduction map between coprime degrees.  Any stronger `d=2` versus
`d=3` overlap statistic must declare a common `M/6M` CRT metric.

For `2|4`, however, there is a canonical inclusion

```text
c mod 2M  ->  2c mod 4M,
mu_4(2c) = 4*mu_2(c).
```

It yields 103,088 exact genus-one quadrisection vertices inherited from the
degree-two high-minimum region:

| degree-two source | count | degree-four minimum | degree-four label |
|---|---:|---:|---|
| genus-one frontier | 63,925 | 32 | genus one, frontier |
| rational bisections | 39,120 | 40 | genus one, deep |
| norm-12 deep holes | 43 | 48 | genus one, deeper |

This is the most important immediate correction to the old degree-four
sample: a height cap of 38 systematically removes all inherited rational
bisection mechanisms.

## Equation complexity and squareclass diversity

For the complete rational-bisection atlas, minimum norm and genus are
constant, but equation cost is not.  The exact priority table has group-law
addition upper bounds from 2 to 17 (median 8), coordinate-input cost from
1,616 to 5,611 bits (median 3,184), and graph degree from 342 to 854.  These
weights give a practical way to choose representatives without pretending
that every one of the 39,120 vertices is equally cheap.

All 39,120 global classes in `QQ(t)^*/QQ(t)^{*2}` are distinct.  The new local
profile reduces each class at every consecutive odd prime from 97 through
163, recording the Gauss valuation parity and squarefree residual class in
`F_p(t)`.  Per-prime entropy varies substantially, so one local place is not a
stable diversity statistic.  The joint signature across the declared prime
block separates all 39,120 classes: every joint bucket is a singleton and its
entropy is `log_2(39120)=15.2556` bits.  This is local squareclass separation,
not proof that two covers give different specialization quotient directions
at a chosen rational fibre.

## Recommended profile

For each degree and genus, retain the following vector of invariants rather
than only its cardinality:

```text
(Aut-orbit entropy,
 covering radius and depth histogram,
 span dimensions mod primes dividing d,
 quotient-separation/intersection entropy,
 small-intersection graph components and clique counts,
 equation-complexity distribution,
 local squareclass entropy over several places,
 cross-degree inherited mass and genuinely new mass).
```

Only after this vector has been calibrated against the rank-25--28 positive
controls and several non-controls should it be compressed into a scalar
`multisection diversity` score.  Equation and local-squareclass coordinates
should remain R17-only until the foundry targets have equation models.

## Complete foundry comparison

The prescribed first comparison is now complete on all `2^17` cosets of R17,
`NS0032-F011`, and `NS0028-F005`.  The exact common-coordinate summary is:

| invariant | R17 | NS0032-F011 | NS0028-F005 |
|---|---:|---:|---:|
| rational vertices | 39,120 | 41,421 | 41,376 |
| rational `Aut(M)` orbits | 39,120 | 20,933 | 41,376 |
| vertex-weighted orbit entropy (bits) | 15.2556 | 14.3488 | 15.3365 |
| all genus-one vertices | 63,968 | 64,064 | 64,082 |
| norm-12 deep cosets | 43 | 315 | 142 |
| rational-pair separation entropy (bits) | 1.5760 | 1.5782 | 1.5665 |
| zero-intersection edges | 8,895,801 | 9,220,380 | 8,987,726 |
| zero-intersection edge density | 0.011626 | 0.010749 | 0.010500 |
| median graph degree | 454 | 450 | 433 |
| graph triangles | 157,553,175 | 162,074,080 | 151,106,959 |
| graph transitivity | 0.116384 | 0.117291 | 0.115480 |
| inherited genus-one degree-four vertices | 103,088 | 105,485 | 105,458 |

All three rational graphs are connected, all rational and genus-one vertex
sets span 17 dimensions over `F_2`, and all three discrete covering radii are
`sqrt(3)`.  Those coordinates do not distinguish this batch.  The remaining
coordinates do, in mutually incompatible orders:

- Raw rational abundance orders the lattices NS0032, NS0028, R17.
- Automorphism-orbit abundance orders them NS0028, R17, NS0032.
- Edge density orders them R17, NS0032, NS0028.
- Raw triangle abundance orders them NS0032, R17, NS0028.
- Pair-separation entropy orders them NS0032, R17, NS0028, while deep-hole
  mass orders them NS0032, NS0028, R17.

The NS0032 symmetry effect is large, not cosmetic.  Its order-four integral
automorphism group induces an order-two action modulo two: the rational set
has 445 fixed vertices and 20,488 two-element orbits.  By contrast, inversion
acts trivially modulo two for the order-two groups of R17 and NS0028, so every
degree-two vertex there is a singleton orbit.

Deep cosets also alter an earlier count-only conclusion.  Counting only the
norm-eight genus-one frontier gives NS0028 `63,940`, R17 `63,925`, and NS0032
`63,749`.  Adding every qualifying norm-12 coset changes the order to NS0028
`64,082`, NS0032 `64,064`, and R17 `63,968`.  The exact inherited degree-four
mass then leaves NS0032 and NS0028 only 27 vertices apart despite their very
different automorphism-orbit profiles.

This is the requested proof of concept: two candidates that differ by only 45
raw rational vertices have radically different symmetry-reduced mechanism
counts, and their graph/depth coordinates do not collapse to the same ranking.
No one-dimensional score is justified yet.

## Specialization-response calibration

The global graph describes the available geometric supply, but a rank jump is
controlled by what happens after specialization.  The complete bisection
specialization certificate identifies every R17 bisection that splits at the
rank-21 mechanism control and the four rank-25--28 controls, together with its
class in the displayed exceptional quotient.  This gives an exact two-stage
picture:

```text
geometric supply S_2 in M/2M
        -> split subset Split_t
        -> exceptional response phi_t(Split_t) in (L_t/M_t) tensor F_2.
```

The first arrow is arithmetic splitting; the second records which new
Mordell--Weil direction the split cover exposes.  Their exact profiles are:

| fibre | displayed quotient rank | split covers | lattice-mask span | exceptional span | visibility fraction | intersection-0 edges | intersection-at-most-1 edges | largest near clique |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| rank-21 mechanism control | 4 | 25 | 17 | 4 | 1.000 | 4 | 73 | 4 |
| rank at least 25 | 8 | 6 | 6 | 5 | 0.625 | 0 | 6 | 3 |
| rank at least 26 | 9 | 3 | 3 | 3 | 0.333 | 0 | 0 | 1 |
| rank at least 27 | 10 | 2 | 2 | 2 | 0.200 | 0 | 0 | 1 |
| rank at least 28 | 11 | 1 | 1 | 1 | 0.091 | 0 | 0 | 1 |

This changes the interpretation of “independent mechanisms.”  At the rank-21
control, 25 split covers span all 17 dimensions of the geometric lattice
quotient but map to only four exceptional directions; four covers even map to
the zero exceptional class.  At rank 25, six lattice-independent masks map to
a five-dimensional exceptional space.  Thus lattice span is an upper supply
coordinate, not arithmetic independence.

Nor is small geometric intersection necessary for independent arithmetic
directions.  The three rank-26 covers and two rank-27 covers give independent
exceptional classes, yet every pair has minimum intersection at least two.
The rank-25 set has six intersection-at-most-one edges and one near clique of
size three, but all fifteen of its pairs also join independent exceptional
directions.  At the rank-21 control, the mutual information between quotient
separation norm and exceptional pair type is only `0.0331` bits.  This small,
selected sample does not estimate a population law, but it directly refutes
identifying graph proximity with exceptional independence.

All 37 split-cover occurrences across the five fibres are different lattice
vertices.  This cross-fibre novelty is real, even though within each fibre the
response can collapse many covers onto few directions.  The appropriate
discovery object is therefore not one scalar but a pair:

```text
(global multisection supply diversity,
 fibre-specific splitting/exceptional-response diversity).
```

Across these historically selected controls, split-cover count and displayed
quotient rank are perfectly oppositely ordered.  That descriptive fact has
`n=5`, is selection-biased, and is not a predictive correlation.  It does show
that maximizing global or specialized cover count alone would optimize the
wrong response for the known high-rank fibres.

## Degree-three stopping boundary

The next nominal step was an exact degree-three graph.  The present complete
degree-three checkpoints retain histograms and checksums, not one membership
bit per coset, so they cannot answer adjacency queries.  Recovering the
rational indicator would replay `64,570,082` inversion representatives for
R17.  A packed indicator on all `3^17=129,140,163` cosets is only about 15.4
MiB, but the small-intersection kernel already contains 2,622 norm-four and
53,344 norm-six difference cosets.  Naive packed shifts would therefore read
roughly 840 GiB before triangle or degree-distribution accounting; a ternary
Fourier convolution needs multi-gigabyte working arrays and an exact rounding
audit.

More importantly, that calculation would refine only the global supply side.
There is no equation-level trisection atlas at the five controls, so it cannot
currently be calibrated against the exceptional-response map that the
degree-two result shows to be essential.  The current 1,025-coset
automorphism-closed degree-three graph remains a correctly labelled sample.
Expanding it now has diminishing value.

Resume the degree-three graph programme when at least one of these becomes
available:

1. a census replay that emits packed per-coset genus labels while CVP work is
   already being performed;
2. an equation-level trisection atlas on R17 that can be specialized at the
   controls; or
3. a substantially larger, prospectively selected set of fibres with exact
   rank/quotient labels.

Until then, the two-stage degree-two calibration is the natural stopping
point: it answers the conceptual question and prevents an expensive global
graph from being mistaken for an arithmetic rank predictor.

## Literature boundary

Garbagnati--Salgado provide the geometric reason to study special
multisections in connection with rank jumps, but they do not define the
finite coset-graph diversity profile used here:

- A. Garbagnati and C. Salgado,
  [*Rank jumps and Multisections of elliptic fibrations on K3 surfaces*](https://arxiv.org/abs/2505.15159).

Automorphism-constrained angle and clique statistics have close analogues in
the spherical-code literature, but that analogy supplies search diagnostics,
not a K3 rank theorem:

- M. Ganzhinov and P. R. J. Oestergard,
  [*Spherical codes with prescribed signed permutation automorphisms inside shells of low-dimensional integer lattices*](https://arxiv.org/abs/2403.04330).

The literature check therefore supports the ingredients, not a claim that
the boxed diversity profile is a known invariant or predicts arithmetic rank.

## Reproduction

Run:

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/analyze_r17_multisection_diversity.py
```

This writes
[`../artifacts/generated-results/elkies-k3-r17-multisection-diversity-v1.json`](../artifacts/generated-results/elkies-k3-r17-multisection-diversity-v1.json).
Its SHA-256 is
`e4f1499a9aad5670202b5948717995d5696116979b9eb7701bc87daa8175c512`.
It pins every input hash, independently MPFR-audits all 43 norm-12
degree-two holes plus a deterministic residue stride, and marks the sampled
degree-three/four graph fields explicitly.

Run the complete three-lattice degree-two comparison with:

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/analyze_r17_multisection_diversity.py \
  --comparison-only
```

This writes
[`../artifacts/generated-results/elkies-k3-r17-foundry-degree2-diversity-comparison-v1.json`](../artifacts/generated-results/elkies-k3-r17-foundry-degree2-diversity-comparison-v1.json).
Its SHA-256 is
`408947ebb3e67048767005a005c3f283cc2bd2b4971e12006716809134a7146c`.
The run recomputes every minimum norm integrally after double-double CVP,
audits a deterministic stride and every deepest coset at 256-bit precision,
and obtains the full automorphism groups with PARI `qfauto`.

After creating that global comparison artifact, run the exact five-control
response calibration with:

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/analyze_r17_multisection_diversity.py \
  --control-calibration-only
```

This writes
[`../artifacts/generated-results/elkies-k3-r17-bisection-control-diversity-calibration-v1.json`](../artifacts/generated-results/elkies-k3-r17-bisection-control-diversity-calibration-v1.json),
with SHA-256
`57945b3431cec5227ff245ae0fa238576529f23ec95bedc9f2667b7d51bc742e`.
Every induced edge, component, and clique was independently recounted by a
separate brute-force check.  The correlations are explicitly descriptive;
the underlying cover, point, and quotient-class certificates retain their own
mathematical status.

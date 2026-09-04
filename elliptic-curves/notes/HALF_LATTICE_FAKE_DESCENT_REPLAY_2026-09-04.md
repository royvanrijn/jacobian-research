# Half-lattice fake-descent replay (2026-09-04)

## Status

The mandatory rank-28 gate passes.  A search that loaded only the published
R17 equation and its seventeen sections recovered a subgroup of exact rank
28 on the fibre (T=-9529/5471).  Only after that output was frozen did the
verification program load the eleven published exceptional points; the
recovered quotient has dimension 11 and equals their full quotient modulo the
specialized generic subgroup.

The replay also changes the interpretation of the method.  The quartics used
here are pointed models birational to the elliptic curve, not nontrivial
2-covering torsors.  They are automatically locally soluble at every place.
Thus the computation establishes that deep half-lattice classes are unusually
effective *search charts*.  It does not establish that record fibres are
explained by simultaneous local solubility of several Selmer classes, and the
present features do not yet predict a rank jump before point search.

An equal-budget selection ablation now separates the two main effects.  On the
three sealed +12 holdouts, the generic top 43 recover mean exact quotient rank
8.33, versus 3.20 averaged over five deterministic random 43-sets, and have
the best pooled rank-per-CPU rate.  Every shallowest-43 arm recovers rank zero
on all seven tested fibres.  Reduction and multi-chart search are still
powerful--random and median charts often work on the easier controls--but
trying many reduced charts alone does not explain the deep-hole result.

The compact machine-readable overview is
[`half_lattice_fake_descent_experiment_matrix_v1.json`](../../artifacts/generated-results/elliptic-curves/half_lattice_fake_descent_experiment_matrix_v1.json).
The equal-budget ablation certificate is
[`half_lattice_search_ablation_summary_v1.json`](../../artifacts/generated-results/elliptic-curves/half_lattice_search_ablation_summary_v1.json).

## Claim ledger

### Exact results

- The seventeen specialized R17 sections are independent: an exact
  finite-reduction certificate has rank 17.
- Every reported quartic point maps by exact rational functions to a point on
  the original specialized curve, and the curve equation is checked exactly.
- In the ablation, every retained candidate has an exact relation in the
  displayed public subgroup.  For every arm on every fibre, the exact
  quotient rank over \(\mathbf Q\) equals the independently retained mod-2 rank.
- On the rank-28 fibre the blind candidates increase the exact
  finite-reduction rank from 17 to 28.
- Post-search fixture verification proves that all 40 candidates outside the
  generic subgroup lie in the published rank-28 group and that their quotient
  span is the full published 11-dimensional complement.
- The complete exact integral CVP census of (M/2M) in the generic R17 height
  form is

  | minimum norm | 0 | 4 | 6 | 8 | 10 | 12 |
  |---:|---:|---:|---:|---:|---:|---:|
  | number of classes | 1 | 1,311 | 26,672 | 63,925 | 39,120 | 43 |

  Hence the 43 generic deepest classes have half-lattice depth 3.  The terminal
  stratum independently reproduces Proposition 8 of
  [Elkies 2026](https://arxiv.org/pdf/2608.25406).

### Numerical lattice evidence

- The specialized canonical-height matrix is evaluated to high precision.
  CVP representatives were recomputed after scaling and rounding it by
  (10^4,10^5,10^6); the top-43 class set is identical at all three scales.
  This is reproducibly stable numerical evidence, not an interval proof of the
  specialized ordering.
- At rank 28, only 22 of the exact generic deepest 43 classes are also in the
  specialized top 43.  Their union has 64 classes.
- Arbitrary subgroups on curves 273 and 302 are less numerically stable: 14--24
  of 43 representatives change between the checked rounding scales.  Their
  depths should not be promoted to exact lattice claims.
- All seven ablation fibres reproduce the specialized top-43 *set* after
  rounding the canonical-height matrix at both (10^5) and (10^6).  This is
  a numerical stability check, not an exact ordering proof.

### Bounded-search results

Unless stated otherwise, each reduced quartic was searched once with PARI
`hyperellratpoints` through reduced-coordinate height 100,000, a 15-second
per-model timeout, and no retry.  A miss is not a proof of absence.

- Searching the 64-class generic/specialized union at rank 28 found 179
  distinct nonbasis points, including 40 outside (M), and certified quotient
  gain 11.
- The fixed generic-deepest-43 search recovers gains (4,8,9,10,9) on the
  known rank (21,25,26,27,28) fibres.  It therefore recovers the full known
  complement through rank 27 and 9 of 11 directions at rank 28.  The
  specialized/generic union supplies the missing two rank-28 directions.
- On three exactly transported rank-29 controls, the same union search recovers
  10 of 12 directions on the 2024 curve, 12 of 12 on ICARM curve 356, and 4 of
  12 on ICARM curve 385.
- With equal per-cover budgets, the generic top-43 arm recovers quotient ranks
  (8,9,10,9) on the +8 through +11 development fibres and (10,12,3) on the
  sealed +12 holdouts.  Five deterministic random arms average respectively
  (8,9,8.8,4.6) and (2.2,5.4,2.0).  The shallowest arm recovers zero in all
  seven cases.

### Heuristic interpretation

- Deep holes are privileged for bounded recovery: the rank-28 quotient is
  found by searching 64 of 131,072 parity classes, and the held-out ablation
  shows a 2.60-fold generic-depth enrichment over the mean random rank.  This
  is search concentration, not a theorem that every exceptional point belongs
  to one canonical half-class.
- Generic and specialized depths are complementary at rank 28: either top-43
  list gives quotient rank 9 and their union gives 11.  Across the sealed
  holdouts, however, generic depth is the stronger selector (25 total
  directions versus 18); specialization-specific CVP is not the main effect.
- Depth is not a calibrated prospective rank predictor.  Raw point counts,
  reduced coefficient size, and modular square density all fail to separate
  records from controls reliably at the tested budget.

## Equal-budget class-selection ablation

The arm rules were frozen after development on the +8 through +11 R17 fibres.
Every non-union arm has exactly 43 covers.  A cover always uses the
specialized-height shortest representative of its class, followed by the same
PARI minimization, reduction, and one `hyperellratpoints` call at height
100,000 with a 15-second timeout.  The five random sets are disjoint initial
blocks of one SHA-256 ordering of the 131,071 nonzero masks.  Shallow and
median arms are fixed SHA-256 samples within their exact generic-depth strata.
All five random sets happen also to be disjoint from both top-43 sets on every
tested fibre; this was measured after selection, not imposed by rejection.
The union is charged for every distinct cover and is normalized back to 43
covers in the table.

Each `G` and `S` cell is `exact quotient rank / rank per CPU second`; because
these arms contain 43 covers, the rank is also rank per 43 covers.  Each `U`
cell is `rank / class count; normalized rank per 43; rank per CPU second`.
`R1...R5` gives all five exact random-arm ranks and their pooled CPU rate.
All CPU figures include parent and completed child user-plus-system time.

| fibre | G: generic top 43 | S: specialized top 43 | U: union | R1...R5 ranks; pooled rate | median 43 | shallow 43 |
|---|---:|---:|---:|---:|---:|---:|
| R17 +8 | 8 / .415 | 8 / .400 | 8/86; 4.000; .204 | 8,8,8,8,8; .446 | 8 / .424 | 0 / 0 |
| R17 +9 | 9 / .600 | 9 / .580 | 9/80; 4.838; .319 | 9,9,9,9,9; .575 | 9 / .541 | 0 / 0 |
| R17 +10 | 10 / .586 | 10 / .636 | 10/75; 5.733; .350 | 9,8,8,10,9; .576 | 8 / .533 | 0 / 0 |
| R17 +11 | 9 / .452 | 9 / .409 | 11/64; 7.391; .354 | 5,3,6,5,4; .224 | 5 / .240 | 0 / 0 |
| alternate-Q80 curve 12 +12 | 10 / .603 | 5 / .309 | 10/72; 5.972; .370 | 0,5,3,1,2; .150 | 1 / .068 | 0 / 0 |
| ICARM 356 +12 | 12 / .724 | 10 / .603 | 12/75; 6.880; .416 | 7,4,4,5,7; .370 | 7 / .465 | 0 / 0 |
| ICARM 385 +12 | 3 / .172 | 3 / .174 | 4/75; 2.293; .135 | 2,0,3,1,4; .123 | 2 / .115 | 0 / 0 |

The sealed holdouts give the cleanest comparison:

| arm | mean exact rank | mean rank per 43 | pooled rank / CPU second |
|---|---:|---:|---:|
| generic top 43 | 8.333 | 8.333 | .494 |
| specialized top 43 | 6.000 | 6.000 | .360 |
| generic/specialized union | 8.667 | 5.049 | .304 |
| five random arms, averaged/pooled | 3.200 | 3.200 | .211 |
| median-depth 43 | 3.333 | 3.333 | .212 |
| shallowest 43 | 0 | 0 | 0 |

This supports four bounded-search conclusions, kept deliberately separate
from theorem claims:

1. Generic deep-hole geometry is a real selector.  It wins decisively on
   curve 12 and curve 356, and its aggregate holdout efficiency is 2.34 times
   the pooled random efficiency.  Curve 385 is an important adverse case: the
   best random arm finds rank 4 while generic depth finds rank 3.
2. Specialized CVP is not the principal selector.  It ties generic depth in
   development, but loses (25-18) in total exact holdout rank.
3. Minimize/reduce plus multi-chart search is independently powerful.  Random
   and median sets already recover the entire quotient on the +8 and +9
   development fibres, but their yield collapses on the harder controls.
4. The union is a recall policy, not an efficient 43-cover policy.  It is
   mandatory for full +11 recovery and adds one direction on curve 385, but it
   is worse per 43 covers and per CPU second than generic depth on the sealed
   holdouts.

There was one fail-closed holdout implementation stop.  Alternate-Q80 curve 12
has 49, not 43, classes in its maximum generic-depth stratum.  An assertion
that this stratum had exactly 43 elements aborted before *any* holdout quartic
was searched.  The only source change removed that assertion; the already
active rule “sort by norm, then mask, take 43,” all seeds, representatives, and
budgets were unchanged.  The repaired search source hash and the pre-search
development hash are both recorded in the compact certificate.  Public
holdout points were first imported only after the final blind artifact was
hashed.

## Blind boundary

The search executable
[`half_lattice_fake_descent_replay.sage`](../cas/half_lattice_fake_descent_replay.sage)
imports only the published R17 model and published section data.  It neither
imports `elkies_rank28.py` nor reads a certificate derived from the exceptional
points.  It freezes its JSON output before the separate executable
[`verify_half_lattice_fake_descent_replay.sage`](../cas/verify_half_lattice_fake_descent_replay.sage)
loads the public fixture.  The search chooses its eleven independent candidates
by exact finite reduction without knowing the target exceptional basis.

The rank-29 and held-out experiments use the same separation.  Their input
builders are the only programs allowed to see the full public fixtures; they
write frozen inputs containing the curve and starting subgroup but omit the
held-out coordinates.  Search and fixture verification are separate programs.

## Reconstructed fake descent

Let

\[
 E:y^2=x^3+Ax+B,
 \qquad P=(x_P,y_P)\in E(\mathbf Q).
\]

Intersect (E) with the line

\[
 y=m(x-x_P)-y_P,
\]

which passes through (-P).  Removing the known intersection (x=x_P), the
other two intersections have a quadratic discriminant.  Writing its square as
(w^2) gives

\[
 C_P:\quad
 w^2=m^4-6x_Pm^2-8y_Pm-3x_P^2-4A.
\]

The exact map back is

\[
 x=\frac{m^2-x_P+w}{2},
 \qquad
 y=m(x-x_P)-y_P.
\]

The third-intersection law shows that a point of (C_P) gives a decomposition
(P=R+S).  Replacing (P) by (P+2T) changes coordinates within the same
class in (M/2M); this is why choosing a short representative *before*
building the quartic matters.

The deterministic pipeline is:

1. enumerate the parity class and compute a shortest representative;
2. form (C_P) from that representative;
3. clear denominators exactly;
4. call PARI `hyperellminimalmodel`, then `hyperellred`;
5. record modular square-sieve densities at fixed small good primes;
6. call `hyperellratpoints` at the declared bound;
7. invert every model change and map the point exactly to (E(\mathbf Q));
8. test membership and independence by exact group law and finite reduction.

For example, the first productive rank-28 class is `0x127e4`, with specialized
representative

\[
(0,0,-1,0,2,-1,1,-1,-1,1,1,0,0,-1,0,0,1).
\]

Its integral quartic has maximum coefficient size 762 bits; minimization and
reduction lower this to 95 bits, after which the bounded search recovers a new
quotient direction.  This is the missing ingredient in the previous naive
line-through-(-P) attempt: translation to a short coset representative and
model minimization/reduction, not a larger search on the original 500-bit
coordinate chart.

This formula matches the procedure Elkies describes under “Finding extra
rational points”: use the rank-17 lattice, search near deep holes in its
half-lattice, transform to quartics, and use `ratpoints`
([Elkies 2007, pp. 10--11](https://arxiv.org/pdf/0709.2908)).  PARI's model
changes and point-search interfaces are documented under
[`hyperellminimalmodel`, `hyperellred`, and `hyperellratpoints`](https://pari.math.u-bordeaux.fr/dochtml/html/Hyperelliptic_curves.html).

## Rank-28 classes

The table lists every searched half-class that produced a point with nonzero
image in the published exceptional quotient.  `q-mask` uses the verification-
only ordered exceptional basis (Q_1,\ldots,Q_{11}); it is not used to select
the class.  “gain” is incremental in specialized-depth search order.  Generic
depth 3 means one of the exact 43 deepest R17 classes.

| half-class | generic depth | specialized rank | specialized depth | q-mask(s) | gain | bits raw -> reduced |
|---|---:|---:|---:|---|---:|---:|
| `0x127e4` | 3 | 1 | 37.378 | `0x460` | 1 | 762 -> 95 |
| `0x05cf7` | 3 | 3 | 36.641 | `0x25c` | 1 | 710 -> 92 |
| `0x08692` | 3 | 7 | 36.245 | `0x008` | 1 | 704 -> 95 |
| `0x0de57` | 3 | 8 | 36.152 | `0x001,0x461` | 1 | 744 -> 94 |
| `0x050a0` | 3 | 9 | 36.075 | `0x461` | 0 | 746 -> 101 |
| `0x19312` | 2.5 | 17 | 35.541 | `0x461` | 0 | 725 -> 93 |
| `0x06226` | 3 | 19 | 35.440 | `0x040` | 1 | 729 -> 94 |
| `0x09765` | 2.5 | 23 | 35.208 | `0x29a` | 1 | 700 -> 93 |
| `0x048ea` | 2.5 | 24 | 35.192 | `0x471` | 1 | 689 -> 93 |
| `0x1655b` | 2.5 | 34 | 34.910 | `0x2ce` | 1 | 711 -> 92 |
| `0x18717` | 2.5 | 38 | 34.843 | `0x001` | 0 | 717 -> 93 |
| `0x18ad2` | 2.5 | 39 | 34.832 | `0x001,0x461` | 0 | 693 -> 93 |
| `0x04223` | 2.5 | 40 | 34.810 | `0x001` | 0 | 691 -> 93 |
| `0x0f687` | 2 | 43 | 34.744 | `0x020` | 1 | 707 -> 92 |
| `0x0b936` | 3 | 115 | 34.331 | `0x471` | 0 | 683 -> 93 |
| `0x18f2e` | 3 | 363 | 33.877 | `0x200` | 0 | 674 -> 93 |
| `0x01b43` | 3 | 1768 | 33.119 | `0x768` | 1 | 650 -> 93 |
| `0x1222e` | 3 | 2872 | 32.828 | `0x123` | 1 | 657 -> 94 |

For every labeled blind point, the verification artifact also stores its exact
relation in (M+\langle Q_1,\ldots,Q_{11}\rangle), quotient coordinates, and a
51-bit auxiliary-good-prime Kummer barcode.  These barcodes separate the known
28-dimensional group in the selected finite-prime sample; they are explicitly
not a complete global Selmer computation.

The productive centers span all 17 dimensions of (M/2M).  Their 153 pairwise
XORs have generic-depth histogram
(1:12, 3/2:47, 2:59, 5/2:35); none is itself deepest.  Thus there is no small
affine cluster of exceptional centers under XOR.

Ten of the 18 productive centers are exact generic deepest holes, but so are 33
of the 46 searched nonproductive centers.  Productive and nonproductive centers
have the same median reduced size (93 bits) and nearly the same median modular
square-density product (0.0333 versus 0.0307).  The hole lists concentrate the
successful search globally, but these elementary features do not rank the
successful classes well within that list.

## Experiment matrix

The following values are exact quotient ranks of the frozen blind candidate
sets, established by the separate verification stage.  “Hits” is the number of
distinct points outside the starting subgroup; it is a bounded search count.

| experiment | classes | median depth | median reduced bits | hits | recovered / target |
|---|---:|---:|---:|---:|---:|
| R17 rank 21, (t=3/8) | 43 generic | 15.17 | 56 | 5,927 | 4 / 4 |
| R17 rank 25, (t=-2/377) | 43 generic | 22.50 | 71 | 1,222 | 8 / 8 |
| R17 rank 26, (t=-308/251) | 43 generic | 25.11 | 75 | 649 | 9 / 9 |
| R17 rank 27, (t=2456/135) | 43 generic | 29.51 | 84 | 262 | 10 / 10 |
| R17 rank 28, (t=-9529/5471) | 43 generic | 34.96 | 93 | 74 | 9 / 11 |
| R17 rank 28, generic/specialized union | 64 | -- | 93 | 179 | 11 / 11 |
| 2024 rank 29 | 72 union | 36.65 | 103 | 94 | 10 / 12 |
| ICARM 356 rank 29 | 75 union | 33.79 | 93 | 239 | 12 / 12 |
| ICARM 385 rank 29 | 75 union | 36.10 | 105 | 77 | 4 / 12 |
| censored R17 (t=4) | 43 generic | 12.02 | 48 | 6,821 | gain 1 |
| censored R17 (t=6) | 43 generic | 13.24 | 52 | 3,974 | gain 1 |
| censored R17 (t=7) | 43 generic | 14.21 | 52 | 3,267 | gain 1 |
| censored R17 (t=8) | 43 generic | 11.44 | 55 | 5,819 | gain 1 |

The four ordinary controls are deliberately *censored*: the table certifies
only the blind gain found, not their full Mordell--Weil ranks.  Even so, it
exposes a serious scoring trap.  Raw quartic hit count is largest on these
controls and decreases along the record sequence while the quotient gain
increases.  Duplicate and generic-subgroup points must be eliminated exactly
before any fibre score is interpreted.

## Why many directions occur on one fibre

The replay recovers eleven independent directions through multiple distinct
half-lattice charts.  This agrees qualitatively with the earlier target-fitted
norm-8 trace-pencil computation, which found eleven distinct quartic
squareclasses through the eleven exceptional points
([R17_RANK28_GENUS_ONE_BISECTIONS_2026-09-02.md](../../elkies-k3/R17_RANK28_GENUS_ONE_BISECTIONS_2026-09-02.md)).
It also confirms that common trace/glue shape must not be identified with one
rank direction.

But the fake-descent quartics cannot supply the desired causal explanation:
their leading coefficient is the square 1, so they have rational points at
infinity and are birational to (E).  Their local-solubility survival fraction
is therefore 100% on every fibre.  The recorded modular conditions are only an
affine square-sieve cost model.

No replacement local signature emerged:

- pairwise XORs of productive centers show no deep-hole closure;
- modular square-density distributions overlap between productive and
  nonproductive classes;
- reduced quartic size does not separate them;
- the exact bad-place residual signature already recorded in
  [`BNF_FREE_RESIDUAL_2SELMER.md`](BNF_FREE_RESIDUAL_2SELMER.md) has rank 15 on
  the generic 17 points and remains rank 15 after all eleven exceptional
  directions are added.

Thus “several distinct coverings become soluble together” remains a plausible
description only for *genuine* 2-covering/Selmer classes that have not yet been
computed.  It is not a conclusion of this fake-descent replay.  A causal test
requires the global intersections and Cassels--Tate information scoped in
[`R17_RECORD_PAIR_HIGHER_2POWER_SELMER_PROGRAM.md`](R17_RECORD_PAIR_HIGHER_2POWER_SELMER_PROGRAM.md).

## Receptivity before point search

The tested pre-search profile contains top-depth counts and distributions,
generic/specialized overlap, reduced coefficient complexity, modular square
density, and their elementary correlations.  The ablation validates generic
depth as a *within-fibre search selector*, but the profile still does not pass
held-out *between-fibre rank-jump prediction*:

- local solubility is vacuous for these pointed charts;
- modular densities overlap and are not monotone in rank;
- coefficient size largely follows parameter height;
- generic/specialized overlap is 22 at rank 28 but only 11 on both ICARM
  rank-29 controls, despite their very different fixed-budget recoveries
  (12/12 and 4/12);
- the deepest-class point count is anticorrelated with quotient gain in the R17
  control matrix.

The operational conclusion is narrower: search generic deep holes first,
minimize and reduce every chart, and use the specialized list or union as a
recall stage when the first 43 do not saturate the intended control quotient.
Prioritize exact quotient gain rather than raw points.  Calling depth a
prospective rank-jump predictor would overstate the evidence.  This agrees
with the presentation-dependence warning in
[`R17_CARRIER_RECEPTIVITY_PROFILE_2026-09-04.md`](../../elkies-k3/R17_CARRIER_RECEPTIVITY_PROFILE_2026-09-04.md).

## Held-out curves 273 and 302

For each curve, the search was given several primitive starting subgroups of
dimensions 12--18 and not the remaining published points.  Exact recovery of
the displayed held-out quotient was measured only afterwards.  An interleaved
rank-17 subset checks sensitivity to point order.  Curve 245 is the known-family
adverse control.

| curve | starting subset | start dim | held-out dim | recovered |
|---|---|---:|---:|---:|
| 273 | prefix | 12 | 18 | 11 |
| 273 | prefix | 15 | 15 | 10 |
| 273 | prefix | 17 | 13 | 10 |
| 273 | prefix | 18 | 12 | 9 |
| 273 | interleaved | 17 | 13 | 9 |
| 302 | prefix | 12 | 19 | 4 |
| 302 | prefix | 15 | 16 | 4 |
| 302 | prefix | 17 | 14 | 3 |
| 302 | prefix | 18 | 13 | 4 |
| 302 | interleaved | 17 | 14 | 4 |
| 245 adverse control | prefix | 12 | 8 | 8 |
| 245 adverse control | prefix | 17 | 3 | 3 |
| 245 adverse control | prefix | 18 | 2 | 2 |

Curve 273 is consistently easier than curve 302 at this budget.  It is not
exceptional relative to the adverse control: curve 245 recovers every held-out
direction in all three tests, and produces thousands of candidates.  Therefore
the experiment detects no distinctive hidden R17/K3 provenance for 273 or 302.
The curve-245 blind finite-reduction field is intentionally null because its
candidate denominators obstructed every small prime in the blind certificate
routine; the displayed ranks come from exact post-search fixture verification.

## Reproduction and certificates

Use the Sage 10.9 Python recorded in the certificates.  The principal commands
are:

```text
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python elliptic-curves/cas/half_lattice_fake_descent_replay.sage
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python elliptic-curves/cas/verify_half_lattice_fake_descent_replay.sage
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python elliptic-curves/cas/replay_r17_generic_deep_holes_matrix.sage
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python elliptic-curves/cas/verify_r17_generic_deep_holes_matrix.sage
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python elliptic-curves/cas/replay_r17_rank21_half_lattice_control.sage
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python elliptic-curves/cas/verify_r17_rank21_half_lattice_control.sage
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python elliptic-curves/cas/replay_half_lattice_rank29_controls.sage
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python elliptic-curves/cas/verify_half_lattice_rank29_controls.sage
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python elliptic-curves/cas/replay_half_lattice_heldout_subgroups.sage
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python elliptic-curves/cas/verify_half_lattice_heldout_subgroups.sage
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python elliptic-curves/cas/replay_half_lattice_search_ablation.sage --phase development
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python elliptic-curves/cas/verify_half_lattice_search_ablation.sage --phase development
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python elliptic-curves/cas/replay_half_lattice_search_ablation.sage --phase holdout
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python elliptic-curves/cas/verify_half_lattice_search_ablation.sage --phase holdout
python3 elliptic-curves/cas/summarize_half_lattice_fake_descent_replay.py
python3 elliptic-curves/cas/summarize_half_lattice_search_ablation.py
```

The large blind JSON files retain every chosen representative, raw and reduced
quartic, inverse model change, modular sieve count, rational search coordinate,
exact mapped point, elapsed search time, and independence certificate.  The
verification JSON files add fixture-relative quotient masks only after the
blind boundary.  Input and executable SHA-256 hashes are embedded in every
certificate.

No change is made to `MATH_STATUS.json`: the exact result is a rank-lower-bound
replay of an already known fibre, while the proposed prospective mechanism and
Selmer-level simultaneity explanation remain unknown.

## Answers to the six questions

1. **Yes.**  The rank-28 quotient is blindly recovered in full from the
   specialized generic rank-17 subgroup.
2. It appears through the 18 productive half-classes listed above; eleven exact
   quotient directions are obtained, with repeated directions across several
   charts.
3. **Yes for bounded recovery, no as a rank predictor.**  The sealed holdouts
   give mean rank 8.33 for generic depth versus 3.20 for random sets, while all
   shallow sets give zero.  Depth does not order productivity perfectly:
   curve 385 has a random arm that beats either deep arm, and the union remains
   necessary for full rank-28 recall.
4. At present, no validated pre-point-search signature distinguishes all
   +8...+12 fibres.  Deep specialized heights track arithmetic complexity, but
   local survival, modular density, coefficient size, overlap, and raw hit count
   fail held-out calibration.
5. The new directions are found in several distinct charts, but these charts
   are pointed and already soluble.  Genuine simultaneous solubility of
   distinct 2-coverings is therefore **not established**.
6. **No provenance signal for 273/302.**  Curve 273 is easier than 302, but the
   adverse known-family control is easier still, so the experiment does not
   reveal a distinctive hidden parent.

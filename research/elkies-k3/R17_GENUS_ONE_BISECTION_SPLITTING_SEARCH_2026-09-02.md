# High-throughput R17 genus-one-bisection splitting search (2026-09-02)

<!-- status-consumer: EC-K3-R17-GENUS1-HIGH-THROUGHPUT-SPLITTING cad3d98ce58c89e7 -->
<!-- status-consumer: EC-K3-R17-NORM12-103B2-INTEGRAL-GLUE 52de13c8443f2b7d -->
<!-- status-consumer: EC-K3-R17-NORM12-103B2-MW-LATTICE-SIEVE aa0d8718eb57de6f -->
<!-- status-consumer: EC-K3-R17-NORM12-103B2-HARD-FIBRE-PRODUCT-H300000 b4fef7ab54b922e0 -->

## Result

The new search removes the rank-28 target fitting from the discovery loop.  It
uses 100 structurally diverse traces selected inside the 1,000
equation-cheapest norm-eight classes and all 43 norm-twelve deep classes.  All
143 compiled branch polynomials are irreducible squarefree quartics over
`QQ`, coprime to their trace denominators and to the degree-24 surface
discriminant.

The production run evaluated 5,474,328 primitive rational parameters: the
complete box

```text
|a| <= 2000,  1 <= b <= 2000,  gcd(a,b)=1
```

together with one million deterministic large-coordinate proposals.  A
separate population contains 900 parameters obtained from multiples on 100
pointed quartic Jacobians.  Twelve Legendre-symbol tables in three disjoint
prime blocks left 77,704 distinct extreme modular collisions.  Every one was
tested exactly, requiring 160,395 integer-square tests.

There is one exact simultaneous split:

```text
t = 1/25
norm-eight class  0x0f6b1
norm-twelve class 0x103b2
```

The norm-eight cover is the neutral member seeded by the sixteenth generic
R17 section at this parameter.  The norm-twelve point is

```text
x = 36075981547811164726251 / 244140625
y = -5898338731136062956741857359589376 / 3814697265625.
```

Exact finite reduction gives baseline rank 17 and combined rank 18 both in a
product of `E(F_p)/2E(F_p)` quotients and independently in a product of
`E(F_p)/3E(F_p)` quotients.  Thus the norm-twelve point is outside the
specialized generic MW17 subgroup.  This proves one quotient direction at
this fibre, hence a rank-at-least-18 specialization; it does not give a new
record rank.

The equation-free follow-up extracts the complete integral character tuple
inside the displayed cover-level rank-18 rational span.  It is the order-two
graph saturation of `R17(2)+<16>`, with determinant `497025024`, minimum
eight, and no roots.  Running each possible eclib saturation prime in an
isolated process proves that the displayed 18-point subgroup is primitive.
No specialization rank upper bound is claimed.

## Integral dissection of `0x103b2`

The census field that was previously unknown is now the exact tuple

```text
(L,A_L,q_L,G,{L_chi},{H_chi})
```

for the displayed cover-level rational span.  Here

```text
G=C2,
L=L_+ orthogonal-sum L_-=R17(2) orthogonal-sum <16>,
rank(L)=18,  det(L)=1988100096,
A_L=(Z/2)^16 + Z/8 + Z/3792,
H_chi=<(tau/2,T/2)> congruent to Z/2,
2R=tau+T.
```

The two projections of the glue generator have quadratic values `6` and `4`
in `QQ/2ZZ`, so their sum is isotropic.  Adjoining it has index two and gives
the full integral saturation inside the declared 18-dimensional rational
span.  The saturated visible lattice has determinant `497025024`, Smith group
`(Z/2)^14 + Z/8 + Z/3792`, minimum eight, no roots, and signed shell counts

```text
norm 8: 2622,   norm 10: 3058,   norm 12: 58254.
```

Those three counts form a unique fingerprint among the 143 selected covers
(126 fingerprints occur).  This is an integral distinction, not the cause of
the specialization: all 143 have the same basic `C2` character dimensions and
half-sum architecture.  At `t=1/25`, 141 other covers do not split.  The only
other split is norm-eight `0x0f6b1`, whose point is the specialized generic
section `P16`, so all 142 other covers fail to add a quotient direction.  The
third cover surviving all twelve local sieve primes, norm-eight `0x0c601`, is
an exact global nonsquare.

The trace representative of `0x103b2` has norm 12 and pinned coordinates

```text
(0,-1,1,1,-1,1,-1,0,1,1,0,1,1,-1,-1,0,-1).
```

Its branch and trace-pole quartics are irreducible with Galois group `S4`.
The four simple branch points avoid both the trace poles and the original
`24I1` discriminant.  Consequently the pulled surface has `48I1`, no
reducible-fibre component groups, and zero local component corrections.  The
exact intersections give `tau^2=24`, `T^2=16`, `R^2=10`, and `R.O=1`.

All eleven rank-28 exceptional directions instead share the single norm-eight
trace `0x1c6bc=-P2-P5`; the mask difference is `0xc50e`, of Hamming weight
seven.  They have the same `R17(2)+<16>` character lattice, the same order-two
graph architecture, and the same saturated determinant and discriminant
form, but not the same integral visible lattice.  Exact integral isometry
testing rejects an isometry, and their shell counts are
`(2638,2910,58886)`.  Their branch quartics likewise avoid the trace poles and
the original surface discriminant, so the local component correction is zero
on both sides; the target differs by having a rigid norm-twelve member and a
quartic `S4` trace-pole field rather than a fitted norm-eight pencil with a
quadratic trace denominator.  Their trace has square 16 after pullback, so
their lifts have height eight and `R.O=0`; the norm-twelve trace raises the
target lift to height ten.

At the specialization, the new point has numerical canonical height
`36.7815147269...` and orthogonal quotient defect `4.52376134659...`; its
shortest integral coset representative requires the nonzero generic correction

```text
(1,1,-2,2,1,-2,1,0,-2,-1,2,1,2,-2,2,-1,1)
```

and has height `22.1020672139...`.  By comparison, each rank-28 exceptional
point is already its own shortest coset representative, with raw height
`47.265917...` to `48.717476...` and quotient defect `16.530291...` to
`27.069127...`.  Thus `0x103b2` sits much closer to the specialized generic
real span, but in a genuinely shifted integral coset.  Exact eclib saturation
at `2,3,5,7,79813,239999` gives index one, proving the displayed specialized
rank-18 subgroup primitive.  This still supplies no rank upper bound and does
not exclude further anti-invariant directions on the cover.

The complete Gram matrices, finite quadratic forms, local polynomials, square
identity, and bridge-height data are stored in
[`elkies-k3-r17-norm12-103b2-mw-glue-v1.json`](../artifacts/generated-results/elkies-k3-r17-norm12-103b2-mw-glue-v1.json).

## Jacobian of the newly pointed `0x103b2` quartic

The point at `t=1/25` turns the norm-twelve branch cover into a pointed
genus-one curve.  After removing a rational square from its published branch
polynomial, an integral model is

```text
s^2 = 52558159476080896*t^4 + 63427547377764064*t^3
      - 88223393949768143*t^2 - 5393571474961890*t
      + 17005267967107473.
```

Its global minimal Jacobian is

```text
y^2 + x*y = x^3
  - 406976193745649728770658795455438*x
  + 3426347203723636438735122309348682280972198597892.
```

The rational torsion subgroup is trivial.  A deterministic
`hyperellratpoints` search of naive height at most `10000` finds 60 signed
affine quartic points.  Their 58 non-base images reduce to 17 independent
points, proving

```text
rank J_103b2(Q) >= 17.
```

eclib's saturation-index bound for this subgroup is

```text
137016286412 = 2^2 * 23 * 37 * 40251553,
```

with Tamagawa candidates `2,3,7,23`.  Isolated saturation at
`2,3,7,23,37` has index one.  The exact PARI 2-descent and saturation at the
remaining candidate `40251553` were stopped after long runs, so no rank upper
bound or claim that these 17 points generate all of `J_103b2(Q)` is made.

The useful change of strategy does not need those two unfinished bounds.
Direct Mordell--Weil lattice sieving evaluates every signed coefficient vector
against all 142 other covers modulo 32 good primes, then replays every local
survivor exactly.  Two complementary exhaustive shells are certified:

```text
coefficients in {-1,0,1}, support <= 8:  9,746,882 vectors
coefficients in {-2,-1,0,1,2}, support <= 5:  6,991,556 vectors
union after removing overlap:             16,496,324 vectors
```

The unit shell leaves four local survivors and the radius-two shell leaves
nine.  In each shell one point is exceptional for the affine inverse map and
one vector gives the known `t=1/25` collision with norm-eight `0x0f6b1`.
Every remaining survivor is an exact rational nonsquare on its locally
surviving cover.  Hence these 16,496,324 structured vectors produce no new
simultaneous split.  This is a bounded lattice-shell result, not a global
rational-point theorem for the pairwise fibre products.

The 17 independent points, exact survivor coordinates, prime histories, and
proof boundaries are stored in
[`elkies-k3-r17-norm12-103b2-mw-lattice-sieve-v1.json`](../artifacts/generated-results/elkies-k3-r17-norm12-103b2-mw-lattice-sieve-v1.json)
and
[`elkies-k3-r17-norm12-103b2-mw-lattice-unit-support8-v1.json`](../artifacts/generated-results/elkies-k3-r17-norm12-103b2-mw-lattice-unit-support8-v1.json).
Their SHA-256 digests are respectively
`1bed775cfd4e71b609c7a61066011b04845333250b3aa2da09fe46e0a2bbf1b0`
and
`d31b50706764c0a67b807bfa5c71e8a49b31a14fc4120abd890422ceae1929ab`.
The rank lower bound concerns the genus-one cover's Jacobian over `QQ`; it is
not a rank statement for the specialized K3 fibre or anti-invariant twist.

The seven covers that occur as false local survivors admit a faster second
filter.  If `0x103b2` and a partner both split at `t`, then the product of
their square coordinates gives a rational point on the degree-eight quotient

```text
z^2 = f_103b2(t) f_partner(t).
```

Each product is squarefree, so this quotient has genus three.  Exact PARI
searches through naive height `300000` find no affine rational point on any of
the seven quotients, and therefore no simultaneous split in that range.  The
certificate is
[`elkies-k3-r17-norm12-103b2-hard-fibre-products-v1.json`](../artifacts/generated-results/elkies-k3-r17-norm12-103b2-hard-fibre-products-v1.json).
This is a bounded necessary-condition search, not a proof that any quotient
has no rational points globally.

## Trace geometry

For a norm-eight trace

```text
tau = (Nx/h^2, Ny/h^3),  deg(h)=2,
```

regular chord slopes form the pencil

```text
M = M0 + lambda*h^2.
```

The search initializes one member without any exceptional point.  It chooses
a distinct small rational fibre and a generic R17 section on that fibre,
solves the one linear incidence equation for `lambda`, and rejects any member
that is reducible, singular, or badly branched.  The resulting rational point
only initializes the pointed binary quartic.  The actual auxiliary parameter
population is then obtained from group-law multiples on its Jacobian using
the exact inverse pointed-quartic map.  Seed parameters are forced distinct,
so repeated generic initialization cannot create an artificial simultaneous
collision.

For every norm-twelve deep trace the exact closest representative has

```text
deg(h)=4,  deg(M0)<8,
q=(M0^4-6*M0^2*Nx-8*M0*Ny-3*Nx^2-4*A*h^4)/h^6,
deg(q)=4.
```

There is no free `lambda` at this degree: `M0` is the unique regular member.
Complete exact construction verifies all 43 quartics.  The bounded
small-height point pass did not point these curves in advance; the hit above
was found by the common parameter sieve.

## Sieve architecture

The old mixed-trace scanner stored exactly 128 cover bits.  The replacement
uses `ceil(number_of_covers/64)` words throughout, so increasing
`--norm8-count` from 100 through 1,000 requires no source change.  Each prime
table stores, for every `(a:b)` modulo `p`, the simultaneous quadratic-residue
mask for every branch quartic.  A prime at which a clearing denominator
vanishes gives that curve no filtering information; prime selection minimizes
this blindness while balancing the `p^2` table cost.

For each parameter the scanner maintains separate intersections in at least
three disjoint prime blocks.  Its ranking key is, lexicographically,

```text
minimum blockwise F2-rank of surviving trace masks,
minimum blockwise surviving-cover count,
sum of blockwise ranks,
sum of blockwise counts,
all-block survivor count.
```

This distinguishes many locally correlated covers from covers spanning many
classes in `R17/2R17`.  It is still only a local diversity heuristic.  For an
exact hit, the lifted formulas construct every point on the specialized
Weierstrass equation, and the common finite-quotient implementation computes
the image of those points modulo the specialized generic MW17 subgroup.

## Literature boundary

Garbagnati--Salgado explain why special multisections are geometrically tied
to rank jumps, but do not supply this simultaneous arithmetic sieve or its
diversity score:

- A. Garbagnati and C. Salgado,
  [*Rank jumps and Multisections of elliptic fibrations on K3 surfaces*](https://arxiv.org/abs/2505.15159).

The binary-quartic/Jacobian step is the standard degree-two genus-one covering
construction.  Magma's genus-one-model handbook gives an independent
software-level description of the same invariant/covariant map:

- [Magma handbook: genus-one models as coverings](https://magma.maths.usyd.edu.au/magma/handbook/text/1594).

These references justify the geometric and covering framework.  They do not
turn a bounded Legendre scan into a completeness theorem.

## Reproduction

Run the stored production profile with:

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/search_elkies_2026_genus_one_bisection_splitting.sage \
  --norm8-count 100 \
  --equation-pool-size 1000 \
  --output artifacts/generated-results/elkies-k3-r17-genus-one-bisection-splitting-search-v1.json \
  --local-directory artifacts/local/elkies-k3/r17-genus-one-bisection-splitting/production-v1
```

The generated artifact has SHA-256
`4745f53993675286173298c9444da023825fd4429e74ede582a1d8c14979d07e`.
The C++ scanner is compiled by the Sage driver.  Passing `--exact-limit 0`
(the default) exact-tests every modular extreme; a positive limit is an
explicit truncated development run and is recorded as such.

Compute the rank-17 subgroup and run the radius-two/support-five lattice shell
with:

```bash
sage -python elkies-k3/scripts/search_r17_norm12_103b2_mw_lattice.sage
```

The complementary unit/support-eight shell is:

```bash
sage -python elkies-k3/scripts/search_r17_norm12_103b2_mw_lattice.sage \
  --max-support 8 --coefficient-radius 1 \
  --output artifacts/generated-results/\
elkies-k3-r17-norm12-103b2-mw-lattice-unit-support8-v1.json
```

Search the seven difficult pairwise genus-three quotients with:

```bash
sage -python \
  elkies-k3/scripts/search_r17_norm12_103b2_hard_fibre_products.sage \
  --height 300000
```

## Boundary

The compact box is exhaustive, but the million large-coordinate proposals,
the trace selection, and the Jacobian multiple ranges are bounded.  The
Legendre ranks are search rankings, not Selmer bounds or Mordell--Weil ranks.
Positive finite-quotient escape proves non-membership in the generic subgroup;
failure to escape the chosen quotient products does not prove dependence.
Nothing here supplies the fifteen quotient directions required for rank 32.

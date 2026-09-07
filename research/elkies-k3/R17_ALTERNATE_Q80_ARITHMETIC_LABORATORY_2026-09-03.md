# Alternate-Q80 arithmetic laboratory (2026-09-03)

<!-- status-consumer: EC-K3-R17-NORM12-11952-INHERITED-COVERS eb333139c25202b2 -->
<!-- status-consumer: EC-K3-R17-NORM12-11952-INHERITED-PRODUCTS b0d303c94466c6f9 -->
<!-- status-consumer: EC-K3-R17-NORM12-11952-ALTERNATE-LAB-1024 c2f6309f8d6cc06d -->

## Result

The direct `norm12-orbit-11952` alternate-Q80 equation is now the first
arithmetic laboratory for bisection characters.  The promoted exact prefix
consists of:

- all 121 published-R17 height-four curves that become degree-two curves over
  the alternate base;
- all 7,260 products of their quadratic characters;
- the complete exact equation-cost ordering of the 39,147 native alternate
  bisection classes; and
- exact equations for only the cheapest 1,024 native classes.

The inherited and native records use different stored parity bases.  After
transport to the direct compiled alternate frame modulo `2M`, the 1,145 raw
records contain two overlaps:

```text
inherited-115d4  = alternate-orbit-11e87  modulo section translation,
inherited-18e8e  = alternate-orbit-11a86  modulo section translation/inversion.
```

In both cases the primitive branch divisor is identical, the displayed
squareclass ratio is `1/4`, and the direct-frame vectors differ by an element
of `2M`.  These are duplicate realizations of one translation class, not
rank-two squareclass collisions.  The combined laboratory therefore has
exactly 1,143 distinct classes.

## The 121 inherited covers

For the alternate fibre class `D=(3,2,w)`, an old height-four section `S_v`
has degree

```text
D.S_v = 5 - <w,v>.
```

Complete exact enumeration gives 121 oriented sections of degree two.  For
each curve, restriction of the compiled alternate coordinate `u=L1/L0` gives
a degree-two relation in the old parameter.  Its discriminant is normalized
in `QQ(u)^*/QQ(u)^*2` while retaining the rational constant squareclass.  All
121 canonical branch polynomials are smooth quadratics coprime to the
alternate `24 I1` discriminant, and the transported point satisfies the
alternate Weierstrass equation over `QQ(u,s)`, `s^2=q(u)`, coefficientwise.

The exact artifact is
[`../artifacts/generated-results/elkies-k3-r17-norm12-11952-inherited-bisection-covers-v1.json`](../artifacts/generated-results/elkies-k3-r17-norm12-11952-inherited-bisection-covers-v1.json),
with SHA-256

```text
6ba1db53ed02de18910456e75dcdad55e87205c7355986080dab1955d83050b9
```

## The 7,260 inherited pair products

All `binomial(121,2)=7,260` products are compared exactly in
`QQ(u)^*/QQ(u)^*2`, including rational constants.  The 121 individual
characters are distinct, all 7,260 pair products are distinct, no product is
a third inherited character, and none formally matches the older comparison
catalogue.  The older rank-28 and `q_103b2` characters live over the published
base; exact Neron--Severi intersections give degrees 8 and 9 over the
alternate base, so a formal coefficient rename cannot supply an alternate
character.

The exact artifact is
[`../artifacts/generated-results/elkies-k3-r17-norm12-11952-inherited-product-characters-v1.json`](../artifacts/generated-results/elkies-k3-r17-norm12-11952-inherited-product-characters-v1.json),
with SHA-256

```text
fe97c92d5ea8609dd31f2e68a9cfe157ee93530d6e0b1e35f4714e59b0d6409d
```

## Cheapest 1,024 native classes

The complete 39,147-class priority table minimizes, in order, the exact group
addition upper bound, support, maximum coefficient, coefficient `L1` norm,
and the section-basis word.  The bounded compiler then applies the unique
regular residual chord to ranks 1 through 1,024 only.  Every trace uses the
finite chart, every resulting branch quadratic is squarefree and avoids the
24 singular fibres, and every lifted section passes both coefficient
identities.

The promoted prefix and its compact squareclass certificate are:

- [`../artifacts/generated-results/elkies-k3-r17-norm12-11952-alternate-bisections-cheapest-1024-v1.json`](../artifacts/generated-results/elkies-k3-r17-norm12-11952-alternate-bisections-cheapest-1024-v1.json),
  SHA-256 `6d58c129e0309c63e0b8421837192f2b0e18f6ac331b5608d74c1288566be7fa`;
- [`../artifacts/generated-results/elkies-k3-r17-norm12-11952-alternate-bisection-collisions-cheapest-1024-v1.json`](../artifacts/generated-results/elkies-k3-r17-norm12-11952-alternate-bisection-collisions-cheapest-1024-v1.json),
  SHA-256 `73ee65a628bef1fdf1d03fbfb7bf027941a543b68343ae028763a8b49a3e2756`.

The native prefix alone has 1,024 distinct branch squareclasses.  After the
two inherited/native overlaps are removed, all 1,143 primitive branch
quadratics are irreducible over `QQ` and pairwise distinct.  Hence distinct
classes share no geometric branch fibre.

## Pair bases, product twists, and triples

Exact Hasse--Minkowski calculations find a rational point on all 1,143
individual conics.  Thus each gives a rational base change with generic rank
at least 18.

Every pair of distinct classes has four disjoint branch points.  Consequently
all

```text
binomial(1143,2) = 652653
```

pair bases are connected genus-one `V4` covers.  Their two new directions
occupy distinct characters and have height matrix

```text
[24  0]
[ 0 24],
```

so the generic rank is at least 19 over each paired base.  This statement does
not assert that every genus-one torsor has a rational point or determine a
base Mordell--Weil rank.

Unique factorization and disjoint branch support make all 652,653 pair-product
squareclasses distinct.  None equals a compiled degree-two section character.
Thus this prefix contains no triple of the three nontrivial characters of one
`V4` group with a section in every character.  A section on an uncatalogued
quartic product twist has not been excluded; that calculation remains
`UNKNOWN`.

On the other hand, every three distinct laboratory characters are independent
over `F_2`.  Hence all

```text
binomial(1143,3) = 248225691
```

triples give connected `(Z/2Z)^3` covers branched at six points.  Their genus
is five, their three pulled directions have height matrix
`diag(48,48,48)`, and character decomposition proves generic rank at least
`17+3=20`.  This is a rank-20 construction over a genus-five base, not the
desired rank-20 construction over a rational or low-rank `V4` base.

The complete bounded classification is
[`../artifacts/generated-results/elkies-k3-r17-norm12-11952-alternate-arithmetic-laboratory-cheapest-1024-v1.json`](../artifacts/generated-results/elkies-k3-r17-norm12-11952-alternate-arithmetic-laboratory-cheapest-1024-v1.json),
with SHA-256

```text
540801bd087a7be7585b9387bea8c2d3c6a5df16405629841e386441cd44217f
```

## Foundry score and boundary

The arithmetic incidence score for this exact prefix is

```text
unique classes                         1143
shared branch-fibre pairs                 0
genus-zero paired bases                   0
catalogued product-character sections     0
three-character V4 closures                0
```

Thus raw bisection abundance and low equation cost do not correlate with
controlled branch incidence in this prefix.  The next useful search is not a
blind extension of the prefix.  It should prioritize either an exact quartic
product-twist section calculation or a foundry frame whose branch-incidence
score is positive.

## Replay

```bash
sage -python \
  elkies-k3/scripts/construct_r17_norm12_11952_inherited_bisections.sage --check

.venv/bin/python \
  elkies-k3/scripts/analyze_r17_norm12_11952_inherited_products.py --check

sage -python \
  elkies-k3/scripts/rank_r17_norm12_11952_alternate_bisection_orbits.sage --check

sage -python \
  elkies-k3/scripts/construct_r17_norm12_11952_alternate_bisections.sage \
  --start 0 --limit 1024 \
  --output artifacts/generated-results/elkies-k3-r17-norm12-11952-alternate-bisections-cheapest-1024-v1.json \
  --check

.venv/bin/python elkies-k3/scripts/hash_bisection_extensions.py \
  --compact \
  --input artifacts/generated-results/elkies-k3-r17-norm12-11952-alternate-bisections-cheapest-1024-v1.json \
  --output artifacts/generated-results/elkies-k3-r17-norm12-11952-alternate-bisection-collisions-cheapest-1024-v1.json

sage -python \
  elkies-k3/scripts/analyze_r17_norm12_11952_alternate_laboratory.sage --check
```

The promoted scope ends at native priority rank 1,024.  The remaining 38,123
native classes and all uncatalogued quartic product-twist Mordell--Weil groups
remain outside this certificate.

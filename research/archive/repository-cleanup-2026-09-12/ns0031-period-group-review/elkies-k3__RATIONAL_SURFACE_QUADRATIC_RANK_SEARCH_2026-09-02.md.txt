# Rational elliptic surfaces and quadratic rank decompositions — 2026-09-02

<!-- status-consumer: EC-K3-RES-QBC-E6A1-RHO19 7103fa2a1a4e7ba2 -->
<!-- status-consumer: EC-K3-NS0031-MARKED-FORMAL-BRANCH b31e99bce4edac0a -->

## Outcome

There is a new one-dimensional rational family of quadratic base changes with
an exact invariant/anti-invariant Mordell--Weil decomposition.  Start from the
rank-one rational elliptic surface

```text
E_c / QQ(u):
y^2 = x^3 + (u-3)x + c^2*u^2 + u - 2,
P = (-1,c*u).
```

For `k in QQ` put

```text
D      = 3*k^2-4,
c      = 2*k/D,
lambda = -(k^2-4)*D/4,
d(u)   = u*(u-lambda).
```

Then the quadratic twist in the convention

```text
d(u)*y^2 = x^3 + (u-3)x + c^2*u^2 + u - 2
```

has the exact section

```text
x = 4*u^2/D^2 + k^2*u/D - 1,
y = 8*u^2/D^3 + 2*(k^2+2)*u/D^2.
```

The cover `u=lambda/(1-t^2)` has squareclass `d(u)`, and its K3 pullback has

```text
fibres             2IV* + I4 + 4I1,
root lattice       2E6 + A3, rank 15,
MW height Gram     diag(1/3,3),
MW rank            2 = 1 invariant + 1 anti-invariant,
Picard rank        19 generically,
abs(disc NS)       36.
```

All equations, sections, the bounded ansatz elimination, heights, and rank
budget are replayed by
[`scripts/certify_rational_surface_quadratic_rank_search.sage`](scripts/certify_rational_surface_quadratic_rank_search.sage).
Its generated certificate is
[`../artifacts/generated-results/elkies-k3-rational-surface-quadratic-rank-search-v1.json`](../artifacts/generated-results/elkies-k3-rational-surface-quadratic-rank-search-v1.json).
Its SHA-256 hash is
`888c7bba059af4cff18a8baa359b60f227575a4835d28b7764f62c42690f8d8c`.
The integral lattice, saturation, transcendental lattice, special fibres, and
first neighbor layer are analysed separately in
[`E6A1_RHO19_K3_DISSECTION_2026-09-02.md`](E6A1_RHO19_K3_DISSECTION_2026-09-02.md).

This is a new equation family in the repository.  It is not being attributed
to the literature cited below.  The underlying `E6+A1` rational-surface chart
is the one-modulus family used by Kimura; the rational section subfamily, the
quadratic branch locus, and the anti-invariant section above are the new exact
construction.

## What “rank 17--19” means here

For this search, `17--19` is the generic **Picard rank** of the K3, not the MW
rank of its displayed base-change fibration.  Shioda--Tate gives

```text
rho(K3) = 2 + root rank + MW rank.
```

A one-dimensional non-isotrivial K3 family cannot have generic Picard rank
20.  Thus a generic MW-rank-17 K3 must be rootless with `rho=19`, as in the
published R17/Golay target.  A rational-surface base-change presentation can
instead store most of the same divisor budget in reducible fibres and display
only one or two sections.  This is exactly the Golay/NS0031 source pattern:
the low-MW source fibration is a construction chart, while a same-surface
neighbour can expose MW rank 17.

For a rational elliptic surface with root rank `R` and MW rank `8-R`, an
unramified quadratic pullback has the visible lower bound

```text
rho >= 2 + 2R + (8-R) = 10+R.
```

Branching at a multiplicative `I_n` fibre replaces the two unramified copies
of `A_(n-1)` by `A_(2n-1)` and adds one more root.  For the selected
`E6+A1` seed, branching at its rational `I2` fibre gives root rank `15` on the
K3.  The invariant generator gives `rho>=18`; one twist section gives
`rho>=19`.  Hence rank sum two is the largest possible generic rank sum in
this fixed root-rank-15 stratum.  It is not claimed to be globally large among
all rational-surface base changes.

## Low-complexity rational-surface seed catalogue

The bounded catalogue prioritizes one-dimensional rank-seven root strata,
because their generic quadratic pullbacks naturally land in the Picard
`17--19` window.

| seed | equation over `QQ(u)` | generic fibres | free MW lattice | rational section status |
|---|---|---|---|---|
| `E7` | `y^2=x^3+(u-3)x+a*u-2` | `III*+3I1` | `[1/2]` | on `a=-c^2-2`, `(-a,c*(a-1))` |
| `D7` | `y^2=x^3+u^2*(u^2-3)x+u^3*(a*u^3+u^2-2)` | `I3*+3I1` | `[1/4]` | geometric rank-one seed; no compact `QQ` generator compiled here |
| `E6+A1` | `y^2=x^3+(u-3)x+c^2*u^2+u-2` | `IV*+I2+2I1` | `[1/6]` | `(-1,c*u)`; selected |
| Golay control | explicit quotient in the prior audit | `I6+I3+3I1` | `[1/2]`, torsion `Z/3` | exact generator and exact anti-invariant direction |
| NS0031 control | explicit quotient over `GF(7)` | `I8+4I1` | rank one by Shioda--Tate | finite-field trace and anti-trace only |

The first three equations are the rank-seven rational-surface charts in
[Kimura’s quadratic-base-change construction](https://arxiv.org/abs/1802.05195).
The lattice entries follow from the rational-surface Shioda--Tate and
discriminant formulas and agree with the classification in
[Oguiso--Shioda](https://rikkyo.repo.nii.ac.jp/records/10006).  The catalogue is
not the full Oguiso--Shioda table.  Herfurtner’s
[four-singular-fibre classification](https://archive.mpim-bonn.mpg.de/id/eprint/860/)
is the natural next source for semistable and mixed-fibre additions.

The higher-dimensional Kimura charts `E6`, `D6`, `D5`, and `A4` have MW
lattices `A2*`, `A1*+A1*`, `A3*`, and `A4*`, respectively, but they were not
promoted into this first arithmetic shortlist: their raw moduli dimensions
are two through four, and the printed `D6/D5` normalizations use
`QQ(sqrt(3))`.  They are the next place to search for rank sums three through
five after imposing rational section-first slices.

## Exact twist and K3 equations

Write `H=1-t^2`.  The quadratic map

```text
u = lambda/H
```

has branch values `0` and `lambda`, because

```text
u*(u-lambda) = (lambda*t/H)^2.
```

After the homogeneous pullback, the K3 equation is

```text
Y^2 = X^3
    + H^3*(lambda-3*H)*X
    + H^4*(c^2*lambda^2+lambda*H-2*H^2).
```

The invariant section is

```text
P0 = (-H^2, c*lambda*H^2).
```

If

```text
p = 4/D^2,
q = k^2/D,
s = 8/D^3,
v = 2*(k^2+2)/D^2,
```

then the anti-invariant section is

```text
P1 = (
  p*lambda^2 + q*lambda*H - H^2,
  lambda*t*(s*lambda^2 + v*lambda*H)
).
```

Literal substitution proves both identities.  Under the deck involution
`t -> -t`, `P0` is fixed and `P1` is negated.  The character decomposition
therefore gives

```text
rank E(QQ(k)(t))
  = rank E(QQ(k)(u)) + rank E^(d)(QQ(k)(u)).
```

The rational surface has root lattice `E6+A1`, trivial torsion, and generator
height `1/6`, so its rank is exactly one.  At `k=1,u=1`, the twist specializes
to

```text
y^2 = x^3 - 49/8*x + 1029/64,
Q = (7/2,-49/8),
```

and exact torsion computation gives infinite order.  Thus the generic twist
section is non-torsion.

The K3 fibres are `IV*` at `t=1,-1`, `I4` at infinity, and four simple `I1`
fibres.  Clearing the component groups by multiplication by twelve and using
compact pole degrees gives

```text
<P0,P0> = 1/3,
<P1,P1> = 3,
<P0,P1> = 0.
```

The cross term vanishes by deck-character orthogonality.  These sections and
the root lattice give `rho>=19`.  If `alpha,beta` are the two residual `I1`
values on the rational surface, a marked configuration invariant of their
four preimages is

```text
Q(lambda)/Q(0)
  = (2*k^2+1)*(9*k^2-20)^2 / (9*(3*k^2+4)^2),
```

which is nonconstant.  Hence the lattice-polarized K3 family is
non-isotrivial.  A non-isotrivial one-dimensional K3 family with nineteen
independent divisor classes has generic Picard rank exactly nineteen.

The integral divisor Gram has discriminant group `Z/3+Z/12`.  Its only
possible proper even overlattices have index three.  At `(k,t)=(1,3)`, none
of `P0`, `P1`, `P0+P1`, or `P0-P1` is three-divisible on the exact good fibre, so the
displayed MW lattice is saturated and `abs(disc NS)=36` is exact.  The full
argument and Gram matrix are in the dissection note linked above.
Shioda--Tate now makes the MW rank exactly two and the twist rank exactly one.

For rational `k`, the clean arithmetic open is especially simple:

```text
k not in {-2,0,2}.
```

The other apparent denominator `3*k^2-4` has no rational zero.  The residual
`I1` discriminant and the branch-smoothness value are

```text
16*(9*k^4+12*k^2+16)^3/(3*k^2-4)^6,
-(2*k^2+1)*(9*k^2-20)^2/(3*k^2-4)^2,
```

so no additional rational exclusion is needed for the generic fibre profile.
Individual rational specializations can still have Picard or MW rank jumps.

## Complete bounded section ansatz

The discovery calculation solved

```text
d(u)*y(u)^2 = x(u)^3 + (u-3)*x(u) + c^2*u^2 + u - 2
```

with

```text
x=p*u^2+q*u+r,
y=s*u^2+v*u+w,
p=z^2, s=z^3,
```

and nonzero leading coefficient and `lambda`.  The constant equation is

```text
(r-2)*(r+1)^2=0.
```

After triangularly eliminating `v,w`, the two resultants are, up to nonzero
scalar and degenerate factors,

```text
r=-1:
(3*lambda*z^2-4*z+1)^2*(3*lambda*z^2+4*z+1)^2,

r=2:
27*lambda^3*z^2-72*lambda^2*z^2+64.
```

The `r=-1` branch gives the clean family above.  The `r=2` branch also has a
rational parameterization, obtained from

```text
3*(h^2-m^2)=16,
lambda=(8-h^2)/3,
c=m/(h^2-8),
h=(3*k^2+16)/(6*k),
m=(16-3*k^2)/(6*k).
```

It again supplies one invariant and one anti-invariant direction, but its
section coefficients grow much faster.  The clean `r=-1` component therefore
dominates it on equation size, parameter height, and arithmetic searchability.

This elimination is complete only in the declared polynomial degree box.  It
does not exclude rational-function sections, higher polynomial degrees, a
second twist direction on another family, or higher rank sums from lower-root
rational surfaces.

## Recovery of the Golay and NS0031 controls

The prior exact audit is retained unchanged in
[`RATIONAL_SURFACE_BASE_CHANGE_AUDIT_2026-09-02.md`](RATIONAL_SURFACE_BASE_CHANGE_AUDIT_2026-09-02.md).

- The rational Golay-chart `3I6+6I1` model is the pullback of an
  `I6+I3+3I1` rational surface under `u=t+1/t`.  Its free rank splits exactly
  as one invariant generator of height `1/2` plus one anti-invariant twist
  direction.  Hidden `3`-torsion and a half-section make its saturated NS
  determinant `20`, so it remains rejected from the determinant-720 target.
- NS0031 model 157 over `GF(7)` is the pullback of an `I8+4I1` rational
  surface under `u=t^2/(t-1)`.  One marked pair gives a nonzero invariant
  trace and a nonzero anti-invariant direction.  This trace calculation is
  finite-field; a later certificate proves a one-parameter formally smooth
  `ZZ_7` marked branch, but no `QQ` source equation or full MW decomposition
  is claimed here.

These controls validate the character split and repeated-fibre diagnostics.
The new `E6+A1` family differs structurally: it is constructed over `QQ` from
a forced twist section, not recognized after an involution happened to appear
on a previously found K3.

## Search priorities

The current Pareto order is:

1. the clean `E6+A1` family above for exact Picard-19 arithmetic experiments;
2. rank-two rational-surface seeds (`E6`, then a rationalized `D6`) with a
   branch at a rational simple singular fibre, targeting rank sum at least
   three and Picard rank `18--19`;
3. semistable rank-one Herfurtner strata, where multiplicative reduction is
   more convenient for specialization sieves;
4. higher-degree or denominator-layer twist sections only after the
   degree-`(2,2)` positive control is used to calibrate the search.

The quadratic-base-change method and injection of rational-surface sections
are standard; see [Kimura](https://arxiv.org/abs/1802.05195).  Explicit
high-geometric-rank K3 constructions and their field-of-definition issues are
treated by [Kumar--Kuwata](https://arxiv.org/abs/1409.2931) and in their
[rank-15 construction](https://arxiv.org/abs/1604.00738).  The distinction is
important here: a geometric rank over `QQbar(t)` is not an arithmetic rank
over `QQ(t)`.  Elkies's current
[rank-17 equation and quadratic-base-change paper](https://arxiv.org/abs/2608.25406)
is the direct literature control for the repository's R17 bisection work; its
rank-18 covers start from a K3, whereas the construction here starts from a
rational elliptic surface.  The current family supplies both displayed
directions over the rational function field itself.

## Reproduction and status boundary

Generate and byte-check the certificate with

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_rational_surface_quadratic_rank_search.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_rational_surface_quadratic_rank_search.sage --check
```

The theorem package proves the displayed generic family and the bounded
ansatz result.  It does not classify all rational elliptic surfaces, prove a
global maximum for the rank sum, or prove that every rational `k`
specialization has exact rank two.  Specializations must be minimized and
their points and independence certified again before arithmetic rank claims
are made.

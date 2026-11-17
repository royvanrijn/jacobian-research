# Rational plane sections of the fixed quintic Keller map

The fixed map has two rational target planes on which every rational fibre
has a rational point. Every other fixed rational affine target plane has only
`O(B^2)` rational-point fibres in projective height `B`. This distinguishes
the plane contributions to ambient counts; it does not determine the full
ambient leading constant.

Retain the normalized map `Phi=L o F0 o A` and target coordinates `(u,b,c)`
from [the fixed-map manuscript](../../unpublished/fixed-map-hasse-failures/sections/02-fixed-map.tex).
Put `pi=-u`. On `pi!=0`, the entire fibre is the derivative-localized root
scheme of

\[
 E_{\pi,b,c}(S)=\pi^5S^5-\frac32\pi^4S^4+\frac32\pi S^3
 +\left(-\frac54\pi-\frac9{32}b\right)S^2
 +\frac9{16}S-\frac9{32}c.                              \tag{1}
\]

A generic rational section of a plane means a source point over its rational
function field. All planes below are defined over `Q`. The height is always
`H([1:u:b:c])`.

## Proposition PS1: two planes with rational sections at every rational target

Every rational target in the planes `u=0` and `c=0` has a rational source
point. On `c=0` there is the polynomial section

\[
 (u,b,0)\longmapsto
 \left(0,\ b-\frac{40}9u,
       -u-\frac32\left(b-\frac{40}9u\right)^2\right).    \tag{2}
\]

**Proof.** For source `x=0`, the source normalization `A` leaves the point
unchanged, and `t=1`, `q=z+(3/2)y^2`. The fixed-map coordinates reduce to

\[
 \Phi(0,y,z)=(-q,y-40q/9,0).
\]

Substitution gives (2), including `u=0`.

On `u=0`, use the boundary computation in
[Lemma FMA1](FIXED_MAP_AMBIENT_LOCAL_CLASSIFICATION.md). An explicit section
on `b!=0` is

\[
 x=\frac4b,\qquad y=-\frac b2,\qquad
 z=\frac{3b^3}{64}\left(
 \frac{10}b+\frac{16}3\left(\frac{3b}{16}\right)^4
              -\frac{16}3\left(\frac{3b}{16}\right)^5-c\right).
                                                               \tag{3}
\]

Indeed the unnormalized source coordinate after `A` is `x/2=2/b`, so `t=0`,
`q=3b^2/32`, and direct evaluation gives `B=b,C=c,Pi=0`. For the remaining
line `u=b=0`, the source point `(c,0,0)` maps to `(0,0,c)`. Thus the assertion
holds at every rational target, not just on the generic open. \(\square\)

## Theorem PS2: all other planes have irreducible generic quintic and only O(B^2) global points

The planes of Proposition PS1 are the only rational affine target planes
admitting a generic rational section. If `P` is any other fixed rational
affine plane, the polynomial (1) is geometrically irreducible of degree five
over its function field, and

\[
 \#\{y\in P(\mathbb Q):H(y)\le B,
                         \Phi^{-1}(y)(\mathbb Q)\ne\varnothing\}
       \ll_P B^2.                                             \tag{4}
\]

**Proof of generic irreducibility.** Work over either `Q` or its algebraic
closure; the same argument applies to both constant fields.

First take a plane `c=alpha pi+beta b+delta`, so `pi,b` are independent
coordinates. Over `K=Qbar(pi)`, write its inverse equation as

\[
 N_\pi(S)-\frac9{32}b(S^2+\beta),\qquad
 N_\pi(S)=\pi^5S^5-\frac32\pi^4S^4+\frac32\pi S^3
             -\frac54\pi S^2+\frac9{16}S
             -\frac9{32}(\alpha\pi+\delta).                    \tag{5}
\]

If a root `rho` of `S^2+beta` is nonzero, it is constant and `N_pi(rho)`
is a polynomial in `pi` with nonzero leading coefficient `rho^5`. It cannot
vanish identically. If `rho=0`, then `beta=0`, and `N_pi(0)` vanishes
identically precisely when `alpha=delta=0`. This last case is exactly the
plane `c=0`. Hence on every other plane of this form,

\[
 \gcd_{K[S]}(N_\pi(S),S^2+\beta)=1.                            \tag{6}
\]

The polynomial (5), regarded as linear in `b`, is primitive over the UFD
`K[S]` and irreducible over `K(S)`. Gauss's lemma makes it irreducible in
`K[S,b]`. Its leading coefficient in `S` is the unit `pi^5`, so it is also
primitive over `K[b]`; a second use of Gauss's lemma gives irreducibility
of degree five over `K(b)`.

For a plane with equation `b=alpha pi+beta`, its coordinates are `pi,c`.
The inverse equation is `N_pi(S)-(9/32)c`, whose coefficient of `c` is a
nonzero constant. The same primitive linear-parameter argument gives
irreducibility of degree five. The remaining planes have `pi=a` constant.
For `a!=0`, work over `Qbar(b)` and again use the scalar coefficient of
`c`. The case `a=0` is the boundary plane `u=0`, already treated in PS1.
These cases exhaust the affine planes. On `pi!=0`, a source point yields a
simple rational root of (1) by the exact reconstruction theorem. Thus an
irreducible generic quintic precludes a generic section.

**Proof of the height bound.** For a plane other than the two sections,
`pi Disc(E)` is a nonzero rational function on the plane: generic
irreducibility in characteristic zero implies generic separability. After
removing its zero locus, the monic equation `E/pi^5` gives a geometrically
integral finite etale cover of degree five of a dense open subset of `P`.
Every rational-point target on this open is the image of a rational point
of that cover, by derivative-unit reconstruction.

Apply Buggenhout--Cluckers--Salberger--Santens--Vermeulen,
[arXiv:2506.13471v2, Theorem 1.1, PDF p. 4](https://arxiv.org/pdf/2506.13471v2),
through the projective closure `Pbar=P^2`. Its hypotheses are an integral
source and a quasi-finite dominant morphism of degree at least five to
projective space. The root cover followed by the open immersion meets all
of them. It gives `O_P(B^2)` image points. A fixed rational linear
isomorphism `P^2 -> Pbar subset P^3` compares the standard height in both
directions by fixed constants, so the exponent is unchanged.

The deleted locus is a proper algebraic curve. In fixed plane coordinates,
its primitive integer representatives lie on a nonzero homogeneous
polynomial in three variables; the elementary polynomial zero bound gives
`O_P(B^2)` such triples. This absorbs all boundary points, including the
intersection with `pi=0`, and proves (4). \(\square\)

The cover theorem is used in its projective form. Its separate affine
leading-form hypotheses are not assumed. The source need not be proper.
The two section planes each have order `B^3` rational-point targets and
therefore exhaust the fixed affine planes with that order. There is no
uniform summation over moving planes or nonlinear surfaces in this argument.

## Corollary PS3: doubled lower coefficient for the ambient local/global counts

Let `G(B)` and `L(B)` count ambient rational targets with a rational source
point and with points over every completion, respectively. Then

\[
 G(B)\ge\frac8{\zeta(3)}B^3+O(B^2),\qquad
 L(B)\ge\frac8{\zeta(3)}B^3+O(B^2).                            \tag{7}
\]

Together with [Theorem ABC3](FIXED_MAP_AMBIENT_COUNT_BOUNDS.md), both counts
have sharp order `B^3`. No ambient leading-constant asymptotic is asserted.

**Proof.** Each coordinate plane has primitive representatives
`[q:m:n]`, `1<=q<=B`, `|m|,|n|<=B`, and `gcd(q,m,n)=1`.
Their exact number is

\[
 P_2(B)=\sum_{d\le B}\mu(d)\lfloor B/d\rfloor
                         (2\lfloor B/d\rfloor+1)^2
       =\frac4{\zeta(3)}B^3+O(B^2).                            \tag{8}
\]

The planes intersect in the rational line `u=c=0`, whose count is
`O(B^2)`. Their union therefore contains
`2P_2(B)-O(B^2)=8B^3/zeta(3)+O(B^2)` rational targets, all globally
soluble by PS1. Global points imply points over every completion, proving
the same lower bound for `L`. \(\square\)

## Exact root-incidence parametrization and remaining scope

For completeness, on `pi!=0` all globally soluble targets are the images of
rational triples `(pi,b,s)` under

\[
 (\pi,b,s)\longmapsto\left(-\pi,b,
   \frac{32}9\pi^5s^5-\frac{16}3\pi^4s^4
    +\frac{16}3\pi s^3-\frac{40}9\pi s^2+2s-bs^2\right),       \tag{9}
\]

restricted to `E'(s)!=0` at that target. Add the plane `u=0` to recover all
rationally soluble targets. This is an image parametrization with possible
multiple rational roots over one target, not an injective target count.
The plane theorem and (9) do not prove that the union of the two planes
accounts for the complete ambient `B^3` coefficient.

## Verification

Run from the repository root with Python 3 and SymPy:

```sh
.venv/bin/python research/arithmetic/verify_fixed_map_plane_sections.py --write-receipt
```

The [checker](verify_fixed_map_plane_sections.py) verifies the two sections,
the generic-plane resultant obstruction and its degenerate case, the scalar
linear-parameter identities, the root-incidence formula, and bounded exact
primitive plane-union counts. The
[receipt](fixed_map_plane_sections_checks.json) pins its source and bounds.
The degree arguments and imported projective theorem are written proof
inputs; bounded controls do not verify the unrestricted plane classification
or the asymptotic. The primary PDF is retained under `sources/` and was read
at the applied theorem. This is contextual internal work, with no formal
verification, independent whole-proof replay, or external review claimed.

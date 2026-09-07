# Curve 302: rational common covers and the branch-first incidence curve

Authorities: `EC-K3-CURVE302-SHARED-BRANCH-GATE`,
`EC-K3-CURVE302-RATIONAL-BRANCH-CARRIER` and
`EC-K3-CURVE302-BRANCH-CARRIER-ARITHMETIC` and
`EC-K3-CURVE302-BRANCH-JACOBIAN-ENLARGEMENT`. The
[parent objective](../elliptic-curves/notes/CURVE302_CONSTRUCTION_RECOVERY.md)
remains open. No alternative parent or generic MW basis is constructed here.

## Why compare one shared branch point?

Let two quadratic covers of the same rational base have distinct two-point
geometric branch sets. Their normalized fibre product has geometric degree
four and Galois group `V4`. If the branch sets share exactly one point, the
union has three points, each with inertia two. Thus

\[
2g-2=-8+3\cdot2=-2.
\]

If both covers have rational unramified lifts of the parameter producing
302, the product has a rational point and this genus-zero normalization is
`P1` over `Q`. Pullback can then put both sections on one surface over a
rational function field. Independence would follow if their specializations,
together with the source sections, were certified independent. A full generic
basis would still require separate rank and saturation proofs.

For disjoint branch sets the same calculation gives
`2g-2=-8+4*2=0`, hence genus one. No nonconstant map from `P1` to that
normalization exists in characteristic zero. This also excludes a common
rational cover of higher degree for those two fixed maps. These are direct
applications of [Riemann–Hurwitz](https://stacks.math.columbia.edu/tag/0C1B).

The checker includes an exact positive control with branch polynomials
`u(u-1)` and `u(u-4)`. Put

\[
r=\frac{k^2-6k+12}{k^2-12},\quad w=kr+3-k,\quad
u=\frac{4r^2-1}{r^2-1},\quad v_2=\frac{w}{r^2-1},\quad v_1=rv_2.
\]

It verifies `v1²=u(u-1)`, `v2²=u(u-4)` and `deg(u(k))=4`.
At `k=1` the point is `(-25/24,35/24,-55/24)`. This is a control of
base-curve geometry, not an elliptic-surface rank calculation.

## Exact result for the 6,417 stored covers

The [certificate](../artifacts/generated-results/elkies-k3-curve302-shared-branch-gate-v1.json)
checks all 31 pinned cubic-pencil checkpoints against their source hashes.
It reconstructs every discriminant from its two residual line quadratics,
verifies the rational unramified lift at `u=0`, and matches every normalized
branch polynomial to the original branch fingerprint.

Every monic quadratic `u²+b*u+c` is irreducible over `Q`. For each one the
certificate records a prime at which the reduced numerator times denominator
of `b²-4c` is a nonsquare. No large integer factorization is used.
Within each pencil the 207 polynomials are distinct. Two distinct irreducible
quadratics over `Q` cannot share an algebraic root. Hence all

\[
31\binom{207}{2}=660951
\]

within-pencil pairs have disjoint branch sets and genus-one fibre products.
There are no rational common covers for these pairs. This does not compare
different source pencils or exclude other pointed lines, multisections or
parents. The older identical-cover grouping alone did not establish this
stronger obstruction.

## All pointed lines for one source: a genus-three gate

The [incidence certificate](../artifacts/generated-results/elkies-k3-curve302-rational-branch-carrier-v1.json)
uses pencil zero, anchored by public points 1–8, and its first basepoint `P`.
It retains the complete symbolic line

\[
X=x_P+s,\qquad Y=y_P+ms.
\]

Dividing the restrictions of `F0,F1` by the known basepoint factor `s` gives
two quadratics `g0(s,m),g1(s,m)`. The branch polynomial for their cover is

\[
\operatorname{disc}_s(g_0+u g_1)=d_a(m)u^2+d_b(m)u+D_4(m).
\]

The checker proves, coefficient by coefficient,

\[
d_b^2-4d_aD_4=16R_8,\qquad R_8=\operatorname{Res}_s(g_0,g_1),
\]

where `R8` is squarefree of degree eight. Its eight rational roots are exactly
the slopes from `P` to the other eight pencil basepoints. Rational branch
points therefore require a rational point on `Y²=R8(m)`, a genus-three curve.
The certificate stores its exact coefficients and all eight rational roots.

A rational residual point on the 302 fibre additionally requires

\[
H^2=D_4(m)=m^4-6x_Pm^2+8y_Pm-3x_P^2-4a,
\]

where `y²=x³+a*x+b` is the stored short model of 302. Its point is recovered
by `s=(-g0[1](m)+H)/(2*g0[2](m))`. Thus this genus-one double cover of the
slope line is birational over `Q` to 302 itself.

The polynomials `R8,D4` are squarefree and coprime. Their normalized `V4`
compositum has twelve branch points and genus `12-3=9`. Equivalently, it is
a degree-two cover of 302 branched at the sixteen distinct points over the
eight roots of `R8`. The checker supplies eight rational representatives
with `Y=0`; their opposite `H` signs give the other eight. These correspond
to degenerate lines through two basepoints and are not new cover candidates.

There is no omitted rational exception: `da` and `D4` both have no projective
root modulo 17, and the vertical line has an irreducible branch quadratic.
Outside the eight degenerate slopes, the simultaneous square conditions are
an exact description of rational pointed lines with rational branch points
and rational residual points on 302. Their global rational solutions remain
**UNKNOWN**. The construction certificate did not run a point search; the
subsequent finite probe is described below.

## Consequence for the next construction attempt

The degree-two map `C9 -> E302` induces pullback and norm on Jacobians with
composition `[2]`. Its pullback has finite kernel, so the certified rank-31
subgroup (`ECR31`) gives `rank J(C9)(Q) >= 31`. The usual Chabauty condition
`rank J(C9)(Q) < genus(C9)=9` cannot hold. Merely compiling this curve does
not make a general genus-nine rational-point search a justified next step.

The genus-three curve `Y²=R8(m)` is a smaller necessary gate whose Jacobian
rank has not been determined. Its arithmetic, followed by the `D4` square
test, is a concrete branch-first route. A pair of resulting covers would
still need a shared branch value, independent section images and a complete
surface construction. Other construction patterns remain eligible; none of
these arguments imposes a generic rank on an alternative parent.

## Arithmetic follow-up: torsion, simplicity and a finite point probe

Authority: `EC-K3-CURVE302-BRANCH-CARRIER-ARITHMETIC`. The
[frozen protocol](../artifacts/generated-results/elkies-k3-curve302-branch-carrier-protocol-v1.json)
compares all 336 ordered choices of three branch points sent to infinity,
zero and one. Square content is removed by trial division at primes below
10,000 followed by a perfect-square check; no general integer factorization
is run. The square-scaled original polynomial has coefficient height 1,318
bits. Chart selection keeps the smallest coefficient-height presentation per
unordered triple, then the best two distinct triples, with heights 1,578
and 1,609 bits. These charts put three branch points at small coordinates;
their coefficients are larger than the square-scaled original's, so no
globally minimal model claim is made.

The control `y²=40320*product(x-i,i=0..7)` has the nondegenerate point
`(8,40320)`. The identical selection and search procedure recovers it in both
chosen control charts at height 64, with exact transports. This verifies the
search and transport machinery; it does not calibrate point heights on the
target curve.

The target bound is `|numerator(z)|<=512`, `1<=denominator(z)<=512`, in
lowest terms. Each chart contains **319,407** reduced rational values.
PARI's results agree with a complete integer-pair replay using eight good
prime sieves and exact integer square tests. The only finite points are the
degenerate points `z=0,1`; each chart's point at infinity is also a known
degenerate branch point. The **638,814** count includes overlap between
charts. There are no nondegenerate genus-three points or new lifts to 302 in
these boxes. This is a finite exclusion only, with no automatic larger search.

The [probe certificate](../artifacts/generated-results/elkies-k3-curve302-branch-carrier-probe-v1.json)
also computes Frobenius at two good primes. Every polynomial is independently
recovered by direct counts over `Fp`, `Fp²` and `Fp³` and Newton identities:

| Prime | Three curve point counts | Jacobian order |
|---:|---|---:|
| 47 | 40, 2,332, 104,248 | 90,176 |
| 53 | 52, 2,988, 148,804 | 147,968 |

Their order gcd is **64**. The eight rational branch points already supply
all `2^(2*3)=64` geometric two-torsion points over `Q`. Prime-to-residue-
characteristic injectivity of torsion under good reduction gives the upper
bound; the possible 47- and 53-primary parts are excluded at the other prime
(`147968 mod47=12`, `90176 mod53=23`). Consequently

\[
J(C_3)(\mathbb Q)_{\rm tors}\cong(\mathbb Z/2)^6.
\]

The Frobenius polynomial at 47 is

\[
T^6-8T^5+93T^4-432T^3+4371T^2-17672T+103823.
\]

Modulo 3 this is `T^6+T^5+T+2`, which is irreducible. The
[Rabin certificate](../artifacts/generated-results/elkies-k3-curve302-branch-jacobian-simple-v1.json)
checks `T^(3^6)=T` modulo this polynomial and supplies Bezout identities for
the gcds with `T^(3^2)-T` and `T^(3^3)-T`. Thus the Frobenius polynomial is
irreducible over `Q`. Any proper positive-dimensional abelian subvariety
defined over `Q` would factor that good-reduction polynomial. Therefore
**the genus-three Jacobian is simple over `Q`**. In particular it has no
elliptic quotient over `Q`, including no factor isogenous to 302. No absolute
simplicity claim or restriction on factors over larger number fields is made.

The arithmetic follow-up initially suggested the following conditional gate.
If its rank were
zero, every rational divisor class would have order dividing two. Choosing
a rational Weierstrass point `W`, a rational point `P` would then satisfy
`2P~2W`. The unique hyperelliptic degree-two linear system forces `P` to be
Weierstrass. Hence a rank-zero proof would close all nondegenerate rational
points in this fixed branch-first problem, not just the two tested boxes.
Conversely, a non-Weierstrass rational point would give a nontorsion divisor
class. **The subsequent construction below disproves the rank-zero premise.**
This conditional argument remains valid but can no longer close this curve.

## A natural divisor proves positive rank

Authority: `EC-K3-CURVE302-BRANCH-JACOBIAN-ENLARGEMENT`. The original
identity `16R8=db²-4daD4` supplies the rational degree-four divisor

```
D4(m)=0,  Y=db(m)/4
```

on `Y²=R8(m)`. Subtract four times the rational Weierstrass point used as
infinity in the first frozen odd-degree chart. This defines a rational
Jacobian point `D`; it does not require an additional rational curve point.

The [generator certificate](../artifacts/generated-results/elkies-k3-curve302-branch-jacobian-generator-v1.json)
transports this divisor and reduces it by Cantor's algorithm. Its reduced
Mumford pair has degrees `(3,2)`. Since the second polynomial is nonzero,
the unique reduced representation differs from its negative. Thus `D` is
not two-torsion. The entire rational torsion group has exponent two, so
**D is nontorsion and rank J(C3)(Q) is at least one**. An independent
finite check also finds `2D != 0` at 47. The rank-zero stopping argument is
therefore unavailable, rather than merely unproved.

For a positive control of Jacobian-to-curve recognition, the checker uses
the non-Weierstrass point `(8,40320)` on the earlier control curve. The
complete finite word-image procedure retains its class at 11. Applying
that procedure to `nD+T`, all 64 rational two-torsion classes `T` and
`|n|<=512`, excludes 65,568 words at 47 and another 24 at 53. The remaining
eight are exactly the Weierstrass classes. This tests a subgroup, not the
whole Jacobian or every rational curve point.

## Recovering a rational half and proving two-saturation

The natural point is not primitive at two. Let `F(z)` be the odd-degree
curve polynomial, `c` its leading coefficient, and `(U,V)` the reduced pair
of `D`. All seven nonzero values `(-c)^3 U(r_i)` at the finite branch roots
are rational squares. Their square roots interpolate a polynomial `S` with
`S²=(-c)^3 U mod F`.

The [enlargement certificate](../artifacts/generated-results/elkies-k3-curve302-branch-jacobian-enlargement-v1.json)
implements [Stoll's halving equations, Proposition 5.1](https://www.mathe2.uni-bayreuth.de/stoll/papers/ratpts-selmer-2016-12-07.pdf):

```
v = w*S mod F,    v = u*V mod U,
deg(u)<=1, deg(v)<=4, deg(w)<=3.
```

The linear system has a one-dimensional kernel. After normalizing `w` to
be monic, it verifies `u²F=v²-(-c)^3 U w²` and `gcd(u,w)=1`, giving the half
`[w,-v/u mod w]`. Exact Jacobian arithmetic verifies its double is `D`.
All 64 torsion translates are checked, and the one with smallest coefficient
bit height is retained as `H_J`. Every one doubles to `D`; none is represented
by a single curve point. This proves the subgroup enlargement

\[
[\langle H_J,J(\mathbb Q)[2]\rangle:
  \langle D,J(\mathbb Q)[2]\rangle]=2.
\]

It also proves the enlarged subgroup **two-saturated in J(Q)**. At 53,
`#J(F53)=147968=2^9*17²` and the two-torsion rank is six. The largest
two-primary invariant-factor exponent is therefore at most `9-5=4`.
The reduction of `H_J` has order `4624=16*17²` and attains that bound.
Every odd multiple of `H_J` plus a two-torsion class still has two-primary
order 16; twice any point over F53 has two-primary order at most eight.
Hence no such odd word has a rational half. If `2X=2nH_J+T`, a nonzero
`T` would give rational four-torsion, already excluded. For `T=0`,
`X-nH_J` is rational two-torsion and belongs to the subgroup. This covers
all possible first index-two enlargements without assuming exact rank one.

The updated finite sieve uses `|n|<=1024` and all 64 torsion translates.
The embedding `nD+T=(2n)H_J+T` preserves every old tested word; every old
finite membership is replayed under it. Of **131,136** words, 131,080 are
excluded at 47 and 48 more at 53. Only the eight Weierstrass classes remain.
These are finite subgroup exclusions, with no unrestricted multiple search.

We now have an explicit **two-saturated rank-one subgroup**, not a full
Mordell-Weil basis of the Jacobian. Its exact rank, other independent
directions, odd-index enlargements and global non-Weierstrass points remain
unknown. No new cover or alternative elliptic parent has been constructed.

## Replay

```
sage -python elkies-k3/scripts/certify_curve302_shared_branch_gate.sage --check
sage -python elkies-k3/scripts/construct_curve302_rational_branch_carrier.sage --check
sage -python elkies-k3/scripts/probe_curve302_branch_carrier.sage --check
sage -python elkies-k3/scripts/certify_curve302_branch_jacobian_simple.sage --check
sage -python elkies-k3/scripts/certify_curve302_branch_jacobian_generator.sage --check
sage -python elkies-k3/scripts/enlarge_curve302_branch_jacobian.sage --check
```

All six checks pass. Their inputs include the retained local checkpoints under
`artifacts/local/elkies-k3/curve302-cubic-pencil-overlap-v1/`. If these are not
present, regenerate them with the original pinned
`search_curve302_cubic_pencil_overlap.sage` script before replay. The new
checks do not enumerate new Mordell-Weil point combinations. The finite probe
replays PARI's point search by integer enumeration and its Frobenius output
by extension-field point counts. This provides independent computational
mechanisms within the checker, not a fully independent reconstruction of
the source pencil or all the geometric arguments.

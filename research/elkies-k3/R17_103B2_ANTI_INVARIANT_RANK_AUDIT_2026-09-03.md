# Anti-invariant rank audit for the R17 cover `0x103b2` — 2026-09-03

<!-- status-consumer: EC-K3-R17-NORM12-103B2-ISOTROPIC-FRAME 47f3a0eb7ee50bcb -->
<!-- status-consumer: EC-K3-R17-POINTED-COVER-JACOBIAN-CONTROL-H10000 4bb087b3a1ebc684 -->

## Outcome

The actual anti-invariant rank is **not yet determined**. The exact result is

\[
1\leq \operatorname{rank}E^{(q_{103b2})}(\mathbb Q(t))\leq 22.
\]

The lower bound is represented by an exact polynomial section of height eight.
The upper bound is the geometric Hodge/Shioda--Tate bound, not a sharp Picard
certificate. Consequently the double cover still has only the proved rank
lower bound `17+1=18`; this audit does not license a `+2` mutation claim.

No Magma calculation is used. All new exact arithmetic uses SageMath,
Singular/msolve, C++, and rational reconstruction.

## The hidden isotropic frame

The lattice test behind the auxiliary genus-one curve is now exact.  In the
pinned decomposition

```text
NS = U + R17(-1),
```

the norm-twelve trace representative

```text
w=(0,-1,1,1,-1,1,-1,0,1,1,0,1,1,-1,-1,0,-1)
```

gives the bisection class

```text
D_103b2=(3,2,w).
```

It satisfies

```text
D^2=0,  D.F=2,  D.O=1,  div(D)=1.
```

Splitting off an exact integral hyperbolic plane containing `D` leaves a
positive rank-17 frame of determinant `948`.  Complete norm-two enumeration
finds no roots.  PARI integral isometry testing identifies this frame with the
published R17 `J2` class and rejects the certified alternate Q80 rootless
frame.  Thus the bisection does hide a second rootless rank-17 fibration at
the lattice level, but it is a new marking in the published frame class, not
the alternate Q80 class.

This gives a compelling geometric explanation for why the individual
pointed-cover Jacobian exposes seventeen rational directions: the quartic is
a member of the isotropic pencil `|D|`.  The certificate deliberately stops
at the exact lattice/J2 statement.  It does not construct the second
Weierstrass equation and its sections over `QQ`, identify its `J1` orbit under
surface automorphisms, or improve the anti-invariant twist-rank upper bound.
The replay and full determinant-one marking are in
[`elkies-k3-r17-norm12-103b2-isotropic-frame-v1.json`](../artifacts/generated-results/elkies-k3-r17-norm12-103b2-isotropic-frame-v1.json).
Its SHA-256 digest is
`efbbebf23bc6608ce83cc0f0cb74947ceeaea6647af8e27343e18b24e4975b95`.

## Ten-cover point-search control

A deterministic uniform sample of ten covers from the 100 other pointed
norm-eight members was run through the same bounded pipeline:

```text
PARI hyperellratpoints(H=10000)
  -> exact pointed-quartic map
  -> eclib Mordell--Weil relation reduction.
```

The seed was `0x103b2`, and the selected production indices were
`4,17,19,25,27,30,55,67,75,96`.  Every control returned exactly the two signed
points above its pointing parameter.  Those two points lie in the exceptional
fibre of the affine point map, so every control had zero nonbase mapped points
and zero rank visible from this search.  Each control still has its separately
certified nontorsion pointed-Jacobian generator, hence a known rank lower bound
of one.  The identical positive regression for `0x103b2` returned 60 signed
points, 58 nonbase images, and relation rank 17.

Therefore ranks `8,12,17` are not common *at this search visibility* in the
sample: `0x103b2` is a sharp `17 versus 0` nonbase outlier.  This remains a
bounded point-height experiment, not an upper bound for any control Jacobian;
their missing independent points may simply be taller.  The selected labels,
denominator-cleared integral models, point counts, and exact rank-reduction
records are in
[`elkies-k3-pointed-cover-jacobian-rank-controls-h10000-v1.json`](../artifacts/generated-results/elkies-k3-pointed-cover-jacobian-rank-controls-h10000-v1.json).
Its SHA-256 digest is
`4b84801a6d599b2402f1808a21bb18d9d45ef7fdc4a0cf219b1be3ba7c4f92ed`.

## Exact twist and known section

Let

\[
E^{(q)}:\quad Y^2=X^3+Aq^2X+Bq^3,
\]

where `A,B` are the compact published R17 coefficients and `q=q_103b2` is
the certified quartic branch polynomial. If `P` is the split bisection on
`s^2=q(t)`, then

\[
T=P-\sigma(P),\qquad
(X,Y)=\bigl(qx(T),q^2[y(T)]_s\bigr)
\]

descends exactly to the twist. Literal rational-function arithmetic gives
`deg X=8`, `deg Y=12`, with denominator one, and the Weierstrass residual is
zero. The cover height is 16, so the twist height is 8. At `p=37` this point
is a reduced isolated point of its fixed-leading polynomial-section block
(Jacobian rank `8/8`). The full coordinates are stored in
[`elkies-k3-norm12-orbit-103b2-twist-section-v1.json`](../artifacts/generated-results/elkies-k3-norm12-orbit-103b2-twist-section-v1.json).

The twist has arithmetic genus `chi=4`. Its four new branch fibres are `I0*`,
while the 24 original nodal fibres remain irreducible. Thus the trivial
lattice has rank `2+4*4=18`. Since `h^{1,1}=10 chi=40`, Shioda--Tate gives
geometric Mordell--Weil rank at most `40-18=22`. Proving rank one is therefore
equivalent to proving Picard rank 19 for this `chi=4` surface.

## Exhaustive `P.O=0` shells

A section disjoint from the zero section has `deg X<=8`, `deg Y<=12`. After
fixing its point on a smooth infinity fibre with nonzero leading `Y`, the high
coefficients recover `Y` recursively from `X`. The remaining twelve
coefficient equations define the complete fixed-leading polynomial-section
scheme.

The C++ enumerator first imposes all fibre-value quadratic-residue conditions.
For small primes it interpolates from eight values; for larger primes it uses
a four-plus-four meet-in-the-middle bitset intersection. Every survivor is
then checked against the literal polynomial identity.

| prime | coverage | interpolants represented | unsigned `X` solutions | result |
|---:|---|---:|---:|---|
| 17 | complete | 227,448,000 | 3 | known point plus two characteristic-specific points |
| 19 | complete | 655,855,200 | 1 | exactly the known signed pair |
| 29 | complete | 79,781,760,000 | 276 | seven full-tangent-rank branches |

The exact `p=19` and `p=29` ledgers are
[`p19-polynomial-section-bruteforce-v1.json`](../artifacts/generated-results/elkies-k3-norm12-orbit-103b2-p19-polynomial-section-bruteforce-v1.json)
and
[`p29-polynomial-section-bruteforce-v1.json`](../artifacts/generated-results/elkies-k3-norm12-orbit-103b2-p29-polynomial-section-bruteforce-v1.json).

The `p=29` full shell restores the leading `X` and `Y` coefficients, giving
13 equations in 10 variables. Seven modular points have tangent rank 10.
High-precision Newton lifting shows:

- the known branch converges through precision `29^800` and reconstructs to
  the exact known rational section;
- the other six reduced branches fail unused equations at valuations one,
  two, or three and do not define `29`-adic sections;
- the remaining 269 singular branches all survive the first lift to `29^2`,
  so they are not eliminated by this calculation.

See
[`p29-hensel-lifts-v1.json`](../artifacts/generated-results/elkies-k3-norm12-orbit-103b2-p29-hensel-lifts-v1.json).

This is a bounded shell theorem only. A second Mordell--Weil generator can
have positive intersection with `O`, and even another rational `P.O=0`
section can have nonintegral reduction at a selected prime. Neither the unique
`F_19` point nor the reduced-branch `29`-adic audit is a global rank upper
bound.

## Heuristic and failed upper-bound routes

A 48-prime Nagao-style pass over `211--491` initially gave block scores
`1.186--1.798` (mean `1.500`), unusually high compared with the earlier
singleton census. Six independent eight-prime blocks over `503--1811` fell
to `0.501--1.569` (mean `1.115`). The signal is therefore compatible with
rank one and is not stable evidence for rank two.

The fixed-leading `p=37` msolve block containing the known point exceeded a
five-minute, roughly 2.9 GB bound. This is a solver limit, not evidence for
another component.

A direct 2-descent does not collapse to rational factor arithmetic: the
2-torsion cubic defines a degree-three cover with genus ten under the 24
simple discriminant branches, so its `S`-class and Jacobian contribution must
be controlled. SageMath and PARI do not currently provide that general
function-field descent.

The remaining clean upper-bound route is a Frobenius polynomial for the
`chi=4` twist surface. Its conductor degree is 32, hence the nontrivial
elliptic `L`-polynomial has degree 28. The open-source controlled-reduction
implementation described by Costa--Harvey--Kedlaya applies to smooth ordinary
projective hypersurfaces, whereas this model is a degree-24 hypersurface in
`P(1,1,8,12)`; a toric/weighted adaptation is required. See
[Costa--Harvey--Kedlaya](https://arxiv.org/abs/1806.00368). Computing the twist
`L`-function only modulo 2 is also too weak here because quadratic twist
characters disappear modulo 2 and the original rank-17 central factor
persists; compare [Boudreau](https://arxiv.org/abs/2110.12156).

## Next decisive computation

The shortest rigorous closeout is one of:

1. implement weighted/toric controlled reduction at two good primes and use
   Shioda--Tate (plus Artin--Tate discriminant classes if both reductions have
   Picard rank 20); or
2. complete a function-field 2-Selmer computation in the genus-ten cubic
   algebra and show the Selmer dimension is one.

Extending polynomial shells alone cannot prove the rank upper bound without
one of these global inputs.

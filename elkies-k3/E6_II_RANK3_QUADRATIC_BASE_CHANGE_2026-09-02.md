# E6 rank-sum-three quadratic base change — 2026-09-02

<!-- status-consumer: EC-K3-RES-QBC-E6-II-RANK3-RHO19 5b10608e230145e9 -->

<!-- status-consumer: EC-K3-RES-QBC-E6-II-Q2-MW4 3aa5084463780acc -->

<!-- status-consumer: EC-K3-RES-D6-RATIONALIZED-SECTION-CHART a94042dd2d76797c -->

## Outcome

The first target is attained by a rational one-parameter family.  Over
`QQ(r)(u)`, put

```text
a =  2*(r^2+r+1)/(r+1),
c = -2*r^2/(r+1),

E_r: y^2 = x^3 + a*u*x + u*(u+c).
```

On the rational open

```text
r not in {-2,-1,-1/2,0,1},
```

the generic fibre profile and two sections are

```text
fibres: IV* + II + 2I1,
P = (1,   u+1),
Q = (r^2, u+r^3).
```

Their height matrix is

```text
[ 2/3  -1/3 ]
[ -1/3  2/3 ],
```

so they generate a rank-two subgroup.  The `E6` root rank is six, and
Shioda--Tate on a rational elliptic surface bounds the geometric MW rank by
two.  Hence the rank over `QQ(r)(u)` is exactly two.

Take the twist squareclass

```text
d(u)=u*(u+c).
```

The twist in the convention

```text
d(u)*y^2=x^3+a*u*x+u*(u+c)
```

has the section `(x,y)=(0,1)`.  Thus

```text
rank E_r(QQ(r)(u)) + rank E_r^(d)(QQ(r)(u)) = 2+1 = 3.
```

The exact replay is
[`scripts/certify_e6_ii_rank3_quadratic_base_change.sage`](scripts/certify_e6_ii_rank3_quadratic_base_change.sage),
with generated certificate
[`../artifacts/generated-results/elkies-k3-e6-ii-rank3-quadratic-base-change-v1.json`](../artifacts/generated-results/elkies-k3-e6-ii-rank3-quadratic-base-change-v1.json).

## K3 pullback and exact Picard rank

Write `H=1-t^2` and use the quadratic map

```text
u=-c/H,
u*(u+c)=(c*t/H)^2.
```

After clearing the standard weights, the K3 equation is

```text
Y^2 = X^3 - a*c*H^3*X + c^2*t^2*H^4.
```

It has the three explicit sections

```text
P = (H^2,     H^2*(H-c)),             invariant,
Q = (r^2*H^2, H^2*(r^3*H-c)),         invariant,
S = (0,       c*t*H^2),                anti-invariant.
```

The cover branches at the type `II` fibre `u=0` and the smooth fibre
`u=-c`.  The K3 fibre profile is

```text
2IV* + IV + 4I1,
```

with root lattice `2E6+A2`, of rank 14.  The MW height matrix is

```text
[ 4/3  -2/3   0  ]
[ -2/3  4/3   0  ]
[ 0      0    2/3 ],
```

of determinant `8/9`.  Consequently `rho>=2+14+3=19`.  The marked fibre
configuration varies with

```text
16*(r^2+r+1)^3 / (27*r^2*(r+1)^2),
```

so the family is nonconstant.  A nonconstant one-dimensional family of K3
surfaces cannot have generic Picard rank 20; therefore its generic Picard
rank is exactly 19 and its displayed MW rank is exactly three.  The character
decomposition makes the rank split exactly `2+1`.

The integral `NS` matrix in the certificate has

```text
abs(disc NS)=24,
Smith invariants 2,2,6.
```

At `(r,t)=(2,2)`, all seven nonzero mod-two combinations of `P,Q,S` fail the
divisibility-by-two test.  Since 24 permits only an index-two saturation
defect, this proves saturation.  Generic torsion is 3-primary from the fibre
component groups, while the same good control fibre has torsion of order two;
specialization therefore eliminates generic torsion.

## Same-NS rootless MW17 search

This success has an exact negative answer to the requested same-NS search,
stronger than a bounded neighbour census.  Rootlessness makes all fibre
corrections in the torsion height formula vanish, so such a fibration has no
MW torsion.  A rootless Jacobian fibration on a Picard-rank-19 K3 would
therefore have an even positive-definite MW frame of rank 17, minimum at least
four, and determinant `abs(disc NS)=24`.  Its Hermite invariant would be

```text
4 / 24^(1/17) = 3.3179596... .
```

Blichfeldt's bound gives

```text
gamma_17 <= (2/pi)*Gamma(2+17/2)^(2/17) = 3.2821242... .
```

The required invariant is larger than the upper bound.  Equivalently, every
rootless even rank-17 lattice must have determinant greater than
`28.8658...`.  Hence this determinant-24 `NS` admits no rootless MW17
fibration.

## Why the next target is genuinely `2+2`

For an `E6` rational surface of rank two, a quadratic cover ramified at one
`I1` fibre and one smooth fibre has K3 root lattice

```text
2E6 + A1,
```

of rank 13.  Generic Picard rank 19 then requires total MW rank four.  In the
deck decomposition this means two independent invariant sections and two
independent anti-invariant sections.  Merely finding two formulas on the
twist does not certify this: a second formula can be a translate of the first
by the invariant `A2*` lattice.  The current rank-four search therefore uses
the height determinant modulo the invariant lattice as its gate.

The smallest polynomial chart tested did not pass that gate.  In its simplest
simultaneous branch, the two candidate conditions reduce to

```text
r^6-9*r^4-27*r^2+27
  = (r^2+3)*(r^4-12*r^2+9) = 0.
```

It has no rational `r`: the quartic gives `r^2=6+/-3*sqrt(3)`.  This is an
obstruction only for that declared low-degree chart, not a global failure of
the rank-four target.  No rank-sum-four family is promoted by this note.

## Complete first quadratic-neighbour shell

The requested immediate low-`q` follow-up is complete.  In the split marking

```text
NS(X)=U+M(-1),     F=e,     O=f-e,
```

the first zero-neutral old-degree-two classes are

```text
D=2e+2f-w,     w.M.w=8.
```

Exact dominant-label enumeration for `2E6+A2`, followed by a closest-vector
enumeration in the rank-three root-orthogonal kernel, gives `268` Weyl
orbits.  Three are divisible degree-one presentations, leaving `265`
primitive classes.  Reduction against the old zero and all finite and affine
fibre components gives

```text
old degree 2: 211,     old degree 1: 42,     old degree 0: 12.
```

The minimum child root rank in the entire shell is thirteen.  It occurs in
exactly six classes, all of type

```text
A6+D7,     geometric MW rank 4,     regulator 6/7,
```

with primitive root lattice and trivial torsion.  All six remain at old
degree two with no physical fixed component.  They also pass the complete
old-section gate.  In the saturated basis `P,Q,S`, the old height matrix is

```text
[ 4/3  -2/3   0  ]
[ -2/3  4/3   0  ]
[ 0      0    2/3].
```

For each candidate, Proposition C2 reduces every potentially negative old
section to the strict finite ellipsoid

```text
(n-delta/2)^T H (n-delta/2) < 2.
```

There are nine points in each ellipsoid, and every tested intersection is
nonnegative.  A remaining negative irreducible root would have old degree
two.  Then `D-C` would be an effective vertical square-zero class `lF`, while
`C^2=(D-lF)^2=-4l=-2` forces `l=1/2`, contradicting primitivity of `F`.
Thus all six classes are nef Jacobian fibrations of exact geometric MW rank
four.

The sparse compiler target is orbit `215`, selected by minimum
`(max_abs,L1,orbit)=(2,12,215)` in the old `NS` basis.  Its divisor is

```text
2E6a_1 + E6a_2 + 2E6a_3 + 2E6a_4 + E6a_5 + E6b_6 + A2_1 + 2S,
```

in the certificate's deterministic component ordering.

This first shell is not R17-like: its best possible MW rank is four.  More
strongly, the determinant-24 Hermite obstruction above rules out a rootless
MW17 fibration in **every** neighbour shell, not just `q=2`.  The six MW4
fibrations are nevertheless useful next compiler targets.  Their equation
models and arithmetic descent over `QQ(r)` remain open; geometric rank four
does not yet imply that all four generators are rational.

The exact replay and generated certificate are
[`scripts/search_e6_ii_rank3_q2_neighbor_candidates.sage`](scripts/search_e6_ii_rank3_q2_neighbor_candidates.sage)
and
[`../artifacts/generated-results/elkies-k3-e6-ii-rank3-q2-neighbor-candidates-v1.json`](../artifacts/generated-results/elkies-k3-e6-ii-rank3-q2-neighbor-candidates-v1.json).
Their SHA-256 hashes are respectively
`9c34734fce6b4cf9ecd9347ae44607d8b7f1333ea4cf8cccf50262d313da5c6b`
and
`f597609d2cc6a3d2c838191195f376ffb9087aa427ea04f8c89fc3ddab37c5ff`.

## Rationalized D6 frontier

The square root in the commonly printed `D6` normal form is a normalization
artifact: it forces an additional `I1` at infinity.  Keeping the leading
coefficient free gives the rational chart

```text
y^2 = x^3
    + u^2*(-u^2+a*u-3)*x
    + u^3*(g*u^3+b*u^2+a*u-2),
```

over `QQ`.  If `-4+27*g^2 != 0`, infinity is smooth and the generic fibre
profile is `I2*+4I1`, with root lattice `D6` and geometric MW rank two.

One polynomial section is rationally parameterized by `z,q` as follows.  Put
`p=q+2` and

```text
y3 = -(64*z^4+q^4)/(64*z^3*q),
a  = (-64*z^4+3*q^4+24*q^3)/(32*z^2*q),
b  = -(-64*z^4*q+3*q^5+256*z^4+12*q^4)/(256*z^4),
g  = y3^2,

x = z^2+p*u,
y = z^3 + 3*z*p*u/2 + 3*(p^2-4)*u^2/(8*z) + y3*u^3.
```

Literal substitution proves the section identity.  With `h=q/z`, equality
of the leading `y` coefficients for two sections in this chart has the
nontrivial necessary branches

```text
h*j*(h^2+h*j+j^2)=64,
h*j*(h^2-h*j+j^2)=-64.
```

The bounded height-30 miss upgrades to an exact obstruction for this chart.
For the first branch put `k=j/h`, `x=h*j`, `X=k`, and `Y=8*k/x`; for the
second use the same variables with `X=-k`.  Both correspondence curves are
birational to

```text
Y^2=X^3+X^2+X.
```

This elliptic curve has rank zero over `QQ` and torsion
`{O,(0,0)}`.  The affine torsion point forces `k=0`, while `O` is an omitted
boundary point.  Hence neither nontrivial branch has a nondegenerate rational
point: two sections cannot coexist inside this particular polynomial marked
chart.  Replay the exact rank computation and the retained height-box
regression with
[`scripts/search_rationalized_d6_rank2_section_chart.sage`](scripts/search_rationalized_d6_rank2_section_chart.sage)
and inspect
[`../artifacts/generated-results/elkies-k3-rationalized-d6-section-chart-search-v1.json`](../artifacts/generated-results/elkies-k3-rationalized-d6-section-chart-search-v1.json).
This is not an obstruction to a second section in a larger `D6` chart.

## Status boundary

Proved here:

- a rational one-parameter `E6` family with exact rank sum three;
- explicit independent invariant sections and one explicit anti-invariant
  section of small degree;
- generic K3 Picard rank 19, saturated `NS`, and determinant 24;
- impossibility of a rootless MW17 fibration in the same `NS`;
- a complete first genuine `q=2` shell with six nef `A6+D7/MW4` frames;
- a rational `D6` chart and one exact polynomial marked-section family;
- the exact two-section obstruction in the stated D6 polynomial chart.

Still open here:

- a rational one-parameter rank-sum-four family;
- two independent anti-invariant directions in either the `E6` or `D6`
  chart;
- equation compilation and arithmetic descent for the six `A6+D7/MW4`
  neighbours.

The `E6` and `D6` surface normal forms are adapted from the rational elliptic
surface charts in [Kimura](https://arxiv.org/abs/1802.05195).  The section
slices, twist, saturation, rank-three K3 family, determinant obstruction, and
rationalized `D6` marked-section chart above are the calculations certified
in this repository.

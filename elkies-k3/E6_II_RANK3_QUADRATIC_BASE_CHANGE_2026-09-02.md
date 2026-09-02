# E6 rank-sum-three quadratic base change — 2026-09-02

<!-- status-consumer: EC-K3-RES-QBC-E6-II-RANK3-RHO19 5b10608e230145e9 -->

<!-- status-consumer: EC-K3-RES-D6-RATIONALIZED-SECTION-CHART 06442a208822d255 -->

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

An exhaustive rational-height box through 30 finds no nontrivial pair on
either branch.  Replay that bounded result with
[`scripts/search_rationalized_d6_rank2_section_chart.sage`](scripts/search_rationalized_d6_rank2_section_chart.sage)
and inspect
[`../artifacts/generated-results/elkies-k3-rationalized-d6-section-chart-search-v1.json`](../artifacts/generated-results/elkies-k3-rationalized-d6-section-chart-search-v1.json).
The box is not a rational-point theorem and does not exclude larger D6
section charts.

## Status boundary

Proved here:

- a rational one-parameter `E6` family with exact rank sum three;
- explicit independent invariant sections and one explicit anti-invariant
  section of small degree;
- generic K3 Picard rank 19, saturated `NS`, and determinant 24;
- impossibility of a rootless MW17 fibration in the same `NS`;
- a rational `D6` chart and one exact polynomial marked-section family;
- the stated bounded D6 necessary-condition search.

Still open here:

- a rational one-parameter rank-sum-four family;
- two independent anti-invariant directions in either the `E6` or `D6`
  chart;
- a same-NS rootless search for a future rank-four success.

The `E6` and `D6` surface normal forms are adapted from the rational elliptic
surface charts in [Kimura](https://arxiv.org/abs/1802.05195).  The section
slices, twist, saturation, rank-three K3 family, determinant obstruction, and
rationalized `D6` marked-section chart above are the calculations certified
in this repository.

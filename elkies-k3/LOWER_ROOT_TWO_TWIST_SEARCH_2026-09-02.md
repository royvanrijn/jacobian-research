# Lower-root simultaneous two-twist search — 2026-09-02

<!-- status-consumer: EC-K3-RES-QBC-E6-II-RANK3-RHO19 5b10608e230145e9 -->

<!-- status-consumer: EC-K3-RES-D6-RATIONALIZED-SECTION-CHART a94042dd2d76797c -->

<!-- status-consumer: EC-K3-RES-D5-TWO-MARKED-TWO-TWIST-POLYNOMIAL ea0496c9566cfdc3 -->

<!-- status-consumer: EC-K3-RES-A4-TWO-POINT-TATE-SLICE-OBSTRUCTION b9729a0a8f2f17be -->

## Outcome

No `2+2` family over `QQ(u)` is promoted.  The simultaneous lower-root pass
nevertheless gives four useful reductions:

1. the old `D6` height-box miss is an exact obstruction for that polynomial
   marked-section chart; and
2. the complete good non-`j=0` `GF(11)` E6 shared-simple-pole census has only
   two mechanisms.  One is a dependent multiple, and the other forces a
   conic with no nondegenerate rational point when the invariant rank-two
   marking is imposed; and
3. a rational section-first D5 seed has two exactly independent invariant
   sections, and its complete `GF(11)` polynomial twist chart has one
   two-section modular candidate; and
4. the simplest two-point A4 Tate slice forces either an extra repeated
   discriminant root (generically an `I2` fibre) or dependence of its marked
   points.

Thus D5 has passed the first modular discovery gate, although it has not
passed rational lifting or the height-determinant gate.  Increasing the old
D6 height box or lifting the raw mod-11 E6 hits would repeat a mechanism
already rejected exactly.

## Rank budget

For an `E6` or `D6` rational elliptic surface, the geometric root rank is six
and the geometric Mordell--Weil rank is two.  A quadratic base change with one
additional `A1` root has K3 root rank thirteen.  A Picard-rank-19 member can
therefore have total Mordell--Weil rank four, with character split

```text
rank E(QQ(u)) + rank E^(d)(QQ(u)) = 2+2.
```

The two twist formulas must have nonzero height determinant.  Merely listing
two points is insufficient: the first E6 modular component below contains the
universal pair `S,-2S`.

## Shared-pole E6 system

Keep the certified rank-two E6 surface in the polynomial marked chart

```text
E: y^2=x^3+a*u*x+u*(u+2-a),
a*(r+1)=2*(r^2+r+1).
```

For a monic quadratic `d=u^2+d1*u+d0` and one shared pole `H=u-h`, impose two
twist sections simultaneously:

```text
x_i=X_i/H^2,       deg X_i <= 2,
y_i=Y_i/H^3,       deg Y_i <= 3,
d*Y_i^2=X_i^3+a*u*X_i*H^4+u*(u+2-a)*H^6.
```

After fixing each `Y_i` monic, degrees seven, six, and five recover its other
three coefficients recursively.  Five residual equations remain per section.
Together with the marked-surface relation, this is eleven equations in eleven
geometric variables.  An inverse variable implements each of the three open
charts covering `X_1 != X_2`.

Direct modular msolve pilots on the unsliced 12-variable saturation charts
timed out at 45 seconds per chart.  This is a solver-size result, not an empty
scheme.  Fibrewise enumeration is much smaller: after fixing the four global
parameters, only `p^3` possible `X` polynomials remain and `Y` is forced.

The replay is
[`scripts/search_e6_shared_pole_two_twist_sections_modp.sage`](scripts/search_e6_shared_pole_two_twist_sections_modp.sage):

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/search_e6_shared_pole_two_twist_sections_modp.sage \
  --prime 11 --enumerate --skip-msolve
```

The generated summary is
[`../artifacts/generated-results/elkies-k3-e6-shared-pole-two-twist-sections-p11-v1.json`](../artifacts/generated-results/elkies-k3-e6-shared-pole-two-twist-sections-p11-v1.json).

## Complete `GF(11)` result

The good non-`j=0` census tests

```text
7,260 global parameter/pole presentations,
9,663,060 candidate X polynomials.
```

It finds eighteen pole presentations with two monic-`Y` sections.  After
canonicalizing the rational functions across different pole presentations,
there are twelve `r`-marked surface/twist records; the pairs of `r` values in
each repeated row define the same `(a,c)` surface.  The six underlying
surface/twist records have only the following mechanisms.

### Dependent `d=B` component

Here

```text
d=u*(u+c),       S=(0,1).
```

The second displayed formula is `-2S`.  It is not a second twist direction.
This is the group-law dependence that a point-counting search would miss.

### Constant-section component

The other records have

```text
d=u^2-h^2
```

and contain a constant twist section `(x,y)=(k,1)`.  Comparing coefficients
gives

```text
a*k+c=0,       k^3=-h^2,       a+c=2.
```

Over `QQ`, write `k=-s^2`, `h=s^3`.  Then

```text
a=2/(1+s^2),       c=2*s^2/(1+s^2).
```

Imposing the rational rank-two E6 marking gives

```text
(1+s^2)*r^2+s^2*r+s^2=0,
disc_r=-s^2*(4+3*s^2).
```

For rational nonzero `s`, this discriminant is negative and cannot be a
rational square.  The case `s=0` makes `d=u^2` a square and is not a quadratic
twist.  Thus this entire simultaneous modular mechanism cannot yield `2+2`
over `QQ(u)`, even though it has valid points over several finite fields.

This conclusion is chart-specific: it does not exclude a different quadratic
`d`, higher section degree, or a non-shared/simple pole divisor.

## Exact D6 chart obstruction

In the rationalized D6 chart from
[`E6_II_RANK3_QUADRATIC_BASE_CHANGE_2026-09-02.md`](E6_II_RANK3_QUADRATIC_BASE_CHANGE_2026-09-02.md),
two copies of the existing polynomial section parametrization must satisfy

```text
h*j*(h^2+h*j+j^2)=64
```

or

```text
h*j*(h^2-h*j+j^2)=-64.
```

Both nontrivial correspondence curves are birational to

```text
Y^2=X^3+X^2+X.
```

That elliptic curve has rank zero over `QQ` and torsion
`{O,(0,0)}`; both torsion points lie on degenerate or omitted boundaries of
the section correspondence.  Hence the declared D6 polynomial chart contains
no nontrivial rational pair.  This is exact for that chart, but says nothing
about a larger D6 rational-function section ansatz.

## Section-first D5 experiment

Kimura's D5 chart becomes rational after replacing the leading pair
`(-1,2/(3*sqrt(3)))` by the rational cusp parametrization `(-3*q^2,2*q^3)`
and taking `q=1`.  In affine coordinate `u` it is

```text
E_(a,b,c): y^2=x^3+u^2*(-3*u^2+a*u-3)*x
                    +u^3*(2*u^3+b*u^2+c*u-2).
```

A low-height section

```text
x=u*(l0+l1*u),       y=u^2*(m0+m1*u)
```

forces `l0=2` or `-1` and

```text
m1^2=(l1-1)^2*(l1+2).
```

Thus `l1=t^2-2`, `m1=+/-t*(t^2-3)` rationalizes the section equation.
Putting two sections on the `l0=-1` branch and choosing

```text
(l1,m1,m0)=(-2,0,1),       (-1,-2,1)
```

gives the exact rational seed

```text
(a,b,c)=(-13,-17,-12),
P=(-u-2*u^2, u^2),
Q=(-u-u^2, u^2*(1-2*u)).
```

Its discriminant has order seven at `u=0`, and the residual finite factor is
squarefree, so the fibre is `I1*` and the other finite fibres are nodal.  At
`u=1`, `P,Q` specialize to a signed basis on

```text
y^2=x^3-19*x-29,       rank E(QQ)=2.
```

Any relation over `QQ(u)` would specialize to a relation there, so this
certifies two independent invariant directions without assuming that all of
the geometric `A3*` Mordell--Weil lattice descends.

On this fixed seed, impose the common pole-free twist chart

```text
d=u^2+d1*u+d0,       d squarefree,
d*y_i^2=x_i^3+f*x_i+g,       deg x_i,deg y_i <= 2.
```

The replay is
[`scripts/search_d5_two_marked_two_twist_polynomial_modp.sage`](scripts/search_d5_two_marked_two_twist_polynomial_modp.sage):

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/search_d5_two_marked_two_twist_polynomial_modp.sage \
  --prime 11
```

The complete `GF(11)` census checks 110 squarefree monic quadratics, 1,331
quadratic `x` polynomials, and 146,410 incidences.  Fourteen incidences occur
on thirteen twists.  Exactly one twist has two distinct-`x` sections:

```text
d=u^2-2*u-2,
(x1,y1)=(-u+7*u^2, 4*u^2),
(x2,y2)=(1-u+3*u^2, 7+8*u+3*u^2)             over GF(11).
```

There is no relation `m*P1+n*P2=0` with `|m|,|n|<=12` in the exact
`GF(11)(u)` group calculation.  This bounded relation check is not an
independence proof.  Analogous complete searches find one two-section twist
over each of `GF(7)` and `GF(13)`, and two over `GF(17)`, so the modular
phenomenon is not isolated to characteristic eleven.

The direct 15-variable inverse-saturated modular Groebner systems are also
exported in three coefficient-pivot charts.  At `GF(11)` all three exceeded a
60-second `msolve` pilot; the fibrewise census is the algorithmic reduction
that completely resolves rational points of the declared finite-field chart.
A search over integral `x` coefficients through height 28 found no rational
pair with a common twist.  That last search is only a bounded lifting pilot.

The generated modular summary is
[`../artifacts/generated-results/elkies-k3-d5-two-marked-two-twist-polynomial-p11-v1.json`](../artifacts/generated-results/elkies-k3-d5-two-marked-two-twist-polynomial-p11-v1.json).

## Exact A4 two-point Tate-slice obstruction

The simplest section-first A4 entrance uses

```text
y^2+a1*x*y+a3*y=x^3+a2*x^2,
a1=1+A*u,       h=H*u,       kappa=B*u+C*u^2,
a2=kappa,       a3=h*(1-h-a1)+kappa.
```

It contains `P=(0,0)` and `Q=(h,h^2)` identically.  The first four
discriminant coefficients vanish automatically, while

```text
[u^4] Delta = -B^2*(H+B)*(A+H-B).
```

On the open chart `B != 0`, an `I5` fibre therefore has two branches,
`H=-B` and `H=B-A`.  On both, exact factorization gives

```text
Delta/u^5 = -(B+(A*B-B^2+C)*u)^2 * cubic(u).
```

Unless `C=B^2-A*B`, the slice has a forced repeated discriminant root,
generically an extra `I2` fibre.  Imposing that condition to remain in the
pure A4 stratum makes the marked points dependent:

```text
H=-B:       P+2*Q=0,
H=B-A:    3*P+2*Q=0.
```

The exact replay is
[`scripts/certify_a4_two_point_tate_slice_obstruction.sage`](scripts/certify_a4_two_point_tate_slice_obstruction.sage),
with generated summary
[`../artifacts/generated-results/elkies-k3-a4-two-point-tate-slice-obstruction-v1.json`](../artifacts/generated-results/elkies-k3-a4-two-point-tate-slice-obstruction-v1.json).
This is an obstruction only for the normalized `r=s=1` two-point Tate slice;
it says that the next A4 chart must restore the more general Bezout data, not
that a two-marked A4 surface is impossible.

## Next experiment

The D5 candidate should next be lifted by slicing the shared coefficient
ideal around its modular points, rather than attacking the unsliced Groebner
system.  A successful lift still needs an exact twist-height determinant.
In parallel, the A4 chart should restore nonconstant coprime Bezout data
`r,s` in the two-point Tate construction, impose the `I5` discriminant jets
afterward, and retain surface moduli in the simultaneous twist equations.

For either route, quotient the obvious `S,-2S` and constant-section components
before modular elimination.  Any lift must pass a twist-height determinant
gate before it is counted as `2+2` or `3+2`.

## Status boundary

Exact here:

- the D6 correspondence-curve obstruction inside the declared polynomial
  marked-section chart;
- construction and complete `GF(11)` enumeration of the declared E6
  shared-simple-pole ansatz;
- the exact algebraic rejection over `QQ` of the two mechanisms containing
  every E6 modular survivor;
- the D5 seed, its fibre calculation, its two invariant sections, and their
  independence by exact specialization;
- the complete `GF(11)` enumeration of the declared D5 polynomial twist
  ansatz; and
- the exact dependence/extra-repeated-fibre dichotomy in the declared A4
  two-point Tate slice.

Not proved here:

- nonexistence of a rational `2+2` family on E6, D6, D5, or A4 surfaces;
- a characteristic-zero lift or twist-height matrix for the D5 modular pair;
  or
- any `3+2` or rank-sum-five construction.

The surface normal forms follow
[Kimura's rational elliptic-surface charts](https://arxiv.org/abs/1802.05195).
The simultaneous two-point viewpoint is compatible with the global `dP2`
model of
[Cvetic--Klevers--Piragua](https://arxiv.org/abs/1303.6970); neither source
asserts the arithmetic rank conclusions tested here.

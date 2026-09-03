# Lower-root simultaneous two-twist search — 2026-09-02

<!-- status-consumer: EC-K3-RES-QBC-E6-II-RANK3-RHO19 5b10608e230145e9 -->

<!-- status-consumer: EC-K3-RES-D6-RATIONALIZED-SECTION-CHART a94042dd2d76797c -->

<!-- status-consumer: EC-K3-RES-D5-TWO-MARKED-TWO-TWIST-POLYNOMIAL ea0496c9566cfdc3 -->

<!-- status-consumer: EC-K3-RES-D5-TWO-MARKED-LOW-SLICE-ELIMINANTS 43d297285eb3655b -->

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
   sections, its complete `GF(11)` polynomial twist chart has one two-section
   modular candidate, and exact eliminants rule out rational points on the
   regular low-section slices through the `p=11,13` candidates; and
4. the simplest two-point A4 Tate slice forces either an extra repeated
   discriminant root (generically an `I2` fibre) or dependence of its marked
   points.

Thus D5 passed the first modular discovery gate.  Exact elimination now kills
the two regular low-section slices through its `p=11,13` survivors: their
selected characteristic-zero points have degrees 88 and 78 over `QQ`, and
neither complete saturated slice has a rational point.  This does not close
the full D5 polynomial chart.  Increasing the old D6 height box or lifting the
raw mod-11 E6 hits would repeat a mechanism already rejected exactly.

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

### Regular low-section slices through the `p=11,13` survivors

The `p=11` and `p=13` pairs admit a smaller exact local chart.  Restrict the
first twist section to

```text
x=u*(A+B*u),       y=C*u^2.
```

The four nonzero coefficient equations factor triangularly as

```text
(A-2)*(A+1)^2=0,
C^2*d0=3*A^2*B-13*A-3*B-12,
C^2*d1=3*A*B^2-3*A-13*B-17,
C^2=(B-1)^2*(B+2).
```

On the two branches containing the modular survivors, put

```text
B=t^2-2,       C=t*(t^2-3).
```

Then the monic twist coefficients are

```text
A=-1:
  d0=1/(t^2*(t^2-3)^2),
  d1=-(3*t^2+1)/(t^2-3)^2;

A=2:
  d0=(9*t^2-56)/(t^2*(t^2-3)^2),
  d1=(6*t^4-37*t^2+27)/(t^2*(t^2-3)^2).
```

The `p=11` point lies on the first branch at `t=8`, and the `p=13` point
lies on the second at `t=2`.  After substituting a general second polynomial
section, its seven coefficient equations form a square system in `t` and the
six section coefficients.  The coefficient Jacobians have full rank seven;
their determinants are `-2 mod 11` and `-3 mod 13`.  Thus both modular points
are reduced isolated points of these slices and each has a unique multivariate
Hensel lift over `Z_11` or `Z_13`.
Exact group arithmetic over each finite function field verifies that the two
sections are neither equal nor opposite and that neither is plus or minus
twice the other.  Their higher coefficients are nonzero.  Thus these local
points are off the repeated, `S,-2S`, and constant-section components.

Lifting to 800 p-adic digits did not produce a rational tuple.  At the
declared checkpoints, coefficientwise rational reconstruction was either
undefined or failed literal substitution in the seven exact equations.  The
one accidental simultaneous reconstruction at 400 digits on the `p=11`
branch had roughly 690-bit coordinates, failed substitution, and did not
persist to 800 digits.  This remains a reproducible negative lifting
experiment; the rigorous nonrationality statement comes instead from the
exact eliminants below.

### Exact local eliminants

Eliminate the six coefficients of the second section and the saturation
inverse, leaving `t`.  Exact block Groebner elimination over `QQ` gives:

```text
A=-1: degree 142, irreducible-factor degrees 2,2,4,4,42,88;
A= 2: degree 128, irreducible-factor degrees 2,4,44,78.
```

Both eliminants are even.  Their full primitive coefficient vectors and all
irreducible factors are pinned in
[`../artifacts/generated-results/elkies-k3-d5-two-marked-two-twist-local-eliminants-v1.json`](../artifacts/generated-results/elkies-k3-d5-two-marked-two-twist-local-eliminants-v1.json).
The residue `t=8 mod 11` selects the irreducible degree-88 factor in the first
line, and `t=2 mod 13` selects the irreducible degree-78 factor in the second.
In both cases the selected root is simple modulo the indicated prime.  Thus
the exact `t`-fields of the two isolated lifts are

```text
K_11 = QQ[t]/(H_88(t)),       K_13 = QQ[t]/(H_78(t)),
```

where `H_88` and `H_78` are the primitive selected factors in that artifact.
In particular neither target is rational.  More strongly, neither complete
saturated slice has a `QQ` point: a rational tuple would give a rational `t`,
but the displayed exact factorizations contain no linear factor over `QQ`.
This closes the two local slices, not the full fourteen-variable D5
coefficient ideal or a different first-section slice.

The replay is
[`scripts/lift_d5_two_marked_two_twist_low_section_slices.sage`](scripts/lift_d5_two_marked_two_twist_low_section_slices.sage):

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/lift_d5_two_marked_two_twist_low_section_slices.sage \
  --digits 800
```

It writes
[`../artifacts/generated-results/elkies-k3-d5-two-marked-two-twist-low-slices-v1.json`](../artifacts/generated-results/elkies-k3-d5-two-marked-two-twist-low-slices-v1.json)
and exact characteristic-zero `msolve` inputs under
`artifacts/local/elkies-k3/d5-two-marked-two-twist-low-slices/`.  The exported
systems invert `x20*t*(t^2-3)`: this removes the tautological repeated-section
component and the denominator-clearing boundary while retaining both target
points.  Run the exact block eliminations with `--run-eliminants`; the first
uses the longer resource envelope.  The `p=7` survivor and the two `p=17`
boundary survivors do not lie in these two low-section slices.

### The `p=7` survivor is obstructed modulo `7^3`

At the ordered `p=7` pair, the full fourteen-equation coefficient Jacobian
has rank thirteen.  Its right kernel in the script's coefficient order is

```text
(1,2,1,6,2,3,6,5,1,0,5,3,6,0),
```

and its left kernel is

```text
(1,4,2,1,4,2,1,0,0,0,0,0,0,0).
```

Solving the first underdetermined Hensel equation gives all seven lifts from
modulo `7` to modulo `7^2`.  For every one of these lifts, pairing the next
Hensel right-hand side with the displayed left kernel gives `1 mod 7`.
Consequently none lifts modulo `7^3`.  In particular there is no solution of
the shared coefficient ideal over `Z_7` reducing to this ordered point, and
therefore no rational coefficient tuple integral at seven which reduces to
it.  This is an exact vertical-characteristic-seven obstruction, not a
bounded rational reconstruction failure.

The twist discriminant and the coefficient pivot `x11-x21` are units at the
point.  Exact group arithmetic over `GF(7)(u)` also checks that the two points
are neither equal nor opposite and that neither is plus or minus twice the
other.  Thus the obstruction is computed off the repeated-section and
`S,-2S` components.
The nonzero higher coefficients also exclude the constant-section component.

### The two `p=17` survivors are a bad-fibre boundary mechanism

Both `p=17` pairs have

```text
d=u*(u+k),
x_i=u*(a_i+b_i*u),
y_i=u*(c_i+e_i*u),       c_i != 0.
```

Before setting `d0=0`, the coefficient of `u^2` in each section identity is
`d0*c_i^2`.  Hence the local unit condition `c_i!=0` forces `d0=0` exactly.
After dividing by `u^3`, each section is governed by

```text
a^3-3*a-2-k*c^2=0,
3*a^2*b-2*k*c*e-c^2-13*a-3*b-12=0,
3*a*b^2-k*e^2-2*c*e-3*a-13*b-17=0,
b^3-e^2-3*b+2=0.
```

The shared eight-equation system in the nine variables
`k,a1,b1,c1,e1,a2,b2,c2,e2` has Jacobian rank eight at both modular hits.
Thus they are smooth points of an expected boundary curve, not isolated
points of the open K3-twist locus.  Exact `GF(17)(u)` group arithmetic again
excludes equality, opposition, and either plus-or-minus-double relation.
Nonzero higher coefficients exclude constant sections as well.

This boundary cannot supply the intended D5 K3 source.  Indeed the short
twist has coefficients `d^2*f,d^3*g`; minimalizing at `u=0` by
`X=u^2*X', Y=u^3*Y'` gives

```text
A_min=(u+k)^2*(-3*u^2-13*u-3),
B_min=(u+k)^3*(2*u^3-17*u^2-12*u-2).
```

These have degrees four and six.  The affine discriminant has degree eleven,
with the twelfth discriminant zero at infinity.  Generically the old `I1*`
becomes `I1` at `u=0`, while the other zero of `d` gives `I0*`; infinity is
also `I1`.  The old discriminant is a unit at `u=-k` for both displayed
modular values.  The minimal twist therefore has `chi=1`: it is a rational
elliptic surface, not the `chi=2` K3 required by the declared
arithmetic-source/NS gate.
Saturating the coefficient ideal by `d0` removes both `p=17` survivors.

For completeness, fixing the displayed residues `k=5` and `k=8` makes each
boundary section system square with nonzero Jacobian determinant.  Their
unique fixed-`k` Hensel lifts were followed for 400 digits; neither yielded a
rational tuple under literal reconstruction.  This last statement is only a
fixed-slice negative experiment.  The exact reason for rejecting these two
survivors from the current search is the rational-surface minimalization, not
the reconstruction failure.

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

The two regular `p=11,13` low-section slices are now closed by exact
nonrationality, the `p=7` point cannot lift even modulo `7^3`, and the two
`p=17` points are removed by the `d0!=0` K3-source gate.  Any remaining D5
work must therefore change the first-section slice or eliminate the larger
coefficient ideal; collecting more primes for these five survivors repeats a
closed mechanism.  A new characteristic-zero lift would still need an exact
twist-height determinant.
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
  ansatz;
- the exact low-section parameterization and full-rank local Jacobians at the
  `p=11,13` D5 survivors, giving unique p-adic lifts in those slices;
- the exact degree-142 and degree-128 `t`-eliminants for those saturated
  slices, their complete rational factorizations, the irreducible degree-88
  and degree-78 factors selected by the two modular points, and consequently
  the absence of a rational point on either complete slice;
- the exact second-order Hensel obstruction eliminating the `p=7` survivor;
- the forced `d0=0` boundary classification and rational-surface
  minimalization eliminating both `p=17` survivors from the K3-source locus;
  and
- the exact dependence/extra-repeated-fibre dichotomy in the declared A4
  two-point Tate slice.

Not proved here:

- nonexistence of a rational `2+2` family on E6, D6, D5, or A4 surfaces;
- a characteristic-zero rational lift elsewhere in the D5 coefficient chart,
  or a twist-height matrix for a D5 modular pair;
  or
- any `3+2` or rank-sum-five construction.

The surface normal forms follow
[Kimura's rational elliptic-surface charts](https://arxiv.org/abs/1802.05195).
The simultaneous two-point viewpoint is compatible with the global `dP2`
model of
[Cvetic--Klevers--Piragua](https://arxiv.org/abs/1303.6970); neither source
asserts the arithmetic rank conclusions tested here.

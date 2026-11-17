# Complete linear carrier-subspace limits of the saved cubic

KCUB29,15 November2026. Fix the integral twelve-variable cubic KCUB3,
its three base coordinates and its nine displayed triangular carrier
outputs C4,...,C12. Impose any homogeneous linear relations among these
nine outputs, with arbitrary coefficients in a characteristic-zero field.
Normalize the relations by highest-pivot row reduction, solve for those
source coordinates, and retain the original base variables, raw free
source coordinates and corresponding raw output rows. No further target
repair or nonlinear source reparametrization is allowed in this class.

Every such graph map is Keller and retains the original integral collision.
Across the entire class, of all dimensions3 through12 and all resulting
degrees, the minimum raw nonlinear occurrence count is exactly13. The
minimum ordinary degree at each dimension is exactly:

| Dimension |3|4|5|6|7|8|9|10|11|12|
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Minimum degree in this class |7|6|5|5|4|4|4|4|4|3|

At the smallest dimensions attaining degrees6,5,4, the corresponding
minimum occurrence counts are15,17,23. Every7-variable quartic in this
class is coefficientwise the saved23-occurrence map. These results
extend the [literal-subset theorem](CARRIER_SUBSET_FRONTIER.md) to arbitrary
simultaneous linear carrier relations, including all cross-weight mixes.
They do not establish any global minimum, any affine-orbit term bound,
any full Pareto classification at intermediate dimensions/degrees, or an
exclusion of nonlinear graphs or subsequent target repairs. The global
ordinary-cubic dimension interval remains[4,12], and all established
headline upper witnesses and the62/12 strict3D profile partition remain
unchanged. The full simplicity goal is open.

## Complete parametrization and Keller inheritance

For each retained carrier index set T, let D be its complement. Every
normalized relation space in its highest-pivot echelon chart has equations

```
Hj = Cj + sum_(k in T, k<j) t_jk*Ck = 0,   j in D.
```

There are binomial(9,|T|) charts in retained dimension3+|T| and512 charts
in total, each with arbitrary scalar parameters, up to20 parameters.
Reverse row reduction over a field covers every relation subspace, including
all boundary strata; the choice of free columns determines T uniquely.
There is no generic nonzero-parameter assumption.

Because Cj=Xj+Pj and Pj uses only earlier source coordinates, Hj is monic
in Xj and independent of later source coordinates. Replacing the deleted
source coordinates by Hj is a triangular polynomial automorphism with
determinant1. Replace the same target carrier rows by their linear Hj
combinations, also with determinant1. The deleted outputs then equal the
new source coordinates. Restricting these identity variables to zero gives
a complementary Jacobian determinant1. All Cj vanish at both original
collision points, so the relations vanish there for every parameter;
projection keeps the unequal base coordinates and therefore a collision.
The checker verifies the actual graph identities, identity linear jets
and both parameter-polynomial source embeddings and collisions in all512
charts. The chain-rule/block-determinant inheritance is a written argument.

The parametrization fixes the raw source complement selected by the
highest-pivot chart. Other, possibly nonlinear, parametrizations of the
same graph are outside this theorem. Additional affine mixing of retained
outputs can change the occurrence count and is not included in that count
bound. Ordinary degree itself is unchanged by invertible affine source
and target changes of an already obtained map.

## Degree certificates

For a proposed degree cap d, expand every retained output and collect its
source-degree-above-d coefficients as polynomials in the chart parameters.
A constant rational linear combination of these coefficient polynomials
equal to1 excludes that chart. All508 excluded charts have such a
certificate;499 have an individual nonzero constant coefficient and nine
need a combination of several coefficients. No Groebner output, primary
row reduction, parameter sampling or finite-field assertion is trusted.

| Dimension | Degree cap tested | Charts | Excluded by unit identity | Surviving chart |
| ---: | ---: | ---: | ---: | --- |
|3|6|1|1|none|
|4|6|9|8|T={5}|
|5|5|36|35|T={5,8}|
|6|4|84|84|none|
|7|4|126|125|T={4,5,9,10}|
|8|3|126|126|none|
|9|3|84|84|none|
|10|3|36|36|none|
|11|3|9|9|none|
|12|3|1|0|the original map|

The three surviving proper charts have fixed nonzero coefficients of
source degrees6,5,4 respectively, as reconstructed below, so their degrees
cannot drop further. For dimension6 the degree4 exclusions give the lower
bound5. Literal subsets already give degree5 in6 variables and degree4
in8 through11 variables. The checker matches all ten saved literal upper
witnesses to parameter-zero specializations of the current reconstructed
maps. This proves both sides of every entry in the degree table.

## The four-variable sextic chart

Here T={5}, U=C5=x5+x2^3, and the parameters t0,...,t6 are the coefficients
of C5 in the deleted C6,...,C12 relations respectively. The sole high-degree
coefficient equation is12*t2=0. Once t2=0, the map differs from the saved
4-variable sextic G0 by

```
G1 = G01 + t0*t1*U^2
G2 = G02 - 12*t1*t6*U^2
G3 = G03 - 4*t1*t5*U^2,
```

and its retained carrier is U. The independent checker reconstructs all
three identities. The coefficient of x1^3*x2^2*x3 in row2 stays12, forcing
degree6. Each nonzero added square has exactly three new nonlinear
monomials, disjoint from the original row supports. Thus the exact count is

```
15+3[t0*t1!=0]+3[t1*t6!=0]+3[t1*t5!=0].
```

In particular the minimum here is15, attained by the original map.

## The five-variable quintic chart

Here T={5,8}, U=C5=x5+x2^3 and V=C8=x8+x2^2. The parameter pairs are

```
t0:(6,5), t1:(7,5), t2:(9,5), t3:(9,8),
t4:(10,5), t5:(10,8), t6:(11,5), t7:(11,8),
t8:(12,5), t9:(12,8).
```

The complete degree5 conditions are

```
t0*t1 = t1*t8 = t1*t6 = 0.
```

Their necessity is independently certified by high-coefficient linear
combinations; direct full-map identities prove sufficiency. Before imposing
them, the differences from G0 are

```
G1-G01 = t0*t1*U^2
G2-G02 = -12*t1*t8*U^2 - 12*t1*t9*U*V
G3-G03 = -4*t1*t6*U^2 - 4*t1*t7*U*V.
```

The U^2 terms disappear under the conditions. Row2 keeps coefficient-12
on x1^3*x3*x8, so degree5 is forced. Each nonzero U*V contributes four
new monomials disjoint from its old row. The exact count is

```
17+4[t1*t9!=0]+4[t1*t7!=0].
```

Hence17 is sharp in this entire quintic chart.

## The seven-variable quartic chart

Here T={4,5,9,10}. The parameter pairs are

```
t0:(6,4),  t1:(6,5),  t2:(7,4),  t3:(7,5),
t4:(8,4),  t5:(8,5),
t6:(11,4), t7:(11,5), t8:(11,9), t9:(11,10),
t10:(12,4),t11:(12,5),t12:(12,9),t13:(12,10).
```

Nineteen exact high-coefficient combinations force t4=t5=0 and all the
coefficients of the three formal products

```
(t2*Y+t3*Z)*(t0*Y+t1*Z),
(t2*Y+t3*Z)*(t6*Y+t7*Z+t8*W+t9*V),
(t2*Y+t3*Z)*(t10*Y+t11*Z+t12*W+t13*V)
```

to vanish. Over a field, the polynomial ring in the independent formal
symbols Y,Z,W,V is a domain. Thus the parameter locus is exactly the union
of these two linear branches:

```
A: t2=t3=t4=t5=0, other ten parameters free;
B: t0=t1=t4=t5=t6=t7=t8=t9=t10=t11=t12=t13=0,
   t2,t3 free.
```

Direct substitution verifies sufficiency and proves that both branches
give exactly the saved7-variable23-occurrence map, with no parameter
remaining. Its row1 coefficient of x1^3*x3 is-1, forcing degree4. There is
therefore no smaller count or new quartic map in this complete chart.
The field/domain step is written, not proof-assistant formalized.

## Thirteen occurrences are forced in the entire class

Collect all nonlinear source coefficients, not only the high-degree ones.
The minimum numbers of parameter-independent nonzero coefficients per
chart in dimensions3 through12 are

```
13,14,13,12,13,15,16,22,29,36.
```

Every chart except one already has at least13 such distinct slots. The
exception is the6-variable chart T={4,9,10}, which has12 fixed slots.
With parameters ordered by its echelon pairs, two additional row3 slots,
for x1*x2^3*x4 and x1^2*x2^2*x3, have coefficients

```
P=8*(t2*t5-1),   Q=12*t2*t5.
```

They are different source monomials and are not among the12 fixed slots.
The exact identity -P/8+Q/12=1 ensures at least one further occurrence at
every parameter. Thus every map in all512 charts has at least13 raw
nonlinear occurrences, without any degree restriction. The original
3-variable13-occurrence map attains this bound. The checker reconstructs
all coefficient groups and the exceptional two-coefficient identity.
This count concerns the fixed retained coordinates only.

## Independent replay and trust boundary

Run

```
.venv/bin/python research/jacobian/cubic-reduction/verify_linear_carrier_subspaces12.py
```

The primary [chart producer](linear_carrier_subspaces12.py) recursively
builds a source embedding and composes the original map. The independent
standard-library Fraction checker instead updates the live map one graph
substitution at a time. It reconstructs every source embedding and output,
every high coefficient, every constant unit, all normal-form identities,
the nineteen quartic generators and both quartic branches, all generic
collisions, every fixed-slot lower bound and all ten positive upper
specializations. The [certificate producer](linear_carrier_certificates12.py)
uses rational row elimination only to discover constant combinations;
the checker merely multiplies them exactly.

- [512-chart coefficient certificates](linear-carrier-subspaces12/certificates.json)
- [Independent receipt and all model hashes](linear-carrier-subspaces12/validation.json)
- [Earlier independently checked positive determinants](carrier-subset-frontier12/validation.json)

The complete mathematical conclusions combine these finite exact checks
with the written echelon-cover, triangular Jacobian, polynomial-domain and
support-count arguments. They are neither formal verification nor external
review. Preliminary weighted-only screens are retained as discovery
history; this theorem uses the full unweighted charts. Primary unresolved
parameter equations are resolved by the subsequent exact certificates,
without erasing their original checkpoints. No bounded search failure is
used as an exclusion.

Further improvement requires leaving the stated raw linear-carrier graph
class: nonlinear graph relations, a changed source parametrization,
subsequent target repair or output-frame sparsification, a changed carrier
realization or a different construction. These routes and all global
simplicity optima remain open.

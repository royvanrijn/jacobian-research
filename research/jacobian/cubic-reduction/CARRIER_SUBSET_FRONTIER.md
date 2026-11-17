# Exact literal carrier-subset frontier of the saved cubic

For the integral twelve-variable KCUB3 map, all512 subsets of its nine
triangular carrier equations have been independently replayed. Setting a
chosen subset of carrier outputs to zero and deleting their source/target
coordinates gives determinant-one polynomial maps with distinct integral
collisions. No target repair is applied.

The nondominated triples (dimension, ordinary degree, expanded nonlinear
monomial occurrences) in this finite family are exactly:

| Dimension | Degree | Nonlinear occurrences | Deleted carrier indices | Explicit coefficients and collision |
| ---: | ---: | ---: | --- | --- |
|3|7|13|4,5,6,7,8,9,10,11,12|[witness](carrier-subset-frontier12/best-n03.json)|
|4|6|15|4,6,7,8,9,10,11,12|[witness](carrier-subset-frontier12/best-n04.json)|
|5|5|17|4,6,7,9,10,11,12|[witness](carrier-subset-frontier12/best-n05.json)|
|7|4|23|6,7,8,11,12|[witness](carrier-subset-frontier12/best-n07.json)|
|12|3|36|none|[witness](carrier-subset-frontier12/best-n12.json)|

Nondominated means that no other member of these512 subsets is no larger
in all three metrics and strictly smaller in at least one. It does not
mean optimal among all counterexamples, all affine frames, or all graph
restrictions. The intermediate constructions give useful explicit upper
witnesses in degrees4,5,6. They do not improve the saved degree3 dimension12,
36-term cubic, raw3D13-term or cubic-linear89-variable witnesses.

## A four-variable integral example

Writing the retained variables as (a,b,c,u)=(x1,x2,x3,x5), put d=3b+ac.
The degree6,15-occurrence map is

```
A = a-a^2*d
B = b+3*a*c+12*a*b*(2*b+a*c)+12*a^2*b^2*d
C = c+8*b^2+8*a*u+6*a*b*c+12*a*b^2*d-8*a^2*u*d
U = u+b^3.
```

Its Jacobian determinant is identically1. The distinct integer points
(0,1,-8,-1) and (-1,1,4,-1) both map to (0,1,0,0).
Expanding and collecting gives nonlinear row counts2,5,7,1, hence15.
The witness file, rather than this factored display, supplies every exact
expanded coefficient used by the checker.

## Graph construction and inheritance

Write K=(K1,K2,K3,C4,...,C12), with each Cj=Xj+Pj depending only on
base variables and earlier carriers. For a subset S, replace Xj by -Pj
recursively for each j in S, and retain the other outputs and variables.
The carrier map changing Xj to Cj only for j in S is triangular with
Jacobian determinant1. In those source coordinates, the deleted outputs
are exactly the deleted identity variables. At their zero fibre, the full
Jacobian has an identity output block; its complementary determinant is
the determinant of the restricted map, hence1. The two KCUB3 collision
points lie on all nine zero-carrier equations. Their projections are still
distinct because their base coordinates differ.

This argument applies to every subset. Independently, the checker computes
the full symbolic determinant of all five displayed frontier maps, rather
than trusting inheritance alone. The retained carrier derivative block is
unit lower triangular; its finite nilpotent inverse is checked in both
orders, and the full3-by3 Schur determinant is expanded exactly to1.
All retained source variables remain symbolic.

## Complete independent replay

Run

```
.venv/bin/python research/jacobian/cubic-reduction/verify_carrier_subsets12.py
```

The primary producer constructs a recursive source embedding and then
composes the full original map. The standard-library Fraction checker
instead substitutes into the live map one deleted row at a time. It
reconstructs all512 restrictions and count records, checks all512 projected
collisions, compares the five selected maps and embeddings, and verifies
their full polynomial determinants and integral coefficients. It computes
the Pareto set afresh from the complete subset table.

- [Primary producer](carrier_subset_frontier12.py)
- [Complete512-case table](carrier-subset-frontier12/checkpoint.json)
- [Independent checker](verify_carrier_subsets12.py)
- [Independent receipt](carrier-subset-frontier12/validation.json)

The least (degree, count) at each fixed dimension is:
3:(7,13),4:(6,15),5:(5,17),6:(5,18),7:(4,23),8:(4,24),
9:(4,25),10:(4,27),11:(4,29),12:(3,36).
Within this family, degree<=4 therefore requires at least7 variables,
degree<=5 requires at least5, and degree<=6 requires at least4.
There is no cubic below12 variables and no map below13 occurrences here.
These are finite-family conclusions only.

## Comparison and remaining routes

The separate [primary comparison](carrier_subset_comparison.py) exhausts
literal carrier-zero subsets of three saved maps: BCR2, KCUB2's13-variable
42-term cubic, and KCUB3. Its [checkpoint](carrier-subset-comparison/checkpoint.json)
is retained as an experiment, without independent replay or global-record
status. The present frontier dominates their reported literal-subset
frontiers. BCR2 is the original six-cubic41-term frame produced by its
saved helper, not the later40-term conjugate. That later frame, general
carrier combinations, target repairs, and other constructions are outside
this comparison. No literature-wide record or novelty claim is made.

The [weighted graph completion theorem](WEIGHTED_GRAPH_COMPLETION.md)
separately excludes one-stage polynomial cubic repair for every
same-weight linear carrier graph. Ordered/coordinated repairs, nonlinear
or cross-weight graphs, different base frames and different constructions
remain available. Global support and special-form dimension minima remain
open; the strict3D profile partition stays62 excluded/12 remaining and
the global3D degree interval remains[5,7].

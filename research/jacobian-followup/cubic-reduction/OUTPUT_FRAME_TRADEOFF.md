# A31-nonlinear-occurrence cubic and the exact output-frame tradeoff

KCUB30,15 November2026. Let K be the integral twelve-variable cubic in
[KCUB3](PROOF12.md). There is an integral degree3 Keller map H with31
nonlinear monomial occurrences,17 linear occurrences and an explicit
integral collision. Its linear part is unipotent and nonidentity. The total
number of monomial occurrences is48, exactly the total for the saved
identity-linear-part map K, which has36 nonlinear and12 linear occurrences.

This is a separate count without identity-linear-part normalization. It
does not improve the established36 count under that normalization, the
12-variable dimension, the13-occurrence raw3D example, or any global
minimum. The construction is a linear output change of the existing map.

With the source coordinates fixed,31 nonlinear occurrences and48 total
occurrences are both sharp among all invertible linear output changes.
The exact Pareto pairs (linear,nonlinear occurrences) are

| Linear | Nonlinear | Total |
| ---: | ---: | ---: |
|12|36|48|
|13|35|48|
|14|34|48|
|15|33|48|
|16|32|48|
|17|31|48|

Every pair is attained by an integral determinant-one output matrix.
The lower bounds also cover affine output changes, since translations
leave linear/nonlinear counts unchanged and only add constant terms here.
By the imported [KCUB26 target-rigidity theorem](TARGET_CUBIC_RIGIDITY.md),
the same minima and Pareto pairs cover every polynomial target change of
degree<=5 whose composition with this fixed K is again a cubic Keller map:
such a target change must be invertible affine. This last corollary uses
KCUB26 as a separately checked premise; the6175-column computation is not
repeated in the new checker.

## Explicit integral map and collision

Define H=M*K by the simultaneous row formulas

```
H2 = K2-3*K7-12*K10
H3 = K3-4*K9-20*K8
H9 = K9+3*K8,
Hj = Kj for all other rows.
```

The matrix M has diagonal1, with the five off-diagonal entries
M_(2,7)=-3, M_(2,10)=-12, M_(3,9)=-4, M_(3,8)=-20, M_(9,8)=3.
It satisfies (M-I)^3=0 and det M=1. Thus H(0)=0, JH(0)=M and det JH=1.
The expanded nonlinear row counts are

```
(4,8,8,2,1,1,1,1,1,1,1,2), total31.
```

The linear rows are the corresponding matrix rows of M. The points

```
p=(0,1,-8,0,-1,0,0,-1,3,0,0,0)
q=(-1,1,4,-1,-1,-1,4,-1,1,-1,1,0)
```

are distinct and both map to (0,1,0,...,0). The
[standalone coefficient witness](output-frame31/witness.json) supplies
all expanded coefficients, matrices and points, not merely these formulas.
The independent checker recomputes the full symbolic determinant in all12
variables using a unit-triangular carrier-block inverse and its3-by3 Schur
determinant. It checks every coefficient and both collision evaluations.

All six Pareto witnesses arise by these successive row operations, where
the last operation uses the already updated ninth row:

```
row2 <- row2-3*row7
row2 <- row2-12*row10
row3 <- row3-8*row8
row9 <- row9+3*row8
row3 <- row3-4*row9.
```

Each step exchanges one nonlinear occurrence for one linear occurrence,
keeping the total48 and determinant1. Every intermediate map is checked.

## Thirty-one nonlinear occurrences are necessary in every output frame

A monomial coefficient is a column vector in the12-dimensional output
space. Invertible output changes preserve nonzero columns, linear
relations and projective distinctness. The saved K has exactly30 distinct
nonlinear source monomials. Thus every output frame needs at least30
nonlinear occurrences, one per column.

Three of these columns, for x1*x2, x1^2*x3 and x1*x10, lie in the span of
output axes e4,e12 and equal

```
3*e4+e12,  e4,  e12.
```

They span a plane and are pairwise nonparallel. Their total support after
any invertible output change is at least4: support3 would put each on one
axis, forcing either three independent axes or two parallel columns.
The other27 columns each need an occurrence, giving31. The displayed H
attains31. This proof allows arbitrary invertible output mixing, including
mixing all the blocks used below, over every characteristic-zero field.

## Forty-eight total occurrences are necessary

Include the12 linear source columns as well. There are42 distinct source
monomials. Partition these columns according to the following original
output row sets. The column sets are disjoint, so their occurrence lower
bounds add, even when a new output matrix mixes their row spaces.

| Original output rows | Columns | Universal occurrence lower bound |
| --- | ---: | ---: |
|{1}|5|5|
|{2,7,10}|13|15|
|{3,8,9}|13|16|
|{4,12}|5|6|
|{5}|2|2|
|{6}|2|2|
|{11}|2|2|
|Total|42|48|

The one-dimensional blocks need one occurrence per nonzero column.
For rows{4,12}, there are two copies of e4, two copies of e12 and one
3*e4+e12. The three projectively different directions force at least one
extra occurrence beyond the five columns, hence6.

For rows{2,7,10}, nine columns are nonzero multiples of h=e2; the four
other columns are

```
u=e7, v=e10, 3*h+u, 12*h+v.
```

If the transformed h has support>=2, those nine columns already cost18,
and all four other columns are nonzero. If h has support1, each pair
(u,3h+u) and (v,12h+v) has total support at least3. Otherwise each member
would have support1, and their difference could not be a nonzero multiple
of the single h axis while remaining independent of h. The bound is
9+3+3=15 in this case, proving15 universally.

For rows{3,8,9}, nine columns are multiples of h=e3. The four light columns
are

```
u=e8, v=e9, 8*h+u-3*v, 8*h+2*v.
```

Again, support>=2 for h gives at least18+4>16. Suppose h has support1.
After deleting its one supporting output coordinate, the light columns
become ubar,vbar,ubar-3*vbar,2*vbar. They need at least5 occurrences:
if all four had support1, the independent ubar and vbar would occupy
different axes and their difference ubar-3vbar would need support2.
In the deleted heavy coordinate, at least two light columns are nonzero.
Any three of the four original light columns are independent in the
3-space span(h,u,v), as the four exact triple minors confirm. A coordinate
functional nonzero on h cannot therefore annihilate three light columns.
The total is at least9+5+2=16. This proves the table and total48.

These are written rank/support arguments. The independent checker
reconstructs all42 columns, their complete block partition, both ninefold
heavy directions, every light vector, all four nonzero triple minors and
the exceptional projective plane. It does not assume primary factorization
or an untested generic matrix.

## Exact Pareto frontier and normalization boundary

For any invertible output matrix, the linear matrix has rank12, so L>=12.
The preceding arguments give N>=31 and L+N>=48. All six integral pairs
in the table attain these bounds. Every other possible pair is dominated
by one of them, proving the entire two-count Pareto frontier in this
fixed-source output class. Constants from affine translations cannot
improve it.

Normalizing H on the left gives M^-1*H=K and returns36 nonlinear terms.
Normalizing its source gives H(M^-1*X)=M*K(M^-1*X), an identity-linear-part
conjugate with74 nonlinear occurrences. Its exact row counts are

```
(7,12,22,6,10,1,3,6,1,1,2,3), total74.
```

The source-normalized coefficients, inverse matrix, derivative and
transported integral collision are also independently reconstructed.
Neither normalization improves the saved36 bound. Arbitrary simultaneous
source/output changes, nonlinear graph parametrizations, higher target
degrees>=6 involving all outputs, other constructions, and all global
sparsity and special-form dimension minima remain open. KCUB29's13 bound
has a different raw graph-class scope and is not promoted to an affine
source/output invariant here.

## Replay and assurance

Run

```
.venv/bin/python research/jacobian-followup/cubic-reduction/construct_output31.py
.venv/bin/python research/jacobian-followup/cubic-reduction/verify_output31.py
```

The first script uses SymPy to construct the six stages and normalized
conjugate. The independent checker uses only standard-library rational
polynomial and integer/rational matrix arithmetic. It reconstructs stage
coefficients by multiplying their output matrices into the raw input,
checks the final map separately from simultaneous row formulas, computes
the full final polynomial determinant, verifies every collision, rebuilds
the inverse-source conjugate and checks every finite premise of the
written universal lower bounds. The base KCUB3 validation is pinned;
the final31-occurrence determinant and collision are also checked directly.

- [Primary constructor](construct_output31.py)
- [Independent checker](verify_output31.py)
- [Independent receipt](output-frame31/validation.json)

The universal output-count deductions and the KCUB26 corollary are written
mathematical arguments, not formal verification or external review. No
new counterexample mechanism, literature-record priority, global minimum,
dimension reduction or total-term improvement is asserted.

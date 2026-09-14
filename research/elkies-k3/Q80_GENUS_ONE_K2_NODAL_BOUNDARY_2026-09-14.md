# Q80 genus-one k2 is confined to a nodal denominator boundary

On the literal direct11952 alternate-Q80 MW17 parent, every solution of
the three genus-one common-quartic cases `k=2,j=0,1,2` must have the following
reduction at `p=131`, in the original parameter coordinate:

```
q_0(t)=t^2+62*t+88,
D=eta*g*d,   g,d monic quadratic,   g_bar=d_bar=q_0,
(g*d)_bar=q_0^2=t^4+124*t^3+90*t^2+39*t+15.                 (1)
```

Here g is the agreement branch factor and d the disagreement factor.
**Both polynomial abscissas have negative Gauss coefficient valuation.**
Thus every integral coefficient solution, and every other branch reduction,
is excluded over Q_131, with no coefficient-height or denominator bound.

**The nodal denominator boundary (1) remains UNKNOWN.** No k2 allocation
is declared empty. All twelve previously surviving genus-one allocations
remain open; this result sharply restricts three of them. The positive
[MW17 parent with two gains and an infinite quadratic base](CORRELATED_QUADRATIC_GAINS_2026-09-12.md)
remains unconstructed. The earlier genus-zero norm-eight replay gap is
neither used nor upgraded here.

If the covering curve has a rational point, eta must moreover be a square
in Q_131, both negative abscissa valuations are even, and its Jacobian has
`v_131(j)=-4*n` for some positive integer n. In particular the existing
positive-rank control base `Y^2=X^3-4*X+1`, which has good reduction at131,
cannot be this base in any k2 presentation on Q80.

## 1. Hypotheses and the additional degree-two character data

Keep the full [common-quartic chart](COMMON_QUARTIC_BRANCH_STRATA_2026-09-14.md):

```
f(x_i)=x_i^3+A*x_i+B=eta*D_0*r_i^2,
deg x_i<=4, deg r_i<=4, D_0 squarefree of binary degree4,
gcd(D_0,Delta)=1,   x_1!=x_2,   deg gcd(D_0,x_1-x_2)=2.       (2)
```

The genus-one curve `C:w^2=eta*D_0` need not have a rational point for the
main theorem. It is geometrically integral. Its two bisection images have
class `2O+4F`, arithmetic genus5 and normalization C; both miss O. In
particular they are degree-one images of their normalizations, not multiple
images of a rational curve. The characteristic-zero branch fibres are
smooth throughout; reduction to a nodal fibre below does not relax this
hypothesis. Work with primitive integral binary forms until (1) permits
monic normalization, so infinity is retained.

For an inherited section T and a selected cubic root e, let `chi_e(T)` be
the regularized squareclass `X_T-e`, using `3*e^2+A` when `T=(e,0)` and
class1 at a pole. The [degree-two reciprocity proof](Q80_DEGREE_TWO_RECIPROCITY_2026-09-14.md)
supplies these classes at every rational and quadratic residue point,
including the nodal orbit. The new check extracts three further exact
properties of that complete table:

| Complete character calculation | Result |
|---|---|
| Rational cubic-root points |110 nonzero, pairwise distinct17-bit codes|
| Nonrational quadratic base orbits with a zero root code |Only the nodal orbit below|
| Unordered rational root-point pairs with XOR code13412 |Only `((35,114),(75,11))`|
| Quadratic root entries with norm code13412 |Only the nodal double root, with its two multiplicities|

The unique quadratic discriminant factor is q_0; the other irreducible
factor degrees are3 and19. All factors are simple. In `F_131[i]`, `i^2=-1`,
the nodal row at `t=100+31*i` has

```
roots       (122+38*i, 122+38*i, 18+55*i),
norm codes  (13412,   13412,     0).                         (3)
```

As polynomials modulo q_0, the double and simple roots are respectively

```
N(t)=29+128*t,        S(t)=73+6*t.                           (4)
```

The checker proves the discriminant factors irreducible and replays the
underlying full character table through its earlier independent checker.
It does not extrapolate from the fully split rows alone.

At either disagreement branch, let e_3 be the unused cubic root. The
anti-invariant difference `P_1-P_2` specializes to e_3 there and to O at
the two agreement branches. Its Kummer representative has odd divisor
only at the complementary roots over d. Therefore reciprocity gives

```
product_over_d Norm(chi_e3(T)) = 1 modulo squares            (5)
```

for every inherited T. This is one character condition, unlike the two
conditions in the earlier genus-zero k0 gate. Each original P_i also
satisfies its character product over all four branches.

For split branch fields, (5) and injectivity of the110 rational codes
force the two unused-root reductions to coincide. For an unramified
quadratic field with nonrational base residue, (3) is the only zero code.
A quadratic field with rational base residue gives a double base reduction;
ramified quadratic fields have rational residue too. Consequently d reduces
either to a double rational base point or to q_0. The nodal case selects
the simple root as e_3 and the double root for both P_i.

## 2. A denominator lemma for genus one

Normalize `v_p(eta)=epsilon` in `{0,1}` and D_0 to be primitive integral.
If an abscissa has valuation `-m<0`, write
`x=p^(-m)*X`, `r=p^(-n)*R` with primitive integral coefficients.
Comparing Gauss valuations in (2) gives `2n=3m+epsilon`, and reduction gives

```
X_bar^3=eta_unit,bar*D_0,bar*R_bar^2.                        (6)
```

Every geometric root of D_0,bar is a root of X_bar. The primitive divisor
`p^m*x_coordinate-X(t)` specializes, as an actual image cycle, to

```
2O + sum_a ord_a(X_bar,h)*F_a.                              (7)
```

Properness of stable maps supplies a genus-one stable-map limit after a
finite DVR extension; proper pushforward identifies its image cycle with
(7). The imported theorem is
[Abramovich--Oort, Theorem2.8 and section2.5](https://arxiv.org/pdf/math/9808074).
This use is the genus-one analogue of the earlier
[genus-zero integrality argument](Q80_GENUS_ZERO_COLLISION_CLOSURE_2026-09-14.md#2-denominators-would-force-an-elliptic-component),
but a genus-one domain can have one elliptic component, so an additional
argument is necessary.

Two distinct smooth fibres in (7) would require two positive-genus domain
components, impossible in arithmetic genus one. Since D_0 factors as two
binary quadratics and the bad-fibre degrees are2,3,19, its bad roots can
only come from q_0. If D_0,bar has at most one smooth geometric root, that
root is rational and its only possible forms, up to a unit, are

```
L_b^4,        L_b^2*q_0,        q_0^2,
L_b=T-b*Z,   L_infinity=Z.                                  (8)
```

For `L_b^2*q_0`, equation (6) forces `X_bar,h` proportional to
`L_b^2*q_0`. The covering curve has a geometric nodal special model
`w^2=L_b^2*q_0`, with rational normalization. Resolving the nodal model
after finite base change adds rational chains; its semistable components
are all rational. Any stable-map model only adds rational components.
It therefore cannot map nonconstantly to the smooth elliptic F_b in (7).
This excludes the middle case, including ramified constant scalars after
adjoining their square roots.

For `L_b^4`, equation (6) forces even multiplicities in X_bar,h. Any second
smooth root is already impossible; a remaining bad quadratic with odd
multiplicity is incompatible with (6). The degree-four form is therefore
proportional to `L_b^4`, and (7) is `2O+4F_b`.

Its unique elliptic domain component would map to F_b with degree four.
Since `p=131` does not divide four, this is an unramified isogeny after
choosing origins. It has four distinct preimages of `O intersect F_b`.
Every such preimage must be a node adjoining a component mapped into O:
the generic image misses O, so the zero locus of the pulled-back Cartier
section for O is entirely vertical and cannot have an isolated zero on
this elliptic component. A rational tree attached there cannot be entirely
contracted in an unmarked stable map; some component of that tree must map
nonconstantly to O, since no rational component can dominate F_b.

The dual graph is a tree once the elliptic component uses the full genus
budget. Hence the four attachment points require four separate rational
trees with positive total mapping degree to O. But (7) has total O degree
only two. This contradiction excludes `L_b^4` as well.

Thus **a negative abscissa valuation is possible only when
`D_0,bar` is proportional to `q_0^2`**. This is a statement about the actual
specialized image and its domain; numerical divisor classes alone would
not prove it. The all-nodal last case of (8) has only rational image
components, so the argument deliberately leaves it open.

## 3. Complete classification with integral abscissas

Suppose both x_i are integral. Then the r_i are integral too. If epsilon
were one, reduction of (2) would give a rational2-torsion section on the
24I1 parent. Its height formula rules this out, so eta is a unit.

If d reduces to `L_b^2`, integrality forces each x_i to select the same
root at both nearby branch points. The two selected roots are distinct
and rational modulo p, so b is one of the sixteen fully split fibres in
the retained table, and the individual norm contribution over d is zero.
The agreement contribution over g must therefore have zero code. The
complete zero-code classification in section1 says that g reduces either
to `L_c^2` with the same selected root at c, or to q_0 with selected root S.

If d reduces to q_0, both x_i select N there, with norm code13412.
The agreement contribution must have that code as well. The complete
13412 classification gives either g reducing to q_0 with root N, or

```
g_bar=(t-35)*(t-75),  x_i(35)=114,  x_i(75)=11.              (9)
```

These exhaust the possibilities, including quadratic branch fields whose
base residues are rational, ramified fields, and infinity. In the
`D_0,bar=q_0^2` case, any integral x_i would reduce to a polynomial section
on the constant twist of E, with ordinate `q_0*r_i,bar`, meeting an I1 node
at N. A section cannot meet that node: the local equation `xy=t` prevents
both section coordinates from vanishing. Thus an integral abscissa is
impossible in this case.

The remaining three types are closed by the following complete finite
coefficient gates.

### Two doubled smooth fibres

Here `D_0,bar=L_b^2*L_c^2`, and b cannot equal c: disagreement at b and
agreement at c would contradict one another after reduction. The reductions
are height-four polynomial points on a constant twist of the parent, with
ordinate `L_b*L_c*r_i,bar`. A nonsquare unit twist has arithmetic rank at
most one by the retained Frobenius certificate, so cannot give distinct
abscissas at b. For the square unit, use the full reduced MW lattice and
the complete1,313-point norm-four roster from the earlier proof.

Exactly **21 unordered section pairs** have a common ordinate zero where
their root choices agree and another where those choices differ. This
enumeration includes the pair whose agreement point is infinity. For each
pair the full simultaneous coefficient system has26 equations in24
variables: five coefficients in each x_i and r_i, and four free coefficients
of the normalized common quartic.

At every pair, the Jacobian modulo131 has rank24, while adjoining the
first error vector from the literal A,B modulo `131^2` gives rank25.
Thus **none lifts even modulo `131^2`**. The errors use the actual rational
parent coefficients; rank of the Jacobian alone is not the obstruction.
At infinity the coefficient of `t^2` is fixed instead of that of `t^4`.
Sign changes of either ordinate just negate its r_i columns and do not
change solvability.

### A doubled smooth fibre and the nodal quadratic

Here `D_0,bar=L_b^2*q_0`, b is a fully split rational fibre, and the
agreement root is S from (4). Every possible x has the unique form

```
x=S+q_0*(z_0+z_1*t+z_2*t^2),      z_i in F_131.             (10)
```

The complete `131^3=2,248,091`-tuple census asks whether `f(x)/q_0` is a
square up to an arbitrary nonzero scalar. There are exactly two abscissas:

| x coefficients, low to high | Scalar | Monic square-root coefficients | Rational square-root zeros |
|---|---:|---|---|
|`61,1,98,17,97`|122|`13,70,57,21,76,1`|none|
|`123,76,86,33,105`|47|`3,32,91,34,32,1`|70|

Neither square root vanishes at any of the sixteen fully split sites.
Thus neither can have the required factor L_b. These are actual finite
section identities, but neither is a solution of this k2 contact problem.

### Four distinct branch residues

The only remaining quartic is `D_0,bar=q_0*(t-35)*(t-75)`. Conditions
`x=N modulo q_0` and (9) determine a unique degree-below-four polynomial
x_0. Every possible x is `x_0+c*D_0,bar`, with `c in F_131`.
For all131 values the exact quotient `f(x)/D_0,bar` fails to be a square
up to scalar. This closes the last integral case.

Together with section2, these calculations prove (1) and show that both
abscissas must have denominators. They do not close the all-nodal boundary.

## 4. The remaining valuation and rational-base conditions

In (1), both g and d are monic integral quadratics with unramified quadratic
root field over Q_131. For an abscissa of valuation `-m_i`, (6) forces

```
X_i,bar=lambda_i*q_0^2,     R_i,bar=mu_i*q_0^2,
lambda_i^3=eta_unit,bar*mu_i^2,
m_i>0,     m_i=epsilon modulo2.                            (11)
```

Since q_0 has no projective F_131 root, the primitive homogeneous value
of g*d is always a unit square modulo131 at a Q_131 base point. Consequently

```
C(Q_131) is nonempty if and only if eta is a Q_131 square.  (12)
```

Sufficiency follows already from either point at infinity after monic
normalization. A rational base point therefore forces epsilon zero and
square unit residue. Locally absorb its square root: then
`m_i=2*a_i>=2`, `x_i=131^(-2*a_i)*X_i`,
`r_i=131^(-3*a_i)*R_i`, and each lambda_i is a nonzero square in F_131.
This local normalization does not discard the global literal scalar.

For a monic quartic `t^4+b*t^3+c*t^2+d*t+e`, put `I=12*e-3*b*d+c^2`.
The [classical binary-quartic invariant formulas, Fisher section7.1](https://arxiv.org/pdf/math/0610318)
give `j=256*I^3/disc(D_0)`. Here `I(q_0^2)=60 modulo131` is a unit.
The discriminants of g and d are units, and

```
disc(g*d)=disc(g)*disc(d)*Res(g,d)^2,
v_131(Res(g,d))=2*n,  n>=1,
v_131(j(C))=-4*n.                                           (13)
```

The middle equality is the norm valuation in the unramified quadratic
field generated by a root of g; g and d are distinct because D_0 is
squarefree. Thus this base has potentially multiplicative reduction.
The retained control `Y^2=X^3-4*X+1` has discriminant `16*229`, a131-unit,
and cannot be isomorphic to such a base over Q. This does not exclude a
different positive-rank genus-one curve satisfying (11)--(13).

## 5. Reproducible evidence and the next gate

The [frozen input](../artifacts/generated-results/elkies-k3-q80-genus-one-two-disagreement-v1/input.json)
binds the literal parent, complete earlier tables, both new implementations
and preserved discovery outputs. The
[Sage producer](scripts/certify_q80_genus_one_k2_nodal_boundary.sage) checks
the characters, all131 interpolants and all21 full lift systems. Its
[C++ census](scripts/certify_q80_nodal_quadratic_census.cpp) checks all tuples
in (10), reconstructing square roots from their leading coefficient.
Sage CPU was0.158 seconds and census CPU0.279 seconds, under respective
40-second and30-second caps; the producer memory cap was4GiB.

The [independent Python checker](scripts/verify_q80_genus_one_k2_nodal_boundary.py)
replays the complete prior arithmetic, proves the degree2,3,19 factors
irreducible, reconstructs all21 matrices from exact rational coefficients,
and repeats the whole nodal census with integer NumPy arrays and square
roots reconstructed from constant coefficients. Its
[receipt](../artifacts/generated-results/elkies-k3-q80-genus-one-two-disagreement-v1/independent-replay.json)
records PASS in2.657 seconds. Seven targeted squareclass, degree,
irreducibility and inconsistent-lift controls pass.

```
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 research/elkies-k3/scripts/verify_q80_genus_one_k2_nodal_boundary.py
python3 -m unittest discover -s research/tests -p 'test_q80_genus_one_k2_nodal_boundary.py' -q
```

One discovery script used an unavailable polynomial CRT method and stopped
before its enumeration. That source is retained; explicit interpolation
resolved the startup issue before freezing. No frozen run failed and no
frozen source was retagged. The written stable-map, parity and local-base
arguments are not formally verified. No older norm-eight assurance changes.

The next k2 gate is the actual characteristic-zero coefficient system near
(1), with the valuation pattern (11), the literal scalar and smooth generic
branch fibres retained. Higher coefficient denominators cannot be discarded
as a bounded-search failure. The other nine genus-one allocations and other
actual MW17 parents remain separate construction routes. No rational
coefficient lift, positive-rank base satisfying the full section system, or
infinite rank-at-least19 source is supplied by this theorem.

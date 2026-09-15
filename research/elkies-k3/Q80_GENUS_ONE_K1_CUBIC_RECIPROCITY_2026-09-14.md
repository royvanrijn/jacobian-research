# Cubic reciprocity excludes the Q80 genus-one k1 chart

On the literal direct11952 alternate-Q80 MW17 parent, all four polynomial
genus-one common-quartic allocations `k=1,j=0,1,2,3` are empty over Q_131,
and hence over Q. This includes every coefficient denominator, literal
scalar, higher contact, infinity and unpointed genus-one base in the chart.

The five `k=0,j=0,...,4` genus-one allocations remain open. Other actual
MW17 parents and higher-height rational abscissas remain possible. The
[MW17 parent with two gains on an infinite quadratic base](CORRELATED_QUADRATIC_GAINS_2026-09-12.md)
is still unconstructed. This proof does not use or upgrade the older
genus-zero norm-eight singular-polynomial replay assurance.

## 1. One rational agreement point and a cubic disagreement divisor

Keep the full [common-quartic chart](COMMON_QUARTIC_BRANCH_STRATA_2026-09-14.md):

```
f(t,x_i)=x_i^3+A(t)*x_i+B(t)=eta*D_0(t)*r_i(t)^2,
deg x_i<=4, deg r_i<=4,
D_0 squarefree of binary degree4, gcd(D_0,Delta)=1,
x_1!=x_2,       deg gcd(D_0,x_1-x_2)=1.                      (1)
```

Thus `D_0=g*d` up to a constant, where g is the linear agreement factor
and d is the cubic disagreement factor. Normalize binary D_0 to be
primitive integral over Z_131 and absorb its scalar into eta, with
`v_131(eta)=epsilon` in `{0,1}`. Neither g nor d is required to be monic
in the original affine chart; their reductions and all local arguments
below are homogeneous.

The g point is Q_131-rational and has a rational projective reduction b.
Every such parent fibre is smooth at131. The reduced discriminant has
simple irreducible factors of degrees2,3,19; write

```
q_0=t^2+62*t+88,
q_3=t^3+6*t^2+81*t+64.                                      (2)
```

These complete fibre facts are independently rechecked from the literal
rational parent. In particular, the agreement residue is never nodal.
The two bisection images have class `2O+4F`, miss O and are birational
images of the genus-one covering curve. These facts follow from (1),
including when infinity is a branch or an abscissa has degree below four.

## 2. Coefficient integrality for k1

Suppose one abscissa has negative Gauss valuation. Write
`x=p^(-m)*X` and `r=p^(-n)*R` with primitive integral coefficients and
m positive. Since A,B are integral, comparing Gauss valuations and then
reducing gives

```
2*n=3*m+epsilon,
X_bar^3=eta_unit,bar*D_0,bar*R_bar^2.                        (3)
```

Every geometric root of D_0,bar is a root of the primitive binary form
X_bar,h. As in the [genus-one denominator argument](Q80_GENUS_ONE_K2_NODAL_BOUNDARY_2026-09-14.md#2-a-denominator-lemma-for-genus-one),
properness of stable maps identifies the actual specialized image cycle as

```
2O + sum_a ord_a(X_bar,h)*F_a.                              (4)
```

The imported properness statement is
[Abramovich--Oort, Theorem2.8 and section2.5](https://arxiv.org/pdf/math/9808074).
Two distinct smooth fibres in (4) would require two positive-genus domain
components, which a connected genus-one domain cannot supply. Since the
rational agreement residue b is smooth, every other root of D_0,bar must
be bad. The complete bad-fibre degrees2,3,19 leave exactly three possible
branch reductions, up to a unit:

```
L_b^4,         L_b^2*q_0,         L_b*q_3,
L_b=T-b*Z,    L_infinity=Z.                                 (5)
```

The first two reductions are excluded by the same written subarguments
in the earlier denominator proof; neither subargument assumes that D_0
factors as two quadratics over Q_131.

For `L_b^2*q_0`, equation (3) forces X_bar,h proportional to that form.
The cover has only rational geometric semistable components: its reduced
nodal model has rational normalization, and semistable resolution adds
rational chains. Scalar twists are harmless after finite extension. Such
components cannot dominate the smooth elliptic fibre F_b in (4).

For `L_b^4`, parity in (3) and the degree-four budget force
`X_bar,h` proportional to `L_b^4`. Its image cycle is `2O+4F_b`.
The unique elliptic domain component must map to F_b by a degree-four
unramified isogeny, with four distinct preimages of `O intersect F_b`.
Each preimage must attach to a component mapped into O: the pulled-back
Cartier section of O has no generic zeros, so its special zeros on the
elliptic component cannot be isolated smooth points of the total domain.
The dual graph is a tree after the elliptic component uses the genus
budget. Each of the four attached trees must contain a component mapping
nonconstantly to O; a wholly contracted unmarked rational tree is unstable,
and no rational component dominates F_b. Total degree at least four onto
O contradicts the degree two in (4). Here p=131 does not divide four.

It remains to exclude `L_b*q_3`. Equation (3) now forces X_bar,h
proportional to `L_b*q_3`, so its zero at b is simple. Pass to a finite
unramified extension that splits the smooth cubic fibre at b, and use
minimal local coordinates if b is infinity. Its three distinct cubic
roots lift to analytic integral functions `rho_1,rho_2,rho_3` on the residue
disc. For each j, the primitive series `p^m*(x-rho_j)` has reduction of
order one. By Weierstrass preparation it has exactly one simple zero in
that disc. The three zeros are distinct because `rho_j-rho_k` is a unit.

Consequently `f(t,x)=product_j(x-rho_j)` has three distinct simple zeros
there. In (1), each odd zero must be a branch zero of D_0. But D_0 has
exactly one zero in this residue disc, since its reduction is `L_b*q_3`.
This is impossible. The preparation statement is the classical complete-DVR
theorem in [Berger, section1 and Corollary1.2](https://perso.ens-lyon.fr/laurent.berger/articles/article33.pdf).
No ordinary-contact assumption is inserted: simplicity follows from the
exact local reduction order.

All possibilities (5) are excluded. Thus **both abscissas are integral**.
Then both r_i are integral as well. If epsilon were one, reduction of
(1) would give a rational2-torsion section on the reduced24I1 parent.
Its height is `4+2*(T.O)>0` for every nonzero section, so such torsion is
impossible. Hence eta is a unit. The same observation prevents an r_i
from reducing identically to zero.

For the remaining k0 chart, the same argument gives a useful necessary
condition. Without the rational agreement point there may be no smooth
root of D_0,bar at all; the only additional degree-four possibility is
`q_0^2`. Thus **any negative abscissa coefficient valuation in the full
genus-one polynomial chart requires D_0,bar proportional to q_0^2**.
That k0 denominator boundary is not excluded here. Its branch field need
not factor as two quadratic fields over Q_131, so the earlier k2 orientation
proof cannot be copied without rechecking that hypothesis.

## 3. The agreement residue cannot collide with disagreement

At the agreement point, both abscissas reduce to the same rational cubic
root e at b. If a disagreement point reduced to b, the integral polynomial
abscissas would have this same residue there as well. The parent fibre is
smooth, so Hensel lifting gives a unique cubic root in that residue class.
They would therefore be equal at that disagreement point, a contradiction.
Thus g_bar and d_bar are coprime, including at infinity.

At a rational residue c of d, the two abscissa residues are distinct
rational roots of its smooth cubic. The fibre must therefore split fully
over F_131. Their unused root has the same rational residue for every
branch point reducing to c, even if its local field is ramified.

For any inherited section T, let chi_e(T) be the regularized class
`X_T-e`, with the derivative `3*e^2+A` at a selected2-torsion zero
and squareclass1 at a pole. The [reciprocity argument](Q80_DEGREE_TWO_RECIPROCITY_2026-09-14.md#1-coupled-norm-conditions-including-zero-ordinates-and-poles)
applies to both original sections and their difference. Let C denote
their17-bit reduced norm code. At the d branches, e_1,e_2 are the roots
chosen by the two sections and e_3 is unused. Necessarily

```
C_d(e_3)=0,
C_d(e_1)=C_d(e_2)=C_b(e),                                  (6)
```

where C_d is the XOR of norm codes over all d factors. At each fibre the
three codes have XOR zero. The complete retained rational roster has
110 distinct nonzero codes C_b(e). The classes specialize on the good
cubic-root-curve model, including nodal residues, zeros and poles.

If d_bar is not squarefree, its degree-three reduction is either
`L_c^3` or `L_c^2*L_a` with distinct rational c,a. The unused-root norm
then equals respectively C_c(e_3) or C_a(e_3): the even multiplicity
cancels and the remaining odd multiplicity does not. Both are nonzero.
This includes ramification in the branch field, since a norm of an
unramified unit class acquires its total local degree as exponent.
Hence **d_bar is squarefree**.

Its three possible factor-degree partitions are now `1+1+1`, `1+2`
and `3`. Only the simultaneous conditions (6), not a single unused-root
zero alone, can close these cases.

## 4. Complete finite norm gates

For the `1+1+1` partition, use the sixteen fully split rational fibres,
three distinct sites and all ordered root choices. All120,960 combinations
fail even the unused-root zero condition.

For `1+2`, combine one of those rational fibres with every nonrational
quadratic base orbit in the retained complete table, preserving the nodal
root multiplicities. There are838,944 ordered choices. Exactly twelve
have zero unused-root norm, but none matches an allowed rational agreement
code. The checker independently finds them by code-indexed joins.

For partition `3`, use

```
F_(131^3)=F_131[theta]/(theta^3+theta+3),
encode(a+b*theta+c*theta^2)=a+131*b+131^2*c.
```

The defining cubic has no F_131 root and is irreducible. All
`(131^3-131)/3=749320` nonrational cubic base orbits are checked, selecting
the least encoded Frobenius conjugate. Exactly125,595 have a smooth cubic
splitting fully over this field; one is nodal. A smooth nonsplit cubic
cannot supply the two distinct selected roots, so no other residue fibre
is eligible.

The unused-root zero condition survives at exactly eight orbits. The
other two root codes coincide, as (6) requires, but **none is one of the
110 rational agreement codes**:

| Cubic base representative | Fibre | Other two norm codes |
|---:|---|---:|
|235367|smooth|89988|
|295766|smooth|65341|
|374594|smooth|76440|
|484763|smooth|87669|
|532125|smooth|51968|
|1050615|nodal|68381|
|1154864|smooth|108754|
|1236534|smooth|123639|

At the nodal row, the simple root has zero code and the double root has
code68381 twice. Coalescing selected roots are retained with multiplicity;
the calculation does not omit singular reductions or require their
characteristic-zero branches to be singular.

Thus no factor-degree partition satisfies (6). Together with sections2--3,
this proves the k1 exclusion over Q_131 for all four j allocations.
The argument needs neither a global nor a local rational point on the
covering genus-one curve.

## 5. Frozen evidence and assurance

The [frozen input](../artifacts/generated-results/elkies-k3-q80-genus-one-k1-reciprocity-v1/input.json)
binds the literal parent, prior complete norm tables, both new census
implementations, the driver, checker and tests. The
[Sage driver](scripts/certify_q80_genus_one_k1_reciprocity.sage) checks
the complete discriminant factorization and composite norms, reduces
the literal basis coordinates and independently verifies all eight special
cubic rows in Sage's extension field.

The [C++ producer](scripts/certify_q80_cubic_norm_census.cpp) uses a direct
`X^(131^3) modulo f` splitting test, explicit degree-three multiplication,
prime-field splitter constants and Frobenius-product norms. The
[independent C++ implementation](scripts/replay_q80_cubic_norm_census.cpp)
uses generic polynomial division, three successive p-Frobenius substitutions,
extension-field splitter constants and determinant norms. The
[Python checker](scripts/verify_q80_genus_one_k1_reciprocity.py) reconstructs
the coefficient input with integer/Fraction arithmetic, verifies all bad
factor degrees, repeats the composite calculation by indexed joins and runs
the complete independent cubic census. The census totals and every
zero-norm row agree exactly. No retained older complete census is repeated.

```
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 research/elkies-k3/scripts/verify_q80_genus_one_k1_reciprocity.py
python3 -m unittest discover -s research/tests -p 'test_q80_genus_one_k1_reciprocity.py' -q
```

Every driver has a30-CPU-second/2GiB limit, each census a90-CPU-second/512MiB
limit, and factorization has at most64 declared splitter trials before
failing closed. The complete discovery runs needed at most8 and10 trials.
All nine targeted failure/control tests pass. The original bounded previews
and complete discovery outputs are preserved; none failed. The final
[execution](../artifacts/generated-results/elkies-k3-q80-genus-one-k1-reciprocity-v1/execution.json)
and [independent replay](../artifacts/generated-results/elkies-k3-q80-genus-one-k1-reciprocity-v1/independent-replay.json)
record the actual resource use.

Independent replay certifies the finite arithmetic. The stable-map
integrality, preparation and reciprocity deductions are written proofs,
not formal verification or external review. The required infinite
rank-at-least19 family is not supplied by this exclusion.

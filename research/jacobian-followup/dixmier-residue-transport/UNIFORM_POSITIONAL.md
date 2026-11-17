# Uniform positional quantum residue beyond the quartic test

14 November 2026. **QRT3.** The QRT2 cocycle has period one on a single
correction category stable under positional polynomial Weyl shears of
arbitrary degree. The proof uses boundary-weight decomposition and finite
coefficient identities over Q(a), with every exceptional integer weight
addressed. It does not infer an all-degree theorem from a finite degree sweep.

The quartic test requested before this extension passes independently:
6125 retained variables, rank 4780, a 1345-dimensional projected lower-lift
family, and 661 current columns. All linear and quadratic dependence
vanishes. Quintic and sextic producer tests also pass and motivated the
uniform calculation; the theorem below supplies their infinite-degree
coverage.

## The category without growing Bernstein cutoffs

Keep the retained seed S,T, the derivation delta, and the 16-term functional
Lambda_plus from [QRT2](COMBINED_COCYCLE.md). All coefficients below belong
to B=Q[X,X^-1,Q], with finite support but no fixed pole order. Declare

\[
s=s_3Z^3+s_2Z^2+s_1Z+s_0,\qquad
t=t_2Z^2+t_1Z+t_0,
\]

where deg_(X,Q)(s_3)<=16 and deg_(X,Q)(t_2)<=15. Here the degree is max(i+j)
on finite Laurent monomials X^i Q^j. The remaining lower-lift coefficients
are unrestricted elements of B. At order four allow

\[
u=u_1Z+u_0,\quad\deg(u_1)\le18,\quad u_0\in B;
\qquad v=v_0,\quad\deg(v_0)\le17.
\tag{1}
\]

This is the increasing union of the degree categories with Z weight w>=4:
lower bounds (3w+16,3),(2w+15,2) and current bounds (w+18,1),(17,0).
Every finite tuple in (1) belongs to one such category. The top-coefficient
bounds are fixed throughout; no unrestricted correction theorem is asserted.

We prove

\[
\Lambda_+d_1(u,v)=0,\qquad\Lambda_+O_5(s,t)=1
\tag{2}
\]

for every compatible lower lift in this category. Thus it has no order-five
lift in the declared parity branch. Odd-hbar corrections, higher Z orders,
other classical seeds and general DC2 remain outside the result.

## Boundary-weight decomposition

Let nu(X,Q,Z)=(1,-1,-2). Exact source identities give

\[
\nu(S)=-2,\quad\nu(T)=-1,\quad
\nu(\delta)=1,\quad\nu(\partial_Z)=2.
\]

The source S,T are homogeneous for this weight, not merely bounded by it.
Every monomial seen by Lambda_plus has weight 12. Put U=XQ and W=X^2Z.
A homogeneous lower-lift sector can be written

\[
s_a=X^a f(U,W),\qquad t_a=X^{a+1}g(U,W),
\tag{3}
\]

with f,g polynomial in U,W and respective W orders at most three,two.
Finite Laurent X support makes the sector decomposition finite for every
actual tuple. The highest-W coefficients obey

\[
a+2j\le10
\tag{4}
\]

for their U^j terms. For a>10 both highest coefficients vanish, so the
orders drop to two,one.

The lower equation has weight a+2; its right side has weight six.
Only a=4 is inhomogeneous. The Moyal-linear part of the period sees
only that sector, while the bilinear period couples only sectors

\[
a+b=8.
\tag{5}
\]

These identities replace an unbounded degree sweep by a central affine
sector, seven finite weight pairs, and one symbolic complementary pair.

## Finite coefficient projection for every integer weight

The Poisson differential decreases an input Q exponent by at most one:
its two monomial output exponents are j+j_seed and j+j_seed-1. Therefore
equations with output Q exponent at most nine depend only on lower inputs
with Q exponent at most ten. Inputs with higher Q exponent cannot alter
this projection. The scalar bilinear period sees Q exponent at most seven
in either input; its Moyal-linear part sees at most nine.

Retain j=0,...,10 in (3) and project the equation to output j<=9. Every
actual homogeneous solution maps to this finite system and all its period
coefficients are retained. The projected systems are supersets of actual
solutions; extending every projected solution to a full solution is unnecessary.
This argument places no bound on a or on the source X pole order.

For a>10 the high-weight system has 55 variables, generic rank 50 and
kernel dimension five. At b=8-a, ignore (4) to obtain a larger low-weight
system with 77 variables, generic rank 70 and kernel dimension seven.
Over Q(a), the cross-period matrix vanishes on these two kernels exactly.

Specialization is certified by full-rank minors of the original matrices.
Up to nonzero rational constants, their determinants are

\[
M_{50}=a^5(a-9)(a-8)^2(a-7)^4
\prod_{r=1}^{6}(a-r)^5\,(a+1)^3(a+2),
\]

\[
M_{70}=a^2(a-11)(a-10)^3(a-9)^5
\prod_{r=3}^{8}(a-r)^7\,(a-2)^6(a-1)^4(a+1).
\tag{6}
\]

Denominators used to clear the kernel bases have no roots outside the
roots of these minors. For every integer a>11, both generic kernels
specialize to complete kernels and the cross-period identity remains valid.
The only exceptional high integer weight is a=11. Its actual systems,
including (4) at the complementary weight -3, have zero cross period
by a separate exact rational calculation.

The remaining pairs are a=-2,-1,0,1,2,3,4 with b=8-a. Their actual
top restrictions (4) are imposed. All seven restricted bilinear matrices
vanish on the complete projected kernels. Every pair of integer weights
in (5) occurs in this list or has a member greater than ten, so this covers
all quadratic lower-lift dependence.

## The affine sector and the constant period

Use the actual retained polynomial order-three lift as an origin in sector
a=4. Its period is exactly one. The projected homogeneous kernel in this
sector has dimension five. Exact calculation annihilates the linear
variation on all five directions. The a=4 bilinear check above annihilates
the quadratic variation. Homogeneity makes every other sector invisible
to the linear part. Hence the period is one on every compatible lower lift.

Current corrections admit a similarly finite check. The necessary weights
are nu(u)=10 and nu(v)=11. For unrestricted u_0 only the eight monomials
X^(j+10)Q^j, j=0,...,7, can meet the period. For u_1Z only
X^(j+12)Q^jZ, j=0,...,3, meet both the weight and degree conditions.
For v_0 only X^(j+11)Q^j, j=0,...,3, are relevant. Their sixteen exact
pairings vanish. This proves the current-order assertion in (2), including
every finite Laurent coefficient.

## Arbitrary positional polynomial Hamiltonians

For any polynomial H(X,Q), QRT1 gives the exact quantum automorphism
alpha_H=exp(-ad_star(H)/hbar), with classical substitution
Z -> Z+delta H. No degree bound on H is now needed.

Classical substitution preserves the highest Z coefficients in (3) and
(1), while all lower coefficients it changes are unrestricted. The induced
order-two seed shifts delta^3(H)*S_ZZZ/24 and its T counterpart have Z
orders at most two,one, respectively. They therefore affect only unrestricted
lower coefficients. At order four, the lower-lift shift involves s_ZZZ,
which is scalar, and t_ZZZ=0. The seed fifth derivative is scalar for S
and zero for T; both sixth derivatives vanish. Thus every induced quantum
shift changes only the unrestricted u_0 coefficient. The inverse shear
has the same properties.

The whole category is consequently stable under exact positional polynomial
shears of arbitrary degree. Lambda_(+,H)=Lambda_plus composed with the
inverse classical substitution is a strong period-one cocycle for every
such changed seed, including its combined native and transported corrections.
Finite words of these positional shears are covered because their potentials
add. This concludes the uniform extension of the quartic test.

## Independent algebraic certificate

[uniform_weight_probe.py](uniform_weight_probe.py) constructs the generic
systems, kernel bases, original-matrix minors, finite exceptions and affine
checks. [uniform_weight_probe.json](uniform_weight_probe.json) retains
both rational and denominator-cleared polynomial bases.

[check_uniform_independent.py](check_uniform_independent.py) imports none
of the producer, legacy algebra, SymPy or FLINT. It reconstructs the
matrices with Fraction polynomial arithmetic. It checks every polynomial
kernel identity, generic basis independence and the zero cross-period
identity. Separate integer Bareiss determinants at 122 integer points
verify the two displayed determinant polynomials: matrix entries are affine
in a, so rank+1 exact evaluations certify each polynomial identity. Complete
division by the listed linear factors certifies all roots. The finite
weight systems, affine variations and current pairings are independently
recomputed by a separate Gaussian elimination and monomial bracket engine.

[uniform_independent.json](uniform_independent.json) records the passed
finite algebraic certificate. The grading, projection relevance and
infinite-class transport above are written proofs, not formal verification.
No independent external mathematical review or novelty priority is claimed.

## First mixed-polarization gate

The exact Weyl shear X -> X-p_x^2, with p_x fixed, is an elementary shear
in the Fourier-conjugate polarization. It leaves the old fiber chart:
for R=2X-3X^2Q,

\[
[X-p_x^2,R]=-4\hbar p_x+12\hbar XQp_x+6\hbar^2Q\ne0.
\]

Every element of the old R-centralizer fiber algebra commutes with R.
Hence the new position does not belong to that algebra. The
[exact gate calculation](mixed_polarization_gate.json) identifies what a
native comparison must add: the transverse momentum omitted from the
three-variable fiber complex. This is not a failure of the Weyl
automorphism or a proof of mixed-polarization nonliftability.

The next problem is therefore a genuinely four-generator correction
complex with a declared native category after this polarization change,
including the transverse direction and all six relations. More positional
degree tests are superseded by QRT3.

## Reproduction

```sh
.venv/bin/python research/jacobian-followup/dixmier-residue-transport/verify_degree_category.py --degree 4
.venv/bin/python research/jacobian-followup/dixmier-residue-transport/check_degree_category_independent.py --degree 4
.venv/bin/python research/jacobian-followup/dixmier-residue-transport/uniform_weight_probe.py
.venv/bin/python research/jacobian-followup/dixmier-residue-transport/check_uniform_independent.py
```

The quartic and uniform inputs, output hashes and coverage are retained in
their receipts. These are short exact calculations with no search campaign.

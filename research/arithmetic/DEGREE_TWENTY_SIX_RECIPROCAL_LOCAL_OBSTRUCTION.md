# Degree-26 reciprocal pencils at the prime three

The mixed-sign reciprocal problem has exact solutions over F3, but none
has a lift modulo nine preserving the reciprocal constraints. A complete
finite census covers both nonzero reciprocal parameters. Consequently a
normalized integral reciprocal degree-26 model cannot have good paired
reduction at three.

The conclusion concerns integrality and nondegenerate reduction in the
specified normalization. It does not exclude rational models, Q3 models
with bad reduction, or nonreciprocal pairings.

## Setting and good reduction

Write N(u) monic of degree thirteen, with N(1)=N'(1)=0, and put

\[
 2uN'(u)-N(u)=(u-1)S(u),\qquad f(z)=N(z^2)/z.                   \tag{1}
\]

The residual critical polynomial S has degree twelve and leading
coefficient 25. Reciprocal pairing uses u->lambda/u with lambda in the base field:
lambda in F3* for the finite census, and lambda in Q3* for the
integral characteristic-zero consequence. In the
fixed-point-free case its coefficient conditions are

\[
 S_j=\lambda^{6-j}S_{12-j}\quad(0\le j\le5),                   \tag{2}
\]

with S_12=25. In particular a=N(0)=S_0=25lambda^6.

Good paired reduction means that modulo three, S is squarefree,
S(0)S(1)!=0, S has no common root with u^2-lambda, N is coprime
to S, and M=N^2/u has precisely six distinct nonzero values at the
S-roots, one for each reciprocal pair. These conditions ensure exactly
thirteen distinct finite branch values, each with two simple critical
points. Since a is a unit, the constant-term identity forces lambda
to be a unit as well.

## Theorem 1. Finite-field feasibility and its complete census

In the setting over F3, with lambda nonzero, there are exactly two
nondegenerate reciprocal paired models. Both have lambda=1.
Their ascending coefficient lists, from constant through leading term,
are

\[
 \begin{split}
 N^{(1)}:\ &(1,1,1,1,0,0,1,1,0,0,1,0,1,1),\\
 N^{(2)}:\ &(1,2,0,2,0,2,0,0,2,0,2,1,2,1).
 \end{split}
\]

The first is the explicit degree-26 pencil

\[
 H(z)=z^{26}+z^{24}+z^{20}+z^{14}+z^{12}+z^6+z^4+z^2,\quad a=1,
\]
\[
 Q_6(y)=y^6-y^5+y^4-y^3-y^2-1,
\]
\[
 \boxed{\operatorname{Disc}_z(H(z)-sz+1)
                        =\bigl(sQ_6(s^2)\bigr)^2
                      \quad\text{in }\mathbb F_3[s].}           \tag{3}
\]

Its Q_6 is squarefree and nonzero at zero. Both models have all
twenty-six finite critical points simple, the desired thirteen
double-transposition branch values, and geometric and arithmetic
monodromy A26. The opposite-sign reciprocal factors have degrees
(8,4) in the first model and (6,6) in the second.

**Complete finite reduction.** The unknown coefficients of N are
N_0,...,N_12, since N_13=1. Synthetic division determines

\[
                  S_j=\sum_{k=j+1}^{13}(2k-1)N_k.
\]

The two zero-branch equations and the six equations (2) are affine
linear in these thirteen unknowns for fixed lambda in F3*.
Row reduction gives:

| lambda | Linear rank | Free coefficients | Complete affine candidates |
|---:|---:|---|---:|
| 1 | 6 | 5,7,8,9,10,11,12 | 2,187 |
| 2 | 8 | 5,7,8,10,11 | 243 |

The checker enumerates every value of these free coefficients.
For lambda=1, exactly 728 candidates pass the squarefreeness,
fixed-point and coprimality tests; exactly two have invariant M,
and both have the required six distinct values. For lambda=2 none
passes all the nondegeneracy conditions.

The value-count test uses the etale algebra A=F3[u]/S.
The reciprocal involution is fixed-point-free on its twelve geometric
points, so its invariant algebra has dimension six. The checker tests
M=lambda-reciprocal and rank six of the seven vectors 1,M,...,M^6.
This is exactly that M separates the six reciprocal pairs; N coprime
to S makes all its values nonzero.

For each of the two surviving models the checker additionally verifies
the complete characteristic polynomial chi_M=Q_6^2, squarefreeness
of Q_6, and the full discriminant identity, including (3).
These are exact polynomial calculations, not evaluations at sampled s.

For monodromy, N(z^2)-sz is irreducible over F3bar(s), since it is
linear in s and N(0)=1. The poles have orders 25 and one. Their
ramification is tame, and the 25-cycle forces primitivity.
The finite critical points have ramification index two, so inertia
contains a double transposition. The primitive double-transposition
theorem gives A26, while the square discriminant gives the reverse
inclusion. The arithmetic groups also lie in A26 by their verified square identities and contain
their geometric groups, so are A26 and their splitting fields have constant
field F3. The precise group theorem is retained in the
[existing primary-source ledger](alternating_hurwitz_priority_ledger.json). ∎

No characteristic-zero lifting assertion follows from this finite-field
existence theorem.

## Theorem 2. Neither finite-field model lifts modulo nine

Neither of the two models in Theorem 1 admits a monic lift N over Z/9
with the normalized double zero, reciprocal residual critical polynomial
and reciprocal paired critical values, even allowing lambda to vary.

**Proof by exact linear obstruction.** For lambda=1 modulo three,
every unit lift of lambda modulo nine is a square. Choose rho with
rho^2=lambda and rho=1 modulo three. Use fourteen variables

\[
                     x=(N_0,\ldots,N_{12},\rho).
\]

There are twenty necessary polynomial constraints over Z3:

1. N(1)=0 and N'(1)=0;
2. the six reciprocal coefficient equations
   C_j=S_j-rho^(12-2j)S_(12-j)=0 for j=0,...,5;
3. twelve remainder coefficients requiring S to divide
   \[
   E(u)=\rho^2u^{24}N(u)^2-
                      \bigl(u^{13}N(\rho^2/u)\bigr)^2.          \tag{4}
   \]

The leading coefficient of S is 25, a unit modulo nine, so these
remainders are defined over Z3. Equation (4) is precisely equality
of M(u) and M(lambda/u) at the critical roots, after multiplying
by units.

For the lift x_0 with coefficients represented by 0,1,2 and rho=1,
let R=F(x_0)/3 modulo three and J be the 20-by-14 constraint
Jacobian modulo three. A lift x_0+3delta would have to satisfy

\[
                              J\delta=-R.
\]

For both models rank J=13, whereas rank(J|-R)=14.
The retained certificate supplies explicit left-null vectors as well:
in the constraint order just listed, the first uses
2C_1+C_2; the second uses C_1+2C_2+C_4+2C_5.
Their derivatives annihilate all fourteen columns of J at their
respective finite points, while their product with -R is one.
Thus the linear equations are inconsistent. ∎

The Jacobian is computed exactly modulo nine by

\[
 J_{\bullet j}
 =\frac{F(x_0+3e_j)-F(x_0)}3\pmod3.
\]

The binomial expansion proves this formula. The checker verifies the
left-null products directly; the conclusion does not depend solely
on a matrix-rank output.

## Corollary 3. A good-reduction obstruction over Z3

There is no normalized monic N in Z3[u] satisfying (1), reciprocal
paired-critical equality, and the good-reduction conditions above.

**Proof.** Its unit residue a forces lambda to be a unit. If lambda
reduces to two, Theorem 1 excludes its nondegenerate reduction.
If it reduces to one, that reduction must be one of the two listed
models. Every unit lambda congruent to one has a square root modulo
nine reducing to one, and Theorem 2 excludes the resulting modulo-nine
lift. ∎

This forces bad reduction at three for any integral normalized rational
reciprocal solution in degree twenty-six. It is not rational
nonexistence. The known degree-ten reciprocal model has rho=2/15,
already illustrating that denominators or nonunit reciprocal parameters
can occur at three. A proposed degree-26 model may likewise fall
outside this good-reduction chart.

## Exact evidence

The [checker](verify_degree_twenty_six_reciprocal_local.py) and
[receipt](degree_twenty_six_reciprocal_local_verification.json) retain
both complete affine enumerations, both model coefficients, both
full discriminant identities, and the exact modulo-nine obstruction
vectors. The computation uses Python's integer arithmetic and pinned
SymPy 1.14.0.

~~~sh
.venv/bin/python research/arithmetic/verify_degree_twenty_six_reciprocal_local.py
~~~

The finite domain and its coverage are the affine reductions above.
The characteristic-zero consequence uses only reduction modulo nine
and the stated integrality/nondegeneracy hypotheses.
Unrestricted tame lifting, nonreciprocal models, and reciprocal models
with bad reduction are not ruled out or constructed here.

This is internal exact verification. No formal verification, independent
external replay, external human review or novelty claim is made.


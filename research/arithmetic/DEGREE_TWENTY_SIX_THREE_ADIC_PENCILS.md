# Degree-26 alternating pencils over the three-adic field

There is a five-parameter family of monic even degree-26 pencils over Q3
with literal square discriminant and full regular A26 monodromy.
The coefficients lie in Z3 and have the good reduction constructed over F3.
Their critical pairs are preserved by an involution of the etale critical
algebra; this involution is not induced by u->lambda/u over Q3.

This establishes positive local existence while the
[reciprocal method has a good-reduction obstruction](DEGREE_TWENTY_SIX_RECIPROCAL_LOCAL_OBSTRUCTION.md).
No rational descent to Q or full-A26 specialization over Q3 is asserted.

## Theorem 1. A five-parameter integral family

Put

\[
 N_0(u)=u^{13}+u^{12}+u^{10}+u^7+u^6+u^3+u^2+u+1
                                       \quad\text{over }\mathbb F_3.
\]

For every independent choice of coefficients N_7,N_9,N_10,N_11,N_12
in Z3 reducing to their values in N_0, there are unique coefficients
N_0,...,N_6,N_8 in the corresponding residue classes such that the monic
polynomial N of degree thirteen has:

* N(1)=N'(1)=0;
* twenty-six simple finite critical points for f(z)=N(z^2)/z;
* precisely thirteen finite branch values, with two critical points over each.

Writing H(z)=N(z^2)-N(0) and a=N(0), one has

\[
 \operatorname{Disc}_z(H(z)-sz+a)
             =\bigl(5^{25}sQ_6(s^2)\bigr)^2                   \tag{1}
\]

for a monic squarefree Q_6 in Z3[y] with unit constant and discriminant.
The geometric and arithmetic groups over Q3(s) are A26, and the
splitting field is regular over Q3.

The uniqueness assertion is in the residue neighborhood of the specified
critical-pairing branch, as constructed below.

## The critical algebra and its fixed involution

Use the exactly reciprocal integer lift

\[
 T(t)=t^{12}+t^9+t^8+t^7+2t^6+t^5+t^4+t^3+1,\qquad
                         A=\mathbb Z_3[t]/T.
\]

Modulo three this is the squarefree residual critical polynomial of N_0.
Consequently A is a finite free, etale Z3-algebra. Its constant coefficient
one makes t a unit, and reciprocity gives the literal involution

\[
                              \tau(t)=t^{-1}.
\]

It has no fixed geometric point modulo three. Since two is a unit,
A splits into its plus and minus eigenspaces, each free of rank six.
The first six coefficient coordinates in the basis 1,t,...,t^11
restrict to an isomorphism from A^- onto Z3^6: their reduction has
rank six, as checked exactly in the retained certificate.

For N in the residue neighborhood of N_0, put R_N(u)=2uN'(u)-N(u).
There is a unique U=U(N) in A, congruent to t modulo three, with

\[
                              R_N(U)=0.                       \tag{2}
\]

Indeed at reduction R_(N0)=(u-1)T, and the derivative at t is
(t-1)T'(t), a unit. Newton iteration in the complete finite free
algebra A gives the root (2), varying integrally and analytically
with the coefficients of N.

Define the eight equations

\[
 G(N)=\left(N(1),N'(1),
       \text{first six coordinates of }(1-\tau)
                                      \frac{N(U(N))^2}{U(N)}\right).
                                                                    \tag{3}
\]

Their reduction vanishes at N_0. The six algebra coordinates vanish
if and only if the entire anti-invariant part vanishes.

## The exact rank calculation and Newton proof

At a critical root, for M_N(u)=N(u)^2/u,

\[
                    M_N'(u)=\frac{N(u)R_N(u)}{u^2}.
\]

Thus critical-root motion makes no contribution to the derivative of (3).
For delta N of degree at most twelve,

\[
 dG(\delta N)=
 \left(\delta N(1),\delta N'(1),
       \text{first six coordinates of }(1-\tau)
                            \frac{2N_0(t)\delta N(t)}t\right)
                                         \quad\text{modulo }3. \tag{4}
\]

In coefficient columns indexed 0,1,2,3,4,5,6,8, its exact minor,
using entries represented by 0,1,2, is

\[
 \begin{pmatrix}
 1&1&1&1&1&1&1&1\\
 0&1&2&0&1&2&0&2\\
 2&1&2&2&1&0&2&0\\
 0&2&2&0&2&2&2&2\\
 2&1&0&2&1&1&2&2\\
 0&1&0&0&2&1&1&1\\
 2&1&2&2&1&1&2&0\\
 0&2&1&2&2&1&0&0
 \end{pmatrix}.
\]

Its integer determinant is 16, hence a unit modulo three. Fix the other
five coefficients as in the theorem. Newton iteration in the selected
eight coefficients converges to the unique solution of (3) in their
specified residue neighborhood.

Here the needed Hensel assertion follows directly from the iteration.
The selected Jacobian minor remains invertible over Z3 throughout this residue ball.
If the current error belongs to 3^m, its Newton correction belongs to
3^m; the integral Taylor remainder makes the next error belong to 3^(2m).
The corrections converge in Z3 and solve all eight equations.
The same Taylor argument applied to two solutions with equal residue,
and their unit Jacobian, proves uniqueness. The root construction (2)
uses the identical argument in the finite free complete algebra A.
Inverses of derivative units expand in convergent integral power series,
so the equations (3) have the required integral Taylor expansions.

This is a written Newton proof, not an application of an uninspected
general tame-Hurwitz lifting theorem.

## Passport, discriminant and monodromy

At a solution, N(1)=N'(1)=0. Good reduction ensures its zero at one has
exact order two, and all remaining critical points are simple and nonzero.
The twelve critical u-roots represented by U(N) remain etale. Equation
(3) pairs their M-values by tau. Modulo three those values are six
different nonzero values, so their pairwise differences and values remain
units. There are exactly six nonzero squared branch values and the zero
branch, giving thirteen double-transposition levels for f.

The six powers of M_N(U) are a basis of A^+ modulo three and hence over
Z3. Thus M_N(U) generates A^+, whose characteristic polynomial is
the monic degree-six Q_6 in (1). Its constant and discriminant stay units.
The even-degree critical-value resultant identity has leading factor

\[
                              25^{25}=5^{50}.
\]

It therefore proves the literal identity (1). This identity concerns
the pencil over Q3(s).

The extension is transitive because N(z^2)-sz is irreducible over
Q3bar(s). The pole of order twenty-five supplies a 25-cycle and forces
primitivity. Finite inertia supplies a double transposition, so the
retained primitive double-transposition theorem gives full geometric A26.
Identity (1) puts the arithmetic group in A26. It contains the geometric
group, hence both equal A26 and their constant-field quotient is trivial.

Finally all coefficients have the good reduction in the earlier local
obstruction theorem. If the critical pairing were u->lambda/u for
lambda in Q3*, that theorem would exclude the model. The pairing furnished
by tau on the etale algebra is therefore genuinely outside that reciprocal
ansatz. No conclusion about rational Q coefficients follows.

## Lemma 2. Uniform tangent rank

Let k be a field of characteristic different from two, let b=2h+1,
and suppose L=2b-1=4h+1 is nonzero in k. Let N be monic of degree b
with N(1)=N'(1)=0. Assume

\[
 2uN'-N=(u-1)S
\]

has b distinct nonzero roots: one at 1 and the other 2h at roots
alpha_1,...,alpha_(2h) of S. Suppose N is nonzero at those S-roots,
and M=N^2/u takes h distinct values there, twice each.
Let iota exchange the two roots in each equal-value pair.

On the b-dimensional space of variations delta N of degree at most
b-1, the linear map

\[
 \delta N\longmapsto
 \left(\delta N(1),\delta N'(1),
                  (1-\iota)\frac{2N\,\delta N}{u}
                         \text{ on the roots of }S\right)      \tag{5}
\]

has rank exactly h+2. The last target summand has dimension h.

**Proof.** Extend k to an algebraic closure; linear rank is unchanged.
Evaluation at the b distinct points 1,alpha_1,...,alpha_(2h)
is a Vandermonde isomorphism on polynomials of degree at most b-1.
The condition delta N(1) and the h weighted pair differences are
independent: the pairs are disjoint and their weights 2N(alpha)/alpha
are nonzero. Their joint rank is h+1.

Consider the Euler variation

\[
                         \delta N=uN'-bN.
\]

Its leading term cancels, including when the characteristic divides b.
It has delta N(1)=0. At a critical S-root, 2uN'=N, so

\[
                       \delta M=(1-2b)M,
\]

which is invariant under iota. But delta N'(1)=N''(1) is nonzero:
the simple critical root at one gives
S(1)=2N''(1)!=0. Hence the derivative coordinate is independent of
the first h+1 conditions. This proves rank h+2. ∎

The special determinant above is an exact control of this general rank
statement. It is not an exceptional generic-rank assumption.

## Corollary 3. Uniform local lifting from a good finite-field model

Let k be a finite field of odd characteristic p with p not dividing
L=4h+1. Suppose a monic polynomial Nbar in k[u] of degree 2h+1
satisfies

\[
 \overline N(1)=\overline N'(1)=0,\qquad
 2u\overline N'-\overline N=(u-1)\overline S,
\]

where Sbar is squarefree of degree 2h, Sbar(0)Sbar(1)!=0,
Nbar is coprime to Sbar, and Nbar^2/u has exactly h distinct
nonzero values at the Sbar-roots, each twice.

Let O be the ring of integers of the unramified p-adic field K with
residue field k. Then this model lifts to an (h-1)-parameter analytic
family of monic N in O[u], with the same normalized paired-simple-critical
passport in characteristic zero. Its pencils satisfy

\[
 \operatorname{Disc}_z(N(z^2)-sz)
                      =L^L\bigl(sQ_h(s^2)\bigr)^2,             \tag{6}
\]

where Q_h is monic of degree h and has unit constant and discriminant.
For n=4h+2>=10, their geometric group is A_n. Their arithmetic group
is A_n with regular splitting field if L is a square in K; otherwise
it is S_n and the constant field is K(sqrt(L)).
In particular the admissible degrees n=m^2+1 have regular A_n over K.

This criterion assumes the finite-field model. It does not construct
such models in every degree or prove rational descent.

**Proof.** The equal-value pairing is unique and Galois invariant.
It therefore defines a fixed-point-free involution iota_0 on the
finite etale algebra A_0=k[u]/Sbar. Choose any monic lift T of Sbar/L
and put A=O[t]/T. This is finite free and etale over O.
The involution lifts by elementary Newton iteration: solve T(v)=0
in A with v reducing to iota_0(t), and define tau(t)=v.
The derivative is a unit because Tbar is squarefree.
The resulting endomorphism reduces to an automorphism, hence is
invertible; alternatively tau^2(t) is a T-root reducing to t and
Newton uniqueness makes it exactly t. Thus tau is an involution.
Since two is a unit, A^- is free of rank h.

For every nearby N, solve R_N(U)=0 in A with U reducing to t.
The derivative at the initial residue root is
L(t-1)Tbar'(t), a unit. Impose N(1)=N'(1)=0 and the h independent
coordinates of (1-tau)N(U)^2/U. Lemma 2 gives rank h+2.
Choose a unit coefficient minor and fix the other h-1 coefficients
arbitrarily in their prescribed residue classes.
The precision-doubling Newton argument already given above proves
existence and uniqueness of the remaining coefficients over O.

All critical-point and critical-value nondegeneracies are units
at reduction, so persist. The invariant algebra A^+ has rank h;
M=N(U)^2/U generates it by reduction and Nakayama.
Its characteristic polynomial gives Q_h and the resultant identity
gives (6).

For n>=10, transitivity, the L-cycle at infinity, and finite
double-transposition inertia give geometric A_n exactly as above.
The discriminant square class in K(s) is L. If L is square, the
arithmetic group is A_n as well. If it is nonsquare, the arithmetic
group is S_n, whose alternating fixed field is the constant
quadratic extension K(s)(sqrt(L)); the geometric group is A_n, so
there is no further constant extension. ∎

## Exact controls

The [checker](verify_degree_twenty_six_three_adic.py) and
[receipt](degree_twenty_six_three_adic_verification.json) recompute (4),
the rank-six anti-invariant projection and the displayed determinant.
They also perform the first full Newton step modulo nine, producing
ascending numerator coefficients

\[
                    (7,1,4,4,6,0,7,1,3,0,1,0,1,1).
\]

The exact critical root in A/9 and all paired-value equations are verified.
The resulting residual critical polynomial is not reciprocal for any unit
lambda modulo nine. This verifies a concrete first approximation; the
infinite family follows from the written Newton argument.

~~~sh
.venv/bin/python research/arithmetic/verify_degree_twenty_six_three_adic.py
~~~

The checker pins and reuses the finite-field polynomial arithmetic in
the [local certificate](verify_degree_twenty_six_reciprocal_local.py).
The full finite-field passport and the imported group theorem are recorded
in the [local note](DEGREE_TWENTY_SIX_RECIPROCAL_LOCAL_OBSTRUCTION.md);
the critical-value resultant identity is retained in the
[all-degree bridge](ALTERNATING_HURWITZ_DISCRIMINANT_BRIDGE.md).

This is internal exact checking with a complete scoped lifting proof.
Formal verification, independent external replay and external human
review are not claimed.


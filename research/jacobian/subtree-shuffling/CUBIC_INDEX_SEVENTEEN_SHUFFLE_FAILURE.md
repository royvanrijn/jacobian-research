# A cubic index-seventeen collision and failure of uniform subtree shuffling

There is an explicit homogeneous cubic H over Q in eighteen variables such
that
\[
 (JH)^{17}=0,\qquad (JH)^{16}\ne0,
\]
while Id+H identifies two distinct rational points. Consequently eventual
cubic constant-vector shuffle membership fails at path length17 and at
every larger path length. The same fixed lengths admit triangular polynomial
automorphisms with unbounded inverse degrees as the dimension grows.
The fixed-length-four question remains open.

The witness is the eighteen-variable homogeneous companion already used in
the retained cubic-linear construction. The contribution here is its exact
nilpotency-index check and the precise shuffle/uniformity consequence.
No new original Keller mechanism, minimum failing length, historical
priority, external review or formal verification is claimed.

## The literal witness and its exact certificate

Use coordinates (X,W,tau)=(x_1,...,x_12,w_1,...,w_5,tau).
Let K=X+Q+C be the retained twelve-variable map at parameter scale11=2;
Q and C are its homogeneous quadratic and cubic parts. Its last seven
cubic components are zero. Put
\[
 H(X,W,\tau)=
 \bigl(\tau Q(X)-\tau^2\iota W,\ C_1(X),...,C_5(X),\ 0\bigr),             \tag{S1}
\]
where iota includes the five W coordinates in the first five positions of
the twelve-component top block.

The accompanying literal input gives all45 nonzero cubic monomials of H,
their coefficients, and these two source points:
\[
\begin{split}
 a={}&(0,1,-8,0,-1,0,0,-1,0,0,0,0,\ 0,0,0,0,-1,\ 1),\\
 b={}&(-1,1,4,-1,-1,-1,4,-1,-1,-1,-2,0,\ -5,12,0,-4,-1,\ 1).
\end{split}
\]
They are distinct and satisfy
\[
 (I+H)(a)=(I+H)(b)=(0,1,0,0,0,0,0,0,0,0,0,0,\ 0,0,0,0,0,\ 1).         \tag{S2}
\]
No floating-point calculation or generic-fibre inference is needed.

**Proposition ExactCubicIndexSeventeenCollision.**
All entries of (JH)^17 are zero as rational polynomials. The entry in row
x_4 and column tau of (JH)^16 is exactly
\[
 288x_1^{10}x_2^4x_3^3\tau^{15}.                                      \tag{S3}
\]
In particular the nilpotency index is exactly17. Since JH is nilpotent,
det(I+JH)=1. Thus I+H is an ordinary polynomial Keller map with the
collision (S2), over every characteristic-zero field.

The standard-library checker differentiates the45 literal cubics and
multiplies the full sparse polynomial matrix successively through power17.
At power16 there are five nonzero entries with ten monomials in total;
(S3) is the entire specified entry. At power17 there are none.
It evaluates every component at both rational points.

It also reconstructs every literal coefficient from the pinned
twelve-variable JSON. Diagonal conjugation x_11 -> x_11/2 changes its
recorded scale11=1 to scale11=2; separating degrees two and three and using
(S1) then reproduces all45 cubics. The input source points are transformed
by the inverse diagonal, then lifted by W=-C(X), tau=1. They reproduce
the two displayed points exactly. This provenance reconstruction is
additional evidence: the direct polynomial powers and collision alone
already establish the proposition.

All identities are over Q, so their zero identities and the nonzero rational
coefficient288 persist under any characteristic-zero extension. The
determinant assertion also follows algebraically: I+JH has polynomial
inverse sum_(j=0)^16(-JH)^j; its determinant is a polynomial unit and
takes value1 at the origin.

## Failure at a specific shuffle length

For fixed p let E_p mean that the constant vector belongs to the cubic
length-p shuffle span at every sufficiently large internal grade.
Use the all-substitutions polarized quotient A_(3,p) from
[the exact subtree-shuffle quotient](SUBTREE_SHUFFLE_AND_INVERSE_DEGREE.md), and let U_k be its weighted inverse
coefficient in internal grade k.

**Theorem CubicShuffleFailureAtSeventeen.**
The classes U_k in A_(3,17) are nonzero at arbitrarily large grades.
In particular E_17 is false. The same statements hold for every p>=17.

Indeed F=I+H has a unique formal inverse at the origin, since its linear
part is the identity. That formal inverse cannot be a polynomial: a
polynomial inverse would make its formal identity with F a global
polynomial identity and would contradict (S2). Thus infinitely many
of its homogeneous inverse terms are nonzero.

Polarize -H into a symmetric ternary product. At every vector a its
multiplication operator is -JH(a)/3, so its seventeenth power is zero.
Polarization in the repeated argument gives precisely the defining
all-substitutions relations of A_(3,17). Consequently each inverse term
is an evaluation of U_k. A universal class that were zero could not have
a nonzero evaluation. The infinitely many nonzero inverse terms therefore
give arbitrarily large k with U_k nonzero. The retained orbit reduction
then makes E_17 false.

For p>=17 the same algebra also satisfies the pth-power identity, and
the same evaluation proves arbitrarily large nonzero U_k in A_(3,p).
Equivalently, the class of maps with JH^17=0 is contained in each of
the later nilpotency classes. There is no implication here from
failure at17 to failure at a smaller length such as4.

## An explicit jet family in linearly growing dimension

The following construction gives much smaller coordinates than a free-algebra
truncation. It is a coefficientwise restriction of scalars over the finite
algebra K[t]/(t^(M+1)); its proof is included directly.

**Lemma FixedSeedJetInverseDegreeCriterion.**
Let H be any homogeneous cubic map on K^n, and let
\[
 G_t(x)=\sum_{j\ge0}t^jG_j(x)
\]
be the formal inverse of F_t(x)=x+tH(x). Then G_0=x, and each nonzero
G_j is homogeneous of ordinary degree2j+1. Write
Y(t)=sum_(j=0)^M y_j t^j and define the polynomial map on n(M+1) coordinates
by the coefficients of
\[
 {\cal F}_M(Y)=Y(t)+tH(Y(t))\pmod {t^{M+1}}.                           \tag{S5}
\]
It is a triangular polynomial automorphism with homogeneous cubic
nonlinear part. Its inverse degree is exactly
\[
 D_M=\max\{2j+1:0\le j\le M,\ G_j\ne0\}.                             \tag{S6}
\]
Consequently the fixed seed I+H has a polynomial inverse if and only if
the degrees of its explicit triangular jet family are uniformly bounded.

The output in grade j is y_j+[t^(j-1)]H(Y), so only strictly lower grades
occur in its nonlinear part. Formal inversion and substitution of the
target series V(t) show that its inverse is G_t(V(t)) modulo t^(M+1).
A coefficient of t^ell G_j(V(t)) has ordinary total degree2j+1 in the
independent target coefficients. This proves the upper bound in (S6).
Setting all positive-grade target inputs to zero exhibits G_j(v_0)
in output grade j, and its degree2j+1 cannot cancel with terms of other
ordinary degrees. This proves equality.

If the D_M are bounded, only finitely many G_j are nonzero, and
specializing t=1 gives a polynomial inverse. Conversely a polynomial
inverse of I+H has only finitely many homogeneous components, which by
the unique formal inverse identity are precisely the G_j. This proves
both directions of the fixed-seed criterion.

If (JH)^p=0, the jet's nonlinear Jacobian has pth power zero. Indeed,
on a tangent series it is multiplication by tJH(Y(t)), up to the
irrelevant sign if the map is written as I minus its nonlinear part.
Its pth power is t^p(JH(Y(t)))^p=0 in the finite scalar algebra.
This uses matrix multiplication, not differentiation of its coefficients.

If JH has exact indexp and M>=p-1, the jet has exact indexp as well:
the block from input grade0 to output gradep-1 in its (p-1)st Jacobian
power is (JH(y_0))^(p-1), a nonzero polynomial matrix.

**Theorem ExplicitUnboundedTriangularJetsAtIndexSeventeen.**
For every M>=16, apply (S5) to the literal eighteen-variable seed.
The resulting map has dimension18(M+1), homogeneous cubic nonlinear
part of exact Jacobian nilpotency index17, and
\[
 2\left\lfloor\frac{M+2}{3}\right\rfloor+1
 \ \le\ \deg {\cal F}_M^{-1}\ \le\ 2M+1 .                            \tag{S7}
\]
In particular its inverse degree tends to infinity.

The upper bound and exact index follow from the lemma. For the lower
bound, the seed's noninjectivity makes G_t an infinite series. If all
G_j in an interval K<=j<=3K-2 were zero, the cubic recursion
\[
 G_j=-\sum_{a+b+c=j-1}\mu(G_a,G_b,G_c)\qquad(j>=1)
\]
would make every later G_j zero: once j>=3K-1 the child-grade sum
j-1 exceeds3(K-1), so some child lies in the already zero tail.
This contradicts nontermination. Hence every such interval contains
a nonzero coefficient. Take K=floor((M+2)/3), so 3K-2<=M, and use (S6).
For M>=16 this also meets the exact-index requirement.

For every fixed p>17 one can obtain exact indexp by direct sum with the
independent cubic chain on p coordinates
\[
 T_{q_0}=0,\quad T_{q_1}=q_0^3,\quad
 T_{q_j}=q_0^2q_{j-1}\quad(2<=j<=p-1).
\]
Its Jacobian has exact indexp, and I-T has inverse degree2p-1.
Thus the direct-sum family has dimension18(M+1)+p, exact indexp,
and inverse degree max(D_M,2p-1), which is unbounded with M.
There is no claim of a minimum dimension or smallest failing index.

## The finite graded reverse implication

For comparison, the retained quotient gives the reverse bridge needed
for the shuffle equivalences below. If U_k is nonzero in A_(3,p),
with k>=p-1, truncate at grade k:
\[
 B_k=A_{3,p}\big/\bigoplus_{j>k}(A_{3,p})_j .
\]
This is a finite-dimensional algebra, because there are finitely many
tree shapes in each grade. The discarded grades form an ideal since
grafting increases the sum of internal grades by one. Its cubic
multiplication map T(y)=mu(y,y,y) gives a triangular I-T with inverse
degree exactly2k+1. The defining length-p relations have grade at
leastp, so the gradep-1 comb survives and gives exact Jacobian indexp.
These are the all-substitutions identities of the retained quotient.

This supplies a second proof of qualitative unboundedness for every
p>=17. The jet family gives the explicit linear dimension and degree
bounds (S7), independently of computing large quotient bases.

## Correct directions and the remaining boundary

The retained orbit reduction and finite graded construction give the
genuine equivalences, for fixed d,p:
\[
 \begin{gathered}
 \text{eventual constant-vector membership}\\
 \Longleftrightarrow\
 \text{a dimension-independent inverse bound for every degree-d map
 I-H with (JH)^p=0}\\
 \Longleftrightarrow\
 \text{a dimension-independent inverse bound even restricted to
 polynomial automorphisms in that class}.
 \end{gathered}                                                        \tag{S4}
\]
The middle assertion includes the existence of polynomial inverses.
For the last-to-first direction, any arbitrarily large nonzero U_k
has the finite graded triangular realization just proved.

A noninjective homogeneous cubic witness with JH^p=0 therefore implies
failure of these equivalent properties. The reverse implication to a
noninjective map is not asserted: failure of a uniform bound is already
witnessed by triangular automorphisms in increasing dimensions.

The accepted cube-zero theorem proves E_p for p<=3. The present witness
proves failure for all p>=17. By monotonicity there is a smallest failing
cubic length, and it lies in[4,17]. This statement does not classify
lengths4through16 or decide whether length4 has a finite upper bound.

The original three-variable degree-seven Keller map is nonhomogeneous.
Its constant determinant alone would not justify nilpotency of its
nonlinear Jacobian or a specific shuffle length. The explicit homogeneous
companion (S1) and the literal identity (S3) are essential for the fixed17
claim. The twelve-variable mixed quadratic/cubic map by itself would not
supply that claim either. The argument concerns cubic arity; it does not
settle a binary fixed-index shuffle question.

## Primary context and reproduction

Bisi, Dyszewski, Gantert, Johnston, Prochno and Schmid,
Random planar trees and the Jacobian conjecture, arXiv2301.08221v3,
Conjecture1.3, Lemma3.3 and Theorem1.4 state the eventual quantifiers,
the fixed-index shuffle identity and the implication from the all-length
assertion to the Jacobian conjecture. Their primary PDF is already
retained under subtree-shuffling/sources/. The quotient and reverse
finite-dimensional realization used above are proved in the retained
ShuffleOrbitReduction and ShuffleInverseDegreeEquivalence.

The polynomial mechanism is Alpöge's. The eighteen-variable companion is
from the already retained economical de Bondt construction
(arXiv2610.24889v1, Sections4.6--4.7), and is used in [the retained cubic-linear construction](../cubic-reduction/POWER89.md).
The separate Thompson24 normalization and its reported index18 are not
premises of this proof.

Run:

    python3 research/jacobian/subtree-shuffling/verify_cubic_index_seventeen_shuffle_failure.py \
      --repo-root . \
      --output /tmp/cubic_index_seventeen_review.json

The literal input, mandatory source ledger and exact receipt accompany
this proof. A coefficient perturbation is rejected by a nonzero
seventeenth power at a rational point; malformed support/degree/index
controls are also included. The finite certificate proves the explicit
polynomial identities and collision. The arbitrarily-high-grade and
unbounded-family conclusions follow from the written argument, not
from scanning a finite inverse prefix.

Updated 20 November 2026.

# Bounded full identity lifts can have unbounded inverses

14 November 2026. **QNI1.** A common finite bound on the forward
images of a formal Weyl automorphism does not imply a common
bound on its inverse. There are exact full examples in every
rank n>=3, agreeing with the identity to any prescribed finite
parameter order, with forward degree at most12 independent
of that order.

This is a rescaling of the retained Keller counterexample,
not a new Dixmier counterexample or a rank-two construction.
It rules out a universal bounded gauge comparison derived
only from formal equivalence and bounded forward images.
In rank two, even the cotangent subclass of that comparison
would imply JC2.

## An all-degree construction from any Keller map

Let K have characteristic zero and let
\[
 f: {\mathbb A}^n\longrightarrow{\mathbb A}^n,\qquad
 \det Df=1,\quad f(0)=0,\quad Df(0)=I.
\]
Any Keller map can be put in this form by affine changes.
For an independent parameter t set
\[
                    f_t(x)=t^{-1}f(tx),\qquad
                    A_t(x)=(Df_t(x))^{-1}.          \tag{1}
\]
Since f has no constant term, f_t is polynomial in t,x.
The inverse Jacobian is polynomial as well, by the adjugate
formula, and
\[
 f_0=x,\qquad A_0=I,\qquad Df_t(x)=Df(tx).           \tag{2}
\]
If deg f=d, then
\[
             \deg_x f_t\le d,\qquad
             \deg_x A_t\le(n-1)(d-1).              \tag{3}
\]
The scalar parameter t contributes no spatial degree.

In the finite PBW algebra over K[h,t] with [p_i,x_j]=h*delta_ij,
define
\[
 X_i=f_{t,i}(x),\qquad
 P_i=\sum_j(A_t)_{ji}(x)p_j.                         \tag{4}
\]
These satisfy every Weyl relation exactly:
\[
 [X_i,X_j]=[P_i,P_j]=0,\qquad [P_i,X_j]=h\delta_{ij}. \tag{5}
\]
Indeed the polynomial vector fields
delta_i=sum_j(A_t)_(ji) partial_(x_j) are dual to df_(t,i).
They satisfy delta_i(f_(t,j))=delta_ij and commute.
Their commutator kills all f_(t,j), and the invertible
Jacobian makes such a vector field zero. First-order
differential-operator commutators are exactly these vector
field brackets. Thus no quantum correction or infinite
recursion is needed for (5).

The resulting endomorphism Phi_(h,t) has common spatial bound
\[
              D=\max(d,(n-1)(d-1)+1).              \tag{6}
\]
If f is not invertible, its generic cotangent Weyl lift is
not surjective. This classical transfer is recalled, for
example, in Bavula,
[Holonomic modules and 1-generation in the Jacobian Conjecture](https://www.numdam.org/item/10.5802/crmath.556.pdf).
The next argument supplies the exact implication needed here.

## Differential order detects generic nonsurjectivity

Work over any characteristic-zero field in which h,t are
nonzero. Differential order in p is preserved by Phi:
the leading p-symbol is substituted through the invertible
matrix A_t(x), and coefficient substitution through f_t
is injective because its Jacobian is a unit. Therefore
\[
 \operatorname{ord}_p\Phi(Q)=\operatorname{ord}_p Q,\qquad
 \Phi(A_n)\cap K[x_1,\ldots,x_n]=K[f_{t,1},\ldots,f_{t,n}].
                                                               \tag{7}
\]
All fields in (7) can be extended to include the parameters.

If Phi were surjective, preimages of the standard positions
would have p-order zero. They would give polynomial inverse
coordinates for f_t. But f_t is conjugate to f by nonzero
scalar dilations, so this would make f invertible.
Conversely, a polynomial inverse of f gives the inverse
cotangent Weyl map. Thus
\[
                 \Phi\text{ surjective}
                       \Longleftrightarrow f\text{ invertible}.
                                                               \tag{8}
\]
This does not assume that a formal inverse is polynomial.

## Converse for assignments preserving the multiplication algebra

Fix a characteristic-zero field \(K\), set \(h=1\), and suppose an
endomorphism \(\Phi\) of the finite Weyl algebra sends each position
\(x_i\) to a polynomial \(f_i(x)\). The Weyl relations give
\([\Phi(p_i),R(f)]=(\partial_i R)(f)\) for every polynomial \(R\).
A nonzero relation of minimum total degree among the \(f_i\) would
therefore give one of smaller degree by differentiation. Thus the
\(f_i\) are algebraically independent. Their differentials form a basis
over \(K(x)\) by separability in characteristic zero.

For the top differential-order symbol \(S_i(x,\xi)\) of \(\Phi(p_i)\),
the order-\((\deg_\xi S_i-1)\) symbol of its commutator with \(f_j\)
is \(\nabla_\xi S_i\cdot\nabla_x f_j\). If the order is at least
two, these symbols all vanish by the Weyl relations. The basis of
differentials then gives \(\nabla_\xi S_i=0\), which contradicts
its positive degree. The relation with \(f_i\) excludes order zero.
Writing
\[
 \Phi(p_i)=\delta_i+b_i,\qquad
 \delta_i=\sum_j A_{ji}p_j,\qquad A=(Df)^{-1}
\]
The first-order coefficient matrix satisfies the polynomial
identity \(Df\,(a_{ji})=I\), so its determinant also shows \(f\) is
Keller.

The derivations \(\delta_i\) commute and are dual to \(df_i\).
The scalar commutator equations say
\(\delta_i b_j=\delta_j b_i\), so \(\sum_i b_i\,df_i\) is a closed
polynomial one-form. It is \(dU\) for a polynomial \(U\) unique up to
constants. In ordinary \(x\)-coordinates, a coefficient vector
homogeneous of degree \(m\) has primitive obtained by contracting
with the Euler vector field and dividing by \(m+1\).
Consequently \(\Phi=\Gamma_U\circ\Phi_f\), where the shear
\(\Gamma_U\) fixes \(x\) and sends \(p_j\) to \(p_j+\partial_j U\),
and \(\Phi_f\) is the cotangent assignment. Conversely these
assignments all satisfy the Weyl relations. Bavula's Proposition 5
classifies the extensions once \(f\) is Keller. The symbol argument
above also recovers this hypothesis from the Weyl assignment itself.
Equation (8) now gives surjectivity precisely when \(f\) is invertible.

In rank two the cotangent assignment has Bernstein degree exactly
\(\deg f\), since its respective position degrees include \(\deg f\)
and its first-order terms have degree at most \(\deg f\). Thus for every
\(d\ge1\), the assertion that all complex plane Keller maps of degree
at most \(d\) invert is equivalent to surjectivity of every rank-two
Weyl endomorphism of Bernstein degree at most \(d\) sending positions
to the multiplication algebra. The same equivalence holds without
a degree cutoff.

For a quantitative inverse bound, suppose \(\Phi\) has Bernstein
degree at most \(D\) and its positional inverse \(g=f^{-1}\) has
degree \(e\). The inverse is \(\Phi_g\circ\Gamma_{-U}\); its
order-zero correction vector is
\[
 -(\nabla U)\circ g=-\bigl((Dg)^{-1}\bigr)^T(b\circ g).
\]
The matrix equality follows by differentiating \(f\circ g=x\).
It has polynomial entries of degree at most \((n-1)(e-1)\) because
\(g\) is Keller, and each \(b_i\circ g\) has degree at most \(De\).
Together with the inverse cotangent term, this gives Bernstein degree
at most
\[
 \max\{e,(n-1)(e-1)+1,De+(n-1)(e-1)\}.
\]
For rank two, \(e\le\deg f\le D\) by the classical polynomial
inverse degree bound, hence the upper bound is \(D^2+D-1\).
[Polynomial inverse-degree theorem as stated by Tsuchimoto, Proposition 4.1](https://ir.library.osaka-u.ac.jp/repo/ouka/all/7472/1700ojm.pdf).

## Arbitrarily high identity contact with fixed forward degree

Choose any integer N>=1 and specialize t=h^N. Then the
finite-PBW map
\[
                       \Phi_h=\Phi_{h,h^N}
\]
satisfies all relations (5), has the same spatial bound (6),
and
\[
                       \Phi_h=\operatorname{id}+O(h^N).
                                                               \tag{9}
\]
Its h-adic completion is an automorphism: the inverse of
the continuous linear map id+O(h^N) is its convergent
Neumann series. As Phi is an algebra map, that inverse
is an algebra map too. Every coefficient of every inverse
image is polynomial, since a fixed h-order uses finitely
many polynomial operations.

Nevertheless, for noninvertible f the generic finite-PBW
map is not surjective by (8). Its unique completed inverse
therefore has no common finite spatial bound.

This failure is visible explicitly on the positions. Let
\[
                    g(x)=x+\sum_{k\ge2}g_k(x)
\]
be the formal inverse at the origin of f, with g_k homogeneous
of degree k. Then
\[
                \Phi_h^{-1}(x)=h^{-N}g(h^Nx)
                      =x+\sum_{k\ge2}h^{N(k-1)}g_k(x). \tag{10}
\]
If only finitely many g_k were nonzero, the two formal inverse
identities would be polynomial identities, making f a
polynomial automorphism. Thus a noninvertible input gives
nonzero coefficients of arbitrarily high spatial degree in
(10). Increasing N postpones the failure to higher parameter
orders and leaves the forward bound unchanged.

## A retained rank-three control of forward degree twelve

Use [F1's certified input](../../foundations/verified/FOUNDATIONAL_GEOMETRY.md).
With u=1+xy, write its outputs as
\[
 a=u^3z+y^2u(4+3xy),\quad
 b=y+3xu^2z+3xy^2(4+3xy),\quad
 c=2x-3x^2y-x^3z.
\]
Take
\[
                              f=(c/2,b,a).          \tag{11}
\]
It satisfies the normalization in (1), has component degrees
(4,6,7), and has the three distinct rational collision points
\[
 (0,0,-1/4),\quad (1,-3/2,13/2),\quad
 (-1,3/2,13/2)
\]
over the common value (0,0,-1/4). The retained F1 theorem
already establishes the determinant and collision; its
provenance is not reclassified by the present result.

For this input the inverse Jacobian column degrees are
(11,9,8), so the six Weyl images in (4) have degrees
\[
                             (4,6,7,12,10,9).        \tag{12}
\]
The common bound12 holds for every N. The scaled collision
points are the displayed points divided by t, with common
image (0,0,-1/(4t)). They belong to K((h)) when t=h^N,
so generic nonsurjectivity is exact.

Identity padding gives the same examples in every n>=3.
This does not descend from rank three to rank two. Higher-rank
Weyl counterexamples are already known; the new use here is
as a fully compatible control for proposed bounded comparison
and inverse arguments.

## Even completed graph isomorphism does not give a bounded gauge

For N>=2 the completed automorphism in (9) is inner by an
h-adic unit:
\[
                       \Phi_h(a)=U a U^{-1},\qquad
                       U\in1+h^{N-1}\widehat A.     \tag{13}
\]
Here is the precise formal input. At a first difference of
order k>=N, the linearized full Weyl relations form a closed
polynomial one-form. Polynomial H1 vanishes, as in
[QFL1](PROOF.md), so that difference is a normalized
Hamiltonian derivation. Remove it by conjugation with
exp(h^(k-1)*G_k). Since k>=2 this exponential is an h-adic
unit. The successive product converges, proving (13).
No finite spatial bound on U, or on the inverse map, follows.

Let M_Phi be the algebraic graph bimodule: its underlying left
module is the regular finite PBW algebra A and its right
action is m*a=m*Phi_h(a). It is left free of rank one and
carries the full right Weyl action. It is a coherent bimodule,
being cyclic over the Noetherian algebra A tensor A^op.

Multiplication m to mU is an isomorphism from its completion
to the completed identity graph, by (13). Generically M_Phi
is not an invertible bimodule: invertibility would make the
right-action map A^op to End_left(M_Phi)=A^op an isomorphism,
forcing Phi to be surjective. The identity graph is invertible.
Thus these algebraizations cannot be compared by bounded
source and target Weyl automorphisms. Surjectivity is preserved
by every such two-sided change.

This is stronger than a failure of finiteness for a selected
lattice: both left modules here are already free, and the full
right actions exist. It does not contradict any uniqueness
theorem requiring prescribed arithmetic support. The common
classical h=0 map does not prescribe that support.

## Two consequences for degree and jet arguments

At t=0 with h=1, (4) is the identity of a fixed Weyl algebra;
for generic t it is not surjective. Hence the automorphism
locus in a bounded-degree endomorphism scheme need not be
open at the identity. The same conclusion holds at the
classical identity point in QFL2's h-relation scheme, after
the degree-nonincreasing change from normal to Weyl symbols.
Its automorphism locus being closed does not imply openness.

For every m, the reduction of (9) modulo h^m is an automorphism,
by a finite Neumann series. Their inverse degrees are unbounded
as m grows, by (10). Thus an inverse-degree theorem over
fields or reduced rings cannot be applied to all these
nonreduced jets. Any fixed inverse bound for all the jets
would give a common bound for the unique completed inverse,
contradicting (8).

## The rank-two comparison gap contains JC2

For fixed n and N>=2, consider the statement P_(n,N):

> Every exact finite-PBW Weyl assignment over K[[h]] equal to
> the identity modulo h^N has a finite-PBW inverse.

For characteristic-zero fields, with the usual passage from
finitely generated coefficient fields to C,
\[
                    {\rm DC}_n\ \Longrightarrow\
                    P_{(n,N)}\ \Longrightarrow\ {\rm JC}_n.
                                                               \tag{14}
\]
For the first implication, DC makes the generic Weyl map an
automorphism. QFL2's valuation argument, applied to its
injective identity reduction, puts its finite inverse
coefficients back in K[[h]]. For the second implication,
apply P to (1)--(9) for an arbitrary normalized Keller map
and use (8). Over C one can state all three assertions over
that field; a generic counterexample over a coefficient
field extension yields one over C by the standard
finitely-generated-field argument of QFL2.

Restricted to the cotangent rescalings (1), the inverse or
bounded-comparison assertion is equivalent to JC_n. If f is
invertible, its normalized inverse g gives the polynomial
inverse family t^{-1}g(tx), regular at zero. Conversely (8)
recovers invertibility of f.

Accordingly, a universal bounded gauge comparison for full
identity lifts is not a routine consequence of formal
Hamiltonian equivalence. It is false in every n>=3 and,
even on the cotangent subclass in n=2, would settle JC2.
No equivalence P_(n,N) to DC_n is claimed, and no assertion
P_(2,N), JC2 or DC2 is proved or refuted here.

## Current consequence

Fixed forward degree, arbitrary finite identity contact and a completed graph
isomorphism do not provide a bounded inverse or bounded comparison gauge.
The complete rank-two realization problem remains open. The all-rank
construction, differential-order argument and inverse nontermination above
are written proofs using the explicit retained Keller input; no formal
verification or independent external review is claimed.

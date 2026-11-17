# Global reconstruction obstructs embeddings with one additional coordinate

14 November 2026. The universal finite-fibre templates leave a genuine
embedding obligation. A classical single-generator criterion supplies
an all-degree restriction on that obligation: a nonautomorphic Keller
map cannot become a closed embedding by adjoining just one polynomial.

This closes the embedded subclass of the earlier universal cotangent
encoding while preserving its exact A/B identities, for arbitrary
polynomial choices of the remaining dual coordinate. It also excludes
all fibre-affine five-coordinate encodings whose dual linear matrix
depends only on the Keller outputs, without a global frame assumption.
Neither result assumes a bound on graph or matrix degree.

## 1. One global additional generator forces a Keller inverse (HC4MYHI0)

Let k have characteristic zero, B=k[z1,...,zm], and
A=k[K1,...,Km] for a polynomial Keller map K, det DK in k nonzero.
If there is h in B with

\[
                              B=A[h],                       \tag{1}
\]

then **A=B**, so K is a polynomial automorphism.

This is the classical Formanek single-generator criterion; see the
[original publication](https://link.springer.com/article/10.1007/BF01210711).
The [source record](sources/reconstruction-import.json) distinguishes
the inspected publisher metadata from the subscription full text,
which is not used as an unproved input. Here is a complete argument.

The invertible Jacobian makes K1,...,Km algebraically independent,
so A is a polynomial UFD. Also Frac(B) is algebraic and finite over
Frac(A). Choose the primitive irreducible polynomial R(Y,T) in A[T]
of h over Frac(A). Gauss's lemma gives

\[
                           B=A[T]/(R).                      \tag{2}
\]

The Jacobian matrix expresses every dz_i as a B-linear combination
of the dK_j. Hence the relative differentials Omega_(B/A) vanish.
From (2), their presentation is

\[
                    \Omega_{B/A}=B/(R_T(K,h))\,dh.          \tag{3}
\]

Thus R_T(K,h) is a unit in B. Polynomial-ring units are k nonzero,
so it is a scalar c nonzero. The polynomial R_T-c is in the kernel
(R) of (2), but has T degree smaller than R. It is therefore zero.
Characteristic zero now gives R=cT+d(Y). Consequently h is in A,
and (1) implies B=A.

The basic free-module and quotient presentations of differentials
used here are [Stacks, Lemmas 10.131.7, 10.131.9 and 10.131.14](https://stacks.math.columbia.edu/tag/00RM),
with inspected primary bytes in the source record. No classification
of Keller maps or finite fibres is used.

Equality A=B gives polynomial expressions z_i=H_i(K). These supply
both compositions of an actual inverse: K(H(Y))=Y follows from
algebraic independence after substituting Y=K(z). This is a ring
generation statement, not just a primitive-element field statement.

For example C[t,t^(-1)]=C[t][t^(-1)] is an unramified monogenic
extension with relation tT-1, but its unit group is larger than C nonzero.
In positive characteristic T^p-T-Y has constant nonzero T derivative
without being linear. These show exactly why trivial source units
and characteristic zero are essential to the proof.

## 2. An HC4 counterexample needs two reconstruction polynomials (HC4MYHI1)

For a polynomial f on A4 with unit full Hessian put g=grad f,
B=k[source], A=k[g1,g2,g3,g4]. If every source coordinate is a
polynomial in g and one **common** additional polynomial h, then
HC4MYHI0 makes g a polynomial automorphism. This needs no Tang theorem.

Equivalently, a nonautomorphic g admits no closed embedding

\[
                   z\longmapsto(g_1(z),...,g_4(z),h(z))
                                  \quad\text{into }A^5.    \tag{4}
\]

There is also a direct graph formulation. Let i:A4->A5 be any actual
polynomial closed embedding, including any graph followed by an
arbitrary source automorphism, and f=Psi composed with i. If global
polynomial recovery identities put all five coordinate functions of
i in k[grad f,h] for one common h, then the embedding gives
B=k[i1,...,i5]=k[grad f,h], and the gradient is invertible under unit
Hessian. For a graph parameterization it suffices to recover its four
retained physical coordinates that generate B. All degrees are free.

In particular, if a polynomial automorphism of the ambient A5 makes
the gradient into the restriction of four ambient coordinate functions,
the fifth is such an h. No HC4 counterexample can have that ambient
coordinate-projection realization. A general graph embedding does
**not** automatically have this property: five source generators versus
four gradient generators is not a relative one-generator theorem.

Under hypothetical nonautomorphy, Tang gives

\[
 f=D(z)+w^T K(z),\qquad J(K)\in k^\times.
\]

Set eta=grad D+DK^T w. This is a polynomial source coordinate change
fixing z, with inverse w=DK^(-T)(eta-grad D). Then

\[
             B=k[z1,z2,eta1,eta2],\qquad A=k[K1,K2,eta1,eta2]. \tag{5}
\]

Thus two additional functions z1,z2 generate B over A, while one
cannot. Over C every hypothetical nonautomorphic HC4 gradient has
**exactly two** as its minimum number of additional algebra generators.
The map (grad f,z1,z2) is an actual closed embedding into A6 over the
gradient target, with polynomial source recovery supplied by (5).
This minimum is over the gradient map; A4 itself of course embeds in
A5. The absolute embedding dimension is not the assertion.

## 3. The universal cotangent substitution cannot be upgraded to an embedding with the same data (HC4MYHI2)

Take any plane Keller pair K=(F(P,Q),G(P,Q)), independent of x,y.
Consider arbitrary polynomials p,q,r in x,y,P,Q satisfying the exact
Meng--Yang linear-form identities

\[
                       A=F,\qquad B=xF+yG.                   \tag{6}
\]

If

\[
          i(x,y,P,Q)=(x,y,p(x,y,P,Q),q(x,y,P,Q),r(x,y,P,Q))  \tag{7}
\]

is a polynomial **closed embedding**, then K is an automorphism.
There is no restriction on the degrees, mixed terms or kernel-direction
corrections to p,q,r preserving (6).

Indeed, base change the closed embedding to the single physical fibre
(x,y)=(0,1). In the original ambient formulas on that fibre,

\[
                         A=p,\qquad B=4p+q,                 \tag{8}
\]

with r the third coordinate. By (6), (p,q,r) on the source plane is
(F,G-4F,h(P,Q)) for some polynomial h. Closed embedding means

\[
                k[P,Q]=k[F,G-4F,h]=k[F,G,h].                 \tag{9}
\]

HC4MYHI0 now gives the plane inverse. The corresponding potential
Psi composed with i is

\[
                       LF^2+MF+N(xF+yG).                   \tag{10}
\]

As in retained [HC4MYGF3](FIRST_SLICE_INVERTIBILITY.md#6-sharp-trace-boundary-the-larger-substitution-class-contains-jc2),
its full Hessian is N^4 J(F,G)^2 for L*N nonzero, and its gradient
is invertible exactly when K is. Therefore the **closed-embedded**
subclass of (6) has an invertible gradient in every degree. No
nonautomorphic Keller pair can be encoded by an actual graph or
closed embedding preserving those same A/B data.

The earlier arbitrary-substitution construction remains valid and still
contains JC2; this theorem identifies its embedding obstruction rather
than silently asserting its parameter inverse. A new graph realization
must change that exact encoding, not merely choose a more complicated
remaining dual function. The result does not identify the residual
coefficient pair in every graph with the physical A/B pair in (6).

More generally, any closed embedded source plane in an ambient A3
whose two independent ambient linear coordinates pull back to a plane
Keller pair has an automorphic pair: complete the two linear rows to
an actual linear frame and apply HC4MYHI0 to the third pullback.

## 4. A target-dependent affine dual frame cannot repair a nonautomorphic base pair (HC4MYHI3)

Work over C. Let K:A2_z->A2_t be Keller, let M(t) be any polynomial
3-by-2 matrix, and v(z) any three-entry polynomial column. If

\[
                   j(z,w)=(K(z),\ M(K(z))w+v(z))            \tag{11}
\]

is a closed embedding A4->A5, then K is an automorphism. Unlike a
chosen global-frame argument, this allows M to drop rank away from
the image of K. A hypothetical rank drop at omitted target points
does not evade the obstruction.

Here is a direct hypersurface proof. The closed image has prime
principal ideal (E(t,d)) because it has codimension one in the UFD
C[t1,t2,d1,d2,d3]. It is smooth, being isomorphic to A4. Let
m(t) be the ordered cross product of the two columns of M.
Immersiveness in the two w directions gives rank M(K(z))=2
everywhere. Hence the three entries of m(K(z)) generate the unit
ideal in C[z]. Differentiating E composed with j in w gives

\[
                    M(K)^T E_d(j)=0.
\]

The unimodular cross product therefore gives one polynomial lambda
in C[z,w] with E_d(j)=lambda*m(K). It can be obtained explicitly by
dotting E_d(j) with a Bezout column for m(K).

If lambda vanished at any source point, E_d(j) would vanish there.
Differentiating E composed with j in z and using invertibility of DK
would then also give E_t(j)=0. This contradicts smoothness of the
image hypersurface. Thus lambda is a polynomial unit, a constant
c nonzero.

For each dual index a, E_(d_a)-c*m_a(t) vanishes on the image and is
therefore a multiple of E. The polynomial E has positive total degree
in d: otherwise E(K)=0 contradicts algebraic independence of K.
The alleged multiples have strictly smaller d degree, so each is zero.
Characteristic zero now gives

\[
                            E=c\,m(t)\cdot d+e(t).          \tag{12}
\]

The generic target fibre of this hypersurface is A2 over C(t).
The generic source fibre of (11) is Spec C(z) times A2 over C(t).
Their algebraic constant fields coincide, so C(z)=C(K): the plane
Keller map is birational. The retained birational Keller case in
[HC4MYOR2](../hc4-relative-row-orbit/PROOF.md#6-properness-and-the-exact-missing-length)
makes it an automorphism. Conversely an automorphic K admits such
an embedding with the constant matrix having rows (1,0),(0,1),(0,0),
and v=0.

This closes an entire output-base encoding mechanism, before imposing
any Hessian condition. The descent of M through K is essential: for
a matrix depending on the source z rather than its outputs, m(z)
cannot be replaced by a target polynomial in the degree argument.
Nonlinear dependence on w is also outside (11). The retained singleton
and unit-ideal graph constructions use different physical-base data
and are not excluded by this theorem.

## Global hypotheses, checks and remaining scope

The equality in (1) must hold globally as polynomial algebras. Equality
of fraction fields, a rational recovery on an open chart, or generation
only modulo the zero-fibre ideal does not suffice. In fact every finite
reduced complex Keller fibre is generated by one separating function
as a finite algebra, so that fibrewise condition is automatic. It does
not remove a denominator or prove any global reconstruction identity.

The first three statements are over characteristic-zero fields with
the stipulated polynomial embedding/recovery data; the geometric
proof of HC4MYHI3 is over C. A polynomial inverse obtained after
coefficient-field extension descends by uniqueness. No use is made of
Hessian invariance under nonlinear source changes.

[verify_reconstruction.py](verify_reconstruction.py) and its
[receipt](receipt-reconstruction.json) check the single-relation
degree mechanism, the indispensable unit/characteristic hypotheses,
cotangent reconstruction, the exact physical fibre (8), linear-frame
and polynomial-kernel controls for (11), and retained source pins.
The arbitrary-degree commutative algebra and hypersurface arguments
are written proofs, separate from finite controls. The classical
single-generator theorem is attributed, not asserted new.

Arbitrary graph embeddings that do not meet these global reconstruction
or encoding conditions, the universal SL2 row equation and both distinct
plane endpoints remain open. No full HC4/JC2 solution, formal
verification, independent replay, external review or novelty priority
is claimed.

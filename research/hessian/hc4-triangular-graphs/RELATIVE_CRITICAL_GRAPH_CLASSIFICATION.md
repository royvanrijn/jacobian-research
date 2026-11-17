# Relative critical schemes classify base-only quadratic graph images

14 November 2026. This note treats actual source parameters and allows
arbitrary polynomial automorphisms of two free variables over a retained
two-variable base. Cotangent alignment is not assumed.

## 1. Two relative critical-scheme obstructions (HC4MYHR0)

Let f in C[s,t,u,v] have nonzero constant Hessian determinant, let
E=span(partial_u,partial_v), and set Z_E=V(f_u,f_v) as a scheme. Write
pi:Z_E->A2_(s,t) for the retained projection and H_E for the 2 by 2
vertical Hessian block. Either of the following conditions forces grad f
to be a polynomial automorphism:

1. pi is an isomorphism of schemes, so Z_E is a global polynomial section.
2. For a nonempty smooth closed curve C in the retained base, Z_E is
   isomorphic **over that base** to C times A1, and H_E has rank one
   everywhere on Z_E.

A generic section, pointwise counts without a global scheme statement,
or an abstract product with no specified base projection is not the
stated hypothesis. The full unit Hessian makes Z_E smooth of dimension
two whenever it is nonempty: the two selected rows of the full Hessian
are everywhere independent.

Assume toward a contradiction that grad f is nonautomorphic. Tang gives
its constant jointly affine plane W. The rank of W's projection to the
retained base is zero, one or two. Base-preserving constant affine
changes reduce to the following exhaustive forms.

### Projection rank zero

Here W=E and f=D(s,t)+uF(s,t)+vG(s,t). Its vertical Hessian is zero
and its relative critical scheme is V(F,G) times A2. It is neither a
section nor a rank-one product as in condition 2.

### Projection rank one

Write

\[
 f=D(t,U)+sA(t,U)+V C(t,U),\qquad J_{t,U}(A,C)\in C^*.       \tag{1}
\]

The relative critical equations are

\[
 C(t,U)=0,\qquad D_U+sA_U+V C_U=0,\qquad
 \det H_E=-C_U^2\quad\hbox{on }Z_E.                         \tag{2}
\]

For condition 1, a zero of C with C_U=0 would have A_U and C_t
nonzero by the Keller identity. Choosing s to kill D_U+sA_U would
then leave a whole V line, contradicting the section. Thus C_U is
nonzero on C=0. Every such zero gives one relative critical point
over each s. The section therefore makes C=0 exactly a graph U=h(t),
with h polynomial: its U coordinate is independent of s because it is
algebraic over C(t). The Keller component C is squarefree, so
C=gamma*(U-h(t)); it is a plane coordinate. A Keller pair with a
coordinate component is invertible, contradicting nonautomorphy.

For condition 2, a zero of C with C_U nonzero would give a point of
Z_E with vertical Hessian rank two, impossible. Hence C_U vanishes
on C=0. Squarefreeness and the Nullstellensatz give C divides C_U,
so the degree in U forces C_U=0 identically. Now C=C(t), and the
Keller determinant -A_U*C_t being constant forces C affine-linear
and A linear in U with constant nonzero coefficient. The pair is
again invertible. The same reasoning applies whenever H_E has rank
zero everywhere on Z_E.

### Projection rank two: the graph and its envelope

Write retained coordinates w=(s,t), quotient coordinates z=(U,V), and

\[
 f=D(z)+w^T K(z),\qquad J_K=DK,\quad\det J_K\in C^*.
\]

Every relative critical point is on the polynomial graph

\[
 w=P(z)=-J_K(z)^{-T}\nabla D(z).
\]

Differentiation gives two exact identities on that graph:

\[
 H_E=-J_K^T DP,\qquad d(D+P^TK)=K^T dP.                   \tag{3}
\]

Under condition 1, P is a polynomial automorphism. Its inverse Z(w)
defines the polynomial envelope c(w)=D(Z(w))+w^TK(Z(w)), with

\[
                \nabla c=K\circ Z,\qquad
                \det\operatorname{Hess}c=\det DK/\det DP\in C^*.
\]

The known two-variable Hessian theorem HC2 makes grad c invertible,
so K is invertible, a contradiction. We explicitly use only the n=2
invertibility assertion of [de Bondt, Theorem 2.1](https://arxiv.org/html/1203.6605v3#S2),
which extends Dillen's earlier result. No four-variable conclusion is
imported from that theorem.

For condition 2, the graph identifies A2_z with C times A1 over the
base. The retained curve formula gives E(C)*tau=tau^2, hence C has
one component, genus zero and one omitted completion point. Thus C is
an actual A1, by the same curve argument retained in HC4MYHK0.
Choose its parameter rho and the free parameter sigma. They form
polynomial coordinates on the actual z plane. Write its closed smooth
embedding in the retained base as iota(rho)=(a(rho),b(rho)). Then
P=iota(rho), and a',b' are coprime because iota is an immersion.

The envelope identity (3) says

\[
 a'K_1+b'K_2=c'(rho),\qquad
 K_\sigma=h(rho,sigma)(b',-a').                            \tag{4}
\]

Here the envelope is independent of sigma, and the syzygy is polynomial
because a',b' are coprime. The Keller determinant in (rho,sigma)
forces h to be a nonzero constant. Integrate in sigma:

\[
 K=h\,\sigma(b',-a')+B(rho),\qquad
 \det DK=h^2\sigma(a''b'-a'b'')+h\det(B',(b',-a')).          \tag{5}
\]

The sigma coefficient vanishes. In characteristic zero the resulting
Wronskian identity makes a',b' proportional to a fixed constant row;
coprimeness makes the common polynomial factor a nonzero constant.
Hence iota is an affine line, K has one constant sigma direction, and
its complementary constant output row is affine-linear in rho by its
unit Jacobian. This is an explicit triangular polynomial automorphism,
the final contradiction. This proves both conditions in every degree.

## 2. The full relative graph family has exactly the old endpoints (HC4MYHR1)

Let Q be any nondegenerate quadratic polynomial on A5 and suppose a
closed embedding, in a constant affine ambient frame, has the form

\[
 j(s,t,u,v)=(s,t,h(s,t),H(s,t,u,v),G(s,t,u,v)).              \tag{6}
\]

No degree bound is placed on h,H,G. Closed embedding implies
C[s,t,u,v]=C[s,t,H,G], since h is already in the base. Thus (H,G)
is a polynomial automorphism over C[s,t], with constant nonzero
relative Jacobian delta. Its polynomial inverse is global.

Put f=Q composed with j and assume its **actual** source Hessian is a
unit. If its gradient is nonautomorphic, then f has either no full
critical point or exactly one reduced full critical point. Consequently
its source orbit is respectively the coordinate/unit-ideal endpoint or
the nondegenerate-quadratic/singleton-zero endpoint of HC4MYHG2.

In fact both conditional endpoints already occur inside (6). Thus
existence of a nonautomorphic unit-Hessian member of this whole family
is exactly the disjunction of those two old plane endpoint problems,
not unrestricted JC2. Neither endpoint is settled.

Every unit-Hessian member, including those with automorphic gradient,
therefore has at most one full critical point. In particular a native
base-only graph with two or more full critical points cannot acquire
unit Hessian by any polynomial automorphism of the two free variables
over the retained base. This closes that entire relative orbit in all
degrees; it does not cover arbitrary mixing of all four source variables.

Here is the proof, including all positions of W. Consider the constant
quadratic rank k=0,1,2 on the two ambient directions H,G. This k is not
the rank of the three-dimensional kernel in HC4MYHO1. Completing
constant squares and making base-dependent translations of H,G gives
three cases, all through actual relative polynomial automorphisms.
Relative gradient ideals transform by their invertible relative Jacobian;
on the relative critical scheme the vertical Hessians transform by
congruence. The extra second-derivative terms vanish there because the
relative gradient is zero. No full source Hessian invariance is assumed.

### k=2

The relative potential is a nondegenerate binary quadratic plus D(s,t).
Its relative critical scheme is a polynomial section: solve the two
linear equations in H,G and then use the relative inverse. HC4MYHR0
condition 1 forces gradient invertibility. There is no endpoint here.

### k=1

The relative potential is

\[
                  cX^2+A(s,t)X+B(s,t)Y+D(s,t),\quad c\ne0,
\]

with (X,Y) a relative polynomial coordinate pair. If B is a nonzero
constant, the potential itself is a relative coordinate and has no full
critical point. Unit Hessian then gives exactly the unit-ideal endpoint.
If B is identically zero the relative critical scheme has dimension
three, impossible under the full unit Hessian. Otherwise, if B has a
zero, the relative critical scheme is

\[
                \{B=0\}\times A^1_Y,\qquad X=-A/(2c).
\]

It is smooth under unit Hessian, so the base curve is smooth and reduced.
The vertical Hessian has rank one there; HC4MYHR0 condition 2 forces
gradient invertibility. This completes the rank-one case without any
assumed alignment of E with W.

### k=0

The relative potential is D(s,t)+A(s,t)X+B(s,t)Y. Its relative critical
scheme is V(A,B) times A2, with vertical Hessian zero there. Under
nonautomorphy consider once more the rank of W's projection to the base.
Rank one is excluded by the rank-one calculation following (2).

In rank two, (3) gives DP=0, so P is constant. The graph relative
critical scheme therefore projects to one reduced base point. Thus
V(A,B) is exactly that reduced point. On its free fibre, the full
Hessian has zero vertical block and cross block D(A,B) times the
relative coordinate Jacobian. Its unit determinant forces D(A,B)
invertible at that base point. The other two gradient equations are
then an invertible affine-linear system in (X,Y), so f has exactly
one full critical point.

It remains to treat rank zero, W=E. In the actual source variables
f is affine in u,v, with Keller coefficient pair K=(F,G). Its relative
gradient and the relative chain rule give

\[
 (F,G)=(A,B)D_{u,v}(X,Y),\qquad (F,G)=(A,B)
                  \quad\hbox{as ideals in }C[s,t].          \tag{7}
\]

The ideal equality follows first after extension to C[s,t,u,v] and
then by contraction. At each zero of (A,B), differentiate the first
identity in the base variables. Because A=B=0, it gives

\[
                     J(F,G)=\delta J(A,B).                 \tag{8}
\]

The left side and delta are fixed nonzero constants. Nondegeneracy of
the ambient Q says that the two coefficients of the isotropic H,G
plane are independent affine-linear functions of the other three
ambient coordinates (s,t,h). If h has zero coefficient in both,
(A,B) is an invertible affine base pair and has one zero. Otherwise
a constant output row annihilates the h coefficient and gives a
nonconstant affine-linear base component. Use that as the first
base coordinate: after constant changes (A,B)=(s,P(s,t)).

If P(0,t) is a nonzero constant, the ideal is the unit ideal. It cannot
vanish identically, as (8) would fail on the whole line. Otherwise its
roots are simple, and (8) says that P_t(0,t) has the **same nonzero
value** at every root. A polynomial with distinct roots and that
property has degree one: P_0'-gamma vanishes at every root but has
smaller degree than P_0. Hence (A,B) is a maximal ideal. In this last
rank position the coefficient pair therefore also has no zero or one
zero. The cotangent critical equations give the same full critical count.

This exhausts all three k and all positions of W. HC4MYHG2 identifies
the two nonautomorphic critical-count strata with the asserted source
orbits. No new properness theorem is obtained.

### The two converse conditional realizations

For a unit-ideal Keller pair F,G, choose aF+bG=1 and set
V=D+uF+vG, R=-bu+av. The retained relative inverse gives the closed
embedding (s,t,1,V-st-R^2,R) into Q=st+XY+Z^2. It has the form (6),
and (V-st-R^2,R) is an actual relative coordinate pair.

For a singleton-zero Keller pair normalized to ideal (s,t), the retained
oriented row theorem writes (F,G)=(s,t)C with C an invertible polynomial
matrix. If (X,Y)^T=C(u,v)^T, then j=(s,t,0,X,Y) is closed and the
nondegenerate ambient quadratic Q=sX+tY+Z^2 pulls back to uF+vG.
These are exactly the old conditional endpoint constructions, not examples
of either nonautomorphic pair.

### A closed example with exactly two affine retained projections

Take H=u+v^2, G=v+H^2, h=0, and
Q=sH+tG-(H^2+G^2)/2+Z^2. This ambient Q is nondegenerate. The two
free pullbacks have degrees two and four, so the only independent affine
retained pullbacks are s,t. Their relative inverse is v=G-H^2,
u=H-v^2. The source potential is cotangent in (s,t) with the invertible
coefficient pair (H,G), and has Hessian determinant one. Its vertical
critical section is v=t-s^2, u=s-(t-s^2)^2. This verifies that the
rank-two criterion applies below the three-affine projection threshold.

## Scope and inputs

This class allows any base polynomial h and any relative automorphism
of the two free variables, with no alignment assumption. It does not
allow arbitrary source mixing that makes the retained base nonlinear
in the actual source coordinates. In particular a native base-only
graph with many critical points is not declared impossible to unitize
under all four-variable source automorphisms. General two-affine
embeddings with three genuinely free-dependent unretained coordinates
remain outside (6).

The [checker](verify_relative_critical_graph.py) and
[receipt](receipt-relative-critical-graph.json) separate exact identities
from these all-degree proofs. Tang, the already retained curve/Hodge
formula and the old endpoint equivalences are reused. The sole explicitly
retained low-dimensional invertibility input here is
[HC2 in de Bondt's primary theorem](sources/relative-critical-graph-import.json).
No general cancellation or embedding rectification theorem is assumed.
Both plane endpoints, the unrestricted row equation and HC4/JC2 remain
open. No formal verification, independent replay, external review or
novelty priority is asserted.

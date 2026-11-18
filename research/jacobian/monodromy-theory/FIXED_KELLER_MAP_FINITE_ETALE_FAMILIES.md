# One fixed Keller map for a rational finite etale family

Every generically separable polynomial family over a rational parameter can
be realized, on a nonempty open parameter set, as the complete fibres of
one polynomial Jacobian-one map of affine four-space. The extra coordinate
records the parameter. A weighted normalization clears the rational
coefficient denominators while keeping the cubic gauge constant.

For the alternating Hasse family this gives one fixed polynomial map
A4_Q -> A4_Q with infinitely many complete fibres having local points
everywhere and no rational point. The generic group of this fixed map is
S301. A fixed three-dimensional map for that product family remains open.
The [degree-ten pencil](FIXED_KELLER_MAP_ALTERNATING_FAMILY.md) already has
the stronger three-dimensional realization.

## Theorem 1. A fixed polynomial compiler in dimension four

Let k have characteristic zero and let P(w,T) in k(w)[T] be generically
separable of constant degree N>=3. There are a polynomial map

\[
 \mathcal F:\mathbb A^4_k\longrightarrow\mathbb A^4_k,
 \qquad \det J\mathcal F=1,
\]

a nonempty open subset U of A1_k, and a rational target curve Gamma on U
such that, for every field extension K/k and w0 in U(K),

\[
 \mathcal F^{-1}(\Gamma(w_0))\simeq
                \operatorname{Spec}K[T]/P(w_0,T).
\]

These are the entire scheme-theoretic fibres, not selected components.
The fixed map has generic degree N and generic arithmetic and geometric
groups S_N.

Proof. Choose tau in k such that the linear and cubic Taylor coefficients
of P(w,tau+X) are nonzero in k(w). Such tau exists: the polynomials
P_T and P_TTT are nonzero over k(w), so only finitely many constants
can be their roots. Write

\[
 P(w,\tau+X)=\sum_{j=0}^N c_j(w)X^j,
 \qquad r(w)=\frac{c_1(w)}{c_3(w)}.
\]

Choose a nonzero polynomial L(w) such that, for every j>=4,

\[
 h_j(w)=\frac{c_j(w)}{c_1(w)}r(w)^j L(w)^{j+1}\in k[w].
\]

For example, take L to be a common multiple of the denominators of
(cj/c1)r^j in reduced form. For N=3 the condition is empty and L=1
suffices. This is a denominator-clearing operation in k[w], without
extracting square roots or changing the parameter field.

In source coordinates (w,x,y,z) put

\[
 t=1+xy,\qquad q=t^2z+y^2(1+3t),\qquad V=tq,
\]
\[
 B=y+3xq+\sum_{j=4}^N jh_j(w)t^2x^{j-2}q^j,
\]
\[
 C=x(5-3t)-x^3z-\sum_{j=4}^N(j-2)h_j(w)(xq)^j,
 \qquad\mathcal F=(w,V,-B/2,C).
\]

All four coordinates are fixed polynomials over k. Their degrees are at most

\[
 \max\left(7,\max_{4\le j\le N}
                         (\deg h_j+6j+2)\right),
\]

where zero profiles are omitted and the inner maximum is omitted for N=3.
Indeed t,q have degrees at most two and five, and the j-th high terms
in B,C have degrees at most deg(hj)+6j+2 and deg(hj)+6j.
At fixed w these are
the retained quadratic-gauge compiler formulas with g1=g3=1 and g2=0.
In the chart t!=0, use S=x/t and Q=y+xq. Then

\[
 G_{w,V}(S)=S+VS^3+\sum_{j=4}^N h_j(w)V^jS^j,
 \qquad D=1-SQ+VS^2=t^{-1},
\]
\[
 B=Q+\frac{G_{w,V}'(S)-1-VS^2}{S},\qquad
 C=2G_{w,V}(S)-BS^2.
\]

The source-chart determinant is t and the plane determinant of (B,C)
is -2D. After normalizing the second output to -B/2, the product is
one. Derivatives in w do not affect it because the first output is w.
The identity holds on a dense chart and hence globally.

For target coordinates (w,v,Y,Z), the inverse polynomial is

\[
 E_{w,v,Y,Z}(S)=G_{w,v}(S)+YS^2-Z/2.
\]

Set

\[
 \boxed{\Gamma(w)=\left(w,
       \frac{c_3}{c_1L^2},\frac{c_2}{c_1L},
                         -\frac{2Lc_0}{c_1}\right).}
\]

Its decisive identity is

\[
 \boxed{E_{\Gamma(w)}(S)
                   =\frac{L(w)}{c_1(w)}P(w,\tau+S/L(w)).}
\]

Indeed h_j times v^j equals (cj/c1)L^(1-j) for each j>=4.
The degrees zero, one, two and three agree by the displayed choices
of Z, the normalized linear coefficient one, Y and v respectively.
This proves the identity for every degree N at once.

Delete the finitely many parameter poles, zeros of L,c1,c3 and the
leading coefficient, and zeros of the discriminant. This defines a
nonempty open U. For w0 in U and every K/k, v is nonzero, so V=v
makes t and q units in the entire fibre algebra. Thus all fibre points
lie in the chart. Moreover

\[
 D=G_{w_0,v}'-BS=E_{\Gamma(w_0)}'
\]

is a unit modulo E, by separability. Put B0=-2Y and
beta=(G'-1-vS²)/S. The inverse in K[S]/E is

\[
 Q=B_0-\beta,\quad t=D^{-1},\quad x=S/D,\quad
 y=Q-vS,\quad q=vD,
\]
\[
 z=D^2\{vD-y^2(1+3/D)\}.
\]

These reconstruct every source equation, as in the
[whole-fibre compiler proof](../../../unpublished/common-arithmetic-fibers/sections/03-realization.tex).
The affine change T=tau+S/L identifies this full algebra with K[T]/P.
Every construction commutes with field extension, so the assertion is
scheme-theoretic, including reducible finite etale fibres.

Over k(w,V), G has degree N and linear coefficient one. The independent
quadratic and constant outputs Y,Z put its inverse in the scope of the
[quadratic-tilt theorem UAQTILT1](../universal-keller-structure/PROOF.md).
It gives arithmetic and geometric S_N and hence generic degree N.
For N=3 the nonzero cubic term is V; for N>=4 the top profile hN is
nonzero. This checks the degree hypothesis in both cases.

## Corollary 2. Arbitrary rational finite etale families

Every finite etale k(w)-algebra A of rank N>=3 admits the preceding
construction after restricting to a nonempty open subset of A1_k.
In particular one fixed map realizes all its separable specializations
on that open subset as entire fibres.

Proof. A finite etale algebra over the infinite field k(w) is monogenic.
One proof chooses a separable closure and an element whose N images under
the distinct embeddings are pairwise distinct. The failures of this
condition form finitely many proper linear subspaces in its underlying
k(w)-vector space. An infinite field is not a union of finitely many
proper linear subspaces, so such an element exists. Its minimal
polynomial has degree N and is separable, giving A=k(w)[T]/P.
Choose an affine open model and invert its presentation denominators
and discriminant. The resulting algebra is the stated finite etale
family. Theorem 1 realizes its complete fibres.

## Corollary 3. Infinitely many Hasse failures in one fixed map

There is one polynomial map

\[
 \mathcal F:\mathbb A^4_{\mathbb Q}\longrightarrow
                          \mathbb A^4_{\mathbb Q},
 \qquad \det J\mathcal F=1,
\]

of generic degree 301 and generic arithmetic and geometric group S301,
and a rational target curve on which infinitely many complete fibres have
points over every completion of Q and no Q-point. Their finite-algebra
splitting fields can be chosen linearly disjoint with group A10.

Proof. Let M/Q(w) be the
[regular degree-ten A10 splitting field](../../arithmetic/RATIONAL_DEGREE_TEN_ALTERNATING_PENCIL.md).
Use the subgroup classes H1,H2,H3,W5 in the
[Hasse theorem](../../arithmetic/ALTERNATING_SUBGROUP_COVERS_AND_HASSE.md)
and set

\[
 A=M^{H_1}\times M^{H_2}\times M^{H_3}\times M^{W_5}.
\]

It is finite etale of rank 301 over Q(w). Apply Corollary 2 and remove
the finite additional presentation exceptions. The infinite integer
specializations in the Hasse theorem, with linearly disjoint A10
closures, survive deletion of that finite set. The whole-fibre
isomorphism preserves their local and global points. The same single
map works for every remaining parameter. S301 describes the ambient
four-dimensional extension; A10 describes the distinguished product
fibre splitting fields.

## Verification and dimensional scope

The [exact checker](verify_fixed_keller_finite_etale_family.py) uses a
rational degree-five family with varying linear and cubic coefficients
and genuine coefficient poles. It checks denominator clearing, the full
two-variable identity, symbolic Jacobian factors and complete quotient
reconstruction at an unbranched parameter. The
[receipt](fixed_keller_finite_etale_family_verification.json) pins the
input formulas, compiler and imported proof sources. These finite
controls supplement the uniform written argument; no external review or
formal verification is claimed.

The dimension-four construction records the family parameter as an
independent output. It does not remove that coordinate or establish a
fixed A3 realization of arbitrary product families. The earlier
three-dimensional degree-ten pencil construction uses its special
coefficient dependence to eliminate this extra coordinate.

# Three affine retained projections suffice in every degree

14 November 2026. The ambient quadratic graph obstruction extends
from four retained affine coordinates to **three**. The two remaining
ambient coordinates may be arbitrary polynomials in all four source
variables. Their degrees and fibre shapes are unrestricted.

**HC4MYHN0.** Let Q be a nondegenerate quadratic polynomial on A5,
j:A4->A5 a polynomial closed embedding, and f=Q composed with j.
Suppose three constant linear combinations of j's five coordinates
have affine-linear pullbacks with independent source linear parts.
If the actual four-variable Hessian of f is a unit, grad f is a
polynomial automorphism.

Through the explicit ambient automorphism Phi of HC4MYHM0, the same
conclusion holds for the original Meng--Yang potential whenever three
normal coordinate combinations pull back affinely and independently
in the actual source parameters. Being merely polynomial coordinates
does not replace that hypothesis.

## 1. Reduction and retained endpoint inputs

All statements are over C. Affine source and ambient changes put

\[
 j(s_1,s_2,s_3,t)=(s_1,s_2,s_3,H(s,t),G(s,t)),\qquad
 C[s,t]=C[s,H,G].                                            \tag{1}
\]

The second equality is exactly the closed-embedding hypothesis.
The only possible two-dimensional kernel V of the ambient retained
projection has quadratic rank 0, 1 or 2. Elementary constant linear
algebra, completing squares and affine shifts of H,G by the retained
coordinates give these exhaustive ambient forms, up to a constant:

\[
 \begin{array}{c|c}
 \operatorname{rank}(Q|_V)&Q\text{ in retained coordinates}\\\hline
 2&q_3(s)+HG,\quad q_3\text{ nondegenerate ternary quadratic},\\
 1&H^2+s_1G+s_2s_3,\\
 0&s_1H+s_2G+s_3^2.
 \end{array}                                                \tag{2}
\]

For rank two use the orthogonal complement of V and split its
nondegenerate binary form as a product. For rank one first split
off its nondegenerate square; its radical pairs with one retained
coordinate, leaving a nondegenerate binary form on the other two.
For rank zero, the totally isotropic V pairs with two retained
coordinates and leaves one nondegenerate square. Cross terms
divisible by the paired retained coordinates are absorbed by affine
shifts of H,G. These changes preserve (1) and the affine source
meaning of all retained coordinates.

Suppose f has a nonautomorphic unit-Hessian gradient. The retained
Tang alternative supplies a constant jointly affine plane W and
constant source coordinates (z1,z2,w1,w2) with

\[
                         f=D(z)+w^T K(z),\qquad J(K)\in C^\times. \tag{3}
\]

Let L be the constant source line ker d(s1,s2,s3), originally the
t direction. Either L is contained in W or L intersects W only at
zero. The proof treats both positions in every row of (2).

We use two already proved endpoints, not a new classification of
relative line embeddings:

- [HC4MYHI0](RECONSTRUCTION_EMBEDDING_OBSTRUCTION.md#1-one-global-additional-generator-forces-a-keller-inverse-hc4myhi0): if C[z1,z2]=C[K1,K2,h] for one common polynomial h, a plane Keller pair K is invertible.
- [HC4MYGR0](CRITICAL_COMPONENTS.md#1-a-plane-keller-fibre-with-an-affine-line-component): an affine-line component of a plane Keller fibre makes the pair invertible. In this proof the line is already straight. Its mate has nowhere-zero polynomial derivative on the line, hence is affine-linear there, so the retained [line-injectivity theorem](https://arxiv.org/abs/alg-geom/9305008) applies directly.

A nonconstant affine-linear component of K also makes K triangular.
In every such case cotangent reconstruction gives a polynomial
inverse to grad f, contradicting (3)'s nonautomorphy assumption.

## 2. Nondegenerate kernel: q3(s)+HG

### The free source line belongs to W

Then f is affine in t, so the product H*G has t degree at most one.
Neither H nor G can be identically zero, since that would leave f
independent of t. By (1) at least one varies with t. After interchanging
the factors, H=h(s) and G=a(s)t+b(s). Generation in (1) forces a(s)
to be a nonzero constant: a relative one-variable coordinate must
be linear with unit leading coefficient.

If h is constant, f_t is constant and the full Hessian is singular.
Otherwise joint affineness along W implies h is independent of the
second W direction. Adapt constant retained coordinates (v,z1,z2)
to that direction, also adjusting t linearly, so

\[
 f=a h(z)t+h(z)b(v,z)+q_3(v,z),\qquad W=\langle\partial_t,\partial_v\rangle.
\]

The second v derivative says h*b_vv+(q3)_vv=0. A nonconstant h
cannot divide a nonzero field constant. Thus (q3)_vv=0 and b_vv=0.
Write b=b0(z)+v*b1(z) and q3=q2(z)+v*ell(z). Nondegeneracy of q3
makes ell nonconstant affine-linear. The Keller coefficient pair is

\[
                            (a h,\ell+h b_1).
\]

Its Jacobian being a unit implies J(h,ell)-c is divisible by h
for a nonzero constant c. But deg J(h,ell)<deg h whenever that
derivative is nonzero. Hence J(h,ell)=c. The pair (h,ell) is
triangular in constant coordinates, so h is a plane coordinate.
A Keller pair with that coordinate as a component is invertible.

### The free source line is transverse to W

Choose constant coordinates so retained s=(w1,w2,z1) and the free
coordinate is z2. The quadratic q3 restricts nontrivially to W:
otherwise a two-plane would be totally isotropic in a nondegenerate
three-dimensional form. Since f is affine in w, H*G has w degree
at most two. Its nonzero quadratic top is the negative of q3's
constant W-quadratic part.

If one factor has w degree zero, that factor divides a nonzero
constant coefficient of this quadratic top, so it is a nonzero
field constant. Then (1) makes the other factor a relative
one-variable coordinate in z2, of the form c*z2+b(s), c in C nonzero.
Consequently f_z2 is constant, again contradicting unit Hessian.

The only remaining degrees are one and one. Write

\[
                   H=a(z)\cdot w+h(z),\qquad
                   G=b(z)\cdot w+g(z).
\]

Their leading product (a(z) dot w)*(b(z) dot w) is a fixed nonzero
quadratic in w, independent of z. Total z degree of a nonzero
polynomial product adds, so both a,b are constant vectors. The
coefficient pair K is an affine-linear function of z1 plus a*g+b*h.
If a,b are dependent, one nonzero output row annihilates both and
gives an affine-linear Keller component. If independent,

\[
                         C[z_1,K_1,K_2]=C[z_1,h,g].
\]

Quotient (1) by w=0 to obtain C[z1,z2]=C[z1,h,g]. Thus one
common additional generator z1 suffices over K, and HC4MYHI0
gives an inverse. This closes the full rank-two row of (2).

## 3. Rank-one kernel: H^2+s1 G+s2 s3

### The free source line belongs to W

On s1=0, f=H0^2+s2*s3 is affine in t. Degree forces H0 independent
of t, so the Keller component f_t vanishes on the whole source
hyperplane s1=0. If ds1 is nonzero on W, that hyperplane projects
onto the whole quotient plane by W, making f_t identically zero,
impossible. Otherwise it projects onto a straight affine line in
the quotient. The nonconstant Keller component f_t vanishes on
that line, so the affine-line endpoint closes the case.

### The free source line is transverse to W

If W is contained in ker ds1, it projects isomorphically onto the
retained (s2,s3) plane. In f|s1=0=H0^2+s2*s3, cancellation of the
W-quadratic part would require a rank-one square to cancel a
nondegenerate binary quadratic. Degree excludes higher W terms,
so this is impossible.

Otherwise choose u spanning W intersect ker ds1 and v in W with
ds1(v)=1. The vector u is not the free line L, so its projection
to (s2,s3) is nonzero. On s1=0, degree in the u direction gives

\[
                              H_0=a w+h(z),                 \tag{4}
\]

where w is the u coordinate and a is a constant, since a squared
cancels the constant quadratic coefficient of s2*s3 in that
direction. If a=0, f_u on the trace is the nonconstant affine-linear
derivative of s2*s3. Joint-affine transport along v extends it
globally as an affine-linear Keller component, closing this subcase.

Let a!=0. Adapt constant coordinates (s,z1,z2,w) so s=s1 and
W=span(partial_s,partial_w). The source scheme

\[
                             \Pi=\{s=0,H_0=0\}
\]

is an actual affine plane with coordinates z, because (4) solves
w=-h(z)/a polynomially. Restrict the closed embedding to Pi.
Its only remaining ambient coordinates are (s2,s3,G), hence

\[
                             C[z]=C[s_2|_\Pi,s_3|_\Pi,G|_\Pi]. \tag{5}
\]

Put q=s2*s3, expressed in the adapted constant source coordinates.
For the Keller pair (A,B)=(f_s,f_w), the first derivatives of the
full potential on Pi give

\[
                            A=G+q_s,\qquad B=q_w.           \tag{6}
\]

Both q_s and q_w are affine-linear functions of the ambient retained
s2,s3. The form q_w is nonzero, because u has nonzero projection
there; A also has coefficient one on G. Thus A,B are two independent
ambient affine-linear coordinates on A3_(s2,s3,G). Complete them
to an affine linear frame by one more coordinate C. Equation (5)
becomes C[z]=C[A,B,C|Pi], and HC4MYHI0 makes (A,B) invertible.
No normal derivative H_s needs to be recovered: its coefficient
2H vanishes exactly on Pi. This closes the rank-one row.

## 4. Totally isotropic kernel: s1 H+s2 G+s3^2

### The free source line belongs to W

The Keller component f_t vanishes on the source plane s1=s2=0.
If W were that plane's direction space, its trace s3^2 would not
be jointly affine. Otherwise the plane projects onto a straight
affine line in the quotient by W. The affine-line endpoint applies.

### The free source line is transverse to W

Consider the rank of (ds1,ds2) on W. Rank zero makes W the direction
space of s1=s2=0, contradicting the trace s3^2. In rank one its
kernel is a nonzero u in W tangent to that plane. Since W misses L,
u moves s3, so D_u^2(s3^2)!=0, again a contradiction.

In rank two, choose the physical pair w=(s1,s2) as the W coordinates
and let z be affine coordinates on their zero plane. The restriction
of the closed embedding gives

\[
                             C[z]=C[s_3^0,H_0,G_0].          \tag{7}
\]

Write K=(f_w1,f_w2) for the Keller pair from (3). Differentiating
at w=0 gives

\[
                  K_i=H_i^0+2s_3^0\,\partial_{w_i}s_3,
                 \qquad (H_1^0,H_2^0)=(H_0,G_0).            \tag{8}
\]

The derivatives of the retained affine s3 are constants. Thus (7)
is C[z]=C[K1,K2,s3^0], another application of HC4MYHI0. This
finishes all three kernel ranks and both positions of W, proving
HC4MYHN0 without degree bounds or a relative line-rectification theorem.

## 5. The extension is strict; two affine projections retain both endpoints (HC4MYHN1)

For example, in rank-two ambient form take

\[
 q_3=s_1s_2+s_3^2,\quad H=s_3+s_2^3,\quad G=t+s_2^2.
\]

The embedding (s1,s2,s3,H,G) has polynomial source recovery
t=G-s2^2. Its affine-linear coordinate combinations have independent
source rank exactly three: the distinct s2^3 and s2^2 terms cannot
cancel between H and G. Yet the pullback has unit Hessian, because
in quotient coordinates (s2,s3) its cotangent pair is
(s2,s3+s2^3). This is a genuine extension of the four-affine graph
criterion. Arbitrarily high distinct powers work the same way.

The two-affine threshold already contains both old open endpoints;
no exclusion of that whole class is proved.

For a singleton-zero pair, normalize the zero to the origin and
its orientation so K=(s,U)C(s,U) with det C=1, using HC4MYOR1.
For source w=(P,Q), put (a,b)^T=C(P,Q)^T and define normal ambient
coordinates for lambda_star+L*T^2+NyZ+2NxW by

\[
             x=s/(2N),\quad y=U/N,\quad T=0,\quad Z=b,\quad W=a.
\]

The pullback is lambda_star+w dot K. Its source recovery uses the
retained affine x,y and C inverse. It is a genuine closed embedding,
and unit Hessian/gradient nonautomorphy are exactly those of the
plane Keller pair, as in HC4MYGP1.

For a unit-ideal pair (F,G), choose aF+bG=1 and set

\[
 V=D(s,U)+PF+QG,\quad R=-bP+aQ.
\]

The actual relative source inverse is
P=a(V-D)-GR, Q=b(V-D)+FR. Now embed by

\[
 x=1/(2N),\quad y=s,\quad T=U,\quad Z=R,\quad
 W=V-\lambda_*-LU^2-NsR.
\]

Its pullback is V, its two retained normal projections y,T are
affine source coordinates, and its inverse follows from the displayed
source formulas. If (F,G) is Keller, V has unit Hessian and its
gradient misses zero, so is nonautomorphic. This is exactly the
existing unit-ideal endpoint, not a constructed example of one.

Applying Phi inverse gives the corresponding genuine Meng--Yang
embeddings. Consequently a universal extension of HC4MYHN0 from
three affine retained projections to two would resolve both distinct
plane endpoints. The whole two-affine class is not asserted equivalent
to either endpoint or to unrestricted JC2.

## Evidence and remaining scope

The [checker](verify_three_affine.py) and
[receipt](receipt-three-affine.json) verify the linear-algebra forms,
product/degree mechanisms, restricted-plane derivative identities,
source recovery, strict-extension control and two conditional endpoint
embeddings. The exhaustive arbitrary-degree case split is the written
proof above, not a finite rank/degree enumeration.

The retained inputs are Tang's reproduced structural alternative,
the fully reproduced classical Formanek criterion, the affine-line
Keller endpoint and the explicit ambient quadratic normalization.
No new external classification, relative A1 embedding rectification,
nonlinear source Hessian invariance or gradient-collision transfer is
assumed. No formal verification, independent replay, external review
or novelty priority is claimed.

Any hypothetical nonautomorphic unit-Hessian realization now has
at most two independent affine-linear retained normal source
coordinates. Arbitrary source reparameterizations, the five-generator
quadratic problem below that threshold, the universal SL2 row equation,
both plane endpoints and unrestricted HC4/JC2 remain open.

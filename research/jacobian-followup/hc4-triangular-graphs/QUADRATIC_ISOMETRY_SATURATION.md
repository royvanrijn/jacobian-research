# Polynomial quadratic symmetries and the invariant embedding frontier

14 November 2026. The visible affine-coordinate count is not an invariant
of a quadratic ambient representation. Polynomial changes of the ambient
coordinates can remove every affine pullback while preserving the actual
source potential exactly. This gives a uniform extension of the existing
obstructions and a necessary correction to the nonlinear-image search gate.

## 1. Saturation under every polynomial quadratic similitude (HC4MYHP0)

Write the nondegenerate ambient potential in centered coordinates as

\[
 Q=\lambda_*+q,\qquad q=xy+zw+t^2.
\]

Let A be any polynomial automorphism of A5 with

\[
                         q\circ A=a q,\qquad a\in\mathbb C^*.       \tag{1}
\]

For a closed embedding j:A4->A5 set j0=A inverse composed with j.
The source parameters are unchanged, and the identity is literal:

\[
                 Q\circ j=\lambda_*+a(Q\circ j0-\lambda_*).        \tag{2}
\]

Consequently the actual source Hessians differ by the constant scalar a,
their determinants by a^4, and their gradients by a. Their gradient maps
are polynomial automorphisms simultaneously. No nonlinear source Hessian
invariance is invoked. Both j and j0 are genuine closed embeddings because
the ambient change has a polynomial inverse.

Conjugate A by the retained ambient normalization Phi (including a constant
linear normalization of its quadratic form). This gives a polynomial
ambient symmetry of the original Meng--Yang potential about lambda_*.
Thus all statements transfer to genuine original-potential embeddings,
again with the same actual source parameters.

It follows that the entire polynomial-similitude saturation of the
[three-affine projection class](THREE_AFFINE_PROJECTIONS.md) is closed:
if j0 has three affine-linear retained projections with independent source
linear parts, unit Hessian of Q composed with j forces an invertible
gradient. There is no degree bound on A, its inverse, or j0.

The same statement applies to every previously closed actual-source
embedding class for this Q. Existence of A is a sufficient certificate;
we give no algorithm deciding its existence for an arbitrary embedding.

For arbitrary source parameterizations, the entire saturation of the
[affine-hyperplane image class](TWO_AFFINE_ENDPOINT_GEOMETRY.md#3-every-affine-hyperplane-image-has-only-the-two-old-endpoint-possibilities-hc4myho2)
still has exactly its old three alternatives. A nonisotropic hyperplane
gives the full nondegenerate Q4 orbit and its singleton-zero endpoint;
an isotropic hyperplane off center gives the coordinate orbit and its
unit-ideal endpoint; an isotropic hyperplane through center has a critical
line and no unit-Hessian member. Formula (2) only rescales the potential
about its ambient critical value, which does not change these conclusions.
Neither endpoint is proved to exist or excluded.

There is no larger polynomial reparameterization of Q hidden in the
fibre-preserving condition. If an ambient polynomial automorphism satisfies
q composed with A=h(q), h in C[T], then h has degree one. Indeed q composed
with A has exactly one critical point, whereas a root rho of h' would make
the whole positive-dimensional level q=rho critical. Thus h=aT+b with
a nonzero. The unique critical point is A inverse(0), and comparison of
critical values gives b=0 and A(0)=0. This proves that every such polynomial
fibre reparameterization is already (1). It does not classify all
automorphisms of an individual quadric level.

## 2. Explicit symmetries of unrestricted degree and word length

For arbitrary h in C[x,z,w], define

\[
 E_h(x,y,z,w,t)=(x,y-2th-xh^2,z,w,t+xh).                    \tag{3}
\]

The unchanged coordinates keep h unchanged. Direct expansion gives
q composed with E_h=q, and E_h inverse=E_{-h}. Interchanging x,y or
interchanging the hyperbolic pairs (x,y),(z,w) gives the other three
families. These are polynomial automorphisms, not substitutions that
might lose a source fibre. Every finite word in them and constant
orthogonal transformations satisfies (1) with a=1. Inverse words are
obtained by reversing the order and changing each h to -h.

All degrees and word lengths are allowed by the proof of Section 1.
We do not assert that these elementary transformations generate every
polynomial symmetry of q, or that every embedding is equivalent to an
affine one under such transformations.

## 3. Zero affine pullbacks already occur in the closed class (HC4MYHP1)

Here is an explicit infinite family, with a written degree argument.
Start with the standard hyperplane embedding

\[
            j0(a,b,c,d)=(a,b,c,d,0),\qquad Q\circ j0=\lambda_*+ab+cd.
\]

Apply successively the x-fixed transformation (3) with h=z^m, the
y-fixed version with coefficient y^n, the z-fixed version with coefficient
z, and the w-fixed version with coefficient w^r. Here m,n,r are arbitrary
positive integers. In source notation put

\[
\begin{aligned}
 Y&=b-ac^{2m}, & T_0&=ac^m,\\
 T_1&=T_0+Y^{n+1}, & X&=a-2T_0Y^n-Y^{2n+1},\\
 T_2&=T_1+c^2, & W&=d-2T_1c-c^3,\\
 T_3&=T_2+W^{r+1}, & Z&=c-2T_2W^r-W^{2r+1}.
\end{aligned}                                                       \tag{4}
\]

Then j=(X,Y,Z,W,T3) is closed, with inverse given by the four inverse
ambient transformations followed by the projection to (a,b,c,d).
Its actual potential is still lambda_*+ab+cd, so its Hessian determinant
is one and its gradient inverse is explicit.

Let A0=2m+1 and B0=(n+1)A0. The five total source degrees are

\[
 \deg Y=A0,\quad \deg W=B0+1,\quad \deg X=(2n+1)A0,\quad
 \deg T_3=(r+1)(B0+1),\quad \deg Z=(2r+1)(B0+1).             \tag{5}
\]

At each step the displayed leading term has strictly greater degree
than the other summands. In particular these degrees are strictly
increasing in the order in (5):

\[
 A0 < B0+1 < (2n+1)A0 < (r+1)(B0+1) < (2r+1)(B0+1).
\]

For the middle inequality use r>=1 and
2(B0+1)-(2n+1)A0=A0+2>0. Thus every nonzero constant linear combination
of the five pullbacks has degree greater than one: the term of greatest
degree cannot cancel. There are **zero** nonzero affine-linear retained
normal projections and no affine equation for the image. Nonetheless
the entire family lies in the closed saturated class. This is an exact
all-degree construction, not a search over m,n,r.

The same four transformations, now with m=n=r=1, also hide both isotropic
hyperplane types. For j0=(epsilon,b,c,d,t), epsilon=0 or 1, use h=z in
the first step. At either epsilon the successive formulae are

\[
\begin{aligned}
 T_0&=t+\epsilon c,&Y&=b-2tc-\epsilon c^2,\\
 T_1&=T_0+Y^2,&X&=\epsilon-2T_0Y-Y^3,\\
 T_2&=T_1+c^2,&W&=d-2T_1c-c^3,\\
 T_3&=T_2+W^2,&Z&=c-2T_2W-W^3.
\end{aligned}                                                       \tag{6}
\]

Their degrees in the order (Y,W,X,T3,Z) are (2,5,6,10,15), so they
also have no affine retained projection or affine image equation.
The potential is epsilon*b+cd+t^2. At epsilon=0 the critical line
still forbids unit Hessian after every source parameterization. At
epsilon=1 its full source orbit still has exactly the open unit-ideal
endpoint. The epsilon=1 displayed parameterization itself has degenerate
Hessian; no unit-ideal Keller pair is constructed.

## Corrected frontier and literature boundary

A new mechanism beyond the established branches must escape their
entire polynomial quadratic-similitude saturations. Merely giving a
nonlinear image hypersurface or eliminating all affine pullbacks does
not do that. For a proposed symmetry certificate, check the polynomial
inverse and the full identity (1), then apply it to the embedding in its
actual source parameters. An arbitrary source automorphism on a closed
three-affine representative need not retain that obstruction; the full
source-orbit statement here is only the separately classified hyperplane
branch.

[Blanc--van Santen, Embeddings of Affine Spaces into Quadrics](https://arxiv.org/abs/1702.00779)
provides primary context for the limits of a general rectification claim:
its Theorem 3 constructs inequivalent affine-plane images in SL2. We do
not transfer that theorem to dimension five, import its classification,
or claim a rectification theorem. The proofs above are direct polynomial
identities and consequences of the retained obstructions.

The [checker](verify_quadratic_symmetry.py) and
[receipt](receipt-quadratic-symmetry.json) verify elementary identities,
literal inverses and small exact degree controls. The unbounded statements
are the written proofs. General embeddings outside these saturated classes,
unaligned two-affine cases, both plane endpoints, the unrestricted row
equation and HC4/JC2 remain open. No formal verification, independent replay,
external review, exhaustive ambient symmetry classification or novelty
priority is asserted.

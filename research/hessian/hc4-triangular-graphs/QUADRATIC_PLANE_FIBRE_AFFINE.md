# A quadratic source plane forces alignment; fibre-affine embeddings retain only the old endpoints

14 November 2026. This strengthens the constant jointly affine-plane test
and classifies a whole two-affine embedding family with arbitrary base
coefficient degrees. All three unretained coordinates may depend on the
two free variables. No cotangent alignment or base-only helper is assumed.

## 1. Two elementary plane Keller facts

We reproduce the exact plane facts needed below rather than importing a
degree search or a larger classification.

**One component affine in one variable.** Let (A(t,U),C(t,U)) be Keller
and C=a(t)U+b(t). If a=0, the unit identity -A_U*b'=c forces A_U and
b' to be nonzero constants, so the pair is triangular. If a!=0 and
the current U degree of A is m>=1, take its leading coefficient a_m(t).
The coefficient of U^m in
J(A,C)=c says

\[
                     a a_m'-m a'a_m=0.
\]

Thus a_m/a^m is a constant gamma. Subtract gamma*C^m from A, without
changing its Jacobian with C. This strictly lowers the U degree.
Repeating in the finite degree of A gives A=P(C)+e(t). Now
J(A,C)=a e'=c makes a and e' nonzero constants. The inverse recovers
t affinely from A-P(C), and then U from C. There is no bound on the
t degree or on the original U degree of A.

**Quadratic plane Keller maps.** A plane Keller map of total degree
at most two is invertible. Normalize its constant and linear parts to
write it as z+H(z), H homogeneous quadratic. The determinant identity
gives tr DH=0 and det DH=0 by their different homogeneous degrees.
The two quadratic components of H are proportional: on z1!=0 write
them as z1^2*a(z2/z1), z1^2*b(z2/z1); their Jacobian zero gives
a*b'-a'*b=0 and hence b/a constant. The zero-component cases are
immediate. Thus H=v h for a constant vector v and homogeneous quadratic
h. The trace identity is D_v h=0. Consequently

\[
                (z+v h)^{-1}=z-v h,
\]

because h is unchanged along v. Undoing the affine normalization gives
the inverse. The quadratic case is classical; Wang's
[primary publication metadata](https://www.sciencedirect.com/science/article/pii/0021869380902331)
identifies its quadratic section. The plane argument here is complete
and uses no uninspected general-dimensional result.

## 2. The unique quadratic two-plane under nonautomorphy (HC4MYHT0)

Let f be a complex four-variable polynomial with nonzero constant
Hessian determinant. Let E be a constant two-plane of source directions.
Suppose f has degree at most two on all E translates, uniformly:
in affine coordinates (s,t,u,v) with E=span(partial_u,partial_v),

\[
 f_{uuu}=f_{uuv}=f_{uvv}=f_{vvv}=0
                  \quad\hbox{as polynomial identities}.       \tag{1}
\]

Then either grad f is a polynomial automorphism or f is jointly affine
along E. Under nonautomorphy, E must be Tang's unique affine plane W.
In particular:

- A genuinely quadratic restriction to a constant two-plane closes the
  unit-Hessian candidate in every degree of the other variables.
- Two distinct constant two-planes satisfying (1) force invertibility.
- For a hypothetical nonautomorphic f, the entire set of constant
  two-planes satisfying (1) is the singleton {W}; it already has all
  its vertical second derivatives zero.

Assume grad f nonautomorphic. Tang supplies W. Its projection rank
to the retained base (s,t) is zero, one or two, and constant affine
changes preserving that base preserve the degree bound along E.

In rank two, the cotangent form is

\[
 f=D(U,V)+sA(U,V)+tC(U,V),\qquad J(A,C)\in C^*.
\]

The E-degree bound makes A,C have ordinary total degree at most two.
The quadratic plane fact gives an inverse, contradicting nonautomorphy.

In rank one, write

\[
 f=D(t,U)+sA(t,U)+V C(t,U),\qquad J_{t,U}(A,C)\in C^*.
\]

Because (U,V) are the free E variables, the E-degree bound forces
deg_U C<=1. The first plane fact again gives an inverse, regardless
of the degree in t. Only rank zero remains, W=E, giving global joint
affineness in u,v. Uniqueness follows also from retained HC4MYHF0.

Condition (1) concerns the full family of parallel translates, not one
quadratic fibre or one direction with a vanishing third derivative.
The result places no bound on the total degree of f. It does not say
every affine-in-E unit-Hessian polynomial has invertible gradient:
that residual coefficient pair still carries the plane Keller problem.

## 3. All fibre-affine two-affine embeddings (HC4MYHT1)

Let Q be any nondegenerate quadratic polynomial on A5 and consider

\[
 j(s,t,u,v)=\bigl(s,t,b(s,t)+A(s,t)(u,v)^T\bigr),            \tag{2}
\]

where A is an arbitrary polynomial 3 by 2 matrix and b is an arbitrary
three-vector of base polynomials. Assume j is a closed embedding.
All degrees in s,t are unrestricted.

Equivalently, A has a polynomial left inverse over C[s,t]. Indeed a
left inverse recovers u,v from the three helpers minus b. Conversely,
differentiate a polynomial source recovery formula with respect to
(u,v) and evaluate at (u,v)=0; it supplies a left inverse of A over
the base ring. Generic matrix rank alone is insufficient.

If f=Q composed with j has unit actual source Hessian, then it has at
most one full critical point. Every nonautomorphic member has the
unit-ideal or singleton-zero plane endpoint already classified. Both
conditional endpoints occur in this fibre-affine class, so its
nonautomorphic unit-Hessian existence problem is exactly their
disjunction. No general two-affine or unrestricted HC4 closure follows.

Here is the all-degree proof. The potential has total free degree at
most two. HC4MYHT0 therefore forces the free plane E to be the actual
cotangent plane whenever the gradient is nonautomorphic. Alignment is
a conclusion here. Thus the quadratic term in u,v vanishes identically:

\[
                         A^T C A=0,                         \tag{3}
\]

where C is the constant Hessian block of Q on the three helper directions.
Its rank r is one, two or three, by ambient nondegeneracy.

For r=3, a nondegenerate three-space has no totally isotropic two-plane:
such a plane would lie in its orthogonal complement of dimension one.
Over C(s,t), A has rank two by its left inverse, contradicting (3).
Thus **every rank-three member has invertible gradient** under unit Hessian.

For r=2 normalize its quadratic part to X*Y, with Z its radical.
The two linear free parts of X and Y have product zero, by (3).
The source polynomial ring is a domain, so one entire linear free part
is zero globally. That helper is a base-only polynomial. For r=1 the
quadratic part is X^2 and its linear free part must be zero as well.
In either case the image is in the retained base-only graph family
of HC4MYHR1, which already allows all relative automorphisms of the
other two helpers and assumes no alignment. Its exact endpoint
classification gives at most one full critical point.

More specifically, rank one forces a coefficient zero at its canonical
base origin by HC4MYHO1; hence a nonautomorphic rank-one member has the
singleton-zero endpoint. Rank two permits either old endpoint. A
member with automorphic gradient always has one full critical point,
which completes the unconditional critical-count assertion.

The reduction in ranks one and two uses a constant linear combination
of helpers, not an unproved polynomial completion of a vector bundle.
No Quillen--Suslin theorem or affine-plane embedding rectification is
needed.

### Both conditional endpoints have fibre-affine realizations

For a unit-ideal Keller pair F,G, choose aF+bG=1, put
V=D(s,t)+uF+vG and R=-bu+av. The embedding

\[
                  j=(s,t,1,V-t^2-sR,R),\qquad
                  Q=sZ+t^2+XY                               \tag{4}
\]

is fibre-affine and has potential V. Recover V=t^2+Y+sZ and use
u=a(V-D)-GZ, v=b(V-D)+FZ. Its three-dimensional helper quadratic
block has rank two. This supplies no example of a unit-ideal Keller pair.

For a singleton-zero Keller pair write (F,G)=(s,t)C(s,t), with C an
invertible polynomial matrix from the retained oriented row theorem.
If (Y,Z)^T=C(u,v)^T, then

\[
                    j=(s,t,0,Y,Z),\qquad Q=X^2+sY+tZ        \tag{5}
\]

is a rank-one fibre-affine embedding with potential uF+vG and literal
polynomial source recovery. These are the two old conditional endpoints,
not a constructed nonautomorphic pair or a properness proof.

## 4. All three helpers can depend on the free variables

This exact unbounded family shows that (2) genuinely reaches beyond
having a base-only helper in the given constant ambient frame. Let
m,n be arbitrary positive integers, a=s^m, b=s^n, and set

\[
\begin{aligned}
 X&=(1-ab)^2u-b^2v,\\
 Y&=v-a^2u,\\
 Z&=(a-ba^2)u+bv.
\end{aligned}                                                   \tag{6}
\]

Then XY+Z^2=uv, so j=(s,t,X,Y,Z) pulls Q=st+XY+Z^2 back to st+uv,
with unit Hessian and invertible gradient. Its global source inverse is

\[
 u=X+2bZ-b^2Y,\qquad
 v=Y+2a(Z-bY)-a^2(X+2bZ-b^2Y).                            \tag{7}
\]

Each helper is affine in (u,v) and genuinely depends on them. Their
source degrees are respectively 2m+2n+1, 2m+1 and 2m+n+1, all distinct
and greater than one. No nonzero constant linear combination of the
helpers is affine-linear or base-only: its largest-degree term cannot
cancel and still depends on u. Thus exactly two independent affine
pullbacks are retained, and no third independent ambient affine
combination is base-only. The image has no affine equation.

The family also has a known ambient quadratic symmetry certificate,
being two elementary helper isometries. The new theorem classifies
all matrices A in (2); it does not require finding such a certificate.
The example illustrates the strict extension of the direct base-only
helper hypothesis, not a claim to escape the earlier symmetry saturation.

## Remaining scope and controls

For a nonautomorphic two-affine embedding outside the classified
fibre-affine family, either its potential is aligned and affine in
the free variables, or its free-variable degree is at least three.
Any argument allowing arbitrary nonlinear helpers must still impose
actual global source generation. A generic full-rank matrix, an
arbitrary substitution or a nonlinear source coordinate change cannot
silently replace that condition.

The [checker](verify_quadratic_plane.py) and
[receipt](receipt-quadratic-plane.json) retain small exact identities
separately from the unrestricted-degree plane lemma, the three Tang
positions and the complete fibre-affine classification. The
[literature context](sources/quadratic-plane-context.json) records
classical degree-two attribution only; both required plane facts are
proved above. Tang, uniqueness and the existing graph/endpoint results
are the mathematical inputs. Both plane endpoints, arbitrary nonlinear
two-affine embeddings and unrestricted HC4/JC2 remain open. No independent
replay, formal verification, external review or novelty priority is claimed.

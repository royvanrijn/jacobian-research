# Nonlinear image graphs and positive-weight quadratic symmetries

14 November 2026. Some nonlinear image equations admit a direct global
quadratic-preserving rectification. This note gives explicit inverses in
unrestricted degree, then proves that every smooth image invariant under
a compatible positive linear weighting is covered, or has a degenerate
critical point. Source parameterizations need not respect the weighting.

## 1. Graphs in an isotropic ideal (HC4MYHU0)

Work over C and write a nondegenerate centered ambient quadratic in
either of the following constant linear coordinate systems. Extra
coordinates U carry the nondegenerate residual quadratic q0.

### Nonisotropic normal

Let

\[
 q=R^2+\sum_{i=1}^k x_i y_i+q_0(U),\qquad
 X=\{R=c+h(x,U)\},\quad h\in(x_1,\ldots,x_k)\,C[x,U].     \tag{1}
\]

The variables U may occur in h with arbitrary degrees; the dual y_i
and the normal R do not occur there. Write h=sum_i x_i a_i(x,U),
which is a literal polynomial ideal certificate. The map

\[
 R'=R-h,\quad y_i'=y_i+2R a_i-h a_i,
 \quad x_i'=x_i,\quad U'=U                                \tag{2}
\]

is a polynomial automorphism preserving the **whole** q. Its inverse is

\[
 R=R'+h,\qquad y_i=y_i'-2R'a_i-h a_i,
 \qquad x=x',\ U=U'.                                     \tag{3}
\]

The h and a_i are unchanged because their variables are fixed. Expanding
(R-h)^2 and sum_i x_i y_i' cancels the cross terms using sum_i x_i a_i=h.
The inverse identities are equally direct; the Jacobian determinant is
one in the order (x,U,R,y). The image of X is exactly R'=c.

This works in every ambient dimension admitting (1), with arbitrary
polynomial degrees. In dimension five it sends every such closed A4
image, under **every source parameterization**, to a nonisotropic affine
hyperplane without changing the actual source potential. Its full source
orbit is the nondegenerate Q4 orbit with critical value lambda_*+c^2.
Nonautomorphic unit-Hessian existence is exactly the old singleton-zero
endpoint. This holds for each fixed graph (1), not merely for their union.

### Isotropic normal

Let

\[
 q=XY+\sum_{i=1}^k z_i w_i+q_0(U),\quad
 X_{\rm image}=\{X=c+h(z,Y,U)\},\quad
 h=\sum_i z_i a_i(z,Y,U).                                \tag{4}
\]

Here Y and U may occur arbitrarily in the coefficients; X and the
dual w_i may not. Then

\[
 X'=X-h,\quad w_i'=w_i+Y a_i,
 \quad Y'=Y,\ z'=z,\ U'=U                               \tag{5}
\]

preserves q and has inverse X=X'+h, w_i=w_i'-Y a_i. It sends the
image graph exactly to X'=c and has Jacobian one. In dimension five,
c=0 gives the rank-three restricted quadratic with a full critical line,
so no source parameterization can have unit Hessian. If c!=0 the
restricted potential is a polynomial coordinate with triangular
completion, and unit-Hessian existence is exactly the old unit-ideal
endpoint. Neither endpoint is constructed or resolved.

The hypotheses are invariantly described by a normal covector and
independent isotropic covectors orthogonal to it in the inverse-Hessian
metric, together with the displayed ideal certificate and independence
from the coordinates moved by the shear. They are sufficient conditions,
not a characterization of all quadratic symmetries or all rectifiable
embeddings. In particular the h dependence cannot be extended to the
dual coordinates in these formulas without another proof.

For each graph the canonical hyperplane parameterization composed with
the written inverse gives a genuine global graph embedding and its
polynomial inverse. Composing instead with the retained conditional
endpoint parameterizations gives the exact endpoint converses. Through
the retained ambient normalization Phi, these are genuine original
Meng--Yang embeddings as well. No nonlinear source Hessian invariance
or ambient embedding rectification theorem is being assumed.

## 2. Every compatible positive-weight invariant image (HC4MYHU1)

Let q be a nondegenerate quadratic form on A5, centered at its critical
point 0, and suppose a diagonal linear action

\[
 z_i\longmapsto a^{w_i}z_i,\qquad w_i>0,
 \qquad q(a^w z)=a^d q(z)                                \tag{6}
\]

has positive rational weights. Scale them to positive integers. Let
j:A4->A5 be a polynomial closed embedding whose image is invariant
under this action. There is no assumed equivariance of j in its actual
source variables.

The image is a smooth prime hypersurface E=0. It contains 0: each orbit
contracts there and the image is closed. Invariance of its principal
ideal makes E weighted homogeneous of some positive degree e. Indeed
E(a^w z) is a generator of the same principal prime ideal, hence a
nonzero scalar times E; comparison with one monomial gives a^e and
forces every monomial to have that weight. Smoothness at 0 gives a
nonzero linear part ell, supported exactly in the weight-e coordinate
space.

Let H=Hess(q) and nu=ell H inverse ell^T. There are exactly two cases.

**If nu=0**, the restriction of H to ker ell is degenerate. At the
source point mapping to 0, ambient grad q is zero, so the source
Hessian is DJ^T H DJ. It is the restricted quadratic form, hence
degenerate. Thus no actual source parameterization has unit Hessian.

**If nu!=0**, then 2e=d. To see this, q's matrix can pair only weight
w_i with weight d-w_i, and its inverse has the same pairing condition.
The inverse-metric norm of ell can therefore be nonzero only at twice
its weight. Split the nondegenerate weight-e block and pair the lower
and higher weight spaces by constant linear algebra, respecting weights:

\[
 q=R^2+\sum_i x_i y_i+q_0(U),\quad
 \mathrm{wt}(R)=\mathrm{wt}(U_j)=e,\quad
 \mathrm{wt}(x_i)<e<\mathrm{wt}(y_i).                     \tag{7}
\]

Normalize E's linear part to R. Every nonlinear monomial of E has
weight e and at least two positive-weight factors. Hence it can use
only variables of weight less than e. Thus

\[
                         E=R-h(x),\qquad h(0)=0.          \tag{8}
\]

The lower weight covectors are mutually isotropic and orthogonal to R.
Equation (8) is exactly (1), with c=0. A compatible ideal certificate is

\[
             a_i=\frac{w_i}{e}\,\partial h/\partial x_i,
                 \qquad h=\sum_i x_i a_i,                \tag{9}
\]

by weighted Euler's identity. Formula (2) is then homogeneous in every
coordinate and commutes with the positive action. It rectifies the image
to R'=0, preserves q globally, and has the explicit polynomial inverse.
No formal-to-polynomial extension or general group-action linearization
theorem is used.

It follows that **every** unit-Hessian parameterization of a compatible
positive-weight invariant image has exactly one full critical point,
with critical value lambda_*. Nonautomorphic members retain exactly the
singleton-zero endpoint of HC4MYHG2. The unit-ideal endpoint does not
occur in this image class because the ambient critical point lies on
the image. For every fixed image in the nonisotropic case, its whole
source orbit is exactly the nondegenerate Q4 orbit; the corresponding
conditional singleton converses are retained. Both the source action
and the actual source parameterization may be nonlinear.

This is a recognition theorem for the quadratic-symmetry saturation:
the required ambient certificate is constructed, not assumed. It works
in all polynomial degrees and for every supplied compatible positive
weight vector. It does not assert that arbitrary invariant images for
nonpositive weights, arbitrary nonlinear C* actions, or all embeddings
admit such a weighting. Given equations in a fixed frame, checking (6)
and weighted homogeneity is a finite set of exact linear weight identities;
no search over all constant frames is asserted.

## 3. A positive linear action on the actual source is stronger

If the actual source variables themselves have positive rational weights
and the actual potential f-lambda_* is weighted homogeneous, unit Hessian
forces its gradient to be a polynomial automorphism. This elementary
fact holds in every dimension. Write its weighted degree as d and its
source weights as alpha_i. Nondegenerate Hessian excludes alpha_i>=d,
since the corresponding derivative would be constant or zero. Thus the
gradient has positive output weights d-alpha_i and vanishes at 0.

The Hessian at 0 is invertible, so there is a unique formal inverse of
the gradient. Equivariance and uniqueness force its i-th coordinate to
have weight alpha_i when the outputs have weights d-alpha_j. There
are only finitely many monomials of a fixed weight for strictly positive
weights. The formal inverse is therefore polynomial, and its two formal
composition identities hold as polynomial identities globally.

In particular a parameterization equivariant for positive **linear**
actions on both ambient space and the actual source closes unconditionally.
Image equivariance alone does not give that conclusion: choosing homogeneous
source coordinates can be a nonlinear change, and unit Hessian is not
invariant under such a change. Section 2 deliberately retains the
singleton-zero endpoint for arbitrary source parameterizations.

## 4. Why both compatibility and strict positivity matter

For q=R^2+xy+zw, the smooth graph R=xy is invariant under positive
weights (wt R,wt x,wt y,wt z,wt w)=(2,1,1,1,1), but q is not homogeneous
for them. Its restricted potential is (xy)^2+xy+zw. It has a critical
point at the origin of value 0 and a critical curve xy=-1/2,z=w=0
of value -1/4. It cannot be isometric to an affine hyperplane restriction,
whose entire critical set has at most one critical value.

If positivity is dropped, the graph R=xy-1 and q above are both
homogeneous for weights (0,1,-1,1,-1), of degree zero. This invariant
image omits 0; its restricted potential has critical values 1 and 3/4,
again including a critical curve. Neither example is a unit-Hessian
counterexample; the critical curves actually forbid unit Hessian after
any source automorphism by HC4MYHG0. They show that invariance alone
does not justify the rectification or orbit classification claimed here.

## Evidence and remaining scope

The [checker](verify_positive_weight_images.py) and
[receipt](receipt-positive-weight-images.json) keep literal isometry,
inverse, weight and critical-point controls separate from the written
unrestricted-degree proof. The [primary context](sources/positive-weight-images-context.json)
records Lamy--Venereau's elementary quadric automorphisms without
importing a tame-generation or classification theorem. The actual
proofs here are direct polynomial identities and the retained hyperplane
and source-orbit endpoint results.

Nonlinear image equations outside these isotropic-ideal graphs and
positive compatible weight systems, both plane endpoints and unrestricted
HC4/JC2 remain open. No formal verification, independent replay, external
review or novelty priority is asserted.

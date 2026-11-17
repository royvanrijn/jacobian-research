# An ambient quadratic normal form and every affine graph direction

14 November 2026. The retained HC5T1 frame straightens the **entire**
scaled Meng--Yang potential, not just one source slice. After an explicit
polynomial ambient automorphism Phi preserving x,y,

\[
 \Psi=\lambda_*+Q_5\circ\Phi,\qquad
 Q_5(x,y,T,Z,W)=LT^2+NyZ+2NxW,\qquad
 \lambda_*=-M^2/(4L).                                      \tag{1}
\]

The quadratic form Q5 is nondegenerate when LN is nonzero. This
function identity does not conjugate gradient maps: Phi is nonlinear.

Every polynomial graph over any affine four-plane, in any constant
ambient direction, has an invertible gradient under actual unit Hessian
when substituted in Q5. Pulling such a **parameterized** graph back
by Phi gives a new all-degree family of genuine Meng--Yang embeddings
with the same conclusion. Neither graph degree nor direction is bounded.

## 1. A literal ambient automorphism (HC4MYHM0)

Put v=xy, u=1+v, and

\[
 a=1-3v+6v^2,\qquad b=y^3(6v^2+15v+10),\qquad
 H=6v^2+15v+4.
\]

The retained Bezout identity is a*u^3-x^3*b=1. In native physical
coordinates (x,y,p,q,r), set

\[
 t=A=u^3p+3xu^2q-x^3r,\quad \rho=q,\quad \sigma=-bp+ar.
\]

This is a determinant-one relative linear change. Let C0 be the
explicit determinant-one normalizer from
[HC4MYOR1](../hc4-relative-row-orbit/PROOF.md#1-the-row-and-an-explicit-polynomial-normalization):

\[
 C_0=\begin{pmatrix}
 (3v^2+6v+2)/2 & x^2(2v+3)\\
 (3y^2/2)(18v^4+63v^3+69v^2+21v-1)&
 36v^5+108v^4+87v^3+3v^2-3v+1
 \end{pmatrix}.
\]

Define (z,w)^T=C0^(-1)(rho,sigma)^T. The exact native row identity,
already retained in the HC5T1 construction, is

\[
                           B=y z+2xw+y^2H t.                \tag{2}
\]

Now put

\[
                   T=t+M/(2L),\qquad Z=z+yH t,\qquad W=w.
\]

These are a relative shear and a constant translation. They give
(1) immediately, and Phi=(x,y,T,Z,W) is an actual polynomial
automorphism with Jacobian determinant one.

For completeness its inverse is completely polynomial. Given
(x,y,T,Z,W), set

\[
 t=T-M/(2L),\quad z=Z-yH t,\quad
 (\rho,\sigma)^T=C_0(z,W)^T,
\]

and recover

\[
 p=a(t-3xu^2\rho)+x^3\sigma,\qquad q=\rho,\qquad
 r=b(t-3xu^2\rho)+u^3\sigma.                                \tag{3}
\]

All formulas are global; x,y,u are never inverted. This normal form
is a consequence of the retained explicit frames, not an imported
Morse lemma or a claim of new priority for those frames. It preserves
functions and closed embeddings. It does not preserve the Hessian of
an arbitrary function under a nonlinear source change, nor transfer
a marked collision of gradient maps automatically.

## 2. Two elementary all-degree HC4 criteria (HC4MYHM1)

The only structural input in this section is the retained reproduced
Tang theorem: under hypothetical nonautomorphy a unit-Hessian
four-variable polynomial has a constant jointly affine two-plane W.
Adapt constant source coordinates (z1,z2,w1,w2) to W. Then

\[
                         f=D(z)+w^T K(z),\qquad J(K)\ne0
                         \text{ constant}.                 \tag{4}
\]

A plane Keller pair with one nonconstant affine-linear component is
triangular in constant coordinates: make that component z1; the unit
Jacobian makes the derivative of the other component in z2 a nonzero
constant. Its inverse and the cotangent gradient inverse are polynomial.

**Quadratic plus one square.** Let q be any polynomial of degree at
most two in four variables and c a nonzero constant. Then

\[
            f=q+cR^2,\quad \det\operatorname{Hess}f\in C^\times
                   \quad\Longrightarrow\quad\nabla f
                   \text{ is a polynomial automorphism}.  \tag{5}
\]

Suppose otherwise and use (4). Regard R as a polynomial in w over
C[z]. If its w degree m is at least two, R^2 has degree 2m>2,
which q cannot cancel to make f affine in w. If its w degree is one,
write R=R0(z)+a1(z)w1+a2(z)w2. Vanishing of the w Hessian of f says
that every a_i*a_j is constant, since the w Hessian of q is constant.
In particular a_i^2 constant makes each a_i a field constant. If the
w degree is zero take both a_i=0.

Consequently the coefficient pair K in (4) has the form

\[
                            K(z)=Az+d+vR_0(z)               \tag{6}
\]

for a constant matrix A and constant vectors d,v. If v=0 it is
affine-linear. Otherwise a nonzero constant row annihilating v
makes one component of an invertible linear target change of K
affine-linear; it cannot be constant because J(K) is a unit.
Thus K is triangular and the gradient invertible, a contradiction.
This proves (5), without any degree bound on R or nondegeneracy
assumption on the quadratic background q.

**One nondegenerate ternary quadratic trace.** If f has unit Hessian
and its restriction to one actual affine source hyperplane is a
nondegenerate quadratic polynomial in three variables, then its
gradient is invertible. Here nondegenerate means rank three of the
quadratic part; affine and constant terms are allowed.

Translate the hyperplane to s=0 and again suppose W exists under
nonautomorphy. If W is contained in ker ds, the ternary quadratic
would be jointly affine along a two-plane. This would be a totally
isotropic two-plane for a nondegenerate three-dimensional bilinear
form, impossible because W would lie in its orthogonal complement
of dimension one.

Otherwise W meets ker ds in a line spanned by u, and has a second
constant vector v with ds(v)=1. Joint affineness gives
D_u(D_u f)=D_v(D_u f)=0. On s=0, D_u f is a nonconstant affine-linear
polynomial, since the ternary quadratic part is nondegenerate.
The transport equation D_v(D_u f)=0 uniquely extends this trace
as the same affine-linear polynomial in translated variables
z-s*v_z. Hence D_u f is globally nonconstant affine-linear.
Because u belongs to W, it is a component of the Keller coefficient
pair in (4). The pair is triangular, again a contradiction.

This is an elementary proof of the needed quadratic-trace case of
the retained trace criteria; no additional Hodge or embedding
classification theorem is imported.

## 3. Every affine graph direction in a nondegenerate quadratic ambient space (HC4MYHM2)

Let Q be any nondegenerate quadratic polynomial on A5. For any
constant point a, constant 5-by-4 matrix B and constant vector v
such that (B,v) is invertible, and any polynomial R on A4, put

\[
                         j(s)=a+Bs+vR(s).                   \tag{7}
\]

This is a genuine polynomial closed embedding: projection along v
recovers s affinely. For f=Q composed with j, unit actual four-variable
Hessian implies a polynomial gradient inverse, in every degree.

Let beta be the bilinear form of Q's quadratic part. If beta(v,v)
is nonzero, write

\[
              f=q(s)+\ell(s)R(s)+cR(s)^2
                =q'(s)+c\bigl(R(s)+\ell(s)/(2c)\bigr)^2,
\]

where q,q' are quadratic and ell is affine-linear. Apply (5).

If beta(v,v)=0, the R coefficient ell(s)=beta(v,a+Bs) plus the
possible original linear term is nonconstant affine-linear. On its
zero hyperplane the R term disappears. The remaining quadratic
part is nondegenerate of rank three: it identifies with the
nondegenerate quotient v^perp/Cv. Equivalently, B*ker(ell_linear)
is a complement of Cv in v^perp, so the restriction has no radical.
Apply the ternary quadratic trace criterion. These two cases
exhaust every constant graph direction and prove the theorem.

Combining (1), (3) and (7), every parameterization

\[
                   i(s)=\Phi^{-1}\bigl(a+Bs+vR(s)\bigr)      \tag{8}
\]

is an explicit genuine Meng--Yang embedding. If its actual pulled-back
potential has unit Hessian, its gradient is invertible. The source
parameter s is the actual affine graph parameter; its degree is free.

For an existing proposed embedding i, a useful sufficient test is:
some constant rank-four projection of the five pullbacks Phi composed
with i is an affine-linear source coordinate system. Then an affine
source change puts it in (8), and the theorem closes the candidate.
The projection must be affine-linear in the actual source coordinates;
being only a nonlinear polynomial coordinate system is insufficient.

## 4. Why full source reparameterization remains different

The image being a graph is not enough for Section 3 if the actual
source parameterization is nonlinear. This distinction is necessary,
not a missing justification to fill by Hessian invariance.

The fixed normal graph T=0 gives

\[
                         \lambda_*+NyZ+2NxW,
\]

a nondegenerate four-variable quadratic. Its full polynomial source
orbit has exactly the distinct singleton-zero Keller endpoint retained
in HC4MYHG2: a nonautomorphic unit-Hessian member would give that
endpoint, and every hypothetical such endpoint is realized in this
orbit. The fixed isotropic-direction graph x=1/(2N) gives

\[
                         \lambda_*+LT^2+NyZ+W,
\]

a polynomial coordinate. Its full source orbit is the coordinate
orbit, with exactly the unit-ideal plane Keller endpoint already
retained in HC4MYGU2/3 and HC4MYHG2. Thus allowing arbitrary source
automorphisms on even these two fixed normal graphs already returns
the two unresolved endpoints. Neither is closed here.

The transverse square equation from HC4MYHL2 is also not a universal
exclusion. When T=0 and D is any constant c with c!=lambda_star,
it has the constant solution h^2=(c-lambda_star)/L over C, for every
plane Keller K. This satisfies the square equation only; it supplies
neither the remaining normal jets nor a global embedding. A further
square-condition campaign alone cannot settle that class.

The remaining general embedding problem is now equivalently a
quadratic representation with source generators: find or exclude
five polynomials x,y,T,Z,W generating the full four-variable source
ring with f-lambda_star=LT^2+NyZ+2NxW and actual unit Hessian.
The linear-projection condition of Section 3 is sufficient to close
such a realization, not automatic from ring generation. The universal
SL2 row condition, arbitrary embeddings lacking that parameterized
graph form, both distinct plane endpoints and unrestricted HC4/JC2
remain open.

The [exact checker](verify_ambient_quadric.py) and
[receipt](receipt-ambient-quadric.json) retain finite algebraic controls
separately from the written all-degree graph argument. The frame
identities, two graph-direction cases, triangular inverses, exact
embedding recovery and source-parameterization limitations are explicit.
No new external classification, formal verification, independent
replay, external review or novelty priority is claimed.

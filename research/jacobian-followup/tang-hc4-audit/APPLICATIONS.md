# Quantitative and recognition consequences of the structural theorem

These consequences use the reproduced external theorem
[HCTANG2](STRUCTURAL_REVIEW.md). They are registered together as HCTANG3.
Their proofs are given here; no novelty priority or external review is
claimed. JC2 and unrestricted HC4 remain unresolved.

## A finite constant-direction certificate

Let S=Hess f have nonzero constant determinant in four variables. Define

\[
 Z_f=\{[v]\in\mathbb P^3:\ v^T\operatorname{adj}(S(X))v=0
                 \text{ as a polynomial in }X\}.
\]

**Theorem.** If Z_f is empty over the algebraic closure, grad f is a
polynomial automorphism. Equivalently, a nonautomorphic gradient necessarily
admits a nonzero **constant** covector v such that adding
t(v.X)^2/2 to f preserves its Hessian determinant for every t.

Indeed, the structural theorem puts a nonautomorphic f into constant-linear
cotangent form R(b)+z A(b)+w B(b). Its Hessian is

\[
 S=\begin{pmatrix}C&DG^T\\DG&0\end{pmatrix},
\]

with DG invertible over the polynomial ring. The base/base block of S^-1
is zero. Every nonzero constant covector supported on the base is therefore
isotropic for S^-1. Transforming it back by the constant source matrix
gives the required constant v. The rank-one determinant identity is

\[
 \det(S+t vv^T)=\det S+t v^T\operatorname{adj}(S)v.
\]

Thus emptiness of the projective scheme is an unconditional certificate of
invertibility. Nonemptiness does not certify noninvertibility. Emptiness
can be checked by extracting all coefficient quadrics in v and saturating
by the irrelevant ideal, or by checking the four affine projective charts
v_i=1 for unit ideals. A bounded search for rational v does not suffice.
The established HC4MR1/HC4MR2 auxiliary-pencil equivalence is consistent
with this result; the external structural theorem supplies recognition
for the nonautomorphic case, which those older statements did not prove.

## Potential degree at most 101

More generally, if every plane Keller pair whose maximum component total
degree is at most m is invertible, then HC4 holds for potentials of degree
at most m+1. In the nonautomorphic cotangent alternative, A=f_z and B=f_w
have degrees at most m; a constant linear source change does not raise
total degree. The plane theorem then contradicts nonautomorphy.

Using only [Moh's established degree-100 range](https://www.math.purdue.edu/~ttm/jacobian)
gives **HC4 through
potential degree 101**. No unrestricted plane theorem is used. This is
a consequence of Tang's degree-preserving structural statement, not of
an abstract equivalence alone. We deliberately do not require the
repository's larger plane frontier or its additional provenance imports.
The characteristic-zero assertion descends from C: embed a finitely
generated coefficient subfield in C and use the unique formal inverse,
whose coefficients are already in that subfield. If it truncates over C,
it truncates there as well.

There is a useful endpoint refinement. If the highest homogeneous part
of f is a nonzero scalar times a proper power H^k, k>=2, then the
cotangent coefficient degrees are at most deg(f)-2. Indeed, a polynomial
proper power that is jointly affine in the two fibres must be independent
of them: its degree in either fibre would otherwise be at least two.
Thus the degree-(D-1) parts of A,B vanish in cotangent coordinates.
**Potentials of degree 102 with a proper-power top are therefore also
invertible**, by the same degree-100 plane import. This refinement uses
the structure of the top, not an improvement of the plane theorem.

## Meng–Yang marked-collision graphs through degree 48

Use exactly the scaled family in
[the retained graph note](../../archive/non-elliptic/HC4_MENG_YANG_GRAPH_OBSTRUCTIONS.md):

\[
 u=1+xy,\quad A=u^3p+3xu^2q-x^3r,
\]
\[
 B=y^2u(4+3xy)p+(y+3xy^2(4+3xy))q+(2x-3x^2y)r,
 \quad\Psi=L A^2+M A+N B.
\]

Substitute any polynomial graph r=R(x,y,p,q) of degree d. The degree of
A_R and B_R is at most max(7,d+3). Hence

\[
 \deg\Psi_R\le\max(14,2d+6).
\]

The full support check gives the same bound directly: the coefficients
of r^0,r^1,r^2 in Psi have maximal base degrees 14,10,6. Thus the
weighted support maximum is max(14,d+10,2d+6)=max(14,2d+6), since the
middle value is the arithmetic mean of the other two. This is a proof
for all d, rather than interpolation from sampled graph degrees.

For the v2 member (L,M,N)=(1,13,2), the marked points are
P+=(1,-3/2,0,0,0) and P-=(-1,3/2,0,0,0), with common gradient
(0,0,-1/2,0,0). If

\[
 R(1,-3/2,0,0)=R(-1,3/2,0,0)=0,
\]

both belong to the graph. The chain rule gives
grad(Psi composed with i_R)=Di_R^T grad Psi. Since the fifth component
of the common gradient is zero, the two induced four-variable gradients
coincide **regardless of the first derivatives of R**.

For d<=47 the substituted potential degree is at most 100 and lies
within the degree-101 HC4 theorem. At d=48 the degree is exactly 102:
the uniquely highest contribution is x^6 R_48^2, a square. The endpoint
refinement applies. In both cases the gradient must be injective if
its Hessian determinant is a nonzero constant, contradicting the collision.
Therefore **no marked-collision graph of degree at most 48 can have
nonzero constant Hessian determinant**.

This closes the collision-preserving quintic-graph gate and all the same
gates through degree 48. It does not prove that every unmarked graph
of these degrees has nonconstant determinant, nor does it settle
arbitrary-degree graphs or nonlinear mixed/coisotropic embeddings.
At d=49 this argument permits potential degree 104, with cotangent
coefficients of degree at most 102, and gives no exclusion. The unique
formal constant-determinant graph series may still exist; it cannot
truncate to degree at most 48 while retaining
the marked collision.

A coarser statement also applies beyond graph parameterizations. For any
polynomial map i from affine four-space to the five-variable source with
degree e, the pullback has degree at most 14e. If distinct a+,a- map to
P+,P- and Di(a+)^T g=Di(a-)^T g for their common gradient g, the chain
rule transfers the collision. Thus no such pullback of degree e<=7 can
have nonzero constant Hessian determinant. Here g=(0,0,-1/2,0,0), so
it suffices that the p-coordinate of i be affine, or more generally have
equal first derivatives at a+ and a-. This does not assume that arbitrary
nonlinear maps automatically transfer the collision. Partial Legendre or
Schur eliminations producing rational expressions are not covered by
this polynomial-pullback degree bound.

The [checker](verify_structural.py) verifies the source pins, universal
line-frame algebra, weight formulas, determinant-update identity, complete
support-degree calculation and marked collision. The all-degree structural
proof is the written review, not an output of that finite checker.

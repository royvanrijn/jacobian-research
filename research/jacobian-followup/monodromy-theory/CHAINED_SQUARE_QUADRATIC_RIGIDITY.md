# The entire chained-square quadratic family is excluded to full rank

The first of the four arbitrary-matrix primitive quadratic families is
now closed in every homogeneous scalar degree. A unique projective
base point forces the scalar factors to be cones. The first trace
then reduces all remaining multipliers to one linear-factor power;
the second trace has an unavoidable nonzero coefficient.

This excludes the whole family, not just selected output matrices.
The three other families and the absolute action problem remain open.

## 1. A unique base point forces scalar cones

Work over algebraically closed characteristic-zero constants and put
\[
 K=(u^2,ux,uy+x^2,uz+y^2),\qquad V=A K,\quad A\in GL_4.
\]
The projective base locus is the single point a=[0:0:0:1].
The projective map phi=[V] is an isomorphism between its source and
target affine charts by KMFI142.

Suppose a nonzero homogeneous scalar F gives a normalized full-rank
input, with h=lambda*u and f=-lambda*u*F as in KMFI125.
KMFI150 gives u not dividing F, and KMFI117 excludes constant F.
For every irreducible factor D of F, KMFI144 gives both
D(p+tV) divisible by D and
\[
                       D(V)=c_D u^{\deg D}D,\quad c_D\ne0. \tag{1}
\]

**KMFI151.** Every such scalar F is independent of z.
Every nonlinear irreducible factor is a projective cone with vertex a;
every linear factor belongs to the span of u,x.

First, phi cannot fix a whole surface pointwise.
Otherwise intersect that surface with a plane avoiding a and take
an irreducible curve component C. The quadratic coordinate sections
of V have no common zero on C, so their map pulls back O(1) to O_C(2).
The map is the identity there, which pulls back O(1) to O_C(1).
Their degrees are respectively 2degC and degC, a contradiction.
Equality with the identity extends from a dense open set because
both morphisms on the surface outside a have separated target.
No irreducibility theorem for a general plane section is needed:
any irreducible curve component suffices.

Let degD>=2. At a general point of D the joining line of p,phi(p)
is therefore defined and lies on D by the divisor homotopy.
The selected-line-family argument of KMFI140 makes this family a
curve and gives a unique selected line through a general point.
Thus phi preserves each general selected projective line L.
Its restriction is generically injective on L's u!=0 part by KMFI142,
so it has degree1. The restricted quadratic sections consequently
have one common base point. Since a is the only base point, L passes
through a. The surface is the closure of a family of lines through a,
hence a cone with that vertex. Its irreducible equation is independent
of z: translation p->p+t*a preserves the cone, so differentiation
would make D divide D_z, of smaller degree, forcing D_z=0.

For linear D=d0*u+d1*x+d2*y+d3*z, compare coefficients in(1).
The x^2 and y^2 coefficients of D(V) are D*A*e2,D*A*e3.
They must vanish. Its uy,uz coefficients are those same expressions,
whereas the right side has coefficients c_D*d2,c_D*d3.
Since c_D!=0, d2=d3=0. Thus linear factors too are independent of z.
All factors of F have now been covered, proving the claim.
This cone reduction uses no bound on individual factor degree.

The finite base locus is essential to the pointwise-fixed-surface
step. For example D=uz-xy and V=(u^2,ux,uy,2uz-xy) have
D(V)=2u^2*D and V=u*p on D, while the base locus consists of two
lines. This is a geometric boundary control, not a normalized
full-rank example or a new action.

## 2. The traces close every arbitrary output matrix

**KMFI152.** For every A in GL4 and every homogeneous nonzero F,
V=A*(u^2,ux,uy+x^2,uz+y^2) gives no normalized input of generic rankDX4.

Assume such an input exists. By KMFI151, F=F(u,x,y), nonconstant
and coprime to u. Use
\[
 N=DV+V\,d\log F-\lambda uI,\qquad N^4=0.                  \tag{2}
\]
The trace equation is V(F)+F*divV=4lambda*u*F.
Its coefficient in z gives exactly
\[
 u(A_{03}F_u+A_{13}F_x+A_{23}F_y)=-A_{03}F.               \tag{3}
\]
If A03!=0, choose u and two linear invariants of the constant
derivative on the left. Then(3) becomes u*partial_u F=-F.
An Euler derivative in one variable has only nonnegative integer
eigenvalues on polynomials, so this is impossible.
Thus A03=0, and
\[
                            A13*F_x+A23*F_y=0.           \tag{4}
\]

If A13=A23=0, invertibility gives A33!=0 and Ae3=A33*e3.
Since F_z=0, column z of N is (A33-lambda)*u*e3.
Nilpotency forces lambda=A33. But then
X=(f,F(V-lambda*u*p)) is independent of z, contradicting full rank.

Otherwise(4) makes F a homogeneous polynomial in the two linear
variables u,t=A23*x-A13*y. Over algebraically closed constants it
splits into linear factors in u,t. It has a non-u factor since it is
nonconstant and coprime to u. KMFI151 requires that factor to belong
to the span of u,x. If A13!=0 its nonzero t coefficient would give
a nonzero y coefficient, impossible. Hence A13=0,A23!=0, and
\[
                                  F=F(u,x).              \tag{5}
\]

We also need A02!=0. If A02=0, then A03=0 makes V_u=uL for a
linear L. Because u does not divide F, the logarithmic form is regular
at the generic point of u0. Row u of N there is
(L|u0,0,0,0). Nilpotency makes its eigenvalue L|u0 vanish.
Thus V_u=a*u^2 with a!=0, a matched chart already excluded to
primitive full rank by KMFI141. This proves A02!=0.
This elementary invariant-critical-hyperplane argument does not assume
a boundary factor in F.

Every linear factor D=d0*u+d1*x of(5) has d1!=0 and satisfies
\[
                                d0*A02+d1*A12=0           \tag{6}
\]
by(1). Since A02!=0, all factors are proportional. Thus
F=constant*(x+t0*u)^s for one integer s>=1.

The simultaneous constant source/output change x_new=x+t0*u
preserves this canonical family. Indeed its inverse substitution gives
\[
 K_{\rm old}=B K_{\rm new},\qquad
 B=\begin{pmatrix}1&0&0&0\\-t0&1&0&0\\t0^2&-2t0&1&0\\0&0&0&1\end{pmatrix}.
\]
The new output matrix is C*A*B, where C performs x_new=x+t0*u.
This is an allowed conjugacy, with arbitrary other parameters retained;
its first-row column2 coefficient remains the same nonzero A02.
After renaming coordinates, F=x^s, and the factor identity forces
V_x=c*u*x with c!=0. The matrix consequently has shape
\[
 A=\begin{pmatrix}
 a&b&d&0\\0&c&0&0\\e&f&g&h\\i&j&k&l
 \end{pmatrix},\qquad c d\ne0.                            \tag{7}
\]
For eta=s*dx/x, the entire trace of(2) is
\[
 \operatorname{tr}N=
 (2a+(s+1)c+g+l-4\lambda)u+b x+(d+2h)y.                  \tag{8}
\]
Hence b=0. The entire rational expression tr(N^2) is polynomial,
and its x^2 coefficient is
\[
                 [x^2]\operatorname{tr}(N^2)
                    =b^2+2(s+2)cd=2(s+2)cd\ne0.           \tag{9}
\]
This contradicts nilpotency. Equations(8)--(9) are whole symbolic
identities, valid for all parameters; no generic parameter division
or scalar degree bound has been used. The proof is complete.

The full-rank hypothesis is sharp. For every positive integer s,
lambda=s+2, A=diag(1,lambda,lambda,lambda) and F=u^s give
f=-lambda*u^(s+1),
g=(-(s+1)u^(s+2),0,lambda*u^s*x^2,lambda*u^s*y^2).
Its P=Dg-p*df is strictly lower triangular and generically regular
nilpotent, and the whole sweep has Jacobian gamma^4. But X is
independent of z and has generic rank3. These are retained
normalized examples, not full-rank inputs.

## 3. Evidence and the unresolved action objective

[Helper](chained_square.py) · [checker](verify_chained_square.py) ·
[model](chained_square_model.json) · [receipt](chained_square_receipt.json).

The written geometry and polynomial proof are unbounded in scalar
degree, factors and all entries of A. Exact controls check canonical
base-point identities, the fixed-surface boundary example, linear-factor
coefficient identities, both complete traces, the allowed conjugacy,
and actual normalized rank3 examples with entire sweep Jacobian.

Only K_A has been closed here. The mixed chain K_B, square shear K_C
and product shear K_D retain arbitrary output matrices and nonconstant
multipliers coprime to u, with individual factor degrees1..3.
Their actual generic groups, higher primitive degrees, nonhomogeneous
inputs and an ordinary-affine exceptional11 realization or absolute
action obstruction remain open. No new occurrence action is claimed.

The previous goal turn was progress: KMFI149--150 excluded every
critical-boundary scalar factor including mixed products. This turn
closes the whole arbitrary-matrix chained-square direction family.
Formal verification, independent whole-proof replay and external review
are not claimed. The full faithful-transitive Keller monodromy
classification remains the active objective.

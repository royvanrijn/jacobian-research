# Dual-substitution graph theorem and the two-divisor obstruction

13 November 2026. Elementary result: HC4MYGT1. Structural corollary:
HC4MYGT2. These extend the standard-chart HC4MYGAP1/HC4MYGALL1 results
to the two triangular parameterizations and all base-independent plane
Keller substitutions in the README. All polynomial
degrees and coefficients are unrestricted. Throughout the field has
characteristic zero and L,N,a,b are nonzero constants.

Use the same five-variable family

\[
 u=1+xy,\quad A=u^3p+3xu^2q-x^3r,
\]
\[
 B=y^2u(4+3xy)p+(y+3xy^2(4+3xy))q+(2x-3x^2y)r,
 \quad\Psi=L A^2+M A+N B.
\]

## Transformation laws and the actual hypothesis

Both dual maps in the README are polynomial automorphisms over K[x,y],
with constant Jacobian determinant ab. Their inverses solve their two
triangular equations successively. Their images with r=R are still
polynomial graph hypersurfaces in the original four coordinates. However,
the gradient/Hessian assertion depends on the chosen parameters.

For any four-variable F and polynomial coordinate map phi,

\[
 \operatorname{Hess}(F\circ\phi)=
 D\phi^T(\operatorname{Hess}F\circ\phi)D\phi
 +\sum_i (F_i\circ\phi)\operatorname{Hess}\phi_i.
\]

The extra term is essential. For example, composing half the sum of four
squares with phi=(x,y,P,Q+P^2), of determinant one, gives Hessian
determinant 1+2Q+2P^2. Thus the new theorem assumes a nonzero constant
determinant for the **actual reparameterized potential**; it does not
infer that condition from the standard graph theorem.

## A useful completed-square identity

The defining formulas give

\[
 x^2 B+(2-3xy)A=u(u+1)p+2x(2u+1)q.
\]

Work temporarily with x inverted and set

\[
 K_0=M-\frac{N(2-3xy)}{x^2},\quad
 \alpha=\frac{Nu(u+1)}{x^2},\quad
 \Gamma=\frac{2N(2u+1)}x,\quad U=A+\frac{K_0}{2L}.
\]

After either parameterization,

\[
 f=L U^2+\alpha p+\Gamma q-\frac{K_0^2}{4L}.
\]

All these expressions are polynomials in y,P,Q over K(x). In y, Gamma
is affine and alpha is quadratic. Localization here is only an algebraic
identity tool; polynomiality of the original R will be used to prove
that the mixed derivatives below cannot vanish.

An affine direction means a constant vector v with D_v^2 f=0 identically.
A jointly affine vector plane must consist entirely of affine directions
and have zero mixed second derivatives as well.

## No affine direction invisible to x can move y

Take v=(0,1,s,t), using scaling to make its nonzero y component one.
Use constant linear coordinates P'=P-sy, Q'=Q-ty, so D_v is differentiation
in y with x,P',Q' fixed. R remains polynomial after this change.

For the first triangular order, p=aP+G(x,y), q=bQ+H(x,y,P),

\[
 f=L U^2+b\Gamma Q'+E(x,y,P'),\qquad
 U_{Q'}=3b x u^2-x^3R_{Q'}.
\]

The terms absorbed into E are independent of Q'. If f were affine in y,
its Q' derivative would also be affine in y. Hence U*U_Q' has y degree
at most one. Since U_Q' is nonzero, polynomial degree additivity forces
U_Q' to be independent of y: either U has degree zero, or it has degree
one and its slope must have zero Q' derivative. But

\[
 U_{Q'y}=x^2(6b u-xR_{Q'y})\ne0,
\]

as the expression in parentheses restricts to 6b at x=0. Contradiction.

For the opposite order, p=aP+H(x,y,Q), q=bQ+G(x,y),

\[
 f=L U^2+a\alpha P'+E(x,y,Q'),\qquad
 U_{P'}=a u^3-x^3R_{P'}.
\]

If f were affine in y, U*U_P' would have y degree at most two, because
alpha is quadratic. However

\[
 U_{P'yy}=x^2(6a u-xR_{P'yy})\ne0.
\]

Thus deg_y U_P'>=2, and deg_y U>=deg_y U_P'>=2, making the degree of
their product at least four. Contradiction. These mixed-derivative
arguments treat every constant s,t and every polynomial R; they are
not checks of selected directions.

It follows in both triangular orders that every affine direction with
zero x component must also have zero y component.

## First order: the remaining pure-dual direction is impossible

Let v=(0,0,c,d) be affine in the first order. On x=y=0 the potential is
L(aP+G(0,0))^2+M(aP+G(0,0)). Restriction commutes with D_v. Its second
derivative is 2La^2 c^2, so c=0. If d is nonzero the Q direction would
have to be affine. Yet

\[
 A_Q=x(3b u^2-x^2R_Q)\ne0.
\]

The completed-square formula makes f have precisely twice the positive
Q degree of A, since the remaining Q term is linear. It is not affine
in Q. Thus the first order has no nonzero affine direction invisible to
x, and consequently no jointly affine two-plane.

## Opposite order: transport between x=0 and u=0

The preceding stronger conclusion about individual directions is false
for the opposite order. A two-plane is still excluded.

For a pure-dual affine direction v=(0,0,c,d), restrict to x=0 and write

\[
 p_0=aP+H_0(y,Q),\quad q_0=bQ+G_0(y),\quad
 f_0=L p_0^2+(M+4Ny^2)p_0+Ny q_0.
\]

Along the v line, any nonconstant polynomial p0 has degree m>=1, giving
f0 degree 2m, since q0 is linear along that line. Thus p0 is constant
along v and a*c+d*H0_Q=0 identically. A nonzero v has d nonzero, so

\[
 H_0(y,Q)=hQ+h_0(y),\qquad h=-ac/d\in K.
\]

If that trace condition fails there are already no nonzero affine
directions invisible to x. If it holds, put T=aP+hQ and S=Q, an
invertible constant linear change. The only possible such direction
is now the S axis. At x=0,

\[
 p_0=T+h_0(y),\qquad (f_S)|_{x=0}=Nby.
\]

Suppose a jointly affine two-plane exists. It must contain the S axis,
because it meets x=0 in a nonzero affine direction, and all of those
directions were just confined to that line. Choose its other generator
v=(1,eta,tau,sigma) in x,y,T,S coordinates. Joint affineness gives
f_SS=0 and D_v(f_S)=0. Thus f_S=C(x,y,T), and

\[
 C_x+\eta C_y+\tau C_T=0,
 \qquad C(0,y,T)=Nby.
\]

The constant-flow identity therefore forces the full polynomial

\[
 C(x,y,T)=Nb(y-\eta x).
\]

Now restrict to u=0, that is y=-1/x, in the domain K[x,x^-1,T,S].
Here A=-x^3R and the same elimination identity becomes

\[
 f=L A^2+(M-5N/x^2)A+\frac{2N}{x}(bS+G(x,-1/x)).
\]

This is a polynomial in S over K(x,T). Its S derivative was required
to be independent of S, so the square term forces A to be independent
of S. Consequently C=2Nb/x on this divisor. The transported identity
instead says C=Nb(-1/x-eta*x). Equality would give eta*x^2+3=0
identically in the free variable x, impossible for constant eta.
This excludes the last jointly affine plane.

This step uses two different divisors and the vanishing mixed derivative
on the whole plane. The existence of one individual affine direction
does not invalidate or establish a cotangent plane.

As a positive control, take the opposite order with

\[
 H=3b x(xy-1)Q,\quad G=0,\quad R=3b y^2u^2Q.
\]

Then A=a u^3P and f_Q=Nb*y*(3u^2-2), so Q really is an affine direction.
Its coefficient is not the transported linear form above. This example
shows why the first-order proof about individual directions cannot just
be copied to the opposite order; the two-divisor plane argument is needed.

## Every base-independent plane Keller substitution

Now let p=p(P,Q), q=q(P,Q) have nonzero constant Jacobian determinant,
and let r=R(x,y,P,Q) be arbitrary polynomial. This hypothesis is on the
Jacobian matrix, and does not require a polynomial inverse map. We prove
that there is no nonzero affine direction with zero x component.

First suppose the direction moves y, and normalize it to
v=(0,1,c,d). Write delta for c*d/dP+d*d/dQ. The x=0 slice is

\[
 f_0=Lp^2+Mp+4Ny^2p+Nyq,
\]

where p,q are independent of y. From (d/dy+delta)^2 f0=0, its y^2 and
y coefficients give

\[
 \delta^2p=0,\qquad \delta^2q=-16\delta p.
\]

The constant coefficient is then

\[
 2L(\delta p)^2+8Np+2N\delta q=0.
\]

Applying delta gives -24N*delta p=0. Consequently delta p=0 and
delta q=-4p. If delta=0 this would force p=0, contrary to its nonzero
Jacobian. If delta is nonzero, its components c,d generate the unit
ideal. The polynomially invertible Jacobian matrix D(p,q) implies that
delta p,delta q generate the same unit ideal: multiplication by its
polynomial inverse preserves that ideal. But here that ideal is (p),
which is not the unit ideal because p is nonconstant. Contradiction.

For a nonzero pure-dual affine direction delta, the y^0 coefficient of
delta^2 f0 says delta^2(Lp^2+Mp)=0. Polynomial degree along a translation
line forces delta p=0. The same Jacobian-matrix argument makes delta q
a nonzero constant, since its principal ideal must be the unit ideal.
Now along that direction A has positive degree: its derivative is
x*(3u^2*delta q-x^2*delta R), nonzero modulo x after removing the factor x.
The completed-square identity has only affine additional dependence on
that direction, so f has twice A's positive degree and is not affine.

Thus this third class has no nonzero affine direction invisible to x,
and hence no jointly affine two-plane. This includes every base-independent
polynomial dual automorphism, regardless of its word length. For example,
p=P+(Q+P^2)^2, q=Q+P^2 is covered, even though it is not one triangular
change in the displayed parameter coordinates. The argument does not
extend automatically when the substitution coefficients depend on x,y.

## The structural consequence

All three classes have no jointly affine constant two-plane. A nonautomorphic
four-variable unit-Hessian gradient would have exactly such a plane by
the external structural alternative reproduced as HCTANG2. Therefore
every unit-Hessian potential in these three graph families has
a polynomial gradient inverse. The elementary no-plane theorem itself
does not use Tang. Polynomial invertibility descends from C to any
characteristic-zero coefficient field by the unique formal inverse.

## What transfers the marked collision

The ambient v2 common gradient is g=(0,0,-1/2,0,0). For any polynomial
parameterization i, its pullback at a marked point is

\[
 Di^T g=-\tfrac12\nabla_{x,y,P,Q}p.
\]

Containing both marked points is therefore insufficient by itself after
a nonlinear parameter change. In the first triangular order, the two
pulled-back gradients coincide exactly when G_x and G_y agree at the
two base points (1,-3/2),(-1,3/2). This is automatic for affine G.
H and R do not enter that condition. In the opposite order the condition
is equality of H_x,H_y,H_Q at the two parameter points, where
Q=-G(x,y)/b and P=-H(x,y,Q)/a. An affine H makes it automatic.

When the graph also contains the marked points and these jet conditions
hold, its Hessian determinant cannot be a nonzero constant, by the
invertibility theorem. Without the jet match the ambient collision need
not transfer; the general invertibility assertion still holds.
For the third class the condition is equality of the two-component
gradient of p(P,Q) at the dual parameter values. A collision of the
dual map alone does not assert equality of those gradients.

The [exact checker](verify.py) verifies the elimination and completed-square
identities, immutable mixed jets, trace and two-divisor algebra, the positive
and negative controls, and the marked-covector formulas. The base-independent
proof uses exact coefficient identities and a polynomial-Jacobian ideal
argument. Polynomial-degree and plane-transport exhaustion are written
arguments, not finite searches. This covers one base-dependent triangular
change in either order and arbitrary base-independent Keller substitutions.
Base-dependent alternating words and nonlinear transformations of x,y
need a separate argument.

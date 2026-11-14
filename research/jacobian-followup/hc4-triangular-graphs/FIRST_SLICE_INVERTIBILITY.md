# All-degree invertibility from the first Meng--Yang slice alone

> The later [automorphic-trace theorem](AUTOMORPHIC_TRACE.md)
> (`HC4MYGV0/1/2`) now closes every y-dependent dual automorphism over
> K[x,y], and every substitution with an automorphic first trace over
> K[y]. Its Hamiltonian pencil argument completes the generic-coordinate
> continuation below. The broader merely-Keller encoding in Section 6
> remains a reduction to JC2; it does not have a certified trace inverse.

14 November 2026. **HC4MYGF1** is an elementary classification of any
affine fibre plane compatible with the prescribed slice. **HC4MYGF2**
uses the reproduced [Tang structural alternative](../tang-hc4-audit/STRUCTURAL_REVIEW.md)
to prove gradient invertibility without assuming that a plane exists.
This removes the second-divisor hypothesis from the earlier
[invertibility extension](TWO_DIVISOR_EXTENSION.md). Its stronger
matching-trace no-plane theorem remains valid.

## Statement

Work over a characteristic-zero field K with LN nonzero. Let p0(P,Q),
q0(P,Q) have a nonzero constant plane Jacobian determinant. Consider
**any** polynomial f(x,y,P,Q) with the complete hyperplane trace

\[
 f(0,y,P,Q)=Lp_0^2+Mp_0+4Ny^2p_0+Nyq_0.             \tag{1}
\]

**HC4MYGF1.** If f has a jointly affine constant vector two-plane W,
then, after an invertible constant linear change of the four variables,

\[
 f=NkYS+XU(Y,T)+V(Y,T),\qquad k\in K^\times.        \tag{2}
\]

If its Hessian determinant is a nonzero constant, necessarily
U=alpha T+beta(Y), with alpha a nonzero constant. The gradient then
has the explicit polynomial inverse below. This elementary implication
uses no Tang theorem or plane Jacobian conjecture.

**HC4MYGF2.** Every polynomial satisfying (1) with nonzero constant
Hessian determinant has a polynomially invertible gradient. Indeed,
HCTANG2 makes a hypothetical nonautomorphic gradient admit W; HC4MYGF1
then contradicts nonautomorphy. There is no degree bound. The theorem
holds over every characteristic-zero coefficient field by the same
finite-coefficient-field and unique-formal-inverse descent as the
[standard graph theorem](../hc4-polynomial-graphs/PROOF.md).

In particular, for the scaled Meng--Yang potential Psi, all substitutions

\[
 p=p_0+xE(x,y,P,Q),\qquad q=q_0+xH(x,y,P,Q),
 \qquad r=R(x,y,P,Q)                               \tag{3}
\]

are covered, for arbitrary polynomials E,H,R. No Jacobian or inverse-map
condition is imposed on the full base-dependent substitution. Only its
x=0 dual trace must be a Keller pair independent of y. Both second
divisor traces may be arbitrary. This closes the different-second-trace
continuation proposed in HC4MYGD1 whenever its first trace has (1).

## 1. What the slice forces on a plane

For completeness, here is the slice argument independently of the
second-divisor calculation. Write delta=c partial_P+d partial_Q.
If an affine direction invisible to x moves y, normalize it to
partial_y+delta. Comparing y coefficients in its second derivative
of (1) gives

\[
 \delta^2p_0=0,\qquad \delta^2q_0=-16\delta p_0,
 \qquad 2L(\delta p_0)^2+8Np_0+2N\delta q_0=0.
\]

Applying delta to the last equation forces -24N delta p0=0. Thus
delta p0=0 and delta q0=-4p0. If delta is zero this contradicts
nonconstancy of p0. Otherwise multiplication by the polynomially
invertible trace Jacobian preserves the component ideal of (c,d),
which is the unit ideal. Its image (0,-4p0) generates a proper ideal.
This is impossible. Every x-invisible affine direction is purely dual.

For a nonzero purely dual affine direction, compare at y=0. If p0
had positive degree m along that translation direction, Lp0^2+Mp0
would have degree 2m, contradicting its vanishing second derivative.
Hence delta p0=0. The same component-ideal argument implies
delta q0=k, a nonzero constant. Choose constant linear dual coordinates
T0,S0 with delta=partial_S0. The trace Keller determinant now forces

\[
 p_0=aT_0+d_0,\qquad q_0=kS_0+h(T_0),\qquad ak\ne0. \tag{4}
\]

All x-invisible affine directions lie on this one line: two independent
ones would annihilate the nonconstant p0. Therefore W meets ker(dx)
in exactly that line. It contains partial_S0 and a second vector
(1,eta,tau,sigma) in x,y,T0,S0 coordinates.

## 2. Transport produces a linear component, not an obstruction

Joint affineness gives f_S0S0=0 and D_v f_S0=0. Thus f_S0=C(x,y,T0),
with

\[
 C_x+\eta C_y+\tau C_{T_0}=0,\qquad C(0,y,T_0)=Nky.
\]

The unique polynomial solution is C=Nk(y-eta x). Set

\[
 X=x,\quad Y=y-\eta x,\quad T=T_0-\tau x,
 \quad S=S_0-\sigma x.
\]

These coordinates make W the span of partial_X and partial_S, and
f_S=NkY. Integration gives f=NkYS+g(X,Y,T). Joint affineness also
gives g_XX=0, proving (2). This proves plane exhaustion; no second
divisor or assumption that every affine direction is absent is needed.

## 3. A unit Hessian forces the remaining coefficient to be triangular

In coordinate order (X,S,Y,T), the Hessian has zero upper-left
two-by-two block. Its off-diagonal block is

\[
 J=\begin{pmatrix}U_Y&U_T\\Nk&0\end{pmatrix}.
\]

Regardless of the bottom-right block,

\[
 \det\operatorname{Hess}f=(\det J)^2=(Nk)^2 U_T^2. \tag{5}
\]

If this is a nonzero constant, polynomiality makes U_T=alpha a
nonzero constant. Hence U=alpha T+beta(Y). For output coordinates
(z_X,z_S,z_Y,z_T) of the gradient, recover successively

\[
\begin{aligned}
 Y&=z_S/(Nk),\\
 T&=(z_X-\beta(Y))/\alpha,\\
 X&=(z_T-V_T(Y,T))/\alpha,\\
 S&=(z_Y-X\beta'(Y)-V_Y(Y,T))/(Nk).                \tag{6}
\end{aligned}
\]

All denominators are nonzero field constants. Inverting the constant
linear source transformation converts this into a polynomial inverse
for the original gradient: gradients transform by the transpose of
that transformation. This proves HC4MYGF1 and the contradiction needed
for HC4MYGF2. No unproved plane properness implication enters (5)-(6).

## 4. Arbitrary word lengths and the boundary that remains

Every elementary dual shear congruent to identity modulo x preserves
the first trace. Every finite composition does too, regardless of
degree or word length. Composing on either side with a base-independent
Keller substitution makes the trace that same Keller pair. Thus all
such base-dependent alternating words satisfy (3), even if their
u=1+xy=0 traces differ. The earlier second-trace square equation remains
a valid necessary condition for a plane, but is no longer an unresolved
unit-Hessian counterexample gate in this class: a surviving plane has
the explicit inverse (6).

The trace hypothesis is substantial. For example p0=P and
q0=Q-4yP give a Keller pair over K[y] but are not independent of y.
The slice becomes LP^2+MP+NyQ and has the affine direction partial_y,
so Step 1 does not extend to that larger trace class. The checker
retains this failure control. It is not a counterexample to HC4.

Words with genuinely y-dependent first traces, nonlinear base changes,
other graph charts, and general mixed or coisotropic reductions remain
outside this theorem. Marked-point collision transfer still requires
the two pulled-back p covectors to agree, as in the earlier theorem.
No assertion about every such word, general HC4 or JC2 is made.

## 5. An actual polynomial unit-Hessian member of (3)

The conclusion is invertibility. The larger class really can contain
a cotangent plane and a unit Hessian, as the following exact control shows.
For a nonzero constant kappa put

\[
 F=LP^2+MP+4Ny^2P+NyQ+\kappa xP,
 \qquad\det\operatorname{Hess}F=(N\kappa)^2.
\]

There are polynomial p,q,R with p(0)=P, q(0)=Q and Psi(x,y,p,q,R)=F.
Here is a reproducible construction using the retained
[HC5T1 completion and normalizer C0](../hc4-relative-row-orbit/PROOF.md).
Write v=xy, u=1+v, pbez=1-3v+6v^2 and
qbez=y^3(6v^2+15v+10). Let

\[
 t=P,\quad z=Q-y(6v^2+15v)P,\quad w=\frac{\kappa}{2N}P,
 \qquad (\rho,\sigma)^T=C_0(z,w)^T.
\]

Then set

\[
\begin{aligned}
 p&=p_{\rm bez}(t-3xu^2\rho)+x^3\sigma,\\
 q&=\rho,\\
 R&=q_{\rm bez}(t-3xu^2\rho)+u^3\sigma.
\end{aligned}
\]

The Bezout identity pbez u^3-x^3 qbez=1 gives A=t=P. The normalized
B row is y z+2x w+y^2(6v^2+15v+4)t, hence
B=4y^2P+yQ+(kappa/N)xP as required. At x=0, C0 has first row (1,0),
so both first traces are the identity. On u=0 instead the Q coefficient
of q is -1/2, so this control is outside the older matching-two-divisor
theorem.

The construction is also a polynomial closed embedding into the
five-variable ambient space: the HC5T1 change and C0 are invertible,
and in their coordinates the image is the graph w=kappa t/(2N).
Its chosen p,q projection need not be invertible. The control prevents
mistaking the new result for a no-plane or unit-Hessian nonexistence
theorem. Its gradient inverse follows directly from (6).

## 6. Sharp trace boundary: the larger substitution class contains JC2

**HC4MYGF3.** If (1) is relaxed to permit a y-dependent first trace
pair that is Keller over K[y], the resulting class of arbitrary
polynomial substitutions (p,q,R) contains the cotangent lift of every
plane Keller map. Over C, universal gradient invertibility in this larger
class is equivalent to JC2 (using HCTANG2 in the JC2-to-HC4 direction).
This is a boundary reduction, not a solution of the larger class.

Here is the explicit encoding. Given any plane Keller pair
(F(P,Q),G(P,Q)) with Jacobian j in K^times, use the same invertible
HC5T1 frames as in Section 5, but now take

\[
 t=F,\quad z=G-y(6v^2+15v+4)F,\quad w=F/2.
\]

Apply the same formulas for rho,sigma,p,q,R. They give

\[
 A=F,\qquad B=xF+yG,\qquad
 f=LF^2+MF+N(xF+yG),                              \tag{7}
\]

with first traces p0=F and q0=G-4yF. Their P,Q Jacobian is j,
a nonzero constant, so this is a Keller pair over K[y]. Its
y dependence is precisely the hypothesis not covered by HC4MYGF2.

In fibre/base order (x,y,P,Q), equation (7) has

\[
 \det\operatorname{Hess}f=N^4j^2.
\]

Set h=LF^2+MF and replace the fibre coordinates (x,y) by

\[
 (x',y')^T=\nabla_{P,Q}h+N D(F,G)^T(x,y)^T.
\]

This is a polynomial source automorphism: its inverse uses the
polynomial matrix (N D(F,G)^T)^(-1) with constant denominator Nj.
After this source change the gradient is the product of (NF,NG)
with an identity map, up to coordinate order. Therefore grad(f) is
a polynomial automorphism if and only if (F,G) is one. An inverse
of a product map restricts at zero in its identity output coordinates
to an inverse of the plane factor. Thus universal invertibility in
this larger substitution class implies JC2. Conversely JC2 and
HCTANG2 imply HC4 and hence invertibility for every member.

This encoding does **not** establish the same equivalence for genuine
dual automorphisms or polynomially embedded graphs. Its full dual
Jacobian is

\[
 j\bigl(p_{\rm bez}(C_0)_{11}+x^3/2\bigr),
\]

a nonconstant polynomial. Nor is a polynomial closed embedding
asserted for an unknown plane Keller pair. The construction is a
polynomial parameterization, with a Keller first trace, inside the
larger substitution class explicitly allowed here. Requiring an
everywhere invertible dual map or an actual embedded graph remains
a separate mathematical condition. In particular, this boundary
reduction does not retire the y-dependent dual-automorphism route.

For that genuine route there is a concrete additional gate. A dual
automorphism over K[x,y] restricts at x=0 to an automorphism over K[y].
Over K(y), its potential trace

\[
 f_0=Lp_0^2+(M+4Ny^2)p_0+Nyq_0
\]

is therefore a coordinate: (p0,f0/(Ny)) is obtained from (p0,q0) by
an explicit triangular shear. Its inverse recovers
q0=f0/(Ny)-Lp0^2/(Ny)-(M+4Ny^2)p0/(Ny).
Thus any claimed surviving plane in this route must be compatible
with this generic-coordinate property, in addition to the whole-plane
equations. A Keller Jacobian alone does not certify that property.
This is a necessary continuation gate, not a classification of the
y-dependent traces or a new JC2 argument.

## Reproduction and assurance

```sh
timeout 60s .venv/bin/python research/jacobian-followup/hc4-triangular-graphs/verify_first_slice.py
```

The [checker](verify_first_slice.py) derives the slice and transport
identities, verifies the determinant and both inverse compositions,
and constructs the actual polynomial unit-Hessian control. It also
checks a different-trace alternating word and the y-dependent trace
failure control, and verifies the universal plane-Keller encoding and
its nonconstant full dual Jacobian. Unbounded direction exhaustion and word length are
written arguments. [Replay receipt](receipt-first-slice.json).
The structural classification is the separately reproduced external
HCTANG2 input, not certified by these local identities. No formal
verification, external review or novelty priority is claimed.

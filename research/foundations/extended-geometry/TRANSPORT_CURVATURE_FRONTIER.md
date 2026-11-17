# Transport curvature and the low-dimensional frontier

## Status

This note tests the proposed common transport mechanism behind plane JC,
the historical Hall-shell route to nonhomogeneous binary GVC, and `HC(4)`.
It proves a scoped all-degree theorem for rank-one binary leading profiles.
It does **not** prove plane JC or `HC(4)`, classify rank-two leading forms, or
promote Hall shells to one algebraic transport object.  Binary GVC is now
proved independently by the
[Hall-envelope theorem](BINARY_GVC_ENVELOPE_CLOSURE.md).

The determinant identities hold over every characteristic-zero field.  The
normalization and affine-latitude conclusions use an algebraically closed
characteristic-zero field.  The symbolic identities have both a structural
checker and an independent direct-Hessian replay.

## 1. The common elementary frame

The known JC(3), Long/GMC(3), and homogeneous GVC(3) constructions all use
pullbacks of

\[
M(\xi,\eta)=
\begin{pmatrix}1+\xi\eta&\eta\\ \xi&1\end{pmatrix}
=E_{12}(\eta)E_{21}(\xi).
\]

Its first column is the unimodular pair

\[
(U,V)=(1+\xi\eta,\xi),\qquad U-\eta V=1,
\]

with affine ratio `U/V=eta+V^(-1)`.  The pullbacks are

\[
(\xi,\eta)=(x,y),\quad (x,a),\quad (x^2,1),
\]

respectively.  Long/GMC(3) and homogeneous GVC(3) also have the same
normalized adjacent jet `[0:1]`; the endpoint-contact exponent changes only
the nonzero scalar moment.  See
[`verified/MARKED_ROOT_MODEL.md`](../verified/MARKED_ROOT_MODEL.md),
[`HOPF_LIFT_CLASSIFICATION.md`](HOPF_LIFT_CLASSIFICATION.md), and
[`THREE_VARIABLE_HOMOGENEOUS_GVC_COUNTEREXAMPLE.md`](THREE_VARIABLE_HOMOGENEOUS_GVC_COUNTEREXAMPLE.md).

Flatness of `M^(-1)dM` is tautological, so a flat `SL_2` frame is not the
wanted invariant.  The first intrinsic differential gate visible in the
rank-one Hessian branch is the projective curvature of its coefficient map.

## 2. Universal dominant-profile formula

Let

\[
X=(x_1,x_2),\qquad U=(u_1,u_2),\qquad
c(X)=\binom{a(X)}{b(X)},\qquad L=c(X)^TU.
\]

Write

\[
J=D_Xc,\qquad
S_U=u_1\operatorname{Hess}(a)+u_2\operatorname{Hess}(b),
\qquad w=J^{-1}c
\]

on `det J != 0`.  For `d>=2` and a nonzero scalar `lambda`,

\[
\boxed{
\det\operatorname{Hess}_{X,U}(\lambda L^d)
=(\lambda d)^4L^{4d-5}(\det J)^2
\left((2d-1)L-(d-1)w^TS_Uw\right).
}
\tag{2.1}
\]

To prove this, put

\[
v=\binom{J^TU}{c},\qquad
K=\begin{pmatrix}S_U&J^T\\J&0\end{pmatrix}.
\]

Then

\[
\operatorname{Hess}(\lambda L^d)
=\lambda dL^{d-2}\left(LK+(d-1)vv^T\right),
\]

while

\[
K^{-1}=\begin{pmatrix}
0&J^{-1}\\J^{-T}&-J^{-T}S_UJ^{-1}
\end{pmatrix},
\quad
\det K=(\det J)^2,
\quad
v^TK^{-1}v=2L-w^TS_Uw.
\]

The matrix determinant lemma gives (2.1).

If the determinant vanishes, coefficient comparison gives

\[
\operatorname{Hess}(c_i)(w,w)=\frac{2d-1}{d-1}c_i,
\qquad
\boxed{\nabla_ww=-\frac d{d-1}w.}
\tag{2.2}
\]

Thus the projective curvature

\[
\kappa(c)=\det(w,\nabla_ww)
\]

vanishes.  A dominant polynomial `c:A^2->A^2` cannot satisfy (2.2): after
reparametrization an integral curve is an affine line `X=p+zv`, but
`D_wc=c` gives

\[
z\frac d{dz}c(p+zv)=-\frac{d-1}{d}c(p+zv),
\]

which has no nonzero polynomial solution.  Hence a constant-amplitude
dominant rank-one profile never has zero four-variable Hessian determinant.

## 3. Arbitrary amplitude forces rank collapse

Consider

\[
G_d(X,U)=\lambda(X)(a(X)u_1+b(X)u_2)^d\ne0.
\]

Adjoin `q` with `q^d=lambda` and put `c_tilde=(qa,qb)`.  If
`D c_tilde` had generic rank two, the preceding argument would apply over
that finite algebraic extension.  Along a general radial integral curve,
reparametrized as `X(z)=p+zv`, the polynomial endpoint coefficients

\[
g_0=(qa)^d=\lambda a^d,
\qquad
g_d=(qb)^d=\lambda b^d
\]

would satisfy

\[
z\frac d{dz}g_i=-(d-1)g_i.
\]

The Euler operator on a polynomial has only nonnegative integral
eigenvalues, contradicting the choice of a point where one endpoint is
nonzero.  Therefore

\[
\operatorname{rank}D\widetilde c\le1.
\]

In characteristic zero, generic differential rank equals the transcendence
degree of the generated function field.  Hence the coefficient algebra

\[
B_d=k[\lambda a^d,\lambda a^{d-1}b,\ldots,\lambda b^d]
\]

satisfies

\[
\boxed{\operatorname{trdeg}_kB_d\le1.}
\tag{3.1}
\]

On a chart where the unscaled `c` is dominant, the alternative is the exact
negative-eigenfunction equation

\[
w(\lambda)=-d\lambda,
\]

and all coefficients `lambda*a^(d-j)*b^j` are `w`-invariant.  This is the
rank-one version of “latitude plus compensator”.

## 4. The collapsed algebra is a polynomial latitude

Let `B` be a finitely generated domain with

\[
k\subsetneq B\subset k[x,y],\qquad \operatorname{trdeg}_kB=1.
\]

Its integral closure `Bbar` in `Frac(B)` is contained in `k[x,y]`, since an
element integral over `B` is integral over the normal ring `k[x,y]`.  The
normal affine curve `Spec(Bbar)` is dominated by a general affine line in
`A^2`, so it is rational.  Also `Bbar^*=k^*`, because both a unit and its
inverse lie in `k[x,y]`.  A normal rational affine curve with no nonconstant
units is `P^1` with exactly one point removed.  Therefore

\[
\boxed{\overline B=k[h]}
\tag{4.1}
\]

for a nonconstant polynomial `h in k[x,y]`.

Apply this to `B_d` and put `R=k[h]`.  Let

\[
g_j=\lambda a^{d-j}b^j\in R.
\]

If `g_0!=0`, write `g_1/g_0=p/q` with coprime `p,q in R`.  Since
`g_d=g_0p^d/q^d` lies in the PID `R`, one has `q^d|g_0`.  Taking

\[
A=q,\qquad B=p,\qquad \mu=g_0/q^d
\]

and restoring the binomial coefficients gives

\[
\boxed{
G_d=\mu(h)(A(h)u_1+B(h)u_2)^d.
}
\tag{4.2}
\]

The endpoint cases are immediate.

## 5. The bordered-Hessian curvature

Define

\[
\mathcal R(h)
=\nabla h^T\operatorname{adj}(\operatorname{Hess}h)\nabla h
=h_x^2h_{yy}-2h_xh_yh_{xy}+h_y^2h_{xx}.
\tag{5.1}
\]

For an arbitrary `Phi=Phi(h,u_1,u_2)`, direct block expansion gives

\[
\begin{aligned}
\det\operatorname{Hess}_{x,y,u_1,u_2}\Phi
={}&\Phi_h^2\det(\operatorname{Hess}h)
       \det(\operatorname{Hess}_{u_1,u_2}\Phi)\\
 &+\Phi_h\mathcal R(h)
       \det(\operatorname{Hess}_{h,u_1,u_2}\Phi).
\end{aligned}
\tag{5.2}
\]

Put

\[
L=A(h)u_1+B(h)u_2,
\qquad W=AB'-BA'.
\]

For `Phi=mu(h)L^d`,

\[
\det\operatorname{Hess}_{u_1,u_2}\Phi=0,
\]

and

\[
\det\operatorname{Hess}_{h,u_1,u_2}(\mu L^d)
=-d^3(d-1)\mu^3W^2L^{3d-4}.
\]

Consequently

\[
\boxed{
\det\operatorname{Hess}_{x,y,u_1,u_2}(\mu L^d)
=-d^3(d-1)\mu^3\mathcal R(h)W^2L^{4d-5}
 (\mu'L+d\mu L_h).
}
\tag{5.3}
\]

Only the projective Wronskian, the first source variation, and the curvature
of `h` remain; all second derivatives of `mu,A,B` cancel.

## 6. Zero curvature means an actual affine latitude

For nonconstant `h in k[x,y]`,

\[
\boxed{
\mathcal R(h)=0
\quad\Longleftrightarrow\quad
h=f(\alpha x+\beta y)
}
\tag{6.1}
\]

for some univariate `f` and `(alpha,beta)!=(0,0)`.

The reverse implication is immediate.  For the forward implication, on
`h_x!=0` use the tangent field

\[
T=h_y\partial_x-h_x\partial_y.
\]

A direct calculation gives

\[
T\left(\frac{h_y}{h_x}\right)=-\frac{\mathcal R(h)}{h_x^2}.
\]

Choose a generic value `c`, so `h=c` is reduced and smooth.  Every
irreducible component then has constant tangent direction and is an affine
line.  Distinct components of this smooth fiber are disjoint, hence the
lines are parallel.  After a linear coordinate change,

\[
h-c=\gamma\prod_i(x-r_i),
\]

with no additional factors or multiplicities.  Thus `h in k[x]` in those
coordinates, proving (6.1).

Combining (3.1), (4.2), (5.3), and (6.1) proves:

> **Rank-one no-virtual-latitude theorem.**  Let
>
> \[
> G_d=\lambda(x,y)(a(x,y)u_1+b(x,y)u_2)^d\ne0,
> \qquad d\ge2.
> \]
>
> If `det Hess(G_d)=0`, its coefficient algebra has transcendence degree at
> most one.  After normalization, write it as in (4.2).  Then either:
>
> 1. `W=0`, so the projective line `[A:B]` is constant;
> 2. `partial_h(mu*L^d)=0`, so the whole profile is source-independent; or
> 3. `h=f(alpha*x+beta*y)`, so the moving profile factors through an actual
>    affine linear latitude.
>
> A zero-Hessian rank-one top cannot carry a genuinely two-dimensional or
> nonlinear virtual latitude.

Rank collapse alone is not enough.  For

\[
c=(y,xy^2),\qquad \lambda=x^2,
\]

the coefficient algebra is `k[xy]`, but `R(xy)=-2xy` and

\[
\det\operatorname{Hess}\bigl(x^2(yt+xy^2m)^2\bigr)
=32x^8y^8(t+mxy)^3(t+2mxy)\ne0.
\]

## 7. Consequence for `HC(4)`

Let

\[
\Psi(X,U)=u_1P(X)+u_2Q(X)+H(X)+G_d(X,U)
          +\text{lower dual-degree terms},
\]

where `G_d` is the highest nonlinear dual-degree term.  The dual-degree
`4d-4` part of `det Hess(Psi)` is exactly `det Hess(G_d)`: every determinant
term uses exactly four dual differentiations, and replacing any top-layer
entry by a lower layer strictly lowers the dual degree.

Thus every constant-Hessian potential whose highest nonlinear dual layer is
a rank-one binary profile passes the theorem above, independently of its
cotangent background, source-only term, and lower dual layers.  The remaining
rank-one problem is reduced to a constant projective line or one affine
source coordinate.  Rank at least two remains open.

## 8. Cohn normal-form calibration

The degree-at-most-two-entry classification of
Chapovskyi--Kozachok--Petravchuk,
[*Decomposition of matrices from `SL_2(K[x,y])`*](https://arxiv.org/abs/2412.03688),
reduces the non-elementary cases, under the equivalences stated there, to
rows of the following forms.  The main theorem above does not depend on this
external classification.

For

\[
c=(\gamma+xy,x^k),\qquad \gamma\ne0,\quad k\ge2,
\]

one obtains

\[
\kappa(c)=
\frac{\gamma k(k-3)+(k-1)(k-2)xy}{k^3}\ne0.
\]

Its radial derivation is triangular on monomials with nonnegative diagonal
weights, so `w(lambda)=-d*lambda` has no nonzero polynomial solution.

For

\[
c=(x^2,\delta+x\psi(y)),\qquad \delta\ne0,\quad \psi'\ne0,
\]

one obtains

\[
\kappa(c)=
-\frac{(2\delta+x\psi)^2\psi''+2\delta x(\psi')^2}
 {8x(\psi')^3}\ne0.
\]

Vanishing would force first `psi''=0` and then `psi'=0`.  The highest
`x`-coefficient of a hypothetical polynomial negative eigenfunction would
satisfy

\[
\psi\lambda_N'+(N+2d)\psi'\lambda_N=0,
\]

hence `lambda_N=C*psi^(-(N+2d))`, again impossible.

Therefore neither Cohn normal-form row can be the dominant highest nonlinear
rank-one dual profile of a four-variable constant-Hessian potential, even
with arbitrary polynomial amplitude and lower layers.  This is a scoped
all-degree exclusion, not a classification of arbitrary `HC(4)` potentials.

On the moment side, the first row has strict one-sided torus weight after
pole clearing.  For the second row, the minimal second-contact condition on
`t^2+rho*x*y=1` forces

\[
\psi(y)=\frac{\rho\delta}{2}y,
\]

and

\[
1-(1-\rho xy)(1+\rho xy/2)^2
=\frac{\rho^2}{4}(xy)^2(3+\rho xy),
\]

again producing a one-sided torus profile.  Higher-contact and multi-profile
variants remain open.

## 9. Implications and remaining frontier

The common hierarchy is now

\[
\text{algebraic exposure}
\longrightarrow
\text{rank collapse}
\longrightarrow
\text{affine-latitude curvature}
\longrightarrow
\text{global holonomy}.
\]

- The historical Hall-shell route to binary GVC is still blocked at exposure:
  its prime-dependent shells have not been promoted to one fixed coefficient
  algebra.  The Hall-envelope proof bypasses this route.
- In the rank-one top sector, `HC(4)` now reaches an actual affine latitude;
  the next task is lower-layer synchronization after making that coordinate
  equal to `x`.
- Plane JC already has a genuine boundary curve; its unresolved stage is
  global endpoint transport and holonomy rather than local curvature.

This replaces the overbroad “flat residual line” proposal by a checkable
rank/curvature/holonomy programme.

## 10. Exact replay and assurance boundary

Run

```bash
python3 scripts/verify_transport_skeleton.py
python3 scripts/verify_transport_curvature_frontier.py
python3 scripts/audit_transport_curvature_frontier_independent.py
```

The first checker replays the elementary frame, exact JC marked root, both
polynomialization identities, and finite instances of the adjacent jet.  The
second checks the structural determinant identities, curvature formulas,
one-latitude formula, Cohn calibrations, and nonlinear-latitude
countercalibration.  The independent checker expands ordinary four-variable
Hessians directly for unrelated degree-two and degree-three profiles and does
not use the block derivation.

The normalization lemma, affine-line integral-curve contradiction, and
zero-curvature classification are mathematical proofs rather than bounded
searches.  The checkers replay their algebraic identities and calibrations;
they are not formal verification or external review.

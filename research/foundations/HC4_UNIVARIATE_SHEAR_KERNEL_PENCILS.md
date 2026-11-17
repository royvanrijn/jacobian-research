# Polynomial shear kernels and univariate reduction for `HC4`

## Status

After `HC4RSD3` closes every affine-in-the-reduced-variables kernel, this
note treats the first degree-unbounded nonlinear stratum, allowing dependence
on both transverse variables.

> **Theorem `HC4RSD4` (transverse polynomial shear-kernel obstruction).**
> Let \(K\)
> have characteristic zero and assume the scalar parent and singular-pencil
> hypotheses
>
> \[
> \det\operatorname{Hess}\Phi=c\in K^\times,
> \qquad
> \det\operatorname{Hess}_x(B+sA)=0.                     \tag{0.1}
> \]
>
> Suppose the normalized kernel line is independent of \(s\) and has, in
> linear coordinates \((x,y,z,w)\), the generator
>
> \[
> v=(P(z,w),1,0,0)^{\mathsf T},
> \qquad P\in K[z,w]\setminus K.                         \tag{0.2}
> \]
>
> Then a linear change of the transverse coordinates makes \(P\) a
> polynomial in \(z\) alone. There are polynomials
> \(p,q,h,k\in K[z]\), constants
> \(\alpha,\eta,\gamma\in K^\times\), and
> \(\xi,\beta\in K\) such that
>
> \[
> \begin{aligned}
> A&=\alpha y+p(z)w+q(z),\\
> B&=x(\eta z+\xi)+y\bigl(\beta-\eta J(z)\bigr)
>    +\frac{\gamma}{2}w^2+h(z)w+k(z),\\
> J'(z)&=P(z).
> \end{aligned}                                           \tag{0.3}
> \]
>
> Every scalar Schur descendant
>
> \[
> \psi_{\kappa,\mu}=B+\frac{\kappa}{2}A^2+\mu A,
> \qquad \kappa\ne0,                                     \tag{0.4}
> \]
>
> has
>
> \[
> \det\operatorname{Hess}\psi_{\kappa,\mu}
> =-\kappa\alpha^2\eta^2\gamma                          \tag{0.5}
> \]
>
> and an explicit triangular polynomial inverse. Thus this entire
> degree-unbounded nonlinear kernel stratum has no collision.

The exact verifier is
[`scripts/verify_hc4_univariate_shear_kernel_pencils.py`](scripts/verify_hc4_univariate_shear_kernel_pencils.py),
and its generated ledger is
[`artifacts/generated-results/hc4_univariate_shear_kernel_pencils.json`](artifacts/generated-results/hc4_univariate_shear_kernel_pencils.json).

## 1. Complete transverse integration

Let \(u=(z,w)\) and put

\[
\delta=P(u)\partial_x+\partial_y,
\qquad r=x-P(u)y.
\]

Then \(\ker\delta=K[u,r]\). If
\(\operatorname{Hess}(f)v=0\), every gradient component of \(f\) lies in
this invariant ring. Write the \(x\)-gradient component as \(X(r,u)\).
For \(u_i\in\{z,w\}\), mixed-partial equality gives

\[
(f_{u_i})_r=X_{u_i}-yP_{u_i}X_r.                         \tag{1.1}
\]

Since \(P\) is nonconstant, one of its derivatives is nonzero, and
invariance of the left side forces \(X_r=0\). The remaining mixed partials
integrate exactly to

\[
\boxed{
f=x a(z,w)+y b(z,w)+G(z,w),
\qquad db=-P\,da.}                                      \tag{1.2}
\]

In particular, \(dP\wedge da=0\). Conversely, (1.2) directly satisfies the
kernel equation. If \(C_f\) denotes the transverse \(2\)-by-\(2\) Hessian
block of \(f\), then

\[
\operatorname{Hess}f=
\begin{pmatrix}
0&0&(da)^{\mathsf T}\\
0&0&-P(da)^{\mathsf T}\\
da&-Pda&C_f
\end{pmatrix},                                          \tag{1.3}
\]

and its rank-one adjugate scalar is

\[
q_f=-(da)^{\mathsf T}\operatorname{adj}(C_f)da.          \tag{1.4}
\]

## 2. The curvature reduction

Apply (1.2) to \(A\). Its kernel-gradient pairing is

\[
D=P a_A+b_A.
\]

The relation \(db_A=-P\,da_A\) gives

\[
dD=a_A\,dP.                                              \tag{2.1}
\]

The bordered unit makes \(D\) a nonzero constant. Because \(dP\ne0\),
(2.1) forces

\[
a_A=0,\qquad b_A=\alpha\in K^\times,
\qquad A=\alpha y+G_A(z,w).                              \tag{2.2}
\]

For \(B\), write its nonconstant coefficient as \(a=a_B\). Generic rank
three gives \(da\ne0\), while \(dP\wedge da=0\). The standard polynomial
common-composite lemma in two variables supplies a polynomial \(H(z,w)\)
and univariate polynomials \(P_0,a_0,b_0\) such that

\[
P=P_0(H),\qquad a=a_0(H),\qquad b=b_0(H).                 \tag{2.3}
\]

Put \(U=dH\). The transverse Hessian block of \(B+sA\) has the form

\[
C_s=\rho\operatorname{Hess}H+\sigma UU^{\mathsf T}
    +\operatorname{Hess}(G_B+sG_A),                      \tag{2.4}
\]

where \(\rho=x a_0'(H)+y b_0'(H)\). For a \(2\)-by-\(2\) matrix the
adjugate is linear, and
\(U^{\mathsf T}\operatorname{adj}(UU^{\mathsf T})U=0\).
Consequently the coefficient of \(x\) in the adjugate scalar (1.4) is a
nonzero polynomial multiple of

\[
Q_H=H_z^2H_{ww}-2H_zH_wH_{zw}+H_w^2H_{zz}.              \tag{2.5}
\]

The bordered unit therefore forces \(Q_H=0\).

The straight-level lemma says that a nonconstant polynomial in two
variables satisfying (2.5) is a polynomial in one linear form. Indeed, on
the open set \(H_w\ne0\), differentiation along the level-curve tangent
\(H_w\partial_z-H_z\partial_w\) gives

\[
\left(H_w\partial_z-H_z\partial_w\right)
\left(\frac{H_z}{H_w}\right)=\frac{Q_H}{H_w^2}.          \tag{2.6}
\]

Thus every smooth component of a generic fiber has constant tangent and is
an affine line. Distinct generic fibers cannot contain nonparallel lines,
so the generic fibers are unions of parallel lines. Hence
\(H=H_0(\ell)\) for a linear form \(\ell\). After a linear transverse
coordinate change, \(P,a,b\) depend only on \(z\). This argument may be
performed after algebraic closure; the unique parallel direction is Galois
invariant and therefore descends to \(K\).

## 3. The univariate bordered unit

Now write the complete potential as

\[
f=x a(z)+y b(z)+G(z,w),\qquad b'=-P a'.                  \tag{3.1}
\]

Its adjugate scalar is \(-a'(z)^2G_{ww}(z,w)\). Applying this to \(A\) and
\(B\), the bordered determinant factors in \(K[s,z,w]\) as

\[
(P a_A+b_A)^2
(a_B'+sa_A')^2
(G_{B,ww}+sG_{A,ww})=c.                                  \tag{3.2}
\]

Every factor is a unit. Equation (2.2) already gives \(a_A=0\) and
\(b_A=\alpha\). The other factors give

\[
a_B'=\eta\in K^\times,\qquad
G_{A,ww}=0,\qquad G_{B,ww}=\gamma\in K^\times.           \tag{3.3}
\]

Integration of (3.3) and \(b_B'=-Pa_B'\) is exactly the normal form (0.3),
with

\[
c=\alpha^2\eta^2\gamma.                                 \tag{3.4}
\]

## 4. Explicit inverse

Put \(S=\kappa A+\mu\), and let
\((F_x,F_y,F_z,F_w)=\nabla\psi_{\kappa,\mu}\). From (0.3),

\[
\begin{aligned}
F_x&=\eta z+\xi,\\
F_y&=\beta-\eta J(z)+\alpha S,\\
F_w&=\gamma w+h(z)+S p(z),\\
F_z&=\eta x-\eta P(z)y+h'(z)w+k'(z)
     +S\bigl(p'(z)w+q'(z)\bigr).
\end{aligned}                                             \tag{4.1}
\]

These recover, in order,

\[
z,\qquad S,\qquad w,\qquad
A=\frac{S-\mu}{\kappa},\qquad
y=\frac{A-pw-q}{\alpha},\qquad x.                       \tag{4.2}
\]

All divisions are by constant units. Thus (4.2) is a polynomial inverse,
proving the theorem.

## 5. Reproduction and frontier

Run:

```bash
.venv/bin/python scripts/verify_hc4_univariate_shear_kernel_pencils.py
# cleanup only: verify committed inputs and historical/current frontier
.venv/bin/python scripts/verify_hc4_univariate_shear_kernel_pencils.py --audit-existing-only
```

The command checks the complete transverse kernel residual, the curvature
coefficient and tangent-slope identities, the univariate adjugate factor,
the bordered and descendant determinants, and every step of the triangular
inverse.

The cleanup-only mode hash-checks the committed ledger and imported equation
helper without importing SymPy, replaying the identities, or rewriting the
artifact. It treats the ledger's broad quasi-translation frontier as
historical: `HC4RSD5` later closes the fixed two-component subcase.

The result contains kernels of arbitrarily high degree and allows arbitrary
dependence on both transverse variables before the unit gate. `HC4RSD5`, in
[`HC4_TWO_COMPONENT_QUASITRANSLATION_KERNELS.md`](HC4_TWO_COMPONENT_QUASITRANSLATION_KERNELS.md),
subsequently extends the argument to every fixed primitive two-component
kernel in a constant support plane; Piola makes its coefficients transverse
before the shear reduction. Fixed kernels with three or four nonlinear
components, and parameter-moving nonlinear kernels, remain open.

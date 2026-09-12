# Complete all-degree closure of the moving rank-two `[2,2]` HC4 branch

## Status and scope

**Current status.**  `HC4RSD60` is the complete all-degree closure of the
rank-two square-zero `[2,2]` stratum.  The former length-three/length-four
handoff in Section 7 is historical: `HC4RSD61--63` close `[3,1]`, and
`HC4MR1` consolidates the reductions; HC4MRA1 and
[HC4MRA2](HC4_NEGATIVE_MOTION_POLYNOMIAL_OBSTRUCTION.md) supply the corrected
maximal-motion closure.

This note combines `HC4RSD56`, `HC4RSD58`, `HC4RSD59`, and the complete
three-variable quasi-translation classification.

> **Theorem HC4RSD60 — complete moving `[2,2]` closure.**
> Let
>
> \[
> S=\operatorname{Hess}\psi,\qquad
> T=\operatorname{Hess}A,\qquad
> N=S^{-1}T,
> \]
>
> over a characteristic-zero field, with
>
> \[
> \det S\in K^*,\qquad
> \operatorname{rank}T=2,\qquad
> N^2=0.
> \]
>
> Then the polynomial potential `psi` is either in a fixed-kernel branch
> already reduced by `HC4RSD42`, or after a constant linear change it has the
> cotangent form
>
> \[
> \psi=yP(x,w)+zQ(x,w)+R(x,w),
> \tag{0.1}
> \]
>
> with
>
> \[
> \det J_{x,w}(P,Q)\in K^*.
> \tag{0.2}
> \]
>
> Consequently the entire rank-two square-zero relative-nilpotent branch is
> exactly `HC2`/triangular or the `JC2` cotangent endpoint.  There is no
> genuinely moving `[2,2]` obstruction specific to `HC4`, in any degree.

The point of the proof is that the apparent parameter dependence in the
three-variable classification does not create a new four-variable geometry.

## 1. Canonical quasi-translation

`HC4RSD56` supplies a constant projective image apex.  Normalize it to

\[
p=e_w.
\]

Then `HC4RSD58` gives

\[
k=S^{-1}p=(\bar k,0),
\]

where

\[
H\bar k=0,
\qquad
H=\operatorname{Hess}_{x,y,z}\psi,
\]

and

\[
D=\bar k\cdot\nabla_{x,y,z}
\]

is a quasi-translation derivation satisfying

\[
D\bar k=0,
\qquad
D(\nabla\psi)=e_w,
\qquad
D(\psi_w)=1.
\tag{1.1}
\]

In particular the components of `bar k` generate the unit ideal, because if

\[
m=\nabla_{x,y,z}\psi_w,
\]

then

\[
m^{\mathsf T}\bar k=1.
\tag{1.2}
\]

## 2. Three-variable quasi-translation normal form

Work temporarily over the coefficient field `L=\overline{K(w)}`.  Corollary
5.8 of de Bondt's *Quasi-translations and singular Hessians* gives, after a
linear change in the ternary variables,

\[
\bar k=(0,b(u)g,a(u)g),
\qquad
 g\in L[u,a(u)v-b(u)t].
\tag{2.1}
\]

Equation (1.2) says the components of `bar k` are unimodular.  Therefore their
common factor `g` is a unit of the polynomial ring, and it may be absorbed into
`a,b`.  Thus

\[
\bar k=(0,b(u),a(u)).
\tag{2.2}
\]

There are now only two cases: the projective ratio `b:a` is constant in `u`,
or it genuinely moves in `u`.

## 3. Constant ratio: HC4RSD59

If `b:a` is constant in `u`, the ternary kernel line is constant over `L`.
The linear frame used to reach (2.2) is allowed to depend on `w`, but
`HC4RSD59` proves that a nonzero constant four-variable Hessian determinant
forces that projective kernel line to be constant in the original coordinates
as well.

Hence this case is already a fixed-kernel packet and reduces through
`HC4RSD42`.

## 4. Moving ratio forces passive linearity

Assume now that `b:a` genuinely varies with `u`.  On the chart `a\ne0`, put

\[
q(u,w)=\frac{b}{a},
\qquad
h=v-q(u,w)t.
\tag{4.1}
\]

The kernel direction is projectively

\[
r=(0,q,1),
\qquad
\partial_r=q\partial_v+\partial_t.
\]

The equation

\[
(\operatorname{Hess}_{u,v,t}\psi)r=0
\tag{4.2}
\]

says that `psi_u`, `psi_v`, and `psi_t` are invariant along `partial_r`.
Integrating the last two equations and using equality of mixed derivatives
gives

\[
\psi=\Phi(u,h,w)+A(u,w)t
\tag{4.3}
\]

up to an irrelevant affine term.

Differentiate (4.3) with respect to `u` at fixed `(v,t,w)`.  Since

\[
h_u=-q_ut,
\]

one gets

\[
\psi_u=\Phi_u-q_ut\Phi_h+A_ut.
\]

Applying `partial_r`, which is simply differentiation with respect to `t` at
fixed `(u,h,w)`, the first component of (4.2) becomes

\[
A_u-q_u\Phi_h=0.
\tag{4.4}
\]

Because `q_u\ne0`, equation (4.4) forces `Phi_h` to depend only on `(u,w)`.
Hence `Phi` is affine in `h`, and after expanding `h=v-qt`,

\[
\boxed{
\psi=vP(u,w)+tQ(u,w)+R(u,w).
}
\tag{4.5}
\]

Thus the genuinely moving ternary quasi-translation branch is already
cotangent-linear in its two passive coordinates.

## 5. Descent through a `w`-dependent frame

It remains to show that the possibly `w`-dependent ternary frame used above
does not turn (4.5) into a new geometry.

Let

\[
X=(x,y,z)^{\mathsf T},
\qquad
Y=T(w)^{-1}X=(u,v,t)^{\mathsf T},
\qquad
B=T^{-1},
\qquad
\Omega=T^{-1}T'.
\tag{5.1}
\]

For (4.5), put

\[
p=P_u,\qquad q_0=Q_u.
\]

The ternary Hessian is

\[
H_Y=
\begin{pmatrix}
* & p & q_0\\
p&0&0\\
q_0&0&0
\end{pmatrix},
\]

with kernel vector

\[
r=(0,q_0,-p)^{\mathsf T}.
\tag{5.2}
\]

Its adjugate is exactly

\[
\operatorname{adj}(H_Y)=-rr^{\mathsf T}.
\tag{5.3}
\]

As in `HC4RSD59`, the pulled-back mixed Hessian column is

\[
n=g_w-H_Y\Omega Y-\Omega^{\mathsf T}g,
\qquad
 g=(vP_u+tQ_u+R_u,P,Q)^{\mathsf T}.
\tag{5.4}
\]

Therefore

\[
\det\operatorname{Hess}_{X,w}\psi
=\det(B)^2(r^{\mathsf T}n)^2.
\tag{5.5}
\]

Since `r^T H_Y=0`,

\[
r^{\mathsf T}n
=r^{\mathsf T}g_w-(\Omega r)^{\mathsf T}g.
\tag{5.6}
\]

Write `s=Omega r`.  The coefficients of the passive variables `v,t` in
(5.6) are respectively

\[
-s_1P_u,
\qquad
-s_1Q_u.
\]

The determinant in (5.5) is a nonzero constant, so `r^T n` is a unit in
`L[u,v,t]`.  Since `(P_u,Q_u)\ne(0,0)`,

\[
s_1=0.
\tag{5.7}
\]

In terms of the first row of `Omega`,

\[
\omega_{12}Q_u-\omega_{13}P_u=0.
\tag{5.8}
\]

But the exceptional case is precisely the case in which the ratio
`Q_u:P_u` is not constant in `u`.  Thus

\[
\boxed{\omega_{12}=\omega_{13}=0.}
\tag{5.9}
\]

Now the first row of `B` evolves only by scalar multiplication:

\[
(e_1^{\mathsf T}B)'
=-e_1^{\mathsf T}\Omega B
=-\omega_{11}e_1^{\mathsf T}B.
\tag{5.10}
\]

Hence the projective active covector line is constant in the original ternary
coordinates.

Choose a constant linear coordinate `x` spanning that covector line.  Then
`u=alpha(w)x`.  Expanding the two passive moving linear forms `v,t` in any
constant complement `(y,z)` and collecting their `x`-parts into `R`, equation
(4.5) becomes, in **constant** coordinates,

\[
\psi=y\widetilde P(x,w)+z\widetilde Q(x,w)+\widetilde R(x,w).
\tag{5.11}
\]

Because `psi` is polynomial, its coefficients

\[
\widetilde P=\psi_y,
\qquad
\widetilde Q=\psi_z
\]

are polynomials as well; no rational descent problem remains.

## 6. Exact cotangent endpoint

For the constant-coordinate form (5.11), direct differentiation gives

\[
\det\operatorname{Hess}\psi
=\left(
\widetilde P_x\widetilde Q_w
-\widetilde P_w\widetilde Q_x
\right)^2.
\tag{6.1}
\]

Since the left side is a nonzero constant, the plane Jacobian

\[
J_{x,w}(\widetilde P,\widetilde Q)
\]

is itself a nonzero constant.  Thus (5.11) is exactly the cotangent lift of a
plane Keller map.

This proves `HC4RSD60`.

## 7. Consequence

The complete rank-two square-zero relative-nilpotent stratum is now closed in
all degrees:

\[
\boxed{
[2,2]\text{ moving frame}
\;\Longrightarrow\;
\text{fixed/triangular/HC2 or exact JC2 cotangent lift}.
}
\]

At this historical stage, moving-frame work only needed to consider the
nilpotency types with a length-three or length-four chain, not `[2,2]`.
The length-three continuation is closed by `HC4RSD61--63`; the final
length-four motion signs use the corrected HC4MRA1/HC4MRA2 proof.

## 8. External input

The only external normal-form theorem is Corollary 5.8 of
M. de Bondt, *Quasi-translations and singular Hessians*, arXiv:1501.05168.
It classifies every three-variable quasi-translation, after a linear change,
as

\[
(0,bg,ag),\qquad a,b\in K[u],\quad g\in K[u,av-bt].
\]

Everything after (2.1), including the unimodular collapse of `g`, the Hessian
integrability dichotomy, and the moving-frame descent, is specific to the HC4
Hessian pencil.

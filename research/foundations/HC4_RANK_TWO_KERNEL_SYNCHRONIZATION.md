# Parameter-moving kernel synchronization in the rank-two `[2,2]` HC4 branch

## Status and scope

**Current status.**  This is an intermediate proof-map note.  The exceptional
quasi-translation form left in Section 4 is closed by `HC4RSD60`, and the
complete square-zero branch is subsumed by `HC4MR1`.

This note continues `HC4RSD58`.  Work over a characteristic-zero field and write

\[
S=\operatorname{Hess}\psi,\qquad
T=\operatorname{Hess}A,\qquad
N=S^{-1}T,
\]

with

\[
\det S\in K^*,\qquad \operatorname{rank}T=2,\qquad N^2=0.
\]

After the constant projective image apex supplied in `HC4RSD56` is normalized to
`e_w`, `HC4RSD58` gives the canonical polynomial kernel section

\[
k=S^{-1}e_w=(\bar k,0),\qquad
H\bar k=0,
\]

where

\[
H=\operatorname{Hess}_{x,y,z}\psi,
\]

and the quasi-translation derivation

\[
D=\bar k\cdot\nabla_{x,y,z}
\]

satisfies

\[
D\bar k=0,\qquad D(\nabla\psi)=e_w,\qquad D(\nabla A)=0.
\]

The remaining issue was that the ternary singular-Hessian classification may be
performed only over `K(w)`, so its normalizing linear frame is allowed to depend
on `w`.

> **Theorem HC4RSD59 — parameter-moving kernel synchronization.**
> Suppose that over `K(w)` the generic rank-two ternary Hessian `H` has a
> constant projective kernel line.  Then that projective kernel line is already
> constant in the original `(x,y,z)` coordinates.  Hence this branch reduces to
> the fixed-kernel classification of `HC4RSD42` and therefore to a triangular
> packet, `HC2`, or the exact `JC2` cotangent endpoint.

The theorem is degree-free and, crucially, allows the `K(w)`-normalizing frame
to move rationally with `w`.

## 1. A moving frame

Let

\[
X=(x,y,z)^{\mathsf T}.
\]

Over `K(w)` choose `T(w)\in GL_3(K(w))` so that in the moving coordinates

\[
Y=T(w)^{-1}X=(u,v,t)^{\mathsf T}
\]

the ternary kernel is the line spanned by `e_t`.  Since
`\operatorname{rank}H=2`, Hessian integrability in the ternary variables gives

\[
\psi(X,w)=\Phi(u,v,w)+a(w)t
\tag{1.1}
\]

up to terms affine in `X`, which do not affect the Hessian.

Put

\[
B=T^{-1},\qquad \Omega=T^{-1}T'.
\tag{1.2}
\]

At fixed `X`,

\[
Y_w=-\Omega Y.
\tag{1.3}
\]

Write

\[
g=\nabla_Y\psi=(\Phi_u,\Phi_v,a)^{\mathsf T},
\]

and

\[
H_Y=\operatorname{Hess}_Y\psi=
\begin{pmatrix}
\Phi_{uu}&\Phi_{uv}&0\\
\Phi_{uv}&\Phi_{vv}&0\\
0&0&0
\end{pmatrix}.
\tag{1.4}
\]

The ternary Hessian in the original variables is

\[
H_X=B^{\mathsf T}H_YB.
\tag{1.5}
\]

## 2. The only relevant mixed derivative

The `Xw` Hessian column, pulled back to the moving frame, is

\[
n
=g_w-H_Y\Omega Y-\Omega^{\mathsf T}g.
\tag{2.1}
\]

Because the third row of `H_Y` vanishes, only the third component of `n`
contributes to the determinant of the full four-variable Hessian.  If

\[
\Omega e_t=(\omega_{13},\omega_{23},\omega_{33})^{\mathsf T},
\]

then

\[
n_3
=a'-\omega_{13}\Phi_u-\omega_{23}\Phi_v-\omega_{33}a.
\tag{2.2}
\]

Furthermore

\[
\operatorname{adj}(H_Y)
=\det\operatorname{Hess}_{u,v}\Phi\;e_te_t^{\mathsf T}.
\tag{2.3}
\]

The standard singular-block determinant identity therefore gives

\[
\det\operatorname{Hess}_{X,w}\psi
=-\det(B)^2
\det\operatorname{Hess}_{u,v}\Phi\; n_3^2.
\tag{2.4}
\]

This is the moving-frame version of the explicit shear identity previously
observed for

\[
u=x-q_1(w)z,\qquad v=y-q_2(w)z.
\]

## 3. Constant determinant kills projective motion

The left side of (2.4) is a nonzero element of `K`.  Work in the polynomial
ring

\[
K(w)[u,v].
\]

The scalar `det(B)` is a unit there.  Hence both

\[
\det\operatorname{Hess}_{u,v}\Phi
\quad\text{and}\quad
n_3
\]

must themselves be units in `K(w)`.

Differentiate (2.2) with respect to `u,v`.  Since `n_3` is independent of
`u,v`,

\[
\operatorname{Hess}_{u,v}\Phi
\begin{pmatrix}
\omega_{13}\\\omega_{23}
\end{pmatrix}
=0.
\]

But the binary Hessian determinant is a nonzero unit.  Thus

\[
\boxed{\omega_{13}=\omega_{23}=0.}
\tag{3.1}
\]

Consequently

\[
\Omega e_t=\omega_{33}e_t.
\]

Since `T'=T\Omega`,

\[
(Te_t)'=T\Omega e_t=\omega_{33}Te_t.
\tag{3.2}
\]

Therefore the projective line `K(w)Te_t` is independent of `w`.
In the original coordinates the ternary kernel direction is constant.

This proves `HC4RSD59`.

## 4. Historical remainder (closed by `HC4RSD60`)

The only rank-two `[2,2]` branch not covered by `HC4RSD59` is the genuinely
nonconstant ternary quasi-translation normal form.  De Bondt's complete
three-variable quasi-translation classification gives, over `K(w)`, after a
linear change,

\[
\bar k=(0,b(u)g, a(u)g),
\qquad
 g=g(u,a(u)v-b(u)t).
\tag{4.1}
\]

The normalization `m^{\mathsf T}\bar k=1` from `HC4RSD58` makes the components
of `\bar k` unimodular, so the common factor `g` is a unit over `K(w)`.  Thus
only

\[
\bar k=(0,b(u),a(u))
\tag{4.2}
\]

remains.  If `a/b` is constant in `u`, this is already `HC4RSD59`.  If the
ratio genuinely varies, Hessian integrability forces the potential to be
linear in the two passive variables; this is the exceptional branch attacked
in the next note.

## 5. External input

The only external classification used in the final reduction statement is
Corollary 5.8 of Michiel de Bondt, *Quasi-translations and singular Hessians*,
arXiv:1501.05168: every three-variable quasi-translation is linearly
conjugate to

\[
(0,bg,ag),
\qquad
 a,b\in K[u],\quad g\in K[u,av-bt].
\]

The moving-frame determinant argument of Sections 1--3 is independent of that
classification.

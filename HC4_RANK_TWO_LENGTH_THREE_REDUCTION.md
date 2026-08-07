# All-degree reduction of the rank-two length-three `[3,1]` HC4 branch

## Status and scope

This note attacks the remaining rank-two nilpotent type after the complete
`[2,2]` closure `HC4RSD60`.

Let

\[
S=\operatorname{Hess}\psi,\qquad
T=\operatorname{Hess}A,\qquad
N=S^{-1}T,
\]

with

\[
\det S=\delta\in K^*,\qquad
\operatorname{rank}T=2,\qquad
N^3=0,\qquad N^2\ne0.
\tag{0.1}
\]

> **Theorem HC4RSD61 — fixed active-plane reduction for `[3,1]`.**
> After a constant linear change of coordinates and removal of affine terms,
>
> \[
> \boxed{A=A(x,w).}
> \tag{0.2}
> \]
>
> In particular
>
> \[
> T=
> \begin{pmatrix}
> K&0\\0&0
> \end{pmatrix},
> \qquad
> K=\operatorname{Hess}_{x,w}A,
> \qquad
> \det K\ne0
> \tag{0.3}
> \]
>
> in the active/passive splitting `(x,w)|(y,z)`.
>
> If
>
> \[
> S=
> \begin{pmatrix}
> E&B\\B^{\mathsf T}&C
> \end{pmatrix}
> \tag{0.4}
> \]
>
> in the same splitting, then the passive Hessian block satisfies
>
> \[
> \det C=0.
> \tag{0.5}
> \]
>
> The square-zero case is exactly `C=0`; hence in the genuine `[3,1]` branch
> `C` has generic rank one.  If `v` spans `ker C` over the fraction field and
> `r=Bv`, then
>
> \[
> \boxed{r^{\mathsf T}\operatorname{adj}(K)r=0.}
> \tag{0.6}
> \]
>
> Thus the complete four-variable moving `[3,1]` problem reduces to a binary
> Hessian metric `K`, a rank-one passive Hessian, and one null coupling vector.

This is an all-degree reduction; no homogeneous or leading-degree assumption
is used.

## 1. The constant apex forces a two-dimensional ternary kernel

The gradient map `nabla A` has Jacobian rank two.  By the same small-rank
Hessian theorem used in `HC4RSD56`, its image has a nonzero constant projective
apex `p`.  At a generic smooth point,

\[
p\in\operatorname{im}T.
\tag{1.1}
\]

Make a constant change so that

\[
p=e_w.
\tag{1.2}
\]

Since `T` is symmetric,

\[
\ker T=(\operatorname{im}T)^\perp.
\]

Therefore every `q in ker T` satisfies

\[
q_w=p^{\mathsf T}q=0.
\tag{1.3}
\]

The kernel of `T` has dimension two, so it is a two-plane contained entirely
in the ternary `(x,y,z)` directions.  Hence

\[
\operatorname{rank}\operatorname{Hess}_{x,y,z}A\le1.
\tag{1.4}
\]

The generic rank is exactly one: if it were zero, the three ternary columns of
`T` would have only their `w` entries available and the full symmetric matrix
could not have rank two with a two-dimensional ternary kernel except in the
already affine/triangular degeneration.

## 2. Rank-one ternary Hessian and synchronization

Work temporarily over `L=\overline{K(w)}`.  The rank-one Hessian
classification gives, after a linear ternary change,

\[
A=F(u,w)+a(w)^{\mathsf T}X+b(w),
\qquad
u=c(w)^{\mathsf T}X,
\tag{2.1}
\]

where `X=(x,y,z)^T`, `c(w)\ne0`, and `F_{uu}\ne0`.  Equivalently, directly in
the original ternary coordinates,

\[
H_3:=\operatorname{Hess}_X A
=F_{uu}\,c c^{\mathsf T}.
\tag{2.2}
\]

Let

\[
m=\nabla_X A_w.
\]

Because `ker T` is the full two-plane `ker H_3`, every vector annihilated by
`H_3` must also be annihilated by `m^T`.  Hence

\[
m\in\operatorname{im}H_3=Lc.
\tag{2.3}
\]

Differentiate

\[
\nabla_XA=F_u c+a
\]

with respect to `w` at fixed `X`.  Since `u_w=c'{}^{\mathsf T}X`,

\[
m
=(F_{uw}+F_{uu}c'{}^{\mathsf T}X)c+F_uc'+a'.
\tag{2.4}
\]

Reduce (2.4) modulo the line `Lc`.  Equation (2.3) gives

\[
F_u[c']+[a']=0.
\tag{2.5}
\]

But `F_{uu}\ne0`, so `F_u` is nonconstant in `u`.  Therefore

\[
[c']=0,\qquad[a']=0
\]

in the quotient `L^3/Lc`.  Thus

\[
c'\parallel c.
\tag{2.6}
\]

The projective active covector line `Lc(w)` is constant.  Returning to the
base field and choosing a constant coordinate `x` on that line, all nonlinear
dependence of `A` is in `(x,w)`; transverse affine terms are constant and may
be removed.  This proves (0.2).

## 3. The passive Hessian is singular

Use the active/passive splitting `(x,w)|(y,z)` and write

\[
S=\begin{pmatrix}E&B\\B^{\mathsf T}&C\end{pmatrix},
\qquad
T=\begin{pmatrix}K&0\\0&0\end{pmatrix}.
\tag{3.1}
\]

The coefficient of `s^2` in

\[
\det(S+sT)=\delta
\]

is the complementary-minor product

\[
[s^2]\det(S+sT)=\det K\,\det C.
\tag{3.2}
\]

Since `det K\ne0`, constancy of the pencil determinant gives

\[
\det C=0.
\]

If `C=0`, then an invertible matrix `S` of the form

\[
\begin{pmatrix}E&B\\B^{\mathsf T}&0\end{pmatrix}
\]

has inverse with zero active-active block.  Consequently

\[
T S^{-1}T=0,
\]

which is equivalent to `N^2=0`.  Therefore the genuine `[3,1]` branch has

\[
\operatorname{rank}C=1.
\tag{3.3}
\]

## 4. The remaining coefficient is one binary null equation

Work over the fraction field and choose a passive basis in which

\[
C=\begin{pmatrix}c&0\\0&0\end{pmatrix},
\qquad c\ne0.
\tag{4.1}
\]

Write the coupling matrix as

\[
B=\begin{pmatrix}b_{11}&b_{12}\\b_{21}&b_{22}\end{pmatrix}
\]

and put

\[
r=\begin{pmatrix}b_{12}\\b_{22}\end{pmatrix}=B e_2.
\tag{4.2}
\]

A direct `4x4` determinant expansion gives

\[
[s]\det(S+sT)
=-c\,r^{\mathsf T}\operatorname{adj}(K)r.
\tag{4.3}
\]

The left side vanishes.  Since `c\ne0`,

\[
r^{\mathsf T}\operatorname{adj}(K)r=0.
\]

In invariant form, if `v` spans `ker C`, then `r=Bv`, proving (0.6).
Notice that `r\ne0`: otherwise the passive kernel vector `(0,v)` would also
annihilate the full metric `S`, contradicting `det S\ne0`.

## 5. Passive-kernel coordinates

For later use, over the active fraction field write a rank-one passive Hessian
in the normal form

\[
\psi=\Phi(h,x,w)+a(x,w)t,
\qquad
h=y-q(x,w)t.
\tag{5.1}
\]

The passive kernel is `(q,1)`, and its coupling into the active directions is
exactly

\[
r=\nabla_{x,w}a-\Phi_h\nabla_{x,w}q.
\tag{5.2}
\]

Substituting (5.2) into (0.6) yields a quadratic polynomial identity in the
nonconstant parameter `Phi_h`.  Since `Phi_{hh}\ne0` in the rank-one branch,
all three coefficients vanish:

\[
\begin{aligned}
(\nabla q)^{\mathsf T}\operatorname{adj}(K)\nabla q&=0,\\
(\nabla a)^{\mathsf T}\operatorname{adj}(K)\nabla q&=0,\\
(\nabla a)^{\mathsf T}\operatorname{adj}(K)\nabla a&=0.
\end{aligned}
\tag{5.3}
\]

Over the two-dimensional active fraction field the null cone of the
nondegenerate binary form `adj(K)` consists of two lines.  Thus `nabla a` and
`nabla q` lie on the same null characteristic whenever both are nonzero.

Equation (5.3) is the residual `[3,1]` obstruction.  It is now a **purely
binary characteristic problem**; the four-variable moving flag has been
eliminated.

## 6. Next target

The next theorem should classify polynomial/rational solutions of

\[
(\nabla q)^{\mathsf T}\operatorname{adj}(\operatorname{Hess}A)\nabla q=0
\tag{6.1}
\]

for a bivariate polynomial `A` with nonzero Hessian determinant.  Geometrically,
(6.1) says that the level curves of `q` are asymptotic curves of the Hessian
metric of `A`.  If every polynomial first integral of such a null foliation is
a function of a polynomial characteristic coordinate, then the remaining
`[3,1]` packet is triangular/`HC2`.

## 7. External input

As in `HC4RSD56`, the constant projective apex for a rank-two polynomial
Hessian is supplied by M. de Bondt, *Polynomial Hessians with small rank*,
arXiv:1609.03904v2.  The rank-one Hessian reduction in Section 2 is the
rank-one case of the same small-rank theory (equivalently the elementary Hesse
classification in rank one).

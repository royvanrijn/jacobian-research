# Rank-two quadratic pivots for HC4

## Status

HC4RSD8 reduces every singular-pencil quadratic scalar pivot to Hessian
rank at most two and closes the passive-rank-zero part of rank two through
HC4RSD5. This note closes the remaining passive-rank-one part.

> **Theorem HC4RSD9 (rank-two quadratic-pivot classification).** Let \(K\)
> have characteristic zero, let \(x=(x_1,\ldots,x_4)\), and put
>
> \[
> \Phi(t,x)=\frac{\lambda}{2}t^2+tA(x)+B(x).
> \tag{0.1}
> \]
>
> Suppose
>
> \[
> \deg A\leq2,\qquad
> \operatorname{rank}\operatorname{Hess}A=2,
> \tag{0.2}
> \]
>
> \[
> \det\operatorname{Hess}_{t,x}\Phi=c\in K^\times,
> \qquad
> \det\operatorname{Hess}_x(B+sA)=0.
> \tag{0.3}
> \]
>
> Then every Schur descendant
>
> \[
> \psi_{\kappa,\mu}
> =B+\frac{\kappa}{2}A^2+\mu A,
> \qquad \kappa\ne0,
> \tag{0.4}
> \]
>
> has a polynomially invertible gradient. In particular, it has no
> nontrivial gradient collision.
>
> More precisely, after scalar extension and constant affine coordinate
> changes, either the pencil kernel lies in a fixed two-plane and HC4RSD5
> applies, or, up to an affine summand in the potential,
>
> \[
> \boxed{
> \begin{aligned}
> A&=xy+w,\\
> B&=xz+\frac{\rho}{2}
>       \bigl(y+h(x)A\bigr)^2
>       +\beta(x)y+\gamma(x)A+\delta(x),
> \end{aligned}}
> \tag{0.5}
> \]
>
> where \(\rho\in K^\times\) and
> \(h,\beta,\gamma,\delta\in K[x]\).

The exact block identities and normal-form calibrations are replayed by
[scripts/verify_hc4_quadratic_rank_two_pivots.py](scripts/verify_hc4_quadratic_rank_two_pivots.py),
which writes
[artifacts/generated-results/hc4_quadratic_rank_two_pivots.json](artifacts/generated-results/hc4_quadratic_rank_two_pivots.json).

## 1. Hyperbolic pivot and passive splitting

By HC4RSD8, the linear part of \(A\) is nonzero on the kernel of its
constant Hessian. Over the algebraic closure, constant affine coordinates
therefore give

\[
A=xy+w.
\tag{1.1}
\]

Use \((x,y)\) as active and \((z,w)\) as passive coordinates. Write

\[
M(s)=\operatorname{Hess}(B+sA)
=
\begin{pmatrix}
K(s)&D\\
D^{\mathsf T}&E
\end{pmatrix}.
\tag{1.2}
\]

The coefficient of \(s^2\) in \(\det M\) is
\(-\det E\), so pencil singularity gives

\[
\det E=0.
\tag{1.3}
\]

The passive part of \(\nabla A\) is \((0,1)\). After using (1.3), the
coefficient of \(s^2\) in the full parent determinant is \(B_{zz}\).
Constancy therefore gives \(B_{zz}=0\), and (1.3) then gives
\(B_{zw}=0\). Hence

\[
B=z\,b(x,y)+C(x,y,w).
\tag{1.4}
\]

The case \(E=0\) is the fixed-support two-component kernel stratum closed
in HC4RSD8 by HC4RSD5. Assume from now on that

\[
C_{ww}\ne0.
\tag{1.5}
\]

## 2. The two determinant faces

Put

\[
p=b_x,\qquad q=b_y,\qquad e=C_{ww}.
\]

Exact block expansion gives

\[
[s]\det M=2epq.
\tag{2.1}
\]

Since \(e\ne0\) and the coefficient ring is a domain, either \(p=0\) or
\(q=0\). Interchanging \(x\) and \(y\), take \(q=0\), so \(b=b(x)\).
If \(b'\) vanished, the \(z\)-row and column of both the pencil and the
parent Hessian would be zero, contradicting (0.3). Thus \(b'\ne0\).

The constant pencil coefficient and the full parent determinant are

\[
\det M=(b')^2
\bigl(C_{yw}^2-C_{yy}C_{ww}\bigr),
\tag{2.2}
\]

\[
\det\operatorname{Hess}\Phi
=(b')^2
\bigl(C_{yy}-2xC_{yw}+x^2C_{ww}\bigr)
+\lambda\det M.
\tag{2.3}
\]

Equations (0.3), (2.2), and (2.3) yield

\[
C_{yw}^2-C_{yy}C_{ww}=0,
\tag{2.4}
\]

\[
(b')^2
\bigl(C_{yy}-2xC_{yw}+x^2C_{ww}\bigr)=c.
\tag{2.5}
\]

Both factors in (2.5) are polynomials, and their product is a unit.
Therefore \(b'\in K^\times\). Rescaling \(z\) and discarding an affine
summand normalizes

\[
b=x.
\tag{2.6}
\]

The directional curvature in parentheses in (2.5) is then a constant
\(\rho\in K^\times\).

## 3. Exact integration

Introduce the polynomial coordinate

\[
v=w+xy=A
\tag{3.1}
\]

and write

\[
\widetilde C(x,y,v)=C(x,y,v-xy).
\]

At fixed \(x\),

\[
\partial_y|_v=\partial_y-x\partial_w.
\]

Consequently, (2.4)--(2.5) become

\[
\det\operatorname{Hess}_{y,v}\widetilde C=0,
\qquad
\widetilde C_{yy}=\rho.
\tag{3.2}
\]

Integrating the second equation gives

\[
\widetilde C
=\frac{\rho}{2}y^2+y f(x,v)+g(x,v).
\tag{3.3}
\]

The first equation is

\[
\rho\bigl(yf_{vv}+g_{vv}\bigr)-f_v^2=0.
\tag{3.4}
\]

Its \(y\)-coefficient gives \(f_{vv}=0\), so

\[
f=\alpha(x)v+\beta(x).
\]

The remaining equation gives

\[
g=\frac{\alpha(x)^2}{2\rho}v^2+\gamma(x)v+\delta(x).
\]

Writing \(h=\alpha/\rho\) yields exactly (0.5).

## 4. Triangular inverse

Let \(F=\nabla\psi_{\kappa,\mu}\), and put

\[
Y=y+h(x)A.
\tag{4.1}
\]

The normal form gives the exact recovery equations

\[
\begin{aligned}
F_z&=x,\\
F_y-xF_w&=\rho Y+\beta(x),\\
F_w&=\rho h(x)Y+\gamma(x)+\kappa A+\mu.
\end{aligned}
\tag{4.2}
\]

Since \(\rho\kappa\ne0\), they recover successively

\[
\begin{aligned}
x&=F_z,\\
Y&=\rho^{-1}\bigl(F_y-xF_w-\beta(x)\bigr),\\
A&=\kappa^{-1}
\bigl(F_w-\rho h(x)Y-\gamma(x)-\mu\bigr),\\
y&=Y-h(x)A,\\
w&=A-xy.
\end{aligned}
\tag{4.3}
\]

Finally \(F_x=z+R(x,y,w)\) for an explicit polynomial \(R\), so \(z\) is
recovered polynomially as well. Thus \(\nabla\psi_{\kappa,\mu}\) is a
polynomial automorphism. Direct differentiation also gives

\[
\det\operatorname{Hess}\psi_{\kappa,\mu}=-\kappa\rho.
\tag{4.4}
\]

## 5. Revised quadratic frontier

HC4RSD7 closes pivot-Hessian rank zero, HC4RSD9 closes rank two, and
HC4RSD8 excludes ranks three and four. At the time of this reduction, the
only quadratic scalar pivots remaining in the singular-pencil reverse-Schur
programme had

\[
\operatorname{rank}\operatorname{Hess}A=1.
\]

Those pivots have normal form \(A=w+q(u)^2/2\) after constant affine
coordinates. Their passive three-variable Hessian block is singular.
HC4RSD10 subsequently completes that low-corank integrability classification
in
[`HC4_QUADRATIC_RANK_ONE_PIVOTS.md`](HC4_QUADRATIC_RANK_ONE_PIVOTS.md)
and proves that every descendant is triangular.

## 6. Reproduction

Run:

~~~bash
.venv/bin/python scripts/verify_hc4_quadratic_rank_two_pivots.py
~~~

The command checks all block determinant faces, the all-degree normal-form
calibration, the descendant Hessian determinant, and the triangular
recovery identities.

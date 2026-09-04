# Rank-one quadratic pivots for HC4

## Status

HC4RSD8 reduces the singular-pencil quadratic scalar branch to pivot-Hessian
ranks one and two, and HC4RSD9 closes rank two. This note closes rank one.

> **Theorem HC4RSD10 (rank-one quadratic-pivot classification).** Let \(K\)
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
> \operatorname{rank}\operatorname{Hess}A=1,
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
> has a polynomially invertible gradient. More precisely, after scalar
> extension, constant affine coordinate changes, rescaling, and removal of
> an affine summand in the potential,
>
> \[
> \boxed{
> \begin{aligned}
> A&=\frac12x^2+w,\\
> B&=xz+\frac{\rho}{2}\bigl(y+h(x)w\bigr)^2
>       +\alpha(x)y+\gamma(x)w+\delta(x),
> \end{aligned}}
> \tag{0.5}
> \]
>
> where \(\rho\in K^\times\) and
> \(h,\alpha,\gamma,\delta\in K[x]\).

Together with HC4RSD7--HC4RSD9, this closes every quadratic scalar pivot in
the identically singular reduced-pencil programme.

The exact block identities and normal-form calibration are replayed by
[scripts/verify_hc4_quadratic_rank_one_pivots.py](scripts/verify_hc4_quadratic_rank_one_pivots.py),
which writes
[artifacts/generated-results/hc4_quadratic_rank_one_pivots.json](artifacts/generated-results/hc4_quadratic_rank_one_pivots.json).
Its nonsingular-quadratic and higher-degree `open_frontier` is stage-local:
`HC4RSD11--16` settle the quadratic cancellation handoff, while `HC4MR1`
consolidates the auxiliary higher-degree reductions but leaves the
negative maximal-motion sign open after
[HC4MRA1](HC4_MOTION_FRAME_TRANSPORT_AUDIT.md).

## 1. The passive three-by-three Hessian

By HC4RSD8, the linear part of \(A\) is nonzero on the kernel of its
constant Hessian. Constant affine coordinates therefore give

\[
A=\frac12x^2+w.
\tag{1.1}
\]

Write \(u=(y,z,w)\), \(a=(0,0,1)^{\mathsf T}\), and

\[
\operatorname{Hess}(B+sA)
=
\begin{pmatrix}
k+s&d^{\mathsf T}\\
d&E
\end{pmatrix},
\qquad
E=\operatorname{Hess}_uB.
\tag{1.2}
\]

The coefficient of \(s\) in the pencil determinant is \(\det E\).
Therefore

\[
\det E=0.
\tag{1.3}
\]

After ordering the full parent Hessian by the active variable \(x\) and
then by \((t,u)\), the coefficient of \(s\) is

\[
\det
\begin{pmatrix}
\lambda&a^{\mathsf T}\\
a&E
\end{pmatrix}
=-a^{\mathsf T}\operatorname{adj}(E)a.
\tag{1.4}
\]

Thus

\[
a^{\mathsf T}\operatorname{adj}(E)a=0.
\tag{1.5}
\]

The bordered-unit corank gate forces the four-by-four pencil in (1.2) to
have generic rank three.

## 2. Passive ranks zero and two

If \(\operatorname{rank}E=0\), one active row and column can raise rank by
at most two, so the pencil rank is at most two. This contradicts the
corank gate.

Suppose \(\operatorname{rank}E=2\). Over the fraction field, write

\[
\operatorname{adj}(E)=qvv^{\mathsf T},
\qquad q\ne0.
\tag{2.1}
\]

Equation (1.5) gives \(v^{\mathsf T}a=0\). Pencil singularity also gives

\[
0=\det
\begin{pmatrix}
k+s&d^{\mathsf T}\\
d&E
\end{pmatrix}
=-d^{\mathsf T}\operatorname{adj}(E)d,
\]

and hence \(v^{\mathsf T}d=0\). Therefore

\[
\begin{pmatrix}0\\v\end{pmatrix}
\]

is a kernel vector of the reduced pencil and is orthogonal to
\(\nabla A=(x,a)\). It remains a kernel vector after adjoining the parent
border, contradicting (0.3). Passive rank two is impossible.

Consequently

\[
\operatorname{rank}E=1.
\tag{2.2}
\]

## 3. Rank-one Hessian integration and the unit frame

Apply the rank-one polynomial-Hessian normal form over \(K(x)\), as recorded
in the canonical proof of HC4T11. After primitive polynomial
factorization, write

\[
E=q\,\ell\ell^{\mathsf T},
\qquad
\ell=(p(x),r(x),h_0(x))^{\mathsf T}.
\tag{3.1}
\]

For an arbitrary passive coupling \(d\), direct determinant expansion gives

\[
\det\operatorname{Hess}\Phi
=q\,\det(a,d,\ell)^2.
\tag{3.2}
\]

The left side is a nonzero constant, so both factors on the right are
polynomial units. Hence \(q=\rho\in K^\times\). Integrating the now
passive-constant Hessian gives

\[
B=\frac{\rho}{2}L^2+b_1(x)y+b_2(x)z+b_3(x)w+\delta(x),
\tag{3.3}
\]

where

\[
L=p(x)y+r(x)z+h_0(x)w.
\tag{3.4}
\]

The \(x\)-passive Hessian column is

\[
d=\rho\bigl((\ell'\mathbin{\cdot}u)\ell+L\ell'\bigr)+b'(x).
\tag{3.5}
\]

Substitution into the determinant frame gives

\[
\det(a,d,\ell)
=\rho L\bigl(p'r-r'p\bigr)+b_1'r-b_2'p.
\tag{3.6}
\]

This polynomial is a unit. Its passive-linear coefficient must vanish:

\[
p'r-r'p=0.
\tag{3.7}
\]

In characteristic zero, (3.7) makes \((p,r)\) a fixed projective
direction. The remaining term in (3.6) is a unit, so their common polynomial
factor is a unit as well. A constant change of \((y,z)\) therefore gives

\[
(p,r)=(1,0),
\qquad
b_2'\in K^\times.
\tag{3.8}
\]

Rescale \(z\), discard its affine summand, and write \(h=h_0\). Then
\(b_2=x\), and (3.3) is exactly (0.5).

## 4. Triangular inverse

Let \(F=\nabla\psi_{\kappa,\mu}\), and put

\[
Y=y+h(x)w,
\qquad
A=\frac12x^2+w.
\tag{4.1}
\]

The normal form gives

\[
\begin{aligned}
F_z&=x,\\
F_y&=\rho Y+\alpha(x),\\
F_w&=\rho h(x)Y+\gamma(x)+\kappa A+\mu.
\end{aligned}
\tag{4.2}
\]

Since \(\rho\kappa\ne0\), recover successively

\[
\begin{aligned}
x&=F_z,\\
Y&=\rho^{-1}\bigl(F_y-\alpha(x)\bigr),\\
A&=\kappa^{-1}
\bigl(F_w-\rho h(x)Y-\gamma(x)-\mu\bigr),\\
w&=A-\frac12x^2,\\
y&=Y-h(x)w.
\end{aligned}
\tag{4.3}
\]

Finally \(F_x=z+R(x,y,w)\) for a polynomial \(R\), so \(z\) is recovered
polynomially. Thus \(\nabla\psi_{\kappa,\mu}\) is a polynomial
automorphism, and direct differentiation gives

\[
\det\operatorname{Hess}\psi_{\kappa,\mu}=-\kappa\rho.
\tag{4.4}
\]

The normal form was obtained after scalar extension. Polynomial
invertibility descends to \(K\): the base-changed coordinate-ring map is an
isomorphism, and faithful flatness reflects that isomorphism.

## 5. Revised reverse-Schur frontier

There is no quadratic scalar-pivot survivor with an identically singular
reduced Hessian pencil:

- rank zero is the affine obstruction HC4RSD7;
- ranks three and four are excluded by HC4RSD8;
- rank two is the triangular classification HC4RSD9;
- rank one is the triangular classification HC4RSD10.

HC4RSD11--HC4RSD16 additionally close every quadratic zero-corner parent
and reduce every nonlinear quadratic nonzero-corner pencil to HC2 or
exactly the JC2 cotangent packet. The remaining scalar mechanisms require a
pencil direction of degree at least three, besides the affine embedding of
direct HC4. Moving matrix-pivot planes and genuinely mixed/coisotropic
reductions also remain open.

## 6. Reproduction

Run:

~~~bash
.venv/bin/python scripts/verify_hc4_quadratic_rank_one_pivots.py
# committed `HC4RSD10` stage artifact only, without symbolic replay:
.venv/bin/python scripts/verify_hc4_quadratic_rank_one_pivots.py --audit-existing-only
~~~

The checker verifies the universal passive block faces, the rank-one frame
determinant, the all-degree normal-form determinant, and the triangular
gradient recovery.

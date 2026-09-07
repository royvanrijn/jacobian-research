# Parameter-moving affine kernel pencils for `HC4`

## Status

This note completes the affine-in-the-reduced-variables scalar branch begun
in [`HC4RSD2`](HC4_AFFINE_MOVING_KERNEL_PENCILS.md). The kernel line may now
move with the pencil parameter.

> **Theorem `HC4RSD3` (parameter-moving affine-kernel obstruction).** Let
> \(K\) have characteristic zero and put
>
> \[
> \Phi(t,x)=\frac{\lambda}{2}t^2+tA(x)+B(x),
> \qquad x\in K^4.
> \]
>
> Assume
>
> \[
> \det\operatorname{Hess}\Phi=c\in K^\times,
> \qquad
> \det\operatorname{Hess}_x(B+sA)=0.                 \tag{0.1}
> \]
>
> If the primitive polynomial generator of the generic kernel line is
> affine in \(x\), then its projective line is independent of \(s\).
> Consequently `HC4RSD1` or `HC4RSD2` applies, and every scalar Schur
> descendant has injective gradient. In particular, no parameter-moving
> affine kernel pencil survives the bordered-unit gate.

No degree bound on \(A\) or \(B\) is assumed. The theorem concerns the
primitive generator itself being affine in \(x\); kernel generators
nonlinear in \(x\) remain outside its scope.

The exact terminal identities are replayed by
[`scripts/verify_hc4_parameter_moving_affine_kernel_pencils.py`](scripts/verify_hc4_parameter_moving_affine_kernel_pencils.py).
Its generated ledger is
[`artifacts/generated-results/hc4_parameter_moving_affine_kernel_pencils.json`](artifacts/generated-results/hc4_parameter_moving_affine_kernel_pencils.json).

## 1. The parameter degree is at most one

Write

\[
M(s)=\operatorname{Hess}_x(B+sA),\qquad g=\nabla A.
\]

The bordered unit forces \(M\) to have generic rank three. As in `HC4RSD2`,
the vector \(\operatorname{adj}(M)g\) is unimodular, and symmetry gives

\[
\operatorname{adj}(M)=\epsilon vv^{\mathsf T},
\qquad
\epsilon\in K^\times,\qquad
g^{\mathsf T}v=\alpha\in K^\times.                       \tag{1.1}
\]

Every entry of \(\operatorname{adj}(M)\) has degree at most three in \(s\).
If the primitive generator \(v\) has parameter degree \(d\), choose a
component of degree \(d\) in (1.1). Its diagonal adjugate entry has degree
\(2d\), so

\[
2d\le3,\qquad d\le1.                                    \tag{1.2}
\]

Degree zero is exactly the parameter-independent branch of `HC4RSD2`.
It remains to exclude

\[
v(s,x)=a_0+L_0x+s(a_1+L_1x).                             \tag{1.3}
\]

## 2. Piola compression pencils

The adjugate Piola identity applied to (1.1) gives

\[
D_vv+(\operatorname{div}v)v=0.                           \tag{2.1}
\]

Applying the affine classification over \(K(s)\), and then taking
coefficients in \(s\), yields

\[
\begin{gathered}
\operatorname{tr}L_0=\operatorname{tr}L_1=0,\\
L_0^2=L_1^2=0,
\qquad L_0L_1+L_1L_0=0,\\
L_0a_0=0,
\qquad L_0a_1+L_1a_0=0,
\qquad L_1a_1=0.                                         \tag{2.2}
\end{gathered}
\]

For every scalar value of \(s\), the pairing in (1.1) makes \(v(s)\)
unimodular. The rank-two affine Piola orbit is therefore impossible, so

\[
\operatorname{rank}(L_0+sL_1)\le1.                       \tag{2.3}
\]

If the two nonzero linear parts are proportional, a translation of \(s\)
makes one pencil member constant in \(x\). At that parameter, the reduced
Hessian has rank three. The constant kernel vector is killed by both
\(\operatorname{Hess}A\) and the recentered \(\operatorname{Hess}B\), and
the coefficient of \(s\) is forced onto the same kernel line. Thus the
primitive line has parameter degree zero.

The exceptional proportional corner \(L_1=0\), where no finite translation
kills the linear part, is treated in Section 3. Otherwise \(L_0,L_1\) are
independent. A two-dimensional linear space consisting entirely of
rank-at-most-one matrices is a compression space: writing
\(L_i=u_i\ell_i^{\mathsf T}\), either the \(u_i\) are proportional or the
\(\ell_i\) are proportional. These are respectively the common-image and
common-covector cases.

In the common-image case write

\[
L_i=u\ell_i^{\mathsf T},\qquad \ell_0,\ell_1
\text{ independent}.                                    \tag{2.4}
\]

The highest kernel coefficient and \(g^{\mathsf T}v_1=0\) imply
\(L_1^{\mathsf T}g=0\), hence \(D_uA=0\). Equation (1.1) then gives
\(D_{a_0}A=\alpha\) and \(D_{a_1}A=0\). Therefore
\(\operatorname{Hess}(A)v_0=0\), and the middle kernel equation gives
\(\operatorname{Hess}(B)v_1=0\). Since
\(\operatorname{Hess}(B)v_0=0\) and \(\operatorname{Hess}(B)\) has rank
three, \(v_0\) and \(v_1\) span the same rational line. Again the primitive
kernel has parameter degree zero.

## 3. The constant-at-infinity corner

Suppose \(L_1=0\) and \(L_0\ne0\). Normalize

\[
v_0=(z,1,0,0)^{\mathsf T},\qquad
v_1=(b_1,b_2,0,b_4)^{\mathsf T}.                          \tag{3.1}
\]

Put \(r=x-yz\). The complete \(v_0\)-kernel potential and the complete
solution of \(g^{\mathsf T}v_0=\alpha\) are

\[
\begin{aligned}
B&=yC(z)+rC'(z)+G(z,w),\\
A&=\alpha y+F(r,z,w).                                    \tag{3.2}
\end{aligned}
\]

At \(s=0\), rank three requires \(C''G_{ww}\ne0\). Direct multiplication
gives

\[
\operatorname{Hess}(B)v_1=
\begin{pmatrix}
0\\0\\
(b_1-b_2z)C''+b_4G_{wz}\\
b_4G_{ww}
\end{pmatrix}.                                           \tag{3.3}
\]

The middle kernel equation and
\(\operatorname{Hess}(A)v_0=-F_r e_z\) force \(b_4=0\) and

\[
F_r=(b_1-b_2z)C''.                                       \tag{3.4}
\]

On the other hand, \(g^{\mathsf T}v_1=0\) says

\[
b_2\alpha+(b_1-b_2z)F_r=0.
\]

Substitution of (3.4) yields

\[
\boxed{b_2\alpha+(b_1-b_2z)^2C''=0.}                    \tag{3.5}
\]

If \(b_2\ne0\), a nonconstant square times a polynomial cannot be the
nonzero constant \(-b_2\alpha\). If \(b_2=0\), a nonzero \(v_1\) has
\(b_1\ne0\), and (3.5) contradicts \(C''\ne0\). Thus this corner is empty.

## 4. The common-covector terminal chart

It remains to take

\[
L_i=u_i\ell^{\mathsf T},\qquad u_0,u_1\text{ independent}. \tag{4.1}
\]

Equation (2.2) puts \(u_0,u_1,a_0,a_1\) in the three-plane
\(W=\ker\ell\). If \(u_1,a_1\) are independent, then \(A\) is invariant in
their two-plane. On the one-dimensional quotient of \(W\), the bordered
pairing becomes

\[
(a+bz)A_r=\alpha.                                        \tag{4.2}
\]

Both factors must be units, so \(b=0\) and \(A\) is affine in \(r\).
Consequently \(\operatorname{Hess}(A)v_0=0\), and the rank-three argument
again collapses \(v_1\) onto \(v_0\).

Hence \(a_1\) is proportional to \(u_1\). Translate the transverse
coordinate so that \(a_1=0\). Joint unimodularity leaves exactly two normal
forms:

\[
v=(sz,z,1,0)^{\mathsf T}
\quad\text{or}\quad
v=(1+sz,z,0,0)^{\mathsf T}.                               \tag{4.3}
\]

For the second, \(g^{\mathsf T}v_1=0\) gives \(D_{e_1}A=0\), while
\(g^{\mathsf T}v_0=\alpha\) becomes \(zD_{e_2}A=\alpha\), impossible for a
polynomial.

For the first, use coordinates \((x,y,r,z)\), put \(q=y-zr\), and write

\[
v=(sz,z,1,0)^{\mathsf T}.                                \tag{4.4}
\]

The bordered pairings and the three kernel coefficients integrate exactly
to

\[
\begin{aligned}
A&=\alpha r+zp'(z)q+h(z),\\
B&=rC(z)+qC'(z)+p(z)x+g(z).                              \tag{4.5}
\end{aligned}
\]

Every member \(B+sA\) is linear in \(x,y,r\). Its Hessian therefore has
the form

\[
\begin{pmatrix}
0_{3\times3}&u\\
u^{\mathsf T}&*
\end{pmatrix},                                           \tag{4.6}
\]

which has rank at most two. This contradicts the rank-three consequence of
the bordered unit and completes the proof.

## 5. Reproduction and next frontier

Run:

```bash
.venv/bin/python scripts/verify_hc4_parameter_moving_affine_kernel_pencils.py
# cleanup only: verify committed inputs and exact boundary
.venv/bin/python scripts/verify_hc4_parameter_moving_affine_kernel_pencils.py --audit-existing-only
```

The command verifies (3.3), both jointly unimodular common-covector normal
forms, the complete terminal integral (4.5), all sixteen \(3\)-by-\(3\)
Hessian minors, and the vanishing bordered determinant.

The cleanup-only mode hash-checks the committed ledger and imported equation
helper without importing SymPy, replaying those identities, or rewriting the
artifact. It also requires the ledger to retain the exact boundary at
primitive kernel generators nonlinear in the reduced variables.

Together, `HC4RSD1`--`HC4RSD3` close every singular scalar pencil whose
primitive kernel generator is affine in the four reduced variables. The
subsequent `HC4RSD4`--`HC4RSD5` results also close every fixed primitive
two-component kernel in a constant support plane, without a degree bound.
The scalar singular frontier now consists of fixed nonlinear kernels with
three or four components and parameter-moving nonlinear generators.

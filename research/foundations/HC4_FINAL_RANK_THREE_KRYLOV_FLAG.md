# Fiberwise Krylov-flag rigidity in the final rank-three `[4]` HC4 stratum

## Status

This note continues `HC4RSD71--73`.  The smooth affine gradient-image case is
already excluded.  Here we spend the new affine orbit formula

\[
\nabla\psi(u,t)=H_0(u)+t\ell(u)
\]

inside the quotient determinant identity of `HC4RSD71`.

> **Theorem HC4RSD74 — triangular Krylov variation.**  On the generic
> Gauss-rank-two locus of the final `[4]` packet, choose quotient coordinates
> and a target tangent frame as in `HC4RSD71`.  Normalize the cyclic quotient
> data so that
> \[
> c=e_1,\qquad M_0c=e_2,\qquad M_0^2c=e_3,
> \qquad r=e_3^{\mathsf T}.
> \]
> If
> \[
> M(t)=M_0+tC
> \]
> is the variation along a collapsed source fiber, then
> \[
> C_{21}=C_{31}=C_{32}=0.
> \]
> Hence
> \[
> C=\begin{pmatrix}
> *&*&*\\0&*&*\\0&0&*
> \end{pmatrix}.
> \]
> In particular the variation preserves the complete cyclic flag
> \[
> \langle c\rangle\subset
> \langle c,M_0c\rangle\subset K^3.
> \]

Combined with HC4RSD72, the induced two-dimensional ruling-direction map has
full rank but triangular first derivative in the canonical Krylov flag.

## 1. Quotient determinant revisited

Let `u=(u1,u2,u3)` be local coordinates on the quotient by the source
quasi-translation and let `t` be its affine fiber parameter.  Thus

\[
F=\nabla A=F(u),
\qquad
x(u,t)=x_0(u)+t k(u).
\]

Put

\[
B_0=F_u.
\]

The strengthened orbit theorem `HC4RSD72` gives

\[
H:=\nabla\psi=H_0(u)+t\ell(u),
\qquad \ell=Sk.
\tag{1.1}
\]

The line `ell` is the Gauss-kernel direction of the developable gradient image.
Therefore its derivative is tangent to `Y`; in a frozen target frame `[B0,n]`
we may write

\[
H_u=B_0(M_0+tC)+nr^{\mathsf T},
\qquad
H_t=\ell=B_0c.                                      \tag{1.2}
\]

All of `M0,C,r,c` are independent of the fiber parameter `t` at the chosen
quotient point.

The quasi-translation flow has Jacobian

\[
J\phi_t=I+tJk,
\qquad
\det(I+tJk)=1,
\]

so the source volume factor in the determinant calculation is also independent
of `t`.

The quotient form of the constant Hessian determinant is therefore

\[
-r^{\mathsf T}\operatorname{adj}
\bigl(I+s(M_0+tC)\bigr)c
=\Gamma s^2,
\qquad \Gamma\ne0,                                  \tag{1.3}
\]

identically in both `s` and `t`.

## 2. Cyclic normalization

At `t=0`, coefficient comparison in (1.3) gives

\[
r^{\mathsf T}c=0,
\qquad
r^{\mathsf T}M_0c=0,
\qquad
r^{\mathsf T}M_0^2c\ne0.
\tag{2.1}
\]

Thus

\[
c,\ M_0c,\ M_0^2c
\]

are a basis.  Rescale the last vector/covector and use this Krylov basis, so

\[
c=e_1,
\qquad M_0c=e_2,
\qquad M_0^2c=e_3,
\qquad r=e_3^{\mathsf T}.                            \tag{2.2}
\]

Then

\[
M_0=
\begin{pmatrix}
0&0&\alpha_0\\
1&0&\alpha_1\\
0&1&\alpha_2
\end{pmatrix}.                                      \tag{2.3}
\]

## 3. The fiber parameter forces triangularity

For a `3 x 3` matrix `M`,

\[
\operatorname{adj}(I+sM)
=I+s(\operatorname{tr}M\,I-M)+s^2\operatorname{adj}M.
\tag{3.1}
\]

The coefficient of `s` in (1.3), using `r^T c=0`, is simply

\[
r^{\mathsf T}(M_0+tC)c=0.
\]

The `t` coefficient gives

\[
C_{31}=0.                                           \tag{3.2}
\]

The coefficient of `s^2` is constant and nonzero.  A direct determinant
calculation in the Krylov basis gives

\[
r^{\mathsf T}\operatorname{adj}(M_0+tC)c
=1+t(C_{21}+C_{32})
+t^2(C_{21}C_{32}-C_{22}C_{31}).                    \tag{3.3}
\]

Using (3.2), independence of `t` gives

\[
C_{21}+C_{32}=0,
\qquad
C_{21}C_{32}=0.
\]

Hence, over a characteristic-zero domain,

\[
\boxed{C_{21}=C_{31}=C_{32}=0.}                    \tag{3.4}
\]

This proves HC4RSD74.

## 4. Geometry of the flag

The first vector `c` represents the selected Gauss-line direction `ell`.
Equation (3.4) says:

1. differentiating `ell` along its own Gauss line changes only its scale;
2. along the first transverse Krylov direction, the derivative of `ell` stays
   inside the osculating plane `span(c,M0 c)`;
3. only the second transverse direction can create the third Krylov component.

Thus the full-rank two-dimensional direction map of HC4RSD72 is not arbitrary:
it comes with a canonical triangular osculating flag.  In the normalized chart
its transverse differential has the form

\[
\begin{pmatrix}
C_{22}&C_{23}\\
0&C_{33}
\end{pmatrix}
\]

modulo the ruling line.  HC4RSD72 implies its determinant is nonzero, so

\[
C_{22}C_{33}\ne0.                                   \tag{4.1}
\]

The remaining singular/focal classification should therefore be run against a
*flagged* two-parameter congruence of Gauss lines, rather than a general
Piontkowski rank-two developable threefold.
# The two-line quartic-denominator Schur packet

## Status

This note proves `HC4NHM4`: the clean generic-corank-one packet with minimal
denominator \(P=x^3y\) and determinant incidence

\[
\det\operatorname{Hess}(h_5)=x^7y^2
\]

is nonempty but rigid. Up to a linear normalization, every such Hessian--Schur
pair belongs to one four-parameter quintic family and its Schur cubics form a
two-dimensional space. The genuinely two-component sections are exactly the
ones with one specified coefficient nonzero.

This is a classification of the leading ternary Hessian--Schur packet, not a
four-variable constant-Hessian completion or an `HC(4)` counterexample.

Replay the exact identities with

~~~bash
.venv/bin/python scripts/verify_hc4_two_line_quartic_denominator_packet.py
~~~

The normalization theorem is
[`HC4_NONREDUCED_HESSIAN_SCHUR_MODULE.md`](HC4_NONREDUCED_HESSIAN_SCHUR_MODULE.md).
The local moving-kernel and constant-kernel gates are in
[`HC4_DIRECT_SEPTUPLE_LINEAR_HESSIAN_GATE.md`](HC4_DIRECT_SEPTUPLE_LINEAR_HESSIAN_GATE.md),
and the one-line quartic-denominator partition is closed in
[`HC4_OCTUPLE_NONUPLE_LINEAR_HESSIAN_GATE.md`](HC4_OCTUPLE_NONUPLE_LINEAR_HESSIAN_GATE.md).

Work over an algebraically closed field of characteristic zero. Put

\[
C=\operatorname{Hess}(h_5),\qquad d=\nabla s_3,
\tag{0.1}
\]

and assume that \(C\) has generic corank one along both lines, that

\[
\det C=\delta x^7y^2,\qquad \delta\ne0,
\tag{0.2}
\]

and that the componentwise-minimal denominator of \(C^{-1}d\) is exactly

\[
P=x^3y.
\tag{0.3}
\]

Thus both \(x\) and \(y\) are essential and

\[
\det C\mid d^{\mathsf T}\operatorname{adj}(C)d.
\tag{0.4}
\]

## 1. The septuple component has constant tangent kernel

The global cleared vector \(e=PC^{-1}d\) has degree three. On the normalized
line \(x=0\), its saturated kernel degree is therefore at most three. One
must not divide globally by the other denominator factor \(y\): the local
residue may have a pole at the intersection point \(x=y=0\).

The degree-one moving-line argument and the degree-two line/conic arguments
of `HC4NHM2` require only that \(x^2\) divide the Hessian determinant. The
defect-free degree-three kernel is excluded by the perfect-square binary
Hessian calculation of `HC4NHM3`. These results eliminate every positive
kernel degree here. The kernel along \(x=0\) is constant.

A constant kernel transverse to \(x=0\) makes the determinant zero once
\(x^7\) divides it. Normalize the remaining tangent kernel to
\(\partial_z\), and use the factorial boundary expansion

\[
h_5=F+xG+\frac{x^2}{2}H+\frac{x^3}{6}J
       +\frac{x^4}{24}K+\frac{x^5}{120}D.
\tag{1.1}
\]

The boundary equations give

\[
F=\alpha y^5,\qquad G=A y^4.
\tag{1.2}
\]

The \(\alpha\ne0\) coefficient ladder loses all \(z\)-dependence and has
zero determinant. Thus \(\alpha=0\), while generic boundary rank forces
\(A\ne0\). The remaining septuple ladder gives

\[
H=h_0y^3,qquad J=j_0y^2,qquad K=By+\Gamma z,
\tag{1.3}
\]

with \(\Gamma\ne0\), and its complete determinant is

\[
\det C=-\frac{\Gamma^2x^7}{108}
       \left(36Ay^2+9h_0xy+j_0x^2\right).
\tag{1.4}
\]

Equation (0.2) forces \(h_0=j_0=0\). Consequently every packet has

\[
\boxed{
h_5=Axy^4+\frac{x^4}{24}(By+\Gamma z)+\frac D{120}x^5,
\qquad A\Gamma\ne0.
}
\tag{1.5}
\]

Conversely,

\[
\det\operatorname{Hess}(h_5)
=-\frac{A\Gamma^2}{3}x^7y^2,
\tag{1.6}
\]

and the Hessian has generic rank two on both lines. Thus (1.5) is the
complete Hessian boundary, not merely a necessary jet.

## 2. Complete Schur-space classification

For (1.5), the Hessian and the entries of its adjugate needed below are

\[
C=
\begin{pmatrix}
\frac D6x^3+\frac{x^2}{2}(By+\Gamma z)&4Ay^3+\frac B6x^3&\frac\Gamma6x^3\\
4Ay^3+\frac B6x^3&12Axy^2&0\\
\frac\Gamma6x^3&0&0
\end{pmatrix},
\tag{2.1}
\]

\[
\begin{aligned}
(\operatorname{adj}C)_{22}&=-\frac{\Gamma^2}{36}x^6,\\
(\operatorname{adj}C)_{23}&=\frac{\Gamma x^3}{36}(24Ay^3+Bx^3),\\
(\operatorname{adj}C)_{13}&=-2A\Gamma x^4y^2,\\
(\operatorname{adj}C)_{33}&=-16A^2y^6+O(x^3).
\end{aligned}
\tag{2.2}
\]

Let \(d=(d_x,d_y,d_z)=\nabla s_3\). If \(d_z\ne0\), its \(x\)-order is at
most two. In the adjugate quadratic form, the unique term of least
\(x\)-order is then \(-16A^2y^6d_z^2\): all cross terms start at least three
orders later. It cannot be divisible by \(x^7\). Hence

\[
d_z=0.
\tag{2.3}
\]

The numerator now reduces exactly to

\[
d^{\mathsf T}\operatorname{adj}(C)d
=-\frac{\Gamma^2}{36}x^6d_y^2.
\tag{2.4}
\]

Divisibility by \(x^7y^2\) gives \(xy\mid d_y\). Since \(d_y\) is
quadratic, gradient integrability and homogeneity give

\[
\boxed{s_3=ax y^2+b x^3.}
\tag{2.5}
\]

Conversely, for every \(a,b\), direct substitution gives

\[
\boxed{
\frac{d^{\mathsf T}\operatorname{adj}(C)d}{\det C}
=\frac{a^2}{3A}x.
}
\tag{2.6}
\]

Thus (2.5) is the complete Schur space.

## 3. The genuinely two-component section

For (2.5), put \(P=x^3y\). The primitive cleared vector is

\[
\boxed{
e=PC^{-1}d=
\begin{pmatrix}
0\\[2mm]
\dfrac{a x^3}{6A}\\[2mm]
\dfrac{12Aa y^3+108Ab x^2y-Ba x^3}{6A\Gamma}
\end{pmatrix}.
}
\tag{3.1}
\]

It satisfies

\[
Ce=Pd,qquad d^{\mathsf T}e
=P\left(\frac{a^2}{3A}x\right).
\tag{3.2}
\]

If \(a\ne0\), then \(e|_{x=0}\ne0\) and \(e|_{y=0}\ne0\), and its entries
have no common factor. Therefore the denominator is exactly \(x^3y\). If
\(a=0\), the minimal denominator drops to \(x\); this is the older
constant-kernel channel and is not part of the clean two-component packet.

We have proved:

> **Theorem `HC4NHM4` -- Two-line quartic-denominator Schur classification.**
> Under (0.1)--(0.4), every packet is linearly equivalent to (1.5) and
> (2.5), with \(A\Gamma a\ne0\). Conversely, every such choice is a
> generic-corank-one Hessian--Schur packet with minimal denominator
> \(P=x^3y\) and determinant \(-A\Gamma^2x^7y^2/3\).

## 4. Prolongation target

Return to the four-variable rank-three degree-five face, with kernel
coordinate \(t\). The nonaligned quartic term is

\[
h_4=t\,s_3+r_4(x,y,z).
\tag{4.1}
\]

The next determinant face is

\[
\det(C)\,\partial_t^2h_3
-d^{\mathsf T}\operatorname{adj}(C)d=0.
\tag{4.2}
\]

Equation (2.6) therefore fixes

\[
\partial_t^2h_3=\frac{a^2}{3A}x,
\qquad
h_3=\frac{a^2}{6A}xt^2+tq_2(x,y,z)+r_3(x,y,z).
\tag{4.3}
\]

Retaining arbitrary ternary \(r_4,q_2,r_3\) and the nondegenerate quadratic
\(h_2\), the next Hessian-determinant face is computed in
[`HC4_TWO_LINE_QUARTIC_DENOMINATOR_PROLONGATION.md`](HC4_TWO_LINE_QUARTIC_DENOMINATOR_PROLONGATION.md).
Its immutable coefficient \(-\Gamma^2a^3/(54A)\) excludes every genuinely
two-component prolongation. The leading module packet classified here remains
nonempty; it simply cannot be completed to a constant-Hessian four-variable
potential.

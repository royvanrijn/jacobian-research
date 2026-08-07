# Fixed-kernel closure for the final rank-three `[4]` HC4 stratum

## Status and scope

Continue `HC4RSD64`.  Let

\[
S=\operatorname{Hess}\psi,
\qquad
T=\operatorname{Hess}A,
\qquad
\det(S+sT)=\delta\in K^*,
\]

and suppose the relative nilpotent `N=S^{-1}T` has Jordan type `[4]`, so
`rank T=3` generically.

> **Theorem HC4RSD65 — fixed top kernel closes globally.**
> If the one-dimensional kernel of `T` is a constant line, then after a
> constant linear coordinate change the packet has the scalar form
>
> \[
> \psi=xP(u_1,u_2,u_3)+Q(u_1,u_2,u_3),
> \qquad
> A=A(u_1,u_2,u_3),
> \tag{0.1}
> \]
>
> and the constant-Hessian condition forces `P` to have a nonzero constant
> direction.  Hence the packet reduces to `HC2` or to the exact `JC2`
> cotangent endpoint.

> **Corollary HC4RSD65a — homogeneous rank-three `[4]` closes in every
> degree.**  If `A` is homogeneous, then the Gordan--Noether theorem makes
> `ker Hess A` a constant line.  Therefore HC4RSD65 applies.

Thus every surviving `[4]` packet must be genuinely nonhomogeneous and have a
genuinely moving top kernel line.

## 1. Normalize the fixed kernel

Assume

\[
\ker T=K e_x.
\]

Since `T=Hess A`, this says

\[
A_x=\text{constant}.
\]

Discarding an affine term, write

\[
A=A(u),\qquad u=(u_1,u_2,u_3).
\]

Then

\[
T=
\begin{pmatrix}
0&0\\
0&H
\end{pmatrix},
\qquad
H=\operatorname{Hess}_u A,
\qquad
\det H\ne0
\tag{1.1}
\]

generically.

Write the metric Hessian in the same splitting as

\[
S=
\begin{pmatrix}
\psi_{xx}&p^{\mathsf T}\\
p&C
\end{pmatrix}.
\tag{1.2}
\]

Then

\[
S+sT=
\begin{pmatrix}
\psi_{xx}&p^{\mathsf T}\\
p&C+sH
\end{pmatrix}.
\]

The coefficient of `s^3` in its determinant is

\[
[s^3]\det(S+sT)
=\psi_{xx}\det H.
\tag{1.3}
\]

The determinant is independent of `s`, and `det H` is nonzero in the
fraction field.  Therefore

\[
\psi_{xx}=0.
\tag{1.4}
\]

Hessian integrability now gives, up to an affine term,

\[
\psi=xP(u)+Q(u),
\tag{1.5}
\]

which proves the normal form (0.1).

## 2. Reuse the global scalar obstruction

For every scalar `s`, put

\[
D_s(u)=Q(u)+sA(u).
\]

Then

\[
\psi+sA=xP(u)+D_s(u).
\]

Its Hessian is

\[
\begin{pmatrix}
0&(\nabla P)^{\mathsf T}\\
\nabla P&\operatorname{Hess}D_s+x\operatorname{Hess}P
\end{pmatrix}.
\tag{2.1}
\]

The determinant is the scalar reverse-Schur block studied in
`HC4RSD52--53`.  Since it is the same nonzero constant `delta`, the exact
algebraic bridge there gives

\[
V(P_{u_1},P_{u_2},P_{u_3})=\varnothing
\tag{2.2}
\]

and

\[
(\nabla P)^{\mathsf T}
\operatorname{adj}(\operatorname{Hess}P)
\nabla P=0.
\tag{2.3}
\]

By the smooth-developable theorem `HC4RSD52`, there is a nonzero constant
vector `v` such that

\[
D_vP=0.
\tag{2.4}
\]

Thus the first nontrivial member of the `[4]` Jordan chain already carries a
second fixed direction.  The fixed-ruling reduction `HC4RSD20` / `HC4RSD53`
then leaves only `HC2` or the exact `JC2` cotangent endpoint.

This proves HC4RSD65.

## 3. Homogeneous corollary

If `A` is homogeneous and

\[
\det\operatorname{Hess}A=0,
\]

the Gordan--Noether theorem in four variables says that `A` is a cone: after
a constant linear coordinate change it is independent of one variable.  In
the generic rank-three case the Hessian kernel is therefore exactly that
constant line.  HC4RSD65 applies immediately.

Consequently no homogeneous moving `[4]` packet survives, regardless of total
degree.

## 4. Remaining frontier

After HC4RSD64--65 the only possible moving nilpotent obstruction has all of
the following features simultaneously:

1. `rank Hess A=3`;
2. `A` is nonhomogeneous;
3. the canonical line `im adj(Hess A)` moves;
4. that line is projectively straight along its own flow;
5. it extends to the divergence-free nested cofactor flag
   `im C3 subset im C2 subset im C1`.

General nonhomogeneous singular-Hessian polynomials in four variables can have
moving kernels, so item 3 cannot be removed by Hesse--Gordan--Noether alone.
The additional cofactor-pencil identities are the structure that remains to
be exploited.

# Direct degree-filtration reduction for unrestricted HC4

> **Archived route.** Superseded by the active direct homogeneous-filtration
> chain; `HC4DF1` is retained for provenance only.

## Status

This note attacks `HC4` without assuming a pre-existing relative-nilpotent pencil.
It extracts such structure directly from the ordinary homogeneous filtration of a
constant-Hessian polynomial.

Let

\[
\Psi=\Psi_2+\Psi_3+\cdots+\Psi_d\in K[x_1,x_2,x_3,x_4],
\qquad
\det\operatorname{Hess}\Psi=\delta\in K^\times,
\]

with each \(\Psi_j\) homogeneous of ordinary degree \(j\).  Write
\(H_j=\operatorname{Hess}\Psi_j\).

> **Theorem HC4DF1 — top cone and first Schur descent.**
> Assume \(d\ge3\) and the generic rank of \(H_d\) is three.  Then, after a
> constant linear change of coordinates \((u,z)=(x_1,x_2,x_3,x_4)\),
>
> \[
> \Psi_d=f(u),
> \]
>
> and the next two homogeneous layers have the forced shapes
>
> \[
> \Psi_{d-1}=z\,g(u)+h(u),
> \]
>
> \[
> \Psi_{d-2}=\frac{z^2}{2}q(u)+z\,r(u)+s(u).
> \]
>
> Moreover
>
> \[
> q\,\det\operatorname{Hess}f
> =
> (\nabla g)^T\operatorname{adj}(\operatorname{Hess}f)\nabla g.
> \tag{DF1}
> \]

Thus unrestricted HC4 begins with the same reverse-Schur mechanism that appears
in the relative-nilpotent analysis, but now it is forced by the ordinary degree
filtration itself.

## Proof

Scale all variables by \(t\).  Since \(\Psi_j(tx)=t^j\Psi_j(x)\),

\[
\operatorname{Hess}\Psi(tx)
=\sum_{j=2}^d t^{j-2}H_j(x).
\]

Its determinant is the constant \(\delta\), hence every positive coefficient in
\(t\) vanishes.

The highest coefficient is

\[
\det H_d=0.
\]

Because \(\Psi_d\) is a homogeneous form in four variables, the
Gordan--Noether theorem in projective dimension three says that a homogeneous
form with vanishing Hessian determinant is a cone.  Under the rank-three
hypothesis its kernel is one-dimensional and constant.  Choose coordinates so
that this kernel is \(e_z\).  Then

\[
\Psi_d=f(u),
\qquad
A:=\operatorname{Hess}_u f
\]

has generic rank three.

In these coordinates

\[
H_d=
\begin{pmatrix}
A&0\\
0&0
\end{pmatrix}.
\]

The next coefficient in the scaled determinant contains three copies of
\(H_d\) and one copy of \(H_{d-1}\).  Since
\(\operatorname{adj}(H_d)\) is supported only in the \(zz\)-entry, it is

\[
\det(A)\,(\Psi_{d-1})_{zz}=0.
\]

Generic invertibility of \(A\) gives

\[
(\Psi_{d-1})_{zz}=0,
\]

hence

\[
\Psi_{d-1}=z\,g(u)+h(u).
\]

Now write the relevant blocks of \(H_{d-1}\) as

\[
H_{d-1}=
\begin{pmatrix}
*&\nabla g\\
(\nabla g)^T&0
\end{pmatrix}.
\]

At the next power of \(t\), the only contributions are one
\(H_{d-2}\) in the kernel--kernel position and the quadratic Schur term from
two copies of \(H_{d-1}\).  Thus

\[
\det(A)\,(\Psi_{d-2})_{zz}
-(\nabla g)^T\operatorname{adj}(A)\nabla g=0.
\]

The right-hand side is independent of \(z\), so
\((\Psi_{d-2})_{zz}\) is independent of \(z\).  Therefore

\[
\Psi_{d-2}=\frac{z^2}{2}q(u)+z\,r(u)+s(u)
\]

and equation `(DF1)` follows.

## Weighted diagonal

The three displayed terms

\[
f(u),\qquad z g(u),\qquad \frac{z^2}{2}q(u)
\]

have one common weighted degree if

\[
\operatorname{wt}(u_i)=1,
\qquad
\operatorname{wt}(z)=2.
\]

Indeed

\[
\deg f=d,
\qquad
\deg g=d-2,
\qquad
\deg q=d-4.
\]

This suggests the general filtration statement:

> **Weighted-descent target.**  After choosing the constant cone direction of
> the top homogeneous part, every successive diagonal term has the form
>
> \[
> \frac{z^k}{k!}a_k(u),
> \qquad
> \deg a_k=d-2k,
> \]
>
> and the determinant equations are precisely the coefficient equations of a
> single weighted reverse-Schur face.

Proving this target to all orders would reduce an arbitrary rank-three HC4
polynomial to the relative-nilpotent/scalar machinery already developed in the
repository.  This is the direct-HC4 bridge currently worth attacking.

## Why this is potentially stronger than the old degree census

The parameter is no longer the total degree.  The top cone produces a fixed
weight-two transverse variable automatically.  The same Schur equation then
propagates down the entire homogeneous filtration.  If the weighted-descent
target holds, all degrees are treated simultaneously.

The rank-\(\le2\) top-Hessian branch should be handled separately by the
small-rank Hessian classifications of de Bondt; the difficult unrestricted
branch is therefore precisely the rank-three cone considered here.

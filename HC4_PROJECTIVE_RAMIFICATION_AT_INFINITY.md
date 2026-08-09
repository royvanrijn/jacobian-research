# Direct projective form of `HC4`: all ramification is at infinity

## Statement

Let

\[
\Psi\in K[x_1,x_2,x_3,x_4],\qquad
\deg\Psi=D>2,\qquad
\det\operatorname{Hess}\Psi=\delta\in K^\times,
\]

and put `m=D-1`, `r=D-2`.

For each component of the gradient, let

\[
G_i(X,T)=T^m\,\Psi_{x_i}(X/T)
\tag{0.1}
\]

be its degree-`m` homogenization.  Consider the projective rational map

\[
\overline F:\mathbb P^4\dashrightarrow\mathbb P^4,
\qquad
[X:T]\longmapsto[G_1:G_2:G_3:G_4:T^m].
\tag{0.2}
\]

> **Theorem HC4-DIR4 — pure projective ramification.**
> The homogeneous Jacobian determinant of the five coordinate forms in
> (0.2) is
>
> \[
> \boxed{
> \det J(G_1,G_2,G_3,G_4,T^m)
> =m\delta\,T^{5r}.
> }
> \tag{0.3}
> \]
>
> Hence the projective extension of the gradient is unramified away from the
> hyperplane at infinity `T=0` wherever it is defined, and its entire
> projective ramification divisor is supported there.
>
> Its base locus is
>
> \[
> \boxed{
> \operatorname{Bs}(\overline F)
> =V(T,\partial_1\Psi_D,\ldots,\partial_4\Psi_D).
> }
> \tag{0.4}
> \]
>
> Thus the direct `HC4` problem is equivalently a problem of controlling the
> base scheme and infinitely-near ramification of an otherwise unramified
> projective gradient map at infinity.

---

## 1. Proof of the Jacobian identity

All five coordinates in (0.2) are homogeneous of degree `m`, so their
Jacobian determinant is homogeneous of degree

\[
5(m-1)=5r.
\]

On the affine chart `T=1`, the Jacobian matrix has block form

\[
\begin{pmatrix}
\operatorname{Hess}\Psi & *\\
0&m
\end{pmatrix}.
\]

Therefore

\[
\left.\det J(G_1,\ldots,G_4,T^m)\right|_{T=1}
=m\det\operatorname{Hess}\Psi
=m\delta.
\]

A homogeneous polynomial of degree `5r` whose dehomogenization at `T=1` is
the nonzero constant `m delta` must equal `m delta T^{5r}`.  This proves
(0.3).

---

## 2. Base locus

At infinity,

\[
G_i(X,0)=\partial_i\Psi_D(X).
\]

The final projective coordinate is zero there as well.  Thus the common zero
scheme is exactly (0.4).

The leading Hessian determinant vanishes, so `Psi_D` is a cone in four
variables by the low-dimensional Hesse--Gordan--Noether theorem.  Consequently
this base locus is always nonempty when `D>2`.

On the generic rank-three branch, after a constant linear change

\[
\Psi_D=f(x_1,x_2,x_3)
\]

and the constant cone vertex

\[
[0:0:0:1:0]
\]

is a base point.  If the ternary projective curve `f=0` is smooth, this is the
**only** base point at infinity.  If `f` is singular, the base locus contains
additional points/structure corresponding exactly to the singular scheme of
`f`.

If the top Hessian rank is at most two, Hesse reduction can be iterated and
the top form depends on at most two constant linear forms; the projective base
locus then has positive dimension.  This is the projective geometry naturally
compatible with the already-known `HC2/JC2` cotangent endpoint.

---

## 3. Relation with the homogeneous filtration

The scaled Hessian polynomial from `HC4-DIR1` is

\[
M(t,x)=\operatorname{Hess}\Psi(tx)
=H_2+tH_3+\cdots+t^rH_D,
\qquad
\det M(t,x)=\delta.
\tag{3.1}
\]

Reversing `t` gives

\[
H_D+tH_{D-1}+\cdots+t^rH_2,
\]

whose determinant is `delta t^{4r}`.  Thus the `t`-adic Schur descent is the
local algebra of the base scheme of (0.2) seen from infinity.

In particular:

- `HC4-DIR2` says a rank-three base point cannot have squarefree ternary
  Hessian ramification;
- `HC4-DIR3` says the first motion away from the cone vertex consumes the
  square-root ramification budget of that ternary Hessian divisor.

So the projective and valuation pictures are two descriptions of the same
obstruction.

---

## 4. Research consequence

The direct route to `HC4` is now:

1. resolve the base scheme (0.4), starting with the rank-three isolated-vertex
   case;
2. use the exact ramification divisor `5r H_infinity` together with the
   gradient/Hessian symmetry to constrain every exceptional divisor;
3. show that an isolated rank-three vertex cannot support the required full
   ramification tower;
4. conclude that the top base locus has positive dimension, hence the leading
   form reduces to at most two active affine directions;
5. identify that residue with the `HC2` or plane-cotangent `JC2` geometry.

This is a direct projective attack on arbitrary `HC4`, not an auxiliary-pencil
reduction.

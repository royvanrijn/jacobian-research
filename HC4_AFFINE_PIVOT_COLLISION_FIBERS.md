# Affine-pivot collision fibers for HC4

## Status

HC4RSD6 gives an exact inverse-Hessian criterion for representing a direct
HC4 candidate by an affine singular scalar Schur pivot. This note imposes
the marked collision and shows that the intersection is empty before any
lower metric-face calculation. The argument then extends from those
singular presentations to every affine zero-corner scalar parent with
constant nonzero bordered Hessian determinant.

> **Theorem HC4RSD7 (affine-pivot collision-fiber obstruction).** Let \(K\)
> have characteristic zero, let
>
> \[
> \psi\in K[x_1,x_2,x_3,x_4],
> \qquad
> H=\operatorname{Hess}\psi,
> \]
>
> and let \(\ell\in K^4\setminus\{0\}\) be a constant covector satisfying
>
> \[
> N_\ell
> =\ell^{\mathsf T}\operatorname{adj}(H)\ell
> \in K^\times.
> \tag{0.1}
> \]
>
> Then the gradient of \(\psi\) is injective on every affine hyperplane
> \(\ell\mathbin{\cdot}x=c\). In particular,
>
> \[
> \nabla\psi(p)=\nabla\psi(q),
> \qquad
> \ell\mathbin{\cdot}p=\ell\mathbin{\cdot}q
> \quad\Longrightarrow\quad
> p=q.
> \tag{0.2}
> \]
>
> Consequently, an affine zero-corner parent
> \(\Phi=tA+B\), with \(A\) affine and constant nonzero bordered Hessian
> determinant, cannot lift a marked nontrivial collision to equal parent
> gradients at a common pivot value. This includes both the singular-pivot
> presentations supplied by HC4RSD6 and the nonsingular exact-remainder
> branch.

The exact identities are replayed by
[scripts/verify_hc4_affine_pivot_collision_fibers.py](scripts/verify_hc4_affine_pivot_collision_fibers.py),
which writes
[artifacts/generated-results/hc4_affine_pivot_collision_fibers.json](artifacts/generated-results/hc4_affine_pivot_collision_fibers.json).

## 1. Adapted affine coordinates

Choose linear coordinates \(y=(u_1,u_2,u_3,r)\) with

\[
r=\ell\mathbin{\cdot}x.
\tag{1.1}
\]

Write \(x=Py\), where the first three columns of \(P\) span
\(\ker\ell\) and \(\ell^{\mathsf T}P=e_4^{\mathsf T}\). The transformed
Hessian is

\[
H_y=P^{\mathsf T}H_xP.
\tag{1.2}
\]

The \(3\)-by-\(3\) tangential block is the Hessian of the restricted
ternary potential

\[
\psi_c(u)=\psi(P(u,c)).
\tag{1.3}
\]

The adjugate covariance identity gives

\[
\det\operatorname{Hess}_u(\psi_c)
=e_4^{\mathsf T}\operatorname{adj}(H_y)e_4
=(\det P)^2
 \ell^{\mathsf T}\operatorname{adj}(H_x)\ell.
\tag{1.4}
\]

Thus every slice \(\psi_c\) has the same nonzero constant Hessian
determinant \((\det P)^2N_\ell\).

## 2. HC3 separates every affine fiber

Suppose \(p,q\) lie on the same fiber of \(\ell\), so in adapted
coordinates

\[
p=(u_p,c),
\qquad
q=(u_q,c).
\]

Equality of the four full gradients implies equality of their first three
adapted components:

\[
\nabla_u\psi_c(u_p)=\nabla_u\psi_c(u_q).
\tag{2.1}
\]

By (1.4), \(\psi_c\) is a three-variable constant-Hessian potential.
The truth of HC3 makes its gradient injective, hence \(u_p=u_q\). The last
coordinate is already equal, so \(p=q\). This proves (0.2).

For a normalized antipodal collision at \(\pm p\), every affine covector
satisfying (0.1) must therefore obey

\[
\ell\mathbin{\cdot}p\ne0.
\tag{2.2}
\]

## 3. Consequence for Schur collision transfer

In the affine presentation of HC4RSD6,

\[
A=\ell\mathbin{\cdot}x+a_0,
\qquad
\Phi(t,x)=tA(x)+B(x).
\tag{3.1}
\]

The conclusion does not actually require the reduced pencil to be
singular. Put \(M=\operatorname{Hess}B\). For any zero-corner affine parent,

\[
\det\operatorname{Hess}\Phi
=-\ell^{\mathsf T}\operatorname{adj}(M)\ell=c\in K^\times.
\tag{3.2}
\]

A quadratic repair replaces \(M\) by
\(M+\kappa\ell\ell^{\mathsf T}\), but the bordered cofactor is invariant:

\[
\ell^{\mathsf T}
\operatorname{adj}(M+\kappa\ell\ell^{\mathsf T})\ell
=\ell^{\mathsf T}\operatorname{adj}(M)\ell=-c.
\tag{3.3}
\]

Thus every repaired affine descendant satisfies (0.1), including the
nonsingular exact-remainder branch.

If two parent gradients agree at a common pivot value \(t_0\), equality of
their pivot components and their four remaining components gives

\[
A(p)=A(q),
\qquad
\nabla B(p)=\nabla B(q),
\qquad\text{hence}\qquad
\ell\mathbin{\cdot}p=\ell\mathbin{\cdot}q.
\tag{3.4}
\]

The two displayed equalities show that every quadratic repair
\(B+\kappa A^2/2+\mu A\) has the same marked collision. But
(3.2)--(3.3) supply (0.1), and (0.2) then forces \(p=q\).
Therefore no nontrivial collision transfers through an affine scalar pivot
from five to four Hessian variables, whether its reduced Hessian is
singular or survives through the exact specialized determinant remainder.

This does not prove that every potential admitting an affine Schur
representation is injective: a collision may join different \(A\)-fibers.
It proves precisely that such a collision cannot be inherited from equal
parent gradients at one pivot value.

## 4. Reproduction and revised frontier

Run:

~~~bash
.venv/bin/python scripts/verify_hc4_affine_pivot_collision_fibers.py
~~~

The command checks the cofactor/slice determinant identity, its adapted
coordinate scaling, and the restriction of the full gradient and Hessian
to the tangential blocks. The final injectivity input is HC3.

The lower metric faces \(N_6,\ldots,N_1\) on the HC4RSD6 rank-three locus
are no longer needed for inherited affine collision transfer: the collision
hyperplane is already empty. They remain relevant only if one wants to
classify affine Schur representations that do not preserve the marked
collision.

For HC4 descent, the live coverage mechanisms are therefore:

1. nonlinear scalar pivots, whose fibers are not affine three-spaces;
2. genuinely mixed source--dual or coisotropic pivots;
3. matrix pivots with moving kernel planes;
4. nonsingular reduced pencils satisfying the specialized exact remainder
   with a nonlinear pivot;
5. or direct degree-five exclusion without an inherited higher-dimensional
   collision.

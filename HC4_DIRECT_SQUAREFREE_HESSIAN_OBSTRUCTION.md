# Squarefree top-Hessian obstruction for direct `HC4`

## Statement

Continue `HC4_DIRECT_HOMOGENEOUS_FILTRATION.md`.  Let

\[
\Psi\in K[x_1,x_2,x_3,w],\qquad
\det\operatorname{Hess}\Psi=\delta\in K^\times,
\]

with degree `D>2`, over a characteristic-zero field.  Suppose the top
homogeneous form has generic Hessian rank three.  By Gordan--Noether, after a
constant linear change,

\[
\Psi_D=f(x_1,x_2,x_3),
\qquad
A_0=\operatorname{Hess}f,
\qquad
\Delta_f=\det A_0\ne0.
\]

> **Theorem HC4-DIR2 — squarefree Hessian discriminant is impossible.**
> If `Delta_f` is squarefree in `K[x_1,x_2,x_3]`, then no such constant-Hessian
> polynomial `Psi` exists.
>
> Equivalently, every rank-three top cone of a hypothetical `HC4`
> counterexample must have a **non-squarefree ternary Hessian determinant**.

This is an all-degree obstruction and uses no auxiliary constant-Hessian
pencil.

---

## 1. The scaled matrix polynomial

Put `r=D-2` and

\[
M(t)=H_D+tH_{D-1}+\cdots+t^rH_2,
\qquad H_j=\operatorname{Hess}\Psi_j.
\]

As in `HC4-DIR1`,

\[
\det M(t)=\delta t^{4r}.
\tag{1.1}
\]

Split off the constant top-kernel direction:

\[
M(t)=
\begin{pmatrix}
A(t)&b(t)\\
b(t)^T&c(t)
\end{pmatrix},
\tag{1.2}
\]

where

\[
A(0)=A_0,
\qquad b(0)=0,
\qquad c(0)=0,
\]

and all three blocks have degree at most `r` in `t`.

The exact block determinant is

\[
\det M(t)=\det A(t)c(t)-b(t)^T\operatorname{adj}(A(t))b(t).
\tag{1.3}
\]

---

## 2. Local divisibility lemma

We need one elementary lemma.

> **Lemma.** Let `A` be a symmetric `3x3` polynomial matrix with nonzero
> squarefree determinant `Delta`.  If a polynomial column `b` satisfies
>
> \[
> \Delta\mid b^T\operatorname{adj}(A)b,
> \tag{2.1}
> \]
>
> then
>
> \[
> A^{-1}b=\frac{\operatorname{adj}(A)b}{\Delta}
> \in K[x_1,x_2,x_3]^3.
> \tag{2.2}
> \]

**Proof.**  Let `pi` be an irreducible factor of `Delta`.  Since `Delta` is
squarefree, the generic rank of `A mod pi` is exactly two.  Hence over the
fraction field of `K[x]/(pi)`,

\[
\operatorname{adj}(A)=\rho vv^T
\]

for a nonzero scalar `rho` and a kernel vector `v`.  Reducing (2.1) modulo
`pi` gives

\[
\rho(v^Tb)^2=0,
\]

so `v^Tb=0`.  Consequently

\[
\operatorname{adj}(A)b=\rho v(v^Tb)=0\pmod\pi.
\]

Thus every irreducible factor `pi` of `Delta` divides every component of
`adj(A)b`.  Squarefreeness gives (2.2).  QED.

---

## 3. Kill the first attempted rotation

Assume first that `b(t)` is nonzero, and let

\[
j\ge1
\]

be its first nonzero order:

\[
b(t)=t^jb_j+O(t^{j+1}),\qquad b_j\ne0.
\tag{3.1}
\]

Because (1.1) vanishes to order `4r`, while `j<=r`, comparison of (1.3) at
orders below `2j` recursively gives

\[
[t^m]c(t)=0\qquad(1\le m<2j).
\tag{3.2}
\]

In particular the homogeneous layer producing `b_j` is affine in `w`.  Thus

\[
b_j=\nabla a_j
\tag{3.3}
\]

for a homogeneous ternary polynomial `a_j`.

At order `2j`, (1.3) gives

\[
c_{2j}\Delta_f
=b_j^T\operatorname{adj}(A_0)b_j.
\tag{3.4}
\]

(If `2j>r`, take `c_{2j}=0`; the same divisibility statement holds.)

By the lemma,

\[
X:=A_0^{-1}b_j
\]

is polynomial.  But homogeneity gives

\[
\deg b_j=D-j-2,
\qquad
\deg A_0=D-2,
\]

so

\[
\deg X=-j<0.
\tag{3.5}
\]

The only polynomial homogeneous vector of negative degree is zero.  Hence

\[
X=0,\qquad b_j=0,
\]

contradicting the definition of `j`.

Therefore

\[
\boxed{b(t)=0.}
\tag{3.6}
\]

---

## 4. Final valuation contradiction

With `b(t)=0`, equation (1.3) becomes

\[
\det M(t)=\det A(t)c(t).
\tag{4.1}
\]

Since `det A(0)=Delta_f` is nonzero, `det A(t)` has `t`-adic valuation zero.
But `c(0)=0` and `deg_t c<=r`, so

\[
\operatorname{ord}_t\det M(t)=\operatorname{ord}_t c(t)\le r.
\tag{4.2}
\]

Equation (1.1), however, requires

\[
\operatorname{ord}_t\det M(t)=4r.
\tag{4.3}
\]

For `r>0` this is impossible.

This proves HC4-DIR2.

---

## 5. Consequence for the direct HC4 programme

A hypothetical four-variable constant-Hessian polynomial of degree larger
than two has only the following possibilities for its top homogeneous form:

1. `rank Hess Psi_D <= 2`; or
2. `rank Hess Psi_D = 3` and the ternary cone base `f` has
   **non-squarefree Hessian determinant**.

The generic rank-three/squarefree case is completely gone in every degree.

This sharply identifies the next algebraic-geometric target: classify ternary
homogeneous forms whose Hessian determinant has repeated irreducible factors,
and feed those factor multiplicities into the same `t`-adic Schur descent.

The expectation is that repeated Hessian factors are precisely where split
variables, composite polar maps, and cotangent/plane-Keller geometry can hide.

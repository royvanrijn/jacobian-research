# Split-metric obstruction in the residual `[3,1]` HC4 branch

## Status

**Current status.**  This is an intermediate proof-map note.  The rational
null-foliation target in Section 5 is closed by `HC4RSD63`, and the complete
rank-two relative-nilpotent branch is subsumed by `HC4MR1`.

Continue `HC4RSD61`.  The entire moving rank-two length-three branch has been
reduced to a bivariate active potential

\[
A=A(x,w),\qquad K=\operatorname{Hess}_{x,w}A,
\qquad \det K\ne0,
\]

and a nonzero active coupling `r` satisfying

\[
r^{\mathsf T}\operatorname{adj}(K)r=0.
\]

In passive-kernel coordinates one has

\[
r=\nabla a-\Phi_h\nabla q,
\]

with `Phi_h` nonconstant in the passive curved coordinate.

> **Theorem HC4RSD62 — split active Hessian metric.**
> In every genuine `[3,1]` packet of `HC4RSD61`, the active Hessian
> determinant is a square up to a field unit:
>
> \[
> \boxed{-\det\operatorname{Hess}A=\eta g^2}
> \tag{0.1}
> \]
>
> for some `eta in K^*` and `g in K[x,w]` after scalar extension if necessary.
> Moreover `nabla a` and `nabla q` lie on the same null characteristic line;
> hence
>
> \[
> da\wedge dq=0.
> \tag{0.2}
> \]
>
> Thus the passive-kernel motion and passive affine correction factor through
> one common one-dimensional active characteristic field.

## 1. The three null identities

`HC4RSD61` gives

\[
\begin{aligned}
(\nabla q)^{\mathsf T}\operatorname{adj}(K)\nabla q&=0,\\
(\nabla a)^{\mathsf T}\operatorname{adj}(K)\nabla q&=0,\\
(\nabla a)^{\mathsf T}\operatorname{adj}(K)\nabla a&=0.
\end{aligned}
\tag{1.1}
\]

Write

\[
K=
\begin{pmatrix}
A_{xx}&A_{xw}\\
A_{xw}&A_{ww}
\end{pmatrix}.
\]

For any nonzero null gradient `nabla q=(q_x,q_w)^T`,

\[
A_{ww}q_x^2-2A_{xw}q_xq_w+A_{xx}q_w^2=0.
\tag{1.2}
\]

## 2. Rational null slope forces square discriminant

On the chart `q_w\ne0`, put

\[
r_q=q_x/q_w\in K(x,w).
\]

Equation (1.2) says that `r_q` is a root of

\[
A_{ww}R^2-2A_{xw}R+A_{xx}=0.
\tag{2.1}
\]

Hence the discriminant

\[
4(A_{xw}^2-A_{xx}A_{ww})
=-4\det K
\]

is a square in `K(x,w)`.  The same conclusion follows on the chart `q_x\ne0`.
Therefore

\[
-\det K=\eta h^2
\tag{2.2}
\]

in the rational function field, with `eta` a field unit.

But `det K` is polynomial.  Since `K[x,w]` is a UFD, a polynomial which is a
square in its fraction field is a square in the polynomial ring up to a unit:
if `f=(p/r)^2` with coprime `p,r`, then `fr^2=p^2` forces `r` to be a unit.
Thus (2.2) may be written as (0.1) with `g` polynomial.

The null equation consequently factors over a quadratic scalar extension as

\[
\bigl(A_{ww}q_x-(A_{xw}+g)q_w\bigr)
\bigl(A_{ww}q_x-(A_{xw}-g)q_w\bigr)=0
\tag{2.3}
\]

(up to the harmless choice of the square root of `eta`).
The residual problem is therefore controlled by one of two rational
characteristic derivations.

## 3. Both moving coefficients use the same characteristic

In a two-dimensional nondegenerate symmetric space, a null line is equal to
its own orthogonal complement.  The first and third equations of (1.1) say
that `nabla q` and `nabla a` are null; the middle equation says they are
orthogonal.  Therefore, whenever both are nonzero,

\[
\nabla a\parallel\nabla q.
\tag{3.1}
\]

Equivalently

\[
J(a,q)=a_xq_w-a_wq_x=0,
\]

which is (0.2).

Thus `a` and `q` are algebraically dependent.  Over the active function field,
Lüroth's theorem supplies a common generator `h` such that

\[
K(a,q)\subseteq K(h)\subseteq K(x,w).
\tag{3.2}
\]

If `a,q` are polynomial and the generator is chosen closed, this is the usual
polynomial-composition form; the rational statement (3.2) is sufficient for
the HC4 reduction.

## 4. Geometric interpretation

The tangent vector to a level curve of `q` is

\[
X_q=(q_w,-q_x).
\]

Equation (1.2) is exactly

\[
X_q^{\mathsf T}(\operatorname{Hess}A)X_q=0.
\]

Thus the level curves of the common characteristic generator are asymptotic
curves of the Hessian metric of `A`.  HC4 requires this null foliation to have
a rational first integral **and** requires the Hessian discriminant to be a
polynomial square.  This is a much thinner locus than a general bivariate
polynomial potential.

## 5. Historical next target (closed by `HC4RSD63`)

With `g^2=-eta^{-1} det Hess A`, the null PDE is first order:

\[
A_{ww}h_x-(A_{xw}\pm g)h_w=0.
\tag{5.1}
\]

The remaining `[3,1]` problem is now the classification of polynomial Hessian
metrics whose null characteristic derivation has a nonconstant rational first
integral.  Every example found so far is triangular after a polynomial active
coordinate change; proving that globally would close `[3,1]`.

# JC2 dual-Schur curvature obstruction

## Scope

Let `K` be a characteristic-zero field and let

\[
F=(P,Q):\mathbb A^2\to\mathbb A^2,
\qquad J(P,Q)=c\in K^\times.
\]

The cotangent doubling is

\[
\Phi(x,y,u,v)=uP(x,y)+vQ(x,y),
\]

with

\[
\det\operatorname{Hess}\Phi=c^2.
\]

This note asks whether the one-variable Schur descent used in the 2026
`JC3 -> HC5` bridge can descend an arbitrary hypothetical `JC2`
counterexample to `HC3`.  Eliminating a dual variable gives an exact answer:
the only obstruction is the curvature of the target pencil.

## 1. Dual Schur matrix

Take `t=v` and `w=(x,y,u)`.  Then

\[
\Phi=tA+B,
\qquad A=Q(x,y),
\qquad B=uP(x,y).
\]

For the Schur parameter `s`, put

\[
M(s)=\operatorname{Hess}_{x,y,u}(uP+sQ).
\]

Explicitly

\[
M(s)=
\begin{pmatrix}
 uP_{xx}+sQ_{xx}&uP_{xy}+sQ_{xy}&P_x\\
 uP_{xy}+sQ_{xy}&uP_{yy}+sQ_{yy}&P_y\\
 P_x&P_y&0
\end{pmatrix}.
\]

Define

\[
\mathcal B(P)
=P_y^2P_{xx}-2P_xP_yP_{xy}+P_x^2P_{yy}
\tag{1.1}
\]

and the mixed curvature

\[
\mathcal B(P;Q)
=P_y^2Q_{xx}-2P_xP_yQ_{xy}+P_x^2Q_{yy}.
\tag{1.2}
\]

A direct determinant expansion gives

\[
\boxed{
\det M(s)=-u\mathcal B(P)-s\mathcal B(P;Q).
}
\tag{1.3}
\]

Thus the dimension-drop hypothesis `det M(s)=0` of the Schur-descent lemma
holds exactly when

\[
\mathcal B(P)=\mathcal B(P;Q)=0.
\]

More generally, after a constant target change one obtains the same statement
with `P` replaced by any nonzero target pencil member `R=aP+bQ`.

## 2. Quadratically repaired descent

For `lambda != 0` and `mu in K`, the usual Schur repair is

\[
\psi_{\lambda,\mu}
=uP+\frac\lambda2Q^2+\mu Q.
\]

Since the bordered determinant of `Phi` equals `c^2`, the rank-one update
identity gives

\[
\boxed{
\det\operatorname{Hess}\psi_{\lambda,\mu}
=-u\mathcal B(P)
-(\mu+\lambda Q)\mathcal B(P;Q)
-\lambda c^2.
}
\tag{2.1}
\]

Hence dual Schur descent produces a three-variable constant-Hessian potential
precisely when the curvature terms vanish.  Because any collision of `F`
lifts to equal gradients of `Phi` at dual coordinates `u=v=0`, such a descent
would preserve the collision and contradict `HC3`.

## 3. Zero curvature is already triangular

> **Theorem JC2DS1.**  If there is a nonzero constant pair `(a,b)` such that
>
> \[
> R=aP+bQ
> \]
>
> satisfies
>
> \[
> \mathcal B(R)=0,
> \tag{3.1}
> \]
>
> then `F` is a polynomial automorphism.

Indeed, because `JF` is invertible,

\[
\nabla R=(a,b)JF
\]

never vanishes.  Equation (3.1) is the zero-curvature equation for the smooth
level curves of `R`.  Thus every generic irreducible fiber is an affine line.
Distinct fibers cannot be nonparallel lines, so after an affine source change

\[
R=\rho(x)
\]

for a univariate polynomial `rho`.  Since `grad R` has no zero, `rho'` has no
root over the algebraic closure, hence `rho` is affine-linear.  Therefore `R`
is an affine coordinate.

After a constant target change take `R=P=x`.  The Keller equation then reads

\[
Q_y=c,
\]

so

\[
Q=cy+h(x),
\]

and `F` is triangular with an explicit polynomial inverse.

Consequently any `JC2` counterexample must satisfy the global obstruction

\[
\boxed{
\mathcal B(aP+bQ)\ne0
\quad\text{for every }[a:b]\in\mathbb P^1(K).
}
\tag{3.2}
\]

Over an algebraically closed base the same statement holds for every constant
projective target direction.

## 4. Interpretation

The homogeneous cubic

\[
C_F(a,b):=\mathcal B(aP+bQ)
\tag{4.1}
\]

is the complete second-order obstruction to descending the cotangent `HC4`
lift along a dual variable.  It packages the curvature numerators of the
entire target pencil.

Thus a hypothetical plane Keller counterexample has to satisfy two simultaneous
properties:

1. every nonzero target pencil member is a polynomial submersion; and
2. every target direction has nonzero curvature polynomial.

The first is automatic from the Keller condition.  The second is the new
Schur obstruction.

This does **not** prove `JC2`: polynomial automorphisms themselves can have
curved constant target directions.  The next target is to combine the cubic
pencil `C_F` with the flatness equations coming from `J(P,Q)=c`, or with the
Newton/Magnus degree-reduction theory, to force a polynomial target change
that lowers the curvature/degree.

# Scalar cancellation dichotomy for reverse HC4 descent

## Status

This note closes scalar exact cancellation as an independent quadratic-pivot
route from a five-variable Hessian collision to HC4. It separates the
nonzero and zero pivot corners.

> **Theorem HC4RSD11 (nonzero-corner pencil equivalence).** Let \(K\) have
> characteristic zero and
>
> \[
> \Phi(t,x)=\frac{\lambda}{2}t^2+tA(x)+B(x),
> \qquad \lambda\in K^\times,
> \quad x\in K^4.
> \tag{0.1}
> \]
>
> Put
>
> \[
> \psi(x)=B(x)-\frac{A(x)^2}{2\lambda}.
> \tag{0.2}
> \]
>
> Then
>
> \[
> \det\operatorname{Hess}_{t,x}\Phi=c
> \quad\Longleftrightarrow\quad
> \det\operatorname{Hess}_x(\psi+sA)=c/\lambda
> \quad\text{for every }s.
> \tag{0.3}
> \]
>
> Equal-gradient pairs of \(\Phi\) are in bijection with equal-gradient
> pairs of the appropriate pencil member \(\psi+sA\). Conversely, every
> constant-Hessian pencil \(\psi+sA\) has the suspension
>
> \[
> \Phi=\frac{\lambda}{2}
> \left(t+\frac{A}{\lambda}\right)^2+\psi.
> \tag{0.4}
> \]
>
> Thus nonzero-corner determinant-term cancellation is exactly the problem
> of classifying four-variable constant-Hessian pencils. Every direct HC4
> candidate embeds by taking \(A\) affine, while nonlinear \(A\) defines a
> smaller structured pencil locus.

> **Theorem HC4RSD12 (zero-corner graph-coordinate obstruction).** Let
>
> \[
> \Phi(t,x)=tA(x)+B(x)
> \tag{0.5}
> \]
>
> have nonzero constant five-variable Hessian determinant. Suppose there is
> a constant direction \(v\in K^4\) with
>
> \[
> D_vA\in K^\times.
> \tag{0.6}
> \]
>
> Then \(\nabla\Phi\) is a polynomial automorphism. In particular \(\Phi\)
> has no gradient collision, independently of whether
> \(\operatorname{Hess}(B+sA)\) is singular or survives by an exact
> specialized determinant remainder.

> Every zero-corner scalar parent with \(\deg A\le2\) satisfies (0.6).
> Consequently no five-variable Hessian collision can pass through a
> quadratic zero-corner scalar pivot.

Together with HC4RSD8--HC4RSD10, the zero-corner conclusion is complete for
quadratic scalar pivots. The singular-pencil descendants are triangular and
every zero-corner quadratic parent is injective. The remaining scalar
quadratic problem is the nonzero-corner nonlinear constant-Hessian pencil
identified by HC4RSD11.

The universal determinant identities are replayed by
[scripts/verify_hc4_scalar_cancellation_dichotomy.py](scripts/verify_hc4_scalar_cancellation_dichotomy.py),
which writes
[artifacts/generated-results/hc4_scalar_cancellation_dichotomy.json](artifacts/generated-results/hc4_scalar_cancellation_dichotomy.json).

## 1. The nonzero corner is a completed-square gauge

Put

\[
g=\nabla A,
\qquad
M(t)=\operatorname{Hess}(B+tA).
\]

The parent Hessian is

\[
\begin{pmatrix}
\lambda&g^{\mathsf T}\\
g&M(t)
\end{pmatrix}.
\tag{1.1}
\]

Fix the pivot-gradient output

\[
y=\partial_t\Phi=\lambda t+A.
\]

Set

\[
s=\frac{y}{\lambda}=t+\frac{A}{\lambda}.
\]

Then direct differentiation of (0.2) gives

\[
\nabla(\psi+sA)=\nabla B+t\nabla A,
\qquad
\operatorname{Hess}(\psi+sA)=M(t)-\lambda^{-1}gg^{\mathsf T}.
\tag{1.2}
\]

The exact block identity

\[
\det
\begin{pmatrix}
\lambda&g^{\mathsf T}\\
g&M
\end{pmatrix}
=\lambda\det(M-\lambda^{-1}gg^{\mathsf T})
\tag{1.3}
\]

proves (0.3), because \((t,x)\mapsto(s,x)\) is a polynomial automorphism.
It also proves collision equivalence: equality of parent gradients first
fixes \(y\), and the other four parent-gradient coordinates are exactly
\(\nabla(\psi+sA)\). A collision of that pencil member lifts by taking

\[
t_i=\frac{y-A(x_i)}{\lambda}.
\tag{1.4}
\]

Conversely, if the whole pencil has constant Hessian determinant, (0.4)
gives

\[
\partial_t\Phi=\lambda t+A,
\qquad
\nabla_x\Phi=\nabla\psi+
\left(t+\frac{A}{\lambda}\right)\nabla A.
\]

A collision of \(\nabla(\psi+s_0A)\) lifts on the level
\(t=s_0-A/\lambda\). Taking \(A\) affine makes the Hessian pencil constant
for every \(\psi\), so this branch contains direct HC4. For nonlinear \(A\),
the all-\(s\) determinant identity is additional structure and remains a
classification target.

## 2. Exact graph-coordinate factorization

Normalize (0.6) by constant affine coordinates and scaling:

\[
x=(u,w),
\qquad
u=(u_1,u_2,u_3),
\qquad
A=w+q(u).
\tag{2.1}
\]

Use the polynomial coordinate

\[
r=A=w+q(u)
\]

and write

\[
C(u,r)=B(u,r-q(u)).
\tag{2.2}
\]

Finally set

\[
\tau=t+C_r(u,r).
\tag{2.3}
\]

In the original variable order \((t,u,w)\), the five parent-gradient
coordinates become

\[
\boxed{
\partial_t\Phi=r,
\qquad
\nabla_u^{\rm old}\Phi
=\nabla_u C(u,r)+\tau\nabla q(u),
\qquad
\partial_w\Phi=\tau.
}
\tag{2.4}
\]

Both changes

\[
(t,u,w)\longleftrightarrow(r,u,\tau)
\]

are triangular polynomial changes with constant Jacobian up to sign. Thus
the Jacobian determinant of (2.4), equivalently the parent Hessian
determinant, is

\[
\boxed{
\det\operatorname{Hess}\Phi
=-\det\operatorname{Hess}_u(C(u,r)+\tau q(u)).
}
\tag{2.5}
\]

This can also be checked directly from the universal Hessian block. If

\[
p=\nabla q,
\quad Q=\operatorname{Hess}q,
\quad H=\operatorname{Hess}_uC,
\quad b=\nabla_uC_r,
\quad d=C_{rr},
\]

then the parent Hessian is

\[
\begin{pmatrix}
0&p^{\mathsf T}&1\\
p&H+bp^{\mathsf T}+pb^{\mathsf T}+dpp^{\mathsf T}+\tau Q&b+dp\\
1&(b+dp)^{\mathsf T}&d
\end{pmatrix},
\tag{2.6}
\]

whose determinant is exactly \(-\det(H+\tau Q)\).

## 3. HC3 makes every graph-coordinate parent injective

If the parent determinant is \(c\in K^\times\), (2.5) says that for every
fixed \(r,\tau\), the ternary potential

\[
C_{r,\tau}(u)=C(u,r)+\tau q(u)
\tag{3.1}
\]

has Hessian determinant \(-c\). By HC3 its gradient is injective over the
algebraic closure.

Suppose two parent gradients agree. Equation (2.4) first gives equal \(r\)
and equal \(\tau\). The middle three equations are then equal gradients of
the same \(C_{r,\tau}\), so HC3 gives equal \(u\). Finally

\[
w=r-q(u),
\qquad
t=\tau-C_r(u,r)
\tag{3.2}
\]

recover the remaining variables. Thus the full five-variable gradient is
injective. Ax--Grothendieck makes it a polynomial automorphism, and faithful
flatness descends the conclusion from the algebraic closure to \(K\).

## 4. Why every quadratic zero-corner pivot is covered

Let

\[
A=\frac12x^{\mathsf T}Qx+a^{\mathsf T}x+a_0,
\qquad Q=Q^{\mathsf T}.
\tag{4.1}
\]

For the zero-corner parent, with \(M(t)=\operatorname{Hess}(B+tA)\),

\[
c=-\nabla A^{\mathsf T}\operatorname{adj}(M(t))\nabla A.
\tag{4.2}
\]

The right side writes the nonzero constant \(c\) as a polynomial
combination of the entries of \(\nabla A=Qx+a\). Hence these affine entries
generate the unit ideal, so \(Qx+a=0\) has no solution over the algebraic
closure. Therefore \(Q\) is singular and \(a\notin\operatorname{im}Q\).
Since \(Q\) is symmetric, there is \(v\in\ker Q\) with

\[
D_vA=v^{\mathsf T}a\ne0.
\tag{4.3}
\]

This is precisely (0.6). No singularity, corank, or exact-remainder
assumption on the reduced pencil was used.

## 5. Revised scalar frontier

The zero-corner quadratic scalar-pivot mechanism is exhausted:

- singular zero-corner descendants are explicitly triangular by
  HC4RSD8--HC4RSD10;
- every zero-corner quadratic parent is collision-free by HC4RSD12, even
  when its reduced pencil is nonsingular;
- every nonzero-corner cancellation is a four-variable constant-Hessian
  pencil by HC4RSD11. Its affine directions are exactly direct HC4, while
  its nonlinear directions remain structured targets.

The live reverse-descent mechanisms are therefore nonlinear
constant-Hessian pencils in the nonzero-corner scalar branch, moving
matrix-pivot planes, and genuinely mixed/coisotropic canonical
transformations. Direct degree-five HC4 classification remains the parallel
route.

## 6. Reproduction

Run:

~~~bash
.venv/bin/python scripts/verify_hc4_scalar_cancellation_dichotomy.py
~~~

The checker verifies the universal \(5\)-by-\(5\) block determinant, the
graph-coordinate Hessian factorization, exact gradient coordinates, a
degree-unbounded nonlinear calibration, and the nonzero-corner suspension
identity.

# Direct HC4 attack by homogeneous filtration

## Purpose

The relative-nilpotent programme isolates a large class of constant-Hessian
pencils, but an arbitrary `HC4` polynomial need not lie on an affine line in
the constant-Hessian locus.  This note starts a direct attack that applies to
**every** four-variable constant-Hessian polynomial.

It also records why the obvious cotangent use of the relative-nilpotent theorem
cannot by itself prove `JC2`.

Throughout let `K` have characteristic zero and let

\[
\Psi\in K[x_1,x_2,x_3,w],\qquad
\det\operatorname{Hess}\Psi=\delta\in K^\times.
\]

Write

\[
\Psi=\Psi_D+\Psi_{D-1}+\cdots+\Psi_2+\text{affine}
\]

for its ordinary homogeneous decomposition, with `D>2`.

---

## 1. Why the direct `JC2` cotangent-pencil attempt is tautological

Let

\[
F=(P,Q):K^2\to K^2,\qquad J(P,Q)=1,
\]

and form its cotangent potential

\[
\Phi(x,y,u,v)=uP(x,y)+vQ(x,y).
\tag{1.1}
\]

For every constant nilpotent matrix `C in Mat_2(K)`, put

\[
A_C=(u,v)\,C\,F(x,y)^T.
\tag{1.2}
\]

Then

\[
\Phi+sA_C=(u,v)(I+sC)F^T
\]

and the standard cotangent determinant identity gives

\[
\det\operatorname{Hess}(\Phi+sA_C)
=\det(I+sC)^2 J(P,Q)^2=1.
\tag{1.3}
\]

Thus **every** plane Keller map carries a whole family of exact
relative-nilpotent HC4 pencils.

If `C` has rank one, `A_C` is a scalar fiber variable times one linear
combination of `P,Q`; consequently `Hess A_C` has a constant fiber-kernel
line.  The fixed-kernel theorem `HC4RSD65` therefore applies, but its terminal
`JC2` cotangent endpoint is exactly (1.1).  Nothing in that argument can
separate an automorphic Keller pair from a hypothetical plane counterexample.

This explains the failure of the previous direct attempt: the construction
embeds `JC2` into the branch that the master theorem deliberately leaves as
its exact endpoint.

---

## 2. Scaling every HC4 polynomial gives a pure-determinant degeneration

Set

\[
r=D-2,
\qquad
M(t)=H_D+tH_{D-1}+\cdots+t^rH_2,
\qquad
H_j=\operatorname{Hess}\Psi_j.
\tag{2.1}
\]

Since

\[
\operatorname{Hess}\Psi(\lambda x)
=\lambda^{D-2}M(\lambda^{-1}),
\]

the constant determinant condition is equivalent to

\[
\boxed{\det M(t)=\delta\,t^{4r}.}
\tag{2.2}
\]

This is stronger than merely saying `det H_D=0`: the complete matrix
polynomial has a determinant with **one root only**, at `t=0`, of multiplicity
`4r`.

Equation (2.2) is the direct replacement for the affine-pencil identity used
in the relative-nilpotent programme.

---

## 3. The top homogeneous Hessian has a constant kernel

The leading coefficient of (2.2) gives

\[
\det H_D=0.
\]

In four variables the homogeneous Gordan--Noether theorem applies.  Hence,
after one constant linear change, the top homogeneous form is a cone.

On the generic rank-three branch we may therefore write

\[
\Psi_D=f(x_1,x_2,x_3),
\qquad
A=\operatorname{Hess}_{x_1,x_2,x_3}f,
\qquad
\det A\ne0
\tag{3.1}
\]

in the fraction field.  Thus

\[
H_D=\begin{pmatrix}A&0\\0&0\end{pmatrix}
\tag{3.2}
\]

with a **constant** top kernel `K e_w`.

This is the crucial advantage over the final moving-kernel problem of the
relative-nilpotent analysis: the degree filtration starts with an affine-fixed
kernel automatically.

---

## 4. First descent equation

Write

\[
\Psi_{D-1}=w\,a(x_1,x_2,x_3)+b(x_1,x_2,x_3)
\tag{4.1}
\]

for the moment.  This form is in fact forced by the determinant.

Indeed, block `M(t)` according to `(x_1,x_2,x_3)|w`:

\[
M(t)=
\begin{pmatrix}
A+tA_1+\cdots & t\,\nabla a+O(t^2)\\
 t(\nabla a)^T+O(t^2) & t c_1+t^2c_2+\cdots
\end{pmatrix}.
\]

The coefficient of `t` in the determinant is

\[
[t]\det M(t)=c_1\det A.
\tag{4.2}
\]

But (2.2) has no term of order `t`, hence `c_1=0`.  Since

\[
c_1=\partial_w^2\Psi_{D-1},
\]

we obtain (4.1).

So the top affine direction propagates through the **first** lower homogeneous
layer for every HC4 polynomial.

---

## 5. Second descent equation: the ternary Schur obstruction

Write the quadratic-in-`w` part of the next layer as

\[
\Psi_{D-2}=\frac{w^2}{2}q(x_1,x_2,x_3)
+w\,c(x_1,x_2,x_3)+d(x_1,x_2,x_3).
\tag{5.1}
\]

The coefficient of `t^2` in the same block determinant, using `c_1=0`, is

\[
[t^2]\det M(t)
=q\det A-(\nabla a)^T\operatorname{adj}(A)\nabla a.
\tag{5.2}
\]

Again (2.2) has no such coefficient.  Therefore every rank-three top cone in
an arbitrary HC4 polynomial satisfies

\[
\boxed{
q\det\operatorname{Hess}f
=(\nabla a)^T
\operatorname{adj}(\operatorname{Hess}f)
\nabla a.
}
\tag{5.3}
\]

The degrees are rigid:

\[
\deg f=D,
\qquad
\deg a=D-2,
\qquad
\deg q=D-4.
\tag{5.4}
\]

Equation (5.3) is a **ternary homogeneous reverse-Schur equation**.  It is the
first genuinely new obstruction in the direct HC4 attack: it does not assume
that the original polynomial lies in any auxiliary constant-Hessian pencil.

Equivalently,

\[
(\nabla a)^T(\operatorname{Hess}f)^{-1}\nabla a=q
\tag{5.5}
\]

in the fraction field, and the left side is forced to be polynomial.

---

## 6. What the new problem really is

A naive propagation claim `a=0` is false.  For example, with

\[
f=x^4+y^4+z^4,
\]

there are nonzero quadratic `a` for which the rational expression in (5.5)
is polynomial (e.g. coordinate-square directions).  These are precisely the
kind of split directions that should reduce the four-variable geometry rather
than contradict it.

The correct all-degree target is therefore the following dichotomy.

> **Ternary Schur dichotomy (target).**  Let `f` be a homogeneous ternary form
> with generically invertible Hessian.  If homogeneous `a,q` of degrees
> `(D-2,D-4)` satisfy (5.3), then either
>
> 1. `f,a` share a constant split direction, allowing degree/variable
>    reduction; or
> 2. `(f,a,q)` has a cotangent/plane-Keller normal form.

If this dichotomy holds and is stable under the subsequent homogeneous layers,
then the degree filtration proves the global structural reduction

\[
\boxed{HC4\Longrightarrow HC2\text{ or the exact }JC2\text{ cotangent case}}
\]

without any relative-pencil hypothesis.

The remaining obstruction would then be literally `JC2`, not a four-variable
moving-flag problem.

---

## 7. Immediate algebraic interpretation

On the Hessian discriminant

\[
\Delta_f=\det\operatorname{Hess}f=0,
\]

a generic point has Hessian rank two and

\[
\operatorname{adj}(\operatorname{Hess}f)=\rho\,vv^T
\]

for its one-dimensional Hessian kernel `v`.  Equation (5.3) implies

\[
(v\cdot\nabla a)^2=0
\qquad\text{on }\Delta_f,
\]

hence

\[
\boxed{v\cdot\nabla a=0\quad\text{on the Hessian discriminant}.}
\tag{7.1}
\]

Thus `a` is constant along the kernel characteristic of the polar/gradient map
of `f` on its ramification divisor.  This is the geometric form of the new
obstruction and suggests attacking (5.3) through the geometry/factorization of
the ternary Hessian discriminant rather than by total degree.

---

## 8. Verification

The companion checker verifies the two universal determinant coefficients in
Sections 4--5 with a completely generic symmetric `3x3` leading block.

Run

```bash
.venv/bin/python scripts/verify_hc4_direct_homogeneous_filtration.py
```

The next research task is to classify the homogeneous solutions of (5.3) by
the factorization type of the Hessian discriminant, starting with squarefree
`Delta_f` and then the repeated-component cases.  This is the direct HC4
attack that was missing from the relative-nilpotent programme.

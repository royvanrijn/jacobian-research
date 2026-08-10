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

`HC4-DIR2` now excludes squarefree `Delta_f` in every degree.  `HC4-DIR3`
closes the first repeated-component stratum

\[
\Delta_f=\ell^2R,
\qquad R\text{ squarefree},
\qquad \gcd(\ell,R)=1,
\]

for every `D>=5`: radical divisibility forces a constant vector, and the two
resulting normal forms have `ell`-multiplicity `m=D-2` or at least `2m`, never
two.  `HC4-DIR3a` applies the same gate to exact multiplicity three: it is
impossible for `D>=6`, while `D=5` has only the additive split top
`f=c*ell^5+h(y,z)`.  `HC4-DIR3b` closes that terminal packet: the first
weighted channel fixes its `w^2*ell` coefficient, and the next is the nonzero
immutable square `-4*B^2*det Hess_(y,z)(h_5)`.  Thus exact linear
multiplicity three is impossible for every `D>=5`.

`HC4-DIR4` closes the generic rank-at-most-one boundary of

\[
\Delta_f=\ell^4R,
\qquad R\text{ squarefree},
\qquad \gcd(\ell,R)=1,
\]

for every `D>=5`.  `HC4-DIR4a` then eliminates the apparent order-two
degree-six split by a quadratic-suspension descent.  The exact rank-two
order-one linear field initially satisfies

\[
\operatorname{Hess}(f)L=\ell^2\nabla a,
\qquad L(a)=\ell^2c.
\]

`HC4-DIR5` proves that its matrix has rank one and reduces the entire exact
quadruple-linear frontier to

\[
D=6,
\qquad f=C\ell^6+h_6(y,z),
\qquad L=\ell\partial_\ell,
\]

up to constant coordinates and rescaling.

`HC4-DIR6` closes this last packet.  The order-one Schur coefficient fixes
the `w^2*ell^2` term; the next passive-top determinant channel fixes the
`w^3` coefficient; and the following channel is the nonzero immutable term

\[
\frac{128}{3375}\frac{\alpha^5}{C^3}
\det\operatorname{Hess}_{y,z}(h_6).
\]

Thus exact linear multiplicity four is impossible for every `D>=5`.

`HC4-DIR7` also closes exact linear multiplicity five on the generic
corank-one boundary.  The multiplicity budget gives `j<=2`; both orders force
the degree-seven split top `f=C*ell^7+h_7(y,z)`.  The order-two weighted
channel dies immediately, while the order-one ladder terminates in the
nonzero square `-9*gamma^2*det Hess_(y,z)(h_7)`.

`HC4-DIR8` closes the complementary lower-rank boundary.  Its two rank-one
jets would need a nonzero second derivative in the first unused transverse
coefficient.  The order-one field equations instead make that coefficient
linear in the transverse variable, or force `m=3` where it has degree one.
Thus exact linear multiplicity five is impossible in every boundary rank.

`HC4-DIR9` begins exact multiplicity six on the generic corank-one boundary.
The half-radical budget gives `j<=3`.  The orders `j=3` and `j=2` both reduce
to `D=8`, `f=C*ell^8+h_8(y,z)` and end respectively in the nonzero channel
multipliers `(24/7)*(alpha^3/C)` and
`-(1875/3136)*(alpha^4/C^2)`.  Only `j=1` survives there, as the exact
quadratic-vector system

\[
\operatorname{Hess}(f)Q=\ell^3\nabla a,
\qquad \ell^3\mid Q(a),
\qquad \deg Q=2.
\]

With `N=Jac(Q)` this is an order-one linear-matrix system, and its boundary
matrix satisfies `rank(N mod ell)<=2`.  `HC4-DIR11` closes rank zero: Hessian
integrability forces `Q=ell^2*u`, and the resulting degree-eight split ladder
ends in the nonzero multiplier
`2187*alpha^7/(4302592*C^5)`.  `HC4-DIR12` then treats rank one.  The axial
case reduces to the earlier linear-field system and the normal image case has
an immediate valuation conflict.  Every survivor is tangent and has
`t=ord_ell(partial_z f)` in `{3,4,5,6}`.  For `t<6` its boundary quadratic is
a square; at `t=6` it is arbitrary, with odd `m` required in the nonsquare
case.  Thus rank one is reduced to four explicit tangent packets, while rank
two remains.  `HC4-DIR13` globally eliminates the auxiliary field from those
rows: `q/ell^2` and `F/ell^m` are functionally dependent.  The quadratic
centralizer split leaves a binary-composite pencil
in `ell` and one linear form, and four primitive conic packets.  Below `t=6`
the composite form is `f=z*G(ell,y)+h(ell,y)`; at `t=6` there is also a
transverse composite orientation.  The primitive parities are
`(t,m)=(3,even),(4,odd),(5,even),(6,odd)`; for `t<6` their quadratic is
`q=y^2+ell*z`.  `HC4-DIR14` closes the invariant composite rows
`t=3,5,6` by the exact binary Hessian valuation and the first two boundary
coefficients.  Only its `t=4` pencil survives, with
`ord_ell(G)=4` and nonzero pure-power boundary value for the binary remainder.
`HC4-DIR15` closes that last pencil: the binary part of the constant field
must annihilate `G`; its zero and nonzero cases respectively contradict
`q mod ell!=0` and the required boundary derivative.  Hence the invariant
composite orientation is empty.
`HC4-DIR16` closes the transverse `t=6` composite orientation as well.  Its
active constant direction would lie in the kernel of the invertible boundary
binary Hessian; after it vanishes, the same impossible logarithmic identity
as in `HC4-DIR15` remains.  Thus all composite rank-one packets are empty.
`HC4-DIR17` closes the primitive rows.  For `t>=4`, both constant vectors lie
in the one-dimensional boundary-Hessian kernel, reducing again to the
impossible pure-quadratic identity.  The sole primitive `t=3` row has
`q=y^2+ell*z` and a nonzero `ell^1` coefficient in `Q(F)`, contradicting
`ell^3|Q(F)`.  Therefore boundary Jacobian rank one is empty, and the generic
sextuple system has boundary Jacobian rank exactly two.
The kernel-locking step is degree-free: on any corank-one Hessian boundary,
two constant directions with directional derivatives in `(ell^2)` must be
parallel.  It can therefore be reused at higher linear multiplicities.
`HC4-DIR18` reduces the remaining rank-two system using the primitive left
kernel of its linear boundary matrix.  That generator has degree at most two,
so either the top vanishes on `ell=0` and the normal component of `Q` lies in
`(ell^2)`, or the boundary binary form has at most three distinct roots.
Over the algebraic closure the latter are the pure-power, two-root monomial,
and three-root `y^a*z^b*(y-z)^c` profiles.  This root-count lemma applies to
any linear rank-two boundary matrix, independently of sextuple multiplicity.
`HC4-DIR19--20` treat the pure-power boundary: exact order first collapses
its second jet to `C*ell^2*y^m`, after which the scalar and matrix field
equations give incompatible coefficients.  Hence the one-root profile is
empty.  `HC4-DIR21` reduces the remaining first jets to two outer monomials
in the two-root case and one scalar class in the three-root case.
`HC4-DIR22` applies the boundary Schur complement.  The three-root scalar
class is nonpolynomial unless it vanishes, while a zero first jet pushes the
first `ell`-dependent term to order eight.  The only other survivors are two
explicit outer second-jet families, requiring the shifted root exponent to
be at least four.

`HC4-DIR10` reduces the lower-rank sextuple boundary to one degree-five
tangent packet at `j=2` and two explicit order-six rank-one jets at `j=1`.

See `HC4_DIRECT_DOUBLE_LINEAR_HESSIAN_GATE.md`.

The next research task is therefore the normal rank-two packet, the two
outer second-jet families, and the order-eight two-/three-root packets,
together with the three lower-rank Hessian packets.
After those come nonlinear repeated factors, linear
multiplicity at least seven, or several distinct repeated factors.
Rank-at-most-two top Hessians remain a separate synchronization problem.  This
is the direct HC4 attack that was missing from the relative-nilpotent
programme.

# HC4 relative-nilpotent master reduction

## Purpose

This note replaces the incremental `HC4RSD17--76` narrative by one proof tree.
It records the strongest statement currently proved for the relative-nilpotent
Hessian-pencil branch and, crucially, isolates the **single remaining global
algebraic bridge**.  It deliberately does not identify local smooth
triangularization with a constant affine flag.

Throughout let `K` be a characteristic-zero field and let

\[
S=\operatorname{Hess}\psi,\qquad
T=\operatorname{Hess}A,\qquad
\det S=\delta\in K^\times,
\]

with

\[
\det(S+sT)=\delta\qquad\text{for all }s.
\tag{0.1}
\]

Put

\[
N=S^{-1}T.
\tag{0.2}
\]

Then `N` is polynomial, nilpotent, `S`-self-adjoint and Hessian-integrable:

\[
N^{\mathsf T}S=SN,\qquad SN=T=\operatorname{Hess}A.
\tag{0.3}
\]

---

## Master theorem

> **Theorem HC4-MR — relative-nilpotent master reduction.**
> Under (0.1), every generic Jordan stratum except one reduces globally to
> either an `HC2` packet or the exact cotangent lift of a plane Keller map.
> The sole residual stratum is the regular rank-three block
> \[
> \operatorname{rank}T=3,\qquad N\sim J_4(0),
> \tag{0.4}
> \]
> with a linearly-independent associated four-variable quasi-translation.
> On this final stratum the complete Jordan flag
> \[
> \ker N\subset\ker N^2\subset\ker N^3
> \tag{0.5}
> \]
> is Frobenius-integrable on the generic constant-rank locus.  Hence `N` is
> locally strictly upper triangularizable there.
>
> What remains is purely global: prove that this polynomial integrable flag
> either admits a constant affine invariant or produces a global plane Keller
> quotient.  No degree-by-degree scalar obstruction remains, and no local
> focal/Frobenius modulus remains.

The theorem is proved in six conceptual steps.

---

## 1. Nilpotent-pencil equivalence

Equation (0.1) is equivalent to

\[
\det(I+sN)=1.
\]

Thus all coefficients of the characteristic polynomial of `N` vanish and

\[
N^4=0.
\tag{1.1}
\]

Conversely nilpotence gives (0.1).  This turns the scalar constant-Hessian
pencil into a finite Jordan-rank problem.

The inverse pencil is automatically finite:

\[
(S+sT)^{-1}
=(I+sN)^{-1}S^{-1}
=(I-sN+s^2N^2-s^3N^3)S^{-1}.
\tag{1.2}
\]

This polynomial inverse-Hessian identity is the common algebraic source of
all cofactor/Krylov constraints below.

---

## 2. Scalar reverse-Schur packets disappear globally

For a residual scalar packet

\[
\Psi(x,y,z,w)=w\,c(x,y,z)+D(x,y,z),
\]

constant nonzero Hessian determinant gives simultaneously

\[
V(c_x,c_y,c_z)=\varnothing
\tag{2.1}
\]

and the three-variable bordered-Hessian equation

\[
(\nabla c)^{\mathsf T}
\operatorname{adj}(\operatorname{Hess}c)\nabla c=0.
\tag{2.2}
\]

Choose one irreducible smooth fiber `c=lambda`.  Equation (2.2) makes its
projective closure developable.  A smooth affine nonplanar tangent developable
would have its edge of regression at infinity, forcing the whole tangent
surface to lie at infinity.  Hence the projective surface is a cone with
vertex at infinity.  The affine fiber is therefore a cylinder.  Since one
irreducible fiber is invariant under a constant translation direction `v`,

\[
c-\lambda\mid D_vc,
\]

and degree forces

\[
D_vc=0.
\tag{2.3}
\]

Therefore the complete scalar reverse-Schur branch has a fixed ruling in every
degree and reduces to `HC2` or the exact `JC2` cotangent endpoint.

The old degree-eight, degree-nine and all-`h=0` calculations remain useful as
independent elementary certificates, but are no longer logically necessary
for this branch.

---

## 3. Rank at most two closes in all degrees

### Rank one

The Hessian direction depends on one affine coordinate.  Hessian
integrability reduces the packet to either an `HC2` suspension or

\[
\psi=zP(x,y)+wQ(x,y)+R(x,y),
\qquad
\det J(P,Q)\in K^\times,
\tag{3.1}
\]

which is the exact cotangent lift of a plane Keller map.

### Rank two: Jordan type `[2,2]`

A constant image apex for the rank-two Hessian produces a primitive kernel
quasi-translation.  The HC4 metric supplies a slice-like normalization and the
cofactor identity

\[
\operatorname{adj}(H)=-\delta\,\bar k\bar k^{\mathsf T}
\tag{3.2}
\]

for the ternary passive Hessian `H`.  The moving ternary kernel either
synchronizes to a fixed direction or has the exceptional quasi-translation
normal form; Hessian integrability sends the latter directly to

\[
\psi=zP(x,w)+yQ(x,w)+R(x,w),
\qquad
\det\operatorname{Hess}\psi=\det J(P,Q)^2.
\tag{3.3}
\]

Thus `[2,2]` is `HC2/JC2`.

### Rank two: Jordan type `[3,1]`

The rank-two apex reduces the nonlinear pencil direction to a fixed active
plane.  Writing the remaining passive geometry as

\[
\psi=\Phi(x,w,h)+a(x,w)z,
\qquad h=y-q(x,w)z,
\]

the full determinant supplies the binary bordered-Hessian equation for every
member of the pencil `a-lambda q`.  Polynomial zero-curvature level curves are
parallel affine lines.  Hence `a,q` share one affine characteristic, and the
complete determinant factorizes as

\[
\det\operatorname{Hess}\psi
=-\bigl(a'-q'\Phi_h\bigr)^2
  \det\operatorname{Hess}_{w,h}\Phi.
\tag{3.4}
\]

A nonzero constant determinant forces `q'=0`; the packet is an `HC2`
suspension.

Therefore every rank-at-most-two relative-nilpotent packet is globally closed.

---

## 4. Rank three: all linearly-dependent kernel packets close

Now suppose

\[
\operatorname{rank}T=3.
\]

The associated singular-Hessian construction gives a primitive polynomial
kernel field

\[
D=k\cdot\nabla,
\qquad
Tk=0,
\qquad
Dk=0,
\tag{4.1}
\]

so `x+t k(x)` is a quasi-translation.

If the components of `k` are linearly dependent, the four-variable
singular-Hessian classification reduces `A`, after a constant affine change,
to a degenerate cotangent form

\[
A=yP(x,w)+zQ(x,w)+R(x,w),
\qquad J(P,Q)=0.
\tag{4.2}
\]

Write

\[
P=p(h),\qquad Q=q(h)
\]

through a closed/generative polynomial `h(x,w)` and put

\[
\Theta=p'q''-q'p''.
\tag{4.3}
\]

The HC4 determinant reduction gives a polynomial identity

\[
(\tau L)(\Theta L+B)=c\in K^\times,
\qquad
\tau=h_x\partial_w-h_w\partial_x.
\tag{4.4}
\]

Both factors in (4.4) are units.  Therefore

\[
\tau L=\alpha\in K^\times,
\qquad
\Theta L+B=\beta\in K^\times.
\]

Applying `tau` to the second identity gives

\[
0=\Theta\tau L=\Theta\alpha,
\]

hence

\[
\boxed{\Theta=0.}
\tag{4.5}
\]

Thus the projective kernel direction is constant.  The supposedly moving
linearly-dependent rank-three branch is impossible; the residue is the
already-closed fixed-direction `HC2/JC2` geometry.

---

## 5. The linearly-independent `[4]` block has no local geometric modulus

The only remaining generic stratum is therefore

\[
\operatorname{rank}T=3,
\qquad N\sim J_4(0),
\]

with the four components of the associated kernel quasi-translation linearly
independent.

Let

\[
F=\nabla A
\]

and let `g` generate the prime relation of the three-dimensional gradient
image

\[
Y=\overline{F(\mathbb A^4)}\subset\mathbb A^4.
\]

Put

\[
n=\nabla g(F).
\tag{5.1}
\]

Then `n` spans `ker T`, is a quasi-translation on the generic locus, and

\[
J_n=(\operatorname{Hess}g)(F)\,T
\tag{5.2}
\]

is nilpotent.

The gradient-image hypersurface is developable.  The HC4 companion gradient
maps every kernel orbit to a cubic-or-lower polynomial curve contained in one
affine tangent hyperplane, and the quotient determinant produces a cyclic
three-dimensional Krylov chain.

The last local obstruction is the Frobenius condition for `ker N^3`.  Choose
an `S`-adapted Jordan chain

\[
Ne_1=0,\quad Ne_2=e_1,\quad Ne_3=e_2,\quad Ne_4=e_3,
\]

with `e1=n` and anti-diagonal `S`.  Let

\[
B=(\operatorname{Hess}g)(F).
\]

Developability makes the Gauss direction `S e1` radical for the second
fundamental form, so

\[
b_{24}=b_{34}=b_{44}=0.
\tag{5.3}
\]

Direct multiplication gives

\[
\operatorname{tr}(BT)=2b_{24}+b_{33}.
\tag{5.4}
\]

Since `J_n=BT` is nilpotent,

\[
\operatorname{tr}J_n=0,
\]

and therefore

\[
\boxed{b_{33}=0.}
\tag{5.5}
\]

For

\[
\lambda=S e_1,
\]

self-adjointness gives

\[
\ker N^3=\ker\lambda.
\]

Hessian symmetry identifies the restriction of `d lambda` to `ker lambda`
with the three coefficients

\[
b_{44},\qquad b_{34},\qquad b_{33}-b_{24}.
\]

Equations (5.3)--(5.5) kill all three, hence

\[
\boxed{\lambda\wedge d\lambda=0.}
\tag{5.6}
\]

Thus `ker N^3` is Frobenius.  The distributions `ker N` and `ker N^2` are
already integrable, so the complete Jordan flag is integrable.  The regular
nilpotent operator is therefore locally strictly upper triangularizable on the
generic constant-rank locus.

This is the genuine dimension-four collapse: after the Gauss-radical equations
there is exactly one remaining Frobenius scalar, and nilpotent trace kills it.
In dimension five several independent scalars remain.

---

## 6. Exact globalization boundary

The preceding section closes the **local differential-geometric** problem.  It
does *not* by itself produce a constant affine flag.

This distinction is essential because ordinary Hessians are tied to the
ambient flat affine coordinates.  A nonlinear coordinate change that
triangularizes `N` does not preserve the simple form

\[
S=\operatorname{Hess}\psi,
\qquad
T=\operatorname{Hess}A.
\]

An exact first-order audit in an `S`-adapted regular-nilpotent frame imposes:

1. Hessian/Codazzi symmetry of `S`;
2. Hessian/Codazzi symmetry of `T`;
3. Frobenius of `ker N^2` and `ker N^3`;
4. quasi-translation normalization `nabla_{e1}e1=0`;
5. constant affine volume from `det S in K^*`.

These equations still permit exactly three transverse projective motions of
the kernel line:

\[
\Gamma^2_{31},\qquad
\Gamma^2_{41},\qquad
\Gamma^3_{41}.
\tag{6.1}
\]

All lower-triangular motions vanish.  Thus the surviving motion is itself
strictly upper triangular, matching the quotient/Krylov geometry, but it need
not be affine-parallel at first order.

Run

```bash
.venv/bin/python scripts/verify_hc4_affine_bridge_first_order.py
```

for the exact `64`-unknown linear audit.

Accordingly the remaining statement is now a single global algebraic problem.

> **Global affine-or-Keller bridge.**  Let a polynomial regular `[4]`
> relative-nilpotent Hessian pencil satisfy (0.1) and suppose its complete
> Jordan flag is Frobenius on the generic locus.  Prove that either
>
> 1. the kernel flag has a constant affine invariant, reducing to the fixed
>    `HC2/JC2` packets; or
> 2. the upper-triangular global twisting descends to a polynomial plane
>    Keller quotient.

This bridge cannot be replaced by the stronger assertion “Frobenius implies
constant affine flag”: the first-order audit disproves that implication at the
formal level.

---

## Consequence for `JC2`

Every plane Keller map

\[
F=(P,Q),\qquad J(P,Q)=c\in K^\times,
\]

has the cotangent potential

\[
\Psi(x,y,z,w)=zP(x,y)+wQ(x,y),
\]

with

\[
\det\operatorname{Hess}\Psi=c^2.
\tag{7.1}
\]

Thus

\[
HC4\Longrightarrow JC2
\]

by the standard cotangent bridge.

The reductions above show considerably more inside the relative-nilpotent
branch: every classified stratum already runs *back* to `HC2` or this exact
plane-cotangent geometry.  If the Global affine-or-Keller bridge is proved,
then the complete relative-nilpotent branch is equivalent, in obstruction
content, to `JC2`.

As of August 2026 this is particularly meaningful: after the dimension-three
Jacobian counterexample and the five-variable Hessian counterexample, the only
open classical low-dimensional statements are precisely `JC2` and `HC4`.

---

## Proof architecture in one line

\[
\boxed{
\text{constant Hessian pencil}
\Rightarrow
\text{self-adjoint nilpotent }N
\Rightarrow
\begin{cases}
\operatorname{rank}T\le2 &\to HC2/JC2,\\
\operatorname{rank}T=3,\ k\text{ dependent} &\to HC2/JC2,\\
\operatorname{rank}T=3,\ k\text{ independent}
&\to\text{developable image}\to\text{Frobenius Jordan flag}\to\boxed{\text{global bridge}}.
\end{cases}}
\]

Everything before the final boxed bridge is degree-free.

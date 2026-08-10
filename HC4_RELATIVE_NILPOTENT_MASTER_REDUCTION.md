# HC4 relative-nilpotent master reduction

## Purpose

This note replaces the incremental `HC4RSD17--80` narrative by one proof tree.
It records the strongest statement currently proved for the relative-nilpotent
Hessian-pencil branch.  The final regular `[4]` globalization is closed by the
affine-plane/flatness calculation and the affine-hyperplane pencil argument;
it does not identify local smooth triangularization with a constant affine
flag.

The consolidated result is registered as `HC4MR1` in `MATH_STATUS.json`.
The labels `HC4RSD41--80` below are local proof-map identifiers rather than
forty additional headline status entries.

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

## Canonical proof map

The master statement consolidates the following active proof notes.  Earlier
parallel routes with conflicting local identifiers are preserved only in the
[HC4 archive](archive/hc4-superseded-branches/README.md).

| Stage | Canonical proof notes |
|---|---|
| Scalar partitions and all-degree tail closure | [degree eight](HC4_DEGREE_EIGHT_CONSTANT_JORDAN_CLOSURE.md), [degree nine](HC4_DEGREE_NINE_SCALAR_CLOSURE.md), [general scalar patterns](HC4_GENERAL_SCALAR_PATTERNS.md), [tail descent](HC4_MINIMAL_EXCESS_TAIL_DESCENT.md), [two-tail closure](HC4_MINIMAL_EXCESS_TWO_TAIL_CLOSURE.md), [minimal-excess completion](HC4_MINIMAL_EXCESS_COMPLETE.md), and [global developable obstruction](HC4_GLOBAL_SMOOTH_DEVELOPABLE_OBSTRUCTION.md) |
| Moving nilpotent frames | [constant-metric closure](HC4_CONSTANT_METRIC_MOVING_NILPOTENT_CLOSURE.md) and [homogeneous cone-pencil classification](HC4_HOMOGENEOUS_CONE_PENCIL_CLASSIFICATION.md) |
| Rank-two Jordan packets | [apex reduction](HC4_RANK_TWO_APEX_REDUCTION.md), [quasi-translation](HC4_RANK_TWO_QUASITRANSLATION.md), [kernel synchronization](HC4_RANK_TWO_KERNEL_SYNCHRONIZATION.md), [complete square-zero closure](HC4_RANK_TWO_COMPLETE_CLOSURE.md), [length-three reduction](HC4_RANK_TWO_LENGTH_THREE_REDUCTION.md), [split metric](HC4_RANK_TWO_LENGTH_THREE_SPLIT_METRIC.md), and [length-three closure](HC4_RANK_TWO_LENGTH_THREE_COMPLETE_CLOSURE.md) |
| Rank-three reduction to the plane endpoint | [cofactor flag](HC4_RANK_THREE_COFACTOR_FLAG.md), [fixed-kernel closure](HC4_RANK_THREE_FIXED_KERNEL_CLOSURE.md), [null-fiber reduction](HC4_RANK_THREE_NULL_FIBER_REDUCTION.md), [one-active-direction obstruction](HC4_RANK_THREE_ONE_ACTIVE_DIRECTION_OBSTRUCTION.md), [affine-generator Wronskian](HC4_RANK_THREE_AFFINE_GENERATOR_WRONSKIAN.md), and [Wronskian-to-JC2 theorem](HC4_RANK_THREE_WRONSKIAN_TO_JC2.md) |
| Final regular `[4]` packet | [developable-image reduction](HC4_FINAL_RANK_THREE_DEVELOPABLE_IMAGE.md), [smooth-chart obstruction](HC4_FINAL_RANK_THREE_SMOOTH_CHART_OBSTRUCTION.md), [Krylov flag](HC4_FINAL_RANK_THREE_KRYLOV_FLAG.md), [Frobenius closure](HC4_FINAL_FROBENIUS_CLOSURE.md), [dimension-four trace/focal identity](HC4_TRACE_FOCAL_DIMENSION_FOUR_MIRACLE.md), and [affine-plane flatness closure](HC4_AFFINE_PLANE_SCHUBERT_BRIDGE.md) |

### Verification boundary

The registered checker replays the exact moving-frame and prolongation
calculations in the final regular `[4]` packet.  It is not an aggregate replay
of every implication in the proof map above.  The scalar, rank-at-most-two,
rank-three reduction, developability, and incidence steps remain written
proofs in their linked canonical notes.  `HC4MR1` currently has neither an
independent end-to-end replay nor external review; those assurance tasks are
separate from the theorem's stated scope.

Then `N` is polynomial, nilpotent, `S`-self-adjoint and Hessian-integrable:

\[
N^{\mathsf T}S=SN,\qquad SN=T=\operatorname{Hess}A.
\tag{0.3}
\]

---

## Master theorem

> **Theorem HC4-MR — relative-nilpotent master reduction.**
> Under (0.1), every generic Jordan stratum reduces globally to either an
> `HC2` packet or the exact cotangent lift of a plane Keller map.  In
> particular, the regular rank-three block
> \[
> \operatorname{rank}T=3,\qquad N\sim J_4(0),
> \tag{0.4}
> \]
> with a linearly-independent associated four-variable quasi-translation has
> no residual moving packet.  On this final stratum the complete Jordan flag
> \[
> \ker N\subset\ker N^2\subset\ker N^3
> \tag{0.5}
> \]
> is Frobenius-integrable, its middle distribution has affine two-plane
> leaves, and the complete flatness prolongation excludes rank-two kernel-line
> motion.  Rank zero is fixed; rank one either makes the middle plane constant
> or produces an affine-hyperplane pencil which again forces a constant affine
> invariant.  The latter alternatives lie in the already closed
> linearly-dependent packet.

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

## 6. Final regular `[4]` globalization

Frobenius alone still does *not* produce a constant affine flag: ordinary
Hessians are tied to the ambient flat affine coordinates, and the exact
first-order audit permits upper-triangular flag motion.  The additional flat
affine equations close that motion as follows.

In an `S`-adapted Jordan frame put

\[
a=\Gamma^3_{4,1}=\Gamma^2_{3,1},\quad
b=\Gamma^2_{4,1},\quad
p=\Gamma^4_{3,3},\quad
q=\Gamma^4_{4,2},\quad
r=\Gamma^3_{4,2}.
\]

The middle distribution `E2=ker N^2` is autoparallel.  Its two transverse
direction derivatives and the projective derivative of `E1=ker N` are

\[
A_3=\begin{pmatrix}0&-(a+q)/2\\0&0\end{pmatrix},\qquad
A_4=\begin{pmatrix}a&r\\0&q\end{pmatrix},
\tag{6.1}
\]

\[
D[e_1]=
\begin{pmatrix}
0&0&a&b\\
0&0&0&a\\
0&0&0&0
\end{pmatrix}.
\tag{6.2}
\]

The complete second-order flatness elimination gives

\[
a(p-q)=0,\qquad 4pa-3a^2-4aq+3q^2=0.
\tag{6.3}
\]

On rank-two motion, `HC4RSD72` makes `pq` a nonzero constant in the canonical
frame.  Prolonging by `d(pq)=0` adds

\[
p(2pa-aq+3q^2)=0.
\tag{6.4}
\]

Equations (6.3)--(6.4), saturated by `a`, generate the unit ideal.  Thus
rank-two motion is empty (`HC4RSD79`).

For lower motion, rank zero is `HC4RSD65`.  Rank one has

\[
a=q=0,\qquad b\ne0,qquad pr=0.
\tag{6.5}
\]

If `r=0`, (6.1) makes `E2` a constant two-plane.  If `p=0`, the affine second
fundamental form of `E3=ker N^3` vanishes, so its leaves are affine
hyperplanes.  Their projective closures form a line in the dual projective
space: the graph of the leaf-hyperplane map equals the incidence variety, and
the incidence degree is therefore one.  A constant pencil direction gives a
constant linear invariant.  A moving direction has conormal
`lambda0+t*lambda1`; differentiating `lambda(t)v(t)=0` shows that the moving
kernel line and its derivative lie in the fixed two-plane
`ker(lambda0) intersect ker(lambda1)`.  Hence `E2` is constant, contradicting
`r!=0`.  This proves `HC4RSD80`.

The exact calculations and the full incidence proof are in
`HC4_AFFINE_PLANE_SCHUBERT_BRIDGE.md`.  Reproduce the local certificates with

```bash
.venv/bin/python scripts/verify_hc4_affine_plane_bridge.py
.venv/bin/python scripts/verify_hc4_affine_plane_prolongation.py
```

Thus the former global affine-or-Keller bridge has no surviving regular `[4]`
packet.

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
branch: every stratum runs *back* to `HC2` or this exact plane-cotangent
geometry.  The complete relative-nilpotent branch is therefore equivalent, in
obstruction content, to `JC2`.

This does not prove unrestricted `HC4` or `JC2`: it closes the
relative-nilpotent Hessian-pencil branch and identifies its only possible
endpoint with the still-open plane Keller problem.

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
&\to\text{Frobenius flag}\to\text{affine }E_2
\to\text{flatness/hyperplane pencil}\to HC2/JC2.
\end{cases}}
\]

The complete reduction is degree-free.

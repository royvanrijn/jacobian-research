# Rank-two apex reduction for the moving `HC4` pencil

## Status and scope

**Current status.**  This is an intermediate proof-map note.  Its
nonhomogeneous `[2,2]` remainder is closed by `HC4RSD58--60`, the `[3,1]`
remainder is closed by `HC4RSD61--63`, and both sit inside the all-degree
relative-nilpotent master reduction `HC4MR1`.  Section 5 records the
historical handoff at `HC4RSD56--57`, not a live frontier.

This note continues `HC4RSD17`, `HC4RSD42`, `HC4RSD54`, and the small-rank
Hessian results of de Bondt.  It separates the two generic rank-two Jordan
types

\[
[2,2]\qquad\text{and}\qquad[3,1].
\]

The square-zero type admits a strong global reduction.

> **Theorem HC4RSD56 — rank-two square-zero apex reduction.**
> Let
>
> \[
> S=\operatorname{Hess}\psi,
> \qquad T=\operatorname{Hess}A,
> \qquad \det S\in K^*,
> \qquad N=S^{-1}T,
> \]
>
> over a characteristic-zero field.  Assume
>
> \[
> \operatorname{rank}T=2,
> \qquad N^2=0.
> \tag{0.1}
> \]
>
> Then the image of `\nabla A` has a nonzero constant projective apex `p`.
> For every such apex,
>
> \[
> k:=S^{-1}p
> \]
>
> satisfies
>
> \[
> Tk=0,
> \qquad p^{\mathsf T}k=0,
> \qquad (S+sT)k=p\quad\text{for every }s.
> \tag{0.2}
> \]
>
> After a constant coordinate change with `p=e_4`,
>
> \[
> (S^{-1})_{44}=0,
> \qquad
> \det\operatorname{Hess}_{x_1,x_2,x_3}\psi=0.
> \tag{0.3}
> \]
>
> Thus the moving `[2,2]` flag is controlled by a ternary singular-Hessian
> family plus the polynomial common kernel section `k`.

> **Corollary HC4RSD57 — homogeneous moving `[2,2]` closes in all degrees.**
> Under HC4RSD56, if `A` is homogeneous, then the pencil is exactly the `JC2`
> cotangent endpoint after a constant linear change:
>
> \[
> A=A(x,y),
> \qquad
> \psi=zP(x,y)+wQ(x,y)+R(x,y),
> \tag{0.4}
> \]
>
> with
>
> \[
> \det J(P,Q)\in K^*.
> \tag{0.5}
> \]
>
> Hence no homogeneous rank-two square-zero moving flag is `HC4`-specific,
> regardless of degree.

## 1. Constant apex from Hessian rank two

The polynomial map `\nabla A` has Jacobian `T` and Jacobian rank two.
De Bondt's small-rank Hessian theorem states that every polynomial Hessian of
rank at most two has a projective image apex.  Thus there is a constant
nonzero vector `p` such that the Zariski closure of `im(\nabla A)` is invariant
under translation by `p`.

At a generic smooth point of this image surface, its tangent plane is
`im T`.  Therefore

\[
p\in\operatorname{im}T
\tag{1.1}
\]

over the fraction field.  Write `p=Tq=SNq`.

Put

\[
k=S^{-1}p=Nq.
\]

Since `N^2=0`,

\[
Tk=SN^2q=0.
\tag{1.2}
\]

Self-adjointness `N^{\mathsf T}S=SN` gives

\[
p^{\mathsf T}k
=(SNq)^{\mathsf T}Nq
=q^{\mathsf T}N^{\mathsf T}SNq
=q^{\mathsf T}SN^2q
=0.
\tag{1.3}
\]

Finally

\[
(S+sT)k=Sk+sTk=p,
\tag{1.4}
\]

so `k` is the same inverse-Jacobian column for **every** member of the pencil.
This proves (0.2).

## 2. Ternary singular-Hessian reduction

Make a constant coordinate change so that the transformed apex is `e_4`.
The scalar identity in (1.3) becomes

\[
(S^{-1})_{44}=0.
\]

Since `det S=\delta\in K^*`, Cramer's rule gives

\[
(S^{-1})_{44}
=\delta^{-1}
 \det\operatorname{Hess}_{x_1,x_2,x_3}\psi.
\]

Hence

\[
\det\operatorname{Hess}_{x_1,x_2,x_3}\psi=0.
\]

Moreover the fourth component of `k=S^{-1}e_4` is zero and the first three
components form a polynomial kernel vector of this ternary Hessian.  Equation
(1.2) says that the same vector also annihilates `Hess A`.

This is substantially narrower than an arbitrary moving null-plane: the
remaining movement lives inside the classical three-variable singular-Hessian
classification.

## 3. Homogeneous `A` supplies two constant apices

For homogeneous `A`, de Bondt's homogeneous small-rank theorem is stronger:
a homogeneous Hessian of rank two has at least two independent projective
image apices.

Consequently `im T` contains a constant two-plane.  Since `rank T=2`,

\[
\operatorname{im}T=L
\]

is itself a constant plane over the fraction field.  Symmetry of `T` gives

\[
\ker T=L^\perp,
\]

also a constant two-plane.

After a constant linear change and removal of affine terms,

\[
A=A(x,y),
\qquad
T=
\begin{pmatrix}
H&0\\0&0
\end{pmatrix},
\tag{3.1}
\]

where `H=Hess A` is generically invertible.

## 4. Square-zero relative nilpotence forces the cotangent block

Write the inverse metric in the same `2+2` splitting as

\[
S^{-1}=
\begin{pmatrix}
K&B\\B^{\mathsf T}&D
\end{pmatrix}.
\tag{4.1}
\]

Since

\[
SN^2=T S^{-1}T,
\]

the condition `N^2=0` gives

\[
0=T S^{-1}T
 =
\begin{pmatrix}
HKH&0\\0&0
\end{pmatrix}.
\]

Generic invertibility of `H` forces

\[
K=0.
\tag{4.2}
\]

Thus

\[
S^{-1}=
\begin{pmatrix}
0&B\\B^{\mathsf T}&D
\end{pmatrix}.
\]

Invertibility of `S^{-1}` forces `B` invertible.  Write

\[
S=
\begin{pmatrix}
P&Q\\Q^{\mathsf T}&R
\end{pmatrix}.
\]

Multiplying `SS^{-1}=I`, the lower-left block is

\[
RB^{\mathsf T}=0,
\]

hence

\[
R=0.
\tag{4.3}
\]

Therefore all second derivatives of `\psi` in the passive variables `z,w`
vanish.  Hessian integrability yields

\[
\psi=zP_0(x,y)+wQ_0(x,y)+R_0(x,y)
\]

up to affine terms.

The determinant of this Hessian is

\[
\det S=\det J(P_0,Q_0)^2
\]

(up to the harmless sign convention from the `2+2` ordering).  Since
`det S` is a nonzero constant,

\[
\det J(P_0,Q_0)\in K^*.
\]

This is exactly the plane-Keller cotangent packet and proves HC4RSD57.

## 5. Historical rank-two handoff (now closed)

The homogeneous square-zero `[2,2]` branch is now closed.  The genuinely new
rank-two possibilities are:

1. nonhomogeneous `[2,2]`, where de Bondt guarantees one constant apex but a
   second apex need not exist;
2. `[3,1]`, where `N^2` has rank one and the square-zero block argument is
   replaced by a rank-one Schur defect.

The first case now carries the canonical common polynomial vector field

\[
k=S^{-1}p,
\qquad
(S+sT)k=p,
\qquad
Tk=0,
\]

and a ternary singular-Hessian metric.  This should be attacked using the
complete three-variable singular-Hessian normal forms rather than another
four-variable coefficient census.

## 6. External input

The small-rank apex statements are from Michiel de Bondt,
*Polynomial Hessians with small rank*, arXiv:1609.03904v2 (2022).  In the
notation of that paper:

- Hessian rank `1<=r<=2` gives an apex at infinity for the image of the
  gradient;
- a homogeneous Hessian of rank `2<=r<=4` gives two independent apices at
  infinity.

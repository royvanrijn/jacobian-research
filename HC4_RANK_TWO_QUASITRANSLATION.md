# Quasi-translation structure in the moving rank-two `[2,2]` HC4 branch

## Status and scope

This note continues `HC4RSD56` after the constant projective apex has been
normalized to

\[
p=e_w.
\]

Write

\[
S=\operatorname{Hess}\psi,
\qquad T=\operatorname{Hess}A,
\qquad N=S^{-1}T,
\qquad \operatorname{rank}T=2,
\qquad N^2=0,
\]

and let

\[
k=S^{-1}p.
\]

`HC4RSD56` gives

\[
Tk=0,
\qquad p^{\mathsf T}k=0,
\qquad (S+sT)k=p.
\tag{0.1}
\]

> **Theorem HC4RSD58 — canonical quasi-translation.**
> In the above branch, write
>
> \[
> k=(\bar k,0),
> \qquad
> D=\bar k\cdot\nabla_{x,y,z}.
> \]
>
> On the generic rank-two locus of the ternary Hessian
>
> \[
> H=\operatorname{Hess}_{x,y,z}\psi,
> \]
>
> one has
>
> \[
> D\bar k=0.
> \tag{0.2}
> \]
>
> Consequently
>
> \[
> D^2x=D^2y=D^2z=0,
> \qquad Dw=0,
> \tag{0.3}
> \]
>
> and, more strongly,
>
> \[
> D(\nabla\psi)=e_w,
> \qquad
> D(\nabla A)=0.
> \tag{0.4}
> \]
>
> In particular
>
> \[
> D(\psi_w)=1,
> \tag{0.5}
> \]
>
> so `D` has the polynomial slice `psi_w`.  The moving `[2,2]` branch is
> therefore a quasi-translation / `G_a`-quotient problem constrained by
> Hessian integrability, rather than an arbitrary moving nilpotent-matrix
> problem.

No claim is made here that an arbitrary locally nilpotent derivation with a
slice in four variables has polynomial invariant ring.  That is too strong in
this generality; the point is that our derivation has the additional Hessian
identities (0.4).

## 1. Block equations

Since `p=e_w` and `p^T k=0`, write

\[
k=(\bar k,0).
\]

Decompose the Hessian metric as

\[
S=
\begin{pmatrix}
H&m\\
m^{\mathsf T}&r
\end{pmatrix},
\tag{1.1}
\]

where

\[
H=\operatorname{Hess}_{x,y,z}\psi,
\qquad
m=\nabla_{x,y,z}\psi_w.
\]

The identity `Sk=e_w` becomes

\[
H\bar k=0,
\qquad
m^{\mathsf T}\bar k=1.
\tag{1.2}
\]

`HC4RSD56` already gives `det H=0`.  Since `S` is invertible and the branch
has generic rank two in the ternary kernel reduction, the generic kernel of
`H` is exactly the line spanned by `\bar k`.

## 2. Differentiate the Hessian-kernel equation along `D`

Apply

\[
D=\sum_{a=1}^3\bar k_a\partial_a
\]

to

\[
H\bar k=0.
\]

This gives

\[
(DH)\bar k+H(D\bar k)=0.
\tag{2.1}
\]

For each component `i`, symmetry of the third derivatives of `psi` gives

\[
\begin{aligned}
((DH)\bar k)_i
&=\sum_{a,b}\bar k_a\,\psi_{iab}\,\bar k_b\\
&=\partial_i\!\left(\frac12\bar k^{\mathsf T}H\bar k\right)
   -\sum_{a,b}(\partial_i\bar k_a)H_{ab}\bar k_b.
\end{aligned}
\]

The second term vanishes because `H\bar k=0`, and the first vanishes because
`\bar k^T H\bar k=0`.  Therefore

\[
(DH)\bar k=0.
\tag{2.2}
\]

Equation (2.1) reduces to

\[
H(D\bar k)=0.
\]

Since the generic kernel of `H` is one-dimensional,

\[
D\bar k=\alpha\bar k
\tag{2.3}
\]

for some rational function `alpha`.

## 3. The normalization `m^T k=1` kills the scaling

Apply `D` to the second equation of (1.2):

\[
D(m^{\mathsf T}\bar k)=0.
\]

Hence

\[
(Dm)^{\mathsf T}\bar k+m^{\mathsf T}D\bar k=0.
\tag{3.1}
\]

Again by symmetry of third derivatives,

\[
(Dm)^{\mathsf T}\bar k
=\bar k^{\mathsf T}(\partial_wH)\bar k.
\tag{3.2}
\]

Differentiate the identity `\bar k^T H\bar k=0` with respect to `w`.  The
terms involving `\partial_w\bar k` vanish because `H\bar k=0`, leaving

\[
\bar k^{\mathsf T}(\partial_wH)\bar k=0.
\tag{3.3}
\]

Thus (3.1) and (2.3) give

\[
0=m^{\mathsf T}(\alpha\bar k)=\alpha.
\]

Therefore

\[
\boxed{D\bar k=0.}
\]

This proves (0.2).

## 4. Gradient identities

Since `k` is the coefficient vector of `D` and `Sk=e_w`, Hessian symmetry
gives

\[
D(\nabla\psi)=Sk=e_w.
\tag{4.1}
\]

Likewise `Tk=0` gives

\[
D(\nabla A)=Tk=0.
\tag{4.2}
\]

In components,

\[
D(\psi_x)=D(\psi_y)=D(\psi_z)=0,
\qquad
D(\psi_w)=1,
\tag{4.3}
\]

and every component of `\nabla A` is a `D`-invariant.

Finally (0.2) says that each coefficient of `D` is itself invariant.  Hence
on the coordinate generators

\[
D^2x_i=D(\bar k_i)=0,
\]

so `D` is a quasi-translation derivation.  In characteristic zero it is
locally nilpotent, with explicit exponential action

\[
(x,y,z,w)\longmapsto(x,y,z,w)+t(k_1,k_2,k_3,0).
\tag{4.4}
\]

Because the coefficients `k_i` are invariant, (4.4) is exact despite the fact
that the direction varies from orbit to orbit.

## 5. What this buys us

The moving direction is now constrained in three simultaneous ways:

1. `k` spans `ker Hess_{xyz} psi`;
2. `k` spans a subdirection of `ker Hess A`;
3. `k` is invariant along its own flow and has the slice `psi_w`.

The remaining nonhomogeneous `[2,2]` problem should therefore be attacked by
combining the ternary singular-Hessian classification with (4.1)--(4.3).
The generic locally-nilpotent-derivation slice problem alone is not enough;
Hessian integrability is the extra rigidity available here.

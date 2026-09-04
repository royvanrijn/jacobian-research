# Affine moving-kernel Hessian pencils for `HC4`

## Status

This note continues
[`HC4RSD1`](HC4_REVERSE_SCHUR_DESCENT.md) beyond constant kernel directions.
It closes the first genuinely moving scalar stratum: the kernel vector may
depend affinely on the four reduced variables, but its normalized line is
fixed as the pencil parameter varies.

> **Theorem `HC4RSD2` (fixed affine moving-kernel obstruction).** Let \(K\)
> have characteristic zero and put
>
> \[
> \Phi(t,x)=\frac{\lambda}{2}t^2+tA(x)+B(x),
> \qquad x\in K^4.
> \]
>
> Assume
>
> \[
> \det\operatorname{Hess}\Phi=c\in K^\times,
> \qquad
> \det\operatorname{Hess}_x(B+sA)=0.                 \tag{0.1}
> \]
>
> Suppose the generic kernel line admits a unimodular generator affine in
> \(x\), and that this normalized line is independent of \(s\). Then either
> the kernel direction is constant, hence is handled by `HC4RSD1`, or an
> affine change of variables puts the family in the form
>
> \[
> \begin{aligned}
> A&=\alpha y+p(z)w+q(z),\\
> B&=yC(z)+(x-yz)C'(z)
>   +\frac{\gamma}{2}w^2+h(z)w+k(z),\\
> C(z)&=\frac{\eta}{2}z^2+\xi z+\zeta,
> \end{aligned}                                             \tag{0.2}
> \]
>
> where \(\alpha\eta\gamma\ne0\). Every scalar Schur descendant
>
> \[
> \psi_{\kappa,\mu}=B+\frac{\kappa}{2}A^2+\mu A,
> \qquad \kappa\ne0,                                      \tag{0.3}
> \]
>
> has determinant
>
> \[
> \det\operatorname{Hess}\psi_{\kappa,\mu}
> =-\kappa\alpha^2\eta^2\gamma,                           \tag{0.4}
> \]
>
> and its gradient is a polynomial automorphism. Thus this stratum has no
> collision and meets no affine-degree-two or affine-degree-three `HC4`
> packet.

The theorem itself does not classify affine \(x\)-kernel vectors whose
projective line varies with \(s\). That branch is subsequently excluded in
full by [`HC4RSD3`](HC4_PARAMETER_MOVING_AFFINE_KERNEL_PENCILS.md).

The verifier is
[`scripts/verify_hc4_affine_moving_kernel_pencils.py`](scripts/verify_hc4_affine_moving_kernel_pencils.py),
and its generated ledger is
[`artifacts/generated-results/hc4_affine_moving_kernel_pencils.json`](artifacts/generated-results/hc4_affine_moving_kernel_pencils.json).

## 1. The unimodular kernel and Piola equation

Work over \(F=K(s)\), put \(R=F[x_1,x_2,x_3,x_4]\), and write

\[
M=\operatorname{Hess}_x(B+sA),\qquad g=\nabla A.
\]

The bordered unit forces \(\operatorname{rank}M=3\). If
\(N=\operatorname{adj}(M)\), then

\[
N=N^{\mathsf T},\qquad \operatorname{rank}N=1,
\qquad g^{\mathsf T}Ng=-c\in F^\times.                    \tag{1.1}
\]

Put \(h=Ng\). Since \(g^{\mathsf T}h=-c\) is a unit, \(h\) is a
unimodular vector. Over the fraction field, every column of \(N\) lies on
the line spanned by \(h\). A row \(\ell\) with \(\ell h=1\) shows that each
column's multiplier is already in \(R\), so \(\operatorname{im}N=Rh\).
Writing \(N=hr^{\mathsf T}\), symmetry and the same splitting give
\(r=\epsilon h\); (1.1) makes \(\epsilon\) a unit. Thus, without any
freeness theorem,

\[
\boxed{N=\epsilon vv^{\mathsf T}},
\qquad \epsilon\in F^\times,                              \tag{1.2}
\]

with \(v\in R^4\) unimodular. The Hessian Piola identity says that every
row of \(N\) is divergence-free. Substitution into (1.2) gives

\[
\boxed{D_vv+(\operatorname{div}v)v=0.}                    \tag{1.3}
\]

Suppose \(v=a+Lx\) and put \(\tau=\operatorname{tr}L\). Constant and linear
coefficients of (1.3) give

\[
(L+\tau I)a=0,
\qquad (L+\tau I)L=0.                                    \tag{1.4}
\]

If \(\tau\ne0\), the second equation makes \(L\) diagonalizable with
eigenvalues \(0,-\tau\). If \(-\tau\) has multiplicity \(r\), then
\(\tau=\operatorname{tr}L=-r\tau\), impossible in characteristic zero.
Thus

\[
\boxed{\operatorname{tr}L=0,\qquad L^2=0,\qquad La=0.}    \tag{1.5}
\]

If \(\operatorname{rank}L=2\), then
\(\operatorname{im}L=\ker L\), so \(a\in\operatorname{im}L\) and
\(a+Lx\) has a zero. This contradicts unimodularity. Rank one similarly
requires \(a\notin\operatorname{im}L\), and affine linear normalization
gives the unique nonconstant orbit

\[
\boxed{v=(z,1,0,0)^{\mathsf T}.}                           \tag{1.6}
\]

## 2. Integrating the moving kernel

Use coordinates \((x,y,z,w)\) and put

\[
\delta=z\partial_x+\partial_y,
\qquad r=x-yz,
\qquad \ker\delta=F[z,w,r].                               \tag{2.1}
\]

The equation \(\operatorname{Hess}(f)v=0\) says that every component of the
gradient of \(f\) lies in \(\ker\delta\). Write those components as
\(P,Q,R,T\in F[z,w,r]\). Mixed-partial equality gives

\[
Q_r=-zP_r,\qquad R_r=P_z-yP_r,\qquad T_r=P_w.              \tag{2.2}
\]

The middle left side is invariant and has no \(y\)-term, so \(P_r=0\).
The remaining \((yz,yw,zw)\) equations then force

\[
P=C'(z),\qquad Q=C(z)-zC'(z),
\]

and integrate the other coordinates through one arbitrary \(G(z,w)\).
Hence the complete potential is

\[
\boxed{f=yC(z)+(x-yz)C'(z)+G(z,w).}                       \tag{2.3}
\]

Direct calculation gives

\[
\operatorname{adj}(\operatorname{Hess}f)
=-C''(z)^2G_{ww}(z,w)\,vv^{\mathsf T}.                    \tag{2.4}
\]

Thus the generic Hessian rank is three precisely when
\(C''G_{ww}\ne0\).

## 3. The pencil and bordered unit

Assume the same normalized kernel (1.6) works for \(f_s=B+sA\). Apply
(2.3) separately to \(A\) and \(B\):

\[
\begin{aligned}
A&=yC_A(z)+(x-yz)C_A'(z)+G_A(z,w),\\
B&=yC_B(z)+(x-yz)C_B'(z)+G_B(z,w).
\end{aligned}                                               \tag{3.1}
\]

Since \(v^{\mathsf T}\nabla A=C_A(z)\), the bordered determinant is

\[
C_A(z)^2
\bigl(C_B''+sC_A''\bigr)^2
\bigl(G_{B,ww}+sG_{A,ww}\bigr)=c.                        \tag{3.2}
\]

This is an identity in the UFD \(K[s,z,w]\), with a nonzero constant on
the right. Each displayed factor is a unit, so

\[
C_A=\alpha\in K^\times,\quad
C_B''=\eta\in K^\times,\quad
G_{A,ww}=0,\quad G_{B,ww}=\gamma\in K^\times.             \tag{3.3}
\]

Integrating (3.3) gives exactly (0.2), with

\[
c=\alpha^2\eta^2\gamma.                                  \tag{3.4}
\]

No degree bound on \(p,q,h,k\) was used.

## 4. Explicit inverse of every descendant

Put \(S=\kappa A+\mu\), and order the gradient coordinates as
\((F_x,F_y,F_z,F_w)\). Equations (0.2)--(0.3) give

\[
\begin{aligned}
F_x&=\eta z+\xi,\\
F_y&=C-zC'+\alpha S,\\
F_w&=\gamma w+h+Sp,\\
F_z&=\eta(x-yz)+h'w+k'+S(p'w+q').
\end{aligned}                                               \tag{4.1}
\]

These coordinates recover, in order,

\[
z,\qquad S,\qquad w,\qquad
A=\frac{S-\mu}{\kappa},\qquad
y=\frac{A-pw-q}{\alpha},\qquad x.                        \tag{4.2}
\]

Every division is by one of the constant units
\(\alpha,\eta,\gamma,\kappa\). Thus (4.2) is a polynomial inverse. It
follows directly that there is no collision and hence no surviving
collision component to reconstruct over \(\mathbb Q\).

## 5. Parameter motion, subsequently closed

The next branch cannot be absorbed into the preceding coordinates. Over
\(K(x)\), the first nonconstant kernel of a symmetric \(4\)-by-\(4\) pencil
of normal rank three has minimal index one:

\[
v(s,x)=v_0(x)+sv_1(x).                                   \tag{5.1}
\]

The elementary calibration

\[
B=xz+\frac12w^2,\qquad A=yz,\qquad v=(-s,1,0,0)           \tag{5.2}
\]

has such parameter motion but bordered determinant \(z^2\), so it fails the
unit gate. The next coefficient system is

\[
\begin{aligned}
\operatorname{Hess}(B)v_0&=0,\\
\operatorname{Hess}(B)v_1+\operatorname{Hess}(A)v_0&=0,\\
\operatorname{Hess}(A)v_1&=0,\\
v_1^{\mathsf T}\nabla A&=0,\\
v_0^{\mathsf T}\nabla A&\in K^\times,
\end{aligned}                                               \tag{5.3}
\]

together with coefficientwise Piola equations, unimodularity, and generic
corank one. The classification of this system and its empty terminal chart
are proved in
[`HC4_PARAMETER_MOVING_AFFINE_KERNEL_PENCILS.md`](HC4_PARAMETER_MOVING_AFFINE_KERNEL_PENCILS.md).

## 6. Reproduction

Run:

```bash
.venv/bin/python scripts/verify_hc4_affine_moving_kernel_pencils.py
# cleanup only: verify committed inputs and later handoffs
.venv/bin/python scripts/verify_hc4_affine_moving_kernel_pencils.py --audit-existing-only
```

The command verifies the affine Piola normal-form calibrations, the integrated
Hessian and adjugate identities, the UFD-normalized pencil, its constant
determinants, and the complete triangular inverse. It writes
`artifacts/generated-results/hc4_affine_moving_kernel_pencils.json`.
The cleanup-only mode instead hash-checks that ledger, the projective-polar
atlas it consumes, and the imported equation helper. It neither imports
SymPy nor recomputes or rewrites the ledger, and it reports that `HC4RSD3--4`
closed the artifact's recorded affine handoffs.

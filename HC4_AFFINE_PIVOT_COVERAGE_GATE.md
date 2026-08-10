# Affine-pivot coverage gate for reverse Schur descent

## Status

The singular-pencil results HC4RSD1--HC4RSD5 classify several large kernel
strata, but they do not show which direct HC4 candidates admit such a
Schur presentation. This note gives an exact coverage criterion for an
affine scalar pivot and intersects it with the first open direct quintic
packet.

> **Theorem HC4RSD6 (affine-pivot coverage gate).** Let \(K\) have
> characteristic zero, let
>
> \[
> \psi\in K[x_1,x_2,x_3,x_4],
> \qquad
> H=\operatorname{Hess}\psi,
> \qquad
> \det H=\delta\in K^\times,
> \]
>
> and fix a nonzero constant covector \(\ell\in K^4\). Put
>
> \[
> N_\ell=\ell^{\mathsf T}\operatorname{adj}(H)\ell.
> \tag{0.1}
> \]
>
> There is an affine-pivot singular Schur presentation
>
> \[
> A=\ell\mathbin{\cdot}x+a_0,\qquad
> B=\psi-\frac{\kappa}{2}A^2-\mu A,\qquad
> \Phi=tA+B,
> \tag{0.2}
> \]
>
> with \(\kappa\in K^\times\), constant nonzero
> \(\det\operatorname{Hess}\Phi\), and identically singular reduced pencil
> \(\operatorname{Hess}(B+sA)\), if and only if
>
> \[
> N_\ell=\nu\in K^\times.
> \tag{0.3}
> \]
>
> In that case \(\kappa=\delta/\nu\),
> \(\det\operatorname{Hess}\Phi=-\nu\), and the Schur descendant is exactly
> \(\psi\).
>
> Suppose further that
> \(\psi=q_2+h_3+h_4+h_5\) is a collision-normalized quintic candidate and
> \(\operatorname{Hess}(h_5)\) has essential rank three. In constant-kernel
> coordinates \((u,t)\), write
>
> \[
> C=\operatorname{Hess}_u(h_5),\qquad
> h_4=t\,s_3(u)+r_4(u),\qquad
> d=\nabla s_3,
> \]
>
> and define
>
> \[
> \Delta=\det C,\qquad
> w=\operatorname{adj}(C)d.
> \tag{0.4}
> \]
>
> If \(\psi\) admits the affine-pivot presentation above, then
> \(\ell=(a,0)\) for a nonzero constant \(a\in K^3\), and
>
> \[
> a^{\mathsf T}w=0.
> \tag{0.5}
> \]
>
> Thus affine-pivot coverage of the open nonsquarefree rank-three packet is
> contained in the closed locus where the three components of \(w\) have a
> nontrivial constant linear relation.

The exact verifier is
[scripts/verify_hc4_affine_pivot_coverage_gate.py](scripts/verify_hc4_affine_pivot_coverage_gate.py),
and its generated ledger is
[artifacts/generated-results/hc4_affine_pivot_coverage_gate.json](artifacts/generated-results/hc4_affine_pivot_coverage_gate.json).

## 1. Exact affine-pivot criterion

Let \(A=\ell\mathbin{\cdot}x+a_0\). Since
\(\operatorname{Hess}A=0\), the reduced pencil in (0.2) is independent of
its parameter:

\[
M=\operatorname{Hess}(B+sA)=H-\kappa\ell\ell^{\mathsf T}.
\tag{1.1}
\]

The rank-one determinant identity gives

\[
\det M=\delta-\kappa N_\ell.
\tag{1.2}
\]

Therefore \(M\) is singular precisely when
\(N_\ell=\delta/\kappa\), which is a nonzero constant. Conversely, if
(0.3) holds, taking \(\kappa=\delta/\nu\) makes \(M\) singular.

The bordered determinant is insensitive to the rank-one correction:

\[
\det
\begin{pmatrix}
0&\ell^{\mathsf T}\\
\ell&H-\kappa\ell\ell^{\mathsf T}
\end{pmatrix}
=-N_\ell.
\tag{1.3}
\]

Hence the parent has constant nonzero Hessian determinant \(-\nu\).
Finally,

\[
B+\frac{\kappa}{2}A^2+\mu A=\psi,
\tag{1.4}
\]

so no extra determinant or reconstruction condition is hidden in the
criterion.

Equivalently, (0.3) says that \(\ell\) has nonzero constant squared norm in
the inverse Hessian metric:

\[
\ell^{\mathsf T}H^{-1}\ell=\frac{\nu}{\delta}.
\tag{1.5}
\]

This is the exact affine-pivot coverage test for a given direct candidate.

The zero-norm branch has a different meaning and is not rejected by this
theorem. If `N_ell=0`, then
`A=(ell.x)^2/2` is itself a rank-one square-zero constant-Hessian pencil
direction. It is the recognition frontend `HC4MR3`, recorded in the
[relative-pencil master note](HC4_RELATIVE_NILPOTENT_MASTER_REDUCTION.md#first-exact-recognition-frontend).
Thus the metric numerator has a useful trichotomy: nonzero constant gives the
affine singular-Schur presentation above, zero gives direct admission to
`PHC4`, and a nonconstant value gives neither conclusion.

## 2. Marked collision transfer

Suppose \(p_+\ne p_-\) and
\(\nabla\psi(p_+)=\nabla\psi(p_-)\). The same pair lifts to equal parent
gradients at a common pivot value exactly when

\[
A(p_+)=A(p_-),
\qquad\text{equivalently}\qquad
\ell\mathbin{\cdot}(p_+-p_-)=0.
\tag{2.1}
\]

Indeed,

\[
\nabla B
=\nabla\psi-(\kappa A+\mu)\ell,
\]

so (2.1) makes the two \(B\)-gradients equal, while
\(\partial_t\Phi=A\). For an antipodal collision at \(\pm p\), the added
condition is simply \(\ell\mathbin{\cdot}p=0\).

The determinant coverage test and the marked-collision hyperplane should
be imposed separately: the first decides whether a singular affine lift
exists, and the second decides whether that lift inherits the chosen
collision.

## 3. Rank-three quintic intersection

Scale the spatial variables by a parameter \(\rho\). In the
constant-kernel coordinates for \(h_5\), the leading Hessian has the block
form

\[
H(\rho)=
\begin{pmatrix}
\rho^3 C+\rho^2E+\cdots&
\rho^2d+\cdots\\
\rho^2d^{\mathsf T}+\cdots&
\rho^2e+\rho f+\cdots
\end{pmatrix}.
\tag{3.1}
\]

Write \(\ell=(a,\tau)\). The degree-nine part of \(N_\ell\) is

\[
[N_\ell]_{\rho^9}=\tau^2\Delta.
\tag{3.2}
\]

Since \(\Delta\ne0\) and \(N_\ell\) is constant, (3.2) forces
\(\tau=0\). Thus any affine pivot is active: it annihilates the constant
kernel of the top quintic.

The degree-eleven face of \(\det H(\rho)=\delta\) is

\[
\Delta e=0.
\tag{3.3}
\]

Hence \(e=D_t^2h_4=0\) and
\(h_4=t\,s_3(u)+r_4(u)\). The next determinant face is the existing
rank-three Schur equation

\[
R=\Delta f-d^{\mathsf T}\operatorname{adj}(C)d=0,
\qquad
f=D_t^2h_3.
\tag{3.4}
\]

Let \(N_7=[N_\ell]_{\rho^7}\) after \(\tau=e=0\). Direct block-adjugate
calculation gives the polynomial identity

\[
\boxed{
\Delta N_7
=R\bigl(a^{\mathsf T}\operatorname{adj}(C)a\bigr)
+\bigl(a^{\mathsf T}\operatorname{adj}(C)d\bigr)^2.
}
\tag{3.5}
\]

Both \(R\) and \(N_7\) vanish for a constant-Hessian candidate admitting
the affine lift. Therefore (3.5), characteristic zero, and the domain
property give

\[
a^{\mathsf T}\operatorname{adj}(C)d=a^{\mathsf T}w=0,
\]

which proves (0.5).

## 4. A finite closed coverage locus

The entries of \(C\) have degree three and those of \(d\) have degree two,
so every component of \(w\) is a ternary form of degree eight. Write their
coefficients in the 45 degree-eight monomials as a \(3\)-by-\(45\) matrix
\(W\). Then

\[
\exists\,a\in K^3\setminus\{0\}:a^{\mathsf T}w=0
\quad\Longleftrightarrow\quad
\operatorname{rank}W\le2.
\tag{4.1}
\]

Consequently affine-pivot coverage is cut out by the \(3\)-by-\(3\) minors
of \(W\), together with the Hessian-integrability, nonsquarefree
discriminant, and Schur equations already defining the direct packet.
This turns the coverage question into a finite elimination problem rather
than a classification of arbitrary nonlinear kernel fields.

The diagonal nonsquarefree top admits a complete leading classification.
Take

\[
C=\operatorname{diag}(u_1^3,u_2^3,u_3^3).
\tag{4.2}
\]

If \(d=\nabla s_3\), Schur divisibility says

\[
u_1^3u_2^3u_3^3\mid
d_1^2u_2^3u_3^3+d_2^2u_1^3u_3^3+d_3^2u_1^3u_2^3.
\tag{4.3}
\]

Reducing (4.3) modulo each coordinate cube gives
\(u_i^3\mid d_i^2\), hence \(u_i^2\mid d_i\). Since each \(d_i\) is
quadratic and the three components are a gradient, this forces

\[
s_3=\frac{\alpha u_1^3+\beta u_2^3+\gamma u_3^3}{3},
\qquad
f=\alpha^2u_1+\beta^2u_2+\gamma^2u_3.
\tag{4.4}
\]

The corresponding Schur vector is

\[
w=(\alpha u_1^2u_2^3u_3^3,\,
   \beta u_1^3u_2^2u_3^3,\,
   \gamma u_1^3u_2^3u_3^2)^{\mathsf T}.
\tag{4.5}
\]

Its three distinguished coefficient channels form
\(\operatorname{diag}(\alpha,\beta,\gamma)\). Therefore the complete
affine-coverage condition on this leading Schur family is

\[
\boxed{\alpha\beta\gamma=0.}
\tag{4.6}
\]

At \(\alpha=\beta=\gamma=1\), equation (3.4) holds but the coefficient
rank is three, so Schur divisibility alone does not imply affine-pivot
coverage. These are leading-face statements, not assertions that any of
the displayed pairs prolongs to a full HC4 candidate.

The zero-norm recognition branch is strictly smaller than (4.6).  The
full lower-face calculation `HC4MR4` proves that every nonaligned member of
the diagonal Schur packet has empty constant-null-covector scheme, even when
one or two of `alpha,beta,gamma` vanish.  Thus (4.6) remains the exact leading
gate for a **nonzero constant** affine metric numerator, but it does not
produce a rank-one zero-norm pencil.  See the
[direct application](HC4_RELATIVE_NILPOTENT_MASTER_REDUCTION.md#application-to-the-direct-diagonal-quintic-packet).

## 5. Reproduction and next target

Run:

~~~bash
.venv/bin/python scripts/verify_hc4_affine_pivot_coverage_gate.py
~~~

The command verifies the universal rank-one and bordered determinant
identities, exact collision transfer, all rank-three homogeneous faces,
identity (3.5), and the full-span diagonal calibration.

The next bounded calculation is now precise: impose

\[
R=0,\qquad
\gcd(\Delta,\partial\Delta)\ne1,\qquad
\operatorname{rank}W\le2
\tag{5.1}
\]

with \(C=\operatorname{Hess}(h_5)\) and \(d=\nabla s_3\), and classify the
resulting constant-span-deficient components before adding lower faces and
the marked collision. Essential ranks one and two require analogous
metric-face calculations. General nonlinear pivots and nonsingular
specialized-determinant cancellation remain outside this theorem.

For inherited collision transfer, this component calculation is
subsequently bypassed by `HC4RSD7` in
[`HC4_AFFINE_PIVOT_COLLISION_FIBERS.md`](HC4_AFFINE_PIVOT_COLLISION_FIBERS.md):
`HC3` makes every affine metric fiber collision-free. The locus (5.1)
remains relevant only for classifying affine Schur representations whose
collisions, if any, lie on different pivot fibers.

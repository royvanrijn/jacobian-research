# HC4 motion-frame correction and the surviving negative sign

## 1. Result and status effect

The final regular Jordan block in
[the HC4 master reduction](HC4_RELATIVE_NILPOTENT_MASTER_REDUCTION.md)
was closed using an unjustified differential equation. The comparison of
the frozen normal frame with the moving Jordan frame introduces the factor
$a^2$. The normalized motion determinant is $pq/a^2$, not $pq$.
Consequently the claimed deduction $d(pq)=0$ in
[the affine-plane bridge, §4](HC4_AFFINE_PLANE_SCHUBERT_BRIDGE.md) does not
follow from the cited frozen-coordinate calculation.

There is a partial repair. Differentiating the actual sign-branch identities
excludes $p=q=a\ne0$ without any constant-motion determinant assumption.
The other sign $p=q=-a\ne0$ survives the next finite prolongation and has an
explicit compatible rational connection jet with $d(pq)\ne0$. The subsequent
[global polynomial-leaf proof](HC4_NEGATIVE_MOTION_POLYNOMIAL_OBSTRUCTION.md)
excludes its realization by the required polynomial endomorphism.

The exact correction and positive-sign closure are registered as **HC4MRA1**.
The audit temporarily made the full master reduction **HC4MR1** and the
equivalence **HC4MR2** partial. **HC4MRA2** now restores them using the
replacement global polynomial argument. The original constant-motion
inference remains withdrawn. This audit does not disprove either proposed theorem, construct
an HC4 counterexample, or prove that the surviving finite jet integrates.
The lower-rank results and the implication
$\mathrm{PHC4}\Rightarrow\mathrm{JC2}$ are not affected by this gap.

## 2. Exact transport from the adapted frame

At a generic point freeze an adapted frame $(e_1,e_2,e_3,e_4)$ with
$S$ the anti-diagonal matrix with entries one and $N=J_4(0)$, so
$T=SN$. Take the collapsed coordinate direction $z=e_1$ and passive
directions $(q_1,q_2,q_3)=(e_2,e_3,e_4)$. Then

\[
 T_{qq}=M=M^{-1}=
 \begin{pmatrix}0&0&1\\0&1&0\\1&0&0\end{pmatrix}.
\]

Use the same connection scalars as the affine-plane bridge:

\[
 a=\Gamma^2_{3,1}=\Gamma^3_{4,1},\quad
 b=\Gamma^2_{4,1},\quad
 p=\Gamma^4_{3,3},\quad q=\Gamma^4_{4,2},\quad
 s=\Gamma^4_{4,3}.
\]

The projective derivative of $e_1$ in the passive directions is

\[
 K=\begin{pmatrix}0&a&b\\0&0&a\\0&0&0\end{pmatrix}.
\]

Since $Te_1=0$ and the ambient derivative of $T$ is symmetric,
the frozen third-derivative matrix $A_{zqq}$ is $-K^TT_{qq}$.
Changing to the Legendre coordinates $\pi=A_q$ gives the graph Hessian

\[
 B=\operatorname{Hess}_{\pi}\mathcal H
   =M(-K^TT_{qq})M
   =\begin{pmatrix}-b&-a&0\\-a&0&0\\0&0&0\end{pmatrix}.
 \tag{2.1}
\]

Thus the Gauss-rank-two open is $a\ne0$. The selected line $Se_1$ is
represented in passive target coordinates by $u=(0,0,1)^T$ at the point.
Its projective derivative, expressed in the $\pi$ directions and modulo
the third target coordinate, is

\[
 U_{\mathrm{proj}}=
 \begin{pmatrix}-q&0\\-s&-p\end{pmatrix},
 \qquad \det U_{\mathrm{proj}}=pq.
 \tag{2.2}
\]

The zero last row used to represent this projective derivative need not
be the last row of the derivative of the actual, unrescaled vector $u$.
Only the induced projective map is used below.

To reach the frozen normal form used in
[HC4RSD72](HC4_FINAL_RANK_THREE_SMOOTH_CHART_OBSTRUCTION.md), make the
constant-at-the-point source change

\[
 q=C\widetilde q,\qquad z=c\widetilde z,\qquad
 C=\begin{pmatrix}1&b/(2a)&0\\0&1&0\\0&0&-a\end{pmatrix},
 \qquad c=-1/a.
 \tag{2.3}
\]

It has determinant $c\det C=1$. The corresponding Legendre and companion
coordinates transform by

\[
 \widetilde\pi=C^T\pi,\qquad \widetilde u=cC^Tu=u.
\]

Therefore

\[
 \widetilde B=cC^{-1}BC^{-T}
 =\begin{pmatrix}0&1&0\\1&0&0\\0&0&0\end{pmatrix},
 \qquad
 \widetilde U=cC^TUC^{-T}.
\]

Since the change is block diagonal relative to the projective line,
the two-dimensional projective determinant transforms as

\[
 \boxed{\det\widetilde U_{\mathrm{proj}}=\frac{pq}{a^2}.}
 \tag{2.4}
\]

The original scalar Hessian determinant is unchanged by (2.3). Hence even
accepting the normalized determinant identity of HC4RSD72 gives
$pq=a^2$ in the adapted unit-metric frame, rather than constancy of $pq$.
This relation was already supplied by the undifferentiated flatness
conditions $p=q$ and $q^2=a^2$.

Constant frame volume cannot remove this factor. Already when $b=0$,

\[
 D=\operatorname{diag}(-1/a,1,1,-a)
\]

satisfies $\det D=1$ and $D^TSD=S$, but $D^{-1}ND\ne N$.
Thus the sign-only centralizer of a frame normalizing both $S$ and $N$
does not apply to the comparison with the frozen normal frame. Freezing
a different matrix at each point does not make its entries constant across
points. No additional polynomial-unit argument for $a$ is established by
the existing proof.

## 3. A valid additional prolongation closes the positive sign

The original curvature elimination gives

\[
 a(p-q)=0,\qquad 4pa-3a^2-4aq+3q^2=0.
\]

On $a\ne0$ this splits into $p=q=a$ and $p=q=-a$. Each is an identity
on its generic branch, so its directional derivatives must also vanish.
For the positive branch adjoin, for every $i$,

\[
 e_i(p-a)=0,\qquad e_i(q-a)=0.
\]

Substitute $p=q=a$ in the 96 curvature equations and eliminate the
68 directional derivative variables. An exact constant left-kernel
combination of the resulting 104 equations is $-4a^2$.
Hence the positive sign is impossible in characteristic zero on $a\ne0$.
The complete coefficient vector is preserved in the replay artifact.
This proof uses the differentiated branch identities and does not use
$d(pq)=0$.

## 4. The negative sign survives this finite prolongation

The same operation for $p=q=-a$ gives the necessary relation

\[
 p_{12}+p_3+p_7=0,
 \tag{4.1}
\]

where $p_0,\ldots,p_{16}$ are the canonical nullspace parameters of the
published first-order constructor, with $a=p_4$, $p=p_0$, and $q=p_8$.
The artifact records the next quadratic compatibility equation obtained by
adjoining (4.1) and its four directional derivatives.

There is a rational finite jet satisfying all 96 curvature equations and
the 12 differentiated negative-branch relations. Take

\[
 p_0=p_8=-1,\qquad p_4=1,\qquad p_j=0
 \quad(j\notin\{0,4,8\}).
\]

Its only nonzero parameter derivatives are

\[
 e_1(p_2)=\tfrac32,\qquad
 e_2(p_0)=-\tfrac32,\qquad
 e_2(p_4)=\tfrac32,\qquad
 e_2(p_8)=-\tfrac32.
\]

Direct substitution gives

\[
 \bigl(e_i(pq)\bigr)_{i=1}^4=(0,3,0,0),\qquad
 \bigl(e_i(pq/a^2)\bigr)_{i=1}^4=(0,0,0,0).
\]

The old augmented system discards this compatible jet precisely because
it adds $d(pq)=0$. A finite jet does not establish all-order integrability,
an algebraic realization, polynomial Hessians, or a collision. Those
remain separate requirements.

At this audit stage the next target was the negative sign with its
actual differential compatibility equations. A proof that its invariant
scale $a$ is constant would recover the old contradiction, but requires
a new argument. The subsequent HC4MRA2 proof does rule out the branch by global
polynomiality, without asserting a local finite-jet contradiction.

## 5. Replay and evidence boundaries

Run from the repository root:

    .venv/bin/python scripts/verify_hc4_motion_frame_transport.py

The [checker](scripts/verify_hc4_motion_frame_transport.py) verifies the
transport by direct symbolic matrices. It reuses only the first-order
frame and curvature constructor from the old prolongation script for the
branch calculations; it does not execute that script's output writer or
import its extra determinant-constancy equations.

The pinned [artifact](artifacts/generated-results/hc4-motion-frame-transport-v1.json)
contains source hashes, the exact frame transition, a literal left-kernel
certificate for the positive sign, and the negative-sign jet. The default
requires byte-identical reproduction; the option --write deliberately
regenerates it.

The original
[prolongation artifact](artifacts/generated-results/hc4_affine_plane_prolongation.json)
is preserved. Its saturated unit ideal remains correct for the
**augmented system that assumes $d(pq)=0$**. It no longer certifies closure
of every geometric input in the final HC4 branch.

The v1 audit artifact retains its historical `status_effect` field from the
intermediate partial stage. Current theorem status is in MATH_STATUS.json
and the HC4MRA2 proof; that field is not a current status authority.

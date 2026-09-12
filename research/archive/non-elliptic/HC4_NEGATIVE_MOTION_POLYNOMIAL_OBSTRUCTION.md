# Polynomiality closes the negative HC4 motion sign

## 1. The repaired step

The [motion-frame audit](HC4_MOTION_FRAME_TRANSPORT_AUDIT.md) showed that the
old maximal-motion argument incorrectly treated $pq$ as constant. The
correct frozen determinant is $pq/a^2$. It also excluded the positive sign
$p=q=a\ne0$ by differentiating the actual branch identities, while leaving
$p=q=-a\ne0$ compatible with a finite connection jet.

The negative sign is excluded by **global polynomiality of $N=S^{-1}T$**,
rather than by a contradiction in that finite jet.

> **Theorem HC4MRA2 (negative-motion polynomial obstruction).** Let
> $S=\operatorname{Hess}\psi$ and $T=\operatorname{Hess}A$ be polynomial
> Hessians in four variables over a characteristic-zero field, with
> $\det S\in k^\times$ and $\det(S+sT)=\det S$. On the regular $J_4(0)$
> locus assume the previously established kernel quasi-translation and
> Frobenius Jordan flag used in the relative-pencil reduction.
> The negative maximal-motion branch $p=q=-a\ne0$ in its normalized
> moving-frame system is empty.

Together with the positive-sign certificate in HC4MRA1, this repairs the
maximal-motion step of
[HC4MR1](HC4_RELATIVE_NILPOTENT_MASTER_REDUCTION.md). It does not prove
unrestricted HC4 or JC2. The other reductions in the master proof map remain
separate written prerequisites, and no independent end-to-end or complete
formal verification is claimed.

## 2. Justifying the normalized kernel frame

Let $k$ be a polynomial quasi-translation field spanning $\ker T$ on the
regular locus, as in the preceding singular-Hessian reduction. Thus
$D_k k=0$. The maps $x\mapsto x+tk(x)$ and $x\mapsto x-tk(x)$ are inverse
polynomial automorphisms over the parameter ring. Their Jacobian
determinant is one, hence $\operatorname{div}k=0$.

Since $T$ is symmetric of rank three, write

\[
 \operatorname{adj}T=\sigma kk^T
\]

over the rational function field, with $\sigma\ne0$. The
[Hessian Piola identity](HC4_RANK_THREE_COFACTOR_FLAG.md) gives row divergence
zero. Expanding it and using $D_k k=0$ and $\operatorname{div}k=0$ gives

\[
 0=\operatorname{div}(\sigma kk^T)=(D_k\sigma)k.
\]

Thus $D_k\sigma=0$. The inverse-pencil identity gives

\[
 N^3S^{-1}=-\frac{\operatorname{adj}T}{\delta},
 \qquad \delta=\det S.
 \tag{2.1}
\]

In an adapted frame with $N=J_4(0)$ and $S$ the unit anti-diagonal matrix,
the left side of (2.1) is $e_1e_1^T$. Consequently

\[
 e_1=fk,\qquad f^2=-\sigma/\delta.
\]

We may choose $f$ in an algebraic extension on the generic locus.
Characteristic zero and $D_k\sigma=0$ give $D_k f=0$, and therefore

\[
 \nabla_{e_1}e_1=0.
 \tag{2.2}
\]

This supplies the affine normalization used by the first-order constructor.
It does not assume that an arbitrary rescaling of a quasi-translation is
again a quasi-translation. Likewise the adapted frame has constant volume:
if $E$ is its matrix, then $(\det E)^2\det S=1$, so
$e_i(\det E)=0$. These justify the kernel and trace conditions of the
published frame calculation.

## 3. Exact identities on the negative sign

The frame computation gives

\[
 a(p-q)=0,\qquad
 4pa-3a^2-4aq+3q^2=0.
\]

On the negative component $p=q=-a\ne0$, differentiating the branch
identities gives the additional relation
$p_{12}+p_3+p_7=0$ in the constructor's nullspace parameters.
The following identities hold there:

\[
 \nabla_{e_1}e_1=0,\qquad
 \nabla_{e_1}e_2=-a e_1,\qquad
 \nabla_{e_2}e_1=\frac52a e_1,
 \tag{3.1}
\]

\[
 \nabla_{e_2}e_2=\eta e_1+\frac12a e_2,\qquad
 e_1(a)=0,\qquad e_2(a)=\frac32a^2,
 \tag{3.2}
\]

where $\eta$ is irrelevant to the argument.

The connection identities are exact substitutions in the first-order
solution. The two derivative identities have short literal certificates
in the 108-equation negative-branch system. Numbering its equations from
zero, $e_1(a)$ equals

\[
 \tfrac12F_0-F_5-\tfrac12F_{10},
\]

and $e_2(a)-\tfrac32a^2$ equals

\[
 -F_9-F_{14}-F_{16}-F_{26}.
\]

Thus (3.2) does not assume $d(pq)=0$ or extrapolate a finite calculation.
The checker verifies these as polynomial identities in all the free
connection and derivative parameters.

## 4. Restriction to one affine two-plane

The previously established middle distribution
$E_2=\ker N^2=\langle e_1,e_2\rangle$ has leaves open in affine two-planes.
Choose a generic leaf meeting $a\ne0$. The projective direction of $e_1$
is constant along the leaf by (3.1). Choose constant independent vectors
$v_1,v_2$ in its direction plane, with $v_1$ spanning that kernel line,
and affine coordinates $(s,t)$ in those directions.

On a suitable local algebraic extension write

\[
 e_1=f(s,t)v_1,\qquad
 e_2=h(s,t)v_1+v(s,t)v_2,
 \qquad fv\ne0.
 \tag{4.1}
\]

The first two identities of (3.1) imply $f_s=v_s=0$. Hence $f$ and $v$
depend only on $t$. The transverse component of the fourth connection
identity and the third identity of (3.1) give

\[
 v'=\frac12a,\qquad
 vf'=\frac52af,
 \tag{4.2}
\]

where a prime denotes the ordinary affine derivative in $t$.
Also $a_s=0$, and the last equation of (3.2) gives

\[
 va'=\frac32a^2.
 \tag{4.3}
\]

Now use the polynomial endomorphism, not the moving frame, to define the
scalar that matters. Since $Nv_1=0$ and $Ne_2=e_1$,

\[
 Nv_2=n(t)v_1,\qquad n=f/v.
 \tag{4.4}
\]

The ambient matrix $N=S^{-1}T$ is polynomial because $\det S$ is a
nonzero constant. Restricting $N$ to the **whole affine plane** therefore
shows that the scalar coefficient in (4.4) is a polynomial in $(s,t)$.
The identities $f_s=v_s=0$ make it a polynomial in $t$ alone. It is nonzero
because $N|_{E_2}$ has rank one on the regular locus.

Neither the chosen frame nor the regular locus needs to extend across
the whole plane. Polynomiality of the ambient matrix is enough: the
matrix-coefficient identities on the nonempty leaf open extend as
polynomial identities to its affine-plane closure. Rank drops at omitted
points do not permit a pole in $n$.

## 5. A polynomial ODE contradiction

Equations (4.2)--(4.3) give

\[
 n'=\frac{2a}{v}n,\qquad
 n''=\frac{6a^2}{v^2}n.
\]

Therefore

\[
 \boxed{2nn''-3(n')^2=0.}
 \tag{5.1}
\]

A nonconstant polynomial $n$ of degree $d\geq1$, with leading coefficient
$c\ne0$, makes the coefficient of $t^{2d-2}$ in the left side equal to

\[
 \bigl(2d(d-1)-3d^2\bigr)c^2=-d(d+2)c^2\ne0
\]

in characteristic zero. Thus every polynomial solution of (5.1) is constant.
Our $n$ is a nonzero constant, so its derivative formula forces $a=0$,
contradicting the selected branch. This proves the theorem.

The compatible finite jet in HC4MRA1 is not contradicted. Nonconstant local
solutions of (5.1) can have a double pole; polynomial truncations can match
arbitrarily many Taylor coefficients of such a function without solving the
ODE identically. The missing condition was global polynomiality on the
affine leaf, precisely the property supplied by $N$.

## 6. Reproduction and proof boundary

Run:

    .venv/bin/python scripts/verify_hc4_negative_motion_polynomial_obstruction.py

The [checker](scripts/verify_hc4_negative_motion_polynomial_obstruction.py)
verifies the literal curvature combinations, the restricted connection,
the cofactor-normalization matrices, the rational differential identity
(5.1), and its symbolic leading coefficient.
The pinned [artifact](artifacts/generated-results/hc4-negative-motion-polynomial-obstruction-v1.json)
records exact source hashes and those identities. The default requires
byte-identical reproduction; --write deliberately regenerates it.

The existence of the preceding quasi-translation and Frobenius flag, the
affine-leaf restriction, and the universal polynomial-degree argument are
explicit written steps. The checker is not an independent replay of the
entire HC4 proof tree. The older augmented-system certificate and the
transport-audit finite jet remain useful regressions and are preserved.

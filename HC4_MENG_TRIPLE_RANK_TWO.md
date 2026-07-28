# Rank-two sextic reduction with simultaneous cubic and quartic layers

## Status

Let

\[
 \psi=q_2+h_3+h_4+h_6
\]

be a collision-normalized Meng potential in four variables, with \(q_2\)
nondegenerate and the other forms arbitrary homogeneous forms of degrees
three, four, and six.

> **Theorem `HC4T21`.**  If
> \(\operatorname{Hess}(h_6)\) has generic rank two and
> \(\det\operatorname{Hess}(\psi)\) is a nonzero constant, then
> \(\nabla\psi\) cannot identify a nonzero antipodal pair.

There is no support restriction.  Together with `HC4T31`, this leaves only
sextic Hessian rank one in the simultaneous cubic--quartic--sextic
coordinate chart; rank zero is `HC4CQ1`.

## 1. The first three faces

After scalar extension, Gordan--Noether makes the two-dimensional kernel
\(U=(t,m)\) of \(\operatorname{Hess}(h_6)\) a constant plane.  Write the
quotient variables as \(X=(x,y)\).  For

\[
 A=\operatorname{Hess}(h_3),\quad
 B=\operatorname{Hess}(h_4),\quad
 C=\operatorname{Hess}(h_6),
\]

the degree-twelve face of

\[
 \det(H_0+\lambda A+\lambda^2B+\lambda^4C)
\]

is

\[
 \det(\bar C)\det(B_U).                              \tag{1.1}
\]

Thus the binary Hessian \(B_U\) is singular.  Its degree-eleven successor
is

\[
 \det(\bar C)\,
 \operatorname{tr}\bigl(\operatorname{adj}(B_U)A_U\bigr). \tag{1.2}
\]

If \(B_U\) has a constant one-dimensional kernel, (1.2) makes that
direction isotropic for \(A_U\).  Hence

\[
 D_vh_6=0,\qquad D_v^2h_4=D_v^2h_3=0,
\]

and the common-direction reduction proved in
[`HC4_MENG_TRIPLE_RANK_THREE.md`](HC4_MENG_TRIPLE_RANK_THREE.md)
excludes the collision.

If \(B_U=0\), the degree-ten face is

\[
 \det(\bar C)\det(A_U).                              \tag{1.3}
\]

The two-variable singular-Hessian classification expresses the nonlinear
part of \(h_3\), over \(K(X)\), through one projective cone.  A moving cone
of degree \(\delta\ge1\) and highest nonlinear dual degree \(d\ge2\) would
need a scalar coefficient of source degree

\[
 3-d-d\delta<0.
\]

Therefore the cone is constant, again giving a common direction and the
same reduction.

## 2. The moving quartic cone

The only remaining alternative in the quartic cone-degree lemma `HC4E46`
is, modulo terms of dual degree at most one,

\[
 h_4=c(X^{\mathsf T}MU)^2,\qquad M\in\operatorname{GL}_2.
\]

Normalize \(M=I\) and put

\[
 L=xt+ym,\qquad k=(-y,x).
\]

Equation (1.2) becomes

\[
 k^{\mathsf T}\operatorname{Hess}_U(h_3)k=0.         \tag{2.1}
\]

Exact linear algebra on all twenty cubic monomials gives rank eight.
Its twelve-dimensional kernel is precisely

\[
 h_3=L(\alpha t+\beta m)
       +\text{terms of dual degree at most one}.     \tag{2.2}
\]

Thus the cubic is forced to align with the moving quartic cone.

This alignment still cannot cancel the later determinant obstruction.  The
dual-degree-four part of the spatial degree-eight face is exactly

\[
 48c^4\det(M)^2(X^{\mathsf T}MU)^4.                 \tag{2.3}
\]

The checker computes (2.3) with a general quadratic Hessian, source-only
sextic block, every compatible lower quartic block, every compatible lower
cubic block, and the aligned cubic (2.2) present.  It is nonzero, so the
moving branch is impossible.  This exhausts all three possibilities for
\(B_U\) and proves the theorem.

## Reproduction

Run:

```bash
.venv/bin/python \
  scripts/verify_hc4_meng_triple_rank_two_reduction.py
```

The checker also replays the exact quartic cone-degree theorem `HC4E46`
and the common-direction reduction `HC4T31`.

The external structural inputs are the
[Gordan--Noether classification](https://arxiv.org/abs/1501.05168),
the de Bondt--van den Essen
[singular-Hessian classification](https://www.sciencedirect.com/science/article/pii/S0021869304004867),
the Hessian conjecture in dimensions at most three, and
[Moh's plane degree bound](https://www.math.purdue.edu/~ttm/jacobian.pdf).

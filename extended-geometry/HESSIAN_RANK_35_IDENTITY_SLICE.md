# A rank-35 nonhomogeneous Hessian-nilpotent witness

## Result and scope

The identity output of the repository's rank-directed 22-variable cubic
collision gives a 21-variable nonhomogeneous nilpotent-Jacobian collision

\[
 F(x)=x+K(x),\qquad \deg K\leq3,\qquad
 \operatorname{rank}_{\mathbb Q(x)}JK=17.
\]

Its 42-variable cotangent potential

\[
 p(x,y)=y^TK(x)
\]

has exact generic Hessian rank

\[
 \boxed{\operatorname{rank}_{\mathbb Q(x,y)}\operatorname{Hess}p=35}.
\]

After the standard complex orthogonal change, this gives a nonhomogeneous
Hessian-nilpotent polynomial \(P\) of degrees two, three, and four whose
gradient map is noninjective.  Thus \(P\) is an HN Vanishing counterexample.

This is an exact rank result only in the **nonhomogeneous** HN class.
It does not change either homogeneous-quartic endpoint:

\[
 \rho_{\mathrm{HN},4}^{\mathrm{hom}}\leq37,\qquad
 n_{\mathrm{HN},4}^{\mathrm{hom}}\leq42.
\]

In particular, the result must not be entered as a rank-35 homogeneous
quartic witness.

## 1. The identity slice

Let

\[
 \widetilde F(q,s)=(q+\mathcal K(q,s),s)
\]

be the homogeneous 22-variable map stored in
[`hessian_rank_reduced_bcw_22_counterexample.json`](../artifacts/generated-results/hessian_rank_reduced_bcw_22_counterexample.json).
Its three rational collision points all have \(s=1\).  Restricting to that
slice gives

\[
 F(q)=q+K(q),\qquad K(q)=\mathcal K(q,1),                 \tag{1.1}
\]

and the three restricted points remain distinct with one common image.  Their
zeroth coordinates are \(0,1,-1\), so the same coordinate-separation
certificate proves noninvertibility.

The correction \(K\) has degrees one, two, and three.  Since the bottom row
of \(J\mathcal H\) is zero, \(JK\) is its upper-left block after \(s=1\).
The exact identity

\[
 (J\mathcal H)^{18}=0
\]

therefore specializes to

\[
 (JK)^{18}=0.                                             \tag{1.2}
\]

Consequently \(\det(I+tJK)=1\) for every scalar \(t\).

## 2. Exact ranks

The Hessian of the cotangent potential has the block form

\[
 \operatorname{Hess}p=
 \begin{pmatrix}
 A&(JK)^T\\
 JK&0
 \end{pmatrix},
 \qquad
 A=\sum_i y_i\operatorname{Hess}K_i.                       \tag{2.1}
\]

Three deterministic good-prime specializations give

\[
 \operatorname{rank}JK=17,\qquad
 \operatorname{rank}A=14,\qquad
 \operatorname{rank}\operatorname{Hess}p=35.              \tag{2.2}
\]

These are lower bounds over the rational function fields.  For the matching
upper bounds, Singular computes:

- four polynomial syzygy generators of \(JK\) which specialize to four
  independent kernel columns;
- thirteen polynomial syzygy generators of \(\operatorname{Hess}p\) which
  contain seven generically independent kernel columns.

Thus the ranks are at most \(21-4=17\) and \(42-7=35\), respectively.
Together with (2.2), both ranks are exact.

Equivalently, the cotangent kernel excess is one:

\[
 35=2\cdot17+1.                                            \tag{2.3}
\]

The circuit therefore already attains the rank-17/excess-one target after
identity slicing.  Its homogeneous source has rank 18 and the same excess
one, which explains the two-rank gap between the sliced value 35 and the
homogeneous value 37.

## 3. Hessian nilpotency and Vanishing failure

For a parameter \(t\), put

\[
 f_t(x,y)=x^Ty+t\,y^TK(x).
\]

The block determinant identity and (1.2) give

\[
 \det\operatorname{Hess}f_t
 =(-1)^{21}\det(I+tJK)^2=-1.                               \tag{3.1}
\]

Make the invertible complex-linear change

\[
 x=u+iv,\qquad y=u-iv
\]

and define

\[
 P(u,v)=\frac12(u-iv)^TK(u+iv).                            \tag{3.2}
\]

After normalizing by the value at \(t=0\), equation (3.1) becomes

\[
 \det(I+t\operatorname{Hess}P)=1.                          \tag{3.3}
\]

Hence every nonconstant characteristic coefficient of
\(\operatorname{Hess}P\) vanishes, so \(\operatorname{Hess}P\) is
nilpotent.  Congruence and the nonzero scalar in (3.2) preserve generic
rank, giving

\[
 \operatorname{rank}\operatorname{Hess}P=35.               \tag{3.4}
\]

The gradient of \(f_1\) is noninjective at the three points
\((q,0)\).  The same source and target changes therefore show that
\((u,v)+\nabla P(u,v)\) is noninjective.  The HN inverse/Vanishing
equivalence then gives

\[
 \Delta^mP^m=0\quad(m\geq1),\qquad
 \Delta^mP^{m+1}\ne0\quad\text{for infinitely many }m.      \tag{3.5}
\]

Thus (3.2) is a 42-variable, degree-at-most-four, nonhomogeneous HN
Vanishing counterexample of exact generic Hessian rank 35.

## 4. Reproduction

Run

```bash
.venv/bin/python scripts/verify_hessian_rank_35_identity_slice.py
python3 scripts/audit_hessian_rank_35_identity_slice_independent.py
```

The first command requires Singular.  It writes the explicit specialized
correction, collision, rank certificates, and good-prime profiles to
[`hessian_rank_35_identity_slice_counterexample.json`](../artifacts/generated-results/hessian_rank_35_identity_slice_counterexample.json).
The second command reconstructs the slice and the three modular block-rank
profiles using only the Python standard library.

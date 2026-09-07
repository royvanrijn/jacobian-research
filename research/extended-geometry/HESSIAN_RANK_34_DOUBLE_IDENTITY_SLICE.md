# A 40-variable rank-34 nonhomogeneous HN witness

## Result and scope

The rank-directed circuit has a second identity output after its homogenizing
coordinate is specialized.  Restricting both identity outputs gives a
20-variable nilpotent-Jacobian collision

\[
 F(x)=x+L(x),\qquad \deg L\leq3,\qquad
 \operatorname{rank}_{\mathbb Q(x)}JL=17.
\]

Its cotangent potential has the exact profile

\[
 p(x,y)=y^TL(x),\qquad
 \boxed{\operatorname{rank}_{\mathbb Q(x,y)}
        \operatorname{Hess}p=34},
\]

in 40 variables.  The cotangent kernel excess is zero:

\[
 34=2\operatorname{rank}JL.
\]

After the standard complex orthogonal change this is a nonhomogeneous
Hessian-nilpotent Vanishing counterexample of degrees two, three, and four.
It improves the repository's nonhomogeneous rank-35 witness while also
returning to the existing ordinary-Laplacian/HN dimension 40.

The homogeneous quartic frontiers remain unchanged:

\[
 \rho_{\mathrm{HN},4}^{\mathrm{hom}}\leq37,\qquad
 n_{\mathrm{HN},4}^{\mathrm{hom}}\leq42.
\]

## 1. The second identity output

Let \(X+K(X)\) be the 21-variable identity slice constructed in
[`HESSIAN_RANK_35_IDENTITY_SLICE.md`](HESSIAN_RANK_35_IDENTITY_SLICE.md).
Exact coefficient comparison gives

\[
 \boxed{K_9-3K_1-K_6=0}.                                  \tag{1.1}
\]

Therefore

\[
 \ell(X)=-3X_1-X_6+X_9
\]

is an identity output:

\[
 \ell(X+K(X))=\ell(X).                                     \tag{1.2}
\]

All three stored collision points satisfy \(\ell=0\).  Put

\[
 X_9=3X_1+X_6                                                \tag{1.3}
\]

and retain the other twenty coordinates.  Equation (1.2) makes this
hyperplane invariant.  The restricted map is \(X+L(X)\), and its three
rational collision points remain distinct.  Their zeroth coordinates are
\(0,1,-1\).

After a rational linear coordinate change, \(JK\) has block form

\[
 JK=
 \begin{pmatrix}
 JL&*\\
 0&0
 \end{pmatrix}
\quad\text{on }\ell=0.                                     \tag{1.4}
\]

The exact nilpotency of \(JK\) therefore implies nilpotency of \(JL\).

## 2. Exact rank and zero kernel excess

The cotangent Hessian is

\[
 M=
 \begin{pmatrix}
 A&(JL)^T\\
 JL&0
 \end{pmatrix},
\qquad
 A=\sum_i y_i\operatorname{Hess}L_i.                       \tag{2.1}
\]

Three deterministic good-prime specializations give

\[
 \operatorname{rank}JL=17,\qquad
 \operatorname{rank}A=12,\qquad
 \operatorname{rank}M=34.                                 \tag{2.2}
\]

Singular supplies the matching characteristic-zero upper bounds:

- three polynomial syzygy generators of \(JL\) specialize to three
  independent kernel columns;
- twelve polynomial syzygy generators of \(M\) contain six independent
  kernel columns.

An exact specialization attains ranks 17 and 34.  Hence both generic ranks
are exact.  In the block-rank formula

\[
 \operatorname{rank}M
 =2\operatorname{rank}JL+\operatorname{rank}(C^TAC),
\]

where the columns of \(C\) span \(\ker JL\), equations (2.2) give

\[
 \boxed{C^TAC=0}.                                          \tag{2.3}
\]

The coefficient equations for a constant kernel vector of \(M\) have full
rank 40.  Thus the 40-variable potential has no further constant linear
direction to remove.

## 3. HN and ordinary-Laplacian consequences

Define

\[
 P(u,v)=\frac12(u-iv)^TL(u+iv).                            \tag{3.1}
\]

As in the first identity slice, nilpotency of \(JL\) gives

\[
 \det(I+t\operatorname{Hess}P)=1.                          \tag{3.2}
\]

Thus \(\operatorname{Hess}P\) is nilpotent.  The three points
\((x,0)\) in the cotangent construction have one gradient image, so
\(I+\nabla P\) is noninjective.  Consequently

\[
 \Delta^mP^m=0\quad(m\geq1),\qquad
 \Delta^mP^{m+1}\ne0\quad\text{for infinitely many }m.      \tag{3.3}
\]

The same 20-variable map supplies the contraction polynomial linear in
twenty dual variables.  The coordinate-zero separation gives an ordinary
Laplacian GVC counterexample in 40 variables.  This matches, rather than
lowers, the repository's current ordinary-Laplacian dimension bound.

## 4. Search boundary

The rank-34 witness was exposed after two bounded diagnostics:

1. the full frozen width-64 circuit census has 140 terminals and a unique
   best first-slice profile \((35,17,1,21)\);
2. 64 neutral low-degree perturbations around that terminal, continued for
   up to ten cleanup steps, give 414 terminals and no first-slice Hessian
   rank below 35.

These searches are not lower bounds.  Their records are
[`identity_slice_hessian_rank_search.json`](../artifacts/generated-results/identity_slice_hessian_rank_search.json)
and
[`identity_slice_local_perturbation_search.json`](../artifacts/generated-results/identity_slice_local_perturbation_search.json).
The exact improvement instead comes from recognizing the constant-kernel
direction as the output identity (1.1).

## 5. Reproduction

Run

```bash
.venv/bin/python scripts/verify_hessian_rank_34_double_identity_slice.py
python3 scripts/audit_hessian_rank_34_double_identity_slice_independent.py
```

The first command requires Singular.  It writes the explicit map, collision,
rank certificates, and HN metadata to
[`hessian_rank_34_double_identity_slice_counterexample.json`](../artifacts/generated-results/hessian_rank_34_double_identity_slice_counterexample.json).
The second command independently reconstructs the relation, elimination,
collision, and modular rank profiles using only the Python standard library.


# Weighted invariant-coordinate Jacobian reduction

This note isolates the reusable differential identity behind the weighted
three-variable constructions.  It converts a three-variable Jacobian equation
into a two-variable determinant in torus invariants.

Work over a characteristic-zero field.  Give the source coordinates
`(x,y,z)` weights

\[
 (1,-1,-k),\qquad k\ge1,
\]

and put

\[
 u=xy,\qquad v=x^kz.
\]

On `x!=0`, every rational function of weight `w` has the form

\[
 x^w A(u,v).
\]

For polynomial maps, the apparent negative powers of `x` must of course
cancel after substituting `u=xy` and `v=x^kz`.

## Determinant lemma

Let

\[
 F_i=x^{w_i}A_i(u,v),\qquad i=1,2,3.
\]

Then

\[
\boxed{
 \det{\partial(F_1,F_2,F_3)\over\partial(x,y,z)}
 =
 x^{w_1+w_2+w_3+k}
 \det
 \begin{pmatrix}
 w_1A_1&(A_1)_u&(A_1)_v\\
 w_2A_2&(A_2)_u&(A_2)_v\\
 w_3A_3&(A_3)_u&(A_3)_v
 \end{pmatrix}.}                                      \tag{1}
\]

### Proof

Use the intermediate coordinates `(x,u,v)`.  Their Jacobian is

\[
 {\partial(x,u,v)\over\partial(x,y,z)}
 =
 \begin{pmatrix}
 1&0&0\\
 y&x&0\\
 kx^{k-1}z&0&x^k
 \end{pmatrix},
 \qquad
 \det=x^{k+1}.                                       \tag{2}
\]

On the other side,

\[
 {\partial F_i\over\partial(x,u,v)}
 =
 x^{w_i}\left(x^{-1}w_iA_i,(A_i)_u,(A_i)_v\right).
\]

Factoring `x^(w_i)` from the three rows and `x^(-1)` from
the first column gives

\[
 \det{\partial(F_1,F_2,F_3)\over\partial(x,u,v)}
 =x^{w_1+w_2+w_3-1}\det M,                           \tag{3}
\]

where `M` is the matrix in (1).  Multiplying (2) and (3) proves the formula.
Since both sides are rational functions agreeing on `x!=0`, the identity
extends polynomially whenever all `F_i` are polynomials.

## Keller consequence

The source volume form has torus weight `-k`.  Thus a nonzero constant
Jacobian is possible only when

\[
 w_1+w_2+w_3=-k.
\]

Under this balance condition, (1) becomes the two-variable equation

\[
\boxed{
 \det
 \begin{pmatrix}
 w_1A_1&(A_1)_u&(A_1)_v\\
 w_2A_2&(A_2)_u&(A_2)_v\\
 w_3A_3&(A_3)_u&(A_3)_v
 \end{pmatrix}
 =\text{constant}.}                                  \tag{4}
\]

For the foundational weights

\[
 k=2,\qquad(w_1,w_2,w_3)=(-2,-1,1),
\]

this is exactly

\[
 \det
 \begin{pmatrix}
 -2A&A_u&A_v\\
 -B&B_u&B_v\\
 C&C_u&C_v
 \end{pmatrix}.
\]

When `C` is nonconstant affine linear, this three-row determinant admits a
further rearrangement.  After an oriented affine change on the invariant
plane and the substitutions `P=C^2A`, `Q=CB`, it becomes the Poisson-square
equation `[P/2,Q]=C^2` at the foundational determinant.  The proof,
orientation factor, explicit pair, and coefficient cascade are in the
[weighted tangent-suspension note](WEIGHTED_TANGENT_SUSPENSION.md).

## What the reduction does and does not prove

Equation (4) is useful for coefficient schemes because all coefficient
comparison takes place in `k[u,v]` rather than `k[x,y,z]`.  It can be used to:

- construct the exact Keller ideal for a declared weighted support;
- compute its Zariski tangent space and quadratic obstructions;
- separate diagonal gauge directions from genuine coefficient deformations;
- test higher-degree weighted seeds with the same source torus.

It does not classify weighted Keller maps without a fixed support and
polynomiality conditions.  In particular, the claimed sixteen-monomial
cubic coefficient scheme cannot be reconstructed from (4) without the exact
three supports, linear normalization, and coefficient labels.

The exact checker is
[`verify_weighted_invariant_jacobian_reduction.py`](../scripts/verify_weighted_invariant_jacobian_reduction.py).

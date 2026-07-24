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

## Quotient form and the contracted divisor

There is a useful quotient interpretation of (4).  Give the target
coordinates `(a,b,c)` weights

\[
 (-k,-1,1)
\]

and let `G=(G_1,G_2,G_3)` be an equivariant polynomial map.  On `x!=0` write

\[
 G_1=x^{-k}A(u,v),\qquad
 G_2=x^{-1}B(u,v),\qquad
 G_3=x\Lambda(u,v).
\tag{5}
\]

The source and target invariant rings are

\[
 k[x,y,z]^{\mathbb G_m}=k[u,v],
 \qquad
 k[a,b,c]^{\mathbb G_m}=k[P,Q],
\]

where

\[
 u=xy,\quad v=x^kz,\qquad P=bc,\quad Q=ac^k.
\tag{6}
\]

Thus `G` descends to the explicit quotient morphism

\[
 \bar G:\mathbb A^2_{u,v}\longrightarrow\mathbb A^2_{P,Q},
 \qquad
 (P,Q)=\bigl(B\Lambda,A\Lambda^k\bigr).
\tag{7}
\]

With the orientation

\[
 \operatorname{Jac}_{u,v}(P,Q)=P_uQ_v-P_vQ_u,
\]

the determinant in (4) satisfies

\[
 \boxed{\operatorname{Jac}_{u,v}(P,Q)
       =-\Lambda^k\det JG.}
\tag{8}
\]

Indeed, direct expansion gives

\[
 \begin{aligned}
 \operatorname{Jac}_{u,v}(P,Q)
 =-\Lambda^k\bigl(&\Lambda\operatorname{Jac}(A,B)
 +B\operatorname{Jac}(A,\Lambda)\\
 &-kA\operatorname{Jac}(B,\Lambda)\bigr),
 \end{aligned}
\]

and the expression in parentheses is exactly the three-row determinant in
(4).  Consequently, if `det JG=K` is a nonzero constant, then

\[
 \operatorname{Jac}_{u,v}(P,Q)=-K\Lambda^k.
\tag{9}
\]

The divisor `D_Lambda=V(Lambda)` is mapped by `bar G` to `(P,Q)=(0,0)`.
When `Lambda` is reduced, (9) says scheme-theoretically that the quotient
Jacobian vanishes to order exactly `k` along every component of
`D_Lambda`.  For `k=2`, this recovers Shaska's order-two quotient-Jacobian
theorem, now with the quotient map and its orientation written explicitly.

The orbit quotient has a separate contraction which should not be confused
with `D_Lambda`.  The coarse quotient maps are

\[
 \pi_X(x,y,z)=(xy,x^kz),\qquad
 \pi_Y(a,b,c)=(bc,ac^k).
\tag{10}
\]

They contract the planes `x=0` and `c=0`, respectively, to their quotient
origins.  By contrast, `D_Lambda` is a divisor *on the source quotient* and
is contracted by the descended plane map (7).

## Explicit refinement for every admissible weighted seed

Specialize to `k=2` and let `H` be any admissible weighted seed.  To avoid a
collision with the target coordinate `v`, write the source invariants as

\[
 u=xy,\qquad v=x^2z.
\]

With the seed constants `a_0,b_0,c`, put

\[
 \gamma=1+a_0u+b_0v,\qquad W=(1+u)\gamma,
\tag{11}
\]

where `b_0c!=0`.  For the associated weighted Keller map `G_H=(A,B,C)`,
the two target invariants and Shaska's descended morphism are

\[
 \boxed{
 \begin{aligned}
 P=BC&=H'(W)+c\gamma,\\
 Q=AC^2&=\frac{W(H'(W)+c\gamma)-H(W)}{c}.
 \end{aligned}}
\tag{12}
\]

Equivalently, `bar G_H` has the family-level factorization

\[
 \mathbb A^2_{u,v}
 \xrightarrow{\ \theta_H\ }
 \mathbb A^2_{W,\gamma}
 \xrightarrow{\ \Psi_H\ }
 \mathbb A^2_{P,Q},
\tag{13}
\]

where

\[
 \theta_H(u,v)=((1+u)\gamma,\gamma),
\qquad
 \Psi_H(W,\gamma)=
 \left(H'(W)+c\gamma,
 \frac{W(H'(W)+c\gamma)-H(W)}{c}\right).
\]

The two Jacobians are

\[
 \det D\theta_H=b_0\gamma,\qquad
 \det D\Psi_H=-c\gamma.
\tag{14}
\]

Therefore

\[
 \boxed{\operatorname{Jac}_{u,v}(P,Q)
       =-b_0c\,\gamma^2
       =-\gamma^2\det DG_H.}
\tag{15}
\]

This identifies the two copies in the order-two vanishing: one is the
Jacobian of the contraction `theta_H`, and one is the critical factor of the
tangent-incidence map `Psi_H`.  The contracted divisor is the reduced affine
line

\[
 D_H=V(\gamma)=V(1+a_0u+b_0v)\cong\mathbb A^1,
\tag{16}
\]

and `theta_H(D_H)=(0,0)`, hence `bar G_H(D_H)=(0,0)`.  Thus the quotient
ramification divisor is `2D_H`, independently of the degree and coefficients
of the admissible seed.

## Stacky stabilizer strata

The coarse planes in (10) suppress stabilizer data.  For the source action
of weights `(1,-1,-k)`, the quotient stack
`[A^3_(x,y,z)/G_m]` has the following stabilizers:

\[
\begin{array}{c|c}
\text{stratum}&\text{stabilizer}\\ \hline
x\ne0\ \text{or}\ (x=0,\ y\ne0)&1\\
x=y=0,\ z\ne0&\mu_k\\
x=y=z=0&\mathbb G_m.
\end{array}
\tag{17}
\]

For the target action of weights `(-k,-1,1)`, the corresponding table is

\[
\begin{array}{c|c}
\text{stratum}&\text{stabilizer}\\ \hline
c\ne0\ \text{or}\ (c=0,\ b\ne0)&1\\
b=c=0,\ a\ne0&\mu_k\\
a=b=c=0&\mathbb G_m.
\end{array}
\tag{18}
\]

For an admissible seed (`k=2`), these tables are compatible stratum by
stratum.  On the contracted source plane `x=0`, weighted polynomiality gives

\[
 C=0,\qquad
 B=-\frac{c}{2+\kappa}\,y,\qquad
 A=b_0(2+\kappa)z+d_Hy^2,
\tag{19}
\]

for a seed-dependent constant `d_H`, where
`\kappa=H''(1)/c!= -2`.  Hence the punctured `mu_2`-line
`x=y=0,z!=0` maps to the punctured target `mu_2`-line
`b=c=0,a!=0`; the origins map to each other; and the remaining points of the
contracted planes have trivial stabilizer.  The divisor `D_H` in (16) is
different: its generic lift lies in `x!=0`, so its generic stabilizer is
trivial even though the descended morphism contracts it.

## Relation to Shaska's grading theorem

Shaska established that positive-weight equivariant Keller maps are
automorphisms in all dimensions, that every equivariant plane Keller map is
an automorphism for every signature, and that in the foundational mixed
signature the Keller condition descends to an order-two quotient-Jacobian
identity.  The results above do not compete with those signature theorems:
they identify the quotient, contraction, multiplicity split, and stabilizer
strata uniformly throughout the admissible weighted-seed family.

The appropriate division of credit is:

> Shaska established the general signature and quotient constraints.  We
> compute the quotient, ramification and normalization data explicitly for
> the weighted Keller family.

See T. Shaska,
[*Graded Keller maps and the Jacobian
Conjecture*](https://arxiv.org/abs/2607.20210), arXiv:2607.20210 (2026).

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

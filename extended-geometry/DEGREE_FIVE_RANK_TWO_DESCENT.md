# Worked example: rank-two symplectic descent in degree five

This note preserves the first explicit rank-two calculation as a worked
degree-five example.  It proves directly that the original one-parameter
stable-moduli family admits a four-dimensional Poisson completion, descending
that family from the universal six-dimensional cotangent lift to two
canonical pairs.  The uniform theorem for every degree `N>=5` is now the
canonical source: [Rank-two symplectic descent](RANK_TWO_SYMPLECTIC_DESCENT.md).

There is a normalization point to settle first.  In the verified family,

\[
 \gamma=1-\frac87xy+x^2z,
 \qquad C_\lambda=x\gamma,                            \tag{1}
\]

so the constant weighted parameter is `a_0=-8/7`, not `-7/8`.  The seed
parameter `lambda` occurs only in `A_lambda,B_lambda`; the third component
`C_lambda` is fixed.

Work over the parameter open `Lambda` of the
[degree-five theorem](DEGREE_FIVE_STABLE_MODULI.md).  In particular
`lambda-1` is a unit.

## 1. A common adapted coordinate system

Let `(x,q,p,zeta)` carry the canonical bracket

\[
 \{p,x\}=\{\zeta,q\}=1.
\]

Put

\[
\begin{aligned}
X&=x,\qquad Q=q,\\
Z&=3x^2p+(2-6xq)\zeta,\\
E&=\frac{1+3xq}{2}p-3q^2\zeta.
\end{aligned}                                        \tag{2}
\]

This is a polynomial coordinate system.  Its inverse is

\[
 p=-6EXQ+2E+3ZQ^2,qquad
 \zeta=-3EX^2+\frac32ZXQ+\frac12Z.                   \tag{3}
\]

The induced brackets are

\[
\begin{gathered}
\{X,Q\}=0,qquad \{X,Z\}=-3X^2,qquad
\{Q,Z\}=6XQ-2,\\
\{E,X\}=\frac{1+3XQ}{2},qquad
\{E,Q\}=-3Q^2,qquad
\{E,Z\}=\frac92QZ.                                   \tag{4}
\end{gathered}
\]

For a shear parameter `k`, define

\[
 W=Z+kQ^2,qquad Y=Q-\frac{XW}{3}.                    \tag{5}
\]

Use the fixed diagonal source change

\[
 (x_{\rm seed},y_{\rm seed},z_{\rm seed})
 =\left(X,\frac{21}{16}Y,-\frac12W\right).           \tag{6}
\]

For the explicit map `F_lambda=(A_lambda,B_lambda,C_lambda)`, put

\[
\begin{aligned}
S_{\lambda,k}&={160\over7(\lambda-1)}
 A_\lambda\left(X,{21\over16}Y,-{W\over2}\right),\\
T_{\lambda,k}&=
 B_\lambda\left(X,{21\over16}Y,-{W\over2}\right),\\
R&=2C_\lambda\left(X,{21\over16}Y,-{W\over2}\right).
\end{aligned}                                        \tag{7}
\]

The last expression is independent of both `lambda` and `k`:

\[
 \boxed{R=2X-3X^2Q.}                                 \tag{8}
\]

Direct chain-rule calculation gives

\[
 \det{\partial(S_{\lambda,k},T_{\lambda,k},R)
 \over\partial(X,Q,Z)}=-1.                           \tag{9}
\]

For the quotient bracket in (4), equation (9) is equivalently

\[
 \{S_{\lambda,k},T_{\lambda,k}\}=1,qquad
 \{R,S_{\lambda,k}\}=\{R,T_{\lambda,k}\}=0.        \tag{10}
\]

## 2. Relative flux and the complete residue obstruction

Let `w_(lambda,k)` be the unique horizontal derivation satisfying

\[
 w(S_{\lambda,k})=w(T_{\lambda,k})=0,qquad w(R)=1,
\]

and let `w_E={E,-}` on `Q(lambda)[X,Q,Z]`.  A completion

\[
 D=E+f(X,Q,Z)                                        \tag{11}
\]

has zero mixed brackets with `S,T` exactly when

\[
 w_{\lambda,k}-w_E=\{f,-\}.                          \tag{12}
\]

The following finite integration constructs the forced localized
Hamiltonian.  First set

\[
 f_0=\int
 {\bigl(w_{\lambda,k}-w_E\bigr)(X)\over3X^2}\,dZ,
                                                               \tag{13}
\]

with zero `Z`-constant, and subtract `{f_0,-}`.  The first two residual
components vanish.  On the third, introduce

\[
 v={1\over X},\qquad \rho=R,qquad
 Q={v(2-\rho v)\over3}.                               \tag{14}
\]

The remaining derivation becomes `3 partial_v`, so one more ordinary
antiderivative produces `h` and

\[
 f_{\lambda,k}=f_0+h.                                 \tag{15}
\]

The exact calculation verifies (12) over the localization `Q(lambda,k)[X,
X^-1,Q,Z]`.  More importantly, its **entire** negative-`X` principal part is

\[
 \mathcal O_{\lambda,k}left(
 {1\over26460X^4}+{Q\over4410X^3}
 +{Q^2\over1176X^2}+{Q^3\over392X}
 \right),                                             \tag{16}
\]

where

\[
 \mathcal O_{\lambda,k}=
 {196(\lambda-1)^2k+27(57\lambda^2-138\lambda+73)
 \over(\lambda-1)^2}.                                \tag{17}
\]

There are no further negative powers.  Consequently the unique pole-free
quadratic shear is

\[
 \boxed{
 \kappa_\lambda=-{27(57\lambda^2-138\lambda+73)
 \over196(\lambda-1)^2}.}                             \tag{18}
\]

For `k=kappa_lambda`, the finite formula (13)--(15) is a polynomial

\[
 f_\lambda\in\mathbb Q[\lambda,(\lambda-1)^{-1}][X,Q,Z].
                                                               \tag{19}
\]

Equations (13)--(15), with the zero integration constants specified there,
are an exact symbolic formula for `f_lambda`.  Its expanded numerator is
large and carries no additional mathematical information; the verifier
constructs it exactly and checks the three components of (12).

## 3. Four-dimensional regression

### Exact specialization — degree-five rank-two descent

For every `lambda in Lambda`, set

\[
 S_\lambda=S_{\lambda,\kappa_\lambda},\qquad
 T_\lambda=T_{\lambda,\kappa_\lambda},\qquad
 D_\lambda=E+f_\lambda
\]

and define

\[
 G_\lambda=(R,T_\lambda,D_\lambda,S_\lambda):
 \mathbb A^4\longrightarrow\mathbb A^4.              \tag{20}
\]

Then

\[
 \{D_\lambda,R\}=1,qquad
 \{S_\lambda,T_\lambda\}=1,                         \tag{21}
\]

and all four mixed brackets vanish.  Hence `G_lambda` preserves the canonical
symplectic form and has Jacobian one.

Moreover, the coordinate change

\[
 (X,Q,Z,E)\longmapsto(X,Y,W,E+f_\lambda)              \tag{22}
\]

is a polynomial automorphism: its inverse first uses

\[
 Q=Y+XW/3,qquad Z=W-\kappa_\lambda Q^2,
\]

and then subtracts `f_lambda`.  Equations (6)--(7), followed by a diagonal
target scaling and permutation, therefore give

\[
 \boxed{G_\lambda\sim_{\rm LR}
 F_\lambda\times\operatorname{id}_{\mathbb A^1}.}    \tag{23}
\]

Thus `G_lambda` has generic degree five and the same fiber schemes as
`F_lambda`.  Since affine four-space has trivial first algebraic de Rham
cohomology, preservation of the symplectic form also makes the lift exact.

### Corollary — uncountable rank-two symplectic moduli

The family contains uncountably many pairwise stably polynomially
left--right inequivalent exact symplectic maps

\[
 \boxed{G_\lambda:\mathbb A^4\longrightarrow\mathbb A^4}
\]

of generic degree five.  Indeed, a stable equivalence between two `G` maps,
together with (23), would give a stable equivalence between the corresponding
`F_lambda`; the degree-five Hessian-root invariant permits at most six
parameters in one stable class.

## Verification and scope

The exact certificate is
[`verify_degree_five_rank_two_descent.py`](../scripts/verify_degree_five_rank_two_descent.py).
It verifies the fixed coordinate, determinant normalization, quotient
brackets, all three relative-Hamiltonian equations, the complete principal
part, uniqueness of (18), polynomiality, and the source automorphism.

This is an internal regression with no recorded external specialist review.
Its role is an explicit worked example for the original degree-five line, not
a separately maintained theorem source.  The
[all-degree theorem](RANK_TWO_SYMPLECTIC_DESCENT.md) proves
completion for every normalized admissible seed in every degree `N>=5`, on
both adapted charts.  Its coarse decorated-normalization corollary transfers
`N-3` stable dimensions to `A^4` with generic degree `N-2`; in degree five
this gives the complete two-dimensional seed surface, with the present line
as the `kappa=-9` specialization.

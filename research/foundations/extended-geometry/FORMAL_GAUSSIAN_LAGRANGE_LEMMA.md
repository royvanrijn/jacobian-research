# Formal Gaussian--Lagrange lemma with constant terms

This note isolates the formal identity used by the weighted Gaussian bridge.
Its point is the case `Phi(0) != 0`: the fixed point then moves away from the
zero section, so a proof stated only in the `z`-adic topology is insufficient.
Everything below is coefficientwise and formal; no analytic expectation of an
exponential is used.

## 1. Statement and conventions

Let `k` be a characteristic-zero field, let

\[
 \Phi=(\Phi_1,\ldots,\Phi_r)\in k[z_1,\ldots,z_r]^r,
 \qquad A\in k[[z_1,\ldots,z_r]],
\]

and impose no condition on `Phi(0)`.  Write

\[
 J\Phi=\left(\frac{\partial\Phi_i}{\partial z_j}\right)_{1\le i,j\le r}.
\]

There is a unique fixed point

\[
 g(u)=u\Phi(g(u))\in u k[[u]]^r.                       \tag{1.1}
\]

Algebraically, let `E` be the `k`-linear Wick functional determined by

\[
 \mathbb E(Z^\beta W^\alpha)
 =\delta_{\alpha\beta}\alpha!.
\]

For `k=C`, this is realized by independent circular complex Gaussian pairs
normalized by

\[
 Z_j=\frac{X_j+iY_j}{\sqrt2},\qquad
 W_j=\frac{X_j-iY_j}{\sqrt2},                         \tag{1.2}
\]

where all `X_j,Y_j` are independent real `N(0,1)` variables.  For a formal
series in `u`, expectation means coefficientwise Wick contraction.

### Theorem 1.1 — constant-term Gaussian--Lagrange identity

In `k[[u]]`,

\[
 \boxed{
 \mathbb E\!\left(A(Z)e^{uW\cdot\Phi(Z)}\right)
 =\frac{A(g(u))}
 {\det\!\left(I-uJ\Phi(g(u))\right)}.}                \tag{1.3}
\]

The denominator has constant term one.  Formula (1.3) is valid whether or not
`Phi` has a constant term.

## 2. The completed rings and coefficientwise finiteness

Use the complete local ring

\[
 S=k[[u,z_1,\ldots,z_r]]
\]

with maximal ideal `(u,z_1,...,z_r)`.  Put

\[
 F_u(z)=z-u\Phi(z).                                    \tag{2.1}
\]

Every `F_{u,i}` belongs to the maximal ideal, even when `Phi_i(0) != 0`.
The continuous endomorphism fixing `u` and sending `z` to `F_u(z)` is an
automorphism: modulo `u` it is the identity, and its inverse

\[
 G(u,y)=y+u\Phi(G(u,y))\in k[[u,y]]^r                 \tag{2.2}
\]

is obtained uniquely, one `u`-coefficient at a time.  Evaluation at `y=0`
is continuous and gives `G(u,0)=g(u)`.

All substitutions in (1.3) are therefore legal.  In particular, every
component of `g` has positive `u`-order, so the coefficient of `u^n` in
`A(g(u))` uses only monomials of `A` of total degree at most `n`.  The
determinant in (1.3) is a unit because it is congruent to one modulo `u`.

The left side is also coefficientwise finite.  The multinomial expansion is

\[
 e^{uW\cdot\Phi(Z)}
 =\sum_{\alpha\in\mathbb N^r}
   \frac{u^{|\alpha|}}{\alpha!}W^\alpha\Phi(Z)^\alpha. \tag{2.3}
\]

For a fixed power `u^n`, only the finitely many multi-indices with
`|alpha|=n` occur.  Wick contraction of such a term extracts one coefficient
of `A Phi^alpha`; hence it is defined even for formal `A`.

For the residue calculation, use the coefficientwise Laurent module

\[
 \mathscr L=
 \left\{\sum_{n\ge0}u^nL_n(z):
 L_n\in k[[z_1,\ldots,z_r]][z_1^{-1},\ldots,z_r^{-1}]
 \right\}.                                            \tag{2.4}
\]

Equivalently, this is the `u`-adic completion of

\[
 S[z_1^{-1},\ldots,z_r^{-1}].
\]

Allowing the lower Laurent bound to depend on the `u`-degree is essential
when the zero section moves.

Each `L_n` has a finite lower exponent in every `z_i`.  Define

\[
 \operatorname{Res}_z L(z)\,dz_1\cdots dz_r
 =[z_1^{-1}\cdots z_r^{-1}]L(z)                       \tag{2.5}
\]

coefficientwise in `u`.  Expressions such as

\[
 (z_i-u\Phi_i(z))^{-b-1}
 =z_i^{-b-1}\sum_{a\ge0}\binom{b+a}{a}
 \left(\frac{u\Phi_i(z)}{z_i}\right)^a               \tag{2.6}
\]

belong to `\mathscr L`: at a fixed `u`-degree only finitely many values of `a`
contribute.  This is the precise completion behind every inverse denominator
below.

## 3. Formal residue change with a moving zero section

We need the change rule for `F_u`, whose value at `z=0` need not vanish as a
series in `u`.  Introduce an auxiliary scalar `s` and set

\[
 F_s(z)=z-su\Phi(z).                                    \tag{3.1}
\]

Over `k[s][[u,z]]`, this is again a continuous automorphism congruent to the
identity modulo `u`.  Pullbacks of Laurent forms are expanded as in (2.6).
Every fixed `u`-coefficient is a finite polynomial in `s`, so differentiation
in `s` and evaluation at `s=0,1` are legal coefficientwise.

Let `omega` be a relative top form in the `z` variables with coefficients in
the module (2.4), and let

\[
 V_s=(\partial_sF_s)\circ F_s^{-1}
\]

be the formal velocity field.  The usual chain rule, valid coefficientwise,
gives

\[
 \frac d{ds}F_s^*\omega=F_s^*(\mathcal L_{V_s}\omega).
\]

Relative top forms are closed.  Cartan's formula therefore gives

\[
 \mathcal L_{V_s}\omega=d_z(\iota_{V_s}\omega).        \tag{3.2}
\]

The formal residue of an exact relative top form is zero.  Indeed, writing an
`(r-1)`-form as a sum with one `dz_i` omitted, its differential is a sum of
terms `partial C/partial z_i`.  To contribute exponent `-1` in `z_i`, the
undifferentiated monomial would have exponent zero, and differentiation then
multiplies it by zero.  This argument applies to every `u`-coefficient.

It follows from (3.2) that `Res_z F_s^* omega` is independent of `s`.  Taking
`s=0` and `s=1` proves the formal change-of-variables rule

\[
 \operatorname{Res}_y C(y)\,dy
 =\operatorname{Res}_z
 C(F_u(z))\det(DF_u(z))\,dz                         \tag{3.3}
\]

for every coefficientwise Laurent top form for which either side is written
using the expansions (2.6).  Thus the moving constant term `-u Phi(0)` causes
no illegal `z`-adic substitution: it is controlled by the `u`-adic
coefficientwise homotopy.

## 4. The coefficient identity

Write

\[
 j(F_u)=\det(DF_u)=\det(I-uJ\Phi).                     \tag{4.1}
\]

For any `B in S` and `beta in N^r`, apply (3.3) to the displayed residue and
use `G(u,F_u(z))=z`:

\[
\begin{aligned}
 [y^\beta]B(G(u,y))
 &=\operatorname{Res}_y
   B(G(u,y))y^{-\beta-\mathbf1}\,dy\\
 &=\operatorname{Res}_z
   B(z)j(F_u)(z)
   \prod_{i=1}^rF_{u,i}(z)^{-\beta_i-1}\,dz.          \tag{4.2}
\end{aligned}
\]

Expanding each inverse by (2.6) gives

\[
 [y^\beta]B(G(u,y))
 =\sum_{\alpha\in\mathbb N^r}u^{|\alpha|}
  \binom{\beta+\alpha}{\alpha}
  [z^{\beta+\alpha}]
  \Phi(z)^\alpha B(z)j(F_u)(z).                       \tag{4.3}
\]

At each `u`-degree this sum is finite: the displayed factor
`u^{|alpha|}` bounds `|alpha|`, and the remaining series have nonnegative
`u`-order.  On the other hand,

\[
 \frac1{\alpha!}[z^\beta]\partial^\alpha C(z)
 =\binom{\beta+\alpha}{\alpha}
  [z^{\beta+\alpha}]C(z).                             \tag{4.4}
\]

Comparing (4.3) and (4.4), coefficient by coefficient in both `u` and `y`,
proves the parameterized Abhyankar--Gurjar identity

\[
 B(G(u,y))=
 \sum_{\alpha\in\mathbb N^r}\frac{u^{|\alpha|}}{\alpha!}
 \partial^\alpha\!\left(
  \Phi(z)^\alpha B(z)j(F_u)(z)
 \right)\bigg|_{z=y}.                                 \tag{4.5}
\]

Now take

\[
 B(z)=\frac{A(z)}{j(F_u)(z)}\in S.                    \tag{4.6}
\]

This is legal because `j(F_u)=1 mod u`.  The Jacobian factors cancel in
(4.5).  Evaluating at `y=0`, and using `G(u,0)=g(u)`, yields

\[
 \sum_{\alpha\in\mathbb N^r}
 \frac{u^{|\alpha|}}{\alpha!}
 \partial^\alpha(A\Phi^\alpha)(0)
 =\frac{A(g(u))}{j(F_u)(g(u))}.                        \tag{4.7}
\]

This proves the formal Lagrange--Good coefficient identity, including
arbitrary constant terms in `Phi`.

## 5. Wick normalization and determinant orientation

For one pair in (1.2), independence of `X,Y` and
`E(e^{aX})=e^{a^2/2}` give, as a formal identity in `s,t`,

\[
\begin{aligned}
 \mathbb E(e^{sZ+tW})
 &=\exp\!\left(\frac{(s+t)^2}{4}\right)
   \exp\!\left(-\frac{(s-t)^2}{4}\right)\\
 &=e^{st}.
\end{aligned}
\]

Independence therefore gives, for multi-indices `alpha,beta`,

\[
 \mathbb E(Z^\beta W^\alpha)
 =\delta_{\alpha\beta}\,\alpha!.                      \tag{5.1}
\]

If `C(Z)=sum_beta c_beta Z^beta`, then

\[
 \mathbb E(W^\alpha C(Z))
 =\alpha!c_\alpha
 =\partial^\alpha C(0),                               \tag{5.2}
\]

where the derivative convention is
`partial^alpha z^beta|_0=alpha! delta_(alpha,beta)`.
Combining (2.3), (5.2), and (4.7) gives (1.3).

Finally, there is no hidden transpose or reciprocal choice in the
determinant.  With rows indexed by components of `Phi` and columns by input
variables,

\[
 DF_u=I-uJ\Phi.
\]

The identity `F_u(G(u,y))=y` gives

\[
 DF_u(G(u,y))DG(u,y)=I.
\]

At `y=0`, therefore,

\[
 \det DG(u,0)
 =\det(I-uJ\Phi(g(u)))^{-1}.                           \tag{5.3}
\]

This is exactly the denominator orientation in (1.3).

## 6. Weighted-bridge corollary and review status

For the two-component polynomial `Phi` constructed in the weighted Gaussian
bridge, the fixed branch is `(g(u),k(u))` and

\[
 \det(I-uJ\Phi(g(u),k(u)))=1.
\]

Theorem 1.1 immediately reduces the bridge's moment calculation to

\[
 \mathbb E(A(Z_1,Z_2)e^{uP_{H,\lambda}})=A(g(u),k(u)).
\]

Thus the determinant-cancelling construction and univariate Lagrange
inversion are the only remaining bridge-specific steps.

This lemma and the weighted bridge are proved and exactly regression-tested in
the repository, but no external specialist review is recorded.  Because the
lemma carries the constant-term topology, formal residue change, determinant
orientation, and Gaussian normalization simultaneously, it is the most useful
single target for such review.

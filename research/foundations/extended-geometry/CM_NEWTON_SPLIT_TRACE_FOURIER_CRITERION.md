# The Fourier scale of the split cubic trace condition

Let \(B/K\) be a polarized abelian threefold with Rosati-compatible
real multiplication by a totally real cubic field \(F\), defined over
\(K\) and integral outside a finite set. At a good absolute degree-one
place \(v\) of characteristic \(p\) split completely in \(F\), its relative Frobenius
polynomial is \(X^2-a_v X+p\) with \(a_v\in\mathcal O_F\). Its three
relative factors are ordinary precisely where the corresponding trace
is a unit. Thus the number of failed factors is
\(\#\{\mathfrak q\mid p:a_v\in\mathfrak q\}\).
[Split trace criterion](CM_NEWTON_TRACE_COFACTOR_SIEVE.md).

The norm test combines all three factors into an integer cofactor of
size \(O(\sqrt p)\). A trace-dual additive test instead detects a single
factor, with frequencies of order \(p^{1/6}\) in normalized
archimedean coordinates. The corresponding archimedean spacing of the
one-failure lattice is exactly of order \(p^{-1/6}\).

## A short additive character for each coefficient prime

**Proposition.** Write \(\mathfrak D^{-1}\) for the codifferent of \(F\).
For each sufficiently large split \(p\) and each \(\mathfrak q\mid p\),
one can choose
\(\lambda_{\mathfrak q}\in\mathfrak q^{-1}\mathfrak D^{-1}\setminus
\mathfrak D^{-1}\) with all three \(|\sigma(\lambda_{\mathfrak q})|
\asymp_F p^{-1/3}\). For every \(a\in\mathcal O_F\),
\[
       a\in\mathfrak q
       \quad\Longleftrightarrow\quad
       \operatorname{Tr}_{F/\mathbb Q}(\lambda_{\mathfrak q}a)\in\mathbb Z.
\]
Every nonzero vector of the dual ideal \(\mathfrak q^{-1}\mathfrak D^{-1}\)
has supremum norm \(\gg_F p^{-1/3}\) in the real embedding coordinates.

In fact, choose a finite set of fractional ideal representatives for
the ideal classes and write \(\mathfrak q=\alpha J\) using one of them.
By Dirichlet's unit theorem, multiply \(\alpha\) by a unit so that
\(|\sigma(\alpha)|\asymp_F p^{1/3}\) at all three embeddings. This
follows by taking the vector of logarithms of absolute values modulo
the full unit lattice in the sum-zero plane after subtracting its mean.
For each representative fix a nonzero \(\xi_J\in J^{-1}\mathfrak D^{-1}\)
and set \(\lambda_{\mathfrak q}=\alpha^{-1}\xi_J\). Its norm is a
nonzero fixed constant divided by \(p\). Since the norms of nonzero
elements of \(\mathfrak D^{-1}\) have a positive lower bound, the chosen
element lies outside that ideal for sufficiently large \(p\).
Every nonzero element of \(J^{-1}\mathfrak D^{-1}\) likewise has norm
bounded below by a positive constant depending on \(J\) and \(F\).
Scaling by \(\alpha^{-1}\) and taking a maximum proves the lower bound
for the dual ideal.

The trace pairing maps the dual ideal to additive characters on
\(\mathcal O_F/\mathfrak q\), with kernel exactly \(\mathfrak D^{-1}\)
by its definition. The character of the chosen element is nontrivial on
an additive group of prime order \(p\), so has trivial kernel there.
For \(x_v=(\sigma(a_v)/\sqrt p)_\sigma\), the phase of the character is
the dot product of \(x_v\) with the frequency vector
\((\sqrt p\,\sigma(\lambda_{\mathfrak q}))_\sigma\). Each coordinate
of the chosen vector has size comparable to \(p^{1/6}\), and no
nonzero representative of such a character has smaller supremum norm
than a fixed multiple of \(p^{1/6}\).

## A sufficient cancellation estimate for a fixed variety

Let \(S(X)\) be the indicated good degree-one split places of norm at
most \(X\), and assume their number tends to infinity. Fix the balanced
choice of each \(\lambda_{\mathfrak q}\) and write
\(e(t)=\exp(2\pi i t)\). Consider the real parts of the sums
\[
 H_m(X)=\sum_{v\in S(X)}\sum_{\mathfrak q\mid p_v}
         e\!\left(m\operatorname{Tr}_{F/\mathbb Q}
                         (\lambda_{\mathfrak q}a_v)\right),
 \qquad m\ge1.
\]

**Corollary.** If \(\operatorname{Re}H_m(X)=o(|S(X)|)\) for every fixed
positive integer \(m\), then the ordinary places have relative natural
density one in the split set. The coefficient primes can be averaged in
these sums and their three individual distributions need not be
estimated.

For any integer \(R\ge1\), the nonnegative periodic function
\(R^{-2}|\sum_{j=0}^{R-1}e(jt)|^2\) is one at every integer. Expanding
it bounds the number of nonordinary places in \(S(X)\) by
\[
 \frac{3|S(X)|}{R}+\frac{2}{R}\sum_{m=1}^{R-1}
        (1-m/R)\operatorname{Re}H_m(X).
\]
Here each nonordinary place contributes at least one failed label.
Taking \(X\to\infty\) with \(R\) fixed and then \(R\to\infty\) proves
the assertion. More generally the limsup of the normalized expression
on the right gives an upper bound for the relative density. This
formulation only asks for cancellation in the real parts of the
selected modes, with the conductor prime in the varying additive
character. The needed cancellation for fixed independent real-cubic
threefolds remains open.

The same test covers coefficient splitting type \((1,2)\) using its
single degree-one coefficient prime. Indeed on making \(\ell\)-torsion
rational for one fixed \(\ell>3\), absolute degree-one reductions have
\(a_v\ne0\) since their rational trace is \(6\bmod\ell\). For all large
unramified \(p\), \(|N_F(a_v)|\le8p^{3/2}<p^2\), so no coefficient prime
of residue degree at least two can divide \(a_v\). The local rank-two
ordinary criterion also holds at these coefficient places. Thus over
this extension one may sum over all good degree-one places, using only
the degree-one coefficient primes in each inner sum. This gives the
same Fejér bound with at most three labels, and no term for inert \(p\).

## Sharp spatial scale and locally realizable data

**Proposition.** Equip the trace cube \([-2,2]^3\) with the supremum
norm. For every large split \(p\), normalized traces of simple
prime-field RM threefolds with exactly one failed coefficient prime
are \(O_F(p^{-1/6})\)-dense in the cube. This holds with the failed prime
prescribed. Conversely, even the union over all three possible failed
primes has covering radius \(\gg_F p^{-1/6}\) in the cube.

Write \(\sigma\) for the Minkowski embedding. The balanced ideal
description above gives the basis matrix of
\(\Lambda_{\mathfrak q}=p^{-1/2}\sigma(\mathfrak q)\) as
\(p^{-1/2}\operatorname{diag}(\sigma_i(\alpha))B_J\), with \(B_J\) one
of a fixed finite set of invertible matrices. Its basis lengths are
\(O_F(p^{-1/6})\), and the inverse matrix has norm \(O_F(p^{1/6})\).
A fundamental parallelepiped supplies the upper covering bound.
Shrink a target toward the interior by \(O_F(p^{-1/6})\) first and,
if necessary, move the selected vector away from zero by one basis
vector. This supplies a nonzero \(a\in\mathfrak q\) with all
\(|\sigma_i(a)|<2\sqrt p\). Since \(|N_F(a)|<8p^{3/2}<p^2\), there
is exactly one failed prime, of valuation one. Also \(a\) generates
\(F\): a rational integer in \(\mathfrak q\) is a multiple of \(p\).
The roots of \(X^2-aX+p\) generate a CM quadratic extension of \(F\).
Its two ordinary split primes have integral slopes and its one
ramified coefficient prime has local degree two and slope \(1/2\).
The Honda–Tate invariants all vanish, giving a simple threefold with
that relative Frobenius polynomial.
[Honda–Tate calculation for the same pattern](CM_NEWTON_SPLIT_POLYNOMIAL_OBSTRUCTION.md#avoiding-fixed-cancellations-on-actual-frobenius-data).

The inverse matrix bound gives \(O_F(p^{1/2})\) points of each lattice
in \([-3,3]^3\), by counting their integer coordinate vectors in a
box of side \(O_F(p^{1/6})\). The union of supremum balls of radius
\(\epsilon p^{-1/6}\) about these points in all three lattices has
volume \(O_F(\epsilon^3)\). For a fixed small \(\epsilon>0\) the union
cannot cover the trace cube. Points of the lattices outside the larger
cube cannot affect this conclusion. This proves the lower bound.

The same upper construction works for any prescribed congruence on
\(a\) modulo a rational integer \(D\) coprime to \(p\), with error
\(O_F(Dp^{-1/6})\), by using the affine coset of \(D\mathfrak q\),
provided \(Dp^{-1/6}\) is sufficiently small in terms of \(F\).
Consequently a bounded family of tests of the trace whose Lipschitz
constants are \(o(p^{1/6}/D)\) cannot distinguish prescribed trace
targets from the approximating one-failure data. At the critical
frequency the additive test is identically one on the corresponding
failed ideal, despite the vanishing spacing of the lattice. This
sharpens the quantitative interpretation of the existing
[joint local realizability result](CM_NEWTON_THREE_BALANCED_GLOBAL.md#simultaneous-archimedean-and-auxiliary-frobenius-data).

## Exact trace phases at 113

In the [level-97 example](CM_NEWTON_SPLIT_FAILURE_WITH_FULL_RESIDUAL_IMAGE.md)
write \(F=\mathbb Q[y]/(f(y))\), \(f(y)=y^3-y^2-2y+1\) and
\(a=8-y-3y^2\) at \(p=113\). If \(\beta\) generates a coefficient
prime, take \(\lambda=(\beta f'(y))^{-1}\). For a polynomial \(C\) of
degree at most two, Lagrange interpolation at the three roots of \(f\)
gives \(\operatorname{Tr}(C(y)/f'(y))=[y^2]C\). Thus multiply \(a\)
by \(p/\beta\) and take the quadratic coefficient after reduction
modulo \(f\), then divide by \(p\). The following table gives all three
exact phases. The displayed \(\beta\)'s have norm \(-113\) and reduce
to zero at the indicated root modulo 113, so each generates the
claimed prime.

| Root modulo 113 | \(\beta\) | \(113/\beta\) | \(113\operatorname{Tr}(\lambda a)\) |
| ---: | --- | --- | ---: |
| 34 | \(-1-y-3y^2\) | \(-57-4y+17y^2\) | 153 |
| 89 | \(-2-y-2y^2\) | \(-42+y+9y^2\) | 104 |
| 104 | \(-1-2y-3y^2\) | \(-49+3y+11y^2\) | 113 |

The unique integral phase corresponds to the failed prime labelled 104.

## Comparison with trace-norm distribution

Suh's proof of Theorem 4.1.1 uses joint Sato–Tate equidistribution
against fixed measure-zero loci when the trace norm divided by the
Weil power is an integer. For one split failure in a cubic field that
quotient is an integer divided by \(\sqrt p\). The Fourier criterion
adapts the equidistribution strategy to a moving frequency in that case.
[Suh, §3.3 and Theorem 4.1.1](https://arxiv.org/pdf/2410.01182).

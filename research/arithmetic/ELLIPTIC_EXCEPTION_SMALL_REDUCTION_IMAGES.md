# Small reduction images and polynomially vanishing absolute losses

The rank-two height pairing bounds how often the rational reduction
image can be small. This gives a polynomially vanishing absolute
local loss at a density-one set of primes, uniformly over the earlier
finite conditions. It does not bound the relative loss after
conditioning on a shrinking surviving set.

## 1. A height bound for all small images

Retain `Gamma=E(Q)` on the existing elliptic slice. Fix the independent
points `P,Q` certified in
[ETC3a](ELLIPTIC_EXCEPTION_COMPLETE_THREE_ADIC_CLASSIFICATION.md).
No full basis or saturation is required here. At a good prime put

\[
 R_p=\operatorname{im}(\Gamma\longrightarrow E(\mathbb F_p)),
 \qquad r_p=|R_p|.
\]

Write `X` for the short-model coordinate `x-3552`, and `H_X` for its
projective rational height. The bounded height comparison and positive
quadratic pairing in
[EEC1](ELLIPTIC_EXCEPTION_COUNT.md) give fixed constants such that

\[
 \log H_X(D)\le q_x(D)+C_0,\qquad
 q_x(uP+vQ)\le C_1(u^2+v^2).                       \tag{1.1}
\]

Constants denoted by `E` below include these fixed, certified points.

**Proposition ESR1 (weighted bound for small reduction images).**
For every real `z>=2`,

\[
 \boxed{\sum_{\substack{p\text{ good}\\r_p<z}}\log p
                                  \ll_E z^2.}      \tag{1.2}
\]

In particular the sum is over a finite set, and

\[
 \#\{p\text{ good}:r_p<z\}\ll_E\frac{z^2}{\log z}. \tag{1.3}
\]

**Proof.** Put `L=floor(sqrt(z))` and consider the box of
`(L+1)^2>z` rational points

\[
 mP+nQ,\qquad 0\le m,n\le L.
\]

If `r_p<z`, two of these points have the same reduction. Their
difference is one of

\[
 D_{u,v}=uP+vQ,\qquad |u|,|v|\le L,
                                  \quad(u,v)\ne(0,0).          \tag{1.4}
\]

Independence modulo torsion makes every difference a nonidentity
point. There are `O(z)` differences, and each has canonical height
`O_E(z)` by (1.1).

Let `b_(u,v)` be the positive reduced denominator of `X(D_(u,v))`
and form the integer

\[
 B_z=\prod_{(u,v)\text{ in }(1.4)}b_{u,v}.
\]

A nonidentity point reducing to the identity at a good prime has a
pole in `X`; hence the collision above makes `p` divide `B_z`.
Every prime in (1.2) is therefore a prime divisor of this one finite
integer. Moreover

\[
 \log B_z\le\sum_{(u,v)}\log H_X(D_{u,v})\ll_E z^2.
\]

Counting each prime only once gives
`sum log p <= log rad(B_z) <= log B_z`, proving (1.2).
Repeating a difference or its negative in the product causes no
problem: it only increases this upper bound.

For (1.3), the primes at most `z` contribute at most `z` to the
cardinality. Every remaining prime contributes at least `log z`
to (1.2). Thus the cardinality is at most
`z+O_E(z^2/log z)=O_E(z^2/log z)`. \(\square\)

The pole-denominator assertion is the integral-model filtration used
in EOP3 and retained from Milne, *Elliptic Curves* (2006), Chapter II,
Theorem 4.1 and Lemma 5.9, printed pp.62--63 and 67:
[primary book](https://www.jmilne.org/math/Books/ectext6.pdf),
[retained copy](sources/milne-elliptic-curves-2006.pdf).
The integer coordinate translation preserves denominators and reduction
at every good prime. No distribution theorem for point orders is
imported into this proof.

## 2. Density-one order growth and uniform absolute loss

**Corollary ESR2 (polynomial loss bound at density-one primes).**
For every fixed real `A>0` and sufficiently large `X`,

\[
 \boxed{\#\left\{p\le X:p\text{ good},\quad
                 r_p<\frac{\sqrt p}{(\log p)^A}\right\}
                     \ll_{E,A}\frac{X}{(\log X)^{2A+1}}.}
                                                               \tag{2.1}
\]

Thus `r_p>=sqrt(p)/(log p)^A` holds on a natural-density-one set
of rational primes. At every sufficiently large prime of that set,
for **every** finite prime set `S` and every global class `a`,

\[
 \boxed{0\le\mu_{S,a}-\mu_{S\cup\{p\},a}
                              \le\frac{3(\log p)^A}{\sqrt p}.} \tag{2.2}
\]

For comparison, for every fixed `0<epsilon<1/2` the same argument
also gives the separate exceptional-set estimate

\[
 \boxed{\#\{p\le X:p\text{ good},\ r_p<p^{1/2-\epsilon}\}
                   \ll_{E,\epsilon}\frac{X^{1-2\epsilon}}{\log X}.}
                                                               \tag{2.3}
\]

Outside that exceptional set the absolute loss is at most
`3p^(-1/2+epsilon)`. Each prime set and bound is independent of
`S,a`. In particular the single-prime coefficients satisfy, along
the density-one set in (2.1),

\[
 \mu_{\{p\},a}=\frac14+
                         O\left(\frac{(\log p)^A}{\sqrt p}\right)
                                                               \tag{2.4}
\]

uniformly in the four classes.

**Proof.** The function `f(t)=sqrt(t)/(log t)^A` is increasing once
`log t>2A`. Apart from finitely many initial primes, every exceptional
prime in (2.1) satisfies `r_p<f(X)`. Apply (1.3) with `z=f(X)`.
Then `z^2=X/(log X)^(2A)` and `log z` is asymptotic to
`(1/2)log X`, proving (2.1). For (2.3), apply (1.3) directly
with `z=X^(1/2-epsilon)`. Both exceptional bounds are
`o(pi(X))`, using the elementary lower bound
`pi(X) >> X/log X`. Finitely many bad or small primes do not
affect density.

For completeness, that elementary prime-count lower bound follows
from `4^n/(2n+1)<=binom(2n,n)`. Legendre's factorial formula gives
`v_p(binom(2n,n))<=floor(log(2n)/log p)`, so
`log binom(2n,n)<=pi(2n)log(2n)`. Comparing the two bounds and
using monotonicity of `pi` proves the estimate for all large `X`.

The intrinsic unconditioned local failure mass is at most `3/r_p`
by [EOP3](ELLIPTIC_EXCEPTION_ONE_PRIME_POSITIVITY.md). At an active
good prime, failure can occur only in the pole residue or the at
most two finite zero residues, and intrinsic Haar measure pushes
forward uniformly to `R_p`. The mass removed from a class-specific
finite-place set is an intersection with that same failure event,
so it is bounded by `3/r_p` regardless of correlations with `S`.
For primes one modulo three it is zero. This proves (2.2).
The identical argument gives its power-epsilon version. Taking `S`
empty, with class mass `1/4`, proves (2.4).
\(\square\)

These are estimates for Haar coefficients, not counting errors
uniform in a growing set of places. The conditional relative loss
would divide (2.2) by `mu_(S,a)`, which can shrink as `S` varies.
Neither the common quotients nor the conditional covariance sums of
[ECM1--ECM3](ELLIPTIC_EXCEPTION_COMMON_QUOTIENT_MOMENTS.md) are bounded
by this argument. It supplies no return recurrence, all-primes
measure value or Hasse-count asymptotic.

The proof is the finite collision box, product of denominators and
height comparison above. It requires no new finite experiment or
checker. The rank, independence, height and local-loss inputs retain
their separate proof sources. No formal verification or external
review is claimed.

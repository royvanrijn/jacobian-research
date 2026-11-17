# A compulsory quadratic correlation in pole returns

The unconditional pole-return sums cannot have variance negligible
compared with their squared means. Rational two-torsion forces a fixed
quadratic correlation at every good prime: the variance is always at
least one seventh of the squared mean. This obstructs that particular
second-moment route, irrespective of reduction-order growth. It does
not obstruct recurrence itself, or automatically apply after fixing
a class or a finite prefix.

## 1. The unavoidable unconditional variance

Retain `Gamma=E(Q)` and `Pi=widehat(Gamma)` with normalized Haar
probability `P`. The rank and exact rational torsion theorem in
[the elliptic arithmetic note](ELLIPTIC_EXCEPTION_ARITHMETIC.md) give

\[
 \Gamma/2\Gamma\simeq(\mathbb F_2)^3.
\]

Let `X_2` be its seven nontrivial quadratic characters, viewed as
functions on `Pi`. At each good odd prime use the actual rational
reduction map and pole event

\[
 \phi_p:\Pi\longrightarrow R_p,\qquad
 r_p=|R_p|,\qquad A_p=\{\phi_p=O\}.
\]

Every `R_p` has even order: the rational point `T_0` of order two
has nonidentity reduction at a good odd prime. Write `X_p` for
the characters of `Pi` factoring through `R_p`.

**Theorem EVO1 (quadratic variance obstruction).** For any finite
set `I` of good odd primes and nonnegative real weights `w_p`, put

\[
 S=\sum_{p\in I}w_p\mathbf1_{A_p},\qquad
 M=\mathbb E S=\sum_{p\in I}\frac{w_p}{r_p}.
\]

Then

\[
 \boxed{\operatorname{Var}(S)\ge\frac{M^2}{7}.}      \tag{1.1}
\]

In particular no unconditioned nonnegative pole-weight sum with
`M` tending to infinity can satisfy the variance hypothesis of ECM3.
For weights `0<=w_p<=1`, its exact common-quotient correlations obey

\[
 \boxed{2\sum_{p<q\in I}
       \frac{w_pw_q(d_{p,q}-1)}{r_pr_q}
                                  \ge\frac{M^2}{7}-M.}        \tag{1.2}
\]

Here `d_(p,q)` is the actual common-quotient order of ECM1.
Furthermore a set of good odd primes with `d_(p,q)=1` for every
distinct pair has at most seven members.

**Proof.** Character orthogonality gives, for `chi in X_2`,

\[
 \mathbb E(\mathbf1_{A_p}\chi)=
 \begin{cases}1/r_p,&\chi\in X_p,\\0,&\chi\notin X_p.
 \end{cases}
\]

Define the nonnegative numbers

\[
 a_\chi=\mathbb E(S\chi)
             =\sum_{\substack{p\in I\\\chi\in X_p}}w_p/r_p.
\]

Every even finite abelian group `R_p` has a nontrivial quadratic
character. Its pullback is one of `X_2`, since `phi_p` is onto.
Thus each summand of `M` contributes to at least one `a_chi`, and

\[
 \sum_{\chi\in X_2}a_\chi\ge M.
\]

The seven characters have mean zero and are orthonormal. Expanding
the nonnegative mean square of
`S-M-sum_(chi in X_2) a_chi chi` proves

\[
 \operatorname{Var}(S)\ge\sum_{\chi\in X_2}a_\chi^2
           \ge\frac17\left(\sum_{\chi\in X_2}a_\chi\right)^2,
\]

where the second inequality is elementary Cauchy--Schwarz.
This proves (1.1) using only a finite quadratic quotient.

By [ECM1](ELLIPTIC_EXCEPTION_COMMON_QUOTIENT_MOMENTS.md),

\[
 \operatorname{Var}(S)=
 \sum_p w_p^2(1/r_p-1/r_p^2)
       +2\sum_{p<q}\frac{w_pw_q(d_{p,q}-1)}{r_pr_q}.
\]

For `0<=w_p<=1` the first sum is at most `M`, giving (1.2).
If `d_(p,q)=1`, the two sets `X_p cap X_2` and `X_q cap X_2`
are disjoint, since a shared nontrivial character would give a
nontrivial common quotient. They are all nonempty subsets of the
seven-element set `X_2`. This proves the final assertion.
\(\square\)

For comparison the same proof at any fixed prime `ell` gives

\[
 \operatorname{Var}(S)\ge
 \frac{(\ell-1)^2}{|\Gamma/\ell\Gamma|-1}
   \left(\sum_{\substack{p\in I\\\ell\mid r_p}}w_p/r_p\right)^2.
                                                               \tag{1.3}
\]

Each image of order divisible by `ell` supplies at least `ell-1`
nontrivial characters of order `ell`; their pole coefficients are
again positive. At two all the primes contribute, giving the
compulsory bound (1.1). At an odd `ell`, (1.3) only obstructs a
small variance if that indicated portion of the mean remains large.

## 2. What changes in a fixed congruence coset

Conditioning cannot be omitted from this conclusion. Pole Fourier
coefficients lose their common nonnegative form after translating
to a nonidentity coset. The following exact refinement records the
affine compatibility instead.

Let `H subset Pi` be an open subgroup and let `C=gamma+H` be a
fixed coset, with normalized conditional Haar probability. Define

\[
 R_p^H=\phi_p(H),\quad t_p=|R_p^H|,\quad
 K_p^H=H\cap\ker\phi_p.
\]

If `-phi_p(gamma)` is not in `R_p^H`, the pole event is impossible
in `C`. Otherwise choose `u_p in H` with
`phi_p(u_p)=-phi_p(gamma)`.

**Proposition EVO2 (affine relative common quotients).** At every
individually possible pole prime,

\[
 \mathbb P(A_p\mid C)=1/t_p.
\]

For two such primes put

\[
 D^H_{p,q}=H/(K_p^H+K_q^H),\quad d^H_{p,q}=|D^H_{p,q}|.
\]

Then

\[
 \boxed{\mathbb P(A_p\cap A_q\mid C)=
 \begin{cases}
 d^H_{p,q}/(t_pt_q),&u_p-u_q\in K_p^H+K_q^H,\\
 0,&\text{otherwise}.
 \end{cases}}                                                  \tag{2.1}
\]

The compatibility condition does not depend on the chosen lifts.
In particular `d^H_(p,q)<=gcd(t_p,t_q)`. For `H=2Pi`, one also
has, with `g=gcd(exp(R_p^H),exp(R_q^H))`,

\[
 d^H_{p,q}\le g^2.                                  \tag{2.2}
\]

**Proof.** Identify `C` with `H` by translation. A pole event is
the coset `u_p+K_p^H`. The two cosets intersect precisely under
the compatibility condition in (2.1). When nonempty their
intersection is one coset of `K_p^H cap K_q^H`. The index identity
of ECM1 inside `H` gives the displayed measure.

The common quotient is a quotient of both `R_p^H,R_q^H`, so its
order divides both orders and its exponent divides `g`. Finally
`2Pi` kills the rational two-torsion and is isomorphic to
`widehat(Z)^2`. A quotient of this group killed by `g` has order
at most `g^2`, proving (2.2). \(\square\)

The relative quotients and compatibility in (2.1) are finite group
data. They can be computed when the corresponding joint quotients
and the coset are certified, using the supplied-input scope of
[FQS2--FQS3](ELLIPTIC_EXCEPTION_FINITE_QUOTIENT_SATURATION.md).
Fixing a full coset modulo `2Pi` removes the specific compulsory
quadratic variation of EVO1. It does not establish any small bound
for the remaining relative quotients.

## 3. A sufficient relative-correlation condition

The relative quotient gives a useful sufficient condition that does
not require asserting independence. Choose only primes where the
pole event is individually possible in `C`, let `0<=w_p<=1`, and put

\[
 M_Y^C=\sum_{p\le Y}w_p/t_p.
\]

**Corollary EVO3 (a conditional pole-recurrence criterion).** If

\[
 M_Y^C\longrightarrow\infty,\qquad
 \sum_{p<q\le Y}\frac{w_pw_q(d^H_{p,q}-1)}{t_pt_q}
                                      =o((M_Y^C)^2),           \tag{3.1}
\]

then `sum_p w_p 1_(A_p)=infinity` almost surely in `C`.
For the pole weights `w_p=b_p` on the active primes of
[ETP4--ETP5](ELLIPTIC_EXCEPTION_TRIANGULAR_PRIME_SIEVE.md), with
zero weights elsewhere, this makes the all-places soluble subset
of `C` Haar-null.

**Proof.** Formula (2.1) bounds the conditional covariance above
by `(d^H_(p,q)-1)/(t_p t_q)`; an incompatible pair has negative
covariance and only improves this bound. The diagonal variance is
at most `M_Y^C`. Thus the conditional variance divided by
`(M_Y^C)^2` tends to zero under (3.1). Apply the elementary
monotone second-moment argument of ECM3. The pole weights are
bounded above by the exact hazards, so ETP5 gives the last claim.
\(\square\)

The weaker signed condition obtained by replacing `d^H_(p,q)-1`
in (3.1) by `epsilon_(p,q)d^H_(p,q)-1` also suffices, where
`epsilon_(p,q)` is the compatibility indicator in (2.1).
That expression is the exact off-diagonal covariance sum. The
unsigned condition (3.1) avoids relying on cancellations, but can
be strictly stronger; affine incompatibility can contribute negative
covariances.

For example, a positive-density set of individually possible primes
with `w_p` bounded below would make the mean diverge, since
`t_p<=#E(F_p)<2p` at large good primes. Trivial pairwise relative
quotients would imply the second condition automatically. Neither
such a prime set nor the weighted relative-correlation estimate is
proved here. The gcd bounds in EVO2 are elementary upper bounds,
not a distribution theorem for their average.

Almost-sure recurrence in a coset also holds after restriction to
any positive-measure finite-place prefix inside that coset. A fixed
global class, however, is generally a union of full mod-two cosets;
it is not permissible to apply EVO1 directly to its conditioned
probability. An all-class or fixed-class conclusion requires the
appropriate coverage by cosets where a recurrence criterion has
actually been proved.

## 4. Scope of the obstruction

The compulsory lower bound concerns **unconditioned pole returns
with nonnegative weights**. Full branch-return hazards contain
translated zero branches, whose Fourier coefficients can have
different signs. The argument does not automatically give a variance
floor for those hazards, for a class-conditioned prefix, or for a
translated congruence coset. It proves failure of one strong
second-moment hypothesis, not failure of recurrence or positivity
of the all-primes measure.

The rank, torsion and good-reduction inputs are already retained
in the elliptic arithmetic and local-image notes. The proof here
is finite-character orthogonality, a mean-square inequality and
exact affine subgroup intersections; it needs no new finite checker.
No formal verification or external review is claimed.

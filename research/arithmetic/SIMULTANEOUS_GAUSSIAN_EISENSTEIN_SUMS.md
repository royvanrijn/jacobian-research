# Two summands with both prime restrictions

The general framework is developed in
[Binary sums with prescribed prime factors](BINARY_SUMS_WITH_PRESCRIBED_PRIME_FACTORS.md).
This case study retains the norm-specific proofs and certificates for modulus 12.
The [quantitative subgroup correlation theorem](SPLITTING_SEMIGROUP_LOCAL_GLOBAL.md)
now proves both the unbalanced every-target statement S12U and balanced S12
below, with the full \(C_{12}+C_{12}\) leading asymptotic and logarithmic
remainder. Its fixed-modulus balanced corollary covers every fixed \(q\ge3\).
The [uniform power-short binary lower bound FBS-U](BINARY_SUMS_WITH_PRESCRIBED_PRIME_FACTORS.md#the-central-binary-conjecture)
remains a stronger open target.

The separate results retained here
include semigroup-scale density-one coverage for fixed power exponents above
\(1/3\), full unweighted variance in fixed-width balanced windows and in
power windows with \(\theta>5/6\),
positive-proportion coverage for every fixed positive power exponent,
a uniform local sieve margin,
a large-prime tail bound, averaged coordinate counts before and after the full
sieve, and an exact finite search. The fixed-fraction positive-proportion
statement is a corollary of Loughran--Matthiesen's correlation theorem; the averaged statements
are consequences of established ideal-counting estimates.
The local calculation includes singular lifts and primes dividing the target.
These earlier statements alone do not supply an every-target theorem; the linked subgroup correlation supplies it by a separate argument.

A separate [one-sided pointwise theorem](ASYMMETRIC_GAUSSIAN_EISENSTEIN_SUMS.md)
now proves a power-error asymptotic for
\(\sum_{a+b=N}\tau(a)^2\mathbf1_{C_{12}}(a)\tau(b)\mathbf1_{A_3}(b)\)
on every sufficiently large \(N\equiv2\pmod6\). Its explicit positive
singular series gives eventual \(C_{12}+A_3\) coverage and many distinct
pairs. The [uniform coefficient and cofactor theorem](ASYMMETRIC_SEMIGROUP_SCALE.md)
now gives matching unweighted counts of order
\(N J_E(N)/(\log N)^{5/4}\) on every sufficiently large such target.
The second weight has a two-character expansion, so its proof does
not supply the four-character versus four-character correlation needed here.
The [general fixed-modulus theory](ASYMMETRIC_FIXED_MODULUS_THEORY.md)
also gives every sufficiently large \(N\equiv2\pmod{12}\)
\(\gg N J_{12}(N)/(\log N)^{3/2}\) pairs with \(a\in C_{12}\),
\(b\in A_3\), and every prime \(7\bmod12\) dividing \(b\) exceeding
\(N^\gamma\), for some sufficiently small fixed \(\gamma>0\).
Here \(J_{12}(N)\) is the product of \((p-1)/(p-2)\) over target
primes \(p\ge5\) outside \(1\bmod12\). There are fewer than
\(1/\gamma\) such excluded factors, with multiplicity.
This does not remove them or provide the raw simultaneous-norm lower bound.
The two-sided targets S12U and S12 now follow from the linked subgroup correlation. Their historical priority audit and the raw coordinate-weight lower bound remain separate obligations.

The general theorem FBS1 gives, for every fixed \(0<\theta<1\), at least
\(c_2X\) targets \(X<N\leq2X\), \(N\equiv2\pmod {12}\), with
\(\gg_\theta X^\theta/(\log X)^{3/2}\) pairs from \(C_{12}\) satisfying
\(|a-b|\leq X^\theta\). Its proof uses a block first moment and Henriot's
corrected upper bound, without a short-interval asymptotic. These are unweighted
pairs; this result does not extend the exact weighted average A12G to shrinking
windows. The stronger [semigroup theorem Dq](PRESCRIBED_PRIME_SEMIGROUP_DENSITY_ONE.md)
gives, for every fixed \(1/3<\theta<1\), the same semigroup-scale lower
bound outside \(O_{\theta,A}(X/(\log X)^A)\) compatible targets for each
fixed \(A>0\). This is derived using prime times small allowed cofactors
and the 2026 short-interval Fourier theorem. The [variance theorem Vq](PRESCRIBED_PRIME_BINARY_VARIANCE.md)
also gives a target-dependent unweighted asymptotic for almost all compatible
targets when \(H=\lambda X\), with fixed \(0<\lambda\le1\). The
[power-window theorem VqPower](PRESCRIBED_PRIME_POWER_WINDOW_VARIANCE.md)
now gives the full main term \(c_{12}^2H(\log X)^{-3/2}S_{12}(N)\) for
\(H=X^\theta\), \(5/6<\theta<1\), with a normalized mean-square
saving \((\log X)^{-1/4}\). The initial sharp four-form energy target
\(\theta=99/100\) is therefore resolved. These are different statements
from the weighted coordinate average A12G. The separate
[D12Power theorem](PRESCRIBED_PRIME_POWER_SAVING_EXCEPTIONS.md) gives
\(r_{12,H}(N)\gg H/(\log X)^{3/2}\), with \(H=X^{99/100}\), outside
\(O(X^{1-10^{-6}})\) compatible targets. This lower-count result comes
from a bounded weighted subset, rather than a full-count asymptotic.
The [raw and fully sifted correlation theorem SimNormPower](SIMULTANEOUS_NORM_POWER_CORRELATIONS.md)
now gives exact coordinate-weight main terms with power-saving mean square.
Outside one \(O(X^{1-7\cdot10^{-6}})\) exceptional set,
\(R(N)\asymp N\mathcal T(N)\) and
\(B_p(N)/R(N)=\beta_p(N)+O(X^{-10^{-6}})\) simultaneously for forbidden
\(p\le X^{10^{-8}}\), including arbitrary target valuations. The fully
sifted weighted proportion is
\(\Gamma(N)=\prod_p(1-\beta_p(N))\ge10-\pi^2\), up to the same error.
Fixed moving balance intervals have the analogous conclusions with a smaller
power error. These weighted statements do not on their own prove every-target coverage.
The older H12 prime-pair baseline has a smaller count scale and the range
\(\theta>2/3\).
For \(2/3<\theta<1\), the linked
[ternary corollary Tq](BINARY_SUMS_WITH_PRESCRIBED_PRIME_FACTORS.md#an-every-target-ternary-baseline)
also represents every sufficiently large \(N\equiv3\pmod {12}\) by three
primes \(1\pmod {12}\), each within \(N^\theta\) of \(N/3\). These are
established baselines; the remaining programme seeks binary estimates using
the whole semigroup: full-count variance in smaller shrinking windows,
a power-saving full-count variance, broader lower-count exception estimates
and a uniform every-target lower bound.

## The target

Put

\[
C=C_{12}=\{c\geq1:p\mid c\Longrightarrow p\equiv1\pmod {12}\},\qquad 1\in C.
\]

**Target S12.** For every sufficiently large integer \(N\equiv2\pmod {12}\),
there are \(a,b\in C\) with

\[
N=a+b,\qquad N/3\leq a,b\leq2N/3.
\]

This is now proved by Corollary FixedModulusBalancedExistence at \(q=12\)
in the [subgroup theory](SPLITTING_SEMIGROUP_LOCAL_GLOBAL.md#balanced-fixed-modulus-families).
The balanced count has leading constant \(A_{12}/3\), target factor
\(J_{12}(N)\), and relative error \(O((\log N)^{-1/32768})\).
The analytic inputs and written proof have internal source/application
review; no external human review or formal verification is asserted.

The conjectured unbalanced version for every positive admissible target agrees
with the finite search below. An effective eventual theorem would need a cutoff
overlapping the verified range, or the intervening range would also have to be
checked. The checked range alone supplies no eventual cutoff.

The [Eisenstein manuscript](../../unpublished/eisenstein-restricted-prime-sums/README.md)
proves binary coverage with only the restriction \(p\equiv1\pmod3\).
Requiring both quadratic splitting conditions selects the single class
\(1\pmod {12}\). The target is not a consequence of that manuscript.

## The first pointwise target and the scaling reduction

**Unbalanced target S12U.** Every sufficiently large \(N\equiv2\pmod{12}\)
has \(N=a+b\) with positive \(a,b\in C_{12}\), without a balance restriction.
This follows from the proved balanced S12. The unrestricted positive
unweighted correlation and its fixed local-state versions are also given
directly by the [quantitative subgroup theorem](SPLITTING_SEMIGROUP_LOCAL_GLOBAL.md).
No converse from the bare S12U existence assertion is used.

**Scaling lemma S12Scale.** Factor uniquely
\[
N=gN_0,\qquad
g=\prod_{p\equiv1\ (12)}p^{v_p(N)}\in C_{12}.
\]
Then \(N_0\equiv2\pmod{12}\) has no prime divisor \(1\pmod{12}\).
A representation \(N_0=a+b\), \(a,b\in C_{12}\), gives
\(N=ga+gb\). It preserves each relative share, since
\(ga/N=a/N_0\), and is injective on value pairs. These facts follow directly
from multiplicative closure of \(C_{12}\).

There is an important quantifier limitation. A large target can have a bounded
reduced core: \(N=14g\) with arbitrarily large \(g\in C_{12}\) always has
core 14. A theorem for sufficiently large reduced cores does not on its own
cover these targets. For unbalanced existence, the finite certificate covers
all admissible cores through \(10^8\). To complete this reduction, the eventual
core theorem must overlap that range or the remaining smaller cores must be
covered by another argument.

For balanced existence, the known exceptional cores are an additional issue:
\[
14,38,50,62,86,230,242
\]
all have no prime factor \(1\pmod{12}\), and none has a thirds-balanced
representation. Scaling their existing witnesses cannot prove balance for
the infinitely many targets with one of those cores. Those multiplier
families require separate treatment even if every larger reduced core is
settled. The new subgroup theorem treats the full target directly and therefore covers these multiplier families without this scaling reduction.

Scaling also changes a power-window estimate. If
\(|a-b|\le N_0^\theta\), its scaled difference is at most
\[
gN_0^\theta=g^{1-\theta}N^\theta.
\]
Thus the reduction alone does not preserve a target-sized power window or
its conjectured semigroup-scale count.

## Unconditional coverage for a positive proportion of targets

For a fixed symmetric window \(J=[\lambda,1-\lambda]\),
\(0<\lambda<1/2\), put
\[
r_J(N)=\#\{(a,b)\in C^2:a+b=N,\ a/N\in J\}.
\]
Symmetry of \(J\) means both summands satisfy the same balance condition.
These are ordered pairs of integer values, without coordinate multiplicities.

**Positive-proportion corollary P12.** There are constants \(c_1,c_2>0\),
depending on the fixed window, such that for all sufficiently large real \(X\),
\[
\#\left\{X<N\leq2X:\ N\equiv2\pmod {12},\quad
r_J(N)\geq c_1\frac{N}{(\log N)^{3/2}}\right\}\geq c_2X. \tag{11}
\]
This includes \(J=[1/3,2/3]\) and the narrower fixed window
\(J=[49/100,51/100]\).

The input is [Loughran--Matthiesen, *Frobenian multiplicative functions and
rational points in fibrations*, arXiv:1904.12845v2, Theorem 1.9](https://arxiv.org/html/1904.12845v2#S1.SS4)
(published in *JEMS* 26 (2024), 4779--4830, Theorem 1.9 on p. 4783).
Here is the complete application. Let \(\xi(n)=\mathbf1_C(n)\), extended by
zero to nonpositive integers. On positive integers it is completely
multiplicative and bounded by one. Outside the exceptional primes 2 and 3,
its prime values detect the
identity Frobenius class in
\(\operatorname{Gal}(\mathbb Q(\zeta_{12})/\mathbb Q)\), a group of order four.
Thus it is Frobenian of prime mean \(m(\xi)=1/4\). The growth and prime-power
bounds in Definitions 2.6--2.7 hold with \(H=1\).

Write
\[
M_1(X)=\sum_{X<N\leq2X}r_J(N),\qquad
M_2(X)=\sum_{X<N\leq2X}r_J(N)^2.
\]
For the first moment apply Theorem 1.9 to the forms \(a,b\) in two variables.
For the second moment, expanding \(a+b=c+d\) and eliminating \(d\) gives the
four forms
\[
a,\quad b,\quad c,\quad a+b-c. \tag{12}
\]
Their nonconstant parts are pairwise linearly independent. To check the domains,
use dilation parameter \(B=2X\). The first-moment convex region in
\([0,1]^2\) is
\[
\tfrac12<u+v\leq1,\qquad
\lambda(u+v)\leq u\leq(1-\lambda)(u+v).
\]
The second-moment region in \([0,1]^3\) adds
\(\lambda(u+v)\leq w\leq(1-\lambda)(u+v)\). The complementary value
\(u+v-w\) then satisfies the same bounds. Both regions have positive volume;
\((3/8,3/8)\) and \((3/8,3/8,3/8)\) are interior points where every form
is positive. The theorem's separate arithmetic positivity condition is also
satisfied: the integer witnesses \((1,1)\) and \((1,1,1)\) make every form
equal to 1. Thus positivity is checked globally, not inferred from local factors.

Theorem 1.9 gives positive constants in the two asymptotics, and in particular
constants \(A,B_0>0\) such that, for all sufficiently large \(X\),
\[
M_1(X)\geq A\frac{X^2}{(\log X)^{3/2}},\qquad
M_2(X)\leq B_0\frac{X^3}{(\log X)^3}. \tag{13}
\]
The logarithmic exponents are respectively
\(2(1/4-1)=-3/2\) and \(4(1/4-1)=-3\).

For the threshold argument put
\(T=(A/4)X/(\log X)^{3/2}\), and let
\(E=\{N:X<N\leq2X,\ r_J(N)\geq T\}\).
There are at most \(2X\) integer targets in this interval once \(X\geq1\).
Consequently the targets outside \(E\) contribute at most half the lower
bound for \(M_1\). Cauchy--Schwarz gives
\[
\#E\geq\frac{\bigl(\sum_{N\in E}r_J(N)\bigr)^2}{M_2(X)}
\geq\frac{A^2}{4B_0}X.
\]
Also \(N/(\log N)^{3/2}\leq2X/(\log X)^{3/2}\) on the interval for large
\(X\). Taking \(c_1=A/8\), \(c_2=A^2/(4B_0)\) proves (11).
Every contributing target is automatically \(2\pmod {12}\).

This corollary proves successful balanced decompositions for a positive
proportion, with many pairs on those targets. It does not establish density
one, every-target coverage, or a practical effective cutoff. Averaging the
target is what makes the forms in (12) independent; fixing \(N\) does not give
the same theorem. More generally, a window should have an overlap of positive
width with its reflection if both summands must lie in it; arbitrary narrower
subintervals can make the balance condition empty.

## Encoding both conditions first

Write \(G(x,y)=x^2+y^2\) and \(E(u,v)=u^2-uv+v^2\).
An integer belongs to \(C\) exactly when it is a positive primitive norm value
in both rings. Indeed a Gaussian primitive norm can contain no prime
\(3\pmod4\), and can contain 2 at most once. An Eisenstein primitive norm is
odd, contains no prime \(2\pmod3\), and contains 3 at most once. Intersecting
the two classifications excludes 2 and 3 and leaves precisely the split primes
\(1\pmod {12}\). Conversely, in each unique factorization ring choose only one
prime from each conjugate pair; the resulting coordinates are coprime.

Use four binary coordinate blocks \(X,U,Z,V\), and count solutions to

\[
G(X)=E(U)=a,\quad G(Z)=E(V)=b,\quad a+b=N. \tag{1}
\]

This is three quadratic equations in eight integer coordinates. Put
\(I=[1/3,2/3]\), and let \(R_I(N)\) be the number of all such ordered coordinate
tuples with \(a/N\in I\), equivalently \(N/3\leq a,b\leq2N/3\). At this stage
the coordinates need not be primitive. Write \(R_{\mathrm{tot}}(N)\) for the
unrestricted count with only \(a,b>0\).

For admissible \(N\), every solution automatically has \((ab,6)=1\).
The common values of \(G\) and \(E\) modulo 4 are 0 and 1; their sum is 2,
so both are 1. The Eisenstein norm has values 0 and 1 modulo 3; again both
summands must be 1. Thus the prime obstructions start at 5.

Let \(\mathcal B=\{p\geq5:p\text{ prime},\ p\not\equiv1\pmod {12}\}\).
For \(p\in\mathcal B\), let \(B_{p,I}(N)\) count the tuples in \(R_I(N)\) with
\(p\mid ab\). A tuple surviving every such obstruction has values in \(C\).
Those values have primitive representations in both rings, even if the surviving
tuple's own coordinates have common split-prime factors.

## The local calculation really is summable

For fixed \(p\geq5\) and nonzero \(N\), define \(\beta_p(N)\) as the limiting
proportion of solutions of (1) modulo \(p^k\) with \(p\mid ab\).
Both numerator and denominator use all residue solutions, including singular
ones. The [local proof](SIMULTANEOUS_NORM_LOCAL_DENSITIES.md) establishes the
existence of the normalized limits and the following bounds.

**Local theorem L12.** For every \(p\in\mathcal B\) and positive integer \(N\),

\[
0\leq \beta_p(N)\leq\frac9{p^2}.
\tag{2}
\]

If \(p\nmid N\), the exact formulas are

\[
\beta_p(N)=
\begin{cases}
\displaystyle\frac{2(3p^2-1)}{p(p^3-p^2+3p+1)},&p\equiv5,7\pmod {12},\\[6pt]
\displaystyle\frac2{p(p-1)},&p\equiv11\pmod {12}.
\end{cases} \tag{3}
\]

These are **not** the raw finite-field fractions. In particular the finite-field
loss is of order \(p^{-3}\) at primes \(11\pmod {12}\), while the lifted loss
is of order \(p^{-2}\). Singular lifts matter but do not destroy summability.
The proof of (2) at primes dividing \(N\) uses Cauchy--Schwarz for the local
norm densities and is uniform in \(v_p(N)\).

Even the coarse bound (2) gives a strict margin, uniformly in the target:

\[
\sum_{p\in\mathcal B}\beta_p(N)
\leq9\sum_{\substack{n\geq5\\(n,6)=1}}\frac1{n^2}
=\pi^2-9<1. \tag{4}
\]

Thus this route has no local mass deficit. Turning these local proportions into
proportions among the integer tuples counted by \(R_I(N)\) remains a global task.

## The large-prime tail can already be controlled

**Tail lemma T12.** For every \(\varepsilon>0\), uniformly in \(N\) and
\(P\geq2\),

\[
\sum_{\substack{p\in\mathcal B\\p>P}}B_{p,I}(N)
\ll_\varepsilon\frac{N^{1+\varepsilon}}P. \tag{5}
\]

To prove this, let \(r_G(n),r_E(n)\) count all ordered integer representations
by the two norms. Their divisor formulas give
\(r_G(n)\leq4\tau(n)\), \(r_E(n)\leq6\tau(n)\). Therefore any fixed pair of
values \(a,b\) contributes at most
\(576\tau(a)^2\tau(b)^2\ll_\varepsilon N^\varepsilon\) eight-tuples.
For each forbidden prime at least one norm is anisotropic modulo \(p\);
\(p\mid a\) forces both coordinates in that norm to be divisible by \(p\),
and hence \(p^2\mid a\). The same holds for \(b\). Counting in the full
positive interval \(1\leq a\leq N-1\) gives at most
\(2\lfloor(N-1)/p^2\rfloor\leq2N/p^2\) relevant pairs, and restricting to
\(I\) cannot increase this count. Summing over integers instead of primes and
using \(\sum_{n>P}n^{-2}\leq1/(P-1)\leq2/P\) proves (5). The floor bound
avoids an unnecessary endpoint error for each prime. It includes singular
tuples and primes dividing \(N\).

This is an unconditional integer-counting bound. To make it negligible relative
to \(R_I(N)\), a sufficiently strong lower bound for \(R_I(N)\) is still needed.

## A sufficient counting theorem, with its exact hypotheses

**Conditional reduction C12.** Suppose there are constants \(c>0\) and
\(\theta>0\), such that for all sufficiently large admissible \(N\),

\[
R_I(N)\geq cN, \tag{H1}
\]

and the summed bad-count estimate holds **inside the same balanced interval**:

\[
\sum_{\substack{p\in\mathcal B\\p\leq N^\theta}}B_{p,I}(N)
\leq R_I(N)\sum_{\substack{p\in\mathcal B\\p\leq N^\theta}}\beta_p(N)
+o(N). \tag{H2}
\]

Then Target S12 holds. Hypotheses (H1) and (H2) are not proved here.

Indeed put \(P=N^\theta\). By (4), dividing (H2) by (H1) gives

\[
\frac1{R_I(N)}\sum_{\substack{p\in\mathcal B\\p\leq P}}B_{p,I}(N)
\leq\pi^2-9+o(1).
\]

By (5) and (H1), the remaining fraction is
\(O_\varepsilon(N^{\varepsilon-\theta})\).
Choose \(0<\varepsilon<\theta\). Both errors tend to zero, so the union of
all bad events occupies at most \(\pi^2-9+o(1)<1\) of the tuples. At least
\((10-\pi^2-o(1))R_I(N)\) tuples survive the union bound. The limiting margin
is approximately 13.04 percent of the balanced candidate count. Their values lie in
\(C\), and balance was built into \(R_I(N)\). This proves the implication.

No estimates for intersections of several bad-prime events are required. One
stronger sufficient version of (H2) is the individual estimate
\[
B_{p,I}(N)\leq\beta_p(N)R_I(N)+Kp^A N^{1-\eta}
\quad(p\leq N^\theta),\qquad \eta>\theta(A+1),
\]
with fixed \(A\geq0\), \(K>0\). Summing its errors gives
\(O(N^{1-\eta+\theta(A+1)})=o(N)\). Every fixed \(\theta>0\) is allowed.
Even the weaker bound \(N^\varepsilon(N/P+\sqrt N)\) would allow this by
choosing \(\varepsilon<\min(\theta,1/2)\); the sharper tail improves the
estimate without enlarging that exponent scope. Moreover, a forbidden prime
larger than \(\sqrt N\) cannot divide either positive summand.
Effective constants, effective validity thresholds, and an
effective error bound replacing the little-oh term would yield an effective
cutoff; the local margin alone does not.

A lower bound for \(R_{\mathrm{tot}}(N)\) does not establish (H1). Even knowing
\(R_I(N)\sim R_{\mathrm{tot}}(N)/3\) would not suffice alongside only a global
bad-count bound \((\pi^2-9+o(1))R_{\mathrm{tot}}(N)\): that bound can exceed
the entire balanced count. Both inputs of C12 concern \(R_I\) and \(B_{p,I}\).

## Pointwise upper bounds and a fixed-prime sieve interface

The [pointwise sieve note](SIMULTANEOUS_NORM_POINTWISE_SIEVE.md) supplies an
unconditional upper bound at the natural target-dependent scale. Put
\[
\mathcal T(N)=
\prod_{\substack{p\mid N\\p\equiv1\ (12)}}(1+9/p)
\prod_{\substack{p\mid N\\p\ge5,\ p\not\equiv1\ (12)}}(1+1/p).
\]
Theorems in that note prove
\[
R_{\mathrm{tot}}(N)\ll N\mathcal T(N),\qquad
B_p(N)\ll \frac N{p^2}\mathcal T(N)
\quad(p\le N^{1/8},\ p\text{ forbidden}),
\]
\[
\sum_{p>Y}B_p(N)
\ll \frac{N\mathcal T(N)}Y+N^{15/16}.
\]
They use Henriot's corrected 2014 New Theorem 5. The exact rescaling of a
bad component leads to dilated correlations
\(\sum_mh(m)h(N-p^2m)\); the proof includes primes dividing the target.
The predicted local singular product is uniformly comparable to
\(\mathcal T(N)\), so this scale has a uniform positive lower bound.
A positive singular product is not a proof of a global representation count.

An alternative sufficient interface is now precise: establish
\[
R_I(N)\gg N\mathcal T(N),\qquad
B_{p,I}(N)/R_I(N)-\beta_p(N)\longrightarrow0
\quad\text{for every fixed forbidden }p.
\tag{H*}
\]
Both assertions remain open. The marginal convergence must hold along all
admissible targets, including varying \(v_p(N)\). The upper tail then permits
fixing a prime cutoff before taking \(N\to\infty\), and gives the same
survival margin \(10-\pi^2\). This route avoids an asymptotic uniform in a
growing prime range. It uses a stronger, naturally normalized lower bound
than bare (H1), and does not prove the original additive-\(o(N)\) assertion
(H2). The note gives an analogous tail for centered \(N^\theta\) windows
when \(\theta>1/2\); their pointwise lower bounds and marginals remain open.

## What the global correlation actually is

The divisor formulas also identify the arithmetic weight explicitly. Set

\[
h(n)=r_G(n)r_E(n)/24.
\]

Let \(K=\mathbb Q(\zeta_{12})\). Its four Dirichlet characters are
\(1,\chi_{-4},\chi_{-3},\chi_{12}\). For \(\Re s>1\), direct Euler-factor
comparison gives

\[
\sum_{n\geq1}\frac{h(n)}{n^s}
=\frac{\zeta(s)L(s,\chi_{-4})L(s,\chi_{-3})L(s,\chi_{12})}
 {L(2s,\chi_{12})},\qquad\chi_{12}=\chi_{-4}\chi_{-3}. \tag{6}
\]

Equivalently, this is exactly \(\zeta_K(s)/L(2s,\chi_{12})\).

With \(t=p^{-s}\), the local factors are

| Prime | Euler factor for \(h\) |
|---|---|
| \(p\equiv1\pmod {12}\) | \((1+t)/(1-t)^3=\sum_{e\geq0}(e+1)^2t^e\) |
| \(p\equiv5,7\pmod {12}\) | \((1+t^2)/(1-t^2)^2=\sum_{j\geq0}(2j+1)t^{2j}\) |
| \(p\equiv11\pmod {12}\), or \(p=2,3\) | \(1/(1-t^2)\) |

There are no omitted ramified-prime corrections. At \(p=2\), the local factor
of \(\zeta_K\) is \((1-t)^{-1}(1+t)^{-1}=(1-t^2)^{-1}\), using
\(\chi_{-3}(2)=-1\); at \(p=3\), the same follows from
\(\chi_{-4}(3)=-1\). In both cases \(\chi_{12}(p)=0\).

Let \(a_K(n)\) count the nonzero integral ideals of \(K\) of norm \(n\).
Multiplying the absolutely convergent series in (6) gives the exact identities

\[
h(n)=\sum_{d^2\mid n}\mu(d)\chi_{12}(d)a_K(n/d^2),\qquad
a_K(n)=\sum_{d^2\mid n}\chi_{12}(d)h(n/d^2). \tag{9}
\]

Their supports coincide: every prime outside \(1\pmod {12}\), including 2
and 3, must occur to an even exponent. Their weights differ, for example
\(h(25)=3\), \(a_K(25)=2\), whereas \(h(121)=1\), \(a_K(121)=2\).
Consequently a correlation theorem for ideal coefficients must be transferred
through (9), including its congruences and square-divisor weights.

In particular,
\[
R_I(N)=576\sum_{\substack{a+b=N\\a/N\in I}}h(a)h(b),\qquad
R_{\mathrm{tot}}(N)=576\sum_{\substack{a+b=N\\a,b>0}}h(a)h(b).
\]
For primitive representations in both rings the normalized weight is instead
\(f(n)=4^{\omega(n)}1_C(n)\), where \(\omega(1)=0\). Each split prime has two
primitive choices in each ring and the units contribute 24. The exact balanced
primitive tuple count is consequently

\[
576\sum_{N/3\leq a\leq2N/3}f(a)f(N-a). \tag{7}
\]

The real density is also explicit: the two binary norm level densities are
\(\pi\) and \(2\pi/\sqrt3\). Integrating their product for the two summands
over \(a\in[N/3,2N/3]\) gives the balanced singular integral
\(4\pi^4N/9>0\). This supplies the archimedean factor, not an integer count.

These identities describe a concrete binary shifted correlation. A one-variable
Dirichlet series does not prove the positivity of its individual additive
coefficients. In particular (6) cannot be substituted for (H1) or (H2).

Nor is the system a generic nonsingular intersection of three quadrics. For
\(Q_1=G(X)-E(U)\), \(Q_2=G(Z)-E(V)\), \(Q_3=G(X)+G(Z)\), its pencil has

\[
\det(\lambda Q_1+\mu Q_2+\nu Q_3)
=\frac9{16}\{\lambda\mu(\lambda+\nu)(\mu+\nu)\}^2. \tag{8}
\]

The determinant discriminant vanishes, and the pencil has minimum rank four.
The affine fibre also has singular boundary surfaces at \(a=0\) or \(b=0\).
Our real balanced region avoids those surfaces but does not repair the failed
hypotheses of a generic circle-method theorem. Rational parametrization gives
denominators; clearing them replaces \(N\) by \(d^2N\), so it does not solve
the integral fixed-target problem.

## Why the one-ring shortcut reaches a sieve barrier

Starting with an Eisenstein shell and excluding prime factors
\(p\equiv7\pmod{12}\) produces a nonsummable local loss.
The Eisenstein form is split modulo this prime. Write the determinant shell
as \(xy+uv=N\) with \(p\nmid N\). It contains
\(p(p-1)(p+1)\) points. The events \(xy=0\) and \(uv=0\) are disjoint,
and together contain \(2(2p-1)(p-1)\) points. Hence the forbidden fraction is
\[
\frac{4p-2}{p(p+1)}\sim\frac4p.
\]
This prime class has density \(1/4\), so its sieve dimension is one.
After imposing component values \(1\pmod{12}\), sieving these primes through
\(\sqrt N\) excludes every such prime factor: a larger one would force a
smaller one by the residue class of the remaining factor. But the classical
lower linear sieve at a natural level \(D\le N\) then has
\(\log D/\log\sqrt N\le2\), where its lower sieve function gives no positive
main term. This diagnoses that route; it is not an impossibility theorem for
stronger analytic input.

A seemingly simpler quadratic replacement is also incorrect.
\(133=11^2+12\) is represented by \(x^2+12y^2\), but
\(133=7\cdot19\) has odd Gaussian-inert exponents, so \(h(133)=0\).
The theta series of that binary form cannot replace the simultaneous weight
\(h\) in the pointwise correlation.

A fixed finite major-arc calculation also needs additional cancellation.
Let \(h_6(n)=h(n)\mathbf1_{(n,6)=1}\). The Frobenian mean theorem
[Loughran--Matthiesen, Lemma 2.8](https://arxiv.org/html/1904.12845v2#S2)
applies both to \(h_6\) and to \(h_6^2\). Their prime means are respectively
\(1\) and \(4\), since the split-prime values are \(4\) and \(16\),
on a prime set of density \(1/4\). Thus
\[
\sum_{n\le X}h_6(n)\asymp X,\qquad
\sum_{n\le X}h_6(n)^2\asymp X(\log X)^3.
\]
For \(F_X(\alpha)=\sum_{n\le X}h_6(n)e(n\alpha)\), Parseval gives the
second quantity as its squared \(L^2\) norm, while \(\|F_X\|_\infty\ll X\).
Any fixed finite collection of rational arcs of total length \(O_Q(X^{-1})\)
therefore captures at most \(O_Q(X)\) of that energy. The complementary
absolute integral remains \(\gg X(\log X)^3\). Bounding the binary minor-arc
coefficient only by this absolute integral cannot prove an \(o(X)\) error.
A successful circle-method approach needs signed cancellation, growing arc
resolution or another stronger input. This rules out that particular
fixed-arc absolute-value estimate, not the circle method itself.

The [character-divisor obstruction](SIMULTANEOUS_NORM_DIVISOR_OBSTRUCTIONS.md)
addresses another natural approximation. Put \(k=h_6*\mu\) and
\(h_D(n)=\sum_{d\mid n,d\le D}k(d)\). Every fixed cutoff \(D=X^\gamma\),
\(0<\gamma<1\), leaves \(L^1\) error \(\gg_\gamma X\) and squared
\(L^2\) error \(\gg_\gamma X(\log X)^3\) on \(X/2<n\le X\), even
restricted to \(n\equiv1\pmod{12}\). At \(n=189805\),
\(h_6(n)=0\) but \(h_{\sqrt n}(n)=4\). The truncated sum is therefore
not a pointwise minorant. These facts rule out the stated absolute
approximation steps, not a refined character-divisor method with cancellation.

## An unconditional global count averaged over targets

The quartic connection does give a balanced count after averaging the target.
The following is a standard consequence of one-variable ideal counting, not a
pointwise substitute for either (H1) or (H2).

**Averaged theorem A12.** Put
\[
\kappa_K=\operatorname*{Res}_{s=1}\zeta_K(s),\qquad
\rho=\frac{\kappa_K}{L(2,\chi_{12})}>0.
\]
As \(X\to\infty\) through positive real values,
\[
\sum_{\substack{X<N\leq2X\\N\equiv2\ (12)}}R_I(N)
=128\rho^2X^2+O(X^{7/4}). \tag{10}
\]
The field is fixed, and the implied constant may depend on it.

The precise imported estimate is [Debaene, arXiv:1611.10103v2,
Corollary 2](https://arxiv.org/html/1611.10103v2): for a fixed degree-four
number field and every \(x\geq1\),
\[
A_K(x):=\sum_{n\leq x}a_K(n)=\kappa_Kx+O_K(x^{3/4}).
\]
It applies to all integral ideals, not only principal ideals. Apply (9) and sum:
\[
\begin{split}
H(x):=\sum_{n\leq x}h(n)
&=\sum_{d\leq\sqrt x}\mu(d)\chi_{12}(d)A_K(x/d^2)\\
&=\rho x+O\left(x^{3/4}\sum_{d\geq1}d^{-3/2}+\sqrt x\right)
=\rho x+O(x^{3/4}).
\end{split}
\]
The \(O(\sqrt x)\) term extends the main coefficient sum to infinity.

Set \(h_6(n)=h(n)\mathbf1_{(n,6)=1}\), and \(H_6(x)=\sum_{n\leq x}h_6(n)\).
Removing the two Euler factors \((1-p^{-2s})^{-1}\) gives exactly
\[
H_6(x)=H(x)-H(x/4)-H(x/9)+H(x/36)
=\rho_6x+O(x^{3/4}),\qquad \rho_6=\tfrac23\rho.
\]
Here \(H(y)=0\) for \(0\leq y<1\); the same estimates can be used on this
bounded range after enlarging their constants. A supported value of \(h_6\)
is automatically \(1\pmod {12}\).

Summing over \(N\) now removes the fixed additive constraint. For a given \(a\),
the allowed \(b\) lie between
\[
L(a)=\max(a/2,X-a),\qquad U(a)=\min(2a,2X-a).
\]
Thus the left side of (10), divided by 576, is exactly
\[
\sum_{X/3\leq a\leq4X/3}h_6(a)
\bigl(H_6(U(a))-H_6(L(a))\bigr).
\]
The strict lower target bound is \(b>X-a\). The other lower bound is
\(b\geq a/2\); it causes no missing endpoint because a supported \(a\) is
odd, so \(a/2\) is not an integer. No extra congruence restriction is needed:
two supported values already sum to \(2\pmod {12}\).

The inner errors contribute \(O(X^{3/4}H_6(4X/3))=O(X^{7/4})\).
Writing \(g=U-L\), its two pieces are
\[
g(a)=\begin{cases}
3a-X,&X/3\leq a\leq2X/3,\\
2X-\tfrac32a,&2X/3\leq a\leq4X/3.
\end{cases}
\]
It vanishes at both endpoints, has total variation \(2X\), and integral
\(X^2/2\). Abel summation using the estimate for \(H_6\) therefore gives
\[
\sum_a h_6(a)g(a)=\rho_6X^2/2+O(X^{7/4}).
\]
Multiplication by \(576\rho_6\) proves (10), since
\(288\rho_6^2=128\rho^2\).

The constant can also be evaluated:
\[
\rho=\frac{\log(2+\sqrt3)}{2\sqrt3},\qquad
128\rho^2=\frac{32}{3}\log^2(2+\sqrt3).
\]
Indeed \(L(1,\chi_{-4})=\pi/4\), \(L(1,\chi_{-3})=\pi/(3\sqrt3)\), and
\[
L(1,\chi_{12})=\int_0^1\frac{1-x^4}{1+x^6}\,dx
=\frac{\log(2+\sqrt3)}{\sqrt3}.
\]
For the last integral, differentiate the logarithm of
\((x^2+\sqrt3x+1)/(x^2-\sqrt3x+1)\). The classical partial-fraction expansion
of \(\pi^2\csc^2(\pi z)\), grouping the four nonzero residue classes modulo
12, gives
\[
L(2,\chi_{12})=\frac{\pi^2}{144}
\left(\csc^2\frac\pi{12}-\csc^2\frac{5\pi}{12}\right)
=\frac{\pi^2}{6\sqrt3}.
\]
These values give \(\kappa_K=\pi^2\log(2+\sqrt3)/36\) and the stated \(\rho\).

The averaged main term is positive. It establishes neither the pointwise lower
bound (H1) nor the bad-prime estimate (H2). Those coordinate-weight inputs remain open even though binary coverage now follows by the separate subgroup theorem.

## The averaged coordinate count after the full sieve

Set \(g(n)=h(n)\mathbf1_C(n)\). Thus
\[
R_I^{\mathrm{good}}(N)=576\sum_{\substack{a+b=N\\a/N\in I}}g(a)g(b)
\]
counts coordinate tuples whose two values lie in \(C\). Their own coordinates
need not be primitive at allowed split primes. Define
\[
D(s)=
\prod_{p\equiv5,7\ (12)}\frac{(1-p^{-2s})^2}{1+p^{-2s}}
\prod_{p\equiv11\ (12)}(1-p^{-2s}),\qquad D=D(1).
\]

**Sieved averaged theorem A12G.** With the same fixed field and interval as A12,
\[
\sum_{\substack{X<N\leq2X\\N\equiv2\ (12)}}R_I^{\mathrm{good}}(N)
=128\rho^2D^2X^2+O(X^{7/4})
=\frac{32}{3}\log^2(2+\sqrt3)D^2X^2+O(X^{7/4}). \tag{14}
\]
This is another consequence of the same ideal-counting input used in A12.

To prove it, first remove 2 and 3 as in \(h_6\), then remove every remaining
forbidden Euler factor. The table following (6) gives
\[
\sum_{n\geq1}\frac{g(n)}{n^s}
=D(s)\sum_{n\geq1}\frac{h_6(n)}{n^s}\qquad(\Re s>1).
\]
There are no further factors at 2 and 3. Writing
\(D(s)=\sum_{d\geq1}\delta_D(d)d^{-2s}\), the local expansions are
\[
\frac{(1-q)^2}{1+q}=1-3q+4q^2-4q^3+\cdots,
\qquad 1-q,
\]
where \(q=p^{-2s}\), for the two respective prime types. At \(\Re s=3/4\)
their absolute coefficient sums are
\(1+3p^{-3/2}+4p^{-3}/(1-p^{-3/2})\) and \(1+p^{-3/2}\).
Their product converges, so
\[
\sum_{d\geq1}\frac{|\delta_D(d)|}{d^{3/2}}<\infty.
\]
In fact \(D(s)\) is absolutely convergent throughout \(\Re s>1/2\), and
\(D>0\).

The exact square-divisor convolution is
\[
g(n)=\sum_{d^2\mid n}\delta_D(d)h_6(n/d^2).
\]
Using the summatory estimate for \(H_6\) from A12 gives
\[
\begin{split}
S_g(x):=\sum_{n\leq x}g(n)
&=\sum_{d\leq\sqrt x}\delta_D(d)H_6(x/d^2)\\
&=\rho_6Dx+O(x^{3/4}).
\end{split}
\]
The error sum is bounded by \(O(x^{3/4}\sum_d|\delta_D(d)|d^{-3/2})\).
Extending the main coefficient sum to infinity costs at most a constant times
\[
x\sum_{d>\sqrt x}\frac{|\delta_D(d)|}{d^2}
\leq x^{3/4}\sum_{d>\sqrt x}\frac{|\delta_D(d)|}{d^{3/2}}.
\]
The balanced-region proof of A12 now applies with \(g\) in place of \(h_6\).
Its support is still odd and \(1\pmod {12}\), so the same endpoint argument
applies. The leading constant is multiplied by \(D^2\), and the error remains
\(O(X^{7/4})\), proving (14).

The [exact interval checker](verify_simultaneous_norm_survival.py) and its
[receipt](simultaneous_norm_survival_checks.json) certify
\[
0.633187<D^2<0.633191.
\]
For reproducibility, it multiplies the squared factors through primes \(10^6\)
with directed integer rounding at scale \(10^{30}\). For each omitted prime
the negative logarithm of its factor is at most \(6/(p^2-1)\). Summing over
all omitted integers yields a tail product at least
\(1-3(10^{-6}+(10^6+1)^{-1})\). Thus the displayed interval is rigorous,
and \(D^2\) rounds to 0.63319 at five decimal places.

Consequently
\[
\frac{\sum_{X<N\leq2X,\ N\equiv2\ (12)}R_I^{\mathrm{good}}(N)}
 {\sum_{X<N\leq2X,\ N\equiv2\ (12)}R_I(N)}
=D^2+O(X^{-1/4}).
\]
Approximately 63.3 percent of the weighted coordinate count survives on this
dyadic average. This is not a fraction of targets, a per-target survival
fraction, or a count restricted to primitive coordinates. The unweighted
successful-pair conclusion is the separate corollary P12.

## Exact finite check

The [checker](verify_simultaneous_norm_sums.py) and its
[receipt](simultaneous_norm_sums_checks.json) cover every target
\(N\equiv2\pmod {12}\), \(2\leq N\leq10^8\): 8,333,334 targets, with last
target 99,999,998. There are no unbalanced exceptions. The balanced exceptions
are exactly

\[
14,\ 38,\ 50,\ 62,\ 86,\ 230,\ 242.
\]

For example \(230=1+229\) and \(242=1+241\) satisfy the unbalanced statement.
All seven unbalanced witnesses are retained in the receipt.

Membership is exact without sieving to \(10^8\): an integer
\(n\equiv1\pmod {12}\) belongs to \(C\) iff it has no forbidden prime factor
\(p\leq\sqrt n\). If all smaller factors are allowed, there can be at most one
remaining forbidden factor, of exponent one. Its residue must then be 1 because
\(n\) and the other factors are 1; this is a contradiction. Thus sieving the
compressed values \(12i+1\) by forbidden primes up to the square root suffices.

Writing \(N=12j+2\), test \(a=12i+1\), \(b=12(j-i)+1\). Up to swapping,
balance is exactly \(\lceil j/3\rceil\leq i\leq\lfloor j/2\rfloor\).
The checker exhausts these indices until it finds a witness; only for a balanced
failure does it search the full interval \(0\leq i\leq\lfloor j/2\rfloor\).
It validates every found balanced witness, independently factors all 8,334
represented integers up to 100,000, rejects known nonmembers, and compares the
vector search with direct scalar enumeration through 5,000.

Replay from the repository root with Python 3.12.3 and NumPy 1.26.4:

```sh
python3 research/arithmetic/verify_simultaneous_norm_sums.py --limit 100000000
python3 research/arithmetic/verify_simultaneous_norm_local.py
python3 research/arithmetic/verify_simultaneous_norm_survival.py
```

The receipt pins checker bytes and deterministic membership/witness hashes. This
is internal computational evidence; no external review or formal proof is claimed.

## The longer-term fixed-modulus question

The central family statement is
[FBS-U](BINARY_SUMS_WITH_PRESCRIBED_PRIME_FACTORS.md#the-central-binary-conjecture).
At modulus 12 it conjectures, for every fixed \(0<\theta<1\) and every
sufficiently large \(N\equiv2\pmod {12}\),
\[
r_{12,N^\theta}(N)\gg_\theta
\frac{N^\theta}{(\log N)^{3/2}},
\]
where the count is unweighted and the two values differ by at most
\(N^\theta\). This remains open. It implies the weaker
[macroscopic family goal FBS-G](BINARY_SUMS_WITH_PRESCRIBED_PRIME_FACTORS.md#the-two-binary-research-campaigns),
whose modulus-12 instance is S12. The two quadratic characters give a specific
norm model here; the local proof and conditional reduction do not supply the
shrinking-window quantitative estimates or the inputs for arbitrary moduli.

## Literature and the remaining obligation

- Ramaré--Viswanadham, [arXiv:2509.23260v1](https://arxiv.org/html/2509.23260v1),
  Introduction and Theorem 1.6: explicitly includes the simultaneous sequence in
  its sieve framework; its stated additive theorem is ternary Gaussian. The
  positive-relative-density subset theorem cannot simply take \(C\), whose
  relative density in the Gaussian sequence tends to zero.
- Malavika--Ramaré, [arXiv:2604.09448v1](https://arxiv.org/html/2604.09448v1),
  Theorems 1.1--1.2 and the introduction's Loeschian discussion: estimates
  quadratic-phase exponential sums and mentions simultaneous norms of sieve
  dimension \(3/4\). This is not the fixed-target binary counting input above.
- Iwaniec--Pomykała, [*Sums and differences of quartic norms* (1993)](sources/iwaniec-pomykala-quartic-norms-1993.pdf):
  the full paper is now retained and its main theorems audited below.
  Theorem 1 is unconditional for differences; Theorem 2 gives conditional
  every-target sums of coprime completely splitting integers, with
  \(N/(\log N)^{3/2}\)-order bounds for the quartic case.
  The latter specializes exactly to \(C_{12}+C_{12}\), assuming the
  stated growing-shift uniformity. It is a direct conditional antecedent.
- Matthiesen, [*Linear correlations amongst numbers represented by positive
  definite binary quadratic forms*](https://arxiv.org/html/1106.4690v3),
  Theorem 1.1: requires pairwise
  affine-independent forms. Repeated arguments \((a,a,N-a,N-a)\) fail this
  hypothesis; it does not establish the needed binary correlation.
- Pierce--Schindler--Wood, [arXiv:1509.04757v2](https://arxiv.org/pdf/1509.04757v2),
  Theorem 1.1: requires at least ten variables and a nonvanishing determinant
  discriminant, and concerns almost all target triples. Neither the hypotheses
  nor that conclusion give the present eight-variable fixed line of targets.
- Dietmann--Elsholtz, [*Sums of two numbers having only prime factors congruent
  to one modulo four*](https://www.math.tugraz.at/~elsholtz/WWW/papers/greaves20.pdf),
  already proves sufficiently-large Gaussian binary coverage using Greaves.
  The new Gaussian paper's explicit all-positive result should be distinguished
  from this earlier asymptotic theorem.

Unrestricted unweighted pointwise counting for S12U and macroscopic balanced S12 is now proved by the [quantitative subgroup theorem](SPLITTING_SEMIGROUP_LOCAL_GLOBAL.md). The following sieve routes retain their own unresolved inputs.
The separate [half-dimensional criterion](ASYMMETRIC_SEMILINEAR_SIEVE.md)
gives another sufficient route: control one signed remainder beyond level
\(N^{1/2}\) for the positive asymmetric sequence. Its local density,
parity reduction and exact kernel compression are established, with
sifted model mass \(\gg N/\sqrt{\log N}\). The required global
remainder bound is open.
A further [critical switching bound](ASYMMETRIC_CRITICAL_SWITCHING.md)
now reduces the sufficient distribution obligation to square-root level
up to logarithms. Its upper bound for the remaining two-large-prime mass
holds on every admissible target, including arbitrary good target factors.
The critical distribution estimate remains open; the two-sided conclusion is now supplied by the separate subgroup-correlation argument.
The balanced inputs (H1)--(H2), or the naturally normalized fixed-prime
interface (H*), remain sufficient routes to S12. None has been established
on every target. SimNormPower establishes both global inputs outside a
power exceptional set, including the full moving balanced count and all
fixed-prime target valuations. Eliminating the remaining exceptions requires
further pointwise information. The local margin, sharper large-prime tail, and averaged
counts before and after the full sieve have complete arguments with the stated
ideal-counting input. The linked general theory now gives unweighted
semigroup-scale density one for \(\theta>1/3\), fixed-width variance,
quantitative full-count variance for \(5/6<\theta<1\), and exact first
moments in every fixed power window. The sharp four-form energy is proved in
that upper range; it remains an obligation for \(0<\theta\le5/6\). The required
every-target power-short binary estimates, an effective cutoff, and the uniform
semigroup-scale family conjecture remain open here. D12Power supplies
power-saving exceptions for the lower count at \(H=X^{99/100}\); a
power-saving full-count asymptotic and elimination of every exception require
further arguments.

### Full Iwaniec--Pomykała theorem audit

The full thirteen-page publisher PDF of [Iwaniec--Pomykała,
*Sums and differences of quartic norms*, Mathematika 40 (1993),
233--245](sources/iwaniec-pomykala-quartic-norms-1993.pdf) has been inspected.
The comparison is with their actual main theorems, not merely their
auxiliary distribution proposition.

**Arithmetic identification.** Their degree-four abelian extension
\(K/\mathbb Q\), of conductor \(\Delta\), determines an index-four
subgroup \(\mathcal H\le(\mathbb Z/\Delta)^\times\). The summands are
integers all of whose prime factors split completely in \(K\).
They count each numerical pair once, independently of how many ideals
have those norms; see pp. 233--234 and the excluded-prime sieve on p. 239.
For
\[
 K=\mathbb Q(\zeta_{12}),\qquad \Delta=12,\qquad \mathcal H=\{1\},
\]
this is exactly \(C_{12}\). Their local sum condition is
\(N\equiv2\pmod{12}\), and the difference condition is
\(N\equiv0\pmod{12}\). The additional parity condition (2) on p. 233
is automatic on these classes.

The simultaneous primitive Gaussian and Eisenstein norm condition is
also exactly \(C_{12}\): primitiveness excludes inert prime factors in
each quadratic field; 2 is inert in the Eisenstein field, 3 in the
Gaussian field, and all remaining allowed primes are \(1\bmod12\).
Conversely those primes generate primitive norms in both class-number-one
rings. Ordinary quartic norms would be a larger set: for example,
\(N_{K/\mathbb Q}(2+i)=25\), with primitive coordinates
\((2,0,0,1)\) in the integral basis \(1,\zeta_{12},\zeta_{12}^2,\zeta_{12}^3\).
That distinction does not remove the historical overlap: the paper
expressly imposes complete splitting of the prime factors.

**Published conclusions, printed p. 234.** Put
\[
 W_\Delta(N)=\prod_{p\mid N}(1-\varepsilon_p/p),\qquad
 \varepsilon_p=
 \begin{cases}
 1,&p\bmod\Delta\in\mathcal H,\\
 -1,&p\bmod\Delta\notin\mathcal H.
 \end{cases}
\]
Their statements concern coprime summands.

| Result | Hypotheses and conclusion |
|---|---|
| Theorem 1 | Unconditionally, every fixed positive admissible difference \(N\) has infinitely many coprime complete-splitting pairs \(a-b=N\). The count with \(a\le X\) has order \(W_\Delta(N)X/(\log X)^{3/2}\) as \(X\to\infty\). The comparison constants depend only on \(\Delta\); no onset uniform in growing \(N\) is supplied by the fixed-shift argument. |
| Theorem 2 | Assuming their (P4) uniformly for all shifts \(|a|<X\), every sufficiently large admissible sum target \(N\) has coprime complete-splitting pairs, of order \(W_\Delta(N)N/(\log N)^{3/2}\). The statement is explicitly conditional. |
| Quantitative scope | These are order bounds, not a leading asymptotic with a specified constant. The authors prove the lower bound and cite a Brun-sieve upper bound. Their theorem does not prescribe a central balance interval. An interval adaptation must be checked before claiming novelty for a balanced refinement. |

The p. 234 remarks also discuss prescribing ideal classes via genus
fields, with additional local conditions. That remark is an antecedent
for the arithmetic motivation; it is not a separately proved
all-order primitive-value leading asymptotic.

In particular their Theorem 2 already gives, under its uniformity
hypothesis,
\[
 \#\{a+b=N:a,b\in C_{12},\ (a,b)=1\}
 \asymp
 \left(\prod_{\substack{p\mid N\\p\equiv1(12)}}(1-p^{-1})\right)
 \left(\prod_{\substack{p\mid N\\p\not\equiv1(12)}}(1+p^{-1})\right)
 \frac{N}{(\log N)^{3/2}},
 \qquad N\equiv2\pmod{12}.
\]
The factor at 2 is fixed, and 3 never divides an admissible sum target.

**The exact analytic hypothesis.** The Proposition on p. 243 assumes
a complex sequence \(\beta_n\), supported on \(n\sim Y\), satisfying
\(|\beta_n|\le\tau(n)\), with no prime factor at most
\(\exp((\log\log n)^2)\). Its (P3) requires
\[
 \sum_{\substack{n\sim Y\\n\equiv l\ (k)}}\beta_n
 -\frac1{\varphi(k)}
       \sum_{\substack{n\sim Y\\(n,k)=1}}\beta_n
 \ll_B Y(\log Y)^{-B}
\]
for every \(B>0\), uniformly over \(k\ge1\) and \((l,k)=1\).
Its other three coefficient sequences are arbitrary bounded sequences. Its (P4) bounds the weighted bilinear-modulus
convolution discrepancy by
\[
 O_{\epsilon,a,B}\bigl(X(\log X)^{-B}\bigr)
\]
under the factor ranges
\[
 X^\epsilon(1+X^{-1}Q^2)R<Y
       <X^{5/6-\epsilon}(Q^2R)^{-2/3}.                 \tag{IP.1}
\]
The published unconditional proposition allows dependence on the
fixed nonzero shift \(a\). Theorem 2 needs that dependence removed
for target-sized shifts. On p. 244 the authors trace it to the
Maass Hecke coefficients and obtain uniformity for \(|a|<X\) under
\(|\lambda_j(a)|\le\tau(|a|)\), their stated Ramanujan hypothesis.
The distinction between fixed shift and growing target is the
reason the unconditional difference theorem does not prove the
sum theorem.

The paper also states an unconditional uniform alternative from
Fouvry, but its factor range yields only a parameter \(1/28\),
where their successful quartic calculation uses \(1/10\).
The numerical discussion on pp. 242--243 reports that even \(1/11\)
does not make that same sieve lower bound positive. Neither a
fixed-shift estimate nor an unspecified polynomial dependence
therefore completes their sum argument. The reported rigorous
positivity computation was not independently replayed in this
source comparison.

**What the later auxiliary statements do and do not supply.**
[Li, arXiv:2602.20917v6](https://arxiv.org/abs/2602.20917v6),
Section 4.1.1, Lemma 4.3, PDF p. 49, reproduces the auxiliary
Proposition in a setting where the nonzero shift is fixed.
Wu's Lemma 2.5, printed p. 222, permits only
\(|a|\le(\log X)^A\). Neither inspected statement directly supplies
the growing-target version assumed in Theorem 2.

A literal substitution \(\beta=\mathbf1_{C_{12}}\) also fails the
auxiliary hypotheses: its entire dyadic mass \(T\) is in the
residue 1 modulo 12, giving discrepancy \(3T/4\), and the allowed
factor 13 violates the rough support condition. These are
obstructions to that substitution. They do not invalidate the
original theorem or distinguish its arithmetic scope; the original
proof constructs different weights.

A relevant modern comparison is [Lichtman,
arXiv:2309.08522v1](https://arxiv.org/abs/2309.08522v1),
Proposition 3.1, p. 9, and Theorem 1.8, pp. 4--5.
It does give a convolution estimate uniform in a growing shift,
but requires triply well-factorable modulus weights, a rough
Siegel--Walfisz factor and its own factor ranges. Its general
Kloosterman estimate retains an explicit shift-dependent spectral
term. A replacement of (IP.1) requires a weight-factorization,
range and positivity argument; quoting Kim--Sarnak or a modern
distribution level alone does not supply one. This checks the
named statements only, not every later refinement of the 1993
method.

**Positioning consequence.** The 1993 main-theorem comparison is
complete at the statement and hypothesis level. The quartic
complete-splitting binary problem, its admissibility conditions,
conditional every-target sum theorem and logarithmic order are
established antecedents. The unconditional difference theorem must
also be credited.

The present [quantitative correlation theorem](SPLITTING_SEMIGROUP_LOCAL_GLOBAL.md)
asserts an unrestricted-pair leading asymptotic with an explicit
constant and logarithmic remainder, using its separately recorded
analytic inputs. Its [coprime corollary](SPLITTING_SEMIGROUP_LOCAL_GLOBAL.md#coprime-summands)
now recovers the older coprime-pair scope by Möbius inversion, with
leading factor \(A_{12}J_{12}(N)P_{12}(N)\) and relative error
\(O((\log N)^{-1/65536})\), where
\(P_{12}(N)=\prod_{p\mid N,\ p\equiv1(12)}(1-p^{-1})\).
The proof explicitly controls the loss from this factor tending to zero.
For the 1993 target weight \(W_{12}(N)=\prod_{p\mid N}(1-\varepsilon_p/p)\),
the factors agree at allowed primes. Since an admissible target is even
and not divisible by three,
\[
 \frac{J_{12}(N)P_{12}(N)}{W_{12}(N)}
 =\frac23\prod_{\substack{p\mid N\\p\equiv5,7,11\ (12)}}
             \frac{p(p-1)}{(p-2)(p+1)}.
\]
Every factor in this product is \(1+O(p^{-2})\), so the ratio has
positive uniform upper and lower bounds. The corollary therefore
recovers their logarithmic order for the same coprime numerical count,
adding a leading constant and remainder through the present analytic
dependencies rather than their hypothesis (P4).
This audit establishes neither a first unconditional result nor
exhaustive historical priority. A later removal of the 1993
uniformity hypothesis within their method and an adapted balanced
statement still require their own verified arguments.

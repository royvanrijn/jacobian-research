# Binary sums with prescribed prime factors

The central research target is a uniform binary lower bound at the scale of the
whole multiplicative semigroup, on every compatible large target. That is the
open conjecture FBS-U below. Positive-proportion coverage, semigroup-scale density one in a range of power
windows, full-count variance in fixed-width windows and a range of power
windows, and an every-target ternary result are established. This note and its linked proofs derive these results from published
tools and keep the remaining binary estimates and their priority assessment
separate.

The [simultaneous Gaussian--Eisenstein case study](SIMULTANEOUS_GAUSSIAN_EISENSTEIN_SUMS.md)
keeps the exact norm identities, lifted local densities, conditional pointwise
criterion and finite certificates for the class \(1\pmod {12}\). The immediate
arithmetic target is every-target existence for \(C_{12}+C_{12}\), initially
without balance. In parallel, the
[Eisenstein power-window theorem](EISENSTEIN_POWER_SHRINKING_WINDOWS.md)
now covers every fixed \(2/3<\theta<1\) for \(q=3\), with explicit
cutoffs and both ordinary and mixed families.
The separate [asymmetric pointwise theorem](ASYMMETRIC_GAUSSIAN_EISENSTEIN_SUMS.md)
now gives a power-error weighted asymptotic for every sufficiently large
\(N\equiv2\pmod6\), with one summand in \(C_{12}\) and the other
in \(A_3\). Its target-dependent main term is explicit and positive.
Its [uniform coefficient extension and cofactor argument](ASYMMETRIC_SEMIGROUP_SCALE.md)
now prove the unweighted order \(N J_E(N)/(\log N)^{5/4}\) on every
sufficiently large compatible target. It proves no every-target
\(C_{12}+C_{12}\) statement.

The [fixed-modulus theory](ASYMMETRIC_FIXED_MODULUS_THEORY.md) now gives an
asymmetric theorem for every fixed \(q\ge3\). Put
\(L=\operatorname{lcm}(2,q)\), choose a nonprincipal quadratic character
\(\chi\bmod L\), and let \(A_\chi\) contain the positive integers whose
prime factors satisfy \(\chi(p)=1\). On every sufficiently large
\(N\equiv2\pmod L\), the unweighted \(C_q+A_\chi\) count has order
\(N J_\chi(N)/(\log N)^{3/2-1/\varphi(q)}\), with the explicit positive
target factor in that proof. A separate clause reaches the full
\(C_q+C_q\) logarithmic scale while allowing only boundedly many excluded
prime factors in the second summand, all above some fixed power of \(N\).
Removing those remaining factors is still necessary for full two-sided
coverage when \(\varphi(q)>2\).

## The central binary conjecture

For a fixed \(q\geq3\), put
\[
C_q=\{n\geq1:p\mid n\Longrightarrow p\equiv1\pmod q\},\qquad
L=\ell_q=\operatorname{lcm}(2,q),\qquad
\kappa_q=2-\frac2{\varphi(q)}.
\]
Let \(r_{q,H}(N)\) count ordered positive pairs \(a,b\in C_q\) with
\(a+b=N\) and \(|a-b|\leq H\).

**Uniform semigroup-scale conjecture FBS-U.** For every fixed \(q\geq3\) and
every fixed \(0<\theta<1\), there exist \(c_{q,\theta}>0\) and
\(N_0(q,\theta)\) such that every \(N\geq N_0(q,\theta)\) with
\(N\equiv2\pmod L\) satisfies
\[
r_{q,N^\theta}(N)\geq
c_{q,\theta}\frac{N^\theta}{(\log N)^{\kappa_q}}. \tag{U}
\]
All parameters are fixed before the target grows. The lower-bound constant is
independent of the individual target. No uniformity in growing \(q\) or varying
\(\theta\) is proposed here. At \(q=12\), the denominator is
\((\log N)^{3/2}\).

This asks for many nearly equal decompositions of every compatible large
integer. It is stronger than an existence assertion and stronger than the
macroscopic balanced family goal FBS-G below. The full family remains a
conjecture. The [Eisenstein cofactor theorem](EISENSTEIN_POINTWISE_SEMIGROUP_SCALE.md)
establishes (U) at \(q=3\) for every fixed \(2/3<\theta<1\), with a
matching upper bound and the correct target-dependent singular factor.


## The general statement

Let \(\mathcal P\) be a fixed Frobenian set of rational primes. Concretely,
outside a finite exceptional set, membership is specified by a union of
conjugacy classes in the Galois group of a fixed finite Galois extension of
\(\mathbb Q\). Membership at the exceptional primes is fixed as well. Suppose
its prime density is \(\alpha>0\), and put
\[
C(\mathcal P)=\{n\geq1:p\mid n\Longrightarrow p\in\mathcal P\},
\qquad \xi(n)=\mathbf1_{C(\mathcal P)}(n).
\]
In particular \(1\in C(\mathcal P)\). Write \(\kappa=2-2\alpha\), so
\(0\leq\kappa<2\). For positive real \(H\), define the unweighted ordered-pair
count
\[
r_{\mathcal P,H}(N)=
\#\{(a,b)\in C(\mathcal P)^2:a+b=N,\ |a-b|\leq H\}.
\]

**Power-window theorem FBS1.** For every fixed \(0<\theta<1\), there are
\(c_1,c_2>0\) and \(X_0\), depending on \(\mathcal P\) and \(\theta\), such
that, with \(H=X^\theta\), every real \(X\geq X_0\) satisfies
\[
\#\left\{N\in\mathbb Z:X<N\leq2X,\quad
r_{\mathcal P,H}(N)\geq c_1\frac{H}{(\log X)^\kappa}\right\}
\geq c_2X. \tag{F1}
\]

Thus a positive proportion of targets have many nearly equal decompositions.
In particular, after changing the constant, these targets have
\(\gg_{\mathcal P,\theta}N^\theta/(\log N)^\kappa\) pairs with
\(a,b=N/2+O(N^\theta)\). All data and \(\theta\) are fixed before
\(X\to\infty\). There is no uniform assertion as the field, exceptional
primes, or window exponent varies. Polylogarithmic windows are outside the
statement.

For a general Frobenian set, (F1) describes a positive proportion of all integer
targets. It neither identifies all admissible residue classes nor gives a
positive proportion in every such class separately. In the congruence examples
below, every successful target automatically lies in the indicated fixed
progression, giving positive-proportion coverage within it as well.

## The imported tools

The indicator \(\xi\) is completely multiplicative on positive integers and
bounded by one. Its restriction to primes is Frobenian of mean \(\alpha\).
The growth and prime-power conditions of Definitions 2.6--2.7 in
[Loughran--Matthiesen, arXiv:1904.12845v2](https://arxiv.org/html/1904.12845v2)
hold with \(H=1\). Their Lemma 2.8, and Lemma 2.4(2), respectively give
\[
\sum_{n\leq x}\xi(n)
=c_{\mathcal P}x(\log x)^{\alpha-1}
+O_{\mathcal P}\bigl(x(\log x)^{\alpha-2}\bigr),
\qquad c_{\mathcal P}>0, \tag{F2}
\]
\[
\sum_{p\leq x}\frac{\xi(p)}p
=\alpha\log\log x+O_{\mathcal P}(1). \tag{F3}
\]
Only a long-interval mean is used for the lower bound.

The upper-bound input is
[Henriot, arXiv:1102.1643v1, Theorem 5 and Corollary 2](https://arxiv.org/html/1102.1643v1#S3)
(Corollary 2 is on preprint p. 7), read together with the
[2014 erratum, pp. 375--377](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/B0BD208D979495FE06B8192209E710EA/S0305004114000280a.pdf/nairtenenbaum_uniform_with_respect_to_the_discriminanterratum.pdf).
The erratum preserves these upper bounds, corrects the sharpness/lower-bound
claims, and includes the leading coefficient in the exceptional-prime factor.
No lower-bound version of Henriot's Theorem 6 is used here.

## A target-uniform short-window upper bound

**Upper-bound lemma FBS2.** For fixed \(\mathcal P,\theta\), there are constants
\(K,C>0\) such that, uniformly for \(X<N\leq2X\) and large \(X\),
\[
r_{\mathcal P,X^\theta}(N)
\leq K\frac{X^\theta}{(\log X)^\kappa}\mathcal S_C(N),
\qquad
\mathcal S_C(N)=\prod_{p\mid N}\left(1+\frac Cp\right). \tag{F4}
\]

Here is the full parameter check. Use
\[
Q_1(t)=t,\quad Q_2(t)=N-t,\quad Q(t)=Nt-t^2,
\qquad \mathcal F(u,v)=\xi(u)\xi(v).
\]
Both linear polynomials and their product are primitive; the product has degree
two, discriminant \(N^2\), coefficient norm at most \(N+1\), and leading
coefficient of absolute value one. Hence the erratum's leading-coefficient
correction adds no exceptional primes. Complete multiplicativity and boundedness
give \(\mathcal F\in\mathcal M_2(1,1,\varepsilon)\) for every
\(\varepsilon>0\), and the auxiliary minimal multiplicative function in the
corollary is \(\mathcal F\) itself.

The allowed \(t\) are contained in the interval \((x,x+y]\) with
\[
x=(N-X^\theta)/2-1,\qquad y=X^\theta+2.
\]
This interval stays in \((0,N)\) for all sufficiently large \(X\). In Henriot's
notation choose short-interval exponent \(\theta/2\), coefficient-size exponent
\(1/2\), and a sufficiently small positive growth parameter. Then
\(x\asymp X\), \(x^{\theta/2}<y\leq x\), and
\(x\geq C_0\|Q\|^{1/2}\), uniformly for the targets. These conditions hold
for every fixed positive \(\theta<1\), however small.

For \(p>2\), the number of roots of \(Q\) modulo \(p\) is two when
\(p\nmid N\) and one when \(p\mid N\). Each linear factor has one root.
Corollary 2 bounds the sum by a constant times
\[
y\Delta(N)
\prod_{2<p\leq x}\left(1-\frac{\rho_Q(p)}p\right)
\prod_{\substack{p\leq x\\p\nmid N}}
\left(1+\frac{\xi(p)}p\right)^2,
\qquad
\Delta(N)\leq\prod_{p\mid N}(1+1/p)^{C_0}. \tag{F5}
\]
The latter estimate is the remark after Henriot's Corollary 1; here bounded
weights and degree two give a fixed exponent. Importantly, the first product
starts at \(p>2\). A fixed divisor at 2 when \(N\) is odd therefore creates no
spurious zero. The finite factor at 2 in the other product is bounded and is
absorbed in the constant. The general upper corollary assumes primitivity; the
extra no-fixed-prime-divisor hypotheses from other statements are not needed.

Away from the prime divisors of \(N\), the baseline local factor is
\((1-2/p)(1+\xi(p)/p)^2\). At a prime divisor of \(N\), comparison with this
factor costs at most
\[
\frac{1-1/p}{(1-2/p)(1+\xi(p)/p)^2}
\leq\frac{p-1}{p-2}\leq1+\frac3p\qquad(p>2).
\]
These factors and \(\Delta(N)\) are bounded by \(\mathcal S_C(N)\) for a
sufficiently large fixed \(C\). Finally, (F3), including its special case
\(\xi\equiv1\), gives
\[
\log\prod_{2<p\leq x}(1-2/p)(1+\xi(p)/p)^2
=-\kappa\log\log x+O_{\mathcal P}(1).
\]
Since \(y\asymp X^\theta\) and \(\log x\asymp\log X\), (F5) proves (F4).
This application is uniform in the varying discriminant \(N^2\), which is why
a fixed-polynomial bound would not suffice.

## Block counting and the two moments

Partition \((X/2,X]\) into \(J\ll X/H\) half-open blocks of length at most
\(H=X^\theta\). Let \(m_j\) count the members of \(C(\mathcal P)\) in block
\(j\). Subtracting the two long-interval means in (F2) gives
\[
\sum_jm_j\gg_{\mathcal P}\frac{X}{(\log X)^{1-\alpha}}.
\]
Every ordered pair from one block has difference at most \(H\) and sum in
\((X,2X]\). Cauchy--Schwarz therefore gives
\[
M_1:=\sum_{X<N\leq2X}r_{\mathcal P,H}(N)
\geq\sum_jm_j^2
\geq\frac{(\sum_jm_j)^2}{J}
\gg_{\mathcal P,\theta}\frac{XH}{(\log X)^\kappa}. \tag{F6}
\]
No estimate for the population of an individual short block was assumed.

**Singular-factor lemma FBS3.** For every fixed \(C>0\) and real \(Y\geq1\),
\[
\sum_{N\leq Y}\mathcal S_C(N)^2\ll_C Y. \tag{F7}
\]
Indeed write \(\mathcal S_C(N)^2=\sum_{d\mid N}b(d)\), where \(b\) is
supported on squarefree integers and
\(b(p)=2C/p+C^2/p^2\). The coefficients are nonnegative, so
\[
\sum_{N\leq Y}\mathcal S_C(N)^2
=\sum_{d\leq Y}b(d)\lfloor Y/d\rfloor
\leq Y\prod_p\left(1+\frac{2C}{p^2}+\frac{C^2}{p^3}\right).
\]
The Euler product converges. This proves (F7) and, by (F4),
\[
M_2:=\sum_{X<N\leq2X}r_{\mathcal P,H}(N)^2
\ll_{\mathcal P,\theta}\frac{XH^2}{(\log X)^{2\kappa}}. \tag{F8}
\]

Put \(Q=H/(\log X)^\kappa\). Choose fixed \(A,B>0\) with
\(M_1\geq AXQ\), \(M_2\leq BXQ^2\) for large \(X\). At threshold
\(T=AQ/4\), the targets with fewer than \(T\) pairs contribute at most
\(2XT=AXQ/2\). The remaining set \(E\) therefore satisfies
\[
\#E\geq\frac{(\sum_{N\in E}r_{\mathcal P,H}(N))^2}{M_2}
\geq\frac{A^2}{4B}X.
\]
This proves FBS1 with \(c_1=A/4\), \(c_2=A^2/(4B)\). Moreover,
\[
\frac{X^\theta}{(\log X)^\kappa}
\geq2^{-\theta}\frac{N^\theta}{(\log N)^\kappa}
\quad(X<N\leq2X),
\]
which gives the stated target-sized formulation. Neither an exact averaged
asymptotic in the thin window nor a short-interval lower bound was used.

## Examples and the role of norm forms

In this table \(H=X^\theta\) with fixed \(0<\theta<1\); the lower bound
applies to a positive proportion of targets, with constants depending on the
fixed prime restriction and window exponent.

| Allowed prime factors | Prime density | Pairs on successful targets |
|---|---|---|
| \(p\equiv1\pmod3\) | \(1/2\) | \(\gg H/\log X\) |
| \(p\equiv1\pmod4\) | \(1/2\) | \(\gg H/\log X\) |
| \(p\equiv1\pmod {12}\) | \(1/4\) | \(\gg H/(\log X)^{3/2}\) |
| \(p\equiv1\pmod q\), fixed \(q\geq3\) | \(1/\varphi(q)\) | \(\gg H/(\log X)^{2-2/\varphi(q)}\) |
| Primes splitting completely in a fixed Galois field of degree \(d\) | \(1/d\) | \(\gg H/(\log X)^{2-2/d}\) |
| Primes splitting in each of \(r\) independent quadratic fields | \(2^{-r}\) | \(\gg H/(\log X)^{2-2^{1-r}}\) |

The first two rows are separate choices. Their union has density \(3/4\),
whereas their simultaneous intersection has density \(1/4\). Independence in
the last row means that the compositum has degree \(2^r\); \(r\) stays fixed.

For \(C_{12}\), successful targets are automatically \(2\pmod {12}\), so
FBS1 gives a positive proportion within that progression as well. In terms of
\(\delta=1-\theta\), it applies for every fixed \(0<\delta<1\). The
summands lie within \(O(N^{1-\delta})\) of \(N/2\).

Prime splitting alone need not characterize representations by a particular
quadratic form: ideal classes and primitivity require separate treatment. In
particular, \(x^2+36y^2\) does not define \(C_{12}\). The number 13 belongs
to \(C_{12}\) but cannot equal that form: nonzero \(y\) contributes at least
36, while 13 is not a square. Conversely \(85=7^2+36\) is represented, even
primitively, but \(85=5\cdot17\notin C_{12}\).

## A known density-one range in the congruence family

The exceptional-set programme must specify the window exponent. For every fixed
\(q\geq3\), write \(C_q=C(\{p:p\equiv1\pmod q\})\) and
\(\ell_q=\operatorname{lcm}(2,q)\). There is already an almost-all result in a
range of power windows, using prime summands alone.

**Known-range corollary Hq.** Fix \(q\geq3\) and \(2/3<\theta<1\). For
every fixed \(A>0\), all but \(O_{q,A,\theta}(X/(\log X)^A)\) targets
\(X<N\leq2X\), \(N\equiv2\pmod {\ell_q}\), have prime representations
\[
N=p_1+p_2,\qquad p_1\equiv p_2\equiv1\pmod {\ell_q},\qquad
|p_1-p_2|\leq X^\theta.
\]
Outside that exceptional set there are \(\gg_{q,\theta}X^\theta/(\log X)^2\)
such ordered prime pairs. These are also \(C_q+C_q\) pairs. We call the
\(q=12\) specialization H12.

This is a corollary of [Halupczok, *Goldbach's problem with primes in arithmetic
progressions and in short intervals*, JTNB 25 (2013), Theorem 3.1,
pp. 339--340](https://www.numdam.org/item/10.5802/jtnb.839.pdf), using the
singular series in her equations (1.1)--(1.2), p. 335. Here are the parameter details.
Put \(Y=X^\theta/4\), and cover \((X/2,X]\) by \(O(X/Y)\) blocks
\((t,t+Y]\), extending the last block if necessary. Apply Theorem 3.1 with
\(X_1=X_2=t\), \(Q_1=Q_2=\ell_q\), and \(a_1=2,a_2=1\). Choose its fixed
\(\varepsilon\) with \(0<\varepsilon<\theta-2/3\). Uniformly for
\(t\asymp X\), all its size conditions hold for large \(X\):
\[
t^{2/3+\varepsilon}\ll Y\leq t,\quad
Y\leq X_2\leq2X_1-Y,\quad
\ell_q\ll Y^{3/2}t^{-1/2},\quad
\ell_q\ll Yt^{-1/2}(\log Y)^{-B}.
\]
Restricting its nonnegative absolute-error sum to \(q_1=q_2=\ell_q\) gives
\[
\sum_{\substack{N/2\in(t,t+Y]\\N\equiv2\ (\ell_q)}}
\left|\sum_{\substack{p_1+p_2=N\\p_1\in(t,t+Y]\\p_1\equiv1\ (\ell_q)}}
\log p_1\log p_2-S(N,\ell_q,1)Y\right|
\ll_{q,A,\theta} Y^2(\log X)^{-A}.
\]
Her singular series satisfies
\(S(N,\ell_q,1)=\mathfrak S(\ell_qN)/\varphi(\ell_q)\geq c_q>0\)
on these targets. Indeed its baseline product is
\(2\prod_{p>2}(1-(p-1)^{-2})>0\), and the additional divisor factors are
at least one. Thus only \(O_{q,A,\theta}(Y/(\log X)^A)\) targets per block can
have weighted prime count below \(c_qY/2\). On every other target, both
\(p_1\) and \(N/2\) lie in the same block, so
\(|p_1-p_2|=2|p_1-N/2|\leq2Y<X^\theta\); also \(p_2\equiv1\pmod {\ell_q}\).
Dividing the weighted lower bound by \(O((\log X)^2)\) counts the pairs.
Summing the exceptions over \(O(X/Y)\) blocks proves Hq.

This is a logarithmic saving, not \(O(X^{1-\eta})\). Allowing arbitrary fixed
\(A\) does not turn it into a power saving. Its prime-pair lower bound also
has \((\log X)^2\), rather than the semigroup scale
\((\log X)^{2-2/\varphi(q)}\), which is \((\log X)^{3/2}\) at \(q=12\).
The smaller-window and power-saving questions require their own priority audits;
density one should not be presented as a new target throughout the whole range.

## An every-target ternary baseline

**Ternary corollary Tq.** Fix \(q\geq3\) and \(2/3<\theta<1\). Every
sufficiently large \(N\equiv3\pmod L\), where
\(L=\operatorname{lcm}(2,q)\), has
\[
N=p_1+p_2+p_3,\qquad p_i\equiv1\pmod L,\qquad
\max_i|p_i-N/3|\leq N^\theta. \tag{T}
\]
In particular, the three primes belong to \(C_q\). This is every-target
coverage in the stated progression, derived from established results.

We can use Halupczok's Theorem 3.1 twice. The first application supplies enough
prime candidates; the second finds a binary remainder outside its exceptional
set. This gives a complete proof from the same input already used for Hq.

**Candidate lemma Tq1.** For fixed \(q\geq3\), \(2/3<\theta<1\) and
\(c>0\), write \(Y=cT^\theta\). For every sufficiently large real \(T\),
\[
\#\{p\in(T,T+Y]:p\text{ prime},\ p\equiv1\pmod L\}
\gg_{q,\theta,c}\frac{Y}{(\log T)^2}. \tag{T1}
\]

To prove this, apply Theorem 3.1 with \(X_1=X_2=T\), length \(Y\),
\(Q_1=Q_2=L\) and \(a_1=2,a_2=1\). All size checks from Hq apply for this
fixed multiple of \(T^\theta\). Define
\[
W_{T,Y}(m)=
\sum_{\substack{u+v=m\\u,v\text{ prime}\\T<u\leq T+Y\\u\equiv1\ (L)}}
\log u\log v.
\]
There are \(2Y/L+O(1)\) targets \(m\equiv2\pmod L\) with
\(2T<m\leq2T+2Y\). Their main terms \(S(m,L,1)Y\) are uniformly at least
\(c_qY\). Taking logarithmic error exponent 4 in Theorem 3.1 and summing gives
\[
\sum_{\substack{2T<m\leq2T+2Y\\m\equiv2\ (L)}}W_{T,Y}(m)
\gg_{q,\theta,c}Y^2.
\]
For a fixed prime \(u\in(T,T+Y]\), at most \(2Y+1\) target values occur;
each determines \(v=m-u\), with both primes of size \(\asymp T\).
Its contribution is at most \(O(Y(\log T)^2)\). Dividing the last display by
this upper bound proves (T1). This weaker-than-PNT candidate count suffices.

Now fix an original target \(N\equiv3\pmod L\), and put
\(T=N/3\), \(Y=N^\theta/10\). Take \(N\) large enough that \(T\geq3Y\).
The candidate lemma supplies \(\gg_{q,\theta}Y/(\log N)^2\) primes
\(p_1\in(T,T+Y]\), \(p_1\equiv1\pmod L\). Discard the possible endpoint
\(p_1=T+Y\), which removes at most one candidate. Their distinct remainders
\[
m=N-p_1\equiv2\pmod L,\qquad 2T-Y<m<2T
\]
are covered by a second application of Theorem 3.1, now with
\[
X_1=T-Y/2,\qquad X_2=T,\qquad Q_1=Q_2=L,
\]
the same length \(Y\), and the same fixed residues. Its binary target range is
\(2T-Y<m\leq2T+Y\). The conditions \(Y\leq X_1\) and
\(Y\leq X_2\leq2X_1-Y\) follow from \(T\geq3Y\); all remaining size
conditions follow as in Hq from \(\theta>2/3\).

The absolute-error bound and \(S(m,L,1)\geq c_q>0\) leave only
\(O_{q,\theta}(Y/(\log N)^4)\) targets in this block with
\(W_{T,Y}(m)<c_qY/2\). Since the candidate-to-remainder map is injective,
these exceptions cannot exhaust \(\gg Y/(\log N)^2\) candidates. Choose one
with a prime representation \(m=p_2+p_3\). Then
\[
T<p_2\leq T+Y,\qquad T-2Y<p_3<T,
\qquad p_2\equiv p_3\equiv1\pmod L.
\]
Together with \(p_1\), their deviations from \(N/3\) are at most
\(2Y<N^\theta\), proving (T). No assertion about just almost all original
ternary targets was substituted for the every-target quantifier.

For comparison, [Matomäki--Maynard--Shao, arXiv:1610.02017v2,
Theorem 1.1](https://arxiv.org/html/1610.02017v2) proves an unrestricted
three-prime almost-equal theorem for every \(\theta>11/20\). That stated
theorem has no prescribed common arithmetic progression. Its stronger exponent
is not imported into Tq, and Tq is not presented as a new ternary Goldbach
breakthrough.

## Semigroup-scale density one and the variance main term

The prime-only result Hq is weaker than the following consequence of recent
short-interval Fourier estimates and an elementary cofactor sieve.

**Semigroup density-one theorem Dq.** Fix \(q\ge3\) and \(1/3<\theta<1\).
There is \(c_{q,\theta}>0\), independent of \(A\), such that for every fixed
\(A>0\), with \(H=X^\theta\),
\[
\#\left\{\begin{array}{c}
X<N\le2X,\quad N\equiv2\pmod L:\\
r_{q,H}(N)<c_{q,\theta}H/(\log X)^{\kappa_q}
\end{array}\right\}
\ll_{q,\theta,A}\frac{X}{(\log X)^A}.\tag{D}
\]
At \(q=12\), the lower bound has denominator \((\log X)^{3/2}\): a factor
\(\sqrt{\log X}\) larger than the prime-pair scale. The
[complete cofactor proof](PRESCRIBED_PRIME_SEMIGROUP_DENSITY_ONE.md) uses
Matomäki--Radziwiłł--Shao--Tao--Teräväinen's published 2026 Theorem 1.1(ii),
not merely an application of Hq. It counts numbers with a unique large prime
factor times a small allowed cofactor, averages Fourier errors before summing
cofactors, and proves a uniform lower bound for the modeled convolution.
The cofactor divisor sum restores the semigroup logarithmic gain even on
highly divisible targets. This is a derived consequence; originality remains
unaudited. The exception bound in Dq is logarithmic.

**Power-saving lower-count theorem D12Power.** For \(q=12\),
\(H=X^{99/100}\), there is a constant \(c>0\) such that
\[
 \#\{X<N\le2X:N\equiv2\pmod {12},\quad
 r_{12,H}(N)<cH/(\log X)^{3/2}\}
 \ll X^{1-10^{-6}}.
\]
The [complete proof](PRESCRIBED_PRIME_POWER_SAVING_EXCEPTIONS.md) uses the
integer-character convolution for \(w=\mathbf1_{C_{12}}\tau_4\),
polynomially uniform short-strip Fourier estimates, a lower sieve and small
cofactors. Sifting makes the retained weight bounded, so the conclusion counts
unweighted pairs. It gives a power-saving exceptional set for a lower bound;
it does not give a power-saving asymptotic for the full count or remove every
exception. Its priority assessment remains open.

The [variance proof and reduction](PRESCRIBED_PRIME_BINARY_VARIANCE.md) also
establish the explicit arithmetic model for \(C_q\). With
\(\alpha=1/\varphi(q)\) and
\(\mathcal B_q=\{p:p\nmid L,\ p\not\equiv1\pmod q\}\), put
\[
\mathfrak S_q(N)=L\mathbf1_{N\equiv2\ (L)}
\prod_{p\in\mathcal B_q}
\begin{cases}p/(p-1),&p\mid N,\\1-(p-1)^{-2},&p\nmid N.\end{cases}
\]
Define \(c_q>0\) by \(\sum_{n\le x}\mathbf1_{C_q}(n)
\sim c_qx(\log x)^{\alpha-1}\), and set
\[
\mathcal M_q(N,H)=c_q^2H(\log X)^{-\kappa_q}\mathfrak S_q(N)
\quad (X<N\le2X).
\]
The product has a uniform positive lower bound on compatible targets.
The character calculation is supplied in that proof; the Euler product is
not inferred just from the local forbidden residues. FBS2's upper-bound
factor \(\mathcal S_C\) remains a different object.

For every fixed \(0<\lambda\le1\), \(H=\lambda X\), theorem Vq proves
\[
\sum_{\substack{X<N\le2X\\N\equiv2\ (L)}}
 |r_{q,H}(N)-\mathcal M_q(N,H)|^2
=o\!\left(\frac{XH^2}{(\log X)^{2\kappa_q}}\right).\tag{V}
\]
This gives a full semigroup-scale asymptotic outside \(o(X)\) targets in
fixed-width windows, including the balanced case. No quantitative rate or
fixed-power shrinking window follows from this application of the fixed-body
Loughran--Matthiesen theorem.

For every fixed \(0<\theta<1\), including \(\theta\le1/3\), theorem Aq
in that proof gives the exact first moment
\[
\sum_{X<N\le2X}r_{\mathcal P,X^\theta}(N)
=c_{\mathcal P}^2\frac{XX^\theta}{(\log X)^{2-2\alpha}}
\{1+O_{\mathcal P,\theta}((\log X)^{-\eta})\}
\]
for some \(\eta>0\) and every fixed Frobenian prime set of density
\(\alpha>0\). Its short-interval input is Mangerel's published Theorem 1.7.
For \(C_q\), fixed-prime projection theorem Pq and reduction Eq show that
the only remaining input for (V) at \(H=X^\theta\) is the sharp upper bound
\[
\sum_{X<N\le2X}r_{q,H}(N)^2
\le\{\mathcal K_q+o(1)\}
 c_q^4\frac{XH^2}{(\log X)^{2\kappa_q}},\qquad
\mathcal K_q=L\prod_{p\in\mathcal B_q}(1+(p-1)^{-3}).\tag{E}
\]
The matching lower bound and both mixed main-term moments are already proved.
The [power-window theorem VqPower](PRESCRIBED_PRIME_POWER_WINDOW_VARIANCE.md)
now proves (E) for every fixed \(5/6<\theta<1\), including
\(\theta=99/100\), and gives
\[
 \sum_{X<N\le2X}|r_{q,H}(N)-\mathcal M_q(N,H)|^2
 \ll_{q,\theta}\frac{XH^2}{(\log X)^{2\kappa_q+1/4}}.
\]
The new proof uses short arithmetic progression means and a Type I/II
factorization covering the entire semigroup. It supplies a full-count
asymptotic outside \(O(X/(\log X)^{1/8})\) targets with relative error
\(O((\log X)^{-1/16})\). Its constants are ineffective. The sharp energy
for \(0<\theta\le5/6\), a power-saving full-count variance and every-target
binary counts remain open. D12Power separately gives power-saving exceptions
for the lower count at \(q=12,\theta=99/100\). Dq still gives a lower bound in its wider range
\(\theta>1/3\), with stronger logarithmic exception estimates; these are
different conclusions. Any proposed power-rate full-count variance needs
a sufficiently refined arithmetic main term, beyond its leading logarithmic
density; the present leading model does not specify its secondary terms.

The independent [Eisenstein angular theorem](EISENSTEIN_SIFTED_ANGULAR_DISTRIBUTION.md)
settles the macroscopic balanced instance of FBS-G at \(q=3\), and also its
mixed dilation \(A+(3\cdot A)\) on targets \(4\pmod6\). It places the first
summand in any prescribed fixed interval of positive relative length for every
sufficiently large compatible target. Its fixed-weight proof supplies no
numerical onset for an arbitrary interval.

The stronger [Theorem EisPower](EISENSTEIN_POWER_SHRINKING_WINDOWS.md) gives,
for every fixed \(0<\delta<1/3\) and every \(n\ge N_\delta\),
\[
\begin{array}{ll}
n\equiv2\pmod6:&n=a+b,\quad a,b\in A,\quad |a-b|\le n^{1-\delta},\\
n\equiv4\pmod6:&n=a+3b,\quad a,b\in A,\quad |a-3b|\le n^{1-\delta}.
\end{array}
\]
The proof gives an explicit formula for \(N_\delta\), at least
\(n^{-\delta}R_o(n)/20\) coordinate witnesses per prescribed orientation,
and \(\gg_{\delta,\nu}n^{1-\delta-\nu}\) distinct value pairs for each
fixed \(\nu>0\). For example, difference at most \(n^{3/4}\) holds for
every admissible \(n\ge2^{2^{40}}\).

The uniform operator bound controls harmonic degree and prime level together.
The exact triangular weight has an absolutely convergent Legendre series with
coefficient sum \(O(h^{-1/2})\). Counting large-prime obstructions only in
the central annulus supplies the additional \(h\) in the main lattice tail.
A separate positive Jackson-kernel construction verifies the proposed
Jackson--Parseval route, which by itself gives the weaker range
\(\theta>6/7\). The stronger theorem gives every fixed \(\theta>2/3\);
the endpoint is open. Its original coordinate-to-value conversion is weaker
than the semigroup-scale lower bound. The explicit cutoffs leave an enormous gap above
the finite balanced certificate. The separate
[value-count theorem EisValueCounts](EISENSTEIN_POINTWISE_VALUE_COUNTS.md)
now gives at least \(\gg_{\theta,\nu}n^\theta/(\log n)^{\log4+\nu}\)
distinct ordered pairs on every sufficiently large admissible target, for
\(2/3<\theta<1\) and \(\nu>0\). It uses a corrected-Henriot moment upper
bound and Hölder's inequality. This improves the former \(n^{\theta-\nu}\)
value-pair bound. The further [cofactor theorem EisSemigroupScale](EISENSTEIN_POINTWISE_SEMIGROUP_SCALE.md)
now proves
\[
 r_H(n)\asymp_\theta\frac{H}{\log n}
 \prod_{\substack{p\mid n,\ p\ge5\\p\equiv2\ (3)}}\frac{p-1}{p-2}
 \qquad(H=n^\theta,\ 2/3<\theta<1)
\]
for every sufficiently large admissible target, in both families and each
mixed orientation. It also covers every fixed central width \(H=cn\),
including the unrestricted count. The small-prime cofactors restore the
missing logarithm in the rough-pair lower sieve, with exact compensation at
split primes dividing the target. This proves the corresponding range of
FBS-U at \(q=3\). The exact leading constant, numerical onset, complete
exception lists, external review and priority remain unresolved.

The [moving-window theorem EisLocationWindows](EISENSTEIN_MOVING_WINDOWS.md)
now gives the same order \(HJ(n)/\log n\) at a prescribed center
\(x\), whenever \(1\le H\le M=\min(x,n-x)\) and
\(H^3\ge n^{1+\varepsilon}M\) for fixed \(\varepsilon>0\). A designated
small summand can therefore lie in \([a n^\theta,b n^\theta]\), for
fixed \(0\le a<b\) and any \(\theta>1/2\), on every sufficiently
large admissible target. Every fixed interior macroscopic interval also
has the expected order of pairs, in either mixed orientation. The nearly
equal threshold remains \(\theta>2/3\).


The further [bounded-factor theorem EisBoundedFactors](EISENSTEIN_BOUNDED_PRIME_FACTORS.md)
gives \(\gg_\theta H/(\log n)^2\) pairs in the same every-target power
windows, with every prime factor of each \(A\)-part exceeding
\(n^{\eta_\theta}\) for a fixed \(\eta_\theta>0\). Thus both parts
have a bounded number of prime factors, counted with multiplicity. Its
composite-level estimate treats all valuations of the target; its elementary
lower sieve does not claim prime or two-prime-factor summands.

## Five complementary research campaigns

Retain the definitions of \(C_q,L,\kappa_q,r_{q,H}\) from FBS-U. A weaker
existence goal remains useful for tracking implications:

**Open family goal FBS-G.** For every fixed \(q\geq3\), every sufficiently
large \(N\equiv2\pmod {\ell_q}\) has
\[
N=a+b,\quad a,b\in C_q,\qquad N/3\leq a,b\leq2N/3.
\]
The threshold may depend on \(q\). The congruence is necessary because every
member of \(C_q\) is odd and \(1\pmod q\); sufficiency for the whole family
is conjectural. The \(q=3\) instance is proved. The \(q=12\) instance is
the open target S12 in the case study. FBS-U implies FBS-G, since its lower
bound tends to infinity and eventually \(N^\theta\le N/3\). FBS-G supplies
neither the shrinking window nor the count in FBS-U.

**Campaign A: pointwise quartic correlation.** A proved intermediate is
[AsymPointwise](ASYMMETRIC_GAUSSIAN_EISENSTEIN_SUMS.md): with
\(U=\tau^2\mathbf1_{C_{12}}\), \(V=\tau\mathbf1_{A_3}\),
\(\sum_{a+b=N}U(a)V(b)=N\mathfrak S(N)+O(N^{1-\eta})\) on every
sufficiently large \(N\equiv2\pmod6\). The positive main term retains
a decreasing target factor at primes \(7\bmod12\). It gives many
unweighted one-sided pairs, with no balance or explicit cutoff claimed.
[AsymSemigroupScale](ASYMMETRIC_SEMIGROUP_SCALE.md) now supplies matching
counts of order \(N J_E(N)/(\log N)^{5/4}\), through a separately proved
small-coefficient family and a positive cofactor construction.
The two-character second weight produces square-root progression moduli;
replacing it by a second modulus-12 weight does not fit that proof.

The [subgroup splitting theorem](ASYMMETRIC_SUBGROUP_SPLITTING.md) gives a
separate mixed quadratic baseline. With \(C_r\) denoting the prime-factor
semigroup for \(1\bmod r\), its R1 specialization is
\[
 \#\{a+b=N:a\in C_4,\ b\in C_3\}
 \asymp \frac{N}{\log N}
       \prod_{\substack{p\mid N\\p\equiv11\ (12)}}\frac{p-1}{p-2}
\]
for every sufficiently large \(N\equiv0,2,6,8\pmod{12}\). The finite
residue density is 3 on these four classes and zero elsewhere. This does
not force both summands into \(C_{12}\): even at \(N\equiv2\pmod{12}\),
\(25+49=74\) has \(25\in C_4\) and \(49\in C_3\), with neither in
\(C_{12}\). The extra prime exclusions are \(5\bmod12\) on the first
side and \(7\bmod12\) on the second. In the tau-weighted R2 local model, either event
has density \((2p-1)/(p^2-p+1)\) when \(p\nmid N\), and zero when
\(p\mid N\). Each family has sieve dimension \(1/2\), so together they
have dimension 1. Starting with \(C_{12}+C_3\) imposes the first
restriction before the final sieve and leaves one half-dimensional family.
This compares the two routes; it does not provide their missing joint
large-coefficient distribution estimate.
The [unweighted coupled sieve](COUPLED_GAUSSIAN_EISENSTEIN_SIEVE.md) has a
half-dimensional Haar model instead: each excluded prime has density
\(1/(p-1)\) away from target divisors. Its bare proposed main
\(g_N(d)R_0(N)\) fails at every positive power level, because power-size
cofactors change the logarithmic density by a nonvanishing factor. This
failure does not refute the weighted R2 model or two-sided coverage. The
same note proves an unconditional two-large-prime upper bound of order
\(\epsilon^{5/4}N J_{12}(N)/(\log N)^{3/2}\); a sufficient sifted lower
bound with a corrected distribution model remains open. Thus the weighted
R2 dimension is not intrinsic to the underlying mixed sumset.
The [semilinear criterion](ASYMMETRIC_SEMILINEAR_SIEVE.md) instead isolates
a sufficient signed sieve remainder with level \(N^{1/2+\sigma}\),
\(\sigma>0\). Its local model, parity reduction and small-kernel tail are
proved; cancellation in its large-coefficient region is open.
The [critical switching theorem](ASYMMETRIC_CRITICAL_SWITCHING.md) improves
the sufficient distribution level to \(N^{1/2}(\log N)^{-B}\): stopping
below \(\sqrt N\) leaves only two large forbidden factors, whose total
weight is \(O(\epsilon^{3/2}X_N\mathcal V_N(z))\) for
\(z=N^{1/2-\epsilon}\). This bound is unconditional and uniform in the
target; its critical-level distribution hypothesis remains open.
The [full auxiliary average](ASYMMETRIC_AUXILIARY_AVERAGE.md) now proves a
weighted absolute correlation estimate through \(c\le N^\xi\), with
\(\xi\) chosen before its smaller error saving. The exact kernel identity
gives an \(L^1\) sieve level \(D=N^{\xi-\delta}\) for every fixed
\(0<\delta<\xi\). This controls a larger part of the critical remainder;
the square-root level is still open.
The [sign-transfer analysis](ASYMMETRIC_SIGN_TRANSFER.md) keeps a further
obligation explicit: transferring the quadratic sign to the second weight
retains the signed target-dependent model term. Two nonprincipal factors
and a prefix bound \(O_\nu(X^{1/2+\nu})\) for every \(\nu>0\) do not remove
its additive main. The
unrestricted coefficient convolution reconstructs the desired modulus-12
second weight. Replacing the remaining prime-supported coefficient by any
fixed-periodic function has full-order mean-square error, even before the
Rosser selector becomes nontrivial. These are specific simplification
obstructions, not a cancellation estimate or a counterexample to binary
coverage. The [well-factorability analysis](ASYMMETRIC_WELL_FACTORABLE_SWITCH.md)
proves that pure well-factorable divisor lower weights cannot replace the
critical Rosser weight, even on infinitely many actual target supports.
It also switches the branch above \(\sqrt D\) exactly to upper Rosser
weights at level \(D/p\), without asserting those upper weights are well
factorable. Fouvry's large-residue theorem controls some boxes satisfying
its precise Siegel--Walfisz hypothesis; neither it nor the switch estimates
the full critical remainder. The main target remains S12U:
every sufficiently large \(N\equiv2\pmod{12}\) is a sum of two members of
\(C_{12}\), without a balance restriction. Seek a pointwise lower bound for
\[
R(N)=576\sum_{a+b=N}h(a)h(b),\qquad
h(n)=r_G(n)r_E(n)/24.
\]
The [pointwise sieve note](SIMULTANEOUS_NORM_POINTWISE_SIEVE.md) identifies
the natural scale \(N\mathcal T(N)\), proves the matching upper bound, and
compares its local product with \(\mathcal T(N)\). A positive local product
is not a lower bound for the integer count on its own. The separate
[weighted correlation theorem SimNormPower](SIMULTANEOUS_NORM_POWER_CORRELATIONS.md)
now proves \(R(N)\asymp N\mathcal T(N)\) outside
\(O(X^{1-7\cdot10^{-6}})\) exceptional targets, with a power-saving weighted
mean square. Its fixed-balance corollary treats the full moving interval.
Every-target H1 and the natural-scale lower bound remain open. Ideal-norm
convolution, spectral methods and the circle method remain candidate
approaches to eliminating the exceptional set. The
[divisor-truncation obstruction](SIMULTANEOUS_NORM_DIVISOR_OBSTRUCTIONS.md)
shows that every fixed-power truncation has full-order \(L^1\) and \(L^2\)
error, even inside the correct residue class; a truncated sum can also be
positive at a value where \(h\) vanishes. Signed cancellation beyond an
absolute divisor approximation is required for that route.

The [scaling lemma S12Scale](SIMULTANEOUS_GAUSSIAN_EISENSTEIN_SUMS.md#the-first-pointwise-target-and-the-scaling-reduction)
removes every allowed prime-power factor from the target and preserves existence
and fixed-ratio balance. Its reduced core need not grow with the target.
An eventual theorem only for large reduced cores therefore needs additional
coverage for bounded cores; a large target can be \(14g\) with \(g\in C_{12}\).
For thirds-balance, the seven recorded exceptional reduced cores require
separate treatment of their multiples. The reduction does not preserve a fixed
power window or the semigroup-scale representation count automatically.

**Campaign B: the pointwise local-to-global sieve.** Independently of a lower
bound, the corrected Henriot theorem now gives
\[
R(N)\ll N\mathcal T(N),\qquad
\sum_{p>Y}B_p(N)\ll N\mathcal T(N)/Y+N^{15/16}.
\]
Consequently a natural-scale lower bound and uniformly valid marginal
asymptotics for each fixed forbidden prime suffice to complete the sieve.
This criterion replaces growing-prime control with fixed-prime marginals
and the natural-scale lower bound. It does not establish the original additive
\(o(N)\) version of H2 on every target. SimNormPower proves actual
bad-prime marginals simultaneously for \(p\le X^{10^{-8}}\), with error
\(O(X^{-10^{-6}})\), outside its common power exceptional set. The fixed
balance error is \(O_I(X^{-10^{-6}/20})\). Both results include arbitrary
prime valuations. The fully sifted coordinate-weight ratio is
\(\Gamma(N)=\prod_p(1-\beta_p(N))\ge10-\pi^2\), with the corresponding
power error. The all-target versions of both global inputs remain open.
This campaign and
Campaign A have the highest priority; neither the upper bound nor the tail
proves every-target coverage by itself.

**Campaign C: quantitative Eisenstein geometry.** EisPower now proves every
fixed power window \(H=n^\theta\), \(\theta>2/3\), in both ordinary and
mixed cases, with explicit onsets. The simultaneous weight/level Hecke bound,
triangle expansion and central-annulus tail have passed renewed internal
adversarial checks. Further goals include the endpoint, an exponent below
\(2/3\), the \(1/2\) scale and an asymptotic with an exact leading constant.
EisSemigroupScale already gives the pointwise semigroup order, with matching
singular-factor upper bounds, in the established window range.
EisLocationWindows extends this to translated windows with localization
cost \(O(\sqrt{M/H})\), and reaches \(\theta>1/2\) for a designated
small summand. This improves endpoint geometry while retaining the central
threshold. The [joint coset theorem EisJointCosetDistribution](EISENSTEIN_BOUNDED_PRIME_FACTORS.md#all-target-valuations-and-joint-cofactor-levels)
now removes the generic residue-space size cost at all target valuations,
including prescribed orientations and prime-square cofactor tests. Its
combined cofactor and sieve level has error proportional to \(BL\),
with fixed logarithmic factors, and approaches the square-root level at
macroscopic width. It does not lower the central exponent or give
every-target simultaneous prime restrictions.
Any proposed polylogarithmic window must also account for arithmetic gaps.
The [individual-moment barrier EisMomentBarrier](EISENSTEIN_POWER_SHRINKING_WINDOWS.md#what-individual-moment-bounds-can-force)
shows that positivity and the present separate Legendre bounds cannot force
an exponent below \(2/3\): an abstract positive measure with an empty
central interval already satisfies those bounds. Further progress needs
additional information, such as joint cancellation across degrees. This
barrier is not an arithmetic counterexample. The stronger
[single-index operator relaxation](EISENSTEIN_SINGLE_INDEX_OPERATOR_BARRIER.md)
retains self-adjoint convolution by a nonnegative measure and the finite
prime-to-index conditional coset means while leaving the same angular gap.
Positivity means preservation of nonnegative functions; the operator has
negative eigenvalues. This continuous model is not an arithmetic shell or
a multiplicative Hecke family: using it as a prime kernel would make the
prime-square recurrence require a negative identity atom. Relations between indices and arithmetic
discreteness remain possible ways beyond this information barrier.
Enlarging a finite scan or reducing an already positive sieve margin would
not settle these questions.

The [central-gap theorem FBSGap](PRESCRIBED_PRIME_CENTRAL_GAPS.md)
quantifies the polylogarithmic obstruction. For each fixed \(q\), infinitely
many compatible targets have no pair with difference at most
\(c_q\log N(\log\log N)^{1-2/\varphi(q)}\). For \(q=12\) this rules
out every-target windows of order \(\log N\); for \(q=3\) it rules out
\(o(\log N)\) and a sufficiently small constant multiple of \(\log N\).
The same note's arithmetic-large-sieve bound shows that a density-one set
fails when \(H=o((\log X)^{\kappa_q})\), including
\(H=o((\log X)^{3/2})\) at \(q=12\). These obstructions leave every
fixed positive power window in FBS-U untouched.

**Campaign D: sharp four-form energy.** The initial \(q=12\),
\(\theta=99/100\) target is resolved by VqPower, which covers every fixed
\(q\ge3\) and \(5/6<\theta<1\). It gives the full target-dependent
main term and the sharp second-moment coefficient, with a logarithmic rate.
D12Power now gives a power-saving exceptional set for a lower count at
\(q=12,\theta=99/100\). The next tasks are to lower the full-variance
exponent to \(\theta\le5/6\), obtain a power-saving full-count variance,
extend the lower-count exception result and develop a separate pointwise
mechanism. A mean-square theorem does not eliminate the last exceptional
target. Campaign D is not a prerequisite for Campaigns A and B.

**Campaign E: proof and priority audit.** Continue independent checks of the
Eisenstein spectral inputs and the power-window proof. Retrieve and audit the
full main theorem of Iwaniec--Pomykała's 1993 *Sums and differences of quartic
norms* before making any priority claim for S12U. Its introduction explicitly
poses the relevant splitting-prime question, but the full theorem has not been
obtained in this audit. The reproduced auxiliary proposition has extra
hypotheses and cannot stand in for it. No every-target modulus-12 result is
claimed here.

For existence in a target-sized window, the auxiliary exceptional count
\[
E_{q,\theta}(X)=
\#\{X<N\le2X:N\equiv2\pmod{\ell_q},\quad
r_{q,N^\theta}(N)=0\}
\]
is weaker than the full-count asymptotic at the same exceptional-set rate.
Dq already proves semigroup-scale density one for every fixed \(q\) and
\(\theta>1/3\), with its window \(X^\theta\le N^\theta\).
Semigroup-scale density one for \(0<\theta\le1/3\) remains unresolved in
this programme; this is not a complete literature-frontier claim.

For arbitrary Frobenian sets, no blanket assertion that local conditions
suffice is made. Loughran--Matthiesen's Remark 1.11 explains global
obstructions in its broader linear-form framework. The every-target conjectures
here concern the stated clean congruence family.

| Stage | Exact scope |
|---|---|
| General Frobenian baselines | FBS1 gives positive-proportion coverage, and Aq an exact first moment, in every fixed power window |
| Semigroup-scale binary density one | Dq gives many unweighted pairs for every fixed \(q\), \(\theta>1/3\), outside \(O_A(X/(\log X)^A)\) targets |
| Fixed-width full-count variance | Vq proves the explicit asymptotic outside \(o(X)\) targets for \(H=\lambda X\) |
| Existing ternary baseline | Tq gives every-target prescribed-AP prime triples for \(\theta>2/3\) |
| Eisenstein power windows | EisPower gives every fixed \(\theta>2/3\), explicit onsets and many coordinate witnesses in both families; \(H=n^{3/4}\) holds above \(2^{2^{40}}\) |
| Every-target Eisenstein counts | EisSemigroupScale proves the order \(H J(n)/\log n\) for \(q=3,2/3<\theta<1\), also for fixed central widths and unrestricted pairs |
| Location-sensitive Eisenstein counts | EisLocationWindows proves the same order when \(H^3\ge n^{1+\varepsilon}\min(x,n-x)\) and \(H\le\min(x,n-x)\); a designated summand may have scale \(n^\theta\) for \(\theta>1/2\) |
| General asymmetric counts | Every fixed \(q\) and chosen quadratic \(\chi\) have the pointwise \(C_q+A_\chi\) order \(N J_\chi/(\log N)^{3/2-1/\varphi(q)}\) on \(N\equiv2\pmod L\) |
| Bounded excluded prime factors | Every-target counts at \(N J_q/(\log N)^{2-2/\varphi(q)}\), with the first summand in \(C_q\) and only boundedly many large excluded prime factors in the second |
| Asymmetric distribution level | A full auxiliary-modulus average proves an \(L^1\) level at every fixed exponent below \(\xi\); the critical square-root level remains open |
| Highest-priority arithmetic target | Prove S12U through the critical sieve remainder or the global norm lower bound and fixed-prime marginals |
| Power-window full-count variance | VqPower proves (E) and quantitative centered variance for every fixed \(q\), \(5/6<\theta<1\); lower exponents and power-saving full variance remain open |
| Power-saving lower-count exceptions | D12Power gives \(H/(\log X)^{3/2}\) pairs outside \(O(X^{1-10^{-6}})\) targets for \(q=12,H=X^{99/100}\) |
| Weighted simultaneous-norm correlations | SimNormPower gives exact raw and fully sifted main terms and actual bad-prime marginals outside a power exceptional set, including fixed balance intervals and all target valuations |
| Main general conjecture | FBS-U asks for the semigroup-scale count on every compatible large target |
| Effective completion | Bridge an effective pointwise cutoff to finite verification; the present Eisenstein onset leaves an enormous gap |

A larger first moment alone does not establish the remaining stages. The
[Iwaniec--Pomykała audit](SIMULTANEOUS_GAUSSIAN_EISENSTEIN_SUMS.md#literature-and-the-remaining-obligation)
remains necessary before claiming novelty for an every-target theorem.

## Verification scope

The proofs above and in the linked variance/cofactor notes use the precisely
identified published inputs. Their applications, character factors, arithmetic
tails and interval geometry received internal AI proof/source checks. The general-family analytic baselines have no effective numerical cutoff
supplied here; EisPower separately gives the explicit onset stated above. No
formal verification or external human review is claimed; the linked notes describe the exact scope of their
diagnostic checks. The
\(C_{12}\) computation remains finite evidence with its original scope in the
case study. Publication and originality are separate from the mathematical
classification of these corollaries.

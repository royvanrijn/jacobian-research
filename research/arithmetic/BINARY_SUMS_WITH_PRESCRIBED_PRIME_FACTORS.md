# Binary sums with prescribed prime factors

The central research target is a uniform binary lower bound at the scale of the
whole multiplicative semigroup, on every compatible large target. That is the
open conjecture FBS-U below. Positive-proportion coverage, semigroup-scale density one in a range of power
windows, fixed-width variance, and an every-target ternary result are proved
baselines. This note and its linked proofs derive these results from published
tools and keep the remaining binary estimates and their priority assessment
separate.

The [simultaneous Gaussian--Eisenstein case study](SIMULTANEOUS_GAUSSIAN_EISENSTEIN_SUMS.md)
keeps the exact norm identities, lifted local densities, conditional pointwise
criterion and finite certificates for the class \(1\pmod {12}\). The immediate
arithmetic target is every-target existence for \(C_{12}+C_{12}\), initially
without balance. In parallel, the
[Eisenstein power-window theorem](EISENSTEIN_POWER_SHRINKING_WINDOWS.md)
now supplies one explicit shrinking exponent for \(q=3\).

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
macroscopic balanced family goal FBS-G below. It remains a conjecture; none of
the baseline results in this note establishes (U).

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
unaudited. The exception bound is logarithmic, not a power saving.

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
The current upper bound (F8) has an unspecified constant and does not prove
(E). Dq proves a lower bound on almost all targets, while (E) would control
deviations of the full count about its target-dependent model.

The independent [Eisenstein angular theorem](EISENSTEIN_SIFTED_ANGULAR_DISTRIBUTION.md)
settles the macroscopic balanced instance of FBS-G at \(q=3\), and also its
mixed dilation \(A+(3\cdot A)\) on targets \(4\pmod6\). It places the first
summand in any prescribed fixed interval of positive relative length for every
sufficiently large compatible target. Its fixed-weight proof supplies no
numerical onset for an arbitrary interval.

The stronger [Theorem EisPower](EISENSTEIN_POWER_SHRINKING_WINDOWS.md) gives,
for every \(n\ge2^{2^{411}}\),
\[
\begin{array}{ll}
n\equiv2\pmod6:&n=a+b,\quad a,b\in A,\quad |a-b|\le n^{99/100},\\
n\equiv4\pmod6:&n=a+3b,\quad a,b\in A,\quad |a-3b|\le n^{99/100}.
\end{array}
\]
Its uniform operator bound controls the harmonic degree and prime level
together. A shrinking Bernstein minorant then supports a nonnegative weighted
union bound. The theorem also gives
\(\gg_\varepsilon n^{99/100-\varepsilon}\) distinct ordered value pairs.
This is weaker than the semigroup-scale count in FBS-U. The enormous explicit
onset does not close the gap above the existing finite balanced certificate;
the complete exception lists remain unresolved. External review and priority
are separate open obligations.

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

**Campaign A: pointwise quartic correlation.** The first target is S12U:
every sufficiently large \(N\equiv2\pmod{12}\) is a sum of two members of
\(C_{12}\), without a balance restriction. Seek a pointwise lower bound for
\[
R(N)=576\sum_{a+b=N}h(a)h(b),\qquad
h(n)=r_G(n)r_E(n)/24.
\]
The [pointwise sieve note](SIMULTANEOUS_NORM_POINTWISE_SIEVE.md) identifies
the natural scale \(N\mathcal T(N)\), proves the matching upper bound, and
compares its local product with \(\mathcal T(N)\). A positive local product
is not a lower bound for the integer count. Neither H1 nor this stronger
natural-scale lower bound is proved. Ideal-norm convolution, spectral methods
and the circle method remain candidate approaches, with the obstacles recorded
in the case study.

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
\(o(N)\) version of H2. Both global inputs remain open, including targets
with varying valuations at a fixed prime. This campaign and
Campaign A have the highest priority; neither the upper bound nor the tail
proves every-target coverage by itself.

**Campaign C: quantitative Eisenstein geometry.** EisPower completes the first
power-shrinking existence target, with exponent \(99/100\), in both ordinary
and mixed cases. Independent verification of its uniform-weight, uniform-level
Hecke argument is essential. Further goals are improved angular discrepancy,
stronger shrinking exponents and, separately, the semigroup-scale count.
The next work is analytic and does not require a large finite scan or a smaller
version of the already positive sieve margin.

**Campaign D: sharp four-form energy.** Fix \(q=12\), initially
\(\theta=99/100\), and prove (E) for \(H=X^\theta\). By Eq this is
equivalent to the full-count centered variance (V). The model, its uniform
positive lower bound, the first moment and the mixed moments are established;
the sharp second-moment upper bound on the thin region remains open. A proof
would give an asymptotic for the full semigroup count outside \(o(X)\) targets,
stronger than Dq's lower bound. Dq does not prove (E), even near \(\theta=1\).

A further quantitative milestone is a power-saving variance that yields
\(O(X^{1-\eta})\) exceptions. Neither the fixed-width Vq nor Dq's arbitrary
fixed logarithmic savings supplies this. A mean-square theorem alone also does
not remove the last exceptional target. Campaign D is not a prerequisite for
Campaigns A and B or for a different pointwise lower-bound mechanism.

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
| Eisenstein power windows | EisPower proves existence for \(H=n^{99/100}\) in both residue classes, above the explicit \(2^{2^{411}}\) onset |
| Highest-priority arithmetic target | Prove S12U by resolving the global lower bound and fixed-prime distribution inputs |
| Sharp variance target | Prove (E), initially at \(\theta=99/100\); seek a power-saving exceptional estimate afterward |
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

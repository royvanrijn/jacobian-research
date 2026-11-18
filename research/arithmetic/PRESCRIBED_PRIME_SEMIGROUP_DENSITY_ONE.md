# Semigroup-scale binary representations in power windows

Fix a modulus \(q\geq 3\). Let \(\mathcal P_q\) be the primes congruent to
\(1\bmod q\), and put
\[
C_q=\{n\geq 1:p\mid n\Longrightarrow p\in\mathcal P_q\},\qquad
L=\operatorname{lcm}(2,q),\qquad
\alpha=\frac1{\varphi(q)},\qquad \kappa=2-2\alpha.
\]
Every element of \(C_q\), including \(1\), is \(1\bmod L\). Let
\(r_{q,H}(N)\) count ordered pairs of positive integers \(a,b\in C_q\)
with \(a+b=N\) and \(|a-b|\leq H\).

**Density-one theorem \(D_q\).** Fix \(1/3<\theta<1\) and put \(H=X^\theta\).
There is a constant \(c_{q,\theta}>0\) such that for every \(A>0\),
\[
\#\left\{X<N\leq 2X:
  N\equiv2\pmod L,\quad
  r_{q,H}(N)<c_{q,\theta}\frac{H}{(\log X)^\kappa}
\right\}
\ll_{q,\theta,A}\frac{X}{(\log X)^A}
\quad (X\to\infty).                                             \tag{D}
\]
The lower-bound constant is independent of \(A\). In particular, at \(q=12\)
the denominator is \((\log X)^{3/2}\). Replacing the width by \(N^\theta\)
in the lower-bound conclusion on this dyadic interval costs at most a fixed
factor, since \(X^\theta\leq N^\theta\leq 2^\theta X^\theta\).

The proof counts a subset in which each summand has an allowed prime factor
larger than its square root. Its smaller cofactor supplies the two powers of
\((\log X)^\alpha\) gained over the scale of pairs of primes. A Fourier
estimate for primes in almost all short intervals transfers this subset count
to a model that removes primes below an explicit cutoff. We give a lower
bound for the model on every compatible target, including targets divisible by
many allowed primes.

## A minorant and its short-interval model

Choose \(\delta>0\) so small that
\[
\delta<\min(1/4,\theta/4),\qquad
\frac{\theta-\delta}{1-\delta}>\frac13,
\]
and set \(D=X^\delta\). Consider a summand interval \(I_x=(x,x+H]\)
with \(X/3\leq x\leq 2X\). Cover these starts by a fixed number of intervals
\([V,2V]\) with \(V\asymp X\), taking \(V=X/4,X/2,X\) and assigning a dyadic
endpoint to just one interval. For \(m\leq D\) put
\[
Y_m=V/m,\qquad h_m=H/m.
\]
There is a fixed \(\varepsilon>0\) for which
\(Y_m^{1/3+\varepsilon}\leq h_m\leq Y_m^{1-\varepsilon}\) for all the
indices and all sufficiently large \(X\). For the lower bound the ratio
\(H/m\) to \((X/m)^{1/3+\varepsilon}\) is smallest at \(m=D\). For the
upper bound the corresponding ratio to \((X/m)^{1-\varepsilon}\) is largest
at \(m=1\). The strict parameter inequalities absorb the fixed factors
\(V/X\).

Write
\[
P(r)=\prod_{p<r}p,\qquad R_Y=\exp((\log Y)^{1/10}),\qquad
\Lambda^\sharp_Y(k)=
 \frac{P(R_Y)}{\varphi(P(R_Y))}\mathbf1_{(k,P(R_Y))=1}.
\]
For \(n\in I_x\) define
\[
\begin{aligned}
g_x(n)&=\frac1{\log(3X)}
  \sum_{\substack{mk=n,\ m\leq D,\ m\in C_q\\ k\equiv1\ (L)}}
           \mathbf1_{\{k\ {\rm prime}\}}\log k,\\
g_x^\sharp(n)&=\frac1{\log(3X)}
  \sum_{\substack{mk=n,\ m\leq D,\ m\in C_q\\ k\equiv1\ (L)}}
           \Lambda^\sharp_{Y_m}(k).
\end{aligned}                                                     \tag{1}
\]
Both functions vanish off \(I_x\). Any prime \(k\) in the first sum is
larger than \(\sqrt n\) for large \(X\), since \(m\leq X^\delta<\sqrt n\).
Such a prime is unique, \(n\in C_q\), and \(k\leq n\leq 3X\). Consequently
\(0\leq g_x(n)\leq\mathbf1_{C_q}(n)\). Mertens' product estimate gives
\(0\leq g_x^\sharp(n)\ll\tau(n)\), where \(\tau\) counts positive divisors.

For finitely supported sequences write
\(\widehat f(\beta)=\sum_n f(n)e^{2\pi i\beta n}\) and
\((f*g)(N)=\sum_n f(n)g(N-n)\). We claim that for every \(K>0\),
\[
\int_{X/3}^{2X}\sup_{\beta\in\mathbb R}
 \left|\widehat g_x(\beta)-\widehat g_x^\sharp(\beta)\right|\,dx
 \ll_{q,\theta,K} XH(\log X)^{-K}.                              \tag{2}
\]

The input is the short-interval estimate for \(\Lambda-\Lambda^\sharp_Y\)
of Matomäki, Radziwiłł, Shao, Tao and Teräväinen. For
\(Y^{1/3+\varepsilon}\leq h\leq Y^{1-\varepsilon}\) and any \(J>0\), the
supremum of the absolute exponential sum with a linear phase is at most
\(h(\log Y)^{-J}\) for \(y\in[Y,2Y]\) outside a set of measure
\(O_{\varepsilon,J}(Y(\log Y)^{-J})\). Their maximal norm permits restriction
to any arithmetic progression in \((y,y+h]\). These statements follow by
taking the fixed unit circle and its exponential character as the nilmanifold
and Lipschitz function in their theorem. [*Higher uniformity of arithmetic
functions in short intervals II. Almost all intervals*, equations (1.1),
(1.4) and Theorem 1.1(ii)](https://doi.org/10.1007/s00222-026-01408-6)

Apply the estimate for each \(m\) at \(Y=Y_m,\ h=h_m,\ y=x/m\), with
\(k\equiv1\bmod L\) and frequency \(m\beta\). The model scale \(Y_m\)
is constant over the dyadic range of \(x\). On the exceptional starts the
trivial absolute bound is \(O(h_m\log X)\), since \(h_m\gg1\), and both
weights are \(O(\log X)\). Changing variables \(dx=m\,dy\) bounds the
mean over \(x\in[V,2V]\) of the supremum for this \(m\) by
\(O_J(h_m(\log X)^{1-J})\). Higher prime powers in \(\Lambda\) have
mean absolute contribution
\(O(h_mY_m^{-1/2}(\log X)^2)\). Indeed there are
\(O(\sqrt{Y_m}\log X)\) powers \(p^j\leq 3Y_m\) with \(j\geq2\), each
of weight at most \(\log(3X)\), and the set of \(y\) for which such a power
belongs to \((y,y+h_m]\) has length at most \(h_m\). Thus the same
arbitrary logarithmic saving holds after replacing \(\Lambda(k)\) by
\(\mathbf1_{\{k\ {\rm prime}\}}\log k\). Sum the integrated errors by the
triangle inequality, use \(\sum_{m\leq D}h_m\ll H\log X\), and divide by
\(\log(3X)\). Taking \(J\) sufficiently large proves (2).

## A uniform lower bound for the model

We prove that for each compatible \(X<N\leq2X\) and each considered
start that satisfies \(|N-(2x+H)|\leq H/2\),
\[
(g_x^\sharp*g_x^\sharp)(N)\gg_{q,\theta}
                   H(\log X)^{-2+2\alpha}.                       \tag{3}
\]
The interval \(I_x\cap(N-I_x)\) has length
\(h_N=H-|N-(2x+H)|\geq H/2\). Put \(z=\exp((\log(3X))^{1/10})\).
Every \(R_{Y_m}\leq z\), and \(\log R_{Y_m}\asymp_\theta\log z\).
Whenever \(k\) is coprime to \(P(z)\), its model weight in (1) is
\(\gg_\theta(\log z)/\log X\).

In the product of the two sums in (1), retain cofactor choices
\(m=du,\ m'=dv\), where \(d\) is a squarefree divisor of \(N\) in \(C_q\),
\(u,v\in C_q\), and
\[
d,u,v\leq D^{1/2},\qquad (u,v)=1,\qquad(uv,N)=1.
\]
Each such pair of cofactors determines \(d=(m,m')\), then \(u,v\). The
quotients \(k,l\) of the summands now satisfy
\[
uk+vl=N',\qquad N'=N/d\equiv2\pmod L,\qquad k\equiv l\equiv1\pmod L.
\]
Because \(u,v,d\equiv1\bmod L\) and \((u,v)=1\), these solutions form a
progression \(k=k_0+vLt,\ l=l_0-uLt\) in \(t\in\mathbb Z\). In fact
\(uk\equiv N'-v\bmod vL\) has a unique solution class modulo \(vL\), and
its reduction modulo \(L\) is \(1\). Restricting \(du k\) to
\(I_x\cap(N-I_x)\) gives an interval of length
\[
J_0=\frac{h_N}{duvL}\gg_q X^{\theta-3\delta/2}.                 \tag{4}
\]

For \(p<z\) let \(\nu_p\) be the number of roots of \(k(t)l(t)\bmod p\).
At \(p\mid L\) there are no roots. For \(p\nmid Luv\) there are two,
coincident precisely if \(p\mid N'\). For \(p\mid uv\), one form is a
nonzero constant modulo \(p\), since \((uv,N)=1\), and the other has one
root. Thus \(\nu_p\leq2\), all primes with \(\nu_p>0\) are odd, and
\(\nu_p=1\) at every \(p\mid N'\) outside \(L\).

Uniformly in these parameters, at least
\(\tfrac12J_0\prod_{p<z}(1-\nu_p/p)\) values of \(t\) in (4) give
\((k(t)l(t),P(z))=1\), for all sufficiently large \(X\).
To see this, truncate inclusion–exclusion at the least odd integer
\(b\geq C\log\log X\), so the truncation is a lower bound. For each
squarefree \(j\) in the truncation the count of simultaneous roots is
\[
J_0\prod_{p\mid j}\frac{\nu_p}{p}
       +O\left(\prod_{p\mid j}\nu_p\right)
\]
by the Chinese remainder theorem. There are at most \(z^b\) such
integers, each contributing an error at most \(2^b\). Their total is
\(X^{o(1)}\). The absolute omitted tail of the product expansion is at most
\(\sum_{i>b}(\sum_{p<z}\nu_p/p)^i/i!\).
Since \(\sum_{p<z}\nu_p/p\leq 2\log\log z+O(1)\), taking \(C\)
sufficiently large makes this \(O((\log X)^{-10})\).
The full product is \(\gg(\log z)^{-2}\) by Mertens' estimate and
\(\prod_{p\geq3}(1-(p-1)^{-2})>0\). Together with (4), these two error
bounds give the asserted count.

After multiplication by the two model weights, the factor
\((\log z)^2\) supplies, up to a constant,
\(\prod_{p<z}(1-1/p)^{-2}\). At \(p\geq3\) with two roots the normalized
factor is \(1-(p-1)^{-2}\). Replacing two roots by one increases it by
\((p-1)/(p-2)\geq1+1/p\). We may therefore use the following lower bound
for the contribution of each retained cofactor pair:
\[
\gg_{q,\theta}
 \frac{H}{duv(\log X)^2}
 \prod_{\substack{p\mid N',\ p<z\\p\in\mathcal P_q}}
        \left(1+\frac1p\right).                                  \tag{5}
\]

## Summing the cofactors

Set \(y=X^{\delta/C_0}\), with \(C_0\) a sufficiently large absolute
constant. Start by summing over coprime \(u,v\) composed of primes
\(p\leq y\) in \(\mathcal P_q\) that do not divide \(N\). Without their
size restrictions their harmonic mass is
\[
\sum_{u,v}\frac1{uv}=
\prod_{\substack{p\leq y,\ p\in\mathcal P_q\\p\nmid N}}
     (1-1/p)^{-2}(1-1/p^2)
\gg_{q,\theta}(\log X)^{2\alpha}
   \prod_{\substack{p\mid N,\ p\in\mathcal P_q\\p\leq y}}(1-1/p)^2.
                                                                    \tag{6}
\]
Here the fixed-progression Mertens estimate is
\(\sum_{p\leq y,\,p\in\mathcal P_q}p^{-1}=\alpha\log\log y+O_q(1)\).
[Loughran–Matthiesen, Lemma 2.4(2)](https://arxiv.org/pdf/1904.12845v2)

Under the probability measure obtained by normalizing this harmonic
mass, either specified coordinate is divisible by \(p\) with probability
\(1/(p+1)\), and its exponent, conditional on being positive, has mean
\(p/(p-1)\). The expected logarithm of each coordinate is at most
\(\sum_{p\leq y}(\log p)/(p-1)\leq 2\log y+O(1)\) by Mertens'
estimate. Markov's inequality and a union bound show that \(u,v\leq
D^{1/2}\) retain at least half the mass for \(C_0\) large. Thus (6) is
also a lower bound up to a constant for the permitted \(u,v\).

Let \(T_N=\{p\in\mathcal P_q:p\mid N,\ p\leq\log X\}\). The sum of
\(1/p\) over the prime divisors \(p>\log X\) of \(N\) is \(O(1)\),
since their number is at most \(\log(2X)/\log2\). Thus the product on the
right of (6) is bounded below by a fixed positive multiple of
\(\prod_{p\in T_N}(1-1/p)^2\). Both \(y\) and \(z\) exceed \(\log X\)
for large \(X\).

Restrict \(d\) to squarefree products of primes in \(T_N\). Every
\(p\in T_N\) not dividing \(d\) still divides \(N'\). Hence the product in
(5) is bounded below by
\[
\prod_{p\in T_N}(1+1/p)\prod_{p\mid d}(1+1/p)^{-1}.
\]
Before the size cutoff, the sum of \(d^{-1}\prod_{p\mid d}(1+1/p)^{-1}\)
equals \(\prod_{p\in T_N}(1+1/(p+1))\). Under this normalized measure,
the mean of \(\log d\) is \(\sum_{p\in T_N}(\log p)/(p+2)=O(\log\log X)\).
At least half the weight therefore satisfies \(d\leq D^{1/2}\) for
large \(X\), uniformly in \(N\). Combining the three factors at each
\(p\in T_N\) leaves
\[
(1-1/p)^2(1+1/p)(1+1/(p+1))=1-3/p^2+2/p^3 .
\]
Their product over any \(T_N\) is bounded below by a positive constant,
since these factors are positive and differ from one by \(O(p^{-2})\).
Summing (5) now proves (3).

## From almost all interval starts to integer targets

Choose \(B>0\). By (2) with \(K=2B\), the set \(E\) of starts
\(x\in[X/3,2X]\) at which the Fourier supremum exceeds
\(H(\log X)^{-B}\) has measure \(O_{q,\theta,B}(X(\log X)^{-B})\).
To sample this continuous estimate at block starts, set \(s=H/2\).
The identity
\[
\int_0^s\#\{j\in\mathbb Z:t+js\in E\}\,dt=|E|
\]
gives a grid of step \(s\) with \(O_{q,\theta,B}((X/H)(\log X)^{-B})\)
exceptional starts. For each grid start \(x\) take the target block
\([2x+H-H/2,\,2x+H+H/2)\). These blocks tile the line. All blocks meeting
\((X,2X]\) have starts in \([X/3,2X]\) for large \(X\). The exceptional
blocks contain \(O_{q,\theta,B}(X(\log X)^{-B})\) targets in total.

On the remaining blocks, Parseval's identity bounds the sum over integers
\(N\) of \(|(g_x*g_x)(N)-(g_x^\sharp*g_x^\sharp)(N)|^2\) by
\[
H^2(\log X)^{-2B}\sum_{n\in I_x}|g_x(n)+g_x^\sharp(n)|^2.
\]
The summand intervals on the grid have bounded overlap. Also
\(\sum_{n\leq 3X}\tau(n)^2\ll X(\log X)^3\): prime-power comparison gives
\(\tau(n)^2\leq\tau_4(n)\), and summing the fourfold divisor count gives
the stated estimate. Consequently the sum of squared differences on the
assigned targets of all these blocks is
\(O_{q,\theta,B}(XH^2(\log X)^{3-2B})\).

For each assigned \(N\), all pairs counted by \(g_x*g_x\) satisfy
\(|a-b|\leq H\), and their weights are bounded by
\(\mathbf1_{C_q}(a)\mathbf1_{C_q}(b)\). Choose \(c_{q,\theta}\) to be
half the constant in (3). At a compatible target violating (D)'s lower
bound on a remaining block the modeled and actual convolutions therefore
differ by at least \(c_{q,\theta}H(\log X)^{-\kappa}\). The last
mean-square estimate bounds their number by
\(O_{q,\theta,B}(X(\log X)^{3+2\kappa-2B})\).
Taking \(B\) sufficiently large in terms of \(A\) proves (D).

The separate [power-window variance proof VqPower](PRESCRIBED_PRIME_POWER_WINDOW_VARIANCE.md)
now gives a full-count asymptotic for \(5/6<\theta<1\). The cofactor
argument here supplies a lower bound over the wider range \(\theta>1/3\),
with stronger logarithmic exception estimates, and does not itself determine
the full count. The separate [D12Power proof](PRESCRIBED_PRIME_POWER_SAVING_EXCEPTIONS.md)
strengthens the exceptional count to \(O(X^{1-10^{-6}})\) for the lower
bound at \(q=12,H=X^{99/100}\), using a bounded weighted rough kernel.
Variance for \(0<\theta\le5/6\) and a lower bound for
every sufficiently large compatible target remain open. Historical priority
for (D) has not been established. The fixed-shape variance and the broader
binary conjecture are recorded in the companion notes.
[Binary variance and its arithmetic main term](PRESCRIBED_PRIME_BINARY_VARIANCE.md)
and [binary sums with prescribed prime factors](BINARY_SUMS_WITH_PRESCRIBED_PRIME_FACTORS.md).

# Prime clusters and prime-counting subadditivity

The second Hardy–Littlewood conjecture asks whether
\(\pi(x+y)\leq\pi(x)+\pi(y)\) for all integers \(x,y\geq2\). The problem
remains open here. Two restricted results concern the width-3159 route.
The Engelsma tuple has no prime translate below the explicit bound in
Section 1. Separately, a sparse deletion of actual primes removes every
potential violation over a range of widths starting at 3159 while preserving
the natural-scale asymptotics of lower-order fixed-shift correlations,
whenever those asymptotics hold for the primes.

## 1. The width-3159 target

Section 2 of [Hensley and Richards, *Primes in intervals*](https://www.impan.pl/en/publishing-house/journals-and-series/acta-arithmetica/all/25/4/100146/primes-in-intervals)
proves that for every \(\varepsilon>0\) and all sufficiently large widths \(y\),
there is an admissible set in \([1,y]\) with at least
\(\pi(y)+(\log 2-\varepsilon)y/(\log y)^2\) elements. Admissibility says
that the offsets omit a residue modulo each prime. The prime-tuple conjecture
would realize each fixed pattern infinitely often. The recent subadditivity
result of [Chahal, Elma, Fellini, Vatwani and Vo, Theorem 1.1](https://arxiv.org/html/2503.02766v1)
uses an error bound in the prime number theorem and applies when
the smaller summand exceeds an error-dependent lower threshold.
[Maynard's cluster theorem, Theorem 3.2](https://arxiv.org/pdf/1405.2593v2)
gives a lower bound of order \(\log y\) for the number of primes in some
width-\(y\) intervals. The required threshold at general large width is
\(\pi(y)+1\sim y/\log y\).

The [Engelsma input](admissible_447_3158.txt) is the [offset list in the MIT narrow-tuple tables](https://math.mit.edu/~primegaps/tuples/admissible_447_3158.txt).
Its 447 offsets occupy \([0,3158]\), so a prime translate \(n+H\) with
\(n\ge3\) would violate subadditivity at \((n-1,3159)\), where
\(\pi(3159)=446\). For a general search at that width the exact statistic is

\[
\sum_{2\le n\le X}\binom{\pi(n+3159)-\pi(n)}{447}.
\]

It is positive precisely when there is a violating width-3159 interval
with start \(2\le n\le X\). It counts all 447-element prime subsets of
those intervals, including patterns other than the Engelsma offsets.

### Exact bound for the specified offsets

The [integer checker](verify_hardy_littlewood_447.py) counts the primes
up to 3159 and computes the permitted start residues modulo every prime
\(p\le447\). Larger primes cannot witness inadmissibility because there
are only 447 offsets. The checker handles starts \(n\le443\) directly,
since a prime divisor of a translated entry could then equal the entry.
For \(n>443\) every tested divisor excludes that start. The 36 primes with
one permitted residue give, by the Chinese remainder theorem,

\[
\begin{aligned}
n&=a+Mq,\qquad q\ge0,\\
a&=226869188698870760778384337119124029399468630472000601300261,\\
M&=264223627364521452551144665505940260576702745525641127244470.
\end{aligned}
\]

The next eight constraints are the first eight nonsingleton primes in
increasing order. They give 54,000 residue classes for \(q\) modulo their
product \(Q=753007512523486741\). These are exactly the possible quotients
in \([0,Q)\) after those eight tests. Every one fails a remaining residue
test. Thus the specified offsets have no prime translate with positive
start strictly below \(a+MQ > 1.98\cdot10^{77}\).

Internal replay used Python 3.12.3 with no third-party dependencies.
From the repository root run
`python3 research/arithmetic/verify_hardy_littlewood_447.py`.
An internal cross-check using an Eratosthenes sieve and the closed-form
CRT sum reproduced the exclusion.
The input SHA-256 is
`ff0d1ee8249430d651a7031b41dc2daf81ae2bc1c7c6b84e2186aa51b3dfa0e6`.

## 2. Deletion and fixed-shift correlations

For a set of primes \(S\), a finite set of distinct integer shifts \(T\),
and a weight \(w_X\) on \([1,X]\cap\mathbb Z\) with \(|w_X|\le1\), put
\[
C_{S,T,w}(X)=\sum_{1\le n\le X}w_X(n)\prod_{h\in T}1_S(n+h).
\]

**Theorem.** Let \(k\ge2\) and \(Y_0\ge k\) be fixed integers and
\(\eta\ge0\) a fixed real number. Write
\(W(t)=\max\{Y_0,\lfloor(\log t)^\eta\rfloor\}\) for \(t\ge2\), using
\(W(t)=Y_0\) when \(\eta=0\). There is a subset \(P^-\) of the primes
with fewer than \(k\) elements in every integer interval \((x,x+y]\)
with \(x\ge2\) and \(1\le y\le W(x)\). For each fixed nonempty \(T\),
uniformly in the weights above,
\[
\left|C_{\mathbb P,T,w}(X)-C_{P^-,T,w}(X)\right|
\ll_{k,Y_0,\eta,T} \frac{X W(X)^{k-1}}{(\log X)^k}.
\]
In particular, for \(r=|T|<k\) and \(\eta<(k-r)/(k-1)\) the difference
is \(o(X/(\log X)^r)\).

**Proof.** Delete all members of each \(k\)-element set of primes
\(\{p_1,\ldots,p_k\}\) with \(p_1<\cdots<p_k\) and
\(p_k-p_1<W(p_1)\). If the remaining primes occupied an interval of
width \(y\le W(x)\) with at least \(k\) points, their least and greatest
among a chosen \(k\) would differ by at most \(y-1<W(p_1)\). This uses
\(p_1>x\) and the monotonicity of \(W\).

To count deletions up to \(X\), each responsible prime tuple has
\(p_1\le X\) and offset set \(H\subset[0,W(X)-1]\) of size \(k\) with
\(\min H=0\). Set \(z=X^{1/100}\). All tuples with \(p_1\le z\) can
account for at most \(z+W(X)\) deleted integers together. For the other
tuples, the anchor \(p_1\) avoids \(\nu_p(H)=\#(H\bmod p)\) residues
for every prime \(p\le z\). Apply [Tao and Teräväinen, Lemma 3.2](https://doi.org/10.1112/jlms.12663)
with this value for \(p\le z\) and zero for \(z<p\le X\). The bound on
the number of such anchors for a fixed \(H\) is
\[
\ll_k X\prod_{p\le z}(1-\nu_p(H)/p)
\ll_k \frac{X}{(\log X)^k}\,\mathfrak S(H),\qquad
\mathfrak S(H)=\prod_p\frac{1-\nu_p(H)/p}{(1-1/p)^k}.
\]
The second estimate uses Mertens' product bound and the fact that for
sufficiently large \(X\), all primes above \(z\) exceed \(W(X)\) and
\(2k\). Their normalized factors are \(1+O_k(p^{-2})\), with a product
bounded away from zero uniformly in \(H\). If a local factor vanishes it
does so below \(z\) and the sifted anchor count is zero.

The fixed-\(k\) Gallagher estimate in [Ford, equation (1)](https://arxiv.org/pdf/1108.3861v2)
gives \(\sum_{J\subset[1,2W],\,|J|=k}\mathfrak S(J)=O_k(W^k)\).
Each set \(H\subset[0,W-1]\) with minimum zero contributes \(W\)
distinct translates to that sum, with its singular series unchanged.
The minimum identifies the translate uniquely, and the factors are
nonnegative. Thus the sum over minimum-zero sets is \(O_k(W^{k-1})\).
Multiplying the anchor count by \(k\) and summing bounds the number of
deletions up to \(X\) by \(O_{k,Y_0,\eta}(X W(X)^{k-1}/(\log X)^k)\),
including the smaller term \(z+W(X)\). For a change in a summand of
\(C_{S,T,w}(X)\), some deleted prime must occur at a specified shift.
Each such prime can change at most \(|T|\) summands, and these primes
are at most \(X+\max(0,\max T)\). This gives the stated bound. \(\square\)

Taking \(k=447\), \(Y_0=3159\) and \(\eta=0\) shows that even fixed-shift
asymptotics through order 446, at errors \(o(X/(\log X)^r)\) at each
order \(r\le446\), would also hold for a prime subset with at most 446
primes in every width-3159 interval. Setting \(\eta>0\) gives the same
obstruction over all widths \(3159\le y\le W(x)\), since
\(\pi(y)\ge446\). For example, all fixed two-point asymptotics survive
for any \(\eta<445/446\). Weights allow a fixed progression or a bounded
restriction on the anchors. The argument leaves open whether the actual
primes have a positive 447th factorial count at width 3159, or a
violation at any other width.

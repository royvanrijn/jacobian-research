# Central gaps for binary sums with prescribed prime factors

The power-window conjecture allows every fixed positive power of the target.
A polylogarithmic version needs a separate qualification: congruences can
exclude every permitted pair in a central interval. The elementary construction
below uses the Chinese remainder theorem and the prime-progression Mertens
estimate. It concerns existence, so coordinate weights cannot remove the
obstruction.

For fixed \(q\ge3\), write
\[
 C_q=\{n\ge1:p\mid n\Longrightarrow p\equiv1\pmod q\},\qquad
 L=\operatorname{lcm}(2,q),\qquad
 \kappa_q=2-\frac2{\varphi(q)}\ge1.
\]
Let \(r_{q,H}(N)\) count ordered positive pairs \(a,b\in C_q\) with
\(a+b=N\) and \(|a-b|\le H\).

**Theorem FBSGap.** For every fixed \(q\ge3\) there is \(c_q>0\) and
infinitely many integers \(N\equiv2\pmod L\) for which
\[
 r_{q,H_q(N)}(N)=0,\qquad
 H_q(N)=c_q\log N\,(\log\log N)^{1-2/\varphi(q)}.
\]
No explicit value of \(c_q\) or starting point is asserted.

For \(q=3,4,6\), this rules out every-target windows \(H=o(\log N)\),
and windows \(c\log N\) for sufficiently small \(c>0\). For
\(\varphi(q)>2\), it rules out every-target windows \(H=O(\log N)\).
In particular, the modulus-12 obstruction has order
\(\log N\sqrt{\log\log N}\). The CRT theorem alone does not exclude every fixed
window \((\log N)^{1+\varepsilon}\), or any positive power \(N^\theta\).

## Prime estimates

Let
\[
 \mathcal B_q=\{p:p\nmid L,\ p\not\equiv1\pmod q\},\qquad
 \beta=1-1/\varphi(q),\quad 2\beta=\kappa_q.
\]
These forbidden primes have density \(\beta\). The fixed-progression
prime number theorem and Mertens estimate give
\[
 \pi_{\mathcal B_q}(y)\sim\beta\frac y{\log y},\qquad
 \sum_{\substack{p\le y\\p\in\mathcal B_q}}\log p\sim\beta y,
 \qquad
 \sum_{\substack{p\le y\\p\in\mathcal B_q}}\frac1p
 =\beta\log\log y+O_q(1).
\]
These are [Loughran--Matthiesen, Lemma 2.4(1)--(3)](https://arxiv.org/html/1904.12845v2#S2)
applied to the Frobenian indicator of \(\mathcal B_q\). Consequently
\[
 V(y)=\prod_{\substack{p\le y\\p\in\mathcal B_q}}(1-2/p)
 \asymp_q(\log y)^{-\kappa_q}.
\]
All these primes are odd because \(2\mid L\), so every factor is positive.

## Pre-sieving symmetric offsets

Fix a large integer \(K\), put \(y=\sqrt K\), and choose residues
\(M\pmod p\) independently and uniformly for the forbidden primes
\(p\le y\). Later impose \(M\equiv1\pmod L\); these conditions are
compatible because \(p\nmid L\).

For an offset \(1\le k\le K\), the pair \(M-Lk,M+Lk\) survives a
prime \(p\) if neither term vanishes modulo \(p\). There are two distinct
forbidden residues for \(M\) unless \(p\mid k\), in which case there
is one. Thus its survival probability is
\[
 V(y)F_y(k),\qquad
 F_y(k)=\prod_{\substack{p\le y\\p\in\mathcal B_q\\p\mid k}}
       \frac{p-1}{p-2}.
\]
The positive squarefree-divisor expansion gives
\[
 \sum_{k\le K}F_y(k)
 =\sum_{d\mid Q_0}\left(\prod_{p\mid d}\frac1{p-2}\right)
    \left\lfloor\frac Kd\right\rfloor
 \le K\prod_{p\mid Q_0}\left(1+\frac1{p(p-2)}\right)\ll K,
\]
where \(Q_0=\prod_{p\le y,p\in\mathcal B_q}p\). The infinite product
in the last bound converges since its logarithm is bounded by a constant
multiple of \(\sum p^{-2}\).

The expected number of surviving offsets is therefore
\(\ll_q K/(\log K)^{\kappa_q}\). Choose one residue system for which
there are at most that many survivors. This is a finite averaging argument:
a residue system can be found, in principle, by enumerating the finitely many
possibilities. No randomness remains in the resulting construction.

## Blocking the remaining pairs

For each surviving positive offset choose a different fresh forbidden prime
\(p>y\), and impose \(M\equiv-Lk\pmod p\). Choose one further fresh
forbidden prime and impose \(M\equiv0\pmod p\), blocking the central
pair \(M+M\).

All fresh primes can be chosen at most \(C_qK\), for a fixed sufficiently
large constant \(C_q\). Indeed their number is
\(O_q(K/(\log K)^{\kappa_q})+1\), whereas the forbidden primes in
\((\sqrt K,C_qK]\) number \(\sim\beta C_qK/\log K\).
Use \(\kappa_q\ge1\) and enlarge \(C_q\) to dominate the constant
when \(\kappa_q=1\).

Let \(Q\) be the product of all small and fresh primes used. The Chinese
remainder theorem combines every imposed condition with \(M\equiv1\pmod L\).
Choose the positive representative \(LQ\le M<2LQ\), and set \(N=2M\).
For large \(K\), \(Q\ge Q_0\) grows exponentially in \(\sqrt K\),
so \(M>LK\) and all the central terms considered are positive.

Every member of \(C_q\) is \(1\pmod L\). A possible pair in the
window \(|a-b|\le2LK\) must consequently be
\(a=M+Lk,b=M-Lk\) with \(|k|\le K\). For \(k\ne0\), either a
small-prime congruence already excluded one term or its assigned fresh prime
excludes one. The offset zero is blocked by its separate fresh prime.
Hence
\[
 r_{q,2LK}(N)=0,\qquad N\equiv2\pmod L.
\]

## Size of the target

The small-prime cost is
\(\log Q_0=O_q(\sqrt K)\). Each fresh prime costs at most
\(\log(C_qK)\), so
\[
 \log N\ll_q\sqrt K+
  \left(\frac K{(\log K)^{\kappa_q}}+1\right)\log(C_qK)
 \ll_q\frac K{(\log K)^{\kappa_q-1}}.
\]
Also \(\log N\ge\log Q_0\gg_q\sqrt K\), ensuring that the constructed
targets tend to infinity. These bounds imply
\(\log\log N\asymp_q\log K\) and
\[
 2LK\gg_q\log N\,(\log\log N)^{\kappa_q-1}.
\]
Choose \(c_q\) below the resulting positive constant. Then
\(H_q(N)\le2LK\) for infinitely many constructed targets, proving the
theorem.

## A density obstruction for smaller polylogarithmic windows

A separate upper-bound sieve gives a stronger conclusion for most targets.
Put \(\beta=1-1/\varphi(q)\), so \(\kappa_q=2\beta\).

**Proposition FBSSmallWindows.** Uniformly for \(1\le H\le X\),
\[
 \#\{X<N\le2X:N\equiv2\pmod L,\ r_{q,H}(N)>0\}
 \ll_q X\left\{(\log X)^{-\beta}
               +\frac H{(\log X)^{\kappa_q}}\right\}.
\]
Thus a density-one set of compatible targets has no representation whenever
\(H=o((\log X)^{\kappa_q})\). For a sufficiently small fixed
\(c_q'>0\), a positive proportion have no representation even at
\(H=c_q'(\log X)^{\kappa_q}\). In particular, the density-one failure
range for \(q=12\) extends to \(H=o((\log X)^{3/2})\).
This is an upper-bound consequence of the classical large sieve, not a new
claim about the sharp threshold.

For the proof write \(N=2M\), \(a=M-k\), \(b=M+k\), where
\(0\le k\le H/2\), after choosing the ordering \(a\le b\). Summing
over these offsets gives an upper bound for successful targets. The case \(k=0\) contributes
\(O_q(X(\log X)^{-\beta})\) by the one-variable Frobenian mean theorem,
[Loughran--Matthiesen, Lemma 2.8](https://arxiv.org/html/1904.12845v2#S2).

For \(k\ge1\), at every forbidden odd prime a successful center must avoid
\(\Omega_p=\{k,-k\}\pmod p\), of size
\(\nu_p=2-\mathbf1_{p\mid k}\). The arithmetic large sieve gives
\[
 \#\{1\le M\le X:M\bmod p\notin\Omega_p\ (p\le Q)\}
 \le\frac{X+Q^2}{\mathcal H_k(Q)},\qquad
 \mathcal H_k(Q)=\sum_{\substack{d\le Q\\d\ \mathrm{squarefree}}}
 \prod_{p\mid d}\frac{\nu_p}{p-\nu_p}.
\]
Use \(\nu_p=0\) at every other prime. The precise imported statement is
[Kowalski, *Amplification arguments for large sieve inequalities*,
arXiv:0911.2020v2, Section 1, (1)--(2), with
\(\Delta=N-1+Q^2\)](https://arxiv.org/html/0911.2020v2#S1).
Since neither summand is divisible by \(p\), the center \(M\) avoids
both residues \(\pm k\pmod p\). No exception for a summand equal to a
sieving prime is needed.

Set \(Q=\sqrt X\), \(y=Q^{1/10}\), and
\[
 B_k(y)=\prod_{\substack{p\le y\\p\in\mathcal B_q}}
       (1-\nu_p/p)^{-1}.
\]
The complete squarefree divisor sum on these primes equals \(B_k(y)\).
Normalize its positive terms to a probability distribution. A prime is
included in its random divisor with probability \(\nu_p/p\), independently
of all the others. Consequently
\[
 \mathbf E\log d=\sum_{\substack{p\le y\\p\in\mathcal B_q}}
                 \frac{\nu_p\log p}{p}
 \le2\log y+O(1)<\tfrac12\log Q
\]
for sufficiently large \(X\), uniformly in \(k\). The estimate
\(\sum_{p\le y}\log p/p=\log y+O(1)\) is the ordinary prime-Mertens
estimate, or Lemma 2.4(3) applied to the constant prime function and partial
summation. Markov's inequality now gives
\(\mathcal H_k(Q)\ge B_k(y)/2\).

The reciprocal product is exactly
\[
 B_k(y)^{-1}=V(y)F_y(k)
 \ll_q(\log X)^{-\kappa_q}F_y(k),
\]
with the same functions as in the CRT argument. The positive divisor expansion
therefore gives
\[
 \sum_{1\le k\le H/2}
 \#\{M\le X:M-k,M+k\in C_q\}
 \ll_q\frac X{(\log X)^{\kappa_q}}
        \sum_{1\le k\le H/2}F_y(k)
 \ll_q\frac{XH}{(\log X)^{\kappa_q}}.
\]
The sieve bound includes all centers in \([1,X]\); restricting to positive
summands and \(X<2M\le2X\) only decreases it. Adding the diagonal case
proves the proposition. Compatible targets number \(X/L+O(1)\), giving
both density conclusions.

## Scope

These are consequences of a classical congruence construction, fixed-modulus
prime estimates and the arithmetic large sieve; no originality claim is made. It introduces no obstruction to
FBS-U's fixed power windows. For the Eisenstein result it shows that an
arbitrarily small positive power cannot simply be replaced by an arbitrarily
small polylogarithmic window. It does not settle the endpoint \(\theta=2/3\)
of EisPower or any new positive-power range.

The congruence averaging, fresh-prime supply, positivity, admissibility,
size inversion and large-sieve denominator estimate received separate internal
derivation checks. No large numerical
search, external review or formal verification is asserted.

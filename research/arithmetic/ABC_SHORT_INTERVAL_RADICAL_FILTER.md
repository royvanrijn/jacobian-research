# Short intervals with small radicals and real-filtered target tests

An interval of sufficiently short positive numerical width contains
few integers of small radical. The proof splits off the bounded number
of individually large primes and uses the existing real strip estimate
only when all the remaining prime powers have small total numerical
weight. This allows a real prefilter in the modular test for popular
target radicals, lowering the required exponent of each independent
prime block by nearly one.

## An elementary bound in relative intervals

Fix \(A>0\), \(0<\tau<1\) and \(\eta>0\). For all sufficiently
small fixed \(L>0\), uniformly for \(K>0\) and \(0\le\ell\le A\), the
number of integers \(n\le X\) satisfying
\(\operatorname{rad}(n)\le X^{L\ell}\) and
\(\lvert n/K-1\rvert\le X^{-\tau}\) is
\[
    O_{A,\tau,\eta,L,\varepsilon}
          \left(X^{(1-\tau/2+\eta)L\ell+\varepsilon}\right).
                                                               \tag{SI1}
\]
As usual constant multiples in the precision can be absorbed
by a smaller fixed allowance in \(\tau\).

Choose \(0<\xi<\min(\eta,\tau)\). The part with
\(n\le X^{\tau-\xi}\) contains at most one candidate for large
\(X\), since the relative interval has absolute length \(O(X^{-\xi})\)
on this part when nonempty. For each remaining \(n\) order
its exact prime powers by decreasing exponent. Take the prefix
\(P\) just before their product first reaches \(n X^{-\tau+\xi}\)
and call the next power \(p^e\). Put \(R=n/(Pp^e)\). Then
\[
 R\le X^{\tau-\xi},\qquad
 \log\operatorname{rad}P
 \le \frac{\log P}{\log n}\log\operatorname{rad}n
 \le (1-\tau+\xi)L\ell\log X .
\]
The first radical inequality uses that the prefix's mean exponent,
weighted by logarithms of its primes, is at least the overall
weighted mean. The second uses
\(\log P<\log n-(\tau-\xi)\log X\) and \(n\le X\).

Fix \(P,e\) and bin \(p,\operatorname{rad}R\) and
\(\operatorname{rad}P\) dyadically. For each \(p\) there is at
most one \(R\) in the strip for large \(X\), since two such
integers would differ by \(O(X^{-\xi})\). For each \(R\) there
is at most one \(p\), since two such positive bases would
differ by \(O(X^{L\ell-\tau})\) using \(e\ge1\). Fix \(L\) small
enough that \(LA<\tau\). If the three binned radical exponents
in base \(X\) are \(b_P,b_p,b_R\), their sum is at most
\(L\ell+o_X(1)\). This uses that the three factors have disjoint
prime supports. Enumerate \(P\) at cost \(X^{b_P+\varepsilon'}\)
by its radical and then the usual uniform fixed-support estimate.
For the other pair, enumerate the smaller of the two possible
dictionaries, at cost \(X^{\min(b_p,b_R)+\varepsilon'}\).
The total exponent is at most
\[
 b_P+\min(b_p,b_R)\le \frac{L\ell+b_P}{2}+o_X(1)
                     \le (1-\tau/2+\xi/2)L\ell+o_X(1).
\]
The logarithmic number of bins and choices of \(e\) have inner
subpower cost. This proves (SI1).
[Uniform enumeration at fixed support](ABC_EXCEPTIONAL_SET_AND_POWERFUL_TRIPLES.md#a-radical-count-contradiction-to-the-claimed-smooth-number-lower-bound)

## Arbitrarily small normalized counts in shorter intervals

Fix \(A>0\) and \(\eta>0\). There exist fixed
\(\theta=\theta(A,\eta)>0\) and \(L_0=L_0(A,\eta)>0\)
such that at every fixed \(0<L<L_0\), uniformly for all
\(K>0\),
\[
 \#\left\{n\le X:\operatorname{rad}(n)\le X^{AL},\
                    \left|\frac nK-1\right|\le X^{-1+\theta}\right\}
                  \ll_{A,\eta,L,\varepsilon}
                          X^{\eta L+\varepsilon}.
                                                               \tag{SI2}
\]
It follows that every interval of length \(X^\theta\) in a
fixed dyadic shell at height \(X\) has the same bound.
One may shrink \(\theta\) when translating to an interval
or absorbing a constant in its length. More generally, if
the centre is at least \(X^{h_0}\) for a fixed \(h_0>0\), the
same conclusion holds for absolute width \(X^{\theta'}\) with
sufficiently small \(\theta'>0\) depending also on \(h_0\).

For the proof put \(b=\min(\eta/16,A/4)\). In each \(n\) call
a prime large when \(p\ge X^{bL}\), and write
\[
                  n=V\prod_{i=1}^m p_i^{e_i},
          \qquad  m\le M=\lfloor A/b\rfloor .
\]
All powers on smaller primes go into \(V\). Each large
prime satisfies \(p_i\le X^{AL}\) and each \(e_i\le1/(bL)\).
Thus the number of choices of \(m\) and its ordered exponents
is bounded at the fixed budget. Order the primes, so this
presentation has exactly one tuple for each complete value.

Choose \(s>0\) with \((M^2+M)s<\eta/4\). Choose a fixed
\(E>0\) sufficiently small for the (RUW1) bound at every
dimension \(1\le m\le M\), with \(C_i=A+1\) and
real precision \(\tau=1/2\). Also take \(E<1/4\).
For \(V\le X^E\), the real strip count costs at most
\(X^{\eta L/4+\varepsilon}\) by (RUW1). It applies for
\(\theta<1/2\). If \(m=0\) in this branch, the relative
interval contains at most one integer because \(V\le X^E\).
Use the lower precision \(h_0/2\) in place of \(1/2\),
and take \(E<h_0/4\), for centres at least \(X^{h_0}\).
[Unequal-height real strip count (RUW1)](ABC_SHORT_HEIGHT_AND_REAL_GAP.md#a-real-strip-with-unequal-root-heights)

For \(V>X^E\), partition its exact prime powers greedily into
groups of radical at most \(X^{2bL}\), using that each prime
has \(\log_X p<bL\). Each completed group other than the
last has radical logarithm greater than \(bL\). Thus the
number of groups is at most \(M+1\). Some group has product
\(F\ge X^{E/(M+1)}\). Choose \(\theta<E/(2(M+1))\).
The bound \(n\le X\) means that every nonempty relative
interval under consideration has absolute width
\(O(X^\theta)\). For a given \(F\) the interval then contains
at most one multiple of \(F\) for large \(X\).
There are at most \(X^{2bL+\varepsilon}\) possible \(F\) by
enumerating \(\operatorname{rad}F\) and then its powers up
to height \(X\). This is below \(X^{\eta L+\varepsilon}\).
For centres at least \(X^{h_0}\) and absolute width
\(X^{\theta'}\), the same argument applies with
\(\theta'\le\min\{h_0/2,E/(2(M+1))\}\).
All choices of constants precede the small fixed budget
and its height limit. No weak convergence of valuation
measures is needed.

In particular for every \(\eta>0\) and all sufficiently large
fixed \(k\), uniformly in intervals of length at most \(X^\theta\)
at height \(X\), the number of \(k\)-full integers is
\(O_{k,\eta,\varepsilon}(X^{\eta/k+\varepsilon})\).
Take \(L=1/k\) and \(A=2\) in a shell with upper endpoint
a fixed multiple of \(X\). The small width exponent is
independent of \(k\) once \(\eta\) has been chosen.
Chan's Theorems 1–3 give uniform bounds in the interval
length, including a separate smooth-number hypothesis for
his unconditional power bound. Trifonov's short-interval
estimates concern squarefull numbers. The result here uses
the small radical budget and sufficiently large \(k\)
to obtain an arbitrarily small coefficient of \(1/k\)
in a short power interval.
[Chan, Theorems 1–3](https://arxiv.org/pdf/2207.08874v1)
· [Trifonov, squarefull numbers in short intervals](https://arxiv.org/abs/2609.17715)

## Coarse complete values and the row count

In a primitive height box suppose the three complete letters
have numerical heights at least \(h_0>0\) and raw radical
caps \(X^{AL}\). Choose the width \(X^\theta\) as in (SI2)
for centres at least \(X^{h_0}/2\). Each interval bin of
one designated letter has at most \(X^{\eta L+\varepsilon}\)
complete values. In balanced individual caps, using (FL1)
for each complete value gives per-bin exponent at most
\(1/6+\eta\) in base \(X^L\), up to arbitrary small allowances.
Thus a near-full balanced family must use at least
\(X^{(1/6-\eta-o(1))L}\) bins on each such letter.

There is also a simultaneous information consequence. On
a uniform family of rows, let \(J_T=\lfloor T/X^\theta\rfloor\).
For every letter \(H(T\mid J_T)\le(\eta L+\varepsilon)\log X+O(1)\).
Since two complete values recover the row, the pair of
corresponding bins loses entropy at most
\(2(\eta L+\varepsilon)\log X+O(1)\). The three bins
satisfy \(J_A+J_B=J_C\) or \(J_A+J_B+1=J_C\).
Keeping one of the two carries gives an exact additive
graph on the bin labels, with at least
\(|\mathcal F|X^{-2\eta L-2\varepsilon}/O(1)\) distinct
edges. Each bin pair lifts to at most the product of the
two single-bin dictionaries. The bin labels themselves
are subject to no asserted radical bound.
[Fixed complete-letter fibre count (FL1)](ABC_FIXED_LETTER_RADICAL_CAPS.md#a-pointwise-count-with-prescribed-prime-pools)

## Consequences for real-filtered popularity tests

Use \(N=X^L\), with \(L\) sufficiently small and fixed,
and original representatives \(A,B,C\le X\) with
\(X/2<C\le X\). Let the radicals of both free letters be
at most \(N^a\), where \(a>0\) is fixed. The test on each
candidate support record includes the real filter
\[
                          |A+B-C|\le X^\theta .
\]
For every \(\eta>0\) one can take a fixed positive \(\theta\)
so that for each target radical \(R_C\), the total number
of eligible free support pairs is \(O(N^{a+\eta+o_X(1)})\).
First fix the exact complete \(C\), at subpower cost on
that target radical. For each pair choose its smaller
summand, breaking ties by role. These values number
\(N^{a+o_X(1)}\) by radical enumeration. Once this summand
is fixed, the larger is in an interval of width
\(2X^\theta\) at height comparable to \(X\). Apply (SI2)
with cap \(a\) and an appropriately smaller width exponent.
The same upper bound holds for triples of representatives
and hence for distinct support pairs.

Apply the modular rejection argument (MP1) only to these
eligible candidates. If a support record is false, its
candidate discrepancies are nonzero integers of magnitude
at most \(X^\theta\). It still has probability at most
\(N^{-\beta+o_X(1)}\) of passing an independent prime block
of total exponent \(\beta\), even if the representatives
change at different test primes. With target radicals
at most \(N\), threshold \(N^u\), and \(j\) independent
blocks, the false contribution to the expected number of
popular target radicals is therefore at most
\[
           O\left(N\min\{1,N^{a+\eta-\beta-u+o_X(1)}\}^{j}\right).
                                                               \tag{SI3}
\]
As in (MP1), the error is in comparison with the true target
tail at half the threshold. For \(\beta>a+\eta-u\), a
fixed sufficiently large \(j\) makes this error smaller by
a fixed power than \(N^{1-u}\). With balanced free caps the
sufficient block exponent can thus be taken arbitrarily close
to \(1-u\) from above, in place of \(2-u\).
In particular a single prime per block can be just above
exponent one for testing occurrence of target radicals and
just above exponent one half for the square-root
multiplicity threshold.
[Independent-block count (MP1)](ABC_MODULAR_TARGET_POPULARITY.md#support-multiplicities-and-independent-blocks)

For the central least-repair problem one may put
\(N=Y^{1/M}\). The target radical is at most \(N\),
the free radical caps are at most \(N^{a_+/a_-}\)
on the balance band \([a_-\lambda,a_+\lambda]\),
and \(C\le N^{1/(a_-\lambda)}\).
For shells \(X\ge N^{1/(2a_-\lambda)}\), the shell budget
\(L=\log N/\log X\) lies in a fixed compact subinterval of
\((0,2a_-\lambda]\). For the short-interval
estimate apply (SI2) at the single budget \(L_*=2a_-\lambda\)
and allowance \(\eta/2\). The caps are at most \(X^{aL_*}\),
and \(L_*/2\le L\) gives \((\eta/2)L_*\le\eta L\).
This ensures uniformity across the shells without an
interchange of budget and height limits.
For smaller \(C\) the unrestricted limiting count with
budget \(3a_+\lambda\) already gives exponent at most
\(a_+/(2a_-)+o(1)\) in \(N\), below one in the balanced
sector. The cutoffs can be imposed independently of the
test primes. Summing the upper shell range or refining
the target multiplicities by shells incurs only a
subpower loss, since the original height is polynomial in
\(N\) at fixed \(\lambda\). All statements concern
primitive quotients.

The real filter is additional input to the modular
incidence question. The prior false Fermat boxes
were uniformly separated from integer equality by a
fixed proportion of the numerical height, so they
do not meet it. A count for popular targets within
the filtered modular incidences is still needed.
[False Fermat-box counts](ABC_THREE_PRIME_FERMAT_NOISE.md#the-three-prime-false-incidence-count)

# Prime tests for popular minimally repaired target radicals

A modular test for the support of one triple and a modular test for
many representations of one target have different thresholds. If
each test block uses primes of total exponent \(\beta\), the
false-positive contribution for popularity threshold \(N^u\) can
be made smaller than the target counting scale with a fixed number
of independent blocks as soon as \(\beta>2-u\). There is also a
lower bound for popular false targets at the critical scale on a
bounded-height Fermat box.

## Support multiplicities and independent blocks

Fix \(D>0\). Let \(N\) tend to infinity and allow positive integer
representatives at most \(N^D\). A support record is a triple of
pairwise coprime squarefree integers \((R_A,R_B,R)\), each at most
\(N\), together with a prescribed eligibility set of triples of
representatives on those exact supports. The sets can include
any height, shape, or least-repair restrictions independent of the
test primes. For each \(R\), let \(m(R)\) count the pairs
\((R_A,R_B)\) for which an eligible triple satisfies \(A+B=C\).
Define \(T(H)=\#\{R:m(R)\ge H\}\).

Fix \(k\ge1\) and positive exponents \(\beta_1,\ldots,\beta_k\),
with \(\beta=\sum_i\beta_i\). In a block, sample one prime in
each \([N^{\beta_i},2N^{\beta_i}]\), uniformly subject to all
primes in that block being distinct. Admit a support record to
the block when, at each sampled prime, some eligible triple
satisfies the equality modulo that prime. The representatives
are allowed to depend on the prime. Write \(m_{\mathcal P}(R)\)
for the number of pairs admitted at a given block \(\mathcal P\).
Take \(j\) independent blocks and let \(P_j(H)\) count the \(R\)
with \(m_{\mathcal P}(R)\ge H\) in all of them. For \(H\ge1\),
\[
 T(H)\le\mathbb E P_j(H)\le
 T(H/2)+O\!\left(
       N\min\{1,N^{\,2-\beta+o(1)}/H\}^{j}
                         \right).                         \tag{MP1}
\]
The implicit terms are uniform in the eligibility sets and in
\(H\), with \(D,k,j,\beta_i\) fixed. The false upper bound depends
only on the support size and polynomial height condition.

To prove this, first fix a false support record. Fixed-support
enumeration gives \(N^{o(1)}\) possible eligible triples,
uniformly at this height. Each of their nonzero discrepancies
\(A+B-C\) has at most \(O_{D,\beta_i}(1)\) prime divisors in the
\(i\)-th test interval. By the prime number theorem and a union
bound on the eligible triples, the probability of acceptance at a
given prime choice of scale \(N^{\beta_i}\) is at most
\(N^{-\beta_i+o(1)}\). Sampling the different primes independently
gives \(N^{-\beta+o(1)}\). Conditioning on distinctness changes this
by at most a bounded factor, since its probability tends to one.
This is the same false-record estimate used in (SO1).

For a fixed target \(R\) there are at most \(N^2\) false support
pairs. Its expected number of false pairs admitted by a block is
at most \(N^{2-\beta+o(1)}\), without any independence between
pairs in that block. If \(m(R)<H/2\), every successful block
requires at least \(H/2\) false pairs. By Markov's inequality the
success probability is at most
\(\min\{1,2N^{2-\beta+o(1)}/H\}\). Independence between blocks
allows this bound to be raised to the \(j\)-th power. Sum over
at most \(N\) possible target radicals. Every true pair passes
every block, giving the other inequality.
[Basic fixed-support modular rejection (SO1)](ABC_MINIMAL_REPAIR_SUPPORT_OCCURRENCE.md#an-average-finite-field-test-equivalent-to-true-occurrence)

For threshold \(H=N^u\) and any \(\beta>2-u\), choose \(j\) with
\(j(\beta+u-2)>u+\sigma\). The error in (MP1) is then
\(O(N^{1-u-\sigma})\). Thus a saving for the expected number of
popular modular targets at that level is equivalent, up to a
factor two in the threshold, to a saving for genuinely popular
targets. For the converse at a fixed exponent one can allow an
arbitrarily small loss, since \(N^{u-\eta}\le N^u/2\) for
fixed \(\eta>0\) and sufficiently large \(N\).

## Application to the repaired height parameter

For a balanced least-repair problem set \(N=Y^{1/M}\).
The repaired target radical is at most \(N\). In a balance
band \([a_-\lambda,a_+\lambda]\) for each original logarithmic
radical relative to \(C\), the two other radicals are bounded
by \(N^{a_+/a_-}\), and every original value is at most
\(N^{1/(a_-\lambda)}\). These assertions follow from
\(R^M\le Y\) and \(C\le R^{1/(a_-\lambda)}\).
One can therefore either slightly enlarge \(N\), or replace
\(2\) in the exponent of (MP1) by \(2a_+/a_-\).
In the second choice the sufficient block threshold from (MP1)
is \(\beta>2a_+/a_--u\). The fixed budgets and the parameters
of the tests are chosen before \(Y\to\infty\). Each block
uses only the original candidates, so quotient representation
counts exclude arbitrary common multipliers.

Uniform fixed-support enumeration makes a true support pair
equivalent to at most \(N^{o(1)}\) primitive quotient witnesses,
uniformly over target radical, up to the exact repaired height
cutoff. Thus one may use support pairs in the multiplicity
criterion after an arbitrarily small adjustment of \(u\).
The upper part \(u\simeq1/2\) of the remaining range in (FL8)
requires block exponent just above \(3/2\) in the balanced
limit. At \(u=0\) this approach still uses block exponent
greater than \(2\) to obtain a strict power error. This does
not bound the expected number of popular modular targets.
[Target multiplicity bound and remaining range (FL5, FL8)](ABC_FIXED_LETTER_RADICAL_CAPS.md#the-remaining-multiplicity-estimate)

## A bounded-height obstruction for one block

Use the dictionaries of the Fermat-box construction
\[
 A=x^d,\quad B=y^d,\quad C=z^d,\qquad
 x,y\in[N,(1+s)N],\quad
 z\in[(1+2s)N,(1+3s)N],\qquad s=(100d)^{-1}.
\]
The three roots are required to be squarefree, mutually
coprime, and units at the sampled primes. In (MP2)–(MP4)
we count only pairs admitted with these unit constraints. This gives
a lower bound for the more permissive test in (MP1). Here \(d\ge3\)
is fixed, so all original heights are \(O_d(N^d)\). The
radicals are exactly \(x,y,z\) and the candidates have
\(\operatorname{rad}(ABC)<C^{3/d}\). For \(M\ge d\), the
least repaired target is \(z^M\). No candidate satisfies
integer equality, since \(2N^d>(1+3s)^dN^d\).

For \(1\le\beta<3/2\), take the prime-scale choices of (FI1).
A proportion bounded below of the tests has
\(\gg N^{3-\beta}\) admitted support records. For
\(3/2\le\beta<38/25\), use the three specific prime scales
of (FI3), for which the same lower bound holds on a
proportion tending to one. In both cases the representatives
are the same at every prime in a block. There is a fixed
\(c>0\) such that on all tests with this lower bound,
\[
 \#\{z:m_{\mathcal P}(z)\ge cN^{2-\beta}\}
       \gg \frac{N^{\,a(\beta)}}{\log^2(2N)},\qquad
             a(\beta)=\min\{1,\,3-3\beta/2\}.              \tag{MP2}
\]
The proof uses a uniform completion bound for a fixed target.
[False records with total test exponent below three halves](ABC_MODULAR_FALSE_INCIDENCES.md)
· [Three-prime false records beyond three halves (FI3)](ABC_THREE_PRIME_FERMAT_NOISE.md#the-three-prime-false-incidence-count)

### The fixed target completion bound

For \(p\nmid dz\), the affine curve
\(x^d+y^d=z^d\) over \(\mathbf F_p\) comes from the smooth
projective Fermat curve. Removing its \(O_d(1)\) points with
a zero coordinate leaves \(O_d(p)\) points. For a nonzero
frequency \((a,b)\) its complete additive Fourier coefficient
is \(O_d(\sqrt p)\). Indeed the linear function \(ax+by\)
has at least one simple pole on the curve at infinity:
the line at infinity has \(d\) distinct intersection
points and the numerator \(aX+bY\) vanishes on at most
one of them. A simple pole rules out a geometric
Artin–Schreier function, so the smooth-curve additive sum
bound applies with all complexity bounded in terms of \(d\).
[Complete additive sums on curves, Bombieri, Section VI, Theorem 5](https://doi.org/10.2307/2373048)

For \(Q=\prod_i p_i\), the Chinese remainder theorem bounds
the Fourier coefficient at \((t_1,t_2)\) by
\(O_{d,k}(\sqrt{Qh})\), where \(h=\gcd(Q,t_1,t_2)\).
If a fixed \(z\) is a nonunit, there are no admitted pairs.
For a unit \(z\), complete the two intervals modulo \(Q\).
Their zero Fourier coefficients are \(O(N)\). On nonzero
frequencies divisible by \(h\) their total absolute mass
in each coordinate is \(O((Q/h)\log(2Q))\), by the geometric
progression bound. Thus the total product mass over the
nonzero vectors with both frequencies multiples of \(h\)
is at most
\[
         O\!\left(\frac Qh\log(2Q)
                  \left(N+\frac Qh\log(2Q)\right)\right).
\]
Multiply by \(\sqrt{Qh}/Q^2\) and sum over proper divisors
of \(Q\). There are at most \(2^k\) divisors and \(Q\gg N\),
so the resulting error is \(O_{d,k}(\sqrt Q\log^2(2Q))\).
The zero mode costs \(O_{d,k}(N^2/Q)\). Dropping the
squarefree and coprimality restrictions only enlarges this upper
bound. Hence for all \(z\),
\[
               m_{\mathcal P}(z)
                     \ll_{d,k} N^2/Q+\sqrt Q\log^2(2Q).
                                                            \tag{MP3}
\]

Now there are \(O(N)\) eligible target radicals \(z\).
Those below \(c N^{2-\beta}\) account for at most half
the \(\gg N^{3-\beta}\) admitted records on taking \(c\)
small enough. Divide the remaining mass by (MP3) to
obtain (MP2). Since \(a(\beta)>\beta-1\) throughout
\(1\le\beta<38/25\), these false records alone make a
one-block popularity bound of order \(N^{1-u}\) fail
at \(u=2-\beta\), allowing constant threshold factors.
The same obstruction holds on an open one-sided
neighborhood \(u<2-\beta\), by the strict exponent
margin.

### Independent repetitions of this obstruction

Choose \(j\) independent blocks of the same test type on the
same Fermat box. Let \(p_z\) be the probability that \(z\)
meets the threshold \(cN^{2-\beta}\) in a block. The positive
proportion assertion and (MP2) give
\(\sum_z p_z\gg N^{a(\beta)}/\log^2(2N)\).
By the power mean inequality over \(O(N)\) possible \(z\),
\[
     \mathbb E P_j(cN^{2-\beta})=\sum_z p_z^j
        \gg_{d,\beta,k,j}
         \frac{N^{\,1-j(1-a(\beta))}}{\log^{2j}(2N)}.
                                                            \tag{MP4}
\]
Thus for \(1\le\beta\le4/3\) the obstruction persists at
the exponent level for every fixed number of independent
blocks. If \(4/3<\beta<38/25\), the displayed lower bound
exceeds the critical scale \(N^{\beta-1}\) whenever
\(\beta<4(j+1)/(3j+2)\). The proof concerns separate
popularity tests on each block. It places no requirement
that the same set of free-support pairs be counted by
different blocks.

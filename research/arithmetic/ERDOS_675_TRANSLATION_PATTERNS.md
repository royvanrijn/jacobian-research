# Erdős #675: a reviewed prime-partition argument and the fixed norm-pattern target

**Project scope.** Deferred stretch goal for later research, kept separate
from the current norm-counting and binary-correlation papers.

**Mathematical status.** The translation property for the fixed set of
sums of two squares remains **open in this work**. The existential
prime-partition question has a valid qualitative Brun-sieve construction,
as in the supplied prior proof posted by Liam Price. This record reviews
that construction, proves an exact reduction of the fixed norm problem,
and separates positive local densities from the missing global existence
assertion. No novelty claim is made for the construction, reductions or
bounded worked examples. The squarefree minimal-translation question is
outside this record.

Let \(E=\{x^2+y^2>0:x,y\in\mathbb Z\}\). For a set
\(A\subseteq\mathbb N\), the translation property means
\[
 \forall N\ge1\ \exists t\ge1\quad
       \mathbf1_A(a+t)=\mathbf1_A(a)\quad(1\le a\le N). \tag{T.1}
\]

## 1. The actual prior proof and the problem's quantifiers

Erdős's [1979 article, Section 5, printed p.69](https://www.renyi.hu/~p_erdos/1979-22.pdf)
asks both the sums-of-two-squares question and whether there **exists**
a partition of the primes into two positive-proportion classes whose
generated sequence has this property. The second assertion is not about
every prescribed partition. The [official problem page](https://www.erdosproblems.com/675)
has an aggregate open label, with a separate
[proof-claims thread](https://www.erdosproblems.com/forum/thread/675/proof-claims).
That aggregate label does not determine the validity of an individual
partial result.

The actual supplied project is *Prime Partitions with Translation*.
The official thread records a claim submitted by Liam Price on
27 July 2026. Its [public read link](https://www.overleaf.com/read/mszmfysmwshf#367bb6)
resolves to the project supplied by the user,
`6a66b64e4a5a746fa9018a16`. The exact source retained as
[a primary source](sources/price-erdos675-prime-partition-2026.tex)
names **GPT 5.6 Sol Pro** as author, and has no internal date.
These are the bytes retrieved on 19 November 2026; the editable project
does not establish that they are identical to its July edition.

Its main convention excludes \(1\). The final remark explicitly
handles the convention including \(1\), by assigning the new prime
factors of \(1+t\) to the allowed class. The proof therefore already
addresses both conventions. It states no result for the fixed set \(E\).

The internal review found no logical gap in its recursive construction,
small-prime locking, disjoint prime assignments, or density control at
intermediate scales. Its explicit numerical bound from Halberstam–Richert
1974, Theorem 10.11, was not checked against the original book pages.
The qualitative bound actually needed follows independently from the
inspected classical lower sieve, as proved next. Thus acceptance of that
numerical formula is unnecessary for the qualitative conclusion below.
This is an internal mathematical review, not external certification or
an update to the official site's status.

## 2. The uniform almost-prime input

**Lemma UniformLinearAlmostPrimes.** For every fixed \(k\ge1\),
there is an integer \(R(k)\) such that for every tuple of positive
integers \(c_1,\ldots,c_k\), and every fixed \(B\), arbitrarily
large positive integers \(u\) satisfy
\[
 \Omega\left(\prod_{i=1}^k(1+c_i u)\right)\le R(k),
 \qquad p\mid\prod_i(1+c_i u)\ \Longrightarrow\ p>B.    \tag{T.2}
\]
Here \(\Omega\) counts factors with multiplicity. The bound
\(R(k)\) is independent of the coefficients and \(B\);
the lower threshold for \(u\) may depend on them.

Put \(F(u)=\prod_i(1+c_i u)\), and let \(\rho(p)\)
be its number of roots modulo \(p\). Since \(F(0)=1\),
\(\rho(p)\le\min(k,p-1)\). For squarefree \(d\), CRT gives
\[
 \#\{U\le u<2U:d\mid F(u)\}
       =U\rho(d)/d+O(\rho(d)),\quad
 \rho(d)=\prod_{p\mid d}\rho(p).                       \tag{T.3}
\]
The usual lower-sieve hypotheses have dimension at most \(k\),
uniformly in the coefficients, and
\[
 V(z):=\prod_{p<z}(1-\rho(p)/p)\gg_k(\log z)^{-k}.
\]
For the last bound, use \(1-\rho(p)/p\ge1/p\) for
\(p\le2k\), and \(1-\rho(p)/p\ge1-k/p\) for larger
primes, followed by Mertens' theorem and the convergent prime-square
correction. The sieve dimension bound follows from \(\rho(p)\le k\).

Apply [Richert 1976, Theorem 11.7, printed pp.141–142](sources/richert-sieve-methods-1976.pdf),
at level \(D=U^{1/2}\) and \(z=U^{1/s}\), with \(s=s(k)\)
a sufficiently large fixed constant. Its relative error can be made
less than \(1/4\). Its weighted remainder is at most
\[
 \sum_{d<D}3^{\omega(d)}\rho(d)
       \ll_k D\log^{3k-1}(2D)=o(UV(z)).                \tag{T.4}
\]
Consequently there are sifted arguments for all sufficiently large
\(U\). If \(C=\max c_i\) and \(U\ge3C\), then
\(F(u)\le(3CU)^k\le U^{2k}\), while every prime factor
is at least \(U^{1/s}\). Thus \(\Omega(F(u))\le2ks\).
Increasing \(U\) makes every prime factor exceed \(B\), proving
the lemma with \(R(k)=\lceil2ks(k)\rceil\).
This proof uses classical sieve input, with no recent ZF or graph theorem.

## 3. The existential prime-partition question

**Theorem PrimePartitionRecurrence.** For every \(0<\delta<1\)
there is a partition of the primes \(P\sqcup Q\) such that
\[
 \pi_P(x)\sim\delta x/\log x,\qquad
 \pi_Q(x)\sim(1-\delta)x/\log x,
\]
and the integers supported on \(P\) have the translation property.
Both conventions for membership of \(1\) can be attained. In the
convention including \(1\), one can enlarge any given baseline
prime set of relative density \(\delta\) by a prime-relative-density-zero
set to achieve the property.

Here is a reconstruction for the convention including \(1\), using
(T.2). Start with any such baseline \(P_0\). Choose integer bounds
\(K_N=\max_{1\le k\le N}R(k)\), let
\(C_N=\sum_{j\le N}K_j\), and choose increasing cutoffs
\(B_N\ge N\) satisfying \(\pi(B_N)\ge N^2 C_N\).
These budgets are available before selecting the stage coefficients.

Suppose the finite set \(S_{N-1}\) of additionally allowed primes
has been chosen. Let \(A_N\) be the integers at most \(N\)
supported on \(P_0\cup S_{N-1}\), including \(1\), and put
\[
 M_N=\operatorname{lcm}(1,\ldots,N)\prod_{p\le B_N}p,
 \qquad F_a(u)=1+(M_N/a)u\quad(a\in A_N).
\]
Every \(M_N/a\) is an integer divisible by every prime at most
\(B_N\). By (T.2), choose \(u_N\ge1\) so that the product
of these forms has at most \(K_N\) prime factors, all exceeding
\(B_N\). Add their distinct prime factors \(T_N\) to the
allowed set, and set \(t_N=M_Nu_N\).

Each old positive position now satisfies \(a+t_N=aF_a(u_N)\),
so remains positive. If an old position was negative, it has an
excluded prime factor \(q\le N\). This factor divides \(M_N\),
hence still divides \(a+t_N\), and no current or future stage can
add it to the allowed primes: all those additions exceed the respective
cutoffs \(B_j\ge j\). Both signs are therefore permanent.

Let \(S=\bigcup_N T_N\) and \(P=P_0\cup S\). For
\(B_N\le x<B_{N+1}\), later stages add no primes at most \(x\),
and the earlier stages added at most \(C_N\) primes altogether. Hence
\[
 \frac{\#(S\cap[1,x])}{\pi(x)}
       \le C_N/\pi(B_N)\le N^{-2}.                    \tag{T.5}
\]
Thus the modification has relative prime density zero, proving the
asserted proportions and recurrence. This density estimate needs no
bound on the largest new prime or on the size of a translation.

For the convention excluding \(1\), Price's retained construction
includes a distinguished form \(F_*(u)=1+M_Nu\), alongside
the forms for positive \(a\ge2\). Its new prime factors can all
be assigned to \(Q\), because a shared prime with \(F_a\)
would divide \(a-1\) while exceeding \(B_N\ge a\).
The positive factors are assigned to \(P\). At free primes, balance
the discrepancy \(d=\pi_P-\delta\pi\) towards the interval
\([-\delta,1-\delta]\). At most \(K_N\) prescribed colors
per stage give \(|d(x)|\le K_N+1\) throughout the stage and
the following free stretch. Choose its starting cutoff with
\(\pi(B_N)>N(K_N+1)\); then \(|d(x)|/\pi(x)<1/N\).
This is the source proof's intermediate-scale density argument. The
qualitative lemma justifies its uniform stage budget, independently of
the unverified explicit numerical choice of \(R(k)\).

This theorem constructs the prime set. It establishes no recurrence
for a prescribed splitting-prime semigroup, or for \(E\), where
prime behavior is fixed and cannot be assigned after the sieve.

## 4. Returns on any prescribed modulus

**Lemma ModularReturn.** If a set has (T.1), then for every \(N\)
and every \(q\ge1\) it has a positive return time for its initial
\(N\)-block divisible by \(q\).

Choose \(q\) positive return times \(r_i\) recursively, making
\(r_i\) return the prefix of length
\(N+\sum_{j<i}r_j\). Every nonempty finite sum of these times
returns the original \(N\)-block: if \(r_i\) is the largest
index in the sum, apply its prefix equality at the argument
\(a+\sum_{j<i,\ j\text{ selected}}r_j\), and proceed inductively.
Among the \(q+1\) partial sums starting at zero, two are congruent
modulo \(q\). Their positive difference is a consecutive sum of
the \(r_i\), hence is the required return.

The hypothesis here is recurrence of **all** initial blocks. A single
return for one fixed block alone does not imply this modular strengthening.

## 5. An exact norm-polynomial reformulation

For this section fix
\[
 A_N=E\cap[1,N],\quad k_N=|A_N|,\quad
 L_N=\operatorname{lcm}(1,\ldots,N),\quad M_N=4L_N^2.
\]
For every \(a\le N\), \(a\mid M_N\), and every prime
\(p\le N\) divides \(M_N/a\). Moreover
\(v_2(M_N)\ge v_2(a)+2\). Define
\[
 Q_N(U)=\prod_{a\in A_N}\left(1+\frac{M_N}{a}U\right)
                                       \in\mathbb Z[U]. \tag{T.6}
\]

**Proposition Erdos675NormReduction.** The sums-of-two-squares
translation property is equivalent to
\[
       \forall N\ge1\ \exists u\ge1:\quad Q_N(u)\in E. \tag{T.7}
\]
For each fixed \(N\) and \(u\ge1\), this polynomial condition
is exactly equivalent to recurrence of the whole initial block with
translation \(t=M_Nu\).

To prove this, for \(p\le N\) one has
\(v_p(a+M_Nu)=v_p(a)\), since the perturbation has larger
valuation. If \(a\notin E\), an inert prime \(p\equiv3\pmod4\)
has odd valuation in \(a\); it is at most \(N\), and its odd
valuation is preserved. Thus every original negative remains negative.

For \(a\in A_N\), the equality
\(a+M_Nu=a(1+(M_N/a)u)\) preserves every inert valuation parity
between the affine form and the endpoint, since \(a\in E\).
At inert primes at most \(N\), all the inner forms are units.
At an inert prime greater than \(N\), at most one endpoint can be
divisible by that prime: a common divisor of two endpoints divides
their nonzero difference of absolute value less than \(N\).
Since the inner factors are units at every prime at most \(N\),
they are in fact pairwise coprime. Consequently the product in (T.6)
has even inert valuations if and only
if every inner factor does. This proves the fixed-progression assertion.
It gives (T.7) \(\Rightarrow\) (T.1); the reverse implication uses
ModularReturn with \(q=M_N\).

The inner affine forms have constant term one, so their product has no
fixed prime divisor: \(Q_N(0)=1\). Their ordinary prime admissibility
is immediate. This observation supplies local compatibility, not prime
values or norm values at a positive integer argument.

## 6. Fully positive local admissibility

**Proposition Erdos675LocalAdmissibility.** For the positive positions
in the preceding progression, every prime has a positive-Haar admissible
set of translation parameters. At an inert prime \(p\le N\)
and at two, every parameter is admissible. At an inert prime \(p>N\),
the exact joint local density is
\[
 \delta_{p,N}=1-\frac{k_N}{p+1}>0,\qquad
 \prod_{\substack{p>N\\p\equiv3\ (4)}}
       \frac{\delta_{p,N}}{(p/(p+1))^{k_N}}
                          \text{ converges positively}. \tag{T.8}
\]

The small inert valuations have already been frozen and are even at
positive positions. At two their valuations and odd units modulo four
are preserved, so the integral two-adic norm condition also remains
satisfied. At split odd primes every local integer is a norm.
For \(p>N\), the affine roots \(-a/M_N\) are distinct modulo
\(p\). Each endpoint has odd valuation with probability
\((1-1/p)\sum_{j\ge0}p^{-2j-1}=1/(p+1)\); these events are
disjoint. This proves the density formula. Its ratio to the product of
single-endpoint densities is \(1+O_N(p^{-2})\), proving the
positive absolute convergence. No infinite simultaneous intersection of
these local conditions on integer parameters follows from this product.

The first question therefore has a particularly explicit admissible
polynomial-norm target (T.7). Arbitrarily many prescribed positive
positions remain to be represented together. The current two-point
comparison theorem does not supply that conclusion. Higher-rank
finite-complexity results do not supply it either merely by adding dummy
variables: the linear parts of these one-variable forms remain parallel.

## Forced prime divisors of a return

**Theorem Erdos675ReturnSize.** Let \(E\) be the positive sums of two squares. Suppose \(t\ge1\) satisfies \(a+t\in E\) whenever \(a\in E\cap[1,N]\). Then every prime \(p\equiv3\bmod4\) with \(200p^{3/2}\le N\) divides \(t\). In particular, for all sufficiently large \(N\),
\[
                         t\ \ge\ \exp(cN^{2/3})
\]
for an absolute constant \(c>0\).

Liu's Lemmas 6 and 7 use two norms in the same nonzero class modulo \(p\) that differ modulo \(p^2\). His Theorem 2 bounds the shift using prime witnesses and Linnik's theorem. We obtain the needed witnesses of size at most \(200p^{3/2}\) using a compactly supported Fourier weight on the norm conic. The norm witnesses need no specified classes modulo \(p^2\) beyond being distinct there. [Liu, Theorem 2 and Remark 12](https://leon2k2k2k.github.io/assets/pdf/erdos/erdos675a.pdf)

Fix \(p\equiv3\bmod4\) and \(r\in\mathbb F_p^\times\). The conic \(x^2+y^2=r\) has \(p+1\) points over \(\mathbb F_p\), by the norm map from \(\mathbb F_{p^2}^\times\). For any nonzero frequency \((u,v)\) its complete sum satisfies
\[
 \left|\sum_{x^2+y^2=r}e_p(ux+vy)\right|\ \le\ 2\sqrt p ,
 \qquad e_p(z)=\exp(2\pi i z/p).
\]
Indeed, insert additive orthogonality and complete the squares in \(x,y\) for each nonzero multiplier \(s\). The product of the two quadratic Gauss sums at \(s\) is the constant \(-p\), so the normalized sum is the negative of
\[
                  \sum_{s\ne0}e_p\!\left(-rs-\frac{u^2+v^2}{4s}\right).
\]
The zero multiplier vanishes since at least one frequency is nonzero. Both coefficients of the Kloosterman sum are nonzero, the second by anisotropy over \(\mathbb F_p\), so the bound follows. This is the prime-modulus case of Tolev's conic calculation. [Tolev, equations (51), (56) and (58)](https://arxiv.org/pdf/1106.3474v2) · [Weil, p. 207](https://www.pnas.org/doi/pdf/10.1073/pnas.34.5.204)

Suppose first that \(p>153^2\), and put \(H=\lceil2p^{3/4}\rceil\). Periodize \(w(x)=(1-|x|/H)_+\) from \(\mathbb Z\) to \(\mathbb F_p\). With the unnormalized Fourier transform, its coefficients are \(H^{-1}|\sum_{j=0}^{H-1}e_p(uj)|^2\). Thus \(\widehat w(0)=H\) and \(\sum_u|\widehat w(u)|=p\), since \(w(0)=1\). Translating the weight leaves the absolute Fourier coefficients unchanged. For any centers \((c_1,c_2)\), Fourier inversion yields the lower bound
\[
 \sum_{x^2+y^2=r} w(x-c_1)w(y-c_2)
          \ \ge\ \frac{p+1}{p^2}H^2-2\sqrt p\ >\ 0 .
\]
Use the centers \((0,0)\) and \((3H,0)\). The first yields a positive norm below \(2H^2\); the second yields one above \(4H^2\) and below \(17H^2\). These bounds use the integer representatives in the support of each weight. They give different integers congruent to \(r\bmod p\), both below \(17H^2\le153p^{3/2}<p^2\). Consequently the two norms differ modulo \(p^2\).

If \(p\le153^2\), choose a modular solution with \(x\not\equiv0\bmod p\), interchanging the coordinates if needed. Lift both coordinates into \([-(p-1)/2,(p-1)/2]\). For a second solution replace \(x\) by \(x-\operatorname{sgn}(x)p\). Its square changes by \(p^2-2p|x|\), which is nonzero modulo \(p^2\). Both norms are at most \((5/4)p^2<200p^{3/2}\). This completes the witness construction for every inert prime.

If \(p\nmid t\), use the two witnesses for \(r=-t\bmod p\). For at least one of them, \(a+t\) has valuation one at \(p\), contradicting \(a+t\in E\). Thus \(p\mid t\) whenever \(200p^{3/2}\le N\). The prime number theorem in arithmetic progressions gives
\[
 \log t \ \ge\ \sum_{\substack{p\le(N/200)^{2/3}\\p\equiv3\bmod4}}\log p
       =\left(\frac{1}{2\cdot200^{2/3}}+o(1)\right)N^{2/3}.
\]
The modulus is fixed at four, as in Liu's argument.

**Corollary Erdos675ReturnValuations.** For a return as in the theorem and every prime \(p\equiv3\bmod4\), an even valuation \(v_p(t)=r\) is impossible if \(200p^{r+3/2}\le N\). An odd valuation \(v_p(t)=r\) is impossible if \(p^{r+1}\le N\). At two, \(v_2(t)\ge\lfloor\log_2N\rfloor-1\).

If \(r\) is even, apply the witness argument to \(t/p^r\) and multiply both positive norm witnesses by \(p^r\), a square. One of the enlarged witnesses gives valuation \(r+1\) in its translate. If \(r\) is odd, \(a=p^{r+1}\) is a square whose translate has valuation \(r\). For the dyadic assertion write \(t=2^r u\) with \(u\) odd. When \(u\equiv1\bmod4\) use \(a=2^{r+1}\), and when \(u\equiv3\bmod4\) use \(a=2^{r+2}\). Both choices are norms, and in each case the odd part of \(t+a\) is \(3\bmod4\), which is incompatible with being a sum of two squares. If \(2^{r+2}\le N\), this contradicts the return.

These bounds apply to the size and divisibility of any return. The existence of a positive return for every initial block in (T.1) remains open.

## 7. An elementary initial-block example

**Proposition Erdos675PellExample.** The initial pattern through
\(N=8\) occurs at infinitely many positive translations. This is a
worked example, with no historical novelty claim.

Let
\[
 z_j+k_j\sqrt3=(1+\sqrt3)(2+\sqrt3)^{504j},\quad j\ge1,
 \qquad m_j=2k_j^2-1,\quad t_j=m_j^2-1.
\]
Then \(z_j^2-3k_j^2=-2\). Coefficientwise,
\((2+\sqrt3)^{18}\equiv1\pmod9\) and
\((2+\sqrt3)^{56}\equiv1\pmod{49}\); hence its power
\(504\) is one modulo \(441\). Thus \(t_j\equiv0\pmod{441}\).
Writing \((z,k,m,t)=(z_j,k_j,m_j,t_j)\), the five positives are
\[
 \begin{array}{lll}
 t+1=m^2,&t+2=m^2+1,&t+4=(m-1)^2+(2k)^2,\\
 t+5=m^2+4,&&t+8=(m-3)^2+(2z)^2.
 \end{array}
\]
The three gaps have \(v_3(t+3)=v_3(t+6)=1\) and
\(v_7(t+7)=1\). The Pell sequence is unbounded, proving the claim.
For a smaller isolated illustration, \(t=288\) already reproduces
the initial eight-bit pattern. Neither this example nor its infinite
family establishes (T.7) for arbitrary \(N\).

## 8. A higher-correlation upper bound, with no existence consequence

**Proposition Erdos675UpperEnvelope.** For every fixed \(N\), with
\(k=k_N\) and \(F_a(u)=1+(M_N/a)u\),
\[
 \#\{1\le u\le X:F_a(u)\in E\text{ for all }a\in A_N\}
                         \ll_N X(\log X)^{-k/2}.       \tag{T.9}
\]
More generally, on a subinterval \(J\subseteq[1,X]\), the upper
bound is \(O_N(|J|(\log X)^{-k/2}+X^{1/20})\).

Here is the actual prime-square sieve budget. At an inert prime \(p>N\),
excluding valuation exactly one at each endpoint excludes
\(\rho(p)=k(p-1)\) residue cells modulo \(p^2\). Write
\(g(p)=k(p-1)/p^2\). The roots are disjoint, and \(g(p)<1\).
For squarefree atoms \(d\), CRT supplies main density \(g(d)\)
and endpoint error \(O(\rho(d))\), where
\(\rho(d)\le k^{\omega(d)}d\). These are necessary norm
conditions, so the abstract Selberg upper sieve applies to these CRT
atoms. This is the quadratic upper-bound identity of Richert's
Theorem 11.5, with the displayed local densities and remainders.

Take \(D=X^{1/100}\), \(z=D^{1/(8k)}\), and use inert
primes \(N<p<z\). Put \(h(p)=g(p)/(1-g(p))\). The
Selberg denominator is
\(G(D)=\sum_{d\le D,\ d\mid\prod p}h(d)\), with squarefree
\(d\). In the full Euler subset law with mass proportional to
\(h(d)\), the inclusion probability of \(p\) is \(g(p)\).
Its mean \(\log d\) is at most \(k\log z+O_N(1)\).
Markov's inequality therefore leaves a fixed positive fraction of the
Euler mass at \(d\le D\). Mertens' theorem in the class three
modulo four gives
\[
 G(D)\gg_N\prod_{N<p<z,\ p\equiv3(4)}(1-g(p))^{-1}
                                  \asymp_N(\log z)^{k/2}.
\]
The optimal Selberg weights obey \(|\lambda_d|\ll_{N,\epsilon}
d^\epsilon\). Consequently their total CRT endpoint error is at most
\[
 \sum_{d,e\le D}|\lambda_d\lambda_e|\rho([d,e])
                     \ll_{N,\epsilon}D^{4+4\epsilon}.
\]
Choose \(\epsilon=1/100\). This is \(O_N(X^{1/20})\),
and the main term is \(O_N(|J|/\log^{k/2}X)\), proving (T.9).
This extends a useful positive envelope to every fixed number of
endpoints; it gives no positive lower bound and hence no solution of
(T.7).

There is also a concrete obstruction to a direct insertion of a third
endpoint into the current marked comparison. On the actual \(N=5\)
progression, \(M_5=14400\), take the split primes
\[
 d=43201,\quad q=115201,\quad p=576001=5q-4.
\]
All are one modulo \(14400\). The raw endpoints are
\(w=p=14400\cdot40+1\) and \(w+4=5q\).
Set \(n=dp\), \(n'=5dq\), and \(R=dq\).
The tentative extra coefficient
\(B(n',Q_0)=\mathbf1_E(n'Q_0-3R)\), with the designated
unshared target mark \(Q_0=q\), equals
\(\mathbf1_E(5q-3)=0\). Swapping target marks replaces it by
\(\mathbf1_E(5d-3)=1\), since
\[
 5q-3=576002=2\cdot7\cdot41143,
 \qquad 5d-3=216002=59^2+461^2.
\]
The target permutation average is therefore \(1/2\), rather than
the required raw value zero. This is a finite control on this proposed
exact encoding. It is not an asymptotic counterexample, nor a proof that
every new multi-endpoint graph construction must fail. The
[graph-control checker](verify_erdos675_graph_controls.py) also checks
the square-atom counts and independent finite controls showing that
positive pair probabilities or bounded positive edge weights alone do
not imply higher intersections or preserve an operator norm.

## 9. Why the nearest norm-value theorems do not supply (T.7)

Browning–Matthiesen–Skorobogatov's
[Theorem 1.2](sources/browning-matthiesen-skorobogatov-pencils-2014.pdf)
treats rational points on systems of norm equations with at least two
homogeneous base variables. Homogenizing our forms as \(v+(M_N/a)u\)
fits that geometry, but produces a rational base parameter \(u/v\).
Clearing denominators rescales the offsets; integrality at the original
offsets requires \(v\mid u\), which the theorem does not supply.

Derenthal–Wei's [Theorem 1.2](sources/derenthal-wei-strong-approximation-1311.3914v2.pdf)
gives strong approximation off infinity, with at least two base variables
and pairwise proportional or affinely independent factors. The slice
\(v=1\) is a different, one-variable system. Its real base component
must also contain balls of arbitrarily large radius. Real norm inequalities
forcing \(v\) into a bounded-width strip fail that requirement, even
if the component is unbounded in another direction. The theorem does
not offer omission of an arbitrary finite inert prime in place of infinity.

Cao–Wei–Xu's [Theorems 1.2–1.3](sources/cao-wei-xu-norm-varieties-1803.11003v5.pdf)
do treat one variable, but their unconditional alternatives require a
norm field to embed in a factor's root field, or its degree to divide that
factor's exponent. Our linear factors have root field \(\mathbb Q\)
and exponent one, so neither alternative holds for \(\mathbb Q(i)\).
Multiplying by squares leaves those exponents odd. Their remaining
cyclic-extension conclusion invokes Schinzel's hypothesis.

Shute's [Theorem 1.1](sources/shute-norm-polynomials-beta-sieve-2025.pdf)
has a quantitative Galois-density hypothesis strictly below \(1/2\).
For the Gaussian norm its bad-cycle proportion \(T(G)\) is \(1/2\).
The theorem's conclusion concerns rational points in any event; it
does not supply a second positive integer parameter from the existing
point \(u=0\).

**Proposition Erdos675NormFieldBarrier.** If a number field \(K\)
contains \(\mathbb Q(i)\), its bad-cycle proportion in Shute's
permutation action is at least \(1/2\).

Let \(\widehat K\) be the Galois closure, and write
\(G=\operatorname{Gal}(\widehat K/\mathbb Q)\),
\(H=\operatorname{Gal}(\widehat K/\mathbb Q(i))\),
\(J=\operatorname{Gal}(\widehat K/K)\). Then \(H\)
is normal of index two, and \(J\subseteq H\). An element outside
\(H\) has only even cycles on \(G/J\): an odd power fixing a
coset would belong to a conjugate of \(J\), contrary to its image
in \(G/H\). Thus at least half of \(G\) has cycle lengths with
nontrivial common divisor, proving the claim. Enlarging the norm field
while retaining Gaussian norms by the tower law cannot repair this
particular theorem's density hypothesis.

There is a straightforward conditional route: Dickson's prime-tuples
conjecture applied to the primitive forms in (T.6) would give simultaneous
prime values. Every coefficient \(M_N/a\) is divisible by four,
so those prime values would be one modulo four and therefore norms.
This proves only an implication from a stronger unproved conjecture.
No use of Dickson or Schinzel is made in the unconditional results here.

## 10. Limits of repairing primes and of a finite squareclass construction

The positive normalized product (T.8) is not positive global Haar measure:
its raw product \(\prod_{p>N,\ p\equiv3(4)}(1-k_N/(p+1))\)
is zero, since the inert reciprocal-prime sum diverges. Finite CRT
compatibility and compactness therefore do not give a positive integer
parameter in (T.7).

The almost-prime lemma also gives no automatic termination after repairing
unwanted inert primes. Forcing previously found primes into square content
\(D^2\) leaves a new primitive residual product to sieve. Its bounded
factor count resets, while the original product's factor count includes
the increasing contribution from \(D^2\). Explicitly, for the single
form \(1+4u\), choose any squarefree product \(D\) of old inert
primes, then fresh distinct inert primes \(p,q\nmid D\), and set
\[
                    u=(D^2pq-1)/4.
\]
This is a positive integer. Every old prime has even valuation, the
residual product has exactly two factors, and two fresh odd inert
valuations remain. Replacing \(D\) by \(Dpq\) permits repeating
the same process. Hence those invariants alone do not imply termination.
This is not a claim that every repair algorithm fails; the single form
itself has many norm solutions.

**Proposition Erdos675FiniteSquareclassObstruction.** For every finite
set \(\mathcal D\) of positive integers, only finitely many integers
\(t\ge1\) have all three of \(t+1,t+4,t+9\) in
\(\bigcup_{d\in\mathcal D}d\mathbb N^2\). If \(1\in\mathcal D\),
that finite union of squareclasses cannot have the global translation
property.

For each such \(t\), choose \(d_1,d_4,d_9\in\mathcal D\)
and integers \(x_1,x_4,x_9\) witnessing the three memberships.
With \(\delta=d_1d_4d_9\), \(Y=x_1x_4x_9\),
\[
 \delta Y^2=(t+1)(t+4)(t+9).
\]
Set \(X=\delta t\), \(Z=\delta^2Y\). These are integral
points on the nonsingular elliptic curve
\[
             Z^2=(X+\delta)(X+4\delta)(X+9\delta).
\]
Siegel's integral-point finiteness applies; an inspected primary proof
is [Noordman, Corollary 4.3, printed p.9](sources/noordman-siegel-2101.07111v2.pdf).
There are only finitely many \(\delta\), so finitely many \(t\).
Global recurrence, however, supplies arbitrarily large returns of any
fixed block: the finite-sums construction in ModularReturn can use
arbitrarily many positive return times. This proves the second assertion.

In particular, a construction confined to one fixed finite collection of
split-prime squareclasses cannot establish arbitrary-prefix recurrence.
The set \(E\) has infinitely many such classes; this proposition
does not rule out its recurrence.

## 11. The current Pell family has no uniform rational identity for the next position

**Proposition Erdos675PellFunctionFieldObstruction.** On the Pell
conic \(z^2-3k^2=-2\), with \(m=2k^2-1\), the rational
function \(m^2+8\) is not a norm from its Gaussian constant-field
extension. Thus the displayed \(N\le8\) family has no rational-function
sum-of-two-squares identity for the additional position \(a=9\).
This assertion does not rule out norm values at infinitely many individual
integer points, or a different family.

Parametrize the conic by
\[
 k=\frac{s^2-2s+3}{s^2-3},\qquad
 z=\frac{-s^2+6s-3}{s^2-3},\quad
 D=(s^2-3)^2,
\]
and set
\[
 \begin{split}
 N(s)&=s^4-8s^3+26s^2-24s+9,\\
 P(s)&=9s^4-16s^3+2s^2+9,\\
 R(s)&=s^4+2s^2-48s+81.
 \end{split}
\]
Then \(m=N(s)/D(s)\) and
\(N(s)^2+8D(s)^2=P(s)R(s)\). Both primitive quartics
remain irreducible over \(\mathbb Q(i)\). Reduction at a Gaussian
prime above five, with the leading coefficients normalized, gives
\[
       s^4+s^3+3s^2+1,\qquad s^4+2s^2+2s+1.
\]
Each has no monic divisor of degree one or two over \(\mathbb F_5\).
The [independent exact certificate](verify_erdos675_constructive_controls.py)
exhausts all sixty such divisor tests and checks the polynomial identities.
Gauss's lemma and reduction prove irreducibility over \(\mathbb Q(i)\).
The quartics are distinct, and neither divides \(D\). The valuation
of \(m^2+8\) at the conjugation-fixed irreducible \(P\) is
therefore one. Every norm \(H\overline H\) has even valuation
there, proving the proposition.

## 12. All positive offsets certified by the Pell family

**Proposition Erdos675PellNormOffsets.** Write \(m=2k^2-1\) in the
rational function field \(K=\mathbb Q(k,z)\) of \(z^2-3k^2=-2\).
For an integer \(a\ge1\), the function \(m^2-1+a\) is a norm from
\(K(i)\) precisely for the positions
\[
 \mathcal S=
 \{1+j^2:j\ge0\}\ \cup\
 \{(1+j^2)^2:j\ge0\}\ \cup\
 \{9j^4-2j^2+1:j\ge0\},                     \tag{T.14}
\]
where \(j\) is integral. For every finite set \(I\) of positive positions
with \(I\cap E\subseteq\mathcal S\), infinitely many positive translations
preserve the membership pattern of \(E\) on \(I\).

Put \(Y=m\), \(U=(Y+1)/2\), \(V=(3Y-1)/2\). Then
\(K=\mathbb Q(Y)(\sqrt U,\sqrt V)\), where the square roots are \(k,z\).
The two squareclasses are independent since \(U,V\) have distinct simple
zeros. The field is rational over \(\mathbb Q\) by the conic
parametrization above, so \(i\notin K\). Let \(A=a-1\).
If \(A=j^2\), we have \(m^2+A=m^2+j^2\). The other two families
have the identities
\[
 \begin{split}
 m^2+j^4+2j^2 &= (m-j^2)^2+(2jk)^2,\\
 m^2+9j^4-2j^2 &= (m-3j^2)^2+(2jz)^2.
 \end{split}                                              \tag{T.15}
\]
In particular every position in (T.14) belongs to \(E\) by evaluating
at \(m=k=z=1\).

Now take \(A>0\) nonsquare. Let \(L=\mathbb Q(\theta)\) with
\(\theta^2=-A\), so \(i\notin L\), and put
\(L'=L(\sqrt{U(\theta)},\sqrt{V(\theta)})\). The substitution
\(Y=\theta+T\) embeds \(\mathbb Q(Y)\) into \(L'((T))\).
Since \(U(\theta),V(\theta)\) are nonzero, the two chosen constant
square roots extend to formal power series square roots of
\(U(\theta+T),V(\theta+T)\) by the binomial expansion at \(T=0\).
They extend the embedding to \(K\).
The image of \(Y^2+A\) is \(2\theta T+T^2\), of order one.
If \(i\notin L'\), the embedding of \(K(i)\) in \(L'(i)((T))\)
intertwines conjugation with coefficientwise conjugation, which
preserves the order in \(T\). Every norm would then have even order.
Consequently a norm requires \(i\in L'\).
The quadratic subextensions
of a multiquadratic extension are exactly those generated by square
roots of the products of the given squareclasses. This follows from
the sign-change automorphisms, with squares omitted if the classes
are dependent. Thus at least one of \(-U(\theta),-V(\theta),-U(\theta)V(\theta)\)
must be a square in \(L\).

For the first possibility, write \(-U(\theta)=(x+y\theta)^2\) with
\(x,y\in\mathbb Q\). The polynomial \(U(Y)+(x+yY)^2\) is divisible
by \(Y^2+A\). Its linear coefficient is \(1/2+2xy=0\), so
\(y\ne0\) and \(h=-x/y=1/(4y^2)\) is a rational square.
Comparison of constant terms gives \(A=h^2+2h\).
Since \(h\) solves a monic integral polynomial it is integral, and
an integral rational square is an integer square. This gives the
second family of (T.14).

The second possibility gives \(3/2+2xy=0\),
\(h=-x/y=3/(4y^2)\), and \(A=h^2-2h/3\).
The rational root theorem applied to \(3h^2-2h-3A\)
gives \(h=u/3\) with \(u\) integral. Now \(u=9/(4y^2)\)
is a square and \(9A=u^2-2u\). If \(3\nmid u\), then
\(u\equiv1\bmod3\) and \(u^2-2u\not\equiv0\bmod3\).
Hence \(u=9j^2\) and \(h=3j^2\), giving the third family.

For the product possibility,
\[
 (x+yY)^2+\frac{3Y^2+2Y-1}{4}=(y^2+3/4)(Y^2+A).
\]
It gives \(xy=-1/4\). With \(s=4y^2>0\), the constant
terms yield \(As^2+(3A+1)s-1=0\). The discriminant
\((3A+1)^2+4A\) is strictly between \((3A+1)^2\) and
\((3A+2)^2\) for integral \(A\ge1\). There is no rational root.
These cases exhaust the possible norms.

For the last assertion, for each \(a\in I\setminus E\) choose an inert
prime \(p_a\equiv3\bmod4\) with \(e_a=v_{p_a}(a)\) odd. Let
\(L_0=\prod_{a\in I\setminus E}p_a^{e_a+1}\), taking the empty
product as one. The unit \(2+\sqrt3\) has finite order modulo \(L_0\)
in \(\mathbb Z[\sqrt3]/L_0\mathbb Z[\sqrt3]\), because its inverse
is \(2-\sqrt3\). For infinitely many positive exponents, its power
is one modulo \(L_0\). The resulting Pell points
\(z+k\sqrt3=(1+\sqrt3)(2+\sqrt3)^j\) have \(k\equiv1\bmod L_0\),
so \(t=m^2-1\equiv0\bmod L_0\), and \(v_{p_a}(a+t)=e_a\) at
every selected negative position. The identities (T.15) and the
square-offset identity preserve the positive positions. These \(t\)
are unbounded.

The norm classification concerns identities in \(K\). Individual
integer specializations may attain additional norm values. In particular,
infinitely many extensions within this Pell family through position 9
and the prescribed-norm problem for arbitrary \(N\) remain open.

## 13. Power substitutions also exclude position 9

**Proposition Erdos675PellPowerCoverObstruction.** Let \(b\in\mathbb Z\),
\(h\ge1\), and take the function field \(K_1=\mathbb Q(X,Y)\) of
\(X^2-3Y^2=1\). Define \(z,k\in K_1\) by
\[
 z+k\sqrt3=(1+\sqrt3)(2+\sqrt3)^b(X+Y\sqrt3)^h .
\]
Then \((2k^2-1)^2+8\) is never a norm from \(K_1(i)\) to \(K_1\).
This includes the substitutions for all arithmetic progressions of the
Pell exponents.

Write \(\varepsilon=2+\sqrt3\). Modulo \(11\), the unit
\(\varepsilon\) has order \(10\): its fifth power is \(-1\), and
\(\varepsilon^2\ne1\). The two base points given by exponents \(1,8\)
satisfy
\[
 (1+\sqrt3)\varepsilon\equiv5+3\sqrt3,\qquad
 (1+\sqrt3)\varepsilon^8\equiv6+3\sqrt3\pmod{11}.
\]
Indeed \(\varepsilon^8\equiv\varepsilon^{-2}=7-4\sqrt3\) here.
At both points \(m=2k^2-1\equiv6\pmod{11}\) and \(m^2+8\equiv0\).
These are simple zeros on the base conic. The derivative of
\(z^2-3k^2+2\) with respect to \(z\) is a unit, so \(k\) is a local
coordinate, and the derivative of \(m^2+8\) is \(8km\ne0\bmod11\).

Write \(h=2^r n\) with \(n\) odd. First use the cover of degree \(2^r\)
with the given shift exponent \(b\). Choose \(e=1\) if \(b\) is odd
and \(e=8\) if \(b\) is even. There is an integer \(j\) with
\(b+2^rj\equiv e\bmod10\). Thus \(X+Y\sqrt3\equiv\varepsilon^j\bmod11\)
maps to one of the two simple zeros. To apply Hensel's lemma, choose
\(\rho\in\mathbb Z_{11}\) with \(\rho^2=3\) and \(\rho\equiv5\bmod11\).
The coordinate \(w=X+\rho Y\) identifies the norm-one conic over
\(\mathbb Z_{11}\) with the invertible elements, using
\(X=(w+w^{-1})/2\) and \(Y=(w-w^{-1})/(2\rho)\).
The base conic has the same coordinate after dividing by
\(1+\rho\), since its norm is \(-2\), a unit here.
In these coordinates the cover is multiplication by a unit times
\(w^{2^r}\), whose derivative at the selected point is a unit.
The zero therefore lifts to a simple zero at some \(w_0\in
\mathbb Z_{11}^{\times}\) of the degree-\(2^r\) pullback.

The full cover substitutes \(w^n\) for \(w\). The polynomial
\(w^n-w_0\) has an irreducible factor of odd degree over \(\mathbb Q_{11}\)
because \(n\) is odd. For a root \(w_1\) of that factor,
\(L=\mathbb Q_{11}(w_1)\) omits \(i\): the irreducible polynomial
\(x^2+1\) over \(\mathbb Q_{11}\) would force an even degree.
Since \(w_1\ne0\), the derivative of \(w^n\) at \(w_1\) is nonzero
in characteristic zero, so the full pullback still has a simple zero.

Finally \(w=w_1+T\) embeds \(K_1\) into \(L((T))\) and maps the
function at position \(9\) to a series of order one. It extends
to an embedding of \(K_1(i)\) into \(L(i)((T))\) that respects
conjugation. Conjugation preserves the order, so every norm has even
order. This proves the obstruction for every \(b,h\).
Individual integer points on these covers may still represent position
\(9\) as a sum of two squares.

## 14. Finiteness on curves above the Pell conic

**Proposition Erdos675PellCurveFiniteness.** Let \(C\) be a smooth
projective geometrically integral curve over \(\mathbb Q\) with
rational functions \(z,k\) inducing a nonconstant map to the projective
completion of \(z^2-3k^2=-2\). If \((2k^2-1)^2+8\) is a sum of two
squares in \(\mathbb Q(C)\), then only finitely many \(P\in C(\mathbb Q)\)
have both \(z(P)\) and \(k(P)\) in \(\mathbb Z\).

For the proof, let \(B\) be the affine Pell conic and let \(C^\circ\)
be its preimage under the induced morphism of smooth projective
curves. The morphism is finite, so \(C^\circ\) is affine and its
coordinate ring is finite over \(\mathbb Q[B]\). Choose finitely
many algebra generators of \(\mathbb Q[C^\circ]\), and for each
choose a monic equation over \(\mathbb Q[B]\). Invert the finitely
many primes in the coefficients of these equations. At every
rational point with \(z,k\) integral, the generator values are
then integral over the resulting localization of \(\mathbb Z\).
They belong to that localization because it is integrally closed
in \(\mathbb Q\). Thus these points form an \(S\)-integral set on
\(C^\circ\) for one fixed finite set \(S\).

Suppose this set is infinite. Siegel's theorem gives genus zero
and at most two geometric points of \(C\setminus C^\circ\).
We use the number-field \(S\)-integral statement, including the
positive-genus case with at most two points at infinity.
[Corvaja–Zannier, Theorem 1 and the following paragraph, p. 268](sources/corvaja-zannier-integral-curves-2002.pdf)

Put \(L=\mathbb Q(\sqrt3)\) and \(v=z+k\sqrt3\). Over \(L\) the
completion of \(B\) is a projective line with coordinate \(v\),
and its two missing points are the zero and pole of \(v\).
Each has a preimage, so \(C\setminus C^\circ\) has exactly two
geometric points. They are \(L\)-rational since the fiber over
each of the two \(L\)-points has a single geometric point and
the characteristic is zero. They are exchanged by conjugation
over \(\mathbb Q\). Pick a rational \(P_0\in C^\circ\) with
\(z(P_0),k(P_0)\) integral. A genus-zero projective curve
with this point is a projective line over \(\mathbb Q\).
Choose a coordinate \(w\) over \(L\) with its zero and pole
above the zero and pole of \(v\), normalized by \(w(P_0)=1\).
For the nontrivial conjugation \(\sigma\) of \(L/\mathbb Q\),
\(w^\sigma w\) is constant by its divisor and equals \(1\)
at \(P_0\). The divisor of \(v\) on \(C\) now gives
\[
 v=v(P_0)w^h\qquad(h\ge1).
\]
Consequently \(X=(w+w^{-1})/2\) and
\(Y=(w-w^{-1})/(2\sqrt3)\) identify \(\mathbb Q(C)\) with the
function field of \(X^2-3Y^2=1\). Indeed both are fixed by
\(\sigma\), and after extending constants to \(L\) they recover
the coordinate \(w=X+Y\sqrt3\).

For an integer point of the Pell conic, its value \(v\) is
\(\pm(1+\sqrt3)(2+\sqrt3)^b\) for some \(b\in\mathbb Z\).
To check this, change its sign and multiply it by a power
of \(2+\sqrt3\) to bring it into
\([1+\sqrt3,5+3\sqrt3)\). These operations preserve integrality
of the coefficients and the norm \(-2\). In this interval,
\(k=(v+2/v)/(2\sqrt3)\) is in \([1,3)\). The integer \(k=2\)
would require \(z^2=10\), so \(k=1\) and the interval selects
\(v=1+\sqrt3\).
We can remove the sign of \(v(P_0)\) by replacing \((z,k)\)
by \((-z,-k)\), leaving \((2k^2-1)^2+8\) unchanged.
The resulting map is exactly a power cover from
Erdos675PellPowerCoverObstruction. Its norm obstruction
contradicts the hypothesis, proving finiteness.

In particular, let \(V\) be the surface given by
\[
 z^2-3k^2=-2,\qquad x^2+y^2=(2k^2-1)^2+8 .
\]
Every proper algebraic subset of \(V\) contains only finitely
many of its integer points. The surface is geometrically integral:
over \(\overline{\mathbb Q}\) its second equation becomes
\((x+iy)(x-iy)=(2k^2-1)^2+8\) over the integral base, with
nonzero right-hand side. On a curve dominating the base the
normalization and the proposition give finitely many integer
points, after discarding the finitely many points where the
normalization is not an isomorphism. On a vertical curve the
base coordinates are fixed, and the real circle bounds the
integer coordinates \(x,y\). Geometric curve components that
are not defined over \(\mathbb Q\) can contribute rational
points only in their finite intersections with distinct
Galois conjugates. Thus any infinite set of integer Pell
specializations with position \(9\) represented would give
Zariski-dense integer points on this surface. The finiteness
assertion on curves leaves the existence of such an infinite
set open.

## 15. Norm identities on a fixed parameter curve

**Proposition Erdos675FixedCurveNormOffsets.** Let \(C/\mathbb Q\)
be a smooth projective geometrically integral curve and let
\(t\in F=\mathbb Q(C)\) be nonconstant. For integers \(1\le a\le A\),
at most \(O_{C,t}(\sqrt A\log(2A))\) of the functions \(t+a\)
are Gaussian norms from \(F(i)\). Fix \(b\ge1\). If \(t+b\)
is nonsquare in \(F(i)\), then at most
\(O_{C,t,b}(\sqrt J\log(2J))\) of the integers \(0\le j\le J\)
give a Gaussian norm \(t+b+j^2\). The number \(r\) of distinct
\(b\)'s for which the entire family \(t+b+j^2\), \(j\ge0\),
consists of norms satisfies \(2^r\mid\deg(t)\). If \(r\ge1\),
the genus of \(C\) is at least \(1+2^{r-2}(r-3)\).
In particular a genus-zero curve certifies at most two
complete families of this form.

Put \(K=\mathbb Q(i)\). For a rational value \(u\), suppose
the fiber of \(t=-u\) on \(C_K\) is the spectrum of a field,
and the morphism is unramified there. Then \(t+u\) has order
one at that unique closed point. Since the fiber descends
to \(\mathbb Q\), complex conjugation fixes this point and
preserves its valuation. Every \(g\bar g\) with \(g\in F(i)\)
has even order there.

We use quantitative Hilbert irreducibility over \(\mathbb Q\)
for monic irreducible polynomials with rational-function
coefficients. At most \(O(\sqrt A\log(2A))\) integers up to \(A\)
can fail to preserve the full arithmetic Galois group,
apart from the finite set of poles and inseparable fibers.
Regularity of the splitting field is unnecessary.
[Zywina, Theorem 1.2, p. 2, attributing the bound to Cohen and Serre](sources/zywina-quantitative-hit-2010.pdf)

Apply this with a primitive element of \(K(C)/\mathbb Q(t)\).
Its minimal polynomial has degree \(2\deg(t)\) and is
irreducible over \(\mathbb Q(t)\). On a nonempty open part
of the base, the algebra of the finite cover is presented
by this polynomial. Indeed the polynomial algebra and the
finite coordinate algebra have the same generic fiber and
become isomorphic after localizing the base. At an integral
specialization that retains irreducibility, the specialized
algebra is a field over \(\mathbb Q\) containing \(K\).
The cover viewed over \(K\) therefore also has a single
closed point in that fiber. Applying this to \(t=-a\), and
discarding the finitely many ramified fibers, proves the
first bound.

For the second bound use the base change \(t=-b-T^2\), with
\(T\) a new rational parameter. Its generic algebra over \(K\)
is a field if \(t+b\) is nonsquare in \(F(i)\). In fact its
function field is \(K(C)(\sqrt{-(t+b)})\), and \(-1\) is
square in \(K\). This field has degree \(\deg(t)\) over \(K(T)\)
and degree \(2\deg(t)\) over \(\mathbb Q(T)\).
The same primitive-element argument at \(T=j\) proves the
bound for \(t+b+j^2\). The finite exceptions also cover the
preimages of all branch values of the original map \(t\).

If all the members of one of these families are norms, the
bound forces \(t+b\) to be square in \(F(i)\). When \(t\)
takes a positive value at a rational point of \(C\), this is
equivalent to \(t+b\) being square in \(F\). Indeed a square
in \(F(i)\) that belongs to \(F\) is a square or the negative
of a square in \(F\). A positive rational value of \(t+b\)
excludes the negative case, since a putative square root
cannot have a pole there. Conversely, if \(t+b=s^2\) with
\(s\in F\), then \(t+b+j^2=s^2+j^2\).

For distinct \(b_1,\ldots,b_r\), the classes of \(t+b_\ell\)
in \(K(t)^\times/K(t)^{\times2}\) are independent, by their
orders at the respective points \(t=-b_\ell\). The same
holds after extending constants to \(\overline{\mathbb Q}\).
Adjoining all their square roots gives a degree \(2^r\)
subextension of \(K(C)/K(t)\) when all the relevant families
are certified. Over \(\overline{\mathbb Q}\) this cover has
ramification index two over each \(t=-b_\ell\) and over
infinity, and is unramified elsewhere. At infinity the
ratios \((t+b_\ell)/(t+b_1)\) are local squares, so the
inertia still has order two. Riemann–Hurwitz gives its
genus \(1+2^{r-2}(r-3)\), and a nonconstant map from
\(C_{\overline{\mathbb Q}}\) to this cover proves the
stated genus bound.

One consequence concerns a fixed finite collection of
curves, with functions \(t_1,\ldots,t_s\) of degrees
\(d_1,\ldots,d_s\). Among the first
\(1+\sum_\ell\lfloor\log_2 d_\ell\rfloor\) positive squares
there is a \(b\) for which every \(t_\ell+b\) is nonsquare
in the respective field after adjoining \(i\). This uses
the degree argument directly for every \(b\) whose function
has a square root. For this \(b\), the norm identities on
the collection miss all but \(O(\sqrt J\log(2J))\) values
of \(j\le J\), although every \(b+j^2\) belongs to \(E\).
Thus a fixed finite list of nonconstant curve identities
cannot give the positive positions for arbitrary initial
blocks. The curves can vary with the block length, and
individual specializations can represent offsets outside
the corresponding function-field norm sets.

**Corollary Erdos675FixedParameterNormOffsets.** The two
counting bounds and the finite-collection obstruction hold
for any geometrically integral parameter variety \(X/\mathbb Q\)
with a nonconstant rational function \(t\). The degree and
genus bounds use the curve whose function field is the
relative algebraic closure of \(\mathbb Q(t)\) in
\(\mathbb Q(X)\).

To see this, let \(L\) be that relative algebraic closure.
It is finite over \(\mathbb Q(t)\), and
\(\mathbb Q(X)/L\) is regular in characteristic zero.
Let \(C\) be the smooth projective model of \(L\).
The dominant rational map \(X\to C\) is a morphism on
a nonempty smooth open of \(X\). Its generic fiber is
geometrically integral by regularity. After restricting
to a nonempty open of \(C\) and its preimage, all fibers
are geometrically integral. We use finite type over a
noetherian integral base and shrink to the open where
the morphism is smooth and surjective.
[Stacks Project, Lemmas 37.27.5 and 37.26.4](https://stacks.math.columbia.edu/tag/0559)
[Geometric reducedness](https://stacks.math.columbia.edu/tag/0578)

For the specializations selected in the proof above, the
fiber on \(C_K\) is an unramified single closed point away
from the finitely many removed values. Its pullback to
the chosen open of \(X_K\) is an integral divisor of
multiplicity one. This follows from flatness and geometric
integrality of the fiber, with \(t+u\) a uniformizer at the
base point. The divisor is fixed by complex conjugation,
so it gives the same odd-order obstruction to a norm
in \(K(X)\). Finally \(K(C)\) is relatively algebraically
closed in \(K(X)\) by regularity after extending constants.
Thus any square root of \(t+b\) in \(K(X)\) already lies in
\(K(C)\), as required in the parabolic and degree arguments.

## Degree and leading coefficient of a polynomial construction

**Theorem Erdos675PolynomialLeadingDivisibility.** Let
\(T\in\mathbb Q[X]\) have degree \(d\ge1\), positive leading
coefficient \(A\), and integer values on \(\mathbb Z\).
Suppose \(T(X)+a\) is a Gaussian norm in \(\mathbb Q(X)\) for
each \(a\in E\cap[1,N]\). Write \(Y=(N/200)^{2/3}\). Then
\[
  \prod_{\substack{p\le Y\\p\equiv3\ (4)}}p\ \mid\ d!A.
\]
The quotient of \(T\) by this product is integer-valued on all
integers. If \(T\in\mathbb Z[X]\), each prime in the product
with \(p>d\) divides every coefficient. A primitive integral
construction, whose coefficients have gcd one, must therefore have
\(d\ge(1-o(1))Y\) as \(N\) grows. At any fixed degree an
integral construction must have
\[
                \log A\ \ge\ (1/2+o(1))Y.
\]

Choose a rational-function norm identity for each prescribed
positive offset. For every sufficiently large integer \(x\),
all the identities specialize without poles and \(T(x)>0\).
A positive rational Gaussian norm that is an integer belongs
to \(E\), since all its inert prime valuations are even. Thus
\(T(x)+a\) belongs to \(E\) at every prescribed positive position.
By Theorem Erdos675ReturnSize, each prime in the product divides
\(T(x)\) for every sufficiently large integer \(x\). Taking the
\(d\)-th forward difference gives \(\Delta^d T(x)=d!A\), an
integer divisible by the product. For the quotient assertion,
shift the parameter origin to a sufficiently large integer and
use Newton interpolation in the binomial basis. All forward
differences at that origin are divisible by the product, and
the binomial basis takes integer values on every integer. If
the coefficients are integral, reduction modulo a prime in the
product gives a polynomial vanishing at every residue. For
\(p>d\) it is the zero polynomial, giving coefficientwise divisibility.

For primitive integral \(T\), all primes in the product must
be at most \(d\). The largest such prime is
\((1-o(1))Y\) by the prime number theorem in the class three
modulo four. At fixed degree the primes with \(p>d\) give
\(\log A\ge\sum_{d<p\le Y,\ p\equiv3\ (4)}\log p
=(1/2+o(1))Y\). More generally for \(d\le Y\) this sum
equals \((Y-d)/2+o(Y)\) uniformly in \(d\), which gives a
tradeoff when both degree and coefficients vary. These conclusions
already follow from the positive positions. They also apply if
one starts with an integral-output polynomial on a fixed
arithmetic progression of parameter values and first substitutes
its affine parametrization into \(T\).

## Quartic polynomial norm identities

**Theorem Erdos675QuarticNormOffsets.** Let \(T\) be an integral quartic
with positive leading coefficient \(A\). The number of integers
\(1\le a\le N\) for which \(T(X)+a\) is a norm from
\(\mathbb Q(i)(X)\) to \(\mathbb Q(X)\) is
\[
 \ll r_2(A)\left(1+\sqrt N+A^{1/8}N^{3/8}\right),
 \qquad r_2(A)=\#\{(u,v)\in\mathbb Z^2:u^2+v^2=A\}.
\]
The implied constant is absolute. A monic integral quartic cannot
certify even the three offsets \(1,2,4\).

Write \(K=\mathbb Q(i)\) and
\(T(X)=AX^4+bX^3+cX^2+dX+e\). If \(T+a\) is a rational-function
norm, it is the norm of a Gaussian integer polynomial
\(Q=\alpha X^2+z_1X+z_0\). For this integrality assertion, split
off the positive integer content of \(T+a\). Its primitive part
is also primitive over \(\mathbb Z[i]\), since a \(\mathbb Z\)-linear
combination of its coefficients is one. In \(K[X]\), a
self-conjugate monic irreducible factor has even exponent in a
norm. Take half of each such exponent and all the factors from
one side of each conjugate pair. By Gauss's lemma, the factors
can be taken primitive over \(\mathbb Z[i]\). The resulting
primitive product and its conjugate multiply to the original
primitive part up to a real unit. Both leading coefficients
are positive, so that unit is one. The content is therefore
itself a rational Gaussian norm. Its inert prime exponents
are even, so it is a Gaussian integer norm by the two-square
criterion. Multiply the primitive product by a Gaussian
integer of that norm.

Fix \(\alpha\in\mathbb Z[i]\) of norm \(A\) and
\(g=\gcd(|\Re\alpha|,|\Im\alpha|)>0\). Write
\[
 z_1/\alpha=B+is,\qquad z_0/\alpha=u+iv,\qquad B=b/(2A).
\]
The relation fixing \(B\) shows that all feasible \(s\) lie in
a coset of \(g^{-1}\mathbb Z\). Indeed the Gaussian integers \(z\)
with \(\Re(z\bar\alpha)=0\) form the lattice
\((i\alpha/g)\mathbb Z\). The coefficient of \(X^2\) gives
\(c/A=B^2+s^2+2u\), with \(Au\in g\mathbb Z\). Hence the
distinct \(s^2\) lie in a coset of \((2g/A)\mathbb Z\).

Use \(X=Y-B/2\) and write
\[
 A^{-1}T(Y-B/2)=Y^4+CY^2+DY+F.
\]
In this coordinate the normalized quadratic has the form
\(Y^2+isY+u'+iv'\), with \(2u'+s^2=C\) and \(2s v'=D\).
If \(s\ne0\), put \(w=s^2\). The constant coefficient is then
\[
             F+a/A=H(w),\qquad
             H(w)=\frac{(C-w)^2}{4}+\frac{D^2}{4w}.
\]
On \(w>0\), \(H''(w)\ge1/2\). On each of its at most two
monotone branches, the \(w\)'s with values in an interval of
length \((N-1)/A\) have diameter at most
\(2\sqrt{(N-1)/A}\). This follows by integrating the inequality
for \(H''\) from the endpoint where the absolute derivative is
smaller. Spacing in \(w\) bounds the number of offsets on that
branch by \(1+\sqrt{AN}/g\). For a fixed sign of \(s\), taking
square roots restricts the \(s\)'s to an interval of diameter
at most \(\sqrt{2}(N/A)^{1/4}\). Their lattice spacing gives the
second bound \(O(1+g(N/A)^{1/4})\) after summing over both signs.
The minimum of these bounds is at most
\(O(1+A^{1/8}N^{3/8})\), using the geometric mean of the two
nonconstant terms.

If \(s=0\), \(D=0\) and the real part \(u'=C/2\) is fixed.
The imaginary part \(v'=v\) lies in a coset of \(g^{-1}\mathbb Z\)
because the real part of \(z_0/\alpha\) is fixed. Now
\(F+a/A=C^2/4+v^2\). On each sign interval the possible \(v\)'s
have diameter at most \(\sqrt{(N-1)/A}\), giving
\(O(1+g\sqrt{N/A})=O(1+\sqrt N)\) offsets. Sum over the
\(r_2(A)\) choices for \(\alpha\). The leading-coefficient obstruction above uses positivity at
all offsets. This estimate counts norm offsets for an arbitrary
integral quartic, including those certifying only some offsets.

For the monic assertion take \(\alpha=1\) by multiplying the
quadratic by a Gaussian unit. Then \(B,s,u,v\) are integers
and \(s^2\equiv c-B^2\bmod2\). If \(s=0\) never occurs,
the positive squares \(s^2\) of this fixed parity are spaced
by at least eight. Each monotone branch has diameter at most
\(2\sqrt3\) when \(a\in[1,4]\). Thus at most two offsets
can occur. If \(s=0\) does occur, put \(u_0=(c-B^2)/2\)
and \(a_0=u_0^2-e\). All feasible \(s\) are even. For
\(s=2k\ne0\), matching the coefficients before translating
gives \(u=u_0-2k^2\) and \(v=Bk\). Thus \(a-a_0\equiv
(Bk)^2\bmod4\). When \(s=0\), \(a-a_0=v^2\) has the same
restriction to zero or one modulo four. The three desired
offsets occupy three distinct residue classes.

**Corollary Erdos675EvenQuarticNormObstruction.** There is a fixed
initial block for which no quartic in \(\mathbb Q[X]\) that becomes
even after a rational translation of \(X\) certifies all positive
positions by Gaussian norm identities. This holds with arbitrary
leading coefficient in that class.

After the translation the norm computation over \(\mathbb Q(i)[X]\)
still applies, with rational coefficients for \(B,s,u,v\).
Polynomial norm factors over this field are obtained by allocating
conjugate factors as above and choosing a leading constant of the
required norm. Evenness gives \(D=0\). Therefore either \(s=0\),
in which case \(F+a/A-C^2/4\) is a rational square, or \(v'=0\),
in which case \(F+a/A\) is a rational square. Each of these two
conditions can hold on no nonconstant four-term arithmetic
progression of offsets, by the theorem on four rational squares.
[Conrad, Theorem 3.4](https://kconrad.math.uconn.edu/blurbs/ugradnumthy/4squarearithprog.pdf)

Choose an arithmetic progression of positive sums of two squares
long enough for the two-colour, four-term case of van der Waerden's
theorem. Such a progression exists by applying Green and Tao's
Theorem 1.2 to the primes one modulo four, which have relative
density one half among primes. If every offset were certified,
colour each position by one of its two square conditions. A
monochromatic four-term subprogression contradicts the preceding
square theorem. The largest element of the chosen progression
gives one fixed initial block for the obstruction.
[Green and Tao, Theorem 1.2](https://annals.math.princeton.edu/wp-content/uploads/annals-v167-n2-p03.pdf)
[Van der Waerden, pp. 212–216](https://cir.nii.ac.jp/crid/1573387449045250304)

For a quadratic \(AX^2+bX+e\), a norm identity instead forces
\(A(e+a)-b^2/4=(\Im(z_0\bar\alpha))^2\) with \(N(\alpha)=A\).
The four norm offsets \(1,5,9,13\) already form a forbidden
progression for that condition. The degree of a nonzero
polynomial norm is even because its order at infinity is even.
These obstructions address identities in the parameter. The
original question concerns integer translations, which can also
come from specializations without such identities.

## 16. Norms on two-squareclass Pell curves

Let \(b<c\) be rational, \(D_1,D_2\) positive rational, and consider
\[
                t+b=D_1 x^2,\qquad t+c=D_2 y^2 .                    \tag{T.18}
\]
The two linear factors have distinct roots, so they give a
biquadratic extension of \(\mathbb Q(t)\) with function field
\(F=\mathbb Q(t,x,y)\). The curve is geometrically integral.

**Proposition Erdos675TwoSquareclassCriterion.**
Fix an offset \(a\), and suppose the curve has a rational seed with
\(t=t_0\) such that \(t_0+a>0\) is a sum of two rational squares.
For \(a\ne b,c\), put
\[
                      U=\frac{a-b}{D_1},\qquad V=\frac{a-c}{D_2}.
\]
Then \(t+a\) is a norm from \(F(i)\) exactly when \(U\), \(V\), or
\(-UV\) is a rational square. For \(a=b\) or \(a=c\), the
corresponding square identity in (T.18) supplies the norm.

For necessity, expand at \(t=-a+s\). The function field embeds in
\(L((s))\), where \(L=\mathbb Q(\sqrt{-U},\sqrt{-V})\). The nonzero
constant terms of the two radicands allow their square roots to be
expanded in \(s\). If \(i\notin L\), conjugation on \(L(i)((s))\)
preserves the order in \(s\), so every Gaussian norm has even order.
But \(t+a=s\) has order one. Thus \(i\in L\), which holds if and only
if one of \(U,V,-UV\) is square, by the quadratic subfields of a
biquadratic field. This also holds when the residue squareclasses
are trivial or dependent, by omitting repeated squareclasses.

If \(U\) is square, \(t+a=D_1(x^2+U)\), and the seed shows
\(D_1\) is a rational norm by taking a quotient of two norms.
The same argument works for \(V\) and at the branch offsets.
If \(-UV\) is square then \(b<a<c\), and
\[
                  t+a=A x^2+B y^2,\qquad
 A=\frac{(c-a)D_1}{c-b},\quad B=\frac{(a-b)D_2}{c-b}.
\]
Since \(AB\) is a square, \(B/A\) is a rational square.
The seed gives \(A\) as a rational norm, and the identity proves
the claim by multiplicativity.

### Obstruction at the three square offsets through 9

**Theorem Erdos675TwoSquareclassSquareAnchorObstruction.**
Suppose (T.18) has a rational point with \(t_0\) a nonnegative integer
and \(t_0+a\in E\) for \(a=1,2,4,5,8,9\). If the curve makes one of
\(t+1,t+4,t+9\) a rational square in its function field, then at least
one of those six positive offsets fails to be a Gaussian norm in
that field. Thus this family of Pell identities cannot certify
the initial block through 9 using any of its three square offsets
as a square.

The nontrivial squareclasses in \(\mathbb Q(t)\) killed by (T.18) are
the two radicands \((t+b)/D_1,(t+c)/D_2\) and their product.
Comparing the valuations at \(t=-b,-c\) shows that making \(t+d\)
square requires a branch at \(d\) with square multiplier, which
we normalize to 1. Call the other branch \(c\), whether it is
above or below \(d\), and write its multiplier as \(D>0\).
At the seed, write \(t_0=m^2-d\) with \(m\) an integer, and
\(w^2=D(t_0+c)\) with \(w\) rational. Multiplying \(D\) by a
rational square scales all the square roots, so we may assume
\(v_3(D)\in\{0,1\}\).

Take first \(d=1\). The offsets \(a=4,8,9\) have
\(a-1=3,7,8\), all nonsquares. They must satisfy one of two
conditions in the proposition:
\[
 k_a^2=D(a-c)\quad (a\ge c),\qquad
 k_a^2=D(a-1)(c-a)\quad (a\le c).                                  \tag{T.19}
\]
At equality either condition uses \(k_a=0\). As \(a\) increases,
the middle-interval condition can switch to the direct condition
at most once. The seed has \(3\nmid m\), because if \(3\mid m\)
then \(m^2+3=t_0+4\) has valuation one at 3.

If all three positions use the direct condition, adding \(w^2\)
to \(k_a^2\) gives \(D(m^2+a-1)\). At \(a=4\) the right side has
valuation \(v_3(D)\). A sum of two 3-adic squares has even
valuation, so \(D\) is a 3-adic unit and \(w,k_4\) are integral at
3. At \(a=9\) the right side is zero modulo 3, forcing
\(w,k_9\equiv0\bmod3\). The equations at \(a=4,8\) then
require \(D\) and \(2D\) both to be nonzero squares modulo 3.

In the mixed case with only 4 in the middle interval, the relation
\(D(4-c)+4D(9-c)=5D(8-c)\) gives
\(k_4^2+15k_8^2=12k_9^2\). Scale to a primitive integer triple.
Reduction modulo 3 gives \(3\mid k_4\). On division by 3, it gives
\(2k_8^2\equiv k_9^2\bmod3\), which forces both remaining
coordinates to be divisible by 3. When 4 and 8 both lie in the
middle interval, the corresponding equation is
\(7k_4^2-15k_8^2=84k_9^2\). First \(7\mid k_8\), and division
by 7 then gives \(k_4^2\equiv5k_9^2\bmod7\). Both coordinates
are divisible by 7, contradicting primitivity again.

If all three positions use the middle-interval condition, the
first two give
\[
 w^2-k_4^2/3=D(m^2+3),\qquad
 w^2-k_8^2/7=D(m^2+7).
\]
In the first equation the two valuation parities are distinct.
If \(v_3(D)=0\), then \(w\) is a unit, \(3\mid k_4\), and
\(D\equiv1\bmod3\). The second equation would have
\(1-k_8^2\equiv2\bmod3\), which is impossible.
If \(v_3(D)=1\), the first equation forces \(3\mid w\).
In the second equation a term of valuation at least two
minus one of even valuation cannot have valuation one.
This exhausts (T.19).

Now let \(d=4\) or \(9\). The other branch must have \(c\le1\)
for the offset \(1\) to pass the real squareclass test. For \(d=4\),
the offsets \(1,2,9\) then require \(k_a^2=D(a-c)\), since
\((a-4)\) is either negative or nonsquare and these offsets
are outside the open middle interval. As before,
\(k_a^2+w^2=D(m^2+a-4)\). The seed condition at \(1\)
forces \(3\nmid m\), since otherwise \(m^2-3\) has valuation
one. At \(a=1\) we get \(v_3(D)=0\), and at \(a=9\) we get
\(3\mid w\). The offsets \(1,2\) would then require \(D,2D\)
both square modulo 3.

For \(d=9\), take \(a=1,2,4,5\), all using the direct
condition. If \(3\nmid m\), then the three seed values
\(m^2-7,m^2-4,m^2-1\) at \(2,5,8\) are all divisible
by 3 and two have valuation one. Thus \(3\mid m\).
The seed value at \(1\) is now a 3-adic unit, so \(D\) is
a unit. All the \(k_a\) and \(w\) in
\(k_a^2+w^2=D(m^2+a-9)\) are integral at 3.
If \(w\equiv0\bmod3\), the offsets 1 and 2 require both
\(D,2D\) square. Otherwise the two square residues at the
classes \(a\equiv1,2\bmod3\) would have to be 0 and 1,
in one order or the other. The pair of offsets in the
zero class gives two \(k_a\)'s divisible by 3. Their squares
differ by \(3D\), an impossibility.

### Necessary prime divisors of the branch difference

**Corollary Erdos675TwoSquareclassBranchDivisors.**
Suppose (T.18) has a rational point with \(t=t_0\) a nonnegative
integer, and certifies \(t+a\) as a Gaussian norm for every
\(a\in E\cap[1,N]\). At every prime \(p\equiv3\bmod4\) with
\(p\ge11\) and \(200p^{3/2}\le N\), both branches are
\(p\)-adically integral and
\[
       b\equiv c\equiv-t_0\bmod p,\qquad
       v_p(D_1)\equiv v_p(D_2)\bmod2,\qquad
       \left(\frac{D_1^*/D_2^*}{p}\right)=-1 ,
\]
where \(D_j^*=D_j p^{-v_p(D_j)}\). In particular the reduced
numerator of \(c-b\) is divisible by the product of these
primes, whose logarithm is
\((1/2+o(1))(N/200)^{2/3}\).

Use the norm witnesses in Theorem Erdos675ReturnSize. Each
nonzero residue modulo \(p\) contains a positive norm at most
\(N\). If \(b\) has negative \(p\)-adic valuation, then at
every integer \(a\) the class of \(U=(a-b)/D_1\) in
\(\mathbb Q_p^\times/\mathbb Q_p^{\times2}\) is \([-1]\).
Indeed \(t_0+b\) is a pole at \(p\), and
\[
         U=-x_0^2\left(1-\frac{t_0+a}{t_0+b}\right).
\]
The parenthesis is a \(p\)-adic square. The norm criterion
would require \(V=(a-c)/D_2\) square for every positive
offset. If \(c\) is also a pole, \(V\) again has class
\([-1]\). If \(c\) is integral, choose a nonzero residue
different from \(c\) at which \(V\) is nonsquare, using
its odd valuation or its nonresidue unit part. Either
case is impossible. The same argument excludes a pole
at \(c\).

Away from \(a=0,b,c\bmod p\), the valuations of \(U,V\) have the
parities of \(v_p(D_1),v_p(D_2)\). If the parities differ, the
product term has odd valuation and the one even-valuation
class must be square at every witness. Half the residue
classes have the opposite unit character before the at
most three exclusions. Suppose the parities agree.
If \(b,c\) are distinct modulo \(p\), each simultaneous
pair of signs of the two unit characters occurs at least
\((p-5)/4\) times away from the two roots. This follows by
expanding the two sign indicators. The sum of the product
of the untwisted characters is \(-1\), by the substitution
\((a-b)/(a-c)\) on \(a\ne c\), and each of the two linear
character sums has absolute value at most one after the
root exclusions. Therefore after excluding \(0\bmod p\)
there is still a nonzero residue where both signs are
negative. None of \(U,V,-UV\) is a square there.
If the roots coincide modulo \(p\) and \(D_1^*/D_2^*\) is
a residue, the two signs are identical, so again they
are both negative at some nonzero residue. This proves
the parity and character requirements.

Finally, if \(t_0+b\) is a unit, then \(t_0+c\) is a unit
of the same residue and (T.18) forces \(D_1^*/D_2^*\)
to be a residue. Hence both \(t_0+b\) and \(t_0+c\) are
zero modulo \(p\). The prime number theorem in arithmetic
progressions counts the resulting divisors of the branch
difference.

## 17. Review and computation scope

The [exact checker](verify_erdos675_reduction.py) independently enumerates
two-square representations in a bounded range, checks locked valuations
and affine coefficients, encloses local densities with exact rational
arithmetic, and verifies the Pell identities and modular returns. It
rejects optimized Python execution, which would disable assertions.
It does not test or prove arbitrary-\(N\) recurrence, a tuple theorem,
or the infinite prime-partition construction. Those written arguments
have separately scoped internal AI review.

The constructive and graph control checkers certify only the stated
polynomial, modular, prime-square-root and finite projection facts. They
are not evidence for arbitrary-prefix recurrence or its negation. The
upper bound and the applicability obstructions leave (T.7) unresolved.

The sources and proof graph retain the fixed \(E\) target as open.
There is no assertion of external human review, formal verification,
independent whole-proof replay, numerical onset or priority. The retained
Price source is credited as prior work; the complete Halberstam–Richert
1974 numerical citation remains outside the completed source audit.

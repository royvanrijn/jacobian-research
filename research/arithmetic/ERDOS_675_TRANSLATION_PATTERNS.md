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
bounded worked examples. The minimal-translation growth question is
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
the simultaneous translation through 9 and the prescribed-norm problem
for arbitrary \(N\) remain open here.

## 13. Review and computation scope

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

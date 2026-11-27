# Erdős #404: factorial sums and unbounded 5-adic cancellation

For an integer \(a\ge1\) and a prime \(p\), let \(f(a,p)\) be the
supremum of \(v_p(\sum_{i\in I}i!)\) over finite sets of distinct
positive indices \(I\) with \(\min I=a\), allowing \(\infty\).
Erdős and Graham ask whether any prime admits an infinite
increasing sequence of indices whose partial factorial sums have
valuations tending to infinity.
[Erdős and Graham, printed pp. 79–80](https://mathweb.ucsd.edu/~ronspubs/80_11_number_theory.pdf)

One sufficient case of [Erdős #404](https://www.erdosproblems.com/404) is
\[
                         f(1,5)=\infty. \tag{404.1}
\]
The exact block calculation below gives \(f(1,5)\ge488272\).

## Uniform lifting obligation

For each precision \(j\), track the reachable residues of sums of **distinct** factorials modulo \(5^j\), with index 1 retained. Find an invariant that (i) permits cancellation of the target residue, (ii) survives the passage from \(j\) to \(j+1\), and (iii) accounts for every newly exposed carry, including at indices around \(5^j\). A finite automaton constructed after fixing a maximum precision is not a uniform invariant. The [existing digit-orbit](image-rigidity/MULTIVARIATE_FACTORIAL_DIGIT_ORBITS.md) and [5-adic automaton](image-rigidity/MULTIVARIATE_PADIC_DIGIT_AUTOMATON.md) notes concern factorial coefficients in polynomial rigidity; their finite residue mechanisms are method analogies, not a proof of (404.1).

If (404.1) is proved, compactness of \(\{0,1\}^{\mathbb N}\) gives an infinite set of indices containing 1 whose 5-adic factorial sum is zero: choose finite cancelling sets at unbounded precisions, take a coordinatewise convergent subsequence, and use that \(v_5(i!)\to\infty\). The limit set cannot be finite, because a finite sum of positive factorials is nonzero. Its increasing enumeration has partial-sum valuations tending to infinity. The same argument applies to any fixed \(a,p\). Conversely, the partial sums of an infinite sequence with valuations tending to infinity imply \(f(a,p)=\infty\) for its starting index \(a\).

## Small completed special case

For every \(a\ge1\) with \(a\equiv6\pmod9\), one has
\[
                         f(a,3)=v_3(a!)+1. \tag{404.2}
\]
After division by \(a!\), the terms at \(a+1,a+2\) are \(7,2\pmod9\); every later term is \(0\pmod9\). The four subset residues are \(1,8,3,1\pmod9\), so valuation one is attained and valuation two is impossible. This is an elementary internal observation; novelty and external review have not been assessed.

## Finite excess in every progression

**Proposition (Arbitrarily large finite factorial-sum excess).**
For each prime \(p\), every arithmetic progression contains a
subprogression on which \(f(a,p)-v_p(a!)\) is a constant finite
integer for positive \(a\). For each \(k\ge0\), these finite
constants are unbounded even within \(a\equiv-2\bmod p^k\).
Along an unbounded sequence of positive \(a\), the finite excess
is at least \((3-o(1))\log_p\log a\).

For \(m\ge1\), consider normalized weights at the integer
argument \(-m\),
\[
 W_j(-m)=\prod_{i=1}^{j}(-m+i)
        =(-1)^j\frac{(m-1)!}{(m-1-j)!}\quad(0\le j<m),
 \qquad W_j(-m)=0\quad(j\ge m).
\]
For \(m\ne2\), no sum of the first \(m\) weights containing
\(W_0(-m)=1\) is zero. To see this, divide by \((m-1)!\)
and write the terms in reverse order. Up to a common sign,
the available terms are \((-1)^\ell/\ell!\) for
\(0\le\ell<m\), with \(\ell=m-1\) required.
If the term at \(\ell=0\) is present, it either cancels the
term at \(\ell=1\), or dominates in modulus the sum of
all terms with \(\ell\ge2\). In the cancellation case at
least one term remains if \(m\ne2\). In a sum whose first
index is \(\ell\ge1\), that first term dominates all
remaining terms in modulus since
\(\sum_{j>\ell}1/j!<1/\ell!\).

Let \(\kappa_m\) be the largest \(p\)-adic valuation of these
nonzero finite sums for \(m\ne2\). Take any \(L>\kappa_m\).
If \(a\equiv-m\bmod p^L\), the products \(W_j(a)\) agree
with \(W_j(-m)\bmod p^L\). For \(j\ge m\) they vanish
modulo \(p^L\) because they contain \(a+m\). The maximum
valuation of normalized positive-index sums containing
\(a!\) is therefore exactly \(\kappa_m\). It suffices to
take \(p^L>3(m-1)!\), since the sum of the moduli of
the negative-start weights is less than \(3(m-1)!\).
Given a progression modulo \(D\), choose any \(-m\) in it
with \(m\ge1\) and \(m\ne2\). The class of \(-m\) modulo
\(\operatorname{lcm}(D,p^L)\) is a subprogression with the
constant excess \(\kappa_m\).

For unbounded finite excess choose \(m=p^s+2\) with
\(s\ge\max(1,k)\). The three indices \(a,a+1,a+3\)
give the normalized sum
\[
 1+(a+1)+(a+1)(a+2)(a+3)=(a+2)^3 .
\]
At \(a=-m\) this has valuation \(3s\), so \(\kappa_m\ge3s\).
Taking \(L\) minimal with \(p^L>3(m-1)!\) and
\(a=p^L-m\) gives \(\log a=O_p(m\log m)\), so
\(\log_p\log a\le s+O_p(\log(s+2))\). These \(a\) are
unbounded and the excess is at least \(3s\).
The two-term lower bound from \(1+(a+1)=a+2\) occurs
in [*A p-adic Reduction for Subsums of Factorials*, Proposition 10.1](https://www.researchgate.net/publication/408593445_A_p-adic_Reduction_for_Subsums_of_Factorials).
The argument here gives large excess on progressions where
the maximum is finite.

The finiteness assertion also holds simultaneously for any
fixed finite set of primes. Use the same \(-m\) and take the
least common multiple of the finitely many moduli and the
original progression modulus. Taking
\(m=2+\prod_{p\in\mathcal P}p^s\) gives simultaneous finite
excesses of at least \(3s\) for all primes in \(\mathcal P\).

The order three at \(-2\) is maximal for any fixed pattern.
For a polynomial \(P(A)=1+\sum_{j\ge1}\epsilon_j W_j(A)\)
with finitely many nonzero \(\epsilon_j\in\{0,1\}\),
\(P(-2)=0\) requires \(\epsilon_1=1\). In that case the
coefficient of \(A+2\) is
\[
  1-\sum_{j\ge2}\epsilon_j(j-2)! .
\]
For this coefficient to vanish exactly one of the terms
\(j=2,3\) must be used, with all others absent. The two
polynomials are \((A+2)^2\) and \((A+2)^3\).
All these polynomials are monic, so a rational root
has to be an integer. The negative-start argument
excludes the other negative integers, and the values
at nonnegative integers are positive. Thus the only
possible rational root is \(-2\).

In particular, the positive integers with \(f(a,p)=\infty\)
cannot fill any residue class modulo a power of \(p\).
This leaves the possible existence of an individual positive
integer with unbounded valuation open.

## Residue obstructions and an exact dyadic class

**Proposition (Dyadic tails miss a residue modulo sixteen).** For
every integer \(a\ge1\), the subsums of \(\{n!/a! : n\ge a\}\)
omit a residue class modulo \(16\). In addition,
\[
 f(a,2)=v_2(a!)+3\qquad(a\equiv4\pmod8).
\]

Only terms before the first zero weight can contribute at a
given precision. The table gives those weights and an omitted
residue. The two rows of precision four separate \(a=4,12\bmod16\).
The weights start at \(1\), with successive multipliers
\(a+1,a+2,\ldots\).

| Start \(a\) | Modulus | Nonzero weights | Omitted residue |
| --- | ---: | --- | ---: |
| \(3\bmod4\) | \(4\) | \(1\) | \(2\) |
| \(1\bmod4\) | \(8\) | \(1,2,6\), in either order after \(1\) | \(4\) |
| \(2\bmod4\) | \(4\) | \(1,3\) | \(2\) |
| \(0\bmod8\) | \(8\) | \(1,1,2,6\) | \(5\) |
| \(4\bmod16\) | \(16\) | \(1,5,14,2\) | \(9\) |
| \(12\bmod16\) | \(16\) | \(1,13,6,10\) | \(9\) |

The omitted class at the displayed modulus excludes all its
lifts modulo \(16\). In the last two rows, the subsums of the
weights after \(1\) omit \(15\bmod16\), so a sum required to
include \(a!\) cannot gain four powers of two. The indices
\(a,a+1,a+3\) give normalized sum \(1+5+2=0\bmod8\)
in both rows. This proves the exact value.

The residue obstruction applies to a strategy that requires
a normalized tail to cover every class. An individual target
residue can still lift through these levels.

**Proposition (Minimum lookback for full residue coverage).**
Let \(v_p(N)=r\), and suppose the preceding block of length
\(W\le p^r\), starting at a positive index, covers every class
modulo \(p^{J_p(N)-J_p(N-W)}\) after division by \((N-W)!\),
where \(J_p(n)=v_p(n!)\). If \(p=2\), necessarily
\(W\ge 2^{r-1}\). For fixed \(p>2\), necessarily
\[
 W\ \ge\ \frac{r-O_p(\log(r+2))}
                   {\log_p 2-1/(p-1)} .
\]
The corresponding leading constant for \(p=5\) is
\(1/(\log_5 2-1/4)>5.53\).

For \(1\le j<W\), \(v_p(N-j)=v_p(j)\), so the precision is
\(r+J_p(W-1)\). Counting subsets and putting
\(s_p(x)\) for the base-\(p\) digit sum of \(x\) gives
\[
 W\log_p 2\ \ge\ r+J_p(W-1)
   =r+\frac{W-1-s_p(W-1)}{p-1}.
\]
At \(p=2\) this forces \(s_2(W-1)\ge r-1\), hence
\(W-1\ge 2^{r-1}-1\). At \(p>2\) the coefficient
\(\log_p2-1/(p-1)\) is positive. If \(W\) is already at least
\(r/(\log_p2-1/(p-1))\) there is nothing to show. Otherwise
\(s_p(W-1)=O_p(\log(r+2))\), giving the bound.

## A local covering criterion at every carry boundary

Let \(N\ge5\) be a multiple of \(5\), let \(A<N\) be a multiple of
\(5\), and put \(J(n)=v_5(n!)\). The *block condition* at \((A,N)\)
means that the subset sums of the integers
\[
  (A+i)!/A!\qquad(0\le i<N-A)
\]
cover \(\mathbb Z/5^{J(N)-J(A)}\mathbb Z\). Assume the sums of subsets
of \(\{2!,3!,\ldots,(A-1)!\}\) cover every residue modulo \(5^{J(A)}\).
For a target residue modulo \(5^{J(N)}\), first realize its residue
modulo \(5^{J(A)}\) with those earlier terms. Divide the difference
by \(5^{J(A)}\). The block sums cover every possible remaining
residue because \(A!/5^{J(A)}\) is a unit. The two sets of factorials
are disjoint, so all residues modulo \(5^{J(N)}\) are then attained.

Consequently it suffices to prove the block condition for **one**
earlier \(A\ge5\) at each multiple \(N>5\). The backwards choice of
\(A\) at each \(N\) avoids any need to reserve disjoint blocks in
advance. Strong induction starts at \(N=5\), where the subset sums
of \(2!,3!,4!\) have all five residues. The resulting full residue
assertion would supply \(-1\) at arbitrarily high precision and hence
(404.1). This condition is sufficient for the full residue route.

There is a finite uniform calculation for the first four jump sizes.
Write \(r=v_5(N)\), and use \(W=N-A\) from the table.

| \(r\) | \(W\) | Precision \(J(N)-J(A)\) | Classes tested |
|---:|---:|---:|---:|
| 1 | 5 | 1 | 1 |
| 2 | 10 | 3 | 4 |
| 3 | 20 | 6 | 100 |
| 4 | 25 | 8 | 500 |

In all four rows \(W\le5^r\). For \(0<j<W\) one has
\(v_5(N-j)=v_5(j)\), so \(J(N)-J(A)=r+v_5((W-1)!)\).
The normalized weights are the integer products
\(\prod_{j=1}^i(N-W+j)\). Their residues therefore depend only on
\(N\) modulo \(5^{r+v_5((W-1)!)}\). Enumerating every class compatible
with \(v_5(N)=r\) covers all positive \(N\) in each row, including
arbitrarily large boundaries. The first row uses the single residue
\(N\equiv0\bmod5\). The
[exact checker](verify_erdos404_local_blocks.py) starts from the
empty sum and updates each bit set by union with its translation
by the next weight. This gives exactly the reachable residues by
induction on the terms. The replay is
`python3 research/factorial/verify_erdos404_local_blocks.py` from the
repository root. The 15-term block at \(N=125\) misses four residues
and serves as a negative control.

The checker also verifies the following individual boundaries.

| \(r\) | \(W\) | Precision | Exact boundaries |
|---:|---:|---:|---|
| 5 | 35 | 12 | \(j5^5\), \(1\le j<25\), \(5\nmid j\) |
| 6 | 40 | 14 | \(5^6,2\cdot5^6,3\cdot5^6,4\cdot5^6\) |

For the size-six jumps it divides the block at \(C=N-25\), using
separate arrays at the lower and upper precision. Set
\(P=5^{J(C)-J(A)}\), \(Q=5^{J(N)-J(C)}\) and
\(U=(C!/A!)/P\), a unit modulo \(Q\). A prefix subsum normalized
by \(A!\) has the form \(b+Pa\) with \(0\le b<P\).
For that low residue \(b\), the upper residues reached by the
whole block, after rescaling by \(U^{-1}\), are the union of
\(U^{-1}a+S\) over the compatible prefix subsums. Here \(S\) is
the normalized suffix subsum set modulo \(Q\). The checker
forms each union by intersecting the translated complements of
\(S\). The intersection must be empty for every \(b\). A 35-term
block at \(N=15625\) fails this test, and the 40-term block passes.
The size-six split uses arrays of length \(Q=5^{10}\), while
enumerating the full prefix subsums with 15 available terms.

The table and the 24 individual boundaries in the Python check
give every residue modulo \(5^{19524}\) using distinct indices
\(2\) through \(78119\).

**Proposition (Full coverage below the first size-nine jump).**
The factorials at indices \(2\) through \(1953119\) have
subsums in every residue class modulo \(5^{488272}\).
In particular \(f(1,5)\ge488272\).

The [word-bitset checker](verify_erdos404_large_blocks.cpp)
implements the same exact split identity for all boundaries
with larger jumps below \(5^9=1953125\).

| Jump size | Width | Suffix width | Boundaries checked |
| ---: | ---: | ---: | ---: |
| 5 | 40 | 35 | 500 |
| 6 | 40 | 25 | 100 |
| 7 | 45 | 25 | 20 |
| 8 | 50 | 35 | 4 |

The suffix support is a cyclic bitset with \(Q\) valid bits.
The checker enumerates all prefix choices, groups their shifts
by the low residue, and intersects the rotated suffix
complements in each group. Taking the complementary
subset pairs low residues \(\ell\) and \(T-\ell\bmod P\),
where \(T\) is the sum of all the normalized weights.
These cosets have identical coverage counts, so testing one
per pair suffices. The width-\(35\) checks at \(15625\) and
\(87500\) are negative controls. To replay from the
repository root, run the Python checker above and

~~~sh
g++ -std=c++20 -O3 -fopenmp research/factorial/verify_erdos404_large_blocks.cpp -o /tmp/erdos404-large
OMP_NUM_THREADS=8 /tmp/erdos404-large 9
~~~

Running the word-bitset checker without the argument \(9\)
certifies the 124 larger jumps below \(5^8\), giving full
coverage through precision \(97648\).
For the full replay use the all-boundary induction with
\(N=1953120\), for which \(J(N)=488272\). Realizing
\(-1\bmod 5^{488272}\) and putting \(1!\) back proves the
lower bound. The all-precision cover would require
controlling jumps of unbounded size. In particular, the
minimum-lookback bound above rules out extending the
finite widths by a rule \(W=5r+O(1)\).

The underlying compactness and integer carry recursion also appear in
[*A p-adic Reduction for Subsums of Factorials*, §§4–5](https://www.researchgate.net/publication/408593445_A_p-adic_Reduction_for_Subsums_of_Factorials).
Its Theorem 5.1 is an exact infinite path criterion. The block
condition above is a sufficient residue test that allows jumps to
be handled by several preceding blocks. The unbounded-jump extension remains open.

## The set of cancelling \(p\)-adic starts

For \(a\in\mathbb Z_p\) define \(W_j(a)=\prod_{i=1}^j(a+i)\) with \(W_0=1\), and set
\[
 \mathcal Z_p=\left\{a\in\mathbb Z_p:
    F_{\epsilon}(a)=\sum_{j\ge0}\epsilon_j W_j(a)=0,\
    \epsilon_0=1,\ \epsilon_j\in\{0,1\}\right\}.
\]
A product of \(j\) consecutive terms in \(\mathbb Z_p\) has valuation at least \(J_p(j)\). For the general setting of factorial polynomial series, see [Anashin and Khrennikov, §3.10.2](https://doi.org/10.1515/9783110203011.1.48). The series therefore converges uniformly in both \(a\) and \(\epsilon\). At a positive integer \(a\), membership in \(\mathcal Z_p\) is equivalent to \(f(a,p)=\infty\), by the compactness argument above after normalization by \(a!\).

**Proposition (Finite root bound for each fixed pattern).**

On the residue disc \(-m+p\mathbb Z_p\), choose a representative \(m\ge1\), \(m\ne2\). Each fixed pattern has at most \(v_p(F_\epsilon(-m))\) zeros on that disc, counted with multiplicity. In particular, for \(p=2\) each pattern has at most three roots on \(\mathbb Z_2\), and all its roots are even.

In \(W_j(-m+px)\), every factor with index congruent to \(m\bmod p\) is divisible by \(p\) as a polynomial in \(x\). Thus the series converges coefficientwise and in the supremum norm on coefficients, since at least \(\lfloor j/p\rfloor\) factors contribute. Its coefficients tend to zero. The coefficient of \(x^i\) is divisible by \(p^i\). At the constant term, the finite negative-start argument gives \(F_\epsilon(-m)\ne0\). The last coefficient of maximal absolute value has index at most \(v_p(F_\epsilon(-m))\). For completeness, if a series \(\sum c_n x^n\) with coefficients tending to zero has last maximal index \(d\) and a root \(z\in\mathbb Z_p\), division by \(x-z\) gives coefficients \(\sum_{i>n}c_i z^{i-n-1}\). They tend to zero and their last maximal index is \(d-1\). This proves the root bound by iteration. For \(p=2\) use \(m=1\) and \(m=4\). The former value is \(1\), and the latter is a subsum with \(1\) required of \(1,-3,6,-6\), whose maximal valuation is \(3\).

**Theorem (Hausdorff dimension of constructed cancelling starts).**

The set \(\mathcal Z_p\) is compact and nowhere dense. Every neighbourhood of \(-2\) contains a compact subset of \(\mathcal Z_p\) of Hausdorff dimension \((p-1)/p\). For \(p=5\), there is such a subset of dimension \(3104/3125\) in every neighbourhood of \(-2\). The dimension is with respect to the metric \(|a-b|_p\).

Compactness follows by projecting the closed zero set in \(\mathbb Z_p\times\{0,1\}^{\mathbb N}\). The excluded negative integers other than \(-2\) are dense, so \(\mathcal Z_p\) has empty interior. To construct subsets, write \(t=a+2\) and
\[
 R_n(t)=\prod_{i=1}^n(t+i).
\]
Fix the first three coefficients of \(F_\epsilon\) to one. Their sum is \(t^2\). At each later offset \(j=n+2\) the term is \((t-1)tR_n(t)\), so it suffices to solve
\[
                   t=(1-t)\sum_{n\in E}R_n(t)                    \tag{404.4}
\]
for a set \(E\) of positive offsets.

Take a power \(L\) of \(p\) and \(1\le D\le J_p(L)\). Suppose that for every sufficiently large multiple \(b=kL\), subsets of the ratios \(R_{b+i}(0)/R_b(0)\) for \(0\le i<L\) cover every residue modulo \(p^D\). Choose one subset for each residue, using the empty subset for residue zero. Fix \(M\ge D\), and let the choice in each block be one digit from \(0,\ldots,p^D-1\). On \(t\in p^M\mathbb Z_p\), the polynomial ratios \(R_{b+i}(t)/R_b(t)\), interpreted as products, retain their residues modulo \(p^D\). The \(L\) factors of the ratio \(R_{b+L}(t)/R_b(t)\) give it valuation at least \(J_p(L)\ge D\).

Use only sufficiently late blocks so all polynomials \(R_n(p^M x)\) in (404.4) have coefficients divisible by \(p^{M+1}\). The right-hand side maps \(p^M\mathbb Z_p\) into itself and is \(p^{-1}\)-Lipschitz there, by subtracting its values in the variable \(x=t/p^M\). Thus iteration starting from \(0\) gives a unique \(t\) for every digit word. A nonzero word cannot have \(t=0\), since its first nonzero block at \(0\) has strictly lower valuation than all subsequent blocks. The negative-start calculation excludes every negative integral start except \(-2\).

For two words first differing at block \(b=kL\), neither root is a zero of \(R_b\). At either root \(t\), the two right-hand sides have difference of valuation
\[
 \gamma_t(kL)+e,\qquad
 \gamma_t(n)=\sum_{i=1}^n v_p(t+i),\qquad 0\le e<D.              \tag{404.5}
\]
The first differing block gives the value in this range by its distinct residue. Subsequent blocks are divisible by \(R_{(k+1)L}(t)\) and cannot change that valuation. The contraction inequality shows that (404.5) is also the valuation of the difference of the two fixed points. Thus the digit words map injectively to the roots. The bound \(\gamma_t(n)\ge J_p(n)\) gives continuity.

For the dimension, if the first \(m\) digits are fixed, the resulting set has diameter at most \(p^{-J_p((k_0+m)L)}\), where \(k_0\) indexes the first usable block. There are \(p^{Dm}\) such sets, so the dimension is at most \((p-1)D/L\). For the lower bound, give the digits independent uniform probabilities. For any \(\delta>0\) and all large positive integers \(i\), roots with \(v_p(t+i)\ge\delta i\) agree on all digits through block starts
\[
                kL\ \le\ \frac{\delta i}{2\log_p i}.
\]
Indeed every factor \(t+j\) at such a start, \(1\le j\le kL<i\), has the valuation of the nonzero integer \(j-i\), at most \(\log_p i\). A first difference there would give a valuation at most \(\delta i/2+D-1\) by (404.5), whereas any two roots in the indicated ball differ by a multiple of \(p^{\lceil\delta i\rceil}\). The probabilities of these events are therefore summable in \(i\), since each is at most \(p^{-c_\delta i/\log i}\) for some \(c_\delta>0\). Counting the tail union of the events and taking \(\delta\) through positive reciprocal integers gives \(v_p(t+i)=o(i)\) for almost every root in this measure. For such \(t\), \(\max_{1\le i\le n}v_p(t+i)=o(n)\). Counting the multiples of \(p^h\) among the \(n\) factors, with an error of at most one for each nonzero level, now gives \(\gamma_t(n)=n/(p-1)+o(n)\).

For almost every \(t\), a ball of radius \(p^{-K}\) around \(t\) forces agreement on \((p-1)K/L+o(K)\) initial digits, by (404.5). Its mass is at most \(p^{-((p-1)D/L)K+o(K)}\). For any \(s<(p-1)D/L\), restrict to a positive-measure set on which the bound \(p^{-sK}\) holds for all \(K\) beyond one fixed threshold. Any sufficiently small ball meeting this set can be centred at a point in the set without changing the ball. Summing the ball bounds over any cover proves that the Hausdorff dimension is at least \(s\), which gives the claimed dimension as \(s\) increases.

For every prime take \(L=p\) and \(D=1\). Inside a block the normalized weights modulo \(p\) are \(0!,1!,\ldots,(p-1)!\). Each is a unit. Adding a new nonzero class to a set of subsum residues strictly enlarges it unless it already equals \(\mathbb F_p\), because a nonzero translation generates the group. Thus every residue is attained.

For \(p=5\), take \(L=3125=5^5\) and \(D=J_5(L-5)=776\). Inside any interval \([b,b+L)\) with \(L\mid b\), every positive multiple of five after \(b\) and before \(b+L\) has jump size at most four. The uniform checks in the [local covering criterion](#a-local-covering-criterion-at-every-carry-boundary) then give full residue coverage through \(b+L-5\) by induction starting at \(b\) with precision zero. The chosen block for each jump starts at or after \(b\) because the first distance to a size-\(r\) boundary is at least \(5^r\), at least the tested width. Also \(J_5(b+L-5)-J_5(b)=J_5(L-5)\), by equality of the valuations at relative indices \(1,\ldots,L-5\). Hence (404.4) constructs a set of dimension \(4\cdot 776/3125\). These dimension bounds concern \(p\)-adic starts. The existence of a positive integer in \(\mathcal Z_p\) remains open.


## How large a correction must be near a prime power

An ordinary block at a multiple \(n\) of \(5\) has some spare
terms. The three weights \((n+1)!,(n+2)!,(n+3)!\), normalized by
\(n!\), have residues \(1,2,1\) modulo \(5\), so their subsums
realize every residue modulo \(5\). The two unused terms can be
selected as one packet with valuation \(v_5(n!)+2\). Indeed
\[
 \frac{n!+(n+4)!}{n!}
   =25+50n+35n^2+10n^3+n^4\equiv25\pmod{125}.
\]
Such packets can add choices at an offset of two. Large jumps
require packets that reach much further.

**Proposition (Size of corrections at large jumps).**
Fix a prime \(p\), write \(J_p(n)=v_p(n!)\), and set \(N=p^r\).
For sufficiently large \(r\), consider a signed sum of distinct
positive factorials whose valuation lies in
\[
 \left[J_p(N-1)+r/3,\ J_p(N-1)+2r/3\right].                  \tag{404.3}
\]
Its nonzero terms have index span \(\gg_p r/\log r\), where
span is largest index minus smallest index plus one. The sum
also contains more than \(\log r/(10\log\log r)\) nonzero terms.
Both statements allow coefficients from \(\{-1,0,1\}\).

Terms with index at least \(N\) can be dropped because their
valuations are at least \(J_p(N)=J_p(N-1)+r\).
Write \(n\) for the smallest remaining index and \(W\) for the
span of the remaining terms. Set \(w=N-n\). The sum equals
\(n!P(n)\), where
\[
 P(Z)=\sum_{j=0}^{W-1} e_j \prod_{\ell=1}^j(Z+\ell),
 \qquad e_j\in\{-1,0,1\},\quad e_0\ne0.
\]
Each factor evaluated at \(n\) is at most \(N\). Thus, if
\(s=v_p(n!P(n))\) lies in (404.3), then
\[
 s-J_p(n)\le (W-1)r+\log_p W .
\]
Counting multiples of \(p\) between \(n\) and \(N-1\) also gives
\(J_p(N-1)-J_p(n)\ge\lfloor(w-1)/p\rfloor\).
It follows that \(w\ll_p(r+1)W\).

The integer \(P(n)\) is congruent to \(P(-w)\) modulo \(p^r\).
If \(W\le c_p r/\log r\) for a sufficiently small constant
\(c_p>0\) and \(P(-w)\ne0\), then
\[
 v_p(P(-w))\le \log_p \bigl(W(w+W)^{W-1}\bigr)<r/3.
\]
The congruence then gives the same valuation to \(P(n)\),
which contradicts (404.3). If \(P(-w)=0\), we must have
\(w\le W\). For \(w\ge W+1\) the absolute values of the
successive factorial products at \(-w\) increase by at least
a factor of two. The largest product with nonzero coefficient
then exceeds the sum of the absolute values of all the others.
When \(P(-w)=0\) the congruence gives \(v_p(P(n))\ge r\), while
\[
 J_p(N-1)-J_p(n)=J_p(w-1)\le (w-1)/(p-1).
\]
The equality holds since \(v_p(N-j)=v_p(j)\) for \(0<j<N\).
For \(w\le W\le c_p r/\log r\) this places the valuation above
the upper end of (404.3). This proves the span bound.

For the number of terms, fix one attainable valuation \(s\)
in (404.3) and take a signed sum attaining it with the smallest
number \(k\) of terms. Its indices satisfy
\(n=n_1<\cdots<n_k<N\). Every proper prefix sum \(S_i\)
has valuation strictly less than \(s\). If a prefix had
valuation \(s\), that prefix alone would suffice. If its
valuation exceeded \(s\), the suffix alone would suffice
by the ultrametric inequality. In particular every such
prefix is nonzero and has \(J_p(n_{i+1})\le v_p(S_i)\),
since all later factorials have valuation at least
\(J_p(n_{i+1})\). Write \(W_i=n_i-n+1\), and write
\(S_i=n!P_i(n)\) using the polynomial construction above.

First use \(|P_i(n)|\le W_i N^{W_i-1}\), together with
\(J_p(n_{i+1})-J_p(n)\ge\lfloor(n_{i+1}-n)/p\rfloor\).
For some constant \(C_p\) this gives
\(W_{i+1}\le C_p(r+2)W_i\) for every \(i<k\).
The analogous final bound is \(w\le C_p(r+2)W_k\), as
for the span argument. Since \(W_1=1\), we obtain
\(w,W_i\le (C_p(r+2))^k\).
This preliminary estimate implies
\(\log_p(w+W_i)=O_p(k\log(r+2))\).

Suppose now that \(k\le\log r/(10\log\log r)\). Put
\(L=C'_p k\log(r+2)\), with \(C'_p\) a sufficiently large
constant. For all sufficiently large \(r\), \(L^k\le r^{1/3}\)
because \(\log L\le(2+o(1))\log\log r\).
We prove \(W_i\le L^{i-1}\) successively, starting at \(i=1\).
If \(P_i(-w)=0\) at a step with this bound, its span forces
\(w\le W_i\). As in the span argument, this would give
\(v_p(S_i)\ge J_p(N)-W_i/(p-1)>s\), which is impossible
for either a proper prefix or the full sum. Otherwise
\[
 v_p(P_i(-w))
 \le \log_p\bigl(W_i(w+W_i)^{W_i-1}\bigr)
 \ll_p k\log(r+2)W_i.
\]
The last bound is less than \(r/3\), by the choice of \(L\)
and \(L W_i\le L^k\). Hence it bounds \(v_p(P_i(n))\)
as well. For \(i<k\) the next index inequality above and
counting multiples of \(p\) give \(W_{i+1}\le L W_i\) after
choosing \(C'_p\) large enough. For \(i=k\) the valuation
bound puts \(s\) below \(J_p(N-1)+r/3\), a contradiction.
This proves the term count.

The proposition rules out making standalone corrections at
the middle digits from a uniformly bounded number of factorials.
It applies also to differences between two candidate fillings.
A full lifting argument must account for corrections of growing
size, or avoid requiring such standalone corrections. The
existence of an infinite carry path for (404.1) remains open.

## Prime 41 and equal low residues

**Theorem (Prime 41 before the first size-eight jump).**
The factorials at indices \(2\) through \(41^8-42\) have
subsums in every residue modulo \(41^{199623130720}\).
In particular \(f(1,41)\ge199623130720\).

For a fixed prime \(p\), write \(J(n)=v_p(n!)\). It suffices
at each boundary \(N\equiv0\bmod p\) to find an earlier positive
multiple \(A\) of \(p\) such that the subset sums from indices
\(A,\ldots,N-1\), divided by \(A!\), cover the modulus
\(p^{J(N)-J(A)}\). In fact, use induction to realize the
desired residue modulo \(p^{J(A)}\) before index \(A\).
After subtraction and division by \(p^{J(A)}\), multiplication
by the unit \(A!/p^{J(A)}\) leaves the entire local target group
available. This is the covering implication from the five-adic
argument. The initial boundary \(p=41\) is covered by indices
\(2,\ldots,p-1\). The subset sums of \(2!,3!\) give four distinct
residues because these two units agree neither with each other
nor with each other's negatives. Each of the remaining \(p-4\)
terms enlarges the support by at least one unless it already
covers the group.

For \(r=v_p(N)\le6\), take \(A=N-p\). The needed precision
is \(r\). The normalized weights
\[
                 w_j(N)=\prod_{i=1}^j(N-p+i)\quad(0\le j<p)
\]
depend only on \(N\bmod p^r=0\) when reduced modulo \(p^r\).
They are consequently also the reductions of the weights at
\(N=p^6\) modulo \(p^6\). The one-block test at prime \(41\)
covers this entire modulus. The same block width works at all
boundaries with jump at most six.

For \(v_p(N)=7\), take \(A=N-2p\), \(C=N-p\) and \(Q=p^7\).
Now \(J(C)-J(A)=1\), \(J(N)-J(C)=7\).
Let \(S\subseteq\mathbb Z/Q\mathbb Z\) be the subset support
from \(C,\ldots,N-1\) normalized by \(C!\). It is the same
at every boundary divisible by \(p^7\). In the prefix block
\(A,\ldots,C-1\) write
\[
 u_j=\frac{(A+j)!}{A!},\qquad
 \lambda=\frac{C!}{pA!},\qquad 0\le j<p.
\]
Choose exactly one term from each of the disjoint pairs of offsets
\[
 (0,1),\ (5,7),\ (6,19),\ (8,22),\ (9,26),\ (11,15),\
 (18,32),\ (23,30).
\]
Each pair has equal values modulo \(p\). The unpaired terms
have subset sums attaining every class modulo \(p\).
Thus any requested low digit can be attained while preserving
all the choices within the pairs. After selecting one fixed
representative from each pair and any unpaired subset with the
required low digit, all pair switches change only the upper residue.
In units normalized by \(C!\), their changes are
\[
                    \lambda^{-1}(u_j-u_i)/p\bmod Q.
\]
The sums of subsets of these eight changes, added to \(S\),
cover the entire upper group in every one of the forty classes
for \(N/p^7\bmod p\) with \(v_p(N)=7\). The changes depend only
on \(N\bmod p^8\), which determines the \(u_j\bmod p^8\) and
the unit \(\lambda\bmod Q\). The upper offset from the chosen
representatives and unpaired subset only translates the group.

The [exact checker](verify_erdos404_prime41_packets.cpp) computes
the one-block supports by cyclic bitset updates. It tests the
paired and unpaired conditions and the forty finite translations
of the complement. The first tested nonzero shift is from pair
\((0,1)\). Its reduction modulo \(p^6\) is shared across the forty
classes. To test the first intersections, the checker groups the
missing residues of \(S\) by their low part modulo \(p^6\). For a
shift with low part \(s\), two missing residues \(a,b\) can have
\(a-b\) equal to the shift only if their low parts differ by \(s\).
A sorted merge finds these pairs of low parts, and the high digits
identify the classes whose first shift equals \(a-b\bmod Q\).
Direct membership tests apply the remaining shifts.
The seventh-power one-block support has \(194716659405\) residues
and misses \(37614476\). The paired two-block test closes those
omissions in the forty boundary classes.

From the repository root, replay with GCC 13.3.0 and OpenMP.
The largest cyclic update uses two bitsets occupying
approximately \(45.4\) GiB in total.

~~~sh
g++ -std=c++20 -O3 -march=native -Wall -Wextra -fopenmp \
  research/factorial/verify_erdos404_prime41_packets.cpp \
  -o /tmp/verify_erdos404_prime41_packets
OMP_NUM_THREADS=12 /tmp/verify_erdos404_prime41_packets
~~~

Induction applies at every boundary through \(N=p^8-p\).
There \(J(N)=\sum_{i=1}^{7}p^i-7=199623130720\).
Selecting residue \(-1\) and putting back \(1!\) proves the
lower bound. A uniform extension through jumps of unbounded
size remains open. The first size-eight boundary requires
an additional argument.

## A long finite dyadic carry

**Proposition (Exact dyadic excess at start 34).**
For \(1\le a\le80\), the value \(f(a,2)\) is finite.
The largest excess \(f(a,2)-v_2(a!)\) in that range occurs at
\(a=34\), where it is \(18412\). The first boundary at which
the carry from \(34!\) cannot lift is \(18448\).

Odd \(a\) have excess zero since the next factor \(a+1\) is
even. For even \(a\), group indices as \(n,n+1\), starting
at \(n=a\). Put \(J(n)=v_2(n!)\). Process two consecutive
factorials at a time, representing earlier partial sums after
division by \(2^{J(n)}\) and by the fixed odd part of \(a!\).
If the unit weight at \(n\) is \(u\), the two weights are
\(u,(n+1)u\). The first slot at \(a\) is required. At the next
boundary \(n+2\), a residue can continue only if its bottom
\(v_2(n+2)\) bits vanish, in which case those bits are divided
out. The new unit weight is
\[
                    (n+1)\operatorname{odd}(n+2)u,
\]
where \(\operatorname{odd}(m)=m/2^{v_2(m)}\).
For cutoff \(B\), store only residues modulo \(2^{J(B)-J(n)}\)
at stage \(n\), preserving every divisibility test through \(B\).
At the first empty successor set, all attained valuations lie
strictly below \(J(n+2)\). Their maximum is \(J(n)\) plus the
maximum valuation of the discarded low remainders. Earlier
rejections have valuation less than \(J(n)\), and every
subsequent factorial has valuation at least \(J(n+2)\).
Consequently the finite failure gives the exact supremum even
when arbitrary later indices are allowed.

The [dyadic checker](verify_erdos404_dyadic_excess.py) tests
the forty even starts through \(80\), with cutoff \(4096\)
except for start \(34\), where the cutoff is \(18448\).
It finds every other maximum excess at most \(1522\).
At \(34\) it gives \(J(18446)=18441\) and maximum valuation
\(18444\), with \(J(34)=32\). Replay from the repository root.

~~~sh
python3 research/factorial/verify_erdos404_dyadic_excess.py
~~~

The exact excess at \(34\) also holds for its residue class
modulo \(2^{18413}\), since the normalized products
\(\prod_{i=1}^j(a+i)\) agree at that precision throughout the
class. The existence of an infinite carry path at some
prime remains open.

## Cancelling starts from variable residue alphabets

**Theorem (Dyadic cancellation dimension and excess density).**
Every \(2\)-adic neighbourhood of \(-2\) contains a compact
subset of \(\mathcal Z_2\) of dimension \(d_2>0.91598\),
where \(d_2\) is the periodic alphabet exponent in (404.6)
below. For fixed \(m\ge1\), the positive integer starts \(a\)
with \(f(a,2)-J_2(a)\ge m\) have a natural density
\(\Delta_2(m)\); the condition includes infinite excess.
There are constants \(c,C>0\) such that
\[
              \Delta_2(m)\ge c(m+1)^{-C}2^{-(1-d_2)m}.
\]
At \(p=41\) the corresponding conclusions hold with
\(d_{41}=1-321/41^8\) in place of \(d_2\), and \(41\) in place
of the base \(2\) in the density estimate.

The dyadic construction uses all the distinct residues available
at each block. Its alphabet changes with the valuation of the next
boundary. Blocks can contribute even when they miss some residues.

### A periodic alphabet criterion

Use the fixed point equation (404.4) in the variable \(t=a+2\).
Fix a block length \(L\), put \(b_k=kL\), and choose
\(1\le D_k\le D_{\max}\). Suppose that on \(t\in p^M\mathbb Z_p\),
where \(M\ge D_{\max}\), the normalized product between consecutive
block starts is divisible by \(p^{D_k}\):
\[
                      p^{D_k}\mid\prod_{i=1}^L(t+b_k+i).
\]
Let \(A_k\) be the attainable subset residues of
\(\prod_{j=1}^i(b_k+j)\), \(0\le i<L\), modulo \(p^{D_k}\),
with the empty product for \(i=0\). Suppose \(|A_k|\ge2\)
and the alphabet sizes are periodic of period \(P_0\). Define
\[
                 d=\frac{p-1}{LP_0}
                            \sum_{0\le k<P_0}\log_p|A_k|.       \tag{404.6}
\]
All of these products give the same normalized weights at
\(t\) modulo \(p^{D_k}\), since \(M\ge D_{\max}\).
Choose one subset for each residue, taking the empty subset
for zero. Use only blocks from a fixed sufficiently large
\(k_0\) onward. The right side of (404.4) then maps
\(p^M\mathbb Z_p\) to itself with Lipschitz constant at most
\(p^{-1}\), as in the fixed-alphabet construction. Each word gives
one root. No such root is a negative integer in the \(t\)-variable:
that would give the excluded negative integer start \(a=t-2\ne-2\).
The zero root in the \(t\)-variable is obtained from the all-zero word.

For two words first differing at \(k\), apply the same contraction
comparison at either root \(t\). The two sums differ there in
valuation \(\gamma_t(b_k)+e\), where
\(\gamma_t(n)=\sum_{j=1}^n v_p(t+j)\) and \(0\le e<D_k\).
The later blocks start at valuation at least
\(\gamma_t(b_k)+D_k\), so they leave this valuation unchanged.
Comparison of the two fixed points gives this exact valuation
for their difference. The map from words is injective.
Agreement through \(r\) successive blocks gives distances
at most \(p^{-J_p(b_{k_0+r})}\). Thus the map is continuous,
its image \(C_p\) is compact, and the usual prefix covers
give dimension at most \(d\). Indeed the number of prefixes
of \(r\) digits is \(p^{d Lr/(p-1)+O(1)}\) by periodicity.

Give the digits independent uniform probabilities. To get
logarithmic control of the metric under this measure, fix
\(i\ge2\) and take \(h_i=\lceil A\log_p(i+1)\rceil\) with \(A\)
a sufficiently large constant. If \(v_p(t+i)\ge h_i\), then
for \(n<i\) the factors \(t+j\) have the same valuations
as \(i-j\), since \(v_p(i-j)< h_i\). Legendre's formula gives
\[
 \gamma_t(n)=J_p(i-1)-J_p(i-1-n)
                  \le \frac{n}{p-1}+\log_p i+1.
\]
Two roots in the ball \(v_p(t+i)\ge h_i\) must therefore
agree on all blocks with \(b_k\le(p-1)h_i/2\), for sufficiently
large \(i\). A first difference there would have valuation
at most \(h_i/2+\log_p i+1+D_{\max}<h_i\).
The ball has measure \(O(p^{-d h_i/2})\) by the prefix count,
and these bounds are summable if \(Ad/2>2\).
By Borel–Cantelli, almost every \(t\) satisfies
\(v_p(t+i)\le C(t)\log_p(i+1)\) for all \(i\).
Choose a fixed \(C_0\) for which the set of such
\(t\) with \(C(t)\le C_0\) has positive measure.

On this set, count the multiples of \(p^e\) among
\(t+1,\ldots,t+n\). Each active level contributes at most
\(np^{-e}+1\), and there are at most \(C_0\log_p(n+1)\)
levels. Thus \(\gamma_t(n)\le n/(p-1)+O(\log(n+1))\)
uniformly on the set. A ball of radius \(p^{-m}\) centred
there enforces agreement on all blocks with
\(b_k\le(p-1)m-O(\log(m+1))\), with the implicit constant
chosen large enough for the first-difference test.
The measure of the restriction to this positive-measure set
therefore assigns each radius-\(p^{-m}\) ball mass
\(O((m+1)^C p^{-dm})\). A ball meeting the set can be
centred there without changing the ball. This bound proves
the lower dimension \(d\) and shows that at least
\(\gg (m+1)^{-C}p^{dm}\) residue classes modulo \(p^m\)
meet \(C_p\).

For a positive start \(a\), divisibility of an anchored sum
after division by \(a!\) through \(p^m\) depends only on
\(a\bmod p^m\). Indeed all \(W_j(a)\) are integer polynomials,
and \(J_p(j)\ge m\) makes \(W_j(a)\) zero at this precision
for every \(a\). Every class meeting \(C_p-2\) has such an
anchored finite sum for every positive representative.
Divide the occupied class bound by \(p^m\) to obtain
the excess-density bound in the theorem. This counts
positive starts with a finite-precision witness, leaving
positive integer membership in \(\mathcal Z_p\) open.

### The two block choices

For \(p=2\), use \(L=16\) and
\[
                     D_k=15+\min(v_2(k+1),5).
\]
Take \(M\ge20\). The first fifteen factors between block starts
have valuation sum \(J_2(15)=11\), since \(t+b_k\) is
divisible by \(16\). The last factor is \(t+16(k+1)\)
with valuation at least \(4+\min(v_2(k+1),5)\).
This proves the required product divisibility.
Each block has at least the two distinct residues zero
and one, and its alphabet depends only on \(b_k\bmod 2^{D_k}\).
The alphabet sizes and the precisions are periodic with
\(P_0=2^{16}\).

The [exact alphabet checker](verify_erdos404_dyadic_dimension.py)
computes the product \(P_2=\prod_{k<P_0}|A_k|\).
For a cap index \(j=\min(v_2(k+1),5)\), the normalized
weights depend only on \(k\bmod 2^{11+j}\), and the
condition on \(v_2(k+1)\) is periodic with this modulus.
Enumerate each allowed class once and give it multiplicity
\(2^{5-j}\) in the full period. The cyclic support update
adds the translate by each successive weight to the current
bit set, starting from the empty sum. The resulting
\(\lfloor\log_2 P_2\rfloor=960483\) is obtained by taking
the bit length of the exact product and subtracting one.
In (404.6), \(d_2=\log_2 P_2/(16\cdot 2^{16})>0.91598\).
The test uses the same exact bit set update as the
five-adic local checker. From the repository root run

~~~sh
python3 research/factorial/verify_erdos404_dyadic_dimension.py
~~~

For \(p=41\), take \(L=p^8\), \(D_k=D=J_p(L-p)\), and
\(M\ge D\). In \([b_k,b_k+L-p]\) no interior boundary
has a carry of size eight. At each size-\(r\) carry,
\(1\le r\le7\), the width \(p\) or \(2p\) used in the
[prime-41 covering proof](#prime-41-and-equal-low-residues)
fits entirely on or after the first index \(b_k\):
the first possible relative boundary of size \(r\)
is \(p^r\). With zero initial precision at \(b_k\),
the block induction proves that sums from this block
cover modulo \(p^{J_p(b_k+L-p)-J_p(b_k)}=p^D\).
The equality follows from the relative index valuation
identity below \(p^8\). Also \(D<J_p(L)\), so
our required next-block divisibility holds for the
nearby \(t\) as well. Thus all the alphabets have size
\(p^D\). Since \((p-1)D=p^8-8p+7\), the dimension formula
gives \(d_{41}=1-321/41^8\).

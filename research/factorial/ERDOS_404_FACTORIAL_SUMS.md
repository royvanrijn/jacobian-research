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

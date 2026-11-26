# Erdős #404: factorial sums and unbounded 5-adic cancellation

For an integer \(a\ge1\) and a prime \(p\), let \(f(a,p)\) be the supremum of \(v_p(\sum_{i\in I}i!)\) over finite sets of distinct positive indices \(I\) with \(\min I=a\); allow the value \(\infty\). [Erdős #404](https://www.erdosproblems.com/404) also asks for an infinite increasing sequence of indices whose partial factorial sums have valuations tending to infinity. The concrete open target here is
\[
                         f(1,5)=\infty. \tag{404.1}
\]
Reported finite witnesses reach valuation at least 20,000, but they have not been imported or replayed in this repository. They do not prove (404.1).

## Uniform lifting obligation

For each precision \(j\), track the reachable residues of sums of **distinct** factorials modulo \(5^j\), with index 1 retained. Find an invariant that (i) permits cancellation of the target residue, (ii) survives the passage from \(j\) to \(j+1\), and (iii) accounts for every newly exposed carry, including at indices around \(5^j\). A finite automaton constructed after fixing a maximum precision is not a uniform invariant. The [existing digit-orbit](image-rigidity/MULTIVARIATE_FACTORIAL_DIGIT_ORBITS.md) and [5-adic automaton](image-rigidity/MULTIVARIATE_PADIC_DIGIT_AUTOMATON.md) notes concern factorial coefficients in polynomial rigidity; their finite residue mechanisms are method analogies, not a proof of (404.1).

If (404.1) is proved, compactness of \(\{0,1\}^{\mathbb N}\) gives an infinite set of indices containing 1 whose 5-adic factorial sum is zero: choose finite cancelling sets at unbounded precisions, take a coordinatewise convergent subsequence, and use that \(v_5(i!)\to\infty\). The limit set cannot be finite, because a finite sum of positive factorials is nonzero. Its increasing enumeration has partial-sum valuations tending to infinity. This implication does not supply the missing cancellations.

## Small completed special case

For every \(a\ge1\) with \(a\equiv6\pmod9\), one has
\[
                         f(a,3)=v_3(a!)+1. \tag{404.2}
\]
After division by \(a!\), the terms at \(a+1,a+2\) are \(7,2\pmod9\); every later term is \(0\pmod9\). The four subset residues are \(1,8,3,1\pmod9\), so valuation one is attained and valuation two is impossible. This is an elementary internal observation; novelty and external review have not been assessed.

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

The table and the 24 individual boundaries give every residue
modulo \(5^{19524}\) using distinct indices \(2\) through \(78119\),
by taking \(N=78120\). This bounded full residue certificate
differs from a single witness for the residue \(-1\).
The next boundary is \(N=78125\), with \(v_5(N)=7\).
A uniform block of fixed width \(W\) cannot work for all \(r\):
once \(5^r>W\), it has at most \(2^W\) subsets and the target has
at least \(5^r\) residues. Thus the covering route needs growing
lookback or a condition that tracks a narrower set of entrance
carries.

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

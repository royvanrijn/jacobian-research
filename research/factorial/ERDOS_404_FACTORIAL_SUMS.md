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

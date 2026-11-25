# Conditional prime-support dictionaries

One can charge the primes of a root by a sequence of conditional
weights. The only restriction is a bound on the sum of the weights
for the next prime, for each already recorded prefix. The resulting
root dictionary bound is uniform over the choice of ordering and
the conditioning rule. At saturation this also constrains
conditional prime incidences, including correlations that are
invisible to unconditioned prime frequencies.

## Coding an ordered prime support

Work in the primitive anatomic box setup of (EB2), with the root
coordinates indexed by their letters and degrees. Use a partition
of the root positions into dictionary blocks \(G\).
One may always take singleton blocks. Alternatively, at positive
numerical heights, each block may be one of the bounded groups of
high-degree roots allowed in the grouped dictionary theorem.
The number of blocks and the truncation are fixed at each budget
before the height limit. The number of positions in one block is at
most a fixed \(m\).
[The coordinate dictionary bound (EB2)](ABC_EFFECTIVE_RADICAL_BUDGET.md#the-one-third-coefficient-for-coordinate-dictionaries)
· [The grouped dictionary bound](ABC_GROUPED_ROOT_DICTIONARY_ENTROPY.md#group-dictionary-count)

A representation of the support of a block is an ordered list
\(s=(z_1,\ldots,z_k)\) of its labelled primes, where \(z_j=(i_j,p_j)\)
records the root position and \(p_j\mid x_{i_j}\). A prime appears
once in its own root. The coordinate coprimality ensures that it
appears in no other root. Exponents are recovered separately.
Fix, for each block and each possible prefix \(P\), next-symbol
weights \(u_G(z\mid P)\ge0\) with
\(\sum_z u_G(z\mid P)\le1\). These weights are allowed to depend
on \(X\). For any representation put
\[
 \begin{split}
 c_G(s)&=\sum_{j=1}^k
       \min\{\log_X p_j,-\log_X u_G(z_j\mid z_1,\ldots,z_{j-1})\},\\
 \Delta_G(s)&=\sum_{j=1}^k
       \bigl(\log_X(p_j u_G(z_j\mid z_1,\ldots,z_{j-1}))\bigr)_+ .
 \end{split}                                                    \tag{CP1}
\]
At zero weight the first summand inside the minimum applies and
the positive discount term is zero. Thus the charge \(c_G\) equals
the logarithmic radical weight of the block in base \(X\) minus
\(\Delta_G\).

For any \(b\ge0\), the block vectors admitting at least one
representation \(s\) with \(c_G(s)\le b\) form a dictionary of
cardinality at most \(X^{b+o_X(1)}\). To see this, give a
transition to \(z=(i,p)\) weight
\(q_G(z\mid P)=\max(p^{-1},u_G(z\mid P))\). From any prefix the sum
of the possible next weights is at most
\[
                  1+m\sum_{p\le X}p^{-1}=O_m(\log\log X).
\]
Iterating this estimate, the sum of the product weights of all
lists of length \(k\) is at most \(O_m(\log\log X)^k\).
Each list with charge at most \(b\) contributes a product
weight at least \(X^{-b}\).
Since the block product is at most \(X^{O(m)}\), it has
\(O_m(\log X/\log\log X)\) distinct prime factors, by comparison
with a factorial. Summing over these lengths gives an
\(X^{o_X(1)}\) factor. The exponents on any one labelled
support likewise cost \(X^{o_X(1)}\) by stars and bars, since
their sum is \(O_m(\log X)\) and their number is
\(O_m(\log X/\log\log X)\).
Both errors are uniform in the transition weights and need no
probability model for the supports.
[Uniform fixed-support enumeration](ABC_EFFECTIVE_RADICAL_BUDGET.md#counting-by-a-sublogarithmic-prime-partition-sum)

Now suppose each row admits representations with total charge
\(\sum_G c_G\le b_0\). Bin the individual charges with total
rounding loss less than any fixed positive multiple of
\(\lambda\). The number of bins is bounded at fixed truncation.
In each use the smaller of the charge cap and the raw
root-height cap for each dictionary exponent. Applying (EB2)
or the grouped theorem gives, for every fixed \(\eta>0\) and
all sufficiently small fixed budgets,
\[
            N(\text{rows admitting total charge }\le b_0)
                     \ \ll X^{b_0/3+\eta\lambda+\varepsilon}.
                                                                  \tag{CP2}
\]
The numerical root-height caps were used unaltered. Since the
true radical remains bounded by \(X^\lambda\), a discount of
at least \(\zeta\lambda\) gives \(b_0\le(1-\zeta)\lambda\), and
thus a strict saving for \(\eta<\zeta/3\).

An example of a prescribed, conditional narrow support is a
set of pairs of primes \(p<q\) with
\(q-p\le p^{1-\nu}\), for some fixed \(0<\nu<1\).
Suppose the pairs are disjoint within each row, and the two
factors of every pair are in the same allowed dictionary
block. Encode the paired factors first, each \(p\) immediately
followed by its \(q\), and then the remaining factors. On even
steps give a candidate \(q\) uniform weight over
all labelled primes in the prescribed interval after the
preceding \(p\). On odd steps use natural prime costs.
The \(q\)-cost on a paired step is at most
\((1-\nu)\log_X p+o_X(1)\), whereas its natural cost is
\(\log_X q\ge\log_X p\). More precisely the error is
\(O_m(1/\log X)\) per step, for a total of \(o_X(1)\).
Thus if the smaller factors of the pairs carry at least
\(\zeta\lambda\) logarithmic radical weight, (CP2) applies
with \(b_0=(1-\nu\zeta)\lambda+o_X(1)\).
The pairs may vary with the row and may use primes outside
every fixed thin prime pool. This example only asserts a
count on the specified sector.

## Conditional occurrence at saturation

For a uniform family \(\mathcal F_X\) in a box, choose one
ordering of each support on each row. The choice may depend
on the row. Once the choice is made, let
\(\pi_G(P)\) be the probability that the list in block \(G\)
has prefix \(P\), and define \(u_G(z\mid P)\) to be the
conditional probability of the next symbol \(z\). A list
ending at \(P\) supplies the remaining probability mass.
Zero-probability prefixes can be given all next weights zero.

Suppose the family saturates the one-third coefficient, so
its count exponent divided by \(\lambda\) tends to \(1/3\)
in the height-first limiting sense. For every fixed
\(\zeta>0\), (CP2) with allowance smaller than \(\zeta/3\)
makes the set of rows with \(\sum_G\Delta_G\ge\zeta\lambda\)
power-smaller for sufficiently late fixed budgets. The
discount is nonnegative and at most \(\lambda\) on every row.
Taking means, this gives
\[
 \mathbb E_{\mathcal F_X}
    \sum_G\sum_{j\le k_G}
       \bigl(\log(p_j u_G(z_j\mid P_j))\bigr)_+
                       =o(\lambda\log X),                     \tag{CP3}
\]
where \(P_j=(z_1,\ldots,z_{j-1})\).
The argument works with frequencies recomputed on any
subfamily that also saturates.

There is a matching bound on negative deviations. Write
\(Z(P)=\sum_z u_G(z\mid P)\) and
\(J_G=|G|\sum_{p\le X}1/p\). For positive \(Z(P)\) the log-sum
inequality on the alphabet of labelled primes gives
\[
 \sum_z u_G(z\mid P)
        \log\frac{1}{p_z u_G(z\mid P)}
       \le Z(P)\log(J_G/Z(P))
       \le Z(P)\log\max\{1,J_G\}+1-Z(P).
\]
It holds also at \(Z(P)=0\) by omitting zero-weight terms.
Sum with coefficients \(\pi_G(P)\) over the prefixes. The sum
of \(\pi_G(P)Z(P)\) is the expected length of the list.
The sum of \(\pi_G(P)(1-Z(P))\) is one, since every list
terminates at a unique prefix. Summing also over the blocks
therefore bounds the expected negative deviation minus
positive deviation by
\[
       O_m(\mathbb E\omega(abc)\log\log\log X)+\#\{G\}
                                                   =o_X(\log X).
\]
Here the fixed number of blocks, the factorial bound on
\(\omega\), and the height limit at each fixed \(\lambda\)
justify the error before dividing by \(\lambda\log X\).
Combining with (CP3) proves
\[
 \mathbb E_{\mathcal F_X}
    \sum_G\sum_{j\le k_G}
        \left|\log(p_j u_G(z_j\mid P_j))\right|
                          =o(\lambda\log X).                  \tag{CP4}
\]
For each fixed \(\mu>0\), the expected logarithmic radical
mass of next symbols whose conditional probabilities fall
outside \([p^{-1-\mu},p^{-1+\mu}]\) is thus \(o(\lambda)\) in
base \(X\). All these statements are uniform in the choice
of one ordering per row and block. The grouping rule limits
which earlier primes can be used as context.

## Conditional presence along a support forest

One can also use conditional presence probabilities, with
no choice of next-symbol distribution from the sample.
Write \(f(z)\) for the probability that the labelled prime
\(z\) occurs in its block and \(f(z\mid y)\) for the probability
that \(z\) occurs there conditional on the occurrence of \(y\).
Give zero weight to the latter term when \(f(y)=0\).
Let \(K_X\ge1\) be a common bound on the number of factors
of \(abc\) on a row, with \(K_X=O(\log X/\log\log X)\).
Then \(\sum_z f(z)\le K_X\) and
\(\sum_z f(z\mid y)\le K_X\) for each parent symbol.

For any rooted forest on the support in a block, order it
with each parent preceding its children. Use the prefix policy
\[
 u_G(z\mid P)=\frac{f(z)}{2K_X}
           +\frac{\sum_{y\in P}f(z\mid y)}
                  {2K_X\max(1,|P|)} .
\]
The two sums each contribute at most \(1/2\) to the total
next weight. A root of the forest may then be charged by
its marginal probability and a child by its probability
conditional on its designated parent, up to an additional
\(O(\log\log X)\) in the unnormalized charge per factor.
The parent designation need not be recoverable from the
prefix alone because the common policy sums over every
available parent.

At saturation, the many-prime sector of (EB4) permits
discarding rows on which
\(\omega(abc)\log\log X\) is at least any fixed positive
multiple of \(\lambda\log X\). The accumulated forest
overhead is therefore negligible in the iterated limit.
Equations (CP1)–(CP2) imply that no forest, even one varying
by row, can secure a fixed positive total discount in
units of \(\lambda\) on more than a power-smaller fraction.
Here the discount on a child is
\((\log_X(p_z f(z\mid y)))_+\), and for a root it uses
\(f(z)\). In particular many correlated children can use
the same parent. This covers concentration patterns
that an argument using disjoint pairs alone could miss.
[Many-prime saving](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)

The conditions concern incidences inside a permitted
dictionary block. They leave open whether the additive
relation forces a positive conditional discount in
every saturated family.


## Conditioning across free letters under a cheap opposite divisor

The coding proof also applies to a block that meets two letters.
Take the groups on the two free letters \(U,V\) allowed in
(CG1), including its bounded-window groups. Their sizes can vary
with the budget in the bounded-window case. They are fixed before
the inner height limit, so the tree estimate still costs
\(X^{o_X(1)}\) for each group. If a prescribed full divisor on
\(T\) has record cost \(f\), a total charge \(b\) on the free groups
then gives the bound \(X^{f+b/2+o(\lambda)+\varepsilon}\).
The dictionaries from the charge bins are upper bounds also
after taking the fixed-divisor slices.
[The two-free-letter dictionary inequality (CG1)](ABC_CROSS_LETTER_ROOT_DICTIONARY_ENTROPY.md#one-divisor-and-two-free-letters)

In a saturated balanced box with arbitrarily cheap positive-weight
full moduli on \(T\), the numerical root caps on \(U,V\) sum to
\((2/3+o(1))\lambda\). If the conditional discount on the free groups
is at least \(\zeta\lambda\) for some fixed \(\zeta>0\), the
charge dips by that amount. Fix an opposite modulus at cost
less than \(\zeta\lambda/4\) and take the coding allowances small
relative to \(\zeta\lambda\). The resulting bound is strictly below
the saturation exponent. Equations (CP3)–(CP4) therefore apply to
these groups, with earlier primes from the other free letter
available as context. The conditional-presence forest argument
can use this same partition.

With two arbitrarily cheap full locks one can charge dictionary
blocks meeting all three letters, by projecting each bounded-price
dictionary to the unfrozen positions after fixing the locks.
Projection does not increase its cardinality. In this case (CG2)
still has coefficient \(1/2\), so this implication by itself
excludes discounts above \((1/3+\zeta)\lambda\), for any fixed
\(\zeta>0\), from saturation. Smaller discounts that require all
three letters simultaneously need an additional count.
[The simultaneous full-divisor bound (CG2)](ABC_CROSS_LETTER_ROOT_DICTIONARY_ENTROPY.md#two-simultaneous-divisors)

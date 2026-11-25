# The diagonal moment limit for three residue dictionaries

Three simultaneous unit congruences, diagonal bounds for every
bounded-order moment of every union of groups, and near-independence
of bounded group marginals on each pair of letters can coexist with
a count equal to the square root of the product of the dictionary
sizes. This example concerns abstract residue maps. Those maps are
assigned to the labels and need not be the reductions of bounded
integer roots.

## A residue model attaining the square-root count

Fix integers \(k\ge3\) and \(H\ge1\). Let \(q>2H\) be prime and
put \(P=\binom{k}{2}\). Take one group for each vertex of the complete
graph on \(k\) vertices. At each incidence \((g,e)\) of a vertex and
an edge choose labels
\[
        x_{g,e,A},x_{g,e,B},x_{g,e,C}\in\mathbb F_q,
        \qquad x_{g,e,C}=x_{g,e,A}+x_{g,e,B}.
\]
The choices at different incidences in a group are independent in
its dictionary. On letter \(W\), the group's label is the vector
\((x_{g,e,W})_{e\ni g}\). Each letter label has \(q^{k-1}\) possible
values. Every pair of letter labels in a group is uniform and
independent in its dictionary. The joint dictionary has \(q^{2(k-1)}\)
entries.

Use the integer encoding
\[
  \beta(t)=\sum_{r=1}^{2H}[t^r]_q B_0^{r-1},\qquad B_0=2Hq+1 ,
\]
where \([\,\cdot\,]_q\) is the least nonnegative representative.
For \(1\le j\le2H\), equality of two sums of \(j\) encoded labels
determines equality of their multisets. Indeed there are no carries,
and Newton's identities in \(\mathbb F_q\) recover the multiset from
its first \(j\) power sums.

Give every edge a fixed orientation. In the modulus belonging to
letter \(T\), choose one radix position for each of the \(2P\) pairs
\((e,W)\) with \(W\ne T\). Insert \(\beta(x_{g,e,W})\) at this
position for the group at the initial end of \(e\), and its negative
at the terminal end. Give the two visible letters opposite signs,
as in a ratio of the other two letters modulo \(T\). Combine the
digits in radix \(R>4H B_0^{2H}+1\). Pass to the cyclic group of
order \(M_T\), where \(M_T>R^{2P}\). Each component is the difference
of two maps on the visible letter labels. An equality of products,
in multiplicative notation, is an equality of these signed sums
modulo \(M_T\). For all the moment comparisons below the integer
difference has absolute value less than \(R^{2P}\). Each signed
digit also has absolute value less than \(R/2\). Thus the modular
equality forces every signed digit to vanish by successive reduction
modulo \(R\).

Select one entry at each vertex. The three product constraints into
the three cyclic groups, at identity targets, now hold exactly when
the two ends of each edge have the same labels on every letter.
In fact, each of a letter's edge positions occurs in the two
modulus components for the other letters, with its own two endpoints
as the only contributors. Thus solving all three constraints means
choosing one Latin triple on every edge, independently. There are
\(q^{2P}\) solutions, whereas the product of all dictionary sizes
is \(q^{4P}\).

The moduli can be interpreted as actual finite unit groups. Take
three distinct fixed odd primes \(\ell_T\) and choose a generator
of the units of \(\ell_T^{E_T}\). Interpret each encoded integer
as an exponent of this generator, with
\(M_T=\varphi(\ell_T^{E_T})>R^{2P}\). At fixed \(k,H\), one may
take \(q\asymp X^{\lambda/(6P)}\) and \(R\asymp_H B_0^{2H}\). The
total raw logarithmic label cost on each letter is then
\((\lambda/3+o_X(1))\log X\), and
\(\log_X(R^{2P})=(2H/3)\lambda+o_X(1)\).
The prime-power modulus can therefore have any prescribed positive
numerical height exponent for all sufficiently small fixed \(\lambda\).
It has only one fixed prime in its radical. The assertions also hold
after reducing to prime-power divisors whose unit group orders exceed
\(R^{2P}\), since the generator maps to a generator under reduction.
This step does not identify these modular maps with residues of small
integers.

## Moments on arbitrary unions and local entropy

Take any union \(J\) of the groups and any \(1\le h\le H\). Its
\(2h\)-th signed moment counts two ordered lists of \(h\) independent
entries from each group in \(J\) with the same signature product.
The equality at each radix position forces matching label
multisets. If an edge has one endpoint in \(J\), its two lists of
\(h\) Latin triples have at most \((h!)^3q^{2h}\) possibilities.
Choose the first list, then order the three fixed multisets. If
both endpoints are in \(J\), exchange the two lists at the
negatively oriented endpoint. The two sides then each have
\(2h\) Latin triples and there are at most
\(((2h)!)^3q^{4h}\) possibilities. These constraints factor over
edges. Consequently the moment is at most
\[
        ((2h)!)^{3P}\left(\prod_{g\in J}|\mathcal D_g|\right)^h .
\]
This includes unions whose size grows with \(k\). Both \(H\) and \(k\)
may be chosen before the height limit. The small-budget threshold
must allow the requisite ambient modulus size, so the example
does not assert uniformity for unrestricted \(H\lambda\).

In the uniform solution family, the Latin triples belonging to
distinct edges remain independent. Consider any partition of the
original label positions \((g,W)\) on two letters into blocks
of at most \(b\) positions. The entropy of a block is the sum of its
contributions from the different edges. On each edge the two
distinct letters are independent, so the only losses from the
raw label costs arise from including the same letter at both
endpoints within a block. Each such inclusion costs \(\log q\),
and there are at most
\(\sum_{\text{blocks }I}\binom{|I|}{2}\le k(b-1)\) of them across
the partition. Thus the fraction of raw pair cost lost is at most
\((b-1)/(2(k-1))\), which tends to zero for fixed \(b\) as
\(k\) increases. For partitions on all three letters, excess
dependence between distinct original groups is also
\(O(bk\log q)\). An edge makes no such contribution unless both
its endpoints occur in the block, and then costs at most
\(3\log q\). The within-group Latin dependence still reduces its
joint cost from three letter costs to two. The original positions
can also be given evenly spaced scaled degree labels on a fixed
compact interval. A group in a bounded number of windows of width
\(\rho\) then has size \(O(1+k\rho)\). The same bound on the loss
fraction tends to zero when \(k\) tends to infinity and \(\rho\)
tends to zero. The assigned degree labels only index the abstract
positions.

The same comparison allows a record \(\mathcal R\) consisting of
any prescribed subset of the original label positions. Conditional
on the record, the edges are still independent. At edge \(e\),
let \(r_e\log q\) be the entropy of the recorded letters, with
\(r_e\in\{0,1,2\}\). The remaining labels of any chosen pair of
letters at either endpoint recover the full edge once the omitted
recorded entries are supplied. Their conditional entropy is
\((2-r_e)\log q\). For any partition \(\mathcal P\) of the remaining
positions of the pair into blocks of size at most \(b\), summing
these endpoint entropies costs at most \(O(bk\log q)\) in duplication
between endpoints within blocks. Adding the record gives
\[
 H(\mathcal R)+\tfrac12\sum_{I\in\mathcal P}H(X_I\mid\mathcal R)
        \ge 2P\log q-O(bk\log q),
\]
where \(X_I\) lists the labels in \(I\). The error bound is uniform
over the prescribed subset and block partition.

The example leaves out the signed additive equation on the
complete integer letters and the requirement that the signatures
are the monomial residues of uniformly small roots. It gives no
counterexample to an arithmetic saving. An argument using only
the displayed moment upper bounds and bounded-marginal entropy
saturation cannot improve the square-root coefficient on the
total joint dictionary cost.

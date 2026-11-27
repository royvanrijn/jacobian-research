# Exact sizes and removable factors in the odd polynomial Hurwitz domains

Fix \(r\ge3\), \(N=2r+1\), and an \(N\)-cycle \(P\). Let \(X_D\) consist
of the ordered \(r\)-tuples of double transpositions with product \(P\).
Let \(X_T\) consist of the corresponding tuples with a 3-cycle in the
first position and double transpositions elsewhere. Write \(d=(N-3)/2\).
The absolute domains have sizes
\[
 |X_D/\langle P\rangle|=N^{r-2}d^r,\qquad
 |X_T/\langle P\rangle|=N^{r-2}d^{r-1}.
 \tag{1}
\]
In the pure domain, the exact fraction of absolute tuples whose last
factor has one exclusive leaf on each transposition is
\[
 \left(1-\frac{2}{N}\right)^{r-2}
 \left(1-\frac{2}{N-3}\right)^{r-1}
 \ \longrightarrow\ e^{-2}.                                          \tag{2}
\]
By an exclusive leaf we mean a sheet moved by this transposition and
fixed by every earlier factor. For \(r=3\), one can define the same
condition on the degree-five residual domain without a monodromy
restriction. The proof below also gives the labeled counts in degree
five when \(r=2\).

For a uniform pure tuple, the number of directly removable
positions has an exact joint-count formula. It has mean
asymptotic to \(e^{-2}r\), is asymptotically normal after centering
and scaling. In each of the pure and marked mixed domains, the
fraction with no removable double transposition is \(O(e^{-cN})\)
for an absolute \(c>0\).

Equation (1) agrees with the Goulden–Jackson cactus count. The derivation
below verifies the cardinalities directly on the fixed-product spaces,
including the free quotient when \(N\) is composite. Equation (2)
gives the size of the directly removable locus in
\(\texttt{GAP-ARITH-HURWITZ-UNIFORM-REMOVABLE-DT-REPRESENTATIVE}\).
The terminal-mass extraction theorem proves that every residual pure-DT orbit
meets it. The later [mixed isolation and connectivity theorem](MIXED_HURWITZ_ISOLATION_AND_CONNECTIVITY.md#theorem-m7-all-degree-mixed-connectivity-extraction-and-restoration) proves the corresponding mixed extraction assertion. [Goulden–Jackson, Theorems 3.1–3.2 and the last line of the proof](https://www-users.cse.umn.edu/~reiner/Classes/Math8680Fall2014Papers/GouldenJacksonOnCacti.pdf)

## The class product as a finite difference

Set \(L=N-1\) and write \(\chi_j\) for the character of the \(j\)-th
exterior power of the standard \(L\)-dimensional representation of
\(S_N\), for \(0\le j\le L\). The generating polynomial on a
permutation \(g\) is
\(\prod_{c}(1-(-s)^{|c|})/(1+s)\), with \(c\) ranging over its cycles.
These characters are distinct irreducibles. Here is an elementary
inner-product check valid in every degree. Counting permutations with
given cycle multiplicities gives the joint cycle index
\[
 \exp\left(\sum_{h\ge1}\frac{z^h}{h}
                  (1-(-s)^h)(1-(-t)^h)\right)
 =1+\frac{(1+s)(1+t)z}{(1-z)(1-stz)} .
\]
Indeed a cycle count vector of mass \(N\) has probability
\(\prod_h(h^{m_h}m_h!)^{-1}\), and expanding the exponential sums these
vectors. For \(N\ge1\), extract \(z^N\) and divide by
\((1+s)(1+t)\) to obtain \(\sum_{j=0}^{N-1}(st)^j\).
This says the exterior-power characters are orthonormal. At \(P\),
\(\chi_j(P)=(-1)^j\). The \(N\) squares of these values exhaust
the column norm \(|C_{S_N}(P)|=N\), so every other irreducible
character vanishes at \(P\).

For classes \(C_i\), the central class sum acts on an irreducible
of dimension \(d_j\) by \(|C_i|\chi_j(C_i)/d_j\). Take the trace of
\(P^{-1}\) times the product of these sums in the regular
representation and divide by \(N!\) to obtain the coefficient of \(P\).
As \(d_j=\binom Lj\), the exact total is
\[
 \frac{\prod_i|C_i|}{N!}\sum_{j=0}^L(-1)^j\binom Lj
                 \prod_{i=1}^r\frac{\chi_j(C_i)}{\binom Lj}.           \tag{3}
\]
For a double transposition \(D\) the generating polynomial is
\((1+s)^L-4s(1+s)^{L-2}\). For a 3-cycle \(T\) it is
\((1+s)^L-3s(1+s)^{L-2}\). Therefore
\[
 \frac{\chi_j(C_i)}{\binom Lj}
 =1-\xi_i\frac{j(L-j)}{L(L-1)},\qquad
 \xi_i=4\ \text{for }D,\qquad \xi_i=3\ \text{for }T.
\]
The product in (3) has degree \(2r=L\) in \(j\).
The alternating binomial sum is the \(L\)-th finite difference,
since \(L\) is even. It is \(L!\) times the leading coefficient.
Using \(|D|=3\binom N4\), \(|T|=2\binom N3\) gives
\[
 \#\{g_1\cdots g_r=P:g_i\in C_i\}
  =\frac1N\prod_i\frac{|C_i|\xi_i}{L(L-1)}
  =N^{r-1}d^{\,\#\{i:C_i=D\}}.                                       \tag{4}
\]
This argument works for any placement of the 3-cycles. It uses ordinary
finite-group character orthogonality, also used in the repository's
[degree-thirteen census](DEGREE_FOURTEEN_HURWITZ_ORBIT.md#an-independent-exact-total-bound).

Every nonidentity power of \(P\) has orbits all of one odd length greater
than one. It cannot preserve a four-element DT support, since that length
does not divide four. Consequently the simultaneous conjugation quotient
of each of the two stated domains is free. Dividing (4) by \(N\) gives (1).

## Counting the directly removable tuples

Compose permutations right to left. Suppose the last factor has
exclusive leaves \(a,b\), let \(R\) be the product of the earlier
factors, and write \(P=R g_r\). Since \(R\) fixes \(a,b\),
the two transpositions must pair \(a\) with \(P^{-1}(a)\) and \(b\)
with \(P^{-1}(b)\). Thus the two leaves are nonadjacent in the cyclic
order of \(P\). Conversely for each nonadjacent unordered pair, these
transpositions define \(g_r\), and \(P g_r\) fixes \(a,b\) and cycles
through all remaining labels in their inherited cyclic order. Indeed,
for a predecessor of either leaf the product skips that leaf and
continues at its successor, and otherwise it follows \(P\). Each
factorization of this residual long cycle by \(r-1\) double
transpositions lifts by fixing \(a,b\) in the earlier factors.

The leaves in each such tuple are determined uniquely. If both
endpoints of a transposition were fixed by every earlier factor,
\(R g_r\) would have that transposition as a two-cycle, contradicting
the long-cycle product. There are \(\binom N2-N=Nd\) choices of leaf
pair. The resulting labeled count is therefore
\[
 Nd (N-2)^{r-2}\left(\frac{N-5}{2}\right)^{r-1}.
\]
Dividing by (4) in the pure case gives (2). Dividing both counts by \(N\)
shows that the same fraction holds after quotienting. For \(N=15\)
the absolute totals are \(212\,576\,400\,000\) pure and
\(35\,429\,400\,000\) marked mixed. The pure count with directly
removable last factor is \(34\,808\,718\,750\). Direct enumeration\nof the labeled domains in degrees 5 and 7 gives pure totals 5 and 392,\nmarked totals 5 and 196, and directly removable counts 0 and 70,\nagreeing with the formulas.

## The joint removable-factor law

In the pure domain let \(I_i\) indicate that factor \(i\) has one
exclusive leaf in each transposition, where an exclusive leaf is
fixed by every other factor in the full tuple. Let \(R_r=\sum_i I_i\).
For a specified set \(J\) of \(k\) positions, put \(q=r-k\) and
\(m=N-2k=2q+1\). For \(0\le k\le r-2\), the exact labeled joint
count and probability are
\[
 \begin{aligned}
 \#\{\text{tuples}: I_i=1\text{ for every }i\in J\}
   &=N\,m^{r-2}q^k(q-1)^q,\\
 p_k:=\mathbb P(I_i=1\text{ for every }i\in J)
   &=\left(\frac{m}{N}\right)^{r-2}
     \left(\frac{q}{r-1}\right)^k
     \left(\frac{q-1}{r-1}\right)^q .
 \end{aligned}                                                        \tag{5}
\]
The joint count is zero for \(k\ge r-1\). The cases \(k=0,1\)
recover (4) and (2). Absolute counts are obtained by dividing the
labeled counts by \(N\).

To prove (5), call the positions in \(J\) marked and move their
factors to the end without changing the relative order of marked
factors. Each swap is
\[
 (g,h)\longmapsto (h,h^{-1}gh),
\]
with \(g\) marked and \(h\) unmarked. It preserves the product.
The unmarked factor fixes each marked exclusive leaf, so both the
swap and its inverse preserve the exclusive leaf of each marked
transposition. This gives a bijection between the joint locus at
the specified positions and the joint locus for a trailing block.

Now delete the two exclusive leaves in each marked factor,
beginning with the last factor and applying the predecessor rule
above. A leaf removed at a given step must have its predecessor in
the final core of \(m\) sheets, after every previous deletion.
Indeed its predecessor is its partner in that factor, and every
other marked leaf is exclusive to a different factor. Therefore
each run of marked leaves following a core sheet has strictly
increasing deletion-step indices, read in the cyclic order of
\(P\). Leaves deleted at the same step occupy two distinct runs.
Conversely these conditions allow deletion in the prescribed
order, since the next leaf in each selected run is adjacent to
its predecessor core sheet. Reinserting the paired transpositions
in reverse order recovers the full tuple. Each newly inserted
transposition joins its designated leaf to a core sheet and
leaves all other marked leaves fixed.

Choose one labeled core sheet as a start, then choose the two
core gaps to receive each of the \(k\) deletion-step indices.
Within each gap the leaves are ordered by their step indices.
The resulting cyclic word determines which original labels are
core sheets and which are leaves, because the long cycle and the
starting label are fixed. Each assignment of the leaf labels is
counted \(m\) times, once per choice of starting core sheet.
There are thus \(N\binom m2^k/m\) assignments. For \(q\ge2\), (4)
gives \(m^{q-1}(q-1)^q\) factorizations of the residual cycle.
Multiplying and using \(\binom m2=mq\) gives (5). For \(q=1\)
one DT cannot produce a 3-cycle. For \(q=0\) there is only one
core gap, so no step can choose two distinct gaps.

The exact falling-factorial moments and probability generating
polynomial follow by summing the joint events:
\[
 \mathbb E(R_r)_k=(r)_k p_k,\qquad
 \mathbb E s^{R_r}=\sum_{k=0}^{r-2}\binom rk p_k(s-1)^k .                \tag{6}
\]
Here \(p_0=1\) and \((r)_k=r(r-1)\cdots(r-k+1)\). A direct
meet-in-the-middle enumeration in degree nine used right-to-left
products and tested each factor against the union of the supports
of the other three. Each specified position was removable in
\(10\,584\) labeled tuples. All six specified pairs of positions
occurred jointly in \(900\) tuples each. The numbers of tuples
with zero, one, and two removable factors were \(22\,113\),
\(31\,536\), and \(5\,400\), with total \(59\,049\).

## Concentration of the removable-factor count

Put \(p=e^{-2}\) and \(v=p-3p^2>0\). In the uniform pure domain,
as \(r\to\infty\),
\[
 \frac{\mathbb E R_r}{r}\longrightarrow p,\qquad
 \frac{\operatorname{Var}R_r}{r}\longrightarrow v,\qquad
 \frac{R_r-\mathbb E R_r}{\sqrt{vr}}
             \ \xrightarrow{\ \mathrm{law}\ }\ \mathcal N(0,1).
                                                                        \tag{7}
\]
For every fixed \(A>0\), the fraction with \(R_r=0\) is
\(O_A(r^{-A})\). These statements hold on labeled and absolute
tuples since the quotient is free and preserves the indicators.

For fixed \(k\), expanding the three logarithms in (5) gives
\[
 \log p_k=-2k+\frac{5k/2-k^2}{r}+O_k(r^{-2}).
\]
Thus \(p_1\to p\) and \(p_2-p_1^2=-2p^2/r+O(r^{-2})\). Inserting
these into \(r p_1(1-p_1)+r(r-1)(p_2-p_1^2)\) gives the variance
in (7). The higher orders of the same expansion control all
fixed cumulants. More precisely, the coefficient of \(r^{-j}\)
in \(\log p_k\) is a polynomial in \(k\) of degree at most
\(j+1\). To see the degree bound, set \(x=1/r\). The three
logarithms have arguments
\[
 1-\frac{kx}{1+x/2},\qquad
 1-\frac{(k-1)x}{1-x},\qquad
 1-\frac{kx}{1-x},
\]
and their multipliers are \(x^{-1}-2,k,x^{-1}-k\),
respectively. Before multiplication the coefficient of \(x^h\)
has degree at most \(h\), which gives the assertion.

For fixed \(j\ge1\), let \(w_j=\Delta^j(\log p_\bullet)(0)\)
be the \(j\)-th forward difference at zero. For all sufficiently
large \(r\) the necessary \(p_i\) are positive. The degree bound
gives \(w_j=O_j(r^{1-j})\), since a \(j\)-th forward difference
annihilates lower degrees. Newton inversion yields
\(\log p_k=\sum_{j=1}^k\binom kjw_j\).

Take \(d\) distinct positions occurring in a cumulant of fixed
total order \(h\), where some positions may occur repeatedly.
Put \(z_i=e^{t_i}-1\) for the \(d\) positions, and give each
subset \(H\subseteq[1,d]\) of size at least two a variable
\(u_H\). Their joint moment generating function is obtained by
putting \(u_H=e^{w_{|H|}}-1\) in
\[
 F(\mathbf t,\mathbf u)=
 \sum_{B\subseteq[1,d]}\prod_{i\in B}p_1 z_i
                   \prod_{\substack{H\subseteq B\\|H|\ge2}}(1+u_H).
\]
Indeed the expansion of \(\prod_i(1+z_i I_i)\) gives the joint
probability \(p_{|B|}\), and Newton inversion factors that
probability as shown. A cumulant with multiplicities is the
corresponding mixed derivative of \(\log F\) at \(\mathbf t=0\).
Through total degree \(h\) in \(\mathbf t\), the coefficients
of \(\log F\) are polynomials in the \(u_H\), since \(F-1\)
has zero constant term in \(\mathbf t\).

For \(d\ge2\), consider a monomial in the \(u_H\) contributing
to the mixed derivative. Regard its subsets \(H\) as a
hypergraph on the \(d\) positions, allowing repeated edges.
Set every variable outside the monomial to zero. If this
hypergraph is disconnected, \(F\) factors over its components,
including isolated vertices, because the summations of \(B\)
inside different components are independent. Then \(\log F\)
is a sum over components, and its mixed derivative involving
all \(d\) positions vanishes. A surviving monomial therefore
has a connected hypergraph on all \(d\) positions. Its edge
sizes satisfy \(\sum_H(|H|-1)\ge d-1\) with multiplicity, since
adjoining an edge of size \(|H|\) reduces the number of connected
components by at most \(|H|-1\). After substitution its size is
\(O_h(r^{1-d})\) by the bound on \(w_j\). There are only
finitely many monomials at fixed \(h\), with bounded coefficients
since \(0\le p_1\le1\). For \(d=1\) the cumulant is bounded as a
cumulant of one Bernoulli variable.

By multilinearity, the cumulant of order \(h\) of \(R_r\) sums
over all \(h\)-tuples of positions. There are \(O_h(r^d)\)
with \(d\) distinct positions. The preceding bound shows that
each fixed cumulant of \(R_r\) is \(O_h(r)\). Consequently the
normalized cumulants in (7) of order at least three tend to
zero. The normalized moments converge to those of the normal
law, which is determined by its moments. For the exceptional
fraction, the moment of order \(2s\) about the mean is \(O_s(r^s)\).
Indeed only partitions with no singleton block contribute
in the moment-cumulant formula for a centered variable, so
there are at most \(s\) cumulant factors, each \(O_s(r)\).
Since \(\mathbb E R_r\asymp r\), Markov's inequality gives
\(\mathbb P(R_r=0)=O_s(r^{-s})\) for every fixed positive integer
\(s\), as claimed.

## The marked mixed domain

Write a mixed tuple as \((g_1,\ldots,g_r)\), where \(g_1\) is
a 3-cycle. Define \(I_i\) for \(i\ge2\) by the same exclusive-leaf
condition on the two transpositions of \(g_i\). Exclusivity tests
every other factor, including \(g_1\).

For a specified set of \(k\) double-transposition positions,
\(0\le k\le r-1\), put \(q=r-k\) and \(m=N-2k\). The stable
movement of marked factors to the end uses only the Hurwitz
moves \(Q_i\) with \(i\ge2\), which act on positions \(i,i+1\).
All marks start at positions at least two and move only to the
right, so \(g_1\) stays first. The same cyclic
core-gap construction gives \(N\binom m2^k/m\) leaf assignments.
The residual \(q\)-factor word begins with a 3-cycle and has
\(q-1\) double transpositions. Its labeled count is
\(m^{q-1}(q-1)^{q-1}\) by (4). At the endpoint \(q=1\), the
residual degree is three, and the only factor is the residual
3-cycle itself. Its count is one, expressed by \(0^0=1\) in the
formula. This endpoint requires only the labeled residual count.

Thus the exact mixed labeled count and joint probability are
\[
 \begin{aligned}
 \#\{\text{mixed tuples with all specified factors removable}\}
   &=N\,m^{r-2}q^k(q-1)^{q-1},\\
 p^{\rm mix}_k
   &=\left(\frac{m}{N}\right)^{r-2}
     \left(\frac{q}{r-1}\right)^k
     \left(\frac{q-1}{r-1}\right)^{q-1}.
 \end{aligned}                                                        \tag{8}
\]
In particular, when all \(r-1\) double transpositions are
required to be removable, the count is \(N3^{r-2}\) labeled and
\(3^{r-2}\) absolute. The generating polynomial is
\[
 \mathbb E s^{\sum_{i=2}^r I_i}=
     \sum_{k=0}^{r-1}\binom{r-1}{k}p^{\rm mix}_k(s-1)^k.                 \tag{9}
\]
For fixed \(k\) and large \(r\), \(\log p^{\rm mix}_k\) differs
from \(\log p_k\) by \(-\log(1-k/(r-1))\). This adds \(k/r\)
to the first correction and preserves the degree bounds on
all higher coefficients. The new first correction
\(7k/2-k^2\) again gives \(p^{\rm mix}_2-(p^{\rm mix}_1)^2
=-2e^{-4}/r+O(r^{-2})\). The connected cumulant argument therefore
gives the same mean and variance limits per \(r\), the same
normal limit with the mixed mean, and \(O_A(r^{-A})\) probability
that all these \(r-1\) indicators vanish.

Independent right-to-left enumeration in degree seven gave
\(196\) mixed tuples, \(70\) for either specified position and
\(21\) for both. Its histogram at zero, one and two removable
factors was \((77,98,21)\). In degree nine the total was
\(19\,683\), with \(5\,292\) for each of three specified positions,
\(900\) for each specified pair and \(81\) for all three. Its
histogram was \((6\,426,10\,719,2\,457,81)\).

## An exponential bound for the exceptional fraction

There is an absolute \(c>0\) such that the proportion with no
removable double transposition is \(O(e^{-cN})\) in each domain.
The argument uses an ordinary labeled tree and a uniform pairing
of its edges.

Expand each double transposition into its two disjoint
transposition edges. In the mixed case expand the 3-cycle into
an ordered pair of transpositions sharing one endpoint. The
expanded graph is connected since the factors generate a group
containing the long cycle, and there are \(2r=N-1\) edge
occurrences. Thus the edges are distinct and form a tree.
Conversely, the product of the edge transpositions of a tree
in any order is a long cycle. Delete a leaf edge \(t\), writing
the old product as \(U t V=(UV)(V^{-1}tV)\). Both \(U\) and
\(V\) fix the leaf. By induction \(UV\) is a long cycle on the
remaining sheets, and \(V^{-1}tV\) inserts the leaf into it.

Sample a uniform labeled tree on the \(N\) sheets, and then
uniformly pair its \(N-1\) edges into \(r\) ordered pairs of
positions, treating the two edges within a pair as unordered.
There are \(N^{N-2}(N-1)!/2^r\) samples. In the pure case,
accept exactly when the two edges in every pair are disjoint.
The accepted samples are in bijection with the pure tuples
over all long-cycle products. Multiplying the labeled count
(4) by \((N-1)!\) gives acceptance probability
\[
                \left(1-\frac{3}{N}\right)^r.                           \tag{10}
\]
For the mixed case, also order the two edges in the first pair
uniformly. Accept when this first pair shares an endpoint and
every later pair consists of disjoint edges. Each mixed tuple
has exactly three preimages, one for each ordered two-edge
presentation of its 3-cycle. On a prescribed three-element
support there are three unordered pairs of distinct
transpositions, and exactly one of the two orders gives the
specified orientation for each pair. The sample space has now
doubled, so (4) gives acceptance probability
\[
 \frac{3}{2(r-1)}\left(1-\frac{3}{N}\right)^r
       = \frac{3}{N}\left(1-\frac{3}{N}\right)^{r-1}.                    \tag{11}
\]
The acceptance probabilities in (10) and (11) are of order
one and \(1/N\), respectively. Both accepted models induce a
uniform tuple over all long cycles. Simultaneous relabeling
maps one fixed-product space bijectively to every other and
preserves removability, so the proportions are the same at a
prescribed \(P\).

A paired double transposition is removable precisely when
each of its two tree edges is incident to a leaf. Indeed a leaf
has degree one in the expanded tree if and only if its endpoint
occurs in just that transposition. Inside a double transposition
the edges are disjoint, and in the mixed case the first two
edges have union equal to the 3-cycle support. For \(N\ge3\)
distinct leaves have distinct incident edges. Write \(L\) for
the number of these leaf edges.

The Prüfer bijection encodes the uniform labeled tree by
\(N-2\) independent uniform sheet labels. At each step delete
the smallest remaining leaf and record its neighbor. Decoding
joins the smallest remaining label absent from the word to
its first entry and ends by joining the final two labels.
In particular each vertex occurs in the word one fewer time
than its degree. Hence \(L\) is the number of absent labels
and \(\mathbb E L=N(1-1/N)^{N-2}\sim N/e\). Changing one word
entry changes \(L\) by at most one. The Doob differences for
revealing the word entries successively thus have conditional
range of width at most one. For a centered variable \(Y\) in
any interval of width one, the second derivative of
\(\log\mathbb E e^{tY}\) is a variance under exponential
tilting and is at most \(1/4\). Integrating, applying this
conditionally to each Doob difference, and optimizing the
exponential Markov bound gives
\[
 \mathbb P(L\le \mathbb E L-u)
         \le \exp\!\left(-\frac{2u^2}{N-2}\right).
\]
Taking \(u=\mathbb E L-N/4\), for all sufficiently large \(N\)
the probability of \(L<N/4\) is at most \(e^{-N/72}\).

Conditional on the tree and the first paired positions, the
remaining edge matching is uniform. Revealing that first
pair removes at most two leaf edges. On \(L\ge N/4\) there
are at least \(L'\ge N/5\) remaining leaf edges for all
sufficiently large \(N\). Explore the partners of unpaired leaf
edges. Provided no two leaf edges have yet been matched,
at any of the first \(\lfloor L'/2\rfloor\) reveals the
probability of a leaf partner is at least \(L'/(2N)\ge1/10\).
Thus the probability of avoiding all leaf pairs among the
remaining positions is at most
\((9/10)^{\lfloor L'/2\rfloor}\le e^{-N/120}\) for all
sufficiently large \(N\). This bounds the unconditioned
probability of the event which, on acceptance, contains
every tuple without a removable factor. Dividing by (10)
or (11) proves the exponential bound, for instance with
\(c=1/240\) in either domain. It also holds after the free
absolute quotient.

An independent enumeration of all \(16\,807\) Prüfer words in
degree seven and all edge matchings gives \(282\,240\) accepted
pure samples and \(423\,360\) accepted ordered mixed presentations
over all long-cycle products. The respective counts with no
removable factor are \(131\,040\) and \(166\,320\). Dividing
by the number \(6!\) of long cycles and, for mixed, by the
presentation multiplicity three gives \(182\) and \(77\)
labeled exceptional tuples at fixed product.

## A sufficient size criterion for a braid orbit

Use the full braid action on the pure factor-type domain and the
distinguished-strand stabilizer action on the marked mixed domain.
An absolute orbit missing the removable locus lies entirely in
the exceptional set with no removable factors. In either
domain the exact absolute size of this set is the size in (1)
times the value at \(s=0\) of the corresponding generating
polynomial (6) or (9). Therefore an orbit of strictly larger
absolute cardinality meets the removable locus. For any fixed \(0<c'<c\), an orbit of relative size at least
\(e^{-c'N}\) meets the locus for all sufficiently large \(N\),
where \(c\) is the absolute exponent in the tree estimate.
This includes every orbit of relative size at least \(r^{-C}\)
for fixed \(C>0\). The separate structural
[Theorem 20](TERMINAL_MASS_HURWITZ_DESCENT.md#theorem-20-every-residual-pure-dt-orbit-admits-extraction)
now proves the every-orbit removable-factor assertion in the pure domain
for odd \(N\ge15\), without an orbit-size hypothesis. Its fixed-projection
argument excludes the last positive terminal-mass minimum. The
degree-thirteen base and two-leaf absorption give pure-DT connectivity in
all residual degrees. These conclusions add no new census or orbit-size
bound. The later [mixed isolation and connectivity theorem](MIXED_HURWITZ_ISOLATION_AND_CONNECTIVITY.md#theorem-m7-all-degree-mixed-connectivity-extraction-and-restoration) establishes marked mixed connectivity, and the
folding/attachment criterion then gives the even full-alternating passport.
[Marked subgroup and braid action](UNMARKED_POLYNOMIAL_HURWITZ_REDUCTION.md#the-spaces-and-the-exact-marked-subgroup)

## The even alternating Laurent passport

Let \(n=2b\ge10\). Fix an \((n-1)\)-cycle with one fixed sheet
as the product of an ordered \(b\)-tuple of double transpositions
generating \(A_n\). This is the full-alternating domain of the
even Laurent passport \(D^b,(n-1,1)\). Its product, graph, and
absolute orbit setting are those of
\(\texttt{ARITH-HURWITZ-UNICYCLE-FIXED-SHEET}\) and
\(\texttt{ARITH-ABSOLUTE-HURWITZ-DEGREE14-AND-UNIFORMITY}\).
In this domain the fraction without a directly removable double
transposition is \(O(e^{-cn})\), with \(c=1/240\) sufficient.
[Cycle folding, Theorem 1](ALTERNATING_HURWITZ_CYCLE_FOLDING.md#theorem-1-the-unicyclic-graph-and-fixed-sheet-trajectory)
· [Even-degree orbit and uniformity](ALTERNATING_HURWITZ_UNIFORMITY_AND_REAL_NODES.md)

For the proof, sample a uniform labeled tree on the \(n\) sheets,
adjoin an independent uniform unordered edge between distinct
sheets, and uniformly pair the resulting \(n\) edge occurrences
into an ordered list of \(b\) pairs. The added edge is
distinguished in this sampling space even if it duplicates a
tree edge. Accept when every paired block consists of disjoint
edges and the resulting tuple is in the stated full-alternating
passport, allowing all products of the prescribed cycle type.
The sample size is
\[
 S_n=n^{n-2}\binom n2\,\frac{n!}{2^b}.                                 \tag{12}
\]

For an accepted tuple the expanded graph is connected and has
one cycle, counting a parallel pair as a cycle of length two.
Its number of sampled preimages equals its cycle length
\(\ell\). Deleting a cycle-edge occurrence leaves a tree and
determines the distinct added occurrence and the ordered pairing.
When the cycle is parallel, the duplicated edge is in two
different blocks because each block has disjoint edges. Either
block can contain the distinguished added occurrence, giving
two preimages with the same underlying tree. Thus the accepted
sampling distribution weights each tuple by \(\ell\in[2,n]\).
A uniform-tuple fraction is at most \(n/2\) times the corresponding
fraction in the weighted distribution.

There are at least
\[
 D_n=n(n-2)!(n-5)(n-1)^{b-2}(b-2)^{b-1}                               \tag{13}
\]
accepted tuples before weighting. For each choice of a fixed sheet
\(z\) and long cycle on the other \(n-1\) sheets, take an odd pure
factorization counted in (4) with \(b-1\) factors. Its group is
\(A_{n-1}\). In fact an odd nontrivial block size is at least three,
and a double transposition moving four sheets cannot permute such
blocks nontrivially. The long-cycle product then gives primitivity.
Jones's Proposition 6.1(3) applies since \(n-1\ge9\).
Split the first factor \(uv\), with \(u,v\) disjoint transpositions,
into \((z\,a)u,(z\,a)v\) in a fixed edge order, where \(a\) is any
of the \(n-5\) old labels outside the support of \(uv\). Coalescing
these two factors recovers the old factorization and the two
doubled-edge blocks recover \(a\). This proves injectivity.
The new group contains \(A_{n-1}\) on the old labels and moves
\(z\). It is transitive and, by orbit–stabilizer and evenness, is
\(A_n\). Thus the acceptance probability is at least
\[
 D_n/S_n =
 \frac{4(n-5)(n-1)^{b-4}(n-4)^{b-1}}{n^{n-1}}
 \asymp n^{-3}.                                                       \tag{14}
\]
[Jones, arXiv:1304.5444v1, Proposition 6.1(3)](https://arxiv.org/pdf/1304.5444v1)
· [Pinned import EXT-PAPER-ALT-PRIORITY-JONES-BEAUVILLE-2013](alternating_hurwitz_priority_ledger.json)

The additional edge destroys at most two of the original tree
leaves. Use the surviving original leaf edges and the Prüfer
concentration and random-pairing estimate above. With
probability \(1-O(e^{-n/72})\) there are at least \(n/4-2\) such
edges, and conditioned on the tree and added edge the probability
of avoiding every paired pair of those edges is at most
\(e^{-n/120}\) for sufficiently large \(n\). On acceptance,
every paired pair of surviving leaf edges gives a removable
factor. Division by (14) and the additional factor \(n/2\) for
removing the cycle-length weighting still leave the bound
\(O(e^{-n/240})\). Conjugacy makes the fraction the same at a
fixed product. The quotient by its centralizer is free because
the tuples generate \(A_n\) and its centralizer in \(S_n\) is
trivial.

For \(n\ge12\), a removable factor gives a reduction to the
full-alternating passport in degree \(n-2\). Move that factor to
the last position by positive Hurwitz swaps. The other factors
retain their permutations and still fix its two exclusive
leaves. The fixed sheet \(z\) cannot be exclusive because a
sheet moved by only one factor is moved by the whole product:
every later permutation after that first move fixes the sheet
and cannot send another sheet back to it. Thus the two leaves
lie on the long cycle. The predecessor rule deletes them from
that cycle and leaves product type \((n-3,1)\). Removing two
leaves from the graph preserves connectedness. The resulting
group is transitive, and its fixed-point stabilizer contains
the long cycle on all other sheets, so it is 2-transitive and
primitive. Jones's same proposition gives the full \(A_{n-2}\).
In reverse, a two-leaf extension of this lower passport preserves
the stated upper passport exactly when its two distinct old
anchors lie on the lower long cycle. Using the lower fixed
sheet as one anchor instead creates a 2-cycle there. If both
anchors lie on the long cycle, the resulting tuple is transitive
and the upper group is full by the same primitive group argument.

It follows that any even orbit of relative size at least
\(e^{-c'n}\) for fixed \(0<c'<1/240\) admits this reduction for
all sufficiently large even \(n\). The sufficiency follows from
the uniform exceptional fraction. The sizes of all even orbits
and the equivalence of all admissible anchor extensions remain
separate obligations.

A finite check of the weighting used independent enumerations
of the transitive degree-six fixed-product tuples and of sampled
trees with one extra edge. At a fixed 5-cycle with one fixed
sheet, \(120\) tuples have cycle length two and \(80\) have
cycle length three. The two lists agree with their individual
tree-sampling multiplicities two and three, respectively,
giving \(480\) presentations. This check concerns the
unicyclic multiplicity at small degree.

## Scope of the conjugate-power mechanism

Féray's [pinned result](https://github.com/1337isnot1337/solved-conjectures/blob/758066055ef49c2c66cb1087918ae62ad75e7ac1/results/feray-conjugate-power-word/PROOF.md)
uses \(\sigma\tau\sigma^2\tau^2=X^2Y^{-1}XY\) for \(X=\sigma\tau\),
\(Y=\tau\), a bijective change of two independent uniform permutations.
Given the conjugacy class of \(X\), the power \(X^a\) is uniform on
its image class by equivariance. Conditional on any particular \(X\)
in that class, \(Y^{-1}X^bY\) is uniform on its image class, a law
independent of the choice of \(X\). This proves the conditional
independence used by the class Fourier argument. The same Schur
averaging mechanism yields the exact coefficient (3) for independent
class variables before conditioning on their product.

The mixing estimate in the external theorem requires few cycles,
including fixed points. A double transposition has \(N-2\) cycles.
The product of two DTs moves at most eight labels, whereas the uniform
law on \(A_N\) assigns probability at most
\(2\sum_{j=0}^8\binom Nj j!/N!\) to that event, tending to zero.
The variation distance for these two classes therefore tends to one.
The exact census and removable-locus density do not imply that a
fixed-product braid orbit intersects the removable locus. In
particular the pure and marked mixed uniform connectivity obligations
remain open.

For verification of the external rate, its two analytic inputs are the
uniform character bound in Teyssier–Thévenin Theorem 1.7 and the
fixed-rank zeta bound (6.3) with rank cutoff one. If both classes have
\(O(\log N)\) cycles, put \(q=(\log\log N+O(1))/\log N\). Schur
averaging and \(\ell^2\) Plancherel bound four times the squared TV
by \(\sum_{d_\chi>1}d_\chi^{-(2-4q)}\). The zeta bound at exponent one
and the \(\ell^s\) norm inequality give \(O(N^{-2}(\log N)^4)\).
For a uniform permutation \(X\), the cycle insertion identity gives
\(\mathbb E 2^{\#\mathrm{cycles}(X)}=N+1\). Powers of fixed nonzero
exponent multiply the cycle count by at most the absolute exponent.
Markov's inequality and parity averaging therefore give
\(O((\log N)^2/N)\) for the exact word \(X^2Y^{-1}XY\).
The two analytic inputs have no role in (1) or (2).
[Teyssier–Thévenin, v4](https://arxiv.org/html/2411.04347v4)

The leading-rate supplement retains the standard character \(F-1\).
For \(N\ge3\), counting marked disjoint cycles gives
\(\mathbb E[(F(X)-1)(F(X^2)-1)]=1\), since
\(F(X^2)=F(X)+2C_2(X)\), \(\operatorname{Var}F=1\), and
\(\operatorname{Cov}(F,C_2)=0\). Its sign twist vanishes when two
unused labels remain. Rank-two character polynomials in \(F,C_2\)
have bounded product moments at \(X,X^2\), so their Fourier
coefficients are \(O(N^{-2})\). At rank cutoff three the same
zeta estimate bounds the high-rank \(L^1\) remainder by
\(O(N^{-3}(\log N)^6)\) on the logarithmic-cycle event.
Using the threshold \(6\log_2 N\), the complement has probability
\(O(N^{-5})\) and discarded low ranks have dimensions \(O(N^2)\).
Thus the leading TV is the mean of \(|F-1|/(2(N-1))\) up to
\(O(N^{-2})\), or \(e^{-1}/(N-1)+O(N^{-2})\) by inclusion–exclusion.
[Pinned leading-rate proof](https://github.com/1337isnot1337/solved-conjectures/blob/758066055ef49c2c66cb1087918ae62ad75e7ac1/results/feray-conjugate-power-word/SHARP-RATE.md)

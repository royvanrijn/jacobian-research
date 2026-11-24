# Prefix complexity of random prime semigroups

Fix an integer base \(b\ge2\), \(0<\gamma<1\), and \(\lambda>0\).
Always admit \(2\) and independently admit each odd prime \(p\)
with probability \(\min(1,\lambda p^{-\gamma})\). Let \(C\) be the
multiplicative semigroup so generated, including \(1\), and put
\(\alpha=1-\gamma\). This is the model of
\(\texttt{ARITH-RANDOM-PRIME-MULTIPLIER-COUNT}\). The membership
language at length \(L\) consists of the zero-padded base-\(b\)
words for \(C\cap[0,b^L)\).
[Multiplier model and counting law](RANDOM_PRIME_SEMIGROUP_LOGARITHMIC_TRANSITIONS.md)

**Theorem.** Almost surely, for all sufficiently large \(L\),
every complete read-once branching program recognizing this
language has width
\[
 W_{\rm nd}\gg_{b,\gamma,\lambda}\frac{b^{\alpha L/2}}{\sqrt L}
 \qquad\text{or, in the deterministic case,}\qquad
 W_{\rm det}\gg_{b,\gamma,\lambda}\frac{b^{\alpha L}}{L^2}.              \tag{1}
\]
The program may depend on \(L\) and on the realized semigroup.
For deterministic programs the minimum width has exponential
rate exactly \(\alpha\), meaning
\(\log_b W_{\rm det}^{\min}/L\to\alpha\).
The bounds (1) also hold for programs that agree with \(C\)
only on the prime inputs in \([b^L/2,b^L)\).

## Complete read-once means a condition on graph paths

The programs are finite directed acyclic graphs with a start and
a single accepting sink. A transition reads a specified digit
position and is enabled by one specified value. Every directed
start-to-sink path queries each of the \(L\) positions exactly
once. A word is accepted when some such path has all its
transitions enabled. A deterministic state specifies one
queried position and at most one successor per possible digit.
The digit order can vary with the path.

Remove vertices on no start-to-sink path. For a retained vertex
\(v\), fix one suffix path to the sink. Every prefix to \(v\)
concatenates with that suffix to form a directed start-to-sink
path, so every prefix reads precisely the positions
complementary to the suffix. Thus \(v\) has a fixed read set
\(J_v\subseteq[1,L]\), and belongs to the layer \(|J_v|\).
Whenever a transition joins \(v\) to \(w\), its queried
position is the unique element of \(J_w\setminus J_v\).
Width is the largest number of retained vertices in a layer.

This is the path property used in the pinned Erdős–Turán
finite-prefix proof. There, disjoint alternative path stretches
give independent digit substitutions whose complementary choices
produce the same two-term sum. The proof counts a word by fixed
endpoint data plus edges that meet every ambiguous stretch.
One must apply the path property to graph paths, including
concatenations of prefixes and suffixes. If one assumes it
only for paths enabled by an input, the read set can vary at a
vertex. For instance, take two positions. The start can reach
one middle vertex by testing digit 0 at either position.
From the middle vertex the sink can be reached by testing digit
1 at either position. The enabled accepting paths on \(01\)
and \(10\) are complete read-once, but their middle read sets
differ. Concatenating the other branches gives inconsistent
two-query paths. The result here uses the syntactic complete
read-once model.
[Pinned Erdős–Turán finite-prefix note, Read-once branching programs](https://github.com/1337isnot1337/solved-conjectures/blob/d2cd4d3738570244ee1830c1d02236d93ecac11c/research/erdos-turan-additive-basis/FINITE-PREFIX-COMPLEXITY.md#read-once-branching-programs)

## Probability of a prescribed prime restriction

Set \(X=b^L\), and let \(M_L\) count the admitted primes in
\([X/2,X)\). For large \(L\) all these primes are odd and
their probabilities equal \(\lambda p^{-\gamma}\).
Ordinary prime number asymptotics and partial summation give
\[
 \mu_L=\mathbb E M_L
 \sim \frac{\lambda(1-2^{-\alpha})}{\alpha}
                 \frac{X^\alpha}{\log X}.                             \tag{2}
\]
Indeed the weighted prime sum on this fixed dyadic interval
has main term
\(\lambda\int_{X/2}^{X} t^{-\gamma}(\log t)^{-1}\,dt\).
Independence gives \(\operatorname{Var}M_L\le\mu_L\).
Chebyshev and \(\sum_L\mu_L^{-1}<\infty\) imply
\(M_L/\mu_L\to1\) almost surely, applying Borel–Cantelli to a
countable sequence of positive error thresholds. No independence
between different length windows is needed.
[Classical prime number theorem, fixed-modulus input with modulus one](../foundations/extended-geometry/sources/soprounov_pnt_arithmetic_progressions.pdf)

A prime is in the semigroup precisely when it was admitted.
Consequently any fixed acceptance language prescribes an
assignment of independent random labels on the prime inputs
in \([X/2,X)\). If its prescribed acceptance set there has
size \(m<\mu_L/2\), it cannot agree on the event
\(M_L\ge\mu_L/2\). If \(m\ge\mu_L/2\), the probability that
the prescribed admitted primes are all admitted is at most
\[
                 \bigl(\lambda2^\gamma X^{-\gamma}\bigr)^m
                       \le \exp(-c X^\alpha)                         \tag{3}
\]
for a fixed \(c>0\) depending on \(\gamma,\lambda\).
Requiring the prescribed rejections as well can only
decrease this probability.

## Counting descriptions of width \(W\)

Give every layer \(W\) named slots, padding unused slots
and fixing the start and sink names. In the nondeterministic
case record \(J_v\) as an \(L\)-bit mask for each of the
\((L+1)W\) slots. For each ordered endpoint pair in
consecutive layers record the subset of the \(b\) possible
digit values carried by its transitions. The endpoint
masks already determine the query position. This uses
\[
              O_b(L^2W+LW^2)
\]
bits for a nondeterministic description, including the
padding. In the deterministic case it suffices to record
the queried position and the \(b\) destinations or missing
transitions at each vertex. This uses
\[
             O_b\!\bigl(LW(\log(2L)+\log(2W))\bigr)
\]
bits. Naming slots and ignoring validity restrictions only
increase the counts of descriptions.

For \(W\le c_1 X^{\alpha/2}/\sqrt L\), the
nondeterministic description length is bounded by a small
constant times \(X^\alpha\), by choosing \(c_1>0\) small
and then taking \(L\) large. The term \(L^2 W\) is
\(o(X^\alpha)\). For
\(W\le c_2 X^\alpha/L^2\) the deterministic description
length has the same bound, by choosing \(c_2>0\) small.
Exponentiating the description counts and summing (3)
shows that the probability of a program of either of
these widths agreeing on the prime inputs, together with
\(M_L\ge\mu_L/2\), is \(O(\exp(-c X^\alpha/2))\).
Another application of Borel–Cantelli proves (1).

A deterministic digit trie has width at most
\(|C\cap[0,b^L)|\) after identifying all accepting leaves
with the sink. Its digit positions are queried in the
same order on every path. The retained multiplier counting
law gives trie width \(O_\omega(b^{\alpha L}L^{\lambda-1})\).
This and (1) give the exact deterministic exponential rate.
The lower bounds continue to hold after any finite change of
the admitted prime set on the same probability-one event,
since the prime inputs in the proof then eventually lie
outside the changed set.

## Grammars of fixed fanout and production rank

For each fixed fanout \(f\) and maximum production rank \(R\),
the minimum size \(G_L\) of a linear multiple context-free
grammar recognizing the length-\(L\) membership language
satisfies almost surely
\[
                       G_L\gg_{b,f,R,\gamma,\lambda}
                            b^{\alpha L}/L .                          \tag{4}
\]
Here size is the number of nonterminals plus productions.
Only the set of words of length \(L\) generated by the
grammar is prescribed.

To count descriptions, first use the pinned projection
normalization to make the grammar nondeleting. This costs
a factor at most \(2^f\) in size and preserves rank and
fanout. It projects each nonterminal onto all nonempty
sets of its components. In a projected rule each remaining
child component is used once, and completely unused
productive children can be dropped. Projecting derivations
and lifting the child projections proves preservation of
the start language. Then remove productions unused in all
length-\(L\) derivations, which preserves the length-\(L\)
language. Every literal digit in a remaining production
appears in some completed word of length \(L\), since
nondeletion retains that occurrence. Thus each production
has at most \(L\) literal digits. Its number of child
component symbols and terminal gaps is \(O_{f,R}(1)\).
Production shape, symbol indices, gap lengths and
literal digits therefore use \(O_{b,f,R}(G(L+\log(2G)))\)
bits to describe a normalized pruned grammar of size
\(O_f(G)\). For \(G\le c_3 X^\alpha/L\) with a suitably
small constant \(c_3>0\), the same union bound using (3)
gives (4). The grammar with one terminal-only production
per accepted word has size \(1+|C\cap[0,b^L)|\). Thus,
in any class allowing these productions, the exponential
rate of the minimum grammar size is also \(\alpha\).
[Pinned projection proof](https://github.com/1337isnot1337/solved-conjectures/blob/d2cd4d3738570244ee1830c1d02236d93ecac11c/research/erdos-turan-additive-basis/README.md#grammar-and-claim)

## Exact Kolmogorov information rate

Let \(V_L=(1_{\{n\in C\}})_{0\le n<X}\), in numerical order,
where \(X=b^L\). For a fixed universal binary decoder, let \(K\)
be either plain or prefix-free Kolmogorov complexity in bits.
Almost surely,
\[
 K(V_L)=
 \left(\frac{\gamma\lambda}{\alpha\log 2}+o(1)\right)X^\alpha .
                                                                    \tag{5}
\]
The same limit holds if the decoder is given \(L\). The base \(b\)
is fixed and known, and \(\gamma,\lambda\) may be arbitrary real
parameters in the ranges at the start of the note. The decoder
uses no information about the realized prime labels except its
binary description. Computation time is unrestricted.

Let \(Y_L\) be the vector of admission labels \(\xi_p\) for odd
primes \(p<X\), in increasing order, and write
\(q_p=\mathbb P(\xi_p=1)\). One can compute \(Y_L\) from \(V_L\)
by extracting its prime entries. Conversely, one can compute
\(V_L\) from \(L\) and \(Y_L\) by factoring each integer in the
interval and using the admission of \(2\). The two vectors
therefore have equal point probabilities at corresponding
outcomes. Write \(p_L(v)=\mathbb P(V_L=v)\).

Write \(S_L=\sum_{2<p<X}\xi_p\log p\) and
\(m_L=\sum_{2<p<X}\xi_p\). The prime number theorem and partial
summation, as in (2), give
\[
 \mathbb E S_L\sim\frac{\lambda}{\alpha}X^\alpha,\qquad
 \mathbb E m_L\sim\frac{\lambda}{\alpha}\frac{X^\alpha}{\log X}.
\]
For the first estimate, apply partial summation to
\(\sum_{p\le t}\log p\sim t\) against \(t^{-\gamma}\).
Independence gives
\[
 \operatorname{Var} S_L\le(\log X)\mathbb E S_L,\qquad
 \operatorname{Var} m_L\le\mathbb E m_L.
\]
Each variance divided by the square of the corresponding mean
is \(O((\log X)/X^\alpha)\), which is summable over \(L\).
Chebyshev and Borel–Cantelli then give \(S_L\sim\mathbb E S_L\)
and \(m_L\sim\mathbb E m_L\) almost surely.

The information \(I_L=-\log_2 p_L(V_L)\) satisfies
\[
 (\log 2) I_L
      =\gamma S_L +O_{\lambda,\gamma}(1+m_L+\mathbb E m_L).
\]
Indeed each sufficiently large admitted prime contributes
\(-\log q_p=\gamma\log p-\log\lambda\), and each such rejected
prime contributes \(-\log(1-q_p)\le2q_p\) once \(q_p\le1/2\).
The remaining finitely many primes contribute a bounded amount
on their support, including zero for any deterministic label.
Thus \(I_L/X^\alpha\to h=\gamma\lambda/(\alpha\log 2)\)
almost surely.

For the lower bound, a fixed decoder has at most \(2^{k+1}\)
binary descriptions of length at most \(k\), and each outcome
with \(I_L\ge(h-\varepsilon)X^\alpha\) has probability at most
\(2^{-(h-\varepsilon)X^\alpha}\). Taking
\(k=\lfloor(h-2\varepsilon)X^\alpha\rfloor\) for
\(0<\varepsilon<h/2\), a union bound for these short descriptions
and this information event gives probability at most
\(2^{1-\varepsilon X^\alpha}\). The calculation also applies
with \(L\) provided to the decoder, because each description at
that fixed \(L\) yields at most one vector. By Borel–Cantelli
and a countable sequence of \(\varepsilon\) tending to zero,
\(\liminf K(V_L)/X^\alpha\ge h\). The lower bound applies to
any fixed partial computable binary decoder, and hence
simultaneously to all such decoders by countability.

For the upper bound let \(\Pi_L\) be the number of odd primes
below \(X\). Give \(L\) and \(m_L\) in a self-delimiting header
using \(O(L)\) bits, then give the lexicographic rank of the
selected subset among the \(m_L\)-element subsets of these
\(\Pi_L\) primes, using exactly
\(\lceil\log_2\binom{\Pi_L}{m_L}\rceil\) bits. The header fixes
the payload length, so this code is prefix-free. One decoder
enumerates the primes from \(L\), recovers the subset from its
rank, and reconstructs the membership vector by factorization.
The decoder does not need \(\gamma\) or \(\lambda\).
Since \(\Pi_L\sim X/\log X\) and \(m_L\sim\lambda
X^\alpha/(\alpha\log X)\), almost surely
\[
 \begin{aligned}
 \log_2\binom{\Pi_L}{m_L}
 &\le m_L\log_2(e\Pi_L/m_L)\\
 &=\left(\frac{\gamma\lambda}{\alpha\log 2}+o(1)\right)X^\alpha .
 \end{aligned}
\]
A fixed universal decoder can run this decoder with a constant
additional number of bits, giving the upper bound in (5).
In particular a uniform algorithm computing membership in \(C\)
would give descriptions from \(L\) alone, so almost surely
\(C\) has no computable membership predicate.

## Scope for the additive questions

These are lower bounds for exact numeral membership
programs and grammars in the repository's random prime
model. They use the independent labels on prime inputs.
The Erdős–Turán bound instead derives program complexity
from bounded counts of representations. The random
semigroup width bound supplies no such count bound.
In particular it neither yields a pointwise asymptotic
for representations nor proves the critical natural
density of \(kC\). An arbitrary set of powerful integers
likewise supplies no independent prime labels, since
every prime is absent from its membership language.
[Random-prime representation and density questions](RANDOM_PRIME_SEMIGROUP_LOGARITHMIC_TRANSITIONS.md)
· [Powerful-integer representation targets](POWERFUL_SUM_GRAPHS.md)

# Growing private rectangle families around one additive edge

A near-full graph on complete primitive abc values has many copies of
any fixed additive rectangle around a typical edge. Polylogarithmically
many copies can share that edge and have private prime packets on every
other complete value. This does not require a complete bipartite
subgraph with both dimensions growing.

## Rooted rectangles with private packets

Put \(n=X^L\), \(\gamma=1/6\), and take a graph of actual positive
primitive triples \(a+b=c\) in \(X/2<c\le X\). Write \(A,B,C\)
for the complete-value dictionaries. Assume
\[
 |A|,|B|,|C|\le n^{\gamma+\kappa},\qquad
 |\mathcal E|\ge n^{2\gamma-\kappa},\qquad
 \operatorname{rad}(W)\le n^{2\gamma+\kappa}\quad(W=a,b,c).
 \tag{BR1}
\]
Fix integers \(r,s\ge2\) and constants \(D>0,\ \xi>0\).
There are \(\kappa_0,L_0>0\) depending only on these parameters
such that for every fixed \(0<\kappa<\kappa_0\) and \(0<L<L_0\),
the following holds for all sufficiently large \(X\). Apart from
at most \(n^{2\gamma-2\kappa}\) anchor edges \((a_0,b_0,c_{00})\),
one can choose \(m=\lfloor(\log X)^D\rfloor\) copies of
\(K_{r,s}\) in the \(A,B\) projection, rooted at that edge.
Every cell is its actual primitive sum. For example, each copy
uses rows \(a_0,a_1,\ldots,a_{r-1}\), columns
\(b_0,b_1,\ldots,b_{s-1}\), and targets \(c_{ij}=a_i+b_j\).
Outside the three anchor values all complete values are distinct,
within and between copies, and have radical greater than
\(n^{1/3-\min(1/12,\xi/2)}\).

For each new complete value \(W\) there is a factor \(P_W\)
formed from whole prime powers of \(W\), such that
\[
                 \log_n\operatorname{rad}(W/P_W)\le\xi.       \tag{BR2}
\]
Within a copy, the support of \(P_W\) is disjoint from the
full supports of all other new values. It is also disjoint
from the full supports of the anchor and all earlier copies.
A prescribed prime pool of cardinality at most \((\log X)^{D+1}\)
can be included among the supports to avoid. The statement is
uniform in this pool, the graph, and each eligible anchor.
All original radical caps continue to apply to every cell.

### How many rectangles contain an anchor

Write \(v=r+s,\ e=rs,\ t=v-2\), and put \(K=3v+4e+3\).
Call an anchor deficient if it is the distinguished first cell of
fewer than \(n^{t\gamma-K\kappa}\) ordered rectangles with
distinct rows and distinct columns. Suppose the deficient-edge
subgraph has at least \(n^{2\gamma-2\kappa}\) edges.
Both of its active parts have size at least \(n^{\gamma-3\kappa}\)
by (BR1), and its edge density is at least \(n^{-4\kappa}\).
For a bipartite graph with part sizes \(M_1,M_2\) and density
\(p\), two applications of the power mean give at least
\(M_1^r M_2^s p^{rs}\) labeled homomorphisms from \(K_{r,s}\).
Indeed averaging the common neighborhood size of \(r\) rows
first gives at least \(M_2p^r\), and averaging its \(s\)-th power
counts the column choices. Thus the deficient-edge subgraph has
at least \(n^{v\gamma-(3v+4e)\kappa}\) such homomorphisms.
Repeated row or column choices cost at most
\(O_{r,s}(n^{(v-1)(\gamma+\kappa)})\), smaller by a fixed
power for small \(\kappa_0\). On the other hand, counting its
injective rectangles by their first cells gives at most
\(n^{2(\gamma+\kappa)}n^{t\gamma-K\kappa}\). The exponents
contradict the lower bound for sufficiently large \(X\).
[Fox–Sudakov, Dependent Random Choice, §2](https://people.math.ethz.ch/~sudakovb/dependent-random-choice.pdf)

### Choosing a new copy

Fix a nondeficient anchor. Throughout the construction the
initial set of candidate copies comprises at least
\(n^{t\gamma-K\kappa}\) distinct assignments to the \(t\)
free row and column variables. Let \(\mathcal P\) contain
the prescribed pool, the full supports of the three anchor values,
and the full supports of every new complete value in the copies
already chosen. This pool has
\(O_{r,s,D}((\log X)^{D+1})\) primes since fewer than \(m\)
copies have been chosen and each value is at most \(X\).
Its cardinality suffices for effective-prime pricing.
The height of the product of these primes need not be polynomial.

Put \(q=t+rs-1\), the number of new complete values in each copy,
and \(\tau=\min(1/8,\xi)/(2q+2)\). Test every new complete
value against \(\mathcal P\), and test every pair of different
new complete values against one another. A test fails at
normalized common radical mass at least \(\tau\).

For a pool of \(O((\log X)^{D+1})\) primes, choose the fixed
price multiplier \(s_0=1-1/(2(D+2))\) on that pool. The
partition sum condition follows from
\(\sum_{j\le O((\log X)^{D+1})}j^{-s_0}=o(\log X)\).
In a fixed-letter fibre, a tested free letter carrying mass
at least \(\tau\) on the pool has effective cap at most
\[
                      2\gamma+\kappa-(1-s_0)\tau .
\]
Equation (GP4) counts it at exponent at most
\(\gamma+\kappa/2-(1-s_0)\tau/2+\eta+\varepsilon/L\).
All three raw caps in (BR1) are retained and give the fixed
raw-budget bound needed by that estimate. Take \(\kappa_0\)
and the allowances sufficiently small that this exponent is
at most \(\gamma-\sigma\), for one fixed \(\sigma>0\)
depending only on \(D,\xi,r,s\). The estimate is uniform in
the fixed complete letter and the prescribed pool.
[Effective fixed-letter count (GP4)](ABC_GROWING_PRIVATE_PRIME_PACKETS.md#high-recurrence-at-the-global-and-fixed-letter-exponents) · [Effective prime pools (EB3–EB5)](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)

For a test against \(\mathcal P\), keep all but one free variable
of the tested value fixed. The new values are of the form
\(a_i,\ b_j,\ a_i+b_j\), with an anchor summand fixed when its
index is zero. Each variable assignment then gives a distinct
triple in one original fixed-letter fibre. For example if the
varying variable is \(a_i\), hold the column \(b_j\) fixed and
use the cell there. If testing the row value itself, any
column in its copy can be used. Hence each such test rejects
at most \(n^{t\gamma+(t-1)\kappa-\sigma}\) assignments, up
to constants and an arbitrarily small fixed allowance.

For a test between two new values, choose a free variable
occurring in just one of their affine forms when one exists.
Hold the other variables fixed and use the prime support of
the other value as the fixed pool. This again gives the same
bound. If two forms use precisely the same free variables,
every pair of distinct forms of that type is incident in
an original cell. For instance the forms \(a_i\) and \(a_i+b_0\)
are coprime by primitivity. If both forms involve two free
variables, their equal variable sets would make them the
same cell and there is no pair to test.

Reject a new value with radical at most
\(n^{1/3-\min(1/12,\xi/2)}\) by the same fixed-letter count
with its raw cap reduced, fixing all but one of its variables.
This costs the same bound after reducing \(\sigma\) if needed.
Equality with any already chosen complete value fixes at least
one free variable once the others are fixed. Equality of distinct
new forms likewise fixes one variable if their variable sets
differ. For equal variable sets their difference is a nonzero
anchor value in the incident case just discussed. All forbidden
numerical coincidences together cost
\(O_{r,s,D}((1+m)n^{(t-1)(\gamma+\kappa)})\).
The finite collection of prime tests and these coincidences
reject fewer than \(n^{t\gamma-K\kappa}\) assignments once
\(\kappa_0\) is sufficiently small compared with
\(\sigma/(K+t)\) and \(\gamma/(K+t)\). The comparison is uniform
at each step. Choose one surviving copy. Strip from each of
its new values the whole powers of primes meeting \(\mathcal P\)
or one of the other new values of that copy. At most \(q\tau\le\xi\)
radical mass is stripped. Repeat with the enlarged pool.
This proves (BR2).

## Multiplicative relations and repaired classes

Each new value has a nonempty packet since its radical exceeds
\(n^{1/4}\) and the stripping cost is at most \(1/16\).
The new values are multiplicatively independent modulo the subgroup of
\(\mathbb Q^\times\) supported on the anchor primes and the
prescribed pool. Indeed in any purported relation select the
last copy with a nonzero integer coefficient. Select a value
with a nonzero coefficient there and a prime in its packet.
This prime divides no other value in that copy and no value
in earlier copies. It therefore gives a nonzero valuation
of the relation. The same reasoning works for any subcollection
on which the corresponding packets are nonempty.

For a prime \(q_0\), suppose the designated values have
normalized mass greater than \(\xi\) on primes whose exact
valuation is not divisible by \(q_0\). Their classes are
linearly independent over \(\mathbb F_{q_0}\) in
\(\mathbb Q^\times/(\mathbb Q^\times)^{q_0}\) modulo the
anchor and prescribed prime subgroup. The packet keeps at
least one of these primes, so the last-copy argument yields
a nonzero valuation modulo \(q_0\). For least repairs of
order \(R\) with \(q_0\mid R\), the same assertion holds
for the repaired multipliers if the mass condition counts
only \(0<v_p(W)<R\) with \(q_0\nmid v_p(W)\).
On such a prime the multiplier has valuation
\(R-v_p(W)\equiv-v_p(W)\bmod q_0\). These hypotheses refer
to the exact original valuations.
[Least target repair classes on translated forms](ABC_ADDITIVE_TRANSLATION_CUBES.md#large-dictionaries-of-weighted-least-repair-classes)

## Why growing common rectangles require additional assumptions

Let \(s(X)\) be any integer tending to infinity. For any actual
complete additive graph with part sizes at most \(X^B\) and
at least \(X^\alpha\) edges for fixed \(B,\alpha>0\),
retain each edge independently with probability
\(p=X^{-(2B+1)/s(X)}=X^{-o_X(1)}\). In any chosen bipartite
projection, the expected number of copies of \(K_{s(X),s(X)}\)
in the retained graph is at most
\[
                  X^{2Bs(X)}p^{s(X)^2}=X^{-s(X)}.
\]
For a fixed marginal degree initially at least \(X^\alpha\),
the probability of retaining less than half its expected
degree is exponentially small in \(pX^\alpha\) by a Chernoff
bound. A union bound over all such vertices and the total
edge count shows that some thinning has no growing biclique,
retains \(X^{-o_X(1)}\) of the original edges, and retains the
exponent of every power-large degree. The same works in all
three projections at once. The numerical labels, primitivity,
and every original radical cap on retained edges are inherited.
Thus if capped near-saturating graphs exist, they can be
thinned to ones with the same count exponent and without
any given growing biclique. This conditional obstruction
does not construct a capped near-saturating graph.

The rooted copies in (BR2) have a fixed number of edges per
new value even when their total size grows. Their combined
packet LCM can have quasipolynomial height and cannot simply
be used as a single polynomial-height input. The construction
uses the polylogarithmic cardinality of the prime pool and
keeps all three radical caps in each fixed-letter count.

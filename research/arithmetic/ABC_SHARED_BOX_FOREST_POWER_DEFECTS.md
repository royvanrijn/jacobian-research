# Joint power defects on shared-value abc boxes

On a shared-value box, a cycle with a cheap power defect can
recover a potential from the other endpoints of that cycle.
Several such recoveries can use the same prime-support record.
A comparison forest gives a simple instance. A large kernel of
the valuation map gives a more flexible way to choose the order
of recovery.

Use the setup of
[shared-value boxes (LB1)](ABC_SHARED_VALUE_RESIDUE_BOXES.md#a-finite-difference-array-theorem).
The three axes, each of size \(m\ge2\), have residue alphabet
bounds \(R_1,R_2,R_3\), with \(R=R_1R_2R_3\). For a specified
pair of opposite residue coordinates, a complete letter has at
most \(K\) possible values. The three arrays are the positive
edge differences of three groups of integer potentials. Orient
each edge from its smaller potential to its larger one, and
write \(L_e\) for its positive length. An ordered box uses
distinct residues within every axis, so its potentials are
distinct there. They are distinct across axes by positivity of
the edge differences. All potential differences have magnitude
\(O_m(X)\).

Fix a prime pool \(P\) whose prime product is at most \(X^B\)
for some fixed \(B\). For a cycle with signed edge weights
\((b_e)\), its defect at order \(q\ge2\) is the product of primes
\(\ell\notin P\) for which
\(\sum_e b_e v_\ell(L_e)\) is nonzero modulo \(q\).
The joint defect of several cycles is the product of the
distinct primes in the union of these supports. The same
convention applies to the joint defect of any subspace of
circulations modulo a prime.

## Recovering a potential with an already recorded support

Suppose a set of potentials is known up to common translation,
a new vertex has unknown potential \(x\), and a circulation
modulo a prime \(p\) is supported on edges among the known
vertices and edges from \(x\) to known vertices. Suppose it has a
nonzero weight on at least one edge incident to \(x\), and that
its defect is supported on a previously recorded set \(S\) of
primes outside \(P\). Its divergence at \(x\) is zero, so
another edge at \(x\) has a nonzero weight. There are only
\(O_m(1)\) choices of the two incident edges.

Add to \(P\cup S\) the prime supports of the fixed edge labels
and of the nonzero differences of the known neighbors of \(x\).
Every prime outside the enlarged pool divides at most one
form incident to \(x\). On either of the two selected forms,
its valuation is a multiple of \(p\), by reduction of the
circulation product modulo \(p\). Thus both canonical
\(p\)-th-power-free coefficients are on the enlarged pool.
The record \(S\) itself is supported on the actual edge labels,
so its prime product is at most \(X^{3m^2}\).

Uniform fixed-support enumeration gives \(X^\epsilon\) choices
for each coefficient, for any fixed positive \(\epsilon\).
For a fixed pair of coefficients, the forms equal \(A y^p\)
and \(A' z^p\). Eliminating \(x\) gives a fixed nonzero signed
sum of these two powers because their affine roots are distinct.
The uniform binary norm estimate gives \(X^\epsilon\) choices
for bounded orders \(p\), with the allowance decreased at each
application. For larger \(p\le\log_2 X\), enumerate one
positive root in \([1,X^{1/p}]\) and pay only \(\log_2 X\)
for all these orders. For \(p>\log_2 X\) the selected forms
themselves are supported on the enlarged pool, so count their
values directly. Consequently there are \(O(X^\epsilon)\)
possible \(x\), uniformly over all orders and circulation
weights, once \(S\) and the known potentials are fixed. The
constants depend on \(m,B,\epsilon\). At every step it suffices
that a suitable circulation exists for the box being counted.
[Binary norm counting with moving coefficients](ABC_RECTANGLE_POWER_CLASS_DEFECT.md#the-fixed-value-bound-with-moving-coefficients)
· [Fixed-support enumeration](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)

## Propagation along a comparison forest

**BF1.** Choose an axis \(i\), and on its \(m\) positions take a
comparison graph with \(c\) connected components, counting
isolated positions as components. An edge of the comparison
graph specifies an oriented four-cycle using its two positions
and two distinct vertices outside axis \(i\). Both outside
vertices may be in the same other axis. The cycles can use
different outside vertices and each can have its own order
\(q\ge2\). If their joint defect is at most \(Z\), the number of
such ordered boxes is
\[
 O_{m,B,\epsilon}\!\left(
 Z X^\epsilon\,\frac{R^m K^{3m-1}}{(K R_i)^{m-c}}
 \right).
 \tag{BF1}
\]
The graph, outside vertices and orders may vary between boxes.

For the proof it suffices to retain a spanning forest of the
comparison graph. Enumerate the support \(S\) of its joint
defect by its squarefree product at cost at most \(Z\). Fix all
vertices of the other two axes and one root of each comparison
tree. The induced graph of these known potentials is connected
since it contains both other axes. The known residues cost at
most \(R^m/R_i^{m-c}\). A spanning tree on the known potentials
has \(3m-(m-c)-1\) edges, so their labels cost at most
\(K^{3m-(m-c)-1}\) and determine the known potentials up to
translation.

Traverse each rooted comparison tree outward. A comparator
cycle for a new vertex uses its known parent and the two
known outside vertices. Its two edges at the new vertex have
unit signed weights, so reducing modulo any prime divisor of
its order meets the recovery hypothesis. Use the already
recorded support at every step. With \(m\) fixed, the choice
of comparison forest, roots and outside vertices costs only
a constant, and the \(X^\epsilon\) allowances can be divided
among the finitely many recoveries. This proves (BF1), including
when the comparison graph is disconnected.

## Recovery from a bound on residual valuation rank

Let \(V\) be the \(3m\) vertices of the complete tripartite graph
and let \(\mathcal C_p\) be its circulation space over
\(\mathbb F_p\). For a prime set \(S\) outside \(P\), define
\[
 \Psi_{p,P\cup S}(b)=
 \biggl(\sum_e b_e v_\ell(L_e)\pmod p\biggr)_{\ell\notin P\cup S}.
\]
The rank in the next statement concerns this whole map from
\(\mathcal C_p\), including any independent cycles that have no
chosen small defect on their own.

**BF2.** Fix an integer \(r\ge0\) and a connected set of anchor
vertices containing \(k_i\) vertices in axis \(i\) and \(k\)
vertices in all. Suppose the remaining vertices admit an
ordering in which each has at least \(r+2\) neighbors among the
anchors and preceding vertices. The number of ordered boxes
for which there exist a prime order \(p\) and some prime set
\(S\) of product at most \(Z\) with
\(\operatorname{rank}\Psi_{p,P\cup S}\le r\) is
\[
 O_{m,B,\epsilon}\!\left(
           Z X^\epsilon K^{k-1}\prod_{i=1}^3 R_i^{k_i}
                         \right).
 \tag{BF2}
\]
Equivalently, a subspace of codimension at most \(r\) in
\(\mathcal C_p\) may be required to have joint defect at most
\(Z\). No basis or list of weights is fixed. The anchor pattern
and recovery order may vary between boxes at bounded cost.

To see this, take \(H=\ker\Psi_{p,P\cup S}\), which has
codimension at most \(r\). Discard primes in \(S\) that divide
no actual edge label and enumerate the remaining support once.
Record the anchor residues and recover their potentials up to
translation by assigning their spanning-tree edge values, at
cost \(K^{k-1}\prod_i R_i^{k_i}\).

At the next vertex \(x\) let \(D\) be its number of already
known neighbors. The known induced graph is connected. On
the \(D\) edges between \(x\) and the known vertices, any
assignment of weights with zero signed divergence at \(x\)
extends to a circulation supported on the known graph
together with \(x\). Indeed the prescribed flows have total
divergence zero on the known vertices. A spanning tree of
their connected graph can route this divergence by successively
balancing its leaves. Consequently projection to the \(D\)
incident edges has image of dimension \(D-1\) on the local
circulation space. Intersecting that space with \(H\) loses
at most \(r\) dimensions. If \(D\ge r+2\), there is still a
circulation in the intersection whose incident weight vector
is nonzero. Its support is entirely on known vertices and
edges from \(x\) to them, and its defect is supported on \(S\).
The preceding recovery count therefore applies. Continue
through the given ordering, absorbing all subpower costs in
\(X^\epsilon\). This proves (BF2).

The equivalence follows by taking a kernel in one direction,
and in the other by taking the union of the outside defect
supports of the subspace. Valuations are evaluated on residues
of the weights, so neither statement depends on choices of
integer representatives.

For example, if \(r=0\) start with one vertex from each of two
axes. Recover one vertex of the third axis first. Any
remaining vertex now has two known neighbors, so \(k=2\)
suffices. For general \(r\), one can take
\(\lceil(r+2)/2\rceil\) anchors per axis whenever that number
is at most \(m\). Every other vertex then has at least \(r+2\)
already known neighbors.

## Consequence at the simultaneous packet scale

In a simultaneous packet tuple, put \(n=X^L\),
\(d=s_1+s_2+s_3\) and \(t=1/6-d\), as in (LB4).
Use \(R_i=n^{s_i+O(\kappa+\eta)}\) and \(K=n^{t+O(\eta)}\).
The ordered lift count is at least
\(n^{md+(3m-1)t-O_m(\kappa+\eta)}\).
Thus for every fixed \(\sigma>0\), with sufficiently small
fixed tolerances and inner allowance, the boxes satisfying
(BF2) with
\[
              Z\le
 n^{\,\sum_i(m-k_i)(t+s_i)-\sigma}
\]
are power fewer than all the forced lifts. In (BF1) the
corresponding threshold is
\(Z\le n^{(m-c)(t+s_i)-\sigma}\). This quantifies the prime
weight needed to reduce the rank of the whole circulation
map or to impose defects along many related cycles. An
endpoint exclusion still needs an arithmetic upper bound on
these joint defect supports for a sufficiently large cohort of
actual boxes.
[Forced actual lifts and their count (LB4)](ABC_SHARED_VALUE_RESIDUE_BOXES.md#consequences-for-simultaneous-packet-strata)

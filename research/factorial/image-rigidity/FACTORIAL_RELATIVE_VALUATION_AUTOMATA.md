# Relative valuation and cancellation in factorial moment automata

Fix a finite support \(\{\nu_h\}\subseteq\mathbb N^d\) and
\(F=\sum_h c_hX^{\nu_h}\) over a finite extension of \(\mathbb Q_p\).
Let \(\pi\) be a uniformizer, write \(w=\operatorname{ord}_\pi\),
and put \(e=w(p)\). Choose a nonzero reference
\[
 H_n=r^n\prod_j(a_j n+b_j)!^{\epsilon_j},
 \qquad r\ne0,\quad a_j,b_j\in\mathbb N,\quad \epsilon_j\in\mathbb Z .
\]
This includes a pure term reference. For \(\sum_h m_h=n\) write
\[
 W(n,m)=\frac{n!}{\prod_h m_h!}\prod_h c_h^{m_h}
                        \prod_{i=1}^d T_i!,\qquad T_i=\sum_h\nu_{h,i}m_h .
\]

1. At each order one can compute the minimum term valuation relative
   to \(H_n\) and sum the actual normalized contributions at any
   fixed precision above that minimum. For fixed support, reference
   and precision this uses a bounded number of graph operations per
   input digit, with auxiliary integer arithmetic of polynomial bit
   complexity.
2. When the reference has slope at most every support column's slope,
   a finite graph decides if those relative term valuations have a
   uniform lower bound. In the bounded case a finite output automaton
   computes the actual normalized moment residues at every order,
   after clearing the bound by a single scalar. It then decides
   all starts of a fixed residue window by finite reachability.

These tests include arbitrary ties between summands. A surviving
residue zero window does not by itself give an exact zero window.

## 1. Valuation and unit data on one graph

Insert \(w(N!)=e(N-s_p(N))/(p-1)\) into the ratio, with \(s_p\) the
base-\(p\) digit sum. Put
\(a_h^*=(p-1)w(c_h)+e|\nu_h|-((p-1)w(r)+e\sum_j\epsilon_j a_j)\)
and \(C=-e\sum_j\epsilon_j b_j\). Then
\[
 \begin{aligned}
 (p-1)w(W(n,m)/H_n)&=\sum_h a_h^*m_h+C\\
 &\quad+e\left(\sum_h s_p(m_h)-s_p(n)-\sum_i s_p(T_i)
                          +\sum_j\epsilon_j s_p(a_j n+b_j)\right).
 \end{aligned}
 \tag{1}
\]
Read the digits of \(n,m_h\) from least to most significant.
The count sum, \(T_i\), and \(a_j n+b_j\) are checked or output by
bounded carries for fixed linear forms. At position \(\ell\) the
contribution is \(A p^\ell+B\), where \(A\) combines count digits
with the \(a_h^*\) and \(B\) is the indicated signed sum of output
digits. A fixed number of zero order digits beyond the canonical
length of \(n\) suffices to clear all carries. Each \(n,m\) has
exactly one track of that length. For existential tests one can
allow all padded lengths, since added padding changes no value in
(1). See [the pure-term valuation test](FACTORIAL_ALL_ORDER_PURE_TERM_VALUATION_TEST.md#factorial-terms-and-effective-slopes).

For the unit data fix precision \(L\) and let
\(J(N)=\prod_{1\le j\le N,\,p\nmid j}j\) in \(\mathcal O/\pi^L\).
Choose \(a\) with \(ae\ge L\). One block of \(p^a\) consecutive
integers has a fixed product of its factors prime to \(p\) in this
finite unit group. Thus \(J(N)\) is periodic, with a period \(p^a\)
times the order of the block product. For \(U(N!)=p^{-v_p(N!)}N!\)
the recurrence
\[
 U((pN+d)!)=U(N!)J(pN+d),\qquad
 v_p((pN+d)!)=N+v_p(N!),\qquad 0\le d<p ,
\]
gives a finite machine reading from most to least significant.
Keep \(v_p(N!)\) modulo the order of \(p/\pi^e\) in the same group.
Multiplying the unit output by \((p/\pi^e)^{v_p(N!)}\) gives the
unit after stripping powers of \(\pi\). Powers of coefficient
units are obtained by keeping exponents modulo their unit orders.
See [Granville, Sections 2–3](https://www.cecm.sfu.ca/organics/papers/granville/).

Reverse the unit machines to run alongside the carry tracks.
Guess their final states, check the reversed transitions, and
end at their original initial states. Each fixed digit track has
one accepting reversed run because the original machines are
deterministic. The guessed states give the unit of \(W(n,m)/H_n\).

## 2. Normalize at the minimum for an individual order

For the fixed order digits propagate through the graph. In each
state keep the least accumulated cost \(a\) and the sum of residues
relative to that cost. Terminal units can be carried as state
labels and multiplied in at the end. To use the scaled costs
directly, adjoin \(\theta\) with \(\theta^{p-1}=\pi\) if \(p>2\)
and use \(\theta=\pi\) otherwise. Work modulo \(\theta^{(p-1)L}\).
A transition adds \(A p^\ell+B\) to the cost. To combine two
incoming pairs \((a,u),(a',u')\), take \(b=\min(a,a')\) and use
\(\theta^{a-b}u+\theta^{a'-b}u'\).
Keep the minimum path cost even if this residue is zero. A common
later cost preserves all discarded relative precision, which
proves the update inductively.

After adding the common \(C\), terminal costs are multiples of
\(p-1\). This gives the actual residue modulo \(\pi^L\) above
the minimum term valuation in the original field. A nonzero
output certifies nonvanishing, including cases with cancellations
in some lower digits. All stored costs have bit length bounded
linearly in digit length. This step allows either sign of \(a_h^*\). For instance, with
\(F=x-y\) and \(H_n=n!\) every term has valuation zero. The
signed sum normalized by \(H_n\) is zero at odd orders and one
at even orders, as the direct sum of alternating signs also gives.

## 3. Exact criterion for a uniform factorial reference

Assume \(a_h^*\ge0\) for all \(h\), so every positional \(A\) is
nonnegative. Consider the subgraph of edges with \(A=0\) on vertices
reachable from the initial state in the full graph and able to
reach acceptance using only \(A=0\) edges. The valuations in (1)
are uniformly bounded below if and only if this subgraph has no
directed cycle of negative total \(B\).

A reachable negative cycle can be repeated and exited by the
prescribed zero-positional suffix. Its fixed prefix contains all
positive positional contributions at fixed digit positions.
Pumping the cycle sends (1) to \(-\infty\).
Conversely, without a negative cycle, any zero-positional suffix
ending at acceptance has additive cost bounded below by a fixed
\(-C_0\), by deleting its nonnegative cycles.
Let \(B_0=\max|B|\). If the last positive positional edge occurs
at index \(\ell\) starting from zero, then the full cost is at least
\[
                         C+p^\ell-B_0(\ell+1)-C_0 .
 \tag{2}
\]
With no positive positional edge it is at least \(C-C_0\).
This proves the equivalence and supplies an effective shift
\(\pi^M\) making every normalized term integral. It also shows
that at each fixed relative precision the counts on columns
with strictly larger slope can have nonzero digits only in a
computably bounded initial segment, independent of the order.
If all column slopes are strictly larger than the reference slope,
every column count is bounded at the fixed precision. Hence the
residue sequence is eventually zero. For a uniform certificate
extending to arbitrarily large starts, only references of exactly
the minimum slope can be useful among the low-slope references.

To compute the residue, read the bounded initial segment and
track its cost exactly. Continue only on edges with \(A=0\),
since (2) permits discarding all later positive positional
contributions at the chosen precision. Let \(h(u)\) be the minimum
additive cost from a vertex \(u\) to acceptance on this suffix
graph, including the empty path at acceptance. This is finite.
Replace each edge \(u\to v\) by the nonnegative cost
\(B+h(v)-h(u)\). Its bounded endpoint corrections allow the resulting
nondecreasing register to be truncated at a fixed threshold
without losing any residue contribution. Include the unit machines
of Section 1 and the valuation shift. The cost on each accepted
track gives a nonnegative integer exponent of \(\pi\).
Keep a vector over the finite coefficient ring counting provisional
states, and sum the weighted terminal outputs. Each input digit
updates that vector linearly, and the fixed padding is included
in the output operation. This is a finite output automaton.

To test a window of length \(w\), combine copies with digit
addition carries for \(n,n+1,\ldots,n+w-1\). The simultaneous
zeros at positive starts constitute a regular language, whose
emptiness is finite reachability. This extends the
[shifted residue automata](FACTORIAL_SHIFTED_MOMENT_AUTOMATA.md)
to reference choices with a certified uniform term bound, without
their special divided-power and constant assumptions.

## 4. An obstruction and an exact stable-perturbation test

The lower-bound condition can fail with a unique least slope.
Take \(p=2\), \(F=x^2y^2+2^D\) with \(D>4\) integral, and
\(H_n=(2n)!^2\). The term with one constant has ratio
\[
                      \frac{2^D n}{((2n)(2n-1))^2}.
\]
At \(n=2^\ell\) its valuation is \(D-\ell-2\).
Using the other pure term as reference also fails, since it has
strictly larger slope. The phase itself has positive ordinary
moments. This example detects failure of the two pure references.

The same obstruction excludes all critical-slope factorial references
for a class of supports, even if fixed negative offsets \(b_j\) are
allowed as long as the factorials are eventually defined. Let \(u\)
be a least-slope exponent column and let a convex combination \(v\)
of other support columns satisfy \(v_i<u_i\) in \(h\ge2\)
coordinates. A rational such combination exists by rational linear
feasibility. Choose an integer \(r_0\ge0\) large enough that
\(a_jr_0+b_j\ge0\) for every reference factor. Clear the combination's
denominators and scale its counts so that their sum \(K\) exceeds
\(r_0\) and \(K(u_i-v_i)>u_i r_0\) in each descending coordinate.

At \(n=p^\ell+r_0\) compare the term with those fixed other counts
and \(n-K\) copies of \(u\) against the pure \(u\) term.
Its multinomial contributes a valuation \(\ell+O(1)\) since the
falling factorial contains \(n-r_0=p^\ell\) exactly once and the
other-count factorial denominators are fixed. In a descending coordinate the lost factorial
interval contains \(u_i p^\ell\) exactly once. All its other
factors have valuations bounded independently of \(\ell\).
It therefore contributes \(-\ell+O(1)\). Increasing coordinates
give bounded changes because their added factors have positive
fixed offsets from \(u_i p^\ell\) when \(u_i\ne0\), and are
fixed integers when \(u_i=0\). Coefficient ratios are fixed
to the power \(K\). The relative valuation over \(p\) is
\((1-h)\ell+O(1)\).

The pure \(u\) term and every critical-slope reference differ
in valuation by \(O(1)\) along this sequence. Their linear
slopes agree and all digit sums in their factorials are bounded.
In fact, for sufficiently large \(\ell\),
\(a_j n+b_j=a_jp^\ell+(a_j r_0+b_j)\) has disjoint digit blocks.
Thus no fixed critical-slope product of linear factorials, with
a scalar exponential, gives a uniform lower term bound for this
geometry. This conclusion is independent of how small the
coefficient of \(v\) is.

For a fixed least-slope face and reference, the negative-cycle
test depends on the distinction between zero and positive
positional weights, independent of the sizes of the positive
gaps. Adjoin columns of strictly larger slope and restrict
acceptance to paths using at least one new column. If this
restricted graph has a reachable negative zero-positional cycle
with a zero-positional exit, increasing the new coefficients'
valuations cannot make every new term disappear uniformly at
any fixed precision. The path pump still applies at every
assignment of positive gaps.

Conversely, in the absence of such a cycle, increasing the
new valuations does suppress all the new terms uniformly to
any chosen precision. Let \(a_*>0\) bound every positive
positional weight below and \(a_{\rm new}\ge a_*\) be the least
gap assigned to a new column. At the last positive position
\(\ell\), any path using a new column has positional cost at least
\(a_{\rm new}+a_*(p^\ell-1)\). If a new digit occurs there,
its cost is at least \(a_{\rm new}p^\ell\). Otherwise
the earlier new digit and the last positive edge together
cost at least \(a_{\rm new}+a_*p^\ell\). The additive lower
bound in (2) still applies. The minimum over \(\ell\) of the
combined bound tends to infinity with \(a_{\rm new}\).
An all-order residue certificate therefore lifts from the
original phase to the enlarged support whenever this termwise
stability test passes and the new valuations are sufficiently large.

## 5. Which supports allow a uniformly dominant pure anchor?

Fix any exponent column \(u\) as an anchor and normalize by its
pure term. There exist coefficient valuations on the other columns
for which the relative term valuations are uniformly bounded below
if and only if both of the following support conditions hold.

1. No convex combination of the other columns is strictly lower than
   \(u\) in two coordinates.
2. In every coordinate where some column is strictly below \(u\),
   the entry \(u_i\) is a power of \(p\), with \(p^0=1\) allowed.

These are finite tests of support alone. The first is a linear
feasibility test for each coordinate pair. When both hold,
sufficiently increasing the other valuations makes the pure anchor
strictly dominate every moment, hence proves all-order nonvanishing.
One can use different descending coordinates for different mixed
terms. Over an arbitrary characteristic-zero field, the same proof
applies to the coefficient ratios whenever they admit the chosen
valuation and satisfy the bounds below. Here \(v_p\) denotes the extension normalized by \(v_p(p)=1\).

If the first condition fails, the fixed-count construction of
Section 4 gives unbounded negative valuation relative to the pure
term for any choice of nonzero coefficients, with no slope
assumption needed for that comparison. If the second fails, let
\(u_i=p^b a\) with \(a>1\) and \(p\nmid a\), and choose a column
\(v\) with \(v_i<u_i\). Take a fixed count \(K\) of \(v\) with
\(K(u_i-v_i)>p^b\), and positive orders \(n\ge K\) tending
\(p\)-adically to \(1/a\). The lost factorial interval in coordinate
\(i\) contains the factor \(u_i n-p^b\), whose valuation tends to
infinity. The multinomial valuation stays bounded since \(1/a\)
is not an integer. A non-descending coordinate contributes only
finitely many factors of the form \(u_j n+l\), \(l>0\), each
with a nonzero rational limit and bounded valuation. Other
descending coordinates only subtract valuations. Thus the
relative valuations are unbounded below for any coefficients.

For sufficiency put \(K=\sum_{v\ne u} k_v\) and consider a term
using the other columns with counts \(k_v\). By the convex
condition the total exponent can be below the pure exponent
in at most one coordinate. Where there is no decrease, all
factorial ratios are integers. If there is a decrease in
coordinate \(i\), its total entry is at least \(u_i(n-K)\).
By the prime-power factorial identity, the ratio in that row
has valuation at least
\[
       -\frac{u_i-1}{p-1}K - v_p\!\left(\frac{n!}{(n-K)!}\right).
\]
The last term is canceled by the falling factorial in the
multinomial. In the absence of a decrease one may likewise
drop that nonnegative multinomial contribution. Let \(B\) be
the maximum of \((u_i-1)/(p-1)\) over coordinates where any
column can be below \(u\), or zero if none exists. Using
\(v_p(k_v!)=(k_v-s_p(k_v))/(p-1)\), the relative valuation is
at least
\[
 \sum_{v\ne u} k_v
       \left(v_p(c_v/c_u)-B-\frac1{p-1}\right)
                  +\frac{\sum_{v\ne u}s_p(k_v)}{p-1}.
\]
Thus \(v_p(c_v/c_u)\ge B+1/(p-1)\) suffices, including equality.
For every mixed term the last sum is strictly positive.
The sole remaining pure anchor has valuation zero after
normalization, proving the assertion.

For planar supports, the first condition says the convex hull
has no point strictly southwest of the anchor. For example, at
\(p=2\), the anchor \((2,4)\) with other exponents \((0,7)\),
\((4,1)\) meets both tests. The two other columns decrease
different coordinates, and all columns lie on the same line
with positive normal \((3,2)\). This support has no
coordinatewise least column.

For unrestricted SFC(2) the remaining step is to supply references
or grouped normalizations giving window separation for arbitrary
coefficients and support. The negative-cycle criterion diagnoses
a failed termwise method and does not certify any exact zero.
No independent replay, formal verification or external review
is claimed.

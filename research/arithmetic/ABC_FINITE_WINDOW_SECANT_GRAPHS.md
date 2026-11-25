# Paired anchor secants concentrated at finitely many scales

Two anchor secants evaluated at the same primitive query are
independent linear forms on its projective direction. A full
divisor on an original query letter fixes their quotient to a
known unit residue. When most of each secant's numerical height
lies in a bounded number of narrow valuation windows, the
wide-cofactor estimate can count this quotient even if its
remaining prime factors have large radical.

## Finite-window secants form a sparse relation

Use the affordable fixed-record setting at height \(X\), with
\(n=X^L\) and \(L\) fixed before the height limit. The original
balanced endpoint truncations supply root caps \(O(L)\) in
base \(X\). Write the positive primitive triples as
\[
 U_j=F u_j,\qquad V_j=G v_j,\qquad
 T_j=\sigma F u_j+\tau G v_j=D w_j .
\]
The fixed positive full factors and signs are common. The
targets are distinct. Assume a fixed positive lower bound for
\(\max(\log_X F,\log_X G)\). Put
\[
 L_i(u,v)=u_i v-v_i u,\qquad
 H_{ij}=\frac{|L_i(u_j,v_j)|}
                  {D\gcd(u_i,u_j)\gcd(v_i,v_j)\gcd(w_i,w_j)} .
\]
These quotients are nonzero integers. Let \(J_i\) be a prime
pool fixed for each anchor, with prime product at most \(X^C\)
for fixed \(C\). Define \(R_{ij}\) as the part of \(H_{ij}\)
on primes outside \(J_iFGD U_iV_iT_i\). These primes also avoid
all three original query letters. Off the excluded pool their
valuations agree with those of the absolute raw secant.
[Quotient determinants and their original supports](ABC_AFFORDABLE_FIBRE_SECANT_VALUATION_WEIGHTS.md#primes-outside-the-endpoint-supports)

Fix \(r_0>0\), \(m\ge1\) and \(0<a<b<\infty\).
There are \(\eta_0,L_0>0\) such that the following holds
for all sufficiently late balanced endpoint boxes with
\(0<L\le L_0\). Let \(N\ge n^{r_0}\) be the number of
vertices. Join \(i\ne j\) when at most \(m\) intervals
in \([a,b]\), each of length at most \(\eta_0\), leave
numerical weight at most \(\eta_0\) uncovered in \(R_{ij}\).
Here a prime power \(p^e\Vert R_{ij}\) has scaled degree
\(Le\) and numerical weight \(e\log_X p\).
The intervals may vary with the edge. Then the number of
directed edges has the bound
\[
                    N^2 n^{-r_0/16+o_X(1)} .                 \tag{SW1}
\]
As usual an arbitrarily small fixed exponent allowance in
the height limit can replace \(o_X(1)\). The constants are
uniform under the fixed positive affordability and exponent
margins. The assertion uses no upper bound on the radical
of the uncovered part.

### A small-cost full divisor after fixing supports

The following observation applies to any subfamily of queries
with exponent bounded below by a fixed positive number.
One may strip and fix the parts of the three letters on a
fixed polynomial-height pool at inner subpower cost, and
dyadically refine the remaining indexed root values. Only a
fixed number of anatomic positions is present at each height
limit. Fixed-support enumeration is uniform over the pool.
Pass to an outer subsequence on which the measures of scaled
degrees of the remaining coordinates converge. Each coordinate
contributes its raw root height exponent divided by \(L\),
and the first moments are bounded by the numerical heights.
These measures are tight by the first-moment bound.

If a letter has positive nonatomic mass, choose a compact
continuity band at positive scaled degrees meeting that
mass. Remove neighborhoods of the finitely many large atoms
in the band with arbitrarily small loss of nonatomic mass.
For any prescribed small \(\mu>0\), take the large-atom
threshold small compared with \(\mu\). Then each retained
original coordinate has root mass less than \(\mu/4\) for
late budgets, by weak convergence. A greedy selection
gives a packet of normalized cost at most \(\mu\) and
positive limiting numerical weight. Its full powers supply
a modulus of fullness \(\gg 1/L\). This packet avoids the
prescribed pool. Freeze its value by pigeonholing its
indexed roots at cost at most \(\mu\) in normalized units.

If all three retained measures are purely atomic, the
query family has arbitrarily small exponent by the original
affordable modulus. To see this, retain a finite list of
shrinking positive atom windows and a shrinking zero window
on the two letters opposite whichever of \(F,G\) has a
positive numerical weight. Freeze the remaining positions
at arbitrarily small normalized root cost. Tightness ensures
that escaping moment, if present, lies in positions of
vanishing normalized mass. In a positive window use the
product of its original roots to a representative degree
and include the discrepancy in a wider correction. The
zero windows also have vanishing numerical height. A
positive-weight part of the original fixed factor remains
after dropping its prime powers of order below
\(\kappa/L\) for a sufficiently small fixed \(\kappa>0\).
The two original opposite values have a fixed unit ratio
modulo this full divisor. Their signed ratio has boundedly
many short roots and a wider factor of vanishing numerical
height. The amplified wide-cofactor estimate counts the
ratio at arbitrarily small normalized cost, and the
primitive values, hence the query, follow from the exact
ratio and the fixed factors. This contradicts the positive
lower exponent for the query subfamily.

Consequently a positive-exponent subfamily supplies a
packet on some original letter at arbitrarily small
normalized cost, with uniform positive weight and fullness
along a limiting subsequence. This uses the same atomic
alternative as the projective packing gap.
[Fixed-support estimate](ABC_EXCEPTIONAL_SET_AND_POWERFUL_TRIPLES.md#a-radical-count-contradiction-to-the-claimed-smooth-number-lower-bound)
· [Packing gap and its atomic case](ABC_AFFORDABLE_FIBRE_SECANT_VALUATION_WEIGHTS.md#setup-and-the-projective-packing-gap)
· [Wide-cofactor release](ABC_WIDE_COFACTOR_CONGRUENCES.md#a-bounded-number-of-small-roots)

### Proof of the sparse relation

Argue by contradiction through budgets \(L\to0\) and window
and tail tolerances \(\eta\to0\). If (SW1) fails, the power
mean for ordered pairs of distinct neighbors produces two
distinct fixed anchors \(i,k\) with a common-neighbor family
of normalized exponent at least \(r_0-2(r_0/16)\), up to
inner allowances. Indeed, for \(E\) edges and query degrees
\(d_j\),
\[
                \sum_j d_j(d_j-1)\ge E^2/N-E.
\]
Divide the lower bound by \(N^2\) to fix the anchor pair.
Their directions are distinct by primitivity and the
distinct-target assumption.

Use one fixed pool \(S\) containing both anchor pools, all
primes of the two anchor triples and \(FGD\), and the primes
dividing the determinants between any two of the two anchor
forms and the three original query forms. The latter forms
are
\[
               F u,\qquad G v,\qquad \sigma F u+\tau G v .
\]
All these determinants are nonzero: the two anchors have
different directions, each anchor is a primitive positive
query pair, and its three original letters are nonzero.
The pool has polynomial product height. Strip and fix the
parts of the original query letters and of both anchor
secants on \(S\), using inner subpower enumeration. Fix the
two signs of the secants as well. The resulting positive
parts \(R_i,R_k\) of the two raw secants are coprime because
any common prime of the values of two linear forms on a
primitive pair divides their coefficient determinant.
They are obtained from \(R_{ij},R_{kj}\) by dropping full
prime powers, so still satisfy the window and tail bounds.

Use the preceding small-cost packet on the original query
letters, taking its normalized cost below \(r_0/4\), and
fix its value \(P\). If its form is \(H\) from the list of
three displayed query forms, then \(H(u_j,v_j)=0\bmod P\).
Over each prime power in \(P\), the coefficient vector of
\(H\) is primitive and both anchor forms are units on the
resulting primitive direction. This uses the exclusion of
the coefficient determinants from the support of \(P\).
Equivalently, the ratio of the anchor forms modulo \(P\)
is the ratio of their coefficient determinants with \(H\),
with the same orientation in numerator and denominator.
The fixed signed parts supported on \(S\) are units at
\(P\). Thus \(R_i/R_k\) has a fixed unit residue modulo \(P\).

For each remaining secant and window let \(y\) be the
product of primes whose exact degree lies in that window.
Take an integer representative \(e_0\) with \(Le_0\) within
\(L\) of the window. At fixed \(L\) there are boundedly
many possible lists of representatives, since all windows
are in \([a,b]\). Fix the lists and partition window
overlaps deterministically. The factor contributed by a
window is
\[
                  y^{e_0}\prod_{p\mid y}p^{v_p(R_i)-e_0}
\]
on the first side, and with \(R_k\) on the second.
All these narrow roots have size \(X^{O_a(L)}\) because
\(Le\ge a\) and each original determinant has height
\(X^{O(1)}\). The total signed discrepancy in numerical
exponent is \(O_a(\eta+L)\). With the uncovered powers
in the wider correction, \(R_i/R_k\) is therefore a signed
monomial in at most \(2m\) narrow roots times a wider
rational of height \(X^{o(1)}\). Both the roots and the
correction are units modulo \(P\).

The exact value of \(R_i/R_k\) determines the query:
restore the two fixed signed factors on \(S\), and take
the ratio of two independent linear forms on \((u_j,v_j)\).
The directions are primitive and positive. To apply
amplification with controlled exact fibres, first color the
primes of each \(R_i,R_k\) by numerator or denominator,
retaining a fixed coloring at inner subpower cost. Each
query has only \(o_X(\log X)\) prime divisors in these
two factors, so averaging over binary colorings suffices.
Among this subfamily, a fixed \(\ell\)-fold product of
exact ratios has only divisor-many ordered decompositions
into the reduced original ratios, for each fixed \(\ell\).
There is no cancellation between the two colored supports.
Each ratio in such a decomposition specifies the original
query.

Apply (GW) to the componentwise products of \(\ell\)
tuples with residue target raised to the power \(\ell\).
Their output multiplicities and exact fibres cost only an
inner subpower factor by the preceding decomposition.
Their short heights still have exponent \(O_{\ell,a}(L)\)
and the wider rational has numerical height \(o(1)\)
for each fixed \(\ell\). The exponent from (GW) before
dividing by \(\ell\) is \(O(L)\) with constant depending
on \(m\) and the positive numerical and fullness margins
of \(P\), but independent of the fixed amplification
order. Consequently the count after the packet lock
has arbitrarily small normalized exponent, upon taking
\(\ell\) sufficiently large before the outer limit.
The initial common-neighbor count after fixing the
packet is bounded below by a positive normalized
exponent. This contradiction proves (SW1).


### Allowing a bounded radical defect

Suppose an edge becomes eligible for the window criterion
after removing a full divisor of \(R_{ij}\) of radical at
most \(n^\beta\), for some fixed \(\beta\ge0\).
For every fixed \(\zeta>0\), sufficiently small window and
tail tolerances now give the edge bound
\[
                    O(N^{3/2}n^{\beta+\zeta}+N) .               \tag{SW2}
\]
In particular this excludes a near-full pair exponent
whenever \(\beta<r_0/2\), with a fixed strict margin.

For an anchor pair as above, first intersect the common
neighborhood with the \(S\)-support and sign records.
Choose one defect divisor for each of the two secants and
strip its \(S\)-part if necessary. For each secant its
radical has at most \(n^\beta\) choices, and a specified
support has only \(X^{o_X(1)}\) integer values of polynomial
height. Thus fixing both exact defect factors costs at
most \(2\beta+o_X(1)\) in base \(n\).
Enlarge \(S\) by their supports. If any remaining powers
in the other anchor form lie on this support, they would
divide both raw forms outside the original \(S\), which
is impossible by the pair determinant. In all cases
further pool-supported parts can be fixed at inner
subpower cost and only decrease the remaining secant
numerical weights. The two fixed defect factors are
now among the factors restored in the exact-ratio
argument. If every height tolerance tending to zero
admitted a common neighborhood of exponent at least
\(2\beta+\zeta\), the refined neighborhood would have
a fixed positive exponent after the two defect records.
Run the small-cost packet and the amplified ratio count
with this new lower exponent. The resulting contradiction
gives uniformly at most \(n^{2\beta+\zeta+o_X(1)}\)
common queries for two distinct anchors at small enough
budgets and tolerances. Substitution in
\(\sum_j d_j(d_j-1)\ge E^2/N-E\), with a smaller
fixed inner allowance if needed, gives (SW2).




### Two prescribed linear forms

The common-neighbor argument also applies to fixed coefficient
vectors unrelated to the original rows. Let
\(\Lambda_1,\Lambda_2\) be integer binary linear forms of
polynomial coefficient height, with nonzero coefficient
determinant. Require each \(\Lambda_i\) to be nonproportional
to each of \(Fu\), \(Gv\), and \(\sigma Fu+\tau Gv\).
A prescribed extra prime pool of polynomial product height
may vary with the two forms but is fixed across the query
family. For a query, let \(R_i\) be the part of the nonzero
absolute value \(|\Lambda_i(u,v)|\) outside that pool.
Suppose for each \(i\) there is a full defect divisor of
\(R_i\) with radical at most \(n^{\beta_i}\), leaving the
same finite-window and numerical-tail criterion as above.
For every fixed \(\zeta>0\), at sufficiently small budgets
and sufficiently small window and tail tolerances the
number of such queries is
\[
                          O(n^{\beta_1+\beta_2+\zeta}).           \tag{SW5}
\]
Here \(m,a,b\), the coefficient and pool height caps, and
the affordability margin are fixed uniformly, and the
arbitrary inner height losses can be absorbed into
\(\zeta\).

Enlarge the pool by the primes of every nonzero
determinant between the two prescribed forms and the
three original forms, including the determinant between
\(\Lambda_1,\Lambda_2\). The three nonproportional
conditions for each prescribed form ensure that a cheap
full modulus on any original query letter makes both
prescribed values units and specifies their ratio modulo
that modulus. Outside the pool the two absolute values
are coprime. Fix signs and pool parts at inner subpower
cost, then record the two exact defect factors at
normalized cost at most \(\beta_1+\beta_2+o_X(1)\).
Their supports can be added to the pool. If (SW5)
failed through arbitrarily small tolerances and budgets,
there would be a positive-exponent query family after
these records. The small-cost packet observation and the
amplified ratio count then contradict this remaining
exponent. Independence of \(\Lambda_1,\Lambda_2\)
makes the exact restored ratio injective on primitive
positive query directions. This is the same use of the
actual original additive equation as in the anchor proof.

### Adaptive exact-power defects of high order

After the simultaneous real and projective peeling, include
the common projective modulus \(M\) in each anchor pool.
Use the consequence of (PS6) giving a fixed upper cutoff
\(\kappa\) for all but an arbitrarily small numerical weight
of the powers outside the pools, on all but a power-small
fraction of pairs. The cutoff is independent of the desired
small tolerance for the discarded numerical weight. Then for
each fixed \(a>0\) and \(\xi>0\), all but a power-small fraction
of the retained pairs satisfy
\[
 \prod_{\substack{p^e\Vert R_{ij}\\ q\nmid e}}p
                        > n^{\,r_0/2-\xi}
                      \quad\hbox{for every integer }q\ge a/L,     \tag{SW4}
\]
where \(r_0\) is a fixed lower bound on the exponent of
the family before peeling. If the right-hand exponent is
nonpositive the nontrivial claim concerns only positive
smaller margins. The integer \(q\) may be chosen separately
on every pair.

Indeed take the exponent loss in peeling smaller than
\(\xi/2\) and apply (SW2) to the resulting fibre with
\(\beta=r_0/2-\xi\). For a violating pair remove the full
powers at \(q\)-defect primes as the defect divisor.
Below the fixed cutoff, the remaining scaled degrees
are in the grid \(Lq,2Lq,\ldots\) and there are at most
\(\lfloor\kappa/a\rfloor\) of them in the permitted range.
They occupy singleton windows. Set the tolerance for
the high-order numerical tail using (SW2) with that
fixed number of windows. This is possible because the
cutoff was selected independently of the tail tolerance.
The remaining tail is only smaller after deleting the
defects. Equation (SW2) now forces the power saving,
without a union bound over the varying orders \(q\).
The same reasoning handles a union of a fixed number
of such grids with spacings bounded below by \(a\).
[Uniform high-order cutoff after peeling](ABC_AFFORDABLE_FIBRE_PROJECTIVE_SECANT_PEELING.md#a-consequence-without-any-secant-power-hypothesis)


### Adaptive defects above a fixed order

For a single exact-power order per edge, the lower cutoff
\(q\ge a/L\) can be replaced by a sufficiently large fixed
integer \(Q_0\). Work in one common dyadic box for both
cofactor coordinates. For every fixed \(0<\xi<r_0/2\)
there are \(Q_0,L_0,\delta>0\) such that, for \(L\le L_0\),
outside at most \(N^2 n^{-\delta+o_X(1)}\) ordered pairs,
\[
 \prod_{\substack{p^e\Vert R_{ij}\\ q\nmid e}} p
                  > n^{\,r_0/2-\xi}
                  \qquad\text{for every integer }q\ge Q_0 .     \tag{SW6}
\]
This uses no numerical upper-order cutoff on the secants.
The anchor pools \(J_i\) can vary as before.

First screen pairs for which \(\operatorname{rad}(R_{ij})
\le n^{3r_0/4}\). For a fixed anchor, the raw determinant
outside its fixed pool \(J_iFGD U_iV_iT_i\) agrees with
\(R_{ij}\) up to sign. Enumerate the outside radical at
cost at most \(n^{3r_0/4}\). For each support, the raw
determinant has \(X^{o_X(1)}\) values by fixed-support
enumeration, uniformly in the pool and support because
their total prime product has polynomial height. For
each signed determinant, the query pairs lie on an
integer affine line with primitive difference direction
\((u_i,v_i)\), hence there are \(O(1)\) of them in the
common dyadic box. This screens a power-small fraction
uniformly in the anchor-dependent pools.
[Fixed-anchor radical count in dyadic cofactor boxes](ABC_SECANT_VALUATION_DEFECT_GRAPHS.md#sparse-pair-bounds-and-adaptive-orders)

Put \(\beta=r_0/2-\xi\). On a remaining violating pair,
the primes whose orders are divisible by the chosen \(q\)
have radical at least \(n^{3r_0/4-\beta}\). Their full
powers have height at most \(X^2\). Therefore
\[
                         q\le \frac{2}{L(3r_0/4-\beta)} .
\]
In particular at fixed \(L\) there are only finitely many
orders in the violated range before the height limit.
This permits choosing one violated order per edge and
pigeonholing the same \(q\) over the whole retained graph.
Taking two anchors by the previous neighbor count now
gives a common-query family with this *same order* for
both secants. Choose \(\delta\) smaller than
\(\min(r_0/16,\xi/8)\). The common-query exponent stays
above \(2\beta\) with a fixed margin.

In this common neighborhood, fix signs and pool factors,
then the full powers on both \(q\)-defect supports. This
costs at most \(2\beta+o_X(1)\) as in (SW2), leaving a
positive-exponent family. Outside the enlarged fixed
pool, both remaining secants are exact \(q\)-th powers,
say \(x^q\) and \(y^q\), with \(x,y\le X^{2/q}\). By
the small-cost packet observation, some original query
letter supplies a full modulus \(P\) of fixed positive
numerical weight on a limiting subsequence, at a cost
smaller than the margin in the query exponent. Fix that
packet. The quotient \((x/y)^q\) has a fixed unit residue
modulo \(P\). At the fixed \(L\) and chosen \(q\), there
are \(X^{o_X(1)}\) unit \(q\)-th roots of this residue:
the kernel has at most \(q\) elements on each odd
prime-power factor and at most \(q^2\) on the single
factor at \(2\). For \(q\) sufficiently large relative to
the positive weight of \(P\), a unit residue determines
\(x/y\) in the indicated height range, since any
nonzero cross difference of two such rationals has
size at most \(2X^{4/q}<P\). The exact quotient of
the two secant forms then determines the query.

If no fixed \(Q_0,L_0\) with (SW6) existed, choose
violations with \(L\to0\) and \(Q_0\to\infty\).
The same-order selection gives \(q\to\infty\) on that
sequence. The fixed positive modulus weight on the
limiting subsequence then makes the reconstruction
argument applicable, contradicting the common-neighbor
exponent. This proves (SW6). The finite-window argument
above additionally accommodates unions of separated
grids without a shared high order.

## A lower bound for diffuse secant cost

Associate to each ordered pair its numerical secant measure
on the compactified degree range \([0,\infty]\):
\[
   \mu_{ij}=\sum_{p^e\Vert R_{ij}}
                    e\log_X p\,\delta_{Le}.
\]
Its mass is bounded by the numerical height of the raw secant.
For \(r_0\) as in (SW1), let \(\mathbb P\) be any iterated
weak limit of the laws of \(\mu_{ij}\) under uniform distinct
ordered pairs. In the height limit the budget is fixed.
One may then take a subsequence of the resulting laws as
\(L\to0\). The space of finite measures of bounded mass
on \([0,\infty]\), and the space of laws on it, are compact
and metrizable. For \(\mathbb P\)-almost every \(\mu\) satisfying
\[
 \mu(\{0,\infty\})=0,\qquad
                 \int_{(0,\infty)}\frac{d\mu(s)}{s}<\infty ,
\]
the nonatomic part \(\mu_{\rm na}\) must have normalized
root cost
\[
                  \int_{(0,\infty)}\frac{d\mu_{\rm na}(s)}{s}
                         \ \ge\ r_0/2.                          \tag{SW3}
\]
If \(\log_n N\to r>0\) in the respective limits, \(r_0\)
may be replaced by \(r\), by taking \(r_0\uparrow r\).

In fact, for a measure satisfying the two hypotheses
and violating (SW3), choose \(\beta<r_0/2\) strictly above
its nonatomic cost. By integrability a finite list of
positive atoms accounts for all but strictly less than
\(\beta\) of its total root cost. Choose an interval
\([a,b]\subset(0,\infty)\) containing those atoms in its
interior, and take the window and tail tolerance from
(SW2) with \(\zeta>0\) so small that
\(\beta+\zeta<r_0/2\). Place windows at the selected
atoms with lengths below the prescribed tolerance and
with endpoints of zero \(\mu\)-mass. Now choose a compact
interval \([a',b']\) containing \([a,b]\), with endpoints
of zero \(\mu\)-mass, so that the numerical mass
outside \([a',b']\) is strictly below the tail tolerance.
The latter interval is chosen after the tolerance. Its
endpoints are used only to choose the discarded powers,
and do not change the short-root bounds for the
original windows.

Remove every prime power of the queried secant in
\([a',b']\) outside the selected windows. Its radical
cost is the integral of \(s^{-1}\) on this fixed
complement of windows. For \(\mu\) the cost is strictly
less than \(\beta\), because every selected atom has
been excluded from the integral. Both strict bounds
persist on a weak neighborhood of \(\mu\). Indeed
the integrands on \([a',b']\) are bounded and their
discontinuities are at endpoints where \(\mu\) has
zero mass. All such pairs meet (SW2), so the probability
of that neighborhood is zero in the iterated limit.
The space is second countable. A countable subcover
of these zero-probability neighborhoods covers all
measures violating (SW3) under its hypotheses.

When the pools also contain the common projective modulus
and the high-order tail has been eliminated by the peeling
bound, one can restrict to limits with no mass at infinity.
An excluded purely atomic law on positive degrees can then
survive only by numerical mass tending to scaled degree zero,
or by divergence of its normalized root cost toward zero.
This separates a diffuse positive-degree obstruction from
large radical appearing in very low scaled orders.
[Full-prime weight after projective peeling](ABC_AFFORDABLE_FIBRE_PROJECTIVE_SECANT_PEELING.md#a-consequence-without-any-secant-power-hypothesis)

## Consequence for the large-radical regime

The criterion includes a bounded number of variable exact-power
layers and windows containing arbitrarily many distinct exact
orders. For example, two such windows can account for the
numerical height even if the tail is squarefree of height
\(X^{o(1)}\) with radical exceeding every fixed power of \(n\) in
the outer limit. The exact-order defect count need not be useful
on that sector. The argument instead uses the quotient of two
secant values and a full divisor on the original query.

After simultaneous real and projective peeling, one may
use anchor pools containing the common projective modulus
and use the bound putting almost all the retained
numerical weight into bounded scaled orders. For every
fixed \(m\) and compact positive scaled-degree range, a
surviving pair family must leave a fixed positive
numerical weight outside any \(m\) sufficiently small
windows in that range. This permits a substantial
contribution from orders near scaled degree zero or
numerical weight spread over sufficiently many positive
scales. Controlling those configurations remains necessary
for an unrestricted small-budget improvement.
[Bounded-order numerical weight after peeling](ABC_AFFORDABLE_FIBRE_PROJECTIVE_SECANT_PEELING.md#a-consequence-without-any-secant-power-hypothesis)
· [Exact-order defect counts](ABC_SECANT_VALUATION_DEFECT_GRAPHS.md)

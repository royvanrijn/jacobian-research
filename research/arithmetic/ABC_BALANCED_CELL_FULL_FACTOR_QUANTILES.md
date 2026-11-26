# Full-factor quantiles in balanced residue cells

A typical balanced residue cell gives many common queries to
several distinct anchors. The same original letter occurs as
an affine function of each query, with a different offset at
each anchor. This lets full powers from repeated copies of that
letter contribute to one numerical congruence. At positive
multiplicity, sufficiently large compact valuation tails are
therefore limited by the cost of prime overlap with the anchor
differences.

## A uniform count for packets on repeated copies of one letter

Put \(n=X^L\), with the height limit taken first at each fixed
\(L\). Work with actual positive primitive \(A+B=C\),
\(X/2<C\le X\). All three original radical caps remain in
force and are bounded in base \(n\). Fix a computed role
\(W\), an anchor role \(U\), and a query role \(V\), all
different. The complete \(U,V\) values specify the row.
There are common signs such that, at fixed \(U=b\),
\[
                        W=\epsilon V+\epsilon' b ,
          \qquad \epsilon,\epsilon'\in\{-1,1\}.
\]
Suppose that the three original radical caps are
\(n^{r_T+\tau}\), the complete \(U,V\) dictionaries both
have size at most \(n^{u+\tau}\), and \(V\le X^{h_V+\tau}\),
where the parameters \(r_T,u,h_V\) are fixed, \(u>0\).
Write
\(\beta=r_W-2u\). A negative \(\beta\) with a fixed margin
already excludes an edge exponent approaching \(2u\)
by the original fixed-letter estimate. We use \(\beta\ge0\).

The queries may lie in one fixed residue class modulo
\(G\le X^K\), where \(G\ge X^{g-\tau}\) for fixed
\(0\le g<h_V\). Write \(H=h_V-g\). The default is
\(G=1,\ g=0\). The modulus \(G\) can vary with height,
and its primes may have any support. A prescribed
interval of length \(X^{h_V+\tau}\) for \(V\) also
suffices in place of the stated numerical cap.

Suppose each row has a packet of full prime powers on \(W\)
with total normalized radical at least \(\beta+\rho\),
where \(\rho>0\) is fixed. Every selected full exponent
\(e\) satisfies \(Le\ge s>0\), and each selected prime has
normalized cost \(\log_n p\le b_0\), with \(s,b_0\) fixed.
The packet may vary with the row and need not be an
admissible grouped dictionary. Only the selected original
full powers will be recorded.

**CQ1.** If a fixed \(t>0\) satisfies
\[
                      t+b_0<u,\qquad st>H ,                 \tag{CQ1}
\]
then there is a fixed \(\delta>0\) such that the row count
is \(O_{L,\varepsilon}(n^{2u-\delta}X^\varepsilon)\) for
sufficiently small fixed \(\tau,L>0\).
The threshold and the saving may depend on the displayed
fixed bounds and their strict margins. The count is
uniform in the integer progression and in the actual
rowwise packet choices.

Choose \(m\) fixed with \(m\rho/2>t\). If the edge count
is at least \(n^{2u-\chi}\) with a small fixed \(\chi>0\),
the star count on the \(U,V\) graph gives a list of
\(m\) distinct complete \(U\) anchors having at least
\[
                     n^{u-m\chi-(2m-1)\tau-o_X(1)}
                                                               \tag{CQ2}
\]
common queries. Indeed a query of degree \(d\ge2m\)
gives \(\gg_m d^m\) ordered distinct anchor lists.
The low-degree queries can account for only
\(O_m(n^{u+\tau})\) edges. Divide the power-mean lower
bound for the lists by the number of possible ordered
anchor lists. This counts genuine original rows, with
the packet property at each of the \(m\) anchors.
[Common-column star count (TM1)](ABC_TWO_MINIMUM_DICTIONARIES_TAIL_PACKING.md#prime-power-packing-with-two-dictionary-caps)

For the chosen anchors \(b_1,\ldots,b_m\), put into a
pool \(P\) all primes of \(G\) and of the nonzero differences
\(b_i-b_j\). The support has a polynomial-height product
uniformly for fixed \(m,K\). At any fixed anchor the number
of queries carrying \(W\)-radical mass at least
\(\beta+\rho/2\) on \(P\) is bounded by
\[
 O\!\left(n^{(r_W-\beta-\rho/2)/2+O(\tau)+\eta}
                     X^\varepsilon\right)
       =O\!\left(n^{u-\rho/4+O(\tau)+\eta}X^\varepsilon\right)
\]
by (FL1), with arbitrary small fixed \(\eta>0\).
The event is empty if the corresponding outside cap is
negative. Choose the star and counting losses smaller
than this fixed margin. Outside \(P\), at least \(\rho/2\)
mass remains in the packet on each affine form \(W_i(V)\).
The outside supports on distinct \(W_i\)'s are disjoint
at each query because their differences are nonzero
fixed integers with prime support in \(P\).
[Fixed-letter bounds with a prescribed pool (FL1)](ABC_FIXED_LETTER_RADICAL_CAPS.md#a-pointwise-count-with-prescribed-prime-pools)

Greedily select full powers from those remaining packets
until their combined radical cost first reaches \(t\).
The overshoot is at most \(b_0\). Their product has
numerical exponent at least \(st\). Write \(V=v_0+Gz\).
On each selected power the affine form has unit slope
in \(z\), so the power determines one residue for \(z\).
The product of these coprime powers exceeds the
parameter interval length \(O(X^{H+O(\tau)})\). For any
record of the selected divisors there is at most one
query in the interval. Their total radical costs at
most \(t+b_0\); enumerate the squarefree prime support,
its allocations to the fixed anchors and the supported
powers at cost \(n^{t+b_0+o_X(1)}\). Both the power values
and the number of anchor allocations have uniform
height-subpower costs on a fixed support. This is less
than (CQ2) by a fixed margin for small enough
\(\chi,\tau\), proving (CQ1).

## Compact high-tail constraints at any positive multiplicity

Use a regular balanced full-exponent cohort as in
(RB1)–(RB4), with common degree exponent \(u>0\).
Pass to an anatomic profile in which \(\nu_W\) is the
limiting exact scaled valuation measure of \(W\),
with total mass \(1/3\). Its numerical height is \(h_W\).
Choose any computed role \(W\) and one of its opposite
roles \(V\) as the query, of numerical height \(h_V\).
A cheap full divisor on \(W\) produces a typical cell
with both opposite complete dictionaries
\(n^{u+o(1)}\) and \(n^{2u-o(1)}\) actual rows. Such a
divisor is available on every role when \(u>0\) by
(FL2), applied in the corresponding original positive
complete fibres. The same argument applies to a cell
which also has a prescribed query progression of
numerical weight \(g\) as above and still has these
two dictionary exponents and this row exponent.
Put \(H=h_V-g\) there, and \(H=h_V\) without the
extra progression.
[Typical balanced cells (RB4)](ABC_BALANCED_DICTIONARY_RESIDUE_CELLS.md#the-actual-residue-cell-sizes)
· [Nonatomic mass in a positive-exponent fibre (FL2)](ABC_FIXED_LETTER_RADICAL_CAPS.md#strict-positive-scale-inequalities-on-the-two-free-letters)

**CQ3.** The following restrictions hold for this computed
and query choice:
\[
 \begin{split}
   (\nu_W)_{\rm na}\big((H/u,\infty)\big)
                              &\le \beta,\
   \nu_W\!\left(
       \left(\frac{\min\{h_W+H,\ \max(h_W,2H)\}}{u},
                       \infty\right)\right)
                              &\le \beta,
                    \qquad \beta=\tfrac13-2u .                \tag{CQ3}
 \end{split}
\]
The intervals here concern the weak limiting measure on
finite scaled degrees. Each restriction may be applied
with either opposite role as query when its own typical
cell has the stated data.

For the first bound, a violating tail gives a compact
continuity region lying above some fixed \(s>H/u\)
with nonatomic mass \(>\beta\). Remove sufficiently
short neighbourhoods of the finitely many atoms of
mass at least a prescribed small threshold. The
nonatomic mass removed can be arbitrarily small.
After weak convergence, every individual prime in
the retained finite union of compact continuity
intervals has small normalized cost. Otherwise a
sequence of such primes would create an atom of at
least that cost at a limiting scaled degree in the
retained compact set. The mass of the retained
packet exceeds \(\beta\) by a fixed margin. Choose
\(t<u\) with \(st>H\), then the atom threshold
small enough for \(t+b_0<u\), and apply (CQ1)
on a typical cell. The packet is supplied by
the original exact valuation positions.

For the first cutoff on total mass, choose a compact
positive band of mass \(>\beta\) at lower order
\(s>(h_W+H)/u\). Every prime in it has cost at most
\((h_W+o(1))/s\), since its full power divides \(W\).
Thus one can again satisfy the two inequalities
in (CQ1). For the other cutoff let
\(s>\max(h_W,2H)/u\).
A full prime power in this band still has cost
less than \(u\) by a fixed margin. A query for
which such a power has numerical weight
greater than \(H+\xi\), for sufficiently small fixed
\(\xi>0\), is reconstructible on one anchor
unless the prime divides \(G\).
To handle that support, first discard queries
having pool overlap at least \(\beta+\rho/2\)
as in (CQ1), and use only powers outside the pool.
For any remaining query with such an outside
power the one-power count is at most
\(n^{(h_W+o(1))/s}\), below the common-query
exponent by a fixed margin. Summing over the fixed
anchor list preserves this loss. On every surviving
query the outside packet still has positive fixed
mass, and its individual prime costs are now at most
\((H+\xi)/s\). Greedy selection and CRT finish because
one can take \(t>H/s\) with
\(t+(H+\xi)/s<u\). This proves the second total-mass
cutoff. Use continuity regions to pass from the bands
to the open tails in (CQ3).

For \(u=1/6\), the trimmed cap \(\beta\) is zero,
and these are the compact support cutoffs from the
minimum-dictionary argument. For intermediate \(u\)
the global complete dictionaries have exponent
\(1/3-u>1/6\), and (RB4) supplies the smaller local
dictionaries needed for the test. The discarded pool
can carry mass \(\beta\), which is why the tail
bounds have that term.
[Earlier compact-support cutoffs at the minimum scale (TM3)](ABC_TWO_MINIMUM_DICTIONARIES_TAIL_PACKING.md#prime-power-packing-with-two-dictionary-caps)

The compact-band proof is at the level of actual prime
powers. For instance, with \(\beta=0\), a prescribed
row sector with a prime of cost bounded below by a
positive constant and scaled degree above the
appropriate finite cutoff by a fixed margin is
excluded at saturation. The prime's scaled degree
is bounded above at such a cost by the height. No
claim is made that degrees form a finite discrete
set, or that a tail of vanishing radical mass has
vanishing numerical weight.

## A two-copy numerical quantile

Take \(G=1\) and retain a compact positive interval
\([a,b]\) of the computed role's degree scale.
Write
\(\mu=(\nu_W)_{\rm na}|_{[a,b]}\).
Let \(F_\mu(t)\) be the largest first moment from
a submeasure of \(\mu\) of mass at most \(t\).
For fixed \(t_1,t_2\ge0\) with \(t_1+t_2<u\), the
same two-anchor argument gives
\[
        F_\mu(t_1)+F_\mu(t_2)-b\beta\le h_V .                 \tag{CQ4}
\]
For example \(t_1=t_2\) tests the pooled top quantile
of two copies of the same original spectrum. The
right side is the numerical width of the actual
query letter, which may be shorter than the
computed letter.

For the proof choose two distinct anchors with a
common-query exponent approaching \(u\) by (CQ2).
Let \(P\) contain the prime support of their nonzero
difference. By (FL1) discard queries on which either
computed term has \(P\)-supported radical mass at
least \(\beta+\eta\), using fixed \(\eta>0\) and then
small chart losses. Fix the exact \(P\)-supported
parts \(J_1,J_2\) at height-subpower cost. Each term
has \(P\)-part of polynomial height on the prescribed
support.

Before dropping the \(P\)-supported primes, each
computed value supplies a compact-order divisor
\(D_i\) of cost at most \(t_i+o(1)\) and numerical
weight at least \(F_\mu(t_i)-o(1)\), with arbitrary
fixed outer allowances. This follows by approximating
a top quantile by finitely many continuity bands,
discarding neighbourhoods of the finitely many
large atoms, and selecting entire small prime
powers from each band. One may approximate the
desired masses from below and use sufficiently
small prime-cost cutoffs to absorb rounding.
The lost \(P\)-part of \(D_i\) has numerical weight
\(l_i\le b(\beta+\eta)\) in base \(X\). The two
outside parts are coprime. Their product together
with \(\operatorname{lcm}(J_1,J_2)\) therefore
supplies numerical weight at least
\[
   F_\mu(t_1)+F_\mu(t_2)-l_1-l_2+\max(l_1,l_2)-o(1)
           \ge F_\mu(t_1)+F_\mu(t_2)-b(\beta+\eta)-o(1).
\]
The fixed supported parts recover at least the
larger of the two discarded contributions. The
outside records cost at most \(t_1+t_2+o(1)\)
in base \(n\). If the strict reverse of (CQ4)
held, take \(\eta\) and the other margins small.
The combined congruences would determine at
most one query per record, contradicting
\(t_1+t_2<u\) and the common-query exponent.
This proves (CQ4).

A fixed small positive numerical relaxation is
also available with (TC4). All queries are units
modulo their combined record by original
primitivity, and the step has polynomial height
and bounded normalized radical cap. If \(h_V\)
is bounded below and the recorded weight is
within a sufficiently small positive margin of
\(h_V\), the resulting progression has an
arbitrarily small normalized occupancy exponent.
Use an allowance smaller than \(u-t_1-t_2\).
[Unit progression counts at a large small-radical step (TC4)](ABC_TARGET_PACKET_CRT_MOMENTS.md#unit-progression-cells-with-a-large-small-radical-step)

The multi-anchor inequality can retain numerical weight
on the fixed prime pool, add a query quantile, and price
the marked divisor through the actual residue classes.
Its finite-copy and limiting forms are in PQ1–PQ5.
[Full powers on the prescribed pool](ABC_RECOVERED_POOL_QUANTILE_STARS.md#pq1-fixed-star-hypotheses-and-the-quantile-test)

## Coverage of the numerical test

For one unequal-height consequence take \(u=1/6\),
a tall computed role \(h_W=1\), and a wholly
nonatomic \(W\)-measure. With \(G=1\), (CQ3) confines
that measure to \([0,6h_V]\).
If its first moment is the full numerical height,
a query height \(h_V<1/2\) is impossible, and
\(h_V=1/2\) would force a point mass at scaled degree
three, contradicting the nonatomic hypothesis.

There is a weaker consequence permitting escaped
first moments. Pigeonholing the original anchor
role in the regular cohort gives a fixed-anchor
fibre of exponent \(1/6\). Its two moving roles are
\(V,W\), so (FW2) gives \(e_W+e_V+\theta\le h_V\)
for a fixed positive \(\theta\); here \(e_T\) is
the escaped first moment on role \(T\).
On the other hand (CQ3) bounds the finite first
moment on \(W\) by \(2h_V\). Hence
\(1=e_W+\int r\,d\nu_W(r)\le3h_V-\theta\).
This excludes \(h_V\le1/3\) in the same wholly
nonatomic limiting profile, without assuming
first-moment convergence.
[Escaped numerical weights in a complete fibre (FW2)](ABC_FIXED_FIBRE_VALUATION_WIDTH.md#full-divisor-locks-leave-a-uniform-amount-of-numerical-width)

These arguments apply for every fixed positive
\(u\). At \(u=0\), even a matching can have the
complete dictionary and cell exponents of (RB4);
the star construction gives no positive common-query
exponent there.

The tests leave tall diffuse candidates at the
level of the limiting inequalities. For example,
take all heights one and all three formal
degree measures equal to
\(\tfrac13\operatorname{Unif}[2,4]\).
With no extra query progression, every cutoff in
(CQ3) is above the support for \(0<u\le1/6\).
The left side of (CQ4) is at most
\(4(t_1+t_2)\le2/3\). This example describes
formal inputs to these particular inequalities.
It does not assert an actual primitive additive
family or compatibility with every other
endpoint constraint.

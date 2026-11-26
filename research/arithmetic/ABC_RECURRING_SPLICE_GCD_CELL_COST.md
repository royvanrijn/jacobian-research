# Recurring splice divisors in residue and log cells

A large average splice gcd in a single cell forces a
fixed divisor of the marked coordinate on many rows.
Across many cells the divisor catalogue has a numerical
cost. The adjustable list exponent of the mixed
real and residue moments pays this cost strictly below
the number of original cells. Both the arithmetic
divisibility and the list exponent are needed.

## Power rarity of cells with large mean extra gcd

Put \(n=X^L\), with the budget \(L>0\) fixed before
taking \(X\) to infinity. Let \(d=1/3\), and consider
a chart of positive primitive rows \(A+B=C\),
with \(C\le X\). Choose distinct coordinate roles
\(W,Z\), use \(W\) as the marked role, and let
\(M\le X\) be a common divisor of all marked values.
The primes in the support of \(W\) are disjoint from
the primes in the support of \(Z\) throughout the
chart. It suffices that every compared value and
its root representation are units at every
divisor of a marked value. Global role colouring
gives this situation at height-subpower cost.

The compared values lie in one exact valuation box
of the type in JRM1, with total root cap \(d+O(\rho)\)
and maximum slot cap at most \(\rho\). Fix a
prime pool \(P\) of polynomial height containing
the primes of \(M\). Every prime \(p\) dividing
a marked value outside \(P\) has \(\log_n p\le\rho\).
The maximum number of chart rows with any fixed
compared value \(Z\) is at most \(n^{u+\rho}\),
with \(u\ge0\).
Here \(\rho>0\) is a chart tolerance. Place each
row in its cell \((Z\bmod M,\ \text{log bin of }Z)\),
with logarithmic bins of width \(\asymp n^{-\xi}\),
where \(\xi\ge0\).
Fix \(s\in[d/2,d)\) and assume
\(\log_n M+\xi=s+O(\rho)\). The constants implicit
in this description are fixed. The definition also
covers pure residue cells by subdividing the
logarithmic range into \(O(\log X)\) bins of
bounded width.

Write \(\gamma=d-s\) and call a cell full for this
statement if its row count \(N\) satisfies
\[
                 n^{d-s-\rho}\le N\le n^{d-s+\rho}.
\]
For ordered distinct rows \(i,j\) of a cell put
\[
         D_{ij}=\frac{\gcd(W_j,Z_i-Z_j)}{M},
           \qquad D_{ii}=1.
\]
The quotient is integral. With \(W=C\) and \(Z=A\),
the numerator is also \(\gcd(C_j,A_i+B_j)\).
For other choices of \(W,Z\), it is the gcd of
the marked value with the corresponding signed
splice.

**RSC1.** Fix \(\eta>0\), put
\(t=\min(\eta/16,\gamma/16)\), and choose
\(0<b<\min(\eta/8,\gamma/4)\).
Assume that
\[
 \Delta=\min(t,\gamma/4)-2u>0,\qquad
 \kappa=\frac{b(d-s-b)-u(s+b)}{s+b-d/2}>0 .
\]
For every \(0<\kappa_0<\min(\kappa,\Delta/2)\),
if the chart tolerance is sufficiently small
in terms of the fixed parameters, then the number
of full cells with
\[
                       \mathbb E_{i,j}\log_n D_{ij}\ge\eta
\]
is \(O(n^{s-\kappa_0})\), for independent uniform
rows \(i,j\) of a cell. The constants can depend on
the fixed budget and exact box data. For \(u=0\)
the gap is positive at every fixed \(d/2\le s<d\).
The conclusion concerns full cells across the
numerical range \(1\le Z,W\le X\). For \(s<d/2\)
the direct single-anchor bound
below gives the corresponding result.
[Joint real and residue list bounds (JRM1–JRM2)](ABC_JOINT_REAL_RESIDUE_MOMENTS.md)
· [Splice gcd and the axis specialization (HPS2–HPS6)](ABC_HYBRID_CELL_PRIVATE_SUPPORT_AND_SPLICE_GCD.md)

### Exceptional powerful excess

Write \(\operatorname{rad}_{P^c}\) for the radical
outside the prescribed pool. Across all cells the
number of ordered distinct pairs satisfying
\[
                  D_{ij}/\operatorname{rad}_{P^c}(D_{ij})
                               \ge n^t
\]
is at most \(n^{2d-s-\Delta+O(\rho)+\varepsilon'}\),
where \(\varepsilon'>0\) can be arbitrarily small
at fixed \(L\). This is the mechanism of HPS3 with
a complete-degree allowance. For clarity, if a
divisor \(E\mid D_{ij}\) has depth at least \(z\)
and outside radical at most \(n^a\), enumerate it
at cost \(n^aX^\varepsilon\). For a fixed choice of
\(E\), the compared values are both units modulo
\(ME\), with identical residues there and the
same log bin. The collision bound JRM2 on distinct
compared values gives at most
\[
                  n^{2d-\min(d,s+z)+O(\rho)}X^\varepsilon
\]
ordered value pairs. Indeed, apply it to the uniform
distribution on all eligible distinct compared
values and multiply by the square of their number.
The original row pairs have the additional factor
\(n^{2u+2\rho}\). The count includes equal compared
values, which causes no problem. Uniform fixed-support
enumeration applies to the divisor since \(E\le X\)
and the given pool has polynomial height.

To see the claimed gap, first suppose the mass of
repeated outside primes of \(D_{ij}\) is at most
\(\gamma/2\). Use as \(E\) the part of \(D_{ij}\)
on the repeated outside primes and on \(P\). Its
depth exceeds its outside radical mass by at least
\(t\), so a fixed fine mesh in the outside mass
gives a saving at least \(\min(t,\gamma/2)\) before
the degree and tolerance charges. Otherwise use a
product of squares of repeated outside primes of
mass between \(\gamma/4\) and \(\gamma/3\), available
since their individual costs are at most \(\rho\).
Its depth is twice its outside mass, which gives a
gap at least \(\gamma/4\). A fixed mesh can have
width \(O(\rho)\). The fine-box gain used here, and
in the final count below, follows directly from
\[
 G_v(x;\mathbf a)\ge \min\biggl(\sum_j a_j,\
                        v(x-\max_j a_j)_+\biggr),
                      \qquad 1\le v\le2 .
\]
For a positive right side \(q\), use
\(h_j=\lceil q/(v a_j)\rceil\) on the positive slots.
Then \(h_j a_j\le q/v+a_j\le x\) and
\(g_{j,x}(1/h_j)=a_j\). If the reciprocal weights
sum to at most \(v\), they give the full gain. If they
sum to more, rescale them to sum \(v\) and use concavity
and \(a_j\ge(q/v)h_j^{-1}\).

Discard full cells in which the exceptional pair
fraction is at least \(n^{-\Delta/2}\). The count
just proved and the lower bound on \(N\) give at
most \(n^{s-\Delta/2+O(\rho)+\varepsilon'}\) such
cells. In any remaining full cell the expected
\(\log_n(D_{ij}/\operatorname{rad}_{P^c}D_{ij})\)
is at most \(t+o_X(1)\). This uses the global upper
bound \(\log_n D_{ij}\le 1/L\) at each fixed budget.

### A divisor shared by many rows of a high-mean cell

Fix a remaining full cell of mean at least \(\eta\).
Take a positive integer \(R\) and
\(m=\lfloor\binom{R+2}{2}/4\rfloor\) independent
uniform anchor rows \(j_1,\ldots,j_m\). There is a
nonzero polynomial \(F(Z,W)\) of total degree at
most \(R\) vanishing at the anchors, with coefficient
vector bounded in absolute value by \(X^R\).
For example, evaluate all polynomials with integer
coefficients between \(0\) and \(\lfloor X^R\rfloor\)
on the anchors. Each evaluation lies in an interval
of integers of length \(O_R(X^{2R})\), so two
different polynomials have the same evaluation
vector for all sufficiently large \(X\). Divide
their difference by its largest power of \(W\).
This preserves vanishing at the anchors and makes
\(F(Z,0)\) nonzero.

For every query row \(i\), \(D_{ij_\ell}\) divides
\(F(Z_i,0)\) by reduction at the anchor arguments,
including the diagonal convention. Except at
most \(R n^{u+\rho}\) query rows, that value is
nonzero and of absolute value at most \(X^{3R}\).
Put \(D^{(\ell)}=\operatorname{lcm}_{h\le\ell}D_{ij_h}\)
and \(D^{(0)}=1\). For any anchors,
\[
 \mathbb E_i\log_n D^{(m)}
            \le 3R/L+Rmn^{u+\rho}/(LN).
\]
Choose \(R\) large enough at the fixed budget that
\(3R/(Lm)<t\). The second term divided by \(m\)
tends to zero since \(u<d-s\) by the stated
strict gap hypotheses and \(\rho\) can be
arbitrarily small. The identity
\[
   \sum_{\ell=1}^m\log_n D_{ij_\ell}
     =\log_n D^{(m)}
          +\sum_{\ell=1}^m
              \log_n\gcd(D_{ij_\ell},D^{(\ell-1)})
\]
shows that the mean of the overlap logs, averaged
over the query, the anchors and a uniform position,
is at least \(\eta-t-o_X(1)\). Passing to the
squarefree part outside \(P\) costs at most another
\(t+o_X(1)\) in expectation, by the exceptional
pair estimate. In detail
\(\operatorname{rad}_{P^c}\gcd(D_{ij_\ell},D^{(\ell-1)})\)
divides the overlap and the quotient of the overlap by
it divides \(D_{ij_\ell}/\operatorname{rad}_{P^c}(D_{ij_\ell})\).

Fix the query, a position, and the previous anchors
so that the expectation over the remaining uniform
anchor \(j\) of this squarefree overlap mass is at
least \(\eta-2t-o_X(1)\). This is possible by
averaging, since \(j\) is independent of these fixed
choices. Every such mass is at most \(1/L\). Thus a
positive fraction, with constants allowed to depend
on \(L,\eta\), of the rows \(j\) have outside overlap
mass at least \(\eta/2\). The overlap is supported
on the fixed union of the previous marked values'
primes. For each selected \(j\), choose greedily a
squarefree divisor of this overlap of mass in
\([b,b+\rho]\). There are height-subpower many
possibilities inside the cell, since the support
is on a product of at most \(m\) integers bounded by
\(X\). It follows that the cell has one fixed
divisor \(E\in[n^b,n^{b+\rho}]\) and
\(N X^{-o_X(1)}\) supporting rows \(j\), all with
\(E\mid D_{ij}\), hence
\[
                            Z_j\equiv Z_i\pmod{ME}.
\]
All these rows also share the original logarithmic
bin. The chosen query and the previous anchor pool
may differ among cells. Only the divisor \(E\) is
catalogued across cells.

### The cost of the cross-cell catalogue

Assign one such \(E\) to every remaining high-mean
full cell. There are at most \(n^{b+\rho}\) numerical
choices of \(E\). Fix one and let \(K_E\) be its
number of cells. Each cell contributes one tested
unit-residue and log-bin pair at modulus \(ME\).
Its selected rows give at least
\(n^{d-s-u-O(\rho)}X^{-o_X(1)}\) distinct compared
values. Distinct base cells have disjoint compared
values because each value determines its base
cell. Thus JRM1 at the new modulus, with the same
logarithmic mesh, gives
\[
 K_E n^{d-s-u-O(\rho)}X^{-o_X(1)}
        \ll K_E^\vartheta n^{d/2+O(\rho)}X^\varepsilon,
      \qquad \vartheta>\frac{d/2}{s+b},
\]
using \(\vartheta<1\). In this range \(s+b<d\), so
\(\vartheta>1/2\). The strict inequality ensures
that the fine-box gain reaches \(d-O(\rho)\)
uniformly in \(E\). The arbitrary list allows
separate residues in different cells and general
moduli \(ME\).

Letting \(\vartheta\) approach \(d/[2(s+b)]\) from
above yields the limiting upper exponent
\[
 b+\frac{(s+b)(s-d/2+u)}{s+b-d/2}
       =s-\frac{b(d-s-b)-u(s+b)}{s+b-d/2}.
\]
For the asserted finite estimate first choose
\(\vartheta\) close enough to this threshold,
then choose \(\rho\) and the counting errors small
enough to preserve \(\kappa_0\). The high-mean
cells discarded earlier also satisfy the asserted
bound. All chart choices, the interpolation
degree, and the counting losses are fixed
before taking the height limit.

## A logarithmic tail in each typical full cell

**RSC2.** Fix a tail depth \(\zeta>0\). In the
setting of RSC1 use \(\eta=\zeta\) in the displayed
parameter conditions. For every
\(0<\kappa_0<\min(\kappa,\Delta/2)\) there is a
fixed \(a>0\) such that, outside
\(O(n^{s-\kappa_0})\) full cells,
\[
          \Pr_{i,j}\bigl(D_{ij}\ge n^\zeta\bigr)
                              \le(\log X)^{-a}.
\]
The parameters and chart tolerance are chosen
before the height limit. Dependence of the
implicit constants on the fixed budget is allowed.

Use the same exceptional-cell exclusion as in
RSC1. Put \(T_{ij}=\operatorname{rad}_{P^c}(D_{ij})\)
if \(D_{ij}\ge n^\zeta\) and the pair is good for
the powerful-excess test at \(t\). Otherwise put
\(T_{ij}=1\). In a remaining full cell with
original tail fraction exceeding \((\log X)^{-a}\),
the kept fraction \(q\) is at least
\(\tfrac12(\log X)^{-a}\) for sufficiently large
\(X\), and \(T_{ij}\ge n^{\zeta-t}\) on kept pairs.
For the anchor construction use
\(R=\lceil(\log X)^{2a}\rceil\) and
\(m=\lfloor\binom{R+2}{2}/4\rfloor\).
The coefficient pigeonhole remains valid with
this growing \(R\). The logarithm of the range
of each evaluation is at most
\(2R\log X+O(\log R)\), while the number of
coefficient vectors has logarithm at least
\(\binom{R+2}{2}R\log X\). Also
\(T_{ij}\mid F(Z_i,0)\). The averaged lcm charge
per anchor is
\[
      O(1/(LR))+R n^{u+\rho}/(LN)=o_X(q),
\]
by the same axis calculation. The overlap logs
therefore have average at least
\((\zeta-t)q-o_X(q)\). They vanish unless the
current anchor pair is kept and are at most
\(1/L\). Consequently a proportion \(\gg_{L,\zeta}q\)
of the query, position and anchor choices have
overlap mass at least \(\zeta/4\). As before,
fix the query, position and previous anchors
so that this proportion holds over the new
anchor rows alone. Every selected overlap now
is squarefree, and again has a divisor of
mass in \([b,b+\rho]\).

The union of the previous anchor supports has
at most \(O(m\log X/\log\log X)\) primes. An
integer at most \(n^{b+\rho}\) has at most
\((b+2\rho)\log n/\log\log X\) distinct prime
divisors for all sufficiently large \(X\).
Both bounds follow from the factorial lower
bound for the product of distinct primes.
Counting subsets of this size in the anchor
pool gives at most
\[
                   n^{\,4a(b+2\rho)+o_X(1)}
\]
possible squarefree divisors in the cell.
Indeed if the pool size \(H\) is at least twice the
size cap \(J\), use
\(\sum_{j\le J}\binom Hj\le(J+1)(eH/J)^J\).
Otherwise \(2^H=n^{o_X(1)}\) suffices.
The retained number of rows for a fixed \(E\)
is thus reduced by at most an exponent
\(4a(b+2\rho)+o_X(1)\) from RSC1.
The factor \(q\) itself is height-subpower.
In the fixed-\(E\) list calculation add this
loss to \(u\). Choose \(a>0\) sufficiently small,
and then the chart and counting allowances
sufficiently small, to preserve the requested
\(\kappa_0\). This proves the per-cell tail bound.

## The direct bound at combined depth below \(d/2\)

For \(s<d/2\), fix a second row \(j\) and enumerate
the divisors of \(W_j\) at height-subpower cost.
If \(D_{ij}\ge n^q\), the first compared value
lies in one cell at modulus \(M D_{ij}\) and the
original log width, of combined depth at least
\(s+q-O(\rho)\). The \(K=1,\vartheta=1\) bound of
JRM1 gives at most
\(n^{d+u-\min(d/2,s+q)+O(\rho)}X^\varepsilon\)
first rows, including the degree charge. There
are at most \(n^{d+u+O(\rho)}\) second rows since
the compared box has at most \(n^{d+O(\rho)}\)
root tuples. The resulting pair exponent is
\[
     2d-s-\min(d/2-s,q)+2u+O(\rho)+\varepsilon'.
\]
Take \(0<q<\eta\). Since \(\log_n D_{ij}\le1/L\),
a full cell of mean at least \(\eta\) has a
positive fraction of pairs in the tail. The
number of high-mean full cells has a power
loss as soon as \(\min(d/2-s,q)>2u\) by a fixed
margin. This is the same single-anchor mechanism
as HPS2 and uses the exact log cells here.

## Consequence for saturation charts

At zero multiplicity in the limiting sense,
every fixed positive mean excess above the
imposed modulus forces a power-small count of
full cells, at every fixed mixed depth \(s<d\).
The full-cell bounds also give power-small row
and within-cell pair mass for those high-mean
cells, compared respectively with the balanced
reference exponents \(d\) and \(2d-s\). Indeed
multiply the cell count by at most
\(n^{d-s+\rho}\) or \(n^{2(d-s+\rho)}\).
The usual alphabet bound and JRM2 show that
full cells carry almost all rows on a regular
saturating chart, after taking its original losses
and fine-layer errors smaller than the cell
tolerance. The small-cell bound is the alphabet
count times the size cap. The large-cell bound
uses the sum of squared occupancies and divides
by the size threshold. For an individual full
cell under the mean threshold, Markov's inequality
bounds the proportion of pairs of excess depth
at least \(\zeta>0\) by \(\eta/\zeta\).
In an RSC2-typical cell, a random ordered packet
of at most \((\log X)^{a/3}\) rows has no pair of
extra depth at least \(\zeta\) with probability
\(1-O((\log X)^{-a/3})\), by the union bound on
ordered pairs.

The typical-cell tail bound leaves open a power
saving on pairs of squarefree excess within each
typical cell. Improving the unrestricted \(abc\)
counting exponent still requires a way to control
the number of original rows, for example through
the third-coordinate cap on an auxiliary family
of genuine splices.

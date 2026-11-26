# Same-role overlap and fresh radicals in hybrid cells

A pair of complete values sharing a positive piece of radical
in a hybrid cell gives a small divisor of many values in that
cell. The divisor is chosen within one fixed anchor support.
Across cells it is paid once, and a mixed list count on the
quotients gives a power saving. This extends the privacy of
complete values beyond the half-depth threshold. Together
with the two-cheapest-letter inequality, it puts physical
shifts and actual splice outputs beyond every fixed radical
cap on typical pairs.

## Power rarity of same-role overlap

Put \(n=X^L\) and \(d=1/3\), taking the height limit at
fixed \(L\). Work in a chart of true positive primitive
rows of height at most \(X\). Let \(Z\) be one original
complete role, with radical cap \(n^{d+O(\rho)}\),
exact valuation roots of total cap at most \(d+O(\rho)\),
and maximum exact-layer cap at most \(\rho\). Here
\(\rho>0\) is a chart tolerance. A fixed positive
coefficient is allowed. Fix a prescribed pool \(P_0\)
supported on a polynomial-height integer and containing
the primes of that coefficient. Every individual prime
on \(Z\) outside \(P_0\) has \(\log_n p\le\rho\).
The last property follows from the fine layers when
they represent all these primes. The maximum complete
\(Z\)-degree in the chart is at most \(n^{u+\rho}\),
where \(u\ge0\).

Let \(M\) be a prescribed polynomial-height modulus
at which \(Z\) and its representation are units. Place
each row in the cell determined by \(Z\bmod M\) and its
logarithmic bin of width \(\asymp n^{-\xi}\), for
\(\xi\ge0\). Take \(s\in[0,d)\) fixed, with
\(\log_n M+\xi=s+O(\rho)\). Call a cell full here if
it contains \(N\) rows with
\[
                     n^{d-s-\rho}\le N\le n^{d-s+\rho}.
\]
For ordered rows \(i,j\) in the same cell define
\[
       I_{ij}=\log_n
             \operatorname{rad}_{P_0^c}\gcd(Z_i,Z_j).
\]

**HPO1.** Fix an overlap threshold \(f>0\). If
\(d/2\le s<d\), choose
\[
             0<b<\min(f,d-s),\qquad q>0,\qquad
     K_*=\frac{b(d-b)/2-s(u+q)}{s-(d-b)/2}>0 .
\]
For every fixed \(0<K_0<K_*\), the number of full
cells in which at least an \(n^{-q}\) fraction of
ordered row pairs have \(I_{ij}\ge f\) is
\(O(n^{s-K_0})\). The chart tolerance and counting
allowances are sufficiently small in terms of the
fixed parameters. Consequently among all full
cells the pairs with \(I_{ij}\ge f\) have exponent
at most \(2d-s-\min(q,K_0)+O(\rho)\).
For \(0\le s<d/2\) the corresponding pair saving,
before allowances, is at least
\[
                    \min(f,(d+f)/2-s)-u
\]
whenever this is positive. At zero limiting complete
multiplicity the overlap tail therefore has a fixed
power loss at every fixed \(s<d\).
The constants are uniform in admissible \(M,P_0\) and
in the bin translates, with the usual dependence on
the budget and fixed root data.
[Mixed residue and real lists (JRM1–JRM4)](ABC_JOINT_REAL_RESIDUE_MOMENTS.md)
· [Local prescribed-pool count (HPS1)](ABC_HYBRID_CELL_PRIVATE_SUPPORT_AND_SPLICE_GCD.md#a-fixed-prime-pool-inside-a-hybrid-cell)

For the first assertion, in any full cell with the
indicated tail, fix an anchor \(i\) having at least
\(Nn^{-q}\) qualifying second rows. For each of those
rows select a squarefree divisor \(E\) of
\(\operatorname{rad}_{P_0^c}\gcd(Z_i,Z_j)\) with
\(n^b\le E\le n^{b+\rho}\). Greedy selection is
possible by the individual-prime cap. There are only
\(X^{o_X(1)}\) possibilities inside this cell, since
all divide the fixed integer \(Z_i\). Some fixed
\(E\) therefore divides \(Z_j\) on at least
\(Nn^{-q}X^{-o_X(1)}\) rows. Assign one such \(E\)
to the cell. Its choices over the entire chart
number at most \(n^{b+\rho}\).

Fix a catalogue entry \(E\) assigned to \(K_E\) cells.
They yield at least
\[
              K_E n^{d-s-u-q-O(\rho)}X^{-o_X(1)}
\]
distinct compared values. Distinct initial cells
have disjoint compared values, and the \(n^{u+\rho}\)
degree loss is paid once. Strip the full
\(E\)-supported part of each selected value and
fix these supported parts slot by slot. This
costs height-subpower for a fixed \(E\), uniformly
in the choice of \(E\); the same prime pool may be
split into the fixed number of slots at subpower
cost. Strip and fix the support before taking
dyadic bins of the residual root ranges. Their
total cap is at most \(d-b+O(\rho)\), with the
original upper bound on each slot up to binning
loss. Indeed \(E\) is outside the fixed coefficient
support, and its primes are removed from the
disjoint variable exact layers. Each fixed
supported-part state sends the original \(K_E\)
residue and log-bin conditions to at most
\(O(K_E)\) such conditions at the same depth.
Division by the known factor permutes the residues
and translates each logarithmic interval.

Choose \((d-b)/(2s)<\vartheta<1\), which is
possible for \(s\ge d/2\). The lower endpoint
exceeds \(1/2\) since \(b<d-s\). The fine-cap gain
in JRM4 gives the upper bound
\[
 K_E n^{d-s-u-q-O(\rho)}X^{-o_X(1)}
             \ll K_E^\vartheta n^{(d-b)/2+O(\rho)}
                  X^\varepsilon
\]
after summing the supported-part states. The
bound uses the actual residual cap as an upper
cap: the quantity \(r-\min(r,2\vartheta s)/2\)
is nondecreasing in \(r\). Let
\(\vartheta\) approach \((d-b)/(2s)\) from above.
The limiting exponent for the number of catalogue
entries weighted by their assigned cells is
\[
       b+\frac{s\{s-(d+b)/2+u+q\}}{s-(d-b)/2}
                               =s-K_* .
\]
For a finite assertion choose \(\vartheta\) close
enough, then absorb all tolerances and inner
counting allowances into \(K_*-K_0\).
In cells not discarded, multiply \(n^{-q}\)
by the full pair population. In discarded cells
use the upper cell-size bound. The total number of
cells is at most \(n^{s+O(\rho)}X^{o_X(1)}\).

For \(s<d/2\), fix an anchor and use its prime
support in (HPS1). The surviving neighbors after
requiring overlap mass at least \(f\) have exponent
at most \(u+\max((d-f)/2,d-f-s)\), with the usual
small allowances. Compare with the lower full-cell
neighbor exponent \(d-s\). An overlap threshold
above the total raw cap gives an empty event.

## Primitive comparison and physical shift outputs

**HPO2.** In the same fine hybrid chart, assume all
three original radical caps are at most
\(n^{d+O(\rho)}\), and the compared role has a fixed
compact positive scaled-degree band \([a_0,a_1]\)
of mass at least \(\gamma>0\) outside \(P_0\) on
each row. Let \(P\) be any prescribed pool of
polynomial height, which can differ from \(P_0\).
For every fixed \(B>0\), the pairs in full cells
with nonzero \(\Delta=Z_i-Z_j\) and
\(\operatorname{rad}_{P^c}|\Delta|\le n^B\) satisfy
\[
 \#\{\text{such pairs}\}
 \ll n^{2d-s-\delta+O(\rho)}
       +n^{d+2u+O(\rho)+\eta}X^\varepsilon .
\]
Here \(\delta>0\) is any fixed overlap loss
available from HPO1 at threshold \(f=\gamma/2\),
with counting allowances taken strictly inside
that loss, and \(\eta>0\) is arbitrary. Thus this
class is power-smaller than \(n^{2d-s}\) if HPO1
has a strict margin and \(2u<d-s\) by a fixed
margin. The height of \(P\) rather than its
normalized radical bounds the pool enumeration.

Indeed discard the high-overlap pairs of HPO1. With
\(g=\gcd(Z_i,Z_j)\), at least \(\gamma/2\) of the
band on \(Z_i\) survives division by \(g\), outside
the anchor's prime support and outside \(P_0\).
The primitive additive comparison with entries
\(Z_i/g,Z_j/g,|\Delta|/g\) therefore has maximum
at least \(X^{a_0\gamma/2}\). Both endpoint caps
are \(d+O(\rho)\) and the outside-\(P\) shift cap
is \(B\). Fix the powers of primes in \(P\) on
all three primitive entries at height-subpower
cost. They enter as arbitrary fixed coefficients.
Use the variable anatomic factorization outside
\(P\) with total variable raw cap at most
\(2d+B+O(\rho)\). Anatomic initial factors can be
fixed within the allotted counting tolerance.
The complete height floor meets the \(n^C\)
condition of (MCB2) for sufficiently small \(L\)
chosen after \(B,a_0,\gamma\) and the allowances.
The stated MCB2 height range permits a common
base \(n=X^L\) throughout. If the outputs are
split into dyadic shells, the number of shells
contributes only a height-subpower factor.
[Bernert–Browning–Lichtman–Teräväinen, Lemma 2.2](https://arxiv.org/html/2410.12234v2#S2)

To account exactly for restoration, let \(b_0\)
be the normalized mass of the outside primes of
\(g\) missing from all three primitive entries.
Every such prime was on both original endpoints.
In a fixed narrow \(b_0\)-bin, use the two endpoint
caps \(d-b_0+O(\rho)\) on the primitive comparison.
The two-cheapest-role bound costs \(d-b_0\) in
the exponent, up to allowances. The missing prime
record costs \(b_0\). Once this record, the primitive
triple and \(P\) are known, all possible \(g\)'s
are supported on a prescribed polynomial-height
pool, so their powers cost height-subpower. This
restores both endpoints exactly. The choice of sign
costs a bounded factor, and the two original
complete-degree costs are at most \(2u+2\rho\).
This proves the displayed estimate.
[Two-role prices for primitive comparisons (ESL2a)](ABC_SPLICE_EXHAUSTED_SUPPORT_LIFT.md#two-letter-costs-after-restoring-the-common-factor)
· [Two cheapest roles with arbitrary fixed coefficients (MCB2)](ABC_MULTIPLE_SHORT_COFACTOR_BOUNDARY.md#bounded-group-dictionaries-with-arbitrary-fixed-coefficients)

## Direct splices and the matched population

**HPO3.** Give the original roles disjoint prime
supports by a fixed role colouring. Mark \(W\) and
take the two actual splice operands from the other
roles \(U,V\) on distinct original rows of one
hybrid cell, with the original signs. Keep
positive splices of maximum complete height
\(\ge X^{h_0}\) for any fixed \(h_0>0\).
Suppose the original operand degrees are at most
\(n^{u_U+\rho},n^{u_V+\rho}\) and the two operands
have outside-\(P\) caps \(r_U',r_V'\) in base \(n\).
For every fixed \(B>0\), the number of such pairs
whose splice output has outside-\(P\) radical at
most \(n^B\) is
\[
              O\!\left(n^{(r_U'+r_V')/2+u_U+u_V
                                  +O(\rho)+\eta}X^\varepsilon\right).
\]
Here \(P\) is any prescribed polynomial-height pool
and \(L\) is sufficiently small in terms of the
fixed cap and height data and the allowances.
No same-role privacy hypothesis is needed for
this output estimate. The operands are coprime
by the role colouring. Fix all \(P\)-supported powers
at height-subpower cost and apply (MCB2) with the
two operands as the charged pair of roles. After
fixing the two complete operands, restore their
original rows at the indicated degree costs.

In a matching of complete labels at exponent \(d\),
the full cells carry essentially all rows after the
initial chart and fine-cap allowances are chosen
smaller than the full-cell tolerance. For details,
the alphabet bound controls the rows in underfull
cells. JRM2 bounds the sum of the squares of cell
sizes by \(n^{2d-s+o(1)}\); division by the large-cell
size threshold controls the rows in overfull cells.
The retained ordered-pair population has exponent
\(2d-s\) up to these allowances, by Cauchy–Schwarz
and the alphabet bound. The repeated complete labels
and diagonal can be omitted in a matching. The
same lower pair exponent holds for any retained
family of full cells carrying \(n^{d-o(1)}\) rows
with degree loss \(u<d-s\) by a fixed margin,
independently of whether this family is obtained
using JRM2.

In the original shell \(X/2<C\le X\) at least one
orientation of each unordered pair gives a positive
splice with maximum complete height comparable to
\(X\). For a sum splice the two positive outputs
sum to \(C_i+C_j\). For a difference splice the
two signed outputs sum to the two genuine positive
marked values, and in the positive orientation
an original target is one of the two operands.
Consequently the same lower exponent applies to
these actual positive, tall splices. At balanced
caps and zero multiplicity, HPO2 and HPO3 therefore
give a fixed power loss on every fixed outside
radical cap \(B\) for shifts and splices, at every
fixed \(s<d\). With positive degree exponents the
comparison requires the strict conditions in HPO2
and, for HPO3,
\[
             (r_U'+r_V')/2+u_U+u_V<2d-s .
\]
These are comparisons against the retained
pair populations on the original three-cap chart.
A prescribed finite-price family of pools can
also be tested by summing their price into the
relevant counting exponent. The physical-shift
result does not give the same estimate for pools
chosen at no price separately for every anchor
at depths \(s\ge d/2\).

The budget threshold can depend on \(B\). The
conclusion says that the normalized outside radical
of these derived outputs exceeds every fixed value
in the height-first limit. It leaves open a bound
on the number of original occupied cells or
original triples. The shallow pair comparison
also leaves the large-multiplicity regime outside
its margins.

## The same-role shift test after recording large primes

**HPO4.** Let \(Z\) be a role with positive nonatomic
scaled valuation mass \(c_Z\) in a balanced limiting
profile. At zero limiting complete multiplicity, for
each fixed \(0<s<c_Z\), one can choose a full-factor
record on \(Z\) and an arbitrarily cheap full lock
on another available role so that the following
conclusion holds on the retained exponent-full
residual chart. Put \(r=d-f\), where \(f\) is the
captured raw mass on \(Z\), including any tail charges.
For actual fixed mixed probes of combined depth \(s\)
at which \(Z\) is a unit, the same-cell pair population
has exponent \(2r-s\) in the height-first convention.
For every prescribed polynomial-height pool \(P\) and
fixed \(B>0\), its distinct-value pairs with physical
shift radical outside \(P\) at most \(n^B\) are power-small.
A fixed gap \(s<c_Z\) and the choice of \(B\) precede
the selection of the small-budget and chart allowances.
The power-saving comparison also works for a
small positive degree when HPO1, HPO2 and JRM6
with \(r\) in place of \(d\) have fixed margins.
For fixed positive degree this uses a lower pair
exponent \(2r-s\) up to allowances and allows an
upper exponent larger by the degree allowances.

To see this, apply JRM5 and choose a mass bin and a
full factor \(F\) removing all varying primes of
normalized size above \(a\), plus an arbitrarily
small valuation tail. Choose the upper bound \(a\)
small enough for the pair estimates, after choosing
the fixed depth and count margins. The original
three caps are kept. The chosen row cohort has
exponent at least \(r-o(1)\), and the tested
quotients \(Z/F\) have raw cap \(r+o(1)\) with
small exact subroot caps. The nonatomic mass of
the limiting measure remains on the quotients,
up to the arbitrarily small tail charge. In
particular there is a compact positive-order
band carrying a fixed positive mass on the
residual part. One may choose its ends to be
continuity points of the original measure. The
finitely many large-prime records remove only an
atomic submeasure in the weak limit, so this band
can be chosen before taking the cap \(a\) small.
The chart alphabet and the rank bounds in JRM4
give the matching lower and upper pair exponents,
and JRM6 makes all but a power-small fraction
of pairs lie in cells at the requisite size.

Use HPO1 on \(Z/F\), with \(r\) replacing \(d\)
throughout its proof. Its count uses only the
total cap of the tested role, the complete degree,
and the full cell sizes. It does not use the other
roles' radical caps to lower that total. The compact
band then ensures a positive primitive height after
dividing the two same-role values by their gcd, as
in HPO2. For the count with radical bound \(B\),
enlarge \(P\) by the prime support of \(F\).
The new pool has polynomial height and does not
invalidate the outside-cap upper bound. Both
complete endpoints have cap \(r+o(1)\) outside this
pool. The same two-role count and support restoration
in HPO2 now costs \(n^{r+2u+o(1)}\) before comparison
with \(n^{2r-s-o(1)}\). Since \(r\ge c_Z-o(1)>s\)
by a fixed margin, this is a power loss at zero
limiting complete multiplicity. In an arbitrary
balanced profile at least two original roles have
nonatomic mass at least \(1/6\), so a chart with
this conclusion is available at any fixed \(s<1/6\).
The price of the first full-factor record is allowed
to be positive on the original radical scale.
[Subroots and pair-mass pruning (JRM5–JRM6)](ABC_JOINT_REAL_RESIDUE_MOMENTS.md#fine-residual-caps-from-small-individual-primes)

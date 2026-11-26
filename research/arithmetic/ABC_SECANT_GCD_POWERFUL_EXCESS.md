# Uniform secant savings for powerful parts of a common gcd

For two genuine small-radical sums, a common secant gcd
specifies congruences on all three original terms. At every
prime of the gcd, at least two of those terms are units.
This gives a uniform power-small exceptional count for
supported or repeated prime powers in the gcd near zero
complete multiplicity, even if the associated prime pool
is prescribed separately for each anchor. The proof uses
the one-class moment at shallow depth and an inexpensive
divisor record for the excessive power.

## A one-class count at low combined depth

Fix \(L>0\) and put \(n=X^L\). Work with an exact valuation
box on one complete integer with positive raw cap \(r\).
The orders, finite coordinate list and cap box are fixed
before the height limit at this budget. Let \(M\) be any
polynomial-height modulus on which the candidate values
are units. For one prescribed residue modulo \(M\) and
one prescribed logarithmic window of width \(O(n^{-\xi})\),
the count is at most
\[
                    O(n^{r-s}X^\varepsilon)
              \quad\text{if}\quad
           0<s\le r/4,\qquad \log_n M+\xi\ge s .
\]
The same applies to a residue alone, by taking \(\xi=0\)
and a subpower list of bins covering polynomial height.
If a layer cap \(a_i\ge s\), the elementary auxiliary
arm of JRM1 for this slot gives the claim, by fixing
the other roots and counting choices in its torsion
progressions and its real interval. Otherwise
\(\sum_i1/\lceil s/a_i\rceil\ge r/(2s)\ge2\).
At each reciprocal weight up to
\(1/\lceil s/a_i\rceil\), the moment gain has slope
\(s\). Allocate weight of sum two in these linear
segments. JRM1 with \(\vartheta=1\) gives gain
\(G_2(s)\ge 2s\). The gain is nondecreasing with
the combined depth, which proves the result.
[Exact-layer joint list moments (JRM1)](ABC_JOINT_REAL_RESIDUE_MOMENTS.md#an-arbitrary-list-of-congruence-and-logarithmic-bins)
· [One-class shallow moment and adaptive pools](ABC_SECANT_GCD_MEAN_OVERLAPS.md#a-one-residue-count-at-low-combined-depth)

## Power-small powerful excess on actual regular charts

Write \(d=1/3\) and let \(S\) be a regular balanced
chart of original positive primitive triples \(A+B=C\) at
height comparable to \(X\), in the height-first convention.
Require each original letter to have raw radical cap
\(d+\tau\) in base \(n\), and each complete value to have
degree at most \(n^{u+\tau}\). Work in the common exact
valuation boxes with raw cap \(d+O(\tau)\), with anatomic
truncations and their losses paid as usual. The assertion
below also applies to a subset of such a chart. No
additional cap is placed on the differences.

Fix \(0<\eta<d/20\), assume \(u<\eta/30\), and choose
all chart and inner count allowances sufficiently small
in terms of \(\eta\). At every sufficiently small fixed
budget, for each anchor \(x=(A,B)\) let \(P_x\) be a
prescribed set of primes supported on an integer at most
\(X^K\) with fixed \(K\). The pool may depend on the
anchor. For a distinct row \(y=(A',B')\) set
\[
 g=g(x,y)=\gcd(A-A',B-B')>0,\qquad
 E_x(g)=
   \prod_{p\in P_x}p^{v_p(g)}
       \prod_{\substack{p\notin P_x\\p\le n^{\eta/12}\\v_p(g)\ge2}}
                                  p^{v_p(g)-1}.
\]
The fixed pool is charged at zero radical price.
The second product measures the additional powers
beyond the radical on small outside primes.

**SGP1.** Uniformly in the anchor and its prescribed
pool, the neighbours with \(E_x(g)\ge n^\eta\)
number
\[
                  O_{K,L,\eta,\varepsilon}
                           (n^{d-\eta/12}X^\varepsilon).
\]
In a near-full chart their ordered pair population
therefore has a uniform positive normalized power
loss. The result allows any limiting valuation
profiles of the original letters. In particular it
requires neither thin exact layers nor a high-degree
floor on the residual factors.

We first record the cost and depth of a suitable
divisor \(Q\mid g\). Let \(b\) be the radical mass of
the primes in the second product in \(E_x(g)\).
If \(b\le\eta\), choose \(Q\) to contain the whole
\(P_x\)-supported part of \(g\) and the entire
prime powers on those outside primes. Then
\(\log_n Q\ge b+\eta\). The outside radical price is \(b\).
If \(b>\eta\), select entire outside prime powers
from that set until the raw mass first reaches
\(\eta/2\). The prime cap gives mass \(b'\) in
\([\eta/2,7\eta/12]\), and \(\log_n Q\ge2b'\) because
each selected valuation is at least two. In this
case drop the pool part. For a fixed upper raw bin
endpoint \(b_+\) of either price, there are at most
\(n^{b_++o_X(1)}\) choices of \(Q\): enumerate the
outside radical, then the powers supported on its
union with the prescribed pool. This is uniform in
\(P_x\) at polynomial height. Use a fixed small
mesh for the raw cost, so replacing the price by
\(b_+\) spends an arbitrarily small prescribed
fraction of the margin.
[Uniform enumeration on fixed support (EB5)](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)

Fix a recorded \(Q\). At each of its prime powers
\(p^j\), at least two of the three anchor complete
values are units modulo \(p^j\) by primitivity.
Since the differences in all three roles are
divisible by \(Q\), the same roles are units on the
query. For each \(T=A,B,C\), let \(Q_T\) be the
product of the powers of \(Q\) at which \(T(x)\)
is a unit. Then
\[
                 \sum_T\log_n Q_T\ge2\log_n Q.
\]
Choose \(T\) with largest unit depth. The
query \(T\) lies in one prescribed unit class
modulo \(Q_T\). If the chosen record has
depth at least \(b+\eta\), one can use a saving
\(s=2(b+\eta)/3\) in the one-class bound. This
is within \(r_T/4\) at late balanced boxes
by \(\eta<d/20\).
If the record has depth at least \(2b'\),
use \(s=4b'/3\), also within this range.
These are lower bounds on the depth, so larger
actual \(Q_T\) cause no change in the argument.
Multiply the complete dictionary count by its
original degree \(n^{u+\tau}\). After summing the
record dictionary, the two margins left by
depth minus raw price are respectively
\[
              \tfrac{2\eta-b}{3}\ge\tfrac\eta3,
                      \qquad
              \tfrac{b'}3\ge\tfrac\eta6 .
\]
They absorb the factor \(n^u\) and the fixed
cap, bin and counting allowances. This proves
SGP1. The anchor's pool can vary arbitrarily
because the record enumeration and the one-class
bound are uniform within each anchor count.

## Relationship to the full common-cofactor problem

The conclusion gives an absolute power loss at fixed
positive excess threshold for the repeated powers on
small primes and for any polynomial-height pool
prescribed at each anchor. It sharpens the necessary
shape of a putative high common-gcd pair population:
its numerical mass must be found in squarefree
outside primes or in primes above the small cutoff,
except for arbitrarily small allowed excess. The
proof gives no uniform power saving for those
remaining squarefree parts, since recording a
squarefree divisor on unknown primes costs its
whole depth, while the guaranteed unit depth in
one complete role can be only two thirds of it.
The full-gcd average and fixed-budget tails concern
different measures of the remaining pairs.
[Mean affine secant gcd with adaptive overlaps](ABC_SECANT_GCD_MEAN_OVERLAPS.md)
· [Full common-gcd tails at fixed budgets](ABC_SECANT_NORMALIZED_GCD_TAIL.md)

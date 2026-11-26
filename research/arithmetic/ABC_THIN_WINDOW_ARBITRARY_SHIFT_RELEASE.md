# Arbitrary shifts with thin valuation windows

A fixed number of powered roots can be released with integral
cofactors longer than any fixed power of the radical scale.
The proof uses a gap among the orders of the few retained roots.
This also extends the arbitrary-coefficient dictionary bound
to groups with a growing number of original positions in a
fixed number of sufficiently thin windows.

## A gap among finitely many orders

**TW1.** Fix \(M,B,\theta,H>0\). There is \(\xi_0>0\)
such that the following holds at sufficiently small fixed
budgets \(L\), with \(n=X^L\). In a positive primitive
three-term sum with maximum \(Q\) satisfying
\[
                  n^{L^{-\theta}}\le Q\le X^H,
\]
suppose each term is a fixed coefficient times an integral
cofactor of size at most \(n^{B L^{-\xi_0}}\) times a
product of powers of positive integers at most \(n^B\).
There are at most \(M\) such roots total and the positive
orders are fixed before the height limit. The count is
\(O(n^\eta X^\varepsilon)\) for every fixed \(\eta>0\) and
all sufficiently small fixed \(L\), uniformly in the
coefficients. Within-letter overlap is allowed.

Let \(R>2(M+2)\) be fixed. Choose \(b_0>\xi_0\) so small
that \(R^{M+1} b_0<\theta\), and put \(b_j=R^j b_0\).
At least one of the \(M+1\) disjoint half-open intervals
\[
             [L^{-b_j},L^{-b_{j+1}}),\qquad 0\le j\le M,
\]
contains none of the at most \(M\) orders. For that \(j\),
absorb all roots of order \(< L^{-b_j}\) into their integral
cofactors, moving every other root into the high list.
The new coordinates are at most \(H_*=n^{B' L^{-b_j}}\)
for a fixed \(B'\). The high orders are at least
\(E=L^{-b_{j+1}}\). For sufficiently small \(L\),
\(Q>2H_*^3\). The normalized exponent in (MCB1) is at most
\[
       O_{M,B'}(L^{-b_j}E^{-1/(2(M+2))})=o_{L\to0}(1).
\]
If the high list is empty, the rational reconstruction
in (MCB1) counts at most one solution. Each exact sum
determines its three complete values. Any extra original
factorizations have height-subpower multiplicity at fixed
orders. The bound follows.
[Boundary release with integral cofactors (MCB1)](ABC_MULTIPLE_SHORT_COFACTOR_BOUNDARY.md#releasing-the-high-roots-and-all-three-cofactors)

## Fixed shifts with a wider integral correction

**TW2.** In (MCB3), fix the bounds on the number of roots,
their \(n\)-height, the opposite raw cap and the polynomial
complete height. Given \(\eta,\theta>0\), there is \(\xi>0\)
such that the same \(O(n^\eta X^\varepsilon)\) bound holds
with the correction permitted to have height \(n^{L^{-\xi}}\)
and with the complete height floor \(Q\ge n^{L^{-\theta}}\),
where \(Q=\max(U,V,|t|)\). All quantifiers use the height
limit first, and the result is uniform in the nonzero signed
shift and fixed coefficients.

For the proof, refine original roots to dyadic sizes at
height-subpower cost. If an opposite block is fixed with
full product \(F\ge n^{L^{-2\xi}}\), the compressed term
is in a prescribed unit class modulo \(F\). Write
\(w=\log_n F\). Trim only the low exact orders on this
block to obtain a \(k\)-full divisor \(D\) with
\(k=\lfloor w/(2(A+1))\rfloor\) and
\(\log_n D\ge w/2\), where \(A\) bounds the opposite raw
mass. In the (GW) parameters the correction height
exponent is \(O(L^{-\xi}/w)=O(k^{-1/2})\). The retained
root heights are \(D^{O(1/k)}\). For an amplified exact
monomial all coordinates divide it since the factors are
positive integers. Its exact fibres and the ordered
decompositions of amplified coordinates have subpower
size in the height limit. Thus (GW) and fixed-order
amplification count the original compressed values at
arbitrarily small normalized cost, uniformly in \(F\),
because \(\log_n D/k\) stays bounded.

If a box had count exponent at least \(\eta\), every
block in this height range would cost at least \(\eta/2\)
to fix for sufficiently small budgets. Set aside positions
of raw cost at least \(\eta/8\) on the opposite letter.
Greedily pack the other positions into boundedly many
blocks of cost at most \(\eta/4\). Their total full product
has height \(n^{O(L^{-2\xi})}\), since each stays below the
release threshold. The set-aside list, together with the
original compressed list, has bounded length \(M\)
depending only on the input bounds and \(\eta\). Their
remaining two integral cofactors are \(n^{O(L^{-2\xi})}\).
The shift has no variables. Choose \(2\xi<\xi_0\) for (TW1)
with this \(M\) and \(\theta\). Applying (TW1) with an
allowance smaller than \(\eta\) gives the contradiction.
Anatomic coefficient and dyadic losses, where needed, may
be charged within arbitrarily small fixed allowances.
[Amplified full-modulus release at a variable scale](ABC_RADICAL_SCALE_REAL_GAP.md#concentration-and-release-at-a-variable-gap)
· [Fixed integral correction (MCB3)](ABC_MULTIPLE_SHORT_COFACTOR_BOUNDARY.md#compressed-fixed-shifts-without-a-unit-shift-modulus)

For example, after fixing the outside record on one variable
letter, let all its remaining orders \(Le\) belong to a
bounded number of intervals of width \(c L^{1-\xi}\)
for fixed \(c>0\). In each interval put the product of the
roots to the minimum original order and move the positive
order differences into one integral correction. The root
heights are \(n^{O(1)}\) and the correction is at most
\(n^{O(L^{-\xi})}\). Taking a slightly smaller window
exponent if necessary to absorb the fixed constant gives
(TW2). The windows can contain a growing number of exact
orders. A weak limiting atom alone gives no such rate.

## Arbitrary coefficients with thin-window dictionaries

**TW3.** In (MCB2) fix a normalized allowance \(\eta>0\),
raw cap \(A\), an upper polynomial height cap, a lower
height parameter \(\theta>0\), and bounds on group size and
number of windows. For some \(\xi>0\) and all sufficiently
small fixed \(L\), the same two-cheapest bound holds
under \(Q\ge n^{L^{-\theta}}\) if each group consists
either of boundedly many original positions or of positions
whose scaled orders \(Le\) lie in a bounded number of
windows of width \(L^{1-\xi}\). Groups belong to a single
role, and keep their actual joint group costs. Arbitrary
fixed coefficients are permitted on all three roles.

Here are the two changes to the concentration proof.
Work along a violating sequence with limiting count
exponent \(d\ge (w-\max_T w_T)/2+\eta\). All exponent
limits are taken in base \(n\) after the height limit.
Set \(S=L^{-2\xi}\). A selected block of groups is large
if its normalized log-height divided by \(S\) tends to
infinity, passing to subsequences to distinguish this
from boundedness. Under a lock from such a block,
apply the (MCB2) character moment argument.
In a thin-window group of raw sum \(b_G>0\), collapse
\(h_G=\lceil A/b_G\rceil\) samples in the diagonal moment
and then combine the positions in each window at their
minimum exact order. Their new roots have height
\(n^{h_G b_G}\le n^{2A}\); the integral correction has
height at most \(n^{2A L^{-\xi}}\). For bounded groups
retain the original collapsed coordinates.
One resulting exact value on a single role has
divisor-many representations, including both the
collapsed original positions and the short-window
representatives. The arbitrary-order rebased (GW)
argument therefore bounds the collapsed residue class
at vanishing normalized cost. Its wide exponent tends
to zero along the sequence, uniformly over groups,
because the modulus has log-height divided by \(S\)
tending to infinity. Hölder gives half the summed
group caps on the opposite roles, just as in (MCB2).

If a letter is forced into boundedly many original groups
and a leftover product with normalized log-height \(O(S)\),
write the retained groups with boundedly many roots as
above. Their integral corrections join the leftover
product. The cofactor still has height \(n^{O(S)}\)
and the roots have height \(n^{O(1)}\). This description
allows vanishing-cost release of such a letter or a pair
of such letters under any large opposite lock. For a pair
use the signed ratio of the variable parts, and first
colour primes by original role at height-subpower loss.
The ratio and all its amplified products then have no
cross-role cancellation. Each positive coordinate divides
the corresponding numerator or denominator. This gives
the requisite divisor bound on exact fibres. If weighing
two positive complete values, primitivity recovers them
from their exact ratio.

For uniformity in \(\xi\), fix \(0<\kappa<\eta/4\) and
call a letter supported if \(w_T>w-2d+2\kappa\).
Every large block on a supported letter costs at least
\(\kappa\) by the moment bound. Greedy packing away from
groups of cost at least \(\kappa/4\) therefore gives the
bounded-group description with an \(O(S)\) remainder.
If a letter costs \(<d-\kappa\) and another letter is
supported, fix the cheaper letter. Every large lock
on the third letter releases the supported letter, so
the third letter concentrates with the same cost floor.
If the fixed letter's variable part was large, its
full lock releases the two other letters. If it was
bounded relative to \(S\), the two concentrated letters
and the fixed letter are counted by (TW1) at vanishing
normalized cost. Both alternatives contradict the cost
of fixing the cheaper letter.

Consequently \(4d\le w+\min_T w_T+4\kappa\): a violation
would make a cheapest letter cost \(<d-\kappa\), with
one of its other two letters supported. This inequality
also excludes three supported letters, since it gives
\(d\le w/3+\kappa\). Two supported letters would leave
the third at cost \(<4d-w-4\kappa\le \min_T w_T\), also
impossible. A most expensive letter is supported since
the failure margin is \(\eta>\kappa\). It is therefore
the unique supported letter.

For the last switches, write \(U,V\) for the other roles
with \(w_V\le w_U\), so \(w_U+w_V\le 2d-2\eta\).
If \(U\) is bounded relative to \(S\), fix \(V\) and
apply (TW1) to the supported letter and \(U\)'s integral
cofactor. If \(U\) is large, first fix \(V\) in a
preliminary count. Every large lock on \(U\) costs at
least \(d-w_V\). Greedy packing with cap \(d-w_V-\kappa\)
finds one large group \(G_U\) of cost
\(\beta_U\ge d-w_V-\kappa\). If every group above this
cost were bounded relative to \(S\), set those
boundedly many groups aside. Every two consecutive
greedy packs of the other groups together exceed the
cap, so only boundedly many cheap packs are needed.
Their individual log-heights would be \(O(S)\), in
contradiction to \(U\) being large.
The cost of fixing its complement on \(U\) is at most
\(w_U-\beta_U\le d-2\eta+\kappa\).
If \(V\) is bounded, fix the complement and use (TW1)
on \(G_U\), the supported letter and the integral
cofactor on \(V\). If \(V\) is large, first fix that
complement in a preliminary count. Every large lock
on \(V\) releases both \(G_U\) and the supported letter,
so costs at least \(d-w_U+\beta_U\). Select a large
group \(G_V\) of cost at least \(d-w_U+\beta_U-\kappa\)
by the same greedy step. Now lock each selected group
in turn and fix the other selected group's complement.
The released ratio has the bounded-root description.
This gives \(d\le\beta_U+w_V-\beta_V\) and
\(d\le\beta_V+w_U-\beta_U\), contradicting
\(w_U+w_V\le 2d-2\eta\).

All concentrated and selected lists in this proof have
a bounded number \(M\) of new roots determined by
\(A,\kappa\) and the group bounds. Choose \(2\xi<\xi_0\)
for the worst such \(M\) in (TW1), so every terminal
bound is valid simultaneously. All lock choices and
the number of original positions are fixed in the
height limit at each budget. The exponent margins
absorb dyadic and prime-colouring losses. This proves
the stated uniformity.
[Bounded-group proof with arbitrary coefficients (MCB2)](ABC_MULTIPLE_SHORT_COFACTOR_BOUNDARY.md#bounded-group-dictionaries-with-arbitrary-fixed-coefficients)

## Consequences for progressions and balanced paths

**TW4.** Given a fixed positive numerical lower height
bound on the tall endpoints and a bounded radical cap
in base \(n\), unit progression blocks of any step \(D\)
and index width at most \(n^{L^{-\xi}}\) have occupancy
\(O(n^\eta X^\varepsilon)\) for some \(\xi>0\) depending
on the prescribed normalized allowance \(\eta\).
Fix one endpoint \(T\), enumerate \(g=\gcd(T,V)\mid T\)
at divisor cost, and use \(g\mid |V-T|/D\) for the unit
condition. The primitive equation has fixed shift
\(T/g\), coefficient \(D\) on the integral index
\(|V-T|/(Dg)\), and the other term \(V/g\) with its
original radical cap. Since \(g\le n^{L^{-\xi}}\),
the primitive maximum has a fixed positive numerical
height floor for small \(L\) if \(\xi<1\). Apply (TW2)
with zero compressed roots. This proof retains any
original constraints on the endpoints.
[Fixed-width predecessor](ABC_MULTIPLE_SHORT_COFACTOR_BOUNDARY.md#arbitrary-span-progression-occupancies)

In (RT8), the negligible path count at any fixed
\(u>0\) likewise holds with primitive quotient at most
\(n^{L^{-\xi}}\). Fix the first original row and its
gcd record and apply (TW2) to the original second row
with its common complete value fixed. Its maximum is
still at least \(X/2\), its marked variable part after
the gcd record is at most \(n^{L^{-\xi}}\), and the
three original radical caps are all retained. Choose the
normalized allowance smaller than \(u\).
On the remaining paths (TW3) permits sufficiently thin
valuation-window groups in (RT7): take the lower height
parameter \(\theta\) smaller than the path cutoff exponent
\(\xi\), then choose the window exponent for (TW3).
Thus no separate low-position charge is needed in that
range even for groups spanning a growing number of exact
orders. These are fixed-positive-multiplicity statements.
They leave the broad continuous profiles and the count
of occupied cells near zero multiplicity open.
[Weighted codegrees and small paths (RT7–RT8)](ABC_BALANCED_FIXED_SHIFT_REVERSE_TRANSFER.md#arbitrary-weighted-shifts-at-radical-scale-height)

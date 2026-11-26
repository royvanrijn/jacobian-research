# Radical amplification on popular additive shifts

A shift shared by many common-target cherries in a near-full
complete-label graph has a larger radical cost than direct radical
enumeration detects. The gain uses a primitive difference triple and
all three original radical caps. In the minimum nonatomic sector, the
same shifts cannot carry a positive numerical weight at arbitrarily
large scaled valuation orders.

## A bound for rich shifts with a radical cap

Let \(n=X^L\). Use positive primitive triples \(a+b=c\) in
\(X/2<c\le X\), and choose two distinct roles \(U,V\).
Represent a triple in signed coordinates \((a,b,-c)\).
For \(h>0\), let \(M_h\) count the unordered pairs of edges
that share the remaining complete letter \(T\), with opposite
signed differences \(h,-h\) in \(U,V\). These are the
cherries used in the translation-cube argument. Put
\(c_0=1/6\), and assume each complete part has at most
\(n^{c_0+\epsilon_0}\) labels and each original letter has
radical at most \(n^{1/3+\epsilon_0}\).

For each fixed \(B_0,\Delta>0\) there are
\(\epsilon_0,L_0>0\) such that at all sufficiently small
fixed \(0<L<L_0\), uniformly in \(0\le b\le B_0\),
\[
 \#\{h>0:M_h\ge n^{1/3-\epsilon_0},\
                    \operatorname{rad}(h)\le n^b\}
     \ll_{L,B_0,\Delta}
             n^{\min\{b,\,1/18+b/3\}+\Delta}.                   \tag{PSA1}
\]
The same bound holds for each of the three choices of the
varied pair. In particular, a shift family of exponent
\(1/6-o(1)\) in the minimum nonatomic reduction cannot have
radical cap \(n^{1/3-\delta}\) for any fixed \(\delta>0\).
The exponent \(1/3\) here is twice the shift dictionary
exponent required by direct radical enumeration.
[Binary shift allocations and cherry extraction](ABC_ADDITIVE_TRANSLATION_CUBES.md#binary-shift-allocations-between-prescribed-letters)

### Removing overlap on the fixed shift support

For any prescribed \(\xi>0\) sufficiently small, discard
cherries where the radical on either varying letter
supported by the primes of \(h\) exceeds \(n^\xi\).
Its two endpoint values have the same divisibility by primes
of \(h\), so either endpoint can be used to test this
condition. Fix the shared original complete \(T\), and
use (FL1) on one original branch. The outside radical cap
on the offending letter is \(n^{1/3+\epsilon_0-\xi}\).
Choose the branch with smaller signed \(U\)-coordinate,
so that branch and \(h\) determine the cherry uniquely.
Summing over the two roles and the at most
\(n^{c_0+\epsilon_0}\) values of \(T\) bounds the discarded
count by \(O(n^{1/3-\xi/4})\) after taking
\(\epsilon_0\), the (FL1) coefficient allowance and the
inner height allowance small enough. The bound is uniform
in \(h\), since \(h\le X\). At least half the cherries
of every rich \(h\) remain for large height.
[Fixed-letter count with prescribed prime pools (FL1)](ABC_FIXED_LETTER_RADICAL_CAPS.md#a-pointwise-count-with-prescribed-prime-pools)

### A tall primitive difference triple

On a retained cherry write the two complete positive values
of its first varied letter as \(u,u+h\), and those of its
other varied letter as \(v,v+h\), with \(u,v>0\). Set
\[
                       g_u=(u,h),\qquad g_v=(v,h).
\]
Each divides \(h\), and they are coprime. If a prime divided
both, it would divide both varied letters on each of the two
original primitive triples. Thus \(g_ug_v\mid h\).
Both larger endpoints are at least \(h\), and at least one
is at least \(X/4\): every original edge has its target and
at least one summand at those heights, and one of those roles
belongs to any chosen pair of varied roles. It follows that
\[
             \frac{u+h}{g_u}\frac{v+h}{g_v}
                       \ge\frac{(u+h)(v+h)}h\ge X/4 .
\]
Choose the first role with its quotient at least \(\sqrt{X}/2\),
in a fixed ordering. In the first role record the triple
\[
                  \left(\frac u{g_u},\frac h{g_u},
                                \frac{u+h}{g_u}\right),
\]
and use the analogous triple in the second role. The recorded
triple is positive and primitive and has height between
\(\sqrt{X}/2\) and \(X\). Using magnitudes makes the same
construction valid when one varied role was the signed
target \(-c\). The multiplier in the chosen role has
radical at most \(n^\xi\).

For fixed \(h\), the selected role and its smaller complete
value have at most \(n^{c_0+\epsilon_0}\) preimages as
cherries, since the other role has at most that many labels.
Indeed the smaller value and \(h\) determine both signed
values in the selected role, whose ordering is fixed by the
role signs and the orientation of the cherry. The other role
determines \(T\) by addition. Each rich shift therefore
yields \(\gg n^{c_0-2\epsilon_0}\) such role and value
records. Conversely, fixing a primitive recorded triple,
its multiplier and its role recovers \(h\) and the smaller
value exactly. Across all shifts the number of multipliers
with radical at most \(n^\xi\) is \(O(n^{\xi+\varepsilon'})\)
for any fixed inner allowance \(\varepsilon'>0\), by
enumeration of the radical followed by supported powers.
[Uniform fixed-support enumeration](ABC_EXCEPTIONAL_SET_AND_POWERFUL_TRIPLES.md#a-radical-count-contradiction-to-the-claimed-smooth-number-lower-bound)

The recorded primitive triple has radical at most
\(n^{2/3+b+2\epsilon_0}\). Apply the limiting one-third
count to its height shell \(Y/2<\max\le Y\). Since
\(X^{1/2}/2\le Y\le 2X\), the required relative budget
\((2/3+b+2\epsilon_0)L\log X/\log(Y/2)\) is \(O_{B_0}(L)\).
Round this budget up to the next point of a finite grid of
mesh \(\rho L\), for fixed small \(\rho>0\) and at the fixed
original \(L\). Its product with \(\log Y/\log X\) is
at most \((2/3+b+2\epsilon_0+\rho+o_X(1))L\).
Thus the one-third estimate at each fixed grid budget and
the dyadic shell summation give exponent
\((2/3+b)/3\) in base \(n\), up to arbitrarily small
fixed margins. Grid size and constants are uniform in
\(b\in[0,B_0]\). Multiply by the number of scalings and
divide by the number of records per shift. Choose \(\xi\)
and all other allowances small compared with \(\Delta\).
This gives exponent \(1/18+b/3+\Delta\). Directly
enumerating the shift radical and then its powers gives
the cap \(b+\Delta\).
[Limiting small-budget coefficient \(1/3\)](ABC_OMITTED_COORDINATE_GRAPH.md#a-limiting-small-budget-coefficient-of-13)

The cap in (PSA1) may also be imposed only outside the primes
of one prescribed integer \(P\le X^{K_0}\), with \(K_0\)
fixed, if the raw radical of \(h\) is at most \(n^{B_0}\).
The bound remains \(n^{\min\{b,1/18+b/3\}+\Delta}\), with
constants also allowed to depend on \(K_0\). The discarded
overlaps and primitive records are unchanged. On each derived
height shell apply (EB4), pricing the primes of \(P\) at zero
and all others at their ordinary logarithmic costs. The raw
relative budget can be rounded up on a fixed grid as above
using \(2/3+B_0+2\epsilon_0\). The actual price in base \(n\)
is at most \(2/3+b+2\epsilon_0\). Since all derived heights
are at least \(\sqrt X/2\), the number of primes in the pool
is \(O_{K_0}(\log Y/\log\log Y)\) at shell height \(Y\).
Together with the prime harmonic sum this verifies the
partition condition of (EB4) uniformly. Its budget allowance
is \(O_{B_0}(\eta L)\) for coefficient tolerance \(\eta\),
which can be chosen small compared with \(\Delta\).
Direct enumeration at exponent \(b\) follows by first
fixing the pool-supported part at subpower cost and then
enumerating the outside radical and its powers.
[Effective radical counts and fixed prime pools (EB4, EB5)](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)

## Valuation mass on a large dictionary of rich shifts

Use minimum nonatomic balanced complete-label boxes in the
height-first limit, with part caps \(n^{c_0+o(1)}\).
Suppose there are \(n^{\tau-o(1)}\) rich shifts with
\(\tau>0\), \(M_h\ge n^{1/3-o(1)}\), and a bounded
normalized radical cap \(B_0\) on each \(h\). After anatomic
refinement, suppose the exact scaled-degree measures of
these shifts converge to \(\nu_h\) on the retained class.
Then
\[
                \tau\le \|\nu_{h,\mathrm{na}}\|,
                    \qquad \tau<\nu_h((0,\infty)).              \tag{PSA2}
\]
Such refinements are available on a subfamily of the same
height exponent. In finitely many exact valuation layers
take dyadic root bins. A tail of orders at least \(D\) has
radical at most \(X^{1/D}\). The truncation can be chosen
arbitrarily large at each fixed budget, so the normalized
tail radical error tends to zero. The bounded first moment
and raw mass then allow weak limits along subsequences.
[Height-first valuation conventions](ABC_SECANT_ENDPOINT_TRANSFER.md#the-height-first-setting)

Use the primitive difference triples in the proof of
(PSA1), allowing the overlap tolerance to decrease to
zero slowly in the budget limit. At each fixed tolerance
the uniform FL1 error discards a power-smaller count.
The recorded triples have exponent at least \(c_0+\tau-o(1)\).
One may fix the selected role and one numerical height
shell at subpower cost. The complete dictionary of their
summand \(u/g_u\) in that role has cap \(n^{c_0+o(1)}\),
because \(u\) comes from one original part and all
multipliers have vanishing normalized radical cap.
Pigeonhole that exact free-opposite value. The resulting
actual primitive fixed-letter fibre has exponent at
least \(\tau\), with \(h/g_u\) free. The shift and its
quotient have the same limiting weak masses because the
valuation changes are supported on \(\operatorname{rad}g_u\)
of vanishing normalized mass.

The quotient height shell has numerical exponent bounded
away from zero. At each fixed \(L\), take an inner subsequence
of chosen shells along which \(\log_XY\) tends to a number
\(H_L\in[1/2,1]\). Use the fixed relative budget unit
\(L/H_L\) at that inner step. Then
\(\log(Y^{L/H_L})=(L+o_X(1))\log X\):
the normalized count exponent and radical masses are
unchanged. All scaled valuation degrees are multiplied
by \(1/H_L\), which has a positive finite limit on an
outer subsequence. This preserves nonatomic mass and mass
at positive scales. Apply (FL2) to the fixed-letter fibre with
the two shifted quotients and the bounded-radical
shift quotient all included. This gives (PSA2).
[Strict positive-scale bound on the two free letters (FL2)](ABC_FIXED_LETTER_RADICAL_CAPS.md#strict-positive-scale-inequalities-on-the-two-free-letters)

## Numerical rigidity at the critical shift radical

In minimum-mass balanced complete-label boxes, suppose the
shift family has exponent \(c_0-o(1)\), each shift is rich
at exponent \(1/3-o(1)\), and its radical is at most
\(n^{1/3+o(1)}\). The shifts in any such class of full
exponent satisfy, in the limiting sense,
\[
        \log_X h\longrightarrow 1,\qquad
        \limsup\sum_{\substack{p^e\parallel h\\ eL\ge R}}
                  e\log_X p=0\quad(R>R_*),                   \tag{PSA3}
\]
where \(R_*\) is the absolute cutoff in (VM2). Here the
assertions exclude a subfamily of full exponent on which
either a fixed positive numerical height gap or a fixed
positive tail weight is imposed. The statement assumes the
critical raw shift cap.

Indeed the original minimum-mass profiles have nonatomic
mass \(c_0\) on each letter. Equations (VM1) and (VM2) give
numerical height one on each and bounded scaled support
with no escaped first moment. For the retained cherries,
a multiplier supported on \(h\) with vanishing normalized
radical mass consequently has numerical weight \(o(1)\).
To see this, split its prime powers at a fixed continuity
cutoff above the bounded support of the original letter.
The weight below the cutoff is bounded by the cutoff
times its normalized radical mass. The weight above the
cutoff tends to zero by convergence of the first moments.

Now the primitive difference triples in (PSA1) have
exponent \(1/3-o(1)\) and radical cap \(n^{1+o(1)}\).
Their numerical height tends to one in base \(X\).
Both source copies are tall and have the same limiting
nonatomic mass \(c_0\), because division by the multiplier
changes only primes of vanishing normalized radical mass.
Apply (VM1) to these two tall letters in the derived
endpoint. Their excess is \(2c_0-1/3=0\), so the third
letter \(h/g\) also has numerical height one. Apply (VM2)
with the same pair. It gives no escaped first moment and
bounded scaled-degree support for \(h/g\), with cutoff
\(R_*\). Adding back \(g\) has numerical cost \(o(1)\)
for the height assertion. To restore the tail, take
\(R_*<R'<R\). For a prime with \(v_p(h)L\ge R\) but
\(v_p(h/g)L<R'\), at least the fraction \(1-R'/R\) of its
weight in \(h\) belongs to \(g\). The total such weight is
\(o(1)\). The remaining primes contribute at most their
\(h/g\) tail above \(R'\) plus the \(o(1)\) multiplier
weight. That tail tends to zero by (VM2); the derived scale
agrees with the original one in the limit, leaving room
between \(R'\) and \(R_*\). If a height-gap or tail
subfamily with fixed positive margins had no uniform
normalized saving for arbitrarily small fixed budgets and
cap tolerances, it would give the violating height-first
sequence just excluded. This gives sector savings with
the width of the raw shift-cap tolerance chosen sufficiently
small relative to the stated margins.
[Height gap and high tails from pair excess (VM1, VM2)](ABC_VARIABLE_MESH_DEGREE_GROUPING.md#bounds-from-the-pair-excess)

There are absolute \(\kappa,\zeta>0\) such that, at the critical
cap, no full-exponent cohort of rich shifts can have normalized
radical mass less than \(\zeta\) below \(3-\kappa\), or less than
\(\zeta\) above \(3+\kappa\), with sufficiently small fixed cap
and richness tolerances. Otherwise take a diagonal inner-then-outer
sequence as these cutoffs approach \(3\). By (PSA1) the shift
measure has mass \(1/3\), and by (PSA3) it has first moment \(1\).
If its support is confined to either half-line ending at \(3\),
it must equal \((1/3)\delta_3\). That contradicts (PSA2) at
\(\tau=1/6\). In particular a full-exponent shift cohort at
the critical cap cannot consist of \(\lceil(3-\kappa)/L\rceil\)-full
integers. This exclusion is for the shifted dictionary, independently
of the original letters' grids.

## Popular shifts in short numerical intervals

Under the cap and richness assumptions of (PSA1), for every
fixed \(\eta>0\) there is a fixed \(\theta>0\) such that
every interval \(I\) of length at most \(X^\theta\) contains
\[
       O_{L,\eta}\!\left(n^{\eta+3\epsilon_0+\varepsilon}\right)
\]
rich shifts, for sufficiently small fixed \(L\) and any
fixed \(\varepsilon>0\). No shift radical cap is needed here.

Assign each cherry to a varied role whose larger complete
value is at least \(X/4\), using a fixed tie rule.
Fix the smaller value \(u\) in that role. When \(h\in I\),
the larger value \(u+h\) is in the interval \(u+I\) in
a shell at height comparable to \(X\). The short-interval
radical count (SI2) gives at most \(n^{\eta+\varepsilon}\)
larger endpoints per anchor, after shrinking \(\theta\)
and the inner allowances if necessary. There are at most
\(n^{c_0+\epsilon_0}\) smaller endpoints and at most that
many choices in the other varied role. The resulting upper
bound for \(\sum_{h\in I} M_h\), divided by the richness
threshold, proves the assertion.
[Arbitrarily small short-interval normalized counts (SI2)](ABC_SHORT_INTERVAL_RADICAL_FILTER.md#arbitrarily-small-normalized-counts-in-shorter-intervals)

The bounds constrain both the radical and the numerical
distribution of rich shifts. The arithmetic graph argument
still allows a full dictionary of shifts whose radicals are
larger than the critical cap.

# Cross-target residue classes from conditional divisor records

A popular complete target gives many distinct witnesses on either summand.
If a powered divisor is prescribed on the other summand, the witnesses
share the target's residue at that divisor. Conditional dictionary
bounds supply lists of such divisors common to many targets.
Their least common multiple can make the residue classes and a short
numerical interval determine the target. This gives a strict upper
bound on the same-target scalar image on specified sectors.

## Divisor records on a complete-target bin

Let \(n=X^L\), and consider eligible positive primitive rows
\(A+B=C\) with \(X/2<C\le X\). The eligibility condition may
include the original small radical cap and any exact least-repair
cutoff \(H_M(C)\le Y\). Restrict to a bin where the multiplicity
\(m(c)\) of every occurring complete target \(c\) lies in
\([n^{u-\tau},n^{u+\tau}]\). The radical caps on the three letters are
\(n^{r_A+\tau},n^{r_B+\tau},n^{r_C+\tau}\), where the three
raw caps are fixed. At each sufficiently small fixed
\(L\), the height limit precedes the passage to small budgets.

Fix \(\ell\ge1\) choices of marks, allowing the same choice more
than once. For choice \(i\), let \(W_i\) be the summand with a
marked full divisor \(F_i\) and let \(V_i\) be the other summand.
Prescribe the following bounds in the eligible row family.

- The complete \(V_i\)-dictionary has size at most \(n^{\alpha_i+\tau}\).
  The complete value \(W_i\) is at most \(X^{h_i+\tau}\).
- The marked divisor is the product of the full prime powers on its
  selected support. Its radical and numerical size satisfy
  \[
       n^{t_i-\tau}\le\operatorname{rad}(F_i)\le n^{t_i+\tau},
       \qquad F_i\ge X^{w_i-\tau}.
  \]
  Each of its nonzero scaled valuations is at most \(S_i\),
  so \(Lv_p(F_i)\le S_i\).

The marked root vector must be a group admissible in (FS1).
One can always take a whole original singleton position.
One can also use a group of the allowed bounded size or shrinking-window
type, or a selected subpacket if splitting the selected and unselected
root factors leaves such an admissible group. The complement on
\(W_i\) is charged at its actual remaining radical cap
\(n^{r_{W_i}-t_i+2\tau}\). This uses full prime powers and disjoint
selected and unselected supports. Prescribing \(F_i\) prescribes
its grouped root record up to \(X^{o_X(1)}\) choices. For exact
original valuation positions the prime factorization of \(F_i\)
determines those positions outright. Alternatively the record can
include the bounded-position allocation of its primes. The latter
has subpower multiplicity on a fixed support at fixed truncation.

One may take \(F_i=1\) with \(t_i=w_i=S_i=0\). That mark
serves only to supply a numerical interval. For nontrivial
marks the stated hypotheses concern the marked group, and the
other original degree positions remain unrestricted. The
(FS1) estimates below use the original numerical root-height
bounds and the prescribed group dictionaries.
[Fixed-shift dictionary count (FS1)](ABC_FULL_SHIFT_DICTIONARIES_AND_CODEGREES.md#a-fixed-shift-count-without-a-radical-cap-on-the-shift)
· [Fixed-letter bound with prescribed prime pools (FL1)](ABC_FIXED_LETTER_RADICAL_CAPS.md#a-pointwise-count-with-prescribed-prime-pools)

## Upper bound for the common target dictionary

Write
\[
 \beta_i=r_{W_i}-2u,\qquad
 \gamma_i=\min(t_i,\beta_i),\qquad
 h_*=\min(1,h_1,\ldots,h_\ell),\qquad
 E=\sum_{i=1}^{\ell}\bigl(\alpha_i-(u-t_i)_++\gamma_i\bigr).
\]
Assume \(\beta_i\ge0\). When any \(\beta_i<0\) by a fixed
margin, (FL1) already excludes such a multiplicity bin for small
enough allowances. A trivial mark uses \(\gamma_i=0\).
If for each nontrivial mark after the first one can choose
\(\zeta_i>\beta_i\) with
\[
         \sum_{i=1}^{\ell} w_i
             -\sum_{\substack{i\ge2\\F_i\ \mathrm{nontrivial}}}
                       S_i\zeta_i
             > h_* ,                                           \tag{TC1}
\]
then, for every fixed \(\eta>0\), the number of occurring complete
targets is
\[
                   O\!\left(n^{E+\eta}X^\varepsilon\right).    \tag{TC2}
\]
Here the positions of the marks are fixed and may be reordered
before applying (TC1). The assertion holds for sufficiently small
fixed \(\tau>0\) and all sufficiently small fixed \(L>0\), depending
on the fixed parameters, the strict margins and the allowed group
data. The inner height losses can be made arbitrarily small at each
such \(L\). In particular \(\sum_i w_i-\sum_{i\ge2}S_i\beta_i>h_*\)
with fixed slack suffices, taking \(S_i=0\) at trivial marks.

One may omit the numerical intervals and use \(h_*=1\) in (TC1).
Keeping the intervals costs only a bounded factor. A trivial mark
on a short summand \(W_i\) supplies that summand's height exponent
in \(h_*\) at target-dictionary cost \(\alpha_i-u\).
A nonempty multiplicity bin necessarily has
\(\alpha_i\ge u-o(1)\), because the opposite values in a single
fibre are distinct.

### Conditional multiplicity of a packet

Let \(m_i(c,F)\) be the number of original eligible rows on
\(c\) on which the selected record for choice \(i\) is \(F\).
Radical and fixed-support enumeration give at most
\(n^{t_i+\tau}X^{o_X(1)}\) records. Thus the rows belonging to
records with \(m_i(c,F)<n^{u-t_i-\theta}\) have count at most
\(n^{u-\theta+\tau}X^{o_X(1)}\), for a fixed \(\theta>0\).

Records with \(m_i(c,F)>n^{u-t_i+\beta_i+\theta}\) form a
dictionary of size at most \(n^{t_i-\beta_i-\theta+\tau}\),
by the upper multiplicity for \(c\). If this exponent is
negative there are no such records for large \(X\). Otherwise,
charge this dictionary on the marked group and the raw remaining
radical cap on its complement. The \(W_i\)-dictionary charge
summed over groups is at most
\[
         t_i-\beta_i-\theta+\tau
               +r_{W_i}-t_i+2\tau
                  =2u-\theta+3\tau
\]
in base \(n\), up to inner allowances. To implement these caps,
fix the exponents and the truncation at each fixed \(L\).
Prime powers above an exact-degree cutoff greater than \(S_i/L\)
are outside the mark. Their support has arbitrarily small
normalized radical at a sufficiently large cutoff, and their
full product can be fixed at arbitrarily small cost. Dyadic root
boxes for the remaining positions give the separate raw caps of
the complement with their sum bounded by its radical cap plus
arbitrarily small allowance. The same procedure applies if the
mark is a specified admissible split of these positions.

The positive target \(c>X/2\) supplies a full divisor of fixed
positive numerical weight and fullness order comparable to
\(1/L\). Indeed remove the powers whose scaled order is below
a sufficiently small fixed positive constant. The original
target radical cap makes the numerical weight removed smaller
than, for example, \(1/4+o(1)\). Thus (FS1), with the complete
target fixed, bounds the heavy-record rows by
\(n^{u-\theta/2+O(\tau)+o(1)}X^\varepsilon\). This is
power-smaller than the fibre for small enough \(\tau\) and
inner allowances. Hence a set of \(n^{u-O(\tau)-o(1)}\)
witnesses is available for choice \(i\) with
\[
        m_i(c,F)\ge n^{(u-t_i)_+-O(\theta+\tau)}
\]
and with upper record multiplicity bounded by the stated
heavy threshold and by \(m(c)\). For a trivial mark this follows
directly from the multiplicity bin. The bounds on \(m_i(c,F)\)
refer to the original eligible witnesses on \(c\), before any
overlap exclusions.

### Many targets with one packet list

Fix \(c\) and select the marks in their prescribed order.
At step \(i\ge2\) let \(P\) be the product of all earlier
selected divisors. It has polynomial height since \(\ell\)
is fixed. Among the current witnesses discard rows with
\[
        \operatorname{rad}\gcd(F_i,P)\ge n^{\zeta_i}.
\]
For a nontrivial mark this event cuts the outside-\(P\)
radical cap on \(W_i\) to \(n^{r_{W_i}-\zeta_i+\tau}\).
It is empty if this cap is below one. Otherwise (FL1)
bounds the discarded rows by
\(n^{(r_{W_i}-\zeta_i)/2+O(\tau)+o(1)}X^\varepsilon\),
uniformly in \(P\). Because \(\zeta_i>\beta_i\), this loss
is power-smaller than \(n^u\). If \(F_i=1\), no overlap
condition is needed.

After the exclusions, divide the remaining count by the
upper record multiplicity at step \(i\). This gives at least
\(n^{(t_i-\beta_i)_+-O(\theta+\tau)-o(1)}\) possible next
records, uniformly in the earlier record list. The trivial
lower bound of one handles \(t_i-\beta_i\le0\).
Counting the resulting lists over all complete targets,
and using the at most \(n^{\sum_i t_i+O(\ell\tau)+o(1)}\)
possible lists, one list is shared by a set \(\mathcal C_*\)
of targets of size at least
\[
       |\{c\}|\,n^{-\sum_i\gamma_i-O_\ell(\theta+\tau)-o(1)} .
                                                               \tag{TC3}
\]
Different choices of marks may use different original rows
on a given target. This step imposes no independence
assumption on the summands or on their prime factors.

### Residue and interval charts for the fixed list

For each fixed \(F_i\) in the list, every \(c\in\mathcal C_*\)
has at least
\(K_i=n^{(u-t_i)_+-O(\theta+\tau)}\) distinct opposite
witnesses \(V_i\equiv c\bmod F_i\). They all lie in
\([c-X^{h_i+\tau},c]\). Partition the opposite dictionary
by residue modulo \(F_i\) and by integer intervals of
length \(\lceil X^{h_i+\tau}\rceil\). At least a fixed
fraction of the \(K_i\) witnesses belongs to one adjacent
interval and to the common residue. The number of
residue-interval pairs containing that many dictionary
values is \(O(n^{\alpha_i+\tau}/K_i)\). Each gives one
specified residue for \(c\) and an interval for \(c\) of
length \(O(X^{h_i+\tau})\).

Let \(Q=\operatorname{lcm}(F_1,\ldots,F_\ell)\).
For a nontrivial step \(i\ge2\), the upper bound on scaled
orders and the overlap exclusion give
\[
 \log_X\gcd\!\bigl(F_i,\operatorname{lcm}(F_1,\ldots,F_{i-1})\bigr)
                         \le S_i\zeta_i .
\]
Thus \(\log_X Q\) is at least the left-hand side of (TC1)
minus \(O(\ell\tau)\). By (TC1) it exceeds the shortest
available interval exponent, or the original height
exponent, by a fixed margin. For any chosen list of
residue-interval pairs, two compatible integers would
differ by a multiple of \(Q\), so at most one target
occurs. Counting these lists gives
\[
        |\mathcal C_*|
              \ll n^{\sum_i(\alpha_i-(u-t_i)_+)
                           +O_\ell(\theta+\tau)+o(1)} .
\]
Combine this with (TC3), choosing \(\theta,\tau\) and all
dictionary allowances small in terms of \(\eta\). This
proves (TC2).

More precise valuation information can replace \(S_i\zeta_i\).
It suffices to bound the maximum numerical weight in \(F_i\)
of the full powers on a subset of its support of normalized
radical at most \(\zeta_i\). The shared prime powers removed
in the least common multiple are bounded by that weight.
This follows directly from the same overlap argument.

## Unit progression cells with a large small-radical step

**TC4.** Fix \(K,q_0,r_0,w_0,\eta>0\). There is a \(\vartheta>0\), depending only on these fixed
constants, with the following property for
all sufficiently small fixed \(L>0\). If a progression step
\(Q\) satisfies \(X^{w_0}\le Q\le X^K\) and
\(\operatorname{rad}(Q)\le n^{q_0}\), then any unit residue
class modulo \(Q\) in an interval of length \(Q X^\vartheta\)
contains at most \(O(n^\eta X^\varepsilon)\) positive
integers \(c\le X^K\) with \(\operatorname{rad}(c)\le n^{r_0}\).
The count is uniform in the interval, step and residue.

Suppose this assertion fails for \(\vartheta\) and \(L\)
arbitrarily small. Work on a sequence with
\(\vartheta\to0,\ L\to0\) and at each fixed pair take the
height limit first. Box the exact valuation layers of
\(c\). At an exact-degree cutoff that is fixed at each
budget and arbitrarily large relative to \(1/L\), the
high tail has arbitrarily small normalized radical and
can be fixed at that cost plus an inner subpower loss.
Fix it and pass to a dyadic box of the remaining roots.
We can keep the total outer discard and record costs
below \(\eta/4\). The scaled valuation measures in
this box are tight by the upper bound on numerical
height. Extract a convergent subsequence.

If the limit has positive nonatomic mass, some compact
positive scaled-degree band carries a positive amount
of it. Select continuity windows in the band containing
a positive nonatomic mass with total mass below
\(\eta/4\). The full product \(D\) of the corresponding
original layers has fixed positive numerical weight and
cost below \(\eta/4\) in the exponent in base \(n\).
These facts follow by bounding the smallest scaled
degree below and enumerating the original root vectors
in the chosen windows. Since \(D\mid c\) and the
progression is a unit class, \((D,Q)=1\). Thus \(D\)
fixes the progression index modulo \(D\), permitting at
most one \(c\) in the interval once \(D>X^\vartheta\).
This counts the chosen box strictly below the presumed
exponent \(\eta\).

Otherwise the limit is purely atomic. Enclose finitely
many positive atoms and the zero atom in slowly
shrinking windows, missing total normalized mass less
than \(\eta/4\). Fix the remaining exact factors at
that cost plus subpower terms. In each positive window
group the squarefree roots at a representative exponent
of order \(1/L\). The small exponent discrepancies and
the zero window give one rational correction of numerical
height \(X^{o(1)}\). The grouped roots have disjoint
supports and normalized height caps bounded in terms
of \(r_0\).

There is a full divisor of \(Q\) of numerical weight at
least \(w_0/2\) and fullness order comparable to \(1/L\).
Indeed for fixed \(\kappa>0\) sufficiently small in
terms of \(w_0,q_0\), removing its prime powers of order
below \(\lfloor\kappa/L\rfloor\) loses at most
\(\kappa q_0<w_0/2\). Divide the unit residue condition
by the fixed factor of \(c\) modulo this full divisor.
The amplified wide-cofactor bound (GW) counts the
remaining values at arbitrarily small normalized
exponent, using the bounded number of grouped narrow
roots and the short rational correction. With no
positive atom, the correction is directly recovered
from the modulus. The selected costs again sum to
less than \(\eta\). The contradiction gives the
claimed fixed \(\vartheta\) and budget cutoff. This proof
uses the full step in the atomic case. Earlier adaptive
locks give the nonatomic case without a fullness
hypothesis on the step.
[Stopping a full divisor (PI4–PI5)](ABC_POWERED_PROGRESSION_OCCUPANCY.md#stopping-a-full-divisor-above-the-block-width)
· [Amplified wide-cofactor bound (GW)](ABC_WIDE_COFACTOR_CONGRUENCES.md#a-bounded-number-of-small-roots)

In (TC2), write \(w_{\mathrm{list}}\) for the left-hand
side of (TC1). Given any fixed positive lower bound
\(w_0<w_{\mathrm{list}}\), there is therefore a relaxation
of the numerical condition to
\[
                     w_{\mathrm{list}}>h_*-\vartheta/2 .
\]
Apply (TC4) to each compatible residue and interval
chart after fixing the packet list. Its least common
multiple has polynomial height, weight at least \(w_0\),
and radical at most \(n^{\sum_i t_i+O(\ell\tau)}\). All
targets are units modulo the step: each marked divisor
is on a summand coprime to that target. Thus each cell
now costs an arbitrarily small normalized exponent
rather than at most one integer. Choose that allowance
as a fixed part of the desired error in (TC2) before
using the corresponding \(\vartheta\), then take
\(\tau\) small enough for both strict numerical margins.
The resulting exponent is still \(E+\eta\). This
weakens the required numerical weight, with the
same need for economical opposite-letter dictionaries.

## Consequence for least-repaired scalar moments

**TC5.** Fix \(k\) and use the scalar records of (SR3) on any chart
where they are defined. The number of target radicals
represented in a bin is at most its number of complete
targets. Keeping one complete representative and at most
\(\lceil n^u\rceil\) eligible witnesses per radical, the
same-target \(k\)-fold image has upper exponent \(E+ku+\eta\)
in base \(n\), with the usual inner height allowance.
The same exponent holds using all witnesses in a bin of
sufficiently small fixed width, after adding \(k\tau\) to
the exponent allowance. The target radical remains subject
to its original cap and to the exact least-repair cutoff.
The count uses the unique least multiplier of each primitive quotient.

At the balanced benchmark this supplies a power-saving
sector of (SR3) whenever \(E<1/3-u\) by a fixed margin.
Indeed the complete-target bin has row exponent at most
\(E+u\) before the arbitrarily small allowances. The
method needs economical complete dictionaries for the
opposite roles used in the marks. It does not show such
dictionaries exist throughout a balanced endpoint.
It also does not force a packet collection satisfying
the numerical inequality (TC1).
[Scalar image and target multiplicity (SR3)](ABC_SCALAR_MASK_TARGET_MOMENTS.md#a-scalar-image-for-weighted-target-occurrence)
· [Least-repair multiplicity and its range (FL5–FL8)](ABC_FIXED_LETTER_RADICAL_CAPS.md#size-of-a-minimally-repaired-target-radical-fibre)

For a numerical example take \(u=0.15\), \(r_A=r_B=1/3\),
and complete dictionary exponents \(0.181\) for \(A\)
and \(0.175\) for \(B\). Suppose \(A\) has height at most
\(X^{0.48+o(1)}\), and mark on \(B\) an admissible packet
of radical cost \(0.06+o(1)\) and numerical weight at least
\(0.50-o(1)\). Use that packet first and a trivial mark
on \(A\) second. The latter supplies the short interval,
and there is no overlap loss. The target exponent is at
most
\[
  (0.181-(0.15-0.06)+(1/3-2\cdot0.15))
       +(0.175-0.15)=0.149333\ldots
            <1/3-0.15 .
\]
Both complete dictionary exponents in the example are
above the minimum \(1/6\). The parameters describe a
sufficient sector for a count saving.

The closest original packing argument fixes complete
anchors and derives several affine forms in one shared
query. The packet argument above fixes full divisors
instead. The heavy residue and interval cells use the
multiplicity through a common target, which permits
intermediate multiplicities and opposite dictionaries
above the minimum scale.
[Mixed packing with fixed complete anchors](ABC_CROSS_ROLE_COMPACT_PACKING.md#a-mixed-packing-inequality-in-actual-additive-graphs)

## A radical-cost limit of the residue-list bound

**TC6.** Work at the balanced limit and let the row exponent
be \(1/3\). Each allowed opposite complete dictionary must
have \(\alpha_i\ge1/6\), as a fixed complete value has fibre
exponent at most \(1/6\) by (FL1) with the balanced caps.
Set \(r_{W_i}=1/3\) in the limiting expression for \(E\).
For every \(u\in[0,1/6]\), the strict exponent inequality
\(E<1/3-u\) implies
\[
                             \sum_i t_i<1/6 .
\]
In particular, if all marked scaled orders are bounded above
by \(S<6\), then \(\sum_i w_i<S/6<1\). Such lists cannot
supply the strict CRT height threshold for profiles with all
three numerical heights one. A short-index relaxation still
needs to cover the fixed difference between the actual
numerical height \(h_*\) and \(S/6\). This restriction on
the charging rule neither constructs nor excludes an
arithmetic near-saturating family.

To prove the implication let \(d=1/6-u\), so that
\(\beta_i=2d\) and the target threshold is \(1/6+d\).
Every summand of \(E\) is at least \(d\). If \(t_i<u\),
it is at least \(d+t_i\). If all marks have \(t_i<u\),
their sum yields
\(\sum_i t_i<1/6-(\ell-1)d\le1/6\).
If some \(t_i\ge u\), its own term is at least
\(1/6+\min(t_i,2d)\). The strict budget requires \(d>0\)
and \(t_i<d\). Its term is then at least \(1/6+t_i\).
Any additional mark would already give at least
\(1/6+d\), so there can be only one mark, with \(t_i<d\).
A fixed exponent gain \(E\le1/3-u-\sigma\) improves both
cases to \(\sum_i t_i\le1/6-\sigma\). The numerical
bound follows from \(w_i\le S t_i\) and the nonnegative
overlap loss. All inequalities are
at the limiting exponents. The same obstruction to a
fixed power saving holds under sufficiently small
perturbations with the numerical gap held fixed.
[Fixed-letter marginal bound (FL1)](ABC_FIXED_LETTER_RADICAL_CAPS.md#a-pointwise-count-with-prescribed-prime-pools)

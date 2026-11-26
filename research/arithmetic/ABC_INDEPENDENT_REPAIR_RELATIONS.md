# Short relations on independently least-repaired letters

Least repair makes the radical of an original letter a single high power
times a cofactor of bounded original height. A small original prime packet
then supplies a modulus that can be much taller than this cofactor. This gives
a counting criterion for repaired ternary relations across the entire
balanced endpoint and a separate bound for repaired relations on forced
additive trades. The inverse repairs that restore the original identities
require large dictionaries of coefficient classes.

## Counting with a modulus from one repaired letter

Fix \(K>0\) before the small-budget limit, put \(n=X^L\) and
\(M=\lceil K/L\rceil\), and take \(X\) large at fixed \(L\).
For a positive integer \(c\le X\) put
\[
 R(c)=\operatorname{rad}(c),\qquad
 \mu_M(c)=\prod_{0<v_p(c)<M}p^{M-v_p(c)},\qquad
 W(c)=\prod_{v_p(c)>M}p^{v_p(c)-M}.
\]
Thus the least target repair has height \(H_M(c)=c\mu_M(c)=R(c)^M W(c)\).
Allow bounded integers \(1\le j\le J\), \(0\le k\le D\) and write
\[
 \Phi_{j,k}(c)=\mu_M(c)^j c^k
       =R(c)^{jM} V_{j,k}(c),\qquad
 V_{j,k}(c)=W(c)^j c^{k-j}.
 \tag{IR1}
\]
Its support is contained in that of \(c\), and \(H(V_{j,k}(c))\le X^{2J+D}\),
where rational height is the maximum of the absolute reduced numerator and
positive denominator. This bound is independent of \(M\). The case
\(j=1,k=0\) uses the least multiplier itself and \(j=k=1\) the least
repaired height. No arbitrary common multiple enters this definition.
[Least target repair](POWERFUL_TRIPLES_COMMON_DIVISOR_BOUND.md#a-reverse-transfer-by-repairing-only-the-target)

Let \(0<\beta<1\) be fixed, and suppose a letter supplies a squarefree
divisor \(d\) of its radical with
\(\beta/2\le\log_n d\le\beta\) and \(v_p(c)<M/2\) for every \(p\mid d\).
Set \(s=\lfloor M/2\rfloor\) and \(N=d^s\). Then \(N\mid \Phi_{j,k}(c)\)
for all the orders above. Recording \(d\) costs \(n^\beta\).
Fix also a height exponent \(B_*>0\). For \(K\) sufficiently large in
terms of \(\beta,B_*,J,D\), and for sufficiently small fixed \(L\), this
modulus has the following counting property. In it consider a unit congruence
\[
          V\prod_{i=1}^m r_i^{\pm j_i M}\equiv u\bmod N
 \tag{IR2}
\]
for \(m=1\) or \(2\), with pairwise coprime squarefree \(r_i\le n^{1/2}\)
and \(H(V)\le X^{B_*}\). The number of tuples in
(IR2) is at most
\[
     O\!\left(n^{(m^2+m)\beta+mB_*/K}X^\varepsilon\right).
 \tag{IR3}
\]
The conclusion holds uniformly in \(d\) and the target unit \(u\).

Indeed \(H(r_i)\le N^{(1/\beta)/s}\), while
\(\log_X N\ge sL\beta/2\ge K\beta/6\) at small fixed \(L\).
Thus \(H(V)\le N^{6B_*/(K\beta)}\). The inverse-product condition of (UW)
holds by taking \(K\) sufficiently large in terms of \(B_*,\beta\), and
the required fullness threshold holds for sufficiently small fixed \(L\).
For its exact-fibre hypothesis, fix a rational value of the left side
of (IR2). This value has polynomial height in \(X\) for fixed \(K\).
The primes of each \(r_i\) on its support can be chosen at height-subpower
cost, by the divisor bound. Off that support, pairwise coprimality makes
the corresponding part of \(r_i\) raised to \(j_iM\) divide the numerator
or denominator of \(V\). That part of \(r_i\) is at most \(n^{B_*/K}\).
After the \(r_i\) are fixed the exact value determines \(V\). Apply (UW),
with chart cost \(N^{(m^2+m)/s+\varepsilon}\), and absorb the arbitrarily
small fixed-height-exponent losses to obtain (IR3).
[Quantitative wide-cofactor estimate (UW)](ABC_QUANTITATIVE_DEGREE_WINDOWS.md#a-product-bound-for-the-admissible-wider-height)

## Ternary relations with an affordable coefficient dictionary

Fix \(J\ge1,D\ge0\) and a coefficient-height exponent \(B>0\). There are
positive \(\gamma,\delta\) and \(K_0\) with this property. For every fixed
\(K\ge K_0\) and all sufficiently small fixed \(L>0\), consider positive
primitive \(a+b=c\) in \(X/2<c\le X\) with
\(\operatorname{rad}(abc)\le n\). For each eligible row suppose some coefficient
triple from a prescribed dictionary of size at most \(n^\gamma\) makes
\[
 \sum_{T\in\{a,b,c\}} A_T\Phi_{j_T,k_T}(T)=0,
 \qquad A_T\in\mathbb Q^\times,\quad H(A_T)\le X^B,\quad
 1\le j_T\le J,\quad 0\le k_T\le D.
 \tag{IR4}
\]
The exponents can be chosen per row. Then the number of eligible rows is
\[
               O_{K,L,J,D,B,\varepsilon}
                   \bigl(X^{(1/3-\delta)L+\varepsilon}\bigr).
 \tag{IR5}
\]
The coefficient dictionary may depend on \(X\). No complete-value dictionary
bound or valuation grid is required.

First suppose the individual original radicals are at most \(n^{1/2}\),
and fix the record of coefficients and exponents. Let \(P\) be the prime
support of the six coefficient numerators and denominators. If one original
letter \(T\) supplies a packet \(d\) as above outside \(P\), reduce (IR4)
modulo \(N=d^s\). The other two original letters, denoted by \(U,V\),
are units by primitivity. The resulting congruence is
\[
 \frac{\Phi_{j_U,k_U}(U)}{\Phi_{j_V,k_V}(V)}
                             \equiv -A_V/A_U\pmod N.
\]
In (IR3) take the two narrow variables \(R(U),R(V)\), which are coprime,
and the wider variable \(V_{j_U,k_U}(U)/V_{j_V,k_V}(V)\), of height at most
\(X^{B_*}\) with \(B_*=2(2J+D)\). Each pair of original radical supports fixes
the two positive complete values at height-subpower cost. Those values give
the third through the original additive equation. Recording \(d\) and summing
(IR3) counts at most
\[
                   O(n^{7\beta+2B_*/K}X^\varepsilon)
 \tag{IR6}
\]
rows for fixed coefficients and orders, uniformly in the role supplying \(d\).
[Uniform supported-power enumeration](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)

To deduce (IR5), the two-letter bound discards a uniformly smaller exponent
if any individual normalized cap exceeds \(1/3+\kappa\), for a small fixed
\(\kappa>0\). Indeed the other two masses then sum to at most \(2/3-\kappa\)
up to arbitrarily small mesh errors. Fix a small positive overlap tolerance
\(\xi\), say \(1/100\). For each coefficient record, (EB5) discards the rows
with normalized radical overlap at least \(\xi\) on \(P\) at exponent
\(1/3-\xi/3+\eta\) in base \(n\), for any fixed small allowance \(\eta>0\).
The pool has polynomial product height, and the relative total budget at
height \(c\) may be padded to \((1+o_X(1))L\), or any fixed slightly larger
multiple of \(L\). Hence this discarded exponent remains below \(1/3\)
after summing the catalogue, by choosing \(\gamma,\eta\) small compared
with \(\xi\).
[Two-letter bound (GD3)](ABC_GROUPED_ROOT_DICTIONARY_ENTROPY.md#group-dictionary-count)
· [Fixed-pool price bound (EB5)](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)

If the rows remaining could saturate, anatomic refinements would give a
height-first balanced endpoint sequence. The total normalized radical
budget in a retained box is \(1+o(1)\) in units of \(L\): a fixed deficit
would contradict the limiting one-third bound, and the original cap gives
the upper bound. The box scale in (GE) is therefore \((1+o(1))L\).
The summed nonatomic mass on the scaled orders \(Lv_p\) is at least
\(1/2\) by (GE). Some letter consequently has nonatomic mass at least
\(1/6\). The bound on its first moment leaves more than \(1/10\) of this mass
below an absolute fixed upper scaled order \(R_0\), allowing small fixed
losses. Choose a compact set of continuity intervals with this much mass
and remove short neighborhoods of its finitely many atoms above the
threshold \(\beta/8\), taking all interval boundaries at continuity
points. By making the removals short enough, a fixed
nonatomic mass greater than \(2\xi\) remains. On the surviving bounded
degree bands no single prime has normalized mass greater than \(\beta/2\)
at sufficiently late boxes. Otherwise a subsequence of such primes at
convergent scaled degrees would give an atom of mass at least \(\beta/2\)
there. The construction uses exact prime valuations on the bounded bands.
The removal of the assigned pool \(P\) costs less than \(\xi\) on every
remaining row, leaving an available packet in \([\beta/2,\beta]\) by
successively adding whole primes. Take \(\beta\) a small fixed constant,
and \(K>2R_0\), to ensure every selected order is below \(M/2\).
[Endpoint diffusion (GE)](ABC_WIDE_COFACTOR_CONGRUENCES.md#necessary-conditions-at-normalized-exponent-13)
· [Actual bounded valuation bands and prime rounding](ABC_FIXED_FIBRE_VALUATION_WIDTH.md#discrete-rounding-on-the-full-three-letter-family)

Choose \(\beta,\gamma\) small and \(K_0\) large so the exponent in (IR6)
plus \(\gamma\) has a fixed margin below \(1/3\). This contradicts the
remaining count on the proposed sequence. The gap can be chosen uniformly over
fixed \(K\ge K_0\). If every candidate gap failed, choose one failing fixed
\(K_\ell\) for each gap tending to zero and then a sufficiently small fixed
\(L_\ell\to0\) at which that gap fails at arbitrarily large heights. The
smallness requirement at \(K_\ell\) for (UW) is imposed before choosing
\(L_\ell\). Both (IR6) and the diffusion contradiction have margins
independent of \(K_\ell\) above \(K_0\). The spectral subsequence, packet choices
and constant allowances are used with the height limit at each fixed
\((K_\ell,L_\ell)\) first.

### Projective power classes and the inverse-repair record

The hypothesis of (IR5) can also be imposed on the coefficient catalogue
after projective normalization and reduction modulo rational \(M\)-th powers.
Normalize a weight vector by its third entry, and choose one normalized
representative from the actual polynomial-height weights in each class.
Each ratio between a normalized entry and its representative has rational
height at most \(X^{4B}\) and is an \(M\)-th power. Its root has rational
height at most \(X^{4B/M}\), so the two entries have at most
\(n^{16B/K+o_X(1)}\) possible lifts together. Enlarge \(K_0\) and take
a smaller catalogue exponent to absorb this cost in (IR5).
Alternatively, a small joint radical cap on all coefficient numerators
and denominators supplies a small exact catalogue by supported enumeration.

The relation \(H_M(a)/\mu_M(a)+H_M(b)/\mu_M(b)=H_M(c)/\mu_M(c)\)
restores the original sum. Its projective coefficient record remains large
even if each ratio of inverse repair weights is taken modulo rational
\(M\)-th powers. In fact any one such record has at most
\[
                         O(n^{3/K}X^\varepsilon)
 \tag{IR7}
\]
primitive preimages, uniformly over all rows of height at most \(X\).

The record gives the valuation vector of the three coefficients at each prime
modulo \(M\) and modulo a common additive residue. If
\(0<v_p(T)<M\), the coefficient on role \(T\) has valuation
\(v_p(T)-M\), nonzero modulo \(M\), and the other two have valuation zero.
Thus the record identifies the unique different slot and recovers every
such prime on its original role. A prime of order at least \(M\) is absent
from the coefficients. The product of these high primes on each letter is
at most \(n^{1/K}\). Record the three missing supports and count the
supported original powers. This proves (IR7). In particular a family of
size \(n^{1/3-o(1)}\) has at least \(n^{1/3-3/K-o(1)}\) such canonical
coefficient classes at large fixed \(K\).

## Fixed shifts of least multipliers and height transforms

For both signs possible in the rational weights, consider
\[
                  A_1\Phi_{j_1,k_1}(c_1)
                  +A_2\Phi_{j_2,k_2}(c_2)=t\ne0.
 \tag{IR8}
\]
Here \(c_i\le X\), \(H(A_i)\le X^B\), and the two orders lie in the
fixed ranges of (IR1). Suppose \(R(c_1)\le n^{1/2}\) and \(c_2\)
provides a packet \(d\) as above coprime to \(c_1\) and disjoint from
the coefficient prime pool. Then the number of such pairs is
\[
                   O(n^{3\beta+(B_*+1)/K}X^\varepsilon).
 \tag{IR9}
\]
Reduce (IR8) modulo \(N=d^s\). The first term is a unit there,
so the target divided by \(A_1\) is a unit. Apply (IR3) with one
narrow variable \(R(c_1)\) and recover \(c_1\) from its support at
height-subpower cost. The exact shift fixes the value of the second
transform. Its prime support identifies every prime of \(c_2\)
of order below \(M\) up to the possible addition of high primes.
Choose the full set of high primes at cost \(n^{1/K}\). Its union
with the transform support is exactly the support of \(c_2\), which
then costs height-subpower to recover. The packet record costs
\(n^\beta\), proving (IR9).

There is also a fixed saving without the packet hypothesis.
For sufficiently small balance tolerance \(\kappa>0\), sufficiently large
fixed \(K\), and sufficiently small fixed \(L\), the count in (IR8)
under both caps \(R(c_i)\le n^{1/3+\kappa}\) is at most
\[
                        O(n^{1/6-\sigma}X^\varepsilon)
 \tag{IR10}
\]
for a fixed positive \(\sigma\). It is uniform in the prescribed shift and
coefficients. To prove it take a large fixed \(R_1>0\) and record the full
original prime powers of both letters at \(Lv_p(c_i)>R_1\).
Their joint radical costs at most \(2/R_1\) in base \(n\).
On the low primes, the transformed order is
\(j_iM+(k_i-j_i)v_p(c_i)\ge j_i(M-R_1/L)\).
Take \(q=\lfloor M-R_1/L\rfloor\), with \(K>R_1\). The remaining
variable parts of the two transforms are \(j_iq\)-full and have numerical
height at most \(X^{j_{\max}(K+L)(1/3+\kappa)+D}\), where
\(j_{\max}=\max(j_1,j_2)\). The frozen high powers enter the fixed
coefficients. Apply (SHF1) at \(q\) with its arbitrary fixed-coefficient
extension. To dispose of a common gcd after clearing denominators,
enumerate its divisor of the fixed nonzero shift and strip the free
parts supported on that divisor at subpower cost. The remaining two
free parts retain their fullness. Reconstruct the original letters from
their recovered transform supports and the frozen high parts.
[Strict fixed-shift saving with arbitrary coefficients (SHF1)](ABC_FIXED_SHIFT_STRICT_HALF_SAVING.md#why-equality-of-the-raw-and-full-height-budgets-is-impossible)

Since \(1\le j_i\le J\), the finitely many applicable half-coefficient
gaps have a positive minimum \(\delta_*\). The normalized bound is at most
\[
 \frac{2}{R_1}
  + \left(\frac{1}{2j_{\max}}-\delta_*\right)
        \frac{j_{\max}(K+L)(1/3+\kappa)+D}{K-R_1-O(L)}
\]
up to arbitrarily small fixed losses. Take \(R_1\) large, then \(K\)
large compared with \(R_1\), and then \(\kappa,L\) sufficiently small.
This has a fixed margin below \(1/6\). Browning–Verzobio's uniform
generalized Fermat estimate treats one fixed exponent per root and their
Theorem 4.1 concerns ternary powerful sums. Here the nonzero shift, the
varying original valuation profile, and the row-dependent repair are
handled by freezing the high valuations and using (SHF1) on the remaining
full factors. [Browning–Verzobio, Theorems 1.2 and 4.1](https://arxiv.org/html/2608.24512v1)

## Application to the forced additive trades

Use the forced eight-edge trade (BT2) and write its four original
targets as \(c_{ij}=a+b+h+it+js\), \(i,j\in\{0,1\}\).
They satisfy the original alternating zero-sum.
Suppose every complete dictionary has size at most \(n^{1/6+\kappa}\)
and every original role radical is at most \(n^{1/3+\kappa}\), with
small fixed \(\kappa>0\). For each fixed set of four nonzero rational
weights of polynomial height and bounded transform orders, the number
of faithful ordered trades that satisfy
\[
                \sum_{i,j\in\{0,1\}}
                         A_{ij}\Phi_{j_{ij},k_{ij}}(c_{ij})=0
 \tag{IR11}
\]
is \(O(n^{5/6-\sigma+4\kappa}X^\varepsilon)\). Choose two signed
terms with the same sign, fix the other two target values, and apply
(IR10) to the chosen pair. The shift in any solution is nonzero because
the chosen terms have the same sign. This costs
\(n^{2(1/6+\kappa)+1/6-\sigma}\) up to height-subpower factors.
Four target values determine \(t,s,a+b+h\). The two actual summand
labels \(a,a+h\) determine the remaining \(a,b,h\), giving at most
\(n^{2(1/6+\kappa)}\) completions. This proves the stated count.
A dictionary of \(n^\gamma\) coefficient vectors adds \(\gamma\) to
the upper exponent. A nonzero output with outside-radical cap \(n^b\)
can also be prescribed by enumeration, at additional cost \(b\) in the
exponent. The output denominator is supported on the coefficient primes,
and all values have polynomial height. The selected same-sign pair
still has a nonzero shift.
[Forced eight-edge trades (BT1)–(BT2)](ABC_BALANCED_LATIN_TRADE_OBSTRUCTION.md#eight-edge-supersaturation-for-actual-additive-values)

In particular the actual least-multiplier zero-sum almost never survives
a faithful trade in a balanced near-full minimum-dictionary graph,
since (BT1) gives \(n^{5/6-o(1)}\) such trades. To see the connection
to repairing the rows, apply the target least multiplier to each of
the eight edges and sum repaired first summands with the signs
of (BT2). In cell \(i,j\), the two repaired first summands differ by
\(h\mu_M(c_{ij})\). The alternating repaired defect is exactly
\[
                       h\sum_{i,j}(-1)^{i+j}\mu_M(c_{ij}).
\]
The second summands give its negative and the targets cancel cellwise.
Thus (IR11) for \(j_{ij}=1,\ k_{ij}=0\) tests survival of the actual
trade with the canonical repair, including the numerical multiplier
on each edge.

On trades where the second queried target has an affordable packet
private to it among the queried targets and outside the weight support,
(IR9) sharpens the upper exponent in (IR11) to
\(2/3+4\kappa+3\beta+(B_*+1)/K\), up to height-subpower losses.
For a balanced minimum nonatomic profile the packet hypothesis holds
on \(n^{5/6-o(1)}\) trades for each fixed weight pool. The private
support estimate (FD2) discards any pairwise radical overlap of
fixed positive normalized mass. For a fixed weight pool the same
proof fixes an incident companion label, extends to a five-label
basis, and applies (FL1) with that prescribed pool, at the same
exponent loss for a large overlap. Then select the packet from
a bounded nonatomic band on the query letter by whole-prime rounding.
Thus in this setting the upper exponent on the surviving private
trades can be made arbitrarily close to \(2/3\).
[Private prime supports on forced trades (FD2)](ABC_PROJECTIVE_RECIPROCAL_GRAPH_OBSTRUCTION.md#the-denominator-of-a-forced-trade-in-a-balanced-capped-graph)
· [Fixed-letter bound (FL1)](ABC_FIXED_LETTER_RADICAL_CAPS.md#a-pointwise-count-with-prescribed-prime-pools)

### The coefficient cost of restoring the target relation

On a subfamily of faithful trades with pairwise target overlap at most
\(n^\xi\), record the vector of the four inverse target multipliers up to
rational scaling and coordinatewise rational \(M\)-th powers.
Fix one value of this projective coefficient record. Enumerate the union
of all nonprivate target primes at cost \(n^{6\xi}\), since it is supported
on the union of the six pairwise gcd radicals. At every other prime in
the coefficient record there is a unique different slot modulo \(M\),
as in (IR7). It gives precisely the private primes at original order
below \(M\), with their roles among the four targets. Choose the missing
private high supports at cost at most \(n^{4/K}\), and then the complete
target values at height-subpower cost. Thus a coefficient record has at
most \(n^{6\xi+4/K+o_X(1)}\) target quadruples. Signs on the four weights
do not affect this reconstruction.

At the balanced minimum-dictionary exponent, (BT1) and the pair overlap
bound provide \(n^{5/6-o(1)}\) trades after taking the count and cap
tolerances sufficiently small relative to \(\xi\). Four target values
leave at most \(n^{1/3+2\kappa}\) trade completions. Consequently the
coefficient class dictionary of the canonical correction on these
trades has exponent at least \(1/2-6\xi-4/K-o(1)\).
This approaches the three-target-parameter exponent \(1/2\) as
\(\xi\) decreases and the fixed repair parameter \(K\) increases.
The record uses the actual target repairs on the primitive quotient.

The bounds require an actual repaired relation or a cap on its output.
Ordinary trade supersaturation gives an original relation. Its canonical
repair can produce a large nonzero defect or require an expensive
coefficient record. The unrestricted exceptional-set exponent remains
open.

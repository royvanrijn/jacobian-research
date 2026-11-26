# Atomic profiles and the price of two compressed records

A bounded-root representation with a small integral correction
forces an atomic valuation profile on the primes outside its fixed
coefficient. This follows from the actual root heights, including
when the roots overlap on the same letter. In particular an
inexpensive compressed record on two fully diffuse letters cannot
be obtained from saturation by discounting the price of their
recorded supports.

## CRP1. A compressed remainder is atomic

Put \(n=X^L\) and take the height limit at each fixed \(L\) before
\(L\to0\). Fix \(A,B,H>0\) and a nonnegative integer \(m\).
Consider positive integers \(T\le X^H\) with
\(\operatorname{rad}(T)\le n^A\), and representations
\[
           T=F u_0\prod_{i=1}^m x_i^{e_i},\qquad
           F,u_0,x_i\in\mathbb Z_{>0},\quad e_i\ge1 .
\]
Assume \(x_i\le n^B\) and that \(\log_X u_0\) tends to zero
in the iterated limits. The root and coefficient supports may
overlap. The orders need not be common to all the integers in
a family. Write
\[
 \nu_{T,L}=\sum_{p\mid T}\log_n p\,\delta_{Lv_p(T)},\qquad
 \mu_{T,F,L}=\sum_{\substack{p\mid T\\p\nmid F}}
                               \log_n p\,\delta_{Lv_p(T)} .
\]
For every subsequence on which \(\nu_{T,L}\to\nu\) and
\(\mu_{T,F,L}\to\mu\) weakly, the measure \(\mu\) is purely
atomic and satisfies \(\mu\le\nu\). In particular the measure
of the primes in the record has nonatomic part \(\nu_{\rm na}\)
in this limit, and its total mass is at least
\(\nu_{\rm na}([0,\infty))\). This applies in particular to
the smaller integral cofactors in (DTC1).
[Two compressed letters (DTC1)](ABC_TWO_COMPRESSED_LETTER_RELEASE.md#two-compressed-terms)

For each prime outside \(F\), put \(\beta_p=\log_n p\),
\(k_{p,i}=v_p(x_i)\), and \(s_p=L v_p(u_0)\). Then
\[
 Lv_p(T)=s_p+\sum_i (Le_i) k_{p,i},\qquad
 \sum_p\beta_p\sum_i k_{p,i}\le mB,\qquad
                         \sum_p\beta_p s_p=o(1).
\]
Thus the primes with \(\sum_i k_{p,i}>K\) have mass at most
\(mB/K\). For a fixed \(\theta>0\), the primes with
\(s_p>\theta\) have vanishing mass. On the other primes the
scaled degree lies within \(\theta\) of a sum of the
\(Le_i\) with nonnegative integer coefficients of total at
most \(K\).

Pass to a further subsequence on which each \(Le_i\) has a limit
in \([0,\infty]\), padding the list by unit roots if needed.
A sum using a diverging order leaves every compact degree range.
The original height bound gives tightness
\(\nu_{T,L}([R,\infty))\le H/R\), including in the presence
of escaping numerical weight. The finite limiting orders
generate a locally finite set of nonnegative integer sums.
Indeed zero orders can be omitted in forming distinct centres,
and any positive limiting order bounds its own coefficient on a
compact interval. Let \(E\) be that set, including zero.

For any continuous function between zero and one with compact
support disjoint from \(E\), take \(\theta\) smaller than its
separation from the relevant centres. For sufficiently late
budgets all the bounded-multiplicity primes with \(s_p\le\theta\)
then contribute zero to its integral. The limiting integral is
at most \(mB/K\), and \(K\) is arbitrary. This proves that
\(\mu\) is supported on \(E\). The inequality \(\mu\le\nu\)
follows directly from the same inequality before taking limits.
Since \(\nu-\mu\) subtracts a purely atomic submeasure, it retains
the full nonatomic part of \(\nu\). The proof imposes no
restriction on numerical weight which escapes to arbitrarily
large degree.

## CRP2. An actual joint-price obstruction on a diffuse pair

Work with actual primitive triples in \(X/2<C\le X\) and bounded
individual radical caps in base \(n\). Use height boxes in which
the original exact valuation measures on the three roles tend to
common limits \(\nu_T\); here common means that the limits hold
for every choice of retained row in the iterated box limits.
Write \(r_T=\nu_T([0,\infty))\) and
\(c_T=(\nu_T)_{\rm na}([0,\infty))\). Suppose a subfamily of
limiting count exponent \(d\) admits joint records \(F_U\mid U\),
\(F_V\mid V\) on two prescribed distinct roles, and their
two remainders satisfy (CRP1) with fixed cap bounds. The joint
record dictionary exponent \(f\) satisfies
\[
         f\ \ge\ d-\frac12
          \min\{a_U+a_V,\ a_U+a_W,\ a_V+a_W\},
 \qquad
         a_U=r_U-c_U,\quad a_V=r_V-c_V,\quad a_W=r_W .        \tag{CRP2}
\]
This is also a lower bound on the joint record entropy in units
of \(\log n\), when one admissible record is selected per row.
The entropy and dictionary assertions use vanishing normalized
allowances in the height-first convention.

Indeed (CRP1) gives these outside-support caps up to arbitrarily
small losses uniformly in the box. Otherwise one could select
violating rows and take weak subsequential limits, contrary to
that assertion. Fix the two record values. Fix also the full
parts of \(U,V\) on their respective recorded prime supports.
For given supports, these additional values cost \(X^\varepsilon\)
with arbitrarily small \(\varepsilon>0\), by the fixed-support
bound. This handles partial recorded powers and roots that share
their support with the recorded factors. The remaining radical
caps on the original three roles are \(a_U,a_V,a_W\) up to
the indicated allowances. Apply (MCB2) to the original primitive
equation, using the full parts as fixed coefficients. Its
two-cheapest bound applies with the original maximum \(C\asymp X\).
The exceptional coefficient and binning losses are taken within
an arbitrarily small normalized allowance at the fixed budget.
Every actual joint record therefore has fibre exponent at most
half the displayed minimum. Counting the fibres gives the
dictionary assertion. For the entropy assertion subtract this
uniform bound from \(H(\mathrm{row})/\log n=d+o(1)\) using the
conditional entropy identity for the selected record.
[Arbitrary-coefficient two-cheapest bound (MCB2)](ABC_MULTIPLE_SHORT_COFACTOR_BOUNDARY.md#bounded-group-dictionaries-with-arbitrary-fixed-coefficients)
· [Uniform counting on prescribed supports](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)

There is a stronger price bound using the atomicity of both residual
measures. In the fixed joint record, the number of possible
complete \(U\)-values has exponent at most \(a_U\) by the
one-variable radical count after the supported powers are fixed.
For each such \(U\), prescribe the \(F_V\)-support on \(V\) in
the original fixed-\(U\) fibre. Its outside measure is purely
atomic by (CRP1), so (FT1) gives a subpower completion count.
Uniformity follows by taking a violating fixed-fibre subsequence
if necessary. The fixed complete value and the pool can vary
with the height and budget in (FT1). Interchanging the two roles
gives a per-record upper exponent \(\min(a_U,a_V)\). In particular
\[
                         f\ge d-\min(a_U,a_V).              \tag{CRP2a}
\]
The two-cheapest bound in (CRP2) remains available with the
third cap. The argument for (CRP2a) uses just the atomicity of
the two outside measures and applies at an arbitrary complete
multiplicity.
[Fixed fibres outside a prescribed prime pool (FT1)](ABC_FIXED_FIBRE_TRIANGLE_RADICALS.md#nonatomic-mass-outside-prescribed-prime-pools)

In a balanced saturation sequence (CRP2a) gives
\(f\ge\max(c_U,c_V)\). If either selected letter is wholly
nonatomic, its residual radical cap is zero in the limiting
exponents, and \(f\ge d=1/3\). The other letter may have a
mixed profile. This gives an actual obstruction to a catalogue
discount on such a pair, uniformly over which records are
selected. It concerns necessary prices in a putative
saturating sequence. Arithmetic realization of that sequence
remains open.
[Joint record extraction criterion (DTC2)](ABC_TWO_COMPRESSED_LETTER_RELEASE.md#a-joint-two-divisor-record-across-the-balanced-range)

## CRP3. The converse at vanishing numerical width

Suppose a fixed-support record leaves purely atomic weak
outside-support measures, for every subsequential limit of
the retained rows in a common profile box. At an arbitrarily
small extra normalized record cost, these remainders admit a
bounded number of powered roots of height \(n^{O(1)}\) and
integral cofactors of height \(X^{o(1)}\).

First record the full powers on the original support at inner
subpower cost. Fix \(\delta>0\), and choose finitely many positive
atoms of the original limit \(\nu\) leaving atomic mass less
than \(\delta\) outside that list and the zero atom. Every
outside-support limit is dominated by \(\nu\), so its atomic
mass has the same finite-list bound. Take disjoint slowly
shrinking windows around the selected atoms and a shrinking
window at zero. Their widths can tend to zero slowly enough
that the residual mass outside them is at most \(\delta+o(1)\)
uniformly in the box. To see this, first fix each positive
window width with continuity boundaries for \(\nu\), apply
subsequential compactness and domination to bound the mass
outside, and then decrease the widths diagonally after the
height limits. Fix all original powers outside these windows.
Their radical cap gives cost at most \(\delta+o(1)\) in base
\(n\), with an inner subpower cost for the full powers. This
includes any high tail with escaping numerical weight.

In each positive window group the squarefree product of its
primes at the minimum original order present. All remaining
orders in the window differ from it by at most the window
width divided by \(L\), up to rounding. The grouped root
has height at most \(n^A\). The positive exponent differences
and the complete zero-window factor make an integral cofactor
of numerical height exponent tending to zero. This proves the
claim. If \(\nu_{\rm na}\) has positive mass, choose a compact
positive continuity band on which it has positive mass.
All that nonatomic mass is covered by the fixed support in
the limit. The full powers of the fixed support in this band
then give a fixed full divisor of positive numerical weight
and fullness comparable to \(1/L\). It may be selected from
the fixed record after recording the full powers.

The construction only guarantees a cofactor of height
\(X^{o(1)}\). A width of order \(L^{1-\xi}\) in scaled
degrees would instead give height \(n^{O(L^{-\xi})}\), which
requires discrete information about the convergence.
[Thin-window integral cofactors (TW2)](ABC_THIN_WINDOW_ARBITRARY_SHIFT_RELEASE.md#fixed-shifts-with-a-wider-integral-correction)

## CRP4. Weak atomicity can miss every stretched width

Take a sequence of integers \(q\to\infty\), put \(L=1/q\),
and fix each budget before its height limit. Set
\(J=2K+1\) with \(K=\lfloor q/\log q\rfloor\), and use
the \(J\) orders \(e_j=3q+j\) for \(-K\le j\le K\).
Give each order its own prime pool in
\([(1-\tau_q)n^{1/(3J)},(1-\tau_q/2)n^{1/(3J)}]\), with
\(\tau_q>0\) fixed and \(\tau_q\sum_j e_j<1/10\).
The pools are disjoint and of equal size up to constants,
by the prime number theorem at fixed \(q\). The integers
\(T=\prod_j p_j^{e_j}\) from these choices lie in \((X/2,X)\),
have radical at most \(n^{1/3}\), and number
\(n^{1/3-o_X(1)}\). The height calculation uses
\(\sum_j e_j=3qJ\). Their common limiting profile is the
single atom \(\frac13\delta_3\), since the scaled spread
\(w_L=KL\) tends to zero.

Fix \(m,B,C>0\) and \(0<\xi<1\). For any representation in
(CRP1) with \(u_0\le n^{C L^{-\xi}}\), the fixed coefficient
in this dictionary must have
\(\operatorname{rad}(F)\ge n^{1/3-o(1)}\).
Indeed choose a fixed multiplicity cutoff \(K_0\) in the
proof of (CRP1). Here
\(\sum_p\beta_p s_p\le C L^{1-\xi}\).
Use a window width
\(\theta_L=(L^{1-\xi}w_L)^{1/2}\). The exceptional cofactor
mass is at most \(C L^{1-\xi}/\theta_L=o(1)\).
Each of the boundedly many remaining windows contains at
most \(O(1+q\theta_L)\) exact orders. Each exact order carries
mass \(1/(3J)+o_X(1)\), so all these windows together have
mass \(O_{K_0,m}((\theta_L+L)/w_L)+o_X(1)=o(1)\).
The remaining outside-support mass is at most \(mB/K_0+o(1)\).
Taking \(K_0\) arbitrarily large proves the assertion.
Moreover a fixed coefficient covering that mass leaves only a
subpower choice of pool labels in this dictionary, in the
iterated exponent sense, since its support identifies the
primes in all but a vanishing fraction of the slots.

With the less restrictive numerical correction condition,
take \(F=1\), use one root \(\prod_j p_j\) at order
\(3q-K\), and put the positive differences into the
cofactor. Its normalized numerical height is
\(O(1/\log q)\) in base \(X\), which tends to zero.
Thus the distinction between weak atomicity and a prescribed
shrinking width persists for actual bounded-height integers
with the original radical cap.
[Exact valuation slots and arithmetic projection dictionaries (RO1)](ABC_RESIDUE_OCCURRENCE_ROOT_PROJECTION_OBSTRUCTION.md#an-arithmetic-dictionary-from-finite-field-evaluations)

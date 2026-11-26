# Prime-power packing in a minimum complete dictionary

A minimum complete-dictionary endpoint has no mass at any finite
scaled valuation degree greater than \(12\), even when its total
nonatomic mass is larger than \(1/2\). Its nonatomic support is
contained in \([0,6]\). Between degrees \(6\) and
\(12\), a positive-mass compact band must contain a prime whose
individual radical mass exceeds an explicit threshold. One can
test several copies of such a prime from different rows together
with smaller prime powers from the same rows.

## Setup and the count of common columns

Put \(c_0=1/6\) and \(n=X^L\), with \(L>0\) fixed before the height
limit. Work with positive primitive triples \(a+b=c\),
\(X/2<c\le X\), and write the relation in signed coordinates
when changing roles. Assume each of the three complete-letter
dictionaries has size at most \(n^{c_0+\kappa}\) and each individual
radical is at most \(n^{2c_0+\kappa}\). The graph of two prescribed
complete letters is simple, since they determine the third letter.

For any fixed positive integer \(m\), a subfamily with
\(E\ge n^{2c_0-\tau}\) and \(\tau+\kappa<c_0\) has \(m\) distinct
complete rows \(a_1,\ldots,a_m\) in any one role with
\[
        \gg_m n^{c_0-m\tau-(2m-1)\kappa}                         \tag{KP1}
\]
common columns in a second role. Indeed if the numbers of distinct
rows and columns are \(A\) and \(B\), respectively, the mean column
degree is at least \(E/B\to\infty\). The number of ordered stars
with distinct rows, summed over columns, is
\(\sum_t (d(t))_m\gg_m E^m/B^{m-1}\), by restricting to
columns of degree at least \(2m\) and applying convexity. Divide
by at most \(A^m\) possible ordered row tuples and use the part
caps. This is the loss in density from selecting the rows. All
target constraints imposed on the subfamily hold on every edge
of the resulting star.

Let \(T_i(t)\) be the signed value in the remaining role at
column \(t\). It has the form \(T_i(t)=\pm t+s_i\), with the
same slope in every row and distinct fixed offsets. For a given
column the targets are positive in absolute value. Set
\[
               P=\prod_{i<j}|s_i-s_j|,
\]
with \(P=1\) if \(m=1\). This is a nonzero integer of height at
most \(X^{O_m(1)}\). A prime not dividing \(P\) cannot divide two
different \(T_i(t)\). For any fixed \(\xi>0\), omit columns for
which, on some row, the normalized radical mass of \(|T_i(t)|\)
supported on \(P\) is at least \(\xi\). With that complete row
fixed, (FL1) applied to the two remaining letters gives at most
\[
       O\bigl(n^{c_0+\kappa/2-\xi/2+\gamma}\bigr)                 \tag{KP2}
\]
such columns per row, with any prescribed \(\gamma>0\) and all
sufficiently small fixed \(L\). The estimate is uniform in the
row values and in \(P\). If the putative outside cap is negative,
there are no such columns. Inner height allowances are included
in \(\gamma\). In particular, for fixed \(\xi,m\), first fix
\(\gamma\) small and then \(\tau,\kappa\) small enough that
\(\xi/2>m\tau+(2m-1/2)\kappa+\gamma\). The omitted columns are
power fewer than those in (KP1).
[Fixed-letter bound with a prescribed prime pool](ABC_FIXED_LETTER_RADICAL_CAPS.md#a-pointwise-count-with-prescribed-prime-pools)
· [Earlier fixed-row argument using outside labels](ABC_DIFFUSE_FACTOR_GRAPH.md#combining-congruences-from-several-fixed-rows)

## Several large primes with smaller-prime contributions

For a full prime power \(p^e\parallel |T_i(t)|\), write
\[
             b(p)=\frac{\log p}{L\log X},\qquad r(p)=Le .
\]
Its numerical weight is \(r(p)b(p)=\log_X(p^e)\).
Suppose a prescribed sector on the target role has the following
properties on every edge of its subfamily.

- The target has a distinguished prime with
  \(b(p)\in[t_-,t_+]\) and \(r(p)\in[s_-,s_+]\),
  where \(0<t_-\le t_+\) and \(0<s_-\le s_+<\infty\).
  A deterministic choice suffices if there are several such primes.
- There are finitely many filler types indexed by \(j\). On each
  applicable row, type \(j\) consists of primes disjoint from the
  distinguished prime and from the other types. Each has
  \(r(p)\in[u_j,U_j]\), with \(0<u_j\le U_j<\infty\), and
  \(b(p)\le b_j\). Their combined radical mass on that row is at
  least \(\rho_j>0\).

Take \(k\) distinguished copies from \(k\) distinct rows in a
star of \(m\ge k\) rows. Require each type \(j\) on \(m_j\)
of those rows, for example on the same \(k\) rows when \(m=k\).
The test below only uses the stated properties on the rows where
their corresponding records are taken. In particular, if the
sector imposes every property on every edge, designate the
applicable rows after selecting the star. For each used type
choose \(q_j>0\) with \(q_j<m_j\rho_j\). Types with \(q_j=0\)
are omitted from the two sums. If
\[
 \boxed{\quad
  k t_++\sum_j(q_j+b_j)<c_0,\qquad
  k s_-t_-+\sum_j u_jq_j>1,
 \quad}                                                        \tag{KP3}
\]
then the sector has a strict uniform count saving below
\(1/3\) in base \(n\), for sufficiently small fixed budgets
and sufficiently tight complete-dictionary and individual radical
caps. The assertion includes the test with no filler, and also
the estimate using only fillers, with \(k=0\). If the sector supplies several prescribed bins of distinguished
primes on every edge, use disjoint lists of rows for the different
bins and sum their respective cost and numerical weight bounds
in (KP3).

For the proof choose \(\xi<t_-\) if \(k>0\), and for each used type
choose it small enough that \(m_j(\rho_j-\xi)>q_j\). Select the
common-column star using (KP1). After the omissions in (KP2),
none of the chosen distinguished primes divides \(P\). For each
type \(j\), remove all primes of \(P\) from the available packets.
On its chosen rows the remaining raw mass is at least
\(m_j(\rho_j-\xi)>q_j\). Accumulate entire prime powers greedily
until the raw cost reaches \(q_j\). The resulting cost is at most
\(q_j+b_j\) and the numerical weight is at least \(u_jq_j\).
Selections for different types have disjoint supports. The
chosen full factors on different rows have disjoint supports
since every chosen prime is outside \(P\).

Let \(M_i\) be the product of all recorded full powers on
row \(i\). Put \(M=\prod_i M_i\) and
\(p_* = k t_+ +\sum_j(q_j+b_j)\). Then
\(\operatorname{rad}(M)\le n^{p_*}\), \(M\le X^m\), and
\(M>X\) by (KP3). Enumerate \(M\) by first specifying its
radical and then the powers on that fixed support. This gives
\(O(n^{p_*+\gamma'})\) possibilities for any \(\gamma'>0\)
at fixed \(L\). At most \(X^{o_X(1)}\) choices of the ordered
factorization \((M_1,\ldots,M_m)\) occur for each \(M\), since
\(m\) is fixed. For a fixed factorization,
\(T_i(t)\equiv0\pmod{M_i}\) for all \(i\) puts \(t\) in one
residue class modulo \(M\) by the Chinese remainder theorem.
The complete column value lies in an interval of length at
most \(X\). Thus at most one column occurs for each record.

The chosen \(\xi\), \(m\) and the strict inequalities in (KP3)
are fixed before the balance and graph tolerances. Take
\(\gamma'\), \(\tau\) and \(\kappa\) small enough that the sum of
the exponent losses in (KP1) and the inner enumeration allowances
is smaller than \(c_0-p_*\). Apply (KP2) with the same fixed choices.
The count of surviving columns now contradicts (KP1). Consequently
for some \(\tau_*,\kappa_*,L_*>0\) depending only on the prescribed
parameters of the test, the sector has fewer than
\(n^{2c_0-\tau_*}\) edges for sufficiently large \(X\) at every
fixed \(0<L<L_*\) and \(0\le\kappa\le\kappa_*\).
All selected blocks may depend on the complete value. The argument
enumerates full integers, so this adaptivity has no extra cost
at a fixed number of rows.
[Enumeration on fixed prime supports](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)

For example, take a distinguished bin near \(t=3/50\),
\(s=15/2\), and one filler type in a compact band starting at
degree \(4\), with mass at least \(1/50\) per row and individual
mass at most \(1/200\). For \(k=m=m_1=2\) and \(q_1=3/100\),
the total raw ceiling is \(31/200<1/6\) and the lower numerical
weight tends to \(51/50>1\) as the distinguished bin shrinks.
The two distinguished copies by themselves contribute weight
\(9/10\), and three copies cost more than \(1/6\).

## A shorter interval for the column parameter

Suppose the common column role lies in one prescribed residue
class modulo an integer \(G\le X^K\), where \(K\) is fixed and
\(G\ge X^\beta\) for some fixed \(0\le\beta<1\). The choice of
\(G\) and of the class can vary with height and budget. Write
\(t=t_0+G y\). The integer \(y\) ranges in an interval of
length at most \(X/G\le X^{1-\beta}\). Include \(G\) in the
pool \(P\) in (KP2). Its height remains polynomial. Each
retained \(M_i\) is then coprime to \(G\), so the congruence
\[
                   T_i(t_0+Gy)\equiv0\pmod{M_i}
\]
determines one class for \(y\) modulo \(M_i\). All the earlier
estimates hold, with the second strict inequality of (KP3)
replaced by
\[
                    k s_-t_-+\sum_j u_jq_j>1-\beta .             \tag{KP4}
\]
The lower bound from (KP1) still counts the actual complete
column values; distinct such values give distinct \(y\).
Only the final reconstruction interval is shortened. If all
three original roles lie in prescribed classes modulo \(G\),
the same width is available in every orientation.

## The compact cutoff and the large-prime alternative

Fix a compact band \([s,S]\), and suppose a target has a
selected subset of primes of normalized radical mass at least
\(\rho>0\) in the band. If each selected prime has individual
mass at most \(b\), then this sector has a strict count saving
whenever
\[
                        (1-\beta)/s+b<c_0.                       \tag{KP5}
\]
Here \(\beta=0\) is always allowed, and positive \(\beta\)
uses the shorter interval above in the chosen orientation.
Use only one filler type, choose \(q>(1-\beta)/s\) with
\(q+b<c_0\), and take \(m\rho>q\). For
\(s>6(2-\beta)\), the bound \(b=1/s\) holds automatically
for every selected prime, since \(p^e\le X\). This proves a
uniform count saving when a fixed compact band above
\(6(2-\beta)\) has any prescribed positive mass.

Consequently, in a height-first endpoint of exponent \(1/3\)
with all the complete dictionaries of exponent at most \(1/6\),
no weak scaled valuation measure can charge
\((6(2-\beta),\infty)\). A positive weak mass there would give
a compact continuity subinterval carrying positive mass in the
original family, so the preceding uniform bound would apply.
The nonatomic support is contained in \([0,6(1-\beta)]\).
Indeed a compact band with positive nonatomic mass strictly
above this cutoff has \(s>6(1-\beta)\). Choose \(b>0\) with
\((1-\beta)/s+b<c_0\). From the band remove small neighborhoods
of the finitely many atoms of mass at least \(b/2\), preserving
a positive nonatomic mass. The remaining compact set can be
taken to be a finite union of continuity intervals. At late
budgets, a prime there of mass greater than \(b\) would, along
a convergent subsequence of its scaled degrees, give an atom
of mass at least \(b\) in this set. This is impossible.
The surviving positive-mass set therefore supplies the filler
in (KP5). These statements concern finite degrees. They leave
open escaped first moment at orders tending to infinity with
the outer budget limit.

For \(6(1-\beta)<s\le6(2-\beta)\), any such fixed compact band
with positive raw mass must, outside a power-smaller subfamily, include a prime
with individual mass greater than
\(c_0-(1-\beta)/s-\epsilon\), for every fixed
\(0<\epsilon<c_0-(1-\beta)/s\). This is a statement about actual
primes in the band, including at atoms of the limiting measure.
For a prescribed bin with limiting mass \(t>0\) and degree \(s>0\),
arbitrarily narrow bin tolerances in (KP3) and (KP4) give,
for each fixed
integer \(k\ge1\),
\[
                         kt<c_0\quad\Longrightarrow\quad kst\le1-\beta. \tag{KP6}
\]
More generally, if the same rows provide filler types with
arbitrarily small individual masses and with the stated positive
mass in compact bands, a strict inequality
\(kt+\sum_jq_j<c_0\), together with \(q_j<m_j\rho_j\),
forces \(kst+\sum_j u_jq_j\le1-\beta\). Finite overshoots are
charged explicitly in (KP3). Both deductions apply to a bin that itself
retains the full count exponent.

### A sharper cutoff in a progression

In a compact band \([s,S]\) with \(s>6\), every prime power
has individual raw mass at most \(1/s<c_0\), so it can be
used as a one-row record. For any fixed \(\eta>0\), consider
the sector with a prime power in that band of numerical
weight above \(1-\beta+\eta\). Its raw mass is at least
\((1-\beta+\eta)/S\). Choose the omitted pool mass \(\xi\)
smaller than this bound. Use one fixed row with many
columns and strip the primes of \(G\) by (KP2). The selected
prime power survives. Its raw enumeration costs at most
\(1/s+o(1)\), and it exceeds the interval length for \(y\).
With the record fixed, its unit congruence for \(y\) yields
at most one column. This sector therefore has a strict
saving.

On the remaining sector, all individual prime powers of the
band have raw mass at most \((1-\beta+\eta)/s\). For
\[
                     s>\max\{6,12(1-\beta)\},
\]
take \(\eta>0\) sufficiently small that
\((2(1-\beta)+\eta)/s<c_0\). The bound (KP5) now rules out
positive mass in the band at full exponent. Thus the
finite-degree support cutoff improves to
\[
                  \operatorname{supp}(\nu_W)
                     \subseteq[0,\max\{6,12(1-\beta)\}].          \tag{KP7}
\]
Here \(W\) is the target in the chosen orientation, and
\(\beta\) is the fixed lower bound on the numerical exponent
of the common column modulus. The one-row exclusion uses a
finite upper endpoint \(S\). It gives no bound on the first
moment escaping to scaled orders tending to infinity.

The original fixed-row CRT argument controls nonatomic support
using outside-label boxes at the minimum nonatomic mass. The
fixed-letter pool discount (KP2) only needs the complete
dictionary caps, and buying actual prime powers tests finite
atoms as well. Thus the compact exclusion applies, for example,
to a critical primitive triangle output once its complete
dictionaries are shown to be minimum, even in the absence of a
minimum nonatomic mass statement for that output.
[Earlier nonatomic support bound (MS11)](ABC_DIFFUSE_FACTOR_GRAPH.md#consequence-for-all-three-nonatomic-supports)
· [Critical primitive triangle transfer](ABC_FIXED_FIBRE_TRIANGLE_RADICALS.md#an-incidence-transfer-at-the-critical-triangle-cap)

## A common modulus from escaping valuation weight

In a height-first minimum-dictionary endpoint, refine the exact
valuation layers and numerical heights so that the weak degree
measures \(\nu_W\) and the limiting numerical heights \(h_W\) exist
on each of the three roles. Write
\[
          e_W=h_W-\int r\,d\nu_W(r),\qquad
          e_{\rm tot}=\sum_W e_W.
\]
There is a refinement with the same count exponent in which all
three complete roles lie in prescribed residue classes modulo
one fixed integer \(G\), with \(\log_XG\) tending to \(e_{\rm tot}\).
The modulus and residue classes are fixed within each resulting
family. They may vary with the height and budget.

First fix a large common continuity cutoff \(R\). For each complete
letter \(W\), take the whole tail
\[
            F_W=\prod_{\substack{p^e\parallel W\\ Le\ge R}}p^e.
\]
It has \(\operatorname{rad}(F_W)\le n^{1/R}\). Indeed each
retained prime contributes at least \(R\) times its raw mass
to the numerical height, and \(F_W\le X\). Enumerate the three
tails and fix them at a total cost of at most \(3/R\) in
normalized exponent, plus an arbitrarily small inner allowance.
The fixed tails are pairwise coprime, because they divide the
three letters of an original primitive row.

To align residues modulo \(F_W\), its own role already has
residue zero. On each of the other two roles, a residue class
modulo \(F_W\) determines the class on its opposite role
bijectively, by the signed equation. If either class contains
fewer than \(H=n^{c_0-\eta}\) available labels, all the pairs of
this kind contribute at most \(2H n^{c_0+\kappa}\) edges.
Each such edge is determined by its two opposite values.
Otherwise, there are at most \(n^{\eta+\kappa}\) class pairs.
Choose one of these pairs with the largest remaining edge
count. This fixes both opposite residues and preserves those
fixed in earlier steps.

For the order of choices, suppose the family before a residue
alignment has count exponent at least \(2c_0-d\). Choose a
fixed \(\eta\) with \(\eta-\kappa>d\), then take the height
large at the fixed budget. The discarded edges have strictly
smaller exponent, and the surviving class has exponent at
least \(2c_0-d-\eta-\kappa\) up to the inner allowance. Apply
this three times, each time putting the next threshold above
the accumulated deficit. For example, if \(\chi>0\) bounds
the inner losses, start with
\(d_0=\tau+3/R+\chi\) and use
\(\eta_j=2(d_{j-1}+\kappa+\chi)\) and
\(d_j=d_{j-1}+\eta_j+\kappa+\chi\).
All thresholds can be below \(c_0\) for sufficiently large
fixed \(R\) and small input losses. The total loss is
\(O(\tau+1/R+\kappa+\chi)\). The Chinese remainder theorem
combines the three sets of classes modulo
\(G=F_AF_BF_C\).

At each fixed continuity cutoff, the tail weight tends to
\(h_W-\int_{[0,R)} r\,d\nu_W(r)\), by convergence of the full
height and of the bounded truncated integral. This quantity
decreases to \(e_W\) as \(R\) increases. Use an increasing
sequence of such cutoffs, and at each cutoff select a late
enough small fixed budget where the tail weights, input
deficits and all chosen tolerances are as close to their
limits as required. The cutoff, residue thresholds and
inner exponent allowances are fixed before the height
limit at that budget. Along this diagonal selection the
count loss tends to zero and the numerical height of \(G\)
tends to \(e_{\rm tot}\). In particular, (KP5) and (KP7) can
be applied for any fixed \(0\le\beta<e_{\rm tot}\) on this
refinement. Since the retained boxes have the original
limits, the resulting finite-degree constraints concern
their original limiting measures.

There is also an absolute \(\theta>0\) for which
\[
                        e_{\rm tot}\le \min_W h_W-\theta       \tag{KP8}
\]
at a minimum-dictionary endpoint. Indeed fix a complete value
on a role different from \(W\); averaging over its dictionary
gives a fibre of exponent at least \(c_0-o(1)\). The complete
\(W\)-values on that fibre are in a progression of step \(G\)
and span at most \(X^{h_W+o(1)}\). If
\(e_{\rm tot}>h_W-\theta\), the progression parameters are
bounded by \(X^{\theta+o(1)}\). Choose \(\theta\) strictly
smaller than the rank-one feature threshold for fixed-fibre
exponent \(r=c_0/2\) and a fixed upper bound on the raw caps.
The fixed-fibre rank-one count contradicts the retained
exponent. This argument allows the step and the first
progression value to vary with height and budget.
[Small affine features in a fixed fibre](ABC_FIXED_FIBRE_TRIANGLE_RADICALS.md#a-small-feature-bound-on-a-single-fibre)

# Uniform private valuation bands on original shared-value boxes

A moving complete letter on a fixed primitive fibre cannot hide its
outside radical at very low valuations. This remains true after fixing
an arbitrary supported part of the letter. There is also a short-interval
count for the corresponding cofactor, uniform in the prescribed
multiplier and the position of the interval. Applied to original boxes,
these estimates give a positive numerical reserve on every edge
outside all the other labels.

## A positive valuation band on a fixed-letter fibre

Fix \(C,r>0\), put \(n=X^L\), and count original primitive rows
\(A+B=C_{\rm row}\) with \(X/2<C_{\rm row}\le X\) and each individual
radical at most \(n^C\). Fix one complete letter \(T\), a different
letter position \(W\), and an arbitrary positive integer \(Q\le X\).
Keep only \(Q\mid W\). Set \(U=W/Q\) and define
\[
    \nu_{W,Q}=\sum_{\substack{p^e\parallel W\\p\nmid Q}}
                     \frac{\log p}{\log n}\,\delta_{Le}.
\]

**NB1.** There are \(0<a<b\), \(\gamma>0\) and \(L_0>0\), depending
only on \(C,r\), such that at every fixed \(0<L<L_0\) the number of
these rows with \(\nu_{W,Q}([a,b])\le\gamma\) is
\(O_{C,r,L}(n^{r/2})\), uniformly in \(T,Q\) and the two roles.

If this fails, test the bands \([1/j,j]\) and thresholds \(1/j\) and
choose violating budgets \(L_j<1/j\). At each fixed violating budget
there are arbitrarily large heights where the exponent is at least
\(r/2-o(1)\). Use the fixed-letter refinement for (FT1), with
\(\operatorname{rad}Q\) as the moving letter's prescribed prime pool.
It is supported on an integer at most \(X\). After the inner
refinements and passage to a subsequence, the outside measure of \(W\)
has nonatomic mass at least \(r/2\). The measure at the sampled
heights is exactly \(\nu_{W,Q}\). For every compact interval in
\((0,\infty)\) its mass tends to zero in the diagonal sequence.
The measures are tight, since their masses are at most \(C\) and their
first moments are at most one. The limit is supported at zero,
contradicting its nonatomic mass. At a fixed budget the \(O\)-constant
can absorb the finite initial range of heights. A finite union over
the choices of roles gives the asserted uniform statement.
[Outside valuation mass on an actual fixed-letter fibre (FT1)](ABC_FIXED_FIBRE_TRIANGLE_RADICALS.md#nonatomic-mass-outside-prescribed-prime-pools)

The \(Q\)-outside part of \(W\) divides \(U\). In particular a row
with \(U\le X^{h_0}\), where \(h_0=a\gamma/2\), has
\(\nu_{W,Q}([a,b])\le h_0/a<\gamma\). If instead we start with a prime
pool supported on \(P\le X^B\) for some fixed \(B\), first record the
entire \(P\)-supported part \(Q\) of \(W\) at cost
\(O_{B,\varepsilon}(X^\varepsilon)\). Its value is at most \(X\),
independently of \(B\). Thus the same valuation band, threshold and
cutoff apply to fibres with any such pool. Only the enumeration
constant depends on \(B\).
[Uniform supported enumeration](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)

## Short cofactor intervals of arbitrary centre

**NB2.** There are \(\theta>0\) and \(L_1>0\) depending only on \(C,r\)
such that for each fixed \(0<L<L_1\), at most
\[
                  O_{C,r,L,\varepsilon}(n^{r/2}X^\varepsilon)
\]
rows in the above fibre have \(U=W/Q\) in any one interval of length
\(X^\theta\). The centre, \(Q\) and \(T\) may vary arbitrarily. An
interval of length \(QX^\theta\) in the complete \(W\)-coordinate
therefore has this bound in the subfamily with \(Q\mid W\).

Count \(U\le X^{h_0}\) by (NB1). If the interval meets the part
\(U>X^{h_0}\), its centre is at least \(X^{h_0}/2\) for
\(\theta<h_0/2\) and large \(X\). The short interval radical count
(SI2) with \(\eta=r/4\) gives at most \(n^{r/4}X^\varepsilon\)
possibilities for \(U\), because \(U\le X\) and
\(\operatorname{rad}U\le n^C\). Each determines \(W\) and then the
row. The argument includes intervals that cross the cutoff and
ones centred at zero.
[Short intervals with a positive-height centre (SI2)](ABC_SHORT_INTERVAL_RADICAL_FILTER.md#arbitrarily-small-normalized-counts-in-shorter-intervals)

## Uniform application to all edges of a box

Use the actual ordered \(m\)-boxes of (LB1), with \(m\ge2\), class
bound \(K\), residue bounds \(R_i\), and \(R=R_1R_2R_3\). Write \(L_e\)
for a complete label and \({\cal P}_{e,P}\) for its private radical
relative to an external polynomial-height prime pool \(P\), as in
(PV1). Let \(M_{e,P}\) be the product of the full exact powers in
\(L_e\) on the support of \({\cal P}_{e,P}\).

**NB3.** For any edge \(e\) between axes \(i,j\), the number of
boxes in which
\[
 \sum_{\substack{p^v\parallel L_e\\ p\mid {\cal P}_{e,P}\\
                  Lv\in[a,b]}}\frac{\log p}{\log n}\le\gamma
\]
is
\[
  O_{m,C,r,B,L,\varepsilon}
   \left(
         X^\varepsilon \frac{R^m K^{3m-2}}{\max(R_i,R_j)}
                              n^{r/2}
   \right),\qquad P\le X^B.                                   \tag{NB3}
\]
In particular, suppose a cohort of boxes has at least
\(R^m K^{3m-1}n^{-r/8}\) members and \(K\max(R_i,R_j)\ge n^r\)
for each axis pair. At all sufficiently small fixed budgets its
members have, with the exception of a fraction \(O(n^{-r/4})\),
normalized radical mass greater than \(\gamma\) on private
primes in the band for **every** edge. Consequently
\(M_{e,P}>X^{a\gamma}\) on every edge of each retained box.
The band and the positive numerical reserve \(a\gamma\) are
independent of \(m\), of the packets and of the outside profiles.

Omit the endpoint on the larger residue axis and fix all other
residues and an anchor spanning tree, at cost at most
\(R^m K^{3m-2}/\max(R_i,R_j)\). The anchor pool consists of the labels
between anchors, the nonzero differences between any two anchored
neighbors of the omitted vertex, and \(P\). Every prime on the
moving edge absent from the pool is private relative to \(P\),
because a common prime with another incident edge divides one
of those neighbor differences. The pool has polynomial height
with a fixed exponent at each \(m\). Fix the entire part \(Q\)
supported on this pool. Complete the edge to a triangle with an
anchored vertex on the third axis. Its opposite edge is a fixed
complete letter and these are actual primitive fibres. The tested
mass for \(\nu_{W,Q}\), with \(W=L_e\), is bounded above by the
private band mass. Now apply (NB1), then restore \(Q\) at subpower
cost. The complete moving value determines the omitted potential.
This proves the displayed estimate. Divide by the cohort bound
and take \(\varepsilon<rL/8\) to sum the exceptional counts over
the finitely many edges. On each retained edge the first moment
on the band exceeds \(a\gamma\).
[Recovery of an omitted potential and its prime pool (PV1)](ABC_SHARED_BOX_PRIVATE_PRIME_SUPPORT.md#a-bound-for-one-position)
· [Actual box lower bound (LB1)](ABC_SHARED_VALUE_RESIDUE_BOXES.md#a-finite-difference-array-theorem)

The same numerical saving applies to a covering by short
cofactor intervals chosen separately from each anchor record.
For example, with a fixed \(Q\) that may depend on the anchors,
\(n^{r/8}\) intervals of length \(QX^\theta\) in \(L_e\) cover at most
\(O(n^{-r/8})\) of the cohort with \(Q\mid L_e\), measured relative
to the whole cohort. Multiply the anchor count by the list size
and (NB2). The conclusion also holds when \(Q\) is the whole
anchor-pool supported part and the lists depend on \(Q\), because
these parts can be enumerated at subpower cost. This gives a
real dispersion constraint on the original cofactor with a
moving interval centre and a possibly large multiplier.

## Consequence for derivative lattice tests

For every fixed \(C,r>0\), choose \(m\) so large that
\[
                  3(m-1)^2a\gamma>3m-2 .
\]
Work with the stated cohort hypotheses for this \(m\). They hold
at the full box exponent in a simultaneous packet stratum with
\(t+\max(s_i,s_j)\ge 2r\) on every pair, taking the fixed
tolerances sufficiently small for this \(m\). Then (NB3) with
\(P=1\) gives every edge a private prime and numerical reserve
\(a\gamma\), away from a power-small fraction. In (DM7) use
\(\omega=a\gamma\) and then take \(L\) sufficiently small.
For these actual boxes, every derivative coordinate set that
contains at least one private prime from each edge satisfies
the family derivative inequality with a fixed positive exponent
of slack, even using every eligible independent prime flag.
The largest allowed short-letter losses are covered by the
choice of \(m\); no lower bound on the numerical height of the
short letter within the shell is used. This conclusion concerns
the stated private-prime coordinate sets and the family
inequality. It leaves coordinate choices that omit private
primes on some edges and other arithmetic constraints open.
[Derivative minor cycle counts and general private weight (DM7)](ABC_SHARED_BOX_DERIVATIVE_MINOR_CYCLES.md#geometric-mean-heights-and-general-private-weight)
· [Simultaneous packet scales (LB4)](ABC_SHARED_VALUE_RESIDUE_BOXES.md#consequences-for-simultaneous-packet-strata)

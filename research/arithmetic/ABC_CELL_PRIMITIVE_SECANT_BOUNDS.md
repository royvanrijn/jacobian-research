# Primitive secant bounds inside simultaneous packet cells

A fixed complete letter gives enough nonatomic mass on either
free letter to exclude a small determinant feature model. This
holds at any fixed positive fibre exponent. Applied to actual
three-packet cells, it bounds short secants at the smaller cell
scale and yields a new radical floor of six times the cell
dictionary exponent. The same test can be used on pairs of
shared-value boxes with matching residue lists.

## Small features on a fixed complete-letter fibre

Let \(n=X^L\), and fix \(C,B>0\), a dimension \(R\ge1\),
and an allowance \(\alpha>0\). There are \(\theta_0,L_0>0\)
such that the following holds at each fixed \(0<L<L_0\).
Take positive primitive triples \(A+B=C_{\rm out}\) in
\(X/2<C_{\rm out}\le X\) with original individual radical
caps \(n^C\). Fix one complete letter \(T\), and take
one of the two free letters, say \(W\). Suppose that each
row is assigned \(Z\in\mathbb Z^R\) with
\(\|Z\|_\infty\le X^{\theta_0}\), and one common real affine
function \(J\) satisfies
\[
               J(Z)\in W\mathbb Z\setminus\{0\},
                    \qquad |J(Z)|\le X^B .
\]
The number of rows is
\[
                           O_{C,B,R,\alpha,L,\varepsilon}
                              (n^\alpha X^\varepsilon).
                                                               \tag{CS1}
\]
The coefficients of the affine function are arbitrary.
The integer-multiple condition is on the values on the
actual primitive rows.

Suppose (CS1) fails uniformly. Extract a sequence of
fixed budgets tending to zero and prospective feature
exponents tending to zero. At each budget take the
height limit on the violating sequence before passing
to the next budget. The normalized count exponent of
the starting fibre is at least \(\alpha-o(1)\).
After subsequences the current affine rank can be fixed
at each of the finitely many stages below.

Fix \(p>0\) with \((R+1)p<\alpha/10\). On the current
family choose a feature witness \(Z_0\). The value
\(J(Z_0)\) is a nonzero polynomial-height integer.
Include its prime support in a prescribed pool for
\(W\); one may include the support of any other
polynomial-height record as well. Fix the exact full
part of \(W\) on this pool at height-subpower cost,
uniformly in the pool. The fixed-letter theorem (FT1)
applies to the current fibre outside this pool. As
long as the retained exponent is at least
\(\alpha-Rp-o(1)\), the residual exact scaled valuation
measure on \(W\) has at least \(\alpha-Rp\) nonatomic
mass in the limit, after the fixed-budget anatomic
refinements.
[Nonatomic mass outside a witness pool (FT1)](ABC_FIXED_FIBRE_TRIANGLE_RADICALS.md#nonatomic-mass-outside-prescribed-prime-pools)

Choose a compact continuity band \([s,S]\subset(0,\infty)\)
with residual nonatomic mass exceeding \(2p\).
Delete small continuity neighborhoods of the finitely
many atoms in the band of mass at least \(p/4\), preserving
this lower bound. The individual prime masses on the
retained part are less than \(p/2\) at all sufficiently
late budgets. Otherwise some sequence of such primes
would give a limiting atom of mass at least \(p/2\)
on the retained compact set. Consequently each row
admits a product \(M'\) of whole prime powers of \(W\)
outside the prescribed pool, of normalized radical cost
between \(p/2\) and \(p\), and numerical weight at least
\(s p/2\) in base \(X\). The supported values of the
products give \(n^{p+o_X(1)}\) possible choices. Fix one
choice, retaining a largest class. The new band and
its positive lower endpoint can be chosen separately
on a subsequence at each stage. Only the cost \(p\)
needs to be fixed across all stages.
[Full prime-power packet selection](ABC_COMPLETE_LETTER_SMALL_FEATURES.md#limiting-setup-and-cost-of-fixed-prime-supports)

If the retained feature vectors span the same affine
space as before, choose a base retained vector and
anchored basis differences with nonzero projected
minor \(\Delta\). Integer Cramer's rule expresses
\(\Delta\) times the old witness evaluation \(J(Z_0)\)
as an integer linear combination of the retained basis
evaluations. The coefficients of \(J\) need not be
rational for this identity. As \(M'\mid J\) on the
retained rows, and \((M',J(Z_0))=1\), this forces
\(M'\mid\Delta\). A minor of rank at most \(R\) has
absolute value at most \(R!(2X^{\theta_0})^R\),
smaller than \(M'\) on the sequence with
\(\theta_0\to0\) at all sufficiently late budgets.
Thus the dimension drops. A drop already caused
by the support or anatomic refinement needs no
further packet.

After at most \(R\) drops \(J\) is constant, so \(W\)
divides a fixed nonzero integer of height \(X^B\).
Each possible \(W\), together with the fixed \(T\),
determines the original row. The divisor bound gives
only \(X^{o_X(1)}\) rows, a contradiction to the
retained positive exponent. This proves (CS1). In
a simultaneous packet cell the original packet product
can be included in each witness pool. The rank drops
then use only the outside radical on the free letter.
[Integer rank descent on fixed-letter fibres](ABC_FIXED_FIBRE_TRIANGLE_RADICALS.md#a-small-feature-bound-on-a-single-fibre)

## The secant inequality at the cell scale

Let \(\mathcal E\) be a family in the same shell and
original radical caps. Write \(E=|\mathcal E|\), and
suppose each complete dictionary has at most \(K\)
entries. View the rows as points \((A,B)\).
For distinct points \(x,y\), define their primitive
direction \(v=(x-y)/\gcd(x_1-y_1,x_2-y_2)\) using a
positive gcd, and set
\[
              H(x,y)=\max(|v_1|,|v_2|,|v_1+v_2|).
\]
When all three entries are nonzero, their absolute
values form a positive primitive additive triple, after
placing its largest value last. Write \(R(x,y)\) for its
radical. Choose a fixed \(B_0>0\) and a prescribed pool
\(P\le X^{B_0}\), enlarging \(B_0\) if necessary to bound
the raw radical cap as well.

For every fixed \(\alpha,\xi>0\) and sufficiently small
fixed \(L\), choose \(\theta>0\) in terms of the
preceding fixed parameters. The number of ordered pairs
that have \(H(x,y)\le X^\theta\), or have larger direction
height, raw radical \(R(x,y)\le n^{B_0}\) and outside
\(P\) radical at most \(n^b\), satisfies
\[
 \ll E+
    E^{3/2}\bigl(K^{1/2}n^{\alpha/2}+n^{b/6+\xi}\bigr)X^\varepsilon
    + E K ,
                    \qquad 0\le b\le B_0.                 \tag{CS2}
\]
The implicit constants are uniform in the family, pool and
cap \(b\). The small-budget threshold can depend on the
fixed normalized allowance \(\xi\). Pairs with a zero direction entry have
height one and are already included.

For short directions, fix distinct anchors \(y,z\).
A query \(x\) off their line with both its directions
short has \(u=v(x,y)\), \(v=v(x,z)\) and
\(d=[u,v]\ne0\). The exact intersection identity is
\[
                           dx=dy+[z-y,v]u .
\]
It expresses \(dA,dB,dC_{\rm out}\) as fixed linear
functions of
\[
                   (d,u_1v_1,u_1v_2,u_2v_1,u_2v_2).
\]
These integer features have absolute value at most
\(2X^{2\theta}\), and each multiple has polynomial
height with a uniform exponent. Choose \(\theta\)
small enough for (CS1). Fix any one complete letter
and use (CS1) on a free one, then sum over at most
\(K\) fixed values. The line through the anchors has
at most \(K\) points since some complete letter is
nonconstant and injective on the line. Thus two
distinct anchors have at most \(O(K n^\alpha X^\varepsilon)\)
common short queries. If \(S\) counts the short pairs
and \(d_s(x)\) counts short anchors at a query, then
\[
          S^2/E\le\sum_x d_s(x)^2
                        \ll S+E^2 K n^\alpha X^\varepsilon.
\]
Solving gives the short term of (CS2). The direction
identity is the same as at the global scale. The local
argument uses the fixed-fibre estimate to avoid a
global near-saturation hypothesis.
[Two short directions at the global scale (SC1)](ABC_PRIMITIVE_SECANT_RADICAL_AMPLIFICATION.md#a-uniform-bound-for-queries-short-in-two-directions)

At height at least \(X^\theta\), the effective radical
count (EB4) gives at most \(O(n^{b/3+2\xi}X^\varepsilon)\)
permissible primitive directions. Indeed one uses the
positive primitive output shell at height \(Y\ge X^\theta\),
pricing the fixed pool at zero. The possible relative raw
budgets are \(O_{B_0,\theta}(L)\). Taking a sufficiently fine
fixed mesh of multiples of \(L\) for the shell conversion
gives effective exponent \(b/3\) in base \(n\) up to the
prescribed allowance.
[Shell conversion and direction count (SC2)](ABC_PRIMITIVE_SECANT_RADICAL_AMPLIFICATION.md#a-pair-inequality-on-the-full-balanced-family)

For \(M_0\) allowed directions, each point is in at most
\(M_0\) lines of those directions. Lines containing fewer
than \(k\) input points therefore contribute at most
\(EM_0 k\) ordered pairs. The Szemerédi–Trotter line
bound \(O(E^2/k^3+E/k)\), summed dyadically for richer
lines of size at most \(K\), gives \(O(E^2/k+EK)\)
further pairs. If \(1\le M_0\le E\), take \(k\) of order
\(\max(2,\sqrt{E/M_0})\); if \(M_0>E\), the resulting
bound \(O(E^{3/2}M_0^{1/2}+EK)\) holds trivially.
This proves the tall term of (CS2).
[Szemerédi–Trotter, Theorem 1](https://trotter.math.gatech.edu/papers/38.pdf)

## Fresh direction radicals on synchronized actual boxes

In a simultaneous packet stratum (SB2), let \(d\) be the
total empirical packet price and \(t=1/6-d\). Suppose
\(t\ge\zeta>0\) for a fixed margin. Each occupied cell
has \(E\ge n^{2t-6\eta}\), \(K\le n^{t+8\eta}\).
Fix \(0<\sigma<6\zeta\). With the tolerances small
enough, (CS2) with \(\alpha=\zeta/4\), \(b=6t-\sigma\),
small \(\xi\), and small inner height allowances gives a fixed power
saving in the fraction of cell pairs that are identical,
short, or pass the outside radical test.
The raw cap \(B_0\) can be any fixed constant.
For example, all but a power-small fraction of cell
pairs have tall primitive direction radical greater
than \(n^{6t-\sigma}\).
[Simultaneous cell bounds (SB2)](ABC_SIMULTANEOUS_COMPOSITE_PACKET_BLOCKS.md#simultaneous-marking-theorem)

Now also assume that all three packet prices are bounded
below by a fixed positive constant. Fix \(m\ge2\) and
consider ordered pairs of actual shared-value arrays
of size \(m\) with the same ordered residue lists.
If \(B_{\mathcal Q}\) counts the arrays on a residue
box \(\mathcal Q\), (LB1) and Cauchy–Schwarz give
\[
 \sum_{\mathcal Q} B_{\mathcal Q}^2
                    \ge n^{md+(6m-2)t-O_m(\kappa+\eta)}
                                                               \tag{CS3}
\]
up to fixed factors and arbitrarily small inner allowances.
Indeed the count of residue boxes is at most
\(n^{md+O_m(\kappa+\eta)}\). This includes repetitions
of an array, whose later rejection is covered by the
identical-row test.
[Actual shared-value arrays (LB1)](ABC_SHARED_VALUE_RESIDUE_BOXES.md#a-finite-difference-array-theorem)

For a specified corresponding cell, fixing a candidate
pair of original rows fixes two spanning-tree edges in
each array: its three labels form a triangle of the
potential graph. The other edges of a spanning tree
cost at most \(K^{6m-6}\) choices, where one may take
\(K\le n^{t+8\eta}\) to bound every complete pair class
simultaneously. The uniform cell bound counts at most
\(K^4 n^{-\delta}\) bad candidate pairs for some fixed
\(\delta>0\). Here identical row pairs alone contribute
at most \(K^2\), which satisfies the same bound after
decreasing \(\delta\). Consequently a power-small
fraction of the synchronized pairs in (CS3) have a
bad corresponding secant in any of their \(m^3\) cells.

## Progression widths and the three pooled valuation tests

For every fixed \(\mu>0\), a cell with \(t\ge\zeta>0\)
can be pruned by a power-small fraction of its rows so
that each complete dictionary meets every ordinary
progression of at most \(X^\theta\) terms in
\[
                          O(n^\mu X^\varepsilon)
                                                               \tag{CS4}
\]
values. The positive \(\theta\) and small-budget threshold
are uniform when the cell tolerances are sufficiently
small in terms of \(\mu,\zeta\). The step and the initial
term of the progression may vary with height and cell.

In fact, a subfamily in which one complete role \(W\)
lies in such a progression can be counted after fixing
any other complete letter \(T\). The progression index
is an integer feature of height at most \(X^\theta\)
and \(J=W\) is affine in it. Use (CS1) with \(R=1\)
and an allowance \(\alpha<\min(\mu,\zeta)/4\).
Sum over the at most \(K\) fixed values of \(T\),
obtaining \(O(Kn^\alpha X^\varepsilon)\) rows.
Iteratively remove vertices of degree below
\(n^{t-\delta}\), with \(14\eta\ll\delta\ll\mu,\zeta\).
The removal costs at most \(3Kn^{t-\delta}\) rows,
a power-small fraction of \(n^{2t-6\eta}\).
Each remaining vertex contributes at least
\(n^{t-\delta}\) rows, so the claimed progression
intersection count follows. The same row count
works for an exact affine description with any
fixed number of small integer features, by (CS1).
In particular, if \(FGM>X^{1-\theta/2}\), every
complete role in the cell lies in a progression
with at most \(X^\theta\) terms for large \(X\).
This gives a power saving over the cell scale
\(n^{2t}\) on that sector. Summing over the packet
and residue records gives a saving over \(1/3\) when
that numerical condition holds on the whole selected
stratum.

There is also a profile-dependent inequality that
uses all three packet weights. Consider a height-first
limiting sequence of actual cells with row exponent
at least \(2t\) and all three dictionary caps at most
\(t+o(1)\), with \(t>0\). Work in original anatomic
boxes throughout. Let \(h_W\) be the limiting numerical
height exponent of complete role \(W\), and let \(w\)
be the limiting \(\log_X(FGM)\). Write \(\rho_W\) for
the weak scaled valuation measure of the complete prime
powers of \(W\) outside the packet product, and \(e_W\)
for its escaped numerical weight, defined as the
limiting outside numerical height minus
\(\int r\,d\rho_W(r)\). Pass to subsequences on these
fixed data. For a finite measure \(\omega\) on
nonnegative scaled degrees, put
\[
 H_\omega(p)=\sup_{\substack{0\le\omega'\le\omega\\
                              \omega'([0,\infty))\le p}}
                    \int r\,d\omega'(r).
\]
Then for every pair of distinct roles \(U,V\),
\[
      w+e_U+e_V+
      H_{(\rho_U)_{\rm na}+(\rho_V)_{\rm na}}(t)
                         \ \le\ \min(h_U,h_V).             \tag{CS5}
\]
The measures \((\rho_W)_{\rm na}\) here are the
nonatomic parts. This uses the valuation spectra
after removing full prime powers, including the
packet on the third role in the numerical term \(w\).

To prove (CS5), choose a complete value of the remaining
role \(T\) whose original fixed-letter fibre has exponent
at least \(t-o(1)\), by the cell row and part caps.
All rows still have the fixed packet and residue data.
The residual profiles on \(U,V\) are unchanged by
this restriction, since in the original boxes their
weak measures are fixed up to vanishing tolerances
and the packet on each role is common. Equation (FT1)
also gives residual nonatomic mass at least \(t\)
on each free letter.
[Nonatomic mass outside the packet pool (FT1)](ABC_FIXED_FIBRE_TRIANGLE_RADICALS.md#nonatomic-mass-outside-prescribed-prime-pools)

If (CS5) fails by a fixed margin, choose \(p<t\) and
a bounded positive degree selection of combined
residual nonatomic mass at most \(p\) and weight
approaching \(H_{(\rho_U)_{\rm na}+(\rho_V)_{\rm na}}(t)\)
within that margin. Integrability permits the bounded
cutoff, and continuity in \(p\) permits the strict
mass margin. Whole powers in compact continuity
bands realize this selection with arbitrarily small
losses. Remove neighborhoods of a finite list of
atoms leaving arbitrarily small unlisted atomic
mass, then round the retained selection at negligible
individual prime mass. In addition retain the entire
tails above a sufficiently large continuity cutoff
on both free roles. The tail radical cost is arbitrarily
small by the first moment bound, and their combined
numerical weight is at least \(e_U+e_V-o(1)\).
The cutoff is above the bounded selection and all
the original packet orders. The selected records
therefore cost less than \(t\) by a fixed margin.
Fix their complete values by ordinary enumeration
of their joint radical, prime allocations and
supported powers.

The selected full powers on \(U,V\) are coprime by
original primitivity and are outside \(FGM\). Forcing
them on their two roles and using the exact relation
with the fixed \(T\) gives one CRT residue for the
shorter of \(U,V\), together with its original
residue modulo \(FGM\). The combined modulus is
larger than the entire range of that letter by the
assumed failure of (CS5). Thus there is at most one
original row per record. This contradicts the pinned
fibre exponent and proves the inequality for that
pair. The three choices of \(T\) can be made independently.
[Fractional upper moments of two degree measures](ABC_WIDE_COFACTOR_CONGRUENCES.md#a-pooled-degree-inequality-and-primitive-shift-fibres)

## Fresh radicals beyond a whole array

Use the synchronized ordered pairs of \(m\)-arrays in
(CS3). Fix margins \(0<\sigma<6\zeta\) and \(\tau>0\).
Put \(b=6t-\sigma\) and
\[
                     r_{\rm tree}=md+(3m-1)t .
\]
A prescribed initial prime pool \(P_0\) may be empty or
have prime product at most \(n^{B_1}\) for a fixed \(B_1\).
For a synchronized pair let \(P_{\rm first}\) be the
union of \(P_0\) and the supports of all complete labels in
the first array. For each diagonal cell \((j,j,j)\) denote
the corresponding primitive direction radical by \(R_j\).
For each spanning tree \(\mathcal T\) of the tripartite
potential graph of the first array, denote the product of
its complete labels by \(D_{\mathcal T}\).

**CS6.** On all but a power-small fraction of the
synchronized pairs, the diagonal directions are tall
and have all entries nonzero. Also, simultaneously for
every tree \(\mathcal T\),
\[
 \begin{split}
  \prod_{\substack{p\mid D_{\mathcal T}\\p\notin P_0}}p
                        &> n^{r_{\rm tree}-\tau},\\
  \prod_{\substack{p\mid \prod_{j=1}^m R_j\\p\notin P_{\rm first}}}p
                        &> n^{mb}.
 \end{split}                                                   \tag{CS6}
\]
Both inequalities are unconditional on the raw direction
radicals. In particular the joint radical of the first
tree labels and the \(m\) primitive direction triples,
outside \(P_0\), exceeds \(n^{r_{\rm tree}-\tau+mb}\).
The field or tree used in a subsequent application may
be chosen after seeing the arrays.

If one only assumes polynomial height of the prime product
of \(P_0\), the second inequality holds on the branch where
all the \(R_j\) have a prescribed raw cap \(n^{B_0}\), with
the exceptional count still measured against all synchronized
pairs in (CS3). In that branch, at each step \(j\) the radical of the current
direction outside \(P_{\rm first}\) and the earlier \(R_i\)
already exceeds \(n^b\).

### Rich lines and a pool fixed by one anchor

First remove synchronized pairs having an identical row,
a short direction, or a line of richness at least
\(k=\lceil n^\omega\rceil\) in any corresponding cell.
Here \(\omega>0\) is fixed and small in terms of \(\sigma,\zeta\),
and richness is measured in the whole ambient cell of (SB2).
In one cell the short-pair estimate in (CS2), the \(E\)
identical pairs and the rich-line bound
\[
                           O(E^2/k+EK)
\]
together give \(O(K^4 n^{-\delta_1})\) pairs for some
fixed \(\delta_1>0\). Indeed \(E\le K^2\), the short pairs
have a power-small fraction by (CS2), and \(K\ge
n^{t-3\eta}\) from the cell edge lower bound. Each prescribed
pair of rows fixes two spanning-tree edges per array, so
there are at most \(K^{6m-6}\) further choices at fixed
residue lists. Summing over lists and cells makes these
exceptions power smaller than (CS3), by taking the
packet tolerances small after fixing \(\omega\).

Now fix the ordered residue lists and the entire first
array. Expose the diagonal rows of the second array in
order. At stage \(j\) put the primes of \(P_{\rm first}\)
and \(R_1,\ldots,R_{j-1}\) into the pool \(P_j\).
This pool is fixed for the current query once the preceding
diagonal rows are known. Its prime product is of polynomial
height uniformly, since the original complete labels and
the primitive direction coordinates are at most \(X\).
At a fixed first row of the current cell, the tall direction
count from (EB4) has at most \(O(n^{b/3+\xi}X^\varepsilon)\)
directions with outside-\(P_j\) radical at most \(n^b\)
and a fixed raw cap \(n^{B_0}\). Each such direction gives
fewer than \(k\) second rows after the rich-line exclusion.
Since \(b/3=2t-\sigma/3\), the bound
\(O(k n^{b/3+\xi}X^\varepsilon)\) is \(O(K^2n^{-\delta_2})\)
for a fixed \(\delta_2>0\) by taking
\(\omega,\xi,\eta,\varepsilon/L\) sufficiently small.

There are at most \(K^{3m-1}\) first arrays at fixed
residue lists. The \(j-1\) prescribed second-array
diagonal rows cost at most \(K^{2(j-1)}\). The bad
current row costs \(O(K^2n^{-\delta_2})\). The \(j\)
diagonal triangles have disjoint potential vertices.
Select two edges per triangle and extend this forest
to a spanning tree of the second array. The remaining
choices cost at most \(K^{3m-1-2j}\). The resulting
exceptional pair bound is \(O(K^{6m-2}n^{-\delta_2})\)
per residue list and stage. Consequently a power-small
fraction of the synchronized pairs can have such a
bad stage. In particular no count of the possible prime
pools from the first array is needed.

The unconditional second inequality of (CS6) follows
from the stated normalized cap on \(P_0\) and the actual
original radical caps. If its outside joint bound fails,
the support of each \(R_j\) is contained in the union
of a squarefree support of size at most \(n^{mb}\),
\(P_0\), and the supports of the \(3m^2\) original labels
of the first array. Each original label has radical at
most \(n^C\) for fixed \(C\). Thus every raw radical
has cap \(n^{mb+B_1+3m^2C}\). Choose a fixed \(B_0\)
above this exponent before applying the stage test.
If all stages passed, the \(m\) disjoint fresh supports
would have product greater than \(n^{mb}\) outside
\(P_{\rm first}\), contradicting joint failure.

### The tree support and the common-field cost

For a prescribed spanning tree, enumerate the squarefree
support outside \(P_0\) of its label product. If its size is
at most \(n^{r_{\rm tree}-\tau}\), there are at most that
many support choices. On the support and the fixed pool
the exact tree labels cost only \(X^{o_X(1)}\) by uniform
supported-value enumeration at polynomial height. Those
labels reconstruct the entire first array and its ordered
residue lists. For each such first array there are at most
\(K^{3m-1}\) synchronized second arrays, by the tree bound.
This gives at most
\(n^{md+(6m-2)t-\tau+O_m(\eta)}X^{o_X(1)}\) bad pairs,
power fewer than (CS3). There are only a bounded number
of spanning trees at fixed \(m\), giving the first
inequality for every tree.
[Tree reconstruction for actual arrays (LB1)](ABC_SHARED_VALUE_RESIDUE_BOXES.md#from-common-differences-to-shared-values)
· [Supported-value enumeration (EB3)](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)

For a common-field application let \(D=FGM\) and write
\(w=\log_X D\) at the current height. Suppose the labels
of a first-array tree \(\mathcal T\), with height caps
\(X^{h_e}\), and the three entries of each diagonal
primitive direction triple become \(Q\)-full principal
ideals over one number field \(F_{\rm lift}\). The
common positive integer order \(Q\) is fixed before the
height limit at fixed \(L\). The original rows in a
corresponding cell have the same three complete residues
modulo \(D\), so their primitive direction height is at
most \(X/D\). The normalized ideal radical of the entire
product of lifted labels and direction entries is
therefore at most
\[
             n^{\,(\sum_{e\in\mathcal T}h_e+3m(1-w))/(LQ)} .
\]
Apply the local discriminant inequality (RL1) to this
product and use (CS6). On the same generic pairs,
\[
 \log_n \operatorname{rd}(F_{\rm lift}) >
       r_{\rm tree}-\tau + m(6t-\sigma)
          -\frac{\sum_{e\in\mathcal T}h_e+3m(1-w)}{LQ}.
                                                               \tag{CS7}
\]
This charges the original tree and the new primitive
directions in the same field, with no duplication of
their counted prime supports. Neither the field nor its
degree has to be prescribed.
[Local discriminant inequality (RL1)](ABC_RAMIFIED_IDEAL_LIFT_COST.md#the-local-cost-and-a-field-attaining-it)
· [Common-field tree count (RL8)](ABC_RAMIFIED_IDEAL_LIFT_COST.md#a-common-field-on-an-actual-shared-value-box)


## A seed triangle and successive one-vertex directions

Use the synchronized pairs and the positive margins of (CS3).
Fix \(0<\sigma<4\zeta\), \(\tau>0\), and write
\[
 b_0=6t-\sigma,\qquad b_1=4t-\sigma,\qquad
 J=3m-3,\qquad B_{\rm dir}=b_0+Jb_1 .
\]
Use the seed cell \((1,1,1)\), followed by the \(J\) cells
\[
 (i,1,1)\ (2\le i\le m),\qquad
 (1,j,1)\ (2\le j\le m),\qquad
 (1,1,k)\ (2\le k\le m)
\]
in this order. Each later cell introduces just one new
potential and uses the two seed potentials on the other axes.
The original complete label joining those two seeds is
therefore fixed once the second-array seed row has been
recorded. Call this the fixed label at that stage.
Write \(R_0,\ldots,R_J\) for the radicals of the primitive
directions at these cells.

**CS8.** Suppose \(P_0\) has prime product at most \(n^{B_1}\)
for a fixed \(B_1\), and define \(P_{\rm first}\) as in (CS6).
Except for a power-small fraction of all synchronized pairs
in (CS3), the directions in these cells are tall and have
no zero entry, and
\[
 \prod_{\substack{p\mid R_0\cdots R_J\\p\notin P_{\rm first}}}p
                    > n^{B_{\rm dir}} .
\]
On the same generic pairs the first inequality of (CS6)
holds simultaneously for every first-array spanning tree.
In particular the joint radical of that tree and these
directions outside \(P_0\) exceeds
\[
             n^{\,r_{\rm tree}-\tau+B_{\rm dir}},
 \qquad B_{\rm dir}=(12m-6)t-(3m-2)\sigma .
 \tag{CS8}
\]
If \(P_0\) is only assumed to have polynomial height, the
same conclusions apply on the branch where all directions
satisfy a prescribed raw radical cap \(n^{B_0}\), with the
exceptional count still relative to all pairs in (CS3).
On that branch the seed contributes fresh radical greater
than \(n^{b_0}\), and the two nonfixed direction entries
of each later stage contribute fresh radical greater than
\(n^{b_1}\). The counted prime supports at different stages
are disjoint.

### Counting the seed and the one-vertex stages

First discard synchronized pairs having an identical row or a
short direction in any of the selected cells. Also discard
those whose seed secant line has ambient cell richness
at least \(k=\lceil n^\omega\rceil\), with fixed small
\(\omega>0\). As in (CS6), (CS2) and the rich-line estimate
\(O(E^2/k+EK)\) give at most \(O(K^4n^{-\delta})\) bad pairs
of rows in any of these tests, where \(\delta>0\) is fixed
and packet tolerances are sufficiently small. The fact
\(K\ge n^{t-3\eta}\) follows from the cell edge lower
bound. A pair of original rows fixes two spanning-tree
edges in each array. Summing the resulting continuation
bound over the finitely many tests and all residue lists
gives a power-small fraction of (CS3). The remaining
directions have no zero entry because any primitive
direction with a zero entry among \(v_1,v_2,v_1+v_2\)
has height one.

Fix the ordered residue lists and the whole first array.
Keep two edges of the second seed triangle. At each later
stage record one edge joining the new potential to a known
seed potential. Those edges determine the new potential,
hence the whole current row. The recorded edges over the
whole sequence form a spanning tree of the second potential
graph, since each new edge attaches a single new vertex to
the connected known graph. The first array costs at most
\(K^{3m-1}\) choices at fixed residue lists, and the second
array costs \(K^2\) at the seed and \(K\) per later stage
before any direction restrictions.

In each test put \(P_{\rm first}\) and all prime supports
of the previously exposed primitive directions into the
prescribed pool. The pool is fixed once the current
prefix is fixed, and has uniform polynomial height. At the
seed, apply the tall direction bound (EB4) of (CS6) with
outside cap \(n^{b_0}\) and raw cap \(n^{B_0}\). This gives
\(O(n^{b_0/3+\xi}X^\varepsilon)\) possible directions at
the fixed first row, with any small fixed allowance \(\xi\).
Each yields fewer than \(k\) second choices by the seed
rich-line test. Because \(b_0/3=2t-\sigma/3\), the bad
seed costs \(O(K^2n^{-\delta_0})\) for some \(\delta_0>0\).

At a later stage let \(\Delta_T\) be the signed difference
in the known complete label between the two arrays, using
the same sign convention as for the primitive direction.
This nonzero integer is fixed before the query. The
absolute value of the corresponding primitive direction
entry \(v_T\) is a divisor of \(|\Delta_T|\), giving
\(X^{o_X(1)}\) possibilities. We may add the support of
\(\Delta_T\) to the stage pool so that this known direction
entry is fully supported on the enlarged pool. For a given
signed primitive direction, \(g=\Delta_T/v_T\) determines
the positive scaling and thus at most one second row.

If the product of distinct primes in the two remaining
direction entries outside the enlarged pool is at most
\(n^{b_1}\), one of those entries has an outside radical
at most \(n^{b_1/2}\), by primitive coprimality. Fix the
value \(|v_T|\), its position in the positive primitive
triple, and a dyadic shell for the direction height. The
direction height is larger than \(X^\theta\) after the
preliminary pruning. Apply (FL1) in that shell to the
fixed complete direction letter and the two outside caps,
with the same fixed \(L\). More explicitly, if the shell
upper endpoint is \(Y\), use the original raw cap
\(n^{B_0}\le Y^{LB_0/\theta}\) and the smaller outside
cap \(Y^{L(b_1/2)\log X/\log Y}\). The prescribed pool
has bounded polynomial height in \(Y\). The exponent
from (FL1) is at most \(b_1/4+\xi\) in base \(n\),
with any prescribed small positive allowance \(\xi\).
The divisor choices, signs and shells are absorbed in
the inner height allowance. Thus the bad current row
costs \(O(n^{b_1/4+\xi}X^\varepsilon)\). This is
\(O(Kn^{-\delta_1})\) for some \(\delta_1>0\), since
\(b_1/4=t-\sigma/4\) and \(\eta,\xi,\varepsilon/L\)
can be chosen sufficiently small.
[Prescribed-pool fixed-letter count (FL1)](ABC_FIXED_LETTER_RADICAL_CAPS.md#a-pointwise-count-with-prescribed-prime-pools)

For each stage the other second-array records cost at
most their unpruned product of \(K^2\) for the seed and
\(K\) for a single-vertex stage. The exception count per
stage and residue list is therefore
\(O(K^{6m-2}n^{-\min(\delta_0,\delta_1)})\), and
summing over the residue lists and stages makes it power
smaller than (CS3). The recorded supports outside each
stage pool are disjoint because all earlier directions
were included in the subsequent pools. All these
supports lie outside \(P_{\rm first}\), even when extra
known difference primes were also removed. This proves
the capped-branch conclusions.

For the unconditional conclusion, if the joint radical
outside \(P_{\rm first}\) is at most \(n^{B_{\rm dir}}\),
each raw direction radical is at most
\(n^{B_{\rm dir}+B_1+3m^2 C}\), where \(n^C\) caps each
complete original radical. Choose a fixed \(B_0\) larger
than this bound's exponent, uniformly in the allowed
\(t\), before applying the capped-branch estimates.
The seed and the later stages then force a disjoint
union of radical greater than \(n^{b_0+Jb_1}\) outside
\(P_{\rm first}\), contrary to the joint failure.
Intersect this result with the generic tree-support
estimate of (CS6), whose proof does not depend on the
choice of direction cells. This gives (CS8).

### The common-field cost with successive directions

Write \(D=FGM\), \(w=\log_X D\). Suppose a first-array
tree has label height caps \(X^{h_e}\). If its labels
and all three entries of each of the \(J+1\) directions
become \(Q\)-full principal ideals in a common field
\(F_{\rm lift}\), the generic pairs in (CS8) satisfy
\[
 \log_n\operatorname{rd}(F_{\rm lift}) >
 r_{\rm tree}-\tau+B_{\rm dir}
  -\frac{\sum_e h_e+3(J+1)(1-w)}{LQ}.
\]
Indeed every direction entry is at most \(X/D\) since
both rows lie in the same complete residue classes
modulo \(D\). Use the \(Q\)-full bound on the
normalized ideal radical of their product and apply
(RL1), as in (CS7). Under the prescribed raw direction
caps, one may use just the three seed direction entries
and the two nonfixed entries of each later direction.
The same bound then holds with \(3(J+1)\) replaced by
\(3+2J\), and fullness of the omitted entries is not
needed on that branch.
[Local discriminant inequality (RL1)](ABC_RAMIFIED_IDEAL_LIFT_COST.md#the-local-cost-and-a-field-attaining-it)

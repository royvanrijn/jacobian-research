# One high power grid in a graph of actual additive values

Two values from the same letter can supply the two powered terms of
a fixed-shift equation. As a result, a dense additive graph can be
counted when only one letter has a high valuation grid. The gain
depends on the sizes of the two other complete-letter dictionaries.
It does not require them to have valuation grids or to have prescribed
limiting degree profiles.

## Difference multiplicities of one grid dictionary

Fix \(f\ge0\) and \(\alpha>1/2\). Set \(n=X^L\) and
\(K=\lceil L^{-\alpha}\rceil\), with \(L\) fixed before the height
limit. Let \(A\subset[1,X]\) be a dictionary of distinct integers
admitting representations
\[
                a=sx^k,\qquad
                \operatorname{rad}(s)\le n^f,\qquad k\ge K .
\]
The order, coefficient, and representation may vary with the value.
For a nonzero integer \(h\) write
\[
              r_A(h)=\#\{(a,a')\in A^2:a-a'=h\}.
\]
For every fixed \(\eta>0\), at all sufficiently small fixed \(L\),
\[
                       \sup_{h\ne0} r_A(h)
                         \ll_{f,\alpha,\eta,L,\varepsilon}
                                X^{(2f+\eta)L+\varepsilon}.
                                                               \tag{OG1}
\]
There is no primitivity assumption on a pair in this difference
count.

Enumerate the two coefficients by their radicals and then by their
powers, at total cost \(n^{2f}X^{o_X(1)}\). After fixing them, the
bases obey
\[
                    s x^k-s'y^{k'}=h .
\]
Fix a cutoff \(J=\lceil L^{-2}\rceil+K\). The orders up to \(J\)
give finitely many choices at each fixed budget. For either
chosen pair of orders, the plane equation is absolutely
irreducible. For example, over \(\mathbb C(x)\) the polynomial
in \(y\) is Eisenstein at a simple zero of \(s x^k-h\).
The coefficient-uniform integral plane-curve bound in the
\(X^{1/K}\)-square gives \(X^{K^{-2}+\varepsilon'}\) possibilities,
with constants allowed to depend on the fixed degrees. When
one order exceeds \(J\), enumerate that whole powered value in
\(O(X^{1/J}\log X)\) ways, including the fixed coefficient if
its base is one. The other complete value is determined by the
shift. Since \(K^{-2}+J^{-1}=o(L)\), this proves (OG1) with the
height limit taken first.
[Browning and Verzobio, Lemma 2.2](https://arxiv.org/html/2608.24512v1#S2.SS1)
· [Fixed-support factor count](ABC_EXCEPTIONAL_SET_AND_POWERFUL_TRIPLES.md#a-radical-count-contradiction-to-the-claimed-smooth-number-lower-bound)

## Exact additive graphs and asymmetric dictionaries

Let \(\mathcal E\) be any set of distinct triples \(a+b=c\) with
complete-label dictionaries \(A,B,C\). Write \(D_A=\sup_{h\ne0}r_A(h)\).
Then
\[
 |\mathcal E|\le
 \min\left\{
 |B|+|C|\sqrt{|B|D_A},\
 |C|+|B|\sqrt{|C|D_A}
 \right\}.                                                     \tag{OG2}
\]
Indeed, two distinct edges through \(b\) give a pair of
equal nonzero differences \(a-a'=c-c'\). The two ordered
pairs of complete labels determine \(b\). Thus the count of
these ordered pairs of incident edges is at most
\(\sum_{h\ne0}r_A(h)r_C(h)\le D_A|C|^2\).
By Cauchy–Schwarz the count is at least
\(|\mathcal E|^2/|B|-|\mathcal E|\). Solving this inequality
gives the first bound. Counting pairs through \(c\), with the
summand sign reversed on the other difference, gives the second.

The role of the grid letter can be arbitrary, by using signed
labels on the three parts. Suppose the other two complete
dictionaries have caps \(|B|\le n^\beta\), \(|C|\le n^\chi\).
Combining (OG1) and (OG2) gives normalized exponent at most
\[
 \min\left\{
          \max(\beta,\ \beta/2+\chi+f),\
          \max(\chi,\ \chi/2+\beta+f)
       \right\}+\eta+o_X(1) .                                  \tag{OG3}
\]
Any fixed positive margin below \(1/3\) gives a uniform strict
saving in this sector. This applies in particular to actual
primitive triples with all three balanced radical caps. If
both of the other dictionaries have exponent at most \(v\),
the bound is \(3v/2+f+\eta+o_X(1)\), which is below \(1/3\)
for \(v<2/9-2f/3\). No restriction on the size of the grid
letter's dictionary is needed in (OG3).

## An arbitrary-order defect on a single letter

Fix a prescribed prime pool \(P\) with prime product at most
\(X^{C_0}\), for fixed \(C_0\). For a complete value \(A\)
and an integer \(q\ge2\) let
\[
             D_{q,P}(A)=
               \prod_{\substack{p\notin P\\q\nmid v_p(A)}}p .
\]
Consider rows in any additive family as in (OG2) for which
the complete letter \(A\) admits some \(q\ge2\) with
\(D_{q,P}(A)\le Z\). The order may vary with the row.

**OG6.** The normalized counts (OG3) and (OG4) hold under
this arbitrary-order condition with \(Z=n^f\). More precisely,
at all \(L\), before using a small-budget reduction, the
selected count is at most
\[
 O_{C_0,\varepsilon}\left(
 X^\varepsilon\min\{ |B|+Z|C|\sqrt{|B|},\
                      |C|+Z|B|\sqrt{|C|} \}\right).
                                                               \tag{OG6}
\]
Here \(|B|,|C|\) are the numbers of complete values in
their roles, and the result also applies after permuting
the signed roles. It needs no lower cutoff on the
admissible orders other than two.

To prove it, partition the \(A\)-dictionary by one
witness order for each value. Orders greater than
\(\log_2 X\) have the same defect, so they use one
representative, and the total number of subfamilies is
\(O(1+\log X)\). For each fixed order \(q\), write
\(A=C_A x^q\) with \(C_A\) its canonical \(q\)-free
coefficient. Enumerate the outside support of \(C_A\)
at cost at most \(Z\), then its exact exponents and
its part on \(P\) at uniform height-subpower cost.
At fixed \(h\ne0\), a pair from this dictionary satisfies
\[
                      C_A x^q-C_{A'}y^q=h.
\]
At fixed coefficients, the coefficient-uniform binary
norm count is subpower for every bounded order. At
larger orders one can choose \(x\) in at most \(X^{1/q}\)
ways and then \(y\) is determined. This yields a bound
\(O_{C_0,\varepsilon}(Z^2 X^\varepsilon)\) for each fixed
order and shift, with a uniform constant and no lower
bound on \(X^{1/q}\). Apply (OG2) separately on these
order subfamilies and sum. The allowances are absorbed
in the displayed \(X^\varepsilon\). It is the sum of
graph counts that is uniform across orders. A shift
multiplicity bound on pairs of different orders is
a separate assertion from (OG6).
[Uniform fixed-value binary form count](ABC_RECTANGLE_POWER_CLASS_DEFECT.md#the-fixed-value-bound-with-moving-coefficients)
· [Fixed-support enumeration](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)

### The count across all cells of one packet tuple

Use the notation \(R_F,R_G,R_M,K,R=R_F R_G R_M\) of
[the finite difference-array theorem](ABC_SHARED_VALUE_RESIDUE_BOXES.md#a-finite-difference-array-theorem).
Keep the same arbitrary-order defect on the \(A\)-role,
with a common pool for the tuple. The pool may depend
on the tuple as long as its polynomial height exponent
is fixed. Then
\[
 \frac{|\mathcal E_{\rm grid}|}{RK^2}\ll_{C_0,\varepsilon}
 X^\varepsilon\min\left\{
 \frac{1}{R_M K}+\frac{Z}{\sqrt{R_G K}},\
 \frac{1}{R_G K}+\frac{Z}{\sqrt{R_M K}}
                  \right\}.                                 \tag{OG7}
\]
This follows without fixing one absolute residue cell.
At a fixed order \(q\), two distinct rows sharing \(C\)
give nonzero opposite summand differences. Both complete
\(A\)-values are in one residue on \(G\), and both complete
\(B\)-values in one residue on \(F\). The number of ordered
pairs of \(B\)-values with the latter equality is at most
\(R_F(R_M K)^2\). Their nonzero difference determines
the difference of the \(A\)-pair. Each \(B\)-pair and
\(A\)-pair determines at most one such ordered pair of
rows. For \(E_q\) edges in that order subfamily,
Cauchy--Schwarz and the shift bound in (OG6) give
\[
    E_q^2/(R_F R_G K)-E_q
                 \ll Z^2 X^\varepsilon R_F (R_M K)^2 .
\]
Solving gives the first entry in (OG7). Counting pairs
with a common \(B\) instead gives the second. In that
orientation the \(C\)-pairs are in one residue on \(F\),
with at most \(R_F(R_G K)^2\) such pairs. Their
difference determines the difference in \(A\),
and the complete \(B\)-dictionary has size at most
\(R_F R_M K\). Sum again over the witness orders.

In a simultaneous packet stratum with empirical prices
\(s_F,s_G,s_M\), let \(d=s_F+s_G+s_M<1/6\) and
\(t=1/6-d\). Take \(R_i=n^{s_i+O(\eta+\kappa)}\) and
\(K=n^{t+O(\eta)}\) from the actual residue and ambient
class bounds. With \(Z=n^f\) and a fixed positive margin,
the \(A\)-grid condition gives a strict saving whenever
\[
                         2f<t+\max(s_G,s_M).
                                                               \tag{OG8}
\]
This also bounds the subfamily when the grid is present
on only some rows. Sum over the at most \(n^{d+3\kappa}\)
packet tuples and choose the tolerances and the inner
height allowance after the margin. The diagonal terms
in (OG7) save because \(t\) is bounded below. The same
test works on either other role, with the prices
permuted. The arbitrary-order condition (OG8) can allow
diffuse scaled degrees on all three letters. For
example, taking a common order two imposes no
atomic structure on the scale of \(Le\) as \(L\to0\).
An extraction of an affordable defect from arbitrary
spectra would be an additional arithmetic step.
[Three-role packet sizes and class pruning](ABC_SIMULTANEOUS_COMPOSITE_PACKET_BLOCKS.md#three-marked-roles)

## Consequences near the minimum nonatomic mass

Use the limiting box convention with balanced radical masses
\(1/3\), exponent \(d\), and nonatomic masses \(c_A,c_B,c_C\).
Write \(\Gamma=\sum_W c_W\), \(c_0=1/6\), and
\(\gamma=\Gamma-1/2\). Keep the grid on one prescribed
complete letter across the selection. One can also impose the
outside defect condition of (OG6) on that letter.

At \(d=1/3\) and \(\gamma=0\), the complete-label reduction
gives an actual additive graph of exponent \(1/3-o(1)\)
whose part exponents are at most \(1/6+o(1)\). Applying
(OG3) bounds the resulting graph by \(1/4+f+o(1)\).
This excludes the limiting endpoint on any one-letter
grid sector with \(f<1/12\).

There is a quantitative neighborhood. Use the cyclic
small-anchor reduction (OA6) in the orientation
\(A\to B\to C\to A\), with \(A\) the grid letter. It yields
a graph of complete values with exponent at least
\(10d-6\Gamma\), up to arbitrarily small normalized
losses. The cap exponents on the other two parts are
\[
           \beta=2c_B+c_C-d,\qquad
           \chi=2c_C+c_A-d .
\]
This uses the pair-fibre bound for outside labels and the joint
record freeze of (OA6); in its stated positive range,
the cheap full locks are automatic. At the balanced
endpoint, the pair inequalities give
\(c_0-\gamma\le c_W\le c_0+\gamma\). Therefore
\(\beta\le c_0+3\gamma\) and
\[
              \beta/2+\chi=\Gamma+\tfrac32c_C-\tfrac32d
                            \le \tfrac14+\tfrac52\gamma
              \qquad(d=1/3).
\]
Comparing the first orientation of (OG3) with the graph's
lower exponent \(1/3-6\gamma\) excludes the endpoint
whenever
\[
                               f+\tfrac{17}{2}\gamma<1/12.
                                                               \tag{OG4}
\]
In this range the diagonal term \(c_0+3\gamma\) is below
\(1/4+(5/2)\gamma+f\), and the lower exponent is positive.
All inequalities have a fixed margin in (OG4). The
conclusion rules out balanced limiting profiles in this
sector. A finite condition guaranteeing the complete
dictionary sizes in (OG3) gives the uniform count directly.
[Complete-label reduction near the minimum](ABC_DIFFUSE_FACTOR_GRAPH.md#the-complete-label-reduction-close-to-minimum-mass)
· [Refined anchor-set reduction (OA6)](ABC_DIFFUSE_FACTOR_GRAPH.md#refining-the-labels-by-small-anchor-sets)

## Real-filtered additive graphs at tall heights

Here impose a short-base condition on the grid letter in
addition to its coefficient radical cap. Fix \(h_*>0\)
and \(C_0>0\). Suppose every dictionary value on each
letter is in \([X^{h_*},X]\), with individual radical caps
at most \(n^A\) for a fixed \(A>0\). Each grid letter value
has a presentation \(a=sx^k\) with
\(\operatorname{rad}(s)\le n^f\), \(x\le X^{C_0 L}\), and
\(k\ge\lceil L^{-\alpha}\rceil\). For example the short-base
condition is automatic if \(k\ge1/(C_0L)\); it is stated
separately because it need not follow from the weaker
lower cutoff \(L^{-\alpha}\) used in (OG1).

Write \(\rho,\beta,\chi\) for complete-dictionary exponents
on \(A,B,C\), respectively. For every fixed \(\eta>0\)
there are fixed \(\theta>0\) and \(L_0>0\) such that the
number of distinct triples in the real filter
\[
                              |a+b-c|\le X^\theta
\]
is at most \(X^{(d_{\rm filt}+\eta)L+\varepsilon}\) for
fixed \(0<L<L_0\), uniformly in the dictionaries, where
\[
 d_{\rm filt}=
 \min\left\{
 \max\!\left(\beta,\frac{\rho+\beta+\chi}{2},\frac\beta2+\chi+f\right),
 \max\!\left(\chi,\frac{\rho+\beta+\chi}{2},\frac\chi2+\beta+f\right)
 \right\}.                                                     \tag{OG5}
\]
In particular at the three minimum complete part caps
\(1/6\), this is again at most \(1/4+f\) before the allowance.
As in (OG2), one may permute the letter roles using signed
labels, with all height and radical conditions on their magnitudes.

Use a positive tolerance \(\xi\) much smaller than \(\eta\).
By (SI2), choose \(h_0\in(0,h_*)\) small enough that every
physical interval of length at most \(4X^{h_0}\) meeting one
of the tall dictionaries contains \(n^{\xi+o_X(1)}\)
possible integer labels. To justify the uniform height
margin one can apply (SI2) with centre height threshold
\(h_*/2\) and choose \(h_0\) below half its admissible
width exponent. Now choose \(0<\theta<h_0/2\).
Let \(W=X^\theta\) and replace each complete value \(a\)
by \(\lfloor a/W\rfloor\), and likewise in the other
parts. Every filtered triple then satisfies exact addition
on its bin labels with one of a bounded number of integer
carries. Keep one carry, translating the target bin labels
by it. Each triple of bins lifts to at most
\(n^{3\xi+o_X(1)}\) triples, so the exact bin graph retains
the original exponent up to \(3\xi+o_X(1)\).

For a nonzero integer difference \(k\) of two grid bin
labels, the original values satisfy
\[
                          |a-a'-kW|<W.
\]
If \(|k|W\ge X^{h_0}\), apply (WP4) with fixed target
\(|k|W\) and minimum numerical height exponent \(h_0\).
After coefficient enumeration, the number of such pairs
is at most \(n^{2f+\xi+o_X(1)}\), uniformly in \(k\).
The relevant order in a base greater than one is
\(O(\log X)\). A base equal to one contributes a fixed
term, and the other term has at most one possible base
per order by spacing, since its height is at least
\(X^{h_*}\) and \(h_*-\theta>C_0L\). If both bases are
one, there is at most one pair for the coefficient record.
All orders for nontrivial bases exceed the fixed auxiliary
degree of (WP4) at sufficiently small budgets.

If \(0<|k|W<X^{h_0}\), the two original grid labels
are within \(2X^{h_0}\). Summing their bin-pair counts
over these differences gives at most \(|A|n^{\xi+o_X(1)}\)
by (SI2). For any one nonzero bin difference there are
at most \(|C|\) target bin pairs. Hence the short part
of the difference collision sum is at most
\(|A||C|n^{\xi+o_X(1)}\), and the large part is at most
\(|C|^2 n^{2f+\xi+o_X(1)}\). Use the Cauchy–Schwarz
argument of (OG2) in the exact bin graph. Repeat with
\(B\) and \(C\) exchanged and absorb \(O(\xi)\) into
\(\eta\). This proves (OG5). This real-filter result
counts approximate edges on arbitrary chosen tall
dictionaries, including the exact additive case under
the stronger short-base condition.
[Short-interval radical bound (SI2)](ABC_SHORT_INTERVAL_RADICAL_FILTER.md#arbitrarily-small-normalized-counts-in-shorter-intervals)
· [Both signs of weighted-power near equality (WP4)](ABC_WEIGHTED_POWER_NEAR_EQUALITY.md#both-signs-and-all-omitted-letter-orientations)

## Scope of the gain

The fixed-shift plane-curve step is the same
coefficient-uniform ingredient used in the fine valuation
grid result. The additional counting step places its two
powered values on one letter, using an actual additive
projection to transfer a bound on that letter's
differences. The two-grid estimate permits an arbitrary
third valuation profile and does not impose small
complete dictionaries on the other letters. Here the
other letters can both have arbitrary valuations, while
their dictionary sizes are essential. General diffuse
families need not meet those size hypotheses or have a
one-letter exact valuation grid. The bounds therefore
leave the unrestricted endpoint open.
[Fine valuation grid coefficient estimate](ABC_FINE_VALUATION_DIVISIBILITY.md#exact-sector-bound)
· [Two valuation grids with an arbitrary third letter](ABC_TWO_VALUATION_GRIDS.md)

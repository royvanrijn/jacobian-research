# Additive translation cubes and least-repair classes

A near-quadratic graph of complete additive labels contains a large
common family for any fixed number of prescribed binary shift
allocations between the letters. The base labels may be artificial
integers, but every allocation recovers an original edge with its
original primitive and radical constraints. On the resulting family,
one can rule out small dictionaries for weighted least-repair
classes in each of the three orientations.

## Binary shift allocations between prescribed letters

Let \(\mathcal E\) be a set of \(E\) triples of signed integers
with \(x_1+x_2+x_3=0\) and with at most \(N\) coordinate
values in each slot. Fix an integer \(m\ge1\) and a list of \(m\)
unordered pairs \(\{i_\ell,j_\ell\}\subset\{1,2,3\}\), with two
different slots in each pair. For \(0<\eta<4^{-m}/4\), suppose
\(E\ge N^{2-\eta}\). For all sufficiently large \(N\) there
are positive integers \(h_1,\ldots,h_m\) and an additive graph
\(\mathcal F\) with at most \(N\) labels in each slot and
\[
                         |\mathcal F|\gg_m N^{2-4^m\eta}.
                                                               \tag{TC1}
\]
For every \(\mathbf u\in\mathcal F\) and every assignment
\(k_\ell\in\{i_\ell,j_\ell\}\), the triple
\[
          \mathbf u+\sum_{\ell=1}^m
                       h_\ell(\mathbf e_{k_\ell}-\mathbf e_3)
                                                               \tag{TC2}
\]
belongs to \(\mathcal E\), where \(\mathbf e_i\) is the \(i\)-th
standard basis vector. The \(2^m\) subset sums of the shifts
are all distinct. In particular, for each base triple the
\(2^m\) triples in (TC2) are distinct. If the original coordinates
are bounded in absolute value by \(X\), then the base
coordinates and shifts are \(O_m(X)\).

### One extraction by counting common shifts

For the proof of a single extraction, permute slots and write
the two selected slots as \(a,b\) and the remaining slot as \(c\).
For \(h\ne0\), let \(M_h\) count pairs of original triples of
the form \((a+h,b,c),(a,b+h,c)\). These correspond one-to-one
to the edges \((a,b,c+h)\) of an auxiliary additive graph.
Let \(A_h\) be the set of \(a\) with both \(a,a+h\) in the
original first part and define \(B_h\) in the same way. Then
\(M_h\le |A_h||B_h|\), both sums \(\sum_h |A_h|\) and
\(\sum_h |B_h|\) are at most \(N^2\), and
\[
 \sum_{h\ne0}M_h=\sum_c \deg(c)(\deg(c)-1)
                  \ge E^2/N-E .
\]
The equality counts ordered pairs of distinct original triples
at \(c\); each determines its oriented shift from the first
coordinate difference.

Assume \(E\ge2N\). For the shifts with
\(M_h<E^4/(16N^6)\) use \(M_h\le\sqrt{M_h|A_h||B_h|}\)
and Cauchy--Schwarz. Their sum is at most \(E^2/(4N)\).
The other shifts contribute at least \(E^2/(4N)\), and
each contributes at most \(N^2\). Thus there are at least
\(E^2/(4N^3)\) shifts for which \(M_h\ge E^4/(16N^6)\).
Interchanging the two original triples gives \(M_h=M_{-h}\),
so at least half of these shifts are positive.
Writing \(k\) for the remaining slot, translate the auxiliary
graph by \(h(\mathbf e_3-\mathbf e_k)\). Its labels
still satisfy the zero-sum equation, and adding either
\(h(\mathbf e_i-\mathbf e_3)\) or \(h(\mathbf e_j-\mathbf e_3)\)
recovers an original triple. Each part of the auxiliary graph
was a subset of a translate of its original part. Hence the
new graph still satisfies the same three part caps.

Repeat with the specified pairs. At each step the current
density exponent deficit is multiplied by at most four, with
a constant depending on the step. Before the last step the
number of eligible shifts is at least a positive fixed power
of \(N\). Exclude the at most \(3^{m-1}\) sums of earlier
shifts with coefficients in \(\{-1,0,1\}\) at each choice.
This guarantees distinct subset sums and allows every
iteration. If two assignments in (TC2) differ, they select
different subsets of shifts in at least one slot, so the
resulting triples differ. At each step the shifts are
differences within a subset of a translate of an original
part. The base coordinates differ from the original ones by
sums of at most \(m\) shifts, proving the height assertion.

The common-label pair argument is also used in the
[eight-edge trade count](ABC_BALANCED_LATIN_TRADE_OBSTRUCTION.md#eight-edge-supersaturation-for-actual-additive-values).
Here only counts of pairs are needed, and the recursion gives
arbitrary fixed lists of allocations. For a sequence with
\(E=N^{2-o(1)}\), the conclusion is
\(|\mathcal F|=N^{2-o(1)}\) for every fixed \(m\).

## Large dictionaries of weighted least-repair classes

Suppose the original triples are the signed representations
\((a,b,-c)\) of positive primitive triples \(a+b=c\), all
letters at most \(X\). Any individual radical caps on the
three original letters remain valid on every original triple
in (TC2). Fix the shifts in (TC1) and fix one base coordinate
\(u_k=z\). Call the other coordinates \(u_i=x\) and
\(u_j=-x-z\). Across all allocations in (TC2), the absolute
values of the two varying letters are affine forms in \(x\),
each with slope \(1\) or \(-1\). Within one letter, combine
repeated copies of the same form.

Except for \(O_m(1)\) choices of \(z\), all remaining forms
are pairwise nonproportional. To see this, the signed
expressions are \(x+\alpha\) and \(-x-z+\beta\), with
\(\alpha,\beta\) among fixed offset lists. Two different
forms on the same letter have different roots. For forms
from different letters to have the same root requires
\(z=\alpha+\beta\). Discard these finitely many base
coordinates in each of the three orientations. At most
\(O_m(N)\) base edges are lost, because fixing a coordinate
leaves at most \(N\) edges. The remaining base graph still
has the bound in (TC1).

For an order \(r\ge2\), define the least target repair
\[
        \mu_r(T)=\prod_{0<v_p(T)<r}p^{r-v_p(T)}.
\]
On a nonexceptional fixed fibre, write \(T_1(x),\ldots,T_s(x)\)
for the distinct positive forms from the two varying letters.
Give them arbitrary integer weights \(e_1,\ldots,e_s\), at
least one of which is nonzero modulo \(r\). Set
\[
 g=\min_{i:r\nmid e_i}\gcd(r,e_i).
\]
Every class of \(\prod_i\mu_r(T_i(x))^{e_i}\) modulo
\((\mathbb Q_{>0}^{\times})^r\) has at most
\[
                           O_{m,\varepsilon}(X^{g/r+\varepsilon})
                                                               \tag{TC3}
\]
members of the fibre, uniformly in the shifts, the fibre,
the order, the weights and the class.

To prove (TC3), take one member \(x_0\) of the class when
there is one. Let \(S\) contain the prime supports of all
\(T_i(x_0)\) and of all nonzero pair determinants of the
affine forms. The product of the primes of \(S\) has
polynomial height with exponent depending only on \(m\).
Outside \(S\) at most one \(T_i(x)\) can be divisible by
any given prime, by subtracting the two affine forms after
multiplying by the alternate slopes. Choose \(j\) with
\(\gcd(r,e_j)=g\) and a prime \(p\notin S\) dividing
\(T_j(x)\). If \(v_p(T_j(x))<r\), equality of the classes
at \(x,x_0\) gives
\[
       r\mid e_j\bigl(r-v_p(T_j(x))\bigr),
       \qquad\text{so}\qquad r/g\mid v_p(T_j(x)).
\]
If \(v_p(T_j(x))\ge r\), the same lower bound \(v_p\ge r/g\)
holds without a congruence. Therefore the product of the
primes of \(T_j(x)\) outside \(S\) is at most \(X^{g/r}\).
Enumerate this squarefree product. For each product there
are \(O_{m,\varepsilon}(X^\varepsilon)\) possibilities
for \(T_j(x)\), by uniform fixed-support counting at
polynomial height. Each one determines \(x\), since the form
has nonzero slope. This proves the bound. No height cap
is needed for the prescribed class representative or the
weights, as the proof uses only the support of one original
realization. Weighted products may include the fixed letter, since its
shifted values are constant on that fibre.
[Uniform fixed-support count](ABC_EXCEPTIONAL_SET_AND_POWERFUL_TRIPLES.md#a-radical-count-contradiction-to-the-claimed-smooth-number-lower-bound)

Consequently, if at most \(J\) weighted repair classes occur
on each fixed fibre in a given orientation, then the retained
base graph has at most
\[
                       O_{m,\varepsilon}(N J X^{g/r+\varepsilon})
                                                               \tag{TC4}
\]
edges. The weights and order can depend on the fixed
coordinate if \(g/r\) in (TC4) is a uniform upper bound for
their quotients. At \(N=X^{(c+o(1))L}\) for small budgets
\(L\) and \(|\mathcal F|=X^{(2c-o(1))L}\), this rules out
an exponent below \(c-g/(rL)\) for \(J\), up to the
extraction losses. In particular, if there is a unit
weight modulo \(r\) and \(rL\to\infty\), every such
dictionary has exponent at least \(c-o(1)\) in the
near-full fibres. The statement is interpreted with the
height limit first and the inner \(\varepsilon\) chosen
smaller than the desired margin times \(L\). The choice
of shifts may depend on height. This applies to the
complete-label minimum-mass reduction with \(c=1/6\);
with nonzero extraction losses one uses (TC1) and (TC4)
directly.
[Complete-label reduction and its losses](ABC_DIFFUSE_FACTOR_GRAPH.md#the-complete-label-reduction-close-to-minimum-mass)
· [Homogeneous repair-star counts](ABC_ADAPTIVE_HOMOGENEOUS_POWER_STARS.md#separable-least-target-repairs-without-affine-realization)

This conclusion constrains repair classes on simultaneously
translated original triples. An upper bound for any one
of the weighted class dictionaries sufficient to contradict
(TC4) remains open in the broad diffuse sector.

## Primitive graphs without three placements of the same shift

For every sufficiently large integer \(M\), take a
three-term-progression-free set \(S\subset[1,M]\) of size
\(M\exp(-O(\sqrt{\log M}))\) as in Behrend's construction.
Keep the positive triples
\[
        (a,b,c)=(a,a+d,2a+d),\qquad
        M\le a<2M,\quad d\in S,\quad (a,d)=1.
\]
All are primitive, and every part has \(O(M)\) labels.
For a fixed \(d\le M\) the interval for \(a\) contains at
least \(\lfloor M/d\rfloor\varphi(d)\) allowed values.
Since \(\varphi(d)/d\ge d^{-o(1)}\) uniformly as \(d\)
grows, the graph has \(M^{2-o(1)}\) edges.
[Behrend, *On sets of integers which contain no three terms in arithmetical progression*](https://pmc.ncbi.nlm.nih.gov/articles/PMC1078964/)

In signed coordinates suppose one base triple allowed all
three placements of some common \(h\ne0\) in the sense of
(TC2). Writing the base as \((a,b,-a-b)\), this would supply
the original triples
\[
     (a+h,b,a+b+h),\qquad
     (a,b+h,a+b+h),\qquad
     (a,b,a+b).
\]
Their \(b-a\) differences are \(b-a-h,b-a+h,b-a\),
three distinct elements of \(S\) in arithmetic progression.
This is impossible. Thus primitivity and numerical additivity
at near-quadratic density suffice for arbitrary binary
allocations but do not force ternary allocations of a common
nonzero shift. The radical budgets of this example are unrestricted.

## Joint radical size on the translated forms

Assume now that the original graph is in the height shell
\(X/2<c\le X\), with individual radical caps
\(n^{1/3+o(1)}\), where \(n=X^L\). Suppose the translated
base graph after removing exceptional fibres has part caps
\(n^{1/6+o(1)}\) and edge exponent \(1/3-o(1)\).
The number of shifts is fixed before the height limit.
On each fixed base-coordinate fibre, let \(P\) be any
prescribed positive integer at most \(X^B\) for a fixed
\(B\). Its choice can depend on the fibre and on the shifts.
For every fixed \(\xi>0\) the number of base edges on
which at least one varying form \(T_j(x)\) has normalized
radical mass at least \(\xi\) on \(P\) is at most
\[
                        n^{1/3-\xi/2+o(1)}                    \tag{TC5}
\]
over all fibres of that orientation. The \(o(1)\) here
allows arbitrary coefficient tolerances in the height-first
convention.

Indeed each varying form is an original letter on one of
the translated triples. Its opposite letter corresponding
to the fixed base coordinate is fixed on the fibre. Remove
the prescribed prime pool from the varying form's radical
cap and apply (FL1) to that original triple. It counts at
most \(n^{1/6-\xi/2+o(1)}\) base points in the fibre, since
the varying affine form is injective. Sum over the bounded
number of forms and the fixed-base coordinate dictionary.
The argument applies simultaneously in all three
orientations.
[Pointwise count with prescribed prime pools (FL1)](ABC_FIXED_LETTER_RADICAL_CAPS.md#a-pointwise-count-with-prescribed-prime-pools)

In particular take \(P\) to include every nonzero pair
determinant and any prescribed polynomial-height pool.
Outside it the distinct affine forms are pairwise coprime.
If each original letter also has radical at least
\(n^{1/3-o(1)}\), the product of the \(s\) distinct varying
forms has outside radical at least \(n^{s/3-o(1)}\) after
discarding smaller-exponent families at any fixed tolerance.
Thus the extraction allows joint radical growth on any
fixed number of translations, even when the shifts and
their common prime pools depend on height.

The same estimate sharpens (TC3) on actual balanced
fibres. With the class-specific pool from its proof, the
outside radical of \(T_j\), for an active weight with
\(\gcd(e_j,r)=g\), is at most \(X^{g/r}\) as well as its
original radical cap. Applying (FL1)
to an original branch containing that form gives class
size at most \(X^{g/(2r)+\eta L+\varepsilon}\), uniformly
in the order and the weights. This uses all the original
radical caps and the fixed original third letter in that
branch.

# Relative coefficient counts on genuine common-query stars

One query with several distinct anchors gives several independently
supported affine values after the anchor-difference primes are fixed.
A comparison of centered coefficients at two queries then tests the
remaining radicals in both directions and on every copy of the original
fibre. The resulting count uses the separate original caps and the
nonatomic part of each conditional degree distribution.

## Complete lifts and conditional coefficient data

Work with \(n=X^L\), positive primitive rows of height
\(X/2<C\le X\), and separate original radical caps bounded in
base \(n\). Use the height limit at each fixed budget before
the limit of decreasing budgets. Fix a list of \(s\ge1\)
distinct complete anchors \(b_1,\ldots,b_s\) in one original
role. A common query is a complete original value \(V\)
that participates in a genuine row with each anchor. The
other complete values are \(W_i=\epsilon V+\epsilon' b_i\)
with the signs fixed by the three original roles. Include
\(W_0=V\). Every \(W_i\) is a free letter in a fixed-anchor
original fibre, including \(W_0\), which can be counted
using any of the anchors.

Prescribe one polynomial-height prime pool for each \(W_i\).
The computed-role pools contain all primes of the
nonzero anchor differences \(b_i-b_j\). They may also
contain actual full-divisor records and previously
recorded support. Fix the complete
supported parts at height-subpower cost as usual, and use
only the outside values in what follows. Outside these pools
the \(W_i\)'s have pairwise disjoint supports at each query:
each \(W_i\) is coprime to \(V\), and any prime common to
two computed values divides the corresponding anchor
difference. The pools may vary with height and with an
already fixed joint record. Their total height exponents
are bounded uniformly for fixed \(s\).

Write \(\ell_i\) for the residual raw radical cap on \(W_i\).
In a stabilized sequence of such stars, let \(m_i\) be the
nonatomic mass of the residual limiting measure on that
letter. Use a chosen nonempty subset \(\mathcal I\) of the
indices with \(m_i>0\). If any free role has \(m_i=0\),
FT1 already excludes a positive query exponent in that
profile. For each \(i\in\mathcal I\) retain actual degree
slots whose measures converge uniformly to the full
nonatomic part of the residual measure for \(i\).
The degree selectors can be taken at continuity
intervals as in APC6. The complement of the retained
slots then converges to the atomic part. The selectors,
the initial record and the pools are fixed at each
retained budget before using the inner height limit.

For \(K\asymp L^{-\alpha_0}\), \(1/2<\alpha_0<1\), use one
common \(k\in[K,2K)\). On each selected slot of exact
degree \(e\) use its centered remainder \(r_e\) modulo
\(k\), with \(|r_e|\le k/2\). Fix any exponent signs
\(\sigma_i\in\{-1,1\}\) for the roles and define
\[
             Q(V)=\prod_{i\in\mathcal I}
                       \prod_{p\ {\rm retained\ on}\ W_i}
                                              p^{\sigma_i r_e}.
\]
Pairwise disjoint supports on one query make this an
exact rational encoding of its retained signed
remainders. Call a query good at this order when the
raw mass distribution of \(|r_e|/k\) on each
selected \(W_i\) is uniformly close to density \(2m_i\)
on \([0,1/2]\). Require also that the sums of positive
and negative remainders on that role, in units of
\(kL\log X\), are each \(m_i/8+o(1)\). On a stabilized
profile the one-slot averaging and second-moment
argument of APC6 show goodness for all but a
vanishing proportion of rows at all but a vanishing
proportion of the orders. For a given order the
bounds below apply to its good queries. If a cover
uses several good orders, their number is finite
at each fixed budget.
[The original fixed-anchor lifts (LFR1)](ABC_BALANCED_ORIGINAL_FIBRE_LIFT.md#counting-a-common-neighbour-in-its-original-fibre)
· [Uniform remainder distributions (APC6)](ABC_ADAPTIVE_POWER_COEFFICIENT_PRICES.md#fine-grid-distances-at-balanced-saturation)

## Counting both directions of a short relative pair

Put
\[
       d_0=\min_{i\in\mathcal I}\{m_i,\ell_i/2\},\qquad
       A_i(e)=\max\{e,\,2e-(\ell_i-m_i)\},\qquad
       \mathcal C(e)=\frac14\sum_{i\in\mathcal I}
                                      \frac{A_i(e)^2}{m_i}
                 \quad (0\le e\le d_0).
\]
For a real \(c\ge0\), let \(\beta(c)\) be the maximum
\(e\in[0,d_0]\) with \(\mathcal C(e)\le c\). The set is
nonempty because \(\mathcal C(0)=0\). The initial
FL1 and FT1 bounds give query exponent at most \(d_0\)
in the stated profile.

**ACQ1.** In any good query subfamily \(\mathcal S\) in
this profile, the count of ordered pairs \((V,V')\)
with \(H(Q(V)/Q(V'))\le X^{c kL}\) is
\[
        O\bigl(|\mathcal S|\,n^{\beta(c)+\eta}X^\varepsilon\bigr)
\]
for every prescribed fixed \(\eta,\varepsilon>0\) and
sufficiently late fixed budgets. As usual the
constants on the stabilized profiles may depend on
the fixed cap data and tolerances. The bound uses
fixed margins in the normalized exponents, so the
inner \(X^\varepsilon\) allowances can be taken
inside any required power margin at each fixed
budget. In particular, for \(h=|\mathcal I|\) and
a common raw cap \(R\) on the selected roles, it
suffices to use
\[
                      \beta(c)\le 2\sqrt{\frac{R c}{h}} .
\]

To prove the count, fix an actual good reference
\(V'\). Let \(E_i(V\mid V')\) be the raw mass of the
retained primes on \(W_i(V)\) outside
\(\operatorname{supp}(Q(V'))\), including those of
zero centered remainder. For any fixed \(v\ge0\),
the condition \(E_i(V\mid V')\le v\) bounds the
normalized count after fixing the reference by
\[
                       \min\left\{v,\frac{\ell_i-m_i+v}{2}\right\}
                  \quad\text{up to a fixed positive allowance.}
                 \tag{ACQ2}
\]
For the second term use the original fixed-anchor
lift for \(W_i\). Adjoin the prime support of the
reference coefficient to the previously
prescribed pool on this role. Its outside raw
mass is at most \(\ell_i-m_i+v+o(1)\). Apply FL1
with the other original free role and fixed
anchor left in the count.

For the first term apply FT1 to that same original
lift with this enlarged prescribed pool. The
outside retained slots have total mass at most
\(v\). The outside unretained slots have limiting
measure dominated by the atomic part of the
original residual profile, since they are
submeasures of the unretained slots before
removing the reference support. They contribute
no nonatomic mass. More explicitly, a sequence
of counterexamples with a fixed positive count
margin can choose the pool and reference at each
inner height. Make the permitted anatomic
refinement on the resulting fixed-anchor
fibres and take the outer weak limits. The
retained outside measure has total mass at most
\(v\), and domination of the unretained weak
limit follows by testing against nonnegative
continuous functions. FT1 then contradicts
a normalized exponent bounded above \(v\)
by the fixed margin. This also gives uniformity
in polynomial-height reference supports along
the stated stabilized profiles. All supported
parts are counted at subpower cost, and all
original radical caps are imposed throughout.
[Outside raw caps (FL1)](ABC_FIXED_LETTER_RADICAL_CAPS.md#a-pointwise-count-with-prescribed-prime-pools)
· [Outside nonatomic caps (FT1)](ABC_FIXED_FIBRE_TRIANGLE_RADICALS.md#nonatomic-mass-outside-prescribed-prime-pools)

Now choose \(e>\beta(c)\) in \([0,d_0]\) with a
fixed margin above \(\beta(c)\). For a sufficiently
small fixed \(\delta_0>0\), for each direction
and each role the pairs with
\(E_i(V\mid V')\le A_i(e)-\delta_0\)
cost at most \(|\mathcal S| n^{e-\delta_1}
X^\varepsilon\) for a fixed \(\delta_1>0\).
Indeed at least one branch of (ACQ2) at that
cutoff is strictly below \(e\).
In the remaining pairs the exclusive support
weights, in units of \(kL\), obey
\[
 2c \ \ge\ \sum_{i\in\mathcal I}
             \frac{ E_i(V\mid V')^2+E_i(V'\mid V)^2}{4m_i}
                           -o(1).
 \tag{ACQ3}
\]
In each direction the centered magnitude mass
on a set of size \(z\) is at least
\(z^2/(4m_i)-o(1)\), by integrating its lower
tail against the uniform CDF bound. Numerator
and denominator absolute weights on primes
exclusive in the two respective directions
cannot cancel. The centered contributions
of zero magnitude cause no difficulty in
this integration. The right side of (ACQ3)
is strictly above \(2c\) when all masses
are at least \(A_i(e)-\delta_0\), because
\(\mathcal C(e)>c\). Thus every close pair
belongs to the exceptional classes counted
by (ACQ2). Take \(e<\beta(c)+\eta\)
by a fixed margin inside the allowance.
When there is no such \(e\le d_0\) the initial
row bound of exponent \(d_0\) already suffices.
This proves ACQ1. The profile-free form follows
from \(A_i(e)\ge e\) and \(m_i\le R\).

For example at the balanced half-cap endpoint
\(e=1/6\) and \(\ell_i\le1/3\), a star of that
exponent must have \(\ell_i=1/3\) and \(m_i\ge1/6\).
It follows that \(A_i(e)=m_i\). For every fixed
gap below \(\sum_i m_i/4\), a power-small fraction
of good pairs fall at or below that radius.
The elementary height bound \(H(Q/Q')\le
H(Q)H(Q')\) gives the matching upper radius
\(\sum_i m_i/4+o(1)\) on good pairs.

## Relative covers and common-query star counts

A cover of \(\mathcal S\) by \(n^{b+o_X(1)}\)
sets of pairwise relative diameter at most
\(X^{c kL}\) implies a query exponent at most
\(b+\beta(c)\), up to any fixed positive
counting allowance. Assign each query to
one set, then the number of ordered close
pairs is at least
\(|\mathcal S|^2/n^{b+o_X(1)}\) by
Cauchy–Schwarz. Apply ACQ1. A ball cover
using radius exponent \(c/2\) also suffices.
A conditional catalogue of full-factor
records of exponent \(f\), paid before
the cover, gives \(f+b+\beta(c)\). The
uniform expression \(2\sqrt{Rc/h}\)
may be used when the profiles differ
among the retained record classes.
The cover and selected coefficient
use the profiles *after* fixing the
actual records.

In an ordinary balanced residue cell let
\(u=1/3-\alpha>0\) be the local complete
dictionary exponent on each side and suppose
the graph reaches edge exponent \(2u\).
For fixed \(s\), (CQ2) gives distinct
\(s\) anchors with common-query exponent
approaching \(u\). It is therefore sufficient
to give the above conditional covers on
the required full-size stars, with numerical
margin
\[
                     f+b+\beta(c)<u.
                      \tag{ACQ4}
\]
For a near-full edge count hypothesis,
the finite star extraction loses at most
a fixed multiple of its edge deficiency
and the dictionary allowances. Thus a
fixed margin in (ACQ4) gives a fixed
edge exponent saving on that sector.
At the balanced individual cap \(R=1/3\),
use all \(s+1\) roles to get the explicit
sufficient condition
\(f+b+2\sqrt{c/(3(s+1))}<u\).
The number of anchors can be chosen before
the budget limit, so this count can use
independently supported observations on
many copies of an original computed role.
It works for every fixed positive remaining
multiplicity \(u-f-b\). Producing such
covers, or directly proving a large close-pair
density on these genuine stars, remains an
extraction problem. The occurrence endpoint
\(u=0\) still needs a separate argument.
[Star extraction (CQ2)](ABC_BALANCED_CELL_FULL_FACTOR_QUANTILES.md#a-uniform-count-for-packets-on-repeated-copies-of-one-letter)
· [Balanced residue cells (RB4)](ABC_BALANCED_DICTIONARY_RESIDUE_CELLS.md#the-actual-residue-cell-sizes)

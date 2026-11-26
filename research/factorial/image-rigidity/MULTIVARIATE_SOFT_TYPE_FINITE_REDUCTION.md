# Finite reductions by distinct soft row types and tail gap budget

Delete equal rows from an ordinary binomial profile. Write
\(d=|b-a|>0\), \(M=\max(a,b)\) for an unequal row \((a,b)\).
Take a computable absolute relative-gap tolerance
\(\eta_0\in(0,1/2]\) from the
[diffuse all-order theorem](MULTIVARIATE_DIFFUSE_ALL_ORDER_RESULTANT_TRANSFER.md).
Call a row hard if \(d/M>\eta_0\), and soft otherwise.
A row type means one ordered pair \((a,b)\), with any
positive multiplicity.

**Theorem.** There is a computable absolute height bound on
all collision profiles with at most one distinct soft type.
In particular this bounds every collision profile with at
most two distinct unequal row types, with arbitrary gaps,
orientations, multiplicities and dimension. All-degree SFC(2)
for this class is decidable by a finite set of exact resultants.

More generally, fix a nonnegative integer \(K\). Impose the
condition that the smaller of the total positive and total
negative absolute gaps among soft rows is at most \(K\).
There is a finite effective reduction of this entire class:
it returns a list of bounded exact collision profiles and
a finite list of families with fixed ordinary cores and
fixed gaps and counts on the remaining height variables.
At each such unresolved family, for its specified adjacent
pair, the twisted core comparators have a common conjugate
pair of odd multiplicity in both. The corresponding real
algebraic locus also has a branch with all remaining
heights tending to positive infinity. These two conditions
are necessary for the unresolved families, and the reduction
returns a finite decision whenever no such family survives.
In particular \(K=0\) treats arbitrary cooriented soft tails.

The first assertion uses the one-unbounded-type height bound,
with distinct types counted after repetitions are collected.
The second combines the same hard-row cutoff with a bounded
absolute tail budget, then applies Hermite splitting and an
exact real-algebraic test to a finite tree of height slices.

## Explicit hard and soft budgets

Let \(J\) be a computable absolute angular cutoff and put
\[
       D=146J-73,\qquad
       W=\lceil 12/\eta_0\rceil,\qquad
       B=\lceil 12/\eta_0^2\rceil .
\]
By the [curvature cutoff](MULTIVARIATE_ABSOLUTE_CURVATURE_12.md)
and the [signed degree cutoff](MULTIVARIATE_ALL_DEGREE_SIGNED_IMBALANCE.md)
every counterexample satisfies \(S<12\) and
\(|\Delta|<D\). It has at least one hard row by the absolute
diffuse theorem. Each hard row contributes \(d^2/M>\eta_0 d\)
to \(S\), so, with repetitions counted,
\[
                     \sum_{\rm hard}d\le W.
\]
Also \(M<(d^2/M)/\eta_0^2<12/\eta_0^2\) for every hard row.
Thus the multiset of hard rows lies in a finite, explicitly
enumerable set. Its absolute gap budget is at most \(W\).
In particular the signed gap of the soft rows satisfies
\[
                    |\Delta_{\rm soft}|\le D+W .                \tag{1}
\]

For a profile with at most one distinct soft type, every
row of height above \(B\) is of that same type. The
[one-unbounded-type theorem](MULTIVARIATE_ONE_UNBOUNDED_TYPE_FINITE_DECISION.md)
therefore applies, with its stated computable height
bound \(H_B\). If the entire profile has at most two distinct
unequal types, at least one is hard at a collision, so it
belongs to this case. After the height bound, \(S<12\)
bounds the number of unequal rows by \(12H_B\), since
the curvature of an unequal integer row of height at most
\(H_B\) is at least \(1/H_B\). This yields the finite
decision, using terminals below \(J\).

For the second assertion, let \(P,Q\) be the positive and
negative absolute soft gap sums. By (1),
\[
             P+Q=|P-Q|+2\min(P,Q)\le D+W+2K .                   \tag{2}
\]
So all remaining soft rows can be labeled individually
and there are finitely many lists of their signed gaps,
including all repetitions. We allow their smaller exponents
to vary independently. This enlarges any family in which
some rows are repeated at the same height. The original
soft constraints \(d/M\le\eta_0\) are simple computable
lower bounds on the individual smaller exponents.

## An effective height-slicing tree

For each enumerated hard core and signed soft gap list fix a
terminal \(N=n+1<J\). The terminal two case is excluded by
real positivity. At a node of the tree some heights have
been fixed and their rows absorbed into the ordinary core
\(C\). The remaining heights are free nonnegative integers
subject to their original lower bounds. Let \(h\) be the
signed gap sum of these remaining rows. The zero-strength
comparators are
\[
                         Q_j^C(j^h z),\qquad j=n,N .
\]
The constant twisting factors use the *remaining* signed
gap, so a row newly absorbed at a deeper node changes
\(h\). Compute \(g\), the number of distinct common roots
of odd multiplicity in both comparators. This is the
degree of the gcd of their odd squarefree parts, computed
exactly over \(\mathbb Q\). At a node with \(g\le1\), the
[fixed-count Hermite bound](MULTIVARIATE_FIXED_GAP_HERMITE_HEIGHT_BOUND.md)
provides a computable cutoff \(L\) so that a collision
requires some remaining height below \(L\). The allowed
bounded values of that height give finitely many child
nodes. A common factor at an ordinary profile must have
degree at least two by real even-order positivity, so this
use of the Hermite bound needs no assumption of simple
comparator roots.

There is a stronger finite test if \(g\ge2\). At positive
smaller exponents \(A_i\) use variables \(t_i=1/A_i\).
Strip factorial endpoints and make the common scale
change of the Hermite theorem. The resulting consecutive
polynomials in \(z\) have rational polynomial coefficients
in the \(t_i\), because the gaps and counts are fixed.
Their leading and constant coefficients are nonzero
for positive heights. Impose the Sylvester minors for
gcd degree at least two. [Real quantifier elimination](https://arxiv.org/pdf/1609.02879)
decides whether this real algebraic locus has points with
all \(0<t_i<\delta\) for every \(\delta>0\). This is a
first-order formula over the rational numbers. If the
answer is negative, a rational positive \(\delta\) for
which the box is empty can be found by exact search and
quantifier elimination. Then some height is bounded,
and we take the corresponding finite slices. Heights zero
can always be included among the bounded slices.

At nodes with just one free height, the ordinary
[affine-ray theorem](MULTIVARIATE_AFFINE_BASE_RAY_ALL_DIRECTIONS.md)
gives an exact bound without a comparator test. At nodes
with no free height one computes the exact integer
resultant. Every splitting step fixes a height and hence
reduces the depth still available. The tree is finite
since the original number of free rows is bounded by (2).
The only leaves where the reduction can leave unbounded
heights have \(g\ge2\) and nonempty real escape as tested
above. Since the comparators have real coefficients and
their even parent has no nonzero real roots for a core
with distinct columns, a positive \(g\) in such a core
already gives a common conjugate pair. If an empty or
equal-column core occurs in the enumeration, the
comparators are proportional to \((1+j^h z)^j\), so
there is no common odd-multiplicity pair and the
Hermite case applies. This proves the asserted necessary
conditions and the finite certificate when there are no
unresolved leaves.

The criterion uses real escape only as a necessary condition.
A real branch at infinity may contain no integer height
profiles. The integer-point problem on those remaining
branches is the step left open by this reduction.

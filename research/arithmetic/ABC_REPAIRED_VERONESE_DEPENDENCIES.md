# Simultaneous power release and interpolation after least repair

An interpolation determinant can vanish through proper
subsums even when the underlying evaluation vectors are
minimally dependent. These cancellations can be handled
by retaining every maximal minor and separating the
prime packets according to where they can be detected.
The argument works for the whole fixed-degree projective
linear system, with no irreducibility assumption on the
curve through the points.

## Power release for a simultaneous system

Fix \(m\) root variables, a number of equations, and their
finite exponent arrays. Each equation is a sum of
nonzero weighted monomials in the \(m\) variables.
Count positive pairwise coprime integer tuples \(x_j\le H\)
that satisfy all the equations without a proper vanishing
subsum in any of them. Suppose that for every variable
\(j\), the variable occurs with two different exponents in at
least one equation, and every positive exponent spread
is at least \(E\). For \(m\ge2\)
the count is
\[
             O\!\left(H^{C E^{-1/(2(m-1))}+\varepsilon}\right).
                                                               \tag{SIH1}
\]
Here \(C\) depends on the number of variables and total
number of monomials. The implied constant may also depend
on the fixed exponent arrays. Both bounds are uniform in
the coefficients. A two-term equation is allowed.

For each equation, divide out every variable whose exponent
is constant on its monomials. The projection onto the
remaining variables satisfies a nondegenerate equation
and every remaining variable has exponent spread at
least \(E\).
Apply (EB1) separately to the visible variables of each
equation, with their own number of variables. If there is
only one such variable, use the ordinary univariate
root bound on the nonzero polynomial. If there are two
terms, fixing their monomial ratio reduces to assigning
the primes of a polynomial-height rational number to
pairwise coprime variables. Each prime has at most \(m\)
possible assignments, and once assigned its valuation in
that variable is determined. This costs \(H^\varepsilon\)
for fixed exponent arrays. An equation without visible
variables costs one empty record. The union of the
visible sets contains all variables. Multiplying the
finitely many list-size bounds, with the allowances
split among the equations, proves (SIH1).
[Coefficient-uniform nondegenerate power release (EB1)](ABC_NONDEGENERATE_POWER_COLLISIONS.md)

## Minimally dependent tuples of repaired points

Use \(n=X^L\), with \(L>0\) fixed before the height
limit. Start with any set of positive primitive rows
\(a+b=c\), \(X/2<c\le X\), and
\(\operatorname{rad}(abc)\le n\). For a fixed \(K>0\)
use \(M=\lceil K/L\rceil\) and the three independent
least repairs
\[
             \mathbf h(a,b,c)=(H_M(a),H_M(b),H_M(c)).
\]
Fix \(d\ge1\) and \(2\le j\le \binom{d+2}{2}\).
For each repaired point use the evaluation vector of all
degree-\(d\) monomials in its three coordinates. A tuple
of \(j\) rows is minimally dependent if these evaluation
vectors are linearly dependent but every proper subset
is independent. Assume that for each two rows, the
ratios of their repaired coordinates in any two prescribed
roles differ. Also assume that the radical of the product
of primes used in at least two of the \(j\) original rows
is at most \(n^\rho\). For every \(\eta>0\) and all sufficiently
small fixed \(L>0\), the number of such tuples is
\[
                      O(n^{\rho+3j/K+\eta}X^\varepsilon).
                                                               \tag{SIH2}
\]
The threshold on \(L\) can depend on \(d,j,K,\eta\).
There is no lower bound on the share of the radical
carried by any letter or degree.

For a complete letter \(U\), let
\(W(U)=\prod_{v_p(U)>M}p^{v_p(U)-M}\).
Then \(H_M(U)=W(U)\operatorname{rad}(U)^M\) and
\(\operatorname{rad}(W(U))\le X^{1/(M+1)}\).
Record all \(3j\) tails at cost \(n^{3j/K}X^\varepsilon\).
Record the shared radical at cost \(n^\rho\) and allocate
its primes to their occupied letters and rows at
subpower cost. This gives a fixed integer \(C_{i,U}\)
multiplying \(x_{i,U}^M\) in each repaired coordinate,
with \(x_{i,U}\) the remaining squarefree private
radical product. The \(3j\) variables are pairwise
coprime and at most \(n\), including the possibility
of units. Their disjointness within each row uses
the coprimality of all three primitive letters.
At fixed records the support of each original letter
is determined once these variables are known. The
fixed-support estimate counts exact original rows at
subpower cost.
[Canonical target repair formula (MPR)](POWERFUL_TRIPLES_COMMON_DIVISOR_BOUND.md#a-reverse-transfer-by-repairing-only-the-target)

Use every \(j\)-row maximal minor of the evaluation
matrix, and partition the signed terms of each minor
into minimal vanishing subsums at the sample. There
are only finitely many choices depending on \(j,d\).
The resulting equations are all nondegenerate by
definition of the blocks. Their exponent spreads in
a private variable are integer multiples of \(M\).
To apply (SIH1) it remains to check that each private
variable has a witness among those blocks.

Suppose a variable \(x_{i,U}\) has no witness. In each
block its exponent is constant among the terms, so
changing that variable alone preserves every vanishing
minor. The other \(j-1\) evaluation vectors are
independent. It follows that the new evaluation
vector at \(i\) always remains in their span. In
projective coordinates the changing point traces
infinitely many points of the line through
\(\mathbf h_i\) and the coordinate vertex of role \(U\).
Every homogeneous degree-\(d\) polynomial vanishing
on the other \(j-1\) points must therefore contain
this line.

There is a nonzero polynomial through the other
\(j-1\) points because there are \(\binom{d+2}{2}\)
coefficients. None of these other points lies on the
line: belonging to it is exactly equality of the
ratio of the two coordinates different from \(U\).
If the chosen polynomial has the line as a factor,
remove its full multiplicity and replace it by the
same power of a different linear form. The resulting
homogeneous polynomial has the same degree, still
vanishes at every other point, and no longer contains
the line. This contradiction supplies a witness for
every private variable. The proof does not require the
other points to lie on a unique degree-\(d\) curve, or
any polynomial through them to be irreducible.

Apply (SIH1) to these blocks, with \(E=M\) and \(H=n\).
At fixed \(L\) all exponents are fixed, and the extra
normalized exponent tends to zero as \(L\to0\).
Using the records and subpower lift above proves (SIH2).

## Prime packets supported on several points

The argument prices fewer shared primes for small circuits.
Let \(s_{d,j}\) be the largest \(s\le\min(j,d+1)\) for
which
\[
                 j-s<\binom{d-s+3}{2}.                     \tag{SIH4}
\]
There is always at least \(s=1\), and every positive
integer \(s'\le s_{d,j}\) satisfies the same inequality
with \(s'\). In (SIH2) it suffices to cap by \(n^\rho\)
the product of primes occurring on more than \(s_{d,j}\)
of the \(j\) original rows. Record these primes and their
role allocations. Split the others into packets by
their exact sets of occupied row-and-role positions.
All variables remain pairwise coprime and no packet
can meet a row in two different positions. There are
at most \(4^j-1\) nonempty patterns.

To check visibility, take a packet meeting \(s\) rows
with \(1\le s\le s_{d,j}\). Varying its variable moves
each affected point on its own coordinate line.
These lines are distinct. For the same variable role
this follows from ratio separation. For different roles
a common line would contain two coordinate vertices,
hence no point with all coordinates positive. None of
the \(j-s\) fixed points is on any of these lines, again
by ratio separation. Pick one moving row \(i\). There
is a nonzero form of degree \(d-s+1\) vanishing on the
\(j-s\) fixed points by (SIH4). Remove any factors of
its \(i\)-th coordinate line and replace them by a
different linear factor. Multiply the result by the
\(s-1\) other moving line equations. This gives a
degree-\(d\) form vanishing on every other point in the
deformation, but not identically on the \(i\)-th line.

If the packet were invisible to every block in every
maximal minor, all the minors would vanish along the
deformation. The other \(j-1\) evaluations stay
independent outside a finite parameter set, since
they were independent at the original tuple.
The constructed polynomial would then vanish on
infinitely many positions of row \(i\), a contradiction.
So all the unrecorded patterns have witnesses and
the same simultaneous power release proves the stated
strengthening.

A shared packet can be invisible when it meets too
many rows. Take the six actual primitive rows
\((1,b,b+1)\) with
\[
                 b=14,\ 20,\ 35,\ 125,\ 2400,\ 4374.
\]
Their products \(b(b+1)\) all have radical \(210\),
with \(\operatorname{rad}(b)\) respectively
\(14,10,35,5,30,6\). For every \(M\ge7\), their repaired
points are \((1,R^M,(210/R)^M)\) and lie on the conic
\(YZ=210^M X^2\). All pairwise coordinate ratios are
distinct. These six points form a degree-two circuit:
the restrictions of quadratic monomials to the conic
with \(X=1\) span \(Y^{-2},Y^{-1},1,Y,Y^2\),
so evaluation on any five of its nonzero distinct
\(Y\)-values has rank five. The prime \(2\) is on the
second coordinates of some rows and the third of the
others. Varying this shared packet replaces \(210^M\)
by a common varying constant in the conic equation.
The determinant remains zero identically in the packet
variable. Grouping its signed terms by the exponent
in that variable yields a partition into minimal blocks
for which the packet is invisible. This illustrates
why shared packets cannot be treated as private ones.

## Generic independence at the small-budget count scale

For any fixed \(2\le k\le\binom{d+2}{2}\), there are
\(K,\delta>0\) such that at all sufficiently small
fixed \(L>0\), the number of ordered \(k\)-tuples
whose repaired evaluation vectors at degree \(d\)
are dependent is
\[
                     O(X^{(k/3-\delta)L+\varepsilon}).       \tag{SIH3}
\]
The result applies throughout the balanced small-budget
endpoint. Here the use of a limiting one-third count
means an arbitrarily small fixed exponent allowance
for sufficiently small fixed budgets.

First exclude repeated repaired-coordinate ratios.
The ratio of two repaired letters determines the
two actual repaired values, since they are positive
and coprime. Their supports fix both original
letters at subpower cost, and the primitive addition
then determines the third one. Thus the corresponding
tuples number at most \(n^{(k-1)/3+\eta}X^\varepsilon\).

Next choose a small fixed \(\rho>0\). If the shared
radical exceeds \(n^\rho\), some pair of rows has
overlap at least \(n^{\rho/\binom{k}{2}}\). Fix
its first row and discount its full prime pool in
the one-third count of the second. This gives a
fixed positive loss in its exponent, uniformly in
the first row. The other \(k-2\) rows use the original
bounds. After the two exclusions, every dependent
tuple contains a minimally dependent subtuple of
size \(j\in[2,k]\). By (SIH2) its count, including
the unused slots, is at most
\[
               n^{(k-j)/3+\rho+3j/K+O_k(\eta)}X^\varepsilon.
\]
First choose \(\rho\) small, then \(K\) large, then the
exponent allowances small. These bounds give (SIH3).
[One-third count after fixed-pool discounts (EB5)](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)

This is a restriction on projective interpolation among
tuples of independently repaired points originating on
the single original additive line. For \(d=1\), the
three-point statement also has an argument based on the
six determinant terms. The simultaneous formulation
handles all fixed \(d\) and minimally dependent sizes up to
the dimension of the linear system, without a case split
over determinant partitions. The count of the original rows
would improve if their additive incidence forced sufficiently
many of these repaired dependencies or inexpensive
coefficient records. No such extraction follows from
interpolation alone.
[Projective incidence for independently repaired letters](ABC_REPAIRED_PRIMITIVE_PROJECTIVE_INCIDENCE.md)

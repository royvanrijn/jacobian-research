# An all-order valuation test for factorial moment windows

Fix a finite nonempty support of distinct exponent columns
\(\nu_1,\ldots,\nu_s\) in \(\mathbb N^d\), including a zero
column if present. For \(f=\sum_h c_h X^{\nu_h}\) with all listed
\(c_h\ne0\), the term of \(\Lambda_d(f^n)\) indexed by
\(m\in\mathbb N^s,\ \sum_h m_h=n\), is
\[
 W(n,m)=\frac{n!}{\prod_h m_h!}\prod_h c_h^{m_h}
                           \prod_{i=1}^d T_i!,
 \qquad T_i=\sum_h \nu_{h,i}m_h .
\]
The term indexed by \(n e_h\) is called pure. A pure term is a
valuation certificate for the moment if its valuation is strictly
smaller than the valuations of all other terms.

**Theorem.** Fix a prime \(p\) and any finite nonempty list of
valuations over \(p\), given by the rational values of each
\(c_h\), with \(v(p)=1\). For a prescribed positive window length
\(w\) there is a terminating decision procedure for the following
property: every \(w\)-window of positive orders contains a moment
with a pure valuation certificate at some place on the list.
If the property fails, the procedure gives a starting order
with a competitor for each pure term at each place in the
window. When it holds for valuations of the actual coefficients,
each such window contains a nonzero exact factorial moment.
In particular \(w=s\) gives the actual-support Strong Factorial
window in that case.

The input support and coefficient valuations are fixed, and the
moment order is unbounded. In the
[fixed-terminal carry graph](MULTIVARIATE_PADIC_DIGIT_AUTOMATON.md)
the exponent columns instead vary. Here the growing part of
each valuation comparison comes from ordinary base-\(p\)
position weights. The carry graph below retains it separately
from the signed digit-sum correction.

## Factorial terms and effective slopes

For each valuation \(v\), put
\[
 \lambda_h=v(c_h)+\frac{|\nu_h|-1}{p-1},
 \qquad |\nu_h|=\sum_i\nu_{h,i}.
\]
Let \(e\) be a common positive denominator for the coefficient
valuations at the chosen place. Legendre's formula, in terms
of the digit sum \(s_p\), gives
\[
\begin{split}
 &(p-1)e\bigl(v(W(n,m))-v(W(n,n e_h))\bigr)\\
 &\quad=\sum_\ell (p-1)e(\lambda_\ell-\lambda_h)m_\ell\\
 &\qquad\quad+e\left(\sum_\ell s_p(m_\ell)-s_p(n)
                  -\sum_i s_p(T_i)+\sum_i s_p(n\nu_{h,i})\right).
                                                               \tag{1}
\end{split}
\]
Here \(e_h\) is the count-vector basis vector, distinct from
the scalar denominator \(e\). Formula (1) follows by inserting
\(v_p(q!)=(q-s_p(q))/(p-1)\) in the multinomial and every row
factor. It also holds if a column is zero.

A pure term with a nonminimum slope cannot be least after a
computable initial range. Compare with any one pure term of
minimum slope. For a gap \(\delta>0\) the difference of their
valuations is \(n\delta\) plus digit-sum differences.
If \(n\) has \(L\) base-\(p\) digits and all exponent entries
are at most \(p^b\) for an integer \(b\ge0\), the size of the
digit-sum correction is at most \(d(L+b)\). This bound uses
nonnegativity of both digit-sum totals, each at most
\((p-1)d(L+b)\). Thus choose \(L_0\) with
\(\delta p^{L-1}>d(L+b)\) for every \(L\ge L_0\).
An explicit search terminates because
\(p^{L-1}/(L+b)\) is nondecreasing and tends to infinity.
Take the maximum cutoff over the finite list of places and
slopes. All windows starting below the cutoff can be checked
by direct finite comparison. For starts at or above the cutoff,
only minimum-slope pure terms need be considered.

## A finite path decision with positional and additive weights

**Lemma.** Let a finite directed graph have an initial state
and accepting states. Each edge bears \(a\in\mathbb N^r\) and
\(b\in\mathbb Z^r\). From initial values \(A=0,\ B\in S\)
update \(A\gets pA+a,\ B\gets B+b\) along a path, where
\(p\ge2\) is fixed and \(S\) is given as a semilinear set in
\(\mathbb Z^r\). It is decidable whether an accepting path
ends with \(A+B\le0\) coordinatewise.

A linear set has the form \(b_0+\sum_j\mathbb N v_j\),
with the periods and base allowed to be signed vectors.
The additive weights of paths in any finite subgraph are
effectively semilinear. One can turn the path language into
a regular expression and use union, sum, and the monoid
closure. For the last operation on a union of linear sets,
choose the subset of sets used at least once. The new base
is the sum of their bases, and the periods consist of
their periods and their bases. Each original period can
be allocated to its first chosen element. Include the
empty subset. This proves the assertion for signed weights
as well.

We use finite linear integer feasibility for the period
multipliers. For completeness, inequalities can be turned
into equations in nonnegative unknowns by adding slack
variables. A system \(C x=d\) of such equations is decided
by a finite carry graph reading digits of \(x\) from least
to most significant. Start the carry at \(-d\), update
\(c\gets(c+C x_{\rm digit})/p\) when this is integral, and
require zero final carry. The magnitude in each row stays
bounded by the maximum of its starting magnitude and the
sum of the absolute coefficients of that row of \(C\).
Hence these feasibility tests need only finite reachability.

Prove the lemma by induction on the number of constrained
coordinates, with the zero-coordinate case given by graph
reachability. Paths whose \(a\)-weights are all zero give
a semilinear test in \(B\) alone. Otherwise choose the
first edge \(f\) whose \(a\)-weight is nonzero. Before it,
\(A=0\). The possible values of \(B\) at its source,
on paths using only zero \(a\)-weights, form an effective
semilinear set. Fix one of its linear components
\(L=b_0+\sum_j\mathbb N v_j\).

Call a coordinate removable if a nonnegative integer
combination of the \(v_j\) is nonpositive everywhere and
strictly negative there. The finite integer test above
finds the set \(I\) of all such coordinates. Summing one
witness per coordinate in \(I\) gives a period combination
\(v_*\) strictly negative on \(I\) and zero on its complement
\(J\). For each *fixed* suffix starting with \(f\), the
constraints on \(I\) can be dropped exactly. Given a feasible
choice in projection to \(J\), lift it to \(L\) and add a
large multiple of \(v_*\). This makes all the dropped
inequalities true for that suffix, without altering the
other inequalities. The revised prefix ends at the same
state with \(A=0\), regardless of its length. If \(J\) is
empty, only reachability of acceptance after \(f\) remains.

If \(a(f)\) vanishes on \(J\), take \(f\), project to \(J\)
and recurse from its target state with offset set
\(L|_J+b(f)|_J\). Since \(a(f)\ne0\) before projection, at
least one coordinate was removed. In the other case some
\(i\in J\) has \(a_i(f)>0\). The real cone generated by
the \(v_j|_J\) meets the nonpositive orthant only at zero.
Otherwise a nonpositive direction strictly negative on some
coordinate of \(J\) can be chosen rational, by rational
linear feasibility, and then multiplied to have integral
multipliers. Add enough \(v_*\) to make the coordinates
on \(I\) nonpositive too. This contradicts the definition
of \(J\). Separate this closed polyhedral cone from
the compact simplex of nonpositive vectors with sum \(-1\).
It gives a vector \(u>0\) on \(J\) with \(u\cdot(v_j|_J)\ge0\).
One can find an integral such \(u\) by linear feasibility
with \(u_i\ge1\), since a real solution can be chosen
rational and scaled.

For a suffix of length \(h\) starting with \(f\), the
contribution of \(a_i(f)>0\) gives \(u\cdot A\ge p^{h-1}\).
Every \(B\) coming from the chosen prefix has
\(u\cdot B\ge u\cdot(b_0|_J)\) before the suffix. Its change
along the suffix is at least \(-M h\) for a fixed
\(M=\max_f |u\cdot(b(f)|_J)|\) over graph edges.
Thus \(A+B\le0\) requires
\[
                p^{h-1}+u\cdot(b_0|_J)-M h\le0 .
\]
The left side is strictly positive beyond an effective
bound, for example once positive with its first difference
positive for all larger \(h\). Enumerate the remaining
finite suffixes and use the linear integer test on the
prefix periods. This proves each branch of the induction,
including sufficiency when coordinates were dropped.
The lemma only uses nonnegative positional digits and
one-sided bounds on all terminal coordinate sums.

The graph lemma concerns a particular affine register
update. General affine Parikh automata permit arbitrary
affine updates and an arbitrary semilinear acceptance
set. Their undecidability result for emptiness does not
supply a decision procedure for this specialization.
See [Cadilhac, Finkel and McKenzie, §4 and Corollary 4.16](https://www.numdam.org/item/ITA_2012__46_4_511_0.pdf).

## Deciding a failed window without bounding its order

For starts at or above the cutoff, failure of the desired
certificate means the following. For each lag
\(0\le r<w\), each place, and each minimum-slope pure
index at that place, there exists a competitor vector at
order \(n+r\), distinct from the pure vector, with
the right side of (1) nonpositive. These are finitely many
simultaneous existential conditions, sharing only \(n\).

Read base-\(p\) digits initially from least to most significant.
Use one track for \(n\) and separate count tracks for the
competitor of each condition. The order \(n+r\) is obtained by
adding \(r\) as an initial carry. Bounded carries enforce
that the counts sum to that order and output the digits
of \(T_i\) and \((n+r)\nu_{h,i}\). For example a fixed linear
form with nonnegative coefficients \(b_q\) uses
\[
                c+\sum_q b_q a_q=p c'+d',
                   \qquad 0\le d'<p
\]
on the input digits \(a_q\). The carry is bounded by
\(\sum_q b_q\) and its starting value. Padding the tracks by
zero digits clears the carries at acceptance. Flags require
a nonzero digit outside the chosen pure count index and
\(n\) at least the cutoff. The latter can also be enforced
by adding a nonnegative slack track to the fixed binary or
base-\(p\) representation of the cutoff. This graph
accepts precisely all padded simultaneous witnesses.

Reverse the graph to read from most to least significant.
Use one pair of weights per condition. The digit increase
of the positional register is the linear combination of
the count digits in the first line on the right of (1).
It is nonnegative, since the reference pure slope is
minimum. The signed additive weight is the second line
using the single digits represented by the transition.
Thus the final register sum is exactly (1), with no
dependence on the padding. Apply the lemma with initially
zero additive offset. If no failed window is accepted, the
certificate holds at every order above the cutoff. If a
failed window exists, one can recover the actual tracks by
enumeration after the existence decision. At a certified
moment the strict minimum cannot cancel, by the
nonarchimedean valuation inequality. A failed certificate
only records term comparisons. The moment sum requires
additional information in those windows.

## A systematic tie in fully shared supports

**Proposition.** At a place over two, suppose every entry of
every column is positive and at least two pure terms have
minimum effective slope. There are arbitrarily late
two-windows without any pure valuation certificate at
this place.

For \(0<a<p^L\), the equality
\(a(p^L-1)=(a-1)p^L+(p^L-a)\) gives
\(s_p(a(p^L-1))=L(p-1)\) by digit complementation.
At \(p^L-p=p(p^{L-1}-1)\) the same sum is
\((L-1)(p-1)\) for \(a<p^{L-1}\).
Thus at both orders all pure terms have the same
row digit-sum correction to their linear slope.
The minimum pure terms tie with one another, and every
higher pure term is beaten by them. At \(p=2\) the orders
\(2^L-2,2^L-1\) are consecutive. This obstruction arises
even if every coefficient is a unit and two minimum
columns have the same total degree with different sorted
exponent entries. It imposes no vanishing on the full
moments.

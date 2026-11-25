# Exact weighted Taylor interpolation on torically sheared supports

This note studies the universal local Taylor-order guarantee in the
first determinant from ST1 and ST9. At a fixed degree stage, finite
exceptional primes may be removed. The comparison below does not
include extra local restrictions on the actual arithmetic sample.

## The local weighted interpolation problem

Let \(I=(I_1,\ldots,I_m)\) and \(j\) be the exponents of a monomial
\(z^Iv^j\). For the congruence \(v\prod z_i^{q_i}\equiv K\bmod p^e\)
introduce \(Y=v\prod z_i^{q_i}\) and \(x=I-qj\), so the same
monomial is \(z^xY^j\), with possibly negative \(x_i\).
Fix a unit chart center for every \(z_i\), and write
\[
 z_i=z_{0i}+p^{b_i}t_i,\qquad Y=K+p^e u.
\]
Here all \(b_i>0\), and \(Y\) is a unit. The binomial expansions
converge on the chart. Up to units common to a Taylor coefficient
vector or common to an original column, the coefficient vector
of \(\prod t_i^{a_i}u^k\) is the restriction to the sheared exponent
set \(S\) of
\[
 \prod_i\binom{x_i}{a_i}\binom{j}{k};
 \qquad\text{its weight is }\sum_i b_i a_i+e k .
\]
The unit factors do not change the order. Clearing the original
evaluation rows by local units also leaves it unchanged.

For a full minor with \(N=|S|\), the least valuation of a nonzero
coefficient of its Taylor expansion is the minimum total weight
of a set of \(N\) independent coefficient vectors. Indeed the
multilinear expansion is indexed by sets of distinct Taylor
multiindices. The alternating evaluation determinants of those
sets have disjoint monomial supports as polynomials in the
independent row variables. No two sets can cancel the coefficient
attached to the same monomial tuple. The nonzero rational
coefficient minors are units at good primes. There is a fixed
finite interpolation problem sufficient for the minimum. In each
variable, powers at or beyond the number of distinct coordinate
values reduce to smaller powers on \(S\), and that reduction
strictly lowers weight. The binomial-to-monomial change is
triangular and weight filtered. Thus one can exclude a finite
set of primes covering the denominators and nonzero minors among
these relevant vectors. This set is fixed at the given degree and
support stage. For larger good primes a nonzero leading
coefficient polynomial can also be realized at some row values
modulo \(p\). Its degree in each row variable is bounded in terms
of the finite interpolation problem. This identifies the best
universal valuation bound in the unconstrained chart model.

Since the local weight is the minimum of linear functions over
interpolation bases, it is concave and homogeneous in the
valuation weights. The same basis problem applies at every
good prime, so joint averaging also upper bounds the gain
available from prime-dependent depths. The algebraic basis
selection is the familiar finite-set interpolation problem;
the issue here is its comparison with primitive real height.
[Lederer, finite point ideals and fibre interpolation](https://arxiv.org/abs/math/0604133)

For a vertical fiber above \(x\), let \(m_x\) be its number of
occupied wider levels. Set \(X_k=\{x:m_x>k\}\). For each nonempty
\(X_k\), choose an interpolation basis \(F_{k,a}(x)\) on \(X_k\).
Then \(j^k F_{k,a}(x)\), over all \(k,a\), give an interpolation
basis on \(S\). To prove independence, in a linear relation first
restrict to fibers on which all coefficients of the highest-degree
polynomial in \(j\) must vanish. In particular the highest
coefficient vanishes on its \(X_k\), whence it is zero by its
interpolation basis. Descend in \(k\). The total number of vectors
is \(\sum_k |X_k|=N\). For one root one may take \(F_{k,a}(x)=x^a\)
for \(0\le a<|X_k|\).

## The broad profile on two occupied wider levels

Take the broad profile \(\nu=(1/3)\operatorname{Unif}[2,4]\)
on all three letters. Allow any bounded signed nonzero integer
multipliers, any asymmetric submeasures on the free letters, and
arbitrary one-root monomial support with exactly two occupied
wider levels. Write \(J\ge1\) for their gap. Let \(n_0,n_1\)
be the numbers of columns on those levels, \(N=n_0+n_1\),
and let \(c\) be the count of repeated horizontal projections in
the sheared support. Thus there are \(N-c\) distinct projections,
and exactly \(c\) of them occur twice. The above interpolation
basis gives the following upper bound for the exact universal
Taylor gain per row on the normalized scale:
\[
 G\le \frac{bL}{2N}\bigl((N-c)(N-c-1)+c(c-1)\bigr)
                  +\frac{wc}{N}.                              \tag{1}
\]
As usual \(b\) is shallow cost, \(w\) is full weight, and \(L\)
is the small common scale. The \(j^1\) vector costs one power
of the full modulus, independently of the wider gap. The bound
already controls arbitrary adaptive depth assignments, since it
comes from one explicit interpolation basis at each good prime.
Excluding finitely many primes at fixed \(L\) has zero normalized
cost on the broad bounded scaled degree interval.

Write \(s_W\) for the retained free mass on side \(W\) and
\(s=(s_A+s_B)/2\). Let
\[
 h=\frac12\sum_W\int|i_W(r)|\,d\nu_W'(r), \qquad
 R=\sum_W\int r\,d\nu_W'(r).
\]
On the sheared support coordinates \((Lx,j)\), the primitive
column mean is the sum of support-minus-centroid expressions
for the two sides, with normals \((\pm i,\pm r)\). The signs
reverse with the side. At any particular normal, move the column
on the unfavorable level to the favorable level at each
duplicated horizontal projection. The sum of column slack
decreases by exactly \(r J c\). On the favorable level there
are now its original \(n_f\) points, plus a copy of the \(c\)
duplicated projections. On the unfavorable level there are
\(n_u-c\) remaining points. Each of these three lists has
distinct horizontal integer coordinates. Measured from the
supporting endpoint in direction \(i\), the sum of their
horizontal slack is at least \(|i|L\,k(k-1)/2\) for a list
of \(k\) points. Integrating and using
\(n_f^2+(n_u-c)^2\ge(N-c)^2/2\) gives
\[
 B\ge\frac{R Jc}{N}
          +\frac{hL}{N}\left(\frac{(N-c)^2}{2}+c^2-N\right).
                                                                    \tag{2}
\]
The placement of the extreme horizontal points need not agree
between the two sides. The argument works for either sign of \(i\).

The root budget of a strict endpoint saving is
\(\mu+b<g=2s-1/3\), and \(b<g\le s\le h\). Rearrangement of
the free profile gives
\[
 R\ge4s+6s^2,\qquad
 w\le4\mu<4g,\qquad
 R-w>6(s-1/3)^2+2/3.
\]
Comparing (1) and (2) yields \(G\le B+hL/2\).
For an outer sequence with a fixed strict root-budget
margin and \(N\to\infty\), the positive slack from \(h-b\)
in the \((N-c)^2\) term offsets the integer end correction.
For bounded \(N\) and \(c>0\), \((RJ-w)c/N\) offsets it
instead. If \(N\) is bounded and \(c=0\), then \(G=O(L)\).
Choose two points on different wider levels. On either free
side the support-minus-mean is at least \(1/N\) times the
absolute difference of their values in each direction.
The difference has the form \(J(r+\alpha i_W)\) for a real
\(\alpha\). For bounded integer multipliers, its sublevel
sets \(|r+\alpha i_W|\le t\) have length \(O(t)\) uniformly
in \(\alpha\), by taking the union of one interval for each
possible multiplier. The density cap and positive retained
mass give a positive lower bound for \(B\). This excludes an improvement from the
exact support-dependent Taylor coefficients on any such
two-level support. If one uses the worst subset mean in ST13,
it is at least the full support's column mean, and in the
two-level proof the same argument can instead be applied directly
to the selected minor support.

For just two levels with disjoint projected horizontal intervals
one can also compute the exact weighted Hilbert function.
Let \(d_0\ge d_1\) be the allowed degrees in \(x\) for \(j^0,j^1\)
at a weight threshold, with \(d_1=-1\) if the second strip is
unavailable. The rank is
\[
 \min\{N,\ d_0+1+\min(d_1+1,n_0,n_1)\}.
\]
When \(d_1<\min(n_0,n_1)\), write the restriction on the levels
as \(A(x)\) and \(A(x)+B(x)\), with \(\deg B\le d_1\).
After \(d_1+1\) derivatives these real polynomials agree.
If \(d_0+d_1+2\le N\), Rolle's theorem gives more roots of their
common derivative than its degree, since the intervals are
disjoint. They then have degree at most \(d_1\), so both are
zero by their respective root counts. This proves injectivity.
For larger domains reduce \(d_0\) to attain dimension \(N\);
the reduced value remains at least \(d_1\). If
\(d_1\ge\min(n_0,n_1)\), interpolate on the smaller slice using
a \(j\)-linear selector, and then use the ordinary interpolation
rank on the larger slice. The \(d_1=-1\) case is immediate.
Thus unequal slice sizes are the only truncation of the mixed
strip in this separated two-level model.

## Fully toric Cartesian supports in any number of roots

A further obstruction applies to supports that, in sheared
coordinates \(x=I-qj\), are Cartesian products of independent
finite integer lists \(X_1,\ldots,X_m,J_0\). The exact weighted
interpolation basis is the product of univariate interpolation
bases. The universal Taylor mean is at most
\[
 \sum_i \frac{L B_i}{2}(\lvert X_i\rvert-1)
                  +\frac w2(\lvert J_0\rvert-1).
\]
For the two free letters a direction in group \(i\) involves
only \(x_i\) and \(j\). Because the support is Cartesian, its
support-minus-mean separates into these two coordinates.
Either one-sided mean distance along a list of \(k\) distinct
integer points is at least \((k-1)/2\). Consequently the
Archimedean mean is at least
\[
 \sum_i Lh_i(\lvert X_i\rvert-1)
                   +\frac R2(\lvert J_0\rvert-1),
 \qquad h_i=\frac12\sum_W\int_{\mathrm{group}\ i}
                              |i_W(r)|\,d\nu_W'(r).
\]
With \(B_i\le2s_i\le2h_i\), as required after freezing cheaper
whole dictionaries, and \(w<R\) as above, this exceeds the
Taylor gain if a vertical error coordinate is used. This is
independent of the number of roots and the shape of each list.
The nesting argument below handles a larger class of non-Cartesian supports.

## Nesting after the toric shear

The Cartesian obstruction extends to arbitrary horizontal sets
nested across the wider levels. More generally, for a nonempty
fiber-count set \(X_k=\{x:m_x>k\}\) and each coordinate \(i\),
let \(T_{k,i}\) be the smaller average distance of its integer
\(x_i\)-coordinates to their minimum or maximum. Write
\[
 T_i=\frac1N\sum_k|X_k|T_{k,i},\qquad
 \bar k=\frac1N\sum_k k|X_k|.
\]
One can choose the interpolation bases above so that their total
gain obeys
\[
                     G\le\sum_i B_i L T_i+w\bar k.               \tag{3}
\]
Indeed for each \(X_k\) translate and independently reflect the
coordinates to put the chosen extreme at zero. This yields
a set of nonnegative integer vectors. Index the interpolation
vectors by those same integer points \(d\) and use
\(\prod_i\binom{U_i}{d_i}\), where \(U_i\) are the translated and
reflected variables. The interpolation matrix is triangular
when ordered by total coordinate sum. The diagonal is one,
and nonzero entries require \(d_i\le U_i\) for all \(i\).
The polynomial has weighted degree at most \(\sum_i b_i d_i\).
Pairing it with \(j^k\) as above and averaging proves (3).

Suppose now that the wider levels can be listed \(j_0,j_1,\ldots\)
so that their horizontal sets are respectively
\(X_0\supseteq X_1\supseteq\cdots\), and \(j_0\) is an extreme
wider level. Orient the vertical coordinate and the two free
letter signs so that \(j_0\) is the minimum, and call the two
sides negative and positive according to their vertical signs.
For group \(i\) let \(h_{i,\pm}\) be the integral of the absolute
multiplier on that side, and \(R_{i,\pm}\) its retained numerical
moment. Thus \(\sum_{i,\pm}R_{i,\pm}=R\).

On the negative side a supporting point attaining the extreme
horizontal coordinate in the direction of its multiplier can
be chosen at \(j_0\), since that level contains the entire
horizontal projection. The horizontal mean distance from this
point is at least the average of the slice minima \(L T_i\).
The vertical mean distance from \(j_0\) is at least \(\bar k\).
Indeed for each fiber the average of \(m_x\) distinct integer
levels above its minimum is at least \((m_x-1)/2\) beyond the
minimum. Averaging over the fibers gives \(\bar k\). Consequently
the negative contribution of group \(i\) is at least
\(h_{i,-}L T_i+R_{i,-}\bar k\).

On the positive side the horizontal estimate \(h_{i,+}L T_i\)
still holds by comparing with the supporting endpoint on each
actual slice \(X_k\). Independently the vertical estimate
\(R_{i,+}\bar k\) follows by comparing with the upper endpoint
on each actual vertical fiber. Take the convex combination of
these estimates with horizontal weight
\[
 \lambda_i=\frac{(B_i-h_{i,-})_+}{h_{i,+}},
\]
interpreting \(0/0\) as zero. The assumption
\(B_i\le2s_i\le h_{i,-}+h_{i,+}\) gives \(0\le\lambda_i\le1\).
The total horizontal coefficient is now at least \(B_i\).
The total vertical coefficient lost by the convex combination
is at most
\[
 \sum_i\lambda_i R_{i,+}\le4\sum_i\lambda_i h_{i,+}
                                      \le4\sum_i B_i,
\]
because \(r\le4\) and all multipliers are nonzero integers.
The root budget is \(\mu+\sum_i B_i<g=2s-1/3\).
The broad moment lower bound above gives \(R>4g\), hence
\[
                       R-4\sum_i B_i>w .
\]
The primitive column mean dominates (3). This rules out all
nested toric supports under this budget, including non-Cartesian
horizontal sets and any number of roots.

When there is only one root, it is enough that an extreme wider
level contains the entire horizontal projection, without nesting
of the other slices. Use (3) with the sets \(X_k\) derived from
fiber sizes. The negative side still compares with the horizontal
global extreme at \(j_0\) and so pays the layer average \(L T_1\).
It also pays its vertical contribution, and the positive vertical
argument is unchanged. The shallow cost is bounded by the
negative side alone, since \(b<g=s_++s_- -1/3\le s_-\le h_-\).
The same reasoning works for several roots if
\(B_i\le h_{i,-}\) for all \(i\).

For lower-rank certificates that use one exact valuation bound
valid for every subset of a specified size \(N\), choose a
subset preserving the property used here. For nested slices take
the first levels in nesting order, with a partial last slice.
For a common extreme take that level first and fill the rest
arbitrarily. For two levels every subset still has at most two
levels. If the chosen subset has a single level, the vertical
Taylor weight is absent and the horizontal interpolation
comparison suffices. The worst Archimedean mean over all
subsets is at least the mean of this one, so this approach
cannot give a positive uniform comparison at smaller rank
either. These obstructions leave general nonnested multiroot
supports and additional arithmetic restrictions on the
evaluation rows open.


[The effective horizontal area bound](ABC_EXACT_TAYLOR_EFFECTIVE_AREA.md)
excludes arbitrary nonnested one-root supports of effective area at
most one on the same broad profile.

## Paying for short real intervals by additional finite-place precision

This comparison concerns a first determinant using the exact
error variable, with the same finite monomial support at all
places. Let \(\mathcal B\) be any selected set of columns,
including an arbitrary prescribed-size subminor, and write \(S\)
for the corresponding sheared exponent set \((I-qj,j)\).
For nonnegative horizontal weights \(\boldsymbol B\) and a
nonnegative vertical weight \(E\), write
\(\mathcal W(S;\boldsymbol B,E)\) for the minimum total weight
of an interpolation basis on \(S\). Here the weight of the
polynomial monomial \(x^\alpha j^k\) is
\(\boldsymbol B\cdot\alpha+Ek\). When some weights vanish, one
can restrict to the fixed finite monomial set obtained by
reducing in each coordinate at its distinct coordinate values.

Suppose, in addition to the finite-place chart, the real
variables obey \(z_i/z_{0i}=1+t_i\) with \(|t_i|\ll n^{-c_i}\).
Keep \(v\) in its usual height box, without an additional
real restriction on it. After taking the centres out of
the monomial columns, the determinant expands in products of
\(t_i\) and the original wider powers \(v^j\). Its generic
Taylor bound gives a discount in the original Archimedean column
weight sum equal to the formal weight
\[
                   \mathcal W(S;L\boldsymbol c,0).
 \tag{AR1}
\]
More precisely this is the least real Taylor weight of a
nonzero coefficient in the unconstrained formal model.
For a given \(j\), the coefficient vectors on its original
columns are restrictions of \(\prod_i\binom{I_i}{\alpha_i}\).
A nonzero determinant term must use precisely the original
number of columns for that \(j\), because vectors from
different \(j\)'s have disjoint column support. The weighted
basis problem therefore splits over \(j\). Vertical degree
is free when \(E=0\), so interpolating in \(j\) realizes that
same direct sum in the sheared problem. On any fixed
\(j\) the shear is a coordinate translation preserving
every horizontal weighted degree filtration. This proves
(AR1), using the independent formal row variables to see
noncancellation among the surviving coefficients. The product
of the extracted centre monomials and the primitive row
multipliers supplies the original Archimedean height bound,
up to an inner subpower loss. Every term uses the same
multiset of wider exponents from the chosen columns. Replacing
\(z_i\) by \(z_{0i}\) in those coefficient bounds changes
them by bounded factors at fixed support and height-box
scales. Counting the intervals themselves uses a normalized
chart cost \(\sum_i c_i+o(1)\) if one pays for the full
logarithmic covers of the \(z_i\).

The finite-place minimum has the superadditivity property
\[
 \mathcal W(S;\boldsymbol B+\boldsymbol D,E)
 \ \ge\
 \mathcal W(S;\boldsymbol B,E)+\mathcal W(S;\boldsymbol D,0)
 \tag{AR2}
\]
for \(\boldsymbol D\ge0\). Each basis has the sum of its two
weights on the right, and each weight is at least its
individual minimum. Thus the same weight discount can be
supplied by finite-place precision on a diffusely supported
full modulus. To make this precise, assume a positive
normalized prime mass on which \(e_p L\) is bounded below
and the maximum individual prime mass tends to zero.
Take the original total shallow costs bounded, with a fixed
number \(m\) of narrow root groups. For any large fixed \(K\),
choose a block of modulus primes with normalized mass
\(1/K+o(1)\), rounded up. At its primes increase the narrow
depths by \(D_i=\lceil Kc_i\rceil\). For each \(i\) the
primes where the original depth is within \(D_i\) of \(e_p\)
have normalized mass \(O(L)\), by the individual lower
bound for \(e_p L\) and the total depth budget. Discard those
primes before making the block choice. This permits all
the depth increases with at most
\(\sum_i c_i+m/K+o(1)\) additional normalized chart cost.
At those primes apply (AR2). Homogeneity in \(\boldsymbol D\)
and the block's prime mass show that the additional
finite-place gain is at least (AR1) after summing log weights.
Finitely many exceptional primes of the fixed support may
be discarded at each inner limit. Any strict margin in the
chart budgets absorbs \(m/K\).

Thus nominal Archimedean subdivisions in the narrow variables
alone cannot make the exact-error Taylor coefficient test
succeed where the flexible finite-place chart did not.
The argument applies to each selected subminor, with no
exchange between the extrema over column sets and bases.
It leaves open cheaper real covers of the actual sample,
simultaneous real restrictions on wider variables, and
further cancellations tied to actual arithmetic evaluations.

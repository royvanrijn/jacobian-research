# An exact count-simplex certificate for affine ordinary tails

Fix any finite ordinary core and finitely many unequal ordinary
tail types with exponent entries \(c_i M+a_i,c_i M+b_i\), where
\(c_i\) are positive integers and \(a_i,b_i\) are integers. The
multiplicities can vary independently with \(M\). There is a
finite rational test of the template alone which, if it passes,
gives a computable \(M_0\) such that every nonempty ordinary tail
of these types satisfies physical SFC(2) in every degree for
\(M\ge M_0\). The same test applies to repeated slopes, different
gaps, and mixed orientations.

The test uses the [absolute angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md)
\(J\), whose numerical value is unevaluated, as an input. Its
conditions involve the literal counts of affine factorial factors
at both orders. It allows templates in which all slopes have the
same \(2\)-adic valuation. The core affects only the height bound.

## The rational test

For each adjacent pair of orders \(n,N\) with \(4\le N<J\),
strip off the factorial at the smaller endpoint in each row.
The remaining factors at order \(j=n,N\) for type \(i\) are drawn
from
\[
 F_{j,i,h}(M)=j c_i M+j\ell_i+h,\qquad
 \ell_i=\min(a_i,b_i),\qquad 1\le h\le j|a_i-b_i|.
\]
At index \(k\), the row contributes the first \(k|a_i-b_i|\)
factors if it is increasing, and the first \((j-k)|a_i-b_i|\)
factors if it is decreasing.

Take as marks one representative \(x(M)\) from each class of
proportional affine factors in these two orders. For a mark and
a vector \(\mu\) in the closed simplex of tail type fractions, set
\[
 Y_{j,k}^{x}(\mu)=\sum_i\mu_i
    \#\{F_{j,i,h}\text{ contributed at }k:\ F_{j,i,h}
                          \text{ is proportional to }x\}.
\]
Let \(s^x_{j,a}(\mu)\), \(1\le a\le j\), be the successive
unit slopes of the lower convex polygon of these ordinates.
Each polygon edge slope is thus repeated by its horizontal
length. Define the nonnegative function
\[
 D_x(\mu)=
  \min_{\substack{a\ne a'\\b\ne b'}}
  \max\left(
    |s^x_{n,a}-s^x_{N,b}|,\ |s^x_{n,a'}-s^x_{N,b'}|
  \right).
\]
The test for this pair is
\[
                \delta_{n,N}:=\min_\mu\max_x D_x(\mu)>0.        \tag{1}
\]
In particular, it suffices that at every \(\mu\) some mark has at
most one slope shared by the two limiting multisets, counting
multiplicity.

All these functions are continuous rational piecewise linear.
Indeed the convex minorant at an integer \(k\) is the minimum of
\(Y_{j,k}\) and the interpolated ordinates at \(k\) on all chords
with integer endpoints bracketing \(k\). Each chord is rational
linear in \(\mu\). Minima, maxima, differences and absolute values
preserve rational piecewise linearity. Partitioning at the
corresponding rational linear equalities and minimizing the final
affine function on each closed cell computes (1) exactly.
The test fails exactly if there is a simplex point where
two disjoint slope matches occur for every mark.
One can also restrict the admissible fractions to any specified
closed rational polytope in the simplex and run the same test there.

## Transfer to actual ordinary arrays

For a fixed pair, write \(x(M)=A M+B\) and any available tail
factor as \(F(M)=A_F M+B_F\). The integer
\(\Delta_F=A B_F-A_F B\) vanishes exactly when \(F\) is
proportional to \(x\). Choose a fixed positive integer \(L\)
divisible by all the nonzero \(|\Delta_F|\) for all marks, and
by all the positive coefficients \(A,A_F\). For \(p\mid x(M)\)
and \(E=v_p(x(M))>v_p(L)\), the identity
\[
                       A F(M)= A_F x(M)+\Delta_F
\]
bounds \(v_p(F)\) by \(v_p(L)\) if \(\Delta_F\ne0\).
It bounds \(|v_p(F)-E|\) by \(v_p(L)\) otherwise.

Write \(r\ge1\) for the total tail count and \(\mu\) for its
vector of fractions. Take \(D=\max_i|a_i-b_i|\), and let
\(B_j\) be \(j!\) times the product of all stripped full core
ramps in order \(j\). Each native coefficient valuation at
order \(j\), divided by \(r E\), differs from \(Y_{j,k}^x(\mu)\)
by at most
\[
           \frac{v_p(B_j)+jD v_p(L)}{E}.
\]
This uses the fact that each stripped core coefficient divides
\(B_j\). Scalar stripping changes neither set of polygon
slopes. The convex minorant map on ordinates is \(1\)-Lipschitz
in the maximum norm by the chord formula. The resulting slopes
change by at most twice the ordinate error.

Assume (1) passes and choose a mark with
\(D_x(\mu)\ge\delta_{n,N}\). Put \(B_*=B_n B_N L^{JD}\).
For any prime with
\(E>\max(1,8/\delta_{n,N})v_p(B_*)\), the test still gives a
strictly positive \(D_x\) for the two actual slope lists,
both divided by the common factor \(rE\). In fact each slope
difference changes by at most \(4v_p(B_*)/E\). If
\(v_p(B_*)=0\), the same conclusion follows from the zero error.
The Newton polygon product rule puts the slope multiset of a
common positive-degree rational factor in both polygons.
Every nonconstant common factor has degree at least two, since
the even ordinary moment is strictly positive on the real axis
for distinct monomials. [Newton polygon product rule, Filaseta and Rocks, Theorem 2.1](https://doi.org/10.4064/cm6761-3-2016)

Consequently a common root forces the chosen \(x(M)\) to divide
\(B_*^{\lceil\max(1,8/\delta_{n,N})\rceil}\). At sufficiently
large \(M\), every mark is larger than that fixed positive
integer. Taking the maximum height over the finitely many
pairs, then using [ordinary low-order coprimality](MULTIVARIATE_TWO_MONOMIAL_NONDIVISIBILITY.md)
and the absolute angular cutoff, proves the certificate.

## A class with overlapping slope valuations

The test passes whenever for each even marked order \(m\) the
following holds. Among types with the same \(v_2(c_i)\), equality
of \((m\ell_i+1)/c_i\) implies equal orientations.
To see this at any simplex point \(\mu\), first restrict to the
types of positive mass. Choose one of maximal \(v_2(c_i)\),
breaking ties by minimizing \((m\ell_i+1)/c_i\), and mark
\(x=m c_i M+m\ell_i+1\). The intercept is odd. The other order
is odd, and comparing the two sides of a putative proportionality
after cross multiplication excludes every type of positive mass
there. It also excludes smaller \(v_2(c_i)\) at order \(m\),
after cancelling the common \(m\) from the equality.

Within the maximal tier at order \(m\), the minimizing intercept
ratio shows that only first factors of tied types can be
proportional to \(x\). They are all of one orientation. The
limiting marked ordinates are zero at exactly one endpoint
and a common strictly positive value everywhere else, with
the odd ordinates all zero. The marked lower polygon has one
nonzero slope repeated \(m\) times. This proves pointwise
positivity in (1), and compactness gives strict positivity
of its minimum.

For example, this covers arbitrary multiplicities of types
\((M-1,M)\) and \((3M,3M-1)\), although their slope valuations
coincide. With both counts \(\lfloor 3\tau M/4\rfloor\) the
tail has zero signed gap and curvature tending to any fixed
\(\tau>0\). More generally, common fixed gaps, zero upper
offsets, and at most one orientation per distinct slope meet
the criterion, regardless of the valuations of those slopes.


## Reflected pairs with arbitrary fixed offsets at one affine slope

The test also passes for any finite list of reflected pairs with
common positive slope \(c\), arbitrary fixed offsets, and arbitrary
independent multiplicities of the pairs. Require equal counts of
the two rows of each pair, equivalently use the simplex of pair
fractions. The positive slope is shared even when upper exponents
of different pairs are unequal.

To prove this take the even order \(m\ge4\). Among pairs of
positive mass choose one with largest midpoint \(cM+\bar a\),
where \(2\bar a\) is an integer. Mark
\(x=mcM+m\bar a+1\). Its intercept is an integer and is
nonzero modulo \(m\) since \(m\ge4\). At the other order \(m'\)
the original affine factors have leading coefficient \(m'c\).
Proportionality to the mark requires their intercept to be
\(m'(m\bar a+1)/m\), which is nonintegral because
\(\gcd(m,m')=1\). Thus all other-order ordinates are zero.

At index \(m/2\), every positive-mass row at order \(m\)
has argument at most \(m(cM+\bar a)\) and contributes no mark.
Away from the center at least one row of the chosen pair has
argument at least \(m(cM+\bar a)+1\), and its smaller-end
factorial precedes the mark. Hence the marked ordinate has
a unique minimum, at the center. A horizontal lower polygon
edge would have both endpoints at a global minimum, so this
polygon has no zero slope. Test (1) passes at every admissible
fraction, and its resulting height bound holds uniformly across
all independent pair multiplicities.

## Restricting the large-count test to the signed-gap hyperplane

For the unconstrained simplex, the large-count certificate only
needs the test (1) on
\[
                       \sum_i \mu_i(b_i-a_i)=0 .
\]
At any physical common root the
[absolute signed imbalance bound](MULTIVARIATE_ALL_DEGREE_SIGNED_IMBALANCE.md)
forces \(|\sum_i\mu_i(b_i-a_i)|\le H/r\), where
\(H=146J-74+|\Delta_C|\). If both orientations are available,
mixing \(\mu\) with a vertex of opposite signed gap puts it on
the hyperplane at distance at most \(2H/r\) in the \(\ell^1\)
norm. The gap at that vertex is at least one in magnitude, so
the required mixing weight is at most \(H/r\).
Each ordinate \(Y_{j,k}^x\) changes by at most the \(\ell^1\)
distance between the two fraction vectors, since each type
contributes at most one factor proportional to \(x\). Thus
\(\max_x D_x\) changes by at most four times that distance.

If the minimum on the balanced simplex is \(\delta>0\),
the actual fraction vector of any physical collision with
\(r\ge16H/\delta\) has \(\max_x D_x\ge\delta/2\).
The prime-depth proof therefore applies for such counts.
With only one orientation the balance slice is empty and
\(r\le H\) at any collision. Either way the remaining smaller
count vectors form a finite set. For each of those vectors
the stripped arrays have polynomial entries in \(M\), so
their ordinary resultant is a polynomial in \(M\). A nonzero
resultant gives an effective upper bound on all its positive
roots. An identically zero resultant gives actual ordinary
common roots for every sufficiently large integer \(M\),
since the stripped arrays then have their stated degrees.
This gives a finite exact decision for the entire large-height
behavior of any template that passes the balanced-slice tests.

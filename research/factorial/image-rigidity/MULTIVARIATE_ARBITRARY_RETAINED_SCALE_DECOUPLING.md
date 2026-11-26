# Exact separation of a populated exponent band from arbitrary retained rows

For a finite ordinary two-monomial profile, remove the rows whose
entries are equal. Let \(d_i=|b_i-a_i|\) be a row's gap, \(M_i\) its
larger exponent, and \(S=\sum_i d_i^2/M_i\). Fix an integer \(J\ge4\)
at least the [absolute angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md).
The value of \(J\) is unevaluated.

Select any group of \(r\ge1\) rows in the profile. For fixed integers
\(A,C\ge1\) and \(B\ge0\), assume that the smaller exponent and gap
of each selected row are
\[
        c_i M+\ell_i,\qquad d_i,\qquad
        1\le c_i\le C,\quad |\ell_i|\le BD,\quad D\le d_i\le AD
\]
for integers \(M,D\ge1\), \(c_i,\ell_i,d_i\). Either orientation
is allowed. All other rows are retained with unrestricted
orientations, heights, gaps and counts.

**Theorem.** Put
\[
\begin{split}
 H&=4J^2 C^2(A+B+1),\qquad K=16JC,\\
 T&=2C(2BK+1)((A-1)K+1),\qquad Q=T(4JT)^T,\\
 D_0&=\max\{32C,\,24(A+B),\,24C^2(J(A+B)+1),\,(24C)^2,\,
                              2^{4000CH}\},\\
 {\cal L}&=\sum_{\substack{i\text{ retained}\\ JM_i>HD}}
                              d_i\log(JM_i).
\end{split}
\]
If \(D\ge D_0\), a physical common root in any two consecutive
orders forces
\[
                 {\cal L}\ \ge\ \frac{rD\log M}{256 C J^2 Q}.       \tag{1}
\]
In particular, fix the band constants and a positive \(S_0\).
Suppose the selected band contributes at least \(S_0\) to \(S\).
Then for \(D\to\infty\), the original profiles are regular in every
degree if the largest retained exponent is \(o(M/D)\). The retained
rows can carry positive curvature.

The factor \(r\) in (1) counts the whole selected group, including
rows omitted by the box selection in the proof. The crucial point
is that an entire box of differing row types can contribute the
same normalized affine mark at all indices outside a short interval.
The proof uses the even order in each consecutive pair.

## Selecting a substantial set of boxes

The [absolute curvature twelve theorem](MULTIVARIATE_ABSOLUTE_CURVATURE_12.md)
gives the result when \(S\ge12\). Suppose \(S<12\). A selected row
gives \(D^2<12(CM+(A+B)D)\). The hypotheses on \(D\) imply
\[
                       M>D^2/(24C).
\]
The [first two adjacent pairs](MULTIVARIATE_TWO_MONOMIAL_NONDIVISIBILITY.md)
and the absolute cutoff leave \(4\le N<J\). Write \(m\) for the
even order in the pair and \(m'\) for the other order, so \(m\ge4\).

Group the selected rows into at most \(T\) boxes by their orientation,
\(c_i\), \(\lfloor K\ell_i/D\rfloor\), and \(\lfloor Kd_i/D\rfloor\).
Put \(b=4JT\) and \(t_j=(r/T)b^{-j}\) for \(0\le j\le T\).
Among the \(T+1\) disjoint intervals \([t_j/b,t_j)\), at least
one contains no positive box count. For that index put \(t=t_j\).
Keep each box of count at least \(t\), and call its rows active.
There is an active box since the largest count is at least \(r/T\).
Each active box has count at least \(t\ge r/Q\); the total count
\(r_0\) in the omitted boxes satisfies \(r_0\le t/(4J)\).
Choose one representative in each active box.

## Marks above the minimum of the representative envelope

Let \(c=2^{\max_i v_2(c_i)}\) where the maximum is over the active
rows. The highest tier consists of active rows with \(v_2(c_i)=v_2(c)\).
For \(0\le k\le m\), let
\[
 e(k)=\max_{i\text{ representative in the highest tier}}
              \frac{m\ell_i+d_i L_i(k)}{c_i},\qquad
 L_i(k)=\begin{cases}k,& b_i>a_i,\\m-k,& b_i<a_i.\end{cases}
\]
Write \(u=\min_k e(k)\), \(w=D/C\), and
\(I=\{k:e(k)<u+w/2\}\). Let \(k_*\) be a minimizing index.
The increasing branches of the envelope increase by at least \(w\)
per step, and the decreasing branches decrease by at least \(w\).
The maximum's sublevel set \(I\) is an interval because it is the
intersection of an interval ending on the right and one starting on
the left. If three consecutive indices lay in \(I\), every increasing
branch at the middle index would be below \(u-w/2\) by comparison
with the right index, and every decreasing branch would be below that
value by comparison with the left index. This contradicts the
definition of \(u\). Thus \(I\) contains \(k_*\) and has width at most
one.

Use the marks
\[
                 x_h=mcM+h,\qquad
                 h\text{ odd},\qquad
                   u+w/4<h/c\le u+3w/8.                          \tag{2}
\]
There are at least \(D/(32C)\) marks and at most \(D\). Indeed the
interval in \(h\) has length \(cw/8\ge D/(8C)\), and the step is
two. Its span is at most \(D/8\). We have
\[
        |h|\le C(J(A+B)+1)D,\qquad x_h>M,
\]
since \(|e(k)|\le J(A+B)D\) and
\(M>D^2/(24C)\ge C(J(A+B)+1)D\).

Strip each row's factorial at its smaller endpoint. Any
remaining factor at order \(j=m,m'\) from a selected row has form
\[
               F=j c_i M+j\ell_i+a,\qquad 1\le a\le j d_i.
\]
Its offset relative to the mark satisfies
\[
        mcF-jc_i x_h=mc(j\ell_i+a)-jc_i h,\qquad
             |mcF-jc_i x_h|\le HD.                               \tag{3}
\]
At the odd order \(m'\), this offset cannot vanish for active rows:
the first term on the right has \(2\)-adic valuation at least
\(v_2(c)+v_2(m)\), while the second has valuation \(v_2(c_i)\).
At the even order one can cancel \(m\) in the putative equality.
The same valuation comparison excludes active rows below the
highest tier. For a highest-tier active row, \(c_i h/c\) is an
integer and the offset vanishes precisely at \(a=c_i h/c-m\ell_i\).
At coefficient index \(k\) this factor occurs exactly when
\[
                   m\ell_i/c_i<h/c\le(m\ell_i+d_i L_i(k))/c_i .
\]

The normalized endpoints in one box differ from those of its
representative by at most \(2JD/K=w/8\). For every highest-tier
active row its normalized start \(m\ell_i/c_i\) is at most \(u+w/8\),
by comparison with the representative endpoint at \(k_*\). At
\(k_*\) all its endpoints are also at most \(u+w/8\), so no active
row contributes a zero-offset factor there. At every \(k\notin I\)
a representative achieving \(e(k)\ge u+w/2\) has all rows of its
box above the marked endpoint threshold \(h/c\le u+3w/8\).
Consequently at least \(t\) active rows contribute a zero-offset
factor at each such index. Each omitted row contributes at most one
zero-offset factor at any index of either order because the
right side of (3) increases strictly in \(a\).

## Prime depth and the retained product

Let \(B_j\) be \(j!\) times the product of all full stripped
ramps of the retained rows at order \(j\). For \(p>HD\) dividing
\(x_h\), put \(E=v_p(x_h)\). Equation (3) shows that the only
selected factors divisible by \(p\) are the zero-offset factors
and that each has valuation \(E\). After stripping all smaller-end
factorials, the ordinates at order \(m'\) are between \(0\) and
\(r_0 E+v_p(B_{m'})\). The absolute value of its slopes is at most
this upper endpoint. If \(s\) is such a slope, then the ordinate of
the order-\(m\) valuation array minus \(ks\), compared with its
value at \(k_*\), is larger at every \(k\notin I\) by at least
\[
        (t-(m+1)r_0)E-v_p(B_m B_{m'}^{\,m})
              \ \ge\ (3t/4)E-v_p(B_m B_{m'}^{\,m}).
\]
If this is positive, only indices in \(I\) can minimize for any
shared slope. The lower polygon then has shared slopes of total
width at most one. The [Newton polygon product theorem](https://doi.org/10.4064/cm6761-3-2016)
bounds the degree of a common rational factor by that width. The
even ordinary moment has no real roots by its positive integral
representation. Its nonconstant rational factors therefore have
degree at least two. A common root forces
\(E\le 4 v_p(B_mB_{m'}^m)/(3t)\).

Distinct marks differ by at most \(D/8\), so each \(p>HD\) divides
at most one mark. If \(P\) is the product of the marks and \(P_{>HD}\)
its part supported on these primes, we obtain
\[
            \log P_{>HD}
               \le \frac{4}{3t}\log (B_mB_{m'}^m)_{>HD}
               \le \frac{4J^2}{3t}{\cal L}.                     \tag{4}
\]
Each retained row contributes at most \(J^2 d_i\) factors to the
middle product, all at most \(JM_i\), and the binomial coefficients
have no primes above \(HD\).

For the reverse estimate let \(k\) be the number of marks. At
each small prime \(p\le HD\) remove a mark of highest valuation
when that valuation is positive. The remaining small-prime part
divides \((k-1)!\): the marks are odd, and their indices form a
single class modulo each \(p^e\) for \(p\ne2\) before removal of a
maximal term. The
[elementary prime count](MULTIVARIATE_LARGE_FACTORIAL_MARK_PRIME_BUDGET.md#uniform-estimates-for-the-certificate)
\(\pi(x)\le9x/\log x\) and the lower bound on \(D\) give
\(\pi(HD)\le D/(256C)\le k/8\). We also have
\(\log k\le\log D\le (2/3)\log M\) since \(M>D^2/(24C)\) and
\(D\ge(24C)^2\). Every mark exceeds \(M\), so
\[
          \log P_{>HD}
             \ge (k-\pi(HD))\log M-\log((k-1)!)
             \ge k\log M/5
             \ge D\log M/(160C).
\]
Combining with (4) and \(t\ge r/Q\) proves (1).

Finally, \(M_i\ge M-BD\ge M/2\) in the selected band. Thus its
curvature is at most \(2rA^2D^2/M\). When this curvature is at
least \(S_0\), we have \(rD\ge S_0 M/(2A^2D)\). For retained rows
of maximum height \(L\), the inequality
\(\sum_{\mathrm{retained}} d_i\le L\sum_{\mathrm{retained}}d_i^2/M_i\)
and \(S<12\) give \({\cal L}\le12L\log(JL)\) for \(L\ge1\).
If \(L=o(M/D)\), the last expression is
\(o((M/D)\log M)\), contradicting (1). An empty retained set has
\({\cal L}=0\). The argument works in the physical variable of
both ordinary moments. Rational coprimality extends by a Bézout
identity to characteristic-zero coefficient fields.


## Bounded gaps from any certified affine template

For a fixed affine template, take the admissible fraction polytope
and the mark gaps \(\delta_{n,N}>0\) of the
[count-simplex criterion](MULTIVARIATE_AFFINE_TAIL_COUNT_POLYTOPE.md).
Let \(r\) count the rows of this template in the profile, leaving
all additional rows unrestricted. Put
\[
                 {\cal L}_0=\sum_{i\text{ retained}}d_i\log(JM_i).
\]
There are explicit constants \(M_0,C_0>0\), depending only on
\(J\) and the template and fraction polytope, such that any
physical common root at height \(M\ge M_0\) forces
\[
                           {\cal L}_0\ge C_0 r\log M.             \tag{5}
\]
The positivity tests for the template still concern only its own
fraction vector.

To prove this keep the notation of the count-simplex criterion,
write \(D\) for the maximal fixed template gap, and take \(L_0\)
divisible by all nonzero affine offsets and positive multipliers
in its proof for every pair below \(J\). Choose the mark \(x(M)\)
with \(D_x(\mu)\ge\delta_{n,N}\). Use \(B_j\) for \(j!\) times the
full ramps of the arbitrarily varying retained rows. For a prime
of mark depth \(E>v_p(L_0)\), the error in the actual normalized
ordinates, divided by \(r E\), has bound
\[
       \left|\frac{v_p(\text{stripped coefficient at }k)}{rE}
                             -Y_{j,k}^{x}(\mu)\right|
              \le \frac{v_p(B_j)/r+jD v_p(L_0)}{E}.
\]
The \(1/r\) on the retained term is important. A chord gives each
lower minorant ordinate as a minimum of weighted averages of
array ordinates. Thus the ordinate-to-minorant map is \(1\)-Lipschitz
in the maximum norm, and the unit slopes are \(2\)-Lipschitz.
With \(C_1=\max_{4\le N<J}\max(1,8/\delta_{N-1,N})\), a possible
common factor forces every prime of \(x(M)\) to satisfy
\[
     v_p(x(M))\le C_1\bigl(v_p(B_nB_N)/r+JD v_p(L_0)\bigr).
\]
Indeed a strict violation ensures the depth hypothesis and changes
each difference of slopes by less than half of the guaranteed gap
\(\delta_{n,N}\). There can then be at most one matched unit slope
in the two polygons, which excludes a common rational factor by
even positivity. The case \(J=4\) is already covered by the cutoffs.

Summing against \(\log p\) and estimating the retained factorial
ramps gives
\[
 \log x(M)
 \le C_1\left(JD\log L_0+\frac{2\log(J!)+2J{\cal L}_0}{r}\right).
\]
For \(M\) large enough, every mark satisfies \(x(M)\ge M\).
Also require \(\log M\ge2 C_1(JD\log L_0+2\log(J!))\).
Using \(r\ge1\) proves (5) with \(C_0=1/(4JC_1)\).

In particular (5) holds for any bounded-gap band at one common
upper exponent \(M\), with arbitrary orientations and multiplicities.
The count-simplex criterion passes for this entire simplex.
At one fraction vector consider the types of positive mass.
If one sign occurs, after joint reversal if necessary an increasing
first-block mark occurs at every positive index in its type and
no mark occurs at zero. No other active row can supply this mark
at zero either. If both signs occur, mark at a nonmultiple
\(x=mM-h\) using the two smallest gaps in the
[common-upper middle-interval lemma](MULTIVARIATE_COMMON_UPPER_EXPONENT_DEEP_MARKS.md#a-common-middle-interval).
The marked count vanishes on an interval of width at most one
and is positive outside it. In either case \(m\nmid h\) implies
there is no proportional mark at the other order, so the two
limiting slope lists share at most one unit slope. Pointwise
positivity and compactness of the simplex give positive
\(\delta_{n,N}\). Here \(m\) can be the upper order.

If a bounded-gap certified band has curvature bounded below
and the other rows have maximum height \(o(M)\), the same
argument as above gives \({\cal L}_0=o(r\log M)\) in the
residual curvature regime \(S<12\). This excludes the
entire profile in all degrees. In particular the other rows
may have arbitrary independent heights and arbitrarily many
row types, even when they still carry positive curvature.


For this last asymptotic conclusion it suffices that the template
passes the tests on its signed-balanced fraction slice. At a
physical collision the
[absolute signed degree bound](MULTIVARIATE_ALL_DEGREE_SIGNED_IMBALANCE.md)
forces the template's mean signed gap to have magnitude at most
\((146J-74+12 L)/r=o(1)\). Over its compact admissible fraction
polytope, every limit of such vectors then lies on the balanced
slice. If that slice is empty, this already excludes an
unbounded collision sequence. If it is nonempty, positivity
of the continuous mark gap on the slice gives a uniform positive
lower bound on a neighborhood of it, and the preceding proof
applies with half the minimum on the slice.


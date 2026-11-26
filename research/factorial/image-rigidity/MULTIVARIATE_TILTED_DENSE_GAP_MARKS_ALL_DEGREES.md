# All-degree coprimality for dense gaps with tilted upper bands

Fix an integer \(B\ge0\). In an ordinary two-monomial profile,
omit equal rows and
select a nonempty group of unequal rows. Let \(X\) and \(M\)
be the largest gap and the largest upper exponent in the group.
Write \(d_i\) for the gap of row \(i\), and \(M-e_i\) for
its upper exponent. Partition the group into finitely many
bands, each of a single orientation. Give each band \(\nu\)
a slope \(\beta_\nu\in[0,B]\) and a density \(\lambda_\nu>0\).
Let \(A_\nu(t)\) count its gaps \(d_i\le t\), including
multiplicity, and write \(\Lambda=\sum_\nu\lambda_\nu\).

Let \(J\ge4\) be at least the
[absolute angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md).
Put
\[
       H=J^2(B+2),\qquad \theta=\frac{1}{64J^3(B+1)^3},
                    \qquad \varepsilon=\frac{\theta}{8J}.
\]
Assume that within each band
\[
 |e_i-\beta_\nu d_i|\le \varepsilon X,\qquad
 \sum_\nu\ \sup_{0\le t\le X}|A_\nu(t)-\lambda_\nu t|
                         \le \theta\Lambda X .
\]
The second inequality concerns the entire group of bands.
No comparison between the densities of the two orientations
is required, and one orientation may be absent.

For a retained row let \(M_i\) be its upper exponent and \(d_i\)
its gap. Put \(\mathcal L=\sum_{\rm retained}d_i\log(JM_i)\),
omitting rows with \(JM_i\le HX\) and rows whose two exponents
are both in \([M-X/(4J),M]\).

**Theorem.** If
\[
 X\ge\max\{144,8J\},\qquad
 \log X\ge576JH,\qquad
                  \mathcal L\le\frac{X\log M}{128J^3},
\]
the original profile satisfies physical SFC(2) in every degree
over every characteristic-zero field. Upper exponents may spread
over an interval of length comparable to the largest gap.
The value of \(J\) is unevaluated.

## Strictly convex counts from a mixture of bands

Fix \(4\le m=N<J\) and \(m'=m-1\). Use the mark
\(x_u=mM-u\) at the integers
\[
               X/4\le u\le X/2,\qquad u\equiv1\pmod m .
\]
At coefficient index \(k\), an increasing row contributes
the factor \(x_u\) when
\(me_i+(m-k)d_i\le u<m(e_i+d_i)\).
For a decreasing row the lower threshold is \(me_i+kd_i\).
Temporarily drop the strict upper bound and call the resulting
count \(Q_k\). As in the
[dense-gap count proof](MULTIVARIATE_DENSE_GAP_MARKS_ALL_DEGREES.md#counting-all-copies-of-a-factor),
the difference from the actual count is the number of full
ramps ending at or before \(u\), independent of \(k\).

In band \(\nu\), let \(h_\nu(k)=m-k\) in the increasing orientation
and \(h_\nu(k)=k\) in the decreasing one. Define
\[
       T(w,q)=\min(X,w/q)\quad(q>0),\qquad T(w,0)=X\quad(w>0).
\]
For \(q=h_\nu(k)+m\beta_\nu\), the band count is between the
values of \(A_\nu\) at \(T(u-m\varepsilon X,q)\) and
\(T(u+m\varepsilon X,q)\). When \(q=0\), all band rows
satisfy the lower threshold since \(u>m\varepsilon X\).
Both arguments are in \([0,X]\). If \(0\le q\le1/8\),
the arguments all equal \(X\), since \(u-m\varepsilon X\ge X/8\).
For \(q>1/8\), the map \(w\mapsto T(w,q)\) is \(8\)-Lipschitz.
Thus in all cases
\[
    Q_k=\sum_\nu\lambda_\nu T(u,h_\nu(k)+m\beta_\nu)+R_k,
                         \qquad |R_k|\le2\theta\Lambda X .
\]

For each band the scalar sequence \(T(u,q_0+k)\), where
\(q_0=m\beta_\nu\), is uncapped at every \(k\ge1\).
Its second difference at an interior index with no cap
in the three adjacent entries is \(2u/[q(q^2-1)]\),
where \(q=q_0+k\le J(B+1)\).
If its entry at zero is capped, its second difference at
index one is at least
\[
              X-2u+\frac{u}{q_0+2}
                      \ \ge\ \frac{u}{q_0+2}
                      \ \ge\ \frac{2u}{J^3(B+1)^3}.
\]
The same lower bound applies in the uncapped case and to
the reversed sequence for the increasing orientation.
The model's second difference is therefore at least
\(2u\Lambda/[J^3(B+1)^3]\ge X\Lambda/[2J^3(B+1)^3]\).
The second difference of the error is at most \(8\theta\Lambda X\).
Thus the actual count is strictly convex on the full grid.
Its minimizing indices span at most one. The comparison is
additive across the bands.

## The native prime comparison

If the total absolute curvature is at least twelve, use the
[absolute curvature theorem](MULTIVARIATE_ABSOLUTE_CURVATURE_12.md).
Otherwise \(M>X^2/12\), since a selected row has gap \(X\).
All marks are greater than \(M\) under the lower bound on \(X\).

Strip the smaller-end factorial from every row. For each selected
row, the offsets from the mark at order \(m\) are \(u-me_i-h\),
\(0\le h<md_i\). At the other order,
\[
         m(m'(M-e_i)-h)=m'x_u+(m'u-mm'e_i-mh),
                       \qquad 0\le h<m'd_i .
\]
The latter offset is nonzero modulo \(m\), and both sorts of
offset have magnitude at most \(HX\). Indeed \(e_i\le(B+1)X\)
and \(u\le X/2\). For any \(p>HX\) dividing \(x_u\),
the only selected factors that can contribute valuations are
the actual copies of \(x_u\) at the marked order.
The retained rows omitted from \(\mathcal L\) also give units
at such a prime. The low rows have every factor at most \(HX\).
In the top strip every order-\(m\) factor is strictly greater
than \(m(M-X/(4J))\ge x_u\). Both offset comparisons then apply
with \(e_i+d_i\le X/(4J)\).

Let \(B_j\) be the product of the full order-\(j\) stripped
ramps of the retained rows counted in \(\mathcal L\).
With \(E=v_p(x_u)\), the marked valuations are \(E\) times
the actual counts plus residual ordinates in \([0,v_p(B_m)]\).
Subtract \(E\) times the minimum count as a constant.
The other polygon has slopes bounded in magnitude by
\(v_p(B_{m'})\). If \(E>v_p(B_m)+m v_p(B_{m'})\), the
[deep-block comparison](MULTIVARIATE_CENTRAL_DEEP_PRIME_DIVISOR_GATE.md#marks-in-arbitrary-blocks)
confines every shared slope to the minimizing interval of width
at most one. Even-degree parity excludes a common factor.
Thus a possible common root forces every such \(E\) to be
at most \(v_p(B_mB_{m'}^m)\).

The progression has span less than \(HX\), count
\(s\ge X/(8J)\), and step \(m\) coprime to its terms.
Each prime above \(HX\) divides at most one mark. Deleting
one maximal-valuation term for each smaller prime leaves
small-prime part dividing \((s-1)!\), by the
[progression product comparison](MULTIVARIATE_LARGE_FACTORIAL_MARK_PRIME_BUDGET.md#a-product-bound-at-a-possible-common-root).
The bound \(\pi(HX)\le9HX/\log(HX)\) and the log hypothesis
give at most \(s/8\) deleted marks. The large-prime part
of the remaining product divides \((B_mB_{m'}^m)_{>HX}\).
Its logarithm is at most \(J^2\mathcal L\le s\log M/16\).
Since \(s\le X\le M^{2/3}\), the whole logarithm is
at most \((2/3+1/16)s\log M\), contradicting the lower bound
\((7/8)s\log M\) from the surviving marks.

Terminal two follows from even positivity, terminal three from
[quadratic nondivisibility](MULTIVARIATE_TWO_MONOMIAL_NONDIVISIBILITY.md#an-inequality-using-the-endpoint-coefficients),
and the absolute cutoff handles \(N\ge J\).
Rational Bézout identities transfer coprimality to every
characteristic-zero field.

For example, take both orientations of
\((M-2d+1,M-d+1)\) for \(d=1,\ldots,X\), with
\(M=\lfloor X^3/9\rfloor\).
There are two bands of unit density and slope one, each with
gap discrepancy at most one and upper-offset error one.
The upper spread is \(X-1\), and the selected lower endpoints
match in the two orientations. The signed gap remains zero
after adding any fixed balanced ordinary core. All factors
of such a core are eventually below the prime cutoff.
The theorem separates these unaugmented integer profiles
for all sufficiently large \(X\).

# A Newton polygon sieve for rational Gaussian bias slopes

Fix an ordinary factorial core, adjacent orders \(n,N=n+1\),
\(n\ge2\), and an integer signed tail gap \(h\). Write \(b_{j,k}\)
for the positive core moment coefficients and put
\(K_n=N\), \(K_N=n\), and \(\Lambda=N/n\). On a prescribed
rational bias slope \(\sigma\in[-1,1]\), the physical Gaussian
reference pair is
\[
 \begin{aligned}
 H_n(q,z)&=\sum_{k=0}^n b_{n,k}q^{K_n k(n-k)}z^k,\\
 H_N(q,\Lambda^h q^\sigma z)&=
    \sum_{k=0}^N b_{N,k}q^{K_N k(N-k)}
                                (\Lambda^h q^\sigma z)^k ,
                       \qquad 0<q<1 .
 \end{aligned}
\]
This corresponds to \(\lambda=\sigma\tau\) and
\(q=e^{-\tau/(2nN)}\) in the
[physical diffuse reference](MULTIVARIATE_DIFFUSE_PHYSICAL_EXCEPTION_POINTS.md#the-finite-set-of-physical-reference-contacts).

For every prime there is an exact test in the rational valuations
of the core coefficients that can exclude the entire positive
reference family at this slope and integer gap. It uses at most
\(n+N-2\) internal clock thresholds, plus affine matching of
edges on the resulting intervals. For nonintegral \(\sigma\),
the test produces a finite menu of possible valuations of \(q\).
Integral slopes can additionally have persistent local matches
on intervals. The criterion makes no genericity assumption
on the core.

## The Newton polygons depend monotonically on the clock valuation

A physical common root at real positive strength is nonreal
because the even reference is strictly positive on the real
axis. Its conjugate is a distinct common root. The two roots
give the same off-diagonal quotient in both parents, so the
[angular eliminant](MULTIVARIATE_ARBITRARY_CORE_DIFFUSE_EXTENSION.md#an-angular-eliminant-for-the-exact-core)
forces \(q\) to be algebraic. Since \(\sigma\) is rational,
\(q^\sigma\) is algebraic as well. For an embedding into
\(\overline{\mathbb Q}_p\), put \(t=v_p(q)\),
\(a_{j,k}=v_p(b_{j,k})\) and \(H_p=h v_p(\Lambda)\).
The Newton ordinates are
\[
    a_{j,k}+K_j k(j-k)t+
       \begin{cases}0,&j=n,\\ k(H_p+\sigma t),&j=N.\end{cases}     \tag{1}
\]
Adding the final affine term changes slopes but preserves the
vertices. By the
[angular clock bounds](MULTIVARIATE_GAUSSIAN_CLOCK_PADIC_EXCEPTIONS.md),
every candidate must satisfy \(-V_p\le t\le P_p\) using
the core bounds defined there.

For each internal index \(k\) of order \(j\) define its
deletion threshold
\[
 \Theta_{j,k}
    =\min_{0\le l<k<u\le j}
       \frac{\displaystyle
          \frac{(u-k)a_{j,l}+(k-l)a_{j,u}}{u-l}-a_{j,k}}
            {K_j(k-l)(u-k)}.                                     \tag{2}
\]
Then \(k\) is a strict lower hull vertex of (1) exactly
when \(t<\Theta_{j,k}\). To see this, an internal point
is a strict lower hull vertex precisely when it lies strictly
below every chord through one index to its left and one
to its right. This follows by comparing the slope from each
left point to this point with the slope from this point
to each right point. The largest left slope is strictly
below the smallest right slope exactly when such a strict
supporting vertex exists. The excess of the concave quadratic
\(K_j k(j-k)\) above the chord of that quadratic is exactly
\(K_j(k-l)(u-k)>0\). Thus each chord condition gives one
strict upper bound on \(t\), and their minimum is (2).
At equality the point can be omitted when forming the strict
vertices. Omitted collinear points do not affect Newton
edge widths or root valuations with multiplicity.

In particular the vertices can only disappear as the clock
valuation increases. Between consecutive values of the
thresholds and interval endpoints their set is constant.
For a lower edge with endpoint indices \(l,u\) and an
upper edge with endpoint indices \(i,k\), the respective
polygon slopes on such an interval are
\[
  \frac{a_{n,u}-a_{n,l}}{u-l}+N(n-l-u)t,\qquad
  \frac{a_{N,k}-a_{N,i}}{k-i}
                  +n(N-i-k)t+H_p+\sigma t .                     \tag{3}
\]
At a threshold one uses the exact hull given by (2).

## The finite test and a restriction on rational slopes

For a fixed prime, partition \([-V_p,P_p]\) by the thresholds
in (2). At the finitely many endpoints compute the overlap of
the lower and upper slope multisets, using edge widths as
multiplicities. In each open interval compare each pair of
active edges by (3). When their affine slopes are distinct
functions, equality holds at most at one explicit rational
\(t\) to be retained if it belongs to the interval.
An identical pair of affine functions contributes its common
width throughout the interval. If \(\sigma\) is nonintegral,
there is no such identical pair because the difference of
the two coefficients of \(t\) is an integer minus \(\sigma\).

Keep a valuation only if the two slope multisets share total
width at least two. If identified affine functions persist on
an interval, first count their total common width and retain
the full interval only when it is at least two. Otherwise
new matches are needed and occur at one of the isolated
candidate values. This is exact for the shared-width test.
The pointwise root valuation rule for Newton polygons makes
the negatives of the slopes the root valuations, with their
respective widths. The two distinct common roots from the
physical reference therefore need at least two units of common
width in every embedding. They may have the same valuation.
If one prime leaves no admissible valuation, the positive
physical reference family is excluded. A surviving value
supplies only local necessary conditions.

There is an additional slope restriction on intervals where
all edges of both parents have width one. A candidate then
requires two distinct simultaneous matches, say lower index
steps \(k,k'\) and upper steps \(\ell,\ell'\). Writing
\(A_{j,k}=a_{j,k}-a_{j,k-1}\), subtraction of their equalities
in (3) gives
\[
 t=\frac{(A_{N,\ell}-A_{N,\ell'})
                     -(A_{n,k}-A_{n,k'})}
          {2[n(\ell-\ell')-N(k-k')]} .                            \tag{4}
\]
Its denominator is nonzero because \(\gcd(n,N)=1\) and
\(0<|k-k'|<n\). At \(t\ne0\) either slope equality then
fixes \(\sigma\) rationally. If \(t=A/B\) is in lowest
terms, the equality has the form
\(\sigma t\in\mathbb Z+\mathbb Zt\). Hence the reduced
denominator of \(\sigma\) divides \(|A|\).
This gives a finite list of admissible rational slopes
for all nonzero valuations in the strict width-one regime,
without prescribing a slope first. A surviving \(t=0\)
does not impose this restriction.

## An exact ordinary core example

Take \(n=3,N=4\), the two shared rows \((7,3),(10,8)\),
\(h=0\) and \(\sigma=1/2\). Its ordinary core curvature
is \(94/35<12\) using the larger exponent in each row.
At \(p=3\) the core coefficient valuations in increasing
index order are
\[
                    (23,20,16,14),\qquad (31,27,26,21,19).
\]
The angular clock interval is \([-3,1/8]\).
The lower deletion thresholds are \(-1/8,1/8\); the upper
ones are \(1/9,-2/3,1/9\).
The endpoints, thresholds and active slope matches from (3)
give the following exact overlap widths.
\[
\begin{array}{c|rrrrrrrrr}
 t&-3&-6/7&-2/3&-1/8&0&1/11&1/9&2/17&1/8\\ \hline
 \text{shared width}&0&1&0&0&1&1&0&1&0
\end{array}
\]
Thus no positive Gaussian physical reference contact on this
rational slope and integer gap survives. The clock interval
alone is nonempty and single common slope matches occur at
four valuations, so it is necessary to use two-root width
in this calculation.

For a fixed core, order and gap, an exclusion by this test
gives actual ordinary stability near any compact positive
strength segment of that rational slope. Indeed the corresponding
finite physical exception list has positive distance from that
segment, and the linked exact coefficient and center errors
tend to zero uniformly in the ordinary tail. Both the
neighborhood and the diffuseness tolerance can be computed
by the [physical exception and transfer theorem](MULTIVARIATE_DIFFUSE_PHYSICAL_EXCEPTION_POINTS.md).
If each relevant order and each of its finitely many
allowed integer gaps admits an exclusion, the same argument
gives an all-degree tolerance on the segment using the
absolute high-order and curvature cutoffs. The sieve can
leave nonempty menus and does not give universal exclusion
for rational slopes.

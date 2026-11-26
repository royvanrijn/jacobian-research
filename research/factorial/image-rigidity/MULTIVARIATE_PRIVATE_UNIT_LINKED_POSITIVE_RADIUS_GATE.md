# An exact positive-radius gate for two private unit rows

Take an ordinary two-monomial support with one private unit
row on each side and arbitrary shared rows. Set
\[
       \Delta=\sum_{\rm shared}(b_i-a_i),\qquad
       V=\sum_{\rm shared}\frac{(b_i-a_i)^2}{\min(a_i,b_i)}.
\]
The matched tanh grid gives \(V=\sum_i2v_i/(1-v_i)\) over
shared Gauss factors, where
\(v_i=(n/T'_i-(n-1)/T_i)^{-1}\). The signed count of these
factors is \(\Delta\).

**Theorem.** At \(N=n+1\), \(n\ge2\), the polynomials
\(Q_n(z)\) and \(Q_N((n/N)^\Delta z)\) are coprime if
\[
                                  V\le\frac{6}{n+2}.           \tag{1}
\]
Here \(Q_m\) are the ordinary factorial moments in the
physical argument. More generally there is an exact rational
positive-radius test for this pair, given below.

Let \(J\ge4\) be at least the
[computable absolute angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md).
If \(V\le6/J\), the displayed compensated pair is coprime in
every degree. The original physical pairs are coprime in every
degree under this curvature cap whenever \(\Delta=0\) or
\(|\Delta|\ge6\). In particular the proof applies to balanced
shared rows of both orientations without a mirror assumption.
The result concerns the specified rescaling and physical
coprimality rather than disjointness under all rescalings.

## A scalar comparison using both orientations

Divide by the constant coefficient in each order. Write the
lower polynomial as \(P(z)\) and the compensated upper one as
\(Q(z)\), with leading coefficients \(a,b>0\), respectively.
Their coefficients are powers of the physical scales and
ordinary factorial ratios, so all are rational. Define
\[
             P^\vee(w)=w^nP(1/w)/a,\qquad
             Q^\vee(w)=w^NQ(1/w)/b,\qquad
             r_*=\frac{a}{b}.
\]
The nonleading coefficients in each comparison satisfy
\[
       [z^k]Q\le[z^k]P,\qquad
       [w^k]Q^\vee\le[w^k]P^\vee
                              \quad(0\le k\le n).             \tag{2}
\]
Indeed, the two private factorials cancel the binomial
coefficient and leave only shared factorial ratios. For
one orientation at fixed \(k\), the
[physical coefficient transport](MULTIVARIATE_BOTH_PRIVATE_PHYSICAL_COEFFICIENT_TRANSPORT.md)
says that compensating a degree-\(m\) coefficient by
\(m^{-\Delta k}\) makes it nonincreasing. On reversing,
\(\Delta\) and the comparison scaling both invert.
The two applications give (2).

Assume that
\[
                                 Q(r_*)\ge P(r_*).             \tag{3}
\]
This is the rational scalar test. To see that it suffices
write
\[
      Q(z)-P(z)=b z^N-\sum_{k=1}^n d_k z^k,\qquad d_k\ge0 .
\]
Its value divided by \(r^N\) on the positive axis is
nondecreasing. Thus the triangle inequality excludes a
nonzero common root of modulus \(r>r_*\), with strictness
unless all \(d_k=0\), in which case exclusion is immediate.
At \(1/r_*\), the reversed difference has the same sign:
\[
 Q^\vee(1/r_*)-P^\vee(1/r_*)
                  =\frac{r_*^{-n}}{a}\,[Q(r_*)-P(r_*)].
\]
Applying the same argument in reverse excludes \(r<r_*\).

If there are unequal shared rows, then \(d_n>0\).
In fact at \(k=n\ge2\) each increasing row contributes
a strictly positive loss to the derivative of its
compensated log coefficient, and the same is true for
each decreasing row. See the explicit sums below.
If equality in (3) gives equality in the triangle
inequality at \(r=r_*\), the term \(d_nz^n\) must align
with \(b z^N\); hence \(z\) is positive real, where the
moments have strictly positive coefficients.
If there are no unequal shared rows the difference
has only its top term. Thus (3) excludes a common root
also on the boundary circle.

## Estimating the scalar test by the matched curvature

For a shared row of smaller exponent \(\ell\), signed
gap \(d\ne0\), and starting exponent \(a_0>0\),
its physical coefficient is a Gamma factorial ratio.
At real order \(m\ge k\), subtract the log derivative
of this ratio from \(dk/m\).
For \(d>0\) the difference is
\(\sum_{h=1}^{dk} h/[m(ma_0+h)]\).
For \(d<0\) it is
\(\sum_{h=0}^{|d|k-1} h/[m(ma_0-h)]\).
In either case it is nonnegative and at most
\[
                          \frac{d^2 k^2+d k}{2m^2\ell}.
\]
It is positive for \(k\ge2\).
Summing, integrating \(m\) from \(n\) to \(N\), and using
\(1-e^{-x}\le x\) for \(x\ge0\), gives
\[
         0\le d_k\le [z^k]P\,
                \frac{V k^2+B_1 k}{2nN},
             \qquad B_1=\sum_{\rm shared}
                               \frac{b_i-a_i}{\min(a_i,b_i)}. \tag{4}
\]
The coefficient difference in the reversed orientation
has the same bound with its own lower coefficients and
with \(B_1\) replaced by \(-B_1\).

Let \(L=a r_*^n=b r_*^N\). If \(L<1\), reverse both
polynomials, which replaces \(L\) by \(1/L\) and does
not change the sign in (3). Thus use an orientation
with \(L\ge1\). Its lower geometric coefficients are
log-convex by log-convexity of each Gamma factor, so
\[
                        [z^k]P\,r_*^k\le L^{k/n}\le L .
\]
Combining this with (4) and the sums of \(k\) and \(k^2\)
gives
\[
        \frac{\sum_{k=1}^n d_k r_*^k}{L}
                     \le \frac{(2n+1)V+3B_1}{12}
                     \le \frac{(n+2)V}{6}.                  \tag{5}
\]
Here \(B_1\) is read in the chosen orientation and
\(|B_1|\le V\). The condition (1) ensures
\(\sum d_k r_*^k\le L=b r_*^N\), which is (3).
The first bound in (5) is an optional sharper
orientation-dependent test.

## The physical argument and the all-degree scope

At \(2\le n<N<J\) the inequality \(V\le6/J\)
implies (1). The first moment has a real root and the
second is strictly positive on the real axis, so they
are coprime under every positive rescaling.
The absolute angular theorem covers \(N\ge J\).
This proves the compensated all-degree result.
When \(\Delta=0\) the compensation is exactly one.

For the large signed drift under the same curvature
cap, write \(g_m\) for the positive endpoint scale of
the unadjusted degree-\(m\) polynomial. The
[linked center bound](MULTIVARIATE_ALL_DEGREE_SIGNED_IMBALANCE.md#an-explicit-linked-radial-condition)
in this two-unit case gives
\[
       \left|\log(g_N/g_n)-\Delta\log(N/n)\right|
                              \le \frac{V}{nN}.
\]
The endpoint-normalized geometric coefficients have
log derivative between \(-V/2\) and \(V/2\) on the
continuous index interval. This follows by integrating
the second derivative bound \(V/m\) between equal
zero log endpoints. The neighboring ratios therefore
put all root moduli within
\([e^{-V/2},e^{V/2}]\) by the Eneström–Kakeya theorem.
An actual common root requires \(|\log(g_N/g_n)|\le V\).
Using \(\log(N/n)>1/N\) then gives
\[
                             |\Delta|<V(N+1/n)<6
                         \qquad (2\le n<N<J,\ V\le6/J).
\]
Thus all \(|\Delta|\ge6\) are excluded as well.
The small nonzero signed drifts in this range need a
different comparison.

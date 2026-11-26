# All-degree ordinary coprimality with private variables on both sides

**Theorem.** Let \(a,b\) be two ordinary exponent columns. Suppose
each column has a positive entry in a row where the other column is
zero. Then their factorial moment polynomials
\[
 P_m(z)=\sum_{k=0}^m\binom mk
               \prod_i((m-k)a_i+kb_i)!\,z^k
\]
are coprime at every pair of consecutive positive orders. The result
holds over every characteristic-zero field. There is no bound on the
dimension, the shared row exponents, their orientations, or the
imbalance of the column totals.

The new ingredient is a first-sector argument using one of the two
coefficient-ratio inequalities that occur in
[endpoint Euclidean descent](MULTIVARIATE_BOTH_PRIVATE_EUCLIDEAN_DESCENT.md).
The inequality needed here always holds at the physical scale after
choosing the column orientation.

## A one-sided coefficient-ratio lemma

Let \(A(z)=\sum_{k=0}^n A_k z^k\) and
\(B(z)=\sum_{k=0}^{n+1}B_k z^k\) have strictly positive
coefficients. If
\[
                  \frac{B_{k+2}}{B_{k+1}}
                         \ge \frac{A_{k+1}}{A_k}
                          \qquad (0\le k<n),                    \tag{1}
\]
then they have no common root of the form \(z=r e^{i\theta}\) with
\(r>0,\ 0<\theta<\pi,\ n\theta<2\pi\).

Indeed, write \(g_k=B_{k+1}/A_k\) and \(u_k=A_k r^k\).
The sequence \(g_k\) is nondecreasing by (1). If \(A(z)=0\),
\(\sum_{k=0}^n u_k\sin(k\theta)=0\). On this index interval,
\(\sin(k\theta)\) is nonnegative until \(k\theta\) exceeds \(\pi\)
and is nonpositive thereafter. If both signs occur, choose \(g_*\)
between the last \(g_k\) with a positive sine and the first
\(g_k\) with a negative sine. Then each term in
\[
              \sum_{k=0}^n g_k u_k\sin(k\theta)
              =\sum_{k=0}^n(g_k-g_*)u_k\sin(k\theta)
\]
is nonpositive. If just one sign occurs the zero sum forces every
sine to be zero, giving the same inequality. But now
\[
          \Im\bigl(e^{-i\theta}B(z)\bigr)
              =-B_0\sin\theta
                 +r\sum_{k=0}^n g_k u_k\sin(k\theta)<0,
\]
which excludes \(B(z)=0\). The inequality also excludes the
conjugate sector. Applying the lemma to the two reciprocal
polynomials proves the same exclusion under the alternative
inequalities \(B_{k+1}/B_k\le A_{k+1}/A_k\) for \(0\le k<n\).

## The physical factorial link

Let \(A=P_n,\ B=P_N,\ N=n+1\). Write \(\Delta=\sum_i(b_i-a_i)\)
and let \(P_a,P_b\) denote the sums of private exponents on the
respective sides. Interchanging the two columns replaces each
polynomial by its reciprocal. This preserves coprimality and the
set of possible common arguments up to sign. Thus take \(\Delta\ge0\).
Both \(P_a\) and \(P_b\) are at least one.

The reverse half of the
[exact Mellin coefficient comparison](MULTIVARIATE_BOTH_PRIVATE_EUCLIDEAN_DESCENT.md#exact-factorial-link-and-the-radius-interval)
gives (1) whenever the normalized first reversed coefficient at the
upper order divided by that at the lower order is at most one.
For these ordinary profiles the bound there gives
\[
    \frac{A_{N,1}^-}{A_{n,1}^-}
       \le (n/N)^{\Delta_{\mathrm{sh}}+P_b-1}
       =(n/N)^{\Delta+P_a-1}\le1 ,
\]
where \(\Delta_{\mathrm{sh}}\) sums the signed gaps on shared
rows. To see the one-sided implication explicitly, reserve one
private factorial slot in each orientation. For the reversed
arrays the quotients of the constant-normalized lower and upper
coefficients are moments of a single positive product of Beta
powers. Their zeroth moment is one and their first moment is
\(A_{n,1}^-/A_{N,1}^-\ge1\). The ratios of successive moments
are nondecreasing by Cauchy–Schwarz, so the reversed lower
coefficient ratios dominate the reversed upper ratios. Taking
reciprocals and reversing the indices gives (1). This step uses
only the reversed comparison.

The first-sector lemma therefore excludes a common physical
root with \(0<|\arg z|<\pi\) and \(n|\arg z|<2\pi\) for every
both-private ordinary profile.

## First-sector extension with one private side

If only one column has private variables and its total degree is
strictly larger, the same first-sector exclusion holds. Orient
that column as \(b\), so \(\Delta=\sum_i(b_i-a_i)\ge1\), and
reverse both polynomial arrays. A private row of exponent \(e\)
in \(b\) becomes a decreasing row. Gauss multiplication writes
its factorial as a geometric factor times
\(\prod_{h=1}^{e}\Gamma(m-k+h/e)\), up to a scalar. The slot
of shape one cancels the \((m-k)!\) in the binomial coefficient.

After constant normalization the reversed coefficients have
the form
\[
  \frac{s^k}{k!}\prod_i
    \frac{\Gamma(m\widetilde a_i+\gamma_i+
                              \widetilde d_i k)}
         {\Gamma(m\widetilde a_i+\gamma_i)} .
\]
Here \(s>0\) does not depend on \(m\), ordinary rows have shape
one, and the remaining slots of the reserved row have
\(\widetilde a_i=1,\ \widetilde d_i=-1,\ 0<\gamma_i<1\).
The sum of the signed gaps \(\widetilde d_i\) is \(1-\Delta\).
The common factor \(s^k/k!\) cancels in the lower-to-upper
coefficient quotient, leaving the same positive Beta moments.

For the first reversed coefficient, each ordinary row
contributes logarithmic derivative at most
\(\widetilde d_i/m\) for real \(m\ge1\). An increasing row
contributes a sum of \(\widetilde a_i/(m\widetilde a_i+h)\)
with \(1\le h\le\widetilde d_i\); a decreasing row contributes
the negatives of those fractions with \(h=0,-1,\ldots,
1-|\widetilde d_i|\). Each residual slot contributes
\(-1/(m-1+\gamma_i)\le-1/m\). Thus the total derivative is
at most \((1-\Delta)/m\le0\), and the first Beta moment is at
least one. The reversed ratios give (1) in the original
orientation. This excludes every one-private common physical
root in the same sector, provided the private side has the
larger total degree.

## Exhaustion of the ordinary both-private cases

If both private totals are one, the
[angle boxes](MULTIVARIATE_PRIVATE_UNIT_ASYMMETRIC_ANGLE_BOXES.md)
place any common nonreal root in the same positive box after
conjugation. The [higher-label theorem](MULTIVARIATE_PRIVATE_UNIT_ALL_HIGHER_LABELS.md)
excludes every label at least two, with the already established
ordinary low terminals used there. In the first box the angle
satisfies \(\pi/n<\theta<2\pi/(n+2)\), so the first-sector argument
excludes it. The even-order polynomial has no real roots. In fact
its value at real \(z\) is the expectation of
\((X^a+zX^b)^m\) for independent unit exponential \(X_i\);
this is strictly positive for even \(m\), since the distinct
monomials cannot cancel identically on the positive orthant.

If both private totals are at least two, use the
[reciprocal binomial bound](MULTIVARIATE_BOTH_PRIVATE_EXPONENT_SUM_TWO.md).
If the smaller total is one and the larger two, use the
[private totals one and two theorem](MULTIVARIATE_PRIVATE_ONE_TWO_ALL_DEGREES.md).
If one total is at least three, the
[square-factorial phase bound](MULTIVARIATE_MIXED_PRIVATE_DEGREE_THREE.md)
and its low-terminal dependencies give the claim. These cases
cover all positive private totals. Finally, the polynomials have
rational coefficients, so complex coprimality gives a rational
Bézout identity that remains valid over every characteristic-zero
field.

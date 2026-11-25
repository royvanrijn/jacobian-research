# Eventual Strong Factorial for two monomials with private variables on both sides

Let \(\Lambda(X^\eta)=\prod_i\eta_i!\) be the ordinary factorial
functional. Fix distinct exponent vectors \(\alpha,\beta\in\mathbb
Z_{\ge0}^s\), and suppose
\[
 D_\alpha=\sum_{\beta_i=0}\alpha_i\ge1,\qquad
 D_\beta=\sum_{\alpha_i=0}\beta_i\ge1.
\]
Then there is \(M=M(\alpha,\beta)\) such that for every \(m\ge M\),
every characteristic-zero field \(K\), and every \((A,B)\ne(0,0)\)
in \(K^2\), at least one of
\(\Lambda((AX^\alpha+BX^\beta)^m)\) and
\(\Lambda((AX^\alpha+BX^\beta)^{m+1})\) is nonzero.

The previous [all-order private-exponent
theorem](MULTIVARIATE_TWO_MONOMIAL_NONDIVISIBILITY.md#all-orders-when-both-private-exponent-sums-are-at-least-two)
requires \(D_\alpha,D_\beta\ge2\) in the overlapping case.
The eventual theorem also covers \(XZ^2\) paired with \(YZ^3\),
and the equal-total-degree pair \(XZ^2W^4\) and \(YZ^5W\).
Both have private exponent sums one.

## Coefficient estimates

Assume \(AB\ne0\), and write \(a_i=\alpha_i\) and \(b_i=\beta_i\).
The coefficients of the moment polynomial in \(t=B/A\) are
\(\binom mk C_{m,k}\), where
\[
 C_{m,k}=\prod_i((m-k)a_i+kb_i)! .
\]
We may drop the coordinates with \(a_i=b_i\), since each such
coordinate contributes a nonzero scalar independent of \(k\).
If all remaining coordinates are private, apply the disjoint-support
consequence in the linked all-order theorem. We now have a coordinate
with \(a_i,b_i>0\) and \(a_i\ne b_i\).

Set \(G_m=(C_{m,m}/C_{m,0})^{1/m}>0\) and \(z=G_m t\).
Write the normalized polynomial as
\[
 P_m(z)=\sum_{k=0}^m c_{m,k}z^k,\qquad
 c_{m,k}=\binom mk\frac{C_{m,k}}{C_{m,0}G_m^k}.
\]
In particular \(c_{m,0}=c_{m,m}=1\).
There is a fixed \(q<1\) such that
\[
 0<c_{m,k}\le q^{\min(k,m-k)}\qquad(0\le k\le m).
 \tag{1}
\]

For a private exponent \(e\), the sequence \((ej)!/(j!)^e\)
is log-convex in \(j\). Its consecutive ratio is
\(\prod_{\ell=1}^e(e-(e-\ell)/(j+1))\), which increases with \(j\).
Thus
\[
 \frac{(ek)!}{((em)!)^{k/m}}\le
 \left(\frac{k!}{(m!)^{k/m}}\right)^e .
\]
The contributions of the private rows of the two columns are
therefore at most the product of \(k!/(m!)^{k/m}\) and
\((m-k)!/(m!)^{(m-k)/m}\). Each quotient is at most one, so the
extra private exponents also preserve this upper bound. The
product cancels \(\binom mk\). Log-convexity of the factorial
bounds each remaining row's normalized contribution by one.

For a varying shared row, reverse \(k\) if needed to write its
factorial as \((ma+dk)!\), with \(1\le a<b\) and \(d=b-a\).
For \(1\le k<m\), the second difference of its logarithm is
\[
 \sum_{j=1}^d\log\!\left(1+\frac{d}{ma+d(k-1)+j}\right)
 \ge\frac{d^2}{mb}.
\]
Indeed \(\log(1+x)\ge x/(1+x)\), and every resulting denominator
is at most \(mb\). Subtracting the quadratic with this second
difference and comparing with the endpoint chord shows that the
normalized contribution is at most
\(\exp[-d^2 k(m-k)/(2mb)]\).
Since \(k(m-k)/m\ge\min(k,m-k)/2\), this proves (1) with
\(q=\exp[-d^2/(4b)]\).

We also need the limit of each fixed coefficient at the two ends.
Put \(\Delta=\sum_i(b_i-a_i)\) and
\[
 G=e^{-\Delta}\frac{\prod_{b_i>0} b_i^{b_i}}
                      {\prod_{a_i>0} a_i^{a_i}}.
\]
Stirling's formula gives \(G_m/m^\Delta\to G\).
For fixed \(k\ge0\) and \(a_i>0\) the factorial ratio
\[
 \frac{(ma_i+(b_i-a_i)k)!}{(ma_i)!}
 \sim (ma_i)^{(b_i-a_i)k}
\]
follows from a finite product, or the inverse product when \(b_i<a_i\).
For \(a_i=0\) the ratio is the fixed factorial \((b_i k)!\).
Together with \(\binom mk\sim m^k/k!\), these formulas give a power
\(m^{(1-D_\beta)k}\) in \(c_{m,k}\).
If \(D_\beta>1\), the limit is zero for all \(k>0\).
If \(D_\beta=1\), the only positive private exponent of the second
column is one and its factorial is \(k!\).
This cancels the \(k!\) in the binomial asymptotic, so the limit is
\(\rho_L^k\) with
\[
 \rho_L=e^\Delta\prod_{a_i b_i>0}(a_i/b_i)^{b_i}.
\]
Indeed the rows with \(b_i=0<a_i\) contribute \(a_i^{-a_i k}\),
which cancels their factor from \(G^{-k}\).
Define \(\rho_L=0\) if \(D_\beta>1\), and define \(\rho_R\) by
swapping the columns. Then \(c_{m,k}\to\rho_L^k\) and
\(c_{m,m-k}\to\rho_R^k\), with the zeroth powers equal to one.
By (1) at \(k=1\), \(0\le\rho_L,\rho_R\le q<1\).

## Angular comparison

Let \(h_m=\lfloor m/2\rfloor\). Split the polynomial exactly as
\[
 P_m(z)=L_m(z)+z^m U_m(z^{-1}),\qquad
 L_m(w)=\sum_{k=0}^{h_m}c_{m,k}w^k,\quad
 U_m(w)=\sum_{j=0}^{m-h_m-1}c_{m,m-j}w^j .
\]
Choose \(1<R<1/q\). The bounds (1) and the fixed-coefficient
limits imply that on \(|w|\le R\), the two polynomials converge
uniformly with their first derivatives to
\((1-\rho_L w)^{-1}\) and \((1-\rho_R w)^{-1}\), respectively.
For the functions and derivatives the majorant sums are
\(\sum_{k\ge0}(qR)^k\) and
\(\sum_{k\ge1}kq^kR^{k-1}\). Both limit functions are bounded away
from zero on the disk. It follows that the two polynomials are
bounded above and bounded away from zero there for large \(m\).

On \(|z|=R^{-1}\) the term \(z^mU_m(1/z)\) then tends to zero
uniformly. By Rouché's theorem \(P_m\) has no zero inside or on
that circle for sufficiently large \(m\). Applying the same
argument to the reciprocal polynomial bounds its zeros outside
as well, so every zero of \(P_m\) has \(R^{-1}<|z|<R\).
At such a zero
\[
 m\log |z|=\log |L_m(z)|-\log |U_m(1/z)| .
\]
The right side is bounded independently of \(m\) and the zero.
Thus the moduli of the normalized zeros tend uniformly to one.

For \(r,r^{-1}<R\), let \(\phi_{m,r}(\theta)\) be the continuous
argument of
\(L_m(re^{i\theta})/U_m(r^{-1}e^{-i\theta})\) on \(0\le\theta\le\pi\),
with \(\phi_{m,r}(0)=0\). This choice is possible because both
polynomials are nonzero on the disk for large \(m\) and their
coefficients are positive. Its derivative equals
\[
 \operatorname{Re}\left(
 \frac{re^{i\theta}L'_m(re^{i\theta})}{L_m(re^{i\theta})}
 +\frac{r^{-1}e^{-i\theta}U'_m(r^{-1}e^{-i\theta})}
       {U_m(r^{-1}e^{-i\theta})}\right).
\]
This follows by differentiating the argument of the ratio.
If \(r\) ranges over the moduli of the zeros of \(P_m\), the
derivative converges uniformly on \(0\le\theta\le\pi\) to one
function independent of \(m\) and of the zero. In fact the
arguments of \(L_m,U_m\) converge uniformly to \(e^{i\theta}\)
and \(e^{-i\theta}\) respectively, and their logarithmic
derivatives converge uniformly on the disk by the established
bounds and derivative convergence.

Suppose the original moment polynomials at \(m,m+1\) share a
complex root \(t\). Both have real positive coefficients, so
conjugating if necessary gives \(0<\theta=\arg t\le\pi\).
Their normalized roots have the same angle \(\theta\) since
\(G_m,G_{m+1}>0\). Let \(r_n=|G_n t|\) for \(n=m,m+1\).
The two root equations give
\[
 n\theta=(2j_n+1)\pi+\phi_{n,r_n}(\theta),\qquad n\in\{m,m+1\},\ j_n\in\mathbb Z .
\]
For sufficiently large \(m\), the two phase derivatives differ
by strictly less than one throughout \([0,\pi]\).
The phases agree at zero, so subtracting the equations gives
\[
 |\theta-2\pi(j_{m+1}-j_m)|<\theta .
\]
For \(0<\theta\le\pi\) the distance to every integer multiple of
\(2\pi\) is at least \(\theta\), a contradiction.

The moment polynomials in \(t\) have rational coefficients.
A rational Bézout identity therefore transfers coprimality from
the complex numbers to every characteristic-zero field. If one
of \(A,B\) is zero, the nonzero monomial has nonzero moments
at every positive order. This completes the proof.

The threshold here depends on the fixed support. The finite prefix in
the overlapping cases with a private exponent sum one, and the supports
with no private variable on one or both sides, are outside this result.

## The left half-plane with exactly two private linear variables

If the private exponent sum of each monomial is one, then
every common root of two consecutive moment polynomials
would have strictly positive real part. This holds in the
original coefficient ratio and at every order, uniformly
over all ordinary shared exponent profiles.

Write the two monomials as \(XU\) and \(YV\), with \(X,Y\)
their respective private linear variables and with \(U,V\)
positive monomials in independent unit-rate exponential
coordinates. Integration in \(X,Y\) cancels the binomial
coefficients because \(\binom nk(n-k)!k!=n!\). Thus the
order-\(n\) polynomial divided by \(n!\) is
\[
 \mathbb E\!\left[\sum_{k=0}^n U^{n-k}(tV)^k\right]
   = A_n(t)-t^{n+1}B_n(t),\qquad
 \begin{split}
 A_n(t)&=\mathbb E[U^{n+1}/(U-tV)],\\
 B_n(t)&=\mathbb E[V^{n+1}/(U-tV)] .
 \end{split}
\]
For \(t=re^{i\theta}\) in the upper half-plane, these integrals
converge since \(|U-tV|\) is bounded below by a positive
multiple of \(U+V\). Both are nonzero. In the positive measure
\(U^n\) times the underlying probability measure, they are
the integrals of \(K(x)=(1-tx)^{-1}\) and \(x^{n+1}K(x)\),
with \(x=V/U>0\). The argument of \(K(x)\) strictly increases
in \(x\) and belongs to \((0,\pi-\theta)\).
Their arguments obey
\(0\le\delta_n=\arg B_n-\arg A_n<\pi-\theta\).
For the first inequality, symmetrize
\(\Im(B_n\overline{A_n})\) as a double integral.
Its integrand is half the product of
\(x^{n+1}-y^{n+1}\) and \(\Im(K(x)\overline{K(y)})\),
which is nonnegative because both factors have the
sign of \(x-y\).

At any upper-half-plane root,
\((n+1)\theta+\delta_n\) is an integer multiple of \(2\pi\).
Consequently a common root would make
\(\theta+\delta_{n+1}-\delta_n\) a multiple of \(2\pi\).
For \(\pi/2\le\theta<\pi\) this number lies strictly between
\(2\theta-\pi\) and \(\pi\). There is no such multiple.
Conjugation gives the lower-half-plane assertion.
On the negative real axis the even-order polynomial of the
pair is positive, by the defining integral of the even
power. The same conclusion therefore holds on the closed
left half-plane. \(\square\)

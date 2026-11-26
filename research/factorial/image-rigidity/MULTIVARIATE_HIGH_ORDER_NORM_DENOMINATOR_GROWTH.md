# High-order factorial differences force scaled norm denominators

Consider an ordinary two-monomial profile written as a core and a
nonempty tail of unequal shared ordinary rows. The core is any
finite profile, allowing equal or private rows. For each tail row
let \(m_i=(a_i+b_i)/2\), \(d_i=b_i-a_i\), \(x_i=d_i/m_i\), and
\(w_i=d_i^2/m_i\). Use
\[
        \eta=\max_i|x_i|\le\tfrac12,\qquad W=\sum_i w_i,\qquad
        B_r=\sum_i w_i|x_i|^{\,r-2}.
\]
Fix adjacent orders and let \(H\) be any nonconstant monic rational
common factor of degree \(d\) in the original physical variable.
Its root product \(c\) is positive and rational. Indeed the even
ordinary moment is strictly positive on the real axis, so the roots
of \(H\) occur in nonreal conjugate pairs, including multiplicities.

Choose either order \(j\), an integer \(2\le r\le j\), and
an integer \(0\le k\le j-r\). Write \(R_{i,j,k}\) for the \(i\)-th tail row's
factorial coefficient multiplier divided by its constant term,
and \(r_{i,j,\ell}=R_{i,j,\ell}/R_{i,j,\ell-1}\) for its
adjacent ratios. Set \(\varepsilon_i=(-\operatorname{sgn}d_i)^r\).
Define an integer vector \(\gamma\) on \(1,\ldots,j\) by
\[
       \gamma_{k+1+\ell}=(-1)^{r-1-\ell}\binom{r-1}{\ell}
              \quad(0\le\ell<r),
       \qquad\gamma_i=0\ \text{otherwise}.
\]
For each orientation choose fixed integer vectors
\(e^{(1,i)},e^{(2,i)}\) of length \(j\), constant within the
orientation, such that
\[
            \sum_{\ell} e^{(1,i)}_\ell
                =\sum_\ell e^{(2,i)}_\ell=d,\qquad
            e^{(1,i)}-e^{(2,i)}=\varepsilon_i\gamma .
\]
They always exist since \(\sum_\ell\gamma_\ell=0\).
For even \(r\) one can take the same two vectors for both
orientations. Write
\(S_t=\prod_{i,\ell}r_{i,j,\ell}^{e^{(t,i)}_\ell}\) for
\(t=1,2\), and \(E=\max_i\|e^{(1,i)}\|_1\).
Finally let \(g_C\) be the first binomial-stripped adjacent
coefficient ratio of the core at order \(j\), so \(g_C=1\) for
an empty core. Define the two rational scaled norms
\[
                         M_i=g_C^d S_i c .
\]

**Theorem.** The difference \(M_1-M_2\) is strictly positive.
If \(D_{\mathrm{diff}}\) is its reduced denominator, then
\[
 \begin{split}
       D_{\mathrm{diff}}\, B_r
        &\ge \frac1{C_{j,r,d,E}(W,\eta)},\\
 C_{j,r,d,E}(W,\eta)
        &=\frac{2^{jd}(r-2)!}{j^{\,r-1}(1-\eta/2)^{r-1}}\,
            \exp\!\left(
              \frac{W E}{2(1-\eta/2)}
            \right).
 \end{split}                                                       \tag{1}
\]
The same lower bound applies to the least common denominator of
\(M_1,M_2\). The estimate is uniform in the choice of the core,
the number of tail rows, the sizes of their exponents, and the
physical common roots. There is no limiting-contact hypothesis.

In particular, on any bounded curvature interval these simultaneous
denominators are \(\gg_{j,r,d,E,W}\eta^{-(r-2)}\).
The exact \(B_r\) in (1) can give a larger requirement, for
example \(B_r\le\eta^{r-4}B_4\) when \(r\ge4\). One can use the
larger order in the adjacent pair for \(r\). If an ordinary
collision is under consideration, the
[absolute curvature bound](MULTIVARIATE_ABSOLUTE_CURVATURE_12.md)
gives \(W<24\) whenever the tail is defined with midpoint
curvatures as above.

## Signed factorial differences

Extend \(\log R_{i,j,k}\) to \(F_i(t)\) for \(0\le t\le j\)
by
\[
               F_i(t)=g(ja_i+t d_i)-g(ja_i),
                   \qquad g(u)=\log\Gamma(1+u).
\]
The differentiated trigamma series gives
\((-1)^r g^{(r)}(u)=(r-1)!\sum_{h\ge1}(u+h)^{-r}\).
Integral comparison then gives
\[
     \frac{(r-2)!}{(u+1)^{r-1}}
          \le (-1)^r g^{(r)}(u)\le \frac{(r-2)!}{u^{r-1}},
                       \qquad u>0 .
\]
For every tail row the argument \(ja_i+t d_i\) is between
\(j m_i(1-\eta/2)\) and \(j m_i(1+\eta/2)\).
The \(i\)-th derivative acquires the factor \(d_i^r\), so
\(\varepsilon_i F_i^{(r)}\) is positive in both orientations.
By repeated application of the fundamental theorem of calculus,
\(\Delta^r F_i(k)\) is the integral of \(F_i^{(r)}\) over the
unit cube in \(r\) variables, with argument \(k\) plus their
sum. Thus
\[
     0<\sum_i\varepsilon_i\Delta^r F_i(k)\le
        \frac{(r-2)!}{j^{\,r-1}(1-\eta/2)^{r-1}}B_r .              \tag{2}
\]
The lower integral comparison also gives
\((r-2)!B_r/[j^{r-1}(1+\eta/2+\eta/j)^{r-1}]\)
as a lower bound, using \(1/m_i\le\eta\).

The identity \(\Delta F_i(\ell-1)=\log r_{i,j,\ell}\) gives
\[
              \log(S_1/S_2)=\sum_i\varepsilon_i\Delta^r F_i(k)>0 .
\]
The same log identity holds for \(M_1/M_2\). The normalization
by \(g_C^d c\) has no effect on this difference of logarithms.

## A uniform bound on the scaled norms

Let \(D_C,D_E\) be the positive endpoint argument scales of the
core and tail at order \(j\). Log convexity of the factorial core
coefficients, before the binomial factor, gives \(g_C\le D_C\).
The parent moment under the argument normalization
\(z\mapsto z/(D_C D_E)\) has endpoints one after its constant
coefficient is removed. All its intermediate coefficients are at
most the corresponding binomial coefficients. Every root of that
normalized polynomial therefore has modulus at most \(2^j\),
by its root equation when the modulus is at least one.
Each original root \(\alpha\) of the common factor consequently
satisfies \(|\alpha|D_C D_E\le2^j\).

Write \(D_{E,i}=\exp(F_i(j)/j)\), so \(D_E=\prod_i D_{E,i}\),
and put \(J_i(t)=F_i(t)-tF_i(j)/j\).
Its curvature satisfies
\(0\le J_i''(t)\le w_i/[j(1-\eta/2)]\) by the trigamma upper bound.
Since \(J_i(j)=J_i(0)=0\), the mean of \(J_i'\) over the entire
interval is zero. Integrating
\(|J_i'(t)-J_i'(s)|\le w_i|t-s|/[j(1-\eta/2)]\) in \(s\) gives
\[
   |J_i'(t)|\le \frac{w_i}{2(1-\eta/2)},\qquad
   \left|\log(r_{i,j,\ell}/D_{E,i})\right|
                       \le\frac{w_i}{2(1-\eta/2)} .
\]
As the exponents of each row in \(S_1\) sum to \(d\),
the normalized root bound now gives
\[
                  M_1\le
             2^{jd}\exp\!\left(
              \frac{W E}{2(1-\eta/2)}
                          \right).
\]
Writing \(L=\log(S_1/S_2)>0\), use
\(M_1-M_2=M_1(1-e^{-L})\le M_1 L\).
Combining with (2) proves
\(0<M_1-M_2\le C_{j,r,d,E}(W,\eta)B_r\).
The numerator of a positive reduced rational is at least one.
This proves (1) and the claim about the simultaneous denominator.

## Exact local budgets and cancellation

At a prime \(p\), let \(s_1,\ldots,s_d\) be the physical valuations
of the roots allocated to the common factor. With
\(V_p=v_p(g_C^d S_2)\) and \(U=S_1/S_2\), one has
\[
        v_p(D_{\mathrm{diff}})
              =\left(-\sum_{i=1}^d s_i-V_p-v_p(U-1)\right)_+ .     \tag{3}
\]
The quantities \(g_C,S_2,U\) are computed from exact factorial
blocks, independent of the candidate factor. The root valuation
sum is the only common-factor datum. The prime caps follow by
maximizing (3) over the shared slope allocations permitted by the
[parent polygon and residual tests](MULTIVARIATE_PHYSICAL_TRUNCATED_NEWTON_FACES.md).
Thus the resulting product must be at least the right side of
(1) divided by \(B_r\).

If \(v_p(U)\ne0\), then \(v_p(U-1)=\min(0,v_p(U))\), so
(3) equals the exponent of the least common denominator of
\(M_1,M_2\). If \(v_p(U)=0\), it is improved by the possible
cancellation \(v_p(U-1)\ge0\). In particular the budget does
not need to charge the trace denominator of a quadratic.
Common reshaping of the two exponent vectors in an orientation
by an integer vector of sum zero leaves \(U\) fixed and can change the budget in
(3). The bound applies separately to each such choice and to
each admissible segment in either parent order.

For an example of strict cancellation in this local budget,
take \(j=4,\ r=4,\ k=0\), an empty core and the one unit-gap
row \((2,3)\). Its four adjacent ratios are \(9,10,11,12\)
and \(U=10^3\cdot12/(9\cdot11^3)\).
The numerator and denominator have equal nonzero residue at
\(p=7\), and \(v_7(U-1)=1\). Thus at this prime the denominator
exponent in (3) loses one unit whenever the uncorrected
simultaneous exponent is positive. This example illustrates
the exact block cancellation, without assuming that this
profile has a candidate factor.

The argument extends the
[rational factorial norm scales](MULTIVARIATE_TRANSVERSE_QUADRATIC_DIRECTIONAL_HEIGHT.md#direct-local-budgets-for-the-denominators)
beyond quadratic contacts. It supplies an effective lower bound
even if the limiting roots and all limiting scaled norms are
rational. It remains necessary to bound the exact local budget
in (3) to exclude an ordinary profile.

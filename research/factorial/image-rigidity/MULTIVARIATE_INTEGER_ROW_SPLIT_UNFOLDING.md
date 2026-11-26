# Integer row splitting unfolds multiple Gaussian contacts

A weak retained row can be split into almost homothetic integer
rows. The leading Stirling terms cancel in the normalized
coefficient arrays. The half-logarithm gives a finite local
differential operator independent of the unknown angular jets.

Consider a physical interior contact of the
[fixed-terminal reference](MULTIVARIATE_FIXED_TERMINAL_COMPACTIFICATION.md#exact-compactified-data)
at adjacent orders \(n,N=n+1\), \(n\ge2\). Write \(T_*>0\) for the
Gaussian strength, \(|\lambda_*|<T_*\) for the bias, \(\Delta_0\)
for the total signed gap, and \(m_n,m_N\) for the multiplicities
of the two contacting nonreal roots. The countable ordinary retained
core has summable curvature. For \(m\ge2\), let \(h_m^-\) and
\(h_m^+\) be the smallest and largest strictly positive zeros of
the probabilists' Hermite polynomial
\[
                  \operatorname{He}_m(X)=e^{-\partial_X^2/2}X^m .
\]

**Theorem.** In each case below the contact is the limit of simple
transverse interior Gaussian contacts of finite ordinary cores at
the same orders and \(\Delta_0\). The strengths, biases and relative
physical root scales converge to those of the original contact.
Every strict total curvature cap persists.

1. If \(m_n=m_N\ge2\) is even, it suffices that infinitely many
   retained rows have absolute gap at least two.
2. If \(m_n=m_N\ge3\) is odd, it suffices that infinitely many
   retained rows have absolute gap at least four.
3. Suppose \(m_n,m_N\ge2\) and
   \[
                  \frac{h_{m_n}^+}{\sqrt n}>
                  \frac{h_{m_N}^-}{\sqrt N}.                 \tag{1}
   \]
   There is a computable bound on the absolute gap, depending only on
   \(n,m_n,m_N\), such that infinitely many retained rows beyond that
   bound suffice. In particular this holds for unbounded retained gaps.

There are no angular velocity or higher complex jet hypotheses in
these tests. Condition (1) holds whenever \(m_n\ge3\) and \(m_N\ge2\),
or when \(m_n=2\) and \(m_N\) is even. The claims concern reference
contacts with a remaining positive Gaussian strength.

## The local operator from splitting one weak row

Choose an integer \(r\ge2\). Select a sequence of rows with absolute
gap \(d\ge r\) and midpoint \(M\), whose curvatures tend to zero.
Summability guarantees this if infinitely many rows meet the gap
bound. Then \(d^2/M\to0\) and \(M\to\infty\). After choosing a
subsequence and reversing both columns if needed, take the gaps
positive. Put \(\eta=d/M\). Split \(d\) into \(r\) balanced positive
integer parts \(d_i\) and set \(p_i=d_i/d\). We have
\(\sum_i p_i=1\) and \(p_i\ge1/(2r-1)\). For a selected row
\((a,b)=(a,a+d)\), choose integer smaller endpoints \(a_i\) by
rounding \(p_i a\), with bounded errors and \(\sum_i a_i=a\).
For instance, round the first \(r-1\) down and use the remainder
for the last. The replacement rows are \((a_i,a_i+d_i)\).
Their midpoints are \(p_i M+\delta_i\), where
\(|\delta_i|\le r\) and \(\sum_i\delta_i=0\). They preserve both
endpoint sums and hence the signed degree difference.

For order \(j=n,N\), write \(s=k/j-1/2\). The selected row contributes
the factorial at \(jx\), where \(x=M(1+\eta s)\). The replacement
rows contribute the factorials at \(jx_i\), where
\(x_i=p_i x+\delta_i\). Analyze the log ratio as a function of
\(\zeta=\eta s\) on \(|\zeta|\le1/4\). Its arguments stay uniformly
in a closed right-half-plane sector with magnitude comparable to \(M\).
The uniform leading
[Stirling expansion, DLMF §5.11(i)–(ii)](https://dlmf.nist.gov/5.11)
gives
\[
 \log\Gamma(1+jx)
     =(jx+\tfrac12)\log(jx)-jx+\tfrac12\log(2\pi)+O_j(M^{-1})
\]
on this disk, with the analytic logarithm real on positive arguments.
In the difference of the leading terms \(jx\log(jx)\), expansion
at \(p_i x\) gives
\[
 \sum_i jx_i\log(jx_i)-jx\log(jx)
  = jx\sum_i p_i\log p_i+j\sum_i\delta_i\log p_i+O_{j,r}(M^{-1}).
\]
The term involving \(\log(jx)+1\) cancels since \(\sum_i\delta_i=0\).
The displayed leading terms are affine in \(\zeta\). The half-log
terms give \(\beta\log(1+\zeta)\) plus a constant and an error
\(O_{j,r}(M^{-1})\), where \(\beta=(r-1)/2\).
The remaining linear and constant Stirling terms are affine in
\(\zeta\) as well. Thus the exact log ratio, modulo affine terms, is
\[
                  \beta\log(1+\zeta)+E_j(\zeta),\qquad
        \sup_{|\zeta|\le1/4}|E_j(\zeta)|=O_{j,r}(M^{-1}).        \tag{2}
\]
All remainders above are analytic on this disk. For the leading
term this follows by Taylor's formula with second derivative
\(j/(p_i x+t\delta_i)\), uniformly bounded for \(0\le t\le1\).
For the half-log term the errors are the differences between
\(\log x_i\) and \(\log(p_i x)\). All quantities are real on real
\(\zeta\) in the disk.

The endpoint normalization subtracts the affine chord of this log
ratio at \(\zeta=\pm\eta/2\), so the large affine terms disappear
exactly. The chord of (2) then changes the reference only by a
nonzero scalar and a real radial dilation with log scale \(O_{j,r}(\eta)\).
This last bound follows from Cauchy's estimate on a smaller disk.
Remove that dilation when comparing angles. The remaining
coefficient multiplier is \((1+\zeta)^\beta e^{E_j(\zeta)}\).

Let \(\xi_j\) be the logarithm of the original root in its endpoint
scale, with common imaginary part \(\vartheta_*\).
Write \(G_j(x)=e^{-jx/2}F_j(T_*,e^x)\), where \(F_j\) is the
original endpoint-normalized reference polynomial. Then
\(G_j(\xi_j+x)=g_j x^{m_j}+O(x^{m_j+1})\), with \(g_j\ne0\).
On the finitely many frequencies \(k-j/2\), a polynomial in the
frequency acts as the corresponding polynomial in differentiation.
For \(T=T_*+\eta^2 y\), the Gaussian clock change contributes
\(e^{\eta^2 y\partial_x^2/(2j)}\) up to a scalar, by direct comparison
of the frequencies. In the coordinate \(x=\xi_j+\eta X\),
the normalized local limits after splitting are
\[
                 g_j e^{y\partial_X^2/(2j)}
                           (1+\partial_X/j)^\beta X^{m_j}.
                                                               \tag{3}
\]
For fractional \(\beta\), the binomial series terminates on the
monomial. To justify this limit, expand both analytic factors in
(2) to any fixed Taylor order \(K>\max(m_n,m_N)\) and truncate
the Gaussian exponential to order \(K\). The remainder on the
original finite frequencies is \(O_{j,r,K}(\eta^{K+1})\) for \(y\)
in a compact interval. It vanishes after division by
\(\eta^{m_j}\). Cauchy's estimate shows that every fixed Taylor
coefficient of \(e^{E_j}-1\) is \(O_{j,r,K}(M^{-1})\), so all
its terms vanish after operator dilation. Apply the finite
operators to \(\eta^{-m_j}G_j(\xi_j+\eta X)\), which tends locally
uniformly with every derivative to \(g_jX^{m_j}\).
This proves (3), including local uniform convergence in \(X\).

For a countable core, truncate after selecting the row so as to
include it and have coefficient error \(o(\eta^{\max(m_n,m_N)})\).
The normalized split multipliers and small clock changes are
bounded on the finite frequencies, so truncation does not change
the limit. At any two fixed clocks in (3) with strict opposite
signs of a ranked angular difference, the exact finite cores
preserve those signs by convergence of the root multisets.
The original roots are isolated in small logarithmic disks,
which can be kept free of boundary roots on the intervening
clock interval. The
[ordinary regularization lemma](MULTIVARIATE_MULTIPLE_CONTACT_ANGULAR_RANK_REVERSAL.md#ordinary-perturbations-with-simple-transverse-angular-contacts)
allows arbitrarily weak additional ordinary rows to enforce
simplicity of nonreal roots and transversality of every contact
on the interval. The sorted angles in the two clusters are
continuous, so a ranked difference changing sign yields the
desired simple transverse contact. No simplicity of the limit
polynomials at the two clocks is needed for the sign transfer.

## Equal multiplicities from a last real-root time

For common multiplicity \(m\), the root variables \(U=jX\) turn
the limit (3) into
\[
             Q_{m,\beta}(U,t)
                =e^{t\partial_U^2/2}(1+\partial_U)^\beta U^m,
                 \qquad t=jy.                                 \tag{4}
\]
The factor \(j^{-m}\) does not affect the zeros. Choose \(r=2\)
and \(\beta=1/2\) if \(m\) is even. At \(t=0\), the polynomial
\(Q_{m,\beta}\) is monic and its value at zero is
\(m!\binom{1/2}{m}<0\). Thus it has at least two real sign-changing
roots, and this remains true for all sufficiently small positive
\(t\).

For odd \(m\ge3\), choose \(r=4\) and \(\beta=3/2\).
The polynomial \(Q_{m,3/2}(-x,0)/m!\) is
\(L_m^{(3/2-m)}(x)\), as follows by expanding
\((1+\partial_U)^{3/2}\). It has at least two positive
sign-changing roots. To see this directly, write
\(\alpha=3/2-m\) and use the Rodrigues expression
\[
       L_m^{(\alpha)}(x)
        =\frac{x^{-\alpha}e^x}{m!}
                 \frac{d^m}{dx^m}\big(e^{-x}x^{m+\alpha}\big).
\]
It is orthogonal to \(1,x\) on \(x>0\) with positive weight
\(e^{-x}x^{-1/2}\). Indeed the two integrals reduce by parts to
the \(m\)-th derivatives of \(x^{m-2}\) and \(x^{m-1}\), which
vanish. At zero each boundary term has a positive power of
\(x\) because \(m+\alpha=3/2\), and at infinity each decays
exponentially. Fewer than two positive sign-changing roots
would allow a polynomial of degree at most one to absorb
all sign changes, contradicting that orthogonality.
Thus \(Q_{m,3/2}(U,0)\) has two negative sign-changing roots.
It also has a positive sign-changing root because
\(Q_{m,3/2}(0,0)=m!\binom{3/2}{m}<0\) and it is monic.
The three roots persist as real sign changes at small positive
\(t\).

For either parity, let \(\tau\) be the supremum of positive
times \(t\) at which (4) has more than \(m\bmod2\) real roots
counted with multiplicity. The preceding argument gives
\(\tau>0\). At large \(t\), scaling \(U=\sqrt t\,Z\) gives the
simple imaginary Hermite roots \(i h\), with \(h\) a zero of
\(\operatorname{He}_m\). Hence for all sufficiently large
\(t\) the only possible real root is the single real root
forced by odd degree, so \(\tau<\infty\).

Choose \(y_->0\) such that \(ny_-\) is a time with extra real
roots and \(Ny_->\tau\), which is possible by the supremum
and \(N>n\). The smallest index with a positive limiting
imaginary part at large \(y\) is
\(i_m=1+\lceil m/2\rceil\) in increasing order. At \(y_-\)
this index has imaginary part zero at order \(n\) and strictly
positive imaginary part at order \(N\), by conjugation and the
real-root counts. At large \(y_+\) the difference between
those same imaginary parts in \(X\) has leading term
\[
                  h_m^-\sqrt{y_+}
                    \big(n^{-1/2}-N^{-1/2}\big)>0.
\]
The signs are opposite. Apply the sign transfer after (3).

## Unequal multiplicities from the large-split limit

For a finite integer \(\beta\ge1\), set \(r=2\beta+1\).
In (3) put \(y=\beta v\), pass to \(Y=X/\sqrt\beta\), and
translate the roots by the real number \(\sqrt{\beta}/j\) in
the \(Y\)-coordinate. On polynomials, the centered operator
has exponent
\[
 \frac{v\partial_Y^2}{2j}
       +\beta\left(\log(1+\partial_Y/(j\sqrt\beta))
                             -\partial_Y/(j\sqrt\beta)\right)
     \ \longrightarrow\ \frac{v-1/j}{2j}\partial_Y^2 .
                                                               \tag{5}
\]
This is coefficientwise convergence of operators on the
fixed-degree polynomial space. Every term above the degree
of the monomial vanishes, and the coefficient of each
remaining power higher than two tends to zero. The real
translations do not affect the angles.

At a fixed \(v_-\) between \(1/N\) and \(1/n\), the pure heat
limit in (5) has real roots at order \(n\) and imaginary
roots at order \(N\). Compare the two full-list ranks that
give the outer positive Hermite root at order \(n\) and the
inner positive Hermite root at order \(N\) on the
positive-heat side. At \(v_-\) the difference is strictly
negative. At sufficiently large \(v_+>1/n\) the same
difference is strictly positive by (1), because the
imaginary displacements there are the respective Hermite
ranks times \(\sqrt{(v_+-1/j)/j}\). A finite integer
\(\beta\) preserves these strict signs in (5). It can be
found by isolating the Hermite zeros, choosing rational
\(v_\pm\), and testing the real and imaginary ranks of
the two algebraic polynomials in (3) for successive
\(\beta\). Root isolation and the strict limiting signs
ensure termination. With that fixed \(\beta\), apply
the sign transfer in (3) using any sequence of weak rows
with gaps at least \(2\beta+1\).

The examples following (1) use ordinary orthogonal
polynomial interlacing. The largest positive Hermite zero
increases with \(m\). The smallest positive zero decreases
along either parity, since \(\operatorname{He}_{2a}(X)\)
and \(\operatorname{He}_{2a+1}(X)/X\) are constant multiples of
\(L_a^{(-1/2)}(X^2/2)\) and \(L_a^{(1/2)}(X^2/2)\), respectively.
These identities follow by comparing polynomial coefficients,
and the zero comparisons follow from interlacing for the
positive Gaussian and Laguerre weights.
Thus \(h_{m_n}^+\ge\sqrt3\) for \(m_n\ge3\),
\(h_{m_N}^-\le\sqrt3\) for \(m_N\ge2\), and
\(h_{m_N}^-\le1\) for even \(m_N\). The extra factor
\(\sqrt{N/n}>1\) gives the strict comparisons stated above.

## Physical scale and curvature

The split preserves the row's signed gap. Its replacement
midpoints tend to infinity and its total curvature tends
to zero together with the original row's curvature.
Indeed \(d_i=p_i d\), \(M_i=p_i M+O_r(1)\), and the \(p_i\)
are bounded below for fixed \(r\). Each relative gap tends
to zero. The retained center changes by \(o(1)\), by the
[curvature bound on its row terms](MULTIVARIATE_FIXED_TERMINAL_COMPACTIFICATION.md#exact-compactified-data).
Truncation and the arbitrarily weak genericity rows change
that center by \(o(1)\) as well. Maintain \(\Delta_0\) by
adjusting the signed gap of the residual Gaussian data.
For equal-angle roots \(z_n,z_N\) of the resulting
endpoint-normalized pair, choose
\[
             \lambda_{\rm fin}
               =2nN\bigl(\Gamma_{\rm fin}-\log|z_N/z_n|\bigr)
\]
in the notation of the fixed-terminal reference. The radial
normalization used in (2) tends to the identity and the
local roots tend to the original roots. Hence
\(\lambda_{\rm fin}\to\lambda_*\), while
\(T=T_*+\eta^2y\to T_*\) at both endpoints and throughout
the interval between them. The strict interior condition
and every strict curvature cap therefore persist.

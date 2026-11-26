# Every-degree coprimality from interlaced coefficient ratios

Let \(Q_m(z)\) be the ordinary factorial moment polynomial of two
distinct exponent columns \(a,b\), with at least one private
variable on each side. Put \(\Delta=\sum_i(b_i-a_i)\), let
\(P_a,P_b\) be the respective private totals, and use
\(\Delta_{\rm sh}\) for the sum of gaps on the shared rows.
For \(m\ge1\), let \(A_{m,1}^+\) be the coefficient of \(z\) after
constant normalization, and let \(A_{m,1}^-\) be the analogous
coefficient for \(z^mQ_m(1/z)\).

**Theorem.** For every \(n\ge1\), \(N=n+1\), the polynomials
\(Q_n(z),Q_N(\lambda z)\) are coprime whenever
\[
                 \frac{A_{N,1}^-}{A_{n,1}^-}
                 \ \le\ \lambda\ \le\
                 \frac{A_{n,1}^+}{A_{N,1}^+}.                  \tag{1}
\]
Their resultant in this order is strictly positive. The interval
contains the full interval
\[
       \left[(n/N)^{\Delta_{\rm sh}+P_b-1},\
                    (n/N)^{\Delta_{\rm sh}-P_a+1}\right],
                                                               \tag{2}
\]
including the signed compensation \((n/N)^\Delta\).
In particular physical SFC(2) holds in every degree if
\(1-P_a\le\Delta\le P_b-1\). This includes **all balanced
both-private ordinary profiles**, with arbitrary private powers
and arbitrary shared rows of either orientation.
The proof uses no terminal cutoff or angular estimate. For
rational \(\lambda\), the assertion holds over every
characteristic-zero field.

More generally two roots of \(Q_n,Q_N\) with the same argument
cannot have their upper-to-lower radius ratio in (1), irrespective
of the angular label. In the unbalanced case the actual physical
test is whether (1) contains one. The guaranteed compensation can
differ from one.

## Endpoint cancellation preserves interlaced ratios

Let \(A(z)=\sum_{k=0}^n A_k z^k\) and
\(B(z)=\sum_{k=0}^{n+1}B_k z^k\) have positive coefficients.
Write \(r_k=A_{k+1}/A_k\) and \(q_k=B_{k+1}/B_k\), and suppose
\[
                     q_k\le r_k\le q_{k+1}
                         \qquad (0\le k<n).                  \tag{3}
\]
Put \(c=B_0/A_0\), \(d=B_{n+1}/A_n\), and define
\[
                  D(z)=\frac{(c+dz)A(z)-B(z)}{z}.
\]
Both endpoint terms in the numerator cancel. At index \(k\)
of \(D\), \(0\le k<n\), set
\(\alpha=cA_{k+1}\), \(\beta=dA_k\), and \(\gamma=B_{k+1}\).
The forward inequalities of (3), multiplied starting at index
zero, give \(\gamma\le\alpha\). The reverse inequalities,
multiplied from the leading endpoint, give \(\gamma\le\beta\).
Thus \(D_k=\alpha+\beta-\gamma>0\), so \(D\) has degree \(n-1\).

For \(0\le k<n-1\), the next coefficient is
\[
                 D_{k+1}=\alpha r_{k+1}+\beta r_k-\gamma q_{k+1}.
\]
In particular
\[
\begin{aligned}
 D_{k+1}-r_kD_k
    &=\alpha(r_{k+1}-r_k)-\gamma(q_{k+1}-r_k)\ge0,\\
 r_{k+1}D_k-D_{k+1}
    &=\beta(r_{k+1}-r_k)-\gamma(r_{k+1}-q_{k+1})\ge0 .
\end{aligned}
\]
Both use \(r_k\le q_{k+1}\le r_{k+1}\) and the corresponding
bound on \(\gamma\). This proves the same interlacing condition
for the pair \(D,A\). At any common root of \(A,B\), the variable
is nonzero since the constant coefficients are positive. Hence
it is also a root of \(D\). Iteration would give a root of a
positive constant, a contradiction.

This comparison also computes the sign without approximating
roots. At each root \(\alpha_j\) of \(A\) one has
\(B(\alpha_j)=-\alpha_j D(\alpha_j)\), and
\(\prod_j\alpha_j=(-1)^n A_0/A_n\). The resultant product
formula gives
\[
                     \operatorname{Res}(A,B)
                         =A_0 A_n\operatorname{Res}(D,A).
\]
Here the reversal of the last two arguments changes no sign,
since \(n(n-1)\) is even. Thus the resultant is a product of
strictly positive endpoint factors, ending with the positive
constant. The endpoint cancellation is the recurrence used
for Laurent biorthogonal polynomials. The interlaced-ratio
positivity proof above is elementary and uses no orthogonality
measure. [Vinet and Zhedanov, Spectral transformations of the Laurent biorthogonal polynomials, II](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/3C4B7EA6A4F71FB4C932C2F2FD12863E/S0008439500014442a.pdf/spectral-transformations-of-the-laurent-biorthogonal-polynomials-ii-pastro-polynomials.pdf)

## Exact factorial link and the radius interval

Reserve one private row in each orientation. Gauss multiplication
on a private row of exponent \(e\) writes its factorial \((ej)!\)
as a positive exponential factor times
\(\prod_{h=1}^e(h/e)_j\). The slot of shape one contributes
\(j!\). The two reserved slots cancel the binomial factorial
denominators. The remaining slots of these two rows have unit
gaps, one nonzero endpoint, and Gamma shapes in \((0,1)\).
Unreserved rows retain shape one. The removed unit slots have
opposite orientations, so the remaining signed total gap is
\(\Delta\).

With \(A_{m,k}^+\) denoting the constant-normalized coefficients,
this gives
\[
 A_{m,k}^+=s^k\prod_i
       \frac{\Gamma(ma_i+\gamma_i+(b_i-a_i)k)}
            {\Gamma(ma_i+\gamma_i)} ,
\]
where the rows are now the residual rows just described and
\(s>0\) is independent of the order. If \(a_i>0\), take
independent \(X_i\sim\mathrm{Beta}(na_i+\gamma_i,a_i)\).
The beta integral, including negative powers with finite moments,
gives
\[
 \frac{A_{n,k}^+}{\lambda^k A_{N,k}^+}
        =\mathbb E\left[
            \left(\lambda^{-1}\prod_{a_i>0}X_i^{\,b_i-a_i}\right)^k
                         \right],\qquad 0\le k\le n.           \tag{4}
\]
Rows with \(a_i=0\) make the same contribution at both orders.
All gamma arguments in the displayed range are positive.
The moments in (4) start at one. Their successive ratios
increase by Cauchy–Schwarz. Under the upper bound in (1), the
first moment is at least one, so the sequence is nondecreasing.
This is exactly \(q_k\le r_k\) for the scaled upper array.
Reversing both arrays changes \(\lambda\) to \(\lambda^{-1}\)
after a scalar normalization. Applying the same argument with
the lower bound in (1) gives \(r_k\le q_{k+1}\). This proves
(3) and hence the asserted coprimality and resultant sign.
[Linked positive Mellin quotient](MULTIVARIATE_BOTH_PRIVATE_LINKED_COEFFICIENT_ORDER.md#positive-mellin-representation-for-the-quotient)

For completeness, each shared row's logarithmic contribution to
\(A_{m,1}^+\) has real-order derivative at most
\((b_i-a_i)/m\) on \(n\le m\le N\). For an ordinary increasing
row it is a sum of \(a_i/(ma_i+h)\le1/m\) over
\(1\le h\le b_i-a_i\). For a decreasing row it is a sum of
\(-a_i/(ma_i-h)\le-1/m\) over \(0\le h<a_i-b_i\).
Every increasing private contribution is independent of \(m\).
The decreasing private contributions sum to at most
\(-(P_a-1)/m\), including the residual slots of fractional
shape, each with derivative \(-1/(m+\gamma_i-1)\le-1/m\).
Consequently
\[
      \frac{A_{N,1}^+}{A_{n,1}^+}
                    \le (N/n)^{\Delta_{\rm sh}-P_a+1},
      \qquad
      \frac{A_{N,1}^-}{A_{n,1}^-}
                    \le (N/n)^{-\Delta_{\rm sh}-P_b+1}.
\]
These are the two inclusions required for (2). Since
\(\Delta=\Delta_{\rm sh}+P_b-P_a\), (2) contains the
usual compensation, and it contains one precisely when
\(1-P_a\le\Delta\le P_b-1\).

For a fixed profile this degree condition is also necessary
for (1) to contain one at every terminal. The first coefficients
are positive rational functions of the order, asymptotic to
positive constants times \(m^{\Delta_{\rm sh}-P_a+1}\) and
\(m^{-\Delta_{\rm sh}-P_b+1}\), respectively. This follows
by counting the linear factors in the first-coefficient
products above. If either exponent is positive its rational
function is eventually strictly increasing, so the corresponding
endpoint of (1) eventually excludes one.

If two roots share a ray, choosing \(\lambda\) to be the
upper radius divided by the lower radius converts the contact
to a common root. Hence all such contacts are excluded
whenever their relative radii lie in (1).

## Reciprocal arrays and circle interlacing

If \(A,B\) are also reciprocal in the same variable, then
\(c=d\) and every descended polynomial is reciprocal.
Normalize the polynomials of the descended chain to be monic.
The recurrence takes the form
\[
                  P_{j+1}=(1+z)P_j-\eta_j zP_{j-1},
                   \qquad \eta_j>0 .
\]
Writing \(P_j=z^{j/2}p_j(\sqrt z+1/\sqrt z)\) gives
\(p_{j+1}=x p_j-\eta_j p_{j-1}\), with \(p_0=1\) and \(p_1=x\).
The finite real symmetric tridiagonal matrices with
off-diagonals \(\sqrt{\eta_j}\) give simple real roots and
strict principal-submatrix interlacing. The strictness follows
also by propagating a hypothetical common eigenvalue through
the recurrence down to the constant. None of the real roots
can have \(|x|\ge2\), since such a value comes from a nonzero
real choice of \(\sqrt z\) and \(P_j(z)>0\) at positive real
\(z\). Mapping \(x\in(-2,2)\) to
\(z=e^{2i\arccos(x/2)}\) gives strict circle interlacing.
In particular this yields a single proof of all-order angular
interlacing for ordinary mirrored both-private profiles with
any private exponents.

The [binomial-weighting obstruction](MULTIVARIATE_BINOMIAL_RATIO_DESCENT_OBSTRUCTION.md)
shows that interlaced ratios and positive Mellin quotients of the
stripped arrays cannot by themselves extend the proof to the
case where a private factorial is unavailable for cancellation.

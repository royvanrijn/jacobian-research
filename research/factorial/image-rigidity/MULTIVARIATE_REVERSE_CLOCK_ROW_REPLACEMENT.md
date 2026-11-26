# A row replacement that reverses the adjacent heat shifts

Use the [fixed-terminal compactified reference](MULTIVARIATE_FIXED_TERMINAL_COMPACTIFICATION.md#exact-compactified-data)
at consecutive orders \(n,N=n+1\). Suppose there is a positive
interior Gaussian contact at \(T_*>0,\ |\lambda_*|<T_*\), with a
double root in the lower order and a triple root at the same
nonreal physical coordinate in the upper order.

**Theorem.** If the retained core contains infinitely many
unequal ordinary rows, the contact is a limit of simple
transverse interior contacts of finite ordinary cores at the
same orders and total degree difference \(\Delta_0\).
The Gaussian strengths, biases, normalized arrays and relative
physical coordinates converge to the original ones. Any strict
total curvature cap persists. With the
[reflected heat unfolding](MULTIVARIATE_MULTIPLE_ROOT_INTERIOR_HEAT_UNFOLDING.md),
this covers all four interior multiplicity pairs in
\(\{2,3\}^2\) for infinite retained cores.

## Uniform compensated row expansion

For a shared row with midpoint \(M\) and signed gap \(d\ne0\),
write \(jL_j(k/j)\) for its endpoint-normalized log multiplier.
Set \(w=d^2/M\) and
\[
 \phi(s)=s(1-s)/2,\quad
 \psi(s)=s(1-s)(2s-1)/12,\quad
 \Xi(s)=((s-\tfrac12)^4-1/16)/12.
\]
Decreasing Gaussian strength by \(w\) on addition of one row
gives the compensated log contribution \(j(L_j+w\phi)\).
The [third-order chord expansion](MULTIVARIATE_THIRD_ORDER_BOUNDARY_AND_INTERIOR_REDUCTION.md#exact-row-expansions-and-the-projected-coefficients)
takes the following form with a uniform remainder:
\[
\begin{split}
 j(L_j+w\phi)(s)
  ={}&j\,\frac{d^3}{M^2}\psi(s)+\frac{d^2}{2M^2}\phi(s)\\
     &+j\,\frac{d^4}{M^3}\Xi(s)-\frac{d^3}{M^3}\psi(s)
                   -\frac{d^2}{6jM^3}\phi(s)
       +O_j\left(\frac{|d|^5}{M^4}\right).                   \tag{1}
\end{split}
\]
Indeed, put \(X=M+d(s-\tfrac12)\) in the exact chord second
derivative \(d^2 j\psi_1(1+jX)\), and expand the trigamma
expression \(X^{-1}-(2jX^2)^{-1}+(6j^2X^3)^{-1}+O_j(X^{-5})\)
through order \(M^{-3}\). If \(|d|/M\le1/2\), each expanded
remainder, after multiplication by \(d^2\), is bounded by
\(O_j(|d|^5/M^4)\) since \(|d|\ge1\).
Integrating with zero endpoint values gives (1).

## Two ordinary replacement rows

Choose a sequence of selected retained rows with \(w=d^2/M\to0\),
which exists by summability and infinitude. Their midpoints tend
to infinity, and they are shared for large \(M\). Fix \(K=2N\) and
define
\[
  \beta_\pm=\frac{K^3(K\pm1)}2,\qquad
  M_+=\beta_+M+K^2,\quad M_-=\beta_-M,\qquad
  d_+=\beta_+d/K,\quad d_-=-\beta_-d/K.
\]
Replace the selected row by the two rows with midpoints \(M_\pm\)
and gaps \(d_\pm\). These are ordinary rows since \(\beta_\pm\)
and \(\beta_\pm/K\) are both even integers. The actual Gaussian
strength is \(T-\Delta w\) at effective strength \(T\), where
\(\Delta w=\sum_\pm d_\pm^2/M_\pm-d^2/M=(K^2-1)w+o(w)\).

Write \(\Delta\) for the sum of the two new values minus the old
one, \(s_+=K^2,\ s_-=0\), and \(q_\pm=d_\pm/(\beta_\pm d)=\pm1/K\).
The two leading nonlinear terms are matched using
\[
 \sum_\pm\beta_\pm q_\pm^3=\sum_\pm\beta_\pm q_\pm^4=1,\qquad
 \sum_\pm q_\pm^3=0,\qquad 2\sum_\pm s_\pm q_\pm^3=1/N.
\]
In fact
\[
\begin{split}
 \Delta(d^3/M^2)&=-\frac{d^3}{NM^3}+O_N(|d|^3/M^4),&
 \Delta(d^4/M^3)&=O_N(|d|^4/M^4),\\
 \Delta(d^3/M^3)&=-\frac{d^3}{M^3}+O_N(|d|^3/M^4),&
 \Delta(d^2/M^2)&=-A+O_N(d^2/M^3),\\[-2pt]
 A&=(1-2/K^2)d^2/M^2>0.
\end{split}
\]
The aggregate coefficient of \(\psi\) in (1) is
\(O_N(|d|^3/M^4)\) in the upper order, because of the
particular shift \(K^2\) in \(M_+\). It is \(O_N(|d|^3/M^3)\)
in the lower order. The coefficient of \(\phi\) in the
compensated difference is \(-A/2+O_N(d^2/M^3)\).
All other terms are \(o(A^{m_j/2})\) on the fixed frequency
sets, for \(m_n=2,\ m_N=3\). In particular
\(|d|^5/M^4=O_N(A^{3/2}w)=o(A^{3/2})\), including when
the selected gaps diverge.

## Transverse crossing and physical scale

Take \(T=T_*+Ay\) and dilate the displacement of the
logarithmic root coordinate by \(\sqrt A\), as in the
[reflected heat proof](MULTIVARIATE_MULTIPLE_ROOT_INTERIOR_HEAT_UNFOLDING.md#exact-reflected-row-compensation).
The \(\phi\) coefficient of the compensated difference gives,
up to a scalar, a quadratic heat increment
\((A/(4j^2)+o(A))\partial_x^2\) in the logarithmic variable \(x\),
since \(\phi(k/j)=1/8-(k-j/2)^2/(2j^2)\).
The change in \(T\) gives \(Ay\partial_x^2/(2j)\).
Thus the rescaled local limiting polynomials, up to nonzero
constants, are
\[
     e^{a_j(y)\partial_X^2}X^{m_j},
                \qquad a_j(y)=\frac{y+1/(2j)}{2j}.            \tag{3}
\]
Here the quadratic term is kept exactly in the coefficient
multiplier. Its exponential can be truncated through any power
exceeding \(m_j\) on the fixed frequency set before dilation.
After dilation the coefficient of \(\partial_X^2\) tends to
\(a_j(y)\), and the rescaled root jet and its derivatives
converge to a nonzero multiple of \(X^{m_j}\). The other log
terms and the exponential truncation give \(o(A^{m_j/2})\)
coefficient errors before dilation. Truncate the retained
core, including the selected row before replacement, with
coefficient errors \(o(A^{3/2})\) and center error tending
to zero. Summability permits these truncations.

For \(y>-1/(2N)\) the positive imaginary displacements in
(3) are
\[
       \sqrt{\frac{y+1/(2n)}{n}},\qquad
       \sqrt{\frac{3(y+1/(2N))}{N}}.
\]
The first square is larger at \(y=-1/(2N)\). The
difference of squares has strictly negative slope \(1/n-3/N\).
The two positive simple branches therefore meet transversely
at a unique point of \(y>-1/(2N)\). The rescaled exact
branches and their first derivatives converge near this
crossing by Rouché continuation and Cauchy differentiation.
Hence the exact angle difference has a transverse zero there
for all sufficiently far terms of the sequence.

At the angular zero take the exact physical bias
\[
    \lambda_{\rm new}=2nN\bigl(\Gamma_{\rm fin}
                      -\log|z_{N,\rm new}/z_{n,\rm new}|\bigr).
\]
Here \(\Gamma_{\rm fin}\) includes the preserved total
\(\Delta_0\). Each replaced row's contribution to
\(\Gamma_{\rm fin}\) tends to zero by the
[compensated center estimate](MULTIVARIATE_FIXED_TERMINAL_COMPACTIFICATION.md#uniform-estimates-for-rows-of-large-size).
For each finite core the residual Gaussian signed gap is
assigned the difference between \(\Delta_0\) and the finite
core's signed gap. Since \(T-\Delta w\to T_*\) and
\(\lambda_{\rm new}\to\lambda_*\), the strict interior
inequality holds eventually. The angular zero and this
positive modulus match give the claimed physical contact.
The replaced curvatures all tend to zero.

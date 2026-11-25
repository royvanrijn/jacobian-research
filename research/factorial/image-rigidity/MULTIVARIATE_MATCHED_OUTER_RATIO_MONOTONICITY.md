# A matched outer-ratio inequality for factorial rows

Let \(n\ge2\) and let two ordinary exponent columns be distinct. On
each unequal row write \(\ell\ge0\) for the smaller exponent and
\(d>0\) for the gap. Let \(u,w\) be the reflected outer means at
orders \(n,n+1\), after positive endpoint normalization and removal
of the binomial factors. Each row contributes
\[
 u_j^2=\prod_{h=1}^d\frac{j\ell+h}{j\ell+(j-1)d+h}
\]
to the squared mean at order \(j\), and \(u,w\) are the products of
these row means at the respective orders.

**Lemma.** For \(R=n^2/(n^2-1)\), one has \(w<u^R\). In particular
\(0<u<1\) and
\[
 \frac{1-w}{1-u}>\frac{1-u^R}{1-u}
     \ge R\left(1-\frac{(R-1)(1-u)}{2u}\right).
\]
The proof also gives complete monotonicity of the rescaled
single-row loss as a function of inverse order on \(0<y\le1/2\).

For a single row, put \(y=1/j\) and \(D(h,s)=\ell+ds+y(h-ds)\).
Integrating each logarithm gives
\[
 -\log u_j = \frac{d(1-y)}2 J(y),\qquad
 J(y)=\sum_{h=1}^d\int_0^1 D(h,s)^{-1}\,ds .
\]
The denominators are bounded away from zero on compact
subintervals of \(0<y\le1/2\), even when \(\ell=0\), since
\(D(h,s)=\ell+(1-y)ds+yh\). Differentiation under the integral gives
\[
 \frac{(-1)^m J^{(m)}(y)}{m!}
    =\sum_{h=1}^d\int_0^1\frac{(h-ds)^m}{D(h,s)^{m+1}}\,ds .
\]
For even \(m\) the terms are nonnegative. For odd \(m\), split
the interval into \(s=(r-b)/d\) with \(1\le r\le d\) and
\(b\in[0,1]\), and pair the \((h,r)\) term with the \((r,h)\)
term at the same \(b\). When \(h\ge r\), put \(a=h-r\).
The denominators \(D_1,D_2\) before and after the swap satisfy
\[
 D_2-D_1=(1-2y)a\ge0 .
\]
The numerators are \((a+b)^m\) and \((b-a)^m\). Their sum
after division by the respective denominators is nonnegative:
if the second numerator is negative then
\((a+b)^m\ge |b-a|^m\) and \(D_1\le D_2\).
The diagonal terms with \(h=r\) are strictly positive for
\(b>0\). In particular \(J\) is strictly decreasing, and
\[
 \frac{1-1/(n+1)}{1-1/n}=R
\]
gives the strict log-ratio inequality at each row.
Summing the row logarithms gives \(-\log w>R(-\log u)\).

The diagonal pairing also supplies a quantitative row bound.
At \(m=1\), the pair sum is
\(a(D_1^{-2}-D_2^{-2})+b(D_1^{-2}+D_2^{-2})\).
Thus \(-J'\) is at least the sum with integrands \(bD^{-2}\)
over all split intervals. At fixed \(h,r\), the factor \(D^{-2}\)
increases with \(b\). Integrating \((b-1/2)D^{-2}\) is
nonnegative by pairing \(b\) with \(1-b\). Cauchy–Schwarz with
total measure \(d\) now gives
\[
          -J'(y)\ge \frac12\sum_{h=1}^d\int_0^1 D(h,s)^{-2}\,ds
                   \ge \frac{J(y)^2}{2d}.
\]
Integration of the derivative of \(1/J\) yields the stronger
row inequality
\[
 -\log u_{n+1}\ge
 R\,\frac{-\log u_n}{1-(-\log u_n)/(d^2(n^2-1))}.
\]
The denominator is positive, as follows from the same integrated
inequality and \(1/J(1/(n+1))>0\).

Finally, Taylor integration gives
\[
 1-u^R=R(1-u)-R(R-1)\int_u^1(s-u)s^{R-2}\,ds .
\]
Here \(1<R\le2\) and \(s^{R-2}\le1/u\) on the integral, which
proves the deficit bound. It strengthens the lower log-ratio
bound in the [matched terminal-ten comparison](MULTIVARIATE_TERMINAL_TEN_PRIVATE_DEGREES.md#shared-outer-mean-boxes)
and also applies to private rows.

## Coupling to the lower-order tanh moments

Write \(M_i=-\log u_{n,i}\), \(U=\sum_i M_i>0\), and
\(\rho=(-\log w)/U\). On row \(i\) the lower-order factors have
\[
 T_{i,h}=\frac{(n-1)d_i}{2(n\ell_i+h)+(n-1)d_i},
 \qquad 1\le h\le d_i.
\]
Let \(V_2=\sum_{i,h} T_{i,h}^2\) and put
\(P=\tfrac12\sum_{i,h}\epsilon_i T_{i,h}^2\) for any row signs
\(\epsilon_i\in\{-1,1\}\). With \(\kappa=(n-1)/(n+1)\), one has
\[
 \rho-R\ge
 \frac{R\kappa^2 V_2^2}{(n^2-1)U^3}
 \ge\frac{4R\kappa^2P^2}{(n^2-1)U^3}.
\]
This relates the extra outer loss in the next order to the same
signed second moment that enters the
[matched tanh transport](MULTIVARIATE_MATCHED_TANH_MOMENTS.md).

In fact \(T_{i,\min}\ge\kappa T_{i,\max}\) by the exact denominators,
and \(M_i=\sum_h\operatorname{artanh}T_{i,h}\ge d_i\kappa
T_{i,\max}\). Consequently
\[
 V_2\le\sum_i M_i T_{i,\max}
      \le \kappa^{-1}\sum_i M_i(M_i/d_i)
      \le\kappa^{-1}U\left(\sum_i (M_i/d_i)^2\right)^{1/2}.
\]
The last step is Cauchy–Schwarz and positivity of the \(M_i\).
Summing the quantitative row inequality with
\((1-x)^{-1}\ge1+x\) on \(0\le x<1\) gives
\[
 \rho-R\ge\frac{R}{(n^2-1)U}\sum_i(M_i/d_i)^2.
\]
Combining these estimates and \(2|P|\le V_2\) proves both bounds.

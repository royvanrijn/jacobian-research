# Paired factorial differences give effective common-factor norm heights

Fix adjacent ordinary moment orders \(n,N=n+1\), \(n\ge2\),
and split the ordinary profile into an arbitrary core and a
nonempty shared tail of unequal rows. For tail row \(i\) write
\(m_i=(a_i+b_i)/2\), \(d_i=b_i-a_i\), and put
\[
 \eta_i=|d_i|/m_i,\qquad \eta=\max_i\eta_i,\qquad
 W=\sum_i d_i^2/m_i,\qquad Z=\sum_i \eta_i^2 .
\]
Assume \(\eta\le 1/2\) and \(W\le W_{\max}\), where
\(W_{\max}>0\) is fixed. There are no bounds on the number
of rows or the exponents in the core. Let \(H\) be a monic
rational common factor of degree \(b>0\) in the original
physical variable and \(c\) the product of its roots.
The even parent has no real roots, so \(c>0\).

For every integer \(2\le r\le n\) there are explicitly
specified positive rational scales \(S_+,S_-\) made from
the literal adjacent factorial blocks of the two parents,
with the following property. The positive rational numbers
\(M_\pm=cS_\pm\) satisfy the uniform, strict comparison
\[
            0<M_+-M_-\le C_{n,r,b,W_{\max}}
                      \bigl(\eta^{r-2}Z+W\eta^r\bigr).       \tag{1}
\]
In particular the reduced denominator of their difference,
and therefore also the least common denominator of \(M_+\)
and \(M_-\), is at least the reciprocal of the right side.
No algebraicity or transversality of a real limiting
contact is required. At \(r=n\), the universal estimate
\(Z\le W\eta\) gives order \(W\eta^{n-1}\) in (1);
if \(Z=O(\eta^2)\), (1) gives order \(\eta^n\).

## Comparing the two factorial grids

Write \(R_{j,k,i}=((j-k)a_i+kb_i)!/(ja_i)!\) and
\(L_{j,i}(k)=\log R_{j,k,i}\).
Use \(\Delta\) for the forward unit step in the index.
Define a discrete functional by using the two extreme
order-\(r\) windows, counting both even if they coincide:
\[
   E_{j,r,i}
     = j^{r-1}\bigl(\Delta^r L_{j,i}(0)
                          +\Delta^r L_{j,i}(j-r)\bigr).
\]
For each unequal row set \(\epsilon_i=(-\operatorname{sgn}d_i)^r\).
The first assertion is the strict inequality
\[
                 \epsilon_i(E_{N,r,i}-E_{n,r,i})>0 .         \tag{2}
\]
It permits odd \(r\) and mixed orientations because the
sign is chosen separately for each ordinary row.

To prove (2), let \(X_{j,r}\) have the distribution of
the sum of \(r\) independent uniforms on
\([-1/(2j),1/(2j)]\), plus the independent random shift
\(\pm(j-r)/(2j)\) with equal probabilities. Applying the
integral form of a forward difference \(r\) times, after
changing the continuous variable to \(k/j-1/2\), gives
\[
  \epsilon_i E_{j,r,i}
    =2|d_i|^r\int_0^\infty t^{r-1} e^{-m_it}
             B_j(t) M_{j,r}(|d_i|t)\,dt ,
 \quad B_j(t)=\frac{1}{j(e^{t/j}-1)},                          \tag{3}
\]
where
\[
 M_{j,r}(u)=
    \left(\frac{\sinh(u/(2j))}{u/(2j)}\right)^r
    \cosh\!\left(\frac{(j-r)u}{2j}\right).
\]
Indeed \(j^{r-1}\) makes each finite difference the average
of the \(r\)-th derivative of
\(j^{-1}\log\Gamma(1+j(m_i+d_i x))\).
The normalized \(r\)-fold integration is over a sum of
independent uniforms, centered as above for the first and
last windows. The derivative formula for \(r\ge2\) follows
by integrating the series for the polygamma function,
\[
 \frac{d^r}{dx^r}\log\Gamma(1+x)
       =(-1)^r\int_0^\infty
                         \frac{u^{r-1}e^{-xu}}{e^u-1}\,du .
\]
After the change of Laplace variable \(t=ju\) this gives (3).
The integrals converge uniformly on the row sampling
interval because its exponents stay positive.

For \(t>0\), \(B_j(t)\) grows strictly with \(j\):
\(y/(e^y-1)\) decreases for \(y>0\), and \(B_j(t)\) is
\(t^{-1}\) times its value at \(y=t/j\).
We also have \(M_{N,r}(u)\ge M_{n,r}(u)\).
For this, set \(x=u/(2N)\), \(h=n-r\),
\(f(x)=\log(\sinh(x)/x)\) and \(g(x)=\log\cosh x\).
The function \(f\) is convex as the log moment generating
function of a uniform variable on \([-1,1]\).
Since \(f(2x)-f(x)=g(x)\), convexity gives
\[
 r\bigl(f((1+1/n)x)-f(x)\bigr)\le (r/n)g(x).
\]
Also \(g\) is convex with \(g(0)=0\), so
\[
 \begin{split}
  g((h+1)x)-g(h(1+1/n)x)
   &\ge (1-h/n)\bigl(g((h+1)x)-g(hx)\bigr)\\
   &\ge (r/n)g(x).
 \end{split}
\]
The two inequalities give the claimed comparison of
the moment generating functions. Their strict positive
kernel factor proves (2).

## The size after cancellation of the common leading term

For \(t>0\) one has
\[
                       1/t-1/j\le B_j(t)\le 1/t .
\]
For the lower bound at \(y=t/j<1\) use
\(e^y\le(1-y)^{-1}\). At \(y\ge1\) it is immediate.
Thus \(B_N-B_n\le1/n\). Each \(X_{j,r}\) has mean zero
and lies in \([-1/2,1/2]\). Taylor's formula for the
exponential gives
\[
 |M_{j,r}(u)-1|\le (u^2/8)e^{|u|/2},
                    \qquad M_{j,r}(u)\le e^{|u|/2}.
\]
In the difference of the integrands in (3), write
\((B_N-B_n)M_{n,r}+B_N(M_{N,r}-M_{n,r})\).
The two bounds above and elementary exponential
integrals at \(m_i-|d_i|/2\ge3m_i/4\) give
\[
 0<\epsilon_i(E_{N,r,i}-E_{n,r,i})
   \le C_{n,r}\left(\eta_i^r+
                      \frac{d_i^2}{m_i}\eta_i^r\right).
\]
Summing and using \(\sum_i\eta_i^r\le\eta^{r-2}Z\) proves
the required bound on the logarithm of \(S_+/S_-\).
The contribution of \(1/t\) common to both grids is
thus measured only against the quadratic-order difference
of the centered quadrature kernels.

## Rational scales and the local test

Write \(r_{j,k,i}=R_{j,k,i}/R_{j,k-1,i}\).
Each row expression in (2) is an integer linear combination
of \(\log r_{j,k,i}\), whose coefficient sum is zero
separately at each order. This follows by taking one
forward difference first, leaving an \((r-1)\)-st
difference of adjacent logs. Write the integer coefficients
as \(e_{j,k,i}^+-e_{j,k,i}^-\), the positive and negative
parts. Let \(b_{j,i}\) be their common sum at order \(j\),
and put \(b_i=b_{n,i}+b_{N,i}\). With \(\Lambda=N/n\)
define
\[
 S_{\pm,i}=\Lambda^{-b_{N,i}d_i}\,
        r_{n,1,i}^{\,b-b_i}
                     \prod_{j,k} r_{j,k,i}^{\,e^\pm_{j,k,i}} .
\]
Take \(S_\pm=(r_{n,1}^{C})^b\prod_i S_{\pm,i}\), where
\(r_{n,1}^{C}\) is the first factorial ratio of the
core after stripping the binomial coefficient. The added
prefactors cancel in \(S_+/S_-\), whose logarithm is
the sum of (2). All exponents are integers, and
\(b_{j,i}\) is bounded in terms of \(n,r\) alone.

For the size bound, denote the endpoint scale of row
\(i\) at order \(j\) by \(D_{j,i}=R_{j,j,i}^{1/j}\).
The normalized adjacent logs
\(\log(r_{j,k,i}/D_{j,i})\) have modulus
\(O_j(d_i^2/m_i)\). Indeed they are differences of
averages of the first derivative in the normalized
continuous variable, whose oscillation is bounded by
the supremum of the second derivative. This second
derivative, after division by \(j\), is
\(j d_i^2\psi'(1+j(m_i+d_i x))\), which is
\(O(d_i^2/m_i)\) on \(|x|\le1/2\) by the elementary
trigamma integral comparison. The
[exact center bound](MULTIVARIATE_POSITIVE_DIFFUSE_EXCEPTION_SIGNED_JETS.md#a-matched-expansion-with-a-second-order-ordinary-remainder)
also gives
\(\log(D_{N,i}/D_{n,i})-d_i\log\Lambda=O_{n,N}(d_i^2/m_i)\).
The definition of \(S_{\pm,i}\) therefore gives
\[
               S_{\pm,i}/D_{n,i}^b\le
                        \exp(C_{n,r,b} d_i^2/m_i).
\]
If \(G_n^C\) is the core endpoint scale, then
\(r_{n,1}^{C}\le G_n^C\) by the convexity of the log
factorial samples, with the binomial coefficients removed.
All roots of the lower parent in the full endpoint
scale \(G_n^C\prod_iD_{n,i}\) have modulus at most \(2^n\),
since the normalized polynomial has endpoints one
and coefficients at most \(\binom nk\). Hence
\(M_\pm\le 2^{nb}\exp(C_{n,r,b}W)\).
Writing \(L=\log(S_+/S_-)>0\), use
\(M_+-M_-=M_+(1-e^{-L})\le M_+ L\).
This and \(W\le W_{\max}\) prove (1).

At each prime let \(s_1,\ldots,s_b\) be the valuations
of the proposed factor roots in the original variable,
counted with multiplicity. The exact denominator of
\(M_+-M_-\) has exponent
\[
                  \left(-\sum_\ell s_\ell
                             -v_p(S_+-S_-)\right)_+.          \tag{4}
\]
The valuation in (4) is computed from the same two
literal products and can capitalize on cancellation
between them. A weaker bound uses
\(\min(v_p(S_+),v_p(S_-))\) in its place.
Shared Newton valuations and
[truncated face compatibility](MULTIVARIATE_PHYSICAL_TRUNCATED_NEWTON_FACES.md)
let one maximize (4) over admissible allocations of
one common factor at each prime. If the product of these
caps is below the lower bound from (1), that allocation
of a common factor is excluded effectively. Unlike the
trace spread, this budget uses only the sum of the selected
root valuations, with one predetermined normalization
in each of the two norms. Every choice of shared tail subset
and every \(r\le n\) gives a separate necessary test.
A uniform upper bound for a covering family of these
tests is still required for unrestricted closure.


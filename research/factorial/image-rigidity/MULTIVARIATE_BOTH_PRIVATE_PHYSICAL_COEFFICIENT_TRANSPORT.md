# Physical coefficient transport after private factorial cancellation

Let an ordinary binomial have private exponent totals \(P_A,P_B\ge1\)
on the two sides and arbitrary shared rows. Orient the argument \(z\)
toward the \(B\) column. Write
\[
 \Delta_{\rm sh}=\sum_{\rm shared}(b_i-a_i),\qquad
 V=\sum_{\rm shared}\frac{(b_i-a_i)^2}{\min(a_i,b_i)}.
\]
At each degree \(m\) normalize the constant coefficient to one
by scalar division only. Let \(C^A_{m,k}\) be the resulting
coefficient at distance \(k\) from that endpoint, and let
\(C^B_{m,k}\) be the same coefficient in the reversed polynomial
normalized at its own constant coefficient. These coefficients
use the physical argument or its reciprocal.

## Compensating powers and derivative bounds

**Theorem.** Define
\[
       \lambda_A=\Delta_{\rm sh}-(P_A-1),\qquad
       \lambda_B=-\Delta_{\rm sh}-(P_B-1).
\]
For \(s=A,B\) the compensated coefficient
\(m^{-\lambda_s k}C^s_{m,k}\) is nonincreasing for each fixed
integer \(k\ge0\), including under real-order interpolation
at \(m\ge k\), \(m>0\). More quantitatively
\[
 0\le \frac{\lambda_s k}{m}-\partial_m\log C^s_{m,k}
 \le \frac{V(k^2+k)}{2m^2}
              +\frac{(P_s-1)k(k+1)}{2m(m-k+1/P_s)} .           \tag{1}
\]
In particular if
\[
             -(P_B-1)\le\Delta_{\rm sh}\le P_A-1,              \tag{2}
\]
the raw coefficients in both orientations are nonincreasing
with the order in the same physical variable. There is no
endpoint-radius correction in (2).

For a shared row with starting exponent \(a>0\), signed gap
\(d\), and \(\ell=\min(a,a+d)\), put \(q=|d|k\).
The log physical ratio at index \(k\) is the log of
\((ma+kd)!/(ma)!\), interpreted using Gamma for real \(m\).
If \(d>0\), its derivative is
\(\sum_{h=1}^q a/(ma+h)\). For \(d<0\) it is
\(-\sum_{h=0}^{q-1} a/(ma-h)\).
Subtracting either derivative from \(kd/m\) gives a
nonnegative sum of \(h/[m(ma+h)]\) or \(h/[m(ma-h)]\).
Both sums are at most
\((d^2k^2+|d|k)/(2m^2\ell)\). This proves the shared bound.

Reserve one private row on each side. A row of exponent \(d\)
supplies \(d\) factors of shape \(j/d\), \(1\le j\le d\),
by grouping the factorial into residue classes. The shape-one
factor in each reserved row cancels the binomial coefficient
up to the degree-only scalar \(m!\). The powers of \(d\)
contribute only a fixed argument scaling. The increasing
private factors at distance \(k\) are independent of \(m\).
On the decreasing side there are precisely \(P_s-1\) factors,
with shapes \(\gamma\in[1/P_s,1]\). Each has raw ratio
\(\Gamma(m-k+\gamma)/\Gamma(m+\gamma)\); its log derivative is
\(-\sum_{h=0}^{k-1}1/(m-h-(1-\gamma))\).
Subtracting that derivative from \(-k/m\) gives a
nonnegative sum, bounded by
\(k(k+1)/(2m(m-k+\gamma))\).
Combine the private and shared contributions to prove (1).

## A positive moment quotient

For any adjacent \(n,N=n+1\), and either orientation \(s\), the ratios
\[
         (N/n)^{\lambda_s k}\frac{C^s_{n,k}}{C^s_{N,k}},
                         \qquad 0\le k\le n,                   \tag{3}
\]
are the moments of one positive random variable, as are
the unadjusted ratios. In the Gamma factors above, write
the positive starting shape as \(na+\gamma\) and signed
slope as \(d\). The factor's ratio across the two orders
is the moment \(\mathbb E X^{dk}\) for
\(X\) of Beta shape parameters \(na+\gamma,a\).
Its moments exist throughout \(0\le k\le n\), since every
initial Gamma argument stays positive up to that index.
Factors with \(a=0\) are independent of the order.
Taking independent products and applying the deterministic
power in (3) proves the claim. In particular (3) is
log-convex in \(k\), and it is nondecreasing since its
first moment is at least one by (1).
This is a physical-variable complement to the
[endpoint-normalized comparison](MULTIVARIATE_BOTH_PRIVATE_LINKED_COEFFICIENT_ORDER.md).

## The smallest mixed private totals

For \(P_A=1\), \(P_B=2\), the signed total degree difference is
\(\Delta=\Delta_{\rm sh}+1\), and
\(\lambda_A=\Delta-1,\ \lambda_B=-\Delta\).
In the signed balance \(\Delta=0\), the \(A\)-oriented log
derivative lies between
\(-k/m-V(k^2+k)/(2m^2)\) and \(-k/m\).
The \(B\)-oriented log derivative loses its \(1/m\) term
exactly; its magnitude is at most
\[
       \frac{V(k^2+k)}{2m^2}
       +\begin{cases}
            \dfrac{k^2}{2m(m-k+1/2)},&
                          \text{one quadratic private row},\\[4pt]
            \dfrac{k(k-1)}{2m(m-k+1)},&
                          \text{two linear private rows}.
         \end{cases}
\]
These last quantities use
\(\sum_{h=0}^{k-1}(h+1-\gamma)=k(k+1-2\gamma)/2\).
Thus the private center drift and the corresponding
part of the right-end coefficient drift cancel before
the actual-root phase comparison. The sign and magnitude
of the remaining signed integer drift must still be
controlled by the modulus equations, and the coefficient
differences alone do not settle the alternating phases.

# Quantitative law of the primitive quotient of a full triple

Use the fixed-order and comparability hypotheses of the
[full-triple expansion](POWERFUL_TRIPLES_COMMON_DIVISOR_BOUND.md#asymptotic-expansion-for-comparable-orders).
Write \(M=\max(s_A,s_B,s_C)\), and let \(q=(u,v,w)\) run over
ordered positive primitive triples \(u+v=w\). Let \(a_q(X)\) count the
triples with quotient \(q\) after division by their common divisor.
The total count is \(A(X)=\sum_q a_q(X)\).

## Summation in total variation

For every admissible truncation \(J\), define
\[
 \mu_j(q)=b_j(\operatorname{rad}(uvw))\,w^{-1/(M+j)}
             \sum_{f\ {\rm allowed}} f^{-1/(M+j)}
 \qquad(0\le j\le J),
\]
using the same factors as in formula (PC). The sign of \(\mu_j(q)\)
is \((-1)^j\), and \(\sum_q \mu_j(q)=C_{\boldsymbol s,j}\).
In fact
\[
 \sum_q\left|a_q(X)-\sum_{j=0}^J
            \mu_j(q) X^{1/(M+j)}\right|
 \ll_{\boldsymbol s,J,\varepsilon} X^{1/(M+J+1)+\varepsilon}.
\]
To see this, apply the coprimality-uniform single-integer expansion
separately to each allowed \(f\) and each \(q\). Its pointwise remainder
is bounded by
\[
  O_{\boldsymbol s,J,\sigma}
     \left(2^{(J+1)\omega(\operatorname{rad}(uvw))}
                      (X/(wf))^\sigma\right)
\]
for \(\sigma>1/(M+J+1)\) sufficiently close to that endpoint.
The [dyadic estimate for summing the common factors](POWERFUL_TRIPLES_COMMON_DIVISOR_BOUND.md#summing-the-quotient-and-common-factors)
proves that the sum of these majorants is \(O(X^\sigma)\).
The triangle inequality applies before any sum over quotients.

In particular \(C_{\boldsymbol s,0}>0\), and the quotient of a
uniformly chosen ordered full triple with \(c\le X\) has law
\(\mathbb P_X(q)=a_q(X)/A(X)\). Its limiting distribution is
\(\mathbb P_\infty(q)=\mu_0(q)/C_{\boldsymbol s,0}\), with
\[
 \|\mathbb P_X-\mathbb P_\infty\|_{\rm TV}
 \ll_{\boldsymbol s,\varepsilon} X^{-1/(M(M+1))+\varepsilon}.
\]
Indeed the \(J=0\) remainder is bounded in \(\ell^1\) by
\(O(X^{1/(M+1)+\varepsilon})\), the scalar normalizing error has the
same bound, and \(A(X)\asymp_{\boldsymbol s}X^{1/M}\).

## Uniform tail in the primitive height

The dyadic summability proof has a fixed negative exponent of the
height scale \(W\), uniformly in each of its two ranges of the
deficient factor \(H\). Consequently, restricting the majorant sum
to \(w>Y\) gives \(O_{\boldsymbol s,J,\sigma}(Y^{-\delta})\) times
\(X^\sigma\) for some \(\delta>0\). The same tail estimate holds for the
coefficient sums since \(1/(M+j)>\sigma\). With \(J=0\) and
\(A(X)\asymp X^{1/M}\) this gives
\[
 \mathbb P_X(w>Y)+\mathbb P_\infty(w>Y)
       \ll_{\boldsymbol s} Y^{-\delta}
\]
uniformly for all sufficiently large \(X\) and \(Y\ge1\).

# Uniform endpoint expansions at growing curvature

Let \(N\to\infty\) along any sequence of ordinary factorial
profiles, in arbitrary finite dimensions. Delete equal rows
and, for the remaining exponent columns \(a,b\), put
\[
 S=\sum_i\frac{(b_i-a_i)^2}{\max(a_i,b_i)},\qquad
 V=\sum_{a_i b_i>0}\frac{(b_i-a_i)^2}{\min(a_i,b_i)}
       +\sum_{a_i b_i=0}\max(a_i,b_i),\qquad
 D=N e^{-S/2}.
\]
If
\[
 (D+\log N)\sqrt{\frac{1+V}{N}}\le10^{-6},             \tag{1}
\]
then the orders \(N-1\) and \(N\) are coprime at every
sufficiently large term of the sequence. This covers all
exponent profiles, including profiles with private variables
and profiles of equal total degrees.

In particular, the conclusion is uniform over profiles with
\[
                      S\ge\log N+2\log\log N+30.
\]
Indeed, \(V\le S(S+2)\). To see this, in a shared row put
\(d=|b_i-a_i|\), \(l=\min(a_i,b_i)\). Its contribution to
\(S\) is \(d^2/(l+d)\ge d/l-1\), since \(l\ge1\). Thus its
contribution to \(V\) is at most \(S+2\) times its
contribution to \(S\). A private row contributes equally
to \(S\) and \(V\). The
[large-curvature criterion](MULTIVARIATE_HIGH_CURVATURE_ANGULAR_SLOPE.md)
already covers \(S\ge2\log N+5\). On the remaining band,
we have \(D\le\sqrt N e^{-15}/\log N\) and
\(\sqrt{1+V}\le S+1\le2\log N+6\), which imply (1)
for large \(N\).
The proof below works on this remaining band, so throughout
it we assume \(S<2\log N+5\). Condition (1) gives
\(D\le10^{-6}\sqrt N\) and hence \(S\to\infty\).

## Normalization and endpoint halves

For \(m=N-1,N\), write
\[
 C_{m,k}=\prod_i((m-k)a_i+kb_i)!,\qquad
 G_m=(C_{m,m}/C_{m,0})^{1/m},\qquad
 c_{m,k}=\binom mk\frac{C_{m,k}}{C_{m,0}G_m^k}.
\]
The original moment polynomial divided by its constant
coefficient has the form \(T_m(G_m t)\), where
\(T_m(z)=\sum_{k=0}^m c_{m,k}z^k\). Its normalized endpoint
coefficients are both one. We use the
[discrete curvature estimate](MULTIVARIATE_HIGH_CURVATURE_ANGULAR_SLOPE.md#an-endpoint-coefficient-estimate)
\[
 c_{m,k}\le\binom mk\exp\left(-\frac{S k(m-k)}{2m}\right).
                                                               \tag{2}
\]
Put \(h=\lfloor m/2\rfloor\) and
\[
 L_m(w)=\sum_{k=0}^h c_{m,k}w^k,\qquad
 U_m(w)=\sum_{k=0}^{m-h-1}c_{m,m-k}w^k,\qquad
 T_m(z)=L_m(z)+z^m U_m(1/z).
\]
The right half is the left half in the opposite orientation,
with a possible shift of one in the final summation index.
For either half \(H_m\), define \(B_m=c_{m,1}\) in its
oriented system when all starting exponents are positive,
and define \(B_m=0\) when a starting exponent is zero.
We will show that
\[
 H_m(w)=e^{B_m w}(1+E_m(w)),\qquad
 \sup_{|w|\le2}|E_m(w)|<\frac1{1000},\qquad 0\le B_m\le 2D,         \tag{3}
\]
uniformly for both orders and both halves. When \(B_m>0\),
its bound follows from (2) at \(k=1\) since
\(m e^{-S/2+S/(2m)}\le2D\) for large \(N\).

Suppose first that a starting exponent is zero and the
corresponding ending exponent is \(e\ge1\). Write
\(f_{m,k}=k!/(m!)^{k/m}\). The normalized factorial from
that row is at most \(f_{m,k}^e\). Indeed
\((ek)!/(k!)^e\) starts at one and is log-convex in \(k\),
since its consecutive ratios are
\(\prod_{j=1}^e(e-(e-j)/(k+1))\).
Moreover \(\binom mk f_{m,k}\le\exp(k)\) by
\(m!\ge(m/\exp(1))^m\).
The single-row curvature bound gives
\(f_{m,k}\le\exp(-k(m-k)/(2m))\). In (2) the other rows
contribute curvature \(S-e\). Therefore, for \(k\le m/2\),
\[
 c_{m,k}\le
 \exp\left(k-\frac{(S-1)k(m-k)}{2m}\right)
 \le \exp\left(-\frac{(S-5)k}{4}\right).
\]
On \(|w|\le2\), a geometric sum gives \(H_m(w)=1+o(1)\),
as required.

## Analytic interpolation at a positive endpoint

Assume now \(a_i\ge1\) at the oriented starting endpoint.
With \(d_i=b_i-a_i\) signed, the coefficients have the
analytic interpolation
\[
 \frac{c_{m,k} k!}{B_m^k}=e^{f_m(k)},\qquad
 \begin{aligned}
 f_m(u)={}&\log\Gamma(m+1)-\log\Gamma(m+1-u)-u\log m\\
 &+\sum_i\bigl[\log\Gamma(ma_i+1+d_i u)-\log\Gamma(ma_i+1)\\
 &\hspace{32mm}
   -u\bigl(\log\Gamma(ma_i+1+d_i)-\log\Gamma(ma_i+1)\bigr)\bigr].
 \end{aligned}                                               \tag{4}
\]
In particular \(f_m(0)=f_m(1)=0\).
Let \(R=\tfrac14\sqrt{m/(1+V)}\). Since \(d_i^2/a_i\le V\),
for \(|u|\le R\) we have
\(|d_i u|\le\tfrac14\sqrt{ma_i}\le ma_i/4\). Also
\(R\le m/4\). All Gamma arguments have positive real part,
so use the analytic logarithm that is real on the
positive axis. The series
\(\psi_1(z)=\sum_{j\ge0}(z+j)^{-2}\) bounds its second
derivative there. Using
\(\sum_{j\ge0}(x+j)^{-2}\le x^{-2}+x^{-1}\) for \(x>0\),
we obtain
\[
 |f_m''(u)|\le\frac{4}{m}\left(1+\sum_i\frac{d_i^2}{a_i}\right)
             \le\frac{4(1+V)}m .
\]
Condition (1) gives \(R\to\infty\). The zeros at \(0,1\)
give \(|f_m'(0)|\le2(1+V)/m\), and integration on the
segment from \(0\) to \(u\) gives
\[
 |f_m(u)|\le\frac{2(1+V)}m(|u|+|u|^2)\le\frac14.
\]
Consequently the coefficients \(\alpha_j\) of
\(e^{f_m(u)}=\sum_{j\ge0}\alpha_j u^j\) satisfy
\(\alpha_0=1\) and \(|\alpha_j|\le2 R^{-j}\).

We next sum the interpolated coefficients with a relative error
on the whole disk \(|w|\le2\). Put
\[
 K=\lceil128(D+\log N)\rceil,\qquad
 J=\lceil16(D+\log N)\rceil .
\]
Write \(\epsilon=10^{-6}\). Since \(m\ge N/2\), (1) gives
\((D+\log N)/R\le4\sqrt2\,\epsilon\). For large \(N\),
\(K\le129(D+\log N)\) and \(J\le17(D+\log N)\), so
\[
               K/R<1024\epsilon,\qquad
               (2B_m+J)/R<128\epsilon.
\]
In particular \(K\) is eventually less than the degree of
either half. The terms with \(k>K\) in the
original polynomial have a negligible relative contribution.
Indeed \(\binom mk\le(em/k)^k\), with natural base \(e\),
and (2) bound \(k^{-1}\log c_{m,k}\) by
\(1-\log w-(S/2)(1-w)\) at \(w=k/m\). This function is
convex on \([K/m,1/2]\). At \(K/m\) it is at most
\(\log(eD/K)+SK/(2m)<-\log8\) for large \(N\), since
\(SK/m=o(1)\). At \(1/2\) it tends to \(-\infty\).
Thus \(c_{m,k}\le8^{-k}\) for \(K\le k\le m/2\), in
either orientation. Dividing the tail bound on \(|w|\le2\)
by \(|e^{B_m w}|\) gives \(O(4^{-K}e^{4D})=o(1)\).

For \(k\le K\) replace \(e^{f_m(k)}\) by
\(P_J(k)=\sum_{j=0}^J\alpha_j k^j\). The pointwise
error is at most \(4(K/R)^{J+1}\). The error after
multiplication by \((B_m w)^k/k!\), summation, and
division by \(e^{B_m w}\) is at most
\[
                  4(K/R)^{J+1}e^{8D}=o(1).
\]
This follows already from \(K/R\le1/2\) and
\(16\log2>8\). We may extend the sum for \(P_J(k)\)
over all \(k\ge0\). For \(j\le J\) and \(k\ge0\),
\((k/R)^j\le e^k\) for large \(N\), by maximizing
\(k^j e^{-k}\). Hence \(|P_J(k)|\le2(J+1)e^k\).
The added tail, divided by \(|e^{B_m w}|\), is at most
\[
 2(J+1)e^{-K+2(e^2+1)B_m}=o(1).
\]
Here the estimate uses \(1\le e^{k-K}\) for \(k>K\), and
sums the resulting exponential series.

For the full exponential sum, the identity
\[
 \sum_{k\ge0}\frac{u^k k^j}{k!}=e^u\,\mathcal T_j(u),\qquad
 \mathcal T_j(u)=\sum_{\ell=0}^j
                   \left\{\begin{smallmatrix}j\\\ell\end{smallmatrix}\right\}
                   u^\ell
\]
follows by expanding \(k^j\) in falling factorials.
The Stirling number of the second kind in this formula
is at most \(\binom j\ell\ell^{j-\ell}\), by choosing
the least representative in each block of a partition and
assigning the remaining elements to blocks. It follows that
\(|\mathcal T_j(u)|\le(|u|+j)^j\). Our full sum divided
by \(e^{B_m w}\) differs from one by less than \(257\epsilon\), because
\[
 \sum_{j=1}^J|\alpha_j\mathcal T_j(B_m w)|
 \le 2\sum_{j=1}^J\left(\frac{2B_m+J}{R}\right)^j<257\epsilon.
\]
Together with the \(o(1)\) tail errors, this proves (3).

## Comparison of the orders

The endpoint coefficients in (3) satisfy
\[
                    B_N-B_{N-1}=o(1)                     \tag{5}
\]
in each orientation. This is immediate for an endpoint
with a zero exponent. At a positive endpoint, the
factorial ratios in \(C_{m,1}/C_{m,0}\) give
\[
 \log(C_{m,1}/C_{m,0})
     =\sum_i d_i\log(ma_i)+O(V/m).
\]
Indeed each row is a product, or reciprocal product, of
\(|d_i|\) factors differing from \(ma_i\) by at most
\(|d_i|\). Their relative changes are eventually below
\(1/2\), using \(d_i^2/a_i\le V=O(\log^2 N)\), so the
logarithmic errors total \(O(\sum_i d_i^2/(ma_i))\).

For the endpoint center, set \(F(x)=x\log x-x\), with
\(F(0)=0\), and \(E(x)=\log\Gamma(x+1)-F(x)\) for
\(x\ge0\). For \(x>0\), convexity of the logarithm
of Gamma and its consecutive differences give
\(0\le E'(x)\le1/x\). Thus a shared row has
\(|E(mb_i)-E(ma_i)|\le|\log(b_i/a_i)|\), and a private
row of degree \(e\) contributes at most \(1+\log(me)\),
by \(E(1)=1\) and \(E(0)=0\). Summing yields, with
\(\Delta=\sum_i d_i\),
\[
 \log G_m=\Delta\log m+\sum_i(F(b_i)-F(a_i))
             +O\left(\frac{V+(2+\log m)S}{m}\right).
\]
Consequently \(B_m/m\) differs from one common positive
scale for the two orders by a relative error
\(O(\log^2 N/N)\). Since \(B_m\le2D\) and
\(D\le10^{-6}\sqrt N\), this proves (5).

By (3), each half is zero-free on \(|w|\le2\) for
large \(N\). Its logarithm real on the positive radius is
\(B_m w+\log(1+E_m(w))\). On \(|w|\le3/2\) the
derivative of the second term has absolute value at
most \(2/999\), by the bound
\(|\log(1+E_m)|\le1/999\) and Cauchy's estimate from the
disk of radius two.

Set \(r_*=\exp(16(D+1)/m)\), which is eventually below
\(3/2\). On \(|z|=r_*\), (3) and \(B_m\le2D\) give
\[
 |L_m(z)|\le e^{3D+1/999},\qquad
 |z^m U_m(1/z)|\ge e^{16(D+1)-2D-1/999}.
\]
The dominating term has all its zeros in \(|z|<r_*\)
because \(U_m\) is zero-free on \(|w|\le2\).
Rouché's theorem for \(T_m\) and then its reciprocal
puts every root radius \(r\) in \((r_*^{-1},r_*)\).

Consider a root of angle \(0<\theta\le\pi\).
The exact root equation and the two half logarithms give
\[
 \begin{aligned}
 m\theta&=(2j_m+1)\pi+\delta_m,\qquad j_m\in\mathbb Z,\\
 \delta_m&=\bigl(B_m^{L}r+B_m^{U}/r\bigr)\sin\theta
                  +\xi_m,\qquad |\xi_m|\le\theta/125 ,
 \end{aligned}
\]
where \(L,U\) label the two endpoint scales. The error
follows by integrating the derivatives of the half logarithm
errors from angle zero along their respective circles.
The angular derivative bound on each circle is \(3/999\).
The radius bound gives
\[
 B_m^{L}|r-1|+B_m^{U}|r^{-1}-1|
    \le128 D(D+1)/m<1/100
\]
for sufficiently large \(N\), by (1).
If a physical root were common to \(N-1,N\), its two
normalized roots would have the same angle by \(G_m>0\).
We may take that angle in \((0,\pi]\) by conjugation and
positivity of the coefficients. By (5) the two angular
corrections differ by at most \(\theta/4\) for large \(N\).
Subtraction would give
\[
 2\pi(j_N-j_{N-1})=\theta-\delta_N+\delta_{N-1}
                     \in(0,2\pi)
\]
for large \(N\), a contradiction. The ordinary moment
polynomials have integer coefficients, so complex
coprimality gives a rational Bézout identity and hence
coprimality over every characteristic-zero field.

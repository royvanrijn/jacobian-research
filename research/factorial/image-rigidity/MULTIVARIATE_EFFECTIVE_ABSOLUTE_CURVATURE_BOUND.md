# An effective absolute curvature bound for consecutive factorial moments

Let \(Q_j(t)=\sum_{k=0}^j\binom jk\prod_i((j-k)a_i+kb_i)!\,t^k\)
for distinct nonnegative integer exponent columns. Remove the equal
rows and put
\[
 S=\sum_i\frac{(b_i-a_i)^2}{\max(a_i,b_i)},\qquad
 V=\sum_{a_ib_i>0}\frac{(b_i-a_i)^2}{\min(a_i,b_i)}
   +\sum_{a_ib_i=0}\max(a_i,b_i).
\]

**Theorem.** If \(S\ge64\), then \(\gcd(Q_{N-1},Q_N)=1\) for every
terminal \(N\ge2\), over every characteristic-zero field.

This makes the large-curvature step of the
[uniform eventual proof](MULTIVARIATE_ALL_DEGREE_REDUCTION.md)
effective. Any open terminal has \(S<64\).

## Reduction and quantitative scales

Write \(n=N-1\). Use the endpoint normalization from the
[large-curvature angular comparison](MULTIVARIATE_HIGH_CURVATURE_ANGULAR_SLOPE.md)
and write its coefficients as \(c_{m,k}\), for \(m=n,N\). Thus both
endpoint coefficients equal one, physical and normalized roots have
the same arguments, and
\[
  c_{m,k}\le\binom mk\exp[-S k(m-k)/(2m)].
  \tag{1}
\]
For \(N\le256\), this gives \(c_{m,k}\le64^{-\min(k,m-k)}\):
for \(k\le m/2\), use \(\binom mk\le N^k\) and
\(N e^{-S/4}\le256e^{-16}<1/64\), and then use symmetry of the bound.
Split the normalized polynomial into halves
\[
 T_m(z)=L_m(z)+z^m U_m(1/z),\quad
 L_m(z)=\sum_{k\le\lfloor m/2\rfloor}c_{m,k}z^k,\quad
 U_m(w)=\sum_{0\le k<m-\lfloor m/2\rfloor}c_{m,m-k}w^k.
\]
On \(|w|\le e^{1/m}\) either half \(H\) satisfies
\(|H-1|<1/15\) and \(|wH'/H|<8/105<1/8\), by summing the
geometric series with \(e^{1/m}/64<1/16\). Rouché on the two
reciprocal circles places both normalized root moduli in
\((e^{-1/m}, e^{1/m})\), since \(e(14/15)>16/15\).
At a root with argument \(0<\theta\le\pi\), the two half arguments
give \(m\theta=(2j+1)\pi+\delta_m\), with
\(|\delta_m|<\theta/4\). The arguments are taken from their real
values on the positive radius. Subtracting these two phase equations
at a common physical root would put an integer multiple of \(2\pi\)
strictly between \(\theta/2\) and \(3\theta/2\), which is impossible.

For \(N\ge256\) the cited angular comparison applies when
\(S\ge2\log N+5\). In the remaining region put
\[
 R=\log N>29.5,\qquad W=1+V,\qquad D=N e^{-S/2},\qquad
 K=\lceil128(D+R)\rceil .
\]
For a shared row \((l,l+d)\), \(l,d\ge1\), its contributions satisfy
\[
 \frac{d^2}{l}\le
 \frac{d^2}{l+d}\left(\frac{d^2}{l+d}+2\right).
\]
Indeed \(1+d/l\le 2+d^2/(l+d)\) is equivalent to
\(d^2-l^2\le d^2 l\), which holds for \(l\ge1\). The same comparison
holds upon reversing the row, and a private row contributes equally
to \(S,V\). Summing gives \(S\le V\le S(S+2)\) and \(W\le(S+1)^2\).
As \(K\le129(D+R)\) and \(n\ge N/2\), it follows that
\[
 \frac{WK}{n}\le258\!\left[(S+1)^2e^{-S/2}
                      +(2R+6)^2 R e^{-R}\right]<2^{-13},
 \qquad \frac{D+1}{n}<2^{-38}.                              \tag{2}
\]
For the first comparison, the two terms in brackets decrease on
\(S\ge64\) and \(R\ge29\), respectively. Their values at these lower
endpoints are less than \(2^{-32}\) and \(2^{-23}\), using
\(e^{-32}<2^{-45}\) and \(e^{-29}<2^{-40}\).
The second comparison follows from \(D/n\le2e^{-32}\) and
\(1/n\le2e^{-29}\). These comparisons use \(\log2<7/10\).

Each half has degree at least \(K\) in this region, by (2).
For \(K\le k\le m/2\), (1) and \(\binom mk\le(em/k)^k\) give
\[
                              c_{m,k}\le8^{-k}.             \tag{3}
\]
Indeed \(1-\log w-(S/2)(1-w)\) is convex for
\(K/m\le w\le1/2\). At the left endpoint it is at most
\(\log(e/128)+SK/(2m)<-\log8\), using (2).
Its value at the right endpoint is \(1+\log2-S/4<-\log8\).
These bounds also apply to the reciprocal half.

## Uniform half-logarithm bounds

Orient the profile so that the half being estimated starts at its
left endpoint. A private row \((0,v)\) at this endpoint gives
\(c_{m,k}\le e^{-(S-5)k/4}\) for \(k\le m/2\), by the
[private-endpoint coefficient estimate](MULTIVARIATE_UNBOUNDED_CURVATURE_EVENTUAL.md#coefficient-halves-and-endpoint-interpolation).
Therefore its half satisfies \(|H_m-1|\le4e^{-59/4}<1/2\)
on \(|w|\le2\), and its logarithm with value zero at \(w=0\)
has modulus at most \(8e^{-59/4}<2^{-18}\).

At a positive starting endpoint, let \(B_m=c_{m,1}\) at integer
order, with the Gamma interpolation \(B_m=m C_{m,1}/(C_{m,0}G_m)\)
for real \(m\in[n,N]\) as in the
[endpoint-logarithm proof](MULTIVARIATE_UNBOUNDED_CURVATURE_EVENTUAL.md#comparison-of-the-endpoint-logarithms).
The normalized coefficient interpolation and connected-hypergraph
bounds in that proof apply unchanged when (2) holds. More precisely,
put \(\beta=8eW/n\). All Gamma arguments needed in the derivative
bounds are at least \(m/2\), or at least \(ma_i/2\) in row \(i\),
and \(m/2>K\). The finite-difference estimates give
\(\max_{2\le r\le K} H_r\le4W/n<1/2\), using
\(2LK/n<1\) for the endpoint ratio bound \(L\le W\) of that proof.
Also \(2\beta K\le1/2\). These are the conditions used in its
bounds (5), (6), and (9).

Here the remaining derivative bound on \(B_m\) is explicit.
After the leading power of \(m\) is removed from the Gamma
ratio \(C_{m,1}/C_{m,0}\), finite products of linear factors
give derivative in logarithm of modulus at most \(2V/m^2\).
Each factor is at \(ma_i+j\) with \(|j|\le|b_i-a_i|\) and
\(ma_i+j\ge ma_i/2\). For the center \(G_m\), write
\[
 \log G_m=\Delta\log m+\sum_i(F(b_i)-F(a_i))+C(m)/m,\quad
 C(m)=\sum_i(E(mb_i)-E(ma_i)),
\]
where \(\Delta=\sum_i(b_i-a_i)\), \(F(x)=x\log x-x\),
\(F(0)=0\), and \(E(x)=\log\Gamma(x+1)-F(x)\).
For \(x>0\), use \(0\le E'(x)\le1/x\) and
\(|E''(x)|\le1/x^2\), as in the endpoint-logarithm proof.
A shared row contributes at most \(|\log(b_i/a_i)|\le d_i^2/\min(a_i,b_i)\)
to \(|C(m)|\), and at most \(2|\log(b_i/a_i)|/m\) to
\(|C'(m)|\), with \(d_i=b_i-a_i\).
For a private exponent \(v\), its contributions are at most
\(1+\log(mv)\le v(2+R)\) and \(1/m\) respectively.
Thus \(|C|\le V+S(2+R)\), \(|C'|\le(2V+S)/m\), and
\[
 \left|(\log B_m)'-1/m\right|
    \le [5V+S(3+R)]/m^2 <1/n .
\]
The last numerator is at most \(51R^2\), since \(S+1<3R\).
Dividing that bound by \(n\) gives at most \(102R^2 e^{-R}<1\)
for \(R>29\). By (1), \(B_n\le D e^{S/(2n)}<2D\),
which now gives \(B_m\le3D\) throughout the real order interval.

Let \(A_m(w)=\sum_{k=1}^K\ell_k(m)w^k\) be the truncated half
logarithm of the cited connected-hypergraph estimate. Its bounds
(8)–(11) give
\[
 |\ell_k(m)|\le B_m p_m^{k-1}/k,\qquad
 |\ell_k'(m)|\le 4B_m p_m^{k-1}/n,\qquad
 p_m=64e^2 WB_m/n<2^{-9}.
\]
For the last inequality, \(e^2<8\), \(B_m\le3D\), and \(128D\le K\)
give \(p_m\le12WK/n\). Consequently \(|A_m|\le8B_m\) on
\(|w|\le4\), \(|A_m|\le4B_m\) on \(|w|\le2\), and
\(\sup_{|w|\le2}|A_N-A_n|\le48D/n\).

At either integer order, \(H_m\) and \(e^{A_m}\) agree in
coefficients through order \(K\). Cauchy's bound for the tail
of \(e^{A_m}\) on \(|w|\le2\), divided by \(|e^{A_m}|\), is
\(2^{-K}e^{12B_m}\le e^{-K/8}\), using \(36D\le9K/32\)
and \(\log2>1/2\).
By (3) the tail of \(H_m\) after the same division is at most
\(4^{-K}e^{4B_m}/3\le e^{-K/8}\).
Thus
\[
       H_m(w)=e^{A_m(w)}(1+E_m(w)),\qquad
       |E_m(w)|\le2e^{-K/8}\le2N^{-16}\qquad(|w|\le2).
\]
Take the logarithm \(\mathcal H_m=\log H_m\) with
\(\mathcal H_m(0)=0\). Its difference from \(A_m\) is bounded
by \(4N^{-16}\). Combining the positive and private starting
endpoints separately, either orientation satisfies
\[
  \sup_{|w|\le2}|\mathcal H_m|\le M:=13(D+1),\qquad
  \sup_{|w|\le2}|\mathcal H_N-\mathcal H_n|<\delta:=2^{-16}.
  \tag{4}
\]
The endpoint support is fixed when comparing the orders, so each
orientation uses the same case at both orders.

## Consecutive angular separation

Let \(r_m=\exp(64(D+1)/m)\), which is less than \(3/2\) by (2).
On the circle \(|z|=r_m\), equation (4) implies that
\(z^m U_m(1/z)\) dominates \(L_m(z)\), since \(64(D+1)>2M\).
The first term is a polynomial whose zeros lie inside the
circle, by the zero-free half disk. Rouché and the reciprocal
argument therefore give \(|\log|z||\le64(D+1)/n\) at every root.

For normalized log radius \(\lambda\) and argument \(\theta\), put
\[
 J_m(\lambda,\theta)=\mathcal L_m(e^{\lambda+i\theta})
                     -\mathcal U_m(e^{-\lambda-i\theta})
\]
using the two half logarithms. If \(t\) is a common physical root,
both normalized roots have the same \(0<\theta\le\pi\) after
conjugation. At the corresponding radii \(\lambda_n,\lambda_N\) the
root equations are
\[
     m(\lambda_m+i\theta)
                =J_m(\lambda_m,\theta)+(2j_m+1)i\pi.
\]
The functions and derivatives in the following bounds are evaluated
where both arguments of the half logarithms have modulus at most
\(3/2\). Cauchy's estimate from (4) gives
\[
 \begin{aligned}
 |J_N-J_n|&\le2\delta,&
 |\partial_\theta(J_N-J_n)|&\le6\delta,\\
 |\partial_\lambda J_m|&\le6M,&
 |\partial_\lambda\partial_\theta J_m|&\le42M .
 \end{aligned}
\]
Indeed each half first derivative is bounded by \(2M\) and its
second derivative by \(8M\), with \(M\) replaced by \(\delta\)
for the order difference.

By the real part of the root equations,
\[
    (n-78(D+1))|\lambda_N-\lambda_n|
          \le 64(D+1)/n + 2\delta .
\]
Using (2) gives \(n|\lambda_N-\lambda_n|<2^{-13}\).
For real \(\lambda\), \(J_m(\lambda,0)\) is real.
Integrating the angular derivative at fixed radius, then the mixed
derivative when changing radius, gives
\[
 \left|\Im[J_N(\lambda_N,\theta)-J_n(\lambda_n,\theta)]\right|
  \le\theta\bigl(6\delta+546(D+1)|\lambda_N-\lambda_n|\bigr)
  <\theta/2 .
\]
Subtracting the imaginary parts of the root equations would again
put a multiple of \(2\pi\) strictly between \(\theta/2\) and
\(3\theta/2\). This proves complex coprimality. The ordinary
polynomial coefficients are rational, so their rational Bézout
identity gives coprimality over every characteristic-zero field.

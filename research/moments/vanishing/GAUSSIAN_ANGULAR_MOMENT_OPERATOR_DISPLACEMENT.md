# Square-root displacement for a Gaussian angular moment deformation

Let \(Z,W\) be a normalized circular Wick pair, and let the independent
ordinary Gaussians \(Y_i\) have positive integer exponents \(k_i\). For
\(s\ge0\), consider
\[
 P_s=W(1+\sqrt{s}\,Z+Z^2)\prod_iY_i^{k_i}.
\]
The even-time subsequence of its unit Wick period, divided by the positive
ordinary Gaussian moments, is
\[
 b_n(s)=[z^{2n}](1+\sqrt{s}\,z+z^2)^{2n}
                   \qquad(n=0,1,\ldots).                         \tag{1}
\]
This is a deformation of the Stieltjes angular slice in
\(\texttt{GAUSSIAN-WICK-BINOMIAL-MOMENT-CLASSIFICATION}\).
[Original Wick evaluation](GAUSSIAN_MONOMIAL_FIBRE_CONVERGENT_CORE.md#1-the-normalized-angular-germ)

For each \(s\), \(b_n(s)\) is the moment sequence of a compact nonnegative
law. The entire family has a simultaneous bounded self-adjoint
realization with a common unit vector and operator displacement at most
\(4\sqrt{s}+s\) from \(s=0\). Conversely, every simultaneous bounded
self-adjoint realization \(T(s)\) on a Hilbert space, with common unit
vector \(v\) and \(\langle v,T(s)^n v\rangle=b_n(s)\) for all \(n\),
satisfies
\[
             \|T(s)-T(0)\|\ \ge\ \frac{3}{32\mathrm e}\sqrt{s},
                         \qquad 0<s\le\tfrac14 .                 \tag{2}
\]
Here the family only needs to be defined at \(0\) and the indicated \(s\).
Thus the exponent \(1/2\) is optimal for the common-vector displacement.
In particular the family admits no norm-Lipschitz realization at zero.

## A matching realization and the lower bound

Write \(z=e^{i\varphi}\) on the unit circle. The coefficient in (1)
is the constant term of
\((z+z^{-1}+\sqrt{s})^{2n}\). With \(\varphi\) uniform on
\([0,2\pi]\), \(X=2\cos\varphi\) has the arcsine law on \([-2,2]\).
Therefore \(b_n(s)=\mathbb E[(X+\sqrt{s})^{2n}]\). On \(L^2\) of the
arcsine law take \(v=1\) and take \(T(s)\) to be multiplication by
\((X+\sqrt{s})^2\). Its displacement from \(T(0)\) is bounded by
\(4\sqrt{s}+s\).

For the universal lower bound put \(m_n=b_n(0)=\binom{2n}{n}\).
Expanding (1) gives the exact coefficient
\[
 \frac{[s^j]b_n(s)}{m_n}
     =\frac{\big(n(n-1)\cdots(n-j+1)\big)^2}{(2j)!}
        \quad (0\le j\le n).                                   \tag{3}
\]
Here the empty product is one. The central binomial bound is
\[
 \frac{4^n}{2\sqrt n}\le m_n\le\frac{4^n}{\sqrt{n+1}}
       \quad(n\ge1),
\]
with the upper bound also at \(n=0\). Both follow by induction using
\(4^{-(n+1)}m_{n+1}/(4^{-n}m_n)=(2n+1)/(2n+2)\).
This ratio lies between \(\sqrt{n/(n+1)}\) and
\(\sqrt{(n+1)/(n+2)}\), as follows by squaring.
The lower bound starts at \(n=1\), with equality.

Let \(0<s\le1/4\) and \(N=\lfloor s^{-1/2}\rfloor\ge2\).
Then \(1/4\le sN^2\le1\). Keeping the linear term of (3) gives
\[
 b_N(s)-m_N\ge sN^2m_N/2\ge 4^N/(16\sqrt N).
\]
For \(0\le j\le N\), the complete series bound from (3) is
\[
 b_{2j}(s)/m_{2j}
       \le\sum_{k\ge0}\frac{(2j\sqrt{s})^{2k}}{(2k)!}
       \le \mathrm e^2.
\]
Hence \(\|T(s)^jv\|\le\mathrm e\,4^j/(2j+1)^{1/4}\), with
the same inequality for \(T(0)\) without the factor \(\mathrm e\).
The noncommutative difference of powers and Cauchy–Schwarz give
\[
 b_N(s)-m_N\le \|T(s)-T(0)\|\,
          \mathrm e\,4^{N-1}
          \sum_{j=0}^{N-1}
           \big((2j+1)(2N-2j-1)\big)^{-1/4}.
\]
The sum is at most \((8/3)\sqrt N\). Indeed, on its first half
\(2N-2j-1\ge N\), and the remaining sum is at most
\(\sum_{k=1}^N k^{-1/4}\le(4/3)N^{3/4}\). Reflect to bound
the second half. Combining with the lower bound for
\(b_N(s)-m_N\) proves
\(\|T(s)-T(0)\|\ge3/(32\mathrm e N)\ge3\sqrt{s}/(32\mathrm e)\).

The transferred step is the use of the first deformation coefficient
at \(N\asymp s^{-1/2}\), controlled by telescoping against only
common-vector moment norms. The Gaussian deformation has a different
central-binomial moment growth from the Apéry family and also supplies
the matching simultaneous realization.
[Pinned Sokal coefficientwise Apéry proof, section on common-vector displacement](https://github.com/1337isnot1337/solved-conjectures/blob/dab1ac8d1ecea8dc34c38932fcd2a4b3f56ce0ef/research/sokal-apery-coefficientwise/OPERATOR-OBSTRUCTION.md#any-common-vector-self-adjoint-path-at-infinity-moves-by-at-least-order-sqrts)

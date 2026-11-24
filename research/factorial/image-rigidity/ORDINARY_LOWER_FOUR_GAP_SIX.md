# All ordinary orders at lower exponent four and gap six

**Theorem.** Let \(K\) be a field of characteristic zero. For every
nonzero binomial \(X^4(A+BX^6)\) in \(K[X]\) and every positive integer
\(n\), its ordinary factorial moments of orders \(n\) and \(n+1\) cannot
both vanish.

The proof combines an effective two-saddle estimate for every terminal
\(N=n+1\ge10{,}000\) with an exact finite-field check for
\(2\le N<10{,}000\). Both orders use the same dilation. The lower order
has the weight \(e^{-s}\), which separates the two saddle factors.
The [relative-gap weighted argument](RELATIVE_GAP_THREE_LOWER_WEIGHTED_SADDLE_TAIL.md)
uses the same distinction between the two orders. Here the physical
ratio \(d/a=3/2\) lies outside that argument's range. The contour and
numerical bounds for this support are certified below.

## 1. The actual pair

When \(B\ne0\), put \(w=A/(BN^6)\) and
\(\psi_w(s)=s^4(w+s^6)e^{-s}\).
With \(m=N,N-1\), the moment \(L((X^4(A+BX^6))^m)\) equals
\[
 B^m N^{10m+1} \int_0^\infty \psi_w(s)^m e^{-(N-m)s}\,ds.
\]
This follows from the ordinary Gamma integral on substituting \(x=Ns\).
The nonzero factors can be discarded. Split the two integrals at
\(c=63/10\), calling the lower piece \(K_m\) and the upper piece \(J_m\).
Both pieces are polynomials in \(w\) for each indicated \(m,N\).
Coefficient-zero phases are nonzero monomials, whose moments are
nonzero over \(K\).

The stationary equation and the negative curvature at a saddle \(r\)
are
\[
 w(r-4)=r^6(10-r),\qquad A(r)=\frac{(r-5)(r-8)}{r^2}.
\]
For \(|w|\le112000\) there is a unique simple upper saddle \(r_U\) in
\(|r-10|<113/100\). For \(|w|\ge88500\) there is a unique simple
lower saddle \(r_L\) in \(|r-4|<11/20\).

To see the first assertion, on the upper circle the squared modulus
of \(r^6(10-r)/(r-4)\) increases with the cosine of the polar angle.
Indeed, with \(R=113/100\), the nonconstant part of its logarithmic
derivative has positive sign because
\[
 10(36+R^2+12R\cos\theta)-(100+R^2+20R\cos\theta)>0.
\]
Its minimum is therefore \((10-R)^6R/(6-R)>112000\).
Rouché comparison with \(r^6(10-r)\) gives one zero counted with
multiplicity. On the lower circle the squared modulus also increases
with the cosine, since, with \(R=11/20\),
\[
 4(36+R^2-12R\cos\theta)-(16+R^2+8R\cos\theta)>0.
\]
The maximum \((4+R)^6(6-R)/R<88500\) gives the comparison with
\(w(r-4)\). In both cases the other factors stay nonzero in the
required disks. These strict comparisons hold on neighborhoods of the
closed parameter regions, and the selected saddles are holomorphic.

## 2. Complete contour bounds

Write \(\lambda_i=\psi_w(r_i)\) for \(i=U,L\). These are nonzero.
Deform the upper contour to height \(\Im r_U\) with a vertical
connector at \(c\). For the lower contour keep \([0,2]\) real and
move \([2,c]\) to height \(\Im r_L\), using both endpoint verticals.
The deformation is valid for the actual integer exponents since
\(\psi_w(s)^m e^{-(N-m)s}\) is entire. On the closing upper vertical
at large \(x\), exponential decay dominates the polynomial factors.
Each horizontal contour contains \(s=r_i+t,\ -2/5\le t\le2/5\).

The [directed interval checker](verify_ordinary_lower_four_gap_six_contours.cpp)
validates these bounds on the full closed parameter ranges.
The bounds on the third derivative hold along the local real
\(t\)-interval, with \(\Phi=\log\psi_w\) locally.
\[
\begin{array}{c|ccc}
 i & \Re A(r_i)\text{ at least} & |A(r_i)|\text{ at most}
   & |\Phi'''(r_i+t)|\text{ at most}\\ \hline
 U& 0.042 &0.14&0.075\\
 L& 0.074 &0.38&0.21
\end{array}
\]
On all other finite pieces up to \(x=20\) in the upper contour,
and on all other pieces of the lower contour, the modulus of
\(\psi_w(s)/\lambda_i\) is at most \(399/400\).
On the upper ray \(x=\Re s\ge20\), the exact derivative gives
\[
 \Re\Phi'(s)
 \le \frac{10}{x}-1+\frac6x
 \frac{112000/x^6}{1-112000/x^6}<-\frac{49}{100}.
\]
It uses the logarithm only on the far ray where the argument
is nonzero. Consequently the modulus bound on the complete tail
gains the factor \(e^{-49(x-20)/100}\).

Here are the certificate formulas. For either saddle \(r\), put
\(v=s/r\), and set
\[
 D=10-r+(r-4)v^6,\qquad
 \frac{\psi_w(s)}{\lambda_i}=\frac{v^4D}{6}e^{r-s},
 \qquad \eta=\frac{6(r-4)v^6}{D}.
\]
On the local interval the certificate proves \(D\ne0\), and
\[
 \Phi'''(s)=\frac{8+20\eta-15\eta^2+2\eta^3}{s^3}.
\]
Thus the local logarithm exists along the interval, continued
from the saddle. The box tests use these two rational expressions
with real exponential bounds for the squared modulus. For small
real \(t=s-r\), cancellation of the constant and linear terms
can be retained by evaluating \(v^4D/6-e^t\) as \(t^2\) times its
Taylor series. For \(|t|\le2\), the remainder of \(e^t/t^2\) after
the tenth Taylor degree is bounded by
\[
 \sum_{j=11}^{\infty}\frac{|t|^{j-2}}{j!}
 \le\frac65\frac{2^9}{11!}<\frac1{60000},
\]
since the successive term ratio is at most \(1/6\).
The disc and parameter tests discard a box only if it misses the saddle disc
or violates \( |r|^{12}|10-r|^2 \le112000^2|r-4|^2\) for the
upper packet, or the reverse inequality with \(88500\) for the
lower packet. Conjugation permits \(\Im r\ge0\) in these checks.
The checker subdivides boxes and requires every leaf to pass one
of the interval comparisons. There is no interpolation from
sampled contour points.

## 3. The weighted and unweighted packets

For \(m=N,N-1\), put \(\ell=N-m\in\{0,1\}\). With \(N\ge10000\),
\[
 (J_m,K_m)_i
 =\lambda_i^m e^{-\ell r_i}
   \sqrt{\frac{2\pi}{m A(r_i)}}\,(1+\varepsilon_{i,m}),
 \qquad |\varepsilon_{i,m}|<\frac9{20}.
\]
The notation on the left denotes \(J_m\) for \(i=U\) and \(K_m\)
for \(i=L\). The Gaussian branch is continued in \(\Re A>0\).

For the error, write \(F(t)=\Phi(r_i+t)-\Phi(r_i)\) and
\(F_2(t)=-A(r_i)t^2/2\). On the local interval
\(|F-F_2|\le M_i|t|^3/6\), where \(M_i\) is the third-derivative
bound in Section 2. The interpolation between \(F\) and \(F_2\)
has real part at most \(-b_it^2\) with \(b_U=0.016\) and
\(b_L=0.023\). This follows by using \(|t|\le2/5\) in the cubic
error estimate. Therefore
\[
 |e^{mF}-e^{mF_2}|\le\frac{mM_i}{6}|t|^3e^{-m b_it^2}.
\]
The amplitude satisfies \(|e^{-\ell t}|\le3/2\) and
\(|e^{-\ell t}-1|\le(3/2)|t|\). Integrating both errors over
the whole real line with the displayed Gaussian majorants and
dividing by the modulus of the leading Gaussian gives at most
\[
 \left(\frac{M_i}{4b_i^2}+\frac{3}{2b_i}\right)
     \sqrt{\frac{|A(r_i)|}{2\pi m}}<\frac{21}{50}.
\]
For the last inequality substitute \(m\ge9999,\ \pi>3\) and the
rational table bounds, then square both sides.

The Gaussian completion past \(|t|=2/5\) has integral at most
\((120/m)e^{-0.00336m}\), by \(\Re A/2\ge0.021\).
The total weighted nonlocal contribution after division by
\(|\lambda_i|^m|e^{-\ell r_i}|\) is at most
\((2400+3/m)e^{-m/400}\). The finite length is less than \(15\),
the amplitude correction is at most \(e^5<160\), and the
upper tail integral is at most \(3/m\) times its boundary
majorant. The full Gaussian leading modulus exceeds
\(3/\sqrt m\), so each of these two additional errors is
less than \(1/1000\) for \(m\ge9999\).
For instance \(e>8/3\) and \(m/400>24\) suffice at \(m=9999\)
for the nonlocal estimate, which decreases thereafter.
Adding the Gaussian, local and contour errors yields \(9/20\).

## 4. Excluding complex common roots

On \(|w|=88500\) one has
\(\log|\lambda_U/\lambda_L|>1/50\), and on \(|w|=112000\) one
has \(\log|\lambda_L/\lambda_U|>1/50\).
For this part the checker uses the rational stereographic samples
\[
 w=\rho\,\frac{\pm(1-t^2)+2it}{1+t^2},\qquad
 t=j/32,\quad 0\le j\le32,\quad \rho\in\{88500,112000\}.
\]
It rigorously encloses each saddle with a radius \(10^{-5}\)
Rouché disk inside the corresponding original saddle disk.
A numerical Newton iteration supplies a center only. The enclosure
inequality compares the constant and higher terms of the shifted
stationary polynomial with its linear term on that small circle.
The computed critical-factor ratio in the required direction
exceeds \(1035/1000\) at every sample.

At a stationary point the logarithmic envelope identity is
\(w\partial_w\log\lambda_i=(10-r_i)/6\). The angular derivative
of \(\log|\lambda_U/\lambda_L|\) has absolute value at most
\((1.13+0.55)/6=0.28\). Consecutive sample arguments differ by
at most \(2/32\), and the closest sample is within \(1/32\).
Conjugation covers the other half of each circle. Using
\(\log(1+x)\ge x/(1+x)\) gives the claimed circle margins.

For \(m=N\), the packet estimates imply \(|K_N/J_N|<1\) on
the smaller circle and \(|J_N/K_N|<1\) on the larger circle.
Indeed the non-exponential curvature factor is less than \(4\)
and the error ratio is less than \(3\). The same strict bounds
hold inside the smaller circle and outside the larger circle
by the maximum modulus principle. The required denominator is
nonzero throughout each region by Section 3. At infinity the
outer ratio is holomorphic since both pieces are polynomials
of degree \(N\) with nonzero leading coefficients. Thus any
common root lies in the annulus between the circles.

At an order \(N\) zero in this annulus, \(|J_N|=|K_N|\)
implies \(|\lambda_L/\lambda_U|^N<12\), hence
\(|\lambda_L/\lambda_U|<2\). If it were also an order \(N-1\)
zero, all four packets would be nonzero and the packet-ratio
formula would give
\[
 1=\left|\frac{J_{N-1}/J_N}{K_{N-1}/K_N}\right|
   <7\cdot2\,e^{\Re(r_L-r_U)}<1 .
\]
The first bound uses \(((1+9/20)/(1-9/20))^2<7\).
The last uses
\(\Re(r_U-r_L)\ge8.87-4.55=4.32\) and \(e^4>16\).
This proves the effective complex result at every \(N\ge10000\).

## 5. Exact prefix and field descent

The [finite checker](verify_ordinary_lower_four_gap_six_prefix.py)
forms the actual monic reciprocal polynomials
\[
 H_m(z)=\sum_{k=0}^m\binom{m}{k}\frac{(4m+6k)!}{(4m)!}z^{m-k}
\]
modulo the prime \(1000003\), verifying coprimality for every
\(2\le N<10000\). Trial division verifies primality and every
factorial argument is at most \(99990\), strictly less than the
prime. The full-degree polynomials therefore reduce to the
recorded pair over the finite field. A nonzero modular resultant
is a nonzero integer resultant. For the tail, absence of complex
common roots also gives coprimality of the two rational polynomials.
At every terminal a rational Bézout identity transfers the assertion
to any characteristic-zero coefficient field.

From the repository root, the two replays are

~~~sh
g++ -std=c++17 -O3 -fno-fast-math -ffp-contract=off \
  research/factorial/image-rigidity/verify_ordinary_lower_four_gap_six_contours.cpp \
  -o /tmp/sfc-gap-six-contours
/tmp/sfc-gap-six-contours
.venv/bin/python research/factorial/image-rigidity/verify_ordinary_lower_four_gap_six_prefix.py
~~~

The contour replay uses binary64 interval bounds rounded outward after
each operation. Its exponential enclosure uses rational powers of two
positive bases and requires no floating-point transcendental estimate.
The prefix replay uses exact FLINT finite-field polynomial arithmetic.
The [receipt](ordinary_lower_four_gap_six_verification.json) records
the internal replay. The proof concerns the single ordinary support
\(\{4,10\}\).


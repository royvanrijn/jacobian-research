# The one-mark test for rough Eisenstein sums

Let \(N\equiv2\bmod6\), \(z=N^{1/3}\), and \(Q(x,y)=x^2-xy+y^2\).
Use ordered pairs \(a+b=N\), \(a,b>1\), with \(ab\) squarefree
and every prime factor of \(ab\) congruent to \(1\bmod6\) and
strictly greater than \(z\). Each summand is prime or a product
of two distinct split primes. Write \(A_{ij}(N)\) for the number
of pairs with \(i\) and \(j\) prime factors, respectively.

The one-mark criterion in
[Rough Goldbach counts](ROUGH_GOLDBACH_SEMIGROUP_PARITY.md#what-the-quadratic-norm-proofs-control)
compares total mass \(R\) with twice the mass \(M\) where the
first summand is a semiprime. For distinct norm values it is
\(R-2M=A_{11}-A_{22}\). For the actual lattice weights, write
\(R^{\mathrm{lat}},M^{\mathrm{lat}}\).
A squarefree split value \(v\) has \(6\cdot2^{\Omega(v)}\)
Eisenstein coordinate witnesses, by unique factorization and
the six units. All these witnesses are primitive, since a
common rational prime divisor of the coordinates would give
a square divisor of \(v\). Thus
\[
 R^{\mathrm{lat}}-2M^{\mathrm{lat}}=144(A_{11}-4A_{22}).
 \tag{1}
\]
The mixed cells cancel in both versions. In particular
\(A_{11}>4A_{22}\) is sufficient for both.

**Theorem (almost-all comparison).** For every fixed \(A>0\), except
for \(O_A(X/(\log X)^A)\) targets \(X<N\le 2X\) with
\(N\equiv2\bmod6\), one has \(A_{11}\gg N/(\log N)^2\) and, uniformly,
\[
 \frac{A_{12}}{A_{11}}=\frac{\log 2}{2}
                  +O_A\!\left(\frac{\log\log N}{\log N}\right),\qquad
 \frac{A_{22}}{A_{11}}=\frac{(\log 2)^2}{4}
                  +O_A\!\left(\frac{\log\log N}{\log N}\right).
 \tag{2}
\]
In particular the marked fractions for values and lattice witnesses
tend to \(\log 2/(2+\log 2)\) and \(\log 2/(1+\log 2)\),
respectively, both below \(1/2\). The weighted margin
\((R^{\mathrm{lat}}-2M^{\mathrm{lat}})/(144 A_{11})\) tends to
\(1-(\log 2)^2\). This is a comparison of the three factor cells
on the same targets. The almost-all existence of two primes in
the prescribed progression is classical.
[Salmensuu, The Goldbach conjecture with summands in arithmetic progressions](https://doi.org/10.1093/qmath/haac008)

**Finite certificate.** For \(20\le N\le20\,000\,000\),
\(N\equiv2\bmod6\), the stronger inequality \(4A_{22}(N)\le
(44/57)A_{11}(N)\) holds with \(A_{11}(N)>0\).
The smaller targets \(2,8,14\) have empty rough shells. The scan starts at 8.
The certificate uses exact integer convolutions and 666 direct
cross-checks. The ratio \(44/57\) is a bound from the calculation,
using rigorous upper bounds for some semiprime counts.

## Circle comparison for the three factor cells

Put \(L=\log X\) and \(e(t)=\exp(2\pi i t)\). Freeze the semiprime
cutoff at either \(Z=X^{1/3}\) or \(Z=(2X)^{1/3}\). The corresponding
semiprime indicators sandwich the moving cutoff, including for a
mixed pair. Define
\[
 P(\alpha)=\sum_{\substack{p\le2X\\p\equiv1\bmod6}}e(\alpha p),
 \qquad
 S_Z(\alpha)=\sum_{\substack{Z<p<r,\,pr\le2X\\p,r\ \mathrm{prime}\\
                          p,r\equiv1\bmod6}} e(\alpha pr).
\]
Unique factorization in \(\mathbb Z\) makes the coefficients of \(S_Z\)
indicators. Deleting small prime summands costs \(O(X^{1/3})\),
uniformly for a fixed target. To pass from the frozen convolutions
to coprime pairs, delete a common prime factor \(q>Z\) dividing \(N\).
There are \(O(1)\) such factors, and each removes at most \(N/q\)
choices of the first summand. This costs \(O(X^{2/3})\), including
any diagonal. In both estimates the other summand is determined by
the first, so no divisor-function loss is needed.

### Minor arcs

Take \(Q_0=L^B\) for a fixed sufficiently large \(B\), to be chosen
for the error needed. The major arcs are the neighborhoods
\(|\alpha-a/q|\le Q_0/X\) modulo one with \(q\le Q_0\) and
\((a,q)=1\). They are disjoint for sufficiently large \(X\).
For every \(K>0\), one can choose \(B\) so that on the complement
\(\mathfrak m\),
\[
 |P(\alpha)|+|S_Z(\alpha)|\ll_K X L^{-K}
 \tag{3}
\]
at both frozen cutoffs.

For \(P\), expand the fixed residue condition as additive characters
modulo 6. For each shifted frequency \(\alpha+b/6\), Dirichlet
approximation with denominator at most \(\lceil X/Q_0\rceil\)
gives \(|\alpha+b/6-a/q|\le1/(q\lceil X/Q_0\rceil)\le q^{-2}\).
If \(q\le Q_0/6\), shifting back gives a major arc for \(\alpha\).
Thus \(Q_0/6<q\le\lceil X/Q_0\rceil\). Apply the
Vinogradov–Vaughan estimate
\((Y/\sqrt q+Y^{4/5}+\sqrt{Yq})(\log(2Y))^{5/2}\)
for the von Mangoldt exponential sum up to \(Y\).
Removing prime powers and partial summation give (3):
the part below \(X L^{-K-2}\) is bounded by its length, and for
larger \(Y\) the displayed terms have the requisite fixed-log saving
for large \(B\).
[Montgomery–Vaughan, Theorem 17.1 and the partial summation in Lemma 18.2](https://personal.science.psu.edu/rcv4/571s25/montgomery-vaughanII.pdf)

For \(S_Z\), omit products up to \(X L^{-K-2}\) at that cost.
On a remaining dyadic block \(p\asymp U,\ r\asymp V\), both lengths
are \(\gg X^{1/3}\) and \(UV\gg X L^{-K-2}\). Cauchy–Schwarz in
\(p\), followed by expansion in \(r_1,r_2\), bounds the squared
block by
\[
 \ll U\left(UV+V\sum_{1\le h\le2V}
                        \min\{U,\|\alpha h\|^{-1}\}\right).
\]
After dropping the outside prime indicator, the admissible \(p\)
for \(r_1,r_2\) lie in one interval, from \(p<r_i\), the dyadic
restrictions and the two product bounds. The geometric-series
bound applies to each nonzero \(h=r_1-r_2\).
With a Dirichlet approximation of \(\alpha\) itself and
\(Q_0<q\le \lceil X/Q_0\rceil\), the sum over \(h\) is
\(\ll(V/q+1)(U+q\log(2q))\). Indeed, in each interval of length
at most \(q/2\) the points \(\alpha h\) modulo one have spacing at
least \(1/(2q)\), so at most a bounded number can use the full
bound \(U\); summing reciprocals for the rest costs \(O(q\log(2q))\).
Dividing the squared block estimate by \((UV)^2\) gives
\[
 \ll V^{-1}+q^{-1}+\frac{\log(2q)}{U}
                      +\frac{q\log(2q)}{UV}.
\]
Taking \(B\) sufficiently large absorbs the \(O(L^2)\) dyadic blocks
and gives (3).

By Parseval, \(\int_0^1|F(\alpha)|^2\,d\alpha\ll X\)
for \(F=P,S_Z\). If \(F,G\) are any two of these sums, Bessel
and (3) give
\[
 \sum_{N\in\mathbb Z}
 \left|\int_{\mathfrak m} F(\alpha)G(\alpha)e(-N\alpha)\,d\alpha\right|^2
 \ll_K X^3 L^{-2K}.
 \tag{4}
\]
At error threshold \(XL^{-3}\), this excludes \(O(X L^{6-2K})\)
targets for each of the five needed convolutions
\(P^2,P S_Z,S_Z^2\) at the two cutoffs.

### Shared major arcs

For \(\alpha=a/q+\xi\) in a major arc, let \(\kappa_q(a)\) be the
mean of \(e(ar/q)\) over units \(r\bmod [6,q]\) with \(r=1\bmod6\).
Let \(V_P,V_{S,Z}\) be Fourier integrals with phase \(e(\xi t)\)
of densities \(v_P,v_{S,Z}\), supported on \([2,2X]\), where
\[
 v_P(t)=\frac1{2\log t},\qquad
 v_{S,Z}(t)=\frac14\int_Z^{\sqrt t}
                    \frac{du}{u\log u\log(t/u)}\quad(t>Z^2),
\]
and \(v_{S,Z}(t)=0\) otherwise. Both exponential sums have the
same arithmetic multiplier,
\[
 P(a/q+\xi)=\kappa_q(a)V_P(\xi)+O_{B,C}(XL^{-C}),\qquad
 S_Z(a/q+\xi)=\kappa_q(a)V_{S,Z}(\xi)+O_{B,C}(XL^{-C})
 \tag{5}
\]
for every \(C>0\). For the prime sum this is Siegel–Walfisz and
partial summation. The permitted reduced classes form a subgroup
of index two in the units modulo \([6,q]\); the product of two
independently uniform members has the same law. For the
semiprime sum apply Siegel–Walfisz in each factor, yielding the
integral of \(e(\xi uv)/(4\log u\log v)\) for
\(Z<u<v,\ uv\le 2X\). Substituting \(t=uv\) gives \(v_{S,Z}\).

For uniform error in the last step, first sum the larger prime
in a fixed reduced progression on \((p,2X/p]\).
Partial summation costs \(O(1+|\xi|X)\).
The error per smaller prime is
\(O((X/p)(1+|\xi|X)L^{-D})\) with arbitrary \(D\).
In the second prime sum, the function of the smaller factor
\(u\in[Z,\sqrt{2X}]\) given by integrating
\(e(\xi uv)/\log v\) on \(u<v\le2X/u\) is
\(O(X/u)\), and its derivative is
\(O((1+|\xi|X)X/u^2)\). Boundary differentiation uses
\(u^2\le2X\). Applying partial summation again gives error
\(O(X(1+|\xi|X)L^{1-D})\) in each choice of classes.
The number of classes and \(1+|\xi|X\) are bounded by fixed
powers of \(L\), so increasing \(D\) gives (5).
[Green, Proposition 4.2 and equation (4.4)](sources/green-roth-primes-2005.pdf)

For any of the three products \(FG\), integration over major arcs
now gives
\[
 \Sigma_{Q_0}(N)(v_F*v_G)(N)+O_B(XL^{-3}),\qquad
 \Sigma_{Q_0}(N)=\sum_{q\le Q_0}\ \sum_{(a,q)=1}
                         \kappa_q(a)^2e(-Na/q).
 \tag{6}
\]
The integrated error in (5) follows from taking \(C\) large in
terms of \(B\). Each density, extended by zero, has bounded
variation uniformly in \(X,Z\). For the semiprime density this
follows by writing it as the product of \(1/(4\log t)\) and a
bounded increasing logarithm on \(Z^2\le t\le 2X\).
Consequently each \(V(\xi)\ll|\xi|^{-1}\), and extending the
truncated Fourier integral to the real line costs
\(O(X/Q_0)\) per factor \(\Sigma_{Q_0}(N)\).
Fourier inversion of the integrable product gives the ordinary
real convolution in (6).

For completeness, the local factor in (6) is uniformly positive
on the chosen progression once \(B\) is large. Write \(q=d m\)
with \(d\) supported on 2 and 3 and \((m,6)=1\). At a prime
power \(2^a\) or \(3^a\) with \(a\ge2\), summing the nontrivial
additive character on the classes prescribed modulo 2 or 3
gives zero. For \(d\mid6\) the contribution from the fixed
residue is the Ramanujan sum \(c_d(N-2)=\varphi(d)\).
On \(m\) the multiplier is \(\mu(m)/\varphi(m)\) by CRT, so
\[
 \Sigma_{Q_0}(N)=\sum_{d\mid6}\varphi(d)
       \sum_{\substack{m\le Q_0/d\\(m,6)=1}}
             u_m(N),\qquad
 u_m(N)=\frac{\mu(m)^2c_m(N)}{\varphi(m)^2}.
\]
For odd primes \(p\ge5\), \(u_p=1/(p-1)\) if \(p\mid N\)
and \(u_p=-1/(p-1)^2\) otherwise.
The absolute mass \(H(N)=\prod_{p\ge5}(1+|u_p|)\) is
\(\ll\log\log X\) by Mertens and a split at \(L^2\).
It is at most a fixed multiple of the positive full product
\(\prod_{p\ge5}(1+u_p)\), since the ratios at \(p\nmid N\)
have convergent Euler product. Also
\[
 \sum_{(m,6)=1}|u_m(N)|\log m
 =H(N)\sum_{p\ge5}\frac{|u_p|\log p}{1+|u_p|}
 \ll H(N)\log L.
\]
For the potentially larger terms \(p\mid N\), sum
\(\log p/p\) up to \(L^2\) by partial summation, and bound the
tail by \(\log N/L^2\). Hence each truncated series in (6)
differs from the full product by \(O(H(N)/B)\), since
\(\log(Q_0/d)=B\log L-O(1)\).
This proves \(1\ll\Sigma_{Q_0}(N)\ll\log L\) and also
justifies the Fourier tail bound used in (6).

### Comparing the densities

The change of variables in the integral gives the exact formula
\[
 \frac{v_{S,Z}(t)}{v_P(t)}
    =\frac12\log\!\left(\frac{\log t-\log Z}{\log Z}\right)
 \quad(Z^2<t\le2X).
\]
On \(X L^{-4}\le t\le2X\), for either frozen cutoff, it equals
\((\log2)/2+O(\log L/L)\). The ratio is uniformly bounded
on the rest of the support. In a convolution at \(X<N\le2X\),
the portions where either argument is below \(X L^{-4}\)
contribute \(O(XL^{-5})\), since the other argument is of order
\(X\) and its density is \(O(1/L)\).
Moreover \((v_P*v_P)(N)\asymp X/L^2\), by the same upper bound
and by restricting the lower bound to \(N/3\le t\le2N/3\).
The convolution ratios in (6) are therefore \((\log2)/2\)
and \((\log2)^2/4\) with error \(O(\log L/L)\).
The same arithmetic factor multiplies all three convolutions.

Take \(2K\ge A+6\) in (4).
Outside the resulting exceptional set each minor-arc error is
\(O(XL^{-3})\). Sandwich the moving cutoff between the two
frozen ones and remove the small factors and shared primes
as above. The common prime-pair main term in (6) is
\(\gg X/L^2\), so division gives (2).
Using the two frozen cutoffs causes no target-by-target
loss of exceptional-set control.

## What remains for the pointwise test

A uniform sufficient estimate is now explicit. Write
\(I_P(N)=(v_P*v_P)(N)\) and \(Z=X^{1/3}\). A lower bound
\[
 \operatorname{Re}\int_{\mathfrak m}
       \big(P(\alpha)^2-4S_Z(\alpha)^2\big)e(-N\alpha)\,d\alpha
 \ \ge -\big(1-(\log 2)^2-\epsilon\big)
                       \Sigma_{Q_0}(N) I_P(N)
 \tag{7}
\]
for some fixed \(0<\epsilon<1-(\log 2)^2\) and every sufficiently
large admissible \(N\) would give \(A_{11}>4A_{22}\).
Indeed the true \(A_{22}\) is at most the count using the
frozen lower cutoff, and the prime endpoint deletion is
power-smaller than the main term. For almost all targets (4)
gives a much smaller error. At a prescribed target the
absolute integral bound available from Parseval is of order
\(X/\log X\), larger than the required scale \(X/(\log X)^2\)
before the singular factor. The proof therefore leaves (7)
open. The finite scan gives no estimate for this integral. The
[full factor-count comparison](EISENSTEIN_ROUGH_FACTOR_DISTRIBUTION.md)
extends the circle calculation to the (N^{1/5}) cutoff and
compares its divisor marks with and without lattice multiplicity.

## Exact scan and replay

Compile and run from the repository root with assertions active.

```sh
g++ -O3 -std=c++20 -Wall -Wextra -pedantic \
    research/arithmetic/verify_eisenstein_one_mark.cpp -o .scratch/verify_eisenstein_one_mark
.scratch/verify_eisenstein_one_mark 20000000
```

The index of a summand \(6i+1\) is \(i\), so convolution at
index \(n\) counts ordered pairs summing to \(6n+2\).
The least-prime-factor sieve supplies the prime indicator and
the direct comparison counts. The prime convolution is exact
after deleting primes at or below the target cutoff and the
diagonal. No pair can have both summands at or below the cutoff
in the scan range. For the semiprime convolution freeze the
cutoff at the lower end of each dyadic slab, generate the
products \(pr\) with \(Z<p<r\), and include shared factors
and the diagonal. This is an upper bound for the target's
semiprime cell. When the upper bound fails the strict
test, the program counts the real pairs directly.
The modular transforms use the prime \(998244353\) and
maximum length \(2^{23}\). The program tests primality and
the order of the generator. Padding removes cyclic wrap,
and every nonnegative convolution coefficient is at most
the number of input positions, which is smaller than the
modulus. No rounding is involved.

The run with GCC 13.3.0 covered 3,333,333 targets, including
\(8\) and \(14\). It used 666 direct cross-checks including
the two fallbacks, and those fallbacks had zero prime pairs
and zero semiprime pairs. The largest ratio \(4U/A_{11}\)
recorded using the exact or upper count \(U\) was \(88/114=44/57\),
at \(N=3458\). Every other scanned target had a strictly positive
margin. A separate enumeration of primitive coordinate vectors
at \(N=20,1964,3458\) checked the multiplicities in (1).
It found \((A_{11},A_{12},A_{22})=(114,10,0)\) at \(3458\),
so the recorded upper count 22 at that target overstates its
actual semiprime cell. The written almost-all proof and finite
scan are internal results with distinct scopes.

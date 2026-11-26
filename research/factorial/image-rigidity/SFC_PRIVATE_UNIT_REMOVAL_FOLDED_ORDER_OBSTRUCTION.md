# Folded likelihood-ratio order does not carry interlacing through private-pair removal

A fixed positive exchangeable core can satisfy both matched
folded likelihood-ratio inequalities at every degree and
still acquire a consecutive common root when its independent
private unit pair is removed. This happens along either of
two natural paths:

- Correlate the two standard complex Gaussian coordinates
  underlying the private pair. Some squared correlation
  \(0<\rho<1\) gives a common root of orders four and five.
- Let the common shape of the two independent private Gamma
  variables increase from one to infinity. Some finite
  \(\gamma>1\) gives a common root of the same orders.

In both cases the common root is on the unit circle. This
obstructs a proof from the folded law inequalities alone.
One example has a fixed exchangeable positive core with
all finite moments. A second uses a single mirrored Gamma
core with positive real, nonintegral exponents.

## A fixed core and its exact order comparisons

Let \(S\) have even density proportional to
\[
                        \exp(-10|s|-s^2/10000)
\]
and put \(U=e^{-S/2}\), \(V=e^{S/2}\). The law of
\((U,V)\) is exchangeable and \(UV=1\). In particular every
geometric mean tilt \((UV)^{j/2}\) gives the same law of
\(S=\log(V/U)\).

Write \(\mu^+\) for the law of \(|S|\), with density \(f\).
Dilation by \(0<t<1\) acts on the variable. The ratio of
the undilated density to the dilated one is
\[
     \frac{f(s)}{t^{-1}f(s/t)}
           = t\exp\left(10(1/t-1)s+(1/t^2-1)s^2/10000\right).
\]
It strictly increases on \(s>0\). Thus for the folded
laws \(\mu_j^+=\mu^+\) the full sandwich
\[
               (j/(j+1))\,\mu_j^+\ \preceq_{\rm LR}\
                           \mu_{j+1}^+\ \preceq_{\rm LR}\ \mu_j^+
                                                               \tag{1}
\]
holds in every degree. The Gaussian factor in the density
makes all exponential moments finite. Therefore the
[private-unit angular theorem for even midpoint laws](MULTIVARIATE_ONE_PRIVATE_MIRRORED_LINKED_PHASE_INTERLACING.md#strict-alternation)
applies at the independent shape-one starting pair.

Let \(M(s)=\mathbb E\cosh(sS)\). For either path, let \(w_{m,k}\)
be the normalized contribution of the private pair, including
the binomial coefficient. After endpoint normalization the
moment polynomial has coefficients
\[
               a_{m,k}=w_{m,k}\,\frac{M(k-m/2)}{M(m/2)}.
                                                               \tag{2}
\]
For the physical path,
\[
      w_{m,k}(\rho)=\sum_{j\ge0}\binom{k}{j}\binom{m-k}{j}\rho^j
\]
by the [mixed Gamma moment formula](SFC_PHYSICAL_MIRRORED_ROOT_CIRCLE.md#strict-coefficient-flattening-for-every-balanced-gap-pair).
For the independent shape path,
\[
        w_{m,k}(\gamma)
          =\binom{m}{k}\frac{(\gamma)_k(\gamma)_{m-k}}{(\gamma)_m}.
\]
Both paths start at \(w_{m,k}=1\) and tend to \(w_{m,k}=\binom mk\)
at the removal end.

## A strict degree-four circle along both paths

Put \(A=M(1)/M(2)\) and \(B=M(0)/M(2)\).
The order-four polynomial is \(1+a z+b z^2+a z^3+z^4\).
Under \(z=e^{i\theta}\), \(x=\cos(\theta/2)\), its centered
value \(e^{-2i\theta}P_4(e^{i\theta})\) equals
\[
               16x^4+(-16+4a)x^2+(2-2a+b).
\]
The discriminant in \(x^2\) is \(16(a^2+8-4b)\).
If the constant and discriminant are positive and \(0<a<4\),
then the two squared roots are positive and have sum
\(1-a/4<1\). Thus the four roots are distinct and lie
in \((-1,1)\).

For the physical path put \(\delta=1-\rho\).
Then \(a=(4-3\delta)A\), \(b=(6-6\delta+\delta^2)B\), and
\[
 \begin{split}
     2-2a+b&=(2-8A+6B)+6(A-B)\delta+B\delta^2,\\
    a^2+8-4b&=(16A^2+8-24B)+24(B-A^2)\delta+(9A^2-4B)\delta^2 .
 \end{split}                                                     \tag{3}
\]
We have \(A>B>0\) since \(S\) is nonconstant, and
\(B>A^2\) since \(M(2)=\mathbb E[2\cosh^2 S-1]>M(1)^2\).
The constant \(2-8A+6B\) equals
\(4\mathbb E(\cosh S-1)^2/M(2)>0\). The two remaining
signs in (3) follow from the exact estimates below.

For the shape path put \(p=1/(\gamma+3)\in[0,1/4]\).
The coefficients now are
\[
    a=(4-12p)A,\qquad b=(6-24p+12p^2/(1-p))B.
\]
At \(\delta=4p\), the physical path has the same \(a\)
and a coefficient of \(B\) at least this large, since
\(12/(1-p)\le16\). Hence \(a^2+8-4b\) is still positive.
Also
\[
      2-2a+b=(2-8A+6B)+24(A-B)p+12Bp^2/(1-p)>0.
\]
Throughout both paths \(0<a<4\) by \(A<1\). Consequently
the order-four roots stay distinct on the unit circle
through removal.

## Exact signs at the removal end

For the comparison without quadratic damping, let \(S_0\)
have even Laplace density \(5e^{-10|s|}\). Its moment generating
function is \(M_0(u)=100/(100-u^2)\) for \(|u|<10\). Multiplying
that density by \(\exp(-s^2/10000)\) gives the law of \(S\) after
normalization. For \(0\le u\le 5/2\) the resulting unnormalized
moment \(\widetilde M(u)\) satisfies
\[
 0\le M_0(u)-\widetilde M(u)
       \le \frac{\mathbb E[S_0^2\cosh(uS_0)]}{10000}
       =\frac{10}{10000}\left(\frac{1}{(10-u)^3}
                                   +\frac{1}{(10+u)^3}\right)
       <\frac1{200000}.
\]
For \(0\le u\le v\le5/2\), the normalized ratio equals
\(\widetilde M(u)/\widetilde M(v)\).
Its difference from \(M_0(u)/M_0(v)\) has absolute value
less than \(1/100000\): after clearing the denominator
\(\widetilde M(v)>1-1/200000\), its numerator is a difference
of two nonnegative errors less than \(1/200000\).

The undamped ratios in (3) are \(A_0=32/33\) and
\(B_0=24/25\). They give
\[
   16A_0^2+8-24B_0=\frac{136}{27225}>\frac4{1000},
             \qquad
   9A_0^2-4B_0=\frac{13984}{3025}>4 .
\]
Using \(0<A,A_0<1\), damping changes the first quantity
by less than \(56/100000\), and the second by less than
\(22/100000\). This proves the strict signs in (3).

At the removal end the degree-five coefficients are
\(1,a,b,b,a,1\), where
\[
                  a=5\,\frac{M(3/2)}{M(5/2)},\qquad
                  b=10\,\frac{M(1/2)}{M(5/2)} .
\]
The centered angular polynomial is
\[
                   2x\bigl(16x^4+(-20+4a)x^2+5-3a+b\bigr).
\]
The discriminant of the factor in \(x^2\) is
\(16(a^2+2a+5-4b)\). For the undamped law,
\(a_0=1875/391\), \(b_0=1250/133\), and
\[
        a_0^2+2a_0+5-4b_0
                      =-\frac{149760}{20333173}<-\frac7{1000}.
\]
The damping changes \(a\) by less than \(5/100000\) and
\(b\) by less than \(10/100000\). As \(a,a_0<5\), the displayed
quantity changes by less than \(1/1000\). The discriminant
therefore stays negative, and degree five has only its
structural \(x=0\) real root at the removal end.

## A common root at an interior parameter

The nonreal degree-five roots persist at correlations
sufficiently close to one and at sufficiently large finite
Gamma shapes, by the strict discriminant at removal.
At the initial pair the real angular polynomials
strictly interlace by (1) and the private-unit theorem.
The four simple order-four angular zeros vary continuously
in \((-1,1)\). If no order-five zero coincided with them
on the path to a parameter with negative discriminant,
the signs of the order-five polynomial at those four roots
would remain the initial alternating signs. Its leading
coefficient stays positive. Those signs force five distinct
real zeros in the five interlacing intervals. This contradicts
the negative discriminant. Thus a common real angular zero
occurs before removal, which gives a unit-circle common root
in the original variable.

## A single mirrored Gamma core with real exponents

The same obstruction occurs for a fixed pair of independent
Gamma variables \(X,Y\) with common unit scale and shape
\(10d\), using
\[
        d=10^{-8},\qquad U=X^{2d}Y^d,\qquad V=X^dY^{2d}.
\]
All moments of \(U,V\) are finite. At order \(m\), the
geometric mean tilt changes the Gamma shapes to
\(d(10+3m/2)\). The folded log ratio is the absolute value
of \(d\log(Y/X)\). By the [Gamma ratio density calculation](MULTIVARIATE_PRIVATE_UNIT_MIRRORED_ROOT_MONOTONICITY.md#an-even-log-ratio-law-and-a-two-term-root-equation),
its density on \(s\ge0\) is proportional to
\[
     e^{-a_m s}(1+e^{-s/d})^{-2a_m d},\qquad a_m=10+3m/2 .
\]
The interpolation formulas (2) use the appropriate
moment function \(M_m\) of this midpoint law in each
degree. The starting private unit pair interlaces by
the [Gamma extension of the private-unit theorem](MULTIVARIATE_ONE_PRIVATE_MIRRORED_LINKED_PHASE_INTERLACING.md).

For \(m=4,5\) the extra factor multiplying the Laplace
density of rate \(a_m\) is between \(1-2 a_m d\) and one.
Indeed \(2^{-2a_m d}\ge1-(2a_m d)\log 2\).
Thus for \(0\le u\le v\le m/2\), the ratio
\(M_m(u)/M_m(v)\) differs from the corresponding Laplace
moment ratio by at most \(2a_md/(1-2a_md)<10^{-6}\).
The latter ratio is at most one, and the normalization
constant cancels. At \(m=4\) the Laplace rate is \(16\),
with benchmark ratios \(A_0=84/85\) and \(B_0=63/64\).
They give
\[
        16A_0^2+8-24B_0=\frac{43}{57800}>\frac7{10000},
                \qquad 9A_0^2-4B_0=\frac{560889}{115600}>4 .
\]
At \(m=5\) the rate is \(35/2\), with benchmark coefficients
\(a_0=375/76\) and \(b_0=500/51\) at removal. They give
\[
                 a_0^2+2a_0+5-4b_0
                       =-\frac{245}{294576}<-\frac8{10000}.
\]
The first two signs change by less than \(56\cdot10^{-6}\)
and \(22\cdot10^{-6}\), respectively. The last changes by
less than \(100\cdot10^{-6}\), using the same coefficient
estimates as above. These signs preserve (3) at order four
and the negative discriminant at order five. The
continuity argument therefore gives an interior adjacent
common root on each private-pair removal path for this
single fixed Gamma core as well.

The small real gap in this example is essential to its
scope. The obstruction applies to the real-exponent Gamma
extension of the private-unit theorem and to the abstraction
using only folded order. The integer-gap circle-location
theorem imposes a separate restriction on ordinary
factorial profiles through both endpoints.

# Gamma quotients and coprimality on the binary factorial boundary

For every integer \(u\ge0\) and every \(m\ge1\), the factorial moments
of \(Y^u(aX+bY)\) in orders \(m,m+1\) have no common projective
coefficient zero over any characteristic-zero field. Interchanging
\(X,Y\) proves the two coordinate axes of
`GAP-FACTORIAL-BINARY-LINKED-BETA-COPRIMALITY`. The general two-positive,
unequal-exponent case remains open.

## 1. Gamma-quotient determinant signs

In the [external proof, Sections 2–7](https://github.com/1337isnot1337/solved-conjectures/blob/758066055ef49c2c66cb1087918ae62ad75e7ac1/results/kvz-gamma-quotient/PROOF.md),
put \(\delta=c-d>0\). Translate each finite increasing grid so that the
denominator is \(x_i+y_j\) and all \(x_i,y_j>0\). Such a split of the
translation \(d\) exists exactly when \(d+x_1+y_1>0\).
With \(G(a,t)=t^{a-1}e^{-t}/\Gamma(a)\), the change of variables
\(r=u+v,\ z=u/(u+v)\), of Jacobian \(r\), gives the decomposition
\[
 \frac{\Gamma(x+y+\delta)}{\Gamma(x+y)}
 =\int_0^\infty\int_0^\infty G(x,u)(u+v)^\delta G(y,v)\,du\,dv .
\]
A direct verification of the signs uses the following facts.

- For increasing \(a_i\), the initial Wronskians of \(e^{a_i t}\) are a
  positive exponential times the positive Vandermonde in the \(a_i\).
  Consequently the gamma-density determinant on ordered grids is
  strictly positive after extracting positive row and column factors.
- For fixed increasing \(v_j\), the order \(h\) initial Wronskian of
  \((u+v_j)^\delta\) has sign
  \(\operatorname{sgn}\prod_{j=0}^{h-1}(-\delta)_j\).
  Differentiation extracts \(\delta^{\underline j}\) from each row,
  \((u+v_j)^\delta\) from each column, and leaves the Vandermonde
  in the decreasing nodes \(1/(u+v_j)\).
- A system with all initial Wronskians nonzero has nonzero ordered
  evaluation determinants of the same signs. Dividing by the first
  function and differentiating reduces the system size, because the
  new Wronskians are the old Wronskians divided by the corresponding
  powers of the first function. Rolle's theorem gives the zero bound
  by induction. Coalescing ordered nodes determines the sign.
- Expand the three determinants in the two continuous Cauchy–Binet
  formulas. All terms are absolutely integrable by
  \((u+v)^\delta\le(1+u)^\delta(1+v)^\delta\). On the strictly ordered
  integration regions the determinant factors have fixed signs, and
  the product is strictly of that sign except in the finite-rank case.

Here \((z)_j=z(z+1)\cdots(z+j-1)\) with \((z)_0=1\). Thus for noninteger
\(\delta>0\) the order-\(h\) sign is precisely the sign of
\(\prod_{j<h}(-\delta)_j\). For an integer \(L=\delta\ge1\) the argument
works through order \(L+1\). Higher minors vanish because the quotient
is a polynomial of degree \(L\) in \(x+y\). In particular, order three
has sign \(\operatorname{sgn}(1-\delta)\), with zero at \(\delta=1\).
For \(c<d\) the Euler beta integral has positive measure and both
alternants \(t^{x_i}\), \(t^{y_i}\) have positive determinants by the
same exponential argument. This gives strict total positivity in that
regime whenever \(c+x_1+y_1>0\).

The external [literal Lean theorem](https://github.com/1337isnot1337/solved-conjectures/blob/758066055ef49c2c66cb1087918ae62ad75e7ac1/lean/Conjectures/KVZGamma.lean)
and [verification record](https://github.com/1337isnot1337/solved-conjectures/blob/758066055ef49c2c66cb1087918ae62ad75e7ac1/docs/verification/lean-kvz-20261121/receipt.json)
are pinned to the same revision. The written proof above rechecks the
strictness, rank and sign independently. The external verification
record reports a strict replay of the frozen closure. This audit uses
the written argument and does not claim a new Lean build.

### A product formula on unit-spaced columns

For increasing \(x_i\), put \(a_i=c+x_i+y_0,\ b_i=d+x_i+y_0>0\), also
assuming \(a_i>0\). Then
\[
 \det_{1\le i\le h,\ 0\le j<h}
       \frac{\Gamma(a_i+j)}{\Gamma(b_i+j)}
 =
 \prod_i\frac{\Gamma(a_i)}{\Gamma(b_i+h-1)}
 \prod_{i<j}(x_j-x_i)
 \prod_{j=0}^{h-1}(d-c)_j.                                         \tag{1}
\]
Indeed, scale each row to obtain evaluation of the polynomials
\((c+x+y_0)_j(d+x+y_0+j)_{h-1-j}\). Their determinant is a scalar
times the Vandermonde since their degrees are at most \(h-1\).
Evaluate at \(x_i=-c-y_0-(i-1)\) to make the polynomial matrix lower
triangular, with diagonal \((-1)^{i-1}(i-1)!(d-c)_{h-i}\). The
signs and factorials cancel those in the Vandermonde at these nodes.
This computation uses only polynomial evaluation after clearing the
Gamma factors. It gives the actual determinant value, including its
size and zeros, on unit-spaced columns.

## 2. An obstruction on the actual linked coefficient matrix

The linked orbit in the
[canonical-density audit](canonical-density/moment_pair_obstructions.md)
has row \(m\) given by the moments of
\(\operatorname{Beta}(mr+1,mu+1)\). After removing the
signed binomial factors, the raw moment matrix is
\[
 M_{m,k}=\frac{(mr+1)_k}{(m(r+u)+2)_k},\qquad
 \frac{M_{m,k+1}}{M_{m,k}}=\frac{mr+1+k}{m(r+u)+2+k}.
\]
For rows \(m_1<m_2\), subtract the consecutive-column ratios. The
numerator is \((m_2-m_1)(r-u(k+1))\) after putting over a common
denominator, so the column-\((k,k+1)\) minor has that sign.
Already \(r=2,u=1\) and \(m_1=3,m_2=4\) give a positive minor for
\(k=0\) and a negative one for \(k=2\), both in the actual range
\(k+1\le m_1\). Positive diagonal rescaling preserves both signs,
and uniform reversal of either index changes both together.

A composition of kernels having a fixed sign per minor order against
a positive measure has a fixed minor sign by Cauchy–Binet, when the
integrals are defined. This excludes such a composition for the
whole linked raw moment matrix in its native orders. It leaves
particular subquotients and the nonlinear Sylvester determinant
as possible receivers. The [shifted-resultant audit](SFC_SHIFTED_RESULTANT_POSITIVITY_AUDIT.md)
poses its positivity question for that different determinant.


The same polynomial determinant argument as in (1) yields a
stronger test on the linked matrix itself. For \(r>0,\ u\ge0\) and
increasing rows \(m_1,\ldots,m_h\), the consecutive-column minor is
\[
\begin{aligned}
 \det_{i,\ 0\le j<h} M_{m_i,k+j}
 ={}&
 \prod_i\frac{M_{m_i,k}}{((r+u)m_i+k+2)_{h-1}}\,
 \prod_{i<j}(m_j-m_i)\,r^{\binom h2}\\
 &{}\times\prod_{i=0}^{h-1}
        \left(1-\frac ur(k+i+1)\right)_{h-1-i}.
\end{aligned}                                                       \tag{2a}
\]
After clearing the indicated row factors the column polynomials
are \((r m+k+1)_j((r+u)m+k+2+j)_{h-1-j}\). Their determinant is a
constant times the Vandermonde. Evaluate at
\(m=-(k+1+i)/r\) to obtain the triangular determinant. Its
\(i\)-th diagonal entry is \((-1)^i i!\) times
\((1-\frac ur(k+i+1))_{h-1-i}\). Dividing by the evaluation
Vandermonde leaves exactly the claimed scalar.

In particular, if \(u/r=\eta>1\) is nonintegral, use
\(h=\lfloor\eta\rfloor+2\). At \(k=0\) exactly one of the
\(\binom h2\) linear factors of the final product in (2a) is positive.
It is the last factor for \(i=0\). Every other factor is strictly
negative, since \(h-1<\eta+1\) and \(i\ge1\) starts at
\(1-2\eta\) or lower. For sufficiently large \(k\) every factor is
negative, so the two minors have opposite signs. Take the same
rows in both minors with \(m_1\ge k+h-1\) to stay inside the actual
coefficient ranges. Reflection swaps \(r,u\). When the ratio of the
larger exponent to the smaller is nonintegral, the reversed
orientation with \(r>u>0\) already fails by the order-two calculation.
Thus neither orientation has a sign-regular raw matrix on these
profiles. For example, at \(r=2,u=3\) and rows \(3,4,5\), the
contiguous minors beginning at columns \(0,1\) are respectively
\(1/14632002\) and \(-7/139605768\).

## 3. All orders when one monomial exponent is zero

Normalize as in [the shifted-window identity](SFC_INITIAL_RESULTANT_AND_SHIFTED_WINDOWS.md).
For \(a=1\) and \(s=b/a\), the actual moment of \(Y^u(X+sY)\) at
order \(m\) is \(m!(mu)!\) times
\[
 P_{m,mu}(s),\qquad P_{k,q}(s)=\sum_{j=0}^k\binom{q+j}{j}s^j.       \tag{2}
\]
This follows by expanding \([z^m](1-z)^{-1}(1-sz)^{-mu-1}\).
The chart \(a=0\) has a nonzero pure monomial moment for \(b\ne0\).
It suffices to prove coprimality over \(\mathbb C\) since the moment
polynomials are rational. When \(u=0\), the difference of the two
polynomials is \(s^{m+1}\), and both have constant term one.

For real \(s\) the even member of the consecutive pair is strictly
positive. The factorial functional is integration against
\(e^{-X-Y}\,dX\,dY\) on the positive quadrant, where the even power
of \(Y^u(X+sY)\) is positive off its zero set. We only need to
exclude common nonreal roots when \(u\ge1\).

For \(q\ge0\) put \(A_{k,q}=(1-s)^{q+1}P_{k,q}\). One has
\[
 \begin{aligned}
 A_{k+1,q}-A_{k,q}
   &=\binom{k+q+1}{k+1}s^{k+1}(1-s)^{q+1},\\
 A_{k,q+1}-A_{k,q}
   &=-\binom{k+q+1}{k}s^{k+1}(1-s)^{q+1}.
 \end{aligned}                                                     \tag{3}
\]
The first adds the terminal coefficient and the second follows by
multiplying \(P_{k,q+1}\) by \(1-s\) and using Pascal's identity.
First advance \(k\) from \(m\) to \(m+1\), then advance \(q\) by \(u\).
With \(q=mu,\ v=1-s\), the result is
\[
 \frac{A_{m+1,q+u}-A_{m,q}}
      {\binom{m+q+1}{m+1}s^{m+1}(1-s)^{q+1}}
  =R(v)=1-(1-v)\sum_{\ell=1}^u a_\ell v^{\ell-1}
  =a_uv^u-\sum_{i=0}^{u-1}(a_{i+1}-a_i)v^i,
 \qquad a_\ell=\frac{(m+q+2)_\ell}{(q+1)_\ell}.                     \tag{4}
\]
All denominators are nonzero at a common nonreal root. The quotient
coefficients are governed by the same unit-shift relation as the
Gamma factors in (1), namely
\(a_{i+1}-a_i=(m+1)a_i/(q+i+1)\).

Set \(\rho=u/(u+1),\ L=a_u\rho^u,\ c_i=(a_{i+1}-a_i)\rho^i\).
Then \(R(\rho w)=Lw^u-\sum_{i=0}^{u-1}c_iw^i\), and
\[
 \frac{c_i}{c_{i-1}}
  =1-\frac{i+1}{(u+1)(q+i+1)}\quad(1\le i<u),\qquad
 (1-1/m)c_0<c_i\le c_0\le(1+1/m)/u.                              \tag{5}
\]
For the lower bound each factor in the ratio product is at least
\(1-1/(m(u+1))\), giving
\(c_i/c_0\ge1-i/(m(u+1))>1-1/m\).
Also \(a_i\rho^i\) is nondecreasing for \(0\le i\le u\), since the
one-step ratio is \(1+(u-i-1)/((u+1)(q+i+1))\) for \(i<u\).
It follows that \(L\ge1\) and
\(R(\rho)=1-u^{-1}\sum_{\ell=1}^u a_\ell\rho^\ell\le0\).
There is exactly one positive root \(w_0\) of \(R(\rho w)\).
The equation \(L=\sum_i c_i w^{-(u-i)}\) has a strictly decreasing
right side from infinity to zero, and gives \(w_0\ge1\).

For \(m\ge4\), take \(x=3/(mu)\). For \(1\le j\le u\),
\((1+x)^{-j}\le(1+jx)^{-1}\le1-jx/(1+ux)\). Sum this bound,
using (5) and \(\sum_{j=1}^u j\ge u^2/2\), to obtain
\[
 \sum_{i=0}^{u-1} c_i(1+x)^{-(u-i)}
 \le (1+1/m)\left(1-\frac{3}{2(m+3)}\right)\le1\le L.
\]
Thus \(1\le w_0\le1+3/(mu)\). Divide \(R(\rho w)\) by \(w-w_0\).
The quotient coefficients \(b_i\) satisfy
\(b_0=c_0/w_0,\ b_{i+1}=(b_i+c_{i+1})/w_0\), so \(b_i>0\) and
\(b_i\le(i+1)c_0\) by induction. For \(i\le u-2\),
\[
 (w_0-1)b_i\le \frac{3(i+1)}{mu}c_0
   <\frac3m c_0\le (1-1/m)c_0<c_{i+1}.
\]
Consequently \(b_{i+1}>b_i\). The quotient is constant for \(u=1\).
Otherwise all its roots satisfy \(|w|<1\) by the
Eneström–Kakeya bound. In the form used here, multiplying a
polynomial with nondecreasing positive coefficients by \(1-w\)
shows that the leading term dominates all the others for \(|w|>1\).
Scaling the argument by the largest ratio of consecutive
coefficients gives a strict disk when the coefficients increase strictly.

For \(P_{m,mu}\) the ratio of a coefficient to the next is
\((j+1)/(q+j+1)\le 1/(u+1)\) for \(0\le j<m\). The same bound
puts every lower root in \(|s|\le1/(u+1)\), hence
\(|v|=|1-s|\ge \rho\). By (4), a nonreal common root would be a
nonreal zero of \(R(v)\) with \(|v|\ge\rho\). All zeros except
the positive \(v=\rho w_0\) satisfy \(|v|<\rho\), a contradiction.

The exact Sylvester resultants for \(m=1,2,3\), respectively, are
\[
 \begin{split}
 &(u+1)^2,\\
 &\frac{(u+1)^3}{4}(17u^3+37u^2+22u+4),\\
 &\frac{(u+1)^4}{432}(44519u^8+235223u^7+518466u^6+622063u^5\\
 &\hspace{38mm}+444511u^4+193890u^3+50496u^2+7200u+432).
 \end{split}                                                       \tag{6}
\]
They are positive for \(u\ge0\). Each is directly checkable from (2).
The Sylvester products use coefficient indices summing to
\(m(m+1)\) and the coefficient of \(s^j\) has degree at most \(j\)
in \(u\). Thus the \(m(m+1)+1\) exact evaluations at successive
integers starting at zero give a full interpolation certificate
for each formula. The [shifted-resultant audit](SFC_SHIFTED_RESULTANT_POSITIVITY_AUDIT.md)
also certifies the full bivariate identities at these orders.
Equations (4)–(6) prove the boundary theorem for every starting order.

## 4. Other factorial obligations

The fixed-shape auxiliary production matrix in the
[cyclic Wronskian note](gap_three_hurwitz_structure.tex) already
has a totally nonnegative bidiagonal factorization on its stated
range. The actual \(H_3\) construction uses cyclic sections of a
rotated Wronskian, multiplying auxiliary minors by the signed
periodic weights \(1,1,0,-1,-1,0\). The positive-measure composition in
Section 1 gives no comparison between the positive and negative
terms. Thus `FACTORIAL-GAP3-HURWITZ-TN`,
`FACTORIAL-GAP3-TERMINAL-MINOR-STRICT` and
`FACTORIAL-LOWER-ONE-TERMINAL-DETERMINANT` stay open.
The explicit eventual cutoff for the ordinary univariate binomial
concerns a different support family.

## 5. Cumulant tails and zero-geometry tests

The [Nuttall proof, Sections 2–6](https://github.com/1337isnot1337/solved-conjectures/blob/758066055ef49c2c66cb1087918ae62ad75e7ac1/results/nuttall-cumulants/PROOF.md)
derives signs for the theta cumulants for all sufficiently large
integration orders, at each fixed minor order. Its transferable
step is the following determinant identity. For a smooth
rapidly decaying function \(\phi\) and its \(m\)-fold tail primitive
\(\Psi_m\), with \(m\ge h\),
\[
 \det[\Psi_m^{(i+j)}(x)]_{i,j=0}^{h-1}
 =
 \frac{\int_{0<t_1<\cdots<t_h}
  \prod_i t_i^{m-h}\,\prod_{i<j}(t_j-t_i)\,
  \det[\phi^{(j-1)}(x+t_i)]_{i,j=1}^{h}\,d\boldsymbol t}
 {\prod_{i=0}^{h-1}(m-i-1)!}.
\]
Integrate by parts \(i\) times in row \(i\). The boundary at zero
vanishes before each integration by parts. The resulting row signs
\((-1)^i\) cancel exactly the sign reversing \(h\) descending
powers of the integration variable in Cauchy–Binet. Both determinants
change sign under the same permutation, leaving a symmetric product.
The \(h!\) cancels on restriction to ordered variables.

In the original theta application, the first theta summand makes
\(\phi(t)\) asymptotic to a positive constant times
\(y^{9/4}e^{-y}\), \(y=\pi e^{4t}\), with a differentiable relative
error \(O(y^{-1})\). The Hankel-Wronskian identity
\(\tau_{h+1}\tau_{h-1}=\tau_h\tau_h''-(\tau_h')^2\) then gives
alternating tail signs: its leading logarithmic second derivative
is \(-16h y\). A determinant factorization evaluates a separated
box \(t_i=\tfrac18\log m+O_h(1)\) with the same sign and size bounded
below by a positive multiple of
\((\tfrac18\log m)^{h(m-h)}e^{-C\sqrt m}\) before the positive scalar.
The possibly opposite-sign region has some \(t_i\) in a fixed
interval and contributes at most
\(C U^m(2\log m)^{m(h-1)}\). The ratio tends to zero because a full
factor \((\log m)^m\) outweighs every fixed-base exponential in
\(m\). These statements use the theta tail and its derivative bounds,
and the compact residual range in \(x\) is uniform because a common
\(x\) multiplies every \(y_i\) by the same positive factor.
The method would need an actual tail-primitive representation with
these quantitative controls for a factorial determinant. The
ordinary-binomial explicit-cutoff obligation has no such identified
representation, and the identity alone gives no cutoff there.

The [Csordas–Vishnyakova proof](https://github.com/1337isnot1337/solved-conjectures/blob/758066055ef49c2c66cb1087918ae62ad75e7ac1/results/csordas-vishnyakova/PROOF.md)
has a directly verifiable two-parameter extension. Let
\(p(z)=z^4+A z^2+B\), where \(B>0\) and \(A\ge0\), and put
\(\varphi(z)=e^{-z^4}\). The roots of \(p\), equivalently their
negatives, have second and fourth power sums \(-2A\) and
\(2A^2-4B\), with odd power sums zero by evenness. Therefore the
product over their shifts is
\[
 \prod_{p(\alpha)=0}\varphi(z+\alpha t)
   =e^{-4z^4}\exp(12A z^2t^2)\exp((4B-2A^2)t^4).
\]
All its Taylor coefficients in \(t\) are nonnegative on the real
\(z\)-axis exactly when \(A^2\le2B\), the necessity following from
the \(t^4\) coefficient at \(z=0\). Yet
\(\varphi(it)\varphi(-it)=e^{-2t^4}\) has a negative coefficient.
A real-rooted approximating polynomial gives a product on the
imaginary axis with only nonnegative coefficients, since each real
root \(r\) contributes \(t^2+r^2\). Compact convergence preserves
each coefficient by Cauchy's formula. This verifies the failed
implication for a full range of test polynomials and excludes a
route from those coefficient inequalities alone to global
Laguerre–Pólya membership.

The [Cardon et al. counterexample](https://github.com/1337isnot1337/solved-conjectures/blob/758066055ef49c2c66cb1087918ae62ad75e7ac1/results/cardon-zero-sector/PROOF.md)
transfers its infinitesimal cubic test to the **actual** factorial
multipliers in the [Grace strip question](ITERATED_GRACE_POLYNOMIAL_ROOT_STRIP.md).
For positive weights \(\gamma_k\), let \(D_t(z^k)=\gamma_k^t z^k\).
On \(z^j(z-1)^3\) first remove the \(z^j\) factor and translate the
cubic to the origin. Its coefficients from degree three to zero
are \(1+O(t),O(t),O(t),t\Delta^3\log\gamma_j+O(t^2)\).
The cubic discriminant is
\[
 -27t^2(\Delta^3\log\gamma_j)^2+O(t^3).
\]
Hence any nonzero third difference forces a nonreal conjugate
pair at all sufficiently small positive \(t\). For the original
Cardon weights \(\log\gamma_k=-k^p\), \(p>2\), the third difference
at zero is strictly negative by integrating the strictly increasing
second derivative of \(x^p\) over translated unit squares. At
\(p=3\) and \(t=\log(16/15)/6\), scaling the variable as in the
external proof gives \(r^4z^3-3rz^2+3z-1\), \(r=15/16\), whose
discriminant is \(-72517275/4294967296\).

For the actual degree-\(n\) factorial factor take
\(\gamma_k=(an+1)_{dk}\), integers \(d\ge1\), \(a\ge0\), and
\(n\ge3\). For \(j+3\le n\) the third difference is strictly negative
by finite products alone:
\[
 \Delta^3\log\gamma_j
 =
 \sum_{\ell=1}^d
 \log\frac{(an+dj+\ell)(an+dj+\ell+2d)}
               {(an+dj+\ell+d)^2}<0 .
\]
Take \(j=n-3\) to make the witness degree exactly \(n\).
The input zeros are \(0,1\), in every positive-width sector about
the nonnegative ray. A nonreal output zero excludes any uniform
sector reduction, by taking initial width below its argument.
The cubic gives the same failure for the inverse fractional
multipliers, since changing the sign of \(\log\gamma\) does not
change the leading squared coefficient. Thus continuous
interpolation by positive fractional powers of the physical
multiplier cannot preserve real-rootedness or reduce sectors
throughout the interpolation. This does not address the direct
finite-order Grace strip or its actual common-root criterion
`GAP-FACTORIAL-POLYNOMIAL-GRACE-STRIP-NONCOLLISION`.



## Verification

The pinned KVZ diagnostic passed all 115 exact rational checks,
and all 20 recorded source hashes matched the pinned revision.
Here, exact rational arithmetic matched (1) on 32 test grids and
(2a) on 120 test grids, including positive, negative and zero
determinants. Direct polynomial expansion matched (4) on 49 pairs
of integer parameters. The determinant proofs above give the
universal identities.

The finite interpolation in (6) was checked at all 23 points by
rational elimination on Sylvester matrices. An independent short
replay using the repository's SymPy environment is

```python
from math import comb
from sympy import Poly, Symbol, resultant

s = Symbol("s")
quotients = (
    (1, [1]),
    (4, [4, 22, 37, 17]),
    (432, [432, 7200, 50496, 193890, 444511, 622063, 518466, 235223, 44519]),
)
for m, (denominator, coefficients) in enumerate(quotients, start=1):
    for u in range(m * (m + 1) + 1):
        lower = Poly(sum(comb(m * u + j, j) * s**j for j in range(m + 1)), s)
        upper = Poly(
            sum(comb((m + 1) * u + j, j) * s**j for j in range(m + 2)), s
        )
        expected = (u + 1)**(m + 1) * sum(
            coefficient * u**j for j, coefficient in enumerate(coefficients)
        )
        assert denominator * resultant(lower, upper, s) == expected
```

# Uniform angular stability of exact mirrored cores with a private unit pair

Fix an ordinary core with a distinguished private unit pair
\((1,0),(0,1)\) and any number of other mirrored rows.
Append arbitrary additional ordinary rows with signed gaps
\(d_i=b_i-a_i\ne0\) and midpoints \(m_i=(a_i+b_i)/2\).
The core remains exact. Put
\[
 W=\sum_i\frac{d_i^2}{m_i},\qquad
 B_\infty=-\sum_i\frac{d_i^3}{m_i^2},\qquad
 D=\sum_i\frac{d_i^4}{m_i^3},\qquad H=|B_\infty|+D,\qquad
 \rho=\frac1{128}.
\]
These sums are over the appended rows only. Let \(J\ge4\) be
any effective absolute terminal cutoff for ordinary angular
separation.

**Theorem.** The extended ordinary profile has disjoint consecutive
root arguments at every order if
\[
        2J\rho^{-5}H^2+33\rho^{-4}D<\frac1{2\cdot10^7 J^2}.
                                                               \tag{1}
\]
Another sufficient condition is \(H<\rho^4/(10^9J^2)\).
In particular there is one absolute computable uniform tolerance
on the relative gaps of the appended rows, independent of the
core and of the sum of appended curvatures. One may use
\[
        \max_i \frac{|d_i|}{\max(a_i,b_i)}
                        <\frac{\rho^4}{216\cdot10^9J^2}.
                                                               \tag{2}
\]
The ordinary result holds over every characteristic-zero field.
The cutoff \(J\) is unevaluated. The tolerance has polynomial
dependence on it. The core need not satisfy any smallness
condition, and the appended rows may have either orientation.
The shape sums in (1) can be read from any matched tanh pair.
For each Gauss factor use
\(u=(n/T'-(n-1)/T)^{-1}\). The exact matched grid gives
\(W=2\sum u\), \(B_\infty=-4\sum\epsilon u^2\) and
\(D=8\sum u^3\), with \(\epsilon\) the row orientation.
[Exact midpoint parameters](MULTIVARIATE_MIDPOINT_INFLATION_POLYNOMIAL_ANGLE_TRANSFER.md#direct-matched-tanh-criteria)

The proof upgrades the
[linked phase comparison](MULTIVARIATE_ONE_PRIVATE_MIRRORED_LINKED_PHASE_INTERLACING.md)
to a uniform angular gap. The exact core with a Gaussian tail
then supplies the reference for the
[midpoint inflation argument](MULTIVARIATE_MIDPOINT_INFLATION_POLYNOMIAL_ANGLE_TRANSFER.md).

## A uniform gap for the exact linked reference

In the private-unit phase formula set \(q=\sin\theta\),
\(c=\cos\theta\), \(\delta=1-c\), \(h=\sin(\theta/2)\), and
\(r_j=R_j(\theta)-1\). The half phase is
\[
     \phi_j(\theta)=\arctan\!\left(\frac{q}{r_j+\delta}\right),
                   \qquad \psi_j(\theta)=j\theta/2-\phi_j(\theta).
\]
This uses \(0<\theta<\pi\). The \(k\)-th positive angle is
labeled by \(\psi_j=(k-1/2)\pi\) and lies between
\((2k-1)\pi/j\) and \(2k\pi/(j+1)\). These facts also hold at
the geometric endpoint with the right angle-box endpoint
included. For two adjacent degrees \(n,N=n+1\), the lower
folded likelihood-ratio comparison gives
\[
    r_N\ge v=r\lambda(r),\qquad r=r_n,\qquad
    \lambda(r)=\exp[-A L],\quad
    A=2+h+\operatorname{arcosh}(1+r),\quad L=\log(1+1/n).
                                                               \tag{3}
\]
The differential argument in the linked phase comparison gives
(3) for all \(0<\theta<\pi\). Its differentiation and expectation
bounds use only \(\delta>0\) and the increasing function
\((\cosh(ts)-1)/(\cosh(ts)-c)\). This also covers \(r=0\)
by using \(r_N\ge0\).

At a positive lower root not on the negative real axis,
write
\[
                         \mathcal P=\theta/2+\phi_n-\phi_N.
\]
We show \(\mathcal P>1/(4\cdot10^6 n)\) for \(n\ge2\).
If \(\theta\ge2\pi/3\), then
\(\phi_N\le(\pi-\theta)/2\) gives \(\mathcal P\ge\pi/6\).
If \(\pi/2\le\theta\le2\pi/3\) and \(r\le4\), then
\(\phi_n\ge\arctan(1/11)\ge1/22\), using \(q\ge1/2\)
and \(\delta\le3/2\). Also \(\theta/2\ge\phi_N\).
If \(r\ge4\), the function \(r\lambda(r)\) is strictly increasing
since its logarithmic derivative is
\(1/r-L/\sqrt{r(r+2)}>0\). Thus \(v\ge4(2/3)^6>1/4\).
Here \(A(4)\le3+\log10<6\). Consequently
\(\phi_N\le\arctan(4/5)\) and
\(\mathcal P\ge\pi/4-\arctan(4/5)\ge1/10\).
This includes the entire positive lower half of order two.

Now \(\pi/n\le\theta\le\pi/2\) and \(n\ge3\).
If \(r\le q/8\), then \(A<6\) and
\[
             \phi_N-\phi_n\le(r-v)/q\le 3/(4n).
\]
The arctangent derivative has absolute value at most \(1/q\).
If \(r\ge8\), then \(v\ge8(n/N)^6\) by monotonicity.
For \(n\ge4\), this gives
\(\mathcal P\ge\theta[1/2-(5/4)^6/8]>\theta/50\).
For \(n=3\), use \(\delta\ge1/2\) and
\(\phi_N\le 1/[8(3/4)^6+1/2]<13/25\), while
\(\theta/2\ge\pi/6>523/1000\).

In the remaining case \(q/8<r<8\), define
\[
              b=\frac{2((r+1)c-1)}{r(r+2)},\qquad
              u=\max(0,rb).
\]
The angle addition identity gives
\(\phi_n+\theta/2=\arctan(q/(\delta+rb))\).
Since \(\delta+rb=(1+c)r/(r+2)>0\), it follows that
\[
           \mathcal P\ \ge\
             \frac{q(v-u)}{(r+\delta)^2+q^2}
             \qquad\text{whenever }v\ge u .                     \tag{4}
\]
This integrates the arctangent derivative from \(u\) to \(v\),
using \(v\le r\). For \(n=3\),
\(b<1/6\) by the degree-three estimate in the linked phase proof.
Thus \(v-u\ge r((3/4)^6-1/6)>r/100\).
With \(q\ge2/n\) and \((r+\delta)^2+q^2\le82\), (4) suffices.

For \(n\ge4\) we use the following margin for every \(b>0\):
\[
                        -n\log b - A >1/40.                    \tag{5}
\]
Indeed this difference increases with \(\theta\) on the positive
\(b\) interval, by the logarithmic derivative from the linked
phase proof. At \(\theta=\pi/n\) we have
\(-n\log b\ge F(r)=2\sqrt{46(r+1)/(5(r+4))}\).
For \(r\le5\) the positive Bernstein bound in that proof gives
\[
 F(r)-\bigl(12/5+\operatorname{arcosh}(1+r)\bigr)
 >\frac{7853/2800}{9\cdot12}>1/40 .
\]
Here the difference of squares is bounded below by the
Bernstein margin divided by \(r+4\le9\), and the sum of
the square roots is below twelve. For \(r\ge5\) the logarithmic
estimate from the same proof gives a margin greater than \(1/10\)
at five and increasing thereafter. At \(\theta=\pi/n\),
\(A\le12/5+\operatorname{arcosh}(1+r)\).
This proves (5).

If \(b>1/8\), then (5) gives \(v-rb>r/(320n)\).
If \(b\le1/8\), use \(\lambda\ge(4/5)^6>1/4\) to obtain
\(v-u\ge r/8\). When \(r\le72/n\), we also have
\(r\le36q\) and \(\delta\le q\). Equation (4) gives the
lower phase bound \(1/(3507200n)\).
When \(r>72/n\) and \(b>0\), the simpler inequality
\[
        -\log b-A L\ge\log(1+r/2)-6/n
                         \ge r/6-6/n >r/12
\]
gives \(v-u\ge(r/4)(1-e^{-r/12})\ge r^2/96\).
It also holds when \(b\le0\) by \(v\ge r/4\).
Because \(q<8r\), the denominator of (4) is at most
\(145r^2\), and (4) is at least \(q/13920\).
All cases prove the claimed uniform phase reserve.

To pass to actual angles, differentiation under the two
resolvent integrals shows
\(0\le(\log r_j)'\le q/\delta=\cot(\theta/2)\)
when \(r_j>0\). The lower sign uses the likelihood-ratio
comparison of the numerator and denominator weights for
\(R_j-1\), since \(1/(\cosh s-c)\) decreases in \(s\).
The upper bound replaces the denominator by \(\delta\).
It follows from the arctangent formula that
\[
              |\psi_j'(\theta)|\le 2j\qquad
                            (\pi/j\le\theta<\pi).
\]
For \(\theta\le\pi/2\), the bound on \(|\phi_j'|\) is
\(1+1/(2q)+\tfrac12\cot(\theta/2)\).
For \(\theta\ge\pi/2\), the three derivative terms are
bounded by \(1,1,1/4\), using \(\delta\ge1\) and
\((1+c)r_j/(r_j+\delta)^2\le1/(4\delta)\).
The geometric endpoint gives the same derivative bound
directly.

At the lower root, \(\psi_N\) exceeds the corresponding
upper label by \(\mathcal P\). Applying the derivative bound
between the already interlaced roots gives an angular distance
at least \(1/(8\cdot10^6 nN)\). The other adjacent gaps follow
directly from the angle boxes and odd structural root at \(-1\).
The case \(n=1\) also follows from its structural root and the
quadratic box. We use the spare margin in the bound
\[
        \text{minimum consecutive circular angle distance}
                               >\frac1{10^7 nN}.                \tag{6}
\]

## Retaining the core and transferring the tail

The reference at degree \(j\) keeps all core coefficients.
For the added rows replace their endpoint-normalized symbol by
a Gaussian with curvature
\[
                   C_j=\sum_i j d_i^2\psi'(1+jm_i).
\]
The Gaussian contributes variance \(C_j/j\) to the midpoint
log-ratio law. The curvatures increase in \(j\) by the positive
matched trigamma increment. Thus
\(C_N/N\ge(n/N)^2 C_n/n\), so the folded Gaussian at order
\(N\) dominates the dilated lower Gaussian in likelihood-ratio
order. The folded convolution kernel preserves the lower
comparison for the core and Gaussian laws successively, by the
[folded convolution proof](MULTIVARIATE_PRIVATE_UNIT_MIRRORED_ROOT_MONOTONICITY.md#adding-the-gaussian-symbol).
Consequently (6) holds for these exact matched reference clocks.

For the reference at each degree \(j\) the angle boxes give a
circular mesh at least \(\pi/j\). The corresponding roots of
\(1+z^j\) can be labeled with angular discrepancies at most
\(\pi/j\). Along linear interpolation in those angular coordinates
the mesh stays at least \(\pi/j\), including the wraparound gap.
At cyclic rank distance \(k\) the sine of the half-angle difference
has modulus at least \(k/j\) throughout. The derivative of its
cotangent bounds the sum of differences from the regular polygon
by \(4\pi j\), using the change of the half-angle at most
\(\pi/j\) and \(\sum_{k\ge1}k^{-2}\le2\).
The regular polygon has zero signed cotangent sum.
Thus the normalized signed force of the reference is at most
\(2\pi\).

The [complex root chart estimate](MULTIVARIATE_MIDPOINT_INFLATION_POLYNOMIAL_ANGLE_TRANSFER.md#a-root-chart-along-a-short-complex-deformation)
therefore applies with \(M=1\), radius \(\rho=1/128\), and
real shift rate \(K=3\). In that estimate the only input
depending on the original Gaussian reference beyond its mesh
is the normalized signed force. Its contribution to the
exponent is now at most \(2\pi\rho\). The constant part of
the exponent is still below \(1/8\), since it is bounded by
\(2\pi\rho+38\rho^2+\rho/4\). The coefficient of the real
shift is at most \(1+1+1/4+(3/2+32)\rho<3\).
All root charts and Cauchy derivative bounds follow with
these constants.

Inflate only the appended midpoints, leaving the core fixed,
and keep the matched Gaussian curvature constant using the
compensating quadratic in the
[midpoint inflation identity](MULTIVARIATE_MIDPOINT_INFLATION_POLYNOMIAL_ANGLE_TRANSFER.md#interpolation-by-inflating-all-midpoints).
That argument differentiates the coefficient multiplier and
does not require the unchanged core to have binomial coefficients.
With \(A_i=m_i/|d_i|\) and \(b_j=|B_\infty|+D/j\), its
one-degree complex log-displacement bound is
\[
  E_j=2\rho^{-3}b_j+S_j,\qquad
  S_j=4\rho^{-4}\sum_i\frac{|d_i|}{(A_i-3)^3}.
\]
It holds if all \(A_i>3\) and \(E_j<\rho/j\).
The [circle projection refinement](MULTIVARIATE_QUADRATIC_SKEW_TANH_ANGLE_TRANSFER.md)
also applies. Its projected root-chart bound uses the same
mesh and the signed force only through \(\rho\,2\pi\).
The power-series exponent in that bound remains below \(1/8\).
Thus the angular displacement is at most
\[
                     2j\rho^{-5}b_j^2+(1+\rho^2/2)S_j .        \tag{7}
\]

## All ordinary degrees

For \(N\ge J\) use the
[absolute ordinary angular theorem](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md).
For \(j\le J\), either smallness test in the theorem forces
each \(A_i>6\), because \(D\ge A_i^{-3}\).
Then \(S_j\le32\rho^{-4}D\).
Under (1) the two terms of \(E_j\) are each less than
\(\rho/(2J)\). For the signed term, use
\(2\rho^{-3}H<\rho^{-1/2}/(\sqrt{10^7}J^{3/2})\)
and \(J\ge4\). For the other term use
\(32\rho^{-4}D<1/(2\cdot10^7J^2)\).
The bootstrap applies and (7) is bounded by the left side of
(1). Two such displacements are strictly smaller than (6).
Under the linear test instead, \(E_j\le34\rho^{-4}H\),
so the same comparison follows directly.

Finally the
[absolute curvature theorem](MULTIVARIATE_ABSOLUTE_CURVATURE_27.md)
already gives angular separation when the total curvature
\(S\ge27\). Otherwise the appended midpoint curvature
\(W\le2S<54\). If \(\eta\) is the relative gap in (2), then
\(\max_i|d_i|/m_i\le2\eta<1\) and
\(H\le W(2\eta+4\eta^2)\le216\eta\).
Condition (2) therefore implies the linear test.

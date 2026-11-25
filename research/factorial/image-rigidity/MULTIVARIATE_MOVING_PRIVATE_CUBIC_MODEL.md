# A uniform moving cubic model for one private linear variable

Orient the exponent columns so the only private row is
\((0,1)\). Delete equal rows and let \(a_i,b_i>0\) index the
other rows. Put \(d_i=b_i-a_i\), \(\Delta=1+\sum_i d_i\),
\[
 S=1+\sum_i\frac{d_i^2}{\max(a_i,b_i)},\qquad
 \sigma=\sum_i d_i^2/a_i,\qquad
 \xi_j=\log(j^\Delta t)+\sum_i d_i\log a_i
\]
on the local branch near the positive boundary. Consider
profiles with bounded \(S\) and \(\sigma\to1\), as
\(n\to\infty\). Their rescaling constants are
\[
 c_2=\tfrac12(1+\sum_i d_i^3/a_i^2),\quad
 \ell=c_2^{-1/3},\quad
 T_j=\ell^2 j^{1/3}(1-\sigma),\quad
 H_j=\ell j^{2/3}\xi_j .
\]

**Theorem.** Uniformly on compact real \(T\)-intervals and
compact complex \(H\)-sets, the normalized moving one-private
period converges to
\[
                F(T,H)=\int_0^\infty e^{-y^3/3-Ty^2/2+Hy}\,dy .
\]
There is a local interpolation with the same limit for its
\(T,H\) derivatives and bounded derivative with respect to
\(\epsilon=j^{-1/3}\). At integer orders its values agree with
the normalized factorial moment up to exponentially small
errors. The normalization never vanishes.

## Exact one-dimensional interpolation

The factor \(k!\) from the private row cancels the corresponding
part of the binomial coefficient. Thus the moment divided by its
constant coefficient is the sum over \(k=0,\ldots,j\) of
\[
 \frac{j!}{(j-k)!}\prod_i
       \frac{\Gamma(ja_i+kd_i+1)}{\Gamma(ja_i+1)}\,t^k .
 \tag{1}
\]
For \(q=k/j\) its analytic interpolation is
\(f_j(q)=A_j(q)\exp(j\Phi(q,\xi_j))\), where
\[
 \Phi(q,\xi)=\xi q+
             \int_0^q\left(\log(1-z)+
                          \sum_i d_i\log(1+d_i z/a_i)\right)dz .
 \tag{2}
\]
This identity defines \(A_j\) by the exact Gamma expression on
a fixed disk about \(q=0\), also for complex \(j\) in a fixed
small sector. Stirling gives
\[
 A_j(q)=A_0(q)+O(j^{-1}),\qquad
 A_0(q)=\frac{\exp(\tfrac12\sum_i\log(1+d_iq/a_i))}
                    {\sqrt{1-q}},
 \qquad A_j(0)=1 .
 \tag{3}
\]
The error is uniform with \(q\)-derivatives on smaller disks,
and \(\partial_j(A_j-A_0)=O(j^{-2})\) by Cauchy. Indeed, bounded
\(S\) bounds the ratios \(a_i/b_i,b_i/a_i\) uniformly. It bounds
\(\sigma\), \(\sum_i1/a_i\) and \(\sum_i|d_i|/a_i\) by a constant.
Thus the shared Gamma arguments stay in fixed sectors and
the sum of the absolute Stirling errors is \(O(1/|j|)\).
The falling Gamma factor in (1) gives the term \(\log(1-z)\)
in (2) and the denominator in (3).

The coefficient \(c_2\) stays above a fixed positive bound.
In fact \(d_i/a_i\ge-R\) for a fixed \(R<1\) on negative rows,
so \(2c_2\ge1-R\sigma\). All higher phase coefficients are
uniformly bounded on a fixed disk. Taylor's formula gives
\[
 \Phi(q,\xi)=\xi q-(1-\sigma)q^2/2-c_2 q^3/3+q^4 E(q)
 \tag{4}
\]
with \(E\) uniformly analytic there.

## From the lattice to a linked analytic integral

At \(\xi=0\) write \(R(q)=(1-q)\sum_i d_i^2/(a_i+d_iq)\).
Then \(\Phi''(q,0)=(R(q)-1)/(1-q)\) and
\[
   R'(q)=-\sum_i d_i^2 b_i/(a_i+d_iq)^2\le-c\sigma
                   \qquad(0\le q\le1),
\]
where \(c>0\) follows from the uniform row ratios.
At any fixed small \(L>0\), (4) gives strictly negative
\(\Phi(L/2,0)\) and \(\Phi_q(L/2,0)\) for all sufficiently
near-critical profiles. The derivative of \(R\) then
shows \(\Phi(q,0)\le-\eta_L<0\) for \(L/2\le q\le1\).
This persists for sufficiently small complex \(\xi\).
The summands of (1) in that range are exponentially
small uniformly. For the falling factorial use the
Riemann sum of \(\log(1-q)\), which differs from \(j\)
times its integral by \(O(\log(j+1))\), including the
final mesh interval. For shared rows, Stirling on the
whole real segment has a bounded logarithmic amplitude
by the same ratio bound as (3).

At integer \(j\), let \(L_j\) be a grid point within \(1/j\)
of \(L\). Sum on \([0,L_j]\) by integrating
\(\pi j\cot(\pi j q)f_j(q)\) on the rectangle of height
\(\pm\nu\), taking half weights at the endpoint residues.
Choose \(\nu\) small and fixed. After division by \(j\)
the constant cotangent terms, \(-i\) above and \(i\)
below, give \(\int_0^{L_j} f_j(q)dq\) by closing
their contours on the real segment. The nonconstant
terms on the horizontal sides have modulus
\(O(e^{-2\pi j\nu})\). Their phase increase above the
real segment is \(O(\nu^2+|\xi|\nu)\), using the real
Taylor coefficients in (2). The real maximum tends
to zero as \(\sigma\to1\), \(\xi\to0\), so these
horizontal errors are exponentially small after
fixing \(\nu\). At \(L_j\) the same is true of the
vertical correction and of the half-endpoint term,
by the strict gap. On the left wall,
\(\cot(i\pi j y)+i=-2i/(e^{2\pi jy}-1)\) for \(y>0\);
the lower half gives the opposite sign. Adding back
the omitted half-summand at zero produces
\[
 \frac1j\sum_{k=0}^j f_j(k/j)
 =\int_0^L f_j(q)\,dq+B_j(\xi)+O(e^{-\eta j}),
 \quad
 B_j(\xi)=\frac1{2j}
       +i\int_0^\nu\frac{f_j(iy)-f_j(-iy)}{e^{2\pi j y}-1}\,dy .
 \tag{5}
\]
In the left side \(f_j\) at the final grid point means the
original summand (1). Moving \(L_j\) to \(L\) costs only
another strict-gap error. The cancellation of the two
values at \(y=0\) makes the left-wall integral convergent.

Formula (5) interpolates analytically in a fixed narrow
sector of complex \(j\), omitting the exponentially
small error. For \(|\arg j|\) small and \(0\le y\le\nu\),
differentiate \(f_j(q)\) on the segment joining \(-iy\)
to \(iy\) and use (2) to bound the numerator by
\[
 C y(1+|j|(|\xi|+y))\,e^{C|j|y(|\xi|+y)} .
\]
Also \(|e^{2\pi j y}-1|\ge e^{2\pi\Re j y}-1\).
Reducing the height and field neighborhoods gives
\(B_j=1/(2j)+O(|\xi|/|j|+|j|^{-2})\).
The same analyticity holds if the linear and quadratic
terms of (4) are replaced by external complex parameters
near their original values while keeping \(A_j\) and
the higher phase terms fixed.

## Bounded moving fields

Set \(\epsilon=j^{-1/3}\) and \(q=\ell\epsilon y\).
After dividing (5) by \(\ell\epsilon\), its integral
part has exponent
\(-y^3/3-Ty^2/2+Hy+\epsilon\ell^4 y^4 E(\ell\epsilon y)\)
and amplitude \(A_j(\ell\epsilon y)\). For small \(L\)
the cubic controls the remainder on the real interval:
the real exponent is at most
\(-c_1y^3+C(y^2+|y|)\) for bounded \(T,H\).
Differentiation in \(T,H,\epsilon\) introduces only
polynomial factors in \(y\) under this bound, by (3)
and (4). The upper endpoint \(y=L/(\ell\epsilon)\)
has an exponentially small contribution and derivative.
The normalized correction \(B_j/(\ell\epsilon)\) is
\(O(\epsilon^2)\) and has bounded \(\epsilon\)-derivative,
by Cauchy on a complex disk of radius proportional to
\(\epsilon\) using the analytic estimate in (5).
Its \(T,H\)-derivatives have the same small bound on
smaller compact sets. This proves the theorem.

**Corollary.** Suppose \(T_n\to0\), both scaled fields
\(H_n,H_{n+1}\) stay bounded, and
\(\Delta=o(n^{1/3})\). Every limiting scaled field of a common-root subsequence
is a multiple zero of \(F(0,\cdot)\).

Indeed, at the same physical \(t\) the exact field drift
is \(\xi_{n+1}-\xi_n=\Delta\log(1+1/n)\). Hence both scaled
fields would tend, along a subsequence, to the same zero
\(H_0\) of \(F(0,H)=\pi\mathrm{Hi}(H)\). Here \(H_0\ne0\) by positivity of the integral at zero.
Suppose it is simple.
[Gil, Segura and Temme, equation (2)](https://ir.cwi.nl/pub/10698/10698D.pdf)
By the implicit function theorem and Rouché's theorem,
the two integer-order roots must lie exponentially close
to the same local branch \(H(T,\epsilon)\) of the
interpolation. Its \(T,\epsilon\) derivatives are uniformly
bounded and \(H\to H_0\). For a fixed profile the physical
branch is \(\xi(\epsilon)=\epsilon^2 H(T(\epsilon),\epsilon)/\ell\),
with \(T(\epsilon)=\ell^2(1-\sigma)/\epsilon\). It satisfies
\[
 \xi'(\epsilon)=\frac{\epsilon}{\ell}
      (2H-T\partial_T H+\epsilon\partial_\epsilon H)
      =\frac{\epsilon}{\ell}(2H_0+o(1)).
\]
Thus its consecutive change is of order \(n^{-5/3}\)
with a nonzero leading coefficient. It cannot equal
the exact drift when \(\Delta=0\), or when
\(|\Delta|\ge1\) and the drift has magnitude at least
\(\log(1+1/n)\).

More generally, for a limit \(T_*\) and bounded fields,
a common-root sequence requires limiting zeros of
\(F(T_*,\cdot)\) at the same imaginary height. If the
limits coincide at a simple zero \(H_*\), the same
comparison applies to nonzero \(\Delta\) and to
\(\Delta=0\) whenever \(2H_*-T_*H_T\ne0\).
Thus the remaining bounded moving window requires
control of the same-height zero sets and the zero
derivatives of the entire family \(F(T,\cdot)\).

The model obeys \(F_{HH}+T F_H-HF=1\) by integration by
parts and \(F_T=-F_{HH}/2\) by differentiation under
the integral. At a simple zero this gives
\(H_T=1/(2F_H)-T/2\), so the adjacent derivative condition
can be expressed using its first \(H\)-derivative alone.
At a double zero, \(F_{HH}=1\) and \(F_T=-1/2\).

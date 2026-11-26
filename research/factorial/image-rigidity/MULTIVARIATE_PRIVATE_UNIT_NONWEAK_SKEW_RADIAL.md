# Absolute skew and integer drift for two private unit rows

Consider an ordinary two-monomial profile with one private unit
row on each side. For each unequal shared row give its signed gap
\(d_i=b_i-a_i\), smaller exponent \(\ell_i\) and midpoint
\(M_i=\ell_i+|d_i|/2\). Put
\[
        V=\sum_i d_i^2/\ell_i,\qquad
        S=\sum_i d_i^2/(\ell_i+|d_i|),\qquad
        \Delta=\sum_i d_i .
\]
Cancel equal multiplicities of reflected rows within each
shared pair, and let \(I\) be the remaining row multiset. Set
\[
 A=\min\left\{\sum_{i\in I}\frac{|d_i|^3}{\ell_i^2},\
   \left|\sum_{i\in I}\frac{d_i^3}{M_i^2}\right|
       +\frac13\sum_{i\in I}\frac{|d_i|^5}{M_i^4}
       +\frac14\sum_{i\in I}\frac{|d_i|^3}{\ell_i^3}\right\}.
 \tag{1}
\]
The cancellation is used only in this definition. The exact
factorial moment polynomials keep every row.

**Theorem.** If \(V>0\) and \(A/V\le10^{-15}\), the ordinary
consecutive moments are coprime in every degree. The signed
integer degree difference and the shared curvature are unrestricted.
At \(V=0\) consecutive geometric sums give the same result.
The constants are independent of the absolute eventual cutoff.

The second bound in (1) combines a signed second matched tanh
moment with unsigned remainders. Indeed, for each Gauss slot
of row \(i\), let \(v_i=|d_i|/(2M_i)\) be computed by
\(v_i^{-1}=n/T_i'-(n-1)/T_i\). The signed sum in (1) is four
times the sum of \(\operatorname{sgn}(d_i)v_i^2\) over the
slots of the unmatched rows.
[Matched tanh coordinates](MULTIVARIATE_MATCHED_TANH_INTEGER_MULTIPLICITY.md#recovering-the-gauss-grid)

## Curvature and signed skew

The [curvature fifteen-eighths theorem](MULTIVARIATE_BOTH_PRIVATE_UNIT_CURVATURE_15_OVER_8.md)
covers \(S\ge15/8\), so take \(0<S<15/8\). Every shared
\(r_i=|d_i|/\ell_i\) is at most two. If \(r_i>2\) then
\(|d_i|\ge3\), and that row contributes
\(|d_i|r_i/(1+r_i)>2\) to \(S\). Consequently \(V\le3S<6\).

After canceling the two private unit factorials, let \(f(x)\)
be the endpoint-normalized log coefficient at real index
\(0\le x\le m\), and \(g(x)=f(m-x)\).
At \(m\ge4\) the row trigamma integral bounds
give
\[
       \frac{V}{4m}\le f''(x),g''(x)\le\frac Vm,\qquad
       |f'|,|g'|\le V/2,\qquad f'(0),g'(0)\le-V/8 .
 \tag{2}
\]
Indeed the row curvature is bounded below by
\(d_i^2/[m(\ell_i+|d_i|)+1]\), at least
\(d_i^2/(4m\ell_i)\) since \(r_i\le2\).
The derivative bounds follow using the two zero endpoint
values. Putting \(t=\min(x,m-x)\), the lower curvature
in (2) gives \(f(x),g(x)\le-Vt/16\).

The direct trigamma skew bound gives
\(m\|f''-g''\|\le\sum_{i\in I}|d_i|^3/\ell_i^2\).
For the second bound in (1), write \(w=x/m-1/2\) and
\(h(t)=t/(e^t-1)\). The scaled curvature of one row is
\[
 m f_i''(x)=d_i^2\int_0^\infty e^{-(M_i+d_iw)u}h(u/m)\,du .
\]
In its difference with the reflection, replacing \(h\) by
one leaves the reciprocal difference
\(-2d_i^3w/[M_i^2(1-d_i^2w^2/M_i^2)]\).
Its leading term is \(-2w d_i^3/M_i^2\).
The absolute remainder is at most \(|d_i|^5/(3M_i^4)\),
using \(|w|\le1/2\) and \(|d_i|/M_i\le1\).
The error from replacing \(h\) by one is at most
\(|d_i|^3/(m\ell_i^3)\).
Indeed \(0\le1-h(u/m)\le u/(2m)\), and the absolute difference
of the two exponential factors is at most
\(2|d_iw|u e^{-\ell_i u}\). Integrate and use
\(2|w|\le1\).
Reflected rows cancel in the exact curvature difference,
so these bounds prove
\[
                       \|f''-g''\|\le A/m .
 \tag{3}
\]

Let \(\epsilon=A/V\).
For \(w_0(x)=(e^{f(x)}-e^{g(x)})/2\), the zero endpoints
of \(f-g\), (2), and (3) give
\[
 |w_0(x)|\le e^{-Vt/16}\frac{A t}{4},\qquad
 |w_0'(x)|\le e^{-Vt/16}\left(\frac A4+\frac{VA t}{8}\right).
\]
The sum of sampled absolute values is at most \(128\epsilon/V\).
The sampled total variation is at most \(72\epsilon\).
These follow respectively by summing on the two halves
with \(\sum_{k\ge1} k e^{-sk}\le s^{-2}\), and by
integrating the derivative on both halves.

## A normalized radial bound

Let \(C_m\) be the average of the normalized coefficient
array and its reciprocal. As in the
[cutoff-free weak radial proof](MULTIVARIATE_PRIVATE_UNIT_ABSOLUTE_WEAK_RADIAL_DRIFT.md#symmetrization-and-coefficient-bounds),
this symmetric polynomial is a geometric moment array.
It has simple unit-circle roots, angular spacing at least
\(\pi/m\), and no positive argument smaller than \(\pi/m\).
If \(A=0\), the actual polynomial is reciprocal and the
radial error is zero.

First suppose \(mV\ge1\). Use the same symmetric positive
moment law \(\mu\), \(Y=|\log T|\), and endpoint resolvent \(E\)
as in that proof. From (2),
\[
 \mathbb E_\mu Y\ge V/8,\qquad
 \mathbb E_\mu Y^2\le V^2/4+V/m\le5V^2/4.
\]
On the extension of the index interval to \([-1,1]\),
both log Mellin components have curvature at most \(2V/m\).
This follows from the smallest factorial argument \(m\ell_i-|d_i|\)
and \(|d_i|/\ell_i\le2\).
Taylor expansion at zero gives
\[
                  \mathbb E_\mu\cosh Y
                      \le e^{V/m}\cosh(V/2)\le1+3V^2.
\]
For \(V\le1\), the middle term is at most \(e^{9V^2/8}\);
convexity in \(V^2\) and \(e^{9/8}<4\) give the claim.
For \(1\le V\le6\), it is at most \(e^{3V/4}\).
The function \((e^{3V/4}-1)/V^2\) has its maximum
on this interval at an endpoint. The sign of its derivative
is the sign of \((3V/4-2)e^{3V/4}+2\), negative at one,
decreasing up to \(4/3\) and increasing afterwards.
At the endpoints use \(e^{3/4}<3\) and \(e^{9/2}<100\).

Let \(\theta\in(0,\pi]\) be a reference root angle and
\(D_\theta(Y)=\cosh Y-\cos\theta\).
The same exact reciprocal resolvent and Hölder bound used in the
weak radial proof give
\[
 \operatorname{Re}E\ge
  \frac{(\mathbb E_\mu Y)^3}
       {2\,\mathbb E_\mu Y(2V+Y)\,\mathbb E_\mu D_\theta(Y)}
 \ge \frac{V}{1536\,\mathbb E_\mu D_\theta(Y)}.
\]
Here \(r^3/(2r+5/4)\ge1/768\) at \(r=\mathbb E_\mu Y/V\ge1/8\).
Since \(\mathbb E_\mu D_\theta(Y)\le3V^2+\theta^2/2\),
this gives \(\operatorname{Re}E\ge1/(5376V)\) when
\(\theta\le V\). For \(V\le\theta\le\pi/2\), Jensen in
the imaginary resolvent gives
\(\operatorname{Im}E\ge1/(14\theta)\).
For \(\theta\ge\pi/2\), \(\operatorname{Re}E\ge1/2\).
Thus \(|E|\ge1/[5376(V+\theta)]\).
The reference phase derivative bound
\((\arg E)'\le\frac12\max(0,\cot\theta)\)
and \(\theta\ge\pi/m\) now give
\[
                       |C_m'(\zeta)|\ge\frac{m}{9000(V+\theta)}
                       \qquad (\zeta=e^{i\theta}).
\]
[Exact reciprocal resolvent and its phase derivative](MULTIVARIATE_PRIVATE_UNIT_ABSOLUTE_WEAK_RADIAL_DRIFT.md#a-reference-slope-at-every-root)

For \(|z-\zeta|\le1/m\), absolute summation of attenuated
coefficients gives \(|C_m(z)|\le132/V\), since \(V<6\).
Variation at most two for the symmetric coefficients gives
the alternative \(44/\theta\) bound after multiplying by
\(1-z\). Thus \(|C_m(z)|\le180/(V+\theta)\).
For the skew perturbation, the sum and variation above give
\(384\epsilon/V\) and \(864\epsilon/\theta\), respectively.
In each variation estimate use \(|z|^m<3\) and
\(|1-z|\ge\theta/4\). Combining gives
\[
              |P_m(z)-C_m(z)|\le\frac{1250\epsilon}{V+\theta}.
\]
For \(\eta=3\cdot10^7\epsilon/m\) and
\(\epsilon\le10^{-15}\), Cauchy's estimate on the half-size
disk bounds the Taylor remainder at radius \(\eta\) by
\(720m^2\eta^2/(V+\theta)\).
On the circle of radius \(\eta\), the linear term has
modulus at least \(3333\epsilon/(V+\theta)\), and the
remainder is at most \(648\epsilon/(V+\theta)\).
The bound is strict for Rouché's theorem. The reference
disks are disjoint, and conjugation handles the lower
half-plane. Thus every actual normalized root satisfies
\[
                          |\log|z||\le6\cdot10^7\epsilon/m .
 \tag{4}
\]

If \(mV\le1\), the
[bounded-\(mV\) radial argument](MULTIVARIATE_PRIVATE_UNIT_WEAK_RADIAL_INTEGER_DRIFT.md#a-quadratic-bound-on-the-radial-displacement)
applies with \(V^2\) as skew size replaced by \(A\).
For clarity, (3) and the slope bound in (2) give variation
at most \(9mA/16\) for \(e^f-e^g\), by the same zero-endpoint
calculation. The radius \(A\) is at most \(1/(4m)\).
The bound on the derivative of the symmetric first-difference
polynomial depends only on \(mV\), so remains valid.
The reference lower modulus is \(5mA/12\) after multiplying
by \(1-z\), against a perturbation below \(3mA/8\).
This gives \(|\log|z||\le2A\le2\epsilon/m\), which
also implies (4).

## Return to the physical variables

For \(N=n+1\), at a common physical root the
[linked endpoint-scale drift](MULTIVARIATE_ALL_DEGREE_SIGNED_IMBALANCE.md#an-explicit-linked-radial-condition)
and (4) give
\[
 |\Delta|\log(N/n)\le
          6\cdot10^7\epsilon(1/n+1/N)+V/(nN).
\]
For \(n\ge7\), multiply by \(N\) and use \(\log(N/n)>1/N\).
The resulting bound is
\[
                   |\Delta|<1.3\cdot10^8\epsilon+6/7<1.
\]
Thus \(\Delta=0\), and the
[balanced Euclidean theorem](MULTIVARIATE_BOTH_PRIVATE_EUCLIDEAN_DESCENT.md)
excludes a common root. The
[terminal-seven theorem and its predecessors](MULTIVARIATE_TERMINAL_SEVEN_ALL_PROFILES.md)
cover \(N\le7\). The resulting rational Bezout identities
transfer coprimality to every characteristic-zero field.

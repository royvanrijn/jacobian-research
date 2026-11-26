# Ordinary reflected rows unfold multiple interior Gaussian contacts

Fix adjacent orders \(n,N=n+1\), with \(n\ge2\), and the total signed
degree difference \(\Delta_0\). Use the
[fixed-terminal compactification](MULTIVARIATE_FIXED_TERMINAL_COMPACTIFICATION.md#exact-compactified-data).
The retained core \(R\) is a finite or countable multiset of ordinary
unequal exponent rows with summable curvature. Its endpoint-normalized
reference polynomials are
\[
 F_j(T,z)=\sum_{k=0}^j a_{j,k}e^{-jT\phi(k/j)}z^k,\quad
 \phi(s)=\tfrac12s(1-s),\qquad j=n,N .
\]
The center in the physical comparison is
\(\Gamma_R-\lambda/(2nN)\), where
\(\Gamma_R=\Delta_0\log(N/n)+\sum_R E(w)\) is the exact retained
center correction. Suppose there is a contact at \(T_*>0\) and
\(|\lambda_*|<T_*\). Its two nonreal root coordinates, before applying
the positive physical scale, have equal argument and positive
modulus ratio. Write their multiplicities as \(m_n,m_N\).

**Theorem.** Each of the following cases gives a sequence of finite
ordinary retained cores with simple transverse interior Gaussian
contacts at the same orders and at the same \(\Delta_0\). In all
cases the normalized reference pairs, signed bias and relative
physical scale converge to the original ones. Any strict total
curvature cap is preserved eventually.

1. If \(2\le m_N\le m_n\le3\), the new Gaussian strengths tend
   to \(T_*\).
2. If \(m_n=m_N=m\in\{4,5\}\), the new Gaussian strengths can tend to
   \(T_*-w\) for any fixed \(0<w<T_*-|\lambda_*|\). The added ordinary
   reflected rows have total variance tending to \(w\).
3. For each \(m\ge4\) there is a computable \(n_0(m)\) such that if
   \(m_n=m_N=m\) and \(n\ge n_0(m)\), the same conclusion holds for
   a choice of arbitrarily small \(w>0\). A finite one-variable root
   test below can also certify the reduction at smaller orders.

These are reductions of multiple contacts in the real compactified
family to transverse contacts of finite ordinary cores. The resulting
contacts retain positive Gaussian strength.

## Exact reflected row compensation

Take \(r\) copies of the reflected unit-gap pair with midpoint \(M\)
and put
\[
 u=M^{-1},\qquad w_M=2ru,\qquad A_M=2ru^2 .
\]
The midpoint is a half-integer, so all added endpoints are positive
integers for \(M\) large. At the new reference transfer \(w_M\) of
Gaussian strength to these rows. The effective strength before that
transfer is denoted by \(T\), so the actual Gaussian strength is
\(T-w_M\).

By the
[third-order chord calculation](MULTIVARIATE_THIRD_ORDER_BOUNDARY_AND_INTERIOR_REDUCTION.md#exact-row-expansions-and-the-projected-coefficients),
the sum of the compensated log multipliers of one reflected pair at
index \(k\) is
\[
 u^2\phi(k/j)+u^3\left(2j\Xi(k/j)-\frac{\phi(k/j)}{3j}\right)
       +O_j(u^4),
 \qquad
 \Xi(s)=\frac{(s-\tfrac12)^4-1/16}{12}.
 \tag{1}
\]
Reflection cancels the two cubic directions. In fact the whole
compensated log multiplier is analytic in \(u\). The rising-product
formula for the increasing row is the endpoint chord of
\[
 \sum_{\ell=1}^k \log(1+(\ell/j-1/2)u).
\]
At order \(u^h\) this produces a polynomial in \(k\) of degree at most
\(h+1\), since it sums a polynomial of degree \(h\) in \(\ell\).
Adding the reflected row and the compensation leaves an even
polynomial in \(k-j/2\). Thus for \(h\ge4\) the degree is at most
\(2\lfloor(h+1)/2\rfloor\).

Write \(\xi_j\) for a logarithm of the given root of \(F_j(T_*,z)\),
with the two imaginary parts chosen equal. The function
\[
 G_j(x)=e^{-jx/2}F_j(T_*,e^x)
\]
has a zero of order \(m_j\) at \(\xi_j\). For its nonzero leading
local coefficient \(g_j\) one has
\[
 A_M^{-m_j/2}G_j(\xi_j+\sqrt{A_M}\,X)\ \longrightarrow\ g_j X^{m_j}
 \tag{2}
\]
locally uniformly on complex \(X\)-disks, including derivatives.
In this logarithmic coordinate a multiplier polynomial in \(k-j/2\)
acts by the corresponding polynomial in differentiation. The clock
change \(T=T_*+A_M y\) contributes \(A_M y\,\partial_x^2/(2j)\)
up to a scalar. This is the trigonometric heat-flow identity after
changing to \(x=\log z\), as in
[Kabluchko, equation (1.4)](https://ahl.centre-mersenne.org/item/10.5802/ahl.227.pdf).
The two highest differential terms of the compensation (1) now give
\[
 a_j(y)\,\partial_X^2 + b_{j,M}\,\partial_X^4,\qquad
 a_j(y)=\frac{y-1/(2j)}{2j},\qquad b_{j,M}=\frac{1}{12j^3w_M}.
 \tag{3}
\]
The scalar contributions do not affect the zeros. The degree-four
term before dilation is \(r u^3\partial_x^4/(6j^3)\), and
\(r u^3/A_M^2=1/(2w_M)\).

If \(w_M\to w>0\), the new logarithmic polynomials divided by
\(A_M^{m_j/2}\) and the nonzero scalar factors converge to
\[
 g_j e^{a_j(y)\partial_X^2+b_j\partial_X^4} X^{m_j},
 \qquad b_j=\frac{1}{12j^3w}.
 \tag{4}
\]
Here and below the convergence is locally uniform in complex \(X,y\).
To justify (4), truncate the log multiplier at a fixed order \(K\)
in \(u\). The error on each of the finitely many frequencies is
\(O(u^{K+1}r)=O(u^K)\), which divided by \(A_M^{m_j/2}\) tends to zero
when \(K>m_j/2\). After the exact scalar is removed, the retained
exponent on each unscaled frequency is \(O(u)\). Thus its exponential
may also be truncated at a fixed power greater than \(m_j\) at that
stage. In the resulting finite operator expansion, the coefficient
of a differential monomial of degree \(d\) from order \(u^h\) in
one pair is \(O(r u^h A_M^{-d/2})\). For \(h\ge4\) this tends to zero
since \(h-1-\lfloor(h+1)/2\rfloor>0\). The degree-two contribution
at \(h=3\) tends to zero as well. The displayed operators in (3) are
the only nonvanishing nonconstant terms. Apply the resulting finite
operator to (2), using holomorphic derivative convergence. On
\(X^{m_j}\) the exponential in (4) terminates.

There is a second limit for \(r\to\infty\) and \(w_M\to0\) when
\(m_j\le3\), for example \(r=\lfloor\sqrt M\rfloor\). Directly on
the finitely many frequencies, the error in retaining just the
quadratic contribution at \(u^2\) in (1) is \(O(r u^3)\). This is
\(o(A_M^{m_j/2})\), since for \(m_j=2\) the ratio is \(O(u)\) and
for \(m_j=3\) it is \(O(r^{-1/2})\). The clock is exact and the
same expansion at (2) now gives the limit
\[
                     g_j e^{a_j(y)\partial_X^2}X^{m_j}.
 \tag{5}
\]
For a countable retained core, first truncate it with coefficient
and center errors \(o(A_M^{\max(m_n,m_N)/2})\). This is possible
by summability and does not change either limit. The compact
holomorphic convergence in (4) and (5) gives analytic continuation
of simple limit root branches on disks in the \(y\)-plane, with
convergence of their first derivatives. This follows by applying
Rouché on an isolating root circle and using a local Cauchy estimate.

## Vanishing variance and multiplicities two or three

For \(y>1/(2n)\), the two limiting polynomials (5) have simple
roots. The largest positive imaginary parts are
\[
 \begin{cases}
   \sqrt{(y-1/(2j))/j}, &m_j=2,\\
   \sqrt{3(y-1/(2j))/j}, &m_j=3 .
 \end{cases}
\]
The upper order has already split when \(y\) first passes
\(1/(2n)\). Its displayed positive imaginary part is therefore
larger than the lower one near that threshold. The difference
of the squares with the lower order first is affine in \(y\).
Its coefficient is strictly positive if \(m_n\ge m_N\) in the
stated range. Hence there is one transverse crossing in the
interior of this \(y\)-interval. The two simple ordinary
perturbation branches give a nearby zero of the exact angular
difference.

## Fixed variance and an effective root comparison

Use \(r=\lfloor w M/2\rfloor\). For \(X=i b_j^{1/4}Y\) the limit
(4), after a nonzero scalar division, becomes
\[
 Q_m(Y,\ell_j(y))
       =e^{-\ell_j(y)\partial_Y^2+\partial_Y^4}Y^m,\qquad
 \ell_j(y)=\sqrt{3jw}\left(y-\frac1{2j}\right)
 \tag{6}
\]
when the multiplicities equal \(m\). All exponentials in (6)
are finite on \(Y^m\).

For \(m\ge4\), there is a first value \(\ell_0>0\) at which all
roots of \(Q_m\) are real. To see this, for large positive
\(\ell\), dilate \(Y=\sqrt{2\ell}Z\) and use convergence to
the ordinary probabilists' Hermite polynomial of degree \(m\).
Its real roots are simple, so nearby real polynomials retain
that property. The set of parameters with all roots real
is closed and is stable under increasing \(\ell\) by the
backward heat flow. One can verify this stability for simple
roots directly. They satisfy
\(Y_i'=2\sum_{h\ne i}(Y_i-Y_h)^{-1}\), which is the gradient flow
of twice the log Vandermonde product on the ordered roots.
The product is nondecreasing and all roots stay bounded on
compact parameter intervals, so collisions among the simple
roots cannot occur. Approximation by simple real roots proves
stability also for repeated roots. For positive time repeated
roots split immediately. Indeed the local proof at a root
of order \(a\) is the same as (2) and (5), with the operator
\(e^{-\partial_Z^2}Z^a\), which has simple real Hermite roots.
Finally \(Q_m(Y,0)\) has a nonzero coefficient below degree \(m\)
and all its terms have nonnegative coefficients and exponents
congruent to \(m\bmod4\). It has no nonzero real root and
does not have all roots real. Stability rules out all negative
\(\ell\) as well. This proves \(\ell_0>0\). At the first value
there is a multiple root, since simple real roots would
remain real for a small negative parameter change.

Choose a multiple root \(R_0\ge0\) at \(\ell_0\), using parity
of (6), and take the largest branch in that root's cluster.
For \(\ell>\ell_0\), let \(R(\ell)\) denote this positive root
with its fixed rank in the ordered list. It is a real-analytic
algebraic function, all roots in that range being simple.
The local Hermite expansion gives
\[
                R(\ell_0+s)=R_0+c\sqrt s+O_m(s),\qquad c>0.
 \tag{7}
\]
At infinity this branch has a nonzero positive Hermite rate
\(R(\ell)\sim C_m\sqrt{\ell}\), \(C_m>0\), since its rank is
in the positive half of the roots.

Set \(D=N/n>1\) and
\[
                 \delta_w=\sqrt{3Nw}\left(\frac1{2n}-\frac1{2N}\right).
\]
At the lower real-root threshold \(x=\ell_n(y)=\ell_0\),
the upper parameter in (6) is \(\sqrt D\,\ell_0+\delta_w\).
The following strict root test is sufficient for the reduction
with some \(w\) arbitrarily close to the one tested:
\[
              D^{-3/4} R(\sqrt D\,\ell_0+\delta_w)>R_0 .
 \tag{8}
\]
Indeed, after factoring out \(b_n^{1/4}\) the two labeled
imaginary parts have difference
\[
                f_\delta(x)=R(x)-D^{-3/4}R(\sqrt D\,x+\delta),
                         \qquad x>\ell_0 .
 \tag{9}
\]
It is negative just above \(\ell_0\) by (8) and positive for
large \(x\) by the common asymptotic rate. The root ranks
used in (9) are simple in both limiting polynomials.

Under the strict test one may vary \(w\) a little, preserving
the two endpoint signs on a fixed compact interval \(x>\ell_0\),
to obtain a transverse crossing in (9). On that compact interval
and a compact small \(\delta\)-interval, finitely many solutions
of \(f_\delta(x)=0\) have \(\partial_\delta f_\delta(x)=0\).
In fact the latter condition requires \(R'(\sqrt D\,x+\delta)=0\),
whose possible arguments form a finite set by analyticity and
nonconstancy of \(R\). For each argument, \(R(x)\) is prescribed
one value and yields finitely many \(x\) in the interval.
At other solutions the zero set is locally the analytic graph
\(\delta=g(x)\). This graph cannot be constant on an interval,
since (9) would vanish identically in \(x\), contradicting the
different rates at infinity. Its stationary points are isolated.
Take a countable cover by such graph neighborhoods. All exceptional
\(\delta\) values therefore form a countable set. Choose one
outside it. The opposite signs give a transverse zero, and
that strict crossing persists for a neighborhood of \(w\).
A rational \(w\) can thus also be used when the strict test
is open to such variation.

For arbitrary fixed \(m\ge4\), test (8) at \(\delta_w=0\).
Equation (7) gives an increase of order \(n^{-1/2}\), while
the dilation factor \(D^{-3/4}\) causes only a loss of order
\(n^{-1}\). More explicitly, for \(n\ge1\),
\(s=\ell_0(\sqrt D-1)\ge\ell_0/(3n)\). For small \(s\),
\(R(\ell_0+s)\ge R_0+\tfrac12 c\sqrt s\), whereas
\(D^{3/4}-1\le 1/n\). Once \(n\) is large enough these
give (8) strictly, and arbitrarily small positive \(w\)
preserve it. The threshold is computable from the algebraic
coalescence point of (6) and a quantitative interval for (7).
For instance the coalescence parameters are roots of the
nonzero discriminant of (6). Polynomial root ordering at
rational parameters on either side and at the algebraic
candidates isolates \(\ell_0,R_0\). A rational interval for
the stated lower bound in (7) can be certified by isolating
the simple local Hermite branch after \(s=v^2\).

For \(m=4,5\), (8) holds at every admissible order even at
\(w=0\) in the test, and the crossing is transverse for each
positive \(w\) in the stated slack interval. In fact the
outer imaginary parts squared of (4) are
\[
 \begin{array}{ll}
 6a_j+\sqrt{24}\sqrt{a_j^2-b_j},&m=4,\\
 10a_j+\sqrt{40}\sqrt{a_j^2-3b_j},&m=5 .
 \end{array}
\]
This follows by writing the two polynomials as a quadratic
in \(X^2\), with an extra \(X\) factor for \(m=5\).
Thus \(\ell_0=1,\sqrt3\), respectively. Put
\(C=6,10\) and \(D_0=\sqrt{24},\sqrt{40}\) in the respective
cases. At \(\delta_w=0\) inequality (8) after squaring
reduces to \(D_0/\sqrt N>C/n\), which holds for \(n\ge2m\).
That last bound follows from nonreal conjugation and the
lower multiplicity \(m\). The outer branch is increasing
in \(\ell\), so every positive \(\delta_w\) preserves the
initial inequality. For \(y\) above both thresholds the
derivative of the lower squared imaginary part is strictly
larger than that of the upper. Indeed \(a'_n>a'_N\) and
\[
 \frac{\ell_0^2 b_n}{a_n(y)^2}>
 \frac{\ell_0^2 b_N}{a_N(y)^2}
\]
there, since these ratios are
\(\ell_0^2/[3jw(y-1/(2j))^2]\). Differentiate the displayed
quadratic roots to obtain the derivative inequality. The
eventual sign is positive by (9), giving the claimed
transverse crossing with no variation of \(w\).

## Physical attainment of each perturbed crossing

At a transverse limiting crossing in (5) or (9), take a
small real interval whose two endpoints have opposite phase
signs and on which the limiting derivative is bounded away
from zero. The continued exact branches, divided by the
angle scale \(\sqrt{A_M}\), preserve both conditions for
large \(M\). Their simple angular zero gives two roots of
the exact reference with equal argument. The root values
approach the given root coordinates because their
logarithmic displacements are \(O(\sqrt{A_M})\).

The added pair has total signed gap zero and its center corrections
cancel exactly by reflection. Truncation preserves the center in the
limit. Keep the total \(\Delta_0\) fixed, adjusting the residual
Gaussian signed gap relative to each finite retained core, and set
\[
       \lambda_M=2nN\bigl(\Gamma_{\rm fin}
                   -\log|z_{N,M}/z_{n,M}|\bigr).
\]
This gives the correct positive modulus ratio and
\(\lambda_M\to\lambda_*\). In the first case \(w_M\to0\). In
the fixed variance cases \(w_M\to w<T_*-|\lambda_*|\).
Thus the new Gaussian strength \(T-w_M\) and the bias are
in the strict interior for large \(M\). Transversality
with respect to \(y\) is equivalent to transversality in
the new Gaussian strength at each fixed \(M\). Finally,
the midpoint and upper-exponent curvatures of the added
rows differ in sum by \(O(w_M/M)\). Compensation preserves
the total curvature limit and hence every strict original cap.

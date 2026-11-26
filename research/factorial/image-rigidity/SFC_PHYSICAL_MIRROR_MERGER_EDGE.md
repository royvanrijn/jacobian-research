# Uniform angular confinement and the single-block merger limit

Consider independent mirrored blocks of common-shape correlated
Gamma pairs. In block \(i\) write \(d_i\) for the positive integer
gap and \(\rho_i\) for the squared correlation. At every degree,
the combined moment roots have circular distance from \(-1\)
at most
\[
                2\sum_i d_i\arcsin\!\sqrt{1-\rho_i}.
                                                               \tag{1}
\]
A block at complete correlation contributes zero to the bound.
The comparison is uniform in the degree and the shared exponents.
For one ordinary block, all positive gaps and shared exponents
have a strictly interlacing Gegenbauer limit at the merger,
after rescaling the angle by \(\sqrt{1-\rho}\).

## Arc bounds under sampling and multiplication

Let \(F(w)=\sum_{j=0}^{dm}\binom{dm}{j}f_j w^j\) have circular
roots lying within an angular distance \(\eta\) of \(-1\).
Its binomial-normalized subsampling
\[
                    \sum_{k=0}^m\binom mk f_{dk} z^k
\]
has no unit-circle root at angular distance greater than
\(d\eta\) from \(-1\). To see this, use the
[circular subsampling identity](SFC_PHYSICAL_MIRRORED_ROOT_CIRCLE.md#passing-to-arbitrary-gaps-and-multiple-independent-pairs).
At a proposed root \(z\) on the unit circle, the evaluation is
a nonzero scalar times the polarization of \(F\) at the
points \(\omega=-1/\zeta\) with \((1+z\zeta^d)^m=0\).
These are precisely the \(d\) solutions of
\(\omega^d=(-1)^{d+1}z\), repeated \(m\) times.
The least angular distance of these points from \(-1\)
equals the angular distance of \(z\) from \(-1\), divided
by \(d\). If this exceeds \(\eta\), an open disk centered
at a positive real number contains every polarization point
and excludes the closed zero arc of \(F\). Indeed the
distance to such a center strictly decreases with angular
distance from \(-1\) on the unit circle. Apply the
Grace–Walsh–Szegő coincidence theorem in that disk to
exclude a zero value.

For two degree-\(m\) circular polynomials with zero arcs of
angular radii \(\eta_1,\eta_2\) about \(-1\), the circular
Schur–Szegő composition has arc radius at most
\(\eta_1+\eta_2\). At a proposed unit root \(z\) farther away,
each polarization point \(-z/\zeta\), where \(\zeta\) is a
root of the second polynomial, has angular distance at most
\(\eta_2\) from \(z\). All points therefore lie farther than
\(\eta_1\) from \(-1\), and the same disk argument applies.
If the radius bound is at least \(\pi\), the claim is
automatic. The unit-gap Gamma formula places all those
root angles within \(2\arcsin\sqrt{1-\rho_i}\) of \(-1\),
since the real angular variable before multiplication by
\(\sqrt{1-\rho_i}\) has modulus at most one. Subsampling
then multiplies each allowable width by \(d_i\).
The composition rule and the
[physical circle-location theorem](SFC_PHYSICAL_MIRRORED_ROOT_CIRCLE.md)
give (1).

## The angular scaling limit for an ordinary block

Take independent standard circular complex Gaussians \(Z,W\)
with \(\mathbb E|Z|^2=1\). Put \(X=|Z|^2\) and
\[
                 Y_\delta=|\sqrt{1-\delta}\,Z+\sqrt\delta\,W|^2,
                 \qquad 0<\delta\le1 .
\]
For an ordinary pair with smaller exponent \(a\ge0\) and
gap \(d\ge1\), put \(b=2a+d\) and
\[
 U_\delta=(XY_\delta)^aX^d,\qquad
 V_\delta=(XY_\delta)^aY_\delta^d,\qquad
 P_m(z;\delta)=\mathbb E[(U_\delta+zV_\delta)^m].
\]
For every fixed integer \(m\ge1\), locally uniformly for
complex \(u\),
\[
 \frac{P_m(-e^{i\sqrt\delta\,u};\delta)}
                      {(-i)^m\delta^{m/2}}
   \ \longrightarrow\
       m!d^m\,\Gamma(1+(b-1)m)
                C_m^{\,1+(b-1)m}\!\left(\frac{u}{2d}\right).
                                                               \tag{2}
\]
Here \(C_m^\lambda\) is the Gegenbauer polynomial.

To prove this, expand the squared modulus in \(Y_\delta\).
In every finite \(L^p\),
\[
 \frac{Y_\delta-X}{\sqrt\delta}
          \ \longrightarrow\ 2\Re(\bar Z W)
                   =\sqrt X\,G ,
\]
where \(G\) is real normal with variance two and independent
of \(X\), by the rotational symmetry of \(W\).
All exponents in \(U_\delta,V_\delta\) are integers. Their
difference, divided by \(\sqrt\delta\), tends in \(L^p\) to
\(-dX^{b-1/2}G\), by factoring \(X^d-Y_\delta^d\).
The term from \(1-e^{i\sqrt\delta\,u}\) tends to \(-iuX^b\).
The convergence is uniform in bounded complex \(u\) in
each finite \(L^p\), by the exponential expansion and
Hölder's inequality for the finitely many products.
It therefore passes to the order-\(m\) moment.

After division by \((-i)^m\), only the even moments of \(G\)
remain. As \(\mathbb E G^{2k}=(2k)!/k!\), the limit is
\[
       \sum_{k=0}^{\lfloor m/2\rfloor}
       \frac{(-1)^k m!\,d^{2k}u^{m-2k}\Gamma(1+bm-k)}
                                      {k!(m-2k)!}.
\]
Factoring \(\Gamma(1+(b-1)m)\) identifies exactly the
Gegenbauer sum on the right side of (2).

The simple real limiting zeros capture all \(m\) roots
in the scaled angle, by Rouché on small disjoint disks,
or by the exact circle location and uniform argument
convergence. The consecutive limiting orders strictly
interlace. For \(b=1\) these are one fixed family.
For \(b>1\), apply the
[adjacent Gegenbauer parameter comparison](MULTIVARIATE_MIRRORED_ROOT_CIRCLE.md#consecutive-interlacing-for-a-single-mirrored-unit-gap)
with \(\lambda=1+(b-1)n\) and increment \(b-1\).
Thus for each fixed consecutive pair of orders there is
a nonempty interval of correlations below one throughout
which the original moments strictly interlace. A sequence of
physical common roots approaching complete correlation for
this fixed single block would have to use unbounded orders.
No order-independent size of that interval is asserted here.

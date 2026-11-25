# The subcritical private half-Gaussian window at escaping fields

Let fixed ordinary factorial profiles \(a,b\) have exactly one
private side, with a unique private linear variable in \(X^b\).
After removing equal rows, write \(d_i=b_i-a_i\),
\(\Delta=\sum_i d_i\), and
\[
 K_a=\prod_{a_i>0}a_i^{d_i},\qquad
 \sigma_a=\sum_{a_i>0}d_i^2/a_i<1 .
\]
Then the actual consecutive moment polynomials are coprime
for all sufficiently large orders.

The [Rice activation argument](MULTIVARIATE_SUBCRITICAL_PRIVATE_ENDPOINT_ACTIVATION.md)
reduces this to excluding common roots with
\(\beta=\log(n^\Delta t K_a)\to0\) and
\(n^{1/2}|\beta|\to\infty\).
The [bounded half-Gaussian windows](MULTIVARIATE_PRIVATE_SUPPORT_SADDLE_GEOMETRY.md#continuous-endpoint-windows)
cover bounded \(n^{1/2}|\beta|\).
Here an exact endpoint-preserving analytic coordinate converts
the local phase to a quadratic. Comparing the endpoint
integrands at the two orders gives a smaller error than
subtracting the two Watson estimates.

## Exact local normalization

Write \(q\) for the conditioning coordinate from the
[private saddle geometry](MULTIVARIATE_PRIVATE_SUPPORT_SADDLE_GEOMETRY.md).
At the \(a\) boundary the single private row has \(x_0=q\),
and all other \(x_i=a_i+d_iq\) are positive at zero.
In the conditional integral use \(q\ge0\), and subtract
the pure boundary phase \(A_a=\sum_{a_i>0}(a_i\log a_i-a_i)\).
The resulting phase near zero is
\[
 \Phi(q,\beta)=
  \log\!\left(1+e^\beta q
        \prod_{a_i>0}(x_i/a_i)^{d_i}\right)
       +\sum_{a_i>0}a_i\log x_i-\sum_i x_i-A_a .
 \tag{1}
\]
It is analytic on a complex bidisk about \(q=\beta=0\).
Put \(c=1-\sigma_a>0\). Then
\[
 \Phi(0,\beta)=0,\quad
 \Phi_q(0,\beta)=e^\beta-1,\quad
 \Phi_{qq}(0,0)=-c .
\]
The local stationary solution \(q_s(\beta)\) is analytic
with \(q_s'(0)=1/c\). Its value
\(S(\beta)=\Phi(q_s(\beta),\beta)\) has
\(S'(\beta)=q_s(\beta)\), because at stationarity
the derivative of the log factor with respect to \(\beta\)
is \(q_s\).

Define, by division at its removable singularity,
\[
 G(q,\beta)=2\,\frac{S(\beta)-\Phi(q,\beta)}
                         {(q-q_s(\beta))^2},\qquad
 \tau(\beta)=q_s(\beta)\sqrt{G(0,\beta)},\qquad
 z=\tau(\beta)+(q-q_s(\beta))\sqrt{G(q,\beta)} .
\]
Use the branch with \(\sqrt{G(0,0)}=\sqrt c\).
This is an analytic coordinate near zero:
\(z(0,\beta)=0,\ z_q(0,0)=\sqrt c\). It gives exactly
\[
  S=\tau^2/2,\qquad \Phi=-z^2/2+\tau z,\qquad
  \tau(\beta)=\beta/\sqrt c+O(\beta^2).                     \tag{2}
\]
The coordinate has real coefficients on real inputs.
On a sufficiently small disk \(\Re\tau'(\beta)>0\);
integration along a vertical line then shows that
\(\Im\tau\) has the sign of \(\Im\beta\).

After multiplication by a nonzero scalar depending on
\(n\), the moment is
\[
 J_n(\beta)=\int_0^L e^{n(-z^2/2+\tau z)}A_n(z,\beta)\,dz
                      +O(e^{-c_0 n}),                    \tag{3}
\]
for a fixed \(L>0\) and small complex \(\beta\).
The amplitudes on a bidisk containing \(|z|\le2L\) have
\[
 A_n=A_0+A_1/n+O(n^{-2}),\qquad A_0(0,0)>0,                 \tag{4}
\]
with uniform errors and derivatives on smaller bidisks.

Here are details about the actual endpoint in (3).
In the conditional coordinates \(y_i=x_i(q)e^{u_i}\),
\(d\cdot u=0\), the factor \(s'(q)\prod_i x_i(q)\)
is analytic and nonzero at zero. The transverse Hessian
\(\operatorname{diag}(x_i(q))\) restricted to this
hyperplane is positive at zero, since a vector supported
only in the single private coordinate is not in the
hyperplane. Restrict the transverse integral to a fixed
small symmetric ball in \(u\). It is then holomorphic
on a complex \(q\)-disk through zero, and has the uniform
Gaussian expansion through order \(1/n\) with relative
\(O(n^{-2})\) error, as in the
[conditional proof](MULTIVARIATE_SHARED_SUPPORT_CONDITIONAL_SADDLE.md#an-exact-conditional-integral).
Indeed on this ball the real part of its quadratic
form is uniformly positive for small complex \(q\);
odd scaled terms integrate to zero.
The original real integral outside a fixed neighborhood
of \(q=0,u=0\) has a strict phase gap at \(\beta=0\):
the pure boundary is then the unique maximum, by
\(\sigma_a<1\) and the real phase derivative.
Compact upper level sets and the one-order integral
bound make this contribution exponentially small
uniformly for nearby complex \(\beta\). Thus the
truncated transverse amplitude is sufficient, without
analytic continuation of the full transverse integral
to \(\Re q<0\).

Use (2) on a fixed small real \(q\)-segment starting
at zero. Its other endpoint maps to a positive real
point at \(\beta=0\), and remains close to it for complex
\(\beta\). Deform the image segment to \([0,L]\)
followed by a short connector to that endpoint, with
\(L\) fixed smaller than the endpoint value at zero.
On the connector \(\Re(-z^2/2+\tau z)\) stays strictly
negative, so its integral is exponentially small.
The amplitude and coordinate Jacobian give (3)–(4).
All scalar normalizations used to obtain (3) are
nonzero and do not change its zeros.

## The escaping estimates

By conjugation assume \(\Im\beta\ge0\), and consider
\(\beta\to0\) with
\[
 M=n|\tau|^2\longrightarrow\infty .
\]
There are two convenient sectors.

If \(7\pi/16\le\arg\tau\le\pi\), rotate the segment of
integration in (3) to angle \(\pi/8\).
There \(\Re(\tau z)\le-c_1|\tau||z|\) and
\(\Re z^2>0\), with a fixed \(c_1>0\).
The outer connecting arc makes an exponentially small
contribution, because \(\Re z^2\) is uniformly positive
on that arc and \(\tau\to0\).
Replacing the Gaussian factor by one and the amplitude
by \(A_0(0,\beta)\) gives
\[
 J_n(\beta)=-\frac{A_0(0,\beta)}{n\tau}
       \left(1+O\left(\frac1{n|\tau|^2}
                       +\frac1{n|\tau|}+\frac1n\right)\right)
                  +O(e^{-c_2 n}).
 \tag{5}
\]
Indeed \(|e^{-nz^2/2}-1|\le n|z|^2/2\) on the ray
and amplitude variation is \(O(|z|+1/n)\).
Their integrals against the linear exponential yield
the stated relative errors. Extending the linear
comparison to infinity costs an exponentially small
tail relative to \(1/(n|\tau|)\), since \(n|\tau|\to\infty\).
This sector therefore has no zeros on the escaping scale.

In the other sector \(\Re\tau\ge c_3|\tau|>0\).
Write the integral in (3) as the full segment
\([-L,L]\) minus the segment \([-L,0]\).
Shift the full segment parallel to the real axis through the
saddle \(z=\tau\), with exponentially small vertical
connectors. This gives
\[
 J_n(\beta)=R_n(\beta)-E_n(\beta)+O(e^{-c_4 n}),\qquad
 R_n=\frac{e^{nS}}{\sqrt n}\,R_0(\beta)(1+O(n^{-1})),
 \qquad R_0(\beta)=\sqrt{2\pi}\,A_0(\tau,\beta),
 \tag{6}
\]
where \(R_0\) is analytic and nonzero on a small disk.
This is a Gaussian estimate with a fixed analytic
amplitude neighborhood, so its relative error
\(O(n^{-1})\) is uniform as \(\tau\to0\).
The endpoint term is exactly
\[
 E_n(\beta)=\int_0^L e^{-n\tau x-nx^2/2}
                            A_n(-x,\beta)\,dx.
\]
Watson comparison as in (5) gives
\[
 E_n(\beta)=\frac{A_0(0,\beta)}{n\tau}(1+o(1))
 \qquad(n|\tau|^2\to\infty),                               \tag{7}
\]
uniformly in this sector. In particular its modulus is
bounded below by \(c_5/(n|\tau|)\).

## Adjacent comparison of the endpoint integrals

At a fixed original \(t\) the consecutive scaled
parameter is
\(\beta'=\beta+\Delta\log(1+1/n)\); put
\(\tau'=\tau(\beta')\). On the escaping scale
\(\tau'-\tau=O(1/n)=o(|\tau|)\), so the second order
remains in a sector with \(\Re\tau'\ge c_6|\tau|\)
whenever the first has \(\Re\tau\ge c_3|\tau|\).
The difference of the two endpoint integrals is
estimated directly, rather than by the separate
error terms in (7).

For \(x\ge0\), interpolate linearly between the two
exponents in the endpoint integrands.
Along that interpolation their moduli are bounded by
\(e^{-c_7 n|\tau|x-nx^2/2}\).
The change in the exponents is bounded by
\[
 |(n+1)\tau'-n\tau|x+x^2/2
              \le C(|\tau|+|\Delta|)x+x^2/2 .
\]
By (4), the change in the amplitude is \(O(n^{-1})\)
if \(\Delta\ne0\), and \(O(n^{-2})\) if \(\Delta=0\).
In the second case the parameters agree so the common
\(1/n\) coefficient cancels to that accuracy.
Integrating these bounds and using (7) gives
\[
 \frac{E_{n+1}(\beta')}{E_n(\beta)}
   =
 \begin{cases}
  1+O(1/(n|\tau|)),&\Delta\ne0,\\[2mm]
  1+O(1/n),&\Delta=0 .
 \end{cases}                                               \tag{8}
\]
In detail the relative moment bounds contribute
\(O((|\tau|+|\Delta|)/(n|\tau|)
       +(n|\tau|)^{-2})\), and
\((n|\tau|)^{-2}=1/(nM)\).
The first bound is \(o(|\tau|)\) and the second is
\(o(|\tau|^2)\) on the escaping scale.

If both orders vanished, then (6)–(7) would give equal
saddle and endpoint ratios, up to exponentially small
relative errors. On the saddle side (6) gives
\[
 \frac{R_{n+1}(\beta')}{R_n(\beta)}
    =\exp\!\left(S+\Delta q_s+O(n^{-1})\right)
                          (1+O(n^{-1})).
\]
Here the increment is the analytic Taylor formula
using \(S'=q_s\), and the analytic Gaussian
amplitudes change by \(O(n^{-1})\).
When \(\Delta\ne0\), \(q_s\) is a nonzero fixed
multiple of \(\tau\) to first order, so the ratio
differs from one by order \(|\tau|\). This contradicts
(8). When \(\Delta=0\), the ratio differs from one
by order \(|S|=|\tau|^2/2\), again contradicting (8).
At a zero all exponentially small errors in (6)
are negligible relative to (7), since \(1/(n|\tau|)\)
is only polynomially small. This rules out all
escaping fields, including angles near the sector
division, where either the nonzero estimate (5)
or the overlapping Gaussian decomposition is valid.

The fixed scaled annulus, the compact complex and
regular real exclusions, the bounded scaled window and
this escaping estimate prove eventual coprimality.
The polynomials have rational coefficients, so the
result transfers to every characteristic-zero field.

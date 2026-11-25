# The critical private cubic endpoint at escaping fields

Let fixed ordinary factorial monomials have exactly one
private side, with one private linear variable in \(X^b\).
Remove equal rows, put \(d_i=b_i-a_i\), and set
\[
 \Delta=\sum_i d_i,\qquad
 K_a=\prod_{a_i>0}a_i^{d_i},\qquad
 \sum_{a_i>0}d_i^2/a_i=1 .
\]
Then the actual consecutive moment polynomials are
coprime for all sufficiently large orders.

The [one-private Rice argument](MULTIVARIATE_ONE_PRIVATE_SIDE_RICE_EVENTUAL.md)
reduces this to excluding common zeros approaching
\(\zeta=n^\Delta t=1/K_a\) with
\(n^{2/3}|\log(\zeta K_a)|\to\infty\).
The [bounded cubic windows](MULTIVARIATE_PRIVATE_SUPPORT_SADDLE_GEOMETRY.md#continuous-endpoint-windows)
cover the complementary scaled fields.
The escaping estimate below uses an explicit
cubic contour and compares adjacent integrands on the
same contour.

## The exact local integral

Set \(\xi=\log(\zeta K_a)\) near the positive endpoint.
In the conditioning coordinate \(q\), subtract the
pure-boundary phase \(A_a=\sum_{a_i>0}(a_i\log a_i-a_i)\).
The exact local phase is
\[
 \Phi(q,\xi)=
  \log\!\left(1+e^\xi q
       \prod_{a_i>0}(1+d_iq/a_i)^{d_i}\right)
       +\sum_{a_i>0}a_i\log(a_i+d_iq)
       -\sum_i(a_i+d_iq)-A_a .
\]
If \(c>0\) is the \(q^2\) coefficient of
\(\log(\zeta(q)K_a)\), where \(\zeta(q)\) is the
stationary map, then
\[
 \begin{gathered}
 \Phi(0,\xi)=0,\qquad \Phi_q(0,\xi)=e^\xi-1,\\
 \Phi(q,\xi)=\xi q-cq^3/3+
                 O(|q|^4+|\xi|\,|q|^2+|\xi|^2|q|).
 \tag{1}
 \end{gathered}
\]
In fact \(c=-\rho''(0)/4>0\) in the
[curvature calculation](MULTIVARIATE_PRIVATE_SUPPORT_SADDLE_GEOMETRY.md#continuous-endpoint-windows).
After multiplication by a nonzero normalization
depending on \(n\), the actual moment has the form
\[
 I_n(\xi)=\int_0^L e^{n\Phi(q,\xi)}A_n(q,\xi)\,dq
                    +O(e^{-c_0n}),\qquad
 A_n=A_0+A_1/n+O(n^{-2}),\quad A_0(0,0)>0.                \tag{2}
\]
The amplitudes and errors are holomorphic and uniform on
a fixed bidisk containing \(|q|\le2L\).
To justify this representation, truncate to small
\(q\ge0\) and a small symmetric ball of the transverse
real coordinates. At \(\xi=0\) the pure boundary is the
unique real global maximum, so its complement has a
strict gap on the real orthant, persisting for nearby
complex \(\xi\) by the one-order integral upper bound.
The transverse Hessian is positive at zero even though
the private coordinate vanishes, because that coordinate
alone is not in the transverse hyperplane.
The bounded transverse integral is holomorphic through
\(q=0\), and its ordinary uniform Laplace expansion
gives (2) after dividing out the leading transverse
power of \(n\). This is the same cutoff argument as for
the [quadratic boundary](MULTIVARIATE_HALF_GAUSSIAN_ENDPOINT_ESCAPING_FIELDS.md#exact-local-normalization).
It does not extend the untruncated transverse contour
to negative real \(q\).

We take \(\Im\xi\ge0\) by conjugation.
Write \(\alpha=\arg\xi\in[0,\pi]\),
\(\rho=|\xi|^{1/2}\), \(\lambda=c^{-1/2}\), and
\(M=n\rho^3\). In the escaping range \(\rho\to0\),
\(M\to\infty\). The local saddle \(q_s\) is analytic
in \(\sqrt{\xi}\), with
\[
 q_s=\lambda\sqrt{\xi}+O(\xi),\qquad
 S(\xi)=\Phi(q_s,\xi)=\tfrac23\lambda\xi^{3/2}+O(\xi^2),
 \qquad S'(\xi)=q_s .                                      \tag{3}
\]
This follows by locally inverting
\(\log(\zeta(q)K_a)=c q^2+O(q^3)\) using an analytic
square root, and then applying (1) and stationarity.
The inverse has real coefficients near zero and
positive derivative there. In particular \(q_s\)
lies on or above the real line in the sector where
the saddle contour below is used.

## The endpoint sector

For \(\alpha\ge3\pi/8\), rotate the finite real
segment in (2) to angle \(\pi/7\). The real parts of
both the linear term and the negative cubic on that
ray have a uniform negative margin:
\(\cos(\alpha+\pi/7)<0\) and \(\cos(3\pi/7)>0\).
For \(L\) small and then \(\xi\) small, (1) gives
\[
 \Re\Phi(e^{i\pi/7}x,\xi)\le
             -c_1(|\xi|x+x^3),\qquad 0\le x\le L .
\]
On the outer connecting arc, the real part is bounded
away from zero negatively, since the negative cubic
dominates uniformly there. The rotated endpoint
integral has the Watson estimate
\[
 -\frac{A_0(0,\xi)}{n(e^\xi-1)}(1+o(1)).                  \tag{4}
\]
Indeed after subtracting its exact linear term
\((e^\xi-1)q\), the phase has error
\(O(|\xi|\,|q|^2+|q|^3)\). Interpolating between the
exact and linear exponents on the ray keeps uniform
linear decay. The integrated relative error from the
quadratic and cubic bounds is respectively
\(O((n|\xi|)^{-1})\) and \(O(M^{-2})\), by integrating
against \(e^{-c_1 n|\xi| x}\).
Amplitude variation contributes
\(O((n|\xi|)^{-1}+n^{-1})\).
Extending the linear integral past \(L\) has
exponentially small relative cost since \(n|\xi|\to\infty\).
This proves (4) and excludes zeros in the endpoint
sector.

## An explicit two-contribution contour

Suppose \(0\le\alpha\le3\pi/8\).
By holomorphy replace \([0,L]\) in (2) by the sum
of the endpoint ray \(0\to L e^{2\pi i/3}\) and
a contour from \(L e^{2\pi i/3}\) to \(L\)
through \(q_s\). For the latter use the ray
from \(i\Im q_s+L e^{2\pi i/3}\) back to \(i\Im q_s\),
then the horizontal segment through \(q_s\) to
\(L+i\Im q_s\), with short connectors at both outer
ends. Those connectors lie in fixed regions of
strict negative real phase by (1).
The endpoint ray obeys the same Watson estimate
(4), because \(\cos(\alpha+2\pi/3)<0\) uniformly
on this sector and its cubic has negative real part.

Here are the saddle bounds on the second contour.
After \(q=\lambda\rho z\), the scaled phase
\(\Phi/(\lambda\rho^3)\) converges uniformly on compact
sets with derivatives to \(e^{i\alpha}z-z^3/3\).
Write the model saddle as \(z_s=x+iy=e^{i\alpha/2}\).
On the model horizontal segment the real phase is
\[
 x^2 v-v^3/3-2 x y^2 \qquad(z=v+iy,\ v\ge0),
\]
with strict maximum at \(v=x\).
On the offset ray \(z=iy+t e^{2\pi i/3}\), \(t\ge0\),
the real phase change from \(iy\) is
\[
 -\bigl(x^2/2+\sqrt3xy\bigr)t
       -(\sqrt3 y/2)t^2-t^3/3 .
\]
Thus all portions outside a fixed model neighborhood
of the saddle have a uniform strict gap on bounded
scaled sets. On unbounded scaled sets up to the fixed
\(q\)-radius, (1) bounds the real phase by
\(-c_2|q|^3+C|\xi||q|\); for \(|q|\ge R\rho\) and
\(R\) sufficiently large this gives a strict cubic
loss also relative to \(S\). Near the actual saddle
on the horizontal segment,
\(\Phi_{qq}(q_s)=-2c q_s+O(\rho^2)\), whose real part
is at most \(-c_3\rho\). On \(|q-q_s|\le c_4\rho\)
Taylor's formula gives a real loss at least
\(c_5\rho |q-q_s|^2\), with \(c_4\) fixed small.
These estimates persist when \(\alpha\) is at either
endpoint of its stated range.

Let \(R_n(\xi)\) denote this actual saddle-contour
integral, and \(E_n(\xi)\) the actual endpoint-ray
integral, both with the amplitudes in (2).
We have \(I_n=R_n+E_n+O(e^{-c_6n})\).
The local Gaussian expansion, on width
\((n\rho)^{-1/2}\), gives
\[
 R_n=e^{nS}\,A_n(q_s,\xi)
      \sqrt{\frac{2\pi}{-n\Phi_{qq}(q_s,\xi)}}\,(1+o(1)).
 \tag{5}
\]
For instance the relative cubic Taylor error is
\(O(M^{-1/2})\) by Gaussian domination on
\(|q-q_s|\le c_4\rho\).
The strict gap on the rest of the bounded scaled
contour and the cubic loss on the tails make their
integrals negligible by comparison to the Gaussian
width. In particular
\(|R_n|\ge c_7 e^{n\Re S}/\sqrt{n\rho}\).
The endpoint bound (4) gives
\(|E_n|\ge c_8/(n|\xi|)\).

## Direct comparison at the consecutive order

For fixed original \(t\), the next scaled parameter is
\(\xi'=\xi+\Delta\log(1+1/n)\).
In the escaping range \(\xi'-\xi=o(|\xi|)\), and the
same chosen contours can be used for the next-order
integrals by holomorphy.
Define the exponent change, as a function of the
integration variable, by
\[
 G_n(q)=(n+1)\Phi(q,\xi')-n\Phi(q,\xi).
\]
Taylor expansion in \(\xi\), using
\(\Phi_\xi,\Phi_{\xi\xi}=O(q)\) on the fixed local
disk, gives
\[
 G_n=\Phi+\Delta\Phi_\xi+O(|q|/n),\qquad
 G_n(q_s)=S+\Delta q_s+O(\rho/n).                          \tag{6}
\]
For \(\Delta=0\) the equality \(G_n=\Phi\) is exact.

On the endpoint ray the real exponent at both orders
has uniform linear and cubic loss. Interpolate between
the two exponents to bound their difference, using
\(|G_n(q)|\le C|q|\) if \(\Delta\ne0\) and
\(|G_n(q)|\le C(|\xi||q|+|q|^3)\) if \(\Delta=0\).
By the analytic amplitude expansion the amplitude
change is \(O(n^{-1})\) or \(O(n^{-2})\), respectively.
Integrate with the endpoint weight
\(e^{-c_9 n|\xi||q|}\) and divide by the lower bound
on \(E_n\). This gives
\[
 E_{n+1}(\xi')/E_n(\xi)=
 \begin{cases}
  1+o(\rho),&\Delta\ne0,\\
  1+o(\rho^3),&\Delta=0 .
 \end{cases}                                               \tag{7}
\]
Indeed the mean endpoint scale is
\((n|\xi|)^{-1}=\rho/M\), so the exponent errors
relative to the packet are \(O(\rho/M)\) and
\(O(\rho^3/M+\rho^3/M^3)\).

For the saddle ratio, evaluate both orders on the
same chosen saddle contour. On its bounded scaled
portion one has
\[
 |G_n(q)-G_n(q_s)|
      \le C(|\Delta|+\rho^2)|q-q_s|.
\]
Against the saddle Gaussian this has relative
mean \(O((|\Delta|+\rho^2)\rho/\sqrt M)\).
On the compact scaled part away from the saddle,
the same bound with the exponential gap contributes
\(O((|\Delta|\rho+\rho^3)\sqrt M e^{-c_{10} M})\).
On the unbounded scaled tails one may instead bound
the increment difference by
\(C[|\Delta|(|q|+\rho)+|q|^3+\rho^3]\).
The cubic loss gives the same little-oh relative
bounds there after \(q=\rho z\).
Since \(q\) stays in a fixed small disk and \(\Delta\)
is fixed, the exponential of the increment difference
minus one is bounded by a fixed multiple of the
increment difference. The amplitude-change estimates
above and the lower bound (5) now yield
\[
 \frac{R_{n+1}(\xi')}{R_n(\xi)}
 =e^{G_n(q_s)}
 \begin{cases}
  1+o(\rho),&\Delta\ne0,\\
  1+o(\rho^3),&\Delta=0 .
 \end{cases}                                               \tag{8}
\]
This argument compares exact contour integrands and
uses no derivative estimate for a bare asymptotic
remainder.

At a common zero the ratios (7) and (8) must agree
up to exponentially small errors, by
\(I_n=R_n+E_n+O(e^{-c_6n})\) and the polynomial lower
bound on the endpoint packet. For \(\Delta\ne0\),
(3) and (6) give a nonzero leading order \(\rho\)
increment \(\Delta q_s\), contradicting (7).
For \(\Delta=0\), the increment is \(S\), whose
modulus has order \(\rho^3\), again a contradiction.
This excludes all escaping sequences and completes
eventual coprimality for the critical private profiles.

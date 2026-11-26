# Compensated angular transport for roots off the circle

The positive Laplace increment of consecutive factorial symbols
also gives a finite root-matching estimate for some noncircular
references. A deterministic jump drift changes only the positive
radial scale. Subtracting it lets the estimate use a second
weighted jump moment even when the increment is asymmetric.

Two geometric hypotheses suffice below. One is a thin
logarithmic annulus. The other is cancellation between clockwise
and counterclockwise neighboring roots. The latter holds for
accurate quantile samples of smooth logarithmic contours and
allows radial width independent of the degree.

## Root geometry and the transfer bound

Let \(P(z)=\sum_{k=0}^m c_k z^k\) have degree \(m\ge2\), nonzero
roots, and angular mesh at least \(\gamma/m\), with
\(0<\gamma\le2\pi\). For each root \(z\) and each other root
\(\zeta\), put
\[
 w_\zeta=\log|z/\zeta|+i(\arg z-\arg\zeta),\qquad
 q_\zeta=\coth(w_\zeta/2),\qquad
 b_z=\frac1{2m}\sum_{\zeta\ne z}q_\zeta .
\]
The branch of the angle difference does not affect \(q_\zeta\).
Let \(M=\pi/\gamma\), \(D=2+4M^2\), and choose
\(0<a_0\le\min(1/32,1/(16M))\).

Use either of the following hypotheses and the constants
\(\mathcal C,\mathcal R,\mathcal B\) it supplies.

- For a thin annulus assume
  \(\max_{\zeta,z}|\log|z/\zeta||\le\alpha/m\). Set
  \[
   \mathcal C=M+1/(2m)+2\alpha M^2,\qquad
   \mathcal R=\alpha M^2,\qquad
   \mathcal B=M(1+\log m)+\alpha M^2 .
  \]
- For paired cancellation enumerate the other roots in circular
  order relative to each \(z\). Pair the \(j\)-th clockwise and
  \(j\)-th counterclockwise roots for \(j<m/2\).
  Suppose the sum of the magnitudes of their paired
  \(q_\zeta+q_{\zeta'}\), plus the magnitude of the lone
  \(q_\zeta\) in even degree, is at most \(A m\) at each root.
  Set
  \[
        \mathcal C=\tfrac12+\pi M/2+A,\qquad
                  \mathcal R=\mathcal B=A/2 .
  \]

In either case put
\[
  K_0=\mathcal C+a_0D/2,\qquad
  K_1=K_0+1/(2a_0)+\mathcal R,\qquad
  B=e^{4a_0^2D+a_0},\qquad b_0=\tfrac12e^{-4a_0^2D}.
\]
For a positive measure \(\nu\) on \(\mathbb R\setminus\{0\}\), let
\[
 \begin{split}
  T &= m\int h^2e^{K_1|h|}\,\nu(dh),\\
  E &=2B\mathcal B T+
       \frac{4B}{a_0}e^{(K_1^2+\mathcal B^2)T/2}
                                      (T+K_1^2T^2).
 \end{split}
\]
If \(T<\infty\) and \(2E/b_0\le a_0\), each root \(z\) of \(P\)
pairs with exactly one root of
\[
 Q(z)=\sum_{k=0}^m c_k z^k
 \exp\left(m\int
       (e^{h(k/m-1/2)}-1-h(k/m-1/2))\,\nu(dh)\right)
\]
of the form \(z e^{y/m}\), with \(|y|\le2E/b_0\).
In particular the angular displacement is at most \(2E/(m b_0)\).
For \(T=0\) the two polynomials agree. An affine change in the
log coefficient multiplier changes only the scalar and positive
radial scaling.

## The off-circle shift estimate

At a reference root define
\[
 F(w)=\frac{e^{-w/2}P(ze^{w/m})}{zP'(z)/m},\qquad
              \widetilde F(w)=e^{-b_z w}F(w).
\]
Factoring the polynomial gives
\[
 F(w)=2m\sinh(w/(2m))\prod_{\zeta\ne z} A_\zeta(w/(2m)),
 \quad A_\zeta(v)=\cosh v+q_\zeta\sinh v
              =\frac{\sinh(v+w_\zeta/2)}{\sinh(w_\zeta/2)}.
\]
Write \(s_\zeta=|\sin((\arg z-\arg\zeta)/2)|\).
For real \(v\), the modulus of
\(\coth(v+w_\zeta/2)\) is at most \(s_\zeta^{-1}\), and the
modulus of its derivative in \(v\) is at most \(s_\zeta^{-2}\).
These follow from
\(|\sinh(x+i\theta)|^2=\sinh^2x+\sin^2\theta\).
At the \(j\)-th circular neighbor
\(s_\zeta\ge\gamma\min(j,m-j)/(\pi m)\). In particular
\(\sum_{\zeta\ne z}s_\zeta^{-2}\le4M^2m^2\).

For the thin annulus, compare \(q_\zeta\) to
\(q_\zeta^0=-i\cot((\arg z-\arg\zeta)/2)\), which deletes its
log-radius difference. The same derivative bound gives
\(|q_\zeta-q_\zeta^0|\le|\Re w_\zeta|/(2s_\zeta^2)\).
Then \(|\Re b_z|\le\mathcal R\), \(|b_z|\le\mathcal B\).
For real \(v\) the difference of the two log shift moduli is
at most \(|\Re w_\zeta|\,|v|/(2s_\zeta^2)\), by integrating
the difference of the real coth terms in the log-radius.
Summing for \(v=h/(2m)\) costs at most \(\alpha M^2|h|\).
The [full sine-mesh estimate](MULTIVARIATE_MATCHED_LAPLACE_MESH_SHARPENING.md#the-exponential-rate-from-the-full-circular-mesh)
bounds the reference product with the \(q_\zeta^0\) by
\(e^{M|h|}\). The gauge costs at most another \(\alpha M^2|h|\).

For paired cancellation, the drift bounds follow directly.
A paired product satisfies
\[
 |A_\zeta(v) A_{\zeta'}(v)|
   \le\cosh^2v\, (1+|q_\zeta q_{\zeta'}|v^2)
                  \exp\big(|q_\zeta+q_{\zeta'}|\,|v|\big),
\]
by expanding the two factors and using \(|\tanh v|\le|v|\).
The lone factor in even degree is at most
\(\cosh v\,e^{|q_\zeta||v|}\).
With angular rank \(j\) in a pair, \(|q_\zeta q_{\zeta'}|
\le(Mm/j)^2\). The logarithm of the product of the factors
\(1+(Mm v/j)^2\) is at most \(\pi M m |v|\) by the
decreasing-integrand comparison on \((0,\infty)\).
Together with the gauge and the sinh prefactor, this gives the
paired growth rate \(\mathcal C\).

The complex perturbation is the same for both hypotheses.
For \(v=h/(2m)\), \(u=y/(2m)\), and \(|y|\le2a_0\),
\[
 \frac{A_\zeta(v+u)}{A_\zeta(v)}
       =\cosh u+q_\zeta(v)\sinh u,\qquad
 q_\zeta(v)=\coth(v+w_\zeta/2).
\]
Here \(|q_\zeta(v)|\le s_\zeta^{-1}\),
\(|q_\zeta(v)-q_\zeta|\le |v|s_\zeta^{-2}\), and
\(|q_\zeta(v)u|\le M a_0\le1/16\). Expanding the logarithm
about one gives an error at most
\(4(1+|q_\zeta(v)|^2)|u|^2\) after its linear term.
The factor \(e^{-b_z y}\) cancels the summed linear terms
\(q_\zeta u\). The remaining linear bound is
\(a_0D|h|/2\), and the sum of the quadratic bounds is at most
\(4a_0^2D\). The prefactor gives the extra constant \(e^{a_0}\)
for complex \(y\). This proves
\[
 |\widetilde F(h+y)|\le B(|h|+2a_0)e^{K_0|h|},
       \qquad |\widetilde F(y)|\ge b_0|y| .
\]
For the lower bound take \(h=0\), so the linear terms cancel
exactly, and use \(|2m\sinh(y/(2m))|\ge|y|/2\).
Cauchy's estimate and Taylor integration along a real segment
give, for \(|y|\le a_0\),
\[
 |\widetilde F'(y)|\le2B,\qquad
 |\widetilde F(y+h)-\widetilde F(y)-h\widetilde F'(y)|
       \le(2B/a_0) h^2 e^{(K_0+1/(2a_0))|h|}.
                                                               \tag{1}
\]

## Centering the positive jumps

First truncate \(\nu\) to \(|h|\ge\epsilon>0\).
Use the compound Poisson sum of intensity \(m\nu\) and subtract
its mean to obtain \(H\). Its log moment generating function is
\(\psi(t)=m\int(e^{t h}-1-t h)\,\nu(dh)\).
This is well defined for complex \(t\) with
\(|\Re t|\le K_1\). Put \(\chi=\mathbb E e^{b_z H}\).
The integral form of the exponential remainder gives
\[
     |\chi|\ge e^{-\mathcal B^2T/2},\qquad
             |\mathbb E[H e^{b_zH}]/\chi|\le\mathcal B T .
\]
For real \(t=\pm K_1\), the same remainder and its derivatives
give \(\psi(t)\le K_1^2T/2\), \(|\psi'(t)|\le K_1T\) and
\(\psi''(t)\le T\). Thus
\[
              \mathbb E[H^2 e^{K_1|H|}]
                \le2e^{K_1^2T/2}(T+K_1^2T^2).
\]
The truncated transformed polynomial in the log coordinate,
with the gauge removed from its prefactor, is
\(\mathbb E[e^{b_zH}\widetilde F(y+H)]\).
Divide by \(\chi\) and use (1) for its difference from
\(\widetilde F(y)\). The tilt in the remainder introduces at
most \(e^{\mathcal R|H|}\), which is included in \(K_1\).
The resulting error is at most \(E\).
Rouché's theorem on \(|y|=2E/b_0\) gives exactly one zero.
The compensated Laplace integral converges coefficientwise,
and \(\chi\) converges to \(e^{\psi(b_z)}\ne0\), so the
strict comparison passes to the limit. The angular disks
are disjoint because \(2a_0<\gamma\).

## Ordinary factorial increments and representative geometries

For \(N=n+1\), use the
[positive ordinary row increment](MULTIVARIATE_MATCHED_LAPLACE_NONLINEAR_TRANSFER.md#positive-curvature-increment-for-ordinary-shared-rows)
for \(L_N-L_n\). A row of smaller exponent \(\ell_i\) and
gap \(d_i\) gives rate \(a_i=\ell_i/d_i+1/2\).
Its Laplace curvature increment density at \(u>0\) is at
most \(u e^{-a_i u}/(2nN)\) at its signed Laplace coordinate.
As \(\nu(dh)=h^{-2}\mu(dh)\), if all \(a_i>K_1\) then
\[
                        T\le\frac{m}{2nN}
                                  \sum_i(a_i-K_1)^{-2}.
\]
For \(m=n\), this transports the actual lower symbol to the
upper symbol on the lower grid. If the profile is fixed and
the geometric hypotheses and a strict rate margin hold with
fixed constants along the orders, the angle bound is
\(O((1+\log n)/n^2)\) in the thin annulus case and
\(O(n^{-2})\) in the paired cancellation case. No equality
of the two column orientations is required.

In particular, the paired hypothesis holds for roots
\[
       z_{m,k}=\exp(u(k/m)+e_{m,k}),\qquad
             |e_{m,k}|\le C_0/m^2,\quad 0\le k<m,
\]
where \(u\) is twice continuously differentiable,
\(u(t+1)=u(t)+2\pi i\), and \(\Im u'(t)\) is bounded
below by a positive constant. For large \(m\) the angular
mesh follows from the derivative bound.
For the paired \(w_+,w_-\) at neighbor rank \(j\), the
second difference gives
\(|w_++w_-|\le(\|u''\|+4C_0)(j/m)^2\).
The identity
\[
      q_++q_-=\frac{\sinh((w_++w_-)/2)}
                    {\sinh(w_+/2)\sinh(w_-/2)}
\]
then bounds every paired sum by a constant, since the
denominator is bounded below by \((\gamma j/(\pi m))^2\).
This verifies the required sum bound.

A fixed annulus and angular mesh alone allow a large gauge
drift. For example, take odd \(m=2r+1\) and roots
\(-\exp(-a|j|/m+i\gamma j/m)\), \(-r\le j\le r\), with
\(a>0\) and \(0<\gamma\le\pi/2\).
The polynomial has positive coefficients by conjugate pairing
in the left half-plane. At its root \(-1\),
\[
 \Re q_j=\frac{\sinh(a|j|/m)}
                 {\cosh(a|j|/m)-\cos(\gamma j/m)}
            \ge C(a,\gamma)\frac{m}{|j|}\qquad(j\ne0).
\]
The inequality follows by bounding the numerator below
linearly and the denominator above quadratically in \(|j|/m\).
Thus \(\Re b_{-1}\) grows at least logarithmically with \(m\).
The additional geometric hypotheses above supply a uniform
drift bound. Establishing such geometry for a broad factorial
family and comparing the two different degree grids are
separate obligations.

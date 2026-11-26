# Polynomial angular transfer by inflation at fixed midpoint curvature

Let \(a,b\) be distinct nonnegative integer exponent columns. Remove
equal rows and write
\[
 d_i=b_i-a_i,\qquad m_i=(a_i+b_i)/2,\qquad A_i=m_i/|d_i|,\qquad
 W=\sum_i\frac{d_i^2}{m_i},\quad
 B_\infty=-\sum_i\frac{d_i^3}{m_i^2},\quad
 D=\sum_i\frac{d_i^4}{m_i^3}.
\]
For each degree \(j\ge2\), take the exact normalized midpoint clock
and set
\[
 t_j=\frac12\sum_i j d_i^2\psi'(1+jm_i),\quad
 \gamma_j=\min\{\pi,\sqrt{32t_j/\pi}\},\quad
 M_j=\frac{\pi}{\gamma_j},\quad
 r_j=\frac1{128M_j^2},\quad K_j=M_j+2 .
\]
Here \(\psi'\) is the trigamma function. When every \(A_i>K_j\), put
\[
 E_j=2r_j^{-3}(|B_\infty|+D/j)
       +4r_j^{-4}\sum_i\frac{|d_i|}{(A_i-K_j)^3}.
 \tag{1}
\]
If \(E_j<r_j/j\), then after a positive rescaling each exact
degree-\(j\) factorial root has a unique Gaussian reference root
\(\zeta\) at time \(t_j\) with
\[
                          |\log(z/\zeta)|\le E_j.                 \tag{2}
\]
The logarithm is the local branch at \(\zeta\). The exact roots
are simple, with no symmetry assumption on the row orientations.

For \(N=n+1\), assuming both one-degree bounds apply, the
ordinary pair has disjoint root arguments whenever
\[
                      E_n+E_N<\frac{t_n}{36nN(1+t_n)}.           \tag{3}
\]
In particular, for any \(w>0\) with \(W\ge w\), set
\(r_w=\min\{1/128,w/(16\pi^3)\}\). One sufficient condition,
with no additional row-rate assumptions, is
\[
              |B_\infty|+D<
                   \frac{r_w^4\,w}{2448 nN(w+4)}.               \tag{4}
\]
For an all-degree sufficient condition, use any computable
\(J\ge4\) at least the
[absolute ordinary angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md)
and replace \(nN\) by \(J^2\) in (4). The tolerance depends
polynomially on the cutoff. The proof of (3) applies directly
in each degree.

## A root chart along a short complex deformation

Fix \(j\ge2\) and a Gaussian reference with normalized time \(t>0\).
Drop subscripts on \(\gamma,M,r,K\) in this section. Its simple roots
\(\zeta_\alpha=e^{i\theta_\alpha}\) have circular mesh at least
\(\gamma/j\), and their signed angular forces satisfy
\[
 \left|\frac1{2j}\sum_{\beta\ne\alpha}
        \cot\frac{\theta_\alpha-\theta_\beta}{2}\right|
 \le\frac{\pi}{4t}.
\]
These are the individual mesh bound and the normalized heat velocity
bound. [Gaussian angular gaps and velocities](MULTIVARIATE_GAUSSIAN_LAPLACE_ANGLE_RESPONSE.md#angular-gaps-and-velocities-of-the-reference)

Suppose a polynomial \(P\) has roots
\(z_\alpha=\zeta_\alpha\exp(w_\alpha/j)\), with \(|w_\alpha|\le r\).
These logarithmic neighborhoods are disjoint by the circular mesh.
At a root \(z=z_\alpha\), define
\[
 Z(w)=\frac{e^{-w/2}P(ze^{w/j})}{zP'(z)/j}.
\]
Factoring the polynomial gives the exact formula
\[
 Z(w)=2j\sinh\frac{w}{2j}
      \prod_{\beta\ne\alpha}\left(
      \cosh\frac{w}{2j}-i q_\beta\sinh\frac{w}{2j}\right),
 \qquad q_\beta=\cot\left(\frac{\theta_\alpha-\theta_\beta}{2}
                    -\frac{i(w_\alpha-w_\beta)}{2j}\right).
 \tag{5}
\]
Write \(q_\beta^0=\cot((\theta_\alpha-\theta_\beta)/2)\).
At cyclic rank distance \(k\ge1\) from \(\alpha\), the sine of
the unperturbed half-gap has modulus at least \(k/(Mj)\).
On the segment joining the two cotangent arguments it has modulus
at least \(k/(2Mj)\), because the perturbation has length at most
\(r/j\), and the sine changes by at most \(2r/j\). Therefore
\[
 |q_\beta-q_\beta^0|\le4M^2rj/k^2,\qquad
                  \sum_{\beta\ne\alpha}|q_\beta-q_\beta^0|
                           \le16M^2rj .                         \tag{6}
\]
The sum uses at most two roots at each cyclic rank distance and
\(\sum_{k\ge1}k^{-2}\le2\).

Let \(A_q(s)=\cosh s-iq\sinh s\). The real-shift bound from
the full circular mesh is
\[
 \prod_{\beta\ne\alpha}|A_{q_\beta^0}(x/(2j))|\le e^{M|x|}
                 \qquad(x\in\mathbb R).
\]
It is obtained by summing the half-log moduli with cyclic distances
and integrating the monotone majorant
\(\log(1+\sinh^2 v/\sin^2(\gamma s/2))\).
[Full-mesh real-shift estimate](MULTIVARIATE_MATCHED_LAPLACE_MESH_SHARPENING.md#the-exponential-rate-from-the-full-circular-mesh)

For \(|y|\le r\) put \(v=x/(2j)\), \(u=y/(2j)\) and
\(b_q(v)=(\tanh v-iq)/(1-iq\tanh v)\), for real \(q\).
Then
\[
 \frac{A_q(v+u)}{A_q(v)}=\cosh u+b_q(v)\sinh u,\qquad
 |b_q(v)|\le\max(1,|q|),\qquad
 |b_q(v)+iq|\le|\tanh v|(1+q^2).
\]
For the reference factors \(|q|\le Mj\). In this range
\[
 \big|\log(\cosh u+b_q(v)\sinh u)-b_q(v)u\big|
                     \le4(1+|b_q(v)|^2)|u|^2 .
\]
Indeed the difference between the argument of the logarithm and
\(1+b_q(v)u\) has modulus at most \(2|u|^2\), and the argument
differs from \(1\) by less than \(1/8\).
Since \(\sum_\beta(1+(q_\beta^0)^2)\le6M^2j^2\), the
linear part of the sum of these logarithms is controlled by the
Gaussian force and \(|\tanh v|\le|v|\). Thus
\[
 \prod_{\beta\ne\alpha}|A_{q_\beta^0}(v+u)|
 \le \exp\left(M|x|+\frac{r\pi}{4t}
                    +\frac{6M^2r|x|}{4}+6M^2r^2\right).
\]
To insert (6), note that \(|A_{q^0}(v)|\ge\cosh v\) and
\(|A_{q^0}(v+u)|\ge |A_{q^0}(v)|/2\).
Hence \(|\sinh(v+u)/A_{q^0}(v+u)|\le 2(|x|+|y|)/j\).
The ratio in the product of moduli introduced by the perturbed
\(q_\beta\) is at most
\(\exp(32M^2r(|x|+r))\).
Finally \(2j|\sinh((x+y)/(2j))|\le
(|x|+r)e^{(|x|+r)/(2j)}\) and \(|x|+r\le e^{|x|}\).
Using \(r=1/(128M^2)\), \(j\ge2\), and
\(r\pi/(4t)\le1/(16\pi^2)\), these bounds give
\[
 |Z(x+y)|\le 2e^{K|x|}\qquad (x\in\mathbb R,\ |y|\le r).
 \tag{7}
\]
For example the coefficient of \(|x|\) before this simplification
is \(M+1+1/(2j)+(6/4+32)M^2r<M+2\).
The remaining exponent is less than \(1/8\).
By Cauchy's estimate
\[
                         |Z^{(k)}(x)|/k!\le2r^{-k}e^{K|x|}.
 \tag{8}
\]

## Interpolation by inflating all midpoints

Fix \(j\) and inflate the midpoint of each row to \(\lambda m_i\),
with \(\lambda\ge1\), preserving the signed gap \(d_i\). Use the
centered grid coordinate \(u=k/j-1/2\), and put
\[
 L_\lambda(u)=j^{-1}\sum_i
               \log\Gamma(1+j(\lambda m_i+d_i u)),\qquad
 C_\lambda=L_\lambda''(0).
\]
Subtract the midpoint tangent and add a centered Gaussian of
curvature \(C_1-C_\lambda\), giving the symbol
\[
 G_\lambda(u)=L_\lambda(u)-L_\lambda(0)-uL_\lambda'(0)
                               +(C_1-C_\lambda)u^2/2 .
\]
Let \(P_\lambda\) have coefficients
\(\binom jk\exp(jG_\lambda(k/j-1/2))\). At \(\lambda=1\)
this is the actual moment polynomial up to a positive coefficient
and argument rescaling. At \(\lambda=\infty\) its roots are the
Gaussian roots at normalized time \(t_j=C_1/2\),
since the row curvatures vanish uniformly on the grid.

The trigamma integral expresses \(L_\lambda''(u)\) as the
Laplace transform of a positive measure \(\mu_\lambda\).
For row \(i\) its measure is the pushforward under
\(h=-\operatorname{sign}(d_i)q\) of
\[
 |d_i|f(q/(j|d_i|))e^{-\lambda A_i q}\,dq,\qquad
                  f(s)=\frac{s}{e^s-1},\quad q>0 .
\]
This follows by differentiating the digamma integral and
substituting \(q=j|d_i|s\). [DLMF 5.9.16](https://dlmf.nist.gov/5.9.E16)
Put \(\nu_\lambda=\sum_i A_i|h|\mu_{\lambda,i}\).
Then \(\partial_\lambda\mu_\lambda=-\nu_\lambda\), and matching
both the midpoint tangent and the curvature gives
\[
 \partial_\lambda G_\lambda(u)=
   -\int\frac{e^{h u}-1-hu-(hu)^2/2}{h^2}\,\nu_\lambda(dh).
 \tag{9}
\]
The row-rate hypothesis makes these integrals converge uniformly
on the real coefficient grid.

For a simple root of \(P_\lambda\), use its chart \(Z\) from (5)
as long as the displacement from every Gaussian root is at most
\(r_j/j\) in complex log. Implicit differentiation of the polynomial
and (9) give
\[
 \frac{\partial_\lambda z}{z}=
   \int\frac{Z(h)-Z(0)-hZ'(0)-h^2 Z''(0)/2}{h^2}\,\nu_\lambda(dh).
 \tag{10}
\]
The factor \(j\) in differentiation is already included in the
root-chart normalization. Expand (10) through the cubic term.
Equation (8), with the fourth-order remainder on the real segment
from zero to \(h\), gives
\[
 \left|\frac{\partial_\lambda z}{z}\right|
 \le2r_j^{-3}\left|\int h\,\nu_\lambda(dh)\right|
        +2r_j^{-4}\int h^2e^{K_j|h|}\,\nu_\lambda(dh).
 \tag{11}
\]

Here the cubic moment is summed with its sign before taking the
absolute value. Since \(0\le1-f(s)\le s/2\), its difference from
\(2B_\infty/\lambda^3\) has absolute value at most
\(3\sum_i A_i^{-3}/(j\lambda^4)\).
Indeed, the leading integral in row \(i\) is
\(-2\operatorname{sign}(d_i)|d_i|/(A_i^2\lambda^3)\).
Replacing \(1\) by \(f\) gives the stated error by the integral of
\(q^3e^{-\lambda A_iq}\). The inequality for \(f\) follows from
\(f'\in[-1/2,0]\), by differentiating
\(f(s)=(s/2)\coth(s/2)-s/2\) and using \(\sinh s\ge s\).
Integrating in \(\lambda\) and using \(|d_i|\ge1\) gives
\[
 \int_1^\infty\left|\int h\,\nu_\lambda(dh)\right|d\lambda
                      \le|B_\infty|+D/j .
\]
Also \(\int_1^\infty\nu_\lambda(dh)d\lambda=\mu_1(dh)\).
Using \(f\le1\) bounds the integrated second term of (11)
by \(4r_j^{-4}\sum_i |d_i|/(A_i-K_j)^3\).
This proves the speed integral bound \(E_j\).

All roots start in the stated disjoint neighborhoods for sufficiently
large \(\lambda\), by convergence to the simple Gaussian roots.
Along continuation towards \(\lambda=1\), there is one root per
neighborhood until a first exit. A collision before exit is impossible
because the degree equals the number of neighborhoods and each carries
one root counted with multiplicity. At a first exit, integrating (11)
from infinity would give log-displacement at most \(E_j<r_j/j\),
a contradiction. This proves (2) and verifies the use of (8)
along the whole interpolation.

## Adjacent physical clocks and a degree bound

The exact midpoint curvatures obey \(C_N(1)>C_n(1)\) and
\(C_N(1)/N<C_n(1)/n\). The first inequality is the positive matched
curvature increment and the second follows from monotonicity of
\(\psi'\) in \(1+jm_i\).
[Exact midpoint clocks](MULTIVARIATE_MATCHED_MIDPOINT_CLOCK_ANGLE_TEST.md#the-gaussian-clocks-and-a-bound-on-root-rays)
Thus the reference pair uses upper raw time \(t_N/N\) in
\([t_n/N,t_n/n]\), where the cross-degree circular distance exceeds
\(t_n/[36nN(1+t_n)]\).
[Polynomial Gaussian physical gaps](MULTIVARIATE_GAUSSIAN_POLYNOMIAL_PHYSICAL_GAPS.md)
Taking the imaginary parts in (2) now gives (3).

For (4), the exact curvature bound \(C_j(1)\ge W/2\) yields
\(t_j\ge w/4\) and \(r_j\ge r_w\). Set \(M_w=(128r_w)^{-1/2}\)
and \(K_w=M_w+2\). The right side of (4) is less than
\((2K_w)^{-3}\), using \(M_w\ge1\). Since
\(D=\sum_i |d_i|/A_i^3\), every row must have \(A_i>2K_w\).
Now \(K_j\le K_w\), so
\(\sum_i |d_i|/(A_i-K_j)^3\le8D\).
It follows that \(E_j\le34r_w^{-4}(|B_\infty|+D)\) for
\(j=n,N\). Under (4) this is less than \(w/[72nN(w+4)]\),
which in turn is less than \(r_w/N\) for \(n\ge2\)
by \(r_w=\min\{1/128,w/(16\pi^3)\}\).
The sum of the two bounds satisfies (3).
At \(N=2\) the first moment has a negative real root and the
second moment is strictly positive on the real line, as the expectation
of the square of two distinct monomials with real coefficients.
With \(nN\) replaced by \(J^2\), the sufficient bound holds for
every remaining degree pair below the absolute cutoff.

## Direct matched tanh criteria

For the Gauss factor indexed by \(h\) in row \(i\), use the
rescaled midpoint parameter \(p_{j,i,h}=jT_{j,i,h}/(j-1)\)
and the sign \(\epsilon_i=\operatorname{sign}(d_i)\). Sums in
this section range over all the factors. Set
\[
 H_j=\sum p_{j,i,h},\qquad P_j^{\rm skew}=\sum\epsilon_i p_{j,i,h}^2,
 \qquad V_{k,j}=\sum p_{j,i,h}^k.
\]
The exact Gauss curvature kernel at \(u=s-1/2\) is
\[
 \int_0^\infty e^{-(1/(2p_{j,i,h})+\epsilon_i u)q}
                 {\cal D}(q/j)\,dq,\qquad
 {\cal D}(v)=\frac{v}{2\sinh(v/2)}.
\]
In particular \(C_j(1)\) is the sum of these kernels at \(u=0\).
One has \(0\le1-{\cal D}(v)\le v^2/24\).
[Exact matched Gauss kernel and bound](MULTIVARIATE_MATCHED_TANH_ORDINARY_ANGLE_LIFT.md#exact-curvature-of-the-assembled-symbol)

There is a direct version of (1) in these moments. Inflate each
individual Gauss midpoint rate \(A=1/(2p_{j,i,h})\) to \(\lambda A\)
in this integral, holding \({\cal D}\) fixed, then subtract the
midpoint tangent and restore the common midpoint curvature
exactly as in (9). The intermediate arrays are defined by these
kernels and agree with the factorial array at \(\lambda=1\).
The equation for the root velocity is still (11). The leading
signed cubic integral of a factor is
\(-2\epsilon_i/(A^2\lambda^3)\), and its error has magnitude at
most \(1/(j^2 A^4\lambda^5)\), by \(1-{\cal D}(q/j)\le q^2/(24j^2)\).
Thus the integrated cubic sum is bounded by
\(4|P_j^{\rm skew}|+4V_{4,j}/j^2\).
The integrated fourth-order bound is the same exponential moment
as above, with rate \(1/(2p_{j,i,h})\) and weight one for each factor.
Consequently, when \(2K_jp_{j,i,h}<1\) for every factor, one may
replace \(E_j\) in (2) and (3) by
\[
 E_j^{\rm tanh}=
  8r_j^{-3}\bigl(|P_j^{\rm skew}|+V_{4,j}/j^2\bigr)
 +32r_j^{-4}\sum_{i,h}
               \frac{p_{j,i,h}^3}{(1-2K_jp_{j,i,h})^3}.          \tag{12}
\]
The same sufficient bound \(E_j^{\rm tanh}<r_j/j\) ensures
continuation inside the logarithmic neighborhoods.

The matching directly controls the upper skew in (12). In the
notation \(\eta_{i,h}=2h/|d_i|-1\), the exact grid formula gives
\[
 p_{j,i,h}^{-1}=2A_i+\eta_{i,h}/j .
\]
For any two orders \(j_0,j\ge2\), let \(q_0=\max_{i,h}p_{j_0,i,h}\)
and \(\alpha=1-q_0|1/j-1/j_0|\). When \(\alpha>0\),
\[
 |P_j^{\rm skew}-P_{j_0}^{\rm skew}|
      \le2|1/j-1/j_0|\alpha^{-3}V_{3,j_0},\qquad
 V_{k,j}\le\alpha^{-k}V_{k,j_0}.                               \tag{13}
\]
Indeed \(p_{j,i,h}=p_{j_0,i,h}/
(1+\eta_{i,h}p_{j_0,i,h}(1/j-1/j_0))\) and \(|\eta_{i,h}|\le1\).
The exponent two uses the mean value theorem. In particular (13)
retains the \(1/(nN)\) factor and one higher power of the moment
when \(j_0=n,\ j=N\).
[Matched tanh transport](MULTIVARIATE_MATCHED_TANH_MOMENTS.md#exact-tanh-transport-and-moment-errors-at-every-order)

An alternative exact input for (1) uses
\(U_{i,h}=(n/T_{N,i,h}-(n-1)/T_{n,i,h})^{-1}=1/(2A_i)\).
It gives \(W=2\sum U_{i,h}\), \(B_\infty=-4\sum\epsilon_i U_{i,h}^2\),
and \(D=8\sum U_{i,h}^3\), without estimating an order drift.
[Exact row-grid extraction](MULTIVARIATE_MATCHED_TANH_INTEGER_MULTIPLICITY.md#recovering-the-gauss-grid)

For a single explicit sufficient test, take the moments at one
order \(j_0\ge2\), put \(w=8H_{j_0}/5\) and set \(r_w\) as in
(4). Then the target pair \(n,N=n+1\), \(n\ge2\), is separated if
\[
       |P_{j_0}^{\rm skew}|+V_{3,j_0}<
               \frac{r_w^4 w}{8000 nN(w+4)} .                  \tag{14}
\]
To check this, let \(M_w=(128r_w)^{-1/2}\) and \(K_w=M_w+2\).
The condition implies \(q_0<1/(16K_w)\). Indeed \(q_0^3\le V_{3,j_0}\)
and \((16K_w)^3r_w^4<1\), by \(K_w\le3M_w\).
From the reciprocal grid at \(j_0\),
\[
 W=2\sum_{i,h} \frac{p_{j_0,i,h}}
                       {1-\eta_{i,h}p_{j_0,i,h}/j_0}
                  \ \ge\ 8H_{j_0}/5=w .
\]
This gives \(r_j\ge r_w\) and \(K_j\le K_w\) for every \(j\ge2\).
Since \(|1/j-1/j_0|\le1/2\), take \(\alpha\ge7/8\).
Then \(p_{j,i,h}\le q_0/\alpha<1/(14K_w)\), so
\(1-2K_jp_{j,i,h}>6/7\). These two denominators
give a product at least \(3/4\), bounding the fourth-order term of
(12) by \(76r_w^{-4}V_{3,j_0}\).
Further, (13) gives
\(|P_j^{\rm skew}-P_{j_0}^{\rm skew}|<(3/2)V_{3,j_0}\), and
\(V_{4,j}/j^2\le V_{3,j_0}/112\) by
\(p_{j,i,h}<1/42\) and \(j\ge2\).
Thus \(E_j^{\rm tanh}\le92r_w^{-4}
(|P_{j_0}^{\rm skew}|+V_{3,j_0})\).
Under (14) this is less than \(w/[72nN(w+4)]<r_w/N\) for
\(j=n,N\), and their sum satisfies (3).

For an all-degree test from one set of tanh data, replace \(nN\)
in (14) with \(J^2\) for an absolute ordinary cutoff \(J\ge4\).
The known eventual result covers the remaining degrees.
The moment bounds need not make (14) small for unrestricted
ordinary rows. A comparison that controls the *difference*
of the two angular responses, including the shared nonquadratic
part of the two symbols, is needed to exploit the full
\(1/(nN)\) matching gain without a smallness hypothesis.

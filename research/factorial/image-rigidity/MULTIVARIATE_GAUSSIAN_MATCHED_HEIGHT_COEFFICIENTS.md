# Matched squared-height coefficients at Gaussian clocks

Let \(n\ge2\), \(N=n+1\), and use the Gaussian reference
\[
 R_j(z;\tau)=\sum_{k=0}^j \binom jk e^{-\tau k(j-k)}z^k .
\]
The reference has simple unit-circle roots for \(\tau>0\).
For a positive angular height \(x_{j,i}(\tau)\) from \(-1\), put
\(\lambda_{j,i}=\sin^2(x_{j,i}/2)\). Let \(b_{j,r}(\tau)\) be the
elementary symmetric function of order \(r\) in these squared heights.
Write \(t=j\tau\) for normalized time, so that the two raw times at
the smallest physical ratio are \(t/n,t/N\).
These definitions and the coefficient evolution below are from the
[physical-clock comparison](MULTIVARIATE_GAUSSIAN_PHYSICAL_CLOCK_COMPARISON.md).

## The linked positive integral

For \(0\le r\le\lfloor n/2\rfloor\) define
\[
 d_{j,r}(t)=\frac{j^r b_{j,r}(t/j)}{\binom{j}{2r}},\qquad
 c_{n,r}(t)=\frac{d_{N,r}(t)}{d_{n,r}(t)} .
\]
For every \(t>0\), \(c_{n,0}=1\) and \(c_{n,r}\) strictly decreases
with \(r\). The explicit bounds are
\[
  \exp\!\left(-\frac{t\,r(2r+1)}{6nN}\right)
        \le c_{n,r}(t)\le1.                                      \tag{1}
\]
In particular the exact elementary symmetric coefficients obey
\[
  \frac{b_{N,r}(t/N)}{b_{n,r}(t/n)}
     =\frac{N}{N-2r}\left(\frac nN\right)^r c_{n,r}(t).
                                                                  \tag{2}
\]
The statement uses only coefficients present in both degrees.

To prove this, the triangular squared-height evolution gives
\[
  \dot d_{j,r}
     =-p_{j,r}d_{j,r}+r(r-\tfrac12)d_{j,r-1},
     \qquad p_{j,r}=r(1-r/j).
\]
Here \(d_{j,0}=1\) and all positive-index initial values are zero.
For ordered times \(0<T_1<\cdots<T_r<t\), put
\(T_{r+1}=t\) and \(u_k=T_{k+1}-T_k\). Integrating the triangular
system gives
\[
 d_{j,r}(t)
  =r!\,(\tfrac12)_r\int_{0<T_1<\cdots<T_r<t}
       \exp\!\left(-\sum_{k=1}^r p_{j,k}u_k\right)\,d\boldsymbol T.
\]
As \(p_{N,k}-p_{n,k}=k^2/(nN)\), the ratio is the positive expectation
\[
 c_{n,r}(t)=
  \mathbb E_{n,r,t}
      \exp\!\left[-\frac1{nN}\sum_{k=1}^r k^2u_k\right],             \tag{3}
\]
where the probability density is proportional to the degree-\(n\)
integrand. In arrival coordinates its density is proportional to
\(\prod_{k=1}^r w_k(T_k)\) where
\(w_k(x)=\exp((p_{n,k}-p_{n,k-1})x)\).

For any positive weights of this product form, the ordered arrival
vector can be coupled to be coordinatewise increasing in the upper
boundary \(t\). Indeed the density of the last arrival \(s\) is
proportional to \(w_r(s)Z_{r-1}(s)\) on \(0<s<t\), where
\(Z_{r-1}(s)\) is the integral with the first \(r-1\) weights and upper
boundary \(s\). Increasing \(t\) only changes this truncation.
Use a common quantile for the two last arrivals, then couple the
previous arrivals by induction. Conditioning on the last arrival of
an \((r+1)\)-vector now couples its first \(r\) arrivals below an
\(r\)-vector with boundary \(t\). The cost in (3) equals
\[
 r^2t-\sum_{k=1}^r(2k-1)T_k.
\]
At length \(r+1\) it has the additional positive term
\((2r+1)(t-T_{r+1})\) when the first \(r\) arrivals are held fixed.
This proves the strict monotonicity.

All \(w_k\) are increasing for \(k\le\lfloor n/2\rfloor\).
The weighted arrival vector stochastically dominates the uniform
ordered arrivals by the same recursive coupling. For the last
arrival, the density ratio to the uniform ordered law is
proportional to \(w_r(s)Z_{r-1}(s)/s^{r-1}\), which is increasing.
Indeed, the ratio \(Z_{r-1}(s)/s^{r-1}\) is the integral of
\(\prod_{k<r}w_k(sv_k)\) over the unit ordered simplex.
After coupling the last arrivals, use boundary monotonicity and
induction on the conditioned earlier arrivals. The cost in (3) is
thus at most its uniform-simplex counterpart in this coupling.
The uniform mean cost is \(t \sum_{k=1}^r k^2/(r+1)=t r(2r+1)/6\).
Jensen's inequality gives (1).

## An interval of clock ratios with all-degree interlacing

For every \(n\ge2\) and every \(\tau>0\), the two Gaussian references
strictly interlace if the upper raw time equals \(\rho\tau\), where
\[
       \frac{C_n}{1+C_n}\le\rho\le1,\qquad
                 C_n=\min\{2n^3,8n^2\}.                            \tag{4}
\]
The inner inequalities already hold for all slower upper clocks.
The proof of the outer inequalities starts with strict interlacing
at small times by the scaled Hermite comparison, since
\(\rho\ge n/N\).

In the notation of the physical-clock comparison, let \(F_{j,i}\)
denote the raw squared-height speed of the \(i\)-th root. The
[speed bound](MULTIVARIATE_GAUSSIAN_PHYSICAL_CLOCK_COMPARISON.md#motion-and-comparison-with-the-equilibrium-roots)
at the lower degree gives
\[
 \frac{F_{n,i}}{1-\lambda_{n,i}}
 \le \frac{\lambda_{n,i}}{2(e^{\tau/2}-1)(1-\lambda_{n,i})}.
\]
Put \(H=4n\) and \(\mu=\cos^2(\pi/(2n))\).
The Hermite starting speed and the monotonicity of
\(\lambda_{n,i}/(1-e^{-\tau/2})\) give
\(\lambda_{n,i}\le H(1-e^{-\tau/2})\).
Here the positive Hermite roots have squares at most \(4n\) by
the row sums of their Jacobi matrix. The equilibrium gives
\(\lambda_{n,i}<\mu\). For \(v=1-e^{-\tau/2}\), the bound above
is therefore at most
\[
 \frac{(1-v)\min(Hv,\mu)}
      {2v[1-\min(Hv,\mu)]}.
\]
It increases until \(v=\mu/H\) and decreases thereafter. Its
maximum is \((H-\mu)/(2(1-\mu))<2n^3\), using
\(1-\mu=\sin^2(\pi/(2n))\ge1/n^2\).

There is a quadratic-in-degree bound from the
[circular angular mesh](MULTIVARIATE_GAUSSIAN_LAPLACE_ANGLE_RESPONSE.md#angular-gaps-and-velocities-of-the-reference).
For \(t=n\tau\le7/16\), one has \(\lambda_{n,i}\le2t\le7/8\)
and \(F_{n,i}<2n\), hence the quotient is less than \(8n^2\).
For \(t\ge7/16\) the mesh \(\delta_0\) is at least
\(\sqrt{14/\pi}/n\). If \(i<k\), then
\[
 \lambda_{n,i}-\lambda_{n,k}
 =\sin\frac{x_{n,i}+x_{n,k}}2\sin\frac{x_{n,i}-x_{n,k}}2
 \ \ge\ \frac{(k-i)^2\delta_0^2}{\pi^2}.
\]
In fact the first sine is at least the second, since
\(0<x_{n,k}<x_{n,i}<\pi\). The squared-height force equation gives
\[
 \frac{F_{n,i}}{1-\lambda_{n,i}}
   \le 2+2\sum_{k>i}\frac1{\lambda_{n,i}-\lambda_{n,k}}
   \le 2+\frac{\pi^4}{3\delta_0^2}
   <8n^2.
\]
The last inequality uses \((22/7)^5<315\) and \(n\ge2\).

At a first outer contact, the exact paired force formula gives
\(F_{N,i}-F_{n,i}\ge 1-\lambda\) at the common squared height.
The contact velocity in lower time is
\(\rho(F_{N,i}-F_{n,i})-(1-\rho)F_{n,i}\).
It is strictly positive by (4) and the preceding bound, which
excludes a first contact.

There is also a restriction on any first failure at the sharp ratio
\(\rho=n/N\). Write \(t=n\tau\), \(v=\lambda/(1-\lambda)\) at
contact, and
\[
 \Delta=\frac14\left[N(1-e^{-tn/N})
                     -n(1-e^{-t(n-1)/n})\right]>0 .
\]
The contact measure has mass \(\Delta\) by the trace identity.
Its support is in \([0,1]\), so
\(\int(\lambda-s)^{-2}\,d\nu(s)\ge
\Delta/\max(\lambda^2,(1-\lambda)^2)\).
Substitution in the contact formula, using
\(2n(e^{\tau/2}-1)\ge t\), gives
\[
 v\ \ge\ t\left(1+
    \frac{2\lambda\Delta}{\max(\lambda^2,(1-\lambda)^2)}\right)
 \ \ge\ t\left(1+2\Delta\min\{1,t(1+t)\}\right).                      \tag{5}
\]
For the second inequality use \(v\ge t\) from the first.
The factor \(2\lambda/\max(\lambda^2,(1-\lambda)^2)
=2v(1+v)/\max(1,v^2)\) is increasing for \(v\le1\) and at
least \(2\) for \(v\ge1\).
Since the contact root lies strictly below its lower equilibrium
height, (5) requires its positive-height index \(i\) to satisfy
\[
 i<\frac12+\frac n\pi
   \arctan\!\frac{1}{\sqrt{t(1+2\Delta\min\{1,t(1+t)\})}} .
\]
This confines the first-contact test to outer indices. It leaves the
sign at those indices unresolved.

## An exact degree lift and its limitation

Write \(D=z\partial_z\), \(P=R_n(z;t/n)\), and
\(\delta=t/(nN)\). The same Gaussian pair has the exact lift
\[
 R_N(z;t/N)=A(z)+zA^*(z),\qquad
 A=e^{-\delta D^2}P=R_n(ze^{-t/N};t/N),\qquad
 A^*(z)=z^n A(1/z).                                               \tag{6}
\]
Indeed the \(k\)-th normalized coefficient in \(A\) equals
\(e^{-t k(N-k)/N}\), and the coefficient of \(z^k\) in \(zA^*\)
supplies the complementary term in Pascal's identity.
The Cayley derivative of the upper polynomial from the cut \(1\)
is consequently \(N(A+A^*)\).

This linked lift does require the actual Gaussian trajectory.
For \(n=3\) and \(P=(1+z)^3\), the lift built from
\(A=e^{-\delta D^2}P\) has value
\[
            2-2(3e^{-\delta}+e^{-9\delta})+6e^{-4\delta}
                =-36\delta^2+O(\delta^3)
\]
at \(z=-1\). A real monic quartic with positive coefficients
and all roots on the unit circle has nonnegative value at \(-1\)
by pairing conjugate factors. The displayed value is negative
for sufficiently small positive \(\delta\).
The input can be taken to have simple circular roots by replacing
it by \((1+z)(1+2(1-\delta^3)z+z^2)\), which changes the value
only by \(O(\delta^3)\). Hence (6) cannot give the desired comparison
through circle-root preservation for arbitrary reciprocal circular inputs.

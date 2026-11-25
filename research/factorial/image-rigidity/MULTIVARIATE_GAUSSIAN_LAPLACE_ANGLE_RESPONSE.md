# Degree-uniform angular response to Laplace curvature directions

For \(m\ge2\) and \(c>0\), write
\[
 G(z)=\sum_{k=0}^m a_kz^k,\qquad
 a_k=\binom mk e^{-c k(m-k)/m},\qquad
 \gamma=\min\{\pi,\sqrt{32c/\pi}\},\qquad
 C=\frac12+\frac{\pi^2}{2\gamma}.
\]
Set \(A=\pi/(4c)\) and \(B=\pi^3/(8\gamma^3)\). Let \(\mu\) be a
real signed measure with
\(\int e^{C|h|}(1+h^2)\,d|\mu|(h)<\infty\).
Let \(L\) be the function with \(L(0)=L(1)=0\) and curvature
\[
                  L''(s)=\int_{\mathbb R}e^{h(s-1/2)}\,d\mu(h).
\]
For every root \(z\) of \(G\), the corresponding local root
\(z(\varepsilon)\) of \(\sum_k a_k e^{\varepsilon m L(k/m)}z^k\)
satisfies
\[
 \left|\left.\frac{d}{d\varepsilon}\arg z(\varepsilon)
                 \right|_{\varepsilon=0}\right|
     \le\int_{\mathbb R} e^{C|h|}(A+B h^2)\,d|\mu|(h).           \tag{1}
\]
The bound has no dependence on \(m\). It measures the angular
response to a normalized log-coefficient direction at a Gaussian
reference. The paired factorial application below concerns the
matched difference of two such directions. A nonlinear estimate
along a deformation is still required for the actual factorial
root arguments.

## Angular gaps and velocities of the reference

The [unitary Hermite circle location](https://arxiv.org/html/2203.05533#S2.SS1)
applies to \(G(-z)\) with parameter \(2c/m\). Its roots are simple,
and their signed displacements \(\theta_i\in(-\pi,\pi)\) from
the negative real axis obey
\[
              \theta_i'(c)=\frac1m\sum_{j\ne i}
                          \cot\frac{\theta_i-\theta_j}{2}.
\]
The factor \(1/m\) comes from the Gaussian time \(c/m\). The
[reference heat-flow proof](MULTIVARIATE_GAUSSIAN_REFERENCE_ANGULAR_INTERLACING.md#root-motion-on-the-upper-semicircle)
also gives the simple Hermite limits at \(c=0\), the reflection
symmetry and the outward motion of every positive displacement.

Let \(\delta(c)\) be the least circular angular gap. Rotate a
minimizing pair to \(0,\delta\). Each other point at \(s\) on
the complementary arc contributes
\(\cot(s/2)-\cot((s-\delta)/2)\) to \(m\) times the gap derivative.
Its magnitude is
\[
         \frac{\sin(\delta/2)}
              {\sin((s-\delta)/2)\sin(s/2)} .
\]
On each half of the complementary arc it decreases toward the
midpoint. The first point on either half is at least one gap
from the adjacent end of the arc. The \(j\)-th point is at least
\(j\delta\) from that end. Placing the points at their respective
shortest distances therefore bounds the negative contributions
by telescoping cotangent differences. With \(r\) points on one
half the resulting lower bound on the gap derivative is
\[
 \frac1m\left(\cot\frac{(r+1)\delta}{2}
                      +\cot\frac{(m-r-1)\delta}{2}\right).
\]
For \(m\delta\le\pi\), put \(a=(r+1)\delta/2\) and
\(b=(m-r-1)\delta/2\). Use
\(\cot a+\cot b=\sin(a+b)/(\sin a\sin b)\),
\(\sin(a+b)\ge2(a+b)/\pi\) and
\(\sin a\sin b\le(a+b)^2/4\). The derivative is at least
\(16/(\pi m^2\delta)\). This holds for every minimizing pair, so
\(\delta^2\) has derivative at least \(32/(\pi m^2)\) wherever
\(\delta\le\pi/m\) and the derivative exists. The least gap is
locally Lipschitz for positive times, starts at zero, and cannot
cross down through \(\pi/m\). Hence
\[
                           \delta(c)\ge\gamma/m.              \tag{2}
\]

The positive displacements satisfy a cooperative system
\(x'=F(x)\), since the derivative of each paired force with
respect to the other positive displacement is positive. For
\(\lambda\ge1\) and \(\lambda x_i<\pi\), the function
\[
        f_\lambda(t)=\lambda\cot(\lambda t/2)-\cot(t/2)
\]
is odd and nonincreasing for \(|t|<2\pi/\lambda\), with removable
value zero at \(t=0\). This follows by differentiation from
\(\sin(\lambda t/2)\le\lambda\sin(t/2)\) for
\(0<t<2\pi/\lambda\). Thus the terms at \(x_i-x_j\) and
\(x_i+x_j\) in the scaling difference sum to a nonpositive
number. The self-reflected term at \(2x_i\) and the odd-degree
central term at \(x_i\) are also nonpositive. It follows that
\(\lambda F(\lambda x)\le F(x)\) componentwise. The comparison
principle for cooperative systems therefore gives
\(\lambda^{-1}x(\lambda^2 c)\le x(c)\). To start the comparison
at positive time, compare the former curve with \(x(c+\eta)\)
for fixed \(\eta>0\). Near zero the former tends to zero and all
components of the latter are positive. Compare forward in the
ordered open chamber and then send \(\eta\) to zero.
Differentiating at \(\lambda=1\), and using outward motion, gives
\[
            0\le x_i'(c)\le x_i(c)/(2c)<\pi/(2c).              \tag{3}
\]
By reflection the absolute bound holds for every signed height.

## Exact radial shifts at a root

For real \(h\) put \(t=h/(2m)\) and
\[
 Z(h)=\frac{e^{-h/2}G(ze^{h/m})}{zG'(z)/m}.
\]
For \(w\ne z\) another root put \(q_w=\cot((\arg z-\arg w)/2)\).
Factoring \(G\) gives exactly
\[
         Z(h)=2m\sinh t
                 \prod_{w\ne z}(\cosh t-iq_w\sinh t).
\]
Write \(s_w=|\sin((\arg z-\arg w)/2)|\). The squared modulus
of the \(w\)-factor is \(1+\sinh^2(t)/s_w^2\). By (2) the
\(j\)-th root counterclockwise from \(z\) has
\(s_w\ge\gamma\min(j,m-j)/(\pi m)\). Use
\(|\sinh t|\le |t|e^{|t|}\) in each factor. At most two roots
occur for each \(\min(j,m-j)\). For \(D\ge0\), the decreasing
integrand bound and integration by parts give
\[
 \sum_{j\ge1}\log(1+D^2/j^2)
    \le\int_0^\infty\log(1+D^2/x^2)\,dx=\pi D .
\]
Apply this with \(D=\pi|h|/(2\gamma)\) to the product. The
prefactor and the extracted exponentials together yield
\[
                           |Z(h)|\le |h|e^{C|h|}.
\]

Each factor has positive real part. The sum of their arguments
is \(-\sum_{w\ne z}\arctan(q_w\tanh t)\). The heat equation and
(3) give \(|\sum q_w|\le m\pi/(2c)\). Also
\(\sum |q_w|^3\le 3(\pi m/\gamma)^3\), since
\(\sum_{j\ge1}j^{-3}\le3/2\). Therefore
\(|\arctan v-v|\le|v|^3/3\) and \(|\tanh t|\le|t|\) show that
the sum of arguments has absolute value at most \(A|h|+B|h|^3\).
The prefactor is real, so
\[
                  |\Im Z(h)|\le h^2 e^{C|h|}(A+B h^2).       \tag{4}
\]
The chord-centered second primitive of \(e^{h(s-1/2)}\) for
\(h\ne0\) is \(h^{-2}\) times that exponential minus its affine
endpoint interpolation. For \(h=0\) the primitive is
\((s^2-s)/2\). At a root of \(G\), any real affine coefficient
direction contributes a real multiple of \(zG'(z)\).
Thus implicit differentiation of the root and (4), including
its continuous \(h=0\) limit, give (1).

## Applying the bound to an exact matched factorial direction

Use the midpoint tanh parameters \(r_i,r_i'\) for the same row
factors at orders \(n,N=n+1\), with orientation \(\epsilon_i\).
The [matched tanh curvature formula](MULTIVARIATE_MATCHED_TANH_ORDINARY_ANGLE_LIFT.md#exact-curvature-of-the-assembled-symbol)
represents the \(i\)-th curvature measure as the pushforward under
\(h=-\epsilon_i u\) of
\[
          e^{-u/(2r_i)}D(u/n)\,du,\qquad
          D(v)=v/(2\sinh(v/2)),\qquad u\ge0.
\]
The upper measure uses \(r_i'\) and \(N\). Let
\(r_*\ge\max_i\{r_i,r_i'\}\), \(d=1-2Cr_*>0\) and
\(V_j=\sum_i r_i^j\). The weighted variation of their matched
difference is at most
\[
 \left(\frac{2A}{d^2}+\frac{48B r_*^2}{d^4}\right)
                   \frac{V_2}{nN-r_*}
 +\left(\frac1{n^2}-\frac1{N^2}\right)
               \left(\frac{2A}{3d^3}+\frac{32 B r_*^2}{d^5}\right)V_3.
                                                               \tag{5}
\]
Indeed, the matched identity bounds
\(|r_i'-r_i|\le r_i^2/(nN-r_*)\).
At fixed order the parameter derivative of the measure density
has magnitude at most \(u e^{-u/(2r)}/(2r^2)\), since \(D\le1\).
After multiplying by \(e^{Cu}(A+Bu^2)\), integrate using
\(1/(2r)-C\ge d/(2r)\). At fixed \(r_i\), the order change is
bounded by \(e^{-u/(2r_i)}(1/n^2-1/N^2)u^2/24\), using the
matched curvature slope bound for \(D\). The same integration
gives the second term of (5).

Equations (1) and (5) bound the angular derivative of the exact
difference of fractional-grid log-coefficient corrections at
any fixed Gaussian reference degree \(m\). Comparing the two
different reference root configurations and controlling the full
deformation remain necessary for a general separation theorem.

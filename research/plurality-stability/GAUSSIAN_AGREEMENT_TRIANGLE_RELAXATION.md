# Angular exactness and a mixed-radius obstruction for Gaussian agreement kernels

For a balanced planar three-cell partition \(A\), its agreement
kernel is \(F_A(x,y)=\sum_{i=1}^3 1_{A_i}(x)1_{A_i}(y)\).
The range, balance and triangle relaxation consists of symmetric
measurable kernels satisfying
\[
  0\le F\le1,\quad
  \int F(x,y)\,d\gamma_2(y)=\frac13,\quad
  F(x,y)+F(y,z)-F(x,z)\le1 .                              \tag{AT1}
\]
The row condition and the inequalities can be imposed almost
everywhere, in \(\gamma_2\), \(\gamma_2^2\) and \(\gamma_2^3\)
as appropriate. The Gaussian objective is
\(\mathcal S_\rho(F)=\mathbb E F(X,Y)\) for standard \(X,Y\)
with cross covariance \(\rho I_2\). Averaging a feasible
kernel over simultaneous rotations and reflections preserves
(AT1) and the objective. This is a convex averaging argument
and does not assume that \(A\) itself has those symmetries.
The graph partition analogue uses the same triangle
constraints and symmetry aggregation.
[van Dam and Sotirov, §4 and §4.1](https://arxiv.org/pdf/1312.0332)

The angular certificate below is exact. The full Gaussian
relaxation fails for every \(\rho\in[9999/10000,1)\), even
with all mixed-radius triangles and the pointwise condition
\(F(x,x)=1\). The counterexample is a fractional agreement
kernel. In particular it violates positive semidefiniteness
and the four-point inequality for three labels.

## The angular inequality and an exact positive certificate

Let \(t\in[0,\pi]\) be circular geodesic distance. Set
\(a=2\pi/3\) and \(F_T(t)=(1-t/a)_+\). This is the simultaneous
rotation average of the agreement kernel of a regular tripod.
For an even angular kernel \(F\) satisfying (AT1) on the circle,
write \(d=1-F\) and \(H(t)=\int_0^t d(u)\,du\).
The triangle inequality on an arc with increments \(u,t-u\)
gives, on integration, \(t d(t)\le2H(t)\).
Consequently \(H(t)/t^2\) is nonincreasing. Also
\(H(\pi)=2\pi/3\), and \(d\le1\) gives \(H(t)\ge t-a/2\).
In particular \(H(a)\ge a/2\). For \(t\le a\) the ratio
inequality now gives \(H(t)\ge t^2/(2a)\). These are precisely
\(H(t)\ge H_T(t)\) on both parts of the interval. Integration
by parts against any decreasing angular weight \(k\) proves
\[
                   \int_0^\pi k(t)(F_T(t)-F(t))\,dt\ge0.       \tag{AT2}
\]
The same proof uses the Stieltjes derivative for bounded
monotone weights. For a smooth weight with \(k'<0\) on
\((0,\pi)\), equality forces \(H=H_T\) and \(F=F_T\) almost
everywhere. In particular this applies to
\(k(t)=e^{\kappa\cos t}\), \(\kappa>0\).

For \(C^1\) weights write
\[
 c_k=\frac2{a^2}\int_0^a t k(t)\,dt,\qquad
 w_k(s)=\frac1{s^3}\int_0^s t^2[-k'(t)]\,dt.
\]
Then
\[
\begin{split}
 \int_0^\pi k(t)(F_T-F)(t)\,dt
   ={}&\int_0^a w_k(s)\int_0^s
      [1-F(u)-F(s-u)+F(s)]\,du\,ds\\
     &+\int_a^\pi[c_k-k(t)]F(t)\,dt .                         \tag{AT3}
\end{split}
\]
To check the identity, put \(G=F_T-F\). The inner triangle
integral is \(2\int_0^s G-sG(s)\). After interchanging the
integrals, its coefficient at \(t<a\) is
\(2\int_t^a w_k(s)\,ds-tw_k(t)\). Its derivative is \(k'(t)\)
and its value at \(a\) is \(k(a)-c_k\). The remaining term
has coefficient \(k-c_k\) on \(G=-F\) for \(t>a\). Finally
\(\int_0^\pi G=0\). Both terms on the right of (AT3) are
nonnegative for decreasing \(k\).

For a radial-angular feasible kernel, the angular average
\(\int F(r,s,t)\,d\mu(r)d\mu(s)\) is feasible for this circle
problem, where \(d\mu(r)=r e^{-r^2/2}dr\). Indeed one can
integrate each triangle over three independently sampled
radii. This proves (AT2) for angular objectives whose radial
factor is independent sampling from \(\mu\) at both endpoints.
The Mehler interaction has a different radial law.

## A balanced mixed-radius counterexample

Let \(R=\sqrt{2\log2}\) be the median Gaussian radius, and let
\(B(x,y)\) indicate that \(|x|,|y|\) lie in different halves
of this median cut. Let \(t=t(x,y)\) be the geodesic angle
between \(x,y\), using any fixed angle convention at the
origin. For \(0\le\epsilon\le1/1000\) set
\[
 a_\epsilon=\frac{a}{1-\epsilon+\epsilon^2/2},
 \qquad
 F_\epsilon(x,y)=\left(1-\frac{t(x,y)}{a_\epsilon}
                              -\epsilon B(x,y)\right)_+ .    \tag{AT4}
\]
Then \(1-F_\epsilon\) is the truncation at one of the sum of
the angular pseudometric \(t/a_\epsilon\) and the cut
pseudometric \(\epsilon B\). Truncating a nonnegative
pseudometric at one preserves its triangle inequality.
Thus \(F_\epsilon\) satisfies the triangle and range
constraints and has pointwise diagonal one.

For fixed \(x\), the radius and angle of \(Y\sim\gamma_2\)
are independent, \(B(x,Y)\) takes values zero and one with
equal probability, and \(t(x,Y)\) is uniform on \([0,\pi]\).
Since \(a_\epsilon<\pi\), the row mean is exactly
\[
       \frac{a_\epsilon}{4\pi}[1+(1-\epsilon)^2]=\frac13.
\]
The angular slope reduction in (AT4) is balanced by a
penalty on pairs across a radial cut. At high correlation
the pairs across the cut are rare enough to give a gain.

### A strict gain on an explicit correlation band

Set \(\alpha=\arccos\rho\), \(\sigma^2=1-\rho^2\) and
\(\delta=1-\rho\), and take \(\epsilon=1/1000\).
On \(t<a\) the inequality \((z)_+\ge z\) gives
\[
 F_\epsilon-F_T
   \ge \epsilon\left[\left(1-\frac\epsilon2\right)\frac{t}{a}-B\right].
\]
On \(t\ge a\), \(F_T=0\) so the contribution to the objective
is nonnegative. Consequently
\[
 \frac{\mathcal S_\rho(F_\epsilon)-S_\rho(T)}{\epsilon}
 \ \ge\
   \left(1-\frac\epsilon2\right)\frac{\mathbb E t}{a}
        -\mathbb P(B=1)-\frac32\mathbb P(t\ge a).           \tag{AT5}
\]
All expectations here use the true correlated Gaussian law.

There are three elementary bounds for (AT5).

- \(\mathbb E t=\alpha\). For each \(x,y\), the fraction
  of uniformly chosen linear directions whose two signs
  disagree is \(t(x,y)/\pi\). Under the Gaussian pair law
  the sign disagreement for any fixed direction is
  \(\alpha/\pi\). For the latter identity, represent the
  two scalar observations by \(G\) and
  \(G\cos\alpha+H\sin\alpha\) with independent standard
  \(G,H\), and use the uniform angle of \((G,H)\).
- \(\mathbb P(B=1)\le\sqrt{\log2/\pi}\,\alpha\).
  Realize the pair as \(X=X_0\), \(Y=X_\alpha\) with
  \(X_u=X\cos u+Z\sin u\) and \(Z\) independent standard.
  At each \(u\), \(X_u\) and its derivative are independent
  standard Gaussians. For a bounded smooth radial
  approximation \(f\) of the median disk indicator,
  \(\mathbb E|f(X_\alpha)-f(X_0)|\le
  \alpha\sqrt{2/\pi}\,\mathbb E|\nabla f(X)|\)
  by integrating along this path.
  A radial transition layer of width tending to zero
  has \(\mathbb E|\nabla f(X)|\to R e^{-R^2/2}=R/2\).
- \(\mathbb P(t\ge a)\le(1-\rho^2)/(4\rho^2)\).
  Since \(a>\pi/2\), condition on \(r=|X|\) and align
  the first axis with \(X\). The event implies
  \(\rho r+\sigma Z_1\le0\). Integrating \(r\) against
  its Rayleigh law and using \(1-e^{-x}\le x\) bounds
  the probability by
  \(\sigma^2\mathbb E[(Z_1^-)^2]/(2\rho^2)=\sigma^2/(4\rho^2)\).

For the remaining arithmetic, \(3/(2\pi)-\sqrt{\log2/\pi}
>7/1000\). For example \(\pi<22/7\), \(\pi>3141/1000\),
and \(\log2<347/500\) give
\(3/(2\pi)>21/44\) and \(\sqrt{\log2/\pi}<4701/10000\).
The log bound follows by bounding the positive tail
of \(2\sum_{j\ge0}3^{-(2j+1)}/(2j+1)\).
Since \(1/a<1/2\),
\[
 \left(1-\frac\epsilon2\right)\frac1a-\sqrt{\frac{\log2}\pi}
       >\frac{27}{4000}.
\]
For \(\rho\ge9999/10000\), the last term of (AT5) is
at most \((19/25)\delta\), and
\(\alpha\ge\sqrt{2\delta}>(7/5)\sqrt{\delta}\).
Substituting and using \(\sqrt{\delta}\le1/100\) gives
\[
        \mathcal S_\rho(F_{1/1000})-S_\rho(T)
                      >10^{-6}\sqrt{1-\rho} >0 .              \tag{AT6}
\]
Thus an exact Gaussian dual certificate using only (AT1)
cannot give the simplex bound across all positive
correlations. The proof makes no quadrature assumption.

## Two necessary conditions that reject this kernel

An ordinary agreement kernel is positive semidefinite
as an integral operator on \(L^2(\gamma_2)\). This also
holds after common rotation averaging. For (AT4), the
angular Fourier coefficient at order \(m\ge1\) for cut
value \(b\in\{0,1\}\) is
\[
         \frac{1-\cos(m a_\epsilon(1-\epsilon b))}
                       {\pi a_\epsilon m^2}.
\]
At \(m=2\), both cosine arguments are in \((\pi,2\pi)\)
for the range of \(\epsilon\) used here. The coefficient
is larger across the radial cut than within a half for
every \(\epsilon>0\). Hence the quadratic form is
negative on the radial half contrast times \(\cos(2\theta)\).
This identifies a single mode that rejects the witness.

A three-label agreement kernel also satisfies
\[
                 \sum_{1\le i<j\le4} F(x_i,x_j)\ge1        \tag{AT7}
\]
because four points include a pair with the same label.
Take four equally spaced angles and alternate the two
radial halves. The two opposite pairs have \(F_\epsilon=0\).
Each adjacent pair has value
\(1/4-\epsilon/4-3\epsilon^2/8>0\). Their sum is
\(1-\epsilon-3\epsilon^2/2<1\). Continuity in the angles
and positive masses of both radial halves show failure
on a positive-measure set of four-point configurations.
Thus either positive semidefiniteness or (AT7) must
exclude this near-one gain. Exactness of a relaxation
with stronger constraints requires a separate argument.
[Symmetry and four-point constraints in graph partition relaxations, §4](https://arxiv.org/pdf/1312.0332)
[Exact finite-cell limits for the Gaussian relaxation](GAUSSIAN_AGREEMENT_TRIANGLE_DISCRETIZATION.md)

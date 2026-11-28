# A parabolic Gaussian score restriction

Fix \(0<\rho<1\) and a nonzero real number \(\kappa\).
For a bounded measurable function \(g\) on the line, set
\[
                    f(x,y)=g(x-\kappa y^2).
\]
If \(T_\rho f\) is constant on an open arc of any parabola
\(x=c+\kappa y^2\), then \(g\) is constant almost everywhere.
In particular, no ordinary stationary Gaussian winner
partition of the plane with two or more positive-mass
labels can have its labels determined by \(x-\kappa y^2\).
This includes invariance under the corresponding parabolic
shear flow, in coordinates with its axis centered on the
Gaussian origin. No maximality assumption is needed.

The proof uses a scalar characteristic function followed
by uniqueness of the Cauchy transform on the two sides of
a smooth curve. The curve is determined by the conditional
Gaussian law of \(x-\kappa y^2\).

## The characteristic curve

Write \(\sigma^2=1-\rho^2\). At observation
\((c+\kappa y^2,y)\) the conditional value of the source
coordinate \(x-\kappa y^2\) has the law of
\[
 X_y=\rho c+\rho\kappa y^2+\sigma Z_1
                         -\kappa(\rho y+\sigma Z_2)^2,
\]
where \(Z_1,Z_2\) are independent standard Gaussian
variables. Direct one-dimensional Gaussian integration
gives, for real \(t\),
\[
 \mathbb E e^{itX_y}=\Psi(t)e^{y^2\Lambda(t)},\qquad
 \Psi(t)=\frac{e^{i\rho ct-\sigma^2 t^2/2}}
                     {(1+2i\kappa\sigma^2 t)^{1/2}},\qquad
 \Lambda(t)=i\kappa\rho t-
                  \frac{i\kappa\rho^2t}{1+2i\kappa\sigma^2t}.
                                                               \tag{PA1}
\]
The square root is taken continuously from \(t=0\).
In particular \(\Psi\) never vanishes on the real line,
and it and its derivatives have Gaussian decay times
polynomial bounds.

Let \(b_0=2\kappa\sigma^2\). The real part of \(\Lambda\)
lies in \([-\rho^2/(2\sigma^2),0]\), and
\[
       \operatorname{Im}\Lambda(t)
             =\kappa\rho t\left(1-\frac{\rho}{1+b_0^2t^2}\right).
\]
Its derivative has the sign of \(\kappa\) and absolute
value at least \(|\kappa|\rho(1-\rho)\). Indeed the
derivative of \(t/(1+b_0^2t^2)\) is at most one.
Thus \(\Gamma=\Lambda(\mathbb R)\) is a smooth proper
graph \(\{\operatorname{Re}z=\psi(\operatorname{Im}z)\}\).
Its complement has two connected components, one to
each side of the graph. The rational formula makes
all derivatives of \(\Lambda\) polynomially bounded.

## Uniqueness on the graph

The Gaussian score is real analytic, so equality on a
parabola arc extends to the whole parabola. Subtract the
constant value from \(g\), obtaining \(g_0\). With
the Fourier convention \(\widehat g_0(t)=\int
e^{-itx}g_0(x)\,dx\) in tempered distributions, (PA1)
gives
\[
              J(w):=\langle\widehat g_0,
                             \Psi(t)e^{w\Lambda(t)}\rangle=0
                                          \quad(w\ge0).        \tag{PA2}
\]
For real \(w=y^2\) this follows from Fourier inversion
of the density of \(X_y\). Equivalently, the characteristic
function in (PA1) is a Schwartz test, so the pairing
gives \(2\pi\mathbb E[g_0(X_y)]\) under this convention.
The same test belongs to the Schwartz class for every
complex \(w\). On any compact set of \(w\)'s all
\(t\)-derivatives of the exponential are bounded by
a polynomial times \(e^{C(1+|t|)}\), which is absorbed
by the Gaussian factor in \(\Psi\). The same bounds
with powers of \(\Lambda\) justify holomorphic
differentiation in the Schwartz topology. Hence \(J\)
is entire and (PA2) holds for every complex \(w\).

For \(z\notin\Gamma\) define
\[
        C(z)=\left\langle\widehat g_0,
                        \frac{\Psi(t)}{z-\Lambda(t)}\right\rangle .
\]
This function is holomorphic on the complement of the
graph by local denominator and derivative bounds.
If \(\operatorname{Re}z>0\), write the reciprocal as
\(\int_0^\infty e^{-u(z-\Lambda(t))}\,du\) to obtain
\(C(z)=\int_0^\infty e^{-uz}J(u)\,du=0\).
The interchange is justified in the Schwartz topology
by \(\operatorname{Re}\Lambda\le0\) and exponential
decay in \(u\) against any fixed polynomial in \(u\).
For \(\operatorname{Re}z<-\rho^2/(2\sigma^2)\), use
\(-\int_0^\infty e^{u(z-\Lambda(t))}\,du\) and \(J(-u)=0\).
Holomorphic continuation gives \(C=0\) on both
components of the complement.

For completeness, define the curve-supported distribution
\[
 \mu(\varphi)=\langle\widehat g_0,
                         \Psi(t)\varphi(\Lambda(t))\rangle
 \qquad(\varphi\in C_c^\infty(\mathbb C)).
\]
It acts only on the restriction of a test to the curve.
Let \(U=(1/(\pi z))*\mu\) as a distribution. This
convolution is well defined here: for a compact smooth
test \(\varphi\), pair \(\widehat g_0\) against \(\Psi(t)\)
times the smooth function
\(\int \varphi(z)/[\pi(z-\Lambda(t))]\,d^2z\).
The kernel is locally integrable and derivatives can
be moved onto the test; the resulting functions have
polynomial bounds before the Gaussian factor. The
distributional identity
\(\partial_{\bar z}(1/(\pi z))=\delta_0\), obtained by
integration by parts off a disk about zero, gives
\(\partial_{\bar z}U=\mu\). Off \(\Gamma\) the distribution
\(U\) equals the holomorphic function \(C/\pi\), so \(U\)
is supported on \(\Gamma\).

In the local coordinates
\(X=\operatorname{Re}z-\psi(\operatorname{Im}z)\),
\(Y=\operatorname{Im}z\), a distribution supported on
the curve admits locally a finite expansion
\(\sum_{j=0}^N u_j(Y)\delta^{(j)}(X)\). This follows
by finite local order and Taylor expansion of a test
in the normal variable: a test with zero normal jets
through that order can be cut off in arbitrarily
small neighborhoods of \(X=0\), where its derivatives
through that order tend uniformly to zero. Thus the
distribution depends only on the finitely many jets,
which give the coefficients \(u_j\). In these same
coordinates
\[
       \partial_{\bar z}
          =\tfrac12\big[(1-i\psi'(Y))\partial_X
                                      +i\partial_Y\big].
\]
If \(U\) has a nonzero highest normal coefficient
\(u_N\), its derivative has nonzero normal coefficient
\((1-i\psi')u_N/2\) at order \(N+1\). But \(\mu\)
acts only on the order-zero jet. Therefore every
coefficient of \(U\) vanishes locally. It follows that
\(\mu=0\). Since \(\Lambda\) is a smooth
embedding, every compactly supported smooth one-variable
test is the pullback of a smooth compactly supported test
in the plane. Consequently the product distribution
\(\Psi\widehat g_0\) vanishes. In the case of (PA1) the
factor \(\Psi\) never vanishes, so local division gives
\(\widehat g_0=0\) and \(g=b\) almost everywhere.

## Consequence for stationary partitions

Suppose the labels of a stationary Gaussian winner
partition are measurable functions of \(u=x-\kappa y^2\),
and all of its used labels have positive Gaussian mass.
Their adjusted Gaussian pair contrasts are nonzero
analytic functions by injectivity of the Gaussian
noise operator. On a generic horizontal line they
restrict to nonzero analytic functions of \(x\), so
their pair equality sets there are discrete. To see
this simultaneously for all pairs, for each contrast
choose some fixed \(x\) at which it is not identically
zero as a function of \(y\). The exceptional horizontal
heights are contained in the union of the discrete
sets of zeros at those choices of \(x\).
Fubini allows choosing such a height at which the
ordinary labels also obey the score rule almost
everywhere.

Along this horizontal line the variable \(u\) differs
from \(x\) by a constant. Every used label occurs on
a set of positive one-dimensional measure, since
the Gaussian law of \(u\) has a positive density
on the whole line. Hence a genuine adjacent label
switch occurs at some \(u=c\). The one-sided label
intervals at that switch correspond to the adjacent
open sides of the whole parabola
\(x=c+\kappa y^2\). The two adjusted scores therefore
tie along it by continuity. Their unadjusted contrast
is \(T_\rho g(x-\kappa y^2)\) for a nonconstant
bounded label contrast \(g\), contrary to the theorem.


## An additional transverse polynomial

The same conclusion holds in higher dimensions for
\(f(x,y,z)=g(x-\kappa y^2-P(z))\), where \(P\) is a real
polynomial on the remaining transverse coordinates.
A constant score on an open part of the graph
\(x=c+\kappa y^2+P(z)\) forces \(g\) to be constant.
This includes every nonflat quadratic graph with
no linear term in the distinguished quadratic direction,
after diagonalization in the transverse block.

Indeed the score equality extends analytically to the
whole graph. At each fixed observation \(z\), the
characteristic formula (PA1) is multiplied by
\[
          e^{it\rho P(z)}
                \big[T_\rho^z(e^{-it P})\big](z).
\]
The \(\Lambda(t)\) factor does not change because the
Gaussian source coordinates \(y,z\) are independent.
For fixed \(z\) the new multiplier is bounded and
has bounded derivatives of any fixed order in \(t\),
by the Gaussian moments of the polynomial. The
same Cauchy argument shows that the product of the
new multiplier and \(\Psi\widehat g_0\) is zero
for every observation \(z\). Given a real frequency
\(t_0\), choose \(z\) at which
\(T_\rho^z(e^{-it_0P})(z)\ne0\). Such an observation
exists by Gaussian convolution injectivity, since
the source function has modulus one. The multiplier
is nonzero in a neighborhood of \(t_0\) at this fixed
\(z\), so local division there makes \(\widehat g_0\)
zero. This holds at every frequency, proving the
assertion. The extension to ordinary stationary
partitions uses the same horizontal-line argument.

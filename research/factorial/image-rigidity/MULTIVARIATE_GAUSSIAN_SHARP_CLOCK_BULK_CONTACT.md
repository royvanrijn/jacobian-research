# Edge confinement for first contacts at the sharp Gaussian clocks

For \(t>0\), let
\[
 Q_j(z;t)=\sum_{k=0}^j\binom jk e^{-t k(1-k/j)}z^k .
\]
The physical lower clock ratio for consecutive degrees \(n\) and
\(N=n+1\) uses the same normalized time \(t\) here, with raw times
\(t/n\) and \(t/N\). The [Gaussian heat theorem](MULTIVARIATE_GAUSSIAN_REFERENCE_ANGULAR_INTERLACING.md)
gives the strictly separated inner inequalities at these times.
Strict outer separation starts at small time by the scaled Hermite
comparison. Put \(\lambda=\sin^2(x/2)\) for an angular height \(x\)
from \(-1\) and order the positive squared heights decreasingly.

**Theorem.** Let \(n\to\infty\) along a sequence of first outer
contacts at times \(t_n\to t_0>0\), with contact squared heights
\(\lambda_n\to\lambda_0\) and positive-height indices \(i_n\).
Then
\[
 \lambda_0=\Lambda(t_0),\qquad \frac{i_n}{n}\ \longrightarrow\ 0 ,
                                                                  \tag{1}
\]
where
\[
 \Lambda(t)=
 \begin{cases}
 \sin^2\!\left(\arcsin\sqrt{t/2}+\sqrt{t(2-t)}/2\right),&0<t<2,\\
 1,&t\ge2 .
 \end{cases}
\]
In particular, for any compact positive normalized-time interval
and any \(\eta>0\), all sufficiently large degrees exclude a first
outer contact of index \(i\ge\eta n\) in that interval. At every
subcritical time, every finite-degree squared height is strictly
below \(\Lambda(t)\), so a first contact can only approach the
edge from within the limiting support.

## Adjacent normalized moment differences

Rotate the roots by \(w=-z\), so that the squared height is
\(s=(2-w-w^{-1})/4\). Define
\[
 u_{j,r}(t)=j^{-1}\sum_w w^r,\qquad
 H_{j,r}(t)=j^{-1}\sum_w \big((2-w-w^{-1})/4\big)^r .
\]
The sums include all \(j\) roots. For \(r\ge1\), \(H_{j,r}\) equals
\(2/j\) times the \(r\)-th power sum of the positive squared heights.
Indeed, the additional root in odd degree has squared height zero.

For \(D=w\partial_w\), the monic rotated polynomial satisfies
\(\partial_t Q=(j^{-1}D^2-D)Q\). Expand its logarithm at infinity
with nonpolynomial part
\(-j\sum_{r\ge1} u_{j,r}w^{-r}/r\).
Comparing coefficients gives the triangular system
\[
 \dot u_{j,r}
   =-r(1-r/j)u_{j,r}-r\sum_{k=1}^{r-1}u_{j,k}u_{j,r-k},
 \qquad u_{j,r}(0)=1.                                              \tag{2}
\]
For each fixed \(r\), this equation defines an analytic extension in
\(h=1/j\) near zero, with derivatives locally bounded for \(t\) in
a compact interval. This follows inductively from its integrating
factor since the sum has only smaller indices.

On the circle, reciprocity gives \(u_{j,-r}=u_{j,r}\). Expanding
\((2-w-w^{-1})^r=(-1)^r w^{-r}(w-1)^{2r}\) now writes \(H_{j,r}\)
as a finite linear combination of the \(u_{j,k}\) and a constant.
It has the same extension \(H_r(h,t)\).
Write \(H_r(h,t)=H_r(0,t)+h G_r(h,t)\) with \(G_r\) analytic
near zero. Then for each fixed \(r\ge1\),
\[
 N H_{N,r}(t)-n H_{n,r}(t)
       =H_r(0,t)+O_r(n^{-2})                                      \tag{3}
\]
locally uniformly in \(t\). The remainder is
\(G_r(1/N,t)-G_r(1/n,t)\), bounded by a fixed derivative times
\(1/(nN)\).

Let \(\mu_t\) be the limiting probability law of squared heights.
Its moments are \(H_r(0,t)\). Existence follows by compactness
on \([0,1]\), and uniqueness follows by polynomial approximation.
The same argument shows weak convergence along varying \(t_n\to t\),
using the uniformity of the analytic moment extension.
For a fixed \(t\), [Kabluchko, Theorems 2.2 and 2.18](https://arxiv.org/abs/2203.05533)
identify the rotated limiting root measure as the free unitary normal
law with parameter \(2t\). Its squared-height image has support
\([0,\Lambda(t)]\) and has no atoms. In the subcritical case
\(t<2\), its density vanishes as a square root at \(\Lambda(t)\).
The mapping from angle to squared height is strictly increasing
on the positive half-circle and has nonzero derivative at the
subcritical edge.

## A zero-free sector at every subcritical degree

For \(0<t<2\) define
\[
 b_*=\arccos\sqrt{t/2},\qquad
 \theta_* = 2b_*-t\tan b_* = \pi-2\arcsin\sqrt{t/2}-\sqrt{t(2-t)} .
\]
Then \(Q_j(z;t)\ne0\) for \(z=e^{i\theta}\) and
\(|\theta|\le\theta_*\), for every \(j\ge1\).
The centered integral has a positive constant \(C_{j,t}\) independent
of \(\theta\) and satisfies
\[
 e^{-ij\theta/2}Q_j(e^{i\theta};t)
   =C_{j,t}\int_{\mathbb R}
        e^{-j x^2/(4t)}\cosh^j\!\left(\frac{x+i\theta}{2}\right)\,dx.
                                                                  \tag{4}
\]
This follows by expanding the hyperbolic cosine power and integrating
each Gaussian term.

Take \(0<\theta<\theta_*\). The function \(2b-t\tan b\) increases
from zero to \(\theta_*\) on \([0,b_*]\), so choose the unique
\(b_0<b_*\) with \(2b_0-t\tan b_0=\theta\), and put
\(Y=t\tan b_0\). For \(x>0\) choose \(y(x)\in[0,Y]\) by
\[
 y(x)=\frac{2t}{x}
     \arctan\!\left(\tanh(x/2)\tan\frac{\theta+y(x)}2\right).        \tag{5}
\]
Its right side minus \(y\) is positive at \(y=0\) and negative
at \(Y\), since \(2\tanh(x/2)<x\) and \(\arctan v\le v\) for
positive \(v\). For all \(0\le y\le Y\) its derivative is at most
\((t/2)\sec^2 b_0-1<0\), by the same bound on \(\tanh(x/2)\)
and the derivative of \(\arctan\).
At \(x=0\) the continuous equation is \(y=t\tan((\theta+y)/2)\),
with solution \(Y\). By the implicit-function theorem \(y(x)\)
extends smoothly and evenly to real \(x\).

At \(z=x+i y(x)\), the one-factor integrand
\(e^{-z^2/(4t)}\cosh((z+i\theta)/2)\) is strictly positive real.
For \(x>0\), the argument of the hyperbolic cosine is the arctangent
in (5), and the phase of the exponential is \(-xy(x)/(2t)\).
At \(x=0\) the hyperbolic cosine is positive since
\(\theta+y(x)\le 2b_0<\pi\). Reflection gives the same positive
integrand at negative \(x\). Deform the line in (4) to the graph
of \(y(x)\). The integrand is entire and the joining vertical
segments vanish in the limit by Gaussian decay, since the contour
height is bounded. The imaginary part of its Jacobian,
\(i y'(x)\), is odd and cancels on symmetric segments. The deformed
integral is therefore the integral of a strictly positive real
function against \(dx\).

At \(\theta=\theta_*\), for each \(x>0\) the equation (5) still has
a strict derivative bound in \(y\), because \(2\tanh(x/2)/x<1\).
The solution is the limit of the solutions for \(\theta<\theta_*\).
All heights are bounded by \(t\tan b_*\), so the Gaussian factor
times \(\exp(j|x|/2)\) gives an integrable majorant independent of
\(\theta\). Dominated convergence gives the same positive integral
at the edge. The value at zero is positive by the real integral,
and conjugation covers negative angles.
Finally \(\cos^2(\theta_*/2)=\Lambda(t)\), which proves strict
finite-degree confinement inside the limiting squared-height edge.

## The limiting first-contact measure

At a first outer contact, pair the ordered positive squared heights
\(x_k\) and \(y_k\) in degrees \(n\) and \(N\). Append a zero to the
lower list if \(n\) is odd. The contact has \(x_i=y_i=\lambda_n\)
and \(y_k\ge x_k\) for every pair. The positive measure
\[
  d\nu_n(s)=\sum_{k\ne i}{\bf1}_{[x_k,y_k]}(s)\,ds
\]
has moment of order \(r\ge0\) equal to
\([N H_{N,r+1}(t_n)-n H_{n,r+1}(t_n)]/[2(r+1)]\).
Therefore (3) gives the weak limit
\[
                    d\nu_\infty(s)=\tfrac12\mu_{t_0}([s,1])\,ds.     \tag{6}
\]
The identification follows from Fubini, and moment convergence gives
weak convergence because the measures have bounded mass on \([0,1]\).

The [exact contact test](MULTIVARIATE_GAUSSIAN_PHYSICAL_CLOCK_COMPARISON.md#a-root-moment-test-for-the-outer-gaps-at-the-physical-clock-ratio)
requires, at a first contact,
\[
 1+2\lambda_n\int\frac{d\nu_n(s)}{(\lambda_n-s)^2}
     \ \le\ \frac{F_{n,i}}{n(1-\lambda_n)}
     \ \le\ \frac{\lambda_n}{t_n(1-\lambda_n)}.                    \tag{7}
\]
The upper bound uses the raw speed estimate and
\(2n(e^{t_n/(2n)}-1)\ge t_n\).
In particular \(\lambda_0\ge t_0/(1+t_0)>0\).

If \(\lambda_0<\Lambda(t_0)\), equation (6) gives a density bounded
below by a positive constant throughout a neighborhood of
\(\lambda_0\). Test weak convergence with
\([(\lambda_n-s)^2+\varepsilon^2]^{-1}\), which tends uniformly to
the corresponding kernel centered at \(\lambda_0\) for each fixed
\(\varepsilon>0\). Then let \(\varepsilon\) decrease to zero.
The integral in (7) must diverge, whereas its right side stays bounded.
This excludes a bulk first contact. The finite-degree confinement
and the continuity of \(\Lambda\) give
\(\lambda_0\le\Lambda(t_0)\), completing the height claim.

At contact there are at least \(i_n\) positive lower heights at
or above \(\lambda_n\). If \(i_n/n\) were bounded below along a
subsequence, weak convergence would put positive mass at or above
\(\lambda_0\) in \(\mu_{t_0}\). More explicitly, use the closed
intervals \([\lambda_0-\varepsilon,1]\) and let \(\varepsilon\)
decrease to zero. This contradicts the height claim and absence of
an atom at the edge. The compact-time statement follows by passing
to a convergent subsequence of times and contact heights.

## The leading edge identity

For \(0<t<2\), set \(L=\Lambda(t)\) and let \(\nu_\infty\) be given
by (6). Its edge resolvent is finite. Direct integration gives
\[
 (1-L)\left[1+2L\int\frac{d\nu_\infty(s)}{(L-s)^2}\right]
      =L(1-L)\int\frac{d\mu_t(u)}{L-u}.                            \tag{8}
\]
The right side is the continuous force in the contact equation,
obtained by replacing the discrete normalized force sum by its
limiting law. Equation (8) is an identity of the limiting
expressions. Convergence of the discrete force at the edge would
require additional near-root estimates, since weak moment
convergence alone does not control its singular kernel. The
identity explains why the leading bulk-resolvent comparison
leaves the sharp edge sign open.

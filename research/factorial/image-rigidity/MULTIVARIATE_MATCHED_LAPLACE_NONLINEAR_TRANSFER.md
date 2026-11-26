# Nonlinear angular transport of a matched factorial symbol

The change of the fractional-grid log-factorial symbol from one
order to the next is a positive Laplace curvature direction.
Grouping by the original shared rows gives the positivity exactly.
As a result, the finite change of symbol on a fixed grid acts by
positive averaging of radial shifts, up to a scalar and a positive
rescaling. The theorem below bounds the resulting root movement
for every degree with a lower bound on the angular mesh.

## Positive curvature increment for ordinary shared rows

Let \(L_j(s)\) be the log-factorial multiplier per degree at order
\(j\), with its endpoints subtracted by affine interpolation and the
binomial coefficients removed. For an increasing row
\((\ell,\ell+d)\), its contribution to the uncentered symbol is
\(j^{-1}\log\Gamma(j(\ell+ds)+1)\). A decreasing row uses \(1-s\).
For each row write \(\epsilon=1\) in the increasing orientation and
\(-1\) in the decreasing orientation, and \(a=\ell/d+1/2\).
The standard trigamma integral, applied after two differentiations
and the substitution \(u=jd\,t\), yields the curvature density
\[
 L_{j,\mathrm{row}}''(s)
    =\int_0^\infty e^{-\epsilon u(s-1/2)}
           d e^{-a u} f\!\left(\frac{u}{dj}\right)\,du,
       \qquad f(t)=\frac{t}{e^t-1}.
\]
This is also the sum of the Gauss-factor densities in the
[matched tanh representation](MULTIVARIATE_MATCHED_TANH_ORDINARY_ANGLE_LIFT.md#exact-curvature-of-the-assembled-symbol).
For \(t\ge0\), \(0\le-f'(t)\le1/2\). The lower inequality follows
from \(e^t(1-t)\le1\). The upper inequality reduces to
\(\sinh t\ge t\). The value at zero is interpreted by continuity.

Put \(N=n+1\). The curvature of \(L_N-L_n\) is the Laplace transform
of a positive measure \(\mu\), obtained by adding the rowwise
increments at the signed coordinate \(h=-\epsilon u\).
The increment density of each such row in \(u\) lies between zero
and \(u e^{-a u}/(2nN)\). In particular \(\nu(dh)=h^{-2}\mu(dh)\)
has a finite weighted first moment whenever the chosen exponential
weight has rate less than every row's \(a\). Up to an affine
function of \(s\), the symbol difference is therefore
\[
       L_N(s)-L_n(s)=\int(e^{h(s-1/2)}-1)\,\nu(dh).
                                                               \tag{1}
\]
Indeed, both sides have the same second derivative. The right side
is convergent at \(h=0\) because the integrand is \(O(|h|)\).
An affine discrepancy in (1), after exponentiation on a grid, gives
only a scalar and a positive rescaling of the argument.

## A positive shift estimate at fixed degree

Let \(P(z)=\sum_{k=0}^m c_k z^k\) have degree \(m\ge2\) and all
roots on the unit circle, with circular angular mesh at least
\(\gamma/m\) for some \(0<\gamma\le\pi\). The following constants
depend only on \(\gamma\):
\[
 \begin{gathered}
 M=\frac{\pi}{\gamma},\quad a_0=\frac{1}{16M},\quad D=2+4M^2,
 \quad C=\frac12+\frac{\pi M}{2},\quad K=C+\frac{a_0 D}{2}+1,\\
 B=\exp(4a_0^2D+a_0),\qquad b_0=\tfrac12\exp(-4a_0^2D),\qquad
 V_m=M(1+\log m).
 \end{gathered}
\]
Suppose \(\nu\) is a positive measure on \(\mathbb R\setminus\{0\}\)
and define
\[
 \beta=m\int |h|e^{K|h|}\,\nu(dh),\qquad
 \rho=\frac{2B}{a_0 b_0}\,\beta\,e^{(K+V_m)\beta}.
\]
If \(\beta<\infty\) and \(\rho\le a_0\), the roots of
\[
       Q(z)=\sum_{k=0}^m c_k
                  \exp\left(m\int(e^{h(k/m-1/2)}-1)\,\nu(dh)\right)z^k
\]
pair with those of \(P\), with angular change at most \(\rho/m\)
for every pair. In fact each root \(z\) of \(P\) has exactly one
root of \(Q\) of the form \(z e^{y/m}\) with \(|y|\le\rho\).
When \(\beta=0\) the two polynomials agree.

To prove the estimate, fix a root \(z\) and put
\[
 F(w)=\frac{e^{-w/2}P(ze^{w/m})}{zP'(z)/m},\qquad
 q_\zeta=\cot\frac{\arg z-\arg\zeta}{2},\qquad
 b=\frac1{2m}\sum_{\zeta\ne z}q_\zeta,\qquad
 \widetilde F(w)=e^{ibw}F(w).
\]
Then \(F(w)=2m\sinh(w/(2m))\prod_{\zeta\ne z}
A_{q_\zeta}(w/(2m))\), where \(A_q(v)=\cosh v-iq\sinh v\).
At the \(j\)-th circular neighbor the absolute sine of the half-gap
is at least \(\gamma\min(j,m-j)/(\pi m)\).
Consequently \(|b|\le V_m\) and
\(\sum_{\zeta\ne z}(2+q_\zeta^2)\le Dm^2\). Here one can use
\(\sum_{j\ge1}j^{-2}\le2\) and the elementary harmonic sum bound.

For real \(h\), put \(v=h/(2m)\). The modulus identity
\[
      |A_q(v)|^2=1+(1+q^2)\sinh^2 v
\]
and \(|\sinh v|\le |v|e^{|v|}\) give
\[
 \left|\prod_{\zeta\ne z} A_{q_\zeta}(v)\right|
       \le \exp\big((C-1/(2m))|h|\big).
\]
Indeed, each neighboring pair contributes at most one factor
\(1+(\pi |h|/(2\gamma j))^2\) after extracting the exponentials.
The sum of the logarithms of these factors is bounded by the
integral of the same decreasing function on \((0,\infty)\), equal
to \(\pi^2|h|/(2\gamma)\).

For complex \(|y|\le2a_0\), put \(u=y/(2m)\) and
\[
          b_q(v)=\frac{\tanh v-iq}{1-iq\tanh v},\qquad
 \frac{A_q(v+u)}{A_q(v)}=\cosh u+b_q(v)\sinh u .
\]
Here \(|b_q(v)|\le\max(1,|q|)\le Mm\) and
\(|b_q(v)+iq|\le|\tanh v|(1+q^2)\).
Since \(|u|\le a_0/m\), the Taylor bounds
\(|\cosh u-1|\le|u|^2\) and \(|\sinh u-u|\le|u|^3\) give
\[
 \left|\log\big(\cosh u+b_q(v)\sinh u\big)-b_q(v)u\right|
                   \le4(1+|b_q(v)|^2)|u|^2 .
\]
For example the argument of the logarithm differs from one by
\(b_q(v)u\) and an error of magnitude at most \(2|u|^2\), so the
ordinary logarithm series applies. In the product \(e^{iby}\)
cancels the summed terms \(-iq_\zeta u\).
Together with the real shift bound, this proves
\[
 |\widetilde F(h+y)|\le B(|h|+2a_0)e^{(K-1)|h|},
 \qquad
 |\widetilde F(y)|\ge b_0|y| .
\]
For the lower bound take \(h=0\), so \(b_q(0)+iq=0\), and
use \(|2m\sinh u|\ge |y|/2\). Cauchy's derivative estimate on
disks of radius \(a_0\) around \(h+y\), followed by integration
along a real segment, now gives
\[
 |\widetilde F(h+y)-\widetilde F(y)|
                    \le (B/a_0)|h|e^{K|h|},\qquad |y|\le a_0 .
                                                               \tag{2}
\]
Here \(2a_0+|h|\le e^{|h|}\).

First truncate \(\nu\) to \(|h|\ge\varepsilon\). Let \(H\) be a
compound Poisson sum with jump intensity measure \(m\nu\) on the
truncation, and let \(J\) be the sum of its absolute jumps.
Expansion of the Poisson series gives
\[
 \mathbb E e^{xH}=\exp\left(m\int(e^{xh}-1)\,\nu(dh)\right),
 \quad
 \mathbb E[J e^{KJ}]\le \beta e^{K\beta}.
\]
For the second identity with inequality, differentiate the
exponential moment of \(J\), using \(e^{K|h|}-1\le K|h|e^{K|h|}\).
Take \(\chi=\mathbb E e^{-ibH}\). Its modulus is at least
\(e^{-|b|\beta}\) by the same series and \(1-\cos t\le |t|\).
The centered shifted polynomial, multiplied by \(e^{iby}\),
is \(\mathbb E[e^{-ibH}\widetilde F(y+H)]\) in the log coordinate
\(y\), up to the common normalization in \(F\). By (2) and
\(|H|\le J\), its difference from \(\chi\widetilde F(y)\),
divided by \(|\chi|\), is at most
\((B/a_0)\beta e^{(K+V_m)\beta}\). On \(|y|=\rho\) this is half
the lower bound \(b_0\rho\) for \(|\widetilde F(y)|\).
Rouché's theorem therefore gives one zero in the disk. All the
bounds are uniform in the truncation. The polynomial coefficients
and \(\chi\) converge as \(\varepsilon\downarrow0\) by the
weighted first-moment bound. The same strict boundary comparison
holds in the limit. The angular intervals of radius \(\rho/m\)
about the distinct reference roots are disjoint, so these disks
account for every root.

For the factorial symbol difference in (1), if all row rates
\(a_i\) exceed \(K\), the quantity in this estimate satisfies
\[
                  \beta\le \frac{m}{2nN}\sum_i\frac1{a_i-K}.
                                                               \tag{3}
\]
In particular for \(m=n\) or \(m=N\), fixed rows and a uniform
mesh constant \(\gamma\) give \(\rho/m=O(n^{-2})\), with explicit
constants above. This is a finite change of log-coefficient symbol
on either fixed grid. It does not compare the two different grids
at one common symbol.

## Actual mirrored reference profiles

For a fully mirrored profile with shared ordinary rows, the
[factorial circle-root theorem](MULTIVARIATE_MIRRORED_ROOT_CIRCLE.md)
writes its order-\(m\) normalized polynomial as a Schur–Szegő
composition of circle-rooted Gegenbauer factors with parameters
\(\alpha=(m\ell+h)/d\), \(1\le h\le d\).
The classical Schur–Szegő composition of two circle-rooted
polynomials preserves the lower bound on the angular mesh of
either factor. To see this, rotate one factor through any angle
strictly between zero and its mesh. The two root sets strictly
interlace on the circle. With phases chosen to make their centered
trigonometric forms real, every real pencil element is circle-rooted.
Grace's theorem carries that pencil through the composition with
the other factor. Composition commutes with rotation. For a
rotation with no shared roots in the output pencil, the two
members strictly interlace. This follows by mapping the circle to
the real line, where a real pencil with only real roots has a
quotient with one fixed imaginary sign in the upper half-plane,
forcing pole residues of the same sign. Repeated output roots
are excluded by using a sufficiently small rotation with no
shared roots. A real perturbation of either sign at a repeated
root with nonzero value of the other pencil member has a nonreal
zero for at least one sign. Therefore the rotated output pair
interlaces for a dense set of shifts up to the input mesh. If the
output mesh were smaller, a shift just above it would cross some
adjacent pairs while leaving others uncrossed, since the sum of
all output gaps is \(2\pi\). The two sets could not alternate.
This proves the mesh bound.

For \(\alpha\ge1\), the transformed Gegenbauer equation for
\(v(t)=(\sin t)^\alpha C_m^\alpha(\cos t)\) is
\[
        v''+\big((m+\alpha)^2+
              \alpha(1-\alpha)\csc^2t\big)v=0 .
\]
On the interval between any consecutive zeros, integration by
parts and the Dirichlet Poincaré inequality give length at least
\(\pi/(m+\alpha)\). This applies also from zero to the first
interior zero because \(v(t)=O(t^\alpha)\) there. Under
\(z=e^{2it}\) the circular mesh is thus at least
\(2\pi/(m+\alpha)\), including the wraparound gap by reflection.
If \(L=\min_i\ell_i/d_i\ge16\), choose one factor with ratio \(L\).
It gives mesh at least \(\gamma/m\) throughout the composition
with \(\gamma=2\pi/(L+3/2)\). Its parameter satisfies
\(\alpha\le mL+1\) and \(\alpha\ge1\). The indicated \(\gamma\)
makes
\[
        K=\tfrac32+(\pi/4+1/16)(L+3/2)+1/(8(L+3/2))
             <L+\tfrac12
\]
for \(L\ge16\), since \(\pi<22/7\).
Hence (3) applies with the actual mirrored order-\(n\) polynomial
as the reference and the upper symbol \(L_N\) placed on the
lower grid. The result gives \(O(n^{-2})\) matching of the two
fixed-grid root arguments for every such fixed profile. In this
application the perturbed coefficients remain real and reciprocal.
Reciprocal conjugation preserves each log disk about a reference
root on the circle, so uniqueness also keeps its paired root on
the circle.
To obtain consecutive coprimality one still needs a sufficiently
strong angular comparison across the two degrees with the
same log-factorial symbol.

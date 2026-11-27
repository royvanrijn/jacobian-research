# A single reflection reduces a three-cell maximizer to a monotone graph

Fix $0<\rho<1$. Suppose an actual balanced three-cell Gaussian
noise-stability maximizer in its active plane is invariant under
a reflection through the origin combined with a transposition of
two labels. Choose coordinates in which the reflection is
$R(x,y)=(-x,y)$ and the transposition exchanges labels $1,2$.
After selecting the signs of the axes, the partition has the form
\[
\begin{split}
 A_3&=\{(x,y):y>h(|x|)\},\\
 A_1&=\{(x,y):x>0,\ y<h(x)\},\\
 A_2&=\{(x,y):x<0,\ y<h(-x)\}.
\end{split}                                                    \tag{RF1}
\]
The height $h$ is finite and analytic on $[0,\infty)$, and
$h'(x)>0$ there, including the right derivative at zero.
Its slope has positive upper and lower bounds, and the graph
has a homogeneous limit under dilation by powers of $\rho$,
with bounded error in height. All three open cell representatives
are connected and there is exactly one
genuine triple junction, at $(0,h(0))$, with three distinct
tangent rays. The reduction applies to any competing global
maximizer at a first tie that has the stated reflection.

## Orientation and Gaussian rearrangement

Write $u=1_{A_3}$, $v=1_{A_1}-1_{A_2}$. Then $u$ is even and
$v$ is odd in $x$, and
\[
 S_\rho(A)=\langle u,T_\rho u\rangle+
       \tfrac12\langle1-u,T_\rho(1-u)\rangle+
       \tfrac12\langle v,T_\rho v\rangle.
                                                               \tag{RF2}
\]
On the folded half-plane $x>0$, let $d\gamma_+(x)=2\,d\gamma_1(x)$.
The even and odd kernels in the $x$ variable, relative to
$d\gamma_+$, are
\[
 K_\pm(x,s)=\tfrac12\big(M_\rho^{[1]}(x,s)
                             \ \pm M_\rho^{[1]}(x,-s)\big).
\]
Here $M_\rho^{[1]}$ denotes the one-dimensional Mehler density
relative to the Gaussian marginal. Both kernels are strictly
positive for $x,s>0$. Replacing $v$ on $x>0$ by $|v|$ and
extending it oddly preserves the cells' masses and the first
two terms in (RF2). The positive odd kernel, multiplied by the
positive Mehler kernel in $y$, strictly increases the last term
if both signs of $v$ occur in the half-plane on positive-measure
sets. Thus, up to exchanging $1,2$, label $1$ occupies the
right half-plane outside $A_3$ and label $2$ occupies the left.

Now keep all the Gaussian fiber masses of $u(x,\cdot)$ fixed.
The part of (RF2) depending on the shapes of these fibers is
the quadratic expression with kernel
$(3K_++K_-)\,M_\rho^{[1]}(y,t)/2$. In fact the linear terms
in the expansion of (RF2) depend only on the fiber masses.
The one-dimensional two-set Gaussian comparison bounds every
pair of fiber contributions by the contribution from two
same-direction half-lines of the respective masses. All such
bounds are attained simultaneously by upper half-lines. When
both fiber masses are strictly between zero and one, equality
forces same-direction half-lines by the two-set equality
theorem. The weight $3K_++K_-$ is strictly positive, so equality
in the integrated comparison forces equality for almost every
pair of nontrivial fibers. Such fibers occur on a set of
positive measure, since otherwise the partition depends only
on $x$, contrary to full active rank. By Fubini all nontrivial
fibers are half-lines of one common orientation. Reverse $y$
if needed. This gives (RF1) with a measurable, extended-real
threshold. The exact inputs here are Theorems 1.1 and 1.3 of
[Mossel and Neeman, *Robust optimality of Gaussian noise stability*](https://ems.press/content/serial-article-files/32017),
applied to ordinary indicators at fixed positive correlation.

## The translation equation fixes the sign of the slope

Use the [maximizer variation conventions](THREE_CELL_GLOBAL_OBSTRUCTIONS.md).
Let $\Sigma_+,\Sigma_-,\Sigma_0$ be the interfaces of $3$ with
$1$, of $3$ with $2$, and of $1$ with $2$, respectively.
Orient the first two out of cell $3$ and the last out of
cell $1$. Then $\Sigma_-$ is the reflection of $\Sigma_+$,
and $\Sigma_0$ lies on $x=0$ with normal $(-1,0)$.
It has positive surface measure by complete adjacency.

In the odd reflection block, a speed $f$ on $\Sigma_+$ has
reflected value $-f$ on $\Sigma_-$, with arbitrary speed on
$\Sigma_0$. Change the sign of the speed on $\Sigma_0$ as
coordinates for this block. In these coordinates the $\Sigma_+$
self kernel is $2M_\rho(p,q)-M_\rho(p,Rq)>0$, since both
points have positive $x$ coordinate. The $\Sigma_0$ self
kernel and the two cross kernels are positive multiples of
$M_\rho$. This follows from the incidence dot products
$d_+\cdot d_-=1$, $d_+\cdot d_0=-1$ and $d_-\cdot d_0=1$.
The folded weighted norm uses multiplicity two on $\Sigma_+$.

The horizontal translation speed belongs to this block and
has quadratic form zero for the first extra-coordinate
form $Q^{(1)}$ by (G1)–(G2). Its gauged speed on $\Sigma_0$
is $+1$. Taking absolute values in the folded coordinates
preserves its weighted norm and strictly increases the
interaction if there are positive and negative sets of
positive surface measure. That would contradict (G1).
Thus the horizontal outward normal component of $A_3$
on $\Sigma_+$ is nonnegative. It is strictly positive by
(G2), since the integral from $\Sigma_0$ is strictly positive
and all other folded contributions are nonnegative.
These operations take place in the weighted space of (G17),
so they also apply to the noncompact translation speed.
The distributional derivative of $u$ in $x>0$ is minus
this outward normal component times surface measure.
It follows that $u$ is nonincreasing in $x$ on almost every
horizontal line. Therefore $h$ has a nondecreasing
extended-real representative.

## Regularity, the axis, and connectedness

The maximizing score rule uses prices
$(\lambda_1,\lambda_2,\lambda_3)$ with $\lambda_1=\lambda_2$.
For example, symmetry permits averaging a valid price vector
with its swapped vector. The odd kernel shows
$T_\rho v(x,y)>0$ for every $x>0$. Thus on the right the
winner is decided by the sign of
\[
 F(x,y)=T_\rho(u-1_{A_1})(x,y)-(\lambda_3-\lambda_1).
\]
For every fixed source horizontal coordinate the integrand
$u-1_{A_1}$ is nondecreasing in its vertical coordinate.
It is nonconstant on a positive-measure set of such fibers.
Differentiating its one-dimensional Gaussian smoothing, or
convolving its nonnegative distributional derivative, gives
$F_y>0$ everywhere. The analytic implicit function theorem
therefore gives an analytic graph on the interval of finite
heights. A jump between two different heights would force
the continuous $F$ to vanish on the vertical interval
between them, contrary to $F_y>0$. On the graph its slope is
strictly positive by the translation identity above, which
holds pointwise by Gaussian integration by parts.

The nondecreasing threshold has $h(0+)>-\infty$. Otherwise,
on every bounded vertical interval both half-planes near the
axis would belong to $A_3$, leaving no $1/2$ interface. Also
$h(0+)<+\infty$, since otherwise $A_3$ would be null. Thus the
finite-height interval is $[0,X)$, and a finite limit at a
finite $X$ would continue the graph by the implicit function
theorem. On the axis $F(0,h(0))=0$ by continuity. At this point
$F_x=-\partial_x T_\rho1_{A_1}<0$: $T_\rho u$ is even, and
the Gaussian derivative for $1_{A_1}$ at $x=0$ integrates the
strictly positive source horizontal coordinate on that cell.
Consequently the graph extends analytically to the endpoint
with $h'(0+)=-F_x/F_y>0$. The reflected graph and the downward
vertical ray meet there. There are no other triple junctions
because $T_\rho v$ has the strict half-plane signs. The graph
description gives the claimed connectedness.

## The graph escapes in both coordinates

Write $p=(x,h(x))$ on $\Sigma_+$ and put $f=N_x>0$,
$w=-N_y>0$ there. Reflection and the vertical and horizontal
translation identities give, with integrals over $\Sigma_+$
unless indicated,
\[
\begin{split}
 a(p)w(p)/\rho&=\int (2M+M_R)w\,d\nu,\\
 a(p)f(p)/\rho&=\int (2M-M_R)f\,d\nu
                        +J_0(p),\qquad J_0(p)=\int_{\Sigma_0}M\,d\nu .
\end{split}                                                    \tag{RF3}
\]
Here $M=M_\rho(p,q)$, $M_R=M_\rho(p,Rq)$, and
$\Sigma_0$ is the downward ray from $(0,h(0))$.
The even block in the first equation adds the reflected
interface and has zero speed along $\Sigma_0$.

We use the following elementary iteration. Suppose a
continuous $m(t)>0$ on a tail satisfies
$m(t)\ge\int m(s)\,dP_t(s)$ with a nonnegative kernel
giving mass at least $1-Ce^{-ct^2}$ to $[\ell t,ut]$ for
each sufficiently large $t$, where $0<\ell<u<1$.
Then $m(t)$ is bounded below by a positive constant on
that tail. Indeed, retain only transitions in $[\ell t,ut]$
and normalize their mass to a probability. Stop at first
entry into $[\ell R,R]$, for large fixed $R$. The entry
occurs in finitely many steps since each retained
transition reduces the parameter by a factor at most $u$.
In reverse order, the parameters before entry are at least
$R,R/u,R/u^2,\ldots$. Thus the product of the retained
masses is bounded below by
$\prod_{j\ge0}(1-Ce^{-cR^2u^{-2j}})>0$, uniformly in the
starting point and the retained path. Iterating the
inequality gives this product times the minimum of $m$
on $[\ell R,R]$.

First suppose $X=\infty$ and $h(x)=o(x)$, which includes
a bounded height. By (RF3), $m=f/w=h'$ satisfies the
iteration inequality with kernel density
$(2M-M_R)w\,d\nu/\int(2M+M_R)w\,d\nu$. In horizontal
coordinate $s$ on the graph, $w\,d\nu=\phi(s)\phi(h(s))ds$.
Write $\sigma=\sqrt{1-\rho^2}$ and
$p_\sigma(t)=\sigma^{-1}\phi(t/\sigma)$. Then
\[
 M w\,d\nu
     =p_\sigma(s-\rho x)p_\sigma(h(s)-\rho h(x))\,ds .
\]
Choose $0<\ell<\rho<u<1$. The unrestricted denominator is
at least $\exp[-o(x^2)]$, by restricting to
$|s-\rho x|\le1$ and using the sublinear heights. The
part of the $M$ integral outside the good interval
$[\ell x,ux]$ is $O(e^{-c x^2})$ by the first Gaussian
factor. On the good interval $M_R/M=e^{-2\rho xs/\sigma^2}$
is exponentially small, and everywhere $0<M_R\le M$.
Therefore the good interval has mass at least
$1-Ce^{-c' x^2}$ for the iteration kernel. This forces
$h'$ to be bounded below by a positive constant for large
$x$, contradicting $h(x)=o(x)$.

For the other coordinate suppose $h$ reaches arbitrarily
large heights and its inverse $k(y)$ is $o(y)$, which
includes finite $X$. Apply (RF3) to $m=w/f=k'(y)$.
The iteration kernel now has numerator $(2M+M_R)f\,d\nu$
and denominator $\int(2M-M_R)f\,d\nu+J_0(p)$.
Parametrization by height $t$ gives
\[
 M f\,d\nu
     =p_\sigma(k(t)-\rho k(y))p_\sigma(t-\rho y)\,dt .
\]
The unrestricted numerator is at least $\exp[-o(y^2)]$
by $|t-\rho y|\le1$, while its mass outside $[\ell y,uy]$
is $O(e^{-c y^2})$, again using $M_R\le M$. Also
$J_0(p)=p_\sigma(-\rho k(y))\Phi((h(0)-\rho y)/\sigma)$
is $O(e^{-c y^2})$. The denominator is at most the
unrestricted numerator plus this last error. Hence the
same good-interval bound holds, forcing $k'$ to be bounded
below on the tail and contradicting sublinear growth.
In particular $X<\infty$ is impossible. With $X=\infty$,
the first argument excludes bounded height, and the two
arguments exclude sublinear growth of either function.

Once the graph escapes in both coordinates, the stationary
score also gives the sharper growth statement. At
$(x,h(x))$, $U_3-U_1$ is a constant in $(-1,1)$ and
$U_2\le\Pr(Y_x<0\mid X_x=x)\to0$. Hence $U_3$ stays
in $[\delta,1-\delta]$ for large $x$, for some $\delta>0$.
The conditional coordinates of $Y$ are independent Gaussians
with means $\rho x,\rho h(x)$ and variance $\sigma^2$.
Choose fixed $R,C$ so that the two coordinate deviations
exceeding $R,C$ have total probability less than $\delta$.
Using the increasing height, the score bound implies
\[
 h(\rho x-R)-C \le \rho h(x)
                      \le h(\rho x+R)+C                     \tag{RF4}
\]
for all sufficiently large $x$. Indeed, violating either
bound puts the entire conditional rectangle on one side
of the graph. In the lower bound put $t=x+R/(1-\rho)$
and subtract $C/(1-\rho)$ from the heights to obtain
$H(\rho t)\le\rho H(t)$. The function $H$ is positive on
a large tail since $h\to\infty$. Iterate until $t$ first
enters a fixed interval $[\rho t_0,t_0]$ on that tail.
This gives a positive lower bound on $H(t)/t$.
In the upper bound use $t=x-R/(1-\rho)$ and add
$C/(1-\rho)$ to the heights. The reverse inequality
$\rho H(t)\le H(\rho t)$ bounds $H(t)/t$ from above by
the corresponding maximum on that fixed interval. Thus
$c_1x\le h(x)\le c_2x$ for all sufficiently large $x$.

## Uniform slopes and the limit along powers of the correlation

Use $\vartheta=f/(f+w)$ on the right graph. Both
$\vartheta$ and $1-\vartheta$ are positive supersolutions
for the same nonnegative kernel
\[
       dP_p(q)=\frac{(2M-M_R)(f+w)(q)\,d\nu(q)}
                         {a(p)(f+w)(p)/\rho}.
\]
Indeed this follows from the two equations of (RF3),
discarding the source $J_0$ and the positive excess
$2M_R$ in the vertical equation.

Let $I(p)=\int M(f+w)d\nu$ on $\Sigma_+$. By (RF4) and
continuity there is a point on the graph within horizontal
distance $R$ and vertical distance $C$ of the conditional
mean $(\rho x,\rho h(x))$. Follow the graph from this point
until one coordinate first increases by one. This segment
stays in an enlarged fixed conditional rectangle and its
integral of $ds+dh$ is at least one. Since
$M(f+w)d\nu=p_\sigma(s-\rho x)p_\sigma(h(s)-\rho h(x))(ds+dh)$,
we obtain $I(p)\ge c>0$ uniformly for large $x$.

For $[\ell x,ux]$ with $\ell<\rho<u$, the part of
$I(p)$ outside the interval is exponentially small in
$x^2$. On the $ds$ part use the Gaussian tail in $s$;
on the $dh$ part use the supremum of the same Gaussian
density outside the interval and integrate the vertical
density. On the good interval the ratio $M_R/M$ is
exponentially small. The integral $J_0$ is also exponentially
small by its Gaussian factor in the horizontal coordinate.
Consequently the denominator in $dP_p$ differs from $2I(p)$
by at most $Ce^{-cx^2}$, and the good interval has kernel
mass at least $1-C'e^{-c'x^2}$. Apply the good-interval
iteration to $\vartheta$ and $1-\vartheta$. Both are bounded
below on the tail. This bounds $h'=\vartheta/(1-\vartheta)$
above and below by positive constants there. The compact
part of the graph has the same type of bounds by positivity
and continuity.

The slope bound in (RF4) gives
$|h(\rho x)-\rho h(x)|\le C_0$ for sufficiently large $x$.
Thus set
\[
 H_\infty(x)=\lim_{n\to\infty}\rho^n h(\rho^{-n}x),\qquad
 H_\infty(\rho x)=\rho H_\infty(x),\qquad
 h(x)=H_\infty(x)+O(1)\quad(x\to\infty).                     \tag{RF5}
\]
The limit is locally uniform for $x>0$ by summing the
successive errors $C_0\rho^n$. The same bound yields the
bounded error at large $x$. The two-sided slope bounds pass
to $H_\infty$, which extends continuously to zero with
$H_\infty(0)=0$. This gives invariance under the one
dilation factor $\rho$ in the limit. A linear limit would
require an additional argument.

## Total positivity in the odd reflection sector

The unrestricted planar second-variation operator in the odd
reflection sector has exactly one zero mode, the rotation.
It has a positive spectral gap after removing rotation and
horizontal translation. This statement concerns the unrestricted
operator on normal speeds. The mass-constrained Hessian is its
compression to the balanced speeds.

Fold the reflected interfaces onto the right graph and the
downward axis ray, giving each point coordinates $(X,Y)$ with
$X\ge0$. Their union is strictly ordered by $Y$, and $X$ is
nondecreasing along this order. On the graph multiply the odd
speed by $\sqrt{2}$, to absorb the multiplicity two. As above
change its sign on the axis ray. The folded kernel is then
$(2M-M_R)$ on graph pairs, $\sqrt{2}M$ on graph-axis pairs,
and $2M$ on axis pairs. After removing positive separable
factors it has the form
\[
 e^{\kappa Y_pY_q}H(\kappa X_pX_q),\qquad
 H(u)=2e^u-e^{-u},\qquad \kappa=\rho/(1-\rho^2).
                                                               \tag{RF6}
\]
Indeed, $H(0)=1$, so the additional $\sqrt{2}$ on each axis
argument is a separable factor. For $u\ge0$ put
$z=e^{-2u}/2$. Then
\[
   (u(\log H)')'=\frac{1-z^2-4uz}{(1-z)^2}>0,
\]
using $z^2\le1/4$ and $4uz\le e^{-1}$. Thus the $H$ factor
has nonnegative two-by-two minors on increasing nonnegative
inputs. The exponential factor has strictly positive
two-by-two minors in the $Y$ order, so the product is strictly
positive in these minors.

The operator with kernel divided by $\sqrt{a(p)a(q)}$ is
compact on $L^2(d\nu)$, modeling the interaction on the
weighted speed space. On the graph $a$ is bounded below:
the first equation of (RF3), (RF4), and the slope bound give
a fixed horizontal interval in a bounded rectangle about
the conditional mean on which the Gaussian integral is
bounded below. On the compact part use the nonvanishing
score derivative. On the axis, $a$ is the positive horizontal
derivative of the $1/2$ contrast. In the conditional product
law at $(0,y)$ this derivative is
\[
 \frac{\rho}{\sigma^2}\mathbb E[|Z_1|
            1_{\{Z_2<h(|Z_1|)\}}],
 \qquad Z_1\sim N(0,\sigma^2),\quad Z_2\sim N(\rho y,\sigma^2).
\]
For $y\le h(0)$ it is bounded below by replacing the event
with $\{Z_2<h(0)\}$. Now the square of the Mehler kernel
times the two Gaussian surface densities has Gaussian decay
in both arguments. The graph slope bound and the straight
axis give a finite Hilbert–Schmidt norm.

The principal eigenfunction is the horizontal translation
speed with eigenvalue $1/\rho$, by (G1)–(G2).
All other eigenvalues are at most one by (G17). There is an
eigenvalue one from the rotation speed: rotating the cells
and evaluating at their moving interfaces preserves their
stationary score differences, and differentiation gives
$K(Jp\cdot N)=a(Jp\cdot N)$ for the right-angle rotation $J$.
This speed is odd and nonzero on the axis. Polynomial growth
of the speed justifies differentiation by Gaussian dominated
integration. The exterior square of the compact interaction
operator has a strictly positive kernel on ordered pairs,
by (RF6). Its top eigenvalue is simple and its eigenfunction
has strict constant sign, by taking absolute values in the
Rayleigh quotient. Its eigenvalues are products of distinct
indexed eigenvalues of the interaction. That interaction is
positive semidefinite by Mehler factorization. Hence the top
product is $1/\rho$ times the second eigenvalue, which is one.
Thus one is simple and the ratio of its eigenfunction to the
principal one is strictly monotone. Compactness gives the
claimed gap below this eigenvalue.

The quotient of rotation by horizontal translation in the
gauged coordinates is $-Y$ on the axis and
$-h(x)-x/h'(x)$ on the graph. It is therefore strictly
decreasing in the folded order. Its graph derivative is
strictly negative as well. To see this, divide the integral
equations for the two eigenfunctions. The quotient at the
observation point is a positive constant times the expectation
of the quotient at the source point with respect to the kernel
tilted by the positive principal eigenfunction.
Differentiating this expectation is the covariance with the
logarithmic derivative of the kernel in the observation
coordinate. Along the graph, the nonseparable part of that
derivative is
$\kappa h'(x)Y_q+\kappa X_q(\log H)'(\kappa xX_q)$.
It is strictly increasing in the source order by (RF6) and
$h'>0$. The covariance is strictly negative, since the quotient
is strictly decreasing and the kernel has full support.
The slope bounds and Gaussian decay justify this differentiation.
In particular
\[
              xh''(x)<h'(x)(1+h'(x)^2)\qquad (x>0).           \tag{RF7}
\]

The even planar sector has only the right graph after folding,
with kernel $2M+M_R$. The same two-by-two argument applies,
using $2e^u+e^{-u}$ whose logarithmic derivative is positive
and increasing for $u\ge0$. The top eigenfunction is vertical
translation. Its only possible unrestricted zero mode is
simple and has a strictly monotone quotient by that
translation, again by the exterior-square argument.

### Circular geometry

The derivative of $x^2+h(x)^2$ has the sign of
$h(x)+x/h'(x)$. This last function is strictly increasing,
starts at $h(0)$ and tends to infinity. If $h(0)\ge0$, every
circle of radius $r>h(0)$ intersects the right graph once.
If $h(0)<0$, the right graph has a single closest point to
the origin, of radius $r_*<|h(0)|$. Each circle with
$r_*<r<|h(0)|$ intersects it twice; circles with
$r>|h(0)|$ intersect it once.

In both cases at almost every radius outside the core
$r\le|h(0)|$, the three circular sections are consecutive
arcs. The fixed label's arc is centered on the upward axis,
and the other two are mirror images. Inside the core, if
$h(0)>0$, circles with $r<h(0)$ contain only labels $1,2$,
each on one semicircle. If $h(0)<0$, circles below $r_*$
contain only label $3$, and those between $r_*$ and $|h(0)|$
contain one arc per reflected label and two arcs of label $3$.
These descriptions follow from the intersection counts and
the labels on the two points of the axis.

## A second reflection

If the maximizing partition also has a different central
reflection that transposes labels, then it is the simplex.
Indeed, if both transpositions are the same, their composition
is a nonidentity planar rotation fixing every moment vector,
contrary to full rank. Otherwise the composition cycles the
labels. Each circular mass must then be $1/3$ because the
rotation preserves radius. The [circular-balance comparison](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#circular-balance-and-compatibility-of-first-angular-moments)
forces the common simplex arcs. A single reflection alone
leaves a curved monotone graph as a possibility.

## The folded objective is not concave in the fiber masses

For a graph profile put $p(x)=1-\Phi(h(x))$ for $x>0$.
In (RF2) the nonlinear part is
$\frac12\iint (3K_++K_-)J_\rho(p(x),p(s))\,d\gamma_+(x)d\gamma_+(s)$,
where $J_\rho(u,v)$ is the joint hit probability of
same-direction Gaussian half-lines of masses $u,v$.
Even on nonincreasing profiles of total mass $1/3$, this
objective has positive second derivatives in some balanced
directions, at every fixed $0<\rho<1$.

For large $M$ take $h=-M$ on an initial half-line interval
$E$ of $\gamma_+$-mass
$p_E=(1/3-\Phi(-M))/(1-2\Phi(-M))$ and $h=M$ on its
complement. Then $\int p\,d\gamma_+=1/3$. Vary these
fiber masses by $t v$, with $v=1$ on $E$ and
$v=-p_E/(1-p_E)$ on its complement. Both signs of small
$t$ give admissible nonincreasing profiles of the same
mass. For a pair of heights $h_x,h_s\in\{-M,M\}$,
differentiation of the Gaussian joint probability gives
\[
 \frac{d^2}{dt^2}J_\rho(p(x)+tv(x),p(s)+tv(s))\bigg|_{t=0}
   =\frac{\phi((\rho h_x-h_s)/\sigma)}{\sigma\phi(M)}
     \big[2v(x)v(s)-\rho(v(x)^2+v(s)^2)\big],
 \qquad \sigma=\sqrt{1-\rho^2}.
\]
For same-height pairs the last factor is positive.
The scalar density for opposite-height pairs divided by
the density for same-height pairs is
$\exp(-2\rho M^2/(1-\rho^2))$. After taking out the common
same-height density, the positive within-block terms
stay bounded below as $p_E\to1/3$. The cross term tends to
zero since $3K_++K_-$ is integrable. Hence the second
derivative is strictly positive for sufficiently large $M$.
The obstruction applies to the fiber-mass objective on its
whole feasible domain. It leaves the stability-constrained
graph equation as a separate comparison.

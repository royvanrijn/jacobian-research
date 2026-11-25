# Actual contours and simple-fold exclusion for fully shared supports

Let \(a,b\) be a fixed ordinary fully shared exponent profile
with distinct columns and one or more varying rows. Use \(w\),
the phases \(\phi_\pm\), and the exact contour from the
[Rice representation](MULTIVARIATE_RICE_CONTOUR_AND_REAL_BOUNDARY.md).
The following concerns compact neighborhoods of interior
points of \(-\pi<\Im w<0\). The single unit-gap row is
already covered by the ordinary univariate theorem, so
assume the sum of the absolute gaps is at least two.

**Theorem.** At every interior simple ramification value
\(w_c\) of the fully shared stationary map, a fixed
physical neighborhood of \(w_c\) excludes consecutive
common factorial roots for all sufficiently large orders,
for both equal and unequal total degrees. This is an
exclusion in the full physical neighborhood, including
unbounded scaled Airy fields.

**Theorem.** For every fixed ordinary fully shared profile
with unequal total degrees and \(\max\rho\ne1\), consecutive
factorial moment polynomials are coprime at all sufficiently
large orders. For equal degrees with \(\max\rho<1\),
every possible sequence of large-order common roots can
accumulate only at nondegenerate two-saddle phase matches
with matching squared Gaussian amplitudes.

## Regular values and nonzero coefficients

The real boundary exclusion and the infinity estimates of the
[supercritical contour proof](MULTIVARIATE_SHARED_SUPERCRITICAL_EVENTUAL_COPRIMALITY.md#descending-paths-and-stirling-bounds-at-infinity)
apply independently of \(\max\rho\). Their hypotheses are
positive integer exponents in each row and total absolute gaps
at least two. They identify all finite cut obstructions and
show that unbounded descending paths have \(h=\Re\phi\to-\infty\).
There are no lower stationary points and the proper upper
map accounts for both upper stationary points with
multiplicities.

At a regular parameter with two equal saddle heights, use a
small quadratic rectangle at each saddle and a common level
\(c\) just below their heights. As in the
[one-saddle isolation](MULTIVARIATE_SHARED_SUPERCRITICAL_EVENTUAL_COPRIMALITY.md#isolation-of-the-larger-saddle-through-transitions),
the union of the two rectangles and \(S_c=\{h<c\}\)
is invariant under descent. A common finite flow time
puts the exact compact loop in this union. Decompose
the resulting local crossings into integer multiples of
the two descending saddle segments plus a compact chain
in \(S_c\). This decomposition persists in a small
parameter neighborhood, with a strict uniform loss
on the exterior chain. It gives the actual Gaussian
expansion at each saddle with its integer coefficient
and a uniform remainder of relative order \(O(n^{-2})\)
after the first relative \(1/n\) correction. This follows
by expanding the analytic amplitude and the Stirling
series by two terms in the exact quadratic coordinate,
with the odd local Gaussian moments zero.

At a point of strict unequal heights, first isolate the
larger saddle by the one-rectangle argument. If its integer
coefficient vanishes, the whole contour can be replaced
by a compact chain below its height and the same argument
now isolates the smaller saddle. It cannot flow back up
to the larger one. These local expansions give zero-freeness
in a small disk if the highest saddle that occurs has no
other active saddle of the same height.

The integer coefficients cannot both vanish in either
construction. Otherwise the contour could be put into
a sublevel below both heights. Then the descending flow
could push it below any prescribed finite level, since
there is no stationary point at or below its starting
level, no accumulation at the finite cut, and unbounded
ends have phase tending to negative infinity.
Choose a bounded real gap, either \((p_L,0)\) or \((1,p_R)\).
At least one such gap is bounded. On its closure the real
phase extends continuously and is bounded below. The
original exact loop winds once about the endpoint in
\([0,1]\) and zero times about the adjacent stationary branch point.
Any closed chain entirely below the gap's minimum height
avoids this connecting segment, so its winding numbers
at the two endpoints must agree. This contradiction
proves at least one nonzero coefficient. The identical
reason applies to a single degenerate upper saddle
after localization to its cubic neighborhood below.

## Localizing the exact contour at a fold

At a simple ramification point \(q_c\) with \(X(q_c)=w_c\),
the first two \(q\)-derivatives of the upper phase vanish
and its third derivative is nonzero. There are no other
stationary points. Use an analytic cubic coordinate \(u\)
near \(q_c\), with analytic coefficients \(C(w),\lambda(w)\)
such that
\[
 \phi_+(q;w)=C(w)+u^3/3-\lambda(w)u,
 \qquad \lambda(w_c)=0,\quad \lambda'(w_c)\ne0.              \tag{1}
\]
For clarity one can construct this form by the analytic
implicit-function argument for a simple cubic turn.
At \(w_c\), absorb the nonzero cubic unit into the local
coordinate \(v\), making the phase \(C(w_c)+v^3/3\).
Solve (1) with \(u=v+g(v,w)\) in the Banach algebra
of series summable with weights \(R^j\) on a sufficiently
small closed disk. The linearized map in \((g,C,\lambda)\)
is \(v^2 g+C-\lambda v\); its inverse is bounded by
removing the constant and linear coefficients and dividing
the remainder by \(v^2\). The analytic implicit-function
theorem applies. On a smaller disk \(u\) is univalent
because its derivative is uniformly close to one.
The mixed derivative \(\phi_{qw}=-1\) gives
\(\lambda'(w_c)=dq/dv(0)\ne0\) by comparing the linear
coefficient in \(v\) after differentiation in \(w\).
This is the local cubic transformation of
[Chester, Friedman and Ursell](https://doi.org/10.1017/S0305004100032655).

At \(w_c\), use a small neighborhood inside the coordinate
disk given by \(|\Re u^3|<b\), \(|\Im u^3|<b\).
It is simply connected, being star-shaped with respect
to zero. In this coordinate along the descent
\(\dot u=-|u'(q)|^2\overline{u}^2\) one has
\(d(u^3)/dt=-3|u'|^2|u|^4\).
The neighborhood can therefore be left only through
\(\Re u^3=-b\), where the height is strictly below
\(c=\Re C(w_c)-b/6\). The union of the neighborhood
and \(S_c\) is invariant. The same finite-time argument
puts the exact loop in this union. In the local
neighborhood its intersection with \(S_c\) has three
components, the cubic inverse images of
\(\Re u^3<-b/2\) in the defining rectangle. Connect
the endpoints of each component to a fixed point on
its central descending ray by paths in that component.
The resulting local contour is a finite integer
combination of paths joining the three rays. Since the
neighborhood is simply connected its local class is
specified by its oriented endpoints in the three
components. If this class were zero the loop could
be put into \(S_c\), which is impossible by the winding
argument in the absence of any remaining saddle below
the fold. Thus the actual local class is nonzero.

Use the canonical coordinate (1) for nearby \(w\) and
attach the same exterior chain to the moving endpoints
of the local pieces by short connecting paths. The
exterior paths and the short pieces have a fixed loss
below \(\Re C(w)\) on a sufficiently small disk.
The local class is the same throughout the disk.
On the fixed central rays of its three cubic decay
sectors one may choose a common physical cutoff radius
within the analytic coordinate disk, with the endpoints
strictly below the fold height. Caps and connectors
used below can be taken at that radius inside the decay
sectors. This is an actual contour identity for both
orders, since the homotopies stay in the exact Rice domain.

Let \(U(\xi)\) be the integral of
\(e^{v^3/3-\xi v}\) along the complete cubic class.
In the basis of one standard Airy contour and its
rotation by \(\omega=e^{2\pi i/3}\), the integer
coefficients are not both zero. The corresponding
functions are \(2\pi i\,\mathrm{Ai}(\xi)\) and
\(2\pi i\,\omega\mathrm{Ai}(\omega\xi)\). In particular
\(U(0)\ne0\), since \(\mathrm{Ai}(0)\ne0\) and
\(1,\omega\) are linearly independent over \(\mathbb R\).
Integration by parts gives \(U''=\xi U\).
Thus no zero of \(U\) is multiple.
[Airy contour integral](https://dlmf.nist.gov/9.5.E4)

## A uniform joint norm for the actual cubic class

Write \(\xi=n^{2/3}\lambda\),
\(J=n^{-1/3}U(\xi)\), and \(K=-n^{-2/3}U'(\xi)\).
These are the complete local integrals with amplitudes
one and \(u\), respectively. For \(|\xi|\) larger
than a fixed \(R_0\), choose the two rotations
\(\omega^j\xi\) with the smallest absolute principal
arguments. Ties may be broken either way. Both arguments
lie in \([-2\pi/3,2\pi/3]\), and their angles differ
by \(2\pi/3\) without crossing the principal cut.
The two corresponding oriented contours form an integer
basis of the relative classes of three cubic ends.
Express the actual class in this basis, with coefficients
\(c_j\) not both zero; only finitely many pairs of
coefficients are needed. Put
\[
 V(\xi)=|\xi|^{-1/4}\sum_j |c_j|
               \exp\!\left(-\Re\frac23(\omega^j\xi)^{3/2}\right).
\]
The sectorial leading Airy and Airy-derivative estimates give
\[
 |U(\xi)|+|\xi|^{-1/2}|U'(\xi)|\ge c_1 V(\xi)
 \quad (|\xi|\ge R_0).
\]
Indeed, after scaling the derivative by a square root, its
two leading terms have opposite relative signs: the derivative
ratio for \(\omega^j\mathrm{Ai}(\omega^j\xi)\) is
\(-\omega^j\sqrt{\omega^j\xi}\) to leading order.
These signs are opposite for the two rotations, which gives
the elementary comparison \(|P+Q|+|P-Q|\ge |P|+|Q|\)
on the two leading values. The asymptotic errors are uniformly
small relative to the sum on these closed sectors for large
fixed \(R_0\). The sectorial bounds are the same as those
used in the [Airy contour envelope](EXPLICIT_AIRY_CONTOUR_ENVELOPE.md#3-explicit-primary-airy-error-bounds-when-xige4).

For each of the two basis integrals choose the compatible
contour in that envelope and rotate it back. Its absolute
integral, with weight \(|c_j|\), is bounded by a constant
times the corresponding term of \(V\). In the physical
coordinate \(u=n^{-1/3}v\), the exits at the chosen fixed
radius lie in closed subsectors of cubic decay for
\(\lambda\) sufficiently small and \(n\) large. Joining
to the original endpoints on the same circle retains a
fixed negative cubic margin. To verify the tail bound
one may continue from those exits by outward rays in the
same closed decay subsectors. Connecting within a sector
at the cutoff keeps \(\Re(u^3/3-\lambda u)\) below a fixed
negative constant. The full integrals of the polynomials
one and \(u\) along these contours still equal \(J,K\).
Their omitted rays have exponentially small absolute
integrals. The compatible contours before truncation
and the fixed-radius compatibility are the
nonnegative-sector bounds in the cited envelope proof.
For \(|\xi|\le R_0\) one can instead use the fixed central
rays, for which the absolute integrals in the \(v\)
coordinate are bounded on the compact disk.

After taking \(|\lambda|\le1\), this gives a bound \(T\)
for the sum of the absolute integrals of the local pieces:
\[
 T\le C_1
 \begin{cases}
 n^{1/3}(|J|+|K|),&|\xi|\le R_0,\\
 |\lambda|^{-1/2}(|J|+|K|),&|\xi|>R_0 .
 \end{cases}                                                \tag{2}
\]
For the first range use the positive minimum of
\(|U|+|U'|\) on the compact disk. In the second range
use \(n^{-1/3}|\xi|^{1/2}=|\lambda|^{1/2}\) and the
lower bound just proved. Any added exponentially small
caps can also be absorbed into (2). In fact both cases give
\[
 |J|+|K|\ge c_2 n^{-2/3}
                     \exp(-C_2 n|\lambda|^{3/2})
\]
on the fixed physical disk, after decreasing \(c_2\) if
needed. For \(|\xi|>R_0\) the bound from one nonzero basis
coefficient is at least
\(c n^{-2/3}|\xi|^{1/4}\exp(-2|\xi|^{3/2}/3)\).
Shrinking that disk absorbs every exterior contribution
and fixed-radius boundary term relative to the joint norm.
The joint norm is never zero, by the Airy equation.

## Comparing both orders on the same cubic contour

At a proposed common root use the same local \(u\) coordinate
for both orders, at the lower-order parameter \(w\).
After removal of nonzero scalars and \(e^{n C(w)}\), the
lower integrand is \(e^{n(u^3/3-\lambda u)} A_{n,0}\).
For the upper integrand evaluated at
\(w'=w-\Delta\log(1+1/n)\), use the same exponential phase
with the amplitude \(A_{n,1}=A_{n,0} B_n\).
Here \(A_{n,0}\) tends to the analytic nonzero amplitude
\(A(q)\,dq/du\), uniformly on the local chart, and
\[
 B_n(u,w)=\exp\bigl(\phi_+(q;w)+\Delta q\bigr)
                   \bigl(1+O(n^{-1})\bigr).
\]
Indeed
\((n+1)\phi_+(q;w')-n\phi_+(q;w)
=\phi_+(q;w)+\Delta(n+1)\log(1+1/n)q\).
The uniform estimates also hold after \(u\) derivatives
on a smaller chart, by Cauchy's formula.

With \(v=\sqrt\lambda\), divide each analytic amplitude by
\(u^2-\lambda\) after interpolation at \(v,-v\):
\[
 A_{n,j}(u)=a_{n,j}+b_{n,j}u+(u^2-\lambda)R_{n,j}(u)
 \qquad (j=0,1).
\]
The interpolation is analytic even at \(\lambda=0\).
The divided functions and their first derivatives are
bounded uniformly on a smaller disk. On a fixed outer
circle the denominator is bounded away from zero, and
the maximum principle and Cauchy's formula give these
bounds on the interior. Integrate the last term by parts
on the truncated local paths, using
\((e^{n(u^3/3-\lambda u)})'=n(u^2-\lambda)e^{n(u^3/3-\lambda u)}\).
The bounds in (2) and the contour loss show that the two
normalized moment values are
\[
 a_{n,j}J+b_{n,j}K+O(T/n)+O(e^{-c_3 n}).
                                                               \tag{3}
\]

For \(\Delta\ne0\), the limiting determinant at the fold
is nonzero. If \(a_*,b_*\) are the value and first \(u\)
derivative of \(A(q)dq/du\) there, the limiting second
row of the matrix in (3) is
\((a_* B(0),b_* B(0)+a_* B'(0))\), where
\(B=e^{\phi_++\Delta q}\).
The determinant is \(a_*^2 B'(0)\), and
\(B'(0)=\Delta(dq/du)(0)B(0)\ne0\) because
\(\phi_u(0,w_c)=0\). The inverse matrix is uniformly
bounded on a small physical disk at large \(n\).
At a common zero (3) would make \(|J|+|K|\) at most
\(O(n^{-2/3})(|J|+|K|)\) plus an exponentially small
relative term, by (2). This is impossible.

For \(\Delta=0\), the determinant at general \(\lambda\)
is
\[
 \frac{A_{n,0}(v)A_{n,0}(-v)}{2v}\,
                         \bigl(B_n(v)-B_n(-v)\bigr).
\]
The limiting difference divided by \(2v\) is
\(-\tfrac23 e^{C(w)}\lambda+O(|\lambda|^4)\),
by the phase values \(C(w)\mp 2v^3/3\).
The finite-\(n\) error here is \(O(n^{-1})\)
uniformly by the \(u\)-derivative bound for \(B_n-B\).
Thus the inverse matrix has norm \(O(|\lambda|^{-1})\)
when \(|\xi|\ge R_1\), after shrinking the physical disk
and taking \(n\) large. Dividing (3) by the joint norm
would then give \(1\le O(1/(n|\lambda|^{3/2}))+o(1)\)
at a common zero. Choose \(R_1\ge R_0\) sufficiently
large to exclude this possibility throughout the
unbounded scaled range.

## Bounded scaled fields for equal total degrees

For the remaining bounded \(\xi\)-range, let
\(\epsilon=n^{-1/3}\) and regard \(\epsilon>0\)
temporarily as real. Pull back \(w\) by
\(\lambda(w)=\epsilon^2\xi\).
Use only the truncated local contour, dropping the exact
exterior chain for this real-order interpolation.
At integer orders the dropped part is exponentially small
relative to the local normalization, uniformly on compact
\(\xi\)-sets. The local normalized function is
the integral after substitution \(u=\epsilon v\) and
division by \(\epsilon\). Its limit is
\(A(q_c)(dq/du)(0)U(\xi)\). The convergence holds with the
\(\xi\)-derivative on compact sets. The \(\epsilon\)-derivative
is bounded there. To check the last assertion, integrate on the
fixed cubic rays in the \(v\) coordinate up to length
\(r_0/\epsilon\). They have a uniform integrable majorant
\(e^{-c|v|^3+C|v|}\).
Differentiating the leading amplitude at
\((u,w)=(\epsilon v,w(\epsilon^2\xi))\)
introduces at most \(O(1+|v|)\). The Stirling error on the
fixed local upper disk is \(O(1/n)\), analytic also in
a sector of complex \(n\) around the positive real axis.
Its \(n\)-derivative is \(O(n^{-2})\) by Cauchy's estimate,
so its \(\epsilon\)-derivative at fixed \(u,w\) is \(O(\epsilon^2)\).
The derivatives through \(u,w\) obey the same polynomial
majorant by the fixed-disk bounds. The moving endpoints
are exponentially small after any polynomial factor
from differentiation. This establishes the derivative
claims for the local interpolation without interpolating
any residue identity away from integer orders.

Each simple zero \(\xi_0\) of \(U\) in a fixed compact
\(\xi\)-range therefore has a local interpolating branch
\(\xi(\epsilon)=\xi_0+O(\epsilon)\) with bounded real-order
derivative, by the implicit function theorem and uniform
convergence. There are only finitely many such branches
on a slightly larger compact range. Rouché's theorem,
with the nonzero derivative on each simple branch, shows
that the actual integer-order roots differ in \(\xi\)
from these branches by an exponentially small quantity.
Outside disjoint small disks about the model roots the
limiting function is bounded away from zero.
Since \(U(0)\ne0\), the derivative of the physical
\(\lambda\)-value along each interpolating branch is
\[
 \frac{d}{d\epsilon}\bigl(\epsilon^2\xi(\epsilon)\bigr)
                         =2\epsilon\xi_0+O(\epsilon^2),
 \qquad \xi_0\ne0 .
\]
Between \(n^{-1/3}\) and \((n+1)^{-1/3}\) the resulting
complex displacement has nonzero leading term of order
\(n^{-5/3}\), dominating the exponentially small errors.
Different branches cannot agree in the physical parameter
at adjacent orders either, by their fixed separation in
\(\xi\) and the fact that the ratio of the two scaling
factors tends to one. With \(\Delta=0\) the physical
parameter \(w\) is the same at common roots, so this proves
exclusion in every fixed bounded scaled field.
Together with the large scaled estimate it proves the
first theorem in a fixed physical disk.

## Consequences away from the real cusp

If the total degrees differ, the
[regular phase comparison](MULTIVARIATE_REGULAR_SADDLE_PHASE_SEPARATION.md)
excludes a limiting common root at any regular two-saddle
scale where both saddles occur and have equal heights.
If only one contributes or their heights differ, the
actual expansions give local zero-freeness. The complex
simple fold is excluded by the preceding proof. A sequence
of possible common roots at arbitrarily large orders
has a limiting scaled root in the fixed annulus by the
[fixed-profile localization](MULTIVARIATE_SHARED_SUPPORT_CONDITIONAL_SADDLE.md#a-fixed-annulus-for-possible-common-roots).
If \(\max\rho\ne1\), every real limit in this annulus
is excluded by the noncusp real-scale argument. The
complex limits have all been excluded above, proving
the unequal-degree assertion.

For equal total degrees, use the same regular localization
and suppose that common roots tend to a regular equal-height
point \(w_0\) where both actual coefficients are nonzero.
Let \(c_L,c_R\) be their integer crossing coefficients
including orientations, and \(A_L,A_R\) the Gaussian
amplitudes with fixed local orientations. At any zero
near \(w_0\) the expansion gives
\[
 e^{nD(w)}
  =-\frac{c_L A_L(w)}{c_R A_R(w)}
                \bigl(1+B(w)/n+O(n^{-2})\bigr),
\]
uniformly, with an analytic \(B\) and nonzero analytic
amplitudes. At a common root, \(w\) is the same for
both orders. Dividing the two equations cancels the
first relative correction and gives
\(e^{D(w)}=1+O(n^{-2})\).
Since \(D'=q_L-q_R\ne0\), it follows that
\(D(w)=D(w_0)+O(n^{-2})\) with
\(D(w_0)\in2\pi i\mathbb Z\). Then \(e^{n D(w)}\to1\),
which forces \(c_L A_L(w_0)=-c_R A_R(w_0)\).

In particular the equality of the squared amplitudes
requires the following rational invariant of the
stationary branches, with their actual contour weights:
\[
 c_L^2 \frac{\prod_i x_i(q_L)}{1-\rho(q_L)}
   = c_R^2 \frac{\prod_i x_i(q_R)}{1-\rho(q_R)}.
\]
Indeed the square of the leading Gaussian coefficient is
proportional to \(A(q_j)^2/X'(q_j)\) with one common
factor. The squared Rice amplitude is proportional to
\(\prod_i x_i(q_j)/(q_j(q_j-1))\) and
\(q_j(q_j-1)X'(q_j)=1-\rho(q_j)\).
This proves the necessary condition in the second
theorem. In the strictly subcritical case all real scales
and the two conjugate fold neighborhoods are already
covered. Therefore possible large-order limits are
confined to the stated regular matches. The phase matches
are isolated away from the fold, because \(D'\ne0\).
Together with the compact annulus and the local exclusions,
this leaves only finitely many possible limit points for a
fixed strictly subcritical equal-degree profile.
The [active-arc classification](MULTIVARIATE_SHARED_SUBCRITICAL_EQUAL_DEGREE.md)
excludes these remaining limits in the strictly subcritical
equal-degree case.

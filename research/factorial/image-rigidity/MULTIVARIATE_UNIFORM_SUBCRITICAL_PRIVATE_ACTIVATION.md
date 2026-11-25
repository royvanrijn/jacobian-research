# Finite contour activation for a varying subcritical private endpoint

Let the only private row of an ordinary exponent profile be
\((0,1)\), and delete equal rows. For the shared rows put
\(d_i=b_i-a_i\) and
\[
 S_{\rm curv}=1+\sum_{\rm shared}d_i^2/\max(a_i,b_i),\qquad
 \sigma=\sum_{\rm shared}d_i^2/a_i,\qquad
 \Delta=1+\sum_{\rm shared}d_i.
\]
Consider varying profiles with \(S_{\rm curv}\) bounded,
\(\sigma\to\sigma_0<1\) and terminals \(n+1\to\infty\).
At a physical parameter \(t\) use the same upper-half-plane
branch for both centered fields
\[
             \xi_j=\log(j^\Delta t)+\sum_{\rm shared}d_i\log a_i,
                   \qquad j=n,n+1 .
\]
**Theorem.** A consecutive common root along such a sequence
must have both \(\xi_n,\xi_{n+1}\to0\). The result includes
\(\sigma_0=0\) and every total-degree difference.

We work in the upper parameter half-plane by conjugation.
The [bounded-private annulus and regular comparison](MULTIVARIATE_UNIFORM_ONE_PRIVATE_REGULAR_COMPLEMENT.md)
give bounded real parts of \(\xi_j\), exclude the negative
real parameter boundary, and force a limit at height zero
to have field zero. We prove exclusion at all interior heights
\(0<\Im\xi<\pi\), including those where a straight left Barnes
line lacks uniform decay.

## Phases and finite descent

Extract the same retained and diffuse limit as in the regular
comparison. Write
\[
 B(q)=\sum_{\rm retained}d_i\log(1+d_iq/a_i)+\tau q,\qquad
 X(q)=\log(1-q)+B(q),\qquad
 \Phi(q;\xi)=\xi q+\int_0^q X(z)\,dz ,
\]
on the upper chart. The lower chart has the conjugate logarithms
and the additional phase term \(-2\pi i q\). The pure \(a\)
boundary has phase zero in each respective chart. The only
saddle is the upper \(q=q(\xi)\) from
\(X(q)=-\xi\), by the
[one-sheet compactness theorem](MULTIVARIATE_UNIFORM_SADDLE_MAP_COMPACTNESS.md).
Define its stationary value \(D(\xi)=\Phi(q(\xi);\xi)\).
Then \(D'(\xi)=q(\xi)\), with \(\Im q>0\), so \(D\) is
injective on the parameter strip. Indeed the difference
quotient between any two points of the strip is the integral
of \(q\) on the segment and has positive imaginary part.
There is a local inverse at the positive endpoint
\(\xi=0,q=0\) since \(X'(0)=\sigma_0-1<0\).
The inverse enters the upper \(q\)-plane for upper \(\xi\),
and \(D(\xi)\to0\) at this endpoint. Injectivity therefore
excludes \(D=0\) at an interior field. Otherwise its open
inverse neighborhood would contain the images of fields
approaching the endpoint, contradicting injectivity.

The right gap is \((1,p_R)\), with \(p_R=+\infty\) if
there is no retained negative row. The upper and lower
phases join through this gap with common real part and
derivatives. Near \(q=0\) the upper phase extends with
derivative \(\xi\), and the lower phase with derivative
\(\xi-2\pi i\). The two descent paths starting at these
two local copies of zero point strictly into their respective
half-planes. They begin at height zero in real phase.
The lower flow cannot return to the upper page.

In this two-chart domain every trajectory staying above a
fixed real-phase level has a stationary accumulation point
unless it enters the sublevel of that level. The finite
real boundary exclusions are those of the
[regular comparison](MULTIVARIATE_UNIFORM_ONE_PRIVATE_REGULAR_COMPLEMENT.md#the-actual-right-saddle-on-compact-field-disks).
There is no logarithmic singularity at \(q=0\), where
the upper velocity into the upper chart is \(\Im\xi>0\).
The remaining finite branch points have the same logarithmic
coefficients and signs as there. The lower vertical speed is
bounded away from zero. At upper infinity the same argument
applies if there are retained shared rows or positive diffuse
mass. For clarity, retained gaps of both signs contribute
at least \(\pi-o(1)\) in total imaginary row logarithms.
If all retained gaps are positive, real descent points left
for large positive real part, so that part of a trajectory
is bounded above. At infinity with real part bounded above,
the private gap and any retained shared gap provide at
least \(\pi-o(1)\) in those row logarithms before subtraction
of the rational arguments. With positive diffuse mass and
no retained shared row, the vertical velocity is positive
at large height from the diffuse term, and at bounded
height and far left from the real-cut argument of the
private row. In the flat limit with no shared or diffuse
term, \(\Phi'_+(q;\xi)=\xi+\log(1-q)\), whose real part
tends to positive infinity as \(|q|\to\infty\). Thus in
every case the speed is bounded away from zero on the
far portion of the trajectory. The finite energy identity
then excludes escape above a finite level, including
escape in finite time. Indeed the total time and then the
length outside a large disk would be finite. A bounded
trajectory away from the excluded real points has a
stationary accumulation point by the same energy identity.

## A finite deformation with two endpoint traces

Use the exact finite Rice rectangle with real sides
\(-c/j\) and \(c_R\), for fixed \(0<c<1\) and
\(1<c_R<p_R\), and any fixed positive top and bottom heights.
The first private Gamma pole is \(-1/j\), so the exact
rectangle encloses the coefficient poles and no other
poles. Keep the central left segment
\[
          \Re q=-c/j,\qquad |\Im q|\le Y_j,\qquad
          Y_j=C\log(j)/j ,
\]
where the positive constant \(C\) will be taken large.
The complementary arc starts at the lower end of this
segment and finishes at its upper end. Its limiting arc
starts and finishes at the two copies of zero.

Assume first \(\Im D(\xi)\ne0\).
Flow the limiting complementary arc forward by descent,
including its two endpoints in the local regular charts.
Neither endpoint trace can approach the saddle. The lower
one stays below, and the upper one has imaginary phase
zero until any crossing to the lower chart. Choose a
level below both zero and \(\Re D(\xi)\). As in the
[finite deformation proof](MULTIVARIATE_UNIFORM_ONE_PRIVATE_REGULAR_COMPLEMENT.md#the-actual-right-saddle-on-compact-field-disks),
a small quadratic rectangle at the saddle, together with
this strict sublevel, is forward invariant. By the preceding
compactness argument, every point of the initial arc enters
this union after finite time. Both endpoint traces themselves
enter the strict sublevel. Compactness of the initial arc
gives one finite time with both properties. The homotopy
identity replaces the initial arc by the two oriented
endpoint traces and the flowed arc. The flowed arc is
an integer multiple of a local descending saddle segment
plus strict-sublevel pieces. This follows by separating the
sublevel into its two descending components inside the
quadratic rectangle, as in the cited finite deformation
argument. The original positively oriented contour includes
the downward central left segment.

All this deformation is finite. The original exact arc
starts at \(-c/j-iY_j\), ends at \(-c/j+iY_j\), and
converges to the limiting arc. Flow it using the same
limiting phase and the same finite time. Near the two
endpoints the phase derivative is regular and its descent
velocity points into the corresponding half-plane. Away
from the endpoint collars the finite deformation stays
in a compact subset of the two half-planes joined through
the right gap, which avoids every exact pole for the
approaching profiles. Thus the exact homotopy is valid,
with the initially retained central segment fixed.
One may move the local saddle segment to the nearby
finite-profile saddle. The strict sublevel and nonzero
Gaussian bounds are unchanged on a smaller field disk.

The saddle integer is locally constant wherever both
endpoint traces enter a strict sublevel of the saddle.
To check this, at such a field choose a common finite
time, a quadratic rectangle with exit faces strictly
below the stopping level, and the resulting finite
crossing decomposition. They transport to nearby fields
with a small deformation of the saddle rectangle.
The oriented crossing number is unchanged. Longer flow
times yield the same number because the two rectangle-side
sublevels are forward invariant and the endpoints are
already in the sublevel. One can compare different
quadratic rectangles and stopping levels by using a common
smaller rectangle, a lower stopping level and longer
finite time. Equivalently the saddle number is the
coefficient of the oriented descending segment in the
relative contour class after keeping the endpoint traces.
For nearby finite profiles the same finite homotopy gives
that integer.

The two endpoint traces beyond their starts are negligible
relative to the central segment after normalization by the
pure \(a\) moment. On their initial pieces with small
\(|\Im q|\) the finite-profile phase obeys
\[
       \Re\Phi_{\pm,j}(q;\xi_j)\le -c_1|\Im q|,
                              \qquad |\Im q|\ge Y_j.
\]
Indeed at the starting point \(-c/j\pm iY_j\), the
negative linear loss is at least a fixed multiple of
\(Y_j\), since the imaginary derivative tends to
\(\Im\xi\) or \(\Im\xi-2\pi\). The quadratic error is
\(O(Y_j^2)\). During the initial limiting flow the
real part of the nearby finite-profile phase decreases
at a rate bounded away from zero, by convergence of
the derivatives and the nonzero endpoint slopes.
The vertical velocity is bounded above in modulus and
bounded away from zero in the outgoing direction
throughout this small neighborhood. This proves the
linear loss on the initial pieces. Beyond that
neighborhood the finite traces stay a fixed amount
below zero in real phase.

To apply this phase bound to the exact kernel, reflect the
private and pinched Gamma factors near zero. Their product
on \(z=jq\) becomes \(\Gamma(-z)\Gamma(1+z)=-\pi/\sin(\pi z)\).
On \(|\Im q|\ge Y_j\), its dominant sine exponential gives
the corresponding upper or lower phase, with an
\(O(j)\) prefactor in the \(dq\) variable. This follows
from the
[exact pinched kernel](MULTIVARIATE_ONE_PRIVATE_SIDE_RICE_EVENTUAL.md#the-pinching-left-contour)
and uniform Stirling for its remaining Gamma ratios on
the fixed local \(q\)-disk. Their arguments have real
parts comparable to \(j\) times their row exponents,
and the summed Stirling errors are \(O(j^{-1})\).
The linear phase loss bounds the integral on either
initial piece by \(O(\exp(-c_1 jY_j))\). The remaining
finite traces have exponential loss by the ordinary
compact Rice estimates. Increasing \(C\) gives a total
\(O(j^{-2})\) error.

On the retained central segment, the same exact pinched
kernel has \(T_j(z)=1+O((1+|z|)^2/j)\) on
\(|z|\le c+C\log j\), uniformly by the Stirling ratio
bounds in the regular comparison. The model
\(\Gamma(-z)\Gamma(1+z)(-e^{\xi_j})^z\) has
integrable exponential decay on the full line at every
fixed interior height. Extending the model integral
from the central segment to that line gives the actual
boundary contribution
\((1-e^{\xi_j})^{-1}+O(j^{-1})\), with the minus sign
of the upward left integral already included. We
have shown a local expansion of the form
\[
 \frac{M_j(t)}{M_j(0)}
     =(1-e^{\xi_j})^{-1}+O(j^{-1})
       +k(\xi)j^{1/2}J(\xi_j)e^{jD_j(\xi_j)}
                              (1+O(j^{-1})),              \tag{3}
\]
where \(D_j,J\) are independent of the order at the
given finite profile, \(D_j\to D\) with derivatives,
\(J\) is bounded away from zero and infinity, and
\(k(\xi)\) is the locally constant saddle integer.
The strict-sublevel contributions in the flowed arc
are exponentially small relative to both the boundary
and saddle scales. The estimate near zero uses only
the retained segment and finite traces.

At \(\Im D=0\) with \(\Re D<0\), choose the stopping
level between \(\Re D\) and zero. The same finite descent
takes the arc entirely below that level, even when an
endpoint approaches the saddle. The boundary term is
then nonzero and dominant. At \(\Im D=0\) with
\(\Re D>0\), neither endpoint can approach the saddle
from its initial value zero. The finite deformation
and (3) still apply locally. There is a single dominant
term, whether the saddle coefficient vanishes or
not. Since \(D=0\) is excluded, no zero can tend to a
field on this separating set.

## Connected sign regions and the two anchors

Write \(\xi=x+i\theta\). Since \(\partial_x\Im D=\Im q>0\),
each sign region at a fixed height is a horizontal
interval, possibly empty. As \(x\to+\infty\), the inverse
near \(q=1\) gives \(q\to1\) and \(\Im D\to\theta>0\).
These limits hold uniformly for compact height intervals
inside the strip. To see the phase limit, the stationary
value has the formula
\[
 D=-\log(1-q)-q-\tau q^2/2+
       \sum_{\rm retained}
            \bigl[a_i\log(1+d_iq/a_i)-d_iq\bigr],          \tag{4}
\]
where the regular sum converges locally near the real
gap endpoints after finite singular terms are removed.
The formula follows by integrating \(X\) and eliminating
\(\xi=-X(q)\). The nonsingular terms at \(q=1\) have
real limit and the inverse of the singular
\(\log(1-q)\) gives the stated imaginary part.
Consequently the positive sign region is connected,
by joining any two heights at a sufficiently large
positive real field.

If \(p_R<\infty\), the inverse near the nearest right
pole gives \(q\to p_R\) as \(x\to-\infty\).
If \(D_R\) is the sum of absolute gaps at that pole,
the stationary equation gives
\(\arg(p_R-q)\to(\theta-\pi)/D_R\), by the singular
logarithm \(-D_R\log(p_R-q)\). Formula (4) therefore
gives
\[
             \Im D\longrightarrow \pi-p_R(\pi-\theta).
\]
Indeed the coefficients \(a_i\) there sum to \(p_R D_R\),
and the imaginary part from \(-\log(1-q)\) tends to
\(\pi\). The local inverse and these limits are uniform
on compact interior height intervals. The negative
sign region projects precisely onto the interval
\(0<\theta<\pi(1-1/p_R)\), and is connected by joining
two heights at a sufficiently negative real field.

If \(p_R=\infty\), all retained shared gaps are
positive. At \(x\to-\infty\), \(q\) must escape,
since no finite pole has \(X\) diverging right in real
part. Its imaginary part stays bounded away from
zero, uniformly over compact interior height intervals.
Otherwise along a sequence with \(\Im q\to0\), the
escaping real part of \(q\) tends either to \(-\infty\)
or \(+\infty\). On the left, \(\arg(1-q)\to0\) and all
remaining contributions to \(\Im X\) are nonnegative,
which contradicts the strictly negative target height.
On the right, \(\arg(1-q)\to-\pi\), and
\[
 0\le\Im B(q)\le(\sigma_0-\tau)\Im q+\tau\Im q .
\]
This contradicts a target bounded away from \(-\pi\).
Hence \(\partial_x\Im D=\Im q\) gives
\(\Im D\to-\infty\) as \(x\to-\infty\), uniformly on
compact height intervals. The negative sign region is
again connected. This argument also covers the flat limit.

For a real field \(x<0\) sufficiently close to zero,
the analytic inverse at \(q=0\) gives \(q(x)<0\),
and \(D(x)>0\) by \(D(0)=0,D'=q\). At a field
\(x+i\theta\) with \(\theta>0\) sufficiently small,
\(\Im D<0\) and \(\Re D>0\). The known
[stable positive-real endpoint expansion](MULTIVARIATE_UNIFORM_ONE_PRIVATE_REGULAR_COMPLEMENT.md#the-positive-real-boundary-and-the-adjacent-comparison)
gives a bounded moment on a complex disk about \(x\).
Thus (3) forces the saddle integer to be zero at this
anchor, and throughout the connected negative sign
region. Near the negative-real parameter boundary,
choose a real field large enough that the gap saddle
has \(\Re D>0\). There \(\Im D\to\pi\) and the
[actual gap saddle expansion](MULTIVARIATE_UNIFORM_ONE_PRIVATE_REGULAR_COMPLEMENT.md#the-actual-right-saddle-on-compact-field-disks)
has nonzero coefficient. Comparing (3) in the neighboring
positive sign region gives a nonzero saddle integer
there, hence throughout that region.

## The active comparison

At any crossing of the imaginary \(q\) axis within the
fiber, write \(q=iy\), \(y>0\), and put
\(\ell_i=a_i/|d_i|>1\) for retained rows. The height
condition gives
\[
       \sum_{\rm retained}|d_i|\arctan(y/\ell_i)+\tau y
                       <\arctan y .
\]
For \(H(u)=u-\arctan u\), the function
\(H(u)/(u\arctan u)=1/\arctan u-1/u\) is strictly
increasing on \(u>0\). Its derivative is positive by
\(\arctan u>u/\sqrt{1+u^2}\), which follows directly
by differentiating the two sides from zero.
Thus the height inequality gives
\[
                   \sum_{\rm retained} a_i H(y/\ell_i)
                                      < H(y).
\]
At \(q=iy\), the diffuse quadratic in (4) is real.
The imaginary part in (4) is bounded above by the
left sum minus \(H(y)\), so \(\Im D<0\) at every
axis crossing. In the connected positive sign region
the saddle thus cannot cross the axis and stays on
the right, by the negative-real-gap anchor.

A limiting interior zero must lie where \(\Im D>0\)
and \(\Re D=0\), by (3) and the sign classification.
At each height the positive sign region is an interval,
and \(\partial_x\Re D=\Re q>0\) on it, so it contains
at most one such point. At a common root the two
limiting fields therefore agree. Their difference is
\(\Delta\log(1+1/n)\), so \(\Delta/n\to0\).
Now the two terms in (3) cancel at each order and
\(\Re q\) has a strictly positive limit. Dividing the
zero equations gives the
[linked order contradiction](MULTIVARIATE_UNIFORM_ONE_PRIVATE_REGULAR_COMPLEMENT.md#the-positive-real-boundary-and-the-adjacent-comparison)
for \(\Delta\ne0\). For \(\Delta=0\), division
forces \(\Re D_j=O(1/n)\), while either zero equation
requires \(n\Re D_j=-\tfrac12\log n+O(1)\).
No interior limit survives, which proves the theorem.

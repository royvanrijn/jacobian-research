# Effective localization in bounded subcritical thermal fields

Use the fully shared ordinary exponent profiles, curvature \(W,\rho\),
quartic center \(q_*,w_*\), and centered fields \(w_j\) of the
[moving quartic window](MULTIVARIATE_SHARED_MOVING_QUARTIC_WINDOW.md).
Let \(\delta=\max_{[0,1]}\rho-1<0\), \(\Delta=\sum_i(b_i-a_i)\)
and
\[
 a_2=\tfrac12F''(q_*),\qquad g=\tfrac1{24}F''''(q_*),\qquad
 r^2=\frac{a_2}{2g},\qquad c=-4gr^4,\qquad
 b_j=\frac{r(w_*-w_j)}c,\qquad M_j=jc .
\]
Both orders use the same profile. Thus \(\Im b_n=\Im b_N\) and
\[
              b_N-b_n=\frac{r}{c}\Delta\log(1+1/n),\qquad N=n+1.
              \tag{1}
\]

**Theorem.** Given positive rational \(W_-\le W_+\), \(B\ge1\)
and \(\epsilon>0\), one can compute positive rational
\(\delta_0\) and \(A\) such that the following holds. Suppose
\[
 W_-\le W\le W_+,\qquad -\delta_0\le\delta<0,\qquad
 n\delta^2\ge A,\qquad |b_n|,|b_N|\le B .
\]
At a consecutive common root, both scaled fields lie in the
same open \(\epsilon\)-disk centered at one of
\[
                          \pm iY_e,\qquad Y_e=\frac2{3\sqrt3}.
\]
There is no imposed bound on the integer \(\Delta\). The
constants are unevaluated.

## Effective finite integrals

The [effective quartic center estimates](MULTIVARIATE_EFFECTIVE_BOUNDED_SHARED_QUARTIC.md#effective-center-and-localization-on-the-real-lattice)
give computable constants on the prescribed \(W\)-band. In
particular \(q_*,1-q_*\) have positive lower bounds and \(-4g\)
has positive lower and upper bounds \(g_-,g_+\). The curvature
maximum \(q_m\) satisfies
\(|q_*-q_m|\le C|\delta|\) and \(\rho'(q_m)=0\).
The uniform bound on \(\rho''\) therefore gives
\[
     |\rho(q_*)-\rho(q_m)|\le C'|\delta|^2 .
\]
For sufficiently small negative \(\delta\) this implies
\(a_2<0\) and two-sided comparability of \(-a_2\) and \(|\delta|\),
using \(2a_2=(\rho(q_*)-1)/(q_*(1-q_*))\).
In particular there are computable positive constants with
\(r^2\asymp|\delta|\) and \(c\asymp\delta^2\).
It suffices to give computable upper bounds on \(r\) and lower
bounds on \(M_n\).

The uniform finite-integral and adjacent-amplitude estimates also
apply in the [escaping-field proof](MULTIVARIATE_EFFECTIVE_SHARED_ESCAPING_QUARTIC_FIELDS.md#effective-finite-contour-near-the-center).
Here are the data needed for the bounded scaled charts.
Choose a fixed small computable physical half-width \(a>0\)
around \(q_*\). On a larger disk the exact phase difference
has the expansion
\[
 \Phi_w(q_*+h)-\Phi_w(q_*)=(w_*-w)h+a_2h^2+gh^4+h^5E(h)
\]
with a computable bound on \(E\) and its Cauchy derivatives.
Choose \(a\) so the last term and its real derivative at
\(\pm a\) are small compared with the quartic and its
derivative. For \(|b|\le B\), \(|w_*-w|\le g_+Br^3\).
The real phase at \(q_*\pm a\) consequently has a fixed
strict loss relative to its value at \(q_*\), and has
strict opposite derivative signs, once \(r\) is small.
As \(F''=(\rho-1)/(q(1-q))<0\), the loss extends to the
whole exterior real lattice. The same bound works on
small neighborhoods of the two endpoints.

For completeness, the amplitude bounds needed here include a
comparison of the actual amplitudes at adjacent orders.
On a fixed sufficiently small complex rectangle about the
physical interval, all \(a_i+d_iq\) have real parts at least
a positive fixed multiple of \(m_i\). Their arguments stay
in a strict sector, as do \(q,1-q\). The same is true of
the rescaled Gamma arguments when the order lies in a fixed
relative complex disk about any sufficiently large positive
real order. Apply Binet's logarithmic Stirling remainder on
these disks. In its positive-real-part integral, the kernel
\[
 \frac{1}{t}\left(\frac1{e^t-1}-\frac1t+\frac12\right)
\]
is bounded by \(1/2\) for \(t>0\). For \(t\ge1\) this follows
from \(e^t\ge1+t\). For \(0<t\le1\), write the parentheses
as \((\coth(t/2)-2/t)/2\) and use
\(x\cosh x-\sinh x=\int_0^x v\sinh v\,dv\).
The resulting remainder is bounded by \(1/(2\Re z)\) for
an argument \(z\) in this half-plane. Summing the row
remainders uses \(\sum_i1/m_i\le W_+\), since each nonzero
\(d_i\) is an integer. Thus the logarithmic remainder is
\(O(j^{-1})\) with a computable constant on each relative
order disk. Cauchy's estimate gives \(O(j^{-2})\) for its
real-order derivative. The leading amplitude and its reciprocal
are bounded uniformly: the relative row logarithms are bounded
by a computable multiple of \(|d_i|/m_i\), whose sum is at most
\(W_+\). These arguments give effective \(q\) derivative bounds
on smaller rectangles as well.

Use the [half-grid cotangent summation](MULTIVARIATE_VARYING_PROFILE_SUBCRITICAL_SADDLE.md#the-lattice-sum-and-its-contour)
on the physical interval with endpoints displaced by at most
\(1/j\) to neighboring half-grid points. The entropy inequality
for the binomial coefficient and the real Stirling bound give
an exponential bound on the omitted lattice summands. The
normalized Stirling amplitudes just described are analytic
where the rectangle is used. Choose a fixed small height
\(\nu>0\) so a vertical displacement changes the real
phase on the walls by less than half the physical gap.
The real increase on either horizontal is bounded by
\(C\nu^2+|\Im w|\nu\) relative to its real-axis value.
The real phase everywhere on the retained interval is at
most its centered value plus \(a|w_*-w|\), as the quadratic
and quartic terms are negative and dominate the remainder.
Choose \(\nu\) and then decrease the upper bound on \(r\)
so the horizontal increase plus this last error is less
than \(\pi\nu\). Thus the nonconstant parts of the
cotangent on the horizontal sides are exponentially small
after normalization. The two constant parts deform to the
finite real interval and add. Moving both endpoints back
to \(q_*\pm a\) also costs an exponential error. All these
steps give the finite-integral formula of the
[moving quartic proof](MULTIVARIATE_SHARED_MOVING_QUARTIC_WINDOW.md#uniform-centers-and-finite-integral-reduction)
with a computable gap uniform on the present fields.

Write the result after division by the nonvanishing centered
factor \(j^{1/2}r e^{j\Phi_{w_j}(q_*)}\) as
\[
 \int_{-a/r}^{a/r} A_j(q_*+rz)e^{M_j\Psi_{r,b_j}(z)}\,dz+\mathcal E_j,
 \qquad
 |\mathcal E_j|\le C r^{-1}e^{-\kappa j},
 \tag{2}
\]
where fixed powers of \(j\) have been absorbed in a smaller
positive computable \(\kappa\). On the physical rectangle,
\[
 \left|\frac{A_N(q)}{A_n(q)}-1\right|\le Cn^{-2}.
 \tag{3}
\]
In (2) the exact scaled phase is
\[
 \Psi_{r,b}(z)=G_b(z)+\frac{r}{-4g}z^5E(rz),\qquad
 G_b(z)=-z^4/4-z^2/2+bz .
 \tag{4}
\]
On every fixed scaled rectangle its difference from \(G_b\),
including the required derivatives, is \(O(r)\) with a
computable constant. On horizontal lines of bounded heights
the far parts admit \(\Re\Psi_{r,b}(x+iy)\le-\kappa_1 x^4\)
for \(|x|\ge R\), up to the physical endpoints, with
computable \(R,\kappa_1>0\). To see this, shrink \(a\) so
\(|rz|\) on the contours makes the quintic correction a
small fraction of \(|z|^4\). Then increase \(R\) to absorb
the real contributions from the bounded heights and the
linear and quadratic terms. There is also a bound
\(|\Psi_{r,b}(x+iy)|\le C(1+|x|)^4\) on the same paths.

## Computing a cover of the nonfold field pairs

Use the model contours of the
[regular thermal proof](MULTIVARIATE_SHARED_SUBCRITICAL_THERMAL_REGULAR_FIELDS.md#the-regular-thermal-contours).
They give the following two classes of local charts.

For every \(b\) off the two imaginary rays starting at the
folds, there is a model simple saddle with an effective
localizing witness on a horizontal line. When \(\Re b>0\)
the positive real saddle branch has strict global maximum
on that horizontal, as follows from the increasing equation
for the imaginary height in the regular thermal proof.
Reflect for \(\Re b<0\). On the imaginary axis between the
folds the branch through \(iv\), with \(|v|<1/\sqrt3\), has
the same strict real quadratic loss. A one-saddle witness
certifies such a field cell and the strict losses on its
contour.

At \(b=iY\) with \(Y>Y_e\) there are two simple model
saddles \(z_R=u+iv\), \(z_L=-u+iv\), where
\(u^2=3v^2-1\) and \(Y=8v^3-2v\). Their model horizontal
can be divided at \(x=0\), with a gap \(u^4/4\) between
that midpoint and each saddle. The two-saddle chart
uses each half-line moved to its own actual saddle
height. Join both halves at a common point \(iv_0\)
on the seam by short vertical segments along \(x=0\).
The seam's strict gap persists. At the physical ends
connect to the original endpoints by vertical segments.
Their strict far loss follows from (4).
The lower imaginary ray has the conjugate chart.

The stationary action at a model saddle \(z\) equals
\(3z^4/4+z^2/2\). At \(iY\) the right-minus-left action
has imaginary part \(4u^3v>0\); the real separation of
the saddles is \(2u>0\). In the upper chart shrink the
field neighborhood so both strictly positive margins
persist. In the lower chart use the absolute imaginary
margin with the opposite sign.

The witnesses can be verified by finite tests. Enumerate
rational field rectangles, rational saddle centers and
isolation disks, and rational widths for the local
quadratic neighborhoods. The simple saddle disks are
checked by linear Rouché applied to \(G'_b\). At center
\(z_0\) and radius \(t\), compare
\[
 |G'_b(z_0)|+3|z_0|t^2+t^3
 \quad\text{against}\quad |G''_b(z_0)|t .
\]
A strictly positive margin on the right accommodates the
\(O(r)\) derivative error in (4) for a computable upper bound
on \(r\).
For the contour inequalities, bound the real second
derivative near each saddle and check the strict loss
on the bounded pieces outside those neighborhoods by
interval subdivision. Use the explicit far tail for
the remaining horizontal parts, and the same bound
on the far vertical connectors. In a two-saddle chart
check the two action/separation margins using the
enclosed saddles. The model inequalities above are
strict, so rational neighborhoods and these tests
produce a witness at every indicated model point.
Use convex field neighborhoods for the two-saddle
charts. Their construction also gives effective \(r\)
bounds and constants uniform on each verified cell.

Consider the compact box of triples
\((\Re b_n,\Re b_N,\Im b_n=\Im b_N)\) in \([-B,B]^3\).
Cover it by one-saddle witnesses in either coordinate,
common convex two-saddle charts in both coordinates,
and the two sets in which both fields belong to the
same open \(\epsilon\)-disk at a fold. These sets give
an open cover. Indeed, if neither field has a
one-saddle witness, both are on the imaginary rays.
Their equal imaginary parts then make them the same
model point. A strictly superfold point has a
two-saddle witness, and a fold point is covered by
the disk alternative. Enumerate finite witness lists
in parallel with rational meshes of the triple box.
At each stage check whether every closed mesh cell
lies in one member of the proposed cover. Inclusion
in the fold alternatives is checked by interval
distance to the two computable algebraic centers. The closed
mesh cells assigned to this alternative have a strict margin
inside the respective \(\epsilon\)-disks for both fields.
Compactness and strict inequalities give termination.
This construction uses only the fixed polynomial
model and uniform bounds on its perturbation.

## Effective comparison on a common two-saddle chart

On each one-saddle chart the ordinary Gaussian estimate
gives a single nonzero leading term, so (2) cannot vanish
for sufficiently large \(M_j\) and small \(r\).
Here and in the two-saddle charts the estimates are
effective by the following local calculation. At the
actual saddle \(z_s\), on the symmetric horizontal
subsegment the real loss is at least \(a_1|x|^2\).
Taylor's bound for the quadratic approximation is
\(C|x|^3\); the pure quadratic has the same strict
real loss after decreasing \(a_1\).
Interpolate the exponential between the actual and
quadratic phases. The error is bounded by a constant
times \(M_j |x|^3 e^{-M_j a_1|x|^2}\).
The amplitude error uses the first moment under the
same Gaussian. The explicit Gaussian integrals and
the strict exterior losses give packet integrals
\(I_{j,s}\) with leading term
\[
 e^{M_j\Psi_{r,b_j}(z_s)}
 A_j(q_*+rz_s)\sqrt{\frac{2\pi}{-M_j\Psi''_{r,b_j}(z_s)}} .
\]
The square root uses the half-plane \(\Re(-\Psi'')>0\).
The relative error is bounded by a computable
constant times \(M_j^{-1/2}\). Each packet also
admits the corresponding absolute Gaussian moment
bounds. The uniform phase and amplitude bounds give
a positive lower bound for the leading modulus after
its stationary exponential and \(M_j^{-1/2}\) are removed.

Let \(D(b)=\Psi_{r,b}(z_R)-\Psi_{r,b}(z_L)\) on
a common convex chart. At a root, the two packets in
(2) have comparable moduli. In fact the exterior
error divided by either leading packet is bounded by
\[
 C r^{-1}\sqrt{M_j}\exp(-\kappa j+C'M_j).
 \tag{5}
\]
Shrink the bound on \(r\) so \(C'c<\kappa/2\).
For \(M_j\ge1\), the inequalities \(g_-r^4\le c\le
g_+r^4\) give \(r^{-1}\le(g_+j)^{1/4}\).
Consequently (5) is bounded by a computable multiple
of \(j^{-2}\), for example by a power-term lower bound
for the positive exponential. Gaussian comparison
and (2) then give \(|\Re D(b_j)|\le C_1/M_j\).

Stationarity gives \(D'(b)=z_R-z_L\). Its real part
has a constant sign and magnitude at least \(d>0\)
on the chart by the verified separation. Integrating
along the real segment in (1) yields
\[
 |b_N-b_n|\le \frac{2C_1}{dM_n},\qquad r|\Delta|\le \frac{4C_1}{d}.
 \tag{6}
\]
This uses \(n\log(1+1/n)\ge1/2\).
The estimate holds without an a priori bound on \(\Delta\).

Now take both upper-order packets along the lower-order
contours, which also give the upper integral in (2)
by analyticity and the same endpoints. Put \(M=M_n\),
\(\psi=\Psi_{r,b_n}\), and \(\omega=r\Delta N\log(1+1/n)\).
The analytic multiplier of the scaled integrand is
\((A_N/A_n)e^{c\psi(z)+\omega z}\), by (1). For each
lower saddle \(z_s\), put \(\Lambda_s=e^{c\psi(z_s)+\omega z_s}\).
Using (3), (6) and stationarity gives
\[
 I_{N,s}=I_{n,s}(\Lambda_s+\varepsilon_s),\qquad
 |\varepsilon_s|\le C_2\bigl(c/M+|\omega|/\sqrt M+n^{-2}\bigr).
 \tag{7}
\]
Indeed the local exponent difference from the stationary
value is bounded by \(C c|z-z_s|^2+|\omega|\,|z-z_s|\).
Integrate the quadratic and linear moments against the
absolute Gaussian bound. On bounded nonlocal portions
use the strict gap together with the factor
\(O(c+|\omega|+n^{-2})\) from the multiplier difference.
On an unbounded scaled horizontal tail use
\[
 |e^\zeta-1|\le|\zeta|e^{|\zeta|},\qquad
 |\zeta|\le Cc(1+|x|)^4+C|\omega|(1+|x|).
\]
The quartic majorant absorbs the exponential on the right
for a computable \(M\), since \(|\omega|\) is bounded by
(6) and \(c\) is small. The remaining polynomial moments
give a computable bound of the form
\(C(c+|\omega|+n^{-2})\sqrt M e^{-\kappa_2 M}\)
after division by the lower packet. The far vertical
segments at the physical endpoints obey the same
bound, using their bounded scaled lengths and their
quartic loss. Thus the nonlocal comparison retains the
small factors in (7). The discrete outer errors are handled
separately by (5), which is also valid for \(\mathcal E_N\)
relative to a lower-order packet.

For \(\Delta\ne0\), one has \(|\omega|\asymp r|\Delta|\).
The stationary exponent difference has real part
\[
    c\,\Re D(b_n)+\omega\Re(z_R-z_L).
\]
Its second term has modulus at least \(d r|\Delta|\).
The first has modulus at most \(C_1/n\), which is smaller
by an arbitrarily large effective factor when \(M\) is
large, since \(r|\Delta|\ge r\) and \(1/n=c/M\le g_+r^4/M\).
Both stationary exponents are bounded by constants from
the chart and (6). Therefore
\(|\Lambda_R-\Lambda_L|\ge d_1 r|\Delta|\) for a
computable \(d_1>0\), using the difference of moduli.

For \(\Delta=0\), the verified imaginary action margin
gives \(|\Im(cD(b_n))|\ge d_2 c\).
Shrink the upper bound on \(c\) so the bounded action
has \(|\Im(cD)|\le1\). Exponentiation then gives
\(|\Lambda_R-\Lambda_L|\ge d_3 c\) with a computable
\(d_3>0\).

At common cancellation, subtract \(\Lambda_L\) times the
lower-order equation (2) from its upper-order counterpart.
The two lower packets are comparable in modulus. Equations
(5) and (7) consequently bound the multiplier difference
by \(C_3(c/M+|\omega|/\sqrt M+n^{-2})\).
For \(\Delta=0\), its ratio to \(c\) tends to zero with
an effective bound in \(M\). For \(\Delta\ne0\), its
ratio to \(r|\Delta|\) also tends to zero effectively,
using \(r|\Delta|\ge r\), \(c\le g_+r^4\), and (6).
This contradicts the respective positive lower bounds.

Take the minimum of the finitely many permitted \(r\)
bounds and the maximum of the finitely many \(M\)
thresholds from the verified cover. The computable
comparisons \(r^2\asymp|\delta|\), \(c\asymp\delta^2\)
give \(\delta_0,A\) in the statement. \(\square\)

# Two isolated ends determine a conical exterior

Let \(A\) be a stationary planar Gaussian three-cell partition at
\(0<r<1\), so each actual pair interface has a constant pair-score
difference. Suppose \(A\) agrees outside a compact set with a
partition \(C\) conical about the Gaussian origin. If two distinct
pairs have bounded isolated graph ends in full forward half-planes,
then \(C\) is the centered regular simplex. The claim does not
require connected cells or an assumption on the third interface.

## The line equations give two reflections of the cone

For one of the two ends, use orthonormal coordinates \((s,v)\)
about the origin, with \(s\) increasing in its outward direction.
The cells above and below the bounded graph \(v=h(s)\) are the
two incident labels. For all large positive \(s\), its whole
transverse slice is outside the exceptional disk. The slice of the
cone at parameter \(s\) is the dilation by \(s\) of its slice at
parameter one. Consequently \(h(s)/s\) is constant for almost
every large \(s\). Boundedness forces that constant to vanish,
and continuity gives \(h=0\) on a tail. This argument applies to
both ends. On either straight ray the pair score tends to zero
at infinity, from the centered binary half-plane profile and
the Gaussian-small probability of a source outside that end's
forward half-plane. The stationary constant is therefore zero.

Fix either pair \(i,j\). Its score contrast is real analytic and
vanishes along its whole supporting line by continuation. Put
\(\sigma^2=1-r^2\). Injectivity of the longitudinal Gaussian
convolution yields, for almost every \(s\),
\[
 \int (1_{A_i}-1_{A_j})(s,v)e^{-v^2/(2\sigma^2)}\,dv=0.
                                                               \tag{CX1}
\]
For example injectivity follows by viewing the transverse average
as a bounded tempered distribution and using the nowhere-zero
Gaussian Fourier multiplier. For \(|s|\) beyond the exceptional
disk, write \(v=|s|u\) and \(\eta=\operatorname{sgn}(s)\) in (CX1).
Conicality gives
\[
    \int (1_{C_i}-1_{C_j})(\eta,u)e^{-s^2u^2/(2\sigma^2)}\,du=0.
\]
For either \(\eta\), the transform is analytic when the coefficient
of \(u^2\) has positive real part. Thus it vanishes for every
positive coefficient. After replacing \(u\) by \(\sqrt t\) on the
positive half-line, injectivity of the Laplace transform makes
the even part of that pair contrast zero. To apply injectivity,
one may first damp the resulting density by \(e^{-t}\) to obtain
a finite measure and then use uniqueness of its Fourier transform.
The contrast takes distinct values on the three labels. Its oddness
means that reflection in the line exchanges \(i,j\) and preserves
the third label in \(C\).
[Related longitudinal slice identity (PG10)](THREE_CELL_PRESSURE_GRADIENT_TAIL.md#exterior-geometry-when-all-three-ends-straighten)

## Compatibility of the two angular reflections

Let \(\theta\) be the smaller angle between the outward directions.
It is at least \(\pi/2\), since an acute angle would put one end
deep inside the other's forward half-plane, where one of its
two labels is absent. Name the common label \(0\), and the other
labels \(+\) and \(-\). In angular coordinates from the forward
bisector the two outward rays have angles \(a\) and \(-a\), where
\(a=\theta/2\). The reflections just obtained are accompanied
by the transpositions \((0,+)\) and \((0,-)\). Their composition
rotates by \(2\theta\) and cycles the three labels. The cases
\(\theta=\pi/2\) and \(\theta=\pi\) are impossible. In those cases
the rotation has order two or one, respectively, and a corresponding
power of the label cycle still has no fixed label.

The forward angular sector \((-a,a)\) belongs to label \(0\).
Indeed each point there lies in at least one defining forward
half-plane, on its common side as determined by the overlap at
the forward bisector. The adjacent arcs \((a,a+\pi/2)\) and
\((-a-\pi/2,-a)\) belong to \(+\) and \(-\), respectively.
Thus a common-label interval can lie only in the forward sector
of width \(\theta\), or in the remaining backward arc of width
\(\pi-\theta\). The third power of the composite rotation preserves
label \(0\), so the forward sector rotated by \(6\theta\) must fit
in these two arcs up to angular null sets. The rotated interval
cannot cross an open noncommon arc. It cannot fit in the backward
arc because \(\theta>\pi-\theta\). Hence it is the original sector
up to null sets, forcing \(6\theta\) to be a multiple of \(2\pi\).
The only such \(\theta\) in \((\pi/2,\pi)\) is \(2\pi/3\).
At that angle the label cycle rotates the common sector into the
two other sectors of width \(2\pi/3\), covering the circle. Thus
\(C\) is the centered regular simplex.

# Effective angular separation for one private side

Consider ordinary two-monomial profiles with exactly one side having
private rows. The exponent of a private row is positive on that side
and zero on the other. Write \(P\) for the total private degree and
\(S\) for the sum of squared unequal row gaps divided by the maximum
exponent in the row, summed over all unequal rows.

**Theorem.** There is a computable absolute terminal cutoff above which
the two consecutive factorial moment polynomials have disjoint root
arguments, even after independent positive rescalings.

The [high-curvature angular proof](MULTIVARIATE_EFFECTIVE_ABSOLUTE_CURVATURE_BOUND.md)
covers \(S\ge64\). With \(S\le64\), the
[regular-band angular theorem](MULTIVARIATE_EFFECTIVE_ONE_PRIVATE_REGULAR_ANGLES.md)
covers \(P\ge2\), and \(P=1\) when the shared curvature
\(\sigma=\sum_{\rm shared}(b_i-a_i)^2/a_i\) stays at least a fixed
positive amount above one. Orient the private row as \((0,1)\)
in the remaining band. Use the centered logarithms \(\xi_j\)
of the [effective one-private comparison](MULTIVARIATE_EFFECTIVE_ONE_PRIVATE_INTERIOR.md).
Their imaginary parts agree at roots of equal argument. Unlike the
shift for a common physical root, their real difference can be arbitrary.

## Interior heights and negative real arguments

For \(P=1,\ 0\le\sigma\le5/4\) and any fixed closed interior height
strip, use the paired equal-height compact cover in the
[effective interior proof](MULTIVARIATE_EFFECTIVE_ONE_PRIVATE_INTERIOR.md#an-effective-cover-of-the-pairs-of-fields).
Its one-term witnesses already exclude individual roots. In a two-term
witness, the two limiting fields coincide and \(\Re q=\Re D'>0\),
where \(D\) is the saddle value relative to the endpoint. The modulus
comparison of the regular angular theorem gives real displacement
\(O((1+\log n)/n^2)\) and division of the zero equations forces
\(|e^{D(\xi_n)}-1|=O((1+\log n)/n)\).

Here every interior balance has \(0<\Im D<\pi\). For \(\sigma>1\)
this is the regular-band phase argument. At the critical limit
\(\sigma=1\), that argument uses \(x_0=0\): the right inverse tends
to \(q=0\) on the upper approach to \(\xi=0\), and the endpoint
phase there is zero. Thus \(\Re D(iy)<0\) for \(0<y<\pi\) by
\(\partial_y\Re D=-\Im q<0\). Any balance with height \(y>0\)
then has \(x=\Re \xi>0\). For \(x>0\) the real entropy phase has
a unique stable positive maximizer, so the upper sheet continues
to its real saddle at \(y=0\), where \(D\) is real. The relation
\(\partial_y\Im D=\Re q>0\) gives the strict lower bound, and
the right-gap value \(\Im D=\pi\) gives the strict upper bound.

For \(\sigma<1\), the
[activation theorem](MULTIVARIATE_UNIFORM_SUBCRITICAL_PRIVATE_ACTIVATION.md#the-active-comparison)
gives \(\Im D>0\) and \(\Re q>0\) on the active balance set.
For fixed real \(x\), once \(\Im D(x+iy)>0\), it cannot return to
zero on increasing \(y\): its derivative \(\Re q\) is positive
throughout the positive sign region by that same theorem.
It therefore stays strictly below its right-gap limit \(\pi\).
Thus the necessary interval margin \(|e^D-1|>0\) exists at
every balance in the closed interior strip. Adding that strict
test to the existing finite search certifies angular exclusion.

Near the negative real ray the
[exact right-line kernel](MULTIVARIATE_EFFECTIVE_ONE_PRIVATE_NEGATIVE_BOUNDARY.md#the-right-line-and-its-uniform-error)
has the factor \((-1)^{j+1}\). The residue orientation is fixed
across orders. Thus its two-packet chart uses the stationary
value \(D(\xi-i\pi)+i\pi\) in (1) of the regular angular theorem,
with the possible remaining constant sign in the amplitude.
The shifted value \(D(\zeta)\) is real at real \(\zeta\) and has
bounded derivatives throughout the uniform strip. Its upper
root-angle convention goes through height \(\pi\). Shrinking
the width, we have \(|\Im D(\xi-i\pi)|\le\pi/2\), so
\(|e^{D(\xi-i\pi)+i\pi}-1|\ge1\). The positive real slope
\(\Re D'>1\) in the
[negative-boundary comparison](MULTIVARIATE_EFFECTIVE_ONE_PRIVATE_NEGATIVE_BOUNDARY.md#consecutive-modulus-separation)
and the regular angular modulus argument now give a computable
threshold independent of the degree difference.

## Bounded quadratic fields

Use the notation of the
[effective bounded quadratic proof](MULTIVARIATE_EFFECTIVE_PRIVATE_QUADRATIC_WINDOW.md).
In particular \(H_j=\sqrt{j/\mu}\,\xi_j\), \(\mu=1-\sigma>0\),
and \(\alpha=(n\mu^3)^{-1/2}\). For \(|H_j|\le B\), the analytic
real-order interpolation \(J_j(H)\) has integer-order residual
\(O(j^{-2})\) and satisfies
\[
 |J_j-G|+|(J_j-G)'|+|j\partial_j J_j|=O(\alpha),
 \qquad G(H)=\int_0^\infty e^{-u^2/2+Hu}\,du .
\]
These bounds hold with computable constants on any fixed enlarged
field cap. At an equal-angle pair the field relation is
\[
     H_N=e^{\lambda/2}H_n+v,\qquad \lambda=\log(N/n),
                  \qquad v\in\mathbb R,\quad |v|\le2B+1.
\]
Let \(B_G(H,v)=\int_0^1G'(H+t v)\,dt\). The compact residual
\(|G|+|B_G|\) has positive infimum by the simplicity and same-height
uniqueness in the cited half-Gaussian proof. This infimum can
be bounded below by interval quadrature as in that source.
At a proposed root \(|G(H_n)|=O(\alpha)\), and at the pair
\[
                    |v B_G(H_n,v)|=O(\alpha |v|+\lambda).
\]
To obtain this, compare the \(J_n\) difference from \(H_n\)
to \(H_n+v\) with \(vB_G\) by the derivative bound, the
real-order change at fixed field by \(|j\partial_jJ_j|\lambda\),
and the remaining field change by the bounded derivative times
\(B(e^{\lambda/2}-1)\). The exact lattice errors contribute
\(O(n^{-2})\). For sufficiently small \(\alpha\) absorb the
\(\alpha|v|\) error using the residual lower bound. This gives
\(|v|\le K_0\lambda\) with computable \(K_0\).

Taylor expansion of \(G\), the same interpolation bounds, and the
lattice error now give
\[
 |H_n G'(H_n)+2(v/\lambda)G'(H_n)|=O(\alpha+\lambda).
\]
This follows by dividing by \(\lambda/2\), using the model
increment \(H_n\lambda/2+v+O(\lambda^2)\); the error in the
field-derivative comparison is \(O(\alpha\lambda)\).
Multiplication by \(\overline{G'(H_n)}\) and taking the imaginary
part remove the \(v\) term. The compact residual
\[
        |G(H)|+\left|\Im H\,|G'(H)|^2\right|
\]
has positive infimum, since \(G\) is positive on the real axis
and \(G'=1+HG\). Its lower bound is computed by interval
quadrature. The bound \(n\mu^3\ge K(B)\) with a computable
\(K(B)\) makes \(\alpha+\lambda\) too small for that residual.
This proves angular exclusion in the bounded quadratic chart.

## Linked estimates on escaping endpoint scales

For the
[cubic](MULTIVARIATE_EFFECTIVE_PRIVATE_CUBIC_FIELD.md#linked-endpoint-and-saddle-ratios),
[quadratic](MULTIVARIATE_EFFECTIVE_PRIVATE_QUADRATIC_CORRIDOR.md#consecutive-comparison-and-the-mixed-scales),
and [thermal](MULTIVARIATE_EFFECTIVE_PRIVATE_THERMAL_ENDPOINT.md#consecutive-thermal-comparison)
endpoint comparisons, use each source's common lower-order scales
\(r,c,M=nc\). They have
\(c=\chi r^\beta\), \(0<\chi\le C\), with \(\beta=2\) or \(3\).
Their assumptions make \(r\) small and \(M\) large.
The saddle value is \(D(b)=\Psi_b(z_s(b))\) in their scaled
fields \(b,b'\), with \(D'=z_s\). Equal physical arguments give
\[
          d=M'(b'-b)\in\mathbb R,\qquad M'=Nc .
\]
The single-order modulus estimates and initial comparison in
those proofs give \(|d|\le1\) for sufficiently small \(r\)
and large \(M\). Their derivation only uses that \(b'-b\)
is real. The same is true of the exact multiplier
\(e^{c\Psi_b(z)+dz}\) on the common lower-order contours.

Equations (7)–(8) of the cubic proof compare the exact endpoint
and saddle packets. At an equal-angle pair with both moments zero,
the ratios agree up to \(O((M/r)e^{-\gamma n})\).
Taking the logs near their leading values gives
\[
 cD(b)+dz_s(b)=O\!\left(
     \frac{c+|d|}{M}+\frac{|d|}{\sqrt M}+n^{-2}
       +\frac{c}{rn}+\frac{c|d|}{r^2 n}
       +\frac Mr e^{-\gamma n}\right).                     \tag{1}
\]
Indeed both ratios are bounded above and away from zero when
\(|d|\le1\), and the chosen leading exponent can be made small
by the same preliminary modulus comparison (its bound for
\(|d|\) tends to zero as \(M\) increases and \(r\) decreases).
The principal logarithms then give (1). The endpoint and saddle
remainders in (1) are precisely the remainder sizes of the linked
contour estimates, including the \(B_j\) difference for the endpoint.

Since \(\Re z_s\ge\kappa>0\), take real parts in (1).
At a root the modulus equation gives
\(|c\Re D|\le C(1+\log M)/n\). Also \(c/r\) and \(c/r^2\)
are bounded in each of the scales. Absorb all terms multiplying
\(|d|\) on the right for sufficiently large \(M,n\).
The physical exponential remainder is \(O(1/n)\) because
\(M/r=nc/r\le Cn\). This gives
\[
                         |d|\le C(1+\log M)c/M.           \tag{2}
\]
Substituting (2) back into the right of (1) gives \(o(c)\)
uniformly as \(M\to\infty\) with computable constants. For
the exponential term divided by \(c\), note that
\(1/r\le(Cn/M)^{1/\beta}\), so \(n e^{-\gamma n}/r\)
is bounded by \(C'n^{3/2}e^{-\gamma n}\).
The other quotients tend to zero directly from \(c/M=1/n\),
the bounds on \(c/r,c/r^2\), and (2).

On the escaping cubic scale the sector-centered pure height
\(D_0=(2/3)b^{3/2}\) has imaginary part bounded away from
zero, with center \(b=e^{i\pi/3}\). On the escaping quadratic
scale the pure height is \(D_0=b^2/2\), with center
\(b=e^{i\pi/4}\). In both sectors the admitted fixed small
width and perturbation in the source proofs, reduced if needed,
give \(\Im D\ge\eta>0\) with effective \(\eta\).
Taking imaginary parts in (1) contradicts this since
\(|d z_s|=o(c)\) by (2).
The mixed scales in the quadratic corridor use single-root
angular confinement and equality of physical imaginary fields,
so they reduce to the bounded quadratic result above as before.

On the thermal equality set away from the real axis the
same positive margin holds. In the notation of its
[pure equality set](MULTIVARIATE_EFFECTIVE_PRIVATE_THERMAL_ENDPOINT.md#the-pure-thermal-equality-set),
write \(z_0=x+iy\), \(a\in\{-1,1\}\). At \(\Re D_0=0\),
\[
   y^2=\frac{x^2(3a+4x)}{3a+12x},\qquad
   \Im D_0=y\left(x(2x+a)-\frac23 y^2\right).
\]
For \(a=1\) the parenthesis is positive since \(x>0\) and
\(y^2\le x^2\). For \(a=-1\), the positive denominator
and \(y^2\ge0\) imply \(x\ge3/4\) and \(y^2\le x^2/3\);
the parenthesis is at least \(x(16x/9-1)>0\).
Thus all thermal balance charts away from real coexistence
have a positive computable imaginary height margin after
shrinking the chart and thermal perturbation bound.
They are excluded by (1)–(2).

## Thermal real coexistence

For \(a=-1\), the only real thermal equality point is
\(b=-3/16\). Take a symmetric disk about that point, with
a smaller concentric disk for comparison. In the representation
(3) of the thermal proof, split the original real integral at
scaled coordinate \(1/8\). The endpoint packet consists of the
integral up to \(1/8\), the left-wall correction, and the
holomorphic physical remainder. The second packet is the rest
of the original real integral. Both are holomorphic in the field
and real for real fields. Indeed \(A_j\) and \(\Phi\) are real
on real variables, the integral defining \(B_j\) pairs conjugates
there, and the physical remainder is the exact normalized
polynomial minus these real terms.

On the larger disk the endpoint leg has uniform linear loss.
Deform the rest of the integral through the continued saddle,
retaining the vertical connector at \(1/8\) in that packet.
The [thermal coexistence contour](MULTIVARIATE_EFFECTIVE_PRIVATE_THERMAL_ENDPOINT.md#effective-contour-neighborhoods)
and [linked ratio bounds](MULTIVARIATE_EFFECTIVE_PRIVATE_THERMAL_ENDPOINT.md#consecutive-thermal-comparison)
hold uniformly on this disk. They also hold below the real
axis by conjugation. For \(|d|\le1\) large \(M\) keeps both fields
needed in the comparison in a larger permitted contour neighborhood.
For each fixed real \(d\), the two exact packet ratios are analytic
functions of the lower field \(b\), using upper field \(b+d/M'\).
They are nonzero by their leading term estimates. On the original
real line the packet ratios are positive, by the real leading
amplitudes and the small relative errors. Define
\[
 \mathcal R(b,d)=
   \log(R_N/R_n)-\log(E_N/E_n)-cD(b)-dz_s(b),
\]
using their real branches. This function is holomorphic in \(b\)
and real when \(b\) is real. The linked estimates bound
\(|\mathcal R|\) by a constant times the right of (1) on the
larger disk, without a zero assumption. The estimates are those
for the same two fixed real packets: the real endpoint segment
does not depend on the field, while the saddle contour and its
connector are deformations of the second real segment.
The exponentially small exact remainder included in the first
packet costs \(O((M/r)e^{-\gamma n})\), as in (1).

At the proposed pair, cancellation gives equality of the exact
ratios. Their real log comparison first gives (2). Now at this
fixed real \(d\), the bound on \(|\mathcal R|\) is \(o(c)\)
throughout the larger disk, by the preceding scale checks.
Cauchy bounds its field derivative by \(o(c)\) on the smaller disk.
Since \(\mathcal R\) is real on the diameter,
\[
                   |\Im\mathcal R(b,d)|\le \Im b\,o(c)
                   \qquad(\Im b>0).
\]
The holomorphic saddle has \(D'=z_s\) with \(\Re z_s\ge\kappa>0\)
on the disk and real stationary value on its diameter.
Thus \(\Im D(b)\ge\kappa\Im b\) and \(|\Im z_s(b)|\le C\Im b\).
At a root pair the log ratio difference is zero, while
\[
   \Im\bigl(cD(b)+dz_s(b)+\mathcal R(b,d)\bigr)
                \ge \Im b\,(\kappa c-C|d|-o(c))>0
\]
for large \(M\) by (2). This contradiction excludes the real
coexistence neighborhood as well. The use of holomorphic,
real-on-real packet errors is needed here because \(\Im b\)
can approach zero with \(M\).

## Combining the uniform thresholds

All the scale reductions of the
[thermal endpoint proof](MULTIVARIATE_EFFECTIVE_PRIVATE_THERMAL_ENDPOINT.md#effective-combination-near-the-physical-endpoint)
use the single-order angle confinement and equality of imaginary
physical fields. Thus they apply to these angular exclusions.
For the bounded cubic remainder use the
[effective bounded cubic angular theorem](MULTIVARIATE_HALF_LINE_CUBIC_ANGULAR_SPEED.md#effective-angular-exclusion-in-bounded-ordinary-private-windows).
For the bounded quadratic remainder use the bounded chart above.
This gives positive rational \(\delta_*,\xi_*\) and a computable
cutoff for angular exclusion when \(0\le\sigma\le1+\delta_*\)
and both \(|\xi_j|\le\xi_*\).

The [positive entry theorem](MULTIVARIATE_EFFECTIVE_ONE_PRIVATE_POSITIVE_ENTRY.md)
is a single-root localization. It places both equal-angle
fields in this neighborhood near the positive real ray for
\(0\le\sigma\le1+\delta\), with a computable \(\delta>0\)
and a computable positive strip width. The negative boundary
and interior comparisons above cover its complementary angles.
Use the regular-band angular theorem with
\(\kappa=\delta/2\), \(K=64\) for the remaining profiles.
All cutoffs in the combination are computable from strict
effective witnesses, and taking their maximum proves the theorem.

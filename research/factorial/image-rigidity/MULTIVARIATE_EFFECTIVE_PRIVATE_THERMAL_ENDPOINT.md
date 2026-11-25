# Effective thermal fields at a one-private linear endpoint

Use the ordinary one-private-linear profiles and fields \(\xi_j\) of the
[moving cubic model](MULTIVARIATE_MOVING_PRIVATE_CUBIC_MODEL.md).
Write \(\mu=1-\sigma\), and let \(c_2\) be its positive cubic
coefficient when \(3/4\le\sigma\le5/4\).
For \(\mu\ne0\), define the thermal scales
\[
 a=\operatorname{sgn}\mu,\qquad r=|\mu|/c_2,\qquad c=c_2r^3,
             \qquad M=nc,\qquad b_j=\frac{r\xi_j}{c}.
\]

**Thermal theorem.** For each rational \(B\ge1\) there are computable
positive rational \(r_T,M_T\) such that consecutive common roots
are impossible whenever \(3/4\le\sigma\le5/4\), \(\mu\ne0\),
\(r\le r_T,\ M\ge M_T\), and \(|b_j|\le B\) for both \(j=n,N=n+1\).
Both signs of \(\mu\) and every integer degree difference are allowed.

**Endpoint consequence.** There are computable positive rational
\(\delta_*,\xi_*\) and an integer \(N_*\) such that no consecutive
common root occurs for \(N\ge N_*\) when
\[
                    0\le\sigma\le1+\delta_*,
                    \qquad |\xi_n|,|\xi_N|\le\xi_* .
                                                               \tag{1}
\]
This supplies an effective exclusion once the physical fields
have entered the stated endpoint neighborhood.
It does not localize arbitrary common roots into that neighborhood.

## The pure thermal equality set

Use the fixed \(L=1/1024\), the explicit positive lower and upper
bounds \(c_*=7/64\) and \(c^*=41/32\) for \(c_2\), and the phase
remainder \(E\) from the
[effective bounded cubic expansion](MULTIVARIATE_EFFECTIVE_BOUNDED_PRIVATE_CUBIC.md#uniform-real-lattice-expansion).
At the thermal scale the phase is
\[
  \Psi_b(z)=P_b(z)+(r/c_2)z^4E(rz),\qquad
                 P_b(z)=bz-a z^2/2-z^3/3 .
\]
By conjugation take \(\Im b_j>0\) at a proposed common root.
The physical logarithm is near a positive real scale, so
a real \(\xi_j\) would mean positive real \(t\) and no root.

For \(b\) in the upper half-plane, write the pure right saddle as
\(z_0=x+iy=(-a+\sqrt{a^2+4b})/2\), with the square root of
positive real part. It has \(y>0,\ 2x+a>0\), and
\[
 \Re b=x^2+a x-y^2,\quad \Im b=(2x+a)y,\quad
                D_0(b)=P_b(z_0)=a z_0^2/2+2z_0^3/3.       \tag{2}
\]
For \(a=1\) temporarily restrict to \(|b|\ge b_0>0\), for rational
\(b_0\), in addition to \(|b|\le B\).
Define \(E_a\) in the closed upper half-plane by the polynomial
relations for \(b,x,y\) in (2), \(\Re D_0=0\), and
\(y\ge0,\ x\ge0,\ 2x+a\ge0\), with the corresponding field caps.

The set is compact, since \(z_0^2+a z_0=b\) bounds \(|z_0|\).
On it, \(x,2x+a,|b|\) have strictly positive computable lower
bounds. Indeed for \(a=1,x=0\) the equality forces \(y=0,b=0\),
which was excluded. For \(a=-1,x=1/2\), one has
\(\Re D_0=-1/24-y^2/2<0\).
For \(a=-1,b=0\) the only candidate is \(z_0=1\), with
\(D_0=1/6\). In particular there are no ramified saddles in
\(E_a\). Its only point on the real line is \(b=-3/16\) for
\(a=-1\), where \(z_0=3/4\). The lower bounds are effective
by real polynomial inequalities on this compact set, for example
by rational bisection with a decision procedure for real closed
fields. If \(E_a\) is empty the subsequent covering argument
uses only the strict regions outside it.

## Effective contour neighborhoods

For the exact integer-order lattice sum use the
[finite cotangent representation and its derivative bounds](MULTIVARIATE_EFFECTIVE_PRIVATE_CUBIC_FIELD.md#effective-lattice-representation-on-a-fixed-small-rectangle).
They give, uniformly with computable constants,
\[
 I_j=\int_0^L A_j(q)e^{j\Phi(q,\xi_j)}\,dq+B_j(\xi_j)
                  +O(e^{-\gamma j}),\quad
 |B_j|+|\partial_\xi B_j|\le C/j,\quad |\partial_j B_j|\le C/j^2 .
                                                               \tag{3}
\]
Their hypotheses are ensured by taking \(r_T\) small:
\(|\mu|=c_2r\), \(|\xi_j|\le Bc_2r^2\), and the uniform
physical real gap beyond \(L/2\) is the one from the bounded
cubic expansion. The same bounds hold on the compact \(q\)
contour neighborhoods used in (3), and the amplitude ratio
at consecutive orders there is \(1+O(n^{-2})\).

For each fixed model field \(b\), the following explicit contours
give strict one-contribution estimates outside \(E_a\), and an
endpoint plus a simple saddle in a neighborhood of \(E_a\).
All long horizontal lines below have bounded imaginary parts,
by (2) and the field cap. For large positive real coordinate \(s\)
and \(rs\le L\), one has \(\Re\Psi_b(s+iv)\le-s^3/8\)
on these lines. Indeed take \(s\ge U=8(B+Y+3)\), where
\(Y\ge2\) bounds the chosen heights. Then
\(|s+iv|\le9s/8,\ \Re(s+iv)^3\ge15s^3/16,\ r|s+iv|\le9L/8\).
The quartic part is bounded by
\(8(9/8)^4(L/c_*)s^3\), and the linear and quadratic terms
contribute at most \(9Bs/8+(9s/8)^2/2\).
These bounds give the asserted strict loss.
The closing segment to the real \(q=L\) has exponential
physical loss by taking \(r_TY\) small in the fixed
rectangle of (3).

For \(\Im b>0\) and \(x>0\), use the vertical leg from zero to
\(iy\), then the horizontal ray through \(x+iy\).
On the vertical leg,
\[
 \Re P_b(iv)=-v\Im b+a v^2/2 .
\]
This is at most a negative constant times \(v\) for \(0<v\le y\).
When \(a=1\), use \(\Im b=(2x+1)y>y\); when \(a=-1\)
both terms are nonpositive. On the horizontal ray the real
derivative is \((x-s)(x+s+a)\) at \(s+iy\).
For \(a=1\) the only maximum is \(s=x\).
For \(a=-1\) the maximum is at \(s=x\) or \(s=0\), and the
latter has strictly negative real height by the vertical loss.
At \(s=x\) the real second derivative is \(-(2x+a)<0\).
The following alternatives therefore have strict estimates.
If \(\Re D_0<0\) the horizontal piece is exponentially small.
If \(\Re D_0>0\) the saddle dominates the rest.
At \(\Re D_0=0\) both contributions are active and all parts
away from their two neighborhoods have fixed loss from both.
For \(a=1,\Im b>0,x\le0\), one has
\(\Re b=x^2+x-y^2<0\); the positive real line then gives a
strictly descending endpoint contribution with no saddle needed.

For real \(b\), when \(a=1,b<0\) the positive real line has
strict linear loss. When \(a=1,b>0\) its positive saddle
is a strict positive-height maximum.
When \(a=-1,b<-3/16\) the real line has strict linear loss
from the identity \(P_b(s)=s(b+s/2-s^2/3)\) and the real
maximum \(3/16\) of \(s/2-s^2/3\).
For \(a=-1,b>-3/16\) its positive saddle is a strict
positive-height maximum. These real-phase single-contribution
bounds persist in field neighborhoods. For a positive-height
saddle, shift to the horizontal line through the continued
saddle. The short vertical leg remains below its positive
height, and the strict loss on the real complement persists
under this small shift.

At \(a=-1,b=-3/16\), keep the real endpoint leg \([0,1/8]\).
On it the real phase derivative is at most
\(-3/16+1/8=-1/16\).
At its end attach the vertical connector to height \(y\),
then the right horizontal line through the continued saddle.
At this field the connector is degenerate and its base real
height is strictly below both maxima. On the rest of the
horizontal segment the saddle is a strict maximum, since
the only other possible stationary point is the local minimum
at \(s=1/4\). These strict inequalities persist near the real
coexistence field. Keep this finite connector inside the
saddle contribution when comparing the two orders. Although
it has fixed loss in the scaled phase, its inclusion in the
linked ratio is needed when \(M\) is large and \(r\) is small
independently.

Here is an explicit way to extract uniform constants from the
contour constructions. Choose a smaller interval about the
model saddle where \(\Re(-P_b'')>0\). For fixed model fields
the indicated real losses on the rest of each finite contour
are strict, allowing rational lower margins and rational
parameter neighborhoods. They can be tested by polynomial
inequalities using (2) and the finite range \(0\le s\le U\).
The endpoint linear estimates include zero by bounding the
phase divided by its endpoint coordinate. For the Gaussian
interval use its nonzero Hessian and Taylor expansion.
Rouché continues every required saddle under the actual
\(O(r)\) perturbation, with an explicit error from the
bounds for \(E,E'\) and the pure derivative's nonzero
factor on a small saddle circle.
The actual contours are the same formulas with the
continued saddle in place of \(z_0\).
All strict margins persist under a computably small
perturbation and the tail beyond \(U\) is controlled above.
Thus rational boxes with certified contour bounds can be
enumerated and their coverage of a closed field set tested
by real polynomial inequalities. One may use larger boxes
for the strict bounds and smaller boxes for the cover.
First cover the compact \(E_a\) by two-contribution boxes,
and then cover the closed remainder of the field cap by the
strict single-contribution boxes. The contour constructions
prove termination of these finite-cover searches. They need
no test of an unevaluated limiting contour or path.

The leading endpoint coefficient on either an initial real leg
or an initial imaginary leg is \(e(b)=-1/b\). The leading
saddle coefficient on the horizontal segment is
\[
              s(b)=A_0(r z_s)\sqrt{2\pi/(-\Psi_b''(z_s))},
\]
using the root on the right half-plane for the negative Hessian.
Both coefficients and their inverses are bounded effectively
on the two-contribution boxes.
The absolute linear and Gaussian integral bounds and their
relative asymptotic errors are those of the
[cubic-field calculation](MULTIVARIATE_EFFECTIVE_PRIVATE_CUBIC_FIELD.md#angular-confinement-of-each-root).
For the finite scaled connector at coexistence use its strict loss
to include it in the relative error for the saddle.
The term \(B_j\) has relative size \(O(c/r)\) at the endpoint.
It follows that for sufficiently small \(r\) and large \(M\)
a root must be in a two-contribution box.
It then satisfies, with \(D(b)=\Psi_b(z_s(b))\),
\[
 jc\Re D(b_j)=-\tfrac12\log(jc)
                +\log|e(b_j)/s(b_j)|+O((jc)^{-1/2}+c/r).
                                                               \tag{4}
\]
These estimates also show that if \(a=1\) and \(|b|\ge b_0\),
the imaginary part of any candidate root is bounded below
by a computable positive constant. Choose the covering boxes
near \(E_1\) to preserve half its positive minimum height.

## Consecutive thermal comparison

At a common root the two scaled fields \(b,b'\) have the same
imaginary part, since \(r,c\) are fixed by the profile.
At both fields the pure saddle has \(\Re z_0\) and
\(\Re(2z_0+a)\) bounded below by the two-contribution bounds.
These bounds persist along the horizontal segment between the
fields. At fixed height \(Y=\Im b>0\), one has
\(y=Y/(2x+a)\) and the real target coordinate is
\[
             x^2+a x-Y^2/(2x+a)^2.
\]
Its derivative in \(x\) is strictly positive when \(2x+a>0\).
Thus \(x\) increases along the horizontal direction of increasing
real field, so its lower bounds at the endpoints suffice.
The actual saddle continues on that segment with
\(\Re D'(b)=\Re z_s(b)\ge\kappa>0\), using the uniform
Rouché estimate and small \(r_T\).
At either root (4) bounds \(\Re D\) in modulus by
\(C(1+\log M)/M\). Hence \(|b'-b|\le C(1+\log M)/M\).
By taking \(M_T\) large, both fields lie in the same larger
two-contribution box of the finite cover, and the analytic
coefficient ratio is Lipschitz along the segment.

Write \(M'=Nc,\ d=M'(b'-b)=rN\Delta\log(N/n)\).
Subtracting (4) and using \(\Re D'\ge\kappa\) gives the
effective preliminary bound \(|d|\le1\) for sufficiently
large \(M\) and small \(r\), as in the
[cubic linked ratio proof](MULTIVARIATE_EFFECTIVE_PRIVATE_CUBIC_FIELD.md#linked-endpoint-and-saddle-ratios).
Explicitly the error terms before absorbing the Lipschitz
contribution \(O(|d|/M)\) are
\(O((1+\log M)/n+M^{-1/2}+c/r)\).

Now use the exact lower-order contours at both orders. The exact
increment in exponent is \(c\Psi_b(z)+dz\), and the amplitude
ratio is \(1+O(n^{-2})\). Thus the endpoint and saddle
comparison equations (7), (8) of the cubic-field proof apply.
The special coexistence connector has fixed negative loss
relative to the saddle. After subtracting the central
multiplier its relative error is at most
\(C(c+|d|)\sqrt M e^{-\gamma_2 M}\).
This is bounded by \(C'(c/M+|d|/\sqrt M)\), retaining the
necessary factor of the increment. The same argument applies
to the bounded horizontal complement. On the long horizontal
tail use the cubic loss and \(|\Psi_b|=O(s^3)\).
The physical remainder has size \(O((M/r)e^{-\gamma n})\)
after division by the endpoint scale.

If \(\Delta\ne0\), then \(|d|\ge r\).
Use \(c/r=c_2r^2\), \(c/r^2=c_2r\), and \(M=nc\);
all ratio errors divided by \(|d|\) are uniformly small
under effective thresholds, as is \(|c\Re D|/|d|\).
The saddle logarithmic modulus difference therefore cannot
equal the endpoint difference, because the term
\(d\Re z_s\) has modulus at least \(\kappa|d|\).
If \(\Delta=0\), (4) gives
\(c\Re D=-(\log M)/(2n)+O(1/n)\), whereas the
remaining ratio errors are \(O(1/n)\). A computable lower
bound on \(\log M\) excludes this case too.
This proves the thermal comparison with the lower field cap
for \(a=1\), and without a lower cap for \(a=-1\).

To remove the lower cap for \(a=1\), take a small fixed
\(b_*>0\) such that \(|b_j|\le b_*\) at both orders implies
the [quadratic corridor](MULTIVARIATE_EFFECTIVE_PRIVATE_QUADRATIC_CORRIDOR.md):
use \(|\xi_j|/\mu^2=|b_j|/c_2\le b_*/c_*\).
Enforce its condition on \(n\mu^3=M c_2^2\) by enlarging \(M_T\).
If either order has \(|b_j|\ge b_*\), its imaginary part is
bounded below by the single-order confinement just proved
with cap \(b_*\). The other field has the same imaginary
part, so both moduli are bounded below by a computable
positive cap. Apply the two-order comparison with this
possibly smaller cap. This proves the thermal theorem.

## Effective combination near the physical endpoint

Let \(r_C,a_C,M_C\) be the constants of the
[escaping cubic corridor](MULTIVARIATE_EFFECTIVE_PRIVATE_CUBIC_FIELD.md).
Choose a rational \(B\ge\max(2,2a_C^{-2})\), and obtain
\(r_T,M_T\) from the thermal theorem.
Choose \(K_1\) so that \(K_1/(c^*)^2\ge\max(M_T,M_C)\).
Take \(\mu_0>0\) small enough that
\(\mu_0\le1/4\) and \(\mu_0/c_*\le r_T\), and impose
\(\xi_*\le c_*r_C^2\).

Suppose first \(|\mu|\le\mu_0\) and \(n|\mu|^3\ge K_1\),
so the thermal scale exists and \(M\ge M_T\).
If both \(|b_j|\le B\) the thermal theorem applies.
Otherwise at the larger field, the cubic scale has
\[
 r_j^C=r\sqrt{|b_j|},\qquad
 \frac{|\mu|}{c_2 r_j^C}=|b_j|^{-1/2},\qquad
 M_j^C=j c |b_j|^{3/2}.
\]
Dropping from that field to its imaginary part loses at
most a factor two by the single-order angle confinement in
the cubic corridor. Because the physical imaginary parts
agree, the other field has thermal modulus at least \(B/2\).
Thus both satisfy the cubic corridor hypotheses, including
the upper bound on \(r_j^C\) by \(\xi_*\), a contradiction.

If \(n|\mu|^3<K_1\), including \(\mu=0\), choose a fixed
rational \(X\) so large that a field
\(|\xi_j|\ge X n^{-2/3}/2\) would satisfy
\[
 \frac{|\mu|}{\sqrt{c_2|\xi_j|}}\le a_C,\qquad
           \frac{j|\xi_j|^{3/2}}{\sqrt{c_2}}\ge M_C .
\]
This is possible using \(|\mu|\le K_1^{1/3}n^{-1/3}\) and
the bounds on \(c_2\). If one field has modulus at least
\(X n^{-2/3}\), its cubic angular confinement forces the other
above \(Xn^{-2/3}/2\), so again both satisfy the cubic corridor
hypotheses. In the remaining case both physical fields are
at most \(X n^{-2/3}\). The cubic fields and temperatures
are then bounded by rational caps depending just on
\(X,K_1\) and the bounds for \(\ell\). The
[effective bounded cubic theorem](MULTIVARIATE_EFFECTIVE_BOUNDED_PRIVATE_CUBIC.md)
applies above its computable cutoff.

Finally if \(0\le\sigma\le1-\mu_0\), take
\(\xi_*\le p_*\mu_0^2\), where \(p_*\) is the quadratic corridor
constant. Then \(|\xi_j|\le p_*\mu^2\), and
\(n\mu^3\ge n\mu_0^3\) eventually meets its computable
threshold. Taking \(\delta_*=\mu_0\) and one cutoff above
all the stated thresholds proves (1).

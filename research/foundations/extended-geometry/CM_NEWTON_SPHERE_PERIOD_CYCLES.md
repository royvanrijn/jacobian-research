# Sphere monodromy cycles at the cyclic limit

The [full-width strip](CM_NEWTON_FULL_WIDTH_STRIP_MOMENTS.md)
relates moment residues to coherent profiles on a cyclic cover.
A nonzero profile can have a nonzero fibre de Rham class even if
the sphere resolvent misses it. The distinction is which fibre
cycles occur in the variation of the sphere.

## The core limit and a relative path for the sphere

**Proposition (compact cycle periods at infinity).**
For the full-width phase
\(f=t^nQ+t^{n-1}B+t^{n-2}D\), set
\(H=t^{s-1}A+t^{s-2}C\) and use the fibre \(f=\varepsilon^{-n}\).
Let \(K\) be a compact subsurface in the \(z\)-sphere avoiding
\(0,\infty\) and the roots of \(P=zQ\). The cover over \(K\) at
\(\varepsilon=0\) is identified with the nearby covers by continuing
the simple roots in \(T=\varepsilon t\). For any closed cycle
transported by this identification, its normalized fibre period
\(\varepsilon^{-(n-s)}\int H\omega/df\) tends, up to sign, to
\[
             \frac{1}{n}\int Y^{-s}A(z)\,\frac{dz}{z},
             \qquad Y^n=Q.
\]

Indeed the equation on the scaled fibre is
\[
         T^nQ+\varepsilon T^{n-1}B
                         +\varepsilon^2T^{n-2}D=1 .
\]
The roots at zero are locally distinct on \(K\), and compactness
makes their continuation uniform. With a fixed sign convention
the normalized differential on a sheet is
\[
 -\frac{T^{s-1}A+\varepsilon T^{s-2}C}
 {nT^{n-1}Q+(n-1)\varepsilon T^{n-2}B+(n-2)\varepsilon^2 T^{n-3}D}
 \,\frac{dz}{z}.
\]
At zero the denominator is \(n/T\) and \(T=1/Y\).
This gives the stated limit with uniform error \(O(\varepsilon)\)
on compact subcharts. The leading correction term vanishes there.
A large counterclockwise loop in the value of \(f\) rotates the
limiting \(T\) sheets by \(e^{2\pi i/n}\) under these identifications.
The original sphere is fixed on that loop since \(f\) is bounded
on it. The fibre cycles for the criterion must come from other
monodromy loops. A limiting projection of a cycle that also has
parts outside the compact continuation requires control of the
periods on those parts.

There is a direct description of the sphere as a relative fibre
integral. In the literal \(x,t\) chart write
\[
 g_\lambda(x,t)=x(f(x,t)-\lambda)
       =P(0)t^n(t^2-1)+x(f_+(x,t)-\lambda),
\]
where \(f_+\) is polynomial. When \(|\lambda|\) is larger than
the maximum of \(|f|\) on the sphere, Rouché's theorem on each
latitude circle \(|x|=\sqrt{(1-t^2)/2}\) shows that \(g_\lambda\)
has exactly one simple interior zero \(\xi_\lambda(t)\), for
\(-1<t<1\). The comparison term is \(-\lambda x\).
The markers are polynomial in \(x,t\), so taking the angular
residue gives
\[
 \mathcal A\!\left(\frac{H}{f-\lambda}\right)=
 \frac12\int_{-1}^{1}
 \frac{H(\xi_\lambda(t),t)}
      {\partial_xg_\lambda(\xi_\lambda(t),t)}\,dt .
\]
For \(t\ne 0,\pm1\) the denominator is \(xf_x\) on the lifted
path. At \(t=0\) the branch \(\xi_\lambda(0)=0\) is smooth in
the plane fibre equation with derivative \(-\lambda\). It
corresponds to an end of the affine fibre after normalization.
The displayed differential extends on this branch. At \(t=\pm1\)
the branch tends to \(x=0\) as the latitude shrinks and is smooth
for large \(|\lambda|\). Thus the sphere integral follows one
specified relative path on the normalized fibre with its two
endpoints over \(t=\pm1\). Tracking this path under monodromy
gives a concrete formulation of the remaining cycle problem.

## A separable conductor-18 phase

**Proposition (a separable conductor-18 phase).**
On \(2xy+t^2=1\), put \(z=x/t^2\) and take
\[
 f=S(t)Q(z),\qquad S(t)=t^{16}(t^2-1),\qquad
 Q(z)=\frac{(z^2+1)^4}{z}.
\]
For \(n=18,s=13\), the two-dimensional space of protected
coherent profiles is
\[
 A(z)=(z^2+1)^2(a_0+a_1z).
\]
For this phase, two sphere monodromy variations, based at a large
negative fibre value, have limiting projections on the normalized
cover \(Y^{18}=Q(z)\) which detect every nonzero profile \(A\).
More precisely, the leading fibre periods for
\(H=t^{12}A(z)+t^{11}C(z)\), with \(\deg C\le5\), are nonzero
constant multiples of
\[
 a_0 I_0+a_1 I_1,\qquad a_0 I_0-a_1 I_1,\qquad
 I_j=\int_0^\infty z^j(z^2+1)^2Q(z)^{-13/18}\,\frac{dz}{z}.
\]
Both \(I_j\) are positive and finite. With \(x=z>0\) the
integrand is \(x^{j-5/18}(1+x^2)^{-8/9}\,dx\).
The substitution \(v=x^2\) gives
\[
 I_0=\tfrac12 B(13/36,19/36),\qquad
 I_1=\tfrac12 B(31/36,1/36).
\]
After dividing out the two nonzero sheet factors the period
matrix has rows \((I_0,I_1)\) and \((I_0,-I_1)\). Its
determinant is \(-2I_0I_1\). In particular, simultaneous
vanishing of all the corresponding weighted sphere moments forces
\(A=0\).

The phase is regular. Its negative-angle part is
\(x^{-1}t^{18}(t^2-1)=-2yt^{18}\), and its other monomials have
nonnegative heights. It is a full-width phase with \(B=0,D=-Q\).
The protected profiles have degree at most five and vanish to
order at least two at both roots of \(z^2+1\).

## Two variations of the sphere

Let \(z_c=1/\sqrt7\), the unique positive critical point of \(Q\),
and let \(t_c=\sqrt{8/9}\). On \(0<z<\infty\), \(Q\) has a strict
minimum at \(z_c\) and tends to infinity at both ends.
Put \(S_c=S(t_c)<0\), \(Q_c=Q(z_c)>0\) and
\(\lambda_-=S_cQ_c\), \(\lambda_+=-\lambda_-\).

Represent the real sphere by a cylinder over \(-1\le t\le1\)
with \(x\) on a circle and \(y=(1-t^2)/(2x)\). Deform the radius
to \(z_ct^2\) for \(|t|\ge\delta\) and \(z_c\delta^2\) for
\(|t|\le\delta\), for small \(\delta>0\). At \(t=\pm1\), cap the
circles by disks in \(y=0\). Interpolating the radii from those
of the real sphere and interpolating the caps gives a homology
between the cycles in the fibre complement at sufficiently large
\(|\lambda|\). The interpolation has compact image.

On the part with \(|t|\ge\delta\), one has \(|z|=z_c\). The
positive coefficients of \(Q\) give \(|Q(z)|\le Q_c\), with
equality only at \(z=\pm z_c\). The maximum of \(|S(t)|\) on
the real interval occurs only at \(t=\pm t_c\). The new caps
have \(f=0\). On the small central cylinder, the inequality
\(|f|\le\delta^{16} Q_c\) follows term by term since
\(16-2j\ge0\) for every angular exponent \(j\) of \(Q\).
Thus the deformed sphere misses the fibre for \(|\lambda|>
|\lambda_-|\). Its first contact along the negative real ray
occurs at \((z_c,\pm t_c)\), and along the positive real ray
at \((-z_c,\pm t_c)\).

At \((z_c,t_c)\) use analytic coordinates
\(Q=Q_c(1+\xi^2)\) and \(S=S_c(1-\tau^2)\), chosen real on the
respective real lines through \(z_c,t_c\). Then
\[
 f-\lambda_-=u^2+v^2,\qquad
 u=\sqrt{-\lambda_-}\,\tau,\qquad
 v=i\sqrt{-\lambda_-}\sqrt{1-\tau^2}\,\xi .
\]
On the cylinder \(\tau\) is real and \(\xi\) is purely
imaginary to first order in the angular coordinate. Its
real tangent plane is therefore \((u,v)\in\mathbb R^2\).
At a nearby negative-ray value \(f-\lambda_-=-\eta\) with
\(\eta>0\), the dual thimble is the disk
\(\{(u,v)\in i\mathbb R^2,\ |u|^2+|v|^2\le\eta\}\).
It meets the sphere transversely at the origin of this
coordinate chart, with local intersection number \(\pm1\).
A small loop about \(\lambda_-\) consequently adds a Leray
tube on its boundary vanishing circle. This is the one-divisor,
ordinary quadratic case of [Pham, §§4.1 and 5, equations
(5.2)–(5.3)](https://www.numdam.org/item/10.24033/bsmf.1628.pdf).
The complement of small neighbourhoods of the two contact
points on the sphere misses every fibre over a sufficiently
small disk about \(\lambda_-\), so the variation is local at
those two points.

Let \(\delta_+\) be the circle from \((z_c,t_c)\). The reflection
\(\rho(z,t)=(z,-t)\) preserves the phase and reverses the
sphere orientation. The two local contributions are therefore
the tube of \(\delta_+-\rho_*\delta_+\), up to a common sign.
The rotation \(\sigma(z,t)=(-z,t)\) preserves the sphere
orientation and sends \(f\) to \(-f\). It gives the corresponding
second variation about \(\lambda_+\). Take the paths from the
two exterior real rays and return to the same large negative
base point along an exterior upper semicircle for the second
variation. The sphere itself stays fixed on that semicircle.
Loops entirely in the exterior, including local monodromy at
infinity, give zero variation of the sphere.

## Transport and limiting periods

For \(\lambda=-R<\lambda_-\), the equation
\(Q(z)=R/|S_c|\) has exactly two positive solutions. Between them
\(S(t)=-R/Q(z)<S_c\). The two roots issuing from \(t_c\) are
complex conjugates. They stay in \(0<\pm\arg t<\pi/18\):
the real positive boundary has \(S(t)\ge S_c\), and
\(S(t)\) has nonzero imaginary part on \(\arg t=\pm\pi/18\).
Traverse the interval on one root and return on the other.
The endpoints are joined where both roots equal \(t_c\).
As \(\lambda\to\lambda_-\) this cycle is the vanishing
circle just described, since the positive \(z\) interval
and its two \(t\) branches give the circle in the imaginary
\((u,v)\)-plane. It thus represents the transported \(\delta_+\).

Put \(\varepsilon=R^{-1/18}\), \(T=\varepsilon t\). The scaled
fibre satisfies \(T^{16}(T^2-\varepsilon^2)Q=-1\).
On each compact subinterval of \(0<z<\infty\), the two paths
limit to \(T=e^{\pm i\pi/18}Q(z)^{-1/18}\). The endpoints
of the \(z\)-interval tend to \(0\) and \(\infty\), which are
both fully ramified in \(Y^{18}=Q\). In this cover use the
identification \(Y=e^{i\pi/18}/T\). The limiting cycle is the
difference between two adjacent lifts of the positive real
ray from \(0\) to \(\infty\).

For the second variation, the rotation \(\sigma\) first gives
the corresponding paths on the negative real axis in the
fibre at \(+R\). For the exterior transport
\(\lambda=R e^{i\theta}\), \(0\le\theta\le\pi\), retain the
central part of that negative ray where \(|Q|\le M/2\), with
\(M=R/|S_c|\). Its two endpoints are continued to the moving
branch points by local inverses of \(q=Q(z)\) near zero and
infinity. In each inverse chart, join \(q=-M/2\) to
\(q=-(M/2)e^{i\theta}\) by the circle of radius \(M/2\), then
join it radially to \(-M e^{i\theta}\).
On the angular and central parts \(|\lambda/q|\ge2|S_c|\),
so both roots of \(S(t)=\lambda/q\) continue without
collision. On the radial tip \(\lambda/q\) is again negative
real, approaching \(S_c\) at the end and closing the lift.
This is a continuous family of cycles along the exterior
transport. In the central part \(\arg\lambda\) increases by
\(\pi\), so the scaled branches acquire the common factor
\(e^{i\pi/18}\) in the limit. The resulting cycle at the
same scaling limit is, up to this common sheet choice, the
difference of adjacent lifts over the negative ray. The
explicit paths transport compact cycles. If a value where
the full nonproper family fails to be locally trivial lies
on a radial base path, one can take a small detour there:
the sphere and the matched paths already continue in a
compact regular neighbourhood along the ray.

For the period calculation use the sign of the fibre differential
for which \(H\omega/df=H\,dz/(zf_t)\) in this chart, with
\(\omega=dz/z\wedge dt\). Multiplying by \(R^{(18-s)/18}\)
and passing to the scaled central fibre gives
\[
 -\frac1{18}T^s A(z)\,\frac{dz}{z}.
\]
Indeed \(f_t=Q S'(t)\),
\(S'(t)=t^{15}(18t^2-16)\), and \(Q T^{18}=-1\).
The \(C\) term has one extra factor of \(\varepsilon\) on
compact subintervals.

The ends of the paths give no additional leading periods.
On the negative radial base ray, put \(w=R/Q(z)\) on the
positive \(z\) matching interval and \(w_0=|S_c|\). On either
\(t\)-branch there is a uniform bound \(|t|\ge c>0\) and
\[
 w^{(18-s)/18}
 \left|\frac{t^{s-1}}{S'(t)}\right|
 \ \le\ K(1-w_0/w)^{-1/2}.
\]
This follows from the simple square-root branch at \(S_c\)
and \(S(t)\sim t^{18}\) as \(|t|\to\infty\) along the
tracked branch. It bounds the absolute scaled integrand by
a constant times
\[
 |A(z)| Q(z)^{-13/18}(1-Q(z)/M)^{-1/2}\,\frac{dz}{z}.
\]
At the two ends use \(q=Q(z)\), where
\(|dz|/|z|\le K|dq|/|q|\). The powers of \(q\) in
\(|A|q^{-13/18}\) are bounded above by \(-13/18\) at zero
and \(5/7-13/18=-1/126\) at infinity. For each \(\alpha>0\)
the integral of \(q^{-1-\alpha}(1-q/M)^{-1/2}\) from
\(q_0\) to \(M\) is \(O(q_0^{-\alpha})\), by splitting
at \(M/2\). On the exterior semicircle the central portion
has \(|w|\ge2w_0\), where the same estimate holds with a
constant independent of the argument, by the large-root
asymptotic and compactness away from the critical values.
In the added tails \(|q|\ge M/2\), the bounded angular legs
and the same square-root estimate on the radial tips give
\(O(M^{-13/18})\) at zero and \(O(M^{-1/126})\) at infinity.
The bounds apply as well to \(C\), whose degree is at most
five. Its factor \(1/t\) is uniformly bounded on all
branches and tends to zero on compact subintervals. Hence
the nuisance contribution has zero limit.

The two sheets of the positive matching path differ by
\(e^{2\pi i/18}\), whose \(s\)-th power is not one.
Reflection in \(t\) multiplies the limiting coherent form
by \((-1)^s=-1\), so the two terms of
\(\delta_+-\rho_*\delta_+\) add. The positive-ray period is
therefore a nonzero constant times \(a_0I_0+a_1I_1\).
For the other variation write \(z=-r\), \(r>0\), on the
negative ray. Then \(dz/z=dr/r\) and
\(A(-r)=(r^2+1)^2(a_0-a_1r)\). The same adjacent sheets and
odd reflection produce a nonzero constant times
\(a_0 I_0-a_1 I_1\). The two leading periods separate the
profiles since \(I_0,I_1>0\).

Finally, the expansion of the sphere resolvent at infinity
has coefficients \(\mathcal A(Hf^k)\). If they all vanish,
its analytic continuations have zero variation. The Leray
residue formula equates each tube integral with the
corresponding fibre period up to a nonzero universal factor.
The two limits then force \(a_0=a_1=0\).

## Odd radial perturbation

Write \(g_b=S(t)Q(z)+bt^{17}\), with \(S,Q\) as above. The
[finite moment calculation](CM_NEWTON_ODD_RADIAL_PERTURBATION.md)
separates the coherent component for every \(b\). The result here
identifies sphere cycles and estimates their periods.

**Theorem (persistence of two detecting sphere variations).**
The two coherent coefficients are detected by sphere-generated
variations for every real \(b\) with \(0<|b|<2Q_c/17\), where
\(Q_c=Q(1/\sqrt7)\). They are also detected for all sufficiently
small complex \(b\). In each case use a variation reached from
each exterior real ray. Their normalized leading periods are
nonzero multiples of \(a_0I_0+a_1I_1\) and \(a_0I_0-a_1I_1\).
Every allowed correction marker has zero normalized limit on
these variations.

### Real parameters and single contacts

It suffices to take \(b>0\), since reflection in \(t\) changes
\(b\) to \(-b\). Use the cylinder and caps of the separable
phase. On the part with \(|z|=z_c\) put \(u=|t|\). Then
\[
 |g_b|\le Q_c u^{16}(1-u^2)+bu^{17}.
\]
The right side has its unique maximum at
\(u=c(b/Q_c)<1\), where \(c(\beta)\) is the positive solution of
\(18c^2-17\beta c-16=0\). This follows by differentiating, and
\(c<1\) is equivalent to \(17b<2Q_c\).
The maximum is strictly greater than the cap value \(b\).
The small central cylinder has arbitrarily small image as its
height interval shrinks. Equality at the maximum requires
\(z=\pm z_c\) and alignment of the two terms in \(g_b\).
Hence there is one contact on the negative ray at
\((z_c,-c(b/Q_c))\), and one on the positive ray at
\((-z_c,c(b/Q_c))\).

At the negative contact the Hessian is diagonal in \(t,z\)
because \(Q'=0\). It is positive on the real \(t\) direction
and negative on the real \(z\) direction. The angular tangent
\(iz_c\) makes the quadratic form positive on the real
tangent plane of the cylinder. In local holomorphic Morse
coordinates the imaginary thimble therefore meets that plane
transversely with intersection number \(\pm1\). The same
quadratic-pinch formula as in the separable case gives one
tube from this contact. The map \((z,t)\mapsto(-z,-t)\) sends
\(g_b\) to \(-g_b\) for every \(b\) and gives the other
variation, up to orientation.

To track the negative-ray circle write \(q=Q(z)>0,\ t=-u\) and
\[
 \varphi_\beta(u)=u^{18}-u^{16}-\beta u^{17},\qquad
 \beta=b/q,\qquad w_c(\beta)=-\varphi_\beta(c(\beta)),\qquad
 L(q)=q\,w_c(b/q).
\]
The derivative \(L'(q)=-S(c(b/q))\) is positive, bounded above,
and bounded away from zero on \(q\ge Q_c\). Thus the matching
interval at \(\lambda=-R\) for \(R>L(Q_c)\) is the positive
\(z\)-interval with \(q\le M=L^{-1}(R)\).
The two branches of \(\varphi_\beta(u)=-R/q\) born at
\(u=c(\beta)\) stay in the sectors \(0<\pm\arg u<\pi/18\).
On the positive real boundary the value \(-R/q\) is below
the minimum. On the other boundary the two terms
\(-u^{16}\) and \(-\beta u^{17}\) have imaginary parts of
the same strict sign in their sum, while \(u^{18}\) is real.
The derivative has no zero in either open sector. These lifts
join at the two endpoints of the matching interval, and
shrink to the local vanishing circle as \(R\downarrow L(Q_c)\).
For \(R>L(Q_c)\) the join points are regular in the fibre,
since \(S(c)\ne0\) and \(Q'\ne0\) at both ends.

On each compact subinterval of the positive ray,
\(R^{-1/18}u\) tends on these branches to
\(e^{\pm i\pi/18}q^{-1/18}\). The sign \(t=-u\) gives a
common sheet factor. For the negative \(z\)-ray use the
image under \((z,t)\mapsto(-z,-t)\). Its polynomial factor
in the limiting coherent integral is \(A(-z)\) at positive
\(z\). This gives the two stated leading periods, provided
the tails have zero limit.

Here is a uniform estimate for the tails. Set \(w=R/q\).
The simple minimum of \(\varphi_\beta\) and its leading
term \(u^{18}\) give
\[
 w^{5/18}\left|\frac{u^{12}}{\varphi_\beta'(u)}\right|
   \le K(1-w_c(\beta)/w)^{-1/2},\qquad |u|\ge c_*>0.
\]
Near the minimum this is the simple square-root estimate,
uniform since \(\beta\in[0,b/Q_c]\) is compact and the
second derivative at \(c(\beta)\) is bounded away from zero.
Away from the minimum use compactness for bounded \(w\), and
\(u^{18}\sim -w\) for large \(w\). Since
\(1-L(q)/R\) is comparable to \(1-q/M\), the same integrable
end bounds as in the separable case apply. Explicitly, near
zero and infinity the angular factors in the normalized
integrand are bounded by \(q^{-13/18}\) and \(q^{-1/126}\),
respectively, times \(K|dq|/|q|\) and the square-root factor.
This holds with either \(A\) or \(C\), since both have degree
at most five. The correction term also has a factor \(1/u\),
bounded on the whole interval and tending to zero on every
compact central piece. Splitting off the ends and then
letting \(R\) grow proves that its normalized limit is zero.

### Persistence for small complex parameters

At \(b=0\) take two fixed based loops enclosing the two
contact values, from sufficiently large opposite real base
values. Each transported sphere representative uses a
compact family of cycles disjoint from its moving fibre.
These families persist for small complex \(b\). Indeed
\(|g_0-\lambda|\) has a positive lower bound on their compact
supports, and \(g_b-g_0=bt^{17}\) converges uniformly there.
The homologies equating the variations at the base with the
two signed tube sums also have compact support and persist
by the same argument. Each tube may be continued using
a fixed small normal radius in the regular compact
neighbourhood of its fibre cycle. The swept tubes remain
disjoint from the perturbed base fibre, by taking \(|b|\)
small compared with that radius and the bound of
\(|t|^{17}\) on these neighbourhoods. Thus the actual sphere
variations at small complex \(b\) are the continued signed
sums of the two circles from each contact pair.

For completeness, the continued circles can be tracked
from the negative base to arbitrarily large \(-R\) by
explicit matching paths. For a circle originating at
\(t\) of sign \(\sigma=\pm1\), put \(t=\sigma u\) and use
\(\beta=-\sigma b/q\) in \(\varphi_\beta\).
Let \(w_*=w_c(0)>0\). Retain the central part of the positive
\(z\)-ray where \(q\le q_0=R/(2w_*)\). For \(|b|\) small,
the nonzero critical values of \(\varphi_\beta\) have
modulus less than \(3w_*/2\) on this part, since
\(|\beta|\le |b|/Q_c\). The two roots at \(b=0\) therefore
continue throughout this central piece and remain distinct.
In either end use the inverse chart of \(q=Q(z)\) near
zero or infinity. Put
\[
 v=q\,w_c(-\sigma b/q).
\]
This is holomorphic at large \(q\), with \(v=w_*q+O(b)\),
and has an inverse there in the sector of the positive ray.
Join \(q_0\) to the inverse at \(v=R/2\) along a short
segment of bounded length, then use the inverse for
\(R/2\le v\le R\). The roots along the latter path satisfy
\(\varphi_\beta(u)=-w_c(\beta)R/v\). They continue from
the two roots at \(b=0\), with a uniform square-root
parametrization at \(v=R\), where they join. To justify
uniformity, on the end \(\beta=O(1/R)\), the two model
branches are simple for \(1<R/v\le2\), and only the
endpoint \(R/v=1\) requires the square-root chart.
The short segment has \(R/q\) near \(2w_*\) and no collision.
This constructs closed cycles continuous in \(b,R\) that
agree with the original matching paths at \(b=0\).
The same construction on the other exterior ray follows
from \((z,t)\mapsto(-z,-t)\).

The central parts have the same adjacent-sheet limits.
The ends have \(|q|\) comparable to \(R\), and
\(|dq|/|q|\le K|dv|/|v|\) on the final pieces.
Their normalized integrals are bounded by \(O(R^{-13/18})\)
near zero and \(O(R^{-1/126})\) near infinity, using the
factor \((1-v/R)^{-1/2}\), which is integrable.
On the central pieces the estimates follow as before,
with no square-root singularity. Hence the correction
periods vanish. The two signs of \(t\) in each signed
sphere variation contribute opposite leading coherent
periods because \(s=13\) is odd. The signed sums therefore
retain the two nonzero sheet factors of the separable
phase. Both exterior base rays belong to the same
resolvent germ at infinity. If all marker moments vanish,
both variations have zero period, forcing \(a_0=a_1=0\).

## Correction principal parts on the cyclic limit

The correction sphere resolvent is identically zero for the
separable seed because \(f(z,-t)=f(z,t)\) and
\(t^{11}C(z)\) is odd in \(t\). For the perturbed phase one
has to control the correction periods.

On the central cover \(Y^{18}=Q(z)\), write the coefficient
of the correction differential as
\(\omega_C=Y^{-12}C(z)\,dz/z\); its core contribution carries
one factor of \(\varepsilon\).
At either root \(\alpha=\pm i\), a local parameter with
\(z-\alpha=w^9\) gives \(Y=w^2\) times a unit in \(w^9\),
so the only possible polar powers are \(w^{-16}dw\) and
\(w^{-7}dw\). At infinity use \(z=w^{-18}\). The contribution
from \(z^j\) begins in order \(83-18j\), giving only
\(w^{-7}dw\) for \(j=5\). At zero use \(z=w^{18}\) and
\(Y\) of order \(-1\); all terms are regular. In particular
the polar terms have meromorphic local primitives.

On a path segment, subtracting the differential of such a
primitive replaces its integral by the difference of the
primitive's endpoint values. Those terms must be retained
when comparing lifts and taking limits. On the real rays
the \(j=5\) correction behaves like
\(\varepsilon x^{-2/3}\,dx\); up to cutoff \(X\) its
integral has size \(O(\varepsilon X^{1/3})\).
The ray integrals converge for degrees at most four.
For the conductor-18 degree-\((5,1)\) system, \(\deg C\le2\)
removes this infinity contribution. The two finite roots
can still have polar terms. The matching paths above
avoid those roots and give a direct bound on their own ends.

The [finite-root calculations](CM_NEWTON_CYCLIC_LIMIT_TAIL_PERIODS.md)
exhibit escaping fibre periods beyond this matching-path family.

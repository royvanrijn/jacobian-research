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
Both \(I_j\) are positive and finite. In particular, simultaneous
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

# No circular exterior arc at a bounded three-cell junction

Fix \(0<\rho<1\). Work in the planar active space of an ordinary balanced
three-cell Gaussian noise-stability maximizer. Let \(A_1,A_2\) be bounded
up to null sets and let \(E\) be the unbounded connected component of the
strict winner of \(A_3\). Suppose \(x_*\) is a trivalent genuine triple
junction with one local region of \(E\) and one of each bounded cell.
Each pair has one incident analytic interface with nonvanishing pair
score gradient, as in [the local regularity result](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#regular-interface-arcs-and-isolated-triple-junctions).
The tangent lines may coincide.

**Theorem.** No open part of a \(1/3\) or \(2/3\) interface lies on a
circle. In particular, a compact theta network with one exterior
component and two bounded winner regions has no circular outer arc.

The proof converts the circular score equation into an exterior-potential
identity for the modified Helmholtz operator. The two bounded labels then
give incompatible second derivatives at the junction. A related
corner obstruction for Helmholtz sources with a continuous limiting
source value is proved by Blåsten. The local argument below treats the
jump between the two analytic source pieces at the junction.
[Blåsten, *Nonradiating Sources and Transmission Eigenfunctions Vanish at Corners and Edges*](https://arxiv.org/abs/1803.10917).

## A circular score determines the exterior potential

Write \(U_l=T_\rho1_{A_l}\), where \(T_\rho\) is the Gaussian noise
operator, and let \(\beta_l\) be the constant prices in the strict winner
rule. Interchanging the bounded labels if needed, the circular arc has
\(U_1-\beta_1=U_3-\beta_3\). As \(U_1+U_2+U_3=1\), it follows that

\[
 T_\rho f=c,\qquad f=2\,1_{A_1}+1_{A_2},
 \qquad c=1+\beta_1-\beta_3>0.
\tag{CI1}
\]

The strict positivity follows directly from the positive Gaussian kernel
and the positive mass of either bounded cell. Let \(p\) and \(R>0\) be
the center and radius of the circle. The real-analytic function
\(\theta\mapsto (T_\rho f)(p+R(\cos\theta,\sin\theta))\) equals \(c\)
everywhere because it equals \(c\) on an open interval. Set
\(\sigma^2=1-\rho^2\), \(y_0=\rho p\), and \(\lambda=\rho R/\sigma^2>0\).
Writing the conditional Gaussian density in (CI1) gives

\[
 \int_{\mathbb R^2} e^{\lambda\omega\cdot(y-y_0)}\,d\mu(y)=C
 \quad(\omega\in S^1),\qquad
 d\mu(y)=e^{-|y-y_0|^2/(2\sigma^2)} f(y)\,dy,\qquad
 C=2\pi\sigma^2 c\,e^{\rho^2R^2/(2\sigma^2)}>0.
\tag{CI2}
\]

Let \((r_y,\phi_y)\) be polar coordinates centered at \(y_0\). The
angular Fourier coefficient of order \(n\in\mathbb Z\) of (CI2) gives

\[
 \int I_{|n|}(\lambda r_y)e^{-in\phi_y}\,d\mu(y)
 =C\,1_{\{n=0\}}.
\tag{CI3}
\]

Here \(I_n\) is the modified Bessel function of the first kind, with
integral representation
\(I_{|n|}(t)=(2\pi)^{-1}\int_0^{2\pi}e^{t\cos\theta}e^{-in\theta}\,d\theta\).

The fundamental solution of \(-\Delta+\lambda^2\) in the plane is
\(G_\lambda(z)=(2\pi)^{-1}K_0(\lambda|z|)\). For a point \(y\) with
\(r_y<r\), the order-\(n\) Fourier coefficient of
\(G_\lambda(y_0+r(\cos\theta,\sin\theta)-y)\) as a function of
\(\theta\) is

\[
 \frac{1}{2\pi}K_{|n|}(\lambda r)I_{|n|}(\lambda r_y)e^{-in\phi_y}.
\tag{CI4}
\]

One can obtain (CI4) directly from the radial equation. The coefficient
is a multiple of \(I_{|n|}(\lambda r)\) for \(r<r_y\) by boundedness
at the center and a multiple of \(K_{|n|}(\lambda r)\) for \(r>r_y\) by
decay at infinity. It is continuous at \(r_y\), and integration of the
radial equation across \(r_y\) gives the derivative jump
\(-e^{-in\phi_y}/(2\pi r_y)\). The Wronskian
\(I_n(t)K_n'(t)-I_n'(t)K_n(t)=-1/t\) gives (CI4). The formula at
\(r_y=0\) follows by continuity. [NIST DLMF, modified Bessel equations and Wronskians](https://dlmf.nist.gov/10.28.E2).

Take \(r\) larger than the distance of the support of \(\mu\) from
\(y_0\). The kernel in (CI4) is then continuous on the circle times
that support, so Fubini and (CI3) compute every Fourier coefficient
of the potential. It follows that

\[
 W=G_\lambda*\mu-CG_\lambda(\,\cdot-y_0)=0
 \quad\text{outside a sufficiently large disk}.
\tag{CI5}
\]

On \(E\setminus\{y_0\}\) this difference solves the homogeneous
modified Helmholtz equation and is real analytic. The punctured open
connected set \(E\setminus\{y_0\}\) is connected, so (CI5) extends to
it. The convolution \(G_\lambda*\mu\) is locally bounded because
\(\mu\) has bounded compactly supported density and \(G_\lambda\)
has only a locally integrable logarithmic singularity. If \(y_0\)
were in \(\overline E\), the equality from \(E\setminus\{y_0\}\)
would contradict the divergence of \(CG_\lambda(\cdot-y_0)\) at
\(y_0\). Hence \(y_0\notin\overline E\) and \(W=0\) on \(E\).

## Second derivatives at the junction

Work in a disk about \(x_*\) avoiding \(y_0\). The potential \(W\) is
\(C^1\) there. Indeed both \(G_\lambda\) and its first derivatives
are locally integrable, and convolutions against a bounded compactly
supported density are continuous by continuity of translation in
local \(L^1\). In the sector of \(A_j\) for \(j=1,2\), its equation is

\[
 (-\Delta+\lambda^2)W=a_j,
 \qquad a_j(x)=m_j e^{-|x-y_0|^2/(2\sigma^2)},
 \qquad (m_1,m_2)=(2,1).
\tag{CI6}
\]

Since \(W=0\) on \(E\), its value and full gradient vanish on each
exterior arc. The analytic \(j/3\) arc extends past the junction
as an analytic curve by the implicit function theorem for the
corresponding pair score. Apply the Cauchy–Kowalevski theorem to
(CI6) with the analytic forcing \(a_j\) and zero value and normal
derivative on this curve. The curve is noncharacteristic because
the principal symbol of \(-\Delta\) is nonzero on real nonzero
normals. This gives an analytic solution \(V_j\) in a neighborhood
of \(x_*\).

At any point of the open \(j/3\) arc sufficiently close to \(x_*\),
the difference \(W-V_j\) solves the homogeneous equation on the
\(A_j\) side and has zero value and normal derivative on the arc.
Extend that difference by zero to the other side in a small disk
about the point. The zero traces cancel the boundary terms in the
distributional equation, so the extension solves the homogeneous
elliptic equation in the disk. Analyticity and vanishing on the
other side imply \(W=V_j\) on an open subset of \(A_j\). Analytic
unique continuation in the connected \(A_j\) sector gives this
equality throughout the sector near \(x_*\).

Use signed Euclidean distance \(d_j(x)\) to the extended analytic
\(j/3\) curve and its orthogonal projection \(q_j(x)\) in a small
tubular neighborhood of \(x_*\). Since \(V_j=\nabla V_j=0\) on this
curve, its Hessian there annihilates the tangent and (CI6) makes its
normal second derivative \(-a_j\). Taylor expansion in signed
distance gives, uniformly as \(x\to x_*\),

\[
 \begin{aligned}
 V_j(x)&=-\tfrac12 a_j(q_j(x))d_j(x)^2+O(|d_j(x)|^3),\\
 |\nabla V_j(x)|&=a_j(q_j(x))|d_j(x)|+O(d_j(x)^2).
 \end{aligned}
\tag{CI7}
\]

At points of the open \(1/2\) arc near \(x_*\), the third adjusted
score is strictly smaller than the other two. Thus these points
lie on neither extended exterior curve, and both \(d_j\) are
nonzero. The common \(C^1\) potential gives \(V_1=V_2\) and
\(\nabla V_1=\nabla V_2\) along this interior arc. In particular the
quotient \(|\nabla V_j|^2/(-2V_j)\) is the same for \(j=1,2\) there.
By (CI7), its limit at \(x_*\) is \(a_j(x_*)\) for each \(j\).
Since \(a_1(x_*)=2a_2(x_*)>0\), the limits differ. This contradiction
proves the theorem, including when the three arcs have a common
tangent.

## Algebraic carriers and the conic consequence

Under the same boundedness assumption, suppose an open \(i/3\) arc
lies in an irreducible complex algebraic curve. Then its polynomial
has even degree \(2m\), and its top homogeneous part is a nonzero
multiple of \((z_1^2+z_2^2)^m\). Consequently no \(i/3\) arc of the
theorem lies in any conic.

To prove the carrier restriction, use the same \(f,c,\sigma\) as in
(CI1), interchanging the bounded labels as needed. The equation on
the open arc can be written

\[
 L(z)=2\pi\sigma^2 c\,
           e^{\rho^2(z_1^2+z_2^2)/(2\sigma^2)},\qquad
 L(z)=\int f(y)e^{\rho z\cdot y/\sigma^2-|y|^2/(2\sigma^2)}\,dy.
\tag{CI8}
\]

Both sides are entire. By the identity theorem, the identity extends
to the connected normalization of the affine curve. The compact
support of \(f\) gives \(|L(z)|\le C_0e^{C_1|\operatorname{Re}z|}\)
with finite constants \(C_0,C_1\). At any point at infinity on the
projective normalization, let \(N\) be the largest pole order of
\(z_1,z_2\). If \(z_1^2+z_2^2\) had pole order \(M>N\), choose
a ray in a local complex parameter for which its leading term is
positive real. The right side of (CI8) would grow as
\(\exp(c_0|t|^{-M})\) or faster for some \(c_0>0\), contradicting the
bound \(|L|\le C_0\exp(C_2|t|^{-N})\) on that ray.

Put \(u=z_1+i z_2\), \(v=z_1-i z_2\). Their largest pole order
is still \(N\), while \(z_1^2+z_2^2=uv\). Thus \(u,v\) cannot both
have a pole at the same point of the projective normalization.
In particular all the values of \(v\) on the curve with a given
large \(|u|\) are uniformly bounded. Otherwise a sequence with
both \(|u|\) and \(|v|\) diverging would accumulate on the compact
normalization at a common pole. The same holds with \(u,v\)
interchanged.

Write the irreducible defining polynomial as
\(P(u,v)=\sum_{k=0}^{m}a_k(u)v^k\). Both coordinates are nonconstant
on the curve since the real arc is nonconstant. Vieta's formulas
and the uniform bound on the roots at large \(|u|\) give that
each \(a_k(u)/a_m(u)\) is bounded at infinity. Therefore
\(\deg a_k\le\deg a_m=:n\), so the bidegree corner \(u^nv^m\)
has nonzero coefficient. The original curve contains an open
real arc. It and its conjugate in \(z_1,z_2\) thus have infinitely
many common points and coincide by irreducibility and Bézout's
theorem. Conjugation interchanges the roles of \(u,v\), so \(n=m\).
The unique term of total degree \(2m\) is a nonzero multiple
of \(u^mv^m=(z_1^2+z_2^2)^m\).

A real conic carrying an open real arc must therefore be an
irreducible quadratic with isotropic quadratic part. Completing
the square shows that it is a circle of positive radius.
A reducible conic can only carry a line arc, which the degree
restriction excludes. The circular obstruction excludes the
remaining conics.

## Algebraic arcs with two noncollinear junctions

**Theorem.** Suppose \(A_1,A_2\) are bounded and the unbounded strict
winner component of \(A_3\) is adjacent to two distinct trivalent
triple junctions. Suppose the three tangent lines at each of those
junctions are distinct. Then no open \(i/3\) arc is contained in a
real algebraic curve, for \(i=1,2\).

Assume such an arc exists and use the irreducible carrier and
notation of (CI8). There is a point at infinity on its projective
normalization where \(u\) has a pole, since \(u\) is nonconstant.
At this point \(v\) is holomorphic by the pole-order argument. We
can take a local parameter \(t\) with \(u=t^{-d}\) for a positive
integer \(d\), by absorbing an analytic root of the leading unit
of \(u\) into the parameter. Write \(v=v(t)\), \(v_0=v(0)\), and
let \(v_l(u)\), \(0\le l<d\), be the values from the \(d\) choices
of \(t\) for large \(|u|\).

Change variables in \(L\) to
\(\zeta=\rho(y_1-i y_2)/(2\sigma^2)-\zeta_0\), choosing \(\zeta_0\)
so that the origin lies in a bounded strict winner. Let \(d\nu\)
be the image of \(f(y)e^{-|y|^2/(2\sigma^2)}dy\) under this map,
and let \(\mathcal E\) be the image of \(E\). The measure \(\nu\)
has compact support and bounded density. Near either junction its
density is a common positive analytic factor times \(2,1,0\)
in the sectors of \(A_1,A_2,E\), respectively. Define

\[
 \alpha=\frac{\rho^2}{2\sigma^2},\qquad
 b=\alpha v_0-\zeta_0,\qquad
 A(u,\zeta)=\frac1d\sum_{l=0}^{d-1} e^{v_l(u)\bar\zeta}.
\]

Divide (CI8) by \(e^{u\zeta_0+v\bar\zeta_0}\) and average over
the \(d\) values of \(v\). The result is

\[
 \int e^{u\zeta} A(u,\zeta)\,d\nu(\zeta)
       =2\pi\sigma^2 c\,e^{bu}H(u),\qquad
 H(u)=\frac1d\sum_l e^{\alpha u(v_l(u)-v_0)-v_l(u)\bar\zeta_0}.
\tag{CI9}
\]

Both averages are single-valued for large \(|u|\), since circling
infinity permutes the \(d\) terms. There are expansions and bounds

\[
 A(u,\zeta)=\sum_{j\ge0} a_j(\zeta)u^{-j},\qquad
 a_0(\zeta)=e^{v_0\bar\zeta},\qquad
 \sup_{\operatorname{supp}\nu}|a_j|\le C_3 C_4^j,\qquad
 |H(u)|\le C_\epsilon e^{\epsilon |u|}
 \quad(\epsilon>0).
\tag{CI10}
\]

The coefficient bound is Cauchy's estimate in \(1/u\) on a fixed
small circle, uniformly on the compact support of \(\nu\).
For the last bound use \(v_l(u)-v_0=O(|u|^{-1/d})\).

### Exterior Borel transforms

For a Laurent series at infinity, apply the transform
\(u^n\mapsto n!/s^{n+1}\) to its nonnegative powers and discard
its negative powers. For \(j\ge1\), the transform of
\(u^{-j}e^{u\zeta}\) is the function

\[
 H_j(s,\zeta)=\frac1{(j-1)!}
       \int_0^\zeta \frac{(\zeta-t)^{j-1}}{s-t}\,dt
 \qquad (|s|>|\zeta|).
\tag{CI11}
\]

Here the integral is along the straight segment. Expand
\(1/(s-t)\) and integrate each power of \(t\) to get the
coefficient \(\zeta^{n+j}n!/(n+j)!\) of \(s^{-n-1}\).
The transform for \(j=0\) is \(1/(s-\zeta)\). Thus the left
transform in (CI9) is

\[
 \int \frac{a_0(\zeta)}{s-\zeta}\,d\nu(\zeta)
    +\sum_{j\ge1}\int a_j(\zeta) H_j(s,\zeta)\,d\nu(\zeta).
\tag{CI12}
\]

The series and the transform agree for sufficiently large \(|s|\)
by (CI10), compact support, and the factorial in (CI11).

The right transform continues along every path avoiding \(0,b\).
In fact, write \(H(u)=\sum_{n\ge0} h_nu^n+\sum_{j\ge1} h_{-j}u^{-j}\).
Cauchy's estimate with radius \(n/\epsilon\) and Stirling's formula
give \(\limsup_{n\to\infty}(n!|h_n|)^{1/n}=0\). The nonnegative part
times \(e^{bu}\) therefore has transform
\(\sum_{n\ge0} n!h_n/(s-b)^{n+1}\), holomorphic off \(b\). The shift
formula follows by expanding \(e^{bu}\) and summing the resulting
binomial series for large \(|s|\). The negative part times \(e^{bu}\)
has transform \(\sum_{j\ge1} h_{-j} H_j(s,b)\), where
\(|h_{-j}|\) grows at most geometrically by Cauchy's estimate on a
fixed large circle. In (CI11), expanding
\(\zeta-t=(\zeta-s)+(s-t)\) shows that \(H_j\) is the sum of a
polynomial in \(s,\zeta\) and
\((\zeta-s)^{j-1}\log(s/(s-\zeta))/(j-1)!\). This gives
continuation along the stated paths. On compact parts of any such
continuation, the polynomial and logarithm formula bounds the
series terms by \(C_5 C_6^j/(j-1)!\). Hence the series converges
locally uniformly and is holomorphic there.

### The derivative at a junction

Choose one of the two junctions whose image \(s_*\) differs from
\(b\). Also \(s_*\ne0\) by the choice of origin. Choose a ray
\(s=s_*+\delta e^{i\psi}\) pointing strictly inside the exterior
tangent sector. For small \(\delta>0\) it lies in \(\mathcal E\).
Continue (CI12) from a point far outside the support along a
fixed piecewise smooth path in \(\mathcal E\) avoiding \(b\) to
this ray, then along the ray toward \(s_*\). Such a path exists
since the unbounded open component stays connected after removing
one point. The fixed initial path is separated from the compact
support and from \(0\). The continued logarithm
\(\log(s/(s-\zeta))\) for each \(\zeta\in\operatorname{supp}\nu\)
therefore has uniformly bounded imaginary part along the fixed
path. On the final straight segment the arguments of \(s\) and
of \(s-\zeta\) each vary by at most \(\pi\). Its real part is bounded
in absolute value by \(C_7+|\log|s-\zeta||\) there. These bounds
justify termwise integration and continuation on a neighborhood
of every point of the chosen path in \(\mathcal E\).

In (CI12) the derivative of the sum over \(j\ge1\) stays bounded
as \(\delta\downarrow0\). Indeed

\[
 \partial_s H_1(s,\zeta)=1/s-1/(s-\zeta),\qquad
 \partial_s H_j(s,\zeta)=
       \frac{\zeta^{j-1}}{(j-1)!\,s}-H_{j-1}(s,\zeta)\quad(j\ge2).
\tag{CI13}
\]

These formulas follow by integration by parts in (CI11) and then
continue analytically. Since \(\nu\) has bounded density and compact
support, the integrals of \(|s-\zeta|^{-1}\) and
\(|\log|s-\zeta||\) are uniformly bounded for \(s\) near \(s_*\).
The polynomial and logarithm formula for \(H_{j-1}\) bounds its
integral in absolute value by \(C_8 C_9^j/(j-2)!\). Combining
this with (CI10) gives a convergent bound for the derivative sum.

The derivative of the first term in (CI12) has a logarithmic
divergence. In polar coordinates \(\zeta=s_*+r e^{i\theta}\)
let \(m(\theta)\) be \(2,1,0\) in the tangent sectors of
\(A_1,A_2,E\), respectively. The contribution from
\(|\zeta-s_*|\le2\delta\) is bounded, because the chosen ray points
strictly inside the exterior sector and its points are at distance
at least \(c_1\delta\) from the support for small \(\delta\), with
\(c_1>0\). On \(2\delta<r<r_0\), replace \((s-\zeta)^{-2}\) by
\((s_*-\zeta)^{-2}\). The integrand error before the area factor
is \(O(\delta/r^3)\) and gives \(O(1)\) in the integral. The
analytic density factor differs from its junction value by
\(O(r)\), as does \(a_0\). Each analytic interface meets the
circle of radius \(r\) at an angle within \(O(r)\) of its tangent
ray. These errors are integrable against \(dr/r\).
Therefore the derivative has the expansion

\[
 -q_*\,\log(1/\delta)\int_0^{2\pi}m(\theta)e^{-2i\theta}\,d\theta
   +O(1),
\tag{CI14}
\]

where \(q_*\ne0\) is the common analytic density factor at
\(s_*\) times \(a_0(s_*)\). In particular the coefficient is
independent of the chosen exterior ray.

The angular integral in (CI14) is nonzero. Integrate over the
three sectors. The jumps of the piecewise constant function
\(m\) are, up to a common sign, \(2,-1,-1\). If the integral
vanished, these jumps would give
\(2e^{-2i\theta_a}=e^{-2i\theta_b}+e^{-2i\theta_c}\)
for the three boundary ray angles. The average of two distinct
unit complex numbers has modulus strictly less than one. The
distinct tangent lines ensure that the three unit complex numbers
here are distinct.

The right transform of (CI9) is holomorphic at the end of the
chosen continuation since \(s_*\) differs from \(0,b\). Its
derivative is bounded there. The left transform has unbounded
derivative by (CI13)–(CI14), contradicting equality after
continuation.

## Algebraic arcs with two trivalent exterior junctions

**Theorem.** Two trivalent exterior junctions suffice to exclude
algebraic outer arcs. Suppose \(A_1,A_2\) are bounded and the
unbounded strict winner component of \(A_3\) is adjacent to two
distinct trivalent triple junctions. Then no open \(i/3\) arc is
contained in a real algebraic curve, for \(i=1,2\).

Under the algebraic arc assumption, choose \(s_*\ne b\) as in
the preceding proof. If this junction has three distinct tangent
lines, the argument in (CI13)–(CI14) applies directly. In the
remaining trivalent case the local winner classification in
(G20) gives three regular analytic curve germs with a common
tangent. On the side of the junction where the middle label wins,
the nonzero function \(q(s)\) in the local classification puts its
two bounding curves on opposite sides of the outer-pair score
equality curve. This distinguishes all three germs. Denote these analytic
germs by \(\Gamma_a,\Gamma_b,\Gamma_c\)
and their holomorphic Schwarz functions near \(p=s_*\) by
\(S_a,S_b,S_c\). Thus \(S_\ell(\zeta)=\bar\zeta\) along
\(\Gamma_\ell\). Such a function exists by extending a regular
real-analytic parametrization \(\chi_\ell\) holomorphically and
setting \(S_\ell=\overline{\chi_\ell}\circ\chi_\ell^{-1}\),
where the overline on a holomorphic function conjugates its
coefficients. In particular both compositions
\(S_\ell\circ\overline{S_\ell}\) and
\(\overline{S_\ell}\circ S_\ell\) are the identity as germs.

Take a small disk centered at \(p\) avoiding \(0,b\), whose
boundary meets each of the three actual arcs once. Continue the
left side of (CI12) from infinity along a fixed piecewise smooth
path in \(\mathcal E\setminus\{b\}\) to a point in a smaller
concentric disk. The continuation of \(\log(s-\zeta)\), viewed
as a function of \(\zeta\), is holomorphic on a neighborhood
of the two local bounded label regions and their boundaries
for \(s\) near that endpoint. Indeed the chosen path to the
endpoint stays a positive distance from the compact support.
Starting with the logarithm at infinity, integrating
\(ds/(s-\zeta)\) along the path produces these holomorphic germs
uniformly in \(\zeta\) near the support. Contributions to
(CI12) from outside the disk continue holomorphically across
\(p\), since every such \(\zeta\) stays a fixed distance from
the smaller disk.

For the local contributions, write the common analytic density
factor as \(\eta(\zeta,\bar\zeta)\). Let \(B_j(z,w)\) be the
holomorphic primitive in \(w\), based at \(\bar p\), of
\(a_j(w)\eta(z,w)\). Here the holomorphic extension of \(a_j\)
in its argument \(w=\bar\zeta\) is supplied by the coefficient
formula in (CI10). The primitives have a uniform bound
\(C_{10}C_{11}^j\) on a common bidisk by the same Cauchy estimate
used in (CI10). Apply complex Green's formula separately on the
two bounded label regions in the disk to the product of \(B_j\)
with the holomorphic kernel \(H_j\) in the integration variable
\(\zeta\), using \(H_0(s,\zeta)=1/(s-\zeta)\). The boundary
densities on the three arcs are \(B_j(\zeta,S_\ell(\zeta))\).
With all arcs oriented away from \(p\), the signed jumps of the
weights \(2,1,0\) are, up to common sign, \(2,-1,-1\).
Relabel the arcs so that their combination is

\[
 F_j(z)=2B_j(z,S_a(z))-B_j(z,S_b(z))-B_j(z,S_c(z)).
\tag{CI15}
\]

The local disk's circular boundary contributes holomorphic terms
near \(p\) because the kernel arguments there stay separated.
On an arc directed from \(p\) to \(q\), the \(j=0\) integral of
an analytic density \(g(\zeta)\) against \(1/(s-\zeta)\) equals
\(g(s)\log(s-p)\) plus a holomorphic germ at \(p\). This follows
by subtracting \(g(s)\) from \(g(\zeta)\) in the integrand and
integrating the constant part. For \(j\ge1\) use the polynomial
and logarithm formula for \(H_j\) after (CI12). Its only term
that can have a branch point at \(p\) is the integral against
\(-(\zeta-s)^{j-1}\log(s-\zeta)/(j-1)!\). Define

\[
 P_j(s,\zeta)=\frac1{(j-1)!}
               \int_s^\zeta (t-s)^{j-1}g(t)\,dt .
\]

Integrating by parts produces the coefficient
\(P_j(s,p)\) of \(\log(s-p)\). The remaining integral is
holomorphic near \(p\), since \(P_j(s,\zeta)/(s-\zeta)\) has
a removable singularity at \(\zeta=s\), and the other endpoint
\(q\) stays separated from \(p\). These formulas apply first near
the chosen exterior endpoint and then as analytic continuations.
Together they express the left side of (CI12) near \(p\) as a
holomorphic germ plus a constant nonzero factor times

\[
 \left[
 F_0(s)+\sum_{j\ge1}\frac{(-1)^j}{(j-1)!}
       \int_p^s(s-t)^{j-1}F_j(t)\,dt
 \right]\log(s-p).
\tag{CI16}
\]

The factor is the common coefficient from Green's formula and
the orientation choice. All the holomorphic remainders, for
\(j\ge1\), are bounded on a sufficiently small disk by
\(C_{12}C_{13}^j/(j-1)!\). For the polynomial part this follows
by expanding \(\zeta-t=(\zeta-s)+(s-t)\) in (CI11), and for the
integrated part it follows by bounding \(P_j\) and
\(P_j/(s-\zeta)\) on the fixed local chart. The logarithms at
the outer boundary and at the fixed endpoint \(q\) have
holomorphic branches uniformly bounded on the smaller disk.
Thus the holomorphic remainders and the coefficient series
converge locally uniformly there.

It remains to find the first nonzero jet in the bracket of
(CI16). Flatten \(\Gamma_a\) by a holomorphic parametrization
\(\chi\) with \(\chi(0)=p\), and write
\(\widetilde S_\ell=\bar\chi^{-1}\circ S_\ell\circ\chi\).
Then \(\widetilde S_a(z)=z\), and the other two germs have
derivative \(1\) at \(0\) because the tangent is common.
Let \(k\ge2\) be the smaller of the vanishing orders of
\(\widetilde S_b-z\) and \(\widetilde S_c-z\).
The sum \(2z-\widetilde S_b(z)-\widetilde S_c(z)\) has a
nonzero jet of order \(\ell\le 2k-1\). Indeed the order-\(k\)
coefficient already works unless the two leading coefficients
are opposite, say \(a\) and \(-a\), with \(a\ne0\).
For either germ, its composition with its coefficientwise
conjugate is the identity. At degree \(2k-1\), with \(d_\ell\)
the coefficient of the individual germ at that degree, this says
\(d_\ell+\bar d_\ell+k|a|^2=0\). The real parts of the
degree-\(2k-1\) coefficients in the two germs therefore agree
and are nonzero. They cannot cancel in the sum.

Undoing the coordinate change gives the same vanishing order
\(\ell\) for \(2S_a-S_b-S_c\) at \(p\). In fact, Taylor expansion
of \(\bar\chi\) about \(z\) multiplies the first nonzero jet
of the sum by \(\bar\chi'(0)\), since all terms quadratic in
\(\widetilde S_b-z,\widetilde S_c-z\) have order at least \(2k\).
The same expansion in the second argument of \(B_j\) about
\(S_a(z)\) shows that every \(F_j\) vanishes to order at least
\(\ell\), and \(F_0\) vanishes to order exactly \(\ell\).
For the last claim, \(\partial_wB_0(p,\bar p)=
e^{v_0\bar p}\eta(p,\bar p)\ne0\).

Each \(j\)-fold primitive in (CI16) vanishes to order at least
\(\ell+j\). Hence the bracket has a nonzero order-\(\ell\) jet.
Its product with \(\log(s-p)\) cannot equal the right transform
of (CI9), which is holomorphic across \(p\). This proves the
algebraic obstruction also at the common-tangent junction.

## Connected semialgebraic maximizers

With the hypotheses of
[connected three-cell networks with algebraic ends](THREE_CELL_SEMIALGEBRAIC_NETWORK_RIGIDITY.md),
every balanced planar global maximizer is the centered regular
simplex when the jump network is semialgebraic. The three-unbounded
and two-unbounded cases are covered by that result. In the remaining
case, the network is a compact theta with two trivalent exterior
junctions. A semialgebraic jump network has empty interior. In any
Boolean description by finitely many polynomial signs, at least one
nonzero defining polynomial vanishes at every point. A point with
all defining signs constant in a neighborhood would be interior.
The product of those polynomials thus vanishes on its exterior
arcs, contradicting the algebraic obstruction. This leaves open
maximizers without connected strict winners or with jump networks
outside the stated semialgebraic class.

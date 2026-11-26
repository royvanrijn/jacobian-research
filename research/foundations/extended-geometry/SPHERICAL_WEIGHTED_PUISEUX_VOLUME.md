# Root contacts in the weighted proper-volume criterion

The multiplicity of a root of the leading angular polynomial gives
a sufficient volume test. The rates at which its branches separate
decide the exact local test. This leads to nonvanishing results even
at the maximum possible root multiplicity in ordinary sphere degree.

Use the literal chart of the sphere \(2xy+t^2=1\) and a rational
angular weight \(\kappa\). For a nonconstant regular phase \(f\), put
\(n=\max(j+\kappa i)\) over its support \(x^it^j\), and assume \(n>0\).
Take a positive integer \(d\) such that \(d\kappa\) is integral. In the
chart \(t=v^{-d}\), \(z=xv^{d\kappa}\), write
\[
 f=v^{-dn}H(v,z),\qquad Q(z)=H(0,z),\qquad
 \omega=-d v^{-d-1}\frac{dz}{z}\wedge dv .
\]
Here \(\omega\) is a nonzero constant multiple of the sphere volume
form, and \(H\) is a polynomial in \(v\) with Laurent-polynomial
coefficients in \(z\).

## The local contact criterion

Let \(z_0\ne0\) be a root of \(Q\) of multiplicity \(r\). In a disk
about \(z_0\), take a finite further ramification \(v=u^e\) so that the
\(r\) nearby roots are holomorphic functions \(z_0+\xi_i(u)\), repeated
with multiplicity. Put \(D=de\) and
\(\nu_{ij}=\operatorname{ord}_u(\xi_i-\xi_j)\), including infinity
for identical functions. In particular \(\nu_{ii}=\infty\).

**Theorem.** The squared volume of \(|f|\le M\) near \((v,z)=(0,z_0)\)
is finite for any fixed \(M>0\) exactly when
\[
             \sum_{j=1}^r\min(D,\nu_{ij})<Dn
                    \quad\text{for every }i.                 \tag{1}
\]
When the largest sum equals \(Dn\), the divergence is logarithmic
after truncating radially in the ramified parameter \(u\).
The test depends only on the contacts \(\operatorname{ord}_v\) obtained
by dividing \(\nu_{ij}\) by \(e\), so a larger common ramification
does not affect it.

Weierstrass preparation gives the local monic degree-\(r\) factor
whose roots tend to \(z_0\). Its squarefree part over the field of
meromorphic germs has nonzero discriminant. After shrinking the
punctured disk, its distinct roots define a finite unramified cover
of that disk. Its monodromy is a permutation of finitely many roots
and becomes trivial after \(v=u^e\) for some \(e\). The roots are
bounded and extend holomorphically across \(u=0\). The multiplicity
of each branch is constant on the punctured disk, and the unit in
Weierstrass preparation has a uniform positive lower bound and
finite upper bound on a small closed bidisk.
[Riemenschneider, Weierstrass preparation](https://www.math.uni-hamburg.de/home/riemenschneider/bookchap03.pdf).

Set \(w=z-z_0\) and, for fixed small nonzero \(u\), choose any closest
root \(\xi_i\) to \(w\). If \(\rho=|w-\xi_i|\), the triangle inequality
and the choice of \(i\) give
\[
 \tfrac12\max(\rho,|\xi_i-\xi_j|)
 \le |w-\xi_j|\le2\max(\rho,|\xi_i-\xi_j|).
\]
Indeed a separation of at least \(2\rho\) supplies the lower bound
by the triangle inequality. A smaller separation is handled by
\(|w-\xi_j|\ge\rho\). Each nonzero separation is comparable to
\(|u|^{\nu_{ij}}\) uniformly in \(\arg u\). Put
\(A_i(a)=\sum_j\min(a,\nu_{ij})\). Since at least the self term
has infinite contact, \(A_i\) increases continuously and strictly
from zero to infinity. Define \(\alpha_i>0\) by
\(A_i(\alpha_i)=Dn\). In the region with nearest root \(i\), the
product of distances is comparable to
\(\prod_j\max(\rho,|u|^{\nu_{ij}})\), with an infinite-contact
power interpreted as zero. At \(\rho=|u|^{\alpha_i}\) this product
is \(|u|^{Dn}\), and increasing the radius by a factor \(R\ge1\)
increases it at least by \(R\) through the self factor. Thus the
sublevel region is covered by disks of radii
\(O_M(|u|^{\alpha_i})\).

Conversely, in any disk \(|w-\xi_i(u)|<c_M|u|^{\alpha_i}\), the
triangle inequality gives the upper bounds by the same maxima,
and the factor \(|w-\xi_i(u)|\) supplies an extra factor \(c_M\).
For small enough \(c_M>0\) the disk lies in the sublevel region.
Therefore the planar area at fixed \(u\) is comparable to
\(|u|^{2\min_i\alpha_i}\). The pulled-back volume form has squared
density comparable to \(|u|^{-2D-2}d^2u\,d^2w\) because \(z_0\ne0\).
Integration gives exactly the condition \(\alpha_i>D\) for every
\(i\), which is (1). At the threshold the radial integral is
logarithmic. The cover away from \(u=0\) is finite, so finiteness
is equivalent in the original weighted chart.

If \(r<n\), condition (1) follows immediately. At \(r=n\), it is
equivalent to having a pair of branches with \(\nu_{ij}<D\).
Indeed the condition that two branches have contact at least
\(D\) is an equivalence relation. If it has more than one class,
every branch has a contact below \(D\) to another branch. At
\(r>n\), the same formula still tests the actual branch separation.

## Consequences for sphere moments

**Corollary.** Suppose the leading angular support has extremes
\(-b,a\) with \(a,b>0\), the full extremes are \(-B,A\), and
\(A-a<n\), \(B-b<n\). If (1) holds at every finite nonzero root
of \(Q\), then infinitely many \(\Lambda(f^m)\) are nonzero.

Use the finite-height and endpoint estimates of the
[weighted proper-volume theorem](SPHERICAL_WEIGHTED_FACE_PROPER_VOLUME_OBSTRUCTION.md#the-bilateral-weighted-criterion).
The new theorem supplies all the remaining finite-angular charts
at infinite height. The same resolution and proper-period argument
then applies. For a characteristic-zero coefficient field, use
an embedding of the finitely generated coefficient field into
\(\mathbb C\). The needed contacts are Puiseux valuations of the
algebraic polynomial \(z^B H(v,z)\) at its roots, so they are
preserved under any such embedding, as are the rational moment identities.

In particular, take \(\kappa=1\) and a crossing leading face of
integer degree \(n\ge2\). The full angular extremes satisfy
\(A,B\le n\) by nonnegative heights and by regularity at negative
angles. The endpoint inequalities therefore hold. For fixed \(f\),
replace it by \(f+\lambda t^{n-1}\). This changes the coefficient
of \(v\) in \(H\) by the constant \(\lambda\). At each nonzero
root \(z_0\), exclude the one value of \(\lambda\) which makes
that coefficient vanish at \(z_0\). If the root multiplicity
is \(r\), the substitution \(v=u^r\), \(z-z_0=u y\), followed by
division by \(u^r\), gives at \(u=0\) a polynomial \(c y^r+c'\)
with \(cc'\ne0\). Its roots are simple, so all distinct branch
contacts in \(v\)-units are \(1/r\). For each branch the sum in
(1), after division by \(D\), is \(1+(r-1)/r<2\le n\). Thus
every parameter outside the explicitly given finite set produces
infinitely many nonzero moments, with the same leading face.

## Regular phases with an angular root of multiplicity \(2n\)

For any \(n\ge2\) and \(1\le B\le n\), take
\(Q(z)=z^{-B}(z-1)^{n+B}=\sum_i q_i z^i\). Lift a term with
\(i\ge0\) to \(q_i x^it^{n-i}\), and a term with \(i=-j<0\) to
\(q_i x^{-j}(t^2-1)^jt^{n-j}\). The latter is a regular monomial
in \(y,t\), because \((t^2-1)/x=-2y\) on the sphere. All deficits
from the face are even. Adding \(t^{n-1}\) gives first correction
\(v\) and supplies the preceding simple-branch argument. These
integer phases are absolutely nonpower, since a proper power
would have multiple Puiseux branches wherever its reduced root
occurs. With \(B=n\), the multiplicity \(2n\) attains the full
angular width bound at weight one.

Even the unperturbed lift works for \(n\ge3\). Its first possible
correction has deficit two, with value at \(z=1\)
\[
 -\sum_{k=0}^{B-1}(B-k)(-1)^{n+B-k}\binom{n+B}{k}
       =(-1)^n\binom{n+B-2}{B-1}\ne0 .
\]
The equality is coefficient extraction from
\(z(1-z)^{n+B}/(1-z)^2\). Use \(v=u^r\), \(z-1=u^2 y\) for
\(r=n+B\). The branch contacts in \(v\)-units are now \(2/r\), and
the test sum divided by \(D\) is \(1+2(r-1)/r<3\le n\).
This retains the parity of the height deficits. The
[critical prime-power examples](SPHERICAL_CRITICAL_FACE_P_POWER_MOMENTS.md)
instead have a cluster at or beyond the cutoff and display divergent
volume at critical multiplicity. Their interval-pole certificates
address precisely a regime that the local volume test can miss.

## Angular endpoints

At an angular endpoint let \(w=z\) or \(w=z^{-1}\) tend to zero.
Let \(C\ge0\) be the pole order of the full angular support in this
direction, so \(G(v,w)=w^C H(v,z)\) is holomorphic. Use the
actual endpoint exponent, so \(G(v,0)\) is not identically zero.
If \(G(0,0)=0\), take its local branches \(\xi_i(u)\) after
\(v=u^e\) as above, and put \(\mu_i=\operatorname{ord}_u\xi_i\).
All \(\mu_i\) are finite. Write \(D=de\) and
\(A_i(a)=\sum_j\min(a,\nu_{ij})\) as before.

**Proposition.** At this endpoint the squared-volume test is exactly
\[
                 A_i(\mu_i+D)<Dn+C\mu_i
                 \quad\text{for every }i .                   \tag{2}
\]
When \(G(0,0)\ne0\), the bounded-phase region avoids the corner.

Let \(\beta_i\) be the unique positive solution of
\(A_i(\beta_i)=Dn+C\mu_i\). In a region where the nearest-root
distance \(\rho=|w-\xi_i(u)|\) is below \(|w|/2\), one has
\(|w|\asymp|\xi_i(u)|\asymp|u|^{\mu_i}\). The bounded-phase
inequality is \(|G|\le M|u|^{Dn}|w|^C\). The same product estimate
as at a finite nonzero root bounds \(\rho\) by
\(O_M(|u|^{\beta_i})\). If \(\beta_i>\mu_i\), the disks of this
radius contribute at most a constant times
\(|u|^{2(\beta_i-\mu_i)}\) to the \(w\)-area weighted by \(|w|^{-2}\).
A disk of radius \(c_M|u|^{\beta_i}\) about \(\xi_i\) supplies
the matching lower bound. If \(\beta_i\le\mu_i\), then
\(A_i(\mu_i)\ge Dn+C\mu_i\). A disk of radius
\(c_M|u|^{\mu_i}\) is in the sublevel region by the same product
upper bound, and has weighted \(w\)-area bounded below by a
positive constant. Either lower bound implies divergence in \(u\)
whenever \(\beta_i-\mu_i\le D\).

It remains to bound the points outside the regions with closest
distance below \(|w|/2\). There each \(|w-\xi_i(u)|\) is comparable to
\(\max(|w|,|u|^{\mu_i})\). For \(q=|w|\), the function
\[
                     q^{-C}\prod_i\max(q,|u|^{\mu_i})
\]
is a power times a constant between successive breakpoints
\(q=|u|^{\mu_i}\). Its minimum is attained at a breakpoint or
at the fixed outer radius. On the innermost interval it is
decreasing or constant as \(q\) increases because \(C\ge0\).
At each breakpoint the exponent of \(|u|\) is
\(A_i(\mu_i)-C\mu_i\), since
\(\min(\mu_i,\nu_{ij})=\min(\mu_i,\mu_j)\).
If all \(\beta_i>\mu_i\), these exponents are strictly below
\(Dn\). The bound \(|G|/|w|^C\le M|u|^{Dn}\) is then impossible
outside the nearest-root regions for small \(u\). The remaining
\(u\)-density is proportional to \(|u|^{-2D-2}d^2u\), so the
upper bounds give convergence exactly when \(\beta_i-\mu_i>D\)
for all \(i\). This is (2).

Consequently, when the full angular support has both signs,
(1) at every nonzero root together with (2) at both endpoints
is necessary and sufficient for finiteness of the squared volume
of each bounded-phase region. The finite-height charts are
already integrable by the weighted proper-volume proof. When
the conditions hold they give infinitely many nonzero moments
even if the leading angular face is axial or one-sided.

For a one-sided example take integers \(n\ge2\), \(1\le L\le n\),
\(1\le B\le n-1\) and
\[
             f=x^L t^{n-L}+y^B t^{n-B-1},\qquad \kappa=1 .
\]
Its top polynomial is \(z^L\). The left endpoint has
\(C=B\) and
\(G(v,z)=z^{B+L}+(-1/2)^Bv(1-v^2)^B\). With \(r=B+L\) and
\(v=u^r\), all branches have \(\mu_i=1\), and all distinct
contacts are one, by their distinct leading coefficients.
Thus \(A_i(1+r)=2r<rn+B\), as required. The right endpoint
has no corner root and there are no finite nonzero roots.
The left loss \(r\) can reach \(2n-1\), beyond the previous
one-sided sufficient bound of loss below \(n\).
At \(n=L=2\), \(B=1\), the phase is \(x^2+y\). Its third
moment is \(2/5\) directly, because the sole angular-constant
term of its cube is \(3x^2y^2=3(1-t^2)^2/4\).

# A proper-volume obstruction from a weighted sphere face

A weighted Laurent face can exclude eventual purity even when the
leading binary form in the usual sphere degree is unstable. The
argument estimates the volume of bounded phase regions at infinity.
At a weighted endpoint, a geometric series over annuli accounts for
the angular logarithmic pole. This resolves the full even-conductor
strip, including arbitrary lower terms and upper-half sources where
the first two arithmetic residues can both vanish.

## The bilateral weighted criterion

On the affine sphere \(2xy+t^2=1\), let \(f\) be a regular phase
over a characteristic-zero field. Use its literal Laurent expression
\(f(x,t)=\sum c_{ij}x^it^j\) from \(y=(1-t^2)/(2x)\), with
\(j\ge0\). Write \(\Lambda\) for angular constant term followed by
the normalized interval integral in \(t\). For a rational
weight \(\kappa\), put
\[
 n=\max_{c_{ij}\ne0}(j+\kappa i),\qquad
 Q(z)=\sum_{j+\kappa i=n}c_{ij}z^i.
\]
Let \(-B,A\) be the extreme angular exponents of the full support,
and \(-b,a\) those of \(Q\).

**Theorem.** Suppose \(n>0\), \(a,b>0\), and
\[
                         A-a<n,\qquad B-b<n.
\]
Suppose also that every nonzero geometric root of \(Q\) has
multiplicity strictly less than \(n\). Then
\(\Lambda(f^m)\ne0\) for infinitely many positive integers \(m\).

### Bounded phase regions at finite height

The finitely many coefficients generate a finitely generated field
over \(\mathbb Q\). Embed this field into \(\mathbb C\). Support,
geometric root multiplicities and moment identities are preserved.
The sphere volume form is a nonzero constant multiple of
\(\omega=(dx/x)\wedge dt\) on the Laurent chart. We show that
its positive squared density has finite integral on \(|f|\le M\)
for each \(M>0\).

We use the following consequence of Weierstrass preparation. If
\(H(v,w)\) is holomorphic near \((0,0)\) and \(H(0,w)\) has
a zero of order \(r>0\), then the planar area of
\(\{w: |w|<\epsilon,\ |H(v,w)|\le\delta\}\) is
\(O(\delta^{2/r})\), uniformly for \(v\) near zero and \(\delta\)
small. On a sufficiently small closed bidisc the unit factor is
bounded away from zero. Hence at least one of the \(r\) linear
factors of the monic polynomial, for each fixed \(v\), must have
magnitude \(O(\delta^{1/r})\) on the sublevel set. The union of
the corresponding disks has the asserted area. If the order is
zero, the function stays bounded away from zero in a neighborhood.
[Riemenschneider, The Preparation Theorem of Weierstraß, Chapter 3](https://www.math.uni-hamburg.de/home/riemenschneider/bookchap03.pdf)

For \(t\) in a compact disk the only ends of the Laurent chart
are \(x=\infty\) and \(x=0\). With \(u=1/x\) at the first and
\(u=x\) at the second, write \(f=u^{-C}H(u,t)\), where \(C=A\)
or \(C=B\), respectively. Both are positive and \(H(0,t)\) is
a nonzero polynomial. Near any zero of this polynomial of order
\(r>0\), the bounded phase region has \(t\)-area
\(O_M(|u|^{2C/r})\). The logarithmic density in \(u\) gives a
convergent integral of \(|u|^{-2+2C/r}\,d^2u\). Away from those
zeros, bounded phase values stay away from \(u=0\). This covers
both ends at bounded \(t\).

### The weighted charts at infinite height

Take a positive integer \(d\) such that \(d\kappa\) is integral.
Over the punctured neighborhood of \(t=\infty\), use the finite
cover \(t=v^{-d}\) and the coordinate \(z=xv^{d\kappa}\). The
cover has degree \(d\) on the Laurent chart, so finiteness of the
pulled-back squared integral implies finiteness below. In these
coordinates
\[
 f=v^{-dn}F(v,z),\qquad F(0,z)=Q(z),\qquad
 \omega=-d\,v^{-d-1}\frac{dz}{z}\wedge dv.
\]
The powers of \(v\) in \(F\) are the nonnegative integers
\(d(n-j-\kappa i)\). At a finite nonzero root of \(Q\) of order
\(r<n\), the sublevel estimate gives \(z\)-area
\(O_M(|v|^{2dn/r})\), and the factor \(|z|^{-2}\) is bounded.
The remaining exponent for \(d^2v\) is
\(-2(d+1)+2dn/r>-2\), so the squared integral is finite there.
On a compact subset away from the roots and angular endpoints,
bounded phase values avoid \(v=0\).

At \(z=\infty\) put \(w=1/z\), \(C=A\), \(c=a\), and
\(G(v,w)=w^A F(v,1/w)\). At \(z=0\) put \(w=z\), \(C=B\),
\(c=b\), and \(G(v,w)=w^B F(v,w)\). In both cases \(G\) is
holomorphic and \(G(0,w)\) has order \(r=C-c\) at zero. If
\(r=0\), bounded phase values avoid the corner since
\(|G(v,w)|\le M|v|^{dn}|w|^C\) there. If \(r>0\), restrict
to a small annulus \(\rho/2<|w|\le\rho\). On the bounded phase
region its \(w\)-area is at most
\(O_M(|v|^{2dn/r}\rho^{2C/r})\) by the uniform sublevel bound.
After multiplication by the logarithmic density \(|w|^{-2}\),
its integral is \(O_M(|v|^{2dn/r}\rho^{2C/r-2})\). The sum
over dyadic annuli converges because \(c>0\) gives \(C>r\).
Finally \(r<n\) makes the remaining integral in \(v\) converge,
as for a finite angular root. These charts exhaust infinite \(t\),
and complete the bounded-phase volume estimate.

### The period conclusion

Resolve the pencil of \(f\) on the smooth projective quadric
compactifying the sphere, with blowups supported on the boundary.
The locus \(Y\) of finite phase values in the resulting smooth
surface is proper over \(\mathbb A^1\) by \(f\) and contains the
original affine sphere. On each relatively compact chart of \(Y\)
the function \(f\) is bounded. The finite squared integral rules
out a pole of the pulled-back volume form along any divisor, since
a transverse pole of order \(k\ge1\) gives a divergent integral
of \(|u|^{-2k}\,d^2u\). The form is therefore regular on \(Y\).
This is the same extension step used in the
[ordinary-degree stability argument](STABLE_LEADING_SPHERE_PHASE_PERIOD_OBSTRUCTION.md)
and the [winding-one cubic argument](WINDING_ONE_CUBIC_PROPER_VOLUME_OBSTRUCTION.md#from-finite-volume-to-nonvanishing).

Apply the [proper-volume lemma](../../../results/primitive-balanced-image-spectrum/proper-volume.tex)
to the compact real sphere cycle inside \(Y\). Normalized integration
there gives \(\Lambda\) on sphere phases. Its exponential period has
constant coefficient one and coefficient \((-1)^m\Lambda(f^m)/m!\)
at order \(m\). Eventual vanishing of the moments would make this
period a nonzero polynomial, excluded by the lemma.

## Axial and one-sided initial faces

**Proposition.** Keep \(A,B>0\) and \(n>0\), and allow either of the
initial endpoint parameters \(a,b\) to be zero or negative. For
each actual endpoint let \(C=A\) or \(C=B\) and \(c=a\) or \(c=b\),
respectively. Let \(w_C\) denote the maximum weight of a term at
that endpoint in the full support, and put \(r=C-c\). Replace the
endpoint hypothesis by the following check at each endpoint with
\(r>0\):
\[
             n-(r-C)_+(n-w_C)>r.
\]
With the same hypothesis at finite nonzero roots, the conclusion of
the weighted theorem still holds. In particular at an axial endpoint
\(c=0\) the condition is simply \(r<n\).

Use the same \(G(v,w)\) at either endpoint. Since the actual endpoint
occurs, the constant coefficient \(G(v,0)\) has a first nonzero term
of order \(k=d(n-w_C)\) in \(v\). When \(r>0\), this order is positive.
Uniformly on a closed smaller bidisc, \(|G(v,w)-G(v,0)|\le C_1|w|\).
Hence for \(|w|<c_0|v|^k\), with a fixed sufficiently small \(c_0\),
one has \(|G(v,w)|\gg|v|^k\). This cannot satisfy
\(|G(v,w)|\le M|v|^{dn}|w|^C\) for small \(v\), since \(dn+kC>k\).
Thus the dyadic sum in the proof stops at radii comparable to
\(|v|^k\). For \(C=r\) it adds at most a logarithmic factor to
\(|v|^{2dn/r}\). For \(C<r\) it is dominated by its smallest
annulus, giving the power
\(|v|^{2(dn-k(r-C))/r}\). For \(C>r\) the previous geometric sum
already applies. In all cases the displayed strict inequality
makes the remaining two-dimensional radial integral converge.

For instance, take \(\kappa=1\) and \(n\ge6\). A phase with initial
monomial \(x^2t^{n-2}\), full angular endpoints \(-1,2\), and
maximum weight \(n-1\) at angle \(-1\) passes this test. At the
negative end \(C=1\), \(c=-2\), \(r=3\), and the left side of the
inequality is \(n-2>3\). There are no nonzero finite roots and
the positive endpoint has \(r=0\). The regular phase
\(x^2t^{n-2}+yt^{n-2}\) realizes this support pattern. Its leading
binary symbol in ordinary sphere degree has a root of multiplicity
\(n+2\).

## Automatic multiplicity bounds for crossing faces

For \(0<\kappa<2\), write \(\kappa=k/d\) in lowest terms. If
\(a,b\ge0\), nonnegative heights imply \(\kappa a\le n\). The
coefficient at negative angle \(-b\) is divisible by
\((1-t^2)^b\) in the literal chart. Its highest height is therefore
at least \(2b\), giving \((2-\kappa)b\le n\). All angular exponents
on the top line are in one residue class modulo \(d\), so every
nonzero root of \(Q\) has multiplicity at most
\[
       \frac{a+b}{d}\le n\left(\frac1k+\frac1{2d-k}\right).
\]
In particular, if \(k\ge2\) and \(2d-k\ge2\), this bound is strictly
smaller than \(n\). Equality in the parenthetical bound could only
occur for \(k=2d-k=2\), which is excluded by lowest terms. At all
these slopes the angular root check is automatic and only the
endpoint checks are needed.

The strict finite-root threshold is sharp for the volume argument.
For \(n\ge2\), the phase
\(f=x^{-1}(t-1)(x-t-1)^n\) is regular on the sphere: its coefficient
of \(x^{-1}\) is divisible by \(1-t^2\). At \(\kappa=1\) its top
angular polynomial is \(z^{-1}(z-1)^n\); both endpoint losses are
zero. On \(x=t+1+w\) with \(w\) bounded and \(|t|\) large, \(f\) is
bounded while \(dx/x\wedge dt=dw\wedge dt/(t+1+w)\) has divergent
squared integral. This example concerns the boundary-volume
method at multiplicity \(n\).

## The full even-conductor strip

In the [full even-conductor strip](CM_NEWTON_FULL_WIDTH_STRIP_MOMENTS.md),
take \(\kappa=2\) and \(Q(z)=z^{-1}P(z)\), with
\(1\le a\le n/2\), \(\deg P=a+1<n\) and \(P(0)\ne0\).
Regularity and the strip give \(B=b=1\). The upper weight bound
gives \(A\le n/2\), so \(A-a<n\). The geometric root hypothesis
follows from \(\deg P<n\). Thus every compatible fixed phase has
infinitely many nonzero unmarked moments. This includes upper-half
source cases with arbitrary odd lower terms and every permitted
geometric configuration of the roots. The same application permits
any angular support in \([-1,\lfloor n/2\rfloor]\) under
\(j+2i\le n\), independently of the lower strip bound, provided
the same top face is fixed.

## Exact primes for cyclically symmetric angular polynomials

An arithmetic witness is available within this geometric conclusion
for a symmetric angular polynomial. Assume the full even-conductor
strip, with \(n\) even and \(P(z)=U(z^R-c)^e\), where
\(U,c\ne0\), \(R\ge2\), and \(a=Re-1\le n/2\). Choose
\(1\le s<n\) with \(\gcd(s,Rn)=1\) and
\[
         1+\left\lfloor\frac{as}{n}\right\rfloor
                       -R\left\lfloor\frac{es}{n}\right\rfloor=R.
\]
For every sufficiently large rational prime
\(p\equiv s^{-1}\pmod{Rn}\), the single moment of order
\(m=(sp-1)/n\) is nonzero. The certifying progression has
natural density \(1/\varphi(Rn)\), with no additional hypothesis
on lower terms.

Indeed \(R\mid m\). Every literal term has weight deficit 0, 1
or 2 by the [full-width support calculation](CM_NEWTON_FULL_WIDTH_STRIP_MOMENTS.md#the-exact-support-bound).
In \(f^m\) the deficit is at most \(2m<2p\). At angular coefficient
zero only deficits 0 and \(p\) give an interval denominator
divisible by \(p\), and the height for deficit \(p\) is odd.
At any good coefficient reduction we obtain
\[
 sp\Lambda(f^m)\equiv[z^m]P(z)^m
        =U^m(-c)^{em-m/R}\binom{em}{m/R}\pmod p.
\]
Writing \(x=s/n\) and \(y=\{ex\}\), the floor condition becomes
\(\lfloor Ry-x\rfloor=R-1\), so \(y>x/R\). For sufficiently large
\(p\), the least base-\(p\) digit of \(em\) is greater than
\(m/R\), which is less than \(p\). The binomial coefficient is
nonzero by Lucas's formula, or directly by factoring powers of
\(1+z\) modulo \(p\). A fixed number-field specialization preserving
\(U,c\ne0\) and the prescribed support gives the same conclusion over
an arbitrary characteristic-zero field.

For the [two-residue cancellation example](CM_NEWTON_FULL_WIDTH_STRIP_MOMENTS.md#exact-checks-and-the-upper-half-obstruction)
with \(n=18\), \(s=13\) and \(P=(z^2-1)^4\), the refinement is
\(p\equiv25\pmod{36}\). It certifies the first moment for every
sufficiently large prime in that refinement and every phase in the
full strip with this top angular polynomial.

For the same fixed phase, at \(p=61\) the order is \(44\) and
\(p\Lambda(f^{44})\equiv5\pmod{61}\). At \(p=43\) the successive
orders \(31,74,117,160\) give normalized residues \((0,0,0,16)\).
In the two added orders every possible angular-constant denominator
lies strictly between \(p^2\) and \(2p^2\), so none has valuation two.
[Exact Laurent replay](verify_cm_newton_full_width_strip_moments.py).

# The quartic–quintic conductor-twelve exclusion

**Theorem QuarticQuinticConductorTwelveExclusion.** Let \(F=(P,Q)\) be an
actual complex polynomial plane Keller map in a source frame \((x,y)\)
where \(P\) is monic of positive degree in \(x\). If its exact coordinate
conductor has global \(x\)-degree twelve, its generic one-place
characteristic sequence cannot be \((5,4)\).

Write \(K=\mathbb C(y)\) and \(J=\operatorname{Jac}_{x,y}(P,Q)\ne0\).
The source coordinate order \(K[P,Q]\subset K[x]\) has normalization
\(K[x]\) and pole semigroup \(\langle4,5\rangle\) in the stated
case. The conductor and its specialization are those of the
[exact coordinate-order theorem](GENERIC_COORDINATE_ORDER_CONDUCTOR_AND_CLOSED_COLLISIONS.md).
In particular each finite source line has a nonconstant specialized
conductor in the noninvertible case under consideration.

## Actual generators and the parameter equation

Reduce the actual generators by the degree-reducing target shears
used in the [cubic–quintic case](CONDUCTOR_EIGHT_CUBIC_QUINTIC_EXCLUSION.md).
In characteristic sequence \((5,4)\) the reduced generator degrees
are four and five. Use leading coefficient scalings and a depressed
source coordinate \(u=x-h(y)\), \(h\in K\). Subtract a multiple of the
quartic from the quintic to remove its degree-four coefficient and
translate both target coordinates to remove their constants. These
ambient target changes have determinant \(\lambda\in K^*\).

Adjoin \(\rho^9=\lambda\), set \(u=\rho z\), and scale the two target
coordinates by \(\rho^{-4}\) and \(\rho^{-5}\), respectively. Over
\(L=K(\rho)\) the target determinant is one, and the generator pair is
\[
 U=z^4+az^2+bz,\qquad V=z^5+dz^3+ez^2+fz .
 \tag{1}
\]
The parameter derivative of the target change, composed with its
inverse, is a polynomial vector field of divergence zero. Write it
as \((-H_V,H_U)\) with \(H\) polynomial. The source change gives
\[
 U_zV_y-U_yV_z=J\rho+\partial_z H(U,V).
 \tag{2}
\]
The left side has \(z\)-degree at most six. As the only positive
degrees at most seven in the exact pole semigroup are four and five,
\(H(U,V)=h_U U+h_V V+\text{constant}\).

At degrees six and five comparison gives
\(D=d-5a/4\) and \(E=e-5b/4\) in \(\mathbb C\). The next two
degrees give
\[
 h_U=-(5a/16+3D/4)b'-(5b/16+E/2)a',
 \qquad h_V=4f'/5-(a/4+3D/5)a' .
\]
With primes denoting \(y\)-derivatives, define
\[
\begin{aligned}
 A&=(a+12D/5)f+5b^2/8+2Eb-5a^3/16-3Da^2/2-9D^2a/5,\\
 B&=(b+8E/5)f-(5a^2/8+3Da/2)b-3Ea^2/4-6DEa/5,\\
 G&=-2f^2/5+(5a/32+3D/8)b^2+Aa/4
       +5a^4/256+Da^3/8+9D^2a^2/40.
\end{aligned}
\]
The residual coefficients in degrees two, one and zero are
\(-A'\), \(-B'\) and \(G'-aA'/4\), respectively. Consequently
\[
 A,B\in\mathbb C,\qquad G'=J\rho.
 \tag{3}
\]
The [exact coefficient replay](verify_conductor_twelve_quartic_quintic_exclusion.py)
checks these identities before imposing \(A'=B'=0\).
Under an automorphism of \(L/K\) sending \(\rho\) to \(\zeta\rho\),
quantities of weight \(j\) transform by \(\zeta^{-j}\). The weights
of \(a,b,f,D,E,A,B,G\) are \(2,3,4,2,3,6,7,8\).

## Matching pairs and poles of the normalized coefficients

At every finite \(y_0\), at least one generic matching pair has
finite distinct limits in the original \(x\)-coordinate. A root
branch of the exact conductor reaches a finite root of the
nonconstant specialized conductor. Monicity of \(P\) makes the
matching coordinate finite. If the limits coincide, the divided
differences of both \(P\) and \(Q\) force \(P_x=Q_x=0\) there,
contradicting \(J\ne0\). This argument takes place in any extension
of the local valuation field containing all generic pairs. Thus
a physical gap \(\rho(z_+-z_-)\) must be a unit at each finite place.

For matching pairs of monic generators of degrees four and five
with regular coefficients, both coordinates are bounded. Otherwise
scale the pair by their largest pole and take the limiting divided
differences of \(z^4\) and \(z^5\). These forms have no simultaneous
nonzero zero. Away from the diagonal this follows from coprimality
of four and five, and on the diagonal it is the derivative test
for these powers. If the limiting generator pair is immersed,
the normalized gap is also a unit, since collapse would force
both derivatives of the limiting generators to vanish. We use
these observations after finite extensions to include all
valuation branches.

Suppose some \(a,b,f\) has a pole at a finite place of \(L\).
Let \(v\) be its integer valuation and \(\ell=v(y-y_0)\). Put
\[
 k=\max\{-v(a)/2,-v(b)/3,-v(f)/4\}>0.
\]
A local extension supplies a scale of valuation \(-k\).
Scale \(z\) by it and the outputs by its inverse fourth and fifth
powers. Write \(\alpha,\beta,\phi\) for the limits of the scaled
\(a,b,f\). From the top weights of \(A,B\) one obtains
\[
 \alpha\phi+5\beta^2/8-5\alpha^3/16=0,\qquad
 \beta(\phi-5\alpha^2/8)=0 .
\]
At least one of \(\alpha,\beta,\phi\) is nonzero. The three cases,
with the corresponding nonzero scaled limit of \(G\), are
\[
\begin{array}{c|c}
(\alpha,\beta,\phi)&\text{scaled limit of }G\\ \hline
 \alpha\ne0,\quad\beta^2=-\alpha^3/2,\quad\phi=5\alpha^2/8
              &-55\alpha^4/256\\
 \alpha\ne0,\quad\beta=0,\quad\phi=5\alpha^2/16
              &-5\alpha^4/256\\
 \alpha=\beta=0,\quad\phi\ne0
              &-2\phi^2/5 .
\end{array}
 \tag{4}
\]
All three limiting pairs are immersed. In the first two cases
normalize \(\alpha=1\) by a constant source scaling and monic
target scalings. Simultaneous derivative vanishing is equivalent to
\[
 4z^3+2z+\beta=0,\qquad z^2+\beta z+4\phi/5=0 .
\]
Their respective resultants are
\(-\beta^2(2\beta-3)(2\beta+3)\) and \(1/4\). The first is nonzero
when \(\beta^2=-1/2\). In the third case \(U_z=4z^3\) and
\(V_z(0)=\phi\ne0\).

The boundedness and immersion tests now give \(v(z_+-z_-)=-k\)
for every matching pair before this local scaling. Equation (4)
gives \(v(G)=-8k\) and (3) gives \(v(\rho)=-8k-\ell\), since
differentiation of a pole with respect to \(y\) reduces its
valuation by \(\ell\). All physical gaps have negative
valuation. This contradicts the finite-pair condition.
It follows that \(a,b,f\) are regular at every finite place.

## Ramified finite places

The regularity just proved makes all normalized matching
coordinates bounded at finite places. A zero of \(\rho\) would
make every physical gap vanish, so \(\rho\) has no finite zero.
At an unramified finite place the regularity of \(G\) implies
the regularity of \(G'\), hence of \(\rho\).

It remains to consider ramification. The Kummer extension
\(L/K\) has degree dividing nine. At a place of ramification
index nine, inertia and the stated weights force
\(D=E=A=B=0\). We treat these homogeneous constants below.
At ramification index three, inertia gives \(D=B=0\) and
the residues of \(a,f\) are zero. Write \(b_*\) for the
residue of \(b\). If \(\rho\) has a pole, a finite unit physical
gap requires a normalized pair to collapse.

If \(b_*=0\), evaluation gives \(A=0\). All \(a,b,f\) vanish
at the place, and \(v(G)\ge2\). Since \(G'=J\rho\) and
\(\ell=3\), a pole forces \(v(G)=2\) and \(v(f)=1\).
If \(E\ne0\), the equation \(B=0\) gives
\[
 f=a^2(5b/8+3E/4)/(b+8E/5).
\]
Its valuation is at least two, a contradiction. If \(E=0\)
the constants are homogeneous.

Suppose instead \(b_*\ne0\). Every common critical point of
the limiting pair \(U=z^4+b_*z,\ V=z^5+(5b_*/4+E)z^2\)
satisfies \(z^3=-b_*/4\) and forces \(E=-5b_*/8\).
Thus \(A=-5b_*^2/8\). Put \(\delta=b-b_*\).
The term linear in \(a\) cancels from \(G\), whose expression is
\[
 G=-2f^2/5+(5a/32)(b^2-b_*^2)+5a^4/256.
\]
Inertia fixes \(b\), so either \(\delta=0\) or \(v(\delta)\)
is a positive multiple of three. Again a pole forces
\(v(G)=2,\ v(f)=1,\ v(\rho)=-1\). The two invariant equations
become
\[
 af+5\delta^2/8-5a^3/16=0,\qquad
 \delta f-a^2(5\delta/8+5b_*/32)=0 .                    \tag{5}
\]
If neither \(a\) nor \(\delta\) vanishes identically, let
\(r=v(a)\ge1\), \(m=v(\delta)\ge3\). Since \(r+1<3r\),
comparison in (5) gives \(r+1=2m,\ m+1=2r\), forcing
\(r=m=1\). This is impossible. If either is zero, (5)
and \(f\ne0\) make both zero.

In that last case \(U=z^4+b_*z,\ V=z^5+(5b_*/8)z^2+fz\).
For a collapsing pair near a common critical point \(z_0\)
write its endpoints \(z_0+h\pm v_0\), where \(z_0^3=-b_*/4\).
The divided difference for \(U\) gives
\(h=-v_0^2/(3z_0)+O(v_0^4)\). Substitution in the divided
difference for \(V\) gives
\(f+5z_0^2v_0^2+O(v_0^4)=0\), with nonzero quadratic
coefficient. These expansions follow by solving the first
polynomial for \(h\), whose derivative at zero is \(12z_0^2\).
Hence the normalized gap has valuation \(1/2\) and its
physical gap has valuation \(-1/2\). A pair with distinct
limiting endpoints also has a physical pole. No finite
unit gap remains.

If the constants are nonhomogeneous, the above excludes all
finite poles and zeros of \(\rho\). Since \(\rho^9\in K\),
\(\rho\) is constant and \(a,b,f\) are polynomials in \(y\).
At infinity any pole among them obeys the same leading
equations (4), so \(\deg_y G=8k\) for
\(k=\max(\deg a/2,\deg b/3,\deg f/4)>0\). A nonconstant
coefficient gives \(k\ge1/4\), hence \(\deg G\ge2\).
If all coefficients are constant then so is \(G\).
Both possibilities contradict \(G'=J\rho\in\mathbb C^*\).

## Homogeneous constants and the specialized generators

Assume \(D=E=A=B=0\). The top equations in (4) now hold
exactly at the generic point, with \(a,b,f\) in place of
\(\alpha,\beta,\phi\). They are not all zero by (3).
In the three cases take respectively a scale \(\tau=b/a\),
\(\tau^2=a\), or \(\tau^4=f\), allowing a finite extension
and nonzero constant scalings. In the variable \(z/\tau\)
the scaled generators are fixed complex polynomials of
degrees four and five with the immersion property of (4).
They have only finitely many matching pairs. Indeed their
function field has index dividing both four and five in
the rational function field of their source, so the map is
birational onto its image. Each pair has a fixed nonzero gap.

The gaps in the original source are then fixed nonzero
multiples of \(w=\rho\tau\). The finite-pair condition makes
\(w\) a unit at every finite place. In the three cases
\(w,w^2,w^4\), respectively, belong to \(K\) by the ratios
or powers defining \(\tau\) in the unscaled coefficients.
They have no finite zeros or poles. Thus \(w\in\mathbb C^*\).
By (4) \(G=g_0\tau^8\) for \(g_0\in\mathbb C^*\).
Equation (3) becomes \(8g_0\tau^8\tau'=Jw\), so
\(\tau^9\) is affine and nonconstant in \(y\).
It follows that \(\lambda^{-1}=\tau^9/w^9\) has a finite zero.

In the centred source coordinate the unscaled generators
are now fixed. Also \(h(y)\) has no finite pole since the
physical matching coordinates have the form \(h(y)\) plus
fixed constants, and one pair stays finite at every finite
parameter. Hence \(h\in\mathbb C[y]\). Let \(R\subset
\mathbb C[u]\) be the fixed order of colength six.
The original \(P,Q\) after centring belong coefficientwise
to \(R[y]\), since generic membership in the fixed vector
subspace \(R\) imposes constant-coefficient linear equations
on each bounded-degree part. The polynomial exact conductor
in the centred coordinate is a scalar polynomial in \(y\)
times the fixed monic conductor of \(R\). This scalar
cannot vanish at any finite \(y\), by exact specialization
of the conductor. At the zero of \(\lambda^{-1}\) the
specialized order is therefore included in \(R\) and has the
same colength, so it equals \(R\).

For every fixed matching pair the original action
\(\int P\,\partial_u Q\,du\) is \(\lambda^{-1}\) times
the fixed canonical action. A polynomial ambient target
change scales the action by its determinant because the
one-form error is exact and the endpoints have equal target
values. The actions continue polynomially at the fixed
endpoints and all vanish at the zero of \(\lambda^{-1}\).
The fixed canonical action vector is nonzero because the
derivative of each physical action is \(-J\) times its nonzero
gap, by differentiation with the endpoint terms cancelling.

We check the ambient degree sequences on the specialized
order \(R\). The
[one-place sequence theorem of Assi and García-Sánchez,
Propositions 2 and 4](sources/assi-garcia-sanchez-1407.0490v1.pdf)
applies as in the cubic–quintic proof. For semigroup
\(\langle4,5\rangle\) the only reduced sequences of
conductor twelve are \((5,4)\) and \((10,4,5)\).
Here is the numerical check. For one gcd drop, the only
two minimal generators are four and five. For two drops,
write \(d=\gcd(r_0,r_1)>1\) and
\(\mu=(r_0/d-1)(r_1/d-1)\ge2\). The conductor formula gives
\(12=d\mu+(d-1)(r_2-1)\). For \(d=2\) this permits
\((6,4,9)\) and \((10,4,5)\), of which the first uses the
nongap \(6\) outside the specified semigroup. For \(d=3\)
it permits \((9,6,4)\), again containing \(6\), and \(d\ge4\)
gives too large a bound because \(r_2\ge4\).
Three or more drops give an initial gcd \(d_2\ge4\).
If \(d_2\ge6\), the lower bound \(2d_2+(d_h-1)(r_h-1)\)
already exceeds twelve. If \(d_2=4\), the next gcd is two,
so \(r_2\ge6\) and the last \(r_3\ge5\). The conductor is
at least \(8+(r_2-2)+(r_3-1)>12\) by the prefix formula.

The sequence \((10,4,5)\) would give a specialized degree-four
generator affine in the fixed quartic \(U\) and, after
rescaling, a degree-ten generator
\[
 W=V^2+g(U)V+h_0(U)
\]
in the fixed generators \(U,V\). Indeed the filtered basis
through degree ten is \(1,U,V,U^2,UV,V^2\).
The one-place approximate root for this sequence is monic
of degree five in the degree-four generator and has
longitudinal degree five. Weights \(4i+10j\) for \(0\le i<5\)
are distinct modulo ten as \(i\) varies. Thus the only
terms of weight twenty are \(U^5\) and a multiple of \(W^2\);
terms involving \(W\) at lower weights are affine in \(W\)
with polynomial coefficients in \(U\). Completing the square
by a lower triangular target shear gives \(W^2-R_0(U)\) of
longitudinal degree five. Consequently this is
\(\gamma V+\mu_0 U+\nu_0\) with \(\gamma\ne0\), since the basis
through degree five is \(1,U,V\).

The monic minimal polynomial of \(V\) over \(\mathbb C(U)\)
has degree four. The putative relation writes it as
\[
 (V^2+g(U)V+h_0(U))^2-R_0(U)-\gamma V-\mu_0U-\nu_0 .
\]
After depressing the quartic in \(V\), its linear coefficient
would be the nonzero constant \(-\gamma\). For the three
fixed shapes of (4), normalize to
\((a,b,f)=(1,\beta,5/8)\), \(\beta^2=-1/2\);
\((1,0,5/16)\); and \((0,0,1)\). Computing
\(\operatorname{Res}_z(U(z)-X,V(z)-Y)\) and depressing it in
\(Y\) gives the linear coefficients
\[
       \beta(80X+41)/1024,\qquad 0,\qquad 0,
\]
respectively. A constant rescaling of source and targets
preserves the distinction between these forms and a nonzero
constant. Therefore the specialized sequence cannot be
\((10,4,5)\).

For the only remaining sequence \((5,4)\), the specialized
reduced generators are related to \(U,V\) by an ambient
invertible affine change, using the basis through degree
five. Their action vector is therefore a nonzero determinant
multiple of the fixed nonzero vector. The ambient changes
from the original specialized generators also have nonzero
determinants. This contradicts vanishing of the original
actions at the chosen specialization and proves the theorem.

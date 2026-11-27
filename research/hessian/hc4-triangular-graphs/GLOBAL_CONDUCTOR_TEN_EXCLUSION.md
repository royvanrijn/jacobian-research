# Exclusion of global directional conductor degree ten

**Theorem KellerDirectionalConductorTenExclusion.** Let \(F=(P,Q)\)
be a noninvertible complex polynomial plane Keller map. In a source
frame \((x,y)\) with \(P\) monic of positive degree in \(x\), the
global \(x\)-degree of its exact directional coordinate conductor is
at least twelve. Individual specialized conductor degrees can
still be smaller.
The [degree-twelve exclusion](GLOBAL_CONDUCTOR_TWELVE_EXCLUSION.md)
raises this global bound to fourteen.

We exclude degree ten below. The evenness of the exact conductor
degree and the previous exclusions through eight are in the
[coordinate-order theorem](GENERIC_COORDINATE_ORDER_CONDUCTOR_AND_CLOSED_COLLISIONS.md)
and the [global octic exclusion](GLOBAL_OCTIC_CONDUCTOR_EXCLUSION.md).
Put \(K=\mathbb C(y)\). At each finite parameter, a generic collision
pair can be selected with distinct finite limiting source
coordinates. The specialized conductor is nonconstant, so one
generic conductor root tends to a finite root. Its mate stays
finite by monicity of \(P\). If their limits coincided, the
divided differences of the original \(P,Q\) would force both
longitudinal derivatives to vanish. This argument works after
further local field extensions, as in the
[cubic–quintic finite-pair argument](CONDUCTOR_EIGHT_CUBIC_QUINTIC_EXCLUSION.md#3-finite-pairs-and-coefficient-poles).

## 1. The characteristic sequences in degree ten

Reduce the generator degrees over \(K\) by triangular target
automorphisms. Apply
[Assi–García-Sánchez, Propositions 2 and 4 and §6](sources/assi-garcia-sanchez-1407.0490v1.pdf)
to the resulting one-place equation as in the
[cubic–quintic classification](CONDUCTOR_EIGHT_CUBIC_QUINTIC_EXCLUSION.md#1-the-characteristic-sequences-of-the-exact-order).
It gives a free symmetric semigroup with reduced characteristic
degrees \(r_0>r_1>\gcd(r_0,r_1),r_2,\ldots,r_h\). The gap count
equals the order's colength by the filtration by longitudinal
degree. The actual conductor degree ten gives colength five,
so the numerical semigroup conductor is ten.

For \(d_i=\gcd(r_0,\ldots,r_{i-1})\) and \(e_i=d_i/d_{i+1}\),
the prefix conductors \(\mu_i\) of the normalized sequences ending
before \(r_i\) are even, and \(d_i\mu_i\) is nondecreasing, with
\(d_2\mu_2\ge2d_2\). The full formula for \(h\ge2\) is
\[
10=d_h\mu_h+(d_h-1)(r_h-1).
\]
For \(h=1\) the equation \((r_0-1)(r_1-1)=10\) with coprime
\(r_0,r_1\) gives \((11,2)\). If \(h\ge3\), strict divisibility of
the \(d_i\) and this formula require \(d_2=4,d_3=2,\mu_2=2\),
hence \(r_0,r_1=12,8\). Freeness gives
\(2r_2\in\langle12,8\rangle\), and \(\gcd(4,r_2)=2\), so \(r_2\ge6\).
But \(d_3\mu_3=8+r_2-2\ge12\), a contradiction. For \(h=2\),
\(d_2\ge4\) already gives a lower bound \(2d_2+(d_2-1)>10\).
The case \(d_2=3\) forces \(r_2=3\), violating the gcd condition.
The two candidates for \(d_2=2\) are \((6,4,7)\) and \((10,4,3)\).
The latter violates the necessary \(2r_2\in\langle10,4\rangle\),
since \(2r_2=6\). Thus only \((11,2)\) and \((6,4,7)\) remain.

## 2. The hyperelliptic branch \((11,2)\)

The [all-degree exclusion of degree-two elements](DEGREE_TWO_ELEMENTS_AND_KELLER_REFLECTION.md)
rules out this sequence. The following period argument gives another
proof of the branch.

Write \(u=x-h(y)\) for the centred longitudinal coordinate. The
canonical generators are \(T=u^2\) and \(V=u g(T)\), where
\(g(T)=\sum_{j=0}^5 a_j T^j\). Choose \(g(u^2)\) to be the
centred actual conductor, rescaling the degree-eleven generator
over \(K\) accordingly. As in
[the octic reflected-pair proof](HYPERELLIPTIC_OCTIC_CONDUCTORS_AND_FOUR_PERIOD_RIGIDITY.md#1-the-centre-and-the-four-reflected-pairs),
monicity of \(P\) and the existence of a finite conductor root
make \(h\) polynomial: a pole of \(h\) would send the reflected
mate of that root to infinity. Thus all \(a_j\) are polynomials.
A finite zero of \(a_0\) would collapse a reflected pair and
force \(P_u=Q_u=0\), so \(a_0\) is a nonzero constant. A persistent
multiple root of \(g\) would make the \(u\)-derivatives of \(P,Q\)
opposite at its reflected pair and make the \(y\)-derivatives
equal there by differentiating the collision. This would make
the constant Jacobian change sign. Hence the five roots of \(g\)
are distinct.

Use the physical actions \(\int P\,dQ\) between matched source
points on a line of fixed \(y\). They scale the canonical actions
by one nonzero \(\Lambda\in K\). For \(u=\pm r\) the canonical
action is \(-4r\sum_{j=0}^5 a_j(r^2)^{j+1}/(2j+3)\).
Differentiation as the integral of \(-2\Lambda u^2g(u^2)\,du\)
removes the moving endpoint terms. The physical action
derivative is \(-2Jr\), with \(J\ne0\) the original Jacobian.
The resulting polynomial identity at the five roots of \(g\)
gives, by division by \(g\),
\[
\frac{(\Lambda a_j)'}{2j+3}-\alpha a_j=\beta a_{j+1}
\quad(0\le j<5),\qquad
\alpha=\frac{(\Lambda a_5)'}{13a_5},\qquad
\beta=-\frac{J}{2a_0}.
\tag{1}
\]

Suppose \(a_5\) has a finite zero of order \(m>0\). Since the
specialized conductor is nonconstant, let \(k\in\{1,2,3,4\}\)
be the largest index with \(a_k\) a unit there. Write \(n\) for
the order of \(\Lambda\). If \(n<1\), avoiding negative powers in
(1) at \(j=0,k\) requires
\(n/3=(n+m)/13=n/(2k+3)\), which is impossible. If \(n>1\),
the left side at \(j=k-1\) has positive order or is zero,
while the right side is a unit. For \(n=1\) the equation
at \(j=k\) gives \(m=(10-2k)/(2k+3)\), one of
\(8/5,6/7,4/9,2/11\), never a positive integer.

Now \(a_5\) is constant. Equation (1) at \(j=0\) shows
\(\Lambda'\) is polynomial, so \(\Lambda\) is polynomial.
It is nonconstant by (1) and \(a_5\ne0\). If its degree is
\(L\ge1\), the leading factor at step \(j\) of (1), assuming
\(\deg a_j=j(L-1)\), is
\[
\frac{j(L-1)}{2j+3}
  +L\!\left(\frac1{2j+3}-\frac1{13}\right)>0 .
\]
Thus \(\deg a_j=j(L-1)\) for all \(j\). It follows that \(L=1\)
and the entire \(g\) is constant.

At the zero of \(\Lambda\), the five fixed pair equalities still
hold. They are independent by interpolation and cut out the
fixed order \(\mathbb C[T,V]\) of colength five. The exact
specialized conductor bounds the colength of the specialized
order by five, so that order equals \(\mathbb C[T,V]\). Its
generators can be reduced afresh to the degrees \(11,2\):
the other sequence \((6,4,7)\) has no element of degree two.
For generators in degrees \(11,2\), triangular shears send them
to the canonical pair with nonzero determinant. The canonical
action vector is nonzero. Otherwise the polynomial
\(\sum_j a_j T^j/(2j+3)\) would vanish at all five roots of \(g\)
and be proportional to \(g\), contrary to its constant and
highest coefficients. Thus the specialized physical action
vector is nonzero. The identity with factor \(\Lambda\) and
the fixed endpoints instead gives zero there by polynomial
continuation. This excludes \((11,2)\).

## 3. Parameter identities for the branch \((6,4,7)\)

In this branch the approximate root has degree seven.
The same weight comparison and square completion as in the
[degree-eight cable proof](GLOBAL_OCTIC_CONDUCTOR_EXCLUSION.md)
give \(W=V^2-U^3\) of degree seven after triangular target
changes. Indeed the approximate root is monic of degree
three in the degree-four generator, with its three coefficients
in the degree-six generator. Their weights modulo six are
distinct, permitting only terms of source degree at most twelve.
Completing the square in the degree-six generator and translating
the degree-four generator removes the quadratic term in that
generator. The remaining affine expression in the latter
has degree at most four. Therefore the seventh degree in the
approximate root is exactly the degree of \(V^2-U^3\).

Make both generators monic, depress the quartic by a source
translation, and subtract a target constant from the sextic.
Then, in the longitudinal parameter \(z\), they have the form
\[
U=z^4+az^2+bz+c,\qquad D=z^2+a/2,\qquad
V=3UD/2-D^3/2+3b^2/8+\kappa z,\qquad \kappa\ne0 .
\tag{2}
\]
Indeed \(V-U^{3/2}\) has top degree one in its Laurent expansion,
and the polynomial part of \(U^{3/2}\) is
\(3UD/2-D^3/2+3b^2/8\).
Let \(\lambda\in K^\times\) denote the target determinant
before scaling. Adjoin \(\rho^{10}=\lambda\), replace the
centred longitudinal variable by \(\rho z\), and scale the
target outputs by \(\rho^{-4},\rho^{-6}\). In (2) use the
new coefficients in \(L=K(\rho)\). The new target determinant
is one and the source Jacobian in \((z,y)\) is \(J\rho\).

As in the [parameter-compatibility proof](CONDUCTOR_EIGHT_CUBIC_QUINTIC_EXCLUSION.md#2-polynomial-parameter-compatibility-for-the-cubicquintic-pair),
the derivative of the actual target change is a polynomial
Hamiltonian vector field, since its divergence is zero.
Hence
\[
U_z V_y-U_y V_z=J\rho+\partial_z H(U,V).
\tag{3}
\]
The left side has degree at most four in \(z\), so
\(\deg H(U,V)\le5\). The semigroup \(\langle4,6,7\rangle\)
permits only a linear combination of \(1,U\) in this degree.
The degree-three coefficient on the left vanishes, so the
\(U\) term also vanishes.

Put \(\Delta=a^2-4c\). The degree-four term in (3) is \(4\kappa'\),
so \(\kappa\in\mathbb C^\times\). Direct differentiation of (2)
now gives the following constants and primitive:
\[
A=b\Delta-\tfrac43\kappa a,\qquad
B=\Delta^2-8ab^2+\tfrac{32}3\kappa b,\qquad
G=\tfrac14 b^3-\tfrac3{16}Aa-\tfrac38\kappa a^2+\tfrac14\kappa\Delta,
\qquad G'=J\rho .
\tag{4}
\]
The degree-two and degree-one terms are \(3A'/4\) and
\(-3B'/32\). Once these vanish, the constant term equals
\(G'\). The constants of the finite algebraic extension with
this derivation are \(\mathbb C\). The correction coefficient is
\(\kappa=\kappa_0/\rho^5\) with \(\kappa_0\in K\).
Therefore \(\rho^5\in K\) and \(L/K\) is trivial or cyclic of
degree five. If \(b\) is identically zero, (4) makes \(a,\Delta,G\)
constant, contradicting \(G'=J\rho\). In general, elimination
of \(a\) in (4) gives
\[
a=\frac{3(b\Delta-A)}{4\kappa},\qquad
\Delta^2-\frac{6}{\kappa}b^3\Delta
   +\frac{6}{\kappa}A b^2+\frac{32}3\kappa b-B=0 .
\tag{5}
\]

## 4. Parameter poles and the physical gaps

For a pair of distinct points with the same \(U,V\), let \(S\)
be their sum and \(g_*\) their difference in \(z\). The divided
\(U\) difference gives their product \((S^3+aS+b)/(2S)\).
The equations cannot have \(S=0\): the first divided
difference would give \(b=0\) and the second would give
\(\kappa=0\). Using (2) for the divided \(V\) difference gives
\[
2S^5+4a S^3-2b S^2+3\Delta S-8\kappa=0,\qquad
g_*^2=-S^2-2a-\frac{2b}{S}
      =-\frac{3b}{S}+\frac{3\Delta}{2S^2}-\frac{4\kappa}{S^3}.
\tag{6}
\]
The actual longitudinal gap is \(\rho g_*\).

Consider a finite parameter place of \(L\) where \(b\) has a pole
of order \(q>0\). Equation (5) allows precisely the two leading
forms in the following table, with \(\sim\) denoting equality
of leading terms. To obtain the root scales in (6), substitute
the indicated scale and retain the highest powers of \(b\).
The resulting face polynomials have simple nonzero roots on
the stated scales, accounting for all five roots after a
local field extension.

| \(\Delta\) branch | Count | Leading scale of \(S\) | Leading squared gap |
|---|---:|---|---|
| \(\Delta= A/b+O(b^{-2}),\ a=O(b^{-1})\) | 3 | \(S^3\sim b\) | \(g_*^2\sim-3S^2\) |
| same | 2 | \(S^2\sim-4\kappa/b\) | \(g_*^2\sim-2b/S\) |
| \(\Delta\sim6b^3/\kappa,\ a\sim9b^4/(2\kappa^2)\) | 2 | \(S/b^2\to \xi=\pm3i/\kappa\) | \(g_*^2\sim(-3/\xi+9/(\kappa\xi^2))b^{-1}\) |
| same | 2 | \(S^2\sim-\kappa/b\) | \(g_*^2\sim-9b^4/\kappa^2\) |
| same | 1 | \(S\sim4\kappa^2/(9b^3)\) | \(g_*^2\sim-27b^4/(2\kappa^2)\) |

For the two roots of order \(b^2\), the second expression for
\(g_*^2\) gives the displayed nonzero coefficient because
\(\xi\ne3/\kappa\). In the small-\(\Delta\) branch, (4) gives
\(G\sim b^3/4\); in the large-\(\Delta\) branch it gives
\(G\sim-243b^8/(32\kappa^3)\). At a finite parameter place of
ramification index \(l\), differentiating with respect to \(y\)
increases the pole order of \(G\) by \(l\). Thus \(\rho=G'/J\)
has pole order \(3q+l\) or \(8q+l\). Every physical gap
then diverges, including the two with normalized squared gap
of order \(b^{-1}\). This contradicts the finite distinct pair
condition. Therefore \(b\) is regular at every finite place.

The monic quadratic (5) makes \(\Delta,a,c,G\) regular there.
The quintic (6) has constant nonzero leading and constant
coefficients, so \(S,S^{-1}\) and \(g_*\) are regular as well.
A finite zero of \(\rho\) would collapse every physical gap.
At an unramified finite place, regularity of \(G\) makes \(G'\)
regular, excluding a pole of \(\rho\). At a ramified finite
place, inertia has order five. It multiplies \(\rho\) by a
primitive fifth root of unity and acts on \(a,b,c,\Delta\) with
weights \(-2,-3,-4,-4\) respectively. The fixed residue field
is \(\mathbb C\), so the regular residues of these coefficients
all vanish. Equation (6) gives \(S^5=4\kappa\) in the residue
field and \(g_*^2=-S^2\ne0\) for every root. A pole of \(\rho\)
would again make every physical gap diverge. Consequently
\(\rho^5\in K\) has no finite zero or pole, and \(\rho\) is a
nonzero constant.

The coefficients now belong to \(K\), and \(b,\Delta,G\) are
polynomials. If \(b\) is constant, (5) makes \(\Delta\) constant,
and \(G'=0\). If \(b\) is nonconstant, apply the same two
leading forms of (5) at infinity. They make \(\deg G\) three
or eight times \(\deg b\). Both contradict \(G'=J\rho\ne0\)
constant. This excludes \((6,4,7)\) and proves the degree-twelve
lower bound.

The [exact identity check](verify_global_conductor_ten_exclusion.py)
replays the coefficient and divided-difference equations used
above.

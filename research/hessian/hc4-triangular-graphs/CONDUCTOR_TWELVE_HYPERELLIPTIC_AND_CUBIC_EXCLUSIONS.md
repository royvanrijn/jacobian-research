# The hyperelliptic and cubic branches of conductor degree twelve

Let \(F=(P,Q)\) be a complex polynomial plane Keller map in a source
frame \((x,y)\) with \(P\) monic of positive degree in \(x\). Write
\(J=\operatorname{Jac}_{x,y}(P,Q)\ne0\), \(K=\mathbb C(y)\), and
\(\kappa\) for the exact directional coordinate conductor. We assume
\(\deg_x\kappa=12\). The [one-place classification](GLOBAL_CONDUCTOR_TWELVE_EXCLUSION.md#characteristic-sequences-of-the-exact-order)
gives the possible reduced characteristic sequences. The two branches
below use the sequences \((13,2)\) and \((7,3)\). The
[semigroup input](sources/assi-garcia-sanchez-1407.0490v1.pdf) is
Assi–García-Sánchez, Propositions 2 and 4.

The [exact conductor theorem](GENERIC_COORDINATE_ORDER_CONDUCTOR_AND_CLOSED_COLLISIONS.md)
gives the actual specialized order on every line, with colength
\(\deg_x\kappa(x,y_0)/2\). At every finite parameter \(y_0\) one can
select a generic collision pair with distinct finite limiting
\(x\)-coordinates. Choose a generic conductor root tending to a
finite root of the nonconstant specialized conductor, then use
monicity of \(P\) to bound its mate, which has the same \(P\)-value.
Equal limits would make the divided differences of both \(P,Q\)
vanish on the diagonal, contrary to \(J\ne0\). This argument applies
at every finite place after a further local extension for the pair.

For a collision pair \(x_-,x_+\), use the action
\(I=\int_{x_-}^{x_+}P\,dQ\) on the source line. Its parameter derivative
is \(-J(x_+-x_-)\). Differentiating the integral gives this formula
because the endpoint values of both coordinates agree. The action in
target coordinates transformed by a polynomial automorphism is the
original action multiplied by the determinant of the automorphism.
At a fixed parameter the difference of the two one-forms, after
this scaling, is closed on the affine target plane and hence exact,
so its integral between the equal target endpoints is zero. We use
these facts over generic parameter opens and apply polynomial
continuation only when the relevant endpoints have finite limits.

## The nodal action on a specialized order

Suppose the actual specialized order has colength one and identifies
two distinct points, centered at \(u=-r,r\) with \(r\ne0\). It is
contained in the kernel of their nonzero difference functional on
\(\mathbb C[u]\). Both subspaces have colength one, so the order is
\(\mathbb C[u^2,u(u^2-r^2)]\). Its semigroup is \(\langle2,3\rangle\).
Reducing the two actual generator degrees by triangular shears
leaves degrees \(3,2\). Indeed the normalized prefix in the
one-place conductor formula has conductor at least two, so a
sequence of length at least three has conductor at least
\(2d+(d-1)(r_h-1)>2\) for \(d\ge2\) and \(r_h\ge2\).
A degree-one element is impossible in this order. For two generators
of degrees \(3,2\), rescaling and a shear give the canonical pair.
The action of that pair from \(-r\) to \(r\) is
\[
 \int_{-r}^{r}u^2\,d\bigl(u(u^2-r^2)\bigr)=8r^5/15.
\]
Every actual generating pair therefore has nonzero action.

## The degree-two branch

**Theorem ConductorTwelveHyperellipticExclusion.** The reduced
sequence \((13,2)\) is impossible under the stated hypotheses.

Use \(u=x-h(y)\) and canonical generators \(T=u^2\) and
\(V=u g(T)\), with \(g=\sum_{j=0}^6 a_jT^j\). These generators follow
from the degree-two element by taking even and odd parts of the
order. The odd factors form an ideal of \(K[T]\). Choose its
generator \(g\) so that \(g(u^2)\) is the centered exact conductor.
As in the [octic reflected-pair argument](HYPERELLIPTIC_OCTIC_CONDUCTORS_AND_FOUR_PERIOD_RIGIDITY.md#1-the-centre-and-the-four-reflected-pairs),
\(h\) is the conductor-root centroid, hence rational, and a finite
pole of \(h\) would send the reflected mate of a finite conductor
root to infinity. Thus \(h\) and the \(a_j\) are polynomials.
A zero of \(a_0\) collapses a reflected pair and contradicts the
Keller identity. Thus \(a_0\) is a nonzero constant. A persistent
multiple root of \(g\) makes the \(u\)-derivatives of \(P,Q\)
opposite at the mates and the \(y\)-derivatives equal, by
differentiating the collision identities. It would change the sign
of the Jacobian. Therefore \(g\) has six simple nonzero generic
roots.

The canonical and physical actions differ by a common
\(\Lambda\in K^\times\). At \(u=\pm r\) with \(g(r^2)=0\), the
canonical action is \(-4r\sum_{j=0}^6 a_j(r^2)^{j+1}/(2j+3)\).
Differentiating its integral form
\(-\int_{-r}^{r}2\Lambda u^2g(u^2)\,du\) cancels the endpoint
terms. The physical action derivative is \(-2Jr\). Dividing the
resulting equation in \(T\) by \(g(T)\) gives
\[
 \frac{(\Lambda a_j)'}{2j+3}-\alpha a_j=\beta a_{j+1}
 \quad(0\le j<6),\qquad
 \alpha=\frac{(\Lambda a_6)'}{15a_6},\qquad
 \beta=-\frac{J}{2a_0}.
 \tag{1}
\]

If \(a_6\) vanishes at a finite parameter, let \(m>0\) be its order
of vanishing, \(n=\operatorname{ord}\Lambda\), and \(k\) the largest
unit index. Nonconstancy of the specialized conductor gives
\(1\le k\le5\). At \(n<1\), the equations at \(j=0,k\) force
\(n/3=(n+m)/15=n/(2k+3)\), impossible since \(m>0\).
At \(n>1\), the equation at \(j=k-1\) has a left side of positive
order or zero and a unit right side. At \(n=1\), the equation at
\(j=k\) forces \(m=(12-2k)/(2k+3)\). The only positive integral
resonance is \(k=1,m=2\). In that case \(a_0+a_1 T\) is the
specialized conductor in the centered frame. Its reflected pair
is finite and distinct. Its canonical action stays finite, while
the physical action tends to zero because \(\Lambda\) vanishes.
The exact specialized order has colength one, so its action is
nonzero by the nodal calculation. The resonance is impossible.

Now \(a_6\) is constant, and \(j=0\) in (1) makes \(\Lambda'\)
polynomial, hence \(\Lambda\) is polynomial. If \(\Lambda\) were
constant then (1) would force \(a_1,\ldots,a_6=0\). With
\(L=\deg\Lambda\ge1\), the leading factor at step \(j<6\) is
\[
 \frac{j(L-1)}{2j+3}+L\left(\frac{1}{2j+3}-\frac1{15}\right)>0,
\]
so induction gives \(\deg a_j=j(L-1)\). Consequently \(L=1\)
and all coefficients of \(g\) are constant.

At the zero of \(\Lambda\), the six fixed reflected pair equalities
cut out the fixed order \(\mathbb C[T,V]\) by interpolation.
The exact specialized conductor gives colength six, so the actual
specialized order equals that fixed order. The only reduced
embedding type for \(\langle2,13\rangle\) is \((13,2)\) by the
one-place classification. Its degree-two generator is affine in
\(T\), and its degree-thirteen generator has odd part a nonzero
constant multiple of \(V\). Triangular changes give the canonical
pair. The canonical action vector is nonzero because
\(\sum_j a_j T^j/(2j+3)\) cannot be proportional to \(g\) by its
constant and leading coefficients. If all canonical actions
vanished, the six simple nonzero roots of \(g\) would imply that
proportionality. Thus the physical actions of the specialized
generators cannot all vanish. Equation (1)'s action identity
with scale \(\Lambda\) instead gives zero at the fixed endpoints.

## The cubic–septic branch

**Theorem ConductorTwelveCubicSepticExclusion.** The reduced
sequence \((7,3)\) is impossible under the stated hypotheses.

Use triangular target reductions, a source translation \(x-h(y)\)
with \(h\in K\), and monic target scalings. In the cubic generator
remove its constant term and depress it by the source translation.
In the septic generator subtract a multiple of the square of the
cubic generator to cancel the degree-six term, then a multiple of
the cubic generator and a constant to cancel degrees three and zero.
Let \(\lambda\in K^\times\) be the determinant of the resulting
target change. Choose \(\rho^{10}=\lambda\), work over
\(L=K(\rho)\), and set \(x-h(y)=\rho z\). Rescale the target
generators by \(\rho^{-3},\rho^{-7}\). The target determinant
becomes one and the normalized generators are
\[
 U=z^3+az,\qquad V=z^7+bz^5+cz^4+dz^2+ez .
\]
The original map in \((z,y)\) has Jacobian \(J\rho\).
The parameter derivative of the target automorphism composed with
its inverse is a polynomial target field of divergence zero, so
one can integrate it in the target variables to obtain
\[
 U_z V_y-U_y V_z=J\rho+\partial_z H(U,V).
 \tag{2}
\]
This is the [polynomial parameter identity](CONDUCTOR_EIGHT_CUBIC_QUINTIC_EXCLUSION.md#2-polynomial-parameter-compatibility-for-the-cubicquintic-pair)
with constant target determinant. The left side has degree at
most seven. Up to degree eight the basis of the order with
semigroup \(\langle3,7\rangle\) is \(1,U,U^2,V\), so these are
the only terms needed in \(H\).

Write primes for differentiation over \(K\) and its finite
extension. The coefficients of \(V,U^2,U\) in \(H\) come from
the terms in \(z^6,z^5,z^2\) and equal
\(3c'/7,\ (ab'-5ba')/6,\ (ad'-2da')/3\), respectively.
Then comparison in degrees \(7,4,3,1\) gives
\[
\begin{aligned}
 B&=b-7a/3,&D&=d-(4a/3+5B/7)c,\\
 E&=e+112a^3/81+10Ba^2/9-2c^2/7,&
 R&=7a^4/81+5Ba^3/27-(2a/7+15B/49)c^2-Ea-6Dc/7
\end{aligned}
\]
with \(B,D,E,R\in\mathbb C\). The constant field of the finite
extension with its extended derivation is still \(\mathbb C\).
The constant term of (2) gives
\[
 G=(4a^3/27+5Ba^2/21-3E/7)c+Da^2/3-2c^3/49,
 \qquad G'=J\rho .
 \tag{3}
\]

### Collision limits and coefficient poles

For a collision pair in \(z\), let \(S\) be the sum and \(p\)
the product of the endpoints. Dividing the difference of \(U\)
by their gap gives \(p=S^2+a\). Dividing that of \(V\) and
substituting for \(p\) gives
\[
 f(S)=S^6+(4a-b)S^4-cS^3+(3a^2-ab)S^2+(d-2ac)S-a^3+a^2b+e=0.
 \tag{4}
\]
The squared normalized gap is \(-3S^2-4a\), and the physical gap
is \(\rho\) times the normalized gap.

Suppose \(a\) has a pole at a finite place. The only possible
leading balance in the constant \(R\) is between \(7a^4/81\)
and \(-2ac^2/7\). Thus \(a,c\) have pole orders \(2q,3q\) with
\(q>0\) and leading coefficients \(A,C\) satisfying
\(C^2=49A^3/162\). Scaling the sum by pole order \(q\) gives the
monic initial polynomial
\[
 S^6+\tfrac53 AS^4-CS^3+\tfrac23 A^2S^2-\tfrac23 ACS+A^3/27.
\]
Its resultant with \(3S^2+4A\) is \(539A^6/3\) after the
leading relation. Hence every collision gap has pole order
\(q\) in the normalized variable. The leading term of \(G\)
is \(11a^3c/81\), of pole order \(9q\). At a finite place,
differentiation of a pole increases its order by the
ramification index. Thus (3) forces \(\rho\) to have a pole,
and every physical gap diverges, a contradiction.

If \(c\) has a pole while \(a\) is regular, scale the sum by
one third the pole order of \(c\). The initial polynomial is
\(S^6-CS^3+2C^2/7\), where \(C\ne0\) is its leading coefficient.
It has no zero at \(S=0\), so the normalized gaps again have
poles. Since \(G\) has leading term \(-2c^3/49\), \(\rho\)
has a pole and all physical gaps diverge. This is also
impossible. Consequently \(a,c\), and then \(b,d,e,G\), are
regular at all finite places. Equation (4) is monic and bounds
the normalized sums and gaps there. A zero of \(\rho\) would
force all physical gaps to vanish. Thus \(\rho\) has no finite
zeros.

### Finite poles of the source scale

A pole of \(\rho\) cannot occur at an unramified finite place
by (3) and regularity of \(G\). The Kummer extension \(L/K\)
has ramification indices dividing ten. Inertia acts by
multiplying \(\rho\) by roots of unity, so \(a,b,c,d,e\)
transform with weights \(2,2,3,5,6\), using inverse powers.
The invariants \(B,D,E,R\) have weights \(2,5,6,8\).
A constant with nontrivial weight for the inertia group must
vanish. A regular coefficient with nontrivial weight has
zero residue.

If the ramification index is two, \(D=0\) and \(c,d\) have
zero residue. Write \(a_*\) for the residue of \(a\).
A finite physical gap at a pole of \(\rho\) requires a common
limiting root of (4) and \(3S^2+4a_*\). Reduction of (4)
modulo the latter polynomial gives
\(E-28a_*^3/81-5Ba_*^2/9\). This must vanish, so the
coefficient of \(c\) in (3) vanishes at \(a_*\).
The invariant \(a-a_*\) has order at least two, and \(c\)
has positive order. Thus \(G\) has order at least three.
Its derivative is regular at a place of index two,
contradicting the pole of \(\rho\).

If the index is five, the residues of \(a,b,c,e\) vanish,
and \(B=E=R=0\). When \(D\ne0\) the limiting sum
polynomial is \(S^6+DS\). It has six simple roots. Only
the sum tending to zero allows a finite physical gap, and
both normalized endpoints in that pair tend to zero.
Every other pair has diverging physical gap. Neither
endpoint of such a pair can have a finite physical limit
by monicity of \(P\).

The pair with small sum is generically transverse. The
divided \(U\)-difference has partial derivative \(-1\) in
the product \(p\), and the resulting divided \(V\)-difference
has a simple root in \(S\). Sum and product are étale
coordinates for distinct ordered endpoints. Thus the
two branches have linearly independent tangent directions.
Neither endpoint has a third mate, since all other
collision pairs diverge. This is an ordinary node, whose
conductor has a simple zero at each endpoint. By
continuity of the roots of the exact (polynomial)
conductor near finite points, the specialized conductor
has degree two. It counts exactly these two finite
simple roots, which have distinct limiting physical
coordinates by the finite-pair argument.

The determinant of the normalized target change is one.
The physical action of the finite pair therefore equals
the \(U\,dV\) action at the normalized endpoints. It tends
to zero because these endpoints coalesce and the
normalized coefficients are regular. The specialized
order has colength one and identifies the two distinct
finite physical endpoints, so its action is nonzero by
the nodal calculation. This excludes a pole of index
five when \(D\ne0\).

### Homogeneous and remaining cases

At index ten all four constants vanish. The same holds
at index five if \(D=0\). Treat the case \(B=D=E=R=0\)
globally. The relation \(R=0\) becomes
\(a(7a^3/81-2c^2/7)=0\). If \(a=0\), choose
\(\tau^3=c\) and put \(w=\rho\tau\), so \(w^3\in K\)
is the coefficient of degree four in the unscaled
septic generator. Here \(c\ne0\) by (3). If \(a\ne0\),
put \(\tau=c/a\) and \(w=\rho\tau\). Then \(w\in K\)
is the ratio of the unscaled degree-four and linear
coefficients in the septic and cubic generators.
The two normalized sum equations after scaling by
\(\tau\) are exactly the initial equations from the pole
analysis. Hence all physical gaps are fixed nonzero
multiples of \(w\). The finite-unit-gap property makes
\(w\) a unit above every finite parameter. A power of \(w\)
is rational, so it is constant, and \(w\in\mathbb C^\times\).

The unscaled reduced generators now have constant
coefficients in \(u=x-h(y)\). Formula (3) gives
\(G=g_0\tau^9\) with \(g_0\in\mathbb C^\times\), so
\(G'=Jw/\tau\) implies \(\tau^{10}\) is nonconstant affine
in \(y\). In particular \(\lambda^{-1}=w^{-10}\tau^{10}\)
has a finite zero \(y_0\). The center \(h\) has no finite
poles, since the centered collision endpoints are fixed
and at least one physical pair has finite limits at
every finite parameter. Thus \(h\) is polynomial.

At \(y_0\) the actual specialized order is contained in
the fixed order of the unscaled reduced generators.
Indeed membership in each finite degree of that fixed
order is a system of constant linear equations on the
coefficients of \(P,Q\) in \(u\), which persists at
\(y_0\). Its colength is at most six by exact conductor
specialization. The fixed order has colength six, so the
two orders agree. The one-place classification gives
only reduced type \((7,3)\) for \(\langle3,7\rangle\).
The basis through degree seven is \(1,U_0,U_0^2,V_0\)
for the fixed unscaled pair \(U_0,V_0\). Thus every
specialized generating pair is related to the fixed
pair by polynomial target changes with nonzero
determinant, by reducing the two degrees and applying
a triangular change within this basis. A fixed canonical
action is nonzero: its generic physical action is
\(\lambda^{-1}\) times that fixed value and its
derivative is \(-J\) times the fixed nonzero gap.
The specialized physical actions all vanish by
polynomial continuation and the zero of \(\lambda^{-1}\).
This contradicts the specialized generator comparison.

Outside the homogeneous case, \(\rho\) has no finite
poles or zeros. Then the rational function
\(\rho^{10}\) is constant, so \(\rho\) is constant and
\(a,c,G\) are polynomials in \(y\). If \(a\) is
nonconstant, the same leading balance at infinity
gives pole order \(9q>1\) for \(G\). If only \(c\)
is nonconstant, the term \(-2c^3/49\) gives degree
greater than one. If both are constant, \(G'=0\).
Each case contradicts \(G'=J\rho\in\mathbb C^\times\).

## Exact controls

The [checker](verify_conductor_twelve_hyperelliptic_and_cubic_exclusions.py)
verifies the coefficient comparisons in (2)–(3), the
collision and inertia equations, the two initial
polynomials and the nonzero resultant in the balanced
pole case. It also checks the hyperelliptic resonance
list and the nodal action. Replay from the repository root.

~~~sh
.venv/bin/python research/hessian/hc4-triangular-graphs/verify_conductor_twelve_hyperelliptic_and_cubic_exclusions.py
~~~

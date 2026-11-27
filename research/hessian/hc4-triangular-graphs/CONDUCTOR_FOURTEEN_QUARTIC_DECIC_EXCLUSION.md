# The quartic–decic sequence \((10,4,7)\) in conductor degree fourteen

**Theorem (The \((10,4,7)\) conductor-fourteen exclusion).** Let
\((P,Q)\) be a complex polynomial plane Keller map in a source frame
\((x,y)\) with \(P\) monic of positive degree in \(x\). If its exact
coordinate-order conductor has global degree fourteen in \(x\), its
generic order cannot have the reduced characteristic sequence
\((10,4,7)\). The proof retains the exact specialized conductor and
allows its degree to drop on individual source lines.

## 1. The two Laurent coefficients and the exact parameter invariants

Write \(K=\mathbb C(y)\), and let \(J=\operatorname{Jac}_{x,y}(P,Q)\ne0\).
The
[coordinate-order theorem](GENERIC_COORDINATE_ORDER_CONDUCTOR_AND_CLOSED_COLLISIONS.md)
identifies the birational generic order with its finite specialized
orders. Its specialized conductor is nonzero and nonconstant on every
source line of a hypothetical noninvertible map. We use the
[characteristic-sequence reduction](CONDUCTOR_FOURTEEN_CHARACTERISTIC_SEQUENCES.md)
of this exact generic order. The
[Assi–García-Sánchez approximate-root theorem](sources/assi-garcia-sanchez-1407.0490v1.pdf),
Propositions 2 and 4, applies to its reduced generators of degrees
four and ten and gives a degree-seven approximate root. Here the
one-place hypothesis and the monic coefficient conditions follow
from the normalization by the source line as in the
[cubic–quintic proof](CONDUCTOR_EIGHT_CUBIC_QUINTIC_EXCLUSION.md).

Use the same square completion in the degree-ten generator as for
[sequence \((10,4,5)\)](CONDUCTOR_TWELVE_CABLE_EXCLUSIONS.md#4-sequence-1045).
Weights \(4i\), \(0\le i<5\), are distinct modulo ten, so the monic
degree-five approximate root in the quartic generator has no terms
of weight greater than twenty. After completing the square in the
decic and translating the quartic, normalize the actual generators
to the form
\[
 W=V^2-(U^5+\alpha U^3+\beta U^2+\gamma U+\delta),\qquad
 \deg W=7 ,
\tag{1}
\]
where the generators \(U,V\) are monic and the quartic is depressed
by a translation \(h(y)\in K\) of the longitudinal source coordinate.
Let \(\lambda\in K^*\) be the determinant of this target change.
As in the degree-twelve proof, choose \(\rho^{14}=\lambda\), write
\(x=h(y)+\rho z\), and rescale the pair by weights four and ten.
The resulting target determinant is one. Continue with the same
letters for the scaled generators and coefficients, writing
\[
 U=z^4+az^2+bz+c,\qquad \Delta=a^2-4c .
\]

Let \(S\) be the monic formal square root of the polynomial in \(U\)
in (1), expanded at \(z=\infty\). Since \(V+S\) has degree ten,
\(V-S=O(z^{-3})\) with a nonzero \(z^{-3}\) coefficient. In particular
\(V=[S]_{\ge0}\). If \(b=0\), then \(S\) and \(V\) are even in
\(z\), and \(W\) cannot have degree seven. Thus \(b\ne0\). The
coefficients of \(z^{-1}\) and \(z^{-2}\) in
\[
 S=U^{5/2}+\tfrac{\alpha}{2}U^{1/2}
                    +\tfrac{\beta}{2}U^{-1/2}+O(z^{-6})
\]
must both be zero. Solving these two coefficient identities gives
\[
 \alpha=-\frac5{64}(3\Delta^2-8ab^2),\qquad
 \beta=-\frac5{256}(\Delta^3+8ab^2\Delta-4b^4),\qquad
 [z^{-3}](V-S)=-\frac5{128}b^3\Delta .
\tag{2}
\]
In particular \(\Delta\ne0\). These identities and the divided
differences below are checked exactly by the
[companion checker](verify_conductor_fourteen_quartic_decic_exclusion.py).

Differentiating the parameter-dependent target automorphism as in the
degree-twelve proof gives the polynomial identity
\[
 U_zV_y-U_yV_z=J\rho+\partial_z H(U,V).
\tag{3}
\]
The contribution of \(S\) to its polynomial part is
\[
 \left[\,U_z(\alpha' U^{1/2}+\beta' U^{-1/2})/2\,\right]_{\ge0}.
\]
The contribution of \(V-S\) is a constant, four times the derivative
of its \(z^{-3}\) coefficient. The left side has degree at most five.
Because the semigroup is \(\langle4,7,10\rangle\), the evaluated \(H\)
lies in the span of \(1,U\) through degree six. The degree-five term
of (3) forces \(\alpha'=0\). The polynomial part of
\(U_z U^{-1/2}\) is \(4z\), so the degree-three and degree-one terms
then force the coefficient of \(U\) in \(H\) and \(\beta'\) both to
be zero. Using (2) at degree zero gives
\[
 \begin{split}
 A&=3\Delta^2-8ab^2\in\mathbb C,\\
 B&=4\Delta^3-A\Delta-4b^4\in\mathbb C,\\
 G&=b^3\Delta,\qquad J\rho=-5G'/32 .
 \end{split}
\tag{4}
\]
We separate the cases \((A,B)\ne(0,0)\) and \((A,B)=(0,0)\).

## 2. Collision equations and the finite-unit-gap condition

For a distinct matched pair in the \(z\)-coordinate let \(s\) be
its sum and \(g\) its gap. The quartic divided difference makes
\(s\ne0\) since \(b\ne0\), and gives the product
\((s^3+as+b)/(2s)\). Its substitution into the decic divided
difference gives \(s^2/32\) times the polynomial below:
\[
 \begin{split}
 F(s)={}&2s^7+8as^5-2bs^4+(5\Delta+8a^2)s^3\\
       &\quad -4abs^2+(10\Delta a-2b^2)s-5\Delta b=0,\\
 g^2={}&-s^2-2a-2b/s .
 \end{split}
\tag{5}
\]
Conversely, each root \(s\) of \(F\) at a generic parameter gives
a matched pair over an algebraic closure by the displayed product.
Indeed \(s\ne0\) since the constant term of \(F\) is nonzero.
If its pair had zero gap, both original derivatives along the
source line would be zero by the two divided differences, contrary
to the nonzero Jacobian. This converse will ensure the existence
of collapsing pairs in some parameter limits.

At any finite \(y_0\), at least one generic matched pair has
finite and distinct limits in the original \(x\)-coordinate,
after passing to its endpoint branches. As in the
[degree-twelve proof](CONDUCTOR_TWELVE_CABLE_EXCLUSIONS.md#1-common-source-and-target-identities),
the specialized conductor is nonconstant and nonzero, so a generic
conductor root reaches a finite root on that fibre. Its mate is
finite by monicity of the original \(P\). If they coincided in the
limit, the divided differences of \(P,Q\) would force
\(P_x=Q_x=0\) there. Thus at least one original gap \(\rho g\) is a
unit at every finite parameter. Likewise, even a single generic
matched pair cannot have both original endpoints finite and equal
in the limit.

## 3. Excluding finite coefficient poles when \((A,B)\ne(0,0)\)

At least one of \(\alpha,\beta\) is nonzero by (2) and (4).
Their scale weights are eight and twelve. The corresponding
unscaled coefficients belong to \(K\), and \(\rho^{14}\in K\),
so \(\rho^2\in K\). Work in \(K(\rho)\) and at a finite parameter
place of ramification index \(l\in\{1,2\}\). Use the normalized
local valuation \(v\) on this field, extending it to the roots
of (5) when necessary.

First suppose \(b\) or \(\Delta\) has a pole. Then
\(b^4\sim\Delta^3\) by (4). With \(k>0\) the pole orders of
\(b,\Delta,a\) are \(3k,4k,2k\), respectively. Here and below
\(\sim\) means equality of nonzero leading terms.
The auxiliary scale \(\sigma=\Delta/b\) has pole order \(k\).
Scale \(z=\sigma w\) and set
\(\widetilde U=\sigma^{-4} U(\sigma w)\) and
\(\widetilde V=\sigma^{-10} V(\sigma w)\). Their limits are
\[
 U_*=w^4+\tfrac38 w^2+w-\tfrac{55}{256},\qquad
 V_*=[U_*^{5/2}]_{\ge0}.
\]
Use the degree-seven function
\[
 \widetilde W=\sigma^{-20}(W+\gamma U+\delta)(\sigma w)
 =\widetilde V^2-\widetilde U^5
   -\alpha\sigma^{-8}\widetilde U^3-\beta\sigma^{-12}\widetilde U^2 .
\]
Its coefficients are bounded from this identity. Its \(w^7\)
coefficient is \(-5b^3\Delta\sigma^{-13}/64\) by (2), and tends
to \(-5/64\), since \(b^4\sim\Delta^3\) and \(\sigma=\Delta/b\).
The limiting \(U_*,V_*\) have derivative resultant \(208000\ne0\).

The coprime degrees four and seven of \(\widetilde U,\widetilde W\)
keep the matched endpoints bounded in the scaled variable \(w\).
In detail, if one endpoint escaped then both would, with ratio
a root of unity for both degrees by their equal values. The ratio
would be one, contradicting the leading quartic divided difference.
Since the limiting \(U_*,V_*\) have no common critical point,
a bounded scaled pair cannot have gap tending to zero.
Thus every normalized gap in \(z\) has pole order \(k\).
But \(G=b^3\Delta\) has pole order \(13k\), and its derivative with
respect to \(y\) has pole order \(13k+l\). Equation (4) makes
every physical gap diverge. The finite-unit-gap condition forbids this.

The only remaining possible pole is a pole of \(a\) with \(b\to0\)
and \(\Delta\to D\) finite. Put \(k=v(b)>0\).
If \(3D^2-A=0\) and \(D\ne0\), the derivative of
\(4\Delta^3-A\Delta-B\) at \(D\) is \(9D^2\ne0\).
The second invariant gives \(\Delta-D=O(b^4)\), so the first
gives \(a=O(b^2)\), contrary to a pole. If \(3D^2-A=0=D\),
then \(A=B=0\), excluded in this section. Thus \(v(a)=-2k\).

When \(D\ne0\), the valuations of the coefficients of \(F\), by
degree \(0,\ldots,7\), are
\[
 k,\ -2k,\ -k,\ -4k,\ k,\ -2k,\ +\infty,\ 0 .
\]
The lower hull has segments from degrees zero to one, one to
three, and three to seven. This gives one sum of order \(3k\),
two of order \(k\), and four of order \(-k\), with multiplicities.
By (5), all gaps have order at least \(-k\). The last face of
\(F\) is \(2s^3(s^2+2a)^2\) in leading terms. Hence for every one
of those four large sums the two leading parts of \(-s^2-2a\)
cancel, while \(b/s\) has positive order. Their gaps have order
strictly greater than \(-k\). Equation (4) gives
\(v(\rho)=3k-l\). If \(2k>l\), all physical gaps collapse, contradicting
the unit gap. The only other case is \(k=1,l=2\). Then all
\(\rho z\) endpoints are bounded and the unit pair forces \(h\)
to be regular at this place. The four large-sum pairs have gaps
tending to zero with finite endpoints in \(x\), again impossible.

When \(D=0\), we have \(B=0,A\ne0\) and \(\Delta\sim-4b^4/A\).
The coefficient valuations are now
\[
 5k,\ 2k,\ -k,\ -4k,\ k,\ -2k,\ +\infty,\ 0 .
\]
The two segments, in degrees zero to three and three to seven,
give three sums of order \(3k\) and four of order \(-k\).
Again every gap has order at least \(-k\).
Here \(v(G)=7k\) and \(v(\rho)=7k-l\), so all physical gaps
collapse. This proves that \(a,b,\Delta\) have no finite poles.

## 4. The parameter scale cannot vanish or have a pole

Still assume \((A,B)\ne(0,0)\). Since \(G\) is regular at every
finite place, a pole of \(\rho=-5G'/(32J)\) would have to occur at
a ramified place \(l=2\). The inertia there sends \(\rho\) to
\(-\rho\), and hence sends \(b\) to \(-b\) because the unscaled
coefficient has weight three. Thus \(b\) vanishes at the place.
This gives \(v(G)=v(b^3\Delta)\ge3\), making \(G'\) regular and
contradicting the pole.

Suppose now that \(\rho\) vanishes at a finite place. If \(\Delta\)
is a unit then \(F\) makes all sums bounded, and \(F(s)/s=0\)
bounds \(b/s\) since its coefficient is \(-5\Delta\) and the
other terms are polynomial in the bounded coefficients and \(s\).
Thus all normalized gaps are bounded and all physical gaps
collapse. If \(\Delta\) tends to zero, then \(b\) must be a unit:
otherwise \(B=0\), \(A\ne0\) and \(a=(3\Delta^2-A)/(8b^2)\)
would have a pole, already excluded. Inertia is then trivial, so
\(l=1\). Write \(q=v(\Delta)>0\). We have \(v(\rho)=q-1>0\).
At the limit \(F\) has a unique simple root at zero, since its
linear coefficient there is \(-2b^2\). The corresponding branch
has \(s\sim-5\Delta/(2b)\) and \(g^2\sim4b^2/(5\Delta)\).
All other roots have nonzero finite limits and give bounded
normalized endpoints. If \(q>2\), even the unbounded normalized
pair has physical gap tending to zero, a contradiction. If \(q=2\),
all \(\rho z\) endpoints are bounded, and the unit pair forces
\(h\) to be regular. Every other matched pair then has both
endpoints finite and equal in the limit, again a contradiction.

It follows that \(\rho^2\in K\) has neither zeros nor poles at
finite \(y\), so \(\rho\) is constant. Now \(b,\Delta,a\in K\)
have no finite poles and are polynomials. If \(b\) or \(\Delta\)
is nonconstant their common pole at infinity has orders \(3k,4k\)
as above, making \(\deg G=13k>1\). If they are constant, then so
is \(G\). Both contradict \(G'\ne0\) constant.

## 5. The homogeneous invariants and the fixed-order determinant

Finally let \(A=B=0\), and put \(\tau=\Delta/b\).
Then \(b=\tau^3\), \(\Delta=\tau^4\), \(a=3\tau^2/8\) and
\(c=-55\tau^4/256\). In \(z=\tau w\) the two generators are
precisely the fixed shape \(U_*,V_*\) from Section 3 multiplied by
\(\tau^4,\tau^{10}\), respectively. Consequently every physical gap
is a fixed nonzero multiple of \(w_0=\rho\tau\). This element
belongs to \(K\): it is the ratio of the unscaled \(\Delta,b\)
coefficients. The unit-gap property makes \(w_0\) a nonzero constant.
Now \(G=\tau^{13}\), so (4) gives \(\tau^{14}\) affine nonconstant
in \(y\). The inverse target determinant
\(\lambda^{-1}=\tau^{14}/w_0^{14}\) has a finite zero \(y_0\).

In the coordinate \(u=x-h(y)\) the unscaled reduced generators
are fixed polynomials, denoted \(U_0,V_0\). The shift \(h\) is
polynomial, because otherwise all original matched endpoints
would escape at its pole while the shifted endpoints are fixed
finite numbers. The actual \(P,Q\) then belong to
\(\mathbb C[y]\otimes A_*\), where \(A_*=\mathbb C[U_0,V_0]\).
This follows by projecting their coefficients to
\(\mathbb C[u]/A_*\). At any specialization their generated order
is contained in \(A_*\), and by the exact conductor theorem its
colength is at most seven. Since \(A_*\) has colength seven,
the specialized order equals it.

Use the
[singular tangent determinant argument](HYPERELLIPTIC_CONDUCTORS_AND_SINGULAR_TANGENT_DETERMINANT.md)
in this fixed order. Lift \(P,Q\) as polynomials over \(\mathbb C[y]\)
in two variables representing \(U_0,V_0\), which is possible by the
coefficientwise inclusion just established. Over \(K\) these lifts
differ from the two components of the inverse target automorphism
by multiples of the fixed irreducible plane equation of \(A_*\).
At any affine singular point of this fixed curve, that equation and
both its partial derivatives vanish. Such a point exists because the
order is not normal. The determinant of the regular lifts at this
point therefore equals \(\lambda^{-1}\) as a rational function of
\(y\). At \(y_0\) the lifts induce the specialized order isomorphism.
The source and image local cotangent spaces both have dimension
two there, and each is its ambient cotangent space because the
respective curve equation has no linear term. Hence the determinant
of the regular lifts is nonzero at \(y_0\), contradicting the zero
of \(\lambda^{-1}\). This excludes the homogeneous case and proves
the theorem.

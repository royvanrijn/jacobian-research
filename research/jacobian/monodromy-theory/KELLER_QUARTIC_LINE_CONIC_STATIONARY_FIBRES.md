# Stationary fibres of a quartic double line and conic

**Theorem KellerQuarticTransverseLineConicExclusion.** In a
noninvertible three-variable Keller map of degree at most five over
a field of characteristic zero, a quartic component cannot have
leading form \(z^2 q\) where \(q\) is a geometrically irreducible
conic transverse to \(z=0\).

**Proposition KellerQuarticTangentLineConicScalarRestriction.**
For a noninvertible map in the same setting, if \(R_4=z^2 q\)
with geometrically irreducible \(q\) tangent to \(z=0\), then
there are affine changes over the algebraic closure making the
component
\[
                       R=z+z^2(x^2+zy+c)
\]
for a scalar \(c\). This is a necessary form only.

## The stationary value identity

Work over \(\mathbb C\) by embedding the finitely generated
coefficient field. The
[fixed-component double-line restriction](KELLER_QUARTIC_DOUBLE_LINE_IRREDUCIBLE_CONIC.md)
gives \(z\mid R_3\). Hence in a frame preserving \(z\) the
component is quadratic in \(x,y\) jointly. Write
\[
                 R=A(z)x^2+B(z)xy+C(z)y^2+d(z)x+e(z)y+f(z).
\]
If its discriminant is nonzero, the
[square discriminant proof](QUINTIC_KELLER_BINARY_QUARTIC_REPEATED_CUBIC_VANISHING.md#the-discriminant-of-a-quadratic-scalar-fibre)
gives \(B^2-4AC=S^2\) with \(S\) polynomial. Here \(S\) can
be taken polynomial because its square is polynomial. In this
square case the coefficient calculation in that proof applies
over \(\mathbb C(z)\) itself. In detail, center the binary
quadratic at its stationary point over \(\mathbb C(z)\) and
diagonalize it as \(R=UV+K(z)\), where
\[
             K=f-\frac{Cd^2-Bde+Ae^2}{4AC-B^2}.
\]
The linear change to \(U,V\) has determinant \(\pm S\).
On \(U\ne0\), the exact fibre residue is
\(\pm dz\wedge dU/(S U)\).
The \(dU/U\) coefficient of the regular Keller primitive is a
rational function of \(z\), polynomial in the fibre value \(t\),
whose derivative in \(z\) at fixed \(t\) is \(1/S\) after a choice
of sign. Evaluating this coefficient on the \(V=0\) line in
the fibre \(t=K(z)\) makes it vanish because the primitive is
regular at \(U=0\). Thus there is a nonconstant polynomial \(H\)
such that
\[
                         \frac{d}{dz}H(K(z))=\frac{1}{S(z)},
 \tag{1}
\]
absorbing the overall sign in \(H\).
The primitive is polynomial in \(t\) before specialising,
since it comes from the two polynomial mates and the substitution
\(V=(t-K(z))/U\). Thus (1) cannot be obtained with a rational
\(H\) having a pole in the finite \(t\)-plane.

## Transverse conic

Normalize \(q=xy+\gamma z^2\) with \(\gamma\ne0\).
The degree bounds give \(\deg A,\deg C\le1\),
\(B=z^2+O(z)\), \(\deg d,\deg e\le2\), and
\(f=\gamma z^4+O(z^3)\).
The discriminant therefore has degree four, so \(S\)
has degree two. The correction term to \(f\) in \(K\) grows
at most quadratically at infinity, since its numerator has
degree at most six and its denominator degree four.
Thus \(K\) has a pole of order four at infinity.
A nonconstant polynomial in \(K\) also has a pole there,
but its derivative cannot equal \(1/S=O(z^{-2})\).
This proves the transverse exclusion.

## Tangent conic

Normalize \(q=x^2+zy\). Now \(A=z^2+O(z)\),
\(\deg B,\deg C\le1,\ \deg d\le2\),
\(e=z^3+O(z^2)\), and \(\deg f\le3\).
A nonzero discriminant has degree at most three, so its
polynomial square root has degree at most one. A nonconstant
linear \(S\) would give a nonzero residue of \(1/S\) at its
root, contrary to (1), which is a rational derivative.
If \(S\) is a nonzero constant, then the term \(Ae^2\) has
degree eight and exceeds the degrees of all other terms
of the numerator of \(K-f\). Hence \(K\) has degree eight,
and (1) cannot hold with a constant nonzero right side.
Thus the discriminant vanishes identically.

From \(B^2=4AC\) and the degree bounds, \(C\) is constant.
If \(C\ne0\), replace \(y\) by the polynomial coordinate
\(Y=y+B(z)x/(2C)+e(z)/(2C)\), preserving \(x,z\).
The scalar becomes \(C Y^2+D(z)x+F(z)\), with
\(D=d-Be/(2C)\) of degree four, since \(B\) then has
degree one. At any root of \(D\), the restriction
of the remainder to that root is a nonconstant
quadratic in \(Y\). This contradicts the degree-independent
[parallel coefficient residue proposition](KELLER_LINEAR_SOURCE_COEFFICIENT_CUBIC_CLASSIFICATION.md#normal-residues-along-a-repeated-line).
The shear preserves polynomiality of the Keller mates,
and no degree bound is needed for that proposition.

It follows that \(C=B=0\), and \(R=e(z)y+A(z)x^2+d(z)x+f(z)\).
At every root of \(e\), the same residue proposition makes the
restriction constant or affine linear in \(x\), so \(A\)
vanishes there. If all restrictions were affine with nonzero slope,
[parallel coefficient completion](KELLER_LINEAR_SOURCE_COEFFICIENT_CUBIC_CLASSIFICATION.md#coordinates-on-the-parallel-support-fibres)
would give inversion in the original degree bounds.
Hence some root \(\xi\) of \(e\) has constant restriction.
The residue proposition then gives
\(A'(\xi)=d'(\xi)=0,\ f'(\xi)\ne0\).
Therefore \(A=(z-\xi)^2\) and \(e=(z-\xi)^3\), since all roots
of \(e\) are roots of \(A\). Also \(d\) is a scalar multiple
of \((z-\xi)^2\). Translate \(z\), complete the square by a constant
translation in \(x\), and absorb the pure cubic into \(y\).
After rescaling the component and the binary variables to make
its linear \(z\)-coefficient one, the remaining scalar is
\(z+z^2(x^2+zy+c)\). This is the asserted form.

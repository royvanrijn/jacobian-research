# The position double locus in the order-\((12,42)\) transfer

In the unsplit order-\((12,42)\) planar profile, the graph lift \((x,P,Q)\) has a double locus over the general position. Its divided differences give two resultants, one in the difference of the two momenta and one in the first momentum. Their leading coefficients track the second resultant and its residue primitive. They determine the determinant of the difference on the double-point algebra, its exceptional ends away from the high roots, and a failure of saturation for the graph lift.

Work over \(\mathbb C\). As in the [target-cylinder calculation](POLYNOMIAL_WEYL_ORDER_TWELVE_FORTY_TWO_TARGET_CYLINDER.md), put \(H=2m\), \(M=5m\), and let \(h=h_0 f^{15}\), with \(f\) the monic squarefree product of the high roots. Write \(B=f g\) for the residue primitive. The [residue equation](POLYNOMIAL_WEYL_ORDER_TWELVE_FORTY_TWO_RESIDUE_MONODROMY.md#the-residue-cover-and-its-blocks) gives \(\deg g=3m\), with simple roots disjoint from those of \(f\). Denote the leading \(X\)-coefficients of \(P,Q,W,Z\) by \(a,b,c,d\), respectively:
\[
             a=h^2,\qquad b=h^7,\qquad c=\delta h^{17/3},\qquad d=K B h^{49/3}.
\]

For a polynomial \(T(x,X)\) set \(T[X,Y]=(T(x,Y)-T(x,X))/(Y-X)\). For fixed \(x\), two distinct momenta represent the same image precisely when both \(P[X,Y]\) and \(Q[X,Y]\) vanish. On the diagonal these are \(P_X,Q_X\), which have no common affine zero by the Keller identity. These are the usual relative double-point equations for a parametrized family. [Nuño-Ballesteros and Peñafort Sanchis, Section 3](https://arxiv.org/pdf/1509.04990)

## The two resultants and the difference norm

**Theorem PlanarFortyTwoDoublePointResultants.**
Over \(\mathbb C(x)\), both
\[
       D_x(s)=\operatorname{Res}_X(P[X,X+s],Q[X,X+s]),\qquad
       C_x(X)=\operatorname{Res}_Y(P[X,Y],Q[X,Y])
\]
have degree \(198\) in their displayed variables. Their leading coefficients are nonzero constant multiples of \(a^{41} b c^2 d\), or equivalently of \(f^{1751}g\). The ordered double-point algebra
\[
            \mathbb C(x)[X,Y]/(P[X,Y],Q[X,Y])
\]
has dimension \(198\). The norm of \(Y-X\) in this algebra is
\[
                              C\, f^{-498} g^{-1}
                 \qquad(C\in\mathbb C^\times).
\]

For the first resultant put \(X=s u\) and let \(s\) tend to infinity. The initial roots of \(s^{-11}P[X,X+s]\) are the eleven simple roots \(u=(\zeta-1)^{-1}\), where \(\zeta^{12}=1\), \(\zeta\ne1\). These lift uniquely in \(1/s\). On each lift \(P(Y)=P(X)\), with \(Y=X+s\). Use the degrees \(12,42,34,99\) along the resulting initial ratio \(Y/X=\zeta\).

For the six roots with \(\zeta^{42}\ne1\), \(Q[X,Y]\) has degree \(41\) in \(s\) and leading coefficient a nonzero constant times \(b\). For the four roots with \(\zeta^{42}=1\), \(\zeta^{34}\ne1\), use
\[
                  Q(Y)-Q(X)=\frac{W(Y)-W(X)}{Q(Y)+Q(X)}.
\]
Here the denominator has degree \(42\) and the numerator degree \(34\). The divided difference therefore has degree \(-9\), with leading coefficient proportional to \(c/b\).

The remaining root is \(\zeta=-1\). Write subscripts \(1,0\) for evaluations at \(Y,X\). Let \(\eta\) denote the constant multiplying \(WP^2Q\) in the second resultant. At a common \(P\)-value \(p\),
\[
       Z_1-Z_0=(Q_1-Q_0)
       \bigl[(W_1^2+W_1W_0+W_0^2-\eta p^2 Q_1)(Q_1+Q_0)
                    -\delta^3p^5-\eta p^2 W_0\bigr].
\]
The bracketed factor has degree \(110=68+42\) and leading coefficient proportional to \(c^2 b\), since the two leading \(W\)-values and the two leading \(Q\)-values coincide at \(\zeta=-1\). The left side has degree \(99\), which is odd. Thus the divided difference of \(Q\) has degree \(-12\) and leading coefficient proportional to \(d/(c^2 b)\).

In the resultant product, the leading coefficient of \(P[X,X+s]\) in \(X\) is \(12a\). The degree in \(s\) is \(6\cdot41-4\cdot9-12=198\). The coefficient product is a nonzero constant times
\[
             a^{41}b^6(c/b)^4d/(c^2 b)=a^{41} b c^2 d.
\]
For \(C_x(X)\), instead fix \(X\) large and solve the first divided difference in \(Y\). The same eleven ratios occur and \(Y-X\) has degree one in \(X\) at every ratio. This gives the same degree and coefficient up to a nonzero constant. The coefficient product equals a scalar times \(h^{350/3}B=f^{1751}g\).

Over \(\mathbb C(x)[s]\), the first divided difference in \(X\) is monic after a unit rescaling. Multiplication by the second in this free algebra has determinant proportional to \(D_x(s)\). Smith normal form over \(\mathbb C(x)[s]\) shows that the double-point algebra has length \(\deg_s D_x\), and that the characteristic polynomial of multiplication by \(s=Y-X\) is \(D_x(s)\) divided by its leading coefficient. The value \(D_x(0)\) is the [polar resultant](POLYNOMIAL_WEYL_ORDER_TWELVE_FORTY_TWO_TARGET_CYLINDER.md#the-polar-norm), a scalar times \(f^{1253}\). Division gives the norm identity.

If \(A^+\) is the exchange-invariant subalgebra, \(Y-X\) is an invertible anti-invariant element. The ordered algebra splits over \(A^+\) as \(A^+\oplus (Y-X)A^+\), and \(A^+\) has dimension \(99\) over \(\mathbb C(x)\). The norm of \((Y-X)^2\) from \(A^+\) has squareclass \(g\), since it agrees with the ordered norm of \(Y-X\) up to sign.

## Ends at the roots of the residue numerator

**Corollary PlanarFortyTwoNonhighDoubleEnds.**
At a finite position outside \(f=0\), pairs of distinct momenta with common \(P,Q\) values can escape to momentum infinity precisely at the roots of \(g\). At each root of \(g\) there is exactly one ordered branch at infinity. It is reduced and smooth in the chart with parameter \(1/(Y-X)\). Exchange of the momenta acts by negation of this parameter. Its projection to the position line has even ramification index. The total multiplicity of the pole of the difference over such a position, as recorded by its norm, is one.

Where \(f\ne0\), the leading coefficient of \(P\) is a unit. If one momentum escapes while \(P(X)=P(Y)\), then both escape with a limiting ratio among the twelfth roots of unity. Ratio \(1\) is impossible for a distinct pair. The divided difference of \(P\) there has its leading nonzero coefficient \(12a\) after dividing by \(X^{11}\), even when \(Y-X=o(X)\). Thus \(Y-X\) tends to infinity and the above eleven ratio charts account for all possibilities. Their respective leading coefficients in \(Q[X,Y]\) are units at the position except possibly at \(\zeta=-1\), and at that ratio the leading coefficient vanishes precisely at a root of \(g\).

Near \(\zeta=-1\) use \(w=1/(Y-X)\). The first divided difference gives a unique analytic \(u=w X\) with \(u(x,0)=-1/2\). Its equation is regular and has simple \(u\)-derivative at the limiting root. In this chart the second equation is equivalent, for \(w\ne0\), to
\[
                  \Xi(x,w)=w^{99}(Z(x,Y)-Z(x,X))=0,
\]
by the bracketed factor of degree \(110\), which is a unit after multiplication by \(w^{110}\). At \(w=0\), \(\Xi\) is a nonzero constant times \(d\). Its position derivative is nonzero at each root of \(g\). Exchange of \(X,Y\) sends \(w\) to \(-w\) and preserves this ratio chart by uniqueness. The odd degree \(99\) makes \(\Xi\) invariant. The implicit equation therefore gives exactly one smooth branch \(x=\chi(w^2)\).

This branch cannot have constant \(x\). On any fixed finite position line the polynomial pair in the momentum is birational onto its image. Otherwise Lüroth's theorem and the unique pole of a nonconstant coordinate give a common polynomial right factor of degree greater than one. More explicitly, the intermediate smooth projective rational curve has one point below momentum infinity, with infinity its entire preimage. Taking that point as infinity makes the intermediate parameter polynomial. A critical point of that polynomial makes both momentum derivatives vanish, contrary to the Keller identity. Thus a curve of distinct double pairs at one fixed position is impossible. The nonconstant \(\chi(w^2)\) gives even position ramification, and the simple pole of \(Y-X\) on this branch agrees with the norm at \(g=0\).

## Saturation for the projection by position

**Theorem PlanarFortyTwoPositionLiftUnsaturated.**
For the lift \((x,P,Q)\), the inverse image of the singular locus of its image closure has a curve component dominating the position line whose image under \((P,Q)\) is unsaturated. That is, some such irreducible component \(L\) has \((P,Q)^{-1}((P,Q)(L))\ne L\) in the affine source.

Remove from \(C_x(X)\) its content as a polynomial in \(X\) over \(\mathbb C[x]\), and call the resulting polynomial \(C_{\mathrm{hor}}\). The content is supported on \(f=0\). Indeed at any other finite position the first divided difference retains its invertible leading coefficient in \(Y\), and birationality on that position line makes the double locus finite. In particular, the \(X\)-leading coefficient of \(C_{\mathrm{hor}}\) still has a simple zero at each root of \(g\).

Over \(f\ne0\) the lift is finite onto its image closure, since \(X\) is integral over the subalgebra containing \(x\) and \(P\). It is generically birational by the fixed-position argument, or by the [position-field computation](POLYNOMIAL_WEYL_ORDER_TWELVE_FORTY_TWO_PLANAR_FIBRES.md#monodromy-of-the-generic-position-projection). The smooth source there is therefore its normalization. An image point with two preimages is singular. At a point with one preimage the inverse function theorem in \((x,P)\) or \((x,Q)\) shows that the finite image is locally a smooth graph, since at least one of the two momentum derivatives is nonzero. Thus the horizontal components of the inverse image of the singular locus are exactly the zero components of \(C_{\mathrm{hor}}\).

Suppose all those components are saturated. Each of their irreducible images has curve closure by quasi-finiteness. The pullback of a reduced equation of that closure has no other irreducible component, since any other component would meet the preimage of the open dense subset in the image. It is reduced at its generic point by the Keller condition. Consequently each factor of \(C_{\mathrm{hor}}\), and hence \(C_{\mathrm{hor}}\) itself, is a polynomial in \(P,Q\), up to constants.

The [graded order algebra](POLYNOMIAL_WEYL_ORDER_TWELVE_FORTY_TWO_PLANAR_FIBRES.md#the-shared-planar-resultant) gives a different leading coefficient for such a polynomial. In order \(198\), the exponent of \(Z\) in a leading monomial is zero or two, by parity. For exponent zero the leading coefficient is proportional to \(h^{33}\); exponent two allows only \(Z^2\), with coefficient proportional to \(B^2h^{98/3}\). Therefore the order-\(198\) coefficient has the form
\[
                 h^{33}\left(\lambda+\mu B^2/h^{1/3}\right),
                    \qquad \lambda,\mu\in\mathbb C .
\]
At a root of \(g\) its order is zero or two, whereas the leading coefficient of \(C_{\mathrm{hor}}\) has order one. The contradiction proves unsaturation. This is the condition on the extra coordinate \(x\) from [Borisov–Gabber–Vasiu, Definition 22.2 and Proposition 27.5](https://arxiv.org/pdf/2609.05746v1). Other possible extra coordinates are outside this calculation.

## Specialization in the outer chart

**Corollary PlanarFortyTwoOuterDoubleAlgebra.**
In the outer chart with \(r=1/x\), the relative ordered double-point algebra is finite free of rank \(198\) over \(\mathbb C[[r]]\). The momentum-coordinate difference \(v_2-v_1\) has determinant of order \(9m\) at \(r=0\). Its vanishing quotient is supported on the common derivative roots \(P_0',Q_0'\) and has total length \(9m\).

Apply the same resultant calculation to the polynomials in the two outer coordinates \(v_1,v_2\). The leading coefficients of \(P,Q,W,Z\) in the outer coordinates are units over \(\mathbb C[[r]]\). The first divided difference gives a monic polynomial over the position ring with the difference adjoined, and the determinant of multiplication by the second gives a degree-\(198\) polynomial with unit leading coefficient. This makes the relative algebra finite over \(\mathbb C[[r]]\). Its generic and special dimensions are both \(198\) by the same calculation, so it is free.

The relation \(X=\Lambda(r)+r^M v\) gives \(v_2-v_1=r^{-M}(Y-X)\). Multiplicativity of determinants makes its norm \(r^{-198M} C f^{-498}g^{-1}\), of order \(498(2m)+3m-198(5m)=9m\). Inserting \(v_1=v_2\) into the divided differences gives \(P_v,Q_v\), so the valuation of the determinant is the length of their common-zero scheme. Support and multiplicity agree with the [local sheet excess](POLYNOMIAL_WEYL_ORDER_TWELVE_FORTY_TWO_TRIANGULAR_FIBRES.md#local-sheet-counts-on-the-outer-curve).

## Vertical content at the high roots

**Corollary PlanarFortyTwoDoubleResultantContent.**
At a high position \(b_*\), let \(d_0,d_1\) be the critical \(P\)-values of the two quadratic characteristic charts, at the zero and nonzero simple roots respectively. Put \(\epsilon_{b_*}=1\) when \(d_0=d_1\) and zero otherwise. The content of \(C_x(X)\) has order \(1252+\epsilon_{b_*}\) at \(b_*\). Thus the \(X\)-leading coefficient of \(C_{\mathrm{hor}}\) is a nonzero constant times
\[
                       g(x)\prod_{b_*}(x-b_*)^{499-\epsilon_{b_*}} .
\]

Evaluate the order of \(C_x(X)\) in \(x-b_*\) with \(X\) generic and finite. This computes the content order. In the resultant product solve \(P(x,Y)=P(x,X)\), dropping \(Y=X\). The [high-root fibre charts](POLYNOMIAL_WEYL_ORDER_TWELVE_FORTY_TWO_PLANAR_FIBRES.md#the-complete-pole-profile-over-a-high-root) give one other affine point of \(Y\), two points on the nonzero simple-root chart and eight sheets in the scalar chart. The conclusions of that chart calculation also apply to the moving value \(P(x,X)\), by lifting the generic simple initial roots with this analytic value on the right.

At the other affine point the difference of the \(Q\)-values and the difference of the momenta are units for generic \(X\). Indeed the opposite root of the quadratic \(P\)-equation gives the negative of the same odd \(Q\)-polynomial in the normalized parameter. On the eight scalar sheets the \(Q\)-value has order \(-1\) in the base parameter \(x-b_*\), by the two stated pole orders and position orders. The denominator \(Y-X\) has order \(-3\) on each. Thus they contribute \(8\cdot2=16\).

For the two nonzero simple-root points \(Y-X\) also has order \(-3\), while the \(Q\)-values are finite. On either quadratic chart the image restrictions are the [polynomials of degrees two and seven](POLYNOMIAL_WEYL_ORDER_TWELVE_FORTY_TWO_PLANAR_FIBRES.md#the-simple-characteristic-roots) with parameter \(v\), \(P=v^2+d_i\), and a prescribed odd \(Q\)-polynomial depending only on \(d_i\). They are birational onto their images by the coprime degrees, and the normalized \(P\)-map has unique critical value \(d_i\). Hence the two image curves differ when \(d_0\ne d_1\), in which case neither of the \(Q\)-values on the nonzero chart can match the affine one at the generic level.

When \(d_0=d_1\), exactly one of the two values matches, with equal normalized parameters \(v\). Along a fixed \(P\)-level the first base derivative of \(Q\) in a chart at a simple root \(\rho\) is a fixed sign times \(\kappa/(2v r_{b_*}'(\rho))\). This follows from \(P=v^2+C_{b_*}(x-b_*)\) and the coordinate change \(x_{b_*}^3v=u(x_{b_*},x_{b_*}^3 Y)\), whose rescaled transverse derivative on the boundary is \(r_{b_*}'(\rho)\). Along the moving level the extra derivatives through \(P\) cancel between the two matching functions, since the restricted \(P,Q\) formulas are the same at \(d_0=d_1\). The fixed-level slopes differ by the [first-face derivative comparison](POLYNOMIAL_WEYL_ORDER_TWELVE_FORTY_TWO_PLANAR_FIBRES.md#monodromy-of-the-generic-position-projection). The matching difference therefore has order one. This gives total order \(6+\epsilon_{b_*}\) from the two simple-root points.

The leading coefficient of the degree-eleven divided difference is a scalar times \(h^2\), of order thirty here. Its power \(41\) contributes \(1230\) to the resultant. The resulting content order is \(1230+16+6+\epsilon_{b_*}=1252+\epsilon_{b_*}\). Comparing with the leading coefficient \(f^{1751}g\) gives the formula.

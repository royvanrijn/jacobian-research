# The high-root ends of the position double locus

The position double locus of the unsplit order-\((12,42)\) planar transfer has length \(198\) over the generic position. At each high root the leading polynomial pairs on the first and scalar faces identify its branches. The two quadratic charts contribute the remaining branches. In particular, the high roots contribute no point fixed by exchange of the momenta.

Use the [position double equations](POLYNOMIAL_WEYL_ORDER_TWELVE_FORTY_TWO_POSITION_DOUBLE_LOCUS.md), and write \(y=x-b\) at a high root. A sheet count includes scheme multiplicity when the generic double algebra is considered before reduction. We give projection indices on the reduced normalization for the explicitly separated branches.

## The first and scalar faces

**Proposition PlanarFortyTwoHighDoubleLeadingPairs.**
The normalized first-face pair
\[
 r(z)=z(z-1)^4(z-q),\qquad A(z)=z(z-1)(z-q)(z-d),\qquad
 2q^2-5q+5=0,\qquad d=(5-q)/3,
\]
has six simple ordered double points with common \(r\)-value nonzero. At the zero \(r\)-value it has the six simple ordered pairs of distinct zeros of \(r\). The normalized scalar pair
\[
 U(z)=z(1+z^7),\qquad V(z)=G(z^7),\qquad
 G(S)=35+280S+560S^2+448S^3+128S^4
\]
has \(84\) simple ordered double points, all with both coordinates nonzero. The finite equation \(G(S)^2=128^2 S(1+S)^7\) has three simple solutions, all disjoint from \(S=0,-1\) and \(G(S)=0\).

The first characteristic roots are \(0,1,q\) after scaling. The leading \(W\)-polynomial on this face is a scalar times \(r^5 A\). Indeed it has the indicated simple factors in \(A\), and the [first-face equation](POLYNOMIAL_WEYL_ORDER_TWELVE_FORTY_TWO_LEADING_FIELDS.md#the-two-high-root-newton-faces) reduces to \(3r A'-2r'A=C r\) after the common factor is removed. The coefficient one order below the top gives \(d=(5-q)/3\), and the exact identity with \(C=-qd\) follows modulo the quadratic of \(q\).

For the divided differences of \(r,A\), elimination of the second parameter gives a nonzero constant times \(z^2(z-1)^2(z-q)^2\) times a squarefree sextic prime to \(r\). The two divided differences and their Jacobian determinant generate the unit ideal over \(\mathbb Q(q)\). Thus all intersections are simple. At \(r=0\), the only points are the six ordered pairs across the three distinct roots. The [exact checker](verify_weyl_forty_two_high_double.py) verifies the resultant, squarefreeness, and the unit ideal over the quadratic field.

For the scalar pair, the [scalar successor equation](POLYNOMIAL_WEYL_ORDER_TWELVE_FORTY_TWO_PLANAR_FIBRES.md#the-complete-pole-profile-over-a-high-root) gives \(G\) by coefficient comparison after \(A_0\) and \(S\) are scaled. A double point cannot have a zero parameter, since \(G(0)\ne G(-1)\). For an ordered pair \(z,z'\) set \(\rho=z'/z\). No nontrivial eighth root of unity can occur from \(U(z)=U(z')\). Away from \(\rho=1\) this equation solves \(S=z^7=-1/H(\rho)\), with \(H(\rho)=1+\rho+\cdots+\rho^7\). The numerator of \(G(\rho^7 S)-G(S)\), after multiplication by \(H^4\), factors as
\[
 -8(\rho-1)(\rho+1)^5(\rho^2+1)^5
   (\rho^6+\cdots+\rho+1)
   (5\rho^6-4\rho^5+\rho^4+6\rho^3+\rho^2-4\rho+5).
\]
The reciprocal sextic is squarefree and disjoint from \(H=0\), \(\rho=0\), and \(\rho^7=1\). After excluding the diagonal and zero denominators, exactly twelve simple ratio roots remain. Each gives seven simple ordered pairs. Their simplicity also follows directly from solving \(z^7=-1/H(\rho)\) and the simple ratio equation. Finally,
\[
      G(S)^2-128^2S(1+S)^7=896S^3+2912S^2+3216S+1225 ,
\]
which is squarefree and nonzero at \(S=0,-1\). It has no root in common with \(G\).

## Distribution at a high position

**Theorem PlanarFortyTwoHighDoubleBranches.**
The ordered double algebra over a punctured neighbourhood of any high root has the following distribution of its \(198\) sheets. The last column specifies the projection indices wherever determined by the limiting equations.

| Limiting charts of the two momenta | Sheets | Projection indices over \(y\) |
|---|---:|---|
| First face with nonzero \(r\)-value | 6 | \(1^6\) |
| Both in the scalar chart at the quadruple root | 84 | \(7^{12}\) |
| One simple chart and the scalar chart in either order | 84 | \(7^{12}\) |
| Both in the same simple chart | 12 | Each cubic root multiplicity twice, for each chart |
| The two distinct simple charts | 12 | Six sheets in each order, described below |

No point over the high root in the normalization of the reduced double curve is fixed by exchange.

Put \(z=y^3 X\) in the first characteristic chart. The rescaled \(y^6P,y^{21}Q,y^{17}W,y^{51}Z\) are regular polynomials in \(z\) at \(y=0\), with leading \(z\)-coefficients units there. For \(Z\), use \(\operatorname{ord}_b B=1\) in its leading coefficient and \(Z\)-degree \(99\). If one of the two scaled momenta tends to infinity, then equality of \(P\) forces both to do so. Their limiting ratio is a twelfth root of unity. The simultaneous equality of \(Q,W,Z\), with degrees \(42,34,99\), forces ratio one. The leading coefficient \(12\) of the divided difference at ratio one excludes distinct pairs there. Thus both scaled momenta stay bounded.

If their common limiting \(r\)-value is nonzero, equality of leading \(P,Q\) gives equality of \(r\), and equality of leading \(W\) gives equality of \(A\). The derivatives \(r',A'\) cannot both vanish there by the first-face differential equation. Thus a coincident limit is impossible by the divided differences of \(P,W\). At the six distinct limits from the proposition, the equations for \(P,W\), rescaled by their respective powers of \(y\), lift simply. Equality of \(Q\) then follows from equality of its square and agreement of its nonzero leading sign. This gives the first row.

Suppose both momenta limit to the quadruple root. Use the [scalar chart](POLYNOMIAL_WEYL_ORDER_TWELVE_FORTY_TWO_PLANAR_FIBRES.md#the-complete-pole-profile-over-a-high-root) with \(y=w^7\), \(\tau^3=X\), \(\tau y=\omega+u\), and \(z_1=u/w^5\). The cube root in \(\tau y\) is chosen analytically at the specified nonzero root. If one scalar parameter is unbounded then so is the other by equality of \(P\). In that sector
\[
                   y^6P\sim c_Pu^8,\qquad
                   y^{17}W\sim c_Wu^{21}
\]
with both constants nonzero. These equivalences follow from the first characteristic orders and the scalar cap: the support of \(y^6P\) in \(y,u\) satisfies \(7i+5j\ge40\), and that of \(y^{17}W\) satisfies \(7i+5j\ge105\). For \(|u|/|y|^{5/7}\to\infty\), terms of smaller \(u\)-degree are negligible by these inequalities, and terms of higher degree vanish relative to the indicated leading power. Equality of \(P,W\) gives a limiting ratio with eighth and twenty-first powers one, hence ratio one. At ratio one the divided difference of \(y^6P\) in \(u\) has initial \(8c_P u^7\) and excludes a distinct pair. The same support estimates apply after taking this divided difference.

In the bounded scalar chart the leading pair is \(U,V\) of the proposition. Its derivatives never vanish simultaneously. Indeed the derivative of \(U\) is \(1+8S\), and at its zero the scalar equation makes \(G'\ne0\). Thus coincident limits are excluded. The \(84\) simple ordered solutions lift in \(w\). The seventh-root deck change acts freely on their nonzero parameters by \(z_1\mapsto\varepsilon^{-5}z_1\). Its orbits give twelve branches over \(y\), each of index seven. An odd-order deck group cannot exchange two distinct parameters within the same ordered-pair orbit. This also shows there is no fixed point here.

For one simple and one quadruple chart, write \(u_s=y^3 v\) in the simple chart. If the scalar parameter on the other side were unbounded, equality of \(P\) would force \(|u_s|\asymp |u|^4\) and \(v\) unbounded. On the simple chart \(y^{17}W\) is regular in \(y,u_s\), with supports \(i+3j\ge17\), since \(W\) is regular in \(y,v\). In the unbounded scalar sector \(|y|=o(|u_s|^{7/20})\). Terms of \(u_s\)-degree at least six are \(o(|u_s|^{21/4})\). For \(j\le5\) their orders in this estimate are at least \(j+7(17-3j)/20\ge114/20>21/4\). These bounds also hold for the analytic series by taking the tail of degree at least six in \(u_s\). Hence \(y^{17}W=o(|u|^{21})\) on the simple side, incompatible with \(W\) on the other side. Here the estimate \(|u_s|\asymp |u|^4\) follows from \(y^6P=u_s^2+O(y^6)\) and \(|u|^8\gg |y|^6\).

In the bounded scalar chart put \(v=w^{-1}V_s\) in the simple chart. The normalized simple outputs tend to \(w^2P=V_s^2\) and \(w^7Q=V_s^7\). For the latter, the regular expansion of \(y^{21}Q\) in \(y,u_s\) has \(i+3j\ge21\); substituting \(y=w^7,\ u_s=w^{20}V_s\) reaches weight \(140\) only at \(i=0,j=7\). Since the scalar outputs have no common zero, \(V_s\ne0\) at a matching value. On the scalar face, having the ratio of the seventh power of \(w^2P\) to the square of \(w^7Q\) equal to one is the three-root equation in the proposition. The normalization of that equation follows from the limit of this ratio at scalar infinity, which is one by the inherited first-face powers. Each root gives seven scalar parameters, and the common \(V_s\) is unique, being \(w^7Q/(w^2P)^3\) at the limit. These \(21\) solutions lift simply because the three roots are simple and both output values are nonzero. They descend to three index-seven branches per choice of simple root and ordering. This gives \(84\) sheets, as in the table.

For a pair limiting to the same simple root, use \(P=v^2+C(y)\). Distinct momenta give opposite parameters and require \(Q(y,v)=Q(y,-v)\). If \(|v|\) approached infinity while \(u_s=y^3 v\to0\), the odd difference of \(y^{21}Q\) would be led by \(2u_s^7\). Indeed terms of \(u_s\)-degree below seven have \(i+3j\ge21\), and higher powers are negligible. Thus \(v\) stays finite. On the limiting chart write \(Q_0(v)=vH_d(v^2)\), where the cubic \(H_d\) has nonzero constant by the [quadratic chart proof](POLYNOMIAL_WEYL_ORDER_TWELVE_FORTY_TWO_PLANAR_FIBRES.md#the-simple-characteristic-roots). Its roots give six ordered sheets. At a multiple root, along constant \(P\) the first \(y\)-derivatives of \(Q\) at the opposite normalized parameters differ by a nonzero multiple of \(1/(v r_b'(\rho))\). The additional derivatives through the moving \(P\)-level vanish here because \(Q'_0(v)=0\). The local double equation is therefore smooth, and its projection index equals the multiplicity of that root of \(H_d\). Simple roots lift directly. Exchange swaps the limiting opposite parameters. Applying this to both simple roots gives the fourth row.

All first-face limits have now been accounted for except pairs between the two different simple roots. The full algebra has length \(198\), so those pairs have length \(12\), equally divided by exchange between the two orders. Their first-face limits are distinct, which rules out exchange-fixed points also in this last sector.

The difference norm also checks the distribution. In row order the orders in the base parameter of the momentum difference are \(-3,-16/7,-3,0,-3\). For the second row, the scalar coordinates differ at their limits, so their momentum difference has order \(w^{-16}\) with \(y=w^7\). Summing with sheet multiplicity gives \(-498\), the previously computed high-root order of the norm.

## Ends between the two simple charts

**Corollary PlanarFortyTwoDistinctSimpleDoubleEnds.**
At a high root let \(d_0,d_1\) be the two critical \(P\)-values of the simple charts, and let \(F_i(p)\) be the polynomial defined there by \(Q_0^2=F_i(P_0)\). If \(d_0\ne d_1\), put \(j=\deg(F_1-F_0)\), so \(0\le j\le2\). In each order of the two charts there is exactly one end with unbounded quadratic momentum parameters, of projection index \(6-2j\). The other \(2j\) sheets have both quadratic parameters finite. If \(d_0=d_1\), all six sheets in each order have finite quadratic parameters and come from the opposite nonzero roots of \(Q_0(v)\); their projection indices are the root multiplicities.

From the [quadratic polynomial part](POLYNOMIAL_WEYL_ORDER_TWELVE_FORTY_TWO_PLANAR_FIBRES.md#the-simple-characteristic-roots), evaluated at \(v=(p-d_i)^{1/2}\) with the positive formal sign, write
\[
           Q_{0,i}(v)=\sqrt p\, A(p)-c/\sqrt p+O(p^{-3/2}).
\]
Here \(c\) is common because both \(d_i\) solve the same critical quartic. Squaring shows \(\deg(F_1-F_0)\le2\). Each \(F_i\) has the unique odd-multiplicity zero \(d_i\), so this difference is nonzero when \(d_0\ne d_1\).

At an unbounded pair write \(L=p^{-1/2}\), choosing the sign so both normalized parameters are \(v_i=L^{-1}(1-C_i(y)L^2)^{1/2}\). Opposite leading signs cannot match \(Q\). The two simple-root charts give, at fixed \(p\),
\[
               \partial_y Q_i\big|_p =
                     \frac{\pm\kappa}{2v_i\,u_{z,i}} ,
\]
with the same sign convention, where \(u_{z,i}\) tends to \(r_b'(\rho_i)\) at the corresponding simple root. Thus the difference of the two \(Q\)-values is
\[
             D(y,L)=D(0,L)+yL\,\Psi(y,L,y^3/L),
                  \qquad \Psi(0,0,0)\ne0.
\]
The formula follows by integrating the derivative at fixed \(p\) from \(0\) to \(y\). The unit is analytic in the three displayed inputs by the analytic charts and the square root of \(1-C_i(y)L^2\). Its constant is nonzero by the [simple-root derivative inequality](POLYNOMIAL_WEYL_ORDER_TWELVE_FORTY_TWO_PLANAR_FIBRES.md#monodromy-of-the-generic-position-projection).

For distinct \(d_i\), \(D(0,L)\) has order \(7-2j\). The implicit equation gives one smooth branch \(y=L^{6-2j}\) times a unit. Substitution of \(y=L^{6-2j}Y\) makes the remaining equation analytic and linear with nonzero derivative in \(Y\) at \(L=0\). All unbounded pairs lie here, since the relation for \(D\) forces this scale. This proves the index and leaves \(2j\) sheets in each order.

For equal \(d_i\) the leading difference vanishes identically, so the analytic unit rules out unbounded pairs. At finite nonzero matching parameters of the same sign, the two fixed-level derivatives still differ, and the common \(Q_0\)-value makes the local difference \(y\) times a unit. At zero parameters use \(Q\) as the local level coordinate instead, since \(Q_0'(0)\ne0\). Then the first derivatives of \(P\) at fixed \(Q\) differ by the same inequality, also ruling out a branch. The only remaining limits have opposite nonzero parameters where \(Q_0(v)=0\). At multiple roots their fixed-\(P\) derivative difference is a nonzero multiple of \(1/r_b'(\rho_0)+1/r_b'(\rho_1)\). Indeed the derivative ratio is \(-(q-1)^4\), and the quadratic for \(q\) implies \((q-1)^4\ne1\). Thus the equations are smooth there with projection index given by the root multiplicity, as for a single simple chart. This accounts for all six sheets in each order.

## A fixed-point consequence at degree sixteen

**Corollary PlanarFortyTwoGlobalExchangeFixedPoints.**
For the degree-sixteen transfer, the normalized projective reduced ordered double curve has respectively \(6,6,8,6,12,6\) exchange-fixed points in the six [index-four cusp profiles](POLYNOMIAL_WEYL_ORDER_TWELVE_FORTY_TWO_INDEX_FOUR_CUSPS.md#all-singular-branch-profiles-at-degree-sixteen).

At finite positions away from the high roots, affine momentum pairs have no fixed points by the Keller identity. The only momentum-infinite ends there are the three roots of the residue numerator \(g\), with one fixed point each by the [nonhigh end calculation](POLYNOMIAL_WEYL_ORDER_TWELVE_FORTY_TWO_POSITION_DOUBLE_LOCUS.md#ends-at-the-roots-of-the-residue-numerator). The high roots contribute none by the theorem above. The outer diagonal contributes \(3,3,5,3,9,3\) by the [index-four double-branch calculation](POLYNOMIAL_WEYL_ORDER_TWELVE_FORTY_TWO_INDEX_FOUR_CUSPS.md#the-ordered-double-locus-at-these-cusps). Adding the finite nonhigh contribution gives the result.

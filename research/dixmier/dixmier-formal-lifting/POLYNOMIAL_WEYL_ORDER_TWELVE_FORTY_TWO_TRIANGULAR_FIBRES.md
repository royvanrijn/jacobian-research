# Triangular target directions of the order-\((12,42)\) transfer

The scalar high-root charts determine the poles in every triangular target direction. For \(Q+\phi(P)\) the poles change when \(\deg\phi\) reaches four and depend on divisibility by seven thereafter. For \(P+\psi(Q)\) they depend on the parity of \(\deg\psi\). These pole lists give the genera and puncture counts of the generic fibres. A resultant identity also gives the total Milnor number of the critical points on the outer boundary in these directions.

Use the [unsplit planar transfer](POLYNOMIAL_WEYL_ORDER_TWELVE_FORTY_TWO_TARGET_CYLINDER.md) over \(\mathbb C\). Write \(H\) for the even number of high roots, let \(f\) be their monic squarefree product, and put \(M=5H/2\). The geometric degree of the transferred pair \(P,Q\) is \(8H\). Its unique affine branch curve has birational outer parametrization \((P_0(v),Q_0(v))\) of degrees \(12,42\). A generic meridian is a single \((M-1)\)-cycle, by the [affine branch calculation](POLYNOMIAL_WEYL_ORDER_TWELVE_FORTY_TWO_RESIDUE_MONODROMY.md#lifting-the-two-residue-permutations).

## Poles on triangular target fibres

**Theorem PlanarFortyTwoTriangularFibrePoles.**
Choose a target coordinate \(T\) from the following table and let \(J\) be the displayed companion. Let \(d=\deg_X T\). On the completed generic \(T\)-fibre the two rightmost columns give, over each high position \(b\), the orders of the poles of \(J\) and the orders of \(x-b\) at those same points. These are all the poles of \(J\). Superscripts in the table indicate repetitions, and \(\phi=0\) is included in the degree-at-most-three row.

| \(T\) | \(J\) | \(d\) | Pole orders of \(J\) | Orders of \(x-b\) |
|---|---|---:|---|---|
| \(Q+\phi(P)\), \(\deg\phi\le3\) | \(P\) | \(42\) | \(2^4\) | \(7^4\) |
| \(Q+\phi(P)\), \(\deg\phi=k\ge4\), \(a=\gcd(k,7)\) | \(P\) | \(12k\) | \(1,(7/a)^a\) | \(k,(7k/a)^a\) |
| \(P+\psi(Q)\), \(\psi\) constant | \(Q\) | \(12\) | \(1,7\) | \(1,7\) |
| \(P+\psi(Q)\), \(\deg\psi=\ell\ge1\), \(c=\gcd(\ell,2)\) | \(Q\) | \(42\ell\) | \((2/c)^{4c}\) | \((7\ell/c)^{4c}\) |

Use the [scalar high-root chart](POLYNOMIAL_WEYL_ORDER_TWELVE_FORTY_TWO_PLANAR_FIBRES.md#the-complete-pole-profile-over-a-high-root). It has \(x-b=w^7\), \(X=\tau^3\), and \(\tau(x-b)=\omega+w^5 z\) for a nonzero \(\omega\). In these variables
\[
                  U=w^2P,\qquad V=w^7Q
\]
are regular at \(w=0\), with specializations
\[
       U(0,z)=A_1 z(A_0+\omega^5 z^7),\qquad
       V(0,z)=A_2 G(\omega^5 z^7).
\]
Here \(A_1,A_2,A_0\) are nonzero, \(G\) has four distinct roots, and \(G(0)G(-A_0)\ne0\). Thus \(U(0,z)\) has eight simple zeros and \(V(0,z)\) has twenty-eight simple zeros. Their zero sets are disjoint. The change of the seventh-root choice is
\[
           (w,z)\longmapsto(\varepsilon w,\varepsilon^{-5}z),
                         \qquad \varepsilon^7=1 .
\]
It acts freely for \(w\ne0\) and permutes the nonzero centres on \(w=0\) in packets of seven.

For \(T=Q+\phi(P)\) with \(\deg\phi\le3\), multiply \(T=T_*\) by \(w^7\). The degree bound makes the equation regular with specialization \(V(0,z)=0\). Its twenty-eight simple roots give four branches after dividing by the seventh-root choice. Each has a parameter \(w\) on one representative chart and has \(P\) of pole order two and \(x-b\) of order seven.

For the higher-degree rows, use the following Newton calculation. Suppose \(Y=0\) is a simple zero of the specialization of \(w^\alpha J\), and use \(Y=w^\alpha J\) as a coordinate with \(w\). Write the other component as \(w^{-\gamma}L(w,Y)\) with \(L(0,0)\ne0\). The level equation for \(w^{-\gamma}L(w,Y)+\theta(J)\), where \(\deg\theta=n\) and \(n\alpha>\gamma\), becomes
\[
 w^{n\alpha-\gamma}L(w,Y)+w^{n\alpha}\theta(w^{-\alpha}Y)-T_*w^{n\alpha}=0.
\]
Let \(g=\gcd(n,\gamma)\). Put \(w=s^{n/g}\) and \(Y=s^{(n\alpha-\gamma)/g}\eta\). After cancellation the initial equation is \(\theta_n\eta^n+L(0,0)=0\), where \(\theta_n\ne0\) is the leading coefficient. Every lower term of \(\theta\) has positive additional order. The \(n\) simple nonzero roots lift, and the \(n/g\)-th roots of unity permute them freely in packets because \(\gcd(n/g,(n\alpha-\gamma)/g)=1\). This gives \(g\) branches, each with order \(n/g\) for \(w\) and pole order \(\gamma/g\) for \(J\). Weierstrass preparation in \(Y\) gives degree \(n\), so these account for the whole local level equation.

For \(T=Q+\phi(P)\), \(k\ge4\), use \((Y,\alpha,\gamma,n)=(U,2,7,k)\) at the eight zeros of \(U(0,z)\). At the seven nonzero zeros the chart change is free on centres. One representative gives \(a\) branches of pole order \(7/a\), with \(x-b\) of order \(7k/a\). Over the fixed centre \(z=0\), the total pole order upstairs is \(a(7/a)=7\). The quotient by the seventh-root choice has degree seven on the punctured branches, so the total downstairs pole order is one. Accordingly there is a single simple pole. The upstairs order ratio of \(x-b\) to the pole order of \(P\) is \(k\), so the downstairs order of \(x-b\) is \(k\).

For \(T=P+\psi(Q)\), \(\ell\ge1\), use \((Y,\alpha,\gamma,n)=(V,7,2,\ell)\) at the twenty-eight zeros of \(V(0,z)\). The four free packets of centres each give \(c\) branches with pole order \(2/c\) and order \(7\ell/c\) for \(x-b\). The constant-\(\psi\) row is the existing [pole computation](POLYNOMIAL_WEYL_ORDER_TWELVE_FORTY_TWO_PLANAR_FIBRES.md#the-complete-pole-profile-over-a-high-root).

The chart at \(w\ne0\) parametrizes actual source points and the indicated roots give distinct branches modulo the seventh-root choice. In every row the pole orders sum to eight over each high position. Since \((T,J)\) is a polynomial target coordinate change of \((P,Q)\), the companion has degree \(8H\) on the generic \(T\)-fibre. The exhibited branches exhaust its poles. The combined pole list on the normalized components is the same at every finite \(T\)-level, by [the cylinder invariance of polar ends](../../jacobian/monodromy-theory/KELLER_DIRECTIONAL_RAMIFICATION_AND_POLAR_ENDS.md#meridians-and-polar-ends).

## Genus and punctures

**Corollary PlanarFortyTwoTriangularFibreTopology.**
Let \(r\) denote the number of poles per high position in the table, counting branches. For the completed generic \(T\)-fibre the genus \(g_T\) and the number \(n_T\) of punctures of its original affine fibre satisfy
\[
          2g_T-2=d(M-2)-(8+r)H,\qquad
          n_T=d+\bigl(d/6+r\bigr)H .
\]
The position projection \(x\) has degree \(d\). It has \(d\) unramified points above \(x=\infty\), and above each high position its only possible ramification is at the scalar poles in the table.

The [directional ramification formula](../../jacobian/monodromy-theory/KELLER_DIRECTIONAL_RAMIFICATION_AND_POLAR_ENDS.md) applies to the polynomial target coordinates \(T,J\). In each row the degree of \(T\) on the outer branch normalization equals \(d\), since the leading degrees \(12,42\) admit no tie under these triangular changes. Thus the finite ramification of the companion on the connected generic fibre has total \(d(M-2)\), and its poles number \(rH\). Riemann–Hurwitz gives the genus formula.

Besides the polar ends, the generic fibre has \(d\) omitted points on the outer boundary. At each high position the nonzero simple-root chart gives \(d/6\) omitted points. The degree on that chart follows from its [boundary restrictions](POLYNOMIAL_WEYL_ORDER_TWELVE_FORTY_TWO_PLANAR_FIBRES.md#the-simple-characteristic-roots), of respective degrees two and seven for \(P,Q\). The other simple-root chart at each high position is affine. These are all the generically visible omitted curves by the [generic \(P\)-fibre calculation](POLYNOMIAL_WEYL_ORDER_TWELVE_FORTY_TWO_OUTER_FIBRES.md#exact-geometric-degree-and-genus). No omitted curve is vertical for either target coordinate, by the same directional result. This counts the punctures.

The position projection has degree \(d\) from the \(X\)-degree and geometric irreducibility of the generic \(T\)-fibre. On the two regular charts over any high position the \(d/6\) simple roots per chart lift with \(x-b\) as parameter. The remaining orders from the table sum to \(2d/3\). On the outer chart the \(d\) simple roots of the restricted \(T\)-level lift with \(1/x\) as parameter, giving the asserted projection data.

For \(H=2\), \(Q\) has genus \(52\) and \(64\) punctures on its generic fibre. More generally all \(Q+\phi(P)\) with \(\deg\phi\le3\) have these same numbers when \(H=2\).

## Polar norms and outer critical points

**Corollary PlanarFortyTwoTriangularPolarNorm.**
For each direction \(T\) of the table, take the norm over \(\mathbb C(x)\) on the polar algebra \(T_X=0\), counting components with their multiplicities. For some \(C\ne0\),
\[
               \operatorname{Norm}(T_x)=C f^{7-5d/2}.
\]
In the outer chart with \(r_0=1/x\), the total intersection multiplicity of \(T_v,T_{r_0}\) along \(r_0=0\) is
\[
                          \left(\tfrac52 d-7\right)H-2(d-1).
\]
This is the sum of the Milnor numbers of the critical points of \(T\) on that boundary.

Use the [polar resultant](POLYNOMIAL_WEYL_ORDER_TWELVE_FORTY_TWO_TARGET_CYLINDER.md#the-polar-norm) for \(P_X,Q_X\), which is a nonzero constant times \(f^{1253}\). Their leading \(X\)-coefficients are scalar multiples of \(f^{30},f^{105}\), and that of \(T_X\) is a scalar multiple of \(f^{5d/2}\). For \(T=Q+\phi(P)\) the derivative \(T_X\) agrees with \(Q_X\) on \(P_X=0\), so the resultant with first member \(P_X\) changes by the leading coefficient of \(P_X\) to power \(d-42\). The resulting exponent is \(30d-7\). For \(T=P+\psi(Q)\) the analogous exponent with first member \(Q_X\) is \(105d-7\), using the power \(d-12\). In both cases the norm of the companion derivative \(J_X\) on \(T_X=0\) has exponent \(5d/2-7\). There \(T_x J_X\) is the nonzero scalar Jacobian up to sign. This gives the norm identity for \(T_x\).

On the outer chart \(X=\Lambda(r_0)+r_0^M v\), the polar equation is \(T_v=0\). The leading coefficient of \(T_v\) is a unit at \(r_0=0\), so its normalized components, counted with multiplicity, have \(d-1\) sheets there. Along the polar one has \(T_x=-r_0^2 T_{r_0}\), since the term with \(T_v\) vanishes. The norm identity gives total order \((5d/2-7)H\) for \(T_x\) there. Subtract \(2(d-1)\) to obtain the intersection multiplicity of the two chart derivatives, using the norm formula for a polar with its component multiplicities. The two derivatives have no common curve through this boundary, so their local intersection lengths are the Milnor numbers.

## Local sheet counts on the outer curve

**Corollary PlanarFortyTwoOuterLocalDegrees.**
At an outer chart point \((0,v_0)\) let \(N_{v_0}\) be the local degree of the finite map germ \((P,Q)\), and put \(e=M-1\). Then
\[
                 N_{v_0}=e+I_{(0,v_0)}(P_v,Q_v).
\]
In particular the excess above \(e\) is supported exactly at the common zeros of \(P'_0,Q'_0\), and it sums to \(9H/2\). If \(s_{v_0}\) is the minimum of their two vanishing orders, then
\[
 s_{v_0}\le N_{v_0}-e\le (M-2)s_{v_0},\qquad
 \deg\gcd(P'_0,Q'_0)\ \ge\ \left\lceil\frac{9H}{5H-4}\right\rceil.
\]
At the smallest \(H=2\), the local degree excess sums to nine and the common derivative has degree at least three.

Use \(r_0,v\) as local source coordinates. The map germ is finite since its fibre is discrete off the boundary and the boundary map is nonconstant. Put \(F=P-P(0,v_0)\), \(G=Q-Q(0,v_0)\), and \(D=F_{r_0}G_v-F_vG_{r_0}\), a nonzero constant times \(r_0^{e-1}\). The germs \(F,F_v\) have no common component. A component of \(F=0\) off \(r_0=0\) is reduced by the Keller identity; if \(F_v\) vanished along it, differentiation in its normalization would make \(r_0\) constant. The restriction \(F(0,v)\) is nonzero, so \(r_0=0\) itself is not a component.

On each branch of \(F=0\), the identity \(F_v\,dG=-D\,dr_0\) gives
\[
                   I(F,D)=I(F,F_v)+I(F,G)-I(F,r_0).
\]
The polar equality \(I(F,F_v)=I(F_{r_0},F_v)+I(F,r_0)-1\) follows by differentiation on the branches of \(F_v=0\) and \(I(r_0,F_v)=I(r_0,F)-1\). This is the transverse-coordinate form of [Teissier's lemma](https://arxiv.org/pdf/1207.1600). Set \(m=I(F,r_0)\) and \(\mu=I(F_{r_0},F_v)\). The two equalities give \(N_{v_0}=I(F,G)=(e-1)m+1-\mu\).

On \(F_v=0\) the determinant identity becomes \(F_{r_0}G_v=D\), so
\[
                   \mu+I(F_v,G_v)=(e-1)(m-1).
\]
Substitution proves the local degree formula. The polar identity in the [outer chart](POLYNOMIAL_WEYL_ORDER_TWELVE_FORTY_TWO_TARGET_CYLINDER.md#the-polar-norm) gives the sum. Finally the local algebra of \((P_v,Q_v)\) has quotient of length \(s_{v_0}\) modulo \(r_0\). The determinant puts \(r_0^{e-1}\) in the ideal, and multiplication by powers of \(r_0\) makes each successive quotient of its \(r_0\)-adic filtration a quotient of the first. This proves both length bounds.

The [outer cusp restrictions](POLYNOMIAL_WEYL_ORDER_TWELVE_FORTY_TWO_OUTER_CUSPS.md)
combine the local excess identity with the rank-one and multiplicity-two
normal forms. They exclude outer charts of differential rank one
everywhere except at the two stated high-root counts.

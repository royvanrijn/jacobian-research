# Averaging the degree-ten pair components for Erdős #364

Let \(C(X)\) count centers \(2\le n\le X\) of three consecutive positive powerful integers. Using the fixed-shift component and lattice estimates of [Plymesser's powerful-products preprint](sources/plymesser-powerful-products-hexagon2611.00602-v1.0.1.pdf) and the published affine curve bounds of [Castryck–Cluckers–Dittmann–Nguyen](sources/castryck-cluckers-dittmann-nguyen-dimension-growth-2020.pdf),

\[
C(X)\ll_\varepsilon X^{\vartheta+\varepsilon},\qquad \vartheta=\frac6{25}-10^{-7}.
\]

The degree-ten components are counted by averaging over pair intervals. A large leading form gives a direct gain in the affine curve count. When both of two leading forms are small, a boundary point determines a short list of reduced lattice directions, which bounds the number of such intervals. A varying rectangle in the contact determinant keeps the few poorly balanced pair windows within this budget. The argument counts all pair points on the selected components, so it applies in particular to those that extend to fundamental triples.

## Pair selection and the contact width

Use the notation and [stability estimate for the old pairwise bound](ERDOS_364_PELL_AND_JOINT_DETERMINANT.md#the-only-maximizer-of-the-pairwise-bound). Work where the old envelope is greater than \(\vartheta\). The tolerance \(77/320-\vartheta=1/1600+10^{-7}\) is less than \(1/1500\). All three pair sums are within \(\kappa=(8/9)\sqrt{1/1600+10^{-7}}\) of \(1/3\), since the middle sum lies between the two adjacent sums in sorted order. Put the three exponents in increasing order, write \(s_* = 103/305\), and select whichever adjacent pair in this ordering has sum \(s\) farther from \(s_*\). The shifts in the original integers are still 1 or 2. If \(\Delta\) is the difference of the selected kernel exponents and \(\tau=s_*-s\), then \(\Delta\le 2|\tau|\): the difference of the two adjacent sums is the span of the three exponents. The bounds \(|\tau|\le\kappa+4/915<1/35\) and \(\Delta\le2\kappa<1/20\) have fixed slack.

Write the selected fixed-shift equation as \(e^3v^2-d^3u^2=h\) with \(h\in\{1,2\}\). The size parameters \(D,E,U,V\) satisfy \(D^3U^2=E^3V^2\), the common pair-window scale. This is comparable to the dyadic height \(X\), with \(O(1/\log X)\) changes in the window exponents. Set

\[
m_0=\frac{57}{305},\qquad \xi=10^{-8},\qquad
m=m_0-\min\{\tau^2/2,1/5000\}+\xi.
\]

There is a bounded-degree pair auxiliary equation in each normalized kernel-ratio interval of width \(X^{-m}\). To prove this, use contact polynomials as in [the quadratic companion, §§5.2–5.3](sources/plymesser-powerful-quadratics-hexagon2611.00699-v1.0.1.pdf), but vary the degree rectangle with \(s\). Let \(b=3(1+\tau)\) and \(a=2(1-\tau)\), and take the integer rectangle of bidegree \((\lfloor bL\rfloor,\lfloor aL\rfloor)\) for a sufficiently large fixed \(L\). The first slope locally has weight \(m\), and the error from its smooth graph has weight 1. The graph in normalized coordinates is \(z=r^{3/2}+O_h(X^{-1})\). Each contact polynomial of bidegree \((2,1)\) or \((1,2)\) can be chosen with weighted order at least \(5m\), by imposing five homogeneous conditions on the six coefficients to cancel the pure first-variable jets. The weight-1 terms cause no loss since \(5m<1\).

Here is the resulting weight calculation for the variable rectangle. For \(a\le b\le2a\), spend \(rL\) factors on the two contacts. The maximal residual area in units of \(L^2\) is

\[
\mathcal A(r)=\begin{cases}
(b-2r)(a-r),&0\le r\le b-a,\\
(b+a-3r)^2/4,&b-a\le r\le(b+a)/3.
\end{cases}
\]

Indeed the residual rectangle for \(kL\) factors of bidegree \((2,1)\) and \((r-k)L\) of bidegree \((1,2)\) has sides \(b-r-k,a-2r+k\). These sides become equal at the maximizing choice \(k=(b-a+r)/2\) when it is feasible, and otherwise the maximizing choice is \(k=r\). Taking integer parts incurs an \(O(L)\) error in the area. The number of jets below weight \(zL\) is \(z^2 L^2/(2m)+O(L)\) at bounded nonnegative \(z\). Multiplication by the chosen contact factors therefore gives a lower bound of \(L^2(\mathcal A(r)-(z-5mr)_+^2/(2m))+O(L)\) on the dimension of the order-\(zL\) subspace. Integrating a nonnegative lower envelope of these dimensions gives a lower bound on the weight of any full rank Taylor coefficient minor, by counting how many of its independent functionals can have weight below each level.

For the integral, write \(j=b-a\) and \(k=2a-b\). The maximizer leaves \(r=0\) at \(z_1=(b+2a)/5\), reaches \(r=j\) at \(z_2=5mj+3k/5\), and reaches \(r=(a+b)/3\) at \(z_3=5m(a+b)/3\). This follows by differentiating the two quadratic expressions for \(\mathcal A(r)-(z-5mr)^2/(2m)\), with \(r\le z/(5m)\). The envelope is \(ab-z^2/(2m)\) up to \(z_1\), adds \((5z-b-2a)^2/(2(25m-4))\) up to \(z_2\), and equals \(9(z_3-z)^2/(2m(50m-9))\) up to \(z_3\). Its integral is

\[
I=\left(\frac{10m}{3}+\frac8{15}\right)j^3
 +\left(\frac{15m}{2}+\frac65\right)j^2k
 +\left(5m+\frac9{10}\right)jk^2
 +\left(\frac{10m}{9}+\frac15\right)k^3.
\]

The weight lower bound is \(I L^3+O(L^2)\). Clearing the rows and undoing the column normalization in the full rectangle costs \((D E)^{\lfloor bL\rfloor/2}(U V)^{\lfloor aL\rfloor/2}\) per row up to a fixed factor. Thus the leading cost of the determinant is \(ab(sb+(1-3s/2)a)L^3/2\) in the exponent of \(X\). The window convention changes coordinate exponents by only \(O(1/\log X)\). The surplus in the leading coefficient is

\[
\frac{305}{18}(m-m_0) + 18\tau^2 +\frac{25}{6}(m-m_0)\tau
 +\left(\frac52-\frac{175m}{6}\right)\tau^2
 +\left(\frac{77}{30}+18s-\frac{935m}{18}\right)\tau^3.
\]

This is uniformly positive. For example, use \(m-m_0\ge-\tau^2/2+\xi\), \(|m-m_0|\le\tau^2/2+\xi\), \(5/2-175m/6>-3\), and \(|77/30+18s-935m/18|<2\) on the indicated intervals. The lower bound has positive coefficients on both \(\tau^2\) and \(\xi\) after using \(|\tau|\le1/35\). Taking \(L\) fixed large enough for the strict surplus, the same Taylor truncation and integer determinant argument as in the preprint gives the auxiliary equation. Cofactors of the bounded-size evaluation matrix have polynomial height in \(X\). Consequently its bounded number of bounded-degree components can be handled by the original pair component estimates with this value of \(m\).

## The two plane projections of a pair component

Let \(\Gamma\) be a geometrically integral \((2,2)\) component, and use the two reduced integral coordinate blocks \(\lambda,\nu\) of the pair lattice. Write \(F\) for the section \(e^3v^2-d^3u^2\), in reduced coordinates. Its boundary on \(\Gamma\) is a divisor of degree 10 since \(\Gamma\) contains points with \(F=h\ne0\). Choose bounded integer linear forms \(\ell_b(\lambda),\ell_a(\nu)\) and for the first plane projection take the three sections

\[
F/h,\qquad \ell_b^2\ell_a^2\lambda_0,\qquad \ell_b^2\ell_a^2\lambda_1 .
\]

For the second plane projection take \(F/h\) and \(\ell_b^3\ell_a\nu_i\) for \(i=0,1\). The linear forms can be chosen separately in the projections to avoid all boundary zeros and give birational maps. For the first projection the field of ratios contains the degree-two kernel-block projection field. Its deck involution fixes the \(\lambda\) ratios and does not fix the \(\nu\) ratio. For \(\ell_a=\nu_0+j\nu_1\), invariance of the resulting ratio to \(F/h\) under the involution gives a polynomial identity in \(j\) of degree at most two. This polynomial is nonzero, since invariance of all its coefficients would imply invariance of \(\nu_1/\nu_0\) by taking a ratio of the linear and constant coefficients. At most two integers are excluded this way and at most ten by the boundary. For the second projection the analogous test uses the cubic \(\ell_b^3\) and excludes at most three integers beyond the boundary exclusions. Choosing from bounded lists works uniformly, including for singular geometrically integral components.

The maps are basepoint-free on the curve and birational, so their plane images have degree 10. The fiber sizes are bounded by 10, by pulling back a line through an image point that does not contain the image. At genuine points \(F=h\) and the affine plane coordinates are integral of height \(O(B)\), where \(B=\max(1,L_1^3T_1^2)\) in the notation of the pair lattice estimates. In each projection the degree-ten part of the primitive plane equation vanishes on the corresponding reduced block ratios of every boundary point, since the ratio of the two nonzero-section coordinates at infinity is precisely that block ratio. This uses the actual affine boundary, without a Weierstrass conversion.

For a primitive degree-ten plane equation \(f\) with highest form of coefficient norm \(H=\|f_{10}\|\), [Castryck–Cluckers–Dittmann–Nguyen, Corollary 4.2.4](sources/castryck-cluckers-dittmann-nguyen-dimension-growth-2020.pdf) bounds the affine integral points by

\[
\ll B^{1/10}H^{-1/100}(1+\log H+\log(2B)) + \log(2B).
\]

Put \(\gamma=75m-14+10^{-5}>0\). If either plane projection has \(H\ge X^\gamma\), the corresponding component contributes \(\ll_\varepsilon X^\varepsilon(1+B^{1/10}X^{-\gamma/100})\). The logarithm in \(H\) is absorbed using the decreasing function \((\log H)H^{-1/100}\) above a fixed threshold; no upper bound on the full plane equation is needed. The pair exponents differ by less than \(1/20<2m/3\), so the usual reference heights satisfy \(L_*^3T_*^2\asymp X^{1-5m/2}\). Summing the lattice moments gives exponent

\[
m+\frac{1-5m/2}{10}-\frac{\gamma}{100}=\vartheta.
\]

## Counting intervals with small leading forms

Let \(\eta=1/8000\). Intervals where \(L_1>X^\eta L_*\) or \(T_1>X^\eta T_*\) are affordable by the ordinary affine bound. The two lattice tails in [Plymesser, Lemma 6.1](sources/plymesser-powerful-products-hexagon2611.00602-v1.0.1.pdf) give \(O(M2^{-2j})\) intervals where either scaled first-minimum reciprocal exceeds \(2^j\), for \(M\asymp X^m\). In bins with indices \(i,j\) for the two reciprocals, the count is \(O(M2^{-2\max(i,j)})\) and the weight from \(B^{1/10}\) is \(O(L_*^{3/10}T_*^{1/5}2^{3i/10+j/5})\), with negative indices capped at zero. Summing the tail saves \(X^{-(3/2)\eta+\varepsilon}\). Its exponent is at most \(1/10+3(m_0+\xi)/4-(3/2)\eta<0.239976442<\vartheta\). In all remaining intervals the ordinary per-component bound is \(\ll_\varepsilon X^{(1-5m/2)/10+\eta/2+\varepsilon}\).

For a degree-ten component with both leading forms of height less than \(X^\gamma\), each of the two reduced block ratios of every boundary point lies in a set of \(O(X^{11\gamma})\) directions, obtained from all nonzero integer binary forms of degree 10 with coefficients bounded by \(X^\gamma\). There is a boundary point far from the interval center in at least one original normalized block. In the ratios \(x=d/e, y=v/u\), the boundary is parametrized by \(x=t^2,y=t^3\). The pullback of a nonzero \((2,2)\) equation has degree at most 10, with no \(t^9\) term. It is nonzero since the boundary curve has bidegree \((3,2)\). If the degree drops, \(t=\infty\) is a boundary point. Otherwise the sum of the roots is zero, so a root has nonpositive real part. On rescaling \(x,y\) by \(D/E,(D/E)^{3/2}\), the same statement holds for the rescaled \(t\). A normalized boundary direction near both positive centers \(\rho,\rho^{3/2}\) would have \(t=y/x\) near \(\sqrt\rho>0\) in those rescaled coordinates. Thus the far direction has chordal distance bounded below by an absolute positive constant in at least one of the blocks, uniformly in \(\rho\) in the compact ratio chart.

Here is the interval count for a fixed direction from the small-form list. In the kernel block its lattice is obtained by scaling \((d,e)\) to \((M(d/D-\rho e/E),e/E)\), up to bounded constants. The square-part case uses \((V,U)\) and center \(\rho^{3/2}\). Denote the two physical scales for the block by \(S\ge S'\). Let \(g_1,g_2\) be the reduced lattice basis and \(P\) the corresponding unimodular basis in the original integer coordinates. The reduced angle is bounded away from zero, the covolume is comparable to \(M/(SS')\), and on the retained intervals

\[
\|g_2\|\ll X^\eta\sqrt{M/(SS')},\qquad
\|g_2\|/\|g_1\|\ll X^{2\eta}.
\]

The two normalized original rows as linear forms in reduced coordinates differ, after multiplying one by the center ratio, by a form of norm \(O(\|g_2\|/M)\). At a far complex direction this difference has size at least a constant times the maximum of their absolute values, by the chordal distance formula. Hence the reduced direction is within \(O(X^{-m+2\eta})\) of the rational kernel direction of the row of \(P\) corresponding to the larger physical scale \(S\). To see the normalization in this distance estimate, the normalized row has norm \(\gg\|g_1\|\) by the covolume and the bound on the norm of the other scaled row. This also holds for the numerator row by reversing the normalized slope, whose center is bounded above and away from zero. The primitive rational direction has height \(O(X^\eta\sqrt{MS/S'})\). Farey separation gives \(O(1+X^{4\eta}S/S')\) possibilities for it in the chordal neighborhood. The same bound applies around a nonreal direction: if the neighborhood contains rational directions, their pairwise separation inside it gives the usual real packing estimate.

A fixed such rational direction determines the corresponding integer row up to sign. It allows only boundedly many intervals. To prove this, divide each reduced basis column by its length. The scaled matrix and its inverse are bounded, so some column has its normalized larger original entry of magnitude \(\gg\|g_j\|\), where \(g_j\) is that column. At that index the ratio of the smaller raw entry to the larger is \(O(S'/S)\). All completions of the fixed primitive row to determinant \(\pm1\) differ by integral multiples of the fixed row after fixing the determinant sign. For either choice of index, the bound on the ratio allows only \(O(1)\) multiples. Finally, the matched column's normalized ratio is within \(O(1/M)\) of the interval center, or of its reciprocal if the numerator row was fixed. For the square-part block the center is the \(3/2\) power of the kernel center. Each of these center maps preserves grid separation up to fixed factors. Thus each full basis allows \(O(1)\) intervals.

In the kernel block \(S/S'\ll X^\Delta\), and in the square-part block \(S/S'\ll X^{3\Delta/2}\). Consequently at most \(O(X^{11\gamma+3|\tau|+4\eta})\) retained intervals can contain a degree-ten component with both leading forms small. Each interval has a bounded number of components. The saving of their summed bound below \(6/25\), before the arbitrary \(\varepsilon\), is at least

\[
\frac{39}{6100}+\frac{3299}{4}\min\{\tau^2/2,1/5000\}
   -3|\tau|-\frac92\eta-\frac{3299}{4}\xi-11\cdot10^{-5}.
\]

In the uncapped case the minimum of the quadratic part after the \(39/6100\) term is \(-18/3299\). In the capped case use \(|\tau|\le1/35\). The resulting upper exponent is less than \(0.239743504<\vartheta\) in both cases. This completes the degree-ten count. Curves with both bidegrees at least two other than \((2,2)\) have lifted degree at least 12 and summed exponent at most \(m+(1-5m/2)/12<0.231284162<\vartheta\) by the ordinary affine bound. A \(\mathbb Q\)-irreducible but geometrically reducible curve contributes only boundedly many rational intersections with its conjugate components. The fixed-ratio and bidegree-one cases are covered by the fixed-shift component estimates, and the exterior boxes by the old pairwise estimates. Summing the dyadic boxes yields the stated bound.

## Scope of the elliptic estimate

The proof closes the averaged counting target. A uniform degree-ten exponent \(162/1625\) on each surviving smooth \((2,2)\) component would also suffice at \(m\downarrow57/305\), since \(57/305+(65/122)(162/1625)=6/25\). That pointwise estimate is a separate open question, as is the joint-cell collision bound.

For example, the smooth pair curve \((x^2+2)y^2-28x+2=0\) contains \((x,y)=(2,3)\) with primitive shift-1 coordinates. The projection to \(x\) has four distinct simple branch points, the two zeros and two poles of \((28x-2)/(x^2+2)\) on \(\mathbb P^1\). Its affine boundary has no rational point: on \((x,y)=(t^2,t^3)\) it is cut by the Eisenstein polynomial \(t^{10}+2t^6-28t^2+2\), with no additional point at infinity. A birational map over \(\mathbb Q\) to a Weierstrass model has a rational pole above the origin, so cannot convert this entire affine chart regularly to the affine Weierstrass curve. This is a pair example. The averaged bound above keeps the boundary in the coefficient and lattice counts.

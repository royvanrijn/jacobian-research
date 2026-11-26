# Uniform small-coefficient counts and Erdős #137

Let \(\mathcal Q\) be the positive powerful integers and put \(\vartheta=6/25-10^{-7}\). For each \(\varepsilon>0\) there is \(\eta>0\) such that, uniformly for positive integers \(f_1,f_2,f_3\le Y^\eta\) and \(Y\ge2\),

\[
\#\{0\le N\le Y:N+r\in f_r\mathcal Q\ (1\le r\le3)\}
 \ll_\varepsilon Y^{\vartheta+\varepsilon}.
\]

This deduction uses the primitive window, fixed-shift component and lattice bounds in [Plymesser, *Counting powerful products of consecutive integers*, §§4–6 and Lemmas 2.3 and 8.2](sources/plymesser-powerful-products-hexagon2611.00602-v1.0.1.pdf), and the degree-ten affine bounds cited in the [averaged pair proof](ERDOS_364_AVERAGED_PAIR_COMPONENT_COUNT.md). The all-width corollary for Erdős #137 also uses [Tao, Lemma 3.1](https://arxiv.org/html/2603.27990v3#S3), the subpower upper bound on the width of a powerful interval product. All these source estimates are explicit inputs.

## Uniformity of the averaged window estimate

Use a primitive pair window at scale \(x\), with kernel exponents \(\beta_1',\beta_2'\), sum \(s\), and imbalance \(\Delta=|\beta_1'-\beta_2'|\). Let \(\tau=103/305-s\). Fix a sufficiently small absolute \(\sigma>0\). The hypotheses for this estimate are

\[
1\le |h'|\le x^\sigma,\qquad
|\tau|\le1/35,\qquad \Delta\le1/20,\qquad \Delta\le2|\tau|+\sigma.
\]

Under these hypotheses the number of primitive points in the window is \(O_\varepsilon(x^{\vartheta+\varepsilon})\), with constants independent of \(h'\) and the window. Use the same \(m=m_0-\min\{\tau^2/2,1/5000\}+10^{-8}\) and degree rectangle as in the averaged pair proof. The normalized graph error is now \(O(x^{-1+\sigma})\), by expanding the square root with relative shift \(h'/(d'^3u'^2)\) of that size. Every Taylor weight \(mi+(1-\sigma)j\) is at least \((1-\sigma)(mi+j)\). Equivalently, after truncation, the rank argument that bounds the weights of the nonzero Taylor minors loses at most the fraction \(\sigma\) in its weight lower bound. The old leading surplus is at least \(16\cdot10^{-8}\) uniformly on the compact range, and the leading Taylor integral there is bounded. Choose \(\sigma\) small enough to leave a fixed positive surplus. The Taylor remainders also have uniform constants because the zero-shift graph is the same smooth function on the compact normalized chart, and the new deviations have a uniform size bound. The integer auxiliary degree is then fixed, after this choice of \(\sigma\) and before \(h'\) grows. Its coefficients have polynomial height in \(x\) uniformly by the same cofactor bound.

The physical strip containment uses \(3|h'|M\le x\), with \(M\asymp x^m\), which holds for large \(x\) at this choice of \(\sigma\). The reduced lattices, their reference heights, and their moment bounds use only these strips and their cell grid. Lemma 8.2 and its proof give the uniform low-bidegree and fixed-ratio cases of Proposition 5.6 on primitive points, also for \(|h'|\) growing in this range. In the two plane maps on a degree-ten component, use first coordinate \(F_{\rm equation}/h'\). It equals 1 on the integer data, so the other two bounded-coefficient sections give integral affine coordinates of height \(O(L_1^3T_1^2)\). The rational coefficients of the first section do not enter the heights of the two counted coordinates. The boundary section is still \(F_{\rm equation}=0\), and in the original pair slopes its parametrization and missing ninth coefficient are unchanged. The fixed-height atlas for the two leading binary forms, the chordal escape, and the Farey and unimodular complement counts therefore have the same uniform constants.

Both balanced reference-height branches hold because \(\Delta\le1/20<2m/3\). Use \(\gamma=75m-14+10^{-5}\) and lattice-tail cutoff exponent \(1/8000\) as in the fixed-shift argument. The large-form contribution is unchanged at exponent \(\vartheta\), and the same tail and degree-at-least-12 bounds apply. The exceptional interval count uses \(3\Delta/2\le3|\tau|+3\sigma/2\), and so loses at most \(3\sigma/2\) in its exponent. The previously explicit upper exponent \(0.239743504\) for this branch has fixed slack below \(\vartheta\). Decrease \(\sigma\) if necessary. All the constants in the averaged window estimate are now uniform over this set of primitive windows and shifts.

## Repaired consecutive terms

Ignore \(N\le Y^\vartheta\), at a cost of at most \(1+Y^\vartheta\). Write \(N+r=f_r d_r^3u_r^2\), with \(d_r\) squarefree and \(F=\max_r f_r\). Use the normalization of [Plymesser, Lemma 8.3](sources/plymesser-powerful-products-hexagon2611.00602-v1.0.1.pdf). For \(i<j\), let \(d=f_i d_i\), \(e=f_jd_j\), \(u=f_j u_i\), \(v=f_iu_j\), and \(h=f_i^2 f_j^2(j-i)\). Divide the two blocks by \(g=\gcd(d,e)\) and \(q=\gcd(u,v)\). The resulting shift is the positive integer \(h'=h/(g^3q^2)\), and \(g,q,h'\le2F^4\). Both reduced products have size between \((N+1)/2\) and \(F^4(N+3)\) because \(g^3q^2\le h\). For \(\eta\) small enough, Lemma 2.3 applies and gives primitive dyadic windows of scale \(x\) with

\[
N+1\ll x\ll F^4(N+3).
\]

Let the raw exponents be \(\beta_r=\log d_r/\log(N+1)\) clipped to \([0,1/3]\). Each actual primitive window exponent for a selected pair differs from the corresponding raw exponent by \(O((1+\log F)/\log(N+1))\). Indeed the changes of kernel height from multiplication by the coefficients and division by \(g\) have logarithmic size \(O(1+\log F)\), and the same bound controls the change of the height scale and clipping. Since \(N>Y^\vartheta\), choose \(\eta\) so small that these discrepancies are within all fixed margins in the averaged window estimate. In particular a pair chosen from the raw exponents by maximizing the distance of the two sorted adjacent pair sums from \(103/305\) satisfies \(\Delta\le2|\tau|+\sigma\) at its actual window parameters. In that window use the actual value of \(s\) to set the rectangle and interval width. This keeps one choice of width for the entire counting key even if several raw boxes give the same primitive window.

For raw envelopes above \(\vartheta\), the [stability estimate](ERDOS_364_PELL_AND_JOINT_DETERMINANT.md#the-only-maximizer-of-the-pairwise-bound) at tolerance \(1/1600+10^{-7}\) gives the required strict margins in \(|\tau|\) and \(\Delta\) before window perturbation. Since \(h'\le2F^4\), shrinking \(\eta\) ensures \(h'\le x^\sigma\). Apply the uniform averaged window estimate. If the raw envelope is at most \(\vartheta\), choose a pair and one branch attaining it and apply the original primitive-window bounds of Lemma 8.2. The elementary exponents \(s,t\) and both arguments of \(\max(m_{\rm old},B_{\rm red})\) are Lipschitz under the window perturbations, with \(m_{\rm old}\) the original width exponent. If the reduced branch was available at the raw parameters, Lemma 6.4 gives \(m_{\rm old}\le B_{\rm red}\) there, and both actual totals \(s,t\) remain at least \(1/20\). Thus the exponents change by only the window error even when the gate changes. Choose the small-shift parameter of Lemma 8.2 in terms of \(\varepsilon\), then choose \(\eta\) even smaller to meet that shift range and absorb the exponent perturbation. With these choices each used window contributes at most \(O_\varepsilon(x^{\vartheta+\varepsilon/4})\).

The counting key is the pair \((i,j)\), \(g,q\), and the dyadic window in the three sizes of Lemma 2.3. There are \(O(F^8(\log Y)^3)\) keys. A primitive solution and its key recover \(N+i=g^3q^2 d'^3u'^2/(f_i^2f_j^2)\), so counting all primitive points on a used window is an upper bound. The sum of their upper bounds is \(O_\varepsilon(F^{O(1)}(\log Y)^3Y^{\vartheta+\varepsilon/4})\). We can work with \(\varepsilon\le1\), so the exponent in \(F^{O(1)}\) is absolute. Taking \(\eta\) small enough in terms of \(\varepsilon\) proves the claimed bound, with no dependence on the particular coefficient triple in the constant.

## Transfer to Erdős #137

By [Tao, Lemma 3.1](https://arxiv.org/html/2603.27990v3#S3), if \(\prod_{j=1}^k(n+j)\) is powerful for \(n\) sufficiently large, then \(k\le K_X=\lceil\exp((\log X)^{7/10})\rceil\) uniformly for \(n\le X\). The bounded initial starts have no effect on the count. For each powerful product write \(f_j=\prod_{p<k,\ v_p(n+j)\ {\rm odd}}p\). The quotient \((n+j)/f_j\) is powerful: the small-prime exponents become even, and a prime \(p\ge k\) divides at most one of the factors and therefore occurs there to exponent at least two. [Plymesser, Lemma 8.1](sources/plymesser-powerful-products-hexagon2611.00602-v1.0.1.pdf) gives three consecutive coefficients bounded by \(\max(30030,k^8)\). For large \(X\) there are at most \(O(K_X^{25})=X^{o(1)}\) choices of starting index and coefficient triple, with all coefficients \(\le K_X^8\). The translated three-term start is at most \(2X\). Apply the small-coefficient bound with \(Y=2X\), using loss \(\varepsilon/2\), and sum the choices to obtain

\[
\#\left\{0\le n\le X:\exists k\ge3,\ \prod_{j=1}^k(n+j)\in\mathcal Q\right\}
\ll_\varepsilon X^{\vartheta+\varepsilon}.
\]

This counts starts across all widths \(k\ge3\). It also bounds starts at each fixed width and the number of nonpowerful integers up to \(X\) contained in such products, since each start accounts for at most \(K_X\) offsets after the same width reduction. The estimate does not decide whether any such product exists.

## Pell repair corollary

In the [Pell descent](ERDOS_364_PELL_AND_JOINT_DETERMINANT.md#pell-descent-and-minimal-repair), a nonfundamental center comes from odd \(k\ge3\), with \(t=T_1<(2X)^{1/k}\) and the three exponent-one defect factors \(f_i\mid\operatorname{rad}(k)\). At each fixed \(k=3,5\) the finitely many defect triples fit the small-coefficient theorem for large \(X\), using start \(t-2\). The remaining indices \(k\ge7\) contribute \(O(X^{1/14}(\log(2X))^2)\) by counting only the repaired middle term, as in the elementary argument. The trace \(t\) determines the squarefree part of \(t^2-1\), so no separate kernel count enters. Thus

\[
C(X)-F(X)\ll_\varepsilon X^{\vartheta/3+\varepsilon}
\ll_\varepsilon X^{2/25+\varepsilon}.
\]

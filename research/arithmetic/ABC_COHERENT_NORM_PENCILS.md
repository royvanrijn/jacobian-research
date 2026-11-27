# Cube classes in coherent quadratic norm configurations

## A pencil through two split primes

Use the [prime-path witnesses](ABC_QUADRATIC_THREE_TORSION_FERMAT_TAIL.md#witness-selection) in an imaginary quadratic field \(K=\mathbb Q(\sqrt{-d})\), with \(d>0\) squarefree. Choose one oriented ideal over each distinct odd unramified split prime \(p_i\) in a fixed coset of \(\mathrm{Cl}(K)[3]\). Write \(\alpha_{ij}=(u_{ij}+v_{ij}\sqrt{-d})/2\) for the coherent generators, with \(\alpha_{ji}=\bar\alpha_{ij}\). The identity \(\alpha_{ki}\alpha_{ij}=p_i^3\alpha_{kj}\) allows a change of anchor to any chosen vertex \(k\). At anchor \(1\), put \(\beta_i=\alpha_{1i}\) and \(\beta_1=p_1^3\). For \(i\ne1\), \(\beta_i\) is primitive in the integral lattice of \(K\): its ideal \((\bar{\mathfrak p}_1\mathfrak p_i)^3\) contains just one orientation over each of its supporting rational primes. Thus \(\beta_i\) cannot be divisible by any rational prime in \(\mathcal O_K\).

Fix vertices \(1,2\) and put \(v=v_{12}\), \(u=u_{12}\). For another vertex \(i\), use \(x=v_{1i}\), \(y=v_{2i}\), \(g=\gcd(|x|,|y|)\). The two anchor vectors are independent, since \(v\ne0\). Their determinant coordinates give
\[
 v\beta_i=x\beta_2-y\beta_1,\qquad
 v^2p_i^3=p_2^3x^2-u xy+p_1^3y^2,\qquad g\mid v.                 \tag{1}
\]
For the second equality, take norms and divide by \(p_1^3\); the cross term uses \(\beta_2+\bar\beta_2=u\). For the divisibility, use the integral coordinates of the primitive \(\beta_i\) in the first equality. The quadratic on the right of (1) has discriminant \(u^2-4(p_1p_2)^3=-d v^2<0\). For a proposed extension, the first equality recovers the only possible generator. Integrality, the specified ideal over \(p_1\), the oriented prime-cube ideal over \(p_i\), and the selected signs can then be checked explicitly.

In particular, fixing the two signed pencil coordinates of vertices \(i,j\) fixes their mutual coefficient:
\[
                 v v_{ij}=v_{1i}v_{2j}-v_{2i}v_{1j}.               \tag{2}
\]
This is the Plücker identity with the two anchor columns kept fixed.

The opposite-product obstruction also bounds the density of an adjusted cube class. Fix vertex \(0\) and put \(w_{ij}=v_{ij}/(v_{0i}v_{0j})\) for the other vertices. Plücker gives \(w_{ij}+w_{jk}=w_{ik}\), so these are differences of rational projective coordinates. If three such differences on a triangle had one cube class, they would give a nonzero rational solution of Fermat's cubic equation. On \(m\) vertices the graph of one adjusted class therefore has at most \(m^2/4\) edges. Indeed, in a triangle-free graph the degrees at the ends of every edge sum to at most \(m\). Summing gives \(\sum\deg^2\le m e\), and Cauchy–Schwarz gives \(4e^2/m\le m e\).

## Uniform bound for one prescribed cube class

**Proposition.** Let all primes lie in \([Z,2Z)\), let \(|v_{ij}|\le V\) with \(V\ge1\), and fix two oriented anchors as above. Fix a class \(\chi\in\mathbb Q^\times/\mathbb Q^{\times3}\). The number of distinct possible completion primes \(p_i\) for which \(x/y\) has class \(\chi\) is
\[
 \ll_\varepsilon (ZV)^\varepsilon
 \left(1+Z^{-1/12}\max\{Z,V^{2/3}\}^{1/3}\right).                  \tag{3}
\]
The implied constant is uniform in \(d\), the anchors and \(\chi\). In particular it is \(O_\varepsilon(Z^{1/4+\varepsilon})\) when \(V^{2/3}\le Z\).

Write \(x=gR r^3,\ y=gS s^3\) with \(R,S\) positive cube-free and coprime, \(r\ne0,\ s>0\), after simultaneously negating \(x,y\) if needed. For each prime in the support of \(\chi\) modulo cubes there are at most two choices of which of \(R,S\) it divides, and then the exponent is determined by \(\chi\). If there is any completion, the product of these primes is at most \(V^2\). Together with \(g\mid v\), this gives \(O_\varepsilon(V^\varepsilon)\) triples \(g,R,S\). The sign is absorbed by \(r^3\).

For each such triple, (1) gives a rational point \([rs:p_i:s^2]\) of height \(O(B)\), \(B=\max\{Z,V^{2/3}\}\), on the plane sextic
\[
 p_2^3R^2 X^6-u R S X^3 W^3+p_1^3S^2 W^6-(v/g)^2Y^3W^3=0.      \tag{4}
\]
In the affine chart \(W=1\), the quadratic in \(X^3\) has two distinct nonzero roots over \(\mathbb C\). The sextic on the right of the resulting equation for \(Y^3\) has six simple roots, so it is no cube in \(\mathbb C(X)\). The cubic in \(Y\) is irreducible over \(\mathbb C(X)\), and homogenization introduces no factor \(W\). The smooth projective model has genus four by Riemann–Hurwitz: the cubic cover of the \(X\)-line is totally ramified at the six simple roots, with no ramification at infinity because the degree is six. A point in this chart determines \(r/s\), hence the coprime \(r,s\) with \(s>0\), then \(p_i\) from \(Y/W=p_i/s^2\).

The two pure coefficients \(p_2^3R^2,p_1^3S^2\) are coprime. A prime in their gcd would have to be \(p_1\mid R\) or \(p_2\mid S\), contrary to the norm primitivity \((p_1,x)=(p_2,y)=1\). Thus (4) is a primitive integral form with coefficient norm at least \(Z^3\). [Walsh, Theorem 1.2](https://arxiv.org/pdf/1308.0574v2) bounds the rational points of projective height at most \(O(B)\) on a primitive irreducible degree-six form by \(O(1+B^{1/3}(\log\|f\|+O(1))\|f\|^{-1/36})\). Bounding the logarithm by an arbitrarily small power of \(\|f\|\) proves (3).

When three vertices \(1,2,3\) are fixed, the condition that any specified two of the three opposite-edge products on \(1,2,3,i\) have the same cube class is a condition of this form. For example the equality of the classes of \(v_{13}v_{2i}\) and \(v_{1i}v_{23}\) prescribes the class of \(v_{1i}/v_{2i}\) by the known ratio \(v_{13}/v_{23}\). The other two comparisons follow by changing the anchor pair among \(1,2,3\). The three products cannot all agree in cube class, by (2) and Fermat for cubes. Thus the union of the three classes of coincidences has the same bound (3). All three can be distinct, and two can coincide in actual ideal witnesses.

For comparison, fixing only a quotient cube class in a box for \(x,y\) gives \(O_\varepsilon(V^{2/3+\varepsilon})\) distinct ratios by the same cube-free decomposition, without using (1). The norm lifts the restricted directions to curves and retains the primitivity needed for the coefficient-sensitive saving. A direct complex incidence bound for the group law of the witnesses would have to contend with the structured alternative in [Raz–Sharir–de Zeeuw, Theorem 1.1](https://arxiv.org/pdf/1504.05012): the relation on the norm-one directions is multiplication on a one-dimensional torus, and local logarithms give the additive form in that alternative.

## Local partitions on four and five vertices

Set \(\epsilon_d=1\) for \(d\equiv3\bmod4\) and \(\epsilon_d=2\) otherwise, and \(D_{ij}=v_{ij}/\epsilon_d\). For each rational prime \(\ell\) and \(e\ge1\), the relation \(i\sim j\) given by equality of vertices or \(\ell^e\mid D_{ij}\) is an equivalence relation, with the partitions nested as \(e\) increases. To see this, reanchor at a vertex of norm prime different from \(\ell\). In the integral basis \(1,(1+\sqrt{-d})/2\) or \(1,\sqrt{-d}\), respectively, the determinant of two anchor vectors is \(p_k^3 D_{ij}\). All vectors, including the rational anchor vector, have a unit coordinate over \(\mathbb Z/\ell^e\mathbb Z\). Their determinant vanishes in this ring exactly when one is a unit multiple of the other: solve for the scalar using a unit coordinate, and use the unit coordinate of the second vector to see the scalar is a unit. This proves the partition claim, including at 2. In particular a prime dividing two coefficients along an adjacent pair to a common depth divides the third edge to that depth. Marginal divisibility data alone do not impose these nested partitions.

## Exact examples and the unresolved count

The [example checker](verify_coherent_norm_pencils.py) constructs oriented principal prime cubes for odd primes at most \(151\) in fields with \(d=23,31,59,83\). It selects the first positive-trace generator in order of its positive imaginary coefficient, verifies the selected roots modulo \(p^3\), and uses their signed products. The run with Python 3.12.3, from the repository root using `python3 research/arithmetic/verify_coherent_norm_pencils.py`, found three distinct opposite product classes in all \(5642\) of those quadruples and checked the local partitions on the first five vertices in each field. A separate example uses \(d=43\), primes \(17,59,127,139\), and cube generators with signed coefficient pairs \((130,8),(904,-10),(2860,18),(3224,90)\). Its opposite products in Plücker order are \(-425251944,-583336000,-158084056\). Their positive cube-free parts are \(72917,72917,19760507\). Thus a two-class collision really occurs. These examples use the identity coset of \(\mathrm{Cl}(K)[3]\). Equations (1)–(4) apply in every coset.

The same checker gives an explicit relaxed rank-one record outside the imaginary quadratic norm locus. Let \(\delta\) be the common path gcd, use \(x=v_{12}/\delta,\ y=v_{23}/\delta,\ z=v_{13}/\delta\), and write \(P_i=p_i^3\). From (1) at anchors \(1,3\), we have \(z^2P_2=y^2P_1+x^2P_3+u_{13}xy\). Thus
\[
 x^2y^2z^2 d\delta^2
      =4x^2y^2P_1P_3-(y^2P_1+x^2P_3-z^2P_2)^2.                  \tag{5}
\]
Take \((p_1,p_2,p_3)=(11,13,19)\) and \((x,y,z)=(389,47,299)\). The full congruences in the prime-path proof hold with \(j=2z\). The lattice with ratios of this prime direction modulo \(x,y,z\), respectively, has only \(0\) and \(\pm(11,13,19)\) in the Euclidean ball of radius \(38\). The right side of (5) is \(-700867140113231117\), so this record has no lift for any \(d>0\). This is an exact rejection of one relaxed record.

The bound (3) controls completions of fixed anchor data in a specified cube class. The [rank-one lattice count](ABC_QUADRATIC_THREE_TORSION_FERMAT_TAIL.md#rank-two-and-the-primitive-coefficient-plane) varies the discriminant and the anchor data, and its per-field cosets can be much smaller than \(Z^{1/4}\) in the range that governs the second moment. A power saving on the full set of occupied rank-one records still requires an aggregate restriction across those pencils, or a different selection of paths. The \(17/14\) second-moment exponent and the low joint three-edge thresholds are unaffected by (3).

A saving to \(V^{3-\sigma}\) on the rank-one count would replace its tail arm by \(X^{5/2-\sigma}/H^{8-3\sigma}\), with the same subpower allowances. For \(0\le\sigma\le1/9\) its crossover with the first moment \(X/H\) occurs at \(H=X^{(3/2-\sigma)/(7-3\sigma)}\). At \(\sigma=1/9\) this meets the unchanged rank-two crossover at \(X^{5/24}\), giving the proposed \(29/24\) ceiling for this one-term improvement.

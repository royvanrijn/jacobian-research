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

The following criterion determines whether the normalized three-edge data lift when the field and prime orientations may be chosen.

### Exact reconstruction from three signed coefficients

**Proposition.** Let \(p_1,p_2,p_3\) be distinct odd primes, \(P_i=p_i^3\), and \(x,y,z\) be nonzero signed, pairwise coprime integers with
\[
 (x,p_1p_2)=(y,p_2p_3)=(z,p_1p_3)=1.
\]
Write \(A=P_1y^2,\ B=P_2z^2,\ C=P_3x^2\). Impose the three full congruences
\[
 x\mid A-B,\qquad y\mid C-B,\qquad z\mid A-C.                       \tag{5}
\]
Define
\[
 u_{12}=\frac{A+B-C}{yz},\qquad
 u_{23}=\frac{B+C-A}{xz},\qquad
 u_{13}=\frac{B-A-C}{xy},\qquad
 E=\frac{4AC-(A+C-B)^2}{x^2y^2z^2}.                                \tag{6}
\]
Then \(u_{ij}\) and \(E\) are integers. There are primitive coherent oriented prime-ideal witnesses in an imaginary quadratic field with \(v_{12}=\delta x,\ v_{23}=\delta y,\ v_{13}=\delta z\), \(\delta>0\), exactly when \(E>0\) and \(\gcd(E,p_1p_2p_3)=1\). Their radicand \(d\) and scale \(\delta\) are unique, given by \(E=d\delta^2\) with \(d\) positive squarefree.

For necessity, the nested local partitions above show that the least valuation of the three \(v_{ij}\) occurs at least twice at every rational prime. Thus division by their common positive gcd \(\delta\) makes the three coefficients pairwise coprime. Apply (1) with the three pairs of anchors to get the trace formulas in (6), which give (5). The edge \(13\) norm gives \(E=d\delta^2\). The primes split unramified, so none divides \(d\). Primitivity on the adjacent edges excludes any of the three primes from \(\delta\).

Conversely (5) and pairwise coprimality make the three trace quotients integral. The numerator \(4AC-(A+C-B)^2\) is symmetric in \(A,B,C\), as is seen by expanding it. It follows that
\[
 \frac{4P_1P_2-u_{12}^2}{x^2}
 =\frac{4P_2P_3-u_{23}^2}{y^2}
 =\frac{4P_1P_3-u_{13}^2}{z^2}=E.                                  \tag{7}
\]
The denominator of \(E\) in lowest terms divides all three squares and is therefore one. When \(E>0\), use its positive squarefree decomposition, and set \(\alpha_{ij}=(u_{ij}+v_{ij}\sqrt{-d})/2\). The norm equations (7) give the integral-basis parity conditions: if \(d\equiv3\bmod4\), then \(u_{ij},v_{ij}\) have the same parity, and otherwise both are even. The trace formulas give \(u_{12}y+u_{23}x=2P_2z\). Expanding the symmetric numerator in (6) also gives \(u_{12}u_{23}-Exy=2P_2u_{13}\). These are the two coefficients of \(\alpha_{12}\alpha_{23}=P_2\alpha_{13}\).

The hypothesis on \(E\) and the incident prime units for \(x,y,z\) make each \(\alpha_{ij}\) prime-primitive. Any prime in its norm must split: an inert prime would contribute an even norm valuation, and a ramified prime of norm valuation three would make the element divisible by the rational prime. Exactly one prime ideal above each of the two rational primes occurs, with exponent three. The coherent product cancels both orientations at \(p_2\), so \(\alpha_{12}\) and \(\alpha_{23}\) use opposite ideals there. The orientations at \(p_1,p_3\) then agree with those forced by \(\alpha_{13}\). These three ideals lie in one coset of \(\mathrm{Cl}(K)[3]\) under the directed convention at the start of the note.

The field-window test for the tail proof is therefore \(\operatorname{sf}(E)\in[X,2X)\), where \(\operatorname{sf}\) removes the maximal square factor. The preceding equivalence concerns existence with suitable orientations. Any fixed selection of primes or witnesses in the tail proof further restricts this set. From (7), \(E\le4\min(P_1P_2/x^2,P_2P_3/y^2,P_1P_3/z^2)\). In particular, if one of \(|x|,|y|,|z|\) is comparable to \(V=Z^3/\sqrt X\), then \(\delta=O(1)\) in the field window. The numerator in (6) is positive precisely when \(\sqrt A,\sqrt B,\sqrt C\) are side lengths of a nondegenerate triangle, by Heron's factorization. It is the remaining real test at the level of three edges.

The checker reconstructs the traces and \(E\) on every generated triple, including a separate \(d=2\) test of the integral-basis condition. The positive-Heron control \((p_1,p_2,p_3)=(3,5,7)\), \((x,y,z)=(1,2,1)\) passes (5) but gives \(E=10475=419\cdot5^2\). Here the lift fails the incident-prime unit test at \(5\).

The same checker gives an explicit relaxed rank-one record outside this real locus. Take \((p_1,p_2,p_3)=(11,13,19)\) and \((x,y,z)=(389,47,299)\). The full congruences in the prime-path proof hold with \(j=2z\). The lattice with ratios of this prime direction modulo \(x,y,z\), respectively, has only \(0\) and \(\pm(11,13,19)\) in the Euclidean ball of radius \(38\). The numerator in (6) is \(-700867140113231117\), so this record has no lift for any \(d>0\). This is an exact rejection of one relaxed record.

The bound (3) controls completions of fixed anchor data in a specified cube class. The [rank-one lattice count](ABC_QUADRATIC_THREE_TORSION_FERMAT_TAIL.md#rank-two-and-the-primitive-coefficient-plane) varies the discriminant and the anchor data, and its per-field cosets can be much smaller than \(Z^{1/4}\) in the range that governs the second moment. A power saving on the full set of occupied rank-one records still requires an aggregate restriction across those pencils, or a different selection of paths. The \(17/14\) second-moment exponent and the low joint three-edge thresholds are unaffected by (3).

A saving to \(V^{3-\sigma}\) on the rank-one count would replace its tail arm by \(X^{5/2-\sigma}/H^{8-3\sigma}\), with the same subpower allowances. For \(0\le\sigma\le1/9\) its crossover with the first moment \(X/H\) occurs at \(H=X^{(3/2-\sigma)/(7-3\sigma)}\). At \(\sigma=1/9\) this gives \(29/24\) for the second moment. The [rank-two divisor bound](ABC_RANK_TWO_NORM_PLANE_DIVISORS.md#consequence-for-a-future-rank-one-saving) lowers the rank-two crossover to \(X^{1/5}\). Reaching the corresponding conditional second-moment exponent \(6/5\) would require \(\sigma=1/4\) uniformly down to that lower threshold.

## Census at the rank-one crossover

The [CRT census](census_coherent_norm_rank_one.cpp) counts prime directions and positive odd normalized coefficients near the scale \(V=Z^{5/4}\) of the \(17/14\) rank-one crossover. Its field cutoff is \(X=\lfloor Z^{7/2}\rfloor\), and \(V=\lfloor Z^{5/4}\rfloor\) in the code. It uses one increasing prime triple in \([Z,2Z)\) per direction. It counts all admissible coefficient triples in the indicated box, with the incident prime units, pairwise coprimality and all three divisibilities (5). The prime-path filter on the two cube classes and the choice of witnesses in a fixed field are separate restrictions.

For each \(x,y\), the first two divisibilities give the square-root lists for \(z/y\bmod x\) and \(z/x\bmod y\). The list for modulus \(1\) is a single trivial class. CRT enumerates \(z\) in the chosen box, then the third divisibility is tested exactly. The rank-one test uses the homogeneous lattice with ratios \(p_1/p_2\bmod x,\ p_3/p_2\bmod y,\ p_1/p_3\bmod z\) and Euclidean ball of squared radius \(12Z^2\). To test rank, the program enumerates the middle coordinate and the two forced residue classes in the enclosing cube, discarding points outside the ball. The ball contains the given prime direction, so absence of a nonparallel vector certifies rank one. For \(E>0\), the field-window test enumerates \(\delta\le\sqrt{E/X}\) and checks whether \(E/\delta^2\) is an integer in \([X,2X)\) divisible by no prime square. All cutoff, lattice and norm comparisons use integer arithmetic. The positive annulus here is \([V,2V]^3\), inclusive at both ends.

| \(Z\) | Congruence incidences | Rank one | Positive \(E\), rank one | Prime units, rank one | In field window, rank one |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 16 | 0 | 0 | 0 | 0 | 0 |
| 24 | 0 | 0 | 0 | 0 | 0 |
| 32 | 2 | 1 | 1 | 1 | 0 |
| 48 | 2 | 1 | 1 | 1 | 0 |
| 64 | 8 | 3 | 3 | 3 | 0 |
| 96 | 48 | 18 | 16 | 16 | 0 |
| 128 | 67 | 35 | 29 | 29 | 3 |
| 192 | 273 | 142 | 116 | 116 | 5 |

In the wider box \([1,2V]^3\) at \(Z=128\), there are \(122720\) congruence incidences, \(162\) rank-one incidences and \(4\) meeting all rank-one and field conditions. At \(Z=192\) in the annulus, the positivity and prime-unit tests together remove \(26\) of the \(142\) rank-one records, and another \(111\) fail the field window. These finite controls identify which tests do the work at these scales. The aggregate rank-one gap requires an estimate over growing \(Z\), with all relevant coefficient boxes.

The census also reports the largest trace lift among prime-unit rank-one records and among field-window rank-one records. The [large-lift count](ABC_COHERENT_NORM_LARGE_LIFT_COUNT.md#exact-replay-near-the-crossover) gives selected outputs and uses the full gcd rather than its radical.

Replay with a C++17 GNU compiler from the repository root:

~~~sh
mkdir -p .scratch/prime-path-census
g++ -O3 -std=c++17 -Wall -Wextra -Werror \n    research/arithmetic/census_coherent_norm_rank_one.cpp -o .scratch/prime-path-census/replay
.scratch/prime-path-census/replay
~~~

The [Cayley cubic realization](ABC_COHERENT_NORM_CAYLEY_CUBIC.md) retains the prime denominators, field-window content and three odd square torsor variables of the normalized triangle.

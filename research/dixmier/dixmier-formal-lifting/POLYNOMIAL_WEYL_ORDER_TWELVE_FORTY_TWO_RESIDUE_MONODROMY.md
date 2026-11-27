# Residue covers and monodromy of the order-\((12,42)\) transfer

For the [unsplit planar profile](POLYNOMIAL_WEYL_ORDER_TWELVE_FORTY_TWO_TARGET_CYLINDER.md), the primitive of the leading Jacobian differential gives a square–cube Belyi map. The monodromy of this residue cover survives in the planar pair on a subset of the sheets. Its possible blocks, together with the absence of a quadratic intermediate cover, determine the geometric monodromy of the pair.

Use the notation of the target cylinder, and put \(H=2m\), \(M=5m\), \(l=M-1\). Work over \(\mathbb C\), after embedding the field of definition as in that note. The conclusion over any algebraically closed characteristic-zero coefficient field follows by invariance of geometric monodromy under algebraically closed field extensions. Write \(h=h_0 f^{15}\) with \(f\) monic and squarefree of degree \(2m\), where \(h_0\ne0\) is a scalar. The second resultant has leading \(X\)-coefficient \(K B h^{49/3}\) and \(B'-(h'/6h)B=1\).

## The residue cover and its blocks

**Theorem PlanarFortyTwoResidueCover.**
The position factors of the exact residue equation are the affine-equivalence classes of the classical square–cube Davenport–Zannier pairs with degrees \(2m,3m\). Equivalently, their dessins are indexed by unrooted plane trees with \(2m\) vertices of valence one or three. The degree-\(6m\) residue cover over \(u\) has the same monodromy as the corresponding rational Belyi map. Any nontrivial proper block size of this permutation group is two or three.

At each root \(b\) of \(f\), the equation for \(B\) first forces \(B(b)=0\) and then gives \(B'(b)=-2/3\). Thus \(B=f g\), \(\deg g=3m\), and
\[
                       2 f g'-3 f' g=2.
\]
The polynomials \(g,f\) are coprime and \(g\) is squarefree. Let \(g_0\ne0\) be the leading coefficient of \(g\). Then
\[
       \phi(x)=\frac{g(x)^2}{g_0^2 f(x)^3},\qquad
       \phi'(x)=\frac{2g(x)}{g_0^2 f(x)^4}.
\]
Its partitions at \(0,\infty,1\) are respectively
\[
              (2^{3m}),\qquad (3^{2m}),\qquad (l,1^{m+1}).
\]
Indeed \(\phi(\infty)=1\) and the derivative has order \(5m\) in \(1/x\), so \(\phi-1\) vanishes to order \(l\) there. The other preimages of \(1\) have nonzero derivative. Equivalently \(\deg(g^2-g_0^2 f^3)=m+1\), the sharp square–cube difference degree. Conversely, for a monic \(f\) and \(g\) in such an extremal pair, \(2fg'-3f'g\) is a nonzero constant by differentiating \(\phi\) and reading the degree at infinity. Rescale \(g\) to make the constant two. Even if squarefreeness is not assumed in the extremal pair, the constant identity forces it for \(f,g\).

This is the classical passport of [Pakovich and Zvonkin, Proposition 2.12 and Example 2.13](https://arxiv.org/pdf/1306.4141v1). Use the dessin over the arc from \(0\) to \(\infty\), and contract its valence-two vertices. The resulting plane graph has trivalent vertices and \(m+1\) monogon faces. Each is bounded by a loop at a distinct vertex, since the graph is trivalent. Removing these loops leaves a connected graph with \(2m\) vertices and \(2m-1\) edges, whose valences are one or three. The inverse construction attaches a loop at every leaf and subdivides every edge by a valence-two vertex. The point of index \(l\) fixes position infinity in the equivalence, so the changes of position are affine.

For the residue cover itself take \(s^2=f(x)\) and \(\lambda=c_0 f^2s\), where \(c_0\ne0\) is chosen so that \(\lambda^6=h\) and \(\lambda^{34}=h^{17/3}\). Its function is
\[
                   u=\frac{K B}{\lambda}=\frac{K g}{c_0 f s}.
\]
In particular \(u^2/u_\infty^2=\phi(x)\), with \(u_\infty=K g_0/c_0\). This is a connected quadratic base change of the \(\phi\)-cover since \(\phi\) is a square times \(1/f\) in \(\mathbb C(x)\). Over the three-punctured target for \(\phi\) the subgroup defining the base change is the kernel of the parity character on loops around \(0,\infty\). A loop around \(\infty\) cubed is in the kernel of the permutation representation and has odd parity for the base change. Multiplying by this loop shows that the restricted representation still has the full original image. The residue cover has one \(l\)-cycle at each of \(u=\pm u_\infty\) and no other finite branch value. Over \(0\) the quadratic base change cancels the index two, and over \(\infty\) it leaves the index three.

A system of blocks for \(\phi\) gives, by Lüroth's theorem, a rational decomposition \(\phi=A\circ T\) with right degree \(b\), where \(1<b<6m\). Since \(l>3m\), \(A\) must ramify at \(T(\infty)\). Every point over this value under \(T\) is then ramified in the \(1\)-fibre of \(\phi\), so that point must be infinity. After a change of the intermediate coordinate, \(T\) is polynomial and \(b\mid l\). Since \(b\mid6m\) as well, \(b\mid 6\). Every finite critical value of \(T\) lies over \(0\) or \(\infty\) for \(\phi\), since the sole ramified point above \(1\) is infinity. If a value lies over \(0\) and has any \(T\)-ramification, the whole \(T\)-fibre must have indices two, by the uniform index two of \(\phi\) there. The corresponding statement over \(\infty\) gives indices three. For \(b=6\) these two types contribute three and four respectively to the finite critical multiplicity of \(T\). They cannot sum to \(b-1=5\). This leaves \(b=2,3\). In each of those cases the same count forces one totally ramified finite fibre, and \(T\) is a power map up to affine changes.

At \(m=1\) the tree is a single edge. In the normalization \(f=x^2-1\), the differential equation gives \(g=(2x^3-3x)/3\). The rational map is
\[
             \phi(x)=\frac{x^2(x^2-3/2)^2}{(x^2-1)^3},
\]
with block size two from the right factor \(x^2\).

## Lifting the two residue permutations

**Proposition PlanarFortyTwoResiduePermutationLifting.**
For the finite normalization of the planar target in \(\mathbb C(x,X)\), the residue monodromy occurs as a subgroup transitive on \(6m\) sheets and fixing all other sheets. Its two finite branch generators are meridians of the unique affine branch curve.

The generic \(P\)-fibre accounts for every horizontal branch curve, as in the target-cylinder argument. The outer chart gives one irreducible image, and the high-root charts over finite \(Q\) are unramified. There are no vertical branch curves. Indeed affine branch curves of the finite normalization are in the nonproper value set of the original Keller map, and [Chau, Theorem 3](https://arxiv.org/pdf/0804.3172v3), gives polynomial parametrizations with both coordinate degrees positive. By the birational outer parametrization and the index of the outer chart, a smooth generic meridian is a single \(l\)-cycle. Isolated exceptional points can be avoided by taking \(P\) large.

Choose a closed disc in the target chart \((t,u)\) of the [target-cylinder proof](POLYNOMIAL_WEYL_ORDER_TWELVE_FORTY_TWO_TARGET_CYLINDER.md#a-target-chart-with-the-residue-variable) containing both finite residue branch values in its interior. For small nonzero \(t\), there are at most two branch values in this disc. In fact \(P_0(v)=t^{-12}\) has twelve solutions with \(t v\) approaching the twelve roots of \(p_{12}(t v)^{12}=1\), where \(p_{12}\) is the leading coefficient of \(P_0\). The limiting pair \((t^{42}Q_0,t^{34}W_0)\) changes by \((\epsilon^{42},\epsilon^{34})\) as the root changes by \(\epsilon\in\mu_{12}\). Each allowed limiting pair occurs twice, since \(\gcd(12,42,34)=2\). The disc uses \((1,\delta)\), hence at most two branches enter it. The finitely many nongeneric \(P\)-values where the outer branch image is singular or its projection is critical are also avoided.

Choose based loops in the disc around \(u_\infty,-u_\infty\), avoiding both values. The inverse lifts of these compact loops under the residue cover have finite \(x\), avoid \(f=0\), and have nonzero derivative. With \(w=1/X\), on the residue double cover solve
\[
                     t=\frac{w}{\lambda}(1+O(w))
\]
for \(w\) as in the target-cylinder proof. Analytic continuation along each compact lifted loop, using the implicit function theorem in a finite cover by local charts, gives the same permutation on \(6m\) distinct lifts for all sufficiently small nonzero \(t\). If two residue points have the same \(x\), their two \(\lambda\)-values are opposite and give different leading \(X\)-values, so they too give distinct sheets. The open source maps injectively into the finite normalization by the quasi-finite birational form of Zariski's main theorem.

The restriction of each loop to these sheets is a nontrivial \(l\)-cycle. Therefore each loop encloses a branch value for small nonzero \(t\), and the bound of two branch values makes each a single meridian. Its actual permutation is a single \(l\)-cycle, so it fixes all sheets outside the lifted subset. The two loops generate the stated subgroup because the two finite meridians already generate the residue monodromy.

## The planar monodromy

**Theorem PlanarFortyTwoFullPlanarMonodromy.**
The transferred planar pair of geometric degree \(16m\) has geometric monodromy \(A_{16m}\) for even \(m\) and \(S_{16m}\) for odd \(m\). In particular the degree-sixteen profile has group \(S_{16}\).

Let \(G\) denote the transitive monodromy of the finite normalization. The affine smooth meridians of its unique branch curve normally generate \(G\) and are conjugate. For normal generation, fill a loop in the affine plane and choose the filling transverse to the smooth branch curve, avoiding its singular points and isolated exceptions. Its intersections express the loop as a product of conjugate meridians. Conjugacy follows by joining smooth points along the irreducible curve.

Suppose \(G\) has a nontrivial proper block system of block size \(b\). A meridian acts nontrivially on the blocks, since otherwise their normal generation would make the action on blocks trivial. A single cycle moving blocks moves all sheets in each moved block, so \(b\mid l\). Also \(b\mid16m\). In particular \(b\) is a power of two since \(\gcd(16m,5m-1)\) divides 16.

Intersect the blocks with the \(6m\) sheets of the lifted residue subgroup. These intersections form a block system for that subgroup's transitive action on its support. If the support lies in one block, the two meridians of the lifted subgroup preserve all blocks, contrary to conjugacy and normal generation. Otherwise each block meeting the support is moved by some element of the subgroup, and therefore contains no sheet outside the support, since those sheets are fixed. The intersection size equals \(b\). The residue theorem gives \(b=2\) or \(3\). Only \(b=2\) is a power of two.

A block of size two would give an intermediate field over which the source field has degree two, hence an involution \(\sigma\) of \(\mathbb C(x,X)\) fixing \(P,Q\). This involution acts on the finite normalization and fixes its outer prime \(E\), the only ramified prime over the affine branch curve. We next determine its action on \(x\).

In the terminal source chart put \(r=1/x\) and write \(X=\Lambda(r)+r^M v\) with the finite Laurent shift \(\Lambda\) from the outer Newton path. The polynomials \(P,Q,Z\) are regular in \((r,v)\). Their \(v\)-degrees are \(12,42,99\), and their leading \(v\)-coefficients are units at \(r=0\), by the outer semigroup calculation. For any divisorial pole of \(x\) in the finite normalization, \(r\) has positive value. Regularity of \(P\) forces \(v\) to have nonnegative value, since the highest \(v\)-term has unit coefficient and would dominate otherwise. If \(v\) has constant residue, both \(P,Q\) have constant residues, which is impossible for a divisor over the finite target. Thus \(v\) has transcendental residue, and this is the terminal divisor \(E\). It follows that \(\sigma(x)\) is polynomial on the affine source, since its only divisorial pole in the normalization is outside that open affine plane.

Whenever \(x\) and \(P\) simultaneously tend to infinity along a curve, the terminal variable \(v\) tends to infinity. The leading \(v\)-terms then give the fixed limit
\[
             \frac{Z^{12}}{P^{99}}\ \longrightarrow\
                  \frac{z_{99}^{12}}{p_{12}^{99}},
\]
where \(z_{99}\) and \(p_{12}\) are the leading coefficients of \(Z_0\) and \(P_0\). On the line \(x=\alpha\), as \(X\) tends to infinity, the same ratio instead tends to \(K^{12} B(\alpha)^{12}/h(\alpha)^2\). This is nonconstant in \(\alpha\), since it has a pole at each root of \(f\). If the polynomial \(\sigma(x)\) depends on \(X\), it tends to infinity along a generic such line. Apply the terminal limit to \((\sigma(x),\sigma(X))\), for which \(P,Z\) are unchanged. The two limits would agree for generic \(\alpha\), a contradiction. Thus \(\sigma(x)\) depends only on \(x\).

A polynomial involution of the position line is affine. If \(\sigma\) fixes \(x\), it would give a nontrivial automorphism of the degree-twelve extension over \(\mathbb C(P,x)\), whose group is [\(S_{12}\)](POLYNOMIAL_WEYL_ORDER_TWELVE_FORTY_TWO_PLANAR_FIBRES.md#monodromy-of-the-generic-position-projection). The point stabilizer is self-normalizing, so that extension has no such automorphism. Thus \(\sigma(x)=-x+\beta\). Over \(\mathbb C(x)\), its action on \(X\) is fractional linear, and the unique pole of \(P\) on the \(X\)-line forces the form \(\sigma(X)=a(x)X+c(x)\). The coefficient of \(X^{12}\) in the identity for \(P\) gives
\[
                  a(x)^{12}=\frac{f(x)^{30}}{f(-x+\beta)^{30}} .
\]
Since twelve does not divide thirty, reflection preserves the simple root set of \(f\) and \(a\) is constant. The midpoint \(x=\beta/2\) cannot be a root, since the other roots occur in pairs and \(2m\) is even. Comparing the \(X^{11}\)-coefficients in \(P\) now gives regularity of \(c(x)\) at the midpoint. On that line the affine involution in \(X\) has a fixed point. Differentiating the Keller pair's invariance there would force the differential of \(\sigma\) to be the identity, although its position derivative is \(-1\). This excludes the quadratic intermediate field.

Thus \(G\) is primitive. Its \(l\)-cycle fixes \(16m-l\ge3\) sheets, so [Jones, Corollary 1.3](https://arxiv.org/pdf/1209.5169v1), gives \(G\supseteq A_{16m}\). The meridian is odd exactly when \(m\) is odd. Normal generation gives the asserted group in both parities.

For both orientations of triangular target change, the
[scalar chart calculations](POLYNOMIAL_WEYL_ORDER_TWELVE_FORTY_TWO_TRIANGULAR_FIBRES.md)
give the poles and topology of the generic fibres and the local sheet
counts on the outer boundary.

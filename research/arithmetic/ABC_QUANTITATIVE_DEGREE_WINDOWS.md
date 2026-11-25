# Quantitative degree windows in small-budget abc boxes

Grouping nearby valuation degrees leaves a rational factor of positive
height. A determinant for one grouped root and one wider cofactor
gives explicit margins that allow this height. The estimate uses the
actual residues visible in a full modulus, so a degree interval can
supply a modulus through many coordinates of individually small mass.

## A quantitative congruence estimate for one grouped root

Use a sequence of anatomic boxes with total logarithmic root mass
\(L\to0\). For each box take \(X\) arbitrarily large at fixed truncation
and fixed degrees. The \(o(1)\) in limiting coefficients refers to the
small-budget limit. An arbitrarily small exponent loss at fixed
truncation is implicit below.

Consider rational units \(z,v\) modulo \(n\) with
\[
 vz^q\equiv K\pmod n,\qquad H(z)\le X^{(h+o(1))L},\qquad
 H(v)\le X^{E+o(1)}.
 \tag{QD1}
\]
Here \(q\) is a positive integer, \(K\) is a unit, and the rational height
\(H\) is the maximum of the absolute reduced numerator and denominator.
Write \(w+o(1)=\log_X n\) with \(w>0\). Let \(k\) be an integer at most
every positive valuation \(e_p=v_p(n)\), with \(kL\to\kappa>0\). Put
\[
 Q_n=\prod_{p\mid n}p^{\lfloor e_p/k\rfloor}.
\]
One may use \(Q=Q_n\) or any divisor \(Q\mid Q_n\) as
the residue modulus. Write \(\log_X Q=(b_Q+o(1))L\)
with \(b_Q\ge b\). Then \(b_{Q_n}\le w/\kappa\). If the modulus
is a product of frozen factors \(x_j^{e_j}\) of disjoint prime support,
with \(e_j\ge k\), taking \(Q=Q_n\) allows \(b\) equal to the limiting normalized
root mass of those factors. Indeed, for \(p\mid x_j\),
\(\lfloor e_p/k\rfloor=\lfloor e_jv_p(x_j)/k\rfloor\ge v_p(x_j)\).
This includes factors with prime powers in their roots.

Choose an integer \(M\ge1\), \(r>0\) and rank fraction
\(0<\alpha<1\), and set \(B=hr+EM\). The following strict
conditions suffice:
\[
 r<\liminf qL,\qquad \alpha(M+1)r<\kappa,\qquad
 \frac{b\alpha(M+1)r}{2}>B,\qquad \frac{2B}{1-\alpha}<w .
 \tag{QD2}
\]
All the estimates below that use \(4B+2E<w\) take
\(\alpha=1/2\) and satisfy this last condition with room.
Then (QD1) has at most \(X^{(b_Q+o(1))L}\) pairs,
so \(Q=Q_n\) gives the upper coefficient \(w/\kappa\). One may use a
fixed amplification order \(\ell\), verify (QD2) with \(h,E\) replaced
by \(\ell h,\ell E\), and obtain the cost \(b_Q/\ell\).
These are statements about the limiting exponent and allow arbitrary
fixed-truncation losses.

### The determinant on a residue class

Group the solutions by the residue of \(z\) modulo \(Q\). For each
\(p\mid n\) set \(b_p=v_p(Q)\le\lfloor e_p/k\rfloor\). The lift
\(v'=Kz^{-q}\) agrees with \(v\) modulo \(p^{e_p}\). On the chosen
\(z\)-residue class, the monomial \(z^i(v')^j=K^j z^{i-qj}\)
has a Taylor expansion in \(z-z_0\) about a unit \(z_0\), with integral
\(p\)-adic coefficients. This follows from the binomial expansion for
an integer exponent, including a negative one. The power of
\(z-z_0\) of order \(s\) is divisible by \(p^{b_p s}\).

Use all monomials of degrees at most \(R=\lfloor r/L\rfloor\) in \(z\)
and \(M\) in \(v\). Let \(D_0=(R+1)(M+1)\) be their number and
\(N=\lfloor\alpha D_0\rfloor\). Any \(N\)-row evaluation minor has
valuation at least \(b_pN(N-1)/2\). Expand its columns by Taylor
order. A determinant term with two identical Taylor vectors vanishes,
so the least total order uses \(0,\ldots,N-1\). The error in replacing
\(v'\) by \(v\) is divisible by \(p^{e_p}\), which is at least
\(p^{b_pN}\) since \(N\le k\) for small \(L\) by (QD2).

Clear the rational denominators in each row to order \(R,M\) in the
two variables. The factor is a \(p\)-adic unit at all primes dividing
\(n\). Each resulting entry has size at most \(X^{B+o(1)}\).
The integral determinant has size at most
\(X^{N(B+o(1))}\), and its divisibility by \(Q^{N(N-1)/2}\)
has exponent per row at least \(b\alpha(M+1)r/2+o(1)\). It therefore vanishes.
The rank is less than \(\alpha D_0\). Put \(G=\alpha B/(1-\alpha)\).
On a maximal set of independent rows, evaluate integer coefficient
vectors of height \(X^{G+\varepsilon}\). The image has fewer than
\(\alpha D_0\) coordinate directions of height
\(O_{L,M}(X^{B+G+\varepsilon+o(1)})\).
Since \(\alpha(B+G)=G\), the number of input vectors exceeds
the number of image vectors for each fixed \(\varepsilon>0\) and
large \(X\). Subtract two inputs with the same evaluations to get
a nonzero integer polynomial \(F(z,v)\), of degrees at most \(R,M\)
and coefficient height \(X^{G+\varepsilon}\), on the residue class.

### Lifting the auxiliary equations

Take each primitive irreducible rational factor of \(F\) in turn. Its
coefficient height is at most a constant depending on the degrees
times \(H(F)\), by Kronecker substitution and the univariate factor
height bound. Denote the factor again by \(F\). Factors supported on
\(zv=0\) have no solutions under consideration. Put
\(J=zF_z-qvF_v\). On a torus component where \(J\) vanishes
identically, the nonzero tangent vector
\(z\partial_z-qv\partial_v\) makes \(vz^q\) constant. In characteristic
zero this would give a component of the irreducible curve
\(vz^q=c\), which has \(z\)-degree \(q>R\). Thus \(F=J=0\)
gives only boundedly many torus points by Bézout. Constants here may
depend on the fixed truncation, including \(q,R,M\).

At any other solution \((z,v)\), the numerator of \(J\) has size
\(X^{B/(1-\alpha)+o(1)}\), by the coefficient bound and the same denominator
clearing as above. Let \(d_p=v_p(J(z,v))\). For primes with
\(2d_p<e_p\), Hensel lifts \(F=0,vz^q=K\) to an exact point
congruent to the given point modulo \(p^{e_p-d_p}\). The Jacobian
determinant is \(z^{q-1}J\), whose valuation is \(d_p\).
At such a prime there are boundedly many exact points, again by
Bézout and \(q>R\). Mark the primes where \(2d_p\ge e_p\) as
discarded. On each retained prime record both its exact local lift
and the value \(d_p\). The lifts and discarded-prime flags have only
boundedly many choices per prime at fixed truncation. The additional
\(d_p\) choices contribute at most \(\prod_{p\mid n}(e_p+1)\), the
divisor count of \(n\). Thus the records number \(X^{o(1)}\) as
\(X\to\infty\) at fixed \(L\).

Two original pairs with the same record have common values of
\(d_p\) and agree modulo the product of \(p^{e_p-d_p}\) over
retained primes. The logarithmic weight lost from
\(n\) is at most \(2\sum_p d_p\log_X p\le2B/(1-\alpha)+o(1)\).
Condition (QD2) leaves a modulus of positive numerical weight,
which exceeds \(2hL\) for small \(L\). Cross-multiplication therefore
recovers the narrow rational \(z\) uniquely. With \(z\) fixed,
the original congruence gives \(v\) modulo the whole modulus \(n\).
Since \(w>2B/(1-\alpha)\ge2E\), rational reconstruction also gives
\(v\) uniquely. Summing over the first residue classes proves the
asserted cost. For amplification, color
primes independently for the numerator and denominator of each
rational coordinate. Some assignment keeps \(X^{-o(1)}\) of the
pairs and prevents cancellation in their \(\ell\)-fold coordinate
products. Divisor bounds give multiplicity \(X^{o(1)}\) at fixed
\(\ell,L\). The products satisfy (QD1) with target \(K^\ell\),
heights \(\ell h,\ell E\), and the same \(q\). Taking \(\ell\)-th
roots gives the amplified exponent.

## Explicit excluded degree windows

Use the balanced limiting boxes with \(d=l_A=l_B=l_C=1/3\) from
the [scaled-degree analysis](ABC_WIDE_COFACTOR_CONGRUENCES.md#necessary-conditions-at-normalized-exponent-13).
The mass of a degree window is its limiting normalized root mass.
Closed windows can be enlarged by a vanishing amount, using continuity
endpoints and a diagonal subsequence, to capture their full limiting
mass. In particular a vanishing omitted root mass can be frozen at
vanishing normalized cost regardless of how much numerical height it
carries.

### One window near degree three

Put \(\eta=1/200\) and let both free letters have their root mass in
\([3-\eta,3+\eta]\) except for combined mass at most \(1/100\).
Suppose the third letter has more than \(4/25\) of nonatomic root
mass in this window. Choose a locked block there of normalized mass
\(\mu=4/25\). To construct the block, remove coordinates of large
individual normalized mass and then use successive small ones until
the mass reaches \(\mu+o(1)\). Any coordinate of positive limiting
individual mass contributes to an atom at a limiting scaled degree.
A diagonal choice of the cutoff decreasing to zero thus leaves
enough mass from the nonatomic part.

The modulus weight is between \(w_-=\mu(3-\eta)\) and
\(w_+=\mu(3+\eta)\). Use \(kL\to\kappa=3-\eta\), below the
actual lowest locked degree. Then \(b\ge\mu\). In the free letters
freeze the omitted factors and group the other roots into a signed
ratio \(z\) of base exponent \(q=\lfloor3/L\rfloor\). We have
\(h=1/3\) and \(E=2\eta/3=1/300\), since each retained root
contributes a power correction of scaled magnitude at most \(\eta+o(1)\).
Take \(M=17\) and \(r=1/6\). Then \(B=101/900\) and the determinant
threshold is \(3/25>B\). The available \(\kappa=599/200\) exceeds
\(3/2\), and \(w_--4B-2E=133/5625\). The congruence comes from
the original signed equation reduced modulo the locked factor.
Its variables, together with the frozen factors, give the ratio of
the two complete terms. Primitivity and the signed equation recover
the triple, with only divisor-many representations at fixed degrees.
The resulting count satisfies
\[
 d\le \frac1{100}+\frac4{25}+\frac{w_+}{\kappa}
   =\frac{19799}{59900}<\frac13.
\]
For example this excludes all three measures nonatomic and supported
in \([3-1/200,3+1/200]\).

### Two windows near degrees two and four

Put \(\eta=1/800\), and suppose each letter has mass \(1/6\)
in each of the windows \([2-\eta,2+\eta]\) and
\([4-\eta,4+\eta]\). The windows may contain atoms. Lock the entire
degree-four window on one letter, at cost \(\mu=1/6\).
Then \(w_\pm=(4\pm\eta)/6\), \(b\ge1/6\), and take
\(\kappa=4-\eta\). Group the free roots with base \(q=\lfloor2/L\rfloor\),
using multiplier one on the degree-two window and two on the
degree-four window. The grouped base ratio has \(h=1/2\) and
the cofactor error has \(E=2\eta/3\). After amplification with
\(\ell=2\), use \(h=1,E=1/600,M=47,r=1/12\). This gives
\(B=97/600<1/6\), the latter being the determinant threshold.
Also \((M+1)r/2=2<\kappa\) and \(w_--4B-2E=79/4800\).
The complete ratio count now gives
\[
 d\le\frac16+\frac{w_+}{2\kappa}
   =\frac{9599}{38388}<\frac13.
\]
All omitted mass tends to zero, so any numerical weight outside the
two windows is allowed by this estimate.

## A quantitative real strip for one grouped root

Consider an anatomic sequence with \(L\to0\). Fix \(h,E,s,r,\tau>0\) and a positive integer \(M\). At each fixed truncation let \(q\) be a positive integer, with \(\liminf qL>r\). Consider positive rational pairs \(z,v\) and a real target \(C>0\) such that
\[
 H(z)\le X^{(h+o(1))L},\qquad H(v)\le X^{E+o(1)},\qquad
 \left|\frac{vz^q}{C}-1\right|\ll X^{-\tau}.
\]
As elsewhere, the small-budget errors are read after taking \(X\) sufficiently large for the fixed truncation. Put \(B=hr+EM\). Suppose
\[
 B<\frac{s(M+1)r}{4},
 \qquad
 \frac{s(M+1)r}{2}<\tau,
 \qquad
 4B<\tau.
 \tag{RSW1}
\]
Then the pair count, uniformly in \(C\), is at most
\[
 X^{(s+o(1))L+\varepsilon}
 \tag{RSW2}
\]
for arbitrary \(\varepsilon>0\) at fixed truncation. The same estimate allows \(E=0\) or \(h=0\) under the displayed strict inequalities.

### Dyadic normalization and the auxiliary polynomial

Split the ranges of \(z,v\) into \(O((\log X)^2)\) dyadic boxes and divide each coordinate by the lower dyadic endpoint of its box. The resulting coordinates \(x,y\) lie in \([1,2]\), and their rational heights increase by at most a factor of two. Indeed, multiplying or dividing by the dyadic endpoint increases one numerator or denominator by at most the original other part times two. In this box the corrected exact curve is
\[
 y=g(x)=c x^{-q}
\]
for a positive \(c\) determined by \(C\) and the box. For a nonempty box, \(c\) stays in a compact interval depending only on the fixed \(q\), since \(x,y\in[1,2]\) and \(y=g(x)+O_q(X^{-\tau})\). Thus \(g\) and its derivatives through any fixed order are uniformly bounded on \([1/2,3]\), with constants allowed to depend on \(q\) and the order.

Take \(R=\lfloor r/L\rfloor\) and the rectangle \(x^iy^j\) with \(0\le i\le R\), \(0\le j\le M\). Write \(D=(R+1)(M+1)\), \(N=\lfloor D/2\rfloor\). Partition the \(x\)-range into cells of length \(\delta=X^{-sL}\). Multiply each row of monomial evaluations by the original reduced coordinate denominators to the powers \(R,M\) to obtain integer entries bounded by \(X^{B+o(1)}\). For a minor of size \(N\) in a cell, expand the evaluations \(x^i g(x)^j\) in the single real coordinate about the center of the cell. A Taylor vector of degree \(a\) has weight \(\delta^a\), while the remainder after degree \(N-2\) has weight \(\delta^{N-1}\). The change from the original \(y\) to \(g(x)\) has this same error bound by the second strict inequality in (RSW1). Repeated Taylor vectors contribute zero to the determinant. Any nonzero term has weight at most \(\delta^{N(N-1)/2}\), since the retained Taylor orders must be distinct and every error vector has weight at least \(N-1\). The determinant is then bounded by \(X^{NB+o(1)}\delta^{N(N-1)/2}\) after clearing denominators. The first strict inequality in (RSW1) gives
\[
 sL(N-1)/2\longrightarrow s(M+1)r/4>B ,
\]
so the integer determinant vanishes for \(X\) sufficiently large. These estimates are uniform in the positive target.

The evaluation rank is less than \(N\le D/2\). Pigeonholing integer coefficient vectors with entries of size \(X^{B+\varepsilon}\) on a maximal set of independent rows gives a nonzero integer polynomial \(F(x,y)\) vanishing at every original pair in this cell. To verify the height, if the rank is \(t\le D/2\), the evaluations of vectors of size \(Q\) take at most \((O_D(Q)X^{B+o(1)})^t\) values on those rows. Taking \(Q=X^{B+\varepsilon}\) for each fixed \(\varepsilon>0\) makes the number of input vectors larger. Subtract two with the same values. Every primitive irreducible rational factor \(G\) of \(F\) has height at most \(c_{R,M}H(F)\), by Kronecker substitution and the univariate factor-height bound. Its degrees in \(x,y\) are at most \(R,M\).

### Counting on an auxiliary factor

Fix a factor \(G\) meeting the original positive coordinate box. Work in the torus and put
\[
 J(x,y)=xG_x(x,y)-qyG_y(x,y).
\]
The curve \(G=0\) has no geometric torus component on which \(yx^q\) is constant. Such a component would equal \(yx^q=c_0\) for \(c_0\ne0\), by irreducibility of this latter curve over an algebraic closure. That component has degree \(q\) in \(x\), while \(R<q\). On the generic smooth points of \(G=0\), the vanishing of \(J\) would make the vector \((x,-qy)\) tangent. This nonzero vector spans the tangent line and annihilates the differential of \(yx^q\), so characteristic zero would make \(yx^q\) constant on that component. Consequently \(G=J=0\) contributes boundedly many torus points by Bézout, with the bound depending only on the degrees at the fixed truncation.

At a remaining original point \((x,y)\), the denominator of \(J(x,y)\) divides the product of the two coordinate denominators to powers \(R,M\). Hence
\[
 |J(x,y)|\ge X^{-B-o(1)}.
\]
Write \(\psi(t)=G(t,g(t))\). The height bound for \(G\) and the exact correction yield
\[
 |\psi(x)|\ll X^{B-\tau+o(1)},\qquad
 \left|\psi'(x)-\frac{J(x,y)}x\right|\ll X^{B-\tau+o(1)},
 \qquad
 |\psi''(t)|\ll X^{B+o(1)}\quad(1/2\le t\le3).
\]
For the middle identity, use \(g'(t)=-qg(t)/t\). In particular \(|\psi'(x)|\ge X^{-B-o(1)}\). For a small fixed \(\xi>0\), the derivative varies by \(O(X^{3B-\tau+\xi+o(1)})\) in the interval of radius \(X^{2B-\tau+\xi}\) about \(x\). By \(4B<\tau\) choose \(\xi\) so this change is strictly smaller than the lower derivative bound. The derivative then has constant sign there, and the radius times its lower bound exceeds \(|\psi(x)|\). The intermediate value theorem gives a zero of \(\psi\) in that interval.

The polynomial \(\psi(t)\) is a nonzero Laurent polynomial. Indeed its exponents \(i-qj\) are distinct for \(0\le i\le R<q\), \(0\le j\le M\). Thus it has at most \(R+qM\) real zeros away from zero. Each resulting exact point \((t,g(t))\) lies within \(X^{2B-\tau+\xi+o(1)}\) in both coordinates of its original pair, since \(g\) is uniformly Lipschitz on the enlarged fixed interval. At a fixed dyadic box, distinct rational coordinates of the relevant heights are separated by \(X^{-2(h+o(1))L}\) in the narrow coordinate or by \(X^{-2E-o(1)}\) in the wider coordinate. Since \(B\ge E\) and \(4B<\tau\), reduce \(\xi\) if necessary so the neighborhood contains at most one original pair. Hence each auxiliary factor contributes boundedly many pairs to the cell. Summing the \(O(X^{sL})\) cells and dyadic boxes proves (RSW2).

### A window on the tall terms with a short height gap

Suppose the limiting numerical height of the short letter is at most
\(4/5\). The two tall root measures may have atoms. If their combined
root mass outside \([3-\eta,3+\eta]\), with \(\eta=1/300\), is at most
\(1/20\), freeze the omitted factors and group the rest with
\(q=\lfloor 3/L\rfloor\). This leaves \(h=1/3\) and \(E=1/450\).
The complete tall ratio is \(1+O(X^{-19/100})\) because the terms
differ by the short term. Frozen factors enter the strip target.
Take
\[
 s=1/4,\qquad r=3/50,\qquad M=11,\qquad \tau=19/100 .
\]
Then \(B=2/45<9/200=s(M+1)r/4\), the Taylor precision requirement
is \(9/100<\tau\), and \(4B=8/45<\tau\). The exact ratio
determines the primitive tall values and then the short term.
Thus \(d\le1/20+1/4=3/10\) in this regime.


## Unequal root heights and the diffuse boundary

### A product bound for the admissible wider height

Fix \(m\ge1\) and individual constants \(C_i\ge1\), and put \(Q=\prod_{i=1}^m C_i\). There is \(c_m>0\), depending only on \(m\), such that for all sufficiently large \(k\), the conclusion of the [bounded-root modulus estimate (GW)](ABC_WIDE_COFACTOR_CONGRUENCES.md#a-bounded-number-of-small-roots) holds under
\[
H(x_i)\le n^{C_i/k},\qquad
H(V)\le n^\rho,\qquad
0\le \rho<c_m/Q .
\tag{UW}
\]
Here \(n\) is \(k\)-full and the congruence is \(V\prod_i x_i^{a_i}\equiv K\pmod n\). As in (GW), the lower cutoff for \(k\) depends on the \(C_i\) but not the fixed signed exponents \(a_i\). Its exponent is \(B_m/k+\varepsilon\), with \(B_m=m^2+m\), multiplied by the exact fiber bound.

Keep the narrow-coordinate dimension induction of (GW). On a geometrically integral narrow component of dimension \(j\), choose as the dominant projection \(j\) of the original coordinate functions. This is possible by taking an algebraically independent subset of the coordinate generators of its function field. For a reducible narrow locus, process its components separately. A nonrational geometric component contributes rational points on a proper intersection with a distinct conjugate and descends the dimension. If a component of bounded degree needs new low-height equations after splitting, apply the bounded-degree evaluation-kernel construction of (GW) to its rational sample. Its output is a subset of that component, and either drops dimension or defines the same irreducible variety. The high-degree covering and bounded-degree sample reset depend on \(\max_i C_i\), but only in their fixed constants and the lower cutoff for \(k\). They leave the original narrow coordinates in use.

At the first mixed determinant, let
\[
M=\lceil A_m Q\rceil,\qquad R=\lceil M\rho k+\sqrt{k}\rceil,\qquad
R_i=\lfloor R/C_i\rfloor
\]
for the \(j\) projected coordinates, with \(A_m\) a sufficiently large fixed constant. Use the rectangular support of degree at most \(R_i\) in each such coordinate and at most \(M\) in the wider variable. The number of columns is
\[
(M+1)\prod_{i\ \mathrm{projected}}(R_i+1)
\ge A_m R^j .
\]
After clearing row denominators, every integer entry has absolute value at most \(n^{(j+1)R/k}\), since \(R_i C_i\le R\) and \(M\rho\le R/k\). For a suitably large constant \(A'_m\), a minor of size \(\lceil A'_mR^j\rceil\) has mean \(j\)-parameter Taylor order larger than \(2(j+2)R\). Its last order is \(O_m(R)\). Over the modulus primes the shallow residue classes gain at least that Taylor sum divided by \(2k\) in the exponent of \(n\), so the minor vanishes for large \(n\). Choose \(A_m\) to make this size at most half the number of columns. The rank defect and kernel pigeonhole argument give an integer auxiliary \(G\) with height \(n^{O_m(R/k)}\) on the same support.

The size estimates now use its individual degrees. Evaluating \(G\), or any of its derivatives, at the rational points has numerator and denominator heights \(n^{O_m(R/k)}\). Indeed
\[
\sum_{i\ \mathrm{projected}} R_i C_i/k+M\rho\le(j+1)R/k .
\]
The other equations \(A_i\) at this stage use only bounded-degree narrow-coordinate relations, with coefficients and evaluations of height \(n^{O_{\mathbf C,m}(1/k)}\). Thus the cleared logarithmic Jacobian minors involving \(G\), the \(A_i\), and the signed monomial have height
\[
n^{O_m(R/k)+O_{\mathbf C,m}(1/k)}.
\]
The same bound holds on the bounded-order wider-derivative strata of \(G\), because differentiation preserves the individual degree bounds. Constants in the chosen linear combinations of minors are absorbed at fixed \(k,\mathbf C\).

The rank-deficient projections and pure-narrow auxiliary cuts descend by the same argument as (GW). In the full-rank branches, a sufficiently small bound on \(R/k\) depending only on \(m\) makes the discarded-prime log weight less than half the whole modulus weight and puts the chart Taylor order \(O_m(R)\) below the available Hensel precision. The additive \(O_{\mathbf C,m}(1/k)\) in the Jacobian height is absorbed by taking \(k\) large. At dimension one the rational reconstruction also holds since \(C_i/k\) tends to zero and \(\rho\) is small. As \(R/k=M\rho+O(k^{-1/2})\), all these conditions follow from \(\rho<c_m/Q\) after making \(c_m\) small in terms of \(m\). The number of chart stages per dimension is the same as in (GW), proving (UW).

For integer narrow factors with disjoint prime supports and exponents \(a_i\) of either sign, the exact-fiber argument from (GW) permits overlap with the wider multiplier. At fixed exact value of the monomial, the mass of each signed narrow factor beyond the ceiling of the matching valuation of that value has its \(|a_i|\)-th power dividing the appropriate numerator or denominator of \(V\). The fiber count is therefore \(n^{\rho\sum_i1/|a_i|+\varepsilon}\). On taking \(\ell\)-fold products after coordinate coloring of the narrow factors and sign coloring of the wider factors, the wide exponent becomes \(\ell\rho\), and the narrow height constants \(\ell C_i\). If the product criterion (UW) holds at these amplified parameters, the count exponent for the original tuples is at most
\[
\frac{B_m}{\ell k}+\rho\sum_i \frac1{|a_i|}+\varepsilon .
\tag{UW1}
\]

### A degree-window inequality

In the limiting scaled-degree measure notation of the balanced note, fix a positive-height full modulus of normalized cost \(\mu\), numerical weight \(w\), and fullness \(kL\to\kappa>0\). On the two free letters choose finitely many disjoint compact continuity intervals of positive normalized degrees. For interval \(I_i\) on letter \(T_i\), choose a representative normalized degree \(q_i>0\), put \(r_i=\nu_{T_i}(I_i)\), and allow also shrinking windows at positive atoms, with the atom mass as \(r_i\). Include the zero windows of both free letters. Write \(r\) for the total of all these retained root masses. Set
\[
E=\sum_i\int_{I_i}|q-q_i|\,d\nu_{T_i}(q),
\]
where exact atom windows have contribution zero. Fix the remaining coordinates. The retained groups give coprime narrow factors with degree \(\lfloor q_i/L\rfloor\) and height exponent \(r_iL+o(L)\). All rounding and zero-window factors combine into a wider signed rational multiplier of height \(X^{E+o(1)}\).

Choose \(C_i>\max(1,\kappa r_i/w)\) and an amplification order \(\ell\). If
\[
\ell E/w < c_m/\prod_i (\ell C_i),
\]
the wider factor satisfies (UW) at the amplified parameters. The complete ratio determines the two coprime values, and (UW1) gives the normalized count bound
\[
d\le l_U+l_V-r+\mu
   +\frac{B_m w}{\ell\kappa}
   +E\sum_i\frac1{q_i}.
\tag{UW2}
\]
All choices on the right are fixed before the small-budget limit. In particular the height constants here reflect the actual chosen modulus weight.

### A restriction at equality of the nonatomic total

Suppose a balanced limiting profile has total nonatomic mass \(c_A+c_B+c_C=1/2\). The three pair inequalities force \(c_A=c_B=c_C=1/6\). Consider any two letters \(U,V\) with finitely many positive atomic scales, and let \(p\) be the total number of their positive atomic scales on the two letters, counting a shared positive scale once per letter. If \(p=0\), this profile is impossible. For \(p\ge1\), their nonatomic measures satisfy, on compact positive-scale intervals away from those atoms, the local mass bound
\[
\nu_{T,\mathrm{na}}(I)\ll |I|^{1/p}\qquad(T=U,V).
\tag{UW3}
\]
For \(p=1\), this entails locally bounded densities on those intervals.

To prove it, free all atom windows on \(U,V\) along with their zero windows. Also choose a small interval \(I\) on one of those letters in a compact scale range positively separated from zero and the finite atom lists. Choose continuity endpoints and let \(t=\nu_{T,\mathrm{na}}(I)>0\), \(\sigma=|I|\). Free it as one additional group at its midpoint. Thus \(m=p+1\), the combined retained cost reduction beyond the atoms is \(t\), and \(E\le \sigma t\). The combined frozen cost on \(U,V\) is \(1/3-t\).

The modulus letter has nonatomic mass \(1/6\). Fix a compact positive range \([a,b]\) containing positive nonatomic mass in its interior. For all sufficiently small \(t\) choose a finite union of continuity intervals in this range of total mass \(\mu\in[t/4,t/2]\). Indeed take an atom-free measurable set of mass \(t/3\), a compact subset of mass at least \(t/4\), and an open neighborhood of total mass below \(t/2\). A finite cover of the compact subset by continuity intervals inside that neighborhood suffices. The degree factors in this union give a modulus of weight \(at/4\le w\le bt/2\). Choose \(\kappa=a/2\), below its fullness.

The height constants for the atom groups can be taken \(O(1/t)\), while for the extra group they are \(O(1)\), uniformly in the location of \(I\) inside the chosen compact positive range. Choose the fixed \(\ell\) large enough that \(B_m w/(\ell\kappa)<t/8\). The product of the amplified height constants is then \(O(t^{-p})\), with a constant independent of the small interval. Since \(\ell E/w=O(\sigma)\), the wider-height criterion holds whenever \(\sigma\le c t^p\) for a sufficiently small positive constant \(c\). Also \(E\sum_i1/q_i<t/8\) for all sufficiently small \(\sigma\), by the fixed positive lower scales. Then (UW2) would give
\[
d\le1/3-t+t/2+t/8+t/8<1/3.
\]
For \(p=0\), any sufficiently short positive-mass interval works, a contradiction. The mass of small intervals in the compact range tends uniformly to zero with their length, since the measure there has no atoms. For \(p\ge1\), the contradiction gives \(t\ll\sigma^{1/p}\), with constants uniform on the compact separations just specified. For discontinuity endpoints, enlarge the intervals by arbitrarily small amounts to continuity endpoints. This proves (UW3).

When \(p=1\), the uniform interval estimate is linear in interval length. It extends to open subsets of each compact range by disjoint intervals and to Borel subsets by regularity, yielding absolutely continuous nonatomic restrictions of locally bounded density. This excludes a singular nonatomic component in this particular finite-atom configuration.

## An explicit positive-degree mass in terms of numerical height

Let \(d=1/3\), fix any letter \(T\), and put
\(\beta_T=\nu_T((0,\infty))\). Then
\[
 \beta_T+c_T>\frac{h_T}{178},\qquad
 \beta_T>\frac{h_T}{356}.
 \tag{QD3}
\]
If \(u\) is the smaller numerical height of the two other letters,
one also has
\[
 \beta_T+c_T\ge\frac{1-u}{3(1+u)} .
 \tag{QD4}
\]
The proof uses prelimit numerical weights and allows escaped
first moments.

Write \(H=h_T>0\), which is positive by the
[short-height bound](ABC_SHORT_HEIGHT_AND_REAL_GAP.md#a-strict-saving-when-the-short-term-is-small),
and let \(U,V\) be the other letters, of heights \(u,1\).
At least one of \(H,u\) is \(1\), since the two tall letters
have height \(1\).

### Freeing a pooled upper quantile

Fix coordinate factors \(A_0,B_0,C_0\) on the three letters.
For any two primitive triples with these fixed factors,
\(ab'-a'b\) is divisible by \(A_0B_0C_0\): divisibility by
\(C_0\) follows by substituting \(c=a+b,\ c'=a'+b'\), and
the three factors are pairwise coprime by primitivity.
The same determinant, up to sign, can be formed from any
pair of letters. Using \(T\) and the letter of height \(u\)
bounds its magnitude by \(X^{H+u+o(1)}\). If the product of
the fixed factors has strictly larger numerical exponent, the
determinant must vanish. Then the ratios agree and primitivity
recovers the entire triple. This is the uniqueness case of
the three-factor lattice estimate.
[Bernert–Browning–Lichtman–Teräväinen, Proposition 4.1](https://arxiv.org/html/2410.12234v2#S4)

Put \(s=\beta_T+c_T\). As \(c_U,c_V\ge1/3-c_T\), the total
atomic mass on each free letter is at most \(c_T\).
Each single prelimit root coordinate on \(U,V\) therefore
has normalized mass at most \(c_T+o(1)\). A coordinate whose
mass exceeds \(c_T\) by a fixed amount would have bounded
scaled degree by the first-moment bound, and a subsequence
would give an atom of such mass.

On \(U,V\) pool the prelimit degree measures into \(\omega_L\).
Its mass \(m_L\) tends to \(2/3\) and its first moment \(J_L\)
to \(1+u\). Let \(F_L(t)\) be the numerical weight of the
highest scaled degrees accounting for root mass \(t\),
allowing a fraction of the last coordinate. For any fixed
\(0<t<1/3-s\), round that last coordinate up and freeze these
highest degrees. This costs at most \(t+c_T+o(1)\). On \(T\),
freeze all degrees above a small positive continuity cutoff
\(\epsilon_0\), at cost at most \(\beta_T+o(1)\) and numerical
weight at least \(H-\epsilon_0/3+o(1)\).
If \(\limsup F_L(t)>u\), choose \(\epsilon_0\) sufficiently
small and take the subsequence attaining the excess. The product
of the three fixed factors then has numerical exponent exceeding
\(H+u\), so the uniqueness count gives \(d\le t+s<1/3\).
It follows that
\[
 \limsup F_L(t)\le u\qquad(0<t<1/3-s).
 \tag{QD5}
\]

### Centering the free pair at its median

Top-quantile averages decrease with their mass, so
\(F_L(m_L/2)\le (m_L/(2t))F_L(t)\). The top half carries at
least half the first moment. Taking limits and letting \(t\)
increase to \(1/3-s\) gives, for \(s<1/3\),
\[
 (1+u)/2\le \frac{u}{1-3s}.
\]
This yields (QD4); the case \(s\ge1/3\) satisfies it
directly.

Now suppose \(s\le s_0=H/178\) and fix
\(t=1/3-s_0-H\cdot10^{-7}\) in (QD5).
For each prelimit pooled measure choose a median \(a_L\).
The absolute deviation about the median is exactly
\(2F_L(m_L/2)-J_L\), by splitting at the median with a
fractional tie if needed. The lower half has weight
\(J_L-F_L(m_L/2)\) and is bounded above by \(a_Lm_L/2\).
The top half has weight at least \(a_Lm_L/2\).
Thus \(a_L\) lies in a compact positive range for small \(L\):
with \(s_1=s_0+H\cdot10^{-7}\), its limiting lower and upper
bounds are \(3(1+u-u/(1-3s_1))\) and \(3u/(1-3s_1)\).
The lower bound exceeds \(2\). Pass to a subsequence with
\(a_L\to a\).

Group every free root into a single signed ratio of base
exponent \(q=\lfloor a/L\rfloor\). The sum of absolute
correction exponents, in numerical units, is at most
\[
 \frac{2u}{1-3s_1}-(1+u)+o(1)
\]
by the median estimate and \(qL-a_L=o(1)\).
The signed sum of the corrections in the ratio tends
to \(u-1\): subtract the two numerical heights, and the
common base exponent times the difference of the two root
masses, which tends to zero. Thus the larger of the
positive and negative correction exponents is at most
\(3us_1/(1-3s_1)+o(1)\).
The rational cofactor can only decrease in height upon
reduction. The pair in (QD1) therefore has \(h=1/3\)
and \(E=(343/20000)H\), since
\(3us_1/(1-3s_1)<(343/20000)H\).

### Counting after locking the high degrees

Use the whole \(T\)-coordinates of scaled degree at least
\(\kappa=(3/2)H\) as a modulus, and \(k=\lceil\kappa/L\rceil\).
Pass to a subsequence for the limiting mass \(\mu\) and weight \(w\).
Then \(\mu\le\beta_T\le s_0\) and \(H-\kappa s_0\le w\le H\).
The bounds hold at an atom of scale \(\kappa\) as well, by using
continuity cutoffs on both sides of \(\kappa\).
The full residue modulus \(Q_n\) has normalized mass at least
\(\max(\mu,w/\kappa-\mu)\ge w/(2\kappa)\), by the root contribution
and the valuation floor loss. This is at least \(353/1068>161/500\).
Choose a divisor \(Q\) of \(Q_n\) by adding copies of whole
root factors up to the multiplicities \(\lfloor e/k\rfloor\)
until its nominal normalized logarithmic mass reaches \(161/500\).
Each increment has mass at most \(\mu+o(1)\). The resulting
record mass \(b_Q\) therefore lies in
\([161/500,161/500+\mu]\) after taking limits.

Apply (QD2) with \(\alpha=543/1000,\ M=6,\ r=(37/100)H\).
We have \(B\le(6787/30000)H\). For \(b=161/500\) the
determinant threshold is \((22642557/10^8)H>B\),
and \(\alpha(M+1)r=(140637/10^5)H<\kappa\). The terminal
condition holds since
\[
 \frac{2B}{1-\alpha}\le\frac{6787}{6855}H
   <\frac{353}{356}H\le w.
\]
The ratio determines the primitive triple. Including the
modulus freeze gives \(d\le161/500+2/178<1/3\), a contradiction.
This proves (QD3); its second inequality follows from
\(c_T\le\beta_T\).

## An explicit diffuse sector with a selective residue modulus

Use normalized fixed anatomic boxes. Suppose a fixed block of coordinates on letter \(T\) has scaled degrees at least \(9/5\), root-mass cost at most \(1/20+o(1)\), and numerical weight at least \(49/50+o(1)\). Assume every coordinate in the block has individual normalized root mass at most \(1/1000+o(1)\). Write \(n\) for the product of their frozen factors. On the other two letters \(U,V\), group all their roots at a common base exponent \(q=\lfloor a/L\rfloor\), for some fixed \(a>9/25\). Suppose this gives a grouped rational root \(z\) of narrow normalized height at most \(1/3+o(1)\) and a rational correction \(v\) of numerical height at most \(1/75+o(1)\). As usual the fixed exceptional coefficients enter the unit target.

Then the normalized count exponent is at most
\[
 \frac{321}{1000}.
 \tag{DS1}
\]
In particular this is a strict small-budget saving on the sector defined by the stated bounds. The conclusion is stable under sufficiently small perturbations of these parameters, with arbitrarily small changes to the numeric count coefficient.

Indeed use \(kL\to\kappa=9/5\). The complete shallow residue modulus has capacity at least \(w/(2\kappa)\ge49/180>27/100\). Choose a selective record of mass \(b\) between \(27/100\) and \(271/1000\), by adding root copies one at a time until the first threshold is reached. Every prime depth in the selected record is at most its original modulus valuation divided by \(k\). The \(271/1000\) bound uses the assumed per-coordinate mass cap. Keep the full \(n\) in the congruence and terminal reconstruction.

Take rank fraction \(\alpha=21/40\), wide degree cap \(M=8\), and \(r=9/25\). The grouped congruence conditions are
\[
 \begin{split}
 B&\le \frac13\frac9{25}+\frac8{75}=\frac{17}{75},\\
 \alpha(M+1)r&=\frac{1701}{1000}<\kappa,\\
 \frac{b\alpha(M+1)r}{2}&\ge\frac{45927}{200000}>\frac{17}{75},\\
 \frac{2B}{1-\alpha}&\le \frac{272}{285}<\frac{49}{50}\le w.
 \end{split}
\]
The \(r<a\) condition follows from the choice of the base. The count exponent is therefore at most the actual record mass \(b\le271/1000\) for each fixing of the block. Each exact pair \(z,v\), with frozen factors restored, gives the full \(U,V\) ratio and hence the primitive triple. Restoring the freeze cost at most \(1/20\) gives (DS1). The usual fixed-truncation and anatomic boxing losses can be absorbed into an arbitrarily small exponent loss.

For use as a direct coordinate test on the two free letters, put
\[
 P_W(a)=\sum_e \beta_{W,e}(eL-a)_+,\qquad
 N_W(a)=\sum_e \beta_{W,e}(a-eL)_+ .
\]
The wider height bound holds if
\[
 \max\{P_U(a)+N_V(a),N_U(a)+P_V(a)\}\le1/75+o(1).
\]
The narrow height is at most the larger free-letter normalized root mass. This includes opposite-sign discrepancies separately, rather than summing their absolute deviations. For free letters of equal normalized root mass and equal numerical height, the two displayed wider contributions agree in the limit.

There is a quantile form of this exclusion at the balanced
endpoint. Keep the stated block on \(T\), let \(u\) be the smaller
height of the free letters, and let \(F\) be their pooled weak
upper quantile as in the [residual quantile formula](ABC_RESIDUAL_QUANTILE_SPECTRUM.md#the-residual-quantile).
If \(e_U,e_V\) are the escaped first moments, then
\[
 F(1/3)+e_U+e_V>u+1/75
\]
is necessary. Indeed, at a pooled prelimit median the absolute
correction equals twice the top-half weight minus the total
moment. The signed imbalance of the two letters is the
difference of their heights, by balance of their root masses.
Thus the larger of the two one-sided correction heights tends
to \(F(1/3)+e_U+e_V-u\). If this is at most \(1/75\),
the bottom-half moment is at least \(1-1/75\).
The median is then at least \(3(1-1/75)\) in the limit,
so a fixed subsequential median satisfies \(a>9/25\).
The above count yields the contradiction.

### An all-nonatomic formal profile in this sector

Give each of \(U,V\) normalized mass \(1/3\), putting fraction \(37/40\) uniformly on \([3-1/200,3+1/200]\) and fraction \(3/40\) uniformly on \([2,4]\). Both first moments are one. Give \(T\) mass \(17/60\) uniformly on \([0,1/50]\), and mass \(1/20\) uniformly on the interval of radius \(1/10\) centered at \(5983/300\). This also has first moment one. All three measures are nonatomic. On \(T\) lock the high interval, whose numerical weight is \(5983/6000>49/50\), whose mass is \(1/20\), and whose limiting individual coordinate mass is zero. On \(U,V\) use the base \(a=3\). The two signed correction exponents are both
\[
 \frac13\left(\frac{37}{40}\frac1{400}+\frac3{40}\frac12\right)
 =\frac{637}{48000}<\frac1{75}.
\]
These measures have no escaped first moment and meet the sector with strict determinant margins.

The earlier one-window bound near degree three uses a local
modulus and a smaller omission allowance on the free letters.
Here the modulus comes from the distant high band, and the freely
grouped letter measures have more than \(1/100\) mass outside the
narrow interval. The example also meets the three-factor upper
quantile constraints. The weight selectable from total root
mass \(1/3\) is at most the whole high \(T\)-weight plus the
top mass \(17/60\) of the opposite pair, which is bounded by
\[
 \frac{5983}{6000}+3\frac{17}{60}+\frac{637}{48000}<2 .
\]
All three measures are atomless, so this verifies the
[upper-third criterion](ABC_RESIDUAL_QUANTILE_SPECTRUM.md#three-letters-without-initial-locks).

## A residue count after locking outliers on the free letters

The [fixed-factor pack bound](ABC_RESIDUAL_QUANTILE_SPECTRUM.md#one-coordinate-rounding-and-the-resulting-inequalities) uses the total numerical weight of all frozen factors. The [one-root congruence count](#a-quantitative-congruence-estimate-for-one-grouped-root) can also exploit the distribution of the factors left free after such locks. In particular one can lock a full modulus on \(T\) and high numerical outliers on its two opposite letters, then take a median of the residual free degrees. The explicit sector below gives the normalized coefficient \(313/1000\) and excludes a formal profile at \(\Gamma=1/2\) satisfying every three-factor pack test, even with fractional atom selection.

### Residual dispersion after arbitrary opposite locks

Freeze prescribed whole original root coordinates on \(U,V\). Let \(\sigma_U,\sigma_V\) be the weak degree measures of their remaining root coordinates, with masses \(m_U,m_V\). Let their full prelimit numerical weights tend to \(u_U,u_V\), and write \(e_i=u_i-\int r\,d\sigma_i(r)\) for their residual escaped first moments. The selected original root sets depend on the fixed boxes rather than on the values of roots within a box. On \(T\), lock a full modulus from a prescribed block.

For a fixed base \(a>0\), group every remaining opposite root at common base exponent \(q=\lfloor a/L\rfloor\), using the unweighted product on each free letter. Their signed ratio is \(v z^q\), with the fixed factors on \(U,V\) reinstated as a known multiplier. The \(U,V\) fixed factors are units modulo the locked \(T\) factor. The narrow height constant in (QD1) can be taken to be \(m_{\max}=\max(m_U,m_V)\) up to the vanishing box errors. The following wider exponent suffices:
\[
 \mathcal E(a)=\frac12\left\{
    \int |r-a|\,d(\sigma_U+\sigma_V)(r)
    +e_U+e_V+
    |(u_U-u_V)-a(m_U-m_V)|
    \right\}.
 \tag{CM1}
\]
Indeed the prelimit sum of absolute correction weights is
\(\int|r-qL|\,d(\sigma_{U,L}+\sigma_{V,L})\). It converges to the first three terms inside braces because
\(|r-a|=r+a-2\min(r,a)\), so only the full moment contributes escape. The imbalance between numerator and denominator contributions is the absolute net signed correction weight. It tends to the final term inside braces. Hence the larger of the two contributions is half their sum plus half their imbalance. This proves (CM1).

When \(m_U=m_V\), select a pooled median of the residual degree measure to minimize (CM1). With unequal masses, one can take a median of the residual pool augmented by an atom of mass \(|m_U-m_V|\) at \((u_U-u_V)/(m_U-m_V)\), restricting to positive bases admissible for the determinant. This follows by rewriting the last term as the absolute deviation from the indicated atom. At a common median and equal residual masses, the minimal wider exponent equals the numerical moment of the top half of the pooled residual root mass, with escape included, minus \(\min(u_U,u_V)\).

More generally one can assign fixed nonnegative integer multipliers \(k_j\) to finitely many residual subblocks \(j\), still with one common exponent \(q\). In \(z\), put each root to its assigned multiplier, and in \(v\) put the rest of its original exponent. The same calculation uses \(\sum_j \int |r-k_j a|\,d\sigma_j\) and signed imbalance \(\left|\sum_j\epsilon_j(u_j-k_j a m_j)\right|\), where \(\epsilon_j\) is the sign of its letter in the ratio. The same escaped moments add to the absolute sum. The narrow height uses the larger of the two weighted root masses \(\sum_{j\in W} k_j m_j\). For pieces with \(k_j>0\), choosing \(a\) amounts to a weighted median on the rescaled measures \(k_j(r/k_j)_*\sigma_j\), with the imbalance atom determined by the signed sums of the numerical weights and weighted masses. The grouped ratio still uses one narrow variable.

### An explicit count with locks on all three letters

Assume a chosen \(T\) block has scaled degrees at least \(2\), normalized root mass at most \(41/1000\), numerical weight at least \(49/50\), and largest individual coordinate mass at most \(1/1000\), all allowing vanishing box errors. Suppose root locks on the opposite letters have combined cost at most \(31/1000\). Require a common base \(a>9/25\) on the residual free pair so that \(m_{\max}\le1/3\) and \(\mathcal E(a)\le11/1000\).

Then the normalized count exponent is at most
\[
 \frac{313}{1000}.
 \tag{CM2}
\]
To prove this, use fullness \(kL\to\kappa=2\) on the frozen factor on \(T\). The complete shallow record capacity is at least \(w/(2\kappa)\ge49/200>6/25\). Add whole root copies to a selective depth record until reaching normalized mass \(6/25\). The overshoot is at most \(1/1000\), so \(6/25\le b\le241/1000\). The complete modulus is still used for the congruence and terminal reconstruction. Take rank fraction \(\alpha=13/25\), wide degree cap \(M=9\), and \(r=9/25\). Then the one-root bounds are
\[
 \begin{split}
 B&\le \frac13\frac9{25}+9\frac{11}{1000}=\frac{219}{1000},\\
 \alpha(M+1)r&=\frac{234}{125}<2,\\
 \frac{b\alpha(M+1)r}{2}&\ge\frac{1404}{6250}>\frac{219}{1000},\\
 \frac{2B}{1-\alpha}&\le\frac{73}{80}<\frac{49}{50}\le w .
 \end{split}
\]
The first-coordinate reconstruction requires only positive retained modulus weight for \(z\), whose numerical height tends to zero. After recovering \(z\), recover \(v\) from the original full modulus. The \(r<a\) condition also holds. The chart cost is therefore at most \(241/1000\). Each counted pair restores the complete primitive triple via the frozen free factors and the signed equation, so the added freezes give a cost at most \((41+31+241)/1000\). This proves (CM2). At fixed budget all loss terms are handled by the usual anatomic boxing; the determinant margins permit sufficiently small fixed perturbations of the data.

### A formal profile with narrow cores and high outliers

Set \(a=510/191\). Each of \(U,V\) consists of the following measure.

- an atom of mass \(1/6\) at \(a\);
- diffuse mass \(91/600\) uniform on \([a-1/10,a+1/10]\);
- diffuse mass \(3/200\) uniform on \([9,11]\).

On \(T\) use these three components.

- an atom of mass \(1/6\) at zero;
- diffuse mass \(1/25\) uniform on \([491/20-1/10,491/20+1/10]\);
- diffuse mass \(19/150\) uniform on \([0,27/95]\).

Each total mass is \(1/3\), each nonatomic mass is \(1/6\), and all three first moments are one. In \(U,V\), the high outlier weight is \(3/20\), and the rest has mass \(191/600\) and weight \(17/20\). On \(T\) the high band has weight \(491/500\), and the low diffuse band has weight \(9/500\). All three supports are bounded and there is no escaped weight.

Lock the high bands on all three letters. The \(T\) block costs \(1/25<41/1000\), has weight \(491/500>49/50\), and its scaled degrees exceed two. Its individual root masses tend to zero because the entire limit on that band is nonatomic. The two opposite freezes together cost \(3/100<31/1000\). Their residual narrow heights are \(191/600<1/3\). At their common median \(a\), the wider exponent is exactly \(91/12000<11/1000\), by equal residual profiles and the uniform core. Thus this formal profile lies strictly inside the coefficients of (CM2), and is excluded at \(d=1/3\).

Before the free locks, the minimum wider exponent over *all common unweighted bases* is
\[
 \frac3{200}(10-a)+\frac{91}{12000}>\frac{11}{100}>\frac1{75}.
\]
Indeed \(a\) is a median of the full profile as well, because the atom at \(a\) covers its half-mass point. The opposed full profiles are equal, so there is no signed imbalance at any common base and the median minimizes the wider exponent.

The formal profile also satisfies all plain three-factor packing tests. Even allowing fractional choices from all atoms, its top total root mass \(1/3\) consists of the three high bands, then the above-median half of the two diffuse cores, then enough mass of the central atoms to fill at scale \(a\). The numerical weight is
\[
 \frac{641}{500}+\frac{79}{300}\frac{510}{191}+\frac{91}{12000}
 =\frac{45673250}{22920000}<2.
\]
Every fixed-factor selection of normalized root mass at most \(1/3\) is therefore below the numerical threshold for three-factor uniqueness. The count (CM2) uses the congruence on the residual pair after the outliers are frozen.

## An explicit lower bound for each letter's atomic weight

At normalized exponent \(d=1/3\) and minimum nonatomic total
\(\Gamma=1/2\), every letter \(T\) satisfies
\[
 \alpha_T:=\int r\,d\nu_{T,\mathrm{at}}(r)>\frac1{500}.
 \tag{AM1}
\]
The [complete-label support bound](ABC_DIFFUSE_FACTOR_GRAPH.md#consequence-for-all-three-nonatomic-supports)
puts the nonatomic measure in \([0,6]\), with mass \(c=1/6\).
The full moment is one by the [minimum-mass conditions](ABC_VARIABLE_MESH_DEGREE_GROUPING.md#consequences-at-the-minimum-total-nonatomic-mass),
so the nonatomic deficit \(\int(6-r)\,d\nu_{T,\mathrm{na}}\)
equals \(\alpha_T\).

Suppose \(\alpha_T\le1/500\). The nonatomic mass of \(T\)
in \([5.95,6]\) is at least \(c-20\alpha_T>0.126\).
Choose two fixed rows of the complete-label graph with
a common column exponent arbitrarily close to \(c\).
The complete \(T\)-terms on those rows have a fixed
nonzero difference \(D\), because the two complete row
values are distinct and the complete column value is the
same on the two rows. This uses only the signed additive
equation, so either addend can be the target \(T\).

Discard the columns in which the outside product on the
second \(T\)-term has radical overlap of mass at least
\(1/10000\) with \(D\). The
[divisor-allocation estimate](ABC_DIFFUSE_FACTOR_GRAPH.md#counting-shared-prime-divisors)
and pair uniqueness bound their exponent by
\(c+\delta-1/10000+o(1)\).
On this row select an outside block of degrees in
\([5.95,6.01]\) of root mass between \(0.083\) and \(0.0831\),
up to limiting errors. The nonatomic mass in the band is
more than enough, and choosing small window loss gives small
individual outside coordinate masses for rounding. Choose
continuity endpoints between \(5.95\) and \(5.96\), and
between \(6\) and \(6.01\). The omitted lower degrees lose
at most \(25\alpha_T\) of the nonatomic mass, so the block
still fits. Freeze the block and strip the full powers
of primes dividing \(D\). Bounded-scale roots are squarefree,
so the retained base mass is at least \(0.0829+o(1)\).
Writing \(n\) for the stripped full modulus, its numerical
weight \(w\) and fullness \(kL\to\kappa\) satisfy
\[
 0.493255\le w\le 0.0831\cdot6.01,\qquad
 \kappa=5.95.
\]
The shallow record mass \(b\) is at least \(0.0829\) because
each surviving prime occurs at depth at least one.
The first \(T\)-term is in a fixed unit residue class
modulo \(n\).

On that first term group the roots of scaled degree at least
three unweighted at base \(q=\lfloor6/L\rfloor\), with all
remaining factors and exponent corrections in \(v\).
This writes the complete term as \(v z^q\) times its known
exceptional coefficient. The narrow height constant is at
most \(c+\alpha_T/3<0.168\). The limiting numerical correction
is at most \(2\alpha_T\): for nonatomic mass at \(r\le6\)
the cost of assigning the coordinate to zero or six is at
most \(6-r\), and for atomic mass its cost is at most \(r\).
At three the two corrections agree, so changing the side
of a coordinate near the cutoff causes no limiting error.
The degree truncation and zero-window errors are absorbed
by convergence of the full first moments. All frozen
coefficients are units at the retained modulus primes.

Apply (QD2) with \(h=0.168\), \(E=0.004\),
\(\alpha=1/2\), \(M=15\) and \(r=0.371\).
The upper bound for \(B\) is \(0.122328\), while the
determinant threshold is at least
\(0.0829\cdot(1/2)\cdot16\cdot0.371/2=0.1230236\).
The Taylor precision is \(2.968<5.95\), the short
exponent parameter is below the base six, and the
terminal threshold \(4B\) is at most \(0.489312\),
strictly less than \(w\). The partner column is determined
after recovering the first complete \(T\)-term, by pair
uniqueness with its fixed row and the divisor-many choices
of outside representations. Every retained prime has
modulus depth below \(2k\), since the original roots are
squarefree and the degrees of the locked block are at
most \((6.01+o(1))/L\). Hence the full shallow record is
just \(\operatorname{rad}(n)\), of mass at most \(0.0831\).
Including the freeze, the good branch costs at most
\[
 0.0831+0.0831=0.1662<c .
\]
The bad branch also saves a fixed amount. This contradicts
the chosen near-full column exponent and proves (AM1).

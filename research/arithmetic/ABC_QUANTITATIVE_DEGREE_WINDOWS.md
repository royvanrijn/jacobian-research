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
Write \(w+o(1)=\log_X n\) with \(w>0\). Choose a positive integer
\(k\) with \(kL\to\kappa>0\), write \(e_p=v_p(n)\), and put
\[
 Q_n=\prod_{p\mid n}p^{\lfloor e_p/k\rfloor}.
\]
One may use \(Q=Q_n\) or any divisor \(Q\mid Q_n\) as
the residue modulus. A prime of \(n\) with \(e_p<k\) contributes
only to the full modulus and terminal reconstruction. Write \(\log_X Q=(b_Q+o(1))L\)
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

At a balanced endpoint with \(d=1/3\) and \(\Gamma=1/2\), let
\[
 \alpha_T=\int r\,d\nu_{T,\mathrm{at}}(r)
\]
be the atomic numerical weight of one letter. Then
\[
 \boxed{\ \alpha_T>1/250\quad\text{for each letter }T.\ }
 \tag{IA1}
\]
This uses [the single-label complete realization and support bound (MS11)](ABC_DIFFUSE_FACTOR_GRAPH.md#complete-values-determined-by-single-outside-labels), and [the quantitative grouped-root count](#a-quantitative-congruence-estimate-for-one-grouped-root) with a shallow residue choice. The count applies to arbitrary atomic lists and degree scales.

At this endpoint each letter has normalized root mass \(1/3\), nonatomic mass \(c=1/6\), numerical height one, and no escaped first moment. By (MS11), \(\nu_{T,\mathrm{na}}\) is supported on \([0,6]\). Consequently its centered deficit is exactly
\[
 \int(6-r)\,d\nu_{T,\mathrm{na}}(r)=\alpha_T .
 \tag{IA2}
\]
If \(\alpha_T\le1/250\), the nonatomic mass in the open band \((119/20,601/100)\) is at least
\[
 c-20\alpha_T \ge \frac16-\frac2{25}=\frac{13}{150}.
 \tag{IA3}
\]
Indeed a nonatomic root at or below \(119/20\) has scaled-degree deficit at least \(1/20\). This band has more mass than the \(33/400\) used below.

### Complete rows and a high band modulus

Use the complete-label graph and take two fixed distinct complete rows on one of the opposite letters, with common outside columns of exponent arbitrarily close to \(c\). The other free complete value is a function of its column label. Thus on the two rows the complete \(T\)-values, denoted \(T_0,T_1\), satisfy
\[
 T_1-T_0=D
\]
for a fixed nonzero integer \(D\) of polynomial height. The same relation holds when the target is originally one of the addends, with \(D\) changing sign if needed. Complete values determine their prescribed labels up to divisor-many representations at fixed truncation, so this does not affect the column exponent.

Let \(\xi=1/10000\). Discard columns in which the outside \(T_1\)-product shares radical of size at least \(X^{\xi L}\) with \(D\). The [fixed-integer divisor allocation bound](ABC_DIFFUSE_FACTOR_GRAPH.md#counting-shared-prime-divisors) gives only \(X^{(c+\delta-\xi+o(1))L}\) such outside labels, where \(\delta\) is the outside-window loss. Their fixed row recovers the corresponding columns up to divisor-many choices. Take all graph construction and window losses sufficiently small compared to \(\xi\), so the remaining near-full columns suffice for a contradiction.

Choose a compact continuity sub-band of \((119/20,601/100)\) of nonatomic mass greater than \(33/400\), possible by (IA3). On the \(T_1\)-letter the outside coefficients retain this diffuse mass. Consecutively choose whole outside roots in that band until their mass is at least \(33/400\) and at most \(827/10000\), allowing vanishing errors. The upper bound follows by making the outside-coordinate mass cap sufficiently small. Freeze these roots. Strip the full powers of primes dividing \(D\) from them. This does not add a choice cost since \(D\) is fixed. Bounded-scale outside roots are squarefree and have disjoint prime support. Hence the unit modulus remaining has weight \(w\) and shallow residue factor \(Q\) satisfying
\[
 \begin{split}
 w&\ge \frac{33}{400}\frac{119}{20}
            -\frac{1}{10000}\frac{601}{100}
      >\frac{49}{100},\\
 b:=\frac{\log_X Q}{L}+o(1)&\ge \frac{33}{400}-\frac1{10000}
      =\frac{103}{1250}.
 \end{split}
\]
For \(Q\) use the product of the unpowered stripped roots, so also \(b\le827/10000\). The congruence is \(T_0\equiv-D\pmod n\) with unit target. The fullness parameter can have \(kL\to\kappa=59/10\), since each retained degree is larger than \(119/20\).

### Grouping the complete opposite target

On \(T_0\) group all nonatomic root factors at the common base exponent \(q=\lfloor6/L\rfloor\). Put the entire atomic factors in a wider rational cofactor together with the signed degree corrections of the grouped nonatomic roots. To realize this on the boxes take shrinking windows at finitely many positive atoms and zero, then increase the finite list slowly so as to capture the full atomic mass and moment. Such a diagonal is possible because there is no escaped numerical weight. The remaining root mass tends to \(c\), and any corrections from unlisted atoms and their degree rounding tend to zero. Thus the grouped narrow root has normalized height \(h\le1/6\).

The wide numerical exponent is at most \(\alpha_T\). Indeed the selected atoms contribute total positive weight \(\alpha_T\) in the limit. By the support condition and (IA2), the diffuse corrections about degree six are nonpositive in the limiting measure and have absolute numerical weight \(\alpha_T\). Reduction bounds the larger of the numerator and denominator exponents by \(\alpha_T+o(1)\). One may instead use finite atom lists with arbitrarily small omitted mass and moment and absorb that error in the fixed determinant margins.

Apply (QD2) in its first-coordinate reconstruction form, with shallow record \(Q\), half-rank, \(M=15\), and \(r=37/100\). Under \(\alpha_T\le1/250\),
\[
 \begin{split}
 B&\le\frac16\frac{37}{100}+15\frac1{250}=\frac{73}{600},\\
 b(M+1)r/4&\ge\frac{3811}{31250}>\frac{73}{600},\\
 (M+1)r/2&=\frac{74}{25}<\frac{59}{10},\qquad r<6,\\
 4B&\le\frac{73}{150}<\frac{49}{100}<w.
 \end{split}
\]
The last condition suffices because the local records first determine the short root, and the original full modulus recovers the wide cofactor. The actual chart price is \(b\), since the residue record is \(Q\). Each exact grouped output gives the full \(T_0\), then the complete column value follows from row zero's fixed additive value. The original column label has only divisor-many possibilities for that complete value. Thus the good-column exponent is at most
\[
 2\cdot\frac{827}{10000}=\frac{827}{5000}<\frac16.
\]
The discarded overlap columns have their own strict saving. This contradicts the complete rectangle with near-full column exponent and proves (IA1).

### A formal profile that this rules out

Give the other two letters each an atom of mass \(1/6\) at degree three plus nonatomic mass \(1/6\) uniform on \([2,4]\). Give \(T\) an atom of mass \(1/6\) at \(3/250\), and nonatomic mass \(1/6\) uniform on \([6-3/125,6]\). All numerical moments are one, \(\Gamma=1/2\), all degrees are at most six, and the measures have no escaped moment. The two-letter fractional packs at mass \(1/6\) have weight at most one and the three-factor packs at mass \(1/3\) have weight at most two. The trimmed support bound (MS12) holds since \(\min(r,s)\ge sr/6\) for \(r\le6\) and \(s\le6\). Its atomic numerical weight is \(\alpha_T=1/500\), so (IA1) excludes the profile.

### Persistence of the bound near minimum mass

The same individual lower bound \(1/250\) persists for
all sufficiently small \(\gamma=\Gamma-1/2\) at \(d=1/3\).
Indeed the [uniform diffuse tail estimate](ABC_DIFFUSE_FACTOR_GRAPH.md#diffuse-tails-and-mesh-discrepancy-near-minimum-mass)
gives \(\int(r-6)_+\,d\nu_{T,\mathrm{na}}=o(1)\) as
\(\gamma\to0\), and the height and escaped moment bounds
give \(h_T=1+o(1)\) and \(e_T=o(1)\). The nonatomic
mass is \(c+O(\gamma)\). If \(\alpha_T\le1/250\),
the nonatomic mass below \(5.95\) is at most
\(20\alpha_T+o(1)\), by integrating the deficit below
six. The mass above \(6.01\) is \(o(1)\), so the
lock band still has mass at least \(13/150-o(1)\).

As in the preceding grouping, the narrow height is at
most \(c+\gamma\). The atomic numerical weight enters
the numerator of the wide factor, and the positive
nonatomic correction above six is \(o(1)\). The
negative nonatomic correction below six is at most
\(\alpha_T+o(1)\), by subtracting its moment from
\(6c_T\). The escaped weight can add at most \(o(1)\).
This gives wider height at most \(1/250+o(1)\).
For example \(h\le1/6+10^{-5}\) and \(E\le0.00401\)
still leave all the displayed determinant and terminal
margins strict. The near-minimum complete-label graph
allows a two-row rectangle with exponent as close
to \(c\) as needed for the \(0.1654\) good-column cost
and the fixed overlap loss. The continued presence of
these margins proves the claim.

### A complementary bound from an opposite high support

At \(\Gamma=1/2\), let \(s_*\) be the larger of the
two opposite nonatomic support endpoints. The complete-label
truncated-moment restriction (MS12) and the same residue
parameters give
\[
  \int r\,d\nu_{T,\mathrm{na}}(r)+\frac{6-s_*}{6}
       > \frac1{250}.
 \tag{IA4}
\]
In the only nontrivial case \(s_*>5.976\), choose \(s<s_*\)
close enough to \(s_*\) and use \(T\) as the common complete
column parameter on many fixed rows. Lock blocks on the
opposite high-support letter, in a fixed band above \(s\),
with total root mass between \(0.0825\) and \(0.0827\).
There are enough rows for the positive nonatomic mass in
this band, however small. Use the fixed shift-difference
product as in the [single-letter mesh count](ABC_DIFFUSE_FACTOR_GRAPH.md#an-explicit-one-letter-mesh-discrepancy)
to discard high-overlap columns and remove in total at
most \(0.0001\) root mass from the selected locks.
The resulting CRT modulus is full to \(kL\to5.9\),
has weight greater than \(0.49\), and has an unpowered
shallow record of mass between \(0.0824\) and \(0.0827\).

Put \(B_T=\int r\,d\nu_{T,\mathrm{na}}(r)\).
By (MS12), \(\int(s-r)_+\,d\nu_{T,\mathrm{at}}(r)\le B_T\),
since the atomic mass is \(c\) and the nonatomic clipped
weight is at most \(B_T\). Group the atomic factors on \(T\)
unweighted at base \(s\) and leave its nonatomic factors in
the wider cofactor. The negative correction is at most \(B_T\),
and its signed total is \(1-s/6\). The larger wide height
is therefore at most \(B_T+(6-s)/6\), and the narrow height
is at most \(c\). If (IA4) fails, this is at most
\(1/250+o(1)\) as \(s\uparrow s_*\). The strict margins in
the preceding one-root calculation apply, and its cost
\(0.1654<c\) contradicts the selected near-full column
count. This proves (IA4).

## A single-letter bound on concentration in short degree intervals

Work at \(d=1/3,\ \Gamma=1/2\) and fix one target letter \(T\). Put \(c=1/6\). Use the unequal-root-height count (UW) with coefficient \(B_m=m^2+m\). Suppose its positive atomic scales can be grouped using \(p\ge1\) fixed positive bases. This includes a list of at most \(p\) positive atoms, using one base per atom. More generally it suffices that the scales belong to the additive semigroup generated by positive bases \(a_1,\ldots,a_p\), with nonnegative integer multipliers. Then
\[
 \nu_{T,\mathrm{na}}(I)\ll |I|^{1/p}
 \tag{IF1}
\]
for positive compact intervals \(I\), uniformly in the patch position for fixed bases and profile. In particular, the nonatomic measure on a single letter whose atoms lie in one lattice has a globally bounded density, irrespective of the opposite scales.

To prove this, take two fixed distinct complete rows in the [complete-label graph](ABC_DIFFUSE_FACTOR_GRAPH.md#complete-values-determined-by-single-outside-labels), with common columns of exponent arbitrarily near \(c\). Their complete targets \(T_0,T_1\) satisfy \(T_1-T_0=D\) for a fixed nonzero integer \(D\) of polynomial height; the same holds in each signed orientation. Fix a compact positive continuity range \([b,B]\) of positive diffuse mass on \(T\). Use the outside coordinates of \(T_1\) there for a modulus, whose cost will be matched to the patch.

Let \(I\) be a patch with nonatomic mass \(t>0\) and width \(\sigma\). On \(T_0\) retain zero and positive atom windows, together with the patch, omitting the shrinking atom windows from the patch coordinates to avoid double counting. At \(\Gamma=1/2\) the whole spectra have bounded support, so a positive semigroup has only finitely many scales in the relevant range. Alternatively one can take finite lists with omitted atomic mass and moment as small as needed and pass to the limiting estimate. The freeze cost is at most \(c-t\), up to these arbitrarily small errors. The \(p\) grouped atomic variables have bounded normalized root heights depending on the bases, because a multiplier at base \(a_i\) is bounded by \(r/a_i\) and the total numerical weight is at most one. The patch is one additional unweighted grouped variable of normalized root height \(t\), at an exponent represented by its center. Thus the wider numerical exponent is at most \(\sigma t+o(1)\).

Set \(m=p+1\) and choose numerical lock weight
\[
 w=\frac{bt}{100(B_m+1)}.
\]
For small \(t\), the outside \(T_1\) coordinates of \([b,B]\) can be accumulated to this weight with vanishing truncation error, since they retain positive nonatomic mass and their individual normalized masses tend to zero when the window losses tend to zero. The lock costs at most \(w/b+o(1)\). Discard columns sharing outside \(T_1\)-radical at least \(X^{\xi L}\) with \(D\), using the fixed-integer divisor allocation bound. Choose \(\xi<w/(4B)\), then all outside and rectangle losses sufficiently small relative to \(\xi,t\). This discarded family has its own strict exponent saving. Strip primes dividing \(D\) from the locked roots. With \(w'\in[3w/4,w+o(1)]\) their weight, the target congruence is \(T_0\equiv-D\pmod n\), with a unit right-hand side. Use fullness parameter \(\kappa=b/2\).

In (UW), the \(p\) normalized narrow height caps relative to the modulus are \(O(1/t)\), and the patch cap is \(O(1)\). Their product is \(O(t^{-p})\), with constants independent of the patch. The relative wider exponent is \(O(\sigma)\). Thus (UW) applies if \(\sigma\le c_0t^p\) for a sufficiently small fixed \(c_0>0\). Exact monomial fibres are divisor-many even if an original root appears in several of the atomic groups: the exact value determines \(T_0\) up to its fixed coefficients, and fixed truncation gives divisor-many original representations. Amplification is unnecessary. The chart cost is at most \(B_mw'/\kappa\le t/50+o(1)\), while the lock cost is at most \(t/100+o(1)\). With the freeze \(c-t\), this gives a strict saving from the near-full column exponent. Therefore \(t\ll\sigma^{1/p}\). Continuity endpoints suffice to prove the bound for all intervals by enlargement. The constants in the height bounds do not depend on the patch position, and all positive patch centers are allowed. Intervals approaching zero use continuity of the nonatomic part there, and the bounded support handles large intervals.

### Folded distributions on shifted meshes

Fix \(a\in(0,1]\), a phase \(\phi\in[-1/2,1/2]\), and windows with \(\|r/a-\phi\|\le\theta\) capturing nonatomic mass \(t>0\) on \(T\). Select finite positive atoms omitting atomic mass at most \(t/16\), with unshifted weighted phase error at most \(\theta t/2\). Group these atoms on the unshifted mesh and the windows on the shifted mesh, omitting the selected atom neighborhoods from the shifted part. Use a first narrow variable of height \(X^{(O(1/a))L}\) and a second of height \(X^{(17t/16+o(1))L}\) carrying the offset exponent \(\lfloor a\phi/L\rfloor\). The wider exponent is at most \(2a\theta t+o(1)\). For these windows take the row modulus at numerical
weight \(w=\epsilon_T t\), with \(\epsilon_T>0\) fixed
and small enough that \(w\) is at most half the
available nonatomic numerical weight of \([b,B]\)
for every \(t\le c\), and that the lock and chart costs
are at most \(t/100\) and \(t/50\), respectively.
This keeps the choice of modulus uniform even for large
captured mass \(t\). The narrow cap product is
\(O(1/(a t))\) and the relative wider exponent is
\(O(a\theta)\). Thus \(t\ll\theta\), uniformly in the phase and mesh, whenever the atomic alignment condition holds. For a lattice spectrum it holds on every refinement of that lattice, so all folded nonatomic measures on those meshes have a common bounded density.

## Partial degree grouping on complete columns

Let \(d=1/3,\ \Gamma=1/2\), and fix a letter \(U\). Choose any prescribed collection of its original root coordinates with limiting root mass \(m\in(0,1/3]\) and weak degree measure \(\sigma\). For example, take all coordinates in a union of fixed continuity windows. Work with the full prelimit group and leave every other factor of \(U\) in the wide cofactor. At base exponent \(q=\lfloor1/(mL)\rfloor\), the narrow root has normalized height at most \(m\). The net signed numerical correction tends to \(1-qLm=0\), so the numerical wider height exponent is
\[
 E= \int (1/m-r)_+\,d\sigma(r)
    =1-\int\min(r,1/m)\,d\sigma(r).
 \tag{PG1}
\]
Indeed only the selected factors contribute negative corrections, while the complement contributes positively. The signed net vanishes, so both sides have the same limiting exponent. The first integrand is bounded and continuous, so this identity holds without using convergence of the first moment of the selected coordinates.

Let \(s_*\) be the maximum nonatomic support endpoint of the other two letters. Then for every integer \(M\ge2\),
\[
 \boxed{\quad
 1-\int\min(r,1/m)\,d\sigma(r)
 \ \ge\ \frac{s_*}{M}
          \left(\frac{829}{40000}-\frac{m}{M+1}\right).
 \quad}
 \tag{PG2}
\]
Only positive right sides are substantive.

To see this, choose an opposite target band of positive nonatomic mass in a compact interval \([s,S]\) with \(S<2s\). Use many fixed complete rows, taking enough to choose their target outside blocks with aggregate root mass between \(0.083\) and \(0.0831\). [The one-letter mesh argument](ABC_DIFFUSE_FACTOR_GRAPH.md#an-explicit-one-letter-mesh-discrepancy) supplies these locks. Strip primes dividing the fixed product of all row differences after discarding the corresponding high-overlap columns. Choose the overlap threshold so that their total root loss is at most \(0.0001\). The remaining target blocks are pairwise coprime, and their congruences specify one unit residue for the complete \(U\)-value. Their full product has numerical weight at least \(0.0829s\), and choosing \(kL\to s\) gives the unpowered shallow record of mass \(b_Q\in[0.0829,0.0831]\).

The grouped congruence now has narrow height \(h\le m\) and numerical wider height \(E\). In (QD2) use half-rank, the given \(M\), and \(r=s/(M+1)\). Since \(s\le6\) by the diffuse support bound and \(m\le1/3\), one has \(r<1/m\) for \(M\ge2\). The Taylor precision is \(s/2<kL\). If
\[
 E<\frac{s}{M}\left(\frac{0.0829}{4}-\frac{m}{M+1}\right),
\]
then \(B=hr+EM<0.0829s/4\), which is strictly below the determinant threshold and gives \(4B<w\). Thus the short coordinate is reconstructed first, and the full modulus recovers the wide factor. Each output determines the complete column value and its associated outside label up to divisor-many choices. The lock and chart costs total at most \(0.1662<1/6\), contradicting the near-full complete rectangle. Let \(s\) approach \(s_*\) to obtain (PG2).

### An explicit excluded two-component profile

Take each opposite letter to have mass \(1/6\) at degree three and nonatomic mass \(1/6\) uniform on \([2,4]\). On \(U\), take
- mass \(1/6\) at \(499/100\);
- nonatomic mass \(1/30\) uniform on \([997/200,999/200]\), centered at \(499/100\);
- nonatomic mass \(2/15\) uniform on \([0,3/100]\).

Every letter has root mass \(1/3\), nonatomic mass \(1/6\), and numerical moment one. All supports are within \([0,6]\), so the support cap, the two- and three-letter packing bounds at mass \(1/6\), and the three-factor packs at mass \(1/3\) hold. The inequality \(\min(r,s)\ge rs/6\) for \(r,s\in[0,6]\) gives (MS12). The nonatomic and atomic numerical moments are bounded away from zero on all letters.

This profile also satisfies the one-letter unshifted phase inequality (DN2). Indeed for \(U\) the low diffuse component alone contributes at least \(1/500\) to \(\int\|r/a\|\,d\nu_U\) uniformly in \(0<a\le1\). For \(a\ge3/50\) its contribution is exactly \(1/(500a)\); otherwise the phase average on its interval is at least \(1/8\). For the other two letters their uniform component on \([2,4]\) supplies more than this lower bound. The required opposite support constants in (DN2) are at most \(6/11500\).

But take \(\sigma\) to be precisely the high part on \(U\), namely its atom and its high diffuse band. Its mass is \(m=1/5\), all its support is below the grouped base \(5\), and its moment is \(499/500\). Thus \(E=1/500\) by (PG1). Use \(s=39/10\), a target band contained in \([39/10,4]\), and \(M=20,\ r=39/210\). Then
\[
 B=\frac15\frac{39}{210}+\frac{20}{500}
    < .0829\frac{39/10}{4}.
\]
The bound (PG2) excludes this formal profile despite the full unshifted phase, support, and packing tests.

### A count saving in a nearby sector

For an explicit open-margin sector, change the low diffuse support of \(U\) to \([0,3/200]\) and the common center of its high atom and high diffuse band to \(999/200\), keeping the high diffuse width \(1/100\). Then all moments still equal one. At base five the high group has \(m=1/5\) and wider height \(E=1/1000\). The following nearby spectral conditions suffice for an explicit count saving. Assume
\[
 |c_W-1/6|\le 1/8000000\quad\hbox{on all letters}.
\]
On one opposite target \(T\), require nonatomic mass greater than \(3/400\) in \([39/10,4]\). On \(U\), require a prescribed group at base five having narrow normalized height at most \(201/1000\) and wider numerical exponent at most \(11/10000\), with all other factors in the wide cofactor. Then
\[
 d<\frac13-\frac1{1000000}.
 \tag{PG3}
\]

Use the
[quantitative complete-label reduction below \(d=1/3\)](ABC_DIFFUSE_FACTOR_GRAPH.md#three-letter-divisor-packing)
with \(d=2c-\epsilon\) and \(\gamma_0=1/8000000\).
Since \(\sum_W c_W-3c\le3\gamma_0\), the complete graph has
edge exponent at least \(2c-35\epsilon-68\gamma_0\) up to
arbitrarily small losses. The global limiting upper bound
gives \(\epsilon\ge0\).
Assume \(\epsilon\le1/1000000\). Take ten fixed complete rows with common column exponent at least
\[
 c-350\epsilon-699\gamma_0-o(1)>c-.000438-o(1).
\]
This follows from the ordered-star count on the complete-label graph with capacities \(c+\gamma_0\), using the conservative \(68\gamma_0\) edge loss. Choose one outside band block per target row of mass \(3/400\) up to the arbitrarily small rounding loss and freeze them at total cost at most \(.0751\). Strip primes dividing the product of all fixed row differences, discarding columns whose radical overlap on any row has root mass at least \(\xi=.0005\). Their exponent is at most \(c+\gamma_0-\xi+o(1)\), below the selected column exponent. The combined root-mass loss is at most \(10\xi=.005\), so the pairwise coprime target modulus gives \(b_Q\in[.07,.0751]\) and \(w\ge.075\cdot3.9-.005\cdot4=.2725\). Use \(\kappa=3.9\).

Apply the half-rank determinant with \(M=16,r=1/5\) on the specified partial \(U\) group. Then \(B\le .201/5+16(.0011)=.0578\), below \(b_Q(17)r/4\ge .0595\). The Taylor precision is \(1.7<\kappa\), \(r<5\), and the retained modulus has positive weight since \(4B=.2312<w\). The actual chart cost is at most .0751, giving total count exponent at most .1502, a contradiction to the common columns. This proves (PG3).

## Selective phase quantiles on complete columns

Work at \(d=1/3,\ \Gamma=1/2\), with \(c=1/6\). Fix a target letter \(U\), and write
\[
 s_*=\max_{T\ne U}\sup\operatorname{supp}\nu_{T,\rm na}.
\]
For every \(a\in(0,1]\) and Borel set \(J\subseteq[0,\infty)\) with \(\nu_U(J)=c+g>c\),
\[
 \boxed{\displaystyle
       \int_J \|r/a\|\,d\nu_U(r)\ \ge\ \frac{s_*g^2}{310}.}
 \tag{SP1}
\]
This tests the discrepancy of retained roots alone, with the remaining roots frozen at a cost below the column exponent. The balance \(g\) is the difference between the retained mass and \(c\).

At a fixed spacing, let \(A_U(a)=\int\|r/a\|\,d\nu_{U,\rm at}(r)\). For every Borel set \(I\) of nonatomic mass \(t>0\), (SP1) gives
\[
 A_U(a)+\int_I\|r/a\|\,d\nu_{U,\rm na}(r)
     \ \ge\ \frac{s_*t^2}{310} .
 \tag{SP2}
\]
Indeed remove the countable set of atoms from \(I\) and add all the atoms back as a separate disjoint piece. More generally, if only a prescribed collection of whole atoms is included, of mass \(c-p\), the same inequality holds with its partial atomic phase discrepancy and with \((t-p)^2\) on the right when \(t>p\). Thus no finite atomic list hypothesis is needed.

In particular
\[
 \nu_{U,\rm na}\{r:\|r/a\|\le\theta\}
 \ \le\ \frac{310\,\theta}{s_*}
          +\sqrt{\frac{310A_U(a)}{s_*}} .
 \tag{SP3}
\]
The same upper bound with \(\theta=\sigma/a\) applies to the mass below scaled degree \(\sigma>0\). Indeed \(\|r/a\|\le\sigma/a\) there. These consequences use \(t^2\le ut+v\Rightarrow t\le u+\sqrt v\). If all positive atoms lie on lattices \(a_j\mathbb Z\) apart from atomic mass \(p_j\), the trimmed form directly gives
\[
 \nu_{U,\rm na}\{r:\|r/a_j\|\le\theta\}
 \ \le\ 2p_j+\frac{1240\theta}{s_*},
\]
since \(t>2p_j\) implies \((t-p_j)^2>t^2/4\). The full discrepancy consequence at \(g=c\) is \(s_*/11160\), while (SP1) also permits excluding arbitrary outlying degree bands.

### Proof of the residual estimate

It suffices to use a finite union \(J\) of bounded intervals with continuity endpoints for \(\nu_U\). Every Borel set can be approximated in \(\nu_U\)-symmetric difference by these unions. To check this, enclose a compact subset inside an open superset whose measure difference is arbitrarily small, then cover the compact set by finitely many inside intervals with endpoints avoiding the countable set of atom locations. The same approximation controls the integral since the phase function is bounded by \(1/2\).

Choose an opposite letter \(T\) and a compact band \([s,S]\) on it of positive nonatomic mass, with \(0<s<S<2s\) and continuity endpoints. By (MS11), \(s\le6\). Use the [complete-label graph](ABC_DIFFUSE_FACTOR_GRAPH.md#complete-values-determined-by-single-outside-labels) to take any fixed number \(m\) of complete rows on the remaining letter, with common \(U\)-columns of exponent arbitrarily close to \(c\). They have distinct fixed complete row values. On the \(i\)-th row, the complete \(T\)-value has form \(\epsilon u+t_i\), where \(u\) is the common complete \(U\)-value, \(\epsilon=\pm1\) is fixed, and the fixed offsets \(t_i\) are distinct.

Select original outside root positions on these rows from the supplying \(T\) band, of total normalized mass \(0.496g\) up to arbitrarily small rounding losses. This can be arranged by prescribing \(0.496g/m\) on each row once \(m\) is large enough, and taking the outside-window errors small enough relative to \(m,g\). Freeze their values at normalized cost at most \(0.497g\). Put \(D=\prod_{i<j}(t_i-t_j)\). Discard the columns on which any outside \(T_i\)-radical has overlap of mass at least \(\xi\) with \(D\), where \(0<\xi<0.0005g/m\) is fixed. The [divisor allocation argument](ABC_DIFFUSE_FACTOR_GRAPH.md#counting-shared-prime-divisors) and row uniqueness bound those columns with a strict saving from \(c\), after taking the graph and extraction losses small enough. Strip primes of \(D\) from the frozen blocks, losing root mass less than \(m\xi\) up to vanishing errors. The resulting full blocks are pairwise coprime because any common prime divides the corresponding shift difference. Their product \(n\) gives a single fixed residue for \(u\bmod n\), which is a unit by primitivity. For \(kL\to s\), each retained valuation is in \([k,2k)\) since the band roots are squarefree and \(S<2s\). The shallow record \(Q_n=\operatorname{rad}(n)\) has mass \(b_Q\le0.497g\), with lower bound \(b=0.495g\), and the modulus weight satisfies \(w\ge0.495gs\). The selected roots were counted before stripping, so no new record charge is required.

On \(U\) freeze the root coordinates outside \(J\), at cost \(c-g+o(1)\) since the total root mass is \(2c\). Any fixed exceptional coefficient can be absorbed into the frozen product at arbitrarily small cost. Every factor of the frozen product is a unit modulo \(n\). For a retained original degree \(r=eL\), assign the nearest nonnegative multiple \(ia\). Take the \(i\)-th power of its original root into the integer narrow root \(z\), and use the common exponent \(q=\lfloor a/L\rfloor\). Put the correction into the wide rational multiplier \(v\). Then the unit congruence is \(v z^q\equiv K\bmod n\), and the normalized narrow and numerical wide heights satisfy
\[
 h\le \frac1a+\frac16\le\frac7{6a},
 \qquad
 E\le a\int_J\|r/a\|\,d\nu_U(r).
\]
For the first estimate use \(i\le r/a+1/2\), total mass \(1/3\), and first moment at most one. For the second, the signed rounding errors have height at most the sum of their absolute weights. The errors from flooring \(a/L\) contribute at most \(L\sum i\beta_e=o(1)\) for fixed \(a\). Weak convergence controls the remaining integral on the continuity set. All errors can be absorbed into the strict margins below.

Suppose the discrepancy is at most \(sg^2/310\). In the [one-root determinant (QD2)](ABC_QUANTITATIVE_DEGREE_WINDOWS.md#a-quantitative-congruence-estimate-for-one-grouped-root) take rank fraction \(\alpha=1/2\), \(M=\lceil19/(ag)\rceil\), and \(r_0=asg/19\). As \(s\le6\) and \(g\le c\), one has \(r_0\le a/19<qL+o(1)\). The Taylor order obeys
\[
 \frac{(M+1)r_0}{2}
 \le \frac{s}{2}\left(1+\frac{2ag}{19}\right)
 \le\frac{29s}{57}<s.
\]
The determinant and reconstruction parameters have the strict margins
\[
 \begin{aligned}
 B=hr_0+EM
   &\le sg\left(\frac7{114}+\frac{19+1/6}{310}\right)
      =\frac{871}{7068}\,sg
      <\frac{99}{800}\,sg,\\
 \frac{b(M+1)r_0}{4}
   &\ge \frac{0.495}{4}\,sg=\frac{99}{800}\,sg,
 \qquad 4B<0.495sg\le w .
 \end{aligned}
\]
Thus the normalized chart exponent is at most \(b_Q\le0.497g\). Each resulting pair \((z,v)\) determines \(u\) given the fixed roots, hence the column up to divisor-many representations at fixed truncation. Total normalized counting cost is at most
\[
 c-g+0.497g+0.497g=c-0.006g.
\]
The discarded bad-overlap columns have their own strict saving. Take the input star and extraction losses smaller than both margins. This contradicts the near-full column exponent. We have proved strict discrepancy \(>sg^2/310\) for every such band. Let \(s\) tend to \(s_*\) to get (SP1).



### An infinite-atom example with a low diffuse component

Let \(V,T\) each have \(c\delta_3+c\,\mathrm{Unif}[2,4]\). Take \(\epsilon=10^{-6}\), and on \(U\) use
\[
 \begin{split}
  \nu_{U,\rm at}&=\sum_{j\ge1}c\,2^{-j}\,
                 \delta_{9/2+\epsilon2^{-j}},\\
  \nu_{U,\rm na}&=\tfrac1{20}\,\mathrm{Unif}[0,\epsilon]
      +\tfrac7{60}\,\mathrm{Unif}\bigl[1,\,23/7-(29/21)\epsilon\bigr].
 \end{split}
\]
The uniform laws each have mass one before multiplying by their coefficients. Every letter has mass \(1/3\), nonatomic mass \(c\), and first moment one. On \(U\) the atomic shift contributes \(\epsilon/18\) to the moment and the low diffuse component adds \(\epsilon/40\), both offset by the decrease \(29\epsilon/360\) in the higher diffuse component.

All degrees are below \(4.501\). Since \(s_*=4\) across the letters, the left side of (MXP1) is at most \(4c+(4.501-4)c<1\). The fractional three-letter pack of mass \(c\) has weight below one, and the pack at mass \(1/3\) below two. The diffuse support bound (MS11) holds. For (MS12), the \(U\)-atomic mass is \(c\) above degree four, and the opposite profiles lie in \([0,4]\) with moment one, so each gives the required compressed moment for \(s\le4\). Every block of mass \(\beta<c/m\) also has \(mw<1\), so (RP1) supplies no common factor. The full phase inequality (DN2) has a margin on all three letters: each higher diffuse uniform interval has length at least two, so its phase average for every \(0<a\le1\) is at least \(1/8\), by integrating the complete periods. Every opposite support endpoint is at most four.

The partial unweighted grouping test (PG2) also holds. On \(U\), for selected mass \(m\le1/4\) the largest possible captured ordinary moment is at most \(3/4+\epsilon/18+121/588\): select all the atoms and at most \(1/12\) of the top mass of the high interval. The average on its upper fraction \(5/7\) is at most \(121/49\). For \(m\ge1/4\), the clipping base \(1/m\le4\) loses at least \(1/12\) from the atomic moment alone. Thus the unweighted residual exceeds \(0.044\), while the right side of (PG2) is at most \(4(0.0829/8)=0.04145\). On \(V,T\), for \(m\le0.23\) the residual is at least \(1-4m\ge0.08\). For \(0.23<m\le0.3\) the omitted mass is at least \(1/30\) at degree at least two, giving residual at least \(1/15\). For \(m>0.3\), clipping at \(10/3\) on the diffuse band gives residual at least \(1/54\). When \(m>0.23\) the right side of (PG2) is negative for \(2\le M\le10\) and at most \(4(0.0829/4)/11<0.008\) for \(M\ge11\).

Nevertheless at \(a=1/2\), the \(U\)-atom list and the low diffuse interval have joint phase discrepancy \(29\epsilon/180\). This follows by summing \(2\epsilon\sum_{j\ge1}c4^{-j}=\epsilon/9\) from the atoms and \(\epsilon/20\) from the low interval. With \(t=1/20\) and \(s_*=4\), the right side of (SP2) is \(1/31000\), a contradiction. Infinitely many distinct atomic scales occur in a fixed compact interval, so the finite-positive-basis hypothesis in (IF1) cannot be used.

### A nearby count saving with one supplying row

Consider a limiting profile of total root mass \(1/3\) on each letter, with \(|c_W-c|\le10^{-7}\) on each letter. Suppose an opposite letter \(T\) has nonatomic mass \(>0.025\) in \([3,4]\). On \(U\) suppose some continuity set \(J\) of root mass at least \(c+0.0499\) satisfies
\[
  \int_J\|2r\|\,d\nu_U(r)\le10^{-6}.
\]
Then \(d<1/3-10^{-6}\).

For if \(d=2c-\varepsilon\) with \(0\le\varepsilon\le10^{-6}\), the [near-minimum complete-label reduction](ABC_DIFFUSE_FACTOR_GRAPH.md#the-complete-label-reduction-close-to-minimum-mass), using the three caps \(c+10^{-7}\), retains edge exponent at least \(2c-35\varepsilon-68\cdot10^{-7}\) up to arbitrarily small losses. This is the same bound used in (MS16), since \(\gamma=\sum c_W-3c\le3\cdot10^{-7}\). Equivalently the simultaneous anchors cost at most \(7\varepsilon+12\cdot10^{-7}\), and the vector factor dictionary step adds at most four times that deficit plus \(8\cdot10^{-7}\). Fix one complete row on the remaining letter, by averaging. Its common \(U\)-columns have exponent at least \(c-35\varepsilon-69\cdot10^{-7}>c-0.000042\), up to arbitrarily small losses.

Take \(g=0.0499\), \(s=3\), and \(a=1/2\) in the preceding count. One supplying row has enough \(T\)-outside roots in \([3,4]\) for root cost \(0.496g<0.025\). There is no row-overlap deletion, and the \(U\)-residue is a unit by primitivity. Use lower and upper record masses \(0.495g\) and \(0.497g\) as above. The freeze outside \(J\) costs at most \(c-g\), while \(10^{-6}<sg^2/310\) supplies the determinant condition. The resulting bound \(c-0.006g=c-0.0002994\) contradicts the fixed-row column count. The displayed hypotheses all hold for the infinite-atom example, using a continuity enlargement of its low interval and atomic cluster.

## Pooled phase quantiles for unrestricted limiting profiles

Use the limiting normalized abc-box notation, with normalized count exponent \(d\le1/3\), scaled degree measures \(\nu_W\), total root masses \(l_A+l_B+l_C\le1\), and numerical heights at most one. The nonatomic masses and the numerical heights need not be equal. There is no complete-label or minimum-\(\Gamma\) hypothesis.

Fix a free pair \(U,V\), its complementary modulus letter \(T\), a spacing \(0<a\le1\), and Borel retained degree sets \(J_U,J_V\). Let
\[
 f=\nu_U(J_U^c)+\nu_V(J_V^c),\qquad g=d-f.
\]
For \(r\in J_W\), assign a nearest nonnegative grid point \(i(r)a\). At a half-integer tie either assignment can be prescribed for that atom. Put \(\sigma_U=+1,\ \sigma_V=-1\), and define
\[
 D_\pm=\sum_{W=U,V}
  \int_{J_W}\Big(\pm\sigma_W\bigl(r/a-i(r)\bigr)\Big)_+\,d\nu_W(r),
 \qquad D=\max(D_+,D_-).
\]
The signs account for the rational height of the correction factor in the \(U/V\) ratio. In particular \(D\le\sum_W\int_{J_W}\|r/a\|\,d\nu_W\).

Suppose \(g>0\) and a compact positive band \([s,S]\), \(S<2s\), has nonatomic \(T\)-mass \(H>0\). Put \(u=\min(g,2H)\). Then
\[
 \boxed{\displaystyle D\ge \frac{s u^2}{400}.}
 \tag{UP1}
\]
If \(l_U,l_V\le1/3\), replace \(400\) by \(310\). For continuity bands separated from their boundary degrees, the comparison can be made strict for each fixed band; the weak form accommodates endpoints and passage to a supremum. The same conclusion holds for any chosen \(0<u\le g\) if \([s,S]\) provides an attainable lock of root mass \(0.496u+o(u)\) through prescribed whole degrees or edge-dependent divisors. A nonatomic band permits this whenever \(u\le2H\).

### Proof of the pooled phase inequality

Approximate the retained sets by finite unions of bounded intervals with continuity endpoints in \(\nu_W\)-measure. The discrepancy is bounded, so this approximation changes the masses and integrals by arbitrarily little. To prescribe a tie assignment, use a shrinking neighborhood at each of the finitely many retained atomic ties in such bounded windows; on that neighborhood use the prescribed adjacent grid point. The limiting signed phase then has the stated integral, since other tie locations are continuity points or have total omitted mass tending to zero. Bands containing boundary atoms can be reduced to smaller continuity bands with nonatomic mass tending to \(H\) and with degree bounds tending to \(s,S\).

Freeze a block of \(T\)-roots in the band with total normalized mass \(0.496u\), up to small rounding losses, always taking its record cost at most \(0.497u\). This is possible since \(0.496u\le0.992H\). To see the nonatomic selection directly in the original boxes, exclude finitely many shrinking atom windows so that every remaining original coordinate is arbitrarily small, while retaining the nonatomic band mass. Bounded-scale original anatomic roots are squarefree and have disjoint supports. The full divisor \(n\) obtained by raising the selected roots to their original degrees is \(k\)-full with \(kL\to s\), and its valuations belong to \([k,2k)\) for small \(L\). Its shallow record \(Q_n=\operatorname{rad}(n)\) has normalized logarithmic size \(b_Q\le0.497u\), with lower bound \(b=0.495u\), and its numerical weight is \(w\ge0.495us\). For the edge-dependent divisor version, fix a dyadic size list first, at subpower cost at fixed truncation. Counting the selected divisors by their sizes gives the same record exponent. The congruence and count are uniform in the divisor choices.

Freeze the \(U,V\) factors outside the retained sets, at cost \(f+o(1)\). Up to this fixed ratio, group the retained factors at exponent \(q=\lfloor a/L\rfloor\), putting \(x_e^{\sigma_W i(eL)}\) in the narrow rational \(z\) and \(x_e^{\sigma_W(e-q i(eL))}\) in the wider rational \(v\). In the ratio congruence \(v z^q\equiv K\bmod n\), the target and all variables are units by primitivity and the signed additive relation. The normalized narrow height and numerical wider height satisfy, up to vanishing errors,
\[
 h\le\max(h_U/a+l_U/2,\ h_V/a+l_V/2)\le\frac{3}{2a},
 \qquad E\le aD.
\]
For the narrow bound, estimate numerator and denominator independently with \(i\le r/a+1/2\). The wider bound uses the independent numerator and denominator sums \(D_+,D_-\). Flooring \(a/L\) contributes at most \(L\sum i\beta_e=o(1)\) by the prelimit total first-moment bound. The bounded phase discrepancy together with the prescribed ties gives the limit by weak convergence. Even a first moment that escapes at infinity causes no problem, since the narrow estimate uses the prelimit first moments and the signed rounding errors are bounded by the phase distance.

If \(D\le su^2/400\), apply (QD2) with rank fraction \(1/2\),
\[
 M=\left\lfloor\frac{24}{a u}\right\rfloor,\qquad r_0=\frac{s}{M+1}.
\]
One has \(a u M\le24<a u(M+1)\). The supplying mass gives \(su\le2\) in the nonatomic case and \(su\le1/0.495\) for adaptive locks, since their numerical weight is at least \(0.495su\). Therefore \(r_0/a<su/24<1\) with room for the floor loss in \(qL\). The Taylor order is exactly \((M+1)r_0/2=s/2<s\). The determinant height parameter satisfies
\[
 B=hr_0+EM\le su\left(\frac1{16}+\frac{24}{400}\right)
 =\frac{49}{400}\,su
 <\frac{99}{800}\,su
 \le\frac{b(M+1)r_0}{4}.
\]
Also \(4B<0.495su\le w\). The normalized chart count is at most \(0.497u\). With the fixed roots, each chart pair determines the actual ratio of \(U,V\) and therefore their two coprime complete values; their relation recovers the third value. The total normalized cost is at most
\[
 f+0.497u+0.497u=f+0.994u< f+g=d,
\]
a contradiction. The slack is at least \(0.006u\). When \(l_U,l_V\le1/3\), use \(h\le7/(6a)\) and replace \(24\) by \(19\). With \(D\le su^2/310\) the height bound becomes \(B/(su)\le 7/114+19/310<99/800\), and the Taylor order remains \(s/2<s\). This proves the sharpened constant.

### A sector away from the minimum nonatomic total

At a formal \(d=1/3\) profile, let \(\nu_T=\frac13\operatorname{Unif}[2,4]\). On both \(U,V\) take
\[
 \nu_W=\tfrac16\delta_3+\tfrac1{20}\operatorname{Unif}[3-\sigma,3+\sigma]
                       +\tfrac7{60}\operatorname{Unif}[2,4],
 \qquad \sigma=10^{-5}.
\]
All masses are \(1/3\), first moments are one, and \(\Gamma=2/3\). The full phase discrepancy of each free measure has a fixed contribution from the broad interval, while the partial sets \(J_U=J_V=[3-\sigma,3+\sigma]\) separate out the low-dispersion component. Take \(a=1\), \([s,S]=[3,4]\), \(H=1/6\). Then \(g\ge1/10\), so \(u\ge1/10\). By the identical free measures and their symmetry about three,
\[
 D=\frac{\sigma}{40}+\frac7{120}\sigma^2
      < 3\cdot10^{-7},
 \qquad
 \frac{s u^2}{310}\ge \frac{3}{31000}.
\]
Thus the profile violates (UP1). In fact any limiting sequence with normalized count \(d\ge 1/3-10^{-3}\), \(l_U,l_V\le1/3\), \(f\le7/30\), \(H\ge1/10\) on this supplying band, and \(D\le10^{-5}\) on the specified retained sets is excluded: \(g\ge0.099\) and the same determinant margins apply. One may use the \(400\) bound and permit \(l_U,l_V\) above \(1/3\) as well. This supplies a fixed coefficient saving in the sector without a near-minimum nonatomic assumption.

## Combining two valuation bands in a full modulus

Use the notation of the pooled phase inequality (UP1). Retain arbitrary prescribed degree sets on a free pair \(U,V\). Let \(f\) be the frozen complement root cost, and write the grouped ratio, up to its fixed multiplier, as \(v z^q\), with \(qL\to a>0\), normalized narrow height \(h\), and numerical wider height \(E\). Select any prescribed blocks on the opposite letter, with root cost \(\mu\), whose full product \(n\) has numerical weight \(w\) and is \(k\)-full for \(kL\to\kappa>0\). Select any divisor \(Q\) of the shallow record \(\prod_{p\mid n}p^{\lfloor v_p(n)/k\rfloor}\), of normalized weight \(b\). The primes used for the shallow record may be a proper subset of the full modulus primes, with lower powers than in the full modulus.

For any \(0<\alpha<1\) and integer \(M\ge1\), put \(D=b\alpha(M+1)/2-h\). If \(D>0\), a sufficient condition for the count exponent bound \(d\le f+\mu+b\) is
\[
 E<
 \frac{D}{M}
 \min\left(a,\ \frac{\kappa}{\alpha(M+1)},\
                   \frac{(1-\alpha)w}{b\alpha(M+1)}\right).
 \tag{MB1}
\]
For the proof, choose \(r\) just above \(EM/D\). The first two upper bounds on \(EM/D\) give \(r<a\) and \(\alpha(M+1)r<\kappa\). The last gives \(hr+EM<(1-\alpha)w/2\) for \(r\) sufficiently close to its lower bound, since at \(r=EM/D\) the left side is \(EM(D+h)/D\). The determinant inequality holds strictly by the choice \(r>EM/D\). Apply (QD2) to the unit-target congruence from the signed relation. Each pair determines the exact coprime terms, and hence the triple, after restoring the frozen factors. The count is \(f+\mu+b\). The statement needs no complete-label graph or hypothesis on \(\Gamma\).

A necessary upper bound on the numerical wider height for any parameter choice in (QD2), for a fixed \(b,w,h\) with \(2h>b\), is
\[
 E<\frac{bw}{32(2h-b)}.
 \tag{MB2}
\]
Indeed its determinant conditions force \(x=b\alpha(M+1)/(2h)>1\) and \(E<(1-\alpha)w(1-1/x)/(2M)\). Also \(1/M\le b\alpha/((2h-b)x)\) because \(x>1\) and \(\alpha<1\). Use \(\alpha(1-\alpha)\le1/4\) and \((x-1)/x^2\le1/4\) to obtain (MB2). This is an upper limit for the parameters in (QD2).

### A phase sector using two bands

Suppose \(a=1\), \(h\le0.656\), \(E\le0.00029\), and \(f\le7/30\). On the modulus letter suppose there is positive nonatomic mass greater than \(0.0332\) in \([3,3.05]\), and a degree-500 atomic window whose root mass is \(0.0002\). This sector has
\[
 d<\frac13-10^{-3}.
 \tag{MB3}
\]
The point values of the atomic mass and degree can be replaced by sufficiently small open tolerances without losing the strict estimates below.

Freeze mass \(0.033+o(1)\) from the nonatomic low band and the entire high atomic window. Freeze cost is at most \(0.0331+0.00021=0.03331\). The full numerical weight is at least \(0.1985\), since the low blocks give \(3(0.033)=0.099\) in the limit, and the high block contributes \(0.1\). Use \(kL\to3\). In the shallow record include the unpowered low radical and the 160th power of the high radical. Indeed every high exponent is \((500+o(1))/L\), and its shallow valuation cutoff is at least 160. Squarefreeness and disjoint original degree positions give record size between \(0.0649\) and \(0.0651\), allowing the usual small losses. Take \(\alpha=1/2\), \(M=80\), \(r=0.0375\). Then
\[
 B\le 0.656(0.0375)+80(0.00029)=0.0478
  <\frac{0.0649\cdot81\cdot0.0375}{4}=0.0492834375,
\]
\[
 (M+1)r/2=1.51875<3,\qquad 4B=0.1912<0.1985 .
\]
The count from (QD2) thus has normalized exponent at most
\[
 f+0.03331+0.0651\le\frac7{30}+0.09841
                      <\frac13-10^{-3}.
\]
Only one modulus letter is needed. The higher block supplies numerical weight and part of the shallow record, and the lower band provides the remainder of the shallow record.

For an exact formal example, let \(\sigma=7/625\). Give both free letters the measure
\[
 \nu_U=\nu_V=\tfrac16\delta_3
          +\tfrac1{20}\operatorname{Unif}[3-\sigma,3+\sigma]
          +\tfrac7{60}\operatorname{Unif}[2,4].
\]
On the modulus letter use
\[
 \nu_T=\tfrac1{5000}\delta_{500}
       +\tfrac1{20}\operatorname{Unif}[3,61/20]
       +\tfrac{4247}{15000}\operatorname{Unif}[r_*-1/10,r_*+1/10],
       \qquad r_*=\frac{44925}{16988}.
\]
Every measure has mass \(1/3\) and first moment one. For the modulus letter the high and middle contributions are \(1/10\) and \(121/800\), with the remaining moment \(599/800\). Its nonatomic total is \(1/3-1/5000\), so the combined nonatomic total is \(2/3-1/5000\), away from \(1/2\).

Retain \(J_U=J_V=[3-\sigma,3+\sigma]\) on the mesh \(a=1\). The grouped root and the two signed remainder weights have exact limiting heights
\[
 h=\frac{13}{20}+\frac7{20}\sigma=0.65392,\qquad
 E=\frac{\sigma}{40}+\frac7{120}\sigma^2
       \in(0.000287,\,0.00029).
\]
The retained excess above the total cost \(f\) at \(d=1/3\) satisfies \(g=1/10+(7/30)\sigma<0.103\). Thus all hypotheses of (MB3) hold. For this fixed retained grouping, either low-only or high-only locks fail even the necessary one-root criterion (MB2), allowing any valid shallow record. Indeed, with \(d>f+\mu+b\), one must have \(\mu+b<g<0.103\). If all modulus primes selected are from the lower scales, \(w\le3.05\mu\), so
\[
 \frac{bw}{32(2h-b)}
 <\frac{3.05(0.103)^2}{4\cdot32(1.3-0.103)}
 <0.000212.
\]
If only the high atom is selected, \(w\le0.1\), so the same expression is at most \(0.1(0.103)/(32(1.3-0.103))<0.000269\). Both are below \(E\). These comparisons are for the prescribed retained grouping.

## The symmetric diffuse profile and the one-group determinant

Consider the formal minimum-nonatomic profile on all three letters
\[
 \nu=\tfrac16\delta_3+\tfrac16\operatorname{Unif}[2,4],
 \qquad d=1/3.
\]
The numerical heights and first moments are all one. For the count based on one grouped power in the free ratio, original-root enumeration for frozen factors, and the determinant conditions (QD2), no choice of retained submeasures, fixed scaled mesh spacing, multiplier assignments, full modulus, shallow residue record, or fixed amplification gives a strict saving below \(1/3\).

The frozen modulus factors are charged their original root mass. Both the retained free factors and the frozen modulus factors may come from fractional selections of the limiting measures.

### Necessary numerical conditions for a successful count

Allow arbitrary retained submeasures on the two free letters, including fractional retention of the mass at degree three. The combined retained mass must be \(1/3+g\) for some \(0<g\le1/3\), since the original frozen complement costs \(1/3-g\). Let \(t\in[g,1/3]\) be the retained nonatomic mass. The retained atomic mass is \(1/3+g-t\).

Let \(qL\to a>0\) be the common power. Group retained roots with arbitrary signed integer multipliers into a rational \(z\), and put the signed residual exponents into \(v\), so that \(v z^q\) is the exact retained ratio. Write \(h\) for the limiting normalized rational height of \(z\) and \(E\) for the limiting numerical rational height of \(v\). The exact retained ratio has numerical height \(m\) satisfying
\[
 m\le ah+E,\qquad
 m\ge\frac{1+3g-t+3t^2}{2}\ge\frac{11}{24}+\frac{3g}{2}.
 \tag{SB1}
\]
The first inequality is the rational height bound on the product. For the second, the retained atomic moment is \(1+3g-3t\), and the least diffuse moment at mass \(t\) is \(2t+3t^2\) because the pooled diffuse density is \(1/6\) on \([2,4]\). The exact numerator and denominator are coprime, so their larger numerical height is at least half their combined numerical weight. To get the last bound, minimize \(-t+3t^2\) at \(t=1/6\).

The numerical wider height \(E\) is at least half the sum of the absolute exponent corrections, since primes in the two original letters and degree coordinates are disjoint. Each correction is at least the distance to a nearest grid point, even if another integer multiplier is chosen. Consequently:
\[
 \begin{array}{ll}
 0<a\le4:&
 E\ge \dfrac{3a t^2}{4(2+a)},\\[4pt]
 4<a<6:& E\ge \dfrac32 t^2,\\[4pt]
 a\ge6:& E\ge \dfrac13+g.
 \end{array}
 \tag{SB2}
\]
For the first bound, the pooled diffuse density of \(r/a\) modulo one is at most \((2+a)/6\). Indeed each phase has at most \(2/a+1\) preimages in \([2,4]\), each contributing \(a/6\). Folding to nearest-integer distance at most doubles the density, so the integral of this distance over retained diffuse mass \(t\) is at least \(t^2/(2(2+a)/3)\). Take half and multiply by \(a\). The same bound holds for any retained submeasure by domination by that density. For \(4<a<6\) the distance in the original degree units from \(r\in[2,4]\) to \(a\mathbb Z\) is at least \(4-r\). Its minimum integral at retained diffuse mass \(t\) is \(3t^2\), by taking the top interval. For \(a\ge6\), every retained degree in the support is at distance at least two from the grid. These arguments pass to limiting submeasures by bounded continuous distance to the fixed grid on compact ranges, with first-moment tightness at this profile.

Let \(\mu\) be the frozen modulus-root mass, \(w\) its numerical weight, and \(b\) the selected shallow record weight. Allow any fixed amplification order \(\ell\ge1\), and put \(b'=b/\ell\) for the effective chart exponent. Saving would require
\[
 \mu+b'<g,\qquad w\le4\mu,\qquad E<w/2\le2\mu.
 \tag{SB3}
\]
The weight bound uses support below or at four and the absence of escaped moment, so also applies to arbitrarily selected divisors of original roots. The terminal condition of (QD2) gives \(E<w/2\), including after amplification.

If \(2h>b'\), the necessary bound (MB2) after replacing \(h,E,b\) by \(\ell h,\ell E,\ell b'\) reads
\[
 E<\frac{b'w}{32\ell(2h-b')}
 \le\frac{g^2}{32(2h-b')}.
 \tag{SB4}
\]
The last comparison uses \(b'w\le4b'\mu\le(b'+\mu)^2<g^2\). If \(2h\le b'\), all upper bounds on \(2h-b'\) used below are already automatic.

### Incompatibility

If \(0<a\le4\), (SB2) and (SB4), or \(2h\le b'\) directly, give
\[
 ah\le \frac a2 b'+\frac{2+a}{48}.
\]
Therefore
\[
 ah+E<2g+\left(\frac a2-2\right)b'+\frac{2+a}{48}
          \le2g+\frac18
          <\frac{11}{24}+\frac{3g}{2},
\]
since \(g\le1/3\). This contradicts (SB1).

If \(4<a<6\), (SB2) and (SB4) give \(2h-b'\le1/48\), and \(E<2\mu\) gives \(\mu>3t^2/4\). Thus \(b'<g-3t^2/4\), and
\[
 ah+E<2g+\left(\frac a2-2\right)b'+\frac{a}{96}
       \le 3g-\frac34 t^2+\frac1{16}.
\]
Subtract the first lower bound on \(m\) in (SB1). The difference is at most
\[
 \frac32 g+\frac12t-\frac94 t^2-\frac7{16}
   \le 2t-\frac94t^2-\frac7{16}\le-\frac1{48},
\]
because \(g\le t\le1/3\). This is again impossible.

If \(a\ge6\), (SB2) gives \(E\ge1/3+g\ge2g\), contrary to \(E<2\mu<2g\). This proves the obstruction for every fixed positive spacing, even using the full range of mixed full/shallow modulus criteria (MB1).

## Fixed differences with a shallow divisor from the rows

Work at the minimum nonatomic endpoint \(d=1/3,\Gamma=1/2\), with \(c=1/6\), using the complete-label graph and its all-scale progression estimate (AC1). Give the letters distinct roles \(A\) (rows), \(B\) (columns), and \(T\) (shifted targets).

The following three conditions exclude the endpoint.

- \(A\) has nonatomic root mass strictly greater than \(0.14\) in a fixed compact band of scaled degrees at least \(4\).
- The supremum of the nonatomic support on \(B\) is at least \(4\).
- The complete \(T\)-measure is supported in \([2.95,3.05]\) and satisfies \(\int|r-3|\,d\nu_T(r)\le0.004\).

In particular this excludes \(\nu_T=c\delta_3+c\,\operatorname{Unif}[2.952,3.048]\) with the prescribed opposite profiles. The fixed row difference supplies the full reconstruction weight. Its selected shallow divisor supplies the Taylor divisibility.

### A determinant using a shallow divisor

The usual one-group determinant proof applies to an arbitrary positive modulus \(n\) if its shallow record \(Q\mid n\) obeys
\[
 v_p(Q)\le\lfloor v_p(n)/k\rfloor\qquad(p\mid Q).
 \tag{FDS1}
\]
Keep all the other hypotheses of (QD2), including the unit condition and the weight \(w=\log_X n\) for terminal reconstruction. This is the only place global \(k\)-fullness is used in the one-group proof. In detail, fix a residue of \(z\bmod Q\) and lift \(v\) at each \(p\mid Q\) to \(Kz^{-q}\). The Taylor determinant uses these primes, and its order \(N<k\) makes \(v_p(n)\ge N v_p(Q)\) on them. Primes outside \(Q\) impose no determinant divisibility requirement. The resulting polynomial vanishes at the original rational sample points. At every prime of the full \(n\), including primes not in \(Q\), the polynomial equation is exact and \(vz^q=K\) holds modulo \(p^{v_p(n)}\). The logarithmic Jacobian, Hensel lift, and stripping and reconstruction arguments therefore apply over all of \(n\). Their numerical budget is still \(2B/(1-\alpha)<w\), because only the log mass in the Jacobian numerator is discarded. The finite-lift and valuation record count remains divisor-many at fixed truncation. In particular there is still just one initial chart charge, \(\log_X Q/L\).

### Selecting fixed rows with a wide difference

Choose a complete-label graph with edge exponent arbitrarily close to \(2c\), and enough regularization to use (AC1) on its complete \(A\)-vertices, with opposite degree bound \(s=3.999\) supplied by \(B\). Fix original outside positions on \(A\) from its indicated band, of aggregate normalized root mass \(p=0.14+o(1)\), where the rounding and outside-window losses can be made arbitrarily small. The nonatomic hypothesis allows this selection with uniformly small original coordinates. The full factor \(F\) from those positions has numerical weight at least \(0.56-o(1)\). Pigeonhole its original root values at cost \(p+o(1)\).

The retained fibre graph has edge exponent at least \(2c-0.14\) up to arbitrarily small losses. Its rows have exponent at most \(c-0.14\), up to the corresponding outside-window error, because these particular outside coordinates have been frozen. Its columns have exponent at most \(c\) up to the initial graph loss. By the ordered two-star count and subtraction of the diagonal, there are at least
\[
 X^{(3c-0.28-o(1))L}
\]
ordered pairs of distinct rows with a common column, up to arbitrarily small graph and extraction losses. The diagonal has a strict saving since the row fibre exponent \(c-0.14\) is positive.

Rows in this fibre are multiples of \(F\). For any fixed one of them, rows within absolute distance \(X^{0.666}\) lie in an interval of length \(2X^{0.666}\) in that same progression. Once the losses are small enough, its effective length exponent for (AC1) is at most \(0.1061\), up to a vanishing error, since \(F\ge X^{0.5599}\). The count of such rows has exponent at most \(0.1061/3.999\), up to arbitrary losses, which is strictly less than \(c-0.14\). Consequently the stars with a close pair have a strict saving from \(3c-0.28\): count a first row, its close second row by (AC1), and any common column by the full column cap. Average the remaining stars over at most \(X^{(2(c-0.14)+o(1))L}\) row pairs, with small graph losses. Thus there are two fixed distinct complete rows with common column exponent arbitrarily close to \(c\) and with numerical difference
\[
 |D|>X^{0.666}.
 \tag{FDS2}
\]
On those rows the two complete targets satisfy \(T_2-T_1=D\), with the fixed signs absorbed into the definition of \(D\). The fixed nonzero \(D\) and either target value determine the complete column value.

### Removing shared prime powers and counting

For each of the two targets let \(S_i\) be the full part supported on primes dividing \(D\). Freeze \(S_1,S_2\) at subpower cost. The bound is uniform for \(0<|D|\le2X\). For fixed \(0<\eta<1\), the number of integers up to \(X\) supported on primes of \(D\) is bounded by \(X^\eta\prod_{p\mid D}(1-p^{-\eta})^{-1}=X^{\eta+o(1)}\). The last equality follows by ordering those primes and using \(\omega(D)=O(\log X)\). Let \(G=\gcd(T_1,T_2)=\gcd(S_1,S_2)\), and denote \(W=\log_X G\). Put \(n=|D|/G\). Since \(T_2/G-T_1/G=D/G\) and the two reduced terms are coprime, their ratio is a unit and is \(1\bmod n\). For \(p\mid n\), the two original \(T_i\) have equal \(p\)-adic valuation: if their valuations were unequal, \(D\) would have their minimum valuation. Hence \(S_2/S_1\) is also a unit modulo \(n\). The residual complete ratio \((T_2/S_2)/(T_1/S_1)\) therefore has a fixed unit target modulo \(n\). By primitivity on both rows, \(F\) is coprime to both targets and divides \(D\), so it divides \(n\).

Strip every prime dividing \(D\) from the original \(T_i\) factors and group all remaining roots unweighted at exponent \(q=\lfloor3/L\rfloor\). Let \(z\) be the quotient of the two residual radical products, and put the signed degree corrections into \(v\). All variables are units modulo \(n\) after stripping. The residual exact ratio is \(vz^q\), and the normalized narrow height and numerical wider height satisfy
\[
 h\le\frac13-\frac{W}{3.05}+o(1),\qquad E\le0.004+o(1).
 \tag{FDS3}
\]
For the first bound, each removed supported part \(S_i\) has numerical weight at least \(W\), giving removed root mass at least \((W-o(1))/3.05\) by the limiting support and the absence of escaped weight at \(\Gamma=1/2\). This remains valid uniformly if a vanishing root subset lies outside the limiting band, since its numerical moment is \(o(1)\). For the second bound, the complete \(T\)-mass and first moment give \(\int(r-3)_+\,d\nu_T=\int(3-r)_+\,d\nu_T\). The ratio correction on either sign side is at most the sum of these two integrals. Removing prime factors decreases each of the individual \(T_i\) correction sums. Floors in \(q\) contribute a total \(o(1)\) by the root-mass bound. Any negligible exceptional factor can be absorbed into the wider term. Fixing \(S_i\) already accounts for the exact removed integer parts. If degree assignments for their primes are recorded as well, there are divisor-many choices at fixed truncation.

Take \(Q=\operatorname{rad}F\), \(kL\to4\), \(b\ge0.1399\), \(b_Q\le0.1401\), \(\alpha=1/2\), \(M=20\), and \(r_0=0.2\). The full primes of \(F\) guarantee (FDS1), since this block lies in original outside positions of scaled degree at least four up to vanishing errors. The determinant height satisfies
\[
 B=h r_0+ME\le\frac{11}{75}-\frac{4}{61}W+o(1)
               \le\frac{11}{75}+o(1)
               < \frac{0.1399(21)(0.2)}{4}.
\]
The Taylor order is \(21(0.2)/2=2.1<4\), and \(r_0<qL+o(1)\). For terminal reconstruction,
\[
 w=\log_X n\ge\max(0.5599,\ 0.666-W).
\]
The difference between this lower bound and \(4B\) is at least
\[
 \max(0.5599,\ 0.666-W)-\frac{44}{75}+\frac{16}{61}W-o(1).
\]
The piecewise linear expression decreases until \(W=0.1061\) and increases thereafter. Its minimum is \(0.5599-44/75+(16/61)(0.1061)>0.001\). Thus the estimate is uniform in the entire possible \(\gcd\)-weight range, including a large shared atomic factor. By the relative-to-shallow (QD2) proof, the chart exponent is at most \(0.1401\). There is no positive column lock charge for \(F\), which was paid in the row fibre before averaging. The two stripped parts are fixed at subpower cost. Each resulting ratio determines \(T_1,T_2\) via \(D\ne0\), then the column. This contradicts the common-column exponent approaching \(c\).


### A profile excluded by the fixed-row determinant

Take
\[
 \nu_A=c\delta_{7/4}+c\,\operatorname{Unif}[4,4.5],\quad
 \nu_B=c\delta_3+c\,\operatorname{Unif}[2,4],\quad
 \nu_T=c\delta_3+c\,\operatorname{Unif}[2.952,3.048].
\]
Each measure has mass \(1/3\), moment one, and nonatomic mass
\(c\). The \(T\)-integral of \(|r-3|\) is exactly \(0.004\), and
both opposite support conditions hold. The mixed trimmed packing
test (MXP1) permits all its locks, because the highest scaled degree
is \(4.5\) and the available normalized lock mass is \(c\), giving
\(4.5c<1\).
[Mixed trimmed packing](ABC_DIFFUSE_FACTOR_GRAPH.md#a-trimmed-packing-bound-using-mixed-additive-moments)

## A primitive-kernel improvement to the one-root determinant

Use the rational-unit congruence \(v z^q\equiv K\bmod n\), with heights \(H(z)\le X^{(h+o(1))L}\), \(H(v)\le X^{E+o(1)}\), and \(\log_X n=w+o(1)\). Let \(Q\mid n\) be the shallow record, with \(\log_X Q=(b_Q+o(1))L\) and \(b_Q\ge b>0\). At all primes in \(Q\), assume \(v_p(Q)\le\lfloor v_p(n)/k\rfloor\), where \(kL\to\kappa>0\). Other primes of \(n\) need no fullness assumption. Use the notation and limiting order of (QD2).

Choose \(M\ge1\), \(r_0>0\), and \(0<\alpha<1\), and put
\[
 B=hr_0+EM,\qquad \Lambda=\frac{b(M+1)r_0}{2}.
\]
The usual one-root conclusion, at most \(X^{(b_Q+o(1))L}\) pairs up to an arbitrary fixed-truncation exponent loss, holds under the strict conditions
\[
 r_0<\liminf qL,\qquad
 \alpha(M+1)r_0<\kappa,\qquad
 B<\alpha\Lambda,\qquad
 4\Lambda\bigl(1-\sqrt{1-B/\Lambda}\bigr)<w .
 \tag{PK1}
\]
For any lower bound \(w_0\le w\) with \(0<w_0\le4\Lambda\), the last condition follows from
\[
 B<\frac{w_0}{2}-\frac{w_0^2}{16\Lambda}.
 \tag{PK2}
\]
These use divisibility of the maximal *nonzero* evaluation minors to improve the auxiliary coefficient height. The determinant forcing and the terminal residue count otherwise use the (QD2) argument.

### Divisibility at the actual rank

In one fixed \(z\bmod Q\) class take the monomial column set \(0\le i\le R=\lfloor r_0/L\rfloor,\ 0\le j\le M\), of cardinality \(D_0=(R+1)(M+1)\). At a rational pair \((z,v)\) multiply its evaluation row by \(\operatorname{den}(z)^R\operatorname{den}(v)^M\), using reduced positive denominators. The resulting rows are integral, with entry sizes at most \(X^{B+o(1)}\). The row factors may differ from point to point. Each is a unit at every modulus prime, so row clearing does not alter the prime-adic valuation bounds on any minor.

Set \(N=\lfloor\alpha D_0\rfloor\). At a prime \(p\mid Q\), abbreviate \(d_p=v_p(Q)\). Replace \(v\) by \(Kz^{-q}\) in the local reference evaluations, with error divisible by \(p^{v_p(n)}\). Taylor expansion in \(z-z_0\) inside the fixed residue class shows that for every \(1\le r\le N\) an \(r\)-row minor is divisible by \(p^{d_p r(r-1)/2}\). The Taylor vectors of orders \(0,\ldots,r-1\) give the smallest total order, since a repetition makes the corresponding determinant contribution zero. Every replacement error is divisible by \(p^{d_p r}\), since \(r\le N<k\) for small \(L\) by the strict precision condition. The same rank-\(r\) divisibility therefore holds for the original integral evaluation rows. At primes not in \(Q\) there is no Taylor requirement. Multiplying over \(p\mid Q\) gives \(Q^{r(r-1)/2}\mid\det\) for every such minor.

At \(r=N\) the divisibility exponent per row tends to \(\alpha\Lambda\), so \(B<\alpha\Lambda\) makes every \(N\)-row minor vanish. In a nonempty chart let \(r<N\) be the actual rank. Choose \(r\) independent rows and form their integer evaluation matrix \(A\) with \(D_0\) columns. Every \(r\)-row maximal minor is divisible by \(Q^{r(r-1)/2}\), while its absolute value is at most \(X^{r(B+o(1))}\). No assumption about genericity of the evaluation ranks is needed.

### Passing from maximal minors to the auxiliary height

For any rank-\(r\) integer matrix \(A\), let \(g_A\) be the positive gcd of its maximal minors. The covolume of its primitive integer kernel is
\[
 \det(\ker A\cap\mathbb Z^{D_0})
       = \frac{\sqrt{\det(AA^\mathsf T)}}{g_A}.
\]
Indeed the wedge of the \(r\) independent integral rows has the maximal minors as coordinates. Dividing by their gcd is the primitive wedge of the saturated row lattice, so its Euclidean norm is the covolume of that lattice. The Hodge dual gives the primitive wedge of the perpendicular integer kernel, with the same norm. Cauchy--Binet bounds the numerator by \(\binom{D_0}{r}^{1/2}\) times the largest minor.

Minkowski's theorem on the integer kernel, of rank \(D_0-r>0\), therefore gives a nonzero integer kernel vector with maximum coefficient at most a constant depending on \(D_0\) times
\[
 X^{\frac{r}{D_0-r}
       (B-bL(r-1)/2+o(1))}.
\]
This vector is in the kernel of the entire chart evaluation matrix, since the selected rows span its row space over \(\mathbb Q\). Put \(x=r/D_0\) and note \(bLD_0/2=\Lambda+o(1)\) as \(L\to0\). The exponent is at most
\[
 \frac{x(B-\Lambda x)}{1-x}+o(1).
\]
This estimate is uniform in actual rank \(r\), since \(x\le\alpha<1\). On \(0\le x<1\), one has
\[
 \sup_x\frac{x(B-\Lambda x)}{1-x}
       =\bigl(\sqrt{\Lambda}-\sqrt{\Lambda-B}\bigr)^2
       =G_* .
\]
For example, putting \(y=1-x\) gives the concave maximum of \((B-\Lambda)/y+2\Lambda-B-\Lambda y\) at \(y=\sqrt{1-B/\Lambda}\). This gives the nonzero auxiliary equation of coefficient height \(X^{G_*+o(1)}\) on the same rectangular degree support.

Apply the usual factor-height, logarithmic derivative and Hensel arguments to this auxiliary. The numerator of the cleared logarithmic derivative on each nonexceptional solution has height at most \(X^{B+G_*+o(1)}\). The exceptional exact-constant components and the intersections with the logarithmic derivative cost boundedly many points as in (QD2), since \(R<q\). The total modulus weight discarded or lost in the remaining Hensel reconstruction is at most \(2(B+G_*)+o(1)\), now using all primes of \(n\). But
\[
 2(B+G_*)=4\Lambda(1-\sqrt{1-B/\Lambda}).
\]
Under (PK1) a positive numerical modulus weight remains, enough to reconstruct \(z\) of height \(X^{hL+o(L)}\). With \(z\) fixed the original congruence reconstructs \(v\), since (PK1) implies \(w>2B\ge2E\). The original residue classes cost \(X^{(b_Q+o(1))L}\). The Hensel and divisor records have subpower cost at fixed \(L\). This proves (PK1), and elementary rearrangement on \(0<w_0\le4\Lambda\) gives (PK2).

For reference, if \(A_0=b(M+1)/2>h\) and one ignores the upper constraints on \(r_0\), the best terminal loss for fixed \(M\) occurs at
\[
 r_0=\frac{EM/A_0}{2\sqrt{h/A_0}(1-\sqrt{h/A_0})},
 \qquad B+G_*=\frac{EM}{1-\sqrt{h/A_0}}.
\]
This follows by differentiating \(2A_0r_0(1-\sqrt{1-h/A_0-EM/(A_0r_0)})\). The rank, degree and precision conditions in (PK1) must still be checked for this choice.

### A wider fixed-difference sector

Use the [shared-row construction](#fixed-differences-with-a-shallow-divisor-from-the-rows)
and its removal of the \(D\)-supported parts, keeping the two
opposite-letter hypotheses. The impossible target sector is now
\[
 \operatorname{supp}\nu_T\subseteq[2.875,3.125],
 \qquad\int|r-3|\,d\nu_T(r)\le0.0104.
 \tag{PK3}
\]
In particular the target \(\frac16\delta_3+\frac16\operatorname{Unif}[2.8752,3.1248]\) is excluded with the same two opposite profiles. The degree band and discrepancy in (PK3) give, after dividing out the frozen \(D\)-supported parts, \(h\le1/3-W/(25/8)+o(1)\), \(E\le0.0104+o(1)\), where \(W\) is the numerical weight of the two targets' common divisor \(G\). This is the same height argument as (FDS3): the common numerical weight has been removed from each target, affine shifts give \(T_2-T_1=D\), and the equal positive and negative whole-\(T\) errors about degree three bound the ratio correction even after stripping. The frozen degree allocations cost at most subpower, and the unit ratio congruence is modulo \(n=|D|/G\), containing the fixed shared row block \(F\).

Use \(Q=\operatorname{rad}F\), \(b=0.1399\), \(b_Q\le0.1401\), \(kL\to4\), \(M=10,\ r_0=0.29,\ \alpha=0.94\). We have \(qL\to3\), \(0.94(11)(0.29)=2.9986<4\), and
\[
 \Lambda=\frac{0.1399(11)(0.29)}2=0.2231405,\qquad
 B\le\frac{29}{300}+0.104-0.0928W+o(1)<0.94\Lambda .
\]
In (PK2) take \(w_0=\max(0.5599,\,0.666-W)\), which is always at most \(0.666<4\Lambda\) and no larger than the full modulus weight. For \(0\le W\le0.1061\), the margin
\[
 \frac{w_0}{2}-\frac{w_0^2}{16\Lambda}
                 -\left(\frac{29}{300}+0.104-0.0928 W\right)
\]
is a concave quadratic in \(W\), so its minimum on that interval is at an endpoint. Its values at zero and \(0.1061\) are greater than \(0.007\) and \(0.001\), respectively. For \(W\ge0.1061\), \(w_0=0.5599\) is fixed and the margin increases. Thus (PK1) holds throughout the gcd-weight range. The chart exponent is at most \(0.1401<c\), and the ratio together with the fixed difference gives at most one column for each counted pair. This contradicts the near-full common-column exponent.

## Newton breakpoints in the residue determinant

The first one-root residue determinant suffices when every selected outside prime appears with a nonzero multiplier in the grouped root. The unique-minimum rule for prime valuations applies directly to the auxiliary equation. For a fixed auxiliary polynomial, its nonzero coefficient primes and finitely many degree breakpoints account for all the selected outside prime factors.

### Valuation breakpoints of the auxiliary equation

Suppose a nonempty chart has a nonzero equation
\[
 F(z,v)=\sum_{(a,j)\in\mathcal S} c_{a,j} z^a v^j=0,\qquad
 0\le j\le M ,
\]
with integer coefficients of polynomial height in \(X\), uniformly at fixed truncation. In \(\mathcal S\) retain only terms with \(c_{a,j}\ne0\). A prime \(p\) from a selected original squarefree root of degree \(e\) on one term of a grouped ratio is assumed to occur in no other original factor of that ratio. If that degree position uses a nonzero grouped-root multiplier \(i\), its valuations have the form
\[
 (v_p(z),v_p(v))=\sigma(i,e-qi),\qquad \sigma\in\{1,-1\}.
\]
One may allow any fixed finite set \(\mathcal I\) of such nonzero multipliers. The case where the original root contributes the same positive valuation factor to both coordinates gives the same argument after dividing by it.

Let \(J\) be a fixed nonzero integer divisible by all the nonzero coefficients. When \(p\nmid J\), the valuation of each term of \(F\) is \(\sigma(ai+j(e-qi))\). For each fixed \(j\), the minimum is attained at the smallest supported \(a\) when \(\sigma i>0\) and at the largest when \(\sigma i<0\). Let \(a_j\) be that extremum, for each \(j\) that appears. The minimum among the resulting affine functions of \(e\) must be attained twice for the sum to vanish. For every two distinct \(j,j'\) their slopes differ, and their equality forces
\[
 e=qi+i\,\frac{a_{j'}-a_j}{j-j'}.
\]
Thus for each \(i,\sigma\) there are at most \(\binom{M+1}{2}\) exceptional original degree indices \(e\). The count is uniform even though the polynomial and the normalized breakpoints \(eL\) may depend on the box and on the chart. In terms of scaled measures these are a bounded list of degree positions. If desired one can pass to a subsequence for the positions in a fixed compact degree interval. For counting, it suffices to use the original outside-coordinate cap \(\delta+o(1)\) at each position, which bounds their combined mass by \(2|\mathcal I|\binom{M+1}{2}\delta+o(1)\) without requiring any uniformity for the locations. This is the useful way to handle chart-dependent breakpoint locations.

The nonzero coefficients were fixed before this count. Their product is nonzero and of polynomial height at fixed truncation, even when some coefficients of the rectangular support vanish. A form with a single supported term has no solution on unit rational coordinates. If some \(j\)-slices are absent, simply omit them. No assumption of generic coefficients or of a full support is needed.

If primes from two factors can coincide only on a fixed polynomial-height integer \(D\), take \(J=D\prod_{\mathcal S}c_{a,j}\). Then the same conclusion holds at every selected outside prime not dividing \(J\). The other primes are supported on this fixed integer and can be counted by divisor allocation. The bound permits the narrow degree of \(F\) to grow as \(L\to0\) with fixed \(M\) and a fixed number of multiplier classes.

### An exclusion using only the shared row factor

At \(d=1/3,\Gamma=1/2\), put \(c=1/6\). Suppose a row letter \(A\) has nonatomic mass strictly above \(0.164\) in a compact band of scaled degrees at least four. This ensures that a continuity subband strictly above four still contains more than \(0.164\). Let \(T\) be the shifted target letter and \(B\) the column letter. The following target sector is impossible:
\[
 \int|r-3|\,d\nu_T(r)\le 0.0224 .
 \tag{NB1}
\]
The hypothesis on the opposite letters is entirely on the row letter \(A\).

Use the complete-label realization, with part caps arbitrarily close to \(c\) and edge exponent arbitrarily close to \(2c\). As usual one can thin to injective label-to-complete-value maps at fixed-truncation divisor loss. For the row letter freeze a prescribed block of outside roots of normalized mass \(0.164\), up to arbitrarily small rounding errors, in the indicated subband. Count its root records at that cost and select one giving the lower edge exponent \(2c-0.164\) up to arbitrary losses. Its full product \(F_0\) divides every retained row term, and the row count in its fibre is at most \(X^{(c-0.164+o(1))L}\) with arbitrarily small extraction loss. The fibre edge count, column cap and this row cap give two distinct fixed rows with common-column exponent arbitrarily close to \(c\) by ordered two-star averaging. Indeed the star exponent is \(3c-2(0.164)\) up to losses, its diagonal has a strict deficit \(c-0.164>0\), and the row-pair dictionary exponent is \(2(c-0.164)\). All these margins absorb the original graph and window losses.

The two complete targets satisfy \(T_2-T_1=D\ne0\) for the fixed rows, with the sign of \(D\) allowed to vary. By the common factor \(F_0\) and primitivity, \(T_2/T_1\equiv1\bmod F_0\) and the ratio is a unit. The fixed common row factor is the full modulus for the determinant. Take \(kL\to4\) and \(Q=\operatorname{rad}F_0\), with shallow size \(b_Q\in[0.1639,0.1641]\) up to vanishing errors. Its primes are in bounded scaled-degree anatomic positions and are squarefree in the original roots.

Group both targets unweighted at \(q=\lfloor3/L\rfloor\). Let \(z\) be the ratio of their radical-root products, and let \(v\) contain the signed exponent corrections and any negligible exceptional coefficients. All variables are units modulo \(F_0\), and \(v z^q\equiv1\bmod F_0\). The normalized narrow height is at most \(1/3+o(1)\). The numerical wider height is at most \(0.0224+o(1)\). The limiting first moment is one and the total \(T\)-mass is \(1/3\), so the absolute error about degree three decomposes into equal positive and negative parts for each target. A numerator of the ratio correction can contain at most the positive part from one target and the negative part from the other. The denominator is analogous. Rounding \(q\) changes these limits by \(o(1)\), and no moment escapes at the minimum nonatomic endpoint.

Use only the first determinant step of (QD2). For \(M=7\), \(r_0=0.5\), \(\alpha=0.99\), and \(R=\lfloor r_0/L\rfloor\), the rectangular support of degree at most \(R\) in \(z\) and \(M\) in \(v\) gives determinant height parameter
\[
 B\le\tfrac16+7(0.0224)+o(1)
      <0.324522=\frac{0.1639\cdot0.99\cdot8\cdot0.5}{2}.
\]
The Taylor order check is \(0.99\cdot8\cdot0.5=3.96<4\). Therefore in every \(z\bmod Q\) chart the rank drops and there is a nonzero polynomial equation \(F(z,v)=0\) with wide degree at most seven. Its coefficients can be bounded by maximal minors of the integral evaluation matrix after clearing row denominators. These minors have polynomial height in \(X\) at fixed truncation, since their sizes are bounded independently of \(X\). The height bound makes the product of the nonzero coefficients an integer of polynomial height, which is the input to the valuation count.

For each chart take the product of its nonzero coefficients and \(D\) as the integer \(J\). A prime of \(T_1\) not dividing \(D\) cannot divide \(T_2\). Every original nonpacked outside prime of \(T_1\) not dividing \(J\) hence has the breakpoint form above, with \(i=1\) and one fixed sign. There are at most \(\binom82=28\) exceptional original degree positions in the chart. The packed coordinate has vanishing normalized root mass, taking the compression index above all bounded scaled degrees used in the argument. The parts of the outside roots in these positions can be enumerated at normalized cost \(O(\delta)+o(1)\). At every other original outside position, all prime factors divide \(J\). Their allocation to original degree positions is divisor-many at fixed truncation. For example, squarefree root parts are chosen from the prime support of \(J\), and the number of allocation patterns is at most \((K+1)^{\omega(J)}=X^{o(1)}\) for fixed number \(K\) of original positions. This bounds the \(T_1\)-outside label count per chart by \(X^{(O(\delta)+o(1))L}\).

The complete-label map fixes \(T_1\) given its outside vector, and the first fixed row and the signed equation fix the complete column. After the injective thinning this is at most one actual column per such vector. Equivalently a fixed complete column has only divisor-many original representations. The shallow charts cost at most \(0.1641+o(1)\) in the normalized exponent, and \(\delta\) can be taken small after \(M\) is fixed. This is strictly below \(c\), contradicting the selected common-column count, and proves (NB1).

### A profile inequality from the first residue determinant

At \(d=1/3,\Gamma=1/2\), let \(A,T\) be distinct letters and
write \(H_A(\kappa)=\nu_{A,\mathrm{na}}((\kappa,\infty))\).
For every \(a,\kappa>0\) and every positive integer \(M\),
\[
 \int \min_{j\ge1,\ j\in\mathbb Z}
 \left(|r-aj|+\frac{\kappa j}{M(M+1)}\right)\,d\nu_T(r)
 \ \ge\ \frac{\kappa}{2M}\min\{1/6,H_A(\kappa)\}.
 \tag{NB2}
\]

Fix \(p<\min\{1/6,H_A(\kappa)\}\) positive. Take the shared
row block of root mass \(p\) in a compact continuity band
strictly above \(\kappa\). The ordered star argument above
supplies two distinct fixed rows with common-column exponent
arbitrarily close to \(1/6\), using the same fibre count
with \(p\) in place of \(0.164\). Their common factor is the
congruence modulus and its radical is the shallow record.
In particular \(kL\to\kappa\) and the shallow exponent tends
to \(p\).

Let \(i(r)\) be a bounded positive integer assignment on a
finite partition into degree intervals. Use continuity
endpoints and choose the same assignment on both shifted
targets. Group the roots with multipliers \(i\) at exponent
\(q=\lfloor a/L\rfloor\). The limiting narrow and wider
height caps are
\[
 h=\int i(r)\,d\nu_T(r),\qquad
 E=\int |r-ai(r)|\,d\nu_T(r).
\]
For the wider height, each side of the signed ratio uses
one positive and one negative correction from the identical
target profiles. The absence of escaped numerical moment
allows these height bounds to pass to the limit.

The Taylor rank drop applies as soon as
\(hr_0+EM<\alpha p(M+1)r_0/2\) and
\(\alpha(M+1)r_0<\kappa\), with \(0<\alpha<1\).
The Taylor coefficients are integral for every integer
monomial exponent, so these conditions permit any \(r_0>0\)
that fits the minor precision. The valuation-breakpoint
count then leaves only subpower-many outside labels per
shallow residue class in the limiting normalized exponent.
Its chart cost is \(p<1/6\). Therefore any strict violation
of \(E+\kappa h/(M(M+1))\ge \kappa p/(2M)\) would exclude the
box, by taking \(\alpha\) close to one and \(r_0\) just below
\(\kappa/(M+1)\).

On the compact limiting degree support the pointwise
minimizers in (NB2) range over finitely many integers.
There are finitely many transitions in the minimum of the
piecewise-linear objectives. At a tie the competing
objectives agree, so shifting the transition away from atoms
costs arbitrarily little in the combined integral. Thus
the permitted assignments approach the displayed minimum.
Finally take \(p\) up to \(\min\{1/6,H_A(\kappa)\}\). If the
latter is zero, the bound follows from nonnegativity.

### A profile with a broader target window

Take
\[
 \nu_A=c\delta_{7/4}+c\,\operatorname{Unif}[4,4.5],\quad
 \nu_B=c\delta_3+c\,\operatorname{Unif}[2,4],\quad
 \nu_T=c\delta_3+c\,\operatorname{Unif}[2.7312,3.2688].
\]
All masses are \(1/3\), all moments are one, and all nonatomic masses are \(c\). The target discrepancy is \(0.2688/12=0.0224\), and the row hypothesis holds. The mixed
trimmed packing bound permits every lock of mass at most \(c\)
because the maximum scaled degree is \(4.5\).
[Mixed trimmed packing](ABC_DIFFUSE_FACTOR_GRAPH.md#a-trimmed-packing-bound-using-mixed-additive-moments)


Under the same opposite profile hypotheses, (NB2) also
excludes \(T=(1/6)\delta_3+(1/6)\operatorname{Unif}[2.72,3.28]\).
With \(a=3,\kappa=4,M=7\), the minimizing integer on this
target is \(j=1\). Its mean absolute discrepancy is
\(7/300\), while (NB2) requires at least \(1/42\) after
subtracting the constant multiplier term.

### The symmetric diffuse profile satisfies the new inequality

Take \(\nu_i=(1/6)\delta_3+(1/6)\operatorname{Unif}[2,4]\) on
each letter. Then \(\kappa p\le1/3\) and \(\kappa p^2\le1/18\)
for \(p\le H_A(\kappa)\). Any violation of (NB2) by an assignment
would have \(h\ge1/3\ge2p\), \(p(M+1)/2>h\), and therefore
\(M\ge4\) and \(E<1/24\). Put \(x=p(M+1)/(2h)>1\).
Using \(1/M\le p/((2h-p)x)\) and
\((x-1)/x^2\le1/4\) in the failed inequality gives
\[
 E<\frac{\kappa p^2}{8(2h-p)},\qquad Eh\le\frac1{216}.
\]
If \(a\le2\), whole periods of the nearest-grid distance in
\([2,4]\) give \(E\ge a/48\). The numerical moment gives
\(ah+E\ge1\), hence \(Eh\ge(1-E)/48>23/1152\), which contradicts
the upper bound. If \(a>2\), the interval of length two on the
diffuse part occupies at most one period. Its folded distance
distribution has unnormalized density at most two, so its distance
integral is at least one before multiplying by the density \(1/12\).
This gives \(E\ge1/12\), again a contradiction.

## Projective heights of sparse monomial supports

Work at \(d=1/3,\Gamma=1/2\) and let \(A,T\) be distinct
letters. Keep the notation \(H_A(\kappa)\) of (NB2).
Use a spacing \(a>0\) and a bounded positive integer
assignment \(i(r)\) on the target, constant on finitely
many intervals with continuity endpoints. For a support
\[
 \mathcal P=\bigcup_{j=0}^{M} (I_j\times\{j\}),
 \qquad B_{\mathcal P}=\int \operatorname{width}_{\mathcal P}
                  \big(i(r),r-a i(r)\big)\,d\nu_T(r),
\]
each \(I_j\) is a finite union of compact intervals.
Empty slices are allowed. The width of \(\mathcal P\) in
direction \(\eta\) is
\(\max_{\mathbf m\in\mathcal P}\mathbf m\cdot\eta-
  \min_{\mathbf m\in\mathcal P}\mathbf m\cdot\eta\).
Write \(V=\sum_j |I_j|\), where \(|I_j|\) is length.
For every such support with \(0<V<\kappa\), one has
\[
 B_{\mathcal P}\ \ge\ \frac V2\min\{1/6,H_A(\kappa)\}.
 \tag{PS1}
\]
A common horizontal translation of the slices leaves
both sides unchanged. Rectangles recover the
fixed-assignment form of (NB2) as \(V\uparrow\kappa\).

### The primitive evaluation rows

Fix \(p<\min\{1/6,H_A(\kappa)\}\) positive. As for (NB2),
freeze outside row factors of mass \(p\) in a compact
continuity band strictly above \(\kappa\). Two distinct
rows in one fibre have near-full common columns, with
the shared full factor giving \(v z^q\equiv1\) modulo
\(F_0\), \(qL\to a\). Its radical \(Q\) has normalized
height \(p+o(1)\) and \(kL\to\kappa\) is an admissible
Taylor depth. Use the same assignment \(i\) on the two
shifted targets \(T_1,T_2\).

Translate the slices if needed to use nonnegative narrow
degrees. At fixed \(L\), the monomials are \(z^u v^j\) with
integer \(u\) and \((Lu,j)\in\mathcal P\). Their number
\(D_0\) satisfies \(LD_0=V+O_{\mathcal P}(L)\).
For each evaluation row, choose primitive integer
coordinates for the rational monomial vector. Its logarithmic
maximum norm is
\[
 \max_{(Lu,j)\in\mathcal P}\log|z^u v^j|
 -\sum_{\ell}\min_{(Lu,j)\in\mathcal P}
             v_\ell(z^u v^j)\log\ell.
 \tag{PS2}
\]
Indeed the sum accounts for the unique positive rational
scaling that makes the row primitive integral. The
archimedean maximum is negligible in the height limit,
since the two targets come from the same dyadic root
box and all degrees are fixed in that limit. Exceptional
coefficients can be fixed in the anatomic extraction.

For an original prime contribution at scaled degree \(r=eL\),
the target supplies the direction
\((i(eL),eL-qL i(eL))\), divided by \(L\) in the valuation
of an evaluation entry. The two targets carry opposite signs
in (PS2). A shared prime causes no additional cost, since
the minimum of a sum is at least the sum of the minima.
The two targets have the same limiting root measure and
no escaped first moment. Their combined contribution to
the second term of (PS2), divided by \(\log X\), is therefore
at most \(B_{\mathcal P}+o(1)\). This is the ordinary
projective height of a monomial vector evaluated using
the actual support.
[Burgos Gil–Philippon–Sombra, canonical toric metric and point height](https://arxiv.org/pdf/1403.4048)

### Taylor minors and the valuation count

In a fixed \(z\bmod Q\) class, the primitive row scalars
are units at the primes dividing \(Q\). The first Taylor
determinant thus has divisibility \(Q^{D_0(D_0-1)/2}\).
Its precision uses \(D_0<k\), by \(V<\kappa\).
The archimedean bound per row is \(B_{\mathcal P}+o(1)\)
in the \(X\)-exponent. If \(B_{\mathcal P}<pV/2\), every
maximal-size determinant vanishes and each chart has a
nonzero integral auxiliary polynomial. Its coefficient
height is polynomial in \(X\) at fixed truncation, by
solving a maximal-rank subsystem using minors.
The Newton breakpoint count applies because \(M\) and
the multiplier list are fixed. Primes on the two targets
coincide only at primes of their fixed difference. Away
from those primes and the nonzero coefficient product,
the auxiliary confines the target's outside primes to
a bounded list of degree positions per chart.
Their normalized cost tends to zero, while the
shallow chart cost is \(p<1/6\). This contradicts the
near-full common columns. Taking \(p\) up to the indicated
minimum gives (PS1).

### A profile satisfying the rectangular inequalities

Put \(c=1/6\) and consider the profile
\[
 \nu_A=c\delta_{7/4}+c\operatorname{Unif}[4,4.5],\qquad
 \nu_B=c\delta_3+c\operatorname{Unif}[2,4],\qquad
 \nu_T=c\delta_3+c\operatorname{Unif}[2.7,3.3].
\]
All masses are \(1/3\) and all moments are one.
For \(\kappa=4\), take \(p=83/500\), \(a=3\), \(i=1\)
and \(M=7\). Center the slices at a common horizontal
coordinate. At \(j=0,7\) take width \(297/1000\), and at
\(j=1,\ldots,6\) take width \(561/1000\). Then
\(V=99/25<4\) and the width in direction \((1,t)\) is
\[
 \max\{561/1000+5|t|,\ 297/1000+7|t|\}.
\]
The transition is \(|t|=33/250\). Integrating the two
linear parts and the atom gives
\[
 B_{\mathcal P}=\frac{1024}{3125}
             <\frac{8217}{25000}=\frac{pV}{2}.
\]
Thus (PS1) excludes the profile in this row-target
orientation. The gap is \(1/1000\).

The profile satisfies (NB2) in every orientation.
To see this, let \(p\le H_A(\kappa)\) and use the
fixed-assignment form with \(h=\int i\,d\nu_T\) and
\(E=\int|r-ai|\,d\nu_T\). We have \(h\ge1/3\),
\(\kappa p\le2/3\), and \(\kappa p^2\le1/9\).
Its right side for \(E\) is at most
\[
 \frac23\max_{M\ge1}\frac{(1/2-2/(M+1))_+}{M}
 =\frac1{42}.
\]
The maximum follows from
\((M-6)(M-7)\ge0\) for integers \(M\).
If \(a\) exceeds the width of the target's diffuse interval,
the interval lies within one period of distance to
\(a\mathbb Z\). Its folded-distance density before
normalization is at most two, giving \(E\ge1/40\) on \(T\)
and \(E\ge1/12\) on \(B\). If \(A\) is the target, the
other two letters have \(\kappa p\le9/20\), so the
rectangular right side is at most \(9/560<1/48\), while
\(E\ge1/48\). For \(a\) at most the width of the target
interval, complete periods give \(E\ge a/48\), and the
numerical moment gives \(ah+E\ge1\).
A failed inequality would have \(E<1/42\). With
\(x=p(M+1)/(2h)>1\), the calculation using
\((x-1)/x^2\le1/4\) gives
\(Eh<\kappa p^2/12\le1/108\), whereas
\(Eh\ge(1-E)/48\). This is a contradiction.
The sparse support thus supplies a strict additional
restriction.

### The broad symmetric profile and sparse supports

For \(\nu_i=c\delta_3+c\operatorname{Unif}[2,4]\) on each
letter, (PS1) holds for all its supports with one grouped
root and positive integer multipliers. Suppose a support
gave \(B_{\mathcal P}<pV/2\). Let \(m\) be the difference
of its largest and smallest occupied wider degrees, and
\(R\) its largest horizontal slice diameter. Then
\(V\le(m+1)R\) and \(B_{\mathcal P}\ge hR\), so
\(h<p(m+1)/2\). As \(h\ge1/3\ge2p\), this requires \(m\ge4\).
Choose monomials in the two extreme wider slices.
Their difference in width gives \(B_{\mathcal P}\ge mE'\)
for \(E'=\int|r-a'i(r)|\,d\nu_T(r)\), with \(a'\) the
spacing shifted by the difference of their scaled narrow
degrees divided by \(m\). The row support gives
\(\kappa p\le1/3\) and \(\kappa p^2\le1/18\), hence
\[
 E'<\frac{\kappa p}{2m}\le\frac1{24},
 \qquad
 E'h<\frac{\kappa p^2(m+1)}{4m}\le\frac5{288}.
\]
A spacing \(a'\le0\) gives \(E'\ge1\). If \(a'>2\) the
folded-distance density on the diffuse interval gives
\(E'\ge1/12\). Finally if \(0<a'\le2\), whole periods
give \(E'\ge a'/48\), and \(a'h+E'\ge1\) gives
\(E'h>23/1152>5/288\). Each case contradicts a strict
determinant margin.

## Two-root Newton count by specialization

Work at \(d=1/3,\ \Gamma=1/2\), and put \(c=1/6\). Use the pair-fibre bound (EN2) for the outside labels on a free pair \(U,V\). On these letters divide the bounded scaled-degree coordinates into two groups, with disjoint prime supports. Assume all atoms of this pair belong to group 1. Let the combined nonatomic masses of groups 1 and 2 be \(t_1,t_2\), so \(t_1+t_2=2c\).

Use one unweighted grouped root in each group, based at degree \(a_i>0\), so \(q_iL\to a_i\). The ratio of group root products is denoted by \(z_i\), and the wider signed correction for the full \(U/V\) ratio by \(v\). Let \(h_i\) be bounds for their normalized narrow heights, and \(E\) a bound for the wider numerical rational height. The usual negligible exceptional factors can be included in the correction or fixed at negligible cost. On the remaining letter freeze bounded-scale roots of cost \(\mu\) and full product \(n\), of fullness \(kL\to\kappa>0\). Let \(Q\) be the shallow modulus used to partition both \(z_i\) coordinates, of normalized weight \(b_Q\), with \(b_Q\ge b\) and \(v_\ell(Q)\le\lfloor v_\ell(n)/k\rfloor\) at its primes.

For any integer \(M\ge1\) and constants \(r_1,r_2>0,\ 0<\alpha<1\), set
\[
 y=\sqrt{2\alpha(M+1)r_1r_2}.
\]
If
\[
 y<\kappa,\qquad
 h_1r_1+h_2r_2+ME < \frac23 b\, y,
 \tag{TN1}
\]
then the normalized count is at most
\[
 \mu+2b_Q+\max(t_1,t_2)
 \tag{TN2}
\]
up to arbitrarily small outside-window and extraction losses. In particular strict inequality of this coefficient from \(2c\) excludes the profile. The Taylor determinant supplies the auxiliary polynomial. Specialization and valuation breakpoints then count the labels.

### First determinant

Primitivity and the sum relation give a unit congruence \(v z_1^{q_1}z_2^{q_2}\equiv K\bmod n\) for each frozen modulus. In each pair of residue classes of \(z_1,z_2\bmod Q\), use rectangular monomial support of degrees at most \(R_i=\lfloor r_i/L\rfloor\) in \(z_i\) and at most \(M\) in \(v\). Write \(D_0=(R_1+1)(R_2+1)(M+1)\) and \(N=\lfloor\alpha D_0\rfloor\). Denominator clearing multiplies each row by a unit modulo \(n\) and gives integer entry heights at most \(X^{h_1r_1+h_2r_2+ME+o(1)}\).

At each \(\ell\mid Q\), substitute \(v=K z_1^{-q_1}z_2^{-q_2}\) in the exact local reference, preserving congruence modulo \(\ell^{v_\ell(n)}\). In the two shallow root variables the signed binomial Taylor coefficients are integral. Distinct Taylor vectors have weighted divisibility at least \(v_\ell(Q)\) times their total orders. Determinants with two occurrences of the same Taylor vector vanish. The \(N\) smallest two-variable total orders have largest order \(\sqrt{2N}+O(1)\) and mean order \(\frac23\sqrt{2N}+O(1)\), by summing the \(s+1\) monomials of order \(s\). The first condition in (TN1) gives largest order below \(k\), so replacement errors of valuation at least \(v_\ell(n)\) can be treated at order \(k\) or higher. Thus the prime-adic divisibility exponent per row in an \(N\)-row minor is at least \((2/3)b y+o(1)\). The second condition in (TN1) forces all such minors to vanish.

Selecting a nonzero integer kernel vector gives a polynomial \(F(z_1,z_2,v)\) in the prescribed rectangle that vanishes on all original sample tuples in this chart. Its coefficients have polynomial height at fixed truncation, by using maximal minors of the evaluation matrix. Polynomial height suffices for the divisor allocations below.

### Specialization and outside labels

Specialize the exact value of \(z_2\). There are only boundedly many exceptional rationals for which \(F(z_1,z_2,v)\) becomes the zero polynomial, because at least one of its coefficient polynomials in \(z_2\) is nonzero. For each fixed exceptional value, its reduced numerator and denominator recover the group-2 original roots up to divisor-many allocations: \(z_2\) is their unweighted ratio and the two free letters have disjoint primes. Consequently the group-2 outside vectors are determined up to subpower choices, while the group-1 outside factors cost at most \(t_1\) up to outside-window errors.

For nonexceptional \(z_2\), enumerate it at normalized cost at most \(t_2\). This is where excluding all atomic factors from group 2 is used. The specialized polynomial has polynomial coefficient height after clearing its \(z_2\) denominators. Take \(J\) to be the product of its nonzero integer coefficients. For a group-1 outside prime not dividing \(J\), its valuations at \(z_1,v\) are proportional, with common signed factor, to \((1,e-q_1)\), where \(e\) is its original degree. Primitivity across letters and disjoint groups ensure that no other original factor changes these valuations. At each exponent \(j\) of \(v\), choose the unique extremal exponent of \(z_1\) giving the smallest valuation. The minimum among the resulting affine functions of \(e\) must tie, so at most \(2\binom{M+1}{2}\) original degree positions are exceptional across the two letter signs. At every other group-1 position the primes lie in \(J\), and divisor allocation counts all their contributions at subpower cost. Each moving exceptional position has outside root cap arbitrarily small. Hence the two free outside vectors are determined up to subpower and arbitrarily small normalized errors.

The pair-fibre bound (EN2) converts either outside-label count to the triple count. The number of shallow charts contributes \(2b_Q\). Include the initial modulus cost \(\mu\) to obtain (TN2).

### An explicit excluded profile

Let
\[
 \phi=\frac{1+\sqrt5}{2},\qquad \rho=\frac1{6+\phi},\qquad
 \delta=\frac3{500},\quad \sigma=\frac1{500},
\]
and put \(r_0=(12-\delta)/(3+\rho)\), \(r_1=r_0+\delta\), and \(r_2=\rho r_0\). Give both free letters the measure
\[
 \nu_U=\nu_V=c\delta_{r_0}
       +\tfrac1{12}\operatorname{Unif}[r_1-\sigma,r_1+\sigma]
       +\tfrac1{12}\operatorname{Unif}[r_2-\sigma,r_2+\sigma].
\]
On the modulus letter put
\[
 \nu_T=c\delta_{657/160}
          +\tfrac1{16}\operatorname{Unif}[4.5,4.6]
          +\tfrac5{48}\operatorname{Unif}[0,0.6].
\]
All masses are \(1/3\), moments are one, and all nonatomic masses are \(c\). In group 1 put the atom \(r_0\) and the interval around \(r_1\) of each free letter, using base \(a_1=r_0\). In group 2 put the interval around \(r_2\) of each free letter, using base \(a_2=r_2\). Then \(t_1=t_2=c\), \(h_1=1/4,\ h_2=1/12\), and
\[
 E=\frac{\delta}{12}+\frac{\sigma}{24}=\frac{7}{12000}.
\]
Freeze normalized root mass \(0.0554\) from the modulus band's nonatomic mass \(1/16\) in \([4.5,4.6]\). Use its radical as \(Q\), with \(b=0.05535,\ b_Q\le0.05545,\ \mu\le0.05545,\ \kappa=4.5\). Set \(M=60,\ r_1'=0.2,\ r_2'=0.6,\ \alpha=0.99\) in (TN1), where the primed variables denote determinant degree parameters rather than the group centers. Then
\[
 y=\sqrt{14.4936}\in(3.806,3.808),\qquad
 B=\tfrac14(0.2)+\tfrac1{12}(0.6)+60(7/12000)=0.135.
\]
The divisibility lower bound \((2/3)(0.05535)y\) exceeds \(0.1404\). The coefficient in (TN2) is at most \(c+3(0.05545)<1/3\), a contradiction.

### Comparison with the rectangular one-root inequality

The fixed-assignment bound from (NB2) holds in every orientation. For any supplying letter of these profiles, \(\kappa p^2\le3.85/144\) and \(\kappa p\le3.85/12\) for any tail cutoff \(\kappa\) and eligible diffuse mass \(p\). This uses the two disjoint support bands and bounds their complete masses. For any positive assignment \(i(r)\), \(h\ge1/3\). A violating \(M\) must satisfy \(M\ge4\), and then \(E<\kappa p/(2M)<0.041\). The rectangular bound also yields \(E<\kappa p^2/[8(2h-p)]\): set \(x=p(M+1)/(2h)>1\), use \(1/M\le p/[x(2h-p)]\), and maximize \((x-1)/x^2\le1/4\). Consequently \(Eh<\kappa p^2/12\le3.85/1728\).

On a free letter, in contrast, \(E/a\ge1/400\) for every \(a>0\). For \(a\le0.004\) use complete mesh periods in the diffuse intervals, and for \(0.004<a\le0.012\) use the minimum average distance on an interval no longer than a period. For \(a>0.012\), write \(d_j=\operatorname{dist}(r_j,a\mathbb Z)\). Since \(d_0+d_1\ge\delta\), and evaluating each uniform-interval integral instead at its center loses at most \(\sigma/2\) times its mass, one has
\[
 E\ge(\delta-\sigma)/12+(d_0+d_2)/12.
\]
This gives the claimed lower bound for \(a\le2/15\). For \(2/15<a\le0.6\), let \(n,m\) be nearest integers to \(r_0/a,r_2/a\), with \(n\ge6\). The algebraic conjugates \(\rho,\rho'\) solve \(41z^2-13z+1=0\). Thus
\[
 |m-\rho n||m-\rho'n|\ge1/41,\qquad \rho'-\rho=\sqrt5/41.
\]
If \(|m-\rho n|<1/(3n)\), the product times 41 is strictly less than \((\sqrt5+41/(3n^2))/3<1\), impossible. It follows that \(d_0+d_2\ge a/(3n)\), since \(0<\rho<1\). As \(n\le r_0/a+1/2\), this yields \(E/a\ge1/(3000a)+a/150>1/400\) by the arithmetic-geometric mean inequality. For \(a>0.6\) the lower cluster is entirely below the first positive assigned mesh point, so \(E/a\ge(1-0.51/a)/12>1/400\). On the modulus letter the broader interval \([0,0.6]\) gives \(E/a\ge(5/48)/8>1/400\) by complete periods when \(a\le0.6\), and its entire mass below the first positive point gives even more for \(a>0.6\). Thus in all cases \(h\ge(1-E)/a\) gives \(Eh>0.959/400>3.85/1728\) whenever \(E<0.041\), a contradiction to a violation.



### Two-root projective row heights

If the free letters have identical limiting degree measures,
the same proof permits sparse horizontal supports
\(\mathcal P=\bigcup_j(P_j\times\{j\})\subset\mathbb R^2\times
\{0,\ldots,M\}\), with \(P_j\) bounded finite unions of
rectangles. Set \(V=\sum_j\operatorname{area}(P_j)\).
Use \(y=\sqrt{2\alpha V}\) in (TN1) and replace its row-height
bound by
\[
 \int\operatorname{width}_{\mathcal P}
     \big(i_1(r),i_2(r),r-a_1i_1(r)-a_2i_2(r)\big)\,d\nu_U(r).
 \tag{TN3}
\]
Here \(i_i\) is the indicator of the corresponding group.
The archimedean logarithms in the projective row formula
tend to zero on the numerical scale, using equality of the
limiting profiles and tightness of their numerical moments.
The two signed finite-place sums give the width integral
as in (PS2). The monomial count has \(L^2D_0\to V\), so the
same Taylor orders and specialization bound apply.

### A count with several specialized roots

The specialization argument extends to \(m\ge2\) disjoint
groups on the two free letters, with all atoms confined
to group one. Let the pooled outside masses be \(t_j\)
with sum \(2c\). Use the analogous congruence with one
grouped ratio \(z_j\) per group and shallow residues
modulo \(Q\) in all \(m\) coordinates.
For a rectangular support put
\(V=(M+1)\prod_{j=1}^m r_j\) and
\(B=\sum_j h_j r_j+ME\). For the projective support when
the two free profiles agree, take \(V\) to be the sum of
the \(m\)-dimensional slice volumes and \(B\) the
corresponding width integral. Put \(y=(m!\alpha V)^{1/m}\).
The smallest Taylor multiindices give largest order
\(y/L+O_m(1)\) and mean order
\(my/((m+1)L)+O_m(1)\). Thus the two conditions
\[
 y<\kappa,\qquad B<\frac{m}{m+1}b y
 \tag{TN4}
\]
give a nonzero auxiliary polynomial on each chart.

Specialize \(z_2,\ldots,z_m\). Their dictionaries have normalized
exponents at most \(t_2,\ldots,t_m\), since they involve only
outside mass. At a specialization where the auxiliary survives,
the Newton breakpoints account for the group-one outside factors
and the cost is \(\sum_{j\ge2}t_j\). The exceptional specializations
lie on a fixed nonzero coefficient polynomial in \(m-1\)
variables. Over finite dictionaries \(S_j\), a nonzero polynomial
of fixed degrees has at most
\(O(\sum_{j\ge2}\prod_{i\ge2,\,i\ne j}|S_i|)\)
zeros. One proves this by choosing a nonzero coefficient
polynomial in the last variable, counting its exceptional
specializations by induction, and using the degree bound
elsewhere. Counting exceptional vectors and then the
group-one outside labels costs at most
\(2c-\min_{j\ge2}t_j\).
In both branches the outside pair-fibre estimate bounds
the number of full triples with an arbitrarily small
normalized loss. Including the modulus and chart costs gives
\[
 d\le \mu+m b_Q+2c-\min_{1\le j\le m} t_j
 \tag{TN5}
\]
whenever (TN4) holds.

### Limits of rank information in several roots

Fix \(m\ge2\) and finitely many nonzero integer multiplier
directions \(\mathbf i\in\mathbb Z^m\). For each direction choose
a primitive nonzero integer vector \(\mathbf u_{\mathbf i}\)
orthogonal to it. Take the Laurent product of the distinct
binomials \(\mathbf z^{\mathbf u_{\mathbf i}}-1\), identifying
duplicates up to Laurent units, and multiply by a monomial to
obtain a polynomial \(G\). At every weight
\((\pm\mathbf i,s)\) the corresponding binomial has both
terms in its initial form. Initial forms multiply, and a
non-monomial Laurent polynomial cannot divide a monomial.
Thus every nonzero multiple of \(G\) has a repeated least
weight for these directions and for every wider slope \(s\).

The multiples of \(G\) in a rectangular monomial space
of narrow degrees \(R_j\to\infty\) and fixed wider degree
\(M\) have dimension fraction tending to one. Multiply the
space with each \(R_j\) reduced by \(\deg_{z_j}G\) by \(G\).
These multiples can in fact be an exact evaluation kernel:
the positive rational points on the union of the split
subtori \(\mathbf z^{\mathbf u_{\mathbf i}}=1\), with arbitrary
positive rational wider coordinate, are Zariski dense there.
Take finitely many evaluations spanning the restriction
dual. This is an abstract rank obstruction. Its arithmetic
sample has no prescribed prime power valuations or height caps.

There is also a dimension cost to forcing unique weights
on the monomial support. If a single fixed multiplier
direction \(\mathbf i\) of target root mass \(\tau>0\)
gives an injective weight on each wider slice, and \(W_L\)
is the largest within-slice integer width in that direction,
there are at most \((M+1)(W_L+1)\) columns. The projective
width bound is at least \(\tau L W_L\). For a growing
number \(D_L\) of columns this is of order at least
\(L D_L\), while the ordinary \(m\)-variable Taylor orders
give divisibility per row of order at most \(L D_L^{1/m}\).
With \(M\) fixed and \(m\ge2\) these bounds cannot force
rank drop. The specialization method handles repeated
weights by counting the exceptional specializations.

### An exclusion using two roots and projective row heights

Keep the profile of the two-root specialization example, changing only \(\delta\) to \(47/5000=0.0094\) in the definitions of \(r_0,r_1,r_2\). Keep \(\sigma=1/500\), the same two free groups, and the same modulus letter and shallow record. The rational wider height now has coefficient \(E=13/15000\), and the normalized grouped-root heights are still \(h_1=1/4,\ h_2=1/12\). The two-root chart budget \(\mu+2b_Q<c\) still holds.

Take \(M=66\). Center each horizontal rectangle at the same horizontal point. For \(|j-33|\le26\), take widths \((0.26,0.655)\); for the other fourteen levels take widths \((0.135,0.61)\). Thus
\[
 \mathcal A=53(0.26)(0.655)+14(0.135)(0.61)
             =\frac{25447}{2500}=10.1788.
\]
The first group's projection width at residual \(u\) is
\(\max(0.26+52|u|,0.135+66|u|)\). Its atom has \(u=0\), and its diffuse residual \(u\) is uniform on \([\delta-\sigma,\delta+\sigma]\), so integrate at breakpoint \((0.26-0.135)/14\). The second group's width is \(\max(0.655+52|u|,0.61+66|u|)\), and its diffuse residual satisfies \(|u|\le0.002\), where the first branch applies. Hence
\[
 B(P)=\frac{5562169}{33600000}<0.165541,\qquad
 \sqrt{2(0.99)\mathcal A}\in(4.4893,4.4894).
\]
With \(b=0.05535\) and \(\kappa=4.5\), one has \(y<\kappa\) and \((2/3)by>0.16565\). Thus (TN1) with the row bound (TN3), together with (TN2) exclude this explicit profile.

For comparison, with the same two grouped roots, no rectangular Taylor criterion (TN1) has enough height room, for any permitted shallow divisor and modulus root selection. For any such choice write \(b=b_Q\) (the largest permitted bound for that selector) and \(C=2/3,\ K_0=\sqrt{2h_1h_2}=1/\sqrt{24}\). The chart budget requires \(\mu+2b<c\). The available numerical modulus weight gives \(\kappa b\le w\le4.6\mu\), since all modulus degrees have scaled support at most \(4.6\) and there is no escaped moment. At a putative rectangular determinant, AM–GM gives
\[
 E<\frac{\kappa}{M}\left(Cb-\frac{K_0}{\sqrt{M+1}}\right).
\]
Indeed for \(y=\sqrt{2\alpha(M+1)r_1r_2}<\kappa\), the narrow height contribution \(h_1r_1+h_2r_2\) is at least \(K_0 y/\sqrt{\alpha(M+1)}\). Put \(x=K_0/(Cb\sqrt{M+1})\in(0,1)\). It follows that
\[
 E<
 4.6\, C^3(c-2b)b^2
 \frac{x^2(1-x)}{K_0^2-C^2 b^2x^2}.
\]
Use \((c-2b)b^2\le c^3/27\) and \(C^2b^2/K_0^2\le2/27\). For \(x\le37/50\) bound the numerator factor by \(4/27\) and the denominator factor below by \(1-(2/27)(37/50)^2\). For \(x>37/50\) use that \(x^2(1-x)\) decreases on this interval, and bound the denominator factor below by \(1-2/27\). The resulting respective upper bounds on \(E\) are \(0.0008660791\) and \(0.0008624588\), both strictly below \(13/15000\). The comparison fixes the two target groups.

## Counting specialization values together with their residue classes

In the setting of (TN4), the count satisfies
\[
 d\le\mu+\max\left(
    b_Q+2c-t_1,\ 2b_Q+2c-\min_{j\ge2}t_j
 \right)
 \tag{TC1}
\]
up to arbitrarily small extraction and truncation losses. For \(m=2\)
these two terms inside the maximum are \(b_Q+t_2\) and \(2b_Q+t_1\).
The bound also applies if some of the outside labels have already been
frozen. In that case \(t_j\) are the retained outside masses, and the
cost \(2c-\sum_j t_j\) of freezing the rest is included in (TC1).
Freezing atomic mass adds its own root cost.

For the specialization where the polynomial survives, choose the
exact \(z_2,\ldots,z_m\) from their full dictionaries and the residue
of \(z_1\). Each exact rational value is a unit at the modulus primes and
determines its own residue. These choices determine
the full chart and hence its fixed polynomial. The breakpoint argument
recovers the group-one outside factors with arbitrarily small loss.
The normalized cost after fixing the modulus and before counting
any frozen outside labels is \(b_Q+\sum_{j\ge2}t_j\).

For the exceptional specialization, let \(S_{i,r}\) be the subdictionary
of exact \(z_i\) in the residue \(r\bmod Q\), and \(S_i\) the full
dictionary, for \(i\ge2\). On every chart
\((r_1,\ldots,r_m)\), use the same polynomial zero estimate as in
(TN5) with these subdictionaries. The bound is a constant times
\(\sum_{j\ge2}\prod_{i\ge2,\,i\ne j}|S_{i,r_i}|\), even though the
coefficient polynomial itself may vary with the chart. Sum over the
charts. For each \(j\) the product is independent of \(r_1,r_j\),
giving at most \(Q^2\) for those choices. For each other \(i\) its
sum is exactly \(|S_i|\). The exceptional vectors consequently have
normalized cost at most
\(2b_Q+\sum_{i\ge2}t_i-\min_{j\ge2}t_j\). Their exact values recover
the corresponding outside labels up to divisor allocations. Adding
\(t_1\) for the remaining outside labels and applying (EN2) proves
(TC1).

### An endpoint profile excluded by the improved chart count

Put
\[
 \rho=\frac{13-\sqrt5}{82}=\frac1{6+\phi},\quad
 d_1=\frac3{50},\quad d_2=\frac8{75},\quad
 \delta=\frac3{250},\quad \sigma=\frac1{500},
 \qquad \phi=\frac{1+\sqrt5}{2}.
\]
Set
\[
 r_0=\frac{75-(9/2)\delta}{17+8\rho},\qquad
 r_1=r_0+\delta,\qquad r_2=\rho r_0 .
\]
Give the two free letters identical measures
\[
 \nu_U=\nu_V=c\delta_{r_0}
        +d_1\operatorname{Unif}[r_1-\sigma,r_1+\sigma]
        +d_2\operatorname{Unif}[r_2-\sigma,r_2+\sigma],
\]
and give the modulus letter
\[
 \nu_T=c\delta_{657/160}
        +\tfrac1{16}\operatorname{Unif}[4.5,4.6]
        +\tfrac5{48}\operatorname{Unif}[0,0.6].
\]
All weights sum to \(1/3\), all first moments are one, and the nonatomic masses are \(c\). In fact \((c+d_1)r_0+d_1\delta+d_2r_2=1\). Group 1 contains the atom at \(r_0\) and the interval around \(r_1\), grouped at \(r_0\); group 2 contains the interval around \(r_2\), grouped at \(r_2\). The resulting parameters are
\[
 t_1=2d_1=0.12,\quad t_2=2d_2=\frac{16}{75},\quad
 h_1=\frac{17}{75},\quad h_2=\frac8{75},\quad
 E=d_1\delta+d_2\sigma/2=\frac{31}{37500}.
\]
Select modulus roots of mass \(0.0592\) from the nonatomic \([4.5,4.6]\) band on \(T\). The band has mass \(1/16>0.0592\). Its radical gives \(b\ge0.05915,\ b_Q\le0.05925,\ \mu\le0.05925\), with arbitrarily small extraction losses and \(\kappa=4.5\).

In the rectangular first determinant take \(\alpha=0.99,\ M=69\) and horizontal degree coefficients \(r_1'=0.26,\ r_2'=0.55\). Then
\[
 y=\sqrt{2(0.99)(70)(0.26)(0.55)}
      =\sqrt{19.8198}\in(4.451,4.453),\qquad
 B=0.17464<\tfrac23(0.05915)y.
\]
The nonexceptional count costs at most \(0.1185+16/75<1/3\), and the exceptional count at most \(3(0.05925)+0.12<1/3\). Thus (TC1) excludes the profile at \(d=1/3\) and \(\Gamma=1/2\).

The old budget \(\mu+2b_Q+\max(t_1,t_2)<1/3\) cannot prove a saving by the rectangular first determinant with these same groups even after optimizing the determinant parameters, full modulus selection and shallow divisor. In fact it requires \(\mu+2b<0.12\) with \(b\le b_Q\), and \(\kappa b\le w\le4.6\mu\). Put \(C=2/3\), \(K^2=2h_1h_2\), and \(x=K/(Cb\sqrt{M+1})\in(0,1)\). AM--GM on the narrow height terms, followed by the Taylor comparison and its precision condition, requires
\[
 E<4.6 C^3(0.12-2b)b^2
       \frac{x^2(1-x)}{K^2-C^2 b^2 x^2}
 \le \frac{4.6 C^3(0.12)^3(4/729)}
              {K^2-C^2(0.12)^2/4}<0.000277 .
\]
This uses \((0.12-2b)b^2\le(0.12)^3/27\) and \(x^2(1-x)\le4/27\). It is incompatible with \(E=31/37500\).

### One-root compatibility of the same profile

In every orientation of these three profiles the maximal supplying tail values satisfy
\[
 \kappa p\le\frac{23}{80},\qquad
 \kappa p^2\le\frac{23}{1280}.
 \tag{TC2}
\]
For \(T\) this follows by considering the \([0,0.6]\) component and then the mass \(1/16\) at degrees at most \(4.6\). For the two identical free letters the lower band is below \(0.55\), and above it the tail mass is at most \(3/50\) at degree at most \(4.17\).

Recall the general necessary product test for a successful sparse projective one-group determinant. If the horizontal volume is \(\mathcal A\), the largest horizontal slice width \(R\), and the gap between extreme wide indices \(J\), then the width height \(B\) satisfies \(B\ge hR\) and \(B\ge J E'\), where \(E'\) is the degree discrepancy at the induced shifted mesh spacing. Also \(\mathcal A\le(J+1)R\). Thus \(B<p\mathcal A/2\), \(\mathcal A<\kappa\) would give \(h<p(J+1)/2\) and \(E'<\kappa p/(2J)\). For the positive multipliers of (PS1), \(h\ge1/3\), so \(J\ge4\).
A nonpositive induced mesh gives \(E'\ge1\) and is incompatible with the width bound.
For a positive mesh (TC2) now gives
\[
 E'<23/640,\qquad E'h<\frac5{16}\,\frac{23}{1280}=\frac{23}{4096}.
 \tag{TC3}
\]

On either free letter we prove \(E'/a\ge0.00588\) for every positive induced mesh \(a\) and every positive assignment. For \(a\le0.004\), integrate complete mesh periods in both diffuse intervals to get \(E'/a\ge c/8\). For \(0.004<a\le0.024\), minimum average distance in a uniform interval of length \(0.004\) gives numerical \(E'\ge0.001c=1/6000\), so \(E'/a>0.00588\). For \(a\ge1\), the lower interval lies below the least positive mesh point and on its own gives the result.

Suppose instead \(0.024<a<1\) and \(E'/a<0.00588\). Both diffuse intervals lie entirely in a single nearest-integer phase cell. Indeed their mean distances to the mesh modulo one are respectively at most \(0.00588/d_1\) and \(0.00588/d_2\), so the center's distance differs by at most \(\sigma/(2a)\le1/24\) from this upper bound; any point in the interval differs in turn by at most \(\sigma/a\le1/12\) from its center. These sums are strictly below \(1/2\), which excludes touching a nearest-cell boundary. Jensen's inequality therefore lower bounds each interval's contribution by its center distance, with no rounding loss.

Put \(x=r_0/a\), and let \(n\) be the nearest integer. The phase discrepancy from the atom implies \(|x-n|<441/12500\). Since \(\delta/a<1/2\), the phase distances of the atom and the nearby upper center sum to at least \(\delta/a\). Using \(c-d_1=d_2\), and \(\rho<1\), yields
\[
 E'/a\ge (d_1\delta/r_0)x+d_2(\|x\|+\|\rho x\|)
       \ge \frac9{52000}\left(n-\frac{441}{12500}\right)
                            +\frac8{75}\|\rho n\|.
\]
Also the first term \((d_1\delta/r_0)x\) before replacing \(x\) by \(n\) implies \(x<34\), hence \(4\le n\le34\). Use the rational interval \(0.13126746<\rho<0.13126747\), obtained by squaring for the expression \((13-\sqrt5)/82\). For each \(n=4,\ldots,34\), replace \(\|\rho n\|\) by the smaller of the two endpoint distances, or zero if the rational interval contains an integer. The minimum of the resulting expression is at \(n=15\), with exact value \(57480433/9750000000>0.00588\). This finite check uses only the displayed endpoints and integer floors. It contradicts the hypothetical small discrepancy.

On \(T\), the diffuse interval \([0,0.6]\) alone gives \(E'/a\ge(5/48)/8\) when \(a\le0.6\) by complete periods, and for \(a>0.6\) it lies entirely below the first positive mesh point and gives more. Thus all orientations have \(E'h\ge (1-E')E'/a\ge(1-23/640)(0.00588)>23/4096\) under the putative \(E'\) bound in (TC3). This proves compatibility with (PS1) in every orientation.

### The broad symmetric profile under the improved chart budget

Let all three profiles be \(c\delta_3+c\operatorname{Unif}[2,4]\)
with \(c=1/6\). In (TN4) choose the same disjoint groups on both
free letters and put all retained atoms in group one. The projective
determinant and chart criterion (TC1) are incompatible for every
fixed \(m\ge2\), even after freezing some diffuse or atomic mass and
choosing a shallow divisor. The width argument below also permits
nonzero integer multipliers on the disjoint groups, conditional on
the same determinant and chart costs.

Write \(d_i\) for the retained diffuse mass in group \(i\) on each
free letter, and \(\beta\) for the omitted atom mass per free letter.
Its total retained group mass is \(\ell_i=d_i+(c-\beta)\mathbf1_{i=1}\).
Let \(b\) be the shallow record weight, which may be used as the
determinant lower bound by taking limits. The modulus has numerical
weight \(w\) with \(\kappa b\le w\le4\mu\). The upper bound follows
from the support bound and convergence of the full first moment.
For a strict saving from (TC1), one needs
\[
 d_1>\beta+x,\qquad d_i>\beta+y_0\ (i\ge2),\qquad
 x=(\mu+b)/2,\quad y_0=(\mu+2b)/2 .
 \tag{TC4}
\]
In particular \(b<c\) and \(b>0\) is required by the first determinant.

Let \(\mathcal A\) be the sum of horizontal slice volumes of the sparse
support and \(J\) its gap in the wider index. For the largest
within-slice width \(X_i\) in each narrow coordinate,
\[
 \mathcal A\le(J+1)\prod_i X_i,\qquad
 B\ge\sum_i h_iX_i,\qquad h_i=\int_{\text{group }i}|i_i(r)|\,d\nu(r).
\]
Each group's width is bounded from below by the projection of two
monomials in a common slice, and the group assignments have disjoint
supports. If \(J>0\), two points in the extreme slices give
\(B\ge J\sum_i E_i\), where \(E_i\) is the integral of
\(\lvert r-\widetilde a_i i_i(r)\rvert\) on group \(i\) for the
corresponding shifted real base \(\widetilde a_i\).

For each group \(h_i\ge\ell_i\) and \(E_i h_i\ge (3/4)\ell_i d_i^2\).
If \(E_i\ge\ell_i\) the latter follows from \(\ell_i\ge d_i\)
and \(d_i\le c\). Otherwise \(a=|\widetilde a_i|>0\) and
\(a h_i\ge2\ell_i-E_i>\ell_i\) since \(r\ge2\).
The discrepancy on diffuse mass is at least its distance to
\(a\mathbb Z\). For \(a\le2\), the folded phase density on
\(\|r/a\|\) from the diffuse interval is at most \((2+a)/6\)
by counting preimages of the two signs. Thus the minimum
distance integral at retained mass \(d_i\) is at least
\(3 a d_i^2/(2+a)\ge(3/4)a d_i^2\). For \(a\ge2\), the phase
interval spans at most one period, so the numerical distance
density is at most \(1/6\). The integral is at least \(3d_i^2\).
This proves the scalar bound in both ranges. Using (TC4) gives
\[
 E\prod_i h_i\ge
 \tfrac34(c+x)y_0^{m-1}\bigl(x^2+(m-1)y_0^2\bigr),
 \qquad E=\sum_i E_i .
 \tag{TC5}
\]

By the two Taylor conditions in (TN4) and the arithmetic-geometric
mean bound on \(\sum h_iX_i\), one needs
\[
 \prod_i h_i < \frac{m!}{(m+1)^m}(J+1)b^m .
\]
If \(J=0\), this fails since \(\prod_i h_i\ge(c+x)y_0^{m-1}>b^m\).
For \(J>0\) use \(JE\le B<mb\kappa/(m+1)\) and
\(\kappa b\le4\mu\) to get
\[
 E\prod_i h_i<\frac{4m m!}{(m+1)^{m+1}}\mu b^m(1+1/J).
 \tag{TC6}
\]
Set \(t=\mu/b\) and \(C_m=(m+1)^m(1+t/2)^{m-1}/m!\).
The bound on \(\prod_i h_i\) gives \(J+1>C_m(c+x)/b\).
Consequently
\((c+x)/(1+1/J)>c+x-b/C_m\ge c\), since \(C_m\ge2\)
and \(x\ge b/2\). Divide (TC5) and (TC6) by
\(b^{m+1}(1+1/J)\). A necessary condition is
\[
 \frac18\left(\frac{t+2}{2}\right)^{m-1}
     \frac{(t+1)^2+(m-1)(t+2)^2}{4}
 < \frac{4m m!}{(m+1)^{m+1}}t .
\]
The left side increases with \(m\ge2\) while the coefficient of \(t\)
on the right decreases. For \(m=2\), the reverse strict inequality
is equivalent to
\(P(t)=54t^3+270t^2-565t+270>0\).
Indeed \(P(5/6)=215/12\), \(P'(5/6)=-5/2\), and
\(P''(t)\ge540\) on \(t\ge0\), so \(P(t)\ge215/12-25/4320>0\).
Thus the improved chart criterion is still compatible with this
formal profile. The calculation concerns the stated determinant
and specialization budget with identical group assignments.

## Direct Newton recovery at arbitrary nonatomic mass

In the notation of (QD1), retain bounded scaled-degree windows on two
free letters \(U,V\), and freeze their complement at normalized cost
\(f\). At each retained exact degree \(e\) on letter \(W\), use a
positive integer multiplier \(i_W(e)\le I_0\), with \(I_0\) fixed
before the small-budget limit. The grouped ratio has one narrow
variable \(z\), base exponent \(q\) with \(qL\to a>0\), and the
wider variable \(v\) carrying exponents \(e-qi_W(e)\). The heights
are bounded by \(h,E\) as in (QD1). Freeze a unit-target full
modulus on the other letter at cost \(\mu\). Choose a shallow
divisor of normalized size \(b_Q\ge b>0\) as in (QD1).
The argument below also allows the sparse projective evaluation
supports of (PS1), provided their row-height and Taylor hypotheses hold.

For a rectangular support take \(B=hr+ME\), where the degrees are
at most \(\lfloor r/L\rfloor\) in \(z\) and \(M\) in \(v\).
Assume just the first Taylor conditions
\[
 y=\alpha(M+1)r<\kappa,\qquad B<by/2.
 \tag{DN3}
\]
In each shallow chart, the determinant proof gives a nonzero
auxiliary polynomial \(F\) of wider degree at most \(M\) and
polynomial coefficient height at fixed truncation. The first
determinant requires neither the terminal Jacobian budget nor
complete-label separation.

For either sign \(W\), let \(\mathcal A_W\) be the sum of the
\(I_0\binom{M+1}{2}\) largest retained atomic masses of \(\nu_W\),
counting fewer if the retained measure has fewer atoms. Only
atoms \(r'\) with \(|r'/i_W-a|\le r\) need be included when the
assignment is fixed at their shrinking windows. If an assignment
can split one atom among several multipliers, use the corresponding
retained submeasures for each multiplier. Then
\[
 d\le f+\mu+b_Q+\mathcal A_U+\mathcal A_V .
 \tag{DN4}
\]
The bounds apply to all limiting profiles with the usual fixed
truncation and vanishing extraction losses. In particular, the
atomic terms vanish if both free letters are entirely nonatomic.

To see the bound on the atomic terms, fix a chart and multiply
the nonzero coefficients of \(F\) to form an integer \(J_F\).
A prime on a retained root with degree \(e\) is on a single free
letter by primitivity. At that prime the valuation of a monomial
\(z^s v^j\), divided by the positive root valuation and with the
letter sign extracted, is \(i_W(e)s+(e-qi_W(e))j\). If the prime
does not divide \(J_F\), a minimum must tie in the original
signed valuations. For each \(j\) just one extremal \(s\) can
participate, because \(i_W(e)>0\), and the choice of extremum
depends only on the letter sign. A tie must occur for two
different wider indices \(j\). Each pair of wider indices
specifies one scalar value of \(e/i_W(e)-q\); since
\(1\le i_W(e)\le I_0\), there are at most
\(I_0\binom{M+1}{2}\) exceptional exact degrees on each letter.
In the rectangle, all have
\(\lvert e/i_W(e)-q\rvert\le\lfloor r/L\rfloor\).

At every other retained position, the prime support lies in
\(J_F\). There are only subpower choices for the integers on that
support, uniformly over the polynomial-height chart coefficients
by the fixed-support bound. The degrees of the exceptional
positions depend on the chart, but their total normalized
mass has upper limit at most the sum of the largest indicated
atoms. Indeed, pass to a subsequence on which the finitely many
exceptional scaled degrees converge in the bounded retained
range. The upper masses of arbitrarily small closed neighborhoods
of those limits shrink to their atomic masses. If several
positions converge to one atom, all their mass together is
bounded by that atom. Fix the actual exceptional roots at this
cost. The full free terms are then determined and recover the
triple. Summing all shallow chart records proves (DN4).

### An entirely diffuse sector beyond terminal lifting

Take all three measures equal to
\(\frac13\operatorname{Unif}[3-3/40,3+3/40]\), with numerical
heights one. On the free pair retain all roots and group at
\(a=3\) with multiplier one. Their symmetric positive and
negative deviations give \(h=1/3\) and \(E=1/80\) for the rational
ratio. The supplying letter has diffuse mass \(1/6\) strictly
above three. Freeze a full band block of mass between \(0.144\)
and \(0.145\) with degrees above three. Use its radical as the
shallow modulus, with \(b=0.144\) and \(\kappa=3\).
Take \(\alpha=99/100,\ M=8,\ r=1/3\). Then
\[
 y=297/100<3,\qquad B=19/90<\tfrac12(0.144)y.
\]
The two atom charges vanish. The bound from (DN4) is at most
\(0.290<1/3\), so a neighborhood satisfying the same strict height
bounds and sufficiently small total atomic mass on the free pair
is excluded, including profiles away from minimum total diffusion.

For the same full free grouping, no rectangle satisfying the old
terminal condition of (QD2) could give a strict saving by choosing
another full or shallow modulus, a fixed amplification, or different \(\alpha,M,r\).
In fact set \(z_0=2h/(\alpha b(M+1))\). The first determinant needs
\(0<z_0<1\). Combining its bound with
\(2B/(1-\alpha)<w\), and optimizing in \(y\), gives
\[
 EM<\tfrac12(1-\alpha)w(1-z_0),\qquad
 hE<\frac{wb}{64(1-b/(2h))}.
\]
For the second inequality use
\(M+1=2h/(\alpha b z_0)\) and
\(\alpha(1-\alpha),z_0(1-z_0)\le1/4\).
If the rectangle is amplified \(\ell\) times, apply this to
\(\ell h,\ell E\) and divide by \(\ell^2\). The same upper bound
holds with \(b/\ell\) in place of \(b\), after dropping a factor
\(1/\ell\). In either case write \(b\) for the effective
shallow cost. The condition \(\mu+b<1/3\) gives \(b<1/3\), and the
degree support gives \(w\le(123/40)\mu\). Thus the right side
is at most \(41/15360<1/240=hE\), a contradiction.
This comparison is for the specified free grouping.

### A broad all-diffuse profile beyond the one-root criterion

There is a barrier even after dropping the terminal requirement.
Give every letter the profile \(\nu=\frac13\operatorname{Unif}[2,4]\).
Consider a one-root projective determinant using the same retained
submeasure and multiplier assignment on both free letters.
Allow arbitrary positive integer multipliers and arbitrary horizontal
slice shapes. Let \(s\) be the retained mass per free letter. A
saving by (DN4) would require
\[
 \mu+b<g:=2s-1/3,\qquad 1/6<s\le1/3.
\]
The modulus weight satisfies \(\kappa b\le w\le4\mu\).
Write \(J\) for the gap between the extreme wider slices and \(R\)
for the largest within-slice width. Let \(h=\int i\,d\nu\) on the
retained mass, and \(E=\int |r-a'i(r)|\,d\nu\) there, where \(a'\)
is the shifted base induced by two monomials in the extreme
slices. The volume \(\mathcal A\) and height \(B\) satisfy
\[
 \mathcal A\le(J+1)R,\qquad B\ge hR,\qquad B\ge JE\quad(J>0).
\]
As in (PS1), the necessary Taylor inequalities are
\(B<b\alpha\mathcal A/2\) and
\(\alpha\mathcal A<\kappa\).
Since \(h\ge s\ge g>b\), these imply \(J\ge2\), and then
\[
 b>\frac{2s}{J+1},\qquad
 E<\frac{2\mu}{J}\le\frac{s}{3},\qquad
 Eh<\mu b(1+1/J)\le\frac{g^2}{3}.
 \tag{DN5}
\]
For the inequality \(E\le s/3\) use
\(\mu<g-b\le s-b\) and \(2(J-1)/(J(J+1))\le1/3\).
For the last inequality, if \(J=2\) then \(b>2g/3\) and
\(\mu b\le (g-b)b\le2g^2/9\). For \(J\ge3\), use
\(\mu b\le g^2/4\) and \(1+1/J\le4/3\).

If \(a'\le0\), positivity of \(i(r)\) gives \(E\ge2s\).
For \(a'\ge2\), folding the interval of length two modulo \(a'\)
gives numerical distance density at most \(1/3\), so \(E\ge3s^2/2\).
With \(h\ge s\) this gives \(Eh\ge g^2/3\). Indeed for
\(s\ge2/9\) use \(g\le s\), and for \(s\le2/9\) use \(g\le s/2\)
and \(s>1/6\).

For \(0<a'\le2\), folding scaled phases has density at most
\((2+a')/3\), so \(E\ge(3/8)a's^2\). The moment of the retained
submeasure is at least \(2s+3s^2\) by the density \(1/6\) on
\([2,4]\). Thus \(a'h+E\ge2s+3s^2\), and \(E\le s/3\) from
(DN5) gives
\[
 Eh\ge\left(\frac{5s}{8}+\frac{9s^2}{8}\right)s^2\ge\frac{g^2}{3}.
\]
For the last inequality put \(t=3s\). The difference after
multiplying by the positive denominator has sign given by
\[
 3t^4+5t^3-32t^2+32t-8
 =(1-t)(-3t^3-8t^2+24t-8)\ge0 .
\]
The last factor is concave on \([1/2,1]\) and positive at its
endpoints. This contradicts (DN5) in all ranges of the shifted
base. The profile remains compatible with this one-root test
throughout the identical-assignment class, even with no atomic
penalty and no terminal lifting requirement.

There is also a chart-budget obstruction to using a direct
multi-group specialization on two entirely diffuse letters.
Put the retained mass \(t_i\) in group \(i\) across the free pair,
so \(\sum_i t_i\le2/3\). In the direct recovery scheme after
specializing \(m-1\) variables, the only free radical cost that can
be saved on a surviving specialization is \(t_1\), and the
exceptional branch saves at most the smallest \(t_i\) with \(i\ge2\)
using the polynomial zero count. Even without any chart or modulus
charge, certifying a strict saving from \(1/3\) by both bounds
would require \(t_i>1/3\) for every \(i\). This is incompatible
with \(\sum_i t_i\le2/3\) for \(m\ge2\). A method counting such
profiles by several group variables needs further arithmetic
information about the exceptional specializations.

## Two first residue determinants and the exact additive relation

Let \(A+B=C\) be primitive positive triples in limiting small-budget boxes.
Fix some root factors, of total normalized cost \(f\), one disjoint full
block \(n_1\mid C\), and another disjoint full block \(n_2\mid B\).
Include the two block costs in \(f\). Once the factors are frozen, express
the two complete odds ratios as
\[
 X_1=A/B=d_1v_1z_1^{q_1},\qquad
 X_2=A/C=d_2v_2z_2^{q_2},\qquad
 X_2=\frac{X_1}{1+X_1}.
\]
Here the constants \(d_i\) are fixed positive rationals from the frozen
factors, and \(q_iL\to a_i>0\) are the two grouped exponents.
On each side group the retained factors at its base \(a_i\)
and collect the exponent corrections in \(v_i\).
At the primes of \(n_i\), assume \(z_i,v_i,d_i\) are units. This is
automatic for groups built from the two original free letters,
by primitivity. Since \(X_1\equiv-1\bmod n_1\) and
\(X_2\equiv1\bmod n_2\), both are one-group congruences as in (QD1).

On side \(i\), assume \(H(z_i)\le X^{(h_i+o(1))L}\),
\(H(v_i)\le X^{E_i+o(1)}\), and choose a shallow record
\(Q_i\) of normalized weight \(b_{Q_i}\ge b_i\) at fullness
parameter \(k_iL\to\kappa_i\). Take fixed \(M_i\ge1\),
\(0<\alpha_i<1\), and \(r_i>0\), and set
\(B_i=h_i r_i+E_i M_i\). Suppose
\[
 \alpha_i(M_i+1)r_i<\kappa_i,\qquad
 \frac{b_i\alpha_i(M_i+1)r_i}{2}>B_i,\qquad
 0<r_i<a_i \qquad(i=1,2),
 \tag{FD}
\]
and require
\[
 \min\left\{
   \frac{4r_1}{a_1-r_1}+\frac{2r_2}{a_2-r_2},
   \frac{2r_1}{a_1-r_1}+\frac{4r_2}{a_2-r_2}
 \right\}<1 .
 \tag{FD1}
\]
Then the normalized count is at most \(f+b_{Q_1}+b_{Q_2}\).
The proof combines first determinants with the exact additive
relation and allows equal-height profiles.

### Initial auxiliary equations and constant coordinates

Apply the first determinant of (QD2) independently on each
shallow residue chart. Its proof uses the two Taylor inequalities in (FD). We obtain nonzero polynomials
\(F_i(z_i,v_i)\) of degrees at most \(R_i=\lfloor r_i/L\rfloor\)
and \(M_i\) in their respective arguments. The coefficients
depend on the charts and frozen values. We only need the
degree bounds here. Pass to their geometric irreducible
factors and the irreducible components of their product
intersected with the real exact ratio relation, in the
algebraic torus. At each fixed \(L\) all degrees and the
number of these components are bounded in terms of the
fixed degree lists and \(q_1,q_2\), uniformly in the coefficients
by Bézout. A product of two plane curves has dimension two.
On a product component where the ratio equation holds
identically, independence of its two coordinates implies
both exact ratios are constant. Such a component contributes
at most one complete triple, because \(A/B\) determines its
primitive numerator and denominator. The same conclusion
holds wherever \(X_1\) is constant. In all other cases the
equation cuts to curves and points.

If \(z_i\) is constant on a component, the original congruence
determines \(v_i\bmod n_i\). The two Taylor conditions in
(FD), and \(b_i\kappa_i\le\log_X n_i+o(1)\), imply
\[
 E_iM_i\le B_i<b_i\kappa_i/2
               \le (\log_X n_i+o(1))/2 .
\]
Thus rational reconstruction determines \(v_i\), and its
complete ratio again determines the triple. This treats
all constant narrow-coordinate branches. A zero-dimensional
component contributes only boundedly many complete ratios.

### Degree of the narrow projection

Work on a remaining smooth projective normalization
\(\widetilde D\) of an irreducible component over
\(\overline{\mathbb Q}\). Let \(H_D(t)\) be the degree
of a nonconstant function \(t\) on this curve, and let it
be zero on constants. Write \(D_i\) for the normalization
of the corresponding irreducible plane projection in
\((z_i,v_i)\), and \(P\) for the normalization of the
narrow plane projection in \((z_1,z_2)\). The latter is
a curve since both narrow coordinates vary. The former
may have \(v_i\) constant. Denote their actual bidegrees
by \(R_i',M_i'\), with \(M_i'\ge1\), \(R_i'\le R_i\),
\(M_i'\le M_i\). The ordinary degree of the projection
to \(z_i\) is \(M_i'\). Hence on \(\widetilde D\)
\[
 H_D(v_i)=\frac{R_i'}{M_i'}H_D(z_i)\le R_iH_D(z_i).
\]
If \(v_i\) is constant this formula simply gives zero.
Let \(D_X=H_D(X_1)=H_D(X_2)\); the equality is Möbius
invariance, and multiplying by a constant has no effect.
By the height inequality for products,
\[
 D_X\ge(q_i-R_i)H_D(z_i).
 \tag{FDH}
\]

Suppose the degree of the narrow plane curve is at most
some fixed \(D_0\). The function field of \(\widetilde D\)
is the compositum of \(K(P)=\overline{\mathbb Q}(z_1,z_2)\)
and \(K(D_1)=\overline{\mathbb Q}(z_1,v_1)\) over
\(\overline{\mathbb Q}(z_1)\). Indeed the exact additive
relation expresses \(v_2\) in terms of \(z_1,z_2,v_1\).
Thus \([\widetilde D:P]\le M_1'\). In characteristic zero
the ramification divisor of \(\widetilde D\to P\)
has degree at most \([\widetilde D:D_1]\) times the degree
of the ramification divisor of \(D_1\to \mathbf P^1_{z_1}\).
For clarity, at a point above two local indices \(e_1,e_2\),
the compositum local index is
\(e=\operatorname{lcm}(e_1,e_2)\), since the complete
residue fields are algebraically closed and characteristic
zero. Its contribution \(e/e_2-1\) to ramification over
\(P\) is at most \((e/e_1)(e_1-1)\), by
\(e/e_1+e/e_2\le e+1\). Sum this over all points
using the pullback degree for the ramification on \(D_1\).

The plane degree of \(P\) bounds its geometric genus in
terms of \(D_0\). Also \(g(D_1)\le(R_1'-1)(M_1'-1)\)
by the arithmetic genus in \(\mathbf P^1\times\mathbf P^1\)
when \(R_1'>0\). When \(R_1'=0\), \(v_1\) is constant,
\(M_1'=1\), and the \(z_1\)-map has no ramification.
In all cases the ramification degree of the \(z_1\)-map is
at most \(2R_1'M_1'\). Riemann–Hurwitz therefore gives
\[
 2g(\widetilde D)-2 \le O_{D_0,M_1}(1)+
             2R_1H_D(z_1).
\]
Apply Riemann–Hurwitz again, now to the function \(X_1\)
and just its fibers over \(0,-1,\infty\).
The function-field abc inequality proved this way is
\(D_X\le2g(\widetilde D)-2+|S|\), where \(S\) is the
union of these fibers. Since \(1+X_1=X_1/X_2\),
\(S\) is contained in the union of the zeros and poles
of \(z_1,v_1,z_2,v_2\), with count at most
\[
 2(R_1+1)H_D(z_1)+2(R_2+1)H_D(z_2).
\]
Consequently
\[
 D_X\le O_{D_0,M_1}(1)
       +(4R_1+2)H_D(z_1)+(2R_2+2)H_D(z_2).
\]
Interchanging the two plane covers gives the same bound with the
coefficients \(4R_i+2\) and \(2R_i+2\) assigned in the reverse order.
By (FD1) one of these bounds, together with (FDH), gives
\(D_X\le O_{D_0,M_1,M_2}(1)+(1-\epsilon)D_X\) with
\(\epsilon>0\). This contradicts \(D_X\ge q_1-R_1\to\infty\).
Thus the degrees of all nonconstant narrow plane
projections tend to infinity uniformly in the
coefficients.

### Counting the large-degree narrow projection

For any prescribed normalized tolerance \(\epsilon>0\),
first take a fixed lower projection degree \(D_0\) large
in terms of \(\epsilon,h_1,h_2\). Then take \(L\) small
enough that the preceding argument excludes this degree
range. The projective height of the narrow point
\([1:z_1:z_2]\) is at most \(X^{(h_1+h_2+o(1))L}\).
Heath-Brown's rectangular projective curve estimate, with
all sides equal, gives the upper exponent \(3/\deg P\) plus
an arbitrary allowance in this projective height. Choose
\(D_0\) so that it counts the narrow points at normalized cost
less than \(\epsilon\), allowing fixed-truncation height losses. The algebraic
curves may be taken over \(\overline{\mathbb Q}\).
If a component is not defined over \(\mathbb Q\), its
rational points lie in the intersection with a distinct
conjugate and contribute only boundedly many by
Bézout. The degree is bounded at fixed \(L\) uniformly
in the coefficients, so all curve constants and the
component count can be absorbed in the inner height
loss. [Heath-Brown, Theorem 3](sources/heath-brown-curves-math0405392v1.pdf)

Given \(z_1,z_2\), the nonvertical \(F_1\) gives at most
\(M_1\) choices of \(v_1\). Its specialization in \(z_1\)
cannot vanish identically, since that would give a
vertical factor of the chosen irreducible component.
Each choice determines \(X_1\) and thus the primitive
triple. Summing over the shallow charts and the frozen
values gives the asserted bound.

### Equal-height diffuse profiles with two locks

Give \(A\) the measure \(\frac13\operatorname{Unif}[3-3/2500,3+3/2500]\)
and each of \(B,C\) the measure
\(\frac13\operatorname{Unif}[3-39/500,3+39/500]\). All first moments
are one and all three measures are nonatomic. Freeze mass
\(\mu=33/400\) from the two extreme ends on each of \(B,C\), with
equal portions on either end. Use these two factors as the moduli
and their radicals as the shallow divisors, so \(b_1=b_2=\mu\).
Take \(\kappa_1=\kappa_2=29/10\), below the least selected degree.
Group both free ratios at \(a_i=3\). The narrow heights are
\(h_i=1/3\). The equal masses on both sides of each deviation
give the wider height bounds
\[
 E_i=\frac{3/2500}{12}
      +\frac{39/500}{12}(1-3\mu)^2
      =0.003780665625 .
\]
For \(\alpha_i=99/100,\ M_i=16,\ r_i=17/100\), the Taylor order
is \(2.8611\) and \(B_i<0.117158<0.118020375=b_i y_i/2\).
The slope condition (FD1) holds since \(6r_i/(3-r_i)<1\).
The total freeze and chart cost is \(4\mu=33/100<1/3\).
Thus the two-curve criterion gives a strict margin in an
equal-height fully diffuse sector.


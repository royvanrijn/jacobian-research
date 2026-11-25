# Small roots and a larger cofactor in abc congruences

A congruence can leave a cofactor whose height is much larger than the
heights of the other roots. Bounding its degree in the auxiliary
polynomials keeps the coefficients small enough for local lifting. The
resulting counts extend the prime-root graph and constrain the distribution
of scaled degrees at the limiting one-third exponent.

## A full-modulus count in three variables

Fix \(C>0\) and a sequence \(\rho_k\to0\) of nonnegative real numbers. For \(k\)-full \(n\), consider positive unit rationals \(x,y,V\) satisfying
\[
Vx^a y^b\equiv K\pmod n,\qquad
H(x),H(y)\le n^{C/k},\qquad H(V)\le n^{\rho_k},
\tag{AW1}
\]
where \(a,b\) are fixed signed integer exponents and \(K\) is a unit. Suppose the exact rational monomial fibers have size at most \(F_0\). For sufficiently large \(k\), independently of \(a,b\), the count is
\[
O_{C,k,a,b,\varepsilon}\!\left(F_0 n^{3/k+\varepsilon}\right).
\tag{AW2}
\]
The constants are uniform in the unit target, and the threshold in \(k\) may depend on the height sequence. The proof uses the same \(T_r\) charts, good-prime Hensel records, and uniform coefficient-height conventions as the [bounded-variable full-modulus lemma (FC1)](ABC_OMITTED_COORDINATE_GRAPH.md#signed-monomial-congruences-over-full-moduli). In this case an asymmetric determinant and a rank defect give a coefficient bound small enough for the wider variable.

## An asymmetric polynomial with controlled coefficients

Take a large fixed constant \(M_0=M_0(C)\), and take
\[
R=\left\lceil M_0\rho_k k+\sqrt{k}\right\rceil.
\]
Then \(R/k\to0\), and we can use polynomial supports whose degrees in \(x,y\) are at most \(R\) and whose degree in \(V\) is a suitably chosen constant at most \(M_0\). All coefficient-height exponents below have the form \(O_C(R/k)\). Constants in degrees and the \(O_C\)-bounds are independent of \(a,b\) and of \(R\).

First lift to the exact surface \(\phi=K\), with \(\phi=Vx^a y^b\), at each \(p^e\Vert n\), by correcting \(V\) while keeping \(x,y\). This preserves the original residues to precision \(e\) because the derivative in \(V\) is a unit. Partition the \(T_r\) charts of dimension two modulo \(p^{\lfloor e/k\rfloor}\), as in (FC1), with \(r=O_C(R)\). This costs \(O(n^{2/k+\varepsilon})\) classes. Old and corrected points give the same Taylor approximations at the determinant precision, since \(r<k\) for large \(k\).

Use the monomials \(x^i y^j V^\ell\) with \(0\le i,j\le R\) and \(0\le\ell\le M\), where \(M\) is a fixed constant to be chosen. There are \(D=(R+1)^2(M+1)\) columns. Multiply each evaluation row by its \(x\)- and \(y\)-denominators to the power \(R\) and its \(V\)-denominator to the power \(M\). The resulting entries are integers of height at most \(n^{(2C+1)R/k}\) if \(M\le M_0\). These row multipliers are units at the modulus primes.

For \(t=\lceil A R^2\rceil\) with fixed \(A\) sufficiently large in terms of \(C\), the sum of the least orders of \(t\) Taylor monomials in two variables exceeds \(2t(2C+2)R\). Their largest order is \(O_C(R)\), so the chart precision can be fixed accordingly. Any \(t\)-by-\(t\) evaluation minor in a class vanishes: its p-adic divisibility over all modulus primes exceeds its real height, using \(\lfloor e/k\rfloor\ge e/(2k)\). Choose \(M\ge 2A+2\) and choose \(M_0\) at least this large. Since \(t\le D/2\), rank is at most \(D/2\). Therefore one can choose a nonzero integral kernel vector with coefficient height \(n^{O_C(R/k)}\). For example, map the box of integer coefficient vectors of radius \(Q\) to the evaluations on a maximal independent row set of size \(r_0\le D/2\). The image has at most \((O_D(Q)n^{(2C+1)R/k})^{r_0}\) points. With \(Q=n^{O_C(R/k)}\) and a large enough absolute multiplier in that exponent, the domain is larger and subtraction gives the vector. The chart count has thus produced an auxiliary polynomial \(G(x,y,V)\) of degree \(O(R)\) in the narrow coordinates and degree \(O(1)\) in \(V\), with coefficient height \(n^{O_C(R/k)}\). Every primitive rational factor of \(G\) has coefficient height at most \(c_{R,M}H(G)\). For instance, make a Kronecker substitution with base larger than the degree in each variable, so exponents in the polynomial and each factor stay distinct, and use the univariate factor-height estimate. The multiplicative constant is independent of the coefficient height. Constants depending on \(R,M,k\) are fixed before \(n\) tends to infinity.

## A second auxiliary and the deficient rank locus

Work on one primitive irreducible rational factor at a time, and denote it by \(G\). On a geometric component of \(G=0\) where \(\phi\) is constant, use \(F_0\). Otherwise choose a nonzero maximal minor \(J\) of the rows \(\nabla G\) and \((a yV,bxV,xy)\) that is coprime to the irreducible rational factor \(G\). All geometric components of a rational irreducible factor have the same constant/nonconstant status, and a nonzero minor on one is generically nonzero on each conjugate. The new polynomial has narrow degree \(O(R)\), wide degree \(O(1)\), and coefficient height \(n^{O_C(R/k)}\). At a rational tuple where it is nonzero its numerator is at most \(n^{O_C(R/k)}\), because the narrow coordinates cost \(n^{C/k}\) each and the wide coordinate occurs only to bounded degree.

If \(G\) depends on \(V\), the locus \(G=J=0\) projects to the narrow plane equation \(A(x,y)=\operatorname{Res}_V(G,J)=0\). The resultant is nonzero by coprimality over \(\mathbb Q(x,y)\). Its Sylvester matrix size is bounded because both wide degrees are bounded. Thus \(A\) has narrow degree \(O(R)\) and height \(n^{O_C(R/k)}\). The locus is handled by the terminal count below with the equations \(A=G=0\). If \(G\) is independent of \(V\), take \(J=xyG_x\) or \(xyG_y\), choosing a derivative not divisible by \(G\). The deficient set has only boundedly many narrow pairs by Bézout. For any fixed narrow pair the congruence determines at most one wide rational coordinate, since \(H(V)^2<n\) for large \(k\).

On \(G=0\) with \(J\ne0\), discard the primes with \(v_p(J)>e/3\). Their prime powers cost total log at most \(O_C(R/k)\log n\). The numerator height just proved applies because coordinates and denominators are units at the modulus primes. At every other prime, lift to \(G=0,\phi=K\) with nonzero rank and precision at least \(2e/3\). The corrected set is a curve in the three-variable box. Uniform \(T_r\) charts in one variable and residues modulo \(p^{\lfloor e/k\rfloor}\) cost \(O(n^{1/k+\varepsilon})\) records, including the bad-prime marker. The same fixed-prime charts as in (FC1) handle the finitely many small primes. Lift errors do not affect the determinant when \(r=O_C(R)=o(k)\).

If \(G\) depends on \(V\), construct a second polynomial \(B(x,y)\) in the two narrow coordinates of bounded degree and height \(n^{O_C(1/k)}\): the usual one-dimensional chart determinant suffices on a sufficiently large fixed square monomial support. It is coprime to \(G\). These children are handled by the terminal count, with \(A=B\) the narrow equation and \(G\) the remaining one.

If \(G\) depends only on \(x,y\), omit a coordinate on which it depends, say \(y\), in the next determinant. Use the monomials \(x^i V^\ell\) with \(0\le i\le R\), \(0\le\ell\le M'\). Force rank at most half the number of these columns by the same argument, now with chart dimension one and minors of size \(\lceil A'R\rceil\), for fixed large \(A'\). The sum of their least Taylor orders is quadratic in the minor size. Good primes contribute at least half the logarithmic prime-power mass once \(k\) is large, so a larger fixed \(A'\) suffices. Take \(M'\ge2A'+2\), and in the first choice of \(M_0\) include this fixed bound as well. The resulting polynomial \(B(x,V)\) has narrow degree \(O(R)\), bounded wide degree and coefficient height \(n^{O_C(R/k)}\), by the same kernel pigeonhole argument. It is coprime to \(G\). If it is independent of \(V\), there are boundedly many narrow pairs and modular uniqueness for \(V\). Otherwise use \(A=G,B=B(x,V)\) in the terminal count.

## The terminal count without eliminating a narrow coordinate

In every terminal branch we have two coprime equations
\[
A(x,y)=0,\qquad B(x,y,V)=0
\]
of narrow degrees \(O(R)\) and heights \(n^{O_C(R/k)}\), with \(B\) of bounded wide degree. Take a primitive irreducible rational factor of \(A\). Its singular points give boundedly many narrow pairs. Since \(A\) and \(B\) are coprime, one coefficient of \(B\) as a polynomial in \(V\) is nonzero modulo the selected factor of \(A\). Where every such coefficient vanishes there are boundedly many narrow pairs by Bézout, and \(V\) is unique by the modulus.

At each other point choose the order \(j\) of the root in the polynomial \(B(x,y,\cdot)\) at \(V\), so \(1\le j\le \deg_V B=O(1)\). The point lies on \(A=0,\ B_j=\partial_V^{j-1}B=0\), and \(\partial_V B_j\ne0\). These two equations have independent gradients because \(\nabla A\ne0\). In this open stratum the components are curves. If \(\phi\) is constant on a geometric component, it gives at most \(F_0\) points. Otherwise the logarithmic Jacobian for \(A,B_j,\phi\) vanishes on only boundedly many points of it. Components outside the open stratum do not enter that argument, and the bound may depend on \(R\) at fixed \(k\).

At every other solution, the nonzero logarithmic Jacobian numerator has height \(n^{O_C(R/k)}\), because narrow degrees and coefficient-height exponents have those bounds and all wide degrees remain bounded. Discard its high-valuation primes at log cost \(O_C(R/k)\log n\), as above. At good primes Hensel lifts to isolated exact common zeros of \(A,B_j,\phi-K\), of which there are boundedly many at each prime by Bézout, uniformly in the target. A record of these isolated lifts and bad primes costs \(O(n^\varepsilon)\). Congruence of two original tuples in one record holds modulo a product of logarithm at least \((1/2-o(1))\log n\), using the corrected precision \(2e/3\). Since \(C/k,\rho_k\to0\), coordinatewise rational reconstruction is unique. This proves (AW2) after the two chart partitions of costs \(n^{2/k+\varepsilon}\) and \(n^{1/k+\varepsilon}\). The exceptional loci use no more chart partitions.

## Amplification when the cofactor shares primes with the roots

Assume now that \(x,y\) have disjoint prime supports and \(|a|,|b|\gg k\). The cofactor \(V\) may share primes with \(x\) or \(y\). For any fixed exact value \(W=Vx^a y^b\), the number of triples in that fiber is at most
\[
O_{k,a,b,\varepsilon}
  \bigl(n^{\,2\rho_k/|a|+2\rho_k/|b|+\varepsilon}\bigr).
\tag{AW3}
\]
To prove this, invert either narrow coordinate as necessary so that \(a,b>0\). Write \(x=x_+/x_-\), \(y=y_+/y_-\), with lowest terms, and similarly \(W=W_+/W_-\) and \(V=V_+/V_-\). For a prime \(p\mid x_+\),
\[
a v_p(x_+)\le v_p(W_+)+v_p(V_-)
\]
because \(p\) divides neither \(x_-\) nor \(y_+y_-\). Split \(x_+\) into the part with exponent \(\min(v_p(x_+),\lceil v_p(W_+)/a\rceil)\) at each \(p\), and the remaining part. The first part divides \(W_+\), so there are divisor-many choices. The \(a\)-th power of the second part divides \(V_-\), so that part is at most \(n^{\rho_k/a}\). Use \(W_-,V_+\) for \(x_-\) and do the corresponding two counts for \(y\). For each \(x,y\) the exact value determines \(V\), and \(H(W)\le n^{O_{k,a,b}(1)}\) allows the divisor factors to be absorbed into \(n^\varepsilon\). This also applies when either narrow root equals one.

Fix an amplification order \(\ell\). Assign each prime one of four colors for the numerator or denominator of \(x\) or \(y\), and independently one of two colors for the numerator or denominator of \(V\). Some coloring retains \(n^{-o(1)}\) of the triples, since the numerator and denominator of each variable have \(o(\log n)\) distinct prime divisors. In the retained \(\ell\)-fold coordinate products, numerator and denominator do not cancel, and the two narrow coordinates have disjoint prime supports. There are divisor-many preimages of each product tuple. Apply (AW2) and (AW3) to the products with height constants \(\ell C,\ell\rho_k\) and target \(K^\ell\). Taking \(\ell\)-th roots bounds the number of original triples by
\[
O_{C,\ell,k,a,b,\varepsilon}
 \bigl(n^{\,3/(\ell k)+2\rho_k/|a|+2\rho_k/|b|+\varepsilon}\bigr).
\tag{AW4}
\]
All constants, exponents and the amplification order are fixed before \(n\to\infty\).

In an anatomic small-budget box, let \(n\) be a full modulus on a locked letter with positive limiting numerical weight and fullness threshold \(k\asymp 1/L\), where \(L\) is the total logarithmic root mass. A fixed number of roots of positive limiting normalized mass has size \(X^{O(L)}\). Each positive-weight root among them has degree comparable to \(1/L\). Thus (AW4) releases two such roots and a cofactor of numerical height \(X^{o(1)}\) in their signed congruence at cost \(o(L)\) in the \(X\)-exponent. The cofactor may be a ratio of products of low factors, and frozen factors can be absorbed into the unit target. When the exact value is the signed ratio of the complete free terms, primitivity recovers both terms from the ratio and then the original equation recovers the third. At fixed truncation each recovered integer has only divisor-many anatomic representations.

## A four-prime graph count

Use the full-modulus counts (FC2) and (AW4) in the [prime-root graph](ABC_EXCEPTIONAL_SET_AND_POWERFUL_TRIPLES.md#a-graph-count-with-two-or-three-prime-roots). In the applications below, first fix a small-budget box and its anatomic truncation, and take the height parameter \(X\) arbitrarily large. All degrees and amplification orders for that box are then fixed. After bounding the box exponent, pass to a sequence of such boxes with \(L\to0\). Errors \(o(L)\) and \(o(1)\) in height and counting exponents refer to this order of limits. Arbitrarily small losses in the \(X\)-exponent at each fixed truncation will be implicit.

Consider a modulus \(n=X^{w+o(1)}\) with a strictly positive limiting weight \(w\), which is \(k\)-full for \(k\asymp 1/L\). At fixed \(n\), consider tuples satisfying
\[
 U\prod_{i=1}^{4}p_i^{a_i}\equiv K V\pmod n .
 \tag{PG4}
\]
All variables and \(K\) are units. The \(p_i\) are distinct primes, each is at most \(X^{O(L)}\), and \(|a_i|\asymp 1/L\). The positive integer cofactors \(U,V\) have height \(X^{o(1)}\) and range over at most \(X^{(f+o(1))L}\) and \(X^{(g+o(1))L}\) values, respectively. Assume each moving cofactor is coprime to all four primes of its own tuple. In the application this follows from the anatomic prime supports after freezing the exceptional coefficients and the small complementary factors of the four roots.

The count of tuples in (PG4), at fixed \(n\), is
\[
 \ll X^{(G(f,g)+o(1))L},\qquad
 G(f,g)=\max\bigl(f,g,\min(f+g/2,g+f/2)\bigr).
 \tag{PG4a}
\]
Indeed, color the primes with four colors, and retain only tuples in which \(p_i\) has color \(i\). Some coloring retains at least \(4^{-4}\) of them. Build a simple bipartite graph on the values \(U,V\), with one prime witness chosen for each edge. For fixed \(U,V\), (FC2) on the four integer prime coordinates gives at most \(X^{o(L)}\) witnesses. Thus it suffices to bound the number of edges \(E\) in the retained graph.

Take ordered two-edge paths with middle vertex \(V\) and distinct endpoints \(U,U'\). Call the path bad if its witnesses share at least two corresponding prime coordinates. Starting with one edge, there are at most six choices for two such coordinates to fix. The second witness has at most two remaining prime coordinates, together with its moving cofactor \(U'\). They have disjoint prime supports within each tuple. Apply (AW4), using a constant coordinate equal to \(1\) if fewer than two prime coordinates are free. The target absorbs the fixed primes and \(V\). This counts at most \(X^{o(L)}\) second edges per first edge. The height hypotheses are uniform because the prime coordinates have height \(n^{O(1/k)}\), the moving cofactor has height \(n^{o(1)}\), and the modulus retains positive limiting weight.

For a good path, fix the two distinct endpoints \(U,U'\). Eliminating the middle cofactor gives
\[
 \prod_{i=1}^{4}(p_i/p_i')^{a_i}\equiv U'/U\pmod n .
\]
A nontrivial quotient is a reduced ratio of two distinct primes, and a trivial quotient is recorded as \(1/1\). The eight integer numerator and denominator coordinates are pairwise coprime, since different prime coordinate positions have different colors. Apply (FC2) to these eight integer coordinates, with signed exponents \(a_i,-a_i\), to count the quotient tuples at cost \(X^{o(L)}\). Every nontrivial quotient determines the two ordered primes. There is at most one common prime in a good path. With one common prime, the first congruence counts this prime and the middle cofactor \(V\) by (AW4) with a constant coordinate equal to \(1\), again at cost \(X^{o(L)}\). With no common prime, the first congruence determines \(V\) modulo \(n\), and hence uniquely because \(V<n\) for sufficiently small budgets and sufficiently large \(X\). Each distinct endpoint pair therefore contributes at most \(X^{o(L)}\) good paths.

All these uses of (FC2) and (AW4) involve a bounded number of coordinates. More explicitly, for each desired normalized loss first choose a sufficiently large fixed amplification order, and then take the budgets sufficiently small. The resulting constants depend on the fixed truncation and on the exponents, while the bounds are uniform in the residue target and the frozen coefficients. This also applies with \(U,V\) interchanged. Accounting for the diagonal paths, the two sums of squared degrees satisfy
\[
 \sum_V\deg(V)^2\ll X^{o(L)}\bigl(|E|+X^{(2f+o(1))L}\bigr),
 \qquad
 \sum_U\deg(U)^2\ll X^{o(L)}\bigl(|E|+X^{(2g+o(1))L}\bigr).
\]
Cauchy–Schwarz first gives the edge exponent \(\max(g,f+g/2)\), and the other side gives \(\max(f,g+f/2)\). Taking the smaller of the two and reinstating witness multiplicity proves (PG4a). The proof counts graph vertices by their numerical cofactor values. At fixed truncation, the divisor bound absorbs the multiplicity of their original anatomic representations into the arbitrarily small \(X\)-exponent loss.

## Strict residual mass at the balanced \(1/3\) endpoint

Use the [balanced-limit notation](ABC_OMITTED_COORDINATE_GRAPH.md#necessary-conditions-at-the-balanced-13-endpoint) \(d=l_A=l_B=l_C=1/3\) and \(\Delta=\delta_A+\delta_B+\delta_C\). We prove that every such limiting profile has
\[
 \Delta>1/3.
 \tag{DR}
\]
The argument uses (D1)–(D3), the residual-mass triangle inequalities, and the bounded-root release (30). Recall that a coordinate of positive limiting numerical weight on \(T\) has normalized mass at least
\[
 t_T=1/3-\delta_U-\delta_V=1/3-\Delta+\delta_T
 \tag{DR1}
\]
when \(t_T>0\). This follows from (D1) with that coordinate as modulus. If also \(\omega_T=0\), fixed finite lists of these positive-weight coordinates capture the whole limiting numerical weight of \(T\). To see the finiteness explicitly, use a diagonal subsequence for the selected largest coordinates in the definition of \(\omega_T\), recording coincidences and limits. Each distinct coordinate sequence of positive weight has mass at least \(t_T\), so there are only finitely many of them. The fixed-degree part has numerical weight tending to zero. The definition \(\omega_T=0\) then gives the full numerical weight on this finite list. The bounded number of selected coordinates have degrees tending to infinity. Write \(P_T\) for their total normalized mass. In particular \(\delta_T+P_T\le 1/3\), as all of these coordinates can be removed in the definition of \(\delta_T\).

First suppose that either \(\Delta<1/3\), or \(\Delta=1/3\) and all three \(\delta_T\) are positive. Then all complementary residual sums are below \(1/3\). By (D2) all numerical heights tend to one, and by (D3) every \(\omega_T\) is zero. Thus \(t_T>0\) and the finite lists just described apply to every letter. Each positive-weight coordinate in these lists has normalized mass with a positive limit, so its degree is \(\asymp 1/L\). The product of all power factors outside a letter's list has numerical height \(X^{o(1)}\), after fixing the exceptional anatomic coefficient at the vanishing truncation loss. By (30),
\[
 P_A+P_B+P_C\le 1-d=2/3.
 \tag{DR2}
\]

Suppose a letter \(T\) has just one positive-weight coordinate. Select one positive-weight coordinate on each of the other two letters \(U,V\), with masses \(\beta_U,\beta_V\). Use the selected power on \(U\) as modulus and freeze the other positive-weight coordinates on \(V\), at normalized cost \(P_V-\beta_V\). The ratio of the two free complete terms, up to the frozen factors, is the signed power ratio of their two remaining narrow roots multiplied by the ratio of their low complements. These complements have height \(X^{o(1)}\). The two remaining roots are coprime, of height \(X^{O(L)}\) and degrees \(\asymp1/L\), so (AW4) releases both roots and the low ratio at cost \(o(L)\). Frozen factors are units at this modulus. All triples counted have pairwise coprime complete terms, so the exact ratio recovers both free terms, and their signed equation recovers the modulus letter as well. The representation count at fixed truncation is divisor-many for any factors recovered only through their complete term. Thus
\[
 d\le\beta_U+P_V-\beta_V.
\]
Swapping \(U,V\) gives \(d\le\beta_V+P_U-\beta_U\). Adding yields \(2d\le P_U+P_V\), contradicting (DR2) since \(P_T>0\) and \(d=1/3\). Hence every letter has at least two positive-weight coordinates.

Choose two on each letter. Their total mass is at least
\[
 2\sum_T t_T=2-4\Delta .
\]
This is at most \(2/3\) by (30), so \(\Delta\ge1/3\). If \(\Delta=1/3\) and the \(\delta_T\) are positive, then \(t_T=\delta_T\), and the sum just counted is exactly \(2/3\). Consequently there are exactly two positive-weight coordinates on each letter and both have mass \(\delta_T\). Now \(3\delta_T=\delta_T+P_T\le1/3\) for each letter, so all three residual masses are \(1/9\). This gives six positive-weight coordinates of mass \(1/9\), two per letter. Their complements have numerical height \(X^{o(1)}\) and normalized mass \(1/9\) per letter.

It remains to consider whether \(\Delta=1/3\) can occur with some \(\delta_T=0\). The residual triangle inequalities give \(\delta_U=\delta_V=1/6\) on the other letters. By (D2), if there is a positive limiting short height gap then \(T\) is the short letter, so \(U,V\) both have positive height in every case. By (D3) their residual numerical weights vanish, and (DR1) makes each of their positive atoms cost at least \(1/6\). Since \(\delta_U+P_U\) and \(\delta_V+P_V\) are at most \(1/3\), each of \(U,V\) has exactly one positive atom of mass \(1/6\), and its complement has numerical height \(X^{o(1)}\).

If \(T\) has limiting numerical height zero, lock the positive atom on \(U\) at cost \(1/6\). In the complete ratio \(V/T\), keep the root of the single atom of \(V\), and treat the ratio of its low complement to \(T\) as the wider rational coordinate. It has height \(X^{o(1)}\). Apply (AW4) with the second narrow coordinate equal to \(1\) and its signed exponent chosen as \(-k\). It counts the complete ratio at vanishing normalized cost, contradicting \(d=1/3\). If \(T\) has positive height, apply (AW4) to \(U/V\) against any positive-weight full modulus from \(T\) of fullness order \(1/L\). Each such modulus must cost at least \(1/3\). The high-block argument forces one positive coordinate on \(T\) with full normalized mass \(1/3\). Indeed, the coordinates of degrees at least \(\lfloor h_T/(2L)\rfloor\), for \(h_T>0\) the limiting height of \(T\), have positive total weight. If their largest normalized mass were below \(1/3\) by a fixed margin, greedy grouping with a cap strictly between that upper bound and \(1/3\) would give boundedly many high blocks. Some block would have positive limiting weight and cost less than \(1/3\). Thus a high coordinate with limiting mass \(1/3\) exists. If its complement retained positive numerical weight, the same high cutoff argument on the complement would supply a positive-weight full modulus at cost \(o(1)\) in normalized units, since the complement's normalized mass tends to zero. The sole atom on \(T\) therefore accounts for its full numerical height. Now lock the atom of \(U\) and count \(T/V\) by (AW4), again at total normalized cost \(1/6\), a contradiction.

We have reduced \(\Delta\le1/3\) to the six-atom profile, with all six normalized masses \(1/9\). We now exclude that profile. Fix \(0<\zeta<1/108\), and work in the boxes along the limiting sequence. If one of the six roots has a divisor \(q\) of normalized logarithmic size \(\log q/(L\log X)\) in \([\zeta,1/9-\zeta]\), use the divisor to the coordinate's original power as modulus. The divisor can be enumerated at normalized cost at most \(1/9-\zeta\). Its numerical weight is bounded below by a positive constant depending on \(\zeta\), since the original degree is \(\asymp1/L\). Freeze one of the two positive-weight roots on each other letter, with combined normalized cost \(2/9+o(1)\). The remaining two roots and the ratio of the low complements are counted by (AW4) at cost \(o(L)\). This is uniform in the divisor throughout the displayed interval. The complete ratio recovers the two free terms and then the third term, with the usual divisor bound for representations. A union bound over the six choices shows that these solutions have normalized exponent at most \(1/3-\zeta+o(1)\).

In all remaining solutions, each of the six roots has a prime factor of normalized logarithmic size greater than \(1/9-\zeta\). Otherwise, since there is no divisor in \([\zeta,1/9-\zeta]\), every prime factor would have normalized size below \(\zeta\). Multiplying these factors with multiplicity up to the first partial product of size at least \(\zeta\) produces a divisor in \([\zeta,2\zeta]\), which is contained in the forbidden interval. The full root has size \(1/9+o(1)\), so the indicated large prime has multiplicity one and its complementary factor has size at most \(\zeta+o(1)\). Freeze the six complementary factors at total normalized cost at most \(6\zeta+o(1)\) and absorb their powers into the corresponding fixed coefficients. All the six large primes are distinct by pairwise coprimality of the original root coordinates. Each remains of positive numerical weight. The low cofactors of the letters have height \(X^{o(1)}\) and range over at most \(X^{(1/9+o(1))L}\) values at fixed truncation, by counting their original root-coordinate boxes. The freezing does not introduce any overlap between a low cofactor and a retained large prime. The exceptional coefficients were fixed separately.

Lock one of the six primes to its original degree as modulus, at normalized cost at most \(1/9+o(1)\). Reduce the signed three-term equation modulo this power. The four prime coordinates on the other two letters and their low cofactors satisfy (PG4), with any fixed component factors absorbed into the unit target. Apply (PG4a) with \(f=g=1/9\), giving normalized cost \(1/6\). The four primes and the two cofactors determine the two complete free terms, hence also the third. The resulting normalized exponent, including all freezing, is at most
\[
 6\zeta+1/9+1/6+o(1)=6\zeta+5/18+o(1)<1/3.
\]
Thus both the divisor families and the surviving prime families have a strict margin below the required \(d=1/3\). All auxiliary parameters, including the finite selections and amplification orders, are chosen at the indicated fixed truncations and used with \(X\) sufficiently large before passage to the small-budget limit. Summing the bounded choices and absorbing fixed-truncation divisor losses preserves the two strict margins. This excludes the six-atom profile and proves (DR).

## A bounded number of small roots

Fix \(m\ge1\), \(C>0\), and \(\rho_k\to0\). Let \(n\) be \(k\)-full, and suppose the positive unit rationals \(x_1,\ldots,x_m,V\) satisfy
\[
H(x_i)\le n^{C/k},\qquad H(V)\le n^{\rho_k},\qquad
\phi(\mathbf x,V)=V\prod_{i=1}^m x_i^{a_i}\equiv K\pmod n
\]
for fixed signed integers \(a_i\) and a unit target \(K\). Assume at most \(F_0\) tuples give any one exact rational value of \(\phi\). For sufficiently large \(k\) in terms of \(m,C,\rho_k\), independently of the exponents \(a_i\), the count is
\[
O_{m,C,k,\mathbf a,\varepsilon}\!\left(
F_0n^{(m^2+m)/k+\varepsilon}\right).
\tag{GW}
\]
The proof uses the \(T_r\) chart and lifting records of (FC1), for both large and finitely many fixed small primes. A step that produces narrow equations with large coefficients is followed by a dimension reduction that controls the coefficients using only the short-coordinate sample.

### Preparing the short-coordinate locus

Use induction on the dimension \(j\) of a narrow algebraic set \(Z\subset\mathbb A^m\) containing the short-coordinate tuples in a branch. Its geometric degree is bounded with the fixed counting parameters, with no restriction on its coefficient height. Split into geometric components. The rational points on a component that is not defined over \(\mathbb Q\) lie in its intersection with a distinct conjugate. Taking the union of conjugate intersections gives a rational algebraic locus of smaller dimension. Split geometrically integral \(\mathbb Q\)-components by degree as follows.

If a \(j\)-dimensional component has sufficiently large degree \(D\), use the Salberger covering step from the proof of the bounded-root release (30). A linear projection with coefficients bounded in terms of \(D,m\) maps it finitely and birationally to a projective hypersurface of the same degree. The projected points have height at most \(O_{D,m}(n^{mC/k})\), by clearing the short-coordinate denominators. The hypersurface covering theorem gives \(O(n^{1/k+\varepsilon})\) auxiliary cuts of proper lower dimension when \(D\) exceeds a fixed threshold depending on \(m,C\): choose the threshold with \(mC(j+1)/D^{1/j}<1\). Pull back the auxiliary forms and descend the narrow dimension. This step uses no height bound on \(Z\). [Salberger, Theorem 0.10](https://londmathsoc.onlinelibrary.wiley.com/doi/full/10.1112/plms.12508) · [Castryck–Cluckers–Dittmann–Nguyen, Proposition 4.3.2](https://rcluckers.perso.math.cnrs.fr/prints/Cast-Cl-Ditt-Nguy.pdf) · [Bounded-root release (30)](ABC_EXCEPTIONAL_SET_AND_POWERFUL_TRIPLES.md#releasing-a-bounded-number-of-high-degree-coordinates)

On a component of bounded degree, replace the narrow variety by one of bounded degree and coefficient height \(n^{O(1/k)}\) using only its rational sample. The fact used is that a projective variety of degree at most \(D\) is cut out set-theoretically by forms of degrees at most \(D\). For a point outside a geometric irreducible component of dimension \(j\), choose a linear projection to \(\mathbb P^{j+1}\) with a fiber through the point disjoint from the component and with image of dimension \(j\). To see the disjointness, project the component from the given point. In the target \(\mathbb P^{m-1}\), a generic plane of dimension \(m-j-2\) misses the image, because the expected intersection dimension is \(-1\). This gives a generic \((m-j-1)\)-plane through the point disjoint from the component. Choose the center as a hyperplane in the latter plane missing the given point. The center is disjoint from the component. Every projection fiber meets the center in a hyperplane, so it cannot contain a positive-dimensional projective subset of the component. Thus the projection is finite there, its image is a projective hypersurface of degree at most \(D\), and its equation separates the given point. Pull back and dehomogenize. Reducible unions use products of the separating forms and their summed degree. When \(j=m-1\), the component itself is a hypersurface and provides the separating equation.

Now take the space of all polynomials of degree at most \(D\) vanishing on the rational sample on the given component. Its basis can be chosen over \(\mathbb Q\) with integer coefficients of height \(n^{O(1/k)}\), by maximal minors of the sample evaluation matrix after clearing short-coordinate denominators. After base change its span contains every polynomial of degree at most \(D\) vanishing on the component, because the evaluation matrix has rational entries. Thus its common zero set is a subset of the original component containing all the sample, and has the claimed degree and coefficient bounds. If its dimension drops, descend immediately. Thus at the current dimension we can assume bounded-degree equations with coefficient height \(n^{O(1/k)}\), independent of \(\rho_k\).

Choose \(j\) bounded integer narrow linear forms \(\mathbf t\) dominant on every top-dimensional component, as in (FC1). Obtain \(m-j\) bounded-degree equations \(A_i\) for the trailing narrow coordinates, with coefficient height \(n^{O(1/k)}\), by elimination, omitting \(\mathbf t\)-only factors and taking squarefree parts. As in (FC1), the omitted factors and the trailing derivative zeros give proper lower-dimensional narrow sets, handled by induction. On the remaining points the rows \(\nabla A_i\) have rank \(m-j\). At dimension zero the narrow points are boundedly many, and each determines at most one wider rational \(V\) modulo \(n\), because \(H(V)^2<n\) for large \(k\).

### A first auxiliary polynomial of bounded degree in the wide variable

The exact system \(A_i=0,\ \phi=K\) in the \(m+1\) variables has rank \(m-j+1\) on the locus with the stated good derivatives: the \(A_i\) are independent of \(V\) and the congruence equation has unit derivative in \(V\). At an original tuple, change only \(V\) at each prime to make \(\phi=K\) exactly, keeping the original \(V\) modulo the full prime power. The reference set at each prime has dimension \(j\). Its \(T_r\) charts, partitioned at depth \(\lfloor e/k\rfloor\) for \(p^e\Vert n\), give \(O(n^{j/k+\varepsilon})\) records.

Choose a sufficiently large fixed \(M_0\) depending on \(m,C\), and put
\[
R=\lceil M_0\rho_k k+\sqrt{k}\rceil .
\]
Use the rectangular support with degree at most \(R\) in each \(t_i\) and degree at most \(M\le M_0\) in \(V\). Its size is \(D=(R+1)^j(M+1)\). Clear denominators in each row by the short projected-variable denominators to the power \(R\) and the \(V\)-denominator to the power \(M\). The integer row heights are \(n^{O_{m,C}(R/k)}\), since \(M\rho_k\le R/k\), and the row multipliers are units modulo \(n\). For a sufficiently large fixed \(A=A(m,C)\), any minor of size \(\lceil A R^j\rceil\) vanishes. Indeed, its \(j\)-parameter Taylor sum divided by its size is a large constant times \(R\). Choose \(A\) to make the divisibility from \(\lfloor e/k\rfloor\ge e/(2k)\) exceed the archimedean bound, and take \(r=O(R)\) above the last order. This is legitimate since \(R=o(k)\). Choose \(M\) large enough that this minor size is at most \(D/2\), and take \(M_0\) to majorize the resulting fixed values for all \(j\le m\). As in the two-short case, the rank defect gives an integer kernel vector with coefficients of height \(n^{O_{m,C}(R/k)}\) by pigeonholing the evaluations on independent rows.

The resulting nonzero \(G(\mathbf t,V)\) cannot vanish identically as a polynomial in \(V\) on any top narrow component, by dominance of \(\mathbf t\). Thus the locus where all its \(V\)-coefficients vanish projects inside a proper lower-dimensional narrow set. Descend there. If \(G\) is independent of \(V\), its vanishing directly descends the narrow dimension. Otherwise, pull \(G\) back to the original coordinates. It has degree \(O(R)\) in the short coordinates, bounded degree in \(V\), and coefficient height \(n^{O(R/k)}\).

### Finite wide fibers and a second, purely narrow auxiliary

For each remaining solution on \(G=0\), use the order \(r_0\) of its \(V\)-root, bounded by \(M\), so it lies on \(G_{r_0}=\partial_V^{r_0-1}G=0\) with \(\partial_VG_{r_0}\ne0\). Now \((A_i,G_{r_0})\) have full rank \(m-j+1\). Restrict to nonzero coordinates, where the monomial is regular, and take closures of the resulting open algebraic strata inside \(Z\times\mathbb A^1\). At each descendant we retain only the original sample tuples in that branch, so the bound on exact fibers continues to hold. Those strata have finite fibers over their narrow projections, hence dimension at most \(j\). Components with constant exact \(\phi\) give \(O(F_0)\) points. On any other geometric component a logarithmic Jacobian minor is generically nonzero, since the differential of \(\phi\) on a nonconstant component cannot vanish in characteristic zero. Choose one bounded-integer combination \(J\) of the minors generically nonzero on all those components. The bound on its integer coefficients may depend on the fixed algebraic degrees at this stage, hence on \(R\), which is harmless for \(n\to\infty\) with \(k\) fixed.

The \(J=0\) locus on these nonconstant components projects to a narrow variety of dimension at most \(j-1\). Its degree is bounded with the fixed parameters, by elimination, and its coefficient height may be arbitrary at this step. Descend on that projection and discard the old wide equation if desired. At the next dimension the degree split and the construction from the rational sample use only the degrees and the heights of the original narrow coordinates. The resulting equations provide the precision bounds.

For the remaining tuples \(J\ne0\). The cleared logarithmic row is of bounded degree in \(V\) and narrow variables. Thus \(J\) has narrow degree \(O(R)\), bounded wide degree and coefficient height \(n^{O(R/k)}\). Its numerator at a solution is at most \(n^{O(R/k)}\), because the \(A_i\) have low height and degree and all terms involving the wide variable have bounded degree in it. Discard primes with \(v_p(J)>e/3\), losing a log mass \(O(R/k)\log n=o(\log n)\). At good primes use a nonzero minor of \((\nabla A_i,\nabla G_{r_0},\nabla\phi)\) and Hensel lift to exact rank-full solutions with precision at least \(2e/3\). As usual the log row and the derivative of the cleared congruence equation agree up to a unit at the original approximate solution, modulo \(p^e\). The exact reference set has dimension \(j-1\).

On this reference set the \(T_r\) determinant is now purely in the \(j\) projected short coordinates. Use large fixed degree, independent of \(R\), to force the auxiliary in each record. The real sizes are \(n^{O(1/k)}\) per bounded single monomial. The retained prime powers have logarithmic weight \((1-o(1))\log n\), which suffices for the usual determinant forcing. The auxiliary polynomial in \(\mathbf t\) is nonzero, so cuts the dimension on every top narrow component. This step has \(O(n^{(j-1)/k+\varepsilon})\) records and descends the narrow dimension. For \(j=1\) the exact references are isolated points, and their records identify the original rationals, since the common reconstruction modulus exceeds \(n^{2\rho_k}\) and \(n^{2C/k}\) for large \(k\). Alternatively use the dimension-zero chart as a list of finitely many exact lifts at each good prime.

### Total cost and amplification

On each path of the recursion a narrow dimension \(j\) gives either the high-degree cover, at cost \(n^{1/k+\varepsilon}\), or the two chart stages, at total cost \(n^{(2j-1)/k+\varepsilon}\), before dimension decreases. The lower-dimensional exceptional branches satisfy the same upper bound. The constant-fiber components contribute \(O(F_0)\), with boundedly many components in each piece. Summing the dimension costs and allowing one covering exponent per dimension gives (GW). All degree, height exponent and Taylor precision requirements before the first narrow normalization depend only on \(m,C\) and the choice \(R\asymp M_0\rho_k k+\sqrt{k}=o(k)\), not on the fixed exponents \(a_i\). Higher degree bounds affect only auxiliary number constants, with a fresh narrow normalization at each descent.

For coprime integer narrow coordinates and signed exponents of absolute value \(\gg k\), the wider rational \(V\) may share prime supports with the narrow coordinates. The exact-fiber estimate generalizes the one in the two-short-coordinate note. If \(W=\phi\) is fixed, for any positive exponent \(a_i\) we have \(x_i^{a_i}\mid W_{\rm num}V_{\rm den}\), and for negative exponent \(x_i^{|a_i|}\mid W_{\rm den}V_{\rm num}\). Cap the contribution from \(W\) at the ceiling of its valuation divided by \(|a_i|\) for each prime. These choices are divisor-many, and the residual base factor to exponent \(|a_i|\) divides the corresponding part of \(V\). Thus the fiber bound is \(n^{\rho_k\sum_i1/|a_i|+\varepsilon}=n^{o(1/k)+\varepsilon}\). Amplify by coloring the narrow prime supports by coordinate, and independently coloring wider primes by numerator or denominator to prevent cancellation. Products keep the narrow coprimality, the wider height exponent is \(\ell\rho_k=o(1)\) at fixed \(\ell\), and the input multiplicity is divisor-many. Dividing the bound (GW) by arbitrarily large fixed amplification shows vanishing normalized cost for all narrow roots and the wider multiplier.

## Necessary conditions at normalized exponent \(1/3\)

Every limiting anatomic profile with normalized exponent \(d=1/3\) has
\[
 l_A=l_B=l_C=1/3,\qquad \Gamma\ge1/2 .
 \tag{GE}
\]
Here \(\Gamma\) is the sum of the nonatomic parts of the root mass distributions on the scaled degrees \(eL\), defined below. In particular the total residual mass \(\Delta\) of the [balanced box analysis](ABC_OMITTED_COORDINATE_GRAPH.md#necessary-conditions-at-the-balanced-13-endpoint) is at least \(1/2\). If the numerical short height tends to zero, the two tall distributions are nonatomic and each has mass \(1/3\), giving \(\Gamma\ge2/3\). The same conclusion holds for the other two letters if all heights are positive and one letter has its numerical weight on finitely many limiting positive scales.

First consider the letter masses. At \(d=1/3\) a supported letter is one of
mass strictly greater than \(1/3\). The [fixed-degree cutoff](ABC_EXCEPTIONAL_SET_AND_POWERFUL_TRIPLES.md#making-the-low-coefficient-short-in-anatomic-boxes) still applies, as its stopping threshold on that letter is positive. As in the [proof of the limiting \(1/3\) bound](ABC_OMITTED_COORDINATE_GRAPH.md#a-limiting-small-budget-coefficient-of-13), if two letters are supported the third has mass strictly less than \(1/3\), so a positive-weight modulus on it or a short-term real strip gives a contradiction. If one letter is supported and both other masses are strictly below \(1/3\), the two exchanges in that proof still have strict margins. Three supported letters have total mass greater than one. In the remaining unbalanced case there is one supported letter \(T\), another letter \(M\) of mass \(1/3\), and a letter \(S\) of mass \(s<1/3\).

If \(S\) has positive numerical height, its high-degree part supplies a
positive-weight full modulus at cost at most \(s\). To count the other two
terms under that modulus, first consider the structure on \(M\). If \(M\)
has positive numerical height, any positive-weight high block on \(M\) of
fullness order \(1/L\) costs at least \(1/3-s\). Otherwise fix \(S\) in full,
lock that block on \(M\), and apply (FC2) to the complete supported letter
\(T\). Choose a fixed \(0<\xi<(1/3-s)/4\), set aside the boundedly many
coordinates on \(M\) of normalized mass at least \(\xi\), and pass to a
subsequence on which their weights have limits. If their complement has
positive numerical weight, keep the degrees of order at least \(1/L\) that
give positive total weight there. The retained coordinates have mass at
most \(\xi+o(1)\). Greedy grouping into blocks of mass at most
\(2\xi+o(1)\) yields boundedly many blocks, one of which has positive
limiting weight. This contradicts the modulus cost \(1/3-s\). The finitely
many selected coordinates with positive limiting weight therefore account
for the height on \(M\). The rest of \(M\) has numerical height \(X^{o(1)}\),
which also holds for the whole letter when \(M\) has height zero.

Now use the modulus on \(S\). Retain the positive-weight coordinates on
\(T,M\) and combine their low cofactors into the wider ratio. There are only
boundedly many small roots and the wider ratio has height \(X^{o(1)}\).
By the amplified bound (GW), the count of the complete ratio costs
\(o(L)\) beyond the modulus. It recovers the triple by primitivity. Hence
\(d\le s\), a contradiction.

If \(S\) has numerical height zero, both \(T\) and \(M\) are tall. Under
any positive-weight full modulus on \(M\), (GW) releases the supported
letter \(T\) and the entire numerical short letter \(S\), with the latter
included in the wider cofactor. Thus every such modulus costs at least
\(1/3\). Applying the same high-block argument on \(M\) forces a single
positive-weight root with full normalized mass \(1/3\), and the other
factors have total mass tending to zero. Freeze those factors. The ratio
of this one root and the entire supported letter \(T\) lies in the real
short strip defined by the numerical short letter \(S\). Their remaining
variables are all of height \(X^{O(L)}\). Equation (SR3) therefore counts
at vanishing normalized cost and gives a contradiction. This proves that
every profile with \(d=1/3\) is balanced.

Call a positive-height letter finite-scale if its numerical weight outside
shrinking neighborhoods of some finite set of positive values of the
normalized degree \(eL\) tends to zero. Such a letter admits finitely many
small roots and a factor of numerical height \(X^{o(1)}\). For each of the
limiting degrees \(r>0\), group the corresponding original root coordinates
as a product, with exponent \(\lfloor r/L\rfloor\). The discrepancy in their
powers has height \(X^{o(1)}\) because their total root mass is \(O(L)\).
The factors outside the selected neighborhoods also have vanishing
numerical weight. The grouped roots have disjoint prime support, and any
prime overlap with the correction factor is allowed by the amplified
bound (GW).

When all three letters have positive numerical height, at most one is
finite-scale. Indeed, if two are finite-scale, under any positive-weight
full modulus on the third their signed ratio can be counted at vanishing
normalized cost. The modulus must cost at least \(1/3\). The same high-block
argument forces a single positive-weight root of full normalized mass
\(1/3\) on the third letter. This letter is then finite-scale as well.
Repeating the operation with the other letters as modulus forces a
full-mass positive root on each. Their normalized masses sum to one, in
contradiction with release (30).

When the short letter has numerical height zero, neither tall letter is
finite-scale. If one were, lock any positive-weight full modulus on the
other tall letter. Release the finite-scale tall term and the whole
numerical short term using (GW). As before the modulus letter must have a
single mass-\(1/3\) positive-weight root, so it is finite-scale too.
Reversing their roles forces both tall letters to have such roots. Freeze
their mass-\(o(1)\) complements. Then (SR3) on the two remaining roots
counts their signed ratio at vanishing normalized cost and recovers the
triple, a contradiction.

For a positive-height letter that is not finite-scale, there are
positive-weight full moduli of arbitrarily small normalized cost. To see
this, suppose every available such modulus costs at least some fixed
\(t>0\). Choose a fixed cutoff \(\xi<t/4\) and set aside the boundedly many
coordinates of normalized mass at least \(\xi\). A complement with positive
limiting weight provides a positive-weight high block of mass at most
\(2\xi+o(1)\) by the same grouping argument used above. This violates the
cost \(t\), so the finitely many coordinates with positive limiting weight
carry the numerical height. Each of these coordinates has a limiting
scaled degree in \((0,\infty)\), making the letter finite-scale.

### Nonatomic mass in the scaled degrees

For each \(T\) take a weak limit \(\nu_T\) of the measures
\[
 \sum_e \beta_{T,e}\,\delta_{eL},\qquad \beta_{T,e}=L_{T,e}/L ,
\]
along a further subsequence. Here \(\delta_r\) denotes a unit point mass
at \(r\). These measures have bounded first moment since
\(\sum_e eL\beta_{T,e}\) is the numerical height exponent of the letter up
to the vanishing truncation loss. Their mass at scales at least \(R\) is
therefore \(O(1/R)\), and \(\nu_T\) has total mass \(l_T\) on
\([0,\infty)\). Let \(c_T\) be its nonatomic mass and put
\(\Gamma=\sum_T c_T\). For each of the fixed deletions defining
\(\delta_T\), the removed mass in a weak limit is a finite sum of point
masses. A fixed-degree coordinate has scale zero, and a deleted coordinate
with scale tending to infinity has mass tending to zero. Thus
\(c_T\le\delta_T\).

Suppose there is a positive-weight full modulus on \(T\) of cost \(\mu\),
with fullness comparable to \(1/L\). For the other two letters \(U,V\),
choose finite sets of positive atoms of their degree measures. Around
each selected scale \(r\), take a window shrinking sufficiently slowly
that the total mass in it tends to the atom \(\nu(\{r\})\). Take also a
slowly shrinking window at zero, whose mass tends to \(\nu(\{0\})\).
One can choose the windows together by first taking fixed disjoint
intervals with endpoints of zero measure, using weak convergence, and
then decreasing the sizes along a diagonal sequence. Freeze all
coordinates outside those windows on \(U,V\). The factors in the zero
window have total numerical weight tending to zero. In each window
around \(r>0\), group the roots into a product with exponent
\(\lfloor r/L\rfloor\), leaving a rational correction of height
\(X^{o(1)}\). This height bound follows by multiplying the maximum
difference \(o(1/L)\) in exponents by the total root mass \(O(L)\).
There are boundedly many grouped roots after the finite sets are fixed,
and their supports are disjoint by primitivity and the anatomic
decomposition. Combine the corrections and the zero-window factors into
the wider ratio. The amplified bound (GW) counts the complete ratio at
vanishing normalized cost beyond the frozen coordinates and the
modulus. If no grouped roots remain, the wider rational ratio is
recovered by rational reconstruction.

Choose successively larger finite sets of positive atoms after taking the
limiting count. The unused atomic mass tends to zero, so the resulting
bound is
\[
 d\le\mu+c_U+c_V .
 \tag{GE1}
\]
This argument also works before imposing balance.

A finite-scale letter has \(c_T=0\), since any mass in a compact interval
of positive scales away from its finitely many concentrations would
carry positive numerical weight. A letter of numerical height zero
likewise has \(c_T=0\), since its weak limit is supported at zero. If all
three letters have positive height and none is finite-scale, each admits
arbitrarily cheap positive-weight moduli. Applying (GE1) with each
choice of modulus letter gives all three pair bounds \(c_U+c_V\ge1/3\),
so that \(\Gamma\ge1/2\). If all heights are positive and one letter is
finite-scale, the other two admit arbitrarily cheap moduli. Equation
(GE1) gives \(c_T=1/3\) on each of those other letters. When the short
height tends to zero, both tall letters admit arbitrarily cheap moduli
and the short letter has \(c_T=0\). Thus (GE1) again gives \(c_T=1/3\)
on each tall letter. More generally, if any \(c_T=0\), the other
two values of \(c\) equal \(1/3\). This follows from the pair bounds when
all three positive-height letters are nonfinite-scale, and from the two
case analyses just given otherwise.

The three uniform distributions of mass \(1/3\) on the interval \([2,4]\)
satisfy these degree-measure constraints and have first moments one. They
illustrate the range of necessary limiting conditions that remain after
the count.

The [quantitative window estimates](ABC_QUANTITATIVE_DEGREE_WINDOWS.md) give additional finite-width exclusions and local mass bounds at equality in \(\Gamma\ge1/2\).

The [short-height sector estimate](ABC_SHORT_HEIGHT_AND_REAL_GAP.md#a-strict-saving-when-the-short-term-is-small) excludes numerical short height zero at the limiting coefficient \(1/3\).

At the same endpoint, the [variable-mesh estimates](ABC_VARIABLE_MESH_DEGREE_GROUPING.md)
bound the short height gap by the tall pair's nonatomic excess above
\(1/3\). They also bound the numerical weight in high scaled degrees,
including escaped first moment on each letter, by the excess on the
complementary pair.

## A quarter inequality for a letter of small radical mass

If a letter \(C\) has normalized root mass \(l_C<d\), both other
letters have mass at most \(1-2d\). Consequently every limiting
small-budget box satisfies
\[
 4d\le 1+\min\{l_A,l_B,l_C\}.
 \tag{LQ1}
\]
Here the fixed-degree cutoff is used in the usual regime \(L\)
comparable to the chosen radical budget. A supported letter is
one with mass strictly greater than \(1-2d\); it has boundedly many
roots of height \(X^{O(L)}\) and a cofactor of numerical height
\(X^{o(1)}\).
[Fixed-degree cutoff](ABC_EXCEPTIONAL_SET_AND_POWERFUL_TRIPLES.md#making-the-low-coefficient-short-in-anatomic-boxes)

To prove the first claim, suppose \(A\) is supported. Freeze the
complete \(C\), at normalized cost at most \(l_C\) by fixing its roots
and the initial exceptional coefficient. Any positive-weight block
on \(B\) formed by full powers of degree at least \(\kappa/L\),
for a fixed \(\kappa>0\), can serve as a modulus. If the block has
normalized cost \(\mu\), the congruence for \(A\) modulo this block
has a unit target given by the fixed \(C\). The amplified (GW)
release of \(A\) gives \(d\le l_C+\mu\). Its output exact value
recovers the triple since \(C\) is fixed. Thus every such \(B\)-block
costs at least \(d-l_C>0\).

If \(B\) has positive height \(h_B\), set aside coordinates of
normalized mass at least a fixed \(\xi<(d-l_C)/4\). If their
complement has a positive numerical contribution \(w\), retain
the part in degrees at least \(w/(2L)\), discarding at most
\(w/2+o(1)\) numerical weight. Its individual normalized root
masses are at most \(\xi+o(1)\); greedy grouping into blocks of
mass at most \(2\xi+o(1)\) gives boundedly many blocks. For
example every completed pair of consecutive blocks has
combined mass at least \(\xi+o(1)\), except for boundary
constants. Some block has positive weight, a contradiction.
Thus boundedly many heavy coordinates of \(B\) carry its height,
and the remainder is a cofactor of numerical height tending
to zero. If instead \(B\) has numerical height zero, its whole
integer is such a cofactor.

If \(C\) has positive height, its high part supplies a full
positive-weight modulus at cost at most \(l_C\). Release
\(A,B\) jointly by (GW) under that modulus. Their exact signed
ratio determines the primitive pair and hence the triple,
giving \(d\le l_C\). If \(C\) has numerical height zero,
\(A,B\) are the tall pair. Their complete ratio, with the
appropriate orientation, is within \(X^{-1+o(1)}\) of one.
Each term consists of boundedly many small roots and a
cofactor of height \(X^{o(1)}\). The unequal-root real strip
count (RUW1) gives normalized cost arbitrarily small, by
taking its chart parameter \(s\) arbitrarily small first
and then letting the cofactor exponent tend to zero.
The exact ratio again determines the pair.
Both conclusions contradict \(d>l_C\).
[Unequal-root real strip](ABC_SHORT_HEIGHT_AND_REAL_GAP.md#a-real-strip-with-unequal-root-heights)

If \(d>(1+l_C)/4\), the known \(d\le1/3\) forces \(l_C<1/3\)
and then \(l_C<d\). The relation
\(l_A+l_B=1-l_C>2(1-2d)\) would force one of the other letters
to be supported. Apply this to a letter of minimum mass to obtain (LQ1).

### Small radical on any one letter

For fixed \(\vartheta\ge0\) and \(\eta>0\), count ordered
positive primitive \(abc\) triples of budget \(\lambda\)
on which at least one letter \(T\) satisfies
\(\operatorname{rad}(T)\le c^{\vartheta\lambda}\). For all
sufficiently small fixed \(\lambda\), the count is
\[
 O_{\lambda,\vartheta,\eta,\varepsilon}\!
 \left(X^{(\min(1/3,(1+\vartheta)/4)+\eta)\lambda+\varepsilon}\right).
 \tag{LQ7}
\]
Indeed, in a sequence violating the bound with
\(\lambda\to0\), choose the same letter \(T\) on a
subsequence. The anatomic boxes have
\(L\le(1+o(1))\lambda\), \(L_T\le(\vartheta+o(1))\lambda\),
and count exponent of order at least \(\lambda\).
The fixed-two-letters bound gives \(L\) of order at least
\(\lambda\). By (LQ1) their exponent is at most
\((L+L_T)/4+o(\lambda)\). The global one-third estimate
gives the other term in the minimum.

The same argument gives \(\max_i l_i\le1-2d\).
If a letter has mass below \(d\), the other two satisfy
this upper bound by the supported-letter exclusion, and
the first does because \(d\le1/3\). If all letter masses
are at least \(d\), the upper bound follows from their
sum being one. Thus for \(d\ge1/3-\xi\), every mass is in
\([1/3-4\xi,1/3+2\xi]\), and the total absolute
deviation from \(1/3\) is at most \(8\xi\).

## A global inequality for a low-mass letter and diffuse degrees

For every limiting small-budget box with \(L\) comparable to its
radical budget, and each permutation \(A,B,C\) of its letters,
\[
 d\le \min_C\max\{l_C,\ c_C+\min_{T\ne C}c_T\}.
 \tag{LQ8}
\]
In particular, if \(l_C=s<1/3\) and \(d=(1+s)/4\), then
\[
 l_A=l_B=(1-s)/2,\qquad
 \min(c_A,c_B)\ge d-c_C\ge(1-3s)/4,\qquad
 \Gamma\ge(1-s)/2.
 \tag{LQ9}
\]
Each of \(A,B\) has positive numerical height and neither is finite-scale.

To prove (LQ8) it suffices to consider \(s=l_C<d\).
Call a letter releasable if it consists of boundedly many roots of
height \(X^{O(L)}\) and a cofactor of numerical height \(X^{o(1)}\).
A finite-scale letter is releasable by grouping its shrinking degree
windows. A letter of numerical height zero is releasable with only
the cofactor. Suppose \(A\) is releasable. Freeze the complete \(C\)
at normalized cost \(s\). If \(B\) has a positive-weight full block of
cost \(\mu<d-s\) and fullness comparable to \(1/L\), (GW) releases \(A\)
at cost \(s+\mu\). This would contradict \(d>s+\mu\). Hence the
high-block stopping argument in the proof of (LQ1) makes \(B\)
releasable too.

If \(C\) has positive height, use a positive-weight full modulus on
\(C\) at cost at most \(s\), and release the complete ratio of \(A,B\)
by (GW). If \(C\) has height zero, \(A,B\) are the two tall letters,
whose relative difference is \(X^{-1+o(1)}\). In this case the strip
count (RUW1) releases their ratio at arbitrarily small cost. The
complete ratio determines the two coprime values and hence the triple.
Both cases contradict \(d>s\). The same argument applies with \(A,B\)
interchanged.

Since neither \(A\) nor \(B\) is releasable, the high-block stopping
argument gives positive-weight full blocks of arbitrarily small fixed
normalized cost on each. Apply (GE1) with the modulus on \(B\) to get
\(d\le c_C+c_A\), and then with the modulus on \(A\) to get
\(d\le c_C+c_B\). In each application freeze the complement of the
selected atom windows on the free letters, including any escaped
first moment. This proves (LQ8). At equality in (LQ1) with \(s<1/3\),
the bounds \(l_A,l_B\le1-2d\) and \(l_A+l_B=1-s\) force the stated
masses. Thus \(\Gamma\ge 2d-c_C\ge(1-s)/2\), proving (LQ9).
For \(s=0\), each opposite letter has root mass \(1/2\) and nonatomic
mass at least \(1/4\).

### A saving when an opposite letter concentrates in finitely many bands

Fix \(0\le\vartheta<1/3\), a positive \(b\), and a positive integer \(K\).
There are \(\delta,\eta,\lambda_*>0\) such that the following triples
of height at most \(X\) have count
\(O_{\vartheta,b,K,\lambda,\varepsilon}
(X^{((1+\vartheta)/4-\delta)\lambda+\varepsilon})\)
for each fixed \(0<\lambda\le\lambda_*\).
Require \(\operatorname{rad}(abc)<c^\lambda\) and
\(\operatorname{rad}(C)\le c^{\vartheta\lambda}\). On a different
letter, require \(\lambda v_p\) for primes of total logarithmic root
mass at least \(((1+\vartheta)/4+b)\lambda\log c\) to lie in
the union of at most \(K\) intervals of length \(\eta\).

Indeed, if the saving fails along vanishing \(\delta,\eta,\lambda\),
(LQ1) and \(L_C\le(\vartheta+o(1))\lambda\) force \(L/\lambda\to1\),
\(l_C\to\vartheta\), and \(d\to(1+\vartheta)/4\) on the selected
boxes. The captured root mass survives the anatomic truncation.
More explicitly, choose compression with \(k\lambda\to\infty\).
The logarithmic radical contribution at degrees at least \(k\) is
at most \((\log c)/k=o(\lambda\log c)\), and the subpower losses
at fixed budget can be taken \(o(\lambda)\) before the small-budget
limit. The other indices retain their exact valuations. If the intervals
vary with the triple, fix the captured subset of exact indices on
each truncated box at finite cost before taking the height limit.
The mass bounds give tightness, so the captured submeasure has a
further weak limit supported on at most \(K\) points. Interval centers
that escape every compact interval contribute zero mass in that limit.
The opposite letter must therefore have atomic mass at least
\((1+\vartheta)/4+b\). But (LQ9) bounds its atomic mass by
\((1-\vartheta)/2-(1-3\vartheta)/4=(1+\vartheta)/4\).


### A saving with an opposite letter near its fullness threshold

Fix \(0\le\vartheta_0<1/3\). For \(0\le\vartheta\le\vartheta_0\),
consider the same radical budget and the cap
\(\operatorname{rad}(C)\le c^{\vartheta\lambda}\).
On another letter \(A\) of height \(h_A=\log A/\log c\),
require the logarithmic radical at scaled degrees
\(\lambda v_p(A)<2h_A/(1-\vartheta)-\tau\) to be at most
\(\tau\lambda\log c\). This includes \(A\) \(p\)-full with
\(p\lambda\ge2h_A/(1-\vartheta)-\tau\).
There are \(\delta,\tau,\lambda_*>0\) depending only on
\(\vartheta_0\) such that at fixed \(0<\lambda\le\lambda_*\) the
count is \(O(X^{((1+\vartheta)/4-\delta)\lambda+\varepsilon})\).

Otherwise let the three allowances tend to zero and pass to
\(\vartheta\to\vartheta^*\le\vartheta_0\) and \(h_A\to h_A^*\).
The same box bounds force equality in (LQ1) with \(l_C=\vartheta^*\)
and \(l_A=(1-\vartheta^*)/2\). If \(h_A^*>0\), the exact degree
positions and the low-degree root bound put the measure of \(A\)
on \([2h_A^*/(1-\vartheta^*),\infty)\). Its total mass and
first-moment upper bound force a single atom at the threshold
with the entire numerical weight. It is then releasable.
If \(h_A^*=0\), it is releasable by its vanishing numerical height.
Both contradict the proof of (LQ8), since \(\vartheta^*<d\).

## Adaptive finite-basis approximation of prime valuations

Let \(T\) be a positive integer in a triple of height \(X_0\),
and fix a small budget parameter \(\lambda>0\).
For fixed \(m\ge1\) and \(H>0\), consider the intrinsic condition
that there exist signed integers \(q_1,\ldots,q_m\) and
\(m_{j,\ell}\), depending on this value of \(T\), with
\[
 \sum_{\ell\mid T}|m_{j,\ell}|\log\ell\le H\lambda\log X_0
       \quad(1\le j\le m),\qquad
 \sum_{\ell\mid T}
 \left|v_\ell(T)-\sum_{j=1}^m q_jm_{j,\ell}\right|\log\ell
      \le\delta\log X_0.
 \tag{AB1}
\]
Here the sums are over prime divisors of \(T\). The exponent list
and the primes used in several groups may vary with \(T\) and
with the height. This condition can be used with the
bounded-root counting theorems without assuming their
constants polynomial in variable exponent parameters.

### Turning adaptive data into finitely many exponent lists

At fixed \(\lambda\), choose \(K=\lceil\lambda^{-2}\rceil\).
On every prime of valuation \(e>K^2\), write
\(e=K\lfloor e/K\rfloor+(e\bmod K)\). These primes supply
one new narrow root of height at most \(X_0^{1/K}\) and
exponent \(K\), with remainder numerical weight at most
\[
 \frac{K\sum_{v_\ell(T)>K^2}\log\ell}{\log X_0}\le1/K.
\]
They will not be included in the other groups.

Among the remaining primes of valuation \(e\le K^2\),
discard from the groups those with \(|m_{j,\ell}|>K^2\)
for some \(j\), or with absolute residual
\(\left|e-\sum_jq_jm_{j,\ell}\right|>K^2\).
The product of the discarded actual prime powers has
numerical weight at most \(mH\lambda+\delta\). In fact
for the first test, sum the \(m\) inequalities in (AB1)
and use \(e\le K^2\le |m_{j,\ell}|\) at a witnessing
index. For the second use \(e\le K^2\) and the residual
sum in (AB1). Double counting primes in both classes
only enlarges this upper bound.

On the retained primes, preserve the exact predicted
exponents \(z_\ell=\sum_j q_j m_{j,\ell}\). The integer
linear system for this list of predictions has coefficients
bounded by \(K^2\) and right sides bounded by \(2K^2\).
It has an integral solution by construction. There is
another integral solution with all entries
\(K^{O_m(1)}\), independent of \(X_0\).
Indeed choose independent rows and an invertible square
minor in their columns of size their common rank.
Let \(\Delta\ne0\) be that determinant. Replace the free
coordinates of an existing integral solution by their
least nonnegative residues modulo \(|\Delta|\).
Solve the pivot coordinates by the adjugate formula.
Each free coordinate change is divisible by \(\Delta\),
so the pivot changes are integers. The independent
rows, and hence all rows of the consistent system, are
still satisfied. Hadamard's bound or the determinant
expansion and the adjugate formula give
\(K^{O_m(1)}\) bounds for both free and pivot coordinates.
In rank zero choose all coordinates zero.

Thus, at fixed \(\lambda\), all resulting exponent
lists are in a fixed finite set, which may be
pigeonholed at constant cost before taking the height
limit. Define positive rational group roots
\(z_j=\prod_\ell\ell^{m_{j,\ell}}\) on the retained
primes, keeping their height bound \(X_0^{H\lambda}\).
One obtains
\[
 T=T_{\rm pack}^{K} U\prod_{j=1}^m z_j^{q'_j},
 \qquad
 H(T_{\rm pack})\le X_0^{1/K},\qquad
 H(U)\le X_0^{\,2\delta+mH\lambda+1/K}.
 \tag{AB2}
\]
For \(U\) account for the retained residual from (AB1),
the discarded actual prime powers, and the packing
remainder. Thus as \(\lambda,\delta\to0\), this is a
bounded number of short roots and a wide cofactor of
vanishing numerical height. The construction still
works after restricting all prime data to any subset
of primes, such as those not dividing a fixed gcd.

### Signed amplification with overlapping groups

Whenever the exact positive rational output determines
the primitive pair or other objects being counted,
keep one representing tuple per exact output.
For modular amplification, colour each prime by its
numerator/denominator sign in that exact output. Every
numerator and denominator has polynomial height in the
height parameter, hence its number of distinct primes
is \(o(\log X)\). Some colouring retains a subpower
fraction of the sample. On that subfamily products of
a fixed number \(\ell\) of outputs have no cancellation.
Each product therefore has only divisor-many
decompositions into the \(\ell\) actual ordered outputs,
by factoring its numerator and denominator. Choose one
coordinatewise product tuple per exact product output.
Short-root and wider height bounds are multiplied in
their exponents by at most \(\ell\), and the residue
target is raised to the \(\ell\)-th power.
Use the unamplified (GW) theorem with exact-fibre bound
one on these tuples. This proves the usual adjustable
amplification cost for the original exact outputs,
even when signed root groups overlap in prime support.
The group supports can overlap. For (RUW1) simply use one tuple per exact
output without amplification.

### Sector with a small radical on a different letter

Fix a positive integer \(m\) and real \(H,\omega>0\). There are
\(\delta_*,\lambda_*>0\), depending only on these constants,
with the following property. For fixed
\(0<\lambda\le\lambda_*\), count primitive positive triples
of height \(c\le X\), radical budget \(\operatorname{rad}(abc)<
c^\lambda\), with \(\operatorname{rad}(C)\le c^{\vartheta\lambda}\)
on one indicated letter and (AB1) with \(X_0=c\) and
\(\delta\le\delta_*\) on a different indicated letter \(A\).
Uniformly for \(0\le\vartheta\le1\), their count is
\[
 O_{m,H,\omega,\lambda,\varepsilon}
      \left(X^{(\vartheta+\omega)\lambda+\varepsilon}\right).
 \tag{AB3}
\]
In particular for \(0\le\vartheta\le\vartheta_0<1/3\) this
gives a strict saving from the quarter coefficient, by
choosing \(\omega<(1-3\vartheta_0)/4\).
It is stronger than merely excluding quarter equality.

Otherwise choose failing small budgets and errors tending
to zero, always taking the height limit first at each fixed
parameter choice. Fix the exponent list in (AB2) at that
height limit. Select ordinary anatomic boxes at losses
\(o(\lambda)\), with total root mass \(L\asymp\lambda\)
by the assumed positive count exponent and the root budget.
Then the selected profiles satisfy \(d-l_C\ge\omega/2\),
since \(D_{\rm box}\ge(\vartheta+\omega-o(1))\lambda\) and
\(L_C\le(\vartheta+o(1))\lambda\). Formula (AB2) makes
\(A\) releasable, with arbitrary overlap allowed by the
preceding amplification. In the proof of (LQ8), fix
the complete \(C\) and use every positive-weight high block
of the opposite \(B\) as modulus. Its cost has to be at
least \(d-l_C>0\), forcing \(B\) to be releasable by the
high-block stopping proof. Under a positive-height full
modulus of \(C\), or in the real short strip if \(C\)
has zero height, release the complete ratio of \(A,B\).
The exact value determines the primitive pair. The
resulting cost is at most \(l_C\), a contradiction.

One special case is the intrinsic \(m\)-center deviation
\[
 \inf_{r_1,\ldots,r_m\ge0}
 \sum_{\ell\mid A}\frac{\log\ell}{\lambda\log c}
       \min_j |\lambda v_\ell(A)-r_j|.
\]
If this tends to zero, assign each prime to a nearest
center and round \(r_j/\lambda\) to an integer exponent.
The group height bound holds with \(H=1\), and rounding
adds at most \(\lambda\) to the original numerical
error. The centers are allowed to diverge. The same
conclusion follows with the squared deviation, by
Cauchy–Schwarz since the total normalized mass is at most
one. In particular (AB3) covers a single exact \(p\)-th
power on \(A\) whenever \(p\lambda\) is bounded below.

## A strict one-third saving with any one complete letter fixed

There exist absolute constants \(\sigma>0\) and \(\lambda_0>0\)
with the following property. For \(0<\lambda<\lambda_0\),
\(\varepsilon>0\), \(T\in\{a,b,c\}\), and \(K\ge1\), the number of
ordered primitive triples
\[
 a+b=c\le X,\qquad \operatorname{rad}(abc)<c^\lambda,\qquad T=K
\]
is
\[
 \ll_{\lambda,\varepsilon}X^{(1/3-\sigma)\lambda+\varepsilon},
 \tag{FT1}
\]
uniformly in \(K\) and in the position \(T\). In particular, if a specified
letter is drawn from a dictionary of at most \(X^{\theta\lambda}\)
complete values for a fixed \(0\le\theta<\sigma\), the count is
\[
 \ll_{\lambda,\varepsilon}
 X^{(1/3-\sigma+\theta)\lambda+\varepsilon}.
 \tag{FT2}
\]

A multiplicative decomposition of one free term gives a large common
neighborhood in a bipartite graph. Prime-power factors of the other
free term on different rows then impose several congruences on the
same column value.

### Passing to a limiting fixed-target box

Write \(d_0=1/3\). If (FT1) fails for every choice of
\(\lambda_0,\sigma>0\), choose counterexample budgets
\(\lambda_j\to0\) and prospective savings \(\sigma_j\to0\).
For each fixed \(j\) the failure occurs at arbitrarily large \(X\), with
a complete value \(K\) that may depend on \(X\). Thus the passage to
limiting boxes must be uniform in \(K\).

Reduce to dyadic shells for the largest letter and use the anatomic
factorization and dyadic root boxes with exceptional coefficients fixed,
at a total loss \(o(\lambda_j)\) in the height exponent. At each fixed
\(j\) the factorization and all degrees are fixed first, and the shell
height tends to infinity before \(j\) tends to infinity. Finite dyadic
box counts, divisor losses, and arbitrary small height-exponent losses
are absorbed at this first limit. This gives a limiting sequence of boxes
whose count is \(X^{D+o(L)}\) in normalized counting exponents, with
\(L\le \lambda_j+o(\lambda_j)\) the total root mass and
\(D\ge(d_0-o(1))\lambda_j\). The elementary count by two complete
letters gives \(D\le 2L/3+o(L)\), so \(L\) is comparable to \(\lambda_j\).
The existing unrestricted \(1/3\) argument rules out a limiting
\(D/L>d_0\). We can therefore assume \(D/L\to d_0\). The necessary
conditions at equality give letter masses \(l_A=l_B=l_T=d_0\), where
\(A,B\) name the two free letters. They also give the degree measures
\(\nu_A,\nu_B,\nu_T\) described in (GE1). Write \(n=X^L\), so the
original count in the box is at least \(n^{d_0-o(1)}\).
A limiting exponent statement means the \(X\)-exponent
statement at each fixed box, followed by the small-budget limit.

The complete value of \(T\) is fixed throughout each box, with
\(1\le T\le X\). Fix the sign convention
\(\epsilon_A A+\epsilon_B B+\epsilon_T T=0\), with each sign in
\(\{1,-1\}\). Consequently either free complete value determines the
other and hence the whole triple, uniformly in the position of \(T\).

### The free letters have wholly nonatomic degree measures

The endpoint conditions preceding (GE1) imply that two finite-scale
letters are impossible when all heights are positive. They also imply
that a tall letter is never finite-scale if the short height vanishes.
Here a finite-scale letter is one whose numerical weight is
concentrated on finitely many limiting positive scaled degrees.
A positive-height letter which is not finite-scale admits
positive-weight full moduli of arbitrarily small normalized root cost,
as in the high-block argument preceding (GE1). A modulus in that
argument has at least \(k\)-fullness with \(k\asymp 1/L\).

Thus if one free letter has zero limiting height, the other is tall
and supplies such cheap moduli. Upon fixing the modulus value on the
tall letter, the zero-height free letter lies in one residue class
because \(T\) is fixed. Its height is eventually below the modulus,
so the intended count is at most the arbitrarily small root cost
of the modulus. If both free letters have positive height and one
is finite-scale, the other has cheap moduli. Group the finite-scale
letter around its finitely many positive scaled degrees and apply
the amplified (GW) under such a modulus. The correction has height
\(X^{o(1)}\) and the target is a unit after fixing any outside
coefficients, by primitivity. This releases the finite-scale free
term at vanishing normalized cost. Again the count is at most the
arbitrarily small modulus cost. These arguments exclude both cases.

We now have cheap moduli on both free letters. With a modulus of cost
\(\mu\) on \(A\), apply the proof of (GE1) to \(B\) alone. The other
complete letter \(T\) is already fixed. Freeze the part of \(B\)
outside the slowly shrinking windows about zero and finitely many
positive atoms of \(\nu_B\). Group the roots in those windows and apply
the amplified (GW) to \(B\)'s monomial and its larger correction. If
there are no positive windows, use rational reconstruction. Its exact
value determines the triple. Increasing the finite set of atoms after
the limiting count yields \(d_0\le \mu+c_B\), with \(c_B\) the
nonatomic mass of \(\nu_B\). The same reasoning on \(A\) gives
\(d_0\le\mu+c_A\). Letting \(\mu\) tend to zero shows
\[
 c_A=c_B=d_0.
 \tag{FT3}
\]

### Common neighborhoods of multiplicative decompositions

Relabel the free letters so that \(B\) is tall and \(h_B=1\).
This is possible because two terms have limiting height one.
Write \(h_A\le1\) for the other free height. Let \(e_B\ge0\) be the
first moment lost at infinity in \(\nu_B\), so
\[
 \int_0^\infty r\,d\nu_B(r)=1-e_B .
\]
By (FT3), \(\nu_B\) is not supported at one point, and there is a finite
band \([s,S]\) of mass \(\tau>0\), with zero-mass endpoints and
\[
 s>(1-e_B)/d_0.
\]
Choose \(q_0\in(0,d_0)\) with \(e_B+s q_0>h_A\). Fix an integer \(m\)
large enough that \(\beta=q_0/m<\tau\). Choose positive
\(\xi,\zeta,\delta\) so small that
\[
 q_0+\zeta<d_0,\qquad
 e_B-\delta+s\bigl(q_0-m(\xi+\zeta)\bigr)>h_A .
 \tag{FT4}
\]
Finally choose \(p>0\) with \(p<\xi\), \(p<d_0-q_0-\zeta\),
and less than the mass of some compact band of \(\nu_A\).
Set \(q=d_0-p\).

Each coordinate in a fixed compact band has mass \(o(1)\) if its degree
measure is nonatomic. Otherwise a sequence of coordinates with mass
bounded below would produce an atom at a limit of their scaled
degrees. Hence a fixed subset of the coordinate positions of \(A\)
can be chosen with normalized root mass \(p+o(1)\), using only that
compact band. For each solution let \(U\) be the product of these
coordinates to their respective powers, including the fixed exceptional
coefficient if desired, and write \(A=UV\). The possible positive
integer values of \(U,V\) number at most \(n^{p+o(1)},n^{q+o(1)}\),
respectively, by the box bounds. Select one factorization per solution
if needed. The graph on the values \(U,V\) has a distinct edge for each
triple, since \(UV\) and the fixed \(T\) determine the triple. It has
at least \(n^{p+q-o(1)}\) edges.

If \(d(V)\) is the number of \(U\)-neighbors of \(V\), the total number
of ordered \(m\)-stars centered on the \(V\)-side, allowing repetitions,
is
\[
 \sum_V d(V)^m\ge
 \frac{(\sum_V d(V))^m}{(\#V)^{m-1}}
 \ge n^{mp+q-o(1)}.
\]
Stars with a repeated \(U\) contribute at most
\(\binom m2(\#U)^{m-1}\#V\le n^{(m-1)p+q+o(1)}\),
which is smaller by a fixed power since \(p>0\). Dividing the count
of distinct stars by \((\#U)^m\) gives distinct fixed integers
\(U_1,\ldots,U_m\) with at least \(n^{q-o(1)}\) common neighbors \(V\).
They may depend on the height. In row \(i\) the complete other free
letter is
\[
 B_i(V)=\eta T+\omega U_iV
\]
for fixed signs \(\eta,\omega\in\{1,-1\}\).
Both \(U_i\mid A_i\) and \(\gcd(A_i,B_i)=1\) hold for every edge.

### Stripping shared primes and retaining escaped weight

Use a fixed subset of \(B\)'s coordinate positions in the band
\([s,S]\) of normalized mass \(\beta+o(1)\), which exists by the same
small-coordinate observation. In row \(i\), write \(R_i\) for
the product of their unpowered roots. These roots are squarefree
and have disjoint supports. To ensure squarefreeness one may use the
specific [BBLT Lemma 2.2](https://arxiv.org/html/2410.12234v2#S2)
factorization with truncation error \(o(L)\): the only possibly
nonsquarefree root has degree \(k=\lceil1/\varepsilon_0\rceil\),
so \(kL\to\infty\) and lies outside every fixed band.

A prime dividing \(B_i(V)\) and \(B_j(V)\) divides
\(T(U_j-U_i)\), by taking \(U_jB_i-U_iB_j\). Put
\[
 J=\prod_{i< j} T(U_j-U_i).
\]
This is nonzero and \(|J|\le X^{O_m(1)}\) uniformly. Call a column bad
if some \(R_i\) shares with \(J\) a squarefree divisor of size at least
\(n^\xi\). For that row choose the divisor of \(\operatorname{rad}(J)\)
and assign its primes to the selected positions whose roots they
divide. There are \(X^{o(1)}\) choices at each fixed truncation by the
uniform divisor and prime-allocation bounds. If the assigned divisors
are \(g_e\) and the root box scales are \(H_e\), the number of multiples
in the relevant interval is \(O(H_e/g_e)\) whenever the interval
contains a positive multiple. Multiplying over all root positions on
\(B\) therefore bounds the number of possible row values by
\(n^{d_0-\xi+o(1)}\), with divisor losses absorbed into arbitrarily
small height exponents at fixed \(j\). Since \(B_i(V)\) determines
\(V\), the same bound holds for bad columns. It is smaller than
\(n^{q-o(1)}\) because \(\xi>p\).

Choose a fixed cutoff \(H>S\), a continuity point of \(\nu_B\), so
large that the total root mass of scaled degrees above \(H\) is less
than \(\zeta\). The numerical weight of that tail has limit
\[
 e_B+\int_{(H,\infty)}r\,d\nu_B(r)\ge e_B.
\]
The mass bound follows from the uniform first moment and the weight
identity follows by bounded convergence below the continuity cutoff
and subtraction from the total height. Denote the product of the tail
powers on row one by \(F\); its root record has at most
\(n^{\zeta+o(1)}\) choices, its radical is at most \(n^{\zeta+o(1)}\),
and its numerical weight is at least \(e_B-\delta\) at sufficiently
small budget and large height. Use \(F=1\) when \(e_B=0\) if desired.

For a good column record the tail root tuple giving \(F\) and the
selected band root tuples on every row. Strip from every recorded
band root all primes dividing \(JF\). The cost in root mass on each
row is at most \(\xi+\zeta+o(1)\) by the good-column condition and
the tail radical bound. The product of the remaining band powers on
row \(i\), denoted \(M_i\), has lower numerical weight
\(s(\beta-\xi-\zeta)+o(1)\). All the \(M_i\) are pairwise coprime
because any common prime would divide \(J\). Each is also coprime to
\(F\). Every prime of \(M_i\) is prime to \(U_i\) by primitivity, and
the same holds for \(F\) and \(U_1\). The recorded factors therefore
impose one residue class for \(V\) modulo \(F\prod_i M_i\). Its
numerical weight is at least
\[
 e_B-\delta+s\bigl(q_0-m(\xi+\zeta)\bigr)+o(1)>h_A.
\]
As \(V\le A\le X^{h_A+o(1)}\), there is at most one \(V\) per
record. There are at most \(n^{\zeta+m\beta+o(1)}\) records. This
exponent is strictly smaller than \(q\), by the choice of \(p\).
Together with the bound for bad columns it contradicts the common
neighborhood of size \(n^{q-o(1)}\). This proves (FT1) by the limiting
argument, and summation gives (FT2).

An unrestricted limiting endpoint therefore requires each of its
three complete-letter dictionaries to have positive normalized
exponent bounded away from zero. Its largest fixed-value fibers also
have a uniform strict saving in normalized exponent. This leaves
unrestricted endpoints with three large complete dictionaries open.

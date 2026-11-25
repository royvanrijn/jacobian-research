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

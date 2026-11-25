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

## A pooled degree inequality and primitive shift fibres

For two chosen letters \(U,V\), take their normalized prelimit root
measure \(\omega\), with one contribution \(\beta_{T,e}\delta_{eL}\)
per box coordinate. Write \(m=\omega([0,\infty))\),
\(J=\int r\,d\omega\), and let \(F(t)\) be the largest numerical
weight of a submeasure of \(\omega\) of mass \(t\). The cut can split
the mass at its final degree. For any median \(a\) of \(\omega\),
\[
 2F(m/2)-J=\int |r-a|\,d\omega(r).
 \tag{PQ1}
\]
This follows by subtracting the lower-half moment from the upper-half
moment. The contributions of the median cancel.

For a sequence of actual boxes tending to \(d=1/3\), there is
an absolute \(\eta>0\) such that on every pair of letters,
\[
 \liminf \bigl(2F(m/2)-J\bigr)\ge \eta.
 \tag{PQ2}
\]
If this failed, a diagonal sequence of actual boxes would have
\(d\to1/3\) and dispersion tending to zero for one fixed pair.
Indeed, each approximating sequence can be taken far enough to
make its budget and box tolerances arbitrarily small, with the
inner height limits at its selected fixed truncation already taken.
By (GE), \(m\to2/3\). The two taller heights tend to one, since
\(c/2\le\max(a,b)\le c\), so \(1-o(1)\le J\le2+o(1)\)
for every pair. Equation (PQ1) gives \(am=J+o(1)\).
Along a subsequence the medians tend to \(a_*>0\), and for every
fixed \(\rho>0\) the numerical weight outside
\((a_*-\rho,a_*+\rho)\) tends to zero. In fact for \(r\)
in this complement and all sufficiently late boxes,
\(r\le (1+2a/\rho)|r-a|\). A slowly shrinking choice of \(\rho\)
then has the same property. Both letters have height tending to
\(a_*/3>0\), and hence both are finite-scale. This violates
the exclusions following (GE), both when the third letter stays
of positive height and when its height tends to zero. This proves
(PQ2) without any assumption on the total nonatomic mass.

There is an associated count for any single prime-resolved box
selection. Choose full divisors on \(U,V\) from their anatomic
roots, whose total normalized record cost is at most \(t\) and
whose joint numerical weight exceeds the shorter free height by
a fixed positive margin. Full means that if a prime is selected,
its whole indexed power is used. For example, for squarefree
root coordinates one may select subproducts of their primes.
Record the chosen root divisors and the indices specifying their
powers. At fixed truncation the number of records is
\(n^{t+o_X(1)}\) with \(n=X^L\). This remains true for adaptive
selections if their combined root product is at most \(n^t\).
Subdivide the chosen divisors by their dyadic sizes and sum over
the subpower number of size patterns. Conditions in the box
and the height margin tolerate arbitrarily small fixed errors
before the inner limit.

For any fixed signed shift \(U\pm V=D_0\), restrict to pairs
with \(\gcd(U,V)=1\). Their recorded divisors are coprime.
The shorter term is zero modulo its own divisor and is a
known residue modulo the other's divisor. Their product exceeds
that term's height bound, giving at most one pair per record.
In particular the shift need not have a small radical or occur
as a third complete value in the original box. For actual
triples write \(H(T)\) for the entropy of the remaining complete
letter on the uniform edge measure, in natural logarithms.
Reconstruction using \(T\) and the selected record also gives
the box inequality
\[
 \frac{\log |\mathcal B|}{\log n}
 \le \frac{H(T)}{\log n}+t+o_X(1).
 \tag{PQ3}
\]
The count holds on any edge subfamily on which the selection works.

Writing \(Y\) for the pair of divisors and their power indices,
there is also the exact identity
\[
 H(T\mid Y)=H(U,V,T\mid Y)\ge \log|\mathcal B|-H(Y)
 \tag{PQ4}
\]
on the subfamily. In particular some fixed record has at
least \(n^{d-t-o(1)}\) distinct admissible targets. Conversely,
once the locks have combined weight above the short height, a
uniform power-saving bound for the number of admissible
targets per record would give the same saving for the edge
count after paying \(t\). This criterion concerns distinct
targets across the pair of numerical progressions determined
by the record.

There exist absolute \(\xi,\omega_0>0\) and a bounded scaled
degree cutoff \(R\) such that the primitive shift count at
every actual limiting \(1/3\) endpoint, on every pair of free
dictionaries satisfying the next prime-size condition, is at
most \(n^{1/3-\xi+o(1)}\), uniformly in the shift. Require
that each prime in an indexed root of scaled degree at most
\(R\) has normalized logarithm at most \(\omega_0\).
By (PQ2), the upper half has numerical weight at least
\(h_{\min}+\eta/3\) at sufficiently small tolerances.
Concavity and \(J\le 2+o(1)\) bound the slope of \(F\)
uniformly on a fixed neighborhood of \(1/3\).
Choose \(t=1/3-\xi\) with absolute small \(\xi>0\) so
\(F(t)>h_{\min}+\eta/4\). Take \(R\) large enough that
the entire tail above \(R\) costs less than \(t/2\).
Retain it in full. Take the anatomic truncation with its
potentially packed coordinate above \(R\); below \(R\) the
roots are squarefree. Select their primes in descending order
of degree until the next prime would exceed the mass \(t\).
This approximates the fractional cut with numerical loss at
most \(R\omega_0\). Choose \(\omega_0<\eta/(8R)\) and apply
the reconstruction bound. In particular, when the two weak
measures are wholly nonatomic the prime-size condition
holds automatically on a bounded degree band, because even
the largest root there has vanishing normalized mass.

The dispersion term retains escaped first moment. In (PQ2) its
limiting expression is
\(e_U+e_V+\int |r-a|\,d(\nu_U+\nu_V)(r)\) along a
converging subsequence of medians. For example, write
\(D_T=\int (3-r)_+\,d\nu_T(r)\) for the moment of the deficit
below the scaled order three used in the
[powerful triple transfer](POWERFUL_TRIPLES_ABC_ENDPOINT_TRANSFER.md#theorem-a-strict-saving-from-endpoint-diffusion).
At balance, \(e_T+\int|r-3|\,d\nu_T=2D_T+h_T-1\le2D_T\).
The median minimizes the absolute-deviation integral, so (PQ2)
gives \(D_U+D_V\ge\eta/2\) for every pair. This connects the
two-letter count directly to the fullness deficit, with a
uniform moment bound across the entire endpoint.

## A quantitative bound for complete-value dictionaries

Let \(N_\lambda(X;T,\mathcal D)\) count ordered positive primitive
triples \(a+b=c\le X\), \(\operatorname{rad}(abc)<c^\lambda\), whose
specified complete letter \(T\in\{a,b,c\}\) belongs to \(\mathcal D\).
There are absolute \(\sigma_0,\lambda_0>0\) such that uniformly over
\(T\) and every dictionary of at most \(X^{\lambda/800}\) values,
\[
 N_\lambda(X;T,\mathcal D)\ll_{\lambda,\varepsilon}
 X^{(1/3-\sigma_0)\lambda+\varepsilon}
 \qquad(0<\lambda<\lambda_0).
 \tag{TD1}
\]
This gives a specific dictionary exponent in (FT2), independent of
the implicit saving in (FT1).

### The two-free-letter weighted quantile

Use the limiting box convention from (FT1), and suppose the box
count exponent is \(d=d_0=1/3\) in base \(n=X^L\), with a complete
\(T\)-dictionary of size \(n^{\theta+o(1)}\), \(0\le\theta<d_0\).
The necessary conditions (GE) give total root mass \(d_0\) on
each letter. Name the other letters \(A,B\) with \(h_A\le h_B=1\).
The existing short-height exclusion keeps all limiting heights
positive.

Both free letters have arbitrarily cheap positive-weight full
moduli. In fact, if one were finite-scale, the other would admit
such moduli by the conclusions preceding (GE1). Fix \(T\) at cost
\(\theta\), lock a cheap modulus on the other letter, and release
the finite-scale letter using the amplified (GW). This would give
\(d_0\le\theta\). A positive-height letter which is not finite-scale
admits the cheap moduli by the same high-block argument as (GE1).
Now use that proof with \(T\) fixed at cost \(\theta\), leaving only
the other free letter in the congruence. This gives
\[
 c_A,c_B\ge d_0-\theta.
 \tag{TD2}
\]
In particular, the total atomic mass on either free letter is at
most \(\theta\), and every individual coordinate on either letter
has mass at most \(\theta+o(1)\). A sequence of coordinates with
mass greater than \(\theta\) by a fixed margin would admit a
bounded subsequence of scaled degrees by the first-moment bound,
producing an atom of mass greater than \(\theta\).

Let \(\nu=\nu_A+\nu_B\) and \(e=e_A+e_B\) be the combined weak
degree measure and escaped numerical weight. For \(0<t<2d_0\)
write \(W_\nu(t)\) for its highest-weight submeasure of mass \(t\),
measured by the integral of the scaled degree. Define
\[
 W^*(t)=e+W_\nu(t).
\]
This is the limiting numerical weight of the fractional top
\(t\)-mass of the two prelimit spectra. For example, the exact
prelimit upper-quantile formula is
\[
 W_{\nu_{\rm pre}}(t)
  =\inf_{b\ge0}\left(bt+\int(r-b)_+\,d\nu_{\rm pre}(r)\right).
\]
The infimum can be restricted to a compact \(b\)-interval for \(t>0\)
fixed, using the bounded first moment. Replacing the unbounded
integrand by \(r-\min(r,b)\), the objectives converge uniformly on
that compact interval by weak convergence, moment convergence,
and the common Lipschitz bound in \(b\). This yields \(W^*(t)\).

Take the full A and B coordinates in descending scaled degree until
a fractional top selection of mass \(t\) is covered. Rounding up
the last root costs at most \(\theta+o(1)\) by (TD2) and can only
increase the weight. After freezing these full factors \(M_A,M_B\)
and the complete \(T\), the signed equation prescribes
\[
 A\equiv0\pmod {M_A},\qquad A\equiv\pm T\pmod {M_B}.
\]
The two moduli are coprime by primitivity for every nonempty
record. Thus when their product has numerical weight \(>h_A\)
there is at most one positive \(A\) per record. If \(W^*(t)>h_A\)
and \(t+2\theta<d_0\), this would contradict the original count.
Consequently, for \(\theta<d_0/2\), continuity gives
\[
 W^*(d_0-2\theta)\le h_A.
 \tag{TD3}
\]
For \(\theta=0\), take \(t\uparrow d_0\) in this argument.
The increasing concave quantile function has increment from \(t\)
to \(d_0\) at most \((d_0-t)W^*(t)/t\). A median \(a\) of \(\nu\)
splits its total mass \(2d_0\) into two halves, splitting a median
atom if necessary. One then has the exact identity
\[
 2W^*(d_0)=h_A+h_B+e+\int|r-a|\,d\nu(r).
\]
Combining this with (TD3) gives the quantitative shared-scale
obstruction
\[
 (1-h_A)+e+\int|r-a|\,d\nu(r)
       \le \frac{4\theta}{d_0-2\theta}.
 \tag{TD4}
\]
At \(\theta=0\) this alone contradicts (TD2), giving a shorter
completion of (FT1) once the two free letters have been shown
wholly nonatomic. At positive \(\theta\) it forces both free
spectra close to one common scale unless the weighted pair lock
has already saved. For example, on the formal profile with both
free spectra \((1/3)\operatorname{Unif}[2,4]\), root rounding costs
\(o(1)\), and the combined highest \(t\)-mass has weight
\(4t-\tfrac32t^2\). The one-row CRT argument alone therefore
excludes target dictionary exponents
\(\theta<(\sqrt{10}-3)/3\) there.

### A first residue determinant in the remaining shared-scale regime

Suppose \(\theta\le1/800\), so the right hand side of (TD4) is at
most \(6/397\). Since the full numerical height of \(B\) is one,
\[
 |a-3|\le 3\left(e_B+\int|r-a|\,d\nu_B(r)\right)\le18/397.
\]
Outside the open band \((14/5,16/5)\), every degree is at least
\(307/1985\) away from \(a\). It follows that \(\nu_B\) puts mass
at most \(30/307\) outside that band. Its nonatomic mass in the
band is at least \(d_0-30/307-1/800>13/80\).

Choose a compact continuity band inside this open band with
nonatomic mass \(>13/80\). Using only coordinates of individually
vanishing normalized mass, choose a whole-coordinate \(B\)-block
of mass \(\mu=13/80+o(1)\) there. This selection is possible by
discarding all coordinates above successively smaller mass
thresholds, since finitely many coordinates above any fixed
threshold can contribute only atom mass in the weak limit.
The full modulus \(N\) from fixing the block has prime powers of
degrees strictly between \(k\) and \(2k\), where \(kL\to14/5\).
The roots in this bounded scaled-degree band are squarefree, by
using the BBLT anatomic construction with truncation error \(o(L)\).
Its sole potentially nonsquarefree root has scaled degree tending
to infinity. Thus the shallow modulus for the first determinant
is exactly \(Q=\operatorname{rad}(N)\), of normalized mass \(\mu\).

Freeze the complete \(T\) and the A coordinates above a fixed
continuity cutoff greater than \(1000\) on the scaled degrees.
The latter cost at most \(1/1000+o(1)\) by the first moment.
Group all remaining A roots into the one root \(z\), with degree
\(q=\lfloor a/L\rfloor\), and write their powers exactly as
\(v z^q\) with a rational correction \(v\). The normalized narrow
height is at most \(1/3\). The wider numerical height exponent is
at most \(6/397\), since
\[
 \limsup \sum_{e\text{ retained}}|eL-qL|\,\beta_{A,e}
   \le \int|r-a|\,d\nu_A(r)
   \le 6/397.
\]
The compact continuity cutoff justifies the first inequality,
and the height of the reduced correction is at most the sum of
the absolute discrepancies. The A-term lies in the unit residue
class determined by \(T\) modulo \(N\), after the fixed coefficient
and fixed A factors are removed. Primitivity makes all these
divisions legitimate.

Use the first determinant in (QD2), or (DNR1), with
\[
 \alpha=99/100,\qquad M=7,\qquad r_0=7/20.
\]
Its normalized Taylor order is
\(y=\alpha(M+1)r_0=693/250<14/5\), and the row height bound is
\[
 B\ \le\ \frac7{60}+\frac{42}{397}
      \ <\ \frac{9009}{40000}
      \ =\ \frac{(13/80)y}{2}.
\]
The width also satisfies \(r_0<a\). All inequalities are strict.
The determinant supplies a nonzero polynomial of bounded wider
degree and polynomial coefficient height on each shallow class.

The one-letter specialization of (DNR2) applies after fixing \(T\).
Explicitly, let \(J_F\ne0\) be the product of the nonzero
coefficients of the auxiliary polynomial. At a prime outside
\(J_F\) in an A-root of original index \(e\), the valuations of
the monomials, after division by the positive root valuation,
are \(u+(e-q)j\), where \(u,j\) are the narrow and wider indices.
For each \(j\) only the smallest \(u\) can achieve the minimum.
A tie of the minimum between distinct \(j\)'s is necessary for
vanishing, and each possible pair admits at most one exact
degree \(e\). There are at most \(\binom{M+1}{2}\) exceptional
exact degrees, located in a bounded scaled-degree interval.
Their total root mass is at most the atomic mass of \(\nu_A\)
in the limiting count, hence at most \(\theta\). This is uniform
in the shallow class: maximize over all exceptional degree
sets of bounded cardinality and pass to convergent scaled
indices. All other A-position roots are supported on \(J_F\)
and give only \(X^{o_X(1)}\) choices at fixed truncation, by
the fixed-support bound for a polynomial-height integer.
The chosen A roots recover A, and the fixed T then recovers
the triple.

The sum of the costs for the T-dictionary, the B freeze and
shallow class, the A cutoff freeze, and the exceptional A
positions is at most
\[
 \theta+2(13/80)+1/1000+\theta
   \ \le\ \frac{657}{2000}<\frac13.
\]
This contradicts the \(d_0\) endpoint count. The standard limiting
box extraction now proves (TD1): a sequence of counterexample
budgets would give \(D/L\to1/3\) and \(L/\lambda\to1\). Here the logarithm of the height
of the selected dyadic shell is \(1-o(1)\) times that of the
ambient height \(X\). Otherwise the unrestricted one-third
estimate at the shell height, with allowance \(o(\lambda)\)
in the exponent along the budget sequence, would put the shell
strictly below the assumed count at ambient height. Thus
the stated size \(X^{\lambda/800}\) yields the excluded
normalized dictionary cap. This extraction fixes the budget,
truncation, and all degrees before the inner height tends to
infinity, including in the determinant and fixed-support losses.

Uniform complete-letter entropy is at least \(1/800\) in
units of \(\log n\) at any limiting endpoint. If it were less
by a fixed margin, Markov's inequality applied to the information
of the complete value would retain a fixed fraction of the edges
in a dictionary of exponent strictly below \(1/800\). That fraction
preserves the endpoint count and contradicts (TD1).

## A strict saving on every bounded-degree complete-letter curve

For every fixed degree bound \(D\) there are
\(\sigma_D,\lambda_D>0\) such that for \(0<\lambda<\lambda_D\),
\[
 \#\{a+b=c\le X_0:\ \gcd(a,b,c)=1,\
        \operatorname{rad}(abc)<c^\lambda,\ F(A,B)=0\}
 \ll_{D,\lambda,\varepsilon}
 X_0^{(1/3-\sigma_D)\lambda+\varepsilon}.
 \tag{ACD}
\]
Here \(a,b,c\) are positive, \(A,B\) are any chosen distinct
complete letters, and \(F\) is any nonzero rational polynomial
of total degree at most \(D\). The bound is uniform in its
coefficients and permits reducible curves and singular contacts.
Equivalently the curve can be given in any coordinates on the
plane \(a+b=c\) related to a pair of complete letters by an
invertible affine linear change.

### Reducing the equation to polynomial height

Factor over \(\mathbb Q\). A factor depending on one variable
fixes at most \(D\) complete values and is covered by (FT1).
A rational irreducible factor which splits geometrically has
\(O_D(1)\) rational points. The Galois action on its components
is transitive, so a rational point belongs to their intersections.
For a geometrically irreducible factor of degree \(d\) containing
more than \(d^2\) points, evaluation at these points in the monomial
basis of total degree at most \(d\) has one-dimensional kernel.
Indeed a second equation of degree at most \(d\) through the
same points would violate Bézout unless it was a scalar multiple.
Take a maximal-rank square minor after deleting one column.
The signed maximal minors of that row basis recover the equation
up to scale as an integer polynomial of height \(X_0^{O_D(1)}\).
The dyadic shell of a violating sequence has logarithmic height
asymptotic to the ambient logarithmic height, by the unrestricted
one-third count as in (TD1), so this remains a polynomial
height bound in the shell. Below \(F\) denotes this new equation,
with positive degrees \(d_A=\deg_A F\) and \(d_B=\deg_B F\).

### Puiseux residue records and fixed prime support

For \(p\)-adic calculations at primes dividing \(B\), regard \(F\)
as a polynomial in \(A\) over \(\mathbb Q(B)\). It is separable
and its constant coefficient is nonzero in \(\mathbb Q(B)\),
by irreducibility and positive degree in both letters. By the
Newton–Puiseux theorem its \(d_A\) roots are Laurent series
\[
 f_j(B)=\sum_{k\ge k_j}a_{jk}B^{k/q_j}\qquad
 (a_{jk_j}\ne0,\quad 1\le q_j\le d_A).
\]
Use their least ramification denominators. The \(q_j\) rotations
of \(B^{1/q_j}\) give distinct roots, which proves the bound
\(q_j\le d_A\).

Use the quantitative Laurent-series form of Eisenstein's theorem,
Bilu–Borichev, Theorem 6.4. For a polynomial of bounded
bidegree over \(\mathbb Q\) the sum of \(\log p\) over primes
where any Puiseux coefficient fails to be integral at some
extension of \(p\) is \(O_D(1+\log H(F))\), with \(H(F)\)
the maximum absolute coefficient after clearing denominators
and removing content. Apply it to \(F(A,B)\) and the reciprocal
\(A^{d_A}F(1/A,B)\) in \(A\), and to the two analogous
polynomials in \(B\). All are separable in their dependent
variable. Outside a fixed prime support of product
\(X_0^{O_D(1)}\), all coefficients of the root series and of
their reciprocals are integral at every extension of \(p\).
The leading coefficient of each root series is therefore a unit.
[Quantitative Eisenstein theorem, Theorem 6.4](https://doi.org/10.1017/S144678871300013X)

At a prime outside that support where \(v_p(B)\ge e>0\)
and \(A\) is a unit, specialize the formal factorization of
\(F\) over the Laurent Puiseux series. For clarity take a common
denominator of the series and choose a root of \(B\) of that
order for this witness in \(\mathbb C_p\). At each prime fix
one embedding of the coefficient field and the needed roots
of unity. With these choices all series converge, since their
coefficients are integral and the chosen root has positive
valuation. The formal factorization specializes to the same
identity over \(\mathbb C_p\). We may omit the \(O_D(1)\)
triples where \(B\) is a zero of the leading coefficient in
\(A\), since each \(B\) has at most \(D\) possible \(A\)'s
on \(F\). Do the same in the opposite direction.
Consequently \(A\) is represented by at least one of the
\(d_A\) evaluated root series. Its leading term has valuation
\(k_j v_p(B)/q_j\), since its leading coefficient is a unit
and every later term has strictly greater valuation. The unit
condition forces \(k_j=0\). On this branch
\[
       A=f_j(s^{q_j}),\qquad B=s^{q_j},\qquad
       v_p(s)\ge e/d_A ,
 \tag{ACD1}
\]
where \(f_j(s^{q_j})\) is the power series in \(s\) from
the formal expansion, with integral coefficients and fixed
constant term. Thus \(A\) values on the same branch differ by
a multiple of \(p^{\lceil e/d_A\rceil}\) in \(\mathbb Z_p\).
The branch index is one of at most \(d_A\) choices, and
the analogous charts for primes dividing \(A\) use \(d_B\).
Each witness chooses a branch index at each locked prime.
The choice of the common root and hence the parameter \(s\)
may vary with the witness. Each record uses the same coefficient
series for all witnesses marked with the same index.

Include the anatomic exceptional coefficients in the fixed bad
support \(S\). At fixed budget and truncation, all its primes
divide one nonzero integer of polynomial height in the shell.
On either letter, if the \(S\)-supported parts of the anatomic
roots have at least \(\xi>0\) normalized mass, divisor allocation
bounds the count by \(n^{1/3-\xi+o(1)}\) in a limiting
one-third box, with \(n=X^L\) at shell height \(X\). This uses
the bounded fibres over each complete value on \(F\).
In a purported endpoint the supported mass can therefore be
taken to tend to zero. Freeze those supported parts at subpower
cost and let \(A',B'\) be the remaining complete integers,
so \(A=K_A A'\) and \(B=K_B B'\). Dyadically box the residual
roots as usual. The two weak degree measures are preserved
and still have mass \(1/3\). Their numerical heights may lose
escaped first moment. All retained primes see \(K_A,K_B\) as
units. Hence (ACD1) gives integral parametrizations of
\((A',B')\) using the original series and the same ramification
precision, without recomputing the bad support.

### Diffusion on both free letters of the curve

A positive-weight full lock on a residual factor of \(B'\)
prescribes subpower many unit classes of \(A'\) modulo
\(\prod_p p^{\lceil e_p/d_A\rceil}\), where \(p^{e_p}\)
ranges over the locked prime powers. Indeed, choose a branch
index for each prime and combine the classes above. There are
at most \(d_A^{\omega(\prod_p p^{e_p})}=X^{o_X(1)}\) records.
The new modulus has at least \(1/d_A\) times the original
numerical weight and fullness at least \(1/d_A\) times the
original fullness. If one residual letter has zero height or
is finite-scale, its nonatomic mass is zero. The other is
wholly nonatomic by (GE1) applied to the original box and
therefore supplies arbitrarily cheap full moduli in a bounded
positive scaled degree range. Reconstruction at zero height,
or the one-letter amplified (GW) estimate at finite scale,
then leaves only subpower possibilities once a cheap modulus
is fixed. Each such complete value gives at most \(D\) points
on the fixed curve. This contradicts the endpoint. Otherwise
both residual letters supply arbitrarily cheap full moduli
by the high-block argument after (GE). The one-letter version
of (GE1), using the new congruences, shows that both residual
degree measures are wholly nonatomic. These steps take degree
bound, budget, truncation and moduli height parameters fixed
in the inner height limit, before passing to limiting profiles.

### A determinant weighted by both residual heights

Let \(h_A,h_B>0\) denote the residual numerical heights.
Choose full locks on \(A',B'\) with numerical weights \(w_A,w_B\)
at combined normalized root cost \(t<1/3\) such that
\[
 W:=w_A/d_B+w_B/d_A >
       \tfrac12(h_A/d_B+h_B/d_A).
 \tag{ACD2}
\]
To obtain the locks, scale the original degree variable by
\(1/d_B\) on \(A'\) and \(1/d_A\) on \(B'\). The pooled upper
half of this weighted root measure has strictly higher finite
moment than its lower half by nonatomicity. Include the escaped
moments at arbitrarily small root cost as in (TD3). The full
upper-half weight thus strictly exceeds the right side of (ACD2).
Reduce the cost to \(t<1/3\) while retaining a margin, and round
the selection to whole roots at vanishing cost, fixing an
arbitrarily large continuity cutoff for bounded indices and
selecting the high tail in full. The available root values
cost \(n^{t+o(1)}\). Fix those locks and their branch indices.
The combined number of branch records is subpower, and the
two locked prime supports are disjoint by primitivity.

Take real upper height caps \(H_A,H_B>0\) arbitrarily close to
\(h_A,h_B\), so \(A'\le X^{H_A}\), \(B'\le X^{H_B}\).
Order monomials by weight \(iH_A+jH_B\), with a monomial
term order to break ties. Let \(x^u y^v\) lead the residual
equation \(F(K_Ax,K_By)\), and put \(K=uH_A+vH_B\).
Since \(K\ge\max(d_A H_A,d_B H_B)\),
\[
 \frac{H_AH_B}{K}\le
       \min(H_A/d_B,H_B/d_A).
\]
By (ACD2) the caps can be chosen so \(W>H_AH_B/K\) by a
fixed positive margin.

For a large fixed \(T\), use the monomials of weight at most
\(T\) that are not divisible by \(x^u y^v\). They lie in the
disjoint strips \(i<u\) and \(i\ge u,\ j<v\). Counting their
lengths and summing their weights yields
\[
 J_T=\frac{K}{H_AH_B}T+O(1),\qquad
 \sum_{\text{basis monomials}}(iH_A+jH_B)
       =\frac{K}{2H_AH_B}T^2+O(T).
\]
The evaluation determinant on any \(J_T\) witnesses in a fixed
record is an integer. By Leibniz its absolute value is bounded
by \(J_T!\) times \(X\) to the sum of the column weights. At
a locked prime, expand each column in the common integral
parameter series in (ACD1). Its valuation is at least
\[
       \frac{e_p}{d_A}\,\frac{J_T(J_T-1)}2
\]
for a \(B'\)-lock, and the same with \(d_B\) for an \(A'\)-lock.
Indeed the first potentially nonzero multilinear term uses
distinct orders, and all row parameters have the stated
valuation lower bound. This calculation takes place in
\(\mathbb C_p\); its fractional lower bound is valid for the
valuation of the integer determinant. Multiplying these
divisibilities yields exponent at least \(W J_T(J_T-1)/2\).
Since \(W>H_AH_B/K\), for large fixed \(T\) and then large
inner height, the divisor exceeds the absolute-value bound.
Every such determinant vanishes.

The monomial basis consequently gives a nonzero auxiliary
polynomial vanishing on all witnesses of this record. It
cannot be divisible by the residual curve equation since
its leading monomial is outside the corresponding monomial
ideal. Bézout gives boundedly many witnesses, with constants
depending on the chosen margin and the fixed parameters.
The record count is therefore \(n^{t+o(1)}\) with \(t<1/3\),
contrary to saturation. As in (FT1), exclusion of every
limiting one-third sequence with \(D\) fixed gives the
uniform \(\sigma_D,\lambda_D\) in (ACD).

### Implication for covers by curves

For fixed \(D\), a family of at most \(X^{\theta\lambda}\)
curves of degree at most \(D\) gives the same strict bound with
\(\sigma_D-\theta\) in place of \(\sigma_D\) for fixed
\(0\le\theta<\sigma_D\). Thus even a sparse matching with large
complete dictionaries must evade every such cover in the
complete-letter plane. A determinant curve in grouped root
variables does not provide a bounded-degree complete-letter
cover without further arithmetic information.

## An absolute curve saving with arbitrary fixed-budget degrees

There are absolute \(\sigma,\lambda_0>0\) such that at every
fixed \(0<\lambda<\lambda_0\) and every finite degree cap \(D\),
\[
 \#\{a+b=c\le X_0:\ \gcd(a,b,c)=1,\
     \operatorname{rad}(abc)<c^\lambda,\ F(A,B)=0\}
 \ll_{\lambda,D,\varepsilon}
 X_0^{(1/3-\sigma)\lambda+\varepsilon}
 \tag{ACG}
\]
uniformly in every nonzero \(F\in\mathbb Q[A,B]\) of degree at
most \(D\). The degree cap may grow arbitrarily as
\(\lambda\to0\). At a fixed budget it is fixed in the inner
height limit.

By (ACD) we may factor \(F\), discard the bounded fibres
fixing a complete letter by (FT1), and use polynomial coefficient
height on a geometrically integral component with positive
bidegrees. Each step holds with constants depending on the
chosen fixed-budget degree. The new input compared with
(ACD) is an exact divisibility constraint on the valuation
orders that any given branch can accommodate.

### Exact divisibility on Puiseux branches

Write \(d_A=\deg_A F\) and let \(q_j\) be the least ramification
denominators of the Laurent Puiseux roots \(f_j(B)\) in \(A\).
In the common variable \(s=B^{1/M}\), where \(M=\operatorname{lcm}_j q_j\),
rotating \(s\) permutes the roots. An orbit containing branch \(j\)
has size \(q_j\), so the sum of these sizes over distinct orbits
is \(d_A\).

Augment the bad support of (ACD) by the primes dividing \(M\)
and the rational ramification primes of all the coefficient
fields of these series. Bilu–Borichev's Theorem 8.5 applies to
Laurent Puiseux series and bounds the sum of their normalized
logarithmic discriminants by \(O_D(1+\log H(F))\). Each
coefficient field has degree at most \(d_A\), since conjugates
of a series are distinct roots of the same equation when their
coefficients differ. Thus the augmented prime product also has
polynomial height. Apply the same construction in the other
orientation. Discard points where a leading coefficient or
discriminant of the univariate equation vanishes at a witness
value. These are boundedly many by the two bounded projection
fibres. For every remaining witness, specialization as in (ACD)
gives distinct evaluated roots at every good prime dividing the
opposite letter, by the discriminant identity.
[Quantitative Eisenstein and discriminant theorems, Theorems 6.4 and 8.5](https://doi.org/10.1017/S144678871300013X)

At a retained \(p\) with \(e=v_p(B)>0\), take a finite unramified
extension \(E/\mathbb Q_p\) containing all the embedded coefficients
of the series, the \(M\)-th roots of unity, and an \(M\)-th root
of the unit part of \(B\). The last inclusion is possible by
Hensel over a suitable residue field because \(p\nmid M\).
Writing \(r=M/\gcd(M,e)\), the extension \(E(s)/E\) has degree
\(r\), as is seen by the value group, or by Eisenstein for
\(p^{1/r}\). It acts on \(s\) by the cyclic subgroup of order
\(r\) in the \(M\)-th roots of unity. If \(f_j(B)=A\) is rational,
this subgroup fixes its evaluated value. Distinctness of the
specialized roots then forces the corresponding formal branch to
be invariant under that subgroup, since the rotations permute
the formal roots. The stabilizer of branch \(j\) has order
\(M/q_j\). It follows that
\[
                        q_j\mid e .
 \tag{ACG1}
\]
For unit \(A\) the branch has order zero as in (ACD). Its
parameter precision is \(e/q_j\), giving both the integer
residue classes of \(A\) and the Taylor divisibility of (ACD).
The series are convergent in this finite local extension,
and its Galois action commutes with evaluation.

### The weighted lock criterion with varying bidegrees

Suppose a sequence of actual one-third endpoint boxes lies on
such curves. Strip and fix the parts of \(A,B\) supported on
the augmented exceptional primes as in (ACD). Any positive
loss of root mass would strictly save by counting that letter's
residual values, then using its bounded curve fibres. Thus
the residual letters \(A',B'\) keep the weak degree measures
of mass \(d=1/3\). Their numerical heights can lose escaped
moment. At least one residual height has positive limit.
Otherwise both weak measures are supported at zero, in
contradiction with (GE1). Write \(h_{A,j},h_{B,j}\) for the
residual height exponents at outer box \(j\), and set
\[
 \Lambda_j=\max(d_{A,j}h_{A,j},d_{B,j}h_{B,j}),\qquad
 s_{A,j}=d_{A,j}/\Lambda_j,\quad s_{B,j}=d_{B,j}/\Lambda_j.
\]
Pass to a subsequence on which \(s_{T,j}h_{T,j}\to u_T\)
for \(T=A,B\), with \(\max(u_A,u_B)=1\).
Scale the numerical weights of a full lock on \(A'\) by
\(s_{A,j}\) and on \(B'\) by \(s_{B,j}\).

A combined scaled lock weight strictly above \(\min(u_A,u_B)\)
at combined root cost strictly below \(d\) suffices. Indeed,
at each fixed budget use the weighted standard monomials of
(ACD) with height caps \(H_{T,j}>h_{T,j}\) whose changes on
multiplication by \(s_{T,j}\) are arbitrarily small.
Write \(K_j\) for the leading monomial's weight. Multiplying
the determinant precision threshold by
\(d_{A,j}d_{B,j}/\Lambda_j\) gives
\[
 \frac{d_{A,j}d_{B,j}}{\Lambda_j}\,
 \frac{H_{A,j}H_{B,j}}{K_j}
 \le \min(s_{A,j}H_{A,j},s_{B,j}H_{B,j}).
 \tag{ACG2}
\]
This uses \(K_j\ge\max(d_{A,j}H_{A,j},d_{B,j}H_{B,j})\).
The normalized precision from locks is at least their summed
scaled numerical weight, since \(q_j\le d_A\) on the
\(B\)-axis and symmetrically for the other axis. The strict
margin in the original precision can shrink with budget
because the determinant size is chosen before the inner height
limit at that budget. The fixed saving in root cost is what
contradicts the endpoint.

### Selecting whole-coordinate locks

By (GE1), either both selected letters have positive
nonatomic mass, or one has zero nonatomic mass and the
other is wholly nonatomic.

First assume \(c_A,c_B>0\). Both residual limiting heights
are positive, since a positive nonatomic part contains
mass on a compact positive-degree interval and the weak
measures survived the support removal. If \(u_A\ne u_B\),
freeze almost all the roots on the letter of greater scaled
numerical height, leaving out a small fixed mass from a
compact band of positive nonatomic mass there. The rates
in that band stay bounded after scaling, because the
letter's residual height tends to a positive number.
The discarded numerical weight can therefore be chosen
smaller than the height margin. This gives a lock of cost
strictly less than \(d\) and scaled weight greater than
\(\min(u_A,u_B)\).

If \(u_A=u_B=1\), the two scaling factors converge to
positive finite numbers. The scaled weak measures thus
retain their positive nonatomic masses. Choose two compact
bands in these nonatomic parts on opposite letters with
every scaled degree on one band strictly greater than every
scaled degree on the other. Such bands exist because
the nonatomic parts cannot both be supported at one point.
From the whole letter carrying the lower band, remove a
small fixed mass of its lower-band roots and add an equal
mass from the opposite higher band. The substitution has
a strict numerical gain. Remove a smaller additional mass
from the lower band to bring the record cost strictly below
\(d\), while keeping the gain. These operations use whole
coordinates: by shrinking and omitting finitely many atom
windows first, enough compact nonatomic mass remains in
each band, with arbitrarily small maximum individual
coordinate mass. The weighted determinant now applies.
This case imposes no growth restriction on the curve degree.

Next suppose \(c_A=0,c_B=d\), so \(B'\) has a positive
residual limiting height. If \(u_B>u_A\), the same
near-full \(B'\)-lock gives the contradiction. It remains to
consider \(u_A=1,\ u_B=u\le1\).
Assume no individual residual A coordinate tends to mass
\(d\), passing to a subsequence with a fixed gap between
its maximal coordinate mass and \(d\). The scaled A
prelimit measures have mass tending to \(d\) and moment
tending to one. Their weak limit exists along a subsequence
by first-moment boundedness. If it loses a positive scaled
escaped moment, lock an arbitrarily far scaled A tail at
arbitrarily small root cost, carrying a fixed positive
numerical gain. Add the top \(B'\)-roots of mass just below
\(d\). Their scaled moment approaches \(u\) as the cost
approaches \(d\), including any B escape. This beats the
required threshold \(u\) at cost strictly below \(d\).

Otherwise the scaled A weak measure has moment one, hence
mean \(1/d=3\). If it has positive mass strictly above
three, take a compact band of positive mass there. If it
does not, the measure equals \(d\,\delta_3\); then take
a slowly shrinking window about three with mass tending
to \(d\). In either case greedily select roots from this
band to have a limiting mass \(p\in(0,d)\) and scaled
weight at least \(3p\): to get such \(p\), reach a small
fixed positive target mass, and use the fixed gap between
the maximum coordinate mass and \(d\) to control the
overshoot. Extract a subsequence for the actual selected
mass. If \(u=0\), this A block alone suffices.

If \(u>0\), the B scaling factor tends to \(u/h_B>0\),
so its scaled weak measure is wholly nonatomic. The
bottom \(p\) of B mass has scaled weight strictly less
than \(up/d\le3p\): its average is strictly lower than
the total finite-moment average of the nonatomic measure,
and escape can only increase the latter on passing to
the full moment. Hence the top \(d-p\) B mass has weight
strictly greater than \(u-3p\), including its escaped
moment. Add a slightly smaller B mass, retaining the
strictness. This uses whole coordinates at vanishing
rounding cost, by nonatomicity. Thus the weighted
determinant applies also in this case.

For these selections, all weak scaled quantiles and
escaped moments are interpreted by the upper-quantile formula used in
(TD3). They can alternatively be
constructed through compact continuity cutoffs and
slowly growing cutoffs for escape. In the last case
scaled A weights, unlike B weights, need not have been
bounded on compact original degree ranges; first-moment
boundedness holds after scaling by definition.

### Releasing a full-mass coordinate by thinning the branches

The only remaining case is a single coordinate on \(A'\)
of normalized mass tending to \(d\). In particular \(c_A=0\)
by tightness of the original degree measures, so \(B'\) is
wholly nonatomic by (GE1). Fix the complements of that
coordinate on \(A'\) at \(o(1)\) cost, writing \(A=Kx^q\)
with \(x\le X^{(d+o(1))L}\) and \(K\) fixed. Any compact
positive-degree band on \(B'\) with positive nonatomic mass
supplies positive-weight locks of arbitrarily small cost.

If \(d_A L\to0\), use any such \(B'\)-lock of weight bounded
below by a positive number. Its branch records prescribe
subpower many unit classes for \(A\) modulo a product of
prime powers of numerical weight at least \(w/d_A\gg L\).
For any nonempty record \(K\) is a unit at these primes.
The unit power congruence for \(x^q\) has at most
\((2q)^{\omega(N)}\) root classes modulo the product \(N\).
This is the local torsion bound for the units of prime powers,
including powers of two. It is subpower at each fixed
budget and truncation. The modulus exceeds the height of
\(x\), so the \(B'\)-lock records count the triples at
cost below \(d\), a contradiction.

Otherwise pass to a subsequence with \(d_A L\ge c_0>0\).
Choose a compact positive continuity band on \(B'\) and
a strictly sub-\(d\) root mass within it whose numerical
weight is at least some \(w_0>0\). For any fixed \(\theta>0\)
there are at most \(1/\theta\) rotation orbits on \(B=0\)
with ramification order \(q_j\ge\theta d_A\), since the
orbit sizes sum to \(d_A\). By (ACG1), every good prime
using one of these orbits has original degree \(e\) with
\(eL\) a positive multiple of a spacing \(q_j L\ge\theta c_0\).
Within the chosen compact band this is a bounded number
of locations. Along a further subsequence these locations
have a finite set of limits, and wholly nonatomicity of
\(B'\) makes the corresponding total root mass tend to zero.
Discard all those original degree positions. The bounded
band can be chosen below the packed anatomic index, and its
retained root factors are squarefree. After the support
stripping each prime in these positions has precisely its
original degree as valuation on \(B\). Thus every prime
in a \(B'\)-lock on the remaining positions has \(q_j<\theta d_A\).
We can choose the lock with weight at least \(w_0/2\) and
cost strictly below \(d\). Its precision is at least
\(w_0/(2\theta d_A)\). Choose \(\theta\) small enough that this
exceeds \(H_B/d_A\) by a fixed factor, with \(H_B\) a numerical
height upper cap on \(B'\) close to its limiting height.
The weighted determinant threshold is at most \(H_B/d_A\),
and the records again give a strict saving. This covers
arbitrary growth of the fixed-budget degrees.

If the stated absolute gap and budget threshold did not
exist, choose budgets tending to zero, gaps tending to zero,
and finite degree bounds witnessing failure. At each fixed
budget and degree bound take an unbounded inner height
sequence, with anatomic losses as small as required, before
passing to the next budget. The unrestricted one-third
bound and the usual shell extraction give \(L/\lambda\to1\)
and \(d\to1/3\). The preceding proof excludes that sequence.
This proves (ACG) with a uniform exponent and the stated
degree-dependent constants. It also excludes covers by at
most \(X^{\theta\lambda}\) such curves for any fixed
\(\theta<\sigma\), by summing the uniform bound.


### Degrees bounded by a power of the logarithm of height

There are absolute \(c_*,\sigma_*,\lambda_*>0\) such that
the bound in (ACG), with \(\sigma_*\) in the exponent,
also holds uniformly for all nonzero equations of degree
at most \((\log X_0)^{c_*}\) at height at most \(X_0\).
The implicit constant may be taken to depend on \(\lambda,\varepsilon\).
This allows the curve to vary with the inner height.

The polynomial-height reduction is quantitative. For degree \(D\),
maximal minors use \(O(D^2)\) columns of monomial values bounded
by \(X_0^D\), so \(\log H(F)\ll D^3\log X_0\). The factor
selection and bounded projection fibres cost polynomial factors
in \(D\). At an endpoint the shell height \(X\) has
\(\log X\sim\log X_0\). Bilu–Borichev Theorems 6.4 and 8.5
give bad-prime log product \(O(D^7\log X)\) in both orientations,
including the reciprocal integrality and field ramification.
In the latter estimate the discriminant degree is \(O(D^2)\),
the factor in Theorem 8.5 is \(O(D)\) times that degree, and
unnormalizing the field discriminants costs at most one factor
\(D\). Primes dividing any Puiseux denominator are at most \(D\).
The excluded discriminant specializations give polynomially
many points in \(D\) and hence subpower loss.

For \(D\le C(\log X)^c\), the number of exceptional primes is
\(O((\log X)^{1+7c})\). Fix \(0<c<1/20\) and set \(\eta=10c\).
By ordering the distinct primes, Euler–Rankin gives
\[
 \log\prod_{p\in S}(1-p^{-\eta})^{-1}
 \ll_\eta (1+\omega(S))^{1-\eta}=o_c(\log X)
\]
since \((1+7c)(1-\eta)<1\). Apply it to the single product
of the unpowered supported parts of the original roots,
bounded by \(X^{L+o_X(1)}\) in the box. Allocation among
the finitely many root positions at fixed truncation is
divisor-many. Thus fixing all supported factors costs at
most \(10c+o_X(1)\) in normalized units. In the overlap
exclusion it suffices for the normalized root mass removed
on a letter to exceed this cost by a positive margin.
In a contradiction sequence with \(c\to0\) the original
weak measures survive as before. For the branch records
on a lock of cost \(t\), its underlying radical is
at most \(X^{(t+o_X(1))L}\). The elementary bound for the
number of its distinct primes gives
\[
  \log D\cdot\omega(N_A N_B)\le (ct+o_X(1))L\log X .
\]
The unit power-congruence records on the sole remaining
coordinate still cost subpower, because its anatomic power
is fixed at the inner height limit. All these estimates
are uniform at fixed inner budget, truncation and \(c\).

For the determinant in (ACG2), use caps
\(H_i=(\Lambda/D_i)(u_i+\epsilon)\), with the prelimit
\(u_i=D_i h_i'/\Lambda\) and small fixed outer \(\epsilon>0\).
At least one ordinary residual height has a fixed positive
lower bound along the endpoint sequence. Thus
\(C_\epsilon/D\le H_i\le2D\) for a positive
\(C_\epsilon\) independent of inner height.
The strict normalized locking margin gives
\(W>(1+\gamma)H_AH_B/K\) with fixed outer \(\gamma>0\).
For the \(B'\)-only argument choose \(H_i=h_i'+\epsilon\)
and use the same inequality via \(H_AH_B/K\le H_B/d_A\).
Put \(C_0=K/(H_AH_B)\ge1/(2D)\). Counting the two
standard-monomial strips by arithmetic progressions and
summing their weights give
\[
 J_T=C_0T+O_\epsilon(D^4),\qquad
 S_T=C_0T^2/2+O_\epsilon(DT+D^6)
\]
for weight cutoff \(T\ge K\), where \(S_T\) is the sum of
column weights. Indeed each strip uses at most \(D\) rays,
with initial weight \(O(D^2)\) and step at least
\(C_\epsilon/D\). With \(T=C_{\epsilon,\gamma}(1+D)^{10}\)
the Taylor divisor beats the product of column bounds by a
fixed multiplicative margin in logarithmic exponents.
The Leibniz factor \(J_T!\) is absorbed because \(J_T\) has
polynomial size in \(D\) and \(\log J_T=O(\log\log X)\).
Bézout then gives polynomially many points per record in \(D\).

If the asserted constants did not exist, choose \(c\to0\)
alongside prospective exponent gaps and budgets tending to
zero. At each fixed inner budget and \(c\), the estimates
above are uniform even for degree varying with height.
Take the unbounded height witnesses beyond the inner
error thresholds. The exceptional support and branch
record losses tend to zero in normalized units along the
outer sequence. The strict lock margins and exact branch
divisibility in (ACG1) therefore exclude the endpoint by
the same whole-coordinate selection and thinning argument.

## Affordable whole-coordinate pair locks at every balanced endpoint

At every actual limiting \(abc\) box with normalized edge
exponent \(d=1/3\), at least one pair of letters admits
two full divisor locks of combined normalized root cost
strictly below \(1/3\) and combined numerical weight
strictly above the shorter numerical height of those
two letters. The construction permits primes of any normalized size. In fact the same
conclusion holds on any pair on which both limiting
nonatomic root masses are positive.

There are absolute positive margins for at least one
pair in actual balanced endpoint sequences. Thus the
affordable-lock hypothesis in the fixed-record target
reduction can be imposed for every endpoint, including
when macroscopic primes prevent rounding the fractional
top quantile. A distinct-target count for the resulting
record fibres is still separately needed.

### Small full blocks in diffuse bands

Name the chosen pair \(U,V\), with normalized root
masses \(l_U=l_V=1/3\) and numerical height limits
\(h_U,h_V\). If a compact positive continuity band
has positive nonatomic mass on a letter, then any
fixed smaller mass can be selected there as a union
of whole original root coordinates, with a vanishing
normalized rounding error. To justify this
independently of atoms in the band, for any fixed
cutoff \(\epsilon>0\) the coordinates of mass at least
\(\epsilon\) are bounded in number. Any weak limit of
their restriction to the band is purely atomic.
Discarding them therefore leaves at least the
nonatomic mass in the band up to vanishing errors.
Then take the cutoff slowly to zero. Greedy addition
of the remaining whole original positions can attain
the desired mass at vanishing overshoot. The same
observation allows deletion of such a prescribed
block from the full product on a letter. The
selected indices are fixed within each sufficiently
fine box, and the record count is bounded by their
combined normalized root cost, with only the usual
fixed-truncation loss.

### Removing or exchanging a diffuse block

If \(h_U>h_V\), start with the whole product of
the original power factors on \(U\). It has cost
\(1/3+o(1)\) and numerical weight \(h_U+o(1)\).
Delete a small fixed positive-mass block from a
compact diffuse band of \(U\), with upper degree
\(a<\infty\). Choose its mass \(s>0\) sufficiently
small that \(a s<h_U-h_V\), with room for
rounding. The residual full divisor on \(U\) then
has mass \(1/3-s+o(1)\) and weight strictly above
\(h_V\). The divisor on \(V\) can be one. The same
argument applies with the letters reversed if
\(h_V>h_U\). If a nontrivial lock on the second
letter is desired, add an arbitrarily small
positive-weight block from its diffuse part,
keeping the strict cost margin.

Now suppose \(h_U=h_V\). There are separated
positive compact continuity bands of nonatomic
mass, one of upper degree at most \(a\) on one
of the letters and one of lower degree at least
\(b>a\) on the other. To see this, if no pair of
positive support points on opposite letters can
be strictly ordered in either direction, both
nonatomic measures must be supported at the same
single point, impossible. Choose the lower band
on \(U\) by interchanging letters if necessary.

Choose \(s>0\) smaller than the nonatomic masses
of both bands and choose \(s'>0\) with
\[
                     (a/b)s<s'<s .
\]
Start again with the whole original power product
on \(U\). Delete a block of root mass \(s+o(1)\)
in its lower band and use the remainder as the
first full divisor. As the second full divisor
take a block of mass \(s'+o(1)\) in the higher
band on \(V\). The combined normalized root cost
is
\[
             \frac13-s+s'+o(1)<\frac13
\]
and its numerical weight is at least
\[
             h_U-a s+b s'+o(1)>h_U=h_V .
\]
Every retained macroscopic prime in \(U\) may be
kept with its entire original power. No
macroscopic prime is split into fractional mass.

This counts the actual full height of the starting
letter, including any escaped numerical moment:
the deleted block lies in a bounded scaled-degree
band, so its numerical weight is precisely the
bounded loss estimated above. Fixed exceptional
coefficients can be included where needed at
vanishing normalized cost. Only the whole root
positions actually selected are enumerated.

### Existence and uniform margin in actual boxes

By (GE1) and its zero-nonatomic alternatives,
every endpoint pair has nonatomic mass at least
\(1/3\) in total. In particular at least two
letters have positive nonatomic mass, so the
preceding construction is always available for
some pair.

To obtain absolute margins for actual endpoint
sequences, suppose there were no absolute
\(\eta>0\) such that some pair has a prescribed
whole-coordinate lock of cost at most
\(1/3-\eta\) and weight at least the shorter
height plus \(\eta\), up to vanishing errors.
Choose a diagonal sequence of actual populated
boxes tending to the endpoint in which every
pair fails these conditions at tolerances
\(\eta_j\to0\). As in (PQ2), fix each inner height,
box and truncation loss before passing to the
outer sequence. Pass to limiting heights and
weak measures. By (GE1) two letters still have
positive nonatomic mass on the limiting sequence.
The preceding construction now chooses fixed
bands and fixed positive strict margins on that
pair. Weak convergence on continuity bands, the
small-coordinate selection argument, and the
convergence of the full numerical heights
realize the same margins on sufficiently late
prelimit boxes. This contradicts the choice
of the sequence.

The argument does not require nonatomic mass to
be lower semicontinuous under a separate weak
limit of profiles. The endpoint hypotheses are
applied to the final diagonal sequence of actual
boxes, where (GE1) again supplies its positive
nonatomic masses.

### Fixed-record and fixed-shift consequence

For the pair \(U,V\) found this way write \(t\)
for the combined record cost and \(Y\) for its
two full divisors, including their prescribed
indices if necessary. The record count is at
most \(n^{t+o(1)}\). On a nonempty record the
divisors are coprime by primitivity. If the
third complete value \(T\) is fixed, the shorter
term is prescribed modulo their product by
the signed relation \(U\pm V=\pm T\).
The product exceeds the shorter height, so
there is at most one primitive pair per record.

It follows in particular that a primitive
fixed-shift count for the chosen pair has
an absolute normalized gap from \(1/3\),
independent of any prime-size assumption.
More importantly, the fixed-record
distinct-target obstruction (PQ4) now applies
on at least one pair in every balanced endpoint
sequence, without a separate affordable-record
hypothesis for the macroscopic-prime sector.
The remaining estimate concerns the number of
actual small-radical targets in that record.


### The exceptional form for an individually prescribed pair

For completeness one can characterize failure on a
prescribed pair, allowing prime-resolved full locks
in this paragraph. Suppose the pair is \(U,V\) and
one of them, say \(U\), has zero nonatomic mass.
The zero-mass alternative from (GE1) says that \(V\)
is then wholly nonatomic. If \(h_U<h_V\), use almost
all of the taller diffuse \(V\) as before. Thus
suppose \(h_U\ge h_V\).

If \(U\) has a positive escaped first moment, choose
its full high-degree tail at arbitrarily small
root cost and fixed positive numerical weight.
Keep almost all of \(V\), deleting a small bounded
diffuse block of strictly larger mass than the
tail cost, but numerically smaller than the tail
weight. This gives strict margins. Hence a failing
pair must have no escaped moment on \(U\).

Write \(l=1/3\). The purely atomic degree measure
on \(U\) now has mean \(h_U/l\ge h_V/l\).
If it is not a single mass-\(l\) atom, it has some
atom of mass \(0<\beta<l\) at a degree
\(r\ge h_V/l\). Keep the entire \(U\) window for
this atom. Its cost and weight tend to \(\beta\)
and \(r\beta\). Let \(B_V(v)\) be the bottom
numerical moment of the diffuse \(V\) measure
at mass \(v\). Its bottom \(\beta\)-mean is
strictly below its full mean, because it is
nonatomic and \(0<\beta<l\). In particular
\[
     B_V(\beta)<\beta h_V/l\le r\beta ,
\]
where the upper bound still holds if \(V\) has
its own escaped moment. Choose a small
\(\epsilon>0\) with \(\beta+\epsilon<l\) and
\(B_V(\beta+\epsilon)<r\beta\) by continuity.
Remove the bottom \(\beta+\epsilon\) of \(V\)
at vanishing overshoot, using bounded continuity
cutoffs and whole small root positions. The
retained \(V\) divisor includes any escaped
moment. Together the \(U\) window and the
\(V\) remainder have cost \(l-\epsilon+o(1)\)
and weight strictly above \(h_V\).

The only case left is
\(\nu_U=l\delta_r\), with \(r=h_U/l\), no escape.
If every individual prime in the shrinking
\(U\) window has normalized logarithm at most
\(l-\gamma\) for a fixed \(\gamma>0\), a positive
proper subfactor of its roots can be chosen
with total mass bounded away from zero and \(l\).
Indeed take one prime of mass at least
\(\gamma/4\) if available, or add primes of
smaller mass until first reaching \(\gamma/4\).
The mass \(\beta\) obtained lies in a compact
subinterval of \((0,l)\). At its original degree
powers the factor weight is \(r\beta+o(1)\).
The same bottom-mass exchange applies, with a
uniform small \(\epsilon\) because
\(r\beta-B_V(\beta)>0\) throughout this compact
interval. If necessary bin the adaptive divisor
size before enumeration. Counting all its
possible prime-resolved records incurs only the
root cost of the selected subfactor, as in (PQ3).

Consequently a failing subsequence of pairs and
prime configurations must have a single prime
on \(U\) whose normalized logarithm tends to
\(l\), and whose numerical power carries the
entire height \(h_U\ge h_V\). Conversely in this
configuration all selections of total cost at
most \(l-\eta\), for fixed \(\eta>0\), must omit
that prime. The other \(U\)-factors have numerical
weight tending to zero and the \(V\) lock has
weight at most \(h_V+o(1)\). Thus a fixed positive
surplus over the shorter height is impossible
on this prescribed pair. The other two letters
are fully diffuse and their pair remains
affordable by the earlier argument.

## Local sparsity of targets over a fixed full-divisor pair

Use the normalized \(abc\) box scale \(n=X^L\). Consider
actual balanced endpoint sequences, and fix two coprime
full divisor records \(F\mid U\) and \(G\mid V\).
Let \(w_F,w_G\) be their numerical weights. Suppose
\[
             w_F+w_G-\min(h_U,h_V)\ge\eta>0
\]
with fixed positive margin. This is the surplus in the
affordable pair-lock reduction. Then for every fixed
normalized allowance \(\epsilon>0\) there is a
positive numerical window exponent \(\delta>0\)
such that every interval of length \(X^\delta\) contains
at most \(n^{\epsilon+o(1)}\) distinct admissible targets
from the fixed record. In the iterated endpoint limit this bound is uniform across
the endpoint sequences and the record choices with
the stated surplus. It does not restrict the number
of separated intervals met by the full target fibre.

### Affine collinearity in one short interval

Write \(U=Fu,\ V=Gv\), and permute the names so \(V\)
has the shorter complete numerical height. In a
fixed interval of target values, fix one reference
point \((u_0,v_0,T_0)\), and write
\(\Delta u_i=u_i-u_0\), \(\Delta v_i=v_i-v_0\),
\(\Delta T_i=T_i-T_0\). Absorb the fixed signs into
the linear equation for \(T\). Then
\[
      F(\Delta u_1\Delta v_2-\Delta u_2\Delta v_1)
           =\pm(\Delta T_1\Delta v_2-
                         \Delta T_2\Delta v_1).
\]
The right side has absolute size at most
\(X^{\delta+h_V-w_G+o(1)}\) up to a constant.
If \(\delta<\eta\), with a fixed margin for the
height tolerances, this is smaller than \(F\).
Consequently all difference pairs are proportional.
The ordinary fixed-record CRT uniqueness ensures
that different pairs give different \(T\).

If at least two pairs occur, the integer points on
their affine line have a primitive integral
direction \((a,b)\). The target values are therefore
contained in one arithmetic progression with
integer step \(Q=\pm F a\pm G b\ne0\). Its index
diameter is at most \(X^\delta/|Q|\). If fewer than
two occur the claimed count is immediate. In the
nontrivial case \(|Q|\le X^\delta\), so \(Q\) has
polynomial height uniformly.

### Arithmetic count on that progression

It suffices first to prove the sequential statement
for intervals of widths \(X^{\delta_j}\) with
\(\delta_j\to0\). Strip and fix the \(Q\)-supported
part of every original target root at subpower
cost by the fixed-support divisor count. For \(|Q|\le X^C\), a number of size \(X^{O(1)}\) with prime support in \(Q\) can be chosen in \(X^{o_X(1)}\) ways, uniformly in \(Q\), by [the supported divisor count](ABC_EXCEPTIONAL_SET_AND_POWERFUL_TRIPLES.md#a-radical-count-contradiction-to-the-claimed-smooth-number-lower-bound).
The exceptional coefficient can be fixed separately
as in the original box. Refine the residual root
scales and pass to a limiting residual degree
measure, whose first moments remain bounded.

If that measure has positive nonatomic mass, choose
a compact positive band containing a positive
amount of it. It supplies a residual full block
of an arbitrarily small fixed root cost and fixed
positive numerical weight. The lock uses whole
small root positions of the residual band after
removing the atomic windows, and is coprime to
\(Q\). The requirement that the target is zero
modulo the lock determines its progression index
modulo the lock. Its positive numerical weight
exceeds \(\delta_j\) eventually. Therefore each
lock record gives at most one target, counted at
the arbitrarily small fixed root cost.

If the residual measure is purely atomic, choose
a finite list of its atoms capturing all but an
arbitrarily small root mass. Keep slowly shrinking
windows at its positive selected atoms and at zero.
Freeze the other positions at arbitrarily small
normalized cost. This includes any escaped
numerical first moment, since it has vanishing
root mass. Group each retained positive window at
its fixed representative scaled degree. The
exponent corrections and the zero window go into
a wider rational factor of numerical height
\(X^{o(1)}\). All fixed factors, including the
stripped supports, multiply the exact odds by a
constant.

The numerical height of the complete target is
bounded below at the endpoint by the existing
short-height exclusion in [the real gap bound](ABC_SHORT_HEIGHT_AND_REAL_GAP.md#a-strict-saving-when-the-short-term-is-small). Taking a value in the
interval as the center, its relative width is
\(X^{-\tau}\) for a fixed \(\tau>0\) once
\(\delta_j\to0\). Apply the unequal-height
real strip [RUW1](ABC_SHORT_HEIGHT_AND_REAL_GAP.md#a-real-strip-with-unequal-root-heights) to the grouped target with
the fixed factor absorbed into the strip center.
Choose its normalized counting tolerance small
first, then use that the wider exponent tends
to zero. Multiplicity at a fixed exact target
is divisor-many at fixed truncation, or can be
made one by choosing a witness representation.
This yields arbitrarily small normalized cost
also in the purely atomic case.

For the uniform positive-width statement, suppose
a fixed \(\epsilon>0\) is contradicted for every
positive width exponent. Choose a diagonal
sequence of actual endpoint boxes and fixed
records at widths with \(\delta_j\to0\), still
violating the bound at that \(\epsilon\), taking
the inner height limits first. The preceding
argument applies to that sequence after fixing
and refining the support records at subpower
cost. This is a contradiction. The two cases
cover all residual spectra, including when
the original target was a macroscopic-prime
atomic letter.


### Numerical interval entropy

For a uniform edge in the box write \(I=\lfloor T/X^\delta\rfloor\)
and let \(Y\) be its affordable pair record of normalized cost \(t\).
A pair \((Y,T)\) determines the edge by (PQ4). If the family has
\(n^{d+o(1)}\) edges, the interval bound yields
\[
 H(I\mid Y)\ge (d-t-\epsilon-o(1))\log n.
 \tag{PQ5}
\]
Indeed \(H(T\mid Y)=H(\mathrm{edge}\mid Y)\ge(d-t-o(1))\log n\),
and \(H(T\mid I,Y)\le(\epsilon+o(1))\log n\). Equivalently, there
must be at least \(n^{d-t-\epsilon-o(1)}\) occupied intervals of
width \(X^\delta\) on at least one record. The quantitative endpoint obstruction is a bound on
the separated intervals met by actual affordable fibres. An
aggregate estimate suffices. If \(\mathcal I_\delta(Y)\) is
the set of occupied interval indices on record \(Y\), then
\[
  \sum_Y|\mathcal I_\delta(Y)|\ \le\ \#\mathrm{edges}\
       \le\ n^{\epsilon+o(1)}\sum_Y|\mathcal I_\delta(Y)|.
\]
To close the reduction it is enough to bound this sum by
\(n^{1/3-\sigma+o(1)}\) for an absolute \(\sigma>0\) with a
window chosen for \(\epsilon<\sigma\). This allows the sizes
of the individual fibres to vary.

The [affine rank argument](ABC_AFFORDABLE_FIBRE_AFFINE_RANK.md)
extends fixed-record sparsity to short windows along polynomial
progressions and other bounded integer parameter systems at
every balanced endpoint. A further full target lock makes a
wider interval fibre affine in a different parameter range.
Its shifted line either provides separate residue classes or
transfers a full modulus to the intercept. The resulting
occupancy is controlled by the smaller residual nonatomic mass.
[Affine count for sparse interval fibres](ABC_SPARSE_AFFINE_CONTENT.md)

The general separated-target count still requires a uniform
cover or a different global bound.

## Progressions in a wholly diffuse target dictionary

A wholly nonatomic target measure gives sparsity in an ordinary
arithmetic progression of arbitrary step with index diameter
at most \(n^C\), for any fixed \(C\). The intersection has size
\(n^{o(1)}\) in the iterated endpoint sense. If the progression
contains two target values of polynomial height, its nonzero
step has polynomial height. Fix and strip its supported target
parts as above. On any subsequence the residual measure is
dominated by the original target measure, so it is nonatomic.
If its mass tends to zero, enumerating all residual coordinates
is subpower in \(n\). Otherwise a fixed positive-weight full
block at arbitrarily small cost determines the progression index
as soon as its size exceeds \(n^C\). Such a block exists by the
same positive-band argument, since \(X^w>n^C\) for every positive
fixed numerical weight \(w\) once \(L\) is sufficiently small.
The argument applies to any subfamily of the dictionary.

There is an additive consequence for a set \(\mathcal S\) of
such target values of size \(n^{r+o(1)}\) with \(r>0\) fixed.
Its doubling tends to infinity in the iterated limit. If the
doubling were bounded on a subsequence, the inverse theorem
would contain \(\mathcal S\) in a proper generalized arithmetic
progression of bounded positive rank \(R\) and size
\(O(|\mathcal S|)\). The longest coordinate has length at
least \(|\mathcal S|^{1/R}\). Fixing the other coordinates gives
\(O(|\mathcal S|^{1-1/R})\) slices of index diameter at most
\(n^{r+o(1)}\). Each meets the target dictionary in \(n^{o(1)}\)
values, a contradiction. By the Balog–Szemerédi–Gowers theorem,
also \(E_+(\mathcal S)=o(|\mathcal S|^3)\), where \(E_+\) counts
ordered equal-sum quadruples. A positive fraction of cubic
energy would give a subset of comparable size with bounded
difference growth and hence bounded doubling by Plünnecke's inequality
applied to \(\mathcal S',-\mathcal S'\). The new subset has the same
progression bound. These conclusions
supply no power bound on either the doubling or the energy.
[Green–Ruzsa, Theorem 1.1](https://arxiv.org/abs/math/0505198)
· [Reiher–Schoen, Theorem 1.1](https://arxiv.org/abs/2308.10245)
· [Petridis, Theorem 1.1](https://arxiv.org/abs/1101.3507)

### Progression rigidity of the three complete letters

The cofactor determinant also works on an arbitrary-step
progression \(T=T_0+Qj\) of index diameter \(H_0\), on a fixed
affordable record. Reduce the signed equation modulo \(F\).
If the progression indices are measured from a chosen solution,
\(G\Delta v\equiv \pm Q\Delta j\bmod F\) with a fixed sign.
Thus \(F\) divides \(\Delta j_1\Delta v_2-\Delta j_2\Delta v_1\).
When \(H_0\ll n^C\) at fixed \(C\), its absolute value is
smaller than \(F\) by the numerical surplus above the shorter
letter. It must vanish, so \(v\) is an affine function of the
target index on the sample. The exact relation then gives the
same conclusion for \(u\). Every moving complete letter in the
sample lies in an arithmetic progression of index diameter
at most \(H_0\). To see this even for a rational slope, put it
in lowest terms; admissible indices occupy one class modulo
its denominator. Reindexing makes the step integral.

Consequently, when at least two of the three original letters
have wholly nonatomic degree measures, the progression bound
holds for the fixed-record targets even if the target itself
has atoms. Two distinct triples have different values in at
least two letters by the signed relation. Thus some wholly
diffuse letter moves and counts the sample injectively. The
bounded-doubling and normalized-energy conclusions then apply
also to this class of fixed-record targets.

## A diffuse target family with no cheap heavy-lock record

This example imposes genuine small-radical, indexed-power
arithmetic on the target values themselves. It isolates
a limitation of trying to control a fixed-record fibre
using only target marginal entropy and the sizes of
its prime coordinates. It does not realize the targets
as sums or differences in a common primitive two-divisor
fibre.

Fix a desired sample exponent \(0<r<1/3\).
At each fixed large integer \(m\), set \(L=1/m\),
\(n=X^L\), and \(J=2m-1\). Use degree positions
\(e=2m+1,\ldots,4m-1\) and give every root the same
normalized box mass \(\beta=1/(3J)\). Set \(H=n^\beta\).
Choose disjoint prime pools in \([H/4,H/2]\) for the
\(J\) positions, each with \(q\) primes, where \(q\)
is a power of two and
\[
                 q=n^{\beta-o_X(1)}.
\]
This is possible by the prime number theorem with
\(m\) fixed in the height limit. Fix a bijection
between each prime pool and the elements of
\(\mathbb F_q\). Use \(J\) distinct field evaluation
points \(a_e\). For \(k=\lfloor r/\beta\rfloor\)
and each polynomial \(f\in\mathbb F_q[z]\) of degree
less than \(k\), form the integer
\[
        T_f=\prod_e p_e\bigl(f(a_e)\bigr)^e ,
\]
where \(p_e\) is the prime assignment for position
\(e\). Take \(m\) large enough that \(1\le k<J\).

The target values are distinct: their prime
factorizations recover the evaluation vectors,
and \(J\ge k\) distinct evaluations determine
\(f\). All lie at the full global height
\(T_f\asymp_m X\), and in particular below \(X\)
for large \(X\), since
\[
   \sum_e \beta=\frac13,\qquad
   \sum_e eL\,\beta=1.
\]
The interval \([H/4,H/2]\) changes the products by
multiplicative constants depending on \(m\), not on
\(X\). Similarly \(\operatorname{rad}(T_f)\asymp_m
n^{1/3}\), with the upper bound \(n^{1/3}\).
The weak indexed root law tends to the wholly
nonatomic \(\frac13\operatorname{Unif}[2,4]\).
The family size is
\(q^k=n^{r_m-o_X(1)}\), where \(r_m=k\beta\to r\).
These claims take the height limit first at each
fixed \(m\).

Every selected full prime factor is one of the
individual \(p_e(f(a_e))^e\). For any fixed record on
\(j\) positions, its fibre within this family has
exactly \(q^{k-j}\) elements when \(j\le k\), by
polynomial interpolation on distinct evaluation
points. For \(j\ge k\) the record determines
\(f\) if it occurs at all. Thus every subset of
\(j\le k\) root coordinates has full entropy
\(j\log q\) on the uniform family. Larger choices
have the full joint entropy \(k\log q\).

The same fibre bound controls adaptive choices.
On any subfamily of size \(n^{r_m-o_X(1)}\), if every
record chooses at least \(j_0\) positions, its record
dictionary has size at least
\[
          n^{\beta\min(j_0,k)-o_X(1)}.
\]
The positions are included in each record, and
no independence of the adaptive choice is needed:
each possible record has at most
\(q^{k-\min(j_0,k)}\) preimages even in the whole
family. In particular a target lock of numerical
weight strictly above \(4r_m\), with fixed margin,
must use more than \(k\) positions, since each
position contributes at most \(4\beta+o_X(1)\).
It specifies at most one target. Its dictionary
therefore has the entire \(n^{r_m-o_X(1)}\)
exponent, even though it may be chosen adaptively.

Thus a conditional-marginal entropy argument alone cannot ensure
a dictionary cost below \(r\) for a target lock whose numerical
weight exceeds \(4r\) by a fixed margin. A primitive additive
realization with common fixed coefficients would require
additional arithmetic.

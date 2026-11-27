# Multiplying the two internal products of a six-cycle

The product \(J_1J_2\) of the internal root ratios on two paths at opposite endpoints fixes a congruence for the two outer roots at the second endpoint, once the two edges at the first endpoint are known. The corresponding weighted root product is strictly less than the target. Divisor recovery then counts every cycle with a prescribed product, including all reciprocal internal products. Equal internal products have a different consequence. The two endpoint root ratios agree, which gives a binary norm equation once that common ratio is prescribed.

Use the [selected actual coprime graph](ERDOS_943_PAIRED_THREE_EDGE_PRODUCTS.md#the-endpoint-equation) at \(n\), with fixed positive coefficients \(\kappa_1,\kappa_2\), kernel caps \(K_1,K_2\), root caps \(R_1,R_2\), and a subpower actual codegree bound. The [critical rectangle](ERDOS_943_CRITICAL_RECTANGULAR_CYCLE_CRITERION.md) has these properties and \(R_1R_2<\sqrt n\) for sufficiently large \(n\). Each of the two paths is oriented from the same left endpoint to the same right endpoint. Denote the common endpoint weighted kernels by \(B_0,B_3\). On path \(i\), write \(X_i,Y_i\) for the endpoint roots, \(J_i=p_i/q_i\) for its reduced internal product, and \(A_i=q_i t_i\) for its positive integral norm coordinate, where \(t_i\) is the second internal ratio. In this notation
\[
         B_0p_i^2 X_i^2+B_3q_i^2Y_i^2
                         =n(p_i^2+q_i^2-A_i^2).                   \tag{1}
\]
The ratios and all selected roots are units modulo \(n\). Counts below fix \(n,\kappa_1,\kappa_2\) and concern simple six-cycles.

## Recovery from a product of two paths

Set \(K=\min(K_1,K_2)\) and \(\Theta(z)=z(1+z^2)\).

**Theorem (Prescribed paired product).** The number of cycles admitting two paths at opposite endpoints with a prescribed value of \(J_1J_2\) is \(O_\varepsilon(n^\varepsilon K^2)\). The same bound holds for a prescribed unit residue of \(J_1J_2\) modulo \(n\), or for one prescribed positive rational squareclass of \(\Theta(J_1J_2)\). No additional height bound on either path product is required. In particular all reciprocal-product cycles in the critical rectangle number \(O_\varepsilon(n^{64/153+\varepsilon})\).

Let \(M=B_0 X_1X_2\) and \(N=B_3Y_1Y_2\). Multiplying the two congruences from (1) gives
\[
                             N^2\equiv M^2(J_1J_2)^2\pmod n.     \tag{2}
\]
Both \(M,N\) are strictly between \(0\) and \(n\). Indeed the two weighted contributions at either endpoint are each below \(n\), and their geometric mean is the weighted root product. Fix the two edges at the first endpoint, thereby fixing \(M\). For each unit residue of \(J_1J_2\) the congruence (2) gives a subpower number of unit square-root classes for \(N\) for every composite or prime target. Each has at most one representative in \(0<N<n\).

For each \(N\), factor \(N=B_3Y_1Y_2\) in divisor-many ways. Keep the weighted cubes appropriate to the second endpoint, and test both ordered roots. An endpoint with one of its selected roots determines its outer edge by the exact target equation: divide the positive complementary summand by its side's fixed coefficient and use the unique expression as a square times a squarefree cube. This recovers the outer edges at the second endpoint, so the two middle edges now have prescribed kernel pairs. The selected graph keeps only one representation per edge. This proves subpower recovery from the two-edge seed, including exact closure and the simple-cycle test.

The first two edges cost \(O_\varepsilon(n^\varepsilon K_2^2)\) when the seed is at the left endpoint, by choosing the two distinct right neighbors and using one actual codegree. Reverse both paths to seed at the other end and cost \(O_\varepsilon(n^\varepsilon K_1^2)\), replacing the prescribed product by its reciprocal. The squareclass prescription is invariant under this reversal since \(\Theta(1/z)=\Theta(z)/z^4\). The actual products \(J_1J_2\) have polynomial height by the root caps, so the [uniform squareclass bound](ERDOS_943_PAIRED_PATH_NORM_CLASSES.md#a-uniform-bound-in-one-squareclass) supplies a subpower list of exact products in one class. This proves the theorem. The geometric mean restriction in (2) supplies the fixed-target bound in place of a count of independent residue witnesses.

In the critical rectangle \(K^2\ll n^{64/153}\). One can therefore sum over a collection of at most \(n^{14/153}\) prescribed product residues or squareclasses within the cycle budget \(n^{26/51}\), up to subpower factors. In particular one can remove the cycles where the squarefree representative of \([\Theta(J_1J_2)]\) is at most \(n^{14/153}\). The relative class \([\Theta(J_2)/\Theta(J_1)]\) is a separate statistic.

## Equal products and the endpoint scale

Suppose \(J_1=J_2=p/q\) and \(R_1R_2<\sqrt n\). The two effective slopes \((p/q)X_i/Y_i\) have the same square modulo \(n\). Cancel \(p/q\), and cross multiply the root ratios. The two resulting squares are both less than \(n\), so \(X_1/Y_1=X_2/Y_2\). Use reduced \(\ell/m\) for this common ratio and write \(X_i=\lambda_i\ell,\ Y_i=\lambda_i m\), with \(\lambda_i=\gcd(X_i,Y_i)\). A simple cycle has distinct roots at each common endpoint, so \(\lambda_1\ne\lambda_2\). Order the paths by \(X_1>X_2\) and put
\[
          \eta=\frac{X_2}{X_1}=\frac{Y_2}{Y_1}
                     =\frac{\lambda_2}{\lambda_1}=\frac{r}{s}<1
\]
in lowest terms. Let \(R=\min(R_1,R_2)\). Then \(s\le R\), since \(s\mid\lambda_1\), and \(r\le R\).

For a specified orientation of the first path, use the factorization
\[
 x_1=\alpha u,\quad y_0=\alpha h w,\quad
 y_1=\beta v,\quad x_2=\beta h z,\qquad
 p=u z,\quad q=v w,\quad A_1=h w z .
\]
Here \(x_0=X_1,y_2=Y_1\), and the three edge pairs on this path are \((x_j,y_j)\) for \(0\le j\le2\). Reversing the two paths interchanges \(p,q\) and interchanges \(w,z\). It preserves \(A_1,A_2\) and \(\eta\). Define the common factor for the chosen orientation to be
\[
                               D=\gcd(h w,q,A_2).                 \tag{3}
\]
This can include a prime power shared beyond \(\gcd(w,A_2)\).

**Theorem (Equal products with endpoint scale or common divisor conditions).** For a prescribed \(\eta\), the number of equal-product cycles is \(O_\varepsilon(n^\varepsilon K^2)\), with no upper cutoff on \(J_1\). For \(H,L\ge1\), the sector \(s-r\le H,\ D\ge L\) in either orientation has bound
\[
  O_\varepsilon\!\left(n^\varepsilon\max(K_1,K_2)^2
                     \left(1+\frac{\min(H,R)}{L^2}\right)\right).  \tag{4}
\]
If \(L=1\) one may use \(K^2\) in place of \(\max(K_1,K_2)^2\). The same bound (4) with \(\min(H,R)\) replaced by \(\sqrt R\) holds for any one prescribed squareclass of \(\eta\) and \(D\ge L\), without a reduced-difference cutoff. At critical caps the equal-product cycle budget therefore holds if \(s-r\le n^{14/153}\), or if \(D\ge n^{19/306}\) in either orientation. A fixed squareclass of \(\eta\) permits the smaller divisor threshold \(D\ge n^{5/612}\).

Fix the first two edges of path 1 as a seed. There are \(O_\varepsilon(n^\varepsilon K_1^2)\) such seeds in left-to-right order, and \(O_\varepsilon(n^\varepsilon K_2^2)\) in reverse order. The seed determines \(X_1,\alpha,u\) and \(k=h w=y_0/\alpha\). Choose \(h\mid k\) and \(v\mid y_1\), fixing \(w,q,\beta\) at divisor cost. The remaining variable on path 1 is \(z\). Equation (1) and equality of the endpoint ratios give
\[
 \begin{aligned}
    A_2^2&=\eta^2 A_1^2+(1-\eta^2)(p^2+q^2),\
    (s A_2)^2-\bigl((s^2-r^2)u^2+r^2 k^2\bigr)z^2
                                &=(s^2-r^2)q^2.                  \tag{5}
 \end{aligned}
\]
For fixed \(r,s\) the coefficient of \(z^2\) and the right-hand side are positive, nonzero and polynomially bounded in \(n\). This indefinite binary norm equation has subpower many solutions in the actual box. If its coefficient is a square, factor the difference of squares. Otherwise embed in the maximal real quadratic order by \(s A_2+z\sqrt{(s^2-r^2)u^2+r^2k^2}\). There are at most divisor-many integral ideals of the prescribed positive norm, as in the ideal count for the [fixed-kernel norm bound](https://www.ulam.ai/research/erdos943.pdf). For one ideal, the quotient of two distinct positive such generators, in the order making it greater than one, is a norm-one unit. Its integral trace is at least three, so the quotient is at least \((3+\sqrt5)/2\). The generators in the box lie between one and a fixed power of \(n\), giving only \(O(\log n)\) choices for each ideal.

A surviving \(z\) gives \(x_2=\beta h z\), hence the last actual edge of path 1. The two endpoints of path 2 and its roots \(X_2=\eta X_1,\ Y_2=\eta Y_1\) are now prescribed. They recover its outer edges, and its middle edge is fixed by the two internal kernel vertices. This proves the stated count for fixed \(\eta\), by choosing the orientation with smaller kernel cap.

For (4), write \(C=(B_0p^2\ell^2+B_3q^2m^2)/n\). The primitive endpoint norm has \(A_i^2+C\lambda_i^2=p^2+q^2\) with this same positive integer \(C\) on both paths. For any prime at which \(D\) is nonunit, \(A_i\) are nonunits and \(p^2+q^2\) is a unit since \(D\mid q\). Thus \(C,\lambda_i\) are units there. Subtract the norm equations and cancel \(C\) and the squared gcd of \(\lambda_1,\lambda_2\) to obtain
\[
                              D^2\mid s^2-r^2.                   \tag{6}
\]
For fixed seed, \(h,v,w,q\) have already been chosen, and \(D\) can be chosen among divisors of \(\gcd(hw,q)\). Choose \(s\mid X_1\) at divisor cost. Equation (6) restricts \(r\) to a subpower number of unit square-root classes modulo \(D^2\). The interval defined by \(0<r<s\), \(s-r\le H\) has length at most \(\min(H,R)\). For one prescribed squareclass of \(\eta\), write \(r=d j^2\) with fixed squarefree \(d\) determined by that class and \(s\). The congruence now restricts \(j\) to a subpower number of unit fourth-root classes and \(j\le\sqrt R\). This proves all the counts using (5) for each chosen ratio. If \(L=1\) the bound on \(D\) imposes no condition, so either orientation can be used. The critical thresholds follow from
\[
      \frac{26}{51}-\frac{64}{153}=\frac{14}{153},\qquad
      \frac{11}{51}-2\frac{19}{306}=\frac{14}{153},\qquad
      \frac{11}{102}-2\frac{5}{612}=\frac{14}{153}.
\]

## Equal products of prime-power norm

**Proposition (Prime-power norm restriction).** In the equal-product setting, suppose \(p^2+q^2=P^e\) is a prime power. Then a simple cycle requires \(P\) odd and \(C=P^{2a}\) for some \(1\le a<e/2\), with \(C\) as in the proof of (6). The two primitive pairs \((A_i/P^a,\lambda_i)\) representing \(P^{e-2a}\) as a sum of two squares must be coordinate swaps.

Indeed \(\gcd(A_i,\lambda_i)=1\), because \(A_i\mid y_{0,i}x_{2,i}\) and \(\lambda_i\) divides both outer endpoint roots. Thus \(P\nmid\lambda_i\). If \(v_P(C)\ge e\), then \(P^{\lceil e/2\rceil}\mid A_i\), contrary to \(0<A_i<P^{e/2}\). Otherwise \(v_P(C)=2a<e\) and \(v_P(A_i)=a\) by comparison of valuations. Divide the norms by \(P^{2a}\). The residual coefficient \(C'=C/P^{2a}\) is prime to \(P\). An even prime-power norm \(p^2+q^2\) with coprime \(p,q>0\) can only be \(2\), in which case \(J=1\) yields \(p=q=A_1=1\) and \(t=1\), impossible for a simple path. So \(P\) is odd.

Modulo \(P^{e-2a}\) the two normalized ratios \((A_i/P^a)/\lambda_i\) are unit square roots of \(-C'\) and agree up to sign. Thus the corresponding difference or sum of their cross products is divisible by \(P^{e-2a}\). Its absolute value is at most \(P^{e-2a}/\sqrt{C'}\) by Cauchy–Schwarz. If \(C'>1\), the sum cannot vanish, so the difference vanishes and the two positive norm pairs coincide. This contradicts the distinct endpoint scales. If \(C'=1\), the difference is strictly less than \(P^{e-2a}\) since both coordinates are positive, and equality in the bound for the sum occurs only at a coordinate swap. Finally if \(a=0\), compare either pair with \((p,q)\), which is already a primitive representation of the same prime power with \(C=1\). The same root and cross-product argument forces \(A_i=p\) or \(q\), which means one of the two consecutive internal ratios is one. Two selected edges at a shared vertex cannot have the same root and different neighbors. This proves the assertion.

## Example and remaining cycle count

At \(n=26497\) take the [two paths with internal products \(14/45\) and \(45/14\)](ERDOS_943_PAIRED_THREE_EDGE_PRODUCTS.md#exact-examples-and-remaining-count). Both endpoint kernels are \(1\). The two outer left endpoint roots are \(82,73\), and the two right endpoint roots are \(73,82\). Thus \(M=N=5986\). From the left seed and \(J_1J_2=1\), (2) allows \(N=5986\) or \(20511\). The former recovers the actual factorizations with right roots \(73,82\); the test counts only completions that pass all six edge equations.

**Corollary (Critical residual after paired products and equal-product scale tests).** From the [effective slope and elliptic dependence residual](ERDOS_943_PAIRED_PATH_ELLIPTIC_DEPENDENCE.md#example-and-residual), delete cycles having a paired product \(J_1J_2\) whose combined squareclass representative is at most \(n^{14/153}\), including every reciprocal pair. Also delete equal-product cycles meeting either of the two unconditional critical thresholds in (4). The critical-budget count on the complement is equivalent to (H) with the same coefficient family.

The two theorems count every removed sector at the pointwise critical budget. The reverse implication follows by subset inclusion. Inside the equal-class sector of the new residual, the two products are equal with larger reduced endpoint differences and smaller common divisors, or they give independent points modulo torsion. The frequency of the unrestricted remaining product classes also remains open.

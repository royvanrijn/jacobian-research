# Matching-ratio counts in the critical six-cycle rectangle

Use the coefficient-shaped graph and target family of [the critical rectangular criterion](ERDOS_943_CRITICAL_RECTANGULAR_CYCLE_CRITERION.md). Fix \(n\) and the two coefficients \(\kappa_1,\kappa_2\). The selected graph has kernel caps \(K_1,K_2\) from (2s), so
\[
   \kappa_i K_i^3\le8n^{32/51},\qquad
   (K_1K_2)^3\le64n^{64/51}.
\]
Every simple six-cycle has two alternating perfect matchings. On either one let \(L\) and \(M\) be the products of the three selected roots on the left and right, respectively, and let \(B,D\) be the products of the three corresponding kernels. Both matchings have the same \(B,D\). The two ratios \(L/M\) agree exactly in the singular case. The estimates below apply whenever **either** matching has the stated ratio property.

## Counts by height and squareclass

**Theorem (Matching height and squareclass).** The number of actual simple six-cycles having one matching ratio of reduced height at most \(Y\ge1\) is
\[
               \ll_\varepsilon n^\varepsilon
                      \bigl(n^{5/17}Y^{8/3}+n^{13/51}Y^2\bigr).     \tag{1}
\]
If additionally the squarefree part of the reduced numerator times denominator is at most \(J\ge1\), the same count is at most
\[
               \ll_\varepsilon n^\varepsilon J^{1/2}
                        \bigl(n^{5/17}Y^{5/3}+n^{13/51}Y\bigr).    \tag{2}
\]
For one specified squareclass, (2) holds with the factor \(J^{1/2}\) replaced by \(1\). Thus height \(Y\le n^{11/136}\) by itself, or \(Y\le n^{11/85}\) in a specified squareclass, gives the critical budget \(O_\varepsilon(n^{26/51+\varepsilon})\). These statements include nonsingular cycles.

Write the reduced matching ratio as \(p/q\). Multiplying the three actual edge congruences on that matching and canceling the common root divisor gives
\[
                     n\mid \kappa_1^3 p^2 B^3+\kappa_2^3 q^2 D^3. \tag{3}
\]
All factors are units modulo \(n\) by coprimality of the sums. Write \(B=g b,D=g d\) with \((b,d)=1\). Positivity in (3) gives
\[
              g\le
             \frac{\kappa_1 p^{2/3} K_1^3+\kappa_2 q^{2/3} K_2^3}
                  {n^{1/3}}.
\]
For fixed \(p,q,g\), \(b/d\) lies in a subpower number of unit cubic residue classes. In one such class the reduced fractions \(b\le K_1^3/g,\ d\le K_2^3/g\) number \(O(1+(K_1 K_2)^3/(g^2n))\) by the [rectangular fraction-packing argument](ERDOS_943_CRITICAL_RECTANGULAR_CYCLE_CRITERION.md#an-unweighted-ratio-with-no-coefficient-loss). Summing \(g\) and using the shaped caps bounds the number of kernel products per ratio by
\[
                    \ll_\varepsilon n^\varepsilon
                                   \bigl(n^{5/17}Y^{2/3}+n^{13/51}\bigr).
\]
Each product pair has subpower many ordered kernel factorizations. A matched labeling fixes the cyclic kernel tuple up to finitely many orderings, and between any two kernels the graph has at most one selected edge. This accounts for actual cycles and gives (1) using \(Y^2\) centers.

For a specified squarefree value \(w\) of \(\operatorname{sf}(pq)\), write \(p=d_1 r^2,\ q=d_2 s^2\) with \(d_1d_2=w\). The number of such centers is at most \(\tau(w)Y/\sqrt w\), since \(r\le\sqrt{Y/d_1},\ s\le\sqrt{Y/d_2}\). Summation over \(w\le J\) gives \(O_\varepsilon(n^\varepsilon YJ^{1/2})\) relevant centers, using the polynomial height cap from the selected roots. This proves (2). There is no use of a coherent labeling in (3).

## Cube ratios and the target factorization

Write \(n=3^e n_0\), with \(3\nmid n_0\). Define \(n_-\) as the product of the full prime powers of \(n_0\) at primes \(p\equiv2\bmod3\), including \(2\). Let \(\mathcal D(n)\) be the unitary divisors of \(n_0\) divisible by \(n_-\), where unitary means \(\gcd(d,n_0/d)=1\). Put \(L_3=Q_3=1\) if \(e=0\), and
\[
                    L_3=3^{\max(1,e-1)},\qquad Q_3=3
\]
if \(e\ge1\). Set
\[
             C(n)=\min_{u\in\mathcal D(n)}
                       \max\!\left(uL_3,\sqrt{(n_0/u)Q_3}\right).
                                                                    \tag{4}
\]
In particular \(C(n)\ge n^{1/3}\).

**Theorem (Matching cube ratios with target factors).** For \(R\ge1\), the number of actual cycles with a matching ratio of the form \((r/s)^3\) with \(r,s\) positive coprime and \(\max(r,s)\le R\) is
\[
             \ll_\varepsilon n^\varepsilon
                       \left(\frac{n^{32/51}R^4}{C(n)}
                                           +n^{13/51}R^2\right).  \tag{5}
\]
There are no such cycles if \(16n^{32/51}R^2<C(n)\). In particular:

- Uniformly in \(n\), the cycle budget holds at \(R\le n^{11/204}\).
- For prime \(n\equiv1\bmod3\), it holds at \(R\le n^{13/136}\).
- For prime \(n\equiv2\bmod3\), the sector is empty for \(R<n^{19/102}/4\).

For such a matching, (3) gives
\[
                 n\mid X^3+Y^3,\qquad
       X=\kappa_1 r^2 b,\qquad Y=\kappa_2 s^2 d
\]
after canceling the unit \(g^3\). Both \(X,Y\) are units modulo \(n\), and \(X+Y\le(\kappa_1r^2K_1^3+\kappa_2s^2K_2^3)/g\). The factorization is
\[
                        X^3+Y^3=(X+Y)(X^2-XY+Y^2).
\]
A prime power of \(n_0\) goes entirely into one of the two factors because simultaneous vanishing at its prime would force \(3Y^2\equiv0\), impossible. If the prime is \(2\bmod3\), cubing is injective in its residue field, so the prime power goes into the linear factor. At \(3\mid n\), \(X+Y\equiv0\bmod3\). Then \(X^2-XY+Y^2=(X+Y)^2-3XY\) has \(3\)-valuation exactly one, so \(L_3\mid X+Y\) and \(Q_3\mid X^2-XY+Y^2\). Thus there is \(u\in\mathcal D(n)\) whose \(uL_3\) divides the linear factor and whose \((n_0/u)Q_3\) divides the quadratic factor. Both are positive and the quadratic factor is at most \((X+Y)^2\). It follows that \(X+Y\ge C(n)\), and hence
\[
                            g\le\frac{16n^{32/51}R^2}{C(n)}.
\]
For \(u\in\mathcal D(n)\), the product of its two required moduli is \(n_0L_3Q_3\ge n\), proving \(C(n)\ge n^{1/3}\). At an odd prime \(n\equiv1\bmod3\) the two options are \(u=1,n\), so \(C(n)=\sqrt n\). At a prime \(2\bmod3\) the only option is \(u=n\) and \(C(n)=n\).

For each \(r,s,g\) the same residue packing as for (1) costs a subpower multiple of \(1+(K_1K_2)^3/(g^2n)\). Summing \(g\) and at most \(R^2\) choices of \(r,s\), then factoring the products as before, proves (5). The argument is a target-sensitive extension of the [singular matching product count](ERDOS_943_COHERENT_PRODUCT_RATIO.md#4-counting-the-positive-product-pairs). It retains the linear and quadratic factors even when different prime powers of a composite target use different factors.

## Full-cycle powers and squareclasses at targets with large sign divisors

Use the ratio of the two matching ratios. With the roots as in [the exact cycle equations](ERDOS_943_SIX_CYCLES.md#exact-equations-for-cycle-labels), this is \(U/V\) where \(U=\prod_i a_i\widetilde c_i\) and \(V=\prod_i\widetilde a_i c_i\). Write \(U/V=P/Q\) in lowest terms and \(H=\max(P,Q)\). In this section the cycle is nonsingular, so \(U\ne V\). Set \(X_0=(2n^{11/51})^6\); both \(U,V\) are at most \(X_0\).

Write \(n=2^e n_0\) with \(n_0\) odd. If \(e=0\) set \(a_2=b_2=1\); otherwise set \(a_2=2^{\max(1,e-1)}\), \(b_2=2\). Define
\[
   \Sigma(n)=\min_{\substack{d\mid n_0\(d,n_0/d)=1}}
                      \max(a_2d,b_2 n_0/d).                         \tag{6}
\]
This is at least \(\sqrt n\). For any odd prime-power factor \(p^k\Vert n\), it is at least \(p^k\). For a prime power target it is at least \(n/2\).

**Theorem (Full-cycle powers at large sign divisors).** The nonsingular cycles with \(U/V\) a rational proper power number
\[
                \ll_\varepsilon n^\varepsilon
                     \left(\frac{X_0}{\Sigma(n)}+\frac{X_0}{n}\right). \tag{7}
\]
The same bound holds for any one prescribed positive rational squareclass, uniformly in that class. In particular both counts meet the critical budget for \(\Sigma(n)\ge n^{40/51}\), and are \(O_\varepsilon(n^{5/17+\varepsilon})\) for prime-power targets.

The exact edge equations give \(n\mid U^2-V^2\), and all roots are units modulo \(n\), so \(n\mid P^2-Q^2\). Each odd primary component of \(n\) divides just one of \(P-Q,P+Q\). At \(2\mid n\), both factors are even and the smaller of their two valuations is exactly one, since \(P,Q\) are odd. The other factor has valuation at least \(\max(1,e-1)\). Assign \(d\) to the factor with this latter valuation, or use either sign when \(n\) is odd. This gives \(H\ge\Sigma(n)/2\). Both signs are nonzero because the cycle is nonsingular.

For a fixed reduced ratio \(P/Q\), there are at most \(X_0/H\) choices of the common scale \(h\) with \(U=hP,V=hQ\). The ordered factorizations of each product into its six roots cost a subpower factor. With roots prescribed, the nonzero determinant in the exact cycle equations determines the kernel contributions uniquely: the left unknowns are \(\kappa_1 b_i^3\), and the right unknowns then follow from their edge equations as \(\kappa_2 d_i^3\). This is the same elimination as in the unweighted formula. Rejecting nonintegral or non-squarefree kernels only reduces the count.

For one fixed squareclass, allocate its squarefree representative as \(d_1d_2=\operatorname{sf}(PQ)\) and write \(P=d_1r^2,Q=d_2s^2\). There are divisor-many allocations when the class occurs, as \(P,Q\le X_0\). The condition \(n\mid P^2-Q^2\) restricts \(r/s\) to a subpower number of unit fourth-root classes. In a dyadic height interval of size \(H_0\), the primitive fractions \(r/s\) have coordinates \(O(\sqrt{H_0})\). The fraction-packing bound gives \(O_\varepsilon(n^\varepsilon(1+H_0/n))\) centers.

For the proper-power result, write \(P=r^j,Q=s^j\) for some \(j\ge2\). For each bounded \(j\), the same argument uses unit \(2j\)-th root classes, and the primitive coordinates are \(O(H_0^{1/j})\). The bound \(O_\varepsilon(n^\varepsilon(1+H_0/n))\) still applies, by the unit groups at prime powers and CRT. To make the proof uniform in \(j\), fix a cutoff \(J\) sufficiently large in terms of \(\varepsilon\). For \(j>J\), the total number of pairs per exponent is at most \(X_0^{2/J}\), a subpower allowance. There are \(O(\log n)\) exponents because at least one of \(r,s\) is at least two. Combining the two exponent ranges, multiplying by \(X_0/H_0\), and summing \(H_0\gg\Sigma(n)\) proves (7) for both sectors.

This is a coefficient-shaped bound using the target's possible sign allocations. The [balanced fixed-ratio count](ERDOS_943_SIX_CYCLES.md#fixed-reduced-ratios-and-proper-powers-for-every-target) also treats balanced cycles when no large sign divisor is available.

## Removing these sets from the cycle problem

**Corollary (Matching-ratio residual equivalence).** From the residual of the [joint continuation criterion](ERDOS_943_JOINT_RECTANGULAR_CONTINUATION.md#the-remaining-cycle-criterion), one may remove all cycles covered by the critical-budget thresholds in (1), (2) or (5). For (2) use any fixed squareclass, or any ranges \(J,Y\) for which its displayed upper bound is within the budget. The pointwise \(O_\varepsilon(n^{26/51+\varepsilon})\) bound for the resulting residual is equivalent to (H) for the same coefficient family, as long as only subpower many such sector tests are imposed. In particular one can use (1), the general cube threshold and the sharper prime cases of (5) with no squareclass choices. For targets with \(\Sigma(n)\ge n^{40/51}\), one may also remove the proper-power full-cycle ratios and any fixed squareclass of full-cycle ratios by (7). The full residual estimate remains open.

Subset inclusion gives one implication. For the reverse implication the uniform bounds count all removed cycles, and the previous exact equivalence finishes the argument.

# Pointwise bounds for sums of two \(k\)-full integers

Let \(r_k(N)\) count the ordered representations of \(N\) as the sum of two positive \(k\)-full integers. An integer is \(k\)-full if every prime that divides it has exponent at least \(k\). Blomer and Schöbel proved \(r_k(N)\ll_{k,\varepsilon}N^{2/(2k+1)+\varepsilon}\) in Corollary 2 of their 2013 paper. This includes the exponent \(2/5\) for square-full numbers, which also appears in Chojecki's March 2026 note. The bound for \(k=2\) obtained by the [two-graph and mixed-cubic argument](ERDOS_943_DUAL_KERNEL_GRAPHS.md) is \(N^{3/10+\varepsilon}\). [Blomer–Schöbel, Corollary 2](https://www.researchgate.net/publication/268024970_Twins_of_powerful_numbers) · [Chojecki, Theorem 1.2](https://www.ulam.ai/research/erdos943.pdf)

**Theorem.** For every fixed integer \(k\ge3\), every \(N\ge1\), and every \(\varepsilon>0\),
\[
                           r_k(N)\ll_{k,\varepsilon}N^{1/(k+1)+\varepsilon}. \tag{1}
\]
This improves the bound in Corollary 2 for every \(k\ge3\). Priority for (1) remains unverified.

## Canonical factors and fixed-coordinate multiplicity

The canonical factorization of a \(k\)-full integer is \(m=\prod_{j=k}^{2k-1}a_j^j\), with the \(a_j\) for \(j>k\) squarefree and pairwise coprime. The factor \(a_k\) may share primes with those factors. For an exponent \(qk+r\ge k\), \(0\le r<k\), a prime belongs only to \(a_k\) with exponent \(q\) if \(r=0\). Otherwise it belongs to \(a_{k+r}\) with exponent one and to \(a_k\) with exponent \(q-1\). This shows both existence and uniqueness. Replacing any one coordinate by a divisor of itself preserves the canonical conditions.

We first consider coprime positive summands in
\[
                \kappa_1\prod_{j=k}^{2k-1}a_j^j+
                \kappa_2\prod_{j=k}^{2k-1}c_j^j=n,                           \tag{2}
\]
where the products are in canonical form and the positive coefficients are fixed. All the factors in both terms are units modulo \(n\), and \((a_j,c_i)=1\) for every \(j,i\). The estimates below are uniform in the positive coefficients.

**Lemma.** For fixed \(e\ge3\), the equation \(P a^e+Q c^e=n\) with positive integers \(P,Q,n\) has \(O_{e,\varepsilon}(n^\varepsilon)\) positive integer solutions, uniformly in \(P,Q\).

For the irreducible binary form, apply Evertse's Theorem 1 with \(S\) consisting of the real place and the primes dividing \(n\). It bounds the number of \(S\)-unit cosets by \((5\cdot10^6e)^{1+\omega(n)}\). At the fixed positive value each coset contributes at most one positive pair by homogeneity. [Evertse, Theorem 1](https://pub.math.leidenuniv.nl/~evertsejh/95-Thue-Mahler.pdf)

For the reducible case, the homogeneous form \(P X^e+QY^e\) has \(e\) distinct linear factors over \(\mathbb C\). Here \(P,Q\) are nonzero and \(P t^e+Q\) has simple roots, and \(Y\) does not divide the homogeneous form. Use Gauss's lemma to write it as the product of two coprime homogeneous integer forms \(G,H\) of positive degrees \(d,f\), placing the content in one factor. Fix the two nonzero integer values \(u=G(a,c),v=H(a,c)\), with \(uv=n\). There are at most \(2\tau(n)\) choices. Then \(t=a/c\) satisfies
\[
                         v^d G(t,1)^f-u^f H(t,1)^d=0 .
\]
The polynomial is nonzero, since \(G(t,1),H(t,1)\) are nonconstant and coprime. Its degree is at most \(df\). Each \(t\) determines at most one positive \(c\) by homogeneity. This gives \(O_e(\tau(n))\) solutions in the reducible case.

## Codegrees for any omitted exponent

Fix an omitted exponent \(e\in\{k,\ldots,2k-1\}\). Vertices on each side are the tuples of coordinates other than the one of exponent \(e\). An edge means that a representation (2) exists for those tuples. Fixing both vertices leaves \(P a^e+Qc^e=n\), so each edge carries \(O_{k,\varepsilon}(n^\varepsilon)\) representations by the lemma.

Suppose two distinct left vertices, with products \(P_1,P_2\) including their common coefficient \(\kappa_1\), share a right vertex with product \(Q\) including \(\kappa_2\). The omitted coordinates \(a_i,c_i\) then satisfy
\[
                  P_i a_i^e+Q c_i^e=n,\qquad i=1,2.
\]
Put \(a_i=g\alpha_i,\ c_i=h\gamma_i\), using the gcd of each pair. Coprimality ensures that \(u/v=\alpha_1\gamma_2/(\alpha_2\gamma_1)\) is already reduced. Eliminating \(Q\) gives
\[
       (u/v)^e\equiv P_2/P_1\pmod n,\qquad
                 g^e(P_1u^e-P_2v^e)=n(\gamma_2^e-\gamma_1^e).               \tag{3}
\]
The unit residue on the right has \(O_k(e^{\omega(n)})\) possible roots, by the local unit groups and CRT. In a dyadic box where the two omitted coordinate scales have product at most \(\sqrt n\), each residue has \(O(1)\) reduced fractions \(u/v\). Indeed \(u,v\le4\sqrt n\), and the vectors for distinct fractions in the same residue class have nonzero determinant divisible by \(n\). Their angles in the first quadrant are separated by \(\gg n/(\max(u,v,u',v'))^2\). This is the angular reconstruction bound, including the constant-width boundary of the height range.

For any one reduced fraction there are divisor-many ways to choose its four factors. If the coefficient of \(g^e\) in (3) vanishes in an actual collision, then \(\gamma_1=\gamma_2=1\), so \(c_1=c_2\). The complete first summands are equal, contradicting distinct vertices by uniqueness of the canonical decomposition after removing the fixed coefficient. Thus (3) determines \(g\). Then \(Qh^e=(n-P_1(g\alpha_1)^e)/\gamma_1^e\) is determined. Divide by \(\kappa_2\) and canonically decompose the result. This recovers the right vertex and \(h\) because replacing the omitted coordinate by its divisor \(h\) preserves the canonical conditions. Hence the maximum codegree on each side is \(O_{k,\varepsilon}(n^\varepsilon)\) in this height range.

If the exponents of the two vertex counts in base \(n\) are \(p,p'\), the number of representations divided by a subpower factor is bounded by
\[
         n^p+n^{p'}+\min(n^{p+p'/2},n^{p'+p/2}).                            \tag{4}
\]
This follows by applying Cauchy–Schwarz to the degrees and bounding their pair count by the corresponding maximum codegree on the other side. Without a height restriction the fixed-pair lemma also gives \(O_{k,\varepsilon}(n^{p+p'+\varepsilon})\) by counting all possible vertex pairs.

## The exponents \(k\) and \(k+1\) suffice

For one side of the box let \(x,y,z\) denote the logarithmic scales of its coordinates of exponents \(k,k+1\) and the sum of the remaining logarithmic scales, respectively. Use \(x',y',z'\) for the other side. With \(t=1/(k+1)\) the size constraints give
\[
 kx+(k+1)y+(k+2)z\le1
\]
and its primed twin. Fixing a complete side counts at most \(n^{x+y+z}\) representations, since its complementary summand and all its coordinates are then fixed. We may therefore assume \(x+y+z\ge t\) and \(x'+y'+z'\ge t\).

The displayed size constraint implies \(x\ge z\). Indeed its left side is \((k+1)(x+y+z)-(x-z)\). Set \(p=y+z\) and \(q=x+z\), with \(p',q'\) defined symmetrically. We have \(p,p'\le t\), and
\[
                               q\le 2(t-p),\qquad q'\le2(t-p').             \tag{5}
\]
To see the first inequality multiply the size constraint by \(2/(k+1)\). The resulting left side minus \(q+2p\) is \((k-1)(x-z)/(k+1)\ge0\).

For the graph omitting \(k+1\), the omitted-coordinate product has logarithmic scale \(y+y'\le 2/(k+1)\le1/2\), as required for (4). For the graph omitting \(k\), if \(x+x'>1/2\) then
\[
                         p+p'\le \frac{2-k(x+x')}{k+1}<t .
\]
The fixed-pair bound covers this case. Otherwise (4) applies to both graphs. The graph omitting \(k\) gives the desired bound unless both \(p+p'/2>t\) and \(p'+p/2>t\), since its linear terms have exponents at most \(t\). In the remaining case \(p,p'>t/2\), so \(q,q'<t\). Furthermore (5) gives \(q+q'/2\le3t-2p-p'<t\) and the symmetric inequality. The graph omitting \(k+1\) gives the desired bound. This proves the uniform estimate in each box, and there are \(O_k((\log(2n))^{2k})\) boxes.

For unrestricted positive \(k\)-full summands \(m_1,m_2\), let \(h=(m_1,m_2)\). Any \(p\mid h\) has exponent at least \(k\) in \(h\). If a prime occurs to positive exponent below \(k\) in \(m_i/h\), it divides \(h\), and its full power in that quotient divides \(h\) too. Extract the product \(\kappa_i\) of these prime powers. Then \(\kappa_i\mid h\), and \((m_i/h)/\kappa_i\) is \(k\)-full. Apply (2) with \(n=N/h\). There are at most \(\tau(N)^3\) choices of \(h,\kappa_1,\kappa_2\); the coefficient-uniform bound and the divisor estimate prove (1). For \(k>2\) the coefficients may contain prime powers up to exponent \(k-1\).

## Square-full kernel pairs at the balanced scale

For \(N\asymp X\) and \(T=X^{1/5}\), let \(K_N(T)\) count ordered pairs of distinct squarefree kernels \(b_1,b_2\in[T,2T)\), coprime to \(N\), for which there is a reduced unit fraction \(u/v\) with \(u,v\le4T^2\) and \(b_1^3u^2\equiv b_2^3v^2\bmod N\).

**Proposition.** On this scale,
\[
                    \frac1X\sum_{X\le N<2X}K_N(T)\ll_\varepsilon T^{1+\varepsilon}.
                                                                            \tag{6}
\]
Also there are infinitely many \(N\) with \(K_N(N^{1/5})\gg N^{1/5}\), even using coprime pairs of kernels.

For the mean, count the \(O(T^6)\) candidate tuples \(b_1,b_2,u,v\). Their differences \(b_1^3u^2-b_2^3v^2\) are nonzero by squarefreeness and \(b_1\ne b_2\), and have size \(O(T^7)\). The number of possible \(N\) for each is bounded by the divisor count \(O_\varepsilon(T^\varepsilon)\). Division by \(X=T^5\) proves (6). In particular the number of \(N\in[X,2X)\) with \(K_N(T)>T^{2-\eta}\) is \(O_\varepsilon(X^{(4+\eta)/5+\varepsilon})\).

For the lower bound use the polynomial identity
\[
 \begin{aligned}
 &x^3(-x^2+14xy+35y^2)^2-y^3(35x^2+14xy-y^2)^2\\
 &\hspace{35mm}=(x-y)^5(x^2-23xy+y^2).
 \end{aligned}                                                               \tag{7}
\]
To define an infinite sequence, take \(r=21000,\ q=r^2-23r+1\), \(L=2^m\) with \(m\ge5\), and \(N=qL^5/256\). Then \(T=N^{1/5}>16L\). Take \(y\in[T,33T/32)\) odd with \((r-1)y\equiv L\bmod q\), and put \(x=y+L\). The congruence is possible since \((r-1,q)=1\); both \(x\equiv r y\bmod q\) and \(y\) are units, since \((r,q)=1\) and \(q\) is odd. Both kernels are in \([T,35T/32)\), and \((x,y)=1\) since \(L\) is a power of two and \(y\) is odd.

There are \(\gg T\) choices for which both kernels are squarefree. In the one progression modulo \(2q\), no prime of \(2q\) divides either kernel. At any other prime \(p\), the proportion lost to \(p^2\) dividing one of them is at most \(2/p^2\), up to an additive count of \(2\). Summing the errors over \(p\le\sqrt{35T/32}\) gives \(O(\sqrt T)\), and \(2\sum_{p\ge3}p^{-2}<1\). The fixed progression has \(\asymp T\) points, so a positive proportion survive.

Let \(U=-x^2+14xy+35y^2\) and \(V=35x^2+14xy-y^2\). They are positive here. Since \(L\) is divisible by 32 and \(y\) is odd, both are \(16\bmod32\). Also neither has a prime factor in common with \(q\). If a prime of \(q\) divided either, it would divide both by (7) and the unit property of \(x,y\). Then it would divide \(U-V=36(y^2-x^2)\), so \(x/y\equiv\pm1\) at this prime. The values \(U/y^2\) would be \(48\) or \(20\), impossible since \(q\) is coprime to \(2\cdot3\cdot5\cdot7\). Now divide \(U,V\) by their gcd \(g\), writing the reduced fraction \(u/v\). Its two coordinates are odd, coprime to \(q\), and at most \(49(35T/32)^2/16<4T^2\). The factor \(x^2-23xy+y^2\) is divisible by \(q\), by \(x\equiv ry\bmod q\). Dividing (7) by \(g^2\), whose 2-part is \(256\) and whose odd part is coprime to \(q\), preserves divisibility by \(qL^5/256\). Thus every surviving pair \((x,y)\) contributes to \(K_N(T)\). This proves the lower bound. In particular a uniform upper bound \(O_\varepsilon(T^{2-\eta+\varepsilon})\) cannot hold for \(\eta>1\).

Pointwise, pairs with \(\gcd(b_1,b_2)\ge T^\rho\) contribute \(O(T^{2-\rho})\) for fixed \(0<\rho<1\), by summing \(O((T/g)^2)\) over \(g\ge T^\rho\). Moreover every witness of height \(H\) requires \(H\gg T g^{3/2}\), since cancelling \(g^3\) modulo \(N\) leaves a nonzero difference of size \(O((T/g)^3 H^2)\). This locates the remaining pointwise problem at small kernel gcd, with rational heights between order \(T\) and order \(T^2\).

**Open question.** Prove \(K_N(T)\ll_\varepsilon T^{2-\eta+\varepsilon}\) for an absolute \(0<\eta\le1\) at \(T\asymp N^{1/5}\). The analogous bounds with fixed interval and height constants would give a balanced-box improvement in the square-full representation problem. Indeed, the [kernel graph common-neighbor argument](ERDOS_943_KERNEL_GRAPH.md) bounds the total neighbor pair count by \(O_\varepsilon(N^\varepsilon K_N(T))\), since each compatible pair contributes at most a subpower number of common neighbors. Cauchy–Schwarz then gives \(O_\varepsilon(N^\varepsilon(T+\sqrt{T K_N(T)}))\) representations in the balanced box, with exponent \(3/10-\eta/10+\varepsilon\) under the proposed estimate.

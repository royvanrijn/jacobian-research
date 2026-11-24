# A power saving for the short-box congruence in Erdős #943

Let \(r(N)\) count ordered representations of \(N\) as the sum of two positive square-full integers. Each square-full integer has a unique expression \(a^2b^3\) with \(b\) squarefree. For a representation
\[
 a^2b^3+c^2d^3=N
\]
the triple \((b,c,d)\) determines at most one positive \(a\), and it satisfies \(c^2d^3\equiv N\bmod b^3\). The congruence can be counted even after dropping the squarefree and positivity-of-the-quotient restrictions. The pointwise \(N^{2/5+\varepsilon}\) benchmark and the dyadic cube-modulus projection are discussed in Chojecki's note. [Erdős #943](https://www.erdosproblems.com/943) · [Chojecki, Theorem 1.2 and §8](https://www.ulam.ai/research/erdos943.pdf)

## Completion and lifting

Write \(M_\ell(q;C,D)\) for the number of solutions to \(c^2 d^3\equiv\ell\bmod q\) in \(C\le c<2C,\ D\le d<2D\), with \((\ell,q)=1\). In the unit case the map
\[
 t=cd,\qquad (c,d)=(\ell^{-1}t^3,\ell t^{-2})\pmod q                \tag{1}
\]
is a bijection from the units modulo \(q\) to all modular solution pairs. The completion bound used below is Chan's Theorem 4. The square-root denominator in its interval-length error applies to the whole numerator \(K+L\). [Chan, Theorem 4 and Lemma 4](https://arxiv.org/pdf/1407.0054)

For an independent check using (1), write \(e_q(x)=\exp(2\pi i x/q)\) and consider
\[
 S_q(A,T)=\sum_{t\bmod q}^{*} e_q(A t^3+Tt^{-2}).
\]
Here and below the star restricts the sum to units. Chinese remaindering reduces the bound to prime powers, with the local frequencies multiplied by units. The prime-power calculation gives
\[
 |S_q(A,T)|\ll_\varepsilon q^{1/2+\varepsilon}(q,A,T)^{1/2}.       \tag{2}
\]

To see this, first divide out any common power of \(p\) from \(A,T,p^m\). This shortens the modulus and introduces exactly the matching number of lifts; if the phase is constant modulo \(p^m\), the bound is immediate. For a primitive frequency pair, put \(f(t)=A t^3+Tt^{-2}\). Its stationary congruence is \(3A t^5-2T\equiv0\), since the derivative has the additional unit denominator \(t^3\). The common valuation of \(3A,2T\) is at most one. After its removal, a stationary congruence in the units either has no solutions or specifies \(t^5\) at the shorter precision. For \(p\ne5\) there are at most five simple roots, and for \(p=5\) the kernel of the fifth-power map on units has size at most five. The removal introduces at most a fixed number of extra lifts at \(2,3\).

For \(m=2s\), sum over \(t=u+p^s v\). The \(v\)-sum forces stationarity modulo \(p^s\) and contributes \(p^s\). For \(m=2s+1\) with \(s\ge1\) and \(p>5\), use \(t=u+p^s(v_0+pv_1)\). The \(v_1\)-sum forces the same congruence, and the \(v_0\)-sum is a quadratic Gauss sum of size \(\sqrt p\), since the derivative at a stationary residue is simple. For \(p\le5\), splitting at \(p^{s+1}\) instead loses at most a fixed factor \(\sqrt p\). For exponent one, the rational-function Weil bound applies at \(p>5\), since each nonzero term has a pole of order less than \(p\). At the three remaining primes the trivial bound suffices. This proves (2), with the product of the fixed constants absorbed in \(q^\varepsilon\).

For intervals of integer lengths \(L_1,L_2\), completion yields
\[
 M_\ell(q;C,D)
 =\frac{\varphi(q)L_1L_2}{q^2}
 +O_\varepsilon\!\left(q^{1/2+\varepsilon}
                  +(L_1+L_2)q^{-1/2+\varepsilon}\right).           \tag{3}
\]
Indeed, the nonzero multiples of a divisor \(h\mid q\) contribute total Fourier magnitude \(O((q/h)\log(2q))\) for either interval. Majorize \(\sqrt{(q,k_1,k_2)}\) by the sum of \(\sqrt h\) over common divisors and isolate the two axes separately. The contributions in the weighted Fourier sum are \(O_\varepsilon(q^{2+\varepsilon})\) off the axes and \(O_\varepsilon((L_1+L_2)q^{1+\varepsilon})\) on the axes with the origin deleted. Equation (2) gives (3). For \(q=1\) the count is the product of the interval cardinalities. The use of integer lengths in (3) also covers intervals whose endpoints are not integers.

At the original balanced scale \(q\asymp B^3\), \(C,D\asymp B\), the completion error \(O(q^{1/2+\varepsilon})\) is too large even for the elementary \(B^{1+\varepsilon}\) bound per modulus. For comparison, the direct total-degree modular determinant in a translated square of side \(L\) also misses this scale. With degree cutoff \(J\), the determinant for monomials \(X^iY^j\), \(i+j\le J\), is divisible by \(q^V\), where
\[
 V=\sum_{k\ge1,\ 5k\le J}\binom{J-5k+2}{2}.
\]
Replace each column by \(X^{i-2k}Y^{j-3k}((c_0+X)^2(d_0+Y)^3-N)^k\) with \(k=\min(\lfloor i/2\rfloor,\lfloor j/3\rfloor)\). The leading monomial makes the transformation triangular. The archimedean bound is \(O_J(L^W)\), with \(W=\sum_{i+j\le J}(i+j)\). The ratio \(V/W\) tends to \(1/10\). Thus this determinant criterion takes \(L<q^{1/10}\) up to exponent margins, requiring order \(B^{7/5}\) subboxes at the balanced scale before counting any points inside them.

## Multiplying two solutions

**Theorem.** Uniformly in the unit target and for \(C,D\ge1/2\),
\[
 M_\ell(q;C,D)
 \ll_\varepsilon (2qCD)^\varepsilon
 \left(\frac{CD}{\sqrt q}+q^{1/4}+\frac{C+D}{q^{1/4}}\right).       \tag{4}
\]
In particular, at \(b,c,d\asymp B\), with \(q=b^3\), the total number of unit triples is \(O_\varepsilon(B^{7/4+\varepsilon})\). This includes all the balanced triples when \(N>2B\) is prime.

Multiply two solutions coordinatewise. Their products lie in \([C^2,4C^2)\times[D^2,4D^2)\), and satisfy the same modular curve with target \(\ell^2\). Each resulting pair has at most \(\tau(c_1c_2)\tau(d_1d_2)\ll_\varepsilon (2CD)^\varepsilon\) ordered preimages. Thus the square of the original count is bounded by a subpower factor times the count in the product box. Apply (3) there and take square roots. When \(C,D\ge1/2\), the lengths including rounding are \(O(C^2)\) and \(O(D^2)\), so the estimate holds in those boundary cases too. At \(q=b^3\asymp B^3\) the three terms in (4) are \(O(B^{1/2})\), \(O(B^{3/4})\) and \(O(B^{1/4})\), respectively. There are \(O(B)\) choices of modulus base.

## Composite targets and shared divisors

**Theorem.** For \(N\ge1\) and \(B,C,D\ge1\), the full number of triples in the half-open dyadic box satisfying \(c^2d^3\equiv N\bmod b^3\) is
\[
 \ll_\varepsilon \tau(N)^3(2BCD)^\varepsilon
 \left(\frac{CD}{\sqrt B}+B^{7/4}+(C+D)B^{1/4}\right).              \tag{5}
\]
This treats all moduli \(b^3\), including those from squarefree \(b\) that share primes with \(N\).

For each \(b\), set \(g=(b^3,N)\). At a prime \(p\mid g\), write \(s=v_p(g)\). For a given pair \(c,d\), truncate their exponents at \(\lceil s/2\rceil,\lceil s/3\rceil\), respectively, and let \(u,v\) be the products of those truncated prime powers. They divide \(g\), so for fixed \(g\) there are at most \(\tau(g)^2\) choices. If \(s<3v_p(b)\), the valuation of \(c^2d^3\) equals \(s\), and the two exponents were not truncated below their actual values. If \(s=3v_p(b)\), their weighted sum can exceed \(s\), but truncation preserves the inequality \(2v_p(u)+3v_p(v)\ge s\). Consequently \(g\mid u^2v^3\), \(uv\ge g^{1/3}\), and, for \(q'=b^3/g>1\), the pair \(c/u,d/v\) satisfies a unit congruence
\[
 (c/u)^2(d/v)^3\equiv
       (N/g)\,(u^2v^3/g)^{-1}\pmod {q'}.                           \tag{6}
\]
The coefficient inverted here is a unit at every prime still dividing \(q'\), by the exact valuations at unsaturated primes. For \(q'=1\) the congruence is vacuous and (4) applies by interval counting. A nonempty branch always has \(C/u,D/v>1/2\).

The choice \(g\mid b^3\) forces \(\prod_{p\mid g}p^{\lceil v_p(g)/3\rceil}\mid b\), so at most \(2B g^{-1/3}\) choices of \(b\) occur in this branch. Applying (4) and summing over those \(b\) gives three terms, up to a subpower factor,
\[
 \frac{CD}{B^{1/2}}\frac{g^{1/6}}{uv},\qquad
 B^{7/4}g^{-7/12},\qquad
 B^{1/4}\left(\frac C u+\frac D v\right)g^{-1/12}.
\]
These are at most the corresponding terms of (5). Summing over \(g\mid N,\ u,v\mid g\) proves (5). In particular the improvement at the balanced scale for \(N\asymp B^5\) is \(O_\varepsilon(N^{7/20+\varepsilon})\) without primality or squarefree-modulus restrictions.

There is a structured lower bound for the raw congruence when \(N\) is composite. For \(N=m^5\), take \(c=d=m\). Every \(b\in[m,2m)\) works, because \(c^2d^3=N\) exactly. Taking prime \(m\) leaves \(\gg m\) choices with \((b,N)=1\), including \(\gg m\) squarefree \(b\). Thus a uniform subpolynomial bound for the raw congruence cannot hold, even on its unit locus. Its associated complementary variable is \(a=0\), outside the positive representation problem.

## Summing the unbalanced representation boxes

**Theorem.** For all \(N\ge1\),
\[
                         r(N)\ll_\varepsilon N^{7/19+\varepsilon}.            \tag{7}
\]
The argument applies to all dyadic ranges, including noncoprime summands.

To provide a separate justification of the two box bounds used with (5), fix \(b,d\) first. Multiplying \(b^3a^2+d^3c^2=N\) by \(b^3\) embeds the positive pairs in a binary norm equation of norm \(Nb^3\) in \(\mathbb Q(\sqrt{-b^3d^3})\). There are at most \(\tau(Nb^3)\) ideals of that norm, from the splitting, ramification and inertia possibilities at each prime. Each principal ideal gives at most six generators. Thus the number of pairs for fixed \(b,d\) is \(O_\varepsilon(N^\varepsilon)\). For fixed \(a,c\) use the binary form \(a^2 X^3+c^2 Y^3\). If it is irreducible, Evertse's uniform Thue–Mahler bound at the places above the primes of \(N\), together with the real place, gives \(O_\varepsilon(N^\varepsilon)\) solutions with fixed value \(N\). The theorem counts \(S\)-unit cosets. Two positive pairs in a single coset with the same positive value would differ by a rational scalar whose cube is one, so the pairs agree. If the form is reducible, the coefficient ratio is a rational cube. Taking out the common cube-free core gives a sum \(X_0^3+Y_0^3=N_0\) of two positive cubes. There are at most \(2\tau(N_0)\) ordered solutions by choosing the divisor \(X_0+Y_0\) and recovering the product \(X_0Y_0\) from \(N_0\). [Evertse, Theorem 1](https://pub.math.leidenuniv.nl/~evertsejh/95-Thue-Mahler.pdf)

Write \(A,B,C,D\) for the four dyadic scales of \(a,b,c,d\), relabeling summands so that \(B\le D\). Necessarily \(A^2B^3,C^2D^3\le N\), and the fixed-pair bounds give \(O_\varepsilon(N^\varepsilon\min(AC,BD))\) representations. Put \(U=BD\). If \(U\le N^{7/19}\) or \(U\ge N^{8/19}\), this bound suffices, in the second case because \(AC\le N/U^{3/2}\). In the intervening range apply (5) and the uniqueness of the positive complementary square root \(a\). Its four terms are bounded as follows:
\[
\begin{aligned}
 CD B^{-1/2}&\le N^{1/2}U^{-1/2}\le N^{6/19},\\
 B^{7/4}&\le U^{7/8}\le N^{7/19},\\
 C B^{1/4}&\le N^{1/2}U^{-5/8}\le N^{41/152},\\
 D B^{1/4}&\le N^{1/4}U^{1/4}\le N^{27/76}.
\end{aligned}
\]
The last bound uses \(D\le N^{1/3}\). Both \(41/152\) and \(27/76\) are less than \(7/19\). The divisor bounds and the number \(O((\log(2N))^4)\) of boxes are absorbed by \(\varepsilon\). This proves (7).

**Open question.** Improve (7) toward subpolynomial growth. The box bounds here meet at \(B,D\asymp N^{4/19}\) and \(A,C\asymp N^{7/38}\). In those boxes, the contribution \(B^{7/4}\) arises from the square-root completion error after taking the square root of the product count. The fixed-square-part bound gives \(AC\asymp N^{7/19}\) at the same place. Further improvement from this route needs a sharper count on that product box or information from the complementary square condition discarded by the projection.

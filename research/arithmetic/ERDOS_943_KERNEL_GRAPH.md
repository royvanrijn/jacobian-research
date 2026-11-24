# A one-third pointwise bound for sums of two square-full integers

Let \(r(N)\) count the ordered representations of \(N\) as a sum of two positive square-full integers. Each square-full integer has a unique expression \(a^2b^3\) with \(b\) squarefree.

**Theorem.** For every positive integer \(N\) and every \(\varepsilon>0\),
\[
                              r(N)\ll_\varepsilon N^{1/3+\varepsilon}.       \tag{1}
\]
The kernel-graph argument was proposed by the user. The exponent improves the earlier [pointwise \(5/14\) bound](ERDOS_943_SIXTH_POWER_CARRY.md). Blomer and Schöbel obtained the \(2/5\) benchmark in 2013. Chojecki's March 2026 note gives the same exponent and discusses the pointwise formulation of Erdős #943. [Blomer–Schöbel, Corollary 2](https://www.researchgate.net/publication/268024970_Twins_of_powerful_numbers) Priority for (1) remains unverified. [Chojecki, Theorem 1.2](https://www.ulam.ai/research/erdos943.pdf)

We work first with positive coprime summands in the coefficient equation
\[
             \kappa_1 a^2 b^3+\kappa_2 c^2 d^3=n,                             \tag{2}
\]
with fixed positive coefficients \(\kappa_1,\kappa_2\) and positive squarefree \(b,d\). All bounds are uniform in the coefficients. Each summand is individually coprime to \(n\), since their gcd is 1 and their sum is \(n\). In particular all the variables and both coefficients are units modulo \(n\) on any nonempty box. The identity \((a,c)=1\) follows from coprimality of the summands.

## Common neighbors and reconstruction

**Proposition.** If \(16A^2 C^2<n\), then for fixed distinct squarefree \(b_1,b_2\) there are \(O_\varepsilon(n^\varepsilon)\) pairs of coprime representations (2) sharing \(d\) with \(a_1,a_2\in[A,2A)\) and \(c_1,c_2\in[C,2C)\). The symmetric bound holds for fixed distinct \(d_1,d_2\) and a common \(b\).

Suppose two such representations share \(d\). Multiplying their equations by \(c_2^2\) and \(c_1^2\), respectively, and subtracting gives
\[
 \kappa_1\bigl(b_1^3(a_1c_2)^2-b_2^3(a_2c_1)^2\bigr)
                         =n(c_2^2-c_1^2).                                  \tag{3}
\]
Thus the reduced positive fraction \(u/v=a_1c_2/(a_2c_1)\) represents a unit solution of
\[
                        (u/v)^2\equiv b_2^3/b_1^3\pmod n.                   \tag{4}
\]
For odd primes a unit square has at most two roots modulo each prime power, and at a power of 2 it has at most four. CRT therefore gives at most \(2^{\omega(n)+1}\) roots modulo \(n\). Both \(u\) and \(v\) are less than \(4AC\). If two fractions in this range represent the same residue, their cross difference is divisible by \(n\) and has absolute value less than \(16A^2 C^2\), so it vanishes. Each root of (4) determines at most one reduced fraction in the range.

Write \(a_i=g\alpha_i,\ c_i=h\gamma_i\) with \((\alpha_1,\alpha_2)=(\gamma_1,\gamma_2)=1\). Also \((\alpha_i,\gamma_i)=1\), for each \(i\). Each factor of \(\alpha_1\gamma_2\) belongs to a pair of factors coprime to each of the two factors of \(\alpha_2\gamma_1\). Hence \(u=\alpha_1\gamma_2,\ v=\alpha_2\gamma_1\) with no cancellation. There are at most \(\tau(u)\tau(v)\) choices for the four reduced coordinates. Dividing (3) by \(h^2\) gives
\[
            \kappa_1 g^2(b_1^3 u^2-b_2^3 v^2)
                         =n(\gamma_2^2-\gamma_1^2).                          \tag{5}
\]
The coefficient of \(g^2\) is nonzero. Equality \(b_1^3 u^2=b_2^3 v^2\) would force \(b_1,b_2\) to have the same prime factors by taking the valuations modulo 2. Thus (5) gives at most one positive \(g\). Then
\[
                 \frac{n-\kappa_1 b_1^3 g^2\alpha_1^2}
                      {\kappa_2\gamma_1^2}=h^2d^3
\]
gives at most one pair \(h,d\) by the unique square-times-squarefree-cube decomposition. This reconstructs both representations. The divisor bounds and \(2^{\omega(n)+1}\) give the claimed subpower estimate, with \(u,v<n\) in this range. Swapping the summands proves the symmetric statement.

## The three box ranges

Let \(A,B,C,D\) be the dyadic lower endpoints for \(a,b,c,d\) and \(U=BD\). For a nonempty box,
\[
                         B,D\le n^{1/3},\qquad AC\le n/U^{3/2}.               \tag{6}
\]
The inequalities also hold with arbitrary positive coefficients in (2).

For fixed \(b,d\), write \(L=\kappa_1b^3,\ M=\kappa_2d^3\). The map \((a,c)\mapsto La+c\sqrt{-LM}\) takes the equation to an algebraic-integer norm equation of norm \(Ln\) in the imaginary quadratic field \(\mathbb Q(\sqrt{-LM})\). An imaginary quadratic field has at most \(\tau(Ln)\) integral ideals of norm \(Ln\), by splitting, ramification and inertia prime by prime. Each principal ideal has at most six generators of this norm. Since \(L\le n\), this gives \(O_\varepsilon(n^\varepsilon)\) representations for any fixed kernel pair. In particular the whole box costs \(O_\varepsilon(n^{1/3+\varepsilon})\) when \(U\le4n^{1/3}\).

For \(4n^{1/3}<U<n^{4/9}\), (6) gives \(16A^2 C^2\le 16 n^2/U^3<n\). Form the simple bipartite graph with kernel vertices \(b,d\) and an edge for a represented kernel pair. Both sides have maximum codegree \(O_\varepsilon(n^\varepsilon)\) by the proposition. If the sides have sizes \(m_1,m_2\), the maximum codegree on the first side is \(K\), and the edge count is \(e\), then
\[
 \frac{e^2}{m_2}-e
 \le 2\sum_{d}\binom{\deg d}{2}
 \le 2K\binom{m_1}{2}.
\]
If either side is empty the edge count is zero. Otherwise solving gives \(e\le m_2+\sqrt K\,m_1\sqrt{m_2}\). The symmetric estimate uses the other maximum codegree. Taking the better square-root term, and then applying the uniform fixed-kernel multiplicity bound, gives a box bound
\[
 \ll_\varepsilon n^\varepsilon(B+D+\min(B\sqrt D,D\sqrt B)).
\]
Here \(B+D\le 2n^{1/3}\), and \(\min(B\sqrt D,D\sqrt B)\le U^{3/4}\le n^{1/3}\).

For \(U\ge n^{4/9}\), fix \(a,c\). Then \(AC\le n^{1/3}\), and the equation to solve is
\[
                              P b^3+Q d^3=n,\qquad
                P=\kappa_1 a^2,\quad Q=\kappa_2 c^2.
\]
If the form \(P X^3+QY^3\) is irreducible over \(\mathbb Q\), use Evertse's Theorem 1 with \(S\) consisting of the real place and the primes dividing \(n\). All the integral pairs are \(S\)-integral and the fixed right-hand side is an \(S\)-unit. The bound is \((15\cdot10^6)^{1+\omega(n)}\) on the number of \(S\)-unit cosets, uniformly in the coefficients. At a fixed nonzero right-hand side each coset contains at most one positive pair, since a common scalar must have cube one. This gives \(O_\varepsilon(n^\varepsilon)\) pairs. The Bombieri–Schmidt primitive Thue bound is another applicable choice because \((b,d)=1\) follows from coprimality of the summands. [Evertse, Theorem 1](https://pub.math.leidenuniv.nl/~evertsejh/95-Thue-Mahler.pdf)

If the form is reducible, \(P/Q\) is a rational cube. On dividing by \(\gcd(P,Q)\) the two coprime remaining coefficients are each integer cubes. The equation changes into \(x^3+y^3=n_0\) for a positive integer \(n_0\mid n\). The sum \(s=x+y\) divides \(n_0\), and \(xy=(s^3-n_0)/(3s)\). Each divisor \(s\) thus gives at most two ordered pairs, which again gives \(O_\varepsilon(n^\varepsilon)\) solutions. The \(O((\log(2n))^4)\) dyadic boxes are absorbed into \(\varepsilon\).

## Arbitrary targets and coprime normalization

For square-full summands \(m_1,m_2\) with sum \(N\), take \(h=(m_1,m_2)\mid N\) and \(n=N/h\). Let \(\kappa_i\) be the product of primes having exponent exactly one in \(m_i/h\). Any such prime divides \(h\), since the original \(m_i\) was square-full. Hence \(\kappa_i\) is a squarefree divisor of \(h\); the positive integer \((m_i/h)/\kappa_i\) is square-full. Its square-times-squarefree-cube decomposition is unique. This puts the quotients \(m_i/h\) in (2) with coprime summands, and recovers the original summands given the normalization data. There are at most \(\sum_{h\mid N}\tau(h)^2\le \tau(N)^3\) choices of \(h,\kappa_1,\kappa_2\). The coefficient-uniform bound for (2) and \(n\le N\) prove (1).

The [carry average](ERDOS_943_SIXTH_POWER_CARRY.md) and the [target-independent collision estimate](ERDOS_943_OFF_DIAGONAL_COLLISIONS.md) are separate open estimates. The graph counts common neighbors among actual representations, using the complementary squares in (3) and (5).

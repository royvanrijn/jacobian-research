# A three-tenths pointwise bound for sums of two square-full integers

Let \(r(N)\) be the number of ordered representations of a positive integer \(N\) as a sum of two positive square-full integers. The kernel graph uses the squarefree cube factors as vertices. Its dual uses the square roots as vertices. Combining the two graphs with a mixed-coordinate cubic bound improves the earlier [one-third estimate](ERDOS_943_KERNEL_GRAPH.md). For the problem formulation and the \(2/5\) benchmark, see the March 2026 note. Priority for the bound below remains unverified. [Chojecki, Theorem 1.2](https://www.ulam.ai/research/erdos943.pdf)

**Theorem.** For every positive integer \(N\) and every \(\varepsilon>0\),
\[
                               r(N)\ll_\varepsilon N^{3/10+\varepsilon}.    \tag{1}
\]

We first count coprime positive summands in
\[
                       \kappa_1 a^2b^3+\kappa_2 c^2d^3=n,\qquad
                         b,d\ \text{squarefree},                           \tag{2}
\]
for fixed positive integer coefficients, with a bound uniform in both coefficients. Every factor of each summand is a unit modulo \(n\). Also \((a,c)=(b,d)=1\) within any one representation. Write \(A,B,C,D\) for the lower endpoints of their positive dyadic intervals and \(R\) for the number of representations in the box. A nonempty box has
\[
                               A^2B^3,\ C^2D^3\le n.                       \tag{3}
\]

## Rational reconstruction and the two graphs

**Lemma.** For \(n\ge2\) and \(H\ge1\), any residue class modulo \(n\) contains at most \(O(1+H^2/n)\) reduced positive fractions \(u/v\) with \(u,v\le H\) and \((v,n)=1\).

For two distinct fractions in the class, \(uv'-vu'\) is a nonzero multiple of \(n\). The vectors \((u,v)\), \((u',v')\) lie in the first quadrant and have Euclidean lengths at most \(\sqrt 2 H\), so the sine of their angular separation is at least \(n/(2H^2)\). The separation itself is at least that large. Packing the angles in an interval of length \(\pi/2\) gives the bound. When \(H^2<n\) the cross-difference bound also gives uniqueness.

A unit has at most \(2^{\omega(n)+1}\) square roots and \(3^{\omega(n)}\) cube roots modulo \(n\). The local root bounds are two for an odd prime power in the square case and four for a power of 2. For cubes there are at most three roots at each prime power. At primes other than 3 this follows by simple-root lifting, and for \(3^k\) the kernel of cubing in the cyclic unit group has at most three elements. CRT proves the bounds. In particular these counts grow at most as \(O_\varepsilon(n^\varepsilon)\).

**Proposition.** The two graph estimates for the box in (2) are
\[
\begin{aligned}
 R&\ll_\varepsilon n^\varepsilon
   \left[B+D+(1+AC/\sqrt n)\min(B\sqrt D,D\sqrt B)\right],\\
 R&\ll_\varepsilon n^\varepsilon
   \left[A+C+(1+BD/\sqrt n)\min(A\sqrt C,C\sqrt A)\right].                     \tag{4}
\end{aligned}
\]

For the first estimate, two distinct \(b_1,b_2\) sharing \(d\) satisfy the exact identity
\[
 \kappa_1\bigl(b_1^3(a_1c_2)^2-b_2^3(a_2c_1)^2\bigr)=n(c_2^2-c_1^2).
\]
Thus \(u/v=a_1c_2/(a_2c_1)\) is a modular square root of \(b_2^3/b_1^3\). Reduce the two root pairs by \(a_i=g\alpha_i,\ c_i=h\gamma_i\). Pairwise coprimalities \((\alpha_1,\alpha_2)=(\gamma_1,\gamma_2)=(\alpha_i,\gamma_i)=1\) give \(u=\alpha_1\gamma_2\), \(v=\alpha_2\gamma_1\) as written in lowest terms. Each fraction therefore gives \(O_\varepsilon(n^\varepsilon)\) choices of the four reduced factors, by the divisor bound and \(u,v<4AC\le4n\). The equation
\[
             \kappa_1 g^2 (b_1^3u^2-b_2^3v^2)=n(\gamma_2^2-\gamma_1^2)
\]
determines \(g\) at most once. The coefficient is nonzero since the parity of the valuations in \(b_i^3\) determines \(b_i\). Substitution in (2) determines \(h^2 d^3\), which determines \(h,d\) by unique factorization with \(d\) squarefree. By the lemma the codegree on the \(b\)-side is \(O_\varepsilon(n^\varepsilon(1+(AC)^2/n))\). The \(d\)-side follows by swapping summands.

For the dual graph, fix distinct \(a_1,a_2\) with a common \(c\). Eliminating the second summand gives the identity
\[
 \kappa_1\left(a_1^2(b_1d_2)^3-a_2^2(b_2d_1)^3\right)=n(d_2^3-d_1^3).
\]
The reduced fraction \(u/v=b_1d_2/(b_2d_1)\) is a unit cube root of \(a_2^2/a_1^2\) modulo \(n\). Write \(b_i=g\beta_i,\ d_i=h\delta_i\), removing the gcd of each pair. Coprimality of the summands gives \((\beta_i,\delta_i)=1\), so \(u=\beta_1\delta_2,\ v=\beta_2\delta_1\) are already reduced. Each \(u/v\) has divisor-many factorizations of this form. The exact equation is now
\[
             \kappa_1 g^3(a_1^2u^3-a_2^2v^3)=n(\delta_2^3-\delta_1^3).        \tag{5}
\]
If the coefficient of \(g^3\) vanishes in a representation, then \(\delta_1=\delta_2=1\). The equation would then give \(a_1^2\beta_1^3=a_2^2\beta_2^3\). Both \(\beta_i\) are squarefree, so equality gives \(\beta_1=\beta_2\) and \(a_1=a_2\), contrary to the choice of distinct vertices. Thus \(g\) is uniquely recoverable. Substitution determines \(h^3c^2\), with \(h\) squarefree, and therefore \(h,c\). Using \(u,v<4BD\) in the angular bound gives \(O_\varepsilon(n^\varepsilon(1+(BD)^2/n))\) for both dual codegrees.

A bipartite graph with sides of cardinality \(m_1,m_2\) and codegree at most \(K\) on each side has
\[
                      e\le m_1+m_2+\sqrt K\min(m_1\sqrt{m_2},m_2\sqrt{m_1})
\]
edges. This follows by summing \(\binom{\deg(v)}2\) on either side and applying Cauchy–Schwarz to the degrees. In both graphs a simple edge has \(O_\varepsilon(n^\varepsilon)\) representations by the fixed-pair bounds below. With \(\sqrt{1+x^2}\le1+x\), the codegree bounds thus give (4).

## Fixed pairs and mixed coordinates

For fixed \(b,d\), set \(L=\kappa_1 b^3,\ M=\kappa_2d^3\) and count \(La+c\sqrt{-LM}\) of norm \(Ln\) in the imaginary quadratic field \(\mathbb Q(\sqrt{-LM})\). Prime-ideal factorization bounds the number of ideals of that norm by \(\tau(Ln)\), and any principal ideal has at most six generators of the specified norm. Consequently \(R\ll_\varepsilon n^\varepsilon BD\), since a nonempty box has \(L\le n\).

For fixed \(a,c\), the binary diagonal cubic \(P b^3+Qd^3=n\) has \(O_\varepsilon(n^\varepsilon)\) solutions, uniformly in the positive coefficients. If the form is irreducible, Evertse's Theorem 1 gives at most \((15\cdot10^6)^{1+\omega(n)}\) cosets for the Thue–Mahler equation over \(\mathbb Q\) with the real place and the primes of \(n\). At a fixed positive value each coset contains at most one positive pair by homogeneity. If reducible, the coefficient ratio is a rational cube. Dividing by the common coefficient and absorbing the two cube coefficients reduces to \(x^3+y^3=n_0\), where \(n_0\mid n\). The divisor \(x+y\) of \(n_0\) determines \(xy\), so each divisor gives at most two ordered pairs. It follows that \(R\ll_\varepsilon n^\varepsilon AC\). [Evertse, Theorem 1](https://pub.math.leidenuniv.nl/~evertsejh/95-Thue-Mahler.pdf)

Fixing \(a,d\) leaves \(P b^3+Q c^2=n\) with \(P,Q>0\). This is an absolutely irreducible affine plane cubic. Indeed, \((n-PX^3)/Q\) has simple complex zeros and is therefore not a square in \(\mathbb C(X)\), so \(Y^2-(n-PX^3)/Q\) is irreducible. Apply the coefficient-uniform Bombieri–Pila bound in a square of side \(2\max(B,C)\). Fixing \(b,c\) gives the symmetric assertion, hence
\[
          R\ll_\varepsilon n^\varepsilon
                \min\left(AD\max(B,C)^{1/3},\ BC\max(A,D)^{1/3}\right).        \tag{6}
\]
The original theorem only requires a square enclosing the rectangle. Its large-side-length threshold depends only on degree and can be absorbed in the constant for the finitely many smaller integer side lengths. [Bombieri–Pila, Theorem 5](https://people.maths.ox.ac.uk/pila/Ovals.pdf)

## Exact optimization

Take logarithms to base \(n\) of the four scales, writing them as \(\alpha,\beta,\gamma,\delta\), and put \(X=\alpha+\gamma,\ Y=\beta+\delta\). By (3),
\[
        2\alpha+3\beta,\ 2\gamma+3\delta\le1,\qquad 2X+3Y\le2.
\]
The five exponent bounds, ignoring the arbitrarily small \(\varepsilon\), are
\[
\begin{gathered}
 Y,\quad X,\\
 \max\left(\beta,\delta,\min(\beta+\delta/2,\delta+\beta/2)+(X-1/2)_+\right),\\
 \max\left(\alpha,\gamma,\min(\alpha+\gamma/2,\gamma+\alpha/2)+(Y-1/2)_+\right),\\
 \min\left(\alpha+\delta+\max(\beta,\gamma)/3,\,
           \beta+\gamma+\max(\alpha,\delta)/3\right).
\end{gathered}
\]
Here \(x_+=\max(x,0)\). If \(X<3/10\) or \(Y<3/10\), one of the fixed-pair bounds gives strictly less than \(3/10\). Otherwise
\[
                             X\le11/20,\qquad Y\le7/15<1/2,
\]
so the dual graph has zero height loss.

First suppose that a kernel exponent is at least \(3/10\). By swapping summands take \(\beta\ge3/10\), so \(\alpha\le1/20\). If \(\gamma<3/10\), the dual bound gives a strict saving because \(\alpha+\gamma/2<1/5\). If \(\gamma\ge3/10\), then \(\delta\le2/15\) and the first mixed bound is at most \(17/60\). Indeed if \(\gamma\ge\beta\), then
\[
       \alpha+\delta+\gamma/3\le1/20+(1-\gamma)/3\le17/60,
\]
and if \(\beta\ge\gamma\), then
\[
       \alpha+\delta+\beta/3\le (1-3\beta)/2+\beta/3+2/15\le17/60.
\]

It remains to take \(\beta,\delta<3/10\). If \(X>1/2\), the kernel graph gives a strict saving because its square-root term has exponent at most \(3Y/4+X-1/2\le X/2\le11/40\). If \(X\le1/2\) and a square-root exponent is at least \(3/10\), take \(\alpha\ge3/10\) by symmetry. Now \(\beta\le2/15\), so the kernel square-root term has exponent at most \(\beta+\delta/2<17/60\), again a strict saving.

We are left with all four exponents below \(3/10\) and both graph height losses zero. For both graph bounds to be at least \(3/10\), their square-root terms must be at least \(3/10\). Each such term is bounded by three-fourths of the sum of its two exponents. Thus \(X,Y\ge2/5\). The inequality \(2X+3Y\le2\) forces \(X=Y=2/5\), and equality in both graph bounds forces \(\alpha=\gamma=\beta=\delta=1/5\). At this point the minimum of the five bounds equals \(3/10\). This establishes the box bound and identifies the unique maximizer of this optimization.

Finally, divide general square-full summands \(m_1,m_2\) by \(h=(m_1,m_2)\mid N\). Every prime occurring with exponent exactly one in a quotient divides \(h\). Extract from each quotient the product \(\kappa_i\) of those primes. The remainder is square-full and can be written uniquely as a square times a squarefree cube, so this is precisely (2) with \(n=N/h\). There are at most \(\tau(N)^3\) choices of \(h,\kappa_1,\kappa_2\), and at most \(O((\log(2N))^4)\) dyadic boxes per choice. Coefficient uniformity and the divisor bound now give (1) after choosing smaller intermediate \(\varepsilon\)'s.

# F-pure thresholds of the Legendre cones

The [Legendre family used for the Gauss–Manin calculation](CM_NEWTON_HORIZONTAL_DENOMINATOR_GROWTH.md) gives a concrete instance of the Hasse-invariant transfer. The threshold computation for this plane cubic is already in [Pagi, Theorem 1.1 and Corollary 3.3](https://arxiv.org/pdf/1706.07309).

## Threshold and congruence classes

For a finite field \(k\) of characteristic \(p>3\) and \(\lambda\in k-\{0,1\}\), put
\[
 f_\lambda=y^2z-x(x-z)(x-\lambda z),\qquad
 H_p(\lambda)=\sum_{i=0}^{m}\binom mi^2\lambda^i,\qquad m=(p-1)/2 .
\]
The affine equation double covers the \(x\)-line with three distinct finite branch points, and the point at infinity \([0:1:0]\) is smooth since the \(z\) derivative there is 1. The projective cubic is thus smooth and adjunction gives trivial canonical class. In Bhatt–Singh's notation the ambient projective dimension is \(n=2\), so both characteristic conditions in Theorem 1.1 hold. The Hasse invariant in that theorem is Frobenius on the line \(H^1(\mathcal O)\). The degree-three hypersurface sequence identifies its scalar, up to trivialization, with the coefficient of \((xyz)^{p-1}\) in \(f_\lambda^{p-1}\) by the local-cohomology description in §4. That coefficient is \(\binom{2m}{m}H_p(\lambda)\): pick exactly \(m\) copies of \(y^2z\), and in the one-variable specialization of the other factor take the coefficient of \(x^{2m}\) in \([x(x-1)(x-\lambda)]^m\), which is \((-1)^m H_p(\lambda)\). The sign from the defining equation cancels it. The binomial prefactor is a \(p\)-adic unit. [Bhatt–Singh, Theorems 1.1 and 4.1, Definition 4.3](https://arxiv.org/pdf/1307.1171v2).

At a point where \(H_p\) is nonzero the order of vanishing on the versal deformation space is zero. At a zero it is positive and at most \(\dim X=1\) by the theorem, hence equals one. Therefore the threshold at the cone vertex is
\[
 \operatorname{fpt}(f_\lambda)=
 \begin{cases}
  1,&H_p(\lambda)\ne0,\\
  1-1/p,&H_p(\lambda)=0.
 \end{cases}
\]
This also determines the global threshold since the divisor outside the vertex is smooth. For \(\lambda=-1\) the coefficient test uses \([x^m](x^2-1)^m\), zero for odd \(m\) and a signed nonzero binomial coefficient modulo \(p\) for even \(m\). Thus the \(j=1728\) member has threshold 1 for \(p\equiv1\bmod4\) and \(1-1/p\) for \(p\equiv3\bmod4\). The \(\lambda=2\) member used in the Gauss–Manin note has the same criterion, by \(x\mapsto z-x\), \(y\mapsto\sqrt{-1}\,y\) over the algebraic closure, which changes the defining equation only by a nonzero scalar. Thresholds are preserved under this finite field extension, as is also immediate from the Frobenius-power monomial membership test.

## Simultaneous comparison of ideals for each rational member

Fix \(\lambda\in\mathbb Q-\{0,1\}\) and any prime \(p>3\) where \(\lambda\) is integral and reduces away from \(0,1\). The multiplier ideal of \(f_\lambda^t\) in characteristic zero is \((f_\lambda^{\lfloor t\rfloor})\) for all \(t\ge0\). Blowing up the origin resolves the cone divisor. The strict transform is smooth and meets the exceptional plane transversely along the smooth cubic. The exceptional divisor has multiplicity 3 and discrepancy 2, so the log canonical threshold is 1. For \(0\le t<1\) the ideal is the unit ideal. For general \(t\), subtracting the integral part of the exponent in the log-resolution formula pulls out \(f_\lambda^{\lfloor t\rfloor}\).

The analogous identity for the test ideals of the reduction holds simultaneously for all \(t\ge0\) if and only if \(H_p(\bar\lambda)\ne0\). In the nonvanishing case the first jump is 1, and principal Skoda iterates the identity from \(0\le t<1\). In the vanishing case the test ideal is proper at the vertex for \(1-1/p<t<1\), while the reduced multiplier ideal is the unit ideal. [Mustaţă–Yoshida, §2, principal Skoda identity](https://arxiv.org/pdf/0706.1124v3).

This comparison concerns the specified principal ideals. The implication in [Mustaţă–Srinivas, Theorem 1.3](https://arxiv.org/pdf/1012.2818) assumes bijectivity of Frobenius on top coherent cohomology on a Zariski-dense set of closed points for every smooth projective characteristic-zero variety and establishes the corresponding general test-ideal comparison conjecture. A further CM/Newton application needs a smooth degree-\((n+1)\) hypersurface in \(\mathbb P^n\), with \(n\ge3\), whose top coherent cohomology line is identified with the proposed arithmetic Frobenius input at good primes. Computing a nonordinary threshold would also need the required order of vanishing on the versal hypersurface space, which a Newton polygon or vanishing along a selected subfamily alone does not determine.

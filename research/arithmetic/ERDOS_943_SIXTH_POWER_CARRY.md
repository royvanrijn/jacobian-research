# Carry and square constraints at the sixth-power modulus

Write \(r(N)\) for the number of ordered representations of \(N\) as a sum of two positive square-full integers. The [cube-modulus product bound](ERDOS_943_CUBE_MODULUS_BOXES.md) gives \(r(N)\ll_\varepsilon N^{7/19+\varepsilon}\). Retaining the complementary squares gives a sixth-power congruence. Its complete Fourier transform can be evaluated on primitive frequencies, and the same calculation identifies a loss in absolute-value completion. A separate bound on the projection modulo the cube improves the full exponent to \(5/14\).

## The carry of a product pair

Put \(q=b^3\). Take two representations with the same \(b\) and put
\[
 C=c_1c_2,\quad D=d_1d_2,\quad S=a_1^2+a_2^2,\quad P=a_1a_2.
\]
The product identity is
\[
 C^2D^3=(N-qa_1^2)(N-qa_2^2)=N^2-NqS+q^2P^2.                  \tag{1}
\]
Assume \((N,b)=1\). For any pair with \(C^2D^3\equiv N^2\bmod q\), its carry class is
\[
 \sigma_N(C,D)=N^{-1}\frac{N^2-C^2D^3}{q}\pmod q.               \tag{2}
\]
A change of either coordinate by a multiple of \(q^2\) preserves the class. The congruence modulo \(b^6\) requires \(S\equiv \sigma_N(C,D)\bmod q\). When \(b^6>2N\), all admissible \(S\) lie in \(0<S<2N/q<q\) for sufficiently large \(N\). Thus at most one integer \(S\) can pass the carry test.

**Theorem.** Suppose \((N,b)=1\) and \((b,30)=1\). For each unit solution pair modulo \(q\), the carries on all its lifts modulo \(q^2\) are exactly uniform. For the complete transform
\[
 \mathcal K(h,k,\ell)=
 \sum_{S\bmod q}\ \sum_{\substack{C,D\bmod q^2\\ C^2D^3=N^2-NqS\bmod q^2}}
 e_{q^2}(hC+kD)\,e_q(\ell S),
\]
at frequencies with \((hk\ell,b)=1\) one has
\[
 |\mathcal K(h,k,\ell)|=
 \begin{cases}
 q^2& h^2 k^3\equiv108N^3\ell^5\bmod q,\\
 0&\text{otherwise}.
 \end{cases}                                                   \tag{3}
\]
In the balanced critical ranges \(X,Y,Z=o(q)\) described below, bounding the nonzero terms of completed rectangular carry counts separately by absolute values leaves a bound of size at least \(q^{1-o(1)}\), even when the carry is required to be a sum of two squares from the prescribed root box.

For the first assertion lift \(C=C_0+q i,\ D=D_0+q j\). Modulo \(q\) the carry has the affine expression
\[
 \sigma_N(C,D)=\sigma_N(C_0,D_0)-N^{-1}
 (2C_0D_0^3 i+3C_0^2 D_0^2j).
\]
One coefficient is already a unit for odd \(b\). Every carry occurs \(q\) times among the \(q^2\) lifts. Also \(\sigma_{N+qv}(C,D)=\sigma_N(C,D)+2v\bmod q\). Averaging over a complete set of target lifts therefore gives the uniform distribution for each fixed product pair. In the critical boxes the natural target interval has length of order \(N\), whereas a complete set of target lifts at fixed \(N\bmod q\) spans \(q^2\). For fixed \(N\) and short coordinate boxes the statement about complete lifts alone supplies no equidistribution.

For the transform, write \(A_h=hN^{-2}\), \(T_k=kN^2\bmod q^2\), and put
\[
 f(t)=A_h t^3+T_k t^{-2},\qquad
 g(t)=N^{-1}(A_h t^3-T_k t^{-2}).
\]
Parametrize by \(t=CD\bmod q^2\). With \(\lambda=N^2-NqS\), the unit coordinates are \(C=\lambda^{-1}t^3,\ D=\lambda t^{-2}\bmod q^2\). Since \((qS)^2=0\bmod q^2\) the phase becomes \(e_{q^2}(f(t))e_q(S(g(t)+\ell))\). Decompose \(t=u+qv\). Orthogonality in \(v\) and \(S\) gives the exact formula
\[
 \mathcal K(h,k,\ell)=q^2
 \sum_{\substack{u\bmod q,\ (u,b)=1\\ f'(u)\equiv0\bmod q\\g(u)+\ell\equiv0\bmod q}}
 e_{q^2}(f(u)).                                                 \tag{4}
\]
Here the value in the sum is independent of the lift of \(u\). For \(b=p>5\) this is stationary phase at \(p^6\) with a simultaneous carry frequency. At unit frequencies, set \(V=A_hu^3,\ W=T_ku^{-2}\). The two stationarity conditions are \(3V=2W\) and \(V-W=-N\ell\bmod q\), so \(V=2N\ell,\ W=3N\ell\). These prescribe \(u^2=T_k/(3N\ell)\) and \(u^3=2N\ell/A_h\). They are compatible precisely when \((u^2)^3=(u^3)^2\), and then determine \(u=u^3/u^2\). This proves (3) directly over the ring and in particular over every squarefree \(b\) coprime to \(30N\) by CRT.

For the loss after completion, take intervals of lengths \(X,Y\) in the product coordinates and a nonnegative carry weight supported in an interval of length \(Z\), of total mass \(W_0\). All these lengths are positive and count integers. Each of the three interval transforms has magnitude at least a fixed multiple of its mass when the corresponding positive frequency is at most \(q^2/(8X), q^2/(8Y), q/(8Z)\), respectively. The same lower bound holds for the weighted transform, by placing all phases in an arc of length at most \(\pi/4\). For each unit \(\ell\) in this range, equation (4) has exactly \(\varphi(q)\) pairs of supported unit residues \((h,k)\bmod q\), one for each \(u\). Each occurs \(\gg (q/X)(q/Y)\) times in the frequency rectangle when \(X,Y=o(q)\). For prime \(b\) the unit count for \(\ell\) is \(\gg q/Z\) as soon as \(q/Z\) grows. For squarefree \(b\), inclusion–exclusion gives that count with a factor \(b^{-o(1)}\) whenever \(q/Z\ge b^\delta\) for fixed \(\delta>0\). Indeed the error is at most \(2^{\omega(b)}\) and the main term is \(\lfloor q/(8Z)\rfloor\varphi(b)/b\). The normalization in Fourier inversion is \(q^{-5}\). The unit terms alone therefore contribute \(\gg b^{-o(1)}qW_0/Z\) to the sum of absolute values. For a carry interval \(W_0=Z\). For its subset of sums \(a_1^2+a_2^2\) with \(a_i\in[A,2A)\), take \(Z\asymp A^2\). This subset has cardinality \(\gg_\varepsilon A^{2-\varepsilon}\), since the number of ordered representations of any \(S\) as two squares is at most \(4\tau(S)\) by Gaussian factorization. The estimate persists with a subpower loss. The claimed loss holds, for instance, at \(B,D_0\asymp N^{2/9},\ A,C_0\asymp N^{1/6}\) with \(b\asymp B,\ X\asymp C_0^2,\ Y\asymp D_0^2,\ Z\asymp A^2\). It concerns a bound by absolute values and gives no lower bound for the actual short-box count.

There is a structured counterexample to the unqualified composite-target carry bound. Take any prime \(b\), put \(N=b^3n\) with \(n\asymp b^{3/2}\), and take \(D=b^2\). Every \(C\asymp b^{3/2}\) then satisfies \(C^2D^3\equiv N^2-Nb^3S\bmod b^6\) for every \(S\). For instance choose \(u=\lfloor b^{3/4}\rfloor,\ n=4u^2,\ S=2u^2\) and \(u^2\le C<2u^2\). Then \(0<C^2D^3<N^2\) and \(S\) is in the desired root scale. This gives \(\gg b^{3/2}\) product pairs and applies to squarefree moduli. It uses the nonunit locus which will be separated below.

## Counting differences in a lifted residue class

Let an \(X\)-by-\(Y\) box here mean two intervals of integer cardinalities \(X,Y\). Its initial position is arbitrary.

**Theorem.** For squarefree \(b\), \((b,30\lambda)=1\) and \(X\le K b^3\), the unit congruence \(x^2 y^3\equiv\lambda\bmod b^3\) in an \(X\)-by-\(Y\) box has at most
\[
 \ll_{K,\varepsilon} b^\varepsilon \min(b,X)\,(1+Y/b^2)\,(1+X/b)^{1/2}        \tag{5}
\]
solutions. In particular, this bounds the unit product pairs surviving the carry test, both per squarefree modulus and after summing over squarefree moduli.

For \(b=1\) the estimate follows by interval counting with \(X\le K\). Otherwise subdivide the \(x\)-interval into at most \(O_K(1)\) intervals of cardinality at most \(b^3\), and the \(y\)-interval into intervals of cardinality at most \(b^2\). There are \(O(1+Y/b^2)\) of the latter. In each resulting box reduce \((x,y)\) modulo \(b\). The number of admissible classes is at most \(\varphi(b)\) by the unit parameter \(t=xy\), and at most \(3^{\omega(b)}X\) by fixing \(x\bmod b\) and counting cube roots modulo squarefree \(b\). This gives \(O_\varepsilon(b^\varepsilon\min(b,X))\) classes.

For each such class \((x_0,y_0)\), choose an integer representative \(y_0\) and the Hensel lift \(x_*\bmod b^3\) solving \(x_*^2 y_0^3\equiv\lambda\bmod b^3,\ x_*\equiv x_0\bmod b\). With \(u=(x-x_*)/b\), \(v=(y-y_0)/b\), the equation in this class is
\[
 u\equiv L v+b Jv^2\bmod b^2,\qquad
 L=-\frac{3x_*}{2y_0},\quad J=\frac{15x_*}{8y_0^2}.                           \tag{6}
\]
Indeed \(x=x_*(1+bv/y_0)^{-3/2}\bmod b^3\) is its unique lift in this class, and the binomial expansion through degree two gives (6). Both \(L\) and \(2J\) are units modulo \(b\). In the trimmed box, \(u\) ranges over an interval of at most \(b^2\) consecutive integers and \(v\) over at most \(b\), with \(|u_1-u_2|\le X/b\).

Consider ordered pairs of distinct points in this class. Their differences \(\alpha=u_1-u_2,\ \beta=v_1-v_2\) satisfy \(\alpha\equiv L\beta\bmod b\). If \(\beta=0\bmod b\), then \(v_1=v_2\), and (6) forces \(u_1=u_2\). Consequently neither \(\beta\) nor \(\alpha\) can be zero for these ordered pairs. For given nonzero \(\alpha\), there are at most two possibilities for \(\beta\), and \((\beta,b)=(\alpha,b)\). Subtracting (6) gives
\[
 \frac{\alpha-L\beta}{b}\equiv J(2v_2\beta+\beta^2)\bmod b.
\]
This has at most \((\beta,b)\) possibilities for \(v_2\) in its interval, and then at most one \(u_2\) by (6). If the class has \(R\) points and \(H=\lfloor X/b\rfloor\), it follows that
\[
 R(R-1)\ll\sum_{1\le|\alpha|\le H}(\alpha,b)
 \le 2\sum_{g\mid b}g\lfloor H/g\rfloor\le 2H\tau(b).
\]
This also covers \(H=0\). Multiplying the resulting bound for \(R\) by the number of classes and blocks gives (5). The contribution of shared divisors of differences is included in \(\tau(b)\), so summing over squarefree bases incurs no additional power. If \(b\) is squarefree and coprime to \(\lambda\) but contains a prime dividing 30, discard this bounded part of \(b\) in the modulus. The same estimate follows with a different \(K\)-dependent constant, since the residual base lies between \(b/30\) and \(b\).

## A uniform exponent of five fourteenths

**Theorem.** For every positive integer \(N\) and \(\varepsilon>0\),
\[
                              r(N)\ll_\varepsilon N^{5/14+\varepsilon}.     \tag{7}
\]

First take coprime summands, allowing fixed positive coefficients
\[
 \lambda_1 a^2 b^3+\lambda_2 c^2 d^3=M
\]
and squarefree \(b,d\). The same bounds hold uniformly in these coefficients. In fact \((\lambda_2 M,b)=1\) because the summands are coprime. The projected unit congruence modulo \(b^3\) has target \(M\lambda_2^{-1}\). Its squared target counts products of two projected solutions. Each product pair has at most \(\tau(c_1c_2)\tau(d_1d_2)\) ordered factorizations, and each original projected solution determines at most one \(a\). These divisor factors are \(O_\varepsilon(M^\varepsilon)\).

Use dyadic scales \(A,B,C_0,D_0\) for \(a,b,c,d\) and swap summands so \(B\le D_0\). Put \(U=BD_0\). The binary quadratic ideal-norm bound with fixed \(b,d\) and the uniform Thue–Mahler bound with fixed \(a,c\) give \(O_\varepsilon(M^\varepsilon\min(U,AC_0))\), as in the cube-modulus note. The same proofs apply to the coefficients \(\lambda_1,\lambda_2\). For the norm equation both positive coefficients are at most \(M\) for a nonempty box, and for Thue–Mahler the prescribed value is \(M\). A reducible binary diagonal cubic reduces to the sum of two positive cubes after removing the common cube-free coefficient. The global size bound \(AC_0\le M/U^{3/2}\) follows by multiplying the two individual inequalities. [Evertse, Theorem 1](https://pub.math.leidenuniv.nl/~evertsejh/95-Thue-Mahler.pdf)

Thus it suffices to take \(M^{5/14}<U<M^{3/7}\). This implies \(C_0^2\le M/D_0^3\le B^3\), so the product box with side cardinalities \(O(C_0^2),O(D_0^2)\) meets the hypothesis of (5), up to a fixed constant. The inequality
\[
 \min(B,C_0^2)(1+C_0^2/B)^{1/2}\ll B^{1/2}C_0
\]
holds by considering \(C_0^2\le B\) and \(C_0^2>B\). Taking the square root of (5), with the divisor multiplicity, and summing over \(O(B)\) bases gives a box bound of
\[
 \ll_\varepsilon M^\varepsilon B^{5/4}C_0^{1/2}(1+D_0/B)
 \ll_\varepsilon M^\varepsilon B^{1/4}D_0 C_0^{1/2}
 \le M^{1/4+\varepsilon}U^{1/4}.
\]
The last term is \(O_\varepsilon(M^{5/14+\varepsilon})\) for the stated \(U\). The number of dyadic boxes costs only a subpower factor.

For general \(N\), write the two square-full summands as \(m_1,m_2\), take \(h=(m_1,m_2)\mid N\), and set \(M=N/h\). Each prime that occurs to exponent exactly one in \(m_i/h\) divides \(h\), because \(m_i\) was square-full. Take \(\lambda_i\) to be the product of those primes. The quotient \((m_i/h)/\lambda_i\) is square-full and has a unique square-times-squarefree-cube expression. This identifies a coprime summand problem of the displayed coefficient form. There are at most \(\tau(N)^3\) choices for \(h,\lambda_1,\lambda_2\), and for each the previous bound applies with \(M\le N\). Absorbing the divisor factor proves (7).

## Exact squares and the pointwise estimate still required

**Proposition.** In the unit case, with \(b^6>2N\), a product pair from two positive representations has a unique carry integer \(S\in(0,2N/q)\). After that carry is chosen, the square conditions
\[
 P^2=(C^2D^3-N^2+NqS)/q^2,\qquad S+2P=u^2,\qquad S-2P=v^2               \tag{8}
\]
with \(P>0,\ u,v\ge0\) recover at most one unordered pair of positive roots \(a_1,a_2\). Together with factor matching and the chosen boxes they are necessary and sufficient for the original product pair to come from two representations.

In fact \(u^2-v^2=4P\) gives \(u>v\) and equal parity. The two roots are \((u+v)/2,(u-v)/2\). The identity \(C^2D^3=(N-qa_1^2)(N-qa_2^2)\) then holds. The two factors are positive because their product is positive and their sum is \(2N-qS>0\). Factor matching consists of checking that \(C=c_1c_2,\ D=d_1d_2\) with \(d_i\) squarefree, and \(N-qa_i^2=c_i^2d_i^3\) for each index. The product identity can redistribute prime exponents between the two factors, making this last check necessary. For example, take the prime \(N=1879\), \(b=7\), \(a_1=1\) and \(a_2=2\). Then \(S=5,\ P=2\), and the two factors \(1536=2^9\cdot3\) and \(507=3\cdot13^2\) have product \(312^2\cdot2^3\). Thus \((C,D)=(312,2)\) passes all three square tests in (8), but each factor has a prime of exponent one. This is a unit example with \(b^6>2N\).

For the coefficient form, the same identity uses the product \(\lambda_2^2 C^2D^3\) and replaces \(q\) in (1) by \(\lambda_1 q\). The normalization by \(h\) in the proof of (7) gives \((\lambda_1,b)=1\), because primes of \(\lambda_1\) have exponent exactly one in the coprime summand. Thus the carry is still unique modulo \(q\) in this application, and the root and factor tests apply with the corresponding coefficients.

For fixed \(b\), denote by \(\mathcal A_b\) the product count with any of the necessary filters above. The factorization bound and Cauchy–Schwarz give a bound of
\[
 M^\varepsilon\left(B\sum_{b\asymp B}\mathcal A_b\right)^{1/2}                \tag{9}
\]
for a box of representations with total \(M\), with an analogous definition of \(\mathcal A_b\) for coprime coefficients. For the exponent \(1/3\) the elementary fixed-pair estimates suffice outside \(M^{1/3}<BD_0<M^{4/9}\). Inside this range the sufficient squarefree average is
\[
 \sum_{b\asymp B}\mathcal A_b\ll_\varepsilon M^{2/3+\varepsilon}/B.           \tag{10}
\]
At \(B,D_0\asymp M^{2/9}\), \(A,C_0\asymp M^{1/6}\), this asks for \(O_\varepsilon(B^{1+\varepsilon})\) product pairs per modulus on average. The bound (5) gives \(O_\varepsilon(B^{5/4+\varepsilon})\) there. At the crossover for (7), \(B,D_0\asymp M^{3/14}\) and \(A,C_0\asymp M^{5/28}\), (5) gives \(O_\varepsilon(B^{4/3+\varepsilon})\); (10) would ask for \(O_\varepsilon(B^{10/9+\varepsilon})\) per modulus on average.

The active ranges for (10) may further be restricted to \(B>M^{4/21}\). To see this, sum the earlier unit product bound over the \(O(B)\) choices of \(b\). Its four terms are
\[
 C_0 D_0 B^{-1/2},\quad B^{7/4},\quad C_0 B^{1/4},\quad D_0 B^{1/4},
\]
up to \(M^\varepsilon\). For \(U>M^{1/3}\), the first and third are at most \(M^{1/2}U^{-1/2}\) and \(M^{1/2}U^{-5/8}\). For \(B\le M^{4/21}\) the second is at most \(M^{1/3}\). If \(D_0 B^{1/4}>M^{1/3}\), fix \(b,c\) in the unit congruence. A specified unit cube has at most \(3^{\omega(b)+1}\) cube roots modulo \(b^3\), so these boxes contribute at most
\[
 \ll_\varepsilon M^\varepsilon B C_0(1+D_0/B^3)
 \ll_\varepsilon M^\varepsilon\left(M^{1/2} B D_0^{-3/2}+M^{1/2}B^{-2}D_0^{-1/2}\right).
\]
The first term is at most \(B^{11/8}\le M^{11/42}\) by \(D_0 B^{1/4}>M^{1/3}\). The second is at most \(M^{1/2}U^{-1/2}B^{-3/2}\le M^{1/3}\). The same unit argument applies for fixed coprime coefficients. Consequently \(b^6>2M\) in every remaining box for sufficiently large \(M\), and uniqueness of the carry does hold throughout the remaining range.


For a direct squarefree average in the balanced one-third box, split the \(Y\)-interval into \(O(1)\) blocks of length at most \(b^2\), as in the proof of (5). There are \(O(B^2)\) admissible modulo-\(b\) classes across all blocks and moduli. Write \(K_{b,\rho}\) for the number of carry-admissible pairs in a class \(\rho\) including its block. The estimate
\[
 \sum_{b\asymp B}\ \sum_{\rho} K_{b,\rho}(K_{b,\rho}-1)
             \ll_\varepsilon B^{2+\varepsilon}                              \tag{11}
\]
would imply \(\sum_b\mathcal A_b\ll_\varepsilon B^{2+\varepsilon}\), since \(K\le 1+\sqrt{K(K-1)}\) and Cauchy–Schwarz applies across the classes.

The candidates for (11) can be enumerated by small differences. For \(b=p>5\) prime put \(\alpha=(C_1-C_2)/p,\ \beta=(D_1-D_2)/p\), with \(0<|\alpha|\ll B^{1/2}\), \(0<|\beta|<p\). Fix them and reduce the difference relation modulo \(p\). On the unit locus
\[
 2D_2\alpha+3C_2\beta\equiv0\bmod p,\qquad
 D_2^5\equiv \frac{9\beta^2 N^2}{4\alpha^2}\bmod p.
\]
This determines at most five base classes. For each one, (6) and its difference equation determine at most one candidate lifted pair in a block, and (2) gives the two actual carries to test. For squarefree \(b\), put \(g=(\alpha,b)=(\beta,b)\). There are at most five choices for the base class at each prime of \(b/g\) and at most \(p\) at each prime \(p\mid g\). After choosing a base class the difference equation allows at most \(g\) lift digits. Thus the number of candidates per difference is at most \(5^{\omega(b)}g^2\). Summing this weight over \(g\mid b\) and the differences, which are multiples of \(g\), bounds the candidate count by \(O_\varepsilon(B^{3/2+\varepsilon})\) per \(b\). This is the same size as the direct class count before testing carries. To obtain (11) this way requires a saving of \(B^{1/2}\) in the count of difference candidates across \(b\) after testing the two carries or the exact squares.

**Open question.** Obtain the short-box carry or exact-square saving (10) on the unit locus, including the unbalanced boxes, to reach the uniform \(1/3\) exponent. The block bound (5) only uses the lift through \(b^3\). At the sixth-power modulus formula (4) gives support and magnitude, while the required cancellation of its phases, or another bound using (8) with factor matching, is open.

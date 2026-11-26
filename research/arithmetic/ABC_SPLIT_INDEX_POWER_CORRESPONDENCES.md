# Joint index power classes under projective correspondences

The split cubic index form of a common-lock lattice gives a curve when multiplied by the index forms of prescribed projective relatives. The branch set is read from the multiplicities of the actual linear factors. The count uses the power class of the product and allows cancellation between rows. Even a single affine shift in an exact target fibre can supply the fourth branch needed for a parity count.

## A joint count from the branch set

Fix \(X\ge2\), a positive integer \(h\), a prime \(q\mid 3h\), and \(C>0\). For each slot \(i\) fix a lattice of the type in (SL1) and its normal-form split cubic
\[
                        f_i=L_{i,1} L_{i,2} L_{i,3}.
\]
The fixed divisors \(F_i,G_i,H_i\) and index \(J_i\) are at most \(X^C\). Use positive primitive original rows in these lattices with complete letters at most \(X\). They may in addition have exact locks and individual original radical caps. Their primitive parameter directions are \([v_i]=[x_i:y_i]\), and \(I_i=f_i(v_i)\) is the index of the corresponding split generator. The parameter height is polynomial in \(X\) uniformly, by the normal form or by inverting the fixed basis. [Split cubic normal form and index (SL1–SL6)](ABC_SPLIT_CUBIC_LOCK_LATTICES.md)

Prescribe integral nonsingular matrices \(A_i\) of coefficient height at most \(X^C\), and require \([v_i]=[A_i v_1]\) for \(A_1=1\). The matrices may vary with \(X\). For an integer \(D\ne0\) and a prime pool \(\mathcal P\), put
\[
              d_q(D;\mathcal P)
               =\prod_{\substack{p\notin\mathcal P\\q\nmid v_p(D)}}p ,
\]
where valuations use absolute values. Write \(B_i\) for the three distinct zero directions of \(f_i\). At a direction \(r\) in the union of the \(A_i^{-1}B_i\), let \(m_r\) be the number of occurrences of \(r\). Let \(m\) be the number of these multiplicities not divisible by \(q\). Set
\[
                             g=(q-1)(m-2)/2.
\]
Assume \(g\ge2\), or \(g=1\) with \(q\in\{2,3\}\). For every \(B\ge1\), the number of original row tuples satisfying the prescribed correspondences and
\[
                      d_q\left(\prod_i I_i;\mathcal P\right)\le B
\]
is
\[
                              O_{h,q,C,\varepsilon}(B X^\varepsilon).
                                                               \tag{IC1}
\]
The prescribed pool has polynomial-height product, with a fixed exponent included in \(C\). One can enlarge it by the primes of the matrix determinants, the fixed lattice coefficients, and the pair resultants of the distinct linear factors. This remains a polynomial-height pool uniformly.

To prove this, a primitive representative of \([A_i v_1]\) is obtained by dividing by the positive gcd of its two entries and choosing a sign. That gcd divides \(|\det A_i|\) by the adjugate identity and primitivity of \(v_1\). Thus away from these determinant primes the index valuations of the tuple sum to the valuation of the homogeneous form
\[
                             G(v_1)=\prod_i f_i(A_i v_1).
\]
Enlarge the pool to include the determinant primes. The signed \(q\)-free representative \(D\) of \(G(v_1)\) has at most \(B X^\varepsilon\) possibilities. First enumerate its outside radical, then assign the bounded exponents on those primes and on the fixed pool. The sign contributes at most a factor two. All relevant quantities have polynomial height. The fixed-support and divisor estimates justify the subpower allowances uniformly.

Over \(\mathbb Q\), collect the proportional linear factors of \(G\), and remove their \(q\)-th powers. The resulting form has degree divisible by \(q\), since \(q\mid3h\), and multiplicities in \(\{1,\ldots,q-1\}\). Consequently its \(q\)-power class is well-defined on a projective parameter. After an integral unimodular change of coordinates we may put all its branch points off infinity, using a bounded list of changes since the branch count is fixed. For each \(D\) the admissible parameters now lift to a curve
\[
                      W^q=d\prod_{j=1}^{m}(Z-r_j)^{e_j},
                                                               \tag{IC2}
\]
where the distinct \(r_j\) and the nonzero \(d\) are rational and have polynomial height. The weighted degree is divisible by \(q\); infinity is unramified. The zeros are totally ramified and Riemann–Hurwitz gives genus \(g\). No input index vanishes.

For \(g\ge2\) work over the fixed field \(\mathbb Q(\zeta_q)\). The curve (IC2) is split, and the rational branch-point differences give the full kernel of \(1-\sigma\) on its Jacobian, where \(\sigma(W)=\zeta_q W\). In detail, principal relations supported on the branch points come from deck eigenfunctions \(W^k R(Z)\), whose branch orders modulo \(q\) are proportional to \((e_j)_j\). Conversely that relation arises from the defining equation, after adjusting by powers of \(Z-r_j\) to cancel the pole at infinity. Thus the degree-zero branch distributions modulo \(q\) give a kernel of rank \(m-2\). This is its full size since the deck trace is zero and the cyclotomic identity gives \([q]=(1-\sigma)^{q-1}\) times an automorphism. Include in \(S\) the primes of \(q\), the numerator and denominator of \(d\), the denominators of all \(r_j\), and the numerators of all nonzero differences \(r_i-r_j\). Away from \(S\), the branch sections are integral and disjoint, \(d\) is a unit, and the tame-cover model has good reduction. Locally at a branch point one removes the unit in the equation after a finite étale extension and normalizes \(W^q=(Z-r_j)^{e_j}\), using \(\gcd(q,e_j)=1\). At infinity the total exponent is divisible by \(q\), so the unit chart is étale. The isogeny Kummer classes are unramified at these good primes. The class group and unit bounds over the fixed base field give rank \(O(1+\omega(S))\), where \(S\) is a polynomial-height product of the exceptional rational primes. Uniform Mordell–Lang for the Abel–Jacobi embedding, using a rational branch point, gives \(X^\varepsilon\) rational parameters per \(D\). [Split cyclic descent](ABC_BINARY_POWER_TWISTS.md#projective-covers-from-conjugate-ratios) · [Gao–Ge–Kühne, Theorem 1.1′](https://arxiv.org/html/2105.15085v4)

If \(q=2,g=1\), move one of the four rational branches to infinity by a rational fractional-linear change of polynomial height. The twist is a split cubic of polynomial coefficient height. Scaling its leading coefficient gives a monic Weierstrass model with rational \(2\)-torsion, and the original directions give points of polynomial height. Dujella's bound counts them at \(X^\varepsilon\) cost per \(D\). If \(q=3,g=1\), there are three rational branches and all three residues \(e_j\bmod3\) agree. Move them to \(0,1,\infty\) and, if necessary, invert the power class to get \(W^3=d Z(Z-1)\) with polynomial-height \(d\). Set \(a=4/d,\ U=aW,\ V=a(2Z-1)\). Then \(V^2=U^3+a^2\) and \((0,a)\) is rational of order three. The \(U\)-height is polynomial, and each \(U\) recovers \(Z\) in boundedly many ways. Dujella applies again. All such height statements are polynomial at fixed \(h,q,C\); in particular the height of a rational \(q\)-th root is the \(q\)-th root of the height. Finally one projective parameter fixes each original row, since its primitive lattice pair is determined up to sign and positivity selects the sign. This proves (IC1). [Dujella, Theorem 1.1 and Corollary 1.3](https://arxiv.org/html/2312.03655#S1)

## Two relatives, exact shifts, and the sharp parity boundary

For \(h=2\), let \(r\) be the number of shared points of \(B_1\) and \(A_2^{-1}B_2\). For \(q=2\) they cancel from the branch set, so the genus is \(2-r\) for \(r=0,1,2\). Thus (IC1) applies exactly on the positive-genus side of these overlap cases. For \(q=3\) every branch remains since the original multiplicities are one and shared branches have multiplicity two. There are always at least three branch points, so (IC1) holds for every prescribed map in the cubic-class count. In the full branch coincidence case it uses the displayed rational \(3\)-torsion model rather than genus at least two. If a dictionary of \(D\) maps and lattice records is prescribed, summing costs \(D\).

A physical shift in a fixed exact target fibre is particularly simple. Fix a lattice \(\Lambda\), a positive target \(C_0\), and a nonzero integer shift \(h_0\) with \(|h_0|<C_0\). Count two primitive rows of target \(C_0\) with \(a_2=a_1+h_0\). In the original direction \(t=a/c\) their map is \(t\mapsto t+h_0/C_0\). The branch sets are \(\{0,1,\infty\}\) and its inverse translate. They share only infinity, giving the genus-one product bound with \(q=2\), even if the individual index kernels are large. The actual shift must satisfy
\[
                                FGJ\mid h_0 .
\]
Indeed \((h_0,-h_0)\in\Lambda\) forces \(FG\mid h_0\), and (SL2) then reduces to \(J\mid h_0/(FG)\) using \(Gu+F=Hv\) with \(v\) a unit modulo \(J\). This is the same primitive axis condition as in the split order. Conversely the divided step gives exactly the differences along the target fibre in \(\Lambda\), before primitivity restrictions. Thus if the radical of \(h_0/(FGJ)\) is at most \(A\), a moving shift can be included at cost \(A X^\varepsilon\) by the one-variable supported-power count. With target \(C_0\) fixed, the pair count is \(O(AB X^\varepsilon)\). The same result works with a fixed complete summand and the corresponding permuted shift step. It can be applied with all three original caps imposed. In a fibre of exponent \(r\) in base \(n=X^L\), the joint shift and parity caps give a power-smaller pair fraction if \(A=n^\alpha,\ B=n^\beta\) and \(\alpha+\beta<2r\) with fixed margin. This requires a cap on the actual divided shifts and the joint actual index parity.

There is a separate Pell bound which requires no shift dictionary. Let a fixed-target fibre have \(T\) rows. For every fixed first row the joint parity cap allows at most \(B X^\varepsilon\) squarefree kernels of the product \(a_2b_2\) of the next row. Indeed, outside the fixed pool of coefficient and target primes the first index class is fixed, and the capped product defect recovers the second class by symmetric difference. The inside prime choices cost a subpower. For a fixed positive squarefree kernel \(D_2\), write \(a_2b_2=D_2 w^2\). Then
\[
                    (C_0-2a_2)^2+D_2(2w)^2=C_0^2 .
\]
The coefficient-uniform binary quadratic norm estimate gives \(X^\varepsilon\) pairs at polynomial height. When instead a complete summand \(A_0\) is fixed, use \((2b_2+A_0)^2-D_2(2w)^2=A_0^2\), now with \(b_2c_2=D_2w^2\). This includes the real Pell and split cases. Thus the pair bound in a fixed-complete-letter fibre is \(O(B\min(A,T)X^\varepsilon)\), with the shift cap used only for the \(A\) alternative. [Browning–Verzobio, Lemma 2.3, citing Vaughan–Wooley](https://arxiv.org/html/2608.24512v1#S2.SS1) · [Uniform elementary Pell recovery in the nonsplit real case](ERDOS_943_PELL_RECOVERED_KERNEL_GCD_WEIGHT.md#1-a-uniform-elementary-norm-lemma)

Two shared roots already admit genuine large families for parity. Use the original lattice \(\mathbb Z^2\) in both slots and \(t\mapsto t/2\), whose inverse branch set shares \(0,\infty\) with the original. Choose coprime positive integers \(u,v\) with \(v\) odd, \(v/4<u<v/2\), and \(v\le\sqrt{X/2}\). Set
\[
                   a=v^2-2u^2,\qquad c=v^2-u^2 .
\]
Then the two primitive rows are \((a,u^2,c)\) and \((a,v^2,2c)\). Indeed \(a\) is odd, and \(\gcd(a,c)=\gcd(v^2,u^2)=1\). Their projective ratios are \(t=a/c\) and \(t/2\), and both heights are at most \(X\). The product of their indices in the split order is \(2a^2 c^2 u^2 v^2\). Its parity defect is one outside the pool \(\{2\}\). The primitive pairs \((u,v)\) in this fixed cone with odd \(v\) number \(\gg X\) by Möbius inversion on the gcd and the positive coprime density. They give distinct ratios since \(u/v\) is positive and \(t=(1-2(u/v)^2)/(1-(u/v)^2)\) is injective there. Both original fibres and the determinant directions are actual arithmetic rows. The example does not impose small caps on their original radicals. It prevents extending the fixed-class subpower claim to the overlapping genus-zero case using only the square class of the index product.

The completion estimates for split covers and for four linear forms are classical inputs of this argument. The additional constraints here are a joint index class of separate original rows, their prescribed projective correspondence, the exact divided shift inside a split-cubic lock lattice, and the branch-overlap cutoff. A full endpoint use would need an arithmetic extraction producing sufficiently many of these cheap joint records. [Secant defect counts with fixed anchors](ABC_SECANT_VALUATION_DEFECT_GRAPHS.md) · [A fourth determinant on a single row](ABC_EUCLIDEAN_GRID_POINT_COUNTS.md#a-parity-count-using-a-fourth-determinant)
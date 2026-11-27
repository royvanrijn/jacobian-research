# Coherent prime paths from three marked norm values

The [cofactor Hessian](ABC_COHERENT_NORM_COFACTOR_CUBIC.md#exact-cubic-and-quadratic-forms) gives three primitive representations of \(np_i\) by one quadratic norm form, where \(n=p_1p_2p_3\). A product condition on the three representing ideals recovers the coherent path. In the range governing the second moment, the prime coset can be selected so that the product condition follows from the representations. Fixing the form and a squared prime puts its representing vector in a lattice of index \(p_i^2\).

## Exact reconstruction and automatic orientation

Let \(\mathcal O_D\) be the imaginary quadratic order of discriminant \(D<0\) in \(K=\mathbb Q(\sqrt{-d})\), \(d>3\) squarefree. Let \(A\) be an invertible integral ideal of \(\mathcal O_D\) and put \(F_A(m)=N(m)/N(A)\). Take distinct odd primes \(p_1,p_2,p_3\) split in \(K\) and prime to \(D\). A primitive vector of \(A\) is one outside \(\ell A\) for every rational prime \(\ell\).

**Proposition (Reconstruction from marked norm representations).** Suppose \(\gamma_i\in A\) are primitive and \(F_A(\gamma_i)=np_i\) for \(i=1,2,3\), where \(n=p_1p_2p_3\). Set \(J_i=(\gamma_i)A^{-1}\). The product condition
\[
                              J_1J_2J_3=(n^2)                       \tag{1}
\]
is equivalent to the existence of one orientation \(\mathfrak p_i\) at each \(p_i\) for which \(J_i=\overline{\mathfrak p}_i^{\,2}\mathfrak p_j\mathfrak p_k\), where \(\{i,j,k\}=\{1,2,3\}\). Under (1) the three marks reconstruct coherent prime-ideal cube witnesses. In a positively oriented basis of \(A\) put \(\Delta_{ij}=\det_A(\gamma_i,\gamma_j)\). The reconstructed coefficients relative to \(\sqrt D\) are
\[
                   x=-\Delta_{12}/p_3,\qquad
                   y=-\Delta_{23}/p_1,\qquad
                   z=-\Delta_{13}/p_2.                               \tag{2}
\]
They give the cofactor order of discriminant \(D\) exactly when \(\gcd(x,y,z)=1\). Every coherent path with cofactor class \([A]\) in that order supplies such marks.

Suppose moreover that the three primes admit orientations whose maximal-order classes lie in a common coset \(c\in C/C[3]\), where \(C=\mathrm{Cl}(K)\), with \(4c\ne0\). If the maximal-order image of \([A]\) is nonzero and is killed by \(3\), then (1) follows from the three primitive norm representations.

Indeed \(J_i\) is integral of norm \(np_i\), coprime to the conductor. Since \(\gamma_i\) is primitive, \(J_i\) is contained in no \(\ell\mathcal O_D\) for a rational prime \(\ell\), so just one orientation at each \(p_j\) occurs in it. At \(p_i\), the exponents in the three \(J_j\) are \(2,1,1\). To give two copies of each orientation in (1), the two single occurrences must use the orientation opposite the double occurrence. This gives the simultaneous pattern.

For automatic alignment, use the proposed maximal-order orientations and write their classes as \(c_i\). The image of each \(J_i\) in \(C/C[3]\) is zero. It is also \((2s_i+s_j+s_k)c\), where the signs are chosen in that row to match the support of \(J_i\). The integer \(2s_i+s_j+s_k\) belongs to \(\{-4,-2,0,2,4\}\); since \(4c\ne0\), it is zero. Thus each row is either \(\overline{\mathfrak p}_i^{\,2}\mathfrak p_j\mathfrak p_k\) or its conjugate. The maximal-order class of the displayed choice is \(c_1+c_2+c_3-3c_i\), independent of \(i\) since the \(3c_i\) coincide. Its two signs are distinct because the maximal-order image of \([A]\) is nonzero \(3\)-torsion, and the maximal-order classes of all \(J_i\) equal its inverse. Hence either all rows are the displayed choice or all rows are conjugate. In either case (1) holds already in \(\mathcal O_D\).

For reconstruction take
\[
 \alpha_{12}=p_1p_2^2\,\gamma_1/\gamma_2,\qquad
 \alpha_{23}=p_2p_3^2\,\gamma_2/\gamma_3,\qquad
 \alpha_{13}=p_1p_3^2\,\gamma_1/\gamma_3 .
\]
For example \(J_1/J_2=(\overline{\mathfrak p}_1\mathfrak p_2)^3/(p_1p_2^2)\). The three elements are integral generators of the appropriate oriented cube ideals and satisfy \(\alpha_{12}\alpha_{23}=p_2^3\alpha_{13}\). Multiplying \(A\) and its marked vectors by \(p_1p_2^2p_3/\gamma_2\) gives the actual cofactor ideal \(I=\overline{\mathfrak p}_1\mathfrak p_2^2\overline{\mathfrak p}_3\) and the three vectors \(p_3\alpha_{12},\,p_1p_2^2p_3,\,p_1\overline{\alpha}_{23}\). This homothety preserves the oriented determinants in the ideal basis. Computing them in the integral order basis and dividing by \(N(I)\) gives (2). For the \(13\) determinant use \(\overline{\alpha}_{12}\overline{\alpha}_{23}=p_2^3\overline{\alpha}_{13}\). Multiplying \(\sqrt D=\sqrt{-d}\sqrt{|D|/d}\) by the common gcd of \(x,y,z\) gives the normalized square-root coordinate of the reconstructed path. Thus its cofactor order has discriminant \(D\gcd(x,y,z)^2\). When the gcd is one, its cofactor ideal is the ideal just reconstructed. The forward direction is the three primitive marks of the cofactor Hessian. For selected rank-one records, apply the signed coefficient windows, unequal cube colors, and the rank-one test to the recovered data.

## A count with one or two prescribed primes

**Proposition (Fixed-class counts and common-prime multiplicity).** Fix \(D\) and the oriented cofactor ideal class in \(\operatorname{Pic}(\mathcal O_D)\). With \(p_i\in[Z,2Z)\), \(Z\ge2\), the number of coherent paths whose prime set contains one prescribed prime \(p\) is
\[
                           O(1+Z^2/\sqrt{|D|}).                     \tag{3}
\]
For two prescribed distinct primes \(p,q\) it is \(O(1+Z/\sqrt{|D|})\). In particular if \(\sqrt{|D|}>32 Z\), there is at most one unordered triple for each choice of the local prime ideal orientations at \(p\) and \(q\) in the squared mark for \(p\). The number of unordered triples containing \(p,q\) is then at most four. Without a prescribed prime the bound is \(O(Z+Z^3/\sqrt{|D|})\).

Fix integral \(A\) in the class. In the norm \(F_A\), the squared mark for \(p\) lies in \(\mathfrak q_p^{\,2} A\) for one of the two primes \(\mathfrak q_p\) over \(p\), and is still primitive in \(A\). Its squared length is at most \(16Z^4\). The Euclidean covolume in the norm \(F_A\) of \(A\) is \(\sqrt{|D|}/2\). Imposing this square factor gives covolume \(p^2\sqrt{|D|}/2\). If the candidate vectors in the resulting disk span a line, there are at most two primitive vectors of \(A\) on it. If they span the plane, successive minima in that disk give \(O(Z^4/(p^2\sqrt{|D|}))\) candidates. Each candidate determines the unordered pair of remaining primes by factoring \(F_A(\gamma)/p^2\). For a fixed order and prime set there are only boundedly many ordered coherent paths, since the orientations and the unit signs specify each generator.

If \(q\) is prescribed too, the squared mark for \(p\) has index \(p^2q\) in \(A\) after both prime orientations in its support are fixed. The same count gives \(O(1+Z/\sqrt{|D|})\). When \(\sqrt{|D|}>32Z\), that covolume exceeds \(16Z^4\), so two candidate vectors cannot span the plane. Their primitive norms coincide and determine the third prime. There are two choices of orientation at each of \(p,q\), proving the stated bound. Summing (3) over \(p\) proves the last assertion.

At the crossover \(V=Z^{5/4}\) the bound is \(O(1+V/Z)\) per prescribed prime and \(O(V)\) per cofactor class throughout the field window. This gives a power saving on any family occupying \(O(V^{2-\sigma})\) pairs \((D,[A])\), for fixed \(\sigma>0\). In balanced boxes \(|D|\asymp X\). The ordinary average of \(3\)-torsion over such orders only gives \(O(X)\) classes, by [Bhargava–Varma, Theorems 2–3](https://arxiv.org/abs/1401.5875). Inserting it into this count gives \(O(XV)\), above \(V^3\). The joint occupancy of the form classes remains an arithmetic counting problem.

## Availability in the second-moment range

**Corollary (Selection outside quotient classes killed by four).** Fix \(\tau>0\) and use the original prime-path construction for \(d\in[X,2X)\), \(H<h_3(-d)\le2H\), \(H\le X^{1/4-\tau}\), and the original prime-range conditions. The per-field lower bound for noncollinear, unequal-color paths holds after requiring \(4c\ne0\). Taking \(0<\tau<1/28\) includes the original second-moment crossover in the automatic orientation range.

The classes over \((C/C[3])[4]\) number \(h_3(-d)|C[4]|\), as multiplication by four is invertible on \(C[3]\). Gauss genus theory gives \(|C[2]|=2^{\omega(D_K)-1}\), where \(D_K\) is fundamental. See [Fouvry–Klüners, Section 1](https://math.uni-paderborn.de/fileadmin-eim/mathematik/AG-Computeralgebra/Publications-klueners/ranks.pdf). Consequently \(|C[4]|\le|C[2]|^2\ll_\varepsilon X^\varepsilon\). The selected primes inject into \(C\) under the original size bound \(Z\le X^{1/3}\). Hence omitting these cosets deletes at most \(O_\varepsilon(HX^\varepsilon)=o(Z/\log Z)\) primes. Use the [noncollinear path selection](ABC_COHERENT_NORM_COFACTOR_CUBIC.md#irreducible-path-selection) on the remaining cosets. This gives the same lower bound, and the image of the cofactor class is nonzero on these paths.

## A nonmaximal-order example

Take \(K=\mathbb Q(\sqrt{-10015})\) with primes \(7,17,23\). In the maximal order, the prime forms \((7,-3,358),(17,7,148),(23,-17,112)\) have common cube class \((8,-1,313)\), whose fourth power is \((40,15,64)\) in reduced form. These identities can be checked by applying PARI/GP's quadratic-form powering function to the displayed forms. The oriented edges in the order \(\mathcal O_{-40060}=\mathbb Z[\sqrt{-10015}]\) are
\[
 \alpha_{12}=688+11\sqrt{-10015},\quad
 \alpha_{23}=7244+27\sqrt{-10015},\quad
 \alpha_{13}=409+20\sqrt{-10015}.
\]
Their cofactor norm form and three primitive marks are
\[
 F(S,T)=59S^2+58ST+184T^2,\qquad (11,-10),\quad(-11,-13),\quad(-23,-10).
\]
The norms are \(19159,46529,62951\), equal to \(2737\) times \(7,17,23\). The determinants in order \(12,23,13\) are \(-253,-189,-340\), giving \(11,27,20\) in (2). One ideal basis giving these coordinates is \(-1026+13\sqrt{-10015},-2711-11\sqrt{-10015}\). These data permit the order and edge reconstruction to be checked directly. The reflected trinomial of \(\Theta=\alpha_{12}\overline{\alpha}_{23}\) is irreducible modulo \(7\), so the maximal-order image of the cofactor class is nonzero.

## Failure of automatic orientation in the identity coset

In \(K=\mathbb Q(\sqrt{-4027})\), use \(p_1=1009,p_2=13,p_3=17\). A prime above \(1009\) is principal by the element \((3+\sqrt{-4027})/2\). The prime cubes above \(13,17\) have generators with doubled coefficient pairs \((69,1),(125,1)\). Their classes have order three and are on different lines: equality up to sign would give a nonrational principal generator of norm \(13\cdot17<4027/4\). Thus the prime classes share the identity coset in \(C/C[3]\), and \(A=\overline{\mathfrak p}_1\mathfrak p_2^2\overline{\mathfrak p}_3\) has nonzero class.

Start with the ideal orientations in (1), then replace the single factor \(\mathfrak p_1\) in \(J_2\) by \(\overline{\mathfrak p}_1\). Each \(J_i\) still has class \([A]^{-1}\), since \(\mathfrak p_1\) and its conjugate are principal. Thus each \(A J_i\) has a generator giving a primitive mark of norm \(np_i\). At \(p_1\) the product now has exponents \(1,3\) instead of \(2,2\). These three particular primitive representations cannot reconstruct a coherent path, despite having the correct norms in one nontrivial cofactor class.

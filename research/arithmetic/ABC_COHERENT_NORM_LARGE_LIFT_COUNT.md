# Counting coherent prime triangles with a large lifted congruence

In a [coherent norm triangle](ABC_COHERENT_NORM_PENCILS.md#exact-reconstruction-from-three-signed-coefficients), write \(P_i=p_i^3\), \(u_{ij}=\operatorname{Tr}(\alpha_{ij})\), and
\[
 (v_{12},v_{23},v_{13})=\delta(x,y,z),\qquad
 \delta=\gcd(v_{12},v_{23},v_{13})>0.
\]
The signed nonzero \(x,y,z\) are pairwise coprime and each is coprime to its two incident primes. At edge \(12\) define \(g_x=\gcd(x,u_{23},u_{13})\), with \(g_y,g_z\) defined by relabeling. The [square node variable](ABC_COHERENT_NORM_CAYLEY_CUBIC.md#exact-dyadic-content-and-a-lifted-congruence) at each edge has square root equal to this gcd if \(E=d\delta^2\) is odd. If \(E\) is even, the square root lies between half the gcd and the gcd. Thus bounds using a lower threshold on these gcds also bound the corresponding square-root sector.

## A uniform sector of the rank-one count

**Theorem (Large lifted congruence bound).** Let \(Z\ge2\) and \(V\ge1\). Count coherent triangles with distinct odd primes in \([Z,2Z)\) and \(1\le |x|,|y|,|z|\le V\), allowing the squarefree radicand and orientations to vary. In a dyadic box \((|x|,|y|,|z|)\asymp(B_x,B_y,B_z)\), the number with \(g_x\ge L\ge1\) is
\[
 \ll_\varepsilon (ZV)^\varepsilon Z^2 B_x
            \left(1+\frac{B_y B_z}{B_x L}\right)
            \left(1+\frac{Z}{B_y B_z}\right).                        \tag{1}
\]
For the occupied [rank-one congruence records](ABC_QUADRATIC_THREE_TORSION_FERMAT_TAIL.md#rank-two-and-the-primitive-coefficient-plane), the sector \(\max(g_x,g_y,g_z)\ge L\) has size
\[
 \ll_\varepsilon (ZV)^\varepsilon
                 \min\left(V^3,\ Z^2 V+\frac{Z^2V^2}{L}\right).       \tag{2}
\]
In particular at \(V=Z^{5/4}\), the sector with a lift of size at least \(V^{3/4}\) contributes \(O_\varepsilon(V^{57/20+\varepsilon})\), a saving of \(V^{3/20}\) up to the subpower allowance. This is uniform in any imposed field window or selected witnesses.

Fix \(p_1,p_2\) and \(x\). There are divisor-many possibilities for \(g=g_x\mid x\), and there are subpower many residue classes of \(u=u_{12}\) modulo \(x^2\), since the norm equation gives
\[
                           u^2\equiv4P_1P_2\pmod{x^2}.
\]
At odd primes of \(x\) the right side is a unit and each prime power has at most two root classes. At \(2^a\Vert x\), write \(u=2u'\); solving \(u'^2\equiv P_1P_2\pmod{2^{2a-2}}\) gives at most eight choices of \(u\bmod 2^{2a}\), including \(a=1\).

The exact [anchor pencil](ABC_COHERENT_NORM_PENCILS.md#a-pencil-through-two-split-primes), using vertices \(1,2\), yields
\[
        x^2p_3^3=P_2z^2-uzy+P_1y^2,\qquad
        xu_{23}=2P_2z-uy.                                          \tag{3}
\]
The second identity also follows from the normalized trace formulas by substituting the first. As \(g\mid u_{23}\), the completion pair \((y,z)\) satisfies \(2P_2z\equiv uy\bmod |x|g\). The chosen class of \(u\bmod x^2\) determines this homogeneous lattice because \(g\mid x\). The homomorphism \((y,z)\mapsto 2P_2z-uy\bmod |x|g\) has image containing the residue of \(2P_2\). Since \(p_2\nmid x\), the index of its kernel is at least \(|x|g/2\). Thus the congruence has a rank-two solution lattice of covolume \(\gg B_xL\).

Scale the two coordinates by \(B_y^{-1},B_z^{-1}\). If the candidate vectors span a line, the primitivity \((y,z)=1\) gives at most two choices. If they span the plane, the two successive minima of the scaled lattice in a fixed ball are bounded, and the lattice count is \(O(B_y B_z/(B_xL))\) by the planar lattice-point bound. This gives \(O(1+B_y B_z/(B_xL))\) choices in all cases.

For a chosen pair, the first equality in (3) gives the separate cubic congruence
\[
                      x^2 p_3^3\equiv P_2z^2+P_1y^2\pmod{|yz|}.     \tag{4}
\]
The coefficient \(x^2\) and the right side are units modulo \(|yz|\). Indeed \((x,yz)=(y,z)=1\), \(p_2\nmid y\) and \(p_1\nmid z\). There are at most three cube-root classes per prime power and hence subpower many classes in (4). In the interval of length \(Z\) these give \(O_\varepsilon((ZV)^\varepsilon(1+Z/(B_yB_z)))\) choices for \(p_3\). Each choice determines the exact signed \(u\) by (3), since \(yz\ne0\), so the field defect \(E=(4P_1P_2-u^2)/x^2\) and its squarefree part \(d\) are determined. The other traces follow from reconstruction. Any remaining choices of orientations and signs cost at most an absolute constant. This proves (1).

For (2), each signed coefficient box already contributes at most \(O_\varepsilon((ZV)^\varepsilon B_x B_y B_z)\) occupied rank-one records by the full-modulus lattice argument. Here \(j=2z\) in the original record construction. If \(B_y B_z<Z\), the rank-one bound is at most a subpower multiple of \(VZ\). Otherwise (1) is at most a subpower multiple of \(Z^2(V+V^2/L)\). Sum the dyadic boxes and allow all three placements of the large lift. This proves (2).

The proof uses the same local cube-root count and lattice-point principle as [Heath-Brown–Pierce, Lemma 4.3 and Section 4.3](https://arxiv.org/pdf/1409.3177). Their lattices count the outer primes from the fixed congruence parameters of a path. Here the lifted lattice counts the two primitive completion coefficients before the third prime is recovered from (4). The [reflected cubic edge bound](ABC_COHERENT_NORM_CUBIC_TRACE_COMPOSITION.md#counting-edges-with-a-prescribed-divisor) counts single edges for a prescribed index modulus. In (1) the divisor is detected jointly by the coefficient of the anchor and the traces of both completing edges.

## Counting two close Heron terms

**Proposition (Close Heron terms).** In the same prime and coefficient boxes, let \(0<\lambda\le1/2\). The number of coherent triangles over all fields for which two of
\[
                       A=P_1y^2,\qquad B=P_2z^2,\qquad C=P_3x^2
\]
differ by at most \(\lambda\) times the larger satisfies
\[
            \ll_\varepsilon (ZV)^\varepsilon
                           (V^2+Z)(1+\lambda Z^2).                 \tag{5}
\]
At \(V=Z^{5/4}\) and \(\lambda=Z^{-1}\), this is \(O_\varepsilon(V^{14/5+\varepsilon})\), including tuples with small lifts.

By relabeling, suppose the close pair is \(A,B\). Fix the signed coefficients \(y,z\). From \(1-\lambda\le A/B\le(1-\lambda)^{-1}\), the ratio \(p_1/p_2\) lies within \(O(\lambda)\) of \((|z|/|y|)^{2/3}\). Here the primes lie in one dyadic interval, so the ratios and their centers in any nonempty case lie in a fixed compact subinterval of the positive reals. The fractions \(p_1/p_2\) are reduced and distinct for distinct ordered pairs, since the primes are distinct. Their spacing is at least \(1/(4Z^2)\) by integral cross multiplication. There are \(O(1+\lambda Z^2)\) possible ordered pairs.

The integer \(A-B\) is nonzero: its two terms have opposite parity of valuation at \(p_1\). The reconstruction divisibility \(x\mid A-B\) therefore supplies subpower many signed choices for \(x\). For each one, use the unit cubic congruence (4) to count possible \(p_3\), and recover the entire triangle as in (3). The factor from the interval for \(p_3\) is \(O_\varepsilon((ZV)^\varepsilon(1+Z/|yz|))\). Summing over \(|y|,|z|\le V\) gives \(O(V^2+Z(1+\log V)^2)\) before the subpower adjustment. This proves (5). The argument counts the sector without imposing rank one or a lower bound on any lift.

## The local projective pattern of a lift

**Proposition (Projective criterion for a lift).** For an odd \(f>1\) dividing \(g_x\), one has \((f,d\delta p_1p_2p_3)=1\). In the ring group
\[
                  A_f=(\mathcal O_K/f\mathcal O_K)^\times/
                                     (\mathbb Z/f\mathbb Z)^\times,
\]
let \(h_f\) be the class of \(\sqrt{-d}\), a nonzero element of order two when \((f,d)=1\). Then the projective residues on \(12,23,13\) are \(0,h_f,h_f\), in additive notation. Conversely this residue pattern, for an odd \(f\) prime to \(d p_1p_2p_3\), holds precisely when \(f\mid g_x\). The other two cases follow by permutation.

If a prime of \(f\) divides \(x\), it cannot be \(p_1,p_2\). If it is \(p_3\) then the incident trace \(u_{23}\) is a unit at that prime by the norm and primitivity. Thus all three edge norms are units at \(f\). A prime dividing \(f\) and either \(d\) or \(\delta\) would give a nonunit norm on edge \(23\) from \(u_{23}\equiv0\bmod f\); thus none does. Now \(v_{12}\equiv0\bmod f\), but \(v_{23},v_{13}\) are units and the corresponding traces vanish. This gives precisely the three projective residues. Conversely \(\alpha_{12}\) is rational modulo \(f\) and the other two are unit multiples of \(\sqrt{-d}\). Since \(f\) is odd and prime to \(d\), these statements read \(f\mid v_{12}\), \(f\mid u_{23},u_{13}\), and \(v_{23},v_{13}\) are units by comparing coordinates. Hence \(\delta\) is a unit and \(f\mid x\).

For a fixed field of discriminant \(D<-4\) and prescribed maximal-order classes \(c_i\) in one coset of \(\mathrm{Cl}(K)[3]\), apply the [exact ring class residue criterion](ABC_COHERENT_NORM_RING_CLASS_RESIDUES.md#ideal-class-lifts-and-exact-residue-patterns). The order-two element \(h_f\) lies in \(3A_f\), so this pattern occurs at fixed \(f\) exactly when \(\theta_{3,f}(c_j-c_i)=0\) for every pair. In that case its limiting proportion among oriented prime triples in these classes is \(|3A_f|^{-2}\) as the prime interval grows. The same criterion and proportion hold for separately assigned, pairwise coprime odd divisors \(f_x,f_y,f_z\) of the three lifts, with \(f=f_x f_y f_z\), by imposing the three permuted patterns on their respective factors and using CRT. For an independent plane of class differences, each good prime of such a divisor that splits in the real quadratic resolvent must split in all four reflected cubic fields, by the [cubic test for the ring class residues](ABC_COHERENT_NORM_RING_CLASS_RESIDUES.md#cubic-fields-detect-the-first-color).

The fixed-conductor pattern and the moving-field bound (2) describe the large-lift sector. Formula (5) also controls close Heron terms. The aggregate rank-one count still includes the sector where all lifts are small and the Heron terms are separated. The retained bound on that sector is \(O_\varepsilon((ZV)^\varepsilon V^3)\).

## Exact replay near the crossover

The [signed witness checker](verify_coherent_norm_pencils.py) checks (3), (4), the anchor trace congruence and the prime units for the odd part of a lift on the retained examples. The [CRT census](census_coherent_norm_rank_one.cpp) gives a separate finite comparison in the positive odd boxes of the [rank-one census](ABC_COHERENT_NORM_PENCILS.md#census-at-the-rank-one-crossover). It computes the three gcds from the normalized coefficients and their trace formulas, for instance \(g_x=\gcd(x,(A-B)/x)\). It tests the threshold by the integer comparison \(\max(g_x,g_y,g_z)^4\ge V^3\).

| \(Z\) and coefficient box | Prime-unit positive rank-one records | Largest lift on those records | Lifts at least \(V^{3/4}\) | Rank-one records in the field window | Largest lift in the window |
| --- | ---: | ---: | ---: | ---: | ---: |
| \(128,\ [V,2V]^3\) | 29 | 23 | 0 | 3 | 1 |
| \(128,\ [1,2V]^3\) | 103 | 543 | 3 | 4 | 5 |
| \(192,\ [V,2V]^3\) | 116 | 35 | 0 | 5 | 7 |

The positive odd census finds no large-lift record in the field window at any retained scale. The wider box contains large lifts before the field test, and its window counts illustrate the need to study the small-lift sector as well. These are finite controls.

The same replay tests close Heron terms by the integer comparison \(Z|A-B|\le\max(A,B)\) and its two permutations. The \(Z=192\) annulus contains one such record among its 116 prime-unit positive rank-one records, before imposing the field window. This is the sole close record among the retained prime-unit positive rank-one boxes. None of the field-window records meets this closeness cutoff.

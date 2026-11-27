# Divisor bounds for the rank-two prime-path count

Use the full congruence moduli of the [quadratic 3-torsion prime-path proof](ABC_QUADRATIC_THREE_TORSION_FERMAT_TAIL.md#full-moduli-including-the-factor-2). Let \(p_1,p_2,p_3\in[Z,2Z)\) be the distinct odd primes of a selected path. The three signed path coefficients are \(O(V)\), and the full moduli \(q_1,q_2,r\) are pairwise coprime. Rank is the dimension spanned by the congruence-lattice vectors in the size-\(O(Z)\) ball used in that proof.

## Counting the rank-two congruences

**Proposition (Rank-two divisor bound).** For \(Z\ge2,V\ge1\), the rank-two records, counted with their prime triples, contribute
\[
 \ll_\varepsilon (ZV)^\varepsilon
    \left(Z^2V+\min\{Z^2V^{4/3},\,Z^3V^{1/2}\}\right).                 \tag{1}
\]
Here a record consists of the signed coefficients and the three local root classes of the full moduli. The estimate also bounds the corresponding selected paths when the field varies, by the exact recovery in the original proof.

Use the [cofactor construction and equation (6)](ABC_QUADRATIC_THREE_TORSION_FERMAT_TAIL.md#rank-two-and-the-primitive-coefficient-plane). Set \((t_1,t_2,t_3)=(q_1,q_2,r)\), and fix dyadic sizes \(Q_i\) for \(t_i\) and \(P\) for the product of the first two successive minima. Write \(Q_*=Q_1Q_2Q_3\), and let \(u\le v\le w\) be the three \(Q_i\) in increasing order. The lattice weight is \(O(Z^2/P)\), with \(P\ll Z^2\), and the three normal coordinates in the corresponding order are \(a_i t_i\). Once the moduli and the bounded signed factors \(b_i\) of the original coefficients are fixed, there are at most \((ZV)^\varepsilon\) local root-class choices, with constants depending on \(\varepsilon\). The trivial bound for the box, with its lattice weights, is therefore \(O_\varepsilon((ZV)^\varepsilon Z^2 Q_*/P)\).

The zero-coordinate case contributes \(O_\varepsilon((ZV)^\varepsilon Z^2w)\) as in the original proof. If the three coordinates are nonzero, \(P\gg w\) and the number of triples \(a_i\) is \(O(P^3/Q_*)\). The cofactor congruence takes the form
\[
        L_1t_1+L_2t_2+L_3t_3=m\,t_1t_2t_3,\qquad m\in\mathbb Z,         \tag{2}
\]
where each \(L_i\) is a nonzero fixed multiple of \(a_i^3\) once the \(b_i\) are fixed. In particular \(|m|\ll P^3/(Q_*u^2)\).

For \(m=0\), the primitive normal reduction and scaled-plane argument of the original proof give \(O_\varepsilon((ZV)^\varepsilon(P^3/Q_*+w^2))\) records in this box. That argument for the exact linear equation applies at every height \(P\). After weighting and taking the trivial bound as well, the contribution divided by \((ZV)^\varepsilon Z^2\) is bounded by a constant times
\[
  \min\{Q_*/P,\,P^2/Q_*\}+w^2/P \ \le\ Q_*^{1/3}+O(w)=O(V).
\]
The first minimum is at most the weighted geometric mean with exponents \(2/3,1/3\).

For \(m\ne0\), permute the indices so that \(t_3\) has size \(u\). With \(a_i,b_i,m,t_3\) fixed, (2) gives
\[
 (m t_3 t_1-L_2)(m t_3 t_2-L_1)=L_1L_2+m L_3t_3^2.                 \tag{3}
\]
When the right side is nonzero, factoring it gives at most subpower many choices of the two factors, each of which determines \(t_1,t_2\) uniquely. There are \(O(u)\) choices for \(t_3\) and \(O(P^3/(Q_*u^2))\) choices for the nonzero integer \(m\). After counting \(a_i\) and weighting, this bounds the contribution by
\[
 O_\varepsilon\left((ZV)^\varepsilon Z^2\frac{P^5}{Q_*^2u}\right).
\]
Combine it with the trivial bound by taking their weighted geometric mean with weights \(3/4,1/4\). The factor after \((ZV)^\varepsilon Z^2\) is at most a constant times
\[
        (Q_*/P)^{3/4}\left(\frac{P^5}{Q_*^2u}\right)^{1/4}
            = P^{1/2}(vw)^{1/4}\ll ZV^{1/2}.                        \tag{4}
\]

When the right side of (3) is zero, \(t_3^2\mid L_1L_2\) and \(m=-L_1L_2/(L_3t_3^2)\). Thus each \(a_i,b_i\) tuple admits at most subpower many \(t_3,m\) pairs in this case. One factor on the left vanishes, fixing one of \(t_1,t_2\), and the other has \(O(w)\) choices. After weighting, the resulting factor is \(P^2w/Q_*\). Its geometric mean with \(Q_*/P\) is \(P^{1/2}w^{1/2}\ll ZV^{1/2}\).

All integers used in divisor bounds have size bounded by a fixed power of \(ZV\), by \(P\ll Z^2\), \(Q_i\ge1\) and the bound on \(m\). Sum over the dyadic boxes and take the better of this bound and the original \(Z^2V^{4/3}\) bound to obtain (1). The lattice counting uses the same successive-minima principle as [Heath-Brown–Pierce, Section 4.3](https://arxiv.org/pdf/1409.3177). Equation (3) uses the additional relation of the three full moduli in one common cofactor.

## Occupied rank-one records with a short second direction

**Proposition (Short second direction).** Let \(\lambda_1,\lambda_2\) be the first two successive minima of the full congruence lattice of an occupied rank-one record. For \(T\ge1\), the number of such records with \(\lambda_1\lambda_2\le T\) is
\[
     \ll_\varepsilon (ZV)^\varepsilon
          \min\{V^3,\,V^2+T^{3/2}V^{1/2}\}.                            \tag{5}
\]
At \(V=Z^{5/4}\), the sector \(\lambda_2\le Z V^{1/30}\) contributes \(O_\varepsilon(V^{59/20+\varepsilon})\).

The shortest direction is the prime-triple line, since the record has rank one in the original ball and is occupied. Choose independent lattice vectors realizing \(\lambda_1,\lambda_2\) and take the cofactor of the entire intersection lattice with their plane. The intersection has the same first two minima. The cofactor congruence (2) still holds, since it uses the local ratios of any two lattice vectors. Thus the same counts apply with one record per lattice and no factor \(Z^2/P\). Here \(P\ll T\), and the dyadic sums cost at most a subpower loss: \(P\le\lambda_1\lambda_2\) up to dyadic constants, and the lattice of integer multiples of its index \(Q\) is contained in the full congruence lattice, so \(P\ll Q^2\ll V^6\).

For \(m=0\) and all cofactor coordinates nonzero, the previous bounds give \(O_\varepsilon((ZV)^\varepsilon(\min\{Q_*,P^3/Q_*\}+w^2))\), which is \(O_\varepsilon((ZV)^\varepsilon(P^{3/2}+V^2))\). With a zero cofactor coordinate the bound is \(O_\varepsilon((ZV)^\varepsilon P w)\), since two of the moduli have product \(O(P)\). In the two cases of (3) with \(m\ne0\), the geometric means used in (4) and the zero-product case, now without the lattice weight, bound the box by \(O_\varepsilon((ZV)^\varepsilon P^{3/2}V^{1/2})\). The term \(Pw\) is at most \(V^2\) for \(P\le V\), and at most \(P^{3/2}V^{1/2}\) for \(P\ge V\). Summing and including the trivial \(V^3\) bound gives (5).

A prime triple in the box has norm \(O(Z)\), so \(\lambda_1\ll Z\). At the stated crossover, \(T\ll Z^2V^{1/30}=V^{49/30}\) for the short second direction. Its term \(T^{3/2}V^{1/2}\) is \(O(V^{59/20})\). This sector bound is uniform under the original field-window and path-selection conditions. It leaves the rank-one records with longer second directions.

## Consequence for a future rank-one saving

In the original tail construction, up to \(X^\varepsilon\) allowances, \(Z=X^{1/2}/H\), \(V=Z^3/X^{1/2}\), and the per-field path weight is \(Z^3 H^2/X\). The two new rank-two terms of (1), after division by this weight, are
\[
         X^{3/2}H^{-4}+X^{3/2}H^{-7/2}.
\]
For \(H\ge1\) the second term dominates the first, and it crosses the first-moment tail \(X/H\) at \(H=X^{1/5}\). The original rank-two term \(X^{11/6}H^{-5}\) remains available for larger \(H\). The original rank-one term \(X^{5/2}H^{-8}\) crosses \(X/H\) at \(H=X^{3/14}\) and still controls the unconditional second-moment exponent \(17/14\).

If the occupied rank-one count saves a factor \(V^\sigma\) uniformly for \(X^{1/5}\le H\le C X^{1597/5000}\), with \(0<\sigma\le1/4\) and the original auxiliary allowances, then
\[
      \sum_{0<|d|\le X}^{*} h_3(d)^2
       \ll_\varepsilon X^{\,1+\max\{1/5,\,(3/2-\sigma)/(7-3\sigma)\}+\varepsilon}.
\]
Indeed the new rank-one tail is \(X^{5/2-\sigma}H^{-(8-3\sigma)}\). Multiply each dyadic tail by \(H^2\) and bound it by the first-moment tail until its crossover. The \(X/H^2\) part of the tail contributes \(O_\varepsilon(X^{1+\varepsilon})\). The same prime selection is available down to \(H=X^{1/5}\), where \(Z=X^{3/10+4\eta}\le X^{1/3}\) for sufficiently small \(\eta\). The reflected real fields use the same estimate by the fixed-scale comparison in the original proof.

At \(\sigma=1/4\) this would give \(6/5\). The existing [rank-one target](ABC_COHERENT_NORM_PENCILS.md#exact-examples-and-the-unresolved-count) asks for a saving only through \(\sigma=1/9\) and starting at \(H=X^{5/24}\). A larger saving and the lower threshold are both required for the \(6/5\) conclusion by this calculation.

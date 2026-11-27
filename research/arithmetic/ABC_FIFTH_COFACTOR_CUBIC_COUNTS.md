# Fifth-power paths with a cubic on the cofactor plane

The [fifth-power prime-path count](ABC_COHERENT_PRIME_PATH_HIGHER_ODD_TORSION.md#the-lattice-count-at-the-fifth-power) attaches a plane cubic to each cofactor normal. Its singular members reduce to linear equations with coefficients that are fifth powers. On the smooth members, the large coefficients of the primitive cubic improve the rational point count.

## Counting a box of cofactor normals

Use the complete moduli \(q_i\), their product \(Q\), and the local root classes of the original count. Fix dyadic sizes \(Q_i\) for \(q_i\), with \(Q_*=Q_1Q_2Q_3\), \(u=\min Q_i\), \(w=\max Q_i\), and \(v\) the remaining size. Write \(P\) for the dyadic magnitude of the integer cofactor normal. If all coordinates are nonzero, the original congruences give
\[
             n=(a_2q_2,a_3q_3,a_1q_1),\qquad
             a_1^5q_1^3+a_2^5q_2^3+a_3^5q_3^3=m q_1q_2q_3.          \tag{1}
\]
Here \(m\) is integral, \(P\gg w\), and each \(|a_i|Q_i\ll P\). Three moduli give at most \((ZV)^\varepsilon\) records, up to a constant depending on \(\varepsilon\), where the moduli are \(O(V)\) and the primes lie in \([Z,2Z)\). This trivial count in a box is \(O_\varepsilon((ZV)^\varepsilon Q_*)\). With a zero normal coordinate and a positive prime triple, two moduli divide both remaining normal components, so \(Q_*\ll P w\), as in the original argument.

### Singular fibers

Singularity of the cubic in (1) implies \(m^3=27(a_1a_2a_3)^5\) from the three partial derivative equations. Put \(t_i=\sqrt[3]{a_i^5}\,q_i\) using real cube roots. The identity
\[
  \sum t_i^3-3t_1t_2t_3
    =(t_1+t_2+t_3)(t_1^2+t_2^2+t_3^2-t_1t_2-t_2t_3-t_3t_1)
\]
shows that \(t_1+t_2+t_3=0\) or the three \(t_i\) are equal. Equality forces all \(a_i\) to have one sign. It is excluded here because the normal is perpendicular to a vector of positive primes.

For \(t_1+t_2+t_3=0\), both \((t_1/t_2)^3\) and \((t_1/t_2)(1+t_1/t_2)\) are rational, the second by expanding \((t_1+t_2)^3\). The ratio therefore has degree at most two and has rational cube. A real cube root of a rational number either is rational or has degree three, so the ratio is rational. The same holds cyclically. It follows that all ratios \(a_i/a_j\) are rational cubes, since the integers \(3\) and \(5\) are coprime. Write \(a_i=s b_i^3\) with \(b_i\) signed nonzero coprime as a triple and \(s\) a common nonzero integer. To see integrality of \(s\), clear the denominators of the ratios of the cube roots and remove the common content from the \(b_i\); their cubes then have gcd one. Thus
\[
                            b_1^5q_1+b_2^5q_2+b_3^5q_3=0.            \tag{2}
\]

Place \(R=\max_i |b_i|^3Q_i\) in a dyadic box, so \(w\ll R\ll P\). There are \(O(R/Q_*^{1/3})\) possible \(b\). For fixed \(b\), scale the \(q_i\) coordinates by \(Q_i^{-1}\). The lattice of (2) in the scaled plane has covolume
\[
       \frac{\|(b_1^5Q_1,b_2^5Q_2,b_3^5Q_3)\|}{Q_*}
                   \gg \frac{R^{5/3}}{Q_* w^{2/3}},
\]
since the unscaled normal is primitive. The scaled determinant formula is the one used in the [quadratic \(3\)-torsion plane count](ABC_QUADRATIC_THREE_TORSION_FERMAT_TAIL.md#rank-two-and-the-primitive-coefficient-plane). If the candidates span a line, primitivity of \((q_1,q_2,q_3)\) gives at most two. If they span the plane, the successive-minima bound in the scaled box is \(O(Q_*w^{2/3}/R^{5/3})\). Multiplying by the number of \(b\) and summing \(R\) gives
\[
          O_\varepsilon\big((ZV)^\varepsilon
                    (P/Q_*^{1/3}+Q_*^{2/3})\big).                  \tag{3}
\]
The factor \(s\) is unnecessary for recovering the moduli or the possible root classes, so it is absent from this count.

For example, the primes \((19,23,29)\) and normalized coefficients \((2,3,5)\) give a singular fiber on the plane with primitive normal \((3,-5,2)\). In (1) one has \(a=(1,1,-1)\) and \(m=-3\), and (2) is \(q_1+q_2-q_3=0\). The full lattice cuts out the whole integer lattice in this plane because its three congruences are reductions of \(3p_1-5p_2+2p_3=0\) modulo \(2,3,5\). The trace formulas give the positive defect \(4569140352316=4\cdot 89\cdot12834663911\), with the two factors after \(4\) distinct primes. In order \(12,23,13\) the trace and imaginary-coefficient pairs are \((6743258,4),(22066828,6),(9429848,10)\). These satisfy the fifth-power norms and the coherent product with middle factor \(23^5\).

### Smooth fibers

Let \(s=\gcd(a_1,a_2,a_3)>0\). For each \(s\) there are \(O(P^3/(Q_*s^3))\) possible \(a\). The primitive equation of the cubic in (1) has coefficient norm
\[
          \gg \big(\max_i |a_i|/s\big)^5\gg (P/(sw))^5.
\]
Indeed the content of the original equation divides \(s^5\), since it divides the gcd of the three pure coefficients. Apply [Walsh, Theorem 1.2, with degree three](sources/walsh-1308.0574v2.pdf) after dividing by the content. The cubic outside the singular case is smooth and hence geometrically irreducible. Since the moduli are positive and primitive, each rational projective point corresponds to at most one ordered integer triple. The count for fixed \(a,m\) is consequently
\[
        O_\varepsilon\big((ZV)^\varepsilon
                     (1+V^{2/3}(sw/P)^{5/9})\big).                 \tag{4}
\]
In the applications \(P\) and the coefficients of (1) are bounded by fixed powers of \(ZV\), so this absorbs the logarithm from Walsh's theorem. Summing the \(a\) count and (4) over \(s\) costs finite sums of \(s^{-3}\) and \(s^{-22/9}\). If \(m\) ranges over at most \(L\) candidates uniformly in \(a\), the smooth record count is at most
\[
    O_\varepsilon\left((ZV)^\varepsilon
          L\,\frac{P^3}{Q_*}\big(1+V^{2/3}(w/P)^{5/9}\big)\right).   \tag{5}
\]

## The rank-two tail term

**Proposition (Refined fifth-power rank-two count).** The number of rank-two coherent fifth-power paths counted in the complete-modulus lattices is
\[
 O_\varepsilon\left((ZV)^\varepsilon Z^2 V^{113/67}\right).
\]
When \(V\ll Z^5\), it is also \(O_\varepsilon((ZV)^\varepsilon Z^{134/45} V^{47/45})\). Both counts apply to either sign of the field.

Each rank-two record has weight \(O(Z^2/P)\) and \(P\ll Z^2\). With that weight, (3) costs \(O(Z^2 V)\), since \(P\gg w\) and \(Q_*^{2/3}/w\le w\). A zero coordinate also costs \(O(Z^2V)\). In the smooth case the bound \(|a_i|Q_i\ll P\) in (1) gives \(L\ll 1+P^5/(Q_*u^2)\). In (5), after putting aside the factor \((ZV)^\varepsilon Z^2\), this gives four terms
\[
 \frac{P^2}{Q_*},\qquad
 \frac{P^7}{Q_*^2u^2},\qquad
 V^{2/3}w^{5/9}\frac{P^{13/9}}{Q_*},\qquad
 V^{2/3}w^{5/9}\frac{P^{58/9}}{Q_*^2u^2}.                           \tag{6}
\]
Take the minimum with \(Q_*/P\) separately for each term. The first term contributes \(O(V)\) by a weighted geometric mean. The third is bounded by \(O(V^{2/3}P^2/Q_*)\) since \(P\gg w\), and its minimum contributes \(O(V^{11/9})\). For the second, use weights \(7/8,1/8\) with \(Q_*/P\), giving \(O(V^{13/8})\). For the fourth, weights \(58/67,9/67\) give
\[
   V^{6/67}w^{5/67}Q_*^{40/67}u^{-18/67}\ll V^{113/67}.
\]
Sum all dyadic boxes for the first assertion.

For the capped bound use weights \(4/5,1/5\) for the fourth term and \(Q_*/P\), cancelling the power of \(u\). The result is
\[
          P^{22/45}(vw)^{2/5}w^{1/9}V^{2/15}
               \ll Z^{44/45}V^{47/45}.
\]
The same weights give \(O(Z^{6/5}V^{4/5})\) for the second term. This is at most the capped target for \(V\ge Z^{10/11}\). Below that threshold the uniform \(V^{13/8}\) estimate is at most the target. The remaining contributions are bounded by it when \(V\ll Z^5\). This proves the capped assertion.

In the signed tail, put \(Z=X^{1/2+4\eta}/H\) and \(V=\max(1,Z^5/\sqrt X)\). Dividing the capped bound by \(Z^3 H^2/X\) gives the rank-two term \(X^{277/90}H^{-36/5}\), up to the original subpower allowances. If \(V=1\) the count is absorbed by \(Z^3\). The new term crosses \(X/H\) at \(H=X^{187/558}\), below the crossover \(X^{12/35}\) with the older rank-two term. The unconditional moment still uses the rank-one crossover \(H=X^{9/26}\).

## Small cofactor quotients in rank one

**Proposition (A rank-one saving at the fifth-power crossover).** Let the rank-one record be occupied by a coherent fifth-power path. Take the plane generated by a shortest lattice direction and a vector realizing the second successive minimum, and use the cofactor of the full lattice in that plane. At \(V=Z^{7/4}\), the records where this cofactor has a zero coordinate, or gives a singular fiber of (1), or gives \(|m|\le V^{1/3}\), number
\[
                         O_\varepsilon(V^{184/63+\varepsilon}).    \tag{7}
\]
The saving is uniform with any additional restrictions from field recovery or path selection.
In balanced boxes \(Q_i\asymp V\), the sector \(\lambda_2\le V^{5/7}\) also has a bound \(O_\varepsilon(V^{176/63+\varepsilon})\), with arbitrary quotient \(m\).

Write the three original minima as \(\lambda_1,\lambda_2,\lambda_3\). Occupation in rank one gives \(\lambda_1\asymp Z\), because the positive prime triple is primitive and lies in the shortest direction. The selected plane has first two minima \(\lambda_1,\lambda_2\). Its cofactor thus has magnitude \(P\asymp\lambda_1\lambda_2\ll\sqrt{ZQ_*}\), by the product of the three minima and \(\lambda_2\le\lambda_3\). This plane still satisfies the local normal relations, so (1) applies. There is only one positive primitive prime direction per record.

A zero coordinate costs \(O_\varepsilon((ZV)^\varepsilon P w)\ll V^{39/14+\varepsilon}\). Singular fibers contribute (3), at most \(O_\varepsilon(V^{2+\varepsilon})\) in this range. For the smooth fibers with \(|m|\le M\), apply (5) with \(L=O(1+M)\). Inserting \(P\ll\sqrt{ZQ_*}\) and \(w\ll V\) yields
\[
    O_\varepsilon\left((ZV)^\varepsilon(1+M)
       \big(Z^{3/2}V^{3/2}+Z^{11/9} V^{17/9}\big)\right).            \tag{8}
\]
At \(Z=V^{4/7}\) the exponents in the parentheses are \(33/14\) and \(163/63\), respectively. Taking \(M=V^{1/3}\) proves (7), which saves \(V^{5/63}\) up to the subpower factors.

For the balanced short-direction sector, \(\lambda_1\asymp Z\) and \(\lambda_2\le V^{5/7}\) give \(P\ll V^{9/7}\). Apply (5) with \(L\ll1+P^5/V^5\), using all possible quotients. The largest term is bounded by \(V^{11/9}P^{67/9}/V^8\ll V^{176/63}\), and the singular and zero-coordinate bounds are smaller in this range.

The remaining rank-one records have nonzero cofactor coordinates and a smooth cubic with \(|m|>V^{1/3}\), and in balanced boxes also \(\lambda_2>V^{5/7}\).

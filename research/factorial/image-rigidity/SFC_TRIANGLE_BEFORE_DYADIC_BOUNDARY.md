# The factorial triangle immediately before a dyadic boundary

For \(F=a xy+b yz+c zx\), the three moments starting at \(2^k-4\) have no common nonzero projective zero for any \(k\ge3\). The leading factorial coefficients force a common zero into a specific set of residues at two. Odd \(k\) has no such residue. For even \(k\), the exact recurrence forces a putative lift to agree modulo two with a root-of-unity configuration. A signed calculation modulo four excludes that configuration's lifts.

Use \(P_m=\Lambda(F^m)/(m!)^2\), \(s=a+b+c,\ h=ab+bc+ca,\ r=abc,\ p_j=a^j+b^j+c^j,\ \Delta=4h-s^2\), and \(W_j=s p_j-2p_{j+1}\). The [exact recurrence](SFC_TRIANGLE_OVERLAP_MOMENTS.md#exact-recurrence) is
\[
                \Delta P_m=c_m rP_{m-1}+W_{m+1},
    \qquad       c_m=\frac{2(2m+1)}m.                              \tag{1}
\]
The orders considered here are the three just before \(Q-1\), where \(Q=2^k\). In particular the window lies outside the scope of the [dyadic-boundary theorem](SFC_TRIANGLE_DYADIC_BOUNDARY.md#windows-containing-2k-1).

## Two leading factorial reductions

Write \(x_a=bc,\ x_b=ca,\ x_c=ab\), and use \(\sum_i\) to sum over the three indexed coefficients \(a_i=a,b,c\).

**Lemma (Leading coefficients just below the boundary).** For \(Q=2^k\), \(k\ge3\), both \(P_{Q-3}\) and \(2P_{Q-4}\) are integral at two. In characteristic two their reductions are
\[
\begin{split}
 \overline{P_{Q-3}}
    &=s^{Q-3}+\sum_{j=2}^{k-1}s^{Q-2^{j+1}+1}\sum_i x_i^{2^j-2},\
 \overline{2P_{Q-4}}
    &=r s^{Q-7}+\sum_{j=3}^{k-1}s^{Q-2^{j+1}+1}
                                      \sum_i a_i x_i^{2^j-3}.
\end{split}                                                       \tag{2}
\]
In these reduction formulas the right-hand sides are read over the residue ring.

For a monomial with exponents \(i_1,i_2,i_3\) of sum \(m\), the contribution at \(q=2^u\) to its coefficient valuation is the [factorial floor difference](SFC_TRIANGLE_DYADIC_BOUNDARY.md#the-backwards-recurrence-at-two)
\[
      \sum_{i<j}\left\lfloor\frac{\rho_i+\rho_j}{q}\right\rfloor
           -\left\lfloor\frac{\rho_1+\rho_2+\rho_3}{q}\right\rfloor ,
 \tag{3}
\]
where \(\rho_i\) is the residue of \(i_i\) between \(0\) and \(q-1\). First take \(m=Q-3\). At \(q=2\) the contribution is zero or two. At \(q\ge4\) with \(q\le Q\), a remainder total of \(q-3\) contributes zero. At total \(2q-3\) a pair fails to carry exactly when its opposite remainder is greater than \(q-3\). At most two remainders can be that large, so the contribution is nonnegative. It vanishes precisely for \((q-1,q-2,0)\) or \((q-2,q-2,1)\), up to permutation. The other possible total is \(3q-3\), with positive contribution.

If a term of valuation zero has a highest level \(q=2^j\) with a carry, its remainders there are those two patterns. They reduce to the required patterns at smaller levels. The next digit, at place \(q\), has no assignment since the total \(2q-3\) already equals the remainder of \(m\) modulo \(2q\). At each larger place with a one-bit, exactly one index receives that bit. The two top patterns together give \(s\sum_i x_i^{q-2}\): the terms with exponents \(q-1,q-2\) contribute \(x_i^{q-2}(s-a_i)\) for the complementary index, and the pattern with two copies of \(q-2\) contributes \(x_i^{q-2}a_i\). The terms without carries give \(s^{Q-3}\). Every valuation-zero coefficient is one modulo two.

For \(m=Q-4\), the \(q=2\) contribution is zero. At \(q=4\) the minimum is \(-1\), attained by the remainders \((2,1,1)\). At \(q\ge8\) the contributions are nonnegative. In fact a negative contribution would require total \(2q-4\) with all three remainders greater than \(q-4\), which is impossible. A term of valuation \(-1\) therefore has the minimum at \(q=4\) and zero at every other level. If it has a highest carry above \(q=4\), at that level \(q=2^j\) precisely two remainders must exceed \(q-4\). The residues modulo four allow only \((q-3,q-3,2)\) and \((q-3,q-2,1)\), up to permutation. They give \(s\sum_i a_i x_i^{q-3}\), by grouping the complementary index as in the first calculation. If there is no higher carry, the initial pattern gives \(rs\), the place \(4\) has no assignment, and the higher digits give \(s^{Q-8}\). Including the higher digit factors in each case gives the second equality in (2).

We also use the [established reduction of \(P_{Q-2}\)](SFC_TRIANGLE_DYADIC_BOUNDARY.md#the-factorial-coefficients-at-a-binary-boundary),
\[
       \overline{P_{Q-2}}
            =s^{Q-2}+\sum_{j=1}^{k-1}s^{Q-2^{j+1}}
                                             \sum_i x_i^{2^j-1}. \tag{4}
\]

## The residue configurations and their lifts

**Proposition (Residues of a putative lower-boundary zero).** Suppose \(P_{Q-4}=P_{Q-3}=P_{Q-2}=0\) with \(abc\ne0\). At every place above two of an algebraic coefficient field, integral projective normalization with at least one unit makes all three coefficients and \(s\) units. In the residue field, after setting \(s=1\), the coefficients must be a permutation of
\[
                     1+\beta,\quad 1+\omega\beta,\quad
                     1+\omega^2\beta,                            \tag{5}
\]
where \(\omega\) is a primitive cube root of unity, \(k\) is even, and
\(\beta^{Q/2+1}=1\), \(\beta^3\ne1\).
Every triple in (5) satisfies the three leading moment equations (2), (4) and the two reduced \(W\) equations. Moreover, in a normalized lift with \(s=3\), every coordinate agrees modulo two with the corresponding coordinate of the Teichmüller lift of (5), formed by taking \(1\) plus each lifted root of unity.

A complex common zero would give one with algebraic coefficients by the weak Nullstellensatz in the chart \(a=1\), with \(bc\) inverted. At such a point (1) gives \(W_{Q-2}=W_{Q-1}=0\). If \(\bar s=0\) and a coefficient is a nonunit, the two other residues must be equal nonzero elements, and (4) reduces to their product raised to \(Q/2-1\), which is nonzero. If all are units and \(\bar s=0\), the [dyadic valuation bound](SFC_TRIANGLE_DYADIC_BOUNDARY.md#the-backwards-recurrence-at-two) gives \(0<\nu(s)\le1/3\) for \(\nu(2)=1\). Divide the \(W\) equations by \(s\) and reduce to get \(p_{Q-2}=p_{Q-1}=0\) over the residue field. The two rows \((1,1,1)\) and \((a,b,c)\) have rank two there, and their kernel is generated by \((a,b,c)\) since \(s=0\). Normalizing the residues to \((1,y,1+y)\) gives \(y^{Q-3}=(1+y)^{Q-3}=1\). Comparing \((1+y)^Q=1+y^Q\) with \((1+y)^3\) yields \(y^2=y\), a contradiction. Thus \(s\) is a unit. A single unit coefficient would make \(W_{Q-2}\) nonzero modulo the maximal ideal. Two units with distinct residues would give two nonzero \(Q-2\) powers with sum zero and also vanishing sum of their \(Q-1\) powers, which forces the residues equal. Two equal unit residues make \(\bar s=0\) when the remaining coefficient is a nonunit. We have proved that all three coefficients are units. The residues are distinct because a repeated pair would make \(p_{Q-2}\) equal the nonzero power of the third residue.

Scale the actual point to \(s=3\), so \(s=1\) in the local ring modulo two. In the following computations work in that quotient ring, which has characteristic two. The two-row kernel argument now applies in the ring, since the minor from two distinct residues is a unit. Its generator is \((1+a_i)_i\), so \(a_i^{Q-2}=\lambda(1+a_i)\) with \(\lambda\) a unit. The cubic \(a_i^3=a_i^2+h a_i+r\) gives
\[
                 a_i^Q=\lambda(h a_i+r).
\]
Summing the three equations gives \(\lambda(h+r)=1\). Put \(u=h/(h+r)\), \(v=r/(h+r)\), so \(u+v=1\), and \(a_i^Q=u a_i+v\).

Define \(Y_i=\sum_{j=0}^{k-1}x_i^{2^j}\) and \(L_m=\sum_i a_i^m Y_i\). At \(s=1\), equations (4) and (2) give respectively \(L_1=0\), \(L_2=r\), by writing the sums in \(x_i\) in terms of \(Y_i\) and using \(x_i=r/a_i\). The Artin–Schreier relation for \(Y_i\) is
\[
            Y_i^2+Y_i=x_i^Q+x_i=v^2a_i^2+v a_i+v^2 h+v .
\]
Indeed \(x_i^Q=u^2x_i+uv(1+a_i)+v^2\) and \(x_i=a_i^2+a_i+h\). From \(L_1^2+L_2=r\) and \(p_2=p_4=1,\ p_3=1+h+r\), it follows that
\[
            r=\sum_i a_i^2(Y_i^2+Y_i)=v^2(1+h)+r .
\]
Since \(v\) is a unit, \(h=1\) modulo two. The elements \(\beta_i=1+a_i\) in that ring satisfy \(\beta_i^3=1+r\) and \(\beta_i^Q=\beta_i/(1+r)\), so \(\beta_i^{Q+2}=1\). They are units. In the residue field the three \(\beta_i\) differ by cube roots of unity, requiring \(3\mid Q+2\), hence even \(k\). Their powers then satisfy \(\beta_i^{Q/2+1}=1\) in the residue field. Conversely, for a triple (5), put \(\gamma=\beta^3\). One has \(s=h=1,\ r=1+\gamma\). The two power equations hold by \(a_i^Q=(a_i+r)/\gamma\), and \(Y_i=\beta_i^Q+\beta_i+k=(r/\gamma)\beta_i\) from \(x_i=\beta_i^2+\beta_i+1\). These give \(L_1=L_0=0,\ L_2=r\), hence (4) and the first equation of (2) vanish. The second equation of (2) at \(s=1\) is \((L_4+r p_3+r^2)/r^3\); this follows by inserting \(Y_i\) and removing its terms of exponents \(1,2,4\). Since \(a_i^4=(1+r)a_i+r\) when \(h=1\), \(L_4=0\), and \(p_3=r\), as required.

For the depth of a possible lift, take the odd-order Teichmüller lifts \(\xi_i\) of the three residue values of \(\beta_i\), and set \(A_i=1+\xi_i\). The \(\xi_i\) are \(\xi,\omega\xi,\omega^2\xi\) and satisfy \(\xi_i^{Q/2+1}=1\). Existence and uniqueness of each lift follow from the simple roots of this odd power equation. Let \(\delta_i=a_i-A_i\); then \(\sum\delta_i=0\). The ring relation \((a_i-1)^{Q+2}\equiv1\pmod2\) gives \(\nu((a_i-1)^{Q/2+1}-1)\ge1/2\), by squaring in characteristic two. The odd exponent and the common residue with \(\xi_i\) give \(\nu(\delta_i)\ge1/2\).

Suppose \(\tau=\min_i\nu(\delta_i)<1\). We have \(h(a)\equiv1\pmod2\) and \(h(A)=3\), so \(\nu(\sum A_i\delta_i)\ge1\), by expanding \(h(A+\delta)\) and using \(\sum\delta_i=0\) and \(\nu(\delta_i\delta_j)\ge1\). Divide the \(\delta_i\) by an entry of valuation \(\tau\) and reduce. The resulting nonzero vector satisfies the two rows \((1,1,1)\), \((\bar A_i)_i\), and is proportional to \((\bar\xi_i)_i\).

Write \(j=Q-2\). At the center \(A\), \(A_i^{Q/2-1}\equiv \xi_i^{-1}\pmod2\) since \(\xi_i^{Q/2}=\xi_i^{-1}\), and squaring gives \(A_i^j\equiv\xi_i^{-2}\pmod4\). Hence \(p_j(A)\equiv p_{j+1}(A)\equiv0\pmod4\) by the sums of the three root-of-unity powers, so \(W_j(A)\equiv0\pmod4\). In expanding \(W_j(A+\delta)\) with sum fixed to \(3\), its linear terms have valuation at least \(1+\tau\), and its terms of degree at least three have valuation at least \(3\tau\). The only possible part at valuation \(2\tau<2\) is \(3\binom j2\sum_i A_i^{j-2}\delta_i^2\). Here \(\binom j2\) is odd. In the leading deviation direction \((\bar\xi_i)_i\), the sum is a nonzero multiple of \(\sum_i\bar A_i^{-2}=(\bar h/\bar r)^2\), since \(A_i^j\equiv\xi_i^{-2}\). This contradicts \(W_j=0\) and proves \(\nu(\delta_i)\ge1\).

## The modulo-four obstruction

**Theorem (The window starting at \(2^k-4\)).** For every \(k\ge3\), the three factorial moments \(P_{2^k-4},P_{2^k-3},P_{2^k-2}\) have no common complex zero with \(abc\ne0\).

Odd \(k\) is excluded by the proposition. For even \(k\), take the normalized center \(A_i=1+\xi_i\) from that proposition, so \(s=h=3,\ r=1+\gamma\), \(\gamma=\xi_i^3\), and \(r\) is a unit. We prove
\[
                            P_{Q-2}(A)\equiv0\pmod4.              \tag{6}
\]
The calculation only needs the coefficients of valuations zero and one in \(P_{Q-2}\).

At \(q=2^u\), the remainder sum for \(m=Q-2\) is \(q-2\) or \(2q-2\). The latter has contribution \(2-\#\{\rho_i=q-1\}\) in (3). At valuation zero, the terms are either carry-free or have a highest carry level \(q=2^j\) with lower remainders \((q-1,q-1,0)\), as in (4). At valuation one there is one level with contribution one. It must be the highest carry level \(q=2^j\), \(j\ge2\), since a higher valuation-zero carry descends with contribution zero at each smaller level. At that level exactly one remainder is \(q-1\); the other two sum to \(q-1\) and reduce to \(q/2-1,0\) at \(q/2\) by the lower valuation-zero condition. Thus the complete pattern is \((q-1,q/2,q/2-1)\). Above the highest carry level the terms have the same independent high digit assignments used for (2).

For the units modulo four, the odd part of \(N!\) has sign given by the sum of \(\lfloor N/2^{u+2}\rfloor\) and the indicators for adjacent binary one-bits in positions \(u,u+1\), over \(u\ge0\). This follows by counting odd factors congruent to three modulo four at each binary level, as in the [dyadic-boundary sign calculation](SFC_TRIANGLE_DYADIC_BOUNDARY.md#the-factorial-coefficients-at-a-binary-boundary). The floor contributions cancel for a valuation-zero coefficient, since (3) at \(q=2\) is zero. For carry-free assignments all adjacent-bit contributions cancel, by comparing the three pair sums in the numerator with the individual indices and their total in the denominator. For a valuation-zero carry at \(q=2^j\) with complementary coordinate \(a_i\), the same cancellation holds inside the lower block and the upper carry-free block. Across the boundary between positions \(j,j+1\), only the total \(m\) and the paired sum of the two high lower remainders have bit \(j\) equal to one. The resulting sign is negative precisely if the first higher bit, at \(2q\), is assigned to \(a_i\). When \(j=k-1\) there is no higher bit and no sign correction.

All calculations that follow evaluate at \(A\), and \(x_i=A_j A_l=r/A_i\) for the complementary pair. One has \(p_{2^v}\equiv-1\pmod4\) for \(v\ge1\). Indeed \(p_1=3\), and \(p_{2d}=p_d^2-2\sum x_i^d\) with \(\sum x_i^{2^v}\equiv h^{2^v}=1\pmod2\). The pure carry-free block is therefore \((-1)^{k-1}=-1\). Before the sign corrections the valuation-zero carry contribution is
\[
       \sum_i x_i^{-1}\sum_{j=1}^{k-1}(-1)^{j+1}x_i^{2^j}.
\]
Here \(x_i=1-\xi_i+\xi_i^2\). For \(j\ge1\), squaring successively modulo four gives
\[
 x_i^{2^j}\equiv
  1+\xi_i^{2^j}+\xi_i^{2^{j+1}}
    +2(\xi_i^{2^{j-1}}+\xi_i^{2^j}+\gamma^{2^{j-1}})\pmod4 .
\]
The alternating sum telescopes to
\[
 1+\xi_i^2+\xi_i^{-2}+2(\xi_i+\xi_i^{-1}+Z),\qquad
                       Z=\sum_{j=0}^{k-2}\gamma^{2^j}.
\]
Using \(x_i^{-1}=(1+\xi_i)/r\) and the sums of the three root-of-unity powers gives an uncorrected carry contribution of \(3+2(1+Z)/r\) modulo four.

It remains to add the sign corrections and the contributions of valuation one. Each carries a factor two, so take their quotients by two in the residue ring. The sign corrections are
\[
          \sum_{j=1}^{k-2}\sum_i x_i^{2^j-1}A_i^{2^{j+1}}
                =\sum_{j=1}^{k-2} r^{2^j-1}p_{2^j+1}.
\]
The valuation-one pattern \((q-1,q/2,q/2-1)\) gives a contribution \(\sum_i A_i^{q-1}x_i^{q/2-1}(s-A_i)\) at a level \(q=2^j\). Higher digit factors are one modulo two. Summing over \(2\le j\le k-1\), using \(s=1\) and \(p_{2^{j-1}}=1\) modulo two, gives
\[
                   \sum_{i=1}^{k-2}r^{2^i-1}(1+p_{2^i+1}).
\]
The power sums cancel between these two corrections. Combining with the carry-free and uncorrected terms, the result modulo four is
\[
         2\left(1+\frac{1+Z+\sum_{j=1}^{k-2}r^{2^j}}r\right)=0.
\]
In the bracket the sum of \(1+Z+\sum_{j=1}^{k-2}r^{2^j}\) equals \(1+\gamma=r\), since \(r=1+\gamma\) and \(k-2\) is even. This proves (6).

At the center, \(p_{Q-1}\equiv0\pmod4\) and \(p_Q\equiv3\pmod4\) follow from \(A_i^{Q-2}\equiv\xi_i^{-2}\pmod4\). Hence \(W_{Q-1}(A)\equiv2\pmod4\). Equation (1) at \(m=Q-2\), with \(\Delta=3\) and \(c_{Q-2}\) a dyadic unit, now gives \(P_{Q-3}(A)\equiv2/r\pmod4\). On the affine plane \(s=3\), the mod-two formula in (2) for \(P_{Q-3}\) is a polynomial in squares. Its partial derivatives in the two free coefficient variables vanish modulo two. Changing the coordinates by \(\delta_i\in2\mathcal O\) with sum zero therefore preserves \(P_{Q-3}\) modulo four by polynomial expansion. This contradicts the putative zero.

## Range and exact checks

The family includes infinitely many starts outside the earlier six prime-power offsets. For \(k\equiv5\pmod{2189738460}\) and \(k>5\), the integers \(2^{k+1}-(7+2t)\) for \(0\le t\le5\) have valuations \(1,1,1,1,2,1\) at \(3,5,53,3,7,47\), respectively, and each exceeds the corresponding prime power. Indeed the valuations hold at \(k=5\), and the stated period is a common multiple of the Euler periods for the six next prime powers. These windows end at \(2^k-2\), so they contain no dyadic boundary \(2^j-1\). For \(k\ge5\) the first two dyadic prefix tests also do not apply: the top prefixes of lengths one and two start at \(Q/2\) and \(3Q/4\), respectively, with their permitted offsets too short to reach \(Q-4\) when \(k\) is large.

Coefficientwise versions of (2), (4) were checked internally in exact arithmetic for \(Q=8,16,32,64,128\). The root-of-unity center residues were also checked at even \(k=4,6,8,10,12\). The all-\(k\) conclusion follows from the written identities. The remaining quadratic-triangle windows stay in the [research DAG](../../MATH_STATUS.json).

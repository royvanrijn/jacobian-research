# Three earlier dyadic windows of the quadratic factorial triangle

For \(F=a xy+b yz+c zx\) with \(abc\ne0\), three consecutive factorial moments cannot all vanish at either start \(2^k-5\) or \(2^k-6\), for any \(k\ge3\). The two power-sum equations from the exact recurrence constrain the same residue triple as the leading factorial coefficients. The resulting finite-field trace equations exclude every residue configuration.

Write \(P_m=\Lambda(F^m)/(m!)^2\), \(Q=2^k\), \(s=a+b+c\), \(h=ab+bc+ca\), \(r=abc\), \(a_i=a,b,c\), \(x_i=r/a_i\) for the product of the two other coefficients, and \(p_j=\sum_i a_i^j\). A three-moment zero at a start \(n\) implies \(W_{n+2}=W_{n+3}=0\) with \(W_j=s p_j-2p_{j+1}\), by the [exact recurrence](SFC_TRIANGLE_OVERLAP_MOMENTS.md#exact-recurrence). A complex common zero with \(abc\ne0\) gives one with algebraic coefficients by using the chart \(a=1\), inverting \(bc\), and applying the weak Nullstellensatz. At a place above two use \(\nu(2)=1\) and normalize the coefficients to be integral with at least one unit.

## The next two leading factorial layers

**Lemma (The layers at \(Q-5\) and \(Q-6\)).** For \(k\ge3\), the polynomial \(P_{Q-5}\) has coefficients integral at two. For \(k\ge4\) the same is true of \(2P_{Q-6}\). Their reductions are
\[
\begin{split}
 \overline{P_{Q-5}}&=s^{Q-5}
     +\sum_{j=3}^{k-1}s^{Q-2^{j+1}+3}\sum_i x_i^{2^j-4},\
 \overline{2P_{Q-6}}&=r^3s^{Q-15}
     +\sum_{j=4}^{k-1}s^{Q-2^{j+1}+1}\sum_i a_i^3 x_i^{2^j-5}.
\end{split}                                                       \tag{1}
\]
On the right-hand sides all operations are in characteristic two.

For monomial exponents \(i_1,i_2,i_3\) of sum \(m\), the contribution to the factorial coefficient valuation at \(q=2^j\) is the sum of the three pair carry indicators minus \(\lfloor(\rho_1+\rho_2+\rho_3)/q\rfloor\), for \(\rho_i=i_i\bmod q\) in \([0,q)\). This is the [factorial floor difference](SFC_TRIANGLE_DYADIC_BOUNDARY.md#the-backwards-recurrence-at-two). At a remainder sum \(q-d\) there are no carries. At sum \(2q-d\) a pair fails to carry precisely when its opposite remainder exceeds \(q-d\). At sum \(3q-d\ge2q\) all pairs carry and the contribution is positive.

For \(m=Q-5\), at \(q=2,4\) the level contributions are zero without a carry and positive otherwise. At \(q\ge8\) a negative contribution would require three remainders at least \(q-4\) summing to \(2q-5\), which is impossible. Thus a unit coefficient having a highest carry at \(q\ge8\) requires two remainders above \(q-5\). Write those two remainders \(q-A,q-B\), with \(1\le A\le B\le4\); the other is \(A+B-5\). The allowed pairs \((A,B)\) are \((1,4),(2,3),(2,4),(3,4),(4,4)\). The sum condition forces \(A+B\ge5\), and \((3,3)\) fails at \(q=2\) by having three odd entries. Each listed pair has zero contribution at every lower level. With \(a_i\) the variable of the complementary coordinate and \(b_i,c_i\) the remaining variables, the sum of the patterns divided by \(x_i^{q-4}\) is
\[
     b_i^3+c_i^3+b_i^2c_i+b_ic_i^2
          +a_i(b_i^2+c_i^2)+a_i^2(b_i+c_i)+a_i^3=s^3 .
\]
Above the highest carry, the bit at \(q\) is absent and each bit of \(m\) at a higher position is assigned without carries. This supplies the remaining factor \(s^{Q-2q}\). The terms with no carries give \(s^{Q-5}\), proving the first formula.

For \(m=Q-6\), the \(q=2\) contribution is always zero. At \(q=4\) the contribution is nonnegative, since the only possible carry total is six and at most two remainders equal three. At \(q=8\) the minimum contribution \(-1\) occurs at \((4,3,3)\) and its permutations, because the total must be ten and every entry greater than two. At \(q\ge16\) a negative contribution would require three remainders at least \(q-5\) summing to \(2q-6\), which is impossible. The terms of valuation \(-1\) therefore use the pattern at eight and have contribution zero at all other levels. If the highest carry is \(q\ge16\), there must be two remainders above \(q-6\). Their residues at eight allow only \((q-5,q-5,4)\) and \((q-4,q-5,3)\) up to permutation. These both satisfy the lower-level conditions and, grouped by the complementary coordinate, give \(s\sum_i a_i^3 x_i^{q-5}\). The higher carry-free bits give the remaining factor \(s^{Q-2q}\). If there is no carry above eight, the pattern \((4,3,3)\) gives \(r^3s\) followed by the factor \(s^{Q-16}\) from the higher bits. Each surviving coefficient of the scaled polynomial is one modulo two. This proves the second formula.

## Unit normalization and trace equations

Suppose there is a common zero at \(n=Q-5\) or \(n=Q-6\). We may assume \(k\ge4\) for the second case, since \(k=3\) there is covered by the [window containing \(3\)](SFC_TRIANGLE_DYADIC_BOUNDARY.md#windows-containing-2k-1). We show that \(s\) and all three coefficients are units at the chosen place.

Suppose first that \(\bar s=0\). If a coefficient is a nonunit, the other two are units with equal residues. The [two-unit valuation test](SFC_TRIANGLE_DYADIC_BOUNDARY.md#the-reduction-with-two-unit-coefficients) gives \(0<\nu(s)\le1/(n+2)\). In fact the nonunit coefficient has this same valuation because the difference of the two units has larger valuation. If all three coefficients are units, the [all-order valuation bound](SFC_TRIANGLE_DYADIC_BOUNDARY.md#the-backwards-recurrence-at-two) gives \(0<\nu(s)\le1/3\) for \(n=Q-5\) and \(0<\nu(s)\le1/4\) for \(n=Q-6\). Dividing the exact \(W\) equations by \(s\) in the three-unit case gives \(p_{n+2}=p_{n+3}=0\) in the residue field. The nullspace of the two rows \((1,1,1),(a,b,c)\) in that field is generated by \((a,b,c)\), since their sum and the sum of their squares are zero. After scaling one residue to \(1\), the other two are \(y,1+y\) with \(y\ne0,1\) and \(y^{n+1}=(1+y)^{n+1}=1\).

At \(n=Q-5\) use the [known leading layer](SFC_TRIANGLE_BEFORE_DYADIC_BOUNDARY.md#two-leading-factorial-reductions)
\[
 \overline{P_{Q-3}}=s^{Q-3}+
                 \sum_{j=2}^{k-1}s^{Q-2^{j+1}+1}\sum_i x_i^{2^j-2}.
\]
Its last sum \(s\sum_i x_i^{Q/2-2}\) has valuation exactly \(\nu(s)<1\). In the two-unit case only the product of those two units contributes to the inner residue sum. In the three-unit case all \(a_i^{Q-4}\) agree and are nonzero by the nullspace condition. Taking square roots in the residue field shows that the three \(x_i^{Q/2-2}\) agree and are nonzero. Their sum is nonzero in characteristic two. All the other terms have a larger power of \(s\), and the difference from the displayed reduction has valuation at least one. This contradicts \(P_{Q-3}=0\).

At \(n=Q-6\), the first reduction of (1) has last term \(s^3\sum_i x_i^{Q/2-4}\). In the two-unit case the inner residue sum is again nonzero. In the three-unit case \(y^{Q-5}=(1+y)^{Q-5}=1\). Comparing \((1+y)^Q=1+y^Q\) with \((1+y)^5=1+y+y^4+y^5\) gives \(y^4=y\). Thus \(y,1+y\) are primitive cube roots of unity and \(Q\equiv5\pmod3\). This makes \(Q/2-4\) a multiple of three, so the three product powers in the inner residue sum each equal one after the scaling. The last term then has valuation \(3\nu(s)<1\) in both unit cases. All other terms have larger valuation, a contradiction to \(P_{Q-5}=0\).

We have proved that \(s\) is a unit. Reduce at two and scale \(s=1\). A single nonzero coefficient residue gives \(p_{n+2}\ne0\). There cannot be just two nonzero residues, since \(p_{n+2}=p_{n+3}=0\) forces them equal and then \(s=0\). A repeated pair among three nonzero residues cancels and leaves a nonzero term in \(p_{n+2}\). Thus the three residues are nonzero and distinct. In the rest of the proof all variables denote their values in the residue field. The two rows \((1,1,1),(a,b,c)\) now have kernel generated by the vector with entries \(1+a_i\). Hence
\[
                        a_i^{n+2}=\lambda(1+a_i)
                               \quad\text{for some }\lambda\ne0.     \tag{2}
\]
The cubic relation for each \(a_i\) is \(a_i^3=a_i^2+h a_i+r\).

For each \(x_i\) define \(Y_i=\sum_{j=0}^{k-1}x_i^{2^j}\) and \(L_v=\sum_i a_i^vY_i\). They satisfy \(Y_i^2+Y_i=x_i^Q+x_i\). Since \(x_i=r/a_i\), the known layers for \(P_{Q-3}\) and \(2P_{Q-4}\) give respectively
\[
          P_{Q-3}=0\ \Longrightarrow\ L_2=r,\qquad
          P_{Q-4}=0\ \Longrightarrow\ L_4=r(1+h).                    \tag{3}
\]
For the first implication, multiply the reduction by \(r^2\) and replace \(a_i^2x_i^{2^j}\) by \(r^2 x_i^{2^j-2}\). The omitted terms at \(j=0,1\) contribute \(r p_1+r^2p_0=r+r^2\), where \(p_0=1\) in characteristic two. For the second implication, multiply the reduction of \(2P_{Q-4}\) by \(r^3\) and account for \(j=0,1,2\). The omitted terms give \(r p_3+r^2 p_2+r^4 p_0\). Since \(p_3=1+h+r\), \(p_2=p_0=1\), this yields (3). Likewise the second identity of (1) gives
\[
          P_{Q-6}=0\ \Longrightarrow\ L_8=r p_7+r^2p_6+r^4.         \tag{4}
\]
Indeed multiply by \(r^5\) and the omitted terms at \(j=0,1,2,3\) give \(r p_7+r^2p_6+r^4p_4+r^8p_0\), where \(p_4=1\). The pure term of the reduction cancels \(r^8p_0\).

## The start at \(Q-5\)

**Theorem (The three-moment window at \(2^k-5\)).** For every \(k\ge3\), the moments \(P_{2^k-5},P_{2^k-4},P_{2^k-3}\) have no common complex zero with \(abc\ne0\).

Multiplying (2) by \(a_i^3\) gives \(a_i^Q=\lambda(h a_i^2+r a_i)\) by the cubic relation. Summation gives \(\lambda(h+r)=1\). Write \(u=h/(h+r),v=r/(h+r)\) so \(u+v=1\) and \(v\ne0\). Multiplying the two corresponding expressions for the complementary coordinates gives \(x_i^Q+x_i=u^2x_i^2+uvr+u x_i\), using \(x_i(1+a_i)=x_i+r\). Equations (3) and \(Y_i^2+Y_i=x_i^Q+x_i\) give
\[
       r^2+r(1+h)=L_2^2+L_4
           =u^2 r^2+uvr+ur(1+h+r).
\]
After cancellation, the difference of the two sides, divided by \(r\), is \(v^2(1+r)\), using \(v h=u r\). Thus \(r=1\).

Take the product of \(a_i^Q=a_i(h a_i+1)/(h+1)\). It gives \((h+1)^3=\prod_i(1+h a_i)=1+h\), hence \(h=0\). Now \(a_i^Q=a_i\), and each residue is a root of the irreducible polynomial \(T^3+T^2+1\) over \(\mathbb F_2\). Therefore \(3\mid k\). The three \(x_i=1/a_i\) are conjugate with sum \(h=0\) in \(\mathbb F_8\), so every block of three terms in each \(Y_i\) sums to zero. This gives \(L_2=0\), contrary to \(L_2=r=1\).

## The start at \(Q-6\)

**Theorem (The three-moment window at \(2^k-6\)).** For every \(k\ge3\), the moments \(P_{2^k-6},P_{2^k-5},P_{2^k-4}\) have no common complex zero with \(abc\ne0\).

The case \(k=3\) was already reduced to the window containing three. For \(k\ge4\), multiply (2) by \(a_i^4\), reduce powers by the cubic, and sum over \(i\). This gives
\[
 D_0 a_i^Q=U a_i^2+V a_i+W,\qquad
 U=h+r,\quad V=h^2,\quad W=hr,\quad D_0=(1+h)(h+r).                 \tag{5}
\]
In particular \(D_0\ne0\), since \(\lambda D_0=1\). The complementary products satisfy
\[
 D_0^2 (x_i^Q+x_i)=U^2 x_i^2+(UV+V^2+D_0^2)x_i+UVr
                +UW(1+a_i^2)+VW(1+a_i)+W^2 .                       \tag{6}
\]
This follows by multiplying the two values of (5) for the indices other than \(i\), using \(x_i(1+a_i)=x_i+r\).

By (3), (4) and \(Y_i^2+Y_i=x_i^Q+x_i\), the sum of \(a_i^8\) times (6) equals \(D_0^2(r^2(1+h)^2+r p_7+r^2 p_6+r^4)\). Subtracting this expression from the sum and reducing by \(a_i^3=a_i^2+h a_i+r\) gives
\[
                             r^4(h^2+r^2+h+1)=0.                    \tag{7}
\]
For the reduction one may use \(p_0=p_1=p_2=p_4=p_8=1\), \(p_3=h+r+1\), \(p_5=h^2+hr+h+r+1\), \(p_6=h^2+r^2+1\), and \(p_7=h^3+h^2r+hr+r^2+h+r+1\). Equivalently reduce the polynomial \(a_i^8\) times the right side of (6) to degree at most two and sum its three coefficients. All these reductions use the same monic cubic.

Let \(t=h+r\). Since \(r\ne0\), (7) gives \(h=t^2+1,\ r=t^2+t+1\) and \(D_0=t^3\). Taking the second and third elementary symmetric functions of (5) yields
\[
        h^Q=t^4+1,\qquad r^Q=(t^8+t^4+1)/t^8 .
\]
To check this calculation, the three images in (5) are the eigenvalues of multiplication by \((W+VA+UA^2)/D_0\) on the quotient by \(f(A)=A^3+A^2+hA+r\). In the basis \(1,A,A^2\), the unscaled columns are
\[
 \begin{pmatrix} W\V\U\end{pmatrix},\qquad
 \begin{pmatrix} Ur\W+Uh\V+U\end{pmatrix},\qquad
 \begin{pmatrix} (V+U)r\Vh+U(h+r)\W+V+U(1+h)\end{pmatrix}.
\]
Substituting \(h=t^2+1,\ r=t^2+t+1\) into its characteristic polynomial after scaling gives the stated coefficients. Now \(h^Q=(t^Q)^2+1\) implies \(t^{24}+t^{20}+t^8+1=0\). This is \(((t+1)(t^2+t+1)(t^3+t^2+1))^4=0\). The middle factor would give \(r=0\).

If \(t=1\), then \(h=0,r=1\), and (5) gives \(a_i^Q=a_i^2\) on the roots of \(T^3+T^2+1\) in \(\mathbb F_8\). Thus \(k\equiv1\bmod3\). The \(x_i\) have trace zero as above, hence \(Y_i=x_i\) and \(L_4=p_3=0\). This contradicts \(L_4=r(1+h)=1\).

In the other case \(t^3+t^2+1=0\) and \(t\in\mathbb F_8\setminus\mathbb F_2\). From \(h^Q=t^4+1\) one gets \(k\equiv1\bmod3\). The following computations are in \(\mathbb F_8[A]/(A^3+A^2+(t^2+1)A+t^2+t+1)\), with \(x=A^2+A+t^2+1\). Repeated squaring and reduction give
\[
\begin{array}{c|l}
  j & A^{2^j}\\hline
  1 & A^2\
  4 & tA^2+(t^2+1)A+t^2+t\
  7 & (t+1)A^2+(t^2+1)A+t^2+t+1\
  9 & A
\end{array}
 \qquad
 \begin{split}
 \sum_{j=0}^{8} x^{2^j}&=0,\
 A^4\sum_{j=0}^{6}x^{2^j}&=A^2+(t^2+1)A+t+1 .
 \end{split}                                                       \tag{8}
\]
Here each power is obtained from the preceding power by squaring and substituting \(A^3=A^2+(t^2+1)A+t^2+t+1\) and \(t^3=t^2+1\). The right-hand side of (5) divided by \(D_0\) equals the row \(j=7\). Among \(k\equiv1\bmod3\), (8) forces \(k\equiv7\bmod9\). Indeed the values for the other two exponents differ as degree-at-most-two polynomials, and could not agree on the three distinct \(a_i\). The nine-term sum in (8) makes \(Y_i\) equal to its first seven terms for this residue class of \(k\). Summing the last displayed expression at all \(a_i\) gives \(L_4=t^2+t+1\), since \(p_2=p_1=p_0=1\). But \(r(1+h)=t^2+t\). This contradicts (3) and proves the theorem.

## The start at \(Q-7\)

**Theorem (The three-moment window at \(2^k-7\)).** For every \(k\ge3\), the moments \(P_{2^k-7},P_{2^k-6},P_{2^k-5}\) have no common complex zero with \(abc\ne0\).

The initial window \(k=3\) is already known. Take \(k\ge4\) and suppose there is a three-moment zero. If \(\bar s=0\) with two unit coefficients, the two-unit valuation test gives \(0<\nu(s)\le1/(Q-5)\). In the first reduction of (1), the last term \(s^3\sum_i x_i^{Q/2-4}\) then has valuation \(3\nu(s)<1\), since the product of the two units is the only nonzero residue product. This is impossible at a zero of \(P_{Q-5}\). If all three coefficients are units with sum a nonunit, the all-order valuation bound again allows division of the exact \(W\) equations by \(s\). Their residue nullspace forces \(y^{Q-6}=(1+y)^{Q-6}=1\) for the normalized residues \(1,y,1+y\). Taking square roots gives \(y^{Q/2}=y^3\) and \((1+y)^{Q/2}=(1+y)^3\). The Frobenius identity now gives \(y+y^2=0\), a contradiction. Thus \(s\) is a unit, and the same argument as before gives three distinct unit residues and (2), with exponent \(n+2=Q-5\).

In this case \(P_{Q-5}=0\) gives \(L_4=r(1+h)\) by the first reduction of (1). Multiplication by \(r^4\) and insertion of the omitted exponents \(j=0,1,2\) give \(L_4=r^4+r p_3+r^2p_2+r^4p_0=r(1+h)\). The moment \(P_{Q-6}=0\) gives (4) as before.

Multiply (2) by \(a_i^5\), reduce by the cubic, and sum. The result again has the form \(D_0 a_i^Q=U a_i^2+V a_i+W\), now with
\[
           U=h^2+h+r,\qquad V=h^2,\qquad W=hr+r^2,\qquad
           D_0=hr+r^2+h+r .
\]
The sum of the powers gives \(\lambda D_0=1\), so \(D_0\ne0\). Formula (6) follows for these values by the same multiplication of the complementary coordinates. Use \(L_4^2+L_8=\sum_i a_i^8(x_i^Q+x_i)\) and the values of \(L_4,L_8\) just obtained. Reducing the right side of (6) times \(a_i^8\) and summing, exactly as in (7), now gives
\[
                              r^5(r+1)=0 .
\]
Since \(r\ne0\), this forces \(r=1\). Substitution into \(D_0\) gives zero, contrary to the nonzero sum multiplier.

## Range and exact checks

The two-coordinate boundary still has the [two-moment conclusion](SFC_TRIANGLE_OVERLAP_MOMENTS.md#consequences-in-characteristic-zero). The new starts include infinitely many outside the earlier prime-power offsets. For \(k\equiv24\pmod{16895838960}\), the eight odd integers \(2^{k+1}-j\) with \(j=9,11,13,15,17,19,21,23\) have valuations \(1,2,1,1,1,1,1,1\) at \(7,3,197,13,3,17,11,3\), respectively. Each valuation follows by evaluating at \(k=24\) and using that the period is a common multiple of the Euler periods at the next prime powers. Each integer is larger than the indicated prime power. The six tests for each start are among these eight. For these \(k\), the three-moment windows contain no order of the form \(2^\ell-1\), and the bands above the first two dyadic prefixes end before \(Q-7\).

Equation (1) was checked coefficientwise in exact arithmetic for \(Q=8,16,32,64,128,256\), omitting \(Q=8\) in its second line. The [exact checker](verify_sfc_triangle_earlier_dyadic_windows.py) also replays the eliminations, finite-field table and progression. Replay from the repository root with Sage 10.9.

~~~sh
sage -python research/factorial/image-rigidity/verify_sfc_triangle_earlier_dyadic_windows.py
~~~

The arguments apply to every \(k\) in their stated ranges. The remaining starts are the open quadratic-triangle problem.

# A two-scale least-repair gap ratio recovers primitive pairs

Let \(n=X^L\), with \(L>0\) fixed before the height limit. Rows \(i=1,2\) are distinct positive primitive triples \(a_i+b_i=c_i\le X\). Keep every original radical cap and any additional least-repair or numerical cutoff.

Fix \(r_->\xi>0\) and \(K>0\) with \(K(r_--\xi)>2\), and set \(M=\lceil K/L\rceil\). Suppose
\[
 n^{r_-}\le \operatorname{rad}(c_i)\le n^{r_+},
 \qquad \operatorname{rad}(\gcd(a_1 b_1 c_1,a_2 b_2 c_2))\le n^\xi.
\]
Write \(Z_m(c)=\prod_{p^e\parallel c}p^{(m-e)_+}\) for the least multiplier making \(c\) \(m\)-full. For \(m=M,M+1,M+2\) put
\[
 U_m=a_1Z_m(c_1)-a_2Z_m(c_2),\quad
 V_m=b_1Z_m(c_1)-b_2Z_m(c_2),\quad
 W_m=U_m+V_m.
\]
For \(m=M,M+1\), put
\[
 \lambda_m=U_{m+1}/U_m-V_{m+1}/V_m,\qquad
 \mathcal R=\lambda_{M+1}/\lambda_M .
\]
All displayed ratios are well-defined and nonzero where required. Assume the shell excludes the exceptional primitive row \((1,1,2)\).

For each exact nonzero rational \(\mathcal R\), the number of ordered row pairs giving that record is, uniformly,
\[
                 O_{r_-,r_+,\xi,K,L,\varepsilon}
                        (n^{\,\xi+2/K}X^\varepsilon).
\]
The theorem needs no assumptions on degree measures, dictionary prices, prescribed shifts, or cheap private-factor extraction.

## Private prime support from the ratio

Let \(Q_{i,m}=Z_{m+1}(c_i)/Z_m(c_i)=\prod_{p^e\parallel c_i,\ e\le m}p\). The private part \(D_i\) of \(Z_M(c_i)\), supported off the entire other row, obeys
\[
 \log_X D_i\ge M L(r_--\xi)-\log_X c_i>1.
\]
Indeed on every private \(p^e\parallel c_i\) use \((M-e)_+\ge M-e\) and sum; removing the shared radical loses at most \(\xi\) in base \(n\). Thus \(D_i>X\) and in particular contains a private deficient prime. At every prime \(p\mid D_1\), write \(e=v_p(c_1)<M\). Then \(U_m,V_m,W_m\) are units for \(m=M,M+1,M+2\): the row-1 contributions vanish modulo \(p\), while the row-2 contributions are units because the prime is absent from that row. Similarly for \(D_2\). In particular the complete differences at all these scales are nonzero.

Primitivity and distinctness give \(h=a_2c_1-a_1c_2\ne0\). The cross identity gives
\[
 \lambda_m=\frac{h Z_m(c_1)Z_m(c_2)(Q_{1,m}-Q_{2,m})}{U_mV_m}
          \quad(m=M,M+1).
\]
At every \(p\mid D_1\), the cross difference \(h\) and the factor \(Q_{1,m}-Q_{2,m}\) are units by primitivity and \(p\nmid c_2\). The same holds at each private deficient \(p\mid D_2\). In particular both displayed gaps are nonzero and at every such prime
\[
 v_p(\lambda_m)=m-v_p(c_i),\qquad v_p(\mathcal R)=1.
\]
Thus every private deficient prime of either target is in the reduced numerator of \(\mathcal R\), to valuation exactly one there.

The shared squarefree radical has at most \(n^\xi\) possible values. Each overfull radical, supported on \(v_p(c_i)\ge M\), has size at most \(X^{1/M}\le n^{1/K}\) and thus at most \(n^{1/K}\) possible values. After these three choices, the support of both targets is contained in their union with that of the fixed numerator of \(\mathcal R\). All these integers have bounded polynomial height in \(X\): for example \(Z_m(c_i)\le\operatorname{rad}(c_i)^m\) and \(Q_{i,m}\le \operatorname{rad}(c_i)\) at the finite scales used. Fixed-support factor enumeration therefore gives \(X^{o_X(1)}\) choices of \(c_1,c_2\), uniformly on the recorded data.

## The fixed-target Pell terminal

For fixed targets write \(x=U_M,\ y=U_{M+1}\), \(w=W_M,\ w'=W_{M+1}\), and put
\[
 T=\frac{Z_{M+1}(c_1)Z_{M+1}(c_2)(Q_{1,M+1}-Q_{2,M+1})}
         {Z_M(c_1)Z_M(c_2)(Q_{1,M}-Q_{2,M})}.
\]
The previous identity cancels the factor \(h\) from the recorded ratio:
\[
       y(w'-y)=\kappa x(w-x),\qquad \kappa=T/\mathcal R\ne0 .
\]
The map from \((a_1,a_2)\) to \((x,y)\) is invertible, since its two linear forms have determinant up to sign \(Z_M(c_1)Z_M(c_2)(Q_{1,M}-Q_{2,M})\ne0\). Completing squares gives
\[
 (2y-w')^2-\kappa(2x-w)^2=w'^2-\kappa w^2.
\]
If the right-hand side is nonzero, this has \(O_\varepsilon(X^\varepsilon)\) solutions under the fixed polynomial bounds. Indeed write \(\kappa=p/q\) reduced and clear denominators to obtain
\[
       (q(2y-w'))^2-pq(2x-w)^2=q(qw'^2-pw^2)\ne0.
\]
If \(pq\) is a square use divisor enumeration. Otherwise remove its square part and use the corresponding quadratic field norm. An algebraic integer of the fixed norm generates an integral ideal of that norm with at most \(X^{o_X(1)}\) possible factorizations. For each principal ideal the generators differ by units. Imaginary quadratic fields have uniformly bounded unit group. For a real quadratic field the two conjugates have polynomial upper bounds and polynomial lower bounds from their fixed nonzero norm. The possible powers of a fundamental unit therefore lie in an interval of \(O(\log X)\) length, since fundamental units above one at a real embedding are uniformly bounded away from one by the integral trace and norm. This bound is uniform for all moving polynomial-height coefficients. The map back to \((a_1,a_2)\) is injective, and \(b_i=c_i-a_i\) gives the original rows.

If the right-hand side vanishes, then \(w,w'\ne0\) and \(\kappa=(w'/w)^2\). The plus sign gives \(y/w'=x/w\). Multiplication by \(w w'\) in their definitions forces \(a_1c_2-a_2c_1=0\), since \(Q_{1,M}\ne Q_{2,M}\), which contradicts distinct primitive rows. The minus sign gives \(y/w'+x/w=1\). At the whole private modulus \(D_1\), the denominators \(w,w'\) are units and both ratios are \(a_2/c_2\) modulo \(D_1\). Hence \(D_1\mid (2a_2-c_2)=a_2-b_2\). Its magnitude is at most \(X\) and \(D_1>X\), so \(a_2=b_2\), giving the excluded primitive row. All degeneracies are therefore absent. Summing the fixed-target bound proves the claim.

## Numerical height of the ratio

Write \(S_0=\mathcal R/T=P_0/Q_0\) in lowest terms, with \(T\)
as in the norm equation. By the ratio identity,
\(S_0=U_M V_M/(U_{M+1}V_{M+1})\).
It is a unit at every private deficient prime. Modulo \(D_i\),
the earlier and later differences reduce to the contributions
of the other row, so
\(S_0\equiv Q_{3-i,M}^{-2}\pmod {D_i}\).
Equality as rationals is impossible: on the private deficient
primes of row \(3-i\), the right side has negative valuation
and the left side is a unit. Hence the nonzero integer
\(Q_{3-i,M}^2 P_0-Q_0\) is divisible by \(D_i\).
With \(H\) denoting the maximum of numerator magnitude and
denominator of a reduced rational, \(H(T)\le n^{3r_+}\) by the
definition of the \(Q_{i,m}\), and \(Q_{i,M}\le n^{r_+}\).
It follows that
\[
 H(\mathcal R)\ \ge\
 \frac{\max(D_1,D_2)}{2n^{5r_+}}
 \ \ge\ \frac12 X^{K(r_--\xi)-1-5r_+L}.
\]
This lower bound rules out a low numerical height for the exact
record itself, independently of its catalogue cardinality.

## Near-full balanced families

In a balanced height shell, \(r_-\) can be \(1/3\) minus a
fixed small cap tolerance. For any fixed first row, EB5 prices its
entire product support at zero. The second rows with shared radical
greater than \(n^\xi\) then have exponent at most
\(1/3-\xi/3\) plus the cap and count tolerances, uniformly for small
fixed \(L\). Thus on a near-full family of exponent \(d\) sufficiently
close to \(1/3\), a full proportion of pairs in exponent satisfy
the overlap cap. The dictionary of \(\mathcal R\) has size at least
\(n^{2d-\xi-2/K-o(1)}\), taking the height limits at fixed \(L\)
before the small-budget limit. This applies to all balanced
dictionary prices.

The successful recovery provides a lower bound on the size of the
actual ratio dictionary. A saving would need an upper bound on a
sufficiently large pair sector, or a reduction to a genuinely cheap
auxiliary. This ratio does not by itself carry the fresh primes of
the [FP5 product difference](ABC_FIELD_FREE_CROSS_PRODUCT_TRANSFERS.md#recovery-from-the-product-difference-alone).
For example the primitive rows \((1,4,5),(1,6,7)\) give \(Z=23\)
in FP1, yet at \(M=2\) the ratio \(\mathcal R=385/1164\) is a
unit at 23. The mod-23 calculation persists whenever
\(M\equiv2\bmod\operatorname{ord}_{23}(5/7)\), since the two
row multipliers are \(5^{M-1},7^{M-1}\) for \(M\ge2\), the
consecutive ratios \(Q_{i,m}\) are constant, and the common scaling
cancels from the relevant expressions modulo 23. The example
addresses only that proposed support implication.

[Least multipliers](ABC_LEAST_MULTIPLIER_SHIFT_COUNTS.md)
· [Effective radicals (EB5)](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)
· [Browning–Verzobio, Lemma 2.3](https://arxiv.org/html/2608.24512v1#S2.SS1).

# Dyadic prefix windows for factorial moments of a triangle

For \(F=a xy+b yz+c zx\), the leading binary layer at orders \((2^t-1)2^k\) restricts the possible reductions of a common zero of three consecutive moments. The layer factors into a power of the sum of the coefficients and a smaller factorial polynomial. At the first two prefixes, the resulting residue equations reduce to gcd tests for pairs of roots of unity that differ by one.

Write \(P_m=\Lambda(F^m)/(m!)^2\), \(s=a+b+c\), \(h=ab+bc+ca\), \(r=abc\), and \(p_j=a^j+b^j+c^j\). Use the [quadratic triangle recurrence](SFC_TRIANGLE_OVERLAP_MOMENTS.md#exact-recurrence)
\[
 \Delta P_m=\frac{2(2m+1)}m\,rP_{m-1}+W_{m+1},\qquad
 \Delta=4h-s^2,\qquad W_j=s p_j-2p_{j+1}.
 \tag{1}
\]

## The leading factorial layer

For \(D=2^e\), \(e\ge0\), put \(F_D=a(xy)^D+b(yz)^D+c(zx)^D\) and
\[
              P^{[D]}_m=\frac{\Lambda(F_D^m)}{2^{2m(D-1)}(m!)^2}.
\]
Thus \(P^{[1]}_m=P_m\). Reduction of an integral polynomial at \(2\) is denoted by a bar.

**Proposition (Dyadic prefix reduction).** Let \(Q=2^k\), \(k\ge2\), and \(L=2^t-1\), \(t\ge1\). In \(\mathbb Z_{(2)}[a,b,c]\), \(2^{k-1}P^{[D]}_{LQ}\) is integral and its coefficientwise reduction modulo two is
\[
         \overline{2^{k-1} P^{[D]}_{LQ}}
                  =\bar r\,\bar s^{\,Q-3}\bigl(\overline{P_{L-1}}\bigr)^Q .
 \tag{2}
\]
In particular the last factor is \(1\) for \(L=1\) and \((\bar s^2+\bar h)^Q\) for \(L=3\). The smaller polynomial \(P_{L-1}\) is integral at two.

First take \(D=1\). At a monomial with exponent triple \(i+j+l=LQ\), the contribution at \(q=2^u\) to the coefficient valuation is
\[
 \sum_{x<y}\left\lfloor\frac{\rho_x+\rho_y}{q}\right\rfloor
                 -\left\lfloor\frac{\rho_1+\rho_2+\rho_3}{q}\right\rfloor ,
 \tag{3}
\]
where \(\rho_1,\rho_2,\rho_3\) are the three remainders. This is the [factorial floor formula](SFC_TRIANGLE_DYADIC_BOUNDARY.md#the-backwards-recurrence-at-two). The contribution at \(q=2\) is zero since \(LQ\) is even. Each \(4\le q\le Q\) contributes at least \(-1\), with equality exactly when the remainders are all positive and sum to \(q\). To see this, a sum of \(2q\) forces all pairs to carry, while at a sum of \(q\) the pair opposite a zero remainder is the only possible carry.

For \(Q<q\le LQ\), the remainder of \(LQ\) modulo \(q\) is \(q-Q\ge q/2\) since \(L=2^t-1\). A contribution of \(-1\) would require a remainder sum of \(q+(q-Q)\) with no pair carrying, so each opposite remainder exceeds \(q-Q\). Summing those three strict bounds contradicts \(q-Q\ge q/2\). A remainder sum of at least \(2q\) always makes every pair carry, and all levels above \(LQ\) vanish. Thus the total valuation is at least \(1-k\).

To attain the minimum, the remainder triple at \(q=4\) is a permutation of \((1,1,2)\). At each doubling up to \(Q\), assign the new digit to exactly one index so that the remainder sum doubles. These choices give \(\bar r\bar s\prod_{u=2}^{k-1}\bar p_{2^u}=\bar r\bar s^{Q-3}\). Write the exponents as these three positive lower remainders plus \(Q(I,J,K)\), with \(I+J+K=L-1\). For \(q=2^vQ>Q\), each lower pair sum is less than \(Q\) and the lower triple sum is \(Q\). The pair floors in (3) are exactly the pair floors of the high remainders of \(I,J,K\) modulo \(2^v\). The extra \(1\) in the triple floor has no effect since the sum of the high remainders is even. Thus the upper levels give exactly the factorial coefficient valuation at \((I,J,K)\) in \(P_{L-1}\). This polynomial is integral by the [dyadic boundary reduction for \(2^t-2\)](SFC_TRIANGLE_DYADIC_BOUNDARY.md#the-factorial-coefficients-at-a-binary-boundary), or directly if \(t=1\). Its valuation-zero monomials give \(\overline{P_{L-1}}(a^Q,b^Q,c^Q)=(\overline{P_{L-1}})^Q\) in the residue field. Every surviving unit coefficient is one modulo two, proving (2).

For \(D=2^e\), the identity
\[
           (2u)!=2^u u!\prod_{j=1}^{u}(2j-1)
\]
iterated \(e\) times shows that each coefficient of \(P_m^{[D]}\) differs multiplicatively from its counterpart in \(P_m\) by an odd integer. Thus the valuations and the leading coefficientwise reductions coincide, giving (2) also for \(D>1\).

## Restrictions on a possible common zero

**Theorem (Residue restriction above a dyadic prefix).** Suppose \(P_n=P_{n+1}=P_{n+2}=0\) with \(abc\ne0\). Let \(n=LQ+d\), \(Q=2^k\), \(k\ge2\), \(L=2^t-1\), \(t\ge1\), and \(0\le d\le Q/2-2\). Work over a number field and at any place above two, scaling to make the coefficients integral with at least one unit. Then all three are units, and their reductions satisfy
\[
                     \bar s\,\overline{P_{L-1}}=0 .                \tag{4}
\]
If \(\bar s=0\), the three residues are distinct, and after scaling one to \(1\) the other two are \(y,1+y\) with
\[
                       y^{n+1}=(1+y)^{n+1}=1.                     \tag{5}
\]

As in the [dyadic boundary argument](SFC_TRIANGLE_DYADIC_BOUNDARY.md#the-backwards-recurrence-at-two), an algebraic common zero follows from a complex one by the weak Nullstellensatz on \(a=1\) with \(bc\) inverted. The proof also applies to any algebraic representative of a common zero. Normalize the valuation by \(\nu(2)=1\). The recurrence gives \(W_{n+2}=W_{n+3}=0\), and the [dyadic two-unit test](SFC_TRIANGLE_DYADIC_BOUNDARY.md#the-reduction-with-two-unit-coefficients) shows that either all coefficients are units, or two have equal unit residues and the third has valuation \(0<\rho\le1/(n+2)\). A single unit or two distinct unit residues contradict the reduced \(W\) equations. In the two-unit case \(\nu(r)=\nu(s)=\rho\), since after scaling to \(a=1\) the test gives \(\nu(b-1)=(n+2)\rho>\rho\). The all-unit case has \(\nu(r)=\rho=0\).

Run (1) backwards from \(P_n=0\) to \(P_{LQ}\). Each term introduced by \(W_{LQ+j+1}\), \(1\le j\le d\), includes the multiplier
\[
       \Delta^{j-1}\prod_{i=1}^j\frac{LQ+i}{2(2(LQ+i)+1)r}.
\]
Since \(i<Q\), its valuation is at least \(\nu(j!)-j(1+\rho)=-s_2(j)-j\rho\), with \(s_2(j)\) the binary digit sum. The \(W\)'s are integral here. The smallest nonnegative integer with \(k-1\) one-bits is \(Q/2-1\), so
\[
                         \nu(P_{LQ})\ge 2-k-d\rho .                \tag{6}
\]
This also holds when \(d=0\), because that moment is zero.

In the two-unit case \(\overline{P_{L-1}}(\bar a,\bar b,0)\ne0\). Indeed, with \(c=0\) the polynomial is the complete homogeneous polynomial in \(a,b\), and at equal residues it equals \((L)\bar a^{L-1}=\bar a^{L-1}\). Therefore the term \(r s^{Q-3}(P_{L-1})^Q\) in (2) has valuation \((Q-2)\rho<1\). The congruence gives
\[
                     \nu(P_{LQ})=1-k+(Q-2)\rho<2-k-d\rho ,
\]
since \((Q-2+d)\rho<1\). This contradicts (6). For three units, a unit value of \(sP_{L-1}\) would give \(\nu(P_{LQ})=1-k\) by (2), also impossible in (6). This proves (4).

When \(\bar s=0\), distinctness of the residues follows from characteristic two and the unit condition. The [all-order dyadic restriction](SFC_TRIANGLE_DYADIC_BOUNDARY.md#the-backwards-recurrence-at-two) gives \(0<\nu(s)\le1/3\), with the bound \(1/4\) for even \(n\). Dividing the two exact \(W\) equations by \(s\) gives \(\bar p_{n+2}=\bar p_{n+3}=0\), since \(2/s\) has positive valuation. Therefore \((\bar a^{n+2},\bar b^{n+2},\bar c^{n+2})\) is in the nullspace of the rows \((1,1,1),(\bar a,\bar b,\bar c)\). The vector \((\bar a,\bar b,\bar c)\) belongs to the same one-dimensional nullspace, since a vanishing sum of the residues has vanishing squared sum in characteristic two. Scaling one residue to \(1\) proves (5). In particular, at such a place the ramification index is at least three, or at least four for even \(n\), and its residue field must contain a root of the gcd for (5).

## The first two prefix tests

For a positive integer \(M\), let \(M_{\mathrm{odd}}\) be its positive odd part and define over \(\mathbb F_2\)
\[
           G_M(T)=\gcd\bigl(T^{M_{\mathrm{odd}}}-1,
                              (1+T)^{M_{\mathrm{odd}}}-1\bigr).
\]

**Corollary (Root-of-unity tests for the first two prefixes).** In the window above \(Q\) the three moments have no common complex zero with \(abc\ne0\) if \(G_{n+1}=1\). In the window above \(3Q\) the same exclusion holds if \(G_{n+1}=1\) and \(G_{3(n+2)}=T^2+T+1\).

For the first window \(L=1\) and (4) forces \(\bar s=0\), so (5) gives the assertion. For \(L=3\) the same argument handles \(\bar s=0\). Otherwise scale the residues to \(s=1\), and use \(h=1\) from (4). In the following argument the coefficients denote their residue values. The cubic with roots \(a,b,c\) is \(T^3+T^2+T+r=(T+1)^3+(1+r)\). The roots cannot all be \(1\), because their power sum \(p_{n+2}\) must vanish when \(s\) is a unit. Thus \(a+1,b+1,c+1\) are the three distinct nonzero cube roots of \(1+r\). In some order
\[
                b+1=\omega(a+1),\qquad c+1=\omega^2(a+1),
\]
where \(\omega^2+\omega+1=0\). The two reduced \(W\) equations again give a one-dimensional nullspace. This time it is generated by \((b+c,a+c,a+b)=(a+1,b+1,c+1)\), so \((b/a)^{n+2}=\omega\), \((c/a)^{n+2}=\omega^2\). The equality \(a+c=\omega(b+c)\) gives
\[
                         \omega b/a+\omega^2 c/a=1 .
\]
Both summands are \(3(n+2)\)-th roots of unity. If \(G_{3(n+2)}=T^2+T+1\), both are cube roots of unity, so \(b/a,c/a\) are distinct nonidentity cube roots of unity. This gives \(s=0\), a contradiction.

An elementary sufficient condition for \(G_M=1\) is \(3\nmid M\) and \(M_{\mathrm{odd}}\mid 2^f+1\) for some \(f\ge1\). On both roots of unity in the gcd, the \(2^f\)-power map would act by inversion. Applying it to \(y+(1+y)=1\) gives \(y^{-1}+(1+y)^{-1}=1\), or \(y^2+y+1=0\), which requires \(3\mid M\). The same proof gives \(G_M=T^2+T+1\) when \(M_{\mathrm{odd}}\mid 2^f+1\) and \(3\mid M\), since the gcd is squarefree by oddness of \(M_{\mathrm{odd}}\).

## Orders outside the earlier windows

The first test covers \(n=(2^v+1)2^u-1\) when \(v\ge2\) is even and \(u\ge0\). Here \(n+1\) has odd part \(2^v+1\), and \(n\) is in the window above \(Q=2^{u+v}\). The exact gcd can also be one without a Frobenius power acting by inversion. For \(M=23\), a Bézout certificate over \(\mathbb F_2[T]\) is
\[
 \left(\sum_{\substack{0\le i\le15\i\ne8}}T^i\right)(T^{23}+1)
              +(T^{15}+T^8+1)((1+T)^{23}+1)=1 .
\]
The order of \(2\) modulo \(23\) is \(11\). This gives \(n=23\cdot2^u-1\) in the window above \(16\cdot2^u\), for every \(u\ge0\). For \(u\equiv8\bmod26311740\), the six odd arguments \(2n+1-2t\), \(0\le t\le5\), have valuation one at the respective primes \(3,61,79,3,7,5\). Indeed this holds at \(u=8\) and \(2^{26311740}=1\) modulo the squares of these primes. These arguments exceed the indicated primes, so none is a prime power. No window in this subsequence contains a dyadic boundary of the form \(2^k-1\).

For the second test take \(n=826\) and \(Q=256\). The two odd exponents are \(827\) and \(621\), with \(G_{827}=1\) and \(G_{621}=T^2+T+1\). The six shifted arguments \(1653,1651,1649,1647,1645,1643\) have the respective factorizations \(3\cdot19\cdot29\), \(13\cdot127\), \(17\cdot97\), \(3^3\cdot61\), \(5\cdot7\cdot47\), and \(31\cdot53\). The window lies between the boundaries \(511\) and \(1023\). 
**Proposition (An infinite second-prefix family).** The three moments starting at \(n=25\cdot2^u-1\), \(u\ge0\), have no common complex zero with \(abc\ne0\).

Here \(Q=8\cdot2^u\) and \(n=3Q+2^u-1\) is in the required band. The odd part of \(n+1\) is \(25\), and \(2^{10}\equiv-1\pmod{25}\), so \(G_{n+1}=1\). To test the other exponent put \(M=3(25\cdot2^u+1)\). If \(y\) is a root of \(G_M\), take \(z=y^{2^u}\), and write
\[
                 yz^{25}=\xi,\qquad (1+y)(1+z)^{25}=\eta
\]
with \(\xi,\eta\) cube roots of unity in \(\mathbb F_4\). Then \(z\ne0,1\) and
\[
       f_{\xi,\eta}(z)=z^{25}(1+z)^{25}+\xi(1+z)^{25}+\eta z^{25}=0
\]
follows by adding the expressions for \(y\) and \(1+y\). This gives nine degree-50 polynomials independent of \(u\). For each one the [exact checker](verify_sfc_triangle_dyadic_prefix.py) computes \(g=\gcd(f_{\xi,\eta},z^{2^{480}}-z)\) by polynomial squaring and verifies \(f_{\xi,\eta}\mid g^2\). This places all possible \(z\) in \(\mathbb F_{2^{480}}\). Put \(Y=\xi/z^{25}\) in the quotient by \(g\), using an inverse obtained from the Euclidean algorithm. For every \(0\le j<480\) the checker verifies that \(\gcd(g,Y^{2^j}-z)\) divides \(z^2+z+1\). These are all Frobenius powers on that finite field. Thus every possible \(z\), and therefore \(y\), is a nonidentity cube root. Both nonidentity cube roots occur since \(3\mid M\). This gives \(G_M=T^2+T+1\) and proves the exclusion.

For \(u\equiv14\pmod{22932289486260}\), the six shifted arguments have valuations \(1,1,1,1,1,2\) at \(191,461,3,43,23,3\), respectively, and exceed the corresponding prime powers. The checker verifies the base values and the period modulo the six next prime powers, so this gives an infinite family outside the earlier prime-power tests. These windows lie strictly between consecutive dyadic boundary orders \(2^k-1\).

Replay from the repository root with Sage.

~~~sh
sage -python research/factorial/image-rigidity/verify_sfc_triangle_dyadic_prefix.py
~~~

The unrestricted three-moment problem remains open. The backwards step uses the \(D=1\) recurrence.

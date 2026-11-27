# Dyadic windows for the quadratic triangle

For \(F=a xy+b yz+c zx\), a window of three consecutive factorial moments containing the order \(2^k-1\) cannot vanish when \(abc\ne0\). The proof uses the exact factorial recurrence both forwards and backwards. At the prime \(2\) the backwards recurrence restricts the discriminant of a possible common zero, and the binary digits of the factorial coefficients give two reductions at the boundary of a dyadic interval.

We use [the quadratic triangle recurrence](SFC_TRIANGLE_OVERLAP_MOMENTS.md#exact-recurrence) with \(P_m=\Lambda(F^m)/(m!)^2\),
\[
 s=a+b+c,\qquad h=ab+bc+ca,\qquad r=abc,\qquad
 \Delta=4h-s^2,\qquad W_j=\sum_{x\in(a,b,c)}(s-2x)x^j .
\]
Three zero moments at \(n,n+1,n+2\) give \(W_{n+2}=W_{n+3}=0\). The [characteristic-zero consequences](SFC_TRIANGLE_OVERLAP_MOMENTS.md#consequences-in-characteristic-zero) exclude \(\Delta=0\), repeated coordinates and zero weights \(s-2x\) when \(r\ne0\). They also give
\[
 ((s-2a)a^j,(s-2b)b^j,(s-2c)c^j)=\lambda(b-c,c-a,a-b),
 \qquad j=n+2,\quad \lambda\ne0.                                  \tag{1}
\]
The support-sized formulation comes from [Edo and van den Essen, Definition 2.7 and Conjecture 2.8](https://arxiv.org/pdf/1304.3956v2). The simplex identity in the recurrence is also used in [Lasserre, Theorem 1.1](https://arxiv.org/pdf/1908.06736).

## The backwards recurrence at two

We may work over \(\overline{\mathbb Q}\). A complex common zero with the required nonvanishings would give an algebraic one by using the chart \(a=1\), inverting the required product, and applying the weak Nullstellensatz over \(\overline{\mathbb Q}\). At a place above \(2\), use a valuation \(\nu(2)=1\) and scale the coefficients to be integral with at least one unit. Denote reduction by an overbar and work in the completed number field when summing series. Write \(s_2(m)\) for the sum of the binary digits of \(m\).

**Proposition (All-order dyadic restriction).** Under the preceding normalization, an algebraic common zero of \(P_n,P_{n+1},P_{n+2}\) with \(abc\ne0\) must satisfy \(\nu(\Delta/r)<1\) at every place above \(2\). If all three coefficients are units there, then \(\nu(s)\le1/3\), and \(\nu(s)\le1/4\) if \(n\) is even.

Put \(B_m=(2m+1)!/(m!)^2\), so \(B_0=1\), \(\nu(B_m)=s_2(m)\), and the recurrence coefficient \(2(2m+1)/m\) is \(B_m/B_{m-1}\). For \(d=\Delta/r\) this gives
\[
 \frac{d^m P_m}{B_m}-\frac{d^{m-1}P_{m-1}}{B_{m-1}}
                       =\frac{d^m W_{m+1}}{\Delta B_m}.            \tag{2}
\]
The coefficient of \(a^i b^j c^k\) in \(P_m\), \(i+j+k=m\), is
\[
             \frac{(i+j)!(j+k)!(k+i)!}{m!\,i!\,j!\,k!}.            \tag{3}
\]
For each power \(q=2^u\), \(u\ge1\), its contribution in the factorial floor formula is
\[
 \sum_{x<y}\Big\lfloor\frac{r_x+r_y}{q}\Big\rfloor
                              -\Big\lfloor\frac{\sum_x r_x}{q}\Big\rfloor, \tag{4}
\]
where \(r_x\) are the three remainders of \(i,j,k\) modulo \(q\), indexed by position. This is at least \(-1\): a total at least \(2q\) makes every pair carry, and a smaller total subtracts at most one. All contributions for \(q>m\) vanish. Thus \(\nu(P_m)\ge-\lfloor\log_2 m\rfloor\) for \(m\ge1\). The formula \(\nu(B_m)=s_2(m)\) follows from \(\nu(m!)=m-s_2(m)\), obtained by summing the factorial floors.

If \(\nu(d)>0\), the first term of (2) tends to zero as \(m\to\infty\). Under a three-moment zero, the sum of the right-hand terms from \(m=n+3\) onwards is therefore zero. Here \(W_{n+4}\ne0\), since three consecutive \(W\)'s would give the impossible grouped Vandermonde equations. Write \(H_\ell(a,b,c)\) for the complete homogeneous polynomial of degree \(\ell\), with \(H_0=1\) and \(H_1=s\). The recurrence of a sum of powers with the three bases \(a,b,c\), together with \(W_{n+2}=W_{n+3}=0\), gives
\[
                        W_{n+4+\ell}=W_{n+4}H_\ell\qquad(\ell\ge0).
\]
Divide the convergent tail by its first term. For \(M=n+3\) the resulting equation is
\[
                 1+\sum_{\ell\ge1}d^\ell
                          \frac{B_M}{B_{M+\ell}}H_\ell=0.         \tag{5}
\]
The valuation of the term of index \(\ell\) is at least
\(\ell\nu(d)-s_2(\ell)+\nu(H_\ell)\), by subadditivity of binary digit sums. If \(\nu(d)\ge1\), all terms with \(\ell\ge2\) have strictly positive valuation. Also \(s\) is a nonunit in that case, because a unit \(s\) gives \(\nu(\Delta)=0\) and \(\nu(r)\ge0\). Thus the term with \(\ell=1\) also has positive valuation, impossible in (5).

When \(a,b,c\) are all units this gives \(\nu(s)<1/2\) and \(\nu(d)=2\nu(s)\). If \(\nu(s)>1/3\), the terms with \(\ell\ge2\) again have positive valuation, since \(s_2(\ell)\le 2\ell/3\). The first term after the constant has valuation at least \(3\nu(s)-1>0\). For even \(n\), \(M\) is odd. Adding an odd \(\ell\) to \(M\) makes at least one carry, so \(s_2(M+\ell)-s_2(M)\le s_2(\ell)-1\) for such \(\ell\). For even \(\ell\), the difference is at most \(s_2(\ell)\). Both bounds are at most \(\ell/2\), so (5) excludes \(\nu(s)>1/4\) when \(n\) is even.

## The reduction with two unit coefficients

Suppose \(W_j=W_{j+1}=0\), \(j\ge3\), with distinct nonzero coefficients and nonzero weights. At a place above \(2\), assume that two coefficients are units with equal residues and the third is a nonunit. After unit scaling and permutation write \(a=1,\ b=1+\xi\), and set \(w=\nu(\xi)>0,\ t=\nu(c)>0\). Then
\[
         w=jt\le1,\qquad w>t,\qquad
                 jt<1\ \hbox{if \(j\) is odd}.                    \tag{6}
\]
Indeed (1) gives \(u=\nu(c+\xi)=\nu(c-\xi)\) and
\(w+u=jt+\nu(2+\xi-c)\). If \(w<t\), then \(u=w\) and \(jt>2w\), impossible. If \(w=t\), the difference of \(c+\xi\) and \(c-\xi\) gives \(u\le1+t\). For \(u<1\) one has \(\nu(2+\xi-c)=u\), making \(w=jt\) impossible. For \(u=1\), the right-hand side is at least \(jt+1>t+1=w+u\). For \(u>1\), it is \(jt+1>1+2t\ge w+u\). Hence \(w>t\) and \(u=t\).

For any exponent \(d_0\) the power sum with weights is
\[
 W_{d_0}=c(1+b^{d_0})+\xi(1-b^{d_0})+(2+\xi-c)c^{d_0}.
\]
At \(d_0=j\) the third term has valuation \(t+w\), and the middle term has valuation at least \(2w\). If \(w>1\), then \(\nu(1-b^j)\ge w>1\), so the first term has valuation \(t+1<t+w\), impossible. Thus \(w\le1\), \(t<1\), and \(\nu(2+\xi-c)=t\), giving \(w=jt\). If \(w=1\) and \(j\) is odd, then \(\nu(b^{j+1}-1)>1\) because \(j+1\) is even and \(\nu(b^2-1)\ge2\). At \(d_0=j+1\) the first term now has valuation \(t+1\), below the middle term of valuation greater than \(2\) and the third term of valuation \(1+2t\). This proves (6).

**Corollary (Restrictions on the coefficient field).** Three consecutive moments cannot vanish at any order if the projective coefficients with \(abc\ne0\) are defined over a number field of degree at most two. Over a cubic field a common zero requires that \(2\) be inert and \(n\equiv3\bmod7\). More generally, at a dyadic place of residue degree one, a common zero requires ramification index \(e\ge n+2\), with strict inequality when \(n\) is odd.

To see this, if exactly one coefficient is a unit then reduction of \(W_j\) is nonzero. If exactly two are units with distinct residues, their sum is a unit, so reduction of the two \(W\) equations would give \(\bar a^j=\bar b^j\) and \(\bar a^{j+1}=\bar b^{j+1}\) for the two unit coordinates. This makes the residues equal. With two equal unit residues, (6) applies and the positive valuation \(t\) is at least \(1/e\). Thus \(e\ge j=n+2\), strictly for odd \(j\). If all coefficients are units with residue degree one, \(W_j\) reduces to \(1\).

At a residue degree two place with \(e\le2\), the two-unit case is excluded. The all-order restriction for three units makes \(s\) a unit, since a positive valuation would be at least \(1/2\). Over \(\mathbb F_4\) each \(p_d=\sum_{x\in(a,b,c)}x^d\) reduces to \(1\), \(\bar s\), or \(\bar s^2\) according to \(d\bmod3\), so \(\bar W_j=\bar s\,\bar p_j\ne0\). These arguments exclude every possible local degree at most two. In a cubic field every decomposition at \(2\) except a single unramified prime of residue degree three has a residue degree one place with \(e\le3\). At that place the bound on \(e\), including strictness at \(j=3\), excludes a common zero. At the unramified residue degree three place all coefficients and \(s\) must be units. In \(\mathbb F_8\) a power sum \(p_d\) of the three units is nonzero at \(d\equiv0,1,2,4\bmod7\), since it equals \(1\) or a Frobenius power of \(\bar s\). The only consecutive pair of exponents outside these classes is \(j\equiv5,\ j+1\equiv6\bmod7\), which requires \(n\equiv3\bmod7\).

## The factorial coefficients at a binary boundary

Put \(Q=2^k\) for \(k\ge2\) and \(p_d=a^d+b^d+c^d\).

**Lemma (Two dyadic boundary reductions).** The polynomials \(P_{Q-1},P_{Q-2}\) have coefficients in \(\mathbb Z_{(2)}\). The first has the coefficientwise reduction
\[
                       P_{Q-1}\equiv\prod_{i=0}^{k-1}p_{2^i}
                                          \pmod4.                \tag{7}
\]
For the second, on \(s=0\) in characteristic two the reduction is
\[
                       \overline P_{Q-2}
                                =(ab)^{Q/2-1}+(bc)^{Q/2-1}
                                  +(ca)^{Q/2-1}.                 \tag{8}
\]

For \(m=Q-1\), the sum of the three remainders in (4) is \(q-1\) or \(2q-1\) at each \(q\le Q\). In the first case the contribution is zero. In the second all three pairs carry and the contribution is two. Thus terms with a carry have valuation at least two, and the valuation-zero terms are precisely the assignments of the set binary digits of \(m\) to the three indices without carries.

To compute the unit coefficients of these terms modulo \(4\), the odd part of \(N!\) is \((-1)^{e(N)}\) modulo \(4\), with
\[
 e(N)=\sum_{h\ge0}
           \left\lfloor\frac{\lfloor N/2^h\rfloor+1}{4}\right\rfloor.
\]
This counts the odd factors congruent to \(3\) modulo \(4\) at each fixed power of \(2\) in the integers up to \(N\). Decompose each summand into \(\lfloor N/2^{h+2}\rfloor\) and the indicator that binary digits \(h,h+1\) of \(N\) are both one. For a carry-free term, the floor parts cancel between numerator and denominator of (3) because the binary supports of \(i,j,k\) are disjoint. For adjacent digits of \(m\) both assigned to the same index, there are two occurrences of the adjacent pair among \(i+j,j+k,k+i\) and two among \(m,i,j,k\). If assigned to distinct indices, there is one occurrence on each side. Hence the signs cancel and each carry-free coefficient is \(1\) modulo \(4\). Summing those assignments proves (7). In particular \(P_{Q-1}\equiv s^{Q-1}\pmod2\).

For \(m=Q-2\), the sum of the remainders is \(q-2\) or \(2q-2\). In the first case the contribution is zero. In the second a pair carries exactly when the opposite remainder is at most \(q-2\). At least one remainder is at most \(q-2\), and the contribution vanishes exactly when the remainders are \(q-1,q-1,0\) in some order. Consequently the coefficients are integral. For valuation zero, a level with sum \(2q-2\) forces the same pattern for every smaller power of \(2\). If the highest such level is \(q=2^j\), where \(1\le j\le k-1\), the two selected indices contain all the first \(j\) digits, the next digit of each index is zero, and the higher digits of \(m\) are assigned without carries. This gives
\[
 \overline P_{Q-2}
 =s^{Q-2}+\sum_{j=1}^{k-1}
       \bigl((ab)^{2^j-1}+(bc)^{2^j-1}+(ca)^{2^j-1}\bigr)
                        s^{Q-2^{j+1}} .
\]
The first term accounts for the assignments without any carry. At \(s=0\) only \(j=k-1\) contributes, giving (8).

## Windows containing \(2^k-1\)

**Theorem (Three consecutive moments across a dyadic boundary).** For \(k\ge1\) and \(n\ge1\) with \(n\le2^k-1\le n+2\), the three moments \(P_n,P_{n+1},P_{n+2}\) have no common complex zero with \(abc\ne0\).

For \(k=1\) use \(P_1=s,\ P_2=s^2-h,\ P_3=s^3-2sh+4r/3\). Take \(k\ge2\), \(Q=2^k\), and a putative common zero as above. By (7) at least one moment has reduction \(\bar s^{Q-1}\), so \(\bar s=0\). There are three distinct unit residues or two equal unit residues and a nonunit third coefficient. In the latter case (6) applies with \(j=n+2\in\{Q-1,Q,Q+1\}\). In each case \((Q-1)t<1\), using the strict bound for the odd exponent \(j=Q-1\). But \(P_{Q-1}(1,1+\xi,c)\) differs from \((\xi+c)^{Q-1}\) by a term of valuation at least one, by reduction modulo \(2\). Since \(w>t\), the latter power has valuation \((Q-1)t<1\), a contradiction.

It remains to take three distinct unit residues. The all-order restriction gives \(u=\nu(s)<1/2\). In \(W_j=W_{j+1}=0\) divide by \(s\) and reduce modulo the maximal ideal. This gives \(\bar p_j=\bar p_{j+1}=0\), since \(\nu(2/s)>0\). In the residue field scale \(a=1\), write \(b=y,c=1+y\), where \(y\ne0,1\). Then
\[
       \bar p_Q=0,\qquad \bar p_{Q+1}=y^Q+y,\qquad
       \bar p_{Q+2}=y^Q+y^2,\qquad
                       \bar r\,\bar p_{Q-1}=\bar p_{Q+2}.
\]
The last identity uses \(bc=h+a^2\) when \(s=0\) in characteristic two and its permutations. With \(n=Q-1\), the two power equations now require \(y=y^2\), impossible. With \(n=Q-3\), they give \(y^Q=y^2\). This window also contains \(P_{Q-2}\). On \(s=0\), the right-hand side of (8) multiplied by \(\bar r\) equals \(\bar p_{Q+1}\): expand using \(c=a+b\) and \(c^{Q/2}=a^{Q/2}+b^{Q/2}\). Vanishing of \(P_{Q-2}\) forces \(y^Q=y\), again impossible.

For \(n=Q-2\), the residue equations give \(y^Q=y\), so \(\bar p_{Q+2}=\bar p_3=\bar r\ne0\). The exact two equations \(s p_Q=2p_{Q+1}\) and \(s p_{Q+1}=2p_{Q+2}\) give
\[
                     \nu(p_{Q+1})=1-u,\qquad
                     \nu(p_Q)=2(1-u)\in(1,2).
\]
We show that the product in (7) has valuation below \(2\). If \(\bar h\ne0\), use
\[
 p_{2d_0}=p_{d_0}^2-2\big((ab)^{d_0}+(bc)^{d_0}+(ca)^{d_0}\big).
\]
At \(d_0=Q/2\) the subtracted term has valuation one. Since \(\nu(p_Q)>1\), this forces \(\nu(p_{Q/2})=1/2\). The same identity backwards halves the valuation at each earlier index \(2^i\), since the subtracted terms have valuation at least one. Thus \(\sum_{i=0}^{k-1}\nu(p_{2^i})<1\).

If \(\bar h=0\), the three residue coordinates after scaling \(a=1\) are the cube roots of unity. Since \(y^Q=y\), \(Q\equiv1\bmod3\). In a finite extension write \(b=\omega+B,\ c=\omega^2+C\), with \(\omega\) reducing to \(y\) and \(\omega^2+\omega+1=0\). Let \(\tau=\min(\nu(B),\nu(C))\le u\). The constant term of \(p_{Q+1}\) is zero and the linear term is \((Q+1)(\omega B+\omega^2 C)\). The higher terms have valuation at least \(2\tau>\tau\). Since \(\nu(p_{Q+1})=1-u>\tau\), the linear combination \(\omega B+\omega^2 C\) has valuation greater than \(\tau\). The matrix of that combination and \(B+C=s\) is invertible over the valuation ring, so \(\nu(B)=\nu(C)=u=\tau\).

In the expansion of \(p_Q\), the constant term is zero and \(B^Q+C^Q\) has valuation \(Qu\): its normalized leading residue is the \(Q\)-th power of the nonzero leading residue of \(s\). Every intermediate coefficient \(\binom{Q}{i}\) has valuation \(k-\nu(i)\), by \(\binom{Q}{i}=Q\binom{Q-1}{i-1}/i\) and the oddness of the latter binomial coefficient. If \(Qu\ge2\), all intermediate terms have valuation at least two, including the term \(i=Q/2\) with coefficient of valuation one. This contradicts \(\nu(p_Q)<2\). Hence \(Qu<2\), and all intermediate terms have valuation strictly greater than \(Qu\), since at \(i=Q/2\) the valuation is at least \(1+Qu/2\) and at the other indices it is at least two. For every \(i<k\), the congruence \(p_{2^i}\equiv s^{2^i}\pmod2\) now gives \(\nu(p_{2^i})=2^i u\), since \(2^i u<1\). The product in (7) has valuation \((Q-1)u<2\). This contradicts \(P_{Q-1}=0\) also in this final case.

These windows include infinitely many starts outside the [initial and first five shifted prime-power tests](SFC_TRIANGLE_OVERLAP_MOMENTS.md#prime-power-windows). For the three choices \(n=Q-3,Q-2,Q-1\), their eight possible odd arguments are \(2^{k+1}-b\), with \(b=1,3,\ldots,15\). At exponent \(k+1=8\), these are \(255,253,251,249,247,245,243,241\). Pick a prime factor of each, say \(3,11,251,3,13,5,3,241\), and keep the valuation at that prime fixed by taking the exponent congruent to \(8\) modulo a common multiple of the orders of \(2\) modulo the corresponding next higher prime powers. At sufficiently large exponents in that progression all eight numbers exceed their fixed prime powers and hence are composite with at least two distinct prime factors. For example, at \(n=1023\) the six arguments are \(2047,2045,2043,2041,2039,2037\). Only \(2039\) is a prime power, and it lies at the fourth shift, whose proved primes are \(11,13,17\). The general three-moment problem outside the proved windows remains open, as does the extension of this recurrence to \(D>1\).

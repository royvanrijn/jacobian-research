# Consecutive factorial moments of a quadratic triangle

Let \(F=a xy+b yz+c zx\), and let \(\Lambda\) take each monomial to the product of its exponent factorials. Set \(P_m=\Lambda(F^m)/(m!)^2\), \(P_0=1\). All coefficient statements below hold over every characteristic-zero field. The three-moment assertions concern \(abc\ne0\); if exactly two coefficients are nonzero the normalized moments are the complete homogeneous polynomials in those coefficients and two consecutive moments suffice.

## Exact recurrence

Put \(s=a+b+c\), \(\Delta=2(ab+bc+ca)-(a^2+b^2+c^2)\) and \(W_j=\sum_{t\in(a,b,c)}(s-2t)t^j\) with multiplicity. For every \(m\ge1\),
\[
       \Delta P_m=\frac{2(2m+1)}{m}abc\,P_{m-1}+W_{m+1}.           \tag{1}
\]

Indeed, on the triangle \(x,y\ge0,\ x+y\le1\), write \(z=1-x-y\) and \(Q=F(x,y,z)\). The beta integrals of \(x^i y^j z^k\) are \(i!j!k!/(i+j+k+2)!\), so \(J_m=\int Q^m=\Lambda(F^m)/(2m+2)!\). This is the homogeneous simplex identity also used by [Lasserre, Theorem 1.1](https://arxiv.org/pdf/1908.06736). For \(\Delta\ne0\) the stationary point of \(Q\) has barycentric coordinates
\[
 (x_0,y_0,z_0)=\Delta^{-1}\bigl(b(s-2b),c(s-2c),a(s-2a)\bigr),
 \qquad Q_0=abc/\Delta .
\]
These formulas follow by substituting into \(x_0+y_0+z_0=1\), the two affine gradient equations and \(Q\). The vector \((x-x_0,y-y_0)\) differentiates \(Q\) to \(2(Q-Q_0)\). Apply the divergence theorem to that vector times \(Q^m\). On \(x=0\) the flux is \(x_0 b^m\int_0^1 t^m(1-t)^m\,dt\), on \(y=0\) it is the analogous term \(y_0 c^m\) times the beta integral, and on \(x+y=1\) the normal factor times the line element yields \(z_0 a^m\) times that integral. Consequently
\[
 2(m+1)J_m-2m Q_0 J_{m-1}=\frac{(m!)^2}{(2m+1)!}
                       (x_0 b^m+y_0 c^m+z_0 a^m).
\]
Substitute the stationary values and clear denominators. The resulting identity is polynomial in the coefficients over \(\mathbb Q\) and extends to \(\Delta=0\) and all characteristic-zero fields.

## Consequences in characteristic zero

If \(P_n=P_{n+1}=P_{n+2}=0\), then
\[
                         W_{n+2}=W_{n+3}=0.                        \tag{2}
\]
Three consecutive \(W_j\) are impossible when \(abc\ne0\). Group the equal nonzero coefficients. A Vandermonde matrix in the distinct values \(t\) shows that three such zeros would require \(m_t(s-2t)=0\) for each value of multiplicity \(m_t\). Two different values cannot both equal \(s/2\), and with only one value \(t\ne0\), the weight is \(3t\ne0\). In particular four consecutive \(P_m\) never vanish. If \(\Delta=0\), (1) and three zero \(P_m\) would give three consecutive zero \(W_j\), so the three-moment claim holds on this locus.

For repeated coefficients the same Vandermonde test with at most two distinct values directly excludes (2). If there are three distinct coefficients and one equals \(s/2\), there are at most two nonzero weights and (2) is again excluded. At three distinct nonzero coefficients, (2) is equivalent to
\[
 \bigl((s-2a)a^{n+2},(s-2b)b^{n+2},(s-2c)c^{n+2}\bigr)
                      =\lambda(b-c,c-a,a-b)                       \tag{3}
\]
for some \(\lambda\ne0\). Over \(\mathbb C\) this also excludes a common modulus, outside the already resolved \(\Delta=0\) case. After scaling to modulus one put \(S=\operatorname{Re}(a\bar b+b\bar c+c\bar a)\). Taking squared absolute values of the first entry of (3) gives \(3-2S+4\operatorname{Re}(b\bar c)=|\lambda|^2(2-2\operatorname{Re}(b\bar c))\), and the same for both other pairs. Thus all three real pairings, hence all three distances, agree. Three distinct points on the unit circle at equal distances are a rotated set of cube roots of unity. There \(s=ab+bc+ca=0\) and \(\Delta=0\).

The initial window also holds: \(P_1=s,\ P_2=s^2-(ab+bc+ca),\ P_3=s^3-2s(ab+bc+ca)+4abc/3\). With a zero coefficient, say \(c=0\), one directly has \(P_m=\sum_{i=0}^m a^i b^{m-i}\), and \(P_{m+1}-aP_m=b^{m+1}\) when both \(a,b\ne0\).

The two equations in (2) can hold on the nondegenerate locus. For \(n=1\), prescribe \(s=1\), \(ab+bc+ca=q\), \(abc=r\), with \(10q^2-6q+1=0\) and \(r=(-1+5q-4q^2)/5\). Take \(a,b,c\) to be the complex roots of \(T^3-T^2+qT-r\). Newton's identities give \(W_3=-1+5q-4q^2-5r=0\) and \(5W_4=(1-4q)(10q^2-6q+1)=0\) after the substitution for \(r\). Here \(r\ne0\), \(\Delta=4q-1\ne0\), and (2) forces distinct coefficients by the preceding Vandermonde argument. But \(P_1=1\). The moment equation supplies an additional constraint on the same projective point.

## Prime-power windows

**Proposition.** If \(n\ge1\) and \(2n+1=p^e\) for an odd prime \(p\) and \(e\ge1\), the three moments \(P_n,P_{n+1},P_{n+2}\) have no common zero with \(abc\ne0\).

A common zero over a characteristic-zero field would give one in \(\overline{\mathbb Q}\) by applying the weak Nullstellensatz to the equations with \(a=1\) and \(bc\) inverted. Work at a prime above \(p\) of a number field containing the coordinates. Homogeneously scale to make all three integral at this place, with at least one unit.

First, for any residue characteristic \(p>2\), a solution of \(W_j=W_{j+1}=0\) for \(j\ge3\) with distinct nonzero coordinates has all three coordinates units under this normalization. A single unit is impossible by reduction of \(W_j\). If there are two units, scale and permute to \(a=1,\ b=1+h,\ c\) nonunit. The equality of the two unit residues follows from the two reduced equations: distinct residues would have zero weights by the two-by-two Vandermonde test, which would require equal residues. Therefore \(w=\nu(h)>0\), \(t=\nu(c)>0\), and \(s-2c=2+h-c\) is a unit. Formula (3) with exponent \(j\), or the nullspace of its two rows, gives
\[
                  w+\nu(c+h)=w+\nu(c-h)=j t .
\]
If \(w=t\), both \(c+h\) and \(c-h\) have valuation \((j-1)t>t\), contradicting their sum \(2c\). For \(w\ne t\), the equation is \(w+\min(w,t)=jt\), forcing \(w>t\). Then
\[
 W_j=c(1+b^j)+h(1-b^j)+(2+h-c)c^j
\]
has its unique term of least valuation in \(c(1+b^j)\), since \(1+b^j\) is a unit. This is also impossible.

For \(m=n-1\), each coefficient of \(P_m\) has \(p\)-valuation at least \(-(e-1)\). In fact, the coefficient of \(a^i b^j c^k\), \(i+j+k=m\), is
\[
       \frac{(i+j)!(j+k)!(k+i)!}{m!\,i!\,j!\,k!}.
\]
Legendre's floor formula, at \(q=p^h\), cancels the quotient parts and leaves the sum of the three pair floors of the residues, minus the floor of their total, all divided by \(q\) inside the floor operation. This difference is at least \(-1\): when the total is below \(2q\) the subtracted floor is at most one, and when it is at least \(2q\) each pair sums to more than \(q\). There are at most \(e-1\) nonzero contributions since \(m<p^e\).

At a common zero, apply (1) at \(m=n\). The factor \(2(2n+1)/n\) has \(p\)-valuation \(e\), so \(W_{n+1}\) reduces to zero. By (2), \(W_{n+2},W_{n+3}\) also reduce to zero. The preceding local argument makes all three residues nonzero. The grouped Vandermonde test excludes two or three distinct residues for every odd \(p\), since the multiplicities are then at most two and any two distinct residues cannot both equal \(s/2\). For a single nonzero residue \(u\), the weight \(3u\) excludes \(p\ge5\). It remains to consider \(p=3\) and three equal residues.

### The equal residue at 3

Here \(j=n+2\) is divisible by 3. After scaling let \(a=1, b=1+h,c=1+k\) with \(h,k\) in the maximal ideal. Use valuation normalized by \(\nu(3)=1\). Let \(L=h+k\), \(Q=(h-k)^2\), and \(G=W_{j+1}-W_j\). Since \(G=0\) is the sum of \((1-h+k)(1+h)^j h\) and \((1+h-k)(1+k)^j k\), the two valuations \(\nu(h),\nu(k)\) agree, say at \(t>0\). Expanding gives
\[
 \begin{aligned}
 G&=L-Q+j(h^2+k^2)+O((h,k)^3),\\
 W_j&=3+(j+1)L+{j\choose2}(h^2+k^2)-jQ+O((h,k)^3).
 \end{aligned}
\]
Here the remainders are integral polynomials of total degree at least three. The first equality and \(G=0\) initially give \(\nu(L)\ge2t\), so \(\nu(h-k)=t\) from \(h-k=2h-L\). Since \(3\mid j\), the same equality now gives \(\nu(L)=2t\). Also \(3\mid {j\choose2}\). The second equality and \(W_j=0\) give \(2t=1\) by comparing the constant, linear and higher terms.

Adjoin \(\pi\) with \(\pi^2=-3\) at the same place, and write \(h=\pi H,k=\pi K\). The two equations
\[
     G/\pi=0,\qquad (W_j-(j+1)G)/\pi^2=0
\]
are polynomials in the integral variables \(H,K\) with integral coefficients. Their reductions are \(H+K=0\) and \((H-K)^2-1=0\). At either solution the reduced Jacobian determinant is \(-4(H-K)\), a unit. There is at most one integral solution in each residue class: the difference of the two polynomial values at two proposed solutions factors as a matrix times the difference of the arguments, with the matrix reducing to this invertible Jacobian. The two explicit solutions are the two orderings \(b,c=(-1\pm\pi)/2\) of the primitive cube roots of unity. They satisfy both equations because \(3\mid j\), and their reductions for \(H,K\) have \(H-K=\pm1\). Thus this branch has \(\Delta=0\), already excluded for a common zero of three consecutive moments.

## Shifted prime-power windows

For polynomials in \(r\), define
\[
\begin{aligned}
 D&=-373248r^5(12r+1)^2,\
 A&=161243136r^5+42456960r^4+3203712r^3+51624r^2-612r+1,\
 B&=-216(559872r^3+201528r^2+23412r+881),\
 C&=-3888(5832r^2+540r-1),\
 U&=A^2-12BD,\qquad J=CD^2,\qquad K=A^3,\qquad
 V=K^{10}U(J/K)\in\mathbb Z[r].
\end{aligned}                                                       \tag{4}
\]

**Proposition (Shifted prime-power windows).** Let \(p\ge5\) be prime, \(e\ge1\), and \(2n-1=p^e\). If the reductions of \(U,V\) in (4) are coprime over \(\mathbb F_p\), the three moments \(P_n,P_{n+1},P_{n+2}\) have no common zero with \(abc\ne0\). This test passes at \(5,7,17\) and at all but finitely many primes.

Suppose there is a common zero. Work over a number field and reduce at a place above \(p\) as in the preceding proposition. Put \(h=ab+bc+ca\), \(r=abc\), \(Q=p^e\), and \(j=n+2=(Q+5)/2\). All three coordinates are units by the preceding valuation test for \(W_j=W_{j+1}=0\). The three residues are distinct: if there are one or two distinct ones, the grouped Vandermonde test applied to these two equations makes every group weight zero. The multiplicities are units and this is impossible, with weight \(3t\) in the single group. Moreover \(W_{n+1}\) is a unit by the three-row Vandermonde test; if it reduced to zero, the three weights on the distinct residues would vanish.

Apply (1) at \(m=n-1\) and \(m=n\). At \(m=n-1\) the coefficient \(2(2m+1)/m\) has valuation \(e\), and \(P_{n-2}\) has valuation at least \(-(e-1)\) by the factorial floor bound, since \(n-2<Q\). At \(m=n\), that coefficient is a unit with residue \(8\), since \(n\equiv1/2\pmod p\). The equality at \(m=n\) makes \(P_{n-1}\) integral at this place. Reduction of the two equalities gives
\[
       \Delta P_{n-1}=W_n,\qquad 8r P_{n-1}+W_{n+1}=0
\]
in the residue field. The three bases of \(W\) are roots of \(T^3-sT^2+hT-r\), so \(W_{n+3}=sW_{n+2}-hW_{n+1}+rW_n\). The two consecutive exact zeros give \(rW_n=hW_{n+1}\). Multiplying the second reduced equality by \(\Delta\) now gives \((\Delta+8h)W_{n+1}=0\). As \(W_{n+1}\) is a unit, the residue point lies on \(s^2=12h\).

If \(s=0\) in the residue field, then \(h=0\). After scaling the three entries are \(1,\omega,\omega^2\) for a primitive cube root \(\omega\). The two equations require \(3\mid j\), because \(W_k=-2(1+\omega^{k+1}+\omega^{2(k+1)})\) there. For \(3\mid j\), these same cube roots give characteristic-zero solutions of the two \(W\) equations. Consider the equations with \(a=1\) and variables \(b,c\). At the reduced cube-root pair, the row of derivatives of \(W_j\) is \((-4,-4)\), since each entry is \(1-2j\) and \(j\equiv5/2\bmod p\). The row for \(W_{j+1}\) is \((-9b,-9c)\), since its entries are \(-2(j+2)b^{j+1},-2(j+2)c^{j+1}\). The determinant \(36(c-b)\) is a unit. In a common finite extension containing the cube roots, the factorization of polynomial differences by the coordinate differences gives uniqueness in each residue class, as in the preceding \(p=3\) argument. Thus the original characteristic-zero point itself is a scaled cube-root triple, already excluded by \(\Delta=0\).

It remains to scale the residue point to \(s=1\), \(h=1/12\), and use \(r\) for its rescaled product. The reduced version of (3) holds with \(\lambda\ne0\) by the two-row nullspace on the three distinct residues. In particular all \(1-2a,1-2b,1-2c\) are nonzero. Squaring (3) at the first coordinate gives
\[
 a^Q=\lambda^2 w(a),\qquad
 w(x)=\frac{2+6x-9x^2}{3x^5(1-2x)^2}.
\]
Here \((b-c)^2=(1-a)^2-4(h-a+a^2)=2/3+2a-3a^2\). The other two coordinates obey the same formula. Thus \(S_1=\sum_{x=a,b,c}w(x)\) is nonzero, since \(\sum a^Q=1\). If \(S_2,S_3\) are the other two elementary symmetric functions of the three \(w\)'s, then
\[
              S_1^2-12S_2=0,\qquad r^Q=S_3/S_1^3 .
\]
The first follows from preservation of \(h=1/12\) under Frobenius and the second from preservation of the product.

The resultant in \(x\) of \(12x^3-12x^2+x-12r\) and \(3T x^5(1-2x)^2-(2+6x-9x^2)\) is \(1152(DT^3+AT^2+BT+C)\). Consequently \(S_1=-A/D,S_2=B/D,S_3=-C/D\), and both \(D\) and \(A\) are nonzero at the point under consideration. The two identities give \(U(r)=0\) and \(r^Q=J(r)/K(r)\). Frobenius also gives \(U(r^Q)=0\), as \(U\) has integer coefficients, so \(V(r)=0\). This contradicts the stated coprimality.

The [exact polynomial checker](verify_sfc_triangle_shifted_windows.py), run with `sage -python research/factorial/image-rigidity/verify_sfc_triangle_shifted_windows.py`, verifies the resultant identity and the coprimality modulo \(5,7,17\). For example \(V\bmod U\) over \(\mathbb F_5\) is \(-r^9+r^7+2r^6+r^4-r^2+2\), coprime to \(U\), whose degree is ten modulo \(5\). In particular \(U,V\) are coprime over \(\mathbb Q\), since a common primitive factor of \(U\) cannot lose degree modulo \(5\) when \(U\) retains its degree. Hence only finitely many primes can fail the test. These windows include infinitely many orders outside the previous prime-power windows: take \(p=5\) and \(e\) a positive multiple of six. Then \(2n+1=5^e+2>3\) has \(3\)-valuation one, so it is not a prime power.

### The second shift

For \(h\) an indeterminate, put \(\rho(h)=(240h^2-40h+3)/128\). In the free rank-three algebra over \(\mathbb Z[1/2][h]\) with generator \(x\) satisfying \(x^3-x^2+hx-\rho(h)=0\), let \(M_f\) denote multiplication by \(f\) in the basis \(1,x,x^2\). Define
\[
\begin{aligned}
 \widehat D T^3+\widehat A T^2+\widehat B T+\widehat C
   &=\det\bigl(T M_{x^7(1-2x)^2}- M_{1-4h+2x-3x^2}\bigr),\
 \widehat U
   &=3\widehat A^4-40\widehat B\widehat D\widehat A^2
     +240(\widehat B\widehat D)^2-128\widehat C\widehat D^2\widehat A,\
 \widehat V
   &=(\widehat A^2)^{53}
          \widehat U\bigl(\widehat B\widehat D/\widehat A^2\bigr).
\end{aligned}                                                       \tag{5}
\]
The degree of \(\widehat U\) is \(53\), so \(\widehat V\) is polynomial.

**Proposition (Second shifted prime-power window).** Let \(p\ge7\), \(p\ne11\), be prime and \(2n-3=p^e\), \(e\ge1\). If the reductions of \(\widehat U,\widehat V\) in (5) are coprime over \(\mathbb F_p\), the three moments starting at \(n\) have no common zero with \(abc\ne0\). This test passes at \(19,29,37,41\) and at all but finitely many primes.

Use the same reduction as for the first shift, now with \(j=(Q+7)/2\). The three coefficients and \(W_{n+1}\) are units, and the residues of the coefficients are distinct. The use of (1) at \(m=n,n-1,n-2\) gives the following reduced equations:
\[
 0=\tfrac{16}3 rP_{n-1}+W_{n+1},\qquad
 \Delta P_{n-1}=8rP_{n-2}+W_n,\qquad
 \Delta P_{n-2}=W_{n-1}.
\]
In detail, the coefficients \(2(2m+1)/m\) for \(m=n,n-1\) reduce to \(16/3,8\), so the two preceding moments are integral at the place in succession. At \(m=n-2\) the factor has valuation \(e\), and \(P_{n-3}\) has valuation at least \(-(e-1)\). The recurrence of \(W\) gives \(rW_n=hW_{n+1}\) and \(r^2W_{n-1}=(h^2-sr)W_{n+1}\). Eliminate the two preceding moments from the reduced equations and cancel the unit \(W_{n+1}\). The resulting relation is
\[
                            3s^4-40s^2h+240h^2-128sr=0.            \tag{6}
\]
If \(s=0\), this gives \(h=0\) and hence the cube-root triple in the residue field. The two \(W\) equations then require \(3\mid j\). Their reduced derivative rows with \(a=1\) are now \((-6,-6)\) and \((-11b,-11c)\) with distinct \(b,c\), giving a unit determinant at the specified primes. The same uniqueness argument forces an actual cube-root triple and excludes a common zero of the moments.

Otherwise set \(s=1\) by scaling in the residue field. Now \(r=\rho(h)\), the coefficients are the roots of the cubic used in (5), and the squared nullspace relation is \(a^Q=\lambda^2 w(a)\) and its two counterparts, where \(w(x)=(1-4h+2x-3x^2)/(x^7(1-2x)^2)\). The determinant in (5) is the product of the three \(T x^7(1-2x)^2-(1-4h+2x-3x^2)\) over the roots of that cubic. This follows by diagonalizing multiplication in the three root coordinates, which are distinct. Thus \(\widehat D\ne0\), and \(\widehat A\ne0\) since \(\sum a^Q=1\). The pair sum and product of the Frobenius-transformed coefficients are
\[
                  h^Q=\widehat B\widehat D/\widehat A^2,\qquad
                  r^Q=\widehat C\widehat D^2/\widehat A^3 .
\]
They satisfy (6) with \(s=1\). Clearing denominators gives \(\widehat U(h)=0\). Since also \(\widehat U(h^Q)=0\), we get \(\widehat V(h)=0\), contrary to coprimality. The same exact checker verifies degree \(53\) and coprimality at the four listed primes. Since the degree survives modulo \(19\), the certificate at \(19\) implies coprimality over \(\mathbb Q\) and thus at all but finitely many primes.

## A residue equation for any fixed shift

Write \(h=ab+bc+ca\), \(r=abc\), \(s=a+b+c\) as above, and let \(t\ge1\) be an integer shift. Define polynomials \(B_{-2}=B_{-1}=0\), \(B_0=1\), and
\[
 B_i=h B_{i-1}-sr B_{i-2}+r^2 B_{i-3}\qquad(i\ge1).
\]
Put \(c_i=8i/(2i-1)\) for \(i\ge1\). The shifted residue equation is
\[
                  \Phi_t(s,h,r)=
          \sum_{i=0}^t \Delta^{t-i}B_i\prod_{j=t-i+1}^{t}c_j=0,    \tag{7}
\]
with an empty product equal to \(1\).

**Proposition (General shifted residue obstruction).** Suppose a three-moment zero exists at \(n\) with \(abc\ne0\). For any \(2n+1-2t=Q=p^e\) with \(p\ge5\), \(p>2t\), \(e\ge1\), its normalized reduction at a place above \(p\) lies on (7). It has three distinct unit coordinates, and its \(Q\)-Frobenius in projective coordinates is given by
\[
 [a^Q:b^Q:c^Q]=
 \left[
 \frac{(b-c)^2}{a^{2t+3}(s-2a)^2}:
 \frac{(c-a)^2}{b^{2t+3}(s-2b)^2}:
 \frac{(a-b)^2}{c^{2t+3}(s-2c)^2}
 \right].                                                        \tag{8}
\]
The curve (7) therefore contains the entire forward orbit of the residue point under this rational map.

Apply the unit and distinct-residue arguments from the first shift to \(W_{n+2}=W_{n+3}=0\). They also give that \(W_{n+1}\) is a unit. Put \(m_0=n-t=(Q-1)/2\). At \(m=m_0+i\) for \(1\le i\le t\), the coefficient \(2(2m+1)/m\) is a unit with residue \(c_i\). Work backwards from \(P_n=0\) in (1) to see that \(P_{m_0},\ldots,P_n\) are integral at the place. At \(m=m_0\) that coefficient has valuation \(e\), and \(P_{m_0-1}\) has valuation at least \(-(e-1)\) by the factorial floor bound. Thus \(\Delta P_{m_0}=W_{m_0+1}\) in the residue field. Substitute this in the next \(t\) steps, clearing powers of \(\Delta\) by multiplication. Using the recurrence of \(W\) with characteristic polynomial \(T^3-sT^2+hT-r\), one has \(B_i=r^i W_{n+1-i}/W_{n+1}\) in the residue field. After canceling the unit \(W_{n+1}\), the expanded last step is (7). Squaring (3) gives (8) because \(2(n+2)=Q+2t+3\), and the reduced weights \(s-2a,s-2b,s-2c\) are nonzero by the distinct-coordinate nullspace formula. Both the map and (7) are defined over the prime field. This proves the orbit condition.

### The third shift at seven

**Theorem (Third-shift windows at seven).** The three moments starting at \(n\) have no common zero with \(abc\ne0\) if \(2n-5=7^e\), \(e\ge1\).

For \(t=3\), (7) multiplied by \(5\) reads
\[
        5\Delta^3+24h\Delta^2+
        128(h^2-sr)\Delta+1024(h^3-2shr+r^2)=0 .
\]
At \(p=7\), its value on \(s=0\) is \(1024r^2\ne0\) because all three coordinates are units. We may therefore normalize to \(s=1\) in the residue field. The reduced relation becomes \(hr+r^2+r+1=0\), so \(h=-r-1-1/r\).

For indeterminate \(r\), put \(h=-r-1-1/r\) in \(\mathbb F_7(r)\). In the cubic algebra given by \(x^3-x^2+hx-r=0\), form the multiplication determinant
\[
       D T^3+A T^2+B T+C
        =\det\bigl(T M_{x^9(1-2x)^2}-M_{1-4h+2x-3x^2}\bigr).
\]
At a possible residue point all roots of the cubic are distinct, the denominators \(x^9(1-2x)^2\) are units, and the sum of the three rational values in (8) is nonzero. Consequently \(D,A\ne0\) there. Just as in the second shift, the two symmetric coordinates of the Frobenius image in the chart \(s=1\) are
\[
                          h'=BD/A^2,\qquad r'=CD^2/A^3.
\]
This image lies on the same curve, so the numerator of \(h'r'+(r')^2+r'+1\) vanishes at \(r\).

Define \(G\in\mathbb F_7[r]\) as the monic numerator of that rational function, with all factors supported on the zeros of
\[
 r\cdot\operatorname{num}(h^2-4h^3-4r-27r^2+18hr)\cdot
               \operatorname{num}(A)\operatorname{num}(D)
               \operatorname{den}(r')
\]
removed. Here numerator and denominator are taken in lowest terms. An admissible \(r\) and its Frobenius image must both be roots of \(G\); the removed factors are excluded by the preceding units and distinctness conditions. The [exact checker](verify_sfc_triangle_shifted_windows.py) computes \(\deg G=144\) and verifies that
\[
                       \gcd\bigl(G,\operatorname{den}(r')^{144}G(r')\bigr)=1
                       \quad\text{in }\mathbb F_7[r].
\]
This excludes the required compatibility \(r^Q=r'\) and proves the theorem for every \(e\).

The third shift includes infinitely many starts outside the three preceding prime-power families. For \(e\equiv60\bmod12180\), the integers \(7^e+6,7^e+4,7^e+2\) have valuations \(1\) at \(29,5,3\), respectively, and exceed those primes. Indeed \(7^{60}\equiv226\bmod29^2\), and \(7^e\equiv1\) modulo \(25\) and \(9\) for these exponents. They cannot be prime powers.

## Frobenius on the individual roots

The symmetric tests for the first two shifts use \(h^Q\) or \(r^Q\) and the residue curve. When a candidate survives those tests, equation (8) still prescribes Frobenius on each separate root of the cubic.

**Proposition (Individual-root refinement of the first two shifts).** In addition to the earlier prime tests, the window \(2n-1=p^e\) is covered for every \(e\ge1\) at \(p=11,29,59,83,359\). The window \(2n-3=p^e\) is covered for every \(e\ge1\) at \(p=13,23,43,53,79,89,293,313\). The same method gives the exclusions in the last column of the following table.

| Shift \(t\) | Prime \(p\) | Exponents still permitted by the residue test |
| ---: | ---: | --- |
| 1 | 13 | \(e\equiv1\pmod3\) |
| 1 | 281 | \(e\equiv1\pmod2\) |
| 2 | 7 or 73 | \(e\equiv0\pmod3\) |
| 2 | 17 | \(e\equiv1,4,5\pmod6\) |
| 2 | 211 | \(e\equiv0\pmod2\) |

The cases \(s=0\) are excluded as in the corresponding shifted proofs. Set \(s=1\), and factor \(\gcd(U,V)\) from (4) for the first shift or \(\gcd(\widehat U,\widehat V)\) from (5) for the second. For an irreducible factor of degree \(d\), use its root as the coordinate \(r\) or \(h\) over \(K=\mathbb F_{p^d}\). The other coordinate follows from the residue curve. Put
\[
 f(x)=x^3-x^2+hx-r,\qquad
 w(x)=\frac{1-4h+2x-3x^2}{x^{2t+3}(1-2x)^2}.
\]
A possible common zero makes \(f\) separable, its three roots units, and all three denominators nonzero. In \(K[x]/(f)\), let \(S_1\) be the trace of multiplication by \(w\); it is the sum over the three roots and is nonzero. The normalized point must satisfy the full identity
\[
                             x^{p^e}=w(x)/S_1\pmod f.              \tag{9}
\]
In fact (8) gives a common scalar and \(\sum a^Q=1\) determines it. The quotient algebra is a product of finite fields. If \(m\) is the least common multiple of the degrees of the irreducible factors of \(f\) over \(K\), then \(x^{p^e}\bmod f\) is periodic in \(e\) with period dividing \(dm\). Thus (9) is checked for a complete set of exponents, for one generator of each parameter factor; its conjugates satisfy the same tests.

The [exact checker](verify_sfc_triangle_shifted_windows.py) uses \(p\)-powering in these quotient algebras. It finds no exponent for the two lists of primes in the proposition, and gives the displayed residues for the other cases. For instance the new second-shift family \(p=13,\ e\equiv0\pmod{60}\) lies outside all the earlier certified windows. In this family \(13^e+4\) and \(13^e+2\) have valuations \(1\) at \(5\) and \(3\), respectively, and \(13^e-2\) is a unit at \(7\). The earlier second-shift symmetric gcd test at \(13\) has the nonconstant factor \(h+1\).

## The fourth shift

**Theorem (Fourth-shift windows).** For \(2n-7=p^e\) the three consecutive factorial moments with \(abc\ne0\) have no common zero at \(p=17\) for all \(e\ge1\), at \(p=11\) for \(3\nmid e\), and at \(p=13\) for odd \(e\).

Use the general offset obstruction with \(t=4\), so \(2(n+2)=Q+11\). For a cubic \(x^3-sx^2+hx-r\), let \(D,A,B,C\) be the four coefficients in decreasing order in \(T\) of
\[
 \det\!\bigl(T M_{x^{11}(s-2x)^2}-M_{s^2-4h+2sx-3x^2}\bigr).
\]
They represent the three symmetric functions \(-A/D,B/D,-C/D\) of the three values on the right of (8). In the chart \(s=0\), normalize to \(r=1\) over the algebraic closure. The unit and separability conditions imply \(D C(-4h^3-27)\ne0\). Frobenius preserves \(s=0\), so \(A=0\). At each of the three primes the equations \(\Phi_4(0,h,1)=0\) and \(A=0\) have no roots after removing \(h=0\) and the zeros of that product, by the exact univariate certificate. At \(h=0\) the residue triple consists of cube roots after projective scaling and forces \(3\mid j=n+2\). In the chart \(a=1\), the derivative rows of \(W_j,W_{j+1}\) are \((-10,-10)\) and \((-15b,-15c)\), since \(j\equiv11/2\bmod p\). Their determinant is a unit at the three primes, so the same residue class contains only the characteristic-zero cube-root solution to the two \(W\) equations. That solution has \(\Delta=0\), already excluded.

If \(s\ne0\), set \(s=1\). The two equations
\[
                 \Phi_4(1,h,r)=0,\qquad \Phi_4(A,BD,CD^2)=0
\]
are necessary by (7) and the symmetric coordinates of the image in (8). Here \(A D r\operatorname{disc}(x^3-x^2+hx-r)\ne0\). Compute the univariate eliminant in \(r\) after saturation by this product. For each irreducible factor of the eliminant, take \(r\) to be its root in the corresponding finite field, take the gcd of the specialized two equations in \(h\), and remove the zeros of the excluded factors. This leaves a finite list of candidate parameters. The scalar coordinate \(r'=CD^2/A^3\) must satisfy the minimal polynomial of \(r\), because it is the Frobenius image of \(r\). Test this on each candidate before checking (9) with \(t=4\).

The [exact checker](verify_sfc_triangle_shifted_windows.py) obtains eliminant degrees \(298,293,298\) at \(11,13,17\), respectively, with no admissible period at \(17\). At \(11\) only \(e\equiv0\bmod3\) survives the root test, and at \(13\) only \(e\equiv0\bmod2\) survives. All surviving parameter tests have linear fibers in \(h\) over the field generated by \(r\); the checker verifies this before running the cubic-root tests. This gives the stated windows. For instance \(p=17\) and \(e\equiv0\bmod1260\) give infinitely many starts outside all the earlier windows. The congruences \(17^e\equiv1\bmod(27\cdot49\cdot25)\) force \(17^e+8,17^e+6,17^e+4\) to have valuations exactly \(2,1,1\) at \(3,7,5\), respectively, and each exceeds the corresponding prime power. Also \(17^e+2\not\equiv0\bmod7\).

### The fifth shift

**Theorem (Fifth-shift windows).** For \(2n-9=p^e\), the three consecutive factorial moments with \(abc\ne0\) have no common zero at \(p=11\) for all \(e\ge1\), and at \(p=13\) for \(4\nmid e\).

Use (7) and (8) with \(t=5\), and put \(j=n+2=(Q+13)/2\). In the chart \(s=1\), use the same two equations and saturation as in the fourth shift with \(\Phi_5\) in place of \(\Phi_4\), and \(x^{13}(s-2x)^2\) in the multiplication determinant. The eliminant degrees are \(533,532\) at \(p=11,13\), respectively. There are no compatible periods at \(11\). At \(13\), the cubic-root test on the remaining fibers permits only \(e\equiv0\pmod4\).

In the chart \(s=0\), normalize to \(r=1\). For \(h=0\) the same cube-root argument applies: the two derivative rows are \((-12,-12)\) and \((-17b,-17c)\) at the reduced cube roots, so the determinant is a unit at both primes. For \(h\ne0\) the equations \(\Phi_5(0,h,1)=0\), \(A=0\) have no roots after the unit and separability exclusions at \(11\). At \(13\) the remaining parameter polynomial is supported on \(h^3-6=0\). In this chart the root test in \(K[x]/(x^3+hx-1)\) checks whether \(x^{p^e}/w(x)\) is constant, as required by (8). The same finite-field period argument applies. The checker finds no compatible exponent at the remaining parameters.

For example \(p=11\) and \(e\equiv10080\bmod318780\) give starts outside all the preceding certified windows. These exponents are multiples of \(1260\), so \(11^e\equiv1\bmod(27\cdot49\cdot25)\). Also \(11^e\equiv404\bmod23^2\). Therefore the five integers \(11^e+10,11^e+8,11^e+6,11^e+4,11^e+2\) have valuations \(1,2,1,1,1\) at \(23,3,7,5,3\), respectively, and exceed the corresponding prime powers.

## Remaining sector

For \(D=1\) there remains the case \(\Delta abc\ne0\) of three distinct coefficients with \(s-2t\ne0\) at every coefficient, at orders outside the proved windows. On that locus (3) together with \(P_n=0\) is equivalent to three consecutive moment zeros, by applying (1) forwards twice. These two power equations and the one moment equation must hold on the same projective point. The proof gives the four-moment result for every order and three-moment results in the initial and shifted prime-power windows. It gives no recurrence for \(D>1\). The [separating-row theorem](SFC_SEPARATING_ROW_AMPLIFICATION.md#3-an-intrinsic-zero-pattern-theorem) handles \(D\ge36\) separately. The original strong condition is the support-sized window in [Edo–van den Essen, Remark 2.9](https://arxiv.org/pdf/1304.3956).

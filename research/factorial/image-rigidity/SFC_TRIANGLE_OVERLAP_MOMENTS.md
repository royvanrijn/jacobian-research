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

## Remaining sector

For \(D=1\) there remains the case \(\Delta abc\ne0\) of three distinct coefficients with \(s-2t\ne0\) at every coefficient, at orders outside the proved windows. On that locus (3) together with \(P_n=0\) is equivalent to three consecutive moment zeros, by applying (1) forwards twice. These two power equations and the one moment equation must hold on the same projective point. The proof gives the four-moment result for every order and a three-moment result for an infinite set of orders. It gives no recurrence for \(D>1\). The [separating-row theorem](SFC_SEPARATING_ROW_AMPLIFICATION.md#3-an-intrinsic-zero-pattern-theorem) handles \(D\ge36\) separately. The original strong condition is the support-sized window in [Edo–van den Essen, Remark 2.9](https://arxiv.org/pdf/1304.3956).

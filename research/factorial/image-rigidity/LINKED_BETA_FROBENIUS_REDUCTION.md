# Frobenius reduction for linked beta moment pairs

For positive integer exponents \(a,b\), put \(d=a+b+1\). Consider the
ordinary factorial moments of \(X^aY^b(X+zY)\), with consecutive
orders \(N-1,N\). Suppose \(N=h p^q\), where \(p>hd\) is prime and
\(q\ge1\). An exact modular coprimality test needs only a fixed
polynomial of degree \(a+b\), and an auxiliary polynomial of degree
\(h\). In particular, for each fixed positive exponent pair, a
finite explicit set of primes is the only exception to the test at
prime-power upper orders. One test at a prime covers all of its
positive powers under the sufficient criterion in Section 4.

## 1. Polynomials and the reduction

After removal of nonzero scalar factors, the native moment
polynomial is
\[
 B_m(z)=\sum_{k=0}^m
   \binom{m(a+1)-k}{ma}\binom{mb+k}{mb}\,z^k .
\]
Indeed the two binomial factors are obtained by dividing the
ordinary factorial coefficient by \(m!(ma)!(mb)!\).
Write \(Q=p^q\), \(n=hQ-1\), and define
\[
\begin{aligned}
 M(Y)&=\sum_{j=0}^h
  \binom{h(a+1)-j-1}{ha-1}\binom{hb+j-1}{hb-1}Y^j,\
 W(Y)&=Y M'(Y),\
 A(Y)&=\frac{h(a+1)M(Y)-W(Y)}{ha},&
 C(Y)&=\frac{hbM(Y)+W(Y)}{hb},\\
 L_b(z)&=\sum_{r=0}^{b-1}(-1)^r\binom{b-1}{r}
                  \frac{a!r!}{(a+r)!}z^r,&
 L_a(z)&=\sum_{r=0}^{a-1}(-1)^r\binom{a-1}{r}
                  \frac{b!r!}{(b+r)!}z^r,\\
 D(z)&=(b+1)z-b,&U(z)&=a+1-az,\\
 H(z)&=z^a D(z)L_b(z)+z^{a-1}U(z)L_a(1/z).
\end{aligned}
\]
All displayed denominators are units modulo \(p\). The degree of
\(H\) over both \(\mathbb Q\) and \(\mathbb F_p\) is \(a+b\), because
its leading coefficient comes from the first summand and is
\((b+1)(-1)^{b-1}a!(b-1)!/(a+b-1)!\). Its constant coefficient is
\((a+1)(-1)^{a-1}b!(a-1)!/(a+b-1)!\), also a unit.

The polynomial \(B_n/Q\) has integer coefficients. The monic gcd of
the reductions of \(B_n/Q\) and \(B_N\) is
\[
 \gcd_{\mathbb F_p[z]}\bigl(H(z),A(z^Q)-zC(z^Q)\bigr).
 \tag{1}
\]
The large exponent in (1) can be computed modulo \(H\) by successive
Frobenius powers. Both reduced moment polynomials retain their
degrees, so a unit gcd in (1) proves characteristic-zero
coprimality by reduction of the Sylvester resultant.

## 2. The two coefficient reductions

Use the base-\(p\) binomial coefficient formula from
[Lucas, *Théorie des fonctions numériques simplement périodiques*](https://www.jstor.org/stable/2369308).
Its reduction follows by expanding \((1+T)^p=1+T^p\) in
characteristic \(p\) through the digits. For the upper order, let
\(U_j'=\binom{h(a+1)-j}{ha}\) and \(V_j'=\binom{hb+j}{hb}\).
The coefficient at \(k=jQ\) reduces to \(U_j'V_j'\), and each
coefficient at \(k=jQ+r\), \(1\le r<Q\), reduces to
\(U_{j+1}'V_j'\). Multiply the sum of geometric blocks by \(1-z\)
and put \(Y=z^Q\). Using \(V_{-1}'=U_{h+1}'=0\) gives
\[
\begin{aligned}
 (1-z)B_N(z)&\equiv A(Y)-z C(Y),\
 [Y^j]A&=U_j'(V_j'-V_{j-1}'),\
 [Y^j]C&=V_j'(U_j'-U_{j+1}').
\end{aligned} \tag{2}
\]
Pascal's identity in the last two lines gives the stated expressions
in terms of \(M,W\).

For the lower order, the contribution to the valuation of
\([z^k]B_n\) at each factorial-floor level \(x=p^i\), \(1\le i\le q\),
is
\[
 \mathbf 1_{k\bmod x<x-a}+\mathbf 1_{k\bmod x\ge b}.
 \tag{3}
\]
For the first factor, the two lower factorial arguments have
residues \(x-a\) and \(x-1-(k\bmod x)\). For the second they have
residues \(x-b\) and \(k\bmod x\). A level contributes a carry
precisely when the residues sum to at least \(x\), by the formula
\(v_p(t!)=\sum_{i\ge1}\lfloor t/p^i\rfloor\). Higher levels give
zero because the largest factorial argument is less than \(hdQ<pQ\).
Since \(x>a+b\), (3) is at least one. At level \(Q\) it equals
one exactly in the two bands \(0\le k\bmod Q<b\) and
\(Q-a\le k\bmod Q<Q\). Those two bands retain a sum of one at
each smaller level. Thus only these bands survive modulo \(p\)
after division by \(Q\).

For the first band, begin at \(k=jQ\), \(0\le j<h\). The first
binomial factor has valuation \(q\) and its unit after division by
\(Q\) is
\[
 -\frac{1}{a}\frac{(h(a+1)-j-1)!}{(ha-1)!(h-j-1)!}.
 \tag{4}
\]
For a verification, the factorial unit
\(F(t)=t!/p^{v_p(t!)}\) satisfies
\(F(t)\equiv(-1)^{\lfloor t/p\rfloor}(t\bmod p)!
F(\lfloor t/p\rfloor)\) by taking full blocks of nonmultiples
of \(p\). Apply the formula through \(q\) levels. The first digit
of the factorial ratio is
\((p-a-1)!/((p-a)!(p-1)!)=1/a\) modulo \(p\).
The next \(q-1\) digits each contribute \(-1\), and the floor
signs contribute \((-1)^q\) by (3). The remaining factorials are
those of (4). The second binomial factor at \(k=jQ\) reduces to
\(\binom{hb+j-1}{hb-1}\) by its digits. Combining yields
\(- (h-j)[Y^j]M/a\).

On advancing \(r<b\) positions in this band the coefficient ratio
is the product of the consecutive ratios
\[
 \frac{[z^{k+1}]B_n}{[z^k]B_n}
   =\frac{(n-k)(nb+k+1)}{(n(a+1)-k)(k+1)}.
\]
At \(k=jQ+i\), \(0\le i<r\), the denominator is a unit, and the
ratio reduces to \(-(b-1-i)/(a+1+i)\). This gives exactly the
coefficients of \(L_b\). Reflection \(k\mapsto n-k\) interchanges
\(a,b\) and deals with the other band. The starting index of that
reflected band in the original variable is \((j+1)Q-1\); its
unit is \(-(j+1)[Y^{j+1}]M/b\). Consequently the complete reduction
is
\[
  \frac{B_n(z)}{Q}\equiv-\left(
    \frac{hM(Y)-W(Y)}{a} L_b(z)
       +\frac{W(Y)}{b z}L_a(1/z)\right),\qquad Y=z^Q.
 \tag{5}
\]
The right side is polynomial after substitution since the negative
powers are multiplied by a factor \(W(Y)\) divisible by \(Y\).
The expressions in (2) and (4) give unit endpoint coefficients for
both reductions. In particular their degrees are unchanged.

## 3. Elimination over the residue field

The upper reduction is nonzero at \(z=0,1\). At zero this follows
from its endpoint coefficient, and at one Vandermonde's identity
gives \(B_N(1)=\binom{N d+1}{N}\equiv\binom{hd}{h}\ne0\).
Also \(H(1)=L_b(1)+L_a(1)=(a+b)/(a+b-1)\) because the beta integral
\[
 L_b(z)=a\int_0^1(1-t)^{a-1}(1-zt)^{b-1}\,dt
\]
gives \(L_b(1)=a/(a+b-1)\), with its reflected analogue for \(L_a\).

The values of \(M,W\) cannot both vanish at \(Y\ne0\). To see this,
the exact consecutive coefficient ratio of \(M\) is
\((h-j)(hb+j)/((j+1)(h(a+1)-j-1))\). It gives the equation
\[
 Y(1-Y) M''+[1-h(a+1)-(1-h+hb)Y]M'+h^2 b M=0.
\]
The argument \(Y=1\) gives the unit \(M(1)=\binom{hd-1}{h}\).
At any other nonzero \(Y\), a double zero of \(M\) would force
all derivatives through degree \(h<p\) to vanish by the equation
and its derivatives, in conflict with the leading coefficient.

At \(z\ne0,1\), equations (2) and (5) give the homogeneous linear
equations in \((M,W)\)
\[
 (az+b)W-hb U(z)M=0,\qquad
 \frac{z L_b(z)}{a}(hM-W)+\frac{L_a(1/z)}{b}W=0.
\]
The determinant is \(-h(zD(z)L_b(z)+U(z)L_a(1/z))\).
The first row is nonzero, because \(az+b\) and \(U(z)\) have
distinct zeros when \(p>d\). The pair \((M,W)\) is nonzero.
Thus on the zero set of the first equation, the second holds
exactly when the determinant vanishes. The factor \(z^{a-1}\)
converts that determinant to a nonzero scalar times \(H\).
The cases \(z=0,1\) already excluded, both gcds have the same roots.
Moreover, the derivative of \(A(z^Q)-zC(z^Q)\) is \(-C(z^Q)\)
in characteristic \(p\). A multiple root would therefore give
\(A=C=0\). The linear map from \((M,W)\) to \((A,C)\) is
invertible because \(p>d\), and even at \(Y=0\) the value \(M(0)\)
is a unit. Hence \(A(z^Q)-zC(z^Q)\) is squarefree. So is the
upper reduction by (2), proving equality of the monic gcds in (1).

## 4. Prime-power upper orders

For \(h=1\), equation (2) becomes \(z^Q D(z)=U(z)\).
Its right-hand side divided by \(D\) is a Möbius map with fixed
points \(1\) and \(r_0=-(a+1)/(b+1)\). They are distinct and
nonzero modulo \(p\) when \(p>d\) and \(p\nmid d+1\). The
remaining case \(p=d+1\) is still covered by the exact test
(1). In the chart
\(x=(z-1)/(z-r_0)\), the Möbius map becomes multiplication by
\(-d\): it fixes zero and infinity, and its derivative at \(z=1\)
is \(-d\). The chart commutes with \(Q\)-th powers over the
residue field. Therefore every upper root other than \(r_0\)
satisfies \(x^Q=-d x\), with \(x\ne0\). If its degree over
\(\mathbb F_p\) is \(l\), the order of \(-d\) in
\(\mathbb F_p^\times\) is exactly \(l/\gcd(l,q)\): it is the
length of the orbit under the \(q\)-th iterate of Frobenius.
When that multiplicative order exceeds \(a+b\), no such root
can satisfy \(H(z)=0\).

The remaining possible root \(r_0\) satisfies the lower reduction
if and only if
\(L_b(r_0)+L_a(1/r_0)\equiv0\pmod p\), by (5) and
\(r_0^{Q-1}=1\). This rational number is positive, since both
polynomials have alternating coefficients and both arguments
are negative. Its denominator is a unit for \(p>d\). Thus a
sufficient condition covering every \(q\ge1\) is
\[
 p>d,\qquad
 p\nmid \operatorname{num}\bigl(L_b(r_0)+L_a(1/r_0)\bigr)
             \prod_{j=1}^{a+b}\bigl((-d)^j-1\bigr).              \tag{6}
\]
The excluded primes are explicitly finite for each exponent pair.
For any admissible cofactor \(h\), (1) also gives an exact smaller
test for primes or exponents where the sufficient condition (6)
is unavailable. A nonunit gcd in (1) records an intersection only
over the residue field.

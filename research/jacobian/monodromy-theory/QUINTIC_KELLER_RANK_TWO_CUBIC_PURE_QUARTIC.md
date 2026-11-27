# Pure quartic in a primitive cubic pencil

**Proposition QuinticRankTwoCubicPureQuarticRestriction.** Let
\(F=(P,Q,R)\) be a complex three-variable Keller map of degree five
with leading forms \(P_5=z^5,\ Q_5=z^2T,\ R_4=z^4\), where \(T\)
is a homogeneous cubic coprime to \(z\) and \(T(X,Y,1)\) is closed.
If \(F\) is noninvertible, then \(R_3(X,Y,0)\) is a nonzero
multiple of \(T(X,Y,0)\), and this binary cubic has a repeated
linear factor.

## Slice and its highest binary degrees

Use the
[pure-quartic formal slice](QUINTIC_KELLER_RANK_TWO_REPEATED_LINE_PENCIL.md#closure-of-the-two-remaining-large-primitive-degrees)
with \(\epsilon=\tau^{-1}\), \(R(\tau X,\tau Y,\tau Z)=\tau^4\)
and \(Z=1+O(\epsilon)\). Write \(\Pi=\tau^{-5}P\),
\(\mathcal Q=\tau^{-5}Q\) on the slice, and \(r=T(X,Y,1)\).
As there, \(\{\Pi,\mathcal Q\}=O(\epsilon^{11})\). Whenever
all lower coefficients of an adjusted \(\Pi\) have been removed,
its next coefficient through order ten belongs to \(\mathbb C[r]\).
A polynomial in \(r\) can be removed at order \(j\) by subtracting
the same polynomial in \(\mathcal Q\) times \(\epsilon^j\).
The adjustment preserves the bracket bound.

Put \(r_j=R_j(X,Y,1)\), \(p_j=P_j(X,Y,1)\), and
\(q_j=Q_j(X,Y,1)\). The order-one coefficient is
\(p_4-5r_3/4=\lambda r+c\).
After subtracting \(\epsilon(\lambda\mathcal Q+c)\), the order-two
coefficient is
\[
 p_3-\tfrac54r_2-cr_3-\tfrac5{32}r_3^2
                           +\lambda(r_3r/4-q_4)\in\mathbb C[r].
 \tag{1}
\]
Here the \(\lambda\)-terms from \(P_4,Q_5\) before expansion are
\(\epsilon\lambda(Z-Z^2)T(X,Y,Z)\). If \(\deg r_3=3\), comparison
in degree six in (1) makes its highest part \(D=R_3(X,Y,0)\)
satisfy \(D(\lambda T(X,Y,0)/4-5D/32)=c_1 T(X,Y,0)^2\).
The ratio of the two cubics is algebraic over \(\mathbb C\), so
\(D\) is proportional to \(T(X,Y,0)\). If \(\deg r_3=2\),
degree five first forces \(\lambda=0\) and then degree four
is impossible. It remains to handle \(\deg r_3\le1\).

We use the binary degree minus \(\epsilon\)-order to track terms.
For example, the contribution of degree \(i+b\) at order \(i\)
is the coefficient of \(\eta^b\) after the formal substitutions
\((X,Y,\epsilon)\mapsto(\eta X,\eta Y,\epsilon/\eta)\).
Under \(\deg r_3\le1\), the highest such part of \(Z\) has degree
minus order zero. Write \(L\) for the linear part of \(r_3\), \(K\)
for the quadratic part of \(r_2\), \(A=\epsilon L\), \(B=\epsilon^2K\),
and \(W=Z^2\) for this part. It satisfies \(W^2+AW+B=1\).
When \(L\ne0\), translate \(z\) to make \(r_3=L\). This does not
affect \(L,K\) or the top binary parts of \(r\) and \(q_4\).
Write \(r_* = T(X,Y,0)\) and \(q_* = Q_4(X,Y,0)\).
The highest part of \(\mathcal Q\), in excess three, is
\(r_* W+\epsilon q_*\). Its lower homogeneous terms have
strictly smaller degree minus order. These rules also apply to
powers of \(\mathcal Q\).

## Linear \(r_3\) and a nonzero first coefficient

Suppose \(L\ne0,\ \lambda\ne0\). Degree four in (1) gives
\(q_*=Lr_*/4\). The adjusted \(\Pi\) has highest possible
degree minus order two, with that part given by
\[
                   \epsilon\lambda r_* (Z-W-A/4).
\]
The coefficient of order three in this part has degree five and
is \(\lambda r_*(K/4-3L^2/32)\). At order two the possible
multiples of \(\mathcal Q\) removed have degree minus order one
and cannot affect this comparison. Degree five is forbidden,
so \(B=3A^2/8\). Now \(Z-W-A/4=-A^3/64-5A^4/1024+O(A^5)\).
At order four, the degree-six term must be a multiple of \(r_*^2\).
Its nonzero value forces \(r_*\) proportional to \(L^3\).
Removing this term by an order-four multiple of \(\mathcal Q^2\)
contributes \(-A^4/128\) inside the parenthesis at the next order,
because \((W+A/4)^2=1-A/2+O(A^2)\). The coefficient
\(-13/1024\) gives a nonzero term of degree seven at order five.
Other removals before order five have degree minus order below two,
so degree seven is impossible.

## Linear \(r_3\) and zero first coefficient

Now suppose \(L\ne0,\ \lambda=0\). Write the nonconstant part
of (1) as \(\gamma r\).
If \(\gamma\ne0\), removing it gives excess one part
\(\epsilon^2\gamma(r_*(1-W)-\epsilon q_*)\).
At order three degree four is forbidden and gives
\(q_*=Lr_*/2\). At order four degree five is forbidden, and with
this value of \(q_*\) the coefficient is
\(\gamma r_*(K/2-L^2/8)\). An order-three removal of a multiple
of \(\mathcal Q\) has excess zero and does not affect it.
Thus \(B=A^2/4\), \(W=1-A/2\), and the entire excess-one part
is zero.

To compute the excess-zero part in this case, write
\(R_2=K+zM+az^2\) with \(M\) linear binary, and
\(Q_4=(L/2)T+zU_3\), with \(U=U_3(X,Y,0)\).
The part of \(Z^2\) in excess minus one is
\(-\epsilon^2 M Z/2\), since \(2W+A=2\) and \(r_3=L\).
Consequently the excess-zero contribution, before removing
subsequent polynomials in \(\mathcal Q\), is
\[
      Z(1+AW/4+B/4+5A^2/32)
                    +\gamma Z(\epsilon^4 Mr_*/2-\epsilon^3U).
 \tag{2}
\]
In fact \(\epsilon^2\gamma T\) from \(\Pi\) and
\(-\epsilon^2\gamma\mathcal Q\) give
\(\epsilon^2\gamma T(1-W-A/2)\) at the top, plus
\(-\epsilon^3\gamma ZU\) below it. After the top vanishes,
the correction in \(Z^2\) gives the other term of (2).
The remaining terms directly use (1) for the quadratic part
of \(p_3\). The leading part of \(\mathcal Q\) is now \(r_*\)
at all orders in excess three. Thus removing any earlier
constant or multiple of a power of \(\mathcal Q\) leaves later
coefficients of (2) unchanged in excess zero. Over the fraction
field in \(X,Y\), write the last two terms of (2) as
\(Z(uA^3+vA^4)\). The vanishing of orders four and five
forces \(u=-11/256,\ v=-7/2048\). But its order-seven
coefficient is then \(-81/524288\), a forbidden degree-seven
term. Hence \(\gamma=0\).

With \(\gamma=0\) the highest possible degree minus order in
\(\Pi\), after the constant adjustments, is zero, and that
part is just \(Z(1+AW/4+B/4+5A^2/32)\). The cubic coefficient
is \(C_3=5L(L^2-8K)/128\). If this vanishes then \(B=A^2/8\)
and the degree-four coefficient \(-5 L^4/1024\) is forbidden.
If \(C_3\ne0\), it must equal \(b_3 r_*\), and the order-three
subtraction contributes \(-b_3\epsilon^3(r_* W+\epsilon q_*)\).
The order-five coefficient now is
\[
                   -L(960K^2-320KL^2+29L^4)/4096 .
\]
Its vanishing gives \(K=k L^2\) for a constant \(k\) satisfying
\(960k^2-320k+29=0\). This implies \(r_*\) proportional to \(L^3\).
The order-four equation fixes \(\epsilon q_*/r_*\) to be a constant
multiple of \(A\). If \(c_j\) are the coefficients of \(A^j\)
in \(Z(1+AW/4+kA^2/4+5A^2/32)\), this makes the top
\(\mathcal Q/r_*\) equal to \(W+(c_4/c_3+1/2)A\).
At order six subtract the resulting
multiple of \(\mathcal Q^2\) to remove the degree-six term.
The order-seven coefficient in \(A\) after both subtractions is
\[
             \frac{5(4k-1)(512k^3-704k^2+200k-23)}{131072}.
\]
The numerator and \(960k^2-320k+29\) are coprime.
Thus a term of degree seven remains. At the other intervening
orders the lower-degree subtractions lie strictly below
excess zero, and cannot change this term. This excludes \(L\ne0\).

## Constant \(r_3\) and the critical locus

If \(r_3\) is constant, translate \(z\) to make \(R_3=0\).
Suppose \(K\ne0\). Now \(Z=1-\epsilon^2r_2/4+O(\epsilon^3)\).
If \(\lambda\ne0\), (1) gives \(q_*=0\) and the order-three
degree-five term after the first subtraction is
\(\lambda K r_*/4\). The order-two subtraction changes it
only by a multiple of \(q_4\), of degree at most three in
this case. Hence \(\lambda=0\).
With \(\gamma\) the coefficient of \(r\) in (1), a nonzero
\(\gamma\) first forces \(q_*=0\) at order three, and then
gives the degree-five term \(\gamma K r_*/2\) at order four
from \(\epsilon^2\gamma(1-Z^2)T\). Since \(Z_1=0\), removing
a multiple of \(\mathcal Q\) at order three cannot cancel it.
It follows that \(\gamma=0\). The order-four contribution
of the excess-zero part \((1-B)^{1/4}(1+B/4)\) is now
\(-5K^2/32\), of degree four. The order-three coefficient
has binary degree at most two, from \(R_1\), the order-two
term of \(c\epsilon Z^4\), and \(P_2\). It is therefore
constant, and its removal has no effect on degree four.
This excludes \(K\ne0\).
For \(K=0\) the quartic \(R=f(z)+a(z)x+b(z)y\) has affine
\(a,b\), so the
[bounded coordinate completion](QUINTIC_KELLER_RANK_TWO_REPEATED_LINE_PENCIL.md#closure-of-the-two-remaining-large-primitive-degrees)
applies.

Finally, suppose the degree-three binary part
\(D=R_3(X,Y,0)\) is squarefree. The binary top parts of
\(R_x,R_y\) are the relatively prime quadrics \(D_x,D_y\).
On each plane \(z=c\) their common zero scheme is finite and
nonempty, by the plane projective intersection theorem and
absence of a common zero on its line at infinity. Take an
irreducible curve in \(R_x=R_y=0\) dominating the \(z\)-line,
and its normalized projective closure. Every point at
infinity has \(z\ne0\) in homogeneous coordinates because
a point with \(z=0\) there would be a common projective root
of \(D_x,D_y\). Thus on the normalization the affine
coordinate \(z\) has a pole at every point at infinity
and \(x/z,y/z\) are regular there. The function
\(R_z=4z^3+\partial_z R_{\le3}\) has a pole at each of these
points, since \(\partial_zR_{\le3}=O(z^2)\).
It therefore has a zero at some affine point of the curve,
giving a critical point of \(R\), a contradiction. This proves
the required repeated factor. The argument extends to
characteristic zero by coefficient-field embedding.

[The exact coefficient checker](verify_quintic_rank_two_cubic_pure_quartic.py)
verifies the rational series coefficients and the polynomial gcd.

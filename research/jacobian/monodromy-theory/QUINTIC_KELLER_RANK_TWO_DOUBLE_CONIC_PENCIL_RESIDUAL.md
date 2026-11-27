# The remaining rank-two double-conic pencil

**Corollary QuinticKellerRankTwoDoubleConicPureFifthPowerResidual.**
Let a characteristic-zero three-variable Keller map have independent
quintic leading parts and cancelling quartic leading form \(q^2\),
where \(q\) is a geometrically irreducible ternary conic. Over the
algebraic closure its quintic leading span must be
\[
                       \langle z^5,\ z^3q+\rho zq^2\rangle
\]
for a linear coordinate \(z\) and scalar \(\rho\). This is a
necessary condition for a Keller completion.

Use the [rank-two pencil representation and valuation identity](QUINTIC_KELLER_RANK_TWO_QUARTIC_PENCIL.md#valuations-along-a-projective-pencil),
with leading quartic \(C=q^2\), primitive degree \(m\), and
\(\deg g+me=5\). The member containing \(q\) can have no
other irreducible factor. Indeed the left side of (2) at \(q\) is
\(10-4\nu_q(g)>0\) since \(\deg g\le4\), and its value at any other
factor of this member is nonpositive. Thus \(m=2\) or \(4\).

If \(m=4\), the member is \(q^2\), and \(g=z\) is linear. The member
containing \(z\) is proportional to \(z^4\) by the same identity,
since every other factor of that member would have zero valuation
in both \(C\) and \(g\). This gives \(\langle z^5,zq^2\rangle\),
already excluded by the
[quartic-argument obstruction](QUINTIC_KELLER_RANK_TWO_DOUBLE_CONIC_QUARTIC_PENCIL.md).

For \(m=2\) the member containing \(q\) is \(q\) itself. Every factor
of \(g\) outside this member belongs to a pencil member. Its order
in \(g\) is proportional to the order in that member because \(C\)
is absent there. If the member is reduced, the proportionality is
integral and its contribution to the degree of \(g\) is even. The
contribution from \(q\) is even as well. Since \(\deg g\) is odd,
the pencil must have a doubled line \(z^2\), contributing an odd
power of \(z\). There is only one doubled line: two would make the
pencil binary and every member reducible, contrary to irreducibility
of \(q\).

If \(e=2\), then \(g=z\) and the quintic span comes from a pair of
coprime binary quadratics in \(z^2,q\) multiplied by \(z\). Some
member of the span has nonzero \(zq^2\) coefficient, since otherwise
the two binary quadratics would have a common \(z^2\) factor. If the
unique direction with zero \(zq^2\) coefficient has a nonzero
\(z^3q\) coefficient, normalize it to \(z^3q+\nu z^5\), and subtract
a multiple of it from the other member to get \(z(q^2+\delta z^4)\).
The [generic](QUINTIC_KELLER_RANK_TWO_DOUBLE_CONIC_GENERIC_PENCIL.md)
and [exceptional](QUINTIC_KELLER_RANK_TWO_DOUBLE_CONIC_EXCEPTIONAL_PENCILS.md)
determinant calculations force \(\delta=\nu=0\).
A transverse \(q\) is excluded by the
[Laurent restriction](QUINTIC_KELLER_RANK_TWO_DOUBLE_CONIC_LAURENT_RESIDUAL.md#exclusion-when-the-conic-meets-the-repeated-pencil-line-transversely),
and a tangent \(q\) by the
[tangent determinant obstruction](QUINTIC_KELLER_RANK_TWO_TANGENT_DOUBLE_CONIC.md).
These exhaust the cases because the restriction of an irreducible
ternary conic to \(z=0\) is nonzero. Thus the direction with no
\(zq^2\) term is \(z^5\). Its complement is
\(zq^2+\mu z^3q\) modulo \(z^5\), with \(\mu\ne0\) by the
quartic-argument obstruction. This yields the stated span.

If \(e=1\), then \(\deg g=3\). The contributions above give
\(g=z^3\) or \(g=z(q+\theta z^2)\) up to scalar. The first gives
\(\langle z^5,z^3q\rangle\). In the second, take
\(Q_5=z^3q+\theta z^5\) and, after subtracting a multiple of this
member, \(P_5=z(q^2-\theta^2 z^4)\). This has the excluded mixed
form, already considered for \(e=2\), with \(\nu=\theta\) and
\(\delta=-\theta^2\).

## Removing the nonzero coefficient

**Theorem QuinticKellerRankTwoDoubleConicPureFifthNonzeroCoefficientExclusion.**
Let \(q\) be a geometrically irreducible ternary conic. Over a field of
characteristic zero there is no Keller completion of
\[
                P_5=z^5,\qquad Q_5=z^3q+\rho zq^2,\qquad R_4=q^2
\]
with \(\rho\ne0\). Consequently the preceding double-conic residual
only has the span \(\langle z^5,z^3q\rangle\).

Work over \(\mathbb C\). Choose coordinates so that the coefficient
of \(x^2\) in \(q\) is nonzero, and normalize that coefficient to one.
The basis \(x^i y^j\), \(i=0,1,\ j\ge0\), is free over
\(\mathbb C[q,z]\) by monic division in \(x\), respecting homogeneous
degrees. At \(z=1\) put \(r=q(x,y,1)\), and use the corresponding basis
over \(\mathbb C[r]\).

Let \(p,t,s\) be the dehomogenizations of the correction
\((P_{5-n},Q_{5-n},R_{4-n})\). When every homogeneous piece of
smaller drop than \(n\) depends only on \(q,z\), the degree-\((11-n)\)
determinant equation gives
\[
       (6r^2+2\rho r^3)p-10r t+(5+10\rho r)s\in\mathbb C[r].
 \tag{3}
\]
In fact, in the coordinates of the [mixed cofactor proof](QUINTIC_KELLER_RANK_TWO_DOUBLE_CONIC_GENERIC_PENCIL.md#the-mixed-homogeneous-equation)
the Euler and \(q\)-derivative columns of the dehomogenized leading
forms are \((5,5(r+\rho r^2),4r^2)\) and \((0,1+2\rho r,2r)\).
Their cross product gives the row in (3). The determinant at this
degree uses at least one derivative tangent to the \(q,z\) levels.
All smaller-drop pieces have zero such derivative, so (3) is the
whole equation by the closed-polynomial criterion in
[Wang–Yan, Proposition 2.2 and Lemma 2.3](https://arxiv.org/pdf/2610.36495v1).
Here \(r\) is irreducible because its degree-two homogenization is
\(q\), and it is closed since a nontrivial polynomial composition
would be reducible over \(\mathbb C\).

For \(n=1\), reduce (3) modulo \(r\). It says that every positive-degree
basis coefficient of \(s\) is divisible by \(r\). Homogenizing in \(z\)
gives \(R_3=qL\) modulo \(\mathbb C[q,z]\), with \(L\) homogeneous
linear. Since the projective conic is smooth, the polarized bilinear
form of \(q\) is nonsingular. A translation of the source can therefore
induce linear correction \(-L/2\) to \(q\), removing this term from
\(R_3\) via \(R_4=q^2\). Such a translation fixes the leading forms,
so (3) still applies. Now for every positive-degree basis monomial
the coefficients \(p_*,t_*\) satisfy
\[
                         t_*=\frac{3r+\rho r^2}{5}\,p_* .
\]
The degree in \(r\) of \(t_*\) is at most one, since it comes from a
quartic and the basis monomial has positive degree. Because
\(\rho\ne0\), this forces both coefficients to be zero. Hence
\(P_4,Q_4,R_3\) depend only on \(q,z\).

For \(n=2,3,4\), reduce (3) modulo \(r\) again. Since \(R_{4-n}\)
has degree at most two, no positive-degree basis coefficient of \(s\)
can be a nonzero multiple of \(r\). Thus \(s\in\mathbb C[r]\).
The same degree comparison for \(p_*,t_*\) then makes them both zero.
This induction makes the full triple polynomial in \(q,z\), so its
Jacobian determinant vanishes. This contradicts the Keller condition.

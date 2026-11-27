# Rank exclusions for the binary quintic Veronese case

**Theorem QuinticVeroneseBinaryRankExclusions.** Under the hypotheses of
[QuinticVeroneseBinaryZReduction](QUINTIC_VERONESE_BINARY_Z_REDUCTION.md),
the coefficient \(A=[z^2]F\) has rank one, and its nonzero affine factor is nonconstant. Thus a linear target change and an affine binary source change put \(A=(0,0,x)\). The final section excludes this remaining form. The conclusions extend to characteristic-zero fields as in the preceding reduction.

## The rank-two coefficient

By the preceding reduction a rank-two coefficient can be put in the form \(A=(x,1,0)\), with \(B=(p,q,r)\) where \(r=r(x)\). Use \(z\mapsto z-q/2\) to complete the square in the second component. The resulting map is
\[
                  F=(xz^2+Pz+S,\ z^2+Q,\ r z+T),
 \qquad P=p-xq,\quad T=C_3-rq/2 .
\]
Here \(C_3\) denotes the third coordinate of the original constant term \(C\) as a polynomial in \(z\). Its top binary form has degree five because the leading parts have independent components. The \(z^3\) coefficient of the determinant is \(2(r'P_y-T_y)\), so \(T=r'P+d(x)\). Primes below are derivatives in \(x\). If \(r=0\), then \(C_3=T=d(x)\) has degree five and a critical point. This contradicts the nonzero determinant.

When \(P\) depends on \(y\) and the map factors through \((x,P,z)\) with entries polynomial in \(P,z\) over \(\mathbb C(x)\), the Jacobian in those coordinates is in \(\mathbb C(x)[P,z]\). Its pullback multiplied by \(P_y\) equals the nonzero constant determinant in \(x,y,z\). Since the only units of \(\mathbb C(x)[y,z]\) are nonzero elements of \(\mathbb C(x)\), the pulled-back determinant is independent of \(y,z\). Substitution of a nonconstant polynomial \(P(x,y)\) preserves positive \(P\)-degree as positive \(y\)-degree, so the determinant in \(x,P,z\) is independent of \(P,z\). The same argument applies with \(Q\) in place of \(P\).

### Constant nonzero \(r\)

Here \(T=d(x)=C_3-rq/2\) has degree five. The \(z^2\) coefficient of the determinant is \(rQ_y+2d'P_y\), and its vanishing gives \(Q=-2d'P/r+J(x)\). In this form the \(z\) coefficient is
\[
 2d'S_y+(4x(d')^2/r+2d''P-rJ')P_y.
\]
Thus
\[
 S=-\frac{d''}{2d'}P^2+
            \left(\frac{rJ'}{2d'}-\frac{2x d'}{r}\right)P+s(x).
\]
If \(P\) were independent of \(y\), every component would be independent of \(y\). Otherwise the map factors through \(x,P,z\) as above. At \(z=0\) its determinant in those coordinates equals
\[
                       r(S_xQ_P-S_PQ_x)-d'P Q_P .
\]
Its coefficient of \(P^2\) is
\[
                  \frac{d'd'''-3(d'')^2}{d'} .
\]
For instance, if \(Q=\gamma P+J\) and \(S=\eta P^2+\text{lower terms}\), this coefficient is \(r(\gamma\eta'-2\eta\gamma')\) with \(\gamma=-2d'/r,\ \eta=-d''/(2d')\). The leading term of \(d'\) has degree four, so the displayed scalar is nonzero. This contradicts the required independence of \(P\).

### Nonconstant \(r\)

Define
\[
 \lambda=x-\frac{r}{2r'},\quad K=r'+\frac{r r''}{r'},\quad
 a=r'''-\frac{3(r'')^2}{r'} .
\]
Work in \(\mathbb C(x)[y,z]\). Since \(T=r'P+d\), the \(z^2\) coefficient of the determinant equals
\[
                  (r-2xr')Q_y+2(r''P+d')P_y+2r'S_y .
\]
Its vanishing gives
\[
 S=\lambda Q-\frac{r''}{2r'}P^2-\frac{d'}{r'}P+s(x).
 \tag{1}
\]
Substituting (1), the coefficient of \(z\) in the determinant is the \(y\) derivative of
\[
     \frac{a}{3}P^3+
       \left(d''-\frac{3d' r''}{r'}\right)P^2
       -\left(\frac{2(d')^2}{r'}+2r's'\right)P
       -\left(KP+\frac{r d'}{r'}\right) Q .
 \tag{2}
\]
These coefficients are obtained by expanding the three derivative columns of \(F\) and collecting the designated powers of \(z\).

If \(P\) is independent of \(y\), (1) gives \(S=\lambda Q+s_0(x)\) and \(T\) depends only on \(x\). Then \(Q\) must depend on \(y\), or the determinant vanishes. In coordinates \(x,Q,z\), the constant coefficient in \(z\) of the determinant has \(Q\)-coefficient \(r\lambda'=r((r')^2+rr'')/(2(r')^2)\). This is nonzero since \((r')^2+r r''=(rr')'\) is nonzero for nonconstant \(r\). It contradicts the required independence of \(Q\).

Suppose \(P\) depends on \(y\). The polynomial \(K\) in the rational-function field is nonzero since \(r'K=(rr')'\). The vanishing of the derivative of (2) gives \((KP+r d'/r')Q\) as a polynomial in \(P\) of degree at most three. Dividing that polynomial by \(KP+r d'/r'\) in \(\mathbb C(x)[P]\) leaves a remainder in \(\mathbb C(x)\). The remainder must vanish because \(KP+r d'/r'\) is nonconstant in \(y\) and divides the remainder in \(\mathbb C(x)[y]\). Thus \(Q,S\) are polynomial in \(P\) over \(\mathbb C(x)\), both of degree at most two. Their quadratic coefficients are
\[
             A_*=\frac{a}{3K},\qquad B_*=\lambda A_*-\frac{r''}{2r'} .
\]
At \(z=0\) the determinant in \(x,P,z\) is
\[
 r(S_xQ_P-S_PQ_x)+P\bigl(r'Q_x-(r''P+d')Q_P\bigr).
 \tag{3}
\]
The coefficient of \(P^3\) is \(2r(A_*B_*'-B_*A_*')+r'A_*'-2r''A_*\).
If \(r=R x^j+O(x^{j-1})\) with \(j\ge2\), then
\[
 A_*=-\frac{j(j-1)}{3}x^{-2}+O(x^{-3}),\qquad
 B_*=-\frac{j^2-1}{3}x^{-1}+O(x^{-2}).
\]
The leading term of the \(P^3\) coefficient is therefore
\[
               \frac{2R j(j-1)(4j^2-1)}{9}x^{j-4},
\]
which is nonzero. The determinant cannot be independent of \(P\).

Finally, when \(r\) has degree one, translate \(x\) and change the first component by a constant multiple of the second to preserve \(A=(x,1,0)\). A constant rescaling of the third component gives \(r=x\). Then \(d=T-P=C_3-p+xq/2\) has degree five, since \(\deg p,\deg q\le3\). Formula (2) and the preceding division show that \(Q=d''P+J(x)\); by (1), \(S\) is also affine in \(P\). Formula (3) now has coefficient \(d'''\) at \(P^2\), from \(P Q_x\), since \(r'=1,\ r''=0\). This is nonzero, again contradicting independence of \(P\). This completes the rank-two exclusion.

## Constant rank-one coefficient

If \(A\) is a nonzero constant vector, equation (2) of the
[binary reduction](QUINTIC_VERONESE_BINARY_Z_REDUCTION.md#the-high-determinant-coefficients)
gives \(xQ_2-yP_2=0\) in the original coordinates. Therefore the leading cubic part of \(B\) is a linear form times \(N\). Apply a linear target change to put \(A=(0,0,h)\) for \(h\in\mathbb C^\times\). Write \(B=(p,q,r)\). The \(z^3\) coefficient of the determinant is \(2h\{p,q\}\), where \(\{\cdot,\cdot\}\) is the binary Jacobian. If the leading cubic part of \(B\) were a nonzero linear form times \(N\) before the target change, then \(p,q\) would have linearly independent leading cubics. Their binary Jacobian would have nonzero top part, by the homogeneous Euler identity. Thus \(B\) has degree at most two.

Now depress \(F_3\) by \(z\mapsto z-r/(2h)\). The map has the form \((P+pz,Q+qz,T+h z^2)\), and \(P,Q,T\) have the same degree-five leading parts as their respective coordinates before depressing, because all modifications have degree at most four. The coefficient of \(z\) in the determinant is
\[
              2h\{P,Q\}+p\{q,T\}-q\{p,T\}.
\]
The last two terms have degree at most seven, while the first has nonzero degree-eight part from the independent quintic forms. This contradiction excludes constant \(h\).

## Leading linear coefficients in the remaining rank-one case

Write the first two components as \(C_1+pz,C_2+qz\) in the
frame \(A=(0,0,x)\). Then \(p,q\) are linearly independent
over \(\mathbb C\) and their degree-three homogeneous parts
are dependent. In fact, if \(p,q\) are dependent, choose \(q=0\)
by a linear change of the first two components. At least one
linear coefficient is nonzero, since otherwise \(F_z\) vanishes
on \(x=0\), so take \(p\ne0\). The \(z^2\) determinant coefficient
then gives \(2x\{p,C_2\}-p(C_2)_y=0\), or
\(\{p^2/x,C_2\}=0\). As in the binary reduction, \(C_2\) has
degree five and no critical point, so \(\mathbb C(C_2)\) is
relatively algebraically closed and contains \(p^2/x\).

If \(x\nmid p\), the reduced denominator of \(p^2/x\) has degree
one. A nonconstant reduced denominator in \(C_2\) would have
degree at least five, and a constant denominator would make
the ratio polynomial. If \(x\mid p\), then \(p^2/x=xt^2\)
for a polynomial \(t=p/x\). This is a nonconstant polynomial
of degree at most five. A polynomial in \(\mathbb C(C_2)\)
is a polynomial in \(C_2\), by the Bézout argument for reduced
numerators and denominators. Thus \(xt^2\) would have to be affine
in \(C_2\), in particular of degree five. This makes \(t\)
nonconstant, but then \(xt^2\) has zero gradient along \(t=0\),
which is impossible for an affine transform of \(C_2\).
This proves independence. Finally the \(z^3\) equation is
\(2x\{p,q\}-p q_y+q p_y=0\). At binary degree five its left side
is \((5/3)x\{p_3,q_3\}\) by Euler's identity, giving the dependence
of the cubic parts.

The [symbolic checker](verify_quintic_veronese_binary_z.py) replays the
mixed determinant identity of the binary reduction and the rank-two
coefficient equations.

**Theorem QuinticVeroneseBinaryZExclusion.** In characteristic zero,
a degree-five three-variable polynomial Keller map cannot have
leading part \(g(x,y)N(x,y)\) with a binary homogeneous cubic
\(g\ne0\) and \(N(x,y)=(x^2,xy,y^2)\).

Work in the last frame \(A=(0,0,x)\) using the previous exclusions
and the embedding of the coefficient field into \(\mathbb C\).
Here \(F_1=C_1+pz\), \(F_2=C_2+qz\), \(F_3=C_3+rz+xz^2\),
with \(p,q\) independent over \(\mathbb C\) of degree at most
three. We first show that \(p,q\) depend on \(x\) alone.

Put \(x=s^2\) temporarily and take \(a=p(s^2,y)/s\),
\(b=q(s^2,y)/s\) in \(R=\mathbb C[s,s^{-1},y]\).
Their binary Jacobian in \(s,y\) is
\((2x\{p,q\}-p q_y+q p_y)/s^3\), zero by the \(z^3\) equation.
By [Gordan's theorem](sources/petravchuk-iena-math-0701588v2.pdf)
as used in Lemma 1, \(\mathbb C(a,b)=\mathbb C(\ell)\).
Suppose \(\ell\) depends on \(y\), and write it as \(U/V\)
with \(U,V\) coprime in \(R\). If a rational function expressing
\(a\) or \(b\) in terms of \(\ell\) has a pole at a finite
\(c\), then \(U-cV\) is a unit in \(R\). Otherwise a prime
factor of \(U-cV\) gives a pole of \(a\) or \(b\), since it
does not divide \(V\). For a pole at infinity, \(V\) must be
a unit by the same argument with \(1/\ell\). Two distinct
poles would make \(U,V\) linear combinations of units and hence
independent of \(y\), since the units of \(R\) are
\(\mathbb C^\times s^{\mathbb Z}\). There is at least one pole,
so there is exactly one. Change the rational parameter to put
that pole at infinity. Then \(\ell\in R\), and \(a,b\) are
polynomials in \(\ell\).

The involution \(s\mapsto -s\) sends \(a,b\) to their negatives,
so on their function field it fixes the unique pole. Its action
on \(\ell\) is affine of order two and nontrivial. After a
translation of \(\ell\) it sends \(\ell\) to \(-\ell\).
Thus \(a,b\) are odd polynomials in \(\ell\). By independence,
the maximum \(m\) of their degrees is at least three. If
\(\ell\) has a pole at \(s=0\) over \(\mathbb C(y)\), one of
\(a,b\) has pole order at least three there, contrary to
their defining expressions. Hence \(\ell=sL(s^2,y)\) for
a polynomial \(L\) which depends on \(y\). The corresponding
original coefficient \(p\) or \(q\) then has degree
\((m+1)/2+m\deg L\ge5\), since the lower powers of \(\ell\)
have strictly lower degree after multiplication by \(s\).
This is impossible. Therefore \(\ell\) is independent of \(y\),
so \(p,q\in\mathbb C[x]\).

Now work over \(\mathbb C(x)\) in coordinates \(u=F_1,y\) and
write
\[
 F_2=\lambda u+W,\qquad F_3=a_2 u^2+R_1u+R_0,\qquad
 \lambda=q/p,\quad W=C_2-\lambda C_1,\quad a_2=x/p^2.
\]
The coefficient functions \(W,R_1,R_0\) are polynomials in
\(y\) over \(\mathbb C(x)\). At fixed \(u\), the binary determinant
of \(F_2,F_3\) in \(x,y\) is \(\kappa/p\), for the nonzero scalar
\(\kappa=\det JF\). Since \(\lambda\) is nonconstant,
the \(u^2\) and \(u\) coefficients successively give
\[
 R_1=\mu W+b,\qquad
 R_0=\frac{\mu'}{2\lambda'} W^2+\frac{b'}{\lambda'} W+c,\qquad
 \mu=\frac{a_2'}{\lambda'},
\]
where \(b,c\in\mathbb C(x)\) and primes denote \(x\)-derivatives.
In particular, \(W\) depends on \(y\), or the determinant
vanishes. Its constant coefficient in \(u\) is
\[
    -\left[
      \left(\frac{\mu'}{2\lambda'}\right)' W^2
       +\left(\frac{b'}{\lambda'}\right)' W+c'
      \right]W_y=\kappa/p .
\]
Both factors are units in \(\mathbb C(x)[y]\). The bracketed
factor is independent of \(y\), so \(\mu'/(2\lambda')=A\) and
\(b'/\lambda'=B\) are constants. Integrating in \(\mathbb C(x)\)
gives \(\mu=2A\lambda+C\), \(b=B\lambda+D\) and
\(a_2=A\lambda^2+C\lambda+E\) with constants \(C,D,E\).
Substitution gives
\[
                 F_3=A F_2^2+C F_1 F_2+E F_1^2+B F_2+D F_1+c(x).
\]
Here \(c(x)\) is polynomial because all the other terms are
polynomial. Subtracting the displayed polynomial expression
in the first two coordinates gives a Keller map whose third
coordinate is \(c(x)\), so \(c'\) is a polynomial unit and
\(\deg c=1\). The homogeneous quadratic part of the expression
in \(F_1,F_2\) is nonzero, since \(F_3\) has \(z^2\) coefficient
\(x\). At \(z=0\), this quadratic evaluated on \(C_1,C_2\)
has degree ten: the degree-five leading forms are independent,
and every linear factor of the binary quadratic evaluates to
a nonzero quintic form. This contradicts \(\deg C_3=5\).

**Corollary QuinticVeroneseLeadingPartExclusion.** In characteristic
zero, a three-variable degree-five polynomial Keller map cannot
have leading part \(g(x,y,z)(x^2,xy,y^2)\) for any nonzero
homogeneous cubic \(g\).

The theorem covers \(\deg_z g=0\).
The [linear](QUINTIC_VERONESE_SQUARE_EXCLUSION.md),
[quadratic](QUINTIC_VERONESE_QUADRATIC_Z_EXCLUSION.md) and
[cubic](QUINTIC_VERONESE_LEADING_RESIDUAL.md) cases
cover all positive \(z\)-degrees.

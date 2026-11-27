# Weight reductions at the square cusp

**Theorem SexticKellerSquareCuspWeight456Reduction.**
Let \(P=xz+y^2\). If a characteristic-zero Keller map has
ordinary highest tuple \((P^3,Px^4,x^6)\), its highest tuple
in source weights \((4,5,6)\) is also \((P^3,Px^4,x^6)\).
In particular the coefficient \(c\) of \(xP^2\) in the ordinary
degree-five correction of the third component vanishes.

**Theorem SexticKellerSquareCuspBinaryEdgeReduction.**
Under the same hypothesis the highest tuple in weights
\((10,15,18)\) has the form
\[
                   (y^6,\ x^4 y^2,\ x^6+k x^3y^2+e y^4).
\]
In particular the degree-five coefficient \(b\) in
\(E_2=b xP^2+x^2 M_3\) also vanishes.

## The initial weight and the transverse terms

The
[degree-five square-cusp constraints](SEXTIC_KELLER_TANGENT_CONIC_CUBIC_OUTER.md#degree-five-divisibility-when-the-cusp-is-at-the-square)
give \(E_2=b xP^2+x^2 M_3\) and \(E_3=c xP^2+x^3 M_2\).
At weights \((4,5,6)\) the highest tuple must therefore be
\[
           (P^3+u z^5,\ Px^4+d x^2z^3,\
                           x^6+c xP^2+v x^3z^2+w z^4)
\]
of weights \(30,26,24\). These are the only tying
terms by the ordinary degree bounds and the displayed
factors. Its determinant is zero, because its
expected determinant weight is \(65>0\).
Compare the following coefficients in the order shown.
\[
\begin{array}{c|r}
 \text{monomial}&\text{coefficient}\\ \hline
 x^6y^7z      & -48v\\
 x^7y^5z^2    & 108d-144v\\
 x^9yz^4      & 108d-60u-48v\\
 x^3y^7z^3    & 72cd-96w
\end{array}
\]
This gives \(u=d=v=w=0\) without a condition on \(c\).

If \(c\ne0\), use \(P,x\) as normalization parameters
of weights ten and four. They are coprime. Their
rational field is relatively closed since at fixed
\(P,x\), one solves rationally for \(z\), leaving \(y\)
free. The pair is onto, including at \(x=0\) where
\(P=y^2\); hence polynomial contraction gives the
polynomial kernel \(k[P,x]\).
The derivative minors of \((P^3,Px^4,x^6+c xP^2)\)
in these parameters, up to sign, are
\[
           12P^3x^3,\quad
           3P^2(6x^5+cP^2),\quad
           x^4(6x^5-7cP^2).
\]
Their gcd is one, because the second excludes \(x\)
and the third excludes \(P\). The
[weighted saturation argument](SEXTIC_KELLER_BINARY_REDUCED_LOWER_PENCILS.md#lifting-the-finite-normalization-at-the-first-new-weight)
therefore excludes a completion. Its determinant order
\(65\) exceeds the largest component weight \(30\).
This proves \(c=0\).

Next keep the weights of \(y,z\) equal to five and
six, and reduce the weight \(s\) of \(x\) from four
to \(10/3\). The first component has highest part
\(y^6\) throughout \(s<4\) in this interval: all its
terms below ordinary degree six had weight strictly
less than \(30\) at \(s=4\), and decreasing \(s\)
cannot increase them. The reference highest parts of
the other two components are \(x^4 y^2,x^6\).
After the preceding eliminations, the possible ties
strictly above \(10/3\) are as follows.
\[
\begin{array}{c|l|l}
 s &\text{second component}&\text{third component}\\ \hline
 23/6 & & yz^3\\
 11/3 & & x^3yz,\ y^2z^2\\
 18/5 & & xz^3\\
 7/2 & x^2yz^2,\ z^4 & y^3z\\
 17/5 & & xyz^2
\end{array}
\]
The table follows by equating \(s i+5j+6l\) to
\(4s+10\) or \(6s\) for a monomial \(x^iy^jz^l\),
using the supports of the degree-five corrections
and all monomials of ordinary degree at most four.
At a tie in just the third component the first two
leading forms fix \(x,y\) up to finitely many choices.
Thus a \(z\)-dependent tie would give a nonzero
leading determinant. At \(s=7/2\), the only possible
leading pair on fixing generic \(y\) is
\[
 B=x^4 y^2+d_1x^2yz^2+d_2z^4,\qquad
 C=x^6+d_3y^3z .
\]
If \(d_3\ne0\), the generic coordinates \(x,C\) on
that plane express \(z\) as \((C-x^6)/(d_3y^3)\).
The zero leading determinant requires \(B\) to be
independent of \(x\) in these coordinates. This is
impossible: its largest \(x\)-degree is \(24\) if
\(d_2\ne0\), \(14\) if only \(d_1\ne0\), and \(4\)
otherwise. Thus \(d_3=0\) and the determinant then
also forces \(d_1=d_2=0\).
At \(s=10/3\) the remaining leading tuple, after
clearing weights, is therefore
\[
                  (y^6,\ x^4y^2+bxy^4,\ x^6+k x^3y^2+e y^4).
 \tag{1}
\]

## Exclusion of the mixed binary term

Suppose \(b\ne0\). Write \(Y=y^2\).
The image field of (1) is \(k(x,Y)\). For the
projective calculation put \(h=x^3/y^2\),
\(\rho=y/x\), and \(f(h)=h^2+k h+e\).
Thus \(x=h\rho^2,\ y=h\rho^3\). Denote the three
entries of (1) by \(A_0,B_0,C_0\). Their ratios
\(C_0^3/A_0^2,\ B_0^3/(A_0 C_0^2)\) are
\[
                    f(h)^3,\qquad
                    R(h)=h(h+b)^3/f(h)^2 .
\]
They generate \(k(h)\). Indeed the possible proper
index divides six and the degree of \(R\), which
is between one and four. Index two would give an
involution on \(h\) preserving \(f^3,R\). It fixes
infinity, the sole pole of \(f^3\), hence is affine
with slope \(-1\). Its invariance of \(f^3\) then
implies invariance of \(f\), by the leading coefficient.
But it would also preserve \(h(h+b)^3\), and could
not fix both distinct finite roots of different
multiplicities.

Index three requires \(R\) to have degree three.
Its numerator has zeros \(0,-b\) of orders one
and three, so this requires cancellation of the
factor \(h\) alone. Thus \(f=h(h+k)\) and \(k\ne b\),
and \(f^3\in k(R)\). Since infinity is its unique
pole and \(R(\infty)=1\), the fibre \(R=1\) can
have no finite point. But the numerator of \(R-1\)
is
\[
              (3b-2k)h^2+(3b^2-k^2)h+b^3 .
\]
It cannot be constant since \(b\ne0\) and the
first two coefficients cannot both vanish. Index
three is also impossible. This proves \(h\) is in
the image field. The entries in (1) are an element
of \(k(h)^\times\) times \(\rho^{18},\rho^{14},\rho^{12}\),
respectively. Their exponents have gcd two, giving
the asserted field \(k(h,\rho^2)=k(x,Y)\).

The two elements \(x,Y\) are integral over the
leading image ring: \(Y^3=A_0\), and \(x\) satisfies
\(x^6+kx^3Y+eY^2-C_0=0\).
Increase the weight of \(z\) from eighteen to
the first tie. By the
[binary normalization lifting argument](SEXTIC_KELLER_BINARY_REDUCED_LOWER_PENCILS.md#lifting-the-finite-normalization-at-the-first-new-weight),
the field generators lift to polynomials of
weights ten and thirty restricting to \(x,Y\) at
\(z=0\). Since the weight of the new variable is
above eighteen, only \(Y\) can lift nontrivially.
The monomial \(x^3z^3\) in the first component ties at
weight twenty, so the first tie is at most twenty.
A monomial of weight thirty involving \(z\) between
weights eighteen and twenty can only be \(xz\) at
weight twenty. Hence the lift is \(H=y^2+\delta xz\),
with \(\delta
e0\). Its expected top determinant weight
is \(175>0\), so the normalization lift applies.
The leading tuple is
binary in the parameters \(x,Z\) and these again
give its finite normalization. Increase the weight
of the unused variable \(y\) from fifteen. The only
possible lifts involving it are a term proportional
to \(xy\) at weight twenty or to \(y\) at weight
thirty in \(Z\). In each case take the lifted
parameter as the new coordinate and continue
increasing the weight of the unused \(y\).
Above thirty neither generator can lift nontrivially,
contrary to existence of a first tie for a Keller
map. Thus the weight-thirty alternative is impossible.

In the weight-twenty case put \(H=y^2+\delta xz\).
The leading tuple is \((H^3,f,g)\), where
\[
              f=x^4H+bxH^2,\qquad
              g=x^6+kx^3H+eH^2 .
\]
The pair \((x,H)\) is onto and has relatively closed
field, by the same argument as for \((x,P)\).
In the parameter ring \(k[x,H]\), the derivative
minors are coprime unless \(k=b/2\).
Indeed \(f_x=H(4x^3+bH)\) and
\(g_x=3x^2(2x^3+kH)\), and the minor of \(f,g\)
is nonzero generically at \(H=0\).
The weighted saturation theorem excludes the coprime
case.

### Descent on the remaining critical divisor

Suppose \(k=b/2\) and put \(D=4x^3+bH\).
The minors in the source ring have common factor
only \(D\), which is irreducible and reduced there
since it is linear in \(z\) with its constant term
coprime to \(x\). After division by \(D\), the two
minors involving \(H^3\) are scalar multiples of
\(H^3,H^2 x^2\), and the remaining one at \(H=0\)
is proportional to \(x^6\).

Use reference triples of the form \(G(\widetilde
x,\widetilde H)\) with the same highest forms.
At a first weighted error loss \(r\) write \(E\)
for the error vector of that weight. The linear
determinant equation puts its normal contraction
in \(S=k[x,H]\). The source is free over \(S\)
with basis \(z^i,yz^i\) for \(i\ge0\), since
\(y^2=H-\delta xz\). Absorb the coefficient of
the basis element \(1\) into the reference, leaving
a tangent error with no scalar coefficient in
this basis. If \(r\ge70\), the last two entries
vanish after this subtraction, so the tangent
equation forces the first to vanish as well.

For \(r<70\), use the generic point of \(D=0\), fix
the actual first component to its generic value, and
take the centered quadratic parameter \(s\) in the
reference second component. This chart is regular
because \(x,H\) have independent differentials there.
At the leading level \(H=H_0\), the critical point
of \(f\) is a simple critical point of the scalar
function, since \(f_{xx}=12x^2H_0\) is nonzero there.
Also \(g_x/f_x=3x^2/(2H_0)\) has nonzero derivative
there. After subtracting the quadratic coefficient of
the third component in the \(s\)-chart, its leading
cubic coefficient is consequently nonzero.

The
[ordinary-cusp pole comparison](SEXTIC_KELLER_BINARY_SIMPLE_CUSP_DIVISORS.md#the-pole-comparison-eliminates-the-local-denominator)
now implies that the restriction of the first error
of the second component on this level, corrected for
the change of the first scalar, is algebraic over the
level constants. The same argument applies with the
free coordinate \(y\): at fixed first two actual
outputs the coefficients of the third through order
less than the determinant order are independent of
\(y\). At an earlier nonzero linear coefficient of
the reference third component, the order-\((p+r)\)
simple pole gives a scalar linear equation for the
restricted error; otherwise its order-\(2r\) simple
pole gives a nonzero quadratic equation because the
cubic coefficient is nonzero. These orders are
available because \(2r<140<175\).

At the leading fixed values \(x=x_0,H=H_0\) the
coordinate ring of the generic critical fibre still
has the free basis \(z^i,yz^i\) over the enlarged
constant field, by the same monic quadratic relation.
The restricted error is a polynomial combination of
these basis elements with no coefficient of \(1\).
Being algebraic over the constants in a rational
function field in \(y\), it is scalar and thus zero.
For the tangent coefficients, the \(H\)-coefficient
is already regular along \(D=0\) by the first output.
The vanishing just obtained makes the remaining
numerator divisible by \(D\), so the \(x\)-coefficient
is also regular there. Elsewhere the derivative
minors solve for both coefficients without a pole.
They are therefore polynomial, and can be absorbed
into \(\widetilde x,\widetilde H\) at that loss.
This descends through all positive weights and
constants, giving a map depending on only two
polynomials and a zero determinant. The contradiction
proves \(b=0\) in (1).

## Returning from the binary edge to the conic

**Theorem SexticKellerSquareCuspWeight101520Reduction.**
In the square-cusp setting above, the highest tuple in
weights \((10,15,20)\) is
\[
             (P^3,\ Px^4,\ x^6+kx^3P+eP^2),
\]
with the same \(k,e\) as on the binary edge.

Increase only the weight of \(z\) from eighteen to twenty.
A monomial involving \(z\) ties if \(10i+15j+wl\) equals
the relevant component weight \(90,70,60\). In the open
interval the possible ties are \(yz^3\) in the second
entry at \(55/3\) and \(yz^4\) in the first entry at
\(75/4\). This uses the divisibilities of both
degree-five corrections and the eliminated coefficients
in the preceding two theorems. At the first tie the other
two entries are independent in \(x,y\). At the second
tie the two unchanged entries \(x^4y^2,\ x^6+kx^3y^2+ey^4\)
are also independent: their binary determinant contains
a nonzero contribution from \(x^6\) of strictly higher
\(x\)-degree than any other contribution.
Thus each tie would give a nonzero weighted determinant
unless its coefficient vanishes.

At weight twenty use the polynomial variables \(x,P,z\)
in the subring \(k[x,y^2,z]\). The leading tuple has the
following general form, by the same exponent test:
\[
\begin{aligned}
 A&=P^3+a Pz^3+a' xz^4,\\
 B&=x^4P+d x^2Pz+d' x^3z^2+f Pz^2+f' xz^3,\\
 C&=x^6+kx^3P+eP^2+g x^4z+g' xPz+g'' x^2z^2+g''' z^3 .
\end{aligned}
\]
The actual source Jacobian of this tuple is \(2y\)
times its determinant in \(x,P,z\). The expected weight
is \(175>0\), so the latter determinant is zero.

In order, the coefficients of \(x^8P^2 z\),
\(x^9Pz^2\), \(x^6P^2z^2\), and \(x^{10}z^3\) give
\(d'=a=f'=a'=0\). Their expressions are
\(36d'\), \(-18a\), \(3ak+15d'g+54f'\), and \(-24a'\),
respectively. After these vanishings compare the following:
\[
\begin{array}{c|r}
 x^7P^3 &18d-12g\\
 x^4P^4 &9dk-12g'\\
 x^5P^3z &6dg+36f-24g''\\
 x^3P^3z^2 &-6dg''+24fg-36g'''\\
 xP^3z^3 &-18dg'''+12fg''\\
 x^2P^4z &-3dg'+18fk .
\end{array}
\]
The first four determine \(g,g',g'',g'''\). Substitution
in the fifth gives \(18(f-d^2/4)^2=0\), and the last
then gives \(kd^2=0\). With \(\gamma=d/2\) and
\(Q=x^2+\gamma z\), the tuple is
\((P^3,PQ^2,Q^3+k x^3P+eP^2)\) and \(k\gamma=0\).

If \(\gamma\ne0\), the
[weighted leading-pair slice](KELLER_LEADING_PAIR_RATIONAL_SLICE.md)
applied to \(P^3,PQ^2\) gives a rational primitive for
the relative area differential on a geometric generic
fibre of \((P,Q)\). Indeed
\(d(P^3)\wedge d(PQ^2)=6P^3Q\,dP\wedge dQ\),
with nonzero constant factor on the generic fibre.
For generic \(P=S,Q=T\), eliminating \(z\) gives
\[
                      y^2=S+(x^3-Tx)/\gamma .
\]
The right side is a squarefree cubic over the geometric
generic field. Also \(dP\wedge dQ\wedge dx
=2\gamma y\,dx\wedge dy\wedge dz\), so the relative
differential is a nonzero scalar multiple of \(dx/y\).
This differential is regular on the smooth projective
model: at a simple finite root, \(y\) is a parameter,
and at infinity \(x,y\) have pole orders two and three.
A rational function with a pole has differential with
a pole in characteristic zero. Thus this nonzero
regular differential cannot have a rational primitive.
The contradiction gives \(\gamma=0\), as claimed.


## Axial pole exclusion with two nonzero coefficients

Write \(A=F_1,\ B=F_2,\ C=F_3\). Use the preceding
\((2,3,4)\) reduction, with third form
\(x^6+kx^3P+eP^2\), and assume \(ke\ne0\).
Write \(C^\flat=C-eA^{2/3}\) when taking a filtered
initial form, using the binomial branch with highest
term \((xz)^2\) for \(A^{2/3}\).
The polynomial
\[
 D=(C^3-e^3 A^2)/(3e^2)
\]
records the same leading information without any
fractional power. If \(A_*=T^3\), \(C_*=eT^2\) and
\(C^\flat\) has strictly lower weight than \(C\), then
\(D_*=T^4(C^\flat)_*\).

**Lemma SexticKellerSquareCuspAxialPoleExclusion.**
Suppose that at weights \((2,3,w)\) with \(w\ge5\) the
highest forms of \(A,C\) are \(T^3,eT^2\), where
\(T=xz\). The second component may be replaced by a
constant multiple of \(C\) subtracted from it. Suppose
its highest form and the normalized third form are
\[
 B_*=Tx^4+H(T)+\frac{bT^m}{x^u},\qquad
 (C^\flat)_*=kTx^3+\frac{cT^n}{x^v},
\]
of weights \(w+10,w+8\), respectively. Here \(H\) is a
polynomial, the displayed powers are integral, \(m,n\ge0\), \(1\le
u<4\) if \(b\ne0\), and \(v>0\) if \(c\ne0\).
Then \(b=c=0\).

Indeed \(T\) is in the field of the four polynomial
highest forms of \(A,B,C,D\) by the ratio of \(T^3\)
and \(eT^2\). Over \(k(T)\), regard the two displayed
functions as rational functions of \(x\). At infinity
they have pole orders four and three. In the intermediate
function field they generate, the place under infinity
has ramification index one by these two orders. Its
residue field is \(k(T)\). Any other place above it
would have a pole of the first function of order at
least four, whereas its only other possible pole is
at zero of order \(u<4\). The degree formula at this
place shows that the intermediate field is all of
\(k(T)(x)\). In particular the four polynomial highest
forms have fraction field \(k(x,z)\).

Increase the weight \(s\) of \(y\) from three to the
first tie in these four polynomials. There must be a
tie by \(s\le3(w+2)\), since at least one of \(A,B,C\)
depends on \(y\) and \(3(w+2)\) is their largest current
weight. At a tie, all three-component leading brackets
vanish. Besides the Keller bracket and the identically
zero bracket for \(A,C,D\), the brackets involving
\(A,B,D\) and \(B,C,D\) are constant multiples of
\(C^2\) and \(A\), respectively. The expected weights
exceed those actual weights by \(4w+22-s>0\).
The algebra of the new highest forms therefore has
transcendence degree at most two. Specialization at
\(y=0\) maps it onto the old algebra of transcendence
degree two, hence gives an isomorphism. This also
identifies their graded fraction fields.

Denote the rational lifts of \(T,x,z\) by \(T',X,Z\).
The equality \((T')^3=A_*\) makes \(T'\) polynomial
by normality of the polynomial source ring. Its
weight is \(w+2\), and it specializes to \(xz\).
Since \(s>3\), no new monomial containing \(y\)
can also contain \(z\). Thus \(T'=xz+f(x,y)\) is
squarefree. It is irreducible if \(x\nmid f\), and
otherwise is \(x\) times a linear polynomial in \(z\)
with coefficient one. At a prime divisor where \(X\)
has a pole, \(T'X^4\) would have negative order and
would be the unique polar term of the lifted \(B_*\).
This uses the bound \(0\le\operatorname{ord}T'\le1\).
So \(X\) is polynomial, and its weight two and
specialization force \(X=x\).

If \(Z=T'/x\) had a pole, it would be at \(x=0\),
where \(T'\) is then a unit. For \(b\ne0\) the
lifted term \(b(T')^m/x^u\) in \(B_*\) makes a
pole there. For \(c\ne0\) the term \(c(T')^{n+4}/x^v\)
in \(D_*\) does the same. Either contradicts
polynomiality. Thus if \(b\) or \(c\) is nonzero,
\(Z\) is polynomial. At a genuine first tie
it must have a nonzero correction \(Z=z+g(x,y)\):
otherwise all four highest forms would be the old
forms in \(x,z\). Use this polynomial shear as a
source change. It removes the tie and preserves the
binary forms at all preceding weights. A subsequent
tie would occur at a strictly larger \(s\). There
are only finitely many possible shear weights
\(2i+sj=w\) with \(i\ge0,j\ge1,s>3\).
Iterating therefore contradicts existence of a tie.
This proves the lemma.

## An affine scalar in the mixed square-cusp sector

**Theorem SexticKellerSquareCuspMixedConstantsAffineReduction.**
In the \((2,3,4)\) square-cusp tuple
\((P^3,Px^4,x^6+kx^3P+eP^2)\), suppose \(ke\ne0\).
After a constant target shear in the second component
and a translation in \(x\), the three complete axial
degrees are \(3,1,2\), with highest \(z\)-terms
\[
               x^3z^3,\qquad a(x)z,\qquad e x^2 z^2 ,
 \qquad \deg a=5 .
\]
Over the algebraic closure, if a completion in this sector
is noninvertible, then
\(a\) has a repeated root at zero. Only at that root
can the constant term \(g(x,y)\) of the second
component restrict to a constant in \(y\). In that
case \(g(0,y)\) is constant and \(g_x(0,y)\) is
a nonzero constant. There are at most two other
distinct roots of \(a\).

The \((2,3,4)\) support gives
\[
\begin{aligned}
 A &=\theta z^4+(x^3+bxy+c x^2+d y+f x+g_0)z^3
                  +A_2z^2+\cdots,\\
 B &=\beta z^3+(\alpha xy+\sigma x^2+\tau y+\lambda x+\mu)z^2
                  +B_1z+B_0,\\
 C &=(e x^2+\eta y+\xi x+\chi)z^2+C_1z+C_0 .
\end{aligned}
\]
Here \(A_2,B_1,C_1\) have \((2,3)\)-highest parts
\(3x^2y^2,x^5,2exy^2+kx^4\), respectively.

At \(w=5\) the \(z\)-leading expansion of the fractional
power gives
\[
 B_*=Tx^4+\alpha xyz^2+\beta z^3,\quad
 (C^\flat)_*=kTx^3+(\eta-2eb/3)yz^2- (2e\theta/3)z^3/x .
\]
To see the fractional terms, write \(A=P^3+(A-P^3)\)
and retain \(P^2+\tfrac23(A-P^3)/P\) at weights
at least thirteen. The only weight-twenty terms
of \(A-P^3\) are \(\theta z^4+bxyz^3\).
On fixing \(T\), the \(y\)-coefficients in the
displayed \(B_*,(C^\flat)_*\) are constant multiples
of \(x^{-1},x^{-2}\), respectively. The coefficient
of \(y\) in their planar bracket forces their product
to vanish. If either vanishes but the other does not,
the remaining derivative in \(x\) gives a contradiction,
since the \(y\)-free parts have nonzero positive powers
\(Tx^4,kTx^3\). Hence \(\alpha=0,\ \eta=2eb/3\).
The axial pole exclusion then gives \(\beta=\theta=0\).

For the higher weights the fractional expansion can be
written uniformly. Put \(a_3=x^3+bxy+cx^2+dy+fx+g_0\).
The coefficients of \(z^2,z\) in \(C^\flat\) are
\[
 R_2=e x^2+\eta y+\xi x+\chi-e a_3^{2/3},\qquad
 R_1=C_1-\frac{2e}{3}A_2 a_3^{-1/3}.
\]
Expand powers of \(a_3\) from its \((2,3)\)-highest
term \(x^3\). Then \(R_1\) has highest part \(kx^4\).
Every monomial of the remaining Laurent series in
\(z\), with \(z\)-exponent \(j\le0\), has
\((2,3)\)-weight at most \(12-4j\). This follows by
expanding \(A^{2/3}\) from \(a_3 z^3\): the coefficients
of the ratios of all other \(A\)-terms to \(a_3z^3\)
have total \((2,3,4)\)-weight at most zero. Thus for
\(w>4\) these terms cannot tie \(kx^4 z\).
The same expansion justifies the coefficient
comparisons in each weight below.

At \(w=6\), the weighted degree-two part of \(R_2\) is
\[
                (\xi-2ec/3)x+\frac{e b^2}{9}y^2/x^2 .
\]
Here \(B_*=Tx^4+\sigma T^2\). On fixing \(T,B_*\),
the free \(y\) forces \(b=0\); the axial pole exclusion
then forces \(\xi=2ec/3\). In particular \(\eta=0\).
Replace \(B\) by \(B-(\sigma/e)C\), cancelling its
\(x^2z^2\) term. At \(w=7\) the only possible new
term of this \(B_*\) is \(\tau y z^2\), and the
only correction to \(kTx^3\) is
\(-2ed\,y z^2/(3x)\). Their \(y\)-coefficients at
fixed \(T\) now have distinct powers \(x^{-2},x^{-3}\).
The same planar bracket argument gives \(\tau=d=0\).

Thus \(a_3\) is univariate and the remaining possible
\(z^2\) coefficient of the sheared second component
is linear in \(x\). At weights eight, ten and twelve
the comparisons have the following forms. In each row
all the earlier corrections have vanished.
\[
\begin{array}{c|cc}
 w & B_*-Tx^4 & (C^\flat)_*-kTx^3\\ \hline
 8 & \text{scalar}\cdot x z^2 & [x^0]R_2\ z^2\\
 10& \text{scalar}\cdot z^2   & [x^{-1}]R_2\ z^2/x\\
 12& 0 & [x^{-2}]R_2\ z^2/x^2 .
\end{array}
\]
With \(z=T/x\), these are exactly the polar terms
of the lemma. All coefficients in the table therefore
vanish. This gives an affine second component and
\[
          a_3(x)^{2/3} - c_2(x)/e=O(x^{-3})
\]
in Laurent series at infinity, where \(c_2=[z^2]C\).
Translate \(x\) to center the monic cubic \(a_3\) as
\(x^3+f_1x+g_1\). The estimate gives
\(c_2/e=x^2+2 f_1/3\), and \(a_3^2-(c_2/e)^3\)
has degree at most one. Its coefficients of \(x^3\)
and \(x^2\) are \(2g_1\) and \(-f_1^2/3\), so
\(f_1=g_1=0\).

It remains to show the second coefficient depends only
on \(x\). At ordinary axial \(z=\infty\), the polynomial
\(D\) now has highest \(z\)-part \(T^5 v\), where
\[
                 u=B_1/x,\qquad
                 v=C_1/x-2e A_2/(3x^2).
\]
The leading bracket of \(A,B,D\) shows that \(u,v\)
are algebraically dependent in \(k(x,y)\): fixing
\(T=xz\), their factors \(T,T^5\) are constant,
and the bracket has expected \(z\)-degree eight versus
actual degree four. Their \((2,3)\)-highest parts are
\(x^4,kx^3\). Over \(k(y)\), the pole orders at
\(x=\infty\) are four and three. At \(x=0\) the pole
order of \(u\) is at most one, and there are no other
possible poles. The same place-degree argument as in
the lemma gives \(k(y)(u,v)=k(y)(x)\).

Over the algebraic closure of \(k\), the curve with
field \(k(u,v)\) is rational: specialize \(y\) to a
scalar, retaining the nonconstant \(u\), to embed its
field in \(k(x)\) and apply Lüroth's theorem. Choose
a parameter with its simple pole at the place under
\(x=\infty\). The degree-one parametrization over
\(k(y)\) then writes it as \(t=a(y)x+b(y)\), and the
coefficient of \(x^4\) in \(u\) makes \(a(y)\) a
nonzero scalar, which can be normalized to one.
If either function has a pole at \(x=0\), then
\(b(y)\) equals one of the fixed finite poles of
that function in \(t\), so \(b(y)\) is constant.
If both are regular there, \(u\) is polynomial,
\(u=u_0(x+b(y))\) with \(u_0\) monic of degree four.
Its \(x^3\) coefficient is \(4b(y)\) plus a
constant. This coefficient is polynomial in \(y\),
and the \((2,3)\) bound on \(u\) makes it constant.
Thus \(B_1=xu\) is a monic univariate polynomial of
degree five, as asserted.

### Restrictions at the roots of the affine coefficient

Put \(R=B=a(x)z+g(x,y)\) in the last frame.
The [normal-residue proposition](KELLER_LINEAR_SOURCE_COEFFICIENT_CUBIC_CLASSIFICATION.md#normal-residues-along-a-repeated-line)
applies without a degree restriction. It says that
\(g(\zeta,y)\) for a root \(\zeta\) of \(a\) is
either affine with nonzero slope, or constant.
The constant case requires a repeated root and
\(g_x(\zeta,y)\) a nonzero constant.

Here the \((2,3)\)-bound of fourteen on \(g\) makes
its \(y\)-degree at most four and its \(y^4\)
coefficient scalar. That coefficient vanishes
by restriction to any root, so the \(y\)-degree
is at most three. If every restriction is affine
with nonzero slope, the
[parallel-fibre coordinate construction](KELLER_LINEAR_SOURCE_COEFFICIENT_CUBIC_CLASSIFICATION.md#coordinates-on-the-parallel-support-fibres)
gives on \(R=t\) the inverse coordinates
\(y=\psi(x,t)+a(x)M,\ z=(t-g(x,y))/a(x)\).
The first has degree at most six in \(x,M\).
The \(y^3,y^2,y,1\) coefficients of \(g\) have
\(x\)-degrees at most \(2,4,5,6\), respectively,
so the quotient for \(z\) has degree at most fifteen.
The two remaining outputs on each fibre have degree
at most \(6\cdot15=90\). Use the planar degree-one-hundred theorem on those fibres
and the fibrewise injectivity argument in the
[bounded coordinate-completion proof](LOW_DEGREE_COORDINATE_KELLER_MAPS_IN_THREE_VARIABLES.md).
The inverse descends to the base field.

If a restriction is constant at \(\zeta\ne0\), the
first and third outputs restricted to \(x=\zeta\)
have nonzero constant planar determinant because
\(R_x(\zeta,y,z)\) is a nonzero constant and
\(R_y=R_z=0\) there. Their \(z\)-degrees are three
and two with nonzero constant highest coefficients.
Such a planar pair is impossible: translate \(z\)
by a polynomial in \(y\) to complete the square
in the quadratic, writing the pair as
\[
             a_0z^3+p_2(y)z^2+p_1(y)z+p_0(y),
                 \qquad b_0z^2+q(y),\quad a_0b_0\ne0.
\]
The constant determinant term requires
\(p_1 q'\) to be a nonzero constant. Both factors
are constant, but the \(z^2\) coefficient gives
\(2b_0 p_1'-3a_0 q'=0\), a contradiction.
The constant restriction can therefore occur
only at zero. There \(g_y\) and \(\partial_x g_y\)
vanish, so the monic degree-four coefficient of
\(y^2\) in \(g\) is divisible by \(x^2\).
It also vanishes at every other root of \(a\).
There can be at most two such distinct roots.

The [exact checker](check_sextic_square_cusp_mixed_axial.py)
replays the first axial supports, the first binomial
coefficients and the critical-plane bracket in Python
with SymPy.


## Inversion in the mixed square-cusp sector

**Theorem SexticKellerSquareCuspMixedConstantsInversion.**
Every characteristic-zero sextic Keller map in the square-cusp
sector with \(ke\ne0\) is invertible.

By the preceding theorem, a noninvertible completion would have
the affine scalar \(R=a(x)z+g(x,y)\) with \(R\) constant on
\(x=0\), \(a'(0)=0\), and \(g_x(0,y)\) a nonzero constant.
Consequently the first and third components restricted to that
plane are a planar Keller pair \((L,M)\).
The axial normal form gives \(\deg_z L\le2\), \(\deg_z M\le1\).
The original \((2,3,4)\)-highest terms \(y^6\) in \(A\) and
\(e y^4\) in \(C\) survive both the target shear in \(R\)
and the translation in \(x\). Positive powers of \(x\)
could not occur in their respective coefficients at those
weights. Hence \(L(y,0)\) has degree six and \(M(y,0)\)
has degree four.

Write \(M=c(y)z+d(y)\). If \(c=0\), its Keller bracket
makes \(d'\) a nonzero constant, contrary to degree four.
If \(c\ne0\), express \(L\) as a polynomial in \(M\) over
\(k(y)\). Its derivative in \(y\) at fixed \(M\) is a
nonzero constant divided by \(c(y)\), so
\(L=H(M)+r(y)\) with \(H\in k[M]\) and \(r\in k(y)\).
The polynomial identity makes \(r\) polynomial, and
\(c r'\) is a nonzero constant. Thus \(c\) is a unit
and \(r\) is affine linear. Now \(L(y,0)=H(d(y))+r(y)\)
has degree at most one or a positive multiple of four,
whereas its degree is six. This excludes the plane and
proves inversion.

## The square-cusp sector with \(e=0\) and \(k\ne0\)

**Theorem SexticKellerSquareCuspZeroEAffineReduction.**
In the \((2,3,4)\) square-cusp frame with \(e=0,k\ne0\),
the third component is affine in \(z\), and the first
has \(z\)-degree three. If the second is not affine
in \(z\), then after translation of \(x\) the highest
\(z\)-terms are \(x^3z^3,\sigma x^2z^2,a(x)z\),
with \(\sigma\ne0\) and \(a\) univariate of degree four.
In that case a noninvertible completion has a repeated
root at zero of \(a\), where the constant term of the
third component is constant in \(y\). All other roots
have affine restrictions with nonzero slope, and there
is at most one such distinct root.

Use the same axial support as in the mixed-sector proof,
putting \(e=0\). At \(w=5\) the first form is \(T^3\),
and the other forms are
\[
        Tx^4+\alpha xyz^2+\beta z^3,\qquad
        kTx^3+\eta yz^2 .
\]
The bracket at fixed \(T=xz\) gives \(\alpha=\eta=0\)
as before. If \(\beta\ne0\), the binary image
generates \(k(T,x)\). Over \(k(T)\) this follows from
the pole orders four and three at infinity, because
the only other pole of the first function has order
three. To descend \(T\), work over an algebraic closure
of the ground field. The possible conjugates \(T'=\omega T\)
are cube-root translates, and the rational-field
identification at fixed \(T\) makes the induced change
of \(x\) fractional linear. The poles force it to preserve
zero and infinity, so \(x'=\zeta x\). Equality of the
two leading positive powers at infinity gives
\(\omega\zeta^4=\omega\zeta^3=1\). Thus \(T,x\) are
already in the binary field.

The same argument applies when the two functions
after fixing \(T\) are \(Tx^4+H(T)+bT^m/x^u\) and
\(kTx^3+cT^n/x^v\), with at least one pole at zero,
\(1\le u<4\) if \(b\ne0\), and nonnegative \(m,n\).
A first transverse tie lifts \(T,x,z\) as rational
functions, because the three highest forms have
rank at most two and specialize onto the old binary
algebra. As in the axial pole lemma, \(T'\) is
squarefree by \(T'^3=A_*\) and its weight \(w+2\).
A pole of the lifted \(x\) in the second form is
impossible, and the lowest pole term in either
remaining form would detect any pole of \(z=T'/x\).
Thus both \(x,z\) lift polynomially. The finite
sequence of shears from that lemma excludes a
genuine tie. In particular \(\beta=0\) at \(w=5\).

At \(w=6\) the possible highest tuple is
\[
 (z^3(x^3+\theta z),\quad x^2z(x^3+\sigma z),
                           \quad kxz(x^3+\rho z)),
 \qquad \rho=\xi/k .
 \tag{1}
\]
These forms generate the full binary field. Indeed with
\(t=z/x^3\) their two weight-zero ratios are
\(A^2/B^3=t^3(1+\theta t)^2/(1+\sigma t)^3\) and
\(B^7/(C/k)^8=t^{-1}(1+\sigma t)^7/(1+\rho t)^8\).
At zero the second has a simple pole and the first a
zero of order three. The only possible other zero of
the first is at \(-1/\theta\), of order two when it
is present. It cannot lie above the same place of
the field generated by the ratios, as its order
would have to be a multiple of three. The degree
at that place is one, so \(t\) is recovered; then
\(B/C\) recovers \(x\).

If \(\theta\ne0\), take rational lifts \(X,Z\) at
a first transverse tie as before. At a prime
divisor put \(n=\operatorname{ord}X,\ m=\operatorname{ord}Z\).
If \(m<0\), polynomiality of the first lift in (1)
forces \(m=3n<0\) and
\(\operatorname{ord}(X^3+\theta Z)\ge-9n\).
Polynomiality of the third lift then requires
\(\rho=\theta\), since otherwise its order is \(7n\).
It would therefore vanish to order at least \(-5n\)
on this prime. Its binary specialization
\(kx z(x^3+\theta z)\) is squarefree, a contradiction.
If \(n<0\) and \(m\ge0\), the second lift forces
\(m\ge-5n\), and the first then vanishes to order at
least \(-12n\). Its binary specialization has
multiplicities at most three, again impossible.
For these comparisons the prime can be represented
by a weighted homogeneous irreducible factor of
the appropriate lifted form. Its specialization
is nonzero and nonconstant because the form has
nonzero specialization and the grading is positive.
Thus both generators lift polynomially, and the
finite-shear argument contradicts \(\theta\ne0\).
For \(\theta=0\), if \(\rho\ne0\), the preceding
pole argument applies with \(B_*=Tx^4+\sigma T^2\);
it gives \(\rho=0\). This leaves \(A\) of \(z\)-degree
three and \(C\) with possible \(z^2\)-coefficient
only the scalar \(\chi\).

If \(\sigma=0\), at \(w=7\) the only possible
correction is \(\tau yz^2\) in the second form.
The independent first and third forms fix \(x,z\)
up to a finite extension, so \(\tau=0\).
At \(w=8\) the remaining ties are \(\lambda xz^2\)
and \(\chi z^2\), in the second and third forms.
At \(w=10\) the only remaining tie is \(\mu z^2\)
in the second. These are the binary poles just
excluded. Thus both components are affine in \(z\)
in this case.

Suppose \(\sigma\ne0\). At the higher weights use
\(B^\flat=B-\sigma A^{2/3}\) and the polynomial
\(D'=(B^3-\sigma^3 A^2)/(3\sigma^2)\), which has
highest form \(T^4(B^\flat)_*\). Here the coefficient
of \(z\) in \(B^\flat\) has \((2,3)\)-highest part
\(x^5\): subtracting \((2\sigma/3)A_2 a_3^{-1/3}\)
has weight at most eight, below ten. The terms
of nonpositive \(z\)-degree have \((2,3,4)\)-weight
at most fourteen and cannot tie \(x^5z\) for \(w>4\).
At \(w=7\), its possible \(y z^2\)-term gives
\(\tau=2\sigma b/3\). At \(w=8\), fixing \(T,C_*\)
forces the \(y^2z^2/x^2\) term \(\sigma b^2/9\)
of \(B^\flat\) to vanish. This gives \(b=\tau=0\).
The two binary poles remaining at \(w=8\) then
give \(\chi=0,\ \lambda=2\sigma c/3\), as
justified below. At \(w=9\) the \(yz^2/x\)-term
forces \(d=0\). The potential poles at weights
ten, twelve and fourteen are now only terms
proportional to \(z^2,z^2/x,z^2/x^2\) in \(B^\flat\).
They force the respective coefficients to vanish,
by the same binary lifting.

For these pole comparisons \(A_*=T^3,\ B_*=\sigma T^2\)
already recover \(T\). Over \(k(T)\), \(C_*\) and
\((B^\flat)_*\) recover \(x\) by the place at infinity:
it has ramification index one from orders three
and four, and no other place can lie above it.
At \(w=8\) the other possible pole of \(C_*\)
has order two; at subsequent ties \(C_*=kTx^3\)
has just its pole at infinity. The four leading
polynomials \(A,B,C,D'\) have rank at most two at
a first transverse tie. This follows from the
same expected-weight inequalities as for \(D\),
now with \(B\) and \(C\) interchanged in the cubic
identity. This justifies their lifts. The squarefree
\(T'\) and the leading term \(kT'X^3\) in \(C_*\)
rule out a pole of \(X\), which therefore equals
\(x\). Either polar correction at \(x=0\) then
rules out a pole of \(Z=T'/x\), by the identities in
\(C_*\) and \(D'_*\). The finite sequence of shears
excludes the tie. This proves the stated vanishing.

Thus \([z^2]B-\sigma a_3^{2/3}=O(x^{-3})\).
The centered-cubic calculation in the mixed-sector
proof gives \(a_3=x^3\), \([z^2]B=\sigma x^2\)
after translation. At \(z=\infty\), the coefficients
\[
 u=[z]B/x-2\sigma A_2/(3x^2),\qquad v=[z]C/x
\]
are algebraically dependent by the bracket of \(A,C,D'\).
Their \((2,3)\)-highest parts are \(x^4,kx^3\),
and \(v\) has at most a simple pole at \(x=0\).
The rational-curve argument of the mixed-sector
proof shows they are univariate: infinity lies
alone above its image since the \(v\)-pole there
has order three, and the only other possible
\(v\)-pole has order at most one. An affine
parameter \(x+s(y)\) on that image therefore has
constant shift, by the fixed finite poles if
present or the \((2,3)\)-bound on \(u\) if absent.
Hence \(a=[z]C=xv\) is a univariate quartic.

For the last assertion apply the
[normal-residue proposition](KELLER_LINEAR_SOURCE_COEFFICIENT_CUBIC_CLASSIFICATION.md#normal-residues-along-a-repeated-line)
to \(C=a(x)z+g(x,y)\). If all root restrictions
have nonzero affine slopes, the parallel-fibre
coordinates give degree at most five for \(y\)
and at most twelve for \(z\). Indeed \(g\) has
\(y\)-degree at most three, with the \(y^3\) coefficient
of degree at most one in \(x\). The other coefficient
bounds are three, five and six, so substitution of
\(y\) gives degree at most sixteen, and division by
\(a\) lowers it by four. The resulting planar degree
at most seventy-two gives inversion by the same
fibrewise argument. A constant root restriction
away from zero is excluded by the planar bracket:
the highest \(z\)-terms of \(A,B\) there are nonzero
constant multiples of \(z^3,z^2\), respectively.
Thus only zero can carry a constant restriction;
it must be repeated. The degree-three coefficient
of \(y^2\) in \(g\) then vanishes twice at zero
and vanishes at every other root, so at most one
other root occurs.

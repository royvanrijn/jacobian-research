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


## Inversion when both lower axial components are affine

**Theorem SexticKellerSquareCuspZeroETwoAffineInversion.**
In the square-cusp sector \(e=0,k\ne0\), any completion
in which the second component is also affine in \(z\)
is invertible.

In this case write the highest axial tuple as
\((a_3(x,y)z^3,b(x,y)z,c(x,y)z)\). Its coefficient
bounds are those of the preceding theorem. In particular
\(a_3\) has \(y\)-degree at most one, \(b=x^5+\cdots\),
and \(c=kx^4+\cdots\), where all the omitted terms
have strictly smaller \((2,3)\)-weight.
Fixing \(c z\) in the highest tuple makes
\[
                    u=b/c,\qquad v=a_3/c^3
\]
algebraically dependent by the leading determinant.
Over \(k(y)\), the pole order of \(u\) at \(x=\infty\)
is one, and the zero order of \(v\) there is nine.
The place below infinity has ramification index one.
No finite place lies above it, since the numerator of
\(v\) has \(x\)-degree three, and a finite zero would
have order at least nine. It follows that
\(k(y)(u,v)=k(y)(x)\).

Over the algebraic closure, specialize generic \(y\)
to embed \(k(u,v)\) in a rational field, as before.
Take a parameter with simple pole at the place under
\(x=\infty\). Over \(k(y)\) it is affine in \(x\), and
the leading coefficient of \(u\) at infinity,
\(1/k\), makes the affine slope a scalar. After scaling,
the parameter is \(x+s(y)\). The constant Laurent
coefficient of \(u\) at infinity involves only
\([x^4]b\) and \([x^3]c\), which are scalars by the
\((2,3)\)-bounds. It therefore makes \(s(y)\)
constant. Thus \(u,v\) are univariate in \(x\).
Since \(a_3=v c^3\) has \(y\)-degree at most one,
\(c\) must be univariate, and so must \(b=uc\).

Now write \(B=b(x)z+g(x,y)\) and \(C=c(x)z+h(x,y)\).
Choose a generic nonzero constant \(\lambda\) and use the
scalar \(C+\lambda B\). Its coefficient \(c+\lambda b\)
has degree five. A repeated root outside the common
roots of \(b,c\) must be a critical point of the
nonconstant rational function \(c/b\). This excludes
only finitely many \(\lambda\). At a common root of
\(b,c\), the restrictions of \(g,h\) cannot both be
constant, since then \(dB,dC\) would be proportional
on that whole plane. At each of these finitely many
roots, avoid the value of \(\lambda\) for which the
restriction of \(h+\lambda g\) is constant, if there
is one. For the resulting scalar the normal-residue
proposition gives an affine restriction with nonzero
slope at every root of its coefficient.

The parallel-fibre coordinates used above now apply.
The \((2,3)\)-weight bound of fourteen on \(h+\lambda g\)
permits a \(y^4\) term only with a scalar coefficient,
which vanishes at every root of the degree-five
coefficient and is therefore zero. The \(x\)-degrees
of the coefficients of \(y^3,y^2,y,1\) are at most
two, four, five and six, respectively. Hence the inverse
coordinates on each fibre have degree at most fifteen,
and the two restricted outputs have degree at most
ninety. The planar bound and the same fibrewise
injectivity argument give inversion and descend it to \(k\).


## The generic-fibre residue in the quadratic-middle case

**Theorem SexticKellerSquareCuspZeroEQuadraticMiddleInversion.**
In the square-cusp sector \(e=0,k\ne0\), any completion
with second component quadratic in \(z\) is invertible.

Suppose the completion is noninvertible.
Use the frame from the affine reduction, writing
\[
 A=x^3z^3+A_2z^2+A_1z+A_0,\qquad
 B=\sigma x^2z^2+B_1z+B_0,\qquad C=a(x)z+h(x,y).
\]
The constant restriction at zero has \(C_x(0,y,z)\)
a nonzero constant and \(C_y=C_z=0\) there. Hence the
restrictions of \(A,B\) to \(x=0\) are a planar Keller
pair, with \(z\)-degrees at most two and one. The proof
of the mixed-sector inversion shows that either
\(B_1(0,y)=0\) and \(B_0(0,y)\) is affine with nonzero
slope, or \(B_1(0,y)\) is a nonzero constant and on
that plane \(A=H(B)+\ell y+r\). In the latter case
\(\deg_y A(0,y,0)=6\) forces \(H\) quadratic and
\(B_0(0,y)\) cubic, since \(\deg H\le2\) and
\(\deg_y B_0(0,y)\le4\).

The \(y^4\)-coefficient of \(B_0\) and the \(y^3\)-coefficient
of \(B_1\) are scalars by the \((2,3)\)-bounds, so both
vanish from this plane calculation. Write \(h=h_2(x)y^2+
h_1(x)y+h_0(x)\): there is no \(y^3\) term since its
coefficient has degree at most one and vanishes twice
at zero by the constant-restriction condition.
The polynomial \(h_2\) has degree three with highest
coefficient \(k\), and \(x^2\mid h_2,h_1\).

From the preceding affine reduction,
\(B_1/x-2\sigma A_2/(3x^2)\) is univariate.
The \(y^2\)-coefficient of \(B_1\) has degree at most
one in \(x\), and vanishes at zero. Its \(x\)-term
is \(2\sigma x\), by the identity with \(A_2\), whose
leading \(y^2\)-coefficient is \(3x^2\). Thus
\[
       [y^2]B_1=2\sigma x,\qquad [y^2]A_2=3x^2,
                    \qquad [y^3]A_2=0 .
\]
Also \([y^4]A_1=3x+\delta\) for a scalar \(\delta\).

On the generic fibre \(C=t\), put \(z=(t-h)/a\).
The two resulting Laurent-rational polynomials in
\(x,y\) have bracket \(\pm\kappa/a(x)\).
Their \(y^6\)-coefficient in \(A\) and \(y^4\)-coefficient
in \(B\) are respectively
\[
    (1-xh_2/a)^3-\delta h_2/a,\qquad
                 \sigma (xh_2/a)(xh_2/a-2).
\]
These are their largest possible \(y\)-degrees.
At \(x=\infty\), the first coefficient tends to zero
and the second tends to \(-\sigma\ne0\). If the first
is nonzero, the leading \(y\)-bracket requires the
fourth power of the first coefficient divided by the
sixth power of the second to be a nonzero constant
in \(x\). This contradicts its limit. Thus the first
is zero identically. Since \(h_2/a\sim 1/x\), comparison
at infinity gives \(\delta=0\), and then \(a=xh_2\).
In particular \(x^3\mid a\), and the \(y^4\)-coefficient
of the fibre polynomial \(B\) is exactly \(-\sigma\).

The scalar residue method used in the
[parallel-coefficient argument](SEXTIC_KELLER_PARALLEL_COEFFICIENT_SCALAR.md#reduction-to-quadratic-mates)
now makes \(dx/a(x)\) rationally exact. Explicitly,
on the fixed \(C\)-fibre the identity
\[
 \partial_x(A B_y B^{-1/4})-\partial_y(A B_x B^{-1/4})
                     =\{A,B\}_{x,y} B^{-1/4}
\]
holds in the formal Laurent series at \(y=\infty\).
The \(y^{-1}\)-coefficient of \(B^{-1/4}\) is a
nonzero scalar, and a \(y\)-derivative has zero
\(y^{-1}\)-coefficient. Taking this coefficient
therefore gives a primitive of \(dx/a(x)\).
A nonzero root of the quartic \(a\) would be simple,
because \(x^3\mid a\), and would give a nonzero
residue. Hence \(a=kx^4\) and \(h_2=kx^3\), giving
\[
        C=kx^3(xz+y^2)+h_1(x)y+h_0(x),\qquad x^2\mid h_1.
 \tag{2}
\]


To finish, give \(x,y\) weights \((-2,3)\) on the
generic fibre (2), and put \(\theta=(t-h_0(0))/k\ne0\).
Since \(x^2\mid h_1\), substituting for \(z\) gives
\[
    (xz+y^2)_*=\theta x^{-3},\qquad \deg_w z\le8 .
\]
Indeed \(xz+y^2=(t-h_0-h_1y)/(kx^3)\), whose numerator
has weight zero and leading term \(t-h_0(0)\).
The identities already proved give
\[
\begin{aligned}
 A-P^3&=(A_2-3x^2y^2)z^2+(A_1-3xy^4)z+(A_0-y^6),\\
 B-\sigma P^2&=(B_1-2\sigma xy^2)z+(B_0-\sigma y^4).
\end{aligned}
\]
The coefficient of \(y\) in \(A_2-3x^2y^2\)
is divisible by \(x^2\), since it equals
\((3x/(2\sigma))[y]B_1\) and \(B_1(0,y)\) is
constant. Its constant coefficient is polynomial,
and it has no higher \(y\)-term. Thus the three
terms of \(A-P^3\) have weights at most \(16,17,15\),
respectively. The linear-in-\(y\) part of
\(B_1-2\sigma xy^2\) is divisible by \(x\);
this coefficient of \(z\) therefore has weight at
most one, and \(B_0\) has \(y\)-degree at most three.
Consequently the highest forms of the two Laurent
polynomials on the generic fibre are
\[
             A_*=\theta^3 x^{-9},\qquad
             B_*=\sigma\theta^2 x^{-6}-\sigma y^4.
\]
Their bracket is \(36\sigma\theta^3x^{-10}y^3\) of
weight twenty-nine, so it is the highest bracket of
the fibre polynomials. The required bracket
\(\pm\kappa/(kx^4)\) has weight eight. This proves the
contradiction and settles the quadratic-middle case.

**Corollary SexticKellerSquareCuspNonzeroKInversion.**
Every characteristic-zero sextic Keller map in the
cubic-outer square-cusp frame with \(k\ne0\) in the
\((2,3,4)\) reduction is invertible.

For \(e\ne0\) this is the mixed-sector inversion.
For \(e=0\) the preceding reduction makes the third
output affine, and the two cases according to the
\(z\)-degree of the second output are settled above.


## Inversion at zero \(k\) and nonzero \(e\)

**Theorem SexticKellerSquareCuspZeroKNonzeroEInversion.**
Every characteristic-zero sextic Keller map in the
cubic-outer square-cusp \((2,3,4)\) frame with
\(k=0\) and \(e\ne0\) is invertible.

Use \(A,B,C\), \(C^\flat=C-eA^{2/3}\), and
\(D=(C^3-e^3A^2)/(3e^2)\) as in the preceding
axial argument. We first extend its pole comparison
to a leading field that contains \(x^2\).

### Lifting \(x^2\) at a transverse tie

Write \(T=xz\) and \(S=x^2\).
Suppose that at weights \((2,3,w)\), \(w>4\), \(A_*=T^3\),
\(C_*=eT^2\), and the two binary forms \(B_*,(C^\flat)_*\)
have the shape
\[
 T x^4+H(T)+b_0T^m/x^u,\qquad G(T,x)+c_0T^n/x^v,
 \tag{3}
\]
where \(B_*\) has weight \(w+10\), \(H,G\) are
polynomials, \(m,n\ge0\), and \(u,v>0\) when the
corresponding coefficient is nonzero.
In the applications below, the field of the four
polynomial forms \(A_*,B_*,C_*,D_*\) is either
\(k(T,x)\) or \(k(T,S)\). In the second case the displayed
terms are even in \(x\), so \(G\) is polynomial in \(T,S\).
In both cases (3) forces \(b_0=c_0=0\).

For the proof use the first increase of the weight \(s\)
of \(y\) that ties in any of the four forms, as in the
[axial pole exclusion](#axial-pole-exclusion-with-two-nonzero-coefficients).
It occurs by \(s\le3(w+2)\), since otherwise the original
three outputs would all be independent of \(y\).
If \(d\) is the weight of \((C^\flat)_*\), then \(D_*\)
has weight \(4(w+2)+d\). The expected weights of
the brackets \((A,B,D)\) and \((B,C,D)\) exceed those
of their actual values by \(3w+14+d-s>0\)
for the positive \(d\) used here. The other two
three-output leading brackets vanish as before.
Thus specialization at \(y=0\) gives an isomorphism
of the new and old binary image rings.
The lift of \(T\) is polynomial by its cube, and
it is squarefree because it is \(xz+f(x,y)\),
with \(f\) containing no \(z\).

In the first field case a pole of the lift of \(x\)
would give a unique pole of \(T x^4\) in the second
output, so \(x\) lifts identically by weight two.
In the square-field case the corresponding argument
for \(T S^2\) gives a polynomial lift of \(S\). Its
only possible change is \(x^2+\nu y\) at \(s=4\).
If \(\nu\ne0\), that prime does not divide the lift
of \(T\), since the latter's \(z\)-coefficient is \(x\).
Any nonzero polar term in \(B_*\) or \(D_*\) would
then have a pole on this prime. So the lift of \(S\)
is \(x^2\). In either case a nonzero polar term forces
\(x\mid T'\) for the lifted \(T'\), by the same
comparison on \(x=0\). The resulting polynomial
\(Z=T'/x\) has highest weight \(w\) and specializes
to \(z\). Its triangular source shear removes the tie.
As before only finitely many shear weights
\(2i+sj=w\) with \(j>0,\ s>3\) are possible, proving
the assertion.

For the field calculations in (3), we will use the
degree at the place below \(x=\infty\), over \(k(T)\).
With pole orders four and either six or two there,
the ramification index is at most two. Any other
pole of the first function at zero with order one or
three cannot lie above the same place, whose pole
order for that function is four or two. If the
degree is two, its involution must preserve infinity
and the only other possible pole of either function,
at zero. This involution is \(x\mapsto -x\).
In particular an odd term excludes degree two when
there is such a finite pole. We also use directly
that the normalized field is \(k(T,S)\) when \(S^2\)
and a nonzero odd function of \(S\) are given.

### Increasing the axial weight

At \((2,3,4)\) use the coefficient notation
\[
\begin{aligned}
 A &=\theta z^4+(x^3+bxy+cx^2+dy+fx+g_0)z^3
                              +A_2z^2+A_1z+A_0,\\
 B &=\beta z^3+(\alpha xy+\sigma x^2+\tau y+\lambda x+\mu)z^2
                              +B_1z+B_0,\\
 C &=(e x^2+\eta y+\xi x+\chi)z^2+C_1z+C_0 .
\end{aligned}
\]
Here the \((2,3)\)-highest parts of \(A_2,A_1,A_0\)
are \(3x^2y^2,3xy^4,y^6\), those of \(B_1,B_0\)
are \(x^5,x^4y^2\), and those of \(C_1,C_0\) are
\(2exy^2,x^6+ey^4\). All lower ordinary terms
satisfy the strict weight bounds from the
\((2,3,4)\)-reduction.

Increase the \(z\)-weight to \(9/2\). The highest
part of \(C^\flat\) is \(x^6+(\eta-2eb/3)yz^2\),
so fixing \(T,B_*=Tx^4\) gives \(\eta=2eb/3\).
At \(w=14/3\) the only possible correction to \(x^6\)
is \(-2e\theta z^3/(3x)\), from the \(z^4\)-term
of \(A-P^3\) in the linear binomial term
\(-2e(A-P^3)/(3P)\).
If it is present, \(T, S^2\), and
\(S^3-2e\theta T^3/(3S^2)\) generate \(k(T,S)\).
The pole comparison gives \(\theta=0\).

At \(w=5\), denote by \(s_0\) the coefficient of
\(y^3z^2\) in \(A_2z^2\). The weight-deficit-one
and weight-deficit-two terms of \(A-P^3\) from
its top weight twenty-one are \(bxy z^3\) and
\(cx^2z^3+s_0y^3z^2\).
Using the linear and quadratic binomial terms,
and expanding \(1/P\) around \(1/T\), gives
\[
 \begin{aligned}
 B_*&=T x^4+\beta z^3+\alpha xy z^2,\\
 (C^\flat)_*&=x^6+(\xi-2ec/3)xz^2+\gamma_0 x^2 y z
          -\tfrac{2e}{3}(s_0-b)y^3z/x
          +\tfrac{e b^2}{9}y^2z^2/x^2,
 \end{aligned}
\]
where \(\gamma_0=[x^2y]C_1\).
If \(\alpha\ne0\), fix \(T,B_*\) and solve for \(y\).
In decreasing order, the coefficients of the third,
second and first powers of \(B_*\) in \((C^\flat)_*\)
would depend on nonzero powers \(x,x^{-2},x^2\),
respectively. Constancy on this level forces them
all to vanish. The remaining expression is
nonconstant in \(x\) since it has the \(x^6\) term,
a contradiction. If \(\alpha=0\) the free \(y\)
directly forces the three coefficients to vanish.
Thus \(\alpha=\gamma_0=b=s_0=\eta=0\).
The remaining poles are \(\beta T^3/x^3\) in the
second output and \((\xi-2ec/3)T^2/x\) in the
normalized third. If either occurs, the pole order
argument gives \(k(T,x)\), and (3) excludes both.
In particular \(\xi=2ec/3\).

The coefficient \(a_3=[z^3]A\) is now
\(x^3+cx^2+dy+fx+g_0\), and \(A\) has \(z\)-degree
three. Use the Laurent expansion
\[
 C^\flat=R_2 z^2+R_1z+R_0+\cdots,\quad
 R_2=[z^2]C-e a_3^{2/3},\quad
 R_1=C_1-\tfrac{2e}{3}A_2a_3^{-1/3}.
\]
Powers of \(a_3\) are expanded from \(x^3\).
The \((2,3)\)-highest part of \(R_0\) is \(x^6\):
the leading terms from \(A_1,A_2\) contribute
\(2ey^4-ey^4=ey^4\) in \(eA^{2/3}\).
The negative powers of \(z\) have weight less than
twelve for \(w>4\), by their \((2,3,4)\) bounds.
Also \(R_1\) has \((2,3)\)-weight at most six:
the weight-eight terms cancel, and the weight-seven
possibility \(x^2y\) was excluded at weight five.
For this last bound the only weight-nine term
allowed in \(A_2-3x^2y^2\) was \(s_0y^3\), now zero.
Both \(R_1\) and its subsequent expansions have
\(y\)-degree at most two once \(d=0\).

At \(w=11/2\) the weight-one term \(-2ed\,y/(3x)\)
in \(R_2\) would give a \(y\)-dependent term in
\((C^\flat)_*\) while \(A_*,B_*\) fix \(x,z\).
Thus \(d=0\). At \(w=6\), the second form is
\(Tx^4+\sigma T^2\). In the third the possible
term proportional to \(y^2z\) again vanishes,
leaving \(x^6+\gamma x^3z+[x^0]R_2\, z^2\).
With \(S=x^2\) these give
\(TS^2+\sigma T^2,\ S^3+\gamma T S+[x^0]R_2\,T^2/S\),
which recover \(S\). Therefore (3) gives \([x^0]R_2=0\).
Subtract \((\sigma/e) C\) from the second output
to cancel its \(x^2z^2\) term. This target shear
does not affect its leading \(z\)-coefficient \(x^5\)
in weights \((2,3)\).

If \(\gamma\ne0\), at \(w=7\) the leading normalized
third output is \(\gamma T x^2\), forcing the possible
\(yz^2\) term of the second output to vanish.
At \(w=8\), the two forms at fixed \(T\) are
\[
 T x^4+b_1T^2/x,\qquad \gamma T x^2+[x^{-1}]R_2\,T^2/x^3.
\]
Nonzero poles give the full field \(k(T,x)\) by
the degree and involution test, and are eliminated
by (3). At \(w=10\) the two forms at fixed \(T\)
are \(TS^2+b_2T^2/S,\ \gamma T S+[x^{-2}]R_2\,T^2/S^2\).
They recover \(S\) by pole orders two and one
at infinity and the bound of one for the other
pole of the second output. Both poles are again
eliminated.

If \(\gamma=0\), at \(w=7\) the normalized third
output can have an \(xy z\) term from the
weight-five part of \(R_1\). At fixed \(T\) the
two \(y\)-coefficients in the second and normalized
third forms are proportional to \(x^{-2}\) and
\(1\), respectively. The same bracket comparison
as at weight five eliminates both. A correction
\([x^{-1}]R_2\,T^2/x^3\) would then give the full
field and a pole, so it vanishes too.
Translate \(x\) to center \(a_3\); the preceding
vanishings \([x]R_2=[x^0]R_2=[x^{-1}]R_2=0\)
give \(a_3=x^3+F x\), as in the mixed-sector
calculation. This translation preserves the
weighted bounds and the identities of highest forms.
At \(w=8\) the normalized third form is
\[
 x^6+\gamma_2 x^2z+\epsilon_2 y^2z/x+
                                [x^{-2}]R_2\,z^2/x^2.
\]
The \(y\)-term vanishes since the second form is
binary. At fixed \(T\), the potential remaining poles
are \(b_1 T^2/x\) in the second and
\([x^{-2}]R_2\,T^2/x^4\) in the normalized third.
If there is a pole and either \(b_1\) or \(\gamma_2\)
is nonzero, the involution test gives \(k(T,x)\).
Otherwise the forms give \(S^2\) and
\(S^3+[x^{-2}]R_2\,T^2/S^2\), recovering \(S\).
Thus the poles vanish by (3).

In both cases the four coefficients of \(R_2\)
through \(x^{-2}\) now vanish. The centered-cubic
calculation gives \(a_3=x^3\), \([z^2]C=e x^2\)
after translation.
When \(\gamma=0,\gamma_2\ne0\), the last possible
pole of the second form at \(w=10\) is \(b_2 T^2/x^2\)
and the normalized third form is \(\gamma_2 T x\).
Thus (3) applies directly. When
\(\gamma=\gamma_2=0\), the weight-three term of
\(R_1\) is proportional to \(y\). Its tie at
\(w=9\) is excluded by fixing \(T,B_*\).
At \(w=10\) the last possible pole \(b_2 T^2/S\)
in the second output accompanies \(S^3+\text{scalar}
\cdot T\) in the normalized third. Over \(k(T)\)
the orders two and three at infinity and the smaller
other pole of the second output recover \(S\), so
(3) applies. This completes the reduction to
an affine middle output with \(B_1\) of highest
\((2,3)\)-part \(x^5\).

### The coefficient of the affine middle output

In the last frame put \(u=B_1/x,\ v=R_1/x\).
The bracket at \(z=\infty\) gives algebraic dependence
of \(u,v\) as in the mixed-sector case, with \(v\)
having \(x\)-degree at infinity at most two.
If \(v\) is nonconstant, use \(L=v\). If \(v\)
is constant, use \(L\) equal to the \(z^0\)-coefficient
of \(C^\flat-v A^{1/3}\). The generic leading
determinant of \(A,B,C^\flat-v A^{1/3}\) gives
algebraic dependence of \(u,L\): the three
\(z\)-degrees are \(3,1,0\) and the bracket is
unchanged. The formula
\[
         L=C_0-\frac{2e}{3}A_1/x+\frac{e}{9}A_2^2/x^4
                              -\frac v3 A_2/x^2
\]
shows its \((2,3)\)-highest part is \(x^6\).
In either case \(u,L\) have no finite poles in
\(x\) other than possibly zero, over \(K=k(y)\).

At infinity \(u\) has pole order four.
When \(L=v\), a nonconstant \(v\) must also be
nonconstant in \(x\), by algebraic dependence with
\(u\). It has either pole order at most two or,
after subtracting its finite value there, zero
order at most two, since its denominator divides
\(x^2\). When using the other \(L\), it has pole
order six at infinity. In both cases the index at
infinity over \(K(u,L)\) is at most two. There is
no other place above it, since the only other
possible pole of \(u\) has order at most one.
Thus the rational-curve parameter from the mixed
case pulls back to a polynomial of degree one or
two in \(x\). Its leading coefficient can be made
one by the constant leading coefficient of \(u\).
If \(u\) or \(L\) has a finite pole, then the
parameter must have fixed scalar value at \(x=0\);
in degree two the two preimages of that value
must coincide there. Hence that parameter is
univariate in \(x\). If neither has a finite pole,
\(u=U(Q)\) is a monic polynomial in the parameter.
In degree one its \(x^3\)-coefficient forces a
constant shift. In degree two the \(x^3,x^2\)
coefficients and their \((2,3)\) bounds give
\[
                    Q=x^2+\alpha x+\beta y+\gamma
\]
with scalar coefficients and \(U\) monic quadratic.

If \(\beta\ne0\), the affine coefficient
\(B_1=xU(Q)\) has a multiplicity-one component
\(x=0\). The residue of \(dx\wedge dy/B_1\) there
is \(dy/U(\beta y+\gamma)\), with an affine pole.
This contradicts the
[linear-source residue restriction](KELLER_QUARTIC_TRIPLE_BOUNDARY_AND_LINEAR_COMPONENTS.md#a-residue-restriction-for-components-linear-in-a-source-variable).
Therefore \(B_1\) is univariate in all cases.

The three highest \(z\)-terms and the weight bound
on \(B_0\) are now exactly those used in the
mixed-sector affine-scalar proof. The normal-residue
restriction and parallel-fibre coordinates give
inversion unless a repeated root produces a
constant scalar restriction. A nonzero such root
would give the impossible planar \((3,2)\)-pair
with constant highest coefficients. At zero the
planar degree mismatch six versus four excludes
the constant restriction as in the mixed-sector
inversion. This completes the proof.


## The first axial forms when both coefficients vanish

**Theorem SexticKellerSquareCuspBothZeroAxialReduction.**
In the square-cusp frame with \(k=e=0\), the highest
tuple at \((2,3,5)\) is \(((xz)^3,x^5z,x^6)\).
At \((2,3,6)\) the highest tuple is
\[
           (x^3z^3,\ x^5z+\sigma x^2z^2,\ x^6+\delta x^3z).
\]
In particular the first component has \(z\)-degree
three and the third has \(z\)-degree at most one.

Use the support from the preceding section, putting
\(k=e=0\). At \(w=9/2\), the \(\eta yz^2\) tie in
\(C\) vanishes by independence of the two binary
forms in \(A,B\). At \(w=5\) the possible other
forms are \(Tx^4+\alpha xyz^2+\beta z^3\) and
\(x^6+\xi xz^2+\gamma_0 x^2 y z\).
The bracket at fixed \(T=xz\) forces
\(\alpha=\gamma_0=0\) by the \(y\)-coefficient
comparison.

Suppose either \(\beta\) or \(\xi\) is nonzero.
The field generated by the forms at \(w=5\), after
adjoining \(T\), is \(k(T,x)\) by the pole-order
argument of the preceding section. Before adjoining
\(T\), the field is \(k(x^3,z/x)\), over the algebraic
closure. Indeed the extension given by \(T^3=A_*\)
has possible conjugates \(T\mapsto\omega T\),
\(\omega^3=1\), which induce fractional linear
changes in \(x\) over \(k(T)\). The pole orders of
the fixed functions at infinity and zero force
\(x\mapsto d x\). Their leading positive powers
give \(\omega d^4=1,\ d^6=1\). The nonzero
pole in either function forces \(d^3=1\).
Conversely simultaneous scaling of \(x,z\) by
a cube root of unity fixes all three forms.
This proves the field claim.

At a first transverse tie lift \(H_1=x^3,\ M=z^3\)
and \(Q=z/x\) as rational functions. The specialization
is an isomorphism of the old and new image rings,
since the leading determinant is zero and the old
binary rank is two. For any prime of the source
write \(n=\operatorname{ord}(H_1)/3\) and
\(m=\operatorname{ord}(M)/3\) for the orders of the
lifts, so \(m-n=\operatorname{ord}(Q)\) is integral.
The first form gives \(n+m\ge0\).
If \(n<0\), the \(x^6\)-term gives a unique pole
in the third form, since the other term \(xz^2\)
has order \(n+2m>0\). Thus \(H_1\) is polynomial.
If \(\beta\ne0\), \(m<0\) gives a unique
pole in the second form. If \(\beta=0,\xi\ne0\),
use the weight-six lift \(H_1\). It either remains
\(x^3\) or acquires \(\nu xy\) at weight four or
\(\nu y\) at weight six. In the changed case \(H_1\)
is squarefree, so \(m<0\) gives \(n+2m<0\).
In the unchanged case \(0\le n\le1\) is integral
and \(m-n\) is integral, with the same conclusion.
Either would give a unique pole in the third form.
Thus \(M\) is polynomial. If \(H_1\) had changed
it would be squarefree, and \(M=H_1 Q^3\) together
with regularity would make \(H_1\mid M\). This is
impossible on specialization to \(x^3,z^3\).
Therefore the lift of \(H_1\) remains \(x^3\), and
\(Z=xQ\) is polynomial by regularity of its cube.
As in the previous axial argument its shear removes
the first tie; only finitely many shear weights
are available. This excludes both \(\beta,\xi\).

At \(w=6\) the only possible \(y\)-term in the
highest tuple is \(y^2z\) in the third output.
The first two binary forms are independent, so it
vanishes. The tuple is
\[
  \bigl(z^3(x^3+\theta z),\
        x^2z(x^3+\sigma z),\
        x^6+\delta x^3 z+\chi z^2\bigr).
 \tag{5}
\]
Put \(t=z/x^3,\ S=x^2,\ f=1+\delta t+\chi t^2\).
The weight-zero ratios \(U=A_* / C_*^2\) and
\(V=B_*^3/C_*^4\) are
\[
         U=\frac{t^3(1+\theta t)}{f^2},\qquad
         V=\frac{t^3(1+\sigma t)^3}{f^4}.
\]
Both vanish to order three at zero. If the
ramification index of \(t=0\) in their generated
field is one, no other point lies above its place:
the only possible other zeros of \(U\) have order
one. Thus the ratio field is \(k(t)\). If the index
is three, \(U\) is a local parameter downstairs,
and \(V/U\) is a power series in \(U\).
Its terms of \(t\)-degree one and two give
\[
       \delta=(3\sigma-\theta)/2,\qquad
       \chi=3(\sigma-\theta)^2/8.
\]
Subtract its constant and the multiple of \(U\)
matching its \(t^3\)-coefficient. The next two
coefficients, up to the nonzero factors
\(-(\sigma-\theta)^2/64\) and
\((\sigma-\theta)^2/64\), are
\[
  9\sigma^2+38\theta\sigma-95\theta^2,\qquad
  27\sigma^3+105\theta\sigma^2-259\theta^2\sigma+31\theta^3.
\]
If \(\sigma\ne\theta\) these cannot both vanish:
for \(\theta=0\) the first suffices, and otherwise
divide by the powers of \(\theta\). The second
polynomial has remainder \(64r-64\) on dividing by
the first in \(r=\sigma/\theta\), and \(r=1\)
does not solve the first. Thus the exceptional
ratio field requires \(\sigma=\theta,\ f=1+\theta t\).
Then \(V=U=t^3/(1+\theta t)\) and the index is three.
In every other case the field of (5) is \(k(t,S)\),
since \(B_*/C_*\) recovers \(S\) once \(t\) is known.

In this last field consider first \(\theta\ne0\).
At a transverse tie let \(n,m\) be the orders of
the rational half-powers \(x,z\), computed from
the respective lifts of \(S,z^2\). If \(m<0\),
the first form requires \(m=3n<0\) and
\(\operatorname{ord}(x^3+\theta z)\ge-9n\).
The second then requires \(\sigma=\theta\).
For the third to be regular, \(f\) must vanish at
\(-1/\theta\). Since the exceptional case was
removed, this requires \(\chi\ne0\), and gives
order at least \(-6n\ge3\) for the third form.
Its binary specialization has multiplicities at
most two, a contradiction. If \(n<0,m\ge0\),
the second form forces \(m\ge-5n\), and the
first has order at least \(-12n\ge6\), again
contradicting its binary multiplicities, at most
three. As in the earlier lift, any homogeneous
prime in either argument specializes nontrivially
with positive weight. Thus the lifts of \(S,z^2\)
are polynomial. The lift of \(T\) is polynomial
by \(T^2=S z^2\). If \(S'=x^2+\nu y\) changed
at weight four, its divisibility into \((T')^2\)
would contradict the coefficient \(x\) of \(z\)
in \(T'\). So \(S'=x^2\), and \(T'/x\) is
polynomial by its square. The finite shear
argument excludes the tie. Hence \(\theta\ne0\)
is possible only in the exceptional ratio case.

If \(\theta=0,\chi\ne0\), we are also in the
field \(k(T,S)\). Use the identities
\(A_*=T^3\), \(B_*=T S^2+\sigma T^2\), and
\(C_*=S^3+\delta T S+\chi T^2/S\) at a first
transverse tie. The squarefree \(T'\) excludes
a pole of \(S'\) by the second form. The last
term of the third forces \(S'=x^2\) and \(x\mid T'\),
by the pole comparison in (3). Both source
coordinates lift, so this tie is likewise impossible.
Thus \(\chi=0\) when \(\theta=0\).

Finally, in the exceptional case with \(\theta\ne0\),
the leading relation \(A_* C_*^2=B_*^3\) survives the
first transverse tie. Unique factorization writes
the new forms as \(D p^3,D p q^2,D q^3\) with
coprime \(p,q\). On specialization \(y=0\),
\(D\) specializes to \(x^3+\theta z\) up to a
unit: its specialization divides this irreducible
polynomial, and a scalar specialization would
force the first form to be a cube. Thus \(p,q\)
specialize to \(z,x\), after constant rescalings.
They have weights six and two. Hence \(q=x\).
For \(3<s<4\), neither \(p\) nor \(D\) can acquire
a \(y\)-term of weight six, so no first transverse
tie occurs there. At \(s=4\) the ordinary term
\(y^6\) in \(A\) has weight twenty-four and must
tie. The displayed factorization at that weight
has \(p=z+\mu xy\), \(D=x^3+\theta z+\nu xy\),
because these are the only extra monomials of
weight six. It would give \(A_*|_{x=0}=\theta z^4\),
contrary to the term \(y^6\) at weight twenty-four.
This excludes the exceptional common factor and
proves the theorem. The [exact checker](check_sextic_square_cusp_mixed_axial.py)
reproduces the weight-five binomial identity and the
coefficient comparisons for the projective ratio field.


## The double-zero square cusp with a nonzero axial slope

**Theorem SexticKellerSquareCuspNonzeroAxialSlopeInversion.**
In the square-cusp sector with \(k=e=0\), write the \((2,3,6)\)-highest
tuple as \((x^3z^3,x^5z+\sigma x^2z^2,x^6+\delta x^3z)\).
If \(\delta\ne0\), the sextic Keller map is invertible.

Work over an algebraic closure. All coordinate changes retained in
the proof are affine. Write
\[
\begin{aligned}
 A&=a_3z^3+A_2z^2+A_1z+A_0,\\
 B&=(\sigma x^2+\tau y+\lambda x+\mu)z^2+B_1z+B_0,\\
 C&=c z+h .
\end{aligned}
\]
The first axial reduction and the \((2,3,4)\) support give
\[
 a_3=x^3+p xy+q x^2+r y+s_1 x+s_0 .
\]
At \((2,3)\) the highest parts of \(A_2,A_1,A_0,B_1,B_0,c,h\)
are \(3x^2y^2,3xy^4,y^6,x^5,x^4y^2,\delta x^3,x^6\),
respectively. In particular \(c\) has \(y\)-degree at most one
and \(h\) has \(y\)-degree at most three. The coefficients of \(y^3,y^2\)
in \(h\) have \(x\)-degrees at most one and two, respectively.

### Regularity of the axial field lift

Put \(T=xz,\ S=x^2\). At weights \((2,3,w)\), \(w>6\), the
highest forms of \(A,C\) are \(T^3,\delta TS\).
If \(\sigma=0\) put \(N=B\).
If \(\sigma\ne0\), put \(B^\flat=B-\sigma A^{2/3}\)
using the binomial branch with highest term \(T^2\).
The leading form of the polynomial
\[
 N=(B^3-\sigma^3 A^2)/(3\sigma^2)
\]
is \(T^4(B^\flat)_*\) when \(B^\flat\) has strictly smaller
weight than \(B\).

We claim that a binary form \(N_*=T^\ell(TS^2+bT^m/x^u)\)
with \(b\ne0\), \(u>0\), \(m\ge0\), \(\ell=0\) or \(4\) as above,
and inner weight \(w+10\), is impossible. After adjoining \(T\),
the generated field is \(k(T,x)\) for odd \(u\) and \(k(T,S)\)
for even \(u\). In the case \(\ell=4\), \(T\) is already in the
field by \(A_*=T^3,\ B_*=\sigma T^2\). For \(\ell=0\), a
nontrivial conjugation over the generated field would send
\(T\) to \(\zeta T\) and \(S\) to \(\zeta^{-1}S\), where
\(\zeta^3=1\). For odd \(u\) it also sends \(x\) to a constant
multiple of \(x\), by its square. In either parity the highest
positive power of \(x\) in \(N_*\) is then multiplied by
\(\zeta^{-1}\), which is incompatible with fixing \(N_*\).
Thus \(T\) belongs to the field in both cases.

Raise the \(y\)-weight to the first tie, including \(N\) when
\(\ell=4\). If the tie weight is \(s\), then \(s\le3(w+2)\),
since at least one original component depends on \(y\).
All the leading triple determinants vanish. For \(A,B,C\)
use the constant determinant and the positive expected weight.
For \(\ell=4\), \(J(A,B,N)=0\), and the other two determinants
with \(N\) are constant multiples of \(B^2\) and \(A\).
The expected weights exceed their actual weights by
\(4w+20-s>0\), using the weights \(3w+6,2w+4,w+6,5w+18\)
for \(A,B,C,N\). Specialization at \(y=0\) maps the new image
ring onto the old one. Both have transcendence degree two,
so this map is an isomorphism. It preserves the grading.

Write primes on the lifts of \(T,S\) and the leading polynomials.
Since \((T')^3=A_*'\), the lift \(T'\) is polynomial. Its
weight and specialization imply \(T'=xz+f(x,y)\).
It is squarefree, and if reducible it is \(xZ\) with
\(Z=z+f/x\). For odd \(u\), the lift \(X\) of \(x\) has
\(C_*'=\delta T'X^2\), so \(X\) is polynomial and its weight
forces \(X=x\). Polynomiality of \(N_*'\) then forces \(x\mid T'\).
For even \(u=2j\) rewrite the identity as
\[
 N_*'=\delta^{-2}(T')^{\ell-1}(C_*')^2+
                   b\delta^j\frac{(T')^{\ell+m+j}}{(C_*')^j}.
\]
At a prime dividing \(C_*'\) but coprime to \(T'\), the second
term has a pole and the first is regular. Thus every prime of
\(C_*'\) divides \(T'\). If \(T'\) were irreducible, \(C_*'\)
would specialize to a power of \(xz\) up to scalar, whereas it
specializes to \(\delta x^3z\). Thus \(T'=xZ\) and
\(C_*'=\delta x^3 Z\), as follows by comparing multiplicities
on specialization. The lifts are \(S'=x^2,\ T'=xZ\).

In either parity the shear \(z\mapsto Z\) removes the first tie.
It preserves the binary specializations and gives the same
identities at any subsequent tie. The possible shear weights
satisfy \(2i+sj=w\) with \(i\ge0,j>0,s>3\), so there are only
finitely many. The first tie would have to increase strictly
after each shear, which proves the claim.

### Axial coefficient comparison

For \(\sigma=0\), at \(w=7\) the term \(\tau yz^2\)
vanishes by the leading determinant at fixed \(T\), since
\(C_*=\delta Tx^2\). The next ties of the quadratic
coefficient in \(B\) are \(\lambda xz^2\) at \(w=8\)
and \(\mu z^2\) at \(w=10\). These are the odd and even
poles excluded above. Thus \(B,C\) are both affine in \(z\).

For \(\sigma\ne0\), the \(z^2\)-coefficient in \(B^\flat\) is
\(\sigma x^2+\tau y+\lambda x+\mu-\sigma a_3^{2/3}\).
The \(z\)-coefficient has highest \((2,3)\)-part \(x^5\),
since the binomial correction \(2\sigma A_2/(3a_3^{1/3})\)
has weight at most eight. The terms of \(z\)-exponent \(j\le0\)
in the binomial series have \((2,3)\)-weight at most \(12-4j\);
this follows by expanding from \(a_3z^3\) at \((2,3,4)\).
Those terms and \(B_0\) cannot tie the weight \(w+10\).

At \(w=7\) the leading determinant of \(A,N,C\) at fixed
\(T\) gives \(\tau=2\sigma p/3\).
At \(w=8\) the weight-two part of the quadratic coefficient is
\((\lambda-2\sigma q/3)x+\sigma p^2y^2/(9x^2)\).
The determinant forces \(p=0\), hence \(\tau=0\), and the
odd pole claim gives \(\lambda=2\sigma q/3\).
At \(w=9\) the term \(-2\sigma r y/(3x)\) forces \(r=0\).
The remaining discrepancies through order \(x^{-2}\) in the
quadratic coefficient tie at \(w=10,12,14\), in order, and
are excluded by the even, odd and even pole cases. Thus
\[
     a_3^{2/3}-(\sigma x^2+\lambda x+\mu)/\sigma=O(x^{-3}).
\]
Translate \(x\) to center \(a_3\). For \(a_3=x^3+f_1x+g_1\),
the estimate says its square differs from
\((x^2+2f_1/3)^3\) by a polynomial of degree at most one.
Comparison at \(x^3,x^2\) gives \(g_1=f_1=0\).
Hence \(a_3=x^3,\ [z^2]B=\sigma x^2\).
The translation preserves all the support bounds above.


### The affine coefficient when \(\sigma=0\)

Write \(B=bz+B_0\). The \(z=\infty\) leading determinant,
at fixed \(cz\), makes \(u=b/c,\ v=a_3/c^3\) algebraically
dependent. Over \(k(y)\) the orders at \(x=\infty\) are a
pole of order two in \(u\) and a zero of order six in \(v\).
For the place below infinity the ramification index is at most
two. Any other point over that place would be a pole of \(u\)
and a zero of \(v\), so it would be a root of \(c\) and a zero
of \(v\) of order at least three. At a root of \(c\), \(v\)
has order at most \(3-3=0\). Hence the place has no other point
above it and the field degree is one or two.

The field \(k(u,v)\) is rational by specialization of \(y\)
at a scalar retaining \(u\) nonconstant, followed by Lüroth.
A parameter \(Q\) with a simple pole at the place below
infinity pulls back to a polynomial in \(x\) of degree one or
two over \(k(y)\). Its leading coefficient is normalized to one
using the constant leading coefficient of \(u\). In degree one,
the \(x\)-coefficient in the Laurent expansion of \(u\) makes
the shift constant, since \([x^4]b,[x^2]c\) are scalars.
In degree two \(u\) has a simple pole in \(Q\) at infinity.
Its \(x\)- and constant Laurent coefficients then give
\(Q=x^2+\alpha x+\beta y+\gamma\) with scalar coefficients,
using that \([x^3]b,[x]c\) are at most affine in \(y\).

If \(\beta\ne0\), write \(v=V(Q)\). The zero of \(V\) at
infinity has order three. Over \(k(x)\) this gives a zero of
order three in \(v\) at \(y=\infty\). Since \(c,a_3\) both
have \(y\)-degree at most one, \(a_3\) is univariate and \(c\)
has \(y\)-degree one. Thus there is one finite pole of \(v\)
in \(y\), of order three, and \(V=\nu/(Q-\eta)^3\) with
scalars \(\nu,\eta\). We obtain \(c=L(x)(Q-\eta)\) for a
linear polynomial \(L\), by the \(x\)-degree of \(c\).
The simple line \(L=0\) has residue of \(dx\wedge dy/c\) a
nonzero scalar multiple of \(dy/(\beta y+\text{constant})\).
This has a nonzero residue at its pole, contrary to the
[linear-source residue proposition](KELLER_QUARTIC_TRIPLE_BOUNDARY_AND_LINEAR_COMPONENTS.md#a-residue-restriction-for-components-linear-in-a-source-variable).
Thus \(u,v\) are univariate. The equality \(a_3=v c^3\)
and its \(y\)-degree bound give that \(c\), and hence \(b\),
are univariate as well.

As in SexticKellerSquareCuspZeroETwoAffineInversion, choose
a generic \(R=C+\lambda B\). Its univariate affine coefficient
has degree five. Its roots outside the common zero set of \(b,c\)
are simple for generic \(\lambda\). On a common root plane
\(B_0,h\) cannot both be constant, by independence of \(dB,dC\).
Thus for generic \(\lambda\) the restrictions of the scalar
term are nonconstant at all these roots, and the
[normal-residue proposition](KELLER_LINEAR_SOURCE_COEFFICIENT_CUBIC_CLASSIFICATION.md#normal-residues-along-a-repeated-line)
makes all the restrictions affine with nonzero slope.
The same inverse fibre degree bound of fifteen applies:
the scalar term has weight at most fourteen and therefore
no \(y^4\)-term, by its restrictions at the roots. The remaining
two outputs on each fibre have degree at most ninety, giving
inversion by the planar degree bound and fibrewise injectivity.

### The affine coefficient when \(\sigma\ne0\)

The \(z=\infty\) leading determinant for \(A,N,C\) now makes
\[
           u=B_1/x-2\sigma A_2/(3x^2),\qquad v=c/x
\]
algebraically dependent. In fact \(N\) has top term \(T^5u\),
and the determinant at fixed \(T\) has expected \(z\)-degree
eight versus actual four. At infinity in \(x\), \(u,v\) have
poles of orders four and two, with leading coefficients \(1,\delta\).
Their only possible finite pole is at zero, where the orders
are at most two and one. Take a parameter \(Q\) of \(k(u,v)\)
with simple pole at the place below \(x=\infty\), again using
specialization and Lüroth.

If \(x=0\) is not above that place, the pullback \(Q\) is a
monic polynomial of degree one or two in \(x\), after scaling
using \(u\). If \(u\) or \(v\) has a finite pole, \(Q(0)\)
is scalar. In degree two the other preimage of that finite
pole would also have to be zero, so \(Q\) is univariate in
this case. If both functions have no finite pole, they are
polynomials in \(Q\). Comparing the \(x^3\)-coefficient
of \(u\) in degree one makes the shift scalar. In degree
two its \(x^3,x^2\)-coefficients give
\(Q=x^2+\alpha x+\beta y+\gamma\) as above. Here
\(v\) is affine in \(Q\), so if \(\beta\ne0\) the simple
component \(x=0\) of \(c=xv\) has residue with a simple
pole in \(y\), which is impossible.

If \(x=0\) also lies above the pole place of \(Q\), then
\(v\) has a simple pole there, is affine in \(Q\), and has
pole order two at \(x=\infty\). Thus
\(Q=x^2+\alpha(y)x+\beta(y)+\gamma(y)/x\), after
normalization, and \(u\) is polynomial of degree two in \(Q\).
The \(x^3,x^2\)-coefficients of \(u\) again make
\(\alpha\) constant and \(\beta\) at most affine in \(y\).
Because \(c=xv\) has \(y\)-degree at most one,
\(\gamma\) is also at most affine. This gives
\(c=G(x)+(\beta_1 x+\gamma_1)y\), with scalar
\(\beta_1,\gamma_1\). If \(\beta_1\ne0\), take the simple
component given by the primitive part of \(c\) as a polynomial
in \(y\) over \(k[x]\). Its function field is \(k(x)\), and
the residue on it is \(\pm dx/(\beta_1 x+\gamma_1)\).
This has no rational primitive, even if its pole lies outside
the affine normalization. Therefore \(\beta_1=0\).

We have shown that \(c\) is univariate or is the graph
coefficient \(G(x)+\gamma_1 y\) with \(\gamma_1\ne0\).
In the graph case the residue \(dx/\gamma_1\) must have a
primitive polynomial in \(h\bmod c\). It follows that
\(h\bmod c=\ell x+\mu_0\) with \(\ell\ne0\). Write
\(h=cs_2+\ell x+\mu_0\). The coordinates \((C,c,z+s_2)\)
recover \(x\) in degree two and \(y\) in degree six.
The polynomial \(s_2\) has \((2,3)\)-highest part
\(x^3/\delta\); its other monomials have weight below six.
Thus the recovered \(z\) has degree at most eight. The
bounded coordinate-completion theorem gives inversion.

### Repeated roots of the univariate coefficient

Put \(a=c(x)\), of degree three. If every \(h(\zeta,y)\) at a
root of \(a\) is affine with nonzero slope, use the same parallel
fibre coordinates as above. At fixed \(C=t\) the recovered
\(y\) has degree at most four and \(z\) has degree at most
ten, using the coefficient bounds on \(h\).
Otherwise the normal-residue proposition gives a repeated root
\(\zeta\) with constant restriction and \(h_x(\zeta,y)\) a
nonzero scalar. The restrictions of \(A,B\) there are a planar
Keller pair. For \(\zeta\ne0\) their \(z\)-degrees are three
and two with nonzero scalar leading coefficients, impossible
by the planar bracket comparison in
SexticKellerSquareCuspMixedConstantsAffineReduction.
So \(\zeta=0\).

On the plane \(x=0\), the same affine-component argument as
in SexticKellerSquareCuspZeroEQuadraticMiddleInversion gives
\(B_1(0,y)=\beta_0\in k\). If \(\beta_0=0\), then \(B_0(0,y)\)
is affine with nonzero slope and \(A_2(0,y)=0,\ A_1(0,y)\)
is a nonzero constant. If \(\beta_0\ne0\), the plane pair has
\(A=H(B)+\ell y+r_0\) for a quadratic \(H\) and
\(\deg_y B_0(0,y)=3\). These degrees follow from
\(\deg_y A_0(0,y)=6\) and \(\deg_y B_0(0,y)\le4\).
In this case \(A_2(0,y)=A_2(0,0)\ne0\), and \(A_1(0,y)\)
has degree at most three. Consequently in both cases
\([y^4]B_0=[y^3]B_1=0,\ [y^4]A_1=3x\).
Since \(c\) is univariate and the two coefficient ratios
\(u,v\) are dependent, \(u\) is univariate. Its identity gives
\[
 [y^2] B_1=2\sigma x,\qquad [y^2] A_2=3x^2,\qquad
 [y^j]A_2=0\ (j\ge3),\qquad x^2\mid [y]A_2 .
\]
Indeed the \(y^2\)-coefficient of \(B_1\) has \(x\)-degree
at most one and vanishes at zero. Also \(h=h_2(x)y^2+
h_1(x)y+h_0(x)\), with \(x^2\mid h_1,h_2\) and
\(\deg h_2\le2\). The coefficient of \(y^3\) vanishes
because it has degree at most one and vanishes twice at zero.

On the generic fibre \(C=t\) substitute \(z=(t-h)/a\).
The resulting planar bracket is \(\pm\kappa/a(x)\).
First suppose that the multiplicity at zero is two. The
other root of \(a\) is simple, so \(h_2=0\) because it
vanishes at that root. Give \(x,y\) weights \((-2,1)\),
and put \(\Theta=(t-h_0(0))/(a/x^2)(0)\), \(Y=x^{1/2}y\)
in a quadratic extension of the rational function field
when writing the initial forms. Then \(z_*=\Theta x^{-2}\).
The top of \(B\) on this fibre is
\[
        x^{-2}\big(\sigma(\Theta^2+2\Theta Y^2)+\beta_0\Theta\big).
\]
If \(\beta_0\ne0\), the top of \(A\) is
\(A_2(0,0)\Theta^2 x^{-4}\). Their bracket has a nonzero
term of weight thirteen, whereas the right side has weight
four. If \(\beta_0=0\), the top of \(A\) is
\[
        x^{-3}\big((\Theta+Y^2)^3+d\Theta^2\big),
               \qquad d=[x]A_2(x,0).
\]
All these expressions follow by substitution of \(z_*\),
using the displayed coefficient identities and the restrictions
of \(A_1\) on \(x=0\). For these last two forms the bracket
would have weight eleven if nonzero, so it must vanish.
Writing them \(x^{-3}U(Y),x^{-2}V(Y)\), vanishing requires
\(U^2/V^3\) constant in \(Y\). This is impossible since
\(U\ne0\) and \(V=\sigma(\Theta^2+2\Theta Y^2)\) is
nonconstant and squarefree over the generic field.

If the multiplicity is three, \(a=\delta x^3\), and write
\(h_2=\rho x^2,\ r=\rho/\delta,\ \Theta=(t-h_0(0))/\delta\).
Use weights \((-2,2)\) with \(Y=xy\). Now
\(z_*=x^{-3}(\Theta-rY^2)\), and the top forms on the fibre
are \(x^{-6}U(Y),x^{-4}V(Y)\), where
\[
\begin{aligned}
 U&=(\Theta+(1-r)Y^2)^3
       + A_2(0,0)(\Theta-rY^2)^2
       + [y^3]A_1(0,y)\,Y^3(\Theta-rY^2),\\
 V&=\sigma(\Theta-rY^2)(\Theta+(2-r)Y^2).
\end{aligned}
\]
The same coefficient bounds give these expressions.
Here \(U\ne0\) by its constant term in \(Y\), since
\(\Theta\) is transcendental over \(k\). Also \(V\) is
nonconstant and squarefree: each nonconstant factor has simple
roots and the two factors have no common root. The expected
bracket weight is twenty versus actual weight six. Vanishing
again forces \(U^2/V^3\) constant in \(Y\), which contradicts
the multiplicity at any root of \(V\).
This excludes the repeated-root obstruction and proves inversion.
The inverse descends to the original field.

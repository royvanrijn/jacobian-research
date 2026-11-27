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
The lift is either \(y^2+\delta z\) at weight
thirty or \(y^2+\delta xz\) at weight twenty,
where \(\delta\ne0\). In all these applications the
expected top determinant weight is positive: the
component weights sum to \(220\), and the moving
weight at a first tie is at most \(90\).

In the former case take \(Z=y^2+\delta z\) as a new
variable of weight thirty. The leading tuple is
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

# Two quadratic mates of the squarefree quartic scalar

**Theorem SexticKellerSquarefreeQuarticTwoQuadraticMatesExclusion.**
In characteristic zero, the squarefree quartic normal form
\(R=x(xy-1)z+b(x,y)\) has no polynomial Keller completion of
ordinary degree at most six with both mate \(z\)-degrees equal to
two. For linearly independent leading \(z\)-coefficients the
exclusion holds for arbitrary \(b\in k[x,y]\).

Work over \(\mathbb C\) via a finitely generated coefficient
field. Put \(v=xy-1,\ a=xv\). Write \(p,q\) for the leading
\(z\)-coefficients of \(P,Q\); both have ordinary degree at
most four. We first prove the assertion when \(p,q\) are
linearly independent. For dependent leading coefficients the
conclusion in the squarefree quartic **normal form** follows
from the [quadratic-affine and planar exclusions](SEXTIC_KELLER_SQUAREFREE_QUADRATIC_AFFINE_MATES.md).

## Independent leading coefficients

At fixed \(r=R\), the leading \(r\)-coefficients of \(P,Q\)
are \(p/a^2,q/a^2\), and their bracket vanishes. By
[Petravchuk–Iena, Lemma 2 and Theorem 1](sources/petravchuk-iena-math-0701588v2.pdf),
they are rational functions of a single closed pencil parameter
\(t\). Write \(t=B_0/D_0\) with coprime polynomials. Every pole
of an outer rational function lies at a parameter whose affine
fibre is supported on \(xv=0\) or is empty. In fact a nonconstant
member of the coprime pencil defines a curve and at the generic
point of each component the parameter specializes to the
corresponding value, with positive multiplicity.

At least one of the divisors \(x=0\) and \(v=0\) is a pole
of \(p/a^2\) or \(q/a^2\), since otherwise both numerators
are divisible by \(a^2\). If only \(x=0\) is a pole, then
\(v^2\) divides both \(p,q\), which contradicts independence
since \(\deg v^2=4\). If only \(v=0\) is a pole, then both
\(p,q\) are divisible by \(x^2\).

First suppose no nonzero pencil member is constant, so no
parameter value has an empty affine fibre. If the two pole
divisors lie at different parameter values, choose these
values to be zero and infinity, respectively. Then
\(t=x^e/v^f\) after rescaling, with \(e,f\in\{1,2\}\) by the
pole orders. The outer functions are Laurent polynomials
in \(t\). A generic member of the leading span has both
poles, so its Laurent expression includes a negative
exponent and a positive exponent. The terms \(a^2t^i\)
have degrees \(6+(e-2f)i\). If the degrees tie, the
leading monomials have distinct \(y\)-exponents \(2-fi\).
Consequently the two opposite exponents cannot both give
degree at most four.

In the other cases put the unique pole value at infinity.
Write \(t=B_0/D_0\) with \(D_0\) a product of positive
powers of the actual pole divisors. The outer functions are
polynomials, of maximum degree \(n\). Each pole order is at
most two, so \(n\le2\). For \(n=1\), the independent leading
span contains \(a^2\), a contradiction. For \(n=2\), the
denominator is squarefree. If both divisors are poles,
\(D_0=a\) after rescaling. The numerator \(B_0\) can be chosen
of degree at most two: degree greater than three gives a
unique excessive degree from \(B_0^2\), and at degree three
the leading forms of \(B_0,a\) must be proportional or no
nonzero binary quadratic in them can have degree below six.
In the proportional case subtract a multiple of \(a\).
The span of allowed polynomials among \(B_0^2,aB_0,a^2\)
is then one-dimensional if \(\deg B_0=2\). If \(\deg B_0=1\),
two independent members would allow one with leading
\(z\)-coefficient \(B_0^2\), excluded by the
[affine-power obstruction](SEXTIC_KELLER_SQUAREFREE_MATE_DEGREES.md#an-affine-power-as-leading-coefficient).

If only \(v=0\) is a pole, \(D_0=v\), and \(p,q\) are in
the span of \(x^2B_0^2,x^2vB_0,x^2v^2\).
The same degree comparison allows \(\deg B_0\le1\) after
subtracting a multiple of \(v\): for degree two the leading
forms of \(B_0,v\) must be proportional, and degree greater
than two makes the square term strictly highest. Thus only
\(x^2B_0^2\) can have degree at most four. This again
contradicts independence.

It remains to consider a pencil with a nonzero constant
member. Choose the closed parameter \(t\) to be a polynomial,
with this member as its denominator. The distinct finite
levels for two pole divisors would make powers of \(x\) and
\(v\) differ by a constant after rescaling, which is impossible
by their leading forms. If both divisors are poles at one
finite level, shift it to zero and write \(t=x^ev^f\) after
scaling. The only possible powers \(a^2t^i\) of degree at
most four have negative exponents, and only \(i\ge
-\min(\lfloor2/e\rfloor,\lfloor2/f\rfloor)\) are allowed.
The degrees are distinct. Independence thus requires
\(e=f=1\), in which case the span contains a constant
leading coefficient, already
[excluded](SEXTIC_KELLER_SQUAREFREE_MATE_DEGREES.md#constant-leading-coefficient).

Finally, if only \(v=0\) is a pole, shift its parameter
value to zero. Closedness of the polynomial generator
forces \(t\) proportional to \(v\), since a higher power
would be a composition. The terms \(a^2 v^i\) available at
degree at most four have \(i=-2,-1\), giving
\[
                        \langle p,q\rangle=\langle x^2v,x^2\rangle.
 \tag{1}
\]

## Excluding the last leading span

By a constant target change assume \(p=x^2 v,\ q=x^2\).
Write \(q_1\) for the \(z\)-coefficient of \(Q\) and \(d\)
for its constant term in \(z\). On \(x\ne0\), use coordinates
\(t=v,x,W=xz+q_1/(2x)\). The coordinate change has determinant
a constant sign times \(x^2\). Thus, for a constant \(\kappa\ne0\),
the determinant in the order \(t,x,W\) of
\[
               R=tW+B,\qquad Q=W^2+D,\qquad
               P=tW^2+EW+H
\]
is \(\kappa/x^2\). All of \(B,D,E,H\) belong to
\(\mathbb C[t,x,x^{-1}]\), since \(y=(t+1)/x\).
In the original coordinates \(B=b-vq_1/(2x)\) and
\(D=d-q_1^2/(4x^2)\).

The coefficients of \(W^3,W^2,W\), in that order, give
\[
 \begin{aligned}
 E&=B+f(t),\\
 H&=\tfrac{t}{2}D+f'B+h(t),\\
 (B+g)D+f''B^2+2h'B&=j(t),\qquad g=f-tf'.
 \end{aligned}
 \tag{2}
\]
For the last line, the \(W\)-coefficient is the \(x\)-derivative
of its left side at fixed \(t\). The integration functions are
polynomials since the constants of \(\partial_x\) in
\(\mathbb C[t,x,x^{-1}]\) are \(\mathbb C[t]\).

If \(B+g=0\), then \(B\) is a function of \(t\). In this
case the constant coefficient of the determinant becomes
\[
            -t D_x(D/2+f''B+h')=\kappa/x^2 .
\]
Integrating in \(x\) gives a quadratic polynomial in \(D\)
equal to a rational function of degree one in \(x\) over
\(\mathbb C(t)\). Degrees of nonconstant rational functions
multiply on composition, so this is impossible.

Otherwise (2) makes \(D\) a rational function of \(t,B\).
If \(B_x=0\), all coefficients in the three transformed
components are independent of \(x\), contradicting the
determinant. Put \(K=B+g\); now \(K\) is a nonconstant
rational function of \(x\) over \(\mathbb C(t)\). Write
\[
                  D=d_1 K+d_0+d_{-1}/K,\qquad d_1=-f'' .
\]
Polynomial division of (2) shows \(d_1,d_0,d_{-1}\in\mathbb C[t]\).
In coordinates \(t,B,W\), the constant determinant coefficient
is
\[
                -K D_t+\tfrac{t}{2}K D_B^2 ,
 \tag{3}
\]
where \(D_t\) holds \(B\) fixed. Indeed the coefficient before
elimination is \((t f'-E)D_t-t D_B(D/2+f''B+h')\).
The \(B\)-derivative of the last line of (2) gives
\(D/2+f''B+h'=-K D_B/2\), proving (3).
Multiplying (3) by \(B_x\) yields \(\kappa/x^2\).

Here \(g'=t d_1\). Substitution in (3) cancels its \(K^{-1}\)
coefficient, which is \(d_{-1}(g'-t d_1)\). The negative-power
part consists only of \(t d_{-1}^2/(2K^3)\). Thus (3) is the
\(B\)-derivative of a polynomial in \(K\) plus
\(-t d_{-1}^2/(4K^2)\). That rational primitive evaluated at
\(K(t,x)\) equals \(-\kappa/x\) plus a function of \(t\).
Its degree in \(K\) must be one by multiplicativity of
rational degree. The double pole forces \(d_{-1}=0\).
Then the primitive is affine in \(K\), with nonzero slope
\(N(t)\in\mathbb C[t]\). Hence
\[
                  B=L(t)/x+J(t),\qquad
                  L,J\in\mathbb C[t],\qquad N L=-\kappa .
\]
The coefficients \(L,J\) are polynomial since \(B\) is
Laurent polynomial in \(x\) over \(\mathbb C[t]\). In
particular \(L\) is a nonzero constant.

Because \(D\) is affine in \(B\) with coefficients in
\(\mathbb C[t]\), and \(B=b-vq_1/(2x)\), \(D\) has at most a
simple pole along \(x=0\) in the original affine plane.
Comparing with \(D=d-q_1^2/(4x^2)\) shows that \(x\) divides
\(q_1\). Now \(B\) is polynomial in \(x,y\). Multiplying
\(B=L(v)/x+J(v)\) by \(x\) and restricting to \(x=0\) gives
\(L(-1)=0\), contradicting the nonzero constant \(L\).

The exclusion for independent leading coefficients did not
use the special form of \(b\). For the squarefree quartic
normal form, proportional leading coefficients are also
excluded: cancel one by a constant target change, then apply
the planar exclusion or the quadratic-affine theorem.

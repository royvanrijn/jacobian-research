# Conic directions in maps affine in a source variable

**Theorem SexticKellerConicDirectionExclusion.**
Let \(F=zA(x,y)+B(x,y)\) be a polynomial map over a
characteristic-zero field, with \(\deg A\le5\).
If the projective image of \(A\) is dense in a conic over
the algebraic closure, the Jacobian determinant of \(F\)
cannot be a nonzero constant.

**Theorem KellerQuadraticConicOddCoefficientExclusion.**
Over a characteristic-zero field let
\[
          F=z^2 A(x)+z B(x,y)+C(x,y),
          \qquad A(x)=(x^2,x,1).
\]
If \(B_2-xB_3\) has odd degree at least three in \(y\)
over \(k(x)\), then \(F\) cannot be Keller.

## The degree bound on the conic parameter

Work over the algebraic closure. The conic is irreducible,
since it is the closure of an irreducible image. Parametrize
it by \(\nu(P,Q)=(P^2,PQ,Q^2)\), where \(P,Q\) are
coprime binary polynomials giving its rational parameter
on the source plane. The three evaluated polynomials
are coprime. After a linear target change,
\(A=G\nu(P,Q)\) for a polynomial \(G\).
Since \(A=F_z\) cannot vanish, \(G\) is a nonzero scalar
and the pair \(P,Q\) has no simultaneous affine zero.
Absorb \(G\) into the target normalization. The bound
on \(A\) gives \(\max(\deg P,\deg Q)\le2\).

If both parts of degree two are independent and coprime
as binary forms, their projective homogenizations have
no common point at infinity. Bézout gives a common affine
point, because the homogenizations have no common
component. Thus independent degree-two parts share
a linear factor. Linear changes of source and conic
parameters make them \(x^2,xy\). Write
\[
       P=x^2+\alpha x+\beta y+\gamma,\qquad
       Q=xy+\delta x+\epsilon y+\eta.
\]
If \(\beta\ne0\), substituting \(y\) from \(P=0\)
into \(Q=0\) gives a polynomial of degree three in \(x\)
with nonzero leading coefficient, hence a common zero.
For \(\beta=0\), every root of \(P\) must equal
\(-\epsilon\), because at any other root \(Q\) is
nonconstant in \(y\). Translating \(x,y\) and rescaling
gives \(P=x^2,\ Q=xy+1\), where the constant in \(Q\)
was nonzero by lack of a common zero.

In the remaining cases a linear change of conic parameters
makes \(Q\) a nonconstant affine polynomial or a nonzero
constant. In the affine case set \(Q=x\). The restriction
of \(P\) to \(x=0\) has no zero, so
\[
                     P=c+x(\alpha x+\beta y+\gamma),
                     \qquad c\ne0 .
\]
Unless \(\beta\ne0\), this pair depends on a single binary
coordinate. For \(Q\) constant normalize \(Q=1\). A binary
quadratic \(P\) has, up to affine changes, one of the
forms \(xy+c,\ x^2+\beta y+c\) with \(\beta\ne0\), or
depends on a single coordinate. This follows by
diagonalizing its quadratic part and translating in its
nondegenerate directions. These alternatives also cover
the affine-linear pairs.

## A square identity on the parameter fibres

Use \(t=P/Q\) and put \(a(t)=(t^2,t,1)\),
\(n(t)=(1,-2t,t^2)=a'(t)\times a(t)\), so \(A=Q^2a(t)\).
On each rational parameter chart \((x,t)\) in the
two-coordinate cases above, use \(\partial_x\) at fixed
\(t\). The coefficient of \(z\) in the determinant
equation of \(F\) gives \(n(t)\cdot B_x=0\) in this chart:
the only contribution is from the normal contraction of
\(\det(B_x,A_t,A)\), since \(A_x\) is proportional to \(A\).
Thus
\[
                    h(t)=B_1-2tB_2+t^2B_3
\]
is a rational function of \(t\). Write \(S=B_2-tB_3\).
The tangent vector \(B_x\) has \(a'(t)\)-coefficient
\(\partial_x S\) in the basis \(a(t),a'(t)\), and
\(n(t)\cdot B_t=h'(t)+2S\). If the original constant
determinant is \(\kappa\ne0\), the constant term in the
chart gives
\[
      \partial_x\big(S+h'(t)/2\big)^2
                  =-\kappa\,\frac{\partial_t y}{Q^2}.
      \tag{1}
\]
Here \(\det(B_x,B_t,A)=\kappa\,\partial_t y\) and
\(\det(a'(t),B_t,a(t))=-n(t)\cdot B_t\).

For \(Q=x,\ \beta\ne0\), one has \(\partial_t y=1/\beta\).
Equation (1) requires a square with a simple pole at \(x=0\)
over \(k(t)\), which is impossible. For \(Q=1,\ P=xy+c\),
the right hand side is \(-\kappa/x\), which has no
rational primitive in \(x\). For \(Q=1,\ P=x^2+\beta y+c\),
the right hand side is a nonzero constant, whose
primitives are linear in \(x\) and have odd pole order
at infinity.

For \(P=x^2,\ Q=xy+1\), one has \(y=x/t-1/x\) and
\(Q=x^2/t\), so (1) gives
\[
                   (S+h'(t)/2)^2=c(t)-\kappa/(2x^2).
\]
If \(c(t)\ne0\), multiplying by \(x^2\) gives a
quadratic with two simple roots over the algebraic
closure of \(k(t)\), again impossible for a square.
If \(c(t)=0\), the square root is a nonzero scalar
multiple of \(1/x\). In the original plane \(S\) is
regular generically on \(x=0\), because \(Q=1\) there.
The function \(t\) has order two on that line, so a
pole of \(h'(t)\) there has even order. Their sum
cannot have a simple pole, excluding this case.

Finally if \(P,Q\) depend on \(x\) alone, \(t'(x)\ne0\).
The coefficient of \(z\) in the determinant gives
\(n(t)\cdot B_y=0\). Then \(h=B_1-2tB_2+t^2B_3\)
depends only on \(x\), and \(S=B_2-tB_3\) satisfies
\[
      \partial_y\big(t'(x)S^2+h'(x)S\big)
                       =\kappa/Q(x)^2 .
\]
This follows by taking the \(a'(t)\)-component of
\(B_y\) and the \(n(t)\)-contraction of \(B_x\).
Completing the square gives a square of a rational
function of \(y\) equal to a polynomial linear in
\(y\) over \(k(x)\) with nonzero slope. This is
impossible by the pole order at infinity.

## A moving conic frame for quadratic dependence

Use the polynomial target frame \((A,A',A'')\) at each
\(x\); its determinant is a nonzero constant. Primes
on functions of \(x\) denote \(x\)-derivatives.
Write \(B=a A+b A'+c A''\), so \(b=B_2-xB_3\).
The \(z^4\) determinant coefficient is, up to the
constant frame determinant, \(2c_y\). Thus \(c\) depends
on \(x\) alone. The triangular source change
\(Z=z+a/2\) removes \(a\) from the linear term
and keeps \(b,c\) there. Use \(z\) for the new
coordinate and write \(C=U A+V A'+W A''\).

In this frame the three Jacobian columns are
\[
 (U_x,\ z^2+b_xz+U+V_x,\ (b+c')z+V+W_x),\quad
 (U_y,\ b_yz+V_y,\ W_y),\quad
 (2z,\ b,\ c).
\]
The \(z^3\) equation gives \(W_y=(b+c')b_y\).
Hence \(W=b^2/2+c'b+w(x)\).
Put \(s=b+c'\) and \(r=V+c''b+w'\).
The \(z^2\) equation, after substituting \(W\), is
\[
                 c U_y+2s V_y+2r b_y=0 .
\]
Its integral is
\[
               cU+2sV+c''b^2+2w'b=u_0(x).
               \tag{2}
\]
The \(z\) equation from the displayed columns is
\[
 (c U_x+2s(U+V_x))b_y+(sb-cb_x)U_y
                           -2(sb_x+r)V_y=0 .
\]
Use (2) and its \(x\)-derivative to obtain
\[
 sb U_y-2r V_y+
           ((2b+c')U-2c''V-c'''b^2-2w''b+u_0')b_y=0.
\]
This integrates to
\[
 sbU-r^2-\tfrac13 c'''b^3+
       ((c'')^2-w'')b^2+(u_0'+2c''w')b=v_0(x).
       \tag{3}
\]
All integrations here assert that a polynomial
annihilated by \(\partial_y\) belongs to \(k[x]\).

If \(c=0\), (2) reads \(2b(V+w')=u_0(x)\),
so \(V=-w'\) and \(u_0=0\) because \(b\) is
nonconstant in \(y\). Then (3) gives
\(b^2(U-w'')=v_0(x)\), and \(U=w''\).
If \(c\ne0\), solve (2) for \(U\) over \(k(x)\).
After substituting into (3) and completing the square,
one obtains
\[
                      (V+L(b))^2=T(b),
                      \qquad L,T\in k(x)[t].
\]
Indeed the only quadratic term in \(V\) on the left of
(3) is \(-r^2\), and after eliminating \(U\) every
coefficient is polynomial in \(b\) over \(k(x)\).
Over the algebraic closure of \(k(x)\), an odd
multiplicity root \(\lambda\) of \(T\) would give
an odd multiplicity root of \(T(b(y))\).
The polynomial \(b(y)-\lambda\) has odd degree and
at least one root of odd multiplicity, and its roots
cannot arise from any other factor of \(T\).
Thus \(T\) is a square over that field, which includes
the zero case, and \(V\) is polynomial in \(b\).
Equation (2) gives the same for \(U\).

In both cases \(U,V,W\) belong to \(\overline{k(x)}[b]\).
The whole second Jacobian column is therefore divisible
by \(b_y\) in \(\overline{k(x)}[y,z]\). Since \(b_y\)
is nonconstant, the determinant cannot be a nonzero
constant.

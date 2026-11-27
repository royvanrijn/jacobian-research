# Cubic mates in the sextic triple-line shear

**Theorem SexticTripleLineShearCubicMateExclusion.** Suppose that a
three-variable Keller map of ordinary total degree at most six over
a characteristic-zero field has a quartic component \(R=g(x)+x^2(z+xy)\).
After the shear writing this scalar as \(r=x+x^2v\), neither
mate has degree three in the remaining free variable \(y\).
A constant invertible combination of the two mates can put their
ordered free degrees in the list \((5,4),(6,4),(6,5)\).

Embed the coefficient field in \(\mathbb C\). Submersion gives
\(g'(0)\ne0\). Subtract and rescale \(R\), and rescale the two
transverse coordinates, so \(g(0)=0,g'(0)=1\). Use
\(v=z+xy+(g(x)-x)/x^2\), as in the
[quadratic-mate exclusion](QUARTIC_AFFINE_SOURCE_COMPONENTS_THROUGH_DEGREE_SEVEN.md).
All ordinary degrees are measured before the shear. Its inverse is
affine in \(y\), so a monomial \(x^i y^j z^k\) has free degree at most
\(j+k\le6\) afterwards. The quadratic-mate result already excludes
either free degree at most two.

## Leading coefficients and equal degrees

Put \(B=\mathbb C[x,v]\), \(k=\mathbb C(r)\), and
\(D=-x^2\partial_x|_r\). For degrees \(m\ge n\) and leading
coefficients \(p,q\), the equations are
\[
             D(P)Q_y-P_yD(Q)=\kappa\ne0,\qquad p^n/q^m\in k .
 \tag{1}
\]
The [support bounds](QUARTIC_AFFINE_SOURCE_COMPONENTS_THROUGH_DEGREE_SEVEN.md#support-after-the-triple-line-shear)
give weight at most six for \({\rm wt}(x)=1,{\rm wt}(v)=2\).
The \(v\)-degree of a leading coefficient for free degree
\(m\) is at most \(\lfloor(6-m)/2\rfloor\). Thus it is at most
affine in \(v\) for degrees three and four, and univariate in
\(x\) for degrees five and six.

A zero or pole of the rational fraction in (1) at \(r=c\)
forces the corresponding coefficient to contain the full
fibre \(r-c\), by squarefreeness and the substituted Bézout
identity. The weight and \(v\)-degree bounds allow at most
one full fibre in either coefficient, with residual factor
univariate of degree at most two. In equal degrees \(p/q\)
itself belongs to \(k\) by relative algebraic closedness.
If the ratio is nonconstant then, by stripping full factors
and taking a constant invertible target combination, the
leading coefficients can be made \(\ell,r\ell\) with
\(\deg\ell\le2\) and \(\ell\) univariate. This case is only
possible for equal degrees three or four. A constant \(\ell\)
is forbidden by the formal residue identity.

We use two consequences of that
[formal residue identity](QUINTIC_KELLER_BINARY_TRIPLE_LINE_SCALAR_EXCLUSION.md#free-variable-degrees),
from [Wang and Yan, Proposition 3.17](https://arxiv.org/pdf/2610.36495v1).
First, an invariant leading coefficient \(a(r)\in\mathbb C[r]\)
must vanish at zero. Multiplying the formal residue primitive
of \(a(r)^{-1/m}\) by that root's reciprocal gives
\(G\in B[a^{-1}]\) with \(DG=c\ne0\), because the coefficientwise
binomial expansion lies in \(B[a^{-1}]\) after the multiplication.
Thus \(G=c/x+H(r)\). If \(a(0)\ne0\), regularity on
\(1+xv=0\) makes \(H\) regular at \(r=0\), and the pole at
\(x=0\) cannot cancel.

Second, if the leading coefficient is \(\ell(x)\) with
\(0<\deg\ell<m\), then \(\ell\) is a power of \(x\) times
a scalar. Put \(w=\ell^{1/m}\). Multiply the formal residue
primitive by \(w\), take the normalized field trace to
\(k(x)\), and rescale by the nonzero constant. The result
\(u\in k(x)\) satisfies
\[
                         u'-\frac{\ell'}{m\ell}u=x^{-2}.
 \tag{2}
\]
Differentiation commutes with trace in characteristic zero.
At a nonzero finite pole of \(u\) of order \(j\), the highest
coefficient on the left is a nonzero multiple of
\(-j-b/m\), where \(b\ge0\) is the vanishing order of \(\ell\)
there. Thus there are no such poles. At zero the same
comparison requires a simple pole of \(u\).
At infinity a highest Laurent exponent \(j\ge0\) would
give the nonzero factor \(j-\deg(\ell)/m\), also impossible.
So \(u\) is a scalar multiple of \(x^{-1}\), and (2) gives
\(\ell'/\ell\) a scalar multiple of \(x^{-1}\).

For the equal-degree case with coefficients \(\ell,r\ell\),
subtract \(r\) times the first mate from the second. This
preserves (1) and polynomiality, and lowers the second free
degree to \(n'<m\). The
[affine-mate exclusion](KELLER_REPEATED_PARALLEL_SCALAR_AFFINE_MATE.md)
ensures \(n'\ge2\). Here \(\ell=\lambda x^b\) with \(b=1,2\)
by (2), so (1) requires \(m\mid b n'\). For \(m=3\) this
is impossible. For \(m=4\) only \(b=2,n'=2\) is possible,
with valuation one of the quadratic leading coefficient
at \(x=0\) over \(k\). This contradicts the
[quadratic-mate square relation](KELLER_REPEATED_PARALLEL_SCALAR_QUADRATIC_MATE.md).
Thus an equal-degree pair always has constant ratio.
Cancellation now preserves the original ordinary bound
and must leave unequal degrees at least three.
In particular, an equal pair of degree three is impossible.


## Two centered cubic calculations

Use a monic cubic mate \(Q=Z^3+cZ+d\) over \(k[x,x^{-1}]\).
Suppose the bracket in \((x,Z)\) is \(\sigma x^{-\epsilon}\)
with \(\sigma\in k^\times,\ \epsilon=0\) or \(1\).

For a monic quartic mate the
[centered quartic-cubic equations](QUINTIC_KELLER_BINARY_TRIPLE_LINE_SCALAR_EXCLUSION.md#the-quartic-and-cubic-free-variable-degrees)
give, for \(U=c+b_0,\ V=d+e_0\) and invariant constants
\(b_0,e_0,B_0,H\),
\[
       UV\in k,\qquad \frac4{27}U^3-B_0U^2-\frac23 V^2
                  =\frac{\sigma}{1+\epsilon}x^{-1-\epsilon}+H .
 \tag{3}
\]
If \(UV\ne0\) the factors are Laurent monomials with
opposite exponents. A positive exponent in either yields
an unmatched positive exponent on the left, and exponent
zero in both gives a constant. Thus \(UV=0\).
If \(V=0\), any pole of \(U\) gives a pole of order
divisible by three, impossible here. Thus \(U=0\).
For \(\epsilon=0\) this gives a square with an odd pole
order one. For \(\epsilon=1\), it forces
\[
                Q=Z^3+c_0(r)Z+d_0(r)+\eta(r)/x,
                 \qquad \eta^2=-3\sigma/4 .
 \tag{4}
\]

For a monic quintic mate \(P=Z^5+aZ^4+bZ^3+eZ^2+fZ+g\),
comparison of bracket powers six through two gives
\(a\in k\) and
\[
 \begin{aligned}
 b&=5c/3+b_1, & e&=4ac/3+5d/3+e_1,\\
 f&=5c^2/9+b_1c+4ad/3+f_1, &
 g&=2ac^2/9+10cd/9+2e_1c/3+b_1d+g_1 ,
 \end{aligned}
\]
with subscripted coefficients in \(k\).
The last two equations make the first of
\[
 \begin{aligned}
 L&=5c^3-36acd-27f_1c-45d^2-54e_1d,\\
 M&=4ac^3+15c^2d+9e_1c^2-18ad^2-27f_1d
 \end{aligned}
 \tag{5}
\]
invariant and the \(D\)-derivative of the second equal
to \(27\sigma x^{-\epsilon}\). These identities are replayed
in the [exact checker](verify_sextic_triple_line_cubic_mate.py).
If one of \(c,d\) has a pole at either Laurent end, \(L\)
forces \(3u=2v>0\) for their pole orders: the terms \(c^3,d^2\)
dominate the mixed and linear terms. The unique dominant
term of \(M\) then comes from \(c^2d\) and has pole order
\(2u+v\ge7\). Its integrated form instead has pole order at
most two at zero and is regular at infinity. Hence \(c,d\)
are invariant, again contradicting the equation for \(M\).

## Unequal degrees with a cubic mate

For \((m,n)=(4,3),(5,3)\), strip the full fibre factors
from (1). The remaining ratio is constant and unique
factorization gives powers \(\ell^m,\ell^3\) up to scalars
and the respective fibre factors. If \(\ell\) is constant,
both leading coefficients are invariant. Weight and the
invariant-leading residue restriction make each a scalar
multiple of \(r\). Divide by the invariant leading
coefficients and center \(Q\) without rescaling the free
variable. This gives the monic calculations with
\(\epsilon=0\), both impossible.
If \(\ell\) is nonconstant, both powers have \(v\)-degree
at most one, so \(\ell\) is univariate and the weight
bound gives degree one. Equation (2) for the cubic leading
root makes \(\ell\) proportional to \(x\). A full fibre
factor times \(x^m\) or \(x^3\) would exceed the weight
bound of six. Thus \(p,q\) are scalar multiples of
\(x^m,x^3\). Normalize their scalar factors and set
\(Z=xy+B_2/(3x^2)\), \(B_2=[y^2]Q\). Then \(Q\) is
monic centered cubic and the bracket in these coordinates
is \(\kappa/x\). The quintic case is excluded by (5).

For the quartic calculation, and also after polynomial
cancellations below, allow the leading coefficient of
\(P\) in \(Z\) to be \(A(r)\in k^\times\).
This coefficient is polynomial in \(r\) because the
\(y\)-leading coefficient \(x^4 A(r)\) is polynomial
and a pole at \(r=c\) would remain on the nonvertical
component of that fibre. Divide by \(A(r)\) in (3).
The coefficients \(c,d\) of \(Q\) in \(Z\) belong to
\(\mathbb C[r,x,x^{-1}]\), because its centered change
in \(y\) only inverts \(x\). In (4) this puts
\(c_0,d_0,\eta\) in \(\mathbb C[r]\) and gives
\(A(r)\eta(r)^2=-3\kappa/4\), so \(\eta\) is a
nonzero constant. At the original \(y=0\), \(Z=B_2/(3x^2)\).
If \(Z\) is regular at \(x=0\), the pole \(\eta/x\)
in (4) persists. If \(Z\) has a pole, its cube has
higher pole order than every other term. Both contradict
polynomiality of \(Q\).

For \((m,n)=(6,3)\), \(p\) is univariate and cannot
contain a full fibre factor. Removing a possible full
factor from \(q\) gives \(q=H(r)\ell(x)\), where \(p\)
is a scalar multiple of \(\ell^2\), and
\(1\le\deg\ell\le3\). If \(q\) has any valuation in
\(k(x)\) not divisible by three, (1) forces the degree
of its mate to be a multiple of three. One may cancel
each such leading term by an invariant multiple of
a power of \(Q\), since the leading ratio belongs
to \(k\). The process preserves the bracket over
\(k(x)\) and would leave a mate independent of \(y\),
whose bracket with \(Q\) is either of degree two or zero.
Thus every valuation of \(q\) is divisible by three.
It follows that \(\ell\) is a cubic power of a linear
form. The weight bound precludes a full fibre factor
in \(q\) times that cubic. Apply the formal residue
identity with the leading cube root \(x-a\) to \(Q\).
It requires a rational primitive of \(1/(x^2(x-a))\)
in \(x\), so a simple pole at \(a\ne0\) is impossible.
Rescale to \(q=x^3,p=x^6\).

Subtract \(Q^2\) from \(P\), lowering the degree of
the resulting polynomial mate to at most five.
The leading coefficient in the centered variable
\(Z=xy+B_2/(3x^2)\) is invariant by the highest bracket
equation. If the new degree is five or four, divide
by this invariant and use (5) or the quartic argument
following (4). If it is three, subtract the invariant
multiple of \(Q\). That invariant is polynomial in
\(r\), since the polynomial \(y\)-leading coefficient is
the invariant times \(x^3\). The affine-mate exclusion
applies below degree two. These cancellations do not
preserve the ordinary degree bound and the last
degree-two case requires the argument below.

## A cubic with a quadratic mate after cancellation

The remaining polynomial quadratic mate has leading
coefficient \(q_2=e(r)x^2\), by (1) with the cubic
coefficient \(x^3\) and relative algebraic closedness.
Here \(e\) is polynomial, since a pole would remain on
the nonvertical part of the corresponding fibre.

Use the
[degree-free quadratic normal form](KELLER_REPEATED_PARALLEL_SCALAR_QUADRATIC_MATE.md),
with \(Z=x(y+B/(2q_2))\) and \(Q_2=eZ^2+K\).
The degree-three odd part gives \(s=\deg F+1=2\).
The [Laurent factorization argument](QUARTIC_AFFINE_SOURCE_COMPONENTS_THROUGH_DEGREE_SEVEN.md#quadratic-mates-in-the-sheared-frame)
applied to \(J(K)-U=\lambda x^{-2}\) gives
\[
             K=\tau(r)+u(r)/x,\qquad F(T)=2c(r)(T-\tau(r)),
                            \qquad c u^2=\lambda\ne0 .
\]
The coefficient of \(y^3\) in the cubic is a nonzero
scalar times \(ce x^3\), so \(ce\) is constant and \(e\)
is proportional to \(u^2\). In particular \(u\in\mathbb C[r]\).
If \(e(0)\ne0\), the coefficients of \(K\) in
\(\mathbb C[r,e^{-1},x,x^{-1}]\) make \(\tau\) regular
and \(u\) a unit at \(r=0\). Then \(K\) has a simple
pole along \(x=0\), whereas \(K=C-B^2/(4q_2)\) can
have only even pole order there. Thus
\({\rm ord}_{r=0}u=b\ge1\).

On \(H=1+xv=0\) the cubic retains its leading
term \(x^3y^3\). Hence \(B|_{H=0}\ne0\), for otherwise
the quadratic becomes independent of \(y\) there
and its bracket with the cubic cannot be constant.
The derivation descends because \(D(H)=xH\).
It follows that \(K\) has pole order \(2b\) along
\(H=0\), all from \(\tau\). Along \(x=0\), \(u/x\)
is regular and \(\tau\) has the same pole order \(2b\).

Choose a branch of \(Z_*^2=(Q_2-\tau)/e\).
The formal odd integral evaluated there is
\(O_*=(4/3)c e Z_*^3\), with pole order \(6b\)
at either component. Its Laurent coefficients in
\(r\) are polynomials in symbolic \(Q_2\) by
binomial expansion at the nonzero scalar leading
coefficient of \(r^{2b}Z_*\).
At fixed \(Q_2\), moving \(K\) from \(\tau\) to
\(\tau+u/x\) changes the odd integral by
\[
            O-O_*=-\frac{c(u/x)^2}{2 e Z_*}
                   \left(1+O\!\left(\frac{u/x}{Q_2-\tau}\right)\right).
 \tag{6}
\]
Indeed the derivative in \(K\) is \(-F(K)/(2 e Z)\)
and the Taylor ratio in (6) has positive order at
both components. Here \(e Z_*\) is a unit along
each component, so (6) is regular on \(H=0\)
and has pole order two on \(x=0\).

Because \(Q_2|_{H=0}\) varies with \(y\), regularity of the
cubic forces the negative Laurent coefficients of its even
part \(\Psi(Q_2,r)\) plus \(O_*\) to vanish as polynomials
in symbolic \(Q_2\). At \(x=0\), choosing the opposite sign
of \(Z_*\) cannot cancel the pole of order \(6b>2\).
With the same sign, the even part plus \(O_*\) is regular
and the pole of (6) remains. This excludes the cubic mate
after all polynomial cancellations.

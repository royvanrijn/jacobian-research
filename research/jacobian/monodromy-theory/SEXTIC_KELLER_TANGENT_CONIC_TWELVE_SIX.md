# The free-degree pair (12,6) for the sextic tangent-conic scalar

**Theorem SexticTangentConicTwelveSixExclusion.** Let a
characteristic-zero three-variable Keller map of ordinary total
degree at most six have a quartic component \(R=g(x)+x^2(xz+y^2)\).
In the localized coordinates \(z=v-y^2/x\), the other two components
cannot have free \(y\)-degrees \((12,6)\) and respective leading
coefficients proportional to \(x^{-6},x^{-3}\).

Work over \(\mathbb C\) by coefficient-field embedding and subtract
\(g(0)\) from \(R\). Put \(g_1=g'(0)\ne0\) by submersion at \(x=0\).
Use the [high-mate setup](SEXTIC_KELLER_TANGENT_CONIC_HIGH_MATES.md):
\(r=g(x)+x^3v\), \(D=-x^3\partial_x|_r\),
\([P,Q]=D(P)Q_y-P_yD(Q)=\kappa\ne0\). The coefficient ring in
\(y\) is \(\mathbb C[r,x,x^{-1}]\). In the original variables put
\(L=xz+y^2=xv\). Every polynomial has the unique expansion
\[
                      \sum_{j,\epsilon} A_{j,\epsilon}(x,L)z^j y^\epsilon,
           \qquad \epsilon\in\{0,1\},
 \tag{1}
\]
because \(y^2=L-xz\). The terms have pairwise distinct free degrees
\(2j+\epsilon\), with leading coefficient
\((-1)^j x^{-j}A_{j,\epsilon}(x,xv)\) for each term.
For original total degree at most six,
\(\mathrm{wt}_{(1,2)} A_{j,\epsilon}\le6-j-\epsilon\), since
the relation is homogeneous in the original grading.

## Reduction to the nonvertical component

Scale the mates to the indicated leading coefficients and set
\(P_1=P-Q^2\). For any positive free degree \(j\) of \(P_1\),
its leading coefficient has \(x\)-valuation \(-j/2\) at fixed
\(r\), by the leading bracket equation with \(Q\). Thus \(j\) is
even and at most ten. In (1) the original \(P\) has
\([y^{11}]P\) a scalar times \(x^{-5}\) by the weight bound
on \(A_{5,1}\). Since \([y^{11}]P_1=0\), the coefficient
\([y^5]Q\) is a scalar times \(x^{-2}\). Center \(Q\)
by translating \(y\) by a scalar multiple of \(x\) at fixed
\(v\), an affine change in the original variables that preserves
\(R\) and the ordinary degree bound.

In the centered variables, \([y^9]P_1=0\). This is immediate
if its free degree is at most eight. If it has degree ten,
its leading coefficient is \(F_0(r)x^{-5}\) with
\(F_0\in\mathbb C[r]\). At fixed \(r\), extend to \(x=t^2\)
and set \(Z=y/t\). The rescaled \(P_1/F_0\) is monic of degree
ten and \(Q\) is monic of degree six, with zero coefficient
of \(Z^5\). The coefficient of \(Z^{14}\) in their bracket
is six times the derivative of \([Z^9](P_1/F_0)\).
This coefficient of \(P_1/F_0\) is therefore invariant under
\(D\). It is also anti-invariant under \(t\mapsto-t\), because
\(Z\mapsto-Z\). The constant field \(\mathbb C(r)\) is fixed
by that involution, so the coefficient vanishes.

In (1) write
\[
                  Q=-z^3+A z^2+C z+D_1 yz+E+F y .
 \tag{2}
\]
Here \(A,C,D_1,E,F\) are polynomials in \(x,L\).
The \(z^5\)-coefficient in (1) for \(P_1\) is \(F_0(r)\)
up to sign, with \(F_0=0\) if the free degree is below ten.
The \(z^5\)-coefficient of \(P\) is affine in \(x\) by
original degree, and that of \(Q^2\) is \(-2A\).
Hence \(A\) is affine in \(x\) plus an invariant polynomial
in \(r\). The \(yz^4\)-coefficient of \(P_1\) vanishes by
\([y^9]P_1=0\). The corresponding coefficient of \(P\)
is affine in \(x\) and that of \(Q^2\) is \(-2D_1\).
Thus \(D_1\) is affine in \(x\).

We have \(F_0(0)=0\). Otherwise (1) on the original plane
\(x=0\) gives \(P_1\) of \(z\)-degree five and \(Q\) of
\(z\)-degree three. That restriction of the original pair
is a plane automorphism by the Keller equation and the
[retained planar bound](LOW_DEGREE_COORDINATE_KELLER_MAPS_IN_THREE_VARIABLES.md).
It remains a plane automorphism after subtracting \(Q^2\).
Its restriction to a generic line \(y=y_0\) contradicts
[Abhyankar–Moh degree divisibility, García Barroso–Płoski Theorem 6.9](sources/garcia-barroso-ploski-1208.0913v1.pdf).

Set \(H=R/x=g_1+xT\), with
\(T=L+(g(x)-g_1x)/x^2\). On \(H=0\), the element \(x\)
is a unit and the localized restriction is \(r=0\),
with coefficient ring \(\mathbb C[x,x^{-1}]\) and the same
nonzero bracket. Now the free degree of \(P_1\) is at most
eight and, if positive, still even. Degree six cancels by
a constant multiple of \(Q\) by the leading bracket equation,
since the constant field is now \(\mathbb C\). Degree two
would have leading valuation \(-1\), excluded by the
[quadratic-mate square relation](KELLER_REPEATED_PARALLEL_SCALAR_QUADRATIC_MATE.md#parity-after-centering).
Degree zero is incompatible with the constant bracket.
The remaining pairs are \((8,6)\) and \((4,6)\).

## The pair (4,6) on \(H=0\)

Rescale the degree-four mate so its lead is \(x^{-2}\).
The [centered nonsquare sextic-quartic calculation](SEXTIC_KELLER_TANGENT_CONIC_TWELVE_FOUR.md#the-centered-quartic-on-the-nonvertical-component)
applies here after interchanging the mates. With the same
notation \(w=t^{-1}\), put \(Z=w(y+s)\), where \(s\) centers
the degree-four mate \(Q_4\), and \(S=Z^2+\alpha x\).
It gives \(Q_4=S^2+\delta Z/t+\beta\). The top
[sextic-quartic coefficient equations](SEXTIC_KELLER_TRIPLE_LINE_SEXTIC_QUARTIC_MATES.md#the-nonsquare-case)
then reconstruct the original sextic as
\[
 Q=S^3+\frac{3\delta}{2t}ZS+\frac{3\delta^2}{8x}
                         +b_0 Q_4+i_0 ,
\]
with all displayed parameters constant. Its \(y^5\)-coefficient
is \(6s x^{-3}\), so \(s=0\) by the original centering of \(Q\).
Replace \(x^{-1}\) by \(-T/g_1\) modulo \(H\).
The term \(x^{-3}y^6\) in \(S^3\) gives a degree-twelve
homogeneous part nonzero modulo \(x\), while every other term
has a polynomial representative of degree at most eight.
Every multiple of \(H\) has its highest homogeneous part
divisible by \(x\). As in the twelve-four calculation,
\((H)\) is the restriction kernel since \(x\) is a unit modulo
\(H\). The degree-twelve term contradicts the original
ordinary degree bound on \(Q\).

## The pair (8,6) on \(H=0\)

Set \(x=t^2\), \(Z=y/t\), and normalize the degree-eight mate
to be monic in \(Z\) using the leading bracket equation.
Write \(S=Z^2+a\). In (2) use \(z=L/x-Z^2\) and \(L=-g(x)/x^2\)
on \(H=0\). This gives
\[
                  Q=S^3+(bZ+c)S+fZ+e ,
 \tag{3}
\]
where
\[
 a=g(x)/x^3+A/3,\qquad b=-tD_1,\qquad c=-C-A^2/3,\qquad
 f=t(F+AD_1/3),\qquad e=E+2A^3/27+AC/3 .
\]
All coefficients are in \(\mathbb C[t,t^{-1}]\). The
coefficients \(b,f\) are odd in \(t\), and \(a,c,e\) are even.
Here \(a\) has a pole of order four at \(t=0\), and \(b\)
is a combination of \(t,t^3\). Also
\[
        \mathrm{ord}_0(c)\ge -4,\qquad
        \mathrm{ord}_0(f)\ge -3,\qquad
        \mathrm{ord}_0(e)\ge -6,\qquad
        \deg_t c\le10 .
 \tag{4}
\]
Indeed \(A\) becomes affine in \(x\) on this component.
The polynomials \(C,F\) have weight at most five, and \(E\)
has weight at most six. On substituting \(L=-g(x)/x^2\),
the exponents of \(x\) in \(L\) range from \(-1\) to at most
two, proving (4). A zero Laurent polynomial has order
\(+\infty\) and highest exponent \(-\infty\).

Use the unique expansion at infinity with \(W=Q^{1/6}=Z+O(Z^{-1})\).
In \(W\) the bracket requires
\[
                D(P_1)|_W=\frac{\kappa_2}{6w W^5 W_Z},
           \qquad w=t^{-1},\quad \kappa_2\ne0 ,
 \tag{5}
\]
with the normalized mate in place of \(P_1\).
Thus the coefficients down through \(W^{-4}\) are constant
and the coefficient of \(W^{-5}\) has derivative
\(\kappa_2/(6w)\). The involution sends \(W\) to \(-W\), so
the odd constant coefficients vanish. Subtracting a constant
multiple of \(Q\) and a constant removes \(W^6,W^0\).
Translating \(Q\) by a constant removes the \(W^2\)-coefficient:
the first correction to \(Q^{4/3}\) under translation is
a constant multiple of \(Q^{1/3}\). This changes only \(e\)
in (3). The resulting mate is the polynomial part in \(Z\)
of \(Q^{4/3}+uQ^{2/3}\) for a constant \(u\).

Let \(G\) be the negative part of this expansion in \(Z\),
and put \(K=c+3u/2\). By (5), \([W^{-1}]G=[W^{-3}]G=0\).
The binomial expansion, replayed by the
[exact checker](verify_sextic_tangent_conic_twelve_six.py),
gives
\[
  be+Kf=b^3/9,\qquad
  ab^3=3b^2f+3bcK-9fe,\qquad
  [W^{-5}]G=b\bigl[b^2(4c+3u)-18f^2\bigr]/243 .
 \tag{6}
\]
For the last coefficient one takes \([Z^{-5}]G+(b/3)[Z^{-2}]G\),
using \(W^{-2}=Z^{-2}-aZ^{-4}-(b/3)Z^{-5}+O(Z^{-6})\).
Integrate (5) for \([W^{-5}]G\) with \(D=-t^5\partial_t/2\).
Its additive constant is zero by anti-invariance. Thus for
some \(\rho\ne0\),
\[
                    b\bigl[b^2(4c+3u)-18f^2\bigr]=\rho t^{-3}.
 \tag{7}
\]

In the Laurent polynomial ring (7) forces \(b\) to be a unit.
Its only possible forms are \(\beta t,\beta t^3\),
\(\beta\ne0\). For \(b=\beta t\), (4) gives
\(\mathrm{ord}_0(b^2(4c+3u))\ge-2\); then (7) forces
\(\mathrm{ord}_0(f)=-2\), impossible for an odd Laurent
polynomial. For \(b=\beta t^3\), (7) forces
\(\mathrm{ord}_0(f)=-3\). The first equation in (6) now
gives \(\mathrm{ord}_0(K)\ge0\), since \(\mathrm{ord}_0(be)\ge-3\).
Eliminate \(e\) from the second equation in (6) by the first,
and then eliminate \(f^2\) by (7). This yields
\[
    a=2f/b+\frac{K(5c+3u/2)}{b^2}
                        -\frac{\rho K t^{-3}}{2 b^5}.
\]
The first two terms on the right have pole order at most
six and \(a\) has pole order four. The last term forces
\(\mathrm{ord}_0(K)\ge12\). Since \(\deg_t K\le10\), it
follows that \(K=0\). Now \(a=2f/b\) has pole order six,
a contradiction. This finishes the pair \((12,6)\).

# A coordinate reduction in the rank-two cubic pencil

**Proposition QuinticKellerRankTwoCubicPencilZeroGamma.** In the cubic-pencil
case of the [rank-two quartic restriction](QUINTIC_KELLER_RANK_TWO_QUARTIC_PENCIL.md#further-restrictions-in-the-cubic-pencil-case),
a map with \(\gamma=0\) is invertible. Consequently a noninvertible
map in that case must have \(\gamma\ne0,\ L\ne0\) and \(T(x,y,0)\)
proportional to \(L^3\). The conclusion holds over characteristic zero
by embedding the finitely generated coefficient field and descending
the inverse by uniqueness of the formal inverse.

The case \(L=0\) is the bounded completion in the preceding note.
Suppose \(L\ne0\), and take \(x=L\). With \(\gamma=0\), the formula
for \(P_3\) in that note gives
\[
 P=z^5+5z^3x+10zx^2+\beta z^4+z^2(ax+by)+\epsilon z^3+P_{\le2}.
 \tag{1}
\]
The Keller condition makes \(P\) critical-point-free. For any explicitly
constructed coordinate completion with inverse degree at most ten,
the [bounded coordinate-completion theorem](LOW_DEGREE_COORDINATE_KELLER_MAPS_IN_THREE_VARIABLES.md)
gives invertibility of \(F\). One remaining form will use direct
planar fibre degrees instead.

## When the \(y^2\) coefficient is nonzero

Completing the square in \(y\) reduces the critical-point problem
to a binary polynomial
\(\widetilde P=A(z)x^2+B(z)x+C(z)\).
Its leading univariate coefficients are \(A(z)=10z+\text{constant}\),
\(B(z)=5z^3+O(z^2)\), \(C(z)=z^5+O(z^4)\).
In particular the rational critical-value function
\(D=C-B^2/(4A)\) has leading term \((3/8)z^5\) at infinity.
At a point with \(A(z)\ne0\), a zero of \(D'\) gives a critical point
of \(\widetilde P\) by taking \(x=-B/(2A)\).

Write \(z_0\) for the unique root of \(A\). If \(B(z_0)\ne0\),
then \(D'\) has a double pole at \(z_0\); multiplying by \(A^2\)
gives a degree-six polynomial with nonzero value at \(z_0\).
It has a root away from \(z_0\), hence a critical point. If
\(B(z_0)=0\), then \(D\) is a polynomial of degree five. Either
\(D'\) has a root away from \(z_0\), or at \(z_0\) the remaining
equation
\(\widetilde P_z=10x^2+B'(z_0)x+C'(z_0)=0\) gives a critical point
because \(\widetilde P_x\) already vanishes for every \(x\).
Thus the \(y^2\) coefficient must be zero.

## When \(P\) is affine in \(y\)

Now write
\[
      P=A(z)x^2+B(z)x+C(z)+y(d x+g(z)),\qquad \deg g\le2.
 \tag{2}
\]
The leading coefficients \(A,B,C\) are as above. If \(d\ne0\), put
\(X=dx+g(z)\) and \(\xi=-g/d\). The part independent of \(y\)
in (2) at \(x=\xi(z)\) is \(K(z)=A\xi^2+B\xi+C\). More generally (2)
becomes \(P=X(y+M(X,z))+K(z)\), where
\[
                      M(X,z)=\frac{2A(z)\xi(z)+B(z)}{d}
                                         +\frac{A(z)}{d^2}X .
\]
In the coordinates \(X,Y=y+M,z\), a root of \(K'\) would give a
critical point at \(X=Y=0\). Hence \(K(z)=\lambda z+\mu\) with
\(\lambda\ne0\). The completion \((X,Y,P)\) has inverse given
successively by \(z=(P-XY-\mu)/\lambda\), \(x=(X-g(z))/d\),
and \(y=Y-M(X,z)\), of degree at most six.

Suppose \(d=0\). If \(g\) is a nonzero constant, use the graph
completion. If \(g=0\), the binary critical-point argument
for \(A x^2+B x+C\) applies directly. For nonconstant \(g\),
a simple root \(z_0\) must satisfy \(A(z_0)=0\), \(B(z_0)\ne0\).
Indeed otherwise one can solve \(P_x=0\) and then choose \(y\)
to solve \(P_z=0\), since \(g'(z_0)\ne0\). Since \(A\) has only
one root, a degree-two \(g\) cannot have two distinct roots.

If \(g\) is linear, translate its root to \(s=0\), so that
\(g=\lambda s\), \(A=10s\), \(B=B(0)+sB_*(s)\), with
\(\lambda B(0)\ne0\). The coordinate \(Y=\lambda y+10x^2+B_*x\)
puts (2) in the form \(P=sY+B(0)x+C(s)\). The inverse of
\((s,Y,P)\) first gives \(x=(P-sY-C(s))/B(0)\), of degree at most
five, and then \(y\) of degree at most ten.

## A repeated root of the coefficient of \(y\)

The remaining case is \(g=\lambda s^2\) after translation, with
\(s=z-z_0\), \(\lambda\ne0\). If \(A(0)\ne0\), fix a generic value
\(u\) of \(P\). The two roots of
\(A(0)x^2+B(0)x+C(0)=u\) continue to holomorphic arcs \(x_+(s),x_-(s)\)
on the fibre with \(y=0\). On the punctured disc in \(s\), the
fibre is the product of that disc with the \(x\)-plane and
has residue volume form \(dx\wedge ds/(\lambda s^2)\), up to sign.
This form has a regular primitive on the entire affine fibre by
the Keller identity, namely a constant multiple of \(Q\,dR\)
restricted to the fibre. Compare it on the punctured product with
\(x\,ds/(\lambda s^2)\). Their difference is closed, and its periods
on the two arcs above a circle in \(s\) are equal because those loops
are homologous in the product. The regular primitive has zero
period on either arc. Thus the residue of
\((x_+(s)-x_-(s))\,ds/(\lambda s^2)\) vanishes, forcing
\((x_+-x_-)'(0)=0\). This is the local residue comparison of
[Wang and Yan, Lemma 3.14](https://arxiv.org/pdf/2610.36495v1),
used here on the two quadratic roots.
The squared difference is \((B^2-4A(C-u))/A^2\).
Its \(s\) derivative at zero has coefficient \(-4A'(0)/A(0)^2\)
of \(u\), nonzero since \(A'=10\). This contradicts vanishing
for generic \(u\).

Thus \(A=10s\). Critical-freeness implies \(B(0)\ne0\):
if \(B(0)=0\), then at \(s=0\) the derivatives in \(x,y\) vanish
and the derivative in \(s\) is a quadratic polynomial in \(x\)
with nonzero leading coefficient. Choose the polynomial
\(\delta(s,u)\) of degree less than two in \(s\) satisfying
\[
       10s\,\delta^2+B(s)\delta+C(s)-u\equiv0\pmod{s^2}.
\]
It exists uniquely with coefficients polynomial in \(u\) by solving
successive coefficients using \(B(0)\ne0\). On each fibre \(P=u\)
use coordinates \(X,s\) via
\[
        x=s^2X+\delta(s,u),\qquad
        y=\frac{u-10s x^2-B(s)x-C(s)}{\lambda s^2}.
\]
These are polynomial substitutions. The inverse
\(X=(x-\delta(s,u))/s^2\) is regular on the whole fibre:
reducing its defining equation modulo \(s^2\) gives
\((x-\delta)(10s(x+\delta)+B(s))=0\), and the second factor
is a unit there. In the fibre ring \(s\) is not a zero divisor,
since the equation, linear in \(y\), is primitive over
\(\mathbb C[x,s]\). Thus these substitutions identify the fibre
with an affine plane. The determinant of the source change in
\(X,u,s\) is \(1/\lambda\), giving a planar Keller pair on each fibre.
At fixed \(u\), \(x,y\) in the substitution have degrees at most
three and five, respectively, so the remaining components have
degree at most twenty-five. The retained planar bound makes each
of these fibre maps injective. The coordinate \(P\) separates
the fibres, giving injectivity of the full complex map and hence
invertibility by the Ax–Grothendieck theorem. This completes
the proposition.

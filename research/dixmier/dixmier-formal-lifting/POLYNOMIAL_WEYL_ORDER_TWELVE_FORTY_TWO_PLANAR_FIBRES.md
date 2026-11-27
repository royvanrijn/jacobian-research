# Quadratic charts and generic-fibre poles in the order-\((12,42)\) profile

At a multiplicity-fifteen root of the leading factor, the
first Newton face has three characteristic roots with
multiplicities \(1,4,1\). The simple roots and the quadruple
root give different parts of the generic fibre. The first
resultant forces the simple roots to give quadratic charts
with regular \(P,Q\). A single quartic controls the critical
values of \(P\) on all those charts. The scalar successor
at the quadruple root gives two poles of \(Q\).

Use \(x\) for position, \(X\) for momentum, and the Poisson
bracket \(\{F,G\}=F_XG_x-F_xG_X\). Throughout let \(k\)
be algebraically closed of characteristic zero. Write
\(\operatorname{in}_{(-1,3)}\) in the order \(x,X\).

## Transfer of the high-root profile

**Proposition PolynomialWeylFortyTwoPlanarHighRootTransfer.**
A proper Weyl pair in the reduced order-\((12,42)\) row
would give a planar Keller pair \(P,Q\) with \(X\)-degrees
\(12,42\), leading coefficients \(h^2,h^7\), and the same
multiplicity-three and multiplicity-fifteen counts in \(h\).
At every multiplicity-fifteen root \(b\), put \(x_b=x-b\).
The transferred pair satisfies
\[
  \operatorname{in}_{(-1,3)}P=R_b^2,\qquad
  \operatorname{in}_{(-1,3)}Q=R_b^7,\qquad
  R_b=x_b^{-3}r_b(x_b^3 X),
\]
where \(r_b\) has degree six and multiplicities \(1,4,1\),
with a simple root at zero. The ratio \(q\) of the nonzero
simple root to the quadruple root satisfies \(2q^2-5q+5=0\),
as in the local two-root calculation. After \(t^3=X\) and
\(x_b=\eta+\omega_b/t\), with \(\omega_b^3\) the nonzero
quadruple root of \(r_b\), the weights
\(w(t)=7,\ w(\eta)=-12\) bound \(P,Q\) by \(2,7\).
The corresponding leading terms have \(\eta\)-degrees
at most \(8,28\), and include the terms of degrees
\(8,28\) from the shifted first face.

The [order-twelve first and scalar faces](POLYNOMIAL_WEYL_ORDER_TWELVE_FORTY_TWO_LEADING_FIELDS.md#the-two-high-root-newton-faces)
supply these finite support conditions for the Weyl pair.
The [two-root equation](POLYNOMIAL_WEYL_BERNSTEIN_FORTY_TWO_REAL_BOUND.md#the-two-roots-and-real-coefficients)
is a finite polynomial condition on the roots.
The momenta and the position translates use polynomial
Weyl coordinates. The ramified substitution uses one
cube root of \(X\) and a Laurent shear.

This uses the
[center reduction and characteristic-zero return](POLYNOMIAL_WEYL_FORTY_TWO_SPLIT_PRIME_F7_TRANSFER.md#algebraic-coefficients-and-split-primes)
with the data just listed. In detail fix finite support
bounds from the hypothetical Weyl pair. Introduce all
roots of \(h\), all characteristic roots at the
multiplicity-fifteen points and the finite coefficients
as variables. Polynomial equations express the factor
identities, the local highest symbols and all the
shifted weight caps. The weight caps in the ramified
Weyl ring use finite normal-order expansions, and
the leading term at the upper \(\eta\)-degree is
fixed by the shifted initial face. Invert the scalars
that must remain nonzero. This gives a finite system
over \(\mathbb Q\), so by the Nullstellensatz a
solution can be taken over a number field. Use large
split primes avoiding the finitely many exceptional
denominators and nonzero conditions.

By
[Lauritzen and Thomsen, Theorem 1.7(iii) and Propositions 1.11–1.12](https://arxiv.org/pdf/1508.02499v1),
the \(p\)-th powers of the reductions are central
and give a nonzero scalar Poisson bracket in the central
variables. Leading symbols at the \(X\)-degree and first
local weights commute with taking these powers, by
positive contraction loss. The center coordinates
remember the same support with coefficients replaced
by their \(p\)-th powers. In particular the coefficient
polynomial \(h\), its prescribed distinct roots and
their multiplicities are preserved in the transferred data.
For the ramified shear the center substitution is
\(\eta^p\mapsto\eta^p+\omega_b^p(t^p)^{-1}\).
The correction from the \((p-1)\)-st derivative of
\(X^{-1/3}\) is zero for \(p>3\), as in the
[right-cap calculation](POLYNOMIAL_WEYL_FORTY_TWO_DIRECT_F7_CERTIFICATE.md#the-transferred-data).
Translation by each \(b\) just translates the position
center variable by \(b^p\).
Positive contraction loss in the shifted weight gives the
same caps. Together with the first shifted face, the cap
forces the indicated maximal \(\eta\)-degrees, since the
maximizers have nonincreasing \(\eta\)-degree as the
penalty on \(\eta\) increases. These are now finite
commutative coefficient equations. Their solvability
for arbitrarily large primes gives a characteristic-zero
solution by the unit-ideal argument in the cited transfer.
The result applies after extending any characteristic-zero
coefficient field, using its finitely generated subfield
of definition.

## The shared planar resultant

**Proposition PlanarFortyTwoSharedHighRootResultant.**
Let \(P,Q\) be a planar pair with the leading \(X\)-degrees,
leading factors and first high-root faces of the transfer
proposition. Assume \(\{Q,P\}=\kappa\ne0\).
After subtracting a polynomial in \(P\) of degree at most
three from \(Q\), there is a single monic cubic \(A(P)\)
and a single constant \(c\) such that
\[
                         W=Q^2-P A(P)^2+2cA(P)
\]
has \(X\)-degree at most \(34\) and weight \(17\) at every
translated high root. In the local expansion with
\(z=x_b^3X\) fixed, choose \(\mathcal U=\sqrt P\) with
leading term \(x_b^{-3}r_b(z)\). Then
\[
                    Q=\mathcal U A(P)-c/\mathcal U+
                                     O(x_b^4).                     \tag{1}
\]
The \(O(x_b^4)\) statement is in \(k(z)((x_b))\).

The order comparison used for the
[first Weyl resultant](POLYNOMIAL_WEYL_ORDER_TWELVE_FORTY_TWO_RESULTANT.md#the-resultant-at-a-multiplicity-fifteen-root)
has the following direct Poisson version.
In the Laurent series field at \(X=\infty\), choose
\(L=P^{1/12}\) with leading term \(h^{1/6}X\).
After subtracting any terms already found, a remaining
term of \(X\)-degree \(j\ge-10\) in \(Q\) has leading
coefficient proportional to \(h^{j/6}\). This follows
from the leading Poisson equation, since the putative
degree \(j+11\) of its bracket with \(P\) is positive.
It can be canceled by a constant times \(L^j\).
For odd \(j\) the coefficient would have nonintegral
order at a multiplicity-fifteen root. The even powers
of \(L\) are rational-coefficient Laurent series, so
the argument inducts over the rational series until
degree \(-11\).

Subtract the nonnegative integer powers of \(P\) and
collect the positive odd half-powers in
\(F=P^{1/2}A(P)\), where \(A\) is monic cubic.
If the largest remaining term of \(Q-F\) has degree
\(j\ge0\) with \(6\nmid j\), then \(Q^2-F^2\) is a
polynomial of \(X\)-degree \(42+j\). At a high root
it has weight at least \(21+j/2\) from its leading
\(X\)-coefficient. Its bracket with \(P\) is \(2\kappa Q\),
so above weight \(17\) it must have a commuting top.
A polynomial top of weight \(w\) commuting with
\(x_b^{-6}r_b^2\) is a multiple of
\(x_b^{-w}r_b^{w/3}\), by the homogeneous Poisson equation.
The simple root of \(r_b\) forces \(w/3\) to be an
integer, and the \(X\)-degree of the top is \(2w\).
Thus \(2w>42+j\) in this case, a contradiction.
The same comparison for a tail with \(-7\le j<0\)
leaves only \(j=-6\). Denote its coefficient by
\(-c\) times that of \(P^{-1/2}\), including \(c=0\)
when it is absent. The correction \(2cA(P)\) makes
the resultant \(W\) have \(X\)-degree at most \(34\).

At every high root \(W\) has weight \(17\).
A higher weight would again force a commuting top
of degree at least \(36\), and a lower weight cannot
supply \(2\kappa Q\) in its bracket with \(P\).
Choose the indicated local branch of \(\mathcal U\).
In the local \(x_b\)-adic expansion \(Q+\mathcal U A(P)\)
has order \(-21\), since both terms have the same
leading symbol \(R_b^7\). The identity
\[
       Q^2-P A(P)^2=W-2cA(P)
\]
first gives \(Q-\mathcal U A(P)=O(x_b^3)\).
Dividing by \(Q+\mathcal U A(P)\) and using this
bound gives (1). The target change by a polynomial
in \(P\) preserves the scalar successor caps, since
its degree is at most three.

**Corollary PlanarFortyTwoDifferentialOrderAlgebra.**
For the pair in the preceding proposition, \(W\) has
\(X\)-degree \(34\) and the \(X\)-degree graded algebra
of \(k[P,Q]\) is
\[
                   k[\tau^6,\tau^{17},\tau^{21},\zeta],
       \qquad \tau=h^{1/3}X^2,\quad \zeta=B\tau^{49}X,
       \qquad B'-\frac{h'}{6h}B=1.
\]
Its graded fraction field has ambient degree \(3n_H\),
where \(n_H\) is the number of high roots.

At weight \(17\) write the top of \(W\) as
\(x_b^{-17}w(z)\).
The bracket with \(R_b^2\) yields
\(3w'r_b-17wr_b'=c_0r_b^6\) with \(c_0\ne0\).
At a root of multiplicity \(m\) in \(r_b\), comparison
of orders forces multiplicity \(5m+1\) in \(w\), since
\(3\operatorname{ord}(w)-17m\) cannot vanish.
The three roots give \(\deg w\ge33\), and degree
\(33\) would leave a degree-\(38\) term on the left
against degree \(36\) on the right.
Hence the degree is \(34\).

The \(X\)-degree and bracket computations for the
[second resultant](POLYNOMIAL_WEYL_ORDER_TWELVE_FORTY_TWO_LEADING_FIELDS.md#the-second-resultant)
now hold with Poisson brackets and commuting products.
In fact the commutators used there become
\(\{W^3,P\}=6\kappa W^2Q\),
\(\{P^5Q,P\}=\kappa P^5\), and
\(\{WP^2Q,P\}=2\kappa P^2Q^2+\kappa WP^2\).
They give degree \(99\) with leading coefficient
\(-(\kappa\delta^2/2)B h^{49/3}\), where
\(\delta h^{17/3}\) is the leading coefficient of \(W\).
The primitive \(B\) is obtained by the same
leading differential equation. It is rational after
undoing the quadratic twist, and has no finite poles
by the pole argument in the
[leading-factor residue proof](POLYNOMIAL_WEYL_LEADING_FACTOR_RESIDUE.md).
The spanning, independence and field calculation in the
[differential-order algebra theorem](POLYNOMIAL_WEYL_ORDER_TWELVE_FORTY_TWO_LEADING_FIELDS.md#the-differential-order-filtration)
use the two monic polynomials in \(Q\) and their leading
symbols. They apply to this commutative algebra as well.

## The simple characteristic roots

**Theorem PlanarFortyTwoSimpleRootQuadraticCharts.**
At either simple root \(\rho\) of \(r_b\), there is
a local coordinate \(u\) with \(u=r_b(z)\) at \(x_b=0\)
such that
\[
                         x_b^6 P=u^2+f_b(x_b).
\]
The critical value series \(f_b\) is divisible by
\(x_b^6\). In the chart \(u=x_b^3 v\), put
\(C_b(x_b)=f_b(x_b)/x_b^6\) and \(d=C_b(0)\).
Then \(P=v^2+C_b(x_b)\), \(Q\) is regular in this
chart, and its restriction is
\[
                 Q(0,v)=
                 \big[\sqrt{v^2+d}\,A(v^2+d)\big]_{\ge0},
                                                                    \tag{2}
\]
where the brackets take the polynomial part of the
Laurent expansion in \(v\) at infinity.
Every \(d\) on these charts is a simple root of the
same quartic
\[
  \mathcal R_A(d)-c,\qquad
  \mathcal R_A(d)=\frac{35}{128}d^4
                +\frac{5a_2}{16}d^3+\frac{3a_1}{8}d^2+\frac{a_0}{2}d,
    \quad A(y)=y^3+a_2y^2+a_1y+a_0 .
                                                                    \tag{3}
\]
In this chart \(C_b'(0)Q_v(0,0)=\kappa/r_b'(\rho)\).
At the nonzero simple root \(\rho\) the curve
\((v^2+d,Q(0,v))\) lies in the nonproper value set
of the planar pair over \(\mathbb C\).

For the quadratic normal form, the critical equation
for \(x_b^6 P\) in \(z\) has a unique solution near
\(\rho\), since \(r_b^2\) has a nondegenerate critical
point there. Subtract its value \(f_b(x_b)\); the
remaining expression is a square times a unit.
Taking the unit square root with the chosen sign
gives \(u\). This is a formal coordinate, convergent
over \(\mathbb C\) by the implicit function theorem
and the unit square root.

Let \(x_b=x_*\) for this calculation.
The series \(x_*^{21}Q\) is regular in \(x_*,u\),
since the weight bound holds before the coordinate
change. On the other hand (1) gives its expansion
modulo \(x_*^{25}\) over \(k((u))\):
\[
 \sum_{j=0}^{3} a_j x_*^{18-6j}(u^2+f_b)^{j+1/2}
                     -c x_*^{24}(u^2+f_b)^{-1/2},\qquad a_3=1 .       \tag{4}
\]
Here half-powers use the branch with leading \(u\).
The change from \(z\) to \(u\) can be made in
\(k((u))[[x_*]]\), since the formal inverse at \(x_*=0\)
maps denominators that are nonzero as functions of
\(z\) to nonzero Laurent series in \(u\).

If \(f_b\) had order \(n<6\), the coefficient of
\(x_*^{4n}u^{-1}\) from the \(j=3\) term in (4)
would be \((35/128)\) times the fourth power of
the leading coefficient of \(f_b\). For \(j<3\) the
first possible term with \(u^{-1}\) has order
\(18-6j+(j+1)n>4n\). The \(c\) term starts at
order \(24>4n\). Thus the nonzero coefficient
cannot cancel, violating regularity in \(u\).
This proves divisibility by \(x_*^6\).

With this divisibility every nonnegative \(u\)-power
in the displayed terms of (4) has combined order
at least \(21\) after \(u=x_*^3v\). Indeed the term
with \(f_b^\ell u^{2j+1-2\ell}\) has order at
least \(18-6j+6\ell+3(2j+1-2\ell)=21\).
Terms of order at least \(25\) in the regular
series for \(x_*^{21}Q\) cannot contribute to
order at most \(21\) in the chart. This gives
regularity of \(Q\) and (2). The coefficient of
\(x_*^{24}u^{-1}\) in (4) is exactly (3) minus \(c\),
so it vanishes.

The coordinate change satisfies
\(dX\wedge dx_* = u_z^{-1}dv\wedge dx_*\).
Indeed \(z=x_*^3 X\), and at fixed \(x_*\),
\(du=x_*^3dv\). At \(x_*=0\), \(u_z=r_b'(\rho)\).
The constant Jacobian in the \(v,x_*\) chart, evaluated
at \(v=x_*=0\), gives the claimed relation.
Differentiation with respect to \(d\) of the expression
in (2) before taking its polynomial part is
\(1/(2v)\) times differentiation in \(v\).
Hence \(2\mathcal R_A'(d)=Q_v(0,0)\ne0\).
This proves simplicity. Finally at \(\rho\ne0\),
the chart at any fixed finite \(v\) has
\(X\sim \rho x_*^{-3}\) as \(x_*\to0\), while
its images converge to the displayed curve.

## The complete pole profile over a high root

**Corollary PlanarFortyTwoHighRootFibrePoles.**
For a planar pair from the transfer proposition, take
the normalized generic \(P\)-fibre and project to \(x\).
Over each multiplicity-fifteen root of \(h\) there
are four unramified points where \(Q\) is regular.
The other points have \(x-b\) of orders \(1,7\) and
\(Q\) of pole orders \(1,7\), respectively.
On the latter two branches the limits of \((x-b)Q\)
have ratio \(-7\) in the order \(1,7\).
In particular the geometric degree of the planar pair
is at least eight times the number of high roots.

The simple-root charts give two unramified points each
at a generic level of \(P\), by \(P=v^2+C_b(x_b)\).
Their \(Q\)-values are finite.
At the quadruple root use \(t^3=X\),
\(x_b=\eta+\omega_b/t\), as in the transfer.
The lattice at weight \((7,-12)\) gives top polynomials
\[
        p=\alpha t^2\eta(A_0+s),\qquad
        q=\beta t G(s),\qquad s=t^{12}\eta^7,\qquad
        \deg G=4 .
\]
The highest coefficients and \(\alpha,\beta\) are nonzero
by the shifted first face. The bracket of the tops is
scalar by the full weight cap and the nonzero Jacobian.
In these variables it gives
\((A_0+8s)G-2s(A_0+s)G'=\text{const}\ne0\).
In particular \(A_0,G(0),G(-A_0)\) are nonzero and
\(G(0)=-7G(-A_0)\).

Put \(x_b=w^7\), \(t x_b=\omega_b+w^5 z_1\).
Then \(s=t^5(t x_b-\omega_b)^7\) tends to
\(\omega_b^5 z_1^7\).
In any monomial \(t^i\eta^j\) the exponent of \(w\)
is \(-7i+12j\); the remaining factor is a regular
unit times \(z_1^j\).
Consequently \(w^2P\) is regular with specialization
\(\alpha\omega_b z_1(A_0+\omega_b^5 z_1^7)\), and
\(w^7Q\) specializes to \(\beta\omega_b G(\omega_b^5 z_1^7)\).
The equation \(P=p_*\) at a generic level lifts
all eight simple roots of the first polynomial.
The seventh-root deck group acts by
\((w,z_1)\mapsto(\varepsilon w,\varepsilon^{-5}z_1)\).
Thus the zero root gives an unramified branch, and
the seven nonzero roots give a single branch of
index seven. On each branch the leading limit of
\(w^7Q\) is nonzero, giving the stated poles and ratio.
These branches account, with the four simple-chart
points, for all twelve sheets of the \(x\)-projection.

The degree of the \(Q\)-map on the normalized generic
\(P\)-fibre is the geometric degree of the planar pair.
It is the sum of all pole orders of \(Q\). This count
can be made after algebraic closure of the generic
\(P\)-level, adding over all components if necessary.
The exhibited poles contribute eight at each high root,
which gives the lower bound.

## Monodromy of the generic position projection

**Theorem PlanarFortyTwoGenericPositionMonodromy.**
For a planar pair from the transfer proposition, the
degree-twelve projection to \(x\) on the normalized
generic \(P\)-fibre has geometric monodromy \(S_{12}\).
Moreover \(k(P,Q,x)=k(x,X)\).

First the generic \(P\)-fibre is geometrically
irreducible. A nontrivial relative algebraic closure of
\(k(P)\) in \(k(x,X)\) would be a function field of
a rational curve, by Lüroth's theorem applied to a
nonconstant restriction of its rational parametrization
from the plane to a general line. Thus \(P=g(H)\)
for rational \(H\in k(x,X)\), \(g\in k(t)\) of degree
greater than one. Choose a pole of \(g\) as infinity
in its source. The denominator of the resulting \(H\)
must be a unit polynomial: any irreducible divisor of
that denominator gives a pole of \(P\). Similarly
\(g\) has no finite pole since every finite level of
the nonconstant polynomial \(H\) is nonempty and
its level polynomial has a divisor. Hence \(P=g(H)\)
with \(g,H\) polynomial. A zero of \(g'\) now gives
a critical point of \(P\), contradicting the constant
nonzero Jacobian. The function-field extension over
\(k(P)\) is consequently regular.

The pole profile gives a seven-cycle in the transitive
degree-twelve geometric group of the position projection.
It forces primitivity. In a nontrivial block system
both the block size and the number of blocks are at most
six. An element of order seven could then act neither
on the blocks nor on any block. By
[Jordan's prime-cycle theorem, as stated in Jones, Theorem 1.1](https://arxiv.org/pdf/1209.5169v1),
the group contains \(A_{12}\).

To find a transposition, fix a high root and examine the
critical values in \(X\) for \(x_*=x-b\) near zero.
The two simple-root charts give the critical values
\(C_\rho(x_*)\), each from a nondegenerate quadratic
point at \(v=0\). Their series are distinct. If their
initial values agree, the nonzero \(Q_v(0,0)\) agrees by
(2), so their slopes differ provided \(r_b'\) at the
two simple roots differs.
Write \(r_b(z)=\gamma z(z-a)^4(z-b_1)\), and
let \(q=b_1/a\). The ratio of these derivatives is
\(-(q-1)^4\). The first-face equation \(2q^2-5q+5=0\)
implies this ratio differs from one. Indeed \(t=q-1\)
satisfies \(2t^2-t+2=0\) and \(t^4=3/4-7t/8\);
the putative \(t^4=-1\) would give \(t=2\).

All other \(X\)-critical values have poles as \(x_*\to0\).
The first polynomial for \(P_X\), in \(z=x_*^3 X\),
is a nonzero multiple of \(r_b r_b'\).
It has one root at each simple root of \(r_b\),
seven at its quadruple root, and two other roots counted
with multiplicity. At these last roots \(r_b\ne0\)
and \(P\) has a pole.
At the quadruple root use the chart of the pole proof.
The derivative of the specialization of \(w^2P\)
is proportional to \(A_0+8\omega_b^5 z_1^7\).
It has seven distinct roots, none a zero of the
specialization of \(w^2P\). All lift to critical
points of \(P\) at fixed \(x_*\), since the coordinate
change has \(X\)-derivative in \(z_1\) nonzero for
\(w\ne0\). These values also have poles.
This accounts for all eleven roots of \(P_X\).

For either simple chart, \(C_\rho'(0)\ne0\). Thus for
\(p_*=C_\rho(0)+\epsilon\), the equation
\(C_\rho(x_*)=p_*\) has a nonconstant inverse series
\(x_*=\psi(\epsilon)\).
The series is algebraic, since it is a branch of the
elimination of \(P_X=0,\ P=p_*\).
It defines an \(x\)-value over the algebraic closure of
\(k(p_*)\). The distinctness and pole comparisons above
exclude any other critical point at the selected branch value.
Equivalently \(P(X,x)-p_*\), at that \(x\)-value,
has exactly one double root and all other roots simple.
Its leading coefficient is nonzero there, and the
quadratic chart with \(C_\rho'(\psi)\ne0\) gives a
simple branch of the position projection. The inertia
element is a transposition. Thus the group is \(S_{12}\).

Finally \(k(P,x)\subset k(P,x,Q)\subset k(x,X)\)
is an intermediate-field chain for the twelve-sheet
projection. Its group over \(k(P,x)\) is already
\(S_{12}\), since the geometric group is \(S_{12}\).
Primitivity leaves no proper intermediate field.
The pole profile shows that \(Q\notin k(P,x)\):
above a single high \(x\)-value there are both poles
and regular points of \(Q\), which is impossible for
a rational function of \(x\) on that fibre.
Hence \(k(P,Q,x)=k(x,X)\).

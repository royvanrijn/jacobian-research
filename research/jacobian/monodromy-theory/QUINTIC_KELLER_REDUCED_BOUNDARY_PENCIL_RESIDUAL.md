# Reduced quartic boundaries and the rank-one quintic pencil

Let \(F=(P,Q,R)\) be a three-variable polynomial Keller map of
ordinary degree five over \(\mathbb C\). The top rank means the
dimension of the constant span of the degree-five parts of its
three components.

## A reduced quartic cannot use a linear-source residual

**Theorem QuinticReducedQuarticHighMultiplicityExclusion.** If
a quartic target component is affine linear in a source variable and
has reduced leading curve, \(F\) is invertible. Consequently, in a
noninvertible map every reduced leading quartic has a required bad
point of plane multiplicity two: it is outside the ordinary node,
cusp and tacnode types, its cubic coefficient \(R_3\) vanishes there,
its full projective quadratic jet has rank one, and the applicable
higher-contact tests fail.

Write the component as \(R=a(x,y)z+b(x,y)\). Suppose first that
\(a\ne0\). In a noninvertible map the
[linear-source classification](KELLER_LINEAR_SOURCE_COEFFICIENT_CUBIC_CLASSIFICATION.md)
gives three possibilities. In the repeated-parallel case
\(R-\mu=\lambda x+x^2H\), so \(x^2\mid R_4\). In the irreducible
cubic case \(a=x^2y+\delta x-1\), the classification gives
\(b\bmod a=\lambda xy+\mu\), which has degree at most two.
Hence \(b_4=a_3 c_1\) for a homogeneous linear \(c_1\), possibly
zero, and \(R_4=a_3(z+c_1)\) is divisible by \(x^2\). In the case
\(a=x(xy-1)\), denote the two affine restrictions of \(b\) by
\(b_i(y)=\lambda_i y+\mu_i\). A representative for \(b\bmod a\)
of degree at most three is
\[
                  (1-xy)b_0(y)+xy\,b_1(y).
\]
It has the correct restrictions on the disjoint components
\(x=0\) and \(xy=1\). Again \(b_4=a_3 c_1\) and \(x^2\mid R_4\).
Thus reduced leading curve excludes every noninvertible case
with \(a\ne0\).

If \(a=0\), the reduced binary leading quartic of \(b\) has
relatively prime homogeneous partial derivatives. Indeed a
common projective zero of the derivatives would be a repeated
root of the binary quartic by Euler's identity. The two cubic
partial derivatives of \(b\), homogenized in three coordinates,
therefore have no common zero at infinity. Bézout gives a common
zero in the affine plane, contradicting the absence of critical
points of \(R=b\).

Now use the
[required bad-point theorem](KELLER_QUARTIC_QUADRATIC_SUSPENSIONS_AND_HIGHER_CONTACT.md#necessary-singular-point-in-a-reduced-leading-curve).
If the quadratic jet at that point vanished, vanishing of the
full linear and quadratic jets would put \(R\) in the preceding
affine-linear source class: after taking the point at
\([0:0:1:0]\), the local degree of a homogenized term from
\(x^iy^jz^k\) is \(4-k\). If the bad point had plane multiplicity
at least three, its rank-one quadratic jet would be a multiple
of the square of the homogenizing variable, giving invertibility
by the
[quadratic-source completion](KELLER_QUARTIC_QUADRATIC_SUSPENSIONS_AND_HIGHER_CONTACT.md#the-rank-one-normal-jet-of-a-quartic-component).
This proves the asserted sharpening.

**Corollary QuinticLowFactorReducedBoundary.** In a noninvertible
map a reduced leading quartic with all irreducible factors of
degree at most two is the union of two irreducible conics, with
a contact of order at least three. In top rank two the contact
has order four and is their only intersection.

A configuration of four lines has only ordinary nodes at its
double points. On a conic with two distinct lines, every double
point is a node or a tacnode, by the bound two on the intersection
multiplicity of a line with a conic. Higher multiplicity cannot
supply the required bad double point in either case. Two distinct
irreducible conics are smooth and meet with total intersection
multiplicity four. Their intersections of multiplicity one and two
are nodes and tacnodes, so one intersection has multiplicity at
least three. The
[rank-two conic pencil restriction](QUINTIC_KELLER_RANK_TWO_QUARTIC_PENCIL.md)
gives the stronger conclusion for top rank two.

## Rank-one pencil

**Theorem QuinticRankOneQuarticPencilResidual.** Suppose the top rank is
one. Choose \(P_5\ne0\) and \(\deg Q,\deg R\le4\). In the noninvertible
case put \(A=P_5,\ C=Q_4,\ D=R_4\). There is a primitive pencil of
coprime homogeneous forms \(p,q\) of degree \(m\), and a common
factor \(g\) with \(C=g C_0(p,q)\), \(D=g D_0(p,q)\). The binary
forms \(C_0,D_0\) are coprime of one degree \(e\) with
\(\deg g+me=4\). The following restrictions hold.

- If \(m=1\), the fixed factor \(g\) is nonreduced.
- If \(m=2\) and the general leading quartic in the cancelling
  pencil is reduced, it is a union of two conics with contact
  order at least three as above.
- If \(m=3\), the cancelling span is \(\langle z^4,zT\rangle\)
  in a linear frame, with \(T\) an irreducible cubic. One has
  \(P_5=z^2(\alpha T+\beta z^3)\), and \(T(x,y,0)\) is a nonzero
  cube of a binary linear form up to scalar.
- If \(m=4\), the top components can be written
  \(P_5=z^5,\ Q_4=z^4,\ R_4=C_*\), where \(C_*\) is irreducible
  of degree four. The binary quartic \(C_*(x,y,0)\) has a root
  of multiplicity at least three. The dehomogenization
  \(Q_3(x,y,1)\) has binary degree at most two. If its degree is
  two, the degree-two part is a scalar multiple of a square of
  a linear form, and \(C_*(x,y,0)\) is proportional to its square.

The [cubic component theorem](QUINTIC_KELLER_CUBIC_COMPONENT.md)
gives independence of \(C,D\). The highest determinant equation is
\(dA\wedge dC\wedge dD=0\). Taking \(u=C/D\) and \(v=A^4/D^5\),
and contracting the resulting equation with the Euler field,
gives \(du\wedge dv=0\). These are rational functions on
\(\mathbb P^2\). Choose a closed generator \(p/q\) of the relative
algebraic closure of \(\mathbb C(u)\) in its function field.
[Petravchuk–Iena, Lemma 2 and Theorem 2](sources/petravchuk-iena-math-0701588v2.pdf)
give this generator and the irreducibility of general members
of its pencil. Common homogeneous factors are removed so that
\(p,q\) are coprime of equal degree. Because \(v\) is algebraic
over \(\mathbb C(u)\), it also belongs to \(\mathbb C(p/q)\).
Writing \(u\) as a ratio of coprime binary forms gives the
factorizations in the theorem.

For an irreducible factor \(f\) of a pencil member \(M\), the
valuation comparison for \(v\), as in the
[rank-two pencil proof](QUINTIC_KELLER_RANK_TWO_QUARTIC_PENCIL.md#valuations-along-a-projective-pencil),
gives
\[
       4\nu_f(A)-5\nu_f(g)=\rho_M\nu_f(M)
\tag{1}
\]
with a constant \(\rho_M\) per member. Any factor of \(A\) or
\(g\) outside every member would have \(4\nu_f(A)=5\nu_f(g)\).
Since \(\deg g\le3\), no such nonzero orders are possible.
General linear combinations of \(C_0,D_0\) have \(e\) simple
roots, since \(C_0/D_0\) is a separable map of projective lines.
The corresponding moving plane members can all be taken
irreducible and coprime to the fixed factor.

### Degrees one and two

For \(m=1\) the moving curves are lines through one point. If \(g\)
is reduced, their generic union with \(g\) has only line and conic
factors, with no possibility of two conic factors, unless \(e=1\)
and \(g\) is an irreducible cubic. The first case contradicts
the reduced boundary corollary. In the second case a general
line of the pencil is transverse to the cubic away from the
base point and avoids the singularities there. This follows
by projecting the cubic from the base point and taking a regular
value on its normalization. At the base point a smooth cubic
has a tangent distinct from the general line. A singular cubic
at that point gives multiplicity at least three on the union.
Thus all double points of the leading quartic are ordinary
nodes or cusps, which also contradicts the required bad point.

For \(m=2\) every factor of a general cancelling quartic has
degree at most two, including the factors of \(g\), which has
degree at most two. The reduced boundary corollary applies.

### Degree three

Now \(e=1\) and \(g\) is linear. Set \(g=z\). Its member \(M\)
must be \(z^3\) up to scalar. Indeed if \(M\) also had another
irreducible factor \(f\), comparing (1) at \(z\) and \(f\) would
give
\[
 (4\nu_z(A)-5)\nu_f(M)=4\nu_f(A)\nu_z(M).
\]
This forces \(4\mid\nu_f(M)\), impossible in a cubic with a
factor \(z\) and a distinct factor. Thus we can choose the
second pencil member to be an irreducible cubic \(T\).

No other pencil member can be a cube of a line. Two such
members would make the pencil binary and all of its cubic
members reducible. Any contribution of \(A\) from a member
other than \(z^3\) therefore has degree divisible by three:
the multiplicities in that member have greatest common divisor
one, and their common proportionality factor in (1) is integral.
Since \(\deg A=5\), this gives \(\nu_z(A)=2\) or \(5\), and
the asserted form \(A=z^2(\alpha T+\beta z^3)\).

For a general member \(T+\lambda z^3\), the union with \(z=0\)
is reduced. Outside the intersections with the line, a singular
point of the cubic is an ordinary node or cusp. At a smooth
point of the cubic on the line, intersection multiplicity one
or two gives an ordinary node or tacnode on the union. An
intersection at a singular point has total multiplicity at least
three. By the required bad double point, some smooth branch
of the cubic and the line must have contact three. This forces
\(T(x,y,0)\) to be a cube, since the intersection divisor is
independent of \(\lambda\) and has degree three.

### Degree four

Here \(g\) is constant. In (1), the contributions to \(\deg A\)
from a member whose multiplicity greatest common divisor is
one or two are even. An odd contribution therefore requires
greatest common divisor four, hence a member \(z^4\) in a
linear frame. We can choose \(Q_4=z^4,\ R_4=C_*\) with \(C_*\)
a general irreducible member. On \(z=1\), let \(r=C_*(x,y,1)\).
It has degree four and generates the relatively algebraically
closed field of the pencil. Since \(A^4/z^{20}\) belongs to that
field, the polynomial \(A(x,y,1)\) belongs to \(\mathbb C(r)\),
hence to \(\mathbb C[r]\), by clearing a coprime univariate
numerator and denominator. Its degree is at most five.
This gives \(A=z(\alpha C_*+\beta z^4)\).

Suppose \(\alpha\ne0\). Put \(S=P-z(\alpha R+\beta Q)\), of
degree at most four, and use the binary bracket in \(x,y\).
The full determinant identity is
\[
 \kappa=S_z\{Q,R\}-Q_z\{S,R\}+R_z\{S,Q\}
                        +(\alpha R+\beta Q)\{Q,R\}.
\]
For \(j=3,2,1\) in descending order, assume that all \(Q_i\)
with \(i>j\) and \(S_i\) with \(i>j+1\) are pure powers of
\(z\). The degree-\(j+6\) equation then is
\[
          \{4z^3 S_{j+1}-(\alpha C_*+\beta z^4)Q_j,C_*\}=0 .
\tag{2}
\]
The \(S_z\)-term has degree at most \(j+5\), and the
\(R_z\)-term has degree at most \(2j+2<j+6\). The remaining terms
of degree \(j+6\) are precisely (2). On \(z=1\), the expression
in the first slot is polynomial in \(r\) by relative closure
and characteristic zero. If \(Q_j(x,y,1)\) were nonconstant,
its product with \(\alpha r+\beta\) would have degree between
five and seven and could not be canceled by \(S_{j+1}(x,y,1)\),
whose degree is at most four. Such a degree is impossible for
a polynomial in \(r\). Hence \(Q_j\) is a pure power of \(z\).
For \(j=3\), \(S_4\) is a constant combination of \(C_*,z^4\);
subtract the corresponding multiple of \(R\) from \(P\).
For \(j<3\), \(S_{j+1}\) has degree less than four and is
already a pure power of \(z\). The induction makes \(Q\) a
univariate quartic, impossible for a Keller component.
Therefore \(\alpha=0\), and rescaling \(P\) gives \(P_5=z^5\).

A general plane member \(C_*+\lambda z^4\) is smooth on \(z\ne0\),
because general level curves of the nonconstant polynomial \(r\)
are smooth by generic smoothness. If \(C_*(x,y,0)\) had root
multiplicities at most two, its simple roots would give smooth
points on that plane member. At each double root, loss of the
linear jet is allowed, but the quadratic coefficient along the
line at infinity is nonzero. Split off that square keeping the
normal coordinate \(z\). The remaining one-variable series has
order at most four for general \(\lambda\), since \(\lambda z^4\)
changes only its quartic coefficient. This gives an ordinary
node, cusp or tacnode. The quartic boundary exclusion forces
a root of multiplicity at least three.

For the constraint on \(Q_3\), set \(x=\tau X,\ y=\tau Y\),
\(z=\tau Z(X,Y,\tau^{-1})\), where \(Z=1+O(\tau^{-1})\) solves
\(Q=\tau^4\) formally. The expansion exists with polynomial
coefficients in \(X,Y\) by formal implicit solution at \(Z=1\),
where the derivative of \(Z^4\) is four. The binary bracket
of the substituted \(P,R\) at fixed \(\tau\) is
\(\pm\kappa\tau^2/Q_z=O(\tau^{-1})\). This follows by differentiating
\(Q=\tau^4\) in \(X,Y\) to remove \(Z_X,Z_Y\).
Use \(q_j=Q_j(X,Y,1)\), \(p_j=P_j(X,Y,1)\).
The coefficient of \(\tau^4\) in \(P\) is \(p_4-5q_3/4\).
It commutes with the leading \(r\) of the substituted \(R\), and
so is affine linear in \(r\). Subtract a multiple of \(R\) from
\(P\) to make this coefficient a constant \(c_0\). Thus
\(P_4=\frac54 zQ_3+c_0z^4\).

The coefficient of \(\tau^3\) in \(P\) is then
\[
                  p_3-\tfrac54 q_2-c_0 q_3-\tfrac5{32}q_3^2 .
\]
To compute it, write \(Z=1+\tau^{-1}Z_1+\tau^{-2}Z_2+\cdots\).
Then \(Z_1=-q_3/4\) and
\(4 Z_2=-6Z_1^2-(\partial_zQ_3)(X,Y,1)Z_1-q_2\).
Substitution and the formula for \(P_4\) give the displayed
expression. It also commutes with \(r\) because the earlier
coefficients of \(P\) are constant in \(X,Y\).
A binary degree-three term in \(q_3\) would give degree six,
which a polynomial in \(r\) cannot have. If \(q_3\) has degree
two, the only degree-four term of the expression is
\(-5/32\) times the square of its top quadratic. It must be
proportional to the top binary term of \(r\), namely
\(C_*(X,Y,0)\). The root of multiplicity at least three on the
latter then forces that quadratic to be a square of a line
up to scalar.

The results transfer to characteristic zero by embedding a
finitely generated coefficient field into \(\mathbb C\), using
geometric factorization and geometric pencil properties over an
algebraic closure. Invertibility conclusions descend by
uniqueness of the formal inverse.

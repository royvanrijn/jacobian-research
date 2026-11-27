# Coprime reduced lower forms in mixed sextic degree profiles

Let \(F=(F_6,F_d,F_e)\) be a characteristic-zero three-variable
Keller map of the indicated component degrees, and write
\((A,B,C)\) for the highest homogeneous parts of these components.
The statements below are over the algebraic closure. In the first
three theorems the pair \(B,C\) is coprime and squarefree. In the
equal-degree profiles this implies their linear independence.

**Theorem SexticKeller654ReducedLowerPencil.**
For \((d,e)=(5,4)\), either \(A,B,C\) are binary in the same two
linear variables, or after linear source changes and constant target
rescalings they have the form
\[
 C=C_4(x,y,z),\qquad B=x(C+\lambda x^4),\qquad
 A=x^2 C+\frac53\lambda x^6,\qquad \lambda\ne0.
\]
In the second case the binary quartic \(C(0,y,z)\) has a root of
multiplicity at least three, and \(C(1,y,z)\) has one of the forms
in the exact area boundary theorem below. In particular \(B\) is
reducible in both alternatives.

**Theorem SexticKeller654NonbinaryFirstCorrection.**
Use the nonbinary frame of the preceding theorem. Denote the
parts immediately below \(A,B,C\) by \(A_5,B_4,C_3\), respectively.
Then \(A_5-2x B_4+x^2C_3\) belongs to \(\langle x^5,xC\rangle\).
Restrict \(G=(B_4-xC_3)(1,y,z)\) to the level
\(C(1,y,z)=-5\lambda\). In either of the two linear-in-\(z\)
forms with two places at infinity, restrict to the open
parametrization with nonzero \(y\) in the normalized frame.
There \(G\) is constant or constant plus one Laurent monomial.
In either of the two one-place forms, the level has a polynomial
affine-line parameter and \(G\) is affine in that parameter.

**Theorem SexticKeller654IrreducibleQuintic.**
For \((d,e)=(5,4)\), suppose \(B\) is irreducible and \(C\) is
coprime to \(B\), allowing repeated factors in \(C\).
Either \(C=x^4,\ A=x^6\) after constant rescalings and a linear
source change, with closed projective parameter \(B/x^5\), or
\(C=H^2,\ A=H^3\) after constant rescalings, with \(H\) quadratic
and closed projective parameter \(B^2/H^5\).

**Theorem SexticKeller655ReducedLowerPencil.**
For \((d,e)=(5,5)\), either \(A,B,C\) are binary in the same two
linear variables, or in linear source and degree-preserving target
frames the leading tuple is
\[
                             (x^6,x^5,H_5(x,y,z)).
\]
The projective parameter \(H_5/x^5\) in the second case is closed,
meaning that its field is relatively algebraically closed in the
projective plane function field.

**Theorem SexticKeller644ReducedLowerPencil.**
For \((d,e)=(4,4)\), either \(A,B,C\) are binary in the same two
linear variables, or the pencil spanned by \(B,C\) contains the
square of a quadratic form. If the degree of its closed projective
parameter is two, one can take the underlying conic pencil to be
\(\langle xz+y^2,x^2\rangle\) or \(\langle xz+y^2,xy\rangle\).
In the second case there are degree-preserving target coordinates
with \(B=p^2,\ C=q^2\) and
\[
 A=a_0p^3+a_1p^2q+a_2pq^2+a_3q^3,\qquad a_1a_2=0,\qquad
 p=xz+y^2,\quad q=xy.
\]

**Theorem KellerSaturatedLeadingPencilExclusion.**
Let \(f_i(p,q)\) be three ordinary homogeneous leading forms of
degrees \(d_i\ge2\), where \(p,q\) are coprime homogeneous forms of
positive degrees and the polynomial kernel of \(\det J(p,q,\cdot)\)
is \(S=k[p,q]\). Give the two variables of \(f_i\) the degrees of
\(p,q\). If the \(2\times2\) minors of the matrix with columns
\((f_{i,p})_i,(f_{i,q})_i\) have greatest common divisor one in
\(S\), no three-variable Keller map has these leading forms.
In particular, in the binary alternative of any of the three profiles
above the three binary derivative minors must have a nonconstant
common factor.

**Theorem QuarticPlaneExactAreaBoundaryReduction.**
Let \(K(y,z)\) be a polynomial of degree four over an algebraically
closed characteristic-zero field whose highest binary quartic has a
root of multiplicity at least three. Suppose its generic affine
level curve is geometrically integral and its relative area
differential \(dy\wedge dz/dK\) has a rational primitive on the
geometric generic curve. Its normalization has genus zero. Up to
affine changes in the plane and an affine rescaling of \(K\), the
polynomial is one of
\[
                 y^3z+b(y),\qquad y^2z+b(y),\qquad
                 z+b(y),\qquad (z+a(y))^2+y,
\]
where \(\deg b\le4\), \(\deg a\le2\). The first two have two places
at infinity and the last two have one.

## Closed parameters and polynomial constants

Work over an algebraically closed field \(k\) of characteristic zero.
The top determinant equation gives \(\det J(A,B,C)=0\). Since
\(B,C\) are independent, \(A\) is algebraic over \(k(B,C)\).
The rational kernel of a nonzero Jacobian derivation is the
relative algebraic closure of its two independent first integrals,
by the characteristic-zero Jacobian criterion.

Let \(L=k(\mathbb P^2)\) be the rational functions of homogeneous
degree zero in \(k(x,y,z)\). For a nonconstant projective ratio \(R\)
below, the relative closure of \(k(R)\) in \(L\) is \(k(t)\) for
a closed parameter \(t=u/v\) with coprime degree-\(m\) homogeneous
forms \(u,v\). This is the generative-function assertion applied to
the rational plane function field.
[Petravchuk–Iena, Lemma 2](sources/petravchuk-iena-math-0701588v2.pdf)

Write \(R=r(t)\) where \(r\) has degree \(s\). A zero or pole of \(r\)
of order \(j\) pulls back to the corresponding member of the
\((u,v)\) pencil, taken to power \(j\). Distinct members are
coprime, so the common homogeneous degree of a reduced numerator
and denominator of \(R\) is \(ms\). A closed parameter of degree one
has two independent linear coordinates \(u,v\).

A homogeneous rational function \(h\) of degree one satisfies
\(k(x,y,z)=L(h)\), since it is a radial linear coordinate times
a nonzero element of \(L\). The relative closure of \(k(R,h)\)
in \(L(h)\) is \(k(t,h)\). Indeed \(L/k(t)\) is regular in
characteristic zero, so adjoining the indeterminate \(h\) preserves
relative closedness, and \(k(t,h)/k(R,h)\) is algebraic. For forms
of the same degree, if \(R=B/C\), every degree-zero rational
function algebraic over \(k(B,C)\) is algebraic over \(k(R)\).
To see this, write \(k(B,C)=k(R,C)\) and use that \(C\)
is transcendental over \(L\).

We will also use polynomial contraction. For an onto pair \(p,q\)
of independent polynomial functions, \(k(p,q)\cap k[x,y,z]=k[p,q]\).
A reduced quotient with a nonconstant denominator would have a
zero of its denominator off the numerator on the target plane,
contradicting the cross-multiplied polynomial identity at a
preimage of that zero.

## Descending through the lower terms

For the saturation theorem, set \(R_0=k[x,y,z]\). The sequence
\(p,q\) is regular because \(p,q\) are coprime. Lift a homogeneous
\(k\)-basis of \(R_0/(p,q)\) to \(R_0\), including \(1\).
Multiplication makes the lifts a homogeneous \(S\)-basis of \(R_0\).
Surjectivity follows by degree induction modulo \((p,q)\).
The regular sequence multiplies the Hilbert series by
\((1-T^{\deg p})(1-T^{\deg q})\), so the source and target of
that multiplication map have equal dimensions in each degree.

Use a reference tuple \(G_i(\widetilde p,\widetilde q)\) with
\(\widetilde p,\widetilde q\) having leading forms \(p,q\), and
\(G_i\) having weighted leading form \(f_i\). Start with the
pure leading tuple. If a nonzero error remains, let \(r\ge1\)
be the smallest degree loss among its three entries and let
\(E\) be its tuple of parts of degrees \(d_i-r\).
Put \(n=f_p\times f_q\), evaluated at \((p,q)\). The coefficient
of degree \(\sum_i d_i-3-r>0\) in the Keller determinant is
\[
                         \det J(p,q,n\cdot E).
\]
Indeed a smaller error or two error entries contribute to lower
degrees, and the reference tuple itself has zero determinant.
Also each entry of \(n\) is killed by the derivation.
It follows that \(n\cdot E\in S\).

For any positive-degree basis lift \(b\) of degree \(h\), its
coefficient vector \(E_b\) in the \(S\)-basis expansion satisfies
\(n\cdot E_b=0\). Over the fraction field of \(S\) write it as
\(a_b f_p+c_b f_q\). For each irreducible prime of \(S\), some
minor is a unit at that prime, so solving there gives
\(a_b,c_b\in S\) by unique factorization. By homogeneity their
degrees are \(\deg p-r-h,\deg q-r-h\), respectively.
Replace \(\widetilde p,\widetilde q\) by
\(\widetilde p+\sum_b b a_b,\widetilde q+\sum_b b c_b\)
and add the coefficient vector of the lift \(1\) to \(G\).
These changes account for \(E\) to degree loss \(r\).
Every term involving two changes has loss at least \(2r\).
Repeating increases the loss until the error is zero. The map
then factors through two polynomials and has zero determinant.
In the binary application take \(p,q\) to be the two linear
coordinates. The remaining source coordinate is free, so the
polynomial kernel is \(k[p,q]\).

## The profile \((6,5,4)\)

Take \(R=B^4/C^5\) and \(h=B/C\). Then
\[
                        B=h^5/R,\qquad C=h^4/R.
\]
Thus \(k(B,C)=k(R,h)\) and its relative closure is \(k(t,h)\).
The cases \(ms=20\) exhaust the possibilities below.

If \(s=1\), \((B,C)\) is onto. For two nonzero prescribed
values \(b,c\), choose a projective point of \(B^4-(b^4/c^5)C^5=0\)
off \(BC=0\). Such points exist because the member is coprime
to \(BC\). Then rescale to the two values, using coprimeness
of 5 and 4. For one prescribed value zero, choose a point on
the corresponding zero curve off the other and rescale.
The origin supplies the remaining point. Consequently
\(A\in k[B,C]\) by contraction. This ring has no nonzero
element homogeneous of degree six.

At a zero of \(r\) of order \(i\), squarefreeness of \(B\)
forces all factors of its member to have multiplicity \(4/i\).
Thus \(i\mid4\) and \(4\mid mi\). At a pole of order \(j\)
one has \(j\mid5\) and \(5\mid mj\) by squarefreeness of \(C\).
Use these restrictions as follows.

- For \(s=2,m=10\), the two pole members are fifth powers
  of coprime quadratics. The zero member is a square and,
  as a nontrivial combination of the pole members, splits
  into five distinct coprime quadratic members of the
  pencil of the two quadratic roots. They would all be
  squares. A pencil spanned by two independent squares has
  rank two at every other point, so this is impossible.
- For \(s=4,m=5\), the pole members are four distinct fifth
  powers of linear forms. The pencil containing two such
  powers contains no third, by the mixed terms in the
  binomial expansion of any other fifth power.
- For \(s=5,m=4\), the pole order is five, so \(C\) is
  a member of the pencil. Some zero order is one, giving
  a fourth power \(x^4\) in the pencil. The only possible
  zero partition is \(4+1\). There cannot be three distinct
  fourth powers, by the same binomial argument. The other
  partition to check is \(2+2+1\), in which a fourth power
  would be a nontrivial sum of two squares of coprime
  quadratics. Over \(k\) this sum factors into two coprime
  nonconstant quadratics, so it cannot be a fourth power
  of a linear form. This leaves \(B=x(C+\lambda x^4)\)
  after rescaling, with \(\lambda\ne0\) by distinctness
  of the three members.
- For \(s=10,m=2\), there are two poles of order five
  and the zero orders are two or four. Order two gives
  a square of a linear form, so there are at most two such
  orders. Their number is odd, hence one. Thus in a pencil
  \(\langle x^2,H\rangle\) there are four distinct members
  \(H-a_i x^2,H-b_i x^2\), \(i=1,2\), with
  \[
   B=x(H-a_1x^2)(H-a_2x^2),\qquad
   C=(H-b_1x^2)(H-b_2x^2).
  \]
- For \(s=20,m=1\), \(B,C\) are binary in the same linear
  parameters. Then \(k(t,h)\) is the rational field of
  those two coordinates, and \(A\) belongs to their
  polynomial ring. This is the binary case.

For the two intermediate cases the relative closure is
\(k(x,C)\) or \(k(x,H)\), respectively. Indeed their second
coordinate divided by the corresponding power of \(x\) is
a closed projective parameter, and \(x\) is a radial
generator. The polynomial kernel is the corresponding
polynomial ring. Each pair is onto since the highest binary
part of \(C(x,y,z)\) or \(H(x,y,z)\), for fixed \(x\), is
nonconstant and independent of the value of \(x\).

The leading pair \(B,C\) has independent differentials
generically on \(x=0\) in both cases, since \(B_x|_{x=0}=C_0\)
and \(C_0=C(0,y,z)\) is nonzero of degree four. Apply the
[rational slice with divisor regularity](KELLER_LEADING_PAIR_RATIONAL_SLICE.md)
and divide by the nonzero determinant. A primitive \(\phi\)
regular generically on that plane would give
\[
                           \{C_0,\phi|_{x=0}\}=1/C_0
\]
in the binary field, up to a nonzero constant, with
\(\{f,g\}=f_y g_z-f_z g_y\). On each generic fibre of
\(C_0\) the relative area differential therefore has a
rational primitive.

To see the restriction on \(C_0\), write \(C_0=y^4c(t)\)
with \(t=z/y\), using a line at infinity in the parameter
avoiding the roots. The relative form along \(C_0=c_0\ne0\)
is proportional to \(y^2\,dt\). Above a root of \(c(t)\)
of multiplicity \(v\), put \(g=\gcd(4,v)\). The ramification
index is \(4/g\), the pole order of \(y\) is \(v/g\),
and the differential order is \((4-2v)/g-1\). A double
root gives a simple pole. If all roots are simple the form
is everywhere regular, also at parameter infinity where
\(y\) has a simple zero. Both a simple pole and an
everywhere regular nonzero form on the smooth complete
curve contradict exactness. This leaves a root of
multiplicity three or four.

For \(m=2\), \(C_0=H(0,y,z)^2\), so the binary restriction
of \(H\) must be a square. In a linear frame
\(H=y^2+xz\). The coefficient of \(z\) after normalizing
the square is nonzero, since otherwise \(H/x^2\) would be
a degree-two function of \(y/x\) and thus not closed.
Let \(\Delta=\det\partial(B,C)/\partial(x,H)\). Its
restriction to \(x=0\) is \(2H^3\), and
\(-y/(x\Delta)\) is a primitive for \(\det J(B,C,\cdot)\).
Any other primitive differs by an element of \(k(x,H)\).
Such an element can have leading pole coefficient on \(x=0\)
only in \(k(y^2)\), since \(H=y^2\) there. Cancellation
would require a nonzero scalar multiple of \(y/(2y^6)\),
which is outside \(k(y^2)\). This excludes \(m=2\).

For \(m=4\), contraction gives \(A=ax^6+bx^2 C\).
The two relevant derivative minors in \(k[x,C]\) are
\(C+5\lambda x^4\) and \(2x(bC+3ax^4)\). If they are
coprime, the saturation theorem applies. Otherwise
\(b\ne0\) and \(3a=5\lambda b\). Rescale \(A\) to
obtain the claimed exceptional shape.

In this case put \(Y=y/x,\ Z=z/x,\ K(Y,Z)=C(1,Y,Z)\).
The closedness of \(C/x^4\) gives geometric integrality of
the generic level of \(K\). The rational primitive for the
leading pair \(B,C\) satisfies
\[
 \{K,\phi(x,Y,Z)\}_{Y,Z}=\frac{1}{x^6(K+5\lambda)}
\]
up to a nonzero scalar, since
\(\det J(B,C,\cdot)=(C+5\lambda x^4)\det J(x,C,\cdot)\)
and \(C=x^4K\). Specializing \(x\) at a nonzero value
avoiding the denominator gives a rational primitive on the
geometric generic level of \(K\), after multiplication by
the constant \(K+5\lambda\) there. This applies the exact
area boundary theorem.

## The first correction on the nonbinary discriminant level

Use the actual quartic component \(F_4\) in the reference tuple
\((x^2F_4+(5/3)\lambda x^6,\ xF_4+\lambda x^5,\ F_4)\).
This tuple has zero determinant since it factors through \(x,F_4\).
Put \(h=C+5\lambda x^4\). The two derivative columns of its
leading tuple with respect to \(x,C\) are
\((2xh,h,0)\) and \((x^2,x,1)\). The weight-loss-one
determinant equation therefore gives
\[
 \det J(x,C,h(A_5-2xB_4+x^2C_3))=0.
\]
Here \(h\) belongs to the constant ring \(k[x,C]\) and is
nonzero. The homogeneous degree-five part of that ring is
\(\langle x^5,xC\rangle\), giving the first assertion.

Now rescale the whole map in an indeterminate \(\epsilon\),
replacing its component of degree \(d\) by
\(\epsilon^d F_d(x/\epsilon,y/\epsilon,z/\epsilon)\).
The determinant is \(\kappa\epsilon^{12}\).
On the fixed quartic level \(c=-5\lambda\) use local
coordinates \(x,t\) near \(x=1\) and a point of the indicated
parametrization. For \(C(1,Y,Z)=Y^j Z+b(Y)\) one uses
\(t=Y=y/x\), with \(t\ne0\) when \(j>0\), and solves for
\(Z=z/x\). For the square form use \(t=Z+a(Y)\) and solve
for \(Y\). Affine normalization of the plane polynomial
transports the level value as well. In both types the defining
derivative for solving the rescaled quartic level is a unit
at the chosen point. Formal implicit solution gives the
restricted rescaled quintic component as
\[
         \widehat B(x,t,\epsilon)
          = xc+\lambda x^5+\epsilon \widetilde G(x,t)+O(\epsilon^2),
          \qquad \widetilde G(1,t)=G(t).
\]
The correction is \(B_4-xC_3\) because the fixed quartic
equation replaces \(xC\) at order one by \(xc-\epsilon xC_3\).

At \(x=1\), the leading derivative in \(x\) is zero and
its derivative in \(x\) is \(20\lambda\ne0\).
Thus one can solve \(\partial_x\widehat B=0\) formally
as \(x=X(t,\epsilon)\) with constant term \(1\).
The remaining derivative
\(\epsilon^{-1}\partial_t\widehat B(X(t,\epsilon),t,\epsilon)\)
has constant term \(G'(t)\) in \(\epsilon\).
If \(G'\) is nonzero and has a root on the parametrized
open set, Weierstrass preparation at that root gives
a zero of the remaining derivative over a Puiseux extension
of the Laurent series field. Both quintic derivatives
then vanish on the quartic level, contrary to the
nonzero determinant at \(\epsilon\ne0\).

On an affine line a polynomial derivative with no root
is a nonzero constant, so \(G\) is affine or constant.
On the indicated multiplicative open set a nonzero
Laurent polynomial with no root is a single Laurent
monomial. Applying this to \(G'\) gives a constant plus
one Laurent monomial for \(G\), or a constant if \(G'=0\).

## An irreducible quintic with a possibly nonreduced quartic

In the \((6,5,4)\) setup the zero fibre of \(r\) has just one
point because \(B\) is irreducible. Its order \(s\) divides four,
so \((s,m)\) is \((1,20),(2,10)\), or \((4,5)\). The case
\(s=1\) is excluded as before by contraction. At a pole of
order \(j\le4\), each factor of the corresponding member has
multiplicity divisible by five, because \(j\) is coprime to
five and the product of pole members is \(C^5\).

If \(s=2,m=10\) and there are two simple poles, their members
are coprime fifth powers of quadratics and the zero member is
the square \(B^2\). The same factorization into five quadratic
squares is impossible. Thus there is one double pole and
\(C=H^2\) with \(H\) quadratic. The ratio \(B^2/H^5\) is
a closed parameter because its square \(B^4/C^5\) has outer
degree two through this pencil. The pair \(B,H\) has relatively
closed field, by its closed projective ratio and the radial
generator \(B/H^2\). It is onto by the same projective
member and radial argument, now for coprime degrees five
and two. Hence \(A\in k[B,H]\), and its degree forces
\(A\) to be proportional to \(H^3\).

If \(s=4,m=5\), every pole member is a fifth power of a linear
form. There are at most two such members, and if there were two,
\(B\) would be a binary form in the two corresponding linear
forms, hence reducible. Thus there is one pole of order four
and \(C=x^4\). The parameter \(B/x^5\) is closed and the
pair \(x,B\) is onto. Consequently \(A=ax^6+bxB\).
In \(k[x,B]\) the derivative minors include \(4x^3\) and
\(6ax^5+bB\) up to sign. The saturation theorem excludes
\(b\ne0\), leaving a pure sixth power.

## The equal-degree lower pencils

Set \(R=B/C\). The closed parameter degree \(m\) divides
the common degree \(d\). If \(m=1\) the tuple is binary
as above, since the pair depends on two linear variables
and generates a subfield of finite index in their field.

For \((6,5,5)\) the remaining degree is \(m=5\), so
\(R\) itself is closed. The rational function
\(\eta=B^6/A^5\) has degree zero, is algebraic over
\(k(B,C)\), and therefore belongs to \(k(R)\).
Viewed on the parameter line, its order at \(R=0\)
is \(1\) modulo five by squarefreeness of \(B\).
Its order at \(R=\infty\) is zero modulo five
by squarefreeness of \(C\) and coprimality.
Since the sum of its orders is zero, there is
another value of \(R\) with order not divisible
by five. Along each geometric component of that
member, \(B\) is nonzero, so the order of \(\eta\)
is divisible by five. All component multiplicities
of that degree-five pencil member are therefore
divisible by five. The member is a pure fifth power.

Take the pencil as \(\langle x^5,H\rangle\).
The parameter \(H/x^5\) is closed, so after adjoining
the radial coordinate \(x\) the constant field of
the lower pair is \(k(x,H)\). Coprimality makes that
pair onto for the same reason as in the mixed case.
Hence \(A=ax^6+bxH\). If \(b\ne0\), its derivative
minor with \(H\) is \(6ax^5+bH\), coprime to the
minor \(5x^4\) of \((x^5,H)\). The saturation
theorem excludes this, leaving a pure sixth power.

For \((6,4,4)\) and \(m=4\) the same argument uses
\(\eta=B^3/A^2\). Its orders at zero and infinity
are odd and even respectively, and some other
member consequently has all factors of even
multiplicity. This is the required squared quadric.

For \((6,4,4)\) and \(m=2\), write \(B,C\) as binary
quadratics in the generators \(p,q\) of the closed
conic pencil. The leading pair slice gives exactness
of the volume differential on their generic fibres.
Fixing \(B,C\) generically fixes \(p,q\) up to
finitely many choices, and the corresponding volume
differential for \(p,q\) differs only by a nonzero
constant on each fibre. On the generic projective
conic of \(q-tp=0\), the affine fibre is the double
cover obtained by adjoining \(1/\sqrt{p_\ell}\),
where \(p_\ell\) is the dehomogenization in a chart.
Division of the radial volume by \(dp\wedge dq\)
gives the plane residue form divided by \(\sqrt{p_\ell}\).
Exactly as in the
[conic differential calculation](SEXTIC_KELLER_INDEPENDENT_LEADING_PENCIL_RESTRICTIONS.md#the-smooth-line-and-conic-cases),
this excludes a double zero and a reduced zero divisor
of \(p\) on the conic, leaving contact order at least
three at a common basepoint.

Replace \(q\) by a pencil member singular at this
basepoint. It is a product of two lines through
that point and at least one is tangent to the
smooth generic conic. Normalizing the tangent
as \(x=0\) and the smooth member as \(p=xz+y^2\)
gives \(q=x^2\) or \(q=xy\).
For \(q=x^2\), the constant field is \(k(x,p)\).
If the outer binary quadratic pencil has no square
member already, its projection still has a square
member over \(k\): the discriminant of a pencil
of binary quadratics is a homogeneous binary
quadratic on the pencil line and has a zero.
For \(q=xy\), use the leading-pair regularity
comparison on \(x=0\) and \(y=0\). The constant
field is \(k(p,q)\), with \(x\) the free fibre
coordinate. If \(\det\partial(B,C)/\partial(p,q)\)
were nonzero generically along \(q=0\),
the rational primitive \(1/(x\det\partial(B,C)/\partial(p,q))\)
has a pole on \(x=0\) and is regular on \(y=0\).
Both divisors are protected and dominate \(q=0\)
in the constant-field plane, so no correction by
a constant can make the primitive regular on both.
Consequently the coefficients of \(p^2,pq\) in
the two outer forms are proportional as column
vectors. Their pencil contains \(q^2\), and the
other quadratic can be taken as a square after
completing the square and changing \(p\) by a
multiple of \(q\). This last change preserves the
normal form of the conic pencil by a shift of \(z\).

In this frame \(B=p^2,\ C=q^2\), and \((p,q)\)
is onto. Its constant field \(k(p,q)\) and
polynomial contraction put \(A\) in \(k[p,q]\).
Writing \(A\) as the cubic specified in the statement,
the minor \(4pq\) and the other two minors
\(2q A_p,\ 2p A_q\), up to signs, are coprime
if \(a_1a_2\ne0\). The saturation theorem forbids
this.

## Exact area on a quartic plane level

Use the geometric generic fibre throughout, with constant
level \(a\) and smooth affine locus by generic smoothness
in characteristic zero. Its projective equation is
\(K_h(s,u,v)-a s^4=0\) where \(K_h\) is the homogenization.
The relative form is nonzero and regular on the affine
curve. A simple root of the leading quartic gives a
smooth point transverse to infinity at which the form is
regular. In suitable affine coordinates \((s,u)\) at the
heavy root, the equation is \(f(s,u)=0\), and the relative
form is \(s\,ds/f_u\) up to a unit. This follows by
writing the original affine coordinates as \(u/s,1/s\)
in this chart. Their area form has a pole of order three,
and the differential of the level has a pole of order four.

For a branch \(i\), let \(n_i\) be the order of \(s\),
\(\delta_i\) the local normalization length of the branch,
and \(I_{ij}\) its intersection multiplicity with branch
\(j\). The order of the relative form is
\[
                n_i-2\delta_i-\sum_{j\ne i}I_{ij}.
                \tag{2}
\]
Here is a local calculation of this order. Use a branch
parameter \(r\) with \(s=r^{n_i}\), and take a monic
Weierstrass equation \(f_i\) over \(k[[s]]\) for \(u\) on
that branch, with the algebraically closed constant field
of the geometric level understood. The order in \(s\) of
the discriminant of the normalized ring \(k[[r]]\) over
\(k[[s]]\) is \(n_i-1\), using the power basis of \(r\).
Passing to the sublattice generated by the power basis of
\(u\) adds \(2\delta_i\) to the discriminant order. This
discriminant is a signed norm of \((f_i)_u\), so
\(\operatorname{ord}_r(f_i)_u=n_i-1+2\delta_i\).
The other factors of \(f\) add their intersection orders
\(I_{ij}\), giving (2) after taking the order of \(s\,ds\).

Write \(\delta=\sum_i\delta_i+\sum_{i<j}I_{ij}\) for the
total normalization length at the heavy point. This sum
follows by adjoining the branches one at a time and using
the exact sequence of the union with its two summands
and their intersection. There is no other singular point:
the affine curve is smooth and the other infinity point,
if present, is simple. The plane quartic has arithmetic
genus three, so the normalization genus is \(3-\delta\).
In particular \(\delta\le3\). A pole of order one excludes
exactness, as does an everywhere regular nonzero form.

First let \(\sum n_i=3\), so there is a second infinity
point of contact one. For one heavy branch, (2) gives order
\(3-2\delta\), leaving only \(\delta=3\). This branch has
multiplicity three. Indeed multiplicity two and contact
three would give local parameter orders two and three,
with normalization length one. For two heavy branches
the contacts are one and two. If both are smooth their
tangents differ, giving \(I=1\) and no pole. Otherwise
the branch of contact two has multiplicity two and costs
at least one in its normalization length. The intersection
has multiplicity at least two, so the genus bound forces
\(I=2\). The smooth branch then has order \(1-2=-1\).
For three branches, all are smooth, and every intersection
has multiplicity one by the genus bound, again giving
simple poles. The only surviving heavy point has
multiplicity three and one branch.

Now let \(\sum n_i=4\). For one branch (2) gives order
\(4-2\delta\), so again \(\delta=3\). Its multiplicity is
two or three, since multiplicity four on an integral
plane quartic would make the curve a cone of lines.
If the multiplicity is two, the infinity line is tangent.
For two contacts one and three, the smooth branch with
contact one is transverse to infinity. If the other branch
is smooth they have distinct tangents, with intersection
one and no pole. Otherwise the genus bound forces a
multiplicity-two singular branch meeting the smooth
branch with intersection two. This gives a simple pole
on the smooth branch. For contacts two and two with both
branches smooth, the intersection is two or three and
both forms are either regular or have simple poles.
If one branch is singular, the genus bound forces one
smooth branch and a multiplicity-two branch of
normalization length one, with intersection multiplicity
two. This leaves a point of multiplicity three and
\(\delta=3\). With three branches, the genus bound forces
three smooth branches and pairwise intersection one.
Some contact is one and gives a simple pole. Four branches
already have total intersection at least six.

Thus the normalization always has genus zero. At a
point of multiplicity three at \([0:0:1]\), the quartic
equation contains no term of \(z\)-degree above one.
Consequently \(K(y,z)=a(y)z+b(y)\), with \(a\ne0\).
On a generic level the curve is the \(y\)-line with
the zeros of \(a\) removed. At each finite zero of
multiplicity \(i\), the corresponding infinity contact
is \(i\), because \(z\) has a pole of order \(i\).
If the contacts are three and one, \(a\) has degree
three and one distinct root, hence is a cube of a
linear factor. If they are two and two, \(a\) has one
double root. In the remaining case the curve has one
place at infinity, so \(a\) is constant. These give the
first three stated forms after affine changes.

In the multiplicity-two case the tangent cone is the
square of the infinity line. The coefficient of \(z^2\)
in \(K\) is therefore a nonzero constant. The leading
binary quartic is a fourth power of \(y\), so the
coefficient of \(z\) has degree at most two in \(y\).
Completing the square writes \(K\), up to a scalar, as
\((z+a(y))^2+b(y)\), with degrees at most two and four
respectively. If \(\deg b\ge3\), a generic level of this
double cover has positive genus. If \(\deg b=2\),
its generic level has two places at infinity, since
the square completion is a polynomial coordinate
change of the affine plane. If \(b\) is constant the
geometric generic fibre is reducible. This leaves
\(\deg b=1\), giving the last stated form.

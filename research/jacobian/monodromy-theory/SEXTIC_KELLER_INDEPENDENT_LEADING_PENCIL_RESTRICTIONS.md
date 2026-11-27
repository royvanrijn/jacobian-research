# Leading pencil restrictions for independent sextic parts

**Theorem SexticKellerIndependentLeadingPencilRestrictions.**
Let a characteristic-zero three-variable Keller map of ordinary
degree six have three linearly independent degree-six parts.
Over the algebraic closure write its leading tuple
\[
                H_i=g\,\Phi_i(p,q),\qquad 6=\deg g+em,
 \tag{1}
\]
where \(p,q\) are coprime of degree \(m\), \(q/p\) is a closed
parameter in the projective plane function field, and the
\(\Phi_i\) are coprime binary forms of degree \(e\ge2\).
Put \(C_t=\{q-tp=0\}\) for the generic projective member and
\(G=g p^{e-2}\). The following are necessary.

- For \(m=1,2\), the zero divisor of \(G\) on \(C_t\) has
  partition \((3,1)\) or \((4)\).
- For \(m=1\), writing \(B=\{p=q=0\}\), one has
  \(e-2+\operatorname{mult}_B(g)\ge3\), or \(g\) has a line
  of multiplicity at least three that does not pass through \(B\).
- For \(m=2\), the conic pencil has a common tangent at a
  fixed basepoint. One can take \(p=xz+y^2\) and \(q\) binary
  quadratic in \(x,y\). If \(e=3\), one can further take
  \(q=x^2\) or \(q=xy\). In the latter case the projective
  outer map \([\Phi_1:\Phi_2:\Phi_3]\) is nonimmersive at
  \([p:q]=[1:0]\). If \(e=2\) and \(g\) is an irreducible
  conic, one can further take \(q=x^2\).
- For \(m=3\), the cubic pencil is generically cuspidal at a
  fixed basepoint. One can take \(p=x^2z+c(x,y)\) and
  \(q=b(x,y)\) with binary cubics \(b,c\), with a repeated factor
  in \(q\).
  If \(x\nmid q\), a pencil change and an affine source
  change allow \(p=x^2z\).

These are conditions on the independent leading parts, with
no restrictions on the lower degree coefficients. Another
condition for unequal component degrees is the
[higher-member basepoint obstruction](KELLER_LEADING_PAIR_RATIONAL_SLICE.md)
when the leading pair is coprime.

**Theorem SexticKellerBirationalBinaryLeadingTuple.**
In (1), suppose \(m=1\), \(e\ge4\), and the outer map
\([\Phi_1:\Phi_2:\Phi_3]\) is birational onto its image.
If \(e=5\) or \(e=6\), such leading parts cannot occur.
If \(e=4\), the Keller map is invertible.

**Theorem SexticKellerCoprimeBinaryLeadingReduction.**
In (1), suppose \(m=1\) and \(e=6\). Either a triangular
source change gives an ordinary degree at most six map
whose degree-six parts span a space of dimension at most
one, or in a polynomial source frame with weights
\((1,1,2)\) the weight-six tuple is a constant target change
of \((P^2,PQ,Q^2)\), where
\[
                     P=xz+y^3,\qquad
                     Q\in\{x^3,x^2y,xy^2\}.
\]
In the second alternative all terms have weight at most six.

**Theorem SexticKellerClosedCubicVeroneseExclusion.**
Take ordinary weights or weights \((1,1,2)\). Let \(P,Q\)
be coprime weight-three forms such that the polynomial
kernel of \(D=\det J(P,Q,\cdot)\) is \(k[P,Q]\).
Then a polynomial map with leading weight-six tuple
\((P^2,PQ,Q^2)\) and no higher-weight terms cannot be Keller.

In the coprime binary leading reduction, this excludes the
two weighted cases \(P=xz+y^3,\ Q=x^2y\) or \(xy^2\).
In the ordinary \(m=3\) case of (1), the only pencils
left are, up to linear source and pencil changes,
\[
 \langle x^2z+y^3,x^3\rangle,\qquad
 \langle x^2z,y^3\rangle .
\]

**Theorem SexticKellerWeightedCubePencilExclusion.**
Take weights \((1,1,2)\) and highest tuple
\((P^2,PQ,Q^2)\) with \(P=xz+y^3\) and \(Q=x^3\).
There is no Keller completion. Consequently every Keller
map in the coprime binary leading reduction can be changed
to ordinary degree at most six with sextic leading rank at
most one.

**Theorem SexticKellerCubicOuterLinePencilReduction.**
Suppose \(m=1,e=3\) in (1). A common factor depending on
the third source direction is impossible. If the common
factor is binary, a triangular source change either lowers
the ordinary degree to at most five or gives a degree at
most six map with sextic leading rank at most one.

**Theorem SexticKellerBinaryQuarticCommonFactorReduction.**
Suppose \(m=1,e=2\) in (1) and the quartic common factor
is binary in \(p,q\). A triangular source change gives an
ordinary degree at most six map with sextic leading rank
at most one, unless in a polynomial source frame the
highest tuple for weights \((2,2,3)\) is
\[
                (P^2,PQ,Q^2),\qquad P=x^2y+z^2,\quad Q=x^3 .
\]
In this exception the binary quartic common factor in
the original leading pencil is a fourth power.

## Reduction to a rational differential on the pencil

Work over \(\mathbb C\) by embedding a finitely generated
coefficient field. The top determinant equation gives
\(\det JH=0\). Apply
[Wang–Yan, Lemma 2.4](sources/wang-yan-2610.36495v1.pdf)
to the independent forms of common degree six. Its factorization
proof takes a parameter on the normalization of the projective
image curve. Replace that parameter by a closed generator of its
relative algebraic closure in \(\mathbb C(\mathbb P^2)\), using
[Petravchuk–Iena, Lemma 2](sources/petravchuk-iena-math-0701588v2.pdf).
Express the new parameter as \(q/p\) using coprime forms of equal
degree. The projective image coordinates are again binary forms
without a common projective zero, now evaluated at \(p,q\).
Their substituted forms have no common factor because \(p,q\)
are coprime. This proves (1) in the closed frame. The degree
relation gives \(m\in\{1,2,3\}\). Closedness gives a geometrically
integral generic pencil member: its function field is the
regular extension \(\mathbb C(\mathbb P^2)/\mathbb C(t)\),
and the pencil has no fixed curve.

The [leading-pair slice theorem](KELLER_LEADING_PAIR_RATIONAL_SLICE.md)
implies exactness of the relative volume differential on a
geometric generic fibre of \((H_1,H_2)\).
In an affine projective chart \(\ell=1\), write source coordinates
as a radial variable \(s\) times the two projective variables.
Use the same symbols for the dehomogenized forms. Then
\[
           H_1=s^6 g p^e\Phi_1(1,t),\qquad H_2=H_1 f(t),
           \qquad f=\Phi_2(1,t)/\Phi_1(1,t).
\]
The ratio \(f\) is nonconstant. At generic fixed \(H_1,H_2\),
\(t\) is fixed after extending constants, and the fibre is a
radical cover of \(C_t\). Let \(\eta\) be the plane residue
of the coordinate two-form divided by \(q-tp\) in this chart.
Since \(dH_1\wedge dH_2=H_1 f'(t)\,dH_1\wedge dt\),
\(\partial_sH_1=6H_1/s\), and
\(dt=d(q-tp)/p\) in the projective variables at that fixed \(t\),
the relative volume differential is a nonzero constant multiple of
\[
                    s^3 p\,\eta = \mathrm{const}\cdot
                             \frac{\eta}{\sqrt{g p^{e-2}}}.
 \tag{2}
\]
The \(s^2\) in the ambient three-form produces the third power
of \(s\) in this expression. The square root on a chosen
component of the cover is supplied by
\(s^6 g p^e=\mathrm{const}\ne0\).

Thus \(\eta/\sqrt G\) must have a rational primitive on the
normalization of the quadratic cover of \(C_t\), taking one
component if \(G\) is square. Indeed if a differential from
a lower field is exact after a finite separable extension,
take the field trace of a primitive and divide by the degree.
The trace commutes with differentiation. All arguments here
can be made after adjoining an algebraic closure of the generic
fibre constants.

## The smooth line and conic cases

If \(m\le2\), \(C_t\) is smooth and \(G\) is a section of
degree \(6-2m\) in plane coordinates, so its zero divisor
on \(C_t\) has degree four. Choose the chart line \(\ell=0\)
transverse to \(C_t\) and disjoint from those zeros. The
residue \(\eta\) is nonzero regular in the affine smooth
chart. At the chart boundary it has order \(m-3\), since the
plane two-form has a pole of order three and the dehomogenized
equation contributes a pole of order \(m\). On the quadratic
cover, \(1/\sqrt G\) has order \(3-m\) there, so (2) is
regular above the chart boundary.

At a zero of multiplicity one of \(G\), the normalized
double cover ramifies and the differential is regular:
a base parameter \(u\) becomes a square, and \(du/\sqrt u\)
is regular. A zero of multiplicity two gives a simple pole
with nonzero residue on each unramified branch. A rational
primitive is therefore impossible if any zero has multiplicity
two, or if all four zeros are simple, since in the latter case
the differential is nonzero and everywhere regular.
This leaves \((3,1)\) and \((4)\).

In a line pencil, the multiplicity at the vertex \(B\) is
\(e-2+\operatorname{mult}_B g\). Away from \(B\), a reduced
irreducible factor of \(g\) not contained in a pencil member
meets a generic line transversely. Indeed the pencil parameter
gives a nonconstant separable function on the normalization
of that factor, and its general fibres are reduced and avoid the
singular points. Distinct factors cannot meet each other away
from \(B\) on a generic line. Thus a zero of order at least
three off \(B\) must come from a factor of \(g\) of that
multiplicity, and this factor has degree one because
\(\deg g\le4\).

In a conic pencil, \(G\) is a plane quadratic. A point of
intersection multiplicity at least three must be a basepoint
by the same separability argument, since no factor of \(G\)
has multiplicity greater than two. At such a basepoint the
tangent of the generic \(C_t\) is fixed. This follows from
tangency to \(G\) if \(G\) is irreducible, tangency to the
underlying line if \(G\) is double, or tangency to one of its
two lines if \(G\) is split and reduced. In the last case one
of the two tangencies holds generically. Choose a smooth
member \(p\) and subtract from the other member a scalar
multiple of \(p\) to make it singular at the basepoint.
Put the point at \([0:0:1]\) and its tangent at \(x=0\).
The resulting pencil has \(p=xz+y^2\) after rescaling
and shifting \(z\) by a binary linear form, and \(q\)
a binary quadric. When \(e=3\), the zero divisor of \(G=p\)
on \(C_t\) is the base intersection. Hence \(q\) must
contain the tangent line \(x=0\): a singular conic at the
basepoint consists of two lines through it, and their
total contact would otherwise be two. This gives \(q=x^2\)
or \(q=xy\) after a linear change and the corresponding
renormalization of \(p\).

## Regularity restrictions in the tangent conic cases

First take \(e=3,\ q=xy,\ p=xz+y^2\). Then
\(\det J(p,q,1/x)=1\) and
\(\mathbb C(x,y,z)=\mathbb C(p,q)(x)\).
If the projective outer map has nonzero derivative at
\([1:0]\), choose two constant combinations \(A,B\) of the
leading components such that
\(\det \partial(A,B)/\partial(p,q)\) is nonzero generically
on \(q=0\). Indeed, on \(p\ne0\) the radial coordinate is
\(p^3\) times a nonvanishing outer coordinate, and one
projective coordinate has nonzero derivative in \(q/p\).
Then \(dA,dB\) are independent at the generic points of both
\(x=0\) and \(y=0\). The function
\[
                \frac{1/x}{\det\partial(A,B)/\partial(p,q)}
\]
is a primitive for their rational Jacobian derivation, up to a
nonzero scalar. It has a pole along \(x=0\), and is regular
along \(y=0\). Both divisors map dominantly to \(q=0\) in
the constant-field plane. The simultaneous regularity in the
leading-pair theorem is therefore impossible.

Next take \(e=2\) with \(g\) an irreducible conic. In coordinates
\(g=xz+y^2\) and \(B=[0:0:1]\), the conics having intersection
order at least three with \(g\) at \(B\) form the space
\(\langle g,x^2,xy\rangle\). This follows by substituting
\((x,y,z)=(-r^2,r,1)\) and killing the coefficients of degrees
zero through two. The entire pencil belongs to this space.
Choose its smooth member as \(p=xz+y^2\) in a new frame by
a shift in \(z\). The resulting pair has
\[
          q=x(\gamma x+\delta y),\qquad
          g=\mu p+x(\nu x+\rho y),\qquad \mu\ne0.
\]
Suppose \(\delta\ne0\). Up to a target change choose
\(A=S=gp^2,\ B=gpq=ST\), \(T=q/p\).
In the chart with \(r=y/x,\ h=z/x\), we have
\(p_\ell=h+r^2,\ q_\ell=\gamma+\delta r\).
At fixed \(T\), the derivative of \(g_\ell\) with respect to
\(r\) is \(L=\mu\delta/T+\rho\). Directly in this chart,
\(\det J(S,ST,r)=6S^2T/(x^3p_\ell)\).
The relation \((x^3 p_\ell)^2g_\ell=S\) then gives
\[
                \det J\left(S,ST,\frac{pg/x}{3S^2(\mu\delta+\rho T)}\right)=1 .
 \tag{3}
\]
Indeed \(pg/x=x^3p_\ell g_\ell\), whose \(r\)-derivative
holding \(S,T\) fixed is \(L x^3p_\ell/2\).

The two divisors \(x=0\) and \(\gamma x+\delta y=0\)
have \(T=0\) with a simple zero. Along each, the restriction
of \(S\) is nonconstant (on the second it has leading term
\(\mu p^3\) as a polynomial in \(z\)). Thus \(dA,dB\) are
independent at both generic points. The primitive in (3) has
a pole along the first divisor and is regular along the second.

If \(g\notin\langle p,q\rangle\), the constant field is
\(\mathbb C(S,T)\). In fact, at fixed \(T\) the affine linear
functions \(g_\ell,q_\ell\) have distinct roots, and the
equation \(x^6=S/(g_\ell p_\ell^2)\) gives a geometrically
integral curve by Eisenstein applied to \(1/x\) at the
simple root of \(g_\ell\). Both divisors dominate the same
curve \(T=0\) in this constant-field plane. If
\(g\in\langle p,q\rangle\), the constant field is
\(\mathbb C(p,q)\) because \(\mathbb C(x,y,z)=\mathbb C(p,q)(x)\)
and \(S,T\) generate a subfield of finite index there.
Both divisors now dominate \(q=0\). In either case (3)
contradicts simultaneous regularity, so \(\delta=0\).

## The cubic pencil

When \(m=3\), \(G\) is constant and (2) is the dualizing
residue of the plane cubic. If the generic cubic is smooth,
the differential is nonzero holomorphic. If it is nodal,
the residue has a simple pole with nonzero residue on each
normalization branch at the node. This follows by taking
local equations of the two transverse branches as coordinates
in the plane. Both cases preclude a rational primitive.
A singular geometrically integral plane cubic with no node
is cuspidal. Its generic singular point belongs to the finite
base locus by Bertini, hence is fixed, and both pencil members
are singular there.

Put that point at \([0:0:1]\). The homogeneous quadratic jets
of the two cubics at the point span a vector space consisting
of squares of linear forms, including zero. Two independent
binary squares have a sum of rank two, so this space has
dimension one. After choosing a basis and a linear coordinate
\(x\) in the tangent direction equation, we get
\(p=x^2z+c(x,y),\ q=b(x,y)\) as stated.

Suppose first \(q(0,1)\ne0\). Set
\(\lambda=c(0,1)/q(0,1)\), and write
\(p-\lambda q=x(xz+h_2(x,y))\).
If \(x^2\) does not divide \(p-\lambda q\), the two
divisors \(x=0\) and \(xz+h_2=0\) are distinct,
the second irreducible since \(h_2(0,1)\ne0\).
Both dominate \(p-\lambda q=0\) in the \(p,q\)-plane.
The forms \(d p,d q\) are independent at their generic
points. At the first this follows from the simple transverse
zero of \(p-\lambda q\) and the nonconstant binary cubic
\(q(0,y)\); at the second it follows from \(p_z=x^2\) and
the nonvanishing gradient of \(q\) generically.

Take the leading pair \(A=q^2,\ B=pq\), which has independent
differentials at both points, since \(q\ne0\) generically there.
As \(\det J(p,q,y/x)=3q\), a primitive for the leading
derivation is \(-y/(6xq^3)\).
It has a pole at the first divisor and is regular at the second.
The constant field is \(\mathbb C(p,q)\) unless \(q\) is
a cube of a linear form: over \(p,q\) the equation is
\(x^3=q/b(1,r)\) for \(r=y/x\), and is geometrically
irreducible unless the binary cubic is a cube.
In the cube case \(q=\ell^3\), the field is
\(\mathbb C(p,\ell)\), since \(\ell\) is independent of \(x\)
and \(x\) is free. In both cases the divisors dominate the
same curve in the constant-field plane, so the two regularity
conditions contradict each other. This forces
\(x^2\mid p-\lambda q\); replacing \(p\) and shifting \(z\)
then gives \(p=x^2z\).

Finally suppose \(q\) is squarefree. Choose two distinct
linear factors \(y-\beta_i x\) of \(q\) with neither equal
to \(x\). On both divisors \(p_z=x^2\) and \(dq\) are
generically nonzero, so for \(A=p^2,\ B=pq\), the leading
pair has independent differentials there. A primitive is
\(\phi_0=y/(6x p^2 q)\), and its constant field is
\(\mathbb C(p,q)\) by the same cubic-cover argument. Both
divisors dominate \(q=0\) and the order of \(q\) is one at
each. For \(\phi_0+h\) with \(h\in\mathbb C(p,q)\) to be
regular at both, the residue of \(qh\) as a rational
function on \(q=0\) would have to be both
\(-\beta_1/(6p^2)\) and \(-\beta_2/(6p^2)\). This is
impossible.

## Binary pencils with a birational outer map

First suppose \(g\) is binary, and use \(p=x,q=y\) so all
three degree-six forms are independent of \(z\). Some
component must have a monomial involving \(z\), since the
map has nonzero determinant. Give \(x,y\) weight one and
\(z\) weight
\[
   w=\min\left\{\frac{6-a-b}{c}:
       x^a y^b z^c\ \hbox{occurs in some component},\ c>0\right\}.
\]
Then \(1<w\le6\), all terms have weight at most six, and
at least one weight-six term involves \(z\). Write \(L_i\)
for the weight-six parts. Since \(18-(2+w)>0\), the
weight-\(18-(2+w)\) determinant equation is \(\det JL=0\).
But \(L(x,y,0)=H(x,y)\), whose image is dense in the
two-dimensional affine cone over the projective outer image
curve. It follows that the closure of the image of \(L\)
is exactly that cone.

By birationality the parameter of this projective outer image
is a rational function of the coordinate ratios of \(L\).
Write its pullback as \(Q/P\) with coprime polynomials. It
is invariant under weight scaling, so \(P,Q\) are weighted
homogeneous of the same degree \(D\). Indeed, after clearing
weight denominators, invariance under each nonzero scaling
and coprimeness imply that each polynomial transforms by
the same scalar, which forces all its monomials to have
the same weight. Neither \(P\) nor \(Q\) vanishes
identically on \(z=0\): if both vanished they would have
the common factor \(z\), and if exactly one vanished the
projective image on \(z=0\) would be constant. Thus \(D\)
is a positive integer, by restriction to the binary plane.

The evaluated forms \(\Phi_i(P,Q)\) have no common
irreducible factor, since \(P,Q\) are coprime and the
binary forms have no common projective zero. Therefore the
projective equality gives \(L_i=G\Phi_i(P,Q)\) for a
weighted homogeneous polynomial \(G\), and \(6=\deg_w G+eD\).
For \(e\ge4\) this gives \(D=1\); hence \(P,Q\) are binary
linear forms. All the \(z\)-dependence of \(L\) is in \(G\).
For \(e=5,6\), the weight of \(G\) is at most one and
this is impossible. For \(e=4\), its weight is two,
so \(w=2\) and \(G=g_2(x,y)+cz\) with \(c\ne0\).
Replace \(z\) by the source coordinate \(Z=G\). This is
a triangular polynomial automorphism, and in the weights
\((1,1,2)\) it preserves the weight filtration. The
weight-six terms are now \(Z\Phi_i(P,Q)\) of ordinary
degree five, and all terms of lower weight have ordinary
degree at most five. Apply
[three-variable inversion through degree five](THREE_VARIABLE_KELLER_INVERSION_THROUGH_DEGREE_FIVE.md)
to this Keller map. Undoing the source automorphism and
descending the unique inverse to the original field gives
the conclusion.

Now suppose \(g\) depends on \(z\). Let \(r=\deg_z g\ge1\)
and let \(\beta(x,y)\) be its leading \(z\)-coefficient.
Here \(e=4\) or \(5\) and \(1\le r\le6-e\).
Suppose some monomial \(x^a y^b z^c\) has \(c>r\), and give
\(z\) weight
\[
   w=\min\left\{\frac{6-r-a-b}{c-r}:
       x^a y^b z^c\ \hbox{occurs in some component},\ c>r\right\}.
\]
The ordinary degree bound and the shape of the sextic
parts give \(w>1\). The highest weighted tuple \(L\) has
weight \(E=6-r+rw\) and includes a term with \(z\)-exponent
greater than \(r\). Every term with \(z\)-exponent less than
\(r\) has weight below \(E\), so \(z^r\) divides every
component of \(L\). The weight-\((3E-2-w)\) equation of
the Jacobian gives \(\det JL=0\). Its affine image closure
is a cone of dimension at most two. Its projective
image equals the projective image closure of \(M=L/z^r\),
which contains the outer image curve by
\(M(x,y,0)=\beta(x,y)\Phi(x,y)\).
The same factorization gives \(M_i=G\Phi_i(P,Q)\) and
\(6-r=\deg_w G+eD\) with \(D\) a positive integer.
Thus \(D=1\), and \(G\) has weight at most one. This again
makes \(M\) independent of \(z\), a contradiction.

The map consequently has degree \(r\) in \(z\). Denote its
leading \(z\)-coefficient by \(A(x,y)\). Its highest
binary parts are \(\beta(x,y)\Phi(x,y)\). The coefficient
of \(z^{3r-1}\) in the Jacobian determinant gives
\(r\det(A_x,A_y,A)=0\). On the chart of any nonzero
component of \(A\), this is the vanishing of the
two-variable Jacobian of its other two coordinate ratios,
up to a nonzero factor. The projective image of \(A\)
therefore lies in a curve and is dense in the outer image
curve by radial scaling. Birationality and coprimeness give
\(A_i=G\Phi_i(P,Q)\) with coprime polynomials in \(x,y\).
For \(D=\max(\deg P,\deg Q)\), the degree identity is
\(6-r=\deg G+eD\), because the degree-\(D\) parts of
\(P,Q\) cannot give zero in all evaluated binary forms.
Thus \(P,Q\) are affine linear with independent linear
parts, by the nonconstant highest projective tuple.
At their common zero \(A=0\). For \(r=1\), this already
gives \(F_z=0\), a contradiction.

For \(r=2\) one has \(e=4\), so \(G\) is constant. Center
the common zero at the origin, making \(A\) homogeneous
quartic, and write the map as \(Az^2+Bz+C\).
The coefficient of \(z^4\) in the determinant is
\[
   \det(A_x,A_y,B)+2\det(A_x,B_y,A)+2\det(B_x,A_y,A)=0.
\]
Its binary degree-six part is \(\det(A_x,A_y,B(0))=0\).
If the constant vector \(B(0)\) were nonzero, this would
force the outer curve to have zero derivative under
projection from \([B(0)]\). To see this, take
\(B(0)=(0,0,1)\) by a constant target change and use Euler's
identity to conclude that \(A_1/A_2\) is constant on the
parameter line. The projective curve would lie in a line,
contrary to independence of the leading tuple. Thus
\(B(0)=0\), and \(F_z=2Az+B\) vanishes on the entire
line above the common zero. This excludes the nonbinary
common factor also for \(e=4\).

## Composite coprime binary sextics

For \(m=1,e=6\) the common factor is constant. The preceding
theorem excludes a birational outer map. Parametrize the
normalization of the projective image. The binary sextic
tuple now has the form \(H_i=\Psi_i(B,C)\), with coprime
binary forms \(B,C\) of degree \(u>1\) and a birational
outer tuple of degree \(v\ge2\). Both tuples are free of
common projective zeros, so \(uv=6\).

Use the first weight threshold of the preceding proof, with
\(z\) of weight \(w>1\). The tuple \(L\) of weight-six
parts has the same projective image as \(H\), so
\(L_i=G\Psi_i(P,Q)\) for coprime weighted forms \(P,Q\)
of degree \(D\). On \(z=0\), their ratio restricts to
\(C/B\). Hence \(D\ge u\). Since \(6=\deg_w G+vD\),
one has \(D=u\), \(G\) constant, and the restrictions of
\(P,Q\) are \(B,C\) up to a common scalar.

If \(u=2\), the only possible \(z\)-dependence is linear
with constant coefficients and \(w=2\). If \(u=3\),
the possibilities are linear with constant coefficients
at \(w=3\), linear with binary linear coefficients at
\(w=2\), and quadratic with constant coefficients at
\(w=3/2\). In either constant-coefficient linear case,
change the pencil basis to make \(Q\) binary and take
\(Z=P\) as a new triangular source coordinate. The weight
filtration is preserved. Now only the term \(Q^v\) in
the weight-six tuple can have ordinary degree six, while
all terms of smaller weight have ordinary degree below
six. This gives the first alternative.

For \(u=3\) the binary quadratics \(\Psi_i\) span all
binary quadratics, so take the leading tuple to be
\((P^2,PQ,Q^2)\). The weighted leading-pair slice for
\((P^2,PQ)\) implies exactness of the volume differential
relative to \(P,Q\) on each geometric generic fibre:
the two volume differentials differ there by an invertible
constant from the two-by-two target change. We use
constants \(p_0,q_0\) for the generic values of \(P,Q\)
and work over the algebraic closure of their field.

For \(w=3/2\), choose \(P=B+cz^2,\ Q=C\) with \(c\ne0\).
In coordinates \(x,t=y/x,z\) the ambient volume form
is \(x\,dx\,dt\,dz\). The identities
\[
           x^3 C(1,t)=q_0,\qquad
           cz^2=\frac{p_0 C(1,t)-q_0 B(1,t)}{C(1,t)}
\]
give the relative differential, up to a nonzero constant,
as
\[
           \frac{dt}{C(1,t)^{1/6}
                     (p_0 C(1,t)-q_0 B(1,t))^{1/2}}.
\]
The two denominator forms have disjoint zero sets, and
the generic binary cubic in the second factor is squarefree
by separability of the basepoint-free pencil. The exponents
at their zeros are all at most \(1/2\), so the form is
regular on the normalized fibre above those points. In
fact a base expression of order \(-a\) with \(a<1\)
gives integral order \(s(1-a)-1>-1\) after ramification
of index \(s\). At a chart boundary chosen off the zero
sets, the two factors have total exponent \(3/6+3/2=2\)
in the pole of \(t\), cancelling the pole of \(dt\).
The form is nonzero and everywhere regular, contradicting
exactness.

For \(w=2\), first take the two binary linear coefficients
of \(z\) to be independent, and use them as \(x,y\), so
\(P=B+xz,\ Q=C+yz\). Put \(R=xC-yB\), a nonzero binary
quartic by coprimeness of \(B,C\). Eliminating \(z\) gives
\(x^3 R(1,t)=q_0-p_0t\). Dividing the ambient volume first
by \(P_z=x\) and then by \(3x^2R(1,t)\) gives a relative
differential proportional to
\[
                 \frac{dt}{R(1,t)^{1/3}(q_0-p_0t)^{2/3}}.
\]
As before it is regular at a chart boundary away from the
zeros, since \(4/3+2/3=2\).
A triple zero of \(R\) gives a simple pole with nonzero
residue, while multiplicities at most two throughout
give an everywhere regular nonzero form. Thus \(R\)
is a fourth power up to a scalar. By simultaneous linear
changes of \((x,y)\) and of the pencil basis \((P,Q)\),
put \(R=\alpha x^4\), \(\alpha\ne0\).
Then \(B=xS,\ C=yS+\alpha x^3\) with \(S\) binary
quadratic. Set \(Z=z+S\). In this frame \(P=xZ,\ Q=yZ+
\alpha x^3\), so only \(Q^2\) can give a term of ordinary
degree six. This again yields the first alternative.

In the remaining \(w=2\) case, choose
\(P=B+xz,\ Q=C\). The relative differential is proportional
to \(dt/C(1,t)^{1/3}\) on \(x^3C(1,t)=q_0\), as follows
from the same two divisions. If \(x\nmid C\), then \(C(1,t)\)
has degree three and the form has a simple pole at \(t=\infty\).
If \(x\) divides \(C\) exactly once and \(C\) is squarefree,
then \(C(1,t)\) is squarefree quadratic. The normalized
cubic cover has genus one by Riemann–Hurwitz, with three
points of full ramification over its two roots and infinity.
The differential is regular at all finite points and has
one double pole at infinity. A rational primitive would
have a single simple pole and hence give a degree-one map
to the projective line, impossible in genus one.

It follows that \(x\mid C\) and \(C\) has a repeated factor.
Up to rescalings and a binary linear change fixing the
line \(x=0\), \(C\) is \(x^3,x^2y\) or \(xy^2\).
Since \(B,C\) are coprime, the coefficient of \(y^3\)
in \(B\) is nonzero. A weight-preserving shift of \(z\)
by a binary quadratic absorbs all the terms of \(B\)
divisible by \(x\), and rescaling gives \(P=xz+y^3\).
This proves the second alternative.

## Cubic Veronese tuples with a closed constant ring

Write \(\nu(p,q)=(p^2,pq,q^2)\) and
\(n=(Q^2,-2PQ,P^2)\). Wang–Yan use the leading adjugate
identity and the constant ring to exclude ordinary
quartic completions of a nonsquare conic pencil.
[Wang–Yan, Proposition 3.9](sources/wang-yan-2610.36495v1.pdf)

Here let the current reference tuple be \(G(\widetilde
P,\widetilde Q)\), with \(\widetilde P,\widetilde Q\)
having weight-three parts \(P,Q\) and \(G\) having
quadratic part \(\nu\). Start with \(\widetilde P=P\),
\(\widetilde Q=Q\) and \(G=\nu\). Write \(E\) for the
difference of the map from the reference tuple, and let
\(E_j\) be its highest weighted part if its weight \(j\)
is positive. The reference determinant vanishes because
it factors through two polynomials. The weight
\(12+j-\sum_i w_i\) determinant coefficient of the
full tuple is
\[
                         2D(n\cdot E_j).
\]
Indeed the terms with two reference Jacobian factors and
one error factor at this weight use just the leading
reference tuple \(\nu(P,Q)\). The two partial columns
of \(\nu\) have cross product \(2n\), and \(D(n)=0\).
Terms with at least two error factors have lower weight
since \(j<6\). The displayed weight is positive, so
the Keller equation gives \(n\cdot E_j\in k[P,Q]\).

Unless \(j=3\), this invariant must be zero, since it has
weight \(6+j\), and \(P,Q\) have weight three. If \(j=3\)
the invariant is a binary cubic in \(P,Q\). Every such
cubic is \(n\cdot M(P,Q)\) for a constant linear map
\(M\) from two coordinates to three. Take the first
entry to supply the \(PQ^2,Q^3\) coefficients and the
third entry to supply the \(P^3,P^2Q\) coefficients.
Absorb \(M(\widetilde P,\widetilde Q)\) into the reference
tuple, so the remaining error of weight \(j\) is tangent
also in this case.

For an error vector \((a,b,c)\) of weight \(j\) satisfying
\(n\cdot(a,b,c)=0\), coprimeness yields \(P\mid a\) and
\(Q\mid c\) by reduction modulo \(P\) and \(Q\).
Writing \(a=2P U,\ c=2Q V\) then gives
\(b=Q U+P V\). If \(j<3\), this forces zero error. If
\(3\le j\le5\), replace \(\widetilde P,\widetilde Q\)
by \(\widetilde P+U,\widetilde Q+V\). Their highest parts
are unchanged because \(U,V\) have weight \(j-3\).
The weight-\(j\) change of the reference is exactly the
tangent vector, since \(2(j-3)<j\) and the linear terms
of \(G\) contribute only at lower weight. Thus the
highest positive error weight strictly decreases. After
absorption of the constant part into \(G\), the original
map factors through \(\widetilde P,\widetilde Q\) and has
zero determinant.

For the applications work over the algebraic closure.
In both weighted cases \(Q=x^2y,\ xy^2\), the generic
curve with \(P,Q\) fixed is rational over \(k(P,Q)\):
use \(x\) as the free parameter for \(x^2y\) and \(y\)
for \(xy^2\). Recover the other binary variable from
\(Q\) and solve for \(z\) from \(P\). Thus the constant
field of \(D\) is \(k(P,Q)\). The map \((P,Q)\) is onto
the affine plane. For \(Q\ne0\) take \(x,y\ne0\) with
the prescribed \(Q\); for \(Q=0\) take \(x\ne0,\ y=0\).
Then solve for \(z\) in both cases.

For the ordinary \(m=3\) case use \(P=x^2z+c(x,y)\)
and \(Q=b(x,y)\) as above. If \(Q\) is not a cube of a
linear form, the generic curve is geometrically integral.
Indeed on \(x\ne0\) it is the cubic cover
\(x^3=Q/b(1,t)\), \(t=y/x\), after solving for \(z\).
The binary cubic has a projective zero of order not
divisible by three, so the right side of the cover
equation is not a cube over the algebraic closure of
the generic constants. The plane \(x=0\) does not
dominate the \(P,Q\)-plane: there both polynomials are
scalar multiples of \(y^3\). The constant field of \(D\)
is again \(k(P,Q)\). This map is onto. For \(Q\ne0\)
choose a direction with \(x\ne0\) and \(b\ne0\) and
scale it to the prescribed value. For \(Q=0\), use a
zero direction of \(b\) with \(x\ne0\), which exists
because \(b\) is not a cube. In either case solve for \(z\).

In these applications surjectivity gives the polynomial
constant ring \(k[P,Q]\). If a reduced quotient
\(A(P,Q)/B(P,Q)\) is polynomial and \(B\) has a
nonconstant factor, choose a point on its zero curve
outside \(A=0\) and lift to the source. The cross-multiplied
equality at the lift is impossible.

It remains in the ordinary cubic pencil to take \(Q=\ell^3\).
If \(x\nmid Q\) the previous cubic reduction gives
\(P=x^2z\), and one can take \(\ell=y\).
If \(x\mid Q\), take \(Q=x^3\). Coprimeness gives
a nonzero \(y^3\) coefficient in \(c\). Shift \(y\) by a
multiple of \(x\) to remove the \(xy^2\) term after
rescaling, and shift \(z\) by a binary linear form to
absorb the terms of \(c\) divisible by \(x^2\). This
gives \(P=x^2z+y^3\).

## Line pencils with a cubic outer map

The outer map is birational onto its image in this case.
Its degree onto the image times the degree of the image
is three, and the image is not a line by independence.
We use coordinates \(p=x,q=y\).

First suppose the cubic common factor is binary. Apply
the first weight threshold of the binary pencil proof.
Its tuple is \(L_i=G\Phi_i(P,Q)\), with
\(\deg_w G+3D=6\), where \(D\) is the common positive
integer weight of \(P,Q\). If \(D=2\), then \(G\)
is constant and the \(z\)-dependence is linear with
constant coefficients at \(w=2\). Change basis to make
\(Q\) binary and use \(P\) as a new source coordinate.
Only the \(Q^3\) term can have ordinary degree six, so
this gives leading rank at most one.

If \(D=1\), take the pencil coordinates to be \(x,y\);
the weight of \(G\) is three. At \(w=3\) one has
\(G=g(x,y)+cz\), \(c\ne0\), so the source coordinate
\(Z=G\) lowers ordinary degree to at most five. At
\(w=2\) one has \(G=g(x,y)+z\ell(x,y)\) for a nonzero
linear form \(\ell\). Set \(\ell=x\). A shift of \(z\)
by a binary quadratic removes the terms of \(g\)
divisible by \(x\). If \(g\) was divisible by \(x\)
all the weight-six terms now have ordinary degree at
most five. Otherwise a rescaling gives
\(G=xz+y^3\). This weighted shape is excluded below.

The only other weight is \(w=3/2\), giving
\(G=g(x,y)+cz^2\), \(c\ne0\). For a generic fibre of
a pair of leading components, the projective line
parameter \(t=y/x\) is fixed after extending constants.
On a component the remaining equation is
\[
                  x^3(cz^2+x^3g(1,t))=a_0,\qquad a_0\ne0.
\]
The two target coordinates are \(A\) and \(A\) times a
nonconstant rational function of \(t\). Dividing the
volume \(x\,dx\,dt\,dz\) by their differentials gives
\(dx/(x^2z)\) up to a nonzero constant on the curve.
Thus with \(u=x^2z\) the differential is \(dx/u\) on
\(c u^2=a_0x-g(1,t)x^7\). The right-hand polynomial
has seven simple roots because \(g(1,t)\ne0\) for generic
\(t\). The differential is regular at the finite branch
points and also at infinity, where \(x\) has a double
pole and \(u\) a pole of order seven. This contradicts
the weighted rational-slice theorem.

Now suppose the common factor has positive \(z\)-degree
\(r\le3\). Its leading \(z\)-coefficient is a nonzero
binary form \(\beta\) of degree \(3-r\). Use the
threshold \(w=\min(6-r-a-b)/(c-r)\) over monomials
\(x^a y^b z^c\) in the map with \(c>r\), if any.
As in the argument for a nonbinary common factor above,
the weight-\((6-r+rw)\) tuple is divisible by \(z^r\).
The projective image of that tuple is the outer image
curve, as seen by restricting the quotient by \(z^r\)
to \(z=0\). Factor the quotient as \(G\Phi_i(P,Q)\)
with \(\deg_w G+3D=6-r\). Hence \(D=1\) and \(G\)
has weight \(3-r\). If \(r\ge2\), this contradicts
the \(z\)-dependence of the quotient. If \(r=1\),
the only possible \(z\)-dependence is at \(w=2\),
and the weight-seven tuple is \(z(g_2(x,y)+cz)\Phi_i(x,y)\)
with \(c\ne0\). A second threshold for any \(z\)-exponent
above two must have \(w>2\), since those monomials have
weight strictly below seven at \(w=2\). At that
threshold the leading tuple is divisible by \(z^2\)
and the quotient has weight three, restricts to
\(c\Phi_i\), and factors through the same outer map.
Therefore the quotient is binary, contradicting the
chosen monomial. This proves that the map has degree
\(s=r\) in \(z\), except that \(r=1\) can give \(s=2\).

Let \(A\) be its leading \(z\)-coefficient. As above,
the \(z^{3s-1}\) Jacobian coefficient gives
\(\det(A_x,A_y,A)=0\), so the projective image of \(A\)
is a curve. Radial scaling shows it is the same outer
image curve, since the highest binary tuple of \(A\)
is a nonzero common factor times \(\Phi_i(x,y)\).
Factor \(A_i=G\Phi_i(P,Q)\) using coprime polynomials
in \(x,y\). The degree identity gives
\(\deg G+3\max(\deg P,\deg Q)=6-s\) for \(s=r\).
For \(s=2>r\) the highest binary degree is three,
so the right side is three instead. In every case
\(P,Q\) are affine linear with independent linear
parts. Center their common zero at the origin.
For \(s=1\) this gives \(F_z=A=0\) at the origin of
the binary plane.

For \(s=2\), write the map as \(z^2 A+zB+C\).
The lowest nonzero homogeneous part \(A_d\) of \(A\)
is a binary scalar times the independent outer tuple
in centered linear coordinates. The binary
degree-\((2d-2)\) part of the coefficient of \(z^4\)
in the determinant is
\(\det((A_d)_x,(A_d)_y,B(0))=0\), since
derivatives of the constant term of \(B\) vanish.
A nonzero \(B(0)\) would force the image of
\([A_d]\) to lie in a line by projection and
Euler's identity, as in the earlier quadratic-in-\(z\)
case. Thus \(B(0)=0\) and \(F_z\) vanishes over the
binary origin.

For \(s=3\), \(A\) is itself an independent
homogeneous binary cubic tuple. Write the map as
\(z^3A+z^2B+zC+E\). Let \(B_j\) denote the binary
degree-\(j\) part of \(B\). In the \(z^7\)
coefficient, the terms linear in \(B_j\) give
\[
 (6-3j)\det(A_x,A_y,B_j)=0.
\]
Indeed on \(x\ne0\) use \(x,t=y/x\) and write
\(A=x^3\alpha(t)\), \(B_j=x^j b_j(t)\). The three
terms have radial factors \(6,0,-3j\) multiplying
\(\det(\alpha,\alpha',b_j)\), respectively, after
the common power of \(x\) is removed. In particular
\(B_0=0\) by the projection argument, and \(B_1\)
is tangent to the outer cone.

If the projective parametrization is immersive, the
minors of \((A_x,A_y)\) have no common factor by
Euler's identity and the absence of common zeros of
the \(A_i\). Every polynomial tangent vector is then
a polynomial combination of these two columns:
solve over the fraction field and use a minor
invertible at each prime to remove denominators.
Hence \(B_1=0\) in this case.
Otherwise put a nonimmersive parameter at \(x=0\).
Two independent target combinations vanish there with
order at least two, and span \(x^2y,x^3\). By shifting
\(y\) by a multiple of \(x\) and removing multiples of
these two components from the first, normalize
\(A=(y^3,x^2y,x^3)\). The same coprime-minor argument
for \((A_x/x,A_y)\) gives
\(B_1=\lambda(0,2y,3x)\) for a constant \(\lambda\).
In the coefficient of \(z^6\), its binary cubic part
uses only \(A,B_1\). It is
\[
 2\det(A_x,(B_1)_y,B_1)
 +2\det((B_1)_x,A_y,B_1)
 +3\det((B_1)_x,(B_1)_y,A)
 =18\lambda^2 y^3 .
\]
It follows that \(B_1=0\) also here. Now the binary
degree-four part of the same coefficient is
\(\det(A_x,A_y,C(0))=0\). Terms quadratic in \(B\)
have degree at least five because \(B\) starts in
degree two. The projection argument gives \(C(0)=0\),
so \(F_z\) vanishes above the binary origin in this
last case as well.

Finally start with weight-six tuple
\((xz+y^3)\Phi_i(x,y)\) for weights \((1,1,2)\).
If there is any term with \(z\)-exponent \(c>1\),
repeat the threshold argument with
\(w=\min(4-a-b)/(c-1)\).
Now \(w>2\): every such term has weight strictly
below six at \(w=2\). The highest tuple has weight
\(4+w\), is divisible by \(z\), and its quotient
by \(z\) has weight four and restricts to
\(x\Phi_i(x,y)\) on \(z=0\). Factoring through the
same projective outer image gives pencil weight
one and common factor weight one. The quotient
is therefore binary, contradicting the chosen term
with \(c>1\). The map is consequently \(zA+B\).
Its highest binary coefficient tuple in \(A\)
is \(x\Phi_i\). The same projective argument for
\(s=1\) gives \(A=G\Phi(P,Q)\) with \(P,Q\) affine
linear with independent linear parts. The
Jacobian vanishes where \(P=Q=0\), as required.

## Binary quartic common factors on a line pencil

Write the ordinary leading tuple as \(g(x,y)\nu(x,y)\)
and take the first weight threshold \(w>1\) in \(z\),
with \(x,y\) of weight one. The weight-six tuple
has the form \(G\nu(P,Q)\) and
\(\deg_w G+2D=6\), where \(D\) is the positive integer
weight of \(P,Q\). We use the same projective
factorization as in the birational binary case.

If \(D=2\), at least one of \(P,Q\) depends on \(z\):
their binary restrictions are \(xh,yh\) with \(h\)
linear, so if both were binary they would have a
common factor. This requires \(w=2\) and constant
linear \(z\)-coefficients in the pencil pair. Make
\(Q\) binary by a pencil change and take \(Z=P\) as
a triangular source coordinate. Now \(G\) has weight
two, \(P=Z\) has ordinary degree one, and \(Q\) has
ordinary degree two. Only \(GQ^2\) can have ordinary
degree six. Terms of weight below six have ordinary
degree at most five.

If \(D=1\), choose \(P=x,Q=y\); \(G\) then has weight
four. At \(w\ge2\), the only possibilities for the
highest \(z\)-degree \(r\) of \(G\) are \(r=1\) for
\(w=2,3,4\), and \(r=2\) for \(w=2\). Suppose a term
with \(z\)-exponent above \(r\) occurs in the full map,
and increase the weight of \(z\) to the first threshold
at which it reaches the terms of \(z\)-exponent \(r\)
in the leading tuple. Dividing the new leading tuple
by \(z^r\) leaves a tuple of weight \(6-rw\), whose
restriction to \(z=0\) is the nonzero leading
\(z\)-coefficient of \(G\nu(x,y)\). Its projective
image is the same conic. The weight of the quotient
is at most four, and the new \(z\)-weight is strictly
greater than two. In any factorization as a common
factor times a Veronese tuple, the two parameter forms
have weight at most two and the common factor has
weight at most two. Thus the quotient is binary,
contradicting the threshold. The full map has
\(z\)-degree \(r\).

For \(r=1\), its coefficient of \(z\) has degree
at most five and projective image in the conic.
Indeed its highest binary part has dense image
in the conic, and the \(z^2\)-determinant equation
forces projective image of dimension at most one.
Apply the
[conic direction exclusion](SEXTIC_KELLER_CONIC_DIRECTION_EXCLUSION.md).
For \(r=2\), write the map as \(z^2A+zB+C\). The
leading binary coefficient \(A\) has degree two,
with highest part \(\nu(x,y)\) up to a scalar.
Its projective image lies in the conic by the
\(z^5\)-equation. Factoring through the conic and
centering the two affine parameter coordinates makes
\(A=\nu(x,y)\) up to a scalar. The degree-two part
of the coefficient of \(z^4\) is then
\(\det(A_x,A_y,B(0))=0\), so \(B(0)=0\)
by independence and Euler's identity. Thus
\(F_z\) vanishes over the binary origin.

For \(D=1\) and \(w<2\), either \(w=3/2\) with
\(G=g(x,y)+z^2\ell(x,y)\), \(\ell\ne0\), or
\(w=4/3\) with \(G=g(x,y)+cz^3\), \(c\ne0\).
The leading-pair differential after fixing \(t=y/x\)
is \(dx/(x^2z)\) in the first case and \(dx/(xz^2)\)
in the second. These follow by differentiating the
first leading value with respect to \(z\) and
dividing \(x\,dx\,dt\,dz\) by the two target
differentials. At \(w=3/2\), set \(u=x^2z\), giving
\(\ell(1,t)u^2=a_0x-g(1,t)x^7\) and differential
\(dx/u\). At \(w=4/3\), set \(u=xz\), giving
\(cu^3=a_0x-g(1,t)x^7\) and differential \(x\,dx/u^2\).
The polynomial on the right has seven simple roots,
since \(a_0\ne0\) and \(t\) is generic. In each
case the differential is regular at the finite
branch points. At infinity the pole orders of
\((x,u)\) are \((2,7)\) and \((3,7)\), respectively.
The differential is regular also there. Exactness
would force this nonzero differential to vanish,
so both weights are impossible.

It remains that \(D=3\), \(G\) constant. On \(z=0\)
the pair restricts to \((xh,yh)\) for a binary quadratic
\(h\), and \(g\) is proportional to \(h^2\). At \(w=3\)
the \(z\)-coefficients are constant, and a triangular
source coordinate as before reduces the leading rank.
At \(w=2\) with independent linear \(z\)-coefficients,
the differential calculation with \(R=xC-yB\) in the
composite binary proof above applies. The case \(R=0\)
is excluded because it would give
\((P,Q)=(x(z+S),y(z+S))\) with a common factor.
Thus \(R\) is a fourth power and the indicated
triangular shift reduces the leading rank. With
dependent \(z\)-coefficients take \(P=B+xz,\ Q=C\).
The same cubic-cover calculation requires \(x\mid C\)
and a repeated factor of \(C\). If \(x\mid h\),
then \(x\mid B\), so a triangular shift gives \(P=xz\)
and only \(Q^2\) has ordinary degree six. Otherwise
\(C\) is proportional to \(xh\), and its repeated
factor comes from \(h\). The weighted pair is then
equivalent to \(P=xz+y^3,\ Q=xy^2\) by a shift
and rescaling, and the closed-constant-ring exclusion
applies.

At \(w=3/2\), put \(P=B+cz^2,\ Q=C\), with
\(B,C=h\) times independent linear forms. Use the
differential
\[
        \frac{dt}{C(1,t)^{1/6}
                      (p_0C(1,t)-q_0B(1,t))^{1/2}}
\]
from the composite binary calculation. If \(h\) is
squarefree, every zero of these denominator factors
has total exponent at most \(2/6+1/2=5/6\),
because the extra linear factors are distinct for
generic \(p_0,q_0\). This gives an everywhere regular
nonzero differential on the normalized curve and is
impossible. If \(h=x^2\) after a linear change and
the remaining linear factor of \(C\) is independent
of \(x\), take \(P=x^3+cz^2,\ Q=x^2y\) by source
and pencil changes. The generic fibre is now
\(cz^2=p_0-x^3\), and its differential is proportional
to \(dx/(x^2z)\). It is regular except for two
double poles above \(x=0\). If it had a rational
primitive, its anti-invariant part under \(z\mapsto-z\)
would be \(z\,r(x)\) with \(r\) rational in \(x\).
It can have at most a simple pole on each point over
\(x=0\) and is regular elsewhere. At a branch point
where \(z=0\), a pole of \(r\) would give a pole
of \(z\,r(x)\), so \(r\) can have only a simple pole
at \(x=0\). Regularity of \(z\,r(x)\) at infinity
requires \(r\) to vanish there to order at least two,
since \(z\) has pole order three and \(x\) pole order
two on the curve. No such nonzero \(r\) exists.

The sole remaining alternative is \(C=x^3\). Then
\(B=x^2y\) up to a pencil change and rescaling.
Rescale \(z\) and clear the weight denominator to
obtain the exceptional tuple in the statement.

## The first correction to the weighted cube pencil

For \(P=xz+y^3,\ Q=x^3\), the polynomial kernel of
\(D=\det J(P,Q,\cdot)\) is \(k[x,P]\).
Indeed \(y\) is a free parameter on the generic fibre
of \(x,P\), with \(z=(P-y^3)/x\). Since \(Q=x^3\),
the Jacobian derivation has the same rational kernel
as \(\det J(P,x,\cdot)\). Surjectivity of \((x,P)\)
over the algebraic closure gives the polynomial kernel
by the previous contraction argument. At \(x=0\) the
values of \(P=y^3\) fill the whole line.

As before put \(n=(Q^2,-2PQ,P^2)\). The determinant
equation at weight thirteen gives
\(n\cdot F_5\in k[x,P]\), of weight eleven, where
\(F_5\) is the tuple of weight-five parts. All the
weight-eleven monomials in \(x,P\) are attained by
contracting \(n\) with vectors whose entries belong
to \(\langle x^5,x^2P\rangle\). Subtract one such
vector from \(F_5\). The remaining tangent part has
the form \((2PU,QU+PV,2QV)\) with \(U,V\) of weight
two, by the same coprime syzygy as above.

Translations of \(x,y\) change \(U\) by multiples
of \(z,y^2\), respectively. Shifting \(z\) by a
binary linear form changes \(U\) by multiples
of \(x^2,xy\). This removes \(U\). The \(x^2\)
part of \(V\), including the contribution from
translating \(x\) in \(Q\), can be absorbed in the
vector with entries in \(\langle x^5,x^2P\rangle\).
Thus the general normalized weight-five tuple is
\[
   \big(a x^5+b x^2P,\quad
      c x^5+d x^2P+PV,\quad
      e x^5+f x^2P+2QV\big),\qquad
          V=v y^2+wxy+\tau z .
\]
The middle coefficient of \(xz^2\) is \(\tau\), and
none of the changes used to remove \(U\) affects
this coefficient. With \(\tau=0\), a diagonal source
change \(x\mapsto\alpha x,\ y\mapsto\beta y,\
z\mapsto\beta^3 z/\alpha\), followed by compensating
target scalings, can normalize \(v=1\) when \(v\ne0\)
or \(w=1\) when \(v=0,\ w\ne0\). This leaves precisely
three charts for \((v,w)\): \((1,w),(0,1),(0,0)\).

The [exact checker](verify_sextic_weighted_cube_pencil.py)
descends over weights \(j=4,3,2,1\). At weight \(8+j\)
of the determinant, the only terms involving an
unknown part of weight at most \(j\) are the terms
linear in \(F_j\), given by \(2D(n\cdot F_j)\).
Expand this linear map on the complete weight-\(j\)
monomial basis in each coordinate. The known terms
come from all triples of previously determined row
weights whose sum is \(12+j\). Row reduction over
\(\mathbb Q\) expresses each \(F_j\) in terms of
respectively \(8,6,3,3\) free coefficients and a
particular solution. The rows with zero left hand side
give the compatibility equations on the coefficients
of the higher parts. The code carries these equations
as an ideal over \(\mathbb Q\); reduction by its
Gröbner basis changes nothing at a common zero.

At the origin the first two columns of the Jacobian
come from \(F_1\) and the last from \(F_2\). In each
of the three charts its determinant has zero
remainder by the final compatibility ideal. This
contradicts the nonzero Keller constant even if the
remaining lower-weight equations are omitted.
The internal exact replay uses SymPy 1.14.0.

If \(\tau\ne0\), the same diagonal change can
normalize \(\tau=1\). The map has degree at most
two in \(z\): a term \(z^3\) would have weight six,
and it occurs in none of the leading components.
The leading tuple gives \(x^2\) as the first
\(z^2\)-coefficient, and the weight-five tuple gives
\(x\) as its second coefficient. To determine the
weight-four contribution, apply the checker with
\(\tau=1\) just to \(j=4\). Its compatibility
equations give the three \(z^2\)-coefficients
of \(F_4\) as \((0,0,1)\). No weight-five term
contributes \(z^2\) to the first or third component.
Thus the full coefficient of \(z^2\) is
\(A(x)=(x^2,x,1)\). For the linear-in-\(z\)
coefficient \(B\), the specified weight-five part
gives \(y^3\) in \(B_2\), with coefficient one,
and the remaining terms of \(B_2\) have \(y\)-degree
at most two. The \(y\)-degree of \(B_3\) is at
most two by the same displayed tuple. Hence
\(B_2-xB_3\) is cubic in \(y\), and the
[quadratic conic obstruction](SEXTIC_KELLER_CONIC_DIRECTION_EXCLUSION.md#a-moving-conic-frame-for-quadratic-dependence)
excludes this last chart. The preceding reductions
now exclude all three weighted pencils of the
coprime binary branch.

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
  \(q=b(x,y)\) with binary cubics \(b,c\).
  If \(x\nmid q\), a pencil change and an affine source
  change allow \(p=x^2z\). If \(x\mid q\), then \(q\)
  has a repeated factor.

These are conditions on the independent leading parts, with
no restrictions on the lower degree coefficients. Another
condition for unequal component degrees is the
[higher-member basepoint obstruction](KELLER_LEADING_PAIR_RATIONAL_SLICE.md)
when the leading pair is coprime.

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

Finally, if \(x\mid q\) and \(q\) is squarefree, its two
other linear factors are distinct from \(x\). The generic
cubic being integral gives \(p(0,1,0)\ne0\), so \(p,q\)
have independent differentials generically along \(x=0\)
and along either of the other factors of \(q\).
Take \(A=p^2,\ B=pq\). Now the primitive is
\(y/(6x p^2 q)\). Its pole orders along the first and the
chosen second divisor are two and one. The constant field
\(\mathbb C(p,q)\) has valuation of \(q\) equal to one along
each divisor and both dominate \(q=0\). Thus a correction in
this field cannot regularize both pole orders.
This excludes squarefree \(q\) divisible by \(x\).

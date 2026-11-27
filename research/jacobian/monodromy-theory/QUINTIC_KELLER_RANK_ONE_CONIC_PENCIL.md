# Conic and cubic arguments in the rank-one quintic pencil

**Proposition QuinticRankOneConicPencilSquareFactor.** Let \(F\) be
a noninvertible three-variable quintic Keller map in characteristic
zero, with top rank one. If the primitive pencil in its cancelling
quartic span has degree two, then the fixed factor is the square of
a line. Over the algebraic closure we may write its leading parts
as
\[
              Q_4=z^4,\qquad R_4=z^2q,\qquad
              P_5=z H_2(z^2,q),
\]
where \(q\) is an irreducible ternary conic and \(H_2\) is a
nonzero binary quadratic.

Use \(A=P_5,\ C=Q_4,\ D=R_4\), and the notation of the
[rank-one pencil reduction](QUINTIC_KELLER_REDUCED_BOUNDARY_PENCIL_RESIDUAL.md#rank-one-pencil).
This gives a primitive conic pencil \(\langle p,q\rangle\) with
generally irreducible members, \(\deg g+2e=4,\ e\ge1\),
and the valuation identity
\[
                  4\nu_f(A)-5\nu_f(g)=\rho_M\nu_f(M)
\tag{1}
\]
for every factor \(f\) of a member \(M\). All factors of \(g\) and
\(A\) occur in pencil members, and distinct members have no common
factor. We work over \(\mathbb C\), with the same finite coefficient
field embedding giving the characteristic-zero statement.

First suppose \(e=2\), so \(g\) is constant. A general cancelling
quartic is a reduced product of two distinct irreducible conics
in the pencil. The [required boundary point](QUINTIC_KELLER_REDUCED_BOUNDARY_PENCIL_RESIDUAL.md#a-reduced-quartic-cannot-use-a-linear-source-residual)
gives intersection multiplicity at least three at one point,
and the [contact-four exclusion](KELLER_QUARTIC_CONTACT_FOUR_CONIC_EXCLUSION.md)
therefore forces intersection multiplicities \(3,1\).

Next suppose \(e=1\), so \(\deg g=2\). If \(g\) is reduced,
a general cancelling quartic is reduced. All of its factors have
degree at most two, so the same boundary arguments force it to be
a product of two irreducible conics with intersection multiplicities
\(3,1\). In this case \(g\) itself is one of the conics, and (1)
puts it in the primitive pencil because it has degree two.
Thus in either of the reduced cases there are two smooth members
of the primitive pencil with intersection multiplicities \(3,1\).

If the pencil contained a doubled line, restricting that member
to either of those smooth conics would give an even intersection
divisor. But the restrictions of all distinct pencil members to a
fixed member give the same divisor. The multiplicities \(3,1\)
therefore rule out every doubled line. All conic members of the
pencil are reduced. On each member (1) then makes the multiplicities
of \(A\) equal on its factors, apart from the possible member \(g\)
when \(e=1\). That member is itself irreducible and of degree two.
Every member contributes an even degree to \(A\), contrary to
\(\deg A=5\).

It follows that \(e=1\) and \(g=z^2\) up to scalar in a linear frame.
The member containing \(z\) must be \(z^2\): a member \(zL\) with
\(L\) a distinct line would make (1) give
\(4(\nu_z(A)-\nu_L(A))=10\).
There is no second doubled line, since two distinct squares would
make the pencil binary and its general conic reducible.
Choose an irreducible member \(q\) and take the cancelling span as
\(\langle z^4,z^2q\rangle\). On all members other than \(z^2\),
the multiplicities in \(A\) are constant along each reduced
member by (1). Their total contribution to \(\deg A\) is even,
so \(\nu_z(A)\) is one, three, or five. Thus after taking out
one \(z\), \(A/z\) is a product of two conic pencil members,
with repetition allowed. This is the asserted form.

## The pure quartic member of the cubic-argument pencil

**Proposition QuinticRankOneCubicPencilAffineLinearMember.** In the
degree-three primitive pencil case of a noninvertible top-rank-one
quintic Keller map, choose the
[leading frame](QUINTIC_KELLER_REDUCED_BOUNDARY_PENCIL_RESIDUAL.md#degree-three)
\(Q_4=z^4,\ R_4=zT\) and take \(x\) so that \(T(x,y,0)\) is a
nonzero multiple of \(x^3\). The member \(Q\) of pure quartic
leading form is affine linear in \(y\). It either depends only on
\(x,z\), or has, after translating \(z\), the form
\[
                  Q=\mu+\lambda z+z^2(z^2+c y+d x+e z+f),
                       \qquad \lambda c\ne0 .
\]

For a general scalar \(s\), the leading curve of \(R+sQ\) is
the union of the line \(z=0\) and the irreducible cubic \(T+sz^3\).
Its unique intersection with the line is \(p=[0:1:0]\).
All singularities away from \(p\) are ordinary nodes or cusps
of the cubic. The
[required double point](QUINTIC_KELLER_REDUCED_BOUNDARY_PENCIL_RESIDUAL.md#a-reduced-quartic-cannot-use-a-linear-source-residual)
must be \(p\); hence the cubic is smooth there. Its tangent is
\(z=0\) and the \(y^2z\) coefficient of \(T\) is some fixed
\(\eta\ne0\). The full projective fibre has vanishing linear jet
and quadratic rank one there for general \(s\).

In the chart \(y=1\), write \(u=x/y,\ v=z/y,\ w\) the normal
coordinate. The quadratic jet restricts to \(\eta v^2\)
on \(w=0\); it must be a square proportional to \((v+\beta(s)w)^2\).
In particular there is no \(uw\) coefficient, and vanishing
of the linear jet excludes a \(y^3\) term of \(R+sQ\).
The coefficient of \(y^2\) in \(R+sQ\) therefore has the form
\[
             \eta z^2+(\alpha+s\alpha')z+(\delta+s\delta')
\]
with discriminant zero for general \(s\). Its \(s^2\) and
\(s\) coefficients force \(\alpha'=\delta'=0\). Since \(Q_4=z^4\),
this proves \(Q\) is affine linear in \(y\).

Write \(Q=a(x,z)y+b(x,z)\), with \(\deg a\le2\).
The [affine-linear coefficient classification](KELLER_LINEAR_SOURCE_COEFFICIENT_CUBIC_CLASSIFICATION.md)
shows that if \(a\ne0\), noninvertibility requires a doubled
affine line in \(a\), and \(b\) constant on that line with
nonzero constant normal derivative. Thus for an affine linear
polynomial \(\ell\),
\(Q-\mu=\lambda\ell+\ell^2 H\), \(\lambda\ne0\).
The leading part of \(\ell\) divides \(z^4\), so it is
proportional to \(z\). Translate and rescale its affine equation
to \(z\). Then \(H_2=z^2\) because \(Q_4=z^4\), and the
coefficient of \(y\) in \(H\) is a nonzero constant because
\(a\ne0\). This gives the stated form.

## Excluding the cubic argument

**Theorem QuinticRankOneCubicPencilExclusion.** A noninvertible
three-variable quintic Keller map of top rank one in characteristic
zero cannot have primitive cancelling-pencil degree three.

Work over \(\mathbb C\). In the preceding coordinates,
\[
            P_5=z^2(\alpha T+\beta z^3),\qquad
            Q_4=z^4,\qquad R_4=zT ,
\tag{2}
\]
where \(T(x,y,0)\) is a nonzero multiple of \(x^3\).
The cubic is smooth at \(p=[0:1:0]\), so its \(y^2z\)
coefficient is nonzero. In particular, the coefficient of
\(y^2\) in \(R\) is \(\eta (z-z_*)^2\), with \(\eta\ne0\),
and \(R\) is quadratic in \(y\).

Suppose first that \(Q\) depends on \(y\). By the preceding
proposition, an affine change in \(y\) and a translation of \(z\)
put it in the form \(Q=z^4+z^2y+\lambda z\), after subtracting
a constant. These changes preserve (2) and the fact that
\(T(x,y,0)\) is a nonzero multiple of \(x^3\).
At \(Q=t\) on \(z\ne0\) substitute
\(y=t/z^2-\lambda/z-z^2\), denoting the resulting functions
by \(P_t,R_t\) in the ring \(\mathbb C(t)[z,z^{-1},x]\).
Their binary Jacobian bracket in \(x,z\) is a nonzero
constant times \(z^{-2}\), by \(Q_y=z^2\).

The degree of \(R_t\) in \(x\) is three. Its leading
\(x\)-coefficient \(r_3\) is an affine linear polynomial in \(z\)
with nonzero slope, from \([x^3](zT)\) and a possible constant
from \(R_3\). There is no \(x^3y\) term in \(R_4\), since
\(R_4\) is divisible by \(z\), and all its other terms involving
\(y\) have \(x\)-degree at most two. The \(x\)-degree of \(P_t\)
is at most four. If it is zero, the bracket cannot be a nonzero
function independent of \(x\), since \(R_{t,x}\) has degree two
in \(x\). If it is \(j>0\) with leading coefficient \(p_j\),
the coefficient of \(x^{j+2}\) in the bracket gives
\[
                     j p_j r_3'-3r_3 p_j'=0.
\]
Thus \(p_j^3/r_3^j\) is independent of \(z\). Taking orders at
the simple zero of \(r_3\) shows \(3\mid j\).
We must have \(j=3\) and \(p_3=c(t)r_3(z)\).

The only monomials of \(P\) containing \(x^3\) have coefficients
of the form \(\alpha' z^2+b z+c_0+d y\), by (2) and the degree
bound. Here \(\alpha'=0\) if and only if \(\alpha=0\).
Thus
\[
                  p_3=\alpha' z^2+b z+c_0+
                                     d(t/z^2-\lambda/z-z^2).
\]
Since \(p_3=c(t)r_3(z)\) has no pole, \(d=0\), and then
\(\alpha'=0\) and \(c(t)\) is constant. Subtract that constant
multiple of \(R\) from \(P\). The bracket on the fibre is
unchanged, but \(P_t\) now has degree at most two in \(x\).
This contradicts the preceding degree argument.

It remains to consider \(Q=Q(x,z)\). It depends on \(x\),
since a univariate quartic \(Q(z)\) has a critical point.
Moreover it is non-composite as a polynomial in \(x,z\):
a proper composition of polynomials would have a critical point
over \(\mathbb C\), because a nonconstant inner polynomial
takes every complex value. The
[closed-polynomial criterion of Wang and Yan](https://arxiv.org/pdf/2610.36495v1)
therefore gives \(\mathbb C(Q)\) as the rational constants of
the binary Jacobian derivation with \(Q\).

Write \(u=z-z_*\) and \(R=\eta u^2y^2+A(x,z)y+B(x,z)\).
If \(P\) were independent of \(y\), the Jacobian would be
a product of \(R_y\) and a polynomial in \(x,z\), impossible
for a nonzero scalar. Thus \(m=\deg_y P>0\). If \(p_m(x,z)\)
is its leading \(y\)-coefficient, the highest equation in
\(y\) reads
\[
                        \{p_m/u^m,Q\}_{x,z}=0.
\]
Consequently \(p_m=u^m h(Q)\) with \(h\) rational.
In fact \(h\) is polynomial. Write it as a ratio of coprime
univariate polynomials; substitution preserves the Bézout
identity, so its denominator after evaluation at \(Q\) must
divide \(u^m\). Any nonconstant polynomial in \(Q\) depends
on \(x\), whereas a divisor of \(u^m\) depends only on \(z\).
It follows that \(h\) is polynomial, and nonzero constant
since \(m+4\deg h\le5-m\) and \(m\ge1\).
In particular \(m\le2\). If \(m=2\), subtract a constant
multiple of \(R\) from \(P\) and reapply this argument.
The case of \(y\)-degree zero is still impossible, so we
obtain \(P=c_1 u y+D(x,z)\) with \(c_1\ne0\).
These operations preserve the degree-five leading part.
Its independence of \(y\) forces \(\alpha=0\) in (2), since
\(T\) has a nonzero \(y^2z\) term. Thus \(D\) has degree five.

Now the coefficient of \(y\) in the determinant is
\[
 c_1\big(u\{A,Q\}_{x,z}-A\{u,Q\}_{x,z}\big)
                           -2\eta u^2\{D,Q\}_{x,z}=0 .
\]
After division by \(u^2\) this gives
\(c_1 A/u-2\eta D=H(Q)\) for a rational \(H\).
The same denominator argument, now with \(u\) in place of
\(u^m\), makes \(H\) polynomial. But \(\deg A\le3\),
\(\deg (uD)=6\), and a nonzero \(uH(Q)\) has degree
\(1+4\deg H\), which cannot be six. The equality
\(c_1A=2\eta uD+uH(Q)\) is impossible.
This finishes the cubic-argument exclusion.

## Consequence for the quintic inversion problem

**Corollary QuinticKellerNoninvertibleQuarticLeadingNonreduced.**
Over a characteristic-zero field, a noninvertible degree-five
three-variable Keller map has dependent degree-five parts. Every
nonzero target-linear combination canceling degree five has degree four
and geometrically nonreduced leading form. In top rank one the primitive
cancelling-pencil degree is one, and both cancelling cubic terms
vanish along the fixed repeated line.

The [independent-leading exclusion](QUINTIC_KELLER_TERNARY_PENCIL_TRANSVERSE_PARAMETER.md#independent-degree-five-leading-parts)
and the [cubic component theorem](QUINTIC_KELLER_CUBIC_COMPONENT.md)
give the degree and dependence assertions. In rank two the
[reduced-quartic inversion theorem](KELLER_QUARTIC_CONTACT_FOUR_CONIC_EXCLUSION.md#reduced-leading-forms-in-top-rank-two)
forces the unique cancelling leading quartic to be nonreduced. In
rank one the [primitive-pencil list](QUINTIC_KELLER_REDUCED_BOUNDARY_PENCIL_RESIDUAL.md#rank-one-pencil)
has degrees one through four. The
[double-line conic exclusion](QUINTIC_KELLER_RANK_ONE_DOUBLE_LINE_CONIC_PENCIL.md),
the cubic-pencil exclusion above and the
[quartic-pencil exclusion](QUINTIC_KELLER_RANK_ONE_QUARTIC_PENCIL_CLOSURE.md)
remove degrees two through four. In degree one every cancelling
quartic contains the nonreduced fixed factor, and the
[line-pencil cubic restriction](QUINTIC_KELLER_RANK_ONE_LINE_PENCIL_CUBIC_RESTRICTION.md)
gives the vanishing of both cubics along its repeated line.

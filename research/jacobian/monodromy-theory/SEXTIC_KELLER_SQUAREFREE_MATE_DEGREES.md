# Source degrees for the sextic squarefree quartic form

Work with the [squarefree quartic form](QUARTIC_AFFINE_SOURCE_COMPONENTS_THROUGH_DEGREE_SEVEN.md)
\(R=az+b\), where \(a=xv\), \(v=xy-1\), and
\[
 b=\lambda_0y+\mu_0+
       xy\bigl((\lambda_1-\lambda_0)y+\mu_1-\mu_0\bigr),
 \qquad \lambda_0\lambda_1\ne0 .
\]
In particular \(\deg b\le3\).

**Theorem KellerSquarefreeQuarticConstantLeadingMateExclusion.**
Over a characteristic-zero field, a polynomial Keller completion
of \(R\) has no mate of positive \(z\)-degree with constant leading
\(z\)-coefficient. This holds without an ordinary degree bound
and for arbitrary \(b\in k[x,y]\).

**Theorem KellerSquarefreeQuarticAffinePowerLeadingExclusion.**
For the same \(a\) and arbitrary \(b\), there is no polynomial
Keller mate of \(z\)-degree \(m>0\) with leading
\(z\)-coefficient a nonzero scalar times \(L^m\), where \(L\)
is affine linear and coprime to \(a\).

**Theorem SexticKellerSquarefreeQuarticMateDegreeBound.**
If \(P,Q,R\) form a characteristic-zero polynomial Keller map
with \(\deg P,\deg Q\le6\), each of \(P,Q\) has \(z\)-degree at
most three. In conjunction with the
[low-degree mate restrictions](KELLER_SQUAREFREE_QUARTIC_FIBRES_AND_MATES.md),
both \(z\)-degrees are positive and at least one is two or three.

**Theorem SexticKellerSquarefreeQuarticCubicLeadingRestriction.**
Under the same ordinary degree bound, the leading \(z\)-coefficient
of any mate of \(z\)-degree three is either \(\alpha a+\beta\)
with \(\alpha\ne0\), or a nonzero scalar times \(x^3\).
In the first case, the other mate can be given \(z\)-degree two
with leading \(z\)-coefficient \(\gamma a+\delta\), \(\gamma\ne0\),
by a constant target change fixing \(R\). In the second case,
after such a change the other mate either has \(z\)-degree two
and leading coefficient \(x^2(\gamma v+\delta)\), or has
\(z\)-degree one and leading coefficient a nonzero multiple of \(x\).

All arguments can be made over \(\mathbb C\) by embedding the
finitely generated field of coefficients. Write \(\{f,g\}\) for
\(f_xg_y-f_yg_x\) with the indicated third coordinate held fixed.

## Constant leading coefficient

**Lemma SquarefreeCubicCoefficientHasNoPlanarJacobianMate.**
For \(a=x(xy-1)\), the polynomial kernel of \(\{a,\cdot\}\)
is \(\mathbb C[a]\), and no polynomial \(h\) satisfies
\(\{a,h\}=c\) for \(c\in\mathbb C^\times\).

At fixed \(a\), the derivation \(\{a,\cdot\}\) is
\(-x^2\partial_x\), since \(y=(a+x)/x^2\). Its rational kernel
is \(\mathbb C(a)\). A rational function in \(a\) with a finite
pole has a pole on a component of the nonempty fibre at that
value. Thus the polynomial kernel is \(\mathbb C[a]\).
If \(\{a,h\}=c\), then \(h=c/x+f(a)\). A pole of \(f\) at
\(a=0\) cannot cancel along \(v=0\), where \(c/x\) is regular.
With \(f\) regular at \(a=0\), the simple pole at \(x=0\)
cannot cancel.

After scaling, write one mate \(P=z^m+\cdots\), \(m>0\).
In the formal Laurent series ring over \(\mathbb C[x,y]\),
use \(w=P^{1/m}\) as coordinate at \(z=\infty\). The binomial
expansion gives \(w=z+O(1)\), and formal inversion near infinity
gives \(z=w+O(1)\) with polynomial coefficients, by successively
solving for each lower coefficient with leading linear multiplier
one. Put dots over polynomials expressed in \(x,y,w\). The
determinant condition in these coordinates is
\[
               \{\dot R,\dot Q\}=\frac{c}{\dot{P_z}}
                    =\frac{c}{m}w^{1-m}+O(w^{-m}),
 \qquad c\in\mathbb C^\times.
\]
Also \(\dot R=aw+O(1)\). If the leading term of \(\dot Q\)
is \(q_k(x,y)w^k\) with \(k>-m\), comparison at \(w^{k+1}\)
gives \(\{a,q_k\}=0\). For \(q_k=\sum_j d_j a^j\), subtraction
of \(\sum_j d_j w^{k-j}\dot R^j\) lowers the Laurent degree and
preserves the bracket. Every coefficient after subtraction is
still polynomial in \(x,y\). Descending finitely to degree
\(-m\) gives a coefficient \(h\in\mathbb C[x,y]\) satisfying
\(\{a,h\}=c/m\), contrary to the lemma.

## An affine power as leading coefficient

Scale \(P\) to have leading coefficient \(L^m\). For a
nonconstant \(L=\ell_x x+\ell_y y+\ell_0\) coprime to \(a\),
the root coordinate \(w=P^{1/m}\) is now \(Lz+O(1)\).
Inversion gives Laurent coefficients in \(\mathbb C[x,y,L^{-1}]\),
and \(\dot R=(a/L)w+O(1)\), \(\dot{P_z}=mLw^{m-1}+O(w^{m-2})\).
Put \(u=a/L\). Solving \(a=uL\) for \(y\) gives
\[
               y=\frac{(1+u\ell_x)x+u\ell_0}{x^2-u\ell_y}.
\]
Thus \(\mathbb C(x,y)=\mathbb C(u,x)\) and the rational kernel
of \(\{u,\cdot\}\) is \(\mathbb C(u)\). Each finite fibre of
\(u\) meets \(L\ne0\) in a curve since \(a-cL\) is nonconstant
and coprime to \(L\). The kernel in \(\mathbb C[x,y,L^{-1}]\)
is therefore \(\mathbb C[u]\).

The same finite descent using \(w^{k-j}\dot R^j\) now makes
\(\{u,h\}=c/(mL)\) for some \(h\in\mathbb C[x,y,L^{-1}]\).
At fixed \(u\),
\[
                   \partial_x h=-\frac{c/m}{x^2-u\ell_y}.
\]
If \(\ell_y\ne0\), the right side has simple poles over
\(\overline{\mathbb C(u)}\) and cannot be a rational derivative.
If \(\ell_y=0\), then \(h=c/(mx)+f(u)\). Along each of \(v=0\)
and \(x=0\), the generic point lies in \(L\ne0\) and \(u\)
has a simple zero. A pole of \(f\) at zero cannot be cancelled
along \(v=0\), and without that pole \(h\) cannot be regular
along \(x=0\). This proves the affine-power exclusion.

## Leading coefficients in the \(R\)-fibre chart

Put \(r=R\), and use \(z=(r-b)/a\). The two mates become
polynomials in \(r\) over \(\mathbb C[x,y,a^{-1}]\), of the
same degrees \(m,n\) as their original \(z\)-degrees. If \(p,q\)
are their leading \(z\)-coefficients, their leading
\(r\)-coefficients are \(p/a^m,q/a^n\). The determinant of
the change of coordinates is \(a\), so at fixed \(r\) the
bracket of the transformed mates is \(\kappa/a\),
\(\kappa\in\mathbb C^\times\). For \(m,n>0\),
\[
                        \{p/a^m,q/a^n\}=0.                 \tag{1}
\]

Call a nonconstant rational function \(u\) closed if \(\mathbb C(u)\)
is relatively algebraically closed in \(\mathbb C(x,y)\).
The kernel of \(\{u,\cdot\}\) is then \(\mathbb C(u)\) by the
characteristic-zero differential criterion for algebraic dependence.
We use [Petravchuk–Iena, Lemma 2 and Theorem 1](sources/petravchuk-iena-math-0701588v2.pdf).
A nonclosed rational function is a composition of a rational
function in \(x,y\) with a rational map of degree at least two on
the projective line. A quotient of coprime, algebraically independent
polynomials is closed when either polynomial is irreducible.

Suppose \(u=p/a^m\) is closed, has a pole of order exactly \(m\)
along one divisor of \(a=0\), and every finite fibre of \(u\)
meets \(a\ne0\) in a curve. Then (1) gives \(q/a^n=f(u)\).
Since \(q/a^n\) is regular on \(a\ne0\), \(f\) is polynomial.
The pole bounds give \(m\deg f\le n\). If \(f=\sum_j c_j T^j\),
subtract \(\sum_j c_j R^{n-mj}P^j\) from \(Q\).
This lowers its \(z\)-degree, preserves the determinant and
keeps it polynomial. Iteration reaches a \(z\)-independent
mate, which is impossible by the no-planar-mate theorem.
The ordinary degree can increase in these subtractions because
the no-planar-mate theorem has no degree bound.

For the closed quotient cases below, the finite fibre assertion
follows in the reduced expression \(u=L/D\), where \(D\) is
supported on \(xv=0\). If \(L\) is nonconstant,
\(L-cD\) is nonconstant for all \(c\) and has no factor in common
with \(D\), since \(\deg D>\deg L\) and \(\gcd(L,D)=1\).
In the applications here both \(x\) and \(v\) divide \(D\), so
the entire finite fibre has no irreducible component supported
on \(a=0\).

The other quotients needed below are
\(u=L/a^d\), with \(L=x\) or \(v\), and \(d\ge2\).
They are closed since \(\mathbb C(x,y)=\mathbb C(u,a)\).
For \(L=x\), solve for \(x=u a^d\), and for \(L=v\), solve
for \(x=1/(u a^{d-1})\). The finite nonzero fibres meet
\(a\ne0\) in curves, by the same reduced-numerator argument
with constant numerator. If a polynomial \(h(x,y)\) satisfies
\(h/a^k=f(u)\) for some \(f\in\mathbb C(T)\) and \(k>0\),
then \(f\) is Laurent polynomial. Its largest exponent is at
most \(\lfloor k/d\rfloor\), since \(u\) has an order-\(d\)
pole on the divisor of \(a=0\) other than \(L=0\).
Thus the possible \(h\) are linear combinations of
\[
          a^{\,k-dj} L^j,\qquad j\in\mathbb Z,\quad dj\le k.
                                                               \tag{2}
\]
For negative \(j\) these are polynomials because \(a^d/L\)
is polynomial. Their distinct ordinary degrees are
\(3k+j(\deg L-3d)\).

## Classification at high \(z\)-degree

Choose \(m=\deg_z P\ge n=\deg_z Q\), and suppose \(m\ge4\).
The constant-leading exclusion removes \(m=6\) and all cases
with constant \(p\). By the no-planar-mate theorem \(n>0\).
For \(m=4,5\) the nonconstant \(p\) has degree at most \(6-m\).

First suppose \(\gcd(p,a)=1\). The polynomials \(p,a\) are
algebraically independent. Otherwise \(p\) belongs to
\(\mathbb C(a)\) since \(\mathbb C(x,y)=\mathbb C(a)(x)\),
and then belongs to \(\mathbb C[a]\), contrary to its
positive degree at most two. If \(p/a^m\) is closed, the
preceding fibre-chart subtraction excludes it.

For a nonclosed \(p/a^m\), represent it as a quotient of
binary forms of common degree \(D\ge2\) evaluated on a
relatively prime polynomial pair \(F,G\). The binary forms
have no common root, so their values are coprime. Their
linear factors evaluate to members of the pencil of \(F,G\).
No nonzero member is constant, since otherwise \(p,a^m\)
would be polynomials in a single pencil generator. The
numerator degree gives \(D=2\), \(\deg p=2\) and \(m=4\).
If the denominator binary form had two distinct roots,
the two distinct polynomial members would be coprime and
nonconstant, and their product would be a constant times
\(x^4v^4\). They would be proportional to \(x^4\) and \(v^4\)
in one order. Every nonzero member of their pencil has
degree at least four, contradicting the degree-two product
in the numerator. Thus the denominator is a square of a
member proportional to \(a^2\). The numerator's two members
are affine linear, and they must be proportional since
otherwise their pencil could not contain \(a^2\).
Consequently \(p\) is proportional to \(L^2\) for an affine
linear \(L\) coprime to \(a\).

If \(x\mid p\), the alternatives are a multiple of \(x\),
a multiple of \(x^2\), or \(xL\) with \(L\) nonconstant
affine linear coprime to \(a\). In the last case
\(p/a^m=L/(x^{m-1}v^m)\) up to a constant. This reduced
quotient is closed by the irreducible numerator criterion.
Its numerator and denominator are algebraically independent:
the denominator would otherwise be in \(\mathbb C[L]\)
by a linear coordinate change, so it could have only linear
irreducible factors. The pole of order \(m\) along \(v=0\)
again gives the fibre-chart subtraction. If \(v\mid p\),
the bound on \(\deg p\) makes \(p\) a multiple of \(v\).

For \(p\) proportional to \(x\) or \(v\), apply (2) with
\(d=m\). Negative \(j\) exceed \(\deg q\le 6-n\), since
their degrees are at least \(3n+3m-\deg p\). If \(n=m\),
the only term allowed in \(q\) is proportional to \(p\);
subtract that multiple of \(P\) from \(Q\) and reduce \(n\).
If \(n<m\), only the \(j=0\) term remains, which forces
\(n=1\), \(q\) proportional to \(a\). Subtracting the same
multiple of \(R\) from \(Q\) gives a planar mate.

The remaining possibility is \(m=4\), \(p\) proportional
to \(L^2\) with \(L\) affine linear and either coprime to
\(a\) or proportional to \(x\). Now \(L/a^2\) is closed.
This follows from the irreducible numerator criterion when
\(\gcd(L,a)=1\), and from the rational coordinates above
when \(L=x\). Thus the kernel of
\(\{p/a^4,\cdot\}\) is \(\mathbb C(L/a^2)\). Formula (2)
applies with \(d=2\), with only nonnegative \(j\) in the
case \(\gcd(L,a)=1\). In the case \(L=x\), negative \(j\)
all exceed the degree bound, having degree at least
\(3n+5>6-n\). The allowed leading coefficients for
\(n=4,3,2,1\) are respectively a multiple of \(L^2\),
none, a multiple of \(L\), and a multiple of \(a\).
After the constant linear reductions at \(n=4,1\),
only the pair \((m,n)=(4,2)\) remains, with \(q\)
proportional to \(L\). For \(\gcd(L,a)=1\), apply the
fibre-chart subtraction with the degree-two mate as
the fixed member, using \(u=L/a^2\) and its order-two
pole. It remains to treat \(L=x\).

## Final centering for \(p=x^2,\ q=x\)

By rescaling the mates, write
\[
 Q=xz^2+c z+d,\qquad P=x^2z^4+e z^3+f z^2+g z+h
\]
with \(\deg c\le5,\deg d\le6,\deg e\le3,\deg f\le4,\deg g\le5\).
In the following subtractions the determinant remains constant.
Reduce the \(z\)-degree to at most three by subtracting \(Q^2\)
from \(P\). The \(z^3\)-coefficient \(e-2xc\) has degree at
most six. By (1) and (2) with \(k=3,d=2,L=x\), that
coefficient equals \(\tau xa\), allowing zero as well.
Subtract \(\tau RQ\). In particular \(2c+\tau a=e/x\)
is a polynomial of degree at most two, and \(\deg c\le3\).

The \(z^2\)-coefficient after these subtractions is
\[
                  f-c^2-2xd-\tau(ac+bx).
\]
Its degree is at most seven. Formula (2) with \(k=2\)
makes it \(\beta a^2+\gamma x\). Subtract
\(\beta R^2+\gamma Q\). The coefficient equation gives
\(\deg(2xd)\le6\) and therefore \(\deg d\le5\).

The \(z\)-coefficient after the last subtraction is
\[
             g-(2c+\tau a)d-\tau bc-2\beta ab-\gamma c .
\]
Its degree is at most seven. Formula (2) with \(k=1\)
allows only a multiple of \(a\), since \(j=-1\) already
has degree eight. Subtracting the corresponding multiple
of \(R\) now gives a planar mate. When a coefficient
vanishes in these steps, the same formula applies with
all its scalar coefficients zero. This excludes \(m=4\)
and completes the proof.

## The cubic leading coefficients

Suppose \(m=3\) and \(\deg p\le3\). If \(p,a\) are algebraically
dependent, the same kernel argument as above gives
\(p\in\mathbb C[a]\), hence \(p=\alpha a+\beta\) with
\(\alpha\ne0\) by the constant-leading exclusion.

Otherwise consider first \(\gcd(p,a)=1\). If \(p/a^3\) is
closed, the polynomial target subtractions in the \(R\)-fibre
chart exclude it. For a nonclosed quotient, the binary forms
in the pencil argument have common degree \(D=2\) or \(3\).
If there are two distinct denominator members, they must be
proportional to pure powers of \(x\) and \(v\). When \(D=2\)
these are \(x^3,v^3\), so the two numerator members have product
of degree at least six. When \(D=3\), one denominator member
would occur twice, which cannot produce the odd multiplicities
in \(a^3\). Three distinct denominator members are impossible
since \(a^3\) has only two prime factors. Hence the denominator
is a pure power. This forces \(D=3\) and the denominator member
to be proportional to \(a\). All three numerator members must
be affine linear, hence proportional to one \(L\) as their
pencil contains \(a\). This makes \(p\) proportional to \(L^3\),
excluded by the affine-power theorem.

For \(v\mid p\), write \(p=vL\) up to a scalar, where \(L\)
is constant or affine linear. If \(L\) is proportional to \(x\)
then \(p\) is proportional to \(a\). If \(L\) is nonconstant
and not proportional to \(x\), the quotient
\(L/(x^3v^2)\) is closed by the irreducible numerator criterion.
The polynomials are independent since the denominator has the
nonlinear factor \(v\), whereas a polynomial in \(L\) has only
linear factors. The quotient has a pole of order three along
\(x=0\), so the polynomial target subtractions apply.

If \(x\mid p\) but \(v\nmid p\), take the maximal power of
\(x\) out of \(p\). For \(p=xL\) with \(\deg L>0\) and
\(\gcd(L,a)=1\), the quotient is \(L/(x^2v^3)\).
Its reduced numerator and denominator are independent. If the
denominator were algebraic over a nonconstant linear \(L\) it
would be a polynomial in \(L\), impossible since \(v\) is
nonlinear. For quadratic \(L\), independence follows too:
if the two polynomials were dependent, they would be
polynomials in a common closed polynomial generator, by
the polynomial form of the closed-generator theorem of
[Wang–Yan, Proposition 2.2](sources/wang-yan-2610.36495v1.pdf).
Its degree divides both two and eight, so is one or two.
For degree one all factors of the denominator would be linear.
For degree two the fibres supporting the factors \(x\) and \(v\)
would either be the same, requiring that generator to have a
fibre of degree at least three, or be distinct. In distinct
fibres the one supported on \(v=0\) is proportional to \(v\),
and the one supported on \(x=0\) is proportional to \(x^2\).
These two polynomials do not differ by a constant after
scaling. Thus independence holds. A nonclosed composition
would have binary degree two. Its denominator cannot be a
square, and two distinct denominator members are proportional
to \(x^2,v^3\). This would force numerator degree at least four,
a contradiction. So the quotient is closed and the order-three
pole on \(v=0\) gives the same exclusion.
For \(p=x^2 L\) with \(L\) nonconstant affine linear and
not proportional to \(x\), the reduced quotient \(L/(xv^3)\)
is closed by the irreducible numerator criterion, with
independence by the nonlinear factor \(v\). This again gives
the order-three pole on \(v=0\). The case \(p\) proportional
to \(x^3\) is retained.

For the three unit quotients with \(p\) proportional to
\(v,x\) or \(x^2\), all are closed. For \(p\) proportional
to \(v\) or \(x\) use the rational coordinates already given.
For \(p=x^2\), \(u=x^2/a^3=1/(xv^3)\) gives
\(\mathbb C(x,y)=\mathbb C(u,v)\). Each quotient has pole
order three along one divisor of \(a=0\) and every nonzero
finite fibre meets \(a\ne0\) in a curve. Apply (2) with
\(d=3,L=p\), after scaling \(p\) as needed. Negative \(j\)
have degree at least \(3n+7>6-n\). At \(n=3\) only the
term proportional to \(p\) can occur, so subtract a multiple
of \(P\) and reduce \(n\). At \(n<3\) the constant term
forces \(n=1\), \(q\) proportional to \(a\), which reduces
to a planar mate and is excluded.

Finally, for \(p=\alpha a+\beta\), equation (1) implies
\(q\in\mathbb C[a]\) and its degree bound leaves \(q\) affine
in \(a\). At \(n=3\) independent leading coefficients would
give a constant nonzero leading coefficient after a linear
target change. Proportional leading coefficients can be
cancelled. At \(n=1\), subtract a multiple of \(R\) to
obtain a planar or constant-leading mate, both excluded.
Thus \(n=2\) and \(q=\gamma a+\delta\) with \(\gamma\ne0\).
For \(p\) proportional to \(x^3\), the rational kernel of
\(\{p/a^3,\cdot\}\) is \(\mathbb C(1/v)\). Equation (2)
with \(d=1,L=x\) describes \(q\), since \(x/a=1/v\).
The terms have degrees \(3n-2j\), \(j\le n\). At \(n=3\)
only \(j=3\) survives, allowing cancellation of the leading
coefficient. At \(n=2\) one obtains \(x^2(\gamma v+\delta)\);
at \(n=1\) one obtains a linear combination of \(a,x\), and
can remove the multiple of \(a\) using \(R\).

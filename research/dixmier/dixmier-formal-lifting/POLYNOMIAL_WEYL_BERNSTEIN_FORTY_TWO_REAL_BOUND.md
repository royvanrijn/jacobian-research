# The degree-forty-two residual and the bound over fields without a square root of minus fifteen

**Theorem PolynomialWeylBernsteinFortyTwoResidual.** Let \(k\) have
characteristic zero. If a nongenerating polynomial Weyl pair with
nonzero scalar commutator has a member of Bernstein degree \(42\),
it has a polynomial presentation over \(k\) with
unordered Bernstein degrees \(\{42,147\}\). In an ordering with
smaller degree first, the crossing edge in Fourier coordinates
\(X,Y\) has direction \((3,-1)\) and a common factor
\[
             R=X L(z),\qquad z=XY^3,\qquad
             L(z)\ \text{proportional to}\ (z-u)^4(z-v),
                    \qquad u v\ne0,\quad v\ne u .
\]
The leading symbols are proportional to \(R^2,R^7\).
Moreover \(u,v\in k\) and \(2(v/u)^2-5(v/u)+5=0\). In particular
\(-15\) is a square in \(k\).
After a fractional shift at a fourth-order root of the common
face polynomial, the next edge has slope \(4/7\) and its
leading bracket must equal the scalar commutator.

**Corollary PolynomialWeylNonsquareMinusFifteenThroughFortySeven.**
If \(-15\) is not a square in a characteristic-zero field \(k\), a
polynomial Weyl pair with nonzero scalar commutator generates when
one member has Bernstein degree at most \(47\), or has nonzero
constant leading differential coefficient and order at most \(47\).

**Corollary PolynomialWeylFormallyRealThroughFortySeven.** In
particular the two generation assertions hold through \(47\)
over every formally real field.

## The crossing and its first fractional shift

Apply the
[Bernstein and leading-degree bounds through forty-one](POLYNOMIAL_WEYL_BERNSTEIN_DEGREE_THROUGH_FORTY_ONE.md)
and the crossing proof there to a member of Bernstein degree \(42\).
The quadratic bracket mate puts its top into one monomial after
a linear Weyl change over \(k\), as verified below. A pure monomial,
oriented along the position
axis, would give differential order below \(42\) and order
plus leading degree below \(42\). The known bounds exclude that
case, whether its coefficient is constant or nonconstant.
In the mixed case orient the smaller exponent in the derivative.
The individual leading-degree bound now saturates, so the
differential leading term is the top monomial. After
fixed-member differential mate descent the orders do not
divide each other: if the smaller reduced partner's order
divided the fixed member's, the same leading-degree argument
would give at least \(84\) for order plus leading degree.
Writing the resulting orders as \(ga,gb\) and the common
leading-factor degree as \(d\), the infinity calculation
gives \(42=a(g+d)\) with \(g+d\ge16\) and \(a\ge2\).
Thus \(a=2,\ g+d=21\), and the source shear gives a
subrectangular pair of Bernstein degrees \(42,21b\)
with \(b\ge3\) odd.

Its crossing is obtained by the arbitrary-pair argument
using Guccione–Guccione–Valqui, Lemma 7.1, in the
degree-forty-one proof. The numerical tests of their
Proposition 7.3 apply to it. Sums at most twenty have
already been excluded for ratio \(2:b\) there. At sum
twenty-one the only row is
\[
          C_0=(6,15),\quad C_1=(1,0),\quad
          (f_1,f_2)=(4,10),\quad(\rho,\sigma)=(3,-1).
\]
This is the same finite corner enumeration as in that
proof, with the normalized sum set to twenty-one.
In particular the possible smaller coordinates before
the direction test are \(u=3,6,7,9\), the integers
\(3\le u<21/2\) with \(\gcd(u,21)>1\).
The alternative \((7,14)\) with \(C_1=(3,2)\)
would have \((f_1,f_2)=(2,4)\) and gives
\(C_2=(4-1/3,4)\), excluded by Proposition 7.3(6).
The corner row has \(C_0=(g,d)\) since its sum
saturates the support rectangle. Thus \(g=6,\ d=15\).
The common factor is \(R=XL(z)\) with \(z=XY^3\),
\(L(0)\ne0\), \(\deg L=5\).

The homogeneous bracket mate \(F=XYG(z)\) of weight two
satisfies \(\{R,F\}=cR\) with \(c\ne0\). Equivalently
\[
                   3zG'+\left(1-2z\frac{L'}L\right)G=c.
                                                               \tag{1}
\]
In particular \(G(0)=c\ne0\), and cancellation at the
roots of \(L\) gives at least one root of \(G\). Comparing
the highest degrees gives \(\deg G=3\). Hence \(L\)
has at most three distinct roots. Its maximal root
multiplicity \(\mu\) belongs to \(\{2,3,4,5\}\).
The same multiplicities occur in \(L(\eta^3)\) at
nonzero roots, where \(\eta=X^{1/3}Y\).

We use [Guccione–Guccione–Valqui II, Proposition 4.3](https://arxiv.org/pdf/1205.6827v1)
as in the
[order-thirty-two successor proof](MONIC_POLYNOMIAL_ORDER_THIRTY_TWO_EXCLUSION.md).
The scalar commutator is normalized to one, both faces
are nonmonomial, their positive weights have nondividing
ratio \(2:b\), their leading bracket vanishes and their
highest \(Y\)-endpoints have negative diagonal difference.
The proposition shifts by a root of maximal multiplicity
and supplies a successor with larger slope \(s\) in
weights \((1,-s)\). Its normalized highest endpoint is
\[
                 C=(1+\mu/3,\mu),\qquad l=3,\qquad
                 \beta=l(\mu-C_x)=2\mu-3.
\]
The operators are in \(W(l)=k[X^{1/l},X^{-1/l}][Y]\)
with derivation \(d/dX\). At another commuting edge the
same hypotheses of Proposition 4.3 hold by its output,
including the negative diagonal difference of each
corner used below. Its further root shifts retain
the ratio \(2:b\).

## Testing the scalar successors

Write \(N=b+2\), an odd integer at least five.
At a scalar successor of \(C\), positivity and preservation
of the weight ratio give
\[
       N(C_x-\mu s)=1-s,\qquad
       3s=\frac{N(3+\mu)-3}{N\mu-1}.
\]
The gcd cancelled from the displayed fraction divides
\(\beta\). Its reduced denominator \(D\) is coprime
to \(\mu\). It divides every deficit from the maximal
\(Y\)-degree \(2\mu\) in \(P\)'s face, and likewise
from \(b\mu\) in \(Q\)'s face. If \(D>2\mu\), the
first face is a monomial, contrary to the successor
being an edge. If \(\beta\le\mu\), then \(D>N-1\),
so neither \(2\mu\) nor \(b\mu\) is divisible by \(D\).
The two face polynomials in \(YX^s\) then both have
zero constant term, and their bracket cannot be scalar.

For \(\mu=2,3\) we have \(\beta\le\mu\).
For \(\mu=4\), \(\beta=5\). If the cancelled gcd
is one, \(D=4N-1>8\). Otherwise it is five and
\(N\equiv4\bmod5\). The first admissible odd \(N\)
is \(9\), and the next is \(19\), for which
\(D\ge15>8\). So a scalar successor at \(\mu=4\)
requires \(N=9,\ b=7,\ D=7\), and the slope
formula gives \(s=4/7\).
For \(\mu=5\), \(\beta=7\). Cancellation of seven
requires \(N\equiv3\bmod7\); the first admissible
odd \(N\) is \(17\), giving \(D\ge12>10\).
Without cancellation \(D=5N-1>10\).
Thus the only possible scalar successor is the
one asserted.

## Testing the commuting successors

At a commuting successor, a single-root common face
allows a root shift preserving \(C,l\), and its
slope has \(ls\in\mathbb Z\) by adjacency of the
coefficients. With \(l=3\), \(1/3<s<C_x/\mu\),
such a slope is possible only for \(\mu=2\), at
\(s=2/3\), and the next successor must have larger
slope. When the common face has multiple distinct
roots, the positive weight equation in the
[successor proof](MONIC_POLYNOMIAL_ORDER_THIRTY_TWO_EXCLUSION.md#the-coefficient-equation-and-the-lattice-table)
gives
\[
  W'-\frac{M}{\mu}\frac{f'}{f}W=c,\qquad
  M=\deg W\ge2,\qquad
  ls=l-\frac{\beta M}{\mu(M-1)},\qquad c\ne0.               \tag{2}
\]
Write \(f=\eta^j J(\eta^D)\) with \(J(0)\ne0\),
where \(D\mid\mu-j\) is the reduced denominator of \(ls\).
In addition to \(M\le\beta\mu+1\), we have
\(M\equiv1\bmod D\). Indeed the projection of \(W\)
to exponents congruent to one modulo \(D\) also solves
the same equation with right side \(c\). This follows
by averaging \(\zeta^{-1}W(\zeta\eta)\) over \(\zeta^D=1\).
The equation forces every polynomial solution to have
degree \(M\), since it vanishes at the multiple distinct
roots of \(f\) and its leading term at infinity
must cancel. Thus that projection has degree \(M\).

With \(l=3,\ \beta=2\mu-3,\ 1/3<s<C_x/\mu\), the
same finite rational test as in the
[order-thirty-two lattice table](MONIC_POLYNOMIAL_ORDER_THIRTY_TWO_EXCLUSION.md#the-coefficient-equation-and-the-lattice-table)
gives these rows before the projection test:

| \(\mu\) | \(M\) | \(s\) | \(D\) | \(j\) |
|---:|---:|---|---:|---|
| 2 | 2 | \(2/3\) | 1 | \(0,1\) |
| 3 | 3 | \(1/2\) | 2 | \(1\) |
| 3 | 4 | \(5/9\) | 3 | \(0\) |
| 4 | 4 | \(4/9\) | 3 | \(1\) |
| 4 | 6 | \(1/2\) | 2 | \(0,2\) |
| 4 | 16 | \(5/9\) | 3 | \(1\) |
| 5 | 5 | \(5/12\) | 4 | \(1\) |
| 5 | 8 | \(7/15\) | 5 | \(0\) |
| 5 | 15 | \(1/2\) | 2 | \(1,3\) |

All rows with \(M=\mu\) force a simple root, incompatible
with multiplier one. The rows \((\mu,M)=(4,6),(5,8)\)
violate the new congruence on \(M\).
For \((3,4)\) the polynomial is \(\eta^3-u_0\) and a
simple-root shift at \(s=5/9\) gives the new endpoint
\(C'=(8/9,1)\) on denominator nine.
The successor has positive weights and slope in
\((5/9,8/9)\). A scalar leading bracket is impossible:
its denominator formula with \(l=9,\mu=1,\beta=1\)
gives \(N-1>b\), so neither face can have a constant
term in its polynomial variable. A commuting edge of
degree one can only be a single-root edge, forcing
\(9s\) integral, preserving both \(C'\) and the lattice.
Thus there are at most two such steps, at slopes
\(2/3\) and \(7/9\), after which the required successor
does not exist. To justify these successive applications
of Proposition 4.3, the endpoints \(C'\) have negative
diagonal difference \(8/9-1\), and the proposition
preserves the remaining positive-weight, edge and
nondividing-ratio hypotheses.

For \((\mu,M)=(4,16)\) the polynomial is
\(\eta(\eta^3-u_0)\) and (2) would make \(1/f^4\)
a rational derivative up to a constant.
At zero the residue of
\(\eta^{-4}(\eta^3-u_0)^{-4}\) is nonzero.
For \((5,15)\) the same reasoning requires zero
residues for \(1/f^3\). If \(j=3\) the polynomial
is \(\eta^3(\eta^2-u_0)\); at zero the coefficient
of \(\eta^8\) in \((\eta^2-u_0)^{-3}\) is nonzero.
If \(j=1\) the polynomial is \(\eta J(\eta^2)\).
Vanishing of the residue at zero forces \(J'(0)=0\),
giving \(f=\eta(\eta^4-u_0)\) up to scale.
At a nonzero root \(\alpha\) the residue is proportional
to \(3f''(\alpha)^2-f'(\alpha)f'''(\alpha)=960\alpha^6\),
which is nonzero. These exclusions exhaust the table.
In the only possible single-root case \(\mu=2\), the
shift at \(2/3\) keeps the corner and denominator and
leaves no row of strictly larger slope.

## The two roots and real coefficients

The required scalar successor has \(\mu=4\). Thus the
first \(L\) has multiplicities four and one, and \(G\)
has the form
\[
                 L=L_0(z-u)^4(z-v),\qquad
                 G=\gamma(z-u)(z-v)(z-w).
\]
The constant and simple-zero comparisons of (1) at
\(z=0,u,v\), after rescaling \(u=1\), give
\[
                   5-5v-5w+4vw=0,\qquad
                   v^2-v-vw+2w=0 .
\]
Equivalently, since \(v\ne1\), the first gives
\(w=5(v-1)/(4v-5)\), and the second gives
\(w=-v(v-1)/(2-v)\). The excluded denominator values
\(4v=5\) and \(v=2\) do not satisfy the original
equations. Equating the expressions yields
\(2v^2-5v+5=0\). Replacing the normalized \(v\)
by \(v/u\) gives the claimed equation.
The slope and quotient test gives \(g=6,d=15,b=7\),
so the derived Bernstein degrees are \(42,147\).

The quadratic homogeneous bracket mate at Bernstein weight
can be chosen over \(k\) by Guccione–Guccione–Valqui,
Theorem 4.1(4). A repeated linear factor would give nilpotent
Hamiltonian action on fixed homogeneous degree. Otherwise its
action on degree-one forms has eigenvalues \(\lambda,-\lambda\)
over an algebraic closure, and on homogeneous degree \(n\)
the eigenvalues are \((i-j)\lambda\), \(i+j=n\). The leading
symbol has eigenvalue \(1\) or \(-1\), so \(\lambda\in k^\times\).
The two independent eigenlines are therefore defined over
\(k\), and can be rescaled to canonical linear coordinates
over \(k\). The axis orientation, mate subtraction and infinity
shear likewise take place over \(k\).

The crossing's common factor is defined over \(k\). Normalize
the coefficients of the leading monomial in the edge symbols
to one. The ratio of the normalized mate symbol to the
\((b-1)/2\)-th power of the normalized first symbol is the
normalized \(R\), by the common-power comparison over the
algebraic closure. Hence the original \(L\) is defined over
\(k\), even though the subsequent shifts may use roots over the
algebraic closure. Its roots of multiplicities four and one
are separately defined over \(k\): the monic gcd of \(L,L'\)
is \((z-u)^3\), and \(v\) then follows from the coefficient
of degree four in \(L\). Thus \(v/u\in k\) and
\[
                         (4v/u-5)^2=-15.
\]

## The nonsquare field bound through forty-seven

Assume \(-15\) is not a square in \(k\).
The residual theorem excludes Bernstein degree \(42\).
Degree \(43\) is excluded by Bernstein mate
descent, because a nondividing pair with smaller degree
\(43\) cannot have a degree gcd at least sixteen.
For constant-leading orders \(42,43\), repeat the
last step of the degree-forty-one proof. A smaller
constant-leading order after a binomial change or a
dividing mate is at most \(21\). At intermediate
weights a nonconstant leading coefficient would have
order plus coefficient degree at most the starting
order, contrary to the Bernstein bound just proved and
the infinity-shear leading-degree inequality.
The weight-one branch would have Bernstein degree
equal to the differential order after depression.
This proves both bounds through \(43\).

For degree \(44\) and then degree \(46\), use the same
axis-orientation, saturation and crossing argument as
in the degree-forty-one proof. The quadratic mate
diagonalizes over \(k\) by the preceding eigenvalue
argument, so the inductive bounds apply in that actual
presentation before any fractional shift. The individual
leading-degree bound at the induction step follows
from the preceding Bernstein bound by infinity shear
and the dividing-partner descent. Thus if the first
member has Bernstein degree \(n=44\) or \(46\), the
subrectangular pair has order ratio \(2:b\) and normalized
vertex sum \(n/2\), since \(n<3\cdot16\).
The crossing exclusions through twenty still apply.
At sum twenty-one the preceding exclusion over \(k\)
uses only the crossing's symbols and \(L\in k[z]\),
so it also applies when the outer Bernstein vertex
has larger sum.

At sum twenty-two the integer conditions of
Guccione–Guccione–Valqui, Proposition 7.3, give
\[
       C_0=(8,14),\quad C_1=(1,0),\quad
       (f_1,f_2)=(4,7),\quad(\rho,\sigma)=(2,-1).
\]
For this additional sum the smaller coordinate \(u\)
must be among \(4,6,8,10\), the integers \(3\le u<11\)
with \(\gcd(u,22)>1\). Substituting the possible
integer multiples of the resulting \(C_0\), and solving
the equal-weight equation with \(0\le s_0<r_0<u\),
gives the displayed row. The factor is \(XL(XY^2)\)
with \(\deg L=7\). Its homogeneous bracket mate
gives the equation
\[
                 2zG'+(1-zL'/L)G=c\ne0 .
\]
It has \(\deg G=3\). Thus \(L\) has at most three
distinct roots, and its maximal root multiplicity
\(\mu\) lies between three and seven.

The first shift of Proposition 4.3 has lattice
denominator \(l=2\) and endpoint \(C=(1+\mu/2,\mu)\),
with \(\beta=\mu-2\). As in the previous cases,
\(\beta\le\mu\) excludes a scalar successor.
No single-root slope is possible, since it would
require integral \(2s\) strictly between \(1\) and
\(1+2/\mu<2\). The multiple-root lattice test, including
\(M\equiv1\bmod D\), leaves the following rows:

| \(\mu\) | \(M\) | \(D\) | \(j\) |
|---:|---:|---:|---:|
| \(3\le\mu\le7\) | \(\mu\) | \(\mu-1\) | \(1\) |
| 5 | 10 | 3 | 2 |
| 6 | 9 | 4 | 2 |
| 7 | 21 | 4 | 3 |

The \(M=\mu\) rows have all roots simple, impossible
with multiplier one in (2). The next row gives
\(f=\eta^2(\eta^3-u_0)\) with multiplier two, excluded
by the
[degree-thirty-two successor calculation](MONIC_POLYNOMIAL_ORDER_THIRTY_TWO_EXCLUSION.md#removing-the-remaining-successors).
The last row gives \(\eta^3(\eta^4-u_0)\) with multiplier
three. The residue of \(1/f^3\) at zero is nonzero
because it uses the nonzero coefficient of \(\eta^8\)
in \((\eta^4-u_0)^{-3}\).

The remaining row has \(\mu=6,\ M=9,\ D=4,\ j=2\),
so \(f=\eta^2(\eta^4-u_0)\) and the slope is \(5/8\).
The maximal root is zero with multiplicity two.
Proposition 4.3 uses the identity shift, so the
operators stay in \(W(2)\) even though the proposition
can be applied in a larger lattice. Its successor
endpoint is \(C'=(3/2,2)\) with \(\beta=1\),
and the slope lies strictly between \(5/8\) and \(3/4\).
It has no integral \(ls\). A scalar leading bracket
is impossible since \(\beta\le\mu\).
The integer bound \(2\le M\le3\) gives only the slopes
\(1/2,5/8\) in (2), neither strictly larger than \(5/8\).
The applications of Proposition 4.3 here use the
positive nondividing ratio and scalar commutator,
nonmonomial faces and zero leading bracket, and the
negative diagonal differences of \(C,C'\).
This excludes sum twenty-two. Sum twenty-three is prime
and cannot satisfy the corner gcd condition.

Thus Bernstein degrees \(44\) and \(46\) are excluded
in succession. Their next odd degrees \(45\) and \(47\)
have no proper divisor at least sixteen, so are excluded
by Bernstein mate descent at each step. For the
corresponding constant-leading orders use the
same weighted three-way reduction as for \(42,43\).
The binomial and intermediate changes have differential
orders at most \(23\). A nonconstant leading coefficient
in the intermediate change contradicts the updated
individual leading-degree bound from the Bernstein step.
This proves both bounds through \(47\).

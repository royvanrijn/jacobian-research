# The next cusp cancellation and graded fields in degree forty-two

**Theorem PolynomialWeylFortyTwoSecondCuspCancellation.**
Use \(P,Q^\ast,S\) in the actual degree-forty-two
[normalized cusp relation](POLYNOMIAL_WEYL_FORTY_TWO_NORMALIZED_CUSP_RELATION.md).
Keep its variables \(t,E,C\) and leading coefficients \(\alpha,\beta\).
Let \(\eta\ne0\) satisfy \(S_0^3=\eta Q_0P_0^5\) for the
lowest Euler components. Then there is a scalar \(a_1\) such that
\[
               T=S^3-\eta Q^\ast P^5-a_1 P^2 Q^\ast S
\]
has \(t\)-degree \(51\), lowest Euler grade \(-147\), and
subrectangular support with corner \((99,246)\). In the
rational Euler algebra, put \(R_0=P_0^3 S_0^{-1}\), of
grade \(-3\). The lowest Euler component of \(T\) is
\[
                       T_0=R_0^{49}(aC+b),\qquad a\ne0.       \tag{1}
\]
Thus the coefficient of \(R_0^{49}\) in its rational
normal form has total pole order at most one.

In the six-root packet description of the finite indicial
polynomial, that pole can only come from an extremal packet.
In step-three coordinates, up to translation, its sorted
indices must be \(\{0,17,21,34,38,55\}\), and its excess
divisor contributes a simple deficit at index \(49\).
Deficits count with multiplicity and can be canceled by
the integer roots \(0,\ldots,146\) of the falling factorial
\((C)_{147}\). Across both packets, at most one deficit
remains after that cancellation.

## Proof of the second cancellation

Write \(P=t^6\alpha H^2+\cdots\),
\(Q^\ast=t^{21}\beta H^7+\cdots\) and
\(S=t^{17}H^5 B+\cdots\) as in the cusp relation.
Let \(b_\infty\ne0\) be the leading coefficient of \(B\).
Choose the formal branch \(H^{2/3}=E^{10}(1+O(E^{-3}))\)
at \(E=\infty\). Its quotient equation gives
\[
 \left(\frac{B}{H^{2/3}}\right)'=\frac{\beta}{\alpha}H^{-2/3},
 \qquad
 \frac{B}{H^{2/3}}=b_\infty-\frac{\beta}{9\alpha}E^{-9}
                          +O(E^{-12}).                       \tag{2}
\]
This follows by expanding \(H^{-2/3}=E^{-10}(1+O(E^{-3}))\)
and integrating term by term.

Begin with \(T_{\mathrm{raw}}=S^3-\eta Q^\ast P^5\).
Its finite grade \(-153\) cancels by the definition of
\(\eta\). The \(t^{51}\)-coefficient is
\(H^{15}(B^3-\eta\beta\alpha^5 H^2)\).
Its coefficient of \(E^{255}\) is zero because that would
be the \((t,E)\)-symbol of the monomial at finite grade
\(-153\) and \(t\)-degree \(51\). Therefore
\(\eta\beta\alpha^5=b_\infty^3\). Equation (2) gives
\[
        \deg_E(B^3-\eta\beta\alpha^5 H^2)=21,\qquad
        \deg_E (T_{\mathrm{raw}})_{51}=246.
\]
In particular the \(t\)-degree is exactly \(51\), and the
corner term produces a nonzero monomial \(X^{99}Y^{246}\)
of Euler grade \(-147\). Normal-order corrections lower
\(t\)-degree, so cannot supply that \(t^{51}E^{246}\) term
from another monomial of lower or equal \(t\)-degree.

Since \([P,S]=-2Q^\ast\), the lowest Euler component of
\([P,T_{\mathrm{raw}}]\) is \(-6Q_0S_0^2\) at grade
\(-165\). In fact the three terms in \([P,S^3]\) contribute
the same commuting product there, and \([P,Q^\ast P^5]\)
has higher lowest grade. If \(T_{\mathrm{raw}}\) has
lowest finite grade \(a<-147\), its leading homogeneous
component \(X^a A(C)\) must commute with \(P_0\).
Comparing one coefficient below the leading degree in
\(f(C+a)A(C)=f(C)A(C-18)\) gives
\(\deg A=-5a/3\), so \(3\mid a\). Its \(t\)-degree
as a polynomial is \(3a+2\deg A=-a/3\le51\).
Since grade \(-153\) has already canceled, the only
possibility below \(-147\) is \(-150\).

The actual polynomial \(P^2Q^\ast S\) has \(t\)-degree
\(50\), and its grade \(-150\) component commutes with
\(P_0\). It is proportional to any other such homogeneous
commuting component by the rational-period equation
for the quotient of coefficient functions. Thus one
subtraction cancels grade \(-150\) if it occurs.
The commutator \([P,P^2Q^\ast S]\) has lowest grade
at least \(-162\), so this correction does not change
the grade \(-165\) component of \([P,T]\).
It leaves the \(t^{51}E^{246}\) corner in place.
The same leading-degree bound rules out any grade
below \(-147\) after this subtraction, and that corner
gives equality. The inequalities \(i-j\ge-147\) and
\(3i-j\le51\) give \(i\le99,\ j\le246\), proving
subrectangularity.

The lowest Euler grades of \(P_0,Q_0,S_0,R_0\) are
\(-18,-63,-51,-3\), respectively, and they commute.
Each of the first three is proportional to the
power of \(R_0\) of the same grade, by uniqueness
up to a scalar of a homogeneous commuting component
in the rational Euler algebra. Write
\(T_0=R_0^{49} A(C)\) with \(A\) rational.
Since
\([P_0,T_0]=-6Q_0 S_0^2\), and \(Q_0 S_0^2\) is
a scalar multiple of \(P_0 R_0^{49}\), one obtains
\[
                       A(C)-A(C-18)=\text{nonzero constant}.
\]
Subtracting a linear particular solution leaves a
rational function of period eighteen, hence a
constant. This proves (1) and its pole bound
because \(T_0\) has a polynomial coefficient.

## The shape of a possible packet deficit

Write \(R_0^{49}=X^{-147}U(C)\). The commutation
equation with \(P_0\), after factoring out
\((C)_{147}\) from \(U(C)\) as a rational function,
gives the excess divisor on each \(3\mathbb Z\)-orbit
as
\[
                \frac{1-z^{49}}{1-z^6}\,\mathcal F(z),        \tag{3}
\]
where \(\mathcal F\) records the twelve excess roots
of \(f\) on that orbit. This follows by projecting the
root-and-pole equation to step-three indices.
Factoring the falling factorial cancels the integer
base part since
\((C-147)_{18}(C)_{147}=(C)_{18}(C-18)_{147}\).
There are no additional divisor contributions on
orbits absent from the right side: their records
would be invariant under a shift by six with finite
support.

For residue \(i\bmod6\) in step-three coordinates, let
\(r_{i,j}\) be the \(j\)-th sorted index of its roots,
as in the cusp relation. The coefficient of (3)
at index \(n\equiv i\bmod6\) equals the number of
\(r_{i,j}\le n\) minus the number of
\(r_{i-1,j}+49\le n\).
All matched differences
\(r_{i,j}-(r_{i-1,j}+49)\) are at most six by the
previous packet bound. They are divisible by six.
Thus a strict reversal is precisely
\(r_{i,j}=r_{i-1,j}+55\), contributing one
deficit at \(n=r_{i-1,j}+49\).
No positive contribution from another rank can
cancel such a deficit in this residue: earlier
ranks have second index at most \(n\), and later
ranks have first index greater than \(n\).

The cusp inequalities give a bound \(55\) by
each of the three length-three paths of offsets
\((3,1,1),(1,3,1),(1,1,3)\) from residue \(i\)
to \(i-1\), with step costs \(21,17\).
Equality \(r_{i,j}-r_{i-1,j}=55\) forces equality
in all three paths. With \(r_{i-1,j}=0\), the six
values in the order \(i,i+1,\ldots,i+5\) are
\(55,38,21,34,17,0\), giving the claimed set.
There cannot be two such rises in one packet,
because each of the other cyclic rises is at
least \(-17\) and the six rises sum to zero.
The only divisor contribution omitted in (3)
is the base \((C)_{147}\), which cancels at most
one deficit at each integer \(0,\ldots,146\).
The total remaining pole order is at most one
by (1).

## The two exact associated graded fields

**Theorem PolynomialWeylFortyTwoGradedFields.**
Let \(\mathcal A\) be the image algebra generated by
the actual pair \(P,Q^\ast\).
Its associated graded subalgebra for \(t\)-degree
(the \((3,-1)\)-filtration) in \(\mathbb C[X,Y]\)
is
\[
                    \mathbb C[R^2,R^7,R^5 J],
         \qquad R=XH(z),\quad J=XY\,D(z),\quad z=XY^3,       \tag{4}
\]
where \(H(z)=(z-u)^4(z-v)\) and \(B(E)=E D(E^3)\).
Its fraction field has degree five in the ambient
symbol field \(\mathbb C(X,Y)\). The degree-five
extension is the base change by \(R\) of a rational
cover with monodromy \(A_5\) and branch cycles of
types \((5),(3,1,1),(3,1,1)\).

For the Bernstein filtration the associated graded
subalgebra instead is
\[
       \mathbb C[M^6,M^{21},M^{17},N],\qquad
       M=X^2Y^5,\quad N=X^{99}Y^{246}.                      \tag{5}
\]
Its fraction field has degree three in the ambient
symbol field. These two degrees concern commutative symbol
fields. For lowest finite Euler grade the associated graded
algebra is generated by \(R_0^6,R_0^{21},R_0^{17},R_0^{49}(aC+b)\).
Its division closure in the rational Euler algebra is generated
by \(C,X^3\) and has index three there. The graded calculations
leave the degree of the actual embedded division algebra open.

For (4), the elements \(P^i(Q^\ast)^j S^k\) with
\(i,k\ge0\), \(j=0,1\), form a basis of
\(\mathcal A\) by triangularity in the ordinary
ordered powers of \(Q^\ast\).
Their top ordered terms are \(P^i(Q^\ast)^{j+2k}\)
with coefficient one, and every smaller mate
power is covered inductively.
The corresponding \(t\)-symbols, up to constants,
are \(R^{2i+7j+5k}J^k\).
They are independent because \(R,J\) are algebraically
independent, and the exponents at fixed \(k\)
are distinct by parity of \(j\).
Thus all of the associated graded algebra consists
of their span, proving (4).
Here \(R=t^3H(E)\) and \(J=t^2 B(E)\), yielding
the displayed expressions in \(X,Y\).
The \(t^2\) coefficient congruence forces
\(B(E)=E D(E^3)\), and its simple zeros at the
roots of \(H\), together with \(\deg B=10\), give
\[
                      D(z)=(z-u)(z-v)L(z),
                      \qquad \deg L=1.
\]
The linear factor does not vanish at \(u\) or \(v\),
since those would be higher zeros of \(B\).
It does not vanish at zero because (1) of the
normalized cusp relation gives \(B'(0)=\beta/\alpha\).

Put
\[
             w=\frac{J^3}{R^2}
               =\frac{z(z-v)L(z)^3}{(z-u)^5}.
\]
This rational function of \(z\) has degree five.
The pair \(R,z\) is algebraically independent.
To describe the degree-five field over
\(\mathbb C(R,J)\), first take the extension
\(\mathbb C(z)/\mathbb C(w)\), base change by the
independent \(R\), and then adjoin the cube root
\(J\) of \(R^2 w\).
The cubic is irreducible over either rational
field with \(R\) adjoined, as the \(R\)-valuation
of \(R^2 w\) is two. Thus base change by that
cubic does not change the degree five.
Moreover \(\mathbb C(X,Y)=\mathbb C(R,J,z)\),
since \(X=R/H(z)\) and \(Y=J/(X D(z))\).
This also proves independence of \(R,J\) used above.

The only pole of \(w\) has order five at \(u\).
Its zeros are simple at \(0,v\) and triple at the
root of \(L\).
At \(z=\infty\), equation (2) gives a nonzero
finite value with a first nonconstant term in
\(z^{-3}\), since \(w=B(E)^3/H(E)^2\) and
\(z=E^3\). Thus the three ramification indices
five, three, three exhaust the Riemann--Hurwitz
sum for a degree-five rational map.
Its transitive monodromy group is in \(A_5\)
and contains both a five-cycle and a three-cycle,
hence is \(A_5\). Indeed a proper subgroup would
have order fifteen or thirty by divisibility;
a subgroup of order fifteen would be cyclic by
Sylow's theorem, and order thirty contradicts
simplicity of \(A_5\).
The cubic base change does not change this group:
the \(A_5\) Galois closure, with the independent
\(R\) adjoined, has no intermediate extension of
degree three, and the cubic base change is Galois
because the ground field contains the cube roots of
unity. Thus their intersection
over \(\mathbb C(w,R)\) is the ground field.

For (5), use the triangular basis
\(P^i(Q^\ast)^j S^k T^\ell\) with
\(i,\ell\ge0,\ j\in\{0,1\},\ k\in\{0,1,2\}\).
Its top ordered powers of \(Q^\ast\) are
\(j+2k+6\ell\), with coefficient one.
Its Bernstein symbols are nonzero scalar multiples
of \(M^{6i+21j+17k}N^\ell\) by the four corner terms.
Their exponent vectors are distinct. The determinant
of \((2,5),(99,246)\) is \(-3\), which fixes
\(\ell\) in an equality of exponent vectors.
Reduction modulo three fixes \(k\), and reduction
modulo six then fixes \(j\).
Thus no terms of a fixed Bernstein degree cancel
as polynomials, giving (5).
Its fraction field is \(\mathbb C(M,N)\) because
\(\gcd(6,21,17)=1\). The determinant of the
exponent lattice gives the degree three.

The lattice index three in (5) agrees with a refinement
of the \(t\)-edge field calculation. Refine \(t\)-degree
by taking the largest \(Y\)-exponent in its leading
symbol. The four elements \(P,Q^\ast,S,T\) have their
Bernstein corner on the \(t\)-top, and those corners
multiply. The same triangular basis used for (5)
has distinct refined leading monomials, so the
refined value group of \(\mathcal A\) has index three
in the ambient exponent lattice.
For the first \(t\)-degree alone the ambient symbol
field and the image symbol field have the same value
group, since the degrees of \(R,J\) are \(3,2\).
Their degree-zero fields are respectively
\(\mathbb C(z)\) and \(\mathbb C(w)\).
The refinement chooses \(z=\infty\), because
\(z\) has first degree zero and positive \(Y\)-exponent.
Its restriction to the image chooses the finite
value \(w(\infty)\), and the order three of
\(w-w(\infty)\) at that place is precisely the
secondary ramification index. Thus the degree-five
field calculation includes the index-three place
that accounts for the refined monomial calculation.

## The finite Euler graded algebra

For lowest Euler grade, the same triangular basis
\(P^i(Q^\ast)^j S^k T^\ell\) has leading term a nonzero
scalar times \(R_0^m F_\ell(C)\), where
\(m=6i+21j+17k+49\ell\) and \(F_\ell\) is a polynomial
of degree \(\ell\) with nonzero leading coefficient.
Indeed each \(T\) supplies one linear factor by (1), and
moving rational Euler monomials past them only shifts \(C\).
At fixed \(m\), different \(\ell\) give different polynomial
degrees. If \(\ell\) also agrees, reduction modulo three
and then modulo six fixes \(k,j\), and \(i\) follows.
Thus no lowest-grade terms cancel on the basis, giving
the claimed generators of the finite graded algebra.

Because \(\gcd(6,21,17)=1\), the graded division closure
contains \(R_0\). Multiplying \(T_0\) by \(R_0^{-49}\)
then gives \(aC+b\), so the closure is generated by \(C,R_0\).
Since \(R_0=X^{-3}\) times a nonzero rational function of
\(C\), it is also the division algebra generated by \(C,X^3\).
In the ambient rational Euler algebra, \(X\) normalizes
this division algebra by the shift on \(C\), so the
left span of \(1,X,X^2\) over it is closed under products.
The three elements are linearly independent by their distinct
exponent residues in the formal skew Laurent series in \(X\)
over \(\mathbb C(C)\). This finite-dimensional span is
a division algebra: right multiplication by a nonzero
element is injective and therefore surjective on the span.
It contains the ambient generators, giving index three.

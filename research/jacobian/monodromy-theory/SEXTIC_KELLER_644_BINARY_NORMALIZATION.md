# Normalization of the reduced binary (6,4,4) profile

**Theorem KellerProjectivelyBirationalBinaryLeadingExclusion.**
Let a characteristic-zero three-variable Keller map have ordinary
component degrees \(d_1,d_2,d_3\) with \(\min_{i\ne j}(d_i+d_j)>2\).
Suppose its highest forms \(A,B,C\) are binary in the same linear
variables \(x,y\), and \(k[x,y]\) is finite over \(k[A,B,C]\).
If the degree-zero subfield of \(k(A,B,C)\) is all of \(k(y/x)\),
there is no such map.

**Theorem SexticKeller644ReducedBinaryNormalizationReduction.**
Suppose a characteristic-zero three-variable Keller map of ordinary
component degrees \((6,4,4)\) has binary highest forms \(A,B,C\),
with the two quartics coprime and squarefree over the algebraic closure.
Then over the algebraic closure the quartic pencil has a member
\(K^2\) with \(K\) quadratic, and after an invertible linear change
on the lower two components, the highest tuple has form
\[
               (aK^3+bKC,\ K^2,\ C).
\]
Here \(K,C\) are coprime, and \(k(K,C)\) is the full fraction field
of the highest tuple. In particular, the projective degree of the
field of the leading tuple is four.

Work over the algebraic closure \(k\). We use the first transverse
weight and finite normalization lift from the
[binary reduced lower pencil proof](SEXTIC_KELLER_BINARY_REDUCED_LOWER_PENCILS.md#lifting-the-finite-normalization-at-the-first-new-weight).

## Projectively birational homogeneous tuples

Put \(E=k(A,B,C)\) and let \(g=\gcd(d_1,d_2,d_3)\).
A product of integral powers, possibly negative, of \(A,B,C\)
has degree \(g\). This product and \(t=y/x\) generate \(k(t,x^g)\),
and all three highest forms lie in this field. Hence \(E=k(t,x^g)\).
For \(g=1\), apply the binary finite birational exclusion.
For \(g>1\), lift \(x^g,y^g,x^{g-1}y\) at the first transverse
tie to \(U,V,M\). They are integral over the image ring by
finiteness. The lifted relation is \(M^g=U^{g-1}V\).
Any positive-weight factor common to \(U,V\) would restrict
to a nonzero positive-degree factor common to \(x^g,y^g\):
it cannot restrict to zero because the products restrict
nontrivially. Thus \(U,V\) are coprime. Unique factorization
forces each of them to be a \(g\)-th power up to scalar,
by \(\gcd(g,g-1)=1\). The roots have weight one and cannot
involve the third variable. Restriction then gives \(U=x^g\),
\(V=y^g\), and \(M=x^{g-1}y\) by the relation and restriction.
The image field lifts without change since these generate \(E\).
This contradicts the defining tie.

## Index two in the projective field

For the \((6,4,4)\) tuple put \(t=y/x\), \(R=B/C\),
\(u=A/C\), and \(\eta=A^2/C^3\). Then \(E=k(u,R,\eta)\)
with degree-zero field \(E_0=k(R,\eta)\), since \(C=u^2/\eta\).
Let \(m=[k(t):E_0]\). Coprimality of \(B,C\) makes
\([k(t):k(R)]=4\), so \(m\in\{1,2,4\}\).
The previous theorem excludes \(m=1\), since coprime binary
forms give finiteness.

If \(m=2\), the involution of the binary projective line can
be diagonalized so that \(E_0=k(t^2)\).
Represent the degree-two rational function \(R\) of \(t^2\)
by binary quadratics in \(x^2,y^2\). They are coprime after
substitution, so \(B,C\in k[x^2,y^2]\).
As \(A^2/C^3\) is invariant under \(y\mapsto-y\), \(A\) is
either even or odd for this involution. In either case
\([k(x,y):E]=4\). Indeed \(u\) is an independent radial
variable over \(k(t)\), and \(u=x^2 f(t)\) for a nonzero
rational function \(f\). Adjoining \(t\) has degree two
and then adjoining \(x\) has degree two, since \(u/f(t)\)
is not a square in \(k(t,u)\).

If \(A\) is odd, \(E\) is contained in the field fixed by
\((x,y)\mapsto(i x,-i y)\), where \(i^2=-1\).
This fixed field is \(k(x^4,y^4,xy)\): the original
field has degree four over it, as seen by adjoining \(x\),
and the four indicated automorphisms preserve it.
The degree count gives equality with \(E\). Lift the
three integral generators at the first tie to \(U,V,W\)
with \(UV=W^4\). The restrictions of \(U,V\) are
\(x^4,y^4\), so they are coprime by the same argument
as above. Both are fourth powers of weight-one
polynomials up to scalars. Again the lifts necessarily
equal their restrictions, a contradiction.

If \(A\) is even, the degree count gives \(E=k(x^2,y^2)\).
At a first tie, the two weight-two generators can change
only at weight two in the third variable, so their
lifts are \(U=x^2+\alpha z,\ V=y^2+\beta z\), with
at least one coefficient nonzero. Suppose \(\alpha\ne0\)
by exchanging \(x,y\) if needed. If \(\beta\ne0\), take
coordinates \(x,y,U\), and fix the generic values of \(U\)
and \(V-(\beta/\alpha)U=y^2-(\beta/\alpha)x^2\).
The [weighted rational slice](KELLER_LEADING_PAIR_RATIONAL_SLICE.md)
for the actual lower pair forces exactness of the
relative area differential on that generic conic.
Indeed the leading lower forms remain independent
functions of \(U,V\), so their differentials multiply
the relative differential by a nonzero constant
on the generic level. But in linear coordinates
\(s=y-\sqrt{\beta/\alpha}\,x,\ t_1=y+\sqrt{\beta/\alpha}\,x\)
the conic is \(st_1=\text{constant}\), with area
differential proportional to \(ds/s\). Its simple
poles contradict exactness.

If \(\beta=0\), retain coordinates \(U,y,x\) of weights
two, one, one. The leading tuple depends only on \(U,y^2\).
At the next increase of the weight of \(x\), the two
normalization generators may only lift to \(U+\gamma x\)
and \(y^2+\delta x\), at weight two in \(x\).
If \(\delta=0\), replace \(U\) by \(U+\gamma x\) and
increase the weight of \(x\) again. No generator of
weight two can then change, a contradiction.
If \(\delta\ne0\), use both lifts as new coordinates of
weight two, with \(y\) still of weight one. The first
tie in \(y\) can only add linear terms in \(y\) to the
two generators at weight two. These can again be
taken as the new coordinates, after which increasing
the weight of \(y\) beyond two cannot change either
generator. This also contradicts the existence of a
tie. All triangular or affine changes in this argument
have nonzero constant Jacobian and preserve the
current weights. Each tie has positive expected
Jacobian degree since the two fixed weights sum to
at most four and the moving weight is at most six.

## Projective index four

Now \(E_0=k(R)\), and \(\eta\) is a rational function
of \(R\). Its order at \(R=\infty\) is odd because \(C\)
is squarefree. All other places where it has odd
order are square members of the quartic pencil.
Indeed along a factor of a member other than \(C\)
the order of \(\eta\) is even. This is the product of
the local member multiplicity with the order on the
parameter line, so every factor of an odd-order member
has even multiplicity.

There are at most three square members. To see this,
choose two distinct ones \(K_1^2,K_2^2\). They are
coprime. Any other member factors as a product
of two distinct coprime members of the quadratic
pencil spanned by \(K_1,K_2\). If the quartic is a
square, both quadratic members must be squares
of linear forms. A binary quadratic pencil contains
at most two such members, since the discriminant is
a quadratic and a pencil of two independent squares
has nonzero discriminant off those two points.
The unordered pair of quadratic squares determines
at most one other quartic member. Since the sum of
the orders of \(\eta\) is zero, there are one or three
square members of odd order.

In the three-member case, write these members as
\(S_i=K_i^2\), \(i=1,2,3\). At every factor of a
\(K_i\), \(A\) has multiplicity at least that of \(K_i\),
by the odd-order calculation. Coprimality of the
three members and the degree of \(A\) give
\(A\) proportional to \(K_1K_2K_3\).
At the first transverse tie, the lifts \(L_{S_i}\)
are still pairwise coprime by their restrictions,
and their product is a square up to scalar by the
lifted identity for \(A^2\). Each \(L_{S_i}\) is
therefore square, with root \(P_i\) of weight two.
The pencil relation has all coefficients nonzero,
so by rescaling the \(P_i\) it has form
\(P_3^2=P_1^2+P_2^2\).
The two factors \(P_1\pm iP_2\) are coprime and hence
are squares of weight-one polynomials up to scalars.
They cannot involve the third variable, so \(P_1,P_2\)
and the entire lifted lower pencil are binary in the
original variables. The same holds for the lift of
\(A\) by the square identity and restriction.
This contradicts the first tie.

In the one-member case, write that member as \(K^2\)
and use \(R_0=K^2/C\) as parameter. The rational
function \(\eta/R_0\) has even orders everywhere, so
is a square in \(k(R_0)\). Therefore \(u/K\), whose
square is \(\eta/R_0\), belongs to \(k(R_0)\).
This gives \(E=k(K,C)\).
The map \((K,C)\) of the binary plane is onto.
For two nonzero values, choose a projective point of
the specified member of \(K^2/C\) away from its zero
and pole divisors, and rescale to set \(K\) to the
first value. The other value is then forced. For one
zero value, use a point of the corresponding zero
divisor, disjoint from the other in the projective
line, and rescale. The origin gives the last value.
Polynomial contraction gives \(A\in k[K,C]\), and
the weight-six part is \(\langle K^3,KC\rangle\).

**Theorem SexticKeller644ReducedBinaryInversion.**
A characteristic-zero three-variable Keller map of ordinary component
degrees \((6,4,4)\) is invertible if the highest tuple is binary and
the two quartic highest forms are coprime over the algebraic closure.

## Inversion after the square-member reduction

Choose two squarefree members of the coprime quartic pencil
by the same [generic member argument](SEXTIC_KELLER_BINARY_REDUCED_LOWER_PENCILS.md#binary-655)
as for the quintic pencil. Apply the normalization reduction above
to write the highest tuple as \((aK^3+bKC,K^2,C)\). The binary ring \(k[K,C]\) is finite
over the image ring, and it generates the image field, so both
generators lift at the first new weight. The same holds in any
subsequent binary frame by transport of integrality and of the
image field. Write \(H\) for the uplift of \(C\).
The leading map keeps weights \((6,4,4)\) throughout the proof.
In these source frames the sum of the weights at a first tie
is at most ten, so the leading determinant vanishes.
Also every weight used is at least one. Thus any Keller map
in the weighted frames below has ordinary degree at most six.

If at a first tie the quadratic generator acquires a term in
the newly introduced variable, that variable has weight two
and the new term is a nonzero scalar times the variable.
Make the lifted generator a source coordinate of weight two.
This triangular change respects the weights. The target
component with highest form equal to the square of that
generator now has ordinary degree at most three: its
highest weighted part is a pure coordinate square, and the
strictly lower weighted parts have integer weights at most
three. At these ties the fixed weights are integers, as is
the moving weight two. The [cubic target theorem](SEXTIC_KELLER_CUBIC_COMPONENT.md)
therefore gives inversion, which descends to the original
map. We may assume the lifted quadratic generator always
equals its old expression.

Initially the source weights are \(1,1,w\) at the first tie.
A nontrivial lift \(H\) then requires
\(w\in\{4/3,3/2,2,3,4\}\) by its weight four.
If \(w\ge2\) and \(H\) is linear in the new variable,
the actual quartic component is affine in that variable:
all terms of exponent at least two would weigh at least four.
The [quartic affine-source inversion theorem](SEXTIC_KELLER_SQUAREFREE_CUBIC_QUADRATIC_MATES.md)
applies in these original linear source coordinates.

The remaining cases use the
[weighted rational slice](KELLER_LEADING_PAIR_RATIONAL_SLICE.md)
for \(K^2,H\). Dividing by the factor \(2K\) on the
generic level gives a rational primitive for the relative
volume differential of \(K,H\). On a generic surface
\(K=s\), this is the residue differential of the plane
area written in coordinates on that surface, with its
coordinate Jacobian factor.

First suppose \(K=xy\), after linear changes and scaling.
On \(K=s\) write \(y=s/x\). The area for the \(H\) level is
\(dx\wedge dz/x\), up to sign, since the derivative of
\(xy\) in \(y\) is \(x\). The polynomial \(C(x,s/x)\)
is a Laurent polynomial with poles of order four at both
\(x=0,\infty\) because \(K,C\) are coprime binary forms.
At \(w=4/3\), \(H=C+c z^3\) with constant \(c\ne0\).
The relative differential is proportional to \(dx/(xz^2)\).
The generic zeroes of \(H-C\) viewed as a Laurent polynomial
on the fixed level of \(H\) are simple in \(\mathbb G_m\);
the differential is regular above them by ramification of
index three. At \(0,\infty\) the pole of \(z\) makes this
form regular. The form is therefore nonzero and regular
on the complete normalization, contrary to exactness.

At \(w=3/2\), \(H=C+a_1(x,y)z^2\) with nonzero
\(a_1\) homogeneous linear. With \(M=a_1(x,s/x)z\) the
level equation becomes \(M^2=a_1(x,s/x)(S-C(x,s/x))\).
Here \(S\) is the generic \(H\)-value. The right hand
side has only simple zeros on \(\mathbb G_m\): the zeros
of the nonzero function \(a_1(x,s/x)\) there are simple
and disjoint from the generic level of \(C\). At \(0,\infty\)
the right hand side has pole order at least three, since
\(a_1(x,s/x)\) has order at most one there. Hence the
relative differential \(dx/(xM)\) is again nonzero and
regular everywhere on the complete curve.

At \(w=2\) with a square term, complete the square by a
weight-preserving shift in \(z\) to write \(H=c Z^2+f(x,y)\),
\(c\ne0\). Now the relative differential is proportional
to \(dx/(xZ)\). At each of \(0,\infty\) where \(f(x,s/x)\)
has a pole the differential is regular, while at an end
where \(f(x,s/x)\) is finite it has a simple pole,
since on a generic level \(Z\) is a unit there.
All other points are regular, including simple branch
points, and the constant case for the Laurent
polynomial has the same simple poles on each component.
Either behaviour precludes exactness.

Now suppose \(K=x^2\) up to scalar. Fix a generic nonzero
\(x\). The slice is the ordinary relative area differential
on the geometric generic plane level of \(H\).
The quartic \(C\) dehomogenized at \(x=1\) has degree
four in \(y\), by coprimality. At the weights \(4/3,3/2\)
we get respectively a cubic in \(z\), or a square in
\(z\) with affine \(y\)-coefficient, added to that quartic.
The [relative area calculations](SEXTIC_KELLER_BINARY_REDUCED_LOWER_PENCILS.md#binary-654)
give an everywhere regular nonzero form, using
\(dy/z^2\) in the first case and \(dy/(a_1 z)\) in the
second. At weight two with a square term, complete its
square. At \(x=1\) the result is \(cZ^2+f(y)\).
The hyperelliptic differential is regular if \(\deg f\ge3\)
and has simple poles if \(\deg f=2\), so
\[
                 H=cZ^2+\gamma x^3y+\delta x^4 .
 \tag{1}
\]
If \(\gamma=0\), increase the weight of \(y\) from one,
with \(x,Z\) of weights one and two. If the quadratic
generator changes the cubic target argument applies.
Otherwise the possible first ties in \(H\) are
\(4/3,3/2,2,3,4\) by the weight of a monomial containing
\(y\). At \(4/3,3/2\) the lifted \(H\) adds a nonzero
cubic or quadratic in \(y\) to \(c Z^2+\delta x^4\)
at \(x=1\), so exactness fails by the same calculation.
At weight two, complete the square in \(Z\), absorbing
any term in \(yZ\). Relative area exactness again
makes the residual degree in \(y\) at most one. If it
is zero, the tuple is again binary and a further tie
must have weight greater than two. The only such
ties have weight three or four, where the only
possible contribution is a linear term in \(y\) with
coefficient a power of \(x\). We are thus left with
\[
 H=cZ^2+\gamma x^{4-j}y+\delta x^4,\qquad
 c\gamma\ne0,\qquad (w_x,w_y,w_Z)=(1,j,2),
 \quad j\in\{1,2,3,4\}.
 \tag{2}
\]
The case \(j=1\) is (1) with nonzero \(\gamma\).

In (2) the actual quartic component \(F_4\) has a nonzero
constant coefficient of \(Z^2\). Its lower-weight coefficient
of \(Z\) has weight below two. Complete the square in the
entire component by a weight-preserving triangular change,
giving \(c Z^2+g(x,y)\). The polynomial \(g\) has no critical
point because the Keller condition makes \(F_4\) a submersion.
If \(j\ge2\), the leading form in (2) and the bound on
all smaller weights show that \(g\) is affine in \(y\).
For \(j=1\) it has the form
\[
 g=a y^3+b(x)y^2+d(x)y+e(x),\quad
 \deg b\le1,\quad \deg d=3,\quad \deg e\le4 .
\]
If \(a\ne0\), depress the cubic by an affine shift in \(y\).
The degree of \(d\) remains three and that of \(e\) remains
at most four. Then the critical equations are
\(3a y^2+d=0,\ d' y+e'=0\).
The polynomial \(3a(e')^2+d(d')^2\) has degree seven,
so has a root. If \(d'\ne0\) there, substitution of
\(y=-e'/d'\) gives a critical point. If \(d'=0\) there,
also \(e'=0\), so solving the first equation still gives
a critical point. Thus \(a=0\).

If \(b\) is nonzero constant, complete the square in
\(y\), reducing the critical equations to the derivative
of \(e-d^2/(4b)\), a polynomial of degree six, and a
linear equation for the shifted variable. This gives a
critical point. If \(\deg b=1\), translate \(x\) so that
\(b=\beta x\). When \(d(0)\ne0\), the expression
\(x^2(e-d^2/(4\beta x))'\) is a polynomial of degree six
with nonzero constant. A root \(x\ne0\) gives a critical
point by \(y=-d/(2\beta x)\). When \(d(0)=0\), completing
the square in \(y\) is a polynomial change and gives
\(\beta xY^2+E(x)\) with \(E=e-d^2/(4\beta x)\) of degree
five. If \(E'\) has a nonzero root it gives a critical
point at \(Y=0\). Otherwise \(E'(0)=0\), and \(x=Y=0\)
is critical. These contradictions give \(b=0\), so \(g\)
is affine in \(y\) also at \(j=1\).

The triangular changes used to obtain the actual
\(c Z^2+g(x,y)\) preserve the weight bounds. Thus the
resulting map has ordinary degree at most six and its
displayed component is either quartic and affine in \(y\)
of degree exactly one (by \(\gamma\ne0\)), or has degree
at most three. The quartic affine-source or cubic target
inversion theorem applies. This finishes the binary profile.

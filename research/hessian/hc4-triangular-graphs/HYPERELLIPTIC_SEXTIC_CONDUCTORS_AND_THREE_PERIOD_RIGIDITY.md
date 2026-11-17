# Hyperelliptic sextic conductors and three-period rigidity

Let K=(P,Q) be an actual complex polynomial plane Keller map in a
source frame where P is monic of positive degree in x. Let kappa
generate the exact coordinate-order conductor and suppose

    deg_x kappa=6.

The generic source-line image C[P(x,c),Q(x,c)] CANNOT belong to
Neumann's hyperelliptic genus-three ambient class.

This excludes one of the THREE realizable genus-three classes in
the primary low-genus table. The intrinsic-degree-four and
intrinsic-degree-six cable classes remain open. No conclusion
that all degree-six conductors are absent is asserted.

The hypothesis concerns the actual plane embedding on a nonempty
Zariski open set of source parameters, up to polynomial target
automorphisms. An abstract curve isomorphism or an unproved
classification of all nodal embeddings is not substituted for it.

## 1. The exact genus-three source input

Neumann, *Complex algebraic plane curves via their links at infinity*,
Invent. Math.98(1989), Section8, printed pages482-483, gives exactly
three realizable genus-three defining-polynomial types, up to
polynomial ambient automorphisms and affine changes of polynomial
value:

    degree4:
      X^4+Y^3+Y*(a*X^2+b*X+c)+d*X^2+e*X;

    degree7:
      X^7+Y^2+a*X^5+b*X^4+c*X^3+d*X^2+e*X;

    degree6 cable:
      (X^3+Y^2)^2+a*X^4+b*X^3+(a^2/4)*X^2+c*X
                       +a*X*Y^2+b*Y^2+d*Y, d!=0.         (1)

The other genus-three knot row in the table is not realizable.
The exact page images, source hashes and applicability are retained
in the [genus-three source audit](sources/neumann-genus-three-source-audit.json).
The table is an ambient polynomial classification, not merely
a topological list of knots.

The [coordinate-order theorem](GENERIC_COORDINATE_ORDER_CONDUCTOR_AND_CLOSED_COLLISIONS.md)
gives exact specialized conductors and normalization C[x] on EVERY
source line, with

    deg_x kappa_c=2*delta(c)>0,

and immersion of that normalization. Thus the generic defect here
is delta=3. The reduced one-place-at-infinity fibre has a good
defining polynomial by corrected Neumann-Rudolph Lemma7.1.
Smoothing its singularities gives a once-punctured genus-three
fibre, as in the retained low-genus applications.

The degree-seven alternative in(1) is the hypothesis being excluded.
Its rational normalized fibre has polynomial normalization

    X=tau^2+gamma,
    Y=tau*f(tau^2), deg f=3, f(0)!=0.                    (2)

Indeed the irreducible equation is quadratic in Y, so
[C(tau):C(X)]=2. The polynomial X therefore has degree two
in the affine-line normalization parameter, and Y has degree seven.
Complete the square in the quadratic X.
The equality Y(tau)^2=G(X(tau)) makes Y odd in tau, and immersion
at tau=0 forces f(0)!=0. Nonzero target scalings absorb leading
coefficients.

For each such fibre its exact order and conductor are

    A=C[tau^2] direct-sum tau*f(tau^2)*C[tau^2],
    conductor=f(tau^2)*C[tau].                           (3)

Indeed the monogenic equation is Y^2=X*f(X)^2 after translating
X by gamma, and its conductor ratio is F_Y/(2tau)=f(tau^2).
Multiple roots of f represent contact of reflected branches.

The order description descends to C(y) without descending a family of
target automorphisms. The midpoint is rational, being the conductor
root centroid. After this rational centering, divide the odd parts of
the actual P,Q by u*kappa over C(y)[u]. Their polynomial remainders
specialize to zero on the nonempty parameter open where the ambient
hypothesis holds. Thus both remainders vanish identically. Consequently

    A_gen subset C(y)[u^2,u*kappa].

Both orders have colength THREE in C(y)[u]: the left side by the exact
generic conductor, and the right side by the even/odd decomposition.
They are therefore equal. This proves the generic order identity used
below by coefficient descent and equal colength, without choosing an
algebraic family of target changes.

## 2. Persistent hyperelliptic contact contradicts the actual Jacobian

Over the generic parameter field translate by the rational
normalization midpoint, writing the centered parameter u.
The exact order is C[u^2,u*f(u^2)]. Suppose f has a repeated root
r^2 in an algebraic extension. Here r!=0 by immersion.

The odd part of every order element is divisible by u*f(u^2).
Its derivative therefore vanishes at u=r and u=-r. For the two
actual target functions this gives

    (P_u,Q_u)(-r)=-(P_u,Q_u)(r).                         (4)

Their endpoint values agree as algebraic functions of y.
Differentiate P(r(y),y)=P(-r(y),y), and use(4). The r' terms
cancel, leaving P_y(r)=P_y(-r). The same holds for Q:

    (P_y,Q_y)(-r)=(P_y,Q_y)(r).                          (5)

Consequently the source Jacobian satisfies J(-r)=-J(r), contrary
to its being one nonzero constant. A rational source shear leaves
that Jacobian unchanged, and all these derivatives are legitimate
in the separable algebraic function field.

Thus f is GENERICALLY SQUAREFREE. Every generic image in the
hypothesized class has exactly three ordinary nodes, paired by
the reflection about its common midpoint.

## 3. Monicity repairs the centroid through every finite degeneration

Let a_6,a_5 be the two highest coefficients of kappa in x. Its
generic root centroid is

    h(y)=-a_5(y)/(6*a_6(y)).

The reflection x -> 2h-x takes a generic conductor point to its
actual same-(P,Q,y) mate.

At a finite parameter c where a_6 vanishes, the specialized conductor
is nonzero of degree two or four. Choose ANY finite conductor root
x_0. It need not be simple. Normalize a nonvertical branch of the
conductor curve through (x_0,c), and clear its parameter unit to write

    y=c+s^alpha, x=x(s) in C[[s]], x(0)=x_0.

Such a branch exists: a vertical component would make kappa_c
identically zero, which the exact source theorem forbids.
This is the elementary curve-normalization/Puiseux construction.

Its generic reflected mate z(s)=2h(c+s^alpha)-x(s) satisfies
P(z(s),c+s^alpha)=P(x(s),c+s^alpha). Since P is MONIC in x,
a negative order of z would make its leading monomial dominate
all lower regular coefficients. Thus z is regular, and h is regular
at c. This repairs every finite possible pole of h, so h is polynomial.

Make the single polynomial shear u=x-h(y). The generic reflection
symmetry makes the entire conductor even:

    kappa=a(y)*u^6+b(y)*u^4+c(y)*u^2+d(y).                (6)

Its constant coefficient d has no zero. Otherwise choose a finite
conductor branch u(s) tending to zero at that parameter. It is
nonzero generically, because f(0)!=0. The reflected mate identities
give

    [P(u(s),y)-P(-u(s),y)]/(2u(s))=0,

and the same for Q. Their limits force P_u=Q_u=0 at u=0,
contradicting the actual Keller Jacobian. Hence

    d in C*, a,b,c in C[y], gcd(a,b,c)=1.                 (7)

The gcd assertion follows because no conductor fibre is constant.
This argument retains repeated special roots and all degree drops.
It uses no classification of those special fibres.

## 4. Three periods give three exact differential equations

Set T=u^2 and g(T)=a*T^3+b*T^2+c*T+d. Let A_i, i=1,2,3,
be its distinct generic roots and choose r_i^2=A_i. A canonical
target pair for the generic fibre is

    p_0=u^2, q_0=u*g(u^2).

Its period at the i-th node is

    integral_(-r_i)^(r_i) p_0*dq_0
       =(8/315)*r_i*C(A_i),
    C(T)=5*a*T^4-9*c*T^2-30*d*T.                         (8)

The three canonical periods do not all vanish. If they did,
the cubic C(T)/T would be a scalar multiple of g(T), since
all A_i are distinct and nonzero. Leading terms force the scalar
to be five, but constant terms would give -30d=5d, impossible.

Polynomial target changes multiply all these periods by their
same constant nonzero Jacobian. Thus the actual periods have
a common algebraic scale Lambda, absorbing the factor 8/315:

    I_i=integral_(-r_i)^(r_i) P*Q_u*du
       =Lambda*r_i*C(A_i).                              (9)

Choose one nonzero canonical period to define Lambda on an open
set; pointwise target equivalence proves all three identities.
No uniform family target automorphism is assumed.

The exact moving-endpoint identity gives I_i'=-2J*r_i, where
J=P_u*Q_y-P_y*Q_u is a nonzero constant. Differentiating(8), or
integrating the canonical moving-endpoint identity, gives modulo g

    d_y C(A_i)+(A_i'/(2A_i))*C(A_i)
       =-(35*a'*A_i^4+45*b'*A_i^3+63*c'*A_i^2)/2.        (10)

Here d'=0 by(7). Consequently the polynomial

    M(T)=2*Lambda'*C(T)
           -Lambda*(35a'*T^4+45b'*T^3+63c'*T^2)+4J

vanishes at all three roots of g and is divisible by g.
Write its degree-one quotient as alpha*T+beta. Comparison of
the leading and constant coefficients gives

    alpha=10Lambda'-35(a'/a)Lambda, beta=4J/d.

Comparison of the other three coefficients yields EXACTLY

    70d^2*a*Lambda'-35d^2*a'*Lambda=-4J*a*c,              (11)

    10a*b*Lambda'+(45a*b'-35b*a')Lambda=-4J*a^2/d,        (12)

    28a*c*Lambda'+(63a*c'-35c*a')Lambda=-4J*a*b/d.        (13)

All coefficient degrees and the degrees of P,Q are unrestricted.

## 5. A zero of a is impossible

Suppose a has a zero of order m>=1.

First let b be a unit there. Put n=ord Lambda in the algebraic
function field, allowing rational orders. In(12) the derivative
terms have leading coefficient proportional to 10n-35m and
possible order m+n-1, whereas the right side has order 2m.

The only cancellation value is n=7m/2, which is greater than
m+1 for every positive integer m and would put the whole left
side above order 2m. Thus no cancellation changes the balance:
equation(12) forces

    n=m+1.

But the left side of(13) then has order at least 2m+ord(c),
while its right side has exact order m because b is a unit.
This is impossible, including c=0 identically.

Otherwise b vanishes, so c is a unit by(7). Eliminate Lambda'
between(11),(13), obtaining

    105d^2*(3a*c'-c*a')Lambda
       =4J*a*(2c^2-5b*d).                               (14)

Its two polynomial factors have orders m-1 and m, respectively.
It therefore forces ord Lambda=1. The leading coefficient of
(13), at order m, is proportional to

    28-35m,

whereas its right side has strictly greater order because b
vanishes. This would give m=4/5, again impossible.

Hence a has no finite zero and is a nonzero CONSTANT.

## 6. Constant a leaves only an impossible isotrivial limit

Equation(11) now gives

    Lambda'=-2J*c/(35d^2).

The derivative of Lambda is polynomial. Its polynomial primitive
differs from Lambda by a constant in a finite algebraic extension
of C(y); that constant belongs to C. This follows by differentiating
a monic minimal polynomial, as in the retained constant-field lemma.
Thus Lambda itself is polynomial.

If Lambda is constant, equations(11),(13) give c=b=0, contradicting
(12). Let n=deg Lambda>=1. Eliminating c and b gives

    c=-35d^2*Lambda'/(2J),
    b=35d^3*[28(Lambda')^2+63Lambda*Lambda'']/(8J^2).

For leading coefficient lambda_n, b has degree 2n-2 and nonzero
leading coefficient

    245d^3*n*(13n-9)*lambda_n^2/(8J^2).

The left side of(12), after division by constant a, has degree
3n-3 and leading coefficient

    1225d^3*n*(10n-9)*(13n-9)*lambda_n^3/(4J^2).

This is nonzero for every positive integer n. Its right side is
constant, so n>=2 is impossible. The only remaining case is n=1,
which makes a,b,c,d all CONSTANT and Lambda affine nonconstant.

The conductor now has six fixed simple roots, because it was
generically squarefree. Generic matching identities extend to every
parameter. Exact specialized defect three identifies EVERY fibre
order with the SAME central three-node order

    A_0=C[u^2] direct-sum u*g(u^2)*C[u^2].                (15)

Every polynomial of A_0 of degree at most SIX is even.

No special fibre can switch to either other genus-three ambient
class in(1). A birational polynomial parametrization of a plane
curve of degree D has maximum coordinate degree D: a general affine
line pulls back to a polynomial of that degree. A polynomial target
change to an intrinsic-degree-four or degree-six model would give
two elements of A_0 of degree at most six. Both would be even and
could not generate its fraction field C(u), a contradiction.

Thus every special fibre remains in the hyperelliptic class.
At least one fixed canonical period in(8) is nonzero and the
corresponding actual period is nonzero on EVERY fibre. Yet(9)
has the nonconstant affine factor Lambda, which has a complex zero.
This final contradiction excludes the entire hypothesized class.

## 7. Exact scope and evidence

The checker pins the primary genus-three source audit and PDFs.
It verifies the period and derivative formulas, polynomial quotient
identity, all three differential equations, both local zero-order
obstructions, and the polynomial-degree obstruction. It includes a
persistent-contact control showing the opposite-Jacobian mechanism.

The result closes the GENERIC HYPERELLIPTIC ambient branch, including
its possible persistent contacts and special degree drops. It does
not exclude genus-three intrinsic quartic or cable embeddings.
No claim that all degree-six conductors or all higher-genus
embeddings are hyperelliptic is made.

The genus-three classification remains the attributed primary-source
input. Finite symbolic controls support the written unrestricted
argument; they do not prove the classification. Independent replay,
formal verification and external human review are not claimed.

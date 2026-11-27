# The leading fields of the order-\((12,42)\) Weyl row

Work over an algebraically closed field \(k\) of characteristic
zero and in the remaining reduced differential orders
\((12,42)\). The operators are normalized as in the
[order-thirty-four resultant](POLYNOMIAL_WEYL_ORDER_TWELVE_FORTY_TWO_RESULTANT.md).
Write \([Q,P]=\kappa\in k^\times\), with leading differential
coefficients \(h^2,h^7\). Set \(\phi=h^{1/3}\in k[x]\).
Its geometric root multiplicities are one and five,
and there is at least one root of multiplicity five.
For the resultant \(W\), choose \(\delta\in k^\times\) so
that its leading coefficient is \(\delta\phi^{17}\).
Let \(B\) be the polynomial solution of
\[
                         B'-\frac{h'}{6h}B=1
\]
from the [residue calculation](POLYNOMIAL_WEYL_LEADING_FACTOR_RESIDUE.md).

## The second resultant

**Theorem PolynomialWeylFortyTwoSecondResultant.**
There is a unique constant \(e\) such that
\[
                 Z=W^3-\delta^3 P^5Q-eWP^2Q
\]
has differential order \(99\). Its leading coefficient is
\[
                         -\frac{\kappa\delta^2}{2} B\phi^{49}.
\]
In particular this result applies also when \(h\) has no
multiplicity-three roots.

Put \(Z_1=W^3-\delta^3 P^5Q\). The two terms have
order \(102\) with the same leading coefficient.
The commutator of \(Z_1\) with \(P\) is
\[
 [Z_1,P]=2\kappa(W^2Q+WQW+QW^2)-\kappa\delta^3P^5 .
\]
It has order \(110\) and leading coefficient
\(6\kappa\delta^2\phi^{55}\).
If \(Z_1\) had order \(j=101\), its would-be
order-\(j+11\) commutator coefficient must vanish.
The leading coefficient would then be proportional
to \(h^{j/6}\), by the differential equation from the
leading symbols. This is impossible because the
multiplicity fifteen in \(h\) makes \(15j/6\)
nonintegral. If \(Z_1\) has order \(100\), the same
comparison makes its leading coefficient proportional
to \(\phi^{50}\). The term \(WP^2Q\) has order
\(100\) with leading coefficient \(\delta\phi^{50}\),
so a unique subtraction kills that order. If the
order is already at most \(99\), set \(e=0\).
The commutator of \(WP^2Q\) with \(P\) has order at
most \(108\), as its terms are \(2\kappa QP^2Q\)
and \(\kappa WP^2\). Thus \([Z,P]\) retains order
\(110\). The upper bound \(\operatorname{ord}[Z,P]\le
\operatorname{ord}Z+11\) makes \(\operatorname{ord}Z=99\).

In \(k(x,\lambda)\) with \(\lambda^2=\phi\), let \(z_{99}\)
be the leading coefficient of \(Z\). The bracket of the
leading differential symbols is
\[
  99 z_{99} (h^2)'-12h^2 z_{99}'
       =6\kappa\delta^2\lambda^{110}.
\]
Dividing out the factor \(\lambda^{111}\) gives
\[
                   (z_{99}/\lambda^{99})'
                        =-\frac{\kappa\delta^2}{2\lambda}.
\]
The primitive \(B/\lambda\) has derivative \(1/\lambda\).
The difference after subtracting its scalar multiple
is constant. The involution \(\lambda\mapsto-\lambda\)
changes the signs of both terms and fixes the constants,
so that difference vanishes. This gives the formula.

## The differential-order filtration

**Theorem PolynomialWeylFortyTwoDifferentialOrderAlgebra.**
For the algebra \(\mathcal A=k\langle P,Q\rangle\) with
the induced differential-order filtration, put
\[
                 \tau=\phi\xi^2,\qquad
                 \zeta=B\tau^{49}\xi .
\]
Inside \(k[x,\xi]\), its associated graded algebra is
\[
                       \operatorname{gr}\mathcal A
                          = k[\tau^6,\tau^{17},\tau^{21},\zeta].
\]
The elements \(\tau,\zeta\) are algebraically independent.
If \(n_H\) counts the geometric multiplicity-fifteen roots
of \(h\), the inclusion of graded fraction fields has degree
\[
                 [\,k(x,\xi):k(\tau,\zeta)\,]=3n_H .
\]
The exact semigroup of differential orders of nonzero elements
of \(\mathcal A\) is \(\langle12,34,42,99\rangle\).

First, the ordered words
\[
               P^i Q^j W^a Z^b,\qquad
                i,b\ge0,\quad 0\le j<2,\quad 0\le a<3
                                                               \tag{1}
\]
span \(\mathcal A\). In normal order with coefficients in
\(k[P]\), \(W\) is monic of degree two in \(Q\) and
\(Z\) is monic of degree six in \(Q\).
Moving \(P\) past \(Q\) lowers \(Q\)-degree in the
correction term. For every normal word \(P^i Q^n\), there
is a unique \(n=j+2a+6b\), and subtraction of (1) leaves
smaller \(Q\)-degree. This proves spanning by induction.

The leading symbols of (1) are nonzero multiples of
\(\tau^{6i+21j+17a}\zeta^b\).
The exponents \(21j+17a\) give all six distinct residues
modulo six. Hence these pairs of monomial exponents
are distinct. To see their independence, observe that
\[
                 \frac{\zeta^2}{\tau^{99}}=S(x),
                   \qquad S(x)=\frac{B^2}{\phi}.
\]
This rational function has a triple pole at each
multiplicity-fifteen root, because \(B\) has a simple
zero and \(\phi\) has a zero of order five there.
It is regular at each multiplicity-three root and at
infinity. At infinity \(\deg B=\deg h/6\).
It follows that \(\deg S=3n_H>0\).
Since \(\tau\) is transcendental over \(k(x)\), the
elements \(\tau,S\), and therefore \(\tau,\zeta\),
are algebraically independent. Thus leading terms of
distinct words at the same differential order cannot cancel.
This proves the graded formula, including the claimed
semigroup.

To calculate the field degree, \(k(x,\xi)=k(x,\tau,\zeta)\)
by \(\xi=\zeta/(B\tau^{49})\). The extension
\(k(x,\tau)/k(S,\tau)\) has degree \(3n_H\).
Adjoining \(\zeta\), a root of
\(U^2-S\tau^{99}\), is quadratic over both fields, since
the order in the independent variable \(\tau\) is odd.
The tower law gives the claim.

## The two high-root Newton faces

**Theorem PolynomialWeylFortyTwoHighRootLeadingFields.**
Translate a multiplicity-fifteen root to \(x=0\), and put
\(X=\xi,\ Y=x\) in symbols. For weights \((1,-s)\)
in \(X,Y\), the leading field of \(\mathcal A\) has the
following degrees in the indicated ambient rational field.
Here \(t^3=X\), and the shift is
\(Y_{\rm old}=Y_{\rm new}+\omega t^{-1}\) at the
multiplicity-four root of the first face.

| Presentation | Weight parameter \(s\) | Ambient field | Degree of ambient field over leading field |
|---|---:|---|---:|
| Unshifted | \(0<s<1/3\) | \(k(X,Y)\) | \(3\) |
| Unshifted | \(s=1/3\) | \(k(X,Y)\) | \(5\) |
| Unshifted | \(1/3<s<1\) | \(k(X,Y)\) | \(1\) |
| Shifted | \(1/3<s<4/7\) | \(k(t,Y_{\rm new})\) | \(5\) |
| Shifted | \(s=4/7\) | \(k(t,Y_{\rm new})\) | \(8\) |

The table is for rational \(s\). A leading field is the
fraction field of the associated graded algebra for that
weight. The degree-five crossing map has pole index five
at the shifted point and index three at its unshifted
infinite parameter. The degree-eight successor map has
index five at its infinite parameter.

At the first weight, use the integral multiple
\((-1,3)\) in \(x,\xi\). As in the first resultant proof,
the leading symbols of \(P,Q,W\) are
\[
             R^2,\quad R^7,\quad V,\qquad
             R=x^{-3}r(z),\quad V=x^{-17}w(z),
                    \quad z=x^3\xi .
\]
The scalings of \(R^2,R^7\) can be made as displayed
because their highest differential coefficients are
the respective powers of the same leading coefficient
of \(h\) at \(x=0\).
The polynomial \(r\) has degree six with multiplicities
\(1,4,1\), including a simple zero at \(z=0\).
We have \(\deg w=34\) and
\[
             3 w'r-17 w r'=c_0r^6,\qquad
             w=r^5 a,\quad \deg a=4,\quad c_0\ne0.
\]
Each root of \(r\) is a simple zero of \(a\).
Thus \(V^3/R^{17}=a^3/r^2=:S_1(z)\) has one pole
of order five, at the quadruple root of \(r\). Its
logarithmic derivative is \(c_0/a\), so it is unramified
off its divisor except at infinity, where \(S_1\)
has finite nonzero value and index three.
That value is \(\delta^3\), by the highest differential
terms of \(R,V\).

The elements \(P^i Q^j W^a\) with \(i,a\ge0\), \(0\le j<2\)
also span \(\mathcal A\), by their distinct leading
\(Q\)-degrees. At the first weight their top symbols
are \(R^{2i+7j} V^a\), which are independent since
\(R,V\) are algebraically independent. Indeed with
\(U=V/R^6=xw/r^6\) we have \(S_1=U^3 R\),
so \(k(R,V)=k(U,S_1)\) and \(k(x,\xi)=k(U,z)\).
The degree is \(\deg S_1=5\).
It follows also that the symbol of \(Z\) at this weight is
\(R^{17}(S_1-\delta^3)\), of weight \(51\) and
differential order \(99\) since \(S_1-\delta^3\)
has order three in \(1/z\).

At the four differential-order leading terms of
\(P,Q,W,Z\), the orders of vanishing at \(x=0\)
are \(30,105,85,246\), respectively.
For \(Z\) use the order-\(99\) leading coefficient
and the simple zero of \(B\).
The resulting monomials attain the maximum at both
\(s=0\) and \(s=1/3\).
Hence they are the unique maximizing monomials for
\(0<s<1/3\). These are powers of \(Y^5X^2\) for
\(P,Q,W\) and the monomial \(Y^{246}X^{99}\) for \(Z\).
The exponent determinant for these last two generators
of the leading field has absolute value three.
Their independence makes (1) a basis with independent
leading symbols of equal weight, by the same residue
comparison modulo six. This proves the degree three.

On the other side without shifting, the symbols of
\(P,Q\) immediately become multiples of \(X^2,X^7\),
and stay so because they use no \(Y\).
The first symbol of \(W\) has lowest \(Y\)-exponent
one, there with term \(YX^6\), since \(w\) has
multiplicity six at \(z=0\). Any monomial of \(W\)
with \(Y\)-degree zero has \(X\)-degree at most five
by the first weight bound. Thus for \(1/3<s<1\)
the top of \(W\) is a multiple of \(YX^6\).
The spanning words \(P^iQ^jW^a\) again have independent
top symbols, and the leading field is \(k(X,Y)\).

Now set \(t^3=X\) and shift \(Y_{\rm old}=Y+\omega t^{-1}\)
with \(\omega^3\) the multiplicity-four root of \(r\).
The shift is in the ramified Weyl ring and preserves
the scalar commutator. The cubic substitution on the
first direction parameter is unramified at this nonzero
root and preserves the pole index five. The
[scalar successor calculation](POLYNOMIAL_WEYL_BERNSTEIN_FORTY_TWO_REAL_BOUND.md#testing-the-scalar-successors)
applies at this finite crossing. Its input is the local
\((6,15)\) edge and ratio \(2:7\), as used in the
[order-twelve crossing argument](POLYNOMIAL_WEYL_ORDER_TWELVE_RESIDUAL.md#the-crossings-at-gcd-six).
Changes of \(Q\) by a polynomial in \(P\) of degree
at most three are below the highest weight on this
positive-weight interval.
The successor at \(s=4/7\) has \(P,Q\) weights
\(2/21,7/21\) and scalar leading bracket.
At the lower end of their shifted first faces the
two symbols are powers of \(t^7Y^4\).
For \(W\) the lower end there is \(t^{38}Y^{21}\),
since \(w\) has order \(21\) at \(\omega^3\).
Indeed in the shifted first symbols \(z=(tY+\omega)^3\),
and the factor \((tY+\omega)^{-17}\) multiplying \(w(z)\)
is a unit at \(tY=0\).
At the successor the weight of \(W\) is \(14/21\)
from \(Q^2-P A(P)^2+2c A(P)\).
The leading \(Q^2\) and \(P^7\) symbols cannot cancel
because the leading \(P,Q\) bracket is nonzero.
The endpoint monomials attain the bounds at both ends
of the interval. By strict interpolation these are the
unique tops inside the interval. The spanning basis with
\(0\le j<2\) again has independent top symbols.
The determinant of the exponent vectors \((7,4),(38,21)\)
in \(t,Y\) is \(-5\), giving the shifted degree five.

At \(s=4/7\), use weights \((7,-12)\) in \(t,Y\).
The top forms have weights two and seven.
As the penalty for \(Y\) increases, the maximizing
\(Y\)-degrees cannot increase. The endpoints from the
shifted first faces thus bound those degrees by eight
and twenty-eight.
Writing \(v=t^{12}Y^7\), the lattice of weights gives,
after nonzero rescalings,
\[
                       p=t^2Y(A_0+v),\qquad q=t G(v),
                       \qquad \deg G=4.
\]
The scalar bracket, computed in \(t,Y\) with the factor
\(1/(3t^2)\) for differentiation in \(X\), implies
\[
                 (A_0+8v)G-2v(A_0+v)G'=c_1,\qquad c_1\ne0 .
\]
Therefore \(A_0 G(0)G(-A_0)\ne0\) and \(G\) has simple
roots. Since \(p,q\) are algebraically independent,
the leading field is \(k(p,q)\) by the ordered
\(P,Q\) basis. Put
\[
               S_2=\frac{p^7}{q^2}=\frac{v(A_0+v)^7}{G(v)^2},
                       \qquad U_2=q/p^3 .
\]
The last two functions generate \(k(p,q)\).
Given \(v,p,q\), one recovers \(t=q/G(v)\) and then \(Y\),
so the ambient field is \(k(v,U_2)\).
The function \(S_2(v)\) has degree eight.
Its logarithmic derivative is \(c_1/[v(A_0+v)G(v)]\).
This has order six at infinity, where \(S_2\) is a
unit, giving index five. Hence the consecutive ramification
indices three and five agree with the adjacent leading
field degrees on the displayed path.

## Geometric monodromy on the high-root path

**Corollary PolynomialWeylFortyTwoLeadingFieldMonodromy.**
The maps \(S_1,S_2\) of the high-root crossing and its
scalar successor have geometric monodromy groups \(A_5\)
and \(A_8\), respectively. For the differential-order
map \(S=B^2/\phi\), the smallest profile without small
roots has geometric monodromy \(S_4\) on six sheets.

The passport of \(S_1\) is \((3,1,1),(5),(3,1,1)\)
over zero, its pole value, and its third value.
The extra triple zero comes from the fourth root of \(a\),
which is simple and lies off \(r=0\).
The transitive group is contained in \(A_5\) since all
branch cycles are even, and has order divisible by
fifteen. There is no subgroup of order fifteen in
\(A_5\): its Sylow five-subgroup would be normal and
would centralize an element of order three, producing
an element of order fifteen. A subgroup of order thirty
would be normal of index two and contain every
three-cycle, which generate \(A_5\).
Thus the group is \(A_5\).

The passport of \(S_2\) is \((7,1),(2,2,2,2),(5,1,1,1)\)
and all branch cycles are again even. A transitive
degree-eight group with a seven-cycle is primitive,
because the block containing its fixed point is invariant
under that cycle and contains either just that point
or the entire orbit of the other seven.
The five-cycle now gives \(A_8\) by
[Jordan's prime-cycle theorem, as stated in Jones, Theorem 1.1](https://arxiv.org/pdf/1209.5169v1).

With two high roots and no small roots, the
[unique leading-factor tree](WEIGHTED_SLICE_TREE_CLASSIFICATION.md#the-order-twelve-profile)
allows normalization \(h=(x^2-1)^{15}\). Then
\(B=x(x^2-1)(2x^2/3-1)\), giving
\[
                    S=\frac{x^2(2x^2/3-1)^2}{(x^2-1)^3}.
\]
The degree-six passport is \((2,2,2),(3,3),(4,1,1)\).
The three blocks of size two come from \(x\mapsto x^2\).
In the cubic quotient the corresponding branch cycles
are two transpositions and a three-cycle, giving \(S_3\).
The sign of each generator on six sheets matches its
sign on the three blocks. Thus the group in
\(C_2^3\rtimes S_3\) lies in the subgroup of order
twenty-four where the sum of the within-block flips
is the parity of the block permutation.
Squaring a four-cycle gives the flip on two blocks.
Its conjugates generate the full even-flip subgroup
of order four, so the group has order twenty-four.
Adding the all-ones flip vector when the block
permutation is odd identifies this group with
\(C_2^2\rtimes S_3\simeq S_4\), acting on the four even
vectors by translations and permutations of coordinates.

These are monodromy groups of the leading rational maps.
Over characteristic zero they can be computed after
embedding the finitely generated coefficient field
into \(\mathbb C\) and using geometric branch cycles.

## The small-root accessory term

**Corollary PolynomialWeylFortyTwoOrderNinetyNineIndicial.**
At each multiplicity-three root, use the integer signed
divisor \(D\) and its normalized gap ideal \(I\) from the
[small-root classification](POLYNOMIAL_WEYL_ORDER_TWELVE_FORTY_TWO_RESULTANT.md#a-third-indicial-shift-at-a-small-root).
Let \(M_j(E)\) be the monic rational function with finite
divisor \(A_jD\), where \(A_j=1+\mathcal S+\cdots+\mathcal
S^{j-1}\) translates by \(0,\ldots,j-1\).
The lowest Euler exponent of \(Z\) is \(-49\), and its
component is a nonzero multiple of
\[
                         x^{-49}M_{49}(E)(E-\theta).
\]
The same \(\theta\in k\) applies at all small roots.
If \(I\) is the full gap set of \(\langle6,17,21\rangle\)
at any root, then \(\theta=24\).
The additional lowest-Euler polynomiality condition
excludes that gap set when \(\theta\ne24\).

Translate a small root to zero and write
\(\phi=\gamma x+O(x^2)\), \(\gamma\ne0\).
The commuting lowest components of \(P,Q,W\) are
\[
 x^{-6}\gamma^6 M_6(E),\quad
 x^{-21}\gamma^{21}M_{21}(E),\quad
 x^{-17}\delta\gamma^{17}M_{17}(E).
\]
The product rule
\[
                        M_i(E-j) M_j(E)=M_{i+j}(E)
\]
holds by the divisor identity \(A_{i+j}=A_j+\mathcal S^j A_i\).
In particular the two exponent-\(-51\) components
of \(W^3\) and \(\delta^3 P^5Q\) cancel.
The Euler exponent of \(WP^2Q\) is \(-50\), so
\(Z\) starts at exponent at least \(-50\).

From the displayed commutator in the second resultant
proof, the lowest Euler component of \([Z,P]\)
has exponent \(-55\) and is
\(6\kappa\delta^2\gamma^{55}x^{-55}M_{55}(E)\).
Indeed the terms from \([WP^2Q,P]\) start at
exponents \(-54\) and \(-29\), and \(P^5\) starts
at \(-30\). If the exponent-\(-50\) component of
\(Z\) were \(x^{-50}z_0(E)\ne0\), it would commute
with \(x^{-6}M_6(E)\). Comparing one degree below
the top in
\(z_0(E-6)M_6(E)=M_6(E-50)z_0(E)\)
would give \(6\deg z_0=600\), above the order bound
\(99\). Hence \(Z\) starts at \(-49\).
Writing that component as \(x^{-49}M_{49}(E)L(E)\)
and using the product rule gives
\[
                      L(E-6)-L(E)=6\kappa\delta^2\gamma^{49}.
\]
A rational periodic function in characteristic zero is
constant, since a pole and all its translates would be
poles, and a nonconstant polynomial is not periodic.
Thus \(L=-\kappa\delta^2\gamma^{49}(E-\theta)\).

To see that \(\theta\) is common, use symmetric
subprincipal coefficients. For a differential operator
of order \(d\) with first two coefficients \(a,b\), this
coefficient is \(b-(d/2)a'\).
Its use in a commutator is the sum of the two Poisson
terms with one leading and one subprincipal symbol.
This follows by expanding the two ordered products
through two derivatives, or by symmetric ordering as in
the [small-root subprincipal calculation](POLYNOMIAL_WEYL_SMALL_ROOT_SPECTRAL_CONSTRAINTS.md#common-spectral-phases-and-the-opposite-fourier-coefficient).
Replace \(D_x\) by \(D_x-S_P\), with
\[
                 S_P=\frac{p_{11}-6(h^2)'}{12h^2} ,
\]
in the rational differential operator ring. This makes
the symmetric subprincipal of \(P\) zero. By the
[small-root mean theorem](POLYNOMIAL_WEYL_SMALL_ROOT_SPECTRAL_CONSTRAINTS.md#exact-subprincipal-coefficients-at-small-roots),
\(S_P\) is regular at every small root, so the
lowest Euler components of the polynomial operators
there are unchanged.

In this gauge the subprincipal of \(Q\) is proportional
to \(\lambda^{41}\), because \([Q,P]\) has order below
\(42+12-2\). Similarly that of \(W\) is proportional to
\(\lambda^{33}\), because \([W,P]=2\kappa Q\) has
order below \(34+12-2\). In each case the equation
for the subprincipal \(b\) of an operator of order \(d\)
is \(b'=(d-1)(h^2)'b/(12h^2)\),
so these claims follow by direct integration.
Both coefficients vanish by the involution of \(\lambda\).
Consequently the order-\(110\) operator \([Z,P]\) has
zero subprincipal, since its three highest products have
zero subprincipal and every other summand has order
at most \(108\). In a symmetric product the first
correction is a bracket of the leading symbols, and
these symbols are all powers of \(\tau\).
The same subprincipal equation for \(Z\) now gives
\(C\lambda^{98}=C\phi^{49}\), with a single constant \(C\).

At the small root \(z_{99}\) has initial term
\(\eta x^{50}\) with \(\eta=-\kappa\delta^2\gamma^{49}\).
The coefficient of \(x^{49}D_x^{98}\) in the lowest
Euler component is thus \(2475\eta+C\gamma^{49}\),
since \(99\cdot50/2=2475\). Passing to \(E\) replaces
\(x^{99}D_x^{99}\) by \((E)_{99}\), whose second
coefficient is \(-4851\). The root sum of the resulting
monic degree-\(99\) Euler polynomial is
\(2376-C\gamma^{49}/\eta\).
The root sum from \(M_{49}(E)(E-\theta)\) is
\(49\cdot48+\theta\), since \(D\) has mass two and first
moment zero. Hence \(\theta=24+C/(\kappa\delta^2)\).

Finally the falling-factorial condition for \(Z\) requires
\(A_{49}(D-[0])+[\theta]\) to be effective.
In the gap-set description of \(D\), its first term
has multiplicities \(C_j-C_{j-49}\) at \(j-|I|\).
Since the conductor is \(50\), the only possible
negative entry occurs at \(j=49\) when \(49\in I\).
The semigroup has twenty-five gaps, and for each pair
\(\{j,49-j\}\) at least one is a gap. Thus each pair
contains exactly one gap. Every gap \(g\) can be reached
from \(49\) by subtracting a sequence of the generators
adding to \(49-g\). Downward closure shows \(49\in I\)
exactly for the full gap set. In that case the unique
negative entry is at \(49-25=24\), proving the claim.

The [planar transfer and generic-fibre calculation](POLYNOMIAL_WEYL_ORDER_TWELVE_FORTY_TWO_PLANAR_FIBRES.md)
carry these high-root caps to a planar pair. They determine
its quadratic charts, its high-root pole profiles and the
monodromy of its generic position projection.

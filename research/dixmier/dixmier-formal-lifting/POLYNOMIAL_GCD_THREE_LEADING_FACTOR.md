# Common leading factors for polynomial Weyl pairs of differential-order gcd three

The [Newton-edge argument](#exclusion-in-all-orders) below excludes all
nondividing differential orders of gcd three. The local factor and degree
restrictions give intermediate tests that also explain how the residue and
Bernstein results apply to this problem.

**Theorem PolynomialGcdThreeLeadingFactor.** Let \(k\) have characteristic
zero. If polynomial differential operators of orders \(3a,3b\), for
coprime \(a,b\ge2\), have nonzero scalar commutator, their leading
coefficients have the form \(\alpha h^a,\beta h^b\), where
\(\deg h\ge6\). Every root of \(h\) has multiplicity one, two,
or a multiple of three at least six. If \(h\) has more than one
distinct root then its degree is a multiple of three.

The [leading-factor residue theorem](POLYNOMIAL_WEYL_LEADING_FACTOR_RESIDUE.md)
gives the common factor, rules out exact triple roots and makes
the degree a multiple of three at least six if there are at least
two distinct roots. We rule out constant \(h\) and a single root
of multiplicity one or two using the leading cubic Hamiltonian.
The same Hamiltonian equation at a higher root restricts its
multiplicity.

## Normalization and weighted descent

Interchange the operators as needed so \(a<b\). Extend the field of constants to an algebraic closure. For constant
\(h\), put \(t=x,\ E=D,\ s=1\). For a single root \(x_0\) of
multiplicity \(m=1\) or \(2\), set \(x=x_0+t^3\) and
\(E=3t^{m}D\). In both cases \([E,t]=t^{1-s}\), where \(s=3-m\)
in the nonconstant cases. Scaling by constants makes the pair
monic in \(E\), with coefficients in \(k[t]\) in the constant
case and in \(k[t,t^{-1}]\) otherwise. Shift \(E\) by a coefficient
function to depress the first operator. Write the resulting
operators again as \(P,Q\) and their scalar commutator as \(\kappa\).

In the nonconstant cases the action
\[
                         (t,E)\longmapsto(\zeta t,\zeta^m E),
                         \qquad \zeta^3=1 ,
\]
fixes the original polynomial operators in these variables. It
commutes with depression: the coefficient of \(E^{3a-1}\) transforms
as \(E\), since \((3a-1)m+m\equiv0\pmod3\). Thus the depressed operators
are invariant.

The formal monic \(3a\)-th root \(L\) of \(P\) gives the finite
decomposition of \(Q\) into constant multiples of \((L^j)_+\).
Indeed their commutators with \(P\) have order at most \(3a-2\).
An order-\(j\) term of any residual with that bound has constant
leading coefficient by the order-\((j+3a-1)\) commutator coefficient,
so one can subtract one term at a time. The recursion for
\(L^{3a}=P\) preserves bounded coefficients at infinity. If all
coefficients of \(P\) were bounded, the same would hold for \(Q\).
The derivation induces zero on the bounded ring's residue field,
contradicting a nonzero scalar commutator there.

Let \(r>0\) be the maximum of the ratios of positive pole degree
at \(t=\infty\) to the corresponding order deficit of \(P\).
Assign weight one to \(t\) and \(r\) to \(E\). The same formal
root recursion bounds \((L^j)_+\) by weight \(jr\), realized by
its monic term. Hence \(P,Q\) have weights \(3ar,3br\). A contraction
loses \(r+s\). The bound \(r\ge1/(3a)\), together with \(b>a\)
and \(s\le2\), gives \((3a+3b-1)r-s>0\).
Their top symbols therefore Poisson-commute with
\(\{\xi,t\}=t^{1-s}\). The weighted Euler identities give the
symbols \(H^a,H^b\), where after a fractional Laurent extension
\[
        H=t^{3r} f(\xi/t^r),\qquad f(\eta)=\eta^3+c_2\eta+c_3 .
\]
The coefficient of \(\eta^2\) vanishes by depression, and at
least one of \(c_2,c_3\) is nonzero by the choice of \(r\).
Consequently \(f\) has at least two distinct roots and has at
least one simple root. The coefficients of \(\xi^{3a-2},\xi^{3a-3}\) in \(H^a\) show that \(jr\) is integral whenever \(c_j\ne0\), for
\(j=2,3\). In particular \(r\ge1/3\). In the nonconstant cases
invariance also gives \(jr\equiv mj\pmod3\) at each such \(j\).

Set \(R=Q^a-P^b\). This has weight less than \(3abr\), while
\([P,R]=a\kappa Q^{a-1}\) has weight \(3b(a-1)r\). As long as this
leading commutator is retained, the weight of \(R\) is at least
\(\theta=(3ab-3a-3b+1)r+s\). Above \(\theta\), its leading symbol
Poisson-commutes with \(H\), so the Euler identities give a
constant ratio of powers with \(H\). Since \(H\) has a simple
factor after adjoining a fractional Laurent variable, valuation
at that factor makes the symbol a constant multiple of an
integral power \(H^i\). The weight range gives
\((a-1)(b-1)\le i<ab\). Every integer in this range is of the form
\(au+bv\) with \(u,v\) nonnegative. To see this, choose \(0\le v<a\)
with \(bv\equiv i\pmod a\); the difference \(i-bv\) is a multiple
of \(a\) at least \(1-a\). Subtract the corresponding constant
multiple of \(P^uQ^v\) from \(R\). If \(v>0\), its commutator with
\(P\) has weight \(3(i-b)r<3b(a-1)r\). Thus the descent ends at
\(\theta\) with leading symbol \(F\) satisfying
\[
                        \{H,F\}=\kappa H^{c},\qquad c=(a-1)(b-1).
\]

## The coefficient equation

Write \(F=t^\theta G(\eta)\) for \(\eta=\xi/t^r\). The Poisson
equation becomes
\[
                    \theta f'G-3r fG'=\kappa f^c .
\]
For \(W=f^{1-c}G\), the result is
\[
                W'-\frac{r+s}{3r}\frac{f'}{f}W
                            =-\frac{\kappa}{3r}.                  \tag{1}
\]
Here differentiation is in \(\eta\). The pole and degree
comparison is as in the leading-factor residue theorem. In
detail, the multiplier of the leading pole term at a finite
point is its negative exponent minus a nonnegative multiple of
the multiplicity of \(f\). It cannot vanish. Thus \(W\) is a
polynomial. It must have a simple zero at every root of \(f\),
and a root of multiplicity \(\mu\) is impossible when
\(\mu(r+s)=3r\). Since \(f\) has at least two distinct roots,
\(M=\deg W\ge2\). Comparing the coefficient of
\(\eta^{M-1}\) in (1) gives
\[
                               M=1+s/r.
\]

If \(m=1\), then \(r\le2\). The integral congruences for \(2r\)
or \(3r\) give \(r=1\) or \(2\). At \(r=1\) a simple root
violates (1). At \(r=2\) the congruence forces \(c_2=0\), so
\(f=\eta^3+c_3\) has three distinct roots but \(M=2\).
If \(m=2\), then \(r\le1\). The congruences give \(r=1/2\)
or \(1\). At \(r=1/2\) a simple root again violates (1).
At \(r=1\) one has \(c_2=0\) and \(M=2\), again impossible.

For constant \(h\), the condition \(r\ge1/3\) and the integral
degree \(M\) give \(r=1/3,1/2,1\). At \(r=1/2\) a simple root
violates (1). At \(r=1/3\) one has \(H=\xi^3+c_3t\).
Take polynomial Weyl coordinates \(Y=E,\ Z=-(E^3+c_3t)/c_3\),
so \([Z,Y]=1\) and \(t=-Z-Y^3/c_3\). The lower terms of \(P,Q\)
after subtracting \((E^3+c_3t)^a,(E^3+c_3t)^b\) have weights
less than \(a,b\), respectively. A monomial \(t^jE^i\) of weight \(j+i/3\)
has order in \(Z\) at most \(j\). Hence the exact differential
orders in the new presentation are \(a,b\), contrary to
the [coprime scalar obstruction](ALL_COPRIME_ALGEBRAIC_STRING_OBSTRUCTION.md).

Finally \(r=1\) gives \(M=2\), so \(f\) has a double root and a
simple root. The Bernstein leading symbol of \(P\) is the
\(a\)-th power of a homogeneous cubic with a repeated and a simple linear
factor. These two factors are independent linear forms and can
be made the two Weyl coordinates by a linear change and scaling.
Then the leading monomials of \(P,Q\) are proportional to
\(x^{2a}\xi^a,x^{2b}\xi^b\), respectively. The iterated adjoint of \(P\)
on \(x\) is nonzero in every order: at step \(j\) it has leading
monomial proportional to \(x^{1+(2a-1)j}\xi^{(a-1)j}\), and the next leading
bracket has multiplier \(a(1+j)\), apart from the leading scalar
of \(P\). This excludes generation of the Weyl algebra, since
the adjoint of an automorphism image of either Weyl generator
is locally nilpotent. But a nongenerating Weyl pair cannot have
Bernstein degree gcd three by
[Guccione–Guccione–Valqui, Corollary 7.4](sources/guccione-guccione-valqui-1111.6100v3.pdf).
This rules out constant \(h\).

## Higher roots

Let \(x_0\) be a root of \(h\) of multiplicity \(m>3\).
Put \(\tau=x-x_0=t^{-e}\), with \(e=3/\gcd(3,m)\), and take
\(E\) proportional to \(\lambda D\) for \(\lambda^3=h\).
Here \([E,t]=t^{1-s}(1+O(t^{-1}))\), with \(s=e(m/3-1)>0\).
The normalized order-zero coefficients before depression
are bounded because substitution of \(D\) in terms of \(E\)
keeps a right factor \(E\) in every positive-order term.

Apply the
[local normalization and reverse-shear comparison](POLYNOMIAL_GCD_TWO_HIGHER_ROOT_EXCLUSION.md#the-quadratic-leading-symbol-at-a-higher-root)
with monic orders \(3a,3b\). Its bounded-coefficient
argument gives \(r>0\) in the depressed pair. A reverse shear
of degree greater than \(r\) would give an unbounded leading
constant coefficient equal to its \(3a\)-th power, so the
reverse shear has degree at most \(r\). The two leading
polynomials in \(\eta=\xi/t^r\) therefore have a common root at
the leading shift of \(\eta\). Weight zero for the first possible
Poisson bracket is impossible: the bracket at the common root
would vanish although the assumed scalar commutator is nonzero.
Thus it is positive and its top symbol vanishes, giving
the leading powers \(H^a,H^b\) with the same depressed cubic
\(f\) as above. Its simple factor suffices for the semigroup
descent and equation (1), which use only \(s>0\) and the leading
derivation symbol. In particular \(M=1+s/r\) is an integer
at least two.

If \(s\ge2r\), the original coefficients satisfy
\[
                  \operatorname{ord}_\tau p_i
                   \ \ge\ i+\frac{(s+r)i-3ar}{e}.
\]
Indeed the coefficient at \(D^i\) in \(E^j\) has valuation at
least \(i+j(m/3-1)\), and the weight bound gives valuation at
least \(-(3a-j)r/e\) for the coefficient at \(E^j\) in the
undepressed \(P\). The sum is smallest at \(j=i\).
Thus \(\tau^i\mid p_i\) for all \(i\ge a\). The
[original polynomial weight comparison](POLYNOMIAL_GCD_TWO_HIGHER_ROOT_EXCLUSION.md#return-to-the-original-polynomial-pair)
applies because the initial weights at the root are in ratio
\(a:b\). At the final derivative-minus-position weight, the order
of \(P\) is at most \(a-1\), so the ratio forces both weights to
be zero. Both operators then preserve \((\tau)\), which is
incompatible with their scalar commutator.

Therefore \(s<2r\), hence \(M=2\) and \(r=s\).
The depressed cubic \(f\) has exactly two distinct roots, and
in particular \(c_2\ne0\). If \(m\) is not a multiple of three,
then \(e=3\) and \(r=s=m-3\). The deck transformation acts by
\((t,E)\mapsto(\zeta t,\zeta^{-m}E)\), since \(\lambda\) is
\(t^{-m}\) times a unit in \(k((\tau))\). It preserves depression.
The coefficient of \(\xi^{3a-2}\) in \(H^a\) is nonzero and
proportional to \(t^{2r}\). Invariance forces
\(2r\equiv-2m\pmod3\), inconsistent with \(r=m-3\). Thus every
higher root has multiplicity a multiple of three, at least six.
It follows that a single root must also have multiplicity at
least six, which finishes the degree bound.

## Multiple roots and perfect cubes

**Corollary PolynomialGcdThreeMultiRootMultiplicities.** Under the
hypotheses above, if \(h\) has at least two geometric roots, it has
at least two roots of multiplicity one or two and at least one root
of multiplicity divisible by three and at least six.

If all multiplicities were divisible by three, \(h\) would be a
cube over the algebraic closure, contrary to the
[perfect-power obstruction](POLYNOMIAL_WEYL_LEADING_FACTOR_RESIDUE.md#perfect-power-leading-factors).
The only remaining multiplicities are one and two. If \(s,t\)
count the simple and double roots, respectively, then \(s+t>0\)
and \(s+2t\equiv0\pmod3\) because \(3\mid\deg h\).
Consequently \(s+t\ge2\). The leading-factor residue theorem
requires some multiplicity greater than three, which must be
a multiple of three at least six. In the one-root case, \(h\)
is a scalar multiple of \((x-x_0)^d\) with \(3\mid d\).

## Normalization at infinity

**Corollary PolynomialGcdThreeInfinityDegreeBound.** Under the hypotheses
of PolynomialGcdThreeLeadingFactor the degree \(d=\deg h\) is divisible
by three and at least fifteen. A polynomial shear of \(D\) gives
\[
                    \deg p_i\le a(d-3)+i,\qquad
                    \deg q_i\le b(d-3)+i .
\]
Here \(p_i,q_i\) are the coefficients of \(D^i\) after the shear.
The Bernstein degrees of the sheared pair are exactly
\(a(d+3),b(d+3)\), and the pair is nongenerating.

The divisibility of \(d\) follows from the theorem, including when
\(h\) has one root. Use the shear \(D\mapsto D-A(x)\), where \(A\)
is the polynomial part of \(p_{3a-1}/(3a p_{3a})\) at infinity, so
\(\deg(p_{3a-1}/p_{3a})\le-1\) at infinity in the sheared pair.
With \(\lambda^3=h\) and \(E=\lambda D\), let \(u=d/3-1>0\).
The leading derivation symbol is proportional to
\(x^{1+u}\partial_x\), which we normalize to have multiplier one.
The subprincipal coefficient of the normalized monic \(P\) prior
to this constant normalization is
\[
        p_{3a-1}/(\alpha\lambda^{3a-1})
                       -\binom{3a}{2}\lambda' ,
\]
of degree at infinity at most \(u\). The shift depressing the
monic operator consequently has degree at most \(u\).

If a lower coefficient of depressed \(P\) exceeds its order
deficit times \(u\), let \(r>u\) be the largest corresponding
ratio. Contractions lose \(r-u\), and the possible leading
commutator weight is positive. Thus the leading symbols are
again powers of a depressed cubic \(H\), with simple factor
and at least two distinct roots. In the descent to (1) now use
\(s=-u\). The threshold \((3ab-3a-3b+1)r-u\) still exceeds
\(3(ab-a-b)r\), so the same semigroup subtractions work. In (1)
the coefficient \((r-u)/(3r)\) is positive. The finite-pole argument
makes \(W\) a polynomial with at least two zeros, but the degree
identity gives \(\deg W=1-u/r<1\). This is impossible.

The weight \(u\) on \(E\) now bounds the depressed \(P,Q\) by
\(3au,3bu\), using the formal-root recursion for \(Q\) since
contractions do not increase weight. The reverse shift preserves
the bounds. In \(E^j\), the degree in \(x\) of the coefficient
at \(D^i\) is at most \(i+ju\), by induction under left
multiplication by \(E\). This gives the stated estimates and the
unique Bernstein leading monomials proportional to \(x^{ad}\xi^{3a}\)
and \(x^{bd}\xi^{3b}\), respectively. The Bernstein degrees are
\(a(d+3),b(d+3)\). Iterated leading commutators of \(P\) with
\(x\) never vanish: at step \(j\) the multiplier of the next
bracket is \(3a+a(d-3)j\) times the leading scalar of \(P\).
Thus the pair is nongenerating, and the Guccione–Guccione–Valqui
bound gives \(d+3>15\). Divisibility by three yields \(d\ge15\).

## Degree-pair exclusions

**Corollary PolynomialGcdThreeDegreePairBounds.** With \(a<b\) and
\(M=d+3\), a pair as above must satisfy \(aM\ge60\) and \(bM\ge105\).

The [Weyl generation bound through 59](POLYNOMIAL_WEYL_BERNSTEIN_FIFTY_SIX.md)
applies to the sheared pair, so \(aM\ge60\). The degrees \(aM,bM\)
are nondividing since \(a,b\ge2\) are coprime. The
[exact degree-pair transfer](POLYNOMIAL_WEYL_CONDITIONAL_THROUGH_FORTY_NINE_FROM_F7.md#two-inputs-from-polynomial-degree-pairs)
would give a complex plane Keller counterexample of these degrees.
[Nguyen, Theorem 3.6](https://doi.org/10.2989/16073606.2025.2482655)
excludes maximum degree at most \(104\).

In particular, taking \(M\) through the multiples of three starting
at \(18\), the full list with \(60\le aM<bM\le104\) consists of
\((72,90),(63,84),(72,96),(60,90),(66,99)\), in order of increasing
\(M\). The inequalities used to enumerate the list are
\(\max(2,\lceil60/M\rceil)\le a<b\le\lfloor104/M\rfloor\) and
\(\gcd(a,b)=1\). Each of the five also lies within Moh's degree-100
bound. The two additional possibilities through maximum degree
\(107\) are \((63,105)\) and \((84,105)\), both with \(M=21\).
They are excluded if one imports the planar bound through \(107\)
from the [GGHV reduction](../../BREAKTHROUGH_TARGET_ASSESSMENT.md#four-variable-hessian-conjecture-and-planar-degree-bounds).
That reduction below \(125\) supplies an existential reduced
representative of degrees \(72,108\). To use its list to exclude
another specific transferred pair of maximum degree at least \(108\)
would require checking that the reduction preserves the relevant
degrees or excludes the transferred shape.

## Exclusion in all orders

**Theorem PolynomialGcdThreeNondividingOrdersExclusion.** Over a
characteristic-zero field, polynomial differential operators of orders
\(3a,3b\), for coprime integers \(a,b\ge2\), cannot have a nonzero
scalar commutator.

Work over the algebraic closure. The one-root and multiple-root cases
both give a root \(x_0\) of \(h\) of multiplicity \(m>3\). For one
root this follows from PolynomialGcdThreeLeadingFactor, and for multiple
roots it follows already from the leading-factor residue theorem. The
[translated-root crossing theorem](POLYNOMIAL_WEYL_TRANSLATED_ROOT_CROSSINGS.md#the-finite-list-for-a-fixed-differential-order-gcd)
then gives \(u=3\). Its only row with \(u=3\) is
\(C_0=(3,6), C_1=(1,0)\) at direction \((3,-1)\), with
\((f_1,f_2)=(2,4)\). It yields \(C_2=(2-1/3,2)\), forbidden
by Guccione–Guccione–Valqui, Proposition 7.3(6).

## Consequences for reduced orders

**Corollary PolynomialWeylReducedOrderFourBound.** In a proper
characteristic-zero polynomial Weyl pair both differential orders
are composite and at least eight, and neither can be nine.
Triangular mate subtraction gives nondividing orders with gcd at
least four.

Use the [preliminary reduced-order bounds](POLYNOMIAL_GCD_TWO_AND_REDUCED_ORDER_BOUNDS.md).
If one order is six or nine, keep that member fixed and repeatedly
subtract powers of it from the other member while the latter's
order is a multiple of the fixed order. The residual order is at
least six and is no longer a multiple of the fixed order. Since
every proper divisor of six or nine is less than six, the two
orders are nondividing. Their gcd is at least three by the
coprime and gcd-two exclusions, and its only remaining possible
value is three, which the theorem excludes. The general triangular
reduction then gives the claimed gcd bound.

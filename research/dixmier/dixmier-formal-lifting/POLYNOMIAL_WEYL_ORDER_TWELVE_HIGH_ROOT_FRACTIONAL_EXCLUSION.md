# The multiplicity-26 Weyl edge and the fractional-power mate

Work over \(\mathbb C\). The case under consideration has reduced
orders \(12,18\), with leading coefficients proportional to \(h^2,h^3\).
Assume there is a multiplicity-26 root of \(h\), translated to zero,
and at least one small root (a root of multiplicity less than six).
Put \(e=\gcd(6,\{\operatorname{mult}_r h\})\), so \(e=1\) or \(2\).

The [translated high-root calculation](POLYNOMIAL_WEYL_ORDER_TWELVE_RESIDUAL.md#the-crossing-at-520-inside-the-gcd-six-support)
gives the edge from \((6,26)\) to \((5,20)\) in Fourier coordinates
of weight \((6,-1)\). In the original differential coordinates use
weight \((-1,6)\) on \((x,y)\) and put \(z=x^6 y\).
The two top symbols are proportional to \(R^2,R^3\), where
\[
                  R=x^{-10}r(z),\qquad
                  r(z)=\rho z^5(z+\beta),\quad \rho\beta\ne0 .
\]
In particular their weights are 20 and 30. The Poisson loss is 5.
These local endpoints use the multiplicity-26 branch of the
order-twelve calculation and do not require a bound on \(\deg h\).
The [next-face argument](#exclusion-by-the-next-newton-face)
excludes the entire stated branch.

## The negative fractional-power coefficients

**Proposition PolynomialWeylOrderTwelveHighRootPureCube.** In the
stated multiplicity-26 branch, after constant scalings, subtraction
of a scalar and a multiple of \(P\) from \(Q\), and translation of
\(P\) by a scalar, the pair satisfies
\[
                 Q=(P^{3/2})_+,\qquad [Q,P]=\kappa\ne0 .
\]
The leading coefficient of \(P\) can be written \(h^2\), and the
positive root in the formal differential expansion is chosen to
have leading coefficient \(h\). The changes preserve the
multiplicity-26 edge.

The [small-root fractional-power formula](POLYNOMIAL_WEYL_SMALL_ROOT_SPECTRAL_CONSTRAINTS.md#fractional-powers-in-the-differential-mate)
gives the powers \(P^{j/(2e)}\) in \(Q\).
Normalize the coefficient of \(P^{3/2}\) to one. Subtract the
\(P\) and constant terms and absorb the \(P^{1/2}\) term by
translation of \(P\): the differential part of
\((P+d)^{3/2}\) is \((P^{3/2})_++(3d/2)(P^{1/2})_+\).
The other binomial terms have negative differential order.
This proves the claim for \(e=1\). For \(e=2\), after also
rewriting the remaining terms in the translated variable,
only \(P^{j/4}\) with \(j=1,3,5\) could remain.

Suppose the largest such \(j\) has nonzero coefficient.
Then \(Q^2-P^3\) has differential order \(18+3j\) and
leading coefficient proportional to \(h^{3+j/2}\).
Its weight is at least \(30+5j>15\) by that monomial
at the high root. Its commutator with \(P\) is
\(2\kappa Q\) of weight 30, so its top symbol must Poisson
commute with \(R^2\).
A homogeneous symbol \(x^{-w}A(z)\) doing so has
\[
                 w A r'-10 A' r=0.
\]
At the simple root of \(r\), this forces \(w\) to be a
nonnegative multiple of 10 and \(A\) proportional to
\(r^{w/10}\). In particular the next possible weight above
\(30+5j\) forces differential order at least
\(18+3(j+1)\), contradicting the order \(18+3j\).
Thus the odd coefficients also vanish.

## The four-simple-root profile is excluded

**Theorem PolynomialWeylOrderTwelveHighRootOddLowMultiplicityExclusion.**
Under the above hypotheses \(e=1\) is impossible.
In particular the degree-thirty boundary profile
\(h=A\,x^{26}(x^4-\lambda_0)\), \(A\lambda_0\ne0\), is excluded.
More generally in the degree-thirty branch with a high root
of multiplicity 26, the remaining multiplicities must be
\((4)\) or \((2,2)\).

Put \(U=P^{1/2}\) after the normalization just proved and
write \(S=U^3-Q\), of negative differential order. Then
\([S,P]=-\kappa\).
If a rational pseudodifferential operator \(T\) has order
\(j>-11\) and scalar commutator with \(P\), its leading
coefficient \(t_j\) must satisfy
\[
               j t_j (h^2)'-12 h^2 t_j'=0 .
\]
Indeed the left side is the coefficient of \(D_x^{12+j-1}\),
which has positive order. A nonzero solution is proportional
to \(h^{j/6}\). For the negative orders \(-10,\ldots,-1\),
rationality when \(e=1\) permits only \(j=-6\).
Thus \(T=S-cU^{-1}\), for some constant \(c\), has order
exactly \(-11\): subtract its order-\(-6\) leading term if
present, and use the nonzero scalar commutator to exclude
order below \(-11\).

The leading coefficient \(t_{-11}\) solves
\[
         -11(h^2)'t_{-11}-12h^2t_{-11}'=-\kappa .
\]
In fact \(t_{-11}=\kappa B/(12h^2)\), where \(B\) is the
[polynomial residue primitive](POLYNOMIAL_WEYL_LEADING_FACTOR_RESIDUE.md)
satisfying \(B'-(h'/6h)B=1\). The difference is proportional
to \(h^{-11/6}\), which is not rational for \(e=1\).
The zero of \(B\) at the high root is simple.

Form the ordinary polynomial differential operator
\[
                        W=Q^2-P^3+2cP .
\]
Its order is 7 and its leading coefficient is
\(-\kappa Bh/6\). To see this, substitute
\(Q=U^3-cU^{-1}-T\). The terms of order 12 cancel and the
highest remaining terms are \(-U^3T-TU^3\).
The coefficient of \(D_x^7\) vanishes to order 27 at
zero, so the weight of \(W\) is at least \(6\cdot7-27=15\).
Its commutator with \(P\) is \(2\kappa Q\).
If its weight were higher than 15 then its top face would
Poisson commute with \(R^2\), and the homogeneous
centralizer calculation above would require differential
order at least 12. Its weight is therefore 15, with
top symbol \(x^{-15}w(z)\), \(\deg w=7\).

The face equation is \(\{x^{-15}w(z),R\}=\kappa R^2\), or
\[
                    15 w r'-10 w' r=\kappa r^2.        \tag{1}
\]
At a root of \(r\) of order \(u=5\) or \(1\), a leading
term of \(w\) of order \(\nu\) has nonzero coefficient
in (1), because \(15u-10\nu\ne0\) for odd \(u\).
Matching orders gives \(\nu=u+1\).
Thus \(w/r\) is a polynomial vanishing at both roots.
Its degree is \(7-6=1\), a contradiction.

For the degree-thirty boundary branch the high root accounts
for 26 and the remaining roots for four. If any of those
multiplicities is odd, their common gcd with 6 and 26 is one.
The only remaining partitions are \((4)\) and \((2,2)\).

## The remaining even profiles

**Corollary PolynomialWeylOrderTwelveHighRootEvenProfile.**
If the multiplicity-26 branch with a small root survives,
then \(e=2\), its moving \(P\)-phase count at every small root
is four, and the local translation stabilizer \(c_i\) has
order two. At degree thirty, this requires the following
extra primitive indicial classes, in addition to the
integer classes \(0,\ldots,k-1\) with \(k=6-u\).

- In the \((2,2)\) partition, each double root has the
  extra classes \(1/2,5/2\bmod 4\).
- In the \((4)\) partition, its quadruple root has
  \(0,0,1,1\bmod 2\), or \(1/2,1/2,3/2,3/2\bmod 2\).

The pure cube proposition still applies for \(e=2\).
For \(S=U^3-Q\) the leading rational term of negative
order above \(-11\) can now have order \(-3,-6,-9\).
Order \(-3\) is impossible: \(Q^2-P^3\) would have order
15 and weight at least 25 from its leading coefficient
proportional to \(h^{5/2}\). Its commuting top face
would require order at least 18. Subtract \(cU^{-1}\)
as above to remove order \(-6\) if present.
The remaining operator must have order \(-9\); otherwise
it has order \(-11\) and the preceding degree-one
contradiction applies. The leading coefficient at
order \(-11\) in that case is still \(\kappa B/(12h^2)\),
since \(h^{-11/6}\) is not rational for \(e=2\) either.

The nonzero coefficient at order \(-9\) in
\(S-cU^{-1}\) is a constant times \(h^{-3/2}\).
In the fixed scalar wave with \(Z^{12}=\lambda\)
it contributes a nonzero term of degree \(-9\) to the
\(Q\) logarithmic derivative equation, after cancellation
of higher terms. Equivalently the fixed spectral phase
has a term of degree \(12-9=3\).
All its higher degrees are multiples of \(6/e=3\) by
the [global phase theorem](POLYNOMIAL_WEYL_SMALL_ROOT_SPECTRAL_CONSTRAINTS.md#common-spectral-phases-and-the-opposite-fourier-coefficient).
Hence its higher phases are in an orbit of size exactly
four. The same count occurs at each small root, and
\(r_P\mid 2c_i\) forces \(c_i=2\).

At a double root, \(k=4\) and the two extra classes must
form a translation-by-two orbit. A pair of integer
classes has even sum modulo 4, but their sum must be
3 by the [sector congruence](POLYNOMIAL_WEYL_SMALL_ROOT_SPECTRAL_CONSTRAINTS.md#indicial-constraints-for-each-monodromy-eigenvalue).
The two classes are half-integral by the double-root
result there. Of the half-integral pairs, only
\(1/2,5/2\) have sum 3 modulo 4.

At a quadruple root, \(k=2\) and translation is by one.
A noninteger sector of mass two would have classes
\(\eta,\eta+1\) and sum \(2\eta+1\equiv 1\bmod2\),
forcing \(\eta\) integral. Thus either all extras are
integer or all four belong to a single noninteger
sector. In that sector (3) gives
\(4\eta+2\equiv2\bmod2\), so it is half-integral.
Translation invariance yields precisely the two lists.

As an additional check in the surviving case, \(W\) from
above has differential order 9 and weight 15. A larger
weight would force a commuting top symbol of order
at least 12. Its leading coefficient gives weight 15
exactly, and (1) has the unique polynomial solution
\[
                  w=rG,\qquad
         G=-\frac{2\kappa}{15\beta^2}\,z(z+\beta)(z-\beta/2).
\]
Indeed the root orders in (1) force the two zeros of
\(G\), and its equation
\(G'-(r'/2r)G=-\kappa/10\) determines the remaining
factor. Both the small-root classes and this high-root
face are necessary conditions for the remaining even
branches. They supply no polynomial Weyl pair.

## A primitive length-two commutator problem

**Proposition PolynomialWeylOrderTwelvePrimitiveCommutatorChain.**
In the surviving degree-thirty even branch, write \(h=H^2\),
absorbing a nonzero constant into \(H\). With the normalization above,
\(\operatorname{lc}P=H^4\), \(\operatorname{lc}Q=H^6\), and
\(\operatorname{lc}W=\delta H^3\) for \(\delta\ne0\). Thus \(P,W\)
have orders \(12,9\), reduced ratios \(4:3\), and common leading
factor \(H\). Its root multiplicities are \((13,2)\) or \((13,1,1)\),
so their gcd with the new order gcd three is one. They satisfy
\[
 \begin{aligned}
 C&=[W,P]=2\kappa Q, & [C,P]&=2\kappa^2,\\
 C^2&=4\kappa^2(W+P^3-2cP),
 & [W,C]&=2\kappa^2(3P^2-2c).
 \end{aligned}                                                   \tag{2}
\]
In particular \(P,W\) generate the same associative algebra as
\(P,Q\). Conversely, polynomial differential operators \(P,W\)
satisfying the second-commutator and squared-commutator equations
in (2), for constants \(\kappa\ne0,c\), recover a polynomial Weyl
pair by \(Q=[W,P]/(2\kappa)\), with exactly the defining identity
\(W=Q^2-P^3+2cP\).

There is also a controlled first binomial resultant. For the rational
pseudodifferential root \(V=P^{1/4}\) with leading coefficient \(H\),
\[
 \begin{aligned}
 \operatorname{ord}(W-\delta V^3)&=7,&
 \operatorname{lc}(W-\delta V^3)&=-\frac{\kappa}{6}BH^2,\\
 B'-\frac{H'}{3H}B&=1,&
 \operatorname{ord}(W^4-\delta^4P^3)&=34,\\
 \operatorname{lc}(W^4-\delta^4P^3)&=
       -\frac{2\kappa\delta^3}{3}BH^{11}.
 \end{aligned}                                                   \tag{3}
\]
Here \(B\) is the polynomial residue primitive already supplied
by the original scalar Weyl pair; no scalar commutator for \(P,W\)
is assumed.

The order-nine leading coefficient was obtained above from the
nonzero order-minus-nine coefficient proportional to \(h^{-3/2}\).
This proves the leading-factor statements. The first line of (2)
is the product rule applied to the defining identity for \(W\).
The squared equation is that identity itself. Since
\([P,Q]=-\kappa\), the same product rule gives
\([W,Q]=\kappa(3P^2-2c)\), proving the remaining equation.
The converse and equality of generated algebras follow by the
displayed formula for \(Q\). The converse alone makes no assertion
about differential orders or root multiplicities; these remain
part of the specified branch.

To prove (3), set \(T=W-\delta V^3\), so \(\operatorname{ord}T<9\)
and \([T,P]=2\kappa Q\). If \(T\) had order eight, its leading
coefficient \(t\) would solve \(8t(H^4)'-12H^4t'=0\), and hence
be proportional to \(H^{8/3}\). This is not rational, already at
the multiplicity-thirteen root. Order below seven cannot produce
the order-eighteen commutator. Thus its order is seven and
\[
                 7t(H^4)'-12H^4t'=2\kappa H^6.
\]
Its rational solution is \(t=-\kappa BH^2/6\). The difference of
two solutions is proportional to \(H^{7/3}\), again not rational.
Expand \(W^4=(\delta V^3+T)^4\). The four words with one \(T\)
have order \(3\cdot9+7=34\) and the same nonzero leading
coefficient; words with two or more \(T\)'s have order at most
32. The term without \(T\) is \(\delta^4P^3\). This proves (3).

### Why the scalar critical-face argument stops here

The closest critical-edge input used in the
[primitive-factor proof](POLYNOMIAL_WEYL_PRIMITIVE_LEADING_FACTOR_OBSTRUCTION.md#no-high-root-in-the-primitive-case)
is Guccione–Guccione–Valqui, *The Dixmier conjecture and the shape
of possible counterexamples*, arXiv:1111.6100v3,
[Lemma 7.1, printed pages 37–38](sources/guccione-guccione-valqui-1111.6100v3.pdf).
Its statement and proof require a scalar commutator. In particular
its zero-leading-bracket conclusion does not extend to \(P,W\).

Indeed write the known top forms at weight \((-1,6)\) as
\(p_0=R^2\), \(q_0=R^3\), \(w_0=x^{-15}w(z)\). Their weights
are 20, 30 and 15, and their already established equation gives
\[
                    \{w_0,p_0\}=2\kappa R^3\ne0,
       \qquad \{\{w_0,p_0\},p_0\}=0.                         \tag{4}
\]
Thus the first commutator is supported on this face, while the
leading part of the second commutator cancels. Moreover \(p_0\)
is not a fourth power: at \(z=-\beta\) its multiplicity is two.
This is exactly where the scalar proof's common-power factorization
would fail for the new ratios \(4:3\).

The first resultant in (3) also survives on the same face. For
the display normalize \(\rho=1\) and put
\(A=-2\kappa/(15\beta^2)\), so \(\delta=A\) there. Its top form is
\[
 x^{-60}A^4z^{24}(z+\beta)^6
       \big((z+\beta)^2(z-\beta/2)^4-z^6\big).                  \tag{5}
\]
The parenthesis has degree four and leading coefficient
\(-3\beta^2/2\). Hence (5) has differential degree 34 and is
nonzero, consistently with (3). The new primitiveness by itself
therefore yields no contradiction through this first resultant.

The [next-face argument below](#exclusion-by-the-next-newton-face)
excludes the even profiles using the polynomial order-nine
operator and its first commutator. The equations in (2)
retain the original pair on the face used there.

The existing [exact checker](verify_weyl_primitive_leading_resultant.py)
also checks (2) in the abstract Weyl algebra over \(\mathbb Q[\kappa,c]\),
the cancellations and degree in (5), and (4) for three nonzero
rational \(\beta\)'s. Run
`python3 research/dixmier/dixmier-formal-lifting/verify_weyl_primitive_leading_resultant.py`
from the repository root. These are exact internal controls for the length-two
reduction and the order-34 resultant.


## Exclusion by the next Newton face

**Theorem PolynomialWeylOrderTwelveMultiplicityTwentySixExclusion.**
A polynomial Weyl pair of reduced differential orders \((12,18)\)
with nonzero scalar commutator cannot have a common leading factor
with a root of multiplicity 26 and a root of multiplicity less than
six. In particular, at common-factor degree thirty, the only
remaining root profile in the \((12,18)\) row is a single root
of multiplicity thirty.

By the odd-multiplicity exclusion we may assume \(e=2\).
Use the normalized pair and its polynomial \(W\) above.
A diagonal Weyl automorphism \(x\mapsto a x,\ y\mapsto a^{-1}y\)
changes \(z\) by a fifth power. Over \(\mathbb C\) it allows
\(\beta=1\); constant rescaling of \(P,Q\) gives \(r=z^5(z+1)\).
Rescale \(W\) to an operator \(B\) with \([B,P]=Q\), so that
its highest symbol at weight \((-1,6)\) is
\(-x^{-15}r z(z+1)(2z-1)/30\).
In particular, at the low-\(z\) endpoints of the known edge,
\(P,B\) contain \(x^{-20}z^{10},x^{-15}z^6/30\), respectively.

Write the symbols using symmetric ordering, so the symbol of a
commutator at weight loss \(w\) is the Poisson bracket in the
convention \(J(P,B)=P_x B_y-P_y B_x\) for \([B,P]\).
For \(w>0\) all further commutator terms have strictly lower
weight. A monomial of \(P\) has the form \(x^{i-20}z^j\)
with \(i,j\ge0\), and the corresponding form for \(B\)
is \(x^{i-15}z^j\) with \(i,j\ge0\).
At weights \((-q,6q-p)\), the losses from the weights of the
two endpoint monomials are \(qi+p(j-10)\) and \(qi+p(j-6)\).
For \(i=0\) those losses are nonnegative, with equality
only at the endpoints.

Suppose a monomial with \(i>0\) gives a ratio
\((10-j)/i\) for \(P\), or \((6-j)/i\) for \(B\), exceeding
\(1/2\). Take the largest such ratio \(q/p\) in lowest terms.
The weights have positive sum \(5q-p\). At these weights the
top symbols of \(P,B\) are
\[
        x^{-20}z^{10}U(t),\qquad x^{-15}z^6 V(t),
        \qquad t=x^p/z^q,
\]
where \(U,V\) are polynomials, \(U(0)=1,\ V(0)=1/30\), and
at least one is nonconstant. The top bracket is
\[
 x^{-30}z^{15}S(t),\qquad
 S=30UV+(20q-10p)tUV'+(6p-15q)tU'V.
\]
This follows by differentiating \(t=x^p/z^q\) in the
coordinates \(x,z=x^6y\), for which the bracket of
\(x^l A(z)\) and \(x^m C(z)\) is
\(x^{l+m+5}(lAC'-mA'C)\).
Since \(S(0)=1\), this is the top symbol of \(Q\).
Its weight is \(30q-15p\), positive and \(3/2\) times
the weight \(20q-10p\) of \(P\). The weight bound on
\([Q,P]\) is \(45q-24p\). A negative bound contradicts the
scalar commutator directly. If it is zero, then \((p,q)=(15,8)\).
Polynomial support gives \(\deg U\le1\) and \(\deg V=0\),
since the two symbols have \(y\)-exponents \(10-qj\) and
\(6-qj\) at index \(j\). It also gives \(\deg S\le1\) from
the bracket formula. Therefore the top symbols of \(P,Q\)
are divisible by \(y^2,y^7\), respectively, and their
bracket cannot be a nonzero scalar at weight zero.
Consequently the bound is positive, and their top
symbols Poisson commute.

The quotient of the square of the \(Q\)-symbol by
the cube of the \(P\)-symbol is \(S^2/U^3\).
Poisson commutation makes its derivative in \(t\) zero,
since \(t\) has weight zero and \(P\) has positive weight
in the toric chart. Its constant value is one by
\(U(0)=S(0)=1\). Thus \(U=R^2,\ S=R^3\) for a polynomial
\(R\) with \(R(0)=1\).
The bracket equation divided by \(R\) becomes
\[
 R^2=30RV+(20q-10p)tRV'-(30q-12p)tR'V.                 \tag{6}
\]
Every root of \(R\) is nonzero. At a root of multiplicity \(e\),
if \(V\) has order \(j\le e\), the derivative
terms in (6) have a nonzero coefficient at order \(e+j-1\),
since \((30q-12p)/(20q-10p)>1\). This order is strictly
less than \(2e\), a contradiction. Hence \(V=RL\) for a
polynomial \(L\) vanishing at each root of \(R\).
If \(R\) is constant, (6) reduces to
\(1=30V+(20q-10p)tV'\), so \(V=1/30\), also impossible
at the first competing face.

For nonconstant \(R\) set \(n=\deg R\) and \(k=\deg L\ge1\).
Divide (6) by \(R^2\) and compare the leading power of \(t\)
at infinity. This gives
\[
              30+(20q-10p)k=(10q-2p)n .
\]
Polynomial support gives \(2qn\le10\) from the \(y\)-exponent
of the last term of the \(P\)-symbol. It follows that
\[
        30+10(2q-p)\le (10q-2p)n\le 50-10p/q,
\]
which forces \(q\le1\) since \(2q-p>0\).
The case \(q=1\) requires \(p=1\) and equality in the
bound, with \(n=5\). But then the last \(P\)-term would
be \(x^{-20} z^{10}(x/z)^{10}=x^{-10}\), which is not a
polynomial. This rules out the first competing face.

It follows that every monomial of \(P\) has weight at
most zero at \((-1,4)\), and every monomial of \(B\)
has weight at most three there. Since \([B,P]=Q\)
and the commutator loses at least three units of weight,
\(Q\) has weight at most zero. The same bound on the
commutator then makes \([Q,P]\) have weight at most
minus three, contradicting the nonzero scalar.
Every multiplicity-26 profile with a small root is
therefore excluded. At degree thirty, a
multiplicity-26 profile has small roots with total
multiplicity four. The characteristic-zero statement follows
by embedding a finitely generated coefficient field, enlarged to
include the specified roots, into \(\mathbb C\).

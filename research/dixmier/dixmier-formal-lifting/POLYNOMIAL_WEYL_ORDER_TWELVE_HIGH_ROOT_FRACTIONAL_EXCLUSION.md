# The multiplicity-26 and multiplicity-30 Weyl edges

Work over an algebraically closed characteristic-zero field.
Consider reduced differential orders \(12,18\) with leading
coefficients proportional to \(h^2,h^3\). Translate a root
of \(h\) of multiplicity \(u\in\{26,30\}\) to zero.
Use the weight \((-1,6)\) on \((x,y)\) and put \(z=x^6y\).
The top symbols are proportional to \(R^2,R^3\), with
\(R=x^{-\omega}r(z)\), \(\omega=36-u\), and
\[
\begin{array}{c|c|c}
 u & \omega & r(z)\\ \hline
 26 & 10 & \rho z^5(z+\beta),\quad \rho\beta\ne0\\
 30 & 6 & zL(z),\quad \deg L=5,\quad L(0)\ne0 .
\end{array}
\]
These are the translated
[(5,20) and (6,30) crossings](POLYNOMIAL_WEYL_ORDER_TWELVE_RESIDUAL.md#the-crossing-at-520-inside-the-gcd-six-support)
in the Fourier coordinates. In both cases \(r\) has degree
six and at least one simple root. The commutator loses weight
five. Write \(e=\gcd(6,\{\operatorname{mult}_r h\})\).

## Normalization to the pure cube

**Proposition PolynomialWeylOrderTwelveHighRootPureCube.**
At either of the stated translated edges, after constant
scalings, subtraction of a scalar and a multiple of \(P\)
from \(Q\), and translation of \(P\) by a scalar, the pair
satisfies
\[
                 Q=(P^{3/2})_+,\qquad [Q,P]=\kappa\ne0 .
\]
The normalization preserves its two top symbols, up to the
stated constant scalings.

Use the
[fractional-power expansion by differential order](POLYNOMIAL_WEYL_LEADING_FACTOR_RESIDUE.md#fractional-powers-by-differential-order)
to write \(Q=\sum_{j=0}^{3e}\gamma_j(P^{j/(2e)})_+\).
Normalize \(\gamma_{3e}=1\).
Subtract the \(P\) and constant terms and absorb \(P^{1/2}\)
by translation of \(P\). For example
\(((P+d)^{3/2})_+=(P^{3/2})_++(3d/2)(P^{1/2})_+\).
When rewriting all the other terms in the translated variable,
their indices only decrease by multiples of \(2e\);
thus they do not generate new indices divisible by \(e\).

Suppose the largest remaining \(j\) with \(\gamma_j\ne0\)
exists. Then \(e\nmid j\), and \(Q^2-P^3\) has differential
order \(18+6j/e\), with leading coefficient proportional
to \(h^{3+j/e}\). In particular its weight is at least
\(\omega(3+j/e)\) by the translated root. Its commutator
with \(P\) is \(2\kappa Q\) of weight \(3\omega\).
Since \(\omega(3+j/e)+2\omega-5>3\omega\),
its top symbol must Poisson commute with \(R^2\).

A homogeneous commuting symbol \(x^{-w}A(z)\) satisfies
\(w A r'-\omega A'r=0\). At a simple root of \(r\)
this forces \(w/\omega\) to be a nonnegative integer.
Then \(A\) is proportional to \(r^{w/\omega}\).
The next possible weight at least \(\omega(3+j/e)\)
therefore forces differential order strictly above
\(18+6j/e\), a contradiction. This proves the proposition.

## The four-simple-root profile is excluded

**Theorem PolynomialWeylOrderTwelveHighRootOddLowMultiplicityExclusion.**
At a multiplicity-26 root, \(e=1\) is impossible.
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

## The even face at multiplicity 26

**Proposition PolynomialWeylOrderTwelveHighRootEvenFace.**
At a multiplicity-26 root with \(e=2\), set \(U=P^{1/2}\)
in the pure cube normalization. There is a constant \(c\)
for which \(U^3-Q-cU^{-1}\) has differential order \(-9\).
The operator \(W=Q^2-P^3+2cP\) then has order nine and
weight fifteen at \((-1,6)\). In the notation of (1)
its top symbol is \(x^{-15}w(z)\), where
\[
        w=rG,\qquad
        G=-\frac{2\kappa}{15\beta^2}\,z(z+\beta)(z-\beta/2).
\]

Put \(S=U^3-Q\). In \(S\) or its successive remainders,
the only possible leading orders in \((-11,0)\) are
\(-3,-6,-9\): their leading coefficients are proportional
to \(h^{j/6}\) by differential-order comparison.
Order \(-3\) is impossible since \(Q^2-P^3\) would have
order fifteen and weight at least twenty-five, by the
multiplicity-26 root. A weight above fifteen requires a
commuting top face, which in this case would have weight
at least thirty and differential order at least eighteen.
Remove order \(-6\), if present, by subtracting \(cU^{-1}\).
Order below \(-11\) would make the commutator zero.
If the remainder has order \(-11\), its leading coefficient
is \(\kappa B/(12h^2)\), since \(h^{-11/6}\) is not
rational for \(e=2\). The order-seven contradiction from
(1) applies without change. Thus the remainder has
order \(-9\) and leading coefficient proportional to
\(h^{-3/2}\).

This makes \(W\) have order nine, and its leading
coefficient gives weight at least fifteen. A larger weight
would require a commuting top face of order at least twelve.
Thus the weight is fifteen and (1) holds. Its root orders
force \(w/r\) to vanish at both \(0,-\beta\).
The equation \(G'-(r'/2r)G=-\kappa/10\) gives the
displayed cubic. Substituting it verifies the equation,
and the difference of two solutions would be proportional
to \(r^{1/2}\), which is not rational.

## The small-root consequence of the even face

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

The even-face proposition supplies a nonzero order-minus-nine
coefficient proportional to \(h^{-3/2}\).
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
with a root of multiplicity 26.

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
This excludes every multiplicity-26 profile. No small-root
spectral information is used in the next-face argument.


## The multiplicity-thirty root forces a sixth power

**Theorem PolynomialWeylOrderTwelveMultiplicityThirtyPower.**
In reduced orders \((12,18)\), the presence of a root of
multiplicity thirty forces \(h\) to be the thirtieth power
of a linear polynomial up to a nonzero scalar.

Use the pure cube normalization at \(u=30\), so
\(\omega=6\). Put \(U=P^{1/2}\), \(S=U^3-Q\).
As \([S,P]=-\kappa\), the order of \(S\) is between
\(-11\) and \(-1\).
For a leading order \(j>-11\), its coefficient is
proportional to \(h^{j/6}\).
If \(j>-6\), \(Q^2-P^3\) has differential order \(18+j\)
and weight at least \(18+j>11\) at the multiplicity-thirty
root. Its commutator with \(P\) is \(2\kappa Q\), so its top
must commute with \(R^2\). This forces its weight to be a
multiple of six and its differential order to be at least
that weight. Since \(18+j\) is strictly between twelve
and eighteen, this is impossible.

Remove a leading order \(-6\) if present by subtracting
\(cU^{-1}\) from \(S\), and denote the resulting order by
\(j\le-7\). Then \(W=Q^2-P^3+2cP\) has order \(18+j\le11\).
Its weight is exactly eleven: anything lower cannot
give \(2\kappa Q=[W,P]\), and anything higher would give
a commuting top of differential order at least twelve.
Write its top as \(x^{-11}w(z)\). Its bracket equation
is
\[
                       6w'r-11wr'=\kappa r^2.              \tag{7}
\]
Here \(r=zL(z)\) has at least two distinct roots, all of
multiplicity at most five. At a root of multiplicity \(m\),
if \(w\) has order \(v\), the leading coefficient on the
left has nonzero multiplier \(6v-11m\).
Matching orders gives \(v=m+1\).
Thus \(\deg w\ge\deg r+2=8\). If \(8\le\deg w<11\),
the left side has degree \(\deg w+5>12\), with nonzero
leading coefficient. It cannot equal the right side.
Therefore \(\deg w=11\), forcing \(j=-7\).

At this leading differential order the rational coefficient
of the remainder is proportional to \(h^{-7/6}\).
Since seven is coprime to six, all multiplicities of \(h\)
are divisible by six. The
[perfect-power residue obstruction](POLYNOMIAL_WEYL_LEADING_FACTOR_RESIDUE.md#perfect-power-leading-factors)
now gives a single root, of multiplicity thirty as claimed.

## The coefficient field at the multiplicity-thirty edge

**Corollary PolynomialWeylOrderTwelveEighteenRationalExclusion.**
There is no proper polynomial Weyl pair over \(\mathbb Q\) with
nonzero scalar commutator and reduced differential orders
\((12,18)\). More generally the same holds over any
characteristic-zero field over which both cubics
\[
    64c^3-216c^2+324c-81,\qquad
    448k^3-1323k^2+3888k-3888
\]
have no root.

Take the order-twelve reduction and the multiplicity-26
exclusion. A multiplicity-thirty root remains, and the
preceding theorem implies that it is the only root of \(h\).
It therefore belongs to the coefficient field, since \(h\)
has characteristic zero. In the translated edge the common
factor \(R\) can be defined over this field by taking the ratio
of the cubed and squared highest symbols. Thus \(L\) and
the scalar equation
\[
                     6zG'+(1-5zL'/L)G=1
\]
can be taken over it, by normalization of the nonzero
constant on the right.

For ratio \(2:3\) the
[(6,30) successor table](POLYNOMIAL_WEYL_ORDER_TWELVE_RESIDUAL.md#the-local-successors-at-630)
requires a maximal root of \(L\) of multiplicity three.
The degree-five polynomial \(L\) has a unique such root,
which is nonzero and defined over the coefficient field.
Scale \(z\) to place it at one. The other multiplicities
are either \((2)\) or \((1,1)\).

In the \((2)\) case, write \(L\) up to scalar as
\((z-1)^3(z-c)^2\). Here \(c\) belongs to the coefficient
field and \(c\ne0,1\). The equation forces
\(G=(z-1)(z-c)(A z^2+B z+C)\).
Comparison from constant term upwards gives
\[
 cC=1,\qquad 7cB+(8c+3)C=0,\qquad
 13cA+(2c-3)B-12C=0,\qquad (4c+9)A+6B=0 .
\]
The value \(4c+9=0\) would force \(B=0\), contrary to the
second equation. The other three equations substituted into the
third give the first cubic in the statement.

In the \((1,1)\) case put
\(L=(z-1)^3(z^2-h_1 z+k)\) up to scalar, with
\(h_1,k\) in the coefficient field and \(k\ne0\).
Now \(G=(z-1)(z^2-h_1 z+k)(A z+B)\).
Coefficient comparison gives \(Bk=-1\) and
\(7kA=(2h_1-8k)B\).
After these substitutions the two remaining equations are
\[
 \begin{aligned}
 h_1^2-4h_1k-9h_1+15k&=0,\\
 16h_1^2-19h_1k+16k^2-21k&=0 .
 \end{aligned}
\]
Their resultant in \(h_1\) is \(7k\) times the second
cubic in the statement. Since \(k\ne0\), this gives a
root in the coefficient field.

Both cubics are irreducible over \(\mathbb Q\): they
have degree three modulo eleven and take nonzero values
at every residue. The coefficient identities, resultant
and residue check are reproduced by the
[order-twelve checker](verify_weyl_order_twelve.py).

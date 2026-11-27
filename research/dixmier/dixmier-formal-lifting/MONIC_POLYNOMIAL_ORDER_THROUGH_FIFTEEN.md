# Constant-leading polynomial Weyl generators through order fifteen

**Theorem MonicPolynomialOrderThroughFifteenGeneration.** Let \(k\)
have characteristic zero and let \(P,Q\in A_1(k)=k[x][D;\partial_x]\)
satisfy \([P,Q]=\kappa\in k^\times\). If \(P\) has nonzero constant
leading coefficient and differential order at most fifteen, then
\(P,Q\) generate \(A_1(k)\).

## Reduction to nondividing orders

Suppose the generated algebra is proper. Use
[PolynomialWeylReducedOrderBounds](POLYNOMIAL_GCD_TWO_AND_REDUCED_ORDER_BOUNDS.md)
to subtract powers of \(P\) from \(Q\) until the order of \(Q\)
is no longer divisible by the order of \(P\), keeping \(P\) fixed.
No member can have order below six or have prime order. In
particular this leaves orders \(6,8,9,10,12,14,15\) for \(P\).
Set \(g\) equal to the gcd of the two remaining orders. The
[coprime scalar obstruction](ALL_COPRIME_ALGEBRAIC_STRING_OBSTRUCTION.md)
and the [gcd-two polynomial exclusion](POLYNOMIAL_GCD_TWO_AND_REDUCED_ORDER_BOUNDS.md)
rule out \(g=1,2\). The leading commutator identity makes the
leading coefficient of \(Q\) constant because that of \(P\) is
constant. For \(g=3\), the
[gcd-three theorem](POLYNOMIAL_GCD_THREE_LEADING_FACTOR.md)
excludes constant common leading factors. Thus the remaining
values of \(g\) are \(4,5,6,7\).

Write the orders as \(ga,gb\), with \(a,b\) coprime. Here \(b\ge2\):
if \(b=1\), then the order of \(Q\) is \(g\), and the preceding
exclusions of orders below six and prime orders leave only \(g=6\).
A pair with a monic member of order six generates, since partner
descent leaves only the gcd-three case with constant leading factor.
Also \(a\ge2\) because the order of \(P\) does not divide that of \(Q\).
We now work over an algebraic closure of \(k\). A proper generated
\(k\)-subalgebra remains proper after this extension, by faithful
flatness applied to its underlying vector space and its quotient
in \(A_1(k)\).

## The common homogeneous symbol

Scale both operators to be monic and depress \(P\) by a polynomial
shift of \(D\) removing the coefficient of \(D^{ga-1}\). Retain
\(\kappa\) for the resulting nonzero scalar commutator. As in the
[gcd-two normalization](POLYNOMIAL_GCD_TWO_LOW_DEGREE_LEADING_FACTOR.md#constant-leading-laurent-normalizations),
the formal \(ga\)-th root \(L=P^{1/(ga)}\) gives \(Q\) as a finite
constant combination of the positive parts \((L^j)_+\).
Their commutators with \(P\) have order at most \(ga-2\), so each
successive leading coefficient in this descent is constant.
If all coefficients of \(P\) were constant, the root recursion
would make those of \(Q\) constant as well, contradicting
the commutator. Give \(x,D\) weights \(1,r\), where \(r>0\) is the
maximum ratio of degree in \(x\) to differential-order deficit
among lower terms of \(P\). Then \(r\ge1/(ga)\).
Root recursion bounds every term of \(L\) by weight \(r\),
so \(P,Q\) have weights \(gar,gbr\).

The first possible Poisson bracket has weight
\((g(a+b)-1)r-1>0\) and must vanish. The weighted Euler identities
give \(q^a=p^b\) for the monic leading symbols, hence
\[
              p=H^a,\qquad q=H^b,\qquad
              H=x^{gr}f(\eta),\qquad \eta=\xi/x^r .
\]
Here \(f\) is monic of degree \(g\), has zero coefficient of
\(\eta^{g-1}\), and is different from \(\eta^g\) by the choice
of \(r\). In particular it has at least two distinct roots.
If \(H\) is a proper power, write \(H=K^e\) with \(e>1\)
maximal, so the multiplicities of the distinct roots of the
monic polynomial defining \(K\) have gcd one. This form
is also available with \(e=1\) and \(K=H\).
Put \(d=g/e\), \(K=x^{dr}f_0(\eta)\), and \(w=dr\).
Then \(f_0\) is monic and depressed, and \(K\) belongs to the
symbol polynomial ring. For the latter assertion, factor
\(H\) in \(k[x,\xi]\). Monicity and evaluation over \(k(x)\)
ensure that its primitive monic root has polynomial coefficients.
Equivalently the coefficients of this root can be recovered
successively at descending \(\xi\)-degrees, dividing only by \(e\).
Here \(d\ge2\), since \(d=1\) with a depressed root would make
\(H=\xi^g\). Thus the only possibilities with \(e>1\) are
\((g,e,d)=(4,2,2),(6,2,3),(6,3,2)\). In particular \(e>1\)
is prime in every such case.

## Descent for the commutator with the primitive symbol

Set \(C=Q^a-P^b\). Its weight is less than \(eab w\), and
\([P,C]=a\kappa Q^{a-1}\) has weight \(eb(a-1)w\).
While preserving this leading commutator, the weight of \(C\)
is at least
\[
                         e(ab-a-b)w+r+1 .
\]
At a strictly larger weight, its leading symbol Poisson-commutes
with \(K\) and is a constant multiple of \(K^i\).
Indeed the Euler identities give a constant ratio of powers
with \(K\). At each distinct root of \(f_0\) the corresponding
valuation forces the exponent ratio times the multiplicity
to be integral. The multiplicities have gcd one, so the power
\(i\) itself is integral. This argument may use a fractional
Laurent extension in \(x\), without changing the finite
operators being subtracted.

At this stage \(e(ab-a-b)+1\le i<eab\). If \(e\mid i\),
then \(i/e\ge(a-1)(b-1)\), so \(i=ea u+eb v\) for
nonnegative \(u,v\). We use here the elementary conductor
bound for two coprime generators: choose \(0\le v<a\)
in the required congruence class modulo \(a\) and the
difference \(i/e-bv\) is a multiple of \(a\) greater
than \(-a\). Subtract a constant multiple of \(P^uQ^v\).
Its commutator with \(P\) has weight at most
\((i-eb)w<eb(a-1)w\), with zero commutator when \(v=0\).
If the descent reaches the lower bound, the leading bracket
there, after division by \(eaK^{ea-1}\), gives
\(\{K,F\}=cK^N\) for its final symbol \(F\), with
\(c\ne0\) and \(N\ge0\).

Otherwise the first exponent \(j\) with \(e\nmid j\) is
coprime to \(e\). Keep this \(C\) and use
\[
                      S=C^e-\gamma P^uQ^v,\qquad
                      ea u+eb v=ej ,
\]
where \(u,v\) exist because \(j\ge(a-1)(b-1)\), and
\(\gamma\) cancels the leading symbol. The weight of \(S\)
is less than \(ejw\), and its commutator with \(P\)
has weight \(((e-1)j+eb(a-1))w\). In fact the leading
symbol of \([P,C^e]\) is \(e\) times the product of the
leading symbols of \(C^{e-1}\) and \([P,C]\); the subtraction
has a lower commutator weight since \(j<eab\).
The resulting lower bound on the weight of \(S\) is
\[
               ((e-1)j+e(ab-a-b))w+r+1 .
\]
At larger weight its leading symbol is again a constant
multiple of \(K^i\). Choose \(0\le \ell<e\) with
\(e\mid(i-\ell j)\). The strict weight bound gives
\(i-\ell j\ge e(ab-a-b)+1\), so
\(i-\ell j=ea u+eb v\) with \(u,v\) nonnegative.
Subtract the corresponding multiple of \(C^\ell P^uQ^v\).
The commutator contribution from \(C^\ell\), if nonzero,
has weight \((i-j+eb(a-1))w\), below the target because
\(i<ej\). The contribution from \(Q^v\), if nonzero,
has weight \((i-eb)w\), also below the target since
\(i<ej\) and \(j<eab\). In both descents the weights
are discrete, and the nonzero leading commutator prevents
crossing the lower bound. The second descent therefore
also gives
\[
                              \{K,F\}=cK^N,\qquad c\ne0,
\]
for a final leading symbol \(F\), with \(N\ge0\).

## The rational equation on the homogeneous face

With \(F=x^\theta G(\eta)\), the last equation gives
\[
        \theta f_0'G-dr f_0G'=c f_0^N,\qquad
        \theta=dNr-(d-1)r+1 .
\]
The weight identity follows from the loss \(r+1\) in
the Poisson bracket. Consequently \(W=f_0^{1-N}G\)
satisfies
\[
             W'-\frac{r+1}{dr}\frac{f_0'}{f_0}W
                             =-\frac{c}{dr}.             \tag{1}
\]
A finite pole of \(W\) cannot occur, since its negative
leading exponent minus a nonnegative multiple of the root
multiplicity is nonzero. Hence \(W\) is a polynomial.
It vanishes at each distinct root of \(f_0\) and
its zero there must be simple to give the nonzero constant
in (1). In particular a root of multiplicity \(\mu\)
is impossible if \(\mu(r+1)=dr\). Since at least two roots
are distinct, comparison at infinity gives
\(\deg W=1+1/r\). Write \(r=1/m\), with \(m\) a positive
integer. Each nonzero lower coefficient of \(K\) at a
deficit \(j\le d\) has degree \(j/m\), and at least one
occurs. Thus \(1\le m\le d\) and
\[
                         f_0(z)=z^p F_0(z^m)
\]
for some \(p\ge0\), chosen so \(F_0(0)\ne0\).

When \(m=1\), the degree of \(W\) is two. The two
distinct roots give \(K\) two distinct homogeneous linear
factors, which can be made position and momentum by a
linear Weyl change and scaling. Thus the Bernstein
degree gcd of \(P,Q\) is \(g\le7\). This contradicts
[Guccione–Guccione–Valqui, Corollary 7.4, printed page 42](sources/guccione-guccione-valqui-1111.6100v3.pdf),
since our pair generates a proper subalgebra.

When \(m=d\), the common symbol \(K\) has the form
\(\xi^d+\sigma x\) with \(\sigma\ne0\). Introduce polynomial
Weyl coordinates \(Y=D,\ Z=-(D^d+\sigma x)/\sigma\),
with \([Z,Y]=1\). The remainders after subtracting the
corresponding powers of \(D^d+\sigma x\) from \(P,Q\)
have weights less than \(ea,eb\). A normal monomial
\(x^jD^i\) of weight \(j+i/d<ea\) has order in \(Z\)
at most \(j<ea\). Hence the exact new orders are \(ea,eb\),
both with constant leading coefficients. For \(e=1,2,3\)
these are excluded respectively by the coprime scalar theorem,
the gcd-two exclusion, or the constant-leading gcd-three theorem.

It remains to take \(2\le m<d\). If \(e>1\), then \(d\le3\);
the only possibility is \(m=2,d=3\), when \(f_0=z^3+\sigma z\)
has a simple root and violates (1). If \(e=1\), two distinct
nonzero root orbits of \(f_0\) under the \(m\)-th roots of unity
would give at least \(2m>m+1=\deg W\) distinct roots. Therefore
\[
                      f_0(z)=z^p(z^m-\sigma)^q ,
                      \qquad \sigma\ne0,\quad q\ge1 .
\]
Here \(p\ge1\), since \(p=0\) and the primitive multiplicities
would require \(q=1\), contrary to \(m<d\). In the polynomial
Weyl coordinates \(Y=D,\ Z=D^m/\sigma-x\), the weight bounds
with \(\operatorname{wt}(Y)=1,\operatorname{wt}(Z)=m\)
give \(P\) weight \(ga\). Thus its order in \(Z\) is at
most \(\lfloor ga/m\rfloor\), and when this order is \(i\)
its leading coefficient has degree in \(Y\) at most \(ga-mi\).

For \(ga=8\) or \(10\) the new order is at most five,
contrary to the reduced order bound. For \(ga=12,14,15\),
the new order is at most seven, so the only possible
remaining nonprime order is six. Its leading coefficient has
degree at most \(ga-12\le3\), whereas keeping an order-six
member of a proper pair and descending its partner gives
orders \(6,3b'\) with odd \(b'\ge3\). The gcd-three theorem
then requires that member's leading coefficient to be
proportional to a square of a polynomial of degree at least
six. This contradiction finishes the proof.

## Primary comparison

Guccione–Guccione–Valqui, [arXiv:1111.6100v3](https://arxiv.org/abs/1111.6100v3),
define \(B\) in their introduction as the minimum of the
total-degree gcd over nongenerating Weyl pairs, and prove
\(B>15\) in Corollary 7.4 over arbitrary characteristic-zero
fields. Their invariant is used only when \(r=1\).
The positive-weight Poisson comparison also appears in
Previato–Rueda–Zurro, [Lemma 5.3](https://doi.org/10.3842/SIGMA.2019.101).
Their fourth-order examples concern commuting operators.

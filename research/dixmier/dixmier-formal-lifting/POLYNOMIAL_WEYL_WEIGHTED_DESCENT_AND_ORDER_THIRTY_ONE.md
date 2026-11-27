# Weighted Weyl descent and the constant-leading order bound

**Theorem PolynomialWeylNonconstantLeadingDegreeBound.** Let \(P,Q\)
be polynomial differential operators over a characteristic-zero field
with nonzero scalar commutator, generating a proper Weyl subalgebra.
Suppose the orders \(ga,gb\) are nondividing, with \(a,b\ge2\) coprime,
and the leading coefficients are \(\alpha h^a,\beta h^b\) with
\(d=\deg h>0\). A polynomial shear gives Bernstein degrees
\(a(g+d),b(g+d)\). In particular \(g+d\ge16\). For any member \(P\)
of a proper pair with nonconstant leading coefficient, whether the
orders divide or not,
\[
                   \operatorname{ord}P+\deg\operatorname{lc}P\ge32.
\]

**Theorem MonicPolynomialOrderThroughThirtyOneGeneration.** In the same
Weyl algebra, a pair with nonzero scalar commutator generates if either
member has nonzero constant leading coefficient and differential order
at most thirty-one. The same conclusion holds for constant-leading
orders \(33,35,37,39\).

The first remaining constant-leading order in this argument is
thirty-two. If a proper pair has such a member, subtracting powers
of it from its mate leaves orders \(32,16b\) with odd \(b\ge3\).
After depression, the relevant weight is the Bernstein weight and
the two leading symbols are powers of a monomial in two independent
linear forms. Its primitive factor has degree \(4,8\), or \(16\).
This is a necessary description of that case in this argument.
The [order-thirty-two exclusion](MONIC_POLYNOMIAL_ORDER_THIRTY_TWO_EXCLUSION.md)
now rules it out.

## A finite weighted descent

We work over an algebraic closure of the coefficient field when
needed. Proper generation is preserved by this extension, since
a nonzero vector-space quotient remains nonzero after extension
of scalars. In the filtered Ore rings below the coefficient
variable \(t\) has weight one, the derivative \(E\) has weight
\(r>0\), and \([E,t]\) has leading term \(t^{1-s}\). The loss of
one contraction is \(r+s>0\), with leading bracket
\(\{\xi,t\}=t^{1-s}\). Formal Laurent coefficients at infinity
are allowed in this calculation. Each nonzero coefficient has
an integral degree in \(t\), and \(r\) is rational, so the weights
form a discrete set.

Assume the top symbols of \(P,Q\) are \(K^{A},K^{B}\), where \(K\)
is monic of positive degree in \(\xi\), with weight \(w\).
By primitive here we mean that the multiplicities of the
distinct roots of \(K\) as a polynomial in \(\xi\), after the
fractional Laurent extension used to split it, have gcd one.
Write \(A=ea,B=eb\) with \(a,b\ge2\) coprime. We show that the
scalar commutator gives a finite element with top symbol \(F\)
satisfying
\[
                          \{K,F\}=cK^N                    \tag{1}
\]
for some \(c\ne0\) and integer \(N\ge0\).

The weighted Euler identities imply that a positive-weight
homogeneous symbol Poisson-commuting with \(K\) has a constant
ratio of powers with \(K\). Valuation at each root of \(K\)
and the gcd-one assumption make this symbol a constant
multiple of an integral power of \(K\). The Poisson factor
\(t^{1-s}\) does not affect this implication.

Start the list of finite operators with \(U_0=P, U_1=Q\)
and top exponents \(k_0=A,k_1=B\) in \(K\). Set
\(d_i=\gcd(k_0,\ldots,k_i)\), \(e_i=d_{i-1}/d_i\) for
\(i\ge1\), and
\[
        F_i^*=\sum_{\ell=1}^i(e_\ell-1)k_\ell-A,\qquad
        F_0^*=-A .
\]
Every multiple of \(d_i\) strictly greater than \(F_i^*\)
is a nonnegative combination of \(k_0,\ldots,k_i\).
For the induction, choose \(0\le v<e_i\) so the difference
after subtracting \(v k_i\) is divisible by \(d_{i-1}\).
This difference is strictly greater than \(F_{i-1}^*\),
and the case \(i=0\) is immediate.

Write \(C_iw\) for the weight of \([P,U_i]\) for \(i\ge1\).
Its leading symbol is a nonzero multiple of \(K^{C_i}\).
Let \(\delta_i=k_i-C_i\). Initially \(C_1=0\) and
\(\delta_1=B\). We maintain strictly decreasing \(\delta_i\)
and \(C_i-A=F_{i-1}^*\). The first cancellation is
\(Q^a-\gamma P^b\), with \(\gamma\) chosen to cancel its top
symbol. In later steps form \(V=U_i^{e_i}-\gamma M\),
where \(M\) is an ordered monomial in earlier operators
of exponent \(e_i k_i\). Such a monomial exists by the
semigroup bound. Indeed a new \(k_i\) is taken strictly
above \(F_{i-1}^*\), and \(e_i k_i\) is a multiple of
\(d_{i-1}\) above the same bound. The weight of \(V\)
is strictly below \(e_i k_i w\). Its commutator with \(P\)
has weight
\[
         C_{i+1}w,\qquad C_{i+1}=(e_i-1)k_i+C_i,\qquad
         C_{i+1}-A=F_i^* .
\]
In fact the commutator of \(U_i^{e_i}\) has the leading
product from the Leibniz rule with coefficient \(e_i\).
Its symbol is a nonzero multiple of \(K^{C_{i+1}}\).
The commutator of an earlier monomial of exponent \(j\)
has weight at most \((j-\delta_{i-1})w\) when \(i>1\),
and \(M\) for \(i=1\) commutes with \(P\).
This uses the strict ordering of the drops \(\delta_i\).

Preserve the leading commutator while lowering the weight
of \(V\). Its weight is at least \(F_i^*w+(r+s)\).
Strictly above this weight the leading symbol Poisson-commutes
with \(K\) and hence is a scalar multiple of \(K^j\),
where \(j>F_i^*\). If \(d_i\mid j\), subtract the matching
monomial in \(U_0,\ldots,U_i\). Its commutator weight
is at most \((j-\delta_i)w<C_{i+1}w\) since \(j<e_i k_i\).
If \(d_i\nmid j\), take the current \(V\) as \(U_{i+1}\),
with \(k_{i+1}=j\). Then \(d_{i+1}<d_i\) and
\(\delta_{i+1}=j-C_{i+1}<\delta_i\).
Each subtraction lowers a discrete weight and each new
operator lowers a positive gcd. The process must reach the
lower bound, where the leading commutator is detected by
\(\{K^A,F\}\). Dividing by \(A K^{A-1}\) gives (1), with
\(N=F_i^*+1\ge0\). Here \(F_1^*=e(ab-a-b)>0\) and the
subsequent \(F_i^*\) are nondecreasing.

## The coefficient equation

In each application, the first possible Poisson bracket
of the original pair lies strictly above the scalar.
It therefore vanishes, and the Euler identities and
monicity give top symbols \(H^a,H^b\). The common factor
\(H\) is monic of degree \(g\) in \(\xi\). Its coefficient
of \(\xi^{g-1}\) is zero by depression of \(P\), and
\(H\ne\xi^g\) because a lower term attains the defining
positive weight of \(P\). Write \(H=K^e\) with \(e\)
maximal and
\[
        K=t^{\mu r}f(\eta),\qquad \eta=\xi/t^r,\qquad
        \mu=g/e .
\]
The polynomial \(f\) is monic, depressed and has at least
two distinct roots, whose multiplicities have gcd one.
The coefficients of \(K\) belong to the same Laurent
symbol ring as \(H\), since they can be recovered by
comparing descending powers of \(\xi\) in \(K^e\) and
dividing by \(e\). In the polynomial case they are
polynomial in \(t\). This realizes the hypotheses of the
finite descent with weight \(w=\mu r\).

Write \(F=t^\theta G(\eta)\) in (1). The explicit
bracket gives
\[
        \theta f'G-\mu r fG'=c f^N,\qquad
        \theta=\mu Nr-(\mu-1)r+s .
\]
Consequently \(W=f^{1-N}G\) is a rational function satisfying
\[
               W'-\frac{r+s}{\mu r}\frac{f'}{f}W
                              =-\frac{c}{\mu r}.              \tag{2}
\]
Because \(r+s>0\), a finite pole of \(W\) would contribute
a nonzero leading term: its negative exponent minus
a nonnegative multiple of the root multiplicity cannot
vanish. Thus \(W\) is a polynomial. It vanishes at every
root of \(f\), hence has degree at least two, and comparing
leading coefficients at infinity in (2) gives
\[
                             \deg W=1+s/r.                \tag{3}
\]

## Nonconstant leading factors and Bernstein degree

The common leading-factor formula follows by taking the
coefficient of differential order \(ga+gb-1\) in the
commutator. The
[leading-factor residue theorem](POLYNOMIAL_WEYL_LEADING_FACTOR_RESIDUE.md)
rules out roots of exact multiplicity \(g\) and gives a
polynomial of degree \(d/g\) vanishing at every distinct root
of \(h\) when there are at least two such roots. Consequently
\(d=g\) is impossible, and \(d<g\) implies that \(h\) has a
single root.

Apply the polynomial shear subtracting from the original
derivative the polynomial part at infinity of \(p_{ga-1}/(ga p_{ga})\).
The resulting ratio \(p_{ga-1}/p_{ga}\) has degree at most
minus one. Choose \(x=t^\epsilon\) at infinity when \(d>g\);
choose \(x=b_0+t^\epsilon\) in the single-root case when
\(d<g\). Here \(\epsilon=g/\gcd(g,d)\) and \(v=\epsilon d/g\).
Choose \(\lambda^g=h\) and take \(E=\lambda D_x\). The
normalized coefficient field is Laurent
at infinity and the leading derivation has
\([E,t]=t^{1-s}(1+O(t^{-1}))\), where \(s=\epsilon-v\).
A constant scaling of \(E\) sets that leading coefficient
to one. Scale the operators to be monic in \(E\).
The \(E^{ga-1}\) coefficient has degree at most \(-s\):
before the harmless constant scales its two terms are
\(\lambda p_{ga-1}/p_{ga}\) and
\(-\binom{ga}{2}d\lambda/dx\). Thus depression in \(E\)
can later be reversed at any weight \(r\ge -s\)
without increasing the weights.

The monic fractional-root recursion shows that the weight of
the mate is \(gbr\) whenever the depressed first member is
bounded by \(gar\) and \(r+s\ge0\). In detail
\([P,(P^{j/(ga)})_+]\) has order at most \(ga-2\).
Successive leading coefficients in the commutator descent
are constant, so \(Q\) is a finite constant combination
of these positive parts. The recursion for each root
coefficient uses only degree-nonincreasing contractions
for the stated weights, and gives weight at most \(jr\) to
the \(j\)-th positive part, attained by its monic term.

If \(d>g\), write \(u=-s>0\). Suppose a lower coefficient
of the depressed first member has degree larger than its
deficit times \(u\). Let \(r>u\) be the maximum ratio of
degree to deficit. The first possible commutator weight is
\((ga+gb-1)r+u>0\), so the coefficient equation above
applies. It gives \(\deg W=1-u/r<1\), a contradiction.
Hence the depressed operators have weights at most \(gau,gbu\)
when the derivative has weight \(u\).

If \(d<g\), then \(s>0\) and the normalized coefficients
are Laurent polynomials. They cannot all be bounded in the
depressed first member. If they were, the same root recursion
would also bound the mate, since the derivation preserves the
bounded ring and induces zero on its residue field.
Reduction would contradict the nonzero scalar commutator.
Choose \(r>0\) by the maximum ratio of degree to deficit.
The first possible commutator weight is
\((ga+gb-1)r-s\ge0\). At equality the denominator
of \(r=s/(ga+gb-1)\) in lowest terms is an integer
\(\ell\ge(ga+gb-1)/(g-1)>a+b\), coprime to \(g\).
Neither \(ga\) nor \(gb\) is divisible by \(\ell\).
Every highest monomial has deficit divisible by \(\ell\),
so both leading symbols are divisible by \(\xi\).
Their bracket vanishes at \(\xi=0\), contrary to a nonzero
scalar at that weight. The inequality is therefore strict,
and (3) gives \(r\le s\).

Let \(R=u\) in the case \(d>g\), and \(R=r\) in the
case \(d<g\). Undo the depression in \(E\), keeping
the bounds \(gaR,gbR\). The coefficient of \(D_x^i\)
in \(E^j\) has \(t\)-degree at most \(\epsilon i-sj\).
This follows by left multiplication by \(E\): differentiating
in \(x\) lowers \(t\)-degree by \(\epsilon\), and multiplication
by \(\lambda\) increases it by \(v=\epsilon-s\).
It follows that each term contributing to \(p_i\) has degree
at most \(\epsilon i+gaR-j(R+s)\). Its maximum for \(j\ge i\)
is bounded by the value at \(j=i\), because \(R+s\ge0\).
Thus total degree is bounded by
\[
                  \deg_x p_i+i
                     \le \frac{gaR+i(2\epsilon-R-s)}{\epsilon}.
\]
The right side is strictly increasing in \(i\): in the first
case \(R+s=0\), and in the second \(R+s\le2s<2\epsilon\).
At \(i=ga\) the displayed value is \(a(g+d)\), attained by
the leading coefficient. The same argument applies to \(Q\),
with \(b\) replacing \(a\). The pair is still proper after
the polynomial shear, so
[Guccione–Guccione–Valqui, Corollary 7.4, printed page 42](sources/guccione-guccione-valqui-1111.6100v3.pdf)
gives \(d+g\ge16\).

For the asserted individual bound, keep a member of order \(i\)
with nonconstant leading coefficient fixed and descend its
partner by subtracting powers of it, as in the
[reduced order theorem](POLYNOMIAL_GCD_TWO_AND_REDUCED_ORDER_BOUNDS.md).
The partner still has positive order because the pair is
proper. If the resulting orders are nondividing, use the
bound just proved and \(i/g\ge2\) to obtain
\(\deg\operatorname{lc}P=(i/g)d\ge32-i\).
Otherwise the new partner has order \(j<i\) dividing \(i\).
Its leading coefficient is nonconstant by the first leading
commutator identity, which also gives
\(\deg\operatorname{lc}P=(i/j)\deg\operatorname{lc}Q\).
Induction on \(i\) now gives \(\deg\operatorname{lc}P\ge32-i\).

The [Bernstein degree improvement](POLYNOMIAL_WEYL_BERNSTEIN_DEGREE_THROUGH_THIRTY_THREE.md)
raises the individual order-plus-leading-coefficient bound to 34.
The bound \(g+d\ge16\) for the shared factor still describes the
nondividing reduction.

## The monic bound and the first residual

Suppose \(P\) has constant leading coefficient and order \(n\).
Reduce the mate by powers of \(P\) while possible.
The leading commutator identity makes its leading coefficient
constant too. If its order divides \(n\), it is a smaller
constant-leading member, so it is excluded by induction when
proving the asserted range. For the remaining nondividing
orders, scale to monic, depress \(P\) by a polynomial shear
and use \(t=x,\ s=1,\ [E,t]=1\). Some lower coefficient
of the depressed first member must be nonconstant by the root
recursion argument. The maximum ratio \(r\) of degree to
deficit satisfies \(r\ge1/(ga)\). The first commutator weight
\((ga+gb-1)r-1\) is positive. The preceding argument
gives (2) and (3), so \(m=1/r\) is an integer with
\(1\le m\le\mu\): a nonzero lower coefficient of \(K\)
at deficit \(j\le\mu\) has degree \(j/m\).

If \(m=1\), then \(W\) has degree two and \(f\) has exactly
two distinct roots. These give independent linear factors
of the homogeneous binary form \(K\). The Bernstein degrees
are the original differential orders \(ga,gb\), so \(g\ge16\)
for a proper pair. Since \(g\) is a proper divisor of \(n\),
this case is excluded for \(n\le31\), as well as for
\(n=33,35,37,39\). At \(n=32\) the only possibility is \(g=16\),
and the primitive degree \(\mu\) divides \(16\). The value
\(\mu=2\) would give two simple roots while the multiplier
\((r+s)/(\mu r)\) in (2) is one, which is impossible at
a simple zero of \(W\). Thus \(\mu=4,8,16\).

If \(m=\mu\), then \(K=\xi^\mu+\sigma x\). Set
\(Y=E,\ Z=-(E^\mu+\sigma x)/\sigma\), so \([Z,Y]=1\).
The strict lower-weight remainders after subtracting the
powers of \(E^\mu+\sigma x\) have orders below \(ea,eb\)
in \(Z\), since a monomial \(x^j E^i\) of weight
\(j+i/\mu<ea\) has order in \(Z\) at most \(j<ea\).
The exact new order \(ea=n/\mu\) has constant leading
coefficient, so the induction applies.

For \(2\le m<\mu\), the polynomial \(f\) has the form
\(z^p F_0(z^m)\). It has at most \(m+1=\deg W\) distinct
roots, so at most one nonzero orbit under the \(m\)-th
roots of unity. Hence \(f=z^p(z^m-\sigma)^q\) with
\(\sigma\ne0,\ q\ge1,\ p\ge1\); the last inequality
uses primitive multiplicities and \(m<\mu\).
Use \(Y=E,\ Z=E^m/\sigma-x\). In this polynomial Weyl
presentation the weights of \(Y,Z\) are \(1,m\) and
the bound on \(P\) is \(n\). Its new order \(i\) is at most
\(n/m\), and its leading coefficient degree at most \(n-mi\).
If this coefficient is constant the induction applies.
Otherwise the individual leading degree bound gives
\(n\ge(m-1)i+32\). For a proper pair \(i\ge6\) and \(i=7\)
is excluded. If \(i=6\), the reduced partner has order \(3b'\),
\(b'\ge3\) odd, and the
[gcd-three infinity bound](POLYNOMIAL_GCD_THREE_LEADING_FACTOR.md)
gives leading degree at least thirty for the order-six member.
This would require \(n\ge42\). For \(i\ge8\) one has \(n\ge40\).
Thus the mixed case is impossible throughout the asserted range.
The induction starts with the
[low-order results](LOW_ORDER_DC1_AUTOMORPHISMS.md).
For the four asserted odd orders above thirty-one, every strictly
smaller order used in the binomial or mixed change is at most
nineteen. A proper divisor in the initial mate descent is at
most thirteen, so those reductions use only the bound through
thirty-one. The same reasoning leaves at order thirty-two
only the Bernstein case described in the statement.

## Primary comparison

The weighted leading-bracket calculation is used in
[Guccione–Guccione–Valqui, Corollary 2.6](https://arxiv.org/pdf/1111.6100v3).
Their Corollary 7.4 concerns the minimum total-degree gcd of
nongenerating pairs. It is applied here only after computing
actual Bernstein degrees. Zheglov's
[arXiv:2410.06959v5, Section 4](https://arxiv.org/pdf/2410.06959v5)
uses a different recursion based on differential order and
its distance invariant. The descent above uses weighted
commutator degrees directly and has no premise from that recursion.

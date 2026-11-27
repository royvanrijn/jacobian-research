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

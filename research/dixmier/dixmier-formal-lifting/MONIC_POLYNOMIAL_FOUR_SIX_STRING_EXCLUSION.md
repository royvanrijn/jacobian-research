# Constant leading coefficient excludes polynomial orders four and six

**Theorem MonicPolynomialFourSixStringExclusion.** Let \(k\) be a
characteristic-zero field, \(D x-xD=1\), and \(P,Q\in k[x][D]\) have
differential orders four and six. If either leading coefficient is
a nonzero constant, then \([P,Q]\) cannot be a nonzero constant.

The proof uses the common square and cube of the highest positive
weight. Squaring the order-six operator and subtracting the cube of
the order-four operator gives a lower-weight operator whose leading
commutator still detects the original scalar. The resulting homogeneous
Poisson equation has no solution for the required binomial.

## Polynomial normalization and the first face

Suppose \([P,Q]=\kappa\in k^\times\). The coefficient of \(D^9\) in
the commutator is \(4p_4q_6'-6q_6p_4'\), so a constant nonzero
leading coefficient of either operator forces the other one to be
constant. Scale both to monic and retain \(\kappa\) for the nonzero
scalar after scaling. The polynomial shear \(D\mapsto D+g(x)\)
depresses \(P\). With unique \(u,v,z\in k[x]\) write
\[
 P=(D^2+u/2)^2+vD+v'/2+z .
 \tag{1}
\]
The case of constant \(v\), including zero, is the retained
[constant-skew exclusion](SKEWADJOINT_PRIME_ORDER_AND_POLYNOMIAL_46_ADMISSION.md).
Assume \(v\) is nonconstant. The same
[fractional-root descent](SKEWADJOINT_PRIME_ORDER_AND_POLYNOMIAL_46_ADMISSION.md)
gives, after subtracting a constant multiple of \(P\) and a scalar,
\[
 Q=B_6+aB_5+bB_3+cB_2+eD,\qquad B_i=(P^{i/4})_+ ,
 \quad a,b,c,e\in k .
 \tag{2}
\]

Give \(x\) weight one and \(D\) weight \(r\), choosing \(r\) to be the
largest \(\deg_x(p_i)/(4-i)\) among the nonzero lower coefficients
\(p_i\) of \(P\) in normal order. If \(u\) is nonconstant then
\(r\ge1/2\). If \(u\) is constant then \(p_1=v\), so \(r\ge1/3\).
Every product contraction lowers weight by \(1+r\). The monic
fractional root \(P^{1/4}=D+\sum_{j\ge1}w_jD^{-j}\) has all terms of
weight at most \(r\). Indeed, compare successively the coefficients of
\(D^{3-j}\) in its fourth power, where the unknown contribution is
\(4w_j\). Each earlier contribution has weight at most \(4r\),
and \(P\) itself has weight \(4r\). It follows that the polynomial
part \(B_i\) has weight exactly \(ir\). Thus \(Q\) in (2) has
weight \(6r\).

Let \(p,q\) be the leading weighted symbols of \(P,Q\) in \(k[x,\xi]\).
The coefficient at the highest possible commutator weight is
their Poisson bracket, with convention \(\{\xi,x\}=1\).
The weight in question is \(9r-1>0\), so \(\{p,q\}=0\).
The two weighted Euler identities imply that both partial derivatives
of \(q^2/p^3\) vanish. The ratio is therefore constant, and comparison
of the top \(\xi\)-coefficients gives \(q^2=p^3\). Unique factorization,
monicity, and the absence of a \(\xi^3\) term in \(p\) give
\[
 p=H^2,\qquad q=H^3,\qquad H=\xi^2+\sigma x^n,
 \qquad \sigma\ne0,\quad n=2r=\deg u .                       \tag{3}
\]
In fact the definition of \(r\) ensures that some lower coefficient
contributes to \(p\), so \(H\) has a nonzero term of \(\xi\)-degree
zero. Since \(H\) has weight \(2r\), its exponent \(n\) is an integer.
The \(\xi^2\) coefficient identifies \(u\)'s leading term.
Comparison of the remaining coefficients in (1) and (3) gives
\[
              \deg v<3n/2,\qquad \deg z<2n .                \tag{4}
\]
The derivatives in (1) have weight strictly below \(4r\) in these
comparisons. This is the usual positive-weight symbol argument,
with the descent in (2) fixing the weight of the actual mate.

## The next commutator and a homogeneous obstruction

Set \(R=Q^2-P^3\). By (3), \(\operatorname{wt}(R)<12r\), while
\([P,R]=2\kappa Q\) has weight \(6r\). We will preserve the leading
symbol of this last commutator as \(R\) is changed. The weight
estimate for any such \(R\) is
\[
            6r\le 4r+\operatorname{wt}(R)-(1+r),
 \qquad\text{so}\qquad \operatorname{wt}(R)\ge3r+1 .
 \tag{5}
\]
If the inequality is strict, the leading weighted symbol of \(R\)
Poisson-commutes with \(H\). A positive-weight homogeneous polynomial
commuting with \(H\) is a scalar multiple of \(H^j\) for an integer
\(j\). To see this, scale the rational weights to integers. The
Euler identities and the vanishing bracket yield a constant ratio
between suitable powers of the polynomial and \(H\). Since \(H\)
is squarefree, unique factorization forces the exponent of \(H\)
to be integral. The weight bounds give \(2\le j\le5\).

For \(j=2,3,4,5\), respectively, subtract the appropriate constant
multiple of \(P,Q,P^2,PQ\). Each subtraction kills the leading
symbol of \(R\). Its commutator with \(P\) has weight below \(6r\),
since the four commutators are \(0,\kappa,0,\kappa P\).
The leading term of \([P,R]\) therefore stays \(2\kappa H^3\).
The weights are discrete and bounded below by (5), so the
subtractions terminate at \(\operatorname{wt}(R)=3r+1\).
Write \(F\) for its leading symbol. At that weight the bracket gives
\[
                  \{H,F\}=\kappa H^2 .                     \tag{6}
\]

We may extend \(k\) algebraically and choose \(s\) with \(s^2=-\sigma\).
On each parametrization of the zero locus of \(H\),
\[
              x=t^2,\qquad \xi=\pm s t^n ,
\]
the derivation \(\{H,-\}\) induces
\(\pm s t^{n-1}\,d/dt\). If a polynomial of positive homogeneous
weight has bracket divisible by \(H\), its two restrictions are
constant and hence zero, since they have positive degree in \(t\)
when nonzero. Therefore \(H\) divides the polynomial. Indeed, the
remainder after division by the monic quadratic \(H\) is
\(A(x)+\xi B(x)\), and the two restrictions force both \(A\) and
\(B\) to vanish.

Apply this to \(F\) in (6), then to \(F/H\), whose bracket is
\(\kappa H\) and whose weight is \(r+1>0\). Thus \(H^2\mid F\)
and the quotient has weight \(1-r\) and bracket \(\kappa\) with
\(H\). For \(r>1\) the quotient cannot be a nonzero polynomial.
For \(r=1\) it is constant. Both cases are impossible.

## The weight \(r=1/2\)

The only remaining positive weight in (3) is \(r=1/2\).
Then \(u\) and \(v\) have degree one, and (4) gives \(\deg z\le1\).
Write \(u=U x+u_0\), \(v=V x+v_0\), with \(UV\ne0\).
For \(P,Q\) of forms (1)–(2), the coefficient of \(x^2\) in
the \(D^2\) coefficient of \([P,Q]\) is \(-15aU^3/32\).
This forces \(a=0\). With \(a=0\), the coefficient of \(x^2\)
in the \(D^1\) coefficient is \(9UV^2/4\), which is nonzero.
The [exact Ore replay](verify_monic_polynomial_four_six_string_exclusion.py)
checks both identities with arbitrary linear coefficients.
This excludes the last weight and completes the proof.

The statement concerns differential operators with actual polynomial
coefficients. For a polynomial pair with nonconstant leading
coefficients, passing to an algebraic coefficient cover to normalize
them does not supply the polynomial hypothesis used in (3)–(6).
The unrestricted \((4,6)\) problem over algebraic function fields
and actual polynomial admission with nonconstant leading coefficients
remain separate obligations.

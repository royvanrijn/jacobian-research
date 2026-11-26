# All-degree separation with private variables on both sides and comparable shared exponents

Let \(a,b\) be distinct ordinary exponent columns, with a
private variable on each side. Suppose that on their shared
variables the exponent lists are comparable coordinatewise.
After exchanging the two monomials if necessary, assume
\[
                         a_i\le b_i\qquad(a_i b_i>0).
\]
There is no order condition on the private rows. Then the
consecutive factorial moments are coprime in every degree,
over every characteristic-zero field.

The same proof applies to any positive base whose degree-\(n\)
left endpoint measure makes either of two specified conditional
expectations nondecreasing. The proof compares imaginary parts
at a proposed common *physical* root. It places no size
restriction on the exponents or the total curvature.

## A conditional monotonicity test at the physical root

Reserve one private row on each side. The
[private factorial grouping](MULTIVARIATE_PRIVATE_UNIT_ASYMMETRIC_ANGLE_BOXES.md#the-geometric-resolvent-without-reflection-symmetry)
cancels the binomial coefficient and leaves positive base
variables \(U,V\). Thus the order-\(m\) polynomial, up to a
nonzero scalar, is
\[
        P_m(z)=\mathbb E[U^m G_m(zT)],\qquad
                  T=V/U,\qquad G_m(w)=\sum_{j=0}^m w^j .
\]
The same \(U,V\) are used at every degree. For the conditional
test this representation can also start with any positive base
with all the relevant mixed moments. Denote expectation under
the tilt by \(U^n/\mathbb E U^n\) by \(\mathbb E_n\).

Suppose the pair at orders \(n,N=n+1\) has a common root.
The [angle boxes](MULTIVARIATE_PRIVATE_UNIT_ASYMMETRIC_ANGLE_BOXES.md)
exclude real roots shared between these degrees and give, after
conjugation, \(z=r e^{i\theta}\), \(r>0\), with the same positive
label \(k\) at both orders. In particular
\[
                       (2k-1)\pi<n\theta<N\theta<2k\pi .
\]
The geometric case \(T\) constant already has no common root.
Set \(q=\sin\theta>0\), \(A=-\sin(N\theta)>0\), and
\(B=-\sin(n\theta)>0\). The lower imaginary part is
\[
 I(s)=\operatorname{Im}G_n(s e^{i\theta})
             =\frac{s[q+s^n(A-Bs)]}{1-2s\cos\theta+s^2},
                                \qquad s>0 .                    \tag{1}
\]
This follows by multiplying the geometric quotient by
\(1-s e^{-i\theta}\).
The bracket in (1) is positive up to one point and
negative thereafter. Its derivative is
\(s^{n-1}(n A-(n+1)B s)\), so it first increases from
\(q>0\) and then decreases to minus infinity.

The lower root equation gives \(\mathbb E_n I(rT)=0\).
For every nondecreasing \(g\ge0\) with
\(\mathbb E_n[g(T)|I(rT)|]<\infty\) it follows that
\[
                           \mathbb E_n[g(T)I(rT)]\le0.          \tag{2}
\]
Indeed subtract a value of \(g\) at the sign change in (1);
both sides of that point give nonpositive products. Equivalently
use bounded truncations of \(g\) first and then dominated
convergence. This also handles jumps in a monotone version of \(g\).

If \(g(T)=\mathbb E_n[U\mid T]\) is nondecreasing, apply (2)
and \(G_N=G_n+w^N\). The imaginary part of
\(\mathbb E_n[U G_N(zT)]\) is strictly negative, because
\(\sin(N\theta)<0\).
If \(g(T)=\mathbb E_n[V\mid T]\) is nondecreasing, use
\(w^{-1}G_N(w)=G_n(w)+w^{-1}\) instead. The upper root equation
would give
\[
             0=z^{-1}\mathbb E_n[U G_N(zT)]
                =\mathbb E_n[V G_n(zT)]+z^{-1}\mathbb E_n U .
\]
Its imaginary part is again strictly negative by (2) and
\(\operatorname{Im}z^{-1}<0\). Both contradictions use the
same physical \(z\) at the two degrees. The integrability in
(2) follows from the order-\(N\) mixed moments.

## Why comparable shared ordinary rows satisfy the test

After reservation, write the base as
\(U=C\prod_i Z_i^{a_i}\), \(V=C'\prod_i Z_i^{b_i}\)
with positive constants and independent Gamma inputs.
The original unreserved inputs are exponential, and the
reserved rows contribute the residual Gamma factors of
shapes \(h/d\), \(1\le h<d\), with exponent one on the
same private side. Shared equal rows can be integrated
separately. Under \(U^n\) the retained inputs are still
independent gamma, with their shapes increased by \(n a_i\).

For each unequal residual row put \(e_i=b_i-a_i\) and
\(S_i=e_i\log Z_i\). Its density under this tilt, up to a
positive constant, is
\[
                 \exp\bigl(\alpha_i s/e_i-e^{s/e_i}\bigr),
                            \quad -\infty<s<\infty,
\]
where \(\alpha_i>0\) is its tilted shape.
The second logarithmic derivative is \(-e^{s/e_i}/e_i^2<0\),
so every \(S_i\) is log-concave, including at negative \(e_i\).
The conditional law of the vector at fixed \(T\) is its law
at fixed \(\sum_i S_i=\log T-\log(C'/C)\).
In the chosen orientation, \(b_i/e_i\ge0\) for every unequal
row. Positive \(e_i\) gives a nonnegative coefficient, and
negative \(e_i\) occurs only on the private \(U\) side, where
\(b_i=0\). Thus \(V\), as a function of the \(S_i\), is
coordinatewise nondecreasing.

[Efron's conditional monotonicity theorem](https://pmc.ncbi.nlm.nih.gov/articles/PMC4847755/)
now makes \(\mathbb E_n[V\mid T]\) nondecreasing in \(T\).
The conditional test proves the assertion. Positive fixed
Gamma shapes and nonnegative real exponents also work if
the two reserved variables have unit exponent and unit
shape. For the ordinary theorem the coefficients are rational,
so the rational Bézout identities give field descent.

## A proof of the conditional monotonicity used here

For completeness, take two independent variables \(X,Y\) with
positive smooth log-concave densities \(f,g\) on the line.
For two sums \(s<t\), the conditional density ratio for \(X\)
is, up to a scalar, \(g(t-x)/g(s-x)\), which is nondecreasing
in \(x\) by log concavity. Thus the conditional law of \(X\)
stochastically increases with the sum. The same argument
applies to \(Y\). Couple the conditional \(X\) values at
\(s,t\) by the same uniform quantile. The \(X\) values are ordered.
The \(Y\) values, recovered by subtracting from their sums,
are also ordered: they are the conditional \(Y\) quantiles
at the reflected uniform quantile. Every coordinatewise
nondecreasing observable therefore increases under this
coupling.

To iterate, the density of \(X+Y\) is log-concave.
For these smooth densities this also follows by differentiation:
the logarithmic derivative of \(f*g\) at \(s\) is the conditional
expectation of the nonincreasing function \((\log g)'(Y)\).
It is nonincreasing by the conditional order just proved.
The log-Gamma densities are bounded with integrable derivatives,
and their convolutions have the same properties, which
justify the differentiation. Group the first two variables
by their sum. The expected observable conditioned on this
partial sum remains nondecreasing in that sum and in the
ungrouped variables. It is independent-group conditional
expectation, so induction proves the general claim.
Truncation handles the nonnegative unbounded observable
used above.

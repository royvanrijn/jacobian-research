# Equal-order cancellation for regular scalar Long repairs

Work in the [exact-first-image scalar Long chart](ALL_COPRIME_ALGEBRAIC_STRING_OBSTRUCTION.md#6-a-specialization-obstruction-without-a-unit-leading-coefficient)
with its regularity and constant-uniformizer hypotheses. Use the prime-order
specialization and coprime-order theorems proved in that source.

## Equal-order cancellation

Suppose a regular finite scalar Long repair has
\(\operatorname{ord}S=\operatorname{ord}T=n\) in the generic chart. Cancel the leading terms and scale by a differential constant
as needed to get a regular operator. The resulting operator has a
composite generic order \(r<n\) with \(\gcd(r,n)>1\). In particular,
\((4,4)\) is excluded, and at \((6,6)\) the cancellation has order four.

Indeed, write \(s_n,t_n\) for the leading coefficients. The coefficient of
\(D^{2n-1}\) in \([S,T]\) is
\(n(s_n\delta t_n-t_n\delta s_n)\), so \(c=t_n/s_n\) is a differential
constant. If \(c\) is regular, take \(U=T-cS\) and use the pair \(S,U\).
Its commutator is \(h\) and its classical orders are \(2,3\), because
\(\overline T\) has higher order than \(\overline S\).
If \(c\) has a pole, take \(U=c^{-1}T-S\) and use the pair \(U,T\).
It has commutator \(-h\), and \(\overline U=-\overline S\) so the
classical orders are again \(2,3\). In both cases \(U\) is regular and its
generic order \(r\) is at least its classical order and strictly below
\(n\). If \(r\) were prime, the prime-order specialization theorem
applied to the pair with \(U\) as its prime-order member would force
\(3\mid 2\) or \(2\mid 3\). Thus \(r\) is composite, and
CoprimeAlgebraicStringExclusion gives \(\gcd(r,n)>1\).

Combining this cancellation with CompositeLongScalarOrders shows that
the smaller of the two generic orders is at least four and the larger
is at least six. If the smaller order is four, the larger is even.
The [selfadjoint power-of-two exclusion](SELFADJOINT_POWER_OF_TWO_STRING_CLASSIFICATION.md)
also rules out the stated algebraic selfadjoint normalization of that
order-four member. These restrictions concern regular finite scalar
repairs with the exact retained first image.

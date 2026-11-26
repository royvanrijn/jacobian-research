# Native Newton edges certify first-order angular separation under weak shared tails

Fix an ordinary factorial binomial core with arbitrary shared and
private rows. Let \(Q_m\) be its moment polynomial, and consider
adjacent orders \(n,N=n+1\), \(n\ge2\). Append any finite nonempty
collection of unequal shared rows of either orientation. Measure the
size of the tail by
\[
 W=\sum_{\rm tail}\frac{(b_i-a_i)^2}{(a_i+b_i)/2}.
\]
The results below permit tails of arbitrary dimension and signed
degree difference.

**Theorem.** Suppose that at a common prime \(p>2nN\) every edge of
both core Newton polygons has horizontal width one. Then the two
augmented moments have disjoint root arguments for every such
ordinary tail with \(0<W<W_0\), for an effective \(W_0>0\) depending
on the core and the order pair. The core itself may be tested
without any assumption of angular separation. There is also a
criterion below for arbitrary binomial edges of width prime to \(p\).

Consequently, let \(J\) be any integer at least the
[absolute ordinary angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md).
If each \(2\le n<n+1<J\) for a single core has either disjoint
arguments already or the stated native polygon certificate, one
common positive tail tolerance gives separation at every order.
This uses the effective but unevaluated absolute cutoff. At orders
one and two the even moment is strictly positive on the real axis.
Positive radial rescaling of either augmented polynomial is allowed.

## Angular response and a local formula at a Newton edge

Put \(D=z\,d/dz\). For a simple nonzero core root \(\alpha\) at
order \(m\), the
[matched shared-tail expansion](MULTIVARIATE_MATCHED_TANH_SIMPLE_CORE_ANGULAR_RESPONSE.md#the-simple-root-response)
gives
\[
 \log\frac{\alpha(W)}{\alpha}
   =W\,\rho_m(\alpha)+O_{Q_m,m}(W^2),\qquad
 \rho_m(\alpha)=\frac12-\frac{1}{2m}
                        \frac{D^2Q_m(\alpha)}{DQ_m(\alpha)} .
 \tag{1}
\]
Here \(\alpha(W)\) is the nearby root after a positive
normalization of the tail endpoints. The remainder is uniform in
the row list. Equivalently the leading log coefficient multiplier
is \(-Wk(m-k)/(2m)\), with all other terms \(O_m(W^2)\); expand
the simple root to obtain (1). The imaginary part is the angular
change to first order.

Embed the common splitting field of the core pair in an algebraic
closure of \(\mathbb Q_p\). For a Newton edge \(E=[a,b]\),
write \(s(E)\) for the negative of its slope and \(u(E)=a+b\).
Assume its residual polynomial is binomial, with \(p\nmid b-a\).
If a root \(\alpha\) has valuation \(s(E)\), we claim
\[
              \frac{D^2Q_m(\alpha)}{DQ_m(\alpha)}
                                  \equiv a+b\pmod{\mathfrak m} ,
 \tag{2}
\]
where \(\mathfrak m\) is the maximal ideal of the valuation ring.
Indeed, scale the variable by an element of valuation \(s(E)\),
and divide the polynomial by a scalar to make its minimum
coefficient valuation zero. The reduction is \(c_a x^a+c_bx^b\)
with both coefficients nonzero, and the root reduces to a
nonzero \(x_0\). At \(x_0\),
\(c_b x_0^b=-c_a x_0^a\), and the reduction of \(DQ_m\) in
the same scale is \((a-b)c_a x_0^a\ne0\).
The reduction of \(D^2Q_m\) is \((a^2-b^2)c_a x_0^a\).
This gives (2) and also proves that such a root is simple.

Suppose a lower root \(\alpha\) and an upper root \(\gamma\)
share an argument. Real arguments cannot be shared, so their
conjugates are different roots of their own polynomials. The
equality of arguments gives the exact equality
\(\alpha/\bar\alpha=\gamma/\bar\gamma\).
For a root \(\xi\) let \(E_\xi\) be its edge in the common
embedding. Equality of angular derivatives in (1) would also give
\(\rho_n(\alpha)-\rho_n(\bar\alpha)
 =\rho_N(\gamma)-\rho_N(\bar\gamma)\).
By multiplying through by \(2nN\) before reducing, (2) yields
the following two necessary conditions for a shared argument
with equal first-order angular response:
\[
 \begin{aligned}
 s(E_\alpha)-s(E_{\bar\alpha})
       &=s(E_\gamma)-s(E_{\bar\gamma}),\\
 N\bigl(u(E_\alpha)-u(E_{\bar\alpha})\bigr)
       &\equiv n\bigl(u(E_\gamma)-u(E_{\bar\gamma})\bigr)\pmod p .
 \end{aligned}                                                     \tag{3}
\]
Thus when all edges in the two polygons are binomial and their
widths are prime to \(p\), it suffices to eliminate (3) over
ordered pairs of distinct edge roots. One may first eliminate
pairs with unequal residual quotients after clearing the same
denominators of all valuations, as in the
[paired Newton test](MULTIVARIATE_PAIRED_MULTI_EDGE_ANGULAR_GATE.md#a-residue-comparison-when-slope-differences-coincide).
The response comparison permits coincident residual quotients.

With width one, the two ordered distinct roots correspond to
different edge indices \(k,l\), and \(u([k-1,k])-u([l-1,l])
=2(k-l)\). The second condition in (3) implies
\[
                        N(k-l)-n(i-j)\equiv0\pmod p
\]
for distinct upper indices \(i,j\). The integer on the left is
nonzero because \(n,N\) are coprime and \(0<|k-l|<n\).
Its absolute value is at most \(N(n-1)+n(N-1)<2nN<p\).
This contradiction proves that every coincident argument has
different first-order angular responses. The proof even
distinguishes the first-order logarithmic speeds of any shared
off-diagonal root quotient between the two core degrees.

There are finitely many root pairs at the fixed orders. A
positive angular gap persists for any pair with distinct initial
arguments. For a coincident pair the uniform \(O(W^2)\) remainder
is smaller than half the linear angular difference whenever \(W\)
is sufficiently small. Hence all pairs separate. The tolerance
can be computed using algebraic root isolation and Taylor bounds
on simple root disks, as in the cited response expansion.
Taking the minimum of the finitely many local tolerances below
\(J\) proves the all-order consequence.

## An exact control for the derivative test

The positive rational polynomials
\[
 (z+101^2)(z^2+z+101),\qquad
 (z^2+101^3z+101^7)(z^2+z+101)
\]
share a nonreal quadratic. Their \(101\)-adic root valuations
are respectively \(0,1,2\) and \(0,1,3,4\), so both Newton
polygons have width one. The even polynomial is strictly positive
on the real line. Thus the same response argument separates
their root arguments under any sufficiently small coefficient
perturbation with the uniform first-order multiplier in (1).
This is an algebraic check of the arithmetic response test on
a pair with an exact initial contact. The control pair is
outside the asserted ordinary factorial core class.

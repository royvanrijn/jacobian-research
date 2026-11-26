# Polynomial height templates with arbitrary multiplicities

An ordinary two-monomial profile may be specified by finitely
many row types whose smaller exponents \(\ell_i(M)\) are fixed
polynomials over \(\mathbb Q\) and whose positive integer gaps
\(d_i\) and orientations are fixed. Assume the endpoints are
ordinary nonnegative integers at the values of the positive
integer parameter \(M\) under consideration. Each row type
has an independently chosen nonnegative integer multiplicity.
Equal rows can be dropped. Additional retained rows are
allowed, with arbitrary heights, gaps and orientations.

For orders \(n,N=n+1\), strip each row's smaller endpoint
factorial. Its possible literal factors at order \(j=n,N\)
are \(j\ell_i(M)+h\), \(1\le h\le j d_i\). A polynomial
factor \(F\in\mathbb Z[M]\) is a separating divisor for type
\(i\) and order \(m\in\{n,N\}\) if:

- \(F\) has at least two distinct complex roots and divides
  \(m\ell_i(M)+h_0\) over \(\mathbb Q[M]\) for some
  \(1\le h_0\le d_i\).
- \(F\) is relatively prime to every other literal factor
  of the template at the two orders.

Occurrences from all copies of the designated factor in
type \(i\) count as one literal factor for this definition.

**Theorem.** Fix a separating divisor and write \(r_i\ge1\)
for the multiplicity of its type. There is a computable
finite prime set \(S\), determined by the template and
orders, and a computable \(\kappa>0\) such that for all
sufficiently large \(M\), a physical common root forces
\[
 \sum_{\substack{p\notin S\\p\mid F(M)}}\
       \sum_{\rho\ {\rm retained}} v_p(L_\rho)\log p
                         \ \ge\ \kappa r_i\log M.              \tag{1}
\]
Here \(L_\rho\) is the least common multiple of the full
stripped ramps at both orders for row \(\rho\).
This bound involves the actual retained rows. In particular
they can be arbitrary fixed ordinary rows, by including
their primes in \(S\). The resulting exact exclusion holds
uniformly over all template multiplicities, orientations
and signed total gaps.

If at every paired order below the absolute angular cutoff
there is a separating divisor of a present type, a common
computable height bound gives ordinary SFC(2) in every
degree with any fixed retained core. For any fixed list of
gaps and polynomial degrees with at least one designated
degree at least two, the separating-divisor requirement on
that designated type holds outside a proper algebraic
subset of polynomial coefficient space. The designated
type must have positive count, and every other count is
arbitrary. In particular this gives all-degree coprimality
above a computable height for generic polynomial height
templates with arbitrary independent multiplicities, without
any signed balance or diffuse curvature assumption.

## Exact prime isolation and the retained charge

Take \(S\) to include the prime factors of the denominators
of all rational polynomial coefficients, those of the
orders' factorials, and the primes of the nonzero rational
resultants of \(F\) with the remaining factors after clearing
denominators. Also include coefficient denominators of the
quotient of the target factor by \(F\). Add \(2\) if needed.
At any prime \(p\notin S\) dividing \(F(M)\), every other
template ramp factor at either order is a unit. Indeed their
polynomial reductions have no common root with the reduction
of \(F\). The target factor has valuation at least \(v_p(F(M))\)
by the polynomial divisibility and the included denominators.
Thus the selected coefficient valuation at the other order
is zero everywhere. Orient the two columns, if needed, so
the marked type is increasing. At order \(m\), index zero
has valuation zero, while every positive index has valuation
at least \(r_i v_p(F(M))\) because the marked factor occurs
in its first block.

For the retained rows use the
[common valuation slope removal](MULTIVARIATE_RETAINED_COMMON_SLOPE_BUDGET.md#arbitrary-retained-rows-at-one-prime).
After one common linear translation their ordinates in each
order have oscillation at most \(R_p=\sum_\rho v_p(L_\rho)\).
At the unmarked order all the translated Newton slopes
therefore have magnitude at most \(R_p\). At the marked
order, after subtracting any of those slopes times the
index, each positive-index ordinate exceeds the index-zero
ordinate by at least
\[
                       r_i v_p(F(M))-(m+1)R_p.
\]
If this is positive, index zero is the unique minimizer,
so there is no shared slope and hence no common factor.
A common root therefore forces
\[
        \sum_{p\notin S,\ p\mid F(M)} R_p\log p
           \ \ge\ \frac{r_i}{m+1}
                          \log\!\left(\frac{|F(M)|}{[F(M)]_S}\right),
                                                               \tag{2}
\]
where \([F(M)]_S\) denotes the product of its prime powers in
\(S\). No full ramp product of a retained row is needed in
this test.

The effective bound on \(S\)-parts from
[Bugeaud, Evertse and Győry, Theorem 2.2](https://arxiv.org/pdf/1708.08290)
states that for an integer polynomial with at least two
distinct complex roots and any fixed finite nonempty \(S\),
\([F(M)]_S\le K |F(M)|^{1-\nu}\) whenever \(F(M)\ne0\), with
computable \(K,\nu>0\). Inserting it in (2) proves (1), since
\(\log|F(M)|=\deg F\log M+O(1)\). A fixed retained core
contributes zero to the left side once its primes are in
\(S\). Low orders supplied by even positivity and the
[absolute angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md)
complete the all-degree statement. The argument is rational
and therefore gives coprimality over any characteristic-zero
field by a Bézout identity.

## The generic polynomial templates

For a designated nonlinear type \(i\), take as \(F\) the
primitive integer multiple of \(N\ell_i(M)+1\) at each low
paired order. Squarefreeness of this polynomial and
coprimality with the factors from the other types are finite
nonvanishing discriminant and resultant conditions. The
other factors within type \(i\) are automatically relatively
prime to it: their affine shifts in \(\ell_i\) are different
since \(1/N\ne h/n\) for \(h\ge1\) and \(1/N\ne h/N\) for
\(h\ne1\).

To see that the exceptional algebraic set is proper, fix
the coefficients of positive degrees. Vary the constant
term of \(\ell_i\) to avoid the finitely many critical
values giving a multiple zero in \(N\ell_i+1\), over all
the low paired orders. Then fix \(\ell_i\) and vary the
constant term of each other polynomial to avoid its
finitely many values causing a ramp factor to vanish at
one of those zeros. Positive rational choices are available
in every sufficiently large interval. With a nonlinear
\(F\) of degree at least two, the hypotheses follow.
The same generic construction can designate several
nonlinear types by choosing constants successively, so
one of them can serve whenever it is present.

For an exact example take a present type with smaller
exponent \(M^2\) and any fixed positive gap. Let every
other moving type have a linear smaller exponent in \(M\)
with positive rational leading coefficient and an arbitrary
fixed positive gap. Then \(N M^2+1\) has no real zeros, hence
no common zeros with their linear factors. It is already
relatively prime to the other factors within the marked row.
The result gives all-degree SFC(2) uniformly over every
independent choice of these counts and with any fixed
retained ordinary core.

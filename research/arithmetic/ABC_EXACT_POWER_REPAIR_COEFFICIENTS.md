# Coefficient records and the effective orders of exact power repair

Let \(a_1+a_2=a_3=c\) be a positive primitive triple. Fix
\(q_1,q_2,q_3\ge1\), set \(Q=\operatorname{lcm}(q_1,q_2,q_3)\), and
write
\[
 g_i=\gcd\!\left(q_i,\operatorname{lcm}_{j\ne i}q_j\right),\qquad
 R_{i,e}=\prod_{\substack{p\mid a_i\\v_p(a_i)=e}}p,\qquad
 m_e=\min_i\log R_{i,e}.
\]
We allow an arbitrary positive integer common multiplier \(m\) and
positive integer coefficient and root choices satisfying
\[
                       ma_i=s_i x_i^{q_i}\qquad (i=1,2,3).
 \tag{EC1}
\]
The height of the repaired exact equation is \(mc\). Paying a coefficient
record by its radical charges \(\log\operatorname{rad}(s_1s_2s_3)\).
For fixed orders these two costs have an exact joint minimum.

## The minimum joint cost with coefficient records

When \(g_i\mid e\), let \(r_i(e)\) be the unique solution in \([0,Q)\)
of
\[
        r_i(e)\equiv-e\pmod {q_i},\qquad
        r_i(e)\equiv 0\pmod {q_j}\quad(j\ne i).
\]
Put \(w_i(e)=r_i(e)/Q\) if \(g_i\mid e\), and \(w_i(e)=1\) otherwise.
Then
\[
 \begin{split}
 &\min_{m,\boldsymbol s,\boldsymbol x\ {\rm satisfying}\ (EC1)}
 \left(\frac{\log(mc)}{Q}+\log\operatorname{rad}(s_1s_2s_3)\right)\\
 &\hspace{35mm}=\frac{\log c}{Q}+
                    \sum_{i,e}w_i(e)\log R_{i,e}\\
 &\hspace{35mm}\ge \sum_{e\ge1}m_e .
 \end{split}                                                     \tag{EC2}
\]
All choices are allowed to depend on the original triple, including
the orders. The equality and inequality hold for each such choice
of orders before any counting limit.

Indeed, primitive support assigns each prime of \(a_1a_2a_3\) to one
letter \(i\) with valuation \(e\). Outside the coefficient radical,
\(v_p(m)\) must satisfy the displayed congruences. Generalized CRT
solves them exactly when \(g_i\mid e\). If solvable, the least local
cost is \(r_i(e)\log p/Q<\log p\), attained with \(p\) absent from
all three coefficients. Putting \(p\) in the coefficient radical
would cost at least \(\log p\). If unsolvable, the prime must be
in that radical. The least cost is then \(\log p\), attained by
setting \(v_p(m)=0\), putting \(e\bmod q_i\) in \(v_p(s_i)\), and
putting \(\lfloor e/q_i\rfloor\) in \(v_p(x_i)\). At a prime outside
\(a_1a_2a_3\), setting all exponents to zero minimizes its cost.
All these local choices combine to give the first equality.

Fix \(e\ge1\). If any \(g_i\nmid e\), then \(\sum_i w_i(e)\ge1\).
If all divide \(e\), then \(\sum_i r_i(e)\equiv-e\bmod Q\), since
this holds modulo each \(q_j\). Nonnegativity of the residues
gives \(e/Q+\sum_i w_i(e)\ge\lceil e/Q\rceil\ge1\).
Thus \(e/Q+\sum_i w_i(e)\ge1\) in either case. Since
\(\log c=\sum_e e\log R_{3,e}\ge\sum_e e m_e\), retaining \(m_e\)
in each of the three nonnegative costs gives the last inequality.

More generally, the joint cost of any single repair is bounded
below by the right side of the equality in (EC2). The exact minimum
records the charges of all prime powers beyond the common overlap.
An added common \(Q\)-th power can only increase the objective.

In base \(X\), let \(\mu_i(e)=\log_X R_{i,e}\), and let
\(L=\sum_{i,e}\mu_i(e)\). The pointwise inequality
\(\sum_i\mu_i(e)-3\min_i\mu_i(e)\le
\sum_{i<j}|\mu_i(e)-\mu_j(e)|\) yields the explicit lower bound
\[
 \frac{\log_X(mc)}{Q}+\log_X\operatorname{rad}(s_1s_2s_3)
 \ \ge\
 \frac{L-\sum_{i<j}\|\mu_i-\mu_j\|_1}{3}.                    \tag{EC3}
\]
This sector includes profile models with essentially uniform
weights on all degree positions in a compact scaled band. Weak
convergence to the same limiting measure alone does not ensure
the microscopic overlap in (EC3).

The unconditioned \(Y^{1/Q}\) benchmark for pairwise coprime orders
already holds on canonical repairs of primitive quotients.
Consequently combining a coefficient-radical enumeration with an
unconditioned weighted-height-only bound at that benchmark cannot
produce the strict endpoint saving on the microscopic balance
sector. This conclusion concerns precisely that addition of the
two costs. A count exploiting the original radical caps, the joint
distribution of coefficient records and repaired heights, or the
additive relation within a restricted record would have more
information than (EC2) uses.
[The canonical repaired-height benchmark](ABC_EXACT_POWER_REPAIR_OBSTRUCTION.md)

## Effective orders for every repair with coefficient records

There is an exact arithmetic alternative to the common multiplier.
For any finite set \(S\) of primes, a repair (EC1) with
all coefficient primes in \(S\) exists if and only if all
\(v_p(a_i)\) with \(p\notin S\) are divisible by \(g_i\).
Equivalently the primitive quotient satisfies
\[
                a_i=t_i y_i^{g_i},\qquad
                \operatorname{supp}(t_i)\subseteq S.        \tag{EC4}
\]
For necessity apply CRT to each prime outside \(S\).
For sufficiency use the same local minimizing choices outside
\(S\), and set \(v_p(m)=0\) on \(S\) while taking the ordinary
remainders of the original valuations modulo \(q_i\).
For (EC4), transfer to \(t_i\) either the full \(S\)-part of \(a_i\)
or only its valuation remainders modulo \(g_i\).
The prime support minimal among valid coefficient records is
exactly the set of original primes with \(g_i\nmid v_p(a_i)\).
Darmon and Granville already use the orders \(g_i\) in their
reduction for general diagonal equations. Here (EC4) retains the
exact allowed coefficient support and applies before height counting.
[Darmon–Granville, §1.1, p. 519](https://www.math.mcgill.ca/darmon/pub/Articles/Research/12.Granville/pub12.pdf#page=7)

Hence common repair with small coefficient radical and high
shared orders already implies a primitive high-divisibility grid
with no greater coefficient-radical cost. The generalized Fermat
bound for large orders applies directly to that grid. Increasing
the \(q_i\) using factors that leave the \(g_i\) unchanged adds no
arithmetic divisibility condition to (EC4).
[Fine-grid count](ABC_FINE_VALUATION_DIVISIBILITY.md#exact-sector-bound)

## The effective orders determine the weighted Fermat curve

Let \(C_{\boldsymbol q}\) be the normalization of
\[
       X_1^{q_1}+X_2^{q_2}=X_3^{q_3}
          \quad\text{in}\quad
       \mathbf P(Q/q_1,Q/q_2,Q/q_3).
\]
The morphism to the line \(a_1+a_2=a_3\) given by taking the three
powers has degree \(D=q_1q_2q_3/Q\). Over \(\mathbf C\) the cover
is Galois and has ramification indices \(g_i\) over \(a_i=0\),
with index \(1\) allowing an unramified point. There is no
ramification elsewhere. In particular
\[
       2\,\operatorname{genus}(C_{\boldsymbol q})-2
                         = D\left(1-\sum_{i=1}^3\frac1{g_i}\right).
                                                               \tag{EC5}
\]
Moreover \(C_{\boldsymbol q}\) and \(C_{\boldsymbol g}\) are
isomorphic over \(\mathbf Q\), as covers of the line.

Here is a direct proof. Over \(\mathbf C\) we may set \(X_3=1\)
on the part with \(X_3\ne0\). The affine curve
\(X_2^{q_2}=1-X_1^{q_1}\) is irreducible, by the valuation at a
simple zero of the right side. The weighted relation on this
chart then quotients by the remaining \((Q/q_3)\)-th roots of
unity. This quotient acts freely when all three coordinates
are nonzero, because \(\gcd_i(Q/q_i)=1\). The generic fibre thus
has cardinality \(q_1q_2/(Q/q_3)=D\). The action of the full
diagonal group on this fibre identifies the deck group with
\[
  (\mathbf Z/q_1\mathbf Z\ \times\ \mathbf Z/q_2\mathbf Z\
       \times\ \mathbf Z/q_3\mathbf Z)/
                       \langle(1,1,1)\rangle .
\]
On a disc around \(a_i=0\) choose one of the nonzero projective
coordinates as gauge. Analytic continuation of the roots around
that point gives the class of the \(i\)-th unit vector. Its order
is the least positive \(k\) congruent to a multiple of
\(\operatorname{lcm}_{j\ne i}q_j\) modulo \(q_i\), namely \(g_i\).
Away from the three coordinates' zeros all roots have local
analytic branches. Riemann–Hurwitz now gives (EC5).

For every prime \(\ell\), sort the three valuations of the orders
as \(\alpha\ge\beta\ge\gamma\). The valuations in the corresponding
positions of \(\boldsymbol g\) are \(\beta,\beta,\gamma\).
Writing \(Q'=\operatorname{lcm}(g_1,g_2,g_3)\), it follows that
\(g_1g_2g_3/Q'=D\). The power map
\[
                       Y_i=X_i^{q_i/g_i}
\]
defines a morphism to the weighted curve of orders \(\boldsymbol g\):
the weights of the three monomials are \(Q/g_i\), proportional
to the target weights \(Q'/g_i\). It commutes with the maps to the
line. Both normalized curves have degree \(D\) over the line.
Thus the induced morphism between the normalizations has degree
one, and it is an isomorphism. The map and its inverse on the
normalizations are defined over \(\mathbf Q\).
The same proof holds with fixed nonzero rational coefficients
in the diagonal equation, using the line with those coefficients.

For pairwise coprime original orders the effective orders are
\((1,1,1)\), the cover has degree one, and the curve is rational.
Positive genus requires actual shared valuation divisibility in
the primitive quotient away from the coefficient support.
Darmon and Granville already give the rationality statement in
the pairwise coprime case and distinguish it from the
primitive generalized Fermat theorem. The effective-order
description above supplies the full cover for arbitrary orders.
[Darmon–Granville, p. 515](https://www.math.mcgill.ca/darmon/pub/Articles/Research/12.Granville/pub12.pdf#page=3)

The obstruction in (EC2) is a pointwise repair-cost inequality.
It neither asserts that microscopically balanced primitive
triples saturate the small-budget count nor bounds their
distinct targets after two affordable divisor records. Those
counts retain the original three radical constraints.

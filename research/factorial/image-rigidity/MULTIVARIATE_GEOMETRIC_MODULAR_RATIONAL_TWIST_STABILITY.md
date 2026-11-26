# All-degree modular row constructions stable under positive rational rescaling

Let \(Q_m(z)\) be the ordinary factorial moment polynomial of
two distinct exponent columns at order \(m\). The constructions
below ensure that \(Q_n(z)\) and \(Q_{n+1}(\kappa z)\) are coprime
for every \(n\ge1\) and every positive rational \(\kappa\).
They allow fully shared profiles, unrestricted curvature,
and unbounded numbers and magnitudes of exponent rows.
A suitable choice of the prime also permits every positive
real algebraic \(\kappa\) of odd degree over \(\mathbb Q\).

Fix \(J\ge4\) at least the
[computable absolute angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md).
Put \(L=\operatorname{lcm}(1,\ldots,J-1)\) and take a prime \(p>J\)
such that no odd prime at most \(J\) divides \(p-1\).
For example one can take \(p\equiv-1\) modulo the product of
the odd primes at most \(J\), by
[Dirichlet's theorem on primes in arithmetic progressions](https://sites.math.rutgers.edu/~alexk/2022S572/Dirichlet1837.pdf).
The cutoff is unevaluated.

Two ordinary constructions at any such prime are as follows.

- Start with two private unit rows, one for each column, and
  use any of the
  [neutral modular row operations](MULTIVARIATE_MODULAR_ROW_OPERATIONS_ALL_DEGREES.md)
  at this \(p\), including common shifts, mixed append blocks,
  independent row inflation and digit packing.
- Use complete exponent-shift orbits as described below.
  Take one marked orbit with gap \(1\bmod L\) and one with
  gap \(-1\bmod L\). Every further orbit can have any gap
  divisible by \(L\). These congruences allow either orientation
  even for the marked orbits. In particular the gap magnitudes,
  their signed sum, and their residues modulo \(p-1\) are
  unrestricted except for the displayed congruences modulo \(L\).

The proof uses reductions which are geometric sums after
possibly different unit scalings at the two degrees. This
permits a weaker gap modulus for the digit-shift orbits than
when the reductions are required to agree with the untwisted
geometric sums.

## The geometric criterion modulo one prime

Write \(G_m(x)=1+x+\cdots+x^m\). Suppose there is a common
nonzero rational \(d\) and, for \(1\le m<J\), a nonzero
rational \(c_m\) such that \(F_m(x)=c_m Q_m(dx)\) is
\(p\)-adically integral coefficientwise, with
\[
                     \overline F_m(x)=a_mG_m(b_m x),
                     \qquad a_m,b_m\in\mathbb F_p^\times.      \tag{1}
\]
The factors \(b_m\) need not agree.

For an adjacent pair with \(N=n+1<J\), the nontrivial
\((n+1)\)-th and \((n+2)\)-th roots of unity in
\(\overline{\mathbb F}_p\) have no ratio in \(\mathbb F_p^\times\).
Indeed, if they have respective orders \(e>1,f>1\), then
\(\gcd(e,f)=1\) and the order of their ratio is \(ef\).
One of \(e,f\) divides an odd integer among \(n+1,n+2\),
so \(ef\) has an odd prime factor at most \(J\). By the
condition on \(p\), such an order cannot divide \(p-1\).
Conversely, if an odd prime \(q\le J\) divides \(p-1\), then
at \(n=q-2\) a nontrivial \(q\)-th root and \(-1\) give a
ratio in \(\mathbb F_p^\times\). Thus the prime condition is
exact for simultaneous separation of these geometric
reductions under every unit rescaling.

All roots of \(F_m\) in \(\overline{\mathbb Q}_p\) have
valuation zero, by (1). The constant term would uniquely
dominate at positive valuation and the leading term would
uniquely dominate at negative valuation. If
\(F_n(x)=F_N(\kappa x)=0\) for rational nonzero \(\kappa\),
both arguments must have valuation zero, so \(\kappa\) is
a \(p\)-adic unit. Reducing gives two roots of
\(G_n,G_N\) whose ratio lies in \(\mathbb F_p^\times\),
since \(\bar\kappa,b_n,b_N\) all belong to that group.
This is impossible. The common \(d\) in the normalization
makes this exactly the comparison of the original \(Q_n\)
and \(Q_N(\kappa x)\). The proof at \(N<J\) covers every
nonzero rational \(\kappa\), including those with nonzero
\(p\)-adic valuation. At \(N\ge J\) the absolute theorem
gives the comparison at every positive \(\kappa\).
Rational coprimality extends to all characteristic-zero fields.

If \(p\equiv-1\) modulo every odd prime at most \(J\), the same
finite-field ratio exclusion holds in \(\mathbb F_{p^f}\) for
every odd \(f\). Indeed \(p^f-1\) is divisible by none of these
odd primes. Let \(\kappa\) be a real algebraic number of odd
degree, and use the number field \(K=\mathbb Q(\kappa)\).
At least one prime of \(K\) above \(p\) has odd residue degree,
since the sum of the products of the ramification indices and
residue degrees is \([K:\mathbb Q]\). Apply the same valuation
argument at that prime. If \(\kappa\) is a unit, its residue
lies in a field of odd degree, so again no root ratio can
equal its residue up to a factor in \(\mathbb F_p^\times\).
For positive \(\kappa\), the absolute theorem continues to
cover \(N\ge J\).

## Constructing the reductions

**Starting from private units.** The original core has
\(Q_m=m!G_m\) and \(p>J\), so its normalized reductions are
\(G_m\) at every relevant order. The coefficient calculations
for the neutral operations preserve these reductions after
scalar normalization and a common diagonal substitution.
They apply at the specified \(p\) because the seed
reductions \(G_m\) already have their full degree and both
endpoints units. The shift and append bounds and the
digit-packing separation are exactly those of the cited
coefficient calculations.
Explicitly, shifting the two core rows gives
\((pu,pu+1),(pv+1,pv)\) with independent nonnegative integers
\(u,v\). Taking both positive gives a fully shared seed with
no requirement that the smaller exponents coincide.

**Starting from complete orbits.** An orbit of base entries
\((A,B)\) and length \(P=p^t\) uses the \(P\) rows
\[
                         (A+i,B+i),\qquad 0\le i<P.
\]
Choose \(P\equiv1\bmod L\) and \((J-1)\max(A,B)<P\).
Such powers exist because \(p\) is coprime to \(L\).
One may take \(A,B\ge1\) in every orbit and choose them
independently across orbits. With signed gap \(e=B-A\),
the [exact digit-orbit reduction](MULTIVARIATE_FACTORIAL_DIGIT_ORBITS.md#exact-orbit-reduction)
contributes the factor
\[
 H_m(ke),\qquad
 H_m(t)=(m!)^{\lfloor t/m\rfloor}(t\bmod m)!,
\]
in the \(k\)-th coefficient in degree \(m\), after normalization.
Here \(0\le t\bmod m<m\), and negative powers of \(m!\) are
interpreted modulo \(p\). For a collection of orbits the
argument substitution in this normalization is
\(x\mapsto(-p)^{-D}x\), where
\(D=\sum e(P-1)/(p-1)\). It is common to every degree.

If \(e\equiv e_0\pmod L\) and \(m<J\), then \(m\mid(e-e_0)\)
and
\[
                    H_m(ke)=
                    (m!)^{k(e-e_0)/m} H_m(ke_0).              \tag{2}
\]
In particular all the changes of gap by multiples of \(L\)
contribute only powers of units of the form \(b_m^k\).
For the two marked gaps use \(e_0=1,-1\). Their other factors
satisfy \(H_m(k)H_m(-k)=k!(m-k)!/m!\) for \(0\le k\le m\).
They therefore cancel the binomial coefficient exactly.
The remaining gap-zero factors contribute only the units
in (2). Thus the reduction is \(G_m(b_m x)\) for some unit
\(b_m\), as in (1). No positivity or coprimality test over
the finite field remains after these calculations.

One can also retain exactly one private unit row and use
only the oppositely oriented marked orbit to supply the
other geometric factor. Extra zero-class orbits can be
appended to any profile satisfying (1) at the same prime:
their common diagonal substitution combines with the
original one and their units \(b_m^k\) preserve (1).

For every profile certified by (1) the
[integer-twist exceptional set](MULTIVARIATE_ARBITRARY_CORE_TWISTED_TAIL_STABILITY.md)
is empty. Its small-curvature stability conclusion therefore
applies after appending arbitrary shared rows, including
tails with either sign of the total exponent gap. This
last corollary needs a small-curvature hypothesis for the
unstructured tail. The modular constructions themselves do
not use such a hypothesis.

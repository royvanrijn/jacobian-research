# Multiplicative spectator units from radical dictionaries

Products of independent draws from a sufficiently large
dictionary of small integers have a nearly uniform unit
distribution at most auxiliary primes. Their Fourier
moments are bounded by the ordinary character large sieve.
This gives an arithmetic supply of independent unit slots
for the projective tree comparison, subject to keeping
those slots independent of the original incidence data.

## Uniform mixing at auxiliary primes

**MSU1.** Fix \(L>0\) and \(n=X^L\), and fix
\(0<a\le h<\beta\) and an integer \(k\ge1\).
Let \(\mathcal D\subset[1,n^h]\) be a dictionary
of \(N\ge n^a\) distinct integers. Put \(Q=n^\beta\)
and let \(\mathcal P\) be any set of primes in
\([Q,2Q]\) with \(|\mathcal P|\ge n^{\beta-o_X(1)}\).
Write \(\mu_q\) for the law of a uniform draw
from \(\mathcal D\) reduced to
\(\mathbb F_q^\times\), and
\[
             \widehat\mu_q(\chi)=
                \frac1N\sum_{x\in\mathcal D}\chi(x).
\]
All elements of the dictionary are units by \(\beta>h\).
Then
\[
 \mathbf E_{q\in\mathcal P}
       \sum_{\chi\ne1}|\widehat\mu_q(\chi)|^{2k}
    \ll_{k,L,h,\beta}
        n^{o_X(1)}
        \frac{Q^2+n^{kh}}{|\mathcal P|\,N^k}
    \le n^{-\delta+o_X(1)},\qquad
       \delta=ka-\max\{\beta,kh-\beta\}.
                                                               \tag{MSU1}
\]
If \(\delta>0\), then for every fixed
\(0<\gamma<\delta\), outside an
\(O(n^{-\delta+\gamma+o_X(1)})\) proportion
of these primes, the law of the product of
\(2k\) independent draws satisfies
\[
 \sup_{u\in\mathbb F_q^\times}
       \left|(q-1)\mu_q^{*2k}(u)-1\right|
                      \le n^{-\gamma}.
                                                               \tag{MSU2}
\]
Here convolution is multiplicative. The bounds
are uniform in the choice of the dictionary and
the subset of test primes.

Write \(b_j\) for the number of ordered
\(k\)-tuples in \(\mathcal D^k\) with product \(j\).
Its length is at most \(n^{kh}\), and
\(\sum_j b_j^2\le N^k n^{o_X(1)}\) by
the fixed-order divisor bound. Moreover
\[
 N^k\widehat\mu_q(\chi)^k=\sum_j b_j\chi(j).
\]
Every nonprincipal character modulo a prime
is primitive. Apply the character large sieve,
in its prime-modulus consequence of the Gauss
identity and the spacing of the fractions \(t/q\)
in Montgomery and Vaughan's Theorem 1. This
bounds the summed squares by
\(O((Q^2+n^{kh})\sum_j b_j^2)\).
Division by \(|\mathcal P|N^{2k}\) gives
(MSU1). The nonprincipal part of multiplicative
Fourier inversion for \(\mu_q^{*2k}\) has
absolute value at most
\(\sum_{\chi\ne1}|\widehat\mu_q(\chi)|^{2k}\)
at every unit. Markov's inequality gives
(MSU2).
[Montgomery–Vaughan, Theorem 1](https://personal.science.psu.edu/rcv4/personal/Publications/large_sieve.pdf)
· [The prime-character consequence used in (MF5)](ABC_MODULAR_JACOBI_INCIDENCE.md#higher-moments-on-an-actual-factorized-letter)

## Actual radical dictionaries and spectator slots

At the balanced occurrence endpoint let \(d=1/3\)
in base \(n\), and take a complete-label matching
with \(n^{d-o(1)}\) rows, with its original caps
and repair cutoff retained. The dictionary of
distinct radicals of any one role has size
\(n^{d-o(1)}\), since a fixed radical of a
complete letter has \(X^{o_X(1)}\) representatives
at polynomial height. Those radicals lie in
\([1,n^{d+o(1)}]\). For any fixed \(\beta>d\)
choose an integer \(k\) with
\[
                           \beta<kd<2\beta .
\]
The interval for \(k\) has length \(>1\).
For sufficiently small balance allowances,
(MSU1) gives a fixed power of uniformity
at most test primes. No Cartesian description
of the complete occurrence dictionary is used.

The same construction gives nearly uniform
products of support radicals raised to an order
\(M\) fixed at each budget. Raising to \(M\)
pushes uniform measure on the units to uniform
measure on their \(M\)-th-power subgroup.
A supremum bound on the relative density
persists under this pushforward by averaging
over fibres. If \(M\) is odd, one may instead
restrict the test primes to a fixed reduced
progression so that \(\gcd(M,q-1)=1\). Impose
\(q\equiv2\pmod p\) for every prime \(p\mid M\),
and \(q\equiv1\pmod2\). The prime number theorem
in this fixed progression gives an eligible
prime set of size \(n^{\beta-o_X(1)}\). The
power map is then a permutation of the full
unit group. The same argument applies to
any fixed finite list of odd orders using
their least common multiple. For even orders
the subgroup observation alone does not
supply full unit slots.

Take any fixed collection of independent
spectator slots, each constructed from a
product of independent draws as in (MSU2).
The dictionaries and the chosen numbers of
draws can differ by slot. A finite union
of the exceptional prime sets leaves a
common collection of good primes, and the
joint relative density on the independent
unit space is \(1+O(n^{-\gamma})\) for a
fixed positive \(\gamma\), chosen below the
smallest applicable exponent. Assign those
units to edges of the bipartite overlap graph
of two tree diagrams and take the incident
products for their leaves. Independent
uniform units push forward to the uniform
law on the prescribed paired-product subgroup,
as in (AT1). The supremum density bound
persists on this subgroup by averaging over
preimages. It yields an \(O_\ell(n^{-\gamma})\)
signed error in (AT2) at fixed depth \(\ell\).
[Overlap law and signed stability (AT1–AT2)](ABC_SIGNED_AUXILIARY_TREE_COMPARISON.md#the-tree-estimate-and-the-law-of-shared-slots)

The averaging uses independent draws from
one-letter dictionaries. It loses the
correlations between the three supports of
a genuine row, and a product of draws need
not itself be a support that occurs. A joint
support-pair product dictionary can already
be at the square-root threshold of (MSU1):
cardinality exponent \(a=d\) and length
exponent \(h=2d\) give
\(kd-\max(\beta,2kd-\beta)\le0\) for every
\(\beta,k\). Joint projective ratios would
require a separate estimate. To use the
spectator supply in occurrence, an arithmetic
history must introduce these slots independently
while retaining its true primitive identities
and a positive statistic on all three
capped letters. This supply alone does not
provide either part of that transfer.
[Two-sided signed tree comparison](ABC_TWO_SIDED_AUXILIARY_TREE_CANCELLATION.md)

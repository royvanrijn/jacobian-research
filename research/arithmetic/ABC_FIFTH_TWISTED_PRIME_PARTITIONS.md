# Genus compatibility and prime partitions for fifth-twisted quadratics

For each quadratic field in the pair
\(\mathbf Q(\sqrt d),\mathbf Q(\sqrt{5d})\), put two oriented split
primes in the same part if their ideal classes have the same image
modulo \(\ell\)-torsion. These parts supply coherent norm witnesses.
The two partitions can have many edges separately while sharing no
edge. An intersection of Hilbert class fields describes the joint
restriction on the prescribed classes and gives realizations by actual
split primes.

## The common class field

**Proposition FJ2.** Let \(d\ne1\) be a signed squarefree integer
prime to 5, \(K_1=\mathbf Q(\sqrt d)\) and
\(K_2=\mathbf Q(\sqrt{5d})\). Write \(H_i\) for the narrow Hilbert
class field, \(F_i\) for its genus field, and \(L=K_1K_2\).
Then
\[
             (H_1L)\cap H_2=F_2=F_1L .                            \tag{1}
\]
Thus among pairs of oriented ideal classes at rational primes
splitting in \(L\), the only compatibility on the two Artin symbols
is equality on this common genus field, after requiring the
symbols to fix \(L\). In particular every pair of square classes
in the two narrow class groups occurs in every sufficiently large
dyadic prime interval for fixed \(d\). Here a square class means
an element of the subgroup \(2\operatorname{Cl}^+(K_i)\)
in additive notation.

The fundamental discriminants are \(D_1\) and \(5D_1\), since 5 is
prime to \(d\) and is 1 modulo 4. The prime discriminants of \(5D_1\)
are those of \(D_1\) together with 5. The genus-field description by
their square roots gives \(F_2=F_1(\sqrt5)=F_1L\), so \(H_2\)
contains \(L\). The genus field is the maximal subfield of \(H_i\)
abelian over \(\mathbf Q\), with fixed subgroup
\(2\operatorname{Cl}^+(K_i)\).
[Stevenhagen, Theorem 2.1 and equations (5)–(8)](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/128E98F4B96AB15B77CDF63E6E27D447/S0305004121000335a.pdf/redei-reciprocity-governing-fields-and-negative-pell.pdf)

The extension \(L/K_1\) ramifies above 5, since \(K_1/\mathbf Q\)
is unramified there and adjoining \(\sqrt5\) is ramified.
Thus \(H_1\cap L=K_1\), and
\[
 \operatorname{Gal}(H_1L/\mathbf Q)
   \simeq \operatorname{Gal}(H_1/\mathbf Q)\times
          \operatorname{Gal}(\mathbf Q(\sqrt5)/\mathbf Q).
\]
Its commutator is \(2\operatorname{Cl}^+(K_1)\), embedded into
\(\operatorname{Gal}(H_1L/L)\), and its fixed field is \(F_1L\).

Let \(E=(H_1L)\cap H_2\). It is Galois over \(\mathbf Q\).
The group \(\operatorname{Gal}(E/L)\) is a quotient of both
\(\operatorname{Gal}(H_1L/L)\) and \(\operatorname{Gal}(H_2/L)\).
Choose an automorphism over \(\mathbf Q\) fixing \(K_2\) and
acting nontrivially on \(K_1\). Conjugation by a lift to the
compositum inverts \(\operatorname{Gal}(H_1L/L)\) and acts
trivially on \(\operatorname{Gal}(H_2/L)\). This follows from
the inversion of ideal classes under quadratic conjugation,
and from the abelian group over \(K_2\). Therefore
\(\operatorname{Gal}(E/L)\) has exponent at most two. Its quotient
map from the first group kills the commutator of
\(\operatorname{Gal}(H_1L/\mathbf Q)\), so \(E\subset F_1L\).
The reverse inclusion was already given by the genus fields.

Two automorphisms of \(H_1L,H_2\) over \(L\) extend to the
compositum exactly when they restrict to the same automorphism
of \(E\). In the first factor a class of \(K_1\) extends
uniquely while fixing \(L\); in the second factor the class
of \(K_2\) must fix \(L\). Square classes fix \(F_i\), hence
satisfy both requirements. Chebotarev applied to the conjugacy
class of the resulting automorphism gives a positive density of
rational primes, all split in \(L\). Choosing a prime of the
compositum with the specified Frobenius orients the primes of
both quadratics. A finite list of desired class pairs can
therefore be realized by distinct rational primes in each
sufficiently large dyadic interval. The bound on that interval
depends on the fields.
[Lenstra, Chebotarev density theorem](https://websites.math.leidenuniv.nl/algebra/Lenstra-Chebotarev.pdf)

## Realizing independently prescribed partitions

**Corollary FJ3.** For any fixed odd \(\ell\ge1\) and finite
bipartite multigraph, there are arbitrarily large negative
prime radicands \(d\) prime to 5 with the following property.
In every sufficiently large dyadic prime interval one can
choose a distinct rational prime for each graph edge, with
oriented prime ideals in \(K_1,K_2\), such that their two
\(\ell\)-torsion coset partitions are exactly the incidences
with the left and right endpoints, respectively. Inside each
part the coherently chosen generators give primitive norm
witnesses with trace and coefficient bounds as in the
[coherent pair construction](ABC_COHERENT_PRIME_PATH_HIGHER_ODD_TORSION.md#coherent-witnesses-with-a-bounded-unit-logarithm).

Take \(d\equiv1\bmod44\) with \(-d\) prime. Infinitely many
such primes exist by Dirichlet's theorem. The prime 11 splits
in both quadratics since \(d\equiv1\bmod11\) and 5 is a square
modulo 11. Choose one ideal \(\mathfrak q_i\) over 11 in
each field. Let \(B\) be the larger of the two vertex counts.
If \(|d|>4\cdot11^{2B\ell}\), the classes of
\(\mathfrak q_i^{2j}\) for \(0\le j<B\) are distinct modulo
\(\operatorname{Cl}(K_i)[\ell]\). Otherwise a positive power
\(\mathfrak q_i^{2s\ell}\) for \(1\le s<B\) would have a
generator \((u+v\sqrt{d_i})/2\) with \(v\ne0\), where
\(d_i=d\) or \(5d\). The principal ideal of a rational
number is invariant under conjugation, which explains the
nonzero coefficient. Its norm would give
\(4\cdot11^{2s\ell}=u^2+|d_i|v^2\ge |d_i|\), a contradiction.

Assign these distinct square classes to the vertices on each
side and choose the pairs of prime ideals by (1) and
Chebotarev. At two distinct primes in the same part the
ratio of the oriented prime ideal classes has order dividing
\(\ell\), so the ideal
\((\overline{\mathfrak p}_i\mathfrak p_j)^\ell\)
is principal and has a bounded generator as in the cited
construction. In an imaginary field the norm gives the
size bound directly. Conversely, if this ideal is principal
then the class ratio lies in \(\operatorname{Cl}(K_i)[\ell]\),
so the primes lie in the same part.

In particular, a simple graph realizes separate witness parts
without any shared edge in the two graphs on these oriented
primes. Complete bipartite graphs make the two part sizes
arbitrarily large. To impose the same parts even when prime
orientations are allowed to vary independently, take the square-class
labels with indices \(1,\ldots,B\) and strengthen the radicand
condition to \(|d|>4\cdot11^{4B\ell}\). The same norm argument
then separates the chosen labels from all their inverses in
the torsion quotient. Any primitive norm witness between two
selected primes uses some orientations and would require their
classes to agree up to inversion in that quotient. This makes
all shared edges absent on the selected primes when the graph
is simple, even with orientations left free.

These are fixed-field realizations with no polynomial bound
for the dyadic intervals as \(d\) varies. They show that a
shared-edge lower bound cannot follow from the two part sizes
and the common genus characters alone. An argument at the
[fifth-power moment scale](ABC_FIFTH_TWISTED_NORM_JOINS.md#effect-at-the-fifth-power-moment-scale)
would have to use the scale itself or another joint constraint.

## Exact four-prime control

Use \(d=-263\). The two fields have class groups \(C_{13}\) and
\(C_6\), so both fifth-torsion subgroups are trivial.
For each field use the integral basis
\(1,\omega_i=(1+\sqrt{d_i})/2\), and write the oriented prime
as \((p,\omega_i-r_i)\). The Artin class coordinates are
relative to ideal class generators
\((2,\omega_1-1)\) and \((7,\omega_2-1)\), respectively.

| \(p\) | \(r_1\) | \(r_2\) | Class in \(C_{13}\) | Class in \(C_6\) |
| --- | --- | --- | --- | --- |
| 1031 | 315 | 27 | 3 | 0 |
| 1439 | 515 | 125 | 3 | 2 |
| 1319 | 556 | 659 | 4 | 0 |
| 1201 | 85 | 114 | 4 | 2 |

The classes follow using `bnfcertify` and `bnfisprincipal` in
PARI/GP on \(x^2-x+66\) and \(x^2-x+329\), for example
`bnfisprincipal(b,idealhnf(b,1031,x-315),0)` in the first field.
The witnesses \((u+v\sqrt{d_i})/2\) are given below, with the
conjugate orientation at the first prime and the designated
orientation at the second. Each has norm \((pq)^5\) and
\(\gcd(v,pq)=1\). The traces and coefficients can be replayed
using the norm and the two reductions at the displayed roots.

| \(d_i\) | \((p,q)\) | \(u\) | \(v\) |
| --- | --- | --- | --- |
| \(-263\) | \((1031,1439)\) | 4293720844298182 | 198044126022888 |
| \(-263\) | \((1319,1201)\) | -4416839867651664 | 278466226246450 |
| \(-1315\) | \((1031,1319)\) | -1746258128548089 | -108754926617975 |
| \(-1315\) | \((1439,1201)\) | 7834403877633839 | -14922688763405 |

The two matchings have no common edge. The same is true on
these four primes after allowing changes of orientation:
the only inverse pairs among the listed class labels occur
inside the original equal-label parts.

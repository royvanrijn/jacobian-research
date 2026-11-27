# Full five-Selmer moments in the negative twist family of X_0(11)

Frey's class-group comparison transfers the quadratic fifth-torsion fourth
moment to the second moment of full five-Selmer groups. The transfer is
classical; the improved exponent comes from the retained
[higher odd torsion moment theorem](ABC_COHERENT_PRIME_PATH_HIGHER_ODD_TORSION.md#signed-tails-and-moments).
All claims here have internal written-proof status. They inherit that
input's assurance limits and make no claim of independent review or priority.

## The curve and the exact twist conditions

Use the model
\[
 E=X_0(11):\quad y^2+y=x^3-x^2-10x-20,
 \qquad P=(5,5).
\]
This is Cremona 11a1, LMFDB 11.a2, rather than every curve in its isogeny
class. Its conductor is 11 and its invariants are
\[
 c_4=496,\qquad c_6=20008,\qquad
 \Delta=-11^5,\qquad j=-2^{12}31^3/11^5.
\]
These data and the model identity are recorded in the
[LMFDB entry](https://www.lmfdb.org/EllipticCurve/Q/11/a/2).
The group law gives \(2P=(16,-61)\), \(3P=(16,60)\),
\(4P=(5,-6)\), and \(5P=O\). The reduction at 5 is good and
\(P\) reduces to an affine point, so it is outside the reduction kernel.
At 11 the model is minimal: \(c_4\) is a unit and \(v_{11}(\Delta)=5\).
The reduction is multiplicative and split, since \(-c_6\equiv1\pmod{11}\).
Thus the conductor exponent there is one; all other primes have good
reduction. In particular \(v_5(j)=0\).

Define the family using signed squarefree radicands:
\[
 \mathcal D_{11}=
 \left\{d<0:\ d\text{ squarefree},\ (d,55)=1,
                       \left(\frac d{11}\right)=-1\right\}.       \tag{1}
\]
An explicit twist model is
\[
 E^{(d)}:\quad y^2=x^3-13392d^2x-1080432d^3.
\]
Only the rational isomorphism class of this model is used.

Frey's theorem, printed p. 651, imposes \(d\equiv3\pmod4\) only when
\(2\mid N_E\). There is no such condition here. His condition at the
torsion prime is also vacuous because \(v_5(j)=0\); either nonzero
residue symbol at 5 is allowed. At the split multiplicative prime 11
his condition is exactly the symbol in (1). His exceptional set
\(\widetilde S_E\) is empty: the only bad prime is 11 and
\(11\not\equiv-1\pmod5\). The corollary on printed p. 652 therefore gives
\[
 h_5(d)\ \mid\ S(d)\ \mid\ h_5(d)^2,
 \quad h_5(d)=|\operatorname{Cl}(\mathbb Q(\sqrt d))[5]|,
 \quad S(d)=|\operatorname{Sel}_5(E^{(d)}/\mathbb Q)|.              \tag{2}
\]
Here Frey's notation denotes the full multiplication-by-five Selmer group,
not a degree-five isogeny Selmer group; see its definition on p. 650.
[Frey, theorem and corollary](sources/frey-selmer-twists-1988.pdf)

In particular (1) includes even radicands and both odd classes modulo 4.
For example \(d=-1,-3,-14,-26\) satisfy (1). Fundamental discriminants
are \(D=d\) if \(d\equiv1\pmod4\), and \(D=4d\) otherwise; hence
\(|d|\le |D|\le4|d|\). The counting parameter below is \(|d|\).
The same exponents hold if it is replaced by \(|D|\).
The family has size \(\asymp X\): the upper bound is immediate and the
lower bound follows, for example, from squarefree integers in the reduced
class \(d\equiv-1\pmod{220}\), by the elementary squarefree sieve.

## Transferred moments

**Corollary (Full five-Selmer moments).** For \(X\ge2\) and every
\(\varepsilon>0\),
\[
 \sum_{\substack{d\in\mathcal D_{11}\\ |d|\le X}} S(d)^2
                  \ll_\varepsilon X^{53/26+\varepsilon}.          \tag{3}
\]
The same argument at adjacent orders gives
\[
 \sum_{\substack{d\in\mathcal D_{11}\\ |d|\le X}}S(d)^{3/2}
       \ll_\varepsilon X^{7/4+\varepsilon},\qquad
 \sum_{\substack{d\in\mathcal D_{11}\\ |d|\le X}}S(d)^{5/2}
       \ll_\varepsilon X^{5/2+\varepsilon}.                       \tag{4}
\]
Indeed (2) gives \(S(d)^s\le h_5(d)^{2s}\). Enlarge the nonnegative
sum to all negative squarefree radicands, and apply the retained moment
theorem at orders \(2s=3,4,5\). For the central order its four terms
have exponents \(2,53/26,9/5,2\), whose maximum is \(53/26\).
There is no equidistribution or density assumption on large torsion in
the restricted family in this argument.

## Large Selmer groups and Sha[5]

**Corollary (Selmer and Sha tails).** Uniformly for \(X\ge2,Y\ge1\),
\[
 \#\{d\in\mathcal D_{11}:|d|\le X,\ S(d)\ge Y\}
   \ll_\varepsilon
      \min\{X,\ X^{53/26+\varepsilon}Y^{-2}\}.                    \tag{5}
\]
The exact sequence
\[
 0\longrightarrow E^{(d)}(\mathbb Q)/5E^{(d)}(\mathbb Q)
 \longrightarrow\operatorname{Sel}_5(E^{(d)}/\mathbb Q)
 \longrightarrow\Sha(E^{(d)}/\mathbb Q)[5]\longrightarrow0
\]
implies the identical bound with \(S(d)\ge Y\) replaced by
\(|\Sha(E^{(d)}/\mathbb Q)[5]|\ge Y\). This concerns finite five-torsion
and requires no finiteness conjecture for the whole Tate-Shafarevich group.
The exact sequence is also displayed in
[James--Ono, equation (1.5)](sources/james-ono-selmer-quadratic-twists.pdf).

To prove (5), each counted summand of (3) is at least \(Y^2\).
For either the five-Selmer dimension or the dimension of \(\Sha[5]\),
the count with dimension at least an integer \(t\ge0\) is consequently
\[
               \ll_\varepsilon
                 \min\{X,\ X^{53/26+\varepsilon}5^{-2t}\}.         \tag{6}
\]
For a fixed \(\eta>0\), (5) gives
\(O_\eta(X^{1-\eta})\) twists with size at least
\(X^{27/52+\eta}\), by taking \(\varepsilon=\eta\).
Thus (5) alone gives a power saving over the family size only above
the threshold exponent \(27/52\). It gives no near-linear second
moment. Using (4) as well adds the terms
\(X^{7/4+\varepsilon}Y^{-3/2}\) and
\(X^{5/2+\varepsilon}Y^{-5/2}\) to the minimum.

## Comparison with the same family and quantity

The classical antecedents already supply the bridge to this family.
The improvement in (3) is a consequence of the improved class-group
moment, with the following precise comparisons.

| Inspected primary statement | Consequence for this application |
|---|---|
| [Frey (1988), pp. 650--652](sources/frey-selmer-twists-1988.pdf), theorem and corollary | Supplies (2) for (1), with constant one, and explicitly names \(X_0(11)\). The Selmer/class-group bridge is not new. |
| [James--Ono (1999), Definition 11 and Proposition 12](sources/james-ono-selmer-quadratic-twists.pdf) | Restates this same double divisibility for the good-twist conditions. Its theorem on trivial Selmer groups (Theorem 14(ii)) has an additional good-at-five hypothesis; it addresses small Selmer groups rather than an upper second moment. No specialization of that additional hypothesis is claimed here. |
| [Česnavičius, arXiv:1301.4724v4, Example 1.12, equations (2)--(3)](sources/cesnavicius-selmer-flat-1301.4724v4.pdf) | Studies this exact curve after number-field base change. For imaginary quadratic \(K\) with 11 inert, it bounds the full five-Selmer dimension between \(2r_5(K)-1\) and \(2r_5(K)+1\). It supplies another classical comparison, not the exponent in (3). |
| [Heath-Brown--Pierce, arXiv:1409.3177v2, Theorem 1.5](sources/heath-brown-pierce-1409.3177v2.pdf) | At \(\ell=5,k=4\), the high branch is \(1+3(2/5)=11/5\). Applying exactly (2) gives \(\sum_{\mathcal D_{11},|d|\le X}S(d)^2\ll_\varepsilon X^{11/5+\varepsilon}\). Equation (3) lowers this transferred exponent by \(21/130\). At class orders three and five the corresponding exponents are \(9/5\) and \(13/5\), compared with (4). |
| [Koymans--Thorner, arXiv:2309.07204, Theorem 1.1](https://arxiv.org/pdf/2309.07204) | For quadratic fields and \(\ell=5,r=4\), its bound with \(O(X)\) fields is \(X^{2+1-4/6+\varepsilon}=X^{7/3+\varepsilon}\). Bounded conversion between radicand and discriminant leaves that exponent unchanged. |
| [Bhargava--Shankar, arXiv:1312.7859, Theorems 1--2](https://arxiv.org/pdf/1312.7859) | The average full five-Selmer size is six over all elliptic curves ordered by height, including families defined by finitely many coefficient congruences. A fixed quadratic-twist family is thin in that ordering, so those statements cannot be restricted to (1) to obtain a moment bound. |

The Česnavičius comparison can also be related directly to twists:
restriction and corestriction for the quadratic extension, with 2
invertible modulo 5, identify
\(\operatorname{Sel}_5(E/K)\) with the direct sum of the Selmer groups
of \(E\) and \(E^{(d)}\) over \(\mathbb Q\). This follows on global
cohomology and on each local Kummer image by the two idempotents
\((1\pm\sigma)/2\). Example 1.12 applied to \(K=\mathbb Q\)
gives \(\dim\operatorname{Sel}_5(E/\mathbb Q)=1\), using both its
bound and parity formula. Thus its imaginary-quadratic bound also gives
\(h_5(d)^2/25\le S(d)\le h_5(d)^2\) on (1). In particular even a
two-sided comparison up to a fixed constant is already available from
this existing work. The moment transfer itself carries no novelty claim.

The checked conclusion is improvement over the explicit
Heath-Brown--Pierce/Frey second-moment bound on the very same family.
These inspections do not establish historical priority among all upper
bounds on five-Selmer moments or tails. No conclusion about all negative
twists of every rational-five-torsion curve follows without checking
their local hypotheses and exceptional ray-class factors.

## Evidence scope

The primary PDFs above are retained with source hashes in
[the structured catalogue](../MATH_STATUS.json). The source inspection
covered Frey's definition, theorem, lemma and corollary; James--Ono's
exact sequence, good-twist conditions and double divisibility; and the
entire Česnavičius example, including its base-change domain. The two
moment comparison theorem statements were checked at the displayed
parameters. The application uses the written torsion moment theorem;
no new Selmer-group computation or torsion experiment is asserted.
Exact rational arithmetic checks of the displayed curve invariants,
torsion multiples and exponent substitutions are finite consistency
checks, not independent replay of the underlying arithmetic theorem.

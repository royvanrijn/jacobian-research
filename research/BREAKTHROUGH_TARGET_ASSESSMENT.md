# Breakthrough target assessment

The objective remains a result with major external mathematical impact. None of
the five unrestricted targets below is resolved by the current repository.
Primary target: JC2. Secondary alternative: every-target C12+C12 coverage.
The other candidates are audited here, not pursued as parallel family-expansion
programmes. Assessment date: 18 November 2026.

## Plane Jacobian conjecture: primary target

**Quantifiers.** For every characteristic-zero field k and every P,Q in k[x,y],
does Jac(P,Q) in k* force k[P,Q]=k[x,y]? One verified noninvertible pair
disproves it.

**External frontier.** There is no single strongest partial theorem: component
degree, generic degree, symmetry and boundary hypotheses are distinct.
Domrina's published general four-sheet theorem excludes generic degree four;
the primary statement and limitations are recorded in
[the plane A4 source audit](jacobian/monodromy-theory/PLANE_A4_MONODROMY_EXCLUSION.md).
Gao's November revision of
[arXiv:2608.00222](https://arxiv.org/html/2608.00222v2) explains actual
counterexamples in dimensions at least three and explicitly leaves the plane
open. This is a preprint. Borisov--Gabber--Vasiu's September
[arXiv:2609.05746v1](https://arxiv.org/html/2609.05746v1), Theorem 28.6,
constructs plane counterexamples to the prime-to-characteristic formulation
over every algebraically closed positive-characteristic field. These do not
by themselves settle characteristic zero.

**Repository evidence.** Polynomial higher-dimensional compilers, finite-fibre
realizations and numerous scoped plane exclusions are retained. Their
dimensions and hypotheses matter. The current work establishes a
[uniform radial lifting mechanism](jacobian/mixed-characteristic-lifting/BGV_ODD_PRIME_LIFTING.md):
the selected odd-prime BGV maps have a compatible determinant-one polynomial
lift modulo every p^n, but the constructed exact p-adic lift is provably
nonpolynomial. Its truncation degrees are necessarily unbounded. The proof
is uniform; finite exact replays only control its formulas.

**Single obstruction on this route.** Construct a bounded-degree polynomial
lift through a different deformation, or prove a universal obstruction to
such algebraization. Higher congruence compatibility alone supplies neither.
The radial construction is excluded by an exact function-field degree
comparison, not merely by a failed coefficient search.

**Impact test.** An actual plane pair or a universal impossibility theorem
would settle the remaining characteristic-zero conjecture. The radial
lifting theorem is useful route evidence, not achievement of that impact.
Further digits or optimal finite-lift degrees are not the next objective.

## Restricted-prime sums: secondary alternative

**Quantifiers.** Let C={n>=1: every prime divisor of n is 1 mod 12}, including 1.
Does there exist N0 such that every integer N>=N0 with N=2 mod 12 has
N=a+b with a,b in C? Requiring N/3<=a,b<=2N/3 is a stronger separate target.

**External frontier.** The directly inspected published inputs include
Loughran--Matthiesen's JEMS 2024 Frobenian correlation theorem, BGLM's GAFA
2023 signed dispersion estimate, Iwaniec's half-dimensional sieve and
Henriot's corrected upper bound. They give substantial averaged or
one-sided tools. The exact full theorem of Iwaniec--Pomykala (1993),
*Sums and differences of quartic norms*, remains an explicit priority and
hypothesis-audit obligation: the accessible primary extract is insufficient
to claim the strongest historical binary theorem. This literature gap
must not be converted into a novelty assertion.

**Repository evidence.** Density-one coverage, quantitative exceptional-set
bounds, local positivity and one-sided C+A3 pointwise results are retained.
The [critical switching criterion](arithmetic/ASYMMETRIC_CRITICAL_SWITCHING.md)
bounds bad survivors but still assumes an unproved signed or absolute
distribution estimate. Current local verifiers passed separately performed
internal replay; they do not verify that analytic estimate.

**Single obstruction.** Prove the signed remainder bound in that criterion
for every admissible target, preserving arbitrary target valuations.
Expanding its companion divisor sum produces weighted combined moduli near
N^(3/4). BGLM's unweighted large-modulus average and small auxiliary modulus
do not cover these terms. Fouvry--Tenenbaum's alternative weighted interface
also misses the needed parameter ranges. Merely absorbing the sieve into
a finite integer-character convolution fails: tau*1_C has half-integer,
rather than integer, Dirichlet-character exponents. Fractional coefficients
retain medium bilinear boxes but lose the existing periodic fallback.

**Impact test.** An unconditional every-target theorem would resolve the
stated binary frontier. More density-one estimates, local identities or
large finite searches do not remove its obstruction.

## Meridional rank

**Quantifiers.** For every classical knot K, is mr(K)=b(K)? The immediate
universal frontier is mr(K)<=3 => b(K)<=3. A certified triple generating
the original knot group together with an independent bridge-four bound
would be a counterexample.

**External frontier.** Boileau--Jang--Weidmann,
[Pacific J. Math. 292 (2018)](https://doi.org/10.2140/pjm.2018.292.61),
prove the three-meridian implication when the double branched cover is a
graph manifold. Recent positive knot families are not a general solution.
The November [Gao--Gu--Li--Wang Whitehead preprint](https://arxiv.org/abs/2611.15673)
is already cited in the repository and remains distinct from a refereed
publication. The detailed [frontier audit](meridional-rank/representations/THREE_MERIDIAN_FRONTIER_AUDIT.md)
records the geometric and source qualifications.

**Repository evidence.** Its prime four-bridge augmentation-gap family
satisfies MRC and has double-cover group rank three, which already
excludes three original meridians. It cannot test the premise above.

**Single obstruction.** Three meridians give a specific deck-inverted
generating pair of the double-cover group. Realizing that marked inversion
by an equivariant genus-two Heegaard splitting is the missing geometric
step. Dropping the marking is false, as T(4,7) demonstrates.

**Impact test.** Proving this universal marked transfer settles the first
meridional-rank frontier; an actual knot counterexample disproves MRC.
Additional augmentation-gap families do neither.

## Gaussian safety-to-purity

**Quantifiers.** For every n>=1 and P in C[x1,...,xn], does
[for every Q there exists M_Q such that E(QP^m)=0 for all m>=M_Q]
imply [E(P^m)=0 for all m>=1], for independent standard real Gaussians?

**External frontier.** This reverses the usual Gaussian Moments Conjecture.
It must not inherit the status of that different implication.
[Long's 2026 preprint](https://arxiv.org/abs/2607.18186) disproves the usual
forward implication in dimensions at least three.
[Wilson's two-variable preprint](https://arxiv.org/abs/2607.23887)
addresses pure moments and one-sided support. A strongest published theorem
for the unrestricted reverse implication was not located; its independent
standing as a longstanding named problem has lower confidence.

**Repository evidence.** Extensive quadratic-fibre, invariant-core and
closure classes have safety-to-purity proofs. The unrestricted actual
original-marker image can be C[t] with a distinguished polynomial vacuum
f(t)=E(exp(tP)); the remaining assertion is f=1.

**Single obstruction.** Exclude a nonconstant polynomial vacuum in the
unweighted Gaussian model. Polynomial tilts do not preserve the desired
conclusion: the known positive tilted model has vacuum 1+t/2.
Thus module finiteness and positivity of tilted integrals are insufficient.

**Impact test.** A general vacuum-rigidity principle could be substantial.
More support or degree classifications alone fail the requested test.

## Strong Factorial Conjecture

**Quantifiers.** For every d>=1, nonzero F in C[x1,...,xd] with s nonzero
monomials, and every n>=1, at least one of
L(F^n),...,L(F^(n+s-1)) is nonzero, where L(x^alpha)=alpha!.

**External frontier.** Edo--van den Essen's
[original paper](https://arxiv.org/abs/1304.3956) links SFC to factorial
vanishing and Furter rigidity. Rocks, J. Algebra 478 (2017), 341--352,
proves powers of linear forms, sums of prime powers and constructions;
its [author preprint](https://people.math.sc.edu/rocks/BradyRocksSFC2015.pdf)
also distinguishes complete cases from partial binomial results.
The audit found no general support-independent consecutive-moment theorem.

**Repository evidence.** Ordinary binomials through degree six and selected
higher support families have internal proofs. The remaining degree-seven
residue families are narrower than the central objective.

**Single obstruction.** Exclude a common nonzero projective zero of s
consecutive homogeneous moment forms uniformly over arbitrary exponent
matrices and every starting n. Existing valuation arguments depend on
support geometry.

**Impact test.** General SFC would imply FC in every dimension and the stated
rigidity consequences. Extending a few terminal-index families would not
have that effect.

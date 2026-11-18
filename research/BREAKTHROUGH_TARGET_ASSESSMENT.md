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
comparison. The stronger [critical-genus obstruction](jacobian/mixed-characteristic-lifting/BGV_CRITICAL_GENUS_AND_RATIONAL_CORRECTIONS.md)
now excludes EVERY dominant rational source correction of the fixed integer
seed, even with an arbitrary dominant rational target correction. Its unique
critical curve has positive genus and ramification index two. This is an
all-degree field-geometric obstruction, not a bounded ansatz search.
The next useful action must change the seed's finite extension or use an
algebraic correction with a proved polynomial landing; another rational
normalization of that seed cannot resolve JC2.

**Impact test.** An actual plane pair or a universal impossibility theorem
would settle the remaining characteristic-zero conjecture. The radial
lifting theorem is useful route evidence, not achievement of that impact.
Further digits or optimal finite-lift degrees are not the next objective.

### Rejected transfer shortcuts: formal boundary data

A different JC2 proposal asks for a global polynomial Darboux pair on
`A_m=C[y,u,v]/(y^m v-u(u-1))`, for some `m>=2`. Its nowhere-vanishing
form is `Omega=dy wedge du/y^m`; the corresponding brackets are
`{y,u}=y^m`, `{y,v}=2u-1`, `{u,v}=m y^(m-1)v`.
The issue is a SINGLE global pair with `{F,G}=1`. Exactness of the form,
a sum of brackets equal to one, and solvability of all boundary jets
do not supply that pair. Such a pair would disprove JC2: the chart
`j(y,s)=(y,sy^m,s(sy^m-1))` pulls `Omega` back to `dy wedge ds`.
If the restricted Keller map were an automorphism, `s` would belong to
`j* A_m`. Its global extension would satisfy `y^m s=u`, impossible on
the omitted line `y=0,u=1`.

Here is an explicit internal countercheck against the boundary-jet
shortcut, uniform in `m`. Work in the y-adic completion of `A_m` and set

```text
q=2u-1,  B=1+4y^m v=q^2,
epsilon=q B^(-1/2),
R=sum_(k>=0) binom(-1/2,k) 4^k v^k y^(mk+1)/(mk+1),
r=epsilon R,
h=(2/3)(sqrt(1+3r)-1),
F=v^2-1-h,  G=v^3-v-(3/2)vh.
```

All square roots have constant term one. Since `epsilon^2=1`,
`d epsilon=0`, and `R_y=B^(-1/2)`, one has
`dr wedge dv=dy wedge dv/q=Omega`. Also
`r=h+3h^2/4`, and direct differentiation gives
`dF wedge dG=(1+3h/2)dh wedge dv=Omega`.
Thus `{F,G}=1` as a formal identity.

On BOTH boundary components `y=0,u=0` and `y=0,u=1`, the pair restricts
to `gamma(v)=(v^2-1,v^3-v)`. This is an immersion because its derivatives
`2v` and `3v^2-1` have no common zero, but
`gamma(1)=gamma(-1)=(0,0)`. For every `N>=1`, choose compatible global
representatives `F_N,G_N` modulo `y^(N+1)`. They exist by the definition
of completion, and satisfy `{F_N,G_N}=1 mod y^N`: the bracket lowers
y-adic order by at most one. Nevertheless this particular formal `F`
is not in `A_m`. On the closed line `v=0,u=0` it becomes
`-1-(2/3)(sqrt(1-3y)-1)`, which is not a polynomial, since a square
polynomial cannot have the simple zero of `1-3y`.

For `m=2` the same check also tests a proposed invariance shortcut.
The conformally symplectic self-map
`H(y,u,v)=(y(2u-1),(2u-1)^2,4v)` has `H*Omega=4Omega`.
Its boundary action is `u -> 1, v -> 4v`. If `H*G` belonged to
`C[F,G]`, evaluation at the boundary points `v=1,-1` would give
the same value. Its actual values are `60,-60`. So even this exact
formal Darboux pair need not generate an H-stable subalgebra.
No assertion about automatic stability of GLOBAL polynomial pairs
follows from this formal countercheck.

The displayed all-order identities are a written argument; exact
truncation checks at `(m,N)=(2,3),(2,5),(3,4),(5,4)` are finite controls
only. They do not exclude another formal pair that terminates, another
global polynomial pair, or JC2 counterexamples. The missing step remains
a GLOBAL polynomial construction or obstruction. Further boundary jets,
or a stability claim based only on these local data, do not remove it.
No new full-conjecture result or independent review is claimed here.

A separate source check also prevents importing a general
rational-component inversion theorem from
[Assi's arXiv:1807.03991](https://arxiv.org/abs/1807.03991):
the author withdrew it in 2019 because the main argument was incomplete.
Any use of a valid rational-fibration inversion theorem must retain its
additional hypotheses.

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

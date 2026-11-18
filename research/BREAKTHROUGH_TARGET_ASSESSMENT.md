# Breakthrough target assessment

The objective remains a result with major external mathematical impact. None of
the five unrestricted targets below is resolved by the current repository.
Primary target: Meridional Rank Conjecture, three-meridian frontier. No secondary target this iteration.
The other candidates are audited here, not pursued as parallel family-expansion
programmes. Assessment date: 18 November 2026.

## Plane Jacobian conjecture

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

### Four-variable Hessian conjecture and planar degree bounds

**Tang's structural alternative.** For a complex four-variable polynomial
$f$ with nonzero constant Hessian determinant, Tang's Theorem 1.1 gives
either an invertible gradient or constant linear coordinates in which

$$
f=R(x,y)+zA(x,y)+wB(x,y),\qquad \det D(A,B)\in\mathbb C^\times.
$$

The repository's v1 reproduction records this alternative as HCTANG2, and
the Ricci-flat HC4 theorem uses it. The extracted v1 and v2 PDF texts agree
after normalizing the version stamps.
[Tang v1](https://hexagonmath.org/2611.00156v1),
[Tang v2](https://hexagonmath.org/2611.00156v2),
[HCTANG2](hessian/tang-hc4-audit/STRUCTURAL_REVIEW.md),
[Ricci-flat HC4](hessian/hc4-triangular-graphs/RICCI_FLAT_HC4_AND_ONE_DIMENSIONAL_JACOBIAN_IMAGES.md).

**Conditional corollary (HC4 through total degree 101).** Assume Tang's
structural alternative. An invertible linear change preserves total degree
and gradient invertibility, since $\nabla(f\circ L)=L^T(\nabla f)\circ L$.
For a hypothetical counterexample of degree at most 101, the alternative
therefore gives $\max(\deg A,\deg B)\le\deg f-1\le100$, since the
$zA$ and $wB$ monomials cannot cancel. The classical planar bound
implies that $G=(A,B)$ has a polynomial inverse. The last two gradient
outputs recover $b=(x,y)$ by $G^{-1}$. The first two recover
$v=(z,w)^T$ by solving

$$
D G(b)^T v=(f_x,f_y)^T-\nabla R(b).
$$

The matrix has a nonzero constant determinant, so its inverse has
polynomial entries. This contradicts the assumed counterexample.
[Tang v2, Lemma 1.2](https://hexagonmath.org/pdf/2611.00156v2),
[Moh's degree-100 result](https://doi.org/10.1515/crll.1983.340.140),
[Moh's statement of the bound](https://www.math.purdue.edu/~ttm/jacobian),
[Wang's revision of Moh's algorithm](https://doi.org/10.11650/twjm/1500407850).

The transfer takes any planar component-degree bound $d$ to potential
degree $d+1$. In particular, Nguyen's published degree-104 theorem
gives a conditional HC4 bound of 105 using Tang. The degree-108 lower
bound claimed in the Guccione–Guccione–Horruitiner–Valqui preprint would
give HC4 through 108, since it excludes planar counterexamples of degree
at most 107. In the reverse direction, cotangent doubling
$f=zA+wB$ gives a potential of degree exactly $1+\max(\deg A,\deg B)$
whose gradient is invertible exactly when $G$ is. Under Tang's
structural alternative, HC4 through degree $d+1$ and JC2 through
degree $d$ are equivalent over $\mathbb C$ for each $d\ge 1$.
The same reduction identifies JC2 as the remaining unrestricted
problem for HC4.

**Conditional degree-pair restriction.** Assume Tang's structural
alternative and Theorem 2.1 of the Guccione–Guccione–Horruitiner–Valqui
preprint. If a complex HC4 counterexample has total degree at most 125,
the base Keller map in Tang coordinates has maximum component degree
at most 124. Its invertibility would give an inverse to the gradient
by the displayed linear system. The planar degree-pair restriction
therefore forces the unordered degrees of $A,B$ to be $72,108$.
Since $zA$ and $wB$ have distinct transverse variables, the potential
has degree at least 109. The term $R(x,y)$ still allows a higher degree
within the stated range.
[Nguyen, Theorem 3.6](https://doi.org/10.2989/16073606.2025.2482655),
[Guccione–Guccione–Horruitiner–Valqui preprint](https://arxiv.org/abs/2204.14178).

**Giraudin's total-degree-five theorem.** The unrefereed preprint gives a
separate low-degree argument using weighted Hessian-determinant layers,
Schur complements, partial Legendre transforms, and exact Gröbner-basis
and radical-membership certificates. Its five ancillary verifiers report
247 checks, with the ternary top handled in a separate preprint. The
dead-direction lemma treats $f=P(u_1,u_2,w)+sR(u_1,u_2)$. If $P_{ww}\ne0$
as a polynomial, constant Hessian
determinant forces $R$ to be affine after a linear change in the base,
and the proof uses fibrewise HC2. If $P_{ww}=0$, write $P=P_0+wP_1$.
Then $\det\operatorname{Hess} f=\det D(R,P_1)^2$, and the proof uses
planar inversion through degree four. The repository already records
the preprint's scope. In the exact quintic Meng–Yang graph exclusion,
the graph has degree five and the pulled-back potential has total
degree sixteen, outside the preprint's hypothesis.
[Giraudin v1, Lemma 4 and Section 6](https://hexagonmath.org/pdf/2611.00161v1),
[ternary-top preprint](https://hexagonmath.org/2611.00134),
[graph-degree comparison](jacobian/affine-complements/NORMAL_FORM_COVERAGE_AND_COMMON_LEADING_DIVISORS.md#8-the-november-preprint-and-the-graph-degree-convention),
[quintic graph exclusion](jacobian/affine-complements/FULL_COMPLEX_QUINTIC_GRAPH_EXCLUSION_AND_CUBIC_NILPOTENCE.md).

**Four-body central configurations.** Tomar's mass map is proper and its
critical loci include folds and cusps. The interval certificates,
reported Lean formalization, and treatment of exceptional parameter
strata may inform certification at the boundary of compactified Keller
maps. A complex Keller map is étale on the affine source. If it were
proper, it would be a finite étale cover of a simply connected target
and hence an automorphism. For HC4, the potential application is to
certification within the existing work on nonproperness at infinity.
[Tomar v1](https://hexagonmath.org/2611.00378v1).

**HC4 research priorities**

1. Independently audit Tang v2, concentrating on the weighted-wall
   termination and the exceptional cases in Sections 3–4. Formalize
   the degree-bound transfer with the structural theorem as an explicit
   hypothesis.
2. Independently reproduce Giraudin's quintic proof. Replay the ancillary
   verifiers, including the separate ternary-top verifier, and check that
   the unary, binary and ternary normal forms cover every complex coefficient
   stratum, including the exceptional branches.
3. Keep global polynomial termination and nonproper planar Keller maps
   as the main breakthrough targets for this programme. The live
   implications are automatic first-velocity membership in the ordinary
   codifferent in a source-finite frame, automatic sixth-moment trace
   perfection, and elimination or realization of the remaining ramified
   polar folds with the full finite target normalization and trace pairing.
   [Codifferent criterion](hessian/hc4-triangular-graphs/CODIFFERENT_PREIMAGES_UNDER_TARGET_DERIVATIONS.md),
   [trace criterion](hessian/hc4-triangular-graphs/EVEN_MOMENT_LATTICES_AND_SIXTH_MOMENT_TRACE_TEST.md),
   [ramified polar folds](hessian/hc4-triangular-graphs/RAMIFIED_FOLD_VALUATIONS_AND_COMMON_COORDINATE_POLES.md).

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

## Restricted-prime sums

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

## Meridional rank: primary target

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

The concrete current test is T(4,5;2,2), the closure of sigma1^4
(sigma1 sigma2 sigma3)^5. Its (1,1)-presentation gives cover rank at most
two, and Himeno's bridge-index PREPRINT gives bridge number four.
The [frontier audit](meridional-rank/representations/THREE_MERIDIAN_FRONTIER_AUDIT.md)
records exact source scopes, its free-by-cyclic model and bounded searches.
Of 595 short-conjugator triples, 286 have finite-quotient nongeneration
certificates and 309 remain unresolved. No three-meridian generation
certificate or new meridional-rank bound has been found. The same audit
gives a trefoil example showing why orbifold generation of a specified
triple does not certify its original lifts.

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

A selected-saddle approach also fails the universal significance test.
Let Z,W be an independent Wick pair with E[Z^a W^b]=delta_ab a!,
and put R=P(X) ZW. For every original marker A(X)Z^a W^b,

```text
E[A Z^a W^b R^m]=delta_ab (m+a)! E[A P^m].
```

Expansion of an arbitrary marker into these monomials proves that P and R
have EXACTLY equivalent safety and purity. The factorials are nonzero in
characteristic zero; the converse uses a=b=0. Over C this stays within
ordinary real Gaussian variables with complex polynomial coefficients,
using Z=(U+iV)/sqrt(2), W=(U-iV)/sqrt(2).

The known necessary condition P(0)=0 for safety makes R vanish to order
at least three. Its selected critical branch g=t C grad R(g), g(0)=0,
is then identically zero and I-t C Hess R(g)=I. Thus an EQUIVALENT
unrestricted problem already has trivial selected critical branch, action
and quadratic response. Multiplication by (ZW)^k raises the vanishing
order arbitrarily with the same equivalence. Correct universal saddle
constraints cannot by themselves remove the nonconstant-vacuum gap.
The reviewed exploratory saddle packet is not promoted to a DAG result.

The vacuum can first be normalized WITHOUT weakening the unrestricted
question. If a safe nonpure P has last nonzero pure moment mu_L, then
R=-P^L/mu_L remains safe and has unit period 1-t. Thus it suffices to
exclude this linear-vacuum case; higher multiplicities are not essential.

Even Wick powers of an ARBITRARY Stein eigenmarker do not supply that
exclusion. Start with the actual
[unweighted seed R](moments/vanishing/POLYNOMIAL_SEEDS_AND_PURE_MARKER_DELAY.md)
with f_R=1-t. Add an independent null pair U,V and set S=(1-U)R, A=1+V.
Exact Wick contraction gives, for EVERY k>=0,

```text
Phi_S[exp(sV)]=f_R(t(1-s))=1-t+st,
Phi_S[A^k]=1+(k-1)t.
```

Since A is affine in a null variable, its Gaussian Wick powers equal its
ordinary powers. Phi_S[A]=1 is a genuine marked Stein eigenclass with
eigenvalue one, but the full Stein orbit of ALL these Wick powers stays
in span{1,t}. Other old markers still have poles, so S is unsafe.
This is an all-order counterexample to the proposed eigenmarker-power
construction, not a safe nonpure Gaussian example.

The same reparametrization preserves the complete actual image and safety:
for an old original marker B, the exact general identity is

```text
Phi_S[B U^a V^b]
  = b!/(b-a)! (-t)^(b-a) D_t^(b-a) Phi_R[B]   if b>=a,
  = 0                                        if b<a.
```

Both inclusions follow because the old image is a time-Weyl module and old
markers are retained. Consequently this construction also equips ANY
hypothetical safe linear-vacuum counterexample with such a bounded-Wick
eigenmarker. A useful global theorem must use further source correlations,
or a distinguished canonical marker, rather than powers of an arbitrary
supplied eigenmarker. Direct exact controls at m<=5,k<=8 confirm the
displayed moment identities and retain an old bad marker with moment m!;
the all-order conclusion uses the contraction formulas above.

The canonical-marker variant remains a separate, unresolved sufficient
lemma: does f=1-t force some Phi[P^(Wick k)] to be nonpolynomial?
The null wrapper does not refute it. For G(X,U) independent of V,
Gaussian heat commutes with evaluation U=0, so
Phi_G[G^(Wick k)]=Phi_(G|U=0)[(G|U=0)^(Wick k)] for every k.
Finite successful examples do not prove this universal assertion.

The same test survives radial padding, unlike the earlier saddle
constraints. Write R=ZW, H_R=exp(partial_Z partial_W), and for a>=1 put

```text
S_a(R)=H_R(R^a)
      =sum_(j=0)^a (a!)^2/[j!^2(a-j)!] R^j.
```

For Q=P R^a, heat factorization gives
Q^(Wick k)=P^(Wick k) H_R^-1(S_a^k). The corresponding m-th correlation
is the old canonical correlation multiplied by

```text
(am)! sum_j [R^j]S_a^k (am)_j,
```

where (n)_j is falling factorial, zero for j>n. This multiplier is strictly
positive for EVERY integer m>=0: all coefficients are nonnegative and the
constant term is (a!)^k. Hence polynomiality of each fixed canonical
response is preserved and reflected. Rescaling Q by 1/a! keeps its vacuum
1-t when the old vacuum is 1-t. The test therefore retains its full
unrestricted content even on phases with arbitrarily high vanishing order.
Exact controls at a<=4,k<=5,m<=6 checked 168 multipliers; the formula
proves the unbounded assertion.

For standard covariance, the classical imaginary-Gaussian formula
([Hu--Yan, Theorem 5.3](https://arxiv.org/html/0901.4911v1)) gives, with
H=exp(Delta/2) and an independent Gaussian copy Y,

```text
Phi_P[P^(Wick k)]
  =E[(HP)(X+iY)^k exp(tP(X))].
```

All expressions are finite polynomial Wick operations before extracting
the time coefficients; no analytic Wick exponential is needed.

A tempting closure shortcut fails. Even when P is safe, one cannot infer
safety of S=P(X)-(HP)(X+iY) from polynomiality of the fixed-k responses.
For an explicit safe input use a Witt pair Z,W and an ordinary Gaussian T,
put q=2ZW+T^2 and P=W-qT-q^2 Z/2. It is the image of W under the
Gaussian-preserving polynomial automorphism exp(q(Z partial_T-T partial_W)).
Direct heat calculation gives

```text
HP=W-(q+5)T-(q^2+14q+35)Z/2.
```

For the coupled phase S, the squared Gaussian gradient at the origin is
35, so its selected critical action starts 35t^2/2. The retained
[critical-action obstruction](moments/vanishing/GAUSSIAN_G_FUNCTIONS_AND_GLOBAL_REGULAR_CORE.md)
therefore makes its UNIT period nonpolynomial. This proves unsafety
without assuming the target reverse implication. Fixed-k polynomiality
cannot be promoted to polynomiality of the coupled diagonal series.
The special hypothesis f_P=1-t, rather than that invalid closure, must
supply any proof of the canonical detection claim.

The ordinary marked Cauchy transform C_B=sum_m E[B P^m]t^m has the same
polynomiality condition as the exponential transform. With
B_0(h)=(h-h(0))/t and L_0=t^2 D_t+t, one has [B_0,L_0]=1 and the exact
unweighted source relation C_((X_i-partial_i)B)=L_0 C_(P_i B) in standard
covariance. The abstract Weyl relation alone also holds for the positive
tilted model with nonpure vacuum. It does not replace these full spatial
relations or prove that the distinguished unit is primitive.

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

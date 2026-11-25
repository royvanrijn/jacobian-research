# Breakthrough target assessment

The objective remains a result with major external mathematical impact.
The restricted-prime binary theorem is established internally under its
stated analytic imports. The assessment below identifies the remaining
targets, including shorter windows and dependent correlations.
Primary target: the unrestricted Strong Factorial Conjecture. MRC is deferred after the current fixed-input checkpoint, at the user's request. One secondary check: unrestricted Gaussian safety-to-purity through the canonical Wick-marker detection obligation described below.
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

**Degree restriction and the residual Laurent systems.** Use the
Guccione–Guccione–Horruitiner–Valqui reduction in its existential form.
If a counterexample of maximum component degree below 125 exists, a
reduced candidate has degrees 72 and 108 and reaches the remaining
\((8,28)\) case of their Section 2. An unordered degree pair in the
reduction describes a chosen representative. Indeed, composing an
unequal-degree pair \((P,Q)\) with the target shear giving \((P+Q,Q)\)
preserves the Keller property, invertibility and maximum degree while
making both component degrees equal to \(\deg Q\) when \(\deg Q>\deg P\).

Assume Tang's structural alternative. A four-variable Hessian
counterexample of potential degree at most 125 gives a planar
counterexample of maximum component degree at most 124. Under the
planar reduction, existence in this range therefore supplies a planar
candidate of degrees 72 and 108. Conversely a planar counterexample
with these degrees gives, by cotangent doubling, a Hessian
counterexample of degree exactly 109, without a structural assumption.
Under both inputs, the existence question through potential degree 125
is equivalent to existence at degree 109, and the bounds through
potential degree 108 and planar degree 107 follow.
[Guccione–Guccione–Horruitiner–Valqui, arXiv:2204.14178v1, Theorem 2.1 and Section 2](https://arxiv.org/abs/2204.14178).

Proposition 4.3 sends the remaining case to two explicit Newton
configurations in \(\mathbb C[x,x^{-1},y]\) with bracket \(x^2\).
Helali reports characteristic-zero unit-ideal identities excluding the
two transcribed coefficient systems. In the harder configuration a
sign split and an involution reduce the exclusion to a forced
specialization followed by exact identities. Strinz gives a different
valuation argument with a separate hypothesis on its alpha-strip
normalization. If the GGHV reduction is exhaustive and the Helali
transcription and exclusions hold for the actual Laurent configurations,
then planar inversion holds through degree 124. Under Tang, the same
inputs give Hessian inversion through potential degree 125.
Independently, any proof of HC4 through potential degree 125 gives
planar inversion through degree 124 by doubling. The coefficient
identities and their normalization have an explicit verification
obligation in the graph.
[Helali's report and exact certificates](https://github.com/bilLkarkariy/jc2-72-108-exact-certificates)
and [Strinz's conditional valuation proof](https://github.com/wstrinz/plane-jacobian-72-108).

[Nguyen, Theorem 3.6](https://doi.org/10.2989/16073606.2025.2482655)
gives the published planar bound through degree 104.

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
For every sufficiently large N=2 mod 12 the repository proves N=a+b with
a,b in C, including the fixed interval N/3<=a,b<=2N/3. The short-window
target asks for every sufficiently large admissible N to have such a pair with
|a-N/2|<=N^(1-epsilon) for some fixed epsilon>0. A separate target asks
for pointwise correlations with three dependent factors.

**External frontier.** The directly inspected published inputs include
Loughran--Matthiesen's JEMS 2024 Frobenian correlation theorem, BGLM's GAFA
2023 signed dispersion estimate, Iwaniec's half-dimensional sieve and
Henriot's corrected upper bound. They give substantial averaged or
one-sided tools. The [full Iwaniec--Pomykała comparison](arithmetic/SIMULTANEOUS_GAUSSIAN_EISENSTEIN_SUMS.md#full-iwaniec--pomykała-theorem-audit)
is now complete. Theorem 1 proves unconditional quartic differences.
Theorem 2 already gives the exact \(C_{12}\) sum problem, with coprime
summands and the expected logarithmic order, conditionally on its
growing-shift convolution estimate. It supplies no leading constant.
The present unrestricted leading asymptotic must be compared against
that direct antecedent; later deconditioning, balanced adaptations
and historical priority have not been exhausted.

**Repository evidence.** The [subgroup correlation theorem](arithmetic/SPLITTING_SEMIGROUP_LOCAL_GLOBAL.md)
gives the unweighted pointwise leading asymptotic for every sufficiently
large admissible target, including fixed macroscopic intervals. The proof
uses the retained ZF and G imports applied in
[PointwiseBinaryCorrelations](arithmetic/PRIMITIVE_NORM_BINARY_CORRELATION.md#6-discharging-the-imported-inputs).
The internal proof status carries the review limits stated in those sources.

**Alternative switching route.** The [critical switching criterion](arithmetic/ASYMMETRIC_CRITICAL_SWITCHING.md)
still needs its signed or absolute remainder estimate. Its companion
divisor sum produces weighted combined moduli near N^(3/4), outside the
applicable BGLM and Fouvry--Tenenbaum estimates. A finite integer-character
convolution misses the half-integer exponents of tau*1_C. Fractional
coefficients retain medium bilinear boxes and lose the existing periodic
fallback. The retained local checks cover the finite identities used by
this route.

**Impact test.** A pointwise lower bound in a power-short window would
extend the established binary theorem. A pointwise dependent three-factor
asymptotic would extend the correlation theory to a different configuration.

## Meridional rank: deferred after the current checkpoint

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

The former test input T(4,5;2,2), the closure of sigma1^4
(sigma1 sigma2 sigma3)^5, is now EXCLUDED as a counterexample:
a direct F_13 braid-colouring certificate proves mr(K)=b(K)=4.
The [frontier audit](meridional-rank/representations/THREE_MERIDIAN_FRONTIER_AUDIT.md)
gives the complete trace argument and four small matrix colours.
They satisfy all nineteen literal braid crossings. An even word has
two determinant-one blocks with traces 2 and 4, exchanged by a meridian.
Any deck-inverted generating pair would preserve every two-dimensional
trace over every field, contradicting this witness. This excludes ALL
meridian triples, including the 309 survivors of the earlier finite
panel. The four-braid supplies the upper bound, so this particular
equality needs neither Himeno's preprint nor a software-derived group.

The closed double cover separately passes interval hyperbolicity after
permanent filling and retriangulation. Its ordinary rank and Heegaard
genus are exactly two. The trace witness excludes cyclicity, and the
(1,1)-presentation supplies the upper bounds. The actual deck action
admits no inverted generating pair, and the least genus of a splitting
with trivial-tangle quotients is three. The same obstruction gives
meridional generating number four in every even-order orbifold quotient.
The universal three-meridian implication remains open.


The further fixed input T(5,19;2,1), up to mirror, is also excluded
from three-meridian generation. Its cover has rank and Heegaard genus
two and its actual deck action fixes EVERY two-dimensional word trace
over every field. Nevertheless a direct seven-sheet branched-cover
certificate proves mr(K)>=4: the cover group has a nontrivial F13
matrix image, whereas three meridians would make it trivial by the
seven-vertex, six-edge Schreier graph argument. All 60 raw relations
are checked from the literal 78-crossing braid, with no software-derived
group or surgery marking in that lower bound. The same obstruction
excludes EVERY inverted generating pair in the actual deck outer class.
Thus trace invariance and ordinary rank two are insufficient.

Vafaee and Krishna--Morton give b(K)=5 using published theorems.
Brittenham--Wu's published classification gives hyperbolicity of both
nonintegral surgeries used to identify the double cover. The exact
remaining distinction is mr(K)=4 versus5. Only an original four-meridian
generation certificate would give a counterexample; none is certified.
The frontier audit records the full small certificate. This is a
falsification of another three-meridian input, not a breakthrough or a
new positive family theorem.

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

Two more shortcuts fail before the vacuum gap. With H=exp(Delta/2),
safety is exactly ord_0 H(P^m)->infinity: Hermite markers recover every
Taylor coefficient. The transported product
A*B=H((H^-1 A)(H^-1 B)) is NOT continuous in this order topology:
for a Wick pair Z,W, both Z^m and W^m tend to zero in order, whereas
their transported product has constant coefficient m!.
Consequently topological nilpotence does not justify treating evaluation
at zero as a multiplicative residue map. Nor does H preserve safety.
For the safe nonlinear phase P=W-qT-q^2 Z/2 above, exact contraction gives
E[(HP)^2]=-30 and the squared gradient of HP at zero is -10.
Its selected action begins -5t^2, so the retained critical-action theorem
makes HP unsafe. Applying the known safe-phase origin condition to HP
would therefore assume a false closure property.

A nonconstant polynomial tilt also cannot be removed by any stable
invertible polynomial coordinate change preserving the COMPLETE moment
functional. Here is the intrinsic obstruction, valid after adjoining
independent Gaussian variables. For D=sum_i A_i partial_i, preservation
of E[rho Q] for every Q is equivalent, by Gaussian integration by parts
and nondegeneracy, to

```text
D(rho)=rho (X.A-div A).
```

Every such D is tangent to each irreducible divisor of rho; hence their
evaluation span drops rank at a smooth point of any nonconstant factor.
For rho=1 in dimension at least two the span is full at EVERY point p:
choose j!=i, R=X_j-p_j, A_i=partial_j R-X_j R and
A_j=X_i R-partial_i R, with other components zero. Then div A=X.A
and D(p)=e_i. A polynomial automorphism conjugating the complete
functionals must identify these rank loci, which is impossible.
Dimension one reduces to dimension two by adjoining a fixed Gaussian.
A constant tilt must be one by evaluation on Q=1. This excludes
coordinate untwisting of the positive tilted model, while leaving
phase-specific suspensions and noninvertible maps outside the argument.

Even the normalized vacuum 1-t does not permit an exact reduction of
the FULL embedded marked image to cubic phases. Independently add its
known polynomial seed to a pure power phase Q_(0,k) from the
[seed construction](moments/vanishing/POLYNOMIAL_SEEDS_AND_PURE_MARKER_DELAY.md).
The unit stays 1-t, while a fixed marker has period (1-t)h_k(t), where
h_k=sum_m (km)!t^m/m!. For m>=1 its m-th raw moment is
(km)!-m(k(m-1))!, asymptotic to (km)! for k>=2.
For every cubic phase in any fixed finite Gaussian dimension and every
fixed polynomial marker A, absolute moments are at most
C^m Gamma((3m+deg A+dimension)/2), up to a fixed factor.
Polynomial time coefficients and finitely many time derivatives preserve
this growth order.
The factorial growth orders differ when k>=2. Thus these normalized
inputs cannot share a cubic phase's full embedded image, even after
finite Gaussian stabilization. These witnesses are unsafe, so this
does not exclude a reduction restricted to safe inputs, a time change,
or an implication that retains only safety and purity.

**Impact test.** A general vacuum-rigidity principle could be substantial.
More support or degree classifications alone fail the requested test.

## Strong Factorial Conjecture: primary target

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
higher support families have internal proofs. The
[ordinary constant-plus-power theorem](factorial/image-rigidity/ALL_UNIVARIATE_TWO_TERM_SFC.tex)
already covers A+BX^d for every positive integer d and every starting index.
Its uniform contour proof and finite certificates do not treat arbitrary
monomial products or supports. The remaining degree-seven residue families
are narrower than the central objective. Full Furter
rigidity already has a credited internal reproduction of liqsweep's
critical-value argument; its factorial transfer is a specific all-dimensional
family, rather than arbitrary supports.

The [unrestricted-mechanism audit](factorial/image-rigidity/UNRESTRICTED_SFC_MECHANISMS.tex)
proves two route obstructions: even a two-term phase can require scalar
homogeneous recurrence order three, and every polynomial automorphism
preserving the entire factorial functional is a variable permutation.
The sharp root-of-unity linear examples also defeat a universal
unique-minimum valuation argument. These are method checks, not an SFC
solution. The exact Gamma-factorization route has two further limits:
arbitrary positive Euler factors can create forbidden zero blocks, while
the actual finite Gamma/Beta blocks lose moment-Gram coercivity at toric
relations in the ordinary limit. The audit gives the exact operator
counterexample and the degeneration for support {1,X,X^2}. Neither is an
SFC counterexample. Finite-block separation and direct control of its
limiting nonvanishing are both missing. Gamma-product Appell simplicity
remains an unproved narrower possibility, with no arbitrary-support
transfer.

Two further direct transfers now have exact counterexamples in the audit.
For F=c0+c1X^2 at orders two and three, no prime and diagonal valuations
produce a regular pair of monomials or rescaled complete homogeneous
forms, despite the original pair being regular. Also
(XYZ)^7(X^5+Y^5+Z^5) has four strict maxima on the positive simplex,
exceeding its three monomials. Thus neither those diagonal targets nor
an ungrouped saddle count supplies the general support-size window.

The proposed raw Gamma-resultant coefficient certificate is now
**falsified**. For the primitive row A={0,1,3,7}, n=1, the exact
resultant has degree 238, order 72 at zero, 162 positive coefficients,
and five negative coefficients at powers 93 through 97. Thus neither
orientation has nonnegative coefficients. Integer-only fraction-free
elimination and a separate polynomial-matrix determinant calculation
agree on all 239 coefficients. The replay is retained in the existing
mechanism checker as --gamma-resultant-counterexample.

This is a counterexample to the proposed coefficient formula, not SFC:
the ordinary specialization R(1) is positive. The earlier {0,1,3},
window 2,3,4 calculation remains a positive calibration only. Connected
row-permutation diagrams still give nonnegative joint cumulants,
but the new calculation shows that their positivity does not survive
resultant elimination.

The row-coalescence identity and primitive positive-product rigidity
lemma remain valid. They constrain that proposed method; they cannot
repair the false coefficient assertion, whose counterexample row is
already primitive. A general SFC argument must prove strict nonvanishing
by a different representation or mechanism. No replacement general
certificate is currently established.

Two different general transfers have now been ruled out. Even for
the fixed support {1,X}, no single nonzero scalar normalization sequence
can turn the moment series of both 1+X and 2+X into algebraic germs.
Thus a direct transfer to the credited polynomial-inverse rigidity theorem
cannot use a common coefficient-independent moment normalization.
Separately, an exact positive-measure example has three regular consecutive
moment forms that cease to be regular after a factorial Borel transform.
Both the original and transformed measures have strict polynomial Gram
positivity. These are proof-route obstructions, not SFC counterexamples.

Even Fock orthogonality and homogeneous holomorphic features do not
replace the actual monomials. An exact three-feature example with
nonzero coefficients has its first five Gaussian moments zero and
sixth moment -259200. Its features mix occupation-number monomials
and have angular cross terms. This falsifies that broader transfer,
while leaving ordinary SFC open.


Positive real Mellin multiplicities do not provide a replacement transfer.
The existing mechanism audit constructs positive variables with moments
(k!)^tau, retaining the canonical log-Gamma negative-jump kernel. There is
a unique tau* in (1/4,1/3) for which two independent copies, together with
a constant feature, have their first three phase moments zero. Increasing
only one fractional multiplicity from 1/4 to tau* changes a regular frame
to that nonregular frame. Thus continuous Mellin positivity and strict
polynomial Gram positivity do not preserve the required window.
The exceptional multiplicity is below one; no actual integer-row inflation
or interpolation already anchored at multiplicity at least one is refuted.
The next transfer argument must retain a condition that excludes this
fractional example.


**Single obstruction.** Exclude a common nonzero projective zero of s
consecutive homogeneous moment forms uniformly over arbitrary exponent
matrices and every starting n. A concrete discrete transfer scheme would
suffice: preserve regularity when adjoining a row e_j to a nonconstant
column, and when identifying two exponent rows by their sum while
columns stay distinct. Independent affine-linear supports provide the
base; token inflation and row identification produce every ordinary
support. Both preservation statements remain unproved for actual factorial
moment frames. The positive-measure counterexample inflates a constant
feature and therefore does not refute the restricted first statement.
The missing step must retain joint integer-toric moment structure:
Stein conjugation of differential operators does not yet yield exactness
of the multiplication Koszul complex after evaluation on the vacuum.
A faithful row-identification interpolation is available through squared
moduli of two correlated complex Gaussians: its moments are
p!q! sum_j binom(p,j)binom(q,j)rho^j, with independent and identified
rows at rho=0 and rho=1. It retains the integer Wick couplings but
does not yet supply strict nonvanishing of the consecutive resultant.

The stronger complex-correlation-disk certificate is false on the actual
integer source features (Y_rho,X Y_rho): their first two forms have
resultant 2(1+4rho-rho^2), vanishing at 2-sqrt(5) in (-1,0).
The merged features X,X^2 remain distinct, and the physical endpoints
have resultants 2 and 8. Thus complex analytic zero-freeness cannot supply
the missing physical row-merge theorem; real rho in [0,1] remains open.

An all-moment alternative was tested through polynomial infinitesimal
symmetries of a simplex. It is now excluded in every degree for planar
polygons: every polynomial derivation preserving the entire uniform
polygon integral has no nonzero constant eigenvalue on rational
functions, even with complex coefficients. More generally, every rational
logarithmic derivative Dr/r that is polynomial has integral zero.
Its polynomial logarithmic primitive vanishes on each edge; Stokes
therefore excludes a nonzero constant cofactor. This closes the entire constant-cofactor Hamiltonian construction
on the triangle, not the general triangle moment problem.
Positive-moment vanishing alone does not provide such an eigenfunction.

**Impact test.** General SFC would imply FC in every dimension and supply
a support-independent consecutive-window principle. The existing rigidity
proof does not supply that principle. Further terminal-index families,
bounded discriminant positivity or a reproof of rigidity would not achieve
the requested breakthrough.

# External mathematical transfers

The external result catalogue is pinned to
[solved-conjectures main, 758066055ef4](https://github.com/1337isnot1337/solved-conjectures/blob/758066055ef49c2c66cb1087918ae62ad75e7ac1/CONTENTS.md).
The research notes are pinned to
[research, dab1ac8d1ece](https://github.com/1337isnot1337/solved-conjectures/tree/dab1ac8d1ecea8dc34c38932fcd2a4b3f56ce0ef/research).
The repository's claim IDs refer to [the current mathematical DAG](MATH_STATUS.json).

| External result | Mathematical transfer and exact DAG target |
| --- | --- |
| Hajdu–Sárközy shifted powers and shifted products, Erdős perturbed squares, and Sárközy smooth fixed gaps | The [determinant and correlation audit](arithmetic/SHIFTED_PRODUCT_AND_CORRELATION_TRANSFER_AUDIT.md) reconstructs the four-cycle exclusion, the summable-coefficient and no-log argument, and the mechanism for eliminating the unknown third-order term in a fixed-gap test. The square-full determinant forces the singular case at large prime targets. The six-cycle and pointwise kernel-pair nodes still require independent bounds. For the Gaussian-norm frontier \(\texttt{GAP-ARITH-CK-AVERAGED-CENTERED-TRIPLES}\), pair limits plus the fourth-moment upper sieve yield triple occurrence for windows longer than \((3+\sqrt3)/6\) mean spacings, with the retained ZF and G inputs. The first moment alone, with the same upper sieve, suffices unconditionally beyond two mean spacings. |
| Karp–Vishnyakova–Zhang, Nuttall, Csordas–Vishnyakova and Cardon et al. | The [factorial audit](factorial/image-rigidity/GAMMA_QUOTIENT_TRANSFER_AUDIT.md) independently derives the all-order gamma-quotient signs and a product determinant on unit-spaced columns. It records the scope of the external Lean receipt. The actual unweighted linked beta array has opposite \(2\times2\) minor signs. An exact contiguous subtraction proves all-order consecutive-moment coprimality for the axis-supported monomial-linear family in \(\texttt{GAP-FACTORIAL-BINARY-LINKED-BETA-COPRIMALITY}\). Distinct positive exponents and the uniform Hurwitz determinant obligations remain open. The zero-geometry comparisons include a negative cubic test for the actual gamma multipliers. |
| Sokal Apéry threshold, Hamburger classification and Ferreira–Simon beta-prime convolution | The [Gaussian audit](moments/vanishing/GAUSSIAN_W_AFFINE_ROTATIONAL_DIVERGENCE.md) collapses a rotation orbit of dominant singularities to one actual-kernel singularity. This covers a symmetry slice beyond \(\texttt{GAUSSIAN-W-AFFINE-SELECTED-DOMINANT-DIVERGENCE}\). It also gives an exact positive-measure classification for normalized coefficients of \(W(1+Z^m)Y\), with separate real-line and nonnegative supports. General nonsymmetric dominant-point cancellation remains open. |
| Hecke bootstrap, CM descent, Artin uniform and average estimates, and Goldbach spectral bridge | The [analytic audit](arithmetic/EXTERNAL_HECKE_ARTIN_GOLDBACH_TRANSFER_AUDIT.md) extends the conditional window interface to the actual angular characters and norm twists, checks the Artin weighted family, and gives an unconditional CM square-norm inclusion and fixed-field splitting density. The gcd-stratified Goldbach projection transfers exactly. The operator-norm interface for shrinking windows loses two powers of the relative width and still requires an estimate for the actual weights. The Eisenstein central-window, sharp four-form energy and switching nodes retain their analytic obligations. |
| Ma–Tang finite-group random subset sums | [Local intensities](arithmetic/RANDOM_PRIME_SEMIGROUP_LOCAL_INTENSITIES.md) proves an unconditional strict improvement of the critical gap-density lower bound and an explicit parity obstruction to any residue-blind subcritical pointwise limit. The exact Walsh certificate in the pinned source uses a binary affine fibre. In the semigroup model the remaining conditions for a critical density are an empirical intensity law and a joint-support bound. |
| Féray conjugate-power words | The [Hurwitz census](arithmetic/HURWITZ_EXACT_CENSUS_AND_REMOVABLE_FACTORS.md) uses conditional class Fourier calculations to give exact absolute totals and to count directly removable last factors, with limiting fraction \(e^{-2}\). The few-cycle mixing hypothesis fails strongly for double transpositions, whose two-factor product stays at total variation \(1-o(1)\) from the parity-uniform law. The exact census targets \(\texttt{GAP-ARITH-HURWITZ-UNIFORM-REMOVABLE-DT-REPRESENTATIVE}\) and the uniform fixed-product fibre connectivity nodes. |
| Stone nonproper-locus obstruction | The degree criterion and the polynomial comparison below locate the missing hypothesis for \(\texttt{GOAL-JC2}\) and \(\texttt{GAP-JACOBIAN-1}\). A known three-dimensional Keller map gives an exact obstruction to inferring generic uniqueness from even a full target line of singleton fibres. |

## Further consequences and remaining estimates

- In the pure and marked mixed Hurwitz domains, the removable-factor
  count has an exact finite probability law and a normal limit. The
  exceptional fraction decays faster than every fixed inverse power.
  For each fixed \(C>0\), any orbit with relative size at least \(r^{-C}\)
  therefore meets the removable locus for all sufficiently large \(r\).
  The unresolved step for \(\texttt{GAP-ARITH-HURWITZ-UNIFORM-REMOVABLE-DT-REPRESENTATIVE}\)
  concerns orbits smaller than these thresholds. [Joint law and orbit criterion](arithmetic/HURWITZ_EXACT_CENSUS_AND_REMOVABLE_FACTORS.md#the-joint-removable-factor-law)
- In native increasing row and column order, a positive-exponent
  raw linked beta matrix is weakly sign regular exactly when \(u/r\)
  is a positive integer. The arbitrary-minor zero criterion gives
  an explicit test for which minors survive at each order.
  The physical moment polynomials also contain binomial row
  factors, and the remaining linked-beta coprimality problem
  requires control of their common complex roots. On the axis,
  the shifted-shape theorem covers \(\Gamma(1)\times\Gamma(\beta)\)
  for every \(\beta>0\).
  [Sign classification and shifted-shape proof](factorial/image-rigidity/GAMMA_QUOTIENT_TRANSFER_AUDIT.md)
- The finite-modulus critical gap bounds have a sharp truncation
  asymptotic along \(\operatorname{lcm}(1,\ldots,\lfloor y\rfloor)\),
  and the all-powers-of-two bound is explicit. The remaining
  critical-density route requires an empirical conditional intensity
  law and control of overlapping unexposed prime supports. The
  quantitative truncation concerns the proposed law itself.
  [Truncation and density criterion](arithmetic/RANDOM_PRIME_SEMIGROUP_LOCAL_INTENSITIES.md#strict-two-adic-improvement-and-sharp-finite-modulus-error)
- Gaussian divergence survives when the inverse equation rotates
  and the numerator breaks rotation symmetry. The full numerator
  separates into the \(m\) residue classes of exponents and each
  component has the same quotient singularity majorant.
  The next case is an inverse with tied dominant points and no such
  rotation of its defining equation.
  [Actual source proof](moments/vanishing/GAUSSIAN_W_AFFINE_ROTATIONAL_DIVERGENCE.md)
- Optimizing the centered Goldbach Hardy test shape has sharp
  relative norm cost \(e^2\pi/(32\Delta^2)\).
  This quantifies the operator-norm route. For the unresolved
  short-window estimate, a direct bound on the moving quadratic
  form could avoid the full norm cost.
  [Sharp cost and moving-window version](arithmetic/EXTERNAL_HECKE_ARTIN_GOLDBACH_TRANSFER_AUDIT.md#sharp-hardy-cost-over-all-supported-test-shapes)

## Transfers from individual proof steps

- Browning–Verzobio combine a determinant cover, a resultant, a proper-power
  test on the normalization, and a degree-sensitive projection. Bounding each
  valuation by the degree of the zero or pole divisor improves their
  coefficient-uniform rectangular count. In the fixed-offset surface argument
  this lowers the leading finite-size loss from \(3r^{-3/2}\) to
  \((2+O(r^{-1}))r^{-3/2}\).
  [Divisor-degree proof and box transfer](arithmetic/POWERFUL_SUMS_SLICE_POWER_TEST.md)
- The generalized Vandermonde in the powerfulness manuscript admits a primewise
  estimate from counting pairs of grid weights in each residue class. This gives
  the smaller all-fixed-offset error \(O_{u,v}(r^{-17/8}(\log r)^{1/8})\).
  The same square grid has a degree-budget barrier toward the target \(r^{-5/2}\).
  [Primewise determinant proof](arithmetic/POWERFUL_MONOMIAL_VALUATION_REFINEMENT.md)
- Mao–Wang–Zheng's symmetry argument fails for unequal linked beta profiles.
  Direct beta moment identities instead give the expected degrees and positive
  leading coefficients for both initial Euclidean remainders. The fixed-shape
  Jacobi ODE gives an individual discriminant product. The varying-parameter
  Jacobi formula applies to the actual consecutive pair, where shared-branch
  comparison and endpoint separation remain to be proved.
  [Exact pivots, source hypotheses and Jacobi receiver](factorial/image-rigidity/LINKED_BETA_UNEQUAL_PROFILE_PIVOTS.md)
- Ordered beta grids and gamma-product density estimates transfer to nonlinear
  Gaussian monomial fibres. The original Wick markers recover the full
  differential rank and the polynomial-time lattice. The raw coefficients
  have strictly positive Hankel minors, even in the range of moment
  indeterminacy. The only globally regular actual submodule in this family
  is the constant line. The inverse kernel also gives unit divergence for
  arbitrary nonconstant spatial coefficients on these monomial fibres,
  using the decay of every inverse branch at infinity.
  [Exact rank, convergence and moment law](moments/vanishing/GAUSSIAN_NONLINEAR_MONOMIAL_FIBRE_FULL_IMAGE.md)

## Further consequences of the transferred proofs

- The determinant depends on the relative Taylor change of the Vandermonde
  quotient. Cauchy–Binet bounds that change using the largest weight, even
  in examples where the valuation of the diagonal quotient is of order
  \(h^3\log r\). The short-congruence loss becomes
  \(O_A((\log r/r)^{1/5})\). Splitting targets by their large-prime part
  improves the unrestricted fixed-offset error to \(O_{u,v}(r^{-11/5})\).
  [Taylor expansion and target split](arithmetic/POWERFUL_MONOMIAL_TAYLOR_REFINEMENT.md)
- Poles at infinity and equal-exponent fibres constrain the proper-power
  degree in the surface projection. Nearby distinct powers have a
  reducible-slice exponent of order \(1/r\), and the equally powered case
  needs only \(W^2+W M^{1/r}\). Optimizing the actual side lengths also
  gives an explicit \(r^{-5/2}\) error coefficient for offsets \((u,0)\).
  [Pole degrees and fixed-offset constants](arithmetic/POWERFUL_SUMS_INFINITY_POLE_BOUND.md)
- Fraction-free division proves positive third and fourth Euclidean
  remainders at every eligible linked-beta order. The exact checks extend
  full chains through starting order nine and consecutive coprimality
  through starting order ten, for every nonnegative integer profile.
  [Later pivots and exact replay](factorial/image-rigidity/LINKED_BETA_LATER_PIVOTS.md)
- Parity normalization and circular residues extend the full convergent
  Gaussian core to every polynomial \(c(Z)\) in \(Wc(Z)\prod_iY_i^{k_i}\).
  A nonzero scalar shift eliminates the globally regular core. The exact
  diagonal calculation covers total active exponent at most two and,
  for constant \(c\) or \(Z^2\mid c\), arbitrary diagonal coefficients.
  [Full-core proof and diagonal descent](moments/vanishing/GAUSSIAN_MONOMIAL_FIBRE_CONVERGENT_CORE.md)

## Incremental upstream audit

The main catalogue is still at its pinned commit. The research branch
has two additions since the prior audit, the Erdős–Turán finite-prefix
complexity argument and the Vizing aggregate-bound construction.
[Pinned commit comparison](https://github.com/1337isnot1337/solved-conjectures/compare/dab1ac8d1ecea8dc34c38932fcd2a4b3f56ce0ef...d2cd4d3738570244ee1830c1d02236d93ecac11c)

- In the Erdős–Turán proof, alternative paths fix the same queried
  positions, and cut edges record all ambiguous stretches. The fixed
  read-set and the grammar projection transfer to the random prime
  semigroup. Independent prime labels give exponential membership
  complexity and an exact almost-sure information rate. The
  representation-count and critical-density problems require separate
  arithmetic estimates. [Prefix proof](arithmetic/RANDOM_PRIME_SEMIGROUP_PREFIX_COMPLEXITY.md)
- Aliabadi–Krop's Proposition 6.1 gives an aggregate minimizing
  segment for the \(2/3\) Vizing bound. The new branch constructs
  actual graphs at every rational point of the segment. Complement
  neighborhoods force a local witness in each block, and the
  Cartesian layer counts charge exceptional fibres to distinct
  centers. The same distinction between feasible aggregate data and
  stronger local information appears in the powerful-triples strict
  cardinality obstruction. There the abstract hypergraph meets the old
  projection bounds, while the arithmetic proof uses congruence classes.
  [Upstream construction](https://github.com/1337isnot1337/solved-conjectures/blob/d2cd4d3738570244ee1830c1d02236d93ecac11c/research/vizing-conjecture/README.md) · [Primary aggregate bound](https://arxiv.org/html/2607.01109v3) · [Powerful-triples obstruction](../unpublished/powerful-number-representations/main.tex)

Four older proof mechanisms also have new applications.

| Upstream method | Local consequence |
| --- | --- |
| Two-difference Motzkin signed dual orientations | The finite two-shift overlap region is determined by one character even with occupied-occupied edges allowed. This gives sharp joint-budget tests for shifted critical gap densities. [Character proof](arithmetic/TWO_SHIFT_OVERLAP_CHARACTER_BOUNDS.md) |
| Vishnyakova's falling-factorial tridiagonal recurrence | An exact quadratic pencil for the actual linked beta polynomials gives strict root annuli and a uniform forbidden sector on positive profiles. Any adjacent collision must satisfy two normalized Rayleigh equalities. [Actual pencil and bounds](factorial/image-rigidity/LINKED_BETA_VISHNYAKOVA_PENCIL.md) |
| Simon's log-Gamma Lévy comparison | Gaussian monomial-fibre coefficient powers remain strict Stieltjes sequences. Their exact raw determinacy threshold is \(\theta\rho=2\), and the leading small-power Hankel determinants are positive. [Mellin and determinant proof](moments/vanishing/GAUSSIAN_MONOMIAL_FIBRE_MOMENT_POWERS.md) |
| Sokal's common-vector operator test | A quadratic Gaussian angular deformation has sharp \(\sqrt{s}\) operator displacement from zero, with an exact arcsine upper realization. [Displacement proof](moments/vanishing/GAUSSIAN_ANGULAR_MOMENT_OPERATOR_DISPLACEMENT.md) |

## Bridges between programmes

| Source and receiver | What transfers |
| --- | --- |
| Random-prime intensities and prescribed-prime binary variance | The excluded-prime factor at order two equals the corresponding prescribed-prime singular-series factor. Its shifted correlation gives the full mixed local constant for two possibly different prescribed-prime semigroups. The known large power-window variance transfers this constant to actual binary counts, uniformly over sublinear power shifts. The short-width diagonal in the four-form energy problem remains open. In the reverse direction, classwise means give a sharp second-moment criterion for the random-prime intensity law. Its missing upper bound is a two-copy additive energy with shared exposed primes. [Local factors, shifted counts and energy criterion](arithmetic/SHIFTED_LOCAL_CORRELATIONS_AND_SEMIGROUP_ENERGY.md) |
| Binary Gamma coprimality and Gaussian vanishing | Blocks of dimensions two and \(d\) realize the two Gamma variables. An eligible marker in the \(d\)-block shifts the second shape and preserves the first. Consecutive zero moments force sharp coefficient-ideal nilpotence over a \(\mathbb Q\)-algebra. Angular separation then gives an explicit original-marker cutoff for a homogeneous family and its specified power of the Laplacian. [Two-block Gaussian and Laplacian proof](factorial/image-rigidity/GAMMA_QUOTIENT_TRANSFER_AUDIT.md#5-two-gaussian-blocks-fixed-markers-and-sharp-nilpotence) |
| Rotational Gaussian divergence and the regular actual-image problem | For \(P_m=W(1+Z^m)Y\), explicit original-marker tails give the full generic rank, a Laurent basis, and the exact polynomial-time lattice. The convergent generic core is rational and the maximal globally regular submodule is the constant line. This excludes this family from the transcendental regular realization search. [Exact actual-image calculation](moments/vanishing/GAUSSIAN_W_AFFINE_ROTATIONAL_DIVERGENCE.md) |
| Polynomial Hurwitz factorizations and the even alternating Laurent passport | The pure census supplies a polynomial acceptance bound for a tree with one added edge. The leaf-pair estimate transfers to the even unicyclic passport. In even degree \(n\ge12\) a removable factor peels to the full alternating passport in degree \(n-2\), and the allowable extensions use two distinct anchors on the lower long cycle. Small-orbit exclusion and equivalence of those anchor extensions are the remaining inputs to this route. [Even Laurent reduction](arithmetic/HURWITZ_EXACT_CENSUS_AND_REMOVABLE_FACTORS.md#the-even-alternating-laurent-passport) |

## Stone's criterion and the location of a singleton fibre

In the [pinned Stone note](https://github.com/1337isnot1337/solved-conjectures/blob/dab1ac8d1ecea8dc34c38932fcd2a4b3f56ce0ef/research/stone-nonproper-locus/README.md),
the complementarity map on each orthant is linear in nonnegative
coordinates with matrix \(B_I\). For a sequence with images converging
and norms escaping, take a fixed orthant along a subsequence and
normalize its nonnegative coordinates by their norms. The limiting
nonzero vector belongs to the nonnegative kernel of \(B_I\).
The limit of the images lies in the closed image cone of that
orthant. In the reverse direction one can add an arbitrary
nonnegative multiple of the kernel vector to any preimage in the
cone. This verifies the identification of the nonproper locus with
the union of the stated polyhedral cones and, in this setting,
with the locus of unbounded fibres.

The other necessary fact in the source is uniqueness at an interior
point of a nonsingular complementary cone. Pivoting on that basis
takes the strictly complementary solution to the solution with
zero complementary input and positive right-hand side. Another
solution would have a nonzero nonnegative input with strictly
negative matrix action on its support, forbidden by semimonotonicity
of the pivoted matrix. Q0 identifies the image with the convex
feasible cone. Along a compact path in its interior avoiding the
nonproper locus all preimages are uniformly bounded, so Brouwer
degree on one large ball is constant. At a unique preimage the
degree is the sign of the corresponding principal minor. At a
strictly positive right-hand side it is \(+1\). If the exceptional
cones intersect the interior only in codimension at least two,
the two-segment general-position argument supplies the path. These
are the two separate hypotheses that yield the degree obstruction.

For complex polynomial Keller maps, the complement of the **nonfinite**
locus carries the finite étale cover of generic degree \(d\). A
single fibre of size one over that complement forces \(d=1\).
The source being open in its finite normalization is essential
here. If the nonfinite locus has **complex** codimension at least two,
the full finite normalization has no height-one ramification. Purity
applies because the normalization is normal, the affine target is
regular, the finite map is quasi-finite, and local dimensions agree by
going-down. The result is a connected finite étale cover of complex
affine space and therefore has degree one. A Keller morphism of
degree one is a birational open immersion by the quasi-finite
form of Zariski's main theorem. Its image meets every target prime
divisor generically: for an irreducible polynomial \(q\) cutting out
the divisor, \(q(F)\) is nonconstant by dominance, and a prime
component of its zero set dominates the divisor by dimension and
quasi-finiteness. The open immersion image thus omits at most a
codimension-two closed set. The coordinate functions of its inverse
extend across that set by normality of the affine target. Both
polynomial composition identities follow on dense opens and hence
everywhere. [Affine inertia and cover proof](jacobian/monodromy-theory/THEORY.md#exact-purity-and-affine-cover-applications)

Real codimension two of a complex exceptional hypersurface can
provide path connectivity and constancy of \(d\), with no conclusion
about its value. Fibre cardinality at an exceptional point need
not give \(d\). For the retained three-dimensional Alpöge–Gallagher
map, put \(U=1+xy\), and write its three outputs as
\[
 A=U^3z+y^2U(4+3xy),\quad
 B=y+3xU^2z+3xy^2(4+3xy),\quad
 C=2x-3x^2y-x^3z .
\]
Its Jacobian determinant is the constant \(-2\).
At a point with \(A=B=0\) and \(U\ne0\), the identity
\(B=y+3xA/U\) gives \(y=0\), then \(z=0\). On \(U=0\) one has
\(B=-2y\ne0\). Thus the fibre of \((0,0,c)\) is the singleton
\((c/2,0,0)\) for every \(c\in\mathbb C\). The arc
\((x,y,z)=(t^{-1},0,2t^2-ct^3)\) escapes with limiting image
\((0,0,c)\). In contrast \((0,1,0)\) has the three distinct preimages
\((0,1,-4),(-2,1,2),(2,-1/2,5/4)\). Their local analytic inverse
branches show \(d\ge3\). These identities and the determinant
were checked in exact symbolic arithmetic. The example is a
control on the inference in arbitrary dimension. The plane
Jacobian obligation and the classification of its actual omitted
divisors remain open. [Retained polynomial control](jacobian/affine-complements/NORMAL_FORM_COVERAGE_AND_COMMON_LEADING_DIVISORS.md#a-mandatory-three-dimensional-control-for-integrality)

## Screened bridges

These dispositions refer to the linked versions and the proposed
receivers in this repository. Each describes the missing hypothesis
or estimate that would have to change before the application is reopened.

| Source | Disposition |
| --- | --- |
| [Rational-Exponential Iterated Integrals and Quantum Polylogarithms](https://github.com/long-mathematics/quantum-polylogarithms-finite-dilogarithms/blob/9f5e6d578c0dd3f649a244fcdd241cdab4176427/papers/quantum_polylogarithms.tex) | Fixed positive irrational parameter and fixed offsets give a functional presentation along a common translation. A mapping to the actual factorial or exponential period packet and an arithmetic specialization theorem are both missing. |
| [Finite Quantum Dilogarithms](https://github.com/long-mathematics/quantum-polylogarithms-finite-dilogarithms/blob/9f5e6d578c0dd3f649a244fcdd241cdab4176427/papers/finite_quantum_dilogarithms.tex) | Finite étaleness concerns the localized raw pentagon solution scheme. A transfer would require equations or a map identifying the actual local moment or period scheme with the relevant pentagon data. |
| [ζ(5) irrationality formalization](https://github.com/long-mathematics/zeta5-irrationality/tree/8fbb4469aec74b88946b711f15f18190cb3c880c) | Defer the rank-correction analogy. The source arithmetic and finite margins concern a specific irrationality construction. They supply no admissible correction of the actual RUW1 mixed evaluation rows with a quantitative minor bound and all denominator costs. The pinned repository also leaves its unconditional Lean endpoint open. |
| [The exact Poisson–binomial comparison in every dimension](https://github.com/long-mathematics/poisson-binomial-comparison/tree/b34cde620016c85e5c7d5461630c0fde3f098d4b) | Its extremum compares count tails of products with ordered Bernoulli parameters and a fixed mean gap. Unexposed prime selectors are independent, while the representation indicators are overlapping products of those selectors. The comparison supplies no weighted overlap bound for those products. |
| [Extremizers of the Expected Missing Mass under an Entropy Constraint](https://github.com/long-mathematics/entropy-constrained-missing-mass/tree/cfc6f0e70596b36c2309cccbb931a1c84f63ad7f) | The objective comes from independent samples of a population law under an entropy constraint. Neither that sampling law nor the constrained objective represents the conditional additive collision count. |
| [Collision-constrained occupancy extrema](https://github.com/long-mathematics/collision-occupancy-extrema/tree/90fc23f63fdec6f762082166f39b4cc22b053cb1) | A collision probability of the form sum of squared population weights constrains independent draws. The semigroup needs weighted reuse of shared primes and the additive equation on each configuration. A scalar collision parameter loses that information. |
| [Partial-Population Inference](https://github.com/long-mathematics/partial-population-inference/tree/519ff7d8ae83ac2455d736fa7348afdfcddcdf69) | The recovery and calibration guarantees use specified sampling kernels. The exposed-prime filtration gives a different conditional law and the arithmetic representation intensity still requires its own estimate. |
| [Super-Coupon Collection under Markovian Replacement](https://github.com/long-mathematics/temporal-super-coupons/tree/01bc2a80095c0c59db36566dd28af55c8b381c08) | The fixed-size subsets are observed under a Markov replacement process. Revealing static weighted primes gives no such transition or uniform replacement law for additive representations. |
| [Bounded-batch species discovery](https://github.com/long-mathematics/bounded-batch-discovery/tree/39d7cbbd29651a1c703cea57aab3bc390651c655) | The model allows within-batch dependence, with independent batches of a fixed bounded size. A prime can occur in many representations across a window, so bounded independent batches for the actual observations would require a new construction. |
| [Missing-mass extremizers and canonical Good–Turing risk](https://github.com/long-mathematics/missing-mass-good-turing/tree/0bf8182b32e057b1e5c7c61c9355606e16bedcf6) | The Poissonized species counts and canonical estimator risk do not identify the conditional intensity or control same-completion overlap in the semigroup. |
| [Deterministic von Neumann fair extractor](https://github.com/long-mathematics/deterministic-von-neumann-fair-extractor/tree/6c5b43349cf96973d434d4cf43d6f92568b6c0ef) | An extraction identity for its input process does not give a switch preserving the arithmetic target and the ratios of unequal prime selector weights. |
| [Finite squarefree-obstruction certificates and bounded prime-degree rigidity for completely multiplicative mixtures](https://rstation.ai/problems/riemann_hypothesis/papers/27/v1) | The certificates constrain squarefree support of the Dirichlet inverse of a finite completely multiplicative mixture. There is no corresponding inverse-support hypothesis for the additive prime-semigroup configurations or the powerful survivor counts. |
| [Large values for cubic conductor pairs and their applications](https://rstation.ai/problems/riemann_hypothesis/papers/87/v1) | The large-value and selected-correlation bounds concern its distinct cubic conductor-pair character family. The actual short-window four-form energy would require a uniform application to its moving weights, including repeated conductors and the short-width diagonal. |
| [Ramification divisors and twisted Nambu brackets of Keller maps](https://arxiv.org/html/2611.34216v1) | Proposition 7.1 treats an assumed Keller counterexample. Sections 7.2–7.3 identify the retained higher-dimensional example with a marked simple root of a binary cubic and compute its cover. The description gives no plane directional-conductor exclusion or classification of the missing boundary components of a plane Keller map. |

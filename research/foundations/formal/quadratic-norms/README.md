# Quadratic norms and splitting semigroups in Lean

This project formalizes parts of
[Pointwise Correlations of Quadratic Form Values and Splitting Semigroups](../../../../results/pointwise-binary-splitting-semigroups/)
and [Exact Sumsets of Primitive Quadratic Norms](../../../../results/eisenstein-restricted-prime-sums/).
It counts numerical values, allows signed or zero coordinates, and
keeps primitive coordinates distinct from coprime summands.

## Build and audit

Run from this directory with Elan. The toolchain and Mathlib are pinned
to v4.32.1. Lake records the exact dependency revisions in the manifest.

```sh
lake exe cache get
lake build
lake env lean Audit.lean
```

The audit prints the transitive axioms of the listed theorems.
The arithmetic and certificate soundness proofs use Lean's standard
foundational axioms `propext`, `Classical.choice`, and `Quot.sound`.
The five full-range scans, the two −11 interval evaluations, and the
−8 interval replay use eight `native_decide` axioms. Lean checks the two finite −11 Euler
products by kernel reduction. Lean checks the implications from these
closed Boolean evaluations to the stated coverage and rational bounds.
The remaining mathematical hypotheses appear in theorem signatures.

## Exact Sumsets of Primitive Quadratic Norms

| Manuscript result | Lean source and scope |
|---|---|
| Primitive-value criterion | [PrimitiveClassification.lean](QuadraticNorms/PrimitiveClassification.lean) proves all four identities `P = C ∪ rC` from integer-coordinate definitions. [NormRoots.lean](QuadraticNorms/NormRoots.lean) proves the Euclidean Bézout and congruence-root bridge. [PrimePowerRoots.lean](QuadraticNorms/PrimePowerRoots.lean) proves lifting and gluing over prime powers. |
| Necessary conditions in the four rows of the exact overview | [LocalObstructions.lean](QuadraticNorms/LocalObstructions.lean) proves the congruence and ramified-square obstructions. [FiniteValues.lean](QuadraticNorms/FiniteValues.lean) proves the three −7 exceptions by exhaustive coordinate search. |
| Packed finite verification | [PackedCertificates.lean](QuadraticNorms/PackedCertificates.lean) proves that the bitset inclusion gives a sum witness for each marked target. It checks the full −7 row through 203 with exact kernel computation. |
| Full Eisenstein finite scans | [EisensteinCertificate.lean](QuadraticNorms/EisensteinCertificate.lean) proves small-prime sieve completeness and mask soundness. [EisensteinReplay.lean](QuadraticNorms/EisensteinReplay.lean) replays all 33,333,333 admissible targets below 10⁸ with the published short-anchor limits. It joins this range to the two large reduced-target obligations at exactly 10⁸. [EisensteinAnchorMenu.lean](QuadraticNorms/EisensteinAnchorMenu.lean) proves that a finite anchor list cannot cover either unrestricted target class. |
| Full −7, −8, and −11 finite scans | [PeriodicCoreCertificate.lean](QuadraticNorms/PeriodicCoreCertificate.lean) proves the periodic membership sieve. [MinusSevenReplay.lean](QuadraticNorms/MinusSevenReplay.lean) proves core coverage through 2²⁴, [MinusEightReplay.lean](QuadraticNorms/MinusEightReplay.lean) proves ordinary and mixed coverage through 2²⁷, and [MinusElevenReplay.lean](QuadraticNorms/MinusElevenReplay.lean) proves core coverage through 2²⁵. |
| Sieve and interval estimates | [SieveCriteria.lean](QuadraticNorms/SieveCriteria.lean) proves the finite unspoiled-candidate argument and the final rational margin for the Eisenstein tail. [OrbitAveraging.lean](QuadraticNorms/OrbitAveraging.lean) proves the finite weighted orbit inequalities for one and two eigenfunctions and the constancy supplied by equivariance. [MinusElevenIntervals.lean](QuadraticNorms/MinusElevenIntervals.lean) checks the rational margins, integer cube roots, and prime cutoffs for all 35 dyadic bins. It recomputes the bounded inert-prime counts and proves the fixed-point rounding bound for the finite portion of the tail. [MinusElevenFrontier.lean](QuadraticNorms/MinusElevenFrontier.lean) proves the divisor limits by compressing each split factorization to the initial split-prime prefix and checking a finite exponent search. |
| −8 signed shells and theta coefficients | [MinusEightShell.lean](QuadraticNorms/MinusEightShell.lean) proves dyadic inclusion–exclusion, the ramified bijection for the even-even part, and the smaller-shell formula at inert primes dividing the target. [MinusEightThetaCoefficients.lean](QuadraticNorms/MinusEightThetaCoefficients.lean) proves stabilization of the level-32 eta coefficients and checks both shell coefficient lists through their Sturm cutoffs. [MinusEightThetaNormalization.lean](QuadraticNorms/MinusEightThetaNormalization.lean) derives the four dyadic counts and exact even selected-shell counts from the modular identities. |
| −8 inert-prime sieve | [MinusEightPrimeBounds.lean](QuadraticNorms/MinusEightPrimeBounds.lean) proves the inert-prime count, inverse-square tail, and local mass bound below 3/16. [MinusEightLattice.lean](QuadraticNorms/MinusEightLattice.lean) and [MinusEightRepresentationCounts.lean](QuadraticNorms/MinusEightRepresentationCounts.lean) prove the signed lattice and representation estimates. [MinusEightBadCounts.lean](QuadraticNorms/MinusEightBadCounts.lean) counts the actual large-prime bad subsets. [MinusEightSieve.lean](QuadraticNorms/MinusEightSieve.lean) derives the selected-shell criteria with normalized large-prime coefficients 8/3, 2, and 4/3. |
| −8 uniform ranges and ordered pairs | [MinusEightEulerBounds.lean](QuadraticNorms/MinusEightEulerBounds.lean) proves the fourth- and sixth-power divisor estimates and the sharper inert-core error bound. [MinusEightAsymptotic.lean](QuadraticNorms/MinusEightAsymptotic.lean) proves the margins at the global, inert-core, and split-prime onsets. It also derives the ordered primitive value-pair bound from the surviving coordinate count and the representation multiplicity. |
| −8 interval certificate and exact sumsets | [MinusEightIntervals.lean](QuadraticNorms/MinusEightIntervals.lean) checks all 107 rows, including the divisor envelopes, cube roots, cutoffs, prime counts, and rounded inverse-square tails. [MinusEightCompletion.lean](QuadraticNorms/MinusEightCompletion.lean) joins the intervals to the finite scan, the uniform ranges, and split-prime descent. [MinusEightHeckeOrbits.lean](QuadraticNorms/MinusEightHeckeOrbits.lean) derives the three exact sumsets and the ordered-pair lower bound from the theta identities, unramified coefficient data, and finite Hecke orbit hypotheses. It proves the dyadic kernel estimate and combines it with the odd operator bound. |
| −11 arithmetic sieve and shell reduction | [MinusElevenPrimeBounds.lean](QuadraticNorms/MinusElevenPrimeBounds.lean) proves the inert-prime count, infinite inverse-square tail, and total local mass bounds. [MinusElevenEulerBounds.lean](QuadraticNorms/MinusElevenEulerBounds.lean) proves the uniform split-divisor and sharp shell error bounds. [MinusElevenSieveMargin.lean](QuadraticNorms/MinusElevenSieveMargin.lean) gives a cutoff with positive margin on every inert core above the finite scan. [MinusElevenShell.lean](QuadraticNorms/MinusElevenShell.lean) proves the odd-shell subtraction identity and the smaller-shell identity at inert primes dividing the target. |
| −11 large inert-prime count | [MinusElevenLattice.lean](QuadraticNorms/MinusElevenLattice.lean) proves the signed-coordinate lattice bound of 12N/5. [MinusElevenPrimitiveCounts.lean](QuadraticNorms/MinusElevenPrimitiveCounts.lean) proves that each root identifies at most two primitive signed representations. [MinusElevenRootCounts.lean](QuadraticNorms/MinusElevenRootCounts.lean) proves multiplicativity and the prime-power bounds for the root counts. [MinusElevenRepresentationCounts.lean](QuadraticNorms/MinusElevenRepresentationCounts.lean) groups representations by their coordinate gcd and proves the split-divisor bound and the large inert-prime bound for the actual odd shell. |
| −11 filtered theta and small-prime estimates | [MinusElevenThetaNormalization.lean](QuadraticNorms/MinusElevenThetaNormalization.lean) derives the odd-shell bounds from the exact theta coefficients and the filtered cusp error. [MinusElevenCuspBounds.lean](QuadraticNorms/MinusElevenCuspBounds.lean) derives that error bound from the unramified Hecke recurrence, multiplicativity, and prime coefficient bounds. [MinusElevenDivisorShellBounds.lean](QuadraticNorms/MinusElevenDivisorShellBounds.lean) proves the uniform bound E < 1 and the small-prime estimate when the inert prime divides the target. [MinusElevenHeckeOrbits.lean](QuadraticNorms/MinusElevenHeckeOrbits.lean) derives the ramified and inert unit-level bad-count estimates from finite operator data. |
| Exact classifications through the scan cutoffs | [ExactSumsets.lean](QuadraticNorms/ExactSumsets.lean) proves bounded ramified restoration and the two −7 primitive patches. [FiniteCoverage.lean](QuadraticNorms/FiniteCoverage.lean) assembles the four finite primitive classifications. With [CoreDescent.lean](QuadraticNorms/CoreDescent.lean), it reduces the unrestricted overview to `BeyondFiniteCoreCoverage`. |
| Eisenstein higher sumsets | [HigherSumsets.lean](QuadraticNorms/HigherSumsets.lean) derives the exact `kC` and `kP` residue formulas from the ordinary and mixed core coverage statements. |

The unrestricted classifications still require coverage above the scan
cutoffs. For −8, the arithmetic and interval arguments are checked.
The remaining inputs are the all-degree theta identities, the
unramified coefficient data, and the unit-level Hecke orbit
construction and spectral hypotheses. For −11, `P11_sum_of_theta_and_hecke_orbits` derives the
unrestricted row from the theta coefficient identity, unramified Hecke
coefficient data, and finite Hecke orbit data at eleven and inert
primes coprime to the reduced target. The remaining work includes
proving the theta identity, establishing the coefficient hypotheses,
and constructing the orbit data with the stated masses, spectral
bound, and return identity. The quantitative lower bounds on ordered
pair counts also remain in the written proof.

## Pointwise Correlations of Quadratic Form Values and Splitting Semigroups

| Manuscript result | Lean source and scope |
|---|---|
| Smilansky’s five ordinary Euler indicators | [OrdinaryIndicators.lean](QuadraticNorms/OrdinaryIndicators.lean) and [ConductorTwoIndicators.lean](QuadraticNorms/ConductorTwoIndicators.lean) prove the signed-coordinate value criteria for x²+dy² with d = 1, 2, 3, 4, 7, including the extra dyadic masks. The latter also proves exact vanishing for x²+4y² at shifts congruent to two modulo four. |
| Subgroup and coset finite images and local solubility | [FiniteImages.lean](QuadraticNorms/FiniteImages.lean) proves the fixed-state finite-image criterion and the equivalence between the joint residue condition and solubility modulo every positive modulus. The result applies to arbitrary nonempty allowed sets of unit classes and masks inside the generated subgroups. [PointwiseCounts.lean](QuadraticNorms/PointwiseCounts.lean) counts the surviving residues at a good prime when both summands require units. |
| Singular series | [PointwiseSingularSeries.lean](QuadraticNorms/PointwiseSingularSeries.lean) proves the normalized good-prime factors, the positive absolutely convergent ordinary product, and uniform bounds for both the ordinary tail and the target-prime tail. |
| Bank record law | [PointwiseBankValuations.lean](QuadraticNorms/PointwiseBankValuations.lean) counts exact valuations below the local modulus exponent on the shifted and reflected lines. It proves the joint product law and exact CRT proportions for events in the finite record space. |
| Bank success and discarded records | [PointwiseBankSuccess.lean](QuadraticNorms/PointwiseBankSuccess.lean) gives the exact success formula, the uniform error bound of 40 times the square-reciprocal mass, and the eventual one-eighth lower bound. [PointwiseBankConcentration.lean](QuadraticNorms/PointwiseBankConcentration.lean) groups independent finite records by their success set and proves the discarded-bank exponent for both size ranges. |
| Bank reconstruction | [PointwiseBankReconstruction.lean](QuadraticNorms/PointwiseBankReconstruction.lean) proves the signed binomial coefficient formula, the threshold polynomial, and the transfer from active-bank correlations to the retained correlation. [PointwiseBankComparison.lean](QuadraticNorms/PointwiseBankComparison.lean) combines this with the discarded-bank bound to control the full comparison from the active-set estimates and positive record envelope. |
| Progression character frame | [PointwiseBankCharacters.lean](QuadraticNorms/PointwiseBankCharacters.lean) proves the conductor criterion, unit phases on products, and cancellation on unit residue classes. [PointwiseBankFrame.lean](QuadraticNorms/PointwiseBankFrame.lean) derives the pointwise and mean kernel bounds, reproduction estimate, and operator norm at most two from the stated off-diagonal cell estimate. |
| Independent prime products | [PointwiseBankProducts.lean](QuadraticNorms/PointwiseBankProducts.lean) proves the exact law for ordered prime draws with repetition, the explicit logarithmic bound, and the tuple-to-product energy estimates. |
| Reflected residue energy | [PointwiseResidueEnergy.lean](QuadraticNorms/PointwiseResidueEnergy.lean) proves the exact character decomposition over all gcd profiles, including cofactor modulus one. It also proves the modulus change under dilation, least-common-multiple refinement, and interval occupancy bounds for arbitrary tails and selections of cofactor characters. |
| Fourier comparison and target grids | [PointwiseFourierComparison.lean](QuadraticNorms/PointwiseFourierComparison.lean) proves the one-period Fourier reduction for the weighted comparison and the continuous Parseval identity for distinct endpoint labels. [PointwiseResidueEnergy.lean](QuadraticNorms/PointwiseResidueEnergy.lean) proves shifted-grid Parseval and integrates a uniform residue-energy bound with the exact grid normalization. |
| C₁₂ Euler identity | [C12EulerProduct.lean](QuadraticNorms/C12EulerProduct.lean) proves the identity for the fourth power of the C₁₂ Dirichlet series on Re(s) > 1, including the separate factors at two and three. |
| Coprime C₁₂ pairs | [PointwiseCounts.lean](QuadraticNorms/PointwiseCounts.lean) proves the exact Möbius inversion. [C12CoprimeWeights.lean](QuadraticNorms/C12CoprimeWeights.lean) verifies the signed main factor, the absolute divisor mass, and the logarithmic moment. [C12CoprimeTransfer.lean](QuadraticNorms/C12CoprimeTransfer.lean) gives the finite error bound. [PointwiseDivisorGrowth.lean](QuadraticNorms/PointwiseDivisorGrowth.lean) bounds the three divisor factors by any positive power of log N and controls the large-divisor remainder. [C12CoprimeAsymptotic.lean](QuadraticNorms/C12CoprimeAsymptotic.lean) derives the uniform coprime asymptotic with saving δ/2 and eventual coprime representations from the ordinary estimate with saving δ. It also proves vanishing outside the admissible class. |
| Basis of exact order thirteen | [HigherSumsets.lean](QuadraticNorms/HigherSumsets.lean) proves the deduction from eventual binary C₁₂ coverage, including the unbounded obstruction to order twelve. |
| Genus and class-group comparison | [Orientation.lean](QuadraticNorms/Orientation.lean) proves the finite abelian group step in orientation saturation. Independent reversals through the order-minus-one bound realize each element of the generated subgroup of doubled classes. |

The analytic estimates supplying the bank comparison hypotheses,
the off-diagonal cell estimate, and the uniform logarithmic saving
in the reflected residue energy remain in the written proof. The CK
and Smilansky asymptotics, the subgroup and coset asymptotics, and the
common-field genus asymptotics and sieve error estimates also remain.
The coprime transfer uses the ordinary C₁₂ pair-count estimate as its
analytic input. The growth bounds use Chebyshev’s prime-counting
estimate and suffice for the saving in the corollary.

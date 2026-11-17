# All-degree HC4 graph obstructions and exact reparameterization endpoints

[PF2RULED152](../../jacobian/affine-complements/RODRIGUES_NILPOTENCE_TRANSFER_AND_HC4_LINEARIZED_TAIL.md)
tests the new Rodrigues strategy with its ACTUAL length17 multiplication.
At fixed HC4 determinant4, a quadratic trace perturbation over eta=b^6
has the exact response eta/(4-13*x^2) in one projected coefficient. Thus
J=(eta) has J^2=0 but arbitrarily late coefficients remain in J, without
increasing powers. 93 exact controls pass. This is a transfer guard, not
a refutation on a smaller remaining trace quotient. The next finite gate
is nilpotence of propagation modulo J, then a full recurrence weight proof.
The required HC4 quotient/ideal certificate and global coverage remain OPEN.

[PF2RULED151](../../jacobian/affine-complements/NORMAL_FORM_COVERAGE_AND_COMMON_LEADING_DIVISORS.md)
makes COVERAGE the next priority. Inside marked quadratic models, valid
target reductions give nondividing ruling degrees>=4 and coupled leading
powers. Both-even degrees force a bad common polynomial G; any odd degree
forces BOTH leading coefficients to have a double zero at u=0. This explains
why good-leading coordinate exclusions do not supply coverage. Global factor
existence, cyclic-image classification and removal of those bad divisors
remain OPEN. 68 exact controls pass, including a 3D generic first-velocity
nonintegrality guard. The preprint audit revalidates 276 retained checks;
quintic graph degree means potential degree16, not degree5. No fresh CAS replay.

[PF2RULED150](../../jacobian/affine-complements/ODD_LEADING_CYCLIC_COORDINATES_AND_RECIPROCAL_BOUNDARY_COFORMS.md)
excludes odd leading coordinates w*A+B, deg_z A=n>=1, deg_z B<=n,
on S2 for every admissible leading a; on MX0 it also covers all contacts
with n>=k-1 and every noncontact root. The reciprocal boundary coframe
keeps zero limiting w, higher contacts, logarithmic equality and all
ramification. PF52 closes the remaining actual n=1,2 line cases.
Together with PF149, all constant-leading coordinates of either parity
are excluded on S2 and MX0 with k<=4. 95 exact controls pass. Other
leading divisors, low odd contact cases and full integrality remain OPEN.

[PF2RULED149](../../jacobian/affine-complements/MONIC_PARAMETER_POLE_BOUNDS_AND_EVEN_LEADING_CYCLIC_TARGETS.md)
removes ALL Newton support bounds for P=B(u,z)+w*A(u,z), deg_z B=n>=1,
deg_z A<n, with nonzero CONSTANT leading coefficient of B. All lower
coefficients are arbitrary. Nonconstant leading c is also allowed when
coprime to M with every root multiplicity<=2*n. A monic parameter pole
bound keeps every generic branch at most logarithmic, excluding rational
mates on MX0 and S2. 83 exact controls pass. PF150 now narrows the odd
leading alternative. Leading divisors outside the stated bounds and
full marked integrality remain OPEN.

[PF2RULED148](../../jacobian/affine-complements/CLOSED_NEWTON_FACES_AND_GENERIC_PARAMETER_DERIVATIVE.md)
excludes EVERY mixed coordinate a*w+F(u,z)+g(u) on quadratic MX0 with
closed support (2*n-1)*j+m*i<=n*m, deg_z F=n>=2, leading coefficient
constant nonzero, and ARBITRARY g. All closed-face coefficients and
multiple leading roots are retained. A generic-parameter derivative
bounds the actual differential at every ramified branch; a genuine MX
triple-root control has index3. 83 exact controls pass. PF149 now removes
this support restriction for its even-leading family, including arbitrary
lower w coefficients. General
coordinates and full marked integrality remain OPEN.

[PF2RULED147](../../jacobian/affine-complements/ADDITIVE_BASE_PERTURBATIONS_AND_RESONANT_CYCLIC_TARGET_EXCLUSION.md)
extends the actual coframe obstruction to a*w+f(z)+g(u), deg f>=2,
with g ARBITRARY in degree and coefficients, on quadratic MX0 and S2.
Every infinity branch is retained, including finite nonzero z limits,
leading cancellations and the full double-root resonant locus. The
required differential stays holomorphic, so no rational source mate
exists. 67 exact controls pass. General multiplicative mixed terms and
full marked integrality/HC4/JC2 remain OPEN.

[PF2RULED146](../../jacobian/affine-complements/HOLOMORPHIC_GENERIC_FIBRES_AND_NONLINEAR_CYCLIC_TARGET_EXCLUSION.md)
excludes EVERY nonlinear separated coordinate a*w+f(z), a!=0, on quadratic
MX0, without target symmetry or any bound on its companion. The actual
generic-fibre coframe is nonzero holomorphic at EVERY infinity valuation,
so it has no rational primitive even in the plane source field. 66 exact
controls pass. PF145's regular linear-z exclusion is already PF52's line
criterion; its critical formulas remain valid. General mixed coordinates
and full HC4/JC2 remain OPEN.

[PF2RULED145](../../jacobian/affine-complements/JACOBIAN_MATE_RIGIDITY_AND_QUADRATIC_CYCLIC_TARGET_EXCLUSION.md)
retains exact critical schemes and mate rigidity. Its regular linear-z
target exclusion already follows from PF52's reached source-line criterion,
also with polynomial coefficients of z. The proofs and 90 controls remain
valid; PF146 addresses the genuinely nonlinear separated class.

[PF2RULED144](../../jacobian/affine-complements/CYCLIC_FIBRE_SCALING_AND_TARGET_FIELD_ESCAPE.md)
excludes EVERY plane target on the MX0 cyclic family whose function field
is preserved by any nontrivial regular fibre scaling. The odd-prime fixed
curve/point obstruction and involution source-line argument allow arbitrary
target degree and nonlinear target action. A remaining target must break
every such scaling. The reached-line criterion also covers all classical S_g.
101 exact controls pass; unrestricted cyclic targets,
PF143 marked regularity and full HC4/JC2 remain OPEN.

[PF2RULED143](../../jacobian/affine-complements/MARKED_QUADRATIC_VELOCITIES_AND_CYCLIC_TRACE_POLES.md)
gives a finite FULL-S-marked quadratic test for a generic first velocity, with
exact charge 1+max ceil((2q+1)/e) over every omitted prime. Regularity is
equivalent to inversion and remains unsupplied. In the surviving MX0 cyclic
family, an actual target coframe forces polar relative first/second moments;
for g>=5 their partial Gram determinant is nevertheless constant and nonzero.
Integral full-normalization markers retain the obstruction in the original
target field. All ramified primes and PF142's mixed minimal retirement remain.
92 exact controls pass; full HC4/JC2 remains OPEN.

[Logarithmic descent and cyclic exact source images](LOGARITHMIC_DESCENT_AND_CYCLIC_EXACT_SOURCE_IMAGES.md), HC4MYMX0, identifies the logarithmic quotient with Pic(V), with closed source forms seeing only torsion. Explicit cyclic surjective plane-source images have exact volume and Pic(V)=Z/g plus a free group of arbitrary rank; no plane target or finite trace-lattice stability is supplied.

[Affine étale primitive descent](AFFINE_ETALE_DESCENT_OF_RATIONAL_PRIMITIVES.md), HC4MYMW0, completes the all-degree residue criterion for every permitted rational slice. It corrects MV0’s independent polynomial-representative gap by splitting the regular cocycle on the actual source relation. Exactness still leaves the target coframe and full trace-lattice problem open.

[PF2RULED142](../../jacobian/affine-complements/PARTIAL_KUMMER_PLANE_CHARTS_AND_CYCLIC_EQUIVARIANT_FACTORS.md) closes MR0/PF139's mixed replica-equivariant alternative at MINIMAL generic degree, for arbitrary source LNDs. A selected open of a branched Kummer cover is an etale plane and gives a smaller Keller comparison whenever h is not a shifted pure power. Every such proper factor of a minimal counterexample is therefore an actual cyclic deck fixed field. MU0's arbitrary-degree conditional order/trace results remain valid, but the mixed degree-five factor now descends from 5b to 2b if realized. No automatic equivariance or full trace-energy landing is proved.

The [trace-free rational-slice theorem](TRACE_FREE_RATIONAL_SLICES_AND_MIXED_SOURCE_IMAGE_GRAPH_MODELS.md), HC4MYMV0, gives a global primitive criterion; MW0 now proves its polynomial source representative automatic from a rational primitive under the actual smooth affine surjective étale hypotheses. It explicitly realizes graph source images of the mixed fibre, then excludes every multi-root graph from an étale plane target by an all-degree exactness obstruction. General rational mixed slices remain open.

The [source-velocity target multiplier ideal](../../jacobian/affine-complements/GENERIC_SOURCE_VELOCITY_MULTIPLIER_IDEAL_AND_TRACE_POLE_CHARGES.md), PF2RULED141, measures MS0/MT0's missing full finite landing exactly: generic velocity pole q+e gives target exponent1+max ceil(q/e), also recovered from the first d power traces. A full source-finite inverse leading coefficient squared clears all four source Jacobian entries and stabilizes S/Lambda in every degree, with divergence accounted for. These rescaled fields are not locally finite under failure; unscaled trace/energy regularity remains open.

The [first mixed-fibre theorem](FIRST_MIXED_FIBRE_SOURCE_IMAGE_AND_QUADRATIC_TRACE_POLE.md), HC4MYMU0, retains the rational slice for h=y²(y−1)³, separates the actual finite order from its primitive order, and forces a quadratic relative trace pole on the multiplicity-three component. An integral intermediate marker detects an original-target trace pole; no Keller realization or monodromy-only exclusion is claimed.

The [all-intermediate missed-line theorem](../../jacobian/affine-complements/COMPLETELY_MISSED_INTERMEDIATE_LINES_AND_SOURCE_IMAGE_EULER_BUDGET.md), PF2RULED140, supplies KV0's previously unspecified deleted-curve geometry: disjoint smooth A1s, full isomorphic power-map packets and zero excess. It gives b=chi(T)+r+A+Eta for the actual maximal source image, while preserving general finite point gaps and the entire-locus distinction. MR0/PF139 mixed fibres therefore require quotient degree b>=r+5. MS0's first-source-derivative finite-lattice tests remain unproved automatic.

The [generic source-translation theorem](GENERIC_SOURCE_TRANSLATION_DERIVATIVES_GENERATE_THE_POLYNOMIAL_RING.md), HC4MYMT0, strengthens field recovery to J_v=C[x,y] and a closed embedding for all but finitely many linear directions in every source frame. Finiteness over the original target and trace-lattice stability remain open.

The [replica degree-descent theorem](../../jacobian/affine-complements/REPLICA_EQUIVARIANT_DEGREE_DESCENT_AND_THE_UNIQUE_MIXED_FIBRE.md), PF2RULED139, combines shared MR0's actual equivariant source image and PF130's general chart construction. In minimal generic degree, every square-free replica has no proper stable intermediate except an actual quadratic deck field with centred odd replica; even-degree square-free replicas have none. General survivors require one mixed normalized fibre of degree>=5 and a repeated replica root. The actual finite trace/energy goal remains open.

The [first-source-derivative criterion](FIRST_SOURCE_DERIVATIVES_AND_FINITE_LATTICE_STABILITY.md), HC4MYMS0, reduces a source LND’s full finite-module/trace-lattice stability to integrality of just deltaP and deltaQ. This integrality remains unproved; finite-field rational resolvents do not supply target-ring regularity or the Mathieu period hypotheses.

The [whole-atlas exactness criterion](../../jacobian/affine-complements/TRANSLATION_ATLAS_EXACTNESS_AND_UNBOUNDED_MOVING_LABEL_FIBRES.md), PF2RULED138, completes the transition-residue test for the actual MQ0 affine translated-source cover. Explicit whole-plane exact-volume free-Ga surfaces have arbitrarily many reduced degenerate fibres, so these properties alone cannot give a fibre-count obstruction. Their compatible etale targets remain unconstructed; actual full trace-dual multiplication/marked-energy remains the goal.

The [replica-degree theorem](REPLICA_DEGREE_AND_SOURCE_DERIVATIVE_FIELD_RECONSTRUCTION.md), HC4MYMR0, proves that a source LND p(y) partial_x can preserve an intermediate field only with index dividing deg(p)+1. Bare translation derivatives of P,Q generate the full source field in every coordinate frame. Global finite-lattice and trace-energy control remain open.

The [separated-variable volume obstruction](../../jacobian/affine-complements/SEPARATED_VARIABLE_VOLUME_AND_SINGLE_ROOT_OBSTRUCTION.md), PF2RULED137, excludes all a(x)y=p(z) models with multiple distinct roots of a, even all repeated, and linear a, for square-free p of degree>=2. This extends the retained monomial period/primitive boundary without classifying general additive surfaces; single-root higher powers and the full trace gate remain open. MP0/MQ0 are credited as the current actual deck/source-image reductions.

The [cyclic source-image theorem](CYCLIC_DECK_SOURCE_IMAGE_AND_THE_UNIQUE_MULTIPLE_FIBRE.md), HC4MYMQ0, eliminates PF136’s possible finite point gaps: the actual plane source surjects étale onto an affine ruled quotient with one irreducible multiple fibre and reduced fibres elsewhere. The remaining finite-target/trace-energy condition is still open.

The [actual translated-chart reduction](../../jacobian/affine-complements/COMMON_OMITTED_ORBITS_AND_ALMOST_SURJECTIVE_DECK_QUOTIENTS.md), PF2RULED136, preserves the actual target while deleting common omitted curve orbits. MP0 now eliminates its noncyclic-deck survivor cases, and MQ0 removes finite point gaps from the actual cyclic source image. Its sealed conditional Euler bounds and the MO0/MA0 scope distinctions are retained.

The [automatic additive-completion theorem](COMMON_UNRAMIFIED_PENCIL_AND_AUTOMATIC_ADDITIVE_COMPLETION.md), HC4MYMP0, combines the retained common unramified source pencil with the whole boundary charts. Every actual canonical etale completion admits source translations, so every generic rational deck group is cyclic in all degrees. The cyclic/trivial-deck finite-target and full trace-energy problems remain open.

The [cotangent deck/trace reduction](../../jacobian/affine-complements/COTANGENT_DECK_NORMAL_FORM_AND_RELATIVE_TRACE.md), PF2RULED135, gives an explicit invariant field and relative trace for every noncyclic generic deck action, using the actual target primitive. A free affine cotangent control has an exact coframe and regular etale target but no whole plane chart. The remaining gate is actual global finite-order/trace-dual regularity across mixed fibres; rational Darboux coordinates do not establish it.

The [Euler-two multiplicity-pair theorem](COMPLETE_RULING_PLANE_CHARTS_AND_THE_LOW_EULER_MULTIPLICITY_PAIR.md), HC4MYMO0, adds a complete-ruling plane-chart criterion. Comparable pairs force universal cyclic relative compression at minimal degree; coprime normalized pairs with both entries at least two and the cyclic finite-target case remain open.

The [finite orbifold-group reduction](../../jacobian/affine-complements/FINITE_ORBIFOLD_GROUPS_AND_UNIVERSAL_COMPLETE_RULINGS.md), PF2RULED134, restricts finite intermediate pi1 and generic deck groups to cyclic, dihedral, A4, S4 or A5. A noncyclic deck group forces an actual equivariant complete A1-ruling on U_N. After the universal orbifold cover, exact volume requires a fibre whose multiplicities all exceed their gcd; the Euler3 case is one mixed coprime two-component fibre. This structures the rigid branch without giving a Ga action or closing the finite-target gap.

The [deck-normalized Euler bound](../../jacobian/affine-complements/DECK_NORMALIZED_EULER_AND_UNRESTRICTED_RADICAL_TOWER_EXCLUSION.md), PF2RULED133, excludes the ENTIRE natural radical-tower family q>=2 in ALL inertia types. A noncyclic generic deck subgroup requires fixed-field degree at least5, while these actions have degree4. Residual ramification lowers the actual deck quotient Euler and does not rescue a low-Euler factor. The index5 natural-S5/zero-excess profile and the full finite-target gap remain OPEN.

The [affine-monodromy deck-index theorem](AFFINE_MONODROMY_DECK_INDEX_AND_UNDECORATED_RADICAL_TOWERS.md), HC4MYMN0, reuses the audited PF132 low-Euler lemma to exclude the FULL natural radical towers H_q for every q>=2, with any inertia. More generally noncyclic generic decks in binary affine sheet actions need at least eight orbits. Cyclic/trivial decks and the full trace-energy problem remain open.

The [Wright/knot tower integration](../../jacobian/affine-complements/LOW_EULER_CYCLICITY_AND_DECORATED_RADICAL_TOWER_EXCLUSION.md), PF2RULED132, combines HC4MYML0 with the Euler budget: the specified lifted-transposition radical towers are excluded in ALL degrees16,32,64,... without minimality. The general low-Euler exact-volume lemma forces cyclic pi1 and improves the noncyclic generic deck cost to b0>=3|G|-1. PF133 above extends the same natural action to all inertia; q=1 and the finite-target gap remain OPEN, and knot representations stay valid.

The [universal inversion-term theorem](UNIVERSAL_ABHYANKAR_GURJAR_NONTERMINATION_AND_MOMENT_TRANSFER.md), HC4MYMM0, gives every plane Keller extension a framing in which every grouped Abhyankar–Gurjar term is nonzero for every polynomial amplitude. Polynomial sums can still result. Source-origin Gaussian safety and binary factorial results therefore need a separate global transfer before they can address HC4.

The [positive Kähler route audit](POSITIVE_KAHLER_COMPLETENESS_AND_POLYNOMIAL_SQUARES_AUDIT.md) separates automatic flatness from unproved completeness. Finite polynomial-square enlargement retains a square Keller map; recent Liouville growth and real-convexity hypotheses are not automatic, even for known nonlinear inverses.

The [complete-ruling plane-chart audit](COMPLETE_RULING_PLANE_CHART_AND_FINITE_TARGET_AUDIT.md) shows that the retained higher-power Danielewski atlas has a complete ruling with its entire nonreduced fibre outside a plane chart. The missing étale target remains essential. ML0 also cannot be transferred to the finite normalization where the pulled-back volume has ramification zeros.

The [exact-volume complete-ruling obstruction](EXACT_VOLUME_AND_COMPLETE_RULING_COMPONENTS.md), HC4MYML0, forces some fibre of every complete-type affine-line ruling of the full canonical completion to have no multiplicity-one component. A two-chart cotangent-torsor residue proves the obstruction. It does not force an additive action or a whole-fibre gcd greater than one.

The [finite-jet, Bézoutian and moduli audit](FINITE_JETS_BEZOUTIAN_DESCENT_AND_MODULI_OPENNESS_AUDIT.md) locates the same global gap in finite diagonal descent and uniform inverse degree. Every source jet retains the full trace defect; every nilpotent parameter jet is invertible already. Connectedness of the fixed-degree Keller moduli supplies no open automorphism neighbourhood. No new HC4 theorem is registered.

The [universal Pascal-energy framing](UNIVERSAL_KELLER_FRAMING_AND_PASCAL_ENERGY_NONTERMINATION.md), HC4MYMK0, puts every plane Keller extension into normalized HC4 forms with exponential energy degrees and all-order Pascal nontermination, preserving the actual target ring and full trace data. Polynomial inverse and Legendre transform coexist with this behavior; nontermination alone does not detect the unresolved gap.

The [Ricci-flat HC4 theorem](RICCI_FLAT_HC4_AND_ONE_DIMENSIONAL_JACOBIAN_IMAGES.md), HC4MYMJ0, proves inversion in every degree when the Hessian metric’s full Levi-Civita Ricci tensor vanishes, using retained Tang. Its global plane lemma makes every one-dimensional Jacobian-matrix image an affine shear. Scalar flatness and affine/Kahler Ricci vanishing do not supply this hypothesis.

The [Hessian connection-pencil theorem](HESSIAN_CONNECTION_PENCIL_AND_INCOMPLETE_AUTOMORPHISMS.md), HC4MYMI0, gives arbitrary-degree polynomial gradient automorphisms whose Hessian metrics are curved and incomplete. Both flat affine endpoints are complete, yet every nonendpoint connection is incomplete; metric completeness and interpolation do not supply the missing inverse theorem.

The [Hodge multiplication audit](HODGE_MULTIPLICATION_AND_NONPROPER_IMAGE_AUDIT.md) rules out a formal purity-to-trace-integrality shortcut: pure complementary Kummer characters multiply to a polar integer character. It also distinguishes canonical Hodge filtration from naive nonproper pushforward. No new HC4 theorem is registered.

The [iterated-deficit transport theorem](ITERATED_DEFICIT_TRANSPORT_AND_FINITE_CENTER_LABELS.md), HC4MYMH0, retains all missing fiber lengths and finite collision corrections. The exact transport reproduces the Euler budget at every iterate; finite-center ramification labels persist. Neither iteration nor stable-image arguments here force trace-pairing regularity.


The [dynamical-degree and positive-thinness theorem](INTEGER_DYNAMICAL_DEGREE_AND_POSITIVE_THINNESS.md), HC4MYMG0, derives integer growth from Favre-Jonsson and gives target automorphisms with arbitrarily fast growth but the identical finite cover. Every omitted prime has positive thinness and lies outside the tight dynamical tree; full trace-pairing regularity remains open.


The [radial-foliation route audit](RADIAL_FOLIATION_AND_PROJECTABILITY_AUDIT.md) proves that a smooth P1 first integral of a nowhere-zero field with constant nonzero divergence reconstructs an omitted-point Keller pair. This classification lead is the same unresolved problem; source projectability is already excluded under failure by retained line and radial rigidity.


The [Picard and quantization route audit](PICARD_AND_QUANTIZATION_ROUTE_AUDIT.md) explains why projective puncture classes remain nontorsion and why a unique quantization of the whole Keller graph does not give an inverse. It also records Dodd v4’s correction of an older existence claim. No new mathematical claim is adopted.


The [isochronous classification source audit](ISOCHRONOUS_CLASSIFICATION_SOURCE_AUDIT.md) records the primary 2023 retraction of the 2022 claimed plane-Jacobian proof. Composed shears also refute its literal abstract linear normal form. No mathematical claim or inversion theorem is newly adopted.


The [generic exact-mark and Liouville control](GENERIC_EXACT_MARKS_AND_SMOOTH_LIOUVILLE_CURVES.md), HC4MYMF0, gives whole-plane submersions with trivial generic connection and both generic marks zero, but nonzero raw volume/Liouville torsion and no polynomial mate. Its two-parameter family also realizes every genus by smooth closed Liouville-exact plane curves.


The [target radial lattice theorem](TARGET_RADIAL_FINITE_LATTICES_AND_EIGENPOLYNOMIALS.md), HC4MYME0, identifies the lifted target Euler field’s locally finite polynomial part exactly with B. A full stable-order discriminant argument supplies the all-degree proof. The binary-affine immersive profile also inherits PF59’s nonpositive normal-degree restriction;the one-cusp profile does not factor through the Orevkov invariant.


The [smooth projective-pencil control](SMOOTH_PROJECTIVE_PENCILS_WITH_ARBITRARY_GENUS.md), HC4MYMD0, retains finite fibers, exactly one omitted target point, and a smooth surjective P1 pencil with every fiber irreducible, yet has arbitrary generic genus. The full constant-Jacobian mark remains essential.


The [Kahler bracket and point-defect audit](KAHLER_POISSON_IMAGE_AND_POINT_DEFECT.md), HC4MYMC0, shows why a finite point defect can produce an exceptional curve in an already normal bracket algebra. Canonical module hulls must precede algebra generation; generic Poisson identities do not close the trace-pairing gap.


The [exact-collision route audit](EXACT_COLLISION_ROUTE_AUDIT.md) retires the queued bare Lagrangian shortcut as already covered by LC0/1. Its exploratory controls preserve the missing-target and finite-equivalence distinctions; no new theorem is adopted.


The [source base-change audit](SOURCE_BASE_CHANGE_AND_RESIDUAL_DIFFERENT.md), HC4MYMA0, retires the different-pullback shortcut: all discriminant defect survives in a residual finite cover. Connected, surjective constant-Jacobian open controls retain sheets over every target point and show why retained sheets alone do not repair it.


The [exact even-moment lattice formula](EVEN_MOMENT_LATTICES_AND_SIXTH_MOMENT_TRACE_TEST.md), HC4MYLZ0/1, includes singular points and gives a uniform sixth-moment trace criterion detecting every ramification index. Its trace form is automatically regular; perfection, equivalently inversion, remains open.


The [trace-linking and half-ramification audit](TRACE_LINKING_HALF_RAMIFICATION_AND_CUBIC_MULTIPLICATION.md), HC4MYLY0/1, shows why a perfect selfdual trace lattice or a zero Witt class can retain a nonzero discriminant module. Regular cubic multiplication data would close the selfdual-lattice gap; their actual-source regularity remains open.


The [Nambu first-derivative lattice and finite-core theorem](NAMBU_SOBOLEV_LATTICE_AND_FINITE_DERIVATIVE_CORE.md), HC4MYLX0/1, identifies each lifted target additive field's finite trace-lattice defect exactly with the discriminant quotient. Its all-orders polynomial Bergman core is only the target ring; arbitrarily long finite prefixes survive under putative failure. The actual global derivative estimate remains open.


The [Bergman trace-dual and multiplication-domain theorem](BERGMAN_TRACE_DUAL_AND_MULTIPLIER_DOMAIN_RIGIDITY.md), HC4MYLW0/1, identifies the actual finite trace dual with polynomial L2 vectors. Hilbert trace retracts automatically; L4 multiplication remains exactly the inverse gap. A nonintegral polynomial multiplier has a nondense domain. Laurent and Runge controls mark the missing whole-plane Keller improvement.


The [integral source-product and complete-replica theorem](INTEGRAL_SOURCE_PRODUCTS_AND_COMPLETE_HAMILTONIAN_REPLICAS.md), HC4MYLV0/1, classifies all divergence-free locally finite stabilizers under Keller failure and shows why complete holomorphic replicas evade that obstruction. Two noncommuting complete fields survive on retained proper exact whole-plane controls; actual finite-target compatibility remains open.


The [Miyanishi parity source audit](MIYANISHI_PARITY_SOURCE_AUDIT.md) records the published odd-plane lead and an unverified quotient smooth-locus bridge. It supplies no centered-even HC4 closure; the connected open-source control explains why local equivariance alone does not settle that bridge.

The [exact-volume balanced source-torus rigidity theorem](EXACT_VOLUME_AND_BALANCED_SOURCE_TORUS_RIGIDITY.md), HC4MYLU0/1, proves that one balanced source action preserving the full finite normalization or trace dual forces inversion, without target-field descent. The required preservation remains open; unbalanced actions and finite-type completions have different conclusions.


The [Hamiltonian stabilizer and boundary-jet theorem](HAMILTONIAN_STABILIZERS_AND_INVISIBLE_BOUNDARY_JETS.md), HC4MYLT0/1, classifies every preserving divergence-free source field. Arbitrary finite source jets and global tangent spanning coexist with a non-Noetherian stabilizer and an obstruction to any preserving unit-Jacobian pair; these local tests do not supply the missing global landing.


The [raw p-curvature and finite Gauss-Manin comparison](RAW_P_CURVATURE_AND_FINITE_GAUSS_MANIN_COMPARISON.md), HC4MYLS0/1, proves that raw p-curvature zero is automatic before a Keller mate and supplies an all-prime geometric failure of finite-model comparison. Katz's published theorem gives an inverse criterion for the finite model; its required arithmetic vanishing remains open.


The [actual irreducible target frame and puncture-extension constraint](IRREDUCIBLE_TARGET_FRAMES_AND_PUNCTURE_EXTENSIONS.md), HC4MYLR0/1, adopts Kaliman's published reduction while preserving the finite trace dual exactly. A counterexample in that frame requires nonsplit infinite affine fibre monodromy; compact semisimplicity does not remove this extension.


The [component discriminant and canonical-boundary unit constraints](COMPONENT_DISCRIMINANTS_AND_CANONICAL_BOUNDARY_UNITS.md), HC4MYLQ0/1, adopt the published all-degree minimal-discriminant inverse criterion and rule out specified canonical puncture configurations in all-irreducible component frames. Neither discriminant minimality nor the required canonical section is proved automatic.


The [submersion connection and volume mark](SUBMERSION_GAUSS_MANIN_AND_MARKED_VOLUME_OBSTRUCTION.md), HC4MYLP0/1, shows that transverse acyclicity is automatic before a Keller mate. Whole-plane geometric examples with unbounded genus and finite monodromy still have nonzero quotients; their volume marks are nonzero and no mate exists. The actual marked Keller implication remains open.


The [full-normalization additive stabilizer](FULL_NORMALIZATION_ADDITIVE_PENCIL_AND_TRACE_DUAL_STABILIZER.md), HC4MYLO0/1, classifies all extending source additive fields in one commuting replica ideal, or none. Two distinct invariant fields preserving the full normalization or finite trace-dual force inversion, including ramified boundary. Existence of those fields remains open.


The [contracting finite-lattice criterion](CONTRACTING_SOURCE_GRADINGS_AND_FINITE_LATTICES.md), HC4MYLN0/1, proves in every dimension and degree that one positive source grading preserving a nonzero finite target module forces properness, hence a Keller inverse. In the plane this makes trace-dual invariance equivalent to regular trace pairing; the required invariance remains open.


The [canonical-completion additive restriction](ADDITIVE_ACTIONS_ON_THE_CANONICAL_ETALE_COMPLETION.md), HC4MYLM0/1, proves that a proper completion has one additive ruling or none. A noncyclic generic rational deck group forces the completion to have no additive action. This concerns deck symmetry, not full monodromy; ramified boundary and full HC4 remain open.


The [common-coordinate and HC4 orbit filter](COMMON_COORDINATE_OUTPUTS_AND_HC4_ORBIT_FILTER.md), HC4MYLK0/1, excludes EVERY gamma-level map from becoming an HC4 gradient under arbitrary polynomial source/target changes, and excludes every global 2+2 mixed graph chart. Generic field degree three already blocks that route; further Darboux calculations for this fixed map are not an HC4 completion gate. Full HC4 remains open.


The [primitive-gamma coordinate construction](PRIMITIVE_GAMMA_COORDINATE_AND_DARBOUX_BARRIER.md), HC4MYLJ0/1, gives explicit polynomial affine-four-space sources at EVERY moment level. This supersedes the gamma source question below. The induced exact symplectic form still lacks standard Darboux coordinates; neither natural position pair can be retained. Full HC4 remains open.


The [cotangent incidence obstruction](COTANGENT_INCIDENCE_HODGE_OBSTRUCTION.md), HC4MYLI0/1, rules out the entire historical primitive-beta source at every moment level and after every affine stabilization. The all-degree Hodge formula also excludes general target-line incidence sources; the special gamma chart passes only this necessary test. Full HC4 remains open.


The [projective differential-generator sieve](PROJECTIVE_DIFFERENTIAL_GENERATORS_AND_MONODROMY.md), HC4MYLH0/1, turns full quadratic differential equations into constant cross-ratios and cyclic projective stabilizers. A primitive noninvertible plane source can then have only odd dihedral or specified degree6/12/20/30 actions. The differential condition is not automatic; full HC4 remains open.


The [uniform integration criterion](UNIFORM_INTEGRATION_AND_EXACT_RELATIVE_QUOTIENT.md), HC4MYLG0/1, proves that one nonzero multiplier integrating the whole source forces inversion. The raw relative quotient has a bijective transverse derivative and exact special-fibre torsion; this acyclicity does not force it to vanish. Full HC4 and automatic uniform integration remain open.


The [trace-averaging flow barrier](TRACE_AVERAGING_COMPLETE_FLOW_BARRIER.md), HC4MYLF0/1, gives an every-odd-degree counterfamily: a complete source LND fixing ramification and preserving exact target volume has an incomplete traced Hamiltonian target field. Trace also fails brackets, flow laws and branch invariance. The controls have etale source Gm times A1; full HC4 and actual whole-plane descent remain open.


The [scalar and one-component descent theorem](SCALAR_DIFFERENTIAL_RIGIDITY_AND_ONE_COMPONENT_DESCENT.md), HC4MYLE0/1, reduces complete-flow inversion to one component derivative plus divergence. It also makes the entire source group preserving one component pencil finite abelian under failure. Existence of a qualifying field or derivative landing, full HC4 and trace regularity remain open.


The [all-conjugate localization theorem](ALL_CONJUGATE_LOCALIZATION_AND_UNIFORM_TRACE_LANDING.md), HC4MYLD0/1, identifies the invariant base as exactly B[1/f] on the proper locus. It distinguishes finite coefficient denominators from uniform trace landing for the whole actual source. Full HC4 and the source-to-trace-regularity implication remain open.


The [exact collision-action theorem](EXACT_COLLISION_ACTION_AND_FIRST_JET_RECONSTRUCTION.md), HC4MYLC0/1, reconstructs all ordered collision chart rings from first chord derivatives, including ramification. The whole pair of actual source-plane subrings recovers the original target by intersection. Primitive action orders may remain singular; full HC4 and the trace-dual integrality implication remain open.


The [global collision-chord theorem](GLOBAL_COLLISION_CHORDS_AND_PENCIL_DEFECT_DROP.md), HC4MYLB0/1, proves that EVERY complete source pencil has generic affine line-image defect>=2 and a strict defect drop somewhere.The actual collision surface retains scheme multiplicities and needs no nontrivial generic deck group.Full HC4 and the universal trace-dual integrality implication remain open.


The [generic deck boundary theorem](GENERIC_DECK_GROUP_AND_UNRAMIFIED_BOUNDARY_COST.md), HC4MYLA0/1, proves |Aut_F(L)| divides gcd(d,1+b0),where b0 counts omitted unramified curves.NONCYCLIC generic deck symmetry requires b0>=2|Aut_F(L)|-1.Sharp published cyclic pseudo-plane quotients retain actual open A2 source charts and exact volume;the required plane target remains missing.Full HC4 is open.


The [complete holomorphic flow theorem](COMPLETE_HOLOMORPHIC_FLOW_DESCENT.md), HC4MYKZ0/1, extends the source-descent obstruction to ALL complete polynomial fields,including non-locally-finite fields with transcendental time maps.The whole nonproperness curve admits no complete polynomial logarithmic field under plane failure.No such field is proved to exist automatically;full HC4 remains open.


The [target-boundary symmetry theorem](TARGET_BOUNDARY_SYMMETRY_AND_FINITE_STABILIZER.md), HC4MYKY0/1, proves that the ENTIRE nonproperness curve of a noninvertible plane Keller map has FINITE FULL polynomial ambient stabilizer.Target Ga/Gm actions lift algebraically,extend over the normalization and preserve all omitted primes;the existing descent obstruction forces inversion.Individual component symmetries do not suffice.Full HC4 remains open.


The [algebraic-flow descent theorem](ALGEBRAIC_FLOW_DESCENT_AND_INTEGRAL_COORDINATES.md), HC4MYKX0/1, excludes EVERY nonzero locally finite source derivation preserving a noninvertible plane Keller target field,including all polynomial LND replicas.It also computes the exact translation multiplier ideal on the full normalization under an integral source coordinate.Normalization actions do not supply target descent;full HC4 remains open.


The [source-line conductor and primitive-collision theorem](SOURCE_LINE_CONDUCTOR_AND_PRIMITIVE_ORDER_COLLISIONS.md), HC4MYKW0/1, identifies the FULL source conductor and its exact specialization to every admissible ruling line.Noninvertibility forces positive even conductor degree on each line.A primitive integral source coordinate in the KNOWN3D Keller map has a singular cubic order at a fully retained etale value,so local inverses cannot justify primitive-hypersurface smoothness.Full HC4 remains open.


The [differential subalgebra classification](DIFFERENTIAL_SUBALGEBRA_CLASSIFICATION.md), HC4MYKV0/1, proves that every target-differential subalgebra inside an actual plane source is one of finitely many normal affine etale rings,selected by its intermediate field and allowed unramified boundary poles.No finite-generation assumption is needed.A generic linear source coordinate generates the whole source under target differentiation.The final HC4 trace-regularity implication remains open.


The [canonical trace-dual reconstruction theorem](CANONICAL_TRACE_DUAL_ETALE_RECONSTRUCTION.md), HC4MYKU0/1/2, proves that B[Lambda] is the ENTIRE affine source etale locus of any finite locally free generically etale algebra. For every actual plane intermediate it also equals the algebra generated by all target derivatives of the finite normalization and its IC quadratic image. Full HC4 remains open:the actual polynomial source must still force Tr(Lambda*Lambda) to be regular. External review is deferred.


The [Kulikov coframe control and Tate filter](KULIKOV_EXACT_COFRAME_AND_TATE_BOUNDARY_FILTER.md), HC4MYKT0/1, show why an actual plane source remains essential:the KNOWN nonplane etale UFD control also has exact coordinates,simple N and quadratic reconstruction.Its four-punctured boundary gives weight4 H2;actual plane completions have weight2 boundary classes.This filters relaxed models,not the already-qualified hard plane packets.


The [trace-zero algebra theorem](TRACE_ZERO_ALGEBRA_AND_FIRST_ENERGY_DEFECT.md), HC4MYKS0/1, identifies the quadratic image with a normal finite-type smooth etale algebra B[N],the maximal partial completion preserving N.Higher products add nothing.Its Picard basis is the single-unramified-lost-sheet defect,and a regular first energy trace removes that defect without proving inversion.


The [quadratic source theorem and published quartic exclusion](QUADRATIC_SOURCE_RECONSTRUCTION_AND_QUARTIC_RETIREMENT.md), HC4MYKR0/1/2, retire ALL generic-degree-four plane Keller realizations by Domrina2000. In every degree>=3,the trace-zero quadratic image recovers the source except for single unramified lost-sheet curve quotients;one boundary prime suffices for complete reconstruction. Actual embedded multiplication and polynomial generators remain unresolved.


Every native polynomial graph eliminating any one of x,y,p,q,r from
the scaled Meng--Yang potential has invertible gradient whenever its
actual four-variable Hessian determinant is a nonzero constant. Graph
degrees are unrestricted. The reparameterization extensions differ:

- **The full plane source has two exact semisimple layers:** [HC4MYKQ0/1/2](SOURCE_MODULE_LAYERS_AND_TRACE_ZERO_KERNEL.md) uses retained Borisov geometry to identify its IC socle and boundary quotient. The trace-zero kernel is semisimple exactly for one residue-one omitted prime per target image;with one boundary and doubly transitive monodromy it is already simple,and the source is indecomposable of length3. The actual extension and polynomial multiplication remain open. This audit also credits PF2BATT1 for the earlier nodal geometric exclusion.
- **The plane trace image is rigid,and nodal inertia cannot overlap:** [HC4MYKP0/1/2](PLANE_TRACE_IMAGE_RIGIDITY_AND_NODAL_INERTIA.md) fixes the complete differential trace image from the target nonproperness curve alone. At a node it has exactly one point defect and no mixed poles,forcing the two meridians to move disjoint sheets. The geometric nodal exclusion already follows from retained PF2BATT1;this supplies a separate trace derivation and explicit support bound. Quartic and degree-six disjoint-support packets remain compatible. The full actions,marked energy regularity and HC4 remain open.
- **The canonical energy is now a single intrinsic HC4 gate:** [HC4MYKO3](CANONICAL_ENERGY_INTEGRALITY_AND_POWER_TRACES.md) proves that G=x dot grad(h)-h is integral over the gradient target iff the gradient has an inverse. Its first d-1 power traces suffice,and it has a pole on every omitted normalization divisor under nonautomorphy. Even a nonnormal monic energy order would close the map. The universal HC4 application uses retained Tang;no automatic integrality is proved.
- **Whole-source power traces give a concrete global inverse gate:** [HC4MYKO0/1/2](WHOLE_SOURCE_POWER_TRACES_AND_DIFFERENTIAL_SPLITTING.md) identifies normalization integrality from the first minimal-degree-minus-one power traces. All source generators,or one generic boundary-detecting linear form,suffice for an inverse. A differential retraction is uniquely normalized field trace. Trace on the finite normalization is insufficient:an actual three-dimensional Keller control has a source trace pole. The missing regularity needs extra plane/gradient structure;HC4 remains open.
- **Global canonical data now exclude every stable complete-intersection normalization:** [HC4MYKN0/1/2](CANONICAL_BOUNDARY_AND_STABLE_COMPLETE_INTERSECTION_OBSTRUCTION.md) proves that every nontrivial embedded Keller intermediate normalization has nontorsion Weil canonical class,in all dimensions and degrees. Thus no global CI presentation works with any number of helpers,even after trivial polynomial/Laurent stabilization. This gives a direct HC4 inverse criterion without Tang. Local Gorensteinness is insufficient;the remaining normalization needs genuinely nontrivial global canonical data. The prior class-group and degree-six observations are credited,and the non-CI remainder stays open.
- **The nonplane horizontal coefficient locus now has a universal inverse:** [HC4MYKM0/1/2](COEFFICIENT_CLOSURE_COMPONENTS_AND_ALL_PRIME_BARRIERS.md) combines a finite integral inverse ideal with the projective plane incidence. Every characteristic-zero component is wholly automorphic or wholly in the plane locus. The integer horizontal closure has an actual inverse outside that locus,including over coefficient rings with nilpotents. Fixed-degree prime/precision certificates exist,but every denominator must contain the square of every prime up to D-1;the all-prime controls pass modulo p squared and still never lift. The ordinary naive ring and plane-Keller components remain open.
- **Arbitrary characteristic-zero potential lifting now reduces exactly to plane lifting:** [HC4MYKL0/1](INTEGRAL_COTANGENT_DESCENT_AND_PLANE_LIFT_EQUIVALENCE.md) combines retained Tang existence,geometric plane uniqueness and a saturated DVR lattice. The plane descends to the coefficient field and has an integral GL4 frame over the same DVR. Thus a noninvertible residue potential lifts iff some residue linear cotangent frame has an actually liftable plane Keller pair;with a unique standard residue plane,the original pair must lift. This does not apply to finite precision or decide HC4/JC2.
- **Nonlinear Frobenius shears now fail a global lifting test:** [HC4MYKK0/1/2](SPECIALIZATION_OF_AFFINE_PLANES_AND_FROBENIUS_SHEARS.md) shows that every noninvertible HC4 good reduction must retain a geometric jointly affine two-plane,using the retained Tang alternative and saturated lattices. Arbitrary polynomial square shears with positive u-degree preserve the characteristic-two Hessian and collision but leave only one affine direction,so none lifts in any degree or ramification. One shear passes the old divided equation,and all have mod4 lifts;those tests are insufficient. General nonlinear source orbits and HC4 remain open.
- **A full-HC4 arithmetic counterexample route is closed in every degree:** [HC4MYKJ0/1/2](RAMIFICATION_INDEPENDENT_LIFTING_OBSTRUCTIONS_FOR_HC4.md) proves sharp lifting thresholds for the characteristic-two plane example and its cotangent potential over every ramified DVR. The plane stops at precision e;arbitrary four-variable potential corrections stop at2e. Thus the unramified potential lifts mod4 but never mod8. Other reductions,nonlinear potential charts and HC4/JC2 remain open.
- **Actual critical type one now recovers the whole source:** [HC4MYKI0/1/2](CRITICAL_GORENSTEIN_FIBRES_AND_FIRST_ORBIT_MODULES.md) closes every h-dependent plinth family locally a complete intersection at its labelled isolated critical values. Gorenstein actual fibres or at most3 generators of ker(delta^2) suffice. Remaining fat components need a genuine critical non-Gorenstein type jump and length at least3;total length is N or at least N+2. The desired rank-two first orbit module itself is not free at a Morse quotient point. Automatic critical type one,new native actions and full HC4 remain open.
- **Weighted periods now recover the entire candidate plinth scheme:** [HC4MYKH0/1/2](WEIGHTED_PERIOD_RECOVERY_AND_THE_GLOBAL_NILPOTENCE_GAP.md) identifies generic Hankel rank with actual scheme length. The leading normalized Morse residue is automatically reduced;its specialization loses precisely the actual nilradical. Saturating that residue embedding would force minimum length,but remains open. Arbitrary-thickness local Morse-volume examples have even compatible nilpotent additive jets while their normalized field is not globally locally nilpotent. This targets the global algebraic action and does not reopen the exhausted ordinary volume periods.
- **Both entire three-point standard quotient classes are excluded:** [HC4MYKG0/1/2](THREE_POINT_TRANSLATION_QUOTIENT_BARRIER_AND_MORSE_COUNTERCONTROL.md) finds fixed affine lines on every regular S_rho fibre for all standard-module coefficient degrees,including collisions. Primitive normalization extends this to every action with the same embedded invariant field,all replicas and h-preserving conjugates. Even pure source alignment requires a different embedded field,independently of E3. An actual primitive tangent-identity quadric action also rules out a blanket Morse-first-jet shortcut;its horizontal plinth component leaves the KF fixed-cylinder gate open.
- **One actual plinth-scheme value now certifies the full source frame:** [HC4MYKF0/1/2/3](FLAT_PLINTH_FAMILIES_AND_ONE_VALUE_RECOGNITION.md) proves fixed-support plinth schemes and their nilradicals finite free over the value line. Native H3 fixes the support count;one genuine ideal fibre of length N forces the whole scheme reduced and recovers the cotangent inverse. The action is a constant primitive Jacobian derivation. Minimum length and the native action remain open;specialized class annihilators can give false positives. A focused line-bundle manuscript audit finds a direct codimension mismatch,without importing its cancellation or coordinate claims.
- **Radicality is now forced by the native finite-cycle basis:** [HC4MYKE0/1/2](FINITE_PLINTH_DESCENT_AND_VANISHING_CYCLE_RECOGNITION.md) weakens the action gate to exact embedded kernel B[h] and exact finite-colength plinth descent I B[h]. A residue pairing gives two local generators without assuming them;every nonreduced support point leaves a primitive generic-fibre cycle outside the finite Morse span. KA exhaustion forces I radical,and KD recovers the full cotangent inverse. Exact ideal descent and the individual native action remain open;E_N is separate.
- **Exact kernel and radical plinth now recover the entire source frame:** [HC4MYKD0/1/2](RADICAL_PLINTH_RIGIDITY_AND_AUTOMATIC_QUOTIENT_RECOVERY.md) removes the source-codimension and one-reached-value hypotheses from JV recovery on A4. Kernel factorial closure supplies codimension;integral H3 detects empty plinth lines and H4 detects extra reached values,so source acyclicity forces exactly one value on every line. One action with exact embedded kernel B[h] and exact radical finite-point plinth now suffices for pure cotangent source alignment. Constructing that action and E2/E3 remain separate open gates.
- **Global polynomial lifts now settle the collective additive geometry:** [HC4MYKC0/1/2](GLOBAL_FINITE_DIFFERENCE_LIFTS_AND_RELATIVE_FLEXIBILITY.md) constructs actual finite-difference shear lifts in every weighted-model degree. Every simple cotangent potential has complete relative ML=C[h],DK=R and a transitive,flexible group of global additive actions on each smooth fibre locus. All native frontiers have the same collective orbit pattern,including collided critical values. One required free action and its precise quotient,source alignment and E2/E3 remain open.
- **Global normal forms now separate the remaining alignment and Keller gates:** [HC4MYKB0/1/2](GLOBAL_ORIENTED_NORMAL_FORMS_AND_SEIFERT_BARRIER.md) reduces every admitted simple cotangent potential with common squared critical orientation to explicit R_(N,c). For S_rho,unit existence is exactly E3 AND an actual source equivalence to the written R_rho;K2 unit would force K2~h1 as well as E2. The nonproper Seifert theorem also forces signed identity homological form and zero positive infinity-relative groups. Actual source alignments,regular links and Milnor-piece embeddings remain open.
- **The complete integral global volume-period test now passes:** [HC4MYKA0/1/2](GLOBAL_INTEGRAL_VANISHING_CYCLES_AND_VOLUME_PERIODS.md) identifies all local vanishing cycles as a unimodular basis of the entire middle homology,including every S_rho collision. Monodromy is identity and every period vector is a common scalar times (s-c_i),exactly the required cotangent pattern. This retires further volume-period calculations within that object;Seifert/thimble/boundary markings and actual critical scheme gluing remain open.
- **The full smooth-family motivic test also passes:** [HC4MYJZ0/1/2](SMOOTH_VALUE_FAMILIES_AND_THE_MOTIVIC_HOMOTOPY_BARRIER.md) gives every normalized smooth value family the same relative unstable A1-type as its required cotangent torsor over an explicit dense value open. The type is a constant wedge of N copies of P1 smash Gm,including the moving S_rho puncture. This retires further unmarked homotopy refinements on that open;actual actions,marked quotient data,critical gluing and source equivalence remain distinct.
- **Every visible coordinate route is now closed for K2 and S_rho:** [HC4MYJY0/1/2](TORUS_CRITICAL_CURVES_AND_ALL_COORDINATE_VISIBILITY.md) excludes a unit source whenever any native K2 coordinate pullback is affine,or any coordinate pullback in the actual S_rho suspension frame is affine. A whole torus Keller fibre with two zeros would force forbidden generic degree two;the S_rho v route instead fails a uniform relative determinant identity. No arbitrary affine alignment or full-source exclusion is asserted.
- **The universal affine frontier now contracts to four requirements:** [HC4MYJX0/1/2](GLOBAL_SIGN_TWISTS_AND_TWO_POINT_KELLER_EQUIVALENCE.md) makes every h_D=2xy+(z^2-1)w+D(z) unit-existence question exactly the existence of an actual two-point plane Keller fibre. An explicit z-preserving source automorphism flips the root orientations and identifies h1 with C2. Thus h0/h1 decisions coincide,and K2 is redundant in the universal test list. What remains is the two old endpoint inversion requirements,absence of a two-point Keller fibre,and the whole S_rho family. Individual K2 source equivalence and all these decisions remain open.
- **The K2 one-coordinate route is now closed in every degree:** [HC4MYJW0/1/2](GENERIC_COORDINATE_SPECIALIZATION_AND_K2_ONE_COORDINATE_CLOSURE.md) excludes every source change with either native x or y pullback affine,even with arbitrary mixing of the other variables. Generic-coordinate specialization rules out all hidden exceptional poles;the complete plane pencil then forbids the required trace coordinate. Any remaining full-source candidate must change both pullbacks nonlinearly. This does not classify every action fixing x or settle the full K2 orbit.
- **Exact action data now recover the global source coordinates:** [HC4MYJV0/1/2](TWO_SLICE_RECOVERY_AND_ORIENTED_COEFFICIENT_MODELS.md) reconstructs the complete polynomial inverse from two global slices and an affine-socle congruence. Reduced plinth holes also recover a full coordinate plane bundle under the stated codimension hypothesis. The remaining fixed-base condition is globally Keller generators in the correct unsquared coefficient orientation; arbitrary nonlinear fibre maps add no further pairs. Explicit all-N auxiliary models encode Keller fibre cardinalities,but no action,pair,N or original graph realization is constructed.
- **The quotient gate is now exact,and all local volume periods are exhausted:** [HC4MYJU0/1/2](QUOTIENT_HOLES_TORSOR_RESIDUES_AND_LOCAL_PERIOD_BARRIER.md) proves the exact radical plinth ideal,its N missing quotient lines with one reached fixed-plane point each,and the whole-value torsor annihilators and affine residues with one common orientation. Every remaining normalized model passes the complete local volume-period test to all orders,including colliding values. A different quotient remains necessary;JV recovers the full inverse under its exact global class hypotheses. No action or source certificate is constructed.
- **The action gate now excludes whole standard quotient fields:** [HC4MYJT0/1/2](PRIMITIVE_QUOTIENTS_PLINTH_AND_GLOBAL_ACTION_CERTIFICATES.md) proves global primitivity and rules out every action sharing a standard K2 embedded invariant field,including invariant replicas and value-preserving source conjugates. A candidate needs a different quotient and a proper nonprincipal plinth ideal. An exact global coordinate-action certificate specifies the complete source inverse and constant coefficient Jacobian needed to produce unit Hessian; no required action or certificate has been constructed.
- **A global additive-action gate replaces the now-uninformative kernel invariants:** [HC4MYJS0/1](GLOBAL_ADDITIVE_ACTIONS_AND_THE_K2_TRANSLATION_BARRIER.md) computes relative ML=C[h] and relative DK=the whole source ring for every current test,using an exact uniform K2 value suspension. A hypothetical unit-Hessian member needs an action free on every noncritical value fibre,with a polynomial kernel and fixed affine planes. Both complete standard K2 translation modules and all value-preserving source conjugates have generic fixed points and are excluded in every degree. Other actions and actual source recovery remain open.
- **Keeping one hyperbolic coordinate is already excluded:** [HC4MYJR0/1](ONE_HYPERBOLIC_COORDINATE_AND_INTRINSIC_CYLINDER_PLANE.md) recovers the exceptional trace plane intrinsically and forces one global free-coordinate scale. It excludes all C[x]- or C[y]-preserving source maps for h0,h1 and all C[P]/C[V]-preserving maps in the exact three-point suspension representative. A constant first normal slope gives the same obstruction for every higher normal jet. K2/h1 equivalences retaining x or y are excluded;JW2 now also closes the K2 one-coordinate unit-Hessian route. Unrestricted source changes remain open.
- **The additional affine tests now normalize to one fixed function and one one-parameter family:** [HC4MYJQ0/1](TWO_NORMAL_MODELS_AND_RELATIVE_RANK_ONE_ALIGNMENT.md) identifies the scalar and length-two cubic families by a global q-preserving shear. All affine outer degrees reduce to the four earlier tests plus K2=w(2xy+z^2-1)+z and K3_rho=K2+2x(rho-z). K3 has an explicit hyperbolic suspension form. Every source map preserving their respective natural bases C[x,w] and C[w,y] is excluded; arbitrary full source changes remain open.
- **All affine quadratic composites now reduce to the four old tests and three explicit resonant cubic families:** [HC4MYJP0/1/2](ALL_AFFINE_PENCILS_TO_THREE_RESONANT_CUBIC_FAMILIES.md) proves positive mixed-plane genus in every higher outer degree and closes every nondegenerate ternary-cone multiplier. Beyond the old coordinate/quadratic and h0,h1 tests, only F_N=q+alpha x+e z+t(lambda q+N-c), q=2xy+z^2, N=0,x^2,2xz, lambda*c*e!=0 remains. These three open full-source families pass the whole-value and common-critical-Hessian controls; no equivalence with the old tests is assumed.
- **Every isotropic-kernel pencil and the entire cubic coordinate-inner family now reduce to the same four fixed tests:** [HC4MYJO0/1](ISOTROPIC_KERNELS_AND_SIMPLE_ZERO_QUARTIC_FRONTIER.md) covers even nonorthogonal offsets, using positive all-degree trigonal genus and an exact binary-cone resonance calculation. Beyond the old coordinate/quadratic orbits and h0,h1, a further affine mechanism must have a simple nondegenerate zero eigenline with nonzero offset and the single quartic form with two nondegenerate ternary quadratics. That quartic family and all four fixed tests remain open.
- **Every inner quadratic with a critical point now reduces to the same four fixed tests:** [HC4MYJN0/1](COMPATIBLE_INNER_CENTERS_AND_QUADRATIC_COORDINATE_FRONTIER.md) proves the needed specific quadratic stabilizations and closes the entire compatible-center and admissible null-kernel branches to the old coordinate/quadratic tests and h0,h1. Beyond those tests, a further affine mechanism must have a genuine degree-two coordinate inner, singular rank-two/rank-three T and no admissible null direction. Explicit cubic/quartic-plus-univariate source forms describe that remaining frontier; they are not asserted solved.
- **Every rank-at-most-one affine pencil now reduces to the same four fixed tests:** [HC4MYJM0/1](ZERO_EIGENVECTOR_SHEARS_AND_RANK_ONE_PENCIL_CLASSIFICATION.md) gives a global zero-eigenvector shear and stabilizes the complete critical-mass obstruction. For a nonzero-eigenvalue rank-one pencil every source orbit is excluded except one invariant quadratic resonance, giving exactly h0 or h1. Together with nilpotent reduction, all low-rank graph/source degrees reduce to the old coordinate/quadratic tests and those two full-source models. Beyond them, additional affine mechanisms must have singular rank-two or rank-three T.
- **Every invertible affine quadratic pencil is now inverse-closed:** [HC4MYJL0/1](JORDAN_POLAR_CURVES_AND_ALL_INVERTIBLE_PENCIL_CLOSURE.md) removes every Jordan, offset-cyclicity and Krylov-pairing qualification in four variables, for all nonlinear outer/source degrees. Plane Jordan polar curves, one complete eigenvector fibre surgery and exact low-genus stationary counts close all remaining nonsemisimple invertible cases. The conclusion is gradient inverse under unit Hessian. Singular pencils, h0,h1 and both actual old endpoints remain open.
- **Suspended univariate multipliers now reduce to four fixed source-orbit tests:** [HC4MYJK0/1](SUSPENDED_MULTIPLIER_ORBITS_AND_BASE_PRESERVING_CLOSURE.md) handles arbitrary C,D in 2xy+w C(z)+D(z), and every compensated graph X=w+A(z). A derivative rigidity excludes squarefree degrees at least three; repeated roots are excluded or inverse-closed. Beyond the two old endpoint orbits only h0=2xy+w(z^2-1) and h1=h0+z remain. Neither can acquire unit Hessian through any source automorphism preserving C[x,z], even with nonlinear base and dual maps. Their entire special-value Hodge profiles pass; their full source orbits remain open.
- **Every invertible semisimple affine pencil is now inverse-closed:** [HC4MYJI0/1](AFFINE_POLAR_CURVES_AND_COMPLETE_INVERTIBLE_SEMISIMPLE_CLOSURE.md) handles every affine offset, isotropic projected offset and repeated eigenvalue in all nonlinear outer/source degrees. An all-offset polar-curve genus calculation and complete Witt/hyperbolic fibre strata close the cases left by the residue bounds. The conclusion is gradient inverse under unit Hessian. Singular and nonsemisimple pencils remain separate; a genuine singular rank-one two-critical-point model passes the present critical/Hodge tests.
- **Simple-spectrum affine pencils now exclude every offset, including noncyclic ones:** [HC4MYJH0/1](SIMPLE_STATIONARY_POLES_AND_NONCYCLIC_KRYLOV_CLOSURE.md) counts the complete finite critical group at every compatible simple determinant pole. The resulting zero reciprocal-Hessian sum excludes unit Hessian throughout every source orbit for invertible simple-spectrum pencils, in all nonlinear outer and admitted graph degrees. A broader orthogonal-Krylov certificate retains the exact active multiplicity budget and inactive-primary hypotheses. Degenerate Krylov spaces and unsupported multiple boundary poles remain separate.
- **Every cyclic affine quadratic-composite graph now closes in all source degrees:** [HC4MYJG0/1](ELLIPTIC_CONE_TRACES_AND_COMPLETE_CYCLIC_CLOSURE.md) resolves the former quadratic infinity resonance by a literal polynomial shear to a genus-one cone trace. Cyclicity prevents singularity at the tangential infinity points; four discriminant branch points give the whole-fibre Hodge obstruction. Together with JF1 this forces gradient inverse under unit Hessian for every admitted cyclic pencil and nonconstant outer graph, including certified quadratic-similitude images. Nonnilpotent noncyclic offsets and the actual old endpoints remain outside this closure.
- **Independent centers now have a global reciprocal-Hessian pencil test:** [HC4MYJF0/1/2](AFFINE_QUADRATIC_PENCIL_RESIDUES_AND_ENDPOINTS.md) excludes entire source orbits when all finite pencil boundaries and infinity are accounted for. Invertible cyclic offsets cover every nonlinear outer degree; all cyclic graph degrees at least two also close. Every nilpotent affine-offset graph family reduces exactly to the two old plane endpoint alternatives for universal inversion. Its former cyclic quadratic double-infinity resonance is now inverse-closed by JG0/1 above. Noncyclic nonnilpotent offsets and both actual endpoints remain open.
- **Every commonly centered quadratic-composite graph image now returns to the old endpoint:** [HC4MYJE0/1](ALL_CENTERED_QUADRATIC_COMPOSITE_IMAGES.md) closes every remaining nonnilpotent Jordan type to gradient inverse under unit Hessian, in arbitrary outer and source degrees. Exact whole-fibre Hodge profiles keep the hyperelliptic/biquadratic genera, all branch-point omissions and the genus-zero degree-two cases. Together with the nilpotent normal form, the entire admitted graph family and its certified quadratic-similitude saturation has precisely the old singleton-zero endpoint; no spectral branch remains inside this family. Independent centers and unrestricted embeddings remain outside it.
- **Nilpotent quadratic-composite images give exactly the old quadratic endpoint:** [HC4MYJD0](NILPOTENT_QUADRATIC_PENCILS_AND_ENDPOINT_ORBITS.md) gives explicit triangular shears for every four-variable nilpotent pencil and arbitrary outer polynomial. Each fixed graph image, over all source parameterizations and certified quadratic similitudes, has precisely the nondegenerate quadratic/singleton-zero orbit. Its former nonzero-Jordan frontier is now inverse-closed by JE0/1 above. The endpoint itself remains open.
- **Semisimple quadratic-composite images are excluded in every source orbit:** [HC4MYJC0/1](SEMISIMPLE_QUADRATIC_PENCILS_AND_ALL_SOURCE_ORBITS.md) proves that q_A+H(q_B), with a common center, nonzero diagonalizable A^(-1)B and nonlinear H, cannot acquire unit Hessian under any polynomial source automorphism. Hence every parameterization of X=P0(q_B), with nonconstant P0, and its entire quadratic-similitude image saturation are excluded in arbitrary degrees. The critical Hessian spectrum supplies the obstruction. Nonsemisimple pencils are essential exceptions: a nonzero nilpotent example has unit Hessian in all outer degrees.
- **All divisor inverse jets are already automatic:** [HC4MYJB0](FORMAL_DIVISOR_RECOVERY_AND_RATIONALITY_GATE.md) reproduces the standard formal-etale thickening theorem here: B/(a^m)=R/(f^m) for every m. Fixed rational recovery or uniformly bounded rational target degrees force the Keller inverse; hypothetical failure requires linearly growing degree for at least one coordinate. A non-Keller etale-open model has every jet and places its trace pole at the finite-image escape. Longer normal-jet searches do not close the global embedding gap.
- **Extra source divisors now have explicit unit and boundary obligations:** [HC4MYJA0/1/2](EXTRA_PRIME_EXCISION_AND_ONE_PLACE_BUDGET.md) gives the whole image equation a(K)=f*g with (f,g)=1, including multiplicities. A selected one-place radical component forces a disjoint component with several punctures and a nonconstant source unit; scalar units on every extra component force the inverse. For polynomial normalization there is at most one one-place extra, only in the all-nodal case, and a finite-image boundary remains compulsory. The full w-dependent recovery obligation remains open.
- **Every descended rank-one trace now closes:** [HC4MYIZ0/1/2](ONE_PLACE_PENCILS_AND_DESCENDED_TRACE_CLOSURE.md) removes the nodal and other multibranch exceptions. Explicit one-place compactification and proper-flat normalization inputs give the intrinsic generic genus G=g0+delta with one puncture, so the complete curve-isomorphism criterion forces degree one. Any positive power of the trace in the Keller target field suffices, with arbitrary helper degrees and the stated affine offset frames. Remaining one-place radical components require extra source primes outside the radical curve; genuinely non-descended traces and higher helper ranks remain open.
- **Earlier contractible-cusp and pencil-defect gates:** [HC4MYIY0/1/2](CONTRACTIBLE_RADICAL_CURVES_AND_PENCIL_DEFECT.md) retains the exact generic-pencil boundary identity and the Lin-Zaidenberg contractible-curve proof. Its former internally multibranch descended frontier is closed by IZ2 above; the retained proof bytes are unchanged.
- **Nonconstant rank-one traces now have a closed radical-curve test:** [HC4MYIX0/1/2](UNBOUNDED_QUADRIC_CLASSES_AND_RADICAL_TRACE.md) makes K a closed embedding on the normal-trace zero curve. Any actual affine-line component closes, including every coordinate-composite normal trace. If the square trace descends to the target, every component must be singular with affine-line normalization under hypothetical failure; smooth or nonparametric components force the inverse. Separately, explicit genuine affine planes realize every integer quadric middle class, so arbitrary planes cannot replace the special IW sections.
- **Every constant-trace intrinsic free projection now closes:** [HC4MYIW0/1/2](QUADRIC_MIDDLE_CYCLE_AND_FREE_PROJECTION_CLOSURE.md) closes the full noncritical rank-one chart with arbitrary polynomial helpers. Two principal-divisor cycles give d copies of one infinite-order integral middle class equal to one copy, forcing the Keller degree d=1. Explicit proper homotopies and a compact real sphere supply the class calculation. Together with IU0/1 this closes every constant trace, including the stated constant affine frames and polynomial offsets. Aligned projection existence and nonconstant traces remain open.
- **Primitive helper recovery must work on source divisors:** [HC4MYIV0/1/2](PRIMITIVE_HELPERS_AND_LEGENDRE_NORMALITY.md) gives the complete noncritical compensation and its gauge freedom, then strengthens the one-helper reconstruction gate to codimension one. For a primitive helper the implicit derivative generates the pulled-back singular ideal; normality, target descent or suitable componentwise saturation forces the Keller inverse. The canonical Legendre energy is primitive, and its hypersurface is normal exactly when the gradient is invertible. Hypothetical cotangent failure already gives collision lines on generic stable source planes. This different image must not inherit normality from the actual smooth embedding. The characteristic-zero weak-type-one mechanism and energy primitivity are explicitly attributed and reproduced.
- **Earlier critical-value and universal-compensation results:** [HC4MYIU0/1/2](FREE_PROJECTIONS_AND_CONSTANT_BASE_TRACE.md) supplies the zero-section rank classification, the critical-value inverse theorem and the universal-compensation immersion obstruction. Its formerly open full noncritical chart is closed by IW1/2 above; the retained proof bytes are unchanged.
- **Intrinsic Nambu and inverse-Hessian directions give new all-degree certificates:** [HC4MYIT0/1/2](INTRINSIC_NAMBU_DIRECTION_RIGIDITY.md) proves that three independent directions with finite actual Legendre-energy orbit force the gradient inverse. One finite direction with nonzero inverse-Hessian square also suffices. Under hypothetical nonautomorphy the finite-direction space and the exact constant inverse-Hessian zero cone are the same linear two-plane; any other cone shape closes the candidate. Zero-square directions give explicit quasi-translations, and pure whole Hodge traces can be mixed with finite-energy directions in a rank-three certificate. No ambient observable or generic-class condition is required. A supplied zero-square two-plane directly reconstructs the actual cotangent coefficient pair by block inversion; a plane-shaped cone alone remains inconclusive, and no extra finite direction is automatic.
- **Dropping the residual direction already contains all of JC2:** [HC4MYIS0/1/2](UNIVERSAL_AFFINE_MODIFICATION_BOUNDARY.md) puts every plane Keller pair in a canonical closed HC4-size graph with primitive multiplier content, smooth geometrically integral observable fibres everywhere and scalar generic units. Its exact class lattice is Z^d modulo the diagonal; the selected source-fixed class vanishes exactly when the original Keller field degree d is one. The literal graph has actual unit Hessian and transfers to the fixed Meng--Yang potential, but its constructed gradient always has an explicit inverse. Thus an arbitrary carried pair must be distinguished from intrinsic nonautomorphic HC4 data. A residual-free universal carried-inverse extension is equivalent to JC2, not a smaller smoothness or coefficient gap.
- **The complete admitted HC4 overlap family now closes:** [HC4MYIR0/1/2](TWO_HELPER_CURVE_PENCIL_CLOSURE.md) closes every two-helper polynomial-direction overlap in every stable dimension and degree, and together with the one-helper theorem every HC4-size overlap with the residual-direction and exact closed-content conditions. The generic helper map is proved flat; each actual source component surjects onto its helper curve by a finite-line argument. Their disjointness and scalar units give a helper-plane pencil of degree [L:F], forced to one by IQ0. This needs no generic ambient filtration, primitive content, curve/action classification or forgetful faithfulness. Multiplier-zero images close in every helper count under the generic residual condition alone. Residual/content failures and compatible arbitrary-graph observables remain outside this theorem.
- **Every one-multiplier-variable overlap now closes:** [HC4MYIQ0/1](ONE_HELPER_OVERLAP_CLOSURE.md) covers (T,M(T,X)Y-P(T,X,Z)) with one X, arbitrary root-variable count and degrees, a surviving polynomial direction on every generic multiplier row, and exact nonconstant closed content levels. Rank-two descent linearizes the source helper; elimination or the finite self-tensor argument then forces the original Keller inverse. Nonprimitive content levels depending only on X are admitted. No generic filtration, monic/Bezout leading term or forgetful-flatness assumption is required.
- **The earlier two-helper reduction is now closed:** [HC4MYIQ2](ONE_HELPER_OVERLAP_CLOSURE.md#3-the-remaining-hc4-overlap-has-two-helpers-and-one-root-variable-hc4myiq2) is retained as the intermediate invariant-curve reduction used by the complete-family theorem above. It is not a remaining open case within the admitted family or a universal HC4 graph normal form.
- **A moving additive action replaces the generic positive filtration:** [HC4MYIP0/1](ADDITIVE_ACTIONS_AND_IMAGE_INVARIANCE.md) forces the carried inverse when a relative locally nilpotent derivation moves the image equation, under faithful flatness, integral height-one fibres and principal normal generic components. Both affine-space quotient unit groups make the component relations scalar, even with a large degree-zero kernel. Generic scalar units are deduced. Joint constants F, including trivial relative Makar-Limanov invariant, suffice; no generic or global positive filtration is required.
- **Every surviving polynomial-direction overlap image is additive-invariant:** [HC4MYIP2](ADDITIVE_ACTIONS_AND_IMAGE_INVARIANCE.md#4-every-surviving-overlapping-polynomial-direction-image-is-vertical-hc4myip2) puts its generic equation in F[X], and primitive multiplier content puts the actual global equation in C[T,V,X]. This covers all overlaps satisfying the residual-direction and exact content conditions, even when no generic filtration is known. An equation outside that invariant subring forces the inverse. All one- and two-helper stationary cases under these hypotheses are closed above; higher-helper cases in larger stable dimensions remain separate.
- **Independent coefficient fields and overlapping root coefficients are now covered:** [HC4MYIO0/1](UNIVERSAL_COEFFICIENT_FIELDS_AND_OVERLAPPING_OUTPUTS.md) permits polynomial coefficients in any independent extension D/C(T), including multiplier-component function fields. Algebraic coefficients themselves are handled before guarded finite-type specialization. A surviving polynomial Z direction on every geometric multiplier component then kills all generic source-fixed classes for arbitrary overlapping P(T,X,Z); row component counts may differ. Coefficients identified with the original source g are excluded.
- **Monic and Bezout leading overlap outputs force the carried inverse:** [HC4MYIO2](UNIVERSAL_COEFFICIENT_FIELDS_AND_OVERLAPPING_OUTPUTS.md#3-a-carried-inverse-and-uniform-generic-filtration-certificates-hc4myio2) combines exact nonconstant closed content levels with a positive generic filtration. Every monic root term with arbitrary overlapping lower coefficients, and P=a(T,X)b(T,Z)+Q with u A+v a=1 over C(T)[X] and deg_Z Q<deg_Z b, has a uniform generic filtration. All degrees and admitted stable dimensions are unrestricted; no global free frame is required. The stronger complete-family theorem above closes the one- and two-helper invariant cases without that filtration.
- **Plane outer irreducibility survives every algebraic parameter extension:** [HC4MYIM0/1](PARAMETER_STABILITY_AND_MOVING_MULTIPLIERS.md) keeps P(f,Z)+g irreducible over L E for every finite E/C(T), equivalently over L kappa_bar, with arbitrarily many polynomial Z variables. Bad horizontal inertia lifts in a Galois tower and stays in the original plane source subgroup while fixing E. Splitting moving multiplier factors over E therefore kills every generic correlated row/root class invariant for disjoint A(T,X)Y-P(T,Z). Full splitting-group preservation is not asserted.
- **Every primitive disjoint moving multiplier now closes without a global frame:** [HC4MYIN0/1](REFLEXIVE_PENCIL_DESCENT.md) covers (T,A(T,X)Y-P(T,Z)) for arbitrary primitive A and arbitrary positive-joint-degree P in disjoint variables. Nonconstant multiplier content is also allowed when each P(c,Z) at its zeros is nonconstant closed. A relative derivation bounds a rank-two pencil lattice; reflexivity on the plane target and faithful flatness supply its global primitive coordinate and saturation. All degrees and admitted stable dimensions are unrestricted. Only a generic filtration is needed. Nonconstant composite content levels remain separate; the next theorem above admits overlapping coefficients under its explicit hypotheses.
- **Earlier global frames remain retained certificates:** [the exact content criterion and leading frames](PARAMETER_STABILITY_AND_MOVING_MULTIPLIERS.md#3-separate-the-global-filtration-from-the-content-fibres-hc4myim2) and [the parameter-unimodular degree-gap frames](PARAMETER_STABILITY_AND_MOVING_MULTIPLIERS.md#4-a-uniform-adapted-frame-beyond-both-original-leading-certificates-hc4myim3) keep their proofs and controls. HC4MYIN0/1 removes the need to find such a frame in the disjoint exact-content family. A constant content level instead forbids a carried Keller pair immediately.
- **Arbitrary fixed multipliers and multivariable suspensions now close:** [HC4MYIL0/1/2](MULTIVARIABLE_SUSPENSIONS_AND_FIXED_MULTIPLIERS.md) proves the carried plane inverse for every (T,A(X)Y-P(T,Z)), with A in C[X] nonzero and independent of T,Z, positive joint Z degree, p,k>=1 and every stable r>=p+k-1. All degrees, variable counts and fixed-factor multiplicities are unrestricted, including ambient/target polynomial automorphism orbits. A generative polynomial transfers outer irreducibility to all geometric factors; each distinct multiplier factor contributes one exact root-lattice copy. Its multiplicity has no effect. The original global free filtration controls every height-one target fibre. HC4 includes arbitrary P(T,Z_1,Z_2) with A=X and arbitrary fixed A(X_1,X_2) with one Z. The moving case is now covered by the exact-content and rank-two descent theorem above; overlapping variables remain separate.
- **Every polynomial triangular output now closes:** [HC4MYIK3](PARAMETER_NORMS_AND_ALL_POLYNOMIAL_OUTPUTS.md#4-every-polynomial-output-in-the-plane-hc4myik3) proves P(f,Z)+g irreducible for every nonconstant P in C(T)[z] and every plane Keller pair. Thus every carried closed stable embedding through (T,xy-P(T,z)), P polynomial of positive z degree and r>=1, forces the Keller inverse. All coefficient, leading-multiplicity, source, embedding and z degrees are unrestricted. Bad inertia supplies polynomial blocks; parametrizable outer critical curves give integral critical values, and parameter norms force block transitivity. No further plane coefficient or factor-degree subdivision is needed. Full splitting-group preservation and a higher-dimensional unrestricted theorem are not claimed.
- **Parameter-cover unit norms replace constant-unit shortcuts:** [HC4MYIK0/1/2](PARAMETER_NORMS_AND_ALL_POLYNOMIAL_OUTPUTS.md) preserves normal monogenic irreducibility over finite algebraic parameter curves in every Keller dimension. A source unit need not be constant, but its finite-free norm is. Integrality of all conjugate monic models makes its target inverse integral, so inertia and a preserved parameter quotient exclude a factor. A direct valuation proof makes integral critical values sufficient for an integral centered monic model. This also completes the no-finite-pole side of every plane binomial, while retaining the earlier stronger group theorem on the finite-pole side.
- **Full-group and higher-dimensional certificates retain their stronger scope:** [relative central kernels and product covers](CENTRAL_GEOMETRIC_QUOTIENTS_AND_BINOMIAL_OUTPUTS.md), [plane rational affine parameters and bad inertia](BRANCH_CURVES_AND_RATIONAL_TRIANGULAR_PARAMETERS.md), [strict projective multiplicity sums](MULTIPLICITY_BOUND_AND_CRITICAL_EXPANSIONS.md), [reduced rational pencils](REDUCED_PENCILS_AND_RATIONAL_MONODROMY.md), [nilpotent target covers](NILPOTENT_COVER_MONODROMY.md) and [regular parameter covers](REGULAR_PARAMETER_MONODROMY.md) remain valid. Their equality, product, multiplicity and dimensional hypotheses are not removed by the plane irreducibility theorem.
- **Earlier coefficient refinements are retained evidence:** [factor resultants](FACTOR_RESULTANT_BARRIER.md), [derivative modifications and normal rescaling](DERIVATIVE_POWER_MODIFICATION.md), [critical supports](CRITICAL_FACTOR_SUPPORT.md), [leading valuations](LEADING_VALUATIONS_AND_CRITICAL_CURVES.md), [leading-fibre units and the global filtered proof](LEADING_FIBRE_UNIT_GATE.md), and [normal monogenic covers](NORMAL_MONOGENIC_COVERS.md) keep their proofs and exact receipts. Their former plane coefficient residuals are superseded by HC4MYIK3; their all-dimensional results remain useful.
- **Normal quotients and missing-sheet classes constrain other routes:** [HC4MYIB0/1/2](NORMAL_QUOTIENT_AND_BOUNDARY.md) forbids nontrivial Galois intermediate fields in every dimension and identifies the plane normalization class group with the free group on missing divisors. Its subnormal class-stabilizer criterion applies beyond nilpotent groups. It does not preserve every arbitrary target cover; nonnormal intermediate fields and other ambient class modules remain distinct.
- **Abelian divisor-class monodromy closes nonfactorial output fibres:** [HC4MYHW0/1/2/3](ABELIAN_CLASSGROUP_OUTPUT_OBSTRUCTION.md) proves that Keller function fields have no nontrivial abelian Galois intermediate field. This makes conjugate image components principal when the generic class group is torsion-free with abelian action; geometrically integral affine height-one fibres suffice for global denominators. The obstruction now covers nondegenerate linear/quadratic output pairs for every r>=1, specified degenerate ranks, and cyclic Danielewski pairs (T,xy-z^m) in every m>=2. The latter have nonfactorial generic fibres outside every older smooth-infinity frame in HC4 size. Nonabelian class actions, torsion and nonintegral special fibres remain distinct boundaries.
- **Ambient output observables may have arbitrary degrees:** [HC4MYHV0/1](SMOOTH_INFINITY_OUTPUT_OBSTRUCTION.md) replaces the linear/quadratic restriction by a smooth scheme-theoretic complete intersection of the two highest forms at infinity. A common filtered basis controls global denominators, and an explicit classical Picard import gives generic factoriality; every closed stable embedding carrying a plane Keller pair then forces that pair invertible. Every ambient extension of a hypothetical nonautomorphic intrinsic pair must fail this infinity condition. Generic factoriality alone, or smooth reduced support, misses the written denominator trap; arbitrary extensions need not admit the required leading forms.
- **Positive-weight invariant images have explicit quadratic rectifications:** [HC4MYHU0/1](POSITIVE_WEIGHT_IMAGE_RECTIFICATION.md) writes polynomial isometries and inverses for nonlinear normal graphs in isotropic-coordinate ideals, allowing arbitrary degrees and coefficients in unchanged variables. Every closed image invariant under a compatible positive linear ambient weighting either forbids unit Hessian at its isotropic tangent normal or admits such a shear to a nonisotropic hyperplane. Under arbitrary source parameterizations the surviving images retain exactly the singleton-zero endpoint. Positive linear weights in the actual source instead force a polynomial gradient inverse; image invariance alone does not give that stronger hypothesis.
- **Quadratic source planes force alignment; all fibre-affine helpers are classified:** [HC4MYHT0/1](QUADRATIC_PLANE_FIBRE_AFFINE.md) proves that under nonautomorphy any constant two-plane with uniformly vanishing cubic derivatives is already Tang's affine plane. A genuinely quadratic plane therefore closes the candidate. This classifies every two-affine embedding j=(s,t,b+A*w), with arbitrary polynomial base coefficients and a global left inverse of A: rank-three helper quadratics close outright, and lower ranks force a base-only helper and the two old endpoints. All unit-Hessian members have at most one full critical point. All three helpers may initially depend on both free variables; alignment is a conclusion.
- **One proper relative Morse component closes a candidate:** [HC4MYHS0/1/2](PROPER_MORSE_CRITICAL_BOUNDARY.md) needs properness and vertical Hessian nondegeneracy on just one nonempty relative-critical component. Each condition separately still contains exactly JC2: the universal Morse model has projection K, while an explicit unimodular-matrix construction gives every K a finite free relative projection with a monomial basis and ramification at zero. The latter projection being proper does not make K proper. Neither potential automatically supplies a compatible ambient graph embedding.
- **Base-only normal graph images are classified without cotangent alignment:** [HC4MYHR0/1](RELATIVE_CRITICAL_GRAPH_CLASSIFICATION.md) allows j=(s,t,h(s,t),H,G) with every relative polynomial automorphism (H,G). A global relative critical section or a rank-one critical curve product forces gradient invertibility. The exhaustive graph classification covers every relative quadratic rank and intrinsic-plane position: all unit-Hessian members have at most one full critical point, and nonautomorphic members are exactly the two old endpoint orbits. Native members with multiple full critical points cannot acquire unit Hessian anywhere in their relative free-variable orbits. Arbitrary mixing of all four source variables remains distinct.
- **The output-base obstruction allows a quadratic ambient observable:** [HC4MYHQ0/1](QUADRIC_OUTPUT_BASE_OBSTRUCTION.md) proves that a closed stable hypersurface embedding carrying a plane Keller pair through a linear and a nondegenerate quadratic ambient observable forces that pair invertible. Both normal types and every r>=2 are covered; generic ambient fibres may be smooth quadrics rather than affine spaces. A normal-monomial filtration makes the pencil and denominator arguments global through singular special levels. This forbids such ambient extensions of a nonautomorphic intrinsic cotangent pair and collapses direct native A/B Keller encodings to the coordinate endpoint. No quadric embedding rectification or resolution of either plane endpoint is claimed.
- **Polynomial quadratic symmetries close genuinely nonlinear ambient images:** [HC4MYHP0/1](QUADRATIC_ISOMETRY_SATURATION.md) extends the three-affine obstruction to its entire polynomial quadratic-similitude saturation in unchanged source parameters. Hyperplane-image saturations under arbitrary source parameterization still retain exactly the two old endpoints or a critical-line obstruction. Four explicit quadratic-preserving shears produce genuine embeddings of unbounded degrees with zero affine retained pullbacks, while their exact source potential remains a unit-Hessian quadratic. New mechanisms must escape those saturated classes; a nonlinear image equation alone does not suffice. No universal rectification or generation theorem for all symmetries is asserted.

- **Two-affine endpoint geometry has global recognition tests:** [HC4MYHO0/1/2](TWO_AFFINE_ENDPOINT_GEOMETRY.md) proves that two constant directional derivatives generating the unit ideal force their unique plane to be the vertical cotangent plane and give an explicit relative coordinate inverse. Cotangent-aligned embeddings with nondegenerate three-dimensional quadratic kernel are exactly the unit-ideal endpoint. Every affine-hyperplane image in normalized coordinates, under arbitrary source parameterization, is exactly a singleton-zero or unit-ideal orbit, or an isotropic-through-center critical-line case forbidding unit Hessian. Nonlinear image hypersurfaces and unaligned two-affine cases remain distinct.
- **Three affine retained normal projections already suffice:** [HC4MYHN0/1](THREE_AFFINE_PROJECTIONS.md) closes every genuine quadratic ambient embedding with three independent affine-linear source projections, allowing arbitrary polynomials in the two remaining ambient coordinates. The proof covers all three quadratic kernel ranks and both positions of the intrinsic plane, using restricted-plane generation and the affine-line Keller endpoint. A unit-Hessian example with exactly three affine projections strictly extends the former graph criterion. Two affine projections already realize both conditional plane endpoints; that class remains open.
- **The entire ambient potential has a quadratic normal form:** [HC4MYHM0/1/2](AMBIENT_QUADRIC_GRAPHS.md) gives an explicit base-preserving polynomial automorphism Phi with Psi=lambda_star+L*T^2+N*y*Z+2*N*x*W. Every polynomial graph in every constant affine direction of this nondegenerate quadratic has invertible gradient under actual unit Hessian. Pulling it back by Phi gives genuine Meng--Yang embeddings in all degrees. Four independent affine-linear retained normal projections suffice to close a candidate. The actual source parameters matter: full source reparameterizations already return both unresolved plane endpoints, and the transverse square condition alone has universal constant-D solutions.
- **The fixed potential closes the middle mixed-projection stratum:** [HC4MYHL0/1/2](FIXED_POTENTIAL_PROJECTION.md) compares the nonlinear one-polynomial boundary at physical x=y=0 with the free affine variable left by eta rank one. It forces a triangular Keller pair, without requiring embedding for that rank. Consequently a nonautomorphic fixed-potential embedding whose physical bases are centered linear gradient combinations must have eta rank two. More general centered gradient-dependent bases transfer the ambient critical value to at most one source critical point. The remaining transverse case has an explicit global Euler square condition; it and arbitrary physical-base data remain open.
- **Every mixed gradient projection has an exact generator profile:** [HC4MYHK0/1](MIXED_PROJECTION_PROFILE.md) gives minimum 4-k+ell additional generators over k gradient combinations, where ell is their intersection dimension with the intrinsic cotangent output plane. The minimum ambient dimension is 4+ell, stable under adjoining free variables. Two-output projections therefore require A4, A5 or A6 according to eta rank two, one or zero. All allowed A5 completions are explicit; they do not automatically match the fixed Meng--Yang potential. Noncoordinate directional derivatives also recognize exactly the intrinsic cotangent plane under nonautomorphy.
- **Every closed hypersurface embedding over a plane Keller base:** [HC4MYHJ0/1](STABLE_HYPERSURFACE_OBSTRUCTION.md) proves that a closed embedding (z,w) to (K(z),H(z,w)) in A5 forces K invertible, with arbitrary nonlinear H. Disjoint generic components form one polynomial pencil; primitive coefficient content removes the target denominators, and smoothness forces a linear relation. This removes every source-dependent matrix and nonlinear-dual escape over the same physical Keller base. In general, a hypothetical nonautomorphic plane Keller map needs exactly r+2 additional algebra generators after adjoining r free source variables. An HC4 cotangent gradient pair cannot extend to an ambient coordinate pair in A5. General graph embeddings need not have that output-base interpretation.
- **Global one-coordinate reconstruction and the exact universal encoding:** [HC4MYHI0/1/2/3](RECONSTRUCTION_EMBEDDING_OBSTRUCTION.md) reproduces the classical Formanek criterion and closes every closed-embedded realization of A=F, B=xF+yG for an independent plane Keller pair, in every degree. A hypothetical HC4 counterexample needs exactly two additional polynomial generators over its full gradient target. Its narrower target-dependent affine-matrix obstruction is now subsumed by HC4MYHJ1.
- **One finite-fibre template family captures the remaining cotangent problem:** [HC4MYHH0/1/2](FINITE_FIBRE_NORMAL_FORM.md) normalizes every hypothetical nonautomorphic unit-Hessian polynomial with n>=1 critical points to D(s)+A_n(s)*P+B_n(s)*U*Q, where A_n=product(s-i), B_n is the inverse of A_n' modulo A_n, and deg D<n. The template degree is n+1; this is no bound on the degree after its source automorphism. For n>=2, acquiring unit Hessian anywhere in its full source orbit is exactly the row problem J((A_n,B_n*U)C)=1 with C in SL2(C[s,U]). Existence for some n>=2 is equivalent to failure of JC2/HC4. The native multi-point rows are non-Keller controls, and no independent Meng--Yang graph embedding or marked collision is inferred. Their global critical/Hodge data already pass the necessary tests, so those data alone do not settle the row problem.
- **An intrinsic plane closes three independent trace directions:** [HC4MYHF0/1](INTRINSIC_AFFINE_PLANE.md) applies to every four-variable polynomial, graph or mixed substitution. Two distinct constant jointly affine two-planes already force gradient invertibility by an elementary argument. Under hypothetical nonautomorphy the plane is unique, and the covectors admitting a pure trace E=t^2+kappa*t are exactly its two-dimensional annihilator. Three linearly independent pure directions therefore close every degree, even with zero or negative coefficients. Two directions determine the only possible plane: failure of any of its three directional Hessian identities also closes the candidate. Two pure directions alone retain an exact JC2 endpoint. [HC4MYHF2](INTRINSIC_AFFINE_PLANE.md#3-parallel-profiles-count-the-exact-residual-fibre-hc4myhf2) identifies complete parallel profiles with the reduced residual fibre, distinguishing empty and singleton ideals without proving properness.
- **Whole value fibres handle arbitrary source mixing:** [HC4MYHG0/1](GLOBAL_VALUE_ORBITS.md) uses the full critical scheme, critical Hessian values and entire three-dimensional value fibres under every actual four-variable source automorphism. A hypothetical nonautomorphic unit-Hessian polynomial has E(f=lambda)=t^3-n*t+n_lambda*t^2, where n counts its full critical points and n_lambda counts those at lambda. Linear-perturbation coefficient defects equal omitted plane-normalization lengths. Every zero-boundary graph y=R(x,p,q,r) or x=R(y,p,q,r), R(0,-)=0, has a degenerate full critical point and therefore no unit-Hessian member in its entire source orbit. The same is true of Q3(z)+b(w) for deg b>=3. No affine source slice need be preserved.
- **Eliminate x or y:** [HC4MYGB0/1/2](BASE_ELIMINATED_GRAPHS.md) closes both native base-eliminated charts by a complete boundary affine-direction classification and the two possible plane-projection ranks. The y chart's last polynomial equation has an uncancellable coefficient in every positive transverse degree. General base/dual mixing remains separate.
- **Every y=R(x,p,q) of positive (p,q) degree, independent of old r:** [HC4MYHA0/1/2/3](CYLINDRICAL_TRACE_CLOSURE.md) proves gradient invertibility under actual unit Hessian after every polynomial retained-three-dual automorphism over K[x]. The exact fixed-base square suspension has generic Euler at least two. A general cylindrical-trace theorem handles the whole exceptional trace, including nonlinear trace equivalence and vertical critical components. Arbitrary mixed p,q degrees and automorphism word lengths are covered. The entire family, including base-only R, has exactly the unit-ideal nonautomorphic endpoint below.
- **Nonlinear old-r-dependent y boundaries:** [HC4MYHC2](HODGE_TRACE_CLOSURE.md#4-nonlinear-old-r-dependent-y-boundaries-hc4myhc2) closes R(0,p,q,r)=P(a*r+b(p,q)), with a a nonzero field constant, and R0=P(C(p,r)) for any nonvertical plane coordinate C. P may be any nonconstant polynomial with arbitrary root multiplicities. All off-slice terms and retained-dual automorphisms are covered. Exact trace strata give a positive Hodge coefficient, including plane coordinates nonlinear in r.
- **Every one-affine-dual-direction y graph and every non-r opposite direction:** [HC4MYHD0/1/2](DIRECTIONAL_BASE_GRAPHS.md) permits w=a(s)*p+b(s)*q+c(s)*r with arbitrary polynomial coefficients and a nonzero row, followed by any polynomial graph R(s,w) of positive w degree and all retained-dual automorphisms. Every direction is closed in y=R(x,w); x=R(y,w) is closed when (a,b) is not identically (0,0). No unimodular row or coefficient/word degree bound is required. A fixed quartic/cubic rank calculation gives a positive whole-fibre Hodge coefficient in every admitted direction.
- **Pure-r opposite endpoint reduction:** [HC4MYHD3](DIRECTIONAL_BASE_GRAPHS.md#4-the-remaining-pure-r-opposite-route-reaches-the-known-endpoints-hc4myhd3) forbids unit Hessian for x=R(y,r) with zero or nonlinear R(0,r). A nonautomorphic member with constant or linear boundary would supply respectively the existing unit-ideal or singleton-zero plane Keller endpoint. This is a necessary reduction; a converse realization inside that graph subfamily and properness are not proved.
- **Every r graph after arbitrary automorphisms of (y,p,q) fixing x:** [HC4MYHB0/1/2](ISOLATED_TRACE_MIXING.md) permits old y to mix with p,q, in all degrees and word lengths. The x-zero trace is equivalent to a ternary quadric: its isolated Milnor number one and generic Euler two contradict the nonautomorphic lower bound three. This criterion needs no generic-base Euler information.
- **One Hodge trace controls arbitrary mixing:** [HC4MYHC0/1/3](HODGE_TRACE_CLOSURE.md) needs a single affine-hyperplane generic value fibre with E=t^2+kappa*t, t=u*v and kappa>0. Unit Hessian then forces gradient invertibility without isolated-point or generic-base data. This closes all H(U)+VW traces, general positive root-budget linear traces, and branched hyperplane trace triples with arbitrary off-slice corrections. The [earlier coordinate-triple theorem](ISOLATED_TRACE_MIXING.md#4-only-a-hyperplane-coordinate-triple-is-needed-hc4myhb3) remains a special case. No global substitution inverse or embedding is inferred.
- **Two traces admit arbitrary plane boundary curves:** [HC4MYHE0/1/2](TWO_TRACE_COMPATIBILITY.md) closes every y graph with R(0,p,q,r)=Z(p,r) of positive r degree and R(a,p,q,r)=b constant at one nonzero a. Boundary genera, singularities and multiplicities are unrestricted; all terms vanishing on both slices and all retained-dual automorphisms are allowed. One pure trace forces every parallel trace of a hypothetical nonautomorphic member to have E=t^2-m*t, m>=0; the boundary trace cannot have that signature. [HC4MYHE3](TWO_TRACE_COMPATIBILITY.md#4-one-place-unibranch-boundaries-need-only-one-trace-hc4myhe3) also closes disjoint unibranch zero curves with A1 normalizations from one trace alone, including every cusp power (r^m-p^n)^k with coprime m,n, without the second-slice condition.
- **Reparameterize a base-only graph x=h(y) or y=h(x):** [HC4MYHG2](GLOBAL_VALUE_ORBITS.md#3-the-two-known-endpoints-are-full-source-orbits-hc4myhg2) now permits every polynomial automorphism of all four source variables. If h(0)=0, unit Hessian is impossible. Otherwise the entire source orbit is the polynomial-coordinate orbit, and unit-Hessian existence is exactly existence of a unit-ideal plane Keller pair, equivalently failure of plane Keller surjectivity. One fixed graph x=1 already realizes that endpoint with the explicit relative inverses in [HC4MYGU1/2/3](RELATIVE_COORDINATE_ENDPOINT.md); arbitrary base/dual mixing creates no new endpoint.
- **Eliminate r or q:** [HC4MYGL1/2](CRITICAL_COORDINATE_JET.md) permits every polynomial base and relative dual automorphism. In fact only an invertible base first jet and an automorphic retained-dual specialization at one base point are required. The critical first jet gives a line on which the residual Keller pair is injective.
- **Eliminate p:** [HC4MYGC3](CURVE_TRACE.md) permits every polynomial base and relative affine dual change. [HC4MYGP1](P_GRAPH_SINGLETON_ENDPOINT.md) identifies the entire nonlinear-dual route exactly with the existing singleton-zero Keller problem. One fixed degree-six embedded graph, with an explicit polynomial parameter inverse, realizes its full oriented row orbit. [HC4MYHG2](GLOBAL_VALUE_ORBITS.md) also identifies the full four-variable source-automorphism orbit of that fixed graph with the nondegenerate quadratic orbit; arbitrary source mixing on it retains exactly the same singleton-zero endpoint. This last extension concerns the fixed graph, rather than every p graph under arbitrary source mixing.

The full q/r first-jet proof imports the reproduced Tang alternative and
[Gwozdziewicz's line-injectivity theorem](sources/critical-jet-import.json).
The narrower [square-axis q proof](Q_GRAPH_SQUARE_AXIS.md) needs only
Tang and an elementary degree lemma. The curve proofs also import the classical
Rentschler criterion from the [retained primary source](sources/automorphic-trace-import.json).
The p orbit identification imports [HC4MYOR1](../hc4-relative-row-orbit/PROOF.md).

```sh
timeout 60s .venv/bin/python research/hessian/hc4-triangular-graphs/verify_source_module_layers.py
```

[Latest replay receipt](receipt-source-module-layers.json) and [source context](sources/source-module-layers-context.json) bind finite-power character controls,boundary trace ranks,the layer budgets,and the explicit three proved imports. The category and topology arguments are written proofs.

The [plane trace-core replay](receipt-plane-trace-core.json) and [source context](sources/plane-trace-core-context.json) bind point-module and Cech controls,actual local Kummer trace factors,and the all-degree moved-sheet inequality. Universal statements are written proofs;no old checker is rerun.

The [canonical-energy replay](receipt-canonical-energy-traces.json) and [source context](sources/canonical-energy-trace-context.json) bind the cotangent energy identity,coefficientwise power traces,target-polynomial energy shifts and the explicit retained Tang dependency.

The [whole-source trace replay](receipt-whole-source-traces.json) and [source context](sources/whole-source-trace-context.json) bind the Newton and local sharpness controls,actual three-dimensional Keller trace poles,and differential trace-pairing identities. The all-degree proofs are written separately.

The [earlier canonical-boundary replay](receipt-canonical-boundary.json) and [source context](sources/canonical-boundary-context.json) bind valuation towers,tame inverse-different lattices,Koszul controls,the actual smooth quadric source chart and characteristic-p countercontrols. The global intermediate-field proof is written separately from these finite checks.

The [earlier coefficient-closure replay](receipt-coefficient-closure.json) and [source context](sources/coefficient-closure-context.json) bind integral inverse prefixes,Pluecker charts,all-prime formulas,explicit nonplane automorphisms and the specialization boundary. No generic elimination or N_D computation was run. Stacks01WC properness is the sole new proved import;Jelonek/Xavier closedness is credited.

The [earlier integral-descent replay](receipt-integral-cotangent-descent.json) and [source context](sources/integral-cotangent-descent-context.json) bind rational-plane/integral-frame controls,cotangent inverse and Hessian congruence identities,and the finite-precision countercontrol. Galois/lattice descent is a written proof using retained Tang and plane uniqueness.

The [earlier specialization replay](receipt-affine-plane-specialization.json) and [source context](sources/affine-plane-specialization-context.json) cover arbitrary-function shear identities,transported collisions,the explicit divided-equation escape,Hasse translation coefficients,even-matrix parity and a saturation rank control. Tang classification remains a retained written input,not a finite-checker conclusion.

The [earlier ramified-lifting replay](receipt-ramified-hc4-lifting.json) and [source context](sources/ramified-hc4-lifting-context.json) cover universal divided coefficients,the full cubic-quartic identity,sharp ramified thresholds and gradient collisions. No old arithmetic checker was rerun and no strong-Keller sufficiency is imported.

The [earlier critical-type replay](receipt-critical-gorenstein-fibres.json) and [source context](sources/critical-gorenstein-fibres-context.json) cover h-dependent critical-family equations,mixed Hessians,socle/generator counts,the genuine nonisolated A4 control and standard first-orbit-module presentation. No Rees classification or automatic critical Gorenstein theorem is imported.

The [earlier weighted-period replay](receipt-weighted-period-recovery.json) and [source context](sources/weighted-period-recovery-context.json) cover weighted residue duality,arbitrary-thickness local-time and analytic-volume identities,Hankel specialization defects and nonterminating pole iterates. No new external theorem or native saturation statement is imported.

The [earlier three-point module replay](receipt-three-point-translation-barrier.json) and [source context](sources/three-point-translation-barrier-context.json) cover the actual quadric kernel/plinth controls,whole three-point curve,puncture and collision identities,and pure alignment action. No new external proved theorem or all-action classification is imported.

The [earlier fixed-support replay](receipt-flat-plinth-families.json) and [source context](sources/flat-plinth-families-context.json) cover fixed-support module algebra,Nambu orientation,primitive cycle labels,specialization false positives and audited manuscript identities/weights. One proved PID module import is explicit;three manuscript sources are context/applicability audits only.

The [earlier finite-plinth replay](receipt-finite-plinth-descent.json) and [source context](sources/finite-plinth-descent-context.json) cover residue/Koszul algebra,the local critical-rank formula,primitive finite-cycle links and an actual nonreduced A4 countercontrol. Classical duality and principal-plinth literature is retained as context;no new external theorem is imported.

The [earlier radical-plinth replay](receipt-radical-plinth-rigidity.json) and [source context](sources/radical-plinth-rigidity-context.json) cover the universal presentation/normal-map algebra,primitive integral boundary blocks and the Euler-balanced negative control. One precisely checked naturality clause is imported from the retained full purity paper.

The [earlier flexibility replay](receipt-relative-flexibility.json) and [source context](sources/relative-flexibility-context.json) cover universal chart/flow identities,polynomial cancellation,derivation algebra,tangent spanning and the actual three-point family. The proof uses no external transitivity theorem.

The [earlier normal-form replay](receipt-global-normal-forms-seifert.json) and [source context](sources/global-normal-forms-seifert-context.json) cover weighted twist inverses,coefficient ideals,the explicit R_rho critical scheme and Seifert block consequences. Three precise proved imports come from the full Neumann--Norbury paper;no tameness is assumed.

The [earlier global-volume replay](receipt-global-vanishing-periods.json) and [source context](sources/global-vanishing-periods-context.json) cover exact residues,native-volume conversion,local node coordinates,the unimodular cycle matrix and collided-value controls. The integral tube construction and completeness argument are written proofs using retained suspension inputs.

The [earlier relative-motivic replay](receipt-motivic-family-barrier.json) and [source context](sources/motivic-family-barrier-context.json) check entire torsor/curve charts,sections and unit identities. Three explicit proved homotopy inputs are imported from the full versioned Asok--Dubouloz--Ostvaer paper;the symbolic checker does not implement homotopy theory.

The [earlier torus-visibility replay](receipt-torus-coordinate-visibility.json) and [source context](sources/torus-coordinate-visibility-context.json) cover the whole relative ideals,torus parameterizations,all-rho determinant class and critical-point controls. The sole new import is Chau's no-affine-line statement for the exceptional curve of a nonsingular polynomial plane map.

The [earlier global-sign-twist replay](receipt-global-sign-twists.json) and [source context](sources/global-sign-twists-context.json) cover the global sign twist and both inverse identities,critical schemes,coefficient-ideal generation,determinant repair and the symbolic cotangent Hessian formula. No new external theorem is imported.

The [earlier coordinate-specialization replay](receipt-coordinate-specialization.json) and [source context](sources/coordinate-specialization-context.json) cover Hamiltonian specialization,actual generic-coordinate inverses,exceptional pole factors,the complete plane pencil,relative critical ideals and affine trace alignment. The sole new classification input is the two-variable Rentschler theorem;the all-degree proof and full-source boundary are in JW0/1/2.

The [earlier two-slice replay](receipt-two-slice-recovery.json) and [source context](sources/two-slice-recovery-context.json) cover the global slice corrections,exact inverse denominator cancellation,normal hypersurface degeneracies,unimodular row completion,unsquared determinant repair and canonical model identities. Five proved algebraic imports are explicit;the all-degree proofs and unresolved action/Keller-generator boundaries are in JV0/1/2.

The [earlier quotient-hole replay](receipt-quotient-holes.json) and [source context](sources/quotient-holes-context.json) cover the raising/lowering obstruction,exact quotient holes,Cech transition and annihilators,common-orientation residues,and the exact all-order period identities. The full proofs and open boundaries are in JU0/1/2.

The [earlier primitive-quotient replay](receipt-primitive-quotients.json) and [source context](sources/primitive-quotient-context.json) cover the local-slice identities,whole-source primitivity including critical fibres,invariant replicas,proper plinth tests and the actual coordinate-action inverse certificate. The all-degree proof and precise open boundaries are in JT0/1/2.

The [earlier global-action replay](receipt-global-additive-barrier.json) and [source context](sources/global-additive-barrier-context.json) cover both global value-presentation composites,four actual polynomial actions,the kernel-chart and fixed-point identities,every module coefficient case,and the boundary between a necessary action and actual unit Hessian. The all-degree proof is in JS0/1.

The [earlier one-coordinate replay](receipt-one-hyperbolic-coordinate.json) and [source context](sources/one-hyperbolic-coordinate-context.json) cover the generic curve units,global cylinder scale,whole primitive-row completion,first-normal slope and all-higher-jet identities,including the nonlinear-coordinate and K2 unequal-slope boundaries. The all-degree proof is in JR0/1.

The [earlier normal-model replay](receipt-two-normal-models.json) and [source context](sources/two-normal-model-context.json) cover all global normalization maps and both inverses,the actual suspension and generic plane fibre,all critical-value collisions,and the entire relative rank-one surface and global section. The exact reduced frontier and natural-base exclusions are in JQ0/1.

The [earlier all-affine replay](receipt-resonant-cubics.json) and [source context](sources/resonant-cubic-context.json) cover the actual mixed-plane quartic branches, complete two-plane divisor strata and singular groups, all resonance source inverses, and every passing critical determinant and whole-value jump of the remaining cubic models. The all-degree reduction and its precise open limits are in JP0/1/2.

The [earlier isotropic-kernel replay](receipt-isotropic-kernels.json) and [source context](sources/isotropic-kernel-context.json) cover the exact trigonal polar and infinity branches, binary-cone profiles, zero-chain moments, every exceptional quadratic critical group and global source maps. The all-degree reduction and the sole remaining quartic frontier are in JO0/1.

The [earlier compatible-center replay](receipt-compatible-inner-centers.json) and [source context](sources/compatible-inner-center-context.json) cover stabilized plane/ternary profile controls, every exceptional plane critical point, actual rank/force fibre maps, compatible-center realizations of all four tests and the global coordinate-inner inverses. The all-degree proof is in JN0/1; the remaining cubic/quartic coordinate-inner forms and actual endpoints stay open.

The [earlier rank-one replay](receipt-rank-one-pencils.json) and [source context](sources/rank-one-pencil-context.json) cover global null-shear inverses, the complete stabilized native critical mass, invariant zero-primary norms, exact finite/infinity data, every resonant multiplier inverse and all four actual graph tests. The all-degree proof is in JM0/1; no bare source-orbit conclusion is assumed to stabilize.

The [earlier Jordan-polar replay](receipt-jordan-polar-curves.json) and [source context](sources/jordan-polar-invertible-context.json) cover Jordan polar and infinity controls, whole open/closed fibre inverses, constant rank profiles, the marked hyperbolic fibre model and every low-genus singular stationary group. The all-degree proof is in JL0/1; the singular spectral boundary and actual endpoints remain open.

The [earlier suspended-multiplier replay](receipt-suspended-multipliers.json) and [source context](sources/suspended-multiplier-context.json) cover derivative rigidity, repeated-root critical schemes, every special-value profile, polynomial normal-form inverses, compensated graphs and nonlinear base/dual orientation transport. The all-degree proof is in JK0/1; the two full-source quadratic models remain open.

The [earlier affine-polar replay](receipt-affine-polar-curves.json) and [source context](sources/affine-polar-semisimple-context.json) cover the polar degree and branch controls, whole Witt strata and inverses, high-primary fibre maps, scalar low-genus exceptions and a genuine singular rank-one two-critical-point normal form. The all-degree proof is in JI0/1; no connectedness or suspension-classification theorem is imported.

The [earlier stationary-pole replay](receipt-stationary-pole-masses.json) and [source context](sources/stationary-pole-mass-context.json) cover the complete nine-critical-point example, 2m-point stationary groups, the orthogonal proper-Krylov split, isotropic-vector distinctions and a genuine multiple-pole boundary residue without finite critical points. The all-degree proof is in JH0/1; no unproved multivariate residue theorem is imported.

The [earlier elliptic-cone replay](receipt-elliptic-cone.json) and [source context](sources/elliptic-cone-context.json) cover exact conic/discriminant identities, generic branch simplicity including chart infinity, distinct puncture counts, full resonant source-map inverses, the whole-fibre Hodge term and negative geometric controls. The all-degree geometric proof is in JG0/1; no new unproved genus theorem is imported.

The [earlier affine-pencil replay](receipt-affine-pencil-residues.json) and [source context](sources/affine-pencil-residue-context.json) cover exact stationarity and critical determinant identities, the full reciprocal-Hessian sum, cyclic Jordan boundaries, the unresolved infinity-resonance control and genuine fixed graph endpoint converses. The all-degree proof is in JF0/1/2; no toric or multivariate residue theorem is imported.

The [earlier centered-Jordan replay](receipt-centered-jordan.json) and [source context](sources/centered-jordan-context.json) cover the exact open-stratum inverses, every zero slice, generic branch divisors, ramification/genus arithmetic, whole-value Hodge profiles and critical-count exceptions. The all-degree and complete-block proof is in JE0/1; the foundations are retained explicit inputs.

The [earlier nilpotent-pencil replay](receipt-nilpotent-pencils.json) and [source context](sources/nilpotent-pencil-context.json) cover every nilpotent block type, divided-difference polynomiality, full triangular inverse identities, the native-Hessian warning and the actual graph potential. The all-outer-degree and exact endpoint-orbit proof is in JD0.

The [earlier spectral-pencil replay](receipt-spectral-pencils.json) and [source context](sources/spectral-pencil-context.json) cover exact critical Hessian factors, scalar coefficient/spectral identities, nonlinear source and ambient inverse certificates, and the genuine nonsemisimple positive boundary. The all-degree proof is in JC0/1; no quartic or matrix-periodicity theorem is imported.

The [earlier formal-divisor replay](receipt-formal-divisor.json) and [source context](sources/formal-divisor-context.json) cover exact inverse-jet compositions, fat-divisor orders, the rational degree threshold and the pole/zero boundary distinction. The formal thickening theorem is attributed and reproduced; the all-order and degree arguments are in JB0.

The [earlier extra-prime replay](receipt-extra-prime-excision.json) and [source context](sources/extra-prime-excision-context.json) cover divisor unit-ideal identities, nonreduced and etale-open controls, cusp/tacnode/nodal budgets and a genuine nonlinear embedding. The section theorem and one-place Milnor identities are explicit inputs; the all-degree arguments are in JA0/1/2.

The [earlier one-place trace replay](receipt-one-place-trace.json) and [source context](sources/one-place-trace-context.json) cover the intrinsic genus identities, normalization-ring reconstruction, a genuine nonlinear nodal-trace embedding and full-preimage controls. Vistoli and both foundational Stacks inputs are explicitly imported; the all-degree argument is in IZ0/1/2.

The [earlier contractible-pencil replay](receipt-contractible-pencil.json) and [source context](sources/contractible-pencil-context.json) cover exact boundary and genus identities, etale nonfinite controls, a genuine nonlinear cusp-trace embedding and the nodal remaining-boundary control. Lin-Zaidenberg classification is an explicit new external input, not a finite-checker conclusion.

The [earlier radical-trace replay](receipt-radical-trace.json) and [source context](sources/radical-trace-context.json) cover the all-degree plane-family inverse identities, exact divisor factors, radical first jets and a genuine nonlinear hyperbolic-trace embedding. The unbounded and geometric conclusions are proved in IX0/1/2 with explicit retained external inputs.

The [earlier middle-cycle replay](receipt-quadric-middle-cycle.json) and [source context](sources/quadric-middle-cycle-context.json) cover the exact section and rotation identities, local divisor charts, transversality and genuine nonlinear embedding controls. Integral cycle-map and Poincare-duality foundations are explicit external inputs; the all-degree proof is in IW0/1/2.

The [earlier primitive-helper replay](receipt-primitive-helpers.json) and [source context](sources/primitive-helper-context.json) cover full compensation and gauge freedom, primitive-helper Jacobian ideals, normality and exact Legendre inverse identities, including non-Keller and nodal negative controls. The characteristic-zero weak-type argument is reproduced with attribution; retained HW0 and normality inputs remain explicit.

The [earlier free-projection replay](receipt-free-projection-trace.json) and [source context](sources/free-projection-trace-context.json) cover helper ranks, actual closed parameter recovery, unrestricted free-variable degree, the exact first jets and a unit-Hessian substitution with a killed source tangent. The one-line theorem and polynomial generative input are retained explicit imports.

The [earlier intrinsic Nambu replay](receipt-intrinsic-nambu.json) and [source context](sources/intrinsic-nambu-context.json) cover the actual Legendre-energy identities, finite orbit certificates, inverse-Hessian cone and nonlinear quasi-translations, including the automorphic plane-shaped-cone control. Nambu-energy identities have recent primary attribution and are reproved; Tang and regular-parameter inputs remain explicit.

The [earlier universal-modification replay](receipt-universal-modification.json) and [source context](sources/universal-modification-context.json) cover the universal graph, exact fibre derivative identity, source-correlated class lattice, helper nonflatness and actual unit-Hessian gradient inverse. Affine-modification terminology has primary context; the needed calculation is reproved directly.

The [earlier two-helper replay](receipt-two-helper-curves.json) and [source context](sources/two-helper-curve-context.json) cover finite-line curve surjectivity, multiplier-zero row descent, generic helper flatness, exact source/helper component correspondence and the complete admitted HC4 overlap closure. The curve-image assertion is proved directly; the helper flatness uses retained Stacks 00R4.

The [earlier one-helper replay](receipt-one-helper-overlap.json) and [source context](sources/one-helper-overlap-context.json) cover generic helper descent, both multiplier cases, finite self tensors, admitted nonprimitive content and actual two-root-variable HC4 compatibility. The self-tensor argument is reproved directly; no new unproved classification is imported.

The [earlier additive-action replay](receipt-additive-actions.json) and [source context](sources/additive-action-context.json) cover derivation degrees with large kernels, the two-sided component argument, exact natural kernels, primitive global descent and actual induced source actions. The degree lemma is reproved directly; no new unproved classification is imported.

The [earlier overlapping-output replay](receipt-overlapping-outputs.json) and [source audit](sources/overlapping-output-import.json) cover algebraic coefficient models, guarded coefficient specialization, independent fields, overlapping row lattices, generic Bezout presentations and actual HC4 controls. The finite-type Nullstellensatz is explicitly imported; all earlier proof and replay bytes are preserved.

The [earlier reflexive-pencil replay](receipt-reflexive-pencil.json) and [source audit](sources/reflexive-pencil-import.json) cover rank-two denominator bounds, the faithful-flatness boundary, exact content fibres and actual HC4 compatibility. The reflexive-surface and miracle-flatness lemmas are explicitly imported.

The [earlier parameter-stability replay](receipt-parameter-stability.json) and [source audit](sources/parameter-stability-import.json) cover parameter-stable root transitivity, moving class rows, exact content fibres and adapted global frames. The Galois-tower inertia-surjectivity lemma is explicitly imported, using unchanged bytes from an existing repository source cache.

The [earlier multivariable replay](receipt-multivariable-suspensions.json) and [source audit](sources/multivariable-suspension-import.json) cover generative factor transfer, exact fixed-multiplier class lattices and the global multivariable suspension proof. The same closed-polynomial implication is explicitly audited over C(T); generic geometric integrality is proved by a direct field argument. No special-fibre classification or coordinate theorem is imported.

The [earlier parameter-norm replay](receipt-parameter-norms.json) and [source audit](sources/parameter-norm-import.json) cover the whole plane polynomial-output proof, the all-dimensional parameter-norm theorem and the direct critical-value scaling lemma. Ordinary one-variable Luroth and structural smoothness inputs are explicit; no Ritt or dynamical good-reduction classification is imported.

The [earlier central-quotient replay](receipt-central-quotients.json) and [source audit](sources/central-quotient-import.json) cover relative central-series descent, algebraic product-cover quotients and finite-pole binomials. The comparison theorem is used through full faithfulness for already algebraic covers; essential surjectivity is not used.

The [earlier branch replay](receipt-branch-parameters.json) and [source audit](sources/branch-parameter-import.json) cover the plane branch, rational-parameter, inertia and complete cubic proofs. They reuse seven explicit retained primary inputs; no new external theorem is imported. The checks are exact controls for the written proofs, not an exhaustive degree search.

The [earlier multiplicity replay](receipt-pencil-multiplicities.json) and its multiplicity
bound reuses the explicit closed-rational equivalence and projective
Riemann--Hurwitz inputs. The inequality, equality counterexample,
coefficient certificate and critical-point expansion are proved directly.
The [context audit](sources/pencil-multiplicity-context.json) keeps inspected
orbifold literature separate: maxima are used here, and no orbifold-cover
classification is imported.
The [reduced-pencil replay](receipt-reduced-pencils.json) retains its exact
earlier controls. The rational-pencil
proof explicitly imports Petravchuk--Iena's closed-rational-function
equivalence and reuses projective Riemann--Hurwitz. Its
[source manifest](sources/reduced-pencil-import.json) records the exact
statements; the point-form, ramification and scaling arguments are written
in the note. Critical points are allowed, but constant finite pencil
members and uncertified exceptional fibres are excluded. The one-exception
rational test additionally requires numerator-denominator independence;
the dependent cubic counterexample is retained explicitly.
The [all-factor replay](receipt-factor-resultants.json) retains its exact
earlier controls. The all-factor proof
uses smoothness, primitive factors, resultant identities and the infinity
valuation directly, with no new external theorem. The
[source audit](sources/factor-resultant-context.json) records the inspected
ruled-surface literature as context and identifies its missing Q-factoriality
and almost-surjectivity hypotheses. Its nonnormal-order conclusion reuses
the retained normal-monogenic theorem.
The [derivative replay](receipt-derivative-modification.json) retains its exact
earlier controls. The derivative
surface proof retains the original A2 source and reuses the explicit qualitative
Jelonek--Lason theorem. Normal monic rescaling additionally imports Stacks
00SB, 031S and 07PF, recorded in its
[source manifest](sources/derivative-modification-import.json).
The [critical-factor replay](receipt-critical-support.json) retains its exact
earlier controls. The critical-factor
support proof uses specialization, UFD factorization and the retained
unit lattice directly. The fixed-model closure reuses the explicit
Jelonek--Lason nonproperness theorem; its
[source context](sources/critical-support-context.json) adds no new
external theorem or S-unit height bound. The
[valuation/critical-curve replay](receipt-valuation-critical.json) accompanies a
proof that explicitly imports the qualitative Jelonek--Lason
nonproperness theorem in its [source record](sources/valuation-critical-import.json).
The root-valuation, normalized-label and cubic critical-curve arguments
are written directly, using the retained Hurwitz input. The
[normal-boundary replay](receipt-normal-boundary.json) accompanies a proof that
explicitly imports DVR multiplicativity and retains the earlier purity,
finite-cover comparison and normalization inputs in its
[source record](sources/normal-boundary-import.json). The boundary-basis
and class-stabilizer arguments are written directly. The
[nilpotent-cover replay](receipt-nilpotent-covers.json) accompanies a proof that
derives the group lemma and common-prime radical-tree bound directly;
general wreath embeddings are retained as [primary context](sources/nilpotent-cover-context.json).
Its full-group conclusion applies specifically to nilpotent covers, while
the general monogenic result below gives root transitivity. The
[leading-fibre replay](receipt-leading-fibre.json) accompanies a proof that
explicitly imports Riemann--Hurwitz and proves its finite-affine-line
constant-field consequence. Its [source audit](sources/leading-fibre-import.json)
excludes the inspected Tian--Sun article's Keller conclusions: affine etaleness
alone does not exclude finite boundary punctures. The
[monogenic-cover replay](receipt-monic-covers.json) accompanies a proof that
records its normality, ramification, purity and finite-cover comparison
[primary inputs](sources/monic-cover-import.json) explicitly; the resultant,
integral-root collision and invariant root-lattice arguments are written in
the note. The [regular-parameter replay](receipt-regular-parameter.json)
uses a proof that explicitly imports the classical noncomposite-to-closed polynomial
[theorem](sources/regular-parameter-import.json); its field denominator step,
finite-cover disjointness and monodromy application are written in the note.
The [abelian class-group replay](receipt-abelian-classgroup.json)
reproduces square-free pullbacks, the cyclic valuation obstruction, class
stabilizer descent and the quadric/Danielewski class groups directly, with
[primary square-free context](sources/abelian-classgroup-context.json)
separate from proof imports. The
[smooth-infinity replay](receipt-smooth-infinity.json) accompanies the
arbitrary-degree output-base proof, which explicitly imports the classical smooth complete-intersection
[Picard theorem](sources/smooth-infinity-import.json); its filtration, global
denominator and disjoint-component arguments are written in the proof. The
[positive-weight image replay](receipt-positive-weight-images.json)
uses direct isometry identities, weighted Euler and a finite-weight formal
inverse argument, with [primary context](sources/positive-weight-images-context.json)
separate from any group-generation or embedding classification. The
[quadratic-plane replay](receipt-quadratic-plane.json)
reproduces both elementary plane Keller inverses and keeps
[classical context](sources/quadratic-plane-context.json) separate from proof imports.
It retains Tang, plane uniqueness and the existing endpoint classification. The
[proper-Morse replay](receipt-proper-morse-boundary.json)
retains only the forward finite-cover comparison and affine proper/finite facts
in its [primary import record](sources/proper-morse-boundary-import.json), beyond
the existing Tang, component and section inputs. The finite-free matrix lemma
is proved by monic division. The [relative critical-graph replay](receipt-relative-critical-graph.json)
proof retains [HC2 explicitly](sources/relative-critical-graph-import.json), together
with the existing Tang, curve/Hodge and endpoint inputs. The
[quadric output-base replay](receipt-quadric-output-base.json)
retains the basic [Nagata criterion](sources/quadric-output-base-import.json)
beyond the retained stable-hypersurface algebra and birational Keller endpoint;
the global filtered primitive lemma is proved in the note. The
[quadratic symmetry replay](receipt-quadratic-symmetry.json) uses
literal ambient symmetry identities and the retained three-affine and hyperplane
results, with [primary literature context](quadratic-symmetry-literature.json)
kept separate from any imported classification. The
[endpoint geometry replay](receipt-two-affine-geometry.json) uses the
retained Tang alternative, source-orbit equivalences and basic Hodge invariant.
The three-projection proof retains the
fully reproduced Formanek criterion, the retained affine-line endpoint and Tang;
it assumes no relative line-embedding rectification theorem.
The ambient quadratic proof reuses the
retained HC5T1 frames and Tang alternative, with no new classification import.
The fixed-potential boundary proof
uses elementary algebra and the retained Tang/stable-hypersurface inputs.
The projection profile retains
[classical plane stable-coordinate context](sources/projection-profile-context.json).
The stable obstruction retains its
[basic field-theory sources](sources/stable-hypersurface-import.json); the classical
reconstruction criterion retains its [attribution record](sources/reconstruction-import.json).
The finite-ideal row theorem has related
[primary literature context](sources/finite-fibre-literature-context.json), with its needed
specialization proved constructively. The basic Hodge--Deligne ring invariant and curve formula have an explicit
[primary-source import](sources/hodge-trace-import.json); the retained Tang and affine-component inputs remain required.
The additional degree-embedding input in the cylindrical-trace theorem is
[Garcia Barroso--Ploski Theorem6.9](sources/degree-embedding-import.json), reusing the retained primary manuscript.
The normalization proof has explicit
[normalization/conductor imports](sources/critical-normalization-import.json).
The opposite-family proof has its explicit
[Chau/Orevkov import manifest](sources/opposite-affine-import.json). The finite
identities and controls do not replace the written all-degree arguments or
independently certify the imported theorems.

Retained results have additional or different strength:

- Intermediate base-chart proofs retain their exact receipts and narrower strengths: [affine-p](RELATIVE_CRITICAL_LOCUS.md), [Newton/separated](NEWTON_CRITICAL_GATE.md), [boundary components](CRITICAL_COMPONENTS.md), [normalization](CRITICAL_NORMALIZATION.md), [global Euler budget](GLOBAL_CRITICAL_BUDGET.md), [first full y closure](Y_GRAPH_CLOSURE.md), and [opposite no-escape/threshold theorem](OPPOSITE_NONLINEAR_FRONTIER.md). The [exceptional-fibre theorem](EXCEPTIONAL_FIBRE_CLOSURE.md) closes all their former positive-p boundaries in both orientations. HC4MYHA2 now extends the y closure to arbitrary mixed p,q dependence. The [opposite affine-p theorem](OPPOSITE_AFFINE_ENDPOINT.md) retains its stronger unit-Hessian nonexistence statement.

- [Automorphic trace](AUTOMORPHIC_TRACE.md), `HC4MYGV0/1/2`: the hyperplane version, including its Hamiltonian endpoint lemma and complete tangent-plane argument.
- [First slice](FIRST_SLICE_INVERTIBILITY.md), `HC4MYGF1/2`: arbitrary x-divisible corrections to a base-independent **Keller** trace pair, without assuming a trace inverse. It also supplies an explicit unit-Hessian embedded control with an affine plane.
- [Plane-Keller encoding](FIRST_SLICE_INVERTIBILITY.md#6-sharp-trace-boundary-the-larger-substitution-class-contains-jc2), `HC4MYGF3`: the larger arbitrary-substitution class with merely Keller y-dependent traces contains every plane Keller cotangent lift. [HC4MYHI2](RECONSTRUCTION_EMBEDDING_OBSTRUCTION.md) now proves that no nonautomorphic pair admits a closed embedding preserving these exact A/B data, irrespective of the remaining dual function. The arbitrary-substitution equivalence remains valid.
- [Matching two divisors](TWO_DIVISOR_EXTENSION.md), `HC4MYGD1/2`: a stronger no-plane conclusion on its stated subclass.
- [Earlier triangular and base-independent theorem](PROOF.md), `HC4MYGT1/2`: retained no-plane arguments and exact marked-covector formulas.

The old-r-independent y chart, both polynomial-coordinate r boundary
classes, every one-affine-direction y graph and every non-r opposite
direction are closed. The remaining pure-r opposite route reduces to
the two existing plane endpoints. The r chart allows full (y,p,q) mixing
fixing x. Compare whole Hodge traces in independent affine directions
before further classification: three pure directions close the candidate;
two determine its only possible plane and three polynomial identities.
The [isolated-trace bound](ISOLATED_TRACE_MIXING.md) and cylindrical-trace
criterion remain available for other signatures. The two-trace closure also
admits arbitrary q-independent r boundaries with one constant nonzero-base
slice, while one-place unibranch boundaries need no second slice. Other
old-r-dependent y boundaries and more general opposite dual dependence, mixing outside
these trace criteria and uncovered embeddings remain outside the closures.
Zero and negative Hodge coefficients are not a new properness theorem.
Base-only graph reparameterizations have the exact unit-ideal endpoint
above; no such Keller pair or plane surjectivity theorem is supplied.
Base and relative dual automorphisms
preserving the fibration no longer reopen either the q or r chart.
The nonlinear p route is the same unresolved singleton-zero properness
problem as the relative Schur row orbit, rather than a new graph-degree
search. Its coefficient ideal is the origin maximal ideal, whereas the
new omitted-point endpoint has unit ideal. No singleton-fibre or
ideal-equality argument proves properness.
No preservation of a Hessian determinant under nonlinear
coordinate change is assumed. No general HC4 or JC2 solution, formal
verification, independent external review or novelty priority is claimed.
Mathematical status remains authoritative in [MATH_STATUS.json](../../MATH_STATUS.json).

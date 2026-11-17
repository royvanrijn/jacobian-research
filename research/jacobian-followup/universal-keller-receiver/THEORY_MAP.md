# Universal Keller theory: constructions, rules and remaining gaps

Updated15 November2026. This is a synthesis and navigation map. Mathematical status is recorded only in [MATH_STATUS.json](../../MATH_STATUS.json), with the canonical proofs linked below. It covers the existing quadratic-gauge compiler; statements about every possible Keller map are explicitly distinguished.

## Current canonical existence: the trace-filter receiver

The [canonical filter proof](../universal-keller-canonical-filter/PROOF.md) now constructs the non-atomic canonical polynomial-bundle receiver in EVERY rank N>=2. For a finite etale algebra A, let V_A=A^3 plus three trivial coordinates. Evaluate the retained cubic in A, take C_A=Tr((3X^2-X)/2), and use the fixed interpolation polynomial lambda_N selecting the integer1 among0,...,2N. One extra cubic block removes every unwanted configuration through an actual empty fibre. The complete marked fibre is canonically Spec(A).

| Requirement | Current result |
|---|---|
| Primitive generator, good line or input-wise basis | None chosen; the formula uses unit, multiplication and trace. |
| Whole source and source/target identification | The same specified vector bundle V_A, of rank3N+3; all geometric fibres are whole affine spaces. |
| Effective gluing and natural S_N symmetry | Constructed over every Z[1/(2N)!]-scheme, including nonreduced bases, and compatible with all algebra isomorphisms and base changes. |
| All subgroup-fixed restrictions | Actual affine-space restrictions of degree3^(s+1) for s label orbits, with exactly the fixed labels in the marked fibre. |
| Degree and canonical volume | Generic degree3^(N+1), coordinate degree4N; volume lift on V_A plus det(V_A)^* has rank3N+4. |
| Atomicity | FALSE for this construction: two canonical equivariant factors have degrees3 and3^N. |
| A natural trivialization over every base, or one fixed all-scheme universal map | Not supplied; the fixed-scheme obstruction remains. |

UAFILTER1, UACANFIB1, UACANFIX1 and UACANVOL1 are the current canonical sources. Earlier references below to an unconstructed whole-affine canonical boundary concern the older degree-N compiler, atomic or stronger fixed-space variants. Non-atomic polynomial-bundle existence and simultaneous fixed-locus gluing are no longer open. The next main task is an atomic coupling outside the output-controlled and rank-one cubic-return classes below, preserving this exact canonical fibre, or an obstruction to that stronger requirement. No independent proof or publication-priority review is claimed.

## Coupling boundaries and atomicity over the ground field

The [feedback-rigidity proof](../universal-keller-feedback-rigidity/PROOF.md) closes two characteristic-zero coupling classes in every polynomial degree.

- UAPARALLEL1 allows every block coefficient to depend on all retained first outputs. Every subset of blocks still gives a polynomial Keller factorization, even with cyclic dependence. The intermediate fields form a Boolean lattice; maximal factor length is the number of blocks.
- UACUBICKER1 identifies the auxiliary cubic translation derivation's rational and polynomial kernels as k(R1,R2) and k[R1,R2], where R1=K2+3K1 and R2=K3-2K1. UASCALARFEED1 then classifies every polynomial scalar control in the trace-filter architecture: the map is Keller exactly when the control lies in k[R1,R2]. Every allowed control is recoverable from the auxiliary outputs, giving a reverse-order polynomial factorization. Arbitrarily high scalar degree cannot remove this factor.
- UAWEILATOM1 evaluates a primitive degree-d seed in a finite etale algebra A. Its maximal stable factor length over k is the number of field components of A, and after K/k it is the number of components of A tensor K. It is k-stably self-atomic exactly when A is a field. For rank greater than one it decomposes over splitting constants, so it is not absolutely stably atomic. The natural rank-N cover is only a selected subconfiguration before filtering, not its complete marked fibre.

These are actual polynomial factors and written general arguments, with selected exact controls. The next atomic construction must leave both closed coupling classes. The joint-feedback theorem below now covers every polynomial joint selector with the same one-direction cubic return. A formal matrix satisfying the feedback determinant identity is not yet a whole-source map with the required canonical fibre. The general canonical absolute-atomic problem remains open.

## All joint selectors and the auxiliary constant ring

The [joint-feedback proof](../universal-keller-joint-feedback/PROOF.md) now closes the remaining selector gate in the retained polynomial-profile, one-constant-return architecture. Let Delta=(DG)^(-1)*mu and let C be its FULL polynomial first-integral ring. Process active compiler blocks (those with Delta g nonzero), leaving inactive blocks untouched; call those intermediate coordinates Z. All original signals remain available.

| Claim | Exact result | Boundary |
|---|---|---|
| UAJOINTKER1 | The lifted derivation D=Delta-W has rational kernel Frac(C)(Z) and polynomial kernel C[Z]. | C is not assumed finitely generated, polynomial or output-recoverable. The precise generic-field Bertini lemma handles each active inverse pencil. |
| UAJOINTSEL1 | A polynomial joint selector beta(P,U) has constant nonzero Jacobian c exactly when beta=(1-c)*t+h(Z,C), with Delta t=1. | One constant return direction and polynomial profiles in U and retained first outputs are required. A polynomial slice does not imply local nilpotence. |
| UAJOINTFACT1 | If C is polynomial, expose its coordinates, apply the time shear, process inactive blocks, then apply the transverse auxiliary map. All are actual whole-affine self-factors. | Stable polynomiality after identities or constants gives the corresponding stable obstruction; its failure is not an atomicity proof. |

This includes every joint selector for the cubic in every constant direction. The displayed maximal atomic cubic chain has length N+1 and degrees {3,d_1,...,d_N}; all maximal chains are not classified. Products of two cubics demonstrate the stronger polynomial-kernel criterion even when those invariants are not recoverable from the original auxiliary outputs. Locally nilpotent directions with a slice also give actual source-automorphism factors.

For N>=2 the compiler stages alone prevent atomicity when the auxiliary is invertible. The complementary endpoint N=1 with an auxiliary automorphism is polynomial left-right equivalent to one known absolutely stable atomic compiler profile; it is not a new stable atomic class.

The [source-only selector theorem](../universal-keller-selector-rigidity/PROOF.md) retains its exact intersection with k[P]: selectors can use active blocks only through their signals. The [earlier rank-one proof](../universal-keller-rank-one-feedback/PROOF.md) retains its output-recovery and geometric generic-line criteria. Their former open joint-selector gate is now classified by UAJOINTSEL1.

For natural S_N symmetry FIXING the auxiliary source, all-active feedback cannot retain the complete transitive N-label fibre for N>1: repeating one block preimage gives a diagonal fixed point. This argument needs no polynomiality of C. Within this one-direction architecture, a remaining auxiliary-only canonical candidate must be all-inactive and have a non-locally-nilpotent direction whose full constant ring lies outside the stable-polynomial criterion. The higher-return extension below now distinguishes factor-preserving descended selectors from actual non-descended escapes. Other compilers and new whole-source profile geometry remain separate departures. Exact marked fibre, compatible descent and absolute atomicity remain required.

## Higher return rank: the exact boundary and marked deformations

The [higher-return proof](../universal-keller-higher-return/PROOF.md) extends the common-kernel theorem while separating it from higher-rank selector classification.

| Claim | Result | Boundary |
|---|---|---|
| UAMULTIKER1 | In every constant return rank, the lifted commuting frame has common rational kernel Frac(C)(Z) and polynomial kernel C[Z]. For the cubic, full rank and quotient rank one, the indicated auxiliary constants are polynomial. | C is the FULL auxiliary common kernel; general finite generation or polynomiality is not assumed. |
| UABASICTIME1 | Selectors PROVED to belong to C[Z,t] give ordinary time Keller-map families and retain actual compiler/auxiliary factors when C is polynomial or stably polynomial. | Higher-rank determinant one does not force descent. Time-map invertibility would contain JC(r) and is unnecessary for the factorization. |
| UAMARKEDDEF1 | Profiles I(sigma)^ell*h_j retain the entire restricted morphism over V(I(a)^ell), universally. For I(T)=T(3T-1)(4T+1), atomic auxiliary maps of every generic degree d>=3 agree with the cubic to any prescribed finite order along the three marked planes. | Generic degree and stable class are not determined by this finite-order data. Canonical trace receiver deformations retain their complete Spec(A) fibre and dimensions, but still factor. |
| UARETURNESC1 | A genuine non-descended rank-two return changes auxiliary degree3 to4 with the three filter fibres unchanged. A full-rank return cancels the auxiliary cubic and recovers a known stable atom. | Neither supplies a new atomic canonical receiver or a new stable class. The fixed-signal rank-two source is a hyperbola complement with nonconstant units, not affine plane. |

The recent Wang-Yan nonlinear nilpotent family has an explicit inverse when used as an ordinary time block; this makes it another descended factor-preserving example. Its use and the distinction between total coordinate degree and generic covering degree are recorded in the [literature audit](../universal-keller-higher-return/LITERATURE_AND_GAPS.md).

Relative to the fixed cubic, the native higher-profile deformation combined with the nonzero trace shift spans all three constant return directions. This is an obstruction only for that specified architecture, not for all rank-two feedback. The next atomic candidate needs a non-descended higher-rank interaction preserving the complete canonical marked scheme and proving full-map atomicity, or another stated departure in compiler or auxiliary geometry. Raising the deformation degree alone leaves explicit factors.

## Source-root corrections: acyclicity replaces a degree search

The [root-feedback proof](../universal-keller-root-feedback/PROOF.md) closes an additional characteristic-zero class beyond descended time maps. Retain every block first output sigma_i and add arbitrary polynomial corrections in all sigma and source roots X_i=x_i*q_i to the other outputs. These roots need not recover from old outputs.

| Claim | Result | Boundary |
|---|---|---|
| UAROOTDAG1 | Keller is equivalent to the exact self equation sigma_i^2 H_i,Xi=2X_i^2 F_i,Xi and an acyclic effective cross-root graph. Independent r_i diagonal variables isolate every principal minor. | All first outputs are retained and corrections lie in k[u,sigma,X]. Ordinary nilpotence is weaker. Positive characteristic and nonreduced-base cycles are explicitly outside. |
| UAROOTFACT1 | Topological processing gives actual polynomial whole-affine stages, each an absolute stable atom of degree at least3. Dependency order ideals give polynomial intermediate models. | A genuinely non-descended rank-two degree4-times4 example still factors. This is a displayed chain, not a classification of all intermediate fields or all maximal chains. |
| UAROOTSYM1 | An invariant DAG gives equivariant orbit factors. Transitivity on ALL blocks makes it edgeless and forces complete fixed-target fibres T0^N, excluding the natural N-label cover for N>1. | With extra auxiliary orbits only the orbit factorization is asserted. The trace receiver shifts its auxiliary first output and is outside this class. |

The matrix lemma is classical and credited to Berson's Lemma5.7/Remark5.8 and its van den Essen reference; the compiler determinant formula, full correction normal form and root-dependency application are the written deductions here. See the [literature boundary](../universal-keller-root-feedback/LITERATURE_AND_GAPS.md). Increasing polynomial root degree or finite-order marked vanishing cannot remove factors inside this class. Further source channels, changed first outputs or different compiler geometry remain possible departures, each still requiring complete marked fibre, effective gluing and absolute atomicity.

## Scalar self-feedback: all source channels and the formal barrier

The [scalar self-feedback proof](../universal-keller-scalar-self-feedback/PROOF.md) replaces a root-channel restriction by a complete one-block affine coefficient-line classification. Let a_j=a_j0(u,sigma)+c_j(u,sigma)*g, with a nonzero effective direction. The actual response depends on g and satisfies a nonlinear first-order equation.

| Claim | Result | Boundary |
|---|---|---|
| UASELFLAUR1 | Laurent degree in r forces controls affine in r; reciprocal Riccati and divisibility classify all Laurent solutions, with explicit monomial-direction exceptions. | Characteristic-zero field/domain arguments. Frozen-response first integrals are not self-consistent coefficients. |
| UASELFCTRL1 | Any polynomial control using ALL own-source coordinates is Keller exactly when it uses only sigma and identity parameters. Even rational controls with whole polynomial output landing give only retained-signal effective profiles. | A rational signal control may cancel a common direction factor. The theorem concerns a single affine coefficient line, not curved families or mutually changing blocks. |
| UASELFFORMAL1 | Explicit rational scalar controls solve Keller exactly, have compatible polynomial reductions over every nilpotent parameter order, and preserve entire marked thickenings. Their sum has a source pole and no finite polynomial family in that line has the same source-dependent tangent. | The completed ring k[P][[lambda]] allows unbounded source degrees. Finite/formal exactness is not whole-source polynomial integration or full-map atomicity. |

The non-atomic canonical receiver inherits these nilpotent/formal auxiliary deformations with the same complete Spec(A) marked fibre, but retains its old factors. The next concrete integration gate is the source denominator 1+c_j*lambda*I(sigma)^ell*sigma^2*X^(2j-2). A curved/multi-direction or genuine compiler/auxiliary construction must remove it on the WHOLE source, retain complete marks and effective gluing, and prove full-map absolute atomicity. Localizing away from the denominator or checking further formal orders does not close that gate. See the [primary-literature comparison](../universal-keller-scalar-self-feedback/LITERATURE_AND_GAPS.md).

## Effective controls and the global polynomial boundary

The [feedback normal-form proof](../universal-keller-feedback-normal-form/PROOF.md) combines all-degree coefficient reduction with a precise account of formal source coordinates.

| Claim | Result | Boundary |
|---|---|---|
| UANATIVE2GEN1 | Every finite polynomial native list is EXACTLY two unique source controls A,B at indices4 and5. All syzygies, the basis wedge sigma^2*X^7 and a full two-control determinant equation are explicit. | Universal with2 invertible, including nonreduced bases. The control degrees, coordinate degree and covering degree remain unbounded. Arbitrary controls are not asserted Keller. |
| UAAFFINER1 | A whole native Keller map whose actual outputs are affine in r=1/t is identical to a retained-profile map, regardless of coefficient count. | One block with retained identities, characteristic zero and the stated native source jets. Nonlinear Laurent-r and mutually changing blocks remain outside. |
| UASIGFIELD1 | Every stable Keller endomorphism normalizing k(sigma) acts on the compiler triple by weighted scaling, with only a vertical Keller family left. Several individual fields allow only block permutation/scaling. | Sigma composed F in k(sigma) is NOT merely first output sigma. General vertical JC(m), the combined multi-signal field and arbitrary source coordinates moving the field are not classified. |
| UAETALEGAUGE1 | Classical formal etaleness makes every nilpotent deformation a unique source coordinate change, with compatible formal limits and marked-thickening preservation. For the explicit scalar family, source-gauge degrees necessarily escape every bound. | This creates no new formal equivalence class or global polynomial integration. Its tangent escapes signal-respecting polynomial coordinate orbits; unrestricted stable inequivalence is not proved. |

Actual retained-profile families with a_d=lambda*I(sigma)^ell can have any generic degree d>=4 and cubic special fibre while remaining formally source-coordinate equivalent to the cubic, fixing the marked thickening. Thus formal triviality supplies no global generic-degree or stable-class constraint.

The native-orbit extension below sharpens this target: Cout must be genuinely nonlinear in r and the pair must escape the generic profile-coordinate reduction, or the first-output/compiler/multi-block architecture must change. Removing the scalar formal denominator must regularize BOTH controls away from the basis degeneracy; merely adding a second response cannot hide it. Exact whole-source landing, complete canonical marks, effective gluing and full-map absolute atomicity remain required. The classical Riemann-Hurwitz/formal-etale inputs and novelty boundaries are recorded in the [literature audit](../universal-keller-feedback-normal-form/LITERATURE_AND_GAPS.md).

## Native profile orbits and an exact old-factor certificate

The [native-orbit proof](../universal-keller-native-orbits/PROOF.md) closes a larger coordinate sector and replaces an unbounded target-degree search by exact field reduction.

| Claim | General result | Boundary |
|---|---|---|
| UAQFSAT1 | For every dominant quasi-finite polynomial affine-space self-map, a rational target function with polynomial pullback is polynomial. | Elementary classical divisor/dimension consequence; non-quasi-finite maps and sources with nonconstant units need not satisfy it. |
| UANATIVEORBIT1 | A native Keller map has Cout affine in r exactly when it is a UNIQUE retained profile followed by b->b+a^5*c^3*h. Generic plane-coordinate affine reducibility and signal-field-respecting stable profile equivalence are equivalent descriptions. | Characteristic zero, one block, retained identities fixed by coordinates. Arbitrary source changes moving the individual signal field are outside. |
| UACUBICDESC1 | Reduce each output in the fixed cubic algebra E=a*S^3+b*S^2+S-c/2. Its S and S^2 coefficients vanish exactly when the output is polynomial in the old cubic outputs. | Fixed dimension three, unrestricted input coefficient sizes. Non-descent excludes only this old factor, not every intermediate map. |

Explicit whole polynomial native examples b->b+lambda*a^5*c^m*I(a)^ell exist for every m>=3 and ell>=1. They preserve the entire marked thickening and all full-receiver factors, so nonlinear-r controls and high marked contact do not by themselves produce a new stable class. The exact orbit landing ideal is (a^5*c^3); source x-order jets alone are weaker, as the non-native quadratic c^2 shear shows.

For outputs polynomial in r, leading coefficients permit a target shear whenever one positive r degree divides the other. Each step lowers the degree sum. Reaching degree1 closes the map into the classified orbit, and the ENTIRE quadratic polynomial-r native sector is identically retained-profile. A surviving degree pair must be non-dividing with both degrees >=2; negative Laurent-r sectors remain separate. These are necessary remaining gates, not atomicity certificates.

The current native candidate must have Cout genuinely nonlinear, survive this finite reduction, remove actual source poles and pass every relevant factor exclusion while preserving complete canonical Spec(A) and effective gluing. The scalar tangent fails old-cubic descent and escapes this restricted orbit, but no finite polynomial integration is constructed. Canonical absolute atomicity and unrestricted source equivalence remain open.

## Profile-factor rigidity and the canonical plane completion

The [plane-completion proof](../universal-keller-plane-completion/PROOF.md) reuses the published HC4MYGL0 line-injectivity theorem, already retained in the HC(4) research, to close a factor branch without further degree fitting.

| Claim | General result | Exact boundary |
|---|---|---|
| UAPROFFACT1 | Rational recovery through ANY retained profile forces the unique global target shear b->b+a^5*c^3*h. A genuinely nonlinear C excludes every such literal first factor. | One block with retained identities. Arbitrary intermediate maps and source changes moving the signal are not classified. |
| UARPLANE1 | Absence of negative r powers is exactly polynomial extension of the generic leaf to the (S,Y) plane, with Jacobian 1-SY and an embedded line. | Generic signal only; special-signal whole-source landing remains a separate requirement. Leading-divisor multiplicities further constrain non-dividing degree survivors. |
| UAPLANEATOM1 | An irreducible Jacobian plus an embedded line forces plane indecomposability; the geometric irreducibility gives the absolute statement. Every native polynomial-r completion satisfies it. | This is plane indecomposability, not stable atomicity of the full Keller map on the OPEN leaf. |

This removes the old UACUBICDESC1 conditional native target-factor-degree-greater-than-one branch: its target map fixes a line and must be invertible. Cubic remainder vanishing now characterizes exactly its degree-three target-shear orbit. Distinct retained profiles also have incomparable literal output fields.

The embedded-curve extension below removes the intermediate-extension hypothesis and closes EVERY signal-retaining whole factor, even with boundary poles and negative Laurent-r outputs. For arbitrary intermediate coordinates and stabilization mixing, a noninvertible target factor must instead have a genuinely NON-PLANE geometric generic retained fibre. No plane Jacobian conjecture is assumed, and the line theorem is applied only to a plane Keller factor, never directly to the nonconstant-Jacobian completion.

The [literature audit](../universal-keller-plane-completion/LITERATURE_AND_GAPS.md) also supplies explicit earlier public attribution for UAQFSAT1's characteristic-zero saturation overlap: Trushin Proposition3.1.1 proves it under the weaker absence-of-contracted-divisors hypothesis. Compiler-application novelty and independent proof review remain unverified.

## Embedded curves close the pole route; arbitrary factors require non-plane surfaces

The [intermediate-surfaces proof](../universal-keller-intermediate-surfaces/PROOF.md) replaces the preceding boundary-extension gate by a curve argument and an exact condition on EVERY stabilized intermediate factor.

| Claim | Result | Remaining boundary |
|---|---|---|
| UARELPLANE1 | A plane Keller target sending a polynomial curve to a closed embedding is invertible. In a whole relative factorization with this curve witness and two free generic coordinates, the target stage is invertible exactly when its geometric generic retained fibre is A2. | Published rectification and line injectivity are imports; a non-plane intermediate fibre is not classified. |
| UASIGFACT1 | Every whole native factorization retaining the signal and all identities has invertible target stage in ALL Laurent-r degrees. | This closes the old intermediate-pole alternative; it is relative atomicity, not full absolute stable atomicity. |
| UAINTERSURF1 | Any factorization after arbitrary identity stabilization and coordinate mixing has actual etale generic surface stages. A proper intermediate is NON-PLANE, smooth, geometrically integral, a complete intersection with free tangent bundle, a whole carried A1 component and geometric unit rank at most1. | Both maps satisfy exact saturated field-intersection identities and miss at most finitely many target points. These are necessary global-factor conditions, not an existence or exclusion theorem for all such surfaces. |

The source curve is i(s)=(u,0,s,a-4s^2), with exact target image (-s/2,0). Its intermediate image remains a polynomial closed embedding even if other intermediate coordinates have boundary poles. This makes a plane Keller target invertible without requiring extension to the completed source plane.

For arbitrary factors, let q=(E_u,E_sigma) on intermediate affine space. Its geometric generic fibre is the surface being tested; no global coordinate completion or plane-form theorem is assumed. Its coordinate ring B is determined by its actual embedded field as B=R intersection Frac(B), and the target plane ring is B intersection K0(b,c). A proper factor must make BOTH whole degrees greater than one. The source surface with the identity first stage is an allowed endpoint, not a proper factor.

The generic source surface itself admits etale self-maps of EVERY positive degree, but their canonical same-signal whole lifts have poles and volume ratio d*r^(2d-2). The open immersion into A2 misses a divisor and fails saturation. These controls prevent generic surface formulas from being mistaken for globally polynomial factors. Complete canonical Spec(A), whole-source landing, effective gluing and absolute stable atomicity remain open, especially for genuinely changing multi-block receivers.

## Primitive units close another entire intermediate branch

The [unit-rigidity proof](../universal-keller-unit-rigidity/PROOF.md) upgrades the preceding unit-rank bound to an exact proper-factor exclusion.

| Claim | Result | Boundary |
|---|---|---|
| UAPRIMUNIT1 | A cofinite dominant etale map from the full hyperbola complement to a smooth affine surface, pulling one target unit back to r primitively, is an isomorphism. | Generic Gm deck transformations extend and fix the retained node. Primitive index, cofinite image and the actual node are hypotheses. |
| UAUNITFACT1 | In ANY whole retained-signal Keller factorization, nonconstant geometric intermediate units occur exactly when the source stage is invertible. | No native jets or Laurent-degree restrictions. A proper native intermediate is unit-free AND non-plane; the latter still uses the curve witness. |

Published root closure is decisive: recovering a positive power of the polynomial t recovers t itself in the whole intermediate polynomial ring. Thus a nonconstant unit cannot hide in a proper finite-index lattice. The earlier generic J_d examples have index d and do not contradict the theorem.

The same criterion identifies the signal pattern. A reducible E_sigma zero fibre occurs exactly at an invertible source prefix, with two irreducible comaximal components pulling back to t and Q. Otherwise EVERY constant level E_sigma-c on whole intermediate affine space is geometrically irreducible. At successive cuts of a factor chain, the unit lattice drops from full Z to zero at the first noninvertible stage and never reappears; later unit-free chain length is not bounded.

The generic intermediate ring is a UFD over k(u,a) by localization. This is arithmetic factoriality. The log-Picard extension below now excludes geometric factoriality for a proper native intermediate. A proper unit-free non-plane surface still needs the exact embedded field, both saturation identities, carried curve, divisor/image conditions and whole-source polynomial realization. Such surfaces have not been excluded or constructed as proper factors. Canonical absolute stable atomicity, finite nonlinear integration, exact complete marks and effective gluing remain open.

## Logarithmic volume forces a geometric Picard jump

The [log-Picard proof](../universal-keller-log-picard/PROOF.md) closes two more general intermediate classes and identifies a sharper residual surface problem.

| Claim | Result | Boundary |
|---|---|---|
| UALOGVOL1 | In every actual whole retained-signal Keller factorization, the intermediate surface has log Kodaira dimension0 exactly when the source factor is invertible; otherwise it has dimension -infinity. | Log pullback is proved on boundary resolutions for nonfinite maps too. The whole Keller volume is r^2 times the source logarithmic volume. |
| UAGEOPIC1 | For native maps, geometric Picard group vanishes exactly when at least one stage is invertible. Every actual retained-signal intermediate has zero invariant geometric Picard subgroup over the retained function field. | The factor criterion imports Miyanishi's negative-log affine-plane characterization and needs the curve witness. Arithmetic UFD alone does not give geometric factoriality. |
| UADANVOL1 | No ordinary smooth Danielewski surface uv=P(z), degree(P)>=2, admits an etale map to A2. Higher powers u^m*v=P(z), m>=2, have explicit global regular volume primitives. | A nonzero compact volume period excludes exponent1. A primitive for higher exponents is not an etale pair or whole polynomial lift; u has no Keller mate for any exponent. |

Thus a proper native intermediate must have constant geometric units, negative log Kodaira dimension and NONZERO geometric Picard group with ZERO retained-field Galois invariants. A specified rational point from the actual retained-base polynomial section makes the line-bundle descent argument possible. The u^m*v=z^d-a models with m>=2 explicitly realize the arithmetic/geometric Picard jump and pass exact volume, but the two required arrows, free tangent frame, exact embedded field and whole-source realization remain unconstructed.

The raw spread a=z^d-u^m*v has critical points over a=0 for every exponent and therefore cannot itself be the whole retained tuple of a Keller factor, which is smooth everywhere. A different global model of the same generic surface is not excluded.

The finite torus-domination classification in Hajra2026 is context only; its finite-map hypothesis is not transferred to our arrows. The obsolete two-surface logarithmic-zero list is not used. Full canonical marked schemes, effective gluing, nonlinear native integration, independent proof/novelty review and absolute stable atomicity remain open.

## Plane charts expose a whole-curve and JC(2) barrier

The [plane-chart proof](../universal-keller-plane-charts/PROOF.md) strengthens the residual negative-log surface problem without a degree sweep.

| Claim | Result | Boundary |
|---|---|---|
| UAPLANECHART1 | A proper native intermediate's WHOLE carried affine line intersects the boundary of EVERY open A2 chart. All its finite/formal neighborhoods are nevertheless the standard plane neighborhoods. | Containment in one chart invokes published plane curve cancellation, then an open-and-closed section inverts the entire target. Formal standardization is classical etaleness, not global chart algebraization. |
| UADANCHART1 | Higher-power Danielewski targets require compatible nonautomorphic Keller pairs on their classical plane atlas, surviving at nondividing ruling degrees>=2. The carried curve meets at least two root components and satisfies deg Z>=m*N-r+1. | Any such target would disprove JC(2), a classical open-plane reduction explicitly documented by Dubouloz2009. A sharp beta-integral contact family fails geometric injectivity for EVERY exponent; other embedded curves are not classified. |

These are conditional requirements on any actual proper intermediate, not JC(2) counterexamples or a constructed etale atlas. The exact primitive from UADANVOL1 solves neither the simultaneous Laurent chart identities nor the whole carried-curve embedding. The raw polynomial spread is already excluded by its critical zero fibre. A different smooth whole model, both arrows, exact saturated fields, complete canonical marks, gluing and absolute stable atomicity remain open.

## Proved descent of affine-type geometric rulings

The [affine-type descent proof](../universal-keller-affine-descent/PROOF.md) addresses the actual Galois-permutation question.

| Claim | Result | Exact limit |
|---|---|---|
| UAMLFACTOR1 | Trivial geometric ML forces a geometric plane; for native/witness maps it is equivalent to whole target-factor invertibility. | Published Daigle classification plus the retained nonexact-volume obstruction; no classification for nonconstant ML is asserted. |
| UAAFFDESC1 | Every affine-type geometric ruling of a proper native intermediate is unique and descends. Centering at the rational point and weighted Galois trace construct the retained-field coordinate; contraction with the canonical volume constructs its additive action. | A proper target with such a ruling would disprove JC(2). Conditional on JC(2), the remaining intermediates are additively rigid and all rulings have COMPLETE base P1. Those Galois-permuted complete-type rulings remain open. |

An exact pair of conjugate polynomial coordinates has trace x+4*x^2*y and generic trace fibre Gm. Averaging different rulings is not a descent theorem. The earlier general ruling-descent gap below is now narrowed to complete type if JC(2) holds; no extra necessary-condition gate is added.

## Ruling descent, product topology and multiple-fibre budgets

The [ruling-descent proof](../universal-keller-ruling-descent/PROOF.md) generalizes the explicit Danielewski chart gate and adds constraints in every actual product-leaf rank.

| Claim | Result | Boundary |
|---|---|---|
| UAARITHRULING1 | ANY A1 ruling over the retained field has base A1, integral arithmetic fibres and reduced geometric fibres. Selecting components gives a finite geometric plane atlas, hence JC(2) counterexamples if the native target factor is proper. | Geometric ruling existence need not descend. Published Galois-invariant/unique-ruling descent DOES apply; conditional on JC(2), a surviving non-descent escape needs nontrivial Galois orbits of inequivalent rulings. The whole curve must meet two components of one split fibre. |
| UAPRODUCTTOPO1 | On an actual U^N times A^h generic leaf, both proper arrows are nonfinite; intermediate pi1 has an abelian subgroup with index DIVIDING source degree. Exact sheet losses obey a composition cocycle. | Source Euler characteristic is1 and pi1=Z^N. Product geometry is a hypothesis, not assumed for changed signals. In dimension2 a source finite completion needs a missing divisor, despite cofinite target image. Euler integrals of nonnegative losses need not be positive. |
| UAORBFIB1 | Every surjective connected-generic geometric fibration to A1 has at most one multiple fibre or exactly two doubles; two doubles require EVEN source-stage degree. | Classical orbifold quotient is imported explicitly. Descended one-block rulings satisfy the stronger reduced-component rule. Other reduced fibres and genuine non-descent geometries are not constructed or excluded. |

The remaining alternatives are an actual compatible JC(2)-counterexample atlas, or a proved ruling-descent failure with the full source, field, curve, volume and global polynomial constraints. A unique ruling is not an unresolved descent case. Complete canonical marks, whole-source landing, nonlinear integration, gluing, costs and absolute stable atomicity remain open.

## The construction pipeline

For a rank-N finite etale algebra A, the compiler starts with a primitive generator theta, rewrites its characteristic polynomial as
\[
 E(S)=S+bS^2+\pi S^3+\sum_{j=4}^N a_j\pi^jS^j-c/2,
\]
and reconstructs the entire Keller fibre by localizing at E'. A squarefree degree-N E gives the whole algebra A.

Five different questions arise:
\[
 \text{generator constraints}
 \ \longrightarrow\ \text{root-free scalar section}
 \ \longrightarrow\ \text{polynomial affine map}
 \ \longrightarrow\ \text{field-versality}
 \ \longrightarrow\ \text{canonical marked descent}.
\]
A result at one step does not prove the next step. Generic monodromy and stable atomicity are additional properties of the constructed map, not substitutes for universality or descent.

## Current universal dimension results

The invariant ktdim_k(S_N) refers to one fixed Keller self-map with a finite-etale fibre family universal over every extension field of k. The table states characteristic-zero results.

| Rank/base | Established upper bound or value | Source |
|---|---:|---|
| N=3, arbitrary k | <=3 | URKC1; full receiver proof also gives the original map in every infinite characteristic other than two |
| N=4, arbitrary k | <=4 | URKC1; a three-dimensional open-source quotient is not an affine A3 construction |
| N=5, arbitrary k | <=4 | UACOMP1 |
| N=6, arbitrary k | <=4 | UASEXT1 |
| N>=7, arbitrary k | <=N-2 | UATRACE1 |
| N>=7 with the binary trace-form test positive | <=N-3 | UADOUBLE1 and UABINARY1 |
| N=7 with sqrt(-14) in k | exactly4 | UAK7MIN1 plus Duncan's external lower bound |

Separately, UAFROB1 computes the positive-characteristic version exactly: ktdim_k(S3)=2 for every characteristic-three field, including finite fields. One fixed Frobenius plane map realizes every cubic etale algebra over every extension field, and UA1DISC1 excludes dimension one. This does not change any characteristic-zero row above.

For every characteristic-zero rank N>=3, UANONLINEARCOST1 now separates the retained fixed universal atomic upper bounds from canonical marked dimension, even for arbitrary nonlinear actions. The explicit fixed upper bounds are N in ranks3,4, four in ranks5,6, and N-2 thereafter; each is strictly below the canonical bound. Thus a dimension-minimizing fixed universal map cannot carry the required canonical marked symmetry at its own dimension. UAFIXDEFECT1 now gives a stronger conclusion for all retained generic-degree-N maps: no number of identity coordinates can make them canonical in characteristic zero. A possible canonical rank-N fibre must instead occur on a genuinely higher-degree map: at least N+3 in even rank and at least N+2 in odd rank. Affine canonical constructions in characteristic zero require at least N+3 in every rank.

For the same characteristic-zero rank-seven bases, UANONLINEARCOST1 now forces dimension at least eight for any natural marked canonical affine-space construction, even with arbitrary nonlinear actions. Thus the fixed minimum four and marked canonical cost differ by at least four; no canonical attainment at eight is asserted.

For a finite base F_Q, UAFINCOUNT1 gives a different, compiler-independent numerical bound: Q^d>=sum_(j=0)^(N-1)p(j). The required complete fibres share one finite source-point supply. Thus every universal S4 map over F2 needs at least three coordinates, and S7 over F2 at least five. UA4CHAR2 now supplies a globally finite atomic four-coordinate quartic map, leaving the S4/F2 minimum in {3,4}; no matching upper bound is asserted for S7/F2. This does not alter any characteristic-zero row. For the regular translation group (C_p)^e, both fixed framed universal and canonical affine-bundle minima are exactly ceil(e/f) over F_(p^f), and one over every infinite characteristic-p base.

The N-3 construction is field-versal over a quadratically closed characteristic-zero base exactly when N's binary expansion has at least three ones. Powers of two and sums of two distinct powers of two are negative cases for this construction. This is not a lower bound excluding a different N-3-dimensional compiler in those ranks.

The general essential-dimension lower bound applies to every field-versal Keller family; construction-specific stable-class defects and toric quotient failures do not supply global dimension lower bounds. Coordinate-degree minima remain a separate problem.

## The rules and their exact boundaries

| Rule | What it establishes | Boundary |
|---|---|---|
| Primitive generator plus discriminant unit | The fibre is an algebra, including disconnected etale algebras | A list of rational preimages is not a full-fibre proof |
| Inverse derivative localization | Repeated inverse roots are deleted | A degree-N polynomial with multiplicities does not give a rank-N etale fibre |
| Pi-zero formula for regular profiles | These complete fibres have rank at most three | Exceptional modular double poles instead have exact rank d or d+2; a full-degree boundary fibre necessarily splits |
| Positive primitive coefficient weight vector | A polynomial monomial scalar section exists | Gcd one is required for generic root-free transport; a coprime pair is sufficient but unnecessary |
| Trace zero | Removes one coefficient uniformly over infinite fields, with no division by N | In this cubic-normalized chart, rank four loses the required cubic coefficient |
| No further invariant linear generator subspace | Beyond trace, a natural homogeneous linear condition cannot keep primitive generators in codimension at least two | More compression must use nonlinear conditions or another receiver |
| First and second traces zero | Exactly the smooth trace-quadric realization criterion | An isotropic scalar radical point is insufficient |
| Rational trace quadric with a smooth point | Dense primitive generators with every denominator nonzero | A general quadric family need not have Zariski-local sections |
| Fixed-field quadratic-form bounds | Often realize every algebra over that one field | These bounds do not persist over arbitrary extensions |
| Binary-block test and odd-degree descent | Universality over all infinite extensions; necessity too when N is invertible | Modular ranks need the nonconstant-vector condition; no complete modular converse is claimed |
| Quadratic tilt with nonzero linear coefficient | One classical theorem gives S_N monodromy in every odd characteristic, for every coefficient-profile section with the top degree retained | It alone does not prove every torsor is a fibre |
| Primitive-monodromy argument | A sufficient condition for absolute stable atomicity | The converse fails in KMFI30; imprimitive monodromy needs an actual polynomial factorization argument |
| Sparse intrinsic Fitting endpoints | Stable equivalence is precisely residual coefficient scaling, even with intermediate zeros | The theorem is characteristic zero |
| Full root-line receiver | Canonical correspondence and polynomial gluing on rank-three vector bundles | It retains the generator line and does not select a section over BS_N |
| Inertia at most mu2 | No section over BS_N for N>=3 while retaining the good generator line | A different receiver discarding that line is not excluded |
| Coefficients replaced by functions of Pi | The rational Jacobian remains one; polynomial landing is classified for all rational profiles | No negative profiles in characteristic zero; modular exceptions require j=2 mod p and pole order at most two |
| Symmetric linear representation lift | Every subgroup-fixed restriction is Keller; a natural marked N-point fibre needs at least N+1 ambient coordinates in every characteristic | At equality a fixed plane would contradict JC(2) in characteristic zero; affine-linear moving components fail with any number of trivial directions |
| Affine canonical source, N>=5 | Dimension at least N when the characteristic divides N, and N+1 otherwise; the orbit barycentre suffices in the rank-unit case | Losing an origin saves at most one coordinate at high rank; arbitrary nonlinear actions remain separate |
| Retained commuting torus | A nontrivial commuting linear torus raises the characteristic-zero canonical linear bound to N+2 unconditionally | This is a marked linear architecture bound, not an unmarked universal dimension bound |
| Globally finite fixed map | Tame cyclic resolvents are constant on the affine target; no odd-characteristic fixed universal S_N map can be finite or proper | A nonfinite completion is allowed; a canonical family on a stack can retain Picard torsion |
| Global affine coarse quotient | Contracts an explicit A1 in every degree | The induced map is not quasi-finite despite its generic degree |
| Target c normalization | Gives a smooth source C=1 with a nonconstant unit | That actual level cannot be affine space, even after identity stabilization |
| Target B normalization with independent seeds | Exact strata and finite-residue counts exclude every level after affine stabilization | The sole affine exception is characteristic three with every retained index divisible by three; coefficient relations require a separate proof |
| Canonical affine-line classification | Exactly the natural groups F_q semidirect C_m with m dividing q-1, in characteristic p; polynomial is a subspace polynomial | The only natural S_N cases are S2 in characteristic two and S3 in characteristic three; composite q loses absolute atomicity in one variable |
| Fixed scalar-complement maps | One plane formula realizes every framed torsor for each F_q semidirect C_m over every extension | Minimum dimension is one for m=1 and two for m>1; the plane is absolutely stably atomic exactly when its multipliers generate F_q over F_p |
| Full-degree canonical symmetric fibre | Its geometric monodromy is normal in S_N and ordinarily primitive, so absolute stable atomicity is automatic | The genuine exception is rank four in characteristic two with group V4; a smaller boundary fibre of a higher-degree map is outside the theorem |
| Scalar-affine weighted deformation | Every actual canonical line model upgrades to a fixed universal plane over the same base; composition carries full subspace flags to polynomial atomic factors | Full F_q constants are unnecessary, but the primitive scalar roots and actual translation subspaces are required |
| Scalar-affine finite-base dimensions | Canonical minimum ceil(e/f); fixed framed upper bound ceil(e/f)+1, exact when f divides e and m>1 | Nondivisible-degree fixed minima retain a one-coordinate gap; coordinate-preserving maps always need the larger endpoint |
| Actual missing sheets at an affine divisor | Inertia embeds in S_(n-r_D); wild inertia needs loss at least p sheets | Rank loss alone is not ramification; omitted unramified branches can cause strict inequality |
| Matrix Lang construction | Canonical affine finite-field spaces in dimension n and fixed universal atomic maps in dimension n*n+n | Periodicity isomorphisms are data; constant AGL_n(F_q) groups require F_q constants |
| Characteristic-two quartics | Canonical affine dimension two; a fixed globally finite atomic universal map exists in dimension four | Canonical planes decompose; canonical absolute-atomicity dimension is 3..6, while fixed F2 universality is 3..4 |
| Finite-field whole-fibre capacity | Every universal rank-N map over F_Q needs Q^d>=sum_(j<N)p(j), by summing rational points over all splitting types | This applies to the finite base itself, not to an arbitrary finite subfield of an infinite base; no generic monodromy or global finiteness is required |
| Normal wild symmetry kernel | Regular deck groups are trivial in characteristic zero and p-groups in characteristic p; target first-jet kernels lie in O_p(G) | Natural S_N for N>=5 retains faithful target tangent symmetry in every characteristic; the first-jet kernel need not satisfy deck-degree divisibility |
| Tame Smith fixed loci | Whole prime-subgroup fixed schemes are smooth connected acyclic; dimensions zero and one have exact fibre-rank alternatives | Wild subgroups are excluded from this argument, and higher acyclic fixed varieties are not automatically affine spaces |
| Arbitrary nonlinear symmetric canonical cost | In characteristic zero: N+1 for prime powers, N for two positive same-prime powers, N-1 generally; every N>=3 three-space route is excluded | Equality in the strengthened cases requires actual acyclic fixed surfaces with a missing-fibre etale morphism; admissible tangent profiles in other ranks are not constructions |
| Nonlinear actions with affine fixed schemes | The whole fixed-locus restriction is again Keller, so empty fixed fibres need dimension at least two | A nonlinear formula is insufficient: published plane and tame three-space linearization transfer the earlier bounds; the known nonlinearizable S3 fixed-line action also fails |
| Even-rank target section | Gives explicit etale maps on affine open charts | The inverse divisor tq remains; filling it needs a polynomial extension proof |

Further uniform results: [all-odd-characteristic atomicity, exact modular boundaries and the equivariant normal form](../universal-keller-structure/PROOF.md).

Canonical references: [scalar compression](../universal-keller-compression/PROOF.md), [trace and boundary fibres](../universal-keller-trace/PROOF.md), [sextic section and sparse stable classification](../universal-keller-sextic/PROOF.md), [trace quadrics and binary criterion](../universal-keller-quadrics/PROOF.md), [full receivers and affine obstructions](PROOF.md).

## Canonical choices versus optimal dimension

The compiler's surviving three-dimensional stable-class parameter counts are:

| Family | Rank range | Parameters of distinct stable classes for each realizable A |
|---|---:|---:|
| Clean untraced coefficient family | N>=5 | N-4 |
| Trace-normalized family | N>=6 | N-5 |
| Double-trace family | N>=7 | N-6 |

These are dimensions of the displayed families of compiler slices. They are not dimensions of a global moduli space of all Keller maps, nor dimensions of the literal fibre of the marked receiver.

The canonical receiver's actual generator-line fibres have dimensions N-1, N-2 and N-3 respectively. Forgetting the three target coordinates accounts for the difference.

In particular rank seven can attain the global universal target-dimension minimum four over the specified bases while a fixed algebra still produces a one-parameter family of distinct three-dimensional compiler classes. Optimal universal dimension and canonical compiler descent are independent questions.

The full receiver has a canonical degree-N root cover, but any natural section over BS_N would need an automorphism image containing the full S_N. Its actual image is trivial or the single opposite-root involution. Adding decorations while retaining that line cannot fix the image.

## The separate costs of a canonical construction

UACANCOST1 gives an all-rank bound for a natural marked S_N fibre in a linear source: d>=N+1 in every characteristic. UAAFFCAN1 now extends the cost to every affine source bundle for N>=5: d>=N if the characteristic divides N, and d>=N+1 otherwise. A wild group order alone is insufficient to evade an origin, since a rank-unit orbit has a fixed barycentre. The point orbit costs N-1 directions beyond the full-group fixed space; the empty fixed fibre needs at least two fixed directions. In characteristic zero equality would produce a nonsurjective plane Keller map. A nontrivial commuting linear torus rules out equality unconditionally by UATORUSCOST1: it either makes that fixed-plane restriction an automorphism, or forces the moving component linear.

Thus at rank seven the known unmarked universal dimension is four over the specified bases, while a canonical marked linear construction requires at least eight coordinates, or nine with a retained commuting torus. Conditional on JC(2), nine is required even without the torus. These bounds do not assert existence at equality or change ktdim.

The modular canonical map Phi_q(x)=x-x^q descends on affine line bundles for every G_(q,m)=F_q semidirect C_m, m dividing q-1. UA1CANCLASS1 proves these are all possible faithful transitive one-dimensional canonical groups, geometrically: no other rank or finite group is missing from that category. The exact base-field criterion requires the multiplier roots and enough F_(p^d) translation constants, where d=ord_m(p); containing all of F_q is sufficient but unnecessary. The source has no origin even after vector-bundle stabilization, and its regular translation monodromy is absolutely atomic exactly for prime q. Over bases containing F_q, restricting the old full affine-group plane supplies the composite-degree canonical upper bound two. The following table retains that full-field constant hypothesis.

| Rank and complement | Canonical marked affine-bundle dimension | With absolute stable atomicity | One fixed framed universal map |
|---|---:|---:|---:|
| prime q, m=1 | 1 | 1 | 1 |
| prime q, m>1 | 1 | 1 | 2 |
| composite prime power q, m=1 | 1 | 2 | 1 |
| composite prime power q, m>1 | 1 | 2 | 2 |

The last column requires a specified G_(q,m) torsor on an open finite-etale whole-cover target chart and does not require absolute atomicity. UAAFFSUB1 gives K_(q,m)(u,x)=(u,x-u^((q-1)/m)*x^q), with exact arithmetic and geometric group G_(q,m). It is absolutely stably atomic precisely when the multipliers generate F_q over F_p. Otherwise a written telescoping polynomial identity supplies nontrivial global Keller factors. A different absolutely atomic framed universal map for those nonprimitive groups is not excluded. A canonical marked plane can instead retain the larger generic group H_q; these requirements must be distinguished.

UASUBDEFORM1 now replaces the binomial by the actual normalized subspace polynomial. If A_W(X)=sum alpha_j X^(B^j), B=p^ord_m(p), then
\[
 \mathcal D_m A_W(u,x)=\sum_j\alpha_j u^{(B^j-1)/m}x^{B^j}
\]
preserves composition and gives a fixed framed universal plane whenever W is an appropriate subspace of the actual base k. This removes the unnecessary full F_q-constant hypothesis from the fixed-map column wherever a canonical line model exists; it does not transfer the canonical atomic plane bound to weaker bases. Full F_B flags give a maximal sequence of degree-B absolutely stably atomic polynomial factors. The explicit finite normalization A_(W_i)(z_i)=t*v_i, u=t^m identifies one retained unramified branch and (q-1)/m omitted index-m branches; for m=1 every omitted branch is unramified. The deformed family is atomic exactly when dim_(F_B)W=1.

Over F_Q=F_(p^f) containing the scalar roots, UASCALARDIM1 handles insufficient constants by using r=ceil(e/f) translation coordinates. The canonical marked minimum is r, while the fixed framed minimum for m>1 is between max(2,ceil(log_Q(q+m-1))) and r+1. The exact value is r+1 when f divides e, and two when r=1; over every infinite base containing the scalar roots it is two. For example (u,x-u*x^4,y-u*y^4) realizes every framed natural degree-16 (C2)^4 semidirect C3 torsor over F4 and all extensions, with exact minimum three, while canonical dimension is two. A retained coordinate forces r+1 in every finite-base case, so it cannot close the remaining unrestricted gap when f does not divide e. See the [actual-base subspace theorem](../universal-keller-subspace/PROOF.md).

Over infinite bases containing F_q, these affine groups have classical essential dimension one. For m>1 their framed Keller minimum two therefore gives a strict gap, explained by the whole affine-space and constant-Jacobian requirements. Arithmetic primitive monodromy is also weaker than absolute atomicity; coefficient promotion and identity stabilization remain different operations.

UAFINITE1 gives another compiler-independent rule: every fixed universal S_N map in odd characteristic is nonproper, and its field-versal finite chart must omit a divisor. A globally finite family has a constant discriminant torsor on affine space. The canonical modular family can instead retain the tame quotient in the stack's torsion character line. See the [cost theorem and exact categories](../universal-keller-costs/PROOF.md).

Beyond linear presentations, UAFIXGEOM1 transfers the ambient bound through polynomial conjugacy: if r trivial identity directions linearize both actions, d+r>=N+1. In characteristic zero every finite plane action is linearizable, so canonical marked planes are excluded for every nontrivial faithful transitive group. Finite tame three-space actions are linearizable too, giving the earlier tame exclusion. UANONLINEARCOST1 now excludes every characteristic-zero natural S_N three-space construction for N>=3, even for nonlinearizable actions. The published nonlinearizable S3 action on A4 still has a whole fixed scheme A1 and fails the same fixed-fibre gate; adding moving representation directions does not change that fact. See the [symmetry classification and nonlinear gates](../universal-keller-symmetry/PROOF.md).

## Fixed tangent symmetry without linearization

UADECKWILD1 extends the regular deck argument to characteristic p: every regular polynomial source deck group is a p-group of order dividing generic degree. For a faithful source G action, the target kernel and its tangent-action kernel lie in the normal p-core O_p(G). High-rank S_N has trivial p-core, so a marked natural S_N fibre, N>=5, forces a faithful target tangent module in every characteristic even for nonlinear polynomial actions. The tangent-kernel order itself need not divide generic degree.

UASMITHFIX1 supplies actual geometry for every tame prime subgroup H: source and target fixed schemes are smooth connected mod-ell acyclic. In dimensions zero and one they are a point or A1. Characteristic-zero fixed-fibre rank must therefore be one in those dimensions; in characteristic p a line can also have positive rank divisible by p. This works with the entire marked fibre and does not require its rank to equal generic degree.

In characteristic zero, UANONLINEARCOST1 evaluates the resulting faithful tangent constraints. Arbitrary polynomial S_N actions cost at least N+1 in prime-power rank, N when N is a sum of two positive powers of the same prime, and N-1 otherwise. This closes all natural symmetric constructions in dimension three for N>=3, including nonlinearizable ones. At equality in either strengthened arithmetic case, a candidate needs an etale map between actual smooth mod-ell acyclic fixed surfaces missing the marked target. Only when both surfaces are A2 does this become a JC(2) counterexample. Outside those rank classes the ordinary standard tangent module passes every prime-subgroup test, so this method cannot itself improve N-1 or provide global gluing. See the [Smith and nonlinear cost proofs](../universal-keller-smith/PROOF.md).

## Generic degree and actual acyclic boundary geometry

UAFIXDEFECT1 separates generic degree M from complete marked fibre rank N. For every subgroup H with smooth connected whole fixed schemes, its generic fixed-restriction degree satisfies
\[
 |\Omega^H|\le d_H\le |\Omega^H|+M-N.
\]
Strict-henselian normalization proves this without flatness: omitted special-fibre scheme length may exceed generic rank, but its number of geometric points does not. A tame prime-power derangement gives acyclic fixed varieties and an empty marked fixed fibre. UAACYCLICRIG1 excludes a birational restriction, so characteristic-zero transitive canonical constructions need M>=N+2 in every dimension, including arbitrary nonlinear actions. All full-degree identity stabilizations remain excluded.

UAACYCLICGAL1 strengthens the geometry. Constant source units and torsion target Picard group make every Galois intermediate field come from an actual finite etale target cover. These hypotheses are automatic for smooth affine rationally acyclic varieties in characteristic zero. Quadratic maps therefore require a nontrivial double cover, equivalently Pic[2], of the target. Mod-2 acyclicity excludes that cover.

For natural S_N in even rank, disjoint transpositions are a 2-derangement. UAPARITYDEF1 therefore raises the necessary generic degree to N+3, with arbitrary polynomial actions. Every affine transitive characteristic-zero construction also requires N+3, in every rank, using its full-group fixed affine spaces. The affine symmetric surplus-three bound persists in every characteristic other than two outside (N,char)=(3,3). For arbitrary-action odd rank in characteristic zero, the retained endpoint N+2 requires nonzero Pic[2] in every relevant fixed target, killed by the quadratic restriction. No canonical attainment is claimed.

The abstract missing-fibre surface gate is attained. UAACYCLICMISS1 computes the full image of the published Dubouloz--Palka cover shears on S_k={u+u^2v=w^k}: at fixed degree k they omit exactly {(a,-1/a,0):epsilon+a Q(a)=0}. Arbitrary finite sets on that curve can be prescribed. The undeformed sample Q=0 is surjective, while the deformations have actual empty fibres. The sharp generic degree for non-surjective maps between mod-ell acyclic affine varieties of dimension at least two is two for odd ell and three for ell=2, geometrically in characteristic zero.

UAACYCLICATOM1 transfers primitive-geometric-monodromy stable atomicity to regular self-maps of the same acyclic affine variety, in every degree. A pluricanonical Jacobian makes every regular factor etale; birational rigidity then finishes the primitive case. Prime degree is a special case. The published shears in prime degree yield infinitely many atomic stable classes at fixed degree, distinguished by the number of connected components omitted after stabilization. Allowing arbitrary different intermediate varieties changes the category and the examples already factor through their finite cover.

UAACYCLICRECV1 supplies a positive canonical cyclic receiver in characteristic zero. Actual regular cyclic scaling makes the formulas descend on the line-bundle hypersurface u+u^2v=phi(w^k), whose complete marked fibre is exactly the input torsor. It includes BS2 over every characteristic-zero field and has minimal dimension two in this rationally acyclic receiver, with stable atomicity in prime degree. A fixed surface map is separately field-versal. The fixed scheme is A1 disjoint union Gm, explaining why rational acyclicity cannot replace mod-prime Smith acyclicity. Its coarse A2 quotient always has nonconstant Jacobian, so forgetting the cyclic cover is not an affine-space Keller construction.

The pseudo-plane fundamental groups persist under identity stabilization, so these sources never become affine spaces. What remains is their realization as actual whole tame fixed schemes, compatibility across all subgroups, effective ambient affine-space gluing and separate ambient generic-monodromy/atomicity proofs. See the [acyclic boundary proof](../universal-keller-acyclic-boundary/PROOF.md), [primary literature](../universal-keller-acyclic-boundary/LITERATURE_AND_GAPS.md) and [normalization budget](../universal-keller-fixed-fibres/PROOF.md).

## Finite self-receivers and canonical volume lifting

UAHOMOSELF1 sharpens the earlier homogeneous PAIR results: every finite self-cover of a characteristic-zero homogeneous variety has abelian regular monodromy. The central fundamental-group extension and Schur's finite-commutator theorem force this even with disconnected isotropy. Thus one fixed full-degree symmetric-universal SELF-MAP is impossible from rank3 onward in that category.

A finite free quotient of a torus is a larger geometry. UAINFRASELF1 proves that all its finite self-covers have an abelian regular normal subgroup in monodromy, using the intrinsic translation lattice, without assuming the map itself is monomial. Primitive monodromy is affine over F_p with irreducible point stabilizer; prime-power degree is necessary. Symmetric universality and canonical full symmetric rank N>=5 are excluded in every dimension and affine identity stabilization.

The symmetric ceiling4 is now attained. UANORMSELF1 constructs a fixed finite universal cubic volume self-map
\[
 (u,s,t)\mapsto(u/9,s^3-3s,t(s^2-1)/3)
 \quad\hbox{on }t^2=u(s^2-4),\ u\ne0.
\]
Its critical fibres include all disconnected cubic algebras. UANORMCYL1 classifies EVERY dominant surface endomorphism and every finite characteristic-zero identity-cylinder map. Fixed cubic universality and full canonical cubic symmetry occupy different exponent directions; no stabilization or finite coefficient coupling makes one map have both properties.

UAQUARTSELF1 constructs an actual fixed universally quartic SELF-MAP on a smooth affine FOURFOLD, over every field outside characteristics2,3. A free S3 base built from the cyclotomic norm-one torus and Gm packages the varying cubic-resolvent norm tori. Squaring the two-dimensional norm torus gives geometric S4 and strong absolute stable atomicity; z->z/2 makes a descended global volume invariant. Coprime norm-cohomology exponents make every input twist have a field point. This is fixed universality; a full canonical marked S4 action on that map remains open.

UAKANLIFT1 gives a separate all-degree functor in every characteristic. Every etale f:X->Y between smooth d-folds lifts naturally on their canonical covolume lines. The tautological d-form has a nowhere-zero exterior derivative, preserved by the lift. Degree, geometric monodromy and the ENTIRE marked zero-section fibre are unchanged, and all compatible group actions lift. Locally the map is (x,w)->(f(x),w/j_f(x)). This supplies volume preservation after receiver geometry exists.

These source spaces retain nonconstant units, and canonical-line lifting preserves unit and characteristic-zero pi1 obstructions. Higher-rank finite self-receivers must escape geometric torus quotients; the original whole affine-space canonical construction still needs its nonfinite higher-degree boundary. See the [proof](../universal-keller-torus-quotient-selfmaps/PROOF.md) and [primary literature](../universal-keller-torus-quotient-selfmaps/LITERATURE_AND_GAPS.md).

## Universal monodromy and all-rank homogeneous receivers

UAUNIVMON1 simplifies the atomicity question for every full-degree fixed universal map. A completely split base-field fibre makes arithmetic and geometric monodromy equal, since the rational-point section of the fundamental-group sequence acts trivially on every label. A full S_N fibre over some extension then makes both groups S_N. Thus full-degree universality already gives primitive monodromy in every characteristic. Use UAFINATOM1 for finite normal-affine factorization, or UAACYCLICATOM1 for nonfinite self-maps of the same acyclic affine variety. Smaller marked boundary rank N inside generic degree M>N does not determine degree-M monodromy.

UAQUOTRECV1 gives a classical special-linear quotient receiver in EVERY rank and characteristic:
\[
 SL_{N+1}/S_{N-1}\longrightarrow SL_{N+1}/S_N.
\]
Left S_N multiplication supplies canonical marked symmetry; effective twisting realizes the entire input algebra over schemes. Hilbert90 makes its fieldwise unmarked maps isomorphic to this one fixed universal cover, including over finite fields. Geometric S_N makes the finite cover strongly absolutely stably atomic. Source and target are different quotient varieties.

UAALTRECV1 replaces S_N by A_N in the quotient pair and supplies full canonical S_N symmetry by an explicit left-right action. Its primitive geometric monodromy is A_N, including the atomic quartic A4 case. Its individual covers have frozen sign resolvents and cannot be one fixed symmetric-universal map.

UAHOMOCOST1 makes the required geometric isotropy precise. For a characteristic-zero homogeneous target L/H, the image of pi1(L) in pi1(L/H) is central. A center-free monodromy group is therefore a quotient of pi0(H). This gives the following exact component minima in the finite affine homogeneous receiver-PAIR category:

| Rank | Canonical S_N | Canonical and strongly absolutely atomic | One fixed full-degree symmetric-universal cover |
|---|---:|---:|---:|
|2|1|1|1|
|3|1|1|6|
|4|1|12|24|
|N>=5|N!/2|N!/2|N!|

Tori attain the low-rank1 entries, alternating quotients the12 and half-factorial entries, and symmetric quotients the factorial entries. These count geometric stabilizer components, not affine-space coordinates or dimension.

UAFINPI1RIG1 explains the remaining bridge: a variety with finite geometric pi1 has no finite etale endomorphism of degree greater than one, and this persists under affine identity stabilization in characteristic zero. The quotient sources and targets have different finite pi1 orders and cannot become isomorphic this way. Thus the finite nonabelian receiver-pair existence gate is closed; the original whole affine-space higher-degree boundary remains open. See the [canonical proofs](../universal-keller-homogeneous-receivers/PROOF.md) and [primary literature](../universal-keller-homogeneous-receivers/LITERATURE_AND_GAPS.md).

## Canonical torus receivers and stable factorization

UATORUSCLASS1 gives an all-rank criterion: a faithful transitive label group G has a marked torus-torsor receiver exactly when the characteristic does not divide the rank N and its action normalizes an abelian regular subgroup of Sym(N). That subgroup need not lie in G. The equivariant permutation-character lattice construction gives dimension at most N-1 and effective torsor gluing. This is the abelian Hopf--Galois normalizer condition applied to this receiver; source and target torus forms are kept distinct unless identified.

For natural symmetric labels the complete classification is:

| Rank | Canonical torus minimum | Excluded characteristic | Absolute atomicity | One fixed symmetric-universal torus map |
|---|---:|---|---|---|
|2|1|2|yes|yes|
|3|1|3|yes|no|
|4|2|2|no|no|
|N>=5|no realization|all|not applicable|no|

The cubic map is cubing on the norm-one torus a^2-ab+b^2=1, with rational source S3 action generated by translation by theta and inversion, where theta^2+theta+1=0. It works over every field of characteristic different from three, including characteristic two. The quartic map squares the two-dimensional product-one torus, with V4 translations and S3 coordinate permutations. Twisting gives exact complete cubic or quartic input-algebra fibres. The fieldwise source torsors are trivial by the coprime norm-cohomology exponents, but this supplies no canonical origin and does not fix the input-dependent torus form.

UATORUSFIX1 proves the separate fixed-map obstruction: cubic sign-resolvent torsors and quartic cubic resolvents are frozen for one fixed torus cover. Extending the base by a variable supplies missing fibre types even when the original field has few extensions. In characteristic zero the obstruction persists for every finite etale endomorphism of a torus times affine space, including arbitrary finite coefficient coupling. A fixed framed cover for a finite abelian group does exist from a quasi-trivial source torus; it may map to a different torus form. Constant C3 already has a one-dimensional fixed universal self-cover on the cyclotomic norm-one torus.

UAFINATOM1 makes every dominant affine factor of a finite map finite; normal intermediates then give finite etale factors and strong prime/primitive atomicity. UATORUSSTABLE1 proves that finite etale torus-cylinder endomorphisms in characteristic zero are torus isogenies times identity up to source automorphism. Stable self-map atomicity reduces exactly to factorization in the Galois endomorphism ring of the unit lattice. In the minimal quartic family a non-split cubic resolvent gives base-field stable self-map atomicity, and a cubic-field resolvent gives the stronger normal-affine atomicity. Every quartic member decomposes after extending constants. Nonisomorphic torus forms cannot become isomorphic by adding affine identity coordinates.

UAPSEUDOSYM1 closes a different route: finite automorphism groups of the original S_k surface are cyclic in characteristic zero, so it cannot support natural S_N for N>=3. A stabilization preserving that surface projection still fails. Generic cyclic degree N>=4 blocks every stabilized full-degree natural S_N lift even with mixing, by the finite-chart normalizer rule. A cubic stabilization that genuinely mixes the surface projection remains open. See the [torus and stable-atomicity proof](../universal-keller-torus-receivers/PROOF.md) and [primary literature](../universal-keller-torus-receivers/LITERATURE_AND_GAPS.md).

## Canonical symmetry, periodic receivers and the quartic exception

If the natural N-point canonical fibre has the full generic degree N, UACANATOM1 gives an invariant finite chart at the marked target. Its full symmetric label action normalizes geometric monodromy, forcing a transitive normal subgroup of S_N. The action is primitive, and the map absolutely stably atomic, except possibly for N=4 in characteristic two with group V4. UAFIXDEFECT1 now makes this implication vacuous in characteristic zero: no such full-degree canonical construction exists in any dimension. The modular implication remains valid. For characteristic-zero boundary candidates M>=N+2 it gives no automatic atomicity, so their generic monodromy and atomicity need separate proofs.

The matrix Lang map (U,X)->(U,X-U*X^[q]) has Jacobian one, full AGL_n(F_q) monodromy and framed every-field universality on det U!=0. The vector map X->X-X^[q] instead descends canonically on affine bundles in dimension n. In basis-free form the receiver consists of an affine bundle with specified compatible Frobenius-periodicity isomorphisms; its difference map glues polynomially and has derivative identity. Zhang's September2026 stack formulation extends this correspondence to arbitrary fibered bases and imperfect fields. It does not choose periodicity data merely from an existence property.

For quartics in characteristic two, S4=AGL2(F2) gives canonical affine-plane dimension exactly two. Every such affine-plane model factors through the two quadratic V4 translation quotients, so its canonical absolute-atomicity minimum is between three and six. A different fixed universal map is explicit:
\[
 (s,t,x,y)\mapsto
 (s,t,x+(1+st)x^2+s y^2,y+t x^2+y^2).
\]
It is globally finite etale of rank four, atomic, and realizes every quartic algebra over every extension field. Two root-free Frobenius matrix charts prove versality even over F2. This leaves the fixed F2 minimum between three and four. The fieldwise normal form is not asserted to glue canonically.

This is also a useful framing control. In characteristic two an A4-framed universal map has a tame C3 quotient and must retain a rank-one ramified divisor by UAINERTIACOST1. The finite S4 map nevertheless realizes every individual A4-induced quartic algebra; it has no global A4 frame. A structured generic family and individual subgroup fibres remain different requirements. See the [affine and Lang proofs](../universal-keller-affine-lang/PROOF.md).

## How the literature changes the next attack

The established generator-scheme and twisting framework removes the need to guess individual generator changes. Brassil--Reichstein's quadratic trace criterion classifies the current N-3 construction, while Duncan supplies an actual matching global lower bound in rank seven.

Cohen's published two-indeterminate-coefficient theorem, together with its correction preserving Theorem1, supplies an older but decisive transfer: arbitrary compiler profiles have full symmetric monodromy and polynomial sections have absolute stable atomicity in every characteristic other than two. The exceptional modular boundary is now completely described by UAMODBOUND1. See the [structural literature map](../universal-keller-structure/LITERATURE_AND_GAPS.md).

The 2026 AI-assisted septic-covariant work offers a different nonlinear constraint, using first and third traces. Its direct rank-seven compiler support has weights six and eight, so a mu2 scalar obstruction remains. The field-extension theorem is proved in that paper; the announced wider etale theorem is not treated here as a completed input.

Quotient-stack theory identifies a genuine canonical correspondence and gives an immediate inertia test for any proposed section. Good-moduli theory also explains why one must check the global coarse morphism, rather than infer its etaleness from the stack morphism. Here an explicit contracted line supplies the failure certificate.

The [coordinate-level and Frobenius extension](../universal-keller-levels/PROOF.md) adds a finite-residue obstruction that ordinary Euler characteristic misses, and an effective modular canonical construction on an affine source bundle. Its lack of origin persists under vector-bundle stabilization; fieldwise trivializations exist but are choices. The [cost extension](../universal-keller-costs/LITERATURE_AND_GAPS.md) transfers Shaska's published 2026 all-signature planar rigidity theorem to the N+2 commuting-torus bound, and uses classical Kummer and tame covering theory to distinguish finite universal maps from canonical stack families. Its one-dimensional canonical Frobenius simplification leaves the fixed-map minimum two for cubic universality in characteristic three.

Zhang's September2026 gerbe formulation of Lange--Stuhler, together with Morales--Sanchez's earlier Frobenius-module construction, supplies the periodic receiver and matrix-frame method; the [Lang literature map](../universal-keller-affine-lang/LITERATURE_AND_GAPS.md) records the exact inputs. Lamy's September2026 survey makes the plane and tame three-variable linearization inputs explicit, and the Freudenburg--Moser-Jauslin action supplies a real nonlinearizable comparison. Kemper--Mattig's older generic polynomials explain why one generic parameter does not give one Keller coordinate. See the [new primary literature map](../universal-keller-symmetry/LITERATURE_AND_GAPS.md), [quadratic map](../universal-keller-quadrics/LITERATURE_AND_GAPS.md) and [receiver map](LITERATURE_AND_GAPS.md).

Haution's2017 Proposition3.2.1 supplies the tame Smith geometry needed to remove the linearizability assumption; the [Smith literature map](../universal-keller-smith/LITERATURE_AND_GAPS.md) records that classical input, the elementary all-rank branching proof and the exact surface limitation.

The retained [KMFI30 proof](../monodromy-theory/SIGNED_ALTERNATING.md) gives actual imprimitive absolutely stably atomic maps. Thus the explicit scalar-complement factorization is a positive certificate, not a presumed converse of primitive-monodromy atomicity.

## Rigid-base reduction and the closed quartic action gate

The [15 November structural extension](../universal-keller-rigid-base-reduction/PROOF.md) computes the full geometric automorphism group of the constructed quartic W as a metabelian group. It cannot contain S4. A1-rigidity forces every cylinder automorphism to project to the base; every finite coupled map reduces to a base map after a source automorphism. Thus the particular W canonical S4 gate is closed on every finite coupled cylinder, while its fixed universality and strong absolute stable atomicity remain.

Every dominant regular map between equal-dimensional geometric torus quotients is finite etale. On any rigid metabelian cylinder, a remaining nonfinite natural S_N datum, N>=4, restricts to a genuine canonical polynomial Keller datum on one affine fibre. The affine coordinate count must itself satisfy r>=c(N). For W quartics this means at least five affine coordinates, total dimension at least nine, and generic degree at least seven times the base degree, with no existence assertion at those bounds.

Whole Sylow subgroups remove the classification-dependent derangement input: every characteristic-zero transitive polynomial label action needs M>=N+2, and every even-degree transitive action needs M>=N+3. In characteristic p the +2 result holds for every transitive group outside p-power rank. See the [source and gap map](../universal-keller-rigid-base-reduction/LITERATURE_AND_GAPS.md).

## All-degree S3 norm families and the stronger canonical boundary

The [holonomy-family proof](../universal-keller-holonomy-families/PROOF.md) extends the universal quartic to degree-n-squared volume self-covers on every Z_m. Specified ((Z/n)^2 semidirect S3) torsors are field-versal exactly away from multiples of three, with mu_n for the constant-group formulation. Every prime p other than three gives a strongly absolutely stably atomic degree-p-squared self-map. These are not arbitrary rank-p-squared universal algebras for p>2.

All dominant endomorphisms of Z_m are classified by one d=1 mod3 exponent, one odd exponent and an integral multiplicity matrix. The additive self-factor length counts prime factors in d,e and det(M). Determinant three gives a degree-nine absolute stable self-atom that factors through a different smooth quotient by two degree-three maps; its intermediate is not affine-stably isomorphic to Z_m.

Repeating standard S3 modules never repairs the finite canonical quartic, even though Z_3 has an actual faithful S4 automorphism subgroup. The free order-three translation still excludes the required marked fibre action. The general projection-block theorem gives a primitive embedding-or-affine-fibre alternative on any A1-rigid base, and a commutative holonomy commutant makes the outer-automorphism kernel metabelian. Use these tests before changing receiver geometry.

## Arbitrary-group free bases, exact isogeny criteria and the scheme boundary

The [general versal-isogeny proof](../universal-keller-versal-isogenies/PROOF.md) constructs a free field-versal torus base for every finite tame group by induced characters. An equivariant self-isogeny then gives a fixed quotient self-cover exactly when its map on every twisted H1 has zero kernel. Coprime degree and a uniform cohomology annihilator suffice. Norm tori have the sharp gcd(n,r) criterion for every transitive group, and signed permutation tori give strong atomic covers in every odd prime-power degree. The norm symmetric family gives every even binary exponent; a cyclotomic Gaussian-period isogeny adds degree8 with AGL1(F8) monodromy. These are specified-group torsors, preserving the full symmetric ceiling4 for torus-quotient self-maps.

Cyclotomic self-atoms correspond to irreducible integers; strong atoms correspond to maximal principal ideals. A nonprincipal ideal changes the quotient geometry in a way no affine identity stabilization can hide. Constant-group claims retain kernel-splitting hypotheses; volume is supplied by the existing canonical-line functor.

The [new literature boundary](../universal-keller-versal-isogenies/LITERATURE_AND_GAPS.md) is crucial for promotion. First--Florence--Rosengarten exclude one torsor over a quasi-compact scheme realizing all smooth affine finite-group families in characteristic zero, with an explicit dimension bound. Pointed base changes, including identity stabilization and canonical-line lifting, preserve the exact realization image. This does not exclude stack receivers retaining inertia or unbounded/dimension-dependent models. Scavia2026 disproves general versality detection by successful Sylow restrictions; negative Sylow degree/norm obstructions remain valid. Kresch--Tschinkel2026 give a legitimate generic-twist test with an actual versal action in their stated setting.

## Exact prime-power degree spectrum through semilinear holonomy

The [semilinear spectrum proof](../universal-keller-semilinear-spectrum/PROOF.md) closes the remaining degree-existence gate. In characteristic zero the nontrivial strong absolute stable atomic degrees of finite geometric torus-quotient self-maps are exactly the prime powers. Every prime power has a volume-preserving, field-versal specified-group version with source dimension at most max(2,2N-3). This preserves the full symmetric self-universality ceiling4 and the nonconstant-unit obstruction to affine-space landing.

A holonomy automorphism changes the ideal obstruction to [P]=(sigma-1)[I]. At prime conductor the ordinary ambiguous-class formula makes sigma-1 bijective, without a class-number-one assumption. Zsigmondy and two direct exceptional constructions then give every degree; an explicit nonprincipal degree2048 integral lattice verifies the new mechanism. The full cyclic free base carries the required holonomy automorphism even when a smaller field base cannot.

Puncturing the canonical covolume line keeps the volume lift in torus-quotient geometry and preserves field realization images. A nonzero covolume at a fixed marked target survives the canonical action exactly when the tangent determinant character is trivial. Prescribed affine groups, sharp dimension bounds and the original coherent whole-affine canonical receiver remain separate open gates.

## Remaining programme, without reviving closed routes

1. Construct a canonical ATOMIC whole-affine-fibre receiver, improve the trace-filter degree/dimension bounds, or prove an obstruction to that stronger architecture. Non-atomic canonical polynomial-bundle existence, same-source/target identification, compatible natural symmetry and simultaneous fixed-locus gluing are now constructed by UACANFIB1 and UACANFIX1. Its explicit two-factor tower cannot be removed by coordinates or identity stabilization. A natural global basis or one fixed all-scheme receiver remains a separate stronger demand. Keep finite quotient limits, smaller marked rank versus generic degree, modular affine minima and pseudo-plane projection mixing within their exact scopes.
2. Find a different affine construction for negative binary-test ranks or bases. The current double-trace compiler is already classified when N is invertible.
3. Decide low-rank A3 universality, especially ranks five and six, and arbitrary-base rank-seven minima. Existing lower-dimensional receivers and open-source maps do not close these questions. For finite bases apply UAFINCOUNT1 first: S4/F2 now has bounds three and four by UA4CHAR2, while S7/F2 has lower bound five without a matching construction. A three-coordinate universal quartic map or a stronger lower bound would close the first gap.
4. Prove a new affine completion or source normalization. UAPROFILE1 classifies rational coefficient profiles and excludes cancellation of the even-rank chart poles in the original source. The global coarse quotient is excluded, and UABLEVEL1 completes the coordinate-level classification for independent promoted seeds: in characteristic zero no Pi, B or C level is stably affine space. Coefficient relations, nonlinear sections and other completions still require new source geometry. UAFINITE1 requires a nonfinite completion for fixed symmetric universality in odd characteristic; filling the chart as a globally finite cover is excluded in every dimension. UAINERTIACOST1 further forces a rank-one ramified divisor for every natural cyclic-affine generic group with nontrivial complement. The finite characteristic-two S4 map lies outside that tame-quotient obstruction.
5. Obtain global coordinate-degree minima or functorial lower bounds applying to every compiler, rather than reuse vertical stable-class counts.
6. Independently review the written mathematical proofs. Exact replay controls are useful evidence but do not formalize the stack, monodromy, or imported arithmetic theorems.
7. For scalar-affine framed universality over F_(p^f) containing the scalar roots, decide the remaining r versus r+1 interval when r=ceil(e/f)>1 and f does not divide e. UASCALARDIM1 settles the divisible-degree family, and excludes all coordinate-preserving attempts at the smaller endpoint. This modular framed gate is separate from the characteristic-zero and symmetric-group problems.

## Atomic cylinder reduction, with categories made explicit

The [atomic-cylinder proof](../universal-keller-atomic-cylinder/PROOF.md) strengthens the general boundary without another degree search.

| Claim | Result | Boundary |
|---|---|---|
| UANORMAFF1 | A nontrivial etale map is absolutely stably atomic through normal affine intermediates exactly when finite with primitive geometric monodromy. | Nonfinite completions give a nontrivial regular factorization; the completion can ramify. |
| UASELFATOM1 | Constant units and torsion canonical bundle suffice for primitive self-map atomicity, without acyclicity, in every characteristic. | The classical Ax self-map lemma concerns the same variety; it does not identify different acyclic varieties. |
| UACYLDICH1 | A nonfinite self-atomic map on any characteristic-zero torus-quotient cylinder has base degree one. | Both base and vertical degree greater than one is an explicit self-factorization. |
| UACYLCAN1 | Its complete marked canonical data restrict to one actual affine fibre, with no metabelian base hypothesis. | Relative restriction gives a pair of affine-space forms, without canonical coordinates or source-target identification. |
| UACYLDEF1 | On the fixed Gm times A3, nonfinite etale degrees are exactly every N>=3, with infinitely many stably inequivalent atomic volume maps in each degree; finite atomic degrees there are primes. | These self-map atoms create no new transitive symmetry or affine-space source; nonfinite normal-affine atoms remain impossible. |

For natural symmetric data the affine coordinate count r itself pays the retained dimension cost. The degree budget remains M>=N+2+indicator(2 divides N). Atomicity on a whole cylinder does not assert atomicity of every special affine fibre. The remaining ATOMIC canonical task is a coupling outside UAPARALLEL1, UAJOINTFACT1, the descended class UABASICTIME1 and the root-polynomial class UAROOTDAG1, with isolated scalar stages subject to UASELFCTRL1, literal retained-profile first factors excluded by UAPROFFACT1, ALL signal-retaining factors excluded by UASIGFACT1, and arbitrary stabilized factors subject to the geometrically unit-free non-plane intermediate surface conditions of UAINTERSURF1, UAUNITFACT1, UALOGVOL1 and UAGEOPIC1: negative log dimension and nonzero geometric Picard group with no retained-field invariant class, with the WHOLE carried curve avoiding containment in every plane chart by UAPLANECHART1 and any Danielewski atlas subject to UADANCHART1, retaining its exact marked fibre and effective gluing; UACANFIB1 already supplies the non-atomic polynomial-bundle version. Moving degree into a torus base cannot solve the atomic problem by itself.

The goal remains active. This map records the resulting theory and the remaining requirements; it does not declare all open problems solved.

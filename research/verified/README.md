# Verified core

This directory contains the stable proof chain:

- [Twelve-variable degree-three reduction](TWELVE_VARIABLE_DEGREE_THREE_KELLER_COUNTEREXAMPLE.md):
  an exact source-graph/target-square reduction of MacFarlane's map, its
  independent sparse determinant replay, and its 19-variable
  cubic-homogeneous parent;
- [Foundational Keller map](FOUNDATIONAL_GEOMETRY.md): exact determinant and collision;
- [Characteristic-two threefold, Mondello plane theorem, and field extension](HUQ_KURUVILLA_CHARACTERISTIC_TWO_AUDIT.md):
  the repository's Huq--Kuruvilla wild-radicial normalization; Mondello's
  external arXiv:2608.02634v1 plane theorem, hidden cubic, actual-target
  irreducibility proof, and polynomial skew-product coordinates; and the
  separately named repository corollary extending the theorem from
  `algebraicClosure(F_2)` to every characteristic-two field;
- [Characteristic-two plane normalization and wild boundary](HUQ_KURUVILLA_PLANE_BOUNDARY_NORMALIZATION.md):
  the finite normalization, primitive-order conductor, unique missing
  Frobenius line, complete generic boundary ledger, and `S_3` monodromy of
  the plane map; this is a repository theorem whose independent second
  normalization implementation remains open;
- [Characteristic-two plane modulo-four obstruction](HUQ_KURUVILLA_PLANE_W2_OBSTRUCTION.md):
  the full de Rham cokernel and nonzero Cartier obstruction proving that no
  polynomially left--right equivalent plane representative has a
  constant-Jacobian lift over `Z/4`, together with the general
  one-coordinate first-order stabilization theorem and an explicit
  compatible polynomial tower over all finite Witt levels; the first stable
  Witt lift has sharp minimum degree 18, and an explicit `W_3` correction has
  degree 25, sharply among extensions of the canonical first digit; an exact
  Boolean certificate excludes unrestricted degree 18 at `W_3`, while a
  directly replayed SAT witness attains degree 19, proving `d_3=19`; the
  preferred frozen determinant-one witness has 440 correction coefficients;
  its fixed-first second-completion space has dimension 804 and exact minimum
  support 160 by incidence-component decoding; at `W_4`, that preferred
  degree-19 representative is obstructed in every degree by a nonzero
  `H^3_dR` class, while the fixed degree-25 representative has exact extension
  degree 52 by a two-coefficient functional and an explicit `8zL` correction;
  the unrestricted next-class-zero search inside the degree-19 ansatz is SAT,
  producing an 818-term representative whose own exact fixed extension degree
  is 52; a joint degree-19 master--subproblem experiment finds only singleton
  codomain-hole obstructions through 64 cuts, while a 599-hole master times
  out rather than proving UNSAT; a shared-minor compiler now reaches the full
  4,340-hole quotient, which also times out, with `z1` and `z3` separately
  SAT and `z0`,`z2` undecided; unrestricted, only `19 <= d_4 <= 52` is proved;

<!-- status-consumer: HKM2W1 904c57385ac0b0dd -->

<!-- status-consumer: HKM2W2 474e0d677133ee23 -->

<!-- status-consumer: HKM2W3 55d99efeee0298af -->

<!-- status-consumer: HKM2W4 6075dd4fbb9cb89c -->

- [Tangent-map core](TANGENT_MAP_CORE.md): the central theorem unifying the
  inverse pencil, plane incidence, Jacobian factor, discriminant normalization,
  reconstruction pole, Hessian Fitting divisor, weighted suspension, and the
  comparison with cancellation suspension;
- [Normalized factorization model](NORMALIZED_FACTORIZATION_MODEL.md): three compact
  propositions giving the polynomial `A^3` source, coefficient--resultant
  étaleness, exact relation to the foundational polynomial, LND/slice proof,
  and unequal-degree extension;
- [Foundational incidence construction](FOUNDATIONAL_INCIDENCE_CONSTRUCTION.md):
  projective normalization, arbitrary hyperplanes, the three contact orbits,
  and the exceptional `(2,1)` affine slice;
- [Cubic marked-root model](MARKED_ROOT_MODEL.md): the marked-root isomorphism
  and its affine-root and root-at-infinity reconstruction charts;
- [Cubic image and nonproperness theorem](IMAGE_AND_NONPROPERNESS.md): exact image, fibers, and nonproperness;
- [Weighted marked-root theorem](WEIGHTED_SEED_THEOREM.md): weighted construction and symmetric
  monodromy.
- [Universal symmetric monodromy](UNIVERSAL_SYMMETRIC_MONODROMY.md): the
  standalone classical proof that every characteristic-zero pencil
  `H(W)-sW+t`, including all exceptional polynomial types, has geometric and
  arithmetic monodromy `S_n`.
- [All-degree rational fibers](ALL_DEGREE_RATIONAL_FIBERS.md): explicit
  integer-root seeds giving a complete regular `N`-point rational fiber and
  `N` nearby real sheets for every `N>=3`.
- [Finite étale Keller fibers](FINITE_ETALE_KELLER_FIBERS.md): every finite
  étale algebra over a characteristic-zero field occurs as a full Keller fiber
  unless its rank is two. Ranks `N>=3` are realized explicitly in `A^3` by
  Jacobian-one maps of coordinate degree at most `6N+2`, compatibly with scalar
  extension. The note includes the degree-two descent, scheme-theoretic
  quotient reconstruction, the collision algebra `A tensor_K A` and its
  rank-`N(N-1)` off-diagonal obstruction, cubic `S_3` normal-closure sheets,
  higher ordered-root configurations reaching the splitting field,
  the exact optimal-quintic decomposition
  `A_5 tensor A_5 = A_5 times (N_6^3 times L_2)`, the updated arithmetic
  chain, and a staged Lean certificate.
- [Universal relative Keller map](UNIVERSAL_RELATIVE_KELLER_MAP.md): packages
  all supplied presentations into one relative Jacobian-one map, compresses
  the map base sharply from `N+1` to `N-3` parameters by moving three inverse
  coefficients into the target, and identifies the universal root fiber over
  the `N`-dimensional incidence open.  Its ordered collision cover splits into
  diagonal rank `N` and off-diagonal rank `N(N-1)`; generic off-diagonal
  monodromy is the `S_N` action on ordered distinct pairs, and the higher
  ordered-configuration tower reaches the full splitting field.  Promoting
  the remaining parameters produces one `Q`-defined map of `A^N` for every
  `N>=3`.  That fixed map is universal for rank-`N` finite-etale fibers over
  characteristic-zero fields, has `S_N` monodromy, and is absolutely and
  stably atomic.  The note also separates the canonical collision operation
  on `BS_N` from the obstructed Tschirnhaus descent of the Keller lift and
  from essential-dimension questions.
- [Generic Tschirnhaus non-descent](GENERIC_TSCHIRNHAUS_NON_DESCENT.md):
  puts the compiler on the clean primitive-presentation groupoid, proves
  that equal stable boundary has codimension `N-4` while projective
  transport has codimension `N-3`, and concludes that a generic
  nonprojective Tschirnhaus change preserves the root algebra but changes
  the intrinsic quadratic-gauge Fitting fingerprint.  The split witness
  `r_i=i`, `u_i=i+i^2` is separated in every rank by one symbolic formula.
  The accompanying all-rank stabilizer theorem exposes the universal
  discriminant intruder `P^2*B^N*C` and proves that arbitrary identity
  stabilization cannot move a physical marked target of any fixed clean
  map for `N>=5`.
  The image of the clean fixed-algebra marked locus in unmarked stable-map
  moduli has geometric dimension `N-4`; this is not the dimension of a
  literal marked receiver fibre, and no total global moduli-fibre dimension
  is claimed.
- [Rank-five ambient and marked Tschirnhaus transition loci](RANK_FIVE_TSCHIRNHAUS_TRANSITION_LOCUS.md):
  computes the ambient and projective loci and their transverse
  intersection, writes the unique canonical seed equivalence explicitly,
  and proves that it carries the selected complete fibre only on the
  root-scaling locus.  Every remaining marked nonprojective lift is reduced
  to a target-orbit problem for the stable self-equivalence group of one
  fixed quintic map.
- [Fixed-quintic stable target stabilizer](RANK_FIVE_STABLE_TARGET_STABILIZER.md):
  proves that the standard marked target orbit is a point in every degree
  and after every number of identity stabilizations.  The unstabilized step
  uses the coordinate-polynomial intruder theorem at `P^2*B^5*C`; Kuroda's
  stable-invariant theorem applied to conjugated stable translations forces
  the arbitrary-stabilization descent.  Vertical automorphisms of the
  identity factors are not classified.  Exact recursive Newton-face pruning
  eliminates quotient degrees eight through twelve, while an all-degree
  graded resolution puts every logarithmic generator in degree seven or
  eight.  Modulo the Koszul submodule, the exact quotient has dimension two,
  degree 296, and support contained in the singular locus.  That singular
  locus has exactly four
  minimal components: prime triple-root and two-double-root curves of
  degrees 17 and 19 and two lines at infinity.  The positive upper Newton hull is
  exactly the edge joining `P^12*C^4` and `P^2*B^5*C`.  The first tie in a
  `P`-zero Koszul ladder fails at target degree fifty, but its degree-55 rung
  cancels at leading order.  An exact depth-nine recursion excludes the
  normalized two-generator continuation of that rung.  Those remaining wall
  calculations concern the full vertical stable group, not the marked
  orbit.  A separate Newton-stratified certificate gives
  `chi(H=h)=246`, so the vanishing-`H^2` cylinder shortcut is unavailable
  but unnecessary.  Directly, the remaining full-group route is the third
  `P`-zero Koszul coefficient and the four classified
  singular-support charts.
- [Clean quadratic-gauge decorated receiver](QUADRATIC_GAUGE_DECORATED_RECEIVER.md):
  constructs the marked quotient mapping to `BS_N` and separates its
  `N-1`-dimensional geometric fibre from the `N-4`-dimensional image in
  unmarked stable-map moduli.  The weight-one coordinate `u_5/u_4` gives a
  global slice, so the clean quotient receiver is represented by an
  explicit scheme and has no residual coefficient-scaling inertia.  The
  receiver span makes precise that the fibre algebra descends while the
  ambient Fitting decoration does not.
- [Stable intruder descent and physical inertia](STABLE_INTRUDER_DESCENT_CRITERION.md):
  packages the arbitrary-stabilization argument into three exact gates:
  an intruder forces descent to the physical polynomial ring, faithful
  normalized boundary data forces physical target identity, and a trivial
  generic deck group forces physical source identity.  It proves trivial
  pointwise physical inertia modulo vertical stabilization gauge, not a
  classification of that vertical group or a global stable-map stack.
- [Keller/Tschirnhaus bridge card in ranks five through seven](KELLER_TSCHIRNHAUS_DESCENT_567.md):
  is the exact low-rank regression for the generic theorem, links the
  canonical all-rank projective and stable-boundary theorems,
  proposes Keller target dimension and Keller coordinate degree as research
  invariants, and pins one exact nonprojective Tschirnhaus experiment whose
  two compiled complete fibres have the same abstract split algebra.  It
  does not assert an ambient equivalence carrying one marked fibre to the
  other.
- [Rank-three collision-framed descent](RANK_THREE_COLLISION_DESCENT.md):
  identifies the cubic off-diagonal collision cover with the full `S_3`
  frame torsor, gives the exact projective and quadratic-Tschirnhaus
  transition cocycles, lifts them to the normalized factorization map after
  target localization, and isolates the remaining global polynomial boundary
  at one explicit denominator.  Within the canonical projective transport,
  the only denominator-free global subgroup is the known scaling torus.
- [Rank-four collision frames and cross ratio](RANK_FOUR_COLLISION_CROSS_RATIO.md):
  identifies the ordered three-root cover with the full `S_4` frame torsor
  and factors the fourth-root interpolation residual by the exact
  cross-ratio defect
  `q_2^2-q_1q_3+q_2q_3e_1+q_3^2e_2`.  It separates this projective
  hypersurface from the primitive-element boundary and writes its cleared
  equation in the universal quartic Keller coordinates; nonprojective Keller
  transport off that hypersurface remains open.
- [All-rank collision frames and projective descent](ALL_RANK_COLLISION_PROJECTIVE_DESCENT.md):
  identifies `Conf_(N-1)` with the full `S_N` frame torsor and proves that
  two primitive root coordinates are projectively related exactly when the
  four columns `1,r,u,r*u` have rank at most three.  The resulting
  determinantal locus is smooth of codimension `N-3` on the primitive open,
  has an explicit normalized coefficient matrix on the universal Keller
  chart, and specializes to the automatic cubic and quartic cross-ratio
  results.  It isolates, but does not construct, the genuinely
  nonprojective Keller transport required off this locus.
- [Rank-four nonprojective Keller descent](RANK_FOUR_NONPROJECTIVE_KELLER_LIFT.md):
  refines the one-point geometric quartic quotient by its exact
  `K*/(K*)^5` ground-field class, constructs an
  arithmetic-neutral primitive quadratic witness, and reduces it to two
  fibers of the fixed map `F_(-124416)`.  It verifies a finite-etale
  straight target line, proves that line has the wrong collision-frame
  sheet partition, and constructs the exact polynomial first-order and
  all-finite-order formal framed lifts.  The `-4` iterate drops the fiber
  cardinality from four to two, proving that the straight target translation
  has no polynomial lift.  The prime discriminant has ordinary degree
  thirteen, so every target self-equivalence through degree twelve is
  exactly in `mu_5`, whose endpoint orbit also fails.  An exact
  logarithmic-boundary reduction and Singular unit-ideal certificate
  exclude endpoint degrees thirteen through eighteen.  A target symmetry
  of degree at least nineteen with the prescribed sheet permutation
  remains open.
- [Universal atomic-map adversarial audit](UNIVERSAL_ATOMIC_MAP_ADVERSARIAL_AUDIT.md):
  attacks polynomiality, the promoted block Jacobian, generic degree,
  coefficient compilation, finite-etale completeness, monodromy, and stable
  atomicity separately; pins connected, split, and product witness targets;
  records every genuine failure boundary; and distinguishes the Lean-verified
  compiler/fiber core from the unformalized monodromy-to-atomicity chain.
- [Universal Keller-fiber multiplicity](UNIVERSAL_KELLER_FIBER_MULTIPLICITY.md):
  over every characteristic-zero field, every finite etale algebra of rank
  at least three is a complete fiber in infinitely many stable classes.
- [Whole-plane stable multiplicity](WHOLE_PLANE_STABLE_MULTIPLICITY.md):
  in every rank at least three, the infinite stable gauge families agree as
  maps on `P=1` and share the full squarefree inverse cover over its
  two-dimensional target plane.  Residue-class subfamilies agree over
  `P^d=1`; the fixed quintic specialization shares one dense quantitative
  family of degree-optimal Hasse failures.
- [Two marked fibers recover the power-shift gauge](TWO_MARKED_FIBER_GAUGE_RECONSTRUCTION.md):
  for the degree-at-least-four power-shift family, `P=1` recovers the
  normalized seed and a second marked fiber at any non-torsion `P=c`
  recovers the exponent.  The universal choice `c=2` works in
  characteristic zero; torsion planes give the exact congruence
  counterexamples, while one transverse marked line reads the exponent as
  a pole order at `P=0`.
- [Finite marked-plane nonreconstruction](FINITE_MARKED_PLANE_NONRECONSTRUCTION.md):
  outside the monomial power-shift locus, polynomial multipliers interpolate
  the value one on any finite set of target planes.  The resulting maps
  agree on every sampled marked inverse cover but have growing numbers of
  degree-drop boundary components, proving that no universal finite number
  of marked fibers recovers an unrestricted stable class in degrees at
  least four.
- [Polynomial-gauge decorated Torelli](POLYNOMIAL_GAUGE_DECORATED_TORELLI.md):
  the intrinsic boundary recovers the base character and the unmarked
  ramified-stratum Fitting divisor recovers the root coordinate, seed, and
  polynomial multiplier up to exactly the ordinary two-scaling action.
  More generally, the full finite normalization morphism plus its
  reconstruction boundary is a complete stable invariant by restriction to
  its distinguished affine open.
- [Universal cubic gauge multiplicity](UNIVERSAL_CUBIC_GAUGE_MULTIPLICITY.md):
  fiber-invisible cubic lifts keep the selected inverse quotient fixed while
  their canonical boundary-component count grows with the lift exponent.
- [Universal power-shifted gauge multiplicity](UNIVERSAL_POWER_SHIFTED_GAUGE_MULTIPLICITY.md):
  one common extra `P`-power on every decoration of degree at least four
  fixes the selected finite-etale fiber and moves the stable Fitting Newton
  area `2N-3+(N-2)m`.
- [Universal quartic gauge multiplicity](UNIVERSAL_QUARTIC_GAUGE_MULTIPLICITY.md):
  power-shifted quartic decorations keep the inverse quotient fixed at
  `P=1` and give pairwise distinct stable Fitting-support indices `2m+5`.
- [Universal multiplicity adversarial audit](UNIVERSAL_MULTIPLICITY_ADVERSARIAL_AUDIT.md):
  stress-tests the generator, clean-torus, invariant, full-fiber,
  power-shifted quartic, Hasse--Minkowski, and selected-root steps and records
  the exact imported stable-normalization dependencies.
- [Universal multiplicity witness cards](UNIVERSAL_MULTIPLICITY_WITNESS_CARDS.md):
  three exact pairwise stably inequivalent presentations of one connected
  field in each of degrees four, five, and six.
- [Low-rank multiplicity boundaries](LOW_RANK_MULTIPLICITY_BOUNDARIES.md):
  an anisotropic biquadratic quartic trace-chord form over
  `Q((a))((b))`, and the exact collapse of all three current cubic
  mechanisms to the foundational stable class.
- [Universal quartic fiber multiplicity](UNIVERSAL_QUARTIC_FIBER_MULTIPLICITY.md):
  every rank-four finite etale algebra over a number field is a complete
  fiber in infinitely many stable classes of determinant-one weighted maps;
  the proof combines the quartic trace-chord quadric, Hasse--Minkowski, and
  weighted selected-root Torelli.
- [Universal quintic fiber multiplicity](UNIVERSAL_QUINTIC_FIBER_MULTIPLICITY.md):
  over every characteristic-zero field, every rank-five finite etale algebra
  is a complete fiber in infinitely many stable quadratic-gauge classes;
  translation moves `a_5^5/(a_3 a_4^6)` nontrivially after choosing a
  generator with nonzero second trace moment.
- [Exact real-sheet spectrum](REAL_FIBER_SPECTRUM.md): every count
  `N,N-2,...,N mod 2` occurs on a nonempty complete regular real target
  chamber, with rational witnesses and an explicit fold-adjacency chain.
- [All-degree scalar vacua](ALL_DEGREE_SCALAR_VACUA.md): Zhu's flat
  unit-volume scalar pullback applied to the complete rational fibers gives,
  for every `N>=3`, an explicit three-scalar theory with exactly `N`
  rational isolated vacua, together with metric incompleteness and the
  varying-multiplicity quantum obstruction.
- [Adelic complete-fiber engineering](ADELIC_FIBER_ENGINEERING.md): weak
  approximation combines any allowed real signature with finitely many
  squarefree local splitting types; one local `N`-cycle gives a complete
  degree-`N` fiber field.
- [Local-to-global Keller fibers](LOCAL_GLOBAL_KELLER_FIBERS.md): arbitrary
  finite étale `Q_p`-algebras, including ramified ones, are combined with a
  real signature and further Frobenius types by coefficient CRT and then
  compiled into a complete determinant-one Keller fiber.  The explicit
  quintic certificate prescribes ramified local factors at both `2` and `3`;
  a universal `2*v_p(Disc)+1` coefficient radius makes the local
  certification automatic.  Reusable arithmetic and Keller compilers provide
  the exact end-to-end construction.  The Keller compiler's optional
  `stable_parameter` emits those infinitely many maps together with the
  boundary-count or Newton-area separation certificate; the portable JSON
  layer and its independent Python and PARI/GP replayers cover both cubic
  and power-shifted outputs.  In every rank at least three, the resulting
  fixed fiber algebra occurs in infinitely many stable map classes.
- [Locally prescribed common fibers](LOCALLY_PRESCRIBED_COMMON_FIBERS.md):
  two fixed stably inequivalent determinant-one maps share infinitely many
  connected fibers with any family-compatible finite collection of local
  algebras and one real signature.  An explicit sextic common field has
  ramified completions at `2` and `3`, signature `(2,2)`, and is inert at `5`.
- [Marked dyadic stable separation](MARKED_Q2_STABLE_SEPARATION.md):
  the connected quintic `T^5+T^3-2T^2+T+1` has marked local action
  `sigma=(1234)(5), tau=x_0=x_1=1` over `Q_2` and is the identical complete
  fiber of two determinant-one maps with stable unit ranks `1` and `2`.
- [Hasse-principle failure for a Keller fiber](HASSE_PRINCIPLE_KELLER_FIBER.md):
  an explicit degree-eight complete regular fiber has points over `R` and
  every `Q_p`, but no rational point.
- [Normal coverings and Hasse-failing fibers](NORMAL_COVERING_HASSE_FIBERS.md):
  arithmetic component stabilizers form a faithful normal covering of the
  splitting-field group, so the component count is at least `gamma(G)` and
  the fiber rank is their subgroup-index sum.  Exact `S_3` and `C_2^2`
  certificates accompany a pinned necessary-candidate transcription of
  Banks' degree-`5` through degree-`10` table and a determinant-one sextic
  Keller realization.
- [Dense multiplicative Hasse family](MULTIPLICATIVE_HASSE_KELLER_FIBERS.md):
  one fixed degree-five map has a Hasse-failing fiber for every noncube
  `a=1 mod 9` supported on primes `1 mod 3`.  Character-filtered
  Selberg--Delange gives order `B/sqrt(log B)` targets of height at most
  `B`; a dependency-free enumeration records all parameters through
  `a=10^6`.
- [Exact geometric-degree spectrum](GEOMETRIC_DEGREE_SPECTRUM.md): the
  spectrum `3,4,5,...` for noninvertible Keller maps of complex affine
  three-space and its stable left--right degree separation.
- [Atomic spectrum and non-generation](ETALE_MONOID_ATOMIC_SPECTRUM.md):
  the exact atomic degree spectrum `3,4,5,...`, degreewise stable atomic
  lower bounds, the forced-atomic/decomposable degree dichotomy, all
  multiplicative degree words in quadratic-gauge atoms, and an explicit
  quartic atom outside the monoid generated by automorphisms and every
  positive-dimensional quasi-torus Keller class.
- [Stable normalization functoriality](STABLE_NORMALIZATION_FUNCTORIALITY.md):
  the construction-independent theorem for normalization, boundary valuations,
  intersections, nilpotents, relative differentials, Fitting ideals, and
  conductor decorations after adjoining identity variables;
- [Quadratic-versus-weighted stable separation](QUADRATIC_WEIGHTED_STABLE_SEPARATION.md):
  complete quadratic boundary exhaustion, intrinsic ordering of the two
  target boundary images, and the stable unit-rank obstruction
  `G_m^2` versus `A^1 x G_m`;
- [Quadratic-gauge stable moduli](QUADRATIC_GAUGE_STABLE_MODULI.md):
  exact stable-orbit classification on the coefficient-torus locus; an
  overlooked independent `P`-scaling makes the quotient dimension `N-4`,
  rather than `N-3`.  The universal discriminant intruder `P^2*B^N*C`
  additionally makes the physical marked stabilizer trivial for every
  clean fixed map with `N>=5`, in arbitrary degree and stabilization
  dimension;
- [Quadratic-gauge/cancellation stable intersection](QUADRATIC_CANCELLATION_STABLE_INTERSECTION.md):
  the two families have exactly one common stable class, the foundational
  cubic; the all-degree separation is certified independently by the
  ramified-stratum Fitting support and by boundary-contact nilpotency;
- [Incidence suspensions through degree four](INCIDENCE_SUSPENSION_DEGREE_FOUR_CLASSIFICATION.md):
  exact marked-line criterion and complete root-preserving,
  `P`-fibration-preserving affine rechart search; the cubic and quartic
  reciprocal charts have unavoidable source poles, leaving only `X=S^2`;
- [Constant-kernel quotient](CONSTANT_KERNEL_QUOTIENT.md): the general
  triangular-extension and fiber-scheme theorem, its GZ-type context, the
  verified 24-to-22 quotient, and the mandatory essential-dimension search
  protocol.
- [Support-saturation principle](SUPPORT_SATURATION_PRINCIPLE.md): the
  equivalent local-cohomology, associated-prime, grade, regular-element,
  and presentation-saturation criteria that extend a defect across its
  possible support.

The [external collision-ideals audit](COLLISION_IDEALS_EXTERNAL_AUDIT.md)
checks the dimension-independent collision-ring and cubic `S_3` interfaces
against a pinned Lean development.  Its credit ledger separately records the
manuscript's attribution to Chloe van der Vlugt, publication of the Lean
project by the GitHub account `what-social-construct`, the manuscript's
AI-assistance disclosure, literature inputs, and the independent provenance
of the explicit three-dimensional counterexample.  It records the exact
conditional boundary of the planar hidden-inertia argument and does not
change the plane Jacobian-conjecture status.

External formal certificates for the foundational map are the pinned
[Lean development](LEAN_FOUNDATIONAL_MAP.md) and the independently authored,
refereed [Archive of Formal Proofs Isabelle/HOL
entry](https://isa-afp.org/entries/Jacobian_Counterexample.html).
Expanded audit and normalization narratives are preserved under
[archive/core-support](../archive/core-support/README.md), while their commands
remain in the public [reproduction guide](../REPRODUCE.md).

Start with [FOUNDATIONAL_GEOMETRY.md](FOUNDATIONAL_GEOMETRY.md). The former
omnibus paper has been retired; these focused theorem notes remain the
canonical sources.

These documents are the primary core references.

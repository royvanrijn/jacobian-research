# Verified core

This directory contains the stable proof chain:

- [Foundational Keller map](FOUNDATIONAL_GEOMETRY.md): exact determinant and collision;
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
  quotient reconstruction, optimal quintic Hasse fiber, updated arithmetic
  chain, and staged Lean certificate.
- [Universal relative Keller map](UNIVERSAL_RELATIVE_KELLER_MAP.md): packages
  all supplied presentations into one relative Jacobian-one map, compresses
  the map base sharply from `N+1` to `N-3` parameters by moving three inverse
  coefficients into the target, identifies the universal root fiber over the
  `N`-dimensional incidence open, and separates the `BS_N` atlas from the
  obstructed Tschirnhaus descent and essential-dimension questions.
- [Universal Keller-fiber multiplicity](UNIVERSAL_KELLER_FIBER_MULTIPLICITY.md):
  over every characteristic-zero field, every finite etale algebra of rank
  at least three is a complete fiber in infinitely many stable classes.
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
  rather than `N-3`;
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

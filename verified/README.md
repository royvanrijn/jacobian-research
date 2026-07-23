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
- [Exact real-sheet spectrum](REAL_FIBER_SPECTRUM.md): every count
  `N,N-2,...,N mod 2` occurs on a nonempty complete regular real target
  chamber, with rational witnesses and an explicit fold-adjacency chain.
- [Adelic complete-fiber engineering](ADELIC_FIBER_ENGINEERING.md): weak
  approximation combines any allowed real signature with finitely many
  squarefree local splitting types; one local `N`-cycle gives a complete
  degree-`N` fiber field.
- [Hasse-principle failure for a Keller fiber](HASSE_PRINCIPLE_KELLER_FIBER.md):
  an explicit degree-eight complete regular fiber has points over `R` and
  every `Q_p`, but no rational point.
- [Exact geometric-degree spectrum](GEOMETRIC_DEGREE_SPECTRUM.md): the
  spectrum `3,4,5,...` for noninvertible Keller maps of complex affine
  three-space and its stable left--right degree separation.
- [Stable normalization functoriality](STABLE_NORMALIZATION_FUNCTORIALITY.md):
  the construction-independent theorem for normalization, boundary valuations,
  intersections, nilpotents, relative differentials, Fitting ideals, and
  conductor decorations after adjoining identity variables;
- [Constant-kernel quotient](CONSTANT_KERNEL_QUOTIENT.md): the general
  triangular-extension and fiber-scheme theorem, its GZ-type context, the
  verified 24-to-22 quotient, and the mandatory essential-dimension search
  protocol.

The external Lean certificate for the foundational map remains in
[LEAN_FOUNDATIONAL_MAP.md](LEAN_FOUNDATIONAL_MAP.md).
Expanded audit and normalization narratives are preserved under
[archive/core-support](../archive/core-support/README.md), while their commands
remain in the public [reproduction guide](../REPRODUCE.md).

Start with [FOUNDATIONAL_GEOMETRY.md](FOUNDATIONAL_GEOMETRY.md). The paper
version is being assembled under
[papers/core-counterexample](../papers/core-counterexample/main.tex).

These documents are the primary core references.

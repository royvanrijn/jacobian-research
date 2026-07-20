# Implementation status against the requested program

## Complete and executed

- Explicit Bass--Connell--Wright/Yagzhev reduction to a 95-dimensional
  cubic-homogeneous map with 148 cubic terms, determinant one, nilpotent
  nonlinear Jacobian, and a stored three-point rational collision.
- Explicit Gorni--Zampieri/Druzkowski pairing to a 510-dimensional
  cubic-linear map `X-(AX)^{*3}`, with `rank(A)=95`, exact pairing matrices,
  determinant-one certificate, and a stored three-point rational collision.
- msolve 0.10.1 installed; F4 over `Q` and large prime fields validated.
- F4SAT component removal validated on a planted reducible system.
- Exact triangular chart `Phi_p`, inverse, Jacobian, and Laurent pullback.
- Exact sparse pole-cancellation matrices and nullspaces.
- Requested rectangular support box processed for canonical representative
  strata through `deg p=4`; sparse-first and three independent bound increases
  also processed.
- Rational/modular rank comparison with bad-prime reporting.
- Stored earliest forbidden pole (`x^-1` in every current stratum).
- Collision and derivative normalization built into nonlinear searches.
- Five named 10–30-coefficient support families solved with exact F4; each has
  a reproducible characteristic-zero input and reduced basis `[1]`.
- Saturated modular runs over three large primes for those families.
- Candidate metrics for exact Jacobian, generic resultant degree, degree-one,
  prime/small degree flags, zero-fiber collisions, obvious line restrictions,
  Newton edges, and a first Puiseux-leading screen.
- Canonical invariant signatures and short-word chart deduplication; depth two
  currently gives 20 sampled invariant classes.
- Three two-divisor unimodular toric charts, with 10- and 15-term cones, solved
  over `Q` and three primes; all six exact systems are empty.
- A greedily reduced 19-equation degree-three contradiction core, exactly
  certified over `Q`.
- The translated non-toric two-divisor chart, with exact simultaneous
  cancellation along `x=0` and `xy+1=0`, and five executed 13–49-term families.
- An exact identification of every symmetric Laurent box in that chart: the
  `[-n,n]^2` kernel is precisely the ordinary degree-`<=n` polynomial space.
- Validation controls: original 3D reconstruction, two weighted-lift seed
  deformations, planted cancellations, polynomial automorphisms,
  characteristic-p fake rejection, F4, and F4SAT.
- Exact ordinary-to-Laurent translation for Newton supports in the translated
  chart, with independent cancellation kernels for the two coordinates.
- Four thin unequal-coordinate Stage D systems: `(3,2)` and `(2,3)` weighted
  bands and both orientations of the primitive `(2,7)` square/cube prototype.
  Each is empty over `Q` and three large prime fields after collision
  normalization, while retaining the identity Keller control.
- The published `(9,27)` terminal obstruction is now an exact regression test:
  eight coefficients are eliminated from the nine displayed Section 5
  equations, reproducing equation (5.9), followed by the checked valuation and
  degree split in (5.11).

## Implemented but deliberately conditional

- HomotopyContinuation.jl is installed and validated on the nonempty,
  zero-dimensional 3D reconstruction residue. It finds one isolated numerical
  solution matching the exact coefficients. It remains deferred for empty or
  positive-dimensional search ideals.
- Hensel lifting/rational reconstruction is implemented and tested, but no
  modular-only candidate has survived to require it.
- Candidate filters run on explicit survivors. Current collision-normalized
  families have none, so there is no field degree or infinity geometry to
  classify beyond controls.

## Still incomplete

- The `p` sweep contains canonical representatives plus exhaustive lower-
  coefficient enumeration over `F_3` and `F_5` and random large-prime samples.
  Rank is constant degree-by-degree in all of them. A characteristic-zero proof
  that rank depends only on `deg p` is not yet stored.
- The full 119+119 coefficient nonlinear ideal for the requested rectangle has
  not been launched. Pole reduction leaves between 39 and 63 parameters per
  coordinate, still too large without sparse/gauge decomposition.
- Canonical signatures are rigorous invariants but not a complete decision
  procedure for equivalence under every affine and triangular automorphism.
- Full Newton–Puiseux expansions and dicritical-component resolution are not
  implemented; the current filter is a one-leading-term screen.
- The explicit Section 5 terminal system for `(9,27)` is transcribed, but the
  earlier reductions to it still depend on cited propositions rather than an
  independent implementation. Both broad Proposition 4.3 `(72,108)` systems
  are now encoded; only the smaller one received a bounded modular run, which
  timed out. Leading-form constraints remain only experimentally encoded.
- The denser total-degree 6/9 edge control has 52 kernel parameters and exceeded
  a 30-second modular F4 budget. It is a scaling datum, not an emptiness result.
- Both Proposition 4.3 reduced `(72,108)` polygons are encoded with all lattice
  coefficients, bracket `x^2`, and nonzero-vertex saturation. Their sizes are
  186/302 and 72/92 variables/equations. The smaller modular F4SAT run timed out
  at 30 seconds; no emptiness claim follows.
- An experimental common square/cube edge substitution reduces those variable
  counts to 163 and 49. The 49-variable modular run also timed out at 30 seconds.
  Its exhaustiveness remains conditional on auditing the cited leading-form
  results.
- Collision normalization is intentionally absent after the Proposition 4.3
  rational transformations; it has not yet been transported into those charts.

These boundaries matter: the stored `[1]` certificates are real finite-family
theorems, while the uncompleted bullets are research implementation tasks, not
negative mathematical results.

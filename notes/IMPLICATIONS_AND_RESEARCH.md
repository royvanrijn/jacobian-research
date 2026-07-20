# Implications and research directions

## Immediate implications

1. The classical Jacobian conjecture is false in dimensions at least three.
2. Any claimed proof for all dimensions contains a false step; the example is a
   compact regression test for locating it.
3. Nonproperness at infinity is the essential geometry. Étaleness prevents
   finite sheets from merging but does not prevent sheets escaping to infinity.
4. Generic fiber degree, monodromy, the nonproperness set, and valuations at
   infinity become primary invariants for classifying Keller maps.
5. The non-normal cubic with `S_3` Galois closure demonstrates the first field
   pattern not covered by birational or Galois-case invertibility results.

## High-priority research questions

- **Independent certification.** Produce short hand-checkable determinant
  proofs, reductions modulo several primes, and verification in two unrelated
  CAS implementations.
- **Provenance.** Recover the actual search prompt, code, random seeds, and
  intermediate candidates used to discover the example. Same-day expositions
  are not substitutes for a citable paper.
- **Minimality.** Minimize maximum coordinate degree, number of monomials,
  coefficient height, and weighted degree under polynomial automorphisms.
- **Normal forms.** Apply affine and triangular changes, tame automorphisms,
  Gröbner/SAGBI reduction, and weighted projective compactification. Determine
  whether degrees `(7,6,4)` or support can be lowered.
- **Infinity.** The normalized inverse graph and all dicritical prime divisors
  are classified in `DICRITICAL_COMPACTIFICATION.md`. Compute a minimal
  explicit blow-up sequence, discrepancies, and the full boundary intersection
  graph, then compare it with the normalized-graph model.
- **Topology.** Root meridians now explicitly give transpositions generating
  `S_3` (see `IMAGE_AND_NONPROPERNESS.md`).  Compute a presentation of the full
  fundamental group of the discriminant complement and its action, and compare
  it with the resolved boundary model.
- **Gradient dynamics.** Resolve the nonhyperbolic equilibrium line in the
  weighted infinity chart for the three-minimum polynomial; determine which
  stable sets at infinity form its basin boundaries. The first exact chart
  calculation is recorded in `GRADIENT_INFINITY.md`.
- **Families.** The universal monodromy, generic discriminant, boundary, and
  contact-partition theorems now replace extrapolation from the 627-seed scan.
  Continue the exceptional-strata calculation in degrees eight and higher,
  make bad-prime exclusions effective, and classify when two seeds give
  polynomially equivalent maps.
- **Known reductions.** The Bass–Connell–Wright/Yagzhev and
  Gorni–Zampieri/Druzkowski reductions have now been executed explicitly,
  yielding certified counterexamples in dimensions 95 and 510 respectively.
  The next task is to optimize the stable-equivalence choices and reduce these
  dimensions and support sizes without losing the transparent certificates.
## Simplification program

Use a multi-objective score

`(max degree, total monomials, max coefficient, weighted support volume)`.

Generate neighbors by affine changes, elementary shears, permutations, diagonal
scalings, and cancellation-aware rational refactorizations. Preserve exact
constant determinant and a certified collision at every step. Generic fiber
degree 3 is invariant under polynomial coordinate changes, so it is also a
useful guard against accidentally switching families.

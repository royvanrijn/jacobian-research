# Open problems for the cancellation programme

The cancellation family has separate construction, arithmetic, boundary, and
rigidity theorems.  Unrestricted parameter equivalence, arithmetic
classification, and universality beyond the current skeleton remain open.

## 1. Finish parameter equivalence

Compute the polynomial target automorphism group preserving the labelled
boundary pair and determine which elements lift to the finite normalization.
The visible quotient is the weighted scaling torus

\[
 (Q,R)\mapsto(uQ,u^{-m}R),
\]

and it fixes every parameter root.  The decisive loose end is the residual
congruence kernel acting trivially on `P=0`.  Proving that its liftable part is
trivial would finish distinct-root left--right equivalence.

## 2. Complete the arithmetic

Prove or disprove irreducibility of `M_{m,r}` for all `m,r`, and classify its
natural Galois group in all degrees.  The complete range `mr<=30`, the full
`m=1` column, and three uniform irreducibility criteria are recorded in
[ARITHMETIC.md](ARITHMETIC.md), but they do not yet imply an all-parameter
theorem.

Determine also the minimal field over which the full symbolic collision is
defined, rather than only the evident compositum of the parameter field and
the collision polynomial's splitting field.

## 3. Compute the stronger upstairs invariant

The target scheme intersection already separates the present examples, but
the next layer should compute the full upstairs--downstairs object:

- all boundary primes over each target stratum;
- their multiple intersections and specialization relations;
- residue-field covers;
- completed local maps, conductors, differents, and valuation filtrations;
- local monodromy around intersections.

This will test whether the formal ladder in
[BOUNDARY_GEOMETRY.md](BOUNDARY_GEOMETRY.md) remains computable after reduced
target intersections cease to distinguish examples.

## 4. Audit the remaining all-parameter inputs

Boundary normalization and exhaustion have a clean-room audit.  The thick
intersection formula now has independent resultant and completed-local-ring
proofs, together with bounded symbolic and standard-library exact-arithmetic
regressions.  The most valuable independent regressions still missing are:

1. the universal cancellation and reconstruction proof;
2. the all-parameter monodromy computation.

Positive invariance regressions should also conjugate a known map by
complicated source and target automorphisms, change primitive elements and
compactification coordinates, and stabilize it.  The reconstructed boundary
object should remain canonically isomorphic.

## 5. Go beyond the skeleton

[RIGIDITY.md](RIGIDITY.md) exhausts additional monomial weights, finitely many
normalized resolvent factors, and arbitrary target-dependent polynomial
derivatives while retaining one triangular reconstruction variable.  A
genuinely broader classification must change that skeleton—for example by
adding another source function, source variable, or independent inverse
variable—and then rebuild the cancellation and boundary analyses.

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

## 3. Prove the all-parameter contact resultant

The full reduced/scheme/formal invariant is now active, and the generic
weighted divisorial profile has executable `(e,f)`, different, completed-DVR,
residue-factor, inertia, and degree-sum extraction.  For cancellation maps,
the generic completed contact over `P=Q=0` is also executable: it retains
`k(R)[[P,Q]]`, its critical-normalization plane `k(R)[[Y,Q]]`, the pullback
`P=(Q-Y)Y^m`, and the separate branch contributions `m^2r` and `mr`.

The full prime-intersection diagram is now executable.  With
`K=K_{m,r}` and the critical coefficient `L=L_{m,r}`, the exact certificate
`Res(K,L) != 0` proves that every geometric `K`-branch has critical contact
length `m`; distinct `K`-branches meet in the common reduced central stratum.
This closes the diagram for the complete `m=1` ladder and for every certified
input, including the regression grid `1 <= m,r <= 5`.

The remaining uniform theorem is purely the all-parameter nonvanishing
problem

\[
 \operatorname{Res}(K_{m,r},L_{m,r})\ne0.
\]

A proof would remove the last per-parameter certificate from the diagram.
Completed finite-stratum conductors, restrictions of the different, valuation
filtrations, and local monodromy are refinements of the now-known incidence
graph rather than missing prime intersections.

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

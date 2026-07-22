# Status and open problems

> **Archived snapshot.**  This file preserves the component-era problem
> ledger.  The current arithmetic overview and frontier are maintained in
> [ARITHMETIC.md](../../cancellation/ARITHMETIC.md) and
> [the active open-problem list](../../cancellation/OPEN_PROBLEMS.md).  The
> fixed-row ramification advance is recorded separately in
> [FIXED_R_NEWTON_RAMIFICATION.md](../../cancellation/FIXED_R_NEWTON_RAMIFICATION.md).

## Current open frontier

The active questions are now sharply separated from completed milestones:

1. classify unrestricted polynomial left--right equivalence among maps
   attached to distinct cancellation parameter roots; target-fixed right-equivalence
   is now completely classified and separates distinct normalized roots;
   compute the enhanced upstairs--downstairs boundary incidence as a candidate
   obstruction to the remaining target automorphisms;
2. prove or disprove irreducibility of `M_(m,r)` outside the known uniform
   criteria, the fixed-`r` density-one theorem, the exact range `mr<=30`, and
   the additional certified cases `(49,1)` and `(97,1)`;
3. classify parameter Galois groups in all degrees and determine minimal
   collision fields; and
4. change the reconstruction skeleton itself, for example by introducing an
   additional source function, source variable, or inverse variable.

The generic-weighted-seed comparison, all allowed tail deformations,
separability, the complete `m=1` irreducibility column, fixed-`r` density-one
irreducibility, and the generalized mechanism inside the current ansatz are
no longer open.

## Archived transfer programme

The former transfer route is no longer part of the active problem list.  Its
finite examples, conditional all-`k` arguments, refuted filtrations, and
scoped cubic obstruction are preserved in
[the archive](../../archive/transfer-program/).

## Settled: `cancellation type (m,1)` versus generic weighted seeds

Fix an algebraically closed characteristic-zero field and `m>1`.  Let
`F_m:A^3 -> A^3` be a polynomial cancellation map with `r=1`, including
an arbitrary allowed tail `h_q(A)+A^2g(A)`.  Let `G_H` be a generic admissible
weighted-seed map whose inverse degree is `n=m+2`.

There do not exist polynomial automorphisms `alpha,beta` of `A^3` such that

\[
 G_H=\beta\circ F_m\circ\alpha.
\]

for `m>1`; the same is true after adjoining the same number of identity
coordinates.

### Boundary-intersection obstruction

The intrinsic data computed in
[RESOLVENT_RAMIFICATION_SIGNATURE.md](RESOLVENT_RAMIFICATION_SIGNATURE.md)
agree:

- generic degree `m+2` and monodromy `S_(m+2)`;
- generic critical partition `1^(m+1)`;
- one discriminant boundary prime with ramification index and sheet loss two;
- `m-1` geometric unramified boundary primes over a second target divisor,
  with total loss `m-1`.

The numerical data alone do not separate the maps.  However, the reduced
intersection of the two intrinsically distinguished target components is
`A^1 disjoint-union G_m` for the cancellation construction and `A^1` for a generic weighted seed.  The
first is disconnected and the second is connected, and polynomial
stabilization preserves that difference.  The complete proof is
[BOUNDARY_INTERSECTION_OBSTRUCTION.md](BOUNDARY_INTERSECTION_OBSTRUCTION.md).

For `m=1`, the explicit linear transformation in
[MASTER_CANCELLATION_CONSTRUCTION.md](MASTER_CANCELLATION_CONSTRUCTION.md)
identifies the cubic cancellation member with the foundational Keller map;
there is no second boundary component to compare.

## Settled: cancellation construction tail deformations

For fixed `(m,r,q)`, the construction permits every tail
`h_q(A)+A^(r+1)g(A)`.  The inverse resolvent, monodromy, and target
boundary-intersection data do not see `g` because the tails are already
source-equivalent.  Replacing `h_q` by `h_q+A^(r+1)g` is exactly
precomposition by the polynomial automorphism

\[
 (x,y,z)\longmapsto(x,y,z+y^{m+1}g(1+xy^m)).
\]

Thus tails do not define distinct polynomial left--right classes.

## Settled inside the ansatz: generalized cancellation

Use the notation of
[GENERALIZED_CANCELLATION_MECHANISM.md](GENERALIZED_CANCELLATION_MECHANISM.md).
For `e=1`, every polynomial leading solution is controlled by

\[
 I_n(q)=\int_0^1u\{1-q(1-u)\}^n\,du.
\]

The generalized cancellation note proves more generally that, for every
fixed `e>=1`, the spectral polynomials

\[
 J_{N,e}(q)=\int_0^1u^e\{1-q(1-u)\}^Ndu.
\]

are pairwise coprime.  Consequently every polynomial leading solution has
`f=c(y-alpha)^m` and `g_0=q(y-alpha)f`, and Hensel uniqueness makes its full
jet a translated/scaled cancellation jet.  The arbitrary `A^(e+1)` tail is removed by
a polynomial source automorphism.  Thus there is no nonmonomial branch, and
no new polynomial-equivalence class, anywhere in the ansatz of that note.

The universal normalized monomial-triangular three-weight class is now
settled in
[THREE_WEIGHT_CANCELLATION_CLASSIFICATION.md](THREE_WEIGHT_CANCELLATION_CLASSIFICATION.md):
all `a,b>=1`, `c>=0`, arbitrary `f`, arbitrary triangular `g`, and every
one-factor polynomial derivative are quantified.  Localized Keller forces
`c=a-1` and the old derivative power before polynomial cancellation is even
imposed; polynomiality then gives only cancellation construction.  Its assumption audit identifies
additional source functions, inverse variables, and nonmonomial `P,Q`
combinations as the genuine ways to leave the classified skeleton.

The alternative minimal relaxation is also settled in
[TWO_FACTOR_RESOLVENT_CLASSIFICATION.md](TWO_FACTOR_RESOLVENT_CLASSIFICATION.md):
a second normalized factor `1-tf_2(Q-Pt)` must coincide with the original
factor or be identically one.  The proof extends to any finite product of
normalized factors, all of whose nontrivial members coalesce.  Polynomiality
therefore again gives only cancellation construction.

The stronger
[TARGET_DEPENDENT_RESOLVENT_CLASSIFICATION.md](TARGET_DEPENDENT_RESOLVENT_CLASSIFICATION.md)
allows an arbitrary polynomial `H(T,P,Q)` as the derivative.  Algebraic
independence of `(s,P,Q)` forces the exact original power
`lambda(1-Tf(Q-PT))^e`.  Thus derivative generalization is completely closed
inside the coordinate skeleton.

**Remaining question.** Can a skeleton with an additional source function,
source variable, or inverse variable escape both Jacobian rigidity and
rising-factorial spectral coprimality? Equivalence among distinct cancellation construction
parameter branches remains a separate classification problem, but the
[target-fixed rigidity theorem](TARGET_FIXED_PARAMETER_RIGIDITY.md) shows
that any such equivalence must use a nonidentity target automorphism moving
the filled `P=0` boundary branch.

The
[target-boundary automorphism theorem](TARGET_BOUNDARY_AUTOMORPHISM_GROUP.md)
computes the restriction of that group to `P=0`: it is the weighted scaling
torus `(P,Q,R)->(lambda^(m+1)P,lambda Q,lambda^(-m)R)`, up to a residual
kernel congruent to the identity on the boundary plane.  The torus fixes
`w=PT/Q` and therefore fixes every parameter root.  The remaining equivalence
problem is to kill, or compute the finite root permutation induced by, the
cover-lifting part of that residual kernel.

## Parameter-polynomial arithmetic

This archived section is superseded by the canonical
[arithmetic overview](../../cancellation/ARITHMETIC.md) and
[active arithmetic frontier](../../cancellation/OPEN_PROBLEMS.md).  Exact
component proofs and finite certificates remain available in
[PARAMETER_IRREDUCIBILITY.md](PARAMETER_IRREDUCIBILITY.md),
[PARAMETER_DISCRIMINANT.md](PARAMETER_DISCRIMINANT.md), and
[PARAMETER_GALOIS_GROUPS.md](PARAMETER_GALOIS_GROUPS.md).

# Status and open problems

## Current open frontier

The active questions are now sharply separated from completed milestones:

1. classify unrestricted polynomial left--right equivalence among maps
   attached to distinct C24 parameter roots; target-fixed right-equivalence
   is now completely classified and separates distinct normalized roots;
2. prove or disprove irreducibility of `M_(m,r)` outside the known uniform
   criteria and the exact range `mr<=30`;
3. classify parameter Galois groups in all degrees and determine minimal
   collision fields; and
4. classify a second independent resolvent factor, now that the first
   one-additional-weight relaxation has been proved to collapse to C24.

The generic-weighted-seed comparison, all allowed tail deformations,
separability, the complete `m=1` irreducibility column, and the generalized
mechanism inside the current ansatz are no longer open.

## Archived all-`k` transfer question

The former C22 programme is no longer active.  Both proposed conceptual
models fail at `k=2`: the ribbon norm kills `X^2`, and the quadratic-remainder
tangent cone omits `X^3`.  The latter proposed special fiber is
positive-dimensional for every `k>=2`; see
[QUADRATIC_REMAINDER_ALGEBRA.md](QUADRATIC_REMAINDER_ALGEBRA.md).

This does not logically disprove finite flatness of the actual transfer ring,
but that narrower question is frozen rather than used as a dependency.  The
superseded programme and its exact historical regressions are preserved in
[the archive](../../archive/transfer-all-k/).

## Settled: `C24_(m,1)` versus generic weighted seeds

Fix an algebraically closed characteristic-zero field and `m>1`.  Let
`F_m:A^3 -> A^3` be a polynomial C24 cancellation map with `r=1`, including
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
`A^1 disjoint-union G_m` for C24 and `A^1` for a generic weighted seed.  The
first is disconnected and the second is connected, and polynomial
stabilization preserves that difference.  The complete proof is
[BOUNDARY_INTERSECTION_OBSTRUCTION.md](BOUNDARY_INTERSECTION_OBSTRUCTION.md).

For `m=1`, the explicit linear transformation in
[MASTER_CANCELLATION_CONSTRUCTION.md](MASTER_CANCELLATION_CONSTRUCTION.md)
identifies C24 with C01; there is no second boundary component to compare.

## Settled: C24 tail deformations

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
jet a translated/scaled C24 jet.  The arbitrary `A^(e+1)` tail is removed by
a polynomial source automorphism.  Thus there is no nonmonomial branch, and
no new polynomial-equivalence class, anywhere in the ansatz of that note.

The first minimal coordinate relaxation is now settled in
[THREE_WEIGHT_CANCELLATION_CLASSIFICATION.md](THREE_WEIGHT_CANCELLATION_CLASSIFICATION.md):
allowing `Q=y+xA^cB` with an independent third weight forces `c=a-1` and the
old derivative power before polynomial cancellation is even imposed.

**Remaining question.** Can a second independent resolvent factor still
admit a finite cancellation operator while escaping the rising-factorial
spectral coprimality argument? Equivalence among distinct C24 parameter
branches remains a separate classification problem, but the
[target-fixed rigidity theorem](TARGET_FIXED_PARAMETER_RIGIDITY.md) shows
that any such equivalence must use a nonidentity target automorphism moving
the filled `P=0` boundary branch.

## Parameter-polynomial arithmetic

The [prime-power Eisenstein theorem](PARAMETER_IRREDUCIBILITY.md) proves
`M_(m,r)` irreducible whenever

\[
 mr+r+1=p^k,\qquad v_p(mr)=k-1.
\]

This includes every pair for which `mr+r+1` is prime and shows that its maps
form one arithmetic Galois-conjugacy class.  It does not settle
irreducibility for arbitrary `(m,r)`.  Classical truncated-binomial theory
does settle the complete column `m=1`.  A second theorem in the same note
proves irreducibility when `mr+1` and `mr+r+2` are prime and `mr+r+2` is
primitive modulo `mr+1`; reduction then gives `Phi_(mr+1)`.  A third theorem
uses the unit-disk transform to cover `binom(mr+r,r)` prime, including every
`r=1` case with `m+1` prime.  The remaining questions are the other
truncated-binomial cases beyond the exact modularly certified range
`mr<=30`, higher Galois groups, and unrestricted polynomial left--right
equivalence between conjugate roots over an algebraic closure. Fixed-target
right-equivalence is already ruled out for distinct roots. The
[closed discriminant formula](PARAMETER_DISCRIMINANT.md) already settles
separability and alternating containment uniformly; square discriminants
occur in an explicit infinite family for every fixed `r`.  The exact
[degree-thirty table](PARAMETER_GALOIS_GROUPS.md) shows that full symmetric
groups coexist with `D_4`, the exceptional degree-six `S_5` action, and
alternating groups `A_12`, `A_16`, `A_17`, and `A_24`, so the higher-group
problem must allow genuine exceptional families.

# Transfer and cancellation: C17--C24

**Workflow state:** C24 is the sole active research branch.  The proposed
C22 transverse filtration is refuted at `k=2`; C18, C21, and all
transfer-dependent all-`k` claims have been archived.  C23's independent
affine-factorization-complement results remain in scope.

## What is the main theorem?

C17, C19, and C20 now supply only exact low-`k` local-algebra examples.
[C23](UNIVERSAL_FACTORIZATION_GEOMETRY.md) independently studies when
universal marked-factorization complements are affine spaces. C24
constructs the master cancellation family,
computes its generic inverse degree and monodromy, and begins its
classification under polynomial left--right equivalence.  The
[generalized cancellation note](GENERALIZED_CANCELLATION_MECHANISM.md)
extracts the determinant argument for arbitrary one-variable input `f(y)`
and reduces polynomiality to a finite functional equation.

## Why is it interesting?

The transfer blocks isolate the local algebra behind collisions of
factorization strata. The cancellation family supplies a broader source of
explicit noninjective Keller maps and asks whether visibly different inverse
resolvents can define equivalent polynomial maps.

## What does it depend on?

The retained low-`k` transfer examples use exact relative Groebner bases. C24
uses its explicit reconstruction formulas, finite normalization, DVR
ramification, and tame branch cycles. Its `(1,1)` member links back to
[C01](../../verified/FOUNDATIONAL_GEOMETRY.md).

## What is fully proved?

C17 remains an exact low-degree theorem.  The independent cubic C23
classification and `(2,3)` obstruction remain internally proved.  The
construction, degree, reconstruction, and monodromy in C24 remain internally
proved.  The `r>=2` ramification obstruction and the
[boundary-intersection obstruction](BOUNDARY_INTERSECTION_OBSTRUCTION.md)
for `r=1,m>1` separate every noncubic C24 member from generic weighted seeds.
The new
[quadratic-remainder audit](QUADRATIC_REMAINDER_ALGEBRA.md) proves that the
proposed C22 associated graded is false.

## What remains only computationally supported?

The explicit `Z_3` and `Z_4` presentations (C19--C20) and displayed small C24
coordinate maps remain exact bounded computations. They are examples, not
evidence for an active all-`k` transfer programme. Parameter arithmetic is
now stronger: three irreducibility criteria and the full `m=1` column are
uniform theorems, every pair `mr<=30` has an exact modular irreducibility
certificate, the discriminant is known uniformly, and Galois groups are
classified exactly through `mr<=30`.

## What is the likeliest failure point?

For transfer theory, the quadratic-remainder tangent cone is now known to
omit `X^3` at `k=2`, while the older Boolean norm omits `X^2`.  For active
C24 work, the risk is
mistaking a primitive-resolvent artifact for an invariant of the polynomial
map. Every proposed obstruction must be identified on the canonical finite
normalization.

## What review expertise is needed?

Finite covers and valuation theory, permutation monodromy, formal and
derived local algebra, invariant theory, and polynomial automorphism/stable
equivalence techniques.

## Completed boundary milestone and next target

The reduced intersection of the intrinsically distinguished target boundary
components is `A^1 disjoint-union G_m` for C24 with `r=1,m>1` and `A^1` for a
generic weighted seed.  This settles the formerly unresolved comparison,
including stable left--right equivalence.  The
[canonical boundary-incidence theorem](CANONICAL_BOUNDARY_INTERSECTION_INVARIANT.md)
now packages, for every characteristic-zero nonproper Keller map, the full
ramification-labelled diagram of target boundary divisors and all their
scheme-theoretic multiple intersections, together with its termwise
reduction.  Its hypothesis audit isolates quasi-finiteness, normality, finite
normalization, and boundary exhaustion.  Unique intrinsic markings are needed
only to extract a named intersection such as the one used here.

The
[upstairs--downstairs enhancement](UPSTAIRS_DOWNSTAIRS_BOUNDARY_INCIDENCE.md)
retains the boundary primes in the finite normalization, their actual finite
maps to target divisors, all upstairs intersections, their scheme images
inside downstairs intersections, and the relative different.  This bipartite
object is functorial and stable under affine cylinders; its C04/C24 vertex
layer is explicit, while the C24 upstairs intersection schemes are the next
geometric computation.

The [monodromy deck-rigidity theorem](MONODROMY_DECK_RIGIDITY.md) identifies
target-fixed automorphisms of a finite normalization with `N_G(H)/H`, or
equivalently with the centralizer of its transitive monodromy action.  A
self-normalizing sheet stabilizer makes the normalization rigid and reduces
comparison of two source models inside it to equality of their distinguished
affine opens.  This separates the general uniqueness argument in C24 from the
parameter-specific pole calculation.

The [boundary-exhaustion certificate](BOUNDARY_EXHAUSTION_CERTIFICATE.md)
separates support containment from local `sum ef` degree accounting.  Its
Theorem 5.1 states explicitly that the displayed C04/C24 divisors are all
components of the canonical normalization boundary and lists all their target
images.  This remains internally proved pending a clean-room compactification
audit.  The full C24 intersection is a disjoint
union of `G_m` and a thickened affine line with nilradical index `m(m+1)`;
the weighted intersection is reduced.  This gives a second stable obstruction
beyond connectedness of the reductions.

The [all-parameter thick-intersection theorem](C24_SCHEME_BOUNDARY_ALL_PARAMETERS.md)
extends the trace formula to every `C24_(m,r)` and the canonical stable
statement to every noncubic pair: the exact nilradical index is `mr(m+1)`,
and the reduced intersection is always
`A^1 disjoint-union G_m`.  Since the generic weighted intersection is reduced,
this gives one scheme-theoretic stable-inequivalence proof for every noncubic
C24 member, replacing the former split between `r=1` and `r>=2`.

The [generalized mechanism](GENERALIZED_CANCELLATION_MECHANISM.md) now
classifies the full two-weight ansatz.  For every `e>=1`, its spectral
polynomials are pairwise coprime; every polynomial leading term is monomial
after translation/scaling; and the unique full jet is C24, while the
`A^(e+1)` tail is removed by a polynomial source automorphism.  Hence a new
constructive family must relax the coordinate skeleton itself. Equivalence
among distinct C24 parameter branches is now narrowed by
[target-fixed rigidity](TARGET_FIXED_PARAMETER_RIGIDITY.md): two branches are
right-equivalent over the identity target exactly when their cancellation
jets agree modulo `A^(r+1)`, so distinct normalized roots are inequivalent
even stably in that category. Any remaining equivalence must use a
nonidentity target automorphism that moves the filled `P=0` boundary branch.
The first relaxation beyond that skeleton is also closed: the
[three-weight classification](THREE_WEIGHT_CANCELLATION_CLASSIFICATION.md)
allows the exponent in `Q=y+xA^cB` to vary independently, but the localized
Jacobian forces `c=a-1` and the old derivative power.  Polynomiality then
reduces to the two-weight theorem, so no additional branch occurs.
The alternative minimal relaxation is closed by the
[two-factor resolvent classification](TWO_FACTOR_RESOLVENT_CLASSIFICATION.md):
a second normalized factor is forced either to coincide with the original
factor or to be identically one.  The same Laurent-unit proof handles any
finite product of normalized factors: all nontrivial factors coalesce.  Thus
both proposed Option B relaxations produce only C24.
More strongly, the
[target-dependent resolvent theorem](TARGET_DEPENDENT_RESOLVENT_CLASSIFICATION.md)
allows an arbitrary polynomial derivative `H(T,P,Q)`.  Algebraic independence
forces it to be exactly the original power.  Option B is therefore complete
inside the coordinate skeleton; any new branch must change reconstruction
itself.
The
[prime-power Eisenstein theorem](PARAMETER_IRREDUCIBILITY.md) proves that all
roots form one arithmetic conjugacy class whenever `mr+r+1=p^k` and
`v_p(mr)=k-1`.  Its cyclotomic companion covers further composite values
when `mr+1` and `mr+r+2` are prime with maximal multiplicative order.  The
unit-disk leading-prime criterion also proves every `r=1` case with `m+1`
prime.  The general truncated-binomial irreducibility and geometric
equivalence questions remain open, but classical truncated-binomial theory
now proves the full `m=1` column irreducible.  Every pair with `mr<=30` is
separately proved irreducible by exact modular degree sieves.  The
[closed discriminant theorem](PARAMETER_DISCRIMINANT.md) proves uniform
separability, gives an all-parameter square criterion, and parametrizes an
infinite even-degree square-discriminant family for every fixed `r`.  The
[degree-thirty Galois table](PARAMETER_GALOIS_GROUPS.md) finds `D_4` at
`(2,2)`, the exceptional transitive `S_5` action at `(2,3)` and `(6,1)`, and
alternating groups at `(4,3)`, `(2,8)`, `(16,1)`, `(17,1)`, `(1,24)`, and
`(12,2)`, so a uniform symmetric-group conjecture is false.  No other branch
in this directory is active.

The archived all-`k` transfer programme is preserved under
[archive/transfer-all-k](../../archive/transfer-all-k/).  It is not an active
dependency chain.  Only the exact low-degree examples in
[TRANSFER_BLOCKS.md](TRANSFER_BLOCKS.md) remain in current navigation.

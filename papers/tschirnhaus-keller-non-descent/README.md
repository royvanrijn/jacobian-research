# Arithmetic descent and geometric non-descent for universal Keller fibres

This directory is the prospectus for a geometric companion paper.  It is not
yet a publication manuscript and is excluded from the paper build.  The
canonical proofs remain in `verified/`; this file records scope and assembly
order rather than duplicating them.

## Proposed headline

For every rank `N>=5`, the clean quadratic-gauge compiler realizes every
primitive finite-etale presentation as a complete Keller fibre.  A generic
Tschirnhaus change preserves the abstract fibre algebra but changes the
intrinsic stable Fitting-boundary fingerprint of the ambient map.  Hence the
fibre algebra descends while this ambient Keller decoration does not.

The precise receiver is the span

\[
 BS_N
 \longleftarrow
 {\mathscr D}_N^{\mathrm{quad},\circ}
 \longrightarrow
 {\mathcal M}_N^{\mathrm{quad}},
\]

not a map from unmarked stable Keller classes to `BS_N`.

## The theorem package already available

1. **Universal complete fibres.**  The relative quadratic-gauge compiler
   produces the prescribed finite-etale algebra as a literal complete
   fibre.
2. **Stable boundary quotient.**  The clean stable quotient is
   `G_m^(N-4)`, with complete fingerprint
   `(I_5,J_6,\ldots,J_N)`.
3. **Generic Tschirnhaus non-descent.**  The equal-boundary relation has
   codimension `N-4`; the projective relation has codimension `N-3`; a
   generic Tschirnhaus arrow is outside both.
4. **Uniform witness.**  The change `i -> i+i^2` is nonprojective and changes
   the displayed boundary in every rank.
5. **Rank-five transition geometry.**  For the presentation
   `(1,2,3,4,5)`, the ambient, projective, and intersection loci have local
   dimensions `4`, `3`, and `2`.
6. **Canonical marked transport.**  On the rank-five ambient hypersurface,
   the canonical coefficient-torus equivalence carries the selected fibre
   only on the one-dimensional root-scaling locus.
7. **Represented decorated receiver.**  The weight-one coordinate
   `lambda=u_5/u_4` gives a global slice for the residual scaling.  Hence
   the clean marked quotient is an explicit scheme of dimension `N-1`,
   its target-forgetting fibres have dimension three, and its image in
   unmarked stable-map moduli has dimension `N-4`.
8. **Stable intruder descent and fixed-map stabilizers.**  The abstract
   criterion separates stable base descent, boundary-faithful physical
   identity, and the source deck-group gate.  For every `N>=5`, the reduced
   discriminant has universal exposed intruder `P^2*B^N*C`.  Kuroda's
   stable-invariant theorem and the coordinate-polynomial theorem make the
   standard marked target orbit of every fixed clean map a point in every
   degree and after every number of identity stabilizations; its physical
   inertia is trivial modulo vertical stabilization gauge.

## Claims deliberately not made

- The equal-boundary and projective loci are not nested; there is no
  quotient of the first by `PGL_2` without additional structure.
- The number `N-4` is not a proved lower bound for global `ktdim` or `kdeg`.
- No finite-type moduli stack of all stable polynomial maps has been
  constructed.  The clean receiver itself is represented by a scheme, but
  this does not construct the surrounding stable-map groupoid.
- The full stable self-equivalence groups of the fixed clean maps are not
  classified: vertical automorphisms of the identity factors may remain.
  Their physical marked stabilizer orbits are nevertheless completely
  determined and cannot move.
- No Tschirnhaus-invariant alternative compiler is asserted.

## Proposed section order

1. finite-etale presentations and Tschirnhaus arrows;
2. the universal quadratic-gauge compiler and its complete fibres;
3. intrinsic Fitting boundary and the stable coefficient quotient;
4. generic arithmetic descent versus ambient non-descent;
5. the all-rank witness and the rank-five transition loci;
6. the marked decorated receiver over `BS_N`;
7. essential dimension, versality, and what the dimension count does not
   prove;
8. the all-rank fixed-map stabilizer theorem and alternative compilers.

The logical center should be the receiver span: it converts the slogan
“same fibre, different ambient map” into two honest morphisms with different
forgetful behaviour.

## Canonical sources

- [`GENERIC_TSCHIRNHAUS_NON_DESCENT.md`](../../verified/GENERIC_TSCHIRNHAUS_NON_DESCENT.md)
- [`RANK_FIVE_TSCHIRNHAUS_TRANSITION_LOCUS.md`](../../verified/RANK_FIVE_TSCHIRNHAUS_TRANSITION_LOCUS.md)
- [`RANK_FIVE_STABLE_TARGET_STABILIZER.md`](../../verified/RANK_FIVE_STABLE_TARGET_STABILIZER.md)
- [`STABLE_INTRUDER_DESCENT_CRITERION.md`](../../verified/STABLE_INTRUDER_DESCENT_CRITERION.md)
- [`QUADRATIC_GAUGE_DECORATED_RECEIVER.md`](../../verified/QUADRATIC_GAUGE_DECORATED_RECEIVER.md)
- [`QUADRATIC_GAUGE_STABLE_MODULI.md`](../../verified/QUADRATIC_GAUGE_STABLE_MODULI.md)
- [`UNIVERSAL_RELATIVE_KELLER_MAP.md`](../../verified/UNIVERSAL_RELATIVE_KELLER_MAP.md)
- [`ALL_RANK_COLLISION_PROJECTIVE_DESCENT.md`](../../verified/ALL_RANK_COLLISION_PROJECTIVE_DESCENT.md)

## All-rank fixed-map gate closed; next receiver gate

For every `N>=5`, let `F_a` be a fixed clean rank-`N` quadratic-gauge map
and let `y` be a standard marked target.  The uniform calculation determines

\[
 \operatorname{Stab}_{\mathrm{st}}^t(F_a)\cdot y.
\]

The orbit is a point in every degree and for arbitrary stabilization.  The
reduced discriminant contains `P^2*B^N*C`, uniquely exposed by the positive
weight `(1,N+1,N)`.  Kuroda's stable-invariant theorem applied to conjugated
stable translations forces both the target automorphism and its inverse to
preserve `k[P,B,C]`; the coordinate-polynomial intruder theorem then makes
that restriction the identity.  Full symmetric monodromy kills the physical
source deck group.

The receiver-side residual quotient is also finished on the clean chart:
`lambda=u_5/u_4` has weight one, and `lambda=1` is a global algebraic
slice.  Thus the clean receiver is an explicit scheme with no
coefficient-scaling inertia.  The stable intruder criterion further shows
that its pointwise physical inertia is trivial modulo vertical
stabilization gauge.  The remaining receiver problem is global: construct
and control that vertical quotient beyond the clean chart, or replace a
section by a Tschirnhaus-compatible correspondence or torsor.

The rank-five calculation below is now an explicit specialization carrying
additional information about the vertical/full-group equations.
Exact recursive Newton-face pruning eliminates the logarithmic spaces
through quotient degree twelve.  The logarithmic module is resolved in every
degree: it has two generators in quotient degree seven, thirteen in degree
eight, eighteen first relations, and six second relations.  Its quotient by
the Koszul submodule has dimension two, degree 296, and support contained in the
discriminant singular scheme.  Exact prime contractions, boundary
saturations, and root-partition exhaustion split that singular scheme into the
degree-17 triple-root curve, the degree-19 two-double-root curve, and two
lines at infinity.  The stronger Newton statement is now
settled: every positive weight exposes only `P^12*C^4`, only
`P^2*B^5*C`, or their common edge.  The edge wall is
`10*w_P-5*w_B+3*w_C=0`.  A `P`-zero Koszul ladder first reaches it at target
degree fifty, where the UFD power condition fails, and admits leading
cancellation at degree fifty-five.  A sparse exact recursion closes the
normalized two-generator continuation of this first cancellable rung at
depth nine.  The remaining `P`-zero Koszul coefficient and the four
classified singular-support charts now belong only to classification of the
full vertical stable group.  They are not a marked non-descent gate.  The
next paper-level step is instead to formulate the exact
correspondence/torsor-valued receiver or construct a genuinely
Tschirnhaus-invariant alternative compiler; the fixed-map kernel cannot
repair the displayed atlas.  An independent rank-five Newton calculation
gives `chi(H=h)=246`,
so the vanishing-`H^2` cylinder shortcut would not have closed this example;
Kuroda's ambient theorem does.

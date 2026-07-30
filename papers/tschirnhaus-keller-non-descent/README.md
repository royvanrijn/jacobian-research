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
7. **Decorated receiver.**  The marked quotient has dimension `N-1`, its
   target-forgetting fibres have dimension three, and its image in unmarked
   stable-map moduli has dimension `N-4`.

## Claims deliberately not made

- The equal-boundary and projective loci are not nested; there is no
  quotient of the first by `PGL_2` without additional structure.
- The number `N-4` is not a proved lower bound for global `ktdim` or `kdeg`.
- No finite-type moduli stack of all stable polynomial maps has been
  constructed.
- The full stable self-equivalence group of the fixed quintic map is not
  classified.  Consequently the current marked theorem concerns the
  canonical equivalence; arbitrary marked transport is reduced to a
  stabilizer-orbit problem.
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
8. the fixed-map stabilizer problem and alternative compilers.

The logical center should be the receiver span: it converts the slogan
“same fibre, different ambient map” into two honest morphisms with different
forgetful behaviour.

## Canonical sources

- [`GENERIC_TSCHIRNHAUS_NON_DESCENT.md`](../../verified/GENERIC_TSCHIRNHAUS_NON_DESCENT.md)
- [`RANK_FIVE_TSCHIRNHAUS_TRANSITION_LOCUS.md`](../../verified/RANK_FIVE_TSCHIRNHAUS_TRANSITION_LOCUS.md)
- [`QUADRATIC_GAUGE_DECORATED_RECEIVER.md`](../../verified/QUADRATIC_GAUGE_DECORATED_RECEIVER.md)
- [`QUADRATIC_GAUGE_STABLE_MODULI.md`](../../verified/QUADRATIC_GAUGE_STABLE_MODULI.md)
- [`UNIVERSAL_RELATIVE_KELLER_MAP.md`](../../verified/UNIVERSAL_RELATIVE_KELLER_MAP.md)
- [`ALL_RANK_COLLISION_PROJECTIVE_DESCENT.md`](../../verified/ALL_RANK_COLLISION_PROJECTIVE_DESCENT.md)

## Next proof gate

Let `F_R` be the compiled map for
`R(T)=prod_(i=1)^5(T-i)` and let `y_R` be its selected target.  The next
strong theorem is obtained by determining

\[
 \operatorname{Stab}_{\mathrm{st}}^t(F_R)\cdot y_R.
\]

If this orbit is a point, canonical rank-five marked non-descent upgrades to
full marked non-descent.  If it is positive-dimensional, its explicit
symmetries describe the residual ambiguity of the decorated receiver.  In
either case the computation is structurally informative and is the first
attack to run before strengthening the paper's marked claims.

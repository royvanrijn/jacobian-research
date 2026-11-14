# Remaining global input for omitted-divisor accounting

This is a completion audit, not a new theorem. The current mathematical
authority remains the PF2OMIT claims in `research/MATH_STATUS.json`.

The Fonollosa control has explicit equations, all omitted primes and
exceptional lengths. For the hypothetical residual immersed-r7 Keller
map, the current data determine necessary numerical and local rows but
do not specify its global finite normalization. The unrestricted
singleton-zero problem has the same missing input. A numerical row is
not an algebra presentation from which an image equation can be eliminated.

## What is already determined

The [residue ledger](RESIDUE_LEDGER.md) gives the three residual rows,
all omitted residue degrees one, and the exact pointwise formulas in
normalization branch counts. The [local proof](SIMPLE_RAMIFICATION_AND_RADIAL.md)
fixes the excess-one point as an ordinary cusp and gives every conditional
singleton radial packet. The [monodromy proof](GLOBAL_MONODROMY.md) gives
S6, the generic deletion rows and the absence of intermediate source fields.

None of these supplies the polynomial parametrization of the second
ramified image B, the optional unramified image H, their coincidences
and intersections with C in fixed target coordinates, or the cusp's
target value. Consequently the branch counts in the pointwise formulas
cannot yet be evaluated at all target points. The individual retained
residue extensions also remain unknown.

Using source A2 again in the finite-normalization Euler calculation
recovers the existing all-prime budget; no additional exclusion was
obtained in this audit. S6 alone does not determine those retained
residue extensions. Neither observation is a realization or an
impossibility proof for a residual row.

## Sufficient concrete input

Either of the following would make the missing calculations well posed.

1. An actual polynomial pair F=(P,Q) on the original source A2, in fixed
   target coordinates, with its constant nonzero Jacobian and degree-six
   function-field extension certified. For an intended terminal-r7
   application, also identify which proved terminal hypotheses apply.
2. A finite normal target algebra B together with the rational source
   coordinates x,y in Frac(B), a certificate that their simultaneous
   regularity open is the original source A2, and the target functions
   P,Q there. A finite cover with the correct numerical rows but without
   this source identification is insufficient.

The source-finite primitive inverse of [PF2OMIT3](ADMISSIBLE_INVERSE.md)
is an alternative computational interface to either package. It needs
the primitive-field, monicity/source-finiteness and normalization
identification certificates. A bare degree-six polynomial, its
discriminant, or the sign double cover does not provide that interface.

## Calculations once the input exists

For a supplied B and source coordinates, form the coordinate denominator
ideals (B:x),(B:y). PF2OMIT1 identifies their product's zero set with the
omitted boundary. Decompose its height-one support and contract each
prime to the target. For every prime record its image equation, omission
flag, valuation ramification index and residue-field degree separately.
Detect ramification on the entire normalization as well, so that a
ramification calculation cannot silently replace the boundary calculation.

For an admissible inverse R(T)=sum a_i T^i, factor its leading coefficient
to obtain the omitted image support. Codimension-one normalization and
Newton/residual data separate the primes above a shared image; a leading
coefficient multiplicity is a pole charge, not a missing-sheet count.
Specialized degree gives the aggregate retained length by PF2OMIT3.
It does not by itself separate colliding normalization points.

At exceptional values compute the local factors of the finite fibre
algebra B/(P-u0,Q-v0). Apply the source-coordinate regularity test to
each support point and record its local Artin length once. For a Keller
source the retained factors have length one. This distinguishes inverse
coordinate collisions, boundary self-image values and actual excess
length, including after the target blowup.

Completion requires actual image equations or certified parametrizations,
their residue data, an exhaustive exceptional-value description, and the
normalized local lengths there, or a proved exclusion of every remaining
row. The present null fields are intentionally retained until then.

## Current input audit

The current machine-readable checkpoints explicitly retain the missing
data:

- `inverse_equations.json`: `actual_r7_inverse_polynomial` and
  `actual_singleton_Keller_inverse_polynomial`;
- `residue_ledger.json`: `actual_ramified_image_B_equation`,
  `actual_unramified_image_H_equation`, `actual_exceptional_value` and
  `actual_Keller_inverse`;
- `global_monodromy_ledger.json`: every `image_equation` and individual
  retained residue list, the actual affine image equations and exceptional
  values, and the unconstructed global representation flags.

These are unresolved inputs, not unsuccessful computations on a supplied
candidate. The earlier optional request for a candidate path or equations
remains unanswered. This audit does not assert that no such candidate
could exist, or that further theoretical progress is impossible.

PF2OMIT9 subsequently uses the concrete marked Wright normalization to
[compute image restrictions and cusp tests](WRIGHT_IMAGE_RESTRICTION.md)
from any proposed finite target pair. It gives complete conditional
quartic value rows and a compulsory inverse-index zero at the cusp.
The target pair itself is still missing; these restrictions do not fill
the actual-image null fields above.

PF2OMIT10 adds [coefficient and index bounds for every singleton packet](SINGLETON_INVERSE_VALUATIONS.md).
It uses the completed local maps to constrain any proposed global inverse;
it does not determine the still-missing affine images or target values.

PF2OMIT11 now supplies [explicit first-image norm equations](IMAGE_NORM_LEDGER.md)
for both complete immersed r7 coefficient charts and computes all first-prime
losses by matrix corank. Its singular-value matching table constrains B,H
and excludes the H=C grouping at thirteen polynomial-chart parameters.
The coefficient parameter and original target coordinate changes, actual
B,H equations and global source normalization are still unresolved.
The earlier null equations are retained as historical checkpoints.

PF2OMIT12 closes the subsequent Laurent rank<=2 gate with an
[exact quartic-cofactor certificate](NO_FOUR_BRANCHES.md). It bounds
every image-parametrization fibre scheme by3, including nonimmersed
parameters, but does not supply B,H or the missing global normalization.

PF2OMIT13 adds [the nonimmersed first-image row](NONIMMERSED_LEDGER.md),
including exact normalized cusp coordinates in both classes and the
matrix corank correction there. The only possible singleton in that
row is the known cusp, conditional on the unknown B passing through it.
The extra cusp on unknown B in the immersed row remains unresolved;
the original target coordinate changes and global source map are absent.

The Laurent affine collision taxonomy is now complete by PF2R7KCOLL1/2.
PF2OMIT14 [reconciles both full catalogues and adds generic-line constraints](GENERIC_LINE_ACCOUNTING.md):
other image degrees are multiples of5; compact generic target/radial
line pullbacks have genus>=3, and removed points, lost sheets and
projection Euler jumps are recorded separately. These do not determine
the actual additional images or justify the old genus-two pencil.

PF2OMIT15 adds an [odd-degree affine intersection obligation](ODD_IMAGE_INTERSECTIONS.md):
for an actual additional image of degree5k with k odd, the reduced
intersection norm with C has degree in[9k,30k]. Shared-image row A
and row C therefore require singleton values; their positions still
need the actual B. Tangency weight is not a count of distinct values,
and the even-k matching cases remain unrestricted by this parity argument.

PF2OMIT16 supplies the [affine pole module and Milnor budget](DISJOINT_IMAGE_GATE.md)
used in the universal result below. Its earlier two-fibre reduction
and pointwise trace witnesses are retained as proof inputs, rather than
an active finite-parameter search. Connected singular support does not
invoke the smooth-set Euler theorem or supply a global source.

PF2OMIT17 closes the [universal disjoint-image gate](UNIVERSAL_IMAGE_INTERSECTION.md).
Every additional polynomial image meets C, so all residual supports
are connected. Shared-image A and cyclic C necessarily have singleton
values; the omitted unramified point in A remains in each packet.
No finite nodal trace gate remains. The actual images, exceptional
values and global finite normalization/source identification are still missing.

PF2OMIT18 determines the [nonlinear image-pencil geometric groups](IMAGE_PENCIL_MONODROMY.md):
S6 in A/B/N, A6 in C, with one geometric source component in all rows.
The generic genus/puncture ledger keeps omitted unramified points,
but needs the actual B/H intersection degrees and its own infinity
partition. The central H_C=0 accounting is supplied below; actual factors still need the image incidences.

PF2OMIT19 gives the [central image-fibre ledger](CENTRAL_IMAGE_FIBRE.md).
Euler/Milnor formulas use distinct normalization roots and actual
singular incidences, forcing source infinity Milnor jump at least46.
Normalized cyclic/cusp points retain their omission flags even when
unramified. Clean-stratum residue partitions and split-section lift
conditions constrain proposals; actual incidence, gluing, source infinity
singularities and global normalization/source remain missing.

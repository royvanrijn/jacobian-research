# Relative Frobenius norms and density-one ordinary reduction

For abelian varieties of adjoint type \(A_1\), the manuscript determines
the possible Newton slopes at every good degree-one prime outside a finite
set after a finite base extension. It computes the multiplicity of slope
\(1/2\) from the local trace factors and proves that the \(p\)-rank exceeds
half the dimension on a density-one set. At most two Hodge-active
geometric factors per rational simple adjoint factor force ordinarity on
that set. A twelvefold with two mixed tensor blocks illustrates this case.
For rationally simple adjoint groups with irreducible CM coefficient
blocks and at most two active embeddings, the paper identifies the
least ordinary-density-one extension and computes the exact ordinary
density over the original field.

The paper also proves density-one ordinary reduction over the original
field for signature \((2,2)\) fourfolds whose full geometric rational
endomorphism algebra is imaginary quadratic. It includes an explicit
genus-four Jacobian. A sharp finite-field bound in higher even rank
determines the range where Pfaffian parity guarantees integrality of the
norm at nonordinary primes exchanging the coefficient embeddings.

The simultaneous relative coefficient norm bounds give the exact
minimum \(p\)-rank of split numerical polygons in every signature.
Projective linkage gives density-one ordinarity for a fixed variety
with three balanced pairs.

[Manuscript PDF](cm-newton-polygon-constraints.pdf) ·
[LaTeX source and replay ZIP](cm-newton-polygon-constraints-source.zip) ·
[Metadata](metadata.json) ·
[Verification record](verification.json) ·
[Checksums](SHA256SUMS.txt)

The research records give the
[adjoint trace and determinant-parity proof](../../research/foundations/extended-geometry/CM_NEWTON_ADJOINT_PROJECTIVE_TRACE_BOUNDS.md),
[two-active-embedding construction](../../research/foundations/extended-geometry/CM_NEWTON_MULTIPLE_ACTIVE_TENSOR_INCIDENCE.md#two-shimura-surfaces),
[fourfold Pfaffian and explicit Jacobian](../../research/foundations/extended-geometry/CM_NEWTON_RANK_FOUR_PFAFFIAN.md),
[extremal polygon and sharp minimum proofs](../../research/foundations/extended-geometry/CM_NEWTON_EXTREMAL_POLYGONS.md),
[three-pair global argument and its residual endomorphism types](../../research/foundations/extended-geometry/CM_NEWTON_THREE_BALANCED_GLOBAL.md),
and [primary-source comparison](../../research/foundations/extended-geometry/CM_NEWTON_POLYGON_LITERATURE.md).
The [companion moment manuscript](../cm-newton-cyclic-cover-moments/)
develops nonvanishing and coherent reconstruction for cyclic covers.
Its [block-rank applications](../../research/foundations/extended-geometry/CM_NEWTON_BLOCK_RANK_MOMENTS.md)
bound the \(a\)-number at specified \(p\)-rank-zero reductions
and recover a two-dimensional Frobenius block from moments.

## Build and replay

From the extracted ZIP root, run

    SOURCE_DATE_EPOCH=1795651200 FORCE_SOURCE_DATE=1 latexmk -pdf -interaction=nonstopmode -halt-on-error -jobname=cm-newton-polygon-constraints main.tex
    python3 -I verification/verify_cm_newton_endpoint_bounds.py --bound 200
    python3 -I verification/verify_cm_newton_vertex_census.py
    python3 -I verification/verify_cm_moment_cover_faces.py
    python3 -I verification/verify_cm_newton_infinite_moment_families.py
    gp -q verification/verify_rank_four_pfaffian_curve.gp

The LaTeX sources use the standard AMS and Hyperref packages. The Python
replays use the standard library. The [dimension-200 summary](../../research/foundations/extended-geometry/cm_newton_endpoint_census_200.json)
records aggregate counts and the SHA-256 of the complete exact minimum
table. Pass `--bound 200 --all-minima` to stream its rows. A separate
integral-breakpoint replay recovers the full polygon multiplicities and
minima for all 1,188 historical signatures through dimension 30. The
PARI/GP replay verifies the finite-field polynomial and endomorphism
certificate used for the explicit genus-four Jacobian.

The minimum formula concerns numerical polygons. The finite-field
sharpness result varies the variety with the residue prime. The adjoint
theorem concerns primes over a finite extension of the original base.
The density of one flat factor at completely split primes with three
independent active labels remains open in the stated cubic setting.

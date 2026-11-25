# Frobenius norm bounds for Newton polygons with CM field action

The manuscript proves simultaneous bounds for relative Frobenius
coefficient norms and an exact closed formula for the least numerical
p-rank in every signature. It determines the four ordinary-forcing
classes, proves density-one ordinarity under projective linkage with
three balanced pairs, and gives infinite cyclic-cover moment exclusions.

[Manuscript PDF](cm-newton-polygon-constraints.pdf) ·
[LaTeX source and replay ZIP](cm-newton-polygon-constraints-source.zip) ·
[Metadata](metadata.json) ·
[Verification record](verification.json) ·
[Checksums](SHA256SUMS.txt)

The research records give the
[extremal polygon and sharp minimum proofs](../../research/foundations/extended-geometry/CM_NEWTON_EXTREMAL_POLYGONS.md),
the [three-pair global argument and its residual endomorphism types](../../research/foundations/extended-geometry/CM_NEWTON_THREE_BALANCED_GLOBAL.md),
the [infinite moment families](../../research/foundations/extended-geometry/CM_NEWTON_INFINITE_MOMENT_FAMILIES.md),
and the [primary-source comparison](../../research/foundations/extended-geometry/CM_NEWTON_POLYGON_LITERATURE.md).

## Build and replay

From the extracted ZIP root, run

    SOURCE_DATE_EPOCH=1795564800 FORCE_SOURCE_DATE=1 latexmk -pdf -interaction=nonstopmode -halt-on-error -jobname=cm-newton-polygon-constraints main.tex
    python3 -I verification/verify_cm_newton_endpoint_bounds.py --bound 200
    python3 -I verification/verify_cm_newton_vertex_census.py
    python3 -I verification/verify_cm_moment_cover_faces.py
    python3 -I verification/verify_cm_newton_infinite_moment_families.py

The LaTeX sources use the standard AMS and Hyperref packages. The Python
replays use the standard library. The [dimension-200 summary](../../research/foundations/extended-geometry/cm_newton_endpoint_census_200.json)
records aggregate counts and the SHA-256 of the complete exact minimum
table. Pass `--bound 200 --all-minima` to stream its rows. A separate
integral-breakpoint replay recovers the full polygon multiplicities and
minima for all 1,188 historical signatures through dimension 30.

The minimum formula concerns numerical polygons. The finite-field
sharpness result varies the variety with the residue prime. The
three-pair global theorem concerns a fixed characteristic-zero variety
with projective linkage. The density question in the remaining
independent-monodromy cases stays open.

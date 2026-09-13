# Latent Mordell--Weil lattice calibration

## Historical calibration

This programme is parked. Its frozen controls recover the known Fermigier
rank-12 primitive closure and a non-oracle R17 rank-16 component before exact
completion. Two independently frozen target sets then failed their own
dimension-recurrence gates. Those failures do not exclude another common
primitive subgroup or identify a new generic family.

The reusable method is [relation-component recovery](../../knowledge/ALGORITHMS.md#method-ec-relation-component-recovery): use exact relation components,
finite quotient ranks, primitive closures and held-out replay to select a
candidate; use exact arithmetic to verify it. Height shells, numerical scores
and selector dimensions remain proposal tools, not generic-rank statements.

The concise [reverse-engineering report](LATENT_LATTICE_REVERSE_ENGINEERING_REPORT.md)
states the control and target boundaries. The complete historical calibration,
all command interfaces and the detailed artifact catalogue are preserved
unchanged in the
[archive](../../archive/elliptic-curves/notes/LATENT_LATTICE_CALIBRATION.md.txt).

The retained [control](../../artifacts/generated-results/elliptic-curves/latent_lattice_graph_walk_calibration_v1.json),
[target-failure](../../artifacts/generated-results/elliptic-curves/latent_lattice_wgxli_frozen_dimension_v1.json)
and [secondary-failure](../../artifacts/generated-results/elliptic-curves/latent_lattice_secondary_frozen_dimension_v1.json)
artifacts remain sufficient to identify the recorded scope. A new target needs
an independently frozen protocol with exact endpoints before any tuning or
source-surface reconstruction.

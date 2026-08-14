# Elliptic-curve reference artifacts

This subdirectory contains compact pinned outputs for the separate
[elliptic-curve programme](../../../elliptic-curves/README.md).  Raw searches
and restart checkpoints belong under the ignored
`artifacts/local/elliptic-curves/` directory.

The JSON manifests are generated and replayed by commands in the programme's
[reproduction catalogue](../../../elliptic-curves/REPRODUCE.md).
`crt_lattice_calibration_v1.json` is an exact low-rank mechanism calibration.
`fermigier_crt_seed_v1.json` is an exact local seed in a high-generic-rank
family, but deliberately has neither a global-conductor result nor a rank
certificate.  Neither of those two artifacts is evidence that a target has
been reached.
`fermigier_rank_certificates_v1.json` exactly certifies twelve generic section
differences and the 22 published E22 points; it proves rank lower bounds but no
upper bound or saturation statement.
`fermigier_rank20_near_miss_v1.json` pins a bounded point search at
`u=28917/20`, an exact 20-point independence certificate, and an exact
conductor below the programme cutoff.  It remains one point short of the
rank-at-least-21 target.

# A 48 term normalized Keller counterexample

The known degree-seven map has a normalized collision frame with **48
nonlinear monomial occurrences**, improving the previous 49-term form.
The three coordinate counts are (6,19,23).

The current support interval is **[10,48]** for normalized three-variable
Keller collisions of raw coordinate degree at most seven. The exact minimum
remains open.

- [Explicit map, proof and two sharp frame minima](PROOF.md)
- [Polynomial coefficients](normalized-48-witness.json)
- [Independent verification](validation.json)
- [Independent Singular replay](singular-validation.json)
- [Discarded field-specialization implementation](INVALIDATED_SPECIALIZATION.md)
- [Specified nonzero-target endpoint obstruction](NONZERO_TARGET_LIMITS.md)
- [Remaining constructive directions](NEXT_FAMILY.md)
- [Scoped metadata and navigation checks](integration.json)

The value 48 is sharp in two specified families: all transverse bases in the
y,z-plane at the chosen collision pair, and all boundary collision pairs
with a lower triangular transverse completion. These do not cover every
affine frame or every Keller map.

From the repository root:

```sh
.venv/bin/python research/jacobian-followup/collision-frame-search/verify.py
```

The default independently reconstructs the map and parameter inputs and
audits the completed exact algebra replays. Full regeneration commands
are in the proof. Bounded beam exploration was used for discovery only.


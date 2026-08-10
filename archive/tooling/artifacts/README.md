# Archived tooling artifacts

These outputs accompany the exploratory searches archived one directory up.
They are not active certificates or inputs to `MATH_STATUS.json`.

The commands embedded in the JSON/text files record the paths used when the
artifacts were generated.  The corresponding scripts now live in
`archive/tooling/`; replay them from the repository root with that path.  The
Malle JSON files were moved byte-for-byte, so their adjacent SHA-256 sidecars
remain valid.

`weighted_seed_scan.json` is the exploratory output of
`archive/tooling/scan_weighted_seeds.py`; it was removed from the active
generated-results directory because it has no theorem-status role.

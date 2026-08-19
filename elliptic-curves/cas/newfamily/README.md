# Newfamily CAS recovery area

This directory is reserved for deterministic replay code for the six-root quartic construction documented in:

- `../../notes/NEWFAMILY_QUARTIC_ROOTS_AND_CONSTANT_SECTIONS.md`
- `../../notes/NEWFAMILY_REPLAY_CHECKLIST.md`
- `../../families/newfamily_symmetric_quartic_roots.json`

## Immediate recovery step

The exploratory work was performed largely in `/tmp` on macOS.  Preserve the surviving scripts verbatim before refactoring them:

```bash
mkdir -p elliptic-curves/cas/newfamily/archive

for f in /tmp/newfamily_*.py /tmp/newfamily_*.cpp /tmp/degree10_*.py; do
  [[ -e "$f" ]] && cp -v "$f" elliptic-curves/cas/newfamily/archive/
done
```

The archive is evidence/provenance, not the final replay interface.

## Intended canonical replay files

The archival scripts should eventually be reduced to these deterministic entry points:

```text
verify_newfamily_construction.py
verify_rich_seed_rank9.py
verify_high_automatic_rank_families.py
verify_constant_section_quadratic_orbit.py
verify_constant_section_degree6_orbit.py
verify_constant_section_degree10_orbit.py
```

Each canonical replay should be self-contained enough to run from the repository root with Sage, print explicit PASS/FAIL markers for its theorem-strength identities, and avoid relying on `/tmp` files.

## Evidence discipline

Keep raw scan logs and large tables under ignored `artifacts/local/elliptic-curves/newfamily/`.  Promote only compact deterministic manifests/certificates into `artifacts/generated-results/elliptic-curves/`.

Do not infer a full Mordell--Weil rank from the rank of the displayed section subgroup.  Do not promote the final degree-ten constant-section classification until its rational-point closure is replayed by a canonical script.

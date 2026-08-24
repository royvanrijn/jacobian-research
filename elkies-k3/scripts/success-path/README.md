# H3 equation-lift success path

This directory is the stable operational entry point for the successful H3
equation-lift route.  Canonical mathematical implementations remain one level
up in `elkies-k3/scripts/`; they are not duplicated here.  `ledger.json` pins
their SHA-256 hashes, exact output statuses, artifact hashes, proof boundaries,
and all pending corridor stages.

List or launch a stage with:

```bash
python3 elkies-k3/scripts/success-path/run_stage.py --list
python3 elkies-k3/scripts/success-path/run_stage.py q42_i8star_physical_marking
```

Verify the scripts and already-produced artifacts without rerunning the long
Sage calculations:

```bash
python3 elkies-k3/scripts/success-path/verify_ledger.py
```

The zero-pole, fast-q6, and identity-halving scripts are retained under
`shortcut_audits`, not as steps in the main route.  They pin three distinct
boundaries: the exact identity shell, the exact q6 degrees `435/703` that
invalidate rational-point transport, and the modular absence of an A11 chord
among the rational halving candidates.  The active construction proceeds by
resolved-surface Riemann--Roch.  A main-route stage may be marked complete here
only after its declared exact artifact and terminal status pass.

This is an operational ledger only.  `MATH_STATUS.json` is the sole
mathematical-status authority, and `STATUS.md` must continue to be generated
from it.

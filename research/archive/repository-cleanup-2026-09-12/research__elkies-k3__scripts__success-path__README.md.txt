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

## Current handoff

- The main equation route is locked through the marked `A11/MW6` child.
- Zero-pole, fast-q6 and identity-halving scripts are shortcut audits, not main
  route stages. They record the exact shell, the bad transport degrees
  `435/703`, and the modular absence of an A11 chord.
- The orbit42 resolved-RR artifact proves the exact A11 child over `QQ`.  The
  exact identity-shell degree fingerprint at the pinned good prime `100003`
  binds it to equation-side orbit64/mapping7; the orientation match retains
  this explicit good-reduction boundary.
- The next exact gate is the q8 lift from that marked A11 frame to `2A5/MW7`.
  `certify_h92_q24_a11_q8_construction_fingerprint.sage` transports the
  historical orbit922 divisor through every root/MW/glue-compatible
  fibre-preserving marking.  It narrows the equation-side nef targets to
  orbits 12 and 2162 and selects orbit12 by the declared minimum-MW-L1 rule.
  The transported divisor still has the historical formula `O+P-2F` and zero
  vertical-root correction.  The exact target-coset follow-up finds the shell
  index five, rejects the old pole-order-four generator as the wrong coset and
  selects `M=(1,0,0,0,0,1)` with
  `P12=M+S6-2*S2-2*S8`.  Coordinates for `M`, hence the characteristic-zero
  section and pencil, remain open.

Simple rule: an output saying `PASS` is evidence, not status. A main-route stage
is complete here only when the declared artifact, terminal status, hashes and
repository status authority all agree.

This is an operational ledger only.  `MATH_STATUS.json` is the sole
mathematical-status authority, and `STATUS.md` must continue to be generated
from it.

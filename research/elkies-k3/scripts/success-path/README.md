# Historical H3 equation-lift prefix

The [ledger](ledger.json) preserves the original q24/orbit42-to-A11 prefix,
its scripts, hashes, exact outcomes and shortcut failures. Its formerly
pending q8 stage is historical. The later q8/orbit376 and q12/orbit5867
endpoint, rational marking, Picard19 and saturated MW17 proof are complete;
use the [current route record](../../ORBIT42_EQUATION_LIFT.md) and
[K3 programme](../../README.md).

The historical target-coset step found index five, rejected the old pole-four
generator and selected the pole-five class. Its coordinates no longer need
discovery for the completed corridor. Keep the retained zero and marking
corrections when replaying this prefix.

List the historical records without running Sage:

```sh
python3 elkies-k3/scripts/success-path/run_stage.py --list
```

`verify_ledger.py` checks the historical hashes and available artifacts. Its
success does not assert a current frontier or an independent theorem replay;
a missing old local file does not authorize rebuilding the chain. `run_stage.py`
still accepts an explicitly selected retained stage for reproduction.

The [original handoff](../../../archive/repository-cleanup-2026-09-12/research__elkies-k3__scripts__success-path__README.md.txt)
and original ledger are preserved. Current mathematical status belongs only
to [MATH_STATUS.json](../../../MATH_STATUS.json).

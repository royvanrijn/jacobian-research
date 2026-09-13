# One-command V3 → transfer autorun

The detached controller in `elliptic-curves/cas/run_v3_transfer_autorun.py` chains the already-running V3 final verification into the fixed cross-parent transfer campaign.

From `research/`, launch exactly once:

```sh
python3 elliptic-curves/cas/run_v3_transfer_autorun.py launch
```

The command returns immediately. The detached controller then:

1. waits for the existing V3 supervisor to publish `replay-M17.json`, `metric-replay-M17.json`, and the final V3 package;
2. requires the independent replay to certify rank at least 31;
3. prepares the frozen determinant-948/native11952 transfer roster once;
4. runs the positive control and its independent replay;
5. releases the three warm rank-27 jobs only if that control gains rank;
6. executes each released case through the existing bounded `v3_transfer_campaign.sage next` supervisor;
7. stops after the fixed four-case roster, after a closed positive-control gate, or on the first preserved method/resource failure.

It never retunes V3, expands a parameter/family roster, retries failed budgets, or writes mathematical status. The existing V3 and transfer evidence remain authoritative.

Inspect without disturbing it:

```sh
python3 elliptic-curves/cas/run_v3_transfer_autorun.py status
tail -f artifacts/local/elliptic-curves/v3-transfer-autorun-v1/autorun.log
```

The controller has a fixed 12-hour wait for the upstream V3 final package. Expiry is an operational stop requiring review, not arithmetic evidence.

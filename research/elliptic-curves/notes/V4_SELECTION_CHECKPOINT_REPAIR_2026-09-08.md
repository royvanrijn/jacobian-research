# V4 selection checkpoint round-trip repair

The `saved V4 selection differs` startup failure is reproducible without a
point search: `ExactParity.solve` returns `minima` as a list of tuples. JSON
serializes both tuples and lists as arrays and reads them back as lists.
Consequently V4's first write followed by its own equality check fails even
when all minima, centres and their order are unchanged. Independent replay
had the same representation mismatch.

V4 now converts a selection to an exact JSON-shaped value before publication
AND before returning it. Only list/tuple representation is normalized; integer
values, coordinate strings, signs, minima ordering and centre ordering are
unchanged. Non-string dictionary keys and unexpected scalar types are rejected.
Genuine mismatches report the first differing field/index. The shared exact-CVP
solver, roster, SHA domain, selector lanes, search bounds and proof checks are
not modified. Preflight now includes the save/read/compare boundary.

## Recover this zero-chart failure

From the repository root, with the old controller stopped:

```sh
git pull --ff-only &&
cd research &&
"${V4_SAGE:-$HOME/.local/bin/sage}" -python -m unittest discover \
  -s elliptic-curves/tests -p 'test_det1092_v4_checkpoint_roundtrip.py' -q &&
python3 elliptic-curves/cas/run_det1092_v4_bootstrap.py repair-resume --hours 24
```

The new `repair-resume` action runs under the controller lock and refuses a live
worker, a modified input, or ANY chart/result/unknown file. It archives the old
prepared V4 directory under `det1092-v4-wide-bootstrap-v2-failed-preparations/`
with a file-hash manifest. The old logs/sessions remain in place. Existing
selection files are copied back byte-for-byte, not discarded or silently
rehash-approved. Fresh preparation binds the repaired source; preflight then
recomputes and checks those selections before point search. The active V4
paths stay unchanged; no V5/new family/version is introduced.

Use ordinary `resume --hours 24` for subsequent resource interruptions. This
special repair is deliberately unavailable after any chart has run. Progress
remains `python3 elliptic-curves/cas/run_det1092_v4_bootstrap.py status`.

## Validation

29 new tests pass in the coding environment, including the real two-dimensional
CVP output, first/repeated publication, unchanged legacy bytes, exact large
integers, deliberate corruption, and mocked preflight/search/replay/interrupt
flows using actual JSON files and worker control flow. Recovery tests verify
byte preservation and refusal of live workers, chart evidence and bad inputs.
The mocked arithmetic is not a mathematical certificate. Native Sage/PARI,
the actual 512-class selection and the user's retained job were not executed
in this environment; those remain runtime verification gates. No rank or
mathematical-status claim is changed.

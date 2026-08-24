# AGENTS.md

This directory inherits the repository rules. `../MATH_STATUS.json` is still
the only authority for what is proved.

## Status in plain language

- One K3 surface can have many elliptic fibrations. Changing the fibration
  moves rank between fibre roots and Mordell--Weil sections; it does not create
  new Neron--Severi classes.
- The generic rank budget is `root rank + MW rank = 17`.
- Equations are certified through `D12/MW5`. The marked lattice route from H3
  to the pinned rootless `MW17` endpoint is exact.
- The next edge is `D12/MW5 --q6 orbit42--> A11/MW6`. A local exact-RR artifact
  currently reports `PASS`, but do not call the edge proved until its canonical
  note, software lock, and `MATH_STATUS.json` entry agree.

## Start here

1. Read `ELKIES_K3_PROCESS_ATLAS.md` for the chronology and reusable lessons.
2. Read all of `RANK_MUTATION_AND_LIFT_THEOREMS.md` before starting a new
   search, lift, route change, or generalization.
3. Read `README.md` for the current frontier.
4. Use `scripts/success-path/ledger.json` for pinned commands and hashes.
5. Check `MATH_STATUS.json` before repeating any claim.

## Maintain the theorem package

For every new exact behavior, first decide whether it is:

- an instance of an existing theorem;
- a stronger theorem with new hypotheses or consequences;
- a counterexample showing that a hypothesis is missing; or
- only an experimental pattern.

Update `RANK_MUTATION_AND_LIFT_THEOREMS.md` when a reusable statement, proof,
proof obstruction, or counterexample is found. Keep the change small and add:

1. the exact hypotheses;
2. a proof or a clearly named missing lemma;
3. the mathematical consequence for the navigation/compiler engine;
4. links to the canonical certificate, note, and primary literature;
5. a warning when a negative result is only bounded or route-specific.

Do not turn one successful calculation into a universal theorem. Do not weaken
an existing hypothesis silently; record why it can be removed. If a result
changes repository status, update `MATH_STATUS.json` and regenerate `STATUS.md`
only after the canonical proof and replay are ready.

## Useful hints

- A fibration is identified by its full marking (`F`, `O`, components,
  sections, chamber and NS transports), not just by an ADE/MW label.
- Saturate glue and resolved local modules before trusting ranks or condition
  counts.
- Track every cover/isogeny degree and every rational-function denominator.
- Keep generic and CM-specialized nodes separate.
- Preserve negative experiments: they are routing information, not theorems.

## Known traps

- Do not reuse the old q8 double-2-cover or missing-`Dx` normalization.
- Do not use the old orbit42 correction-one or `P.O=2` profile; the selected
  class has correction `3`, `P.O=3`, and divisor `D42=O+P+V`.
- The proposed fast q6 point transport has degrees `435/703`, so it is not a
  section transport.
- The rational-halving scan found no A11 chord, but that is not a global
  non-existence theorem.
- Do not reopen q24/orbit85, zero-pole, or archived Hensel searches unless a
  new question explicitly requires a regression.

For documentation or ledger changes, run:

```bash
python3 elkies-k3/scripts/analyze_process_ledger.py --check-document
python3 elkies-k3/scripts/success-path/verify_ledger.py
```

Long Sage calculations are separate targeted replays; report them if skipped.

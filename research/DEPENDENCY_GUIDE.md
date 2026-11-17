# Maintaining mathematical implications

The graph lives inside [MATH_STATUS.json](MATH_STATUS.json), alongside the existing
principal results. [Scope and gaps](SCOPES.md), [dependency diagrams](DEPENDENCIES.md),
[status](STATUS.md) and topic guides are generated views. There is one authored
status record.

Use this graph to trace the consequences of a new lemma, locate the unresolved
premises of a proposed theorem, or find claims affected by a changed proof. The
initial graph covers principal results and identified gaps across every topic.
Most source proofs are bundled leaves with **partial dependency coverage**: their
internal lemmas and imported theorems have not all been expanded into separate
nodes. A green node means its declared route is supported by current, internally
accepted evidence. It does not certify that the argument is correct or complete.

## Read the graph

A node states one exact proposition. State its category, base field, characteristic,
dimension, quantifiers, parameter restrictions and exceptional cases in its scope.
Our own results and existing research can both supply leaves. An imported theorem
needs its precise version and theorem locator, plus the hypotheses used here.

An implication has a nonempty list of premises and one conclusion:

```text
A AND B --proof of the implication--> C
C       --proof of the implication--> D
```

All premises of one rule are required. Different incoming rules are alternative
routes: any valid route can support the conclusion. An arrow must have its own
proof; related subject matter or an attractive research plan does not establish
an implication. Equivalent statements need one acyclic direction in this graph,
with the equivalence documented in the source, rather than a dependency cycle.

The engine distinguishes these effective states:

| State | Meaning |
|---|---|
| proved | A current declared proof route has all premises and no unresolved assumptions. This is an internal evidence classification. |
| conditional | A route is supported under its listed assumptions. |
| open | No currently available route establishes the exact statement. |
| stale | A supporting source changed or disappeared, blocking the available route. An unaffected alternative can still establish the node. |
| falsified | The exact recorded proposition has a counterexample; it cannot support a positive premise. |
| unmapped | A principal result has no declared proof route. Validation rejects this for the current catalogue. |

Proving a counterexample theorem is a positive result about a *negated* conjecture.
Give that theorem its own scope and node. A falsified premise disables a rule; it
does not prove that rule's conclusion false.

Dependency coverage (`complete`, `partial`, `unmapped`) describes how fully the
proof dependencies have been modelled. It is separate from effective state.
The authored principal-result state also stays visible: a partial programme may
contain proved restricted theorems. Its graph status concerns only its stated
restricted scope.

Formal verification, independent replay, external review and publication retain
their existing records and scope. They never propagate along mathematical arrows.
In particular, an imported Lean axiom can leave the formal theorem conditional
even when a separate written mathematical proof is available.

## Inspect consequences and gaps

Run these commands from the repository root:

```sh
python3 research/scripts/research.py gaps
python3 research/scripts/research.py gaps --topic hessian --kind mathematical
python3 research/scripts/research.py explain HESSIAN-MOMENT-TRACE
python3 research/scripts/research.py impact HESSIAN-MOMENT-TRACE
python3 research/scripts/research.py graph --format dot > /tmp/research-dependencies.dot
python3 research/scripts/research.py graph --format json > /tmp/research-dependencies.json
```

Use `show ID` to inspect a raw node and `explain ID --json` for recursive proof
provenance, exact source pins, all routes and blockers. `impact` reports structural
descendants, including currently blocked alternatives; it does not assert that
every listed claim must lose its support. Use `impact --source PATH` with an
evidence path relative to `research/` to trace a source change. The diagrams omit
fresh unconditional proof-source leaves for readability and retain their count
on each rule diamond. The JSON and DOT exports include the whole graph.

Gaps are classified as mathematical, formalization, verification or workflow.
Review and publication tasks do not become mathematical assumptions. A full
conjecture can remain open with no known sufficient reduction in the graph;
do not invent an arrow merely to connect every goal.

## Try a hypothetical improvement

Choose an unresolved node from `gaps` and substitute its identifier for `ID`:

```sh
python3 research/scripts/research.py simulate --prove ID
python3 research/scripts/research.py simulate --prove ID --json
python3 research/scripts/research.py simulate --stale ID
python3 research/scripts/research.py simulate --reopen ID
```

For a current example, run `simulate --prove FACTORIAL-GAP-FULL-RANK` to inspect
the consequence for the lower-one higher-gap family.

Multiple overrides may be supplied. Simulations are read only. An assumed proof
leaves `ID` among downstream assumptions, so conclusions remain explicitly
hypothetical and conditional. A simulation does not replace a proof or change any
assurance flag.

## Record an actual proof

1. Check the exact statement and every hypothesis against the canonical argument.
   Follow [the mathematics guide](../working-guide/MATHEMATICS.md), including
   imported results and the scope of any finite or formal check.
2. Add or update the canonical proof and retained evidence. For a new result,
   create a narrowly scoped node. Keep a wider unresolved goal separate.
3. Update its graph node or incoming implication in `MATH_STATUS.json`. A proved
   leaf requires `sources`; a proved implication requires `proof`. Each record
   contains a repository-relative `path` (relative to `research/`), the file's
   64-character `sha256`, and a precise `locator`. Multiple indispensable source
   files can be listed on a leaf. An implication whose proof imports additional
   prerequisites must model them as premises or explicitly record partial
   coverage and its remaining audit boundary.
4. Represent dischargeable assumptions as premise nodes. Change an open leaf to
   `proved` only after the proof is checked, or add a proved route into that node.
   Add both premises and the proof of their implication. Do not change a broad
   conjecture merely because a restricted family is now proved.
5. Update the topic's exact scope notes, gap mapping and classified open questions.
   Reassess its authored state separately; the renderer never rewrites it.
6. Run the relevant mathematical checks, then:

   ```sh
   make render-navigation
   make check
   ```

The engine recomputes every reachable result on each command. A completed premise
can support several levels of consequences without manually updating their graph
statuses. Generated pages are checked for drift, so CI rejects summaries left
behind after a dependency update.

For retained TeX and Lean sources, `make check` also detects unpinned local
input/import dependencies and retained Lean toolchain/configuration files. Add
them to the source leaf after checking their role; imported packages retain the
formal project's existing dependency pins. This closure check covers literal
local imports, not every possible macro or dynamically loaded dependency.

A source change invalidates its old pins and affected routes. Use `impact --source`
to identify consumers, inspect the change, rerun the necessary mathematical checks,
and update pins only after revalidating the claims. The renderer can display stale
evidence; `make check` rejects it. Hash replacement alone supplies no mathematical
evidence.

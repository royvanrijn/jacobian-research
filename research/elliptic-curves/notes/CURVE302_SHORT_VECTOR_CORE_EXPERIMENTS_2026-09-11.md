# Curve302 short-vector/core experiments

These experiments test the strongest interpretation left after the closure work:
the fourteen seeded V3 trajectories may be repeatedly recovering a distinguished
low-height **integral sublattice** of the displayed quotient
`D/M17 ~= Z^14`, rather than following a predictive adaptive gradient.

They consume only the sealed `curve302-closure-structure-v1` outputs
(`quotient-relations.json`, `trajectories.json`, `REPORT.json`). No point search,
Mestre scoring, V3 execution, rank search or new curve construction is launched.

Run from `research/`:

```sh
python3 elliptic-curves/cas/run_curve302_short_vector_core.py run
python3 elliptic-curves/cas/run_curve302_short_vector_core.py status
python3 elliptic-curves/cas/run_curve302_short_vector_core.py check
```

Use `--source` if the passed closure-structure evidence is not in the canonical
local folder. The default output is
`artifacts/local/elliptic-curves/curve302-short-vector-core-v1/`.

## Exact performance fixes after the first real run

The first attached run stopped without a receipt and is retained as UNKNOWN.
The detached `curve302-short-vector-core-v2` run completed and sealed exact
enumeration: **1,288,441 primitive directions**, 5,789,129 nodes, covering all
180 acquisitions. Filtration then reached its 3,600-second cap; its timeout
receipt and the complete enumeration remain preserved. No basin result or
completed full check came from that run.

The optimized implementation retains the same mathematical definitions and
limits:

- Incremental rational echelon reduction replaces repeated general matrix
  rank calls. Filtration consumes the complete final norm shell, then stops
  after reaching rank14: saturation is Z^14 at every later shell.
- A single shell scan computes all observed rank intervals with exact ties.
- Basin annihilators are cleared to primitive integer rows. Optional NumPy
  int64 dot products run only when the exact bound
  `sum(abs(a[j])*max(abs(v[j]))) <= 2^63-1` proves every intermediate fits.
  Unsafe cases and environments without NumPy use arbitrary-size Python
  integers. No floating-point containment or cutoff is introduced.
- Enumeration obtains each leaf's norm from the exact accumulated LDL sum,
  instead of evaluating all196 matrix terms again. Full deterministic `check`
  still re-enumerates the entire ball and compares the TSV byte-for-byte.

Seventeen tests pass, including comparison with exact matrix ranks and rational
basin tests, overflow fallbacks, exact ties and boundary vectors, a full
synthetic pipeline/check, successful enumeration reuse and rejection of a
corrupt donor before output creation.

From the repository root, the fresh optimized run is:

```sh
sage -python research/elliptic-curves/cas/run_curve302_short_vector_core.py run \
  --reuse-enumeration research/artifacts/local/elliptic-curves/curve302-short-vector-core-v2 \
  --folder research/artifacts/local/elliptic-curves/curve302-short-vector-core-v3
sage -python research/elliptic-curves/cas/run_curve302_short_vector_core.py check \
  --folder research/artifacts/local/elliptic-curves/curve302-short-vector-core-v3
```

`--reuse-enumeration` is for a fresh folder only. It verifies the donor plan,
stage seal, output hashes, identical source bindings, exact observed norm bound
and the current node/vector limits. Imported files are copied and rehashed;
the new plan and import receipt preserve their producer's code hashes. The
failed donor is never resumed, edited or relabelled. Sage Python is used here
because system Python lacks SymPy; the installed SymPy is1.14.0.

## Experiment 1: complete primitive-vector vocabulary

Let `Q` be the exact rational 14-dimensional Schur form stored by the previous
closure experiment (exact for its entrywise-rounded height metric). Compute

```
H = max Q(v)
```

over the 180 already acquired primitive quotient vectors, then enumerate **all**
primitive integer directions modulo sign satisfying `Q(v) <= H`.

Enumeration is exact Fincke--Pohst recursion using a rational LDL decomposition.
Integer interval endpoints are obtained with integer square roots; floating point
cannot exclude a boundary vector. `--max-directions` and `--max-nodes` are
fail-closed resource guards: hitting either makes the experiment UNKNOWN rather
than returning a truncated vocabulary.

This removes the 49/180 and 70/180 vocabulary-coverage limitation of the previous
next-move experiment. Every acquisition must occur in the complete vocabulary.

## Experiment 2: intrinsic saturation filtration

Order the complete primitive directions by exact `Q`-norm. At every norm shell
where the rational span grows, form

```
Lambda(h) = Sat_Z < v primitive : Q(v) <= h >.
```

Smith/Hermite arithmetic returns a canonical saturated integral basis. This gives
an intrinsic short-vector filtration depending only on the quotient height form.

Independently reconstruct every one of the **194** actual trajectory prefix
lattices. Each is required to have saturation index one. At every attained
quotient dimension, intersect the prefix lattices of all runs present there to
recover the exact common integral core. The experiment then asks when the
intrinsic `Lambda(h)` first contains, or exactly equals, each observed common
core.

A close match would support the statement that V3 is exposing the successive
short-vector saturation filtration of the exceptional quotient. A mismatch is a
clean falsification of that stronger explanation.

## Experiment 3: exact acquisition ranks and basin multiplicity

Every one of the 180 acquisitions receives its exact static-height rank interval
inside the complete primitive vocabulary, including exact norm ties. There are no
out-of-vocabulary events.

For each dimension at which the common integral core grows, inspect every actual
prefix one step earlier. For every complete-vocabulary short direction `v` not
already in that prefix, test whether

```
Sat(prefix + Z v)
```

contains the next common core. Since the inputs are saturated, this is exactly a
rational-span containment test. Report basin fractions both:

- through the actual next acquisition's norm; and
- through the experiment's complete global bound `H`.

This distinguishes a simple **energy mechanism** (the core is represented by a
few exceptionally short vectors) from a **multiplicity/basin mechanism** (many
different short extensions saturate into the same underlying core).

## Interpretation boundary

The quotient form is the already published rounded-height Schur form, not an exact
canonical-height matrix. The experiment is retrospective and uses the 180 known
acquisitions only to choose the complete enumeration bound and to measure the
observed common cores. It does not infer a V3 transition probability, prove that
a short lattice vector has a rational point/chart of bounded complexity, or prove
a self-propagation theorem.

The useful theorem-shaped next target, if the experiment is positive, is:

> A distinguished low-height saturated sublattice of `D/M17` has a large basin of
> short primitive extensions under the half-lattice geometry; once enough of that
> sublattice is acquired, later exceptional directions admit substantially cheaper
> half-lattice representatives.

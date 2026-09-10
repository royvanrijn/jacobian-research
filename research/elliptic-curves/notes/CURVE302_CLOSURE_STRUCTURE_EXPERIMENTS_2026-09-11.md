# Curve302 closure-structure experiments

This experiment asks a narrower question than another rank search: **why does almost any one of the fourteen displayed exceptional directions on curve302 unlock most or all of the other thirteen under adaptive half-lattice search?**

The fourteen directions are independent modulo the recovered generic `M17`, so the target is not an ordinary Mordell--Weil linear relation. The experiment instead looks for a closure law, short quotient near-relations, and repeated actual transition subspaces.

The single controller is:

```sh
sage -python elliptic-curves/cas/run_curve302_closure_structure.py prepare
sage -python elliptic-curves/cas/run_curve302_closure_structure.py run
sage -python elliptic-curves/cas/run_curve302_closure_structure.py status
sage -python elliptic-curves/cas/run_curve302_closure_structure.py check
```

The default evidence folder is `artifacts/local/elliptic-curves/curve302-closure-structure-v1/`.

## Experiment 1: exhaustive closure law

Consumes the immutable `2^14` exceptional-subgroup landscape. For the global minimax and historical-tail thresholds it evaluates the closure of every subset, extracts inclusion-minimal trigger hyperedges `S -> j`, all inclusion-minimal full generators, critical singleton-to-direction arrival thresholds, and checks the closure axioms exhaustively. It also tests matroid exchange and anti-exchange on closed states rather than assuming either structure.

This is exact combinatorics of the retained finite-atlas costs. It does not promote those costs to exact CVP, chart runtime, or rank probabilities.

## Experiment 2: quotient near-relations

Reconstructs the same 31-point rounded canonical-height metric used by the existing landscape and verifies its hash. It takes the Schur complement modulo the generic `M17`, giving a 14-dimensional positive-definite quotient form, then:

- computes the best pairwise discrete projected half-lattice relation;
- compares pairwise projected unlock against the actual retained-atlas unlock matrix;
- enumerates the shortest primitive support-2/3/4 quotient combinations;
- searches each small minimal closure trigger for a bounded symbolic near-relation of the form `2 e_j - sum k_i e_i`;
- transports the exact ten-dimensional strict Kummer kernel into the same 14-dimensional quotient coordinates.

All quotient norms are exact for the **entrywise 1e-6 rounded** height metric, not exact canonical-height identities. A short nonzero vector is a near-relation, not a Mordell--Weil dependency.

## Experiment 3: actual V3 trajectory reconciliation

Requires all fourteen completed raw seeded V3 runs: the two original seeded-amplifier runs plus the twelve-seed universality panel. `prepare` freezes the seed packet, terminal record, and every stage-final independent-point audit used by this experiment. Missing or changed raw evidence aborts before analysis.

For every new point acquired by V3, a high-precision height solve proposes coordinates in the fixed 31-point basis. The candidate is then checked by exact elliptic-curve group arithmetic; numerical recognition alone is never accepted. The resulting 14-dimensional quotient vectors are reduced mod 2 to obtain a canonical quotient-subspace signature at every attained rank.

The output therefore shows whether different starting seeds converge to the same intermediate quotient subspaces, which exact quotient-vector families recur across runs, and whether those vectors lie in the transported ten-dimensional strict Kummer subspace.

## Fail-safe boundary

The controller launches **no point search** and modifies no existing campaign. Inputs and local raw evidence are hash-frozen at `prepare`. Each experiment runs as a separately supervised stage with a receipt and output seal. A failed, censored, or unreceipted interrupted stage is retained as `UNKNOWN`; it is not silently rerun in place. Use a new folder after investigating such a failure.

The GitHub Actions smoke runs synthetic policy tests plus a full real-data recomputation of the fourteen singleton minimax thresholds from the committed 16,384-state landscape. CI deliberately does not fabricate the local raw V3 directories needed by experiment 3.

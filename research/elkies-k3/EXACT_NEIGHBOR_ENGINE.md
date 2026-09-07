# Reusable exact elliptic-neighbour lattice engine

[`scripts/exact_neighbor_engine.sage`](scripts/exact_neighbor_engine.sage)
turns a supplied primitive isotropic divisor and explicit section/component
data into an exact child fibration.  Its degree-independent entry point is
`degree_q_neighbor(ns, divisor, old_fiber, curves, expected_old_fiber_degree=q)`;
the stable compatibility wrapper
`degree_two_neighbor(ns, divisor, old_fiber, curves)` pins `q=2`.

```sage
load("elkies-k3/scripts/exact_neighbor_engine.sage")
result = degree_q_neighbor(
    ns, divisor, old_fiber, section_component_data,
    expected_old_fiber_degree=q,
)
child = result["minimized_child_frame"]
transport = result["transport"]
```

The result also records `reduced_divisor`, the ordered
`fixed_component_sequence`, the unminimized `child_frame`,
`child_root_data`, the reduced `child_mw_height`, a canonical row-major frame
digest, and the child simple-root classes lifted back into the parent
Neron--Severi basis.  `transport_parent_vector_to_child` and
`transport_marked_parent_vectors` preserve marked sections/components through
the neighbor; `lift_child_frame_vector` performs the inverse geometric lift.

The input curves are ordered `(name, class)` pairs for effective `(-2)`
curves.  The engine first removes every fixed supplied curve with negative
intersection, then verifies that the surviving divisor has the declared
positive old-fiber degree.  It splits a primitive hyperbolic plane integrally, returns the raw child
frame and its unimodular Neron--Severi transport, and makes the child ready for
another neighbor: it chooses lexicographic simple roots and LLL-reduces the
saturated Mordell--Weil quotient.

The root/MW adaptation is explicit about torsion.  For a primitive root
lattice it reports `ROOT_MW_MINIMIZED`; a rootless child receives a full-frame
LLL reduction.  For a nonprimitive root lattice it retains the exact child
and returns `PARTIAL_NONPRIMITIVE_ROOT_LATTICE`, its Smith invariants, and no
fictional saturated MW height.  A torsion/glue-aware quotient can then be
supplied by a downstream routine.

The output distinguishes what the calculation establishes from any external
nef argument.  In particular, `nonnegative_on_supplied_curves` means exactly
that, not that those curves exhaust all effective `(-2)` curves.  A caller
must retain its exact section/CVP and multisection proof when that is needed
to certify global nefness.  A proof reference can be preserved verbatim in
the certificate's `proof_metadata`.

## Serialized runner and certificates

The versioned input and output contracts are
[`../schemas/elkies_k3_exact_neighbor_input.schema.json`](../schemas/elkies_k3_exact_neighbor_input.schema.json)
and
[`../schemas/elkies_k3_exact_neighbor_certificate.schema.json`](../schemas/elkies_k3_exact_neighbor_certificate.schema.json).
The standalone runner accepts the input contract and writes the certificate:

```sh
sage elkies-k3/scripts/run_exact_neighbor_engine.sage \
  --input path/to/neighbor-input.json --output path/to/certificate.json
```

Certificate hashes are SHA-256 over canonical sorted JSON, while every frame
and transport has its own canonical row-major SHA-256.  This records the
actual Sage version without presenting its LLL representatives as intrinsic.

Build or verify the two pinned example certificates (the command refuses to
replace changed contents):

```sh
sage elkies-k3/scripts/build_exact_neighbor_engine_certificates.sage
```

The resulting artifacts are
[`../artifacts/generated-results/elkies-k3-exact-neighbor-q80-first-q4.json`](../artifacts/generated-results/elkies-k3-exact-neighbor-q80-first-q4.json)
(`b6ea4c8b421cf782bf57416935b20bb3424118c0531a236c4e66548bc07895c3`)
and
[`../artifacts/generated-results/elkies-k3-exact-neighbor-h3-d13-q24.json`](../artifacts/generated-results/elkies-k3-exact-neighbor-h3-d13-q24.json)
(`093a0e1b7fe8a1ef93cfffaa758762ae7e7ff83278ee74d439e1ff4ea052c01c`).

The q80 first `q=4` degree-two pencil is the regression example:

```sh
sage elkies-k3/scripts/verify_exact_neighbor_engine.sage
```

It reconstructs the known first child with root data `(13,164,20)`, i.e.
`D9+A4/MW4`, from the source `E6+D5+A3/MW3` frame.  The same checker then
feeds the pinned section/component data for the first rank-growing H3 arrow
`D13/MW4 --q24--> D12/MW5` to the engine and obtains `(12,264,4)`.
The existing
[`scripts/analyze_q80_rootless_first_neighbor.sage`](scripts/analyze_q80_rootless_first_neighbor.sage)
is the companion global-nef certificate; this engine reuses its wall layer but
does not replace its exact MW/CVP and bisection argument.

The complete H3 `D13`-to-rootless replay now calls the same engine for every
primitive `U`-split while retaining its pinned root-adapted bases and original
artifact digest:

```sh
sage elkies-k3/scripts/verify_h3_d13_to_mw17_path.sage
```

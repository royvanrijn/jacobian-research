# NS0024 edge-1 compiler preparation — 2026-09-01

## Outcome

The equation-independent preparation for the first preferred NS0024 edge is
complete.  The exact lattice/source marking identifies

```text
A3+A4+A6/MW4 --q4/orbit1, old-fibre degree 2--> A1+A2+A4+D5/MW5
```

with the minimum-pole source section `P3` and the literal divisor

```text
D = O + P3 + 2F - C2 - 2C3 - C4
```

on the normalized split `I5` fibre at `t=1`.  Here `P3` meets component 1.
The complete old-generic-fibre ambient is

```text
<1, t, t^2, m>,       m=(y+y(P3))/(x-x(P3)).
```

The new split-multiplicative toric adapter evaluates all three supported
components and both multiplicity layers on `C3`.  It produces the exact
condition matrix over the source coefficient field.  Edge 1 is admitted only
when that matrix has profile `4 -> 2`, the same kernel compiles to a degree-four
binary quartic, and the child Jacobian has

```text
I1* + I5 + I3 + I2 + 7I1,
root type D5+A4+A2+A1,
root rank 12, root determinant 120, Euler number 24.
```

This is compiler preparation, not an equation result.  No finite-field MW4
family was available to this work, so no NS0024 edge has yet been promoted at
equation level.

## What was fixed by the abstract marking

The source frame alone has four abstract MW quotient coordinates.  The exact
minimum-pole basis certificate replaces that ambiguity by the following
equation-facing basis:

| section | `P.O` | components on `I7,I5,I4` |
|---|---:|---|
| `P1` | 0 | `(1,0,0)` |
| `P2` | 0 | `(2,1,3)` |
| `P3` | 0 | `(2,1,1)` |
| `P4` | 1 | `(1,1,1)` |

The q4/orbit1 MW projection is `(-1,0,0,0)`, which is exactly the stored
quotient coordinate of `P3`.  Subtracting `O+P3+2F` from the selected fibre
leaves root-frame vector `(1,0,2,1,0,...,0)`.  In the effective `I5` chain
order `(1,0,2,3)` this is the vertical correction `(0,-1,-2,-1)`.  The new
fibre degrees on `C1,...,C4` are `(0,0,2,0)`, and the affine degree is zero.
There is no vertical correction on `I7`, `I4`, or any irreducible fibre.

The machine-readable replay is
[`../artifacts/generated-results/elkies-k3-lattice-foundry-ns0024-edge1-compiler-preparation.json`](../artifacts/generated-results/elkies-k3-lattice-foundry-ns0024-edge1-compiler-preparation.json).

## Construction methods reused

Three previous exact q4 mechanisms were compared before selecting the adapter.

1. The physical q4/orbit1584 calculation showed that a degree-two q4 divisor
   with a connected middle-double correction can reduce the standard
   four-dimensional chord ambient to a two-plane by one resolved local jet.
2. The q4_6855 characteristic-zero replay identified the middle-double module
   as value plus first-jet conditions along the analytic nodal centre.
3. The q4_1938 replay supplied the robust split-multiplicative toric model when
   the horizontal itself specializes to the singular node and its oriented
   resolved component matters.

NS0024 needs the third geometry: `P3` meets the `I5` node in component 1, so a
coarse node value or ADE lookup cannot certify the module.  The adapter uses
the analytic centre `c`, selects the tangent orientation from the exact
component-1 valuations of `P3`, constructs toric arcs on `C2,C3,C4`, and
expands every ambient function through the requested component multiplicity.
This retains the connected-module cancellations and does not count blow-up
centres as independent scalar rows.

The local adapter lives in
[`scripts/elliptic_neighbor_compiler.sage`](scripts/elliptic_neighbor_compiler.sage)
as `split_multiplicative_toric_chord_condition`.  Its independent regression
uses the universal Tate-normal-form 5-torsion curve: the marked point meets
component 1, the exact correction is `-C2-2C3-C4`, and the resulting matrix
has rank two over `GF(101)`.

## Incoming family contract

The compiler consumes a JSON document validated by
[`../schemas/elkies_k3_lattice_foundry_ns0024_mw4_family.schema.json`](../schemas/elkies_k3_lattice_foundry_ns0024_mw4_family.schema.json).
Its essential shape is:

```json
{
  "schema": "elkies-k3.lattice-foundry-ns0024-mw4-family-modp.v1",
  "status": "PASS_EXACT_MODULAR_NS0024_MW4_FAMILY_MARKING",
  "prime": 101,
  "parameters": ["u"],
  "surface": {
    "A_coefficients_low_to_high": ["..."],
    "B_coefficients_low_to_high": ["..."]
  },
  "sections": {
    "P1": {"x": "...", "y": "..."},
    "P2": {"x": "...", "y": "..."},
    "P3": {"x": "...", "y": "..."},
    "P4": {"x": "...", "y": "..."}
  },
  "marking": {
    "minimum_basis_sha256": "<sha256>",
    "normalized_supports": {"I7": "0", "I5": "1", "I4": "infinity"},
    "section_profiles_I7_I5_I4": {
      "P1": [1, 0, 0],
      "P2": [2, 1, 3],
      "P3": [2, 1, 1],
      "P4": [1, 1, 1]
    },
    "section_intersection_gram": [
      [-2, 1, 2, 1],
      [1, -2, 0, 1],
      [2, 0, -2, 1],
      [1, 1, 1, -2]
    ]
  },
  "proof_boundary": {"proved": "...", "not_proved": "..."}
}
```

Expressions are exact Sage expressions over `GF(p)`, an optional exact finite
extension `GF(p^d)`, or a one-variable function field over either constant
field.  A finite extension is declared by an optional record such as
`"extension": {"generator": "z", "modulus": "z^2 + 1"}`; the compiler
checks irreducibility rather than trusting the declaration.  The names `t`
and `U` are reserved for the old and new elliptic bases.  A single model is
represented by an empty `parameters` array, while a generic modular family
uses one parameter.

The recovery search may first produce an isolated marked point either as a
compact P4 record over `GF(p^d)` plus an oriented prime-field MW3 seed, or as
one joint record carrying the surface and all four sections over `GF(p^d)`.
The independent adapter
[`scripts/adapt_lattice_foundry_ns0024_mw4_point_for_edge1.sage`](scripts/adapt_lattice_foundry_ns0024_mw4_point_for_edge1.sage)
joins those records and replays all four curve equations, the `I7+I5+I4+8I1`
profile, every component label, and the complete section-intersection Gram
matrix before granting the compiler input status.  Thus finite-extension
points need no manual transcription.

A direct prime-field `MW4SEED` hit can be converted without hand editing by
[`scripts/convert_lattice_foundry_ns0024_mw4_seed_to_point.py`](scripts/convert_lattice_foundry_ns0024_mw4_seed_to_point.py).
The converter is lossless only; it does not bypass the adapter's independent
source-marking replay.

For a joint zero-dimensional slice, the preferred extractor is
[`scripts/extract_lattice_foundry_ns0024_joint_rur_point.sage`](scripts/extract_lattice_foundry_ns0024_joint_rur_point.sage).
It decodes every irreducible eliminant factor in its own arbitrary-degree
residue field and invokes the independent joint source verifier before this
adapter.  The older quadratic fixed-MW3 assumption is absent from that path.

The minimum-pole certificate also enumerates the q4-containing basis frontier.
Its resolved component-depth recommendation has absolute profiles
`(6,0,0),(2,1,1),(4,2,0),(6,4,3)` and Gram matrix
`[[-2,1,1,1],[1,-2,0,2],[1,0,-2,2],[1,2,2,-2]]`; in that basis the same
abstract q4/orbit1 horizontal is named `P2`.  The adapter and compiler accept
either pinned basis variant, replay absolute labels, and select the horizontal
from the bound marking.  This avoids confusing absolute component labels with
the relative multipliers returned when `P4` is used as a local generator.

The compiler independently checks:

- short-Weierstrass K3 degree bounds and nonsingularity;
- exact discriminant orders `(7,5,4)` at `0,1,infinity`;
- a squarefree residual discriminant of degree eight;
- all four displayed section identities and their pole profile `(0,0,0,1)`;
- the finite-node incidence fingerprint, all `I7/I5/I4` component labels,
  and the complete section-intersection Gram matrix forced by the marking;
- the resolved component-1 orientation of the bound q4 horizontal (`P3` in
  the original basis, `P2` in the resolved recommendation) at the split `I5`;
- the complete `4 -> 2 -> 2` resolved-RR calculation;
- exact square removal, binary-quartic degree four, finite minimization, the
  infinity boundary, target root data, and Euler number 24.

The output is validated by
[`../schemas/elkies_k3_lattice_foundry_ns0024_edge1_compilation.schema.json`](../schemas/elkies_k3_lattice_foundry_ns0024_edge1_compilation.schema.json).

## Reproduce

Prepare or check the equation-independent handoff:

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/prepare_lattice_foundry_ns0024_edge1_compiler.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/prepare_lattice_foundry_ns0024_edge1_compiler.sage --check
```

Run the compiler regression:

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/verify_elliptic_neighbor_compiler.sage
```

Compile an incoming certified family:

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/compile_lattice_foundry_ns0024_edge1_modp.sage \
  --input artifacts/generated-results/<mw4-family>.json \
  --output artifacts/generated-results/<mw4-family>-edge1.json
```

For a compact marked point emitted by the current residue-algebra recovery,
first form the certified compiler input and then compile it:

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/adapt_lattice_foundry_ns0024_mw4_point_for_edge1.sage \
  --point artifacts/generated-results/<mw4-point>.json \
  --output artifacts/generated-results/<mw4-source>.json

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/compile_lattice_foundry_ns0024_edge1_modp.sage \
  --input artifacts/generated-results/<mw4-source>.json \
  --output artifacts/generated-results/<mw4-source>-edge1.json
```

The same two stages can be dispatched in one fail-closed command:

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/run_lattice_foundry_ns0024_edge1_handoff.sage \
  --input artifacts/generated-results/<mw4-point>.json \
  --source-output artifacts/generated-results/<mw4-source>.json \
  --edge-output artifacts/generated-results/<mw4-source>-edge1.json
```

For the older split compact-point format, add `--seed
artifacts/generated-results/<mw3-seed>.txt`.  For an input already in the
certified family schema, omit both `--seed` and `--source-output`.  Add
`--check` to replay the complete handoff without rewriting either artifact.

Then rerun the same command with `--check`.  Do not update
`MATH_STATUS.json` merely because a modular family compiles: characteristic
zero, source identity, Picard rank, and the equation-side effective child zero
remain separate gates.

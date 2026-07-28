# Proof-carrying arithmetic compilation

The local-to-global modules now have a certificate boundary.  A JSON
specification supplies monic local models, coefficient neighborhoods, a real
coefficient box, an irreducibility prime, a translation, and any available
local action witnesses.  It may also supply a nonnegative `stable_parameter`.
The compiler emits:

- each local model, its discriminant valuation, and the universal precision
  `2*v_p(discriminant)+1`;
- the denominator-aware coefficient CRT moduli and residues;
- the synthesized global polynomial and rational real-root isolating
  intervals;
- a good-prime irreducibility witness;
- the selected admissible translation, target, and inverse-polynomial
  identity;
- the determinant-one map degrees, expanded term counts, and a canonical
  SHA-256 hash of the fully expanded map;
- when `stable_parameter` is present, the cubic boundary-count or
  higher-degree Fitting-support/Newton-area record, together with the `P=0`
  boundary ledger;
- exact local factor and Frobenius-action certificates where the
  specification provides them; and
- for the minimal gauge, the names of the generated Lean specializations.

The first specification is
[`ramified_quintic.json`](specifications/ramified_quintic.json).  It prescribes
a tame ramified cubic and unramified quadratic over `Q_2`, a tame ramified
quadratic and unramified cubic over `Q_3`, inert reduction at `5`, and
Frobenius cycle type `(1,2,2)` at `7`.  Its compiled polynomial is

```text
T^5 - (9855/30241)T^4 + (163265/30241)T^3
    + (190/30241)T^2 + (113214/30241)T - 7266/30241.
```

The checked artifact is
[`arithmetic_keller_quintic.json`](../artifacts/generated-results/arithmetic_keller_quintic.json).
Its expanded determinant-one map has SHA-256
`a67bb7ddc8f0516e5e2ed236695fe896fffcb4125dd40158d2d82be995836a6b`.

Two further artifacts exercise the optional stable-family block:

- [`arithmetic_keller_quintic_stable_m2.json`](../artifacts/generated-results/arithmetic_keller_quintic_stable_m2.json)
  retains the identical ramified quintic, target, and inverse polynomial.
  Its power shift `m=2` has Fitting support
  `[(0,0),(1,2),(6,3),(7,4)]`, normalized area `13`, and expanded-map
  SHA-256
  `e06f20c245099f261f556131c75f01df45a4f1d9491bf89d984119365f3539f6`.
- [`arithmetic_keller_cubic_stable_n7.json`](../artifacts/generated-results/arithmetic_keller_cubic_stable_n7.json)
  compiles the connected cubic `T^3-T-1`.  Its family parameter `k=3`
  gives exponent `n=7`, boundary-component count `7`, and expanded-map
  SHA-256
  `b7ad8dfaca5f3ac4d119a46e49e459de1c4ddb0e7c1f0a77b196159dc8cfe2e2`.

## Independent verification

The dependency-free verifier
[`verify_arithmetic_keller_certificate.py`](../scripts/verify_arithmetic_keller_certificate.py)
uses only Python's standard library.  It implements rational polynomial
arithmetic, resultants, finite-field irreducibility, Sturm isolation, sparse
three-variable arithmetic, and the Jacobian calculation itself.  It imports
neither SymPy nor `jcsearch`.

The PARI/GP verifier
[`verify_arithmetic_keller_certificate.gp`](../scripts/verify_arithmetic_keller_certificate.gp)
reads the same JSON.  `jq` is used only as a JSON parser; all arithmetic,
factorization, Sturm, map reconstruction, and Jacobian work is redone in
PARI/GP.  Both verifiers reconstruct the fully expanded map and independently
recompute its canonical hash.  When the optional stable block is present,
they also recompute its parameter conversion, support, Newton area or
boundary count, and boundary-prime ledger before selecting the corresponding
cubic or power-shifted reconstruction formula.

```bash
python3 scripts/verify_arithmetic_keller_certificate.py
gp -q -f scripts/verify_arithmetic_keller_certificate.gp
python3 scripts/verify_arithmetic_keller_certificate.py \
  artifacts/generated-results/arithmetic_keller_quintic_stable_m2.json
ARITHMETIC_CERTIFICATE=artifacts/generated-results/arithmetic_keller_quintic_stable_m2.json \
  gp -q -f scripts/verify_arithmetic_keller_certificate.gp
python3 scripts/verify_arithmetic_keller_certificate.py \
  artifacts/generated-results/arithmetic_keller_cubic_stable_n7.json
ARITHMETIC_CERTIFICATE=artifacts/generated-results/arithmetic_keller_cubic_stable_n7.json \
  gp -q -f scripts/verify_arithmetic_keller_certificate.gp
```

These are exact algebra and map replays.  The stable functoriality that makes
the recorded boundary count or Newton area a separating invariant is the
written theorem in the corresponding canonical multiplicity note.

## Formal algebra-to-Keller layer

The compiler also generates
[`GeneratedArithmeticQuintic.lean`](../formal/finite-etale-keller/FiniteEtaleKeller/GeneratedArithmeticQuintic.lean).
It proves a denominator-cleared Bézout identity and hence squarefreeness,
checks the selected translation, and instantiates the existing supplied-
translation fiber equivalence, Jacobian-one theorem, inverse-polynomial
identity, and automatic page-one certificate.

```bash
cd formal/finite-etale-keller
lake env lean FiniteEtaleKeller/GeneratedArithmeticQuintic.lean
```

The Lean file formalizes the algebra-to-Keller implication for the explicit
polynomial.  The local-field stability theorem and the arithmetic input
claims are replayed by the two exact verifiers; they are not silently claimed
as Lean formalizations.  The two stable-lift artifacts deliberately omit a
`lean_instantiation` block: the present Lean development formalizes the
minimal gauge, not the fiber-invisible cubic or common power-shift formulas.

## Regeneration

Regeneration is explicit because the JSON and Lean outputs are pinned
certificates:

```bash
.venv/bin/python scripts/compile_arithmetic_keller_certificate.py
.venv/bin/python scripts/compile_arithmetic_keller_certificate.py \
  --stable-parameter 2 \
  --certificate artifacts/generated-results/arithmetic_keller_quintic_stable_m2.json
.venv/bin/python scripts/compile_arithmetic_keller_certificate.py \
  --spec arithmetic/specifications/connected_cubic_stable_n7.json \
  --certificate artifacts/generated-results/arithmetic_keller_cubic_stable_n7.json
```

Use `make verify-arithmetic-compilation` for all seven non-mutating checks and
`make refresh-arithmetic-compilation` only when intentionally updating the
three generated JSON files or the legacy Lean file.

The pinned generation environment uses Python 3 with SymPy `1.14.0` from
`requirements.txt`.  The second replay was checked with PARI/GP `2.17.4` and
`jq`; the formal replay uses the Lean toolchain pinned under
`formal/finite-etale-keller/`.  The current whole-file SHA-256 values are
recorded in `artifacts/generated-results/README.md`.

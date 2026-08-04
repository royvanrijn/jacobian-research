# Lean audit of the GVC manuscript

This package is a partial formalization of
[`../../papers/generalized-vanishing-two-variables/main.tex`](../../papers/generalized-vanishing-two-variables/main.tex).
It is pinned to Lean `v4.32.1` and the matching Mathlib release.

Build it with:

```sh
lake build
```

## What is checked

- `Definitions` gives a coefficientwise semantics for arbitrary
  constant-coefficient differential symbols, proves agreement with
  Mathlib's formal partial derivative on coordinate symbols, and states GVC
  in that model.
- `CuspIdentity` proves `x*C = rho^3 - t^2*A^2` over every commutative ring.
- `EndpointCoefficients` proves the all-order adjacent-coefficient step:
  it constructs the normalized polynomial primitive, proves the order
  `2m+1` endpoint flatness from the actual derivative
  `(1-(1+u)^2)^(2m)`, and then proves that the pure coefficient is zero while
  the neighboring mixed coefficient is `c_m`.  It also proves nonvanishing
  of the displayed exact rational scalar.
- `ConcreteWitness` defines the manuscript's literal `rho`, `A`, `C`, `P`,
  `Delta`, `Lambda`, and `Q` in `Q[x,y,t]`.  From the explicitly named
  Reynolds/phase bridge it derives the pure identities, the exact mixed
  value, and failure of GVC for this pair.
- `Envelope` checks the intermediate-value step which turns a positive
  continuous envelope gap and a later nonpositive value into a common
  threshold.
- `BinaryReduction` states the common-threshold support inequalities and
  proves the binary theorem from two separately named obligations: envelope
  closure and common-threshold terminality.

There are no `sorry`, `admit`, or explicit `axiom` declarations.
`ConcreteCounterexampleBridge` and `BinaryEnvelopeBridge` are structures of
hypotheses; Lean does not construct values of them in the current package.

## What is not yet verified

This package is not a complete formal proof of the paper.  The following
load-bearing arguments remain to be formalized:

1. the beta evaluation of `c_m`, the Reynolds--apolar identity, and quadric
   phase extraction connecting the concrete ternary polynomials to
   `endpointKernel`;
2. the exposing-point consequence of Duistermaat--van der Kallen and the
   good-prime shifted-ray separation theorem;
3. binary Hall localization, no-reversal of the moving Newton intervals, and
   the full common-threshold cutoff;
4. the arbitrary cusp-profile family, unused-variable padding, and the final
   all-dimension classification.

Accordingly, `MATH_STATUS.json` should continue to record
`formal_verification: false` for the paper's headline claims until those
bridges are inhabited by proofs.

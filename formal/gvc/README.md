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
- `TopContraction` proves that multiplication of differential symbols is
  composition of their actions and that equal-degree homogeneous symbols
  contract homogeneous polynomials to the coefficientwise apolar pairing.
- `CuspIdentity` proves `x*C = rho^3 - t^2*A^2` over every commutative ring.
- `EndpointCoefficients` proves the all-order adjacent-coefficient step:
  it constructs the normalized polynomial primitive, proves the order
  `2m+1` endpoint flatness from the actual derivative
  `(1-(1+u)^2)^(2m)`, and then proves that the pure coefficient is zero while
  the full neighboring ladder has coefficient
  `choose (m-1) (ell-1) * c`.  The same proof is instantiated for an
  arbitrary rational endpoint profile.  A formal polynomial integration
  argument proves the beta evaluation
  `c_m = 2^(2m) (2m)! / (4m+1)!!`, not merely its nonvanishing, and the file
  also proves nonvanishing of the displayed exact rational scalar.
- `ConcreteWitness` defines the manuscript's literal `rho`, `A`, `C`, `P`,
  `Delta`, `Lambda`, and `Q` in `Q[x,y,t]`.  It proves that `P` and `Lambda`
  are nonzero homogeneous polynomials of degree/order twelve.  Its
  normalized algebraic Reynolds functional is the coefficientwise top
  contraction, so the differential identity itself is proved rather than
  assumed.  It factors the final argument through an explicitly named
  quadric phase-extraction bridge.
- `ReynoldsExpansion` expands `Delta^k` coefficientwise, evaluates every
  normalized diagonal weight by a formal beta integral, and proves the
  algebraic Reynolds--phase identity for arbitrary ternary polynomials.
- `QuadricPhase` constructs the Laurent restriction
  `x -> z`, `y -> (1-t^2)z^-1`, proves the paper's exact formulas for
  `rho`, `A`, `C`, and `P`, and identifies the Laurent constant term of every
  homogeneous input with the Reynolds phase polynomial.
- `PhaseKernel` performs the even-phase coefficient extraction for `P^m`
  and `x^2 P^m`, proves the coefficientwise algebraic change of height
  variable, identifies the result with `endpointKernel`, and constructs a
  value of `ConcreteCounterexampleBridge`.  Thus the pure identity, exact
  mixed scalar, and rational GVC counterexample are unconditional Lean
  theorems.
- `Padding` proves that injective renaming commutes with the differential
  action, that GVC descends along a variable embedding (so the binary theorem
  implies the unary theorem), and that the ternary counterexample carries to
  every dimension `n >= 3`.
- `BaseChange` proves the analogous compatibility with coefficient-ring
  maps.  It transports the pure identity, exact mixed scalar, and
  counterexample from `Q` to every characteristic-zero field and composes
  this with unused-variable padding.  `VerifiedCounterexample` supplies the
  constructed bridge and proves unconditionally that GVC fails over every
  characteristic-zero field in every finite dimension `n >= 3`.
- `ProfileFamily` defines the manuscript's literal
  winding--profile--radial family, proves the degree/order formula
  `2 * (6r + 3e + h)`, and checks the minimal and radial specializations.
  It also derives the pure identity, every exact multiplier-ladder scalar,
  and failure of GVC from a bridge containing only the remaining
  multivariate phase-extraction equalities and the paper's nonzero-moment
  hypothesis.
- `FactorialValuation` proves over the rational `p`-adic valuation both the
  coordinatewise lower bound and the exact prime-block formula in Lemma 3.2,
  combines them into its two-coordinate factorial product, and proves the
  floor-sum estimate used for the non-Frobenius terms in Proposition 4.1.
- `Envelope` checks the intermediate-value step which turns a positive
  continuous envelope gap into a common threshold.  It also proves that the
  strictly negative slope on the final affine envelope pieces supplies the
  required later nonpositive value.
- `BinaryReduction` states the common-threshold support inequalities and
  proves the binary theorem from two separately named obligations: envelope
  closure and common-threshold terminality.

There are no `sorry`, `admit`, or explicit `axiom` declarations.
`verifiedConcreteCounterexampleBridge` constructs the concrete bridge.
`ProfileFamilyBridge` and `BinaryEnvelopeBridge` remain structures of
hypotheses for the arguments not yet formalized.

## What is not yet verified

This package is not a complete formal proof of the paper.  The following
load-bearing arguments remain to be formalized:

1. the exposing-point consequence of Duistermaat--van der Kallen, transfer
   of the rational `p`-adic factorial calculation to an unramified
   number-field prime, and the full good-prime shifted-ray separation
   theorem;
2. binary Hall localization, no-reversal of the moving Newton intervals, and
   the full common-threshold cutoff;
3. the multivariate phase bridge for the arbitrary cusp-profile family and
   the positive half of the final dimension classification.  The
   counterexample's coefficient base change and unused-variable padding are
   now checked.

Theorem 8.1 and the negative `n >= 3` half of Theorem 10.1 are now fully
Lean-verified.  `MATH_STATUS.json` records `formal_verification: true` for
the canonical concrete-counterexample entry `GVC3HC`; it remains false for
the paper's binary, full-profile, and complete-classification entries.

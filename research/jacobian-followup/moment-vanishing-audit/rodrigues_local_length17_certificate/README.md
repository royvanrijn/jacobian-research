# Exact Rodrigues local-algebra replay

Requirements: Python 3 and SymPy. All coefficients and matrix operations use
exact integers and rational numbers; no numerical approximation is used.

From this directory run:

```sh
python verify_certificate.py > verification.log 2>&1
```

The verifier imports and reruns the four calculation scripts in sequence:

1. `rod_local_taylor_quartic.py`: sphere moments, the first-ten Jacobian,
   formal corrections through degree three, and fourth-order residuals.
2. `rod_local_taylor_quintic.py`: corrections and residuals through degree six.
   The original analysis-time guard has been removed; the sixth-order step
   always runs.
3. `rod_local_taylor_septic.py`: an independently checked one-sum marked-moment
   evaluator, the sixth implicit correction, and seventh-order residuals.
4. `rod_local_full_prefix28.py`: residual moments 15 through 28 and incremental
   exact ideal ranks modulo `(a,b)^8`.

The verifier then checks the transverse determinant, the weighted homogeneous
initial ideals, the two reduced nonzero coefficients, and independently built
32-by-36 and 50-by-36 rational matrices. It constructs multiplication by `a`
and `b` in the truncated quotient and checks its Hilbert function and exact
maximal-ideal nilpotence index. Direct row-membership checks include every
ordinary degree-seven monomial, every monomial of weight at least thirteen
within the truncation (weights `a=3,b=2`), and nonmembership of `b^6`.

Output files:

- `verification.log`: run output and assertion status.
- `certificate_data.json`: exact residual coefficients, both rational matrices
  with row labels, weighted footprints, a quotient basis, multiplication
  matrices, and the Hilbert function. Rational values are JSON strings.

In the JSON matrices, a row label `[n,u,v]` means the coefficients of
`a^u*b^v*R_n` modulo `(a,b)^8`. Columns use the explicit `monomials` array.
Polynomial coefficient records have the form `[[a_exponent,b_exponent],value]`.

Scope: these scripts certify the exact finite calculation. The mathematical
argument must additionally justify the formal transverse interpretation,
the passage by Nakayama's lemma from checked degree-seven membership to the
completed local algebra, and the filtered tail bound making the uncomputed
moments redundant. These are local calculations and do not assert a global
classification of pure phases.
